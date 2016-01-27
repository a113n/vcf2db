"""
Take a VCF and create a gemini compatible database
"""
import sys

import itertools as it
import re
import zlib
import cPickle
import time
from collections import defaultdict
import multiprocessing

import numpy as np
import sqlalchemy as sql
from pedagree import Ped
import geneimpacts
from sqlalchemy import String, Float, Integer, Boolean
import cyvcf2

import cProfile
import StringIO
import pstats
import contextlib

def grouper(n, iterable):
    iterable = iter(iterable)
    piece = list(it.islice(iterable, n))
    while piece:
        yield piece
        piece = list(it.islice(iterable, n))

@contextlib.contextmanager
def profiled():
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    s = StringIO.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('time')
    ps.print_stats(60)
    # uncomment this to see who's calling what
    # ps.print_callers()
    print s.getvalue()

def set_column_length(e, column, length, saved={}):
    table = column.table
    c = column.table.columns[column.name]
    if c.type.length >= length:
        return
    c.type.length = length
    column.type.length = length
    if saved.get((table.name, c.name), 0) < length:
        sys.stderr.write("changing varchar field '%s' to length %d\n" %
                                     (c.name,  length))
    saved[(table.name, c.name)] = c.type.length
    if e.dialect.name.startswith("postgres"):
        e.execute('ALTER TABLE %s ALTER COLUMN %s TYPE VARCHAR(%d)' %
                            (table.name, c.name, length))
    elif e.dialect.name == "mysql":
        e.execute('ALTER TABLE %s MODIFY %s VARCHAR(%d)' %
                                (table.name, c.name, length))

#import blosc
#def pack_blob(obj):
#    if obj is None: return ''
#    return buffer(blosc.compress(obj.tostring(), obj.dtype.itemsize, clevel=5, shuffle=True))

def pack_blob(obj, _none=zlib.compress(cPickle.dumps(None, cPickle.HIGHEST_PROTOCOL))):
    if obj is None: return _none
    return zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 1)

def clean(name):
    """
    turn a vcf id into a db name
    """
    return name.replace("-", "_").replace(" ", "_").strip('"').strip("'").lower()

def info_parse(line,
        _patt=re.compile("(\w+)=(\"[^\"]+\"|[^,]+)")):
    """
    >>> ret = info_parse('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">')
    >>> assert ret == {"ID": "AC", "Number": "A", "Type": "Integer","Description": '"Allele count in genotypes, for each ALT allele, in the same order as listed"'}, ret
    """
    assert line.startswith("##INFO=")
    stub = line.split("=<")[1].rstrip(">")
    return dict(_patt.findall(stub))


type_lookups = {
        "Integer": sql.Integer(),
        "Float": sql.Float(),
        "Flag": sql.Boolean(),
        "Character": sql.String(1),
        "String": sql.String(5),
        }

class VCFDB(object):
    gt_cols = ("gts", "gt_types", "gt_phases", "gt_depths", "gt_ref_depths",
               "gt_alt_depths", "gt_quals")

    effect_list = ["CSQ", "ANN", "EFF"]
    _black_list = []

    def __init__(self, vcf_path, db_path, ped_path=None, blobber=pack_blob,
            black_list=None):
        self.vcf_path = vcf_path
        if not db_path.startswith(("sqlite:", "mysql", "postgres")):
                db_path = "sqlite:///" + db_path
        self.db_path = db_path
        self.engine = sql.create_engine(db_path, poolclass=sql.pool.NullPool)
        self.impacts_headers = {}
        self.metadata = sql.MetaData(bind=self.engine)

        self.blobber = blobber
        self.ped_path = ped_path
        self.black_list = list(VCFDB._black_list) + list(VCFDB.effect_list) + (black_list or [])
        self.pool = multiprocessing.Pool(2)

        self.vcf = cyvcf2.VCF(vcf_path)
        # we use the cache to infer the lengths of string fields.
        self.cache = it.islice(self.vcf, 10000)
        self.create_columns()
        self.samples = self.create_samples()
        self.load()
        self.index()

    def _set_variant_properties(self, v, d):
        d['type'] = v.var_type
        d['sub_type'] = v.var_subtype
        d['call_rate'] = v.call_rate
        d['num_hom_ref'] = v.num_hom_ref
        d['num_het'] = v.num_het
        d['num_hom_alt'] = v.num_hom_alt
        d['aaf'] = v.aaf
        # TODO:
        #d['hwe'] = v.hwe
        #d['inbreeding_coef'] = ??
        #d['pi'] = ??

    def _load(self, iterable, create, start):

        variants = []
        keys = set()
        i = None
        for i, v in enumerate(iterable, start=start):
            d = dict(v.INFO)
            for c in self.gt_cols:
                # named gt_bases in cyvcf2 and gts in db
                arr = v.gt_bases if c == "gts" else getattr(v, c, None)
                if arr is not None:
                    arr = arr[self.sample_idxs]
                #d[c] = self.blobber(arr)
                d[c] = arr

            d['chrom'], d['start'], d['end'] = v.CHROM, v.start, v.end
            d['ref'], d['alt'] = v.REF, ",".join(v.ALT)

            d['qual'], d['filter'], d['vcf_id'] = v.QUAL, v.FILTER, v.ID
            d['variant_id'] = i
            self._set_variant_properties(v, d)

            # TODO: just save required keys outside.
            keys.update(d.keys())

            variants.append(d)
            # http://docs.sqlalchemy.org/en/latest/faq/performance.html
            if not create and (i % 10000) == 0:
                self.insert(variants, keys, i)
                variants = variants[:0]

        if len(variants) != 0:
            self.insert(variants, keys, i, create=create)

        return i

    def load(self):
        self.t0 = self.t = time.time()

        i = self._load(self.cache, create=True, start=1)
        self.cache = []
        #with profiled():
        self._load(self.vcf, create=False, start=i+1)

    def check_column_lengths(self, dicts, cols):
        change_cols = defaultdict(int)
        for name, c in cols.items():
            l = c.type.length
            for d in dicts:
                if len(d.get(name) or '') > l:
                    change_cols[c.name] = max(change_cols[c.name], len(d.get(name)))
        return dict(change_cols)

    def insert(self, variants, keys, i, create=False):
        ivariants, variant_impacts = [], []
        te = time.time()
        for variant, impacts in self.pool.imap(gene_info, ((v,
                     self.impacts_headers, self.blobber, self.gt_cols, keys) for
                     v in variants), 10):
            variant_impacts.extend(impacts)
            ivariants.append(variant)
        te = time.time() - te

        variants = ivariants
        vlengths = vilengths = {}

        if create:
            self.create(variants, variant_impacts)

        elif self.engine.dialect.name != "sqlite":
            vlengths = self.check_column_lengths(variants, {c.name: c for c in self.variants_columns if
                c.type.__class__.__name__ == "String"})

            vilengths = self.check_column_lengths(variant_impacts, {c.name: c for c in
                self.variant_impacts_columns if c.type.__class__.__name__ ==
                "String"})

        self._insert(vlengths, variants,
                     vilengths, variant_impacts)

        fmt = "%d variant_impacts:%d\teffects time: %.1f\tchunk time:%.1f\t%.2f variants/second\n"
        vps = i / float(time.time() - self.t0)
        sys.stderr.write(fmt % (i, len(variant_impacts), te, time.time() - self.t, vps))
        self.t = time.time()

    def _insert(self, vlengths, v_objs, vilengths, vi_objs):

        for name, clen in vlengths.items():
            col = self.variants.columns[name]
            set_column_length(self.engine, col, clen)

        self.__insert(v_objs, self.metadata.tables['variants'].insert())

        for name, clen in vilengths.items():
            col = self.variant_impacts.columns[name]
            set_column_length(self.engine, col, clen)

        self.__insert(vi_objs, self.metadata.tables['variant_impacts'].insert())


    def __insert(self, objs, stmt):

        tx = time.time()
        # (2006, 'MySQL server has gone away'
        # if you see this, need to increase max_allowed_packet and/or other
        # params in my.cnf (or we should detect and reduce the chunk size)
        if len(objs) > 6000:
            for group in grouper(5000, objs):
                g = list(group)
                try:
                    self.engine.execute(stmt, g)
                except:
                    with self.engine.begin() as trans:
                        for o in g:
                            trans.execute(stmt, o)
                    raise
        else:
            try:
                self.engine.execute(stmt, objs)
            except:
                with self.engine.begin() as trans:
                    for o in objs:
                        trans.execute(stmt, o)
                raise
        return time.time() - tx

    def create_columns(self):
        self.variants_columns = self.get_variants_columns()
        self.variant_impacts_columns = self.get_variant_impacts_columns()

    def create(self, dvariants, dvariant_impacts):
        # update the lengths of the string columns based on the variants that
        # we've seen so far
        v_cols = {c.name: c for c in self.variants_columns if c.type.__class__.__name__ == "String"}
        self._create(dvariants, v_cols)

        vi_cols = {c.name: c for c in self.variant_impacts_columns if c.type.__class__.__name__ == "String"}
        self._create(dvariant_impacts, vi_cols)

        self._create_tables()

    def _create(self, dicts, cols):
        exclude_cols = set()
        for name, col in cols.items():
            if name in exclude_cols: continue
            for d in dicts:
                try:
                    if col.type.length < len(d.get(name) or ''):
                        #col.type.length = int(1.618 * len(d[name]) + 0.5)
                        col.type.length = int(1.2 * len(d[name]) + 0.5)
                except:
                    print name, col.type.length
                    raise
                if col.type.length > 48:
                    col.type = sql.TEXT()
                    exclude_cols.add(name)
                    break

    def _create_tables(self):
        self.variant_impacts = sql.Table("variant_impacts", self.metadata, *self.variant_impacts_columns)
        self.variant_impacts.drop(checkfirst=True)

        self.variants = sql.Table("variants", self.metadata, *self.variants_columns)
        self.variants.drop(checkfirst=True)

        self.variants.create()
        self.variant_impacts.create()
        self.create_vcf_header_table()

    def create_vcf_header_table(self):
        h = self.vcf.raw_header
        t = sql.Table("vcf_header", self.metadata,
                      #sql.Column("vcf_header", sql.TEXT(len(h)))
                      sql.Column("vcf_header", sql.TEXT)
                      )
        t.drop(self.engine, checkfirst=True)
        t.create()
        self.engine.execute(t.insert(), [dict(vcf_header=h)])

    def get_variant_impacts_columns(self):
        return [sql.Column("variant_id", Integer,
                           sql.ForeignKey("variants.variant_id"), nullable=False),
                ] + self.variants_gene_columns()

    def index(self):
        sys.stderr.write("indexing ... ")
        t0 = time.time()
        sql.Index("idx_variants_chrom_start", self.variants.c.chrom, self.variants.c.start).create()
        sql.Index("idx_variants_exonic", self.variants.c.is_exonic).create()
        sql.Index("idx_variants_coding", self.variants.c.is_coding).create()
        sql.Index("idx_variants_impact", self.variants.c.impact).create()
        sql.Index("idx_variants_impact_severity", self.variants.c.impact_severity).create()
        sys.stderr.write("finished in %.1f seconds...\n" % (time.time() - t0))
        sys.stderr.write("total time: in %.1f seconds...\n" % (time.time() - self.t0))

    def create_samples(self):
        p = Ped(self.ped_path)
        samples = self.vcf.samples
        cols = ['sample_id', 'family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
        idxs = []
        rows = []
        for i, s in enumerate(p.samples(), start=1):
            idxs.append(samples.index(s.sample_id))
            assert s.sample_id in samples
            if i == 0:
                cols.extend(s.attrs)
            rows.append([i, s.family_id, s.sample_id, str(s.paternal_id), str(s.maternal_id),
                '1' if s.sex == 'male' else '2' if s.sex == 'female' else '-9',
                '2' if s.affected is True else '1' if s.affected is False else '-9',
                ] + s.attrs)

        scols = [sql.Column('sample_id', Integer, primary_key=True)]
        for i, col in enumerate(cols[1:], start=1):
            vals = [r[i] for r in rows]
            l = max(len(v) for v in vals)
            scols.append(sql.Column(col, String(l)))

        t = sql.Table('samples', self.metadata, *scols)
        t.drop(checkfirst=True)
        t.create()

        self.engine.execute(t.insert(), [dict(zip(cols, r)) for r in rows])

        # track the order to pull from the genotype fields.
        self.sample_idxs = np.array(idxs)

    def get_variants_columns(self):
        columns = self.variants_default_columns()
        columns.extend(self.variants_calculated_columns())
        columns.extend(self.variants_gene_columns())
        columns.extend(self.variants_sv_columns())
        columns.extend(self.variants_info_columns(self.vcf.raw_header))
        columns.extend(self.variants_genotype_columns())
        return columns

    def variants_default_columns(self):
        return [
            sql.Column("variant_id", Integer(), primary_key=True),
            sql.Column("chrom", String(10)),
            sql.Column("start", Integer()),
            sql.Column("end", Integer()),
            sql.Column("vcf_id", String(12)),
            #sql.Column("anno_id", Integer()),
            sql.Column("ref", sql.TEXT()),
            sql.Column("alt", sql.TEXT()),
            sql.Column("qual", Float()),
            sql.Column("filter", String(10)),
           ]

    def variants_gene_columns(self):
        # all of these are also stored in the variant_impacts table.
        return [
            sql.Column("gene", String(20)),
            sql.Column("transcript", String(20)),
            sql.Column("is_exonic", Boolean()),
            sql.Column("is_coding", Boolean()),
            sql.Column("is_lof", Boolean()),
            sql.Column("is_splicing", Boolean()),
            sql.Column("exon", String(8)),
            sql.Column("codon_change", sql.TEXT()),
            sql.Column("aa_change", sql.TEXT()),
            sql.Column("aa_length", String(8)),
            sql.Column("biotype", String(50)),
            sql.Column("impact", String(20)),
            sql.Column("impact_so", String(20)),
            sql.Column("impact_severity", String(4)),
            sql.Column("polyphen_pred", String(20)),
            sql.Column("polyphen_score", Float()),
            sql.Column("sift_pred", String(20)),
            sql.Column("sift_score", Float()),
            ]


    def variants_calculated_columns(self):
        return [
            sql.Column("type", String(8)),
            sql.Column("sub_type", String(20)),
            sql.Column("call_rate", Float()),
            sql.Column("num_hom_ref", Integer()),
            sql.Column("num_het", Integer()),
            sql.Column("num_hom_alt", Integer()),
            sql.Column("aaf", Float()),
            sql.Column("hwe", Float()),
            sql.Column("inbreeding_coef", Float()),
            sql.Column("pi", Float()),
           ]

    def variants_sv_columns(self):
        return []

    def variants_genotype_columns(self):
        return [sql.Column(name, sql.LargeBinary()) for name in self.gt_cols]

    def update_impacts_headers(self, hdr_dict):
        "keep the description so we know how to parse the CSQ/ANN fields"

        desc = hdr_dict["Description"]
        if hdr_dict["ID"] == "ANN":
            parts = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]
        elif hdr_dict["ID"] == "EFF":
            parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
        elif hdr_dict["ID"] == "CSQ":
            parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
        else:
            raise Exeception("don't know how to use %s as annotation" % hdr_dict["ID"])
        self.impacts_headers[hdr_dict["ID"]] = parts

    def variants_info_columns(self, raw_header):
        "create Column() objects for each entry in the info field"
        for l in (x.strip() for x in raw_header.split("\n")):
            if not l.startswith("##INFO"):
                continue

            d = info_parse(l)
            if d["ID"] in self.effect_list:
                self.update_impacts_headers(d)
                continue

            id = clean(d["ID"])
            if d["ID"] in self.black_list or id in self.black_list:
                continue

            if id == "id":
                id = "idx"
            c = sql.Column(id, type_lookups[d["Type"]], primary_key=False)
            if id.endswith(("_af", "_aaf")) or id.startswith(("af_", "aaf_", "an_")):
                c = sql.Column(id, Float(), default=-1.0)
            yield c

def gene_info(d_and_impacts_headers):
    # this is parallelized as it's only simple objects and the gene impacts
    # stuff is slow.
    d, impacts_headers, blobber, gt_cols, req_cols = d_and_impacts_headers
    impacts = []
    for k in (eff for eff in ("CSQ", "ANN", "EFF") if eff in d):
        if k == "CSQ":
            impacts.extend(geneimpacts.VEP(e, impacts_headers[k], checks=False) for e in d[k].split(","))
        elif k == "ANN":
            impacts.extend(geneimpacts.SnpEff(e, impacts_headers[k]) for e in d[k].split(","))
        elif k == "EFF":
            impacts.extend(geneimpacts.OldSnpEff(e, impacts_headers[k]) for e in d[k].split(","))
        del d[k] # save some memory

    top = geneimpacts.Effect.top_severity(impacts)
    if isinstance(top, list):
        top = top[0]
    keys = ('gene', 'transcript', 'is_exonic', 'is_coding', 'is_splicing',
            'is_lof', 'exon', 'codon_change', 'aa_change', 'aa_length',
            'biotype', 'top_consequence', 'so', 'effect_severity',
            'polyphen_pred', 'polyphen_score', 'sift_pred', 'sift_score')


    for k in keys:
        d[k] = getattr(top, k)

    d['impact'] = top.top_consequence
    d['impact_so'] = top.so
    d['impact_severity'] = top.effect_severity

    for c in gt_cols:
        d[c] = blobber(d[c])

    # add what we need.
    u = dict.fromkeys(req_cols)
    u.update(d)
    d = u

    # TODO: check "exonic" vs is_exonic"
    # TODO: check top_consequence
    gimpacts = []
    for impact in impacts:
        #gimpacts.append({k: getattr(impact, k) for k in keys})
        gimpacts.append(dict(variant_id=d['variant_id'],
                             gene=impact.gene, transcript=impact.transcript,
                             is_exonic=impact.is_exonic, is_coding=impact.is_coding,
                             is_splicing=impact.is_splicing, is_lof=impact.is_lof,
                             exon=impact.exon, codon_change=impact.codon_change,
                             aa_change=impact.aa_change, aa_length=impact.aa_length,
                             biotype=impact.biotype, top_consequence=impact.top_consequence,
                             so=impact.so, effect_severity=impact.effect_severity,
                             polyphen_pred=impact.polyphen_pred,
                             polyphen_score=impact.polyphen_score,
                             sift_pred=impact.sift_pred,
                             sift_score=impact.sift_score))
    return d, gimpacts

if __name__ == "__main__":

    import doctest
    doctest.testmod()

    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("VCF")
    p.add_argument("db")
    p.add_argument("ped")
    a = p.parse_args()

    v = VCFDB(a.VCF, a.db, a.ped)