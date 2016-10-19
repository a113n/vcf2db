import pandas as pd
import csv

mart_export1 = pd.read_csv("Homo_sapiens.GRCh38.86.mart_export.gz", header = 0, sep='\t', skip_blank_lines=True)
mart_export2 = pd.read_csv("Homo_sapiens.GRCh38.86.cds_length.mart_export.gz", header = 0, sep='\t', skip_blank_lines=True)
pheno_export = pd.read_csv("HMD_HumanPhenotype.rpt", header = None, usecols = [1,5], names=['EntrezGene ID','mph'], sep='\t', skip_blank_lines=True)
rvis_export = pd.read_csv("GenicIntolerance_v3_12Mar16.txt", header = 0,  usecols = [0, 1, 2],  sep='\t', skip_blank_lines=True)
hgnc_export = pd.read_csv("HGNC.txt", header = 0,  usecols = [0,3,4], sep='\t', skip_blank_lines=True)

#Parse HGNC alias and previous names into a dict of comma separated values
hgnc_dict = {}
infile = open("HGNC.txt")
reader = csv.reader(infile, delimiter='\t')
header = reader.__next__()
for row in reader:        
    if not row[0] in hgnc_dict:
        hgnc_dict[row[0]]={"Synonyms":set(),"HGNC Symbol":row[1]}
    
    if row[3] != "":
        hgnc_dict[row[0]]["Synonyms"].add(row[3])
    if row[4] != "":
        hgnc_dict[row[0]]["Synonyms"].add(row[4])

hgnc_export = {"HGNC ID":[], "HGNC Symbol":[], "Synonyms":[]}
for key,value in hgnc_dict.items():
    hgnc_export["HGNC ID"].append(key)
    hgnc_export["HGNC Symbol"].append(value["HGNC Symbol"])
    hgnc_export["Synonyms"].append(",".join(sorted(value["Synonyms"])))
    
hgnc_export = pd.DataFrame.from_dict(hgnc_export)

# Dataframe postprocessing
mart_export2["Protein_length"] = mart_export2["CDS Length"]/3-1
mart_export2 = mart_export2.round({"CDS Length": 0, "Protein_length": 0})
mart_export2 = mart_export2.drop('Ensembl Gene ID', 1)

pheno_export = pheno_export.dropna().drop_duplicates(subset="EntrezGene ID")

# Left joining tables
merged_df = pd.merge(mart_export1,mart_export2,how='left',on='Ensembl Transcript ID')
merged_df = pd.merge(merged_df,hgnc_export,how='left',left_on='HGNC ID(s)',right_on="HGNC ID")
merged_df = pd.merge(merged_df,pheno_export,how='left',on='EntrezGene ID')
merged_df = pd.merge(merged_df,rvis_export,how='left',left_on='HGNC Symbol',right_on="GENE")
merged_df["hgnc_flag"] = merged_df["HGNC ID(s)"].apply(lambda x: 0 if pd.isnull(x) else 1)
merged_df["is_hgnc_id"] = merged_df["HGNC ID(s)"].apply(lambda x: "None" if pd.isnull(x) else x)

# Output
merged_df.to_csv('detailed_gene_table_v86', sep="\t", index=False,
                 header=["Chromosome","Gene_name","Is_hgnc","Ensembl_gene_id","Ensembl_transcript_id","Biotype",
                              "Transcript_status","CCDS_id","HGNC_id","CDS_length","Protein_length",
                              "Transcript_start","Transcript_end","strand","Synonyms", 
                              "Rvis_pct","entrez_gene_id","mammalian_phenotype_id"],
                 columns=["Chromosome Name","Associated Gene Name","hgnc_flag",
                          "Ensembl Gene ID","Ensembl Transcript ID","Transcript type",
                          "Status (gene)","CCDS ID","is_hgnc_id","CDS Length","Protein_length",
                          "Transcript Start (bp)","Transcript End (bp)","Strand","Synonyms",
                          "%ALL_0.01%","EntrezGene ID","mph"])
merged_df.to_csv('summary_gene_table_v86', sep="\t", index=False,
                 header=["Chromosome","Gene_name","Is_hgnc","Ensembl_gene_id",
                         "HGNC_id","Synonyms", "Rvis_pct","Strand","Mammalian_phenotype_id"],
                 columns=["Chromosome Name","Associated Gene Name","hgnc_flag",
                          "Ensembl Gene ID","is_hgnc_id","Synonyms","%ALL_0.01%",
                          "Strand","mph"])
