#!/usr/bin/python
import os
from tqdm import tqdm
from urllib.request import urlopen
"""
This script is for updating the files below:
    alias files
    genes_<genome>.bed
    chrom.sizes.<genome>
"""
def get_num_lines(file_path):
    with open(file_path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

##########################################################################
###  HG38
###  -- Update hg19/alias_human.txt
###  -- Update hg38/alias_human.txt
###  -- Update hg38/genes_hg38.bed
##########################################################################

alias = []
genes = []
file_path = "hg38/gencode.v32.annotation.gtf"
print("Parsing "+file_path)
with open(file_path) as gtf_hg38:
    for line in tqdm(gtf_hg38, total=get_num_lines(file_path)):
        if line.startswith("#"):
            pass
        else:
            line = line.replace('"', "")
            line = line.replace(';', "")
            l = line.split()
            if l[2] == "gene":
                id_symbol = l.index("gene_name")
                id_ensembl = l.index("gene_id")
                if "hgnc_id" in l:
                    tag_HGNC = l[l.index("hgnc_id")+1]
                else:
                    tag_HGNC = ""
                if "havana_gene" in l:
                    tag_Havana = l[l.index("havana_gene")+1]
                else:
                    tag_Havana = ""

                alias.append([l[id_symbol + 1], l[id_ensembl + 1],
                              tag_HGNC, tag_Havana])
                genes.append([l[0], l[3], l[4], l[id_symbol + 1], "0", l[6]])

res_alias = list(set(tuple(g) for g in alias))
res_genes = list(set(tuple(g) for g in genes))
# Save alias to HG38
alias_file = "hg38/alias_human.txt"
with open(alias_file, "w") as f:
    for g in res_alias:
        print("\t".join(g), file=f)

# Save alias to HG19
alias_file = "hg19/alias_human.txt"
with open(alias_file, "w") as f:
    for g in res_alias:
        print("\t".join(g), file=f)

# Save genes to HG38
genes_file = "hg38/genes_hg38.txt"
with open(genes_file, "w") as f:
    for g in res_genes:
        print("\t".join(g), file=f)

##########################################################################
###  HG19
###  -- Update hg19/genes_hg19.bed
##########################################################################

genes = []
file_path = "hg19/gencode.v19.annotation.gtf"
print("Parsing "+file_path)
with open(file_path) as gtf_hg19:
    for line in tqdm(gtf_hg19, total=get_num_lines(file_path)):
        if line.startswith("#"):
            pass
        else:
            line = line.replace('"', "")
            line = line.replace(';', "")
            l = line.split()
            if l[2] == "gene":
                id_symbol = l.index("gene_name")
                genes.append([l[0], l[3], l[4], l[id_symbol + 1], "0", l[6]])

res_genes = list(set(tuple(g) for g in genes))

# Save genes to HG19
genes_file = "hg19/genes_hg19.txt"
with open(genes_file, "w") as f:
    for g in res_genes:
        print("\t".join(g), file=f)


##########################################################################
###  MM10
###  -- Update mm9/alias_mouse.txt
###  -- Update mm10/alias_mouse.txt
###  -- Update mm10/genes_mm10.bed
##########################################################################

alias = []
genes = []
file_path = "mm10/gencode.vM23.annotation.gtf"
print("Parsing "+file_path)
with open(file_path) as gtf_mm10:
    for line in tqdm(gtf_mm10, total=get_num_lines(file_path)):
        if line.startswith("#"):
            pass
        else:
            line = line.replace('"', "")
            line = line.replace(';', "")
            l = line.split()
            if l[2] == "gene":
                id_symbol = l.index("gene_name")
                id_ensembl = l.index("gene_id")
                if "hgnc_id" in l:
                    tag_HGNC = l[l.index("hgnc_id")+1]
                else:
                    tag_HGNC = ""
                if "havana_gene" in l:
                    tag_Havana = l[l.index("havana_gene")+1]
                else:
                    tag_Havana = ""

                alias.append([l[id_symbol + 1], l[id_ensembl + 1],
                              tag_HGNC, tag_Havana])
                genes.append([l[0], l[3], l[4], l[id_symbol + 1], "0", l[6]])

res_alias = list(set(tuple(g) for g in alias))
res_genes = list(set(tuple(g) for g in genes))
# Save alias to MM10
alias_file = "mm10/alias_mouse.txt"
with open(alias_file, "w") as f:
    for g in res_alias:
        print("\t".join(g), file=f)

# Save alias to MM9
alias_file = "mm9/alias_mouse.txt"
with open(alias_file, "w") as f:
    for g in res_alias:
        print("\t".join(g), file=f)

# Save genes to MM10
genes_file = "mm10/genes_mm10.txt"
with open(genes_file, "w") as f:
    for g in res_genes:
        print("\t".join(g), file=f)


##########################################################################
###  MM9
###  -- Update mm9/genes_mm.bed
##########################################################################

genes = []
file_path = "mm9/gencode.vM1.annotation.gtf"
print("Parsing "+file_path)
with open(file_path) as gtf_mm9:
    for line in tqdm(gtf_mm9, total=get_num_lines(file_path)):
        if line.startswith("#"):
            pass
        else:
            line = line.replace('"', "")
            line = line.replace(';', "")
            l = line.split()
            if l[2] == "gene":
                id_symbol = l.index("gene_name")
                genes.append([l[0], l[3], l[4], l[id_symbol + 1], "0", l[6]])

res_genes = list(set(tuple(g) for g in genes))

# Save genes to MM9
genes_file = "mm9/genes_mm9.txt"
with open(genes_file, "w") as f:
    for g in res_genes:
        print("\t".join(g), file=f)

##########################################################################
###  Chromosome Size
##########################################################################

chrom_size = {
    "hg19": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
    "hg38": "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    "mm9": "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes",
    "mm10": "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
}

for genome in chrom_size.keys():
    with open(os.path.join(genome, "chrom.sizes." + genome), "w") as f:
        data = urlopen(chrom_size[genome])
        for line in data:
            l = line.decode('utf-8').strip()
            if "_" not in l:
                print(l, file=f)

