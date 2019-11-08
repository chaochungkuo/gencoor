from gencoor.util import GenomeConfig
from gencoor.coordinates import GenCoor, GenCoorSet
import numpy
from tqdm import tqdm

class Annotation:
    def __init__(self, str name, str genome):
        self.name = name
        self.header = ["chrom", "start", "end", "symbol", "score", "strand",
                       "annotation source", "feature type", "gene_type",
                       "gene_id", "hgnc_id", "havana_gene",
                       "transcript_name", "transcript_id", "transcript_type", "havana_transcript",
                       "exon_number", "exon_id", "protein_id"]
        self.table = []

        GConfig = GenomeConfig(genome)
        if GConfig:
            self.load(gtf_file_path=GConfig.get_gtf())
        else:
            self.load(gtf_file_path=genome) # Add error when no such file

    def __len__(self):
        return numpy.size(self.table,0)
        
    def load(self, str gtf_file_path):
        """
        0 chromosome name	chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M} or GRC accession a
        1 genomic start location	integer-value (1-based)
        2 genomic end location	integer-value
        3 genomic strand	{+,-}
        4 genomic phase (for CDS features)	{0,1,2,.}
        5 annotation source	{ENSEMBL,HAVANA}
        6 feature type	{gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
        7 gene_id	all	ENSGXXXXXXXXXXX.X b,c _Xg	all
        8 transcript_id d	all except gene	ENSTXXXXXXXXXXX.X b,c _Xg	all
        9 gene_type	all	list of biotypes	all
        10 gene_status e	all	{KNOWN, NOVEL, PUTATIVE}	until 25 and M11
        11 gene_name	all	string	all
        12 transcript_type d	all except gene	list of biotypes	all
        13 transcript_statusd,e	all except gene	{KNOWN, NOVEL, PUTATIVE}	until 25 and M11
        14 transcript_name d	all except gene	string	all
        15 exon_number f	all except gene/transcript/Selenocysteine	integer (exon position in the transcript from its 5' end)	all
        16 exon_id f	all except gene/transcript/Selenocysteine	ENSEXXXXXXXXXXX.X b _Xg

        """
        cdef str line
        cdef int num_lines
        cdef list l
        cdef str symbol
        cdef str ensembl
        cdef str gene_type
        cdef str hgnc_id
        cdef str havana_gene
        cdef str transcript_name
        cdef str transcript_id
        cdef str transcript_type
        cdef str havana_transcript
        cdef str exon_number
        cdef str exon_id
        cdef str protein_id
        cdef list new_row


        def check_label(str label, list l):
            if label in l:
                return l[l.index(label) + 1]
            else:
                return ""
        print("Loading annotation file: "+ gtf_file_path)
        num_lines = sum(1 for line in open(gtf_file_path, 'r'))
        with open(gtf_file_path) as f:
            for line in tqdm(f, total=num_lines):
                if not line.startswith("#"):
                    line = line.replace('"', "")
                    line = line.replace(';', "")
                    l = line.strip().split()
                    symbol = l[l.index("gene_name")+1]
                    ensembl = l[l.index("gene_id")+1]
                    gene_type = l[l.index("gene_type") + 1]
                    hgnc_id = check_label("hgnc_id", l)
                    havana_gene = check_label("havana_gene", l)
                    transcript_name = check_label("transcript_name", l)
                    transcript_id = check_label("transcript_id", l)
                    transcript_type = check_label("transcript_id", l)
                    havana_transcript = check_label("havana_transcript", l)
                    exon_number = check_label("exon_number", l)
                    exon_id = check_label("exon_id", l)
                    protein_id = check_label("protein_id", l)

                    new_row = [l[0], l[3], l[4], symbol, l[5], l[6], l[1], l[2], gene_type,
                               ensembl, hgnc_id, havana_gene,
                               transcript_name, transcript_id, transcript_type, havana_transcript,
                               exon_number, exon_id, protein_id]
                    self.table.append(new_row)
        self.table = numpy.array(self.table)

    def output_regions(self, table, str name, str gene_name="symbol"):
        cdef int gn_idx
        res = GenCoorSet(name)
        if gene_name == "symbol":
            gn_idx = 3
        elif gene_name == "gene_id":
            gn_idx = 9
        for r in table:
            res.add(GenCoor(chrom=r[0], start=int(r[1]), end=int(r[2]), name=r[gn_idx], strand=r[5]))
        return res

    def filter(self, dict query, str name):
        """

        :param query:
        :return:
        """
        res = self.boolean_index(query)
        return self.output_regions(res, name)

    def boolean_index(self, query):
        filter_critiria = []
        for key, values in query.items():
            if values is not list:
                values = [values]
            for v in values:
                filter_critiria.append(self.table[:, self.header.index(key)] == v)
        index = numpy.all(filter_critiria, axis=0)
        return self.table[index, :]

    def get_genes(self):
        q = {"feature type": "gene"}
        return self.filter(q, self.name+"_genes")

    def get_transcripts(self):
        q = {"feature type": "transcript"}
        return self.filter(q, self.name+"_transcript")

    def get_all_promoters(self, length=1000):
        promoters = self.get_transcripts()
        promoters.relocate(mode='5end as 3end', width=length, inplace=True)
        promoters.name = self.name+"_promoters"
        return promoters


