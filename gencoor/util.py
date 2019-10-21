import configparser
import os
from .exceptions import GenomeUndefinedError

data_path = os.path.expanduser(os.getenv("GENCOORDATA", os.path.join(os.getenv("HOME"), "gencoor_data")))

class GenomeConfig():
    def __init__(self, genome):
        self.genome = genome
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(data_path, 'data.config'))
        self.config.read(os.path.join(data_path, 'data.config.user'))

        if genome not in self.config:
            print(f"The input genome ({genome}) is not defined in the configure file \
            ({os.path.join(data_path, 'data.config.user')})")
            raise GenomeUndefinedError

    def get_fasta(self):
        return self.config[self.genome]["genome"]

    def get_gtf(self):
        return self.config[self.genome]["annotation"]

    def get_chromosome_sizes(self):
        return self.config[self.genome]["chromosome_sizes"]

    def get_gene_alias(self):
        return self.config[self.genome]["gene_alias"]

    def get_genes(self):
        return self.config[self.genome]["genes"]


