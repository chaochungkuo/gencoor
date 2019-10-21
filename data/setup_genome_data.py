#!/usr/bin/env python
"""
Usage:
    setup_genome_data.py mm9 [--fasta] [--gtf]
    setup_genome_data.py mm10 [--fasta] [--gtf]
    setup_genome_data.py hg19 [--fasta] [--gtf]
    setup_genome_data.py hg38 [--fasta] [--gtf]

Download the latest genome FASTA file or annotation GTF according to the given genome.

Options:
    --fasta                  Download the genome FASTA file [default: False]
    --gtf                    Download the annotation GTF file [default: False]
"""
from docopt import docopt
import sys
import os
import gzip
import ftplib
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed
from tqdm import tqdm
import urllib.request
from contextlib import closing

latest_fasta = {
"mm9":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/NCBIM37.genome.fa.gz",
"mm10":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz",
"hg19":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz",
"hg38":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz"
}
latest_gtf = {
"mm9":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz",
"mm10":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz",
"hg19":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
"hg38":"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"
}


class TqdmUpTo(tqdm):
    """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""

    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)  # will also set self.n = b * bsize


def download(url, target_dir):
    base_name = os.path.basename(url).replace(".gz", "")
    tmp_file = os.path.join(target_dir, base_name+".gz")
    goal_file = os.path.join(target_dir, base_name)

    eg_link = url
    with TqdmUpTo(unit='B', unit_scale=True, miniters=1,
                  desc=tmp_file) as t:  # all optional kwargs
        urllib.request.urlretrieve(url, filename=tmp_file,
                           reporthook=t.update_to, data=None)
    # Unzip the file
    gzfile = gzip.GzipFile(tmp_file, 'rb')
    s = gzfile.read()
    gzfile.close()
    output = open(goal_file, 'wb')
    output.write(s)
    output.close()
    os.remove(tmp_file)


if __name__ == '__main__':
    arg = docopt(__doc__)
    genome = ""
    for gg in latest_fasta.keys():
        if arg[gg]:
            genome = gg

    if arg['--fasta']:
        download(latest_fasta[genome], genome)
    if arg["--gtf"]:
        download(latest_gtf[genome], genome)
