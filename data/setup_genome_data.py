import argparse

parser = argparse.ArgumentParser(description='Setup the genome data under the defined directory.')
parser.add_argument('-g','--genome', dest='genome', type=str, help='the given genome name (such as hg38, mm10...).')
parser.add_argument('-a', '--all', dest='all', action='store_true', default=False,
                    help='Download both the genome (FASTA) and annotation (GTF). (default: False)')
parser.add_argument('--fasta', dest='fasta', action='store_true', default=False,
                    help='Download the genome (FASTA). (default: False)')

args = parser.parse_args()
print(args.accumulate(args.integers))