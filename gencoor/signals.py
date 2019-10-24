"""A class to calculate the signal over the given GenCoorSet."""
import pysam
import sys
import pyBigWig
from collections import defaultdict
from tqdm import tqdm
from gencoor.coordinates import GenCoorSet, GenCoor

class SignalProfile:
    def __init__(self, regions, bin=200, step=50):
        self.regions = regions
        self.bin = bin
        self.step = step
        self.cov = {}

    def bam_read_pair_generator(self, bam, total, region_string=None):
        """Generate read pairs in a BAM file or within a region string. Reads are added to read_dict until a pair is found.
        """
        print("Iterating all the paired reads...")
        read_dict = defaultdict(lambda: [None, None])
        fetch_reads = bam.fetch(region=region_string)
        pbar = tqdm(total=total)
        for read in fetch_reads:
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                else:
                    yield read_dict[qname][0], read
                del read_dict[qname]
            pbar.update(1)
        pbar.close()
        return read_dict

    def bam_total_reads(self, filename):
        idxstats= pysam.idxstats(filename).split("\n")
        idxstats = [l.split() for l in idxstats if l]
        return sum([int(l[2]) + int(l[3]) for l in idxstats])

    def bam_detect_fragment_size(self, filename):
        """Detect the average length of the fragment from the paired-end reads. (Only used for paired-end bamfiles)"""
        print("Detecting fragment size...", end="", flush=True)
        length_list = []
        total_reads = self.bam_total_reads(filename)
        bam = pysam.AlignmentFile(filename, "rb")
        for read1, read2 in self.bam_read_pair_generator(bam, total=total_reads):
            if read2.is_reverse:
                l = read2.reference_start + read2.infer_query_length() - read1.reference_start
            else:
                l = read1.reference_start + read1.infer_query_length() - read2.reference_start
            length_list.append(l)
        print(" OK")
        return sum(length_list)/len(length_list)

    def bam_get_paired_reads_positions(self, filename):
        """Return the position of left end and right end of the fragment (including read1 and read2). (Only used for paired-end bamfiles)"""
        print("Iterating the paired reads position...", end="", flush=True)
        pos_left = []
        pos_right = []
        total_reads = self.bam_total_reads(filename)
        bam = pysam.AlignmentFile(filename, "rb")
        for read1, read2 in self.bam_read_pair_generator(bam, total=total_reads):
            if read2.is_reverse:
                pos_left.append(read1.reference_start)
                pos_right.append(read2.reference_start + read2.infer_query_length())
            else:
                #XXXXXXXXXXXXXXXXXXXXXXx
                pos_left.append(read1.reference_start)
                pos_right.append(read2.reference_start + read2.infer_query_length())
                l = read1.reference_start + read1.infer_query_length() - read2.reference_start
            length_list.append(l)
        print(" OK")
        return sum(length_list) / len(length_list)

    def load_bam(self, filename, label, extension_size=-1):
        """Load a BAM and calcualte the coverage on the defined genomic coordinates. If extenstion is not defined, the length of the paired reads will be calculated. extension is only used for single-end sequencing."""
        self.cov[label] = {}
        # Detect fragment size
        # average_frag_size = self.bam_detect_fragment_size(filename)
        buffer_extesion = 500
        bamfile = pysam.AlignmentFile(filename, "rb")
        for r in self.regions:
            win1 = r.start
            win2 = r.start + self.bin
            self.cov[label][r] = []
            while win2 < r.end:
                reads_in_win = bamfile.fetch(r.chrom, win1-buffer_extesion, win2+buffer_extesion)
                for r in reads_in_win:
                    if not r.is_proper_pair or r.is_secondary or r.is_supplementary:
                        continue
                read_start = [r.reference_start for r in reads_in_win if r]
                c = sum([1 for read in ])
                self.cov[label][r].append(c)
                win1 += self.step
                win2 += self.step

    def load_bigwig(self, filename, label):
        self.cov[label] = {}
        bw = pyBigWig.open(filename)
        for r in self.regions:
            win1 = r.start
            win2 = r.start + self.bin
            self.cov[label][r] = []
            while win2 < r.end:
                c = bw.stats(r.chrom, win1, win2, type="mean")[0]
                if not c:
                    c = 0
                self.cov[label][r].append(c)
                win1 += self.step
                win2 += self.step

if __name__ == '__main__':
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=100000, end=150000, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr3_random", start=1, end=80, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test2", strand="."))
    signal_profile = SignalProfile(regions=genset1)
    signal_profile.load_bam(filename="/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/1_Alignment/bam_filtered/Y6.bam",
                            label="test")
    signal_profile.load_bigwig(
        filename="/media/work_disk/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/2_Peak_Calling/THOR_res_masked/CTCF_ChIPSeq-s1-rep1.bw",
        label="test_bw")
    print(signal_profile.cov)
