"""A class to calculate the signal over the given GenCoorSet."""
import pysam
import os
import pyBigWig
from collections import defaultdict
from tqdm import tqdm
from gencoor.coordinates import GenCoorSet, GenCoor
import pprint

class SignalProfile:
    def __init__(self, regions, bin=200, step=100):
        self.regions = regions
        self.bin = bin
        self.step = step
        self.cov = {}
        self.file_path = {}

    def bam_read_pair_generator(self, bam, chrom, start, end):
        """Generate read pairs in a BAM file or within a region string. Reads are added to read_dict until a pair is found.
        """
        read_dict = defaultdict(lambda: [None, None])
        fetch_reads = bam.fetch(chrom, start, end)
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
        return sum(length_list)/len(length_list)

    def bam_count_paired_reads(self, bam, chrom, start, end):
        """Return the position of left end and right end of the fragment (including read1 and read2). (Only used for paired-end bamfiles)"""
        c = 0
        for read1, read2 in self.bam_read_pair_generator(bam, chrom, start, end):
            c += 1
        # print(c)
        return c

    # def bam_get_paired_reads_positions(self, bam, chrom, start, end):
    #     """Return the position of left end and right end of the fragment (including read1 and read2). (Only used for paired-end bamfiles)"""
    #     pos_left = []
    #     pos_right = []
    #     for read1, read2 in self.bam_read_pair_generator(bam, chrom, start, end):
    #         if read2.is_reverse:
    #             pos_left.append(read1.reference_start)
    #             pos_right.append(read2.reference_start + read2.infer_query_length())
    #         else:
    #             pos_left.append(read2.reference_start )
    #             pos_right.append(read1.reference_start + read1.infer_query_length())
    #     return pos_left, pos_right

    def load_bam(self, filename, label):
        """Load a BAM and calcualte the coverage on the defined genomic coordinates. If extenstion is not defined, the length of the paired reads will be calculated. extension is only used for single-end sequencing."""
        print("Loading BAM file " + filename)
        self.file_path[label] = filename
        self.cov[label] = {}
        # Detect fragment size
        # average_frag_size = self.bam_detect_fragment_size(filename)
        # buffer_extesion = 300
        bam = pysam.AlignmentFile(filename, "rb")
        pbar = tqdm(total=len(self.regions))
        for r in self.regions:
            try:
                win1 = r.start
                win2 = r.start + self.bin
                self.cov[label][r] = []
                while win2 < r.end:
                    c = self.bam_count_paired_reads(bam, r.chrom, win1-self.step, win2+self.step)
                    self.cov[label][r].append(c)
                    win1 += self.step
                    win2 += self.step
            except:
                pass
            pbar.update(1)
        pbar.close()

    def load_bigwig(self, filename, label):
        print("Loading BigWig file " + filename)
        self.file_path[label] = filename
        self.cov[label] = {}
        bw = pyBigWig.open(filename)
        pbar = tqdm(total=len(self.regions))
        for r in self.regions:
            try:
                win1 = r.start
                win2 = r.start + self.bin
                self.cov[label][r] = []
                while win2 < r.end:
                    c = bw.stats(r.chrom, win1, win2, type="mean")[0]
                    # print(c)
                    if not c:
                        c = 0
                    self.cov[label][r].append(c)
                    win1 += self.step
                    win2 += self.step

            except:
                pass
            pbar.update(1)
        pbar.close()

    def normalize_by_scaling_factors(self, factors):
        """Rescale the coverage by the given scaling facotrs. The input factors is a dictionary with label as the key and scaling factor as the value."""
        for k, v in factors.items():
            for r in self.regions:
                self.cov[k][r] = [j * float(v) for j in self.cov[k][r]]

    # def
if __name__ == '__main__':
    # genset1 = GenCoorSet(name="Test_set")
    # genset1.add(GenCoor(chrom="chr1", start=100000, end=150000, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr3_random", start=1, end=80, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test2", strand="."))
    genset1 = GenCoorSet(name="Test_set")
    bedfile = os.path.join(os.getenv("HOME"), "gencoor_data/hg38/genes_hg38.bed")
    genset1.load(filename=bedfile, filetype="BED")
    genset1.relocate(mode="5end as 3end", width=1000)
    signal_profile = SignalProfile(regions=genset1[1000:1002])
    signal_profile.load_bam(filename="/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/1_Alignment/bam_filtered/Y6.bam",
                            label="bam")
    signal_profile.load_bigwig(
        filename="/media/work_disk/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/2_Peak_Calling/THOR_res_masked/CTCF_ChIPSeq-s1-rep1.bw",
        label="bw")
    sc = {"bam": 10, "bw": 0.5}
    signal_profile.normalize_by_scaling_factors(factors=sc)

    pp = pprint.PrettyPrinter(indent=4)

    pp.pprint(signal_profile.cov["bam"])
    pp.pprint(signal_profile.cov["bw"])
    # print(signal_profile.cov["test_bw"].values()[50:60])
