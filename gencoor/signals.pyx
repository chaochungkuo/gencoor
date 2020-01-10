"""A class to calculate the signal over the given GenCoorSet."""
import pysam
import os
import pyBigWig
from collections import defaultdict, OrderedDict
from tqdm import tqdm
from gencoor.util import GenomeConfig
from gencoor.coordinates import GenCoorSet, GenCoor
import pprint
from natsort import natsorted
import sys
import numpy as np
import re

class SignalProfile:
    def __init__(self, regions, str genome, int bin=200, int step=100):
        self.regions = regions
        self.regions.sort()
        self.genome = genome
        self.bin = bin
        self.step = step
        self.cov = OrderedDict()
        self.file_path = OrderedDict()
        self.scaling_factor = OrderedDict()

    def load_files(self, file_dict, disable_progressbar=False):
        """Load each file according to the given labels and paths"""
        cdef str p
        cdef int i
        cdef str l
        for i, l in enumerate(list(file_dict.keys())):
            p = file_dict[l]
            vlist = list(self.file_path.values())
            if p in vlist:
                klist = list(self.file_path.keys())
                self.file_path[l] = p
                self.cov[l] = self.cov[klist[vlist.index(p)]]
            else:
                if p.lower().endswith(".bam"):
                    self.load_bam(filename=p, label=l, disable_progressbar=disable_progressbar)
                elif p.lower().endswith(".bigwig") or p.lower().endswith(".bw"):
                    self.load_bigwig(filename=p, label=l, disable_progressbar=disable_progressbar)
                elif p.lower().endswith(".bedgraph"):
                    self.load_bedgraph(filename=p, label=l, disable_progressbar=disable_progressbar)

    def get_chrom_size_tuples(self):
        cdef str line
        genome = GenomeConfig(genome=self.genome)
        res = []
        with open(genome.get_chromosome_sizes()) as f:
            for line in f:
                l = line.strip().split()
                res.append((l[0], int(l[1])))
        res = natsorted(res, key=lambda x: x[0])
        return res

    def bam_read_pair_generator(self, bam, str chrom, int start, int end):
        """Generate read pairs in a BAM file or within a region string. Reads are added to read_dict until a pair is found.
        """
        cdef str qname
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

    def bam_total_reads(self, str filename):
        idxstats= pysam.idxstats(filename).split("\n")
        idxstats = [l.split() for l in idxstats if l]
        return sum([int(l[2]) + int(l[3]) for l in idxstats])

    def bam_detect_fragment_size(self, str filename):
        """Detect the average length of the fragment from the paired-end reads. (Only used for paired-end bamfiles)"""
        print("Detecting fragment size...", end="", flush=True)
        cdef int l
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

    def bam_count_paired_reads(self, bam, str chrom, int start, int end):
        """Return the position of left end and right end of the fragment (including read1 and read2). (Only used for paired-end bamfiles)"""
        cdef int c = 0
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

    def load_bam(self, str filename, str label, disable_progressbar=False, str progressbar_mode="region"):
        """Load a BAM and calcualte the coverage on the defined genomic coordinates. If extenstion is not defined, the length of the paired reads will be calculated. extension is only used for single-end sequencing."""
        print("Loading BAM file ..." + filename[-30:])

        cdef int win1
        cdef int win2
        cdef int c

        self.file_path[label] = filename
        self.cov[label] = OrderedDict()
        # Detect fragment size
        # average_frag_size = self.bam_detect_fragment_size(filename)
        # buffer_extesion = 300
        bam = pysam.AlignmentFile(filename, "rb")
        if progressbar_mode=="region":
            pbar = tqdm(total=len(self.regions), disable=disable_progressbar)
        for r in self.regions:
            # try:
            if progressbar_mode=="window":
                print("Loading BAM file ..." + filename[-30:] + "\t" + r.chrom)
                pbar = tqdm(total=len(r)//self.step+1, disable=disable_progressbar)
            win1 = r.start
            win2 = r.start + self.bin
            self.cov[label][str(r)] = []
            cont_loop = True
            while cont_loop:
                start = win1
                # start = win1 - self.step
                # if start < 0:
                #     start = 0
                # print([start, win2])
                c = self.bam_count_paired_reads(bam, r.chrom, start, win2)
                # print(c)
                self.cov[label][str(r)].append(c)
                win1 += self.step
                win2 += self.step
                if win1 > r.end:
                    cont_loop = False
                if progressbar_mode=="window":
                    pbar.update(1)
            # except:
            #     pass
            if progressbar_mode=="region":
                pbar.update(1)
            elif progressbar_mode=="window":
                pbar.close()

    def load_bigwig(self, str filename, str label, disable_progressbar=False):
        print("Loading BigWig file ..." + filename[-30:])

        cdef int win1
        cdef int win2
        cdef int c

        self.file_path[label] = filename
        self.cov[label] = OrderedDict()
        bw = pyBigWig.open(filename)
        pbar = tqdm(total=len(self.regions), disable=disable_progressbar)
        for r in self.regions:
            try:
                win1 = r.start
                win2 = r.start + self.bin
                self.cov[label][str(r)] = []
                while win2 < r.end:
                    c = bw.stats(r.chrom, win1, win2, type="mean")[0]
                    # print(c)
                    if not c:
                        c = 0
                    self.cov[label][str(r)].append(c)
                    win1 += self.step
                    win2 += self.step

            except:
                pass
            pbar.update(1)
        pbar.close()

    def load_bedgraph(self, str filename, str label, disable_progressbar=False):
        print("Loading BedGraph file ..." + filename[-30:])
        cdef int i
        cdef int j
        cdef int win1
        cdef int win2
        cdef float c
        cdef int ind_step

        self.file_path[label] = filename
        self.cov[label] = OrderedDict()
        bg = GenCoorSet(name=label)
        bg.load(filename=filename, filetype="BedGraph")
        bg.sort()
        chroms, starts, ends, names, scores, strands, datas = bg.to_lists()
        i = 0
        j = 0
        loop_continue = True
        pbar = tqdm(total=len(self.regions), disable=disable_progressbar)

        while loop_continue:
            r = self.regions[j]
            if GenCoor.overlap_pairing(r.chrom, r.start,
                                       r.end, r.strand,
                                       chroms[i], starts[i], ends[i], strands[i]):
                win1 = r.start
                win2 = r.start + self.bin
                init_cov = False
                ind_step = 0
                if r not in self.cov[label].keys():
                    init_cov = True
                    self.cov[label][str(r)] = []
                while win2 < r.end:
                    if GenCoor.overlap_pairing(r.chrom, win1, win2, r.strand,
                                               chroms[i], starts[i], ends[i], strands[i]):
                        c = float(scores[i])
                    else:
                        c = float(0)
                    if init_cov:
                        self.cov[label][str(r)].append(c)
                    else:
                        self.cov[label][str(r)][ind_step] += c
                    win1 += self.step
                    win2 += self.step
                    ind_step += 1
                i += 1
            elif r > bg[i]:
                i += 1
            elif r < bg[i]:
                # self.cov[label][str(r)] = [0] * (len(r)//self.step)
                j += 1
                pbar.update(1)
            if i == len(bg):
                loop_continue = False
            if j == len(self.regions):
                loop_continue = False

        pbar.close()

    def normalize_by_scaling_factors(self, factors):
        """Rescale the coverage by the given scaling facotrs. The input factors is a dictionary with label as the key and scaling factor as the value."""
        cdef str k
        cdef float v
        cdef float j

        for k, v in factors.items():
            for r in self.regions:
                self.cov[k][str(r)] = [j * float(v) for j in self.cov[k][str(r)]]

    def norm_bakcground(self, lower_bound=20, upper_bound=80, bin=10000, genome=""):
        """Normalize the coverage by removing the bins with values between the given lower bound and upper bound in percentage."""
        # Getting background
        # if genome:
        #     ref_back = GenCoorSet(name="background")
        #     ref_back.get_chromosomes(genome=genome)
        #     ref_back.filter_chromosome(chrom="chr22")
        #     # ref_back.standard_chromosome()
        #     # print(len(ref_back))
        #     ref_back.segmentize(width=bin, inplace=True)
        #     # print(len(ref_back))
        #     sig = SignalProfile(regions=ref_back, bin=bin, step=bin)
        #     sig.load_files(labels=self.file_path.keys(), paths=list(self.file_path.values()),
        #                    disable_progressbar=False)
        # else:
        # sig = SignalProfile(regions=self.regions, bin=self.bin, step=self.step)
        # sig.load_files(labels=self.file_path.keys(), paths=list(self.file_path.values()),
        #                disable_progressbar=False)

        # Get array
        cov_array = self.cov2array()
        # print(cov_array)
        for i, label in enumerate(self.cov.keys()):
            cov = cov_array[i,:]
            # print(cov)
            low_b = np.percentile(cov[cov > 0], lower_bound)
            upp_b = np.percentile(cov, upper_bound)
            # print([low_b, upp_b])
            self.scaling_factor[label] = np.sum(cov[(cov > low_b) & (cov < upp_b)])

        # Get scaling factors
        base = min(list(self.scaling_factor.values()))
        for k,v in self.scaling_factor.items():
            self.scaling_factor[k] = base/v
        print(self.scaling_factor)
        self.normalize_by_scaling_factors(self.scaling_factor)

    def norm_library_size(self):
        """Normalize the coverage by simply the library size (total number of reads) of each sample. It is only applicable when all the input files are BAM files."""
        cdef int i
        cdef str label
        cdef float min_count

        read_counts = {}
        for i, label in enumerate(self.cov.keys()):
            read_counts[label] = self.bam_total_reads(self.file_path[label])
        min_count = float(min(list(read_counts.values())))
        for i, label in enumerate(self.cov.keys()):
            self.scaling_factor[label] = min_count/read_counts[label]
        print(self.scaling_factor)
        self.normalize_by_scaling_factors(self.scaling_factor)



    def cov2array(self):
        """Return a dictionary with labels as the keys and the arrays of the coverage as the values."""
        cdef str lab
        cdef str k
        res = {}
        for lab, d in self.cov.items():
            # print(d)
            res[lab] = []
            for k, r in d.items():
                res[lab].append(r)
            res[lab] = np.array(res[lab])
            # r = []
            # for v in list(d.values()):
            #     r += v
            # res.append(r)
            # print(res[lab])
        return res

    def cov2bigwig(self, cov, filename):
        """Generate BigWig file from coverage profile"""
        cdef int i
        cdef float s
        # Get Bedgraph
        bg = [[], [], [], []]
        for r in cov.keys():
            l = re.split(":|-| ", r)
            # print(l)
            for i, s in enumerate(cov[r]):
                bg[0].append(l[0]) # chrom
                bg[1].append(int(l[1]) + i * self.step) # start
                bg[2].append(int(l[1]) + (i+1) * self.step) # end
                bg[3].append(float(s)) # score
        bw = pyBigWig.open(filename, "w")
        bw.addHeader(self.get_chrom_size_tuples())
        bw.addEntries(chroms=bg[0], starts=bg[1], ends=bg[2], values=bg[3])
        bw.close()

    def cov2bedgraph(self, cov, filename):
        """Generate BedGraph file from coverage profile"""
        # Get Bedgraph
        with open(filename, "w") as f:
            for r in cov.keys():
                l = re.split(":|-| ", r)
                for i, s in enumerate(cov[r]):
                    print("\t".join([l[0],
                                     str(int(l[1]) + i * self.step),
                                     str(int(l[1]) + (i+1) * self.step),
                                     str(s)]), file=f)

    def coverages2bigwigs(self, directory):
        """Generate BigWig files for each sample"""
        if not os.path.exists(directory):
            os.makedirs(directory)
        for lab in self.cov.keys():
            print(lab)
            self.cov2bedgraph(cov=self.cov[lab], filename=os.path.join(directory, lab+".bedgraph"))
            self.cov2bigwig(cov=self.cov[lab], filename=os.path.join(directory, lab+".bw"))

    def minus_coverage(self, cov_dict):
        """Reduce the coverage by the given coverage dictionary which has the same labels as its keys."""
        for lab in self.cov.keys():
            for r in self.cov[lab].keys():
                self.cov[lab][r] = [x1 - x2 for (x1, x2) in zip(self.cov[lab][r], cov_dict[lab][r])]
                self.cov[lab][r] = [0 if i < 0 else i for i in self.cov[lab][r]]

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
    genset1.relocate(mode="5end as 3end", width=2000)
    genset1.list = genset1.list[0:100]
    signal_profile = SignalProfile(regions=genset1, genome="hg38", bin=200, step=50)
    signal_profile.load_bam(filename="/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/1_Alignment/bam_filtered/Y6.bam",
                            label="bam")
    # signal_profile.load_bigwig(
    #     filename="/media/work_disk/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/2_Peak_Calling/THOR_res_masked/CTCF_ChIPSeq-s1-rep1.bw",
    #     label="bw")

    # signal_profile.coverages2bigwigs(directory="./")
    # genset1 = GenCoorSet(name="Test_set")
    # genset1.load(filename="/media/work_disk/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/3_Peaks_analysis/bedgraph/GSE40279/GSE40279_pearson_hyper250.bed", filetype="BED")
    # genset1.relocate(mode="center as center", width=2000)
    # signal_profile = SignalProfile(regions=genset1)
    # signal_profile.load_bedgraph(
    #     filename="/media/work_disk/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/3_Peaks_analysis/bedgraph/GSE40279/GSE40279_pearson.bedgraph",
    #     label="bedgraph")
    # sc = {"bam": 10, "bw": 0.5}
    # signal_profile.normalize_by_scaling_factors(factors=sc)
    signal_profile.coverages2bigwigs(directory="/home/joseph/source_code/gencoor/test/diffpeaks")
    # pp = pprint.PrettyPrinter(indent=4)

    # pp.pprint(signal_profile.cov["bam"])
    # pp.pprint(signal_profile.cov["bw"])
    # pp.pprint(signal_profile.cov["bedgraph"])
    # print(signal_profile.cov["bedgraph"])
    # for k, v in signal_profile.cov["bedgraph"].items():
    #     if sum(v) != 0:
    #         print(v)
