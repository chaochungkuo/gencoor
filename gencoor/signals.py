"""A class to calculate the signal over the given GenCoorSet."""
import pysam
import pyBigWig

from gencoor.coordinates import GenCoorSet, GenCoor

class SignalProfile:
    def __init__(self, regions, bin=200, step=50):
        self.regions = regions
        self.bin = bin
        self.step = step
        self.cov = {}

    def load_bam(self, filename, label):
        self.cov[label] = {}
        samfile = pysam.AlignmentFile(filename, "rb")
        for r in self.regions:
            win1 = r.start
            win2 = r.start + self.bin
            self.cov[label][str(r)] = []
            while win2 < r.end:
                c = sum([1 for read in samfile.fetch(r.chrom, win1, win2)])
                self.cov[label][str(r)].append(c)
                win1 += self.step
                win2 += self.step

    def load_bigwig(self, filename, label, strand_specific=False):
        self.cov[label] = {}
        bw = pyBigWig.open(filename)
        for r in self.regions:
            win1 = r.start
            win2 = r.start + self.bin
            self.cov[label][str(r)] = []
            while win2 < r.end:
                c = bw.stats(r.chrom, win1, win2, type="mean")
                self.cov[label][str(r)].append(c)
                win1 += self.step
                win2 += self.step

if __name__ == '__main__':
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr3_random", start=1, end=80, name="test", strand="."))
    # genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test2", strand="."))
    signal_profile = SignalProfile(regions=genset1)
    signal_profile.load_bam(filename=)
