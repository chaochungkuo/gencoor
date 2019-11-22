#!/usr/bin/env python
"""
Usage:
    main.py bam2bw <bam_file_path> <bw_file_path>
    main.py diffpeak <config_file> <output_directory>

DiffPeak is used for detecting differential peaks from multiple ChIP-Seq data.

Options:

"""

from docopt import docopt
import configparser
import sys
from gencoor.signals import SignalProfile
from gencoor.coordinates import GenCoorSet

class DiffPeakConfig():
    def __init__(self, filepath):
        self.config = configparser.ConfigParser(allow_no_value=True)
        self.config.read(filepath)
        self.files_dict = self.output_files_dict()

    def output_files_dict(self):
        res = {}
        for lab in self.config["Condition 1"]:
            res[lab] = self.config["Condition 1"][lab]
        for lab in self.config["Condition 2"]:
            res[lab] = self.config["Condition 2"][lab]
        return res

    def genome(self):
        return self.config["Parameters"]["genome"]

    def get_inputs(self):
        res = {}
        if "Input" in self.config.sections():
            for lab in self.config["Input"]:
                res[lab] = self.config["Input"][lab]
        return res


if __name__ == '__main__':
    arg = docopt(__doc__)
    if arg["diffpeak"]:
        config = DiffPeakConfig(filepath=arg["<config_file>"])
        ref_back = GenCoorSet(name="background")
        ref_back.get_chromosomes(genome=config.genome())
        #
        # genset1 = GenCoorSet(name="Test_set")
        # genset1.load(
        #     filename="/projects/epi_aging_signature/exp/CTCF_ChIPSeq_analysis_UKA/3_Peaks_analysis/windows/66_CpGs.bed",
        #     filetype="BED")
        # genset1.relocate(mode="center as center", width=2000)
        bin = 1000000
        step = 500000
        sig = SignalProfile(regions=ref_back, genome=config.genome(), bin=bin, step=step)
        sig.load_files(file_dict=config.files_dict)

        # Normalization
        # sig.norm_bakcground(genome="hg38")

        # Input
        if "Input" in config.config.sections():
            sig2 = SignalProfile(regions=ref_back, genome=config.genome(), bin=bin, step=step)
            sig2.load_files(file_dict=config.files_dict)
            sig.minus_coverage(sig2.cov)
        sig.coverages2bigwigs(directory=arg["<output_directory>"])

        # Detecting diff. peaks



    # elif arg["diffpeak"]:

    print(arg)

