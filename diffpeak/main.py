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

    def bin_size(self):
        return int(self.config["Parameters"]["bin_size"])

    def step_size(self):
        return int(self.config["Parameters"]["step_size"])

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
        sig = SignalProfile(regions=ref_back, genome=config.genome(),
                            bin=config.bin_size(), step=config.step_size())
        sig.load_files(file_dict=config.files_dict)

        # Normalization
        # sig.norm_bakcground(genome="hg38")
        sig.norm_library_size()

        # Input
        if "Input" in config.config.sections():
            sig2 = SignalProfile(regions=ref_back, genome=config.genome(),
                                 bin=config.bin_size(), step=config.step_size())
            sig2.load_files(file_dict=config.get_inputs())
            sig.minus_coverage(sig2.cov)
        sig.coverages2bigwigs(directory=arg["<output_directory>"])

        # Detecting diff. peaks



    # elif arg["diffpeak"]:

    print(arg)

