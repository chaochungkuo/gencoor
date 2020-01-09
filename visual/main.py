#!/usr/bin/env python
"""
Usage:
    main.py lineplot <matrix_file_path> <output> [--col] [--row] [--c] [--average] [--scol] [--srow]

lineplot is used for generating lineplots between the given regions and the signals.

Options:

"""

from docopt import docopt
from gencoor.experiment_matirx import ExpMatrix
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import configparser
import sys
from gencoor.signals import SignalProfile
from gencoor.coordinates import GenCoorSet



if __name__ == '__main__':
    arg = docopt(__doc__)
    if arg["lineplot"]:
        exp = ExpMatrix()
        exp.read(path=arg['<matrix_file_path>'])
        # exp.print()

        # tags for sorting
        tags_rows = exp.get_all_tags(arg["--row"])
        tags_cols = exp.get_all_tags(arg["--col"])
        tags_colors = exp.get_all_tags(arg["--c"])

        # Plotting
        fig = plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(len(tags_rows), len(tags_cols),
                               width_ratios=[1, 1])

        for i, row in enumerate(tags_rows):
            for j, col in enumerate(tags_cols):

                plt.subplot(gs[i*len(tags_cols)+j])
                xyz = np.array(np.random.random((100, 3)))
                plt.scatter(xyz[:, 0], xyz[:, 1])

        plt.show()
        plt.savefig(arg["<output>"])

    # elif arg["diffpeak"]:

    print(arg)

