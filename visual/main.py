#!/usr/bin/env python
"""
Usage:
    main.py lineplot <matrix_file_path> <output_directory> [options]

lineplot is used for generating lineplots between the given regions and the signals.

Options:
--col <col-header>  Define the header in the matrix for distinct columns
--row <row-header>  Define the header in the matrix for distinct rows
--c <color-header>  Define the header in the matrix for distinct colors
--scol              Sync the scale according to the columns.
--srow              Sync the scale according to the rows.
--test              Take only 10 regions for each file for quick testing.

"""

from docopt import docopt
from gencoor.experiment_matirx import ExpMatrix
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec
import matplotlib.cm as cmx
import os
from scipy.interpolate import make_interp_spline, BSpline
import seaborn as sns
import configparser
import sys
from gencoor.signals import SignalProfile
from gencoor.coordinates import GenCoorSet



if __name__ == '__main__':
    arg = docopt(__doc__)
    print(arg)
    if arg["lineplot"]:

        ext = 1000
        exp = ExpMatrix()
        exp.read(path=arg['<matrix_file_path>'])
        # exp.print()

        # tags for sorting
        tags_rows = exp.get_all_tags(arg["--row"])
        tags_cols = exp.get_all_tags(arg["--col"])
        tags_colors = exp.get_all_tags(arg["--c"])

        # cmap = plt.cm.get_cmap('Spectral')
        # colors = cmap(np.arange(len(tags_colors)))
        # cmn = np.linspace(0, 256, len(tags_colors))
        # colors = [cmap(c) for c in cmn]
        cmap = plt.get_cmap('jet')
        cNorm = colors.Normalize(vmin=0, vmax=1)
        colors = [cmap(i) for i in np.linspace(0, 1, len(tags_colors))]
        #
        # jet = cm = plt.get_cmap('jet')
        # cNorm = colors.Normalize(vmin=0, vmax=1)
        # scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        print(colors)
        print(tags_rows)
        print(tags_cols)
        print(tags_colors)

        # Plotting
        fig = plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(len(tags_rows), len(tags_cols))
        # print(exp.tags)
        for i, row in enumerate(tags_rows):
            for j, col in enumerate(tags_cols):
                plt.subplot(gs[i * len(tags_cols) + j])
                # labels
                if j == 0:
                    plt.ylabel(row)
                if i == 0:
                    plt.title(col)

                for k, color in enumerate(tags_colors):
                    beds = exp.filter_by_tags([row, col, color, "regions"])
                    signals = exp.filter_by_tags([row, col, color, "signals"])
                    # print([row, col, color] + beds + signals)
                    for bed in beds:
                        regions = GenCoorSet(name=bed)
                        regions.load(filename=exp.get_file(bed))
                        if arg["--test"]:
                            regions.list = regions.list[0:50]
                        regions.relocate(mode='center as center', width=ext*2)
                        sig = SignalProfile(regions, genome="hg19", bin=200, step=100)
                        for signal in signals:

                            sig.load_bigwig(filename=exp.get_file(signal),
                                            label=signal,
                                            disable_progressbar=True)
                        arrs = sig.cov2array()
                        # print(arrs)
                        for s in signals:
                            mean_a = np.mean(np.array(arrs[s]), axis=0)
                            x = np.linspace(-ext, ext, num=len(mean_a))
                            y = mean_a
                            x2 = np.linspace(-ext, ext, 2*len(mean_a))
                            spl = make_interp_spline(x, y, k=3)
                            y2 = spl(x2)
                            plt.plot(x2, y2, marker='', color=colors[k], linewidth=2, alpha=0.4,
                                     label=color)

                        # print(sig.cov.keys())


                # xyz = np.array(np.random.random((100, 3)))
                #
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        fig.savefig(os.path.join(arg["<output_directory>"], "output.pdf"),
                    bbox_inches='tight')

    # elif arg["diffpeak"]:



