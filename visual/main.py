#!/usr/bin/env python
"""
Usage:
    main.py lineplot <matrix_file_path> <output_file> [options]
    main.py heatmap <matrix_file_path> <output_file> [options]

lineplot is used for generating lineplots between the given regions and the signals.

Options:
--col <col-header>  Define the header in the matrix for distinct columns
--row <row-header>  Define the header in the matrix for distinct rows
--c <color-header>  Define the header in the matrix for distinct colors
--srow              Sync scale of y-axis across the rows.
--scol              Sync scale of y-axis across the columns.
--average           Show the average of the lines under same grouping.
--test              Take only 10 regions for each file for quick testing.
--genome <genome>   Define the genome, e.g. mm9, hg19... etc.
--bin <bin>         Define the bin size [default: 200]
--step <step>       Define the step size [default: 100]
--ext <ext>         Define the width of the window in the view (bp) [default: 2000]
--cores <cores>     Define the cores for parallelprocessing [default: 1]


"""

from docopt import docopt
from gencoor.experiment_matirx import ExpMatrix
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from matplotlib import gridspec
# import matplotlib.cm as cmx
import os
from scipy.interpolate import make_interp_spline
    # , BSpline
# import seaborn as sns
# import configparser
# import sys
from gencoor.signals import SignalProfile
from gencoor.coordinates import GenCoorSet


def reshape_axes(axes, n_row, n_col):
    try:
        res = axes.reshape((n_row, n_col))
    except:
        res = np.array([[axes]])
    return res


def share_y(arg):
    if arg["--srow"] and arg["--scol"]:
        sharey = "all"
    elif arg["--srow"] and not arg["--scol"]:
        sharey = "row"
    elif not arg["--srow"] and arg["--scol"]:
        sharey = "col"
    else:
        sharey = "none"
    return sharey

def show_uniq_legend(ax, handles, labels):
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              loc='center left', bbox_to_anchor=(1, 0.5))

def smooth_plot(mean_array, extend):
    x = np.linspace(-extend, extend, num=len(mean_array))
    x2 = np.linspace(-extend, extend, 2 * len(mean_array))
    spl = make_interp_spline(x, mean_array, k=3)
    y2 = spl(x2)
    return x2, y2

def color_limit(tags_colors):
    if len(tags_colors) < 10:
        res = plt.get_cmap("Set1").colors
    else:
        cmap = plt.get_cmap('jet')
        res = [cmap(i) for i in np.linspace(0, 1, len(tags_colors))]
    return res

def bed_sig2arr(bedname, sig_names, exp, arg):
    regions = GenCoorSet(name=bedname)
    regions.load(filename=exp.get_file(bedname))
    if arg["--test"]:
        regions.list = regions.list[0:50]
    regions.relocate(mode='center as center', width=2*int(arg["--ext"]))
    sig = SignalProfile(regions, genome=arg["--genome"],
                        bin=int(arg["--bin"]), step=int(arg["--step"]),
                        cores=int(arg["--cores"]))
    for signal in sig_names:
        if exp.get_file(signal).endswith(".bw") or exp.get_file(signal).endswith(".bigwig"):
            sig.load_bigwig(filename=exp.get_file(signal), label=signal,
                            disable_progressbar=False, verbal=False)
        elif exp.get_file(signal).endswith(".bam"):
            sig.load_bam(filename=exp.get_file(signal), label=signal,
                            disable_progressbar=False, verbal=False)
    res = sig.cov2array()

    return res


def set_yaxis(n_row, n_col, axes, arg):
    ymax = []
    ymin = []
    for i in range(n_row):
        ymax.append([])
        ymin.append([])
        for j in range(n_col):
            ymax[i].append([])
            ymin[i].append([])
            bottom, top = axes[i, j].get_ylim()
            ymax[i][j] = top
            ymin[i][j] = bottom
    ymax = np.array(ymax)
    ymin = np.array(ymin)
    if arg["--srow"] and arg["--scol"]:
        ymax = np.amax(ymax)
        ymin = np.amin(ymin)
        for i in range(n_row):
            for j in range(n_col):
                axes[i, j].set_ylim(bottom=ymin, top=ymax)
    elif arg["--srow"] and not arg["--scol"]:
        ymax = np.amax(ymax, axis=1)
        ymin = np.amin(ymin, axis=1)
        for i in range(n_row):
            for j in range(n_col):
                axes[i, j].set_ylim(bottom=ymin[i], top=ymax[i])
    elif not arg["--srow"] and arg["--scol"]:
        ymax = np.amax(ymax, axis=0)
        ymin = np.amin(ymin, axis=0)
        for i in range(n_row):
            for j in range(n_col):
                axes[i, j].set_ylim(bottom=ymin[j], top=ymax[j])

    # Set the panel aspect ratio
    for i in range(n_row):
        for j in range(n_col):
            axes[i, j].set_aspect(1. / axes[i, j].get_data_ratio())


if __name__ == '__main__':
    arg = docopt(__doc__)
    print(arg)
    linewidth = 2
    alpha = 0.8

    if arg["lineplot"]:

        exp = ExpMatrix()
        exp.read(path=arg['<matrix_file_path>'])

        # tags for sorting
        tags_rows = exp.get_all_tags(arg["--row"])
        tags_cols = exp.get_all_tags(arg["--col"])
        tags_colors = exp.get_all_tags(arg["--c"])
        # print(tags_rows)
        # print(tags_cols)
        # print(tags_colors)

        colors = color_limit(tags_colors)

        # Plotting
        n_row = len(tags_rows)
        n_col = len(tags_cols)

        fig, axes = plt.subplots(ncols=n_col, nrows=n_row, constrained_layout=True,
                                 figsize=(n_col*3+2, n_row*3), sharex="all")
        axes = reshape_axes(axes, n_row, n_col)
        legend_handles = []
        legend_labels = []

        for i, row in enumerate(tags_rows):
            for j, col in enumerate(tags_cols):
                # print([row, col])
                # labels
                if j == 0:
                    axes[i, j].set_ylabel(row)
                if i == 0:
                    axes[i, j].set_title(col)

                for k, color in enumerate(tags_colors):
                    beds = exp.filter_by_tags([row, col, color, "regions"])
                    signals = exp.filter_by_tags([row, col, color, "signals"])
                    print([row, col, color] + beds + signals)

                    if len(beds) == 0 or len(signals) == 0:
                        continue

                    if arg["--average"]:
                        sum_coverage = []
                    for bed in beds:

                        arrs = bed_sig2arr(bed, signals, exp, arg)

                        for s in signals:
                            # print([row, col, color, bed, s])
                            mean_array = np.mean(np.array(arrs[s]), axis=0)
                            x, y = smooth_plot(mean_array, 2 * float(arg["--ext"]))
                            if not arg["--average"]:
                                axes[i, j].plot(x, y, marker='', color=colors[k],
                                                linewidth=linewidth, alpha=alpha, label=color)
                            else:
                                sum_coverage.append(y)
                    if arg["--average"]:
                        y_average = np.mean(np.array(sum_coverage), axis=0)
                        axes[i, j].plot(x, y_average, marker='', color=colors[k],
                                        linewidth=linewidth, alpha=alpha, label=color)

                    handles, labels = axes[i, j].get_legend_handles_labels()
                    legend_handles += handles
                    legend_labels += labels

                if i == 0 and j == n_col - 1:
                    show_uniq_legend(axes[i, j], legend_handles, legend_labels)
                else:
                    axes[i, j].legend().set_visible(False)
        set_yaxis(n_row, n_col, axes, arg)
        fig.savefig(arg["<output_file>"], bbox_inches='tight')

    elif arg["heatmap"]:

        exp = ExpMatrix()
        exp.read(path=arg['<matrix_file_path>'])

        # tags for sorting
        tags_rows = exp.get_all_tags(arg["--row"])
        tags_cols = exp.get_all_tags(arg["--col"])

        # Plotting
        n_row = len(tags_rows)
        n_col = len(tags_cols)

        fig, axes = plt.subplots(ncols=n_col, nrows=n_row, constrained_layout=False,
                                 figsize=(n_col * 3 + 2, n_row * 3), sharex="all")
        axes = reshape_axes(axes, n_row, n_col)

        for i, row in enumerate(tags_rows):
            for j, col in enumerate(tags_cols):



                beds = exp.filter_by_tags([row, col, "", "regions"])
                if len(beds) > 1:
                    print("There are more than one BED files sharing the same tags, only the first one will be used.")
                signals = exp.filter_by_tags([row, col, "", "signals"])
                if len(signals) > 1:
                    print("There are more than one BED files sharing the same tags, only the first one will be used.")

                arrs = bed_sig2arr(beds[0], [signals[0]], exp, arg)
                a = arrs[signals[0]]
                hm = axes[i, j].imshow(a, cmap='hot', interpolation='None')
                # Y ticks
                axes[i, j].get_yaxis().set_ticks([])
                # X ticks
                x_label_list = ['-'+arg["--ext"], '0', arg["--ext"]]
                xmin, xmax = axes[i, j].get_xlim()
                axes[i, j].set_xticks([xmin, int(0.5 * (xmax - xmin)), xmax])
                axes[i, j].set_xticklabels(x_label_list)

                # labels
                if j == 0:
                    axes[i, j].set_ylabel(row)
                if i == 0:
                    regions = GenCoorSet(name=beds[0])
                    regions.load(filename=exp.get_file(beds[0]))
                    axes[i, j].set_title(col+" ("+str(len(regions))+")")


        cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
        fig.colorbar(hm, cax=cbar_ax)

        set_yaxis(n_row, n_col, axes, arg)
        fig.savefig(arg["<output_file>"], bbox_inches='tight')

    # elif arg["boxplot"]:
