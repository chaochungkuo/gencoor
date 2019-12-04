#!/usr/bin/env python
"""
Usage:
    coordinates.py resize <length> <input_BED_file_path> <output_BED_file_path>
    coordinates.py split_strand <input_BED_file_path> <output_BED_directory_path>

    Blablabla

Options:

"""

from docopt import docopt
from gencoor.coordinates import GenCoorSet
import os

if __name__ == '__main__':
    arg = docopt(__doc__)
    if arg["resize"]:
        gc = GenCoorSet(name="input")
        gc.load(arg["<input_BED_file_path>"], filetype="BED")
        gc.relocate(mode='center as center', width=int(arg["<length>"]))
        gc.save(arg["<output_BED_file_path>"], filetype="BED")

    elif arg["split_strand"]:
        name = os.path.basename(arg["<input_BED_file_path>"]).split(".")[0]
        print(os.path.join(arg["<output_BED_directory_path>"], name + "_" + "+" + ".bed"))
        gc = GenCoorSet(name="input")
        gc.load(arg["<input_BED_file_path>"], filetype="BED")
        res = gc.split_by_strands()
        for k, g in res.items():
            g.save(os.path.join(arg["<output_BED_directory_path>"],
                                name + "_" + k + ".bed"),
                   filetype="BED")
