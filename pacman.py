#!/dls_sw/apps/python/anaconda/1.7.0/64/bin/python
"""
usage: pacman [-h|--help] [FILE or DIR] [OPTIONS]

Plot I24 chip results.

Methods:
    FILE                        Plot results from a hitfinding file
    DIR
    file=/path/to/file.out      Alternate method to plot from file
    dir=/path/to/dir

Additional options can be passed:
    binding=[alpha|shot]    Binding type, if appropriate
    col=3                   Which column from the input file to read. Accepts
                            numbers (for columnar .out files) or names (for new
                            json-style files) [default: 3 or total_intensity]
    chiptype=1              The chip type. [default: 1]
    blocks=                 Limit to specified blocks e.g. A1,A2,A4
    ms=8                    Plotting marker size [default: 8]
    cmap=terrain            Matplotlib cmap. [default: terrain]
    dpi=200                 Matplotlib save resolution, in DPI [default: 200]
    xlim=min,max            Matplotlib x range to plot
    ylim=min,max            Matplotlib y range to plot
    zlim=min,max            Matplotlib z range to plot
    alpha=1                 Matplotlib marker plotting alpha [default: 1]
    method=[new|old]        Plotting method. "old" plots markers for each well,
                            which causes overlap and bad results when zooming.
                            "new" causes images with a pixel per well to be
                            drawn for each block. [default: new]
"""

from __future__ import absolute_import, print_function

import json
import os
import sys
import time
from collections import namedtuple

import numpy as np
from matplotlib import pyplot as plt

from Chip_StartUp_v5 import (
    address_to_index,
    get_alphanumeric,
    get_format,
    get_shot_order,
    get_xy,
)

try:
    from typing import Dict
except ImportError:
    pass


def hits_scrape(filename, column_choice, shot_order_addr_list, bound=False):
    # type (str, Union[str,int], Any, bool) -> Dict
    with open(filename) as f:
        raw_line_data = f.readlines()

    # Decide what kind of file this is... old columnar, or json stream
    try:
        data = [json.loads(line) for line in raw_line_data]
        print(
            "COLUMN CHOICES:\n"
            + "\n".join(
                "    " + x
                for x in data[0].keys()
                if isinstance(data[0][x], (int, float))
            )
        )
        column_choice = column_choice or "total_intensity"
        # Do a sanity check that the hit file doesn't contain more hits than we
        # have addresses. Less is probably better - a partial chip might be in
        # the same address order?
        if "expected_total" in data[0] and int(data[0]["expected_total"]) > len(
            shot_order_addr_list
        ):
            sys.exit(
                "Error: {} expected hits in {} but only {} expected addresses for chip".format(
                    data[0]["expected_total"], filename, len(shot_order_addr_list)
                )
            )
        # This data file isn't in hit-order, because the hitfinding service
        # returns asynchronously. So use the image index as a lookup into
        # the address list. This will also solve problems of missing hits/
        # skipped images
        hits_dict = {}
        for hit in data:
            hits_dict[shot_order_addr_list[hit["file-pattern-index"]]] = hit[
                column_choice
            ]

    except ValueError:  # json.JSONDecodeError:
        # We don't have a jsony file, we have a columnar one
        columns = raw_line_data[1].split("|")
        column_choice = 3 if column_choice is None else column_choice

        print("COLUMN CHOICES")
        for i, item in enumerate(columns[2:10]):
            print("Column", i + 1, item)

        hits_dict = {}
        for i, line in enumerate(raw_line_data[3:]):
            # End of table indicated by ------
            if line.startswith("-"):
                break
            line.rstrip("\n")
            listy = line.split("|")
            entry = [x.lstrip().rstrip() for x in listy[1:]]
            if not entry[column_choice]:
                continue
            addr = shot_order_addr_list[i]
            if bound:
                hits_dict[addr] = i
            else:
                hits_dict[addr] = entry[column_choice]
    return hits_dict


def det_dist_file_scrape(fid, shot_order_addr_list):
    print(shot_order_addr_list[:10])

    det_dict = {}
    f = open(fid)
    multi = 0
    for line in f.readlines():
        entry = line.split()
        img_num = int(entry[0])
        # try:
        #   x = det_dict[img_num]
        #   multi +=1
        # except:
        if entry[1] == "None":
            val = None
        else:
            val = float(entry[1])
        if val == 0:
            val = None
        det_dict[img_num] = val
    print(multi)
    hits_dict = {}
    for i in range(len(shot_order_addr_list)):
        addr = shot_order_addr_list[i]
        if i in list(det_dict.keys()):
            val = det_dict[i]
        else:
            val = None
        hits_dict[addr] = val
    return hits_dict


ChipFormat = namedtuple(
    "ChipFormat",
    [
        "blocks_x",
        "blocks_y",
        "wells_x",
        "wells_y",
        "well_distance",
        "block_distance_x",
        "block_distance_y",
    ],
)


class ChipHitPlotter(object):
    """
    Handles tracking and plotting of blocks.

    Blocks are drawn as a single bitmap image - so each well is distinct
    and there are absolutely no overlaps.
    """

    def __init__(self, chip_type):
        # type: (str) -> None
        self.chip_type = chip_type
        self.metrics = ChipFormat(*get_format(chip_type))

        # Create the array of block contents
        self.blocks = np.full(
            (
                self.metrics.blocks_y,
                self.metrics.blocks_x,
                self.metrics.wells_y,
                self.metrics.wells_x,
            ),
            -1,
        )
        for col in range(self.metrics.blocks_x):
            for row in range(self.metrics.blocks_y):
                self.blocks[row, col][:5, :5] = 0.4e6

    def set_from_dict(self, hits_addr_dict):
        # type: (Dict[str, float]) -> None
        """Fill out all the wells from an address dictionary"""
        for addr, val in hits_addr_dict.items():
            self.set_well(addr, val)

    def set_well(self, addr, value):
        # type: (str, float) -> None
        """Set the value of a well.

        Args:
            addr: The address string of the well e.g. A1_aj
            value: The value to set the cell to
        """
        blockR, blockC, windowR, windowC = address_to_index(addr)
        self.blocks[blockR, blockC, windowR, windowC] = value

    def draw(self, fig, **kwargs):
        """
        Draw a chip.

        Args:
            fig: The target to draw on e.g. a figure, or pyplot
            kwargs: Any matplotlib imshow drawing argument
        """
        vmin = kwargs.pop("vmin", np.max(np.min(self.blocks), 0))
        vmax = kwargs.pop("vmax", np.max(self.blocks))

        block_width = (self.metrics.wells_x - 1) * self.metrics.well_distance
        block_height = (self.metrics.wells_y - 1) * self.metrics.well_distance
        # Edge-to-edge well coverage
        well_extents_x = self.metrics.wells_x * self.metrics.well_distance
        well_extents_y = self.metrics.wells_y * self.metrics.well_distance
        image_artists = []
        for col in range(self.metrics.blocks_x):
            for row in range(self.metrics.blocks_y):
                # Work out position of top-left corner of well
                x = (
                    col * (block_width + self.metrics.block_distance_x)
                    - self.metrics.well_distance / 2
                )
                y = (
                    row * (block_height + self.metrics.block_distance_y)
                    - self.metrics.well_distance / 2
                )
                image = self.blocks[row, col]
                image_masked = np.ma.masked_where(image < 0, image)
                image_artists.append(
                    fig.imshow(
                        image_masked,
                        extent=[x, x + well_extents_x, y, y + well_extents_y],
                        vmin=vmin,
                        vmax=vmax,
                        origin="lower",
                        **kwargs
                    )
                )
        return image_artists


def make_plot_arrays(hits_addr_dict, chip_type):
    x_list, y_list, z_list = [], [], []
    for addr, val in hits_addr_dict.items():
        x, y = get_xy(addr, chip_type)
        x_list.append(x)
        y_list.append(y)
        z_list.append(val)

    X = np.array(x_list)
    Y = np.array(y_list)
    Z = np.array(z_list)

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()
    return x, y, z


def run_fromdir_method(*args):
    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
        if k.startswith("dir"):
            path = v
            print("\n\n\nThis is the path:", path)

    addr_list = get_shot_order(chip_type="1")
    hits_dict = {}
    if path:
        if os.path.isdir(path):
            fid_list = os.listdir(path)
            int_list = [fid for fid in fid_list if fid.startswith("int-0")]
            print(int_list[:10])
            shotnum_list = [
                int(fid.split("_")[1].rstrip(".pickle")) for fid in int_list
            ]
            print("Number of Integrated:", len(shotnum_list))
            for x in shotnum_list[:20]:
                print(x, end=" ")

            for i, addr in enumerate(addr_list):
                if i in shotnum_list:
                    hits_dict[addr] = 1
                else:
                    hits_dict[addr] = 0
            x, y, z = make_plot_arrays(hits_dict, chip_type="1")
            return x, y, z
        else:
            sys.exit("Error: {} is not a valid path".format(path))
    else:
        raise SyntaxError("Unknown directory\n\n\n")
        return 0


def run_fromfile_method(
    filename, col=None, binding=None, chiptype="1", blocks=[], **kwargs
):
    if isinstance(blocks, str):
        blocks = blocks.split()
        for block in blocks:
            print(block)

    if binding:
        if binding.startswith("alpha"):
            addr_list = get_alphanumeric(chiptype)
        elif binding.startswith("shot"):
            addr_list = get_shot_order(chiptype)
        bound = True
    else:
        addr_list = get_shot_order(chiptype)
        bound = False

    print("Length of address list:%s\n" % len(addr_list))
    new_addr_list = []
    if len(blocks) > 0:
        for addr in addr_list:
            block = addr.split("_")[0]
            if block in blocks:
                new_addr_list.append(addr)
        addr_list = new_addr_list
        print("After:            ", len(addr_list))

    if filename.endswith("out"):
        hits_addr_dict = hits_scrape(filename, col, addr_list, bound)
        plotter = ChipHitPlotter(chiptype)
        plotter.set_from_dict(hits_addr_dict)
        x, y, z = make_plot_arrays(hits_addr_dict, chiptype)
        return x, y, z, plotter
    elif filename.endswith(".txt"):
        hits_addr_dict = det_dist_file_scrape(filename, addr_list)
        plotter = ChipHitPlotter(chiptype)
        plotter.set_from_dict(hits_addr_dict)
        x, y, z = make_plot_arrays(hits_addr_dict, chiptype)
        return x, y, z, plotter

    sys.exit("Unrecognised file extension: Expect .txt or .out")


def plot(x, y, z, plotter, *args):
    alfa = 1
    col = 3
    mrksz = 8
    pixels = 200
    cmap_choice = "terrain"
    xlim_min, xlim_max = (-1, 25)
    ylim_min, ylim_max = (-1, 25)
    zlim_min, zlim_max = (None, None)
    method = "new"

    for arg in args:
        k = arg.split("=")[0]
        v = arg.split("=")[1]
        if k.startswith("file"):
            fid = v
        if k.startswith("dir"):
            fid = v
        elif k.startswith("col"):
            col = str(v)
        elif k.startswith("ms"):
            mrksz = int(v)
        elif k.startswith("alpha"):
            alfa = float(v)
        elif k.startswith("cmap"):
            cmap_choice = str(v)
        elif k.startswith("dpi"):
            pixels = int(v)
        elif k.startswith("xlim"):
            xlim_min, xlim_max = v.split(",")
        elif k.startswith("ylim"):
            ylim_min, ylim_max = v.split(",")
        elif k.startswith("zlim"):
            zlim_min, zlim_max = v.split(",")
            zlim_min = float(zlim_min)
            zlim_max = float(zlim_max)
        elif k.startswith("method"):
            method = v

    raw_fid = fid.split("/")[-1]
    print("\nraw_fid", raw_fid)
    if raw_fid == "":
        raw_fid = fid.split("/")[-2]

    out_fid = "%s_%s_%smap.png" % (
        raw_fid.split(".")[0],
        col,
        time.strftime("%Y%m%d_%H%M%S"),
    )
    fig = plt.figure(figsize=(10, 10), facecolor="0.75", edgecolor="w")
    ax1 = fig.add_subplot(111, aspect=1, facecolor="0.5")
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97)

    # plt.axis('off')
    ax1.set_xlim(float(xlim_min), float(xlim_max))
    ax1.set_ylim(float(ylim_min), float(ylim_max))

    if not (zlim_min, zlim_max) == (None, None):
        a_list = []
        for entry in z.tolist():
            if float(entry) > zlim_max:
                a_list.append(zlim_max)
            elif float(entry) < zlim_min:
                a_list.append(zlim_min)
            else:
                a_list.append(entry)
        z = np.array(a_list)

    # plt.scatter(x, y, c=z, s=mrksz, marker='s',vmin=95, vmax=98.5, alpha=alfa, cmap=cmap_choice)

    if method == "old":
        plt.scatter(
            x,
            y,
            c=[float(zcol) for zcol in z],
            s=mrksz,
            marker="s",
            vmin=zlim_min,
            vmax=zlim_max,
            alpha=alfa,
            cmap=cmap_choice,
        )
    else:
        plotter.draw(plt, cmap=cmap_choice)

    ax1.invert_yaxis()
    plt.colorbar()
    path = "."
    plt.savefig(
        path + out_fid,
        dpi=pixels,
        facecolor="0.15",
        bbox_inches="tight",
        pad_inches=0,
        transparent=True,
    )
    print("\n\nPLOTTING KEYWORDS\n", 30 * "-")
    print("col=%s" % col)
    print("ms=%s" % mrksz)
    print("cmap=%s" % cmap_choice)
    print("dpi=%s" % pixels)
    print("xlim=%s,%s" % (xlim_min, xlim_max))
    print("ylim=%s,%s" % (ylim_min, ylim_max))
    print("zlim=%s,%s" % (zlim_min, zlim_max))
    print("alpha=%s" % alfa)
    print(30 * "-", "\n\n")
    print("Imaged saved as: %s%s" % (path, out_fid))
    plt.show()


def run_dosecolumns(*args):
    path = "/dls/i24/data/2017/nt14493-54/processing/find_spots/hector_pacman/"
    chip_name = "hector"
    num_of_doses = 5
    list_of_lists = []
    for x in range(num_of_doses):
        list_of_lists.append([])

    if os.path.isdir(path):
        fid_list = sorted(os.listdir(path))
        list_of_files = [fid for fid in fid_list if fid.startswith(chip_name)]
        for fid in list_of_files:
            print(fid)
            f = open(path + fid, "r")
            for i, line in enumerate(f.readlines()[3:]):
                # print i, i % num_of_doses
                if line.startswith("-"):
                    break
                else:
                    list_of_lists[i % num_of_doses].append(line)
            f.close()

    f = open(path + fid, "r")
    header = f.readlines()[:3]
    f.close()

    for i in range(num_of_doses):
        new_fid = chip_name + ".dose%s" % (i + 1)
        print("Making:", new_fid)
        g = open(new_fid, "a")
        for line in header:
            print(line.rstrip())
            g.write(line)
        for line in list_of_lists[i]:
            g.write(line)
        g.close()
    print()


def main(args=None):
    # Manual parsing because currently rather different from normal argparse
    args = args or sys.argv[1:]
    if not args:
        print(__doc__.strip().splitlines()[0])
        sys.exit(1)
    if "-h" in args or "--help" in args:
        print(__doc__.strip())
        sys.exit()

    allowed_keyword_list = [
        "file",
        "dir",
        "dosecolumns",
        "alpha",
        "binding",
        "blocks",
        "chiptype",
        "cmap",
        "col",
        "dpi",
        "ms",
        "xlim",
        "ylim",
        "zlim",
        "method",
    ]
    # Convert the list of a=b c=d arguments to a dictionary by splitting on "="
    arg_dict = {key: val for key, val in [x.split("=", 1) for x in args if "=" in x]}
    # Handle any positional, including identifying duplications
    positionals = [x for x in args if "=" not in x]
    if len(positionals) > 1 or (
        positionals and ("file" in arg_dict or "dir" in arg_dict)
    ):
        sys.exit("Error: Can only specify one file or directory location")
    if positionals:
        positional = positionals[0]
        args.remove(positional)
        if os.path.isfile(positional):
            arg_dict["file"] = positional
            args.append("file={}".format(positional))
        elif os.path.isdir(positional):
            arg_dict["dir"] = positional
            args.append("dir={}".format(positional))
        else:
            sys.exit("Error: Could not identify {} as a file or dir".format(positional))

    # Validate all passed arguments are in our allowlist
    for key in arg_dict.keys():
        print("Keyword argument: %s=%s" % (key, arg_dict[key]))
        if key not in allowed_keyword_list:
            sys.exit(
                "Unknown arg in args: {}\n    Allowed keywords: {}".format(
                    key, ", ".join(allowed_keyword_list)
                )
            )
    if "file" in arg_dict and "dir" in arg_dict:
        sys.exit("Error: Must specify file or dir, not both")

    if "file" in arg_dict:
        method = "fromfile"
    elif "dir" in arg_dict:
        method = "fromdir"
    elif "dosecolumns" in arg_dict:
        method = "dosecolumns"
    else:
        sys.exit("Error: Cannot determine method. Please pass file= or dir=")

    print("\n\nCalled PACMAN method = %s\n" % method)
    print("\n\nDATA KEYWORDS\n", 30 * "-")
    print("file=filename.out")
    print("blocks=A1,A2,A3,A4 ...................,,,,.. Blocks MUST be in block order")
    print("          Total number of images in file must not exceed 200*num_of_blocks")
    print(30 * "-", "\n")

    plotter = None

    if method == "fromfile":
        x, y, z, plotter = run_fromfile_method(arg_dict["file"], **arg_dict)
    elif method == "fromdir":
        x, y, z = run_fromdir_method(*args)
    elif method == "dosecolumns":
        run_dosecolumns(*args)
        print("Exiting")
        return 0
    else:
        raise SyntaxError("Unknown method")

    plot(x, y, z, plotter, *args)
    print("EOP")


if __name__ == "__main__":
    main()
