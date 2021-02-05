#!/usr/bin/env python3
"""
usage: pacman.py [-h|--help] [-o OUTPUT] [OPTIONS] [FILE OR DIR]

Plot I24 chip results.

positional arguments:
  FILE OR DIR           Plot results from a hitfinding file (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        The output filename to write plot to (default: None)
  --show                Always show plot, even if specifying output filename.
                        (default: False)
  --binding {alpha,shot}
                        Binding type, if appropriate (default: None)
  --column COLUMN       Which column from the input file to read. Accepts
                        numbers (for columnar .out files) or names (for new
                        json-style files) [default: 3 or total_intensity])
                        (default: 3)
  --chiptype CHIPTYPE   The chip type (default: 1)
  --blocks BLOCKS       Limit to specified blocks e.g. A1,A2,A4 (default:
                        None)
  --ms MS               Plotting marker size (default: 8)
  --cmap CMAP           Matplotlib cmap (default: terrain)
  --dpi DPI             Matplotlib save resolution (default: 200)
  --xlim XLIM           Matplotlib x range to plot (default: (-1, 25))
  --ylim YLIM           Matplotlib y range to plot (default: (-1, 25))
  --zlim ZLIM           Matplotlib z range to plot (default: None)
  --alpha ALPHA         Matplotlib marker plotting alpha (default: 1)
  --method {new,old}    Plotting method. "old" plots markers for each well,
                        which causes overlap and bad results when zooming.
                        "new" causes images with a pixel per well to be drawn
                        for each block. (default: new)

For legacy compatibility, any option can be passed in as optionname=value
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Chip_StartUp_v5 import (
    address_to_index,
    get_alphanumeric,
    get_format,
    get_shot_order,
    get_xy,
)


def hits_scrape(filename, column_choice, shot_order_addr_list, bound=False):
    # type (str, Union[str,int], Any, bool) -> Dict
    """Read in hits from a hit file."""

    # Minor fudge: this is set like this so that we get a nice "default"
    # in the automatic help string....
    if column_choice == "3 or total_intensity":
        column_choice = None

    with open(filename) as f:
        raw_line_data = f.readlines()

    # Decide what kind of file this is... old columnar, or json stream
    try:
        column_choice = column_choice or "total_intensity"
        data = [json.loads(line) for line in raw_line_data]
        # Correct any case sensitivity issues
        column_choice = [
            x for x in data[0].keys() if x.lower() == column_choice.lower()
        ][0]
        print(
            "COLUMN CHOICES:\n"
            + "\n".join(
                (" -> " if x == column_choice else "    ") + x
                for x in data[0].keys()
                if isinstance(data[0][x], (int, float))
            )
        )
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

    except json.JSONDecodeError:
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


class ChipHitPlotter:
    """
    Handles tracking and plotting of blocks.

    Blocks are drawn as a single bitmap image - so each well is distinct
    and there are absolutely no overlaps.
    """

    def __init__(self, chip_type: str) -> None:
        self.chip_type = chip_type
        self.metrics = get_format(chip_type)

        # Create the array of block contents
        self.blocks = np.full(
            (
                self.metrics.blocks_y,
                self.metrics.blocks_x,
                self.metrics.wells_y,
                self.metrics.wells_x,
            ),
            -1.0,
        )

        # Vestigial code to colour the top of each block
        # for col in range(self.metrics.blocks_x):
        #     for row in range(self.metrics.blocks_y):
        #         self.blocks[row, col][:5, :5] = 0.4e6

    def set_from_dict(self, hits_addr_dict: Dict[str, float]) -> None:
        """Fill out all the wells from an address dictionary"""
        for addr, val in hits_addr_dict.items():
            self.set_well(addr, val)

    def set_well(self, addr: str, value: float) -> None:
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
                        **kwargs,
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


def run_fromdir_method(path):
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
            sys.exit(f"Error: {path} is not a valid path")
    else:
        raise SyntaxError("Unknown directory\n\n\n")
        return 0


def run_fromfile_method(
    filename, column=None, binding=None, chiptype="1", blocks=None, **kwargs
):
    blocks = blocks or []
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

    print(f"Length of address list:{len(addr_list)}\n")
    new_addr_list = []
    if len(blocks) > 0:
        for addr in addr_list:
            block = addr.split("_")[0]
            if block in blocks:
                new_addr_list.append(addr)
        addr_list = new_addr_list
        print("After:            ", len(addr_list))

    if filename.suffix == ".out":
        hits_addr_dict = hits_scrape(filename, column, addr_list, bound)
        plotter = ChipHitPlotter(chiptype)
        plotter.set_from_dict(hits_addr_dict)
        x, y, z = make_plot_arrays(hits_addr_dict, chiptype)
        return x, y, z, plotter
    elif filename.suffix == ".txt":
        hits_addr_dict = det_dist_file_scrape(filename, addr_list)
        plotter = ChipHitPlotter(chiptype)
        plotter.set_from_dict(hits_addr_dict)
        x, y, z = make_plot_arrays(hits_addr_dict, chiptype)
        return x, y, z, plotter

    sys.exit("Unrecognised file extension: Expect .txt or .out")


def plot(x, y, z, plotter, options):
    default_output = not options.output
    if not options.output:
        options.output = f"{options.file.stem}_{options.column}_{time.strftime('%Y%m%d_%H%M%S')}map.png"

    fig = plt.figure(figsize=(10, 10), facecolor="0.75", edgecolor="w")
    ax1 = fig.add_subplot(111, aspect=1, facecolor="0.5")
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.03, top=0.97)
    ax1.set_xlim(options.xlim)
    ax1.set_ylim(options.ylim)

    if options.zlim:
        z = np.clip(z, *options.zlim)

    if options.method == "old":
        plt.scatter(
            x,
            y,
            c=[float(zcol) for zcol in z],
            s=options.ms,
            marker="s",
            vmin=options.zlim[0],
            vmax=options.zlim[1],
            alpha=options.alpha,
            cmap=options.cmap,
        )
    else:
        assert plotter
        plotter.draw(plt, cmap=options.cmap)

    ax1.invert_yaxis()

    # Make the colorbar match the height of the image
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.07)
    plt.colorbar(cax=cax)

    plt.tight_layout()
    plt.savefig(
        options.output,
        dpi=options.dpi,
        # facecolor="w",
        bbox_inches="tight",
        pad_inches=0,
        transparent=False,
    )

    print(f"Imaged saved as: {options.output}")

    # If we asked to or are writing a generated filename, show
    if default_output or options.show:
        plt.show()


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Plot I24 chip results.",
        epilog="For legacy compatibility, any option can be passed in as optionname=value",
        usage=__doc__.strip().splitlines()[0][6:].strip(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "file",
        metavar="FILE OR DIR",
        type=Path,
        help="Plot results from a hitfinding file",
        nargs="?",
    )
    parser.add_argument("-o", "--output", help="The output filename to write plot to")
    # dir= and file= can be passed via option form but are deprecated
    parser.add_argument("--dir", "--file", help=argparse.SUPPRESS, type=Path)
    parser.add_argument(
        "--show",
        help="Always show plot, even if specifying output filename.",
        action="store_true",
    )
    parser.add_argument(
        "--binding", choices=["alpha", "shot"], help="Binding type, if appropriate"
    )
    parser.add_argument(
        "--column",
        default="3 or total_intensity",
        help="Which column from the input file to read. Accepts numbers (for columnar .out files) or names (for new json-style files) [default: 3 or total_intensity])",
    )
    parser.add_argument("--chiptype", default="1", help="The chip type", type=str)
    parser.add_argument("--blocks", help="Limit to specified blocks e.g. A1,A2,A4")
    parser.add_argument("--ms", default=8, type=int, help="Plotting marker size")
    parser.add_argument("--cmap", default="terrain", help="Matplotlib cmap")
    parser.add_argument(
        "--dpi", default=200, type=int, help="Matplotlib save resolution"
    )

    def _range(opt):
        """Parse an argument string as two separate floats"""
        a, b = opt.split(",")
        return float(a), float(b)

    parser.add_argument(
        "--xlim", type=_range, default=(-1, 25), help="Matplotlib x range to plot"
    )
    parser.add_argument(
        "--ylim", type=_range, default=(-1, 25), help="Matplotlib y range to plot"
    )
    parser.add_argument("--zlim", type=_range, help="Matplotlib z range to plot")
    parser.add_argument(
        "--alpha", default=1, type=float, help="Matplotlib marker plotting alpha"
    )
    parser.add_argument(
        "--method",
        default="new",
        choices=["new", "old"],
        help="""Plotting method. "old" plots markers for each well, which causes overlap and bad results when zooming. "new" causes images with a pixel per well to be drawn for each block.""",
    )

    # Pre-convert anything option=like to an --option=like so that we
    # remain (mostly) compatible with the old way of running pacman.
    args = args or sys.argv[1:]
    args = [f"--{x}" if ("=" in x and not x.startswith("-")) else x for x in args]
    options = parser.parse_args(args)

    # Handle redundancy of file/file=/dir=
    if options.file and options.dir:
        sys.exit("Error: Must specify file or dir, not both")
    elif options.dir:
        print(
            "Warning: file= and dir= options are deprecated. Pass path directly instead."
        )
        options.file = options.dir
    elif not options.file:
        parser.print_usage()
        sys.exit("Error: Must specify path to process")

    if not options.file.exists():
        sys.exit(f"Error: {options.file} is not a file or a directory")

    plotter = None

    if options.file.is_file():
        x, y, z, plotter = run_fromfile_method(options.file, **vars(options))
    else:
        x, y, z = run_fromdir_method(options.file)

    plot(x, y, z, plotter, options)


if __name__ == "__main__":
    main()
