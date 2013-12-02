#!/usr/bin/env python
'''
###############################################################################

Bam Coverage

A program for generating coverage information for BAM sequence
alignment files.

Author: Bernie Pope (bjpope@unimelb.edu.au)
Licence: BSD

usage: bamcover [-h] [--version] --coords COORDS [--log FILE] bams [bams ...]

Generate coverage information for BAM files

positional arguments:
  bams             bam files containing mapped reads

optional arguments:
  -h, --help       show this help message and exit
  --version        show program's version number and exit
  --coords COORDS  TSV coordinates file (1-based) for region of interest
  --log FILE       Log progress in FILENAME.

###############################################################################
'''

from argparse import ArgumentParser
import logging
import sys
import pysam
import csv
from bx.intervals.intersection import IntervalTree
import pkg_resources
import matplotlib.pyplot as plt
import os


# The version number of the program as specified in the package.
VERSION = pkg_resources.require("bamcover")[0].version

# The maximum number of characters to take from the prefix of an input BAM
# file to use in the legend of the graph.
MAX_BAM_PREFIX = 10

# default name of log file if none is specified on the command line
DEFAULT_LOG_FILE = "bamcover.log"


def parse_args():
    'Parse the command line arguments for the program.'
    parser = ArgumentParser(
        description="Generate coverage information for BAM files")
    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + VERSION)
    parser.add_argument(
        '--coords', required=True, type=str,
        help='TSV coordinates file (1-based) for region of interest')
    parser.add_argument(
        'bams', nargs='+', type=str,
        help='bam files containing mapped reads')
    parser.add_argument(
        '--log', metavar='FILE', type=str, default=DEFAULT_LOG_FILE,
        help='Log progress in FILENAME.')
    return parser.parse_args()


def get_coords(coords_filename):
    '''Read the input coordinates file in TSV format and check
    that the start and end coordinates are integers. Log any
    coordinates which are skipped.'''
    result = []
    with open(coords_filename) as coords_file:
        reader = csv.reader(coords_file, delimiter='\t')
        for row in reader:
            if len(row) >= 3:
                chrom, start, end = row
                if start.isdigit() and end.isdigit():
                    result.append((chrom, int(start), int(end)))
            else:
                logging.warn("Skipping invalid coordinate: {}".format(row))
    return result


def bam_name_legend(bam_filename):
    '''Take the first MAX_BAM_PREFIX characters from the prefix of the input BAM
    file name to be used in the legend of the graph.'''
    basename = os.path.basename(bam_filename)
    return basename[:MAX_BAM_PREFIX]


def plot_coverage(coords, bams):
    '''Given the name of a DNA coordinates firl and a list of bam file names,
    plot the read aligment coverage for each bam file for each coordinate.
    One graph per coordinate will be generated. The coverage for each
    BAM file for a given coordinate will be plotted on the same graph.
    The coordinates file should be in TSV format.'''
    coords = get_coords(coords)
    for chrom, start, end in coords:
        logging.info("processing coord {} {} {}".format(chrom, start, end))
        # Start plotting the graph and generate a name for the output file
        graph_filename = start_graph(chrom, start, end)
        coords_range = range(start, end+1)
        for bam_filename in bams:
            # interval tree tracks the start and end mapped coordinates
            # of each read in the bam file that lies within our region
            # of interest.
            interval_tree = IntervalTree()
            with pysam.Samfile(bam_filename, "rb") as bam:
                logging.info("processing bam file {}".format(bam_filename))
                # Collect all the reads from the BAM file which lie in
                # the region of interest.
                # fetch uses 0-based indexing. Our input coordinates are
                # in 1-based coordinates.
                reads = bam.fetch(chrom, start-1, end-1)
                # Insert the start and end of each aligned read into the
                # interval tree.
                for read in reads:
                    if len(read.positions) > 0:
                        # Add 1 to convert from 0-based to 1-based coordinates
                        first_pos = read.positions[0] + 1
                        last_pos = read.positions[-1] + 1
                        interval_tree.add(first_pos, last_pos, None)
            # For each base position in our region of interest,
            # count the number of reads which overlap this position.
            # This computes the coverage for each position in the region.
            counts = [len(interval_tree.find(pos, pos))
                      for pos in coords_range]
            # Plot the coverage information for this bam file
            legend_text = bam_name_legend(bam_filename)
            plot_graph(counts, coords_range, legend_text)
        # Close the drawing of the graph for this set of coordinates
        end_graph(graph_filename)


def start_graph(chrom, start, end):
    '''Begin plotting of a graph for a given coordinate.'''
    plt.ylabel("Coverage")
    plt.xlabel("Position")
    plt.title("Read alignment coverage\n{}:{}-{}".format(chrom, start, end))
    return "{}.{}.{}.png".format(chrom, start, end)


def plot_graph(counts, coords_range, bam_name):
    '''Plot the coverage data for a given bam file.'''
    plt.plot(coords_range, counts, label=bam_name)


def end_graph(filename):
    '''Begin plotting of a graph for a given coordinate.'''
    plt.legend(title="Sample")
    plt.savefig(filename)
    plt.close()


def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOGFILE.'''
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    # Log the command line that was used to run the program
    logging.info('command line: {0}'.format(' '.join(sys.argv)))


def main():
    '''Program entry point.'''
    args = parse_args()
    start_log(args.log)
    plot_coverage(args.coords, args.bams)


if __name__ == '__main__':
    main()
