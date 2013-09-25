#!/usr/bin/env python
'''
Bam Coverage

A program for generating coverage information for BAM sequence
alignment files.

Author: Bernie Pope (bjpope@unimelb.edu.au)
Licence: BSD

Usage: 
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

version = pkg_resources.require("bamcover")[0].version

def parse_args():
    parser = ArgumentParser(description="Generate coverage information for BAM files")
    parser.add_argument(
    '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument(
    '--coords', required=True, type=str, help='DNA coordinates to generate coverage for')
    parser.add_argument(
        'bams', nargs='+', type=str, help='bam files containing mapped reads')
    parser.add_argument( '--log', metavar='FILE', type=str,
        help='Log progress in FILENAME, defaults to stdout.')
    return parser.parse_args() 

def get_coords(coords_filename):
    result = []
    with open(coords_filename) as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) >= 3:
                chr, start, end = row
                if start.isdigit() and end.isdigit():
                    result.append((chr, int(start), int(end)))
    return result

def bam_name_legend(bam_filename):
    basename = os.path.basename(bam_filename)
    return basename[:10]

def plot_coverage(coords, bams):
    coords = get_coords(coords)
    for chr, start, end in coords:
        graph_filename = start_graph(chr, start, end)
        coords_range = range(start, end+1)
        for bam_filename in bams:
            interval_tree = IntervalTree()
            with pysam.Samfile(bam_filename, "rb") as bam:
                logging.info("processing bam file {}".format(bam_filename))
                reads = bam.fetch(chr, start, end + 1)
                for read in reads:
                    if len(read.positions) > 0:
                        first_pos = read.positions[0]
                        last_pos = read.positions[-1]
                        interval_tree.add(first_pos, last_pos, None)
            counts = [len(interval_tree.find(pos, pos)) for pos in coords_range]
            plot_graph(counts, coords_range, bam_name_legend(bam_filename))
        end_graph(graph_filename)
            

def start_graph(chr, start, end):
    plt.ylabel("Coverage")
    plt.xlabel("Position")
    plt.title("Read alignment coverage\n{}:{}-{}".format(chr, start, end))
    return "{}.{}.{}.png".format(chr, start, end)
            

def plot_graph(counts, coords_range, bam_num):
    plt.plot(coords_range, counts, label=bam_num)


def end_graph(filename):
    plt.legend(title="Sample")
    plt.savefig(filename)
    plt.close()


def main():
    args = parse_args()
    if args.log is None:
        logfile = sys.stdout
    else:
        logfile = args.log
    logging.basicConfig(
        filename=args.log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))
    plot_coverage(args.coords, args.bams)


if __name__ == '__main__':
    main()
