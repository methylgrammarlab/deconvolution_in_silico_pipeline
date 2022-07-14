'''
This script takes in an epiread file and counts
how many CpGs each reads covers, then returns the mean.
'''
import json
from epiread_tools.naming_conventions import *
METHYLATION_COLUMN = 7

def count_coverage(epireads):
    n_lines = 0
    total_data = 0
    for epiread in epireads:
        with open(epiread, "r") as infile:
            for line in infile:
                n_lines += 1
                meth = line.split("\t")[METHYLATION_COLUMN]
                for m in meth:
                    if (methylation_state[m] == METHYLATED) or (methylation_state[m] == UNMETHYLATED):
                        total_data += 1
    return total_data/n_lines

def main(epireads, outpath):
    res = count_coverage(epireads)
    with open(outpath, "w") as outfile:
        json.dump(res, outfile)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('epireads', action="store", nargs="+", type=str)
parser.add_argument('--outfile', type=str)
args = parser.parse_args()
main(args.epireads, args.outfile)
