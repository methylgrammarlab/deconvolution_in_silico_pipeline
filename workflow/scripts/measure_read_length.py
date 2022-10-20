'''
This script takes in an epiread file and counts
how many CpGs each reads covers, then returns the mean.
'''

# ###################################################
#
# MIT License
#
# Copyright (c) 2022 irene unterman and ben berman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###################################################

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
