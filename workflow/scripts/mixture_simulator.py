'''
calculate number of reads to get from each cell type
total reads = number of regions (TIMs) * coverage
draw number per cell with p=alpha
save to file
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

import numpy as np
import argparse
import json
from collections import Counter
import re
#%%

def main(args):
    with open(args.lengths, "r") as infile:
        factor = json.load(infile)
    total_reads = int(args.coverage * args.n_cpgs / factor)
    assert total_reads > 0
    z = np.random.choice(args.cell_types, p=args.alpha, size=total_reads)
    res = dict([(x, 0) for x in args.cell_types])
    res.update(dict(Counter(z)))
    with open(args.outfile, "w") as outfile:
        json.dump(res, outfile)

parser = argparse.ArgumentParser()

parser.add_argument("lengths", type=str) #json with mean # cpgs covered by read
parser.add_argument("coverage", type=float)
parser.add_argument("n_cpgs", type=int)
parser.add_argument('alpha', type=str)
parser.add_argument('cell_types', type=str)
parser.add_argument('outfile', type=str)

args = parser.parse_args()
args.alpha = np.array([float(x) for x in re.split(',', args.alpha.strip("[]'").replace(" ",""))])
args.alpha = args.alpha / args.alpha.sum()
args.cell_types = re.split(',', args.cell_types.strip("[]").replace(" ",""))

main(args)

