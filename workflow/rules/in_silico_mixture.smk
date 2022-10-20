
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
import numpy as np
import pandas as pd

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

def count_lines(wildcards):
    with open(wildcards) as foo:
        lines = len(foo.readlines())
    return lines

rule measure_length: #TODO: replace with something better?
    input:
        epireads=expand("interim/{cell_type}_mixture_epipaths.epiread", cell_type=config["cell_types"])
    output:
        expand("interim/{name}_mean_cpgs_per_read.json", name=config["name"])
    shell:
        """python3 workflow/scripts/measure_read_length.py {input.epireads} --outfile {output}"""


rule calculate_reads: #number to sample from each cell type
    input:
        cpgs=expand("results/{name}_atlas_over_regions.txt", name=config["name"]), #number of cpgs to cover
        lengths=expand("interim/{name}_mean_cpgs_per_read.json", name=config["name"]) #TODO: where is this?
    params:
        coverage=lambda wildcards: str(runs.depth[runs.index == int(wildcards.param_id)].values[0]),
        alpha = lambda wildcards: str(runs.true_alpha[runs.index == int(wildcards.param_id)].values[0]),
        cell_types=str(config["cell_types"])
    output:
        expand("interim/{name}_{{param_id}}_rep{{instance_id}}_read_number.json", name=config["name"])
    run:
        n_cpgs = count_lines(next(shell("echo {input.cpgs}", iterable=True)))
        shell("""python3 workflow/scripts/mixture_simulator.py {input.lengths} {params.coverage}"""+\
              """ {n_cpgs} '{params.alpha}' '{params.cell_types}' {output}""")
    #this works: python3 mixture_simulator.py 10 584 '[0.04761905, 0.0952381 ,
    # 0.14285714, 0.19047619, 0.23809524, 0.28571429]' '[Alpha, Beta, Delta, Duct, Acinar, Endothel]' try.json

rule sample_from_epiread:
    input:
        epiread="interim/{cell_type}_mixture_epipaths.epiread",
        counts=expand("interim/{name}_{{param_id}}_rep{{instance_id}}_read_number.json", name=config["name"])
    params:
        cell_type="{cell_type}"
    output:
        "interim/{cell_type}_{param_id}_rep{instance_id}_random_reads.epiread"
    run:
        with open(input.counts[0], "r") as infile:
            counts = int(json.load(infile)[params.cell_type])
            if counts > 0:
                shell("""bedtools sample -i {input.epiread} -n {counts} | sort -k1,1 -k2,2n > {output}""")
            else:
                shell("touch {output}")

rule make_mixture:
    input:
        expand("interim/{cell_type}_{{param_id}}_rep{{instance_id}}_random_reads.epiread", cell_type=config["cell_types"])
    output:
        expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread", name=config["name"])
    run:
        shell("""sort -k1,1 -k2,2n -m {input} > {output}""")

rule zip_output:
    input:
        "{path}.epiread"
    output:
        "{path}.epiread.gz"
    shell:
        """bgzip -c {input} > {output}"""


rule index_output:
    input:
        "{path}.epiread.gz"
    output:
        "{path}.epiread.gz.tbi"
    shell:
        """tabix -p bed {input}"""

