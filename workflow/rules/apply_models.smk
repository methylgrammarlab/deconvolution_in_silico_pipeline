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
import pandas as pd
import json


runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index


def get_regions_file(wildcards):
    if len(config["regions_file"]):
        return config["regions_file"]
    else:
        return expand("results/{name}_tims.txt", name=config["name"])[0]

rule sort_regions_file:
    input:
         regions=get_regions_file
    output:
        temp("sorted_regions_file.bed")
    shell:
        """sort -k1,1 -k2,2n {input.regions} > {output}"""

rule merge_regions_file:
    input:
        "sorted_regions_file.bed"
    output:
        "merged_sorted_regions_file.bed"
    shell:
        """bedtools merge -i {input} > {output}"""

rule create_run_config:
    input:
        mixture=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz", name=config["name"]),
        index=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz.tbi", name=config["name"]),
        atlas=expand("results/{name}_atlas_over_regions.txt", name=config["name"]),
        lambdas= expand("results/{name}_lambdas.bedgraph",name=config["name"]),,
        thetas= expand("results/{name}_thetas.bedgraph",name=config["name"]),
        regions=get_regions_file #not merged
    params:
        run_config = lambda wildcards: runs[runs.index == int(wildcards.param_id)].iloc[0, :].to_dict()
    output:
        temp("interim/{param_id}_rep{instance_id}_run_config.json")
    run:
        basic_config = {"bedfile":True, "header":False,"cpg_coordinates":config["cpg_file"]}
        basic_config.update(params.run_config)
        basic_config["epiread_files"] = input.mixture
        basic_config["epiformat"] = config["epiformat"]
        basic_config["atlas_file"] = input.atlas[0]
        basic_config["genomic_intervals"] = input.regions
        basic_config["lanbdas"] = input.lambdas[0]
        basic_config["thetas"] = input.thetas[0]
        with open(output[0], "w") as outfile:
            json.dump(basic_config, outfile)

rule run_model:
    input:
        "interim/{param_id}_rep{instance_id}_run_config.json"
    output:
        temp("interim/{param_id}_rep{instance_id}_{model}.npy")
    shell:
        """deconvolution --model {wildcards.model} -j {input} outfile={output}"""

rule write_output:
    input:
        expand("interim/{param_id}_rep{instance_id}_{model}.npy", param_id = param_ids, instance_id = np.arange(config["reps"]), model=config["models"])
    output:
        expand("results/{name}_alpha_estimates.tsv", name=config["name"])
    run:
        #open each input
        for filename in input:
            alpha, i = np.load(filename, allow_pickle=True)
            with open(output[0], "a+") as outfile: #dump to file
                outfile.write(filename+"\t"+str(list(alpha))+"\t"+str(int(i))+"\n")
