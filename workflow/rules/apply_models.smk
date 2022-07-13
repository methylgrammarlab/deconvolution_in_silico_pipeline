import numpy as np
import pandas as pd


runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

def get_atlas_file(wildcards):
    if len(config["atlas_file"]):
        return config["atlas_file"]
    else: #no user supplied atlas
        return expand("results/{name}_atlas.bedgraph", name=config["name"])

def get_regions_file(wildcards):
    if len(config["regions_file"]):
        return config["regions_file"]
    else:
        return expand("results/{name}_tims.txt", name=config["name"])

rule sort_regions_file:
    input:
         regions=get_regions_file
    output:
        temp("sorted_regions_file.bed")
    shell:
        """sort -k1,1 -k2,2n {input.regions} > {output}"""

rule atlas_over_regions: #intersect with processed tim file
    input:
        regions="sorted_regions_file.bed", #should be sorted
        atlas=get_atlas_file
    output:
        expand("results/{name}_atlas_over_regions.txt", name=config["name"])
    shell:
        """bedtools intersect -a {input.atlas} -b {input.regions}  -u -header -sorted > {output}"""

rule run_model:
    input:
        mixture=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz", name=config["name"]),
        index=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz.tbi", name=config["name"]),
        atlas=get_atlas_file,
        regions=get_regions_file #not merged
    params:
        instance="{instance_id}",
        run_config=lambda wildcards: runs[runs.index == int(wildcards.param_id)].iloc[0,:].to_dict()
    output:
        temp("interim/{param_id}_rep{instance_id}_{model}.npy")
    shell:
        """deconvolution --model {wildcards.model} {run_config}""" #todo: make sure config works like this


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
