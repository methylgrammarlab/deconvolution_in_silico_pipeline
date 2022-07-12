import numpy as np
import pandas as pd


runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

def get_atlas_file(wildcards): #TODO: change to import from one source
    if len(config["atlas_file"]):
        return config["atlas_file"]
    else: #no user supplied atlas
        return expand("results/{name}_atlas.bedgraph", name=config["name"])

rule run_model:
    input:
        mixture=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz", name=config["name"]),
        index=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz.tbi", name=config["name"]),
        atlas=get_atlas_file
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
