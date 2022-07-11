import numpy as np
import pandas as pd

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

rule atlas_over_tims: #intersect with processed tim file
    input:
        tims=expand("results/{name}_processed_tims.txt", name=config["name"]),
        atlas=expand("results/{name}_atlas.bedgraph", name=config["name"])
    output:
        expand("results/{name}_atlas_over_tims.txt", name=config["name"])
    shell:
        """bedtools intersect -a {input.atlas} -b {input.tims}  -u -header -sorted > {output}"""

rule run_model:
    input:
        mixture=expand("data/{name}/{{param_id}}_rep{{instance_id}}_data.npy", name=config["name"]),
        atlas=expand("data/{name}/{{param_id}}_rep{{instance_id}}_metadata_{{model}}.npy" ,name=config["name"])
    params:
        instance="{instance_id}",
        run_config=lambda wildcards: runs[runs.index == int(wildcards.param_id)].iloc[0,:].to_dict()
    output:
        temp("interim/{param_id}_rep{instance_id}_{model}.npy")
    shell:
        """deconvolution --model {wildcards.model}}""" #todo: add params


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
