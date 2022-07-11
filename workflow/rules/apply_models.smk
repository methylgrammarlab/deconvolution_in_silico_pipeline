import numpy as np
import pandas as pd

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

rule parse_mixture_for_celfie:
    input:
        mixture=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz", name=config["name"]),
        index=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz.tbi", name=config["name"]),
        tims=expand("results/{name}_tims_summed.txt", name=config["name"]),
        cpg_coords=config["cpg_file"]
    output:
        data=expand("data/{name}_{{param_id}}_rep{{instance_id}}_celfie_data.npy", name=config["name"]),
        cov=expand("cov/{name}_{{param_id}}_rep{{instance_id}}_celfie_coverage.json", name=config["name"])
    shell:
        """python3 scripts/mixture_epiparser.py {input.tims} {input.mixture} {input.cpg_coords} {output.data}"""+\
        """ {output.cov}  --celfie"""


rule parse_mixture_for_celfie_plus:
    input:
        mixture=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz", name=config["name"]),
        index=expand("results/mixtures/{name}_{{param_id}}_rep{{instance_id}}_mixture.epiread.gz.tbi", name=config["name"]),
        tims=expand("results/{name}_processed_tims.txt", name=config["name"]),
        cpg_coords=config["cpg_file"]
    output:
        data=expand("data/{name}_{{param_id}}_rep{{instance_id}}_celfie_plus_data.npy", name=config["name"]),
        cov=expand("cov/{name}_{{param_id}}_rep{{instance_id}}_celfie-plus_coverage.json", name=config["name"]),
    shell:
        """python3 scripts/mixture_epiparser.py {input.tims} {input.mixture} {input.cpg_coords} {output.data}"""+\
        """ {output.cov}"""

rule merge_coverage:
    input:
        expand("cov/{name}_{param_id}_rep{instance_id}_{model}_coverage.json", name=config["name"],
        param_id = param_ids, instance_id = np.arange(config["reps"]), model=config["models"])
    output:
        expand("results/mixtures/{name}_mixture_coverages.tsv", name=config["name"])
    run:
        # open each input
        for filename in input:
            with open(filename, "r") as infile:
                cov = json.load(infile)
            with open(output[0], "a+") as outfile:  # dump to file
                outfile.write(filename + "\t" + str(cov) + "\n")

rule parse_celfie_plus_atlas: #beta values-based
    input:
        tims = expand("results/{name}_processed_tims.txt", name=config["name"]),
        atlas = expand("results/{name}_atlas_over_tims.txt", name=config["name"]),
        cpg_coords=config["cpg_file"]
    output:
        expand("results/{name}_celfie_plus_atlas.bedgraph", name=config["name"])
    shell:
        """python3 scripts/atlas_epiparser.py {input.tims} {input.atlas} {input.cpg_coords} {output}"""


rule run_celfie:
    input:
        data=expand("data/{name}_{{param_id}}_rep{{instance_id}}_celfie_data.npy", name=config["name"]),
        metadata=expand("results/{name}_tims_summed.txt", name=config["name"]) #atlas from tim output
    params:
        instance = "{instance_id}",
        num_iterations = lambda wildcards: int(runs.num_iterations[runs.index == int(wildcards.param_id)].values),
        random_restarts = lambda wildcards: int(runs.random_restarts[runs.index == int(wildcards.param_id)].values)
    log:
        "logs/{param_id}_rep{instance_id}_celfie.log"
    output:
        "interim/{param_id}_rep{instance_id}_celfie.npy"
    shell:
        """python3 scripts/celfie_step.py {input.data} {input.metadata}  """ + \
        """ {params.num_iterations} {params.random_restarts} {output}"""



rule run_celfie_plus:
    input:
        data=expand("data/{name}_{{param_id}}_rep{{instance_id}}_celfie_plus_data.npy", name=config["name"]),
        metadata=expand("results/{name}_celfie_plus_atlas.bedgraph", name=config["name"])
    params:
        instance = "{instance_id}",
        num_iterations = lambda wildcards: int(runs.num_iterations[runs.index == int(wildcards.param_id)].values),
        random_restarts = lambda wildcards: int(runs.random_restarts[runs.index == int(wildcards.param_id)].values)
    log:
        "logs/{param_id}_rep{instance_id}_celfie_plus.log"
    output:
        "interim/{param_id}_rep{instance_id}_celfie-plus.npy"
    shell:
        """python3 scripts/celfie_plus.py {input.data} {input.metadata}  """ + \
        """ {params.num_iterations} {params.random_restarts} {output}"""

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
