
import json
import numpy as np
import pandas as pd

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

# epiread file paths for mixture. can be different from atlas (e.g. holdout samples)
LONG_PATH = {}
SAMPLES_PER_TYPE = {}
for cell_type, paths in config["mixture_epipaths"].items():
    SAMPLES_PER_TYPE[cell_type] = []
    for path in paths:
        short = path.split("/")[-1].split(".")[0]
        LONG_PATH[short] = path
        SAMPLES_PER_TYPE[cell_type].append(short)

def get_atlas_file(wildcards): #TODO: change to import from one source
    if len(config["atlas_file"]):
        return config["atlas_file"]
    else: #no user supplied atlas
        return expand("results/{name}_atlas.bedgraph", name=config["name"])

def get_samples_per_type(wildcards):
    return ["small/" + wildcards.cell_type +"_"+ x + "_small_verified.epiread" for x in SAMPLES_PER_TYPE[wildcards.cell_type]]

def count_lines(wildcards):
    with open(wildcards) as foo:
        lines = len(foo.readlines())
    return lines

rule measure_length: #TODO: replace with something better?
    input:
        epireads=expand("interim/{cell_type}.epiread", cell_type=config["cell_types"])
    output:
        expand("interim/{name}_mean_cpgs_per_read.json", name=config["name"])
    script:
        "scripts/measure_read_length.py"

def get_regions_file(wildcards):
    if config["regions_file"]:
        return config["regions_file"]
    else:
        return expand("results/{name}_tims.txt", name=config["name"])

rule atlas_over_tims: #intersect with processed tim file
    input:
        regions=get_regions_file,
        atlas=get_atlas_file
    output:
        expand("results/{name}_atlas_over_regions.txt", name=config["name"])
    shell:
        """bedtools intersect -a {input.atlas} -b {input.regions}  -u -header -sorted > {output}"""

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
        shell("""python3 scripts/mixture_simulator.py {input.lengths} {params.coverage}"""+\
              """ {n_cpgs} '{params.alpha}' '{params.cell_types}' {output}""")
    #this works: python3 mixture_simulator.py 10 584 '[0.04761905, 0.0952381 ,
    # 0.14285714, 0.19047619, 0.23809524, 0.28571429]' '[Alpha, Beta, Delta, Duct, Acinar, Endothel]' try.json

rule assert_intersection: #make sure tabix and bedtools agree
    input:
        epiread= lambda wildcards: LONG_PATH[wildcards.sample],
        regions=get_regions_file
    output:
        "small/{cell_type}_{sample}_small_verified.epiread"
    shell:
        """bedtools intersect -u -a {input.epiread} -b {input.regions}  | sort -k1,1 -k2,2n > {output}"""

rule merge_bioreps:
    input:
        get_samples_per_type
    output:
        "interim/{cell_type}_mixture_epipaths.epiread"
    run:
        shell("""sort -m -k1,1 -k2,2n {input} | sed 's/$/\t{wildcards.cell_type}/' > {output}""")
        #adding source to each row

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

