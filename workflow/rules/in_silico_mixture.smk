'''
For each TIM, given alpha, sample from epiread files (tabix the region, randomly select reads)
    - tabix out all TIMs from epiread, to have small files
    - bedops --everything epireads of each cell type to one medium file
    - calculate number of reads to get from each cell
    - bedtools sample -i input.bed -n 1000
    - merge to new epiread file
'''
#from TIM output need to remove header and slop tims
import json
import numpy as np
import pandas as pd

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index

LONG_PATH = {}
SAMPLES_PER_TYPE = {}
for cell_type, paths in config["epipaths"].items():
    SAMPLES_PER_TYPE[cell_type] = []
    for path in paths:
        short = path.split("/")[-1].split(".")[0]
        LONG_PATH[short] = path
        SAMPLES_PER_TYPE[cell_type].append(short)

def get_samples_per_type(wildcards):
    return ["small/" + wildcards.cell_type +"_"+ x + "_small_verified.epiread" for x in SAMPLES_PER_TYPE[wildcards.cell_type]]

def count_lines(wildcards):
    with open(wildcards) as foo:
        lines = len(foo.readlines())
    return lines


rule measure_length:
    input:
        epireads=expand("interim/{cell_type}.epiread", cell_type=config["cell_types"])
    output:
        expand("interim/{name}_mean_cpgs_per_read.json", name=config["name"])
    script:
        "../scripts/measure_read_length.py"

rule calculate_reads: #number to sample from each cell type
    input:
        cpgs=expand("results/{name}_atlas_over_tims.txt", name=config["name"]), #number of cpgs to cover
        lengths=expand("interim/{name}_mean_cpgs_per_read.json", name=config["name"])
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
        tims=expand("results/{name}_processed_tims.txt", name=config["name"])
    output:
        "small/{cell_type}_{sample}_small_verified.epiread"
    shell:
        """bedtools intersect -u -a {input.epiread} -b {input.tims}  | sort -k1,1 -k2,2n > {output}"""

rule merge_bioreps:
    input:
        get_samples_per_type
    output:
        "interim/{cell_type}.epiread"
    run:
        shell("""sort -m -k1,1 -k2,2n {input} | sed 's/$/\t{wildcards.cell_type}/' > {output}""")
        #adding source to each row

rule sample_from_epiread:
    input:
        epiread="interim/{cell_type}.epiread",
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

