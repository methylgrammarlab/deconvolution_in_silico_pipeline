
'''
process TIM output files for creating mixtures
and running celfie-plus

'''

rule slop_tims:
    input:
        # tims=expand("results/{name}_tims_summed.txt", name=config["name"]),
        tims=expand("results/{name}_raw_tims.txt", name=config["name"]),#should have no header!
        genome=config["genome"]
    params:
        slop=config["slop"]
    output:
        temp("slopped_tims.bed")
    run:
        if config["slop"] > 0:
            shell("""bedtools slop -b {params.slop} -i {input.tims} -g {input.genome} > {output}""")
        else:
            shell("""cp {input.tims} {output}""")

rule sort_tims:
    input:
        "slopped_tims.bed"
    output:
        temp("sorted_slopped_tims.bed")
    shell:
        """bedtools sort -i {input} > {output}"""

rule merge_tims: #so no CpG is read twice
    input:
        "sorted_slopped_tims.bed"
    output:
        expand("results/{name}_processed_tims.txt", name=config["name"])
    shell:
        """bedtools merge -i {input} > {output}"""

rule atlas_over_tims: #intersect with processed tim file
    input:
        tims=expand("results/{name}_processed_tims.txt", name=config["name"]),
        atlas=expand("results/{name}_atlas.bedgraph", name=config["name"])
    output:
        expand("results/{name}_atlas_over_tims.txt", name=config["name"])
    shell:
        """bedtools intersect -a {input.atlas} -b {input.tims}  -u -header -sorted > {output}"""

rule sum_tims:
    input:
        tims=expand("results/{name}_processed_tims.txt", name=config["name"]),
        atlas=expand("results/{name}_atlas.bedgraph", name=config["name"]),
        genome=config["genome"]
    output:
        expand("results/{name}_tims_summed.txt", name=config["name"])
    run:
        n_cols = len(config["cell_types"])*2
        cols = ",".join([str(x) for x in range(4, n_cols+4)])
        shell("""bedtools map -b {input.atlas} -a {input.tims} -o sum -c {cols} -null 0 -g {input.genome} > {output}""")

# rule generate_tims:
#     input:
#         expand("results/{name}_atlas.bedgraph", name=config["name"])
#     params:
#         t=len(config["cell_types"]),
#         slop=int(config["slop"])*2
#
#     output:
#         raw=expand("results/{name}_raw_tims.txt", name=config["name"]),
#         summed=expand("results/{name}_tims_summed.txt", name=config["name"])
#     shell:
#         """sbatch scripts/tim.sh -i {input} -o {output.raw} -s {output.summed} -w {params.slop} """+\
#         """-n {config[n_tims]} -t {params.t} -d 15 -e 1"""
