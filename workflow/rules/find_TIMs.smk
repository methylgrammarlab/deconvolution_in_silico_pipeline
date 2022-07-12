
def get_atlas_file(wildcards):
    if len(config["atlas_file"]):
        return config["atlas_file"]
    else: #no user supplied atlas
        return expand("results/{name}_atlas.bedgraph", name=config["name"])

rule generate_tims:
    input:
        get_atlas_file
    params:
        t=len(config["cell_types"]),
        slop=int(config["slop"])*2

    output:
        raw=temp(expand("results/{name}_raw_tims.txt", name=config["name"])),
        summed=temp(expand("results/{name}_tims_summed.txt", name=config["name"]))

    shell:
        """sbatch scripts/celfie_scripts/tim.sh -i {input} -o {output.raw} -s {output.summed} -w {params.slop} """+\
        """-n {config[n_tims]} -t {params.t} -d 15 -e 1"""

rule remove_header:
    input:
        tim_file = expand("results/{name}_raw_tims.txt", name=config["name"])
    output:
        temp("tims_no_header.bed")
    run:
        shell("tail -n +2 {input.tim_file} | cut -f 1,2,3 > {output}")

rule slop_tims:
    input:
        tims="tims_no_header.bed", #make sure this has no header
        genome=config["genome"]
    params:
        slop=config["slop"]
    output:
        temp("slopped_tims.bed")
    shell:
        """bedtools slop -b {params.slop} -i {input.tims} -g {input.genome} > {output}"""

rule sort_tims:
    input:
        "slopped_tims.bed"
    output:
        expand("results/{name}_tims.txt", name=config["name"])
    shell:
        """bedtools sort -i {input} > {output}"""


