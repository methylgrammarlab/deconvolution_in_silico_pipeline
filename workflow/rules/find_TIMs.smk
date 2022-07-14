
def get_atlas_file(wildcards): #TODO: move/remove
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
        """sbatch workflow/scripts/celfie_scripts/tim.sh -i {input[0]} -o {output.raw} -s {output.summed} -w {params.slop} """+\
        """-n {config[n_tims]} -t {params.t} -d 15 -e 1"""

#sbatch workflow/scripts/celfie_scripts/tim.sh -i results/test_atlas_.bedgraph -o results/test_raw_tims_cl.txt -s results/test_raw_tims_summed_cl.txt -w 500 -n 100 -t 6 -d 15 -e 1

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
        expand("results/{name}_tims.txt", name=config["name"])
    shell:
        """bedtools slop -b {params.slop} -i {input.tims} -g {input.genome} > {output}"""

rule atlas_over_regions: #intersect with processed tim file
    input:
        regions="sorted_regions_file.bed", #should be sorted
        atlas=get_atlas_file
    output:
        expand("results/{name}_atlas_over_regions.txt", name=config["name"])
    shell:
        """bedtools intersect -a {input.atlas} -b {input.regions}  -u -header -sorted > {output}"""