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

from bimodal_detector import AtlasEstimator

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
        regions="merged_sorted_regions_file.bed", #should be sorted
        atlas=get_atlas_file
    output:
        expand("results/{name}_atlas_over_regions.txt", name=config["name"])
    shell:
        """bedtools intersect -a {input.atlas} -b {input.regions}  -u -header -sorted > {output}"""


#only for epiread models
rule create_epistate_atlas:
    input:
        regions = "merged_sorted_regions_file.bed"
    output:
        lambdas = expand("results/{name}_lambdas.bedgraph", name=config["name"]),
        thetas = expand("results/{name}_thetas.bedgraph", name=config["name"])
    run:
        epiread_files = []
        labels = []
        for cell_type, v in config["atlas_epipaths"].items():
            for path in v:
                epiread_files.append(path)
                labels.append(cell_type)
        basic_config = {"genomic_intervals":input.regions, "cpg_coordiantes":config["cpg_file"],
                        "cell_types":config["cell_types"],
                  "epiread_files":epiread_files, "labels":labels, "outdir":"results", "name":config["name"],
                "epiformat":config["epiformat"], "header":False, "bedfile":True, "parse_snps": False,
        "get_pp":False, "walk_on_list": False,"verbose" : False}
        runner = AtlasEstimator(config)
        runner.run()