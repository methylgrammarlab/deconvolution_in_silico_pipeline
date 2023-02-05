

##################################################
# epiread file paths for atlas_from_epiread
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


def unpack_epipaths(epipaths):
    '''
    read config epiread paths
    :param epipaths: config epiread paths
    :return: short name to path, cell type to short name
    '''
    LONG_PATH = {}
    SAMPLES_PER_TYPE = {}
    for cell_type, paths in epipaths.items():
        SAMPLES_PER_TYPE[cell_type] = []
        for path in paths:
            short = path.split("/")[-1].split(".")[0]
            LONG_PATH[short] = path
            SAMPLES_PER_TYPE[cell_type].append(short)
    return LONG_PATH, SAMPLES_PER_TYPE

LONG_PATH = {}
if "atlas_epipaths" in config:
    atlas_long_paths, atlas_samples_per_type = unpack_epipaths(config["atlas_epipaths"])
    LONG_PATH.update(atlas_long_paths)
if "mixture_epipaths" in config:
    mixture_long_paths, mixture_samples_per_type = unpack_epipaths(config["mixture_epipaths"])
    LONG_PATH.update(mixture_long_paths)

def get_samples_per_type(wildcards):
    if wildcards.target == "atlas":
        return atlas_samples_per_type[wildcards.cell_type]
    elif wildcards.target == "mixture":
        return mixture_samples_per_type[wildcards.cell_type]
    else:
        raise ValueError("target wildcard unknown")

rule merge_regions_file: #if not merged, overlapping regions will get duplicated
    input:
        regions="sorted_regions_file.bed"
    output:
        expand("results/{name}_merged_regions_file.bed", name=config["name"])
    shell:
        """bedtools merge -i {input.regions} > {output}"""

rule cut_regions_from_epireads:
    input:
        epiread= lambda wildcards: LONG_PATH[wildcards.sample],
        regions=expand("results/{name}_merged_regions_file.bed", name=config["name"])
    output:
        "interim/small/{cell_type}_sample_{sample}_small_verified.epiread"
    shell:
        """bedtools intersect -u -a {input.epiread} -b {input.regions}  | sort -k1,1 -k2,2n > {output}"""

rule merge_bioreps:
    input:
        lambda wildcards: expand("interim/small/{{cell_type}}_sample_{sample}_small_verified.epiread", sample=get_samples_per_type(wildcards))
    output:
        "interim/{cell_type}_{target}_epipaths.epiread"
    run:
        shell("""sort -m -k1,1 -k2,2n {input} | sed 's/$/\t{wildcards.cell_type}/' > {output}""")
        #adding origin to each row

###############################################################################

rule epiread_to_bedgraph: #only for epiread
    input:
        epiread="interim/{cell_type}_atlas_epipaths.epiread.gz",
        index="interim/{cell_type}_atlas_epipaths.epiread.gz.tbi",
        cpg_file=config["cpg_file"],
        regions=expand("results/{name}_merged_regions_file.bed", name=config["name"])
    output:
        "interim/merged/merged_{cell_type}_from_epiread.bedgraph"
    shell:
        """epireadToBedgraph -A --cpg_coordinates={input.cpg_file} --epiread_files={input.epiread} --outfile={output} -b --genomic_intervals={input.regions}"""


###############################################################################
#only for bedgraph

def get_bedgraph_paths(wildcards):
    return config["bedpaths"][wildcards.cell_type]

rule merge_and_sort_files:
    input:
        get_bedgraph_paths
    output:
        temp('interim/merged/merged_{cell_type}_from_bedgraph.bedgraph')
    shell:
        """sort -m -k1,1 -k2,2n {input} > {output}"""


###############################################################################
#common
def epiread_or_bedgraph(wildcards):
    if len(config["regions_file"]):
        return f'interim/merged/merged_{wildcards.cell_type}_from_epiread.bedgraph'
    else:
        return f'interim/merged/merged_{wildcards.cell_type}_from_bedgraph.bedgraph'

rule filter_excluded_regions: #filter with whitelist, optional
    input:
        sorted_file=epiread_or_bedgraph,
        include_regions=config["include_list"]
    output:
        'interim/filtered/filtered_{cell_type}.bedgraph'
    shell:
        """bedtools intersect -u -a {input.sorted_file} -b {input.include_regions} > {output}"""

def use_include_list(wildcards):
    if len(config["include_list"]): #user-supplied regions to include
        return f'interim/filtered/filtered_{wildcards.cell_type}.bedgraph'
    else: #skip filtering step
        return epiread_or_bedgraph(wildcards)

rule format_merged_file:
    input:
        use_include_list
    output:
        'interim/formatted/formatted_{cell_type}.bedgraph'
    shell:
        """awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4*$5,$5}}'  {input} > {output}"""

rule separate_strands:
    input:
        bedgraph='interim/formatted/formatted_{cell_type}.bedgraph',
        cpgs=config["cpg_file"],
        genome=config["genome"]
    output:
        plus='interim/stranded/plus_strand_{cell_type}.bedgraph',
        minus='interim/stranded/minus_strand_{cell_type}.bedgraph'
    run:
        shell("""bedtools intersect -v -sorted -a {input.bedgraph} -b {input.cpgs} -g {input.genome} > {output.minus} """)
        shell("""bedtools intersect -u -sorted -a {input.bedgraph} -b {input.cpgs} -g {input.genome} > {output.plus} """)

rule shift_minus_strand:
    input:
        minus='interim/stranded/minus_strand_{cell_type}.bedgraph',
        genome=config["genome"]
    output:
        'interim/stranded/shifted_minus_strand_{cell_type}.bedgraph'
    shell:
        """bedtools shift -i {input.minus} -s -1 -g {input.genome} > {output}"""

rule cpg_filter: #only keep data aligned to cpg file
    input:
        bedgraph='interim/stranded/shifted_minus_strand_{cell_type}.bedgraph',
        cpgs = config["cpg_file"],
        genome=config["genome"]
    output:
        'interim/stranded/shifted_filtered_minus_strand_{cell_type}.bedgraph'
    shell:
        """bedtools intersect -u -sorted -a {input.bedgraph} -b {input.cpgs} -g {input.genome} > {output} """

rule merge_sorted_files:
    #for input that is individually sorted
    input:
        minus='interim/stranded/shifted_filtered_minus_strand_{cell_type}.bedgraph',
        plus='interim/stranded/plus_strand_{cell_type}.bedgraph'
    output:
        'interim/collapsed/collapsed_{cell_type}.bedgraph'
    shell:
        """sort -m -k1,1 -k2,2n {input.minus} {input.plus} > {output}"""

rule aggregate_and_merge:
    input:
        'interim/collapsed/collapsed_{cell_type}.bedgraph'
    output:
        'interim/aggregated/aggregated_{cell_type}.bedgraph'
    shell:
        'bedtools merge -d 0 -o sum -c 4,5  -i {input} > {output}'


rule split_atlas_methylation_and_coverage:
    input:
        'interim/aggregated/aggregated_{cell_type}.bedgraph'
    output:
        methylation='interim/methylation/{cell_type}.bedgraph',
        coverage='interim/coverage/{cell_type}.bedgraph'
    run:
        shell("""awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4}}' {input} > {output.methylation}""")
        shell("""awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$5}}' {input} > {output.coverage}""")

rule unite_atlas_samples:
    input:
        methylation=expand('interim/methylation/{cell_type}.bedgraph', cell_type=config["cell_types"]),
        coverage=expand('interim/coverage/{cell_type}.bedgraph', cell_type=config["cell_types"])
    output:
        methylation = expand('interim/combined_methylation_{name}.bedgraph', name=config["name"]),
        coverage = expand('interim/combined_coverage_{name}.bedgraph', name=config["name"])
    run:
        methylation_header = "\t".join(["CHROM", "START", "END"] + [x+"_METH" for x in config["cell_types"]]) + "\n"
        coverage_header = "\t".join(["CHROM", "START", "END"] + [x+"_COV" for x in config["cell_types"]]) + "\n"
        with open(output.methylation[0], 'w') as outfile:
            outfile.write(methylation_header)
        with open(output.coverage[0], 'w') as outfile:
            outfile.write(coverage_header)
        shell("""unionBedGraphs -i {input.methylation} >> {output.methylation}"""),
        shell("""unionBedGraphs -i {input.coverage} >> {output.coverage}"""),

rule create_interim_atlas_files:
    input:
        methylation = expand('interim/combined_methylation_{name}.bedgraph', name=config["name"]),
        coverage = expand('interim/combined_coverage_{name}.bedgraph', name=config["name"])
    output:
        pasted = expand('interim/pasted_methylation_coverage_{name}.bedgraph', name=config["name"])
    run:
        shell(""" paste {input.methylation} {input.coverage} > {output.pasted}""")

rule create_final_atlas:
    input:
        expand('interim/pasted_methylation_coverage_{name}.bedgraph', name=config["name"])
    output:
         expand("results/{name}_atlas.bedgraph", name=config["name"])
    run:
        n = len(config["cell_types"])
        cols = list(range(1,4))+list(chain.from_iterable([(x, x+n+3) for x in range(4,4+n)]))
        cols = ['$'+str(x) for x in cols]
        cols = ",".join(cols)
        shell(""" awk -v FS='\t' -v OFS='\t' '{{print {cols} }}' {input} > {output} """)
