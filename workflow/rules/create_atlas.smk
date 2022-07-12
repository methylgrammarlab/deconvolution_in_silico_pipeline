

rule epiread_to_bedgraph: #only for epiread
    input:
        epiread="interim/{cell_type}.epiread.gz", #from in_silico_mixture.smk
        index="interim/{cell_type}.epiread.gz.tbi",
        cpg_file=config["cpg_file"],
        regions=config["user_regions"]
    output:
        'interim/merged/merged_{cell_type}_from_epiread.bedgraph'
    shell:
        """epireadToBedgraph --cpg_coordinates {input.cpg_file} --epireads"""+\
        """ {input.epiread} --outfile {output} -b {input.regions}"""

rule sort_files: #only for epiread
    input:
        'interim/merged/merged_{cell_type}_{method}.bedgraph'
    output:
        'interim/sorted/sorted_{cell_type}_{method}.bedgraph'
    shell:
        """sort -k1,1 -k2,2n {input} > {output}"""


def get_bedgraph_paths(wildcards): #only for bedgraph
    return config["bedpaths"][wildcards.cell_type]

rule merge_and_sort_files: #only for bedgraph
    input:
        get_bedgraph_paths
    output:
        temp('interim/sorted/sorted_{cell_type}_from_bedgraph.bedgraph')
    shell:
        """sort -m -k1,1 -k2,2n {input} > {output}"""


###############################################################################
#common
def epiread_or_bedgraph(wildcards):
    if config["regions_file"]:
        return f'interim/sorted/sorted_{wildcards.cell_type}_from_epiread.bedgraph'
    else:
        return f'interim/sorted/sorted_{wildcards.cell_type}_from_bedgraph.bedgraph'

rule filter_excluded_regions: #filter with whitelist, optional
    input:
        sorted_file=epiread_or_bedgraph,
        include_regions=config["include_list"]
    output:
        'interim/filtered/filtered_{cell_type}.bedgraph'
    shell:
        """bedtools intersect -u -a {input.sorted_file} -b {input.include_regions} > {output}"""

def use_include_list(wildcards):
    if config["include_list"]: #user-supplied regions to include
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


rule split_methylation_and_coverage:
    input:
        'interim/aggregated/aggregated_{cell_type}.bedgraph'
    output:
        methylation='interim/methylation/{cell_type}.bedgraph',
        coverage='interim/coverage/{cell_type}.bedgraph'
    run:
        shell("""awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4}}' {input} > {output.methylation}""")
        shell("""awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$5}}' {input} > {output.coverage}""")

rule unite_samples:
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

rule create_interim_files:
    input:
        methylation = expand('interim/combined_methylation_{name}.bedgraph', name=config["name"]),
        coverage = expand('interim/combined_coverage_{name}.bedgraph', name=config["name"])
    output:
        pasted = expand('interim/pasted_methylation_coverage_{name}.bedgraph', name=config["name"])
    run:
        shell(""" paste {input.methylation} {input.coverage} > {output.pasted}""")

rule stitch_samples:
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
