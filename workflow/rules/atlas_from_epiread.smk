'''
create atlas for celfie from "methylation_extract" bedgraphs
takes in bedgraphs formattted chrom start end %meth coverage
final is chrom start end A_meth A_cov B_meth B_cov

'''

from itertools import chain

rule epiread_to_bedgraph:
    input:
        epiread="interim/{cell_type}.epiread.gz",
        index="interim/{cell_type}.epiread.gz.tbi",
        cpg_file=config["cpg_file"],
        regions=expand("results/{name}_processed_tims.txt", name=config["name"]), #TODO: verify
    output:
        'interim/merged/merged_{cell_type}.bedgraph'
    shell:
        """python3 scripts/epireadToBedgraph.py {input.cpg_file} {input.epiread} {output} -b {input.regions}"""

rule sort_files:
    input:
        'interim/merged/merged_{cell_type}.bedgraph'
    output:
        'interim/sorted/sorted_{cell_type}.bedgraph'
    shell:
        """sort -k1,1 -k2,2n {input} > {output}"""

rule filter_excluded_regions: #filter with whitelist
    input:
        sorted_file='interim/sorted/sorted_{cell_type}.bedgraph',
        include_regions=config["include_list"]
    output:
        'interim/filtered/filtered_{cell_type}.bedgraph'
    shell:
        """bedtools intersect -u -a {input.sorted_file} -b {input.include_regions} > {output}"""

rule format_merged_file:
    input:
        'interim/filtered/filtered_{cell_type}.bedgraph'
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

