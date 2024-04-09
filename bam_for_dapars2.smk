FILES = glob_wildcards("TRIMMED_READS/{name}.trim.bam")

# extract the {name} values into a list
NAMES = FILES.name

rule all:
    input:
        # use the extracted name values to build new filenames
        "DaPars2/mapping_wig_location_with_depth.txt",
       # "DaPars_data/",
        expand("WIG/{name}.wig", name=NAMES),
	    expand("DaPars2/{name}.mapping_wig_location_with_depth.txt", name=NAMES),

rule bam_to_wig:
    input:
        "TRIMMED_READS/{name}.trim.bam"
    output:
        "WIG/{name}.wig"
    shell:
        "bedtools genomecov -ibam {input} -bga -split -trackline > {output}"

rule mapped_reads:
    input:
        "TRIMMED_READS/{name}.trim.bam",
        "WIG/{name}.wig"
    output:
        "DaPars2/{name}.mapping_wig_location_with_depth.txt",
    shell:
        """
        test=$(samtools view -c {input[0]})
        echo -e {input[1]}'\t'$test>> {output[0]}
        """
rule mapping_wig_location_with_depth:
    input:
        expand("DaPars2/{name}.mapping_wig_location_with_depth.txt", name=NAMES),
    output:
        "DaPars2/mapping_wig_location_with_depth.txt"
    shell:
        "cat {input} >> {output}"


rule DaPars2:
    input:
        expand("Dapars2_configure_file"),
        expand("chr.txt")
    output:
        "DaPars_data/"
    shell:
        "python DaPars2-master/src/DaPars2_Multi_Sample_Multi_Chr.py {input[0]} {input[1]}"   

