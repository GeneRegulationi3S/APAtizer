FILES = glob_wildcards('RAW_BAM/{name}.bam')

# extract the {name} values into a list
NAMES = FILES.name

rule all:
    input:
        # use the extracted name values to build new filenames
        expand("SORTED_BAM/{name}.sorted.bam", name=NAMES),
        expand("SORTED_QC/{name}.sorted_fastqc.zip", name=NAMES),
        expand("SORTED_QC/{name}.sorted_fastqc.html", name=NAMES),
        expand("TRIMMED_READS/{name}.trim.bam", name=NAMES),
        expand("TRIMMED_READS/{name}.trim.bam.bai", name=NAMES),
        expand("TRIMMED_QC/{name}.trim_fastqc.zip", name=NAMES),
        expand("TRIMMED_QC/{name}.trim_fastqc.html", name=NAMES),
        expand("TRIMMED_htseq/{name}.trim.htseq.txt", name=NAMES),
        "DaPars2/mapping_wig_location_with_depth.txt",
        "DaPars_data/",
        "TRIMMED_htseq/FILTERED/",
        expand("WIG/{name}.wig", name=NAMES),
        expand("DaPars2/{name}.mapping_wig_location_with_depth.txt", name=NAMES)

rule samtools_sort:
    priority: 20
    input:
        "RAW_BAM/{name}.bam",
    output:
        "SORTED_BAM/{name}.sorted.bam",
    shell:
        "samtools sort -@ 4 -o {output} {input} -T SORTED_BAM/"

rule fastqc:
    priority: 19
    input: rules.samtools_sort.output,
    output:
        zip = "SORTED_QC/{name}.sorted_fastqc.zip",
        html = "SORTED_QC/{name}.sorted_fastqc.html",
    params:
        path = "SORTED_QC/",
    shell:
        "fastqc -t 4 {input} -o {params.path}"

rule picard:
    priority: 18
    input: rules.samtools_sort.output
    output:
        "TRIMMED_READS/{name}.trim.bam",
    shell:
        "gatk MarkDuplicates --REMOVE_DUPLICATES true -I {input} -O {output} -M marked_dup_metrics.txt -TMP_DIR TRIMMED_READS/"

rule trimmed_fastqc:
    priority: 17
    input: rules.picard.output,
    output:
        zip = "TRIMMED_QC/{name}.trim_fastqc.zip",
        html = "TRIMMED_QC/{name}.trim_fastqc.html",
    params:
        path = "TRIMMED_QC/",
    shell:
        "fastqc {input} -o {params.path}"

rule index_bam:
    priority: 16
    input: rules.picard.output,
    output:
        "TRIMMED_READS/{name}.trim.bam.bai",
    shell:
        "samtools index -@ 4 {input}"

rule htseq_count:
    priority: 15
    input: rules.picard.output,
    output:
        "TRIMMED_htseq/{name}.trim.htseq.txt",
    shell:
        "htseq-count -f bam -r name -s no {input} src/hg38.ncbiRefSeq.gtf > {output}"

rule protein_coding_filter:
    priority: 14
    output:
        directory("TRIMMED_htseq/FILTERED")
    shell:
        "Rscript src/protein_coding_filter.R"

rule bam_to_wig:
    priority: 13
    input:
        "TRIMMED_READS/{name}.trim.bam"
    output:
        "WIG/{name}.wig"
    shell:
        "bedtools genomecov -ibam {input} -bga -split -trackline > {output}"

rule mapped_reads:
    priority: 12
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
    priority: 11
    input:
        expand("DaPars2/{name}.mapping_wig_location_with_depth.txt", name=NAMES),
    output:
        "DaPars2/mapping_wig_location_with_depth.txt"
    shell:
        "cat {input} >> {output}"

rule DaPars2:
    priority: 10
    input:
        expand("src/Dapars2_configure_file"),
        expand("src/chr.txt")
    output:
        directory("DaPars_data")
    shell:
        """
        sed -i "7s#Aligned_Wig_files=.*#&$(ls WIG/*.wig | tr '\n' ',' | sed 's/,$//')#" {input[0]} && \
        mkdir DaPars_data && \
        python3 src/DaPars2_Multi_Sample_Multi_Chr.py {input[0]} {input[1]} && \
        mv DaPars_data_chr*/*temp* DaPars_data/ && \
        rm -r DaPars_data_chr* && \
        sed -i '7s/.*/Aligned_Wig_files=/' {input[0]}
        """
