FILES = glob_wildcards('../RAW_BAM/{name}.bam')

# extract the {name} values into a list
NAMES = FILES.name

rule all:
    input:
        # use the extracted name values to build new filenames
        expand("../SORTED_BAM/{name}.sorted.bam", name=NAMES),
        expand("../SORTED_QC/{name}.sorted_fastqc.zip", name=NAMES),
        expand("../SORTED_QC/{name}.sorted_fastqc.html", name=NAMES),
        expand("../TRIMMED_READS/{name}.trim.bam", name=NAMES),
        expand("../TRIMMED_READS/{name}.trim.bam.bai", name=NAMES),
        expand("../TRIMMED_QC/{name}.trim_fastqc.zip", name=NAMES),
        expand("../TRIMMED_QC/{name}.trim_fastqc.html", name=NAMES),
        expand("../TRIMMED_htseq/{name}.trim.htseq.txt", name=NAMES),
        expand("../WIG/{name}.wig", name=NAMES),
        expand("../DaPars2/{name}.mapping_wig_location_with_depth.txt", name=NAMES)

rule samtools_sort:
    priority: 20
    input:
        "../RAW_BAM/{name}.bam",
    output:
        "../SORTED_BAM/{name}.sorted.bam",
        zip = "../SORTED_QC/{name}.sorted_fastqc.zip",
        html = "../SORTED_QC/{name}.sorted_fastqc.html"
    params:
        path = "../SORTED_QC/"
    shell:
        "samtools sort -o {output[0]} {input} -T ../SORTED_BAM/ && fastqc {output[0]} -o {params.path}"

rule picard:
    priority: 18
    input: rules.samtools_sort.output[0]
    output:
        "../TRIMMED_READS/{name}.trim.bam",
    shell:
        "gatk MarkDuplicates --REMOVE_DUPLICATES true -I {input} -O {output} -M ../marked_dup_metrics.txt -TMP_DIR ../TRIMMED_READS/"

rule htseq_count:
    priority: 14
    input: rules.picard.output,
    output:
        "../TRIMMED_READS/{name}.trim.bam.bai",
        "../TRIMMED_htseq/{name}.trim.htseq.txt",
        zip = "../TRIMMED_QC/{name}.trim_fastqc.zip",
        html = "../TRIMMED_QC/{name}.trim_fastqc.html"
    params:
        path = "../TRIMMED_QC/"
    shell:
        "fastqc {input} -o {params.path} && samtools index {input} && htseq-count -f bam -r name -s no {input} src/annotations/mm10.ncbiRefSeq.gtf.gz > {output[1]}"

rule bam_to_wig:
    priority: 12
    input:
        "../TRIMMED_READS/{name}.trim.bam"
    output:
        "../WIG/{name}.wig"
    shell:
        "bedtools genomecov -ibam {input} -bga -split -trackline > {output}"

rule mapped_reads:
    priority: 11
    input:
        "../TRIMMED_READS/{name}.trim.bam",
        "../WIG/{name}.wig"
    output:
        "../DaPars2/{name}.mapping_wig_location_with_depth.txt",
    shell:
        """
        test=$(samtools view -c {input[0]})
        echo -e {input[1]}'\t'$test >> {output}
        """
