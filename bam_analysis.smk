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
        expand("TRIMMED_QC/{name}.trim_fastqc.zip", name=NAMES),
        expand("TRIMMED_QC/{name}.trim_fastqc.html", name=NAMES)

rule samtools_sort:
    input:
        "RAW_BAM/{name}.bam",
    output:
        "SORTED_BAM/{name}.sorted.bam",
    shell:
       "samtools sort -o {output} {input} -T /YOUR_WORKING_DIRECTORY/SORTED_BAM/"

rule fastqc:
    input:rules.samtools_sort.output
    output:
        zip = "SORTED_QC/{name}.sorted_fastqc.zip",
        html = "SORTED_QC/{name}.sorted_fastqc.html",
    params:
        path="SORTED_QC/",
    shell:
        "fastqc {input} -o {params.path}"    

rule picard:
    input:rules.samtools_sort.output,
    output:"TRIMMED_READS/{name}.trim.bam",
    shell:
        "gatk MarkDuplicates --REMOVE_DUPLICATES true -I {input} -O {output} -M marked_dup_metrics.txt -TMP_DIR /YOUR_WORKING_DIRECTORY/TRIMMED_READS/"

rule trimmed_fastqc:
    input: rules.picard.output,
    output:
        zip = "TRIMMED_QC/{name}.trim_fastqc.zip",
        html = "TRIMMED_QC/{name}.trim_fastqc.html",
    params:
        path="TRIMMED_QC/",
    shell:
        "fastqc {input} -o {params.path}"
