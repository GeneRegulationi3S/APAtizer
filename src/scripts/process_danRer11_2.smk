FILES = glob_wildcards('../RAW_BAM/{name}.bam')

# extract the {name} values into a list
NAMES = FILES.name

rule all:
    input:
        # use the extracted name values to build new filenames
        "../DaPars2/mapping_wig_location_with_depth.txt",
        "../DaPars_data/",
        "../TRIMMED_htseq/FILTERED/"

rule mapping_wig_location_with_depth:
    priority: 100
    input:
        expand("../DaPars2/{name}.mapping_wig_location_with_depth.txt", name=NAMES),
    output:
        "../DaPars2/mapping_wig_location_with_depth.txt",
    shell:
        "cat {input} >> {output}"

rule DaPars2:
    priority: 20
    input:
        expand("src/dapars_configs/Dapars2_configure_file_danRer11"),
        expand("src/dapars_configs/chr_danRer11.txt"),
    output:
        directory("../DaPars_data")
    shell:
        """
        sed -i "7s#Aligned_Wig_files=.*#&$(ls ../WIG/*.wig | tr '\n' ',' | sed 's/,$//')#" {input[0]} && \
        mkdir ../DaPars_data && \
        python3 src/dapars_configs/DaPars2_Multi_Sample_Multi_Chr.py {input[0]} {input[1]} && \
        mv ../DaPars_data_chr*/*temp* ../DaPars_data/ && \
        sed -i '7s/.*/Aligned_Wig_files=/' {input[0]}
        """

rule protein_coding_filter:
    priority: 10
    output:
        directory("../TRIMMED_htseq/FILTERED")
    shell:
        "Rscript src/filter_files/protein_coding_filter_danRer11.R"
