# APAtizer
A web app for gene expression and alternative polyadenylation analysis of TCGA RNA-Seq data.

# Workflow
![image](https://github.com/brss12/APAtizer/assets/121204829/ea7c789d-907a-42bc-b331-94387a9b4325)

## 1. Pre-processing BAM files
The script [bam_analysis.smk](bam_analysis.smk) employs the snakemake workflow to create fastqc reports, sort and remove the duplicates from the raw TCGA BAM files. The script requires the raw BAM files to be placed in a folder called **RAW_BAM** and, upon running, creates two new folders called **SORTED_BAM** and **TRIMMED_READS**. In the **SORTED_BAM** folder, is where the sorted bam files along with fastqc reports are placed and in the **TRIMMED_READS** is where the de-duplicated bam files along with fastqc reports are placed.

Depending on the size and ammount of BAM files, we recommend running the script in a HPC environment. Run the following script on the same directory of the **RAW_BAM** folder. 
```shell
snakemake -s bam_analysis.smk --cores [number_of_cores] --wait-for-files
```
Indexing the de-duplicated BAM files from the previous step is a good practice. The script [index_script.sh](index_script.sh) can be executed by running the following command and specifying as argument the number of threads to use.
```shell
./index_script.sh [number_of_threads]
```
## 2. Creating the DaPars2 txt files
Alternative polyadenylation analysis in APAtizer can be done by using one of two algorithms, APAlyzer and DaPars2. To use the APAlyzer algorithm, the de-duplicated bam files obtained in the previous step are used as input. On the other hand, for the DaPars2 algorithm, the user needs to run the first part of the DaPars2 analysis in the terminal and APAtizer only deals with the final manipulations and calculations to obtain the final gene lists.

The script [dapars2.smk](dapars2_files/dapars2.smk) runs the first part of the DaPars2 analysis. The script must be executed in the same directory of the **TRIMMED_READS** folder. The following command runs the script.
```shell
snakemake -s dapars2.smk --cores [number_of_cores] --wait-for-files
```
After this, run the following command to make a folder called **DaPars_data_final** and move all of the output files to it.
```shell
mkdir DaPars_data_final/ && mv DaPars_data_chr*/*temp* DaPars_data_final/
```
## 3. Creating htseq files for DGE analysis
The htseq files are used in the differential gene expression (DGE) analysis in APAtizer. These files are created via htseq-count by running the script [htseq_script.sh](htseq_script.sh). Also, those files are filtered and manipulated by some extra commands as listed below.
```shell
./htseq_script.sh && mkdir --parents TRIMMED_htseq/FILTERED && mv TRIMMED_READS/*htseq* TRIMMED_htseq/ && mv filter_folder/ TRIMMED_htseq/ && Rscript protein_coding_filter.R
```
