#!/bin/bash

# Prompt the user to select a genome version
echo "Please select the genome version used in the creation of the BAM files:"
echo "1. mm9"
echo "2. mm10"
echo "3. hg19"
echo "4. hg38"
echo "5. dm6"
echo "6. danRer11"
read -p "Enter the number corresponding to your choice: " choice

# Map the user's choice to the genome version
case $choice in
    1)
        echo "Running Snakemake to create inputs for BAM files aligned to mm9..."
        snakemake -s src/scripts/process_mm9_1.smk --cores 4 --wait-for-files
        snakemake -s src/scripts/process_mm9_2.smk --cores 1 --wait-for-files
        ;;
    2)
        echo "Running Snakemake to create inputs for BAM files aligned to mm10..."
        snakemake -s src/scripts/process_mm10_1.smk --cores 4 --wait-for-files
        snakemake -s src/scripts/process_mm10_2.smk --cores 1 --wait-for-files
        ;;
    3)
        echo "Running Snakemake to create inputs for BAM files aligned to hg19..."
        snakemake -s src/scripts/process_hg19_1.smk --cores 4 --wait-for-files
        snakemake -s src/scripts/process_hg19_2.smk --cores 1 --wait-for-files
        ;;
    4)
        echo "Running Snakemake to create inputs for BAM files aligned to hg38..."
        snakemake -s src/scripts/process_hg38_1.smk --cores 4 --wait-for-files
        snakemake -s src/scripts/process_hg38_2.smk --cores 1 --wait-for-files
        ;;
    5)
        echo "Running Snakemake to create inputs for BAM files aligned to dm6..."
        snakemake -s src/scripts/process_dm6_1.smk --cores 4 --wait-for-files
        snakemake -s src/scripts/process_dm6_2.smk --cores 1 --wait-for-files
        ;;
    6)
        echo "Running Snakemake to create inputs for BAM files aligned to danRer11..."
        snakemake -s src/scripts/process_danRer11_1.smk --cores 4 --wait-for-files
        snakemake -s src/scripts/process_danRer11_2.smk --cores 1 --wait-for-files
        ;;
    7)
    	echo "Running Snakemake to create inputs for BAM files aligned to ce11..."
     	snakemake -s src/scripts/process_ce11_1.smk --cores 4 --wait-for-files
      	snakemake -s src/scripts/process_ce11_2.smk --cores 1 --wait-for-files
        ;;
    *)
	echo "Invalid choice. Please run the script again and choose a valid option." exit 1
        ;;
esac
