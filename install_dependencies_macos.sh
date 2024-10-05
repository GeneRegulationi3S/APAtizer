#!/bin/sh

echo "\nINSTALLING TOOLS (samtools, fastqc, snakemake and htseq-count)...\n"
brew install samtools fastqc snakemake python3 && pip install HTSeq
echo "\nTOOLS INSTALLED.\n"
