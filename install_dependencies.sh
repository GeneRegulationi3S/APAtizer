#!/bin/sh

echo "\nINSTALLING TOOLS (samtools, fastqc and htseq-count)...\n"
sudo apt install samtools fastqc && pip install HTSeq
echo "\nTOOLS INSTALLED.\n"

