#!/bin/bash

threads=$1

for file in $(ls TRIMMED_READS/*.bam);
do
  samtools index -@ $threads $file;
done

echo "index run ended"
