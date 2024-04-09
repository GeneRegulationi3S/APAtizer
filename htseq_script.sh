#!/bin/bash

for file in $(ls TRIMMED_READS/*.bam); do
  htseq_file="${file%.bam}.htseq.txt"
  if [ ! -e "$htseq_file" ]; then
    htseq-count -f bam -r name -s no "$file" hg38.ncbiRefSeq.gtf > "$htseq_file"
    echo "Processed $file"
  fi
done

echo "htseq run ended"
