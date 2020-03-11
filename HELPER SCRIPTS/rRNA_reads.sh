#!/bin/bash

# get rRNA reads info

set -e
if [ $# -eq 0 ]; then
  echo "Usage : ./rRNA_reads.sh bam_file_1 [bam_file_2, ...]"
  exit 1
else
  for i in $*
   do
   echo ""; echo "    bam file : ${i}"; echo ""
   samtools idxstats $i | cut -f 3 | tail -n +2 | awk '{s+=$1}END{print s}'
   echo "    done. Next!"
   done
echo "all done !"
fi
#./scripts/rRNA_reads.sh spike-in/*spike_in_reads_sorted.bam_fwd.bam ~/Desktop/these_Carlo/R_terminal/D-seq/fastq2/mapped/*spike_in_reads_sorted.bam.fwd.bam