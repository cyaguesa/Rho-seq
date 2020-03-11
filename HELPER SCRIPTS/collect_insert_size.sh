#!/bin/bash

# Collect insert sizes using picard


set -e
if [ $# -eq 0 ]; then
  echo "Usage : ./collect_insert_size.sh bam_file_1 [bam_file_2, ...]"
  exit 1
else
  for i in $*
   do
   echo ""; echo "    bam file : ${i}"; echo ""
   #picard CollectInsertSizeMetrics I=genome/mapped/${i}_accepted_hits.bam O=inner_size/${i}.txt H=inner_size/${i}.pdf W=400
   picard CollectInsertSizeMetrics I=../fastq2/mapped/${i}_accepted_hits.bam O=inner_size/${i}.txt H=inner_size/${i}.pdf W=400  
   echo "    done. Next!"
   done
echo "all done !"
fi
#./scripts/collect_insert_size.sh A560T1 A560T2 A560T3 A560T4
