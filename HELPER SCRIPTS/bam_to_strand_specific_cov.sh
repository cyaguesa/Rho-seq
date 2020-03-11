#!/bin/bash
if [ $# -eq 0 ]; then
  echo "Usage : mapped_to_coverage_strand_specific.sh bam_file_1 [bam_file_2, ...]"
else
  for i in $*
    do
    echo ""; echo "bam file : ${i}"; echo ""
    echo "pre-processing forward strand : ..."
    samtools index $i
    samtools view -b -f 128 -F 16 $i > ${i}_fwd1.bam && samtools index ${i}_fwd1.bam
    echo "     1/3..."; samtools view -b -f 80 $i > ${i}_fwd2.bam && samtools index ${i}_fwd2.bam
    echo "     2/3..."; samtools merge -f ${i}_fwd.bam ${i}_fwd1.bam ${i}_fwd2.bam
    samtools index ${i}_fwd.bam ; echo '    3/3' ; echo ' ';    echo "pre-processing reverse strand : ..."
    samtools view -b -f 144 $i > ${i}_rev1.bam && samtools index ${i}_rev1.bam
    echo "     1/3..."; samtools view -b -f 64 -F 16 $i > ${i}_rev2.bam && samtools index ${i}_rev2.bam
    echo "     2/3..."; samtools merge -f ${i}_rev.bam ${i}_rev1.bam ${i}_rev2.bam
    samtools index ${i}_rev.bam ; echo "     3/3." ; echo ""; echo "computing coverage..." ; samtools depth ${i}_rev.bam > ${i}_rev.depth; samtools depth ${i}_fwd.bam > ${i}_fwd.depth ; echo ""; echo " cleaning up..."
    rm ${i}_rev1.bam; rm ${i}_rev2.bam; rm ${i}_fwd1.bam; rm ${i}_fwd2.bam; 
    rm ${i}_rev1.bam.bai; rm ${i}_rev2.bam.bai; rm ${i}_fwd1.bam.bai; rm ${i}_fwd2.bam.bai; echo ""; echo "Done. Have a good day !"
    done
fi

