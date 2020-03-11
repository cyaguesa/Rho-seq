#!/bin/bash

# Compute D-ratio (ratio between the number of reads stopping at one position and the number of reads covering a position)
# for the spike-in (forward strand only).

set -e
if [ $# -eq 0 ]; then
  echo "Usage : spike-in_D-ratio.sh bam_file_1 [bam_file_2, ...]"
  exit 1
else
  for i in $*
   do
   bedtools genomecov -ibam spike-in/${i}_sorted.bam_fwd.bam -strand + -d -5 | sort -k 1,1 -k2,2n -r > spike-in/cov/${i}.fwd.5end
   samtools sort -l 0 -@ 1 -m 3G -n spike-in/${i}_sorted.bam_fwd.bam | bedtools bamtobed -i stdin -bedpe | awk 'BEGIN {OFS="\t"};{m=$2;M=0;for(i=2;i<=6;i++)if((i != 4)){if(($i<m))m=$i;if(($i>M))M=$i};print $1,m,M}' | sort-bed --max-mem 5G - | bedtools genomecov -i stdin -d -g spike-in/spike_in_rDNA_chromInfo.txt | sort -k 1,1 -k2,2n -r > spike-in/cov/${i}.fwd_fragment.cov
   paste spike-in/cov/${i}.fwd_fragment.cov spike-in/cov/${i}.fwd.5end | awk 'NR==62{print $3,$6,$6/$3}' >> spike-in/cov/${i}.spike_in_cov
   paste spike-in/cov/${i}.fwd_fragment.cov spike-in/cov/${i}.fwd.5end | awk 'BEGIN {OFS="\t"};{print $1,$2,$3,$6,($6+1)/($3+1)}' > spike-in/cov/${i}.fwd_cov_stop_ratio.cov

   done
fi

#./scripts/bam_to_strand_specific_cov.sh spike-in/*sorted.bam
#./scripts/spike_in_D-ratio.sh D045T11_spike_in_reads D045T12_spike_in_reads D045T13_spike_in_reads D045T14_spike_in_reads D045T15_spike_in_reads D045T16_spike_in_reads D045T17_spike_in_reads D045T18_spike_in_reads D045T19_spike_in_reads D045T20_spike_in_reads D045T21_spike_in_reads D045T22_spike_in_reads D045T23_spike_in_reads D045T24_spike_in_reads D045T25_spike_in_reads D045T26_spike_in_reads D045T27_spike_in_reads D045T28_spike_in_reads D045T29_spike_in_reads D045T30_spike_in_reads

