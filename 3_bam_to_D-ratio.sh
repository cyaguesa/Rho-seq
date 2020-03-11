# !/bin/bash

# Compute the ratio between fragment stopping at a site and fragment coverage.

# This code requires samtools, rscript and bedtools installed and in your PATH.
# The scripts splitbed_to_splitcov.R and bedpe_helper.R should be executable (chmod +x).
# when used on spike-in/rDNA, chromInfo = spike_in_rDNA_chromInfo.txt, otherwise,  chromInfo = chrom_summary.txt

cd /Users/carlo/Desktop/these_Carlo/R_terminal/D-seq/2018_D_seq_pombe_single_mutants/genome/
chromInfo="chrom_summary.txt" # path the chromosome information
splitbed_to_splitcov="../scripts/splitbed_to_splitcov.R" #path to helper script 1
bedpe_helper="../scripts/bedpe_helper.R" #path to helper script 2

set -e
if [ $# -eq 0 ]; then
  echo "Usage : ban_to_D-ratio.sh bam_file_1 [bam_file_2, ...]"
  exit 1
else
  for i in $*
   do
   echo ""; echo "bam file : ${i}"; echo ""
   echo "split bam : ..."
   samtools view -b -f 82 -F 256 mapped/${i} > mapped/${i}.fwd1.bam     # 1. alignments of the first in pair if they map to the reverse strand + properly paired + not secondary alignment
   samtools view -b -f 130 -F 272 mapped/${i} > mapped/${i}.fwd2.bam    # 2. alignments of the second in pair if they map to the forward strand + properly paired + not secondary alignment
   samtools view -b -f 66 -F 272 mapped/${i} > mapped/${i}.rev1.bam    # 3. alignments of the first in pair if they map to the forward strand + properly paired + not secondary alignment
   samtools view -b -f 146 -F 256 mapped/${i} > mapped/${i}.rev2.bam     # 4. alignments of the second in pair if they map to the reverse strand + properly paired + not secondary alignment
   samtools index mapped/${i}.fwd1.bam
   samtools index mapped/${i}.fwd2.bam
   samtools index mapped/${i}.rev1.bam
   samtools index mapped/${i}.rev2.bam

   echo "done ... compute 5'end : ..."
   bedtools genomecov -ibam mapped/${i}.fwd1.bam -d -5 | sort -k 1,1 -k2,2n -r > cov/${i}.fwd1.5end
   bedtools genomecov -ibam mapped/${i}.fwd2.bam -d -5 | sort -k 1,1 -k2,2n -r > cov/${i}.fwd2.5end
   bedtools genomecov -ibam mapped/${i}.rev1.bam -d -5 | sort -k 1,1 -k2,2n -r > cov/${i}.rev1.5end
   bedtools genomecov -ibam mapped/${i}.rev2.bam -d -5 | sort -k 1,1 -k2,2n -r > cov/${i}.rev2.5end

   echo "done ... compute extended coverage : ..."
   samtools merge mapped/${i}.fwd_.bam mapped/${i}.fwd1.bam mapped/${i}.fwd2.bam
   samtools sort mapped/${i}.fwd_.bam > mapped/${i}.fwd.bam
   rm mapped/${i}.fwd_.bam
   samtools index mapped/${i}.fwd.bam
   samtools merge mapped/${i}.rev_.bam mapped/${i}.rev1.bam mapped/${i}.rev2.bam
   samtools sort mapped/${i}.rev_.bam > mapped/${i}.rev.bam
   rm mapped/${i}.rev_.bam
   samtools index mapped/${i}.rev.bam
   bedtools genomecov -ibam mapped/${i}.fwd.bam -d -split | sort -k 1,1 -k2,2n -r > cov/${i}.fwd_split.cov
   bedtools genomecov -ibam mapped/${i}.rev.bam -d -split | sort -k 1,1 -k2,2n -r > cov/${i}.rev_split.cov
   samtools sort -n mapped/${i}.fwd.bam | bedtools bamtobed -i stdin -bedpe | bedpe_helper | sort -k 1,1 | bedtools genomecov -i stdin -d -g $chromInfo | sort -k 1,1 -k2,2n -r > cov/${i}.fwd_fragment.cov
   samtools sort -n mapped/${i}.rev.bam | bedtools bamtobed -i stdin -bedpe | bedpe_helper | sort -k 1,1 | bedtools genomecov -i stdin -d -g $chromInfo | sort -k 1,1 -k2,2n -r > cov/${i}.rev_fragment.cov

   echo "done ... compute split intron coverage : ..."
   samtools view -H mapped/${i}.fwd.bam > temp.H
   samtools view mapped/${i}.fwd.bam | awk '$6 ~/N/' > temp.sam
   cat temp.H temp.sam | samtools view -Sb - > temp.bam
   rm temp.sam
   samtools sort -n temp.bam | bedtools bamtobed -i stdin -split | splitbed_to_splitcov | sort -k 1,1 | bedtools genomecov -i stdin -d -g $chromInfo |  sort -k 1,1 -k2,2n -r > temp0

   samtools view mapped/${i}.rev.bam | awk '$6 ~/N/' > temp.sam
   cat temp.H temp.sam | samtools view -Sb - > temp.bam
   rm temp.sam
   samtools sort -n temp.bam | bedtools bamtobed -i stdin -split | splitbed_to_splitcov | sort -k 1,1 | bedtools genomecov -i stdin -d -g $chromInfo |  sort -k 1,1 -k2,2n -r > temp0

   echo "done ... merging, computing D-ratio and merging: ..."
   paste cov/${i}.fwd_fragment.cov cov/${i}.fwd_split_introns.cov cov/${i}.fwd1.5end cov/${i}.fwd2.5end | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,$6,$9,$12}' > cov/${i}.fwd.cov
   paste cov/${i}.rev_fragment.cov cov/${i}.rev_split_introns.cov cov/${i}.rev1.5end cov/${i}.rev2.5end | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,$6,$9,$12}' > cov/${i}.rev.cov
   awk '{print(($6+0.01)/($3+0.01-$4))}'  cov/${i}.fwd.cov > cov/${i}.fwd.Dratio
   awk '{print(($6+0.01)/($3+0.01-$4))}'  cov/${i}.rev.cov > cov/${i}.rev.Dratio
   paste cov/${i}.fwd.cov cov/${i}.fwd.Dratio | awk 'BEGIN {OFS="\t"}; {print $0}' > cov/${i}.fwd_final.cov
   paste cov/${i}.rev.cov cov/${i}.rev.Dratio | awk 'BEGIN {OFS="\t"}; {print $0}' > cov/${i}.rev_final.cov

   echo "done ... cleaning up : ..."
   #rm cov/*.5end
   #rm cov/*.Dratio
   #rm cov/*d.cov
   #rm cov/*v.cov
   echo "Done. Have a good day!"
   done
fi

# ./scripts/bam_to_D-ratio.sh A782T01_accepted_hits.bam
# i=A782T01_accepted_hits.bam
#
# ./scripts/bam_to_D-ratio.sh D045T11.bam D045T12.bam D045T13.bam D045T14.bam D045T15.bam D045T16.bam D045T17.bam D045T18.bam D045T19.bam D045T20.bam D045T21.bam D045T22.bam D045T23.bam D045T24.bam D045T25.bam D045T26.bam D045T27.bam D045T28.bam D045T29.bam D045T30.bam
#./scripts/bam_to_D-ratio.sh D045T28.bam D045T29.bam D045T30.bam
#paste genome/cov/*.fwd.Dratio > genome/cov/all_single_all_rep_fwd_0.01.Dratio
#paste genome/cov/*.rev.Dratio > genome/cov/all_single_all_rep_rev_0.01.Dratio
#paste genome/cov/*.fwd.cov | awk 'BEGIN {OFS="\t"}; {print $6,$12,$18,$24,$30,$36,$42,$48,$54,$60,$66,$72,$78,$84,$90,$96,$102,$108,$114,$120,$3-$4,$9-$10,$15-$16,$21-$22,$27-$28,$33-$34,$39-$40,$45-$46,$51-$52,$57-$58,$63-$64,$69-$70,$75-$76,$81-$82,$87-$88,$93-$94,$99-$100,$105-$106,$111-$112,$117-$118}' > genome/cov/all_single_all_rep_fwd.cov
#paste genome/cov/*.rev.cov | awk 'BEGIN {OFS="\t"}; {print $6,$12,$18,$24,$30,$36,$42,$48,$54,$60,$66,$72,$78,$84,$90,$96,$102,$108,$114,$120,$3-$4,$9-$10,$15-$16,$21-$22,$27-$28,$33-$34,$39-$40,$45-$46,$51-$52,$57-$58,$63-$64,$69-$70,$75-$76,$81-$82,$87-$88,$93-$94,$99-$100,$105-$106,$111-$112,$117-$118}'  > genome/cov/all_single_all_rep_rev.cov

