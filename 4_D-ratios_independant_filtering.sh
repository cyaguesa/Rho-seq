#!/bin/bash

# Find sites where that could support RT stops using independant filtering
# This code requires bedtools installed and in your PATH.
# Get the sequence of those sites and the gene they overlap

# genome annotation are annot/Schizosaccharomyces_pombe.ASM294v2.26_no_CDS.gff | annot/Schizosaccharomyces_pombe.ASM294v2.26.gff3 | annot/Schizosaccharomyces_pombe.fasta 

#####################
# parameters
#####################
cov_expression_cutoff=100 # the site should be sufficiently transcribed
averagedDratio=0.018 # the site should support sufficient stop
sd_Dratio=0.02 # there should be some variability between 
rep1="/Users/carlo/Desktop/these_Carlo/R_terminal/D-seq/fastq2/cov" # Path to replicate 1 cov
rep2="/Users/carlo/Desktop/these_Carlo/R_terminal/D-seq/duplicate_pombe/genome/cov" # Path to replicate 2 cov
annot="/Users/carlo/Desktop/these_Carlo/data_raw/annotation/pombe" # path to annot
#####################

## getting data and prefiltering

paste $rep1/A560T1_accepted_hits.bam.fwd_final.cov $rep1/A560T2_accepted_hits.bam.fwd_final.cov $rep1/A560T3_accepted_hits.bam.fwd_final.cov $rep1/A560T4_accepted_hits.bam.fwd_final.cov $rep2/A782T01_accepted_hits.bam.fwd_final.cov $rep2/A782T02_accepted_hits.bam.fwd_final.cov $rep2/A782T03_accepted_hits.bam.fwd_final.cov $rep2/A782T04_accepted_hits.bam.fwd_final.cov | awk 'BEGIN {OFS="\t"};{ print $0, ($3+$10+$17+$24+$31+$38+$45+$52)/8, ($7+$14+$21+$28+$35+$42+$49+$56)/8}' > cov/full_data_fwd.cov
awk '($57 > 50) && ($58 > 0.018) { print $0 "\t+"}' cov/full_data_fwd.cov  > cov/full_data_I_filtered.cov
cut -f 1-3,6,7,10,13,17,20,24,27,31,34,38,41,45,48,52,55-60,14,21,28,35,42,49,56 cov/full_data_I_filtered.cov > $rep2/full_data_I_filtered.cov2
cut -f 1-3,6,7,10,13,17,20,24,27,31,34,38,41,45,48,52,55-60,14,21,28,35,42,49,56 cov/full_data_fwd.cov > cov/full_data_fwd_clean.cov

paste $rep1/A560T1_accepted_hits.bam.rev_final.cov $rep1/A560T2_accepted_hits.bam.rev_final.cov $rep1/A560T3_accepted_hits.bam.rev_final.cov $rep1/A560T4_accepted_hits.bam.rev_final.cov $rep2/A782T01_accepted_hits.bam.rev_final.cov $rep2/A782T02_accepted_hits.bam.rev_final.cov $rep2/A782T03_accepted_hits.bam.rev_final.cov $rep2/A782T04_accepted_hits.bam.rev_final.cov | awk 'BEGIN {OFS="\t"};{ print $0, ($3+$10+$17+$24+$31+$38+$45+$52)/8, ($7+$14+$21+$28+$35+$42+$49+$56)/8}' > $rep2/full_data_rev.cov
awk '($57 > 50) && ($58 > 0.018) { print $0 "\t-"}' $rep2/full_data_rev.cov  > $rep2/full_data_I_filtered_rev.cov  
cut -f 1-3,6,7,10,13,17,20,24,27,31,34,38,41,45,48,52,55-60,14,21,28,35,42,49,56 $rep2/full_data_I_filtered_rev.cov > $rep2/full_data_I_filtered_rev.cov2
cut -f 1-3,6,7,10,13,17,20,24,27,31,34,38,41,45,48,52,55-60,14,21,28,35,42,49,56 cov/full_data_rev.cov > cov/full_data_rev_clean.cov

## [R-process here]

## Getting sequences of prefiltered data 

echo "done... Converting into bed...";  
awk 'BEGIN {OFS="\t"};{ print $1,($2-2),($2-1),".",".",$29}' $rep2/full_data_I_filtered_fwd.cov3 > $rep2/full_data_fwd.cov3.bed
awk 'BEGIN {OFS="\t"};{ print $1,($2),($2+1),".",".",$29}' $rep2/full_data_I_filtered_rev.cov3 > $rep2/full_data_rev.cov3.bed
cat $rep2/full_data_fwd.cov3.bed $rep2/full_data_rev.cov3.bed > $rep2/full_data.cov3.bed
bedtools getfasta -s -tab -fi $annot/Schizosaccharomyces_pombe.fasta -bed $rep2/full_data.cov3.bed -fo $rep2/cov3.seq

## [R-process here]

## analysis of sites

   echo "done... Converting into bed...";  
   awk 'BEGIN {OFS="\t"};{ print $1,($2-3),$2,".",".",$3}' $rep2/full_data_I_filtered_fwd.cov4 > $rep2/full_data_fwd.Dsites.bed
   awk 'BEGIN {OFS="\t"};{ print $1,($2-1),($2+2),".",".",$3}' $rep2/full_data_I_filtered_rev.cov4 > $rep2/full_data_rev.Dsites.bed
   cat $rep2/full_data_fwd.Dsites.bed $rep2/full_data_rev.Dsites.bed > $rep2/full_data_ALL.Dsites.bed

   echo "done... Match genes and sequence...";  
   bedtools intersect -s -wa -wb -a $annot/Schizosaccharomyces_pombe.ASM294v2.26_no_CDS.gff -b $rep2/full_data_ALL.Dsites.bed| egrep "Name=gene|Name=tRNA_gene|Name=ncRNA_gene|Name=snoRNA_gene" > $rep2/full_data_ALL.Dsites.txt
   bedtools intersect -s -wa -wb -a  $annot/Schizosaccharomyces_pombe.ASM294v2.26.gff3 -b $rep2/full_data_ALL.Dsites.bed > $rep2/full_data_features.Dsites.txt
   bedtools getfasta -s -tab -fi $annot/Schizosaccharomyces_pombe.fasta -bed $rep2/full_data_ALL.Dsites.bed -fo $rep2/full_data.seq
   echo "done."
