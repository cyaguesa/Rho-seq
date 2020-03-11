#!/bin/bash

# This script maps the remaining reads on the genome
# It requires tophat2 and samtools installed and in the $PATH.

################  PARAMETERS
set -e #stops script at the first error
cd /Users/carlo/Desktop/these_Carlo/R_terminal/D-seq/2018_D_seq_pombe_single_mutants
genome_index="~/Desktop/these_Carlo/data_raw/annotation/pombe/bowtie_index_pombe" #path to tophat 2 genome index basename
#############################

for i in $(ls filtered/*final*.fastq.gz | rev | cut -c 27- | rev | uniq)
  do
  echo ""; echo "fastq file : ${i}"
  echo "... done... mapping..."
#--max-multihits 20 = default setting
  tophat2 --rg-sample ${i} --rg-id ${i} -r 200 --mate-std-dev 150 -N 2 --read-edit-dist 2 --min-intron-length 5 --max-intron-length 3000 --no-discordant --no-coverage-search -p 4 --library-type fr-firststrand --output-dir ${i} $genome_index ${i}.filtered_final.1.fastq.gz ${i}.filtered_final.2.fastq.gz
  samtools index ${i}/accepted_hits.bam
  samtools flagstat ${i}/accepted_hits.bam > ${i}/accepted_hits.bam.flagstat
  echo "... done. Good job"
  done

