#!/bin/bash

# This script quality-trim the D-seq paired-end raw reads
# and maps the remaining reads on the spike-in sequence
# It requires bowtie2 and samtools installed and in the $PATH.
# Path to trimmomatic should be manually set.

################  PARAMETERS
set -e #stops script at the first error
cd ~/Desktop/these_Carlo/R_terminal/D-seq/2018_D_seq_pombe/raw
javapath="/Users/carlo/Documents/binaries/Trimmomatic-0.36/trimmomatic-0.36.jar" #path to trimmomatic
mkdir -p ../filtered #create trimming output directory
mkdir -p ../spike-in #create spike-in output directory
spike_index="../spike-in/spike_in_T7_rRNA" #path to bowtie 2 spike index basename
###############################


for i in $(ls *.fastq.gz | rev | cut -c 13- | rev | uniq)
  do
  echo ""; echo "fastq file : ${i}"
  echo "... trimming..."
  java -jar $javapath PE -threads 4 -phred33 ${i}.R1.fastq.gz ${i}.R2.fastq.gz ../filtered/${i}.R1.filtered.paired.gz ../filtered/${i}.R1.filtered.unpaired.gz ../filtered/${i}.R2.filtered.paired.gz ../filtered/${i}.R2.filtered.unpaired.gz ILLUMINACLIP:/Users/carlo/Documents/binaries/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:20
  echo "... done... aligning..."
  bowtie2 -x $spike_index -1 ../filtered/${i}.R1.filtered.paired.gz -2 ../filtered/${i}.R2.filtered.paired.gz --no-discordant --no-mixed --phred33 -p 4 --un-conc ../filtered/${i}.filtered_final.fastq --no-unal | samtools view -uS - > ../spike-in/${i}_spike_in_reads.bam
  echo "... done... gzipping .fastq"
  gzip ../filtered/${i}.filtered_final.*.fastq
  echo "... done... sorting and indexing .bam"
  samtools sort ../spike-in/${i}_spike_in_reads.bam -o ../spike-in/${i}_spike_in_reads_sorted.bam
  samtools index ../spike-in/${i}_spike_in_reads_sorted.bam
  rm ../spike-in/${i}_spike_in_reads.bam
  echo "... done."
  done
