#!/bin/bash

# This script uses hisat2 to align paired-end reads.
# It sorts the mapped reads by name and coordinates.
# instructions on https://github.com/cyaguesa/SL-quant

# REQUIREMENTS (see installation instruction if unmet):      VERSION USED 
# - bowtie2 installed and in your path                       2.4.4            
# - samtools installed and in your path                      1.14
# - bedtools installed and in your path                      2.30.0           
# - bedops installed and in your path                        2.4.40

# Carlo Yague-Sanz, 2022


# options
reference="reference/spike-in"  # basename of the bowtie2 index
sample="data/rho_seq_exemple"   # basename for raw read data (fastq)
coverage_threshold=50           # coverage threshold
Dratio_threshold=0.05           # D-ratio threshold


#Align toy exemple data to spike-in reference
bowtie2 -x $reference -1 ${sample}_R1.fq.gz -2  ${sample}_R2.fq.gz --no-discordant --no-mixed --phred33 -p 4 --no-unal | samtools view -uS - > ${sample}_aligned.bam

#Sort and index
samtools sort ${sample}_aligned.bam -o ${sample}_aligned_sorted.bam
samtools index ${sample}_aligned_sorted.bam

#Split mapped reads according to strand and orientation. Here we only process the forward strand.
samtools view -b -f 82 -F 256 ${sample}_aligned_sorted.bam > ${sample}.fwd1.bam && \
samtools index ${sample}.fwd1.bam
samtools view -b -f 130 -F 272 ${sample}_aligned_sorted.bam > ${sample}.fwd2.bam && \
samtools index ${sample}.fwd2.bam
samtools merge ${sample}.fwd.bam ${sample}.fwd1.bam ${sample}.fwd2.bam && \
samtools index ${sample}.fwd.bam && \
rm ${sample}.fwd1.bam* && \
rm ${sample}.fwd2.bam*

#Calculate the number of RT stop events per position, i.e., the number of R2 reads ending in that position.
bedtools genomecov -ibam ${sample}.fwd.bam -strand + -d -5 | \
sort -k 1,1 -k2,2n -r > ${sample}.fwd.RTstop

#Calculate the extended fragment coverage, i.e., the read coverage from the left-most coordinate to the right-most coordinate of a read pair.
samtools sort -l 0 -@ 1 -m 3G -n ${sample}.fwd.bam | \
bedtools bamtobed -i stdin -bedpe | \
awk 'BEGIN {OFS="\t"};{m=$2;M=0;for(i=2;i<=6;i++)if((i != 4)) {if(($i<m))m=$i;if(($i>M))M=$i};print $1,m,M}' | \
sort-bed --max-mem 5G - | \
bedtools genomecov -i stdin -d -g ${reference}.info | \
sort -k 1,1 -k2,2n -r > ${sample}.fwd.coverage

#Merge coverage and RTstop files across conditions and calculate D-ratio. A pseudocount of 1 is added to the numerator and denominator to avoid division by zero.
paste ${sample}.fwd.coverage ${sample}.fwd.RTstop | awk 'BEGIN {OFS="\t"};{print $1,$2,$3,$6,($6+1)/($3+1)}' > combined.Dratio

#Filter out sites with average fragment coverage < $coverage_threshold and average D-ratio < $Dratio_threshold
awk -v x=$coverage_threshold -v y=$Dratio_threshold 'BEGIN {OFS="\t"} ; ($3 > x && $5 > y && $2 > 1) {print $0}' combined.Dratio > combined.Dratio.filtered

#Find upstream nucleotide in the reference transcriptome file and keep only sites with a T/U upstream. Note that contrary to the coverage files, the intermediary bed file is zero-based for the start coordinate.
awk 'BEGIN {OFS="\t"}; {print $1, $2-2, $2-1}' combined.Dratio.filtered > filtered.bed
bedtools getfasta -tab -fi $reference.fa -bed filtered.bed -fo filtered.seq
paste combined.Dratio.filtered filtered.seq | awk 'BEGIN {OFS="\t"}; ($7 == "T") {print $0}' > combined.Dratio.filtered2
head combined.Dratio.filtered2
#The combined.Dratio.filtered2 file should have column like this:
#chrom | position | coverage | stop | Dratio | chrom | (position-2) | (position-1) | upstream_sequence
# meaning that after the unique "U" of the spike-in sequence in position 42, there is a RT-stop site in (at position 43) that passes the thresholds.

