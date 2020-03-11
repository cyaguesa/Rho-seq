#!/usr/local/bin/Rscript

# This script parse .bedpe files (output of 'bedtools bamtobed -bedpe')
# in bed files with coordinates corresponding to the leftmost and rightmost
# coordinates of the read pair.

# It is designed to work in a pipe such as:
# bedtools bamtobed -i .bam -bedpe | ./bedpe_helper.R > .bed

# improvement : It should be possible, and faster, to do the same thing using awk.
# awk '{m=$2;M=0;for(i=2;i<=6;i++)if((i != 4)){if(($i<m))m=$i;if(($i>M))M=$i};print $1,m,M}'

f <- file("stdin")
options(scipen = 999)
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  line=strsplit(line, "\t")[[1]]
  write(paste(line[1],(min(as.numeric(line[c(2,3,5,6)]))),(max(as.numeric(line[c(2,3,5,6)]))), sep="\t"),file="")
}
