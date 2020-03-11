#!/usr/local/bin/Rscript

# This script parse "split" bed files  (output of 'bedtools bamtobed -split') 
# into bed files of the regions within read coordinate that are not directly
# covered by the read sequence (that is, regions corresponding to introns).

# This coverage can then be substracted from fragment coverage to obtain physical coverage.

# It is designed to work in a pipe such as:
# bedtools bamtobed -i .bam -split | ./splitbed_to_splitcov > .bed

library(IRanges)
options(scipen = 999)
f <- file("stdin")
open(f)
current="not_an_ID"
current_ID="not_an_ID"
start=NULL
end=NULL
startA=NULL
endA=NULL
while(length(line <- readLines(f,n=1)) > 0) {
  line=strsplit(line, "\t")[[1]]
  readID=substr(line[4], 1, nchar(line[4])-2)
  if(readID==current_ID){
    if(line[4]==current){
        chrom=line[1]
        start=c(start,line[2])
        end=c(end,line[3])
    }
    if(line[4]!=current){
      if(!is.null(start)){
        startA=c(startA,end[1:(length(end)-1)])
        endA=c(endA,start[2:(length(start))])
      }
      current=line[4]
      start=line[2]
      end=line[3]
    }
  }
  if(readID!=current_ID){
    current_ID=readID
    if(!is.null(start)){
      startA=c(startA,end[1:(length(end)-1)])
      endA=c(endA,start[2:(length(start))])
    }
    current=line[4]
    start=line[2]
    end=line[3]
    if (!is.null(startA)){
      answer=reduce(IRanges(as.numeric(startA), as.numeric(endA)))
      startA=NULL;endA=NULL
      for (j in c(1:(length(answer)))){
        write(paste(chrom,start(answer)[j],end(answer)[j], sep="\t"),file="")
      }
    }
  }
}
if(!is.null(start)){
  startA=c(startA,end[1:(length(end)-1)])
  endA=c(endA,start[2:(length(start))])
}
if (!is.null(startA)){
  answer=reduce(IRanges(as.numeric(startA), as.numeric(endA)))
  startA=NULL;endA=NULL
  for (j in c(1:(length(answer)))){
    write(paste(chrom,start(answer)[j],end(answer)[j], sep="\t"),file="")
  }
}
