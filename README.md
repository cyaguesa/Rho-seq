# Rho-seq

*Transcription-wide distribution of dihydrouridine (D) into mRNAs reveals its requirement for meiotic chromosome segregation.*

_Raw data on GEO: [GSE145686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145686)_

## Summary

Here we have developed Rho-seq, an integrated pipeline detecting a range of modifications through differential modification-dependent Rhodamine labeling. Using Rho-seq, we report that the reduction of uridine to dihydrouridine by the Dus reductase family targets tRNAs in E. coli but expands to mRNAs is yeast. The modified mRNAs are enriched for cytoskeleton related encoded protein. We show that the a-tubulin encoding mRNA nda2 undergoes dihydrouridination, which affects its protein expression level. The absence of the modification onto the nda2 mRNA strongly impacts meiosis by inducing a metaphase delay or by completely preventing the formation of spindles during meiosis I and meiosis II, eventually resulting in low gamete viability. Collectively these data show that the codon specific reduction of uridine within mRNA is required for proper meiotic chromosome segregation and gamete viability.


## Analysis

-  _1_mapping_spike.sh:_ Reads were quality- and adaptor-trimmed using trimmomatic 0.36 with options `PE ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 HEADCROP:1 SLIDINGWINDOW:4:15 MINLEN:20`. The filtered reads were then mapped on the spike-in sequence with bowtie2 (options `--no-discordant --no-mixed --phred33 -p 4 --un-conc`).

-  _2_mapping_genome.sh:_ The reads that did not map on the spike-in sequence were remapped on the organism genome using also bowtie2 for the E.coli libraries or tophat2 for the S. pombe libraries (options `--mate-std-dev 150 -N 2 --read-edit-dist 2 --min-intron-length 5 --max-intron-length 3000 --no-discordant --no-coverage-search -p 4 --library-type fr-firststrand`).

- _3_bam_to_D-ratio.sh:_ Properly paired mapped reads were split by strand. Then bedtools (options `genomecov -ibam -d -5`) was used to calculate the number of RT termination events for each strand-specific position as the number of 5'end of R2 reads that map at this position. In parallel, the paired reads were computationally extended/merged into fragments that covered the genome from the left end of the left-most read to the right end of the right-most reads with the notable exeption that the fragments were not extended in identified introns. Those spliced fragments where then used to compute a fragment coverage reflecting the number of time a reverse transcriptase passed through a specific position. Then the ratio between the number of RT termination events and fragment coverage at each position (the D-ratio) was calculated.

- _4_D-ratios_independant_filtering.sh_: After independent pre-filtering of the positions unlikely to have significant differences in D-ratio accross conditions (see manuscrit), the number of RT termination events out of the number of RT readthrough were modelised in generalized linear model of the binomial family (function glm in R) according to the effect of the treatment (R+ vs R-), the effect of the strain (wild-type vs delta4dus), and their interaction. The positions were the interaction between the wild-type and the delta4dus lead to a statistically significant (FDR<0.1) increase of the proportion of RT termination events were retained as putative dihydrouridilation sites.
