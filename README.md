# SplitVision
SplitVision -SV a software for extracting and analysing SV breakpoints
SplitVision accepts an SV vcf file or bedpe file, as well as a bam file. SplitVision extracts split reads(supplementary alignments) bridging the SV breakpoints. SplitVision also detects breakpoint homology and de novo insertions at the breakpoints


command:

python splitvision.py --vcf input.vcf --bam input.bam

A tab separated table will be printed to stdout
