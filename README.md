# SplitVision
SplitVision -SV a software for extracting and analysing SV breakpoints
SplitVision accepts an SV vcf file or bedpe file, as well as a bam file. SplitVision extracts split reads(supplementary alignments) bridging the SV breakpoints. SplitVision also detects breakpoint homology and insertions at the breakpoints.

# Run

First activate the conda environment

    source activate SplitVision

To analyse the breakpoints in a vcf:

    python splitvision.py --analyse --vcf input.vcf --bam input.bam --fa ref.fasta --db repeats.db --working_dir output

To analyse the breakpoints in a bedpe file:

    python splitvision.py --analyse --vcf input.tab --bam input.bam --fa ref.fasta --db repeats.db --working_dir output

To generate a database of repeats:

    python splitvision.py --db --tab repeatmask.tab --prefix output

The repeatmask.tab file may be downloaded via USCS tablebrowser. the output file is a sqlite database containing all the repeats of the input file.

The software will create an excel file containing the positions and features of each juction. These properties are found based on split reads.

# Install

Install the SplitVision conda environemnt:

    ./create_conda_env.sh


The conda environment contains the following tools:

    ABYSS
    ClustalW2
    cd-hit
    xlwt(python module)
    freebayes
    vt

SplitVision requires these tools aswell:

    bwa
    samtools

These are not installed thtough the conda install script. If you do not have samtools and bwa, these could be installed using conda aswell.

The reference fasta needs to be indexed using bwa.

# Algorithm
Splitvision extracts split reads (supplementary + primary alignment) from the breakpoint junctions found in the vcf or bed file. Splitvision will search for these reads within the distance given by the "padding" parameter. These split reads are collected into a fasta file, and clustered using cd-hit. If multiple clusters are formed, the biggest cluster(i.e the cluster containing most reads) will be analysed. The reads of that cluster is analysed using ClustalW. Using the output of clustalW, a consensus sequence of the breakpoint is generated. The consensus sequence is aligned to the reference using BWA mem.

If no split reads are found, ABYSS will atempt to produce a contig spanning the  breakpoint junction. This is done by extracting and asssebling all reads within the paddding distance of the breakpoints.

This sequence is analysed for microhomology and breakpoint insertions, and the sequence itself is printed to the excel file.
Additionally, freebayes is used to call SNVs located within the user set snp_distance. These snvs are decomposed and normalised using vt.
