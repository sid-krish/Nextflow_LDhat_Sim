#!/usr/bin/env python
import sys
import pysam


inputFile = sys.argv[1]
numSequences = int(sys.argv[2])
numSitesInAln = int(sys.argv[3])
phaseType = int(sys.argv[4])

# inputFile = "../Output/s_3_m_0.1_r_0.1/s_3_m_0.1_r_0.1_Aligned.sam"
# numSequences = 100
# numSitesInAln = 250000
# phaseType = 1

samfile = pysam.AlignmentFile(inputFile)

with open('LDhat_formated.fa', 'w') as fileOut:
    # Add required first line for LDhat

    # From Manual: Full sequence data should be aligned and in a modified FASTA format,
    # with the first line detailing the number of sequences/genotypes,
    # the number of sites in the alignment and a flag (1 or 2) that
    # details whether the data is haplotype/phased (1) or genotype/unphased (2).

    fileOut.write(f"{numSequences} {numSitesInAln} {phaseType}\n")
    for read in samfile.fetch():
        fileOut.write(f">{read.query_name}\n")
        fileOut.write(f"{read.query_sequence}\n")

