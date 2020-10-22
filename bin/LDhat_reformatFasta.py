#!/usr/bin/env python
import sys

from Bio import SeqIO

inputFile = sys.argv[1]
numSequences = sys.argv[2]
numSitesInAln = sys.argv[3]
phaseType = sys.argv[4]

records = list(SeqIO.parse(inputFile, "fasta"))
# print(">entry_" + records[0].id)  # first record
# print(records[0].seq)

with open('LDhat_reformated.fa', 'w') as fileOut:
    # Add required first line for LDhat

    # From Manual: Full sequence data should be aligned and in a modified FASTA format, 
    # with the first line detailing the number of sequences/genotypes, 
    # the number of sites in the alignment and a flag (1 or 2) that 
    # details whether the data is haplotype/phased (1) or genotype/unphased (2).

    # Values used: sample size (how many genomes), genome size (size of each genome), 1 (for haplotype)
    # For the number of sites in alignment, it only requires one number so full genome size was used.
    fileOut.write(f"{numSequences} {numSitesInAln} {phaseType} \n")
    # Modify fasta header names
    for i in records:
        fileOut.write(">genome_" + i.id + '\n')
        fileOut.write(str(i.seq) + '\n')

