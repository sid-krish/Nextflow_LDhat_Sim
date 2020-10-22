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
    # Add header line for LDhat
    fileOut.write(f"{numSequences} {numSitesInAln} {phaseType} \n")
    for i in records:
        fileOut.write(">genome_" + i.id + '\n')
        fileOut.write(str(i.seq) + '\n')

