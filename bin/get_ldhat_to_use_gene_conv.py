#!/usr/bin/env python
import sys

locs_file = sys.argv[1]

# locs_file = "../locs.txt"

with open(locs_file, 'r') as fileIn:
    locs_file_list = fileIn.readlines()

# change from crossover/linear (L) to gene conversion (C) mode
old_first_line = list(locs_file_list[0])
old_first_line[-2] = 'C'
new_first_line = ''.join(old_first_line)
locs_file_list[0] = new_first_line

with open("locs_C.txt", 'w') as fileOut:
    fileOut.writelines(locs_file_list)
