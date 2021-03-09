#! /usr/bin/env python
import sys

run_result = sys.argv[1]
rho_rate =  sys.argv[2]
sample_size =  sys.argv[3]
genome_size =  sys.argv[4]

# run_result = "pairwise_outfile.txt"
# rho_rate = 15
# sample_size = 20
# genome_size = 30000

with open(run_result, 'r') as f:
    lines = f.readlines()

max_rho, max_lk = lines[4].split(':')

max_rho = max_rho.strip()
max_lk = max_lk.strip()

max_rho = max_rho[26:]
max_lk = max_lk.strip("Lk =")

with open("processed_results.csv", 'w') as file:  # open in write mode (create new file)
    file.write(f"rho,sample_size,genome_size,max_rho,max_lk\n")
    file.write(f"{rho_rate},{sample_size},{genome_size},{max_rho},{max_lk}\n")

