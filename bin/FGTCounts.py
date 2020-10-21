#!/usr/bin/env python
import pandas as pd
import sys

# Add recom rate, mutation rate, seed value, recom rate and count per run
# All of this will be collected to create a collected csv at the end

FGT_results = sys.argv[1]

seedVal = int(sys.argv[2])

mutation_rate = float(sys.argv[3])

recombination_rate = float(sys.argv[4])

df = pd.read_csv(FGT_results)
numRows = df.shape[0]

with open ("summaryCounts.csv", 'w') as f:
    f.write(f"seed_val,mut_rate,recom_rate,count\n")
    f.write(f"{seedVal},{mutation_rate},{recombination_rate},{numRows}\n")