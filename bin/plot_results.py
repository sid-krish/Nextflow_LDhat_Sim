#! /usr/bin/env python
import sys
import seaborn as sns
import pandas as pd


collected_results = sys.argv[1]

# collected_results = "../Output/Results/collected_results.csv"

df = pd.read_csv(collected_results)

df["sample_size,genome_size"] = df.apply(lambda x: str(x["sample_size"])+','+str(x["genome_size"]), axis="columns")

df = df.sort_values(by=["sample_size,genome_size"])

sns.set_theme(style="whitegrid", palette="flare")
# Flare and crest palettes are nice

g = sns.FacetGrid(df, col="rho")
g.map_dataframe(sns.barplot, x="max_rho", y="sample_size,genome_size")
g.set_axis_labels("estimated_rho", "sample_size,genome_size")

g.fig.suptitle('LDhat Pairwise : rho vs estimated_rho', y=1.05)
g.savefig("results_plot.png", dpi=500)
