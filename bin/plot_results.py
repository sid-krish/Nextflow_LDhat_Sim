#! /usr/bin/env python
import sys
import seaborn as sns
import pandas as pd

collected_results = sys.argv[1]

# collected_results = "../Output/Results/collected_results.csv"

df = pd.read_csv(collected_results)

df["sample_size,genome_size"] = df.apply(lambda x: str(x["sample_size"])+','+str(x["genome_size"]), axis="columns")

df = df.sort_values(by=["sample_size,genome_size"])

plot = sns.barplot(data=df, x="rho", y="max_rho", hue="sample_size,genome_size")

plot.set_title("Pairwise : rho(simulated) vs max_rho(estimated)")

fig = plot.get_figure()

fig.savefig("results_plot.png")
