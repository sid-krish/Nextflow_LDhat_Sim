#!/usr/bin/env python

import os
import numpy as np

"""
TODO: Parametric sweep + Replicates

####
Sweep 1: Recombination rate simulations
####

Rho: 0.1, 0.5, 2.5, 25, 50
Theta: 0.01 (fixed like PIIM)
Genome_size: 100,000; 250,000; 500,000; 750,000; 1,000,000
mean_cov: .5, 1, 2.5, 5, 10
Samples: 10 (currently fixed)

####
3 Replicates with different seed values
####

Seed: 123, 456, 789

"""


def parametric_sweep(rho, theta, genome_size, sample_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, genome_size, sample_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Generate simulated datasets
    for rho, theta, genome_size, sample_size, seed in sweep_1_combinations:
        # subprocess.run(["nextflow", "run", "sim_LDhat.nf",
        #                 "--rho_rates", f"{rho}",
        #                 "--mutation_rate", f"{theta}",
        #                 "--genome_sizes", f"{genome_size}",
        #                 "--sample_sizes", f"{sample_size}",
        #                 "--seed", f"{seed}"])

        cmd = f"nextflow run sim_LDhat.nf --rho_rates {rho} --mutation_rate {theta} --genome_sizes {int(genome_size)} --sample_sizes {int(sample_size)} --seed {int(seed)}"
        os.system(cmd)

    return None


if __name__ == '__main__':

    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.01, 0.025, 0.05, 0.075, 0.1]
    theta_sweep_1 = [0.01]
    genome_size_sweep_1 = [10000, 25000, 50000, 75000, 100000]
    sample_size_sweep_1 = [10]
    seed_sweep_1 = [123, 456, 789]

    parametric_sweep(rho_sweep_1, theta_sweep_1, genome_size_sweep_1, sample_size_sweep_1, seed_sweep_1)
