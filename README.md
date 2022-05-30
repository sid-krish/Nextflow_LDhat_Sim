<h1 align="center">Nextflow_LDhat_sim - Bacterial Sequence Simulation Pipeline</h1>
  <p align="center">
    Bacterial sequence simulator with recomnination and mutation

- [About](#about)
  - [Built With](#built-with)
- [Getting Started](#getting-started)
  - [Set up using conda](#set-up-using-conda)
  - [Set up using docker](#set-up-using-docker)
- [Quick Start and Output](#quick-start-and-output)
- [Pipeline Options](#pipeline-options)
  - [sim_LDhat.nf](#sim_ldhatnf)
- [Issues and Contributing](#issues-and-contributing)
- [Contact](#contact)

<!-- ABOUT -->
## About
Bacterial sequence simulator based on msprime. Designed to simulate sequences with mutation and gene conversion type recombination. Originally designed to generate simualated sequences to help test LDHat https://github.com/auton1/LDhat, the nextflow pipeline for LDhat pairwise in gene conversion mode https://github.com/sid-krish/Nextflow_LDhat and the reimplementation of LDhat pairwise in gene conversion mode https://github.com/sid-krish/Rhometa_Full_Genome.

### Built With
* [Nextflow](https://www.nextflow.io/)
* [Msprime](https://tskit.dev/msprime/docs/stable/intro.html)


<!-- GETTING STARTED -->
## Getting Started

Nextflow_LDhat_sim is designed to be run on linux and requires nextflow to be installed. 
Dependencies are resolved either via conda or docker images. Support for HPC, docker, singularity, AWS and many other systems are provided via nextflow.

While it is possible to resolve the dependencies using conda for running on macOS, its recommended that this option be used on linux systems for which it has been extensively test.
If running on macOS it recommended that docker be used with the provided image, in which case it is similar to running in a linux environment.

It is also possible to install and run the program on Windows via [wsl](https://docs.microsoft.com/en-us/windows/wsl/install).

### Set up using conda
Instructions for installing nextflow and dependencies via conda
1. Clone the repo
   ```sh
   git clone https://github.com/sid-krish/Nextflow_LDhat_Sim.git
   ```
2. Install the conda package manager: [Miniconda download](https://conda.io/en/latest/miniconda.html)
3. Install nextflow
   ```sh
   conda install -c bioconda nextflow
   ```
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
   Disable the use of docker by setting the docker option to false. Disabling the use of container engines will cause conda packages to be used by default:
   ```sh
   docker {
       enabled = false
   }
   ```
5. The pipeline is now ready to run, and all dependencies will be automatically resolved with conda.

### Set up using docker
Instructions for installing nextflow and using the provided docker image for dependencies
1. Clone the repo
   ```sh
    git clone https://github.com/sid-krish/Nextflow_LDhat_Sim.git
   ```
2. Install nextflow [Nextflow install](https://www.nextflow.io/index.html#GetStarted)
3. Install docker desktop [Docker install](https://docs.docker.com/desktop/linux/).
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
5. In the sim_gen.nf file comment the lines related to conda, for instance:
   ```
   // conda 'conda-forge::msprime=1.1.1 conda-forge::gsl'
   ```
6. Ensure docker is running.
7. The pipeline is now ready to run, all the required dependencies are present in the docker image, that the pipeline is preconfigured to use.


<!-- QUICK START AND OUTPUT -->
## Quick Start and Output
The pipeline is preconfigured to generate 2 small multi-sequence fasta files for testing purposes. This following section shows how to generate them.

For this example the command to run is:
```sh
nextflow run sim_LDhat.nf
```

By default, running the command will output the files to 'Msp_Sim_Output'.

The simulation parameters are reflected in the filenames. Note the recomibation and mutation rates are the unscaled values. This will be discussed further in the pipeline options section.


<!-- PIPELINE OPTIONS -->
## Pipeline Options
### sim_LDhat.nf
In the workflow section of the sim_LDhat.nf script, the following parameters can be adjusted, these are the available options for the pipeline.

```
params.recom_rates = [0.005, 0.01] // Unscaled recombination rate. This value is scaled as such: 2 * ploidy (1) * effective_population_size (1) * unscaled_recombination_rate * tract_length
params.mutation_rates = [0.005] // Unscaled mutation rate. This value is scaled as such: 2 * ploidy (1) * effective_population_size (1) * unscaled_mutation_rate
params.genome_sizes = [50000]
params.sample_sizes  = [20]
params.seeds = [123]

params.tract_len = 1000
```

It is important to note that the recom_rates and mutation_rates inputs are unscaled values, the final simulated values will be the scaled values following the formulation shown next to the parameters as above.

Additionaly, in the sim_LDhat.nf script where the input to options are in the form of a list, multiple values can be provided, for instance:

```
params.recom_rates = [0.005, 0.01]
```

<!-- ISSUES AND CONTRIBUTING -->
## Issues and Contributing
If you have any issues please open an issue with the details and steps for reproducing the issue. If you have any questions please open a issue with the tag "question" or alternatively email one of the authors from the contact section.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".


<!-- CONTACT -->
## Contact
Sid Krishnan - sidaswar.krishnan-1@student.uts.edu.au, sid.kr15n@gmail.com \
Aaron Darling - aaron.darling@uts.edu.au \
Matt DeMaere - matthew.demaere@uts.edu.au