#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique

    input:
        each rho
        each theta
        each genome_size
        each sample_size
        each seed

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val("rho_${rho}_theta_${theta}_genome_size_${genome_size}_sample_size_${sample_size}_seed_${seed}")

    script:
    """
    """

}


process MSPRIME {
    publishDir "Msp_Sim_Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho),
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed),
            val(fn_modifier)

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("msp_out.fa")
             
    script:
    """
    msp_sim_fa.py ${sample_size} ${genome_size} ${rho} ${params.recom_tract_len} ${seed} ${theta}
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    params.rho_rates = [0.005] // unscaled r values. rho = 2 . p . N_e . r . tractlen
    params.theta_rates = [0.005] // unscaled u values. theta = 2 . p . N_e . u
    params.genome_sizes = [25000]
    params.sample_sizes  = [50]
    params.seeds = [1]
    
    params.recom_tract_len = 1000
    // params.effective_pop_size = 1

    rho_rates = Channel.from(params.rho_rates)
    theta_vals = Channel.from(params.theta_rates)
    genome_sizes = Channel.from(params.genome_sizes)
    sample_sizes = Channel.from(params.sample_sizes)
    seed_vals = Channel.from(params.seeds)

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need

    RATE_SELECTOR(rho_rates, theta_vals, genome_sizes, sample_sizes, seed_vals)

    MSPRIME(RATE_SELECTOR.out)

}