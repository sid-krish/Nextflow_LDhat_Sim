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


process MS {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        val rho_rate
        val sample_size
        val seed
        val genome_size
        val path_fn_modifier

    output:
        path "trees.txt", emit: trees_txt
  
    script:
//     To carry out simulations with gene conversion but no crossing-over, one uses 
//     the -r option with  equal to zero, and the -c option. In this case value following
//     the -c option is the value of 4N0g, rather than the ratio g=r.
    """
    ms ${sample_size} 1 -T -seeds ${seed}  -t ${params.mutation_rate} -r 0 ${genome_size} -c ${rho_rate} ${params.recom_tract_len} > trees.txt
    """
}


process FAST_SIM_BAC {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

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
            path("trees.txt")
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${seed} -T -t ${theta} -r ${rho} ${params.recom_tract_len} > trees.txt
    """
}


process CLEAN_TREES {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}
    
    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("trees.txt")


    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("cleanTrees.txt")

    script:
    """
    clean_trees.py trees.txt
    """
}


process SEQ_GEN {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("cleanTrees.txt")

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("seqgenOut.fa")

    script:
    // 1 partiion per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process LDHAT_REFORMAT_FASTA{
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path(seqgenOut)

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("LDhat_reformated.fa")

    script:
    """
    LDhat_reformat_fasta.py ${seqgenOut} "${sample_size}" "${genome_size}" 1
    """
}


process LOOKUP_TABLE_LDPOP {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("LDhat_reformated.fa")

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("LDhat_reformated.fa"),
            path("lookupTable.txt")

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores 4 -n ${sample_size} -th ${theta} -rh ${params.ldpop_rho_range} --approx > lookupTable.txt
    """
}


process LDHAT_CONVERT{
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("LDhat_reformated.fa"),
            path("lookupTable.txt")

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("lookupTable.txt"),
            path("locs.txt"),
            path("sites.txt")

    script:
        // -2only: Specifies that only polymorphic sites with exactly two alleles
        // will be analysed and outputted Although only those sites with two alleles
        // are analysed in pairwise and interval, outputting all segregating sites
        // may be of interest and can be used to estimate a finite-sites estimate of
        // Watterson’s theta per site within pairwise
        """
        convert -seq LDhat_reformated.fa -2only
        """

}


process SWITCH_TO_GENE_CONVERSION_MODE{
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("lookupTable.txt"),
            path("locs.txt"),
            path("sites.txt")

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("lookupTable.txt"),
            path("sites.txt"),
            path("locs_C.txt")

    script:
        """
        get_ldhat_to_use_gene_conv.py locs.txt
        """

}


process LDHAT_INTERVAL{
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        path lookup_table_file
        path freqs_forLDhatInterval
        path locs_forLDhatInterval
        path sites_forLDhatInterval
        val path_fn_modifier

    output:
        path "bounds.txt", emit: bounds_forLDhatStat
        path "rates.txt", emit: rates_forLDhatStat
        path "new_lk.txt", emit: new_lk_txt
        path "type_table.txt", emit: type_table_txt
        path "intervalOut.txt", emit: interval_out_txt

    script:
        // pre-generated lookup tables available from the ldhat github page was used. One that matches number of samples and mutation rate was selected.
        // the arguments -its, -samp, -bpen use recommened values given in the manual for the interval program.

        // The information printed on screen was useful so decided to save that also.
        """
        interval -seq sites.txt -loc locs.txt -lk lookupTable.txt -its 1000000 -samp 2000 -bpen 5 > intervalOut.txt
        """

}


process LDHAT_INTERVAL_STAT{
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        path rates_forLDhatStat
        val path_fn_modifier

    output:
        path "res.txt" 
        path "statOut.txt"

    script:
        // The information printed on screen was useful so decided to save that also.
        """
        stat -input rates.txt > statOut.txt
        """

}


process LDHAT_PAIRWISE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("lookupTable.txt"),
            path("sites.txt"),
            path("locs_C.txt")

    output:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt")

    script:
        // uses pexpect to handle unavoidale prompts
        """
        run_pairwise_with_pexpect.py ${params.recom_tract_len} sites.txt locs_C.txt lookupTable.txt > pairwise_stdOut.txt
        """

}


process PAIRWISE_PROCESS_OUTPUT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt")

    output:
        path "processed_results.csv", emit: processed_results_csv

    script:
        """
        pairwise_process_output.py pairwise_outfile.txt ${rho} ${sample_size} ${genome_size}
        """

}

process PLOT_RESULTS{
    publishDir "Output/Results", mode: "copy"

    input:
        path collectedFile

    output:
        path "rho_comparison.png"
        path "max_lk_comparison.png"

    script:
        """
        plot_results.py collected_results.csv
        """

}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    params.seed = 123
    params.mutation_rate = 0.01
    params.recom_tract_len = 500
    params.ldpop_rho_range = "101,100"
    params.effective_pop_size = 1
    params.rho_rates = 0.001
    params.sample_sizes  = 20
    params.genome_sizes = 25000

    rho_rates = Channel.from(params.rho_rates)
    theta_vals = Channel.from(params.mutation_rate)
    genome_sizes = Channel.from(params.genome_sizes)
    sample_sizes = Channel.from(params.sample_sizes)
    seed_vals = Channel.from(params.seed)

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need

    RATE_SELECTOR(rho_rates, theta_vals, genome_sizes, sample_sizes, seed_vals)

    // MS(RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, params.seed, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(MS.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    FAST_SIM_BAC(RATE_SELECTOR.out)

    CLEAN_TREES(FAST_SIM_BAC.out)

    SEQ_GEN(CLEAN_TREES.out)

    LDHAT_REFORMAT_FASTA(SEQ_GEN.out)

    LOOKUP_TABLE_LDPOP(LDHAT_REFORMAT_FASTA.out)

    LDHAT_CONVERT(LOOKUP_TABLE_LDPOP.out)

    SWITCH_TO_GENE_CONVERSION_MODE(LDHAT_CONVERT.out)

    // LDHAT_INTERVAL(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, LDHAT_CONVERT.out.freqs_txt, LDHAT_CONVERT.out.locs_txt, LDHAT_CONVERT.out.sites_txt, RATE_SELECTOR.out.path_fn_modifier)

    // LDHAT_INTERVAL_STAT(LDHAT_INTERVAL.out.rates_forLDhatStat, RATE_SELECTOR.out.path_fn_modifier)

    LDHAT_PAIRWISE(SWITCH_TO_GENE_CONVERSION_MODE.out)

    PAIRWISE_PROCESS_OUTPUT(LDHAT_PAIRWISE.out)

    // collectedFile = PAIRWISE_PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    // PLOT_RESULTS(collectedFile)
}