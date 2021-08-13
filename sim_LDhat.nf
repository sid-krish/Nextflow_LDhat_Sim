#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    maxForks 1 // Run sequentially

    input:
        each rho_rate
        each sample_size
        each genome_size

    output:
        val "${rho_rate}", emit: rho_rate
        val "${sample_size}", emit: sample_size
        val "${genome_size}", emit: genome_size

        val "rho_${rho_rate}_sam_${sample_size}_gen_${genome_size}/rho_${rho_rate}_sam_${sample_size}_gen_${genome_size}", emit: path_fn_modifier

    script:
    """
    """

} 


process MS {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 

    input:
        val rho_rate
        val sample_size
        val seed
        val genome_size
        val path_fn_modifier

    output:
        path "trees.txt", emit: trees_txt
  
    script:
    """
    ms ${sample_size} 1 -T -seeds ${seed}  -t ${params.mutation_rate} -r 0 ${genome_size} -c ${rho_rate} ${params.recom_tract_len} > trees.txt
    """
}


process FAST_SIM_BAC {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        // path "rho_calc.txt", emit: rho_rho_calc_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${params.seed} -T -t ${params.mutation_rate} -r ${rho_rate} ${params.recom_tract_len} > trees.txt
    #calc_rho.py ${params.effective_pop_size} ${rho_rate} ${genome_size}
    """
}


process CLEAN_TREES {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}
    maxForks 1

    input:
        path trees
        val path_fn_modifier


    output:
        path "cleanTrees.txt", emit: cleanTrees_txt

    script:
    """
    clean_trees.py trees.txt
    """
}


process SEQ_GEN {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path cleanTrees
        val genome_size
        val path_fn_modifier

    output:
        path "seqgenOut.fa", emit: seqgenout_fa

    script:
    // 1 partiion per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${params.seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process LDHAT_REFORMAT_FASTA{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path seqgenOut
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "LDhat_reformated.fa", emit: ldhat_reformated_fa

    script:
    """
    LDhat_reformat_fasta.py ${seqgenOut} "${sample_size}" "${genome_size}" 1
    """
}


process LOOKUP_TABLE_LDPOP {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val sample_size
        val path_fn_modifier

    output:
        path "lookupTable.txt", emit: lookupTable_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores 4 -n ${sample_size} -th ${params.mutation_rate} -rh ${params.ldpop_rho_range} --approx > lookupTable.txt
    """
}


process LDHAT_CONVERT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path fasta_forLDhatConvert
        val path_fn_modifier

    output:
        path "freqs.txt", emit: freqs_txt
        path "locs.txt", emit: locs_txt
        path "sites.txt", emit: sites_txt
        path "convertOut.txt", emit: convert_out_txt

    script:
        // The information printed on screen was useful so decided to save that also.
        // -2only, only output sites with exactly two alleles
        """
        convert -seq LDhat_reformated.fa -2only > convertOut.txt
        """

}


process SWITCH_TO_GENE_CONVERSION_MODE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path locs_txt
        val path_fn_modifier

    output:
        path "locs_C.txt", emit: locs_C_txt

    script:
        """
        get_ldhat_to_use_gene_conv.py locs.txt
        """

}


process LDHAT_INTERVAL{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path lookup_table_file
        path locs_C_txt
        path sites_txt
        val path_fn_modifier

    output:
        path "pairwise_freqs.txt", emit: pairwise_freqs_txt
        path "pairwise_outfile.txt", emit: pairwise_outfile_txt
        path "pairwise_stdOut.txt", emit: pairwise_stdOut_txt

    script:
        // uses pexpect to handle unavoidale prompts
        """
        run_pairwise_with_pexpect.py ${params.recom_tract_len} sites.txt locs_C.txt lookupTable.txt > pairwise_stdOut.txt
        """

}


process PAIRWISE_PROCESS_OUTPUT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path pairwise_outfile_txt
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "processed_results.csv", emit: processed_results_csv

    script:
        """
        pairwise_process_output.py pairwise_outfile.txt ${rho_rate} ${sample_size} ${genome_size}
        """

}

process PLOT_RESULTS{
    publishDir "Output/Results", mode: "copy"

    maxForks 1

    input:
        path collectedFile


    output:
        path "rho_comparision.png", emit: rho_comparision_png
        path "max_lk_comparision.png", emit: max_lk_comparision_png

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
    
    // precomputed likelihood table
    // lookup_Table = Channel.fromPath("$baseDir/lookupTable.txt")
    
    // trees = Channel.fromPath("$baseDir/trees.txt")
    // fasta = Channel.fromPath("$baseDir/simbac.fasta")

    rho_rates = Channel.from(0.1) // For fastsimbac use this for recom rate (it doesn't accept rho)
    sample_sizes = Channel.from(30)
    genome_sizes = Channel.from(25000)

    RATE_SELECTOR(rho_rates, sample_sizes, genome_sizes)

    // MS(RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, params.seed, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    FAST_SIM_BAC(RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(MS.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(trees, RATE_SELECTOR.out.path_fn_modifier)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    LDHAT_REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // LDHAT_REFORMAT_FASTA(fasta, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    LOOKUP_TABLE_LDPOP(RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.path_fn_modifier)

    LDHAT_CONVERT(LDHAT_REFORMAT_FASTA.out.ldhat_reformated_fa, RATE_SELECTOR.out.path_fn_modifier)

    SWITCH_TO_GENE_CONVERSION_MODE(LDHAT_CONVERT.out.locs_txt, RATE_SELECTOR.out.path_fn_modifier)

    // LDHAT_INTERVAL(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, LDHAT_CONVERT.out.freqs_txt, LDHAT_CONVERT.out.locs_txt, LDHAT_CONVERT.out.sites_txt, RATE_SELECTOR.out.path_fn_modifier)

    // LDHAT_INTERVAL_STAT(LDHAT_INTERVAL.out.rates_forLDhatStat, RATE_SELECTOR.out.path_fn_modifier)

    LDHAT_PAIRWISE(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, SWITCH_TO_GENE_CONVERSION_MODE.out.locs_C_txt, LDHAT_CONVERT.out.sites_txt, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_PROCESS_OUTPUT(LDHAT_PAIRWISE.out.pairwise_outfile_txt, RATE_SELECTOR.out.rho_rate, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    collectedFile = PAIRWISE_PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    PLOT_RESULTS(collectedFile)
}