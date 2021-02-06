#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process FAST_SIM_BAC{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1 // Run sequentially
    time '2h' // Should fastsimbac freeze, skip
    // errorStrategy 'ignore'// Should fastsimbac freeze, skip. Without this line program will stop

    echo true
    
    input:
        each recom_rate
        each mutation_rate
        each seed

    output:
        path "trees.txt", emit: trees_txt
        val "${recom_rate}", emit: r_val
        val "${mutation_rate}", emit: m_val
        val "${seed}", emit: s_val
             
    script:
    """
    fastSimBac ${params.sampleSize} ${params.genomeSize} -s ${seed} -T -t ${mutation_rate} -r ${recom_rate} ${params.recom_tract_len} > trees.txt
    """
}


process CLEAN_TREES{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}
    maxForks 1

    input:
        path trees
        val recom_rate
        val mutation_rate
        val seed

    output:
        path "cleanTrees.txt", emit: cleanTrees_txt

    script:
    """
    cleanTrees.py trees.txt 
    """
}


process SEQ_GEN{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1

    input:
        path cleanTrees
        val recom_rate
        val mutation_rate
        val seed

    output:
        path "seqgenOut.fa", emit: seqgenout_fa

    script:
    // 1 partiion per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${params.genomeSize} -z ${seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process LDHAT_REFORMAT_FASTA{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1

    input:
        path seqgenOut
        val recom_rate
        val mutation_rate
        val seed

    output:
        path "LDhat_reformated.fa", emit: ldhat_reformated_fa

    script:
    """
    LDhat_reformat_fasta.py seqgenOut.fa "${params.sampleSize}" "${params.genomeSize}" 1
    """
}


process LOOKUP_TABLE_LDPOP{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1 
    
    input:
        val recom_rate
        val mutation_rate
        val seed

    output:
        path "ldtable_log.txt", emit: log_txt
        path "lookupTable.txt", emit: lookupTable_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores 4 -n ${params.sampleSize} -th ${mutation_rate} -rh ${params.ldpop_rho_range} --approx --log ldtable_log.txt > lookupTable.txt
    """
}


process LDHAT_CONVERT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1

    input:
        path fasta_forLDhatConvert
        val recom_rate
        val mutation_rate
        val seed

    output:
        path "freqs.txt", emit: freqs_txt
        path "locs.txt", emit: locs_txt
        path "sites.txt", emit: sites_txt
        path "convertOut.txt", emit: convert_out_txt

    script:
        // The information printed on screen was useful so decided to save that also.
        // -2only, only output sites with exactly two alleles
        """
        convert -seq  LDhat_reformated.fa -2only > convertOut.txt
        """

}


process LDHAT_INTERVAL{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1

    input:
        path lookup_table_file
        path freqs_forLDhatInterval
        path locs_forLDhatInterval
        path sites_forLDhatInterval
        val recom_rate
        val mutation_rate
        val seed

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


process LDHAT_STAT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "s_${seed}_m_${mutation_rate}_r_${recom_rate}/s_${seed}_m_${mutation_rate}_r_${recom_rate}_${filename}"}

    maxForks 1

    input:
        path rates_forLDhatStat
        val recom_rate
        val mutation_rate
        val seed

    output:
        path "res.txt" 
        path "statOut.txt"

    script:
        // The information printed on screen was useful so decided to save that also.
        """
        stat -input rates.txt > statOut.txt
        """

}

// Params and Channels for workflow
// Note: Channels can be called unlimited number of times in DSL2
params.genomeSize = '250000'
params.meanFragmentLen = '150'
params.sampleSize = '50'
params.recom_tract_len = '500'
params.ldpop_rho_range = '10,100'

// precomputed likelihood table
// lookup_Table = Channel.fromPath("$baseDir/lookupTable.txt")

// recom_rates = Channel.from(0,0.0001,0.001,0.01,0.1)
// mutation_rates = Channel.from(0.001, 0.01, 0.1)
// seed_values = Channel.from(1,2,3,4,5,6,7,8,9,10)

recom_rates = Channel.from(0.1)
mutation_rates = Channel.from(0.1)
seed_values = Channel.from(12345)

workflow {
    // A process component can be invoked only once in the same workflow context.
    FAST_SIM_BAC(recom_rates, mutation_rates, seed_values)

    CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)

    LDHAT_REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)

    LOOKUP_TABLE_LDPOP(FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)

    LDHAT_CONVERT(LDHAT_REFORMAT_FASTA.out.ldhat_reformated_fa, FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)

    LDHAT_INTERVAL(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, LDHAT_CONVERT.out.freqs_txt, LDHAT_CONVERT.out.locs_txt, LDHAT_CONVERT.out.sites_txt, FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)

    LDHAT_STAT(LDHAT_INTERVAL.out.rates_forLDhatStat, FAST_SIM_BAC.out.r_val, FAST_SIM_BAC.out.m_val, FAST_SIM_BAC.out.s_val)
}