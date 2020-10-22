#! /usr/bin/env nextflow

params.genomeSize = '100000'
params.meanFragmentLen = '150'
params.sampleSize = '100'

// Rates the user wishes to change
recom_rates = [0.001]
mutation_rates = [0.001]
seed_values = [1]

// likelihood table
lt_file = Channel.fromPath("$baseDir/lk_n100_t0.001") // n=100, theta=0.001 per site

process fastSimBac{
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$seed"+"_m_"+"$mutation_rate"+"_r_"+"$recom_rate"+"/"+
                                                "s_"+"$seed"+"_m_"+"$mutation_rate"+"_r_"+"$recom_rate"+"_"+"$filename" }

    maxForks 1 // Run sequentially
    time '2h' // Should fastsimbac freeze, skip
    errorStrategy 'ignore'// Should fastsimbac freeze, skip. Without this line program will stop
    
    input:
        each recom_rate from recom_rates
        each mutation_rate from mutation_rates
        each seed from seed_values

    output:
        file "trees.txt" into trees
        val seed into sOut1, sOut2, sOut3, sOut4, sOut5, sOut6, sOut7, sOut8, sOut9, sOut10, sOut11
        val mutation_rate into mOut1, mOut2, mOut3, mOut4, mOut5, mOut6, mOut7, mOut8, mOut9, mOut10, mOut11
        val recom_rate into rOut1, rOut2, rOut3, rOut4, rOut5, rOut6, rOut7, rOut8, rOut9, rOut10, rOut11
        
    script:
        """
        fastSimBac "${params.sampleSize}" "${params.genomeSize}" -s "${seed}" -T -t "${mutation_rate}" -r "${recom_rate}" 500 > trees.txt
        """
}


process cleanTrees{
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$sOut1"+"_m_"+"$mOut1"+"_r_"+"$rOut1"+"/"+
                                                "s_"+"$sOut1"+"_m_"+"$mOut1"+"_r_"+"$rOut1"+"_"+"$filename" }

    maxForks 1

    input:
        file trees
        val sOut1
        val mOut1
        val rOut1

    output:
        file "cleanTrees.txt" into cleanTrees

    script:
        """
        cleanTrees.py trees.txt
        """
}

// Using all trees
process seqGen{
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$sOut2"+"_m_"+"$mOut2"+"_r_"+"$rOut2"+"/"+
                                                "s_"+"$sOut2"+"_m_"+"$mOut2"+"_r_"+"$rOut2"+"_"+"$filename" }

    maxForks 1

    input:
        file cleanTrees
        val sOut2
        val mOut2
        val rOut2

    output:
        file "seqgenOut.fa" into seqgenOut

    script:
        // 1 partiion per tree
        // program crashes if seq length is not as the one set for fastsimbac
        """
        numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
        seq-gen -m HKY -t 4 -l "${params.genomeSize}" -z "${sOut2}" -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
        """
}


process reformatFasta{
    // reformats headers of fasta enteries into a format suitable for LDhat
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$sOut3"+"_m_"+"$mOut3"+"_r_"+"$rOut3"+"/"+
                                                "s_"+"$sOut3"+"_m_"+"$mOut3"+"_r_"+"$rOut3"+"_"+"$filename" }

    maxForks 1

    input:
        file seqgenOut
        val sOut3
        val mOut3
        val rOut3

    output:
        file "LDhat_reformated.fa" into fasta_forLDhatConvert

    script:
        //  Assumption setting it to haploid (1) causes the estimator to use 2ner

        // Additional notes in python script
        """
        LDhat_reformatFasta.py seqgenOut.fa "${params.sampleSize}" "${params.genomeSize}" 1
        """
}


process LDhat_convert{
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$sOut4"+"_m_"+"$mOut4"+"_r_"+"$rOut4"+"/"+
                                                "s_"+"$sOut4"+"_m_"+"$mOut4"+"_r_"+"$rOut4"+"_"+"$filename" }

    maxForks 1

    input:
        file fasta_forLDhatConvert
        val sOut4
        val mOut4
        val rOut4

    output:
        file "freqs.txt" into freqs_forLDhatInterval
        file "locs.txt" into locs_forLDhatInterval
        file "sites.txt" into sites_forLDhatInterval

    script:

        """
        convert -seq LDhat_reformated.fa
        """

}


process LDhat_interval{
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$sOut5"+"_m_"+"$mOut5"+"_r_"+"$rOut5"+"/"+
                                                "s_"+"$sOut5"+"_m_"+"$mOut5"+"_r_"+"$rOut5"+"_"+"$filename" }

    maxForks 1

    input:
        file lt from lt_file
        file freqs_forLDhatInterval
        file locs_forLDhatInterval
        file sites_forLDhatInterval
        val sOut5
        val mOut5
        val rOut5

    output:
        file "bounds.txt" into bounds_forLDhatStat
        file "rates.txt" into rates_forLDhatStat
        file "new_lk.txt"
        file "type_table.txt"

    script:
        // pre-generated lookup tables available from the ldhat github page was used. One that matches number of samples and mutation rate was selected.
        // the arguments -its, -samp, -bpen use recommened values given in the manual for the interval program.
        """
        interval -seq sites.txt -loc locs.txt -lk lk_n100_t0.001 -its 1000000 -samp 2000 -bpen 5
        """

}


process LDhat_stat{
    publishDir "Output", mode: "copy", saveAs: { filename -> "s_"+"$sOut6"+"_m_"+"$mOut6"+"_r_"+"$rOut6"+"/"+
                                                "s_"+"$sOut6"+"_m_"+"$mOut6"+"_r_"+"$rOut6"+"_"+"$filename" }

    maxForks 1

    input:
        file rates_forLDhatStat
        val sOut6
        val mOut6
        val rOut6

    output:
        file "res.txt"
        file "statOut.txt"

    script:
        // The information printed on screen was useful so decided to save that also.
        """
        stat -input rates.txt > statOut.txt
        """

}

