#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                  Chromatin-state discovery pipeline   
========================================================================================

@#### Authors
Federica Gervasoni <federica.gervasoni@ifom.eu>
Raoul JP Bonnal <raoul.bonnal@ifom.eu>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:

- STEP1: Binarization process
- STEP2: Model learning
- STEP3: Compare models
----------------------------------------------------------------------------------------

Commmand line example:
	nextflow run main_chromHMM.nf --align "alignfiles_chromHMM.txt" --chromHMMconfig 'cellmarkfiletable_organoids.txt' --genome 'hg38' --numstates '8,10,12' -w work 

*/

params.align="alignfiles_chromHMM.txt"

chromHMMconfig = file(params.chromHMMconfig)

bam_to_chromHMM = Channel.create()

Channel
    .fromPath(params.align)
    .splitCsv(header: true)
    .map { sample -> 
    println("Read data to chromHMM")
    println([sample.path])
    [file(sample.path)] }
    .into {bam_to_chromHMM}

/*
 * STEP1: Binarization process
 */

process binarize {
    cpus 4
    memory '10 GB'

    input:
    file '*' from bam_to_chromHMM.collect()

    output:
    file 'OUTDIR/*' into in_learnmodel

    script:
    """
    java -Xmx4g -jar ${params.chromhmm.jar} BinarizeBam ${params.chromhmm.chromsizes}/hg38.txt . ${chromHMMconfig} OUTDIR
    """
}

numStates = params.numstates.toString().tokenize(', ')

/*
 * STEP2: Model learning
 */

process learnmodel {
    cpus 4
    memory '10 GB'
    publishDir "${params.outdir}/chromhmm/", mode: 'copy'

    input:
    file '*' from in_learnmodel.collect()
    val(numstates) from numStates

    output:
    file 'OUTDIR/*' into learned_models
    file("OUTDIR/emissions_${numstates}.txt") into emissions_files

    script:
    """
    java -Xmx4g -jar ${params.chromhmm.jar} LearnModel -p ${task.cpus} -gzip . OUTDIR ${numstates} hg38
    """
}

in_comparemodels = ( numStates.size() > 1 ? emissions_files : Channel.empty() )

/*
 * STEP3: Compare models
 */

process comparemodels {
    cpus 4
    memory '10 GB'
    publishDir "${params.outdir}/chromhmm/comparison", mode: 'copy'

    input:
    file '*' from in_comparemodels.collect()

    output:
    file("${maxStates}_comparemodels*") into states_comparison

    script:
    maxStates = numStates.max()
    """
    java -Xmx4g -jar ${params.chromhmm.jar} CompareModels emissions_${maxStates}* . ${maxStates}_comparemodels
    """
}















    








































