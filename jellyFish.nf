#!/usr/bin/env nextflow

params.fastqDir = false
params.outDir = false
params.help   = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        jellyfish.nf --coolDir ~/path/to/coolFiles/location --outDir ~/path/to/where/to/save/files
    
    Mandatory arguments:
    
        --fastqDir [path]             Name of the a direcoty with cool files (can be local or valid S3 location).
        --outDir [path]             Path to a diectory to save result files (can be local or valid S3 location).

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.fastqDir && params.outDir)) {
        exit 1, "--fastqDir and --outDir are required arguments. Use --help to get the full usage." 
}


process jellyfish {
    label 'cpu'
    container "mblanche/jellyfish"
    
    publishDir "${params.outDir}/jellyfish",
	mode: 'copy'
    
    input:
    tuple id, path(fq) from Channel.fromFilePairs("${params.fastqDir}/*_{1,2}.fastq.gz")
    
    output:
    tuple id, path("*.jf") into histo_ch
    
    script:
    """
    jellyfish count -C -m 21 -s 1000M -t ${task.cpus}  <(zcat ${fq})  -o ${id}.jf
    """
}


process jellyfish_histo {
    label 'cpu'
    container "mblanche/jellyfish"
    
    publishDir "${params.outDir}/jf_histo",
	mode: 'copy'
    
    input:
    tuple id, path(jf) from histo_ch
    
    output:
    path("*.histo")
    
    script:
    """
    jellyfish histo -t ${task.cpus}  ${jf} > ${id}.histo
    """
}
