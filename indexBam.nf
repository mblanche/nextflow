#!/usr/bin/env nextflow

params.bamDir = false

params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        indexBam  --bamDir ~/path/to/bam/location
    
    Mandatory arguments:
        --bamDir [path]              Name of the a direcoty with bam files and their index to compute coverages (can be local or valid S3 location.
    
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.bamDir) {
        exit 1, "--bamDir is a required arguments. Use --help to get the full usage." 
}


process index {
    label 'index'
    tag '_${id}'
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "${params.bamDir}",
	mode: "copy"

    input:
    path(bam) from Channel.fromPath("${params.bamDir}/*.bam")

    output:
    path("*.bam.bai")
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    samtools index -@${task.cpus} ${bam}
    """
    

}
