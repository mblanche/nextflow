#!/usr/bin/env nextflow

params.bamDir = false
params.outDir = false
params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        plotECDF  --bam ./path/to/lib1,./path/to/lib2,./path/to/lib3,... 

         or

        bam2bw  --prefix myECDF --bam ./path/to/lib1,./path/to/lib2,./path/to/lib3,... --outDir path/to/results.pdf

    Mandatory arguments:
        --bamDir [path]              Name of the a direcoty with bam files and their index to compute coverages (can be local or valid S3 location.
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.bamDir) {
        exit 1, "--bamDir is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.bamDir).getParent() + "/bigwigs"
} else {
    outDir = params.outDir
}



Channel
    .fromPath("${params.bamDir}/*.bam")
    .set{bam_ch}

process index {
    label 'index'
    tag '_${id}'
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "${params.bamDir}",
	mode: "copy"

    input:
    path(bam) from bam_ch

    output:
    tuple id, path(bam), path("*.bam.bai") into bamNidx_ch
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    samtools index -@${task.cpus} ${bam}
    """
    

}

process bam2bw {
    label 'batch'
    tag "_${id}"
    cpus 20
    memory '150 GB'
    container 'mblanche/r-cov'
    
    publishDir "${outDir}",
	mode: 'copy'
    
    input:
    tuple val(id), path(bam), path(idx) from bamNidx_ch
    
    output:
    path "*.bw" into bigwig_ch
    
    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}

