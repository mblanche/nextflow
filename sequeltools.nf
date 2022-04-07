#!/usr/bin/env nextflow


params.design = false
params.outDir = false
params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        sequeltools  --design ~/path/to/design.txt

         or

        sequeltools  --design ~/path/to/design.txt --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --design [path]              Path to a file with the location of bam files to process, either local or S3.
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

    """.stripIndent()
}

if (!params.design) {
        exit 1, "--design is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = "${PWD}/SequelToolQC/"
} else {
    outDir = params.outDir
}

if (params.help){
    helpMessage()
    exit 0
}

process qc {
    echo true
    label 'sequelQC'
    tag '_${id}'
    cpus 48
    memory '100 GB'
    container 'mblanche/sequeltools'
    
    publishDir "${outDir}",
	mode: "copy"
    
    input:
    path(bam) from Channel
	.fromPath(params.design)
	.splitText()
    
    output:
    tuple id, path("${id}_QC") into out_ch
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    SequelTools.sh -t Q -k -n ${task.cpus} -o ${id}_QC -u <( echo ${bam})
    """
}

