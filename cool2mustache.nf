#!/usr/bin/env nextflow

params.coolDir = false
params.outDir = false
params.resolution = 5
params.help   = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        cool2mustanche  --coolDir ~/path/to/coolFiles/location

         or

        cool2mustanche  --coolDir ~/path/to/coolFiles/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --coolDir [path]             Name of the a direcoty with cool files (can be local or valid S3 location).
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save result files (can be local or valid S3 location).
        --resolution [integer] Resolution in kb to build look for loop. Default: 5kb

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.coolDir) {
        exit 1, "--coolDir is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.coolDir).getParent() + "/mustache"
} else {
    outDir = params.outDir
}


Channel
    .fromPath("${params.coolDir}/*.cool")
    .set{mustache_mcool_ch}

process mustache {
    tag "_${id}"
    cpus 24
    memory '48 GB'
    container "mblanche/mustache"

    publishDir "${params.outDir}",
	mode: 'copy'

    input:
    tuple path(mcool), val(res)  from mustache_mcool_ch
	.combine(Channel.from(1000,4000,16000))
    
    output:
    tuple id, path("*.tsv") into mustache_out_ch

    script:
    id = mcool.name.take(mcool.name.lastIndexOf('.'))
    """
    touch ${id}_${res}kb_loops.tsv 
    mustache -p ${task.cpus} \
	-f ${mcool} \
	-r ${res} \
	-o ${id}_${res}kb_loops.tsv
    """
}

