#!/usr/bin/env nextflow

params.mcoolDir = false
params.outDir = false
params.resolutions = false
params.help   = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        mustache  --mcoolDir ~/path/to/coolFiles/location

         or

        mustache  --mcoolDir ~/path/to/coolFiles/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --mcoolDir [path]           Name of the a direcoty with mcool files (can be local or valid S3 location).
    
    Facultative arguments
        --outDir [path]             Path to a diectory to save result files (can be local or valid S3 location).
        --resolutions [integers]    Comma-seperated list of resolutions in kb. Default: [5,10]

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.mcoolDir) {
        exit 1, "--mcoolDir is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.mcoolDir).getParent() + "/mustache"
} else {
    outDir = params.outDir
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10]
}

Channel
    .fromPath("${params.mcoolDir}/*.mcool")
    .set{mustache_mcool_ch}


process mustache {
    cpus 24
    memory '48 GB'
    container "mblanche/mustache"
    
    publishDir "${outDir}",
	mode: 'copy'
    
    input:
    tuple path(mcool), val(res)  from mustache_mcool_ch
	.combine(Channel.from(1000,4000,16000))
    
    output:
    tuple id, path("*.tsv") into mustache_2_merge_ch
    
    script:
    id = mcool.name.toString().take(mcool.name.toString().lastIndexOf('.'))
    """
    mustache -p ${task.cpus} \
	-f ${mcool} \
	-r ${res} \
	-o ${id}_${res}kb_loops.tsv
    """

}
