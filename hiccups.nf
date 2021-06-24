#!/usr/bin/env nextflow

params.hicDir = false
params.outDir = false
params.resolutions = false
params.help   = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        hiccups  --coolDir ~/path/to/coolFiles/location

         or

        hiccups  --coolDir ~/path/to/coolFiles/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --hicDir [path]             Name of the a direcoty with cool files (can be local or valid S3 location).
    
    Facultative arguments
        --outDir [path]             Path to a diectory to save result files (can be local or valid S3 location).
        --resolutions [integers]    Comma-seperated list of resolutions in kb. Default: [5,10,25]

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.hicDir) {
        exit 1, "--hicDir is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.hicDir).getParent() + "/hiccups"
} else {
    outDir = params.outDir
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10,25]
}

Channel
    .fromPath("${params.hicDir}/*.hic")
    .set{hiccups_ch}

process hiccups {
    label 'gpu'
    accelerator 1
    cpus 6
    memory '30 GB'
    container "mblanche/hiccups-gpu"

    publishDir "${outDir}",
	mode: 'copy'
    
    input:
    tuple path(hic), val(res)  from hiccups_ch.first()
	.combine(Channel.from(resolutions.collect{it*1000}.join(',')))
    
    output:
    tuple id, path("${id}_loops") into hiccups_out_ch
    

    script:
    id = hic.name.toString().take(hic.name.toString().lastIndexOf('.'))
    """
    java -Xmx24000m \
	-jar /juicer_tools.jar \
	hiccups \
	--threads ${task.cpus} \
	--ignore-sparsity \
	-m 500 \
	-r ${res} \
	-k KR \
	${hic} \
	${id}_loops
    """
}
