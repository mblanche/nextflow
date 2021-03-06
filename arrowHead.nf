#!/usr/bin/env nextflow

params.hicDir = false
params.outDir = false
params.resolutions = false
params.help   = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        arrowHead  --coolDir ~/path/to/coolFiles/location

         or

        arrowHead  --coolDir ~/path/to/coolFiles/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --hicDir [path]             Name of the a direcoty with cool files (can be local or valid S3 location).
    
    Facultative arguments
        --outDir [path]             Path to a diectory to save result files (can be local or valid S3 location).
        --resolutions [integers]    Comma-seperated list of resolutions in kb. Default: [5,10]

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
    outDir = file(params.hicDir).getParent() + "/arrowHead"
} else {
    outDir = params.outDir
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10]
}

Channel
    .fromPath("${params.hicDir}/*.hic")
    .first()
    .set{arrowHead_ch}

process arrowhead {
    cpus 12
    memory '40 GB'
    container "mblanche/juicer"

    publishDir "${outDir}",
	mode: 'copy'
    
    input:
    tuple path(hic), val(res) from arrowHead_ch
	.combine(Channel.from(resolutions))

    output:
    tuple id, path("${id}_${res}kb") into arrowhead_ch

    script:
    id = hic.name.toString().take(hic.name.toString().lastIndexOf('.'))
    bpRes = res.toInteger() * 1000
    """
    mkdir -p ${id}_${res}kb && touch ${id}_${res}kb/${bpRes}_blocks.bedpe
    java -Xmx24000m \
	-jar /juicer_tools.jar \
	arrowhead \
	--threads ${task.cpus} \
	--ignore-sparsity \
	-r ${bpRes} \
	-k KR \
	${hic} \
	${id}_${res}kb
    """
}
