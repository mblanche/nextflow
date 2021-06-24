#!/usr/bin/env nextflow

params.validPairsDir = false
params.outDir = false

params.help = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        runJuicer  --validPairsDir ~/path/to/.valid.pairs.gz

         or
        runJuicer  --validPairsDir ~/path/to/.valid.pairs.gz --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --validPairsDir [path]      path to a direcoty with .valid.pairs.gz files (can be local or valid S3 location).
    
    Facultative arguments
        --outDir [path]             Path to a diectory to save result files (can be local or valid S3 location).

    runJuicer.ng will covert every file ending in .valid.pairs.gz in the specified directory to a .hic file. If no --outDir option is specified,
    the files will be saved in the parent direcoty of --validPairsDir in a direcoty name hicFiles
    
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.validPairsDir) ){
    exit 1, "--validPairsDir is required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.validPairsDir).getParent() + "/hicFiles"
} else {
    outDir = params.outDir
}

Channel
    .fromPath("${params.validPairsDir}/*.valid.pairs.gz",
	      checkIfExists: true)
    .set{pairs_ch}


process chrm_size {
    echo true
    tag "_${id}"
    cpus 4
    memory '24 GB'
    container 'mblanche/pairtools'
    
    input:
    path(pairs) from pairs_ch
    
    output:
    tuple id, path(pairs), path("*.px2"), path("*.tsv") into juicer_ch
    
    script:
    id = pairs.name.toString().take(pairs.name.toString().indexOf('.',1))
    """
    pairix ${pairs}

    pairix -H -f ${pairs} \
	| awk -v OFS='\t' '/^#chromsize/  {print \$2,\$3}' \
	| sort -V -k1,1 \
	> chr_size.tsv
    """
}

process juicer {
    tag "_${id}"
    cpus 24
    memory '150 GB'
    container 'mblanche/juicer'
    
    publishDir "${outDir}",
	mode: 'copy'
    
    input:
    tuple id, path(pairs), path(idx), path(chr_sizes)  from juicer_ch
    
    output:
    tuple id, path("*.hic") into arrowHead_ch
 
    script:
    """
    java -Xmx96000m -Djava.awt.headless=true \
	-jar /juicer_tools_1.22.01.jar pre \
	--threads ${task.cpus} \
	-j ${task.cpus} \
	-k VC,VC_SQRT,KR,SCALE \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    """
}
