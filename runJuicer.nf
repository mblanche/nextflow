#!/usr/bin/env nextflow

params.expDir = false
params.expName = false

params.validPairDir = false

params.help = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        runJuicer --validPairsDir path/to/Dir  --expDir prodEpi --expName myExpName
    
    Mandatory arguments:
        --expDir [nameOfDir]         Name of the base directory in ref_push where the data will be saved.
        --expName [nameOfDir]        Name of the experiment directory where the data will be saved.
        --validPairsDir [str]        Path to a directoy of fastq files.


    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.expDir && params.expName && params.validPairsDir) ){
    exit 1, "--expDir, --expName and --validPairsDir  are required arguments. Use --help to get the full usage." 
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
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/hicFiles",
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
