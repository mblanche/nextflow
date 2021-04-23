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
    .fromPath("${params.coolDir}/**.cool")
    .into{cool_ch;cool_ch2}

process get_chr {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    input:
    path cool from cool_ch

    output:
    stdout into chr_ch
    
    script:
    id = cool.name.take(cool.name.lastIndexOf('.'))
    """
    cooler dump -t chroms -c name ${cool} |grep -v chrY | awk -v f=${id} -v OFS=, '{print f,\$0}'
    """
}

process mustache {
    cpus 24
    memory '48 GB'
    container "mblanche/mustache"

    chr_ch
	.splitText()
	.splitCsv()
	.set{left}
    
    cool_ch2
	.map{cool ->
	    id = cool.name.take(cool.name.lastIndexOf('.'))
	    return(tuple(id,cool))
	}
	.cross(left)
	.map{file,chr ->
	    return(tuple(file[0],file[1],chr[1]))
	}
	.set{results}
    
    input:
    tuple id, path(cool), val(chr) from results
    
    output:
    tuple id, path("*.tsv") into mustache_2_merge_ch
    
    script:
    """
    mustache -p ${task.cpus} \
	-f ${cool} \
	-r ${params.resolution}kb \
	-ch ${chr} \
	-o ${id}_${chr}_${params.resolution}kb_loops.tsv
    """
}



process merge_mustache {
    echo true
    cpus 4
    memory '8 GB'
    container "ubuntu:20.04"
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple id, path(mustache) from mustache_2_merge_ch
	.groupTuple()
    
    output:
    tuple id, path("*.tsv") into loops_ch
    
    script:
    """
    cat <(head -n 1 ${mustache[0]})  <(tail -q -n +2  ${mustache}) > ${id}_${params.resolution}kb_loops.tsv
    """
}
//##
