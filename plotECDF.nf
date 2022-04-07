#!/usr/bin/env nextflow

params.bam = false
params.outDir = false
params.prefix = 'myResult'
params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        plotECDF  --bam ./path/to/lib1,./path/to/lib2,./path/to/lib3,... --outDir path/to/results/

         or

        plotECDF --prefix myResult --bam ./path/to/lib1,./path/to/lib2,./path/to/lib3,... --outDir path/to/results/

	or
        
        plotECDF --bam \$(find \${PWD}/to/bam/ -name '*.bam' |xargs |tr ' ' ,) --outDir path/to/results/
	
    Mandatory arguments:
        --bam [path]             Comma delimited list of paths to bam files
        --outDir [path]          Path to a diectory to save the bigwig coveage files (can be local or valid S3 location. Default 
    
    Facutativ argument:
        --prefix [string]        A string to name the _ecdf.pdf file. Default: myResult

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.bam && params.outDir)) {
    exit 1, "--bam and --outDir area required arguments. Use --help to get the full usage." 
}

process index {
    label 'index'
    tag '_${id}'
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    path(bam) from Channel
	.from(params.bam)
	.splitCsv()
	.flatten()
    
    output:
    path(bam) into bam_ch
    path("*.bam.bai") into idx_ch
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    samtools index -@${task.cpus} ${bam}
    """
}
process plotECDF {
    echo true
    label 'batch'
    tag "_${id}"
    cpus 24
    memory '48 GB'
    container 'mblanche/r-prod-qc'
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    path(bam) from bam_ch.collect()
    path(idx) from idx_ch.collect()
    
    output:
    path "*.pdf"
    
    script:
    """
    plotECDF.R ${task.cpus} ${params.prefix} ${bam}
    """
}

