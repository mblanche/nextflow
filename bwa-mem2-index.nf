#!/usr/bin/env nextflow

params.design = false

if ( !(params.design && params.outDir) ){
    exit 1, "--design and --outDir are required arguments"
}

process index_genome {
    label 'large'
    
    container 'mblanche/bwa-mem2'
    
    publishDir "bwa-mem2-index/${id}",
	mode: 'copy'
    
    input:
    tuple id, path(fasta) from Channel
	.fromPath(params.design)
	.splitCsv(header: false)
    
    output:
    tuple val(id), path(fasta), path("*") into bwa_index
    
    script:
    """
    bwa-mem2 index -p ${id} ${fasta}
    """
}
