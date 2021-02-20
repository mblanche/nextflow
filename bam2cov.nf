#!/usr/bin/env nextflow
/*
* TODO: Update this for the current script
*/


params.expDir = 'capture'
params.expName = 'Twist-panGenome-2'

Channel
    .fromFilePairs("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bam/*.{bam,bam.bai}", flat: true)
    .set{bam_ch}

process deeptools_bw {
    tag "dt_cov_${id}"
    echo true
    cpus 48
    memory '100 GB'
    container 'mblanche/deeptools'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bigwigs",
    	mode: 'copy'
    
    input:
    tuple val(id), path(bam), path(idx) from bam_ch
    
    output:
    path "*.bw" into bw_out_ch

    script:
    """
    bamCoverage -p ${task.cpus} -bs 1 -b ${bam} -o ${id}.bw
    """
}

