#!/usr/bin/env nextflow

params.expDir = 'capture'
params.expName = 'Twist-PG-newTaq'

Channel
    .fromFilePairs("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bam/*.{bam,bam.bai}",flat: true)
    .set{bam_ch}


process bam2bw {
    label 'batch'
    tag "_${id}"
    cpus 20
    memory '100 GB'
    
    container 'mblanche/r-cov'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bigwigs",
    mode: 'copy'
    
    input:
    tuple val(id), path(bam), path(idx) from bam_ch
    
    output:
    path "*.bw" into bigwig_ch
    
    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}

