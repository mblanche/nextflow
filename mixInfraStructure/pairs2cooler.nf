#!/usr/bin/env nextflow

params.expDir = 'capture'
params.expName = 'tmp2'
params.noJuicer = false

Channel
    .fromPath("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bam/*.bam")
    .map { file -> tuple(file.name.toString().replaceFirst(/.bam/,''),file) }
    .set{bam_ch}

Channel
    .fromFilePairs("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/**.{valid.pairs.gz,valid.pairs.gz.px2}",flat: true)
    .set { pairs_ch }



process make_chr_size {
    container 'mblanche/bwa-samtools'

    input:
    tuple val(id), path(bam) from bam_ch

    output:
    tuple val(id), path('chr_size.tsv') into chrSizes_ch

    script:
    """
    samtools view -H ${bam} | \
	awk -v OFS='\t' '/^@SQ/{split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' \
	> chr_size.tsv
    """
}

chrSizes_ch
    .join(pairs_ch)
    .set { ch2 }

process cooler_cload {
    tag "_${id}"
    label 'batch'
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles2/${id}",
    	mode: 'copy'
    
    input:

    tuple val(id), path(chr_sizes), path(pairs), path(idx) from ch2
    
    output:
    path "*.cool" into cooler_ch

    script:
      """
    cooler cload pairix \
	-p ${task.cpus} \
	${chr_sizes}:1000 \
	${pairs} \
	${id}.cool
    """
}

process cooler_zoomify {
    tag "_${id}"
    label 'batch'
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles2/${id}",
    	mode: 'copy'
    
    input:
    path cooler from cooler_ch

    output:
    path "*.mcool" into mcool_ch
    
    script:
    id = cooler.name.toString().tokenize('.').get(0)
    """
    cooler zoomify --balance -p ${task.cpus} ${cooler}
    """
}
