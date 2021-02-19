#!/usr/bin/env nextflow

params.expDir = 'prodEpi'
params.expName = 'tmp'
params.noJuicer = false
Channel
    .fromPath("/home/ubuntu/ebs/ref_push/${params.expDir}/${params.expName}/bam/*.bam")
    .map { file -> tuple(file.name.toString().replaceFirst(/.bam/,''),file) }
    .set{bam_ch}

Channel
    .fromFilePairs("/home/ubuntu/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/**.{valid.pairs.gz,valid.pairs.gz.px2}",flat: true)
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

process juicer {
    label 'bigTask'
    cpus 48
    memory "100 GB"
    container 'mblanche/juicer'
    
    publishDir "result/",
    	mode: 'copy'
    
    input:
    tuple val(id), path(chr_sizes), path(pairs), path(idx) from ch2
    
    output:
    path "*.hic" into juicer_out_ch
    path chr_sizes into chrSizes_ch2

    when:
    !params.noJuicer
    
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


