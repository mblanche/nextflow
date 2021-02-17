#!/usr/bin/env nextflow

params.expDir = 'capture'
params.expName = 'tmp'
params.genome = 'hg38'

Channel
    .fromPath("/mnt/ebs/ref_push/${params.expDir}/${params.expName}/bam/*.bam")
    .map { file -> tuple(file.name.toString().replaceFirst(/.bam/,''),file) }
    .set{bam_ch}

Channel
    .fromFilePairs("/mnt/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/**.{valid.pairs.gz,valid.pairs.gz.px2}",flat: true)
    .set { pairs_ch }

bam_ch
    .join(pairs_ch)
    .view()
    .into {ch1;ch2}

process make_chr_size {
    container 'mblanche/bwa-samtools'

    
    input:
    tuple val(id), path(bam), path(pairs), path(idx) from ch1


    output:
    path 'chr_size.tsv' into chrSizes_ch

    script:
    """
    samtools view -H ${bam} | \
	awk -v OFS='\t' '/^@SQ/{split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' \
	> chr_size.tsv
    """
}


process juicer {
    echo true
    cpus 48
    memory '100 GB'
    container 'mblanche/juicer'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/hicFiles",
    	mode: 'copy'
    
    input:
    tuple val(id), path(bam), path(pairs), path(idx) from ch2
    path chr_sizes from chrSizes_ch
    
    output:
    path "*.hic" into juicer_out_ch
    
    script:
    """
    java -Xmx16000m -Djava.awt.headless=true \
	-jar /juicer_tools_1.22.01.jar pre \
	--threads ${task.cpus} \
	-j ${task.cpus} \
	-k VC,VC_SQRT,KR,SCALE \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    """
}


