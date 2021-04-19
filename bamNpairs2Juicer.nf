#!/usr/bin/env nextflow

params.expDir = 'working'
params.expName = 'DHockemeyer_cellLine'

patt = 'CTTCCTTC'

Channel
    .fromPath("/mnt/ebs/ref_push/${params.expDir}/${params.expName}/bam/*.bam")
    .map { file -> tuple(file.name.toString().replaceFirst(/.bam/,''),file) }
    .set{bam_ch}

Channel
    .fromFilePairs("/mnt/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/**.{valid.pairs.gz,valid.pairs.gz.px2}",flat: true)
    .set{pairs_ch}

bam_ch
    .join(pairs_ch)
    .set{ch1}

process make_chr_size {
    container 'mblanche/bwa-samtools'

    input:
    tuple val(id), path(bam), path(pairs), path(idx) from ch1
    
    output:
    tuple id, path('chr_size.tsv'), path(pairs), path(idx) into cooler_sort_ch, juicer_ch
    
    script:
    """
    samtools view -H ${bam} | \
	awk -v OFS='\t' '/^@SQ/ && !(\$2 ~ /:(chr|"")M/) {split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' | \
	sort -V -k1,1 \
	> chr_size.tsv
    """
}

process cooler_sort {
    echo true
    tag "_${id}"
    label "batch"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    input:
    tuple id, path(chr_sizes), path(pairs), path(idx) from cooler_sort_ch
    
    output:
    tuple id, path(chr_sizes), path("*.gz"), path("*.px2") into pairs_ch_cooler
    
    script:
    """
    cooler csort \
	-c1 2 -p1 3 -c2 4 -p2 5 \
	-i pairix \
	-p ${task.cpus} \
	--out ${id}_sorted.valid.pairs.gz \
	${pairs} \
	${chr_sizes}
    """
}

process cooler_cload {
    echo true
    tag "_${id}"
    label "batch"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
	mode: 'copy'
    
    input:
    tuple id, path(chr_sizes), path(pairs), path(idx) from pairs_ch_cooler
    
    output:
    tuple id, path("*.cool") into cooler_ch

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
    label "batch"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
    	mode: 'copy'
    
    input:
    tuple id, path(cooler) from cooler_ch

    output:
    path "*.mcool" into mcool_ch
    
    script:
    """
    cooler zoomify --balance -p ${task.cpus} ${cooler}
    """
}


process juicer {
    tag "_${id}"
    label "batch"
    cpus 24
    memory '24 GB'
    container 'mblanche/juicer'

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/hicFiles",
    	mode: 'copy'

    input:
    tuple id, path(chr_sizes), path(pairs), path(idx) from juicer_ch
    
    output:
    path "*.hic" into juicer_out_ch

    script:
    """
        java -Xmx24000m -Djava.awt.headless=true \
	-jar /juicer_tools_1.22.01.jar pre \
	--threads ${task.cpus} \
	-j ${task.cpus} \
	-k VC,VC_SQRT,KR,SCALE \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    """
}
