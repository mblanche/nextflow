#!/usr/bin/env nextflow

/*
* TODO: Update this for the current script
*/


params.expDir = 'capture'
params.expName = 'tmp'
params.genome = 'hg38'


Channel
    .fromFilePairs("/mnt/ebs/ref_push/${params.expDir}/${params.expName}/fastqs/*_R{1,2}*.fastq.gz", flat: true)
    .map { prefix, file1, file2 -> tuple(prefix.split('_')[0], file1, file2) }
    .groupTuple()
    .set{fastqs_ch}

Channel
    .fromPath("/mnt/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
    .ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    .set { bwa_index }

process bwa_mem {
    cpus 48
    memory '100 GB'
    
    container 'mblanche/bwa-samtools'
    
    input:
    set id, file(R1s), file(R2s) from fastqs_ch
    file index from bwa_index.first()
    
    output:
    path "${id}.sam" into sam4chrSize,
	sam4pt
    
    script:
    """
    bwa mem -5SP -t ${task.cpus} \
	${index}/${params.genome} \
	<(zcat ${R1s}|head -n 100000) \
	<(zcat ${R2s}|head -n 100000) \
	>${id}.sam
    """
}

process make_chr_size {
    container 'mblanche/awscli'

    input:
    path sam from sam4chrSize

    output:
    path 'chr_size.tsv' into chrSizes_pt,
	chrSizes_cooler, chrSizes_juicer
    
    """
    awk -v OFS='\t' '/^@SQ/{split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' ${sam} > chr_size.tsv
    """
}

process pairtools_parse {
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    path sam from sam4pt
    path chr_sizes from chrSizes_pt
    
    output:
    path "*.pairsam" into pairsam_ch

    script:
    id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    pairtools parse \
	--min-mapq 40 \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in \$(nproc) --nproc-out \$(nproc) \
	--chroms-path ${chr_sizes} \
	${sam} \
	> ${id}.pairsam
    """
}

process pairtools_sort {
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    path sam from pairsam_ch
    
    output:
    path "*_sorted.pairsam" into sorted_ps_ch

    script:
    id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    pairtools sort --nproc \$(nproc) $sam >${id}_sorted.pairsam
    """
}

process pairtools_dedup {
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    path sam from sorted_ps_ch
    
    output:
    path "*_dedup.pairsam" into dedup_ps_ch
    path "*_unmapped.pairsam" into unmapped_ps_ch
    path "*_pairtools.stats" into ps_stats
    
    script:
    id = sam.name.toString().tokenize('_').get(0)
    """
    pairtools dedup --nproc-in \$(nproc) --nproc-out \$(nproc) \
	--mark-dups \
	--output-stats ${id}_pairtools.stats  \
	--output ${id}_dedup.pairsam \
	--output-unmapped ${id}_unmapped.pairsam \
	${sam}
    """
}

process pairtools_split_dedup {
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'
        
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}",
    	mode: 'copy',
    	saveAs: {filename ->
    	if (filename.endsWith('.valid.pairs.gz')) "coolerFiles/${id}/$filename"
	else if (filename.endsWith('.valid.pairs.gz.px2')) "coolerFiles/${id}/$filename"
    }
    
    input:
    path sam from dedup_ps_ch
    
    output:
    path "*.bam" into pt_bam_ch, samtools_bam_ch
    path "*.valid.pairs.gz" into pairs_ch_cooler, pairs_ch_juicer
    path "*.px2" into pairs_px2_cooler, pairs_px2_juicer

    script:
    id = sam.name.toString().tokenize('_').get(0)
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_PT.bam  \
	--output-pairs ${id}.valid.pairs.gz  \
	${sam}
    pairix ${id}.valid.pairs.gz
    """
}

process bam_sort {
    cpus 48
    memory '150 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/bam", mode: 'copy'
    
    input:
    path bam from samtools_bam_ch
    
    output:
    path "*.bam" into bam_forindex_ch
    
    script:
    id = bam.name.toString().replaceFirst(/_PT.bam/,'')
    """
    samtools sort -m 2G \
	-@ ${task.cus} \
	-o ${id}.bam \
	${bam}
    """
}

process bam_index {
    cpus 48
    memory '150 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/bam", mode: 'copy'
    
    input:
    path bam from bam_forindex_ch
    
    output:
    path "*.bai" into bam_samtools_index_ch
    
    script:
    id = bam.name.toString().tokenize('.').get(0)
    """
    samtools index -@${task.cpus} ${bam}
    """
}

process pairtools_split_unmapped {
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'
        
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/unmapped",
    	mode: 'copy'

    input:
    path sam from unmapped_ps_ch
    
    output:
    path "*_unmapped.bam" into unmapped_bam_ch
    path "*_unmapped.valid.pairs.gz" into unmapped_pairs_ch

    script:
    id = sam.name.toString().tokenize('_').get(0)
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_unmapped.bam  \
	--output-pairs ${id}_unmapped.valid.pairs.gz  \
	${sam}
    """
}

process cooler_cload {
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
    	mode: 'copy'
    
    input:
    path pairs from pairs_ch_cooler
    path idx from pairs_px2_cooler
    path chr_sizes from chrSizes_cooler
    
    output:
    path "*.cool" into cooler_ch
    
    script:
    id = pairs.name.toString().tokenize('.').get(0)
    """
    cooler cload pairix \
	-p ${task.cpus} \
	${chr_sizes}:1000 \
	${pairs} \
	${id}.cool
    """
}

process cooler_zoomify {
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
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

process juicer {
    cpus 48
    memory '100 GB'
    container 'mblanche/juicer'

    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/hicFiles",
    	mode: 'copy'

    input:
    path pairs from pairs_ch_juicer
    path idx from pairs_px2_juicer
    path chr_sizes from chrSizes_juicer
    
    output:
    path "*.hic" into juicer_out_ch
    
    script:
    id = pairs.name.toString().tokenize('.').get(0)
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

workflow.onComplete {
    
    def mainDir = file('$PWD/work')
    println "done!"
    if (params.cleanup) {
	println "Deleting directory " + mainDir.name
	mainDir.deleteDir()
    }
}

