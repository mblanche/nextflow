#!/usr/bin/env nextflow
/*
* TODO: Update this for the current script
*/


params.expDir = 'capture'
params.expName = 'tmp'
params.genome = 'hg38'

params.noPairTools = false
params.noJuicer = false
params.noCooler = false
params.noCoverage = false


Channel
    .fromFilePairs("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/fastqs/*_R{1,2}*.fastq.gz",flat: true)
    .map { prefix, file1, file2 -> tuple(prefix.split('_')[0], file1, file2) }
    .groupTuple()
    .set{fastqs_ch}

Channel
    .fromPath("${HOME}/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
    .ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    .set { bwa_index }

process bwa_mem {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
	
    input:
    tuple val(id), file(R1s), file(R2s) from fastqs_ch
    path index from bwa_index.first()
    
    output:
    path "*.sam" into sam4chrSize, sam4pt
    
    script:
    """
    bwa mem -5SP -t ${task.cpus} \
	${index}/${params.genome} \
	<(zcat ${R1s}) \
	<(zcat ${R2s}) \
    > ${id}.sam
    """
}


process make_chr_size {
    tag "_${id}"
    label 'local'
    echo true
    
    container 'mblanche/bwa-samtools'

    input:
    path sam from sam4chrSize

    output:
    path 'chr_size.tsv' into chrSizes_pt, chrSizes_cooler, chrSizes_juicer
    
    when:
    !params.noPairTools

    script:
    id = id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    samtools view -H ${sam} | \
	awk -v OFS='\t' '/^@SQ/{split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' > chr_size.tsv
    """
}

process pairtools_parse {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'

    input:
    path sam from sam4pt
    path chr_sizes from chrSizes_pt
    
    output:
    path "*.pairsam" into pairsam_ch

    when:
    !params.noPairTools

    script:
    id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    pairtools parse \
	--min-mapq 40 \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--chroms-path ${chr_sizes} \
	${sam} \
	> ${id}.pairsam
    """
}

process pairtools_sort {
    tag "_${id}"
    cpus 14
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    path sam from pairsam_ch

    output:
    path "*_sorted.pairsam" into sorted_ps_ch

    when:
    !params.noPairTools

    script:
    id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    pairtools sort --nproc ${task.cpus} $sam >${id}_sorted.pairsam
    """
}

process pairtools_dedup {
    tag "_${id}"
    cpus 14
    //memory '100 GB'
    container 'mblanche/pairtools'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    path sam from sorted_ps_ch

    output:
    path "*_dedup.pairsam" into dedup_ps_ch
    path "*_unmapped.pairsam" into unmapped_ps_ch
    path "*_pairtools.stats" into ps_stats

    when:
    !params.noPairTools
    
    script:
    id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    pairtools dedup --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--mark-dups \
	--output-stats ${id}_pairtools.stats  \
	--output ${id}_dedup.pairsam \
	--output-unmapped ${id}_unmapped.pairsam \
	${sam}
    """
}

process pairtools_split_dedup {
    tag "_${id}"
    cpus 14
    //memory '100 GB'
    container 'mblanche/pairtools'

    input:
    path sam from dedup_ps_ch
    
    output:
    path "*.bam" into bam_parts_ch
    path "*.valid.pairs.gz" into pairs_parts_ch

    when:
    !params.noPairTools
    
    script:
    id = sam.name.toString().take(sam.name.toString().lastIndexOf('.'))
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_PT.bam  \
	--output-pairs ${id}.valid.pairs.gz  \
	${sam}
    """
}

process merge_pairs {
    tag "_${id}"
    cpus 14
    //memory '100 GB'
    container 'mblanche/pairtools'

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
    	mode: 'copy'
    
    input:
    tuple val(key), path(pairs) from pairs_parts_ch
    	.map { file ->
	    def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
	}
	.groupTuple()

    output:
    path "*.valid.pairs.gz" into pairs_ch_cooler, pairs_ch_juicer
    path "*.px2" into pairs_px2_cooler, pairs_px2_juicer

    when:
    !params.noPairTools
    
    script:
    id = key.tokenize('_').get(0)
    """
    pairtools merge -o ${id}.valid.pairs.gz --nproc ${task.cpus}  ${pairs}
    pairix ${id}.valid.pairs.gz
    """
}


process merge_bam {
    tag "_${id}"
    echo true
    cpus 48
    //memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple val(key), path(bam_part) from bam_parts_ch
	.map { file ->
	    def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
	}
	.groupTuple()

    when:
    !params.noPairTools

    output:
    path "*.bam" into merged_bam_sort_ch

    script:
    id = key.tokenize('_').get(0)
    """
    samtools merge -@ ${task.cpus} ${id}_PT.bam ${bam_part}
    """
}

process bam_sort {
    tag "bam_sort_${id}"
    cpus 48
    memory '150 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bam", mode: 'copy'
    
    input:
    path bam from merged_bam_sort_ch
    
    output:
    path "${id}.bam" into bam_bw_ch
    path "${id}.bam.bai" into idx_bw_ch
    
    when:
    !params.noPairTools
    
    script:
    id = bam.name.toString().replaceFirst(/_PT.bam/,'')
    """
    samtools sort -m 2G \
	-@ ${task.cus} \
	-o ${id}.bam \
	${bam}

    samtools index -@${task.cpus} ${id}.bam
    """
}

process pairtools_split_unmapped {
    tag "_${id}"
    cpus 14
    //memory '100 GB'
    container 'mblanche/pairtools'
        
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/unmapped",
    	mode: 'copy'

    input:
    path sam from unmapped_ps_ch
    
    output:
    path "*_unmapped.bam" into unmapped_bam_ch
    path "*_unmapped.valid.pairs.gz" into unmapped_pairs_ch

    when:
    !params.noPairTools

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
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
    	mode: 'copy'
    
    input:
    path pairs from pairs_ch_cooler
    path idx from pairs_px2_cooler
    path chr_sizes from chrSizes_cooler
    
    output:
    path "*.cool" into cooler_ch

    when:
    !params.noCooler && !params.noPairTools

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
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/${id}",
    	mode: 'copy'
    
    input:
    path cooler from cooler_ch

    when:
    !params.noCooler && !params.noPairTools

    output:
    path "*.mcool" into mcool_ch
    
    script:
    id = cooler.name.toString().tokenize('.').get(0)
    """
    cooler zoomify --balance -p ${task.cpus} ${cooler}
    """
}


process juicer {
    tag "_${id}"
    cpus 24
    memory '24 GB'
    container 'mblanche/juicer'

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/hicFiles",
    	mode: 'copy'

    input:
    path pairs from pairs_ch_juicer
    path idx from pairs_px2_juicer
    path chr_sizes from chrSizes_juicer
    
    output:
    path "*.hic" into juicer_out_ch

    when:
    !params.noJuicer && !params.noPairTools
    
    script:
    id = pairs.name.toString().tokenize('.').get(0)
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

process bam2bw {
    tag "_${id}"
    cpus 16
    
    container 'mblanche/r-cov'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bigwigs",
    	mode: 'copy'
    
    input:
    path bam from bam_bw_ch
    path idx from idx_bw_ch
        
    output:
    path "*.bw" into bigwig_ch

    when:
    !params.noCoverage && !params.noPairTools
    
    script:
    id = bam.name.toString().replaceFirst(/.bam/,'')
    """
    bam2bw ${bam} ${id}.bw ${task.cpus} | tee ${id}.bw
    """
}

workflow.onComplete {
    // will probably have to do something else with S3 work dir
    println "Working directory is " + workDir
    if (params.cleanup) {
	println "Deleting directory " + workDir
	workDir.deleteDir()
    }
}
