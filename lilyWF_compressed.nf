#!/usr/bin/env nextflow
/*
* TODO: Update this for the current script
*/


params.expDir = 'prodEpi'
params.expName = 'tmp3'
params.genome = 'hg38'

params.noPairTools = false
params.noJuicer = false
params.noCooler = false
params.noCoverage = false


Channel
    .fromFilePairs("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/fastqs/*_R{1,2}*.fastq.gz",flat: true)
    .groupTuple()
    .set{fastqs_ch}

Channel
    .fromPath("${HOME}/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
    .ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    .set { bwa_index }



process bwa_mem {
    tag "bwa_mem_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    set id, file(R1s), file(R2s) from fastqs_ch
    file index from bwa_index.first()
    
    output:
    path "${id}.bam" into sam4chrSize, sam4pt
    
    script:
    """
    bwa mem -5SP -t ${task.cpus} \
	${index}/${params.genome} \
	<(zcat ${R1s} | head -n 40000) \
	<(zcat ${R2s} | head -n 40000) |\
    samtools view -@ ${task.cpus} -Shb - \
	>${id}.bam
    """
}

process make_chr_size {
    tag "chr_size_${id}"
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
    tag "pt_parse_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'

    publishDir "result", mode: "copy"
    
    input:
    path sam from sam4pt
    path chr_sizes from chrSizes_pt
    
    output:
    path "*.pairsam.gz" into pairsam_ch

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
	> ${id}.pairsam.gz
    """
}


process pairtools_sort {
    tag "pt_sort_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    path sam from pairsam_ch

    output:
    path "*_sorted.pairsam.gz" into sorted_ps_ch

    when:
    !params.noPairTools

    script:
    id = sam.name.toString().replaceFirst(/.pairsam.gz/,'')
    """
    pairtools sort --nproc ${task.cpus} $sam \
	>${id}_sorted.pairsam.gz
    """
}

/*
process pairtools_dedup {
    tag "pt_dedup_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    path sam from sorted_ps_ch

    output:
    path "*_dedup.pairsam.gz" into dedup_ps_ch
    path "*_unmapped.pairsam.gz" into unmapped_ps_ch
    path "*_pairtools.stats" into ps_stats

    when:
    !params.noPairTools
    
    script:
    id = sam.name.toString().replaceFirst(/_sorted.pairsam.gz/,'')
    """
    pairtools dedup --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--mark-dups \
	--output-stats ${id}_pairtools.stats  \
	--output ${id}_dedup.pairsam.gz \
	--output-unmapped ${id}_unmapped.pairsam.gz \
	--cmd-in gzip --cmd-out gunzip \
	${sam}
    """
}

process pairtools_split_dedup {
    tag "pt_split_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    path sam from dedup_ps_ch
    
    output:
    path "*.bam" into bam_parts_ch
    path "*.valid.pairs.gz" into pairs_parts_ch

    when:
    !params.noPairTools
    
    script:
    id = sam.name.toString().replaceFirst(/_decup.pairsam.gz/,'')
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_PT.bam  \
	--output-pairs ${id}.valid.pairs.gz  \
	--cmd-in gzip --cmd-out gunzip \
	${sam}
    """
}

process merge_pairs {
    tag "pt_merge_${id}"
    cpus 48
    memory '100 GB'
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
    tag "bam_merge_${id}"
    echo true
    cpus 48
    memory '100 GB'
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
    samtools merge -@ ${task.cpus} ${id}.bam $bam_part 
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
    tag "pt_split2_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/pairtools'
        
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/unmapped",
    	mode: 'copy'

    input:
    path sam from unmapped_ps_ch
    
    output:
    path "*_unmapped.bam" into unmapped_bam_ch
    path "*_unmapped.valid.pairs.gz" into unmapped_pairs_ch

    script:
    id = sam.name.toString().replaceFirst(/_unmapped.pairsam.gz/,'')
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_unmapped.bam  \
	--output-pairs ${id}_unmapped.valid.pairs.gz  \
	--cmd-in gzip --cmd-out gunzip \
	${sam}
    """
}

process cooler_cload {
    tag "cooler_${id}"
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
    tag "zoomify_${id}"
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
    tag "juicer_${id}"
    cpus 48
    memory '100 GB'
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
    java -Xmx96000m -Djava.awt.headless=true \
	-jar /juicer_tools_1.22.01.jar pre \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    
    ## -j ${task.cpus} \
    ## -k VC,VC_SQRT,KR,SCALE \
    """
}


process deeptools_bw {
    tag "dt_cov_${id}"
    echo true
    cpus 48
    memory '100 GB'
    container 'mblanche/deeptools'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/bigwigs",
    	mode: 'copy'
    
    input:
    path bam from bam_bw_ch
    path idx from idx_bw_ch
    
    when:
    !params.noCoverage && !params.noPairTools
    
    script:
    id = bam.name.toString().replaceFirst(/.bam/,'')
    """
    bamCoverage -p ${task.cpus} -bs 1 -b ${bam} -o ${id}.bw
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
*/
