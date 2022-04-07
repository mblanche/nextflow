#!/usr/bin/env nextflow

params.outDir = false
params.bamDir = false

params.resolutions = false
params.ABresolutions = false
params.mapQ = 0

params.help = false

def helpMessage() {
    log.info"""
    Usage:
    
    The typical command for running the pipeline is as follows:

        bam2cool --bamDir ~/path/to/bamDir --outDir ~/ebs/ref_push/prodEpi/myExpName

    Mandatory arguments:
        --bamDir [path]              Path of a directory containing the bam file(s) to be cooled, can be local or an S3 accessible bucket.
        --outDir [pathOfDir]         Path of where to publish the data, can be local or an S3 accessible bucket.
    
    Additional parameters:
        --mapQ [int]                 Quqlity score to filter reads. Integer between 0 and 60. Default: 0.
        --resolutions [integers]     Arrowhead and hiccup resolutions. Comma-seperated list of resolutions in kb for loops finding. Default: [5,10]
        --ABresolutions [integers]   ABcomp resolutions. Comma-seperated list of resolutions in kb. Default: [32,64,128]
    
""".stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.bamDir && params.outDir) ){
            exit 1, "--bamDir and --outDir are required arguments. Use --help to get the full usage." 
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10]
}

if (params.ABresolutions){
    ABresolutions = params.ABresolutions.toString().split(/,/,-1)
} else {
    ABresolutions = [32,64,128]
}

process chr_size {
    tag "_${id}"
    label "mezzo"    
    container 'mblanche/bwa-samtools'

    publishDir "${params.outDir}/chr_size_res",
	mode: 'copy'
	
    input:
    path(bam) from Channel.fromPath("${params.bamDir}/*.bam")
    
    output:
    tuple id, path(bam), path("*.tsv") into pairtools_parse_ch, cooler_chrsize_ch
    
    script:
    id = (bam.name.toString().split(/\./))[0]
    """
    samtools view -H ${bam} | \
	awk -v OFS='\t' '/^@SQ/ && !(\$2 ~ /:(chr|"")M/) {split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' | \
	sort -V -k1,1 \
	> chr_size.tsv
    """
}

process pairtools_parse {
    tag "_${id}"
    label "cpu"
    container 'mblanche/pairtools'
    
    input:
    tuple id, path(bam), path(chr_sizes) from pairtools_parse_ch

    output:
    tuple id, path("*.pairsam.gz") into pairsam_ch

    script:
    """
    pairtools parse \
	--min-mapq ${params.mapQ} \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--chroms-path ${chr_sizes} \
	--output ${id}.pairsam.gz \
	${bam} 
    """
}


process pairtools_sort {
    tag "_${id}"
    label "cpu"
    container 'mblanche/pairtools'
    
    input:
    tuple id, path(pairsam) from pairsam_ch

    output:
    tuple id, path("*_sorted.pairsam.gz") into sorted_ps_ch
    
    script:
    """
    pairtools sort \
	--nproc ${task.cpus} \
	--output ${id}_sorted.pairsam.gz \
	$pairsam
    """
}


process pairtools_dedup {
    echo true
    tag "_${id}"
    label "cpu"
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/pairtools_stat",
	mode: 'copy',
	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    tuple id, path(pairsam) from sorted_ps_ch

    output:
    tuple id, path("*_dedup.pairsam.gz") into dedup_ps_ch
    path("*.stats")
    
    script:
    """
    echo redoing
    pairtools dedup --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--mark-dups \
	--output-stats ${id}_pairtools.stats  \
	--output ${id}_dedup.pairsam.gz \
	${pairsam}
    """
}


process pairtools_split_dedup {
    tag "_${id}"
    label "cpu"
    container 'mblanche/pairtools'

    publishDir "${params.outDir}/valid_pairs",
	mode: 'copy'

    input:
    tuple id, path(pairsam) from dedup_ps_ch
    
    output:
    tuple id, path("*.valid.pairs.gz") into pairs_ch
    
    script:
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-pairs ${id}_PT.valid.pairs.gz  \
	${pairsam}
    """
}

process index_pairs {
    tag "_${id}"
    label "mezzo"
    container 'mblanche/pairtools'

    publishDir "${params.outDir}/valid_pairs",
	mode: 'copy'

    input:
    tuple id, path(pairs) from pairs_ch
    
    output:
    tuple id, path(pairs), path("*.px2")  into pairs_ch_cooler
    
    script:
    """
    pairix ${pairs}
    """
}


/*
process cooler_cload {
    tag "_${id}"
    label "cpu"
    container 'mblanche/cooler'

    input:
    tuple id, path(pairs), path(idx) from pairs_ch_cooler

    path(chr_sizes) from cooler_chrsize_ch
    
    output:
    tuple id, path("*.cool") into balance_cooler_ch
        
    script:
    """
    cooler cload pairix \
	-p ${task.cpus} \
	${chr_sizes}:1000 \
	${pairs} \
	${id}.cool
    """
}

/*
process balance_cooler {
    tag "_${id}"
    label "large"
    container 'mblanche/cooler'
    
    publishDir "${params.outDir}/coolerFiles",
    	mode: 'copy'
    
    input:
    tuple id, path(cooler) from balance_cooler_ch

    output:
    tuple id, path(cooler) into zoomify_cooler_ch 
    
    script:
    """
    cooler balance --force -p ${task.cpus} ${cooler}
    """
}

process cooler_zoomify {
    tag "_${id}"
    label "cpu"
    container 'mblanche/cooler'
    
    publishDir "${params.outDir}/coolerFiles",
    	mode: 'copy'
    
    input:
    tuple id, path(cooler) from zoomify_cooler_ch

    output:
    tuple id, path("*.mcool") into mustache_mcool_ch, abcomp_mcool_ch
    
    script:
    """
    cooler zoomify --balance -p ${task.cpus} ${cooler}
    """
}
*/
