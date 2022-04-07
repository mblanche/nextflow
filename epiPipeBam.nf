#!/usr/bin/env nextflow


params.bamDir = false
params.outDir = false
params.help   = false

params.genome = 'hg38'

params.resolutions = false
params.ABresolutions = false

params.capture = false
params.hichip = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        bam2bw  --bamDir ~/path/to/bam/location

         or

        bam2bw  --bamDir ~/path/to/bam/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --bamDir [path]              Name of the a direcoty with bam files and their index to compute coverages (can be local or valid S3 location.
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

    Alignment:
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.

    Arrowhead parameters:
        --resolutions [integers]     Comma-seperated list of resolutions in kb for loops finding. Default: [5,10]
        --ABresolutions [integers]   Comma-seperated list of resolutions in kb. Default: [32,64,128]

    Sub-Steps:
        --capture                    Runs only steps required for HiC capture. Will not permoform AB, TAD and loop calls.
        --hichip                     Runs only steps required for HiChIP. Will not permoform AB, TAD and loop calls.
    """.stripIndent()
}

if (!params.bamDir) {
        exit 1, "--bamDir is a required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.bamDir).getParent()
} else {
    outDir = params.outDir
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

if (!params.genome =~ /hg19|hg38|mm10|dm3/){
    exit 1, "Only hg38, mm10 and dm3 genomes are currently offered"
} else {
    Channel
	.fromFilePairs("${HOME}/ebs/genome/nextflow/${params.genome}/*.{amb,sa,pac,ann,bwt,fa}", size: -1, checkIfExists: true)
	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
	.set { bwa_index }
    Channel
	.fromPath("${HOME}/ebs/genome/nextflow/${params.genome}/${params.genome}.fa", checkIfExists: true)
	.ifEmpty { exit 1, "Genome not found: ${params.genome}" }
	.set { abcomp_genome_ch }
    
}

if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromPath("${params.bamDir}/*.bam")
    .set{bam_ch}


process bamFilter {
    label 'bamFilter'
    tag "_${id}"
    cpus 8
    memory '32 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    path(bam) from bam_ch.first()

    output:
    tuple id, path("${id}_filtered.bam") into bamFilt_ch
    //tuple id, path("${id}_filtered.bam"), path("${id}_filtered.bam.bai") into bam2bwIn_ch
    
    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    samtools view -h -f1 -F 3328 ${bam} | samtools sort -n -O BAM -@ ${task.cpus} -o ${id}_filtered.bam -
	##samtools index ${id}_filtered.bam
    """

}


/*
process bam2bw {
    tag "_${id}"
    cpus 20
    memory '150 GB'
    container 'mblanche/r-cov'
    
    publishDir "${outDir}/bigwigs",
    	mode: 'copy'
    
    input:
    tuple id, path(bam), path(bai) from bam2bwIn_ch
        
    output:
    tuple id, path("*.bw") into bigwig_out_ch

    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}
 */


process bam2pairs {
    label 'bam2pairs'
    tag "_${id}"
    cpus 16
    memory '64 GB'
    container 'dovetailg/pairtools'
    
    publishDir "${outDir}/validPairs",
    	mode: 'copy'
    
    input:
    tuple id, path(bam) from bamFilt_ch
    
    output:
    tuple id, path("*.gz"), path("*.px2") into pairs_chrSize_ch

    script:
    """
    bam2pairs ${bam} ${id}
    #mv ${id}.bsorted.pairs.gz ${id}.valid.pairs.gz
    #mv ${id}.bsorted.pairs.gz.px2 ${id}.valid.pairs.gz.px2
    """
}





process chr_size {
    tag "_${id}"
    cpus 8
    memory '32 GB'
    container 'dovetailg/pairtools'
    
    input:
    tuple id, path(pairs), path(idx) from pairs_chrSize_ch
    
    output:
    tuple id, path(pairs), path(idx), path("*.tsv") into pairs_ch_cooler, pairs_ch_juicer
    
    script:
    """
    pairix -H -f ${pairs} \
	| awk -v OFS='\t' '/^#chromsize/  {print \$2,\$3}' \
	| sort -V -k1,1 \
	> chr_size.tsv
    """
}

process cooler_sort {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    input:
    tuple id, path(pairs), path(idx), path(chr_sizes) from pairs_ch_cooler
    
    output:
    tuple id, path("*.pairs.sorted.txt.gz"), path("*pairs.sorted.txt.gz.px2"), path(chr_sizes) into sort_cload_ch
    
    script:
    """
    cooler csort \
	-c1 2 -p1 3 -c2 4 -p2 5 \
	-i pairix \
	-p ${task.cpus} \
	--out ${id}.pairs.sorted.txt.gz \
	$pairs \
	$chr_sizes
    """
}

process cooler_cload {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    input:
    tuple id, path(pairs), path(idx), path(chr_sizes) from sort_cload_ch
    
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
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${outDir}/coolerFiles",
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
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${outDir}/coolerFiles",
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

process juicer {
    tag "_${id}"
    cpus 24
    memory '150 GB'
    container 'mblanche/juicer'
    
    publishDir "${outDir}/hicFiles",
    mode: 'copy'
    
    input:
    tuple id, path(pairs), path(idx), path(chr_sizes) from pairs_ch_juicer
    
    output:
    tuple id, path("*.hic") into arrowhead_ch, hiccups_ch

    script:
    """
    java -Xmx96000m -Djava.awt.headless=true \
	-jar /juicer_tools.jar pre \
	--threads ${task.cpus} \
	-j ${task.cpus} \
	-k VC,VC_SQRT,KR,SCALE \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    """
}

if (!(params.capture|params.hichip)){
    process arrowhead {
	tag "_${id}"
	cpus 12
	memory '40 GB'
	container "mblanche/juicer"
	
	publishDir "${outDir}/arrowHead",
	    mode: 'copy'
	
	input:
	tuple id, path(hic), val(res) from arrowhead_ch
	    .combine(Channel.from(resolutions))
	
	output:
	tuple id, path("${id}_${res}kb") into arrowhead_out_ch
	
	script:
	bpRes = res.toInteger() * 1000
	"""
	mkdir -p ${id}_${res}kb && touch ${id}_${res}kb/${bpRes}_blocks.bedpe
	java -Xmx24000m \
	    -jar /juicer_tools.jar \
	    arrowhead \
	    --threads ${task.cpus} \
	    --ignore-sparsity \
	    -r ${bpRes} \
	    -k KR \
	    ${hic} \
	    ${id}_${res}kb
	"""
    }
    
    process hiccups {
	tag "_${id}"
	label 'gpu'
	accelerator 1
	cpus 6
	memory '32 GB'
	container "mblanche/hiccups-gpu"
	
	publishDir "${outDir}/hiccups/",
	    mode: 'copy'
	
	input:
	tuple id, path(hic), val(res)  from hiccups_ch
            .combine(Channel.from(resolutions.collect{it*1000}.join(',')))
	
	output:
	tuple id, path("${id}_loops") into hiccups_out_ch
	
	script:
	"""
	java -Xmx24000m \
	    -jar /juicer_tools.jar \
	    hiccups \
	    --threads ${task.cpus} \
	    --ignore-sparsity \
	    -m 500 \
	    -r ${res} \
	    -k KR \
	    ${hic} \
	    ${id}_loops
	"""
    }
    
    process mustache {
	tag "_${id}"
	cpus 24
	memory '48 GB'
	container "mblanche/mustache"
	
	publishDir "${outDir}/mustache",
	    mode: 'copy'
	
	input:
	tuple id, path(mcool), val(res)  from mustache_mcool_ch
	    .combine(Channel.from(1000,4000,16000))
	
	output:
	tuple id, path("*.tsv") into mustache_2_merge_ch
	
	script:
	"""
	touch ${id}_${res}kb_loops.tsv 
	mustache -p ${task.cpus} \
	    -f ${mcool} \
	    -r ${res} \
	    -o ${id}_${res}kb_loops.tsv
	"""
    }
    
    process ABcomp {
	tag "_${id}"
	cpus 1
	memory '12 GB'
	container "mblanche/fan-c"
	
	publishDir "${outDir}/AB_comp",
	    mode: 'copy'
	
	input:
	tuple id, path(cool), val(resKB) from abcomp_mcool_ch
    	    .combine(Channel.from(ABresolutions))
	
	path(genome) from abcomp_genome_ch.first()
	
	output:
	tuple id, path("*.bed"), path("*.ab") into fanc_out_ch
	
	script:
	res = resKB.toInteger() * 1000
	"""
	fanc compartments \
	    -f \
	    -v ${id}_eigenV_${resKB}kb.bed \
	    -d ${id}_AB_${resKB}kb.bed \
	    -g ${genome} \
	    ${cool}@${res} \
	    ${id}_${resKB}kb.ab
	"""
	
    }
}
*/
