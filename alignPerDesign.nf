#!/usr/bin/env nextflow

params.help   = false

params.design = false
params.outDir = false

params.genome = 'hg38'
params.mapQ = 40

params.resolutions = false
params.ABresolutions = false

params.chicago = false
params.hichip = false

params.panel = 'hs_pc_1'
params.chicagoRes  = "5,10,20"

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        alignPerDesign.nf  --design ~/path/to/bam/location/design.csv --fastqDir ~/path/to/bam/location

         or

        alignPerDesign.nf  --design ~/path/to/bam/location/design.csv --fastqDir ~/path/to/bam/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

        --design [path]              Path to a design files in csv. First col: id, second col: replicate, third col: path to R1, fourth col: path to R2
    
    Alignment:
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.
        --mapQ [int]                 Quqlity score to filter reads. Integer between 0 and 60. Default: 40.

    Additional parameers:
        --resolutions [integers]     Arrowhead and hiccup resolutions. Comma-seperated list of resolutions in kb for loops finding. Default: [5,10]
        --ABresolutions [integers]   ABcomp resolutions. Comma-seperated list of resolutions in kb. Default: [32,64,128]

    Sub-Steps:
        --chicago                    Runs only steps required for HiC capture. Will not permoform AB, TAD and loop calls.
        --hichip                     Runs only steps required for HiChIP. Will not permoform AB, TAD and loop calls.

    
    Chicago parameters:
        --resolutions [integers]     Comma-seperated list of resolutions for computing genomic bins. Default: [5,10,20]
        --panel [string]             Id of the panel to use. Valid id: hs_pc_1, mm_pc_1. Default: hs_pc_1.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.outDir && params.design)) {
    exit 1, "--outDir and --design is are required arguments. Use --help to get the full usage." 
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
	.fromFilePairs("${HOME}/ebs/genome/nextflow/bwa-mem2-index/${params.genome}/*.{0123,amb,ann,bwt.2bit.64,pac,fa}", size: -1, checkIfExists: true)
	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
	.set { bwa_index }
    
    Channel
	.fromPath("${HOME}/ebs/genome/nextflow/${params.genome}/${params.genome}.fa", checkIfExists: true)
	.ifEmpty { exit 1, "Genome not found: ${params.genome}" }
	.set { abcomp_genome_ch }
    
}

// Creating fastq channel from design files
Channel
    .fromPath(params.design)
    .splitCsv(header: false)
    .view()
    .set{ fqs_ch }


// Creating Chicago reolution channel
Channel
    .from(params.chicagoRes.toString())
    .splitCsv(header: false)
    .flatten()
    .map{it.toInteger() * 1000}
    .set{res_ch}

process bwa_mem {
    tag "_${prefix}"
    label 'cpu'
    container 'mblanche/bwa-mem2'
    
    input:
    tuple id, val(rep), path(R1s), path(R2s) from fqs_ch.first()
    tuple index, path(index_files) from bwa_index.first() 
    
    output:
    tuple id, val(rep), val(prefix), path("*.bam"), path("*.tsv") into  pairtools_parse_ch
    
    script:
    prefix = R1s.name.toString().replaceFirst(/.fastq.+/,"")
    """
    bwa-mem2 mem -5SP -t ${task.cpus} \
    	${index} \
    	<(zcat ${R1s}) \
    	<(zcat ${R2s}) \
	|samtools view -@ ${task.cpus} -Shb -o ${prefix}.bam - \
	&& samtools view -H ${prefix}.bam | \
	awk -v OFS='\t' '/^@SQ/ && !(\$2 ~ /:(chr|"")M/) {split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' | \
	sort -V -k1,1 \
	> chr_size.tsv
    """
}


process pairtools_parse {
    tag "_${prefix}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, val(rep), val(prefix), path(sam), path(chr_sizes) from pairtools_parse_ch

    output:
    tuple id, val(rep), val(prefix), path("*.pairsam.gz") into pairsam_part_ch

    script:
    """
    pairtools parse \
	--min-mapq ${params.mapQ} \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--chroms-path ${chr_sizes} \
	--output ${prefix}.pairsam.gz \
	${sam} 
    """
}

process pairtools_merge_lane {
    tag "_${prefix_out}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'
    
    input:
    tuple id, val(rep), val(prefix), path(sam) from pairsam_part_ch
	.groupTuple(by: [1,0])
    
    output:
    tuple id, val(prefix_out), path("${prefix_out}*pairsam.gz") into pairsam_ch

    script:
    prefix_out="${id}-${rep}"
    if (sam.sort().size() >1) {
	"""
	pairtools merge -o ${prefix_out}.pairsam.gz --nproc ${task.cpus} ${sam}
	"""
    } else {
	"""
	ln -sf ${sam} ${prefix_out}_ML.pairsam.gz
	"""
    }
    
}

process pairtools_sort {
    tag "_${prefix}"
    cpus 14
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, val(prefix), path(sam) from pairsam_ch

    output:
    tuple id, val(prefix), path("*_sorted.pairsam.gz") into sorted_ps_ch
    
    script:
    """
    mkdir -p tmp 
    pairtools sort --tmpdir ./tmp  \
	--nproc ${task.cpus} \
	--output ${prefix}_sorted.pairsam.gz \
	$sam 
    """
}

process pairtools_dedup {
    tag "_${prefix}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    tuple id, val(prefix), path(sam) from sorted_ps_ch

    output:
    tuple id, val(prefix), path("*_dedup.pairsam.gz") into dedup_ps_ch
    tuple id, val(prefix), path("*_unmapped.pairsam.gz") into unmapped_ps_ch
    path("*_pairtools.stats") into ps_stats_ch

    script:
    """
    pairtools dedup --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--mark-dups \
	--output-stats ${prefix}_pairtools.stats  \
	--output ${prefix}_dedup.pairsam.gz \
	--output-unmapped ${prefix}_unmapped.pairsam.gz \
	${sam}
    """
}


process pairtools_stats_merge {
    tag "_${id}"
    cpus 1
    memory '4 GB'
    container 'mblanche/pt-stats'

    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    path(stats) from ps_stats_ch
	.collect()
    
    output:
    path('pairtoolsStats.csv') into merged_stats_ch
    
    script:
    """
    pairtoolsStat.sh ${stats} > pairtoolsStats.csv
    """
}

process pairtools_split_dedup {
    tag "_${prefix}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, val(prefix), path(sam) from dedup_ps_ch
    
    output:
    tuple id, val(prefix), path("*.bam") into bam_parts_ch
    tuple id, val(prefix), path("*.valid.pairs.gz") into pairs_parts_ch, pairs_parts_ch_test

    script:
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${prefix}_PT.bam  \
	--output-pairs ${prefix}_PT.valid.pairs.gz  \
	${sam}
    """
}

process merge_bam {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, val(rep), path(bam_part) from bam_parts_ch
	.groupTuple()

    output:
    tuple id, path("*.bam") into merged_bam_sort_ch

    script:
    bam_files = bam_part.sort()
    if (bam_files.size() >1) {
	"""
	samtools merge -@ ${task.cpus} ${id}_MB.bam ${bam_part}
	"""
    } else {
	"""
	ln -s ${bam_part} ${id}_MB.bam
	"""
    }
}

process merge_pairs {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/validPairs",
    	mode: 'copy'
    
    input:
    tuple val(id), val(rep), path(pairs) from pairs_parts_ch
	.groupTuple()
    
    output:
    tuple id, path("*.valid.pairs.gz"), path("*.px2") into pairs_chrSize_ch

    script:
    pair_files = pairs.sort()
    if (pair_files.size() >1) {
	"""
	pairtools merge -o ${id}.valid.pairs.gz --nproc ${task.cpus}  ${pairs}
	pairix ${id}.valid.pairs.gz
	"""
    } else {
	"""
	ln -s ${pairs} ${id}.valid.pairs.gz
	pairix ${id}.valid.pairs.gz
	"""
    }
}

process bam_sort {
    tag "bam_sort_${id}"
    cpus 48
    memory '150 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "${params.outDir}/bam",
	mode: 'copy',
	pattern: "${id}.bam"
        
    input:
    tuple id, path(bam) from merged_bam_sort_ch
    
    output:
    tuple id, path("${id}.bam"),path("${id}.bam.bai") into bam_bigwig_ch,
	bam_capStats_ch,
	bam_cleanBam_ch,
	bam_mapFile_ch

    script:
    """
    samtools sort -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} 

    samtools index -@${task.cpus} ${id}.bam
    """
}

process chr_size {
    tag "_${id}"
    cpus 4
    memory '24 GB'
    container 'mblanche/pairtools'
    
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


process cooler_cload {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    input:
    tuple id, path(pairs), path(idx), path(chr_sizes) from pairs_ch_cooler
    
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

process balance_cooler {
    tag "_${id}"
    cpus 48
    memory '100 GB'
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
    cpus 48
    memory '100 GB'
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

process bam2bw {
    tag "_${id}"
    cpus 20
    memory '175 GB'
    
    container 'mblanche/r-cov'
    
    publishDir "${params.outDir}/bigwigs",
    	mode: 'copy'
    
    input:
    tuple id, path(bam),path(idx) from bam_bigwig_ch
        
    output:
    tuple id, path ("*.bw") into bigwig_out_ch

    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}

process juicer {
    tag "_${id}"
    cpus 24
    memory '150 GB'
    container 'mblanche/juicer'
    
    publishDir "${params.outDir}/hicFiles",
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


if (!(params.chicago|params.hichip)){
    process arrowhead {
	tag "_${id}"
	cpus 12
	memory '40 GB'
	container "mblanche/juicer"
	
	publishDir "${params.outDir}/arrowHead",
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
	memory '30 GB'
	container "mblanche/hiccups-gpu"
	
	publishDir "${params.outDir}/hiccups/",
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
	
	publishDir "${params.outDir}/mustache",
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
	
	publishDir "${params.outDir}/AB_comp",
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

if (params.chicago){
    switch(params.panel) {
	case "hs_pc_1":
	    baits ="/mnt/ebs/genome/nextflow/${params.panel}/baits_v1.0.bed.bgz"
	    probes ="/mnt/ebs/genome/nextflow/${params.panel}/probes_v1.0.bed.bgz"
	    genome = "hg38"
	    break
	    
    case "mm_pc_1":
	    baits ="/mnt/ebs/genome/nextflow/${params.panel}/baits_mm10_v1.0.bed.bgz"
	    probes ="/mnt/ebs/genome/nextflow/${params.panel}/probes_mm10_v1.0.bed.bgz"
	    genome = "mm10"
	    break
	    
	default:
	    exit 1, "Wrong capure panel;. Has to be hs_pc_1 or mm_pc_1."
	    break
    } 

    
    process capture_Stats {
	label 'capStats'
	tag "capStats"
	cpus 25
	memory '48 GB'
	container 'mblanche/r-cap-stats'
	
	publishDir "${outDir}/captureStats", mode: 'copy'
	
	input:
	bam_capStats_ch
	    .multiMap { id, bam, idx ->
		bam: bam
		idx: idx
	    }
	    .set{ result }
	
	path(bams) from result.bam.collect()
	path(idx) from result.idx.collect()
	path(probeFile) from Channel.fromPath(probes)
	
	output:
	tuple path("*.pdf"), path("*.csv") into capStats_ch
	
	script:
	"""
	capStats ${probeFile} ${task.cpus} ${bams}
	"""
    }
     
    process make_mapFiles {
	tag "_${id}"
	cpus 1
	memory '8 GB'
	container 'mblanche/chicago:3.3'

	publishDir "${outDir}",
	    saveAs: {filename -> filename.endsWith('.rmap') ? filename : null},
	    mode: 'copy'
	
	input:
	tuple path(baitFile), val(id), path(bam), path(idx), val(res) from Channel.fromPath(baits)
	    .combine(bam_mapFile_ch.first())
	    .combine(res_ch)
	
	output:
	tuple val(res), path("*.rmap"), path("*.baitmap") into mapFiles_ch
	
	script:
	"""
	prep4Chicago ${baitFile} ${res} ${bam}
	"""
    }

    process cleanUpBam {
	label 'cleanUp'
	tag "_${id}"
	cpus 48
	memory '100 GB'
	container 'mblanche/bwa-samtools'

	input:
	tuple id, path(bam), path(idx) from bam_cleanBam_ch
	
	output:
	tuple id, path("*-cleanedUp.bam") into cleanBam_ch

	script:
	"""
	samtools view -@ ${task.cpus} -Shu -F 2048 ${bam} \
	    | samtools sort -n -@ ${task.cpus}  -o ${id}-cleanedUp.bam -
	    """
    }

    process make_design {
	tag "_${id}"
	cpus 1
	memory '8 GB'
	container 'mblanche/chicago:3.3'
	
	input:
	tuple val(res), path(rmap), path(baitmap) from mapFiles_ch

	output:
	tuple val(res), path(rmap), path(baitmap), path("${res}kDesingFiles*") into design_ch
        
	script:
	"""
	python3 /makeDesignFiles_py3.py \
	    --minFragLen 75 \
	    --maxFragLen 30000 \
	    --maxLBrownEst 1000000 \
	    --binsize 20000 \
	    --rmapfile ${rmap} \
	    --baitmapfile ${baitmap} \
	    --outfilePrefix ${res}kDesingFiles
	"""
    }

    process run_Chicago {
	tag "_${id}"
	cpus 1
	memory '182 GB'
	container 'mblanche/chicago:3.3'

	errorStrategy { task.exitStatus in 1 ? 'ignore' : 'terminate' }

	publishDir "${outDir}/${id}_${res}",
	    mode: 'copy'
	
	input:
	tuple id, path(bam), val(res), path(rmap), path (baitmap), path(designFiles) from cleanBam_ch
	    .combine(design_ch)

	output:
	tuple id, path("${id}_${res}_chinput"), path("${id}_${res}bp") into chicago_ch

	script:
	"""
	bam2chicago.sh ${bam} ${baitmap} ${rmap} ${id}_${res}_chinput
	
	runChicago \
	    --design-dir . \
	    --cutoff 5 \
	    --export-format interBed,washU_text,seqMonk,washU_track \
	    ${id}_${res}_chinput/${id}_${res}_chinput.chinput \
	    ${id}_${res}bp
	"""
    }

    
}
