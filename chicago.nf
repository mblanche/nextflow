#!/usr/bin/env nextflow

params.bamDir = false
params.outDir = false
params.help   = false
params.libraryID = false

params.resolutions = "5,10,20"

params.panel = 'hs_pc_1'

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        chicago.nf  --bamDir ~/path/to/bam/location
    
    or
    
        chicago.nf  --bamDir ~/path/to/bam/location --panel mm_pc_1

    Mandatory arguments:
        --bamDir [path]              Name of the a direcoty with bam files and their index to compute coverages. Can be local or valid S3 location.
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save the chicago result files. Can be local or valid S3 location. Default: Parent directory of the bam directory
        --panel [string]             Id of the panel to use. Valid id: hs_pc_1, mm_pc_1. Default: hs_pc_1.

    Restriting to a set of library:
        --libraryID [str]            Comma seperated list of library prefixes. 

    Chicago parameters:
        --resolutions [integers]     Comma-seperated list of resolutions for computing genomic bins. Default: [5,10,20]

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.bamDir) {
    exit 1, "--bamDir is a required argument. Use --help to get the full usage." 
}

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


if(!params.outDir){
    outDir = file(params.bamDir).getParent() + "/chicago"
} else {
    outDir = params.outDir
}

if (params.libraryID){
    Channel
	.fromPath("${params.bamDir}/*.bam")
	.map{ bam ->
	    for ( item in params.libraryID.split(/,/, -1) ){
		if ( bam.name =~/${item}/ ){
		    return(bam)
		}
	    }
	}
	.set{bam_ch}
} else {
    Channel
	.fromPath("${params.bamDir}/*.bam")
	.set{bam_ch}
}

Channel
    .from(params.resolutions)
    .splitCsv(header: false)
    .flatten()
    .map{it.toInteger() * 1000}
    .set{res_ch}


process index_bam {
    label 'index'
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'

    input:
    path(bam) from bam_ch
    
    output:
    tuple id, path(bam), path("*.bai") into bam_capStats_ch, bam_cleanBam_ch, bam_mapFile_ch

    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    samtools index -@${task.cpus} ${bam}
    """
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
    container 'mblanche/chicago'

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
    container 'mblanche/chicago'
    
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
    container 'mblanche/chicago'

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
	--design-dir .  \
	--cutoff 5 \
	--export-format interBed,washU_text,seqMonk,washU_track \
	${id}_${res}_chinput/${id}_${res}_chinput.chinput \
	${id}_${res}bp
    """
}
