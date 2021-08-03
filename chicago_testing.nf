#!/usr/bin/env nextflow

params.bamDir = false
params.baits  = false
params.outDir = false
params.help   = false

params.resolutions = false

params.genome = 'hg38'

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        chicago.nf  --bamDir ~/path/to/bam/location --baits ~/path/to/baits.bed

         or

        chicago.nf  --bamDir ~/path/to/bam/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --bamDir [path]              Name of the a direcoty with bam files and their index to compute coverages. Can be local or valid S3 location.
        --baits [path]               Path to the bed files containing the baits. Can be local or valid S3 location.
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save the bigwig coveage files. Can be local or valid S3 location.
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.

    Chicago parameters:
        --resolutions [integers]     Comma-seperated list of resolutions for computing genomic bins. Default: [5,10,20]

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.bamDir && params.baits)) {
    exit 1, "--bamDir and --baits are required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
    outDir = file(params.bamDir).getParent() + "/chicago"
} else {
    outDir = params.outDir
}

if (!params.genome =~ /hg19|hg38|mm10|dm3/){
    exit 1, "Only hg38, mm10 and dm3 genomes are currently offered"
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10,20]
}



process cleanUpBam {
    label 'index'
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'

    input:
    path(bam) from Channel
	.fromPath("${params.bamDir}/*.bam")
    
    output:
    tuple id, path("*-cleanedUp.bam") into cleanBam_ch

    script:
    id = bam.name.toString().take(bam.name.toString().lastIndexOf('.'))
    """
    samtools index -@${task.cpus} ${bam} \
	&& samtools view -@ ${task.cpus} -Shu -F 2048 ${bam} \
	| samtools sort -n -@ ${task.cpus}  -o ${id}-cleanedUp.bam -
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
    tuple path(baits), val(genome), val(res) from Channel.fromPath(params.baits)
	.combine(Channel.from(params.genome).first())
	.combine(Channel.from(resolutions).map{ it * 1000 })
    
    output:
    tuple val(res), path("*.rmap"), path("*.baitmap") into mapFiles_ch
    
    script:
    """
    prep4Chicago ${baits} ${res} ${genome}
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
    memory '16 GB'
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
