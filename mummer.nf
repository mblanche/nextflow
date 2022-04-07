#!/usr/bin/env nextflow

params.fasta = false
params.outDir = false
params.genome = 'hg38'

params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        mummer  --fasta ./path/to/hap1,./path/to/hap2

    Mandatory arguments:
        --fasta [path]              Comma delimited list of paths to bam files
        --outDir [path]             Path to a diectory to save the bigwig coveage files (can be local or valid S3 location). Default basedir of fasta file

    Facultative arguments
    
        --genome [path or def]      Id of existing genome or path to a fasta file. Default: hg38

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.fasta && params.outDir) ) {
        exit 1, "--fasta and --outDir are required arguments. Use --help to get the full usage." 
}

if (!params.genome =~ /hg19|hg38|mm10|dm3|hg38_hla|susScr11|rn6|GRCg6a/){
    Channel
	.fromPath("${params.genome}",checkIfExists: true)
	.ifEmpty { exit 1, "Genome file ${params.genome} does not exist." }
	.set { mm_ref }
    
} else {
    Channel
	.fromPath("${HOME}/ebs/genome/nextflow/${params.genome}/${params.genome}.fa/", checkIfExists: true)
	.set { mm_ref }
}


process mummer {
    label 'cpu'
    tag "_${id}"
    container 'mblanche/mummer4'
    
    publishDir "${params.outDir}",
    mode: 'copy'
    
    input:
    path(fasta) from Channel
	.from(params.fasta)
	.splitCsv()
	.flatten()

    path(ref) from mm_ref.first()
    
    output:
    tuple val(id), path("*.mums") into paf_ch
 
    script:
    id = (fasta.name.toString().split(/\./))[0]
    """
    nucmer --maxmatch -l 80 -c 100 -t ${task.cpus} -p ${id} ${ref} ${fasta} 
    """
}

/*
process plotly {
    label 'median'
    tag "_${id}"
    container 'mblanche/minimap2'
    
    publishDir "${params.outDir}/plots",
    mode: 'copy'
    
    input:
    tuple id, path(paf) from paf_ch
    
    output:
    path "*"
    
    script:
    """
    pafCoordsDotPlotly.R \
	-i ${paf} \
	-o ${id}.plot \
	-m 2000 \
	-q 500000 \
	-k 10 \
	-s -t -l -p 12
    """
}
*/
