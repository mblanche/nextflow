#!/usr/bin/env nextflow

params.fastaDir = false
params.outDir = false
params.genome = 'hg38'

params.help   = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        mininmap  --fastaDir ./path/to/fasta --outDir path/to/results

    Mandatory arguments:
        --fastaDir [path]           Path to a directory containing the fasta files to align (can be local or valid S3 location)..
        --outDir [path]             Path to a diectory to save the bigwig coveage files (can be local or valid S3 location).

    Facultative arguments
    
        --genome [path or def]      Id of existing genome or path to a fasta file. Default: hg38

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.fastaDir && params.outDir) ) {
        exit 1, "--fastaDir and --outDir are required arguments. Use --help to get the full usage." 
}



if (!params.genome =~ /hg19|hg38|mm10|dm3|hg38_hla|susScr11|rn6|GRCg6a/){
    Channel
	.fromPath("${params.genome}")
	.set { mm_ref }
    
} else {
    Channel
	.fromPath("${HOME}/ebs/genome/nextflow/${params.genome}/${params.genome}.fa/", checkIfExists: true)
	.set { mm_ref }
}



process minimap {
    echo true
    label 'cpu'
    tag "_${id}"
    container 'mblanche/minimap2'
    
    publishDir "${params.outDir}",
    mode: 'copy'
    
    input:
    path(fasta) from Channel
	.fromPath("${params.fastaDir}/*.{fasta,fa,fasta.gz,fa.gz}")
    
    path(ref) from mm_ref.first()
    
    output:
    tuple val(id), env(chrN), path("*.paf") into paf_ch
    
    script:
    id = (fasta.name.toString().split(/\./))[0]
    """
    chrN=\$(grep -c '^>' ${ref})
    minimap2 -x asm5 -t ${task.cpus} ${ref} ${fasta}  > ${id}.paf
    """
}



process plotly {
    label 'median'
    tag "_${id}"
    container 'mblanche/minimap2'
    
    publishDir "${params.outDir}/plots",
    mode: 'copy'
    
    input:
    tuple id, val(chrN), path(paf) from paf_ch
    
    output:
    path "*"
    
    script:
    """
    pafCoordsDotPlotly.R \
	-i ${paf} \
	-o ${id}.plot \
	-m 2000 \
	-q 500000 \
	-k ${chrN} \
	-s -t -l -p 12
    """
}
