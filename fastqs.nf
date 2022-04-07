#!/usr/bin/env nextflow

params.outDir = false

params.biosample = false
params.genewizMap = false
params.fastqDir = false

params.libraryID = false

params.genome = 'hg38'

params.resolutions = false
params.ABresolutions = false

params.HiChIP = false

params.help = false

params.mapQ = 40
params.capture = false
params.hichip = false
    
params.copyFastq = false
params.config = 'default'

def helpMessage() {
    log.info"""
    Usage:
    
    The typical command for running the pipeline is as follows:

        fastqc.nf --fastqDir ~/ebs/ref_push/prodEpi/myExpName/fastqs/ --outDir ~/ebs/ref_push/prodEpi/myExpName

    Mandatory arguments:
        --outDir [pathOfDir]         Path of where to publish the data, can be local or an S3 accessible bucket

    Input (minimum one of these, can be mixed):
        --fastqDir [str]             Path to a directoy of fastq files. Can be S3 path.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !params.fastqDir ){
    exit 1, "You need to use --fastqDir path/to/bam/files. Use --help to get the full usage." 
}

if ( !(params.outDir) ){
    exit 1, "--outDir is a required arguments. Use --help to get the full usage." 
}
	      
	      


process FASTQC {
    cpus 8
    container 'mblanche/fastqc'

    publishDir "${params.outDir}/fastqc",
	mode: 'copy',
        saveAs: { filename ->
        filename.endsWith('.zip') ? "zips/$filename" : filename
    }

    input:
    path(reads) from Channel
	.fromPath("${params.fastqDir}/*.fastq.gz",
		  checkIfExists: true)
    
    output:
    path '*.{zip,html}' into ch_fastqc_reports_mqc
    
    script:
    """
    fastqc -q -t $task.cpus ${reads}
    """

}

process MULTIQC {
    container 'mblanche/fastqc'
    
    publishDir "${params.outDir}/multiqc",
	mode: 'copy'
    
    input:
    path ('fastqc/*') from ch_fastqc_reports_mqc.collect().ifEmpty([])

    output:
    path '*multiqc_report.html' into ch_multiqc_report
    path '*_data'

    script:
    """
    multiqc . 
    """
}
