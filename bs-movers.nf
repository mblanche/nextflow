#!/usr/bin/env nextflow

params.outDir = false

params.biosample = false

params.help = false
params.config = 'default'

def helpMessage() {
    log.info"""
    Usage:
    
    The typical command for running the pipeline is as follows:

        ba-movers.nf --biosample bisample1,biosample2,biosample3  --outDir ~/ebs/ref_push/prodEpi/myExpName

    Mandatory arguments:
        --outDir [pathOfDir]         Path of where to publish the data, can be local or an S3 accessible bucket
        --biosample [str]            Comma seperated list of BaseSpace biosamples.

    basespce token:
        --config [config]            if using a basespce config file different from ~/.basename/default.cfg pass ~/.basespace/<config>.cfg

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.biosample) ){
    exit 1, "You need to use either --biosample aBaseSpace-biosample or --fastqDir path/to/bam/files or --genewizMap mappingFiles.csv. Use --help to get the full usage." 
}

if ( !(params.outDir) ){
    exit 1, "--outDir is a required arguments. Use --help to get the full usage." 
}

Channel
    .from(params.biosample)
    .splitCsv()
    .flatten()
    .set { biosample_ch }

/////
// Get BS tokent from config file
////
/// TODO: Add some logic if things are not working
bsConfig = file("${HOME}/.basespace/${params.config}.cfg", checkIfExists:true)

for ( line : bsConfig.readLines() ) {
    if (m = line =~ /apiServer.+=(.+)/) {
        host = m[0][1].replaceAll(' ','')
    } else if ( m = line =~ /accessToken.+=(.+)/){
        token = m [0][1].replaceAll(' ','')
    }
}

    
process get_bs_files {
    cpus 1
    memory '1G'
    container 'mblanche/basespace-cli'
    
    input:
    val bs from biosample_ch
    
    output:
    stdout into bs_id_ch
    
    script:
    """
    findNewestBS.sh ${bs} ${token} ${host}
    """
}


process download_bs {
    echo true
    label "movers"
    cpus 4
    memory '4G'
    container 'mblanche/basespace-cli'
    queue 'moversQ'
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple id, val(bs) from bs_id_ch
	.splitCsv(header: false)
    
    output:
    tuple bs, file("*.fastq.gz") into fastqs_ch
    
    script:
    """
    bs file download --api-server ${host} --access-token ${token} -i ${id} -o . 
    """
}
