#!/usr/bin/env nextflow

params.outDir = false

params.biosample = false

params.help = false

def helpMessage() {
    log.info"""
    Usage:
    
    The typical command for running the pipeline is as follows:

        ba-movers.nf --biosample bisample1,biosample2,biosample3  --outDir ~/ebs/ref_push/prodEpi/myExpName

    Mandatory arguments:
        --outDir [pathOfDir]         Path of where to publish the data, can be local or an S3 accessible bucket
        --biosample [str]            Comma seperated list of BaseSpace biosamples.

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
    bs biosample content -n ${bs} -F Id -F FilePath -f csv | \
	awk 'BEGIN{OFS =","} \
	NR == 1 {print "biosample", \$0} \
	NR > 1  {print "${bs}", \$0}'
    """
}

process download_bs {
    label "movers"
    cpus 4
    memory '4G'
    container 'mblanche/basespace-cli'
    queue 'moversQ'
    
    publishDir "${params.outDir}/fastqs",
	mode: 'copy'
    
    input:
    tuple bs, val(oriFname), val(newFname), val(id) from bs_id_ch
	.splitCsv(header: true)
	.map { row -> tuple(row.FilePath,row.biosample, row.Id )}
	.groupTuple()
	.map{if (it[2].size() >1){
		x = []
		fname = it[0]
		id = fname.take(fname.indexOf('_'))
		suff = fname.substring(fname.indexOf('_')+1)
		bs = it[1][0]
		
		for (i in (1..it[2].size())) {
		    x.add([bs,fname,"${id}_dupName${i}_${suff}",it[2][i-1]])
		}
		return(x)
	    } else {
		fname = it[0]
		bs = it[1][0]
		id = it[2]
		return([bs,fname,fname,id])
	    }
	}
	.flatten()
	.collate(4)
    
    output:
    tuple bs, file("*.fastq.gz") into fastqs_ch
    
    script:
    """
    bs file download -i ${id} -o . 
    
    if [ "${oriFname}" != "${newFname}" ];then 
    mv ${oriFname} ${newFname}
    fi
    """
}
    
