#!/usr/bin/env nextflow

params.design = false
params.outDir = false
params.help   = false
params.k = false
params.K = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    hlaAssembly.nf  --design --outDir path/to/results.pdf


    Mandatory arguments:
        --design [file]           Path to a csv file with specie,path/to/ccs. Multiple ccs can be used by assigning same specie to different file. 
        --outDir [path]           Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.design&& params.outDir)) {
        exit 1, "--bamDir and --outDi are required arguments. Use --help to get the full usage." 
}

if (params.k) {
    k = params.k.toString().split(/,/,-1)
}  else {
    k = 501
}

if (params.K) {
    K = params.K.toString().split(/,/,-1)
}  else {
    K = 5001
}

Channel
    .from(k)
    .set{k_ch}

Channel
    .from(K)
    .set{K_ch}

Channel
    .fromPath(params.design)
    .splitCsv()
    .map { id, ccs ->
	return(tuple (id,ccs,'--reads ' + file(ccs).name))
    }
    .groupTuple()
    .set{CCS_ch}


process jla {
    label 'cpu'
    tag '_${id}'
    container 'dovetailg/lja:alpha'

    errorStrategy { task.exitStatus == 1 ? 'ignore' : 'retry' }
    
    publishDir "${params.outDir}",
	mode: "copy"

    input:
    tuple val(k), val(K), val(sp), path(ccs), val(cmd) from k_ch
	.combine(K_ch)
	.combine(CCS_ch)
    
    output:
    path("${id}")
    tuple id, path("${id}/assembly.fasta") into fasta_busco_ch, fasta_n50_ch
        
    script:
    reads = cmd[0]
    id = "${sp}_k${k}_K${K}"
    """
    lja -k ${k} -K ${K} -t ${task.cpus} -o ${id} ${reads}
    """
}


process busco {
    label 'cpu'
    tag '_${id}'
    container 'ezlabgva/busco:v5.3.1_cv1'
    
    input:
    tuple id, path(fasta),val(recNum) from fasta_busco_ch
	.map{id,fasta ->
	    recNum = fasta.countFasta()
	    return(tuple(id,fasta,recNum))
	}

    output:
    tuple id, path("${id}/short_summary.specific.*.txt") into busco_report_ch
    
    when:
    recNum > 0
    
    script:
    """
    busco -m genome -i ${fasta} -f -o ${id} -l eukaryota_odb10 --cpu ${task.cpus} --quiet
    """
}


process busco_reports {
    label 'mezzo'
    tag '_${id}'
    container 'ubuntu:20.04'
    
    publishDir "${params.outDir}/busco",
	mode: "copy"

    input:
    tuple id, path(busco) from busco_report_ch.collect()
    
    output:
    path("*.csv")
    
    script:
    """
    busco_report.sh ${busco} > busco_report.csv
    """
}


process n50 {
    label 'mezzo'
    tag '_${id}'
    container 'mblanche/n50'
    
    publishDir "${params.outDir}/n50",
	mode: "copy"
    
    input:
    tuple id, path(fasta), val(recNum) from fasta_n50_ch
    	.map{id,fasta ->
	    recNum = fasta.countFasta()
	    return(tuple(id,fasta,recNum))
	}

    output:
    path("*.json")

    when:
    recNum > 0
    
    script:
    """
    n50.py -a ${fasta} -o ${id}_n50.json
    """
}
