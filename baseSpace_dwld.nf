#!/usr/bin/env nextflow

/*
* TODO: Update this for the current script
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run baseSapce_dwld.nf  --design design.csv --genome hg38 -profile docker
    
    Mandatory arguments:
    
    NextFlow config:  
    -profile                      Configuration profile to use. Can use multiple (comma separated)
                                  Available: standard, awsbatch. Default: standard
    Generic
      --genome                      Name of genomes reference. Available: hg38, mm10
    
    References                      If not specified in the configuration file or you wish to overwrite any of the references
    
    Alignments
      --keepMito                    Reads mapping to mitochondrial contig are not filtered from alignments
      --keepDups                    Duplicate reads are not filtered from alignments
      --keepMultiMap                Reads mapping to multiple locations are not filtered from alignments
      --skipMergeReplicates         Do not perform alignment merging and downstream analysis by merging replicates i.e. only do this by merging resequenced libraries
      --saveAlignedIntermediates    Save the intermediate BAM files from the alignment step - not done by default
    
    Other
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

params.destDir = 'capture'
params.expName = 'tmp'
params.genome = 'hg38'

params.bpProject = 'Capture'
params.biosample = 'PC_I14_29-51'

Channel
    .from(params.biosample)
    .splitCsv()
    .flatten()
    .set { biosample_ch }

Channel
    .fromPath("/mnt/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
    .ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    .view()
    .set { bwa_index }


process findFile {
    container 'mblanche/basespace-cli'

    input:
    val biosample from biosample_ch
    
    output:
    stdout out_ch
    val biosample into bs_ch
    
    """
    bs list dataset \
	--like-type=illumina.fastq.v1.8 \
	--input-biosample=${biosample} \
	--project-name=${params.bpProject} \
	--like-type=illumina.fastq.v1.8 \
	-f csv
    """
}

out_ch
    .splitCsv(header: true)
    .set { out_ch }

process download {
    cpus 4
    container 'mblanche/basespace-cli'
    
    publishDir "/mnt/ebs/ref_push/${params.destDir}/${params.expName}/fastqs", mode: 'copy'
    
    input:
    val bp_id from out_ch.flatMap { it["Id"] }

    output:
    val bp_id into id
    file '*_R1*.fastq.gz' into R1s_ch
    file '*_R2*.fastq.gz' into R2s_ch
    
    """
    bs download dataset \
     	--id=${bp_id} \
 	-o .
    """
}

R1s_ch
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
    }
    .groupTuple()
    .set{ R1_groups_ch }

R2s_ch
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
    }
    .groupTuple()
    .set{ R2_groups_ch }

process bwa_mem {
    cpus 48
    memory '100 GB'
    
    container 'mblanche/bwa-samtools'
    
    publishDir "/mnt/ebs/ref_push/${params.destDir}/${params.expName}/bam", mode: 'copy'
    
    input:
    set id, file(R1_files) from R1_groups_ch
    set id2, file(R2_files) from R2_groups_ch
    file index from bwa_index.first()
    
    output:
    path "*.bam" into sam_ch
    
    script:
    """
    bwa mem -t \$(nproc) \
	-5SP ${index}/${params.genome} \
	<(zcat $R1_files|head -n 100000) \
	<(zcat $R2_files|head -n 100000) \
	| samtools  view -@ \$(nproc) -Shb - \
	| samtools sort -m 2G  -@ \$(nproc) -o ${id}.bam - 
    """
}

process index_bam {
    cpus 48
    memory '100 GB'
    
    container 'mblanche/bwa-samtools'
    
    input:
    path bam from sam_ch
    
    output:
    path bam into sam_ch2

    """
    samtools index -@ \$(nproc) ${bam}
    """
}

sam_ch2.view()

