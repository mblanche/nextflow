#!/usr/bin/env nextflow

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


//.fromPath("/mnt/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
Channel
    .fromPath("/mnt/ebs/genome/${params.genome}", checkIfExists: true)
    .ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    .view()
    .set { bwa_index }

Channel
    .fromPath( './tmp', type: 'dir' )
    .set { tmpDir_ch}

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
    container 'mblanche/bwa-samtools'
    
    publishDir "/mnt/ebs/ref_push/${params.destDir}/${params.expName}/bam", mode: 'copy'
    
    input:
    set id, file(R1_files) from R1_groups_ch
    set id2, file(R2_files) from R2_groups_ch
    file index from bwa_index.first()
    path tmpDir from tmpDir_ch.first()
    
    output:
    path "*.bam" into sam_ch
    
    script:
    """
    bwa mem -t \$(nproc) \
	-5SP ${index}/${params.genome} \
	<(zcat $R1_files|head -n 100000) \
	<(zcat $R2_files|head -n 100000) \
	| samtools  view -@ \$(nproc) -Shb - \
	| samtools sort -m 2G -T ${tmpDir} -@ \$(nproc) - 
	> ${id}.bam
    """
}

sam_ch.view()

