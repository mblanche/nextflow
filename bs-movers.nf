params.destDir = 'capture'
params.expName = 'tmp'
params.genome = 'hg38'

params.bpProject = 'Capture'
params.biosample = 'PC_I14_29-51'

/*
// PC_T33_29-80-2,PC_T33_29-80-3,PC_T34_29-86-2,PC_T34_29-86-3,PC_T35_29-88-2,PC_T35_29-88-3,PC_T36_29-92-2,PC_T36_29-92-3,PC_T37_29-96-2,PC_T37_29-96-3,PC_T38_29-100-2,PC_T38_29-100-3,PC_T39_29-104-2,PC_T39_29-104-3,PC_T40_29-105-2,PC_T40_29-105-3
*/

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
    queue 'movers-Q'
    
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
