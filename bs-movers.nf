params.expDir = 'capture'
params.expName = 'tmp'

params.bpProject = 'Capture'
params.biosample = 'PC_I14_29-51'


Channel
    .from(params.biosample)
    .splitCsv()
    .flatten()
    .set { biosample_ch }


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
    .view()
    .set { out_ch }

process download {
    cpus 4
    container 'mblanche/basespace-cli'
    queue 'moversQ'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/fastqs", mode: 'copy'
    
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

