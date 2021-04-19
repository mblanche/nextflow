params.expDir = 'capture'
params.expName = 'tmp'

params.bpProject = 'Capture'
params.biosample = 'PC-OM-GM-4P-1'

Channel
    .from(params.biosample)
    .splitCsv()
    .flatten()
    .set { biosample_ch }

process findFile {
    cpus 2
    memory '2G'
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
    memory '4G'
    container 'mblanche/basespace-cli'
    queue 'moversQ'
    
    publishDir "/mnt/ebs/ref_push/${params.expDir}/${params.expName}/fastqs", mode: 'copy'
    
    input:
    val id from out_ch.flatMap { it["Id"] }

    output:
    tuple id, path('*_R1*.fastq.gz'), path('*_R2*.fastq.gz') into fastqs_ch
    
    """
    bs download dataset \
     	--id=${id} \
 	-o .
    """
}

fastqs_ch
    .map { id, file1, file2 -> tuple(file1.name.toString().split('_')[0], file1, file2) }
    .groupTuple()
    .view()
