
/*
nextflow run test.nf --expName Twist-PG-plexing --expDir capture --biosample $(echo PC-OM-GM-4P-{1..4} PC-OM-GM-8P-{1..8} PC-OM-GM-2P-{1..2}|tr ' ' ,)
*/
params.expDir = 'capture'
params.expName = 'tmp'

params.bpProject = 'Capture'

params.biosample = false
params.bamDir = false
params.validPairsDir = false
params.fastqDir = false


params.genome = 'hg38'

params.noPairTools = false
params.noJuicer = false
params.noCooler = false
params.noCoverage = false

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
    val bs from biosample_ch.first()

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


process get_bs_dataset {
    echo true
    label 'movers'
    cpus 4
    memory '2G'
    container 'mblanche/basespace-cli'
    
    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/fastqs",
	mode: 'copy'
    
    input:
    tuple bs, id, filePath from bs_id_ch
	.splitCsv(header: true)
	.map {row -> tuple(row.biosample,row.Id,row.FilePath) }

    output:			
    path "*.fastq.gz" into fastq_ch
    
    script:
    """
    bs file download -i ${id} -o .
    
    new_fq=\$(echo ${filePath}| perl -pe 's/.+?(_.+)/${bs}\$1/')

    if [[ ${filePath} != \$new_fq ]]
    then
    mv ${filePath} \$new_fq
    fi
    """
}

fastq_ch.view()
