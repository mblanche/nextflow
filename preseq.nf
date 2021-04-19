



params.expDir = 'capture'
params.expName = 'Twist-PG-plexing'

params.fastqDir="${HOME}/ebs/ref_push/capture/Twist-PG-plexing/fastqs"
params.genome = 'hg38'

Channel
    .fromPath("${HOME}/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
    .ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    .set { bwa_index }

Channel
    .fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
		   chekIfExists: true,
		   flat: true)
    .set{fastqDir_ch}


process bwa_mem {
    tag "_${id}"
    cpus 48
    memory '48 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple val(id), file(R1s), file(R2s) from fastqDir_ch
    path index from bwa_index.first()
    
    output:
    tuple id, path("*.bam") into bam_part_ch
    
    script:
    """
    bwa mem -5SP -t ${task.cpus} \
    	${index}/${params.genome} \
    	<(zcat ${R1s}) \
    	<(zcat ${R2s}) \
	|samtools view -@ ${task.cpus} -Shb - > ${id}.bam
    
    """
}

process merge_lane {
    echo true
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'
    
    input:
    tuple id, path(sam) from bam_part_ch
	.map {id, file ->
            def key = id.tokenize('_').get(0)
            return tuple(key, file)
	}
	.groupTuple()

    output:
    tuple id, path("*.bam"), path("*.bai") into bam_ch

    script:
    if (sam.sort().size() >1) {
	"""
	samtools merge -@ ${task.cpus} - ${sam} \
	    |samtools sort -@ ${task.cpus} -o ${id}_MB.bam \
	    && samtools index -@ ${task.cpus} ${id}_MB.bam
	"""
    } else {
	"""
	ln -s ${sam} ${id}_MB.bam
	samtools index -@ ${task.cpus} ${id}_MB.bam
	"""
    }
    
}


process preSeq {
    tag "${id}"
    cpus 48
    memory '150 GB'
    container 'mblanche/preseq'
    maxRetries 4

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}",
	mode: 'copy'
    
    input:
    tuple id, path(bam), path(idx) from bam_ch
    
    output:
    tuple id, path("${id}_cmplx_stat") into preseq_ch
    
    script:
    if (task.attempt == 3){
	"""
	touch "${id}_cmplx_stat"
	"""
    } else {
	"""
	preseq lc_extrap -B -P -e 2.1e9 -s 1e8 -o ${id}_cmplx_stat ${bam} 
	"""
    }
}

process merge {
    tag "${id}"
    cpus 1
    memory "2 GB"

    container 'ubuntu:20.04'

    publishDir "${HOME}/ebs/ref_push/${params.expDir}/${params.expName}",
	mode: 'copy'
    
    input:
    tuple id, path(cmplx) from  preseq_ch.collect()

    output:
    path("*.csv") into res
    
    script:
    """
    echo "Library,Complexity_@400M" > complexity.csv
    for f in ${cmplx};do
    echo \${f%_cmplx_stat},\$(grep '^400000000.0' \$f|cut -d \$'\t' -f2) >> complexity.csv
    done
    """
}

