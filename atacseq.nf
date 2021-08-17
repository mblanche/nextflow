#!/usr/bin/env nextflow


params.fastqDir = false
params.design = false
params.outDir = false
params.help   = false
params.genome = 'hg38'
params.libraryID = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        atacseq.nf  --design ~/path/to/bam/location/design.csv --fastqDir ~/path/to/bam/location

         or

        atacseq.nf  --design ~/path/to/bam/location/design.csv --fastqDir ~/path/to/bam/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

        --design [path]              Path to a design files in csv. First col: id, second col: replicate, third col: prefix of fastq
    
    Alignment:
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.

    """.stripIndent()
}

if (!params.outDir || !params.design ) {
        exit 1, "--fastqDir and --design is are required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
     outDir = file(params.fastqDir).getParent()
 } else {
     outDir = params.outDir
}

if (params.help){
    helpMessage()
    exit 0
}

if (!params.genome =~ /hg19|hg38|mm10|dm3/){
    exit 1, "Only hg38, mm10 and dm3 genomes are currently offered"
} else {
    Channel
	.fromFilePairs("${HOME}/ebs/genome/nextflow/${params.genome}/*.{amb,sa,pac,ann,bwt,fa}", size: -1, checkIfExists: true)
	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
	.set { bwa_index }
}

Channel
    .fromPath(params.design)
    .splitCsv(header: ['id','rep','lib','fqDir'])
    .map {row ->
	R1s=file("${row.fqDir}/${row.lib}_*R1*.fastq.gz")
	out=[]
	R1s.each { R1 ->
	    def m = R1.name.toString() =~ /(.+)_R1/
	    def prefix = m[0][1]
	    def R2 = file("${row.fqDir}/${prefix}*_R2*.fastq.gz")
	    if (R2.size() > 1) {
		return( "R2 has more than one entry")
	    } else {
		out.add(tuple(row.id,row.rep,R1,R2[0]))
	    }
	}
	return(out)
    }
    .flatten()
    .collate(4)
    .set{fq_ch}

process NGmerge {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/ngmerge'
    
    input:
    tuple id, val(rep), path(R1), path(R2)  from fq_ch

    output:
    tuple id, val(rep), val(prefix), path("*merged_1.fastq.gz"), path("*merged_2.fastq.gz") into fastqs_ch
    
    script:
    prefix = R1.name.toString().replaceFirst(/.fastq.+/,"")
    """
    NGmerge -a -1 ${R1} -2 ${R2} -o ${prefix}_merged -v -n ${task.cpus}
    """
}

process bwa_mem {
    tag "_${id}"
    cpus 36
    memory '24 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, val(rep), val(prefix), file(R1s), file(R2s) from fastqs_ch
    tuple index, path(index_files) from bwa_index.first() //This is an hack to make sure all files are in staged area
    
    output:
    tuple groupID, path("*.bam") into bam_lane_ch
    
    script:
    groupID = "${id}-rep${rep}"
    """
    bwa mem -5SP -t ${task.cpus} \
    	${index} \
    	<(zcat ${R1s}) \
    	<(zcat ${R2s}) \
	|samtools view -@ ${task.cpus} -Shb -o ${prefix}.bam - \
    """
}


process bam_merge_lane {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, path(bam) from bam_lane_ch
	.groupTuple()
    
    output:
    tuple id, path("*_ML.bam") into bam_ML_ch, bam_sort4bw_ch

    script:
    if (bam.sort().size() >1) {
	"""
	samtools merge -@ ${task.cpus} ${id}_ML.bam ${bam}
	"""
    } else {
	"""
	ln -sf ${bam} ${id}_ML.bam
	"""
    }
    
}

process bam_sort {
    tag "bam_sort_${id}"
    cpus 48
    memory '150 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "${outDir}/bam",
	mode: 'copy',
	pattern: "${id}.bam"
        
    input:
    tuple id, path(bam) from bam_ML_ch
    
    output:
    tuple groupID, path("${id}.bam") into bam_final_ch

    script:
    groupID = id.replaceFirst(/-rep.+/,"")
    """
    samtools sort -n -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} 
    """
}


process bam_sort_4bw {
    tag "bam_sort_${id}"
    cpus 48
    memory '150 GB'
    container 'mblanche/bwa-samtools'
    
    publishDir "${outDir}/bam",
	mode: 'copy',
	pattern: "${id}.bam"
        
    input:
    tuple id, path(bam) from bam_sort4bw_ch
    
    output:
    tuple groupID, path("${id}.bam"),path("*.bai") into bam_bigwig_ch

    script:
    groupID = id.replaceFirst(/-rep.+/,"")
    """
    samtools sort -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} 

    samtools index -@ ${task.cpus} ${id}.bam
    """
}

process bam2bw {
    tag "_${id}"
    cpus 20
    memory '175 GB'
    
    container 'mblanche/r-cov'
    
    publishDir "${outDir}/bigwigs",
    	mode: 'copy'
    
    input:
    tuple id, path(bam),path(idx) from bam_bigwig_ch
        
    output:
    tuple id, path ("*.bw") into bigwig_out_ch

    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}

process peak_calling {
    echo true
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/ngmerge'

    publishDir "${outDir}/mergedFQ",
	mode: 'copy'
    
    input:
    tuple id, path(bam) from bam_final_ch.groupTuple()

    output:
    tuple id, path("*.bed") into peaks_ch
    
    script:
    """
    Genrich  -t \$(echo ${bam}|tr ' ' ,)  -o ${id}_peaks.bed  -j  -y  -r  -e chrM  -v
    """
}

