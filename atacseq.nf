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
        --fastqDir [path]            Path to a direcoty with bam files and their index to compute coverages (can be local or valid S3 location.
                                     Mulitple paths can be specified by seperating them with commans

        --design [path]              Path to a design files in csv. First col: id, second col: replicate, third col: prefix of fastq
    
    Facultative arguments
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

    Restriting to a set of library:
        --libraryID [str]            Comma seperated list of library prefixes. Need to be used in conjunction with --fastqDir.
    
    Alignment:
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.

    """.stripIndent()
}

if (!params.fastqDir || !params.design ) {
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



line.split(/,/, -1).each{println "Hello $it"}

println(params.fastqDir 

//Channel.fromPath([ params.fastqDir ]).view()

/*


if (params.libraryID){
    Channel
	.from(params.libraryID)
	.splitCsv()
	.flatten()
	.map {id -> file("${params.fastqDir}/${id}_*R1*.fastq.gz")}
	.flatten()
	.map {R1 ->
	    def m = R1.name.toString() =~ /(.+)_R1/
		def id = m[0][1]
	    def R2 = file("${params.fastqDir}/${id}*_R2*.fastq.gz")
	    if (R2.size() > 1) {
		return( "R2 has more than one entry")
	    } else {
		return(tuple(id, R1, R2[0]) )
	    }
	    }
	.set{fq_ch}
} else {
    Channel
	.fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
		       chekIfExists: true,
		       flat: true)
	.set{fq_ch}
    
}


process NGmerge {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/ngmerge'

    publishDir "${outDir}/mergedFQ",
	mode: 'copy'
    
    input:
    tuple id, path(R1), path(R2)  from fq_ch

    output:
    tuple id, path("*merged_1.fastq.gz"), path("*merged_2.fastq.gz") into fastqs_ch
    
    script:
    prefix = R1.name.toString().replaceFirst(/.fastq.+/,"")
    """
    NGmerge  -a  -1 ${R1}  -2 ${R2} -o ${prefix}_merged -v -n ${task.cpus}
    """
}

process bwa_mem {
    tag "_${id}"
    cpus 36
    memory '24 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, file(R1s), file(R2s) from fastqs_ch
    tuple index, path(index_files) from bwa_index.first() //This is an hack to make sure all files are in staged area
    
    output:
    tuple id, path("*.bam") into bam_lane_ch
    
    script:
    """
    bwa mem -5SP -t ${task.cpus} \
    	${index} \
    	<(zcat ${R1s}) \
    	<(zcat ${R2s}) \
	|samtools view -@ ${task.cpus} -Shb -o ${id}.bam - \
    """
}

process bam__merge_lane {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, path(bam) from bam_lane_ch
	.map {id, file ->
            def key = id.tokenize('_').get(0)
            return tuple(key, file)
	}
	.groupTuple()

    output:
    tuple id, path("*_ML.bam") into bam_ML_ch

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
    
    publishDir "${params.outDir}/bam",
	mode: 'copy',
	pattern: "${id}.bam"
        
    input:
    tuple id, path(bam) from bam_ML_ch
    
    output:
    tuple id, path("${id}.bam") into bam_final_ch

    script:
    """
    samtools sort -n -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} 
    """
}

process peak_calling {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/ngmerge'

    publishDir "${outDir}/mergedFQ",
	mode: 'copy'
    
    input:
    tuple id, path(bam) from bam_final_ch

    output:
    tuple id, path("*.bed") into peaks_ch
    
    script:
    """
    Genrich  -t ${bam}  -o ${id}_peaks.bed  -j  -y  -r  -e chrM  -v
    """
}
*/
