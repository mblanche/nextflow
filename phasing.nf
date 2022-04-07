
params.help   = false
params.design = false
params.outDir = false
params.genome = 'hg38'


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        phasing.nf  --design ~/path/to/bam/location/design.csv --outDir ~/path/to/results/dir
    
    Mandatory arguments:
        --outDir [path]              Path to a diectory to save the bigwigfiles (can be local or valid S3 location).

        --design [path]              Path to a design files in csv. First col: id, second col: replicate, third col: path to R1, fourth col: path to R2
    
    Alignment:
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.



    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.outDir && params.design)) {
    exit 1, "--outDir and --design is are required arguments. Use --help to get the full usage." 
} else {
    outDir = params.outDir
}


if (!params.genome =~ /hg19|hg38|mm10|dm3/){
    exit 1, "Only hg38, mm10 and dm3 genomes are currently offered"
} else {
    Channel
    	.fromFilePairs("${HOME}/ebs/genome/nextflow/bwa-mem2-index/${params.genome}/*.{0123,amb,ann,bwt.2bit.64,pac}",
    		       size: -1,
    		       checkIfExists: true)
    	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
    	.set { bwa_index }
    
    Channel
	.fromFilePairs("${HOME}/ebs/genome/nextflow/${params.genome}/*.{fa,fa.fai}", checkIfExists: true, flat: true)
	.ifEmpty { exit 1, "Genome not found: ${params.genome}" }
	.set {genome_ch }
    
}

process bwa_mem {
    tag "_${prefix}"
    label 'cpu'
    container 'mblanche/bwa-mem2'
    
    publishDir "${params.outDir}/bam",
	mode: 'copy'

    input:
    tuple id, val(rep), path(R1s), path(R2s) from Channel
	.fromPath(params.design)
	.splitCsv(header: false)
    
    tuple index, path(index_files) from bwa_index.first() 
    
    output:
    tuple id, path("*.bam"), path("*.bai") into bam_dv_ch,bam_hairs_ch
    path("*.bam") into bam_chr_ch
    
    script:
    prefix = R1s.name.toString().replaceFirst(/.fastq.+/,"")
    """
    bwa-mem2 mem -5SP -t ${task.cpus} \
    	${index} \
    	${R1s} \
	${R2s} \
	|samblaster -r \
	|samtools view -@ ${task.cpus} -Shb - \
	|samtools sort -m 2G -@ ${task.cpus} -o ${prefix}.bam - \
	&& samtools index ${prefix}.bam 
    """
 }


process dv_simple {
    label 'cpu'
    container "google/deepvariant:1.3.0"
    
    publishDir "${params.outDir}/deepVariant",
	mode: 'copy'
    
    input:
    tuple id, path(bam), path(idx) from bam_dv_ch
    tuple ref, path(fa), path(fai) from genome_ch.first()
    
    output: 
    path "${id}.vcf" into hairs_vcf_ch,hc2_vcf_ch
    path "*.gvcf"
    path "${id}_raw.vcf"
    
    script:
    """ 
    /opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=${fa} \
	--reads=${bam} \
	--output_vcf=${id}_raw.vcf \
	--output_gvcf=${id}.gvcf \
	--intermediate_results_dir /intermediate_results_dir \
	--num_shards=${task.cpus} \
	&& cat ${id}_raw.vcf | grep -E '^#|0/0|CHROM|1/1|0/1|1/0|0/2|2/0' -w > ${id}.vcf
    """
}

process extractChr {
    tag "_${prefix}"
    label 'mezzo'
    container 'mblanche/bwa-mem2'
    
    input:
    path(bam) from bam_chr_ch
    
    output:
    path("*.txt") into chrInfo_ch
    
    script:
    """
    samtools view -H ${bam} \
	|awk -v OFS='\t' '/^@SQ/ && !(\$2 ~ /:(chr|"")M/) {split(\$2,chr,":");split(\$3,ln,":");print chr[2]":1-"ln[2]}' \
	|sort -V -k1,1 \
	> chrInfo.txt
    """
}

process extractHairs {
    tag "_${id}"
    label 'median'
    container 'mblanche/hapcut2'
        
    input:
    tuple val(region), val(id), path(bam), path(idx), path(vcf) from chrInfo_ch
    .splitText(){it.trim()}
    .combine(bam_hairs_ch)
    .combine(hairs_vcf_ch)

    output:
    tuple id, val(chr), path("*.hairs") into hairs_hc2_ch
    
    script:
    chr= region.split(":",2)[0]
    """
    extractHAIRS \
	--bam ${bam} \
	--VCF ${vcf} \
	--hic 1 \
	--region ${region} \
	--out ${id}_${chr}.hairs
    """
}

process hapcut2 {
    tag "_${id}"
    label 'alto'
    container 'mblanche/hapcut2'
    
    publishDir "${params.outDir}/hapcut_res",
	mode: 'copy'
    
    input:
    tuple id, val(chr), path(hairs), path(vcf) from hairs_hc2_ch
	.groupTuple()
        .combine(hc2_vcf_ch)
    
    output:
    path("*")
    
    script:
    """
    cat ${hairs} > all.hairs && \
	HAPCUT2 \
	--hic 1 \
	--fragments all.hairs \
	--VCF ${vcf} \
	--output ${id} \
	--outvcf 1
    """
}
