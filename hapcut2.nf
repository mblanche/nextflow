
params.help = false
params.bam = false
params.idx = false
params.vcf = false
params.outDir = false


process extractChr {
    tag "_${prefix}"
    label 'mezzo'
    container 'mblanche/bwa-mem2'
    
    publishDir "${params.outDir}/hairs",
	mode: 'copy'
    
    input:
    path(bam) from Channel.fromPath("${params.bam}")
    
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
    label 'mezzo'
    container 'mblanche/hapcut2'
    
    publishDir "${params.outDir}/hairs",
	mode: 'copy'
    
    input:
    tuple val(region), path(bam), path(idx), path(vcf)  from  chrInfo_ch
	.splitText(){it.trim()}
	.combine(Channel.fromPath("${params.bam}"))
	.combine(Channel.fromPath("${params.idx}"))
	.combine(Channel.fromPath("${params.vcf}"))

    output:
    path("*.hairs")
    
    script:
    id = bam.name.toString().replaceFirst(/.bam.+/,"")
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

