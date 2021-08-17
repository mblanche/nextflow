#!/usr/bin/env nextflow


params.fastqDir = false
params.design = false
params.outDir = false
params.help   = false
params.genome = 'hg38'
params.libraryID = false
params.mapQ = 40
params.bsDesign = false

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        alignPerDesign.nf  --design ~/path/to/bam/location/design.csv --fastqDir ~/path/to/bam/location

         or

        alignPerDesign.nf  --design ~/path/to/bam/location/design.csv --fastqDir ~/path/to/bam/location --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

        --design [path]              Path to a design files in csv. First col: id, second col: replicate, third col: prefix of fastq
    
    Alignment:
        --genome [str]               Name of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.

    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

if (!params.outDir || !(params.design ||params.bsDesign)) {
        exit 1, "--outDir and --design is are required arguments. Use --help to get the full usage." 
}

if(!params.outDir){
     outDir = file(params.fastqDir).getParent()
 } else {
     outDir = params.outDir
}


if (!params.genome =~ /hg19|hg38|mm10|dm3/){
    exit 1, "Only hg38, mm10 and dm3 genomes are currently offered"
} else {
    Channel
	.fromFilePairs("${HOME}/ebs/genome/nextflow/${params.genome}/*.{amb,sa,pac,ann,bwt,fa}", size: -1, checkIfExists: true)
	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
	.set { bwa_index }
}


if (params.design){
    Channel
    .fromPath(params.design)
    .splitCsv(header: true)
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
		    out.add(tuple(row.sample,row.rep,row.lib,R1,R2[0]))
	    }
	}
	return(out)
    }
    .flatten()
    .collate(5)
    .set{fq_ch}

} else {
    Channel
	.empty()
	.set{fq_ch}
}

if (params.bsDesign){
    process get_bs_files {
	cpus 1
	memory '1G'
	container 'mblanche/basespace-cli'
	
	input:
	tuple bs,val(sample),val(rep) from Channel
	    .fromPath(params.bsDesign)
	    .splitCsv(header: true)
	    .map{tuple(it.bs, it.sample, it.rep)}
	
	output:
	stdout into bs_id_ch
	
	script:
	"""
	bs biosample content -n ${bs} -F Id -F FilePath -f csv | \
	    awk -v bs=${bs} -v sample=${sample}, -v rep=${rep} \
	    'BEGIN{OFS =","} \
	    NR == 1 {print "bs","sample","rep",\$0} \
	    NR > 1  {print bs,sample,rep,\$0}'
	"""
    }
    
    process download_bs {
	label "movers"
	cpus 4
	memory '4G'
	container 'mblanche/basespace-cli'
	queue 'moversQ'
	
	publishDir "${params.outDir}/fastqs",
	    mode: 'copy'
	
	input:
	tuple bs, val(oriFname), val(newFname), val(id) from bs_id_ch
	    .splitCsv(header: true)
	    .map { row -> tuple(row.FilePath,row.biosample, row.Id )}
	    .groupTuple()
	    .map{if (it[2].size() >1){
		    x = []
		    fname = it[0]
		    id = fname.take(fname.indexOf('_'))
		    suff = fname.substring(fname.indexOf('_')+1)
		    bs = it[1][0]
		    
		    for (i in (1..it[2].size())) {
			x.add([bs,fname,"${id}_FC${i}_${suff}",it[2][i-1]])
		    }
		    return(x)
		} else {
		    fname = it[0]
		    bs = it[1][0]
		    id = it[2]
		    return([bs,fname,fname,id])
		}
	    }
	    .flatten()
	    .collate(4)
	
	output:
	tuple bs, file("*.fastq.gz") into bs_ch
	
	script:
	"""
	bs file download -i ${id} -o . 

	if [ "${oriFname}" != "${newFname}" ];then 
	    mv ${oriFname} ${newFname}
	fi
	"""
    }
} else {
    Channel
	.empty()
	.set{bs_ch}
}

process bwa_mem {
    tag "_${id}"
    cpus 12
    memory '24 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, val(rep), val(lib), file(R1s), file(R2s) from fq_ch
	.mix(bs_ch)
    
    tuple index, path(index_files) from bwa_index.first() //This is an hack to make sure all files are in staged area
    
    output:
    tuple id, val(rep), val(lib), val(prefix), path("*.bam"), path("*.tsv") into  pairtools_parse_ch
    
    script:
    prefix = R1s.name.toString().replaceFirst(/.fastq.+/,"")
    """
    bwa mem -5SP -t ${task.cpus} \
    	${index} \
    	<(zcat ${R1s}|head -n 8000) \
    	<(zcat ${R2s}|head -n 8000) \
	|samtools view -@ ${task.cpus} -Shb -o ${prefix}.bam - \
	&& samtools view -H ${prefix}.bam | \
	awk -v OFS='\t' '/^@SQ/ && !(\$2 ~ /:(chr|"")M/) {split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' | \
	sort -V -k1,1 \
	> chr_size.tsv
    """
}

process pairtools_parse {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, val(rep), val(lib), val(prefix), path(sam), path(chr_sizes) from pairtools_parse_ch

    output:
    tuple id, val(rep), val(lib), path("*.pairsam.gz") into pairsam_part_ch

    script:
    """
    pairtools parse \
	--min-mapq ${params.mapQ} \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--chroms-path ${chr_sizes} \
	--output ${prefix}.pairsam.gz \
	${sam} 
    """
}

process pairtools_merge_lane {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'
    
    input:
    tuple id, val(rep), val(lib), path(sam) from pairsam_part_ch
	.groupTuple(by: [1,0])
    
    output:
    tuple id, val(rep), path("*.pairsam.gz") into pairsam_ch

    script:
    if (sam.sort().size() >1) {
	"""
	pairtools merge -o ${id}-${rep}.pairsam.gz --nproc ${task.cpus} ${sam}
	"""
    } else {
	"""
	ln -sf ${sam} ${id}_ML.pairsam.gz
	"""
    }
    
}

process pairtools_sort {
    tag "_${id}"
    cpus 14
    memory '100 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, val(rep), path(sam) from pairsam_ch

    output:
    tuple id, val(rep), path("*_sorted.pairsam.gz") into sorted_ps_ch
    
    script:
    """
    mkdir -p tmp 
    pairtools sort --tmpdir ./tmp  \
	--nproc ${task.cpus} \
	--output ${id}-${rep}_sorted.pairsam.gz \
	$sam 
    """
}

sorted_ps_ch
    .view()
/*
process pairtools_dedup {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    tuple id, val(rep), path(sam) from sorted_ps_ch

    output:
    tuple id, val(rep), path("*_dedup.pairsam.gz") into dedup_ps_ch
    tuple id, val(rep), path("*_unmapped.pairsam.gz") into unmapped_ps_ch
    tuple id, val(rep), path("*_pairtools.stats") into ps_stats_ch

    script:
    """
    pairtools dedup --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--mark-dups \
	--output-stats ${id}_pairtools.stats  \
	--output ${id}_dedup.pairsam.gz \
	--output-unmapped ${id}_unmapped.pairsam.gz \
	${sam}
    """
}

process pairtools_split_dedup {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, val(rep), path(sam) from dedup_ps_ch
    
    output:
    tuple id, val(rep), path("*.bam") into bam_parts_ch
    tuple id, val(rep), path("*.valid.pairs.gz") into pairs_parts_ch, pairs_parts_ch_test

    script:
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}-${rep}_PT.bam  \
	--output-pairs ${id}-${rep}_PT.valid.pairs.gz  \
	${sam}
    """
}




process merge_bam {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'
    
    input:
    tuple id, val(rep), path(bam_part) from bam_parts_ch
	.groupTuple()

    output:
    tuple id, path("*.bam") into merged_bam_sort_ch

    script:
    bam_files = bam_part.sort()
    if (bam_files.size() >1) {
	"""
	samtools merge -@ ${task.cpus} ${id}_MB.bam ${bam_part}
	"""
    } else {
	"""
	ln -s ${bam_part} ${id}_MB.bam
	"""
    }
}

process merge_pairs {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/validPairs",
    	mode: 'copy'
    
    input:
    tuple val(id), val(rep), path(pairs) from pairs_parts_ch
	.groupTuple()
    
    output:
    tuple id, path("*.valid.pairs.gz"), path("*.px2") into pairs_chrSize_ch

    script:
    pair_files = pairs.sort()
    if (pair_files.size() >1) {
	"""
	pairtools merge -o ${id}.valid.pairs.gz --nproc ${task.cpus}  ${pairs}
	pairix ${id}.valid.pairs.gz
	"""
    } else {
	"""
	ln -s ${pairs} ${id}.valid.pairs.gz
	pairix ${id}.valid.pairs.gz
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
    tuple id, path(bam) from merged_bam_sort_ch
    
    output:
    tuple id, path("${id}.bam"),path("${id}.bam.bai") into bam_bigwig_ch

    script:
    """
    samtools sort -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} 

    samtools index -@${task.cpus} ${id}.bam
    """
}

process chr_size {
    tag "_${id}"
    cpus 4
    memory '24 GB'
    container 'mblanche/pairtools'
    
    input:
    tuple id, path(pairs), path(idx) from pairs_chrSize_ch
    
    output:
    tuple id, path(pairs), path(idx), path("*.tsv") into pairs_ch_cooler, pairs_ch_juicer
    
    script:
    """
    pairix -H -f ${pairs} \
	| awk -v OFS='\t' '/^#chromsize/  {print \$2,\$3}' \
	| sort -V -k1,1 \
	> chr_size.tsv
    """
}


process cooler_cload {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'

    input:
    tuple id, path(pairs), path(idx), path(chr_sizes) from pairs_ch_cooler
    
    output:
    tuple id, path("*.cool") into balance_cooler_ch
        
    script:
    """
    cooler cload pairix \
	-p ${task.cpus} \
	${chr_sizes}:1000 \
	${pairs} \
	${id}.cool
    """
}

process balance_cooler {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${params.outDir}/coolerFiles",
    	mode: 'copy'
    
    input:
    tuple id, path(cooler) from balance_cooler_ch

    output:
    tuple id, path(cooler) into zoomify_cooler_ch 
    
    script:
    """
    cooler balance --force -p ${task.cpus} ${cooler}
    """
}

process cooler_zoomify {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/cooler'
    
    publishDir "${params.outDir}/coolerFiles",
    	mode: 'copy'
    
    input:
    tuple id, path(cooler) from zoomify_cooler_ch

    output:
    tuple id, path("*.mcool") into mustache_mcool_ch, abcomp_mcool_ch
    
    script:
    """
    cooler zoomify --balance -p ${task.cpus} ${cooler}
    """
}

process bam2bw {
    tag "_${id}"
    cpus 20
    memory '175 GB'
    
    container 'mblanche/r-cov'
    
    publishDir "${params.outDir}/bigwigs",
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

process juicer {
    tag "_${id}"
    cpus 24
    memory '150 GB'
    container 'mblanche/juicer'
    
    publishDir "${params.outDir}/hicFiles",
    mode: 'copy'
    
    input:
    tuple id, path(pairs), path(idx), path(chr_sizes) from pairs_ch_juicer
    
    output:
    tuple id, path("*.hic") into arrowhead_ch, hiccups_ch

    script:
    """
    java -Xmx96000m -Djava.awt.headless=true \
	-jar /juicer_tools.jar pre \
	--threads ${task.cpus} \
	-j ${task.cpus} \
	-k VC,VC_SQRT,KR,SCALE \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    """
}
*/
