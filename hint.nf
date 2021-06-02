#!/usr/bin/env nextflow

params.outDir = false

params.biosample = false
params.genewizMap = false
params.fastqDir = false

params.libraryID = false

params.genome = 'hg38'

params.resolutions = false
params.ABresolutions = false

params.HiChIP = false

params.help = false


def helpMessage() {
    log.info"""
    Usage:
    
    The typical command for running the pipeline is as follows:

        epiPipe --biosample bisample1,biosample2,biosample3  --outDir ~/ebs/ref_push/prodEpi/myExpName

         or

        epiPipe --fastqDir ~/ebs/ref_push/prodEpi/myExpName/fastqs/ --outDir ~/ebs/ref_push/prodEpi/myExpName

    Mandatory arguments:
        --outDir [pathOfDir]         Path of where to publish the data, can be local or an S3 accessible bucket

    Input (minimum one of these, can be mixed):
        --biosample [str]            Comma seperated list of BaseSpace biosamples.
        --fastqDir [str]             Path to a directoy of fastq files.
        --genewizMap [str]           CSV file descrbibing the LIMS id to S3 path on GeneWiz bucket. Field: LIMS_ID,LIBRAY_ID,S3_path_R1,S3_path_R2

    Restriting to a set of library:
        --libraryID [str]            Comma seperated list of library prefixes. Need to be used in conjunction with --fastqDir.
    
    Arrowhead parameters:
        --resolutions [integers]     Comma-seperated list of resolutions in kb for loops finding. Default: [5,10]
        --ABresolutions [integers]   Comma-seperated list of resolutions in kb. Default: [32,64,128]

    References:
        --genome [str]               Name of the genome to use. Possible choice: hg38, mm10, dm3. Default: hg38.
    
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.biosample || params.fastqDir || params.genewizMap) ){
    exit 1, "You need to use either --biosample aBaseSpace-biosample or --fastqDir path/to/bam/files or --genewizMap mappingFiles.csv. Use --help to get the full usage." 
}


if ( !(params.outDir) ){
    exit 1, "--outDir is a required arguments. Use --help to get the full usage." 
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10]
}

if (params.ABresolutions){
    ABresolutions = params.ABresolutions.toString().split(/,/,-1)
} else {
    ABresolutions = [32,64,128]
}

if (!params.genome =~ /hg38|mm10|dm3/){
    exit 1, "Only hg38, mm10 and dm3 genomes are currently offered"
} else {
    Channel
	.fromPath("${HOME}/ebs/genome/nextflow/${params.genome}", checkIfExists: true)
	.ifEmpty { exit 1, "BWA index not found: ${params.genome}" }
	.set { bwa_index }
    Channel
	.fromPath("${HOME}/ebs/genome/nextflow/${params.genome}/${params.genome}.fa", checkIfExists: true)
	.ifEmpty { exit 1, "Genome not found: ${params.genome}" }
	.set { abcomp_genome_ch }
    
}

if (params.biosample){
    Channel
	.from(params.biosample)
	.splitCsv()
	.flatten()
	.set { biosample_ch }
} else {
    Channel
	.empty()
    	.into { biosample_ch;fastqs_ch }
}

if (params.fastqDir){
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
	    .set{fastqDir_ch}
    } else {
	Channel
	    .fromFilePairs("${params.fastqDir}/*_{R1,R2}*.fastq.gz",
			   chekIfExists: true,
			   flat: true)
	    .set{fastqDir_ch}
	
    }
} else {
    if (params.libraryID){
	exit 1, "--libraryID needs to be used in conjuction with --fastqDir. Use --help to get details"
    } else {
	Channel
	    .empty()
	    .set{fastqDir_ch}
    }
}

if (params.genewizMap) {
    def i = 0
    Channel
	.fromPath(params.genewizMap)
	.splitCsv()
	.map{row ->
	    i=i+1
	    tuple(row[0]+"-rep"+i,file(row[2]),file(row[3]))
	}
	.into{keyMapping_ch;genewiz_ch}
    
    process copy_mapping {
	cpus 1
	memory "2 GB"
	container 'mblanche/basespace-cli'
	
	publishDir "${params.outDir}",
	    mode: 'copy'
	
	input:
	val mapping from keyMapping_ch
	    .map { l = [it[0],it[1].toString(),it[2].toString()]
		  l.join(",")
	    }
	    .flatten()
	    .collect()

	output:
	path ("*.csv")  into outMapping

	script:
	"""
	echo "id,R1_fastq,R2_fastq" > keyMapping.csv
	for l in ${mapping}; do
	echo \${l::-1} >> keyMapping.csv
	done
	"""
    }

    process copy_manifest {
	cpus 1
	memory "2 GB"
	container 'mblanche/basespace-cli'
	
	publishDir "${params.outDir}",
	    mode: 'copy'
	
	input:
	path manifest from Channel.fromPath(params.genewizMap)
	
	output:
	path manifest into outManifest
	
	script:
	"""
	cat ${manifest}
	"""
    }
    
}else{
    Channel
	.empty()
	.set{genewiz_ch}
}

if (params.biosample){
    process get_bs_files {
	cpus 1
	memory '1G'
	container 'mblanche/basespace-cli'
	
	input:
	val bs from biosample_ch
	
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
	cpus 1
	memory '1G'
	container 'mblanche/basespace-cli'
	
	publishDir "${params.outDir}/fastqs",
	mode: 'copy',
	    saveAs: {filename -> filename.endsWith('.error') ? filename : null}
	
	input:
	tuple file_name, bs, id from bs_id_ch
	    .splitCsv(header: true)
	    .map { row -> tuple(row.FilePath,row.biosample, row.Id )}
	    .groupTuple()
	    .map { tuple(it[0], it[1].unique()[0],it[2]) }
	
	output:			
	path "bs_data_set.csv" optional true into bs_dataset_ch
	path "*.error" optional true into error_ch
	
	script:
	if (id.size() >1){
	    """
	    echo found more than one identical file names for biosample ${bs} > ${bs}.error
	    """
	} else {
	    """
	    bs list dataset \
    		--like-type=illumina.fastq.v1.8 \
    		--input-biosample=${bs} \
    		--like-type=illumina.fastq.v1.8 \
    		-f csv | \
		awk 'BEGIN{OFS =","} \
		NR == 1 {print "biosample", \$0} \
		NR > 1  {print "${bs}", \$0}' > bs_data_set.csv 
	    """
	}
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
	tuple bs, val(bs_id) from bs_dataset_ch
            .splitCsv(header: true)
	    .unique()
	    .map{row -> tuple(row.biosample, row.Id)}
	
	output:
	tuple bs, file("*_R1*.fastq.gz"), file("*_R2*.fastq.gz") into fastqs_ch
	
	script:
	"""
	bs download dataset \
     	    --id=${bs_id} \
 	    -o . 
	"""
    }
}

process hint_pre {
    tag "_${id}"
    cpus 48
    memory '48 GB'
    container 'suwangbio/hint'
    
    input:
    tuple val(id), file(R1s), file(R2s) from fastqs_ch
    	.mix(fastqDir_ch)
	.mix(genewiz_ch)

    path index from bwa_index.first()
    
    output:
    tuple id, path("*.bam"), path("*.tsv") into  pairtools_parse_ch
    
    script:
    """

    hint pre -d <(zcat ${R1s|head -n 400000}),<(zcat ${R2s|head -n 400000}) \
	-i ${index}/${params.genome} \
	--refdir /path/to/refData/hg19 \
	-g hg19 \
	--informat fastq \
	--outformat cooler \
	-n test \
	-o /path/to/outputdir \
	--pairtoolspath /path/to/pairtools \
	--samtoolspath /path/to/samtools \
	--coolerpath /path/to/cooler\
    """
}

process pairtools_parse {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, path(sam), path(chr_sizes) from pairtools_parse_ch

    output:
    tuple id, path("*.pairsam.gz") into pairsam_part_ch

    script:
    """
    pairtools parse \
	--min-mapq 40 \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--chroms-path ${chr_sizes} \
	--output ${id}.pairsam.gz \
	${sam} 
    """
}

process pairtools_merge_lane {
    tag "_${id}"
    cpus 14
    memory '50 GB'
    container 'mblanche/pairtools'
    
    input:
    tuple id, path(sam) from pairsam_part_ch
	.map {id, file ->
            def key = id.tokenize('_').get(0)
            return tuple(key, file)
	}
	.groupTuple()

    output:
    tuple id, path("*.pairsam.gz") into pairsam_ch

    script:
    if (sam.sort().size() >1) {
	"""
	pairtools merge -o ${id}.pairsam.gz --nproc ${task.cpus} ${sam}
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
    tuple id, path(sam) from pairsam_ch

    output:
    tuple id, path("*_sorted.pairsam.gz") into sorted_ps_ch
    
    script:
    """
    mkdir -p tmp 
    pairtools sort --tmpdir ./tmp  \
	--nproc ${task.cpus} \
	--output ${id}_sorted.pairsam.gz \
	$sam 
    """
}

process pairtools_dedup {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    tuple id, path(sam) from sorted_ps_ch

    output:
    tuple id, path("*_dedup.pairsam.gz") into dedup_ps_ch
    tuple id, path("*_unmapped.pairsam.gz") into unmapped_ps_ch
    tuple id, path("*_pairtools.stats") into ps_stats_ch

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

process pairtools_stats_merge {
    tag "_${id}"
    cpus 1
    memory '4 GB'
    container 'mblanche/pt-stats'

    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    path(stats) from ps_stats_ch
	.map{it[1]}
	.collect()
    
    output:
    path('pairtoolsStats.csv') into merged_stats_ch
    
    script:
    """
    pairtoolsStat.sh ${stats} > pairtoolsStats.csv
    """
}

process pairtools_split_dedup {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'

    input:
    tuple id, path(sam) from dedup_ps_ch
    
    output:
    tuple id, path("*.bam") into bam_parts_ch
    tuple id, path("*.valid.pairs.gz") into pairs_parts_ch, pairs_parts_ch_test

    script:
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_PT.bam  \
	--output-pairs ${id}_PT.valid.pairs.gz  \
	${sam}
    """
}

process merge_bam {
    tag "_${id}"
    cpus 48
    memory '100 GB'
    container 'mblanche/bwa-samtools'

    input:
    tuple id, path(bam_part) from bam_parts_ch
	.map {id, file ->
	    if ( id.contains("-rep") ){
		def key = id.tokenize('-rep').get(0)
		return tuple(key, file)
	    } else {
		return( tuple(id,file) )
	    }
	}
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

process pairtools_split_unmapped {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/unmapped",
    	mode: 'copy'

    input:
    tuple id, path(sam) from unmapped_ps_ch
    
    output:
    path "*_unmapped.bam" into unmapped_bam_ch
    path "*_unmapped.valid.pairs.gz" into unmapped_pairs_ch

    script:
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${id}_unmapped.bam  \
	--output-pairs ${id}_unmapped.valid.pairs.gz  \
	${sam}
    """
}

process merge_pairs {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'mblanche/pairtools'
    
    publishDir "${params.outDir}/validPairs",
    	mode: 'copy'
    
    input:
    tuple val(id), path(pairs) from pairs_parts_ch
    .map {id, file ->
	if ( id.contains("-rep") ){
	    def key = id.tokenize('-rep').get(0)
	    return tuple(key, file)
	} else {
	    return( tuple(id,file) )
	}
    }
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

process arrowhead {
    tag "_${id}"
    cpus 12
    memory '40 GB'
    container "mblanche/juicer"

    publishDir "${params.outDir}/arrowHead",
	mode: 'copy'
    
    input:
    tuple id, path(hic), val(res) from arrowhead_ch
	.combine(Channel.from(resolutions))

    output:
    tuple id, path("${id}_${res}kb") into arrowhead_out_ch

    script:
    bpRes = res.toInteger() * 1000
    """
    mkdir -p ${id}_${res}kb && touch ${id}_${res}kb/${bpRes}_blocks.bedpe
    java -Xmx24000m \
	-jar /juicer_tools.jar \
	arrowhead \
	--threads ${task.cpus} \
	--ignore-sparsity \
	-r ${bpRes} \
	-k KR \
	${hic} \
	${id}_${res}kb
    """
}

process hiccups {
    tag "_${id}"
    label 'gpu'
    accelerator 1
    cpus 6
    memory '30 GB'
    container "mblanche/hiccups-gpu"

    publishDir "${params.outDir}/hiccups/",
	mode: 'copy'
    
    input:
    tuple id, path(hic), val(res)  from hiccups_ch
        .combine(Channel.from(resolutions.collect{it*1000}.join(',')))

    output:
    tuple id, path("${id}_loops") into hiccups_out_ch

    script:
    """
    java -Xmx24000m \
	-jar /juicer_tools.jar \
	hiccups \
	--threads ${task.cpus} \
	--ignore-sparsity \
	-m 500 \
	-r ${res} \
	-k KR \
	${hic} \
	${id}_loops
    """
}

process mustache {
    tag "_${id}"
    cpus 24
    memory '48 GB'
    container "mblanche/mustache:1.1.3"

    publishDir "${params.outDir}/mustache",
	mode: 'copy'

    input:
    tuple id, path(mcool), val(res)  from mustache_mcool_ch
	.combine(Channel.from(1000,4000,16000))
    
    output:
    tuple id, path("*.tsv") into mustache_2_merge_ch

    script:
    """
    touch ${id}_${res}kb_loops.tsv 
    mustache -p ${task.cpus} \
	-f ${mcool} \
	-r ${res} \
	-o ${id}_${res}kb_loops.tsv
    """
}

process ABcomp {
    tag "_${id}"
    cpus 1
    memory '12 GB'
    container "mblanche/fan-c"
    
    publishDir "${params.outDir}/AB_comp",
	mode: 'copy'

    input:
    tuple id, path(cool), val(resKB) from abcomp_mcool_ch
    	.combine(Channel.from(ABresolutions))
    
    path(genome) from abcomp_genome_ch.first()

    output:
    tuple id, path("*.bed"), path("*.ab") into fanc_out_ch
    
    script:
    res = resKB.toInteger() * 1000
    """
    fanc compartments \
	-f \
	-v ${id}_eigenV_${resKB}kb.bed \
	-d ${id}_AB_${resKB}kb.bed \
	-g ${genome} \
	${cool}@${res} \
	${id}_${resKB}kb.ab
    """
    
}
