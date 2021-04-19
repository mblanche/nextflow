

params.expDir = 'working'
params.expName = 'DHockemeyer_cellLine'


Channel.from(10,25).set{res}

Channel
    .fromPath("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/hicFiles/*.hic")
    .first()
    .set{arrowHead_ch}


process arrowhead {
    echo true
    cpus 24
    memory '24 GB'
    container "mblanche/juicer"

    publishDir "arrowhead",
	mode: 'copy'
    
    input:
    tuple path(hic), val(res)  from arrowHead_ch
        .combine(Channel.from(10,25))

    output:
    path "*.bed" into arrowhead_ch
    path "*.log" into arrowhead_log_ch
    
    script:
    id = hic.name.toString().take(hic.name.toString().lastIndexOf('.'))
    """
    java -Xmx24000m -Djava.awt.headless=true \
	-jar /juicer_tools_1.22.01.jar arrowhead \
	--threads ${task.cpus} \
	--ignore-sparsity \
	-r \$(( ${res} * 1000 )) \
	-k KR \
	${hic} \
	${id}.tads.${res}kb.bed &> ${id}_${res}kb.log ||
	touch ${id}.tads.${res}kb.bed
    """
}

 

