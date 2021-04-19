

params.expDir = 'working'
params.expName = 'DHockemeyer_cellLine'

Channel
    .fromPath("${HOME}/ebs/ref_push/${params.expDir}/${params.expName}/coolerFiles/**.cool")
    .view()
    .set{cool_ch}

process mustache {
    echo true
    cpus 24
    memory '100 GB'
    container "mblanche/mustache"

    publishDir "mustache",
	mode: 'copy'
    
    input:
    path(cool) from cool_ch
	.first()

    //output:
    //path "*.tsv" into loops_ch
    //path "*.txt" into time_ch
    
    script:
    id = cool.name.toString().take(cool.name.toString().lastIndexOf('.'))
    """
    mustache -p ${task.cpus} -f ${cool} -r 5kb -o ${id}_5kb_loops.tsv -ch chr21
    """
}
