

Channel
    .fromPath("${HOME}/ebs/ref_push/working/DHockemeyer_cellLine/hicFiles/DHAM001-A.hic")
    .set{hic_ch}


process hiccups {
    echo true
    accelerator 1
    container "mblanche/hiccups"
    
    input:
    path hic from hic_ch
    
    script:
    id = hic.name.toString().take(hic.name.toString().lastIndexOf('.'))
    """
    java \
	-Xms512m \
	-Xmx4G \
	-jar /juicer_tools.jar \
	hiccups \
	-m 1024 \
	-r 5000,10000 \
	
	--threads 0 \
	--ignore-sparsity \
	${hic} \
	${id}.loops
    """
}

