
params.mcoolDir=false
params.bedpe=false
params.outDir=launchDir
params.suffix=false
params.help = false

params.res=8000
params.size=1000

def helpMessage() {
    log.info"""
    Usage:
    
    The typical command for running the pipeline is as follows:

        hicPlotMatrix --mcool ~/path/to/mcoolFile.mcool --bedbe ~/path/to/chr1-chr2.bepde --outDir ~/ebs/ref_push/prodEpi/myExpName
    
    Mandatory arguments:
        --mcoolDir [path]            Path to a hic mcool file.
        --tloc [path]                Path to a bedpe containing pairs of genomic coordiate to plot.
    
    Optional arguments:
        --outDir [pathOfDir]         Path of where to publish the data, can be local or an S3 accessible bucket. Default: current dir "./"
        --suffix [string]            Extra suffix to add the file name to differentiate it from other run
    
    ##NOT USED FOR NOW
    Additional parameters:
        --mapQ [int]                 Quqlity score to filter reads. Integer between 0 and 60. Default: 0.
        --resolutions [integers]     Arrowhead and hiccup resolutions. Comma-seperated list of resolutions in kb for loops finding. Default: [5,10]
        --ABresolutions [integers]   ABcomp resolutions. Comma-seperated list of resolutions in kb. Default: [32,64,128]
    
""".stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if ( !(params.mcoolDir && params.bedpe) ){
            exit 1, "--mcoolDirDir and --bedpe are required arguments. Use --help to get the full usage." 
}

////////////////////////////////////////////////////////////////////////////////
// Creating channel first out of processes, make code cleaner
////////////////////////////////////////////////////////////////////////////////
// mcool channel
Channel.fromPath("${params.mcoolDir}/*.mcool",checkIfExists: true)
    .set{ mcool_ch }

// Process the bedpe into a channel
def i=1
Channel
    .fromPath("${params.bedpe}")
    .splitText()
    .splitCsv(sep:"\t",header:false)
    .map{row ->
	out = tuple(i,row[0],row[1],row[3],row[4],row[6])
	i=i+1
	return(out)
    }
    .set{ bedpe_ch }
////////////////////////////////////////////////////////////////////////////////

process plotMatrix {
    tag "_${id}"
    label "median"    
    container 'mblanche/hicexplorer'

    input:
    tuple iter, val(chr1), val(start1), val(chr2), val(start2), val(name), path(mcool)  from bedpe_ch
	.combine(mcool_ch)
	    
    output:
    tuple val(id), val(iter), path("*.png") into pngMerge_ch
    
    script:
    r1Start = start1.toInteger() - params.size * 1000
    r1End = start1.toInteger() + params.size * 1000
    r2Start = start2.toInteger() - params.size * 1000
    r2End = start2.toInteger() + params.size * 1000
    id = (mcool.name.toString().split(/\./))[0]
    """
    hicPlotMatrix -m ${mcool}::/resolutions/${params.res} \
	-out ${iter}.png \
	--title "${name}" \
	--region ${chr1}:${r1Start}-${r1End} \
	--region2 ${chr2}:${r2Start}-${r2End}
    """
}

process mergePNG {
    tag "_${id}"
    label "mezzo"    
    container 'mblanche/hicexplorer'
    
    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    tuple id, path(png) from pngMerge_ch
	.map{ col  -> tuple(col[0],[col[1],col[2]]) }
	.groupTuple()
	.map { id,file ->
	    file.sort{ a,b -> a[0] <=> b[0] }
	    return tuple(id,file) 
	}
	.transpose()
	.flatten()
	.collate(3)
	.groupTuple()
	.map{id,num,files -> tuple(id,files) }
    
    
    output:
    path("*.pdf")
    
    script:
    extraSuff = params.suffix ? "_" + params.suffix : ''
    fname = id + "_" + params.size + "k" + extraSuff + ".pdf"
    """
    convert ${png} ${fname}
    """
}
