#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3){
    cat("Usage: bam2cov.R bamfile bigwig_file_name.bw cpu_number\n")
    quit(save="no")
} else {
    library(GenomicAlignments)
    library(rtracklayer)
}

bam <- BamFile(args[1])
bigwig_file <- args[2]
ncores <- args[3]

cov <- as(mclapply(seqlevels(bam),function(chr){
    params <- ScanBamParam(which=GRanges(chr,IRanges(1,seqlengths(bam)[chr])))
    aln <- readGAlignments(bam,param=params)
    coverage(aln)[[chr]]
},mc.preschedule=TRUE,mc.cores=ncores),"RleList")

names(cov) <- seqlevels(bam)

gr <- as(cov,"GRanges")

export(gr,bigwig_file)



