#!/usr/bin/env Rscript
library(GenomicAlignments)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

bam <- BamFile(args[1])
ncore <- args[2]

chrs <- seqnames(bam)

if (dir.exists(file.path(expDir,"genome"))) {
    chrs <- read.table(file.path(expDir,"genome","coolerFiles","chr_size.txt"))
} else {
    chrs <- read.table("~/ebs/ref_push/QCscripts/coolerDir/hg38.txt")
}

cov <- as(mclapply(seqlevels(bam),function(chr){

    params <- ScanBamParam(which=GRanges(chr,IRanges(1,seqlengths(bam)[chr])),
                           flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                            isDuplicate=FALSE),
                           mapqFilter=40)

    aln <- readGAlignments(bam,param=params)
    coverage(aln)[[chr]]
},mc.preschedule=TRUE,mc.cores=system("nproc",inter=TRUE)),"RleList")

names(cov) <- seqlevels(bam)

cov <- cov[chrs[,1]]
gr <- as(cov,"GRanges")


bw.dir <- file.path(dirname(dirname(path(bam))),'bigwigs')
bw.f <- file.path(bw.dir,sub("\\.bam","\\.bw",basename(path(bam))))
dir.create(bw.dir,showWarnings=FALSE)

export(gr,bw.f)



