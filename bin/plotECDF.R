#!/usr/bin/env Rscript
library(GenomicAlignments)
library(ggplot2)
library(parallel)

args = commandArgs(trailingOnly=TRUE)

cores <- args[1]
prefix <- args[2]
bams <- BamFileList(args[-(1:2)])


allWidths <- lapply(bams,function(bam){
    cat(paste("working on",basename(bam),"\n"))
    
    ## Take the longest chromosome of the first BAM
    chr <- names(sort(seqlengths(bam),decreasing=TRUE))[1]
    chr.l <- sort(seqlengths(bam),decreasing=TRUE)[1]
    
    bin.size <- 1e5
    nbreaks <- ceiling(chr.l/bin.size)
    ends <- c(seq(from=bin.size+1,length.out=nbreaks-1,by=bin.size),as.integer(chr.l))
    widths <- do.call(c,mclapply(ends,function(end){
        params <- ScanBamParam(which=GRanges(chr,IRanges(end-bin.size,end)),
                               mapqFilter=40,
                               flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                                isDuplicate=FALSE))        
        aln <- readGAlignmentPairs(bam,param=params)
        aln.width <- width(unique(granges(aln,on.discordant.seqnames='drop')))
    },mc.preschedule=FALSE,mc.cores=cores))
    return(widths)
})

d.f <- data.frame(lib=rep(sub(".bam", "",names(bams)),sapply(allWidths,length)),
                  data=do.call(c,allWidths))

p <- ggplot(d.f,aes(data,color=lib))
p <- p+stat_ecdf()
p <- p+scale_x_log10()
p <- p+labs(title="Cummulative distribution of linkage",
            subtitle="Computed largest chromsome",
            x="Linkage Distance [log10(bp)]",
            y="Fraction of Reads")


out.file <- file.path(paste0(prefix,"_ecdf.pdf"))
pdf(out.file)
print(p)
dev.off()
