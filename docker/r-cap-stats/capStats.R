#!/usr/bin/env Rscript
args  <- commandArgs(trailingOnly=TRUE)

##args <- c('probes_v1.0.bed.bgz', 48 ,'PC-bGurd-4.bam', 'PC-bGurd-1.bam', 'PC-bGurd-3.bam', 'PC-bGurd-2.bam')

args <- c("probes_v1.0.bed.bgz", 25, "PC-Gurdon-WG-5.bam", "PC-Gurdon-WG-6.bam")

if (length(args) < 3) {
    cat("Usage: capStats probes.bed ncores lib1.bam lib2.bam lib3.bam ...\n")
    quit(save="no")
}

library(GenomicAlignments)
library(rtracklayer)
library(ggplot2)

probes <- import(args[1])
ncpus <- args[2]
bams <- args[-(1:2)]

################################################################################
## Resize the probes by 200 bp on both side
################################################################################
probes.rs <- resize(probes,width= width(probes) + 2*150,fix='center')
seqlevelsStyle(probes.rs) <- 'UCSC'

################################################################################
## process the bam file in parallel across all chromosome
## Couunt the number of reads overlaping the probes
################################################################################
chrs <- rep(seqlevels(probes.rs),length(bams))
bams.f <- rep(bams,length(seqlevels(probes.rs)))

res <- mclapply(seq_along(bams.f),function(i){

    bam <- BamFile(bams.f[[i]])
    chr <- chrs[i]
    chr.l <- seqlengths(bam)[chr]
    
    params <- ScanBamParam(flag=scanBamFlag(isPaired=TRUE,
                                            isSecondaryAlignment=FALSE,
                                            isSupplementaryAlignment=FALSE),
                           which=GRanges(chr,IRanges(1,chr.l)))
    ga <- readGAlignments(bam,param=params)
    o.l <- findOverlaps(ga,probes.rs)
    
    cov <- coverage(ga)[[chr]]
    r <- ranges(probes.rs[seqnames(probes.rs) == chr])
    
    return(
        list(readsInProbe = length(queryHits(o.l)),
             nReads = length(ga),
             ol.probes = countOverlaps(probes.rs,ga),
             meanProbeCovs = mean(unlist((Views(cov,r)))),
             meanChrCovs = mean(unlist(Views(cov,gaps(r,1,chr.l))))
             )
    )
},mc.preschedule=FALSE,mc.cores=ncpus)

################################################################################
### Massging the returns to compute some stats
################################################################################
reads <- do.call(rbind,tapply(res,bams.f,function(res){
    r <- sum(do.call(c,lapply(res,'[[','readsInProbe')))
    total <- sum(do.call(c,lapply(res,'[[','nReads')))
    return(c(nreads=total,
             onProbes=r,
             onProbesFrac=r/total))
}))

meanCovProbe <- tapply(res,bams.f,function(res){
    mean(do.call(c,lapply(res,'[[','meanProbeCovs')))
})

meanCovChr <- tapply(res,bams.f,function(res){
    mean(do.call(c,lapply(res,'[[','meanChrCovs')))
})

stats <- cbind(data.frame(lib=sub("\\.bam","",rownames(reads))),
    reads,
    meanCovProbe,
    meanCovChr)

stats$enrichment <- stats$meanCovProbe/stats$meanCovChr
    
write.csv(stats,row.names=FALSE,quote=FALSE,"capture_stat.csv")

################################################################################
### Ploting boxplot distributions of read count per probes
################################################################################
probes <- tapply(res,bams.f,function(res){
    p <- colSums(do.call(rbind,lapply(res,'[[','ol.probes')))
})

d.f <- data.frame(lib=rep(sub("\\.bam","",names(probes)),sapply(probes,length)),
                  counts=do.call(c,probes))

p <- ggplot(d.f,aes(lib,counts)) +
    geom_boxplot() +
    labs(title="Number of reads per probe in each library",
         x="Library",
         y="Counts of reads per probe")

pdf("DistReadProbe.pdf")
print(p)
dev.off()



