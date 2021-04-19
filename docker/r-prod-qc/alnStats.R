#!/usr/bin/env Rscript
library(GenomicAlignments)
library(R.utils)

args  <- commandArgs(trailingOnly=TRUE)

exp.name <- args[1]
NCORES <- as.numeric(args[2])
cmplx_stats <- unlist(strsplit(args[3],','))
bams <- BamFileList(unlist(strsplit(args[4],',')))

bams <- bams[order(sub("\\.bam$","",path(bams)))]

################################################################################
### From the bam file, geting the number of reads ussed in the aligment
### ie: number of R1 reads...
################################################################################
readPair.count <- as.data.frame(do.call(rbind, mclapply(bams,function(bam){
    total.read.flags <- scanBamFlag(isFirstMateRead=TRUE,
                                    isSecondaryAlignment=FALSE,
                                    isSupplementaryAlignment=FALSE)
    aligned.flags <- scanBamFlag(isPaired=TRUE,
                                 isFirstMateRead=TRUE,
                                 isUnmappedQuery=FALSE,
                                 hasUnmappedMate =FALSE,,
                                 isSecondaryAlignment=FALSE,
                                 isSupplementaryAlignment=FALSE)
    counts <- do.call(c,mclapply(c(total.readPairs=1,total.alnPairs=2),function(i){
        if(i==1){
            total.readPairs <- countBam(bam,param=ScanBamParam(flag=total.read.flags))$records
        } else if (i==2){
            total.alnPairs <- countBam(bam,param=ScanBamParam(flag=aligned.flags))$records
        }
    },mc.preschedule=FALSE,mc.cores=2))
},mc.preschedule=FALSE,mc.cores=floor(NCORES/2))))

readPair.count$f.aligned <- readPair.count$total.alnPairs/readPair.count$total.readPairs
    
################################################################################
### Running a loop over the individual bam files
################################################################################
aln.stats <- as.data.frame(do.call(rbind,lapply(bams,function(bam){
    lib.name <- sub(".bam","",basename(bam))
    cat(paste("Working on",lib.name,"\n"))
    
################################################################################
### Grabbing the complextity at a given read depth (from precomputed preSeq)
################################################################################
    cmplxStat.f <- grep(paste0('^',lib.name,"_"),cmplx_stats,value=TRUE)
    seq.depth <- 4e8
    if(file.exists(cmplxStat.f) & countLines(cmplxStat.f) > 1){
        cmplx.s <- read.table(cmplxStat.f,header=TRUE)
        preSeq <- cmplx.s$EXPECTED_DISTINCT[cmplx.s$TOTAL_READS == seq.depth]
    } else {
        preSeq=NA
    }
    
################################################################################
### To make it faster, computing metrics only on longest Scafold/Chromosome
################################################################################
    f <- order(seqlengths(bam),decreasing=TRUE)[1]
    chrs <- seqnames(seqinfo(bam))[f]

################################################################################
### creating the mapping for the parallel step
################################################################################
    steps <- do.call(rbind,lapply(chrs,function(chr){
        bin.size <- 1e6
        chr.l <- seqlengths(bam)[chr]
        nbreaks <- ceiling(chr.l/bin.size)
        starts <- seq(length.out=nbreaks,by=bin.size)
        ends <- c(starts[-1],as.integer(chr.l))
        dd <- data.frame(chr=rep(chr,nbreaks),
                         start=starts,
                         end=ends)
    }))

################################################################################
### Computing the aligment stats in parallel, doing the reduce step in line
################################################################################
    aln.stats.file <- colMeans(do.call(rbind,mclapply(seq(nrow(steps)),function(i){
        d <- steps[i,]
        
        params <- ScanBamParam(which=GRanges(d$chr,IRanges(d$start,d$end)),
                               flag=scanBamFlag(isSecondaryAlignment=FALSE),
                               what='mapq',
                               mapqFilter=0)
        aln <- readGAlignmentPairs(bam,param=params)

        mq.gt.10 <- mcols(first(aln))$mapq >=10 & mcols(second(aln))$mapq >=10
        mq.gt.40 <- mcols(first(aln))$mapq >=40 & mcols(second(aln))$mapq >=40
        mq.gt.50 <- mcols(first(aln))$mapq >=50 & mcols(second(aln))$mapq >=50
        mq.eq.60 <- mcols(first(aln))$mapq ==60 & mcols(second(aln))$mapq ==60
        
        ## working on aln with mapq>40
        aln <- aln[mq.gt.40]

        ## Marking dups and removing them
        dups <- duplicated(paste(seqnames(first(aln)),start(first(aln)),seqnames(second(aln)),end(second(aln)),sep='-'))
        aln <- aln[!dups]

        ## flaging trans reads
        aln.grl <- grglist(aln)
        seqlevels(aln) <- as.character(unique(seqnames(unlist(aln.grl))))
        chrs.sub <- as.character(seqnames(unlist(aln.grl)))
        idx <- rep(seq(aln),elementNROWS(aln.grl))
        discordant <- tapply(chrs.sub,idx,function(x) length(unique(x)) != 1)

        ## Converting cis read to granges
        aln.gr <- granges(aln, on.discordant.seqnames="drop",use.names=FALSE)

        aln.count <- length(aln)
        
        return(c(f.mq.gt.10=sum(mq.gt.10)/length(mq.gt.10),
                 f.mq.gt.40=sum(mq.gt.40)/length(mq.gt.40),
                 f.mq.gt.50=sum(mq.gt.50)/length(mq.gt.50),
                 f.mq.eq.60=sum(mq.eq.60)/length(mq.eq.60),
                 f.dups=sum(dups)/length(dups),
                 f.inter.chr = sum(discordant)/aln.count,
                 f.1to350 =     sum(width(aln.gr) <= 350)/aln.count,
                 f.350.1kb =    sum(width(aln.gr) > 350 & width(aln.gr) <= 1e3)/aln.count,
                 f.1kb.10kb =   sum(width(aln.gr) > 1e3 & width(aln.gr) <= 1e4)/aln.count,
                 f.10kb.100kb = sum(width(aln.gr) > 1e4 & width(aln.gr) <= 1e5)/aln.count,
                 f.100kb.1mb =  sum(width(aln.gr) > 1e5 & width(aln.gr) <= 1e6)/aln.count,
                 f.gt1mb =      sum(width(aln.gr) > 1e6)/aln.count
              ))
    },mc.preschedule=TRUE,mc.cores=NCORES)),na.rm=TRUE)

    res <- c(aln.stats.file,preSeq)
    names(res)[length(res)] <- paste0('complexity@',seq.depth/1e6,'M')
    return(res)
})))

res <- cbind(readPair.count,aln.stats)

cat("Writing the Report")
report <- file.path(paste0(exp.name,".csv"))
write.csv(res,report)
