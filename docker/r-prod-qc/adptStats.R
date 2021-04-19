#!/usr/bin/env Rscript
library(ShortRead)

args = commandArgs(trailingOnly=TRUE)

expName <- args[1]
fqs <- args[-1]

## What fraction of reads have a ligated junction
stats <- as.data.frame(do.call(rbind,mclapply(fqs,function(fl){

    adpt <- DNAString('GGTTCGTCCATCGATC')
    
    m <- c(adpt,subseq(reverseComplement(adpt),5,length(adpt)))
    patt <- subseq(m,length(adpt)-7,width=12)
    half.adp.patt <- subseq(reverseComplement(adpt),5)

    half.adp.stat <- mm1 <- mm2 <- l <- 0
    
    f <- FastqStreamer(fl, 1e6)
    
    while (length(fq <- yield(f))) {
        r <- sread(fq)
        nreads <- length(r)
        mm <- vmatchPattern(patt,r)
        mm1 <- c(mm1,sum(elementNROWS(mm) > 0)*2/length(r))
        filt <- elementNROWS(mm) > 0
        start.with.adpt <- sapply(relist(start(unlist(mm[filt])) == length(adpt)-5,mm[filt]),any)
        if (length(start.with.adpt) > 0){
            mm2 <- c(mm2,sum(start.with.adpt)/length(start.with.adpt))
        }
        l <- l + length(r)
        mm.half.adp <- vmatchPattern(half.adp.patt,r,max.mismatch=2)
        filt.half.adp <- elementNROWS(mm.half.adp) > 0
        half.adp.cnt <- sum(start(unlist(mm.half.adp[filt])) == 1)
        half.adp.stat <- c(half.adp.stat,half.adp.cnt*2/length(r))
    }
    close(f)
    c(length=l,mm1=mean(mm1),mm2=mean(mm2),half.adp.stat=mean(half.adp.stat))

},mc.preschedule=FALSE,mc.cores=system('nproc',intern=TRUE))))

libs <- factor(sub("_S.+","",basename(fqs)))
dd <- data.frame(lib.id=levels(libs),
                 nReadPairs=tapply(stats$length,libs,sum),
                 fracWithHalfAdp=tapply(stats$half.adp.stat,libs,mean)*100,
                 fracWithLig=tapply(stats$mm1,libs,mean)*100,
                 startWithBridge=tapply(stats$mm2,libs,mean)*100)

dd <- dd[order(dd$lib.id),]


outFile <- paste0(expName,"_adptStats.csv")
write.csv(dd,outFile)

