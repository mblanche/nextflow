#!/usr/bin/env Rscript
library(GenomicAlignments)
library(rtracklayer)
library(ggplot2)
library(ggpmisc)


args = commandArgs(trailingOnly=TRUE)

base.dir <- args[1]
prot <- args[2]
prefix <- args[3]
libs <- args[-(1:3)]

prot2encode <- c(CTCF="ENCFF017XLW",
                 H3K27ac="ENCFF361XMX",
                 H3K4me3="‘ENCFF859HPC",
                 YY1="ENCFF223MUF",
                 PolII="ENCFF455ZLJ",
                 RAD21="ENCFF618UTS",
                 SMC3="ENCFF599ZQS",
                 mCTCF="ENCFF187ZKI",
                 mH3K27Ac="ENCFF135ATY")


type <- c(CTCF="TF",
          H3K27ac="histone",
          H3K4me3="histone",
          YY1="TF",
          PolII="TF",
          RAD21="TF",
          SMC3="TF",
          mCTCF="TF",
          mH3K27Ac="histone")


dir <- file.path("~/ebs/ref_push",base.dir)

bams <- BamFileList(sapply(libs,function(x) list.files(dir,paste0(x,"\\.bam$"),recursive=TRUE,full=TRUE)))

exp.dir <- dirname(dirname(path(bams[[1]])))

encodeExp <- prot2encode[prot]

encode.f <- file.path('~/ebs/ref_push/QCscripts/encodeFiles',paste0(encodeExp,'.bed.gz'))
cols <- c('signalValue'="numeric",
              'pValue'='numeric',
          'qValue' = 'numeric',
          'peak'='integer')

encode <- import(encode.f,extraCols=cols)

width <- 2e3

ROI <- resize(resize(shift(encode,mcols(encode)$peak),1),width,'center')
mcols(ROI)$name <- paste0("id",1:length(ROI))
qs <- quantile(ROI$signalValue,seq(0.2,.8,0.2))
ROI <- ROI[mcols(ROI)$signalValue > qs[4] ]


stat <- lapply(bams,function(bam){
    exp <- sub("\\.bam","",basename(path(bam)))
    
    cat(paste("Working on",exp,"\n"))
    
    total.read.flags <- scanBamFlag(isFirstMateRead=TRUE,
                                    isSecondaryAlignment=FALSE,
                                    isSupplementaryAlignment=FALSE)

    total.readPairs <- countBam(bam,param=ScanBamParam(flag=total.read.flags))$records
    
    res <- mclapply(seqlevels(bam),function(chr){
        
        params <- ScanBamParam(which=GRanges(chr,IRanges(1,seqlengths(bam)[chr])))
        aln <- readGAlignments(bam,param=params)
        
        aln.gr <- granges(aln)
        ## Shift the reads 73 bp so that they are over the nucleosome diad
        if (type[prot] == 'histone') {
            aln.gr[strand(aln.gr) == '+'] <- shift(aln.gr[strand(aln.gr) == '+'],73)
            aln.gr[strand(aln.gr) == '-'] <- shift(aln.gr[strand(aln.gr) == '-'],-73)
        }
        
        cov <- coverage(aln.gr)
        
        sub.roi <- ROI[seqnames(ROI) == chr]
        sub.roi <- sub.roi[start(sub.roi) >= 1 & end(sub.roi) <= seqlengths(bam)[chr]]
        sub.cov <- Views(cov[[chr]],ranges(sub.roi))

        sub.roi.peaks <- resize(sub.roi,400,'center')
        read.in.peaks <- suppressWarnings( sum(countOverlaps(sub.roi.peaks,aln)) )
        
        if (length(sub.cov) == 0){
            return(list(mat=vector(),
                        count=0))
        } else {
            sub.mat <- list(mat=as(sub.cov,"matrix"),
                            count=read.in.peaks)
        }
    },mc.preschedule=FALSE,mc.cores=detectCores())

    cov.mat <- do.call(rbind,lapply(res,'[[','mat'))
    
    rel.cov.mat <- cov.mat/rowSums(cov.mat)
    
    rel.cov <- colMeans(rel.cov.mat,na.rm=TRUE)*width

    count.peak <- sum(sapply(res,'[[','count'))

    obs.exp <- (count.peak/total.readPairs)/(sum(width(resize(ROI,400,'center')))/sum(seqlengths(bam)))


    return(list(rel.cov=rel.cov,obs.exp=obs.exp))
})


dd <- lapply(stat,'[[','rel.cov')

d.f <- data.frame(cov = unlist(dd),
                  x = rep(seq(width),length(dd)),
                  lib = rep(names(dd),sapply(dd,length))
                  )


d.t <- data.frame(lib=names(stat),
                 'obs/exp'=sprintf("%.1f",sapply(stat,'[[','obs.exp')))

rownames(d.t) <- NULL

p <- ggplot(d.f,aes(x,cov,color=lib))+geom_line()

p <- p+scale_x_continuous(breaks=c(0,width/2,width),labels=c(-width/2,"peak",width/2))
p <- p+labs(title=paste("Relative read coverage around \nhigh occupancy",prot,"binding sites"),
            x=paste("Base pair around",prot,"binding sites"),
            y="Relative Read Coverage")


ylim <- layer_scales(p)$y$range$range

p <- p + annotate(geom = "table",
                  x = 25,
                  y = ylim[2],
                  label = list(d.t), 
                  vjust = 1,
                  hjust = 0)


dir.out <- file.path(exp.dir,'QCplots')

dir.create(dir.out,showWarnings=FALSE)

f.name <- file.path(dir.out,paste0(prefix,"_",prot,"_profile.pdf"))


pdf(f.name)
print(p)
dev.off()
cat(paste("Done with",prot,"\n"))

