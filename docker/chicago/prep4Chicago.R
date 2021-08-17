#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

##args <- c('baits_v1.0.bed.bgz',5000,'PC-UNI1946-cleanedUp.bam')

if (length(args) != 3){
    cat("Usage: prep4Chicago.R baitFile.bed.bgz 5000 myCleanedUp.bam")
    quit('no')
}

library(rtracklayer)
library(Rsamtools)

baits <- import(args[1])
bin.size <- as.integer(args[2])
bam.f <- BamFile(args[3])

## Make sure the bam file and the baits uses same seqlevel style (UCSC or Ensembl)
seqlevelsStyle(baits) <- seqlevelsStyle(bam.f)

## Resize the baits region to include shoulders of read coverage
## Merge overlaping resized baits
baits.r <- reduce(resize(baits,width(baits)+2*120,'center'))

## Bin the genomes in block of bin.size
bins <- do.call(c,lapply(seqlevels(bam.f),function(chr,s.i){
    bin.starts <- seq(by=bin.size,length = (seqlengths(s.i)[chr] %/% bin.size) +1)
    bins <- GRanges(chr,
                    IRanges(bin.starts,c(bin.starts[-1]-1,seqlengths(s.i)[chr])),
                    seqinfo=s.i)
},s.i=seqinfo(bam.f)))


## Break the bins overlaping with the baits
OE.raw <- disjoin(c(bins,baits.r))

## finds with broken bins overlaps with resized baits
ol <- findOverlaps(OE.raw,baits.r)

## merge the baits with bins not overlaping baits
OE <- c(OE.raw[-queryHits(ol)],
        baits.r)
OE <- sort(OE)

## Add a unique id to the bins
mcols(OE)$name <- seq_along(OE)

## Add an extra columns with the concatenated original bait ids
baits.ol <- findOverlaps(OE,baits)

ids <- tapply(baits[subjectHits(baits.ol)]$name,
                 queryHits(baits.ol),
                 paste,
                 collapse=",")

f <- as.integer(levels(factor(queryHits(baits.ol))))
mcols(OE)$id <- ''
mcols(OE)$id[f] <- ids

## Write the rmap and baitmap files
write.table(data.frame(seqnames(OE),start(OE),end(OE),OE$name),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE,
            sep="\t",
            file.path(paste0(bin.size/1000,"kb.rmap")))


br <- OE[OE %in% baits.r]

write.table(data.frame(seqnames(br),start(br),end(br),br$name,br$id),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE,
            sep="\t",
            file.path(paste0(bin.size/1000,"kb.baitmap")))






