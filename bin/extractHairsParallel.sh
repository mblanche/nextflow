#!/usr/bin/env bash

bam=$1; shift
vcf=$1; shift

for chrInfo in $(samtools view -H ${bam} \
		     |awk -v OFS='\t' '/^@SQ/ && !($2 ~ /:(chr|"")M/) {split($2,chr,":");split($3,ln,":");print chr[2]":1-"ln[2]}' \
		     |sort -V -k1,1); do

    fname=$(basename $bam)
    chr=$(echo $chrInfo|cut -d: -f1)
    
    extractHAIRS \
	--bam ${bam} \
	--hic 1 \
	--VCF $vcf \
	--region $chrInfo \
	--out ${fname%.bam}_${chr}.hairs &
    
done

wait
