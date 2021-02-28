FROM mblanche/r-base

USER root

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'install.packages("BiocManager")' \
	-e 'BiocManager::install(version = "3.12")' \
	-e 'BiocManager::install(c("GenomicAlignments","rtracklayer"))'

COPY bam2bw.R /usr/bin/bam2bw

RUN chmod 755 /usr/bin/bam2bw


