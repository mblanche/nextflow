FROM mblanche/r-base

USER root

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'install.packages("BiocManager")' \
	-e 'BiocManager::install(version = "3.12")' \
	-e 'BiocManager::install(c("GenomicAlignments","rtracklayer"))'

