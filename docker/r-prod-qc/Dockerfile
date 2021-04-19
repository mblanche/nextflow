FROM mblanche/r-base

USER root

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'install.packages(c("R.utils","BiocManager"))' \
	-e 'BiocManager::install(version = "3.12")' \
	-e 'BiocManager::install(c("GenomicAlignments","ShortRead"))'

COPY alnStats.R /usr/bin/alnStats
COPY adptStats.R /usr/bin/adptStats

RUN chmod 755 /usr/bin/alnStats /usr/bin/adptStats


