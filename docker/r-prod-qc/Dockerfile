FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && apt-get install -y \
    software-properties-common \
    dirmngr \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    gnupg2 \
    wget

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update -y \
    && apt-get install -y r-base

RUN apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'install.packages("BiocManager")' \
    -e 'BiocManager::install()'

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'install.packages("R.utils")'

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'BiocManager::install("GenomicAlignments")'

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'BiocManager::install("ShortRead")'

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'BiocManager::install("ggplot2")'

COPY alnStats.R /usr/bin/alnStats
COPY adptStats.R /usr/bin/adptStats

RUN chmod 755 /usr/bin/alnStats /usr/bin/adptStats


