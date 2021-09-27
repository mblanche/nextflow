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
    python3 \
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
    -e 'BiocManager::install("rtracklayer")'

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'install.packages("Rsamtools")'

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && mv bedtools.static.binary /usr/local/bin/bedtools \
    && chmod a+x /usr/local/bin/bedtools

RUN wget https://raw.githubusercontent.com/dovetail-genomics/chicago/main/chicagoTools/bam2chicago.sh \
    && chmod 755 bam2chicago.sh \
    && mv bam2chicago.sh /usr/local/bin/

RUN wget https://raw.githubusercontent.com/dovetail-genomics/chicago/main/chicagoTools/makeDesignFiles_py3.py

COPY prep4Chicago.R /usr/bin/prep4Chicago

RUN chmod 755 /usr/bin/prep4Chicago

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'BiocManager::install(c("argparser","Chicago"))'

RUN wget https://raw.githubusercontent.com/dovetail-genomics/chicago/main/chicagoTools/runChicago.R \
    && chmod 755 runChicago.R \
    && mv runChicago.R /usr/local/bin/runChicago
