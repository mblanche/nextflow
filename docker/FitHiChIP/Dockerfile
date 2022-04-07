FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

# base dependencies
RUN apt-get update -y \
    && apt-get install -y \
    python3 \
    python3-pip \
    software-properties-common \
    dirmngr \
    build-essential \
    samtools \
    git \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    gnupg2 \
    wget

# python stuff
RUN pip install pysam bx-python numpy scipy py2bit pyBigWig matplotlib deeptools networkx

# r base
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update -y \
    && apt-get install -y r-base

# Install required R packages
RUN R -e "install.packages(c('BiocManager', 'optparse', 'ggplot2', 'data.table', 'splines', 'fdrtool', 'parallel', 'tools', 'plyr', 'dplyr', 'RColorBrewer'), quitely=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e 'BiocManager::install()'
RUN R -e "BiocManager::install('GenomicRanges');"
RUN R -e "BiocManager::install('edgeR');"

## Install miniconda.
RUN wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O ~/anaconda.sh
RUN bash ~/anaconda.sh -b -p /usr/local/anaconda
RUN rm ~/anaconda.sh
ENV PATH $PATH:/usr/local/anaconda/bin

## Install all HiC-Pro dependencies using conda
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /usr/local/anaconda/envs/HiC-Pro_v3.1.0/bin:$PATH

## Install HiCPro
RUN cd /tmp && \
    echo "master.zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && \
    unzip hicpro_latest.zip && \
    cd HiC-Pro-master  && \ 
    make configure prefix=/ && \
    make install && \
    cd .. && \
    rm -fr HiC-Pro*
    ENV PATH="/HiC-Pro_3.1.0/bin/:${PATH}"
    ENV PATH="/HiC-Pro_3.1.0/bin/utils/:${PATH}"

# Install bowtie2
RUN conda install -c bioconda bowtie2

# install MACS2
RUN conda install -c bioconda macs2

# install bedtools
RUN conda install -c bioconda bedtools

# Get FitHiChIP
RUN apt-get -y install git
#RUN cd / && git clone https://github.com/ay-lab/FitHiChIP
# this version has messages suppressed to stop nextflow from freaking out
RUN cd / && git clone https://github.com/torchij/FitHiChIP.git
RUN cd /FitHiChIP && sed -i 's/\/home\/sourya\/packages\/HiCPro\/HiC-Pro_2.9.0\//\/HiC-Pro-2.11.4\//g' configfile_*

# Cleanup
RUN rm -rf /*tar*
RUN chmod -R 777 /root
RUN conda install -y -c bioconda iced 
RUN conda clean --all -y
RUN apt -y install bc
