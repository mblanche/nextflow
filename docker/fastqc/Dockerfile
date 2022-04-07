FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    default-jre-headless \
    wget \
    unzip \
    build-essential \
    cpanminus \
    python3 \
    python3-pip \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN cpanm FindBin

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && rm -rf fastqc_v0.11.9.zip

RUN chmod 755  /FastQC/fastqc \
    && ln -f -s /FastQC/fastqc /usr/local/bin/fastqc

RUN pip install multiqc
