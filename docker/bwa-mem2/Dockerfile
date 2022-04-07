FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    build-essential \
    samtools \
    git \
    zlib1g-dev \
    seqtk \
    wget

RUN apt-get --purge -y remove \
    && apt-get clean \
    && apt-get  purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*


RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 -O  - \
    |tar -xj 

ENV PATH="/bwa-mem2-2.2.1_x64-linux:${PATH}"

RUN git clone git://github.com/GregoryFaust/samblaster.git \
    && cd samblaster \
    && make -j \
    && mv samblaster /usr/local/bin/. \
    && rm -rf ../samblaster

COPY splitUG /usr/bin/splitUG
