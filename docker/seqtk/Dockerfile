
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
    bedtools \
    seqtk \
    && git clone https://github.com/lh3/bwa.git \
    && cd /bwa \
    && make -j \
    && mv bwa /usr/local/bin \
    && rm -rf /bwa \
    && apt-get --purge -y remove \
    && apt-get clean \
    && apt-get  purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*





