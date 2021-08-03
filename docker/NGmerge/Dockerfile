FROM ubuntu:21.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    build-essential \
    git \
    zlib1g-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN git clone https://github.com/jsh58/NGmerge.git \
    && cd NGmerge \
    && make -j \
    && mv NGmerge /usr/local/bin


RUN git clone https://github.com/jsh58/Genrich.git \
    && cd Genrich \
    && make -j \
    && mv Genrich /usr/local/bin
