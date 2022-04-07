
FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    git \
    wget \
    curl \
    build-essential

RUN  apt-get clean \
    && apt-get  purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz \
    && tar xvf mummer-4.0.0rc1.tar.gz \
    && cd mummer-4.0.0rc1 \
    && ./configure \
    && make -j \
    && make install \
    && ldconfig
