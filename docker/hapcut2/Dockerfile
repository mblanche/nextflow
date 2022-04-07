FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    build-essential \
    libhts-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    git 

RUN apt-get --purge -y remove \
    && apt-get clean \
    && apt-get  purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*


RUN git clone https://github.com/vibansal/HapCUT2.git \
    && cd HapCUT2 \
    && make -j \
    && make install \
    && rm -rf ../HapCut2






