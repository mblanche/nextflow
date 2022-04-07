FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    python3 \
    python3-pip \
    && pip3 install \
    argh \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

COPY n50.py /usr/local/bin

RUN chmod 755 /usr/local/bin




