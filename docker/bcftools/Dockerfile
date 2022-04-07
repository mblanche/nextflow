FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    bcftools

RUN apt-get --purge -y remove \
    && apt-get clean \
    && apt-get  purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*







