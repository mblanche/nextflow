FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
	&& apt-get install -y \
	build-essential \
	libcurl4-openssl-dev \
	libxml2-dev \
	gnupg2 

RUN echo deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/ | tee -a /etc/apt/sources.list \
	&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
	&& apt-get update -y \
	&& apt-get install -y r-base

RUN apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*
