FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
	&& DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y \
	python3 \
	python3-pip \
	&& pip3 install \
	numpy \
	scipy \
	py2bit \
	pyBigWig \
	pysam \
	matplotlib \
	deeptools \
	&& apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*






