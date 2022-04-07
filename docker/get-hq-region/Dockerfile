FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
	&& DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y \
	wget \
	python3 \
	python3-pip \
	tabix \
	&& pip3 install \
	numpy \
	scipy \
	py2bit \
	argparse \
	pyBigWig \
	pysam \
	matplotlib \
	deeptools \
	&& apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && mv bedtools.static.binary /usr/local/bin/bedtools \
    && chmod a+x /usr/local/bin/bedtools

COPY get_HQ_region_bed.py /usr/bin/

RUN chmod 755 /usr/bin/get_HQ_region_bed.py
