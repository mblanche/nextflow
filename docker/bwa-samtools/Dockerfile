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
	build-essential \
	&& git clone https://github.com/lh3/bwa.git \
	&& cd /bwa \
	&& make -j \
	&& mv bwa /usr/local/bin \
	&& rm -rf /bwa \
	&& apt --purge -y remove \
	build-essential \
	git \
	&& apt clean \
	&& apt purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*





