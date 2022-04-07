FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
	&& DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y \
	python3 \
	python3-pip \
	git \
	samtools \
	tabix \
	libbz2-dev \
	lz4 \
	zlib1g-dev \
	build-essential \
	autotools-dev \
	automake \
	&& pip3 install \
	numpy \
	&& pip3 install click \
	cython \
	nose \
	&& pip3 install \
	pairtools \
	&& git clone https://github.com/nh13/pbgzip.git \
	&& cd pbgzip \
	&& sh autogen.sh \
	&& ./configure \
	&& make -j \
	&& make install \
	&& cd / \
	&& rm -rf /pbgzip \
	&& git clone https://github.com/4dn-dcic/pairix \
	&& cd pairix \
	&& make -j \
	&& cp bin/* /usr/local/bin/ \
	&& cd / \
	&& rm -rf /pairix \
	&& apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*

RUN git clone https://github.com/4dn-dcic/bam2pairs \
	&& mv /pairix/util/bam2pairs/bam2pairs /usr/bin/ \
	&& rm -rf pairix