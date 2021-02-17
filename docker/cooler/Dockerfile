FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
	&& DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y \
	python3 \
	python3-pip \
	git \
	&& pip3 install cooler \
	&& git clone https://github.com/4dn-dcic/pairix \
	&& cd pairix \
	&& make -j \
	&& cp bin/* /usr/local/bin/ \
	&& rm -rf /pairix \
	&& apt clean \
	&& apt purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*






