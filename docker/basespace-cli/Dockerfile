FROM ubuntu:20.04

RUN apt-get update -y \
	&& DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y \
	wget \
	gnupg2 \
	&& wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /usr/local/bin/bs \
	&& chmod 755 /usr/local/bin/bs \
	&& apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*






