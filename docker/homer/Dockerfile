FROM mblanche/r-base

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    build-essential \
    wget \
    zip \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*


RUN mkdir - p /homer \
    && cd /homer \ 
    && wget http://homer.ucsd.edu/homer/configureHomer.pl \
    && perl configureHomer.pl -install

ENV PATH="/homer/bin:${PATH}"

RUN wget -P /usr/local/bin  http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
    && chmod 755 /usr/local/bin/bedGraphToBigWig
