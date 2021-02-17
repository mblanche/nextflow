FROM mblanche/awscli

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt update \
    && apt install \
    build-essential \
    samtools \
    git \
    zlib1g-dev \
    build-essential \
    -y --no-install-recommends \
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





