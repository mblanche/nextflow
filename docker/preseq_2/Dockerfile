FROM mblanche/basespace-cli

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    build-essential \
    wget \
    zlib1g-dev \
    && wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz \
    && tar xvf preseq-3.1.2.tar.gz \
    && cd /preseq-3.1.2 \
    && ./configure \
    && make -j \
    && mv preseq /usr/local/bin \
    && rm -rf /preseq* \
    && apt-get --purge -y remove \
    build-essential \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*





