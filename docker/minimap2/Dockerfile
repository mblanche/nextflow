
FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    git \
    curl \
    wget \
    zlib1g-dev \
    software-properties-common \
    dirmngr \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    gnupg2

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    libssl-dev \
    pandoc

RUN apt-get --purge -y remove \
    && apt-get clean \
    && apt-get  purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update -y \
    && apt-get install -y r-base


RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'install.packages(c("devtools","optparse", "ggplot2", "plotly"))'


RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    | tar -jxvf - \
    && mv /minimap2-2.24_x64-linux/k8 /minimap2-2.24_x64-linux/paftools.js minimap2-2.24_x64-linux/minimap2 /usr/bin 


RUN git clone https://github.com/tpoorten/dotPlotly.git \
    && cp /dotPlotly/pafCoordsDotPlotly.R /usr/bin
