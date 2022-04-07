FROM ubuntu:20.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && apt-get install -y \
    software-properties-common \
    dirmngr \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    gnupg2 \
    wget \
    samtools \
    git \
    zlib1g-dev \
    python3 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3 10

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update -y \
    && apt-get install -y r-base



RUN apt-get clean \
	&& apt-get purge \
	&& rm -rf /var/lib/apt/lists/* /tmp/*

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'install.packages("BiocManager")' \
    -e 'BiocManager::install()'

RUN git clone https://github.com/ISUgenomics/SequelTools.git \
    && cd SequelTools/Scripts \
    && chmod +x *.sh *.py *.R 

ENV PATH="${PATH}:/SequelTools/Scripts"
