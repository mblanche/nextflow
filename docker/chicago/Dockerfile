FROM mblanche/r-base

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    wget \
    python3 \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && mv bedtools.static.binary /usr/local/bin/bedtools \
    && chmod a+x /usr/local/bin/bedtools

RUN wget https://raw.githubusercontent.com/dovetail-genomics/chicago/main/chicagoTools/bam2chicago.sh \
    && chmod 755 bam2chicago.sh \
    && mv bam2chicago.sh /usr/local/bin/

RUN wget https://raw.githubusercontent.com/dovetail-genomics/chicago/main/chicagoTools/makeDesignFiles_py3.py

RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
	-e 'install.packages(c("R.utils","BiocManager"))' \
	-e 'BiocManager::install(version = "3.12")' \
	-e 'BiocManager::install(c("rtracklayer","Rsamtools"))'

COPY prep4Chicago.R /usr/bin/prep4Chicago

RUN chmod 755 /usr/bin/prep4Chicago


RUN /usr/bin/R -e 'options(Ncpus = parallel::detectCores())'\
    -e 'BiocManager::install(c("argparser","Chicago"))'

RUN wget https://raw.githubusercontent.com/dovetail-genomics/chicago/main/chicagoTools/runChicago.R \
    && chmod 755 runChicago.R
    && mv runChicago.R /usr/local/bin/runChicago
