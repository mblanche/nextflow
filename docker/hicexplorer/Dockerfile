FROM ubuntu:20.04
USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y\
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    python3 \
    python3-pip \
    imagemagick

RUN sed -i '/disable ghostscript/d' /etc/ImageMagick-6/policy.xml \
    && sed -i '/pattern="PS"/d' /etc/ImageMagick-6/policy.xml \
    && sed -i '/pattern="PS2"/d' /etc/ImageMagick-6/policy.xml \
    && sed -i '/pattern="PS3"/d' /etc/ImageMagick-6/policy.xml \
    && sed -i '/pattern="EPS"/d' /etc/ImageMagick-6/policy.xml \
    && sed -i '/pattern="PDF"/d' /etc/ImageMagick-6/policy.xml \
    && sed -i '/pattern=XPS"/d' /etc/ImageMagick-6/policy.xml

RUN  apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*


RUN pip3 install hicexplorer








