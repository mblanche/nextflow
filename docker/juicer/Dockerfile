FROM mblanche/awscli

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    default-jre-headless \
    wget \
    && wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar \
    && apt clean \
    && apt purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*



