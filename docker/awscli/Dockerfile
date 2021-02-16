FROM ubuntu:20.04

RUN apt update  \
    && apt install -y \
    unzip \
    curl \
    && curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && apt clean \
    && apt purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* ./aws




