FROM mblanche/awscli

RUN apt update  \
    && apt install -y \
    wget \
    gnupg2 \
    && wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O /bin/bs \
    && chmod 755 /bin/bs \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV BASESPACE_API_SERVER="https://api.basespace.illumina.com"
ENV BASESPACE_ACCESS_TOKEN="dbd2961a36a74bae9f207f4d82dea3f2"




