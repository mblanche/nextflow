FROM nvidia/cuda:8.0-devel

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    wget \
    default-jre-headless \
    && wget -O /juicer_tools.jar https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar


## java -Xms512m -Xmx2048m -jar /juicer_tools.jar hiccups -m 500 -r 5000,10000 --threads 0 -k KR --ignore-sparsity /mnt/DHAM001-A.hic test.loops
