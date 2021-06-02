FROM ubuntu:21.04

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
    build-essential \
    wget \
    zip \
    git \
    libssl-dev\
    libjsoncpp-dev \
    zlib1g-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz \
    && tar -zxvf cmake-3.20.0.tar.gz \
    && cd cmake-3.20.0 \
    && ./bootstrap --parallel=$(nproc) \
    && make -j \
    && make install

RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2 \
    && tar xvf eigen-3.3.9.tar.bz2 \
    && cd eigen-3.3.9 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make install

RUN git clone git://github.com/pezmaster31/bamtools.git \
    && cd bamtools \
    && mkdir build \
    && cd build \
    ##&& cmake -DBUILD_SHARED_LIBS=ON .. \
    && cmake .. \
    && make install


RUN git clone https://github.com/dixonlab/hic_breakfinder.git \
    && cd hic_breakfinder \
    && ./configure CPPFLAGS="-I/usr/local/include/bamtools/  -I/usr/local/include/eigen3/" LDFLAGS="-L/usr/local/lib" \
    && make -j 
    ##make install
