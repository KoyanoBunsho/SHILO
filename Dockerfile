FROM ubuntu:24.04

RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        emacs \
        git \
        libeigen3-dev \
        libssl-dev \
        sudo \
        unzip \
        wget \
        zlib1g && \
    rm -rf /var/lib/apt/lists/*

# SET path to compilers.
# https://stackoverflow.com/questions/17275348/how-to-specify-new-gcc-path-for-cmake
ENV CC=/usr/bin/gcc \
    CXX=/usr/bin/g++

# Please use below directory to install cpp libraries.
WORKDIR /usr/src

CMD ["bash"]
