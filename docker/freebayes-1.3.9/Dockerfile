FROM ubuntu:24.04

MAINTAINER Ryan Golhar <ryan.golhar@bms.com>

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    bc \
    samtools \
    vcftools \
    parallel \
    git \
    cmake \
    pkg-config \
    libseqlib-dev \
    libvcflib-tools \
    libvcflib-dev \
    libhts3 libhts-dev \
    libsmithwaterman-dev \
    libfastahack-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    liblzma-dev \
    meson \
    ninja-build \
    libwfa2-0 libwfa2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /opt

# Clone freebayes repository and checkout version 1.3.9
RUN git clone --recursive https://github.com/ekg/freebayes.git \
    && cd freebayes \
    && git checkout v1.3.9

# Build freebayes
WORKDIR /opt/freebayes
RUN meson build/ --buildtype release && \
    cd build && \
    ninja && \
    ninja test && \
    cp freebayes /usr/local/bin/freebayes

# Default command
ENTRYPOINT ["freebayes"]
CMD ["--help"]