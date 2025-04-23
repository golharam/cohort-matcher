FROM ubuntu:20.04

# Metadata
LABEL maintainer="NGS Bioinformatics"
LABEL description="Dockerfile for freebayes 1.2.0 - a Bayesian genetic variant detector"
LABEL version="1.2.0"

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    cmake \
    pkg-config \
    python3 \
    python3-pip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    autoconf \
    automake \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /opt

# Clone freebayes repository and checkout version 1.2.0
RUN git clone --recursive https://github.com/ekg/freebayes.git \
    && cd freebayes \
    && git checkout v1.2.0 \
    && git submodule update --recursive

# Build freebayes
RUN cd /opt/freebayes \
    && make \
    && make install

# Add freebayes to PATH
ENV PATH="/opt/freebayes/bin:${PATH}"

# Create a directory for data
RUN mkdir -p /data

# Set working directory to /data
WORKDIR /data

# Default command
ENTRYPOINT ["freebayes"]
CMD ["--help"]