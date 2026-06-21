FROM debian:stable-slim

# Install build dependencies for Buildroot SDK and UPX
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    file \
    upx-ucl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
