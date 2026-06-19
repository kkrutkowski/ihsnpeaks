FROM debian:stable-slim

# Install building tools and dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    make \
    wget \
    tar \
    bash \
    file \
    upx-ucl \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /app

