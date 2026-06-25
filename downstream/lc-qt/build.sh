#!/bin/sh
set -e

echo "🐳 Building Docker compilation image..."
docker build -t lc-qt-builder .

echo "🏗️ Compiling lc-qt statically inside Alpine..."
docker run --rm -v "$(pwd)":/app -v "$(pwd)/../../include":/external_include lc-qt-builder bash -c "

    # Delete buildroot if it is a git repo, so we can replace with SDK
    if [ -d "buildroot" ] && [ -d "buildroot/.git" ]; then
        echo '🧹 Removing Git repository buildroot to replace with prebuilt SDK...'
        rm -rf buildroot
    fi

    if [ ! -d "buildroot" ]; then
        echo '📥 Downloading prebuilt Buildroot SDK...'
        wget -q https://github.com/fifbroman8/buildroot/releases/download/v1.1.0/x86_64-buildroot-linux-musl_sdk-buildroot.tar.gz
        tar xf x86_64-buildroot-linux-musl_sdk-buildroot.tar.gz
        mv x86_64-buildroot-linux-musl_sdk-buildroot buildroot
        rm -f x86_64-buildroot-linux-musl_sdk-buildroot.tar.gz
    fi
    
    echo '⚙️ Relocating toolchain SDK...'
    ./buildroot/relocate-sdk.sh
    
    echo '🧹 Disabling dynamic libraries in sysroot...'
    find buildroot/x86_64-buildroot-linux-musl/sysroot -name \"*.so\" -delete
    
    echo '🏗️ Compiling lc-qt statically using the Buildroot toolchain...'
    cd /app
    rm -rf build
    mkdir -p build && cd build
    cmake .. \
        -DCMAKE_TOOLCHAIN_FILE=/app/buildroot/share/buildroot/toolchainfile.cmake \
        -DEXTERNAL_INCLUDE_DIR=/external_include
    make -j\$(nproc)
"






echo "✅ Compilation finished. Checking file type of build/lc-qt:"
file build/lc-qt
