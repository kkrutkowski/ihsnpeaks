#!/bin/sh
set -e

echo "🐳 Building Docker compilation image..."
docker build -t lc-qt-builder .

echo "🏗️ Compiling lc-qt statically inside Alpine..."
docker run --rm -v "$(pwd)":/app lc-qt-builder bash -c "
    if [ ! -d \"x86_64-buildroot-linux-musl_sdk-buildroot\" ]; then
        echo '📥 Downloading Buildroot SDK...'
        wget -q https://github.com/fifbroman8/buildroot/releases/download/v1.1.0/x86_64-buildroot-linux-musl_sdk-buildroot.tar.gz
        tar xf x86_64-buildroot-linux-musl_sdk-buildroot.tar.gz
        rm -f x86_64-buildroot-linux-musl_sdk-buildroot.tar.gz
    fi
    
    echo '⚙️ Relocating toolchain SDK...'
    ./x86_64-buildroot-linux-musl_sdk-buildroot/relocate-sdk.sh
    
    echo '🧹 Disabling dynamic libraries in sysroot...'
    find x86_64-buildroot-linux-musl_sdk-buildroot -type d -name \"sysroot\" -exec sh -c 'for d; do rm -f \"\$d\"/usr/lib/*.so; done' _ {} +
    
    echo '🏗️ Building static executable...'
    rm -rf build
    mkdir -p build && cd build
    cmake .. \
        -DCMAKE_EXE_LINKER_FLAGS='-static' \
        -DCMAKE_TOOLCHAIN_FILE=/app/x86_64-buildroot-linux-musl_sdk-buildroot/share/buildroot/toolchainfile.cmake
    make -j\$(nproc)
"





echo "✅ Compilation finished. Checking file type of build/lc-qt:"
file build/lc-qt
