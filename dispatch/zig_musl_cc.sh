#!/bin/sh
set -eu

tool=$(basename "$0")
case "$tool" in
    aarch64-*|arm64-*)
        target=aarch64-linux-musl
        ;;
    x86_64-*)
        target=x86_64-linux-musl
        ;;
    *)
        target=${IHSNPEAKS_ZIG_TARGET:-}
        if [ -z "$target" ]; then
            echo "zig_musl_cc.sh: cannot infer target from '$tool'" >&2
            exit 2
        fi
        ;;
esac

set -- "$@"
translated="-Dconstexpr=const"
for arg do
    case "$target:$arg" in
        aarch64-linux-musl:-march=armv8-a)
            ;;
        aarch64-linux-musl:-march=armv8-a+simd)
            ;;
        aarch64-linux-musl:-march=armv8.2-a+sve)
            translated="$translated -mcpu=generic+sve"
            ;;
        aarch64-linux-musl:-march=armv8.5-a+sve2)
            translated="$translated -mcpu=generic+sve2"
            ;;
        *)
            quoted=$(printf "%s\n" "$arg" | sed "s/'/'\\\\''/g")
            translated="$translated '$quoted'"
            ;;
    esac
done

eval "exec zig cc -target '$target' $translated"
