"""Build hwloc static library from source for cross-compilation release targets.

Shared by all dispatch/build_release*.py scripts. Downloads and extracts the
hwloc tarball once per build_root, then runs `./configure && make && make install`
using the given target cross-compiler. Returns paths to the installed headers
and static library.

The hwloc configure script's compiler-sanity-check sometimes fails under exotic
cross-compilers (e.g. zig cc wrapper scripts). We work around this by writing a
trivial `CC` wrapper that pretends to behave like a normal gcc when configure
probes it, while still passing the real cross-compilation flags through.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

HWLOC_VERSION = "2.14.0"
HWLOC_TARBALL_URL = f"https://download.open-mpi.org/release/hwloc/v2.14/hwloc-{HWLOC_VERSION}.tar.gz"


def _download_hwloc_tarball(tarball: Path) -> None:
    tarball.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading hwloc {HWLOC_VERSION} -> {tarball}", flush=True)
    urllib.request.urlretrieve(HWLOC_TARBALL_URL, tarball)


def _extract_hwloc(tarball: Path, build_root: Path) -> Path:
    src_dir = build_root / f"hwloc-{HWLOC_VERSION}"
    if src_dir.exists():
        return src_dir
    print(f"Extracting {tarball} -> {build_root}", flush=True)
    with tarfile.open(tarball) as tf:
        tf.extractall(build_root, filter="data")
    if not src_dir.exists():
        raise SystemExit(f"Expected hwloc source dir at {src_dir} after extraction")
    return src_dir


def ensure_hwloc(
    root: Path,
    build_root: Path,
    cc: list[str] | str,
    cflags: list[str] | None = None,
    ldflags: list[str] | None = None,
    host: str | None = None,
    stamp_name: str = ".hwloc_built",
) -> tuple[Path, Path]:
    """Build hwloc static library for the given target cross-compiler.

    Parameters
    ----------
    root : Path
        Repository root (used as working directory when running configure).
    build_root : Path
        Per-variant build directory. Contains hwloc-source, hwloc-build, and
        the installed prefix.
    cc : list[str] | str
        Target cross-compiler. Either a command list (["zig", "cc", "-target", ...])
        or a single tool name ("x86_64-alpine-linux-musl-gcc").
    cflags : list[str], optional
        Extra CFLAGS to pass through to configure.
    ldflags : list[str], optional
        Extra LDFLAGS to pass through to configure.
    host : str, optional
        Target triple to pass to configure's --host flag (e.g.,
        "x86_64-apple-darwin" or "aarch64-linux-musl"). Enables cross-compilation
        mode so configure skips runtime tests.
    stamp_name : str
        Name of stamp file inside build_root. Per-arch builds should use
        distinct stamp names so they don't collide.

    Returns
    -------
    (include_dir, libhwloc_path)
    """
    stamp = build_root / stamp_name
    include_dir = build_root / "hwloc-install" / "include"
    lib_dir = build_root / "hwloc-install" / "lib"
    libhwloc = lib_dir / "libhwloc.a"

    if stamp.exists() and libhwloc.exists() and include_dir.exists():
        return include_dir, libhwloc

    tarball = build_root / f"hwloc-{HWLOC_VERSION}.tar.gz"
    if not tarball.exists():
        _download_hwloc_tarball(tarball)
    src_dir = _extract_hwloc(tarball, build_root)

    hwloc_build = build_root / "hwloc-build"
    hwloc_prefix = build_root / "hwloc-install"
    hwloc_build.mkdir(parents=True, exist_ok=True)
    hwloc_prefix.mkdir(parents=True, exist_ok=True)

    # Build the CC wrapper. Configure probes the compiler with things like
    # `cc -E` (cpp mode) and `a.out` execution tests. We wrap the real CC so
    # that: (a) the wrapper itself is a shell script executable, (b) it
    # forwards all args to the real cross-compiler, (c) we set CFLAGS/
    # LDFLAGS explicitly so configure's flag probes see them.
    cc_cmd = cc if isinstance(cc, list) else [cc]
    wrapper_path = build_root / "cc-wrapper.sh"
    quoted_args = " ".join(f"'{arg}'" if " " in arg or "'" in arg else arg for arg in cc_cmd)
    wrapper_path.write_text(
        f"#!/bin/sh\nexec {quoted_args} \"$@\"\n",
        encoding="utf-8",
    )
    wrapper_path.chmod(0o755)

    cflags_env = " ".join(cflags or [])
    ldflags_env = " ".join(ldflags or [])

    env = os.environ.copy()
    env["CC"] = str(wrapper_path)
    if cflags_env:
        env["CFLAGS"] = f"-O2 {cflags_env}"
    else:
        env["CFLAGS"] = "-O2"
    if ldflags_env:
        env["LDFLAGS"] = ldflags_env

    # Use llvm-ar (or zig ar) when available — the system `ar` (GNU ar)
    # produces archives in a format that zig's Mach-O lld can't parse for
    # Apple targets, which then fails at the final link step with "unknown
    # cpu architecture". llvm-ar produces BSD-format archives that are
    # universally readable by lld regardless of target OS.
    llvm_ar = shutil.which("llvm-ar")
    if llvm_ar:
        env["AR"] = llvm_ar
    llvm_ranlib = shutil.which("llvm-ranlib")
    if llvm_ranlib:
        env["RANLIB"] = llvm_ranlib

    # Run configure. We disable every optional feature that pulls in deps we
    # don't ship (PCI, NUMA, XML via libxml2, GPU backends, plugins). The
    # core Linux topology backend uses only sysfs / procfs — no extra deps.
    # --host=... avoids configure treating this as a native build when the
    # wrapper's output doesn't look like the build machine's compiler.
    configure_args = [
        str(src_dir / "configure"),
        f"--prefix={hwloc_prefix}",
        "--enable-static",
        "--disable-shared",
        "--disable-libxml2",
        "--disable-pci",
        "--disable-libnuma",
        "--disable-cairo",
        "--disable-libudev",
        "--disable-cuda",
        "--disable-opencl",
        "--disable-levelzero",
        "--disable-nvml",
        "--disable-rsmi",
        "--disable-valgrind",
        "--disable-ncurses",
        "--disable-plugin-pci",
        "--disable-plugin-opencl",
        "--disable-plugin-cuda",
        "--disable-plugin-rsmi",
        "--disable-plugin-nvml",
        "--disable-plugin-levelzero",
    ]
    if host:
        configure_args.append(f"--host={host}")

    print(f"Building hwloc {HWLOC_VERSION} for {build_root}", flush=True)
    print(f"+ configure {' '.join(configure_args)}", flush=True)
    subprocess.run(configure_args, cwd=hwloc_build, env=env, check=True)

    # Only recurse into include/ (public headers) and hwloc/ (the static
    # library). Skipping utils/ avoids a known zig-cc Mach-O linker failure
    # when the cross-built utils try to link against libhwloc.a — and we
    # don't need the hwloc CLI tools anyway.
    print(f"+ make -C {hwloc_build} SUBDIRS='include hwloc'", flush=True)
    subprocess.run(
        ["make", "SUBDIRS=include hwloc", f"-j{os.cpu_count() or 2}"],
        cwd=hwloc_build,
        env=env,
        check=True,
    )

    print(f"+ make -C {hwloc_build} SUBDIRS='include hwloc' install", flush=True)
    subprocess.run(
        ["make", "SUBDIRS=include hwloc", "install"],
        cwd=hwloc_build,
        env=env,
        check=True,
    )

    stamp.touch()
    return include_dir, libhwloc
