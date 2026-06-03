#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from variants import repo_root


MAX_TWIDDLE_REUSE = "8"


@dataclass(frozen=True)
class ArmVariant:
    name: str
    symbol: str
    cflags: tuple[str, ...]
    width: int
    feature: str
    aliases: tuple[str, ...] = ()

    @property
    def build_name(self) -> str:
        return self.symbol.replace("_", "-")


ARM_VARIANTS: tuple[ArmVariant, ...] = (
    ArmVariant("arm64", "arm64", ("-march=armv8-a",), 128, "generic", aliases=("aarch64",)),
    ArmVariant("arm64-neon", "arm64_neon", ("-march=armv8-a+simd",), 128, "neon", aliases=("neon", "asimd")),
    ArmVariant("arm64-sve128", "arm64_sve128", ("-march=armv8.2-a+sve", "-msve-vector-bits=128"), 128, "sve"),
    ArmVariant("arm64-sve256", "arm64_sve256", ("-march=armv8.2-a+sve", "-msve-vector-bits=256"), 256, "sve"),
    ArmVariant("arm64-sve512", "arm64_sve512", ("-march=armv8.2-a+sve", "-msve-vector-bits=512"), 512, "sve"),
    ArmVariant("arm64-sve2-128", "arm64_sve2_128", ("-march=armv8.5-a+sve2", "-msve-vector-bits=128"), 128, "sve2"),
    ArmVariant("arm64-sve2-256", "arm64_sve2_256", ("-march=armv8.5-a+sve2", "-msve-vector-bits=256"), 256, "sve2"),
    ArmVariant("arm64-sve2-512", "arm64_sve2_512", ("-march=armv8.5-a+sve2", "-msve-vector-bits=512"), 512, "sve2"),
)

ARM_DISPATCH_ORDER = (
    "arm64-sve2-512",
    "arm64-sve512",
    "arm64-sve2-256",
    "arm64-sve256",
    "arm64-sve2-128",
    "arm64-sve128",
    "arm64-neon",
    "arm64",
)


def run(cmd: list[str], cwd: Path) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def first_existing_tool(candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        if shutil.which(candidate):
            return candidate
    return candidates[0]


def resolve_tool(tool: str, role: str) -> str:
    resolved = shutil.which(tool)
    if resolved:
        return resolved
    raise SystemExit(
        f"Required {role} tool '{tool}' was not found. "
        "Install the aarch64 musl cross toolchain or pass --cc/--objcopy/--strip explicitly."
    )


def detect_stdflag(cc: str, cwd: Path) -> str:
    for stdflag in ("-std=gnu23", "-std=gnu2x"):
        probe = subprocess.run(
            [cc, stdflag, "-x", "c", "-", "-fsyntax-only"],
            input="constexpr int x = 1; int main(void){return x - 1;}\n",
            text=True,
            cwd=cwd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if probe.returncode == 0:
            return stdflag
    return "-std=gnu11"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def try_generate_proxy_scaling(root: Path, build_root: Path, host_cc: str, stdflag: str, name: str, cflags: list[str]) -> Path | None:
    proxy_dir = build_root / "scaling-proxy" / name
    proxy_dir.mkdir(parents=True, exist_ok=True)
    nufft_obj = proxy_dir / "nufft1.o"
    scaling_gen = proxy_dir / "scaling_gen"
    scaling_header = proxy_dir / "scaling.h"
    common_flags = [stdflag, "-D_GNU_SOURCE", "-DHAS_MIMALLOC=0", f"-DMAX_TWIDDLE_REUSE={MAX_TWIDDLE_REUSE}", "-pthread"]
    include_flags = [f"-I{root / 'include'}", f"-I{root / 'src/nufft'}", f"-I{proxy_dir}", f"-I{root}"]

    try:
        run(
            [
                host_cc,
                *common_flags,
                "-O3",
                "-ffast-math",
                "-fno-sanitize=all",
                *cflags,
                *include_flags,
                "-c",
                str(root / "src/nufft/nufft1.c"),
                "-o",
                str(nufft_obj),
            ],
            root,
        )
        run(
            [
                host_cc,
                *common_flags,
                "-O3",
                "-ffast-math",
                "-fno-sanitize=all",
                *cflags,
                *include_flags,
                str(root / "src/nufft/scaling.c"),
                str(nufft_obj),
                "-static",
                "-lm",
                "-o",
                str(scaling_gen),
            ],
            root,
        )
        run([str(scaling_gen), str(scaling_header)], root)
    except subprocess.CalledProcessError:
        print(f"Scaling proxy {name} could not be generated.", flush=True)
        return None

    return scaling_header


def x86_proxy_scaling_sources(root: Path, build_root: Path, host_cc: str, stdflag: str) -> dict[int, Path]:
    proxies = {
        128: try_generate_proxy_scaling(root, build_root, host_cc, stdflag, "x86-64-v2", ["-march=x86-64-v2", "-mtune=generic"]),
        256: try_generate_proxy_scaling(root, build_root, host_cc, stdflag, "x86-64-v3", ["-march=x86-64-v3", "-mtune=generic"]),
        512: try_generate_proxy_scaling(root, build_root, host_cc, stdflag, "x86-64-v4", ["-march=x86-64-v4", "-mtune=generic"]),
    }
    missing = [str(width) for width, path in proxies.items() if path is None]
    if missing:
        raise SystemExit(
            "Unable to generate required x86 proxy scaling bucket(s) "
            + ", ".join(missing)
            + ". ARM64 cross release on x86 requires runnable x86-64-v2/v3/v4 proxy builds."
        )
    return {width: path for width, path in proxies.items() if path is not None}


def arm_proxy_scaling_sources(root: Path, build_root: Path, host_cc: str, stdflag: str) -> dict[int, Path]:
    arm128 = try_generate_proxy_scaling(root, build_root, host_cc, stdflag, "arm128", ["-march=armv8-a+simd"])
    if arm128 is None:
        raise SystemExit("Unable to generate ARM 128-bit scaling proxy.")

    arm256 = try_generate_proxy_scaling(root, build_root, host_cc, stdflag, "arm256", ["-march=armv8.2-a+sve", "-msve-vector-bits=256"])
    if arm256 is None:
        arm256 = arm128
        print(f"Scaling proxy arm256: reusing shorter proxy {arm256}", flush=True)

    arm512 = try_generate_proxy_scaling(root, build_root, host_cc, stdflag, "arm512", ["-march=armv8.2-a+sve", "-msve-vector-bits=512"])
    if arm512 is None:
        arm512 = arm256
        print(f"Scaling proxy arm512: reusing shorter proxy {arm512}", flush=True)

    return {128: arm128, 256: arm256, 512: arm512}


def variant_wrapper_source(variant: ArmVariant) -> str:
    return f"""/* Generated by dispatch/build_release_arm64.py. */
#define main __attribute__((visibility("default"))) ihsnpeaks_main_{variant.symbol}
#include "src/main.c"
#include "src/nufft/nufft1.c"
"""


def support_function(variant: ArmVariant) -> str:
    if variant.feature == "generic":
        return "supports_arm64"
    if variant.feature == "neon":
        return "supports_arm64_neon"
    return f"supports_{variant.symbol}"


def dispatcher_source(variants: list[ArmVariant]) -> str:
    by_name = {variant.name: variant for variant in variants}
    ordered = [by_name[name] for name in ARM_DISPATCH_ORDER if name in by_name]
    prototypes = "\n".join(f"int ihsnpeaks_main_{variant.symbol}(int argc, char **argv);" for variant in variants)
    table_rows = "\n".join(
        f'    {{"{variant.name}", "{variant.symbol}", ihsnpeaks_main_{variant.symbol}, {support_function(variant)}}},' for variant in ordered
    )
    alias_rows: list[str] = []
    for variant in variants:
        for alias in (variant.name, variant.symbol, variant.build_name, *variant.aliases):
            alias_rows.append(f'    {{"{alias}", "{variant.name}"}},')
    aliases = "\n".join(alias_rows)

    return f"""/* Generated by dispatch/build_release_arm64.py. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/auxv.h>
#include <sys/prctl.h>

#ifndef HWCAP_ASIMD
#define HWCAP_ASIMD (1UL << 1)
#endif
#ifndef HWCAP_SVE
#define HWCAP_SVE (1UL << 22)
#endif
#ifndef HWCAP2_SVE2
#define HWCAP2_SVE2 (1UL << 1)
#endif
#ifndef PR_SVE_GET_VL
#define PR_SVE_GET_VL 51
#endif
#ifndef PR_SVE_VL_LEN_MASK
#define PR_SVE_VL_LEN_MASK 0xffff
#endif

{prototypes}

typedef int (*ihsnpeaks_entry_fn)(int argc, char **argv);

typedef struct {{
    int asimd;
    int sve;
    int sve2;
    int sve_bits;
}} ihsnpeaks_cpu_features;

typedef struct {{
    const char *name;
    const char *symbol;
    ihsnpeaks_entry_fn entry;
    int (*supported)(const ihsnpeaks_cpu_features *features);
}} ihsnpeaks_variant_entry;

typedef struct {{
    const char *alias;
    const char *canonical;
}} ihsnpeaks_variant_alias;

static ihsnpeaks_cpu_features ihsnpeaks_detect_cpu(void) {{
    ihsnpeaks_cpu_features features = {{0}};
    unsigned long hwcap = getauxval(AT_HWCAP);
    unsigned long hwcap2 = getauxval(AT_HWCAP2);
    features.asimd = (hwcap & HWCAP_ASIMD) != 0;
    features.sve = (hwcap & HWCAP_SVE) != 0;
    features.sve2 = (hwcap2 & HWCAP2_SVE2) != 0;
    if (features.sve) {{
        long vl = prctl(PR_SVE_GET_VL);
        if (vl > 0) features.sve_bits = (int)((vl & PR_SVE_VL_LEN_MASK) * 8L);
    }}
    return features;
}}

static int supports_arm64(const ihsnpeaks_cpu_features *features) {{
    (void)features;
    return 1;
}}

static int supports_arm64_neon(const ihsnpeaks_cpu_features *features) {{
    return features->asimd;
}}

{''.join(f'''
static int supports_{variant.symbol}(const ihsnpeaks_cpu_features *features) {{
    return features->{variant.feature} && features->sve_bits >= {variant.width};
}}
''' for variant in variants if variant.feature in {'sve', 'sve2'})}

static const ihsnpeaks_variant_entry ihsnpeaks_variants[] = {{
{table_rows}
}};

static const ihsnpeaks_variant_alias ihsnpeaks_aliases[] = {{
{aliases}
}};

static const char *canonical_variant_name(const char *name) {{
    size_t n = sizeof(ihsnpeaks_aliases) / sizeof(ihsnpeaks_aliases[0]);
    for (size_t i = 0; i < n; ++i) {{
        if (strcmp(name, ihsnpeaks_aliases[i].alias) == 0) return ihsnpeaks_aliases[i].canonical;
    }}
    return name;
}}

static void print_dispatch_variants(const ihsnpeaks_cpu_features *features) {{
    size_t n = sizeof(ihsnpeaks_variants) / sizeof(ihsnpeaks_variants[0]);
    fputs("Compiled ihsnpeaks dispatch variants:\\n", stdout);
    for (size_t i = 0; i < n; ++i) {{
        const ihsnpeaks_variant_entry *entry = &ihsnpeaks_variants[i];
        printf("  %s%s\\n", entry->name, entry->supported(features) ? "" : " (unsupported on this CPU)");
    }}
}}

int main(int argc, char **argv) {{
    ihsnpeaks_cpu_features features = ihsnpeaks_detect_cpu();
    const char *forced = getenv("IHSNPEAKS_DISPATCH");
    size_t n = sizeof(ihsnpeaks_variants) / sizeof(ihsnpeaks_variants[0]);

    if (forced && (strcmp(forced, "help") == 0 || strcmp(forced, "list") == 0)) {{
        print_dispatch_variants(&features);
        return 0;
    }}

    if (forced && forced[0] != '\\0') {{
        const char *canonical = canonical_variant_name(forced);
        for (size_t i = 0; i < n; ++i) {{
            const ihsnpeaks_variant_entry *entry = &ihsnpeaks_variants[i];
            if (strcmp(canonical, entry->name) == 0) {{
                if (!entry->supported(&features)) {{
                    fprintf(stderr, "Requested dispatch variant '%s' is not supported on this CPU.\\n", forced);
                    return 2;
                }}
                return entry->entry(argc, argv);
            }}
        }}
        fprintf(stderr, "Unknown dispatch variant '%s'. Set IHSNPEAKS_DISPATCH=help to list compiled variants.\\n", forced);
        return 2;
    }}

    for (size_t i = 0; i < n; ++i) {{
        const ihsnpeaks_variant_entry *entry = &ihsnpeaks_variants[i];
        if (entry->supported(&features)) return entry->entry(argc, argv);
    }}

    fputs("No supported ihsnpeaks dispatch variant is available.\\n", stderr);
    return 2;
}}
"""


def compile_variant(
    root: Path,
    build_root: Path,
    cc: str,
    objcopy: str,
    stdflag: str,
    base_cflags: list[str],
    variant: ArmVariant,
    scaling_source: Path,
) -> Path:
    variant_dir = build_root / variant.build_name
    variant_dir.mkdir(parents=True, exist_ok=True)
    scaling_header = variant_dir / "scaling.h"
    wrapper_c = variant_dir / "ihsnpeaks_variant.c"
    variant_obj = variant_dir / "ihsnpeaks_variant.o"
    include_flags = [f"-I{root / 'include'}", f"-I{root / 'src/nufft'}", f"-I{variant_dir}", f"-I{root}"]

    print(f"Scaling proxy for {variant.name}: {scaling_source}", flush=True)
    shutil.copyfile(scaling_source, scaling_header)
    write_text(wrapper_c, variant_wrapper_source(variant))
    run(
        [
            cc,
            *base_cflags,
            *variant.cflags,
            *include_flags,
            "-c",
            str(wrapper_c),
            "-o",
            str(variant_obj),
        ],
        root,
    )
    run([objcopy, "--localize-hidden", str(variant_obj)], root)
    return variant_obj


def default_target_tool(name: str) -> str:
    machine = platform.machine().lower()
    if machine in {"aarch64", "arm64"}:
        return name
    return first_existing_tool((f"aarch64-alpine-linux-musl-{name}", f"aarch64-linux-musl-{name}"))


def main(argv: list[str]) -> int:
    root = repo_root()
    parser = argparse.ArgumentParser(description="Build the static Linux ARM64 release dispatcher.")
    parser.add_argument("--cc", default=os.environ.get("CC", default_target_tool("gcc")))
    parser.add_argument("--objcopy", default=os.environ.get("OBJCOPY", default_target_tool("objcopy")))
    parser.add_argument("--strip", default=os.environ.get("STRIP", default_target_tool("strip")))
    parser.add_argument("--host-cc", default=os.environ.get("HOST_CC", "cc"))
    parser.add_argument("--build-dir", default=str(root / "build/release/linux-arm64-musl"))
    parser.add_argument("--output", default=str(root / "dist/ihsnpeaks-linux-arm64"))
    args = parser.parse_args(argv)

    cc = resolve_tool(args.cc, "target compiler")
    objcopy = resolve_tool(args.objcopy, "target objcopy")
    strip = resolve_tool(args.strip, "target strip")
    host_cc = resolve_tool(args.host_cc, "host compiler")
    build_root = Path(args.build_dir).resolve()
    output = Path(args.output).resolve()
    build_root.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    host_stdflag = detect_stdflag(host_cc, root)
    host_machine = platform.machine().lower()
    if host_machine in {"x86_64", "amd64"}:
        scaling_by_width = x86_proxy_scaling_sources(root, build_root, host_cc, host_stdflag)
    elif host_machine in {"aarch64", "arm64"}:
        scaling_by_width = arm_proxy_scaling_sources(root, build_root, host_cc, host_stdflag)
    else:
        raise SystemExit(f"release-linux-arm64-musl requires x86_64 or arm64 host/proxy scaling, got {platform.machine()}")

    variants = list(ARM_VARIANTS)
    print("Release variants: " + ", ".join(variant.name for variant in variants), flush=True)
    stdflag = detect_stdflag(cc, root)
    base_cflags = [
        stdflag,
        "-D_GNU_SOURCE",
        "-DHAS_MIMALLOC=0",
        f"-DMAX_TWIDDLE_REUSE={MAX_TWIDDLE_REUSE}",
        "-Ofast",
        "-fno-sanitize=all",
        "-ffunction-sections",
        "-fdata-sections",
        "-fvisibility=hidden",
        "-pthread",
    ]

    variant_objects = [
        compile_variant(root, build_root, cc, objcopy, stdflag, base_cflags, variant, scaling_by_width[variant.width]) for variant in variants
    ]

    dispatcher_c = build_root / "dispatch_main.c"
    dispatcher_o = build_root / "dispatch_main.o"
    write_text(dispatcher_c, dispatcher_source(variants))
    run(
        [
            cc,
            stdflag,
            "-D_GNU_SOURCE",
            "-O2",
            "-fno-sanitize=all",
            "-ffunction-sections",
            "-fdata-sections",
            "-march=armv8-a",
            "-pthread",
            "-c",
            str(dispatcher_c),
            "-o",
            str(dispatcher_o),
        ],
        root,
    )
    run(
        [
            cc,
            "-static",
            "-pthread",
            "-Wl,--gc-sections",
            str(dispatcher_o),
            *(str(obj) for obj in variant_objects),
            "-lm",
            "-o",
            str(output),
        ],
        root,
    )
    run([strip, "-s", str(output)], root)
    print(f"Release binary: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
