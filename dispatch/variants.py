from __future__ import annotations

import platform
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Variant:
    name: str
    symbol: str
    cflags: tuple[str, ...]
    required_flags: tuple[str | tuple[str, ...], ...]
    aliases: tuple[str, ...] = ()

    @property
    def build_name(self) -> str:
        return self.symbol.replace("_", "-")


VARIANTS: tuple[Variant, ...] = (
    Variant(
        name="x86-64",
        symbol="x86_64",
        cflags=("-march=x86-64", "-mtune=generic"),
        required_flags=(),
    ),
    Variant(
        name="x86-64-v2",
        symbol="x86_64_v2",
        cflags=("-march=x86-64-v2", "-mtune=generic"),
        required_flags=("cx16", "lahf_lm", "popcnt", "sse4_1", "sse4_2", "ssse3"),
    ),
    Variant(
        name="x86-64-v2+avx",
        symbol="x86_64_v2_avx",
        cflags=("-march=x86-64-v2", "-mavx", "-mxsave", "-mtune=generic"),
        required_flags=("cx16", "lahf_lm", "popcnt", "sse4_1", "sse4_2", "ssse3", "avx", "xsave"),
        aliases=("x86-64-v2-avx", "x86-64-v2+avx+xsave", "x86-64-v2-avx-xsave"),
    ),
    Variant(
        name="x86-64-v3",
        symbol="x86_64_v3",
        cflags=("-march=x86-64-v3", "-mtune=generic"),
        required_flags=(
            "cx16",
            "lahf_lm",
            "popcnt",
            "sse4_1",
            "sse4_2",
            "ssse3",
            "avx",
            "xsave",
            "avx2",
            "fma",
            "bmi1",
            "bmi2",
            "f16c",
            ("abm", "lzcnt"),
            "movbe",
        ),
    ),
    Variant(
        name="x86-64-v4",
        symbol="x86_64_v4",
        cflags=("-march=x86-64-v4", "-mtune=generic"),
        required_flags=(
            "cx16",
            "lahf_lm",
            "popcnt",
            "sse4_1",
            "sse4_2",
            "ssse3",
            "avx",
            "xsave",
            "avx2",
            "fma",
            "bmi1",
            "bmi2",
            "f16c",
            ("abm", "lzcnt"),
            "movbe",
            "avx512f",
            "avx512bw",
            "avx512cd",
            "avx512dq",
            "avx512vl",
        ),
    ),
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def cpu_flags() -> set[str]:
    flags: set[str] = set()
    try:
        text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return flags

    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep and key.strip().lower() in {"flags", "features"}:
            flags.update(value.split())
            break
    return flags


def requirement_met(requirement: str | tuple[str, ...], flags: set[str]) -> bool:
    if isinstance(requirement, tuple):
        return any(flag in flags for flag in requirement)
    return requirement in flags


def host_supported_variants() -> list[Variant]:
    machine = platform.machine().lower()
    if machine not in {"x86_64", "amd64"}:
        raise SystemExit(f"release-linux-x86_64-musl requires x86_64, got {platform.machine()}")

    flags = cpu_flags()
    supported: list[Variant] = []
    for variant in VARIANTS:
        if all(requirement_met(requirement, flags) for requirement in variant.required_flags):
            supported.append(variant)

    if not supported or supported[0].name != "x86-64":
        supported.insert(0, VARIANTS[0])
    return supported
