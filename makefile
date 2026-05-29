.PHONY: all native check_compiler install clean

CC ?= cc
AR ?= ar

ifeq ($(origin CC),default)
ifneq ($(origin cc),undefined)
    CC := $(cc)
endif
endif

GCC_MIN_VERSION := 14
CLANG_MIN_VERSION := 19
ICX_MIN_VERSION := 2025

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR := $(MAKEFILE_DIR)build/native
NUFFT_SRC := $(MAKEFILE_DIR)src/nufft/nufft1.c
NUFFT_HDR := $(MAKEFILE_DIR)src/nufft/nufft1.h
NUFFT_OBJ := $(BUILD_DIR)/nufft1.o
SCALING_GEN := $(BUILD_DIR)/scaling_gen
SCALING_HEADER := $(BUILD_DIR)/scaling.h

CFLAGS_BASE := -D_GNU_SOURCE -DMI_OVERRIDE=1 -march=native -flto -fno-sanitize=all -I$(MAKEFILE_DIR)include -I$(MAKEFILE_DIR)src/nufft -I$(BUILD_DIR) -static
LDFLAGS_BASE := -Wl,--gc-sections
LDLIBS := -lm

CC_RESOLVED := $(shell readlink -f $(shell command -v $(CC) 2>/dev/null) 2>/dev/null || command -v $(CC) 2>/dev/null)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)

# x86-64 feature-level detection is intentionally retained as build metadata
# for future multi-dispatch/release selection. The current source build still
# uses -march=native while the PSWF NuFFT transition settles.
UNAME_M := $(shell uname -m 2>/dev/null)
CPU_FLAGS := $(shell lscpu 2>/dev/null | grep -oP 'Flags:\s*\K.*' | head -n 1)
ifeq ($(strip $(CPU_FLAGS)),)
    CPU_FLAGS := $(shell grep -m1 -E '^(flags|Features)[[:space:]]*:' /proc/cpuinfo 2>/dev/null | cut -d: -f2)
endif

FEATURES_X86_64 := sse sse2 lm cmov cx8 fpu fxsr mmx syscall
FEATURES_X86_64_V2 := cx16 lahf_lm popcnt sse4_1 sse4_2 ssse3
FEATURES_X86_64_V2_AVX := avx xsave
FEATURES_X86_64_V3 := avx2 fma bmi1 bmi2 f16c abm movbe
FEATURES_X86_64_V4 := avx512f avx512bw avx512cd avx512dq avx512vl

missing_features = $(filter-out $(CPU_FLAGS),$(1))
has_features = $(if $(strip $(call missing_features,$(1))),0,1)

X86_64_SUPPORTED := 0
X86_64_V2_SUPPORTED := 0
X86_64_V2_AVX_SUPPORTED := 0
X86_64_V3_SUPPORTED := 0
X86_64_V4_SUPPORTED := 0
X86_DISPATCH_LEVEL := generic
X86_DISPATCH_VERSION := 0

ifneq ($(filter x86_64 amd64,$(UNAME_M)),)
    X86_64_SUPPORTED := $(call has_features,$(FEATURES_X86_64))
    ifeq ($(X86_64_SUPPORTED),1)
        X86_DISPATCH_LEVEL := x86-64
        X86_DISPATCH_VERSION := 1
    endif
    ifeq ($(call has_features,$(FEATURES_X86_64) $(FEATURES_X86_64_V2)),1)
        X86_64_V2_SUPPORTED := 1
        X86_DISPATCH_LEVEL := x86-64-v2
        X86_DISPATCH_VERSION := 2
    endif
    ifeq ($(call has_features,$(FEATURES_X86_64) $(FEATURES_X86_64_V2) $(FEATURES_X86_64_V2_AVX)),1)
        X86_64_V2_AVX_SUPPORTED := 1
        X86_DISPATCH_LEVEL := x86-64-v2+avx
        X86_DISPATCH_VERSION := avx
    endif
    ifeq ($(call has_features,$(FEATURES_X86_64) $(FEATURES_X86_64_V2) $(FEATURES_X86_64_V2_AVX) $(FEATURES_X86_64_V3)),1)
        X86_64_V3_SUPPORTED := 1
        X86_DISPATCH_LEVEL := x86-64-v3
        X86_DISPATCH_VERSION := 3
    endif
    ifeq ($(call has_features,$(FEATURES_X86_64) $(FEATURES_X86_64_V2) $(FEATURES_X86_64_V2_AVX) $(FEATURES_X86_64_V3) $(FEATURES_X86_64_V4)),1)
        X86_64_V4_SUPPORTED := 1
        X86_DISPATCH_LEVEL := x86-64-v4
        X86_DISPATCH_VERSION := 4
    endif
endif

ifneq ($(findstring gcc,$(CC_RESOLVED)),)
    CC_TYPE := gcc
    CC_VERSION_NUMBER := $(shell $(CC) -dumpversion)
    MIN_VERSION := $(GCC_MIN_VERSION)
    OPTFLAGS := -Ofast
else ifneq ($(findstring clang,$(CC_RESOLVED)),)
    CC_TYPE := clang
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=clang version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(CLANG_MIN_VERSION)
    OPTFLAGS := -O3 -Wno-nan-infinity-disabled
else ifneq ($(findstring icx,$(CC_RESOLVED)),)
    CC_TYPE := icx
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=icx version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(ICX_MIN_VERSION)
    OPTFLAGS := -O3 -Wno-nan-infinity-disabled
else
    CC_TYPE := unknown
    CC_VERSION_NUMBER := 0
    MIN_VERSION := 0
    OPTFLAGS := -O3
endif

CC_MAJOR_VERSION := $(firstword $(subst ., ,$(CC_VERSION_NUMBER)))
ifeq ($(CC_MAJOR_VERSION),)
    CC_MAJOR_VERSION := 0
endif

ifeq ($(shell test $(CC_MAJOR_VERSION) -lt $(MIN_VERSION) && echo true || echo false),true)
    STDFLAG := -std=gnu11
else
    STDFLAG := -std=gnu23
endif

CFLAGS := $(STDFLAG) $(OPTFLAGS) $(CFLAGS_BASE)

all: native

check_compiler:
	@echo "Compiler command: $(CC)"
	@echo "Compiler: $(CC_VERSION)"
	@echo "C standard: $(STDFLAG)"
	@echo "Host machine: $(UNAME_M)"
	@echo "x86 dispatch level: $(X86_DISPATCH_LEVEL)"
	@echo "x86 candidates: x86-64=$(X86_64_SUPPORTED) x86-64-v2=$(X86_64_V2_SUPPORTED) x86-64-v2+avx=$(X86_64_V2_AVX_SUPPORTED) x86-64-v3=$(X86_64_V3_SUPPORTED) x86-64-v4=$(X86_64_V4_SUPPORTED)"
	@echo "CFLAGS: $(CFLAGS)"

$(BUILD_DIR):
	@mkdir -p $@

$(NUFFT_OBJ): $(NUFFT_SRC) $(NUFFT_HDR) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -DMAX_TWIDDLE_REUSE=8 -c $< -o $@

$(SCALING_GEN): $(MAKEFILE_DIR)src/nufft/scaling.c $(NUFFT_OBJ) $(NUFFT_HDR) | $(BUILD_DIR)
	$(CC) -std=gnu11 -O3 -ffast-math -Wall -Wextra -I$(MAKEFILE_DIR)src/nufft -DMAX_TWIDDLE_REUSE=8 -march=native -mtune=native -static $< $(NUFFT_OBJ) $(LDLIBS) -o $@

$(SCALING_HEADER): $(SCALING_GEN)
	$(SCALING_GEN) $@

native: $(NUFFT_OBJ) $(SCALING_HEADER)
	@echo "Compiling ihsnpeaks with PSWF NuFFT"
	$(CC) $(CFLAGS) $(MAKEFILE_DIR)src/main.c $(NUFFT_OBJ) $(LDFLAGS_BASE) $(LDLIBS) -o ihsnpeaks
	strip -s ihsnpeaks
	@echo "Compilation complete"

install: native
	@if [ "$$(id -u)" -ne 0 ]; then \
		echo "Error: 'make install' must be run as root."; \
		exit 1; \
	fi
	cp ./ihsnpeaks /usr/local/bin/ihsnpeaks

clean:
	@echo "Cleaning up..."
	@rm -rf $(MAKEFILE_DIR)build
