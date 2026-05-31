.PHONY: all native check_compiler install clean format release release-linux-x86_64-musl clean-docker

CC ?= cc
AR ?= ar
MIMALLOC ?= 0

ifeq ($(origin CC),default)
ifneq ($(origin cc),undefined)
    CC := $(cc)
endif
endif

GCC_MIN_VERSION := 14
CLANG_MIN_VERSION := 19
ICX_MIN_VERSION := 2025

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ALLOCATOR_BUILD := system
ifeq ($(MIMALLOC),1)
    ALLOCATOR_BUILD := mimalloc
endif
BUILD_DIR := $(MAKEFILE_DIR)build/native/$(ALLOCATOR_BUILD)
RELEASE_IMAGE ?= ihsnpeaks-release-linux-x86_64
RELEASE_CONTAINER ?= ihsnpeaks-release-extract
RELEASE_BIN := $(MAKEFILE_DIR)dist/ihsnpeaks-linux-x86_64
MIMALLOC_COMPAT_DIR := $(BUILD_DIR)/compat
MIMALLOC_COMPAT_HEADER := $(MIMALLOC_COMPAT_DIR)/mimalloc/mimalloc.h
MIMALLOC_OVERRIDE_HEADER := $(MIMALLOC_COMPAT_DIR)/ihsnpeaks-mimalloc-override.h
NUFFT_SRC := $(MAKEFILE_DIR)src/nufft/nufft1.c
NUFFT_HDR := $(MAKEFILE_DIR)src/nufft/nufft1.h
TRIG_HDR := $(MAKEFILE_DIR)src/utils/trig.h
NUFFT_OBJ := $(BUILD_DIR)/nufft1.o
SCALING_GEN := $(BUILD_DIR)/scaling_gen
SCALING_HEADER := $(BUILD_DIR)/scaling.h

ifeq ($(MIMALLOC),1)
MIMALLOC_HEADER_PATH := $(shell for d in /usr/local/include /usr/include "$$HOME/include"; do [ -d "$$d" ] || continue; p=$$(find "$$d" -maxdepth 3 -type f -name mimalloc.h -print -quit 2>/dev/null); if [ -n "$$p" ]; then printf '%s\n' "$$p"; break; fi; done)
MIMALLOC_STATIC_LIB := $(shell for d in /usr/local/lib /usr/lib "$$HOME/lib"; do [ -d "$$d" ] || continue; p=$$(find "$$d" -maxdepth 4 -type f -name libmimalloc.a -print -quit 2>/dev/null); if [ -n "$$p" ]; then printf '%s\n' "$$p"; break; fi; done)
MIMALLOC_DYNAMIC_LIB := $(shell for d in /usr/local/lib /usr/lib "$$HOME/lib"; do [ -d "$$d" ] || continue; p=$$(find "$$d" -maxdepth 4 -type f -name libmimalloc.so -print -quit 2>/dev/null); [ -n "$$p" ] || p=$$(find "$$d" -maxdepth 4 -type l -name libmimalloc.so -print -quit 2>/dev/null); [ -n "$$p" ] || p=$$(find "$$d" -maxdepth 4 -type f -name 'libmimalloc.so.*' -print -quit 2>/dev/null); [ -n "$$p" ] || p=$$(find "$$d" -maxdepth 4 -type l -name 'libmimalloc.so.*' -print -quit 2>/dev/null); if [ -n "$$p" ]; then printf '%s\n' "$$p"; break; fi; done)
MIMALLOC_HEADER_DIR := $(shell h='$(MIMALLOC_HEADER_PATH)'; if [ -n "$$h" ]; then if [ "$$h" != "$${h%/mimalloc/mimalloc.h}" ]; then printf '%s\n' "$${h%/mimalloc/mimalloc.h}"; elif [ "$$h" != "$${h%/mimalloc.h}" ]; then printf '%s\n' "$${h%/mimalloc.h}"; fi; fi)
MIMALLOC_HEADER_STYLE := $(shell h='$(MIMALLOC_HEADER_PATH)'; if [ -n "$$h" ]; then if [ "$$h" != "$${h%/mimalloc/mimalloc.h}" ]; then printf nested; elif [ "$$h" != "$${h%/mimalloc.h}" ]; then printf flat; fi; fi)
MIMALLOC_DYNAMIC_LIB_DIR := $(dir $(MIMALLOC_DYNAMIC_LIB))

HAS_MIMALLOC_HEADER := 0
MIMALLOC_HEADER_CFLAGS :=
MIMALLOC_HEADER_DEP :=
ifeq ($(MIMALLOC_HEADER_STYLE),nested)
    HAS_MIMALLOC_HEADER := 1
    MIMALLOC_HEADER_CFLAGS := -I$(MIMALLOC_HEADER_DIR)
    MIMALLOC_HEADER_DEP := $(MIMALLOC_OVERRIDE_HEADER)
else ifeq ($(MIMALLOC_HEADER_STYLE),flat)
    HAS_MIMALLOC_HEADER := 1
    MIMALLOC_HEADER_CFLAGS := -I$(MIMALLOC_COMPAT_DIR) -I$(MIMALLOC_HEADER_DIR)
    MIMALLOC_HEADER_DEP := $(MIMALLOC_COMPAT_HEADER) $(MIMALLOC_OVERRIDE_HEADER)
endif

MIMALLOC_STATIC_LINK := 0
ifneq ($(strip $(MIMALLOC_STATIC_LIB)),)
    MIMALLOC_STATIC_LINK := $(shell printf '#include <mimalloc.h>\nint main(void){void *p=mi_malloc(16); mi_free(p); return 0;}\n' | $(CC) $(MIMALLOC_HEADER_CFLAGS) -x c - -x none -static $(MIMALLOC_STATIC_LIB) -o /tmp/ihsnpeaks_mimalloc_static_probe >/dev/null 2>&1 && echo 1 || echo 0)
endif
MIMALLOC_DYNAMIC_LINK := 0
ifneq ($(strip $(MIMALLOC_DYNAMIC_LIB)),)
    MIMALLOC_DYNAMIC_LINK := $(shell printf '#include <mimalloc.h>\nint main(void){void *p=mi_malloc(16); mi_free(p); return 0;}\n' | $(CC) $(MIMALLOC_HEADER_CFLAGS) -x c - -x none -L$(MIMALLOC_DYNAMIC_LIB_DIR) -lmimalloc -o /tmp/ihsnpeaks_mimalloc_dynamic_probe >/dev/null 2>&1 && echo 1 || echo 0)
endif

HAS_MIMALLOC := 0
ALLOCATOR_LDLIBS :=
STATIC_LDFLAGS := -static
ifeq ($(HAS_MIMALLOC_HEADER),1)
    ifeq ($(MIMALLOC_STATIC_LINK),1)
        HAS_MIMALLOC := 1
        ALLOCATOR_LDLIBS := $(MIMALLOC_STATIC_LIB)
        ALLOCATOR_NAME := mimalloc static
    else ifeq ($(MIMALLOC_DYNAMIC_LINK),1)
        HAS_MIMALLOC := 1
        ALLOCATOR_LDLIBS := -L$(MIMALLOC_DYNAMIC_LIB_DIR) -lmimalloc
        STATIC_LDFLAGS :=
        ALLOCATOR_NAME := mimalloc dynamic
    else
        ALLOCATOR_NAME := system
    endif
else
    ALLOCATOR_NAME := system
endif
else
HAS_MIMALLOC_HEADER := 0
MIMALLOC_HEADER_CFLAGS :=
MIMALLOC_HEADER_DEP :=
HAS_MIMALLOC := 0
ALLOCATOR_LDLIBS :=
STATIC_LDFLAGS := -static
ALLOCATOR_NAME := system
endif

ifeq ($(HAS_MIMALLOC),1)
    ALLOCATOR_CFLAGS := -DHAS_MIMALLOC=1 $(MIMALLOC_HEADER_CFLAGS) -include $(MIMALLOC_OVERRIDE_HEADER)
else
    ALLOCATOR_CFLAGS := -DHAS_MIMALLOC=0
endif

CFLAGS_BASE := -D_GNU_SOURCE $(ALLOCATOR_CFLAGS) -march=native -flto -fno-sanitize=all -I$(MAKEFILE_DIR)include -I$(MAKEFILE_DIR)src/nufft -I$(BUILD_DIR)
LDFLAGS_BASE := -Wl,--gc-sections
LDLIBS := -lm $(ALLOCATOR_LDLIBS)

CC_RESOLVED := $(shell readlink -f $(shell command -v $(CC) 2>/dev/null) 2>/dev/null || command -v $(CC) 2>/dev/null)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)

# x86-64 feature-level detection is intentionally retained as build metadata
# for future multi-dispatch/release selection. The current source build still
# uses -march=native while the PSWF NuFFT transition settles.
UNAME_M := $(shell uname -m 2>/dev/null)
CPU_FLAGS := $(shell lscpu 2>/dev/null | sed -n 's/^Flags:[[:space:]]*//p' | head -n 1)
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
	@echo "Allocator: $(ALLOCATOR_NAME)"
	@echo "x86 dispatch level: $(X86_DISPATCH_LEVEL)"
	@echo "x86 candidates: x86-64=$(X86_64_SUPPORTED) x86-64-v2=$(X86_64_V2_SUPPORTED) x86-64-v2+avx=$(X86_64_V2_AVX_SUPPORTED) x86-64-v3=$(X86_64_V3_SUPPORTED) x86-64-v4=$(X86_64_V4_SUPPORTED)"
	@echo "CFLAGS: $(CFLAGS)"

$(BUILD_DIR):
	@mkdir -p $@

$(MIMALLOC_COMPAT_HEADER): | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	@printf '#include <mimalloc.h>\n' > $@

$(MIMALLOC_OVERRIDE_HEADER): $(MIMALLOC_COMPAT_HEADER) | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	@printf '#ifndef IHSNPEAKS_MIMALLOC_OVERRIDE_H\n#define IHSNPEAKS_MIMALLOC_OVERRIDE_H\n#include <mimalloc/mimalloc.h>\n#define malloc(n) mi_malloc(n)\n#define calloc(c,n) mi_calloc(c,n)\n#define realloc(p,n) mi_realloc(p,n)\n#define free(p) mi_free(p)\n#define strdup(s) mi_strdup(s)\n#define strndup(s,n) mi_strndup(s,n)\n#define realpath(f,n) mi_realpath(f,n)\n#define aligned_alloc(a,n) mi_aligned_alloc(a,n)\n#define posix_memalign(p,a,n) mi_posix_memalign(p,a,n)\n#define reallocarray(p,c,n) mi_reallocarray(p,c,n)\n#endif\n' > $@

$(NUFFT_OBJ): $(NUFFT_SRC) $(NUFFT_HDR) $(TRIG_HDR) $(MIMALLOC_HEADER_DEP) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -DMAX_TWIDDLE_REUSE=8 -c $< -o $@

$(SCALING_GEN): $(MAKEFILE_DIR)src/nufft/scaling.c $(NUFFT_OBJ) $(NUFFT_HDR) | $(BUILD_DIR)
	$(CC) -std=gnu11 -O3 -ffast-math -Wall -Wextra -I$(MAKEFILE_DIR)src/nufft -DMAX_TWIDDLE_REUSE=8 -march=native -mtune=native -static $< $(NUFFT_OBJ) $(LDLIBS) -o $@

$(SCALING_HEADER): $(SCALING_GEN)
	$(SCALING_GEN) $@

native: $(NUFFT_OBJ) $(SCALING_HEADER) $(MIMALLOC_HEADER_DEP)
	@echo "Compiling ihsnpeaks with PSWF NuFFT"
	$(CC) $(CFLAGS) $(MAKEFILE_DIR)src/main.c $(NUFFT_OBJ) $(STATIC_LDFLAGS) $(LDFLAGS_BASE) $(LDLIBS) -o ihsnpeaks
	strip -s ihsnpeaks
	@echo "Compilation complete"

release:
	docker build -f $(MAKEFILE_DIR)Dockerfile.release -t $(RELEASE_IMAGE) $(MAKEFILE_DIR)
	@docker rm -f $(RELEASE_CONTAINER) >/dev/null 2>&1 || true
	docker create --name $(RELEASE_CONTAINER) $(RELEASE_IMAGE) >/dev/null
	@mkdir -p $(dir $(RELEASE_BIN))
	docker cp $(RELEASE_CONTAINER):/work/dist/ihsnpeaks-linux-x86_64 $(RELEASE_BIN)
	docker rm $(RELEASE_CONTAINER) >/dev/null
	@echo "Release binary: $(RELEASE_BIN)"

release-linux-x86_64-musl:
	python3 $(MAKEFILE_DIR)dispatch/build_release.py --build-dir $(MAKEFILE_DIR)build/release/linux-x86_64-musl --output $(RELEASE_BIN)

clean-docker:
	@docker rm -f $(RELEASE_CONTAINER) >/dev/null 2>&1 || true

install: native
	@if [ "$$(id -u)" -ne 0 ]; then \
		echo "Error: 'make install' must be run as root."; \
		exit 1; \
	fi
	cp ./ihsnpeaks /usr/local/bin/ihsnpeaks

clean:
	@echo "Cleaning up..."
	@rm -rf $(MAKEFILE_DIR)build
format:
	clang-format -i ./src/*.c
	clang-format -i ./src/*.h
	clang-format -i ./src/utils/*.h
	clang-format -i ./src/nufft/*.c
	clang-format -i ./src/nufft/*.h
	clang-format -i ./include/*.h
	clang-format -i ./include/klib/*.h
