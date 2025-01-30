# Define the minimum required versions
GCC_MIN_VERSION := 14
CLANG_MIN_VERSION := 19
ICX_MIN_VERSION := 2025

#Default CFLAGS
CFLAGS := -std=gnu23 -D_GNU_SOURCE -DMI_OVERRIDE=1 -static -march=native -flto -fno-sanitize=all -Wl,--gc-sections -I../include -lm -L../lib
FFTW_CONFIGURE_FLAGS := --quiet --enable-single --disable-double --disable-fortran
# --enable-neon --enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-fma
# gcc-only: --enable-generic-simd128 --enable-generic-simd256

# Determine the compiler and its version
CC := $(CC)
CC_VERSION := $(shell $(CC) --version | head -n 1)

# Debug: Print CC and CC_VERSION
# $(info CC is set to: $(CC))
# $(info CC_VERSION is: $(CC_VERSION))

# Resolve 'cc' to the actual compiler
ifeq ($(notdir $(CC)),cc)
    CC_RESOLVED := $(shell readlink -f $(shell which $(CC)))
    $(info Resolved 'cc' to: $(CC_RESOLVED))
else
    CC_RESOLVED := $(CC)
endif

# Extract the version number based on the compiler
ifneq ($(findstring gcc, $(CC_RESOLVED)),)
    CC_TYPE := gcc
    CC_VERSION_NUMBER := $(shell $(CC) -dumpversion)
    MIN_VERSION := $(GCC_MIN_VERSION)
    CFLAGS += -Ofast
else ifneq ($(findstring clang, $(CC_RESOLVED)),)
    CC_TYPE := clang
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=clang version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(CLANG_MIN_VERSION)
    CFLAGS += -O3
else ifneq ($(findstring icx, $(CC_RESOLVED)),)
    CC_TYPE := icx
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=icx version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(ICX_MIN_VERSION)
    CFLAGS += -Wno-nan-infinity-disabled -O3
else
    $(error Unsupported compiler: $(CC))
endif

# Extract the major version
CC_MAJOR_VERSION := $(firstword $(subst ., ,$(CC_VERSION_NUMBER)))

# Debug: Print CC_TYPE, CC_VERSION_NUMBER, and CC_MAJOR_VERSION
$(info Compiler type: $(CC_TYPE))
$(info Compiler version: $(CC_VERSION_NUMBER))
#$(info CC_MAJOR_VERSION is: $(CC_MAJOR_VERSION))
$(info flags passed to the compiler are:$(CFLAGS)")

# Target to check the compiler version
check_compiler:
	@echo "Checking compiler version..."
	@if [ -z "$(CC)" ]; then \
	    echo "Error: CC is not set."; \
	    exit 1; \
	fi
	@if [ "$(CC_TYPE)" = "gcc" ]; then \
	    echo "Detected GCC $(CC_VERSION_NUMBER)"; \
	elif [ "$(CC_TYPE)" = "clang" ]; then \
	    echo "Detected Clang $(CC_VERSION_NUMBER)"; \
	elif [ "$(CC_TYPE)" = "icx" ]; then \
	    echo "Detected ICX $(CC_VERSION_NUMBER)"; \
	else \
	    echo "Error: Unsupported compiler: $(CC). Please use one of the following supported compilers and versions:"; \
	    echo "\t- GCC (minimum version $(GCC_MIN_VERSION).0.0)"; \
	    echo "\t- Clang (minimum version $(CLANG_MIN_VERSION).0.0)"; \
	    echo "\t- ICX (minimum version $(ICX_MIN_VERSION).0.0)"; \
	    exit 1; \
	fi
	@if [ $(CC_MAJOR_VERSION) -lt $(MIN_VERSION) ]; then \
	    echo ""; \
	    echo "Error: Unsupported compiler: $(CC_TYPE)-$(CC_VERSION_NUMBER). Please use one of the following supported compilers and versions:"; \
	    echo "\t- GCC (minimum version $(GCC_MIN_VERSION).0.0)"; \
	    echo "\t- Clang (minimum version $(CLANG_MIN_VERSION).0.0)"; \
	    echo "\t- ICX (minimum version $(ICX_MIN_VERSION).0.0)"; \
	    exit 1; \
	fi
	@echo "Compiler check passed: $(CC_TYPE)-$(CC_VERSION_NUMBER)"

# Default target
all:
	@echo ""
	@echo "Building with $(CC_TYPE)-$(CC_VERSION_NUMBER)"
	$(MAKE) clean

# Clean target
clean:
	@echo "Cleaning up..."
