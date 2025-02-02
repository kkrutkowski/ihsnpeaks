# Define the minimum required versions for C23 support
GCC_MIN_VERSION := 14
CLANG_MIN_VERSION := 19
ICX_MIN_VERSION := 2025

# Default CFLAGS.
CFLAGS := -D_GNU_SOURCE -DMI_OVERRIDE=1 -static -march=native -flto -fno-sanitize=all -Wl,--gc-sections -I../include -lm -L../lib
FFTW_CONFIGURE_FLAGS := --quiet --enable-single --disable-double --disable-fortran
FFTW_COMPILER_FLAGS := -march=native -fno-sanitize=all
# --enable-neon --enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-fma
# gcc-only: --enable-generic-simd128 --enable-generic-simd256

# Check if CC is set and points to a valid executable
ifeq ($(shell which $(CC) >/dev/null 2>&1 && echo 1 || echo 0), 0)
    $(error Invalid compiler: '$(CC)' is not a valid executable)
endif

# Determine the compiler and its version
# (Assumes CC is already set externally, e.g. via environment or command line)
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

# Extract the version number based on the compiler type
ifneq ($(findstring gcc, $(CC_RESOLVED)),)
    CC_TYPE := gcc
    CC_VERSION_NUMBER := $(shell $(CC) -dumpversion)
    MIN_VERSION := $(GCC_MIN_VERSION)
    CFLAGS += -Ofast
else ifneq ($(findstring clang, $(CC_RESOLVED)),)
    CC_TYPE := clang
    # Expect version string like "clang version 15.0.0 ..."
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=clang version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(CLANG_MIN_VERSION)
    CFLAGS += -O3
else ifneq ($(findstring icx, $(CC_RESOLVED)),)
    CC_TYPE := icx
    # Expect version string like "icx version 2025.1.0 ..." (for example)
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=icx version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(ICX_MIN_VERSION)
    CFLAGS += -Wno-nan-infinity-disabled -O3
else
    $(error Unsupported compiler: $(CC))
endif

# Extract the major version (the first number before the first dot)
CC_MAJOR_VERSION := $(firstword $(subst ., ,$(CC_VERSION_NUMBER)))

# Debug: Print CC_TYPE, CC_VERSION_NUMBER, and CC_MAJOR_VERSION
$(info Compiler type: $(CC_TYPE))
$(info Compiler version: $(CC_VERSION_NUMBER))
$(info Compiler major version: $(CC_MAJOR_VERSION))
$(info Flags passed to the compiler: $(CFLAGS))

#----------------------------------------------------------------------
# Set the C standard flag based on compiler version:
#
# If the compiler's major version is less than the required minimum
# for C23, then use C11 (gnu11); otherwise, use C23 (gnu23).
#----------------------------------------------------------------------
ifeq ($(shell test $(CC_MAJOR_VERSION) -lt $(MIN_VERSION) && echo true || echo false),true)
    $(info Warning: $(CC_TYPE)-$(CC_VERSION_NUMBER) does not fully support C23. Falling back to gnu11.)
    CFLAGS := $(CFLAGS) -std=gnu11
else
    CFLAGS := $(CFLAGS) -std=gnu23
endif

#----------------------------------------------------------------------
# Targets
#----------------------------------------------------------------------

# Target to check the compiler version.
# In the case where the compiler version is below the C23 requirement,
# a warning is printed, but the build continues (using gnu11).
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
	    echo "Warning: $(CC_TYPE)-$(CC_VERSION_NUMBER) doesn't support C23. Falling back to C11 (gnu11)."; \
	fi
	@if [ $(CC_MAJOR_VERSION) -lt $(MIN_VERSION) ]; then \
	    echo "Warning: $(CC_TYPE)-$(CC_VERSION_NUMBER) doesn't fully support C23. Using C11 standard (-std=gnu11) instead."; \
	else \
	    echo "Compiler check passed: $(CC_TYPE)-$(CC_VERSION_NUMBER) supports C23."; \
	fi

download:
	#@wget …

fftw:
	@wget https://github.com/kkrutkowski/ihsnpeaks/releases/download/beta-1.1.0/fftw-3.3.10_debloated.tar.xz -O /tmp/fftw-3.3.10_ihsnpeaks.tar.xz
	@tar -xf /tmp/fftw-3.3.10_ihsnpeaks.tar.xz -C /tmp
	# Unpacked to: /tmp/fftw-3.3.10/

all: check_compiler
	@echo ""
	@echo "Building with $(CC_TYPE)-$(CC_VERSION_NUMBER)"
	$(MAKE) clean
	# ... add further build commands here ...

clean:
	@echo "Cleaning up..."
	@rm -rf /tmp/fftw-3.3.10/
	@rm -f /tmp/fftw-3.3.10_ihsnpeaks.tar.xz
