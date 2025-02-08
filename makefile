.PHONY: check_compiler download fftw mimalloc clean native all

# Define the minimum required versions for C23 support
GCC_MIN_VERSION := 14
CLANG_MIN_VERSION := 19
ICX_MIN_VERSION := 2025

# Get the directory where the Makefile is located
MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Default CFLAGS.
CFLAGS := -D_GNU_SOURCE -DMI_OVERRIDE=1 -static -march=native -flto -fno-sanitize=all -Wl,--gc-sections -I../include -lm -L../lib
FFTW_CONFIGURE_FLAGS := --quiet --enable-single --disable-double --disable-fortran
FFTW_COMPILER_FLAGS := -march=native -fno-sanitize=all
# --enable-neon --enable-avx --enable-avx2 --enable-avx512 --enable-fma
# gcc-only:

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

# Features organized by microarchitecture versions
FEATURES_V1 := lm cmov cx8 fpu fxsr mmx syscall sse2
FEATURES_V2 := cx16 lahf_lm popcnt sse4_1 sse4_2 ssse3 avx
FEATURES_V3 := avx2 bmi1 bmi2 f16c fma abm movbe xsave
FEATURES_V4 := avx512f avx512bw avx512cd avx512dq avx512vl

# Combine all features into a single list
FEATURES := $(FEATURES_V1) $(FEATURES_V2) $(FEATURES_V3) $(FEATURES_V4)

# Get CPU flags using lscpu
CPU_FLAGS := $(shell lscpu | grep -oP 'Flags:\s*\K.*')

# Loop through each feature and check if it’s supported.
$(foreach feature, $(FEATURES), \
    $(if $(filter $(feature),$(CPU_FLAGS)), \
        $(eval CFLAGS += -m$(feature)) \
        $(info $(shell echo $(feature) | tr 'a-z' 'A-Z') is supported, enabling $(shell echo $(feature) | tr 'a-z' 'A-Z')), \
        $(info $(shell echo $(feature) | tr 'a-z' 'A-Z') is not supported, disabling $(shell echo $(feature) | tr 'a-z' 'A-Z')) \
    ) \
)

# Initialize version to 0 (generic build)
version := 0

# Loop over version groups in descending order.
$(foreach ver,4 3 2 1, \
    $(if $(strip $(foreach f, $(FEATURES_V$(ver)), $(filter $(f),$(CPU_FLAGS)))), \
         $(if $(filter 0,$(version)), $(eval version := $(ver))) \
    ) \
)

ifeq ($(version),0)
  FFTW_CONFIGURE_FLAGS += --enable-generic-simd128 --enable-generic-simd256
else ifeq ($(version),1)
   FFTW_CONFIGURE_FLAGS += --enable-sse --enable-sse2
else ifeq ($(version),2)
  FFTW_CONFIGURE_FLAGS += --enable-avx
else ifeq ($(version),3)
  FFTW_CONFIGURE_FLAGS += --enable-avx2 --enable-fma
else ifeq ($(version),4)
  FFTW_CONFIGURE_FLAGS += --enable-avx2 --enable-avx512
else
  $(error Unknown version: $(version))
endif


#----------------------------------------------------------------------
# Set the C standard flag based on compiler version:
#
# If the compiler's major version is less than the required minimum
# for C23, then use C11 (gnu11); otherwise, use C23 (gnu23).
#----------------------------------------------------------------------
ifeq ($(shell test $(CC_MAJOR_VERSION) -lt $(MIN_VERSION) && echo true || echo false),true)
    $(info Warning: $(CC_TYPE)-$(CC_VERSION_NUMBER) does not fully support C23. Falling back to gnu11.)
    CFLAGS += $(CFLAGS) -std=gnu11
else
    CFLAGS += $(CFLAGS) -std=gnu23
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
	#@wget https://github.com/kkrutkowski/ihsnpeaks/releases/download/beta-1.1.0/……… -O $(MAKEFILE_DIR)/ihsnpeaks
	#@chmod +x $(MAKEFILE_DIR)/ihsnpeaks

install:
	@if [ -x ./ihsnpeaks ]; then \
		if [ "$$(id -u)" -ne 0 ]; then \
			echo "Error: 'make install' must be run as root."; \
			exit 1; \
		fi; \
		cp ./ihsnpeaks /usr/local/bin/ihsnpeaks; \
		/usr/local/bin/ihsnpeaks generate; \
	else \
		echo "Error: executable not found."; \
		echo "Please run 'make download' or 'make native' to build it first."; \
		exit 1; \
	fi

fftw:
	@mkdir -p $(MAKEFILE_DIR)/lib
	@wget https://github.com/kkrutkowski/ihsnpeaks/releases/download/beta-1.1.0/fftw-3.3.10_debloated.tar.xz -O /tmp/fftw-3.3.10_ihsnpeaks.tar.xz
	@tar -xf /tmp/fftw-3.3.10_ihsnpeaks.tar.xz -C /tmp
	# Unpacked to: /tmp/fftw-3.3.10/
clean:
	@echo "Cleaning up..."
	@rm -rf /tmp/fftw-3.3.10 /tmp/fftw-3.3.10_ihsnpeaks.tar.xz || true
