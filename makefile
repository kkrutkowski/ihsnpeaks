# Define the minimum required versions
GCC_MIN_VERSION := 14
CLANG_MIN_VERSION := 18
ICX_MIN_VERSION := 2025

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
else ifneq ($(findstring clang, $(CC_RESOLVED)),)
    CC_TYPE := clang
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=clang version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(CLANG_MIN_VERSION)
else ifneq ($(findstring icx, $(CC_RESOLVED)),)
    CC_TYPE := icx
    CC_VERSION_NUMBER := $(shell $(CC) --version | grep -oP '(?<=icx version )\d+\.\d+\.\d+' | head -n 1)
    MIN_VERSION := $(ICX_MIN_VERSION)
else
    $(error Unsupported compiler: $(CC))
endif

# Extract the major version
CC_MAJOR_VERSION := $(firstword $(subst ., ,$(CC_VERSION_NUMBER)))

# Debug: Print CC_TYPE, CC_VERSION_NUMBER, and CC_MAJOR_VERSION
$(info Compiler type: $(CC_TYPE))
$(info Compiler version: $(CC_VERSION_NUMBER))
#$(info CC_MAJOR_VERSION is: $(CC_MAJOR_VERSION))

# Check if the compiler version meets the minimum requirement
ifeq ($(shell [ $(CC_MAJOR_VERSION) -ge $(MIN_VERSION) ] && echo 1 || echo 0),0)
    $(info )
    $(info Error: Unsupported compiler. Please use one of the following supported compilers and versions:)
    $(info _        - GCC (minimum version $(GCC_MIN_VERSION).0.0))
    $(info _        - Clang (minimum version $(CLANG_MIN_VERSION).0.0))
    $(info _        - ICX (minimum version $(ICX_MIN_VERSION).0.0))
    $(info Detected $(CC_TYPE)-$(CC_VERSION_NUMBER) is not supported.)
    $(error Unsupported compiler)
endif

# Default target
all:
	@echo ""
	@echo "Building with $(CC_TYPE)-$(CC_VERSION_NUMBER)"

# Clean target
clean:
	@echo "Cleaning up..."
