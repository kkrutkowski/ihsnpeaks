.PHONY: debug native clean

CC ?= gcc
MIMALLOC ?= 0
ROOT := ../..
ALLOCATOR_BUILD := system
ifeq ($(MIMALLOC),1)
    ALLOCATOR_BUILD := mimalloc
endif
BUILD_DIR := $(ROOT)/build/native/$(ALLOCATOR_BUILD)
NUFFT_OBJ := $(BUILD_DIR)/nufft1.o
SCALING_HEADER := $(BUILD_DIR)/scaling.h

native:
	$(MAKE) -C $(ROOT) native

debug:
	$(MAKE) -C $(ROOT) $(abspath $(NUFFT_OBJ)) $(abspath $(SCALING_HEADER))
	$(CC) -std=gnu23 -D_GNU_SOURCE -DMAX_TWIDDLE_REUSE=8 -mavx -fsanitize=address -march=x86-64-v2 -O2 -g main.c $(NUFFT_OBJ) -Wl,--gc-sections -lm -I$(ROOT)/include -I$(ROOT)/src/nufft -I$(BUILD_DIR) -o ihsnpeaks

clean:
	$(MAKE) -C $(ROOT) clean
