CC = gcc
CXX = nvcc
CFLAGS = -O3 -march=native -Wall -Wextra -fopenmp
NVCCFLAGS = -O3 -arch=sm_120 -Xcompiler -O3,-march=native
LIBS = -lm -lcudart -lgomp

# Include paths
INCLUDES = -I../../src/nufft -I../src/nufft

all: 1_nufft

# Compile the CPU NUFFT dependency
nufft1.o: ../../src/nufft/nufft1.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Compile the GPU NUFFT code
cunufft.o: src/nufft/cunufft.cu
	$(CXX) $(NVCCFLAGS) -c $< -o $@

# Build the test binary
1_nufft: tests/1_nufft.c cunufft.o nufft1.o
	$(CXX) $(NVCCFLAGS) -Xcompiler -fopenmp $(INCLUDES) $^ -o $@ $(LIBS)

clean:
	rm -f *.o 1_nufft src/nufft/*.o tests/*.o

run: 1_nufft
	./1_nufft

.PHONY: all clean run
