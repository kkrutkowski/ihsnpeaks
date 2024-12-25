all:
	musl-gcc -std=gnu2x -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast main.c -o ihsnpeaks ../lib/libmimalloc.a -Wl,--gc-sections -lm #-mx32 -flto
	strip -s ihsnpeaks
	#gcc -Ofast main.c -o ihsnpeaks -fsanitize=address -Wl,--gc-sections

