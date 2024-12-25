all:
	musl-gcc -std=gnu2x -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast main.c -o ../lib/libmimalloc.a ihsnpeaks -Wl,--gc-sections -lm #-mx32 -flto
	strip -s ihsnpeaks
	#gcc -Ofast main.c -o ihsnpeaks -fsanitize=address -Wl,--gc-sections

