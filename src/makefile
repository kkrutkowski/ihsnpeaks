all:
	musl-gcc -std=gnu2x -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm #-mx32
	#gcc -std=gnu2x -DMI_OVERRIDE=1 -D_GNU_SOURCE -fsanitize=address -march=x86-64-v3 -Ofast -flto ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm #-mx32
	strip -s ihsnpeaks
