all:
	#musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm -lfftw3f #-mx32
	#musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib/libfftw3f.a  #-mx32
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a
	#gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -fsanitize=address -march=x86-64-v3 -Ofast -flto ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib/libfftw3f.a #-mx32
	strip -s ihsnpeaks
