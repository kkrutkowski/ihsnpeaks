all:
	#musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm #bugged FFT library, temporarly unavailable
	#musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a
	#gcc -std=gnu23 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f
	gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -fsanitize=address -march=x86-64-v3 -Ofast -flto main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a #-mx32
	strip -s ihsnpeaks

v1:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a
	strip -s ihsnpeaks
v2:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v2 -mavx -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a
	strip -s ihsnpeaks
v3:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a
	strip -s ihsnpeaks
