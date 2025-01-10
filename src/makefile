all:
	#musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all ../lib/libmimalloc.a main.c -o ihsnpeaks -Wl,--gc-sections -lm #bugged FFT library, temporarly unavailable
	#musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a
	#gcc -std=gnu23 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f
	gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -fsanitize=address -march=x86-64-v3 -Ofast -flto main.c -o ihsnpeaks -Wl,--gc-sections -lm -L/usr/local/lib -lfftw3f ../lib/libmimalloc.a #-mx32
	strip -s ihsnpeaks

v1:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks-x86-64-v1 -Wl,--gc-sections -lm -L../lib -lfftw3f-v1 -lmimalloc
	strip -s ihsnpeaks-x86-64-v1
v2:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v2  -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks-x86-64-v2 -Wl,--gc-sections -lm -L../lib -lfftw3f-v2 -lmimalloc
	strip -s ihsnpeaks-x86-64-v2
avx:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -mavx -march=x86-64-v2  -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks-x86-64-avx -Wl,--gc-sections -lm -L../lib -lfftw3f-avx -lmimalloc
	strip -s ihsnpeaks-x86-64-avx
v3:
	musl-gcc -std=gnu23 -DMI_OVERRIDE=1 -D_GNU_SOURCE -static -march=x86-64-v3 -Ofast -flto -fno-sanitize=all main.c -o ihsnpeaks-x86-64-v3 -Wl,--gc-sections -lm -L../lib -lfftw3f-v3 -lmimalloc
	strip -s ihsnpeaks-x86-64-v3
rename:
	mv ihsnpeaks-x86-64-avx ihsnpeaks

#FFTW3.10 bootstrap flags;
#v1;    ./configure --enable-float --enable-sse  --disable-fortran
#v2;    ./configure --enable-float --enable-sse --enable-sse2 --disable-fortran
#avx;    ./configure --enable-float --enable-avx  --disable-fortran
#v3;    ./configure --enable-float --enable-avx --enable-avx2 --enable-fma --disable-fortran
