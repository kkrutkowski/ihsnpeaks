x86:
	ln -s ./mufft.x86.h ./mufft.c
	gcc -c -mavx ./mufft.c -Ofast -o mufft_x86.o
	rm ./mufft.c

generic:
	ln -s ./mufft.generic.h ./mufft.c
	gcc -c ./mufft.c -Ofast -o mufft_generic.o
	rm ./mufft.c
