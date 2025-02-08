# ihsnpeaks
A minimalistic C23 implementation of high-performance multi-harmonic Rayleigh's Z-test and a circularly-independent F test, designed as a replacement for [FNPEAKS and MPEAKS](http://helas.astro.uni.wroc.pl/deliverables.php?active=fnpeaks)

## Installation
### Release builds for x86
The release builds targetting the x86 architectures for the Linux operating system can be found at [the releases section of this GitHub repo](https://github.com/kkrutkowski/ihsnpeaks/releases)
Additionally, their installation can be automated by running the makefile. The basic CLI installation of the release build (compiled with gcc-14 compiler, on top of the [MUSL standard library](https://musl.libc.org/) and [mimalloc memory allocator](https://github.com/microsoft/mimalloc)) can be performed in 3 commands as shown below
```
git clone https://github.com/kkrutkowski/ihsnpeaks && cd ihsnpeaks
make download
sudo make install
```
### Building from source
The makefile installation requires running the makefile with CC set to a gnu11-compatible compiler, although the gnu23 compatibility is expected to result in a noticeable performance increase over the older compilers. To install the application on the systems with no release candidates available please run;
```
git clone https://github.com/kkrutkowski/ihsnpeaks && cd ihsnpeaks
make fftw native clean
sudo make install
```
## Credits
