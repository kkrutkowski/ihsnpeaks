# ihsnpeaks
A minimalistic C23 implementation of high-performance multi-harmonic Rayleigh's Z-test and a circularly-independent F test, designed as a replacement for [FNPEAKS and MPEAKS](http://helas.astro.uni.wroc.pl/deliverables.php?active=fnpeaks)

## Installation
### Dependencies
The makefile requires [lscpu](https://man7.org/linux/man-pages/man1/lscpu.1.html), [wget](https://www.gnu.org/software/wget/), [readlink](https://www.gnu.org/software/coreutils/manual/html_node/readlink-invocation.html#readlink-invocation) and a gnu11-compatible C compiler. 

While these binaries should come preinstalled on most modern Linux distributions they may be absent on some of the other POSIX-compliant operating systems. In that case, manual compilation may be required. In that case the compiler has to support at least the gnu99 C standard (with a fallback to ANSI C99 planned for the full 1.0 release).

Additionally the [**photview**](https://github.com/kkrutkowski/ihsnpeaks/blob/main/src/photview) extension (written in the python3 language) relies on [NumPy](https://pypi.org/project/numpy/), [AstroPy](https://pypi.org/project/astropy/), [PyQt6](https://pypi.org/project/PyQt6/) and [PyQt6-Charts](https://pypi.org/project/PyQt6-Charts/).
### Release builds for x86
The release builds targetting the x86 architectures for the Linux operating system can be found at [the releases section of this GitHub repository](https://github.com/kkrutkowski/ihsnpeaks/releases)
Additionally, their installation can be automated by running the makefile. The basic CLI installation of the release build (compiled with gcc-14 compiler, on top of the [MUSL standard library](https://musl.libc.org/) and [mimalloc memory allocator](https://github.com/microsoft/mimalloc)) can be performed in 3 commands as shown below;
```
git clone --depth=1 https://github.com/kkrutkowski/ihsnpeaks && cd ihsnpeaks
make download
sudo make install
```
### Building from source
The makefile installation requires running the makefile with CC set to a gnu11-compatible compiler, although the gnu23 compatibility is expected to result in a noticeable performance increase over the older compilers. To install the application on the systems with no release candidates available please run;
```
git clone --depth=1 https://github.com/kkrutkowski/ihsnpeaks && cd ihsnpeaks
make fftw native clean
sudo make install
```

Alternatively running
```
git clone --depth=1 https://github.com/kkrutkowski/ihsnpeaks && cd ihsnpeaks
make
sudo make install
```
Will make the makefile choose the preferred installation method based on the availability of a suitable release candidate for Linux operating systems, depending on the target's architecture.

## Planned features
* Prewhitening
  * Transit detection and outlier filtering
* MAST .fits photometry support
* Multi-target archive support
  * Custom .phot64 binary archive
  * Tape archive input (.dat + .fits content)
  * Multiband periodograms
  * LZ4 runtime decompression
  * Archive creation
* GUI results evalutation tool

## Credits
