# Makefile for compiling and compressing lc-qt statically with Buildroot musl

.PHONY: all clean build-all

all: build/lc-qt

build/lc-qt: build-all
	@echo "🗜️ Compressing executable with UPX and cleaning up..."
	docker run --rm -v "$(shell pwd)":/app lc-qt-builder sh -c " \
		strip -s /app/build/lc-qt; \
		upx -1 /app/build/lc-qt; \
		mv /app/build/lc-qt /app/lc-qt-temp; \
		rm -rf /app/build /app/x86_64-buildroot-linux-musl_sdk-buildroot; \
		mkdir -p /app/build; \
		mv /app/lc-qt-temp /app/build/lc-qt; \
		chown -R $(shell id -u):$(shell id -g) /app/build; \
		chmod 755 /app/build/lc-qt"
	@echo "🎉 Done! Statically compiled and compressed executable is ready at: build/lc-qt"


build-all:
	@echo "🐳 Building Docker image and running static compile..."
	chmod +x build.sh
	./build.sh

clean:
	rm -rf build x86_64-buildroot-linux-musl_sdk-buildroot
