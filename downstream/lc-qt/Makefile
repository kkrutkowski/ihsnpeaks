# Makefile for compiling and compressing lc-qt statically with Buildroot musl

.PHONY: all clean build-all

all: lc-qt

lc-qt: build-all
	@echo "🗜️ Compressing executable with UPX..."
	docker run --rm -v "$(shell pwd)":/app lc-qt-builder sh -c " \
		/app/buildroot/bin/x86_64-buildroot-linux-musl-strip /app/build/lc-qt; \
		strip -s /app/build/lc-qt; \
		upx /app/build/lc-qt; \
		chown $(shell id -u):$(shell id -g) /app/build/lc-qt; \
		chmod 755 /app/build/lc-qt; \
		mv /app/build/lc-qt /app/lc-qt; \
		rm -rf /app/build"
	@echo "🎉 Done! Statically compiled and compressed executable is ready at: ./lc-qt"


build-all:
	@echo "🐳 Building Docker image and running static compile..."
	chmod +x build.sh
	./build.sh

clean:
	rm -rf build

deep-clean: clean
	@echo "🧹 Removing entire buildroot cache directory..."
	rm -rf buildroot

docker-clean:
	@echo "🐳 Removing builder Docker containers and image..."
	docker rm -f $$(docker ps -a -q --filter ancestor=lc-qt-builder) 2>/dev/null || true
	docker rmi -f lc-qt-builder 2>/dev/null || true


