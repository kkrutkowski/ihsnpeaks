# lc-qt — Project Guide

## Project Overview

`lc-qt` is a clean-room rewrite of Prof. Grzegorz Pojmański's **lc** light-curve analysis tool (1997–2018), reimplemented with a Qt6 GUI frontend and **ihsnpeaks** as the backend periodogram engine. The project is a work in progress.

The application renders a dark-themed UI replicating the original `lc` tool's layout: raw/phased light-curve plots, object list management, classification system with customizable labels, period search (AoV, IHS, GB, BLS), and period modify controls. Plot widgets are currently mock `QFrame` subclasses with painted grids — no actual data rendering or ihsnpeaks backend integration yet.

**Parent project:** This is a downstream subproject of `ihsnpeaks` — located at `ihsnpeaks/downstream/lc-qt/`. It uses `json.h` from the parent project (`ihsnpeaks/include/json.h`) for configuration persistence.

## Architecture

- **Qt6 Widgets:** Uses `QMainWindow`, `QFrame`, `QGroupBox`, `QLineEdit`, `QPushButton`, `QCheckBox` — classic widget-based desktop UI, no QML.
- **Dark theme:** Inline `setStyleSheet()` calls define a graphite-dark palette (background `#242429`, accent `#3b82f6`).
- **MockPlotWidget:** A `QFrame` subclass that paints a grid and title — placeholder for future chart rendering.
- **ClassificationDisplay:** A read-only `QLineEdit` with `Qt::NoFocus` that displays the current classification as "N — label". Installed as an `eventFilter` on `QApplication` to intercept 0-9 key presses globally when no text field has focus.
- **Static linking:** Build produces a fully static musl binary, optionally UPX-compressed.
- **AUTOMOC:** CMake `CMAKE_AUTOMOC ON` is enabled for Qt Meta-Object Compiler support (required for `Q_OBJECT` in dialog classes).

## Classification System

The main window has a split "Log Files" row with two `QGroupBox`es:
- **Statistics** (left): "Class stats" button → opens `ClassificationStatsDialog`
- **Classification** (right): "Customize labels" button, "Current:" label, and read-only `ClassificationDisplay` box. Approving by clicking Enter (standard or NumPad) copies the classification to the "Type" input box.

### ClassificationDisplay (`src/main.cpp`)
- Read-only `QLineEdit` with `Qt::NoFocus` policy (never steals focus)
- Installed as `eventFilter` on `QApplication` — intercepts 0-9 key presses globally to change classification, and Enter/Return key presses to approve the classification (copying text to "Type" box)
- Only catches keys when the focused widget is NOT a `QLineEdit` (so typing in text fields works normally)
- Displays current classification as "N — label"

### CustomizeLabelsDialog (`src/windows/customize_labels.h` / `.cpp`)
- QDialog with 10 editable label rows (0-9), all editable
- Defaults: 0="nonvar", 1="var", 2-9="unknown"
- **NumPad navigation checkbox** at top (default off): when enabled, hides slots 2, 4, 6, 8 (numpad arrow keys), leaving only 0, 1, 3, 5, 7, 9. Visible rows are compacted to the top of the grid (no gaps).
- Checking (but not unchecking) this option dynamically closes and automatically reopens the dialog (preserving user label edits) to clean up layout presentation.
- Emits `labelsChanged()` signal on OK (or when automatically accepted for NumPad navigation change), which triggers `saveConfig()` and display refresh

### ClassificationStatsDialog (`src/windows/classification_stats.h` / `.cpp`)
- QDialog showing a table: number (0-9) | label text | count (all 0 for now)

## Configuration Persistence

On startup, `loadConfig()` reads `.lc-config.json` from the working directory using `json.h` (sheredom's header-only JSON parser from `../../include/json.h`). The config file uses key-value pairs:

```json
{
  "label0": "nonvar",
  "label1": "variable",
  "label2": "variable",
  ...
  "label9": "variable",
  "numpad_nav": "false"
}
```

`saveConfig()` is called whenever the Customize Labels dialog is accepted.

## Building

### Prerequisites

- Docker (required — the entire build runs inside containers)
- The build downloads a prebuilt Buildroot SDK (musl + Qt6 static) from `https://github.com/fifbroman8/buildroot/releases/download/v1.1.0/`

### Build

```sh
make          # Full build: Docker image + static compile + strip + UPX compress
```

This runs `build.sh`, which:
1. Builds the `lc-qt-builder` Docker image (Debian slim + build-essential, cmake, upx-ucl)
2. Downloads and extracts the Buildroot musl SDK into `buildroot/` (if not already present)
3. Relocates the SDK via `buildroot/relocate-sdk.sh`
4. Removes `.so` files from sysroot to enforce static linking
5. Mounts `../../include` as `/external_include` inside Docker (for `json.h`)
6. Runs CMake with the Buildroot toolchain file and `-DEXTERNAL_INCLUDE_DIR=/external_include`
7. Strips and UPX-compresses the final binary, output as `./lc-qt`

Other Make targets:
```sh
make clean         # Remove build/ directory
make deep-clean    # Remove build/ and buildroot/ SDK cache
make docker-clean  # Remove builder containers and image
```

### Manual build (inside Docker)

```sh
docker build -t lc-qt-builder .
docker run --rm -v "$(pwd)":/app -v "$(pwd)/../../include":/external_include lc-qt-builder bash -c "
    ./buildroot/relocate-sdk.sh
    find buildroot/x86_64-buildroot-linux-musl/sysroot -name '*.so' -delete
    mkdir -p build && cd build
    cmake .. -DCMAKE_TOOLCHAIN_FILE=/app/buildroot/share/buildroot/toolchainfile.cmake -DEXTERNAL_INCLUDE_DIR=/external_include
    make -j\$(nproc)
"
```

## Running

```sh
./lc-qt
```

Displays the Qt6 GUI window. Press 0-9 keys (when no text field is focused) to set classification. Edit labels via Customize Labels dialog (persists to `.lc-config.json`).

## File Structure

```
lc-qt/
├── src/
│   ├── main.cpp                    — Qt6 application entry point: UI layout, MockPlotWidget,
│   │                                  ClassificationDisplay, loadConfig/saveConfig, main()
│   └── windows/
│       ├── customize_labels.h      — CustomizeLabelsDialog declaration
│       ├── customize_labels.cpp    — CustomizeLabelsDialog: 10 editable labels, NumPad nav toggle
│       ├── classification_stats.h  — ClassificationStatsDialog declaration
│       └── classification_stats.cpp — ClassificationStatsDialog: count-per-label table
├── CMakeLists.txt                  — CMake config (C++17, AUTOMOC, static linking, Qt6 Widgets,
│                                       external include path for json.h)
├── Makefile                        — Orchestrates Docker-based static build + UPX compression
├── build.sh                        — Docker build script: image build, SDK download, cmake compile
├── Dockerfile                      — Debian slim with build-essential, cmake, upx-ucl
├── README.md                       — Project description (work in progress)
├── QWEN.md                         — This file
├── lc-qt                           — Pre-built output binary (git-tracked)
├── .lc-config.json                 — Runtime config (labels + numpad_nav, created on first save)
└── buildroot/                      — Prebuilt Buildroot musl SDK (downloaded, not committed)
    ├── relocate-sdk.sh             — Toolchain relocation script
    ├── bin/                        — Cross-compiler binaries
    ├── share/                      — CMake toolchain file
    └── x86_64-buildroot-linux-musl/
        └── sysroot/                — Static Qt6 libraries + musl libc
```

## Development Conventions

- **Language:** C++17 (`CMAKE_CXX_STANDARD 17`, required).
- **Build system:** CMake (invoked via Docker with Buildroot toolchain). Makefile is a thin wrapper for the Docker workflow.
- **Static linking:** The final binary must be fully static (`-static`, `--gc-sections`, `.so` files stripped from sysroot). No dynamic linking.
- **Qt6 only:** Uses `Qt6::Widgets`. No Qt Quick/QML.
- **AUTOMOC:** Required for any class using `Q_OBJECT` — already enabled in CMakeLists.txt.
- **Style:** Dark-themed UI via inline stylesheets — no external CSS/QSS files.
- **Headers as implementation units:** Dialog windows use separate `.h`/`.cpp` pairs in `src/windows/`.
- **No test framework:** No tests exist yet.
- **The `buildroot/` directory** is downloaded by the build process and should not be committed to git (it is a large prebuilt SDK). Do not use `make deep-clean` to remove it — it is already gitignored.
- **json.h dependency:** The sheredom `json.h` header-only parser lives in the parent project at `ihsnpeaks/include/json.h`. It is mounted into Docker builds at `/external_include/` via `build.sh`. The CMake `EXTERNAL_INCLUDE_DIR` variable controls the path (defaults to `../../include` for local builds, overridden to `/external_include` in Docker).
