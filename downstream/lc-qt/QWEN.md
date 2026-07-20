# lc-qt — Project Guide

## Project Overview

`lc-qt` is a clean-room rewrite of Prof. Grzegorz Pojmański's **lc** light-curve analysis tool (1997–2018), reimplemented with a Qt6 GUI frontend and **ihsnpeaks** as the backend periodogram engine. The project is a work in progress — currently closer to a functional mock-up than a full lc replacement.

The application renders a dark-themed UI replicating the original `lc` tool's layout: raw/phased light-curve plots, object list management, classification system with customizable labels, period search (AoV, IHS, GB, BLS), and period modify controls. The **Raw Light Curve** plot is functional — it loads and displays `.dat`/`.fits` files using upstream ihsnpeaks C readout functions, with detrending and proper time-axis preservation. All other plot widgets remain mock placeholders.

**Parent project:** This is a downstream subproject of `ihsnpeaks` — located at `ihsnpeaks/downstream/lc-qt/`. It uses `json.h` from the parent project (`ihsnpeaks/include/json.h`) for configuration persistence, and the upstream `src/utils/readout.h` (plus its dependency chain) for light-curve file I/O.

## Architecture

- **Qt6 Widgets:** Uses `QMainWindow`, `QFrame`, `QGroupBox`, `QLineEdit`, `QPushButton`, `QCheckBox` — classic widget-based desktop UI, no QML.
- **Dark theme:** Inline `setStyleSheet()` calls define a graphite-dark palette (background `#242429`, accent `#3b82f6`).
- **LightCurvePlotWidget** (`src/windows/lightcurve_plot.h/.cpp`): A `QFrame` subclass that renders actual scatter data with auto-scaling, grid, and smart axis tick formatting.
- **MockPlotWidget:** A `QFrame` subclass that paints a grid and title — placeholder for phased plot, period search, and period modify charts.
- **C bridge** (`src/lc_readout.h` / `src/lc_readout.c`): A single C translation unit (compiled with `-std=gnu23`) that wraps upstream ihsnpeaks header-only libraries (`readout.h`, `fast_convert.h`, `qfits.h`, `sds.h`, `fdist.h`) and exposes a clean C API (`lc_load_dat`, `lc_load_fits`, `lc_detrend`, `lc_free`) to the C++ frontend.
- **ClassificationDisplay:** A read-only `QLineEdit` with `Qt::NoFocus` that displays the current classification as "N — label". Installed as an `eventFilter` on `QApplication` to intercept 0-9 key presses globally when no text field has focus.
- **Static linking:** Build produces a fully static musl binary, optionally UPX-compressed.
- **AUTOMOC:** CMake `CMAKE_AUTOMOC ON` is enabled for Qt Meta-Object Compiler support (required for `Q_OBJECT` in dialog and widget classes).

## Light Curve Loading

### Supported formats
- **`.dat`** — whitespace-delimited text (time, magnitude, error per row), loaded via upstream `read_dat`
- **`.fits`** — FITS binary table (TESS/OGLE style), loaded via upstream `read_fits`

### Loading flow
1. User selects a file via **File → Load Light Curve...** dialog, or a path is passed as `argv[1]` (auto-loads on startup).
2. The C bridge stats the file and allocates exactly-sized buffers (newline count for `.dat`, table row count for `.fits`).
3. After reading, `lc_detrend` performs weighted linear regression detrending (via upstream `linregw_buffer`), then restores the original time offset and adds the mean magnitude back.
4. The detrended data (with original timescale) is displayed in the Raw Light Curve plot.

### Display
- Y-axis tick precision: `floor(-log10(ySpan / 100))` decimal places (resolves 1/100 of visible amplitude).
- X-axis tick precision: integers when values exceed 1000, otherwise up to 3 decimal places.
- No y-axis unit label; x-axis labeled "Time [d]".

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

On startup, `loadConfig()` reads `.lc-config.json` from the working directory using `json.h` (sheredom's header-only JSON parser from `../../include/json.h`). The config uses a JSON array format:

```json
{
  "labels": ["nonvar", "var", "unknown", ...],
  "numpad_nav": "false",
  "files": [
    {"path": "/path/to/OGLE-BLAP-035.dat", "n": 2240}
  ]
}
```

Legacy `"label0"`...`"label9"` key-value format is still read for backward compatibility. `saveConfig()` always writes the new array format.

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
5. Mounts `../../include` as `/external_include` and `../../src` as `/upstream_src` inside Docker
6. Runs CMake with the Buildroot toolchain file, `-DEXTERNAL_INCLUDE_DIR=/external_include`, and `-DUPSTREAM_SRC_DIR=/upstream_src`
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
docker run --rm -v "$(pwd)":/app -v "$(pwd)/../../include":/external_include -v "$(pwd)/../../src":/upstream_src lc-qt-builder bash -c "
    ./buildroot/relocate-sdk.sh
    find buildroot/x86_64-buildroot-linux-musl/sysroot -name '*.so' -delete
    mkdir -p build && cd build
    cmake .. -DCMAKE_TOOLCHAIN_FILE=/app/buildroot/share/buildroot/toolchainfile.cmake -DEXTERNAL_INCLUDE_DIR=/external_include -DUPSTREAM_SRC_DIR=/upstream_src
    make -j\$(nproc)
"
```

## Running

```sh
./lc-qt                  # Launch with empty plot
./lc-qt path/to/file.dat # Launch and auto-load a light curve
```

Displays the Qt6 GUI window. Press 0-9 keys (when no text field is focused) to set classification. Edit labels via Customize Labels dialog (persists to `.lc-config.json`). Load light curves via File menu or command-line argument.

## File Structure

```
lc-qt/
├── src/
│   ├── main.cpp                    — Qt6 application entry point: UI layout, light curve loading,
│   │                                  ClassificationDisplay, loadConfig/saveConfig, main()
│   ├── lc_readout.h                — C bridge public API (lc_data_t, lc_load_dat/fits/detrend/free)
│   ├── lc_readout.c                — C bridge implementation (gnu23, wraps upstream readout.h)
│   └── windows/
│       ├── lightcurve_plot.h       — LightCurvePlotWidget declaration (Q_OBJECT)
│       ├── lightcurve_plot.cpp     — Scatter plot renderer with auto-scaling and smart ticks
│       ├── customize_labels.h      — CustomizeLabelsDialog declaration
│       ├── customize_labels.cpp    — CustomizeLabelsDialog: 10 editable labels, NumPad nav toggle
│       ├── classification_stats.h  — ClassificationStatsDialog declaration
│       └── classification_stats.cpp — ClassificationStatsDialog: count-per-label table
├── CMakeLists.txt                  — CMake config (C23 + C++17, AUTOMOC, static linking, Qt6 Widgets,
│                                       upstream include paths isolated to C bridge only)
├── Makefile                        — Orchestrates Docker-based static build + UPX compression
├── build.sh                        — Docker build script: image build, SDK download, cmake compile
├── Dockerfile                      — Debian slim with build-essential, cmake, upx-ucl
├── README.md                       — Project description (work in progress)
├── QWEN.md                         — This file
├── lc-qt                           — Pre-built output binary (git-tracked)
├── .lc-config.json                 — Runtime config (labels array + files list + numpad_nav)
└── buildroot/                      — Prebuilt Buildroot musl SDK (downloaded, not committed)
    ├── relocate-sdk.sh             — Toolchain relocation script
    ├── bin/                        — Cross-compiler binaries
    ├── share/                      — CMake toolchain file
    └── x86_64-buildroot-linux-musl/
        └── sysroot/                — Static Qt6 libraries + musl libc
```

## Development Conventions

- **Languages:** C23 (`-std=gnu23` for the upstream C bridge) and C++17 (`CMAKE_CXX_STANDARD 17`).
- **Build system:** CMake (invoked via Docker with Buildroot toolchain). Makefile is a thin wrapper for the Docker workflow.
- **Static linking:** The final binary must be fully static (`-static`, `--gc-sections`, `.so` files stripped from sysroot). No dynamic linking.
- **Qt6 only:** Uses `Qt6::Widgets`. No Qt Quick/QML.
- **AUTOMOC:** Required for any class using `Q_OBJECT` — already enabled in CMakeLists.txt.
- **Upstream include isolation:** The upstream `src/` include paths are applied ONLY to `lc_readout.c` via `set_source_files_properties` COMPILE_FLAGS — never via `target_include_directories` — to prevent moc from parsing conflicting upstream headers (notably `strings.h`).
- **Style:** Dark-themed UI via inline stylesheets — no external CSS/QSS files.
- **Headers as implementation units:** Dialog windows and widgets use separate `.h`/`.cpp` pairs in `src/windows/`.
- **No test framework:** No tests exist yet.
- **The `buildroot/` directory** is downloaded by the build process and should not be committed to git (it is a large prebuilt SDK). Do not use `make deep-clean` to remove it — it is already gitignored.
- **json.h dependency:** The sheredom `json.h` header-only parser lives in the parent project at `ihsnpeaks/include/json.h`. It is mounted into Docker builds at `/external_include/` via `build.sh`. The CMake `EXTERNAL_INCLUDE_DIR` variable controls the path (defaults to `../../include` for local builds, overridden to `/external_include` in Docker).
- **Upstream C dependency:** The C bridge (`lc_readout.c`) compiles the entire upstream header-only dependency chain (`readout.h` → `common.h`, `fast_convert.h`, `qfits.h`, `sds.h`, `fdist.h`, `convolution.h`, `nufft1.h`) in a single translation unit. These libraries use non-static function definitions and must not be included elsewhere.
