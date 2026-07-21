# lc-qt — Project Guide

## Project Overview

`lc-qt` is a clean-room rewrite of Prof. Grzegorz Pojmański's **lc** light-curve analysis tool (1997–2018), reimplemented with a Qt6 GUI frontend and **ihsnpeaks** as the backend periodogram engine. The project is a work in progress — currently closer to a functional mock-up than a full lc replacement.

The application renders a dark-themed UI replicating the original `lc` tool's layout: raw/phased light-curve plots, object list management, classification system with customizable labels, period search (AoV, IHS, GB, BLS), and period modify controls. The **Raw Light Curve** plot is functional — it loads and displays `.dat`/`.fits` files using upstream ihsnpeaks C readout functions, with detrending and proper time-axis preservation. The **Period Search** (Negative Log-Likelihood spectrum) plot is functional — it computes and displays IHS, AoV, GB, and BLS periodograms using the upstream ihsnpeaks computation pipeline via a multithreaded C bridge. All other plot widgets remain mock placeholders.

**Parent project:** This is a downstream subproject of `ihsnpeaks` — located at `ihsnpeaks/downstream/lc-qt/`. It uses `json.h` from the parent project (`ihsnpeaks/include/json.h`) for configuration persistence, and the upstream `src/` headers (readout, process, convolution, aov, nufft1, kthread) for the computation backend.

## Architecture

- **Qt6 Widgets:** Uses `QMainWindow`, `QFrame`, `QGroupBox`, `QLineEdit`, `QPushButton`, `QCheckBox` — classic widget-based desktop UI, no QML.
- **Dark theme:** Inline `setStyleSheet()` calls define a graphite-dark palette (background `#242429`, accent `#3b82f6`).
- **LightCurvePlotWidget** (`src/windows/lightcurve_plot.h/.cpp`): A `QFrame` subclass that renders actual scatter data with auto-scaling, grid, and smart axis tick formatting.
- **SpectrumPlotWidget** (`src/windows/spectrum_plot.h/.cpp`): A `QFrame` subclass that renders the NLL periodogram as a line plot with peak-preserving downsampling (≤8 pts/pixel), not-a-knot cubic spline interpolation when sparse, 5% y-margin (starting at 0), and a progress percentage overlay in the upper-right corner during computation.
- **C bridge** (`src/lc_periodogram.h` / `src/lc_periodogram.c`): A single C translation unit (compiled with `-std=gnu23`) that wraps the entire upstream ihsnpeaks computation pipeline (`kthread.h`, `process.h`, `readout.h`, `convolution.h`, `aov.h`, `nufft1.h`, etc.) and exposes a clean C API to the C++ frontend: file I/O (`lc_load_dat`, `lc_load_fits`, `lc_detrend`, `lc_free`), periodogram computation (`lc_compute_periodogram_ctx`), and progress monitoring (`lc_progress_*`).
- **Persistent computation context** (`lc_compute_ctx_t`): Caches the thread pool, NuFFT plans, and worker buffers between runs. Rebuilds only when data size, method, fmax, grid mode, or nterms change.
- **CustomizePeriodSearchDialog** (`src/windows/customize_period_search.h/.cpp`): QDialog with period search parameters: harmonics, oversampling, min/max frequency, zoom factor, search radius, upsampling factor (4/3 or 2/1), oversmoothing factor (GBLS/BLS), and number of bins (BLS).
- **ClassificationDisplay:** A read-only `QLineEdit` with `Qt::NoFocus` that displays the current classification as "N — label". Installed as an `eventFilter` on `QApplication` to intercept 0-9 key presses globally when no text field or spin box has focus.
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

## Period Search (Spectrum Generation)

### Overview

The "Period Search" section computes and displays Negative Log-Likelihood (NLL) frequency spectra using the upstream ihsnpeaks computation pipeline. Four methods are available: **AoV**, **IHS**, **GB** (GBLS), and **BLS**. Computation is fully multithreaded using a shared `kt_forpool_t` thread pool.

### Button mapping

| Button | periodogramMethod | gbEvalMode | mode | fstep formula |
|--------|-------------------|------------|------|---------------|
| IHS | PERIODOGRAM_IHS | — | 2 | `1/(nterms × oversampling × span × 0.5)` |
| AoV | PERIODOGRAM_AOV | — | 2 | same as IHS |
| GB | — | GB_EVAL_GBLS | 5 | `gb_direct_frequency_step(n, span, params)` |
| BLS | — | GB_EVAL_BLS | 5 | `bls_direct_frequency_step(span, params)` |

### Threading model

- **Thread pool**: Created with all logical CPUs (SMT enabled) via `sysconf(_SC_NPROCESSORS_ONLN)`.
- **GB/BLS**: Use all logical threads (embarrassingly parallel per-frequency evaluation with independent scratch spaces).
- **IHS/AoV**: Workset count capped at physical core count (`max_slices` parameter) to avoid SMT contention on NuFFT-heavy workloads.
- **Progress**: NuFFT-block-level atomic progress updates for IHS/AoV; per-frequency atomic updates for GB/BLS. Progress percentage is drawn in the upper-right corner of the spectrum plot (font 12 bold). Time estimation is shown in the status bar.
- **Cancellation**: Stop button sets an atomic cancel flag. GB/BLS workers stop at chunk boundaries; IHS/AoV run to completion (fast) and discard the result.

### Computation flow

1. `PeriodogramTask` (QObject on a QThread) builds `lc_periodogram_config_t` from dialog settings.
2. `lc_compute_periodogram_ctx()` uses the persistent context (`lc_compute_ctx_t`):
   - Rebuilds cached resources (plans, buffers) only when data size, method, fmax, grid mode, or nterms change.
   - Copies data into the primary buffer, runs `preprocess_buffer`.
   - Computes frequency grid using the **actual time span** (`x[n-1] - x[0]`), not absolute `x[n-1]`.
   - Runs the parallel sweep (IHS/AoV via NuFFT blocks, GB/BLS via direct grid).
   - Converts power→NLL in parallel (`aov_likelihood_from_r2` for AoV, `correct_ihs_res` for IHS, direct for GB/BLS).
3. Result (`fmin`, `fstep`, `nll[]`) is emitted to the GUI thread and rendered by `SpectrumPlotWidget`.

### Customize Period Search dialog

Accessible via the "…" button or the "Customize Period Search..." menu action. Parameters:
- **Number of harmonics** (`nterms`, default 3)
- **Oversampling** (default 5)
- **Min frequency** (QLineEdit with "auto" placeholder; empty = auto `2/span`)
- **Max frequency** (default 24)
- **Zooming factor** (default 4, stored but unused)
- **Search radius** (default 0.1, stored but unused)
- **Upsampling factor** (combo: "4/3" → pswf43, "2/1" → pswf21)
- **Oversmoothing factor** (default 0.2 → `gbAlpha` for GBLS, `blsMinRelWidth` for BLS; BLS max width hardcoded 0.5)
- **Number of bins** (default 10 → `blsWidthCount`, BLS only)

All settings persist to `.lc-config.json`. Decimal delimiter is forced to "." (`QLocale::c()`).

### Spectrum rendering

- **Y-axis**: Starts at 0, 5% headroom above max, not inverted (higher NLL = higher on screen).
- **Dense data** (≥1 pt/px): Peak-preserving downsample to ≤8 pts/pixel (keeps local maxima).
- **Sparse data** (<0.75 pt/px): Not-a-knot cubic spline interpolation at ~6 samples/px.
- **Empty state**: Shows "No data" centered.

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
  ],
  "period_search": {
    "nterms": 3,
    "oversampling": 5.0,
    "fmin": 0.0,
    "fmax": 24.0,
    "zoom_factor": 4.0,
    "search_radius": 0.1,
    "pswf": 43,
    "oversmoothing": 0.2,
    "nbins": 10
  }
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
│   │                                  ClassificationDisplay, PeriodogramTask worker, period search
│   │                                  wiring, loadConfig/saveConfig, main()
│   ├── lc_readout.h                — C bridge public API: lc_data_t, lc_load_dat/fits/detrend/free
│   ├── lc_periodogram.h            — C bridge public API: lc_periodogram_config_t, lc_compute_ctx_t,
│   │                                  lc_compute_periodogram_ctx, lc_progress_*
│   ├── lc_periodogram.c            — C bridge implementation (gnu23): wraps upstream kthread.h +
│   │                                  process.h + readout.h; multithreaded IHS/AoV/GB/BLS sweeps,
│   │                                  parallel NLL conversion, persistent context, physical core count
│   └── windows/
│       ├── lightcurve_plot.h       — LightCurvePlotWidget declaration (Q_OBJECT)
│       ├── lightcurve_plot.cpp     — Scatter plot renderer with auto-scaling and smart ticks
│       ├── spectrum_plot.h         — SpectrumPlotWidget declaration (line plot + progress overlay)
│       ├── spectrum_plot.cpp       — NLL spectrum renderer: peak-preserving downsample, cubic spline
│       ├── customize_period_search.h  — CustomizePeriodSearchDialog + PeriodSearchSettings struct
│       ├── customize_period_search.cpp — Dialog: harmonics, oversampling, fmin/fmax, pswf, GBLS/BLS params
│       ├── customize_labels.h      — CustomizeLabelsDialog declaration
│       ├── customize_labels.cpp    — CustomizeLabelsDialog: 10 editable labels, NumPad nav toggle
│       ├── classification_stats.h  — ClassificationStatsDialog declaration
│       └── classification_stats.cpp — ClassificationStatsDialog: count-per-label table
├── CMakeLists.txt                  — CMake config (C23 + C++17, AUTOMOC, static linking, Qt6 Widgets,
│                                       lc_backend OBJECT library, build-time scaling.h generation)
├── Makefile                        — Orchestrates Docker-based static build + UPX compression
├── build.sh                        — Docker build script: image build, SDK download, cmake compile,
│                                       /src symlink for kthread.h cross-mount include
├── Dockerfile                      — Debian slim with build-essential, cmake, upx-ucl
├── README.md                       — Project description (work in progress)
├── QWEN.md                         — This file
├── lc-qt                           — Pre-built output binary (git-tracked)
├── .lc-config.json                 — Runtime config (labels, files, period_search settings)
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
- **AUTOMOC:** Required for any class using `Q_OBJECT` — already enabled in CMakeLists.txt. `main.cpp` includes `main.moc` at the end for the in-file `PeriodogramTask` Q_OBJECT class.
- **Backend OBJECT library:** `lc_backend` holds `lc_periodogram.c` + `${UPSTREAM_SRC_DIR}/nufft/nufft1.c`. Both compiled with `-std=gnu23 -D_GNU_SOURCE -DMAX_TWIDDLE_REUSE=8`. Dispatch-ready: parameterized `LC_BACKEND_MARCH` for future multi-microarchitecture support.
- **Build-time scaling.h:** Generated by compiling+running upstream `scaling.c` + `nufft1.c` as `scaling_gen` (static musl binary runs inside the build container). Emits tuned PSWF constants into `${CMAKE_BINARY_DIR}/scaling.h`.
- **Single C translation unit:** All upstream headers (`kthread.h`, `process.h`, `readout.h`, `convolution.h`, `aov.h`, `nufft1.h`, header-only libs) are included ONLY in `lc_periodogram.c`. Non-static definitions in these headers would cause duplicate-symbol link errors if included elsewhere.
- **Upstream include isolation:** The upstream `src/` include paths are applied ONLY to the backend sources via `set_source_files_properties` / OBJECT library compile options — never via global `target_include_directories`.
- **Style:** Dark-themed UI via inline stylesheets — no external CSS/QSS files.
- **Headers as implementation units:** Dialog windows and widgets use separate `.h`/`.cpp` pairs in `src/windows/`.
- **No test framework:** No tests exist yet.
- **The `buildroot/` directory** is downloaded by the build process and should not be committed to git (it is a large prebuilt SDK). Do not use `make deep-clean` to remove it — it is already gitignored.
- **json.h dependency:** The sheredom `json.h` header-only parser lives in the parent project at `ihsnpeaks/include/json.h`. It is mounted into Docker builds at `/external_include/` via `build.sh`. The CMake `EXTERNAL_INCLUDE_DIR` variable controls the path.
- **build.sh symlink:** Creates `/src -> /upstream_src` inside the container so `kthread.h`'s relative include (`../../src/params.h`) resolves across the separate Docker mounts.
