# lc-qt — Project Guide

## Project Overview

`lc-qt` is a clean-room rewrite of Prof. Grzegorz Pojmański's **lc** light-curve analysis tool (1997–2018), reimplemented with a Qt6 GUI frontend and **ihsnpeaks** as the backend periodogram engine.

The application renders a dark-themed UI replicating the original `lc` tool's layout: raw/phased light-curve plots, object list management, classification system with customizable labels, period search (AoV, IHS, GB, BLS), zoomed spectrum inspection with MMB peak finding, period scrolling / keyboard modification, and phased model overlays. The **Raw Light Curve** plot loads and displays `.dat`/`.fits` files using upstream ihsnpeaks C readout functions, with detrending and proper time-axis preservation. The **Period Search** (Negative Log-Likelihood spectrum) plot computes and displays IHS, AoV, GB, and BLS periodograms using the upstream ihsnpeaks computation pipeline via a multithreaded C bridge.

**Parent project:** This is a downstream subproject of `ihsnpeaks` — located at `ihsnpeaks/downstream/lc-qt/`. It uses `json.h` from the parent project (`ihsnpeaks/include/json.h`) for configuration persistence, and the upstream `src/` headers (readout, process, convolution, aov, nufft1, kthread) for the computation backend.

## Architecture

- **Qt6 Widgets:** Uses `QMainWindow`, `QFrame`, `QGroupBox`, `QLineEdit`, `QPushButton`, `QCheckBox` — classic widget-based desktop UI, no QML.
- **Dark theme:** Inline `setStyleSheet()` calls define a graphite-dark palette (background `#242429`, accent `#3b82f6`).
- **Main Window** (`src/windows/main_window.h/.cpp`): Main window implementation managing UI layout, toolbar, menus, light curve loading, classification persistence, and UI-worker signal routing.
- **LightCurvePlotWidget** (`src/windows/lightcurve_plot.h/.cpp`): A `QFrame` subclass that renders actual scatter data with auto-scaling, grid, and smart axis tick formatting.
- **PhasedLightCurveWidget** (`src/windows/phased_lightcurve.h/.cpp`): Renders phase-folded light curves duplicated over phase `[0, 2]` with analytical or smoothed model curves overlaid (`lc_model.c`).
- **SpectrumPlotWidget** (`src/windows/spectrum_plot.h/.cpp`): A `QFrame` subclass that renders the NLL periodogram as a line plot with peak-preserving downsampling (≤8 pts/pixel), not-a-knot cubic spline interpolation when sparse, 5% y-margin (starting at 0), and a progress percentage overlay in the upper-right corner during computation.
- **ZoomedSpectrumWidget** (`src/windows/zoomed_spectrum.h/.cpp`): High-resolution zoom window around the selected peak / pivot frequency with interactive FOV scaling and MMB local peak picking.
- **Computation Workers** (`src/period_worker.h/.cpp`):
  - `PeriodogramTask` (`QObject` on a `QThread`): Runs full multithreaded periodogram computation in the background via `lc_compute_periodogram_ctx`.
  - `PhasedModelWorker` (`QObject` on a dedicated `QThread`): Recomputes analytical/smoothed phased models asynchronously as the user scrolls the pivot frequency.
- **C bridge** (`src/lc_period.h` / `src/lc_period.c`): A single C translation unit (compiled with `-std=gnu23`) that wraps the upstream ihsnpeaks computation pipeline (`kthread.h`, `process.h`, `readout.h`, `convolution.h`, `aov.h`, `nufft1.h`, etc.) and exposes a clean C API to the C++ frontend: file I/O (`lc_load_dat`, `lc_load_fits`, `lc_detrend`, `lc_free`), periodogram computation (`lc_compute_periodogram_ctx`), phased model evaluation (`lc_compute_phased_model`, `lc_compute_phase_offset`), and progress monitoring (`lc_progress_*`).
- **Phased Model Evaluation** (`src/lc_model.c`): Single-frequency model generator for IHS, AoV (Szego orthogonal polynomials), GBLS (convolution smoothing), and BLS (boxcar fit), included into `lc_period.c`.
- **Persistent computation context** (`lc_compute_ctx_t`): Caches the thread pool, NuFFT plans, and worker buffers between runs. Rebuilds only when data size, method, fmax, grid mode, or nterms change.
- **CustomizePeriodSearchDialog** (`src/windows/customize_period_search.h/.cpp`): QDialog with period search parameters: harmonics, oversampling, min/max frequency, statistic (`raw` vs `bayes`), zoom factor, search radius, upsampling factor (4/3 or 2/1), oversmoothing factor (GBLS/BLS), and number of bins (BLS).
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
- **Object List** (left): Input box, Open button, and Entry No. field
- **Classification** (right): "Customize labels" button, "Class stats" button, "Cur:" label, and read-only `ClassificationDisplay` box. Approving by clicking Enter (standard or NumPad) copies the classification to the "Type" input box.

## Period Search (Spectrum Generation)

### Overview

The "Period Search" section computes and displays Negative Log-Likelihood (NLL) / Bayesian Relative Evidence frequency spectra using the upstream ihsnpeaks computation pipeline. Four methods are available: **AoV**, **IHS**, **GB** (GBLS), and **BLS**. Computation is fully multithreaded using a shared `kt_forpool_t` thread pool.

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
- **Progress**: NuFFT-block-level atomic progress updates for IHS/AoV; per-frequency atomic updates for GB/BLS. Progress percentage is drawn in the upper-right corner of the spectrum plot. Time estimation is shown in the status bar.
- **Cancellation**: Stop button sets an atomic cancel flag. GB/BLS workers stop at chunk boundaries; IHS/AoV run to completion (fast) and discard the result.

### Computation flow

1. `PeriodogramTask` (QObject on a QThread) builds `lc_periodogram_config_t` from dialog settings.
2. `lc_compute_periodogram_ctx()` uses the persistent context (`lc_compute_ctx_t`):
   - Rebuilds cached resources (plans, buffers) only when data size, method, fmax, grid mode, or nterms change.
   - Copies data into the primary buffer, runs `preprocess_buffer`.
   - Computes frequency grid using the **actual time span** (`x[n-1] - x[0]`), not absolute `x[n-1]`.
   - Runs the parallel sweep (IHS/AoV via NuFFT blocks, GB/BLS via direct grid).
   - Converts power→NLL in parallel (Pochhammer rising factorial expansions for AoV/GB, `correct_ihs_res` for IHS, or Bayesian evidence).
3. Result (`fmin`, `fstep`, `nll[]`, detected peaks) is emitted to the GUI thread, rendered by `SpectrumPlotWidget`, and passed to `ZoomedSpectrumWidget`.

### Customize Period Search dialog

Accessible via the "…" button or the "Customize Period Search..." menu action. Parameters:
- **Number of harmonics** (`nterms`, default 3)
- **Oversampling** (default 5)
- **Min frequency** (QLineEdit with "auto" placeholder; empty = auto `2/span`)
- **Max frequency** (default 24)
- **Statistic** (`raw` for classical NLL, `bayes` for Bayesian evidence ratio)
- **Zooming factor** (default 4)
- **Search radius** (default 0.1)
- **Upsampling factor** (combo: "4/3" → pswf43, "2/1" → pswf21)
- **Oversmoothing factor** (default 0.2 → `gbAlpha` for GBLS, `blsMinRelWidth` for BLS; BLS max width hardcoded 0.5)
- **Number of bins** (default 10 → `blsWidthCount`, BLS only)

All settings persist to `.lc-config.json`. Decimal delimiter is forced to "." (`QLocale::c()`).

## Configuration Persistence

On startup, `loadConfig()` reads `.lc-config.json` from the working directory using `json.h` (sheredom's header-only JSON parser from `../../include/json.h`).

```json
{
  "labels": ["nonvar", "var", "unknown", ...],
  "numpad_nav": "false",
  "files": [
    {"path": "/path/to/OGLE-BLAP-035.dat", "n": 2240}
  ],
  "targets": [
    {"path": "/path/to/OGLE-BLAP-035.dat", "class": 1}
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
    "nbins": 10,
    "scroll_rate": 1.0,
    "auto_center": "true",
    "display_frequency": "false",
    "statistic": "bayes"
  }
}
```

`saveConfig()` writes the array format back whenever settings or classifications change.

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

## File Structure

```
lc-qt/
├── src/
│   ├── main.cpp                    — Application entry point: QApplication setup, main window launch
│   ├── period_worker.h/.cpp        — Multithreaded PeriodogramTask & PhasedModelWorker QObject workers
│   ├── lc_readout.h                — C bridge public API: lc_data_t, lc_load_dat/fits/detrend/free
│   ├── lc_period.h                 — C bridge public API: lc_periodogram_config_t, lc_compute_ctx_t,
│   │                                  lc_compute_periodogram_ctx, lc_compute_phased_model, lc_progress_*
│   ├── lc_period.c                 — C bridge implementation (gnu23): wraps upstream kthread.h +
│   │                                  process.h + readout.h; multithreaded IHS/AoV/GB/BLS sweeps,
│   │                                  Pochhammer rising factorials / Bayesian conversions, persistent context
│   ├── lc_model.c                  — Phased model curve generation (Szego polynomials, convolution, boxcar),
│   │                                  included into lc_period.c
│   └── windows/
│       ├── main_window.h/.cpp      — MainWindow class: layout, menus, wiring, key filters, config I/O
│       ├── lightcurve_plot.h/.cpp  — Scatter plot renderer with auto-scaling and smart ticks
│       ├── phased_lightcurve.h/.cpp— Phased scatter plot with model fit overlay
│       ├── spectrum_plot.h/.cpp    — NLL spectrum renderer: peak-preserving downsample, cubic spline
│       ├── zoomed_spectrum.h/.cpp  — Interactive high-res FOV zoom around selected pivot frequency
│       ├── customize_period_search.h/.cpp — Dialog: harmonics, oversampling, fmin/fmax, pswf, GBLS/BLS, stats
│       ├── customize_labels.h/.cpp — CustomizeLabelsDialog: 10 editable labels, NumPad nav toggle
│       ├── classification_stats.h/.cpp — ClassificationStatsDialog: count-per-label table
│       └── period_scroll.h/.cpp    — PeriodScrollDialog: scroll rate & fine tuning
├── test_phased_verify/             — Standalone verification test harness for phased light curve rendering
│   ├── CMakeLists.txt              — Test harness build config
│   └── main_verify.cpp             — Automated verification test suite
├── CMakeLists.txt                  — CMake config (C23 + C++17, AUTOMOC, static linking, Qt6 Widgets,
│                                       lc_backend OBJECT library, build-time scaling.h generation)
├── Makefile                        — Orchestrates Docker-based static build + UPX compression
├── build.sh                        — Docker build script: image build, SDK download, cmake compile,
│                                       /src symlink for kthread.h cross-mount include
├── Dockerfile                      — Debian slim with build-essential, cmake, upx-ucl
├── README.md                       — Project description (work in progress)
├── AGENTS.md                       — This file
├── lc-qt                           — Pre-built output binary (git-tracked)
├── .lc-config.json                 — Runtime config (labels, files, period_search settings)
└── buildroot/                      — Prebuilt Buildroot musl SDK (downloaded, not committed)
```

## Development Conventions

- **Languages:** C23 (`-std=gnu23` for the upstream C bridge) and C++17 (`CMAKE_CXX_STANDARD 17`).
- **Build system:** CMake (invoked via Docker with Buildroot toolchain). Makefile is a thin wrapper for the Docker workflow.
- **Static linking:** The final binary must be fully static (`-static`, `--gc-sections`, `.so` files stripped from sysroot). No dynamic linking.
- **Qt6 only:** Uses `Qt6::Widgets`. No Qt Quick/QML.
- **AUTOMOC:** Required for any class using `Q_OBJECT` — enabled in CMakeLists.txt.
- **Backend OBJECT library:** `lc_backend` holds `lc_period.c` + `${UPSTREAM_SRC_DIR}/nufft/nufft1.c`. Both compiled with `-std=gnu23 -D_GNU_SOURCE -DMAX_TWIDDLE_REUSE=8`. Dispatch-ready: parameterized `LC_BACKEND_MARCH` for future multi-microarchitecture support.
- **Build-time scaling.h:** Generated by compiling+running upstream `scaling.c` + `nufft1.c` as `scaling_gen` (static musl binary runs inside the build container). Emits tuned PSWF constants into `${CMAKE_BINARY_DIR}/scaling.h`.
- **Single C translation unit:** All upstream headers (`kthread.h`, `process.h`, `readout.h`, `convolution.h`, `aov.h`, `nufft1.h`, header-only libs) are included ONLY in `lc_period.c` (and its `#include "lc_model.c"`). Non-static definitions in these headers would cause duplicate-symbol link errors if included elsewhere.
- **Upstream include isolation:** The upstream `src/` include paths are applied ONLY to the backend sources via `set_source_files_properties` / OBJECT library compile options — never via global `target_include_directories`.
- **Style:** Dark-themed UI via inline stylesheets — no external CSS/QSS files.
- **The `buildroot/` directory** is downloaded by the build process and should not be committed to git (it is a large prebuilt SDK). Do not use `make deep-clean` to remove it — it is already gitignored.
- **json.h dependency:** The sheredom `json.h` header-only parser lives in the parent project at `ihsnpeaks/include/json.h`. It is mounted into Docker builds at `/external_include/` via `build.sh`. The CMake `EXTERNAL_INCLUDE_DIR` variable controls the path.
