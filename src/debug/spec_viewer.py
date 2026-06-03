#!/usr/bin/env python3.12

import sys
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
from scipy.interpolate import CubicSpline


# Activate spline mode when original visible samples become sparse
# relative to horizontal screen resolution.
OVERZOOM_POINTS_PER_PIXEL = 0.75

# Rendered spline samples per horizontal screen pixel.
SPLINE_SAMPLES_PER_PIXEL = 6.0

# Extra original points taken on both sides of the visible range when
# constructing the local spline. This reduces edge artifacts.
SPLINE_MARGIN_POINTS = 6


def sanitize_xy_columns(x, y):
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    if len(x) == 0:
        raise ValueError("no finite x-y pairs found")

    # Required for searchsorted and spline construction.
    if np.any(np.diff(x) < 0):
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]

    # CubicSpline requires strictly increasing x.
    unique = np.empty(len(x), dtype=bool)
    unique[0] = True
    unique[1:] = np.diff(x) > 0

    x = x[unique]
    y = y[unique]

    if len(x) < 2:
        raise ValueError("need at least two unique x-values")

    return x, y


def load_spectrum_file(path):
    data = np.loadtxt(path, dtype=np.float64)

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("expected at least two numeric columns")

    x = data[:, 0]
    columns = []

    for col_idx in range(1, data.shape[1]):
        y = data[:, col_idx]
        col_x, col_y = sanitize_xy_columns(x, y)
        columns.append((col_idx - 1, col_x, col_y))

    return columns


def visible_indices(x, xmin, xmax):
    i0 = np.searchsorted(x, xmin, side="left")
    i1 = np.searchsorted(x, xmax, side="right")

    i0 = max(0, min(i0, len(x)))
    i1 = max(0, min(i1, len(x)))

    return i0, i1


def evaluate_local_cubic_spline(x, y, xmin, xmax, width_px, i0, i1):
    """
    Build a local cubic spline on the visible data plus a small margin,
    then evaluate it densely only inside [xmin, xmax].

    This avoids building a global spline for the entire million-point file.
    """
    n = len(x)

    if i1 <= i0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    j0 = max(0, i0 - SPLINE_MARGIN_POINTS)
    j1 = min(n, i1 + SPLINE_MARGIN_POINTS)

    xv = x[j0:j1]
    yv = y[j0:j1]

    if len(xv) < 2:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    a = max(float(xmin), float(x[0]))
    b = min(float(xmax), float(x[-1]))

    if not np.isfinite(a) or not np.isfinite(b) or b <= a:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    # The default bc_type is "not-a-knot".
    # For local splines this is usually better than "natural", because
    # natural boundary conditions can artificially flatten the local edges.
    spline = CubicSpline(xv, yv, extrapolate=False)

    n_samples = max(64, int(width_px * SPLINE_SAMPLES_PER_PIXEL))

    xp = np.linspace(a, b, n_samples, dtype=np.float64)
    yp = spline(xp)

    finite = np.isfinite(xp) & np.isfinite(yp)
    xp = xp[finite]
    yp = np.maximum(yp[finite], 0.0)
    return xp, yp


def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} file.tsv [file2.tsv ...]")
        return 1

    dataset_groups = []

    for path in sys.argv[1:]:
        try:
            columns = load_spectrum_file(path)
        except Exception as exc:
            print(f"failed to load {path!r}: {exc}")
            return 1

        file_label = Path(path).name
        multi_column = len(columns) > 1
        group = []

        for iteration, x, y in columns:
            label = f"{file_label} iter {iteration}" if multi_column else file_label
            group.append(
                {
                    "path": path,
                    "label": label,
                    "x": x,
                    "y": y,
                    "y_max": float(np.nanmax(y)),
                    "curve": None,
                    "mode": None,
                }
            )
            print(f"loaded {label}: {len(x):,} unique finite points")

        group_y_max = max(dataset["y_max"] for dataset in group)
        dataset_groups.append((group_y_max, group))

    dataset_groups.sort(key=lambda item: item[0], reverse=True)
    datasets = [dataset for _, group in dataset_groups for dataset in group]

    app = QtWidgets.QApplication([])

    pg.setConfigOptions(antialias=False)

    title = " | ".join(dataset["label"] for dataset in datasets)
    win = pg.GraphicsLayoutWidget(show=True, title=title)
    win.resize(1200, 700)

    plot = win.addPlot()
    plot.showGrid(x=True, y=True)
    plot.setLabel("bottom", "x")
    plot.setLabel("left", "y")
    plot.addLegend(offset=(10, 10))

    # Mouse wheel/drag affects x-axis only.
    # Y-axis is automatically adjusted.
    plot.setMouseEnabled(x=True, y=False)

    for idx, dataset in enumerate(datasets):
        color = pg.intColor(idx, hues=max(3, len(datasets)))
        dataset["curve"] = plot.plot([], [], pen=pg.mkPen(color=color, width=1), name=dataset["label"])
        dataset["curve"].setZValue(idx)

    vb = plot.getViewBox()

    state = {
        "updating": False,
    }

    def autoscale_y(y_chunks):
        y_max = None
        for y_values in y_chunks:
            if len(y_values) == 0:
                continue
            candidate = np.nanmax(y_values)
            if not np.isfinite(candidate):
                continue
            if y_max is None or candidate > y_max:
                y_max = candidate

        if y_max is None:
            return

        if y_max <= 0:
            plot.setYRange(0.0, 1.0, padding=0.0)
        else:
            plot.setYRange(0.0, 1.2 * float(y_max), padding=0.0)

    def update_plot():
        if state["updating"]:
            return

        state["updating"] = True

        try:
            xmin, xmax = vb.viewRange()[0]

            if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
                return

            width_px = max(64, int(plot.width()))
            y_chunks = []

            for dataset in datasets:
                x = dataset["x"]
                y = dataset["y"]
                curve = dataset["curve"]
                i0, i1 = visible_indices(x, xmin, xmax)
                visible_count = i1 - i0

                if visible_count <= 0:
                    curve.setData([], [])
                    dataset["mode"] = None
                    continue

                points_per_pixel = visible_count / width_px
                overzoomed = points_per_pixel < OVERZOOM_POINTS_PER_PIXEL

                if overzoomed:
                    xp, yp = evaluate_local_cubic_spline(
                        x=x,
                        y=y,
                        xmin=xmin,
                        xmax=xmax,
                        width_px=width_px,
                        i0=i0,
                        i1=i1,
                    )

                    curve.setDownsampling(auto=False)
                    curve.setClipToView(False)
                    curve.setData(xp, yp)

                    # Use spline-rendered values for autoscaling so spline
                    # overshoot is not clipped.
                    y_chunks.append(yp)
                    dataset["mode"] = "local cubic spline"

                else:
                    curve.setData(x, y)
                    curve.setClipToView(True)
                    curve.setDownsampling(auto=True, method="peak")

                    y_chunks.append(y[i0:i1])
                    dataset["mode"] = "raw"

            autoscale_y(y_chunks)

        finally:
            state["updating"] = False

    vb.sigXRangeChanged.connect(lambda *_: update_plot())

    xmin = min(float(dataset["x"][0]) for dataset in datasets)
    xmax = max(float(dataset["x"][-1]) for dataset in datasets)
    plot.setXRange(xmin, xmax, padding=0.0)
    update_plot()

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
