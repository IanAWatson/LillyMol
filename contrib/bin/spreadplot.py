#!/usr/bin/env python3
"""Plot sel,distance CSV trajectories written by nplotnn -S.

No molecule or fingerprint parsing is performed. Matplotlib is loaded only when
rendering; SciPy is required only for optional monotonic spline interpolation.
"""

from __future__ import annotations

import argparse
import csv
import math
from contextlib import nullcontext
from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, localcontext
from pathlib import Path
import sys
from typing import Sequence


@dataclass(frozen=True)
class Trajectory:
    source: str
    label: str
    selected: list[int]
    distance: list[float]

    def fraction(self, index: int) -> float:
        # nplotnn positions are zero-based, including in reduced trajectories.
        return (self.selected[index] + 1) / (self.selected[-1] + 1)


def parse_thresholds(value: str) -> list[float]:
    """Parse one comma-separated --dist argument."""
    try:
        fields = value.split(",")
        if any(not field.strip() for field in fields):
            raise ValueError
        result = [float(field) for field in fields]
        if any(not math.isfinite(x) or not 0 <= x <= 1 for x in result):
            raise ValueError
        return result
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--dist requires comma-separated finite distances between 0 and 1") from exc


def crossing_point(values: Sequence[float], cutoff: float) -> int | None:
    """First array index at or below cutoff in a non-increasing trajectory."""
    first, last = 0, len(values)
    while first < last:
        middle = (first + last) // 2
        if values[middle] <= cutoff:
            last = middle
        else:
            first = middle + 1
    return first if first < len(values) else None


def crossing_points(values: Sequence[float], delta: float) -> list[int]:
    """Markers at 1-delta, 1-2*delta, ...; each observation is marked once."""
    if not math.isfinite(delta) or not 0 < delta <= 1:
        raise ValueError("--point must be greater than zero and at most 1")
    interval = Decimal(str(delta))
    result = []
    with localcontext() as context:
        # Keep the regular grid exact even when the requested interval is too
        # small for subtraction from the float 1.0 to change its value.
        context.prec = max(50, -interval.adjusted() + 25)
        end = int((Decimal(1) / interval).to_integral_value(rounding=ROUND_CEILING))
        multiple = 1
        while multiple < end:
            threshold = Decimal(1) - multiple * interval
            first, last = 0, len(values)
            while first < last:
                middle = (first + last) // 2
                if Decimal(str(values[middle])) <= threshold:
                    last = middle
                else:
                    first = middle + 1
            if first == len(values):
                break
            if not result or result[-1] != first:
                result.append(first)
            # Jump over every additional threshold crossed at this observation.
            below = int((Decimal(1) - Decimal(str(values[first]))) // interval) + 1
            multiple = max(multiple + 1, below)
    return result


def read_trajectory(filename: str, label: str = "Selected") -> Trajectory:
    """Read either the full or plateau-reduced nplotnn CSV, preserving sel."""
    selected, distance = [], []
    context = nullcontext(sys.stdin) if filename == "-" else open(
        filename, newline="", encoding="utf-8-sig")
    try:
        with context as input_file:
            reader = csv.DictReader(input_file, strict=True)
            fields = reader.fieldnames or []
            if "sel" not in fields or "distance" not in fields:
                raise ValueError(f"{filename}: input must contain sel,distance columns")
            if len(set(fields)) != len(fields):
                raise ValueError(f"{filename}: duplicate CSV column names")
            for row in reader:
                try:
                    if None in row:
                        raise ValueError("too many CSV fields")
                    position = (row["sel"] or "").strip()
                    if not position.isascii() or not position.isdecimal():
                        raise ValueError("sel must be a non-negative integer")
                    position = int(position)
                    value = float(row["distance"])
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        raise ValueError("distance must be finite and between 0 and 1")
                    if selected and position <= selected[-1]:
                        raise ValueError("sel must be strictly increasing")
                    if distance and value > distance[-1]:
                        raise ValueError("distance must be monotonically non-increasing")
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"{filename}: CSV line {reader.line_num}: {exc}") from exc
                selected.append(position)
                distance.append(value)
    except csv.Error as exc:
        raise ValueError(f"{filename}: malformed CSV: {exc}") from exc
    if len(selected) < 2:
        raise ValueError(f"{filename}: at least two trajectory rows are required")
    return Trajectory(filename, label, selected, distance)


def threshold_rows(trajectories: Sequence[Trajectory],
                   thresholds: Sequence[float]) -> list[dict]:
    """Build the same threshold results used in plot legends and CSV output."""
    rows = []
    for trajectory in trajectories:
        for threshold in thresholds:
            index = crossing_point(trajectory.distance, threshold)
            rows.append(dict(
                file=trajectory.source, label=trajectory.label, threshold=threshold,
                status=("not_reached" if index is None else
                        "at_start" if index == 0 else "crossed"),
                sel="" if index is None else trajectory.selected[index],
                fraction_selected="" if index is None else trajectory.fraction(index)))
    return rows


def plot_points(trajectory: Trajectory, plateaus: bool = False):
    """Use first positions at distinct distance levels for visual interpolation.

    Keep the last position too, so the displayed extent is unchanged. The input
    observations remain intact for threshold calculations and marker placement.
    """
    if plateaus:
        return trajectory.selected, trajectory.distance
    indices = [0]
    indices.extend(i for i in range(1, len(trajectory.distance))
                   if trajectory.distance[i] != trajectory.distance[i - 1])
    if indices[-1] != len(trajectory.distance) - 1:
        indices.append(len(trajectory.distance) - 1)
    return ([trajectory.selected[i] for i in indices],
            [trajectory.distance[i] for i in indices])


def make_plot(trajectories: Sequence[Trajectory], thresholds: Sequence[float], *,
              colors: Sequence[str] = (), title: str | None = None,
              lwd: float = 4, point: float | None = None, spline: int | None = None,
              normalize: bool = False, logx: bool = False, plateaus: bool = False,
              width: float = 10, height: float = 6, dpi: int = 150):
    """Create a headless figure. Crossings always use original observations."""
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        from matplotlib.ticker import FuncFormatter
    except ImportError as exc:
        raise RuntimeError(
            "Matplotlib is required: uv pip install matplotlib "
            "(or python3 -m pip install matplotlib)") from exc

    figure = Figure(figsize=(width, height), dpi=dpi, layout="constrained")
    FigureCanvasAgg(figure)
    axes = figure.subplots()
    for number, trajectory in enumerate(trajectories):
        x, y = plot_points(trajectory, plateaus)
        if spline is not None:
            try:
                import numpy as np
                from scipy.interpolate import PchipInterpolator
            except ImportError as exc:
                raise RuntimeError(
                    "--spline requires SciPy: uv pip install scipy "
                    "(or python3 -m pip install scipy)") from exc
            interpolation = PchipInterpolator(x, y)
            x = np.linspace(x[0], x[-1], spline)
            y = interpolation(x)
        x = [(value + 1) / (trajectory.selected[-1] + 1) for value in x] if normalize else x
        # A zero selection position cannot be displayed on a logarithmic axis.
        visible = [(a, b) for a, b in zip(x, y) if not logx or a > 0]
        color = colors[number % len(colors)] if colors else (
            "green" if len(trajectories) == 1 else None)
        line, = axes.plot([a for a, _ in visible], [b for _, b in visible],
                          label=trajectory.label, color=color, linewidth=lwd)
        if point is not None:
            indices = crossing_points(trajectory.distance, point)
            marker_x = [trajectory.fraction(i) if normalize else trajectory.selected[i]
                        for i in indices]
            markers = [(a, trajectory.distance[i]) for a, i in zip(marker_x, indices)
                       if not logx or a > 0]
            axes.scatter([a for a, _ in markers], [b for _, b in markers],
                         color=line.get_color(), zorder=3, label="_nolegend_")

    rows = threshold_rows(trajectories, thresholds)
    for threshold in thresholds:
        descriptions = []
        for row in rows:
            if row["threshold"] != threshold:
                continue
            if row["status"] == "not_reached":
                location = "not reached"
            elif row["status"] == "at_start":
                location = "at start"
            else:
                location = f'selected {row["sel"]:,}'
                if normalize:
                    location += f' ({row["fraction_selected"]:.1%})'
            prefix = f'{row["label"]}: ' if len(trajectories) > 1 else ""
            descriptions.append(prefix + location)
        axes.axhline(threshold, color="grey", linestyle="--", linewidth=1,
                     label=f"Distance {threshold:g}: " + "; ".join(descriptions))

    # The initial spread distance is usually an uninformative placeholder.
    # Deliberately base the y-limit on the SECOND observation in each curve.
    ymax = max([t.distance[1] for t in trajectories] + list(thresholds))
    axes.set_ylim(0, ymax if ymax > 0 else 0.01)
    axes.set_ylabel("Distance")
    axes.set_xlabel("Fraction selected" if normalize else "Selection position (zero-based)")
    if logx:
        axes.set_xscale("log")
    else:
        axes.set_xlim(left=0)
    if normalize:
        axes.set_xlim(right=1)
        axes.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{100 * x:g}%"))
    else:
        axes.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    if title:
        axes.set_title(title)
    axes.grid(alpha=0.2)
    axes.legend()
    return figure


def positive_float(value: str) -> float:
    try:
        result = float(value)
        if not math.isfinite(result) or result <= 0:
            raise ValueError
        return result
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero") from exc


def positive_int(value: str) -> int:
    try:
        result = int(value)
        if result < 1:
            raise ValueError
        return result
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Plot full or reduced sel,distance CSV files from nplotnn -S.",
        epilog="Example: spreadplot.py -S spread --dist 0.1,0.2 "
               "--summary thresholds.csv trajectory.csv")
    result.add_argument("input", nargs="+", help="nplotnn CSV files; - reads stdin")
    result.add_argument("-S", required=True, metavar="STEM", help="output filename stem")
    result.add_argument("--format", choices=("png", "svg", "pdf"), default="png",
                        type=str.lower, help="image format (default: png)")
    result.add_argument("--dist", action="append", type=parse_thresholds, default=[],
                        metavar="D[,D...]", help="horizontal distance thresholds; repeatable")
    result.add_argument("--summary", metavar="CSV",
                        help="threshold summary CSV; requires --dist; - writes stdout")
    result.add_argument("--label", action="append", default=[], metavar="TEXT",
                        help="curve label; repeat once per input when plotting several files")
    result.add_argument("--color", action="append", default=[], metavar="COLOR",
                        help="curve colour; supply one for all curves or one per input")
    result.add_argument("--title")
    result.add_argument("--lwd", type=positive_float, default=4, help="line width (default: 4)")
    result.add_argument("--point", type=positive_float, metavar="DELTA",
                        help="markers at distances 1-DELTA, 1-2*DELTA, ...; at most 1")
    result.add_argument("--plateaus", action="store_true",
                        help="retain the exact input polyline instead of joining distance levels")
    result.add_argument("--spline", type=positive_int, metavar="N",
                        help="plot N monotonic interpolation samples (at least 2; needs SciPy); "
                             "use at least as many samples as input datapoints")
    result.add_argument("--normalize", action="store_true",
                        help="plot fractions (sel+1)/(last_sel+1) to compare collection sizes")
    result.add_argument("--logx", action="store_true",
                        help="logarithmic x-axis; omit sel=0 from the unnormalized plot")
    result.add_argument("--width", type=positive_float, default=10, help="width in inches")
    result.add_argument("--height", type=positive_float, default=6, help="height in inches")
    result.add_argument("--dpi", type=positive_int, default=150, help="PNG resolution (default: 150)")
    return result


def check_outputs(inputs: Sequence[str], outputs: Sequence[str]) -> None:
    """Do not overwrite inputs or alias the image and summary output."""
    resolved_inputs = [Path(name).resolve() for name in inputs if name != "-"]
    resolved_outputs = [Path(name).resolve() for name in outputs if name != "-"]
    previous = []
    for path in resolved_outputs:
        for other in resolved_inputs + previous:
            if path == other or (path.exists() and other.exists() and path.samefile(other)):
                raise ValueError(f"output '{path}' would overwrite an input or another output")
        previous.append(path)


def main(argv: Sequence[str] | None = None) -> int:
    argument_parser = parser()
    args = argument_parser.parse_args(argv)
    thresholds = list(dict.fromkeys(value for group in args.dist for value in group))
    if args.summary and not thresholds:
        argument_parser.error("--summary requires --dist")
    if args.input.count("-") > 1:
        argument_parser.error("standard input (-) can only be read once")
    if args.label and len(args.label) != len(args.input):
        argument_parser.error("supply one --label per input")
    if args.color and len(args.color) not in (1, len(args.input)):
        argument_parser.error("supply one --color or one per input")
    if args.point is not None and args.point > 1:
        argument_parser.error("--point must be at most 1")
    if args.spline is not None and args.spline < 2:
        argument_parser.error("--spline must be at least 2")
    image = args.S + "." + args.format
    try:
        check_outputs(args.input, [image] + ([args.summary] if args.summary else []))
        trajectories = []
        for number, filename in enumerate(args.input):
            label = (args.label[number] if args.label else
                     "Selected" if len(args.input) == 1 else
                     "stdin" if filename == "-" else Path(filename).stem)
            trajectories.append(read_trajectory(filename, label))
        if args.spline is not None:
            for trajectory in trajectories:
                count = len(trajectory.selected)
                if args.spline < count:
                    print(f"Warning: {trajectory.source}: --spline {args.spline} is "
                          f"smaller than the {count} input datapoints; use at least "
                          f"--spline {count} to avoid losing detail.", file=sys.stderr)
        figure = make_plot(
            trajectories, thresholds, colors=args.color, title=args.title,
            lwd=args.lwd, point=args.point, spline=args.spline,
            normalize=args.normalize, logx=args.logx, plateaus=args.plateaus,
            width=args.width, height=args.height, dpi=args.dpi)
        figure.savefig(image, format=args.format, dpi=args.dpi)
        figure.clear()
        if args.summary:
            context = nullcontext(sys.stdout) if args.summary == "-" else open(
                args.summary, "w", newline="", encoding="utf-8")
            with context as output:
                writer = csv.DictWriter(output, fieldnames=[
                    "file", "label", "threshold", "status", "sel", "fraction_selected"])
                writer.writeheader()
                writer.writerows(threshold_rows(trajectories, thresholds))
    except (ValueError, OSError, RuntimeError, ImportError, OverflowError) as exc:
        argument_parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
