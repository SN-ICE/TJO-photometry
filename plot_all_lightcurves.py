#!/usr/bin/env python3
"""Every tjophot light curve in one figure: one panel per band, all objects.

    python3 plot_all_lightcurves.py
    python3 plot_all_lightcurves.py --bands g r i z --subtracted-only
    python3 plot_all_lightcurves.py --targets SN2026acd TDE2025bri

A sample-overview figure in the style of the Foundation DR1 light-curve plot:
apparent magnitude against MJD, one stacked panel per band, every object drawn
on top of every other, each object's own points joined by a thin dashed line so
that individual light curves stay readable inside the pile.

Lives in the pipeline root so that `tjophot run --clean` cannot remove it, and
so a rerun after reducing a new object picks that object up automatically.
Everything is resolved relative to this file:

    runs/<TARGET>/<TARGET>_lightcurve.csv     every one found is included

Writes runs/all_lightcurves.{pdf,png} (--outdir to put it elsewhere).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"

# Every band the pipeline might report, bluest first. Panels follow this order.
WAVELENGTH_ORDER = ["u", "U", "B", "g", "V", "r", "R", "i", "I", "z", "y", "C"]

# Bands covering essentially the same part of the spectrum share a colour slot;
# they rarely appear together (a target is reduced in SDSS *or* Johnson
# filters). Each band gets its own panel anyway, so a shared colour costs
# nothing here -- this only keeps the palette identical to plot_tjo_vs_poise.py.
COLOR_SLOT = {"u": 0, "U": 0, "B": 1, "g": 2, "V": 3,
              "r": 4, "R": 4, "i": 5, "I": 5, "z": 6, "y": 6, "C": 6}

# Checked with the dataviz validator (adjacent pairs, light surface): lightness
# band, chroma floor, CVD separation and contrast all pass. Band identity is
# carried by the panel and its label as well as by colour, never by colour
# alone.
RAMP = ["#a855f7", "#1d4ed8", "#0891b2", "#15803d",
        "#d97706", "#be123c", "#6b21a8"]

MARKER = {"u": "o", "U": "o", "B": "s", "g": "o", "V": "^", "r": "s",
          "R": "v", "i": "^", "I": "D", "z": "D", "y": "v", "C": "X"}

INK = "#1a1a1a"
MUTED = "#6b7280"
JOIN = "#9aa0a6"        # the thin dashed line tying one object's points
SURFACE = "#fcfcfb"
SHADE = "#6b7280"       # alternating calendar-year bands

MJD_EPOCH = date(1858, 11, 17)      # MJD 0


def mjd_of_year(year: int) -> float:
    """MJD at 00:00 UTC on 1 January of `year`."""
    return float((date(year, 1, 1) - MJD_EPOCH).days)


def year_spans(mjd_lo: float, mjd_hi: float):
    """[(year, start_mjd, end_mjd), ...] covering the range, ends included."""
    first = (MJD_EPOCH + timedelta(days=int(np.floor(mjd_lo)))).year
    last = (MJD_EPOCH + timedelta(days=int(np.ceil(mjd_hi)))).year
    return [(y, mjd_of_year(y), mjd_of_year(y + 1))
            for y in range(first, last + 1)]


def read_all(targets, subtracted_only: bool, flag_exclude: set[str]):
    """{band: {target: [(mjd, mag, err), ...]}} from every run directory."""
    out: dict = defaultdict(lambda: defaultdict(list))
    for csv_path in sorted(RUNS.glob("*/*_lightcurve.csv")):
        target = csv_path.name[: -len("_lightcurve.csv")]
        if targets and target not in targets:
            continue
        with open(csv_path) as fh:
            for row in csv.DictReader(fh):
                if row["is_limit"].lower() == "true" or not row["mag"]:
                    continue
                if subtracted_only and not row["method"].startswith("difference"):
                    continue
                if flag_exclude & set(filter(None, row["flags_str"].split("|"))):
                    continue
                out[row["band"]][target].append(
                    (float(row["mjd"]), float(row["mag"]),
                     float(row["mag_err"] or 0.0)))
    for band in out:
        for target in out[band]:
            out[band][target].sort()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", nargs="+", default=None,
                    help="restrict to these objects (default: every run found)")
    ap.add_argument("--bands", nargs="+", default=None,
                    help="restrict to these bands (default: everything present)")
    ap.add_argument("--min-epochs", type=int, default=1,
                    help="skip a band with fewer detections than this")
    ap.add_argument("--subtracted-only", action="store_true",
                    help="keep only template-subtracted epochs")
    ap.add_argument("--exclude-flags", nargs="+",
                    default=["bad_subtraction", "no_zeropoint"],
                    help="drop epochs carrying any of these flags "
                         "(pass none to keep everything)")
    ap.add_argument("--mjd-offset", type=float, default=None,
                    help="x axis is MJD minus this (default: a round number "
                         "just below the earliest epoch)")
    ap.add_argument("--gap-days", type=float, default=60.0,
                    help="break an object's line across gaps longer than this; "
                         "0 joins everything")
    ap.add_argument("--mag-range", nargs=2, type=float, default=None,
                    metavar=("BRIGHT", "FAINT"))
    ap.add_argument("--errorbars", action="store_true",
                    help="draw magnitude error bars (off: too dense by default)")
    ap.add_argument("--no-year-shading", action="store_true",
                    help="do not shade the odd calendar years")
    ap.add_argument("--year-alpha", type=float, default=0.1,
                    help="opacity of the odd-year shading (default 0.1)")
    ap.add_argument("--year-labels", choices=("top", "inside", "none"),
                    default="top",
                    help="where to write the year numbers: 'top' above the "
                         "first panel (default), 'inside' along the top of it "
                         "(may collide with the title), or 'none'")
    ap.add_argument("--panel-height", type=float, default=1.75)
    ap.add_argument("--title", default="TJO/MEIA2",
                    help="label for the top-left of the first panel")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--name", default="all_lightcurves",
                    help="output file stem")
    args = ap.parse_args(argv)

    flag_exclude = {f for f in args.exclude_flags if f.lower() != "none"}
    data = read_all(set(args.targets) if args.targets else None,
                    args.subtracted_only, flag_exclude)
    if not data:
        print("no light curves found under runs/", file=sys.stderr)
        return 1

    present = set(data)
    if args.bands:
        present &= set(args.bands)
    bands = [b for b in WAVELENGTH_ORDER if b in present]
    bands += sorted(present - set(bands))
    bands = [b for b in bands
             if sum(len(v) for v in data[b].values()) >= args.min_epochs]
    if not bands:
        print("nothing left to plot", file=sys.stderr)
        return 1

    points = [p for b in bands for v in data[b].values() for p in v]
    mjds = np.array([p[0] for p in points])
    mags = np.array([p[1] for p in points])

    offset = (args.mjd_offset if args.mjd_offset is not None
              else float(np.floor(mjds.min() / 1000.0) * 1000.0))
    pad = 0.03 * (mjds.max() - mjds.min() or 1.0)
    xlim = (mjds.min() - offset - pad, mjds.max() - offset + pad)
    if args.mag_range:
        ylim = tuple(args.mag_range)
    else:
        ylim = (np.floor(mags.min() - 0.6), np.ceil(mags.max() + 0.6))

    # One colour per band, same slot rule as plot_tjo_vs_poise.py. Bands sharing
    # a slot (i and I, say) are nudged to the nearest free colour -- they sit in
    # separate panels, so this is only tidiness.
    used, colour_of = set(), {}
    for band in bands:
        want = COLOR_SLOT.get(band, len(RAMP) - 1)
        slot = next((s for s in range(want, len(RAMP)) if s not in used),
                    next((s for s in range(want, -1, -1) if s not in used), want))
        used.add(slot)
        colour_of[band] = RAMP[slot]

    fig, axes = plt.subplots(
        len(bands), 1, sharex=True,
        figsize=(9.2, args.panel_height * len(bands) + 1.4))
    axes = np.atleast_1d(axes)
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(hspace=0.0)

    spans = year_spans(xlim[0] + offset, xlim[1] + offset)

    tally: list = []
    for ax, band in zip(axes, bands):
        ax.set_facecolor(SURFACE)
        colour, marker = colour_of[band], MARKER.get(band, "o")
        n_obj = n_pt = 0

        for target, rows in sorted(data[band].items()):
            arr = np.array(rows, dtype=float)
            x, y = arr[:, 0] - offset, arr[:, 1]
            n_obj += 1
            n_pt += len(arr)

            # One object, one line. Broken across a long gap so that a second
            # season is not tied to the first by a straight segment that no
            # measurement supports.
            if len(arr) > 1:
                yl = y.copy()
                if args.gap_days > 0:
                    yl[np.flatnonzero(np.diff(x) > args.gap_days)] = np.nan
                ax.plot(x, yl, ls="--", lw=1.4, color=JOIN, alpha=0.85, zorder=2)

            if args.errorbars:
                ax.errorbar(x, y, yerr=arr[:, 2], fmt="none", ecolor=colour,
                            elinewidth=0.6, alpha=0.6, capsize=0, zorder=3)
            ax.plot(x, y, ls="none", marker=marker, ms=5.0, mfc=colour,
                    mec="black", mew=0.5, zorder=4)

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        # Alternating calendar years, so a bump can be placed on the calendar
        # without arithmetic. Odd years carry the tint; even ones are the bare
        # surface. zorder 0 keeps it behind the joining lines and the markers,
        # and axvspan spans the full height whatever the y limits end up being.
        if not args.no_year_shading:
            for year, lo, hi in spans:
                if year % 2:
                    ax.axvspan(lo - offset, hi - offset, color=SHADE,
                               alpha=args.year_alpha, lw=0, zorder=0)

        ax.invert_yaxis()
        ax.tick_params(direction="in", top=True, right=True, which="both",
                       colors=INK, labelsize=10)
        for spine in ax.spines.values():
            spine.set_color(INK)
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(AutoMinorLocator())

        # Marker plus dashed line, exactly as the panel draws them.
        ax.legend(handles=[Line2D([], [], marker=marker, ms=7, mfc=colour,
                                  mec="black", mew=0.5, ls="--", lw=1.4,
                                  color=JOIN, label=f"${band}$")],
                  loc="upper left", frameon=False, handletextpad=0.4,
                  borderaxespad=0.7,
                  bbox_to_anchor=(0.0, 0.80) if ax is axes[0] else None,
                  labelcolor=colour, fontsize=13)
        tally.append((band, n_obj, n_pt))

    # One number per year, centred on the part of that year that is actually
    # visible, so the first and last years are labelled where they are seen
    # rather than off the edge. A sliver too narrow for the text is skipped.
    if args.year_labels != "none":
        top = args.year_labels == "top"
        width = xlim[1] - xlim[0]
        for year, lo, hi in spans:
            left, right = max(lo - offset, xlim[0]), min(hi - offset, xlim[1])
            if right - left < 0.045 * width:
                continue
            axes[0].annotate(
                str(year), xy=(0.5 * (left + right), 1.012 if top else 0.975),
                xycoords=("data", "axes fraction"),
                ha="center", va="bottom" if top else "top",
                fontsize=10, color=MUTED, clip_on=not top)

    axes[0].annotate(args.title, xy=(0.012, 0.93), xycoords="axes fraction",
                     ha="left", va="top", fontsize=15, color=INK)
    axes[-1].set_xlabel(f"MJD − {offset:.0f}", fontsize=13, color=INK)
    fig.supylabel("Apparent Brightness (mag)", fontsize=13, color=INK, x=0.035)

    outdir = args.outdir or RUNS
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / args.name
    fig.savefig(f"{out}.pdf", facecolor=SURFACE, bbox_inches="tight")
    fig.savefig(f"{out}.png", dpi=200, facecolor=SURFACE, bbox_inches="tight")
    for band, n_obj, n_pt in tally:
        print(f"  {band:>2}  {n_obj:3d} object{'s' if n_obj != 1 else ' '}"
              f"  {n_pt:5d} epochs")
    print(f"wrote {out}.pdf and {out}.png  "
          f"({len(bands)} bands, "
          f"{len({t for b in bands for t in data[b]})} objects, "
          f"{len(points)} epochs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
