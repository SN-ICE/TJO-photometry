#!/usr/bin/env python3
"""Plot a tjophot light curve against POISE and ZTF photometry of the same object.

    python3 plot_tjo_vs_poise.py <TARGET> [options]
    python3 plot_tjo_vs_poise.py TDE2025bri --mjd-max 61100
    python3 plot_tjo_vs_poise.py SN2026acd

Lives in the pipeline root so that `tjophot run --clean` cannot remove it.
Everything is resolved relative to this file:

    runs/<TARGET>/<TARGET>_lightcurve.csv    tjophot output   (required)
    runs/<TARGET>_poise.txt                  JD band mag err  (optional)
    runs/<TARGET>_detections.csv             ZTF alerts       (optional)
    runs/<TARGET>_non_detections.csv         ZTF limits       (optional)

The figure is written to runs/<TARGET>_tjo_vs_poise.{pdf,png} -- beside the
inputs rather than inside runs/<TARGET>/, so that it too survives --clean.
Use --outdir to put it somewhere else.

Bands are whatever the files contain, stacked in wavelength order with constant
offsets. Nothing about a particular object is hard-coded.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"

# Every band the pipeline or POISE might report, bluest first. The stack, the
# offsets and the labels all follow this order.
WAVELENGTH_ORDER = ["u", "U", "B", "g", "V", "r", "R", "i", "I", "z", "y", "C"]

# Bands that measure essentially the same part of the spectrum share a colour
# slot; they normally do not appear together (a target is reduced in SDSS *or*
# Johnson filters). When they do -- POISE i alongside TJO I, say -- the second
# one is pushed to the next free colour by assign_colours().
COLOR_SLOT = {"u": 0, "U": 0, "B": 1, "g": 2, "V": 3,
              "r": 4, "R": 4, "i": 5, "I": 5, "z": 6, "y": 6, "C": 6}

# Checked with the dataviz validator (adjacent pairs, light surface): lightness
# band, chroma floor, CVD separation and contrast all pass. The green/orange
# pair sits in the 6-8 CVD band, which is legal here because band identity is
# also carried by the vertical offset, the direct labels and the marker shape --
# never by colour alone.
RAMP = ["#a855f7", "#1d4ed8", "#0891b2", "#15803d",
        "#d97706", "#be123c", "#6b21a8"]
FALLBACK = "#475569"

ZTF_FID = {"1": "g", "2": "r", "3": "i"}
ZTF_MAX_ERR = 0.5      # sanity cut; the alert stream is already 5-sigma

INK = "#1a1a1a"
MUTED = "#6b7280"
GRID = "#e5e7eb"
SURFACE = "#fcfcfb"


# --------------------------------------------------------------------------
def pretty_target(target: str) -> str:
    """'SN2026acd' -> 'SN 2026acd', 'TDE2025bri' -> 'TDE 2025bri'."""
    for k, ch in enumerate(target):
        if ch.isdigit():
            return f"{target[:k]} {target[k:]}" if k else target
    return target


def assign_colours(bands: list[str]) -> dict[str, str]:
    """One colour per band, stable for a given band whenever the slot is free."""
    used, out = set(), {}
    for band in bands:
        slot = COLOR_SLOT.get(band, 6)
        while slot in used and slot < len(RAMP) - 1:
            slot += 1
        if slot in used:
            out[band] = FALLBACK
            continue
        used.add(slot)
        out[band] = RAMP[slot]
    return out


def mjd_to_date(mjd: float) -> str:
    return (datetime.date(1858, 11, 17) +
            datetime.timedelta(days=float(mjd))).strftime("%Y-%m-%d")


# --------------------------------------------------------------------------
def read_tjo(path: Path, lo: float, hi: float) -> dict:
    """Adopted measurement per epoch, split by how it was measured.

    'difference_*' means a template was subtracted; plain 'aperture'/'psf'
    still contains the host galaxy, which matters a great deal for anything
    sitting on its nucleus.
    """
    out: dict = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            band, mjd = row["band"], float(row["mjd"])
            if not (lo <= mjd <= hi):
                continue
            rec = out.setdefault(band, {"diff": [], "direct": [], "limit": []})
            if row["is_limit"].lower() == "true":
                if row["limit_mag"]:
                    rec["limit"].append((mjd, float(row["limit_mag"])))
                continue
            if not row["mag"]:
                continue
            key = "diff" if row["method"].startswith("difference") else "direct"
            rec[key].append((mjd, float(row["mag"]), float(row["mag_err"] or 0.0)))
    return out


def read_poise(path: Path, lo: float, hi: float) -> dict:
    out: dict = {}
    if not path.exists():
        return out
    with open(path) as fh:
        next(fh)
        for line in fh:
            if not line.strip():
                continue
            jd, band, mag, err = line.split()
            mjd = float(jd) - 2400000.5
            if lo <= mjd <= hi:
                out.setdefault(band, []).append((mjd, float(mag), float(err)))
    for band in out:
        out[band].sort()
    return out


def read_ztf(det_path: Path, lim_path: Path, lo: float, hi: float) -> dict:
    """ZTF alert photometry: detections, and the non-detections that constrain.

    'mag' is the PSF magnitude on the ZTF difference image -- the transient
    alone, which is what POISE and the template-subtracted TJO points also
    measure. ('mag_corr' adds the reference flux back in and so carries the
    host, like the un-subtracted TJO epochs.)

    Only non-detections *before* the first detection in their own band are
    kept. The rest come from shallow frames taken while the transient was
    comfortably detected on either side -- their limits are brighter than the
    source, so an arrow drawn there constrains nothing and reads as a
    contradiction.
    """
    out: dict = {}
    if det_path.exists():
        with open(det_path) as fh:
            for row in csv.DictReader(fh):
                band = ZTF_FID.get(row["fid"])
                mjd = float(row["mjd"])
                if band is None or not (lo <= mjd <= hi):
                    continue
                if row.get("dubious") == "true" or row.get("isdiffpos") != "1":
                    continue
                err = float(row["e_mag"])
                if err <= ZTF_MAX_ERR:
                    out.setdefault(band, {"det": [], "limit": []})["det"].append(
                        (mjd, float(row["mag"]), err))
    for band in out:
        out[band]["det"].sort()

    if lim_path.exists():
        with open(lim_path) as fh:
            for row in csv.DictReader(fh):
                band = ZTF_FID.get(row["fid"])
                mjd = float(row["mjd"])
                if band is None or not (lo <= mjd <= hi):
                    continue
                rec = out.setdefault(band, {"det": [], "limit": []})
                first = rec["det"][0][0] if rec["det"] else np.inf
                if mjd < first:
                    rec["limit"].append((mjd, float(row["diffmaglim"])))
    for band in out:
        out[band]["limit"].sort()
    return out


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="target name, e.g. SN2026acd")
    ap.add_argument("--bands", nargs="+", default=None,
                    help="restrict to these bands (default: everything present)")
    ap.add_argument("--mjd-min", type=float, default=-np.inf)
    ap.add_argument("--mjd-max", type=float, default=np.inf,
                    help="drop epochs past this MJD, e.g. an isolated late season")
    ap.add_argument("--offset", type=float, default=0.8,
                    help="magnitude step between band tracks (default 0.8)")
    ap.add_argument("--gap-days", type=float, default=15.0,
                    help="break the POISE guide line across gaps longer than this")
    ap.add_argument("--label-gap", type=float, default=0.40,
                    help="minimum separation between band labels, mag")
    ap.add_argument("--no-overlap", action="store_true",
                    help="do not shade the TJO/POISE overlap")
    ap.add_argument("--outdir", type=Path, default=None,
                    help="where to write the figure (default: the runs directory)")
    args = ap.parse_args(argv)

    target = args.target
    lc = RUNS / target / f"{target}_lightcurve.csv"
    if not lc.exists():
        print(f"no light curve at {lc}", file=sys.stderr)
        return 1

    lo, hi = args.mjd_min, args.mjd_max
    tjo = read_tjo(lc, lo, hi)
    poise = read_poise(RUNS / f"{target}_poise.txt", lo, hi)
    ztf = read_ztf(RUNS / f"{target}_detections.csv",
                   RUNS / f"{target}_non_detections.csv", lo, hi)

    present = set(tjo) | set(poise) | set(ztf)
    if args.bands:
        present &= set(args.bands)
    bands = [b for b in WAVELENGTH_ORDER if b in present]
    bands += sorted(present - set(bands))          # anything unrecognised, last
    if not bands:
        print("nothing to plot", file=sys.stderr)
        return 1

    colour_of = assign_colours(bands)
    mid = (len(bands) - 1) / 2.0
    offset_of = {b: (k - mid) * args.offset for k, b in enumerate(bands)}

    # ---- x range, with room on the left for the band labels ---------------
    xs = [x[0] for v in tjo.values() for k in v for x in v[k]]
    xs += [x[0] for v in poise.values() for x in v]
    xs += [x[0] for v in ztf.values() for k in v for x in v[k]]
    span = max(xs) - min(xs) or 1.0
    xlim = (min(xs) - 0.16 * span, max(xs) + 0.05 * span)

    fig, ax = plt.subplots(figsize=(10.4, 0.55 * len(bands) + 3.4))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    drawn, labels = set(), []
    for band in bands:
        colour, shift = colour_of[band], offset_of[band]
        t = tjo.get(band, {"diff": [], "direct": [], "limit": []})
        z = ztf.get(band, {"det": [], "limit": []})

        p = np.array(poise.get(band, []), dtype=float)
        if len(p):
            # Break the guide line across a season gap: a straight segment
            # drawn over empty months would read as a measured decline.
            ly, lx = (p[:, 1] + shift).copy(), p[:, 0].copy()
            ly[np.flatnonzero(np.diff(lx) > args.gap_days)] = np.nan
            ax.plot(lx, ly, "-", color=colour, lw=1.0, alpha=0.45, zorder=2)
            ax.errorbar(p[:, 0], p[:, 1] + shift, yerr=p[:, 2], fmt="o", ms=4.5,
                        color=colour, mec=SURFACE, mew=0.6, ecolor=colour,
                        elinewidth=0.9, capsize=0, zorder=3)
            drawn.add("poise")

        d = np.array(t["diff"], dtype=float)
        if len(d):
            ax.errorbar(d[:, 0], d[:, 1] + shift, yerr=d[:, 2], fmt="s", ms=6.5,
                        mfc="none", mec=colour, mew=1.5, ecolor=colour,
                        elinewidth=0.9, capsize=0, zorder=4)
            drawn.add("tjo_diff")

        n = np.array(t["direct"], dtype=float)
        if len(n):
            ax.errorbar(n[:, 0], n[:, 1] + shift, yerr=n[:, 2], fmt="D", ms=5.0,
                        mfc="none", mec=colour, mew=1.0, ecolor=colour,
                        elinewidth=0.7, alpha=0.55, capsize=0, zorder=4)
            drawn.add("tjo_direct")

        lim = np.array(t["limit"], dtype=float)
        if len(lim):
            ax.plot(lim[:, 0], lim[:, 1] + shift, "v", ms=6.0, mfc="none",
                    mec=colour, mew=1.0, alpha=0.55, zorder=4)
            drawn.add("tjo_limit")

        zd = np.array(z["det"], dtype=float)
        if len(zd):
            ax.errorbar(zd[:, 0], zd[:, 1] + shift, yerr=zd[:, 2], fmt="*", ms=9,
                        mfc="none", mec=colour, mew=1.0, ecolor=colour,
                        elinewidth=0.8, capsize=0, zorder=5)
            drawn.add("ztf_det")

        zl = np.array(z["limit"], dtype=float)
        if len(zl):
            ax.plot(zl[:, 0], zl[:, 1] + shift, "1", ms=12, color=colour,
                    mew=1.5, alpha=0.75, zorder=5)
            drawn.add("ztf_limit")

        # Label anchor: near the bright end of the track, not its first epoch.
        # Surveys start on different dates and a light curve is usually rising
        # steeply at whichever one comes first, so the left-hand edge is a poor
        # guide to where the track reads. The 5th percentile sits near peak
        # without being dragged by a single contaminated point. Drawn after the
        # loop, once the whole set can be de-collided.
        mags = np.concatenate([a[:, 1] for a in (p, d, n, zd) if len(a)]) \
            if any(len(a) for a in (p, d, n, zd)) else np.array([])
        if len(mags):
            labels.append([band, colour, shift,
                           float(np.percentile(mags, 5)) + shift])

    # Bands can sit within a magnitude of one another and the surveys start on
    # different dates, so the offsets alone do not guarantee the labels clear
    # each other or even come out in stack order. Walk them in wavelength order
    # and push each down by the least that works.
    for k in range(1, len(labels)):
        labels[k][3] = max(labels[k][3], labels[k - 1][3] + args.label_gap)
    for band, colour, shift, y in labels:
        sign = "+" if shift >= 0 else "−"
        ax.annotate(f"${band}$ {sign}{abs(shift):.1f}",
                    xy=(xlim[0] + 0.01 * span, y), color=colour, fontsize=11,
                    fontweight="semibold", va="center", ha="left", clip_on=False)

    # ---- shade the TJO/POISE overlap when it is only part of the range ----
    if poise and tjo and not args.no_overlap:
        p_all = [x[0] for v in poise.values() for x in v]
        t_all = [x[0] for v in tjo.values() for k in v for x in v[k]]
        olo, ohi = max(min(p_all), min(t_all)), min(max(p_all), max(t_all))
        if ohi > olo and (ohi - olo) < 0.75 * (xlim[1] - xlim[0]):
            ax.axvspan(olo, ohi, color="#0891b2", alpha=0.05, zorder=0)
            ax.annotate(f"overlap\n{olo:.0f}–{ohi:.0f}",
                        xy=((olo + ohi) / 2, 0.985),
                        xycoords=("data", "axes fraction"), ha="center",
                        va="top", fontsize=9, color=MUTED)

    ax.invert_yaxis()
    ax.set_xlim(*xlim)
    ax.set_xlabel("MJD", fontsize=11, color=INK)
    ax.set_ylabel("apparent magnitude + offset", fontsize=11, color=INK)
    ax.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=10)

    step = max(1.0, round((xlim[1] - xlim[0]) / 7 / 5) * 5)
    ticks = np.arange(np.ceil(min(xs) / step) * step, xlim[1], step)
    top = ax.secondary_xaxis("top")
    top.set_xticks(ticks)
    top.set_xticklabels([mjd_to_date(t) for t in ticks], fontsize=9, color=MUTED)
    top.spines["top"].set_color(GRID)
    top.tick_params(colors=MUTED)

    catalogue = [
        ("poise", dict(ls="-", marker="o", ms=5, color=MUTED, mec=SURFACE,
                       mew=0.6, lw=1.0, alpha=0.8, label="POISE")),
        ("tjo_diff", dict(ls="none", marker="s", ms=6.5, mfc="none", mec=INK,
                          mew=1.5, label="TJO, template-subtracted")),
        ("tjo_direct", dict(ls="none", marker="D", ms=5, mfc="none", mec=INK,
                            mew=1.0, alpha=0.6,
                            label="TJO, no subtraction (host included)")),
        ("tjo_limit", dict(ls="none", marker="v", ms=6, mfc="none", mec=INK,
                           mew=1.0, alpha=0.6, label="TJO 3$\\sigma$ upper limit")),
        ("ztf_det", dict(ls="none", marker="*", ms=9, mfc="none", mec=INK,
                         mew=1.0, label="ZTF, difference magnitude")),
        ("ztf_limit", dict(ls="none", marker="1", ms=12, color=INK, mew=1.5,
                           alpha=0.75,
                           label="ZTF pre-discovery limit (5$\\sigma$)")),
    ]
    handles = [Line2D([], [], **kw) for key, kw in catalogue if key in drawn]
    # Below the axes: pre-discovery limits tend to sit in the lower-left corner,
    # and a legend there would cover them.
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.085),
              ncol=min(3, len(handles)), fontsize=9.5, frameon=False,
              labelcolor=INK, handletextpad=0.6, columnspacing=2.2)

    surveys = [name for name, data in (("TJO/MEIA2", tjo), ("POISE", poise),
                                       ("ZTF", ztf)) if data]
    if len(surveys) > 1:
        surveys = [", ".join(surveys[:-1]) + " and " + surveys[-1]]
    ax.set_title(f"{pretty_target(target)} — {surveys[0]} photometry",
                 fontsize=14, color=INK, pad=26, loc="left")
    fig.text(0.5, 0.012,
             "Magnitudes as reported by each survey, in its own natural system. "
             "No cross-system colour terms applied.",
             fontsize=8.5, color=MUTED, ha="center")

    outdir = args.outdir or RUNS
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{target}_tjo_vs_poise"
    fig.tight_layout(rect=(0.025, 0.035, 1.0, 1))
    fig.savefig(f"{out}.pdf", facecolor=SURFACE)
    fig.savefig(f"{out}.png", dpi=200, facecolor=SURFACE)
    print(f"wrote {out}.pdf and {out}.png  ({len(bands)} bands: {' '.join(bands)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
