"""Judging a difference image where it matters: at the sources.

The reduced chi-square the kernel solve returns is measured on the stamps used
to fit the kernel, which are bright field stars. It says the PSF match worked
*there*. It does not say the subtraction is clean at the transient, and a frame
can subtract well across the field and badly on a bright host -- which is how
three SN2022xkq epochs came through with errors of 0.7 to 1.1 mag while every
chi-square on that target stayed below 1.8.

The test that catches those is the shape of what is left behind. A residual
from a photometric mismatch is a monopole: the star over- or under-subtracts by
some fraction and leaves a blob of one sign. A residual from a *registration*
error is a dipole -- positive on one side, negative on the other -- because the
same flux has been subtracted a fraction of a pixel away from where it sits.
Dipoles are the signature no flux rescaling can remove, and counting them over
the field stars is a direct measure of whether the alignment is good enough.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .utils import get_logger

log = get_logger("tjophot.diffquality")


@dataclass
class DifferenceQuality:
    verdict: str = "unknown"          # 'pass' | 'downgrade' | 'fail' | 'unknown'
    dipole_fraction: float = float("nan")
    n_sources: int = 0
    residual_sigma: float = float("nan")   # bright-star residual, in noise units
    target_dipole: bool = False
    message: str = ""

    def as_dict(self) -> dict:
        return {"diff_verdict": self.verdict,
                "diff_dipole_fraction": self.dipole_fraction,
                "diff_quality_nsources": self.n_sources,
                "diff_residual_sigma": self.residual_sigma,
                "diff_target_dipole": bool(self.target_dipole)}


def _noise(diff: np.ndarray) -> float:
    finite = diff[np.isfinite(diff)]
    if finite.size < 16:
        return float("nan")
    med = float(np.median(finite))
    return float(np.median(np.abs(finite - med))) * 1.4826


def _is_dipole(stamp: np.ndarray, sigma: float, nsigma: float,
               antisymmetry: float) -> Optional[bool]:
    """Positive and negative peaks, both significant, on opposite sides.

    Returns None when the stamp cannot be judged (all masked, no noise
    estimate), which is not the same as judging it clean.
    """
    if not np.isfinite(sigma) or sigma <= 0:
        return None
    good = np.isfinite(stamp)
    if good.sum() < 9:
        return None
    values = np.where(good, stamp, 0.0)
    hi = np.unravel_index(int(np.argmax(values)), values.shape)
    lo = np.unravel_index(int(np.argmin(values)), values.shape)
    pos, neg = float(values[hi]), float(values[lo])
    if pos < nsigma * sigma or neg > -nsigma * sigma:
        return False
    # Opposite sides of the stamp centre, and comparable in size: a bright
    # source beside a faint negative is not a dipole, it is a detection.
    cy, cx = (np.array(stamp.shape) - 1) / 2.0
    v1 = np.array([hi[0] - cy, hi[1] - cx])
    v2 = np.array([lo[0] - cy, lo[1] - cx])
    if float(v1 @ v2) >= 0:
        return False
    ratio = min(abs(pos), abs(neg)) / max(abs(pos), abs(neg), 1e-12)
    return bool(ratio >= antisymmetry)


def assess(diff, positions, fwhm_px: float, cfg,
           target_xy: Optional[tuple] = None) -> DifferenceQuality:
    """Score a difference image on the shape of its residuals.

    *positions* are (x, y) of field sources on the difference stamp.
    """
    q = DifferenceQuality()
    if not bool(cfg.get("subtraction.quality.enabled", True)):
        q.verdict = "pass"
        return q

    data = np.asarray(diff.difference, dtype=float)
    sigma = _noise(data)
    if not np.isfinite(sigma) or sigma <= 0:
        q.message = "no noise estimate"
        return q

    nsigma = float(cfg.get("subtraction.quality.dipole_nsigma", 5.0))
    anti = float(cfg.get("subtraction.quality.dipole_antisymmetry", 0.3))
    half = max(int(round(1.5 * float(fwhm_px))), 3)
    ny, nx = data.shape

    def stamp(x, y):
        x0, x1 = int(round(x)) - half, int(round(x)) + half + 1
        y0, y1 = int(round(y)) - half, int(round(y)) + half + 1
        if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny:
            return None
        return data[y0:y1, x0:x1]

    verdicts, residuals = [], []
    for x, y in np.atleast_2d(np.asarray(positions, dtype=float)):
        sub = stamp(x, y)
        if sub is None:
            continue
        got = _is_dipole(sub, sigma, nsigma, anti)
        if got is None:
            continue
        verdicts.append(got)
        residuals.append(float(np.nanmedian(np.abs(sub))) / sigma)

    q.n_sources = len(verdicts)
    if residuals:
        q.residual_sigma = float(np.median(residuals))
    if target_xy is not None:
        sub = stamp(*target_xy)
        if sub is not None:
            q.target_dipole = bool(_is_dipole(sub, sigma, nsigma, anti))

    min_sources = int(cfg.get("subtraction.quality.min_sources", 8))
    if q.n_sources < min_sources:
        q.message = f"only {q.n_sources} sources to judge on"
        return q

    q.dipole_fraction = float(np.mean(verdicts))
    warn = float(cfg.get("subtraction.quality.dipole_warn", 0.05))
    fail = float(cfg.get("subtraction.quality.dipole_fail", 0.20))
    if q.dipole_fraction >= fail:
        q.verdict, q.message = "fail", (
            f"{100 * q.dipole_fraction:.0f}% of sources leave a dipole -- the "
            f"reference is misregistered, not merely mismatched in flux")
    elif q.dipole_fraction >= warn or q.target_dipole:
        q.verdict = "downgrade"
        q.message = ("the transient itself leaves a dipole" if q.target_dipole
                     else f"{100 * q.dipole_fraction:.0f}% of sources leave a dipole")
    else:
        q.verdict = "pass"
    return q
