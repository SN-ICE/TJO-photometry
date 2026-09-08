"""Detector configurations: which camera, read out how.

A single telescope's archive is not a single instrument. TJO's frames span at
least three configurations:

* MEIA2 read out as a centred 1048 x 1048 window   (6.3 arcmin, 0.361 "/px)
* MEIA2 read out full frame, 2048 x 2048          (12.3 arcmin, 0.361 "/px)
* MEIA3, 3108 x 3108                              (20.8 arcmin, 0.402 "/px,
                                                   gain 1.94, read noise 8.3)

Two of those differ only in how much of the same chip was digitised; the third
is a different camera with different gain, read noise and plate scale. Both
kinds of change matter, and they matter in the same place, so the pipeline
treats them uniformly: a *detector configuration* is the tuple (instrument,
stored array shape, binning), and anything measured on the detector -- fringe
patterns, illumination residuals -- is keyed by it.

That is not fussiness. A fringe pattern is fixed in detector coordinates and
has no meaning on an array of a different size, and a windowed readout is not
a rescaling of the full frame: it is a crop, at an origin the header does not
record (TJO writes OFFSET1 = OFFSET2 = 0 for both). Stretching a 2048-pixel
pattern across a 1048-pixel frame would put the ripples in the wrong place and
add structure rather than remove it.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Optional

from .utils import slugify


@dataclass(frozen=True)
class Detector:
    """How one frame was read out. Hashable, so it can key a dict."""

    instrument: str
    ny: int
    nx: int
    binning: int = 1

    @property
    def token(self) -> str:
        """Filesystem-safe key, e.g. ``meia2_2048x2048``, ``meia3_3108x3108_b2``."""
        name = slugify(self.instrument).lower() or "unknown"
        token = f"{name}_{self.nx}x{self.ny}"
        if self.binning != 1:
            token += f"_b{self.binning}"
        return token

    @property
    def label(self) -> str:
        """Short human-readable form for logs and tables."""
        text = f"{self.instrument} {self.nx}x{self.ny}"
        if self.binning != 1:
            text += f" bin{self.binning}"
        return text

    def field_arcmin(self, pixel_scale_arcsec: float) -> float:
        """Width of the field in arcmin, given the frame's pixel scale."""
        return max(self.nx, self.ny) * float(pixel_scale_arcsec) / 60.0

    def matches_shape(self, shape) -> bool:
        return tuple(int(v) for v in shape[:2]) == (self.ny, self.nx)


def _int(value, default: int = 0) -> int:
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return default


def identify(header, instrument_map: Optional[dict] = None) -> Detector:
    """The detector configuration a FITS header describes.

    ``instrument_map`` optionally renames the raw INSTRUME value, for archives
    where the same camera was written under more than one spelling.
    """
    raw = str(header.get("INSTRUME", "") or "unknown").strip()
    if instrument_map:
        raw = str(instrument_map.get(raw, raw)).strip()
    ny = _int(header.get("NAXIS2"), 0)
    nx = _int(header.get("NAXIS1"), 0)
    binning = _int(header.get("XFACTOR"), 1) or 1
    return Detector(instrument=raw or "unknown", ny=ny, nx=nx, binning=binning)


def group_by_detector(entries: Iterable, header_of) -> dict:
    """Split an iterable of anything into ``{Detector: [entry, ...]}``.

    ``header_of`` maps an entry to its FITS header.
    """
    groups: dict = {}
    for entry in entries:
        header = header_of(entry)
        if header is None:
            continue
        groups.setdefault(identify(header), []).append(entry)
    return groups


# --------------------------------------------------------------------------
# sizes that should follow the frame
# --------------------------------------------------------------------------
def field_radius_arcmin(nx: int, ny: int, pixel_scale_arcsec: float) -> float:
    """Radius of the smallest circle containing the whole frame, in arcmin."""
    return math.hypot(float(nx), float(ny)) / 2.0 * float(pixel_scale_arcsec) / 60.0


def resolve_size(value, frame_arcmin: float, fraction: float = 1.0,
                 lo: Optional[float] = None, hi: Optional[float] = None,
                 default: float = 9.0, step: float = 1.0,
                 rounding: str = "up") -> float:
    """Read a config field that may be a number or the word ``auto``.

    A fixed arcmin figure is right only for one field of view. TJO's three
    readouts span 6.3 to 20.8 arcmin, so a catalogue radius tuned for the
    middle one fetches four times too much sky for the smallest and misses
    half the field of the largest. ``auto`` scales the figure to the frame
    actually being measured; an explicit number still wins, so a value that
    was tuned by hand keeps working.
    """
    if value is None:
        return float(default)
    if isinstance(value, str):
        if value.strip().lower() not in ("auto", "frame"):
            return float(value)
        size = float(frame_arcmin) * float(fraction)
        # Downloads are cached under a key that includes the size, so a size
        # that drifts with every frame would mean a fresh download per frame.
        # Snapping to whole arcmin keeps one cache entry per detector.
        if step and step > 0:
            size = (math.ceil(size / step) if rounding == "up"
                    else math.floor(size / step)) * float(step)
        if lo is not None:
            size = max(size, float(lo))
        if hi is not None:
            size = min(size, float(hi))
        return size
    return float(value)


def frame_arcmin(image=None, detector: Optional[Detector] = None,
                 pixel_scale: Optional[float] = None, default: float = 12.3) -> float:
    """Width of a frame in arcmin, from an image or from a detector + scale."""
    if image is not None:
        try:
            return max(image.data.shape) * float(image.meta.pixel_scale) / 60.0
        except Exception:
            pass
    if detector is not None and pixel_scale:
        return max(detector.nx, detector.ny) * float(pixel_scale) / 60.0
    return float(default)


def stamp_arcmin(cfg, image=None, detector: Optional[Detector] = None,
                 pixel_scale: Optional[float] = None) -> float:
    """Working stamp size around the target, honouring ``auto``.

    Three quarters of the frame keeps a margin for the reprojection edges and
    for dithering, and the clamp stops a 6 arcmin subframe asking for a stamp
    wider than itself or a 21 arcmin frame dragging down a huge template.
    """
    return resolve_size(cfg.get("subtraction.cutout_arcmin", "auto"),
                        frame_arcmin(image, detector, pixel_scale),
                        fraction=0.75, lo=4.0, hi=12.0, default=9.0,
                        step=1.0, rounding="down")
