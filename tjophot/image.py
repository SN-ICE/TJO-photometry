"""Image ingest: FITS -> a fully-described science frame.

``ScienceImage`` is the single object the rest of the pipeline works with. It
owns the pixel data, the mask, the WCS, an uncertainty map and the metadata
extracted through the instrument header mapping.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u
from photutils.background import Background2D, MedianBackground

from .utils import first_header_value, get_logger

log = get_logger("tjophot.image")


@dataclass
class ImageMeta:
    """Everything the pipeline needs to know about a frame, instrument-neutral."""

    path: Path
    object_name: str
    band: str
    system: str                  # 'ab' or 'vega'
    mjd: float
    exptime: float
    gain: float                  # e-/ADU
    readnoise: float             # e-
    saturation: float            # ADU
    pixel_scale: float           # arcsec/px
    target_coord: Optional[SkyCoord] = None
    airmass: Optional[float] = None
    header_filter: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def epoch_id(self) -> str:
        return f"{self.mjd:.5f}_{self.band}"


def _sexagesimal_coord(ra_value, dec_value) -> Optional[SkyCoord]:
    """TJO writes RA in hours and Dec in degrees, both sexagesimal strings."""
    if ra_value is None or dec_value is None:
        return None
    try:
        if isinstance(ra_value, str) and (":" in ra_value or " " in ra_value.strip()):
            return SkyCoord(str(ra_value), str(dec_value), unit=(u.hourangle, u.deg))
        return SkyCoord(float(ra_value), float(dec_value), unit=(u.deg, u.deg))
    except Exception as exc:  # pragma: no cover - malformed headers
        log.warning("could not parse target coordinates %r %r: %s", ra_value, dec_value, exc)
        return None


class ScienceImage:
    """A single reduced TJO exposure, ready for photometry."""

    def __init__(self, data: np.ndarray, header: fits.Header, meta: ImageMeta,
                 wcs: Optional[WCS] = None):
        self.data = np.asarray(data, dtype=float)
        self.header = header
        self.meta = meta
        self.wcs = wcs
        self.mask = np.zeros(self.data.shape, dtype=bool)
        self._background = None      # Background2D result
        self.bkg_map: Optional[np.ndarray] = None
        self.bkg_rms: Optional[np.ndarray] = None
        self.error: Optional[np.ndarray] = None
        self.fwhm_px: Optional[float] = None
        self.fwhm_arcsec: Optional[float] = None
        self.astrometry = None       # AstrometryResult, once the WCS is refined
        self.qa: dict = {}

    # ------------------------------------------------------------------ load
    @classmethod
    def from_file(cls, path, cfg, instrument: dict, pattern=None) -> "ScienceImage":
        path = Path(path)
        with fits.open(path, memmap=False) as hdul:
            hdu = next((h for h in hdul if getattr(h, "data", None) is not None), hdul[0])
            data = np.asarray(hdu.data, dtype=float)
            header = hdu.header.copy()

        hmap = instrument.get("header", {})
        defaults = instrument.get("defaults", {})

        raw_filter = str(first_header_value(header, hmap.get("filter", ["FILTER"]), "")).strip()
        band = instrument.get("filters", {}).get(raw_filter)
        if band is None:
            band = instrument.get("filters", {}).get(raw_filter.replace("SDSS ", ""), raw_filter)
            log.warning("filter %r not in the instrument map; using %r", raw_filter, band)
        system = instrument.get("systems", {}).get(band, "ab")

        mjd = first_header_value(header, hmap.get("mjd", ["MJD-OBS"]))
        if mjd is None:
            jd = first_header_value(header, hmap.get("jd", ["JD"]))
            mjd = float(jd) - 2400000.5 if jd is not None else float("nan")

        binning = float(first_header_value(header, hmap.get("binning_x", ["XFACTOR"]), 1) or 1)

        wcs = None
        try:
            candidate = WCS(header)
            if candidate.has_celestial:
                wcs = candidate.celestial
        except Exception as exc:
            log.warning("%s: unusable WCS (%s)", path.name, exc)

        if wcs is not None:
            pixel_scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.0)
        else:
            pixel_scale = float(defaults.get("pixel_scale", 1.0)) * binning

        meta = ImageMeta(
            path=path,
            object_name=str(first_header_value(header, hmap.get("object", ["OBJECT"]), path.stem)).strip(),
            band=band,
            system=system,
            mjd=float(mjd),
            exptime=float(first_header_value(header, hmap.get("exptime", ["EXPTIME"]), 1.0)),
            gain=float(first_header_value(header, hmap.get("gain", ["GAIN"]), defaults.get("gain", 1.0))),
            readnoise=float(first_header_value(header, hmap.get("readnoise", ["RDNOISE"]),
                                               defaults.get("readnoise", 0.0))),
            saturation=float(cfg.get("image.saturation_adu", defaults.get("saturation_adu", 60000))),
            pixel_scale=pixel_scale,
            target_coord=_sexagesimal_coord(
                first_header_value(header, hmap.get("target_ra", ["OBJRA"])),
                first_header_value(header, hmap.get("target_dec", ["OBJDEC"])),
            ),
            airmass=first_header_value(header, hmap.get("airmass", ["AIRMASS"])),
            header_filter=raw_filter,
            extra={"quality": first_header_value(header, hmap.get("quality", ["IMGQS"]))},
        )

        image = cls(data, header, meta, wcs)

        if pattern is not None:
            from . import patterns as pattern_module
            fits_ok, why = pattern_module.geometry_matches(pattern, header, data.shape)
            if not fits_ok:
                log.warning("%s: not applying the %s pattern -- %s",
                            path.name, pattern.band, why)
                pattern = None

        if pattern is not None:
            from . import patterns as pattern_module
            try:
                image.data, info = pattern_module.apply_pattern(image.data, pattern, cfg)
                image.qa.update({f"pattern_{k}": v for k, v in info.items()})
                image.qa["pattern_applied"] = True
                log.debug("%s: pattern removed (fringe scale %.1f ADU)",
                          path.name, info.get("fringe_scale_adu", 0.0))
            except Exception as exc:
                log.warning("%s: pattern correction failed (%s); using the raw frame",
                            path.name, exc)
                image.qa["pattern_applied"] = False
        else:
            image.qa["pattern_applied"] = False

        image.build_mask(cfg)
        return image

    # ------------------------------------------------------------------ prep
    def build_mask(self, cfg) -> np.ndarray:
        """Mask non-finite, saturated, dead and edge pixels."""
        mask = ~np.isfinite(self.data)
        mask |= self.data >= self.meta.saturation
        for bad in cfg.get("image.bad_value_mask", []) or []:
            mask |= self.data == float(bad)

        margin = int(cfg.get("image.edge_margin_px", 0) or 0)
        if margin > 0:
            edge = np.zeros_like(mask)
            edge[:margin, :] = True
            edge[-margin:, :] = True
            edge[:, :margin] = True
            edge[:, -margin:] = True
            mask |= edge

        # Saturated stars bleed: grow the saturation mask slightly.
        if np.any(self.data >= self.meta.saturation):
            from scipy.ndimage import binary_dilation
            sat = binary_dilation(self.data >= self.meta.saturation, iterations=2)
            mask |= sat

        self.mask = mask
        self.qa["masked_fraction"] = float(mask.mean())
        return mask

    def estimate_background(self, cfg) -> None:
        """2D background and its RMS via photutils.Background2D."""
        box = int(cfg.get("image.background.box_size", 128))
        filt = int(cfg.get("image.background.filter_size", 5))
        nsig = float(cfg.get("image.background.sigma_clip", 3.0))
        try:
            bkg = Background2D(
                self.data,
                box_size=box,
                filter_size=filt,
                mask=self.mask,
                sigma_clip=SigmaClip(sigma=nsig),
                bkg_estimator=MedianBackground(),
                exclude_percentile=30.0,
            )
            self._background = bkg
            self.bkg_map = bkg.background
            self.bkg_rms = bkg.background_rms
        except Exception as exc:
            log.warning("%s: Background2D failed (%s); falling back to a constant sky",
                        self.meta.path.name, exc)
            _, median, std = sigma_clipped_stats(self.data, mask=self.mask, sigma=nsig)
            self.bkg_map = np.full(self.data.shape, median, dtype=float)
            self.bkg_rms = np.full(self.data.shape, std, dtype=float)

        self.qa["sky_adu"] = float(np.nanmedian(self.bkg_map))
        self.qa["sky_rms_adu"] = float(np.nanmedian(self.bkg_rms))

    def build_error_map(self) -> np.ndarray:
        """Per-pixel 1-sigma uncertainty in ADU.

        sigma^2 = sky_rms^2 (sky shot noise + read noise, measured empirically)
                  + source counts / gain  (Poisson noise of the source itself)
        """
        if self.bkg_rms is None:
            raise RuntimeError("call estimate_background() before build_error_map()")
        source = np.clip(self.data - self.bkg_map, 0.0, None)
        var = self.bkg_rms.astype(float) ** 2 + source / max(self.meta.gain, 1e-6)
        self.error = np.sqrt(var)
        return self.error

    # ------------------------------------------------------------- accessors
    @property
    def data_sub(self) -> np.ndarray:
        """Background-subtracted data."""
        if self.bkg_map is None:
            raise RuntimeError("background not estimated yet")
        return self.data - self.bkg_map

    def pixel_of(self, coord: SkyCoord) -> tuple[float, float]:
        if self.wcs is None:
            raise RuntimeError(f"{self.meta.path.name} has no WCS")
        x, y = self.wcs.world_to_pixel(coord)
        return float(x), float(y)

    def contains(self, coord: SkyCoord, margin: int = 0) -> bool:
        try:
            x, y = self.pixel_of(coord)
        except Exception:
            return False
        ny, nx = self.data.shape
        return (margin <= x < nx - margin) and (margin <= y < ny - margin)

    def center_coord(self) -> SkyCoord:
        ny, nx = self.data.shape
        return self.wcs.pixel_to_world(nx / 2.0, ny / 2.0)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"<ScienceImage {self.meta.path.name} band={self.meta.band} "
                f"mjd={self.meta.mjd:.4f} fwhm={self.fwhm_arcsec}>")
