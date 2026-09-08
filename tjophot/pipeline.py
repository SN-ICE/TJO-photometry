"""Orchestration: from a directory of FITS files to a calibrated light curve.

The run is deliberately memory-frugal. Images are never all held at once:

  phase 0  read headers only and build an index of the epochs
  phase 1  build (or load from cache) one reference image per band
  phase 2  process one science image at a time, writing a row per method
  phase 3  assemble, plot, write
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
import astropy.units as u

from . import astrometry as astrom
from . import catalogs as cat
from . import lightcurve as lc
from . import nameresolve
from . import photometry as ph
from . import plots
from . import zeropoints
from . import patterns as pat
from . import templates as tpl
from . import subtraction as sub
from .calibration import ZeroPoint, calibrate_image
from .detection import detect_sources, filter_by_profile, measure_fwhm, select_psf_stars
from .image import ScienceImage
from .utils import ensure_dir, first_header_value, get_logger, slugify

log = get_logger("tjophot.pipeline")


# --------------------------------------------------------------------------
@dataclass
class EpochIndex:
    path: Path
    mjd: float
    band: str
    object_name: str
    exptime: float


def index_images(paths: Iterable[Path], instrument: dict,
                 target: Optional[str] = None) -> list[EpochIndex]:
    """Read headers only and build the epoch list."""
    hmap = instrument.get("header", {})
    filters = instrument.get("filters", {})
    index: list[EpochIndex] = []
    for path in sorted(paths):
        try:
            header = fits.getheader(path)
        except Exception as exc:
            log.warning("cannot read %s: %s", path, exc)
            continue
        raw = str(first_header_value(header, hmap.get("filter", ["FILTER"]), "")).strip()
        band = filters.get(raw, filters.get(raw.replace("SDSS ", ""), raw))
        mjd = first_header_value(header, hmap.get("mjd", ["MJD-OBS"]))
        if mjd is None:
            jd = first_header_value(header, hmap.get("jd", ["JD"]))
            mjd = float(jd) - 2400000.5 if jd is not None else float("nan")
        name = str(first_header_value(header, hmap.get("object", ["OBJECT"]), "")).strip()
        if target and name and name.lower() != target.lower():
            continue
        index.append(EpochIndex(Path(path), float(mjd), band, name,
                                float(first_header_value(header, hmap.get("exptime", ["EXPTIME"]), 1.0))))
    return sorted(index, key=lambda e: (e.band, e.mjd))


def find_images(cfg, target: str) -> list[Path]:
    """Locate the FITS files for *target* under ``paths.data_root``."""
    root = Path(cfg.data_root).expanduser()
    patterns = cfg.get("paths.glob", ["*.fits", "*.fit", "*.fits.gz"])
    candidates: list[Path] = []
    for directory in (root / target, root):
        if not directory.is_dir():
            continue
        for pattern in patterns:
            candidates.extend(sorted(directory.glob(pattern)))
        if candidates:
            break
    return sorted(set(candidates))


def resolve_target_coord(cfg, index: Sequence[EpochIndex], instrument: dict) -> tuple:
    """SN position. Returns (coord, was_explicit).

    An explicit ``target.ra``/``target.dec`` wins. Otherwise the telescope
    pointing (OBJRA/OBJDEC) is used as a starting guess -- on TJO this is
    typically about an arcsecond from the supernova, which is enough to bias a
    1.2 x FWHM aperture, so it should be refined (see ``refine_target_coord``).
    """
    ra = cfg.get("target.ra")
    dec = cfg.get("target.dec")
    if ra is not None and dec is not None:
        if isinstance(ra, str) and (":" in ra or " " in str(ra).strip()):
            return SkyCoord(str(ra), str(dec), unit=(u.hourangle, u.deg)), True
        return SkyCoord(float(ra), float(dec), unit="deg"), True

    hmap = instrument.get("header", {})
    ras, decs = [], []
    for entry in index[: min(len(index), 20)]:
        header = fits.getheader(entry.path)
        rv = first_header_value(header, hmap.get("target_ra", ["OBJRA"]))
        dv = first_header_value(header, hmap.get("target_dec", ["OBJDEC"]))
        if rv is None or dv is None:
            continue
        try:
            coord = SkyCoord(str(rv), str(dv), unit=(u.hourangle, u.deg))
        except Exception:
            continue
        ras.append(coord.ra.deg)
        decs.append(coord.dec.deg)
    if not ras:
        raise RuntimeError("cannot determine the target position; set target.ra / target.dec")
    coord = SkyCoord(float(np.median(ras)) * u.deg, float(np.median(decs)) * u.deg)
    log.info("starting from the telescope pointing (%s)", coord.to_string("hmsdms", precision=1))
    return coord, False


def refine_target_coord(coord: SkyCoord, cfg, index: Sequence[EpochIndex], instrument: dict,
                        n_epochs: int = 5, patterns: Optional[dict] = None,
                        astrometry_ref=None) -> tuple:
    """Centroid the transient near *coord* on a few epochs and adopt the median.

    Returns (coord, offset_arcsec). The refinement is rejected if it wants to
    move further than ``target.refine_radius_arcsec``, which stops it locking
    onto a neighbouring star when the supernova is undetected.
    """
    from photutils.centroids import centroid_quadratic

    max_shift = float(cfg.get("target.refine_radius_arcsec", 3.0))
    positions, offsets = [], []

    # Sample across the whole campaign rather than taking the first few frames,
    # which would all come from one band and one night.
    stride = max(1, len(index) // max(n_epochs * 4, 1))
    for entry in list(index)[::stride]:
        if len(positions) >= n_epochs:
            break
        try:
            image = prepare_image(entry.path, cfg, instrument,
                                  pattern=patterns.get(entry.band),
                                  astrometry_ref=astrometry_ref)
            if not image.contains(coord, margin=40):
                continue
            fwhm = image.fwhm_px or 4.0
            x0, y0 = image.pixel_of(coord)
            half = int(np.ceil(max(2.0 * fwhm, max_shift / image.meta.pixel_scale)))
            xi, yi = int(round(x0)), int(round(y0))
            stamp = image.data_sub[yi - half:yi + half + 1, xi - half:xi + half + 1]
            if stamp.size == 0 or not np.isfinite(stamp).all():
                continue
            peak = float(np.nanmax(stamp))
            if peak < 10.0 * float(np.nanmedian(image.bkg_rms)):
                continue                       # nothing convincing to centre on
            cx, cy = centroid_quadratic(stamp)
            if not (np.isfinite(cx) and np.isfinite(cy)):
                continue
            dx, dy = cx - half, cy - half
            shift = float(np.hypot(dx, dy) * image.meta.pixel_scale)
            if shift > max_shift:
                continue
            positions.append(image.wcs.pixel_to_world(x0 + dx, y0 + dy))
            offsets.append(shift)
        except Exception as exc:
            log.debug("position refinement skipped %s: %s", entry.path.name, exc)

    if len(positions) < 2:
        log.warning("could not refine the target position; using it as given")
        return coord, float("nan")

    refined = SkyCoord(np.median([c.ra.deg for c in positions]) * u.deg,
                       np.median([c.dec.deg for c in positions]) * u.deg)
    offset = float(coord.separation(refined).arcsec)
    log.info("refined the target position by %.2f\" using %d epochs -> %s",
             offset, len(positions), refined.to_string("hmsdms", precision=2))
    return refined, offset


# --------------------------------------------------------------------------
def prepare_image(path: Path, cfg, instrument: dict, pattern=None,
                  astrometry_ref=None) -> ScienceImage:
    """Load a frame and get it ready for photometry.

    Order matters: the detector pattern is removed before the background is
    measured (otherwise the fringes end up in the background model), and the
    WCS is refined before anything places an aperture through it.
    """
    image = ScienceImage.from_file(path, cfg, instrument, pattern=pattern)
    image.estimate_background(cfg)
    image.build_error_map()
    measure_fwhm(image, cfg)

    if astrometry_ref is not None and str(cfg.get("astrometry.mode", "linear")).lower() \
            not in ("off", "none"):
        try:
            sources = detect_sources(image, cfg, fwhm_px=image.fwhm_px)
            result = astrom.refine_wcs(image, sources, astrometry_ref, cfg)
            image.astrometry = result
            if result.ok and result.wcs is not None:
                image.wcs = result.wcs
                image.qa.update({"astrom_rms_arcsec": result.rms_after,
                                 "astrom_rms_before_arcsec": result.rms_before,
                                 "astrom_nstars": result.n_stars,
                                 "astrom_shift_arcsec": result.shift_arcsec,
                                 "astrom_mode": result.mode})
            else:
                image.qa["astrom_message"] = result.message
        except Exception as exc:
            log.warning("%s: WCS refinement failed (%s)", path.name, exc)

    return image


def load_patterns(cfg, bands: Sequence[str]) -> dict:
    """Detector patterns for the bands that have one stored."""
    if not bool(cfg.get("image.pattern.enabled", True)):
        return {}
    calib_root = _calib_root(cfg)
    out = {}
    for band in bands:
        pattern = pat.load_pattern(calib_root, band)
        if pattern is not None:
            out[band] = pattern
            log.info("%s: using stored detector pattern (%d frames, illumination %.2f%% pp, "
                     "fringe %.2f%% pp)", band, pattern.n_frames,
                     100 * pattern.illumination_amplitude, 100 * pattern.fringe_amplitude)
    if not out:
        log.info("no stored detector patterns; run 'tjophot build-pattern' to make them")
    return out


def _calib_root(cfg) -> Path:
    from .config import PIPELINE_ROOT
    value = cfg.get("image.pattern.calib_root", "calib")
    path = Path(value).expanduser()
    return path if path.is_absolute() else (PIPELINE_ROOT / path)


def load_astrometry_reference(cfg, coord: SkyCoord, cache_root: Path,
                              refresh: bool = False):
    """Reference catalogue for the WCS refinement, with fallbacks."""
    mode = str(cfg.get("astrometry.mode", "linear")).lower()
    if mode in ("off", "none"):
        return None

    providers = [cfg.get("astrometry.catalog", "gaia_astrometry")]
    providers += list(cfg.get("astrometry.fallback_catalogs", []) or [])
    radius = float(cfg.get("astrometry.search_radius_arcmin", 10.0))
    local_path = cfg.resolve_path(cfg.get("calibration.local_catalog"))

    for provider in providers:
        try:
            table = cat.get_catalog(provider, float(coord.ra.deg), float(coord.dec.deg),
                                    radius, cache_root, refresh=refresh,
                                    local_path=local_path)
            if table is not None and len(table) >= int(cfg.get("astrometry.min_stars", 8)):
                log.info("astrometric reference: %s (%d sources)", provider, len(table))
                return table
        except Exception as exc:
            log.warning("astrometric reference %s unavailable: %s", provider, exc)
    log.warning("no astrometric reference catalogue; keeping the header WCS")
    return None


def build_references(cfg, index: Sequence[EpochIndex], instrument: dict, coord: SkyCoord,
                     target: str, cache_root: Path, bands: Sequence[str],
                     refresh: bool = False, patterns: Optional[dict] = None,
                     astrometry_ref=None) -> dict:
    """One reference image per band, built once and reused for every epoch."""
    references: dict[str, Optional[tpl.ReferenceImage]] = {}
    max_images = int(cfg.get("subtraction.tjo_stack.max_images", 8))
    mjd_min = cfg.get("subtraction.tjo_stack.mjd_min")
    mjd_max = cfg.get("subtraction.tjo_stack.mjd_max")

    for band in bands:
        order = cfg.get(f"subtraction.reference_sources.{band}") or \
            ["file", "tjo_stack", "panstarrs"]
        siblings: list[ScienceImage] = []
        if "tjo_stack" in order:
            if mjd_min is None and mjd_max is None:
                log.info("band %s: no template MJD window set, skipping tjo_stack "
                         "(set subtraction.tjo_stack.mjd_min/mjd_max to use the "
                         "target's own late-time epochs)", band)
            else:
                selection = [e for e in index if e.band == band
                             and (mjd_min is None or e.mjd >= float(mjd_min))
                             and (mjd_max is None or e.mjd <= float(mjd_max))]
                if len(selection) > max_images:
                    selection = selection[:max_images]
                for entry in selection:
                    try:
                        siblings.append(prepare_image(
                            entry.path, cfg, instrument,
                            pattern=(patterns or {}).get(band),
                            astrometry_ref=astrometry_ref))
                    except Exception as exc:
                        log.warning("template epoch %s failed: %s", entry.path.name, exc)

        try:
            references[band] = tpl.get_reference(
                _FakeImage(band), cfg, cache_root, coord, siblings, target, refresh=refresh)
        except Exception as exc:
            log.warning("reference for band %s failed: %s", band, exc)
            references[band] = None
        siblings.clear()
    return references


class _FakeImage:
    """Minimal stand-in so templates.get_reference can be called before any
    science frame is loaded (it only needs the band and the MJD to exclude)."""

    class _Meta:
        def __init__(self, band):
            self.band = band
            self.mjd = float("nan")

    def __init__(self, band: str):
        self.meta = self._Meta(band)


# --------------------------------------------------------------------------
def process_image(image: ScienceImage, cfg, coord: SkyCoord, cache_root: Path,
                  reference: Optional[tpl.ReferenceImage], outdir: Path,
                  target: str, refresh: bool = False) -> list[dict]:
    """All measurements for one frame: returns one row per photometric method."""
    band = image.meta.band
    fwhm = image.fwhm_px or 4.0
    ann = cfg.get("photometry.annulus_scale", [5.0, 8.0])
    annulus = (ann[0] * fwhm, ann[1] * fwhm)
    ap_radius = float(cfg.get("photometry.aperture_scale", 1.2)) * fwhm

    common = dict(
        target=target, mjd=image.meta.mjd, band=band, system=image.meta.system,
        fwhm_arcsec=image.fwhm_arcsec, sky_adu=image.qa.get("sky_adu"),
        sky_rms_adu=image.qa.get("sky_rms_adu"), exptime=image.meta.exptime,
        airmass=image.meta.airmass, aperture_radius_px=ap_radius,
        quality=str(image.meta.extra.get("quality", "")), image=image.meta.path.name,
        astrom_rms_arcsec=image.qa.get("astrom_rms_arcsec", np.nan),
        astrom_nstars=image.qa.get("astrom_nstars", np.nan),
        pattern_applied=bool(image.qa.get("pattern_applied", False)),
        fringe_scale_adu=image.qa.get("pattern_fringe_scale_adu", np.nan),
    )

    sources = detect_sources(image, cfg, fwhm_px=fwhm)
    psf_stars = select_psf_stars(image, sources, cfg, fwhm,
                                 max_stars=int(cfg.get("detection.max_stars_for_fwhm", 40)),
                                 min_snr=float(cfg.get("calibration.min_snr", 20.0)))
    psf_stars = filter_by_profile(image, psf_stars, fwhm, cfg)
    apcor = ph.curve_of_growth(image, psf_stars, fwhm, cfg)
    zp = calibrate_image(image, sources, fwhm, cfg, cache_root,
                         apcor_mag=apcor.apcor_mag, refresh=refresh)

    common.update(
        apcor=apcor.apcor_mag, apcor_err=apcor.apcor_err,
        zp=zp.zp, zp_err=zp.zp_err, color_term=zp.color_term,
        color_term_err=zp.color_term_err, cal_provider=zp.provider,
        cal_nstars=zp.n_stars, cal_rms=zp.rms,
    )

    limit = ph.empirical_limit(image, fwhm, zp.zp, cfg, apcor.apcor_mag) if zp.ok else np.nan

    cal_flags = 0
    if reference is not None:
        lo, hi = reference.mjd_range
        if np.isfinite(lo) and np.isfinite(hi) and lo - 1e-6 <= image.meta.mjd <= hi + 1e-6:
            cal_flags |= lc.FLAG_IN_TEMPLATE
            log.warning("%s is itself part of the reference image; its difference is null "
                        "by construction", image.meta.path.name)
    if not zp.ok:
        cal_flags |= lc.FLAG_NO_ZEROPOINT
    else:
        if zp.n_stars < 2 * int(cfg.get("calibration.min_stars", 5)):
            cal_flags |= lc.FLAG_FEW_CAL_STARS
        if np.isfinite(zp.rms) and zp.rms > float(cfg.get("calibration.max_rms", 0.25)):
            cal_flags |= lc.FLAG_POOR_CALIBRATION
    common["flags_base"] = cal_flags

    sn_xy = image.pixel_of(coord)
    rows: list[dict] = []

    # ---------------------------------------------------------------- direct
    tab = ph.measure_apertures(image, [sn_xy], ap_radius, cfg, annulus_px=annulus)
    sky_sigma = ph.blank_aperture_noise(image, ap_radius, cfg) \
        if bool(cfg.get("photometry.empirical_errors", True)) else float("nan")
    aper_err = ph.combine_errors(float(tab["flux"][0]), float(tab["flux_err"][0]),
                                 sky_sigma, image.meta.gain)
    rows.append(_photometry_row(common, "aperture", float(tab["flux"][0]),
                                aper_err, zp, apcor.apcor_mag,
                                image, x=sn_xy[0], y=sn_xy[1],
                                flags=int(tab["flags"][0]), limit=limit))

    epsf, n_psf = ph.build_epsf(image, psf_stars, fwhm, cfg)
    psf_corr = 0.0
    if epsf is not None and bool(cfg.get("photometry.psf.enabled", True)):
        psf_corr = ph.psf_to_aperture_correction(image, epsf, psf_stars, fwhm, cfg,
                                                 apcor_mag=apcor.apcor_mag)
        fit = ph.psf_photometry_at(image, epsf, [sn_xy], fwhm, cfg)
        if len(fit):
            psf_err = ph.combine_errors(float(fit["flux"][0]), float(fit["flux_err"][0]),
                                        sky_sigma, image.meta.gain)
            rows.append(_photometry_row(common, "psf", float(fit["flux"][0]),
                                        psf_err, zp, psf_corr, image,
                                        x=float(fit["x"][0]), y=float(fit["y"][0]),
                                        dx=float(fit["dx"][0]), dy=float(fit["dy"][0]),
                                        flags=int(tab["flags"][0]), limit=limit))

    # ------------------------------------------------------------ difference
    diff = None
    if reference is not None and bool(cfg.get("subtraction.enabled", True)):
        diff = sub.subtract(image, reference, cfg, coord, fwhm_px=fwhm)
        if diff.ok:
            rows.extend(_difference_rows(diff, image, cfg, coord, common, zp, limit,
                                         science_epsf=epsf, science_apcor=apcor.apcor_mag,
                                         science_psf_corr=psf_corr))

    if bool(cfg.get("output.write_qa_plots", True)):
        plots.epoch_qa(Path(outdir) / "qa" /
                       f"{slugify(image.meta.path.stem)}.png",
                       image, zp, apcor, diff=diff, sn_xy=sn_xy, cutout=coord,
                       stamp_px=int(float(cfg.get("output.stamp_arcsec", 60.0)) /
                                    max(image.meta.pixel_scale, 1e-3) / 2))

    preference = str(cfg.get("output.primary_method", "auto"))
    best = lc.choose_primary(rows, preference)
    if best is not None:
        best["primary"] = True
    return rows


def _difference_rows(diff, image, cfg, coord, common, zp, limit,
                     science_epsf=None, science_apcor: float = 0.0,
                     science_psf_corr: float = 0.0) -> list[dict]:
    """Aperture and PSF photometry on the difference image.

    When the reference was the frame that got convolved, the difference keeps
    the science PSF exactly, so the PSF model and aperture correction already
    measured on the full science frame are reused -- far more reliable than
    rebuilding them from the handful of stars inside the stamp.
    """
    rows: list[dict] = []
    stamp = _DifferenceImage(diff, image)
    fwhm = _matched_fwhm(diff, image)
    ann = cfg.get("photometry.annulus_scale", [5.0, 8.0])
    annulus = (ann[0] * fwhm, ann[1] * fwhm)
    ap_radius = float(cfg.get("photometry.aperture_scale", 1.2)) * fwhm

    x, y = diff.wcs.world_to_pixel(coord)
    sn_xy = (float(x), float(y))

    max_chi2 = float(cfg.get("subtraction.max_chi2", 2.5))
    bad_subtraction = lc.FLAG_BAD_SUBTRACTION if (
        np.isfinite(diff.chi2) and diff.chi2 > max_chi2) else 0
    ref_band = diff.meta.get("reference_band")
    if ref_band and ref_band != image.meta.band:
        bad_subtraction |= lc.FLAG_CROSS_BAND_TEMPLATE

    diff_meta = dict(
        ref_source=diff.meta.get("reference_source"),
        ref_band=diff.meta.get("reference_band"),
        ref_nimages=diff.meta.get("reference_n_images"),
        diff_mode=diff.mode, diff_chi2=diff.chi2,
        diff_kernel_sum=diff.kernel_sum, diff_nstamps=diff.n_stamps,
    )

    # Aperture correction measured on the matched image, which shares the PSF
    # of the difference image.
    matched = _MatchedImage(diff, image)
    reuse_science_psf = (diff.mode == "convolve_reference" and science_epsf is not None)

    matched_stars, apcor = None, None
    if reuse_science_psf:
        apcor_mag = float(science_apcor)
    else:
        from .detection import detect_sources as _detect, select_psf_stars as _select
        try:
            matched_sources = _detect(matched, cfg, fwhm_px=fwhm)
            matched_stars = _select(matched, matched_sources, cfg, fwhm, max_stars=30,
                                    min_snr=20.0)
            matched_stars = filter_by_profile(matched, matched_stars, fwhm, cfg)
            apcor = ph.curve_of_growth(matched, matched_stars, fwhm, cfg)
            if apcor.n_stars < int(cfg.get("photometry.min_apcor_stars", 5)):
                apcor = None      # too few stars to trust; fall back to no correction
        except Exception:
            matched_stars, apcor = None, None
        apcor_mag = apcor.apcor_mag if apcor is not None else 0.0
    scale = diff.flux_scale if np.isfinite(diff.flux_scale) and diff.flux_scale > 0 else 1.0

    tab = ph.measure_apertures(stamp, [sn_xy], ap_radius, cfg, annulus_px=annulus)
    sky_sigma = ph.blank_aperture_noise(stamp, ap_radius, cfg) \
        if bool(cfg.get("photometry.empirical_errors", True)) else float("nan")
    diff_err = ph.combine_errors(float(tab["flux"][0]), float(tab["flux_err"][0]),
                                 sky_sigma, image.meta.gain * scale)
    rows.append(_photometry_row({**common, **diff_meta}, "difference_aperture",
                                float(tab["flux"][0]) / scale,
                                diff_err / scale, zp, apcor_mag, image,
                                x=sn_xy[0], y=sn_xy[1],
                                flags=int(tab["flags"][0]) | bad_subtraction,
                                limit=limit))

    if not bool(cfg.get("photometry.psf.enabled", True)):
        return rows

    if reuse_science_psf:
        epsf, psf_corr = science_epsf, float(science_psf_corr)
    elif matched_stars is not None and len(matched_stars):
        epsf, _ = ph.build_epsf(matched, matched_stars, fwhm, cfg)
        psf_corr = ph.psf_to_aperture_correction(matched, epsf, matched_stars, fwhm, cfg,
                                                 apcor_mag=apcor_mag) if epsf is not None else 0.0
    else:
        epsf, psf_corr = None, 0.0

    if epsf is not None:
        fit = ph.psf_photometry_at(stamp, epsf, [sn_xy], fwhm, cfg)
        if len(fit):
            dpsf_err = ph.combine_errors(float(fit["flux"][0]), float(fit["flux_err"][0]),
                                         sky_sigma, image.meta.gain * scale)
            rows.append(_photometry_row({**common, **diff_meta}, "difference_psf",
                                        float(fit["flux"][0]) / scale,
                                        dpsf_err / scale, zp, psf_corr,
                                        image, x=float(fit["x"][0]), y=float(fit["y"][0]),
                                        dx=float(fit["dx"][0]), dy=float(fit["dy"][0]),
                                        flags=int(tab["flags"][0]) | bad_subtraction,
                                        limit=limit))
    return rows


def _matched_fwhm(diff, image) -> float:
    value = diff.meta.get("reference_fwhm_px", np.nan)
    science = image.fwhm_px or 4.0
    if diff.mode == "convolve_science" and np.isfinite(value):
        return float(value)
    return float(science)


class _StampBase:
    """Adapter presenting a stamp array with the ScienceImage interface."""

    def __init__(self, data, variance, wcs, image, mask=None):
        self.data = np.nan_to_num(data, nan=0.0)
        self.wcs = wcs
        self.meta = image.meta
        self.header = image.header
        self.mask = (~np.isfinite(data)) if mask is None else mask
        self.bkg_map = np.zeros_like(self.data)
        with np.errstate(invalid="ignore"):
            rms = np.sqrt(np.abs(variance))
        median_rms = float(np.nanmedian(rms[np.isfinite(rms)])) if np.isfinite(rms).any() else 1.0
        self.bkg_rms = np.where(np.isfinite(rms) & (rms > 0), rms, median_rms)
        self.error = self.bkg_rms.copy()
        self.fwhm_px = image.fwhm_px
        self.fwhm_arcsec = image.fwhm_arcsec
        self.qa = {"sky_adu": 0.0, "sky_rms_adu": median_rms}

    @property
    def data_sub(self):
        return self.data

    def build_error_map(self):
        return self.error

    def pixel_of(self, coord):
        x, y = self.wcs.world_to_pixel(coord)
        return float(x), float(y)

    def center_coord(self):
        ny, nx = self.data.shape
        return self.wcs.pixel_to_world(nx / 2.0, ny / 2.0)


class _DifferenceImage(_StampBase):
    def __init__(self, diff, image):
        super().__init__(diff.difference, diff.variance, diff.wcs, image)


class _MatchedImage(_StampBase):
    def __init__(self, diff, image):
        super().__init__(diff.reference_stamp, diff.variance, diff.wcs, image)


def _photometry_row(common: dict, method: str, flux: float, flux_err: float,
                    zp: ZeroPoint, apcor_mag: float, image, x=np.nan, y=np.nan,
                    dx=np.nan, dy=np.nan, flags: int = 0, limit: float = np.nan) -> dict:
    """Turn a flux measurement into a calibrated row."""
    common = dict(common)
    flags = int(flags) | int(common.pop("flags_base", 0))
    exptime = max(image.meta.exptime, 1e-6)
    snr = flux / flux_err if (np.isfinite(flux_err) and flux_err > 0) else np.nan

    mag = np.nan
    mag_err = np.nan
    is_limit = True
    if zp.ok and np.isfinite(flux) and flux > 0 and np.isfinite(snr) and snr >= 3.0:
        inst = -2.5 * np.log10(flux / exptime) + apcor_mag
        mag = inst + zp.zp
        mag_err = float(np.hypot(1.0857362 * flux_err / flux, zp.zp_err))
        is_limit = False

    return lc.new_row(
        **common, method=method, mag=mag, mag_err=mag_err, is_limit=is_limit,
        mag_offset=apcor_mag, zp_source=("fit" if zp.ok else ""),
        flux=flux / exptime, flux_err=flux_err / exptime, snr=snr,
        x=x, y=y, dx=dx, dy=dy, flags=flags, limit_mag=limit,
    )


# --------------------------------------------------------------------------
def run(cfg, target: str, bands: Optional[Sequence[str]] = None,
        max_images: Optional[int] = None, refresh: bool = False,
        skip_subtraction: bool = False, append: bool = False,
        skip: int = 0) -> dict:
    """Full pipeline for one target. Returns a dict of output paths."""
    started = time.time()
    instrument = cfg.get("instrument_config", {})
    outdir = ensure_dir(Path(cfg.work_root) / target)
    cache_root = ensure_dir(cfg.cache_root)

    paths = find_images(cfg, target)
    if not paths:
        raise RuntimeError(f"no FITS files found for {target} under {cfg.data_root}")
    full_index = index_images(paths, instrument, target=None)
    if bands:
        full_index = [e for e in full_index if e.band in bands]
    # Templates are always built from every available epoch, even when the run
    # itself is restricted with --max-images.
    index = full_index[int(skip):] if skip else list(full_index)
    if max_images:
        index = index[:max_images]
    if not index:
        raise RuntimeError("no epochs left after filtering")

    present = sorted({e.band for e in index})
    log.info("%s: %d images, bands %s, MJD %.3f-%.3f", target, len(index), ",".join(present),
             min(e.mjd for e in index), max(e.mjd for e in index))

    coord, explicit = resolve_target_coord(cfg, full_index, instrument)

    mode = str(cfg.get("target.resolve", "auto")).lower()
    if mode == "always" or (mode == "auto" and not explicit):
        found = nameresolve.resolve(target, cfg)
        if found is not None:
            if explicit:
                moved = coord.separation(found.coord).arcsec
                log.info("%s moves the position by %.2f\"", found.source.upper(), moved)
            coord, explicit = found.coord, True
            if bool(cfg.get("target.resolve_write", True)):
                try:
                    nameresolve.update_target_file(target, found)
                except Exception as exc:
                    log.warning("could not update the target file: %s", exc)
        elif not explicit:
            log.info("name resolution found nothing; falling back to the telescope pointing")

    patterns = load_patterns(cfg, present)
    astrometry_ref = load_astrometry_reference(cfg, coord, cache_root, refresh=refresh)

    refine = str(cfg.get("target.refine", "auto")).lower()
    offset = float("nan")
    if refine == "always" or (refine == "auto" and not explicit):
        coord, offset = refine_target_coord(coord, cfg, full_index, instrument,
                                            patterns=patterns, astrometry_ref=astrometry_ref)
    elif not explicit:
        log.warning("using the telescope pointing as the SN position without refinement; "
                    "set target.ra / target.dec for precise forced photometry")
    log.info("SN position: %s", coord.to_string("hmsdms", precision=2))

    references: dict = {}
    if bool(cfg.get("subtraction.enabled", True)) and not skip_subtraction:
        references = build_references(cfg, full_index, instrument, coord, target, cache_root,
                                      present, refresh=refresh, patterns=patterns,
                                      astrometry_ref=astrometry_ref)

    rows: list[dict] = []
    failures: list[dict] = []
    for k, entry in enumerate(index, 1):
        log.info("[%d/%d] %s", k, len(index), entry.path.name)
        try:
            image = prepare_image(entry.path, cfg, instrument,
                                  pattern=patterns.get(entry.band),
                                  astrometry_ref=astrometry_ref)
            if not image.contains(coord, margin=int(cfg.get("image.edge_margin_px", 24))):
                raise RuntimeError("the target position falls outside this frame")
            rows.extend(process_image(image, cfg, coord, cache_root,
                                      references.get(entry.band), outdir, target,
                                      refresh=refresh))
        except Exception as exc:
            log.error("%s failed: %s", entry.path.name, exc)
            failures.append({"image": entry.path.name, "error": str(exc)})
        finally:
            image = None

    if not rows:
        raise RuntimeError("no photometry could be measured")

    frame = lc.to_frame(rows)

    if bool(cfg.get("calibration.transfer_zeropoints", True)):
        frame, transfer_report = zeropoints.transfer(frame, cfg)
        if transfer_report.get("recovered"):
            frame = lc.reselect_primary(frame, str(cfg.get("output.primary_method", "auto")))
    else:
        transfer_report = {}

    paths_out = lc.write(frame, outdir, target, append=append,
                         preference=str(cfg.get("output.primary_method", "auto")))
    if append:
        import pandas as pd
        frame = pd.read_csv(paths_out["all"])

    if bool(cfg.get("output.write_qa_plots", True)):
        plots.lightcurve(outdir / f"{target}_lightcurve.png", frame, target)
        plots.run_diagnostics(outdir / f"{target}_diagnostics.png", frame, target)

    summary = {
        "target": target,
        "config_files": cfg.get("_config_files", []),
        "n_images": len(index),
        "n_rows": len(frame),
        "n_failures": len(failures),
        "failures": failures,
        "bands": present,
        "coord": coord.to_string("hmsdms", precision=3),
        "position_refined_by_arcsec": None if not np.isfinite(offset) else round(offset, 3),
        "position_from_config": bool(explicit),
        "elapsed_s": round(time.time() - started, 1),
        "zeropoint_transfer": {
            "recovered": len(transfer_report.get("recovered", [])),
            "unrecovered": len(transfer_report.get("unrecovered", [])),
            "band_offsets": transfer_report.get("offsets", {}),
        },
        "references": {
            band: (None if ref is None else {
                "source": ref.source,
                "band": ref.source_band or ref.band,
                "path": str(ref.path) if ref.path else None,
                "pixel_scale_arcsec": round(float(ref.pixel_scale), 4),
                "n_images": ref.n_images,
            })
            for band, ref in references.items()
        },
    }
    with open(outdir / "run_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    cfg.dump(outdir / "config_used.yaml")

    log.info("done in %.1f s -- %s", summary["elapsed_s"], outdir)
    print("\n" + lc.summarise(frame) + "\n")
    return {**paths_out, "outdir": outdir, "summary": summary, "frame": frame}
