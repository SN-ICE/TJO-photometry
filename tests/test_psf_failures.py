"""Failed PSF fits must remain flagged in measurements and calibration."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from tjophot import diffquality, lightcurve as lc, photometry as ph, pipeline
from tjophot.calibration import ZeroPoint
from tjophot.config import Config
from tjophot.image import ImageMeta, ScienceImage


def noise_image(data):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10., 20.]
    wcs.wcs.crpix = [data.shape[1] / 2, data.shape[0] / 2]
    wcs.wcs.cdelt = [-1 / 3600., 1 / 3600.]
    meta = ImageMeta(Path("test.fits"), "SN2026test", "g", "ab", 60000.,
                     1., 4., 1., 60000., 60000., 1.)
    image = ScienceImage(data, fits.Header(), meta, wcs)
    image.bkg_map = np.zeros_like(data)
    image.bkg_rms = np.ones_like(data)
    image.error = np.ones_like(data)
    image.fwhm_px = image.fwhm_arcsec = 3.
    return image


def psf_fit(raw_flags=0, error=10., flux=1000.):
    return Table({"x": [50.], "y": [50.], "flux": [flux], "flux_err": [error],
                  "dx": [0.], "dy": [0.], "psf_flags": [raw_flags],
                  "qfit": [0.12], "cfit": [0.03]})


@pytest.mark.parametrize("raw_flags", [2, 8, 16, 64, 128, 256])
def test_photutils_failure_bits_are_explicitly_rejected(raw_flags):
    assert ph.psf_fit_failed(psf_fit(raw_flags)).tolist() == [True]


@pytest.mark.parametrize("raw_flags,flux", [(0, 1000.), (1, 1000.), (4, -10.), (32, 1000.)])
def test_photutils_diagnostics_are_not_automatically_failures(raw_flags, flux):
    assert ph.psf_fit_failed(psf_fit(raw_flags, flux=flux)).tolist() == [False]


@pytest.mark.parametrize("raw_flags,error", [(8, 10.), (16, np.nan), (0, 0.)])
def test_psf_failures_survive_both_measurement_paths_and_csv(monkeypatch, tmp_path,
                                                          raw_flags, error):
    image = noise_image(np.zeros((101, 101)))
    coord = image.wcs.pixel_to_world(50., 50.)
    cfg = Config({"output": {"write_qa_plots": False},
                  "photometry": {"empirical_errors": True},
                  "calibration": {"injection": {"enabled": False}}})
    for name in ("detect_sources", "select_psf_stars", "filter_by_profile"):
        monkeypatch.setattr(pipeline, name, lambda *a, **k: Table())
    monkeypatch.setattr(ph, "curve_of_growth", lambda *a, **k: SimpleNamespace(apcor_mag=0., apcor_err=0.))
    zp = ZeroPoint("g", 25., .01, 0., 0., .5, 20, .01, "test", "ab")
    monkeypatch.setattr(pipeline, "calibrate_image", lambda *a, **k: zp)
    monkeypatch.setattr(ph, "empirical_limit", lambda *a, **k: 23.)
    monkeypatch.setattr(ph, "blank_aperture_noise", lambda *a, **k: 10.)
    monkeypatch.setattr(ph, "measure_apertures", lambda *a, **k: Table({
        "flux": [1000.], "flux_err": [10.], "flags": [0]}))
    monkeypatch.setattr(ph, "build_epsf", lambda *a, **k: (object(), 10))
    monkeypatch.setattr(ph, "psf_to_aperture_correction", lambda *a, **k: 0.)
    monkeypatch.setattr(ph, "psf_photometry_at", lambda *a, **k: psf_fit(raw_flags, error))
    diff = SimpleNamespace(
        ok=True, difference=image.data, reference_stamp=image.data,
        variance=image.error**2, wcs=image.wcs, mode="convolve_reference",
        flux_scale=1., chi2=1., kernel_sum=1., n_stamps=10, meta={},
    )
    monkeypatch.setattr(pipeline.sub, "subtract", lambda *a, **k: diff)
    monkeypatch.setattr(diffquality, "assess", lambda *a, **k: SimpleNamespace(
        verdict="pass", message="", as_dict=lambda: {}))
    monkeypatch.setattr(diffquality, "host_residual", lambda *a, **k: {"sigma": 0., "n": 10})
    reference = SimpleNamespace(mjd_range=(np.nan, np.nan))
    rows = pipeline.process_image(image, cfg, coord, tmp_path, reference, tmp_path, "test")
    frame = lc.to_frame(rows)
    psf = frame[frame.method.isin(["psf", "difference_psf"])]
    assert len(psf) == 2
    assert (psf["flags"].astype(int) & lc.FLAG_BAD_PSF).ne(0).all()
    assert psf.psf_flags.eq(raw_flags).all()
    # A finite empirical error must not rehabilitate a failed fit.
    assert np.isfinite(psf.flux_err).all()
    assert frame.loc[frame.primary, "method"].tolist() == ["difference_aperture"]
    direct = [r for r in rows if not r["method"].startswith("difference")]
    assert lc.choose_primary(direct)["method"] == "aperture"
    paths = lc.write(frame, tmp_path, "test")
    stored = pd.read_csv(paths["all"])
    stored_psf = stored[stored.method.str.endswith("psf")]
    assert stored_psf.psf_flags.eq(raw_flags).all()
    assert stored_psf.psf_qfit.eq(.12).all()
    assert stored_psf.psf_cfit.eq(.03).all()
    assert stored_psf.flags_str.str.contains("bad_psf").all()


def test_failed_psf_fits_cannot_set_the_aperture_correction(monkeypatch):
    image = noise_image(np.zeros((101, 101)))
    stars = Table({"xcentroid": [50.], "ycentroid": [50.]})
    monkeypatch.setattr(ph, "measure_apertures", lambda *a, **k: Table({
        "flux": [500.], "flags": [0]}))
    monkeypatch.setattr(ph, "psf_photometry_at", lambda *a, **k: psf_fit(8))
    assert ph.psf_to_aperture_correction(image, object(), stars, 3., Config({})) == 0.


def test_failed_psf_cannot_be_restored_by_the_host_guard():
    difference = lc.new_row(method="difference_psf", mag=19., diff_chi2=1.,
                            flags=lc.FLAG_BAD_PSF | lc.FLAG_DIPOLE_RESIDUAL)
    direct = lc.new_row(method="aperture", mag=17., flags=0)
    assert lc.choose_primary([difference, direct]) is direct


def test_nonlinear_core_is_distinct_from_a_failed_psf():
    repaired = lc.new_row(method="psf", mag=17., flags=lc.FLAG_NONLINEAR_CORE)
    aperture = lc.new_row(method="aperture", mag=18., flags=0)
    assert lc.describe_flags(repaired["flags"]) == "nonlinear_core"
    assert lc.choose_primary([repaired, aperture]) is repaired
    repaired["flags"] |= lc.FLAG_BAD_PSF
    assert set(lc.describe_flags(repaired["flags"]).split("|")) == {"nonlinear_core", "bad_psf"}
    assert lc.choose_primary([repaired, aperture]) is aperture
