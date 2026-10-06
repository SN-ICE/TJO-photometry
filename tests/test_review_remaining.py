"""Regression coverage for discovery, quality selection and noise calibration."""

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from photutils.aperture import CircularAperture, aperture_photometry
from scipy.ndimage import gaussian_filter

from tjophot import lightcurve as lc, photometry as ph, pipeline
from tjophot.config import Config
from tjophot.image import ImageMeta, ScienceImage


def test_discovery_includes_header_matches_when_filenames_also_match(tmp_path):
    entries = {
        "TJO_2026test_g.fits": "2026test",
        "renamed.fits": "SN 2026test",
        "nested/exposure.fits": "AT2026test",
        "TJO_2026test_wrong.fits": "2026other",
        "unrelated.fits": "2026other",
    }
    for filename, target in entries.items():
        path = tmp_path / filename
        path.parent.mkdir(exist_ok=True)
        fits.PrimaryHDU(np.zeros((4, 4)), fits.Header({"OBJECT": target})).writeto(path)
    (tmp_path / "broken.fits").write_text("not a FITS file")
    cfg = Config({"paths": {"data_root": str(tmp_path)}})
    found = pipeline.find_images(cfg, "SN2026test")
    assert {p.relative_to(tmp_path).as_posix() for p in found} == {
        "TJO_2026test_g.fits", "renamed.fits", "nested/exposure.fits"}
    assert len(found) == 3


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


@pytest.mark.parametrize("correlation", [0., 1.3])
def test_blank_noise_and_depth_recover_untrimmed_gaussian_scatter(correlation):
    rng = np.random.default_rng(123)
    data = rng.normal(size=(800, 800))
    if correlation:
        data = gaussian_filter(data, correlation)
    image = noise_image(data)
    radius, count = 5., 6000
    rng = np.random.default_rng(42)
    positions = rng.uniform(15, 785, (count, 2))
    sums = aperture_photometry(data, CircularAperture(positions, r=radius))["aperture_sum"]
    expected = float(np.std(sums))
    cfg = Config({"photometry": {"aperture_scale": 1.}, "image": {"edge_margin_px": 0}})
    noise = ph.blank_aperture_noise(image, radius, cfg, n_apertures=count)
    limit = ph.empirical_limit(image, radius, 25., cfg, n_apertures=count)
    inferred = 10. ** ((25. - limit) / 2.5) / 3.
    assert noise == pytest.approx(expected, rel=0.07)
    assert inferred == pytest.approx(expected, rel=0.07)


def test_blank_noise_rejects_sources_and_partly_masked_apertures():
    data = np.random.default_rng(12).normal(size=(600, 600))
    baseline = noise_image(data.copy())
    image = noise_image(data.copy())
    yy, xx = np.mgrid[:600, :600]
    # Sparse bright positive and negative sources, as in a difference frame.
    for i, x in enumerate(range(60, 590, 100)):
        for y in range(60, 590, 100):
            image.data += (-1.) ** i * 1000 * np.exp(-((xx-x)**2 + (yy-y)**2) / 8.)
    image.mask[:, :180] = True
    image.data[:, :180] = np.nan
    cfg = Config({})
    expected = ph.blank_aperture_noise(baseline, 5., cfg, n_apertures=5000)
    actual = ph.blank_aperture_noise(image, 5., cfg, n_apertures=5000)
    assert actual == pytest.approx(expected, rel=0.08)


def test_blank_noise_returns_unknown_without_enough_complete_apertures():
    image = noise_image(np.ones((100, 100)))
    image.mask[:] = True
    assert np.isnan(ph.blank_aperture_noise(image, 5., Config({})))


@pytest.mark.parametrize("method,bad", [
    ("difference_psf", lc.FLAG_MASKED_PIXEL),
    ("difference_psf", lc.FLAG_OFF_IMAGE),
    ("difference_psf", lc.FLAG_BAD_PSF),
    ("difference_aperture", lc.FLAG_NONLINEAR_CORE),
])
@pytest.mark.parametrize("warning", [0, lc.FLAG_DIPOLE_RESIDUAL | lc.FLAG_BAD_SUBTRACTION])
def test_host_guard_never_restores_invalid_target_measurements(method, bad, warning):
    difference = lc.new_row(method=method, mag=19., flags=bad | warning, diff_chi2=1.)
    direct = lc.new_row(method="psf", mag=17., flags=0)
    assert lc.choose_primary([difference, direct]) is direct


def test_host_guard_requires_dipole_evidence():
    difference = lc.new_row(method="difference_psf", mag=19.,
                            flags=lc.FLAG_BAD_SUBTRACTION, diff_chi2=1.)
    direct = lc.new_row(method="psf", mag=17., flags=0)
    assert lc.choose_primary([difference, direct]) is direct
    difference["flags"] |= lc.FLAG_DIPOLE_RESIDUAL
    assert lc.choose_primary([difference, direct]) is difference


def test_host_guard_can_keep_a_valid_psf_with_reconstructed_core():
    difference = lc.new_row(method="difference_psf", mag=19., diff_chi2=1.,
                            flags=lc.FLAG_NONLINEAR_CORE | lc.FLAG_DIPOLE_RESIDUAL
                            | lc.FLAG_BAD_SUBTRACTION)
    direct = lc.new_row(method="psf", mag=17., flags=0)
    assert lc.choose_primary([difference, direct]) is difference


def test_empirical_limit_on_small_image_falls_back_without_crashing():
    image = noise_image(np.random.default_rng(7).normal(size=(20, 20)))
    cfg = Config({"image": {"edge_margin_px": 32}})
    assert ph.empirical_limit(image, 3., 25., cfg) == ph.limiting_magnitude(image, 3., 25., cfg)
