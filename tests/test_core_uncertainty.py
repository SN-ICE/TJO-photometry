"""Nonlinear-core uncertainty survives measurement, CSV and calibration transfer."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from tjophot import lightcurve as lc, pipeline, zeropoints
from tjophot.calibration import ZeroPoint
from tjophot.config import Config


def measured_row(method, core, calibrated=False, flux=10000.):
    zp = ZeroPoint("r", 24. if calibrated else np.nan, .01, 0., 0., .5,
                   20, .01, "test", "ab", ok=calibrated)
    common = dict(target="test", band="r", system="ab", image="missing.fits",
                  mjd=60003.5, zp=zp.zp, zp_err=zp.zp_err, cal_provider="test",
                  core_syst_mag=core, core_nonlinear_px=4, core_masked_px=4,
                  flags_base=lc.FLAG_NONLINEAR_CORE | (0 if calibrated else lc.FLAG_NO_ZEROPOINT))
    return pipeline._photometry_row(common, method, flux, 10., zp, 0.,
                                     SimpleNamespace(meta=SimpleNamespace(exptime=10.)),
                                     flux_err_sys=5.)


def campaign(row):
    rows = []
    for night in range(4):
        for band, zp in (("g", 25.), ("r", 24.)):
            if night == 3 and band == "r":
                continue
            rows.append(lc.new_row(target="test", image=f"{night}_{band}.fits", band=band,
                                    mjd=60000.5 + night, method="psf", zp=zp, zp_err=.01,
                                    zp_source="fit", cal_provider="test", flags=0))
    return lc.to_frame(rows + [row])


@pytest.mark.parametrize("method", ["psf", "difference_psf"])
@pytest.mark.parametrize("core", [0., .07])
def test_saved_core_uncertainty_overrides_new_config_after_csv(tmp_path, method, core):
    row = measured_row(method, core)
    frame = campaign(row)
    paths = lc.write(frame, tmp_path, "test")
    loaded = pd.read_csv(paths["all"])
    assert loaded.loc[loaded.image == "missing.fits", "core_syst_mag"].iloc[0] == core
    cfg = Config({"photometry": {"nonlinear_core_syst_mag": .2}})
    out, report = zeropoints.transfer(loaded, cfg)
    assert len(report["recovered"]) == 1
    recovered = out.loc[out.image == "missing.fits"].iloc[0]
    expected = np.sqrt((1.0857362 * np.hypot(1., .5) / 1000.)**2 + recovered.zp_err**2 + core**2)
    assert recovered.mag_err == pytest.approx(expected)
    assert recovered.core_syst_mag == core
    assert recovered.snr == pytest.approx(1000.)
    assert not recovered.is_limit
    again, report = zeropoints.transfer(out, cfg)
    assert not report["recovered"]
    pd.testing.assert_frame_equal(out, again)


@pytest.mark.parametrize("configured,expected", [(None, .05), (.08, .08)])
def test_legacy_core_uncertainty_is_estimated_with_warning(caplog, configured, expected):
    frame = campaign(measured_row("psf", .07)).drop(columns="core_syst_mag")
    cfg = Config({"photometry": {"nonlinear_core_syst_mag": configured}}) if configured is not None else None
    out, _ = zeropoints.transfer(frame, cfg)
    row = out.loc[out.image == "missing.fits"].iloc[0]
    assert row.core_syst_mag == expected
    assert row.mag_err >= expected
    assert "no stored nonlinear-core uncertainty" in caplog.text


def test_unflagged_legacy_rows_do_not_gain_core_uncertainty():
    row = measured_row("psf", 0.)
    row["flags"] = lc.FLAG_NO_ZEROPOINT
    frame = campaign(row).drop(columns="core_syst_mag")
    out, _ = zeropoints.transfer(frame)
    assert out.loc[out.image == "missing.fits", "core_syst_mag"].iloc[0] == 0.


def test_core_uncertainty_is_kept_even_for_nondetections():
    out, _ = zeropoints.transfer(campaign(measured_row("difference_psf", .07, flux=-10.)))
    row = out.loc[out.image == "missing.fits"].iloc[0]
    assert row.is_limit and np.isnan(row.mag_err)
    assert row.core_syst_mag == .07


def test_measured_calibrated_error_includes_the_saved_systematic():
    row = measured_row("psf", .07, calibrated=True)
    expected = np.sqrt((1.0857362 * np.hypot(10., 5.) / 10000.)**2 + .01**2 + .07**2)
    assert row["mag_err"] == pytest.approx(expected)
    assert row["core_syst_mag"] == .07
