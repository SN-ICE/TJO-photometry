# tjophot

[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

A supernova photometry pipeline for the Telescopi Joan Oro (TJO), built on
astropy and photutils. It takes a directory of reduced TJO frames and produces a
calibrated, multi-band light curve with per-epoch diagnostics.

Nothing in it is specific to one supernova: point it at any target directory and
it will work out the bands, the epochs, the seeing, the zero points and the
reference images by itself. It handles both TJO cameras, MEIA2 and MEIA3, and
picks the right one per frame.

`snname` below stands for the target's name -- its name on TNS, which is also
the name of its directory under `paths.data_root` and of its file in `targets/`.
Substitute your own.

---

## Quick start

```bash
cd ~/Desktop/TJO/pipeline
pip install -e .                         # -e matters: it keeps config/ and calib/ live

tjophot build-pattern snname --extra-dirs ~/Desktop/TJO/finalfiles   # once per instrument
./full-run.sh snname                  # clean, resolve, fetch, measure
```

`QUICKSTART.md` covers this in more detail, including how to add a new target.

Outputs land in `runs/snname/`:

| file | contents |
|---|---|
| `snname_lightcurve.csv` | one row per epoch: the adopted measurement |
| `snname_photometry_all.csv` | every method for every epoch, for cross-checking |
| `snname_lightcurve.png` | the light curve |
| `snname_diagnostics.png` | seeing, zero point and depth over the campaign |
| `qa/<image>.png` | per-epoch panel: stamps, curve of growth, zero-point fit, kernel |
| `templates/template_<band>.png` | reference and science frame on a common North-up grid |
| `run_summary.json`, `config_used.yaml`, `run.log` | exactly what was run, including which config files were merged (API keys redacted) |

Without installing, every command works as `python -m tjophot.cli ...` run from
inside the pipeline directory.

---

## What it does, step by step

**0. Detector pattern.** If a pattern has been measured for the band (see
*Detector patterns* below) it is removed before anything else, because fringes
left in the frame end up in the background model and in every sky annulus.

**1. Ingest.** Headers are read through an instrument description
(`config/instruments/tjo_meia2.yaml`), so nothing in the code hard-codes TJO
keywords. `FILTER = 'SDSS g'` becomes the canonical band `g`; `B`/`V` are
reported in Vega magnitudes, `ugriz` in AB. Saturated, dead and edge pixels are
masked, saturation is grown by two pixels to catch bleeding, and a 2D background
(`photutils.Background2D`) plus a per-pixel error map
(sky RMS squared + source counts / gain) are built.

**2. Seeing.** Bright, isolated, unsaturated stars are fitted with 2D Gaussians
and the median gives the FWHM, iterated once so the detection FWHM and the
measured FWHM agree. Every aperture and fitting box in the pipeline is expressed
in units of this FWHM, so the same configuration works in 1.5-arcsecond and
4-arcsecond seeing.

**2b. Astrometry.** The header carries an astrometry.net solution, which is
close but not tight. Detected stars are matched to the reference catalogue and
the WCS is re-fitted -- translation only, full TAN, or TAN plus SIP -- with
outlier rejection and two matching passes. The refit is thrown away if it moves
the frame further than `astrometry.max_shift_arcsec`, leaves a worse residual
than `astrometry.max_rms_arcsec`, or is simply worse than what it started from,
so the step can only help. On SN2026acd it takes the residuals from about
0.33 px to 0.18 px. This matters because every aperture, PSF fit and difference
stamp is placed through the WCS.

**3. Photometry.** Two independent measurements are always made:

- *Aperture*: a circular aperture of 1.2 x FWHM with a sigma-clipped annulus
  background, plus an aperture correction to 4 x FWHM measured from the field's
  own curve of growth.
- *PSF*: an empirical PSF (`EPSFBuilder`) built from stars that pass a
  point-source test, fitted at the SN position with a bounded recentring. PSF
  fluxes are put on the aperture-corrected system with a
  PSF-to-aperture offset measured on the same stars, so one zero point applies
  to both.

**4. Difference imaging.** A reference image is resampled onto the science
frame, PSF-matched with an Alard & Lupton kernel (a basis of Gaussians modulated
by polynomials, solved by linear least squares on stamps around bright stars,
with a low-order polynomial absorbing the differential background) and
subtracted. Aperture and PSF photometry are then repeated on the difference.
Implemented in numpy/scipy only -- no HOTPANTS or other external binary.

If the reference has worse seeing than the science frame, the pipeline detects
it and convolves the *science* frame instead (matching a blurry reference to a
sharp science image would be a deconvolution), then divides the measured flux by
the kernel sum to return to science-frame counts.

**5. Calibration.** Field stars are matched to a reference catalogue and

    m_cat = m_inst + ZP + C x (colour - pivot)

is fitted by weighted least squares with sigma clipping, where `m_inst` is the
aperture-corrected instrumental magnitude *per second*. Both `ZP` and the colour
term `C` are free by default; `calibration.color_term: fixed` uses values you
supply instead (see `tjophot colorterms`).

**6. Depth and errors.** The limiting magnitude is measured empirically, from
the scatter of hundreds of blank apertures, not from an analytic formula. The
same trick sets the noise floor on the SN measurement, which matters most on
difference images: resampling and PSF matching correlate neighbouring pixels, so
summing the per-pixel variances underestimates the aperture noise by a large
factor. Anything below 3 sigma is reported as an upper limit rather than a
magnitude.

---

## Detector patterns

TJO z-band frames carry a strong fringe pattern -- interference of night-sky
emission lines inside the thinned CCD -- of about 5 per cent of the sky
peak-to-peak with a wavelength near 100 px, plus a ~1 per cent illumination
residual present in every band. The second is a flat-field residual, not
vignetting: measured on stacked, source-masked frames the response is flat to
+-1 per cent right out to the corners, with only the outermost row and column
deviating. There is no shadowed border on MEIA2 to trim. Both patterns are
fixed in *detector* coordinates, so they can be measured once and removed from
every frame:

```bash
tjophot build-pattern snname --extra-dirs ~/Desktop/TJO/finalfiles
```

This median-combines source-masked, sky-normalised frames in detector
coordinates and writes `calib/pattern_<band>.fits` plus a PNG preview. TJO
dithers by about 50 px, which is enough for the median to reject stars, but
adding frames of other fields with `--extra-dirs` makes a cleaner pattern and
the command warns when everything comes from a single field.

The two components are corrected differently because they are physically
different things. The smooth component is a sensitivity residual -- it scales
source flux as well as sky -- so the frame is *divided* by it. The ripples are
sky emission that never went through the source's optical path, so a scaled copy
is *subtracted*, with the amplitude fitted on each frame: fringe strength
follows the sky emission lines, and moonlight raises the continuum without
adding fringes. Before subtracting, the pipeline checks on the frame itself
that the fringe term actually flattens it, and skips it otherwise -- under a
bright twilight gradient the fit can latch onto structure that is not the fringe
pattern.

**What it is worth.** On SN2026acd the z-band 3-sigma limiting magnitude
improves by **1.2 to 1.4 magnitudes** and the SN's own photometric errors drop
by roughly a third to a half (0.048 to 0.031 mag over the campaign; 0.037 to
0.021 on the six epochs re-measured against the current pattern). The gain is large because the fringe
amplitude, a few tens of ADU per pixel, is comparable to a faint source once
integrated over an aperture, and the sky annulus samples a different fringe
phase than the aperture does. Bands other than z gain little (r: +0.03 mag) and
are never made worse.

Patterns are tied to a readout geometry: a frame with different binning, a
different window offset or a different size is left uncorrected, with a warning,
rather than having the pattern stretched onto it.

## Reference catalogues

Configured by `calibration.catalogs`, tried in order until one produces a usable
fit for the band in question.

| provider | bands | notes |
|---|---|---|
| `gaia_astrometry` | -- (positions only) | Gaia DR3 positions for the WCS refinement. Separate from `gaia_synth` because an astrometric fit wants every well-measured position, not just the brighter subset that has synthetic photometry. |
| `gaia_synth` | U B V R I, u g r i z | **default.** Gaia DR3 synthetic photometry (GSPC) through the ESA Gaia TAP service. All-sky and internally consistent, so *B, V and griz come from one system* -- the reason it is first. Limited to G < 17.65. |
| `panstarrs` | g r i z y | PS1 DR2 mean objects via the MAST API. Deeper, but Dec > -30 and no B/V. |
| `apass` | B V g r i | AAVSO APASS DR9 via VizieR. The B/V fallback where GSPC is thin. |
| `local` | whatever the file has | a CSV/ECSV of your own with `ra, dec, <band>, e_<band>`. |

Results are cached under `cache/catalogs/`, so a field is downloaded once and
every later run is offline. `tjophot fetch-catalog` does that deliberately --
useful before travelling, or on a machine behind a firewall.

There are **no colour transformations between systems** anywhere in the
pipeline. B and V are calibrated against catalogues that genuinely measure B and
V; the residual bandpass mismatch between TJO's filters and the standard system
is absorbed by the fitted colour term. Nothing is converted with published
polynomial relations, because those add a systematic you cannot see.

---

## Reference images for subtraction

Set per band in `subtraction.reference_sources`, tried in order:

- **`file`** -- an explicit path in `subtraction.reference_files.<band>`. Use
  this whenever you have a genuine pre- or post-explosion image; it is always
  the best option.
- **`tjo_stack`** -- a median stack of the target's *own* TJO epochs, chosen by
  the MJD window in `subtraction.tjo_stack`. Same telescope, same filters, same
  pixel scale, no cross-instrument colour issues. It requires that window to be
  set explicitly, and it refuses to run otherwise, because a stack that still
  contains the supernova will quietly subtract signal from the light curve.
- **`panstarrs`** -- a PS1 DR2 stacked cutout, downloaded on first use and
  cached, linearised from its asinh (`BSOFTEN`/`BOFFSET`) storage. griz(y) and
  Dec > -30 only, so B and V cannot use it. This is the default for g, r, i and
  z, and needs no configuration: run `tjophot fetch-template <TARGET>` to fetch
  and inspect the cutouts on their own, or just let `tjophot run` do it.
  The request is sized from `subtraction.cutout_arcmin` (1.25x, capped by
  `subtraction.panstarrs.max_size_px`) -- about 2700 px and 29 MB per band at
  the default 9 arcmin stamp.

When no reference is available for a band, the pipeline says so and reports
direct PSF/aperture photometry for that band instead of failing.

Difference images are checked against their own noise model. The reduced
chi-square of the difference appears as `diff_chi2` in the output; above
`subtraction.max_chi2` (2.5) the epoch is flagged `bad_subtraction` and the
pipeline falls back to direct photometry when choosing what to report.

---

## Reading the output

`snname_lightcurve.csv` carries one row per epoch. `method` says how the
adopted value was measured, in order of preference:

`difference_psf` > `difference_aperture` > `psf` > `aperture`

`is_limit` marks non-detections, whose depth is in `limit_mag`. `flags_str`
spells out the quality bits:

| flag | meaning |
|---|---|
| `masked_pixel` | a masked pixel fell inside the aperture |
| `saturated` | saturated pixel inside the aperture |
| `off_image` | the aperture ran off the frame |
| `bad_subtraction` | the difference image failed its own noise model |
| `few_cal_stars` | the zero point rests on very few stars |
| `no_zeropoint` | the frame could not be calibrated at all |
| `poor_calibration` | the zero-point fit scatter exceeds `calibration.max_rms` |
| `in_template` | this epoch is part of the reference image, so its difference is null by construction |

Two further columns record the new preparation steps: `astrom_rms_arcsec` is the
residual of the refitted WCS against the reference catalogue (blank if the
refinement was skipped or rejected), and `pattern_applied` says whether a
detector pattern was removed from that frame.

`snname_photometry_all.csv` keeps every method for every epoch. Comparing
`aperture`, `psf` and `difference_*` for the same epoch is the fastest sanity
check there is: in an uncrowded field with a faint host they should agree to a
few hundredths, and where they do not, the host galaxy is the usual reason.

---

## Configuration

Four layers, merged in order, each overriding the last:

1. `config/default.yaml` -- every knob, documented inline
2. `targets/<NAME>.yaml` -- per-target settings (position, template window,
   reference files)
3. `--config my.yaml` -- overrides for this run
4. command-line flags (`--bands`, `--aperture-scale`, `--catalogs`, ...)

`config_used.yaml` in the output directory records the fully merged result, so
you can always see exactly what a run used.

A target file is where the SN position belongs:

```yaml
target:
  name: snname
  ra: "12:12:12.79"      # sexagesimal hours, or decimal degrees
  dec: "+13:13:18.2"
```

Without it the pipeline starts from the telescope pointing (`OBJRA`/`OBJDEC`)
and then refines it: the transient is centroided on the first few epochs and the
median position is adopted, provided the shift stays within
`target.refine_radius_arcsec`. On SN2026acd that correction is 1.1 arcsec --
three pixels, and enough to bias a 1.2 x FWHM aperture noticeably -- so it is
worth having, but it is still no substitute for the real coordinates. The
adopted position and the size of the refinement are recorded in
`run_summary.json`. Set `target.refine: always` to refine even when you have
supplied coordinates, or `off` to force photometry exactly where you said.

---

## Commands

```bash
tjophot list snname -v                 # epochs, bands, MJD coverage
tjophot fetch-catalog snname           # download and cache the catalogues
tjophot run snname                     # the full pipeline
tjophot run snname --bands g r         # a subset of bands
tjophot clean snname                   # remove this target's previous outputs
tjophot resolve snname --write         # look the position up on TNS
tjophot fetch-template snname          # fetch and inspect the reference images
tjophot build-pattern snname           # measure the detector patterns
tjophot run snname --no-subtraction    # direct photometry only
tjophot run snname --no-astrometry     # keep the header WCS
tjophot run snname --no-pattern        # ignore the stored detector patterns
tjophot run snname --refresh           # ignore cached catalogues/templates
tjophot run snname --skip 40 --append  # resume a long run, merging into the table
tjophot colorterms snname              # fit the instrument colour terms
```

`--append` merges a run into an existing table instead of replacing it, keyed on
(image, method), so you can redo one band or resume after an interruption
without losing the rest. `--skip N` starts from the Nth epoch.

`tjophot colorterms` runs the calibration over many frames and prints the
weighted mean colour term per band, already formatted as YAML to paste into a
config. Once you trust those numbers, switching `calibration.color_term` to
`fixed` makes every epoch's zero point a one-parameter fit, which is noticeably
more stable on frames with few calibration stars.

---

## Adding another instrument

Copy `config/instruments/tjo_meia2.yaml`, adjust the header keyword lists, the
`filters` map and the `systems` map, and set `instrument:` to the new file's
name. No Python changes are needed unless the data model itself differs (for
example multi-extension files).

---

## Layout

```
pipeline/
  config/default.yaml            all defaults, documented
  config/instruments/            header maps, one file per instrument
  targets/                       per-target settings
  tjophot/
    config.py        layered YAML configuration
    patterns.py      detector patterns: build, store, apply
    image.py         FITS ingest, pattern removal, masking, background, error map
    astrometry.py    WCS refinement against a reference catalogue
    detection.py     source detection, seeing, point-source selection
    photometry.py    apertures, curve of growth, EPSF, forced PSF, depth
    catalogs.py      reference catalogue providers and disk cache
    calibration.py   matching and the zero-point fit
    templates.py     reference image acquisition (file / TJO stack / PS1)
    subtraction.py   Alard & Lupton PSF matching and the difference image
    lightcurve.py    table assembly, flags, CSV output
    plots.py         QA figures
    pipeline.py      orchestration
    cli.py           command line
  calib/             measured detector patterns, one FITS + PNG per band
  tests/             offline validation (see below)
```

---

## Validation

`tests/make_synthetic_catalog.py` measures the field on one epoch and writes a
reference catalogue with zero points and colour terms of your choosing. Running
the pipeline against it must recover exactly those numbers, which exercises the
matching, the fit, the clipping and the whole output path without any network
access:

```bash
python tests/make_synthetic_catalog.py snname \
    --data-root ~/Desktop/TJO/finalfiles --mjd 61121
python -m tjophot.cli run snname --data-root ~/Desktop/TJO/finalfiles \
    --config config/test_local.yaml --work-root runs_test
```

On SN2026acd this recovers the injected zero points to 0.01-0.02 mag and the
injected colour terms to a few thousandths.

---

## Notes on the TJO data

Things worth knowing, measured on the SN2026acd campaign:

- Transparency varies a great deal between nights. The same field stars are
  1.36 mag fainter in g on MJD 61097 than on MJD 61121, and the offset is
  wavelength-dependent (1.36 in g, 1.16 in r, 1.07 in i) -- thin cloud or high
  airmass, not a data problem, but it is why every frame gets its own zero point
  and why a run-wide zero point would be badly wrong.
- The field is sparse: 20-70 detections per frame at 5 sigma. Frames taken in
  poor seeing can drop below the five stars needed for an empirical PSF, and the
  pipeline then reports aperture photometry for that epoch.
- The z band fringes at about +-3 per cent of the sky; g, r, i and V show the
  same fringe family five times weaker, and every band sits in a ~1 per cent
  vignetting bowl. `build-pattern` handles all of it.
- Some frames were taken in bright twilight or moonlight (sky up to 22000 ADU
  against a normal 250-400). The measured limiting magnitude tracks this
  faithfully, so cutting on `limit_mag` is an effective quality filter.

---

## Licence and citation

BSD 3-Clause, see `LICENSE`. If it contributes to a publication, `CITATION.cff`
has the metadata; please also cite astropy and photutils, which do the heavy
lifting underneath.
