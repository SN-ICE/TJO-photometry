# What has been tested, and how

All of the following was run against the real SN2026acd data set (155 frames,
6 bands, MJD 61097-61213) on a machine with no access to astronomical catalogue
servers. Everything that does not need the network was therefore exercised end
to end; the two network paths are marked below.

## Calibration: closed-loop test

`tests/make_synthetic_catalog.py` measures the field on one epoch, adds a chosen
zero point and colour term per band, and writes the result as a `local`
catalogue. Calibrating that same epoch against that catalogue must return the
injected numbers.

| band | injected ZP | recovered | injected C | recovered |
|---|---|---|---|---|
| B | 21.500 | 21.522 | 0.06 | 0.089 |
| g | 22.000 | 21.992 | 0.04 | 0.044 |
| i | 21.600 | 21.600 | 0.01 | 0.001 |
| r | 22.000 | 21.992 (median over frames) | -0.02 | 0.004 |

Fit residuals on that epoch are 0.003-0.010 mag. (The B colour term is fitted
against B-V while the test injects against g-r, which is why it differs most.)

Across all 103 epochs of the four-band demonstration run the median recovered
zero points were B 21.55, V 21.82, g 21.98, r 22.00 against injected values of
21.50, 21.80, 22.00, 22.00 -- the residual spread being real night-to-night
extinction, which is exactly what a per-frame zero point is for.

## Photometry

- **EPSF fidelity.** The empirical PSF's width, measured on the real pixel grid,
  matches the frame's measured seeing: 4.77 vs 4.77 px in g, 5.23 vs 5.31 in r,
  8.99 vs 9.05 in i.
- **PSF vs aperture.** On field stars, PSF flux over aperture(4 x FWHM) flux is
  0.98-1.00 with 1-5 per cent scatter, and the residual offset is removed by the
  measured PSF-to-aperture correction.
- **Method agreement on the SN.** Over the r-band campaign,
  `difference_aperture - difference_psf` has a median of +0.005 mag with a MAD
  of 0.030 -- two independent measurements on the same difference image agreeing
  at the 3 per cent level.

Three real bugs were found and fixed by these checks: photutils' `quartic`
smoothing kernel collapses the EPSF to near a delta function at oversampling 2
(now disabled by default, with a width sanity check that rejects a bad build),
and extended sources passing the brightness and isolation cuts were broadening
the PSF (now removed by a per-candidate width test). Separately, the pipeline's own
logging setup created a plain logger named ``astroquery`` before astroquery had
installed its Logger subclass, which made every astroquery import fail with
``'Logger' object has no attribute '_set_defaults'`` -- that is, every catalogue
query would have failed on the first machine with a working network connection.
Noisy libraries are now silenced with a log filter, which never creates a logger.

## Difference imaging

Reduced chi-square of the difference against its own noise model is 0.93-1.05
across the r-band campaign -- the subtraction leaves nothing but photon noise
behind. Frames where it does not (bad seeing, bright twilight) come out at
chi2 > 2.5 and are flagged `bad_subtraction`.

The error model was corrected during this work: resampling and PSF matching
correlate neighbouring pixels, so summing per-pixel variances underestimated the
aperture noise by a factor of about five. Errors are now floored by the scatter
of blank apertures measured on the difference image itself, after which the
difference-image and direct-aperture signal-to-noise agree.

## Target position

The telescope pointing (`OBJRA`/`OBJDEC`) sits 1.1 arcsec -- three pixels --
from the actual supernova in these frames, enough to bias a 1.2 x FWHM aperture.
The pipeline now centroids the transient on several epochs and adopts the median
position, and records the size of that correction in `run_summary.json`.

## Detector pattern correction

Built from 46 z-band frames spanning four fields and 16 months, and tested on
frames held out of the build:

| held-out frame sky | large-scale sky scatter before | after | flatter by |
|---|---|---|---|
| 354-609 ADU (6 frames) | 1.5-1.7% | 0.38-0.53% | 3.0-4.2x |
| 1280, 2555 ADU (moonlit) | 3.9%, 1.8% | unchanged | 1.0x |

The moonlit frames are not a failure: their large-scale structure is a moon
gradient, which is not a fixed detector pattern, and the guard described below
correctly declines to subtract there.

**Worth, in magnitudes.** Six SN2026acd z-band epochs, same frames, three
treatments:

| pattern | median mag_err | 3-sigma limit |
|---|---|---|
| none | 0.037 | 18.13 |
| 8 frames, a different field | 0.021 | 19.39 |
| 46 frames, four fields | 0.021 | 19.49 |

The correction is worth 1.4 magnitudes of depth and halves the error on the
supernova. Which frames it is built from is worth another 0.1 mag -- so having
a pattern matters far more than having a perfect one.

**How long a pattern lasts.** The amplitude is fitted per frame, so only the
shape has to hold. Patterns built 2.5 months apart from the same field correlate
at r = 0.969, against a split-half noise floor of 0.928-0.988 measured *within*
a single window: the drift is smaller than the noise of measuring it. A pattern
spanning 2026 corrects frames from 10 months earlier about as well as one
field's pattern corrects another field's frames. Nothing in 16 months of data
shows the shape moving, so rebuild after an instrument intervention -- a
warm-up, a re-coating, a filter change -- rather than on a schedule.

**Build from several fields.** A pattern measured on one field absorbs that
field's own sky structure and subtracts it back out of every frame of that
field. On held-out SN2026fvx frames a same-field pattern scored 0.194 in
residual large-scale scatter against 0.276 for a four-field pattern -- it looks
better while removing real signal. `build-pattern` warns when all its input
frames are of one field; the warning is worth heeding.

End to end on eight z-band epochs, with and without the correction:

| | 3-sigma limit | SN mag error | calibration rms |
|---|---|---|---|
| without pattern | 17.63 | 0.048 | 0.076 |
| with pattern | 19.05 | 0.031 | 0.064 |

**+1.2 magnitudes of depth in z**, and the epochs that gained most were the ones
that were worst (0.118 -> 0.037 mag, 0.096 -> 0.043 mag). The measured
magnitudes themselves agree to a few millimagnitudes, as they must for a bright
source. In r the same machinery gains +0.03 mag and never loses more than 0.05.

Bright-star repeatability in z improves by 7 per cent (0.0159 -> 0.0147 mag) and
its systematic excess over the formal errors by 11 per cent. The modest figure
for bright stars, against 1.2 magnitudes at the faint end, is the expected
signature of an *additive* sky artefact.

One guard was added after testing: an early version subtracted fringes from a
bright-twilight r frame (sky 22600 ADU) and made it half a magnitude *shallower*,
because the fit had latched onto the twilight gradient. The correction now
measures, on each frame's own block grid, whether the fringe term actually
flattens it, and skips it when it does not.

## WCS refinement

Tested against a consensus catalogue built from the median positions of stars
detected across the r-band campaign, so the numbers measure internal
consistency rather than absolute accuracy:

| | residual |
|---|---|
| header (astrometry.net) solution | 0.121" = 0.33 px |
| after a linear refit | 0.064" = 0.18 px |

Nine of nine frames improved; the fitted solutions moved the frame by 0.06-0.13
arcsec. Against the reference catalogue used in the demonstration run the
residuals go from 0.20-0.53" to 0.11-0.23". The absolute improvement against
Gaia should be larger, since the header solutions carry systematic errors that
an internal comparison cannot see.

## Not yet exercised against a live server

- `gaia_synth`, `panstarrs` and `apass` catalogue queries. The ADQL, the MAST
  query and the VizieR query are written against the documented schemas and the
  column mapping and caching around them are tested, but no response has been
  parsed from the real services. **Run `tjophot fetch-catalog SN2026acd` first**
  and check the reported source counts and bands before trusting a science run.
- The `gaia_astrometry` query used by the WCS refinement. In the demonstration
  run the refinement fell back to the local catalogue, which exercises the same
  code path.
## Pan-STARRS templates

Now exercised against real cutouts (griz, 2700 x 2700 px, 0.25"/px, taken
2012-2013, downloaded on the user's machine):

- Difference images against the PS1 g stack give a reduced chi-square of
  0.92-1.02 across eight epochs -- the subtraction is consistent with its own
  noise model, with a reference from a different telescope, a different pixel
  scale and a 14-year baseline.
- The cutouts arrive already linear: `fitscut` output carries no
  `BSOFTEN`/`BOFFSET`, so the asinh conversion is a no-op here. It is kept for
  the case where a raw stack is supplied instead, and now refuses to run twice.

Two problems this turned up:

- **Skycell edges.** PS1 stacks are tiled and the cutout service cuts from one
  tile only. This field sits 3 arcmin from a tile boundary, so a single-tile
  cutout left 22 per cent of the 9-arcmin stamp blank. The neighbouring
  skycells are now probed and mosaicked in; the combiner is exact where grids
  match and reprojects where they do not.
- **The convolution direction was being chosen by noise.** `_stamp_fwhm` took
  the science frame's FWHM as both its detection kernel and its fitting-box
  size, so on a deep reference the fit wandered onto blends and the "measured"
  reference width simply tracked the science seeing (ratio 0.83-1.17 across
  frames whose seeing ranged over a factor of two). Two of eight epochs were
  therefore convolving the sharp science frame down to a blurrier "reference"
  that was not actually blurrier. Measured properly, the PS1 stamp gives a
  constant 3.79-3.90 px regardless of the science frame, and all eight epochs
  now correctly convolve the reference.

## The gain model

The header's `GAIN = 1.04` e-/ADU with `RDNOISE = 9.0` e- does not describe
the frames TJO delivers. That much has held up; the first account of *why*,
written here earlier, did not, and is corrected below.

**What was measured.** Pixel variance (from second differences, which remove
gradients and fringe-scale structure) against sky level, on 42 MEIA2 frames
in six filters spanning a lunar cycle, so the sky runs from 135 ADU on a dark
night to 6886 ADU under the full moon. One line fits all of them to a few per
cent:

    sigma^2 = S / 1.95  -  21 ADU^2  +  (0.0041 S)^2

- **Gain 1.95 +- 0.03 e-/ADU**, and the same in every filter fitted on its
  own: B 1.95, V 1.96, g 1.94, r 1.95. A camera property, not a sky effect.
- **Offset -21 +- 2 ADU^2.** A clean bias-subtracted frame gives +rn^2/g^2,
  never a negative number. Negative means the "sky" carries an additive
  pedestal that makes no shot noise: 41 ADU with no read noise, 83 ADU with
  the header's 9 e-. The sky alone cannot split the two, and nothing in the
  pipeline needs it split -- local backgrounds remove the pedestal from every
  flux, and the noise prediction depends only on the line.
- **0.41 per cent pixel-to-pixel flat noise**, visible only above ~4000 ADU.
  Over a stellar footprint it averages down to ~0.05 per cent of the flux, so
  it is left out of the error model.

**Cross-check.** Two B frames taken 3.4 minutes apart give the temporal noise
at the sky with the flat cancelled: 39-44 ADU^2 at 134-139 ADU. Shot noise
alone at the header gain would be 133 ADU^2, so 1.04 cannot describe a single
frame whatever the read noise. The same pair shows 4-16 ADU^2 of fixed-pattern
noise on top, varying from night to night -- calibration masters, most likely.
It is real noise in one frame and the measured `bkg_rms` already carries it.

**What was wrong before.** The first version solved each frame's own sky for
the gain: sigma^2 = sky/g + (rn/g)^2 at a single level. At one level a gain
and a pedestal are degenerate, so with the pedestal ignored the answer came out
2.9-3.6 e-/ADU on dark nights and ~1.9 under the moon, falling with sky
brightness -- which is also why it fell from B to z. The "average of about
three sub-exposures" reading drawn from that was an artefact: averaging scales
the slope and the intercept together and can never make the intercept
negative. The full-well check now agrees with no averaging at all: stars clip
near 64000 ADU, which is ~125000 e- at 1.95.

**What the code does now.** The per-frame solve is gone. The instrument file
carries the measured `noise_model` (gain, offset_adu2), which wins over the
header. A negative offset is represented as a pedestal with zero read noise,
which reproduces the same line at every level: the error map and
`combine_errors` use the gain; L.A.Cosmic is handed the frame with the
pedestal taken off, so its threshold tracks the true noise at every level;
and the sky-noise check predicts with the pedestal removed. On the 42 frames
that check reads 0.997 median (0.94-1.17), against ~0.47 with the header.

**Where it matters.** Only the source's own shot noise. The sky term of every
error bar comes from the measured `bkg_rms`, so faint-source errors, the
difference chi-square and the limiting magnitudes were never affected. The
header model overstated a bright source's shot-noise variance by 1.9x; the
per-frame solve understated it by up to 1.75x on dark nights and less under
the moon. Neither moves a magnitude, and on a bright target the zero-point
term dominates the error regardless.

**Not measured:** MEIA3. Its instrument file still uses the header values,
and the sky-noise check will say so if they are wrong. Measuring it needs
frames across a range of sky levels -- moonlit nights give the leverage.
