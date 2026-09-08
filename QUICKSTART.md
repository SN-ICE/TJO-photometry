# Running tjophot

## One-time setup

```bash
cd ~/Desktop/TJO/pipeline
pip install -e .
```

`-e` matters: it installs in place, so the code keeps reading `config/`,
`targets/` and `calib/` from this directory. That also means editing a config
file takes effect immediately, with no reinstall.

If you would rather not install anything, every command below works as
`python -m tjophot.cli ...` *run from inside `~/Desktop/TJO/pipeline`*.

Python 3.9 or newer. `pip install -e .` pulls in numpy, scipy, astropy,
photutils, astroquery, reproject, matplotlib, pandas and PyYAML. If you keep a
conda environment for this kind of work, activate it first.

If pip warns that the `tjophot` script is not on your PATH, either add
`~/.local/bin` (or your environment's `bin`) to PATH, or just use
`python -m tjophot.cli` instead — they are the same thing.

## Running it

Throughout this file `snname` stands for the target's name -- the name it has
on TNS, which is also the name of its directory under `paths.data_root` and of
its file in `targets/`. Substitute your own everywhere it appears.

From a clean slate, for one target:

```bash
./full-run.sh snname
```

Extra flags are split sensibly: options every subcommand understands
(`--config`, `--data-root`, `--work-root`, `--bands`, `--log-level`) are
forwarded to *all* the steps, and anything else goes to the final `tjophot run`.
Each step echoes the command it is about to run, so you can see exactly what
happened. It is a thin wrapper over five commands you can equally run by hand:

```bash
# 1. clear out anything left from a previous run of this target
tjophot clean snname

# 2. look the position up on TNS and record it in targets/snname.yaml
tjophot resolve snname --write

# 3. cache the reference catalogues for the field          (internet, once)
tjophot fetch-catalog snname

# 4. fetch the reference images for difference imaging     (internet, once)
tjophot fetch-template snname

# 5. measure the light curve
tjophot run snname
```

Steps 2-4 are cached, so they cost nothing on a re-run. Step 1 only removes
`runs/snname/` -- the downloads and the detector patterns survive, which is
what you want when you are iterating on the photometry.

`tjophot run --clean` folds step 1 into step 5, so the shortest honest
re-run once the field is set up is:

```bash
tjophot run snname --clean
```

**Detector patterns** are a property of the camera, not of the target, so they
sit outside this recipe. Build them once and forget them:

```bash
tjophot build-pattern snname --extra-dirs ~/Desktop/TJO/finalfiles
```

They live in `calib/` and are applied automatically. Rebuild only when a
readout configuration appears that has none -- see "More than one camera"
below, and `tjophot list` will tell you.

**Starting genuinely from scratch**, including re-downloading everything:

```bash
tjophot clean snname --all      # add -n first to see what would go
```

`--all` also deletes this field's cached catalogues and reference images;
`--templates` and `--catalogs` pick one or the other. `-n` / `--dry-run` lists
what would be removed without touching anything, and the command refuses to
delete outside `runs/` and `cache/`.

A full run over about 150 frames takes roughly 20 minutes.

## Pan-STARRS templates

Difference imaging needs a reference image of the field without the supernova
in it. For g, r, i and z at Dec > -30 that comes from Pan-STARRS DR2, and it is
already the default -- `tjophot run` will download a cutout the first time it
needs one. `tjophot fetch-template` just does that step on its own so you can
check it before committing to a run:

```
  g: panstarrs  2700x2700 px, 0.250"/px, 100% valid, median 0.3, range -3.1 to 892.4
     cached at  cache/templates/ps1_183.05315+13.22175_g-sdss_2700.fits
     preview    runs/snname/templates/template_g.png
```

**Look at the previews.** Each one shows the full cutout, then the template and
one of your own science frames cut around the *same* sky coordinate at the same
angular size, with the target marked in both. The supernova should be plainly
inside the marker on the science panel and absent from the template panel. If
it is in the template, that "template" contains the thing you are measuring and
the photometry will be biased.

Both stamps are resampled onto one common North-up, East-left grid centred on
the target, with a compass drawn in the corner, so they can be compared with
each other and with any finder chart directly. This matters more than it
sounds: TJO's WCS has both CD diagonal terms negative, so raw TJO frames are
North-*down*, and drawing them as stored makes the two panels mirror images of
one another.

Each cutout is 2700 x 2700 PS1 pixels (about 29 MB) and takes a minute or so.
The size follows `subtraction.cutout_arcmin`; shrink that, or set
`subtraction.panstarrs.max_size_px`, if the service is slow.

PS1 stacks are tiled into skycells and the cutout service only ever cuts from
one tile, so a field near a tile edge comes back with a blank strip -- the
SN2026acd field is 3 arcmin from one, and a single-tile cutout loses 22 per cent
of the stamp. The pipeline probes the neighbouring skycells and fills the gap
(`subtraction.panstarrs.mosaic`, on by default). The `% valid` figure printed by
`fetch-template` is what to watch.

B and V have no Pan-STARRS counterpart of their own, and no CCD survey covers
them all-sky. Three ways round that, best first:

1. **Take your own.** A set of deep TJO B and V images of the field, once the
   supernova is gone, is the ideal template: same telescope, same filters, same
   pixel scale, no colour term at all. Point `reference_files` at them.
2. **Download a public-archive template.** Any survey published as a HiPS can
   be rendered to a FITS cutout by the CDS `hips2fits` service, which is
   all-sky and returns exactly the grid asked for -- so, unlike a tiled cutout
   service, never with a blank strip through it. Configured per band:

   ```yaml
   subtraction:
     reference_sources:
       B: [file, hips, panstarrs, tjo_stack]
     hips:
       B: {survey: dss2-blue}      # or a full id: CDS/P/DSS2/blue
       V: {survey: dss2-red}
   ```

   Shorthands: `dss2-blue`, `dss2-red`, `dss2-nir`, `sdss9-u/g/r/i/z`,
   `ps1-g/r/i/z/y`. Anything else takes a full HiPS id; browse them at
   <https://aladin.cds.unistra.fr/hips/list>. `pixel_scale` is optional and
   defaults to the survey's own resolution.

   DSS2 blue and red are the closest wide-field bandpass matches to B and V,
   and are the default. They are photographic, though: their response is not
   proportional to flux, and a single kernel scaling cannot match a non-linear
   plate to a linear CCD frame across a bright galaxy. Sometimes a deeper,
   linear survey in the *wrong* filter wins anyway --
   `subtraction.reference_band: {B: g, V: r}` borrows the PS1 stacks, and those
   epochs come out flagged `cross_band_template`.

   Which is better is a question about your field, not an argument to have:

   ```bash
   tjophot fetch-template snname --compare --bands B V
   ```

   builds a difference against each candidate on one of your own frames and
   prints the reduced chi-square and residual scatter for each. Lower is
   better; near 1 means the difference is consistent with photon noise. You can
   add more candidates to the comparison by listing them:
   `--compare sdss9-g dss2-red`.
3. **No template.** The pipeline reports plain PSF photometry and says so. Fine
   while the supernova is much brighter than the host, useless once it is not.

For option 1, point the target file at the images:

```yaml
subtraction:
  reference_files:
    B: ~/Desktop/TJO/refs/snname_ref_B.fits
    V: ~/Desktop/TJO/refs/snname_ref_V.fits
```

or let the pipeline stack late-time epochs itself with an MJD window under
`subtraction.tjo_stack`. Epochs that end up inside their own template are
flagged `in_template`.

A same-filter template that still contains the supernova is *worse* than a
clean cross-band one: it subtracts the signal you are trying to measure, and
nothing downstream can recover it. That is why PS1 is listed ahead of
`tjo_stack` for B and V in the SN2026acd target file.

## More than one camera

A TJO archive is not one instrument. Frames in this one come from three
different readout configurations:

| | array | scale | field | gain | read noise |
|---|---|---|---|---|---|
| MEIA2, windowed | 1048 x 1048 | 0.361 "/px | 6.3' | 1.04 | 9.0 |
| MEIA2, full frame | 2048 x 2048 | 0.361 "/px | 12.3' | 1.04 | 9.0 |
| MEIA3 | 3108 x 3096 | 0.402 "/px | 20.8' | 1.94 | 8.3 |

The first two are the same chip read out through different windows; the third
is a different camera. A single campaign can straddle either change -- SN2025few
switches from the 1048 window to full frame at MJD 61013, mid-light-curve.

The pipeline works this out per frame, from `INSTRUME` and the array size, and
everything that depends on it follows automatically:

* **gain, read noise, plate scale and saturation** come from the frame's own
  instrument file, chosen through the `instruments:` map in `config/default.yaml`;
* **detector patterns** are stored per readout, as
  `calib/pattern_<band>_<detector>.fits`, because a fringe pattern is fixed in
  detector coordinates and is meaningless on an array of another size. TJO
  records `OFFSET1 = OFFSET2 = 0` even for a windowed readout, so where the
  window sits on the chip is not recoverable from the header -- which is why a
  full-frame pattern is not cropped onto a subframe, but measured separately;
* **catalogue search radii and stamp sizes** default to `auto` and are derived
  from the frame's own field of view, snapped to whole arcmin so the download
  cache stays small. A fixed 9 arcmin radius covers a 12.3' frame but wastes
  most of its stars on a 6.3' one and misses half of a 20.8' one.

`tjophot list` shows the split, and says which patterns are still missing:

```
2 detector configurations:
  detector                    N    MJD first     MJD last  bands
  MEIA2 1048x1048            81   60786.9465   61011.1862  B,V,g,i,r,z
  MEIA2 2048x2048            40   61017.1079   61134.0660  B,V,g,i,r
```

Every epoch records the camera it came from in the `detector` column of the
output tables, so a step in the light curve can be checked against the swap.

Building patterns for a new readout is the same command; it groups the frames
itself and writes one pattern per band per configuration. Feed it the whole
archive so each pattern is measured over several fields:

```bash
tjophot build-pattern snname --extra-dirs ~/Desktop/TJO/finalfiles
tjophot build-pattern snname --extra-dirs ~/Desktop/TJO/finalfiles \
    --detectors meia2_1048x1048        # or just one of them
```

Adding a camera the pipeline has not seen means two edits: a new file in
`config/instruments/`, and a line in the `instruments:` map naming it. Without
them the run still works, but it warns and falls back to the default
instrument's numbers, which is exactly the silent-wrong-answer case worth
avoiding.

## What you get

Everything lands in `runs/snname/`:

| file | what it is |
|---|---|
| `snname_lightcurve.csv` | the light curve: one row per epoch |
| `snname_photometry_all.csv` | every method for every epoch |
| `snname_lightcurve.png` | the light curve, plotted |
| `snname_diagnostics.png` | seeing, zero point and depth over the campaign |
| `qa/<image>.png` | per-epoch panel: stamps, curve of growth, zero-point fit |
| `run.log`, `run_summary.json`, `config_used.yaml` | what was actually run |

Start with `run.log` and `snname_diagnostics.png` when something looks wrong.

## The SN position

The pipeline looks the target name up and writes the answer into the target
file, so this only has to happen once:

```bash
tjophot resolve snname --write
```

It queries the Transient Name Server first, then falls back to SIMBAD. A run
does the same automatically whenever `ra`/`dec` are missing from the target
file (`target.resolve: auto`); set `always` to re-check every run, or `off` to
leave it alone. Only the `ra` and `dec` lines are rewritten -- every comment and
every other setting in the file is preserved -- and a comment recording where
the position came from and when is written above them.

TNS needs a bot key. Create one at <https://www.wis-tns.org/bots>, then:

```bash
cp config/tns.yaml.example config/tns.yaml
# edit config/tns.yaml with the bot's api key, numeric id and name
```

`config/tns.yaml` is loaded automatically on every command -- no `--config`
needed -- and is gitignored. `TNS_API_KEY`, `TNS_BOT_ID` and `TNS_BOT_NAME`
work as well, with the file taking precedence. Without credentials it falls
back to SIMBAD, which needs none but is coarser and lags TNS for very recent
objects.

Every command prints which configuration files it merged, in order:

```
INFO tjophot.config  configuration: default.yaml + tns.yaml + targets/snname.yaml
```

If a setting is not taking effect, that line is the first thing to check. The
same list is recorded in `run_summary.json`, and `config_used.yaml` holds the
fully merged result with the API key redacted.

You can of course still write the position by hand:

```yaml
target:
  ra: "12:12:12.79"
  dec: "+13:13:18.2"
```

Without it the pipeline starts from the telescope pointing and centroids its way
to the transient, which works but is not forced photometry at a known position.
The TNS coordinates are what you want here.

## Useful variations

```bash
tjophot run snname --bands g r i          # a subset of bands
tjophot run snname --max-images 6         # a quick sanity run
tjophot run snname --no-subtraction       # direct photometry only
tjophot run snname --no-pattern           # skip the detector-pattern removal
tjophot run snname --no-astrometry        # keep the header WCS
tjophot run snname --refresh              # ignore cached catalogues/templates
tjophot run snname --bands z --append     # redo one band, merge into the table
tjophot run snname --skip 80 --append     # resume after an interruption
tjophot colorterms snname                 # fit the instrument colour terms
```

Every flag has a config equivalent; `--config my.yaml` merges a file of your own
over the defaults and the target file.

## A new supernova, start to finish

Nine steps. Steps 1-3 are the only ones that need thought; the rest is
mechanical, and 4-8 are what `full-run.sh` does for you.

### 1. Put the frames where the pipeline looks

```bash
mkdir -p ~/Desktop/TJO/finalfiles/snname
# copy or link the reduced frames in
tjophot list snname
```

`list` is the check that matters: it prints one row per band with the frame
count and the MJD range. If it says `0 images`, the frames are in the wrong
place or `paths.data_root` points somewhere else. Files may sit either in
`<data_root>/snname/` or directly in `<data_root>/`; the subdirectory is
tried first and wins.

Nothing else about the instrument needs setting up. The detector patterns in
`calib/` are a property of the camera, not the target, and are picked up
automatically -- rebuild them only if the CCD, filters or binning change.

### 2. Write the target file

```bash
cp targets/TEMPLATE.yaml targets/snname.yaml
```

Edit `target.name`. Leave `ra`/`dec` out unless you already know them -- the
next step fills them in. The rest of the template is the subtraction setup,
which the next section covers; the defaults are right for a northern field.

### 3. Get the position

```bash
tjophot resolve snname --write
```

Queries TNS, falls back to SIMBAD, and rewrites only the `ra`/`dec` lines of
the target file, leaving every comment intact and recording where the position
came from. This is worth doing explicitly rather than letting the run do it,
because it is the one input that silently ruins everything downstream: forced
photometry at the wrong position gives you a light curve of empty sky. Check
the printed coordinates against the TNS page before going on.

TNS needs a bot key (see "The SN position" above). Without one you get SIMBAD,
which is coarser and lags TNS badly for objects discovered in the last months.

### 4. Cache the reference catalogues

```bash
tjophot fetch-catalog snname
```

Read the star counts it prints, per band. The pipeline needs at least a handful
of usable calibrators per frame; a field with two APASS B stars will produce B
magnitudes that are technically calibrated and practically worthless. Gaia
synthetic photometry covers griz well almost everywhere. B and V come from
APASS, which thins out at high galactic latitude -- if the counts look thin
there, expect `few_cal_stars` and large `zp_err` in those bands, and plan to
lean on `recalibrate`.

### 5. Get the templates and *look at them*

```bash
tjophot fetch-template snname
```

Then open `runs/snname/templates/template_*.png`. Each preview puts the
template and one of your own science frames side by side, on a common
North-up, East-left grid, with the target marked in both. Two things to
confirm, per band:

- the supernova is plainly inside the marker on the science panel;
- it is **absent** from the template panel.

A template containing the supernova subtracts away the signal you are trying to
measure, and no amount of later processing recovers it. This is the single
worst failure mode in the whole pipeline and the preview is the only thing that
catches it.

Also watch the `% valid` figure. PS1 cutouts come from tiled skycells and a
field near a tile edge can come back with a blank strip; the pipeline mosaics
the neighbours automatically, but if `% valid` is still well below 100 the
stamp has a hole in it.

For B and V there is no Pan-STARRS counterpart, so decide between the
photographic DSS2 template and a cross-band PS1 one by measuring it on your own
field rather than arguing about it:

```bash
tjophot fetch-template snname --compare --bands B V
```

It differences each candidate against one of your frames and prints the reduced
chi-square. Lower wins; near 1 means the residual is consistent with photon
noise. Set the winner in the target file.

If the field is south of Dec -30, Pan-STARRS does not reach it at all: drop
`panstarrs` from `reference_sources` and give each band a HiPS survey that does
cover it (SkyMapper, DECaLS).

### 6. A short trial run

```bash
tjophot run snname --max-images 6
```

Six frames, a couple of minutes. Then look at two or three panels in
`runs/snname/qa/`. Each one shows the science stamp, the template, the
difference, the curve of growth and the zero-point fit. What you are checking:
the difference stamp is flat apart from the supernova, the curve of growth
flattens out, and the zero-point fit is a tight line rather than a scatter
plot. If those look right on six frames they will look right on two hundred.

### 7. The full run

```bash
./full-run.sh snname
```

This repeats steps 3-6 in order and then measures everything. Steps 3-5 are
cached, so the repetition is free. Roughly 20 minutes for 150 frames.

To re-run later without re-downloading anything, `tjophot run snname
--clean` is the short form.

### 8. Read the output

Everything is in `runs/snname/`. In the order worth looking at them:

1. `snname_diagnostics.png` -- seeing, zero point and depth across the
   campaign. Bad nights are obvious here and nowhere else.
2. `snname_lightcurve.png` -- filled circles are directly calibrated,
   open circles are zero points transferred from another band that night.
3. `snname_lightcurve.csv` -- one row per epoch. Always read `flags_str`
   alongside `mag`; `zp_err` is the honest uncertainty and the flags say how
   the number was arrived at.
4. `run.log` when any of the above surprises you.

The run finishes by transferring zero points to any frame that could not fit
one, so `recalibrate` is normally already done. Run it by hand only after you
change a calibration setting.

### 9. Later: replace the borrowed templates

Once the supernova has faded, take your own deep B and V frames of the field.
Same telescope, same filters, same pixel scale, no colour term -- strictly
better than anything downloaded. Point the target file at them and re-run:

```yaml
subtraction:
  reference_files:
    B: ~/Desktop/TJO/refs/snname_ref_B.fits
    V: ~/Desktop/TJO/refs/snname_ref_V.fits
```

Or let the pipeline stack late-time epochs itself with an MJD window under
`subtraction.tjo_stack` -- but only once the supernova is genuinely gone, not
merely faint. At 4 magnitudes below peak it still contributes a couple of per
cent of the peak flux, which is enough to bias exactly the late-time points you
built the template for.

### The short version

```bash
mkdir -p ~/Desktop/TJO/finalfiles/snname   # frames here
tjophot list snname                        # check they are found
cp targets/TEMPLATE.yaml targets/snname.yaml
$EDITOR targets/snname.yaml                # set target.name
tjophot resolve snname --write             # position from TNS
tjophot fetch-catalog snname               # check the star counts
tjophot fetch-template snname              # LOOK AT THE PREVIEWS
tjophot run snname --max-images 6          # trial; check runs/*/qa/
./full-run.sh snname                       # the real thing
```

## Frames that came back with no zero point

A frame that finds too few catalogue stars gets `no_zeropoint` and no
magnitude, even though its flux was measured fine. Two things help.

**Fix the cause.** The catalogue-error cut is per survey
(`calibration.max_cat_err`), because Gaia synthetic photometry is good to
millimagnitudes while APASS B errors run 0.05-0.15. A single tight cut throws
APASS away: on the SN2026acd field it kept 2 of 27 B stars, where 0.12 keeps 21.

**Recover what is left, after the fact:**

```bash
tjophot recalibrate snname -n     # report only
tjophot recalibrate snname        # rewrite the tables and plots
```

This does *not* substitute the band's average zero point. The zero point tracks
atmospheric transparency and spans 1.4-2.6 magnitudes within a band over a
campaign, and calibration tends to fail on exactly the poor nights where an
average is furthest from the truth. What is stable is the *difference* between
two bands on the same night, since the extinction largely cancels: on SN2026acd
those offsets hold to 0.02-0.10 mag. So a missing zero point is predicted from
every other band observed that night and the predictions combined.

How far apart those predictions land is the useful part. Agreement means the
night's extinction was grey and the transfer is trustworthy; disagreement means
it was not, and that spread goes into the uncertainty. Tested by hiding each
fitted zero point in turn and predicting it: median error 0.030 mag, 90th
percentile 0.154, with 93 per cent of cases inside the quoted 2 sigma.

Recovered epochs are flagged `zp_transferred`, so they are easy to exclude.
Sort them by `zp_err`: below about 0.10 the night was grey and the value is
solid; above 0.25 treat it with real caution. The original table is kept as
`*_photometry_all.csv.bak`. Runs do this automatically at the end
(`calibration.transfer_zeropoints`).

## If something goes wrong

**"cannot find .../config/default.yaml"** — the package was installed away from
its config files. Reinstall with `pip install -e .` from the pipeline directory.

**"no FITS files found for X"** — check `tjophot list X` and `paths.data_root`.

**Catalogue queries fail** — check the network, then rerun `tjophot
fetch-catalog` and look at the source counts it prints. A run will fall back
through `gaia_synth` to `panstarrs` to `apass`, and reports `no_zeropoint` in
the output if none of them works.

**A run dies partway** — rerun with `--skip N --append` to pick up where it
stopped.
