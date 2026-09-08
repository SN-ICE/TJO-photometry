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

From a clean slate, for one target:

```bash
./full-run.sh SN2026acd
```

Extra flags are split sensibly: options every subcommand understands
(`--config`, `--data-root`, `--work-root`, `--bands`, `--log-level`) are
forwarded to *all* the steps, and anything else goes to the final `tjophot run`.
Each step echoes the command it is about to run, so you can see exactly what
happened. It is a thin wrapper over five commands you can equally run by hand:

```bash
# 1. clear out anything left from a previous run of this target
tjophot clean SN2026acd

# 2. look the position up on TNS and record it in targets/SN2026acd.yaml
tjophot resolve SN2026acd --write

# 3. cache the reference catalogues for the field          (internet, once)
tjophot fetch-catalog SN2026acd

# 4. fetch the reference images for difference imaging     (internet, once)
tjophot fetch-template SN2026acd

# 5. measure the light curve
tjophot run SN2026acd
```

Steps 2-4 are cached, so they cost nothing on a re-run. Step 1 only removes
`runs/SN2026acd/` -- the downloads and the detector patterns survive, which is
what you want when you are iterating on the photometry.

`tjophot run --clean` folds step 1 into step 5, so the shortest honest
re-run once the field is set up is:

```bash
tjophot run SN2026acd --clean
```

**Detector patterns** are a property of the instrument, not of the target, so
they sit outside this recipe. Build them once and forget them:

```bash
tjophot build-pattern SN2026acd --extra-dirs ~/Desktop/TJO/2023ixf/ph
```

Rebuild only if the camera changes -- a new CCD, a new filter, different
binning. They live in `calib/` and are applied automatically.

**Starting genuinely from scratch**, including re-downloading everything:

```bash
tjophot clean SN2026acd --all      # add -n first to see what would go
```

`--all` also deletes this field's cached catalogues and reference images;
`--templates` and `--catalogs` pick one or the other. `-n` / `--dry-run` lists
what would be removed without touching anything, and the command refuses to
delete outside `runs/` and `cache/`.

A full run over SN2026acd's 155 frames takes roughly 20 minutes.

## Pan-STARRS templates

Difference imaging needs a reference image of the field without the supernova
in it. For g, r, i and z at Dec > -30 that comes from Pan-STARRS DR2, and it is
already the default -- `tjophot run` will download a cutout the first time it
needs one. `tjophot fetch-template` just does that step on its own so you can
check it before committing to a run:

```
  g: panstarrs  2700x2700 px, 0.250"/px, 100% valid, median 0.3, range -3.1 to 892.4
     cached at  cache/templates/ps1_183.05315+13.22175_g-sdss_2700.fits
     preview    runs/SN2026acd/templates/template_g.png
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
   tjophot fetch-template SN2026acd --compare --bands B V
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
    B: ~/Desktop/TJO/refs/SN2026acd_ref_B.fits
    V: ~/Desktop/TJO/refs/SN2026acd_ref_V.fits
```

or let the pipeline stack late-time epochs itself with an MJD window under
`subtraction.tjo_stack`. Epochs that end up inside their own template are
flagged `in_template`.

A same-filter template that still contains the supernova is *worse* than a
clean cross-band one: it subtracts the signal you are trying to measure, and
nothing downstream can recover it. That is why PS1 is listed ahead of
`tjo_stack` for B and V in the SN2026acd target file.

## What you get

Everything lands in `runs/SN2026acd/`:

| file | what it is |
|---|---|
| `SN2026acd_lightcurve.csv` | the light curve: one row per epoch |
| `SN2026acd_photometry_all.csv` | every method for every epoch |
| `SN2026acd_lightcurve.png` | the light curve, plotted |
| `SN2026acd_diagnostics.png` | seeing, zero point and depth over the campaign |
| `qa/<image>.png` | per-epoch panel: stamps, curve of growth, zero-point fit |
| `run.log`, `run_summary.json`, `config_used.yaml` | what was actually run |

Start with `run.log` and `SN2026acd_diagnostics.png` when something looks wrong.

## The SN position

The pipeline looks the target name up and writes the answer into the target
file, so this only has to happen once:

```bash
tjophot resolve SN2026acd --write
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
INFO tjophot.config  configuration: default.yaml + tns.yaml + targets/SN2026acd.yaml
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
tjophot run SN2026acd --bands g r i          # a subset of bands
tjophot run SN2026acd --max-images 6         # a quick sanity run
tjophot run SN2026acd --no-subtraction       # direct photometry only
tjophot run SN2026acd --no-pattern           # skip the detector-pattern removal
tjophot run SN2026acd --no-astrometry        # keep the header WCS
tjophot run SN2026acd --refresh              # ignore cached catalogues/templates
tjophot run SN2026acd --bands z --append     # redo one band, merge into the table
tjophot run SN2026acd --skip 80 --append     # resume after an interruption
tjophot colorterms SN2026acd                 # fit the instrument colour terms
```

Every flag has a config equivalent; `--config my.yaml` merges a file of your own
over the defaults and the target file.

## A new supernova

```bash
mkdir ~/Desktop/TJO/finalfiles/SN2027xyz      # put its frames here
cat > targets/SN2027xyz.yaml <<'EOF'
target:
  name: SN2027xyz
  ra: "01:23:45.6"
  dec: "-12:34:56.7"
EOF
tjophot fetch-catalog SN2027xyz
tjophot run SN2027xyz
```

No new pattern build is needed — `calib/` already has them.

For difference imaging on a new target, Pan-STARRS templates are used
automatically for g, r, i and z (Dec > -30). B and V have no PS1 counterpart, so
either point `subtraction.reference_files.B` / `.V` at a real reference image in
the target file, or set an MJD window under `subtraction.tjo_stack` once the
supernova has faded. Without a template the pipeline reports direct PSF
photometry for those bands and says so.

## Frames that came back with no zero point

A frame that finds too few catalogue stars gets `no_zeropoint` and no
magnitude, even though its flux was measured fine. Two things help.

**Fix the cause.** The catalogue-error cut is per survey
(`calibration.max_cat_err`), because Gaia synthetic photometry is good to
millimagnitudes while APASS B errors run 0.05-0.15. A single tight cut throws
APASS away: on the SN2026acd field it kept 2 of 27 B stars, where 0.12 keeps 21.

**Recover what is left, after the fact:**

```bash
tjophot recalibrate SN2026acd -n     # report only
tjophot recalibrate SN2026acd        # rewrite the tables and plots
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
