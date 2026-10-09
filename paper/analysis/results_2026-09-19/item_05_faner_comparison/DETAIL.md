# Item 5 — quantitative comparison against the superheated-hexane experiments

Fills **declared missing result 5** as far as the declared domain allows
(`sec_missing.tex` item 5; `sec_validation.tex` §The Faner experiments;
`sec_results.tex` §Comparisons against the two literature anchors). It does not
lift the three model-side boundaries; it poses everything that can be posed
without them, and counts and names every element that stays refused.

## Method

`run_item05_faner.py` (wall time 0.10 s). Part A re-runs the committed
instrument `tools/analysis/faner2008_fig421_mass_biot.py` unchanged. Part B
applies that instrument's own `identify` to the Faner 2019 Figure 2 traces, the
only dimensional loading curves in that paper. Part C inverts the frozen
n-hexane saturation surface at 101 325 Pa by bisection and scores it against
the digitized solid-temperature plateaus, and inverts the measured
constant-rate specific demand into the external heat coefficient the rig
requires. `make_item05_figure.py` draws `item05_faner.pdf/.png`.

## (a) Faner 2008 Figure 4.21, 136 °C — the falling-rate form identified

| series | K/m (1/s) | samples | X_c crossed (s) | last X (kg/kg) | front fraction | observed traverse (s) | unhindered (s) | Bi_m | reading band | declared 1e5 predicts |
|---|---|---|---|---|---|---|---|---|---|---|
| soybean | 1.0995e-2 | 10 | 57.7 | 0.1006 | 0.7965 | 242.4 | 8.0 | 269.1 | 211.7–377.0 | 87 081 s = 359× observed |
| sunflower | 1.0583e-2 | 10 | 61.8 | 0.1042 | 0.8069 | 238.4 | 8.0 | 280.7 | 205.5–377.3 | 82 094 s = 344× observed |

Margin sweep (0.03/0.05/0.08 kg/kg): 252.7/269.1/300.3 and 268.6/280.7/293.4.
Alternate 334.5 K anchors: 206.1 and 200.1. Kill criterion K3: bands overlap
(211.7–377.0), both inside [1, 1e4], series form identifiable — **true**.

## (b) Faner 2019 Figure 2 — the same identification on the only invertible curves

Both whole traces are **refused** by the instrument, verbatim: `soybean: the
last sample 0.007415 is not inside (X_e, X_c)` and the same text for sunflower.
Both traces dry below the model's declared equilibrium anchor 0.0105 kg/kg.
Truncated to the samples above that anchor, the identification runs:

| material, gas T | anchors (X_c, X_e) | samples kept | K/m (1/s) | X_c crossed (s) | observed traverse (s) | unhindered (s) | Bi_m | declared 1e5 / observed |
|---|---|---|---|---|---|---|---|---|
| soybean, 120 °C | 0.20, 0.0105 (source X_c) | 20 | 7.264e-3 | 45.5 | 144.5 | 26.0 | 9.66 | 8 487× |
| soybean, 120 °C | 0.1888, 0.0105 | 20 | 7.264e-3 | 48.3 | 141.7 | 24.5 | 10.19 | 8 123× |
| soybean, 120 °C | 0.2017, 0.0254 | 15 | 7.264e-3 | 45.0 | 94.9 | 24.2 | 6.18 | 12 046× |
| sunflower, 100 °C | 0.20, 0.0105 | 23 | 5.096e-3 | 65.1 | 154.8 | 36.8 | 7.19 | 10 594× |
| sunflower, 100 °C | 0.1888, 0.0105 | 23 | 4.990e-3 | 68.6 | 151.3 | 35.4 | 7.39 | 10 365× |
| sunflower, 100 °C | 0.2017, 0.0254 | 19 | 5.096e-3 | 64.5 | 115.5 | 34.5 | 4.97 | 14 104× |

**The substantive finding.** Under one and the same closed-form series law, the
two sources identify mass Biot numbers a factor of roughly 30 apart: 269–281
from the 2008 thesis traces and 5–10 from the 2019 Figure 2 traces. The 2008
traces stop at a front fraction of 0.80 (the front has barely moved) while the
2019 traces run to 0.15–0.21, so the two data sets constrain different parts of
the same curve and are not consistent under a single Bi_m. The declared 1e5 is
two to four orders of magnitude above every identification. The 2019 dataset
carries no per-point reading bands, so no reading band is propagated for (b);
the anchor set is the sensitivity axis instead.

## (c) The zero-parameter constant-rate solid temperature, and the duty

Inverting the frozen n-hexane saturation curve at 101 325 Pa gives
**341.8645 K = 68.7145 °C**, reproducing the paper's value exactly.
Scored against every digitized solid-temperature sample taken while the loading
is above X_c = 0.20:

| material | plateau samples | plateau mean | spread | prediction − plateau | gas mean over the window | driving ΔT |
|---|---|---|---|---|---|---|
| soybean, 120 °C | 4 | 68.809 °C | 1.003 °C | −0.095 °C | 120.260 °C | 51.450 K |
| sunflower, 100 °C | 5 | 68.634 °C | 0.251 °C | +0.081 °C | 101.823 °C | 33.189 K |

Both deviations are inside the ±0.7 °C read uncertainty of that figure. The
selection rule here (every sample while X > X_c) differs from the three and four
markers quoted in `sec_validation.tex`, and gives 68.809 and 68.634 rather than
68.851 and 68.621; the spreads, 1.003 and 0.251 °C, match the ones the paper
already discloses.

**What the constant-rate period requires of the film.** Using only the measured
specific demand, the frozen Δh_vap at the predicted boiling point
(3.34929e5 J/kg), the source's own equivalent sphere diameter and apparent
particle density, and the measured gas-minus-solid difference, the external
heat coefficient the experiment demands is **16.18 W m⁻² K⁻¹** (soybean) and
**15.98 W m⁻² K⁻¹** (sunflower): specific duty 2433.0 and 1706.8 W per kg dry meal,
external area 2.9223 and 3.2184 m² per kg.

## The register of refused elements — 33, none skipped

`item05_refused_register.csv`. By class: 12 data-side not invertible (the six
Figure 3 sunflower and six Figure 4 soybean normalized experimental traces,
whose ordinate is (X − X_e)/(X₀ − X_e) with neither X₀ nor X_e printed per
run); 12 not data (the authors' own model lines on those figures); 6 instrument
refusals (the two Figure 2 whole traces at three anchor sets each); 3 model-side
declared-domain boundaries (pure-hexane atmosphere, water-free charge,
thin-layer film coefficient).

## Claim boundary

Nothing here is a trajectory of this model against those curves; the three
model-side boundaries stand. (a) and (b) are identifications of a *closed-form
falling-rate law's* one free group from digitized traces, not solver runs; (c)
is a phase-equilibrium prediction plus an inversion of the measured data, not a
rate comparison. The 16 W m⁻² K⁻¹ is what the rig requires, not a coefficient
for the model: the geometry is a thin layer, not a packed bed, and it may not
be substituted for the correlation of Table S1. Not physically qualifying, not
plant predictive.
