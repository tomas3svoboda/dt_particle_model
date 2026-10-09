# Item 35: the Faner (2008) comparison re-scored against the version-2 digitization

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

Recorded 2026-09-27 (machine date). Brief: re-score every consumer of Fig. 4.21 of
Faner (2008) against `paper/analysis/datasets/faner2008_thesis_95_108_v2/` (version
2), re-march only where a march input changed. Every committed folder was only read;
nothing under `paper/01_particle_jfpe/`, `docs/`, `src/`, `scripts/`, `tools/` or any
other result item was written. No git operation was run by this item.

## 1. Which inputs derive from the curves, and which consumers needed what

Measured from the drivers (`v2common.py`, `item05_measured_side.py`,
`window_quantities.py`; every version-1 value recomputed and checked against the
committed files by float equality before any version-2 value was taken):

| input | where it enters | v1 | v2 | march input? |
|---|---|---|---|---|
| X0 (first point, 0.58 s) | stage-1 duration t_c, the transfer area (via 1+X0), the tail start | 0.7350938 | 0.7350938 (unchanged, byte-equal) | **yes, unchanged** |
| last measured loading | march target (hit time only) | 0.1041772 / 0.1005939 | unchanged | recorded only |
| constant-rate window and slope (X > 0.1888 + 0.05) | item 5: measured rate, duty, demanded h, every predicted-over-measured ratio | sun 10 pts, 0.010583 1/s; soy 10 pts, 0.010995 1/s | sun **11 pts, 0.009830 (-7.12 %)**; soy 10 pts, **0.010546 (-4.09 %)** | no (measured side) |
| measured crossing of X = 0.20 | item 5 crossing ratios; the scored window start (t >= crossing) | 57.699 / 54.930 s | **60.144 / 58.034 s (+2.445 / +3.104 s)** | no (measured side) |
| plateau Tp mean up to the crossing | item 5 demanded h | 69.134 / 69.980 C | 69.134 / **69.729 C** | no |
| falling-rate window points and reading bands | every window RMS, in-band count, span | 15 + 15 points | sun **14** (rule) or 15; soy 15 with 10 values changed | no |
| first sample of the window | item 28 / surface-record leverage window, hence the 95 % factor | 60.04 / 60.03 s | sun **70.11 s** (rule) / 60.04 s (15 pt); soy 60.03 s | **the band factor is a march input: checked, unchanged** |

**Consequence.** No march input moved: X0 is byte-identical, and the one derived
input that could have moved, the 95 per cent prediction factor on which the band
edges of items 27, 28, 30 and 31 are marched (item 28's F = 3.1401 sunflower,
3.1323 soybean at the boundary read; item 24's 3.1651 / 3.1698 at the outer-cell
read), is the maximum leverage over the window, reached at about 297 s, after the
window's new start; it is unchanged to the float at every level and both window
variants (`outputs/window_quantities.json`, `factors_that_moved`: only the
outer-cell factors of the sunflower `source_eq426` and `source_low` closures move,
3.1711 -> 3.1610 and 3.1703 -> 3.1585, neither the input of any march). **Every
consumer was therefore re-scored, none re-marched for a changed input.**

**The only marches run here** are twelve item 27 thesis cells (60/4, 120/4, 240/4,
60/8 and the two 240/4 band edges on each trace) whose per-step files the committed
folder did not keep (`large_dumps.sha256`). They were re-run through the committed
driver, redirected to `remarch/`, only to recover their trajectories: all twelve
reproduce the committed per-step CSV byte for byte (sha256 and size) and every
committed summary leaf except wall time (`log_remarch.txt`; 36 to 156 s each, 12:27
to 12:46 local, no host sleep).

| consumer | what changed on its measured side | treatment |
|---|---|---|
| item 5 rig (2026-09-22), stage 1 | slope -7.1 / -4.1 %; crossing +2.4 / +3.1 s; soybean plateau -0.25 C | re-computed on the measured side (predicted rates and t_c unchanged) |
| item 5 rig, falling-rate window (`rerun_traces.csv`, 88 thesis traces) | points, bands, window start | re-scored |
| item 5 surface record (2026-09-24) | window first sample (sunflower) | recomputed from its per-step record |
| item 5 space-resolved (2026-09-21c) | points, bands, window | re-scored (stored traces) |
| items 24, 27, 28, 29, 30, 31 | points, bands, window start | 300 thesis runs re-scored from their per-step records (12 via the byte-identical re-march); every committed analysis script re-run unchanged on the re-scored runs |
| item 32 | points, bands, window | 48 runs re-scored (stored traces and steps) |
| the identified scalar 1.1654e-9 m2/s (item 5) | its fitting data | its two stored traces re-scored; **not re-identified** (a new fit is outside this brief) |

## 2. The sunflower 60.04 s point, and which window the paper carries

The version-2 sunflower point at 60.0409 s reads X = 0.20027, 0.0003 (0.1 px) above
0.20; the version-2 crossing is 60.144 s. **The committed rule, applied exactly as
written (every measured sample with t >= the measured crossing), excludes it: the
sunflower window is 14 points from 70.11 s, span 0.06975 kg/kg.** Under the unchanged
rule the paper carries the 14-point window, and the headline column below is that
one. The 15-point window (the 60.04 s point kept; span 0.09609) is reported beside
it at every number, because whether the point is in or out is decided below reading
resolution. The soybean window is the same 15 points in both (span 0.07912 ->
0.09138 because its first point moved from 0.1797 to 0.1920).

## 3. Verification, run before any version-2 value was read

* **Re-scoring identity:** 5,850 checks exact over 300 runs of items 24 to 31 (every
  leaf of both committed scorings, decimated "as printed" and full-step, the measured
  crossing, and items 30/31's stored decimated traces element by element;
  `outputs/rescore_identity.json`).
* **Analysis harness identity:** the committed `analyze_item24/27/29/30/31.py`, run
  unchanged in a disposable shadow copy of the results tree with nothing replaced,
  reproduce all 31 committed analysis outputs byte for byte
  (`outputs/shadow_v1/HARNESS_IDENTITY.json`).
* **Item 5 measured side:** 707 checks exact (slopes, windows, crossings, plateaux,
  duties, demanded h, every film case's ratios, the circularity probe, all 88 stored
  tail residuals); the identified-scalar traces (stored at 10 significant digits)
  reproduce 0.1584 / 0.2645 to the printed digits.
* **Window statistics:** 744 checks exact (surface record and item 28).
* **Space-resolved:** 928 checks exact over 52 traces (items 5-dxr and 32).
* **Re-march:** 12 of 12 per-step files byte-identical to the committed sha256.

## 4. The headline, before and after

Window RMS as a fraction of the window span; "as printed" = the committed decimated
scoring at 60 cells and the reported step; in band = points inside the reading band.

| quantity | sunflower v1 | **sunflower v2 (rule, 14 pts)** | sunflower v2 (15 pts) | soybean v1 | **soybean v2** |
|---|---:|---:|---:|---:|---:|
| declared sphere (0.7215 / 0.666 mm), as printed | 13.3 % (10 of 15) | **17.2 % (9 of 14)** | 12.6 % (9 of 15) | 19.7 % (6 of 15) | **15.6 % (8 of 15)** |
| converged (order one to observed) | 13.04-13.05 | **16.67-16.70** | 12.52 | 18.89-18.92 | **15.17-15.19** |
| numerical uncertainty of the printed value | +0.24 to +0.26 pts | **+0.51 to +0.53** | +0.05 | +0.79 to +0.82 | **+0.43 to +0.45** |
| mesh orders (r = 1, 2, 4) | 1.91-1.98 | **2.00-2.00** | 1.67-1.92 | 2.01-2.04 | **2.05-2.12** |
| time orders (60, 120, 240 cells) | 1.10-1.11 | **1.10-1.11** | 1.08-1.10 | 1.16-1.19 | **1.15-1.16** |
| 95 % range (full-step) | 5.1-40.4 | **5.4-52.6** | 6.2-37.0 | 8.9-50.3 | **7.7-42.2** |
| band minimum (parabola), D multiplier | 1.45 | **1.46** | 1.44 | 1.56 | **1.50** |
| continuum mean signed residual / span | +0.111 | **+0.158** | +0.095 | +0.178 | **+0.134** |
| 0.885 mm sphere (volume-equivalent, 2019 d), as printed | 24.4 (2) | **31.7 (1 of 14)** | 22.4 (2) | 36.6 (0) | **30.3 (1)** |
| 0.885 mm converged | 24.0 | **31.2** | 22.2 | 35.8 | **29.7** |
| R_d = 0.655 mm (2019 d), as printed | 7.9 (15) | **9.9 (14 of 14)** | 8.1 (14) | 18.7 (7) | **14.7 (9)** |
| thesis volume-equivalent sphere (psi = 1) | 29.1 (2) | **37.9 (1)** | 26.7 (2) | 37.5 (0) | **31.0 (1)** |
| crossing ratio, simulated over measured (every sphere) | 0.89 | **0.86** | (not a crossing) | 1.01 | **0.95** |

**Every ladder re-establishes its orders.** Under the rule window the declared-sphere
ladder is second order in the mesh (2.00 sunflower, 2.05 to 2.12 soybean) and first
order in the step (1.10 to 1.16), mesh-independent within the time-step error
(largest 60-to-240 change 0.0012 and 0.0008 of the span against a 240/1-to-240/4
change of 0.0016 and 0.0014); items 27 and 30 likewise (0.885 mm: mesh 1.97 to 2.02,
time 1.02 to 1.16 over all four conditions; R_d: mesh 2.01 to 2.11, time 1.07 to
1.16). No rejecting outcome of items 27, 30 or 31 changes state (item 27's R3 still
fires on its horizon-remainder clause, a march fact untouched by scoring). In the
15-point variant the sunflower mesh changes are 0.0003 to 0.0004 of the span, below
item 31's 0.0005 materiality threshold, where the observed order (1.67 to 1.92) is
not tested; at R_d one sunflower mesh chain (r = 4) does not form in that variant for
the same reason (changes of order 1e-5).

**What the move is.** The corrected version-2 points lie 0.008 to 0.028 higher on the
early falling-rate segment; under the unchanged rule the sunflower window loses its
first point and its span shrinks by 20 per cent, so the same trajectory's residual
rises on sunflower and falls on soybean. The two traces trade places: soybean is now
the closer trace at the declared spheres. The march still dries more slowly than both
traces (continuum mean signed residual positive on both).

## 5. Item 5: the constant-rate leg

| quantity (Whitaker, charge-mass basis unless stated) | v1 | **v2** |
|---|---|---|
| duty, simulated over measured: 2008 sun, 2008 soy, 2019 soy, 2019 sun | 0.98, 0.88, 0.91, 0.82 | **1.06, 0.92**, 0.91, 0.82 |
| range printed "0.82 to 0.98" | 0.82-0.98 | **0.82-1.06** |
| crossing, simulated over measured, same four | 0.89, 1.01, 0.97, 1.12 | **0.86, 0.95**, 0.97, 1.12 |
| range printed "0.89 to 1.12" | 0.89-1.12 | **0.86-1.12** |
| dry-meal reading: duty / crossing ranges | 0.51-0.61 / 1.45-1.74 | **0.53-0.61 / 1.45-1.68** |
| demanded h (bed area), 2008 sun, soy, W/m2K | 97.1, 85.5 | **90.2, 81.7** |
| Whitaker over the demand, 2008 sun, soy (2019: 0.909, 0.873) | 0.998, 0.885 | **1.074, 0.927** |
| range printed "0.87 to 1.00" | 0.87-1.00 | **0.87-1.07** |
| Bird and Bradshaw-Myers over the demand (printed "1.8 to 2.6") | 1.78-2.55 | **1.78-2.74** |
| demanded h, dry reading, 2008 soy, sun | 148.4, 168.5 | **141.8, 156.5** |
| particle-area demand (S13.2, range over the four conditions) | 16.18-18.65 | **16.18-17.82** |
| specific duty, 2008 sun, soy, W/kg dry | 3544.7, 3682.6 | **3292.2, 3532.1** |
| charge-reading flux definition against the plotted slope, closest miss, sun, soy | 10.2 %, 4.8 % | **10.0 %, 0.7 %** |

**Stated because it touches a written rejecting outcome.** Item 5's M2 wrote before
its arithmetic that a failure of the three printed bed correlations to bracket the
measured demand would withdraw the stage-1 posing. At version 1 Whitaker was below
the demand at all four conditions and Bird and Bradshaw-Myers above. At version 2 the
2008 sunflower demand (90.2) lies **below all three** (Whitaker 96.9, 1.074 times the
demand); the other three conditions still bracket. This record does not rule on the
consequence; it is a Class-B question for the lead or the owner. The duty ratio on
that trace moves from "slightly low" (0.98) to "slightly high" (1.06).

## 6. The numbers the paper prints from these items, before and after

One row per number the particle paper prints from items 5, 24, 27, 28, 29, 30, 31
and 32 (the ids are those of `paper/01_particle_jfpe/check_numbers_2026-09-23.py`;
multi-number rows list their values in the order of the "quantity" column). "v1
recomputed" is this item's own recomputation from the committed runs and was checked
against the printed digits; where it differs in the last digit the printed rounding
is noted in DETAIL. Rows marked "unchanged" do not depend on the Fig. 4.21 curves
(radii, identity counts, ledgers, model-only quantities, the Faner 2019 conditions,
refusal times) or are consequences of the version-1 identification, which was not
re-run. Source of every value: `outputs/printed_numbers_v1_v2.csv`, built by
`build_printed_table.py` from `outputs/`.

| id | printed | quantity | v1 recomputed | **v2, rule (14 sunflower pts)** | v2, 15 pts | note |
|---|---|---|---|---|---|---|
| F1 | 68.7145 | constant-rate solid temperature (model) | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (saturation curve) |
| F2/F3 | -0.137, +0.093 | deviation from the Faner 2019 Fig. 2 plateaux | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (Faner 2019 figure) |
| F4/F5 | 0.19935, 1.257 | critical loading, GAB activity | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (model) |
| F6/F7 | 45.9, 131.3 | frozen pair Re_eps and h | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (correlation) |
| tab:basis | 0.82-0.98 | duty, simulated over measured, Whitaker, 4 conditions, charge basis | 0.82, 0.98 | **0.82, 1.06** | 0.82, 1.06 |  |
| 5.1 | 0.98, 0.88 | duty ratio, Whitaker: 2008 sunflower, 2008 soybean | 0.98, 0.88 | **1.06, 0.92** | 1.06, 0.92 |  |
| 5.1 | 0.91, 0.82 | duty ratio, Whitaker: 2019 soybean, 2019 sunflower | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (Faner 2019 curves) |
| tab:basis | 0.89-1.12 | crossing, simulated over measured, Whitaker, 4 conditions | 0.89, 1.12 | **0.86, 1.12** | 0.86, 1.12 |  |
| tab:fanermarch | 0.89, 1.01 | crossing ratio, Whitaker: sunflower, soybean (every sphere) | 0.89, 1.01 | **0.86, 0.95** | 0.86, 0.95 |  |
| B1/B2 | 0.51-0.61 | duty ratio, Whitaker, dry-meal reading, 4 conditions | 0.51, 0.61 | **0.53, 0.61** | 0.53, 0.61 |  |
| B3/B4 | 1.45-1.74 | crossing ratio, Whitaker, dry-meal reading, 4 conditions | 1.45, 1.74 | **1.45, 1.68** | 1.45, 1.68 |  |
| 5.1 | 0.87-1.00 | Whitaker h over the h the measured duty demands, 4 conditions | 0.87, 1.00 | **0.87, 1.07** | 0.87, 1.07 |  |
| 5.1 | 1.8-2.6 | Bird and Bradshaw-Myers h over the demanded h, 4 conditions | 1.78, 2.55 | **1.78, 2.74** | 1.78, 2.74 |  |
| 5.1 | 2.5 | spread of the three correlations (h ratio at one state) | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| RESULT_05 | 97.1, 85.5 | h demanded by the measured duty, bed area, W/m2K: sunflower, soybean | 97.1, 85.5 | **90.2, 81.7** | 90.2, 81.7 |  |
| RESULT_05 | 0.998, 0.885 | Whitaker over the demand: sunflower, soybean | 0.998, 0.885 | **1.074, 0.927** | 1.074, 0.927 |  |
| RESULT_05 | Whitaker below at 3 of 4; the three bracket the demand at all four | does Whitaker lie below the demand (sunflower, soybean) [journal conditions unchanged: below] | 1, 1 | **0, 1** | 0, 1 | 1 = below the demand, 0 = above; with 0 at sunflower all three comparators lie above the demand there |
| S10.8 | 129.0-168.5 | demanded h on the dry-meal reading, thesis traces (soybean, sunflower) | 148.4, 168.5 | **141.8, 156.5** | 141.8, 156.5 | 129.0 is a journal condition (unchanged); the thesis values are shown |
| F19/N7 | 16.18-18.65 | particle-area demand h, S13.2 (range over 4 conditions; 18.65 = 2008 soybean) | 16.18, 18.65 | **16.18, 17.82** | 16.18, 17.82 |  |
| R10 | 17.228 | particle-area demand h, the live S13.2 value (2008 sunflower in v1 is 17.73; 17.228 is a journal condition) | 17.732 | **16.469** | 16.469 | shown: 2008 sunflower; the printed 17.228 belongs to a Faner 2019 condition and is unchanged |
| R9 | 15.98 | 2026-09-19 particle-area demand (superseded record) | unchanged | **unchanged** | unchanged | item05_constant_rate.csv of 2026-09-19; not re-scored (superseded construction) |
| F20/F21 | 1671.4, 2.9231 | duty low (journal), frozen particle area | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| S13.2 | 3544.7, 3682.6 | specific duty W/kg dry: 2008 sunflower, soybean | 3544.7, 3682.6 | **3292.2, 3532.1** | 3292.2, 3532.1 |  |
| M2/M3 | 1 to 10 % | flux definition reproduces the plotted slope (closest |implied/measured - 1|, charge reading): sunflower, soybean | 10.2, 4.8 | **10.0, 0.7** | 10.0, 0.7 |  |
| F9/F10 | 9.9, 20.6 | window RMS/span, superseded outer-cell read (lit_whitaker, surface law), sunflower, soybean | 9.9, 20.6 | **12.7, 16.4** | 9.7, 16.4 |  |
| F9 | 14 of 15 | in band, same, sunflower | 14 | **13** | 13 |  |
| F17 | 16.3 | soybean window, source coefficient (source_central) | 16.3 | **14.0** | 14.0 |  |
| F18 | 9.7, 15 of 15 | sunflower window, frozen pair at the rig state | 0.0968, 15.0000 | **0.1175, 14.0000** | 0.0987, 14.0000 |  |
| A27/A28 | 0.141, 0.11 | surface edge, borrowed leg (frozen_1.0): soybean, sunflower | 0.1413, 0.1102 | **0.1354, 0.0994** | 0.1354, 0.1179 |  |
| F11 context | 15.84, 26.45 | identified scalar (1.1654e-9, fitted on v1 sunflower points): fitted-trace and held-out residual | 15.84, 26.45 | **20.59, 21.50** | 14.54, 21.50 | traces stored at 10 significant digits; the identification itself is not re-run (a fit) |
| F11/F12 | 1.165e-9, 2.91 | identified scalar and its ratio to PHY-019 | unchanged | **unchanged** | unchanged | a fit to the version-1 sunflower window; NOT re-identified (would be a new fit) |
| F13-F16 | 0.410, 1.507, 3.513e-9, 0.959 | law regression | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (Cardarelli data) |
| N15 | 96.8 | surface record: lowest duration fraction of the window inside the measured loading range, ten closures | 96.8 | **96.8** | 96.8 |  |
| N16 | 0.29-0.74 | surface record: leverage range over the window (low end over all twelve closures, high end over the ten particle-mass closures, as the record states it) | 0.29, 0.74 | **0.29, 0.74** | 0.29, 0.74 |  |
| N18 | 11.53 | surface record: volume mean over the measured maximum, largest, ten closures | 11.53 | **11.53** | 11.53 |  |
| B5/N20 | 3.17-3.18 | surface record: 95 % factor at the window's maximum leverage, the six rows the record tabulates (frozen_1.0, lit_whitaker, source_central) | 3.165, 3.184 | **3.165, 3.184** | 3.165, 3.184 | over all ten particle-mass closures: 3.165-3.196 (v1), 3.159-3.196 (v2 rule); only sunflower source_eq426 and source_low move, neither a band input |
| N17/N19 | 0.30, 0.72 | largest fitted leverage; temperature extrapolation | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| B6/B7, S39/S40 | 0.770, 0.697 | space-resolved D(X(r)), paper configuration (0.885 mm, borrowed pair): soybean, sunflower | 0.770, 0.697 | **0.682, 0.899** | 0.682, 0.639 |  |
| K1/K2 | 24.4, 36.6 | 0.885 mm sphere, reported cell 60/1, as printed: sunflower, soybean | 24.4, 36.6 | **31.7, 30.3** | 22.4, 30.3 |  |
| K10 | 2, 0 of 15 | in band, as printed: sunflower, soybean (of 14 / 15 under the v2 rule) | 2, 0 | **1, 1** | 2, 1 |  |
| K3/K4 | 24.0, 35.8 | 0.885 mm converged (continuum, order one to observed) | 24.0, 35.8 | **31.2, 29.7** | 22.2, 29.7 |  |
| K5/K6 | 1.96-2.03 | 0.885 mm mesh orders, all four conditions, three time levels | 1.96, 2.03 | **1.97, 2.02** | 1.97, 2.03 |  |
| K7/K8 | 1.05-1.16 | 0.885 mm time orders, all four conditions | 1.05, 1.16 | **1.02, 1.16** | 1.02, 1.15 |  |
| K9 | 0.4, 0.8 | 0.885 mm reported level above continuum, points (as printed minus order-one continuum) | 0.4, 0.8 | **0.6, 0.5** | 0.2, 0.5 |  |
| K11-K14 | 33.9, 13.2 | one-sigma edges (D/1.507, D x 1.507), sunflower, full-step | 33.9, 13.2 | **44.1, 17.0** | 31.1, 12.6 |  |
| K11-K14 | 46.6, 24.2 | one-sigma edges (D/1.507, D x 1.507), soybean, full-step | 46.6, 24.2 | **39.0, 19.7** | 39.0, 19.7 |  |
| K15-K19 | sun 47.6 (1), 14.7 (8), min 5.7 (15) at F^0.75 | item 27 band at F=3.165, sunflower, full-step: D/F, D x F, D x F^0.75 | 47.6, 14.7, 5.7 | **61.9, 18.3, 6.1** | 43.4, 14.2, 6.8 |  |
| K15-K19 | sun in band 1, 8, 15 | same, points in band | 1, 8, 15 | **0, 8, 14** | 1, 8, 14 |  |
| K16/K19 | soy 61.4 (0), 10.7 (12) | item 27 band, soybean, full-step: D/F, D x F | 61.4, 10.7 | **51.8, 9.9** | 51.8, 9.9 |  |
| K16/K19 | soy in band 0, 12 | same, points in band | 0, 12 | **1, 11** | 1, 11 |  |
| K20 | 0.42, 0.51 | item 27 band widths over the span (at least) | 0.42, 0.51 | **0.56, 0.43** | 0.37, 0.43 |  |
| K21-K23 | 2.78-3.38, 9.7-9.8 | threshold (first-stride loading change over Xc) | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (model only) |
| K24 | 54,658 | identity count | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| K25/K26 | 66.6, 60.2 | journal conditions, boundary read | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (Faner 2019 curves) |
| K27 | 1.51, 1.28, 1.15 | outer-cell over face loading | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (model only) |
| K28/K29 | 19.7, 30.8 | outer-cell read at the finest level 240/4, full-step: sunflower, soybean | 19.7, 30.8 | **25.4, 25.6** | 18.4, 25.6 |  |
| K30 | 490 | identity count | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| A1-A3, L12/L13 | 3.1401, 3.1323 | 95 % factor at the boundary read (window max leverage, 60/1) | 3.1401, 3.1323 | **3.1401, 3.1323** | 3.1401, 3.1323 | the band marches' input; unchanged, so no band re-march |
| A4-A7 | 47.4, 14.4, 61.2, 10.4 | item 28 band edges at 0.885 mm (full-step): sun D/F, sun D x F, soy D/F, soy D x F | 47.4, 14.4, 61.2, 10.4 | **61.7, 17.9, 51.7, 9.6** | 43.3, 14.0, 51.7, 9.6 |  |
| tab:fanermarch | 5.7-47.4, 10.4-61.2 | 0.885 mm 95 % range as tabulated: sun (item 27 minimum, item 28 slow edge), soy (item 28 edges) | 5.7, 47.4, 10.4, 61.2 | **6.1, 61.7, 8.5, 51.7** | 6.8, 43.3, 8.5, 51.7 | lowest marched point inside the band on each trace; the sweep is item 27's (F=3.165) |
| A8/A9 | 0.20-0.69 | leverage range over the window at the boundary read, 60/1, both traces | 0.20, 0.69 | **0.20, 0.69** | 0.20, 0.69 |  |
| A10/A12 | 1.032, 7.920e-4 | lowest boundary loading in the sunflower window over the measured minimum; its value | 1.031635, 0.000792 | **1.031635, 0.000792** | 1.031635, 0.000792 |  |
| A11 | 33,248 | identity count | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| H2/H3 | 0.677, 0.601 | radius factor minimising the window RMS (parabola), sunflower, soybean | 0.677, 0.601 | **0.672, 0.612** | 0.681, 0.612 |  |
| H9/H10 | 2.15, 2.78 | constant factor on D minimising the RMS at 0.885 mm | 2.15, 2.78 | **2.18, 2.64** | 2.14, 2.64 |  |
| H11/H12 | sun 0.198-0.282, soy 0.326-0.395 | Ea +/- 16.1 kJ/mol at 0.885 mm, as printed: sun +, sun -, soy +, soy - | 0.198, 0.282, 0.326, 0.395 | **0.257, 0.366, 0.267, 0.329** | 0.182, 0.259, 0.267, 0.329 |  |
| H13 | 0.0093 | D/R^2 cross-grid interpolation agreement (largest |difference|) | 0.0093 | **0.0146** | 0.0066 |  |
| H4-H8, H14, H15, H18 | 44.81, 44.8, 4.294, 1.27, 32.65, 3.5e-13, 7.607e-4, 1.102/1.017 | implied Ea from the identified scalar; path identity; datum; radius factors | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves, or a consequence of the version-1 identification (not re-identified) |
| H16 | -0.080, -0.065 | Ea-pair slope per unit ln factor at 0.885 mm: sunflower, soybean | -0.080, -0.065 | **-0.103, -0.058** | -0.073, -0.058 |  |
| H17/L1/L2 | 0.49-0.55 | RMS per unit ln R at the reported point (as printed, one- and two-step, both traces) | 0.49, 0.55 | **0.47, 0.66** | 0.44, 0.47 |  |
| H19/H20 | 0.458; 0.093 (14) | soybean radius sweep: factor 1.2; minimum marched (factor 0.6) and its in-band count | 0.458, 0.093, 14.000 | **0.383, 0.077, 14.000** | 0.383, 0.077, 14.000 |  |
| L3/L4 | 11.4, 15.3 | soybean at the 0.575 and 0.620 mm spheres (bracketing psi 0.67), as printed | 11.4, 15.3 | **8.7, 11.9** | 8.7, 11.9 |  |
| L5/L6 | 0.575, 0.620 | marched radii | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| G1-G4 | 7.9 (15), 18.7 (7); 0.0791, 0.1868 | R_d = 0.655 mm reported cell, as printed | 0.0791, 0.1868 | **0.0993, 0.1473** | 0.0811, 0.1473 |  |
| G1-G4 | in band 15, 7 | R_d reported cell, points in band as printed | 15, 7 | **14, 9** | 14, 9 |  |
| G5-G8 | 7.8, 17.8 (0.0777-0.0779, 0.1784-0.1786) | R_d converged (order one, observed) | 0.0777, 0.0779, 0.1784, 0.1786 | **0.0940, 0.0943, 0.1428, 0.1429** | 0.0820, 0.0820, 0.1428, 0.1429 |  |
| G9 | 2.01-2.13 (mesh) | R_d mesh orders, both traces | 2.00, 2.13 | **2.01, 2.11** | 1.61, 2.57 |  |
| G10/G11 | 1.10-1.33 (time) | R_d time orders, both traces | 1.10, 1.33 | **1.07, 1.16** | 1.07, 1.87 |  |
| G34/G35 | +0.0013, +0.0083 | R_d numerical uncertainty (as printed minus continuum, low end) | 0.0013, 0.0082 | **0.0050, 0.0044** | -0.0009, 0.0044 |  |
| G12-G15 | 5.2-40.0, 9.1-49.6 | R_d 95 % range (full-step, lowest marched to the far edge) | 5.2, 40.0, 9.1, 49.6 | **5.4, 51.4, 7.6, 41.6** | 6.4, 36.9, 7.6, 41.6 |  |
| G16/G17/G48 | D x 1.20, D x 1.51; 0.676 | R_d band minimum (parabola) D multiplier; as a radius factor 0.74/sqrt(mult) | 1.199, 1.510, 0.676 | **1.209, 1.451, 0.673** | 1.191, 1.451, 0.678 |  |
| G36 | 0.35, 0.40 | R_d band width over the span | 0.35, 0.40 | **0.46, 0.34** | 0.30, 0.34 |  |
| G49/G50 | 0.0524 (15), 0.367 (2) | R_d band table, sunflower: lowest marched, slow edge (full-step) | 0.0524, 0.3667 | **0.0542, 0.4771** | 0.0642, 0.3355 |  |
| G37/G38 | 3.137; 1.036, 1.048 | R_d information-only factor (sunflower); lowest boundary loading over measured minimum | 3.137, 1.036, 1.048 | **3.137, 1.036, 1.048** | 3.137, 1.036, 1.048 |  |
| G41 | +0.053, +0.167 | R_d continuum mean signed residual | 0.053, 0.167 | **0.082, 0.124** | 0.043, 0.124 |  |
| G44/G45 | 0.00037, 0.00130; 0.00070, 0.00199 | R_d largest 60-to-240 change; 240/1-to-240/4 change | 0.00037, 0.00130, 0.00070, 0.00199 | **0.00100, 0.00081, 0.00155, 0.00138** | 0.00003, 0.00081, 0.00026, 0.00138 |  |
| G30/G31 | -0.033 to +0.027 | R_d population shifts (decimated, minus single), both traces | -0.033, 0.027 | **-0.035, 0.033** | -0.028, 0.023 |  |
| G42/G43 | 0.0546, 0.2136 | R_d population (a1) sunflower, (b) soybean (decimated) | 0.0546, 0.2136 | **0.0639, 0.1708** | 0.0631, 0.1708 |  |
| G39 | 0.029 | R2 distance outside the single-sphere range (largest) | 0.029 | **0.033** | 0.025 |  |
| G18-G20, G25-G29, G32, G33, G40, G46, G47 | 0.655, 1.31, 1.443, 1.332, 1.95, 0.9, 0.32, 33,330, 68, 268.6/292.9, 5.5e-10, 2.25 | radii, diameters, Fourier number, counts, model times and coefficients | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| G21/G22, L7/L8 | 29.1, 37.5 | thesis volume-equivalent sphere (psi = 1), as printed | 29.1, 37.5 | **37.9, 31.0** | 26.7, 31.0 |  |
| tab:fanermarch | 2 of 15, 0 of 15 | thesis volume-equivalent sphere, in band | 2, 0 | **1, 1** | 2, 1 |  |
| I1/I2/I33/I34, G23/G24 | 13.3 (10), 19.7 (6) | DECLARED SPHERES, reported cell 60/1, as printed: sunflower, soybean | 13.3, 19.7 | **17.2, 15.6** | 12.6, 15.6 |  |
| I1/I2 | 10, 6 of 15 | declared spheres, points in band as printed (of 14 / 15 under the v2 rule) | 10, 6 | **9, 8** | 9, 8 |  |
| I5/I6/I37/I38 | 13.0, 18.9 (0.1304-0.1305, 0.1889-0.1892) | declared spheres, converged (order one, observed) | 0.1304, 0.1305, 0.1889, 0.1892 | **0.1667, 0.1669, 0.1517, 0.1519** | 0.1252, 0.1252, 0.1517, 0.1519 |  |
| I39/I40 | +0.0024 to +0.0026, +0.0079 to +0.0082 (0.2, 0.8 points) | declared spheres, numerical uncertainty (as printed minus continuum) | 0.0024, 0.0026, 0.0079, 0.0082 | **0.0051, 0.0053, 0.0043, 0.0045** | 0.0005, 0.0005, 0.0043, 0.0045 |  |
| I35/I36 | 0.1308, 0.1898 | declared spheres, finest cell 240/4, full-step | 0.1308, 0.1898 | **0.1675, 0.1523** | 0.1254, 0.1523 |  |
| I7/I8 | 1.91-1.98 (sun), 2.01-2.04 (soy) | declared spheres, mesh orders at r = 1, 2, 4 | 1.91, 1.98, 2.01, 2.04 | **2.00, 2.00, 2.05, 2.12** | 1.67, 1.92, 2.05, 2.12 |  |
| I9 | 1.10-1.11 (sun), 1.16-1.19 (soy) | declared spheres, time orders at 60, 120, 240 cells | 1.10, 1.11, 1.16, 1.19 | **1.10, 1.11, 1.15, 1.16** | 1.08, 1.10, 1.15, 1.16 |  |
| I41/I42 | 0.00070, 0.00132; 0.00097, 0.00200 | declared spheres, largest 60-to-240 change; 240/1-to-240/4 change | 0.00070, 0.00132, 0.00097, 0.00200 | **0.00118, 0.00084, 0.00161, 0.00140** | 0.00040, 0.00084, 0.00061, 0.00140 |  |
| I43 | +0.111, +0.178 | declared spheres, continuum mean signed residual | 0.111, 0.178 | **0.158, 0.134** | 0.095, 0.134 |  |
| I10-I13 | 5.1-40.4, 8.9-50.3 | declared spheres, 95 % range (full-step) | 5.1, 40.4, 9.0, 50.3 | **5.4, 52.6, 7.7, 42.2** | 6.2, 37.0, 7.7, 42.2 |  |
| I48 | 5.0-40.6, 9.3-50.5 | declared spheres, 95 % range (decimated) | 5.0, 40.6, 9.3, 50.5 | **5.5, 52.8, 7.7, 42.4** | 6.1, 37.0, 7.7, 42.4 |  |
| I14/I15 | D x 1.45, D x 1.56 | declared spheres, band minimum (parabola, full-step), D multiplier | 1.45, 1.56 | **1.46, 1.50** | 1.44, 1.50 |  |
| I44-I47 | 0.0508 (15), 0.0896 (14); 0.404 (1), 0.334 (5) | declared spheres band table: lowest marched sun, soy; slow edge sun; fast edge soy | 0.0508, 0.0896, 0.4045, 0.3335 | **0.0537, 0.0773, 0.5262, 0.3013** | 0.0624, 0.0773, 0.3697, 0.3013 |  |
| I44-I47 | in band 15, 14; 1, 5 | same, points in band | 15, 14, 1, 5 | **14, 14, 0, 3** | 14, 14, 1, 3 |  |
| I49/I16 | 3.138, 3.130; 1.034, 1.047 | declared spheres, information-only factor; lowest boundary loading over measured minimum | 3.138, 3.130, 1.034, 1.047 | **3.138, 3.130, 1.034, 1.047** | 3.138, 3.130, 1.034, 1.047 |  |
| I17-I20, I59, I60 | sun 10.0-17.7 (Ea+ 0.100 (12)), soy 18.0-22.4 (Ea- 0.224 (3)) | declared spheres, Ea +/- 16.1 kJ/mol, as printed: sun +, sun -, soy +, soy - | 0.100, 0.177, 0.180, 0.224 | **0.130, 0.227, 0.142, 0.182** | 0.095, 0.166, 0.142, 0.182 |  |
| I59/I60 | 12; 3 | in band: sun Ea+, sun Ea-, soy Ea+, soy Ea- | 12, 5, 7, 3 | **11, 4, 8, 4** | 12, 4, 8, 4 |  |
| L10/L11 | 0.069, 0.137 | mean signed residual at Ea + 16.1 (full-step): sunflower, soybean | 0.069, 0.137 | **0.097, 0.098** | 0.057, 0.098 |  |
| I61/I62 | 0.073 (15), 0.355 (0) | constant matched to Ea+ (sun), to Ea- (soy), as printed | 0.073, 0.355 | **0.084, 0.293** | 0.080, 0.293 |  |
| S38 table | 1.5: 0.050 (15), 0.095 (13); 44.81: 0.198, 0.259; transported: 0.131, 0.192 | declared spheres: constant 1.5; the 44.81 kJ/mol case; the transported identification (sun, soy each) | 0.050, 0.095, 0.198, 0.259, 0.131, 0.192 | **0.052, 0.076, 0.258, 0.219, 0.170, 0.155** | 0.062, 0.076, 0.182, 0.219, 0.122, 0.155 |  |
| I50-I52 | -0.072, -0.042; -0.190, -0.243; 0.38, 0.17 | declared spheres levers (decimated): Ea pair; matched constant pair; ratio | -0.072, -0.042, -0.190, -0.243, 0.381, 0.171 | **-0.091, -0.039, -0.258, -0.196, 0.355, 0.197** | -0.067, -0.039, -0.162, -0.196, 0.417, 0.197 |  |
| I53/I54 | 1.2e-14, 5.5e-14; 2.1e-16, 5.8e-16 | D/R^2 identity, |RMS difference| (decimated) | 1.1e-14, 3.7e-16, 5.2e-14, 7.8e-16 | **1.5e-14, 5e-16, 4.6e-14, 1.1e-16** | 1e-14, 3.1e-16, 4.6e-14, 1.1e-16 |  |
| I25/I26 | -0.0329, +0.0269 (0.033, 0.027) | declared spheres population shifts (decimated), range over both traces | -0.0329, 0.0269 | **-0.0406, 0.0322** | -0.0281, 0.0235 |  |
| I63/I64 | 0.1025 (14), 0.2240 (3) | population (a1) sunflower, (b) soybean (decimated) | 0.1025, 0.2240 | **0.1315, 0.1797** | 0.0998, 0.1797 |  |
| I66 | +0.017 | item 30 unscaled (b) against the soybean species-own single sphere | 0.017 | **0.015** | 0.015 |  |
| I67/L9 | 0.67 (reproduces 7.9 %) | sunflower sphericity at which the radius sweep gives item 30's R_d residual | 0.67 | **0.67** | 0.67 | linear interpolation of item 29's sweep (as printed); approximate by construction |
| I3/I4, I21-I24, I27-I32, I55-I58, I65, I68 | 0.7215, 0.666, 2.85, 2.43, 32.3, 27.4, 0.3, 17,043, 82, 0.81525, 0.75254, 1.9e-12, 4.293, 7.74e-10, 32.34, 2.854, 1.1017, 1.95 | radii, transported identification (a consequence of the version-1 identification), counts, ledgers | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves |
| S1/S2 | 0.736 (2), 0.875 (0) | space-resolved, sunflower at 0.7215 mm, 60 cells: Whitaker, borrowed (decimated) | 0.735, 0.875 | **0.956, 1.129** | 0.672, 0.801 |  |
| S1/S2 | in band 2, 0 | same, points in band | 2, 0 | **1, 0** | 1, 0 |  |
| S10/S11 | 0.564, 0.538 | space-resolved, 0.885 mm, Whitaker leg, 60 cells: sunflower, soybean | 0.564, 0.538 | **0.734, 0.473** | 0.516, 0.473 |  |
| S20/S21 | -0.37 to -0.82 | space-resolved mean signed residual where the march covers the window (decimated) | -0.37, -0.82 | **-0.34, -1.07** | -0.34, -0.76 |  |
| S22/S23 | 0.95-1.0 at the window start; 0 to 0.42 at its end | space-resolved volume fraction above the measured maximum in the window (covering marches): max of window maxima; max of window minima | 1.00, 0.42 | **1.00, 0.42** | 1.00, 0.42 |  |
| S3-S9, S12-S19, S24-S38 | 88.7, 147.1, 3732, 131.34, 96.89, 75.71, 300.1, refusal times, 97-101, 15-61, ledgers | refusal times, coefficients, ratios, ledgers, identity | unchanged | **unchanged** | unchanged | independent of the Fig. 4.21 curves (march outputs) |

## 7. Sentences the paper would carry (not applied; the paper is untouched)

Under the unchanged window rule (14 sunflower points):

* Abstract / Table 2 / Sec. 5.2: "on each sample's declared volume-to-surface
  sphere the falling-rate window residual is 17.2 and 15.6 per cent of the window
  span (sunflower, soybean; 9 of 14 and 8 of 15 points in the reading band),
  converged 16.7 and 15.2 per cent, numerical uncertainty 0.5 and 0.4 percentage
  points; second order in the shell mesh and first order in the step; the march dries
  more slowly than both traces."
* Comparator: "the 0.885 mm volume-equivalent sphere gives 31.7 and 30.3 per cent
  (converged 31.2 and 29.7); the thesis's own volume-equivalent spheres 37.9 and
  31.0 per cent; the 2019-diameter volume-to-surface sphere 9.9 and 14.7 per cent."
* Band: "carried at the law's 95 per cent prediction factor the residual ranges over
  at most 5.4 to 52.6 and 7.7 to 42.2 per cent, its minimum inside the band at
  D x 1.46 and 1.50."
* Temperature: "an activation energy inside the source's +/- 16.1 kJ/mol puts the
  residual between 13.0 and 22.7 per cent (sunflower) and 14.2 and 18.2 per cent
  (soybean), short of the band minimum."
* Population: "a declared 1 to 3 mm population moves the residual by -0.041 to
  +0.032 of the span."
* Sphericity: "sunflower 9.9 to 37.9 per cent between sphericity 0.67 and 1.0;
  soybean bracketed by the 0.575 and 0.620 mm spheres at 8.7 and 11.9 per cent and
  31.0 per cent at unit sphericity"; the residual moves by 0.47 to 0.66 of the span
  per unit log radius at 0.885 mm.
* Constant rate: "Whitaker's correlation gives 0.82 to 1.06 of the measured
  constant-rate duty and 0.86 to 1.12 of the measured transition time (dry-meal
  reading 0.53 to 0.61 and 1.45 to 1.68); at the rig's state it gives 0.87 to 1.07
  of the coefficient the measured duty demands, Bird's and Bradshaw and Myers's 1.8
  to 2.7 times it." The sentence "the three bracket the demand at all four
  conditions" is **false at version 2** (2008 sunflower: all three above) and must
  not be carried without the owner's ruling on item 5's M2 outcome.
* Table 3 crossing row: 0.86 and 0.95.

If the 15-point window were ruled instead, the sunflower column reads 12.6 per cent
(9 of 15), converged 12.5, band 6.2 to 37.0, 0.885 mm 22.4 per cent; the soybean
column is identical.

## 8. Claim and non-claim

**Claimed.** The version-2 digitization changes no march input (X0 byte-identical;
the 95 per cent band factor unchanged at every level and window); every consumer
of Fig. 4.21 in items 5, 24, 27, 28, 29, 30, 31 and 32 was re-scored from its stored
trajectories by the committed scorer and analysis code, after reproducing every
committed number exactly at version 1 (5,850 + 707 + 744 + 928 checks and 31 byte-
identical analysis outputs); the ladders re-establish second order in the mesh and
first order in the step at every geometry; the numbers of section 6 are the
version-2 values under the unchanged window rule, with the 15-point alternative
beside them.

**Not claimed.** That either window variant is the right one (the membership of the
60.04 s point is below reading resolution; the rule as written excludes it); a
re-identification of the scalar 1.1654e-9 m2/s or of anything derived from it (the
44.8, 32.3 and 27.4 kJ/mol, the factors 4.294, 2.85, 2.43), which are version-1 fits
and would need a new fit; a ruling on the bracketing outcome of item 5's M2; that
the corrected digitization is the source's instrument record (it is a 141 dpi
reading with declared envelopes); the run-to-table mapping of Tables 4.5/4.6 row 06;
anything about the Faner 2019 conditions (unchanged); calibration, validation,
physical qualification or plant prediction.

Detail, per-run tables and the pipeline: `DETAIL.md`.
