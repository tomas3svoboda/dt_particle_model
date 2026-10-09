# Item 31 — the Faner march at each species' own volume-to-surface sphere

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

## 0. Rejecting outcomes, written 2026-09-27 before any march of this item except the identity check

Brief (fifth referee read, majors M1 and M2): the Faner falling-rate march
(both thesis traces at 136 C; film case `lit_whitaker`, temperature case
`sensible`, the carried law read at the boundary value, item 27) at each
species' own volume-to-surface sphere, R_d,sun = 0.74 x 1.95 mm / 2 =
0.7215 mm for the sunflower trace and R_d,soy = 0.74 x 1.80 mm / 2 =
0.6660 mm for the soybean trace (section 1). The film transfer area per kg
dry is held fixed at the committed sample-holder value for every radius.

* **Q0 (identity first).** Item 30's two single reported-level runs at these
  radii (`conv_thesis_aris_<key>`, its conventions table: 13.3 and 19.7 per
  cent) are reproduced through this item's driver by exact float identity,
  every summary leaf except tag and wall time and every recorded step value.
  A single mismatch stops the script. *Run before this section was written:
  **17,043 of 17,043 checks exact over 2 marches** (`identity.json`); these
  two runs are the ladder's 60/1 cells.*
* **R1 (the ladder at either species-own radius is not second order in the
  mesh, or moves by more than the time-step error).** Item 30's R1 exactly:
  on each trace at its own radius and each time level r in {1, 2, 4}, the
  mesh chain 60/120/240 of the full-step window RMS over the span. R1 fires
  if any mesh chain does not form, or forms with an observed order below 1.5
  while its 60-to-240 change exceeds 0.0005 of the span; or if the largest
  60-to-240 change over the three time levels exceeds the 240/1-to-240/4
  change.
* **R2 (a level refuses inside the scored window).** Fires if any march of
  this item (ladder, band edges and sweep, temperature alternatives,
  similarity checks, population classes) ends before the last scored sample
  (300.1 s). Refusals after the window are recorded verbatim and do not fire
  R2. Nothing (floor, tolerance, retry budget, horizon, level) is relaxed.
* **R3 (the temperature alternatives at the species-own radii contradict
  item 29's ordering of the levers).** Item 29's ordering, per unit change of
  the logarithm at the reported point: radius (0.49 to 0.55 of the span)
  greater than a constant factor on D (0.25 to 0.27, and exactly minus
  one half of the radius's) greater than an activation energy pivoted at
  95 C (0.065 to 0.080, a third to a quarter of the constant factor's). R3
  fires on a trace if (i) the D/R^2 identity fails: the residual at
  (R_own, law) differs from that at (0.885 mm, law / rf^2), or the residual
  at (R_own, Ea + 16.1 kJ/mol) from that at (0.885 mm, same alternative /
  rf^2), by more than 1e-9 of the span (rf = R_own / 0.885 mm is not a power
  of two, so float identity is not expected; item 29 measured 3e-15 at
  rf = 0.5); or (ii) the Ea +/- 16.1 kJ/mol pair moves the residual, per unit
  change of ln(factor on the law at 136 C), by as much as or more than the
  constant pair matched to the same two factors at 136 C (1.694 and 0.590,
  marched at the same radius), in magnitude; or (iii) the two pairs move the
  residual in opposite directions.
* **Predicted by the D/R^2 similarity before the runs, not thresholds:** the
  constant-factor minimum at the species-own radius sits at item 29's minimum
  transported, D x (rf_own / rf_min)^2 = (0.8153 / 0.677)^2 = 1.45
  (sunflower) and (0.7525 / 0.601)^2 = 1.57 (soybean), inside the band; the
  identified constant 1.165e-9 m2/s at 0.885 mm is the same identification as
  1.165e-9 x rf^2 at R_own, factors 4.294 x 0.6646 = 2.85 and 4.294 x 0.5663
  = 2.43 on the law at 136 C, positions about 0.92 and 0.78 in the band (inside),
  implied Ea on item 29's pivot about 32 and 27 kJ/mol; the 44.8 kJ/mol case
  (factor 4.29 at 136 C) is 4.29 / rf^2 = 6.46 and 7.58 times the law at
  0.885 mm in D/R^2 terms, beyond item 29's constant-factor extension.
* **Reported without a threshold:** per cell the window residual (decimated
  "as printed" and full-step), points in the reading band, crossing time
  (closed form; identical at every cell), strides, ledgers; orders; the
  continuum estimate by item 27's two-factor construction and the numerical
  uncertainty of the reported level; the band edges and the interior minimum;
  every temperature alternative's factor at 136 C, position in the band,
  residual and in-band count; the population mixtures.

**The population at the species-own geometry, declared here.** Item 30's
five weightings are kept, with every class diameter multiplied by
d_thesis / 1.77 (1.1017 sunflower, 1.0169 soybean) before the psi = 0.74
scaling, so that (a1) and (a3) have their mass-mean at the species' own
1.95 and 1.80 mm and the relative spread (1.0 to 3.0 mm about 1.77 mm) is
item 30's. This is a declaration: the thesis prints no size distribution for
either sample, and the 1 to 3 mm statement is Faner 2019's industrial soybean
meal. For reference, item 30's unscaled mixtures (the absolute 1 to 3 mm
classes) are compared with the species-own single sphere as well, from item
30's own runs (no new march).

---

## 1. The declared radii, and what ran

**Derivation (declared, not fitted).** With sphericity psi = pi d_v^2 / A
(d_v the volume-equivalent diameter, A the particle surface), the sphere with
the particle's own volume-to-surface ratio has radius 3 V/A = psi d_v / 2
(item 30 section 1). Here d_v is each thesis sample's own printed
equivalent-sphere diameter, 1.95 mm (sunflower) and 1.80 mm (soybean), so

* **R_d,sun = 0.74 x 1.95 mm / 2 = 0.7215 mm** (radius factor 0.81525 on the
  frozen 0.885 mm sphere),
* **R_d,soy = 0.74 x 1.80 mm / 2 = 0.6660 mm** (radius factor 0.75254).

psi = 0.74 is the mean sphericity Faner et al. (2019) measured on industrial
soybean meal. **For the sunflower sample it is an assumption, declared**: no
held source prints a sphericity for the thesis's sunflower meal (nor, strictly,
for the thesis's soybean sample; 0.74 belongs to the 2019 meal). What changes
against item 30's R_d = 0.655 mm is the diameter only: each trace now uses its
own sample's printed diameter instead of the 2019 meal's 1.77 mm, so the
sunflower trace no longer borrows the soybean meal's diameter (referee audit 5,
M1). The radius factor is item 30's float expression for its `thesis_aris`
convention, so the reported level is item 30's conventions-table row (13.3 and
19.7 per cent) by exact identity. The film transfer area per kg dry is held at
the committed sample-holder value (0.534 and 0.638 m2/kg) for every radius and
class; the constant-rate leg and the crossing do not depend on the radius.

**What ran** (82 marches, one at a time; logs `log_*.txt`):

* Q0: **17,043 of 17,043 checks exact over 2 marches** (`identity.json`),
  item 30's `conv_thesis_aris_<key>`; these are the ladder's 60/1 cells.
* The 3 x 3 (60, 120, 240 cells x time level 1, 2, 4) at each species-own
  radius: 18 marches (identity plus 12 further cells run by the previous agent
  and reused, complete and deterministic; the four 240/2 and 240/4 cells run
  on completion).
* The 95 per cent edges at item 28's factors 3.1401 and 3.1323, 60/1: 4; the
  interior sweep at item 30's fractions (-0.75 to +0.75 of ln F): 18.
* Item 29's temperature alternatives at the species-own radii, 60/1: 20 (per
  trace Ea_law +/- 16.1 kJ/mol, the 44.81 kJ/mol case, the Ea implied by the
  identified constant transported to the species-own radius, constant factors
  1.5, 2.0, 2.5, 3.0, and the constant pair matched to the Ea pair's factor at
  136 C). The law itself (Ea at the source's centre) is the ladder's 60/1.
* D/R^2 similarity checks at 0.885 mm: 4.
* Population classes at the species-own scaling (section 0), 60/1: 18.

Ledgers over all 82 marches: every accepted stride closes to at most 1.3e-14
of the particle inventory, every whole march to at most 1.9e-12, every remap
to at most 1.1e-14; every stride's read equals the face datum. **Every march
covers the scored window.** Three refusals, all after it: the sunflower
quarter-step cells 60/4, 120/4 and 240/4 end at 2399.895 s, 0.10 s short of
the 2400 s horizon (verbatim `RuntimeError: dry-shell Newton line search
failed inside the admissible storage domain`, seven retries), the
horizon-remainder class of items 24, 27, 29 and 30, more than 2,000 s after the
last scored sample (300.1 s); their scored numbers are kept. One cell, soybean
240/1, recorded 3,672 s of wall time against 54.6 s for the sunflower 240/1
cell; its output is complete (484 strides, t_max, 0 retries) and the march is
deterministic, so the wall time is a host stall, not a solver event.

Residuals are the window RMS as a fraction of the window span; "as printed" is
the committed decimated scoring, "full-step" item 27's refinement metric; in
band = measured points (of 15) inside the reading band.

## 2. The refinement ladder at the species-own radii

| trace (radius) | cells / time level | as printed | in band | full-step | in band | crossing, s |
|---|---|---:|---:|---:|---:|---:|
| sunflower (0.7215 mm) | 60/1 (reported) | **0.1330** | **10** | 0.1325 | 11 | 51.5704 |
| sunflower | 60/2, 60/4 | 0.1320, 0.1315 | 11, 11 | 0.1318, 0.1315 | 11, 11 | same |
| sunflower | 120/1, 120/2, 120/4 | 0.1323, 0.1314, 0.1310 | 11 | 0.1320, 0.1313, 0.1310 | 11 | same |
| sunflower | 240/1, 240/2, 240/4 | 0.1321, 0.1313, 0.1309 | 11 | 0.1318, 0.1312, **0.1308** | 11 | same |
| soybean (0.6660 mm) | 60/1 (reported) | **0.1971** | **6** | 0.1931 | 6 | 55.2361 |
| soybean | 60/2, 60/4 | 0.1923, 0.1913 | 6, 7 | 0.1917, 0.1911 | 7, 7 | same |
| soybean | 120/1, 120/2, 120/4 | 0.1962, 0.1913, 0.1902 | 6, 7, 7 | 0.1920, 0.1907, 0.1900 | 6, 7, 7 | same |
| soybean | 240/1, 240/2, 240/4 | 0.1959, 0.1911, 0.1900 | 6, 7, 7 | 0.1918, 0.1904, **0.1898** | 7 | same |

The crossing time is one float at every cell (51.570405232184214 and
55.236094599373544 s, identical to items 24, 27, 29 and 30).

**Observed orders** (full-step RMS, item 24's rule): **mesh 1.98, 1.93, 1.91
(sunflower) and 2.01, 2.02, 2.04 (soybean)** at time levels 1, 2, 4; **time
1.10 to 1.11 (sunflower) and 1.16 to 1.19 (soybean)** at 60, 120, 240 cells;
diagonal 1.47 and 1.48. Every chain forms. (Decimated scoring: mesh 1.99 to
2.00 and 1.96 to 2.06; time 0.97 to 1.22 and 2.13 to 2.18, the soybean
decimated time chain carrying the output sampling.)

**Mesh-independent within the time-step error.** The largest 60-to-240 change
is 0.00070 (sunflower) and 0.00132 (soybean) of the span, below the
240/1-to-240/4 change, 0.00097 and 0.00200.

**Converged values** (item 27's two-factor construction from 240/4, order one
and observed orders):

| trace | reported, as printed | reported, full-step | 240/4 | **continuum** | **numerical uncertainty of the reported value** (as printed minus continuum) |
|---|---:|---:|---:|---:|---:|
| sunflower | 0.1330 (10/15) | 0.1325 (11/15) | 0.1308 | **0.1304 to 0.1305** | **+0.0024 to +0.0026** (full-step +0.0020 to +0.0021) |
| soybean | 0.1971 (6/15) | 0.1931 (6/15) | 0.1898 | **0.1889 to 0.1892** | **+0.0079 to +0.0082** (full-step +0.0039 to +0.0042; the rest output sampling) |

The mean signed residual stays positive in the continuum (+0.111 and +0.178 of
the span): at the species-own radii the march dries more slowly than both
traces. Time to the last measured loading, continuum 326.3 and 303.0 s after
the crossing (item 30 at R_d: 268.6 and 292.9 s).

Against item 30's R_d = 0.655 mm (the 2019 meal's diameter for both traces):
sunflower 13.3 per cent (10 of 15) against 7.9 (15 of 15), converged 13.0
against 7.8; soybean 19.7 (6 of 15) against 18.7 (7 of 15), converged 18.9
against 17.8. The diameter differs by 10 per cent for sunflower (1.95 against
1.77 mm) and by 2 per cent for soybean (1.80 against 1.77 mm).

## 3. The 95 per cent band at the species-own radii (60/1)

| trace | D/F (low edge) | minimum inside the band (marched) | parabola through the minimum | central line | D x F (high edge) | **95 % range** |
|---|---:|---|---|---:|---:|---|
| sunflower, F = 3.1401 | 0.404 (1/15) | 0.0508 (15/15) at D x F^0.3 = D x 1.41 | 0.050 at D x 1.45 (fraction 0.324) | 0.133 | 0.312 (0/15) | **at most 5.1 to 40.4 %** |
| soybean, F = 3.1323 | 0.503 (0/15) | 0.0896 (14/15) at D x F^0.4 = D x 1.58 | 0.089 at D x 1.56 (fraction 0.389) | 0.193 | 0.334 (5/15) | **at most 8.9 to 50.3 %** |

Full-step scoring, as items 27, 28 and 30 printed the band (decimated: 5.0 to
40.6 and 9.3 to 50.5 %; decimated minima 0.050 and 0.093 at D x 1.45 and
1.57). As at R_d, both edges lie on the far side of a minimum inside the band;
the fast edge D x F is too fast (mean signed -0.29 and -0.28 of the span). The
minima sit where section 0 predicted from item 29's radius minima through the
D/R^2 similarity: (0.8153/0.677)^2 = 1.45 and (0.7525/0.601)^2 = 1.57.

Sweep (full-step RMS, in band, by fraction of ln F): sunflower -0.75: 0.348
(2), -0.50: 0.283 (2), -0.25: 0.211 (3), +0.10: 0.100 (14), +0.20: 0.070
(15), +0.30: 0.051 (15), +0.40: 0.058 (15), +0.50: 0.087 (14), +0.75: 0.191
(7); soybean -0.75: 0.437 (0), -0.50: 0.363 (1), -0.25: 0.281 (2), +0.10:
0.158 (9), +0.20: 0.125 (12), +0.30: 0.099 (12), +0.40: 0.090 (14), +0.50:
0.104 (13), +0.75: 0.202 (8). Band width 0.35 and 0.41 of the span. For
information only: item 28's leverage construction on the species-own
reported-level steps gives F = 3.138 and 3.130 (within 0.1 per cent of item
28's), the boundary loading inside the measured range at every window stride
(lowest 1.034 and 1.047 times the measured minimum).

## 4. The temperature alternatives at the species-own radii (60/1)

Item 29's construction exactly (pivot 95 C, multiplier along the whole
temperature path; constant factors by `ScaledCase`). Position = ln(factor on
the law at 136 C) / ln F. RMS as printed (in band); full-step in
`temperature_table.csv`, within 0.005.

| alternative | Ea, kJ/mol | factor at 136 C | position | sunflower RMS (in band) | item 29 at 0.885 mm | soybean RMS (in band) | item 29 at 0.885 mm |
|---|---:|---:|---:|---:|---:|---:|---:|
| the law (Ea at the source's centre) | 0.31 | 1.000 | 0.00 | 0.133 (10) | 0.244 (2) | 0.197 (6) | 0.366 (0) |
| Ea + 16.1 | 16.41 | 1.694 | +0.46 | **0.100 (12)** | 0.198 (2) | **0.180 (7)** | 0.326 (0) |
| Ea - 16.1 | -15.79 | 0.590 | -0.46 | 0.177 (5) | 0.282 (3) | 0.224 (3) | 0.395 (0) |
| implied by the identified constant at 0.885 mm | 44.81 | 4.293 | +1.27 | 0.198 (7) | 0.157 (7) | 0.259 (4) | 0.268 (5) |
| implied by the identified constant transported to the species-own radius (1.165e-9 rf^2) | 32.34 (sun), 27.45 (soy) | 2.854 (sun), 2.431 (soy) | **+0.92, +0.78** | 0.131 (9) | | 0.192 (7) | |
| constant 1.5 | | 1.5 | +0.35 | **0.050 (15)** | 0.134 (10) | **0.095 (13)** | 0.247 (1) |
| constant 2.0 | | 2.0 | +0.61 | 0.127 (9) | 0.058 (15) | 0.139 (11) | 0.159 (9) |
| constant 2.5 | | 2.5 | +0.80 | 0.213 (6) | 0.068 (15) | 0.228 (7) | 0.103 (12) |
| constant 3.0 | | 3.0 | +0.96 | 0.290 (0) | 0.126 (9) | 0.311 (5) | 0.099 (13) |
| constant matched to Ea + 16.1 | | 1.694 | +0.46 | 0.073 (15) | | 0.098 (14) | |
| constant matched to Ea - 16.1 | | 0.590 | -0.46 | 0.273 (2) | | 0.355 (0) | |

**The D/R^2 identity (verified).** (R_own, law) against (0.885 mm, law /
rf^2): RMS equal to 1.2e-14 (sunflower) and 5.5e-14 (soybean) of the span,
same strides (483, 484), same in-band count, stride lengths equal to 2.7e-15 s,
mean-loading path equal to 8.3e-13 relative. (R_own, Ea + 16.1) against
(0.885 mm, Ea + 16.1 times 1/rf^2): 2.1e-16 and 5.8e-16, strides 481 and 482,
path 4.4e-13 and 3.2e-13. The Arrhenius multiplier does not depend on D's
magnitude, so the similarity carries through the temperature alternatives.

**Lever slopes** (secant over the matched pair, per unit change of ln(factor at
136 C), as printed): Ea pair -0.072 (sunflower) and -0.042 (soybean) (item 29
at 0.885 mm: -0.080 and -0.065); matched constant pair -0.190 and -0.243;
the radius by the similarity +0.38 and +0.49 per unit ln R. Ea over constant:
0.38 and 0.17. The constant pair's secant is flatter than item 29's local
slope because at the species-own radii its upper point (1.694) lies past the
residual's minimum on the sunflower trace and near it on the soybean trace.

**Which of item 29's statements hold at the species-own radii, and which
change:**

* **Holds: the march depends on R and D only through D/R^2** (above), so the
  radius is not an independent lever.
* **Holds: the ordering of the levers**, radius > constant factor on D >
  activation energy pivoted at 95 C, the Ea pair moving the residual in the
  same direction as the matched constant pair and by 0.38 (sunflower) and
  0.17 (soybean) of it per unit log factor (item 29: about a third to a
  quarter; soybean now about a sixth, a secant effect near the minimum).
* **Changes in magnitude: "the activation energy anywhere inside the source's
  95 per cent half-width does not reach the data".** At the species-own radii
  the +16.1 kJ/mol edge takes sunflower from 0.133 (10) to 0.100 (12 of 15 in
  band) and soybean from 0.197 (6) to 0.180 (7); the band minima are 0.050 and
  0.089. It closes 40 per cent (sunflower) and 16 per cent (soybean) of the gap
  between the law and the minimum and brings neither trace to 15 of 15, so the
  statement still holds in the sense of not reaching the minimum or the whole
  reading band; the range inside +/- 16.1 kJ/mol is 10.0 to 17.7 per cent
  (sunflower) and 18.0 to 22.4 per cent (soybean), not item 29's 19.8 to 28.2
  and 32.6 to 39.5.
* **Holds with new numbers: the constant factor that minimises the residual
  lies inside the 95 per cent band**: at D x 1.45 (sunflower) and 1.56
  (soybean), positions 0.32 and 0.39 (item 29 at 0.885 mm: 2.15 and 2.78,
  positions 0.67 and 0.90).
* **Changes: the identified constant outside the law's band (item 29's S6).**
  The 44.81 kJ/mol activation energy is a property of the identification at
  0.885 mm and stays outside (factor 4.29, position +1.27); marched at the
  species-own radii it overshoots (0.198 and 0.259, worse than the law's 0.133
  and 0.197). The same identification transported to the species-own radius by
  the exact similarity, 1.165e-9 x rf^2 = 7.74e-10 (sunflower) and
  6.60e-10 m2/s (soybean), is a factor 2.854 and 2.431 on the law at 136 C,
  positions **+0.92 and +0.78, inside the band**; the activation energy it
  implies on item 29's pivot is **32.3 and 27.4 kJ/mol**, still outside the
  +/- 16.1 kJ/mol half-width about the law's 0.31 kJ/mol. (The referee's
  arithmetic at R_d = 0.655 mm, factor 2.35 and about 26 kJ/mol, is the same
  construction at the other radius.) Marched as Arrhenius alternatives these
  two give 0.131 (9) and 0.192 (7), about the law's own residual; as a constant
  D they are, by the similarity, the paper's constant identification (15.8 and
  26.5 per cent at 0.885 mm).

## 5. The population at the species-own geometry (60/1, psi = 0.74 on every class)

Classes are item 30's declared diameters times d_thesis / 1.77 (1.1017
sunflower, 1.0169 soybean), then psi-scaled; the mixture is the
dry-mass-weighted mean loading, scored with the same window residual.
Self-checks: the species-mean class is the ladder's 60/1 march, and the
weight-one mixture reproduces the single sphere's full-step and decimated
scores exactly.

| distribution (declared, section 0) | sunflower mass-mean d, mm | RMS as printed (in band) | minus single | soybean mass-mean d, mm | RMS as printed (in band) | minus single |
|---|---:|---:|---:|---:|---:|---:|
| single sphere at the species-own radius | 1.950 | 0.1330 (10) | 0 | 1.800 | 0.1971 (6) | 0 |
| (a1) tails 20 %, mean held | 1.950 | 0.1025 (14) | **-0.0305** | 1.800 | 0.1641 (8) | **-0.0329** |
| (a2) literal 10 / 80 / 10 | 2.001 | 0.1171 (12) | -0.0159 | 1.847 | 0.1798 (7) | -0.0172 |
| (a3) fine tail 10 %, mean held | 1.950 | 0.1082 (12) | -0.0248 | 1.800 | 0.1703 (7) | -0.0268 |
| (b) five equal-mass bins, midpoints | 2.203 | 0.1575 (4) | **+0.0245** | 2.034 | 0.2240 (3) | **+0.0269** |
| (b') five equal-mass classes, endpoints | 2.203 | 0.1285 (11) | -0.0045 | 2.034 | 0.1922 (7) | -0.0049 |

(Full-step shifts within 0.0005 of these, same signs.) **Signs, per
distribution and species, from the numbers:** (a1), (a2), (a3) and (b') lower
the residual on both traces (by 0.030, 0.016, 0.025 and 0.005 sunflower;
0.033, 0.017, 0.027 and 0.005 soybean); (b) at the bin midpoints raises it on
both (+0.025, +0.027). The pattern is item 30's: at the species-own radii the
single sphere dries too slowly; the finest classes (the scaled 1.0 and 1.2 mm
classes, mean signed -0.29 and -0.15 sunflower, -0.28 and -0.12 soybean) dry
faster than the data and pull the mixture towards it; (b) at the midpoints has
its finest class at 1.2 mm rather than 1.0 mm and a mass-mean above the
species diameter (2.20 and 2.03 mm). The shifts, at most 0.033 of the span, are 4 to 13 times the
numerical uncertainty of the reported level and at most a quarter of the
digitization half-width (0.142, 0.157).

For reference, item 30's unscaled mixtures (absolute 1 to 3 mm classes, from
item 30's runs, no new march) against the species-own single sphere:
sunflower -0.078, -0.068, -0.075, -0.030, -0.058 for (a1), (a2), (a3), (b),
(b') (all lower); soybean -0.043,
-0.027, -0.037, **+0.017**, -0.015 ((b) at the midpoints higher, the rest
lower).

## 6. Verdict against the rejecting outcomes

* **Q0 — passed.** 17,043 of 17,043 checks exact over 2 marches.
* **R1 — does not fire.** Every mesh chain forms, orders 1.91 to 1.98
  (sunflower) and 2.01 to 2.04 (soybean), all above 1.5 (the sunflower
  60-to-240 changes, 0.00063 to 0.00070, exceed 0.0005, so the order test is
  the one that applies, and it passes); the largest 60-to-240 change is below
  the 240/1-to-240/4 change on both traces (0.00070 < 0.00097;
  0.00132 < 0.00200).
* **R2 — does not fire.** All 82 marches cover the scored window; three
  refusals at 2399.895 s, after it.
* **R3 — does not fire.** (i) The D/R^2 identity holds to at most 5.5e-14 of
  the span (threshold 1e-9); (ii) the Ea pair moves the residual by 0.38 and
  0.17 of the matched constant pair per unit log factor, below it; (iii) both
  pairs lower the residual as the factor rises, the same direction.

## 7. Sentences the paper could carry (not applied)

* "With each thesis sample's own equivalent diameter (1.95 mm sunflower,
  1.80 mm soybean) at the same sphericity, 0.74, the volume-to-surface sphere
  has radius 0.722 and 0.666 mm; the sphericity was measured on the 2019
  industrial soybean meal and is assumed for the sunflower sample."
* "At these radii the falling-rate window residual is 13.3 and 19.7 per cent of
  the window span (10 and 6 of 15 points in the reading band), converged 13.0
  and 18.9 per cent, second order in the mesh and first order in the step,
  numerical uncertainty 0.2 and 0.8 percentage points; the march dries more
  slowly than both traces."
* "For the sunflower trace, the 7.9 per cent at R_d = 0.655 mm takes both the
  sphericity and the diameter from the 2019 soybean meal; with the sunflower
  sample's own diameter it is 13.3 per cent. For soybean the two diameters
  agree to 2 per cent and the residuals are 18.7 and 19.7 per cent."
* "Carried at the law's 95 per cent prediction factor the residual ranges over
  at most 5.1 to 40.4 and 8.9 to 50.3 per cent, with its minimum inside the
  band at D x 1.45 and 1.56."
* "At these radii an activation energy inside the source's +/- 16.1 kJ/mol
  puts the residual between 10.0 and 17.7 per cent (sunflower) and 18.0 and
  22.4 per cent (soybean), short of the band minimum; the identified constant
  transported to these radii by the D/R^2 similarity is 2.9 and 2.4 times the
  law at 136 C, inside the law's 95 per cent band, and implies 32 and 27
  kJ/mol, outside the source's half-width. The 44.8 kJ/mol figure belongs to
  the volume-equivalent sphere."
* "A declared 1 to 3 mm population scaled to each sample's diameter moves the
  residual by -0.033 to +0.027 of the span: down whenever a 1.0 mm class
  carries 10 to 20 per cent of the mass, up only for the five-bin midpoint
  quadrature, whose finest class is 1.2 mm and whose mass-mean lies above the
  sample's diameter."

## 8. Claim and non-claim

Claimed: at the declared species-own volume-to-surface radii (0.7215 and
0.6660 mm; psi = 0.74 from Faner 2019's soybean meal, assumed for sunflower;
diameters from the thesis), nothing fitted, the Faner window residual with the
boundary-value read is second order in the mesh and first in time,
mesh-independent within the time-step error, converged to 0.130 and 0.189 of
the span with the reported level 0.0024 to 0.0026 and 0.0079 to 0.0082 above;
the 95 per cent band gives at most 5.1 to 40.4 and 8.9 to 50.3 per cent with
interior minima at D x 1.45 and 1.56; item 29's temperature alternatives and
the D/R^2 identity as in section 4; the population shifts of section 5.

Not claimed: that 0.74 is the sphericity of the sunflower sample or of the
thesis soybean sample; that either sample is a sphere of these radii or that
the volume-to-surface sphere is the correct kernel geometry (item 30 section
9's caveats apply unchanged: the equivalence is exact only at short times,
the window reaching a Fourier number of about 0.3, and the convention family
was proposed after item 29 showed that smaller spheres fit better); any size
distribution for the thesis samples (the declared distributions scale a
one-line statement about the 2019 meal); a preferred temperature dependence
(the transported identification is a consequence of the similarity, not an
activation-energy measurement); anything about the space-resolved solve, the
other film cases or the journal conditions; calibration, validation,
qualification or plant prediction.

## 9. Adversarial paragraph

The species-own radius removes the borrowed diameter from the sunflower trace
but keeps a borrowed sphericity on both: psi = 0.74 was measured on a 2019
industrial soybean meal, not on either thesis sample, and sunflower meal
(hull fragments, a different flake form) need not share it. Because the
residual depends on psi only through D/R^2, a sunflower sphericity of 0.67
would reproduce item 30's 7.9 per cent at the thesis diameter
(0.67 x 0.975 mm = 0.655 mm), and psi = 1 gives the thesis volume-equivalent
sphere's 29.1 per cent; the comparison cannot tell these apart, and item 29
section 5's measurement (a printed flake thickness or sphericity for the
thesis samples) is still the one that would settle it. The mixing of sources
is also asymmetric: R_d takes the 2019 meal's psi and d for both traces, the
species-own radius combines the 2019 psi with the thesis d; for the sunflower
sample neither is a single-source geometry. Neither convention is preferred
by this measurement; it states what each gives.
