# Item 29 — the Faner march's sensitivity to the particle radius and to the temperature dependence of D

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

## 0. Rejecting outcomes, written 2026-09-26 before any march of this item ran

Brief: a sensitivity measurement on declared alternatives, both thesis traces
(sunflower and soybean, 136 C), item 27's reported level (60 cells, reported
step, boundary-value read). Two levers: the sphere radius (0.5 to 1.2 times
0.885 mm, film area per kg dry held fixed) and the temperature dependence of
the law (Arrhenius alternatives pivoted at 95 C, and a plain constant factor as
the reference scale), and their cross grid.

* **Q0 (identity first).** The reported level through this item's driver
  (radius passed explicitly as 0.885e-3 m, the Arrhenius wrapper at the law's
  own activation energy) reproduces item 27's two thesis runs
  `refine_<key>_n060_t1` by exact float identity, every summary leaf and every
  recorded step value. A single mismatch stops the script.
* **S1 (the radius lever does not bracket the data).** Fires on a trace if the
  window RMS residual is monotone over the radius factors 0.5 to 1.2, so that
  no crossing radius can be interpolated inside the declared ladder. A
  declared extension (smaller factors) is then run and marked as such.
* **S2 (the diffusion-time scaling fails).** If the march depended on R and D
  only through R^2/D, the residual at radius factor f_R and D factor 1 would
  equal the residual at radius 1 and D factor f_R^-2. Fires if, at any radius
  factor where both are measured (the D-axis interpolated linearly in ln f_D),
  they differ by more than 0.01 of the window span. The expected breakers are
  named in advance: the particle temperature path runs in real time with a
  radius-independent film heat input per kg, the stride rule's absolute floor
  (2 ms) and ceiling (5 s), and the output decimation.
* **S3 (a lever cannot be ranked).** The two levers are compared as the change
  of the window residual per unit change of ln R and of ln D at the reported
  point (central differences on the cross grid). Fires if the two local slopes
  change sign or rank between the one-step and the two-step differences, so
  that "stronger per unit log change" is not a local statement.
* **S4 (a march refuses inside the window).** Any march that does not cover the
  last scored sample is reported verbatim and left unscored; nothing is relaxed
  (floor, tolerance, retry budget, horizon) to obtain a number.
* **S5 (the hand value of the implied activation energy does not hold).** The
  brief's hand value, about 32.6 kJ/mol, for the activation energy that takes
  the law at 95 C to the identified constant 1.165e-9 m2/s at 136 C. Fires if
  the exact computation, with the law evaluated at the loading the march
  reads, differs from it by more than 1 kJ/mol.
* **S6 (the identified constant is outside the law's 95 per cent band).** Fires
  if the factor on the law at 136 C that alternative (b) implies exceeds the
  95 per cent prediction factor at the boundary read (3.140 sunflower, 3.132
  soybean, item 28).

Reported without a threshold: per march the window residual (as printed, the
committed decimated scoring; and full-step), the points inside the reading
band, the stride counts (to the last scored sample, inside the window, whole
march to the 2400 s horizon), the ledgers.

---

## 1. What ran

* **The march** is item 27's instrumented march at the reported level (60
  cells, reported step, boundary-value read, `lit_whitaker` film, `sensible`
  temperature case), imported read only. Two inputs are varied: `radius_m`,
  which the committed `march_tail` accepts and item 27 never passed, and a
  multiplier on the carried law. **The film heat-transfer area per kg dry is
  held fixed** at the committed 0.534 (sunflower) and 0.638 (soybean) m2/kg:
  it is the sample-holder area over the dry mass, not a particle surface, and
  the constant-rate leg (closed form, crossing at 51.570 and 55.236 s) does not
  depend on the radius.
* **Q0: 16,520 of 16,520 checks exact** over the two reported-level runs
  (radius passed explicitly as 0.885e-3 m, the Arrhenius wrapper at the law's
  own activation energy): every summary leaf and every recorded step value of
  item 27's `refine_<key>_n060_t1`.
* **62 marches** one at a time, 9 to 16 s each: radius ladder 0.5, 0.6, 0.7,
  0.85, 1.0, 1.2 (and 0.55, 0.65, 0.75, a declared refinement near the
  minimum); four Arrhenius alternatives; constant factors 1.5 to 3.0 (and 3.5,
  4.0, 5.0, a declared extension because the soybean minimum lay at the 3.0
  end); the cross grid radius {0.7, 0.85, 1.0, 1.2} x factor {1/1.5, 1, 1.5, 2}.
  Every accepted stride closes to at most 1.2e-14 of the particle inventory,
  every march to at most 1.5e-12. Every march covers the scored window; four
  refuse 0.18 s short of the 2400 s horizon (the law +16.1 kJ/mol and the
  factor 1/1.5, both traces; verbatim `RuntimeError: dry-shell Newton line
  search failed inside the admissible storage domain`, seven retries), the
  horizon-remainder class items 24 and 27 recorded, more than 2,000 s after
  the last scored sample.

Residuals are the window RMS as a fraction of the window span, "as printed"
(the committed decimated scoring behind the paper's 24.4 and 36.6); in band is
out of 15 points. Full-step values are in the CSV tables and move nothing
below by more than 0.005.

## 2. The radius

| radius factor | radius, mm | sunflower RMS | in band | soybean RMS | in band | strides, crossing to last scored sample (whole march) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.443 | 0.238 | 5 | 0.186 | 9 | 58 / 58 (490) |
| 0.55* | 0.487 | 0.162 | 8 | 0.119 | 12 | 57 / 56 (488) |
| 0.60 | 0.531 | 0.099 | 12 | **0.093** | **14** | 56 / 55 (487) |
| 0.65* | 0.575 | 0.057 | 15 | 0.114 | 12 | 55 / 54 (486) |
| 0.70 | 0.620 | **0.055** | **15** | 0.153 | 9 | 54 / 53 (485) |
| 0.75* | 0.664 | 0.086 | 14 | 0.195 | 7 | 53 / 52 (484) |
| 0.85 | 0.752 | 0.156 | 4 | 0.272 | 1 | 52 / 51 (483) |
| **1.00** | **0.885** | **0.244** | **2** | **0.366** | **0** | 50 / 49 (481) |
| 1.20 | 1.062 | 0.330 | 2 | 0.458 | 0 | 49 / 49 (481) |

\* declared refinement near the minimum. Stride counts sunflower / soybean.
The first stride scales with R^2 (0.80 s at 0.5 against 3.21 s at 1.0); the
5 s ceiling binds at every radius (2e-3 tau is at least 16 s), so the count
grows only from 481 to 490 over the whole march and from 50 to 58 up to the
last scored sample.

**Crossing radii.** The residual minimum by a parabola in ln R through the
three lowest points: **0.677 x 0.885 = 0.599 mm (sunflower)**, RMS 0.051, and
**0.601 x 0.885 = 0.532 mm (soybean)**, RMS 0.093; the zero of the mean
signed residual at 0.675 (0.597 mm) and 0.583 (0.516 mm). As diameters 1.20
and 1.06 mm, against the thesis's printed equivalent-sphere diameters 1.95
and 1.80 mm (the march's 1.77 mm is the frozen radius of the manuscript, not
the thesis's): 0.61 and 0.59 of the thesis sphere. At the thesis's own
equivalent-sphere radius (factors 1.102 and 1.017) the law gives 0.290 and
0.374 (interpolated), worse than the reported 0.244 and 0.366.

**The march depends on R and D only through D/R^2.** With the film area per
kg held fixed, the dry-shell equation, its Dirichlet datum, the loading-based
read and the per-kg heat balance are invariant under R -> kR, D -> k^2 D, and
the stride rule (tau = R^2 cap / D) maps onto itself. Measured: radius 0.5
against D x 4 at radius 1.0 gives the same residual to 1e-16 (sunflower) and
3e-15 (soybean), the same 490 strides, the same stride lengths to 1.0e-15 and
the same mean-loading path to 3.5e-13 relative at every stride. Across the
whole grid, every radius/factor cell against the factor axis interpolated at
the same D/R^2: at most 0.0093 of the span apart, all of it the linear
interpolation's error on a curved function (largest next to the minimum).
S2 does not fire. The radius is therefore not an independent lever: the
residual is a function of ln D - 2 ln R.

## 3. The temperature dependence of D

Pivot 95 C (368.15 K): the top of the fitted range, the data nearest the
136 C application, so every alternative agrees with the law where its data
end and differs only in the 41 K extrapolation. (The fit's centroid in 1/T
would be the covariance-minimising pivot; 95 C is chosen so that no
alternative changes the law at the hottest measured temperature.) Each
alternative multiplies the law by exp(-(Ea_alt - Ea_law)/Rg (1/T - 1/368.15))
along the whole temperature path, from 68.7 C at the crossing (a factor below
1 when Ea_alt exceeds the law's 0.306 kJ/mol) to about 134 C at the window's
end.

**(b), the activation energy implied by the identified constant, computed
exactly** (`implied_ea.json`): at the loading of the 136 C boundary datum,
7.607e-4 kg/kg, the law gives 2.713e-10 m2/s at 136 C and 2.686e-10 at 95 C;
taking it to 1.165e-9 at 136 C needs **Ea = 44.81 kJ/mol**. The hand value
**32.65 kJ/mol is reproduced by its own construction**, Rg ln(1.165e-9 /
4.0e-10) / (1/368.15 - 1/409.15), but that construction takes the Table S2
reference value 4.0e-10 as the law at 95 C; the law at the loading the march
reads is 2.69e-10. S5 fires. The number is not a property of the data alone:
read at the 95 C boundary loading instead (3.34e-3 kg/kg, where the law gives
1.110e-9), the implied Ea is 1.5 kJ/mol. The exact value and the hand
construction were both marched.

| alternative | Ea, kJ/mol | factor on law at 136 C | at 68.7 C | position in 95 % band | sunflower RMS (in band) | soybean RMS (in band) |
|---|---:|---:|---:|---:|---:|---:|
| the law (Ea at the source's centre) | 0.31 | 1.000 | 1.000 | 0.00 | 0.244 (2) | 0.366 (0) |
| (a) centre + 16.1 | 16.41 | 1.694 | 0.667 | +0.46 | 0.198 (2) | 0.326 (0) |
| (a) centre - 16.1 | -15.79 | 0.590 | 1.498 | -0.46 | 0.282 (3) | 0.395 (0) |
| (b) implied, exact | 44.81 | **4.294** | 0.327 | **+1.27** | 0.157 (7) | 0.268 (5) |
| (b) hand construction | 32.65 | 2.883 | 0.444 | +0.93 | 0.158 (8) | 0.286 (3) |
| (c) constant 1.5 | | 1.5 | 1.5 | +0.35 | 0.134 (10) | 0.247 (1) |
| (c) constant 2.0 | | 2.0 | 2.0 | +0.61 | **0.058 (15)** | 0.159 (9) |
| (c) constant 2.5 | | 2.5 | 2.5 | +0.80 | 0.068 (15) | 0.103 (12) |
| (c) constant 3.0 | | 3.0 | 3.0 | +0.96 | 0.126 (9) | **0.099 (13)** |
| constant 3.5 (ext.) | | 3.5 | 3.5 | +1.10 | 0.184 (7) | 0.136 (11) |
| constant 4.0 (ext.) | | 4.0 | 4.0 | +1.21 | 0.238 (5) | 0.186 (9) |
| constant 5.0 (ext.) | | 5.0 | 5.0 | +1.41 | 0.335 (0) | 0.284 (6) |

Position = ln(factor at 136 C) / ln F with F the item 28 factor at the
boundary read, 3.1401 (sunflower) and 3.1323 (soybean); inside the band is
|position| <= 1 (positions for soybean differ from sunflower's in the third
decimal). **S6 fires**: the exact implied Ea puts the law at 4.29 times itself
at 136 C, outside the 95 per cent band (+1.27). Its residuals, 15.7 and 26.8
per cent, sit next to the paper's constant identification (15.8 and 26.5), as
they should, since both carry about 1.13e-9 m2/s through the window
(geometric mean over the marched strides to the last scored sample).

**Findings.** (i) Within the source's own 95 per cent interval on Ea the
temperature lever puts the residual between 0.198 and 0.282 (sunflower) and
0.326 and 0.395 (soybean); it does not reach the data. (ii) An activation
energy buys less than a constant factor of the same size at 136 C, because
with the pivot at 95 C it lowers D over the cooler early window: per unit
ln(factor at 136 C) the +/-16.1 kJ/mol pair moves the residual by -0.080
(sunflower) and -0.065 (soybean), a third to a quarter of a constant factor's
-0.25 and -0.27; per kJ/mol, -0.0026 and -0.0021 of the span. (iii) The
constant factors that minimise the residual, **2.15 (sunflower) and 2.78
(soybean)** by parabola in ln f (mean-signed zeros 2.20 and 2.95), lie
**inside the 95 per cent band** (positions 0.67 and 0.90); the law's own
prediction band does not discriminate them.

## 4. The two levers compared (cross grid)

Residual change per unit change of the logarithm at the reported point,
three-point differences on the grid (one step: 0.85 and 1.2 in R, 1/1.5 and
1.5 in D; two steps: 0.7 and 1.2 in R, 1/1.5 and 2 in D):

| trace | per ln R (one / two step) | per ln D (one / two step) | ratio R / D |
|---|---:|---:|---:|
| sunflower RMS | +0.508 / +0.493 | -0.253 / -0.247 | -2.005 / -1.994 |
| soybean RMS | +0.546 / +0.537 | -0.272 / -0.269 | -2.005 / -1.999 |
| sunflower mean signed | +0.490 / +0.491 | -0.245 / -0.246 | -2.002 / -2.000 |
| soybean mean signed | +0.530 / +0.528 | -0.265 / -0.264 | -2.003 / -1.999 |

**Per unit log change the radius is twice as strong as D, with the opposite
sign, and exactly so**: the ratio is -2.00 to within 0.6 per cent, which is
the D/R^2 similarity of section 2, not an empirical coincidence. S3 does not
fire. Per unit log change of D at 136 C, the constant factor is 3 to 4 times
stronger than an activation-energy change pivoted at 95 C. Ranking per unit
log change: radius (0.49 to 0.55 of the span) > constant D (0.25 to 0.27) >
Ea pivoted at 95 C (0.065 to 0.080).

## 5. What one thesis page could settle, and what it could not

The held thesis pages print only the equivalent-sphere diameters (1.95 mm
sunflower, 1.80 mm soybean) and the holder-thickness bound: no flake
thickness and no size distribution. A smaller sphere here is a proxy for a
shorter diffusion path, a flake's half-thickness, which the spherical kernel
cannot represent directly.

* **A printed flake thickness would fix the geometric half of D/R^2** from a
  measurement independent of the drying curves, so that the law's magnitude
  is tested at a stated geometry with no free lever left. The crossing spheres
  (0.60 and 0.53 mm) correspond to slab half-thicknesses of about R/3 to R/2
  (equal surface-to-volume ratio, or equal slowest decay mode), 0.20 to 0.30
  and 0.18 to 0.27 mm, flakes 0.4 to 0.6 and 0.35 to 0.53 mm thick (derived
  here, not measured). A printed thickness well above that range would say
  the law is too slow at 136 C by about the factor 2.2 to 2.8 of section 3;
  one inside it would say the gap is geometric.
* **It could not settle** the split cleanly with this kernel: the sphere to
  slab mapping is uncertain by R/3 against R/2, a factor 2.25 in D/R^2, as
  large as the gap itself, and a slab kernel would be needed to remove it. Nor
  could it settle the temperature dependence (Ea is bounded by the source,
  not measured; the 41 K extrapolation stays), the evaluation convention (the
  space-resolved solve's 0.70 and 0.77 of the span are untouched), or the
  loading dependence.
* **A size distribution** would turn the single march into a mass-weighted
  sum of marches over D/R^2 and test the tail's shape (fines speed the early
  window, coarse particles slow the tail). It could not separate geometry
  from D either, since every member depends on D/R^2 alone.

## 6. Verdict against the rejecting outcomes

* **Q0 — passed.** 16,520 checks exact.
* **S1 — does not fire.** Both minima lie inside the declared 0.5 to 1.2
  ladder (sunflower at 0.68, soybean at 0.60); the 0.55, 0.65, 0.75 rungs are
  a refinement, not an extension. (The constant-factor axis needed its
  declared extension to 3.5 to 5.0 to bracket the soybean minimum.)
* **S2 — does not fire.** Exact D/R^2 similarity: the residual at radius 0.5
  equals that at D x 4 to 3e-15, the loading path to 3.5e-13; interpolated
  cross-grid differences at most 0.0093 (interpolation error).
* **S3 — does not fire.** Slopes stable between one and two steps; ratio
  -2.00.
* **S4 — does not fire.** Every march covers the window; four
  horizon-remainder refusals 0.18 s before 2400 s, recorded.
* **S5 — fires.** Exact implied Ea 44.81 kJ/mol against the hand 32.65; the
  hand value used the 4.0e-10 reference value as the law at 95 C.
* **S6 — fires.** The identified constant sits at 4.29 times the law at 136 C
  at the boundary loading, position +1.27, outside the 95 per cent band.

## 7. Claim and non-claim

Claimed: at the reported level with the boundary-value read, the window
residual is a function of D/R^2 alone (exactly, with the film area per kg
held fixed), so the radius is the same lever as D with twice the log
sensitivity; the residual is smallest at a sphere 0.68 (sunflower) and 0.60
(soybean) times the frozen 0.885 mm, or equivalently at 2.15 and 2.78 times
the law, both inside the law's 95 per cent band; the source's own interval on
Ea (+/-16.1 kJ/mol, pivoted at 95 C) moves the residual by at most 0.05 of the
span and does not reach the data; the identified constant implies Ea = 44.8
kJ/mol at the marched loading (not 32.6), outside the source's interval and
putting the law outside its band.

Not claimed: a fitted radius, factor or activation energy (all are declared
alternatives; the minima are reported, not adopted); that the particles are
0.6 mm spheres or 0.4 to 0.6 mm flakes (the kernel cannot represent a slab
and the thesis prints no thickness); anything about the space-resolved solve,
the other film cases, the journal conditions or finer discretizations (item
27 puts the reported level within 0.4 and 0.8 points of its continuum value);
a better or validated agreement with the data; physical qualification or
plant prediction.
