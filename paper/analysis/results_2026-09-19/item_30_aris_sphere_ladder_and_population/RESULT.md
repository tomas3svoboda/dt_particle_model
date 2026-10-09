# Item 30 — the Faner march at the volume-to-surface equivalent sphere: refinement ladder, band, population

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

## 0. Rejecting outcomes, written 2026-09-26 before any march of this item except the identity check

Brief: the Faner falling-rate march (both thesis traces, sunflower and soybean
at 136 C; film case `lit_whitaker`, temperature case `sensible`, the carried
law read at the boundary value, item 27) at the sphere with the particle's
own volume-to-surface ratio, R_d = psi R_eq = 0.74 x 0.885 mm = 0.6549 mm,
derived from the source's measured sphericity (section 1), not fitted. The
film transfer area per kg dry is held fixed at the committed sample-holder
value for every radius and every population class.

* **Q0 (identity first).** Item 29's reported-level cells at radius factors
  0.70 and 0.75 (both traces) are reproduced through this item's driver by
  exact float identity, every summary leaf except tag and wall time and every
  recorded step value. A single mismatch stops the script. *Run before this
  section was written, as the brief orders: 33,330 of 33,330 checks exact over
  4 marches (`identity.json`).*
* **R1 (the residual at R_d is not converged in the mesh within the time-step
  error).** On each trace and each time level r in {1, 2, 4}, the mesh chain
  60/120/240 of the window RMS over the span (full-step scoring, item 27's
  metric). R1 fires if any chain either does not form or forms with an
  observed order below 1.5 while its 60-to-240 change exceeds 0.0005 of the
  span (item 27's threshold); or if the largest 60-to-240 change over the
  three time levels exceeds the time-step change at 240 cells (240/1 to
  240/4), i.e. the residual moves across the ladder by more than the
  time-step error.
* **R2 (the population moves the residual outside the single sphere's range
  by more than the digitization band).** The single sphere's range at R_d is
  the interval spanned by its nine ladder cells and its two continuum
  estimates (order one, observed orders), and, separately reported, its
  95 per cent band range. The digitization band is the traces' own reading
  band: half-width 0.01246 kg/kg on every scored point, 0.142 (sunflower) and
  0.157 (soybean) of the window span. R2 fires if, for any declared
  distribution and either trace, the population's window RMS (60/1,
  full-step and decimated) lies outside the ladder range by more than that
  half-width.
* **R3 (a level refuses inside the scored window).** Fires if any march of
  this item (ladder, band edges, population classes, conventions) ends before
  the last scored sample (300.1 s). Refusals after the window are recorded
  verbatim and do not fire R3. Nothing (floor, tolerance, retry budget,
  horizon, level) is relaxed to obtain a number.
* **Reported without a threshold:** per cell the window residual (decimated
  and full-step), points inside the reading band, the crossing time (closed
  form; it must be identical at every cell), stride counts and ledgers; the
  observed orders in mesh and time; the continuum estimate by item 27's
  two-factor construction and the numerical uncertainty of the reported-level
  value (reported level minus the continuum interval); the prediction factor
  that item 28's construction would give on the R_d reported-level steps (for
  information; the band edges use item 28's factors as instructed).

A note on the population declaration, written with the rejecting outcomes:
distribution (a) as briefed (three classes at 1.0, 1.77 and 3.0 mm, mass-mean
1.77 mm, the 1.0 and 3.0 mm classes 10 per cent of the mass each) is
over-determined: with 10, 80 and 10 per cent the mass-mean is 1.816 mm, and
with the mass-mean held at 1.77 mm the tails cannot both carry 10 per cent
(0.77 w_1 = 1.23 w_3). Three weightings of the same three marched classes are
therefore reported, all declared here: (a1) mass-mean 1.77 mm with the two
tails together 20 per cent (w = 0.1232, 0.8000, 0.0768), the primary; (a2) the
literal 10/80/10 (mass-mean 1.816 mm); (a3) the 1.0 mm tail at 10 per cent and
the mean held (w = 0.1000, 0.8374, 0.0626). Distribution (b), five classes
uniform in mass on 1 to 3 mm, is quadratured at the bin midpoints 1.2, 1.6,
2.0, 2.4, 2.8 mm (primary) and, as a declared variant, at 1.0, 1.5, 2.0, 2.5,
3.0 mm with equal weights; its mass-mean is 2.0 mm, not 1.77, which is a
property of the brief's declaration and is stated with the result. Every
class diameter is multiplied by psi = 0.74. The mixture is the dry-mass-weighted
mean loading, X_pop(t) = sum w_i X_i(t), scored by the same window residual.

---

## 1. The derived sphere, and what ran

**Derivation (source, not fit).** Faner et al. (2019), as resolved in
`docs/GT_PS2_PACKET_A2A_SOURCE_RESOLUTION.md` section 1: soybean meal from an
industrial extractor, mean equivalent diameter 1.77 mm, more than 95 per cent
of particles between 1 and 3 mm, mean sphericity 0.74. Frozen ruling PHY-034
carries one sphere of R_eq = 0.885 mm and names the 1 to 3 mm population as a
mandatory model-form sensitivity. With sphericity psi = pi d_v^2 / A (d_v the
volume-equivalent diameter, A the particle surface), the particle's
volume-to-surface ratio is V/A = psi d_v / 6, and the sphere with the same V/A
(Aris's characteristic length) has radius 3 V/A = psi R_eq:
**R_d = 0.74 x 0.885 mm = 0.6549 mm** (radius factor 0.74 on the frozen
sphere, passed to the march exactly as item 29 passes its factors).

**What ran** (68 marches, one at a time, about 25 minutes; logs `log_*.txt`):

* Q0: item 29's 0.70 and 0.75 reported-level cells, **33,330 of 33,330 checks
  exact** over 4 marches; the only input changed below is the radius.
* The 3 x 3 at R_d (60, 120, 240 cells x time level 1, 2, 4), both thesis
  traces: 18 marches.
* The 95 per cent edges at R_d, 60/1, at item 28's factors 3.1401 (sunflower)
  and 3.1323 (soybean), `ScaledCase(ln_shift = -/+ ln F)` exactly as items 27
  and 28 used them: 4 marches; and a declared interior sweep (fractions -0.75
  to +0.75 of ln F): 18 marches, added after section 0 and before the edges
  were read, because item 29's radius minima map (by the exact D/R^2
  similarity) to multipliers 1.2 and 1.5 at R_d, inside the band.
* Population classes at psi x {1.0, 1.2, 1.5, 1.6, 1.77, 2.0, 2.4, 2.5, 2.8, 3.0}
  mm, 60/1: 20 marches. **The film transfer area per kg dry is held fixed at
  the committed sample-holder value (0.534 and 0.638 m2/kg) for every class**,
  so every class shares the closed-form constant-rate leg and the crossing
  time; only the falling-rate tail differs by class.
* The source's (thesis) equivalent spheres, 1.95 mm (sunflower) and 1.80 mm
  (soybean), and their psi-scaled spheres, 60/1: 4 marches.

Ledgers over all 68 marches: every accepted stride closes to at most 1.0e-14
of the particle inventory, every whole march to at most 1.8e-12, every remap
to at most 8.8e-15. Every stride's read equals the face datum. **Every march
covers the scored window.** One refusal, after it: the band sweep at -0.50 ln F
(sunflower) refuses 0.19 s short of the 2400 s horizon (verbatim
`RuntimeError: dry-shell Newton line search failed inside the admissible
storage domain`, seven retries), the horizon-remainder class of items 24, 27
and 29, more than 2,000 s after the last scored sample; its scored numbers are
kept.

Residuals are the window RMS as a fraction of the window span; "as printed"
is the committed decimated scoring behind the paper's 24.4 and 36.6 per cent,
"full-step" is item 27's refinement metric. In band = measured points (of 15)
inside the reading band.

## 2. The refinement ladder at R_d

| trace | cells / time level | as printed | in band | full-step | in band | crossing, s |
|---|---|---:|---:|---:|---:|---:|
| sunflower | 60/1 (reported) | **0.0791** | **15** | 0.0791 | 15 | 51.5704 |
| sunflower | 60/2, 60/4 | 0.0786, 0.0784 | 15, 15 | 0.0786, 0.0784 | 15, 15 | same |
| sunflower | 120/1, 120/2, 120/4 | 0.0786, 0.0783, 0.0781 | 15 | 0.0788, 0.0783, 0.0781 | 15 | same |
| sunflower | 240/1, 240/2, 240/4 | 0.0785, 0.0783, 0.0781 | 15 | 0.0787, 0.0782, **0.0780** | 15 | same |
| soybean | 60/1 (reported) | **0.1868** | **7** | 0.1826 | 7 | 55.2361 |
| soybean | 60/2, 60/4 | 0.1818, 0.1807 | 7, 7 | 0.1812, 0.1806 | 7, 7 | same |
| soybean | 120/1, 120/2, 120/4 | 0.1859, 0.1808, 0.1797 | 7 | 0.1815, 0.1802, 0.1795 | 7 | same |
| soybean | 240/1, 240/2, 240/4 | 0.1856, 0.1806, 0.1794 | 7 | 0.1813, 0.1799, **0.1793** | 7 | same |

The crossing time is one float at every cell (51.570405232184214 and
55.236094599373544 s, identical to items 24, 27 and 29): closed form, no
discretization error, radius independent.

**Observed orders** (full-step RMS, item 24's order rule): **mesh 2.00, 2.05,
2.13 (sunflower) and 2.01, 2.03, 2.05 (soybean)** at time levels 1, 2, 4;
**time 1.29 to 1.33 (sunflower) and 1.10 to 1.14 (soybean)** at 60, 120, 240
cells; diagonal 1.56 and 1.44. Every chain forms.

**Mesh-independent to within the time-step error.** The largest 60-to-240
change is 0.00037 (sunflower) and 0.00130 (soybean) of the span, below the
240/1-to-240/4 change, 0.00070 and 0.00199.

**Converged values** (item 27's two-factor construction from 240/4; order one
and observed orders):

| trace | reported, as printed | reported, full-step | 240/4 | **continuum** | **numerical uncertainty of the reported value** (as printed minus continuum) |
|---|---:|---:|---:|---:|---:|
| sunflower | 0.0791 (15/15) | 0.0791 | 0.0780 | **0.0777 to 0.0779** | **+0.0013 to +0.0014** (full-step +0.0012 to +0.0013) |
| soybean | 0.1868 (7/15) | 0.1826 | 0.1793 | **0.1784 to 0.1786** | **+0.0082 to +0.0084** (full-step +0.0039 to +0.0042; the rest output sampling) |

The mean signed residual stays positive in the continuum (+0.053 and +0.167 of
the span): at R_d the march still dries more slowly than both traces, much
less so than at 0.885 mm (+0.22 and +0.35). Time to the last measured loading,
continuum 268.6 and 292.9 s after the crossing (0.885 mm: 493 and 537 s).

## 3. The 95 per cent band at R_d (60/1)

| trace | D/F (low edge) | minimum inside the band (marched) | parabola through the minimum | central line | D x F (high edge) | **95 % range** |
|---|---:|---|---|---:|---:|---|
| sunflower, F = 3.1401 | 0.367 (2/15) | 0.0524 (15/15) at D x F^0.2 = D x 1.26 | 0.050 at D x 1.20 | 0.079 | 0.400 (0/15) | **at most 5.2 to 40.0 %** |
| soybean, F = 3.1323 | 0.496 (0/15) | 0.0913 (15/15) at D x F^0.4 = D x 1.58 | 0.090 at D x 1.51 | 0.183 | 0.350 (4/15) | **at most 9.1 to 49.6 %** |

Full-step scoring, as items 27 and 28 printed the band (decimated: 5.2 to
39.8 and 9.4 to 49.8 %). The residual is not monotone across the band: at R_d
both edges lie on the far side of a minimum, and the fast-drying edge D x F is
now **too fast** (mean signed -0.38 and -0.30 of the span), unlike at 0.885 mm
where the data sat at the fast edge. The minima at multipliers 1.20 and 1.51
are item 29's radius minima (0.68 and 0.60 of R0) seen through the D/R^2
similarity: 0.74/sqrt(1.20) = 0.676 and 0.74/sqrt(1.51) = 0.602. The band width,
0.35 and 0.40 of the span, is below the span: item 24's O8 does not fire.
For information only: item 28's leverage construction applied to the R_d
reported-level steps gives F = 3.137 (sunflower) and 3.129 (soybean), within
0.1 per cent of item 28's factors, the boundary loading inside the measured
range at every window stride (lowest 1.036 and 1.048 times the measured
minimum).

## 4. The population at the reported level (psi = 0.74 on every class)

The mixture is the dry-mass-weighted mean loading of the separately marched
classes, scored with the same window residual. Self-checks: the 1.77 mm class
is the ladder's 60/1 march (same float radius, same scores), and the
weight-one mixture reproduces the single sphere's full-step and decimated
scores exactly.

| distribution (declared, section 0) | mass-mean d, mm | sunflower RMS as printed (in band) | minus single | soybean RMS as printed (in band) | minus single |
|---|---:|---:|---:|---:|---:|
| single sphere at R_d | 1.77 | 0.0791 (15) | 0 | 0.1868 (7) | 0 |
| (a1) 1.0 / 1.77 / 3.0, w 0.123 / 0.800 / 0.077 | 1.77 | 0.0546 (15) | **-0.025** | 0.1542 (9) | **-0.033** |
| (a2) literal 10 / 80 / 10 | 1.816 | 0.0653 (15) | -0.014 | 0.1697 (7) | -0.017 |
| (a3) fine tail 10 %, mean held | 1.77 | 0.0585 (15) | -0.021 | 0.1602 (9) | -0.027 |
| (b) five equal-mass bins, midpoints 1.2 ... 2.8 | 2.0 | 0.1030 (14) | **+0.024** | 0.2136 (4) | **+0.027** |
| (b') five equal-mass classes, 1.0 ... 3.0 | 2.0 | 0.0752 (15) | -0.004 | 0.1819 (7) | -0.005 |

(Full-step values within 0.005 of these, same signs; `population_table.csv`.)

* **Distribution (a) moves the residual more** than (b') and about as much as
  (b), in the opposite direction: **down** by 0.025 and 0.033 of the span,
  because the 1 mm fines (psi-scaled radius 0.37 mm) dry fast and the march at
  R_d is on the slow side of the data, so they pull the mixture towards it; the
  3 mm class (1.11 mm) costs little because it is 8 per cent of the mass. (b)
  at the bin midpoints moves it **up** by 0.024 and 0.027, because its
  mass-mean is 2.0 mm, not 1.77; the same five-class weight at the endpoints,
  which puts 20 per cent on the 1.0 mm fines, nearly cancels (-0.004, -0.005).
  The direction is set by the fine tail, not by the mean alone.
* The population shifts (at most 0.033 of the span) are 2.5 to 25 times the
  numerical uncertainty of the reported level (0.0013 and 0.008) and at most a
  quarter of the digitization half-width (0.142 and 0.157): resolvable by the
  march, not by the data.
* **What a printed histogram would settle:** the weights, and above all the
  mass fraction in the fine tail near 1 mm, which carries the sign of the population
  correction (on the mean signed residual the classes up to 1.6 mm, sunflower,
  and up to 1.2 mm, soybean, dry faster than the data, the larger classes
  slower; per-class scores in `analysis.json`). It
  would not separate geometry from D (every class depends on D/R^2 alone,
  item 29), nor settle whether each class has the mean sphericity 0.74, nor
  apply directly to the thesis samples (the histogram statement is Faner
  2019's industrial soybean meal).

## 5. The conventions side by side (60/1, as printed; continuum where marched)

| convention | sunflower: d, mm | RMS (in band) | soybean: d, mm | RMS (in band) |
|---|---:|---:|---:|---:|
| volume-equivalent sphere, Faner 2019 (the paper, 0.885 mm) | 1.77 | **0.244** (2); continuum 0.240 | 1.77 | **0.366** (0); continuum 0.358 |
| thesis volume-equivalent sphere | 1.95 | 0.291 (2) | 1.80 | 0.375 (0) |
| psi x thesis sphere (mixes the 2019 psi with the thesis d) | 1.443 | 0.133 (10) | 1.332 | 0.197 (6) |
| **psi x Faner 2019 sphere, R_d = 0.655 mm** | 1.310 | **0.079** (15); continuum 0.078 | 1.310 | **0.187** (7); continuum 0.178 |
| population (a1) at psi = 0.74 | 1.77 mass-mean | 0.055 (15) | same | 0.154 (9) |
| item 29's residual minimum (reported, not adopted) | 1.20 | 0.051 | 1.06 | 0.093 |

The thesis spheres were marched here (item 29 interpolated 0.290 and 0.374;
marched 0.291 and 0.375). The volume-equivalent conventions bracket the
comparison from the slow side (24 to 29 and 37 to 38 per cent), the
volume-to-surface conventions give 8 to 13 and 19 to 20 per cent; no
convention from the sources reaches item 29's minimum (0.60 and 0.53 mm).

## 6. Verdict against the rejecting outcomes

* **Q0 — passed.** 33,330 checks exact over 4 marches.
* **R1 — does not fire.** Every mesh chain forms at order 2.00 to 2.13; the
  sunflower 60-to-240 changes (at most 0.00037) are below 0.0005 anyway; the
  largest 60-to-240 change is below the 240/1-to-240/4 change on both traces
  (0.00037 < 0.00070; 0.00130 < 0.00199).
* **R2 — does not fire.** Every declared distribution lies within 0.029 of
  the single sphere's range (0.0777 to 0.0791 and 0.1784 to 0.1868), against
  digitization half-widths of 0.142 and 0.157 of the span.
* **R3 — does not fire.** All 68 marches cover the scored window; one refusal
  0.19 s before the 2400 s horizon, after the window.

## 7. Sentences the paper could carry (not applied)

* "The representative sphere is the volume-to-surface equivalent of the
  measured particle, R = psi R_eq = 0.74 x 0.885 mm = 0.655 mm, with psi = 0.74
  and R_eq from the same source (Faner et al. 2019); the volume-equivalent
  sphere (0.885 mm) is reported alongside."
* "At R = 0.655 mm the falling-rate window residual is 7.9 and 18.7 per cent of
  the window span (sunflower, soybean; 15 and 7 of 15 points inside the reading
  band), converged values 7.8 and 17.8 per cent, observed orders 2.0 in the
  mesh and 1.1 to 1.3 in the step, numerical uncertainty 0.1 and 0.8
  percentage points. The march still dries more slowly than both traces."
* "Carried at the law's two-sided 95 per cent prediction factor (3.14 and 3.13)
  the residual ranges over at most 5.2 to 40.0 and 9.1 to 49.6 per cent; both
  band edges lie on the far side of a minimum inside the band."
* "A population over 1 to 3 mm at the same sphericity moves the residual by at
  most 3.3 percentage points, downward when the fine tail carries 10 to 12 per
  cent of the mass and upward when the mass-mean is 2.0 mm; the sign is set by
  the fine tail near 1 mm, whose mass fraction no held source prints."
* "With the volume-equivalent sphere the residual is 24.4 and 36.6 per cent
  (thesis spheres 29.1 and 37.5); the choice of equivalent sphere is the
  largest single model-form lever in the comparison, larger than the law's
  one-sigma scatter and comparable to its 95 per cent band."

## 8. Claim and non-claim

Claimed: at R_d = psi R_eq = 0.655 mm, derived from Faner 2019's measured
sphericity and diameter and not fitted, the Faner window residual with the
boundary-value read is second order in the mesh and 1.1 to 1.3 in time,
mesh-independent within the time-step error, converged to 0.078 (sunflower)
and 0.178 (soybean) of the span, with the reported level 0.0013 and 0.0083
above; the 95 per cent band at item 28's factors gives 5.2 to 40.0 and 9.1 to
49.6 per cent; declared 1 to 3 mm populations at psi 0.74 move the residual by
-0.033 to +0.027 of the span; the conventions table above.

Not claimed: that the particles are spheres of 0.655 mm or that the
volume-to-surface sphere is the correct kernel geometry (section 9); that the
thesis samples have Faner 2019's sphericity or size distribution (the thesis
prints only 1.95 and 1.80 mm); a fitted radius, factor or distribution (the
distributions are declarations bracketing a one-line source statement, not a
histogram); a better or validated agreement (the lower residual follows from a
geometry convention proposed after item 29 showed that smaller spheres fit
better, which section 9 addresses); anything about the journal conditions,
the space-resolved solve, the other film cases or the temperature dependence;
calibration, validation, qualification or plant prediction.

## 9. Adversarial paragraph

The strongest argument against the volume-to-surface sphere as the a priori
convention is that V/A equivalence is exact only in the short-time limit: for
any convex shape the early loss per unit volume is (A/V) 2 sqrt(D t / pi)
times the surface driving difference, so matching V/A matches the early
window, but the late window is governed by the slowest eigenmode, which
depends on the shape and its smallest dimension, not on V/A. A slab of
half-thickness L has V/A = L and the equal-V/A sphere radius 3L; the slab's
slowest mode decays at (pi/2)^2 D/L^2 = 2.47 D/L^2, the sphere's at
pi^2 D/(9 L^2) = 1.10 D/L^2, a factor 2.25 apart (item 29's R/3 to R/2
ambiguity), while a compact particle (psi near 1) sits near the
volume-equivalent sphere instead. The scored window here is not short-time:
with the window's geometric-mean coefficient at R_d (5.5e-10 and 5.7e-10
m2/s) the Fourier number D t / R^2 at the last scored sample is about 0.32
and 0.33, past the sphere's short-time regime, so the equivalence is used
where it is approximate; a flake of the same psi would dry faster late in the
window than the R_d sphere, a compact one slower. A second argument is
procedural: the convention was proposed after item 29 showed that the
residual falls with the radius; a source-derived number that moves towards
the minimum is still a convention chosen with the data in view, and the paper
must present it as one of the bracketing conventions (section 5), not as a
refinement. A third: psi = 0.74 is a mean over Faner 2019's industrial soybean
meal, applied here to both thesis samples and to every size class. **The
measurement that would refute the convention** is a march of the measured
shape itself (a slab or spheroid kernel at a printed flake thickness and
aspect), or the same meal dried at two sieve cuts: if the true-shape march at
the measured geometry gives a window residual differing from the R_d sphere's
by more than the digitization half-width (0.14 to 0.16 of the span), or if
two sieve cuts do not collapse on the (V/A)^2 time scale in their early
window, the volume-to-surface sphere is not an adequate equivalent for this
comparison.
