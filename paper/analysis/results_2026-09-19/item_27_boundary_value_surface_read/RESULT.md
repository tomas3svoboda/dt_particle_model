# Item 27 — the surface law read at the boundary value of the loading

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

## 0. Rejecting outcomes, written 2026-09-26 before any march of this item ran

Brief: the principled repair of item 24's finding, as a measurement. Item 24
found the Faner falling-rate march first order in mesh and time, the printed
window residuals 9.9 and 20.6 per cent becoming 19.7 and 30.8 per cent at the
finest level and about 23 to 25 and 35 to 36 per cent extrapolated, and named
the cause: the "surface" law is read at the outer cell centre, where the
loading is 1.51 times the surface value at the reported mesh. Here the law is
read at the boundary value of the loading (the Dirichlet datum at r = R, the
value the D22 boundary-anchored line returns at the face) and item 24's 3 x 3
(60, 120, 240 cells; time level 1, 2, 4) is re-run at all four conditions.
Nothing else changes.

* **Q0 (identity first).** With the committed read switched in, item 24's four
  reported-level cells (60/1) and its four reported-level 95 per cent edges are
  reproduced by exact float identity, every summary leaf and every recorded step
  value. A single mismatch stops the script and nothing below is measured.
* **R1 (the first-order mesh dependence is not removed).** For each condition
  and each time level r in {1, 2, 4}, the mesh chain 60/120/240 of the window
  RMS over the span (full-step scoring). R1 fires if, anywhere, the 60-to-240
  change exceeds 0.0005 of the span (half the last printed digit, item 24's O3)
  and the chain either forms with an observed order below 1.5 or does not form.
* **R2 (the converged value moves).** The converged value under the new read is
  the finest corner (240/4) carried to the continuum by the same two-factor
  construction item 24 used (mesh chain at r = 4, time chain at n = 240; at order
  1 and at the observed orders; a chain that does not form contributes its
  finest difference as an uncertainty, not a correction). R2 fires if, at any
  condition, it lies farther from item 24's continuum interval (between item
  24's order-1 and observed-order estimates) than item 24's finest-level
  difference, `|item 24 order-1 estimate - item 24 240/4|` (0.0367, 0.0429,
  0.0160, 0.0124 of the span for thesis sunflower, thesis soybean, journal
  soybean, journal sunflower).
* **R3 (a level refuses that item 24 accepted).** R3 fires if any 3 x 3 level
  whose item 24 trace covered the scored window does not cover it under the new
  read (all 36 did in item 24), or if any level item 24 marched to the 2400 s
  horizon refuses before it. Refusals are recorded verbatim; no tolerance,
  floor, retry budget, horizon or level is changed to obtain a chain.
* **Reported without a threshold, as measurements:** per cell the window
  residual, the crossing time (closed form; it must be identical to item 24's
  float at every level), the ledgers; the observed orders in mesh and time;
  whether the residual is mesh-independent to within the time-step error (the
  largest 60-to-240 change against the 1-to-4 change at 240 cells); the
  threshold width (does the first-stride loading change still scale with the
  square root of the first stride); the 95 per cent prediction band at the
  reported mesh and step (item 24's factor, read from
  `rerun_2026-09-24_surface/surface_loading.json`), with item 24's O8 test (a
  residual interval wider than the window span is not resolvable).

---

## 1. What changed in the read, and what ran

* **The read.** The committed march evaluates the "surface" law once per stride
  at `ds.storage(c[-1], T, p)/rho`, the loading of the outer cell centre, half
  a cell inside r = R. Here it is evaluated at `ds.storage(c_s, T, p)/rho`, the
  loading of the Dirichlet datum `c(R) = c_s` that the same stride is marched
  against. That is the D22 boundary-anchored value at the face: the line
  through the exact boundary datum and the outer cell centre, evaluated at R,
  returns the datum (with a Dirichlet condition the outer cells set the
  gradient, not the face value). Checked at every accepted stride of every
  march: the read equals the face loading to the float, and the floating-point
  anchored line at R equals `c_s` exactly. One substituted line in the
  re-bound `march_tail`; nothing else changes (stride rule, levels, `tau`,
  tolerances, retries, horizon, film, temperature case, law, scoring).
* **Q0: 54,658 of 54,658 checks exact over 8 marches.** With the committed read
  switched in, item 24's four reported-level cells and four reported-level 95
  per cent edges reproduce every summary leaf and every recorded step value.
  The only difference below is the read.
* **The full 3 x 3** (60, 120, 240 cells x time level 1, 2, 4) at all four
  conditions under the new read, 36 marches, all scored; the 95 per cent band
  (edges at 60/1 and 240/4, one-sigma edges, a six-point sweep inside the
  band), 24 marches; a declared extension (480/1, 960/1, 60/8 on the thesis
  traces), 6 marches. 74 marches, one at a time, about 65 minutes. Ledgers:
  every accepted stride closes to at most 8.4e-15 of the particle inventory,
  every whole march to at most 3.3e-13, every remap to at most 1.3e-14.

## 2. The numbers

Window RMS as a fraction of the window span. "As printed" is the committed
decimated scoring at 60/1 (the scoring behind the paper's 9.9 and 20.6); the
other columns are full-step.

| trace | item 24: printed / finest 240/4 / continuum | **new read: 60/1 as printed** | 60/1 full-step | 240/4 | **continuum (order 1 to observed)** | 60/1 as printed minus continuum |
|---|---|---:|---:|---:|---:|---:|
| thesis sunflower 136 C | 0.099 / 0.197 / 0.234 to 0.246 | **0.244** | 0.243 | 0.241 | **0.2402 to 0.2404** | +0.004 |
| thesis soybean 136 C | 0.206 / 0.308 / 0.351 to 0.362 | **0.366** | 0.363 | 0.359 | **0.3579 to 0.3582** | +0.008 |
| journal soybean 120 C | 0.605 / 0.646 / 0.662 to 0.664 | 0.666 | 0.665 | 0.664 | 0.6636 to 0.6637 | +0.003 |
| journal sunflower 100 C | 0.554 / 0.585 / 0.598 to 0.599 | 0.602 | 0.601 | 0.599 | 0.5987 to 0.5988 | +0.003 |

**Observed orders** (DETAIL section 4): **mesh 1.96 to 2.03** at every time
level and condition (contraction ratios 3.90 to 4.08; 2.00 and 2.01 on the
extended chain 120/240/480); **time 1.05 to 1.16** (1.02 and 1.03 on r = 2, 4,
8); diagonal 1.42 to 1.56. Item 24: mesh 0.77 to 0.97, time 0.89 to 1.09. The
leading mesh error was the read; with it removed the cell-centred scheme is
second order in the mesh, and the time error (backward Euler with the
coefficient frozen per stride) stays first order.

**Mesh-independent to within the time-step error.** The largest 60-to-240
change at any time level is 0.0011 and 0.0018 of the span (thesis sunflower,
soybean; 0.0007 journal), smaller than the change from the reported to the
quarter stride at 240 cells, 0.0012 and 0.0019 (0.0008, 0.0010 journal), at
all four conditions.

**The reported level under the new read against item 24.** 24.4 and 36.6 per
cent (as printed; 24.3 and 36.3 full-step) at 60 cells and the reported step,
against item 24's finest 19.7 and 30.8 (+4.6 and +5.4 points) and its
extrapolated 23.4 to 24.6 and 35.1 to 36.2 (sunflower inside its interval,
soybean 0.4 points above). The reported level now sits within 0.4 and 0.8 points of its
own continuum value (0.3 and 0.5 on the full-step scoring, the rest output
sampling), all of one sign, reported level above the continuum. Item 24's
printed 9.9 and 20.6 were 14 and 15 points below it.

**What it means for the comparison.** The fit to the thesis traces is worse
than the printed numbers, not better: the march under the continuum read dries
more slowly than the data (mean signed residual +0.22 and +0.35 of the span;
2 of 15 and 0 of 15 points inside the reading band as printed, 1 of 15
soybean full-step; time to the last measured loading 493 and 537 s after the
crossing). This is the continuum limit item 24 predicted; item 24's printed
line owed its closeness to a first-order read error of about +50 per cent in
the coefficient.

**Crossing times.** Identical to item 24's float at every level (51.570 and
55.236 s thesis, 44.252 and 72.998 s journal): closed form, no discretization
error.

## 3. Verdict against the rejecting outcomes

* **Q0 — passed.** 54,658 checks exact.
* **R1 — does not fire.** Every mesh chain at every time level and condition
  forms at order 1.96 to 2.03 (threshold: below 1.5 with a change above 0.0005
  of the span). The first-order mesh dependence is removed.
* **R2 — does not fire.** The new continuum values lie inside item 24's
  continuum interval on thesis sunflower (0.2402 to 0.2404 in 0.2338 to 0.2459),
  thesis soybean (0.3579 to 0.3582 in 0.3514 to 0.3616) and journal soybean;
  journal sunflower lies 0.0002 above it (tolerance 0.0124). Both reads have
  the same continuum limit, as they must; item 24's extrapolation was right to
  within its own spread, and the new one has about 1/35 to 1/60 of that spread.
* **R3 — fires, on its second clause only.** No 3 x 3 level becomes unscored
  (all 36 cover the window). But two cells that reached the 2400 s horizon in
  item 24 now refuse 0.213 s short of it on the horizon-remainder stride:
  thesis sunflower 60/4 and thesis soybean 120/4 (verbatim: `RuntimeError:
  dry-shell Newton line search failed inside the admissible storage domain`,
  tail t = 2399.787 s, after seven retries that cannot enlarge the remainder).
  Conversely thesis sunflower 120/4, refused there in item 24, now reaches the
  horizon. It is the same horizon-remainder refusal class item 24 recorded at
  nine 3 x 3 cells, more than 2,000 s after the last scored sample, moved between cells;
  it touches no scored number. Also recorded: 480/1, refused at the first
  stride in item 24, now marches on both thesis traces; 960/1 refuses at the
  second stride (tail t = 3.214 s).

## 4. The threshold width under the new read

The loading removed over the first stride (item 24's switch measure) is 2.78,
3.26, 3.36 and 3.38 per cent of Xc at 60, 120, 240 and 480 cells (reported
step): it converges in the mesh to about 3.4 per cent, a third of item 24's 9.7
to 9.8 per cent, because the first stride is now marched at the face
coefficient (6.39e-9 against 5.57e-8 m2/s). **It still scales with the square
root of the first stride** on the resolved meshes (ratio 1.42 and 1.44 per
halving at 240 cells; 1.47 and 1.60 at 120); at 60 cells the ratio approaches 2
because 60 cells do not resolve the first-stride profile at the smaller
strides. The width is still a first-stride width, not a cell width, and still
not under one per cent of Xc at the reported step.

## 5. The 95 per cent prediction band under the new read

Item 24's factor F (3.165 and 3.170, read from `surface_loading.json`), D/F and
D x F at every state, at 60/1:

| trace | one-sigma edges (D/1.507, D x 1.507) | 95 % edges (D/F, D x F) | minimum marched inside the band | **95 % range** |
|---|---|---|---|---|
| sunflower 136 C | 33.9 %, 13.2 % | 47.6 % (1/15), 14.7 % (8/15) | 5.7 % at D x F^0.75 (15/15) | **at most 5.7 to 47.6 %** |
| soybean 136 C | 46.6 %, 24.2 % | 61.4 % (0/15), 10.7 % (12/15) | 10.7 % at D x F | **at most 10.7 to 61.4 %** |

Widths at least 0.42 and 0.51 of the span: item 24's O8 does not fire. At the
finest level both edges now form (47.2 and 14.9; 60.9 and 10.6 per cent), each
within 0.5 points of 60/1; item 24 could not form the band at 240/4. The
central line lies at the slow edge of what the law's 95 per cent band allows:
the data sit inside the band, toward D x F.

## 6. Sentences the paper could carry (not applied)

* Method: "The loading-dependent law is evaluated once per stride at the
  loading of the surface boundary value (the Dirichlet datum), not at the
  outermost cell; with that read the window residual is second order in the
  shell mesh and first order in the step."
* Sec. 5.1, replacing the falling-rate numbers: "The window RMS residuals are
  24 and 36 per cent of the window span (sunflower, soybean; 24.4 and 36.6 per
  cent at the reported discretization, converged values 24.0 and 35.8 per cent,
  numerical uncertainty below one percentage point). The march dries more
  slowly than both thesis traces."
* Replacing the one-sigma sentence: "Carrying the regression's
  one-standard-deviation scatter gives 13.2 to 33.9 and 24.2 to 46.6 per cent;
  at the two-sided 95 per cent prediction factor of 3.17 the residual ranges
  over at most 5.7 to 47.6 and 10.7 to 61.4 per cent, and the measured curves
  lie inside the band toward its fast-drying edge."
* The crossing ratios 0.89 to 1.12 carry no discretization error (unchanged).
* Threshold (Faner march): "the loading removed while the surface state is
  taken up converges in the mesh to about 3 per cent of Xc at the reported step
  and scales with the square root of the first step; it is a time-step width,
  not a measured physical sharpness."
* Any sentence that relies on the printed 9.9 and 20.6 per cent, or on 14 of
  15 sunflower points inside the reading band, is withdrawn by this measurement.

## 7. Claim and non-claim

Claimed: the boundary-value read reproduces item 24 bit for bit when switched
back (54,658 checks); under it the experiment-facing window residual is second
order in the mesh and first order in time, mesh-independent to within the
time-step error, converged to 24.0 and 35.8 per cent of the span on the thesis
traces (66.4 and 59.9 journal) with the reported level 0.4 and 0.8 points above;
it shares its continuum limit with the committed read (item 24's extrapolation
confirmed); the crossing times are exact; the switch width is a first-stride
width of about 3 per cent of Xc; the 95 per cent band gives the ranges above.
R3 fired on two horizon-remainder refusals after the window.

Not claimed: a better agreement with the data (the converged comparison is
worse than the printed one); a space-resolved D(X(r)) (the certified shell
still takes one scalar per stride; the boundary value is the continuum limit
of the committed surface convention, not the variable-coefficient solution);
convergence past 480 cells (960 refuses at the second stride); anything about
the moving-front solver; a fit, calibration or validation (nothing was
fitted); physical qualification or plant prediction.
