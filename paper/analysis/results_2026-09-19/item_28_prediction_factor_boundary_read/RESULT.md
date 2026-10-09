# Item 28 — the 95 per cent prediction factor re-derived at the boundary-value read

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

Date 2026-09-26. Asked for by the fourth referee read of the particle
manuscript (`paper/01_particle_jfpe/REFEREE_AUDIT_4_2026-09-26.md`, item M2):
the two-sided 95 per cent prediction factor on the diffusivity law that item 27
carried through its band (3.165 sunflower, 3.170 soybean) was derived by
`item_05_faner_comparison/rerun_2026-09-24_surface/run_surface_loading.py` at
the loading of the **outer cell centre**, the read item 27 replaced. The
referee's concern: the boundary loading is 1/1.51 of the outer-cell loading at
the reported mesh, the outer-cell window minimum is 1.05e-3 kg/kg against a
measured minimum of 7.677e-4, so the boundary read might sit below the
measured loading range, at a larger leverage and a larger factor, and "at most
5.7 to 47.6 per cent" would then not bound the band.

Nothing in any other folder was edited or written. Item 27's driver was
imported read only, with its `RUNS` output directory re-pointed at this
folder's `runs/`.

## 1. The construction (unchanged) and its identity check

Exactly the surface script's: design row `x = [1, -1/(Rg T), ln(X/0.01)]` in
the law's own parameterization, leverage `h = x'(X'X)^-1 x` with the recorded
`cov_unscaled` of the soybean fit (n = 17), factor
`F(h) = exp(t95 sigma sqrt(1 + h))` with t95 = 2.145 and sigma = 0.410374, and
F carried = F at the **largest leverage over the accepted strides in the scored
window** (absolute time between the first measured sample at or after the
measured crossing, 60.04 and 60.03 s, and the last, 300.12 s). The loading is
the one the law was read at in that stride, and T the stride's model
temperature.

* `(X'X) C = I` on the seventeen points to 6.0e-10.
* **Identity:** applied to the committed outer-cell steps
  (`surface_loading_steps.csv`, `lit_whitaker`), the construction reproduces
  the committed window step counts (48, 48), leverage extremes
  (0.33014-0.71330, 0.32357-0.71772) and factors **3.1651154925439124 and
  3.169818627063594 by exact float equality**. The only change below is the
  loading and temperature the leverage is evaluated at.

## 2. The factor at the boundary-value read

Item 27's reported-level carried line (60 cells, time level 1; the read equals
the face datum to the float at every stride, checked again here):

| trace | F, outer-cell read (item 24/27) | **F, boundary read** | change | leverage in window (outer / boundary) | largest leverage at |
|---|---:|---:|---:|---|---|
| sunflower 136 C | 3.1651 | **3.1401** | **-0.79 %** | 0.33-0.71 / **0.20-0.69** | 296.9 s, X = 7.92e-4 kg/kg, 134.4 C |
| soybean 136 C | 3.1698 | **3.1323** | **-1.18 %** | 0.32-0.72 / **0.20-0.68** | 295.6 s, X = 8.01e-4 kg/kg, 134.0 C |

Largest fitted leverage 0.3046 (factor 2.733); 9 and 10 of the 48 window
strides are at or below it. The factor is **lower** at the boundary read, not
higher: the window's largest leverage falls from 0.713 and 0.718 to 0.690 and
0.683, because at the boundary read the high-loading opening strides (up to
2.19e-2 and 3.11e-2 kg/kg at the outer cell) are gone and the largest leverage
now sits at the window's last strides, at the lowest loading and the highest
temperature.

The factor is insensitive to the discretization: 3.14011, 3.14014, 3.14015,
3.14015 (sunflower) and 3.13234, 3.13238, 3.13239, 3.13239 (soybean) at 60,
120, 240 and 480 cells (time level 1); 3.1417 and 3.1345 at the half step. On
the band-edge marches themselves it is 3.1448 and 3.1349 (sunflower D/F,
D x F) and 3.1373 and 3.1269 (soybean).

**The soybean change exceeds one per cent** (the brief's threshold), so the
band edges were re-marched (section 4).

## 3. Does the boundary loading leave the measured range? (measured)

Measured loading range of the law: 7.677e-4 to 1.5616e-2 kg/kg.

| trace, reported level | boundary loading over the window, kg/kg | min / measured min | strides inside the measured range | duration inside | first stride below the measured minimum |
|---|---|---:|---:|---:|---|
| sunflower 136 C | 7.920e-4 to 1.022e-2 | 1.032 | **48 of 48** | **1.000** | tail time 710.4 s (absolute about 762 s), after the window |
| soybean 136 C | 8.011e-4 to 1.444e-2 | 1.043 | **48 of 48** | **1.000** | tail time 780.4 s (absolute about 836 s), after the window |

Inside the scored window the boundary loading **never falls below the measured
minimum** at any stride, at any of the seven refinement levels read (60/1 to
480/1, 60/2 to 240/2) or on the 60/1 band edges; the lowest window value is
3.2 and 4.3 per cent above the measured minimum. It falls below the measured
minimum only after the window closes (the last scored sample is at 300.1 s),
at tail times 710.4 and 780.4 s on the carried line (430 to 838 s on the
edges), reaching 7.63e-4 kg/kg by the 2400 s horizon; no scored number depends
on those strides. The outer-cell read of the committed record held 47 and 46
of 48 strides inside (98.3 and 96.8 per cent of the duration, sunflower and
soybean), its misses
being above the measured maximum at the window's opening. The temperature is
unchanged in kind: 5 and 7 of the 48 window strides lie inside the fitted 50
to 95 C.

## 4. The band edges re-marched at the new factor (reported level)

Item 27's driver (`run_item27.run_march`, `read_boundary`, item 24's
`ScaledCase(ln_shift = -/+ ln F)`), 60 cells, time level 1.

* **Identity first:** item 27's four reported-level 95 per cent edges,
  re-marched here at item 27's factor, reproduce item 27's summaries and every
  recorded step value: **33,248 of 33,248 checks exact over 4 marches**
  (`identity27_*` in `runs/`).
* **The four edges at the new factor**, all to the 2400 s horizon with no
  retry; every stride's read equals the face datum; step ledgers at most
  7.8e-15 and whole-march closure at most 2.7e-13 of the particle inventory.

Window RMS as a fraction of the span; item 27 printed the band on the
full-step scoring (RMS and band count), so that is what is carried; the
decimated scoring is beside it.

| trace | edge | ln shift | item 27 (F outer-cell), full-step | **item 28 (F boundary), full-step** | band count, full-step (decimated) | decimated RMS | mean signed, full-step |
|---|---|---:|---:|---:|---:|---:|---:|
| sunflower | D/F | -1.1443 | 0.4756 | **0.4744** | 1 (1) of 15 | 0.4752 | +0.443 |
| sunflower | D x F | +1.1443 | 0.1469 | **0.1439** | 8 (8) of 15 | 0.1426 | -0.129 |
| soybean | D/F | -1.1418 | 0.6141 | **0.6121** | 0 (0) of 15 | 0.6135 | +0.589 |
| soybean | D x F | +1.1418 | 0.1071 | **0.1043** | 12 (13) of 15 | 0.1060 | -0.029 |

**The band at the boundary-read factor.** Sunflower: the edges give 47.4 and
14.4 per cent; item 27's sweep minimum, 5.7 per cent at +0.75 of its ln F
(ln shift +0.864, 15 of 15 in band), lies inside the new band (|0.864| <
1.144), so the band's range is **at most 5.7 to 47.4 per cent** (item 27: 5.7
to 47.6). Soybean: the edges give 61.2 and 10.4 per cent; the residual is not
monotone in the multiplier (item 27's sweep: 10.8 per cent at +0.865, 10.7 at
+1.154, and here 10.4 at +1.142), so the minimum lies inside the band near its
fast-drying edge and was not refined further: **at most 10.4 to 61.2 per
cent** (item 27: 10.7 to 61.4, where "monotone to D x F" is corrected by this
measurement). Widths at least 0.42 and 0.51 of the span, below the full span,
so item 24's O8 test (a residual interval wider than the window span is not
resolvable) does not fire. The measured curves still lie inside the band
towards its fast-drying edge (mean signed residual, full-step, -0.129 and
-0.029 of the span at D x F, +0.443 and +0.589 at D/F).

Not re-marched: the finest-level (240/4) edges at the new factor (item 27
formed them at the old factor within 0.5 points of 60/1), and item 27's
declared sweep (its interior shifts all lie inside the new band). The
one-sigma scatter band (factor 1.507, no leverage) is unchanged: 13.2 to 33.9
and 24.2 to 46.6 per cent.

## 5. Verdict

* The referee's concern does **not** materialize: at the boundary read the
  law is read inside its measured loading range at every stride of the scored
  window (48 of 48, both traces), the largest leverage is lower (0.690 and
  0.683 against 0.713 and 0.718), and the factor is lower, **3.140 and 3.132**
  against 3.165 and 3.170 (-0.8 and -1.2 per cent).
* Because the soybean change exceeds one per cent, the band printed at the
  boundary read is now the one marched at the boundary-read factor: **at most
  5.7 to 47.4 per cent (sunflower) and 10.4 to 61.2 per cent (soybean)**.
* The factor, the leverage and the coverage are functions of where the law is
  read; they are properties of the one-scalar surface convention, not of the
  law alone.

## 6. Claim and non-claim

Claimed: the construction reproduces the committed factor to the float; at
item 27's boundary read the factor is 3.1401 and 3.1323, the leverage 0.20 to
0.69 and 0.20 to 0.68, the window loading inside the measured range at every
stride; the band edges at that factor give the ranges above, reported level
only.

Not claimed: a new law, a refit or any change to sigma, t95 or the covariance;
convergence of the new edges (the factor itself moves by at most 3e-5 over
60 to 480 cells, but the edges were marched only at 60/1); coverage of the
41 K temperature extrapolation, which no band covers; anything about a
space-resolved coefficient; calibration, validation, qualification or plant
prediction.

## 7. Files

* `derive_factor_item28.py` — the factor, the leverage and the coverage;
  writes `factor_item28.json` and `factor_item28_steps.csv` (the 96 window
  strides of the two reported-level carried lines with the boundary and
  outer-cell loadings and the leverage). No march.
* `run_band_item28.py` — identity re-march of item 27's four edges, then the
  four edges at the new factor; writes `band_item28.json`, `log_band28.txt`
  and `runs/` (eight summaries and their per-step CSVs).
* `MANIFEST.sha256`.

Run order (from this folder, pins set: PYTHONHASHSEED, OMP, OPENBLAS, MKL,
NUMEXPR = 1, PYTHONPATH = the repository's `src`):
`derive_factor_item28.py`, then `run_band_item28.py` (about 100 s of marching,
one at a time).
