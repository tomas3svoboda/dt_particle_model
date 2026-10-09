# Item 26 — space-time refinement of the dry-particle uptake benchmark

## 0. Rejecting outcomes (written 2026-09-26 before any refined level ran)

These were fixed while only the reported level (m1 k1) was running, and they
are applied below exactly as written.

- **R1 — identity.** If the reported level does not reproduce every committed
  float of item 13, of item 05's converted benchmark and of the form-mapping
  fit, the refinement is not a refinement of the printed pipeline: stop and
  report, no uncertainty is claimed.
- **R2 — refusal.** Any typed refusal of the certified step on a refined level
  is reported verbatim; that level is withheld from every chain it belongs to,
  and a printed number whose chain is broken carries no numerical uncertainty
  from this item.
- **R3 — ledger.** A run whose largest per-step inventory-against-face-flux
  ledger exceeds 1e-9 of the equilibrium inventory, or whose cumulative ledger
  exceeds 1e-8, is not admissible evidence.
- **R4 — materiality.** A printed number is called numerically resolved at its
  printed precision only if the reported-to-finest difference is below half a
  unit of its last printed digit. If the difference in fractional uptake at a
  measured time exceeds the digitization's declared ±0.010 band, or a printed
  ratio moves by more than 5 per cent, the comparison is numerically limited
  and the paper's "nor is it a discretization artefact" sentence is refuted.
- **R5 — orders.** An observed order is reported only for a three-level chain
  whose two differences have the same sign and contract (ratio above 1);
  otherwise "no asymptotic order" is written and the uncertainty is the plain
  reported-to-finest difference, never an extrapolation.
- **R6 — counts.** If an integer in-band count differs between levels, the
  printed count is not numerically resolved and its range over levels is
  printed instead.
- **R7 — isotherm.** If the isotherm source prints no uncertainty for its fit,
  the identification branch stops there; no uncertainty is invented.

## 1. Outcome

**The printed uptake comparison is not discretization-limited, but several of
its printed digits are not resolved.** Every run is admissible (R3), no run
ends in a refusal, every in-band count at the declared coefficients is
identical on every level, and the largest change in simulated fractional
uptake at any measured time is 0.0053, half the digitization's declared
±0.010 box (R4 materiality not triggered; the paper's "nor is it a
discretization artefact" stands). The error is time-step dominated (backward
Euler, observed order 1.0); the mesh contributes an order of magnitude less
(observed order 2.00). The upper ends of the 7.6 to 18.9 factor and of the
8.53 to 19.54 retardation move by 1.3 and 1.4 per cent, and the converted
benchmark's "35 of 100" is 34 or 35 depending on the time step (R6).

## 2. What ran

- **Identity first.** The reported level (60 cells for the item-13 and item-05
  marches, 240 for the form-mapping master curves, each at its own step
  controller) reproduced the committed outputs to float identity:
  **128,514 checks, 0 failures** (`outputs/identity_check.json`): all 34
  item-13 residual rows x 24 numeric fields and status plus every point of its
  committed thinned traces; the 17 converted rows of item 05's
  `uptake_benchmark_residuals.csv` and every point of their committed traces;
  the 8 form-mapping master curves point for point and its 13 fit rows.
- **All series, no subset.** Every series the paper scores: 17 series x 2
  declared coefficients (contract A, item 13); 17 series at the converted
  coefficient (contract B, item 05) in two variants, *pipeline* (the
  retardation re-measured at the same level, so `D* = r D_table` is recomputed
  there, as the paper's procedure would) and *fixed* (the printed `D*` held, so
  only the re-march is refined); and the 13-series form-mapping fit on 8
  master curves (contract C).
- **Levels.** Mesh x1, x2, x4 (uniform cells); time x1, 1/2, 1/4 (every step of
  the reported controller split into 2 or 4 equal backward-Euler sub-steps,
  the controller itself unchanged, so level k1 is the reported march). Run:
  m1k1, m2k1, m4k1, m1k2, m1k4, m2k2, m4k4, i.e. the diagonal and both
  single-axis chains. 515 marches (59 at the reported level, 76 at each
  refined one), one at a time, 9,799 s wall. Not run: m2k4
  and m4k2.
- **Refusal (R2).** One sub-step at m4k4 (figure 3.17, a_h 0.423, tabulated
  coefficient; first step, sub-stride 0.0340 s) was refused by the certified
  step, verbatim `RuntimeError: dry-shell Newton has no admissible
  positive-density update; reject, do not clamp` (`outputs/retry_probe.json`);
  the contract retried it at twice the stride and the march completed, so that
  run is no longer a pure step-quartering. It and the converted run its
  retardation feeds are withheld from every chain at m4k4. Figure 3.17 is the
  duplicate of 3.13, so no independent span is affected. **Departure from R2 as
  written, stated:** R2 says "that level" is withheld; it was applied to the
  affected run (and its dependent converted run) only, not to the other 74 runs
  of m4k4, because the refusal was recovered by the contract's own retry rule
  and ended no march. Read literally, R2 removes m4k4 and the finest level
  becomes m1k4 (time) and m4k1 (mesh); by `outputs/printed_numbers_by_level.csv`
  every verdict in section 4 is the same on m1k4, the finest values differing
  by at most one unit of the last digit shown (frozen time-ratio upper end
  18.253 against 18.243, retardation 8.516–19.282 against 8.513–19.271, form
  factor 18.655 against 18.654).
  No other retry and no terminal refusal on any level.
- **Ledgers (R3).** Inventory change against the Dirichlet-face flux on the
  converged state, per step: largest 1.70e-10 of the equilibrium inventory,
  cumulative at most 3.71e-9, on every level. Admissible.

## 3. Observed orders (R5)

Per-series three-level chains, ratio 2 (`outputs/analysis.json`,
`outputs/per_series_refinement.csv`):

| quantity | mesh chain (k1) | time chain (m1) | diagonal |
|---|---|---|---|
| A: t63, time ratio, retardation (34 rows) | 2.00 (34/34) | 0.97–1.05 (34/34) | 1.02–1.33 (33/33) |
| A: rms residual | 2.00–2.04 | 0.99–1.04 | 0.86–1.80 |
| A: time to 0.5 and 0.9 uptake | 2.00–2.01 | 0.99–1.06 | 1.01–1.52 |
| B fixed D*: t63, time ratio | 2.00–2.01 | 1.00–1.06 | 1.04–1.15 |
| B: rms residual | 2.00–2.03 | 0.59–1.50 | 0.99–2.27 |
| B pipeline: D* | 2.00 | 0.97–1.05 | 1.02–1.33 |
| C: fitted ratio (13 series) | 2.00 (9; 4 unchanged to 5 digits) | **no asymptotic order**: 0.11–1.87 on 12, one not contracting | same |
| C: fit rms | 2.00 | 0.93–1.05 | 0.91–1.05 |

The marches converge at the formal orders (second in space, first in time).
The fitted form-mapping coefficient shows no clean order in time (it is an
argmin over a handful of points on an interpolated master curve), so for it
only the plain reported-to-finest difference is stated, and that is a lower
bound on the reported level's error.

## 4. The numerical uncertainty of each printed number

Reported (m1k1, = printed) → finest (m4k4) → formal-order Richardson estimate
(`q* = q11 - 2(q11 - q12) - 4/3 (q11 - q21)`, valid where section 3 shows the
orders; for C's fitted ratio shown for scale only).

| printed | reported | finest m4k4 | Richardson | verdict |
|---|---|---|---|---|
| factor 7.6 to 18.9 (13 series) | 7.621–18.898 | 7.594–18.654 | 7.592–18.53 | upper end −0.24 (−1.3 %) at least; prints 7.6–18.7 |
| 7.6 to 18.6 (without 3.17) | 7.621–18.617 | 7.594–18.425 | 7.592–18.37 | prints 7.6–18.4 |
| soybean 7.6 to 17.9 | 7.621–17.920 | 7.594–17.683 | 7.592–17.57 | prints 7.6–17.7 |
| fit rms 0.040 to 0.115 | 0.0403–0.1146 | 0.0408–0.1139 | 0.0409–0.1137 | prints 0.041–0.114; still the analytic sphere's 0.041–0.112 |
| declared time ratio 0.58 to 18.3 (frozen) | 0.5758–18.295 | 0.5738–18.243 | 0.5732–18.23 | prints 0.57–18.2 |
| declared time ratio 6.1 to 24.4 (tabulated) | 6.131–24.426 | 6.105–24.286 | 6.096–24.24 | prints 6.1–24.3 |
| declared rms (frozen / tabulated) | 0.096–0.635 / 0.445–1.193 | 0.097–0.632 / 0.444–1.191 | 0.098–0.631 / 0.444–1.191 | third digit only |
| 11 and 3 of 100 in the declared box; 0 in the two-reading band | 11, 3, 0 | 11, 3, 0 | — | **resolved on every level** |
| "a factor of 32" (frozen band) | 31.77 | 31.79 | 31.80 | resolved |
| retardation 8.53 to 19.54 | 8.534–19.541 | 8.513–19.271 | 8.506–19.19 | upper end −0.27 (−1.4 %); prints 8.51–19.27 |
| converted rms 0.041 to 0.297 | 0.0405–0.2970 | 0.0414–0.2968 (pipeline) | 0.0416–0.2967 | resolved at printed precision |
| converted time ratio 0.57 to 2.19 | 0.5656–2.1946 | 0.5554–2.1958 (pipeline); 0.5529–2.1832 (fixed) | 0.552–2.196 | prints 0.56–2.20 (pipeline) |
| band "contracts to 3.9" | 3.880 | 3.954 (pipeline); 3.948 (fixed) | 3.98 | **3.9 to 4.0**, not resolved |
| converted, 35 of 100 in the declared box | 35 | 34 at every k2 and k4 level (pipeline); 35 on every level (fixed) | — | **34 to 35** (R6) |
| converted, 2 in the two-reading band | 2 | 1 at k4 (both variants) | — | **1 to 2** (R6) |
| converted coefficients D* (panel d) | 2.620e-10–7.193e-9 | 2.608e-10–7.173e-9 | 2.604e-10–7.168e-9 | −0.5 to −1.4 % |

Materiality (R4): the largest change of simulated fractional uptake at any
measured time, reported against finest, is 0.0045 (declared coefficients),
0.0053 (converted, fixed D*), 0.0047 (converted, pipeline); the master curves
move by at most 0.0031. No printed ratio moves by 5 per cent (largest 2.6 per
cent, a converted t63). **The comparison is not numerically limited; the third
digits of the upper ends are.** The Table S16 per-figure spans at every level
are in `outputs/printed_numbers_by_level.csv` and
`outputs/printed_spans_summary.txt`.

## 5. The identification's remaining uncertainty: the isotherm (R7)

The storage capacity behind the conversion is the frozen Cardarelli and
Crapiste (1996) native-soybean GAB (`paper/analysis/datasets/cardarelli1996/`).
That source prints the five GAB constants (its Table 2) with **no standard
error, confidence interval or fit statistic** (checked in the source text);
its dataset record states the authors print no pointwise uncertainty either,
and the bounds held there are digitization bounds, not a fit uncertainty. The
thesis's own per-temperature refit (Tabla 2.6,
`cardarelli1998_thesis_ch4_and_eq224/`) is likewise printed without one ("the
source supplies no uncertainty, no replicate structure and no error budget").
**The branch stops: no isotherm fit uncertainty is printed, so none is
propagated and the movement of the 7.6 to 18.9 factor under it is not
measured.** What exists is a second fit by the same author whose 50 °C
constants differ from the 1996 global fit (k within 0.4 to 0.9 per cent, C by
12 to 32 per cent, per that dataset record); that is a fit-to-fit divergence,
not an uncertainty, and running it would be a sensitivity case for the owner to
ask for, not this measurement.

## 6. Sentences the paper could carry

- "Refined to four times the mesh and a quarter of the time step, the uptake
  marches converge at second order in space and first in time; the simulated
  fractional uptake at every measured time moves by at most 0.005, half the
  digitization's declared reading box, and every in-band count at the declared
  coefficients is unchanged."
- "The factor between the storage-consistent and the source's coefficient is
  7.6 to 18.7 at the finest level (7.6 to 18.9 at the reported one), and the
  retardation 8.51 to 19.27 (8.53 to 19.54)."
- "At the converted coefficients 34 to 35 of the hundred points fall inside the
  declared box depending on the time step, and the time-constant band
  contracts to 3.9 to 4.0."
- "The isotherm behind the conversion is printed without a fit uncertainty, so
  the conversion's dependence on it is not quantified here."

## 7. Claim and non-claim

Claim: the numerical uncertainty of the printed dry-particle uptake numbers,
measured on every scored series at seven space-time levels of the three
contracts that print them, from a float-identical reproduction. Non-claim: no
physical validation, no identification uncertainty (section 5), no change to
any printed number, coefficient, law or tolerance; one temperature, one
radius, digitized points; `physically_qualifying: false`,
`plant_predictive: false`. Levels m2k4 and m4k2 not run.

Detail: `DETAIL.md`. Files: `MANIFEST.sha256`.
