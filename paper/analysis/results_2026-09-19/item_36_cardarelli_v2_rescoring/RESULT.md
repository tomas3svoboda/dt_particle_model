# Item 36 — the Cardarelli (1998) uptake benchmark re-scored and re-fitted on the version-2 digitization

Date 2026-09-27. Consumers read-only: item 13, item 05's converted benchmark
(rerun 2026-09-21), the form mapping (docs/evidence/form_mapping_2026-09-15),
item 26. Data: `paper/analysis/datasets/cardarelli1998_thesis_ch3_v2/`.
Nothing in the paper, in docs/ or in any other item was written.

## 1. Outcome

The version-2 points move every printed number that compares the model with
the measured points, and none of the numbers that do not. The form-mapping
factor becomes **8.4 to 20.1** (8.4 to 20.0 refined) on the independent set,
against the printed 7.6 to 18.6 (18.4). The converted benchmark gets
visibly better: rms 0.033 to 0.146, 63 % time ratios 0.86 to 2.12, the band
**2.5** instead of 3.9 to 4.0, and **52 of 124** points in the box (51 to 52
refined) instead of 34 to 35 of 100. The retardation 8.51 to 19.27 (refined) does
not move at its printed precision. At the declared coefficients **22 of 124 and 3 of 124** points are
in the box, where the paper prints 11 and 3 of 100. The frozen upper end is
**17.5** where the paper prints 18.3. The old 13-series set on the version-2
points gives the same ranges except where the Figura 3.19 arm sets an end
(section 3, last column), so each change comes from the points, not from the
set.

## 2. Which inputs come from the curves (brief step 1)

| input | from the curves? | consequence |
|---|---|---|
| measured points (t, W/W_inf), their per-point box, the readings band | yes | every residual, count and fit **re-scored / re-fitted** |
| measured t63 (linear between bracketing points, item 13 rule) | yes | every time ratio re-scored; `t63_table.csv` reproduces the dataset record's t63_v2 to its printed digits for all 18 series |
| the series set (which series exist with at least 2 points; duplicates) | yes | independent set = 3.8, 3.9, 3.13: 12 series, 124 points; 3.19 open circle now posable (9 points) but a re-plot, never counted |
| **measured window of each series** (last measured time; the step controller of contracts A and B caps the stride at window/400 inside it) | **yes** | a **march input**, so contracts A and B were **re-marched** at the reported level and at item 26's finest level m4k4 (see below). The windows move by up to a factor 4.5 (the reassigned 1.06 h and 2.0 h corner markers, the added late markers), yet the effect on the marches is small: on the independent set the simulated t63 moves by at most 0.71 %, a single retardation by at most 0.089 and an rms_rel by at most 0.0016 at the reported level (all rows: 1.4 %, 0.275, on the re-plot 3.17 a_h 0,813 tabulated, whose window fell from 3822 to 1309 s); the retardation span's ends move by 0.014 and 0.001 (`analysis.json`, `window_effect_reported_level`) |
| D_eff (Tablas 3.6/3.7), activities, T = 50 C, r_p = 0.67 mm, capacity | no | unchanged |
| retardation r and D* = r D_table | only through the window | the printed values of r and D* are unchanged at printed precision; D* re-derived per level (pipeline) and also held at the printed value (fixed) |
| form-mapping master curves (D t master, 240 cells, no window) | no | reused as committed at the reported level; re-marched at m4k4; only the **fit** is redone |

So: **re-scored** = the existing committed marches, scored against the v2
points (section 3, note 5); **re-fitted** = the form-mapping
identification (item 13's "converted coefficients" are not a fit, they are
r D_table, so there is nothing to re-fit there); **re-marched** = contracts A
and B with the v2 windows at m1k1 and m4k4, and the master curves at m4k4.
Previous runs untouched.

## 3. Every printed number, before and after

Before = printed (v1 points, 13-series set of 100 points; refined value in
brackets from item 26 m4k4). After = v2 points, **independent 12 series / 124
points**, re-marched with v2 windows at the reported level (refined m4k4 in
brackets). Last column = v2 points on the **old 13-series set** (134 points),
reported level, for attribution. Converted = pipeline D*.

| printed number | before (printed) | v1 points, independent subset 12 / 92 | **after, independent 12 / 124** | v2, old 13 / 134 |
|---|---|---|---|---|
| form-mapping factor (paper's 9 series without 3.17) | 7.6–18.6 (7.6–18.4) | 7.6–18.6 (8 series, 56 pts) | **8.43–20.15 (8.36–19.96)**, 8 series, 76 pts | 8.43–20.15 (paper's 9); 8.43–20.25 with 3.17 |
| same, 13-series fit incl. 3.17 (supplement) | 7.6–18.9 (7.6–18.7) | — | as above | 8.43–20.25 (8.36–20.05) |
| soybean factor | 7.6–17.9 (7.6–17.7) | same | **8.67–20.14 (8.62–19.84)** | same |
| fit rms, present form | 0.040–0.115 (0.041–0.114) | 0.040–0.115 | **0.013–0.098 (0.013–0.097)** | 0.013–0.099 |
| fit rms, analytic sphere, same markers | 0.041–0.112 | 0.041–0.112 | **0.015–0.085** | 0.015–0.091 |
| retardation (all rows) | 8.53–19.54 (8.51–19.27) | same | **8.52–19.54 (8.51–19.27)** | same |
| declared-coefficient time ratio, frozen | 0.58–18.3 (0.57–18.2) | same | **0.589–17.51 (0.586–17.46)** | 0.579–17.51 |
| declared-coefficient time ratio, tabulated | 6.1–24.4 (6.1–24.3) | same | **7.77–23.49 (7.74–23.40)** | 7.64–23.49 |
| frozen time-constant band factor ("32") | 31.77 (31.79) | 31.77 | **29.7 (29.8)** | 30.2 |
| declared rms_rel, frozen / tabulated | 0.096–0.635 / 0.445–1.193 | same | **0.081–1.225 / 0.342–1.236** (raw rms 0.061–0.456 / 0.305–0.496) | same |
| in declared box, frozen / tabulated / readings band | 11 / 3 / 0 of 100 | 9 / 3 / 0 of 92 | **22 / 3 / 10 of 124 (22 / 3 / 9)** | 24 / 3 / 11 of 134 |
| converted rms_rel | 0.041–0.297 (resolved) | same | **0.033–0.146 (0.034–0.145)** | same |
| converted 63 % time ratio | 0.57–2.19 (0.56–2.20) | same | **0.862–2.125 (0.862–2.119)** | same |
| converted band ("3.9") | 3.88 (3.9–4.0) | 3.88 | **2.46 (2.46)** | 2.46 |
| converted, in declared box | 35 of 100 (34–35) | 35 of 92 | **52 of 124 (52 pipeline, 51 fixed D\*: 51–52)** | 53 of 134 |
| converted, in readings band | 2 (1–2) | 2 | **14 (14)** | 14 |
| converted D* (panel d) | 2.620e-10–7.193e-9 | same | 2.620e-10–7.193e-9 (2.608e-10–7.173e-9) | same |

Notes. (1) The readings band is now the min/max over the up to three
readings (v1, the audit, reading 3), not version 1's two automatic readings
that agreed to 2e-4; the "0 in the two-reading band" statement has no
like-for-like successor, and its v2 analogue is 10 (9 refined), not
resolved in the sense of item 26 R6. (2) The declared box is version 2's
per-point half-width (the larger of the v1 box and the readings' half-range);
the v1 per-figure box gives the same counts (22, 3; converted 52). (3) The
frozen rms_rel upper end 1.225 is Figura 3.9 a_h 0.80, whose loading range
fell from 0.95 to 0.37 because its t = 0 marker is now a 3.9 a_h 0.60 triangle.
The raw rms of that series barely moves (0.42 to 0.46), so the rms_rel upper
end is a normalization effect and should not be printed as a
worsening. (4) What sets each new end: form factor 3.13 a_h 0,212 (8.43) and 3.13
a_h 0,813 (20.15), with 3.9 a_h 0.80 at 20.14; converted ratio 3.13 a_h 0,621 (0.862,
the old 0.57 end was 3.9 a_h 0.60, whose t63 fell from 0.0589 to 0.0311 h) and
3.8 a_h 0.60 (2.12); frozen upper end 3.9 a_h 0.80; tabulated lower end 3.13
a_h 0,212 (on the old set 3.19's arm, 7.64).
(5) Re-scoring only, i.e. the committed marches at their committed inputs
(v1 windows, printed D*) scored against the v2 points, independent set:
frozen / tabulated time ratios 0.588–17.51 / 7.79–23.55, box counts 22 / 3,
readings band 10, retardation 8.53–19.54 (unchanged by construction);
converted rms 0.033–0.145, ratios 0.862–2.116, band 2.46, box 52, readings
band 15; the form-mapping re-fit is identical to the column above, because the
master curves carry no window. Every re-scored span end is within 0.5 %
of the re-marched one (largest: the converted upper end 2.116 against 2.125,
the tabulated upper end 23.55 against 23.49), and every count is equal, except
the converted readings band, which is 15 against 14; this is the window effect of section 2.

## 4. Sample-description wording correction

"set of spheres" / "the soybean set is a packed cylinder read as a set of
spheres" should read: **the Figura 3.8 (sunflower) runs are 10 mm wire-mesh-cylinder
runs (d_c = 0.01 m) that the source tabulates twice, as a cylinder and,
refitted with the sphere solution at d_p = 1.34 mm, as a set of spherical
particles; the paper uses the latter row. For soybean (Figura 3.9) the same
construction is inferred from Tabla 3.7, not shown by a cross-plot.**
Figura 3.13 is isolated particles (a meal layer on a basket), mean radius 0.67 mm.

## 5. Independent set

Figuras 3.8, 3.9 and 3.13: **12 series, 124 points** (47 + 41 + 36, origins
excluded). Figura 3.17 (all four series) and both arms of 3.19 are re-plots
(3.17 and 3.19's particle arm of 3.13's runs, 3.19's cylinder arm of 3.8's).
They are scored and reported, carry `independent = no`, and are never
counted. The form mapping's committed 13 series held 8 independent ones
(3.13 x 4, 3.9 x 4); on version 2 they carry 76 fitted points (one 3.9 a_h 0.60
marker at t < 0 is dropped by the fit's own filter). The paper's 13 series and
100 points were 12 series and 92 points independent on version 1.

## 6. Sentences the paper can carry

- "The data are twelve independent uptake series at 50 °C from that thesis
  with 124 measured points. A further six series re-plot the same runs and are
  scored but not counted."
- "At either declared coefficient the simulated curves stay apart from the
  measured ones; 22 and 3 of the 124 points fall inside the declared
  uncertainty box."
- "With the coefficient fitted per series the present form reproduces each
  series as well as the analytic sphere does (rms 0.013 to 0.098 against
  0.015 to 0.085), with a coefficient 8.4 to 20.1 times the tabulated value
  (8.4 to 20.0 refined)."
- "Over the simulated cases the retardation is 8.52 to 19.54 (8.51 to 19.27
  refined)."
- "At the converted coefficients the residuals fall to 0.033 to 0.146 and the
  63 % time ratios to 0.86 to 2.12; the time-constant band, a factor of 30 at
  the reference coefficient alone, contracts to 2.5, and 52 of the 124 points
  (51 to 52 refined) fall inside the declared box, against 22 and 3 at the
  declared coefficients."
- Frozen time ratios "0.59 to 17.5 (0.59 to 17.5 refined)"; tabulated "7.8 to
  23.5 (7.7 to 23.4)".

## 7. Health

627 float-identity checks at the reported level (every existing A and B march
re-driven with its committed inputs reproduces item 13's and item 05's
committed floats; the v1 form-mapping fit reproduces s2_partA_fit.csv, and item
26's m1k1 fit) and 52 at m4k4 (the re-marched masters give item 26's m4k4 fit
to the float): **0 failures**. No march refused. One m4k4 march retried a
reported step: 3.19 particle filled triangle, tabulated coefficient, sub-stride
0.096057 s, verbatim `RuntimeError: dry-shell Newton has no admissible
positive-density update; reject, do not clamp`. The contract doubled the
stride and the march completed. Under item 26's R2 that run and the converted
run its retardation feeds are withheld at m4k4. It is a re-plot, so no
independent number is affected. Ledgers: largest per step 3.4e-11 (m1k1) and
8.5e-12 (m4k4) of the equilibrium inventory, cumulative at most 3.7e-9.

## 8. Claim and non-claim

Claim: the printed uptake-benchmark numbers re-scored and the form mapping
re-fitted on the version-2 reading of Cardarelli's figures, with the same
methods, float-identical at the committed inputs, at the reported and the
finest refinement level. Non-claim: no new method, no change to any
coefficient, law or tolerance. Nothing is plant data, calibration or
qualification: one temperature, one radius, points read from a scanned figure.
The Figura 3.9 a_h 0.80 t63 (0.0184 h) sits about 6 px from the axis, so
ratios built on it are resolved to about ±15 % by the scan. Not re-scored
here: the supplement's measured-law seeded comparison ("20 of the 100 points"),
which is another consumer of the same points. `physically_qualifying: false`,
`plant_predictive: false`.

Detail: `DETAIL.md`. Files: `MANIFEST.sha256`.
