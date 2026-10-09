# Item 37: the journal-condition comparisons and the floor statements re-scored against the version-2 Faner et al. (2019) digitization (2026-09-27)

Inputs, read only: `paper/analysis/datasets/faner2019_v2/` (DATASET_RECORD.md,
curves.csv, adjudication_per_point.csv), item_05, item_22, item_23, item_29,
item_30, item_31 and item_34. Nothing was fitted. Nothing outside this folder
was written: no paper file, nothing under `docs/`, no other result item.
Every number below was printed by a script in this folder (see DETAIL.md).

`fitted: false. physically_qualifying: false. plant_predictive: false. qsc10_complete: false.`

## 1. item_22: the F32 sweep re-scored against the version-2 Figure 2 traces

**What changes and what does not.** item_22 scores only the two Figure 2 loading
traces. It does not use Figures 3 and 4, so the corrected series identities
there do not reach it. The runner (`tools/farm/faner_sensitivity_sweep.py`)
predicts the constant rate from the job's film coefficient, the fixed gas
temperature (100 or 120 °C), the saturation temperature, the latent heat and
the layer geometry. The initial-loading axis is the fixed grid 0.40 to 0.55. It
is not read from the curve. The curve enters only the measured slope over
(0.2388, X0_grid]. None of the runner's inputs changed, so nothing was
re-marched. The 48 stored reports were re-scored with the runner's own
`score_grid_point`.

**Identity first.**
- Re-scoring with the version-1 curve reproduced all 48 stored reports: 11,424
  leaves, 0 mismatches.
- item_22's `analyze_f32.py` on those outputs reproduced its committed
  `f32_summary.json`: 503 leaves, 0 mismatches.
- The acceptance band is item_05's thesis-run band and does not depend on this
  dataset. It is unchanged: signed −0.1179 to −0.0184, magnitude cap 0.1179.

**Which points moved (3 of 50).**
- Soybean t = 0: 0.4906 → 0.4924.
- Soybean 9.884 s: 0.4264 → 0.4209.
- Sunflower t = 0: 0.4993 → 0.5016.

The other changed points lie below the window's lower bound of 0.2388 and no
window can see them: soybean 40.154 s (0.2250) and both 240 s final loadings.

**The effect of the sunflower t = 0 point.** It now reads 0.5016, which is above
the X0 = 0.50 grid value, so the 0.50 window drops t = 0. It goes from 6 points
to 5 and becomes identical to the 0.45 window. The 0.55 window keeps all 6
points.

**Measured window rates (1/s), v1 → v2:**

| trace | X0 = 0.40 | 0.45 | 0.50 | 0.55 |
|---|---|---|---|---|
| soybean | 6.0016e-3 (unchanged) | 7.3936e-3 → 7.1200e-3 (−3.7 %) | 7.2643e-3 → 7.2606e-3 | as 0.50 |
| sunflower | 4.7411e-3 → 4.7920e-3 | 4.8504e-3 → 4.8840e-3 | 4.9904e-3 → 4.8840e-3 (−2.1 %, t = 0 dropped) | 4.9904e-3 → 5.0481e-3 |

**Residual grid.**
- Largest shift of any grid residual: 0.080 on soybean (90 of 120 rows move)
  and 0.040 on sunflower (all 120 rows move).
- Range, soybean: −0.7551 to 1.5716 → −0.7506 to 1.5716.
- Range, sunflower: −0.7840 to 0.9400 → −0.7864 to 0.9194.

| | before (v1, committed) | after (v2) |
|---|---|---|
| soybean, best in whole sweep (v, X0, L) | (0.20, 0.40, 4) +0.0132 | unchanged, (0.20, 0.40, 4) +0.0132 |
| soybean, best in printed box | (0.30, 0.45, 4) **+0.0437** | (0.25, 0.45, 4) **−0.0264**, now in the signed band |
| sunflower, best in whole sweep | (0.10, 0.45, 2) −0.0035 | (0.20, 0.40, 3) +0.0081 |
| sunflower, best in printed box | (0.20, 0.45, 3) **−0.0040** | same coordinate, **−0.0109** |
| soybean in signed band / within cap (box) | 9 / 20 (4 / 8) | 9 / 18 (4 / 7) |
| sunflower in signed band / within cap (box) | 8 / 14 (3 / 5) | 6 / 13 (2 / 4) |

**Band width per lever** (median residual span across the lever, other levers
fixed), v1 → v2:
- Soybean: gas velocity 0.688 → 0.702; initial loading 0.170 → 0.157; layer
  depth 1.116 → 1.117; sample mass 0 → 0.
- Sunflower: gas velocity 0.615 → 0.610; initial loading 0.034 → 0.034; layer
  depth 0.968 → 0.957; sample mass 0 → 0.

The mass span stays exactly zero, as item_22 found. The ordering of the levers
is unchanged: layer depth, then velocity, then initial loading.

**Does one coordinate serve both traces?**
- On the signed band: **no**, before or after. No coordinate lies in the band on
  both traces, in the box or in the whole sweep.
- On the magnitude cap in the printed box: before, four common coordinates;
  after, three: (0.20, 0.50, 3), (0.30, 0.45, 4) and (0.30, 0.50, 4). The lost
  one is (0.20, 0.45, 3), where the soybean residual moves 0.0966 → 0.1387,
  beyond the cap.
- On the whole-sweep magnitude cap: 9 → 7 common coordinates.
- The three box coordinates that remain still carry opposite signs:
  - soybean +0.117, +0.084 and +0.063;
  - sunflower −0.011, −0.058 and −0.058.
  item_22's statement stands unchanged: the two traces do not share one
  condition within the signed band.

**Flips of band membership.**
- Sunflower (0.10, 0.50, 2) and (0.20, 0.50, 3) leave the signed band (−0.031
  → −0.010). The cause is the t = 0 point dropped from the 0.50 window.
- Soybean (0.10, 0.45, 2) and (0.20, 0.45, 3) and sunflower (0.25, 0.50, 3)
  leave the magnitude cap.
- No row enters either set.

## 2. item_23 and item_05: every ratio at the final loading 7.415e-3, restated as an interval

**Version-2 final loadings** (DATASET_RECORD §3.2):
- Sunflower: 0.0105, flag inferred, plausible range 0.0075 to 0.0125.
- Soybean: 0.0080 ± 0.006, i.e. 0.0020 to 0.0140, flag merged.

The two traces no longer share one value. Floors are item_23's committed
`item23_floors.csv`. The version-1 ratios were recomputed first and matched
item_23's columns exactly.

**Ratios.** Central value first, then (low end of the range to high end).

| trace | final loading m | m / arm-free floor | m / floor with arm (declared Ω_w) | log-gap closure |
|---|---|---|---|---|
| sunflower 100 °C, before | 0.007415 | 2.79 | 1.77 | 44.7 % |
| sunflower 100 °C, after | 0.0105 (0.0075 to 0.0125) | **3.96 (2.83 to 4.71)** | **2.50 (1.79 to 2.98)** | **33.4 % (44.2 to 29.6)** |
| soybean 120 °C, before | 0.007415 | 6.09 | 3.35 | 33.2 % |
| soybean 120 °C, after | 0.0080 (0.0020 to 0.0140) | **6.57 (1.64 to 11.50)** | **3.61 (0.90 to 6.32)** | **31.8 % (121 to 24.5)** |

- **The coefficient bracket (Ω_w low to high).** Sunflower m / F_arm: 2.37 to
  2.62 at the central value, 1.70 to 3.12 over the range. Soybean: 3.40 to 3.82
  at the central value, 0.85 to 6.69 over the range.
- **item_23's "33 and 45 per cent".** This becomes 32 per cent (soybean, 24.5
  to above 100) and 33 per cent (sunflower, 29.6 to 44.2). The 45 per cent is
  now the low-loading end of the sunflower interval, not its centre.
- **The measured-VLE arm (sunflower).** m / floor goes from 1.50 to 2.13 (1.52
  to 2.54). The closure goes from 60.2 % to 45.0 % (59.6 to 39.9).
- **item_23 R3 ("the floor with the arm exceeds the measured loading").**
  - Sunflower: clear over the whole range. The largest floor, 4.42e-3, is below
    0.0075.
  - Soybean: **not clear at the low end.** At 0.0020 the floor with the arm
    (2.22e-3, bracket 2.09e-3 to 2.35e-3) exceeds the loading, so the closure
    exceeds 100 %. The soybean ratio is therefore not bounded away from 1 by
    this reading. It is bounded away from 1 only for m > 2.35e-3.
- **Why item_23 excluded the mole basis.** Its argument was that the mole-basis
  floor exceeds the measured loading.
  - Sunflower: that floor is 9.86e-3. It exceeds only the low end (0.0075), not
    the central 0.0105 or the high 0.0125.
  - Soybean: that floor is 3.10e-3. It exceeds only the low end (0.0020).
  - The mole-basis exclusion therefore no longer follows from this dataset
    alone. It holds only at the low end of each interval, and the exclusion
    needs its other grounds.
- **item_05's own older arm (its 2026-09-21 floors).**
  - m / sorbate floor: sunflower 2.79 → 3.96 (2.83 to 4.71); soybean 6.09 →
    6.57 (1.64 to 11.50).
  - m / corrected floor (ideal): sunflower 1.91 → 2.70 (1.93 to 3.22); soybean
    3.98 → 4.29 (1.07 to 7.51).
  - m / corrected floor (low): sunflower 2.23 → 3.16 (2.26 to 3.77); soybean
    4.62 → 4.99 (1.25 to 8.73).

**item_05's journal-condition inputs, re-derived with item_05's own rules**
(v1 identity against `rerun_faner.json`: 0 mismatches at 1e-14 relative):

| quantity | sunflower 100 °C, v1 → v2 | soybean 120 °C, v1 → v2 |
|---|---|---|
| X0 (first curve point) | 0.499279 → 0.50164 | 0.490628 → 0.492364 |
| last loading | 0.007415 → 0.0105 | 0.007415 → 0.0080 |
| constant rate over X > 0.2388 (1/s) | 4.9904e-3 → 5.0481e-3 (6 pts) | 7.2643e-3 → 7.2606e-3 (4 pts) |
| X = 0.20 crossing (s) | 65.08 → 65.06 | 45.46 → 45.86 |
| Tp plateau mean to crossing (°C) | 68.63 → 68.51 | 68.81 → 69.11 |
| Tg mean to crossing (°C) | 101.82 → 101.55 | 120.26 → 119.44 |
| h at constant rate, bed area (W/m²K) | 91.61 → 92.93 | 86.52 → 88.31 |
| Whitaker predicted / measured | 0.823 → 0.815 | 0.906 → 0.908 |
| Whitaker constant-rate duration (s) | 73.00 → 73.46 | 44.25 → 44.46 |

item_05's plateau and gas means include the t = 0 marker. The soybean t = 0
readings moved on the axis (Tp 68.68 → 69.74 °C, Tg 119.57 → 117.87 °C). The
record's §9 values (68.91 and about −0.5 °C) leave that marker out, which is
why they differ from these.

## 3. Floors at published compositions (from item_34) and the replacement sentence

**How the table is built.** Total storage uses item_06's construction. The arm
is item_23's weight-basis oil arm at w_o = 0.0195 and Ω_w = 5.2. The class is
the published 100 to 500 ppm. The floors are in ppm of dry meal; y is the
n-hexane mole fraction.

| composition | source | T (°C) | a_h | floor, no arm | floor, with arm | vs class |
|---|---|---:|---:|---:|---:|---|
| sparged-section entry vapour, 0.03 kg/kg (y 0.0064) | Cardarelli 1998, p. 149 | 105 | 0.00233 | 11.6 | 20.3 | below |
| same, at 99 °C | Cardarelli 1998, p. 149 | 99 | 0.00271 | 14.7 | 24.8 | below |
| Fig. 5.3 vapour, y 0.0038 to 0.0149 | Cardarelli 1998, Fig. 5.3 | 105 | 0.0014 to 0.0054 | 6.8 to 26.8 | 12.0 to 47.1 | below |
| model vapour, 10.3 to ~11 g/kg (y 0.0022 to 0.0023) | Cardarelli et al. 2002, Fig. 2 | 105 | 0.00079 to 0.00084 | 3.9 to 4.2 | 6.9 to 7.3 | below |
| thesis model vapour (y 0.0021) | Cardarelli 1998, Fig. 5.12 | 105 | 0.00078 | 3.9 | 6.8 | below |
| last trays, y 0.003 to 0.03 (discharge statement) | the paper's declared grid | 100 to 110 | 0.0010 to 0.0123 | 4.5 to 66.0 | 8.1 to 112.4 | below; with arm, up to inside (0.03 at 100 °C only) |
| last trays, y 0.1 | the paper's declared grid | 100 to 110 | 0.032 to 0.041 | 152 to 222 | 273 to 378 | inside |
| vessel-top exit vapour, y 0.612, at its own 75 °C | Cardarelli 1998, Fig. 5.1 | 75 | 0.504 | 4834 | 6928 | above; not a gas the dry meal meets |
| same paired with 105 °C (the paper's 1216.4 / 2085.5) | withdrawn pairing | 105 | 0.222 | 1216 | 2085 | above; mis-paired, withdrawn |

- **Every published composition of the gas the sparged meal meets** gives 3.9 to
  26.8 ppm without the arm and 6.8 to 47.1 ppm with it. All are below the class.
- **Replacement for the vessel-gas sentence.** This is item_34's sentence
  verbatim:

  > At the vapour composition Cardarelli (1998, p. 149) gives for the gas at the
  > entry of the sparged section, 0.03 kg of n-hexane per kg of vapour (mole
  > fraction 0.0064) at 105 °C, the same construction floors at 11.6 ppm without
  > the oil-solution arm and 20.3 ppm with it, and across that source's Fig. 5.3
  > (0.004 to 0.015) at 6.8 to 47.1 ppm; at the model vapour of Cardarelli et
  > al. (2002), at most about 10 g per kg (0.0022), it floors at 3.9 to 6.9 ppm.
  > All lie below the published class, so equilibrium does not bind at any
  > published composition of the gas the meal meets in the sparged section. The
  > 88.3 per cent n-hexane printed for the vessel's exit vapour (Cardarelli
  > 1998, Fig. 5.1) is the gas leaving above the pre-desolventizing stages at
  > 75 °C and bears on no dry-meal floor.

- **The discharge statement.** Add this after it:

  > At the last trays of the declared grid, 100 to 110 °C, a vapour of n-hexane
  > mole fraction 0.003 to 0.03 floors at 4.5 to 66.0 ppm without the arm and
  > 8.1 to 112.4 ppm with it, reaching the class only with the arm at 0.03 and
  > 100 °C; a vapour of 0.1 floors inside the class, at 152 to 222 ppm without
  > the arm and 273 to 378 ppm with it.

- **Withdrawn everywhere:** the 1216.4 and 2085.5 ppm values and the phrase
  "the vapour composition reported for the industrial vessel".
- No plant number of this project appears.

## 4. Step 3: the thesis sunflower trace at the 2019 sunflower sphere, 0.670 mm

**The march.** Thesis sunflower trace, 136 °C, boundary-value read, at
R = 0.74 × 1.81 mm / 2 = 0.6697 mm, 60 cells, time level 1. The driver is
item_30's `run` with `rf_of_diameter(1.81)`, the function item_31 drives. The
scoring basis is item_30's thesis points, which keeps the value comparable
with the 0.655 mm column.

**Identity first.** item_30's 0.655 mm cell was re-marched here and matched in
8,529 checks with 0 mismatches.

| sphere | R (mm) | window residual (rms / span) | in band | mean signed |
|---|---:|---:|---:|---:|
| ψ × 2019 soybean d, 1.77 mm (item_30's column, applied to both traces) | 0.6549 | 0.0791 | 15 of 15 | +0.0607 |
| **ψ × 2019 sunflower d, 1.81 mm (this item)** | **0.6697** | **0.0911** | **14 of 15** | +0.0746 |
| ψ × thesis sunflower d, 1.95 mm (item_31) | 0.7215 | 0.133 | 10 of 15 | — |

- **Status:** t_max. 484 steps, 0 retries, and the trace covers the scored
  window.
- **Ledgers:** 3.8e-15 per step, and a total closure of −6.9e-15 relative to
  the initial inventory.
- **Full-step score:** 0.0908 (14 in band).
- **Consistency check:** the value lies inside item_29's radius sweep at
  0.664 mm (0.086) and 0.752 mm (0.156).
- **Corrected column, per species:** "ψ × Faner 2019 sphere" should read
  0.670 mm for sunflower (9.1 per cent, 14 of 15) and 0.655 mm for soybean
  (18.7 per cent, 7 of 15, item_30).
- **Not done here:** the 0.885 mm "volume-equivalent" row is also the 2019
  soybean diameter. The 2019 sunflower volume-equivalent sphere, 0.905 mm, was
  not marched.
- Re-scoring against the version-2 thesis digitization belongs to item_35.

## 5. Step 4: the sphericity statement

**What the article says.** It measured ψ = 0.74 on both of its own meals: "The
average sphericity shape factor (ψ) was 0.74 for both sunflower and soybean
meals" (Faner et al. 2019, p. 5).

**What is no longer assumed:** that 0.74 is borrowed from soybean for
sunflower.

**What remains assumed:** that the 2008 thesis samples share the 2019 meals'
sphericity. The thesis prints no sphericity or flake thickness.

**What the radius sweep says.** The residual depends on ψ only through D/R²
(item_29 and item_31: 0.49 to 0.55 of the window span per unit ln R).
- **Sunflower.** At the thesis diameter, ψ from 0.67 to 1.0 moves the residual
  from 7.9 to 29.1 per cent (item_31 §9).
- **Soybean.** The residual is bracketed by the marched 0.575 and 0.620 mm
  spheres at 11.4 and 15.3 per cent and reaches 37.5 per cent at ψ = 1 (item_29,
  item_30).
- **Rank.** This is the largest single model-form uncertainty on the headline.

**Corrected sentence:**

> The sphericity 0.74 was measured by Faner et al. (2019) on both their
> sunflower and their soybean meal; that the 2008 thesis samples share it is
> assumed, and it is the largest single model-form uncertainty on the
> headline: between ψ = 0.67 and 1.0 the sunflower residual runs from 7.9 to
> 29.1 per cent of the window span, and the soybean residual reaches 37.5 per
> cent at ψ = 1.

**Sentences to replace in the paper** (DATASET_RECORD §9; not edited here):
- "measured on industrial soybean meal";
- "assumed for the sunflower sample";
- "borrows both factors from the soybean meal".

**Also flagged:** the supplement's "Sample mass is logged every 2 s" is
unsupported. The article says only "at second intervals" (p. 2), and the
interval is not printed.

## 6. Claim and non-claim

**Claimed:**
- The item_22 re-scoring on the version-2 Figure 2 points, identity-checked.
- The item_05 inputs re-derived with item_05's rules, identity-checked.
- The item_23 and item_05 ratios as intervals over the version-2 final-loading
  ranges, identity-checked against the committed version-1 columns.
- The floors table, selected from item_34's committed CSV without recomputation.
- The 0.670 mm march, with its identity check.
- The sphericity sentence.

**Not claimed:**
- That the sunflower final loading is known better than 0.0075 to 0.0125. It
  is inferred, not read.
- Any new floor computation.
- A re-scoring of Figures 3 and 4. item_22 does not use them. Their consumer is
  the t50 and t90 validation script.
- A re-scoring against the version-2 thesis digitization (item_35).
- Calibration, validation, qualification or plant prediction.
- Any statement about a plant.
