# Detail: the F32 collection, its verification, and the full tables (2026-09-25)

Companion to `RESULT.md`. Nothing here is a new claim; this is the method, the
verification, and the tables `RESULT.md` summarizes.

## 1. Reproducing it

Three scripts, run in this order from the repository root under the pinned
environment (`PYTHONHASHSEED=1`, the four thread caps at 1,
`PYTHONDONTWRITEBYTECODE=1`, `PYTHONIOENCODING=utf-8`, `PYTHONPATH=src`,
interpreter `.venv/Scripts/python.exe`):

```
.venv/Scripts/python.exe -B -s <this folder>/collect_f32.py
.venv/Scripts/python.exe -B -s <this folder>/analyze_f32.py
.venv/Scripts/python.exe -B -s <this folder>/make_figure.py
```

`collect_f32.py` reads the extracted farm return (read-only; its location is the
module constant `DEFAULT_PACKAGE`, overridable with the environment variable
`F32_PACKAGE_ROOT`), verifies every artifact against its receipt, places the 48
sweep outputs in `outputs/`, and writes the evidence sidecar. Once `outputs/`
exists the other two need nothing but the tree.

## 2. What was verified at collection

For each of the 48 job folders:

* the folder holds exactly one `F32_faner_sweep_*.json`;
* its sha256 and byte count equal the `artifacts` entry of the job's own
  `receipt.json` (48 of 48 matched; no mismatch, no retry, no repair);
* `farm_stdout.txt` and `farm_stderr.txt` match their receipt digests
  (every `farm_stderr.txt` is empty, sha256 `e3b0c442…` at 0 bytes);
* the receipt's `job.job_id` equals the folder name and the `suffix` is unique
  across the campaign;
* every path the receipt's `export_paths` and `extra_paths` name resolves in the
  working tree (0 missing).

Campaign-level: one commit across all 48
(`bad3513ef877910642e19df8e36e603f940c8710`), one host (TIASERVER), one
interpreter (CPython 3.14.5 / numpy 2.4.6 / scipy 1.18.0), `outcome: done` and
`exit_code: 0` on all 48, `grid_scored_n == grid_n == 20` on all 48, 0 typed
refusals over 960 grid points.

Receipts copied into `docs/evidence/farm_f32_2026-09-25/receipts/` carry one
mechanical substitution and no other change: the directory prefix a farm worker
is told to write its artifact into is replaced by the token
`<worker-artifact-dir>`, because this collection places the outputs in this item
folder instead and a repository-rooted literal with nothing behind it is refused
by `tests/test_repository_link_integrity.py` clause (c). The three path
components the token stands for are recorded, spelled apart, in
`COLLECTION.json`; each job row there carries both `receipt_sha256_original` and
`receipt_sha256_sanitized`, so the untouched receipt is still pinned.

## 3. The score, in full

From `tools/farm/faner_sensitivity_sweep.py`:

```
area_per_kg_dry      = 1 / (rho_env * (1 - eps) * layer_diameters * d_p)
predicted_rate [1/s] = h * (T_gas - T_sat) * area_per_kg_dry / dh_vap
measured_rate  [1/s] = -lsq_slope(t, X) over { (t, X) : X_lo < X <= X_0 }
window_residual      = (predicted_rate - measured_rate) / measured_rate
X_lo                 = 0.1888 + 0.05 = 0.2388 kg/kg
```

`h` is `rig_state.predict(state, "whitaker_4_18")` — Whitaker (1972),
`rig_state.BEST_GROUNDED_INDEPENDENT` — with the Reynolds number recomputed from
the swept velocity and the gas viscosity inverted from the thesis's own printed
Table 4.7, and every other element of the state (voidage, Prandtl number,
specific conductance) interpolated in temperature by `rig_state.py` unchanged.
`T_sat` is `single_component_pore_gas.saturation_temperature(101325 Pa)` =
341.86451 K = 68.7145 °C, and `dh_vap` is the certified hexane equation of
state's, 334 929.467 J/kg. `rho_env` is 1160 (soybean) and 1023 (sunflower)
kg/m³ and `d_p` is 1.80e-3 and 1.95e-3 m, both from `rig_state.py`.

`X_lo` is the frozen identification instrument's critical loading plus its
constant-rate margin, cited by the runner as two float literals rather than
imported. The window is therefore the *same* window for every point of a given
(trace, X₀) pair, which is why the measured rate takes only three distinct values
per trace.

## 4. The per-job table

| job | material | mass (kg) | velocity (m/s) | h (W m⁻² K⁻¹) | Re_eps | scored | residual min | residual max | in signed band | wall (s) | rc |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| F32-001-f32soym015v005 | soybean | 0.015 | 0.05 | 46.601 | 93.0 | 20/20 | −0.7551 | −0.0950 | 1 | 99.065 | 0 |
| F32-002-f32soym015v010 | soybean | 0.015 | 0.10 | 69.615 | 186.1 | 20/20 | −0.6342 | +0.3519 | 1 | 97.371 | 0 |
| F32-003-f32soym015v015 | soybean | 0.015 | 0.15 | 88.174 | 279.1 | 20/20 | −0.5367 | +0.7123 | 3 | 95.583 | 0 |
| F32-004-f32soym015v020 | soybean | 0.015 | 0.20 | 104.343 | 372.2 | 20/20 | −0.4517 | +1.0263 | 0 | 93.545 | 0 |
| F32-005-f32soym015v025 | soybean | 0.015 | 0.25 | 118.946 | 465.2 | 20/20 | −0.3750 | +1.3099 | 4 | 92.177 | 0 |
| F32-006-f32soym015v030 | soybean | 0.015 | 0.30 | 132.417 | 558.3 | 20/20 | −0.3042 | +1.5716 | 0 | 90.579 | 0 |
| F32-007-f32soym030v005 | soybean | 0.030 | 0.05 | 46.601 | 93.0 | 20/20 | −0.7551 | −0.0950 | 1 | 89.038 | 0 |
| F32-008-f32soym030v010 | soybean | 0.030 | 0.10 | 69.615 | 186.1 | 20/20 | −0.6342 | +0.3519 | 1 | 87.509 | 0 |
| F32-009-f32soym030v015 | soybean | 0.030 | 0.15 | 88.174 | 279.1 | 20/20 | −0.5367 | +0.7123 | 3 | 85.976 | 0 |
| F32-010-f32soym030v020 | soybean | 0.030 | 0.20 | 104.343 | 372.2 | 20/20 | −0.4517 | +1.0263 | 0 | 84.000 | 0 |
| F32-011-f32soym030v025 | soybean | 0.030 | 0.25 | 118.946 | 465.2 | 20/20 | −0.3750 | +1.3099 | 4 | 82.625 | 0 |
| F32-012-f32soym030v030 | soybean | 0.030 | 0.30 | 132.417 | 558.3 | 20/20 | −0.3042 | +1.5716 | 0 | 81.285 | 0 |
| F32-013-f32soym045v005 | soybean | 0.045 | 0.05 | 46.601 | 93.0 | 20/20 | −0.7551 | −0.0950 | 1 | 79.955 | 0 |
| F32-014-f32soym045v010 | soybean | 0.045 | 0.10 | 69.615 | 186.1 | 20/20 | −0.6342 | +0.3519 | 1 | 78.564 | 0 |
| F32-015-f32soym045v015 | soybean | 0.045 | 0.15 | 88.174 | 279.1 | 20/20 | −0.5367 | +0.7123 | 3 | 76.913 | 0 |
| F32-016-f32soym045v020 | soybean | 0.045 | 0.20 | 104.343 | 372.2 | 20/20 | −0.4517 | +1.0263 | 0 | 75.369 | 0 |
| F32-017-f32soym045v025 | soybean | 0.045 | 0.25 | 118.946 | 465.2 | 20/20 | −0.3750 | +1.3099 | 4 | 73.705 | 0 |
| F32-018-f32soym045v030 | soybean | 0.045 | 0.30 | 132.417 | 558.3 | 20/20 | −0.3042 | +1.5716 | 0 | 72.378 | 0 |
| F32-019-f32soym060v005 | soybean | 0.060 | 0.05 | 46.601 | 93.0 | 20/20 | −0.7551 | −0.0950 | 1 | 70.788 | 0 |
| F32-020-f32soym060v010 | soybean | 0.060 | 0.10 | 69.615 | 186.1 | 20/20 | −0.6342 | +0.3519 | 1 | 69.354 | 0 |
| F32-021-f32soym060v015 | soybean | 0.060 | 0.15 | 88.174 | 279.1 | 20/20 | −0.5367 | +0.7123 | 3 | 67.783 | 0 |
| F32-022-f32soym060v020 | soybean | 0.060 | 0.20 | 104.343 | 372.2 | 20/20 | −0.4517 | +1.0263 | 0 | 66.302 | 0 |
| F32-023-f32soym060v025 | soybean | 0.060 | 0.25 | 118.946 | 465.2 | 20/20 | −0.3750 | +1.3099 | 4 | 64.922 | 0 |
| F32-024-f32soym060v030 | soybean | 0.060 | 0.30 | 132.417 | 558.3 | 20/20 | −0.3042 | +1.5716 | 0 | 63.665 | 0 |
| F32-025-f32sunm015v005 | sunflower | 0.015 | 0.05 | 48.945 | 98.9 | 20/20 | −0.7840 | −0.3178 | 0 | 62.050 | 0 |
| F32-026-f32sunm015v010 | sunflower | 0.015 | 0.10 | 73.138 | 197.7 | 20/20 | −0.6772 | +0.0194 | 2 | 60.447 | 0 |
| F32-027-f32sunm015v015 | sunflower | 0.015 | 0.15 | 92.651 | 296.6 | 20/20 | −0.5910 | +0.2914 | 0 | 58.844 | 0 |
| F32-028-f32sunm015v020 | sunflower | 0.015 | 0.20 | 109.655 | 395.5 | 20/20 | −0.5160 | +0.5284 | 2 | 56.787 | 0 |
| F32-029-f32sunm015v025 | sunflower | 0.015 | 0.25 | 125.014 | 494.3 | 20/20 | −0.4482 | +0.7425 | 0 | 55.141 | 0 |
| F32-030-f32sunm015v030 | sunflower | 0.015 | 0.30 | 139.182 | 593.2 | 20/20 | −0.3856 | +0.9400 | 4 | 53.760 | 0 |
| F32-031-f32sunm030v005 | sunflower | 0.030 | 0.05 | 48.945 | 98.9 | 20/20 | −0.7840 | −0.3178 | 0 | 52.441 | 0 |
| F32-032-f32sunm030v010 | sunflower | 0.030 | 0.10 | 73.138 | 197.7 | 20/20 | −0.6772 | +0.0194 | 2 | 50.671 | 0 |
| F32-033-f32sunm030v015 | sunflower | 0.030 | 0.15 | 92.651 | 296.6 | 20/20 | −0.5910 | +0.2914 | 0 | 44.703 | 0 |
| F32-034-f32sunm030v020 | sunflower | 0.030 | 0.20 | 109.655 | 395.5 | 20/20 | −0.5160 | +0.5284 | 2 | 43.288 | 0 |
| F32-035-f32sunm030v025 | sunflower | 0.030 | 0.25 | 125.014 | 494.3 | 20/20 | −0.4482 | +0.7425 | 0 | 41.429 | 0 |
| F32-036-f32sunm030v030 | sunflower | 0.030 | 0.30 | 139.182 | 593.2 | 20/20 | −0.3856 | +0.9400 | 4 | 39.981 | 0 |
| F32-037-f32sunm045v005 | sunflower | 0.045 | 0.05 | 48.945 | 98.9 | 20/20 | −0.7840 | −0.3178 | 0 | 38.733 | 0 |
| F32-038-f32sunm045v010 | sunflower | 0.045 | 0.10 | 73.138 | 197.7 | 20/20 | −0.6772 | +0.0194 | 2 | 37.553 | 0 |
| F32-039-f32sunm045v015 | sunflower | 0.045 | 0.15 | 92.651 | 296.6 | 20/20 | −0.5910 | +0.2914 | 0 | 36.041 | 0 |
| F32-040-f32sunm045v020 | sunflower | 0.045 | 0.20 | 109.655 | 395.5 | 20/20 | −0.5160 | +0.5284 | 2 | 34.493 | 0 |
| F32-041-f32sunm045v025 | sunflower | 0.045 | 0.25 | 125.014 | 494.3 | 20/20 | −0.4482 | +0.7425 | 0 | 32.917 | 0 |
| F32-042-f32sunm045v030 | sunflower | 0.045 | 0.30 | 139.182 | 593.2 | 20/20 | −0.3856 | +0.9400 | 4 | 30.844 | 0 |
| F32-043-f32sunm060v005 | sunflower | 0.060 | 0.05 | 48.945 | 98.9 | 20/20 | −0.7840 | −0.3178 | 0 | 29.413 | 0 |
| F32-044-f32sunm060v010 | sunflower | 0.060 | 0.10 | 73.138 | 197.7 | 20/20 | −0.6772 | +0.0194 | 2 | 27.956 | 0 |
| F32-045-f32sunm060v015 | sunflower | 0.060 | 0.15 | 92.651 | 296.6 | 20/20 | −0.5910 | +0.2914 | 0 | 26.541 | 0 |
| F32-046-f32sunm060v020 | sunflower | 0.060 | 0.20 | 109.655 | 395.5 | 20/20 | −0.5160 | +0.5284 | 2 | 25.134 | 0 |
| F32-047-f32sunm060v025 | sunflower | 0.060 | 0.25 | 125.014 | 494.3 | 20/20 | −0.4482 | +0.7425 | 0 | 23.771 | 0 |
| F32-048-f32sunm060v030 | sunflower | 0.060 | 0.30 | 139.182 | 593.2 | 20/20 | −0.3856 | +0.9400 | 4 | 22.144 | 0 |

The four rows at each (material, velocity) are identical in every score column,
which is the degeneracy of section 3 of `RESULT.md` read off the raw table. The
wall times fall monotonically with the job index because all 48 were admitted to
a 32-processor host together and finished within twelve seconds of each other;
the per-job compute is under two seconds, as the package README's own local
smoke run measured.

## 5. The full residual surface

Window residual at the source's printed **lower** initial loading, 0.45 kg/kg.
Rows are gas velocity in m/s, columns bed depth in particle diameters. Cells
inside the signed acceptance band (−0.117872 to −0.018422) are marked `*`;
cells inside the printed depth range (3 to 5) are the middle three columns.

**soybean, 120 °C**

| v \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| 0.05 | −0.2654 | −0.5103 | −0.6327 | −0.7062 | −0.7551 |
| 0.10 | +0.0974 | −0.2684 | −0.4513 | −0.5610 | −0.6342 |
| 0.15 | +0.3900 | −0.0734 `*` | −0.3050 | −0.4440 | −0.5367 |
| 0.20 | +0.6448 | +0.0966 | −0.1776 | −0.3421 | −0.4517 |
| 0.25 | +0.8751 | +0.2500 | −0.0625 `*` | −0.2500 | −0.3750 |
| 0.30 | +1.0874 | +0.3916 | +0.0437 | −0.1650 | −0.3042 |

**sunflower, 100 °C**

| v \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| 0.05 | −0.3332 | −0.5554 | −0.6666 | −0.7333 | −0.7777 |
| 0.10 | −0.0035 | −0.3357 | −0.5018 | −0.6014 | −0.6678 |
| 0.15 | +0.2623 | −0.1584 | −0.3688 | −0.4951 | −0.5792 |
| 0.20 | +0.4940 | −0.0040 | −0.2530 | −0.4024 | −0.5020 |
| 0.25 | +0.7032 | +0.1355 | −0.1484 | −0.3187 | −0.4323 |
| 0.30 | +0.8963 | +0.2642 | −0.0519 `*` | −0.2415 | −0.3679 |

Window residual at the printed **upper** initial loading, 0.50 kg/kg (identical
to 0.55, since no digitized point lies above 0.50):

**soybean, 120 °C**

| v \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| 0.05 | −0.2523 | −0.5015 | −0.6262 | −0.7009 | −0.7508 |
| 0.10 | +0.1169 | −0.2554 | −0.4415 | −0.5532 | −0.6277 |
| 0.15 | +0.4147 | −0.0569 `*` | −0.2926 | −0.4341 | −0.5284 |
| 0.20 | +0.6741 | +0.1161 | −0.1629 | −0.3303 | −0.4420 |
| 0.25 | +0.9084 | +0.2723 | −0.0458 `*` | −0.2366 | −0.3639 |
| 0.30 | +1.1246 | +0.4164 | +0.0623 | −0.1502 | −0.2918 |

**sunflower, 100 °C**

| v \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| 0.05 | −0.3519 | −0.5679 | −0.6759 | −0.7407 | −0.7840 |
| 0.10 | −0.0315 `*` | −0.3543 | −0.5157 | −0.6126 | −0.6772 |
| 0.15 | +0.2269 | −0.1821 | −0.3865 | −0.5092 | −0.5910 |
| 0.20 | +0.4521 | −0.0320 `*` | −0.2740 | −0.4192 | −0.5160 |
| 0.25 | +0.6554 | +0.1036 | −0.1723 | −0.3378 | −0.4482 |
| 0.30 | +0.8431 | +0.2287 | −0.0785 `*` | −0.2628 | −0.3856 |

The X₀ = 0.40 surface, which discards measured points, is in
`f32_trace_surface.csv` with the rest; the two tables above are the
printed-admissible loadings only.

Read down a column: the residual rises monotonically with velocity at every
depth, on both traces, and crosses zero once. Read across a row: it falls
monotonically with depth. That single crossing per line is why each trace meets
the band somewhere, and why that alone is not evidence of anything.

## 6. The implied footprint against the printed holder

`implied_footprint_area_m2 = area_per_kg_dry × sample_mass`, which is the bed
volume divided by the layer depth — the footprint the swept (mass, depth) pair
implies. The printed holder cross-section is 0.011 m² (both sources, Sec. 2.2
of the 2019 article; `rig_state.BED_AREA_M2`).

**soybean** (m²)

| mass (kg) \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| 0.015 | 0.01142 | 0.00761 | 0.00571 | 0.00457 | 0.00381 |
| 0.030 | 0.02283 | 0.01522 | **0.01142** | 0.00913 | 0.00761 |
| 0.045 | 0.03425 | 0.02283 | 0.01713 | 0.01370 | 0.01142 |
| 0.060 | 0.04567 | 0.03045 | 0.02283 | 0.01827 | 0.01522 |

**sunflower** (m²)

| mass (kg) \ depth | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| 0.015 | 0.01061 | 0.00707 | 0.00531 | 0.00424 | 0.00354 |
| 0.030 | 0.02122 | 0.01415 | **0.01061** | 0.00849 | 0.00707 |
| 0.045 | 0.03184 | 0.02122 | 0.01592 | 0.01273 | 0.01061 |
| 0.060 | 0.04245 | 0.02830 | 0.02122 | 0.01698 | 0.01415 |

Only the ratio mass/depth enters, so the diagonal `(0.015, 2)`, `(0.030, 4)`,
`(0.045, 6)` is one footprint. The article's own approximate 30 g at four
particle diameters — the middle of its own printed three-to-five range — gives
1.038 and 0.965 times the printed holder area. This is a consistency check on
the sweep's geometry, not an identification of the run: any of the three diagonal
cells reproduces the same footprint, and the article prints its mass as
approximate.

## 7. Where the item-5 march sits

| trace | item-5 Whitaker h (W m⁻² K⁻¹) | bracketing sweep velocities (m/s) | bracketing sweep h | interpolated velocity (m/s) | thesis printed approach velocity (m/s) |
|---|---:|---|---|---:|---|
| soybean 120 °C | 78.649 | 0.10, 0.15 | 69.615, 88.174 | **0.1243** | 0.117 to 0.129 |
| sunflower 100 °C | 80.013 | 0.10, 0.15 | 73.138, 92.651 | **0.1176** | 0.117 to 0.129 |

The interpolation is linear in the sweep's own two tabulated coefficients and is
labelled as such; it is a placement of the item-5 march on this axis, not a
measurement of the coefficient at that velocity. The item-5 transfer area on the
charge reading, 0.5466 and 0.5497 m²/kg dry, corresponds in the sweep's depth
convention to 2.785 and 2.574 particle diameters, which is below the printed
three-to-five range; the two areas are built differently (item 5 recovers the
holder from the source's own Eq. 4.25, the sweep imposes a depth) and the
difference is reported rather than reconciled.

## 8. Caveats carried forward

* **The 0.40 kg/kg loading column discards measured points** and, on the soybean
  trace, leaves a two-point window — one point above the runner's own refusal
  threshold, a secant rather than a regression. The soybean best-overall
  coordinate sits there. The printed-box row of `RESULT.md` §4, which excludes
  it, is the one to quote.
* **The transfer-area convention item 5 left open is not reopened here.** The
  sweep uses one construction (depth × envelope density × solid fraction)
  throughout; the charge-versus-dry-meal question of item 5 §M3 is orthogonal to
  it and unaffected by this record.
* **Whitaker's printed Reynolds bound is a disclosure, not a clamp**, and the
  sweep evaluates at the voidage Reynolds number while the thesis's own
  attribution table prints the bound without the voidage subscript. That
  disclosure is carried by the runner in every output and is unresolved here, as
  in item 5.
* **Nothing about the falling-rate leg is measured by this campaign.** F32 scores
  the constant-rate leg only. The falling-rate reading of the journal traces,
  and the thirtyfold shell-coefficient gap the supplement's §S10.9 reports, are
  untouched.
* **The sweep does not reach beyond 0.30 m/s.** Both jointly admitted
  coordinates inside the printed box sit at or near that edge, so whether a
  still higher velocity would tighten or reverse the sign pattern is not
  measured. That would be a new campaign, not a re-reading of this one.
