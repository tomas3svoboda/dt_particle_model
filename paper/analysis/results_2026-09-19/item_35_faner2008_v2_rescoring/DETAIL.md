# Item 35: detail

## 1. Pipeline (run in this order; every step a hard stop on an identity failure)

Environment on every call: PYTHONHASHSEED, OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS = 1,
PYTHONDONTWRITEBYTECODE = 1 (enforced by `v2common.py`, so the committed folders
receive no bytecode), PYTHONIOENCODING = utf-8, PYTHONPATH = `src`, interpreter
`.venv\Scripts\python.exe`, working directory this folder.

| step | script | what it does | log |
|---|---|---|---|
| 0 | `remarch_item27.py` | re-runs the 12 item 27 thesis cells whose per-step files were not kept, through `run_item27.run_march` (imported read only, output redirected to `remarch/`); accepts only byte-identical per-step CSVs and exact summaries | `log_remarch.txt` |
| 1 | `rescore_runs.py` | rebuilds every thesis run's decimated and full-step trajectory from its per-step record, reproduces both committed scorings at v1, scores at v2 (both windows) | `log_pipeline.txt`, `outputs/rescore_identity.json`, `outputs/rescore_runs.csv` |
| 2 | `item05_measured_side.py` | item 5's measured side (slope, crossing, plateau, duty, demanded h, film-case ratios, circularity probe, 88 stored tail traces, identified-scalar traces) | `outputs/item05_*.{json,csv}` |
| 3 | `window_quantities.py` | the leverage / factor / in-range statistics of the surface record and of item 28 over the window, v1 identity then v2 | `outputs/window_quantities.json` |
| 4 | `rescore_spaceresolved.py` | items 5-dxr and 32 re-scored with `dxr_common.window_score` | `outputs/spaceresolved_*` |
| 5 | `shadow_analyses.py` | copies the results tree to a temporary directory, replaces the measured points, bands and window start (rig files, surface-record window) and the run scores, runs the committed `analyze_item24/27/29/30/31.py` unchanged; mode v1 must be byte-identical to the committed outputs | `outputs/shadow_{v1,v2_rule,v2_p15}/`, `outputs/log_shadow_*.txt` |
| 6 | `flatten_diff.py` | every analysis.json leaf, v1 / v2 rule / v2 15-point | `outputs/analysis_leaves_v1_v2.csv` (9,735 leaves; 2,536 change under the rule) |
| 7 | `build_printed_table.py` | the printed-number table of RESULT section 6 | `outputs/printed_numbers_v1_v2.{csv,md}` |

`log_pipeline.txt` is the console of steps 1 to 7 in one pass (12:55 to 12:58 local);
step 0 ran 12:27 to 12:46. No step sleeps; no wall time includes a host sleep.

## 2. Window variants

`v2common.window_start`: v1 and `v2_rule` use the measured crossing of X = 0.20
(linear interpolation of the points), exactly as `run_rerun_faner.py` computes it and
every scorer applies it (t >= start). `v2_p15` starts the window at the first sample
at or after the version-1 crossing: 60.0408676 s sunflower (15 points), 60.030346 s
soybean (the same 15 points as the rule). In `v2_p15` the run JSONs' `crossing_ratio`
is still t_c over the version-2 crossing; the crossing ratios in RESULT come from the
rule mode only.

| window | sunflower points / span kg/kg / start s | soybean points / span / start |
|---|---|---|
| v1 | 15 / 0.087668 / 57.699 | 15 / 0.079121 / 54.930 |
| v2 rule | **14 / 0.069746 / 60.144** | 15 / 0.091376 / 58.034 |
| v2 15-point | 15 / 0.096093 / 60.041 | 15 / 0.091376 / 60.030 |

## 3. The declared-sphere ladder, every cell (item 31 runs re-scored)

As printed (decimated) / full-step, in band in brackets.

| trace | cells/level | v1 | v2 rule | v2 15-point |
|---|---|---|---|---|
| sunflower | 60/1 | 0.1330(10) / 0.1325(11) | 0.1721(9) / 0.1703(10) | 0.1257(9) / 0.1264(10) |
| sunflower | 60/2, 60/4 | 0.1320, 0.1315 | 0.1695, 0.1687 | 0.1260, 0.1258 |
| sunflower | 120/1, 120/2, 120/4 | 0.1323, 0.1314, 0.1310 | 0.1711, 0.1686, 0.1678 | 0.1251, 0.1257, 0.1255 |
| sunflower | 240/1, 240/2, 240/4 | 0.1321, 0.1313, 0.1309 | 0.1708, 0.1684, 0.1676 | 0.1250, 0.1256, 0.1254 |
| sunflower | 240/4 full-step | 0.1308 | 0.1675 | 0.1254 |
| soybean | 60/1 | 0.1971(6) / 0.1931(6) | 0.1562(8) / 0.1545(8) | same as rule |
| soybean | 60/2, 60/4 | 0.1923, 0.1913 | 0.1539, 0.1532 | same |
| soybean | 120/1, 120/2, 120/4 | 0.1962, 0.1913, 0.1902 | 0.1554, 0.1533, 0.1525 | same |
| soybean | 240/1, 240/2, 240/4 | 0.1959, 0.1911, 0.1900 | 0.1552, 0.1531, 0.1524 | same |
| soybean | 240/4 full-step | 0.1898 | 0.1523 | same |

Every other cell of every item: `outputs/rescore_runs.csv` (columns
`<mode>_<dec|full>_<rms_over_span|in_band|mean_signed_over_span|rms_kg_kg>`,
window points, span and start per mode, and the sha256 of the per-step file scored).
Every derived quantity: `outputs/shadow_<mode>/<item>/` (the committed analysis
scripts' own tables) and `outputs/analysis_leaves_v1_v2.csv`.

## 4. Rounding and attribution notes on the printed-number table

* "1.8 to 2.6" (Bird and Bradshaw-Myers over the demand): the recomputed v1 maximum
  is 2.548 (2008 sunflower, Bradshaw-Myers), which the paper prints as 2.6; v2 2.743.
* "17.228" (R10) and "16.18" (N7) and "1671.4" (F20) and "129.0" (dry demand, low)
  belong to Faner 2019 conditions and are unchanged; the 2008 values are shown in
  their rows.
* "0.29-0.74" (N16) is the surface record's own statement: the low end is reached
  on a dry-area closure (0.287), the high end over the ten particle-mass closures;
  the row is computed the same way and is unchanged. "3.17-3.18" (B5/N20) is the six
  tabulated rows; over all ten particle-mass closures 3.165 to 3.196 (v1) and 3.159
  to 3.196 (v2 rule).
* "0.736" (S1): recomputed 0.7355 (decimated), printed rounded up.
* I53/I54 (the D/R^2 identity): all four differences stay below 1e-13 of the span at
  v2 (1.5e-14, 5e-16, 4.6e-14, 1.1e-16); the row lists them in the table's order.
* I67/L9 (sphericity reproducing item 30's R_d value): linear interpolation of item
  29's as-printed radius sweep on the branch above its minimum, psi = rf x 0.885 /
  0.975; 0.67 at v1 and at v2 (where the R_d value itself is 9.9, not 7.9).
* K15-K19: item 27's own band (F = 3.165) is kept as recorded; the paper's
  0.885 mm range combines item 27's interior sweep minimum (sunflower) with item 28's
  edges; the "tab:fanermarch" row reproduces that combination.

## 5. Cells not re-scored

44 run JSONs keep their version-1 scores inside the shadow tree because no per-step
record exists for them (`outputs/shadow_v2_rule/STALE_INPUTS.json`): item 24's 32
identity, 480/960-cell and band_fine cells and item 27's six outer-read identity
cells. None is behind a number the paper prints (item 24's printed 19.7 / 30.8 are its
240/4 cells, which have per-step records and were re-scored); the item 24 continuum
and extension entries of `shadow_v2_*/item_24_*` that depend on them are not used.

## 6. Files

* Scripts: `v2common.py`, `remarch_item27.py`, `rescore_runs.py`,
  `item05_measured_side.py`, `window_quantities.py`, `rescore_spaceresolved.py`,
  `shadow_analyses.py`, `flatten_diff.py`, `build_printed_table.py`.
* `remarch/item_27_boundary_value_surface_read/runs/`: the 12 re-marched cells
  (summary JSON, per-step CSV, `remarch_log.jsonl`).
* `outputs/`: every table named above.
* `MANIFEST.sha256`: sha256 of every file in this folder except itself.
