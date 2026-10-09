# Item 32 — DETAIL

## Code path

* `run_item32.py identity | declared | extension`, `analyze_item32.py`.
* The march is `dxr_common.march_tail_dxr` (item 5, `rerun_2026-09-21c_dxr/`,
  imported read only) calling `dry_shell_variable_coefficient.step` (read only).
  The script replaces `dxr_common.vc` by a `Recorder` proxy whose `step` calls
  the real `step` with the same arguments and returns its result unchanged,
  then records, from public functions of the module (`total_hexane`,
  `cell_coefficients`, `harmonic_face_coefficient`): the inventory before and
  after, the implicit-state Dirichlet-face outflow `dt * A_R * D_face / (R - r_N)
  * (c_N - c_s)` (the discrete flux the residual contains), the per-cell loading
  and coefficients, and the solver strategy. Only accepted strides are
  recorded; refused attempts raise before the recorder.
* Arguments passed that the committed driver does not pass: `radius_m`
  (`rf * 0.885e-3` with `rf = 0.74 * d_thesis / 1.77`, item 30's float
  expression, so the radii equal item 31's), `n_cells` (60, 120), and in the
  extension `dt_max_s` (5/4, 5/16). All are arguments of the committed function.
* Legs. Borrowed: `rerun_2026-09-21/rerun_faner.json`, `film_cases["1.0"]`
  (h = 131.34 W/m2K; crossing 38.04 s sunflower, 31.84 s soybean). Whitaker:
  `rerun_2026-09-22_rig/rerun_residuals.csv`, `lit_whitaker`, `sensible`,
  `law_log_loading_surface` row, stage-1 fields only (h = 96.89 and 75.71
  W/m2K; crossing 51.57 and 55.24 s). The rig's frozen anchors were checked
  equal to the sibling's, and the measured window start equal. Area per kg
  dry 0.5337 (sunflower) and 0.6377 (soybean) m2/kg on both legs, every
  radius.
* Law: `dxr_common.soybean_law()` (the carried law, not the held one), the law
  of the paper's printed space-resolved values.
* Scoring: `dxr_common.window_score` (committed rule) on the decimated trace
  (every fifth stride, as committed) and on every accepted stride.
  Cross-checked against item 31's scores to the float.

## Environment

PYTHONHASHSEED=1, OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1,
PYTHONDONTWRITEBYTECODE=1, PYTHONIOENCODING=utf-8, PYTHONPATH=<repo>/src,
`.venv\Scripts\python.exe`; one march at a time. Wall: identity 24 marches
about 4 min; declared 16 marches 1.4 to 10.7 s each; extension 32 marches
1.4 to 5.8 s each.

## Order of work

1. Identity (log_identity.txt): 3,732 of 3,732 exact.
2. RESULT.md section 0 written; declared phase run (log_declared.txt). The
   first pass stopped with a `ValueError` in the script's own window
   statistics (an empty list, because the soybean 0.885 mm 120-cell march
   refused before the window); the statistics were guarded and the phase
   re-run; the five marches already written were skipped, not re-marched.
   Nothing in the march changed.
3. Section 0b written; extension run (log_extension.txt).
4. `analyze_item32.py` (log_analyze.txt, analysis.json, runs_table.csv).

## Files

* `runs/<tag>.json` per march (summary, both scores, ledgers, window
  statistics), `runs/<tag>_steps.csv` (every accepted stride, repr floats),
  `runs/<tag>_trace.csv` (the scored decimated trace, crossing included).
  Tags: `<leg>_<condition>_<r0885|rown>_n<cells>` and `ext_<c4|c16>_...`.
* `runs_table.csv`: one row per march, including the state at the last
  accepted stride of every refused march.

## The refusal, as recorded

Every refusal is the certified inner Newton's `line search failed inside the
admissible storage domain` after seven retries (the committed ladder retries at
twice the stride, re-capped by the ceiling). At the committed 5 s ceiling the
last accepted stride of a refused march has the outer cell at 0.003 to 0.012
kg/kg and the surface datum at 0.001 to 0.009 kg/kg, the coefficient varying
by 15 to 61 across the particle; at the 5/16 s ceiling the march refuses two to
seven strides after the crossing, at coefficient ratios of 1 to 7 and outer
loadings of 0.03 to 0.14 kg/kg. The item does not diagnose why the step fails
(the module is read only and the brief is a measurement).
