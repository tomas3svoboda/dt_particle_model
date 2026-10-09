# Measurement: the voidage-Reynolds definition against the coded film closure (2026-09-22)

Read-only measurement. No source edited, no recommendation made. Grounds:
`paper/analysis/datasets/faner2008_thesis_4_5_3/DATASET_RECORD.md` and
`source_reading.md` (Equation 4.23, `Re_epsilon = D_EE rho_g v_g / (mu_g (1-eps))`,
divides by the solid fraction).

## Sites that compute a voidage Reynolds number or evaluate Coletto B7

| Site | Formula as coded | Divides by | eps source |
|---|---|---|---|
| `src/dtdc_simulator/core2/cell_engineering_feasibility.py:2088-2092` (`_evaluate_binary_no_inert`, the live BINARY_NO_INERT gas-film closure) | `reynolds_void = W_g * D / (A_t * mu_g * bed_void)` | `eps` | `area_law.bed_void_fraction` |
| `src/dtdc_simulator/core2/qsc_k2_native_film_authority.py:97-102` (`_accepted_bulk_film_scalars`, QSC K=2 native-film authority) | `reynolds = rho_g * v * D / (mu_g * bed_void)` | `eps` | `cell_inputs.film_area_law.bed_void_fraction` |
| `src/dtdc_simulator/core2/particle/bed_film.py:172-214` (`coletto_b7_b10_pair`, the canonical B.7-B.10 oracle; called by `particle/engineering_foundation.py:352` `source_correlated_film_pair`, which is the paper's frozen anchor) | `reynolds = rho_g * v * D / (mu_g * eps)` | `eps` | `PackedBedFilmState.bed_void_fraction` |
| `src/dtdc_simulator/core/dt_solver.py:64` (retired core, design-note docstring) | `Reeps = rho_V * uV * 2rP / (mu_V * eps_b)`, tagged `[DERIVED]`: no primary source in hand at the time, an independent re-derivation converged on the same divide-by-eps form at "moderate" stated confidence | `eps` | design constant |

All four divide by the voidage `eps`; Equation 4.23 divides by the solid
fraction `(1-eps)`. The Nu-to-h conversion (Eq. 4.21, `h = Nu_eps * lambda *
(1-eps) / (D*eps)`) is coded correctly at every site checked; only the
Reynolds argument fed into the correlation differs.

`_evaluate_positive_carrier` (`cell_engineering_feasibility.py:3170+`) does
**not** compute Re_eps: it reads a frozen declared constant,
`gas_side_heat_coefficient_w_m2_k = 106.5` (RS-2 sealed prior, `schema_compatibility_inert`
migrated legacy value), untouched by this question.

## Quantified effect

Bed voidage is not per-tray state in this model: it is the single frozen
constant `eps = 0.40` at every site above (`dtdc_stack.py:118
REFERENCE_BED_VOID`, `qsc10_k4_a0_genesis_authority.py:65
QSC10_K4_A0_FIXED_BED_VOID_FRACTION`, law id
`FIXED_VOID_FRACTION_0P40_MEMORYLESS_NO_COMPACTION_HISTORY_V1`,
`qsc10_k4_rs2_sealed_prior_graph_supplier.py:575,590`,
`particle/thermal_oracle.py:61`). Because the definitional gap enters only
through `eps`, the ratio is the same multiplicative factor at every tray,
every step, everywhere these sites fire:

| Quantity | Ratio (source convention / as-coded) |
|---|---|
| `Re_eps` | `eps/(1-eps) = 0.6667` |
| `Nu_eps`, `h_Q` (correlation exponent 0.579 only) | `0.6667^0.579 = 0.7908` |
| `h_Q` combined with the printed-prefactor gap (`0.6941/0.6949`) | `0.7898` |

Concrete case, the paper's frozen anchor (`source_correlated_film_pair`,
"Coletto source prior", `eps=0.40`, committed in
`paper/analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21/rerun_faner.json`
`part2_frozen_anchors.coletto_film_pair`):

| | As coded | Under Eq. 4.23 |
|---|---|---|
| `Re_eps` | 45.93 | 30.62 |
| `h_Q` (W m-2 K-1) | 131.34 | 103.74 |

Both values sit below the source's own hexane-family Re_eps range, 189.06 to
275.29 (`tables.csv`, 12 rows), and the full 36-run domain, 18.80 to 275.29.
The corrected Re_eps (30.62) is further below that domain than the coded
value (45.93). Because `eps=0.40` is fixed and below the source's own
implied bed-voidage span (0.531-0.734, all runs), the 0.790 ratio applies
uniformly across the production march's trays; no per-tray Reynolds
extraction from `docs/evidence/d11_n_class_face_2026-09-21/lanes/…op200s4pc4.json`
changes this, since the ratio is `eps`-only and `eps` does not vary by tray
or time step anywhere in the checked code.

## Surfaces touched

`eps=0.40` is pinned at the A0 genesis authority (named packing law), the
RS-2 sealed prior graph supplier, the top-level `dtdc_stack` reference
constant, and the particle module's default — all frozen/pinned QSC-10 K4
surfaces, not a typed `DECLARED_ENGINEERING_ASSUMPTION` row. `COLETTO_B7_NUSSELT_PREFACTOR`
/ `REYNOLDS_EXPONENT` are P1EF-era frozen parameters (`bed_film.py` docstring:
"Gate 1g must not vary... h_Q and h_M... does not qualify the derived
voidage-corrected Reynolds convention"). The RS-2 sealed prior's
`106.5` is a separate frozen legacy constant unaffected by this question.

## Statement

Every site in this codebase that evaluates the Coletto B7 correlation divides
the voidage Reynolds number by `eps`; the digitized source divides by
`(1-eps)`. Bed voidage is a single frozen `0.40` constant everywhere in the
model rather than a per-tray quantity, so the two conventions differ by one
fixed multiplicative factor on the film coefficient, 0.790, at every site and
every tray; the frozen anchor's own numbers (Re_eps 45.93 vs 30.62, h_Q 131.34
vs 103.74 W m-2 K-1) are the concrete instance. Both readings sit below the
correlation's own hexane-family Reynolds range (189-275).
