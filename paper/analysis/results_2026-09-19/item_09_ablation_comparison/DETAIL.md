# Item 9 — non-isothermal against the isothermal-ablation oracle

Fills **declared missing result 9** (`sec_missing.tex` item 9;
`sec_validation.tex` §The Faner experiments) as far as the ladder allows: the
comparison over a *complete* trajectory still needs the admissible ladder of
item 1, which is refused at every level above N = 48 and at most cells below
it, but the comparison over every accepted step that exists at this HEAD is
performed here and is a definite number.

## Method

`run_item09_ablation.py`. The oracle is eq:faneroracle,
s = R((X_h − X_e)/(X_c − X_e))^(1/3), evaluated through the committed
authority `core2.sorption_interface.front_radius_m` on the live X_c(T),
X_e(T, P) pair at the pressure-feasible activity ceiling. For every accepted
step of every completed job in `item_01_refinement_ladder/cells/`, the total
hexane loading is formed from the step's own hexane inventory and the
particle's dry mass, X_c and X_e are evaluated at the step's own interface
temperature, and the oracle front is compared with the solver's. The harness
stores the front as the wet volume fraction z = (s/R)³; both coordinates are
reported. `pytest` over the four oracle suites (`test_core2_thermal_oracle`,
`test_core2_gate1g_limit_oracles`, `test_core2_falling_rate_closed_form_oracle`,
`test_core2_front`) and `scripts/qualify_gate1g_limit_oracles` were also run.

## The comparison

Fifteen accepted steps across four completed jobs (composition-only at N = 12,
Δt = 0.075 and at N = 48, Δt = 0.0375; combined at N = 12 and N = 24,
Δt = 0.075); **twenty-three jobs refused** and are listed in
`item09_ablation.json`. Every scored step is in the RECEDING_FRONT regime.

| quantity | range over the fifteen steps |
|---|---|
| relative difference in s/R | **9.43 % to 10.93 %** |
| relative difference in the wet volume fraction z | **25.71 % to 29.33 %** |

The sign is the same on every step: the resolved non-isothermal front sits
**above** the isothermal-ablation oracle, that is the resolved core retains
more wet volume than the geometric ablation identity implies at the same
inventory. The deviation grows monotonically along each job (9.75 → 9.91 →
10.57 % on the N = 12 combined job) and is mesh-insensitive at the one
available pair (9.75/9.91/10.57 % at N = 12 against 9.71/9.89/10.54 % at
N = 24 — a mesh effect of 0.03 points against a trajectory effect of 0.8).
The one job that ran at a finer time step, composition-only at N = 48 and
Δt = 0.0375 s, resolves the same drift over twice as many steps and reaches the
same place: 9.43, 9.67, 9.93, 10.21, 10.53, 10.86 % over its six steps, against
9.75, 9.91, 10.57 % at N = 12, Δt = 0.075 s over the same 0.225 s. The drift is
therefore a property of the trajectory, not of the discretization.
Per-step table: `item09_step_comparison.csv`.

## The single-state oracles, re-run

`scripts/qualify_gate1g_limit_oracles` at HEAD: `numerical_oracles_passed:
true`, `physically_qualifying: false`. The ideal-binary chain-rule recovery
closes to `1.27e-16` relative and the complete coupled API flux to `2.47e-15`,
while the finite-pressure relative coefficient offset at the smallest water
level is **0.037992** — the open discrepancy `sec_verification.tex` records as
3.819 %; at this HEAD it reads **3.799 %**, and the paper's figure should be
re-quoted. The recovered-pressure limit likewise passes with all component
energy ledgers closing and all algebraic total fluxes recovering zero.

The four oracle suites: **75 passed, 0 failed** (26.1 s).

## What is still missing

A complete trajectory. The horizon available here is the harness's three
0.075 s macro steps, 0.225 s in all, over which the front moves from z = 0.216
to about z = 0.193 — roughly a tenth of the wet volume. A comparison over a
full recession needs the admissible ladder of item 1, and the ladder is refused
at 23 of the 27 jobs that started, with
`AdaptiveSameCellMacrostepError: same-cell failed-leaf subdivision exhausted
its documented minimum dt 0.000292969 s with exact macro rollback`, and at
N = 96 with `candidate crosses a master face; use the exact event chart, not
mesh continuation`.

## Claim boundary

Fifteen steps over at most 0.225 s of a stress schedule, not a drying trajectory. The
loading used in the oracle is formed from the step's hexane inventory over the
particle's dry mass (4/3)πR³ρ_dm,p, which is the convention the dry-shell
storage uses; it is stated because it is a choice. The oracle carries no
experimental authority — it is the authors' own front relation used as an
oracle, as `sec_validation.tex` states. Not physically qualifying, not plant
predictive.
