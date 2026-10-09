# Item 33 — the refinement ladder re-run on the D20 tree (F30, run locally)

`physically_qualifying: false. plant_predictive: false. nothing fitted, nothing marched by this item.`

**Fills** the brief's collection of the locally run F30 campaign against the
ladder passage of `sec_results.tex` (Sec. 6.1) and the supplement's
`sup:stressconsolidation`, which report the pre-D20 ladder of item_01
(`../item_01_refinement_ladder/RESULT.md`, rerun 2026-09-20: 36 jobs, 5
complete, two adjacent mesh pairs, no three-level chain). **Outcome: the ladder
still forms no three-level chain on any axis, scenario or observable set; four
trajectories complete, one adjacent pair forms; every refusal is the D21
departure-tangent class, reached through the failed-leaf floor at N ≤ 48 and
directly, in the mesh bootstrap, at N = 96.**

## 1. What ran

The 36 jobs of `tools/farm/campaigns/F30_refinement_ladder_d20.json`
(N ∈ {12, 24, 48, 96} × Δt ∈ {0.075, 0.0375, 0.01875} s × {composition only,
temperature only, combined}, one job per process, one mesh, one step and one
scenario per job) through `scripts/qualify_phase1_particle_engineering.py` at
the frozen engineering case (`nominal_log_center`, film multiplier 1.0,
k = 0.24 W/m/K, R = 0.885 mm), on the workstation, from HEAD `b9a48f03`
(the launcher's commit; the jobs record source provenance, combined sha256
`c8ab779c…` on all 27 output files), eight at a time with the pinned
environment, 2026-09-26 21:37 to 2026-09-27 11:09. Horizon of every job: the
stress schedule's three square segments, 0.225 s. Evidence placed in
`docs/evidence/f30_local_2026-09-27/` (logs, ledger, launcher,
`COLLECTION.json`, `MANIFEST.sha256`); scored by `analyse_f30.py` beside this
file, which imports no solver.

**Exit codes.** 0 = `PASS_ENGINEERING_NUMERICAL_CASE`; 2 =
`FAIL_CLOSED_ENGINEERING_NUMERICAL_CASE`; 1 = an uncaught exception, no output
file. **Tally: 3 × 0, 24 × 2, 9 × 1.** The exit code is *not* the trajectory
verdict, in both directions:

* `f30n48d01875comb` exits **0 with a refused trajectory**: its refusal matches
  a registered R7 wall of the harness (`R7_WALL_SIGNATURES`, combined, N = 48,
  Δt = 0.01875 s, "subdivision exhausted its documented minimum dt"), which the
  runner accounts as a verified fail-closed exclusion, not an evidence child.
* `f30n12d075comp` and `f30n24d075comp` exit **2 with complete trajectories**:
  every per-job gate is true, the release contract passes, determinism is bit
  exact, and the campaign's structural pass fails on one term alone — the
  composition-only pulse endpoint (hot, rich) is flagged
  `particle_dry_storage_chart_compatible: false`, a property of the scenario's
  boundary, not of the solve. No composition-only job can exit 0 on this tree.
* `f30n48d01875comp` is also a registered R7 wall and exits 2 for the same
  composition-only reason.

So **four trajectories complete** (the runner's own `completed`): composition
only N = 12 and N = 24 at Δt = 0.075 s, temperature only N = 12 at 0.075 s,
combined N = 48 at 0.0375 s.

## 2. The 36 jobs

"Last accepted t" is the last committed time on the requested grid (for a
refused same-cell job, the input time of the rolled-back macrostep); "furthest
uncommitted leaf" is how far the provisional failed-leaf subdivision got before
exact macro rollback (never committed). "D21 tangent text inside" is the
departure-tangent refusal found in the job's own attempt record (for N = 96, in
its traceback), with the scaled residual it printed against the 2.000e-11
contract. Ledger: worst normalized cumulative ledger residual over accepted
steps (complete jobs) or over the uncommitted provisional leaves (refused
jobs). Wall: the launcher's, which for six N = 96 jobs includes a host sleep
(Sec. 5).

| job | N | dt (s) | scenario | exit | runner status | trajectory | refusal (typed) | last accepted t (s) of 0.225 | furthest uncommitted leaf (s) | D21 tangent text inside | worst ledger (cum.) | wall s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| f30n12d075comb | 12 | 0.075 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.17285 | no admissible root 1.061e-10 | 1.47e-13 | 360.1 |
| f30n12d075comp | 12 | 0.075 | comp | 2 | FAIL_CLOSED | complete | - | 0.225 | - | - | 2.89e-13 | 129.8 |
| f30n12d075temp | 12 | 0.075 | temp | 0 | PASS | complete | - | 0.225 | - | - | 1.59e-13 | 151.6 |
| f30n12d0375comb | 12 | 0.0375 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.17402 | no admissible root 0.007844 | 6.15e-14 | 431.1 |
| f30n12d0375comp | 12 | 0.0375 | comp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.18047 | no admissible root 0.007805 | 6.11e-14 | 441.7 |
| f30n12d0375temp | 12 | 0.0375 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.18076 | reverse-core seed | 7.24e-14 | 850.8 |
| f30n12d01875comb | 12 | 0.01875 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.1875 | 0.18779 | no admissible root 0.05041 | 4.24e-14 | 434.9 |
| f30n12d01875comp | 12 | 0.01875 | comp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.1875 | 0.19746 | reverse-core seed | 4.00e-14 | 1272.5 |
| f30n12d01875temp | 12 | 0.01875 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.1875 | 0.19746 | reverse-core seed | 1.05e-13 | 516.7 |
| f30n24d075comb | 24 | 0.075 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16465 | no admissible root 9.998e-11 | 1.08e-13 | 555.0 |
| f30n24d075comp | 24 | 0.075 | comp | 2 | FAIL_CLOSED | complete | - | 0.225 | - | - | 1.02e-13 | 94.0 |
| f30n24d075temp | 24 | 0.075 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.17285 | reverse-core seed | 2.14e-14 | 561.2 |
| f30n24d0375comb | 24 | 0.0375 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.15908 | no admissible root 0.01271 | 7.24e-14 | 696.8 |
| f30n24d0375comp | 24 | 0.0375 | comp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.15938 | no admissible root 0.01281 | 6.51e-14 | 534.6 |
| f30n24d0375temp | 24 | 0.0375 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16055 | no admissible root 0.01272 | 7.31e-14 | 621.3 |
| f30n24d01875comb | 24 | 0.01875 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16787 | no admissible root 0.01777 | 1.34e-14 | 663.9 |
| f30n24d01875comp | 24 | 0.01875 | comp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.16875 | 0.17227 | no admissible root 0.04683 | 1.43e-14 | 929.7 |
| f30n24d01875temp | 24 | 0.01875 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.16875 | 0.17256 | no admissible root 0.04639 | 4.73e-14 | 1105.6 |
| f30n48d075comb | 48 | 0.075 | comb | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16465 | reverse-core seed | 1.60e-14 | 596.0 |
| f30n48d075comp | 48 | 0.075 | comp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16172 | reverse-core seed | 1.02e-14 | 456.7 |
| f30n48d075temp | 48 | 0.075 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16553 | no admissible root 0.01024 | 2.25e-14 | 897.3 |
| f30n48d0375comb | 48 | 0.0375 | comb | 0 | PASS | complete | - | 0.225 | - | - | 1.90e-14 | 310.3 |
| f30n48d0375comp | 48 | 0.0375 | comp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.15791 | reverse-core seed | 5.32e-14 | 812.4 |
| f30n48d0375temp | 48 | 0.0375 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.16348 | reverse-core seed | 2.09e-14 | 919.6 |
| f30n48d01875comb | 48 | 0.01875 | comb | 0 | PASS | refused | same-cell floor (R7 wall) | 0.15 | 0.16025 | reverse-core seed | 2.94e-14 | 1361.4 |
| f30n48d01875comp | 48 | 0.01875 | comp | 2 | FAIL_CLOSED | refused | same-cell floor (R7 wall) | 0.15 | 0.16201 | reverse-core seed | 6.62e-14 | 1306.8 |
| f30n48d01875temp | 48 | 0.01875 | temp | 2 | FAIL_CLOSED | refused | same-cell floor | 0.15 | 0.15996 | reverse-core seed | 1.19e-14 | 1470.3 |
| f30n96d075comb | 96 | 0.075 | comb | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 3501.8 |
| f30n96d075comp | 96 | 0.075 | comp | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 3526.8 |
| f30n96d075temp | 96 | 0.075 | temp | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 3552.2 |
| f30n96d0375comb | 96 | 0.0375 | comb | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 44006.0 |
| f30n96d0375comp | 96 | 0.0375 | comp | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 44184.5 |
| f30n96d0375temp | 96 | 0.0375 | temp | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 44175.7 |
| f30n96d01875comb | 96 | 0.01875 | comb | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 43473.8 |
| f30n96d01875comp | 96 | 0.01875 | comp | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 43967.0 |
| f30n96d01875temp | 96 | 0.01875 | temp | 1 | none (no output) | refused | CutEventMacrostepError, N=96 bootstrap | 0 | - | no admissible root 5.279e-11 | - | 43919.6 |

On every complete job the maximum scaled residual is 1.13e-11 to 9.34e-11
(criterion 1e-10), the worst step and cumulative ledger residuals 1.9e-14 to
2.9e-13, restart and serial replay bit exact, and each commits exactly one
exact interior-face event (atomic face-event recovery) and no accepted
subdivision. On all 23 refused same-cell jobs the failure class is
`same_cell_nonlinear_minimum_dt_exhausted`, no failed state is committed,
every rejected attempt preserves its rollback identity, and the uncommitted
leaves' ledgers stay at or below 1.5e-13. **Every one of those 23 refusals
sits at or after t = 0.15 s**, the start of the return segment (the exact
return to the hot, lean baseline after the pulse); the nine N = 96 refusals
come before any accepted time.

## 3. The chains

The acceptance criterion (Sec. 3.4 of the main text, S3.4) has two halves:
adjacent differences must contract (ratio ≤ 0.75, which needs three completed
levels on one axis for one scenario), and the final two levels must agree to
1 % in normalized front position, component inventories, integrated component
fluxes and event time, and to 0.1 K in local temperature (the harness's 22
integrated convergence observables, scored as in item_01).

**Of the 30 possible three-level chains (per scenario: six mesh chains 12-24-48
and 24-48-96 at each Δt, four time chains at each N), none forms; one has two
completed levels, eight have one, 21 have none. No observed order exists
anywhere.**
Where each chain breaks (full list `chains.csv`):

| chain | completed levels | breaks at | refusal there |
|---|---|---|---|
| composition only, mesh 12-24-48, Δt 0.075 | **2** | N = 48 | floor, reverse-core tangent seed |
| composition only, mesh 24-48-96, Δt 0.075 | 1 | N = 48; N = 96 | floor (reverse-core); N = 96 bootstrap tangent |
| composition only, time at N = 12 | 1 | 0.0375; 0.01875 | floor (no admissible root 7.8e-3; reverse-core) |
| composition only, time at N = 24 | 1 | 0.0375; 0.01875 | floor (no admissible root 1.28e-2; 4.68e-2) |
| temperature only, mesh 12-24-48, Δt 0.075 | 1 | N = 24; N = 48 | floor (reverse-core; no admissible root 1.02e-2) |
| temperature only, time at N = 12 | 1 | 0.0375; 0.01875 | floor (reverse-core; reverse-core) |
| combined, mesh 12-24-48, Δt 0.0375 | 1 | N = 12; N = 24 | floor (no admissible root 7.8e-3; 1.27e-2) |
| combined, mesh 24-48-96, Δt 0.0375 | 1 | N = 24; N = 96 | floor (no admissible root 1.27e-2); N = 96 bootstrap |
| combined, time at N = 48 | 1 | 0.075; 0.01875 | floor (reverse-core); floor at the registered R7 wall (reverse-core) |
| the other 21 chains | 0 | every level | 11 on the floor alone; 10 through N = 96, also on the bootstrap tangent |

**The one adjacent pair: composition only, N = 12 → 24 at Δt = 0.075 s.**
15 of the 22 observables meet their criterion and 7 do not
(`pair_differences.csv`):

| observable | criterion | difference | met |
|---|---|---|---|
| pulse water flux time integral | 1 % | **8.043 %** | no |
| return water flux time integral | 1 % | 2.415 % | no |
| baseline water flux time integral | 1 % | 1.905 % | no |
| baseline hexane flux time integral | 1 % | 1.643 % | no |
| cumulative boundary hexane out | 1 % | 1.013 % | no |
| max instantaneous radial temperature spread | 0.1 K | **0.4093 K** | no |
| peak particle temperature | 0.1 K | 0.3990 K | no |
| pulse hexane flux time integral | 1 % | 0.864 % | yes |
| cumulative boundary energy out | 1 % | 0.775 % | yes |
| return hexane flux time integral | 1 % | 0.347 % | yes |
| final front z (wet volume fraction) | 1 % | 0.268 % | yes |
| final total hexane | 1 % | 0.0895 % | yes |
| final, minimum front position over radius | 1 % | **0.0894 %** | yes |
| final total energy | 1 % | 0.0198 % | yes |
| final pore vapour water | 1 % | 0.0122 % | yes |
| final retained and total water | 1 % | 0.0078 % | yes |
| final interface temperature | 0.1 K | 0.0019 K | yes |
| maximum front position, minimum material temperature, free liquid water | — | 0 | yes |

The pattern is the pre-D20 12 → 24 pair's (combined scenario, 14 of 22, worst
8.161 % on the same pulse water flux integral and 0.4529 K on the same radial
spread): the front and the integral inventories inside their criteria by an
order of magnitude, the finite-window flux integrals and the two peak
temperatures outside. It is a different scenario, so this is a like pattern,
not a like-for-like number.

## 4. Against the pre-D20 ladder, cell by cell

The pre-D20 record is item_01's rerun of 2026-09-20 (its `cells/`, all three
scenarios per process). Its source provenance differs from this run's in
exactly **two of 65 tracked files**,
`src/dtdc_simulator/core2/particle/bed_film.py` and
`src/dtdc_simulator/core2/particle/engineering_foundation.py`, both last
changed in `4cbedea9` (the commit that carries D20 level 0); the runtime
(CPython 3.14.7, numpy 2.4.6, scipy 1.18.0) and the coefficient case are
identical. Full table `comparison_pre_d20.csv`.

| verdict | cells | which |
|---|---|---|
| refused pre-D20, **completes now** | 3 | temperature only N = 12 @ 0.075; composition only N = 24 @ 0.075; combined N = 48 @ 0.0375 |
| completed pre-D20, **refuses now** | 4 | combined N = 12 and N = 24 @ 0.075 (floor, no admissible root at 1.06e-10 and 1.00e-10); composition only N = 48 @ 0.0375 (floor, reverse-core); composition only N = 96 @ 0.075 (bootstrap tangent) |
| completes on both | 1 | composition only N = 12 @ 0.075 — **not bit-identical**: front position moves 0.045 %, total hexane 0.63 %, the finite-window flux integrals 6.2 to 13.2 %, cumulative boundary energy out 30.9 %, radial temperature spread 11.035 → 10.292 K |
| refuses on both | 28 | — |

**The refusal class moved at N = 96 only.** At N ≤ 48 every refused job, on
both trees, is the same-cell floor with one departure-tangent refusal in its
attempt record; the tangent's text flips between its two forms in six cells
(no admissible root ↔ reverse-core seed), and the one caloric-bracket seed of
the pre-D20 ladder (temperature only, N = 96) is gone. At N = 96 the pre-D20
row had composition only @ 0.075 complete, combined and temperature only @
0.075 on the floor at t = 0.15 s, and the six refined-step jobs refusing at the
shared baseline prefix (`RuntimeError`: time continuation requires a cut
candidate); now all nine refuse **earlier**, inside the mesh bootstrap's
48 → 96 master-face route (the 2026-09-20 route that let N = 96 run at all),
where the exact face arrival's departure tangent finds no admissible
primary-drainage root at q_h = 0.3782, scaled residual **5.279e-11** against
2.000e-11 — bit-identical in all nine jobs, because that bootstrap step is
built at the segment duration and is independent of the requested Δt. **No
N = 96 job reaches any accepted time; none reaches further than its pre-D20
cell.** The pair structure therefore changes from two pairs (combined 12 → 24,
composition only 12 → 96) to one (composition only 12 → 24).

## 5. The N = 96 walls

The three Δt = 0.075 s cells took 3,502 to 3,552 s; the six refined-step cells
43,474 to 44,185 s. The six ran through a host sleep: the Windows System log
records sleep at 2026-09-26 23:09:50 and resume at 2026-09-27 10:25:52 (the
clock stepped by 40,560 s), and the launcher's wall includes it. Awake time is
therefore about 2,900 to 3,600 s, the same as the 0.075 s trio, which finished
before the sleep. **The 12-hour walls are not solver time spent before the
refusal**: their stdout is empty, the refusal is raised in
`_bootstrap_first_steps` before any campaign step (the traceback), and it is
bit-identical to the 0.075 s cells'. The finer steps bought nothing: the code
path up to the refusal does not depend on Δt.

## 6. What this ladder proves and does not

* The complete moving-front trajectory **does not form a three-level chain on
  this tree**, in mesh or in time, for any scenario; so the contraction half of
  the criterion is still untested, not failed, and no observed order exists.
* The coarse levels stop on the same-cell failed-leaf floor
  (2.930e-4 s, exact macro rollback); inside each of those 23 jobs the attempt
  after the initial same-cell refusal is the exact-face dispatch
  (`accepted_history_uncertified_face_dispatch`), whose departure tangent refuses
  — 12 with no admissible primary-drainage root (two near the literal,
  1.00e-10 and 1.06e-10; ten at 7.8e-3 to 5.0e-2) and 11 with a right-limit
  seed requesting reverse-core motion. The fine level (N = 96) stops on the
  exact-face tangent itself. **This is the D21 class, not a different
  obstruction**: the same refusal texts, the same bimodal residuals (a
  near-literal group and a 1e-2 plateau five decades above), the same
  reverse-core group that prints no residual.
* What is new is **where** the N = 96 tangent refuses (the mesh bootstrap,
  not the campaign) and **which group** it falls in: 5.279e-11 is the
  near-literal group, 0.027 of D21's tightest bound (1.938e-9), not the
  outward-core plateau. But it sits at q_h = 0.378, outside the 0.217 to 0.241
  band D21 built that bound on, so whether it is within the arrival's own
  resolution at its own state is **not** established here.
* Even if every tangent refusal under D21's bound were admitted (the two
  near-literal N = 12/24 combined jobs and the N = 96 bootstrap), no chain with
  a level below N = 96 would form: each still needs a reverse-core job or a
  plateau job (every chain through N = 48 at Δt = 0.075 s does). The three time
  chains lying wholly at N = 96 are undetermined, since no N = 96 job ran past
  the bootstrap (item_01 also recorded single-mesh N = 96 invocations dying in
  the campaign's determinism audit on the pre-D20 tree). Admission of a tangent
  would in any case not by itself complete a job.
* It does **not** show that the D20 correction caused the changed cells. The
  two trees differ in two particle source files; which change moves which cell
  is not measured.

## 7. Sentences the supplement could carry

Beside the existing refinement passage of `sup:stressconsolidation` (nothing in
the main text is edited by this item):

> A re-run of the full ladder at the present state of the formulation
> (36 jobs, item 33) forms no three-level chain either. Four trajectories
> complete — composition only at N = 12 and 24 and temperature only at N = 12,
> all at Δt = 0.075 s, and the correlated scenario at N = 48, Δt = 0.0375 s —
> admitting one adjacent pair, composition only from 12 to 24 cells, on which
> 15 of the 22 observables meet their criterion, the worst again the pulse
> water flux integral (8.04 %) and the radial temperature spread (0.409 K),
> with the front position at 0.089 %. Every refusal is the departure-tangent
> refusal described above: at N ≤ 48 reached through the failed-leaf floor, and
> at N = 96 inside the mesh continuation itself, where the tangent misses the
> contract at 5.28 × 10⁻¹¹ and that level reaches no accepted time.

A one-sentence version, if the page is tight:

> Re-run at the present state of the formulation, the ladder again forms no
> three-level chain (four of 36 jobs complete, one adjacent pair), and every
> refusal is the same departure-tangent refusal, reached through the
> failed-leaf floor at N ≤ 48 and directly in the mesh continuation at N = 96.

## 8. Claim and non-claim

**Claim:** 36 jobs ran from `b9a48f03` on the workstation; exit codes 3 × 0,
24 × 2, 9 × 1; four trajectories complete, every contract gate the structural
pass reads true on each; one adjacent pair, 15 of 22 criteria met; no three-level chain in
mesh or time; every refusal carries the D21 departure-tangent text; against
the pre-D20 ladder 3 cells gained, 4 lost, 1 changed without changing verdict,
28 unchanged in verdict, the refusal class moved only at N = 96.

**Non-claim:** not a convergence demonstration, not an observed order, not a
contraction ratio; no tolerance, budget, seed, floor or refusal moved, no
criterion re-scored; the R7 wall PASS is the harness's registered exclusion,
not a trajectory; not a statement that D20 is right or wrong or that it caused
any cell's change; not physically qualifying, not plant predictive, no
calibration, no QSC-10 content.

Files: `analyse_f30.py` (the analysis), `jobs_table.csv`, `chains.csv`,
`pair_differences.csv`, `comparison_pre_d20.csv`, `tangent_refusals.csv`,
`analysis.json`, `collection_rows.json`, `DETAIL.md`, `MANIFEST.sha256`.
