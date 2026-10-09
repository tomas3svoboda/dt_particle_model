# Item 22, 2026-09-25: which of the Faner (2019) unprinted conditions actually moves the comparison

The F32 farm campaign swept the four experimental quantities the 2019 journal
article does not print per run, and scored the model's constant-rate leg against
the source's own digitized drying curves over that box. This record collects the
return (48 jobs, all `rc=0`), assembles the measurement, and states what it does
and does not establish.

`fitted: false. physically_qualifying: false. plant_predictive: false.`
`production_wired: false. qsc10_complete: false.`

**The short answer.** Of the four unprinted quantities, **two do not matter and
two do.** The sample mass leaves the score bit-identical, so it cannot be the
source of the disagreement. The initial loading is pinned to within digitization
error by the trace's own first point, so the article's printed 0.45 to 0.50
kg/kg range is not a real degree of freedom either. What remains is the **bed
depth** and the **gas velocity**, and they are exactly the two the article prints
least well: a range of three to five particle diameters for the first, nothing at
all for the second. A single condition set does **not** serve both traces on the
signed reading of the acceptance band, and serves both on the magnitude reading
only at two coordinates, both of which demand a gas velocity above the one the
companion thesis rig itself ran at.

---

## 1. What ran

| | |
|---|---|
| campaign | `F32`, `tools/farm/campaigns/F32_faner_sensitivity_sweep.json` |
| runner | `tools/farm/faner_sensitivity_sweep.py` |
| commit | `bad3513ef877910642e19df8e36e603f940c8710` (all 48 jobs, one commit) |
| host | farm 2, TIASERVER, AMD EPYC 9474F, 32 logical processors, 128 GB RAM |
| interpreter | CPython 3.14.5, numpy 2.4.6, scipy 1.18.0 |
| jobs | 48 submitted, **48 done, exit code 0 on every one**, 0 failed, 0 outstanding |
| wall | 22.144 s to 99.065 s per job, 2923.533 s summed; all 48 started within 89 s and ended within 12 s of each other, so the summed wall is queue overlap on a 32-processor host, not per-job compute |
| ran | 2026-09-24T20:54:26Z to 2026-09-24T20:56:17Z |
| grid | 20 points per job, **960 in total, 960 scored, 0 typed refusals** |
| return | `results_TIASERVER_20260925T023602Z.zip`, sha256 `b28eddf4fee4b86f84ec0a1a03cd909141b77e353075a394081ddac6e17d2b4e` |

The 48 jobs are the outer product of the **series** (the two journal traces:
soybean at 120 °C, sunflower at 100 °C) with the outer **condition** axes
(sample mass 0.015, 0.030, 0.045, 0.060 kg; gas velocity 0.05, 0.10, 0.15, 0.20,
0.25, 0.30 m/s): 2 × 4 × 6 = 48. Each job ran the two inner axes itself —
initial loading 0.40, 0.45, 0.50, 0.55 kg/kg and bed depth 2, 3, 4, 5, 6 particle
diameters — for 20 grid points each. Every job's receipt, argv, outcome and
artifact digest is in `docs/evidence/farm_f32_2026-09-25/COLLECTION.json`; the
per-job outcome table is `f32_job_summary.csv`.

## 2. The score, and the band it is read against

The score is the runner's own, not one invented here
(`tools/farm/faner_sensitivity_sweep.py`, its module docstring and
`score_grid_point`):

> The score, "window residual", is
> `(predicted_constant_rate_1_s - measured_constant_rate_1_s) / measured_constant_rate_1_s`
> where the predicted rate is `h * (T_gas - T_boil) * area_per_kg_dry / dh_vap`
> (the same constant-rate formula `run_rerun_faner.py` uses) and the measured
> rate is minus the least-squares slope of the windowed (time, loading) points.

`h` is the Whitaker (1972) coefficient — the item-5 record's own best-grounded
independent comparator — evaluated at the swept velocity and otherwise at the
material's printed hexane-vapour family state. The window runs from the frozen
identification instrument's critical loading plus its margin, 0.1888 + 0.05 =
0.2388 kg/kg, up to the swept initial loading. A window with fewer than two
points is a typed refusal; none fired.

**The acceptance band is the item-5 record's own thesis-run outcome at exactly
this closure**, read off `item_05_faner_comparison/rerun_2026-09-22_rig/`'s
`rerun_film_band.csv` rather than retyped. `RESULT_2026-09-22.md` states it:

> Marched at Whitaker, the constant-rate demand is predicted to **0.98, 0.88,
> 0.91 and 0.82** of the measurement and the constant-rate leg ends at **0.89 to
> 1.12** of the measured transition time. Nothing is fitted to any drying curve.

The first two of those four are the Faner (2008) thesis runs, the conditions for
which the source does print a per-run table. In the F32 score's units:

| reading | definition | value |
|---|---|---|
| thesis-run band, **signed** | `window_residual` between the two thesis runs' own values | **−0.117872 ≤ r ≤ −0.018422** |
| thesis-run band, **magnitude** | the same, read as a size only | **\|r\| ≤ 0.117872** |
| 2008 sunflower 136 °C | the better thesis run | −0.018422 (0.9816 of the measurement) |
| 2008 soybean 136 °C | the worse thesis run | −0.117872 (0.8821) |
| 2019 soybean 120 °C, item-5 march | the journal trace as item 5 left it | −0.093890 (0.9061) |
| 2019 sunflower 100 °C, item-5 march | the journal trace as item 5 left it | −0.176679 (0.8233) |

The signed band is one-sided and negative: on **both** thesis runs, where the
conditions are printed, this closure under-predicts the measured constant rate.
That sign is part of the band, and section 5 is where it bites.

## 3. Two of the four swept quantities cannot carry the disagreement

**Sample mass is exactly degenerate.** Across the four masses, at every one of
the 12 (material, velocity) pairs, the twenty window residuals are
**bit-identical**. This is not an approximation and not a small effect: the
transfer area per kg dry meal of a layer of a given depth is
`1 / (rho_env (1 - eps) depth)`, in which the mass cancels, so the score cannot
see it. The 48 jobs therefore carry **240 distinct scored conditions**, not 960.
The runner says so in its own docstring and the return confirms it.

**Initial loading is pinned by the trace itself.** The sweep runs 0.40 to 0.55
kg/kg, but the source's own digitized curve starts at 0.490628 (soybean) and
0.499279 (sunflower) kg/kg, so 0.50 and 0.55 select the *same* window and 0.40
discards real measured points. Only three distinct windows exist per trace:

| trace | X₀ = 0.40 | X₀ = 0.45 | X₀ = 0.50 and 0.55 |
|---|---:|---:|---:|
| soybean 120 °C | 2 points, 6.0016e−3 s⁻¹ | 3 points, 7.3936e−3 | 4 points, 7.2643e−3 |
| sunflower 100 °C | 4 points, 4.7411e−3 | 5 points, 4.8504e−3 | 6 points, 4.9904e−3 |

Both digitized start values sit inside the article's own printed initial-loading
range of 0.45 to 0.50 kg/kg, which is therefore corroborated, not open.

**What the axes are worth.** Residual span along one axis with the other two
held fixed, over all such lines (median of the spans):

| axis | swept range | soybean | sunflower |
|---|---|---:|---:|
| bed depth | 2 → 6 particle diameters | **1.116** | **0.968** |
| gas velocity | 0.05 → 0.30 m/s | **0.688** | **0.615** |
| initial loading | 0.40 → 0.55 kg/kg | 0.170 | 0.034 |
| sample mass | 0.015 → 0.060 kg | 0.000 (exactly) | 0.000 (exactly) |

The signed band is 0.099450 wide. In band widths, the median spans above are
11.2 and 9.7 for the bed depth, 6.9 and 6.2 for the gas velocity, 1.7 and 0.3
for the initial loading, and exactly 0 for the sample mass.

## 4. Where each trace meets the band

The score varies smoothly and monotonically through zero along both live axes, so
that *some* grid point lands in the band is a consequence of continuity and is
not the measurement. The measurement is **where** the crossing sits relative to
what the source prints, and whether the two traces put it in the same place.

| trace | 120 points | residual span | in signed band | within magnitude cap | best condition (v m/s, X₀ kg/kg, depth) | its score |
|---|---:|---|---:|---:|---|---:|
| soybean 120 °C | all | −0.7551 to +1.5716 | 9 | 20 | 0.20, 0.40, 4 | **+0.01317** |
| soybean 120 °C | printed box only | — | 4 of 36 | 8 of 36 | 0.30, 0.45, 4 | **+0.04370** |
| sunflower 100 °C | all | −0.7840 to +0.9400 | 8 | 14 | 0.10, 0.45, 2 | **−0.00353** |
| sunflower 100 °C | printed box only | — | 3 of 36 | 5 of 36 | 0.20, 0.45, 3 | **−0.00401** |

"Printed box" is the sub-box the article's own text admits: bed depth 3 to 5
particle diameters and initial loading 0.45 to 0.50 kg/kg (`conditions.csv`,
Sec. 2.2). Gas velocity has no printed value to restrict.

The soybean best-overall point sits on the two-point window (X₀ = 0.40), the
weakest measured slope in the set and one point above the runner's own refusal
threshold; the printed-box row, which excludes it, is the one to quote.

**Both traces are brought well inside the band by some swept condition** — to
1.3 % and 0.4 % of the measured constant rate, against 1.8 % and 11.8 % on the
two thesis runs where the conditions are printed. Neither number is a fit and
neither is evidence of agreement at the source's real conditions; both say the
under-determined box is wide enough to contain the measurement, which is the
problem, not the result.

## 5. Does one condition set serve both traces? No on the signed band, barely on the magnitude

| test | soybean coordinates | sunflower coordinates | common (v, X₀, depth) | common (v, depth) |
|---|---:|---:|---|---|
| signed band, whole sweep | 9 | 8 | **none** | **none** |
| magnitude cap, whole sweep | 20 | 14 | 9 | (0.10, 2), (0.20, 3), (0.30, 4) |
| signed band, printed box | 4 | 3 | **none** | **none** |
| magnitude cap, printed box | 8 | 5 | 4 | **(0.20, 3), (0.30, 4)** |

The four jointly admitted coordinates inside the printed box, with each trace's
score:

| gas velocity (m/s) | initial loading (kg/kg) | bed depth (diameters) | soybean | sunflower |
|---:|---:|---:|---:|---:|
| 0.20 | 0.45 | 3 | **+0.0966** | **−0.0040** |
| 0.20 | 0.50 | 3 | +0.1161 | −0.0320 |
| 0.30 | 0.45 | 4 | **+0.0437** | **−0.0519** |
| 0.30 | 0.50 | 4 | +0.0623 | −0.0785 |

**The signs differ at every one of them.** Wherever the sweep brings both
journal traces within the size of the thesis-run disagreement, it over-predicts
one and under-predicts the other, whereas on both thesis runs it under-predicts.
No swept condition reproduces the thesis-run signature on both journal traces at
once. That is the operative form of the under-determination statement: the two
journal traces do not merely want unprinted conditions, they want *different*
ones.

## 6. Three independent cross-checks the return passes

1. **The item-5 march lands where it should on the velocity axis.** The item-5
   Whitaker coefficient at the rig's state is 78.649 (soybean) and 80.013
   (sunflower) W m⁻² K⁻¹; the sweep brackets each between its v = 0.10 and
   v = 0.15 values, and linear interpolation in the sweep's own tabulated
   coefficient puts them at **0.1243 and 0.1176 m/s** — both inside the thesis's
   printed approach-velocity range of **0.117 to 0.129 m/s**, which the sweep
   never saw. The sweep's velocity axis contains the companion rig's own
   velocity, and reproduces its coefficient there.
2. **The implied footprint recovers the printed sample holder.** The holder's
   printed cross-section is 0.011 m². The sweep's implied footprint depends only
   on mass over depth; 30 g — the article's own approximate sample mass — at
   **4 particle diameters**, the middle of the printed 3-to-5 range, gives
   0.011417 m² (soybean, 1.038 × printed) and 0.010612 m² (sunflower, 0.965 ×).
   The quantity that cannot move the score is the one that closes the geometry.
3. **The correlation is never read outside its printed bound.** Voidage Reynolds
   number over all 960 points: 93.0 to 593.2, inside Whitaker's printed 10 to
   10000 at every point (`outside_printed_whitaker_reynolds_bound` false, 0 of
   960); the thesis's own hexane-vapour runs sit at 189 to 275, which the sweep's
   v = 0.10 to 0.15 band reproduces. The saturation temperature the runner used,
   68.71 °C at 101 325 Pa from the certified equation of state, matches the
   article's printed boiling temperature of 68.7 °C, which was not an input.

## 7. The figure

`item22_f32_residual_vs_velocity.pdf` / `.png`, built by `make_figure.py`,
caption in `item22_f32_residual_vs_velocity.caption.txt`:

> Window residual of the constant-rate leg against the swept gas velocity, for
> the two drying curves of Faner et al. (2019). Columns: soybean at 120 °C
> (a, c) and sunflower at 100 °C (b, d); the top row spans the whole swept range
> and the bottom row the same data around the acceptance band. Each line is one
> bed depth in particle diameters at the source's printed lower initial loading
> of 0.45 kg/kg, with the ribbon reaching its printed upper value of
> 0.50 kg/kg. The grey strip is the acceptance band, which is the residual this
> same closure reaches on the two Faner (2008) thesis runs, and the black
> diamond is the item-5 march of the journal trace placed on this velocity axis
> by its own film coefficient. Sample mass does not appear because it leaves the
> residual bit-identical.

## 8. Sentences the supplement could carry

Offered for §S10.9 (`sup:fullmarch`), which already says the journal comparison
is under-determined by its own source. These make that statement quantitative.
They are drafts for the lead; no `.tex` file is edited by this record.

> A sensitivity sweep over the four experimental quantities the 2019 article
> does not tabulate per run — sample mass, gas velocity, bed depth and initial
> loading, 960 scored grid points at the independent film closure — separates
> them. The sample mass cannot enter: the transfer area per unit dry mass of a
> layer of given depth is independent of it, and the swept residuals are
> bit-identical across a fourfold change. The initial loading is fixed to within
> digitization error by the first point of the source's own curve, 0.491 and
> 0.499 kg/kg, both inside the range the article prints. Only the bed depth and
> the gas velocity remain open, and they are what the comparison is sensitive
> to: over their swept ranges they move the constant-rate residual by roughly
> ten and six times the width of the thesis-run band respectively, against under
> two for the initial loading and exactly nothing for the sample mass.

> Each journal trace can be brought inside the thesis-run band by some point of
> that box, but not by the same point. Within the depth and loading ranges the
> article prints, no coordinate places both traces inside the band on its signed
> reading, and the coordinates that place both within its magnitude — a gas
> velocity of 0.20 m/s at three particle diameters, or 0.30 m/s at four —
> over-predict the soybean trace while under-predicting the sunflower one,
> whereas both thesis runs are under-predicted. The two published journal curves
> therefore demand different unprinted conditions from one another, which is why
> the gas velocity of the Figure 2 runs is the single most useful quantity to
> request from the authors.

The second paragraph strengthens the concrete ask already drafted in
`FANER_EMAIL_DRAFT_2026-09-24.md`: of the eight quantities that draft requests,
this measurement ranks the gas velocity and the per-run bed depth first and
second, and shows that the sample mass and the initial loading would add nothing
to this particular comparison.

## 9. Claim and non-claim

**Claimed.** That 48 farm jobs ran to `rc=0` at one commit and were collected
with every artifact digest verified against its receipt; that the runner's own
score over the declared 4-D box takes the values tabulated here; that the sample
mass axis is exactly degenerate and the initial-loading axis effectively so;
that the bed depth and gas velocity axes each move the score by several band
widths; that both traces meet the band somewhere in the box and that no common
coordinate meets it on the signed reading; and the three cross-checks of
section 6.

**Not claimed.** This is **not a fit and not a calibration.** Every coordinate
scored was handed to the runner by the campaign file; none was searched for.
Naming the best-scoring grid point is an identification of *which unprinted
quantity the score is sensitive to*, not a statement that the experiment ran at
that point, and the existence of a crossing follows from continuity rather than
from agreement. Nothing here is authenticated plant data, calibration, physical
qualification, production release or QSC-10. It does not supersede
`item_05_faner_comparison/RESULT_2026-09-22.md`, does not move any constant,
tolerance, gate, budget or seed family, and does not complete drying-rate
validation. No `.tex` file, no pinned kernel and no existing record was edited.
No plant, company or site appears anywhere in this record or its evidence.

---

Files: `DETAIL.md` (method, verification and the full tables),
`collect_f32.py`, `analyze_f32.py`, `make_figure.py`, `outputs/` (the 48 sweep
JSONs as they came back, 504 416 bytes in total), `f32_grid.csv` (960 rows),
`f32_trace_surface.csv` (240), `f32_job_summary.csv` (48),
`f32_common_conditions.csv`, `f32_footprint.csv`, `f32_summary.json`, the
figure and its caption. Evidence sidecar:
`docs/evidence/farm_f32_2026-09-25/`. Collection record:
`docs/GT_PS2_F32_COLLECTION_RECORD_2026-09-25.md`.
