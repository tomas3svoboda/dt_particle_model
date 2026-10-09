# Item 39: the journal conditions' constant-rate range ends on version 2, item 29's outcome S2 re-evaluated, and item 5's bracketing pre-registration quoted (2026-09-27)

`fitted: false. physically_qualifying: false. plant_predictive: false. qsc10_complete: false.`

Brief: close three loose ends left by `paper/01_particle_jfpe/CORRECTIONS_FOLDIN_2026-09-27.md`
section 1 (items 2, 3 and 5). Only this folder was written. No paper file, no
`docs/` file and no other result item was edited, and no git operation was run.
Every number below was printed by a script in this folder (see DETAIL.md).

---

## 1. The four-condition constant-rate ranges with the journal conditions on version 2

**Method.** item 5's own module (`item_05_faner_comparison/rerun_2026-09-22_rig/run_rerun_faner.py`)
was imported read only. Its readers were re-pointed at the version-2 curves.
Its stage-1 block (steps A and H, the measured window, what the rig requires,
the twelve film cases) was run verbatim. Stage 2 and the floor were left out
because they enter no range.

**Identity first.**
- Version 1, all four conditions: **1,466 leaves exact** by float equality
  against the committed `rerun_faner.json`.
- `conditions.csv`, v1 against v2: the 15 keys stage 1 reads are identical
  (sample mass 30 g, holder area 0.011 m2, Fig. 2 gas temperatures, diameters,
  porosities, densities, oil contents, critical loading). Only the curves move.
- Thesis conditions on version 2: 122 leaves equal item 35's values exactly.
- Journal conditions on version 2: 18 inputs agree with item 37 to 1e-14.
- One caveat, stated because it is real. The current tree's
  `source_correlated_film_pair()` returns 103.736 W/m2K (Re_eps 30.62), not the
  131.337 (Re_eps 45.93) the 2026-09-22 record was built on. That value enters
  only the borrowed `frozen_1.0` case and each case's film multiplier. Neither
  appears in any printed range. The committed value was used, which is why the
  identity holds. This change in `src/` is outside this item and was not
  traced here.

**What moves at the journal conditions (v1 -> v2).** The Whitaker, Bird and
Bradshaw-Myers coefficients do not depend on the curves and stay the same.

| quantity | 2019 soybean 120 C | 2019 sunflower 100 C |
|---|---|---|
| X0 | 0.490628 -> 0.492364 | 0.499279 -> 0.50164 |
| measured constant rate (1/s) | 7.2643e-3 -> 7.2606e-3 | 4.9904e-3 -> 5.0481e-3 |
| measured crossing (s) | 45.457 -> 45.864 | 65.085 -> 65.058 |
| specific duty (W/kg dry) | 2433.0 -> 2431.8 | 1671.4 -> 1690.8 |
| demanded h, bed, charge reading (W/m2K) | 86.52 -> 88.31 | 91.61 -> 92.93 |
| demanded h, dry-meal reading | 128.97 -> 131.79 | 137.35 -> 139.55 |
| demanded h, frozen particle area | 16.18 -> 16.53 | 17.23 -> 17.50 |
| Whitaker / demand | 0.909 -> 0.891 | 0.873 -> 0.861 |
| Bird / demand; Bradshaw-Myers / demand | 2.084 -> 2.042; 2.290 -> 2.244 | 1.777 -> 1.751; 2.199 -> 2.168 |
| duty ratio, Whitaker, charge | 0.9061 -> 0.9076 | 0.8233 -> 0.8152 |
| crossing ratio, Whitaker, charge | 0.9735 -> 0.9695 | 1.1216 -> 1.1291 |
| duty ratio, Whitaker, dry-meal | 0.6079 -> 0.6082 | 0.5491 -> 0.5429 |
| crossing ratio, Whitaker, dry-meal | 1.4511 -> 1.4468 | 1.6816 -> 1.6955 |
| slope match, closest miss, charge | 0.88 % -> 0.72 % | 1.13 % -> 0.13 % |
| the three bracket the demand | yes -> yes | yes -> yes |

**The ranges, before and after.** The thesis columns are item 35's version-2
values. "Printed now" is what the paper carries today: thesis on version 2,
journal on version 1.

| printed range | before (v1 everywhere) | printed now | **after (v2 everywhere)** | printed digits change? |
|---|---|---|---|---|
| duty, Whitaker, charge | 0.823-0.982 | 0.823-1.057 | **0.815-1.057** | no: 0.82 to 1.06 |
| crossing, Whitaker, charge | 0.894-1.122 | 0.857-1.122 | **0.857-1.129** | **yes: 0.86 to 1.12 -> 0.86 to 1.13** |
| duty, dry-meal reading | 0.508-0.608 | 0.530-0.609 | **0.530-0.609** | no: 0.53 to 0.61 |
| crossing, dry-meal reading | 1.451-1.745 | 1.451-1.682 | **1.447-1.695** | **yes: 1.45 to 1.68 -> 1.45 to 1.70** |
| Whitaker / demand | 0.873-0.998 | 0.873-1.074 | **0.861-1.074** | **yes: 0.87 to 1.07 -> 0.86 to 1.07** |
| Bird and Bradshaw-Myers / demand | 1.777-2.548 | 1.777-2.743 | **1.751-2.743** | no: 1.8 to 2.7 |
| demanded h, charge (W/m2K) | 85.5-97.1 | 81.7-91.6 | **81.7-92.9** | **yes, top end 91.6 -> 92.9** |
| demanded h, dry-meal (W/m2K) | 129.0-168.5 | 129.0-156.5 | **131.8-156.5** | **yes: 129.0 -> 131.8** |
| demanded h, particle area, S13.2 | 16.18-18.65 | 16.18-17.82 | **16.47-17.82** | **yes: 16.18 -> 16.47** (the low end is now the 2008 sunflower, 16.469) |
| specific duty (W/kg dry) | 1671.4-3682.6 | 1671.4-3532.1 | **1690.8-3532.1** | **yes: 1671.4 -> 1690.8** |
| slope match, closest miss, charge | 0.9-10.2 % | 0.7-10.0 % | **0.1-10.0 %** | **yes: "1 to 10 %" -> "0.1 to 10 %"** |

Per-condition duty ratios (Sec. 5.1): 1.06, 0.92, 0.91 and 0.82. These are
unchanged at two decimals (2019 sunflower 0.8233 -> 0.8152). The journal
demands 86.5 and 91.6 become 88.3 and 92.9. Item S13.2's journal value 17.228
becomes 17.504 (sunflower 100 C); the soybean 120 C value becomes 16.531.

**Bracketing on the executed reading (demand).** The three correlations
bracket the demand at both journal conditions, before and after. At the 2008
sunflower condition the result "all three above" (Whitaker at 1.074 times the
demand) comes from item 35 and is unchanged. So on this reading the outcome
fires at one of four conditions, both as printed now and after this item.

---

## 2. item 29's rejecting outcome S2 on the version-2 scoring

**S2 as written** (`item_29_faner_geometry_and_temperature_sensitivity/RESULT.md`,
lines 23 to 30):

> If the march depended on R and D only through R^2/D, the residual at radius
> factor f_R and D factor 1 would equal the residual at radius 1 and D factor
> f_R^-2. Fires if, at any radius factor where both are measured (the D-axis
> interpolated linearly in ln f_D), they differ by more than 0.01 of the window
> span.

**Re-evaluation.** Item 29's own construction was applied to item 35's
re-scored cells. Version 1 reproduces item 29's `analysis.json`, and version 2
reproduces item 35's shadow analysis: 114 checks exact.

| trace | v1 | **v2, rule window (14 points; the paper's window)** | v2, 15-point window |
|---|---|---|---|
| sunflower, largest difference (at f_R = 0.65) | 0.00926 | **0.01463** | 0.00648 |
| soybean, largest difference (at f_R = 0.60) | 0.00707 | 0.00652 | 0.00652 |
| S2 by its letter | does not fire | **FIRES (sunflower)** | does not fire |

**What the firing is.** Linear interpolation in ln f_D sits between the D
factors 2.0 and 2.5, and the sunflower residual has its minimum between them.
S2's pre-registration named other expected breakers: the real-time
temperature path, the stride floor and ceiling, and the output decimation.
Item 29's RESULT attributed its v1 value (0.0093) to the interpolation but did
not measure that. This item measured it directly. Each radius rung was re-marched against its
exact partner, radius 1.0 at the D factor f_R^-2. That gives seven new pairs
per trace (14 marches); item 29 had marched only 0.5 <-> 4.0. The redirected
driver first reproduced item 29's committed `r0.65` cells: **16,716 checks
exact**.

- Every pair takes the same number of strides.
- Stride lengths agree to 7.1e-15 s.
- The mean-loading path agrees to 5.4e-13 relative.
- The residual difference is at most **2.8e-14** of the span on v2, rule
  window (2.6e-14 full-step), 2.0e-14 on v1 and 1.6e-14 on v2 with 15 points.
- At f_R = 0.65 sunflower (v2, rule window) the rung and its partner both give
  0.06033. The linearly interpolated value is 0.0750, so the whole 0.0146 is
  interpolation error.
- None of the named breakers acts at this resolution: decimated and full-step
  scores agree between partners alike.

**Verdict.** S2 **fires by its letter** at version 2 under the paper's
window rule, on sunflower only, by 0.0046 above its threshold. It does **not**
fire on what it was written to test: the D/R^2 similarity holds to 3e-14 at
all eight radius rungs on both traces.

**What this means for the sentences the paper takes from item 29:**
- **"The residual is a function of D/R^2 alone, exactly."** Stands. It is now
  measured at eight exact pairs instead of one.
- **Radius lever per unit log change** (0.47 to 0.66 at v2) and the **sweep
  minima** (0.672 and 0.612 times 0.885 mm). These are computed directly on
  the radius axis without the D-axis interpolation, so the firing does not
  touch them. Both are already version-2 values in item 35.
- **The ratio -2.00 against D.** This comes from central differences on the
  grid, not the interpolation. It is unaffected.
- **The supplement's verdict line is wrong at version 2.** It says S2 "did not
  fire", marked as evaluated on version 1. Its similarity sentence already
  prints 0.0146 as "the interpolation's own error" without saying that 0.0146
  is above S2's own 0.01.

**Corrected sentence** (supplement, item 29's verdict paragraph; replaces the
S2 clause):

> S2 (the $D_\mathrm{eff}/R_p^2$ similarity, tested against the factor axis
> interpolated linearly in $\ln f_D$, differs by more than 0.01 of the span)
> fires at version 2 on the sunflower trace, 0.0146 at the radius factor 0.65
> next to the minimum, and not on the soybean trace (0.0065); marched exact
> partners at all eight radius factors reproduce the radius runs to
> $3\times10^{-14}$ of the span, the same strides and the same loading path to
> $5\times10^{-13}$, so the firing is the linear interpolation's error on a
> curved function and the similarity itself holds (item 39).

Similarity sentence: replace "to within 0.0146 of the span, the
interpolation's own error" with "to within 0.0146 of the span, above the 0.01
S2 was written with; exact partners marched at every radius factor agree to
$3\times10^{-14}$, so the difference is the interpolation's error (item 39)".

---

## 3. item 5's pre-registered bracketing outcome, quoted

### The text as recorded

**(a) The pre-registration.** `paper/01_particle_jfpe/LIMITS_AUDIT_2026-09-22.md`,
line 101: the M2 row of the ranked table, whose column header on line 98 reads
"Rejecting outcome, written first". The working-tree file is modified and
uncommitted; sha256 3708ff9c...0da0. The row's rejecting-outcome cell, verbatim:

> If the three correlations at the rig's state do not bracket the source's measured coefficients, the area or the group convention is wrong and the stage-1 posing is withdrawn

The same row's "Effect" cell names the leg:

> `sec_validation` 5.1's "That agreement is not a second prediction" gains an independent column; limits 7 and 8 are answered rather than declared

Limits 7 and 8 (same file, lines 38 and 39) both read "Bears on a claim: yes:
the constant-rate leg".

**(b) The execution record.** `paper/analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig/DETAIL.md`,
lines 109 to 111:

> It was written first: "if the three correlations at the rig's state do not
> bracket the source's measured coefficients, the area or the group convention is
> wrong and the stage-1 posing is withdrawn." They do bracket it.

The record then tests it against "the coefficient the measured duty demands on
the printed footprint". On lines 116 to 118 it adds:

> the tabulated mean lies inside the
> three-correlation band at the two thesis conditions and 3.8 % and 4.1 % below
> its low edge at the two journal conditions.

**(c) The summary record.** `paper/analysis/results_2026-09-19/item_05_faner_comparison/RESULT_2026-09-22.md`,
lines 84 to 86:

> **The three bracket the demand at all four conditions**, Whitaker below and the
> other two above, so the rejecting outcome — that a failure to bracket would
> withdraw the stage-1 posing — did not fire.

No separate pre-registration file exists in `rerun_2026-09-22_rig/`. The audit
row in (a) is the text written first.

### What the numbers say under each reading (measured here)

| reading | 2008 sunflower | 2008 soybean | 2019 soybean | 2019 sunflower | fires? |
|---|---|---|---|---|---|
| (i) the words written first: "the source's measured coefficients", the tabulated same-rig central value, which does not depend on the digitization | 97.27 inside 96.89-247.45 | 94.25 inside | **75.73 below 78.65** (3.7 %) | **76.67 below 80.01** (4.2 %) | **yes, at both journal conditions, on v1 and v2 alike** |
| (ii) as executed: the demand the plotted curve implies, v1 | 97.12 inside | 85.50 inside | 86.52 inside | 91.61 inside | no |
| (ii) as executed, v2 | **90.20 below 96.89** (Whitaker 1.074x) | 81.70 inside | 88.31 inside | 92.93 inside | **yes, at the 2008 sunflower condition** |

Under reading (i), the matched same-rig band at the 2019 sunflower condition,
73.9 to 78.6, lies wholly below Whitaker's 80.0. At the 2019 soybean condition
the band, 58.2 to 103.9, straddles it.

Under reading (ii), the digitization correction alone moved the 2008
sunflower demand by 7.1 % (97.1 -> 90.2). The firing margin is 7.4 %. The
criterion is therefore being decided at the resolution of the source's figure.

---

**The following is my own reading, not the record's.**

**What "withdrawing the stage-1 posing" would remove.** Stage 1 is the
constant-rate leg as posed: the film coefficient read at the rig's state,
applied on the sample-holder footprint per kg dry meal, running to the model's
critical loading. Withdrawing it would remove:
- the abstract's sentence "0.82 to 1.06 of the measured constant-rate duty and
  0.86 to 1.12 of the transition time";
- Table 2's second row, where Whitaker's correlation is labelled a measured
  trajectory;
- Sec. 5.1's duty, crossing and coefficient ratios, and the dry-meal-reading
  ranges;
- the matching sentences in the discussion and conclusions;
- the S10.3 and S10.8 material that is presented as a result rather than a
  disclosure.

The falling-rate comparison would lose its time origin, because every marched
tail is grafted on at the stage-1 end time t_c. That affects the 17.2 and
15.6 %, the crossing row 0.86 and 0.95 of Table 3, and item 29's sweep. None
of these could be printed as they are. They would need either a re-posed
origin, for example starting the tail at the measured crossing (a new,
different comparison), or withdrawal with the leg. What would survive:
- the property check of the constant-rate solid temperature;
- the uptake identification;
- the floor statements;
- the rig-state reconstruction as a disclosure.

**What a declared departure would have to say.** It would have to:
1. state which reading of "the source's measured coefficients" governs, and
   that the executed reading (ii) is not the reading written first (i);
2. state that under (i) the outcome had already fired on 2026-09-22, at both
   journal conditions, and was reported but not called firing, and that under
   (ii) it fired at version 2 at the 2008 sunflower condition;
3. give the reason for keeping the posing anyway, in numbers, without
   loosening the test after the fact. For example: the misses are 3.7 to 7.4 %
   against a three-correlation spread of a factor of 2.5; the reading (ii)
   margin is smaller than the move the re-digitization alone produced; and the
   consequence the audit named, a wrong area or group convention, is tested
   separately and bounded (M3, M7);
4. say that the departure was decided after the outcome was seen, by an owner
   ruling (Class B, a written consequence of a physics comparison), in a dated
   record;
5. say that the constant-rate sentence carries it wherever it is printed;
6. leave the pre-registration text and this item unaltered.

---

## 4. Claim and non-claim

**Claimed:**
- item 5's journal-condition stage 1 recomputed on version 2 by its own code,
  after 1,466 exact version-1 leaves;
- the printed four-condition ranges before and after;
- S2's letter fires at version 2 on sunflower (0.0146 > 0.01), while the
  similarity it tests holds to 3e-14 at eight marched exact pairs;
- the pre-registration and its execution quoted with file and line.

**Not claimed:**
- a ruling on the M2 consequence (the owner's);
- a re-march of any stage-2 journal trace (none enters a printed range);
- any trace of why the current tree's frozen film pair differs from the
  committed one;
- calibration, validation, qualification or plant prediction.
