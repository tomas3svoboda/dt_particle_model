# Item 23 — the oil-solution arm at a source-pinned activity coefficient

**Executes** item **B1** of `paper/01_particle_jfpe/SUBMISSION_TODO_2026-09-24.md`
and mitigation **M6** / limit **17** of
`paper/01_particle_jfpe/LIMITS_AUDIT_2026-09-22.md`.

**Outcome: the limit closes at its source.** The activity coefficient of
n-hexane in vegetable oil, which this project did not hold, is now pinned from
four independent published determinations, and the oil molar mass leaves the
construction entirely. All four rejecting outcomes declared in
`DECLARATION_2026-09-25.md` are clear. The arm roughly doubles and still does
not close the residual-floor gap: the measured final loading is **3.35** and
**1.77** times the floor at the two journal conditions instead of 6.09 and
2.79, which closes **33** and **45 per cent** of the gap in the logarithm
instead of 24 and 37.

**Two things the closure does not buy, and they are the new limit.** First,
every condition at which the arm is evaluated lies **above the concentration
limit below which the sources state Henry's law applies** — 1.5 to 2.5 times
above it — so the pinned coefficient is extrapolated there, and the one
measured isotherm that covers a journal condition says the extrapolation reads
**0.67 of** the measured equilibrium. Second, the overlap with the total-meal
isotherm is unchanged in kind and larger in size: **23 to 72 per cent** of the
measured retention at that isotherm's own conditions, against 17 to 53 before.
The arm therefore stays the *upper* bound of a bracket whose lower bound is
the arm-free floor.

No plant data enter. No parameter is fitted. Nothing under `src/` was edited.

**Method.** `run_item23_oil_arm.py` (about ten seconds) over
`oil_activity_coefficient.py` beside it, which extends — and does not edit —
`item_05_faner_comparison/rerun_2026-09-21/oil_solution_arm.py`. The four
sources are digitized in `paper/analysis/datasets/hexane_in_oil_activity/`
and checked against the owner's stated reading in `source_reading.md` there.
Outputs: `item23_oil_arm.json` and seven CSV files in this folder.

---

## 1. The coefficient

Four independent determinations, each reduced to the weight basis using **only
quantities its own source prints**, evaluated across the desolventizer band by
each source's own temperature form:

| Ω_w^∞ | 100 °C | 105 °C | 110 °C | 115 °C | 120 °C |
|---|---:|---:|---:|---:|---:|
| Smith & Wechter 1950, static vapour pressure | 5.798 | 5.662 | 5.533 | 5.411 | 5.295 |
| King & List 1990, inverse gas chromatography | 5.804 | 5.789 | 5.775 | 5.761 | 5.748 |
| Belting et al. 2014, dilutor | 4.741 | 4.709 | 4.679 | 4.649 | 4.621 |
| Belting et al. 2015, static cell via its own UNIQUAC | 4.856 | 4.812 | 4.769 | 4.727 | 4.684 |
| **mean** | **5.300** | **5.243** | **5.189** | **5.137** | **5.087** |
| full spread | 20.1 % | 20.6 % | 21.1 % | 21.6 % | 22.2 % |

Envelope over the band **4.621 to 5.804**, mean **5.191**, full spread
**22.8 per cent**.

| declared | value |
|---|---|
| Ω_w^∞ | **5.2**, bracket **4.6 to 5.9** |
| as a Henry constant, `p_hex = H_w · w_hex` | **14.4 atm = 1.46 MPa** at 105 °C (bracket 12.7 to 16.3 atm) |
| as a mole-basis coefficient | **0.51** with the sources' own oil molar masses (870 to 876 g/mol) |
| oil molar mass | **no longer declared** — it leaves the construction |

The disagreement is not scatter about a common value. It is a systematic split
into two camps, the two older determinations at 5.3 to 5.8 and the two Belting
determinations at 4.6 to 4.9, and nothing in the four sources decides between
them. The declared value splits it and the bracket covers both.

### Two corrections to the print, both measured

**King & List Table 2 at 123.4 °C.** The printed 4.77 disagrees with three
other printed quantities of the same paper: Table 3's 0.572 implies 5.775,
Table 4's 24.5 atm implies 5.709, and Table 5's heat of mixing
0.14 ± 0.06 kcal/mol is returned by 5.77 (van't Hoff slope 0.141 kcal/mol) and
not by 4.77 (0.818 kcal/mol, R² 0.70). Read as 5.77 the paper's weight-basis
coefficient is flat to 3.4 per cent over its whole 65 K range, which is what
its own "athermal" reading of the alkanes says it should be.
`item23_kinglist_internal_consistency.csv` carries all four checks.

**Belting 2015 Table 11's c column.** It is placed differently by two text
extractors and neither placement reproduces the source's own data. Of the 361
assignments of the eighteen printed c values to the two rows, one reproduces
the printed VLE table to 1.92 per cent (printed 1.69) and the printed
infinite-dilution coefficients to 3.36 per cent (printed 1.23); the runner-up
is 2.10 and 2.37. The reconstruction is used only as the fourth method and
never as a primary pin.

---

## 2. The rejecting outcomes

| | criterion | measured | verdict |
|---|---|---|---|
| **R1** | any of the four methods outside 4.6 to 5.9 over 100–120 °C | envelope 4.621 to 5.804 | **clear** |
| **R2** | `H_w` at the declared value differs from the printed Henry constants by more than the four-method spread (22.8 %) | −8.7 % against 14.2 atm at 100.9 °C; −13.5 % against 15.5 atm at 102.2 °C | **clear** |
| **R3** | the arm plus the arm-free floor exceeds the measured final loading at either journal condition | largest floor 4.423e-3 against 7.415e-3 kg/kg | **clear** |
| **R4** | the arm exceeds 100 % of the measured total-meal retention at that isotherm's own conditions | 72.0 % at the declared value, 81.7 % at the bracket's high arm | **clear** |

R2 is the sharpest of the four and deserves its number in full: the pinned
coefficient reads **below** both printed Henry constants, by about half the
four-method spread. That is the expected consequence of splitting a bimodal
disagreement whose two printed Henry constants both sit in the upper camp.

---

## 3. The arm, and the floor

`item23_floors.csv`. The arm-free floors are recomputed from the certified
dry-shell storage functional and reproduce
`item_05_faner_comparison/rerun_2026-09-21/rerun_floor.csv` to 1e-15 kg/kg
before anything is built on them.

| condition | *T* | `a_h` | arm-free floor | **floor with the arm** | bracket | measured / floor | **log gap closure** |
|---|---:|---:|---:|---:|---|---:|---:|
| journal, soybean | 120 °C | 0.2535 | 1.2169e-3 | **2.2161e-3** | 2.0922e-3 – 2.3540e-3 | 6.09 → **3.35** | **33.2 %** (30.0 – 36.5) |
| journal, sunflower | 100 °C | 0.4114 | 2.6546e-3 | **4.2010e-3** | 4.0038e-3 – 4.4225e-3 | 2.79 → **1.77** | **44.7 %** (40.0 – 49.7) |
| thesis, soybean | 136 °C | 0.1786 | 7.6072e-4 | 1.4544e-3 | 1.3695e-3 – 1.5485e-3 | 132.2 → 69.2 | 13.3 % |
| thesis, sunflower | 136 °C | 0.1786 | 7.6072e-4 | 1.4010e-3 | 1.3227e-3 – 1.4879e-3 | 136.9 → 74.4 | 12.4 % |
| industrial vessel | 105 °C | 0.2219 | 1216.4 ppm | **2085.5 ppm** | 1978.4 – 2204.6 ppm | — | — |

All loadings are kg n-hexane per kg dry meal. The residual oil is 1.95 per
cent by mass for the soybean rows and 1.80 for the sunflower ones, as the
sources report. The two 136 °C rows are **outside the declared 100–120 °C
domain** and are printed for continuity only; nothing in the manuscript rests
on them, and at those conditions the floor is under 1.5 per cent of the
measured loading either way.

The arm closes its ledger to exactly 0.0 at all three activities, vanishes to
exactly 0.0 at zero residual oil, and is linear in the oil fraction to one
unit in the last place (measured slope spread 1.6e-16 relative).

---

## 4. The new limit: the arm runs above the Henry region

Smith's Fig. 4 prints the concentration above which his own Henry constant
stops applying. Converted to an activity through his own constant:

| *T* | Henry limit, `a_h` | the arm's condition | ratio |
|---|---:|---:|---:|
| 100 °C | 0.163 | 0.4114 | **2.53 ×** |
| 105 °C | 0.150 | 0.2219 | **1.48 ×** |
| 120 °C | 0.104 | 0.2535 | **2.43 ×** |
| 136 °C | 0.098 | 0.1786 | 1.83 × |

**Every condition is above it.** The pinned coefficient is therefore
extrapolated at every condition the paper reports, and the two linear forms of
that extrapolation — weight basis and mole basis — diverge there. Belting
2015's measured isotherm at 373.15 K settles which way each errs, because it
covers the 100 °C journal condition exactly, at its own temperature and its
own activity, with no model of any kind (`item23_extrapolation_test.csv`):

| `a_h` at 373.15 K | weight basis / measured | mole basis / measured |
|---:|---:|---:|
| 0.050 (inside the Henry region) | 0.83 | 0.91 |
| 0.150 (at the limit) | 0.81 | 1.10 |
| 0.2535 | 0.76 | 1.41 |
| **0.4114 (the journal condition)** | **0.67** | **2.89** |

The measured mole fraction at the 100 °C journal condition is **0.5615**,
against 0.439 from the weight-basis reading and 0.802 from the mole-basis one.

Read straight off that measured isotherm, the arm at the 100 °C condition is
**2.2729e-3 kg/kg**, the floor **4.9275e-3**, the measured loading **1.50**
times it, and the log gap closure **60.2 per cent**. At the same pinned
coefficient the mole-basis reading gives an arm of 7.207e-3 kg/kg, a floor of
**9.86e-3 — above the measured final loading of 7.415e-3** — which is the
falsification R3 exists to catch, and it is why the weight basis is the form
carried.

**What this means for the manuscript.** The pinned arm is a *lower* reading of
the oil term at these activities and the measured equilibrium is 1.46 times
larger. That does not disturb the arm's role as the *upper* bound of the floor
bracket, because that bracket is about the overlap with the total-meal
isotherm and not about the coefficient. It does mean the arm is not tight, and
it is stated so.

---

## 5. The overlap, unchanged in kind and larger in size

`item23_overlap.csv`. Measured against the frozen total-meal retention
isotherm at that isotherm's own measured conditions (50 to 95 °C, activities
0.179 to 0.8):

| | previous, ideal solution | now, at the pinned coefficient |
|---|---:|---:|
| arm / measured native retention | 17 – 53 % | **23 – 72 %** |
| at the bracket's high arm (Ω_w = 4.6) | — | up to **81.7 %** |

The total-meal isotherm was fitted on meal that already carried its residual
oil, so this fraction is what the arm could be duplicating. It does not reach
100 per cent, so R4 is clear and the arm may be added; it is now within a
third of doing so, which is why the corrected floor stays an upper bound and
no single corrected number is presented as the model's floor.

---

## 6. The sentences the paper can carry

Each is checked against `item23_oil_arm.json` by `check_numbers_2026-09-23.py`.

- The arm's activity coefficient is now source-pinned: the weight-fraction
  activity coefficient of n-hexane in refined vegetable oil at infinite
  dilution is **5.2, bracketed 4.6 to 5.9** by four independent published
  determinations over 100 to 120 °C, equivalently a Henry constant
  `p = H_w w` with `H_w` = **14.4 atm at 105 °C**; the oil molar mass leaves
  the construction.
- At the two journal conditions the arm raises the floor from
  1.217e-3 to **2.216e-3** kg/kg at 120 °C and from 2.655e-3 to **4.201e-3**
  kg/kg at 100 °C, so the measured final loading of 7.415e-3 kg/kg is **3.35**
  and **1.77** times above it instead of 6.09 and 2.79. In the logarithm that
  closes **33** and **45 per cent** of the gap and does not close it.
- At the industrial-vessel vapour composition the same construction floors at
  **2085.5 ppm** instead of 1216.4.
- The coefficient is measured; the *extrapolation* is not. All four conditions
  lie 1.5 to 2.5 times above the concentration limit below which the sources
  state Henry's law applies, and the one measured isotherm that reaches a
  journal condition puts the pinned arm at **0.67** of the measured
  equilibrium there.
- The arm remains the upper bound of a bracket whose lower bound is the
  arm-free floor, because the total-meal isotherm was measured on meal that
  already held its residual oil; at that isotherm's own conditions the arm is
  **23 to 72 per cent** of the measured retention, which sizes the possible
  overlap. That overlap, and not the coefficient, is now the unmeasured term.

---

## 7. For the lead, not for the paper: the tower model's own oil closure

Read-only; `src` was not edited. `item23_tower_oil_closure.csv` and §5 of
`DETAIL.md`. `src/dtdc_simulator/core2/props/sorption.py` carries PHY-048's
`OilIsotherm`, a power law `q_o(a) = 0.9635 · a^2.7036` in kg hexane per kg
oil, with **no temperature dependence**, composed as
`W_total = α_o·W_native + β_o·q_o(a)` with
`β_o = (w_o − 0.0195)/(1 − 0.0195)`.

Two findings, either of which is a Class-B question for the owner:

1. **At the reference oil fraction the closure contributes nothing.**
   `β_o` is exactly zero at `w_o = W_O_REF = 0.0195`, so `q_o` is a *signed
   off-anchor correction* and not a standing hexane-in-oil equilibrium. The
   tower model at its own reference composition carries no oil-solution term
   at all.
2. **Where it does contribute, it disagrees with the four sources by orders of
   magnitude at low activity.** At 105 °C, inside the Henry region, the
   sources' weight-basis equilibrium exceeds the power law by **506 ×** at
   `a_h` = 0.01, 156 × at 0.02, 33 × at 0.05, 10 × at 0.10 and 5.2 × at 0.15.
   The two cross near `a_h` ≈ 0.41 and the power law is the larger above it.
   **Every point inside the Henry region exceeds the four-method spread**, by
   one to nearly three orders of magnitude. The low-activity end is exactly
   where a desolventizer's last trays sit.

The absence of temperature dependence is the one thing the sources support:
they place the heat of mixing near zero and the weight-basis coefficient moves
4 per cent over 100 to 120 °C.
