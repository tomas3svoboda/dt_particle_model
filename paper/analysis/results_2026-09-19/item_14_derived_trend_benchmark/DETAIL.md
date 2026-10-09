# Item 14 — derived-trend benchmarking, steps 3 and 4 (2026-09-20)

`physically_qualifying: false. plant_predictive: false.`
No parameter was fitted, tuned, smoothed or blended. **A refused comparison is a
recorded result.** Not one pre-registered number was revised; `PRE_REGISTRATION.md`,
`derive_expectations.py` and `item14_expectations.json` are unmodified, and the
scorer asserts bit-for-bit reproduction of twelve of their numbers before it scores
anything.

Scripts, both run on this workstation under the pinned environment
(`PYTHONHASHSEED=1`, single-threaded BLAS/OMP/MKL/NumExpr):
`run_item14_with_water.py` (863.9 s, 120 certified dry-shell marches) and
`score_item14.py`. Data: `item14_simulated.json`, `item14_trend_table.csv`,
`item14_band_widths.csv`, `item14_traces.csv` (7 649 rows),
`item14_refused_register.csv`. Figure: `item14_figure.pdf` / `.png`.

---

## 0. A claim this pass made, checked, and corrected

The first execution of `run_item14_with_water.py` asserted in its docstring and in
`certified_lane.headline` that one step above `a_h = 1` the coupled-pore authority
refuses with *"activate the mobile-liquid topology, do not clamp"* — while its own
JSON recorded `free_hexane_surface_refusal: null`. The probe was defective: it
stepped `y_hexane` by `+1e-8` from the **Raoult** seed, where the certified state's
hexane activity is only 0.9977, so nothing was ever pushed above 1 and the
authority never spoke. The probe was rebuilt to start **on the certified plateau**,
where the certified state's hexane activity is 1 to within `1e-15`, and to step the
composition from there. The authority then does say it, and the sentence in section
2 below is that measurement and not a paraphrase. The correction is recorded here
rather than made silently, because the first text asserted a refusal that had not
been captured.

---

## 1. The anchor, and the two selection rules

External conditions are part C's `faner2008_soybean_136C` condition, read out of
its own `rerun_faner.json`: gas 409.15 K, 101 325 Pa (a READING — the source prints
no pressure), frozen P1EF geometry (`R = 0.885 mm`, `rho_dm,p = 1159.65`,
`eps_p = 0.141`, `w_o = 0.0195`), measured `X_0 = 0.7350938` kg/kg dry, measured
constant rate `1.0995122e-2` 1/s, measured `X_c = 0.20` crossing 54.930 s, last
measured loading 0.1005939 kg/kg at 300.119 s. **Only the water changes.**

The brief's step 1 asks for the declared interval END closest to the measurement,
chosen by a stated rule, never a fitted value. Both rules are recomputed in the
runner from part C's own CSVs:

| interval | rule | ends measured by part C | selected |
|---|---|---|---|
| film `COMPACT_FILM_STRESS_MULTIPLIERS` (0.5, 1.0, 2.0) | argmin of `\|ln(predicted/measured)\|` on the constant-rate demand | 3.5073 / 7.0147 / 14.0293 | **`h x 0.5` = 65.6687 W m⁻² K⁻¹** |
| `D_INTERVAL` from `run_item06_tail.py`, six ends | argmin of `\|ln(t_sim/t_meas)\|` for the time to the last measured loading, at the selected film end, sensible case | 2088.82 / 792.63 / 701.69 / 599.59 / **403.01** / 216.58 s against 300.119 s | **`phy019_frozen_4.0e-10`** |

Each is an argmin over a three- or six-point declared set. **Neither selection
touches a single scored number**: every trend below is a ratio against the anchor,
a temperature, or an equilibrium loading, and `h` cancels algebraically out of all
three. The declared film interval contributes exactly zero width to every row of
the trend table, and that is algebra, not an approximation. Part C's own section 6
finds the raw film band is not even the right comparison basis (mesh footprint
against particle area, ratio 7.8–9.7); that area ratio is a comparison basis there,
is never fed back as a coefficient there, and is not used here either.

---

## 2. Step 3, first preference: the certified coupled lane, probed three ways

**Adding the water RE-OPENS the certified lane at the state level.** All three
sorbed charges are admitted where the water-free anchor is refused.

| case | `coupled_pore.evaluate_equilibrium` at the case's plateau | two-equation certified plateau solve | `sphere.initialize_qualified_feed` at `X_0` |
|---|---|---|---|
| anchor, water free | **REFUSED** | **REFUSED** | **REFUSED** |
| sorbed 0.05 | admitted | admitted | admitted |
| sorbed 0.10 | admitted | admitted | admitted |
| sorbed 0.20 | admitted | admitted | admitted |
| free water film | **REFUSED** | **REFUSED** | **REFUSED** |

Verbatim, the three distinct refusals:

```
CoupledPoreTopologyError: coupled binary pore state requires 0 < y_hexane < 1
```
```
CoupledPoreTopologyError: dry-pore water activity must lie in (0,1], got
1.0196361579160484; activate the external/free-water topology, do not clamp
```
```
ValueError: feed retained water must lie on the qualified positive-moisture
Luikov branch [0.023, 0.235670904] inside the evidence cap
```

The free-film refusal is the correct one and not a defect: a free liquid film is by
definition not a dry-pore state, and the authority names the topology that owns it.

**The certified two-equation plateau solve.** For each sorbed charge, solve on
`coupled_pore.evaluate_equilibrium` ALONE for the `(T, y_h)` at which the certified
state's hexane activity is 1 (free hexane present) and its retained-water loading
is the charge. This carries the full binary virial fugacities and the frozen Luikov
and GAB branches; **it shares no arithmetic with the Raoult construction of
`derive_expectations.py`**, so its plateau is an independent determination.

| charge | pre-registered plateau (°C) | certified solve (°C) | difference (K) | certified `y_h` | Raoult `y_h` | residuals |
|---|---|---|---|---|---|---|
| 0.05 | 66.500879 | **66.507248** | **+0.0064** | 0.93465309 | 0.93235260 | `a_h - 1` = −4.4e-16, `W - 0.05` = 2.1e-15 |
| 0.10 | 64.448405 | **64.425080** | **−0.0233** | 0.87660260 | 0.87288923 | −8.9e-16, 6.1e-15 |
| 0.20 | 62.896497 | **62.850849** | **−0.0456** | 0.83433299 | 0.82994086 | −5.6e-16, 4.0e-14 |

Against a pre-registered depression of 2.21 / 4.27 / 5.82 K and a pressure-reading
band of ±0.33 K, the certified authority lands **0.006 to 0.046 K** from the
pre-registration. The whole residue is the declared Raoult approximation, and it is
visible in the composition, not the temperature: the certified vapour is 0.0023 to
0.0044 richer in hexane than Raoult.

**Second preference, and why.** The constant-rate LEG still cannot be marched on
that lane. A surface carrying free hexane sits exactly on `a_h = 1`; stepping the
composition one part in 10⁸ above it from the certified plateau gives

```
HexaneSupersaturationTopologyError: dry-pore hexane activity must lie in (0,1],
got 1.0000000099261663; activate the mobile-liquid topology, do not clamp
```

and the mobile-liquid (cut / moving-front) lane takes **no ordinary macrostep at
the P1EF production coefficients** — 50 probes, zero accepted, measured in
`docs/GT_PS2_PARTICLE_DRYING_VALIDATION_2026-09-03.md` section 6 — and is
additionally outside this pass's declared scope because another agent holds its
modules under edit. So part C's two-stage construction is used, with the water
entering only through the plateau temperature, the vapour composition there, the
duty split between the two latent heats, and the hexane partial pressure at the
surface. Every one of those comes from a frozen closure.

For the record, part C's pure-hexane carrier floor at 136 °C is reproduced:
`y_h` refused at 0.85700736 on the qualified Luikov lower-bound activity. The
surface vapour compositions of this item are 0.9324 / 0.8729 / 0.8299 / 0.7901, so
two of them sit above that 136 °C floor — but **every case is evaluated at its own
plateau**, where the floor is lower, and all three sorbed charges are admitted
there.

---

## 3. Step 3, the posed two-stage construction

Retained-water state: `sorption.water_activity` on the qualified Luikov band for
0.05, 0.10 and 0.20 (all three loadings inside `[0.023, 0.235670904]`),
`dry_limit_moisture_branch.dry_limit_water_activity(0.0) = 0.0` exactly for the
anchor, and `a_w = 1` by definition for a free film, which uses no isotherm.

| case | `a_w` | plateau (°C) | Δ (K) | `y_water` | rate ratio | model's own `X_c` | `t_c` ratio (own `X_c`) | `t_c` ratio (fixed `X_c`) | water-departure temperature (°C) |
|---|---|---|---|---|---|---|---|---|---|
| anchor | 0 | 68.714513 | — | 0 | 1 | 0.1993532 | 1 | 1 | — |
| 0.05 | 0.256035 | 66.500879 | −2.2136 | 0.06765 | 0.929703 | 0.2000712 | 1.074171 | 1.075612 | **143.223** |
| 0.10 | 0.527195 | 64.448405 | −4.2661 | 0.12711 | 0.869643 | 0.2007334 | 1.146934 | 1.149897 | **118.956** |
| 0.20 | 0.756499 | 62.896497 | −5.8180 | 0.17006 | 0.827125 | 0.2012320 | 1.204768 | 1.209008 | **107.985** |
| free film | 1 | 61.397575 | −7.3169 | 0.20994 | 0.788197 | 0.2017118 | 1.263133 | 1.268719 | **99.974** |

The fixed-`X_c` column reproduces the pre-registration exactly. The own-`X_c`
column is the model's own behaviour and it is the scored one: `X_c` RISES as the
plateau falls, so part of the slower rate is cancelled by having less free hexane
to remove. The water-departure temperature is where `a_w Psat_w(T) = P` on the
frozen IAPWS-95 surface; for 0.05 it is ABOVE the 136 °C gas and for the other
three it is below it, which is the pre-registration's own refusal expressed as a
temperature.

**The tail.** `particle.dry_shell` marched from `X_c` at each case's own plateau to
2400 s, at the surface activity `a_h = (P − p_w)/Psat_h(T_p)`, in two temperature
cases (`bounding_gas`, `sensible`) and two water constructions, over all six ends
of the declared `D_INTERVAL`: 120 marches. `fixed_charge_activity` holds
`p_w = a_w,charge Psat_w(T_p)` and REFUSES once `p_w >= P`.
`water_departure_switch` is a DECLARED CONSTRUCTION OF THIS PASS: the water is
either still on the particle at its charge activity or it has left, because the
model's own coexistence law forbids it to stay. It is a sharp switch and not a
capped activity, because capping at the coexistence ceiling is the degenerate limit
the derivation itself refuses to call a residual (`p_w = P`, hexane partial
pressure zero).

DECLARED LIMITATION of the `sensible` rows: the only latent sink integrated is the
solver's own hexane evaporation, because the certified dry shell carries no water
component. The water's own departure enthalpy is not charged against the particle,
so a with-water sensible tail heats at least as fast as the anchor's. The bias is
toward the anchor, and the `bounding_gas` rows — which integrate no temperature —
are free of it. Every scored tail number below is a `bounding_gas` row.

Frozen `D_eff`, loading at 2400 s:

| case | fixed activity, bounding | fixed activity, sensible | departure switch, bounding | departure switch, sensible |
|---|---|---|---|---|
| anchor | 9.552661e-3 | 9.995639e-3 | — | — |
| 0.05 | 5.988671e-3 | 6.391865e-3 | 5.988671e-3 | 6.391865e-3 |
| 0.10 | **REFUSED at start** | **REFUSED at 23.63 s**, 120.72 °C | 9.867115e-3 | 1.030625e-2 |
| 0.20 | **REFUSED at start** | **REFUSED at 16.34 s**, 111.10 °C | 9.981678e-3 | 1.043615e-2 |
| free film | **REFUSED at start** | **REFUSED at 11.15 s**, 100.36 °C | 1.009242e-2 | 1.056152e-2 |

The three `sensible` refusals land within one stride of the analytic departure
temperatures 118.956 / 107.985 / 99.974 °C, measured independently on the frozen
surface. That is the model finding its own coexistence limit by marching into it.

---

## 4. Step 4 — the band, and whether each test is weak

**The band.** The dominant declared uncertainty on every absolute number is the
operating pressure, which the source never prints. Atmospheric is a reading; the
pre-registration records `dP/dT ~ 3011 Pa/K` at the plateau, so **±1 kPa**, declared
as a reading, is propagated through EVERY expected number by re-evaluating
`derive_expectations.py`'s own functions at 100 325 and 102 325 Pa. Part C declared
no narrower band. The film interval enters nothing (section 1). The declared
`D_eff` interval enters only the marched-tail rows, which are run at all six of its
ends and report the spread.

**Is the band wide enough that the test is weak?** No — it is the opposite. The
half-width as a fraction of the predicted effect (`item14_band_widths.csv`):

| observable | band half-width / predicted effect |
|---|---|
| 1 plateau change | 0.52 % – 0.57 % |
| 2 constant-rate ratio | 0.079 % – 0.13 % |
| 3 crossing-time ratio | 0.100 % – 0.14 % |
| 4 residual floor ratio (0.05 only) | 0.87 % |

A ±1 kPa reading error moves the expectation by less than one part in a hundred of
the effect it predicts. **The band is narrow, so the magnitude test is strict.**
What is weak in some rows is not the band but the INDEPENDENCE, and the table says
per row which it is: the `simulated_is_independent_of_the_derivation` column is
`True` only where the simulated side comes from the certified coupled-pore
authority or from the model's own `X_c`, and `False` for the two free-film rows
(no certified determination exists, so the simulated side IS the derivation's own
construction and those two rows agree exactly and trivially) and for the
equilibrium-floor row (the same frozen isotherm at the same activity on both
sides).

### The full trend table

| water | observable | expected | ±1 kPa band | simulated | dev. of predicted effect | independent | agreement |
|---|---|---|---|---|---|---|---|
| 0.05 | 1 plateau Δ (°C) | −2.213634 | [−2.226261, −2.200957] | **−2.207265** | +0.29 % | yes | **FULFILLED** |
| 0.05 | 2 rate ratio | 0.929703 | [0.929614, 0.929795] | 0.932870 | +4.51 % | yes | sign only |
| 0.05 | 3 `t_c` ratio | 1.075612 | [1.075506, 1.075715] | 1.074171 | −1.91 % | yes | sign only |
| 0.05 | 4 residual ratio | 0.167319 | [0.159997, 0.174497] | **0.167319** | 0.00 % | no | **FULFILLED** |
| 0.05 | 4b marched residual | — | — | 0.626911 (D band [0.1747, 0.9609]) | — | no | supporting, sign agrees |
| 0.10 | 1 plateau Δ (°C) | −4.266108 | [−4.289484, −4.242633] | **−4.289432** | −0.55 % | yes | **FULFILLED** |
| 0.10 | 2 rate ratio | 0.869643 | [0.869508, 0.869784] | 0.874973 | +4.09 % | yes | sign only |
| 0.10 | 3 `t_c` ratio | 1.149897 | [1.149711, 1.150076] | 1.146934 | −1.98 % | yes | sign only |
| 0.10 | 4 residual ratio | **REFUSED** | — | **REFUSED** | — | no | **FULFILLED (refusal predicted and reproduced)** |
| 0.10 | 4b marched residual | no prediction | — | 1.032918 (D band [1.0000, 1.0329]) | — | no | not scored |
| 0.20 | 1 plateau Δ (°C) | −5.818016 | [−5.848975, −5.786920] | −5.863663 | −0.78 % | yes | sign only |
| 0.20 | 2 rate ratio | 0.827125 | [0.826970, 0.827286] | 0.833473 | +3.67 % | yes | sign only |
| 0.20 | 3 `t_c` ratio | 1.209008 | [1.208772, 1.209234] | 1.204768 | −2.03 % | yes | sign only |
| 0.20 | 4 residual ratio | **REFUSED** | — | **REFUSED** | — | no | **FULFILLED (refusal predicted and reproduced)** |
| 0.20 | 4b marched residual | no prediction | — | 1.044911 (D band [1.0000, 1.0449]) | — | no | not scored |
| film | 1 plateau Δ (°C) | −7.316938 | [−7.354819, −7.278882] | **−7.316938** | 0.00 % | no | **FULFILLED** |
| film | 2 rate ratio | 0.788197 | [0.788033, 0.788368] | **0.788197** | 0.00 % | no | **FULFILLED** |
| film | 3 `t_c` ratio | 1.268719 | [1.268443, 1.268982] | 1.263133 | −2.08 % | yes | sign only |
| film | 4 residual ratio | **REFUSED** | — | **REFUSED** | — | no | **FULFILLED (refusal predicted and reproduced)** |
| film | 4b marched residual | no prediction | — | 1.056504 (D band [1.0000, 1.0565]) | — | no | not scored |

**16 of 16 scored trends fulfil the expected SIGN. 8 of 16 also land inside the
±1 kPa band.** Every one of the eight magnitude misses has a single named cause,
and neither cause is a new physical claim:

* observable 2, +3.7 % to +4.5 % of the effect — the certified vapour composition
  is 0.0023 to 0.0044 richer in hexane than Raoult. That is the size of the
  derivation's **first declared approximation** (Raoult on both species; no
  activity coefficient for the hexane–water liquid pair exists in the repository),
  measured rather than argued;
* observable 1 at 0.20, −0.78 % — the same non-ideality, 0.046 K;
* observable 3, −1.91 % to −2.08 % of the effect on all four contents — the
  derivation holds `X_c` fixed by construction; the model's own PHY-050 critical
  loading RISES as the plateau falls, 0.1994 to 0.2017 kg/kg, which shortens the
  constant-rate leg and cancels about a fiftieth of the predicted delay. The
  fixed-`X_c` restatement reproduces the pre-registered value exactly, so this row
  measures a **missing term in the derivation, not a disagreement in the model.**

### Row 4b, and why it is not scored

The primary residual row is an equilibrium statement; 4b is a kinetic one and is
reported beside it, never folded into the count. For 0.05 the water survives to the
gas temperature (departure at 143.2 °C), the tail runs WITH water on the surface,
and its 2400 s loading falls to 0.627 of the anchor's — the pre-registered
direction, over a declared `D_eff` spread of [0.175, 0.961]. For 0.10, 0.20 and the
free film the pre-registration makes **no residual prediction at all** — its
prediction is the refusal — and the marched tail reproduces the mechanism behind
that refusal: the water departs at 119.0 / 108.0 / 100.0 °C, below the 136 °C gas,
so by the time the hexane tail runs the particle is water-free and the residual
returns to the anchor's. The ratios 1.033 / 1.045 / 1.057 are **not a water
effect**: they are the higher critical loading those charges start the tail from
(0.2007 / 0.2012 / 0.2017 against the anchor's 0.1994). Where the `D_eff` end makes
the tail slow, the ratio is exactly 1.000 — the two runs are then the same run.

---

## 5. The refused register — 17 elements, none skipped

`item14_refused_register.csv`. By class:

| class | count | content |
|---|---|---|
| certified coupled lane, state-level refusal | 2 | the anchor at `y_h = 1`; the free film at `a_w = 1.0196` |
| certified coupled lane, no plateau determination | 2 | the same two cases have no certified plateau |
| certified coupled lane, feed-entry refusal | 2 | the same two cases refuse `initialize_qualified_feed` |
| certified coupled lane, mobile-liquid topology required | 1 | the free-hexane surface, one step above `a_h = 1` |
| certified march lane, measured elsewhere | 1 | the cut / moving-front march at production coefficients |
| model-side law: the state does not exist | 3 | the equilibrium residual floor at 136 °C for 0.10, 0.20, free film |
| model-side law: the hexane tail cannot be posed while the water is still there | 6 | the `fixed_charge_activity` tails of those three charges, in both temperature cases |

The last nine are **the pre-registration's own prediction coming true**: its section
on the coexistence ceiling says the meal cannot hold more than 0.058262 kg/kg of
sorbed water at 136 °C and 1 atm, so for 0.10 and 0.20 "the model's own laws say
the water leaves before the hexane tail begins ... The refusal is the prediction."
Nine independent refusals, from three different code paths, say exactly that.

The three `fixed_charge_activity` `bounding_gas` refusals fire at the first
evaluation; the three `sensible` ones fire after the particle has marched up to its
own coexistence temperature, verbatim, e.g.

```
refused at t = 23.6310 s of the tail: with-water tail surface refused: the retained
water demands a partial pressure of 108626.7 Pa at 394.3042070758377 K, at or above
the total 101325.0 Pa; no n-hexane partial pressure remains, so the hexane tail
cannot be posed at this particle temperature
```

---

## 6. Errors found in the pre-registration

**None.** Twelve of its numbers were recomputed from its own script and are
bit-identical; its four declared gaps are all live and two of them are quantified
above rather than contradicted; its refusal predictions all fired. One
**incompleteness**, not an error, is recorded: its crossing-time derivation holds
`X_c` fixed, and the model does not — section 4, observable 3. That is stated in
the derivation itself ("so its ratio is the inverse of the rate ratio at the same
`X_0` and `X_c`"), so it is a declared scope of the derivation and the measured
size of the omission, about 2 % of the predicted effect, is this item's addition.

---

## 7. Claim boundary

This is a **consistency demonstration** between the model's coupled response and
the thermodynamics it is built on, against a pure-hexane measured anchor. **It is
not a validation against measured wet-meal drying curves.** The Cardarelli
chapter-3 uptake curves of item 13 — meal carrying its own moisture — are the
closest measured wet-meal comparison in this repository and must be
cross-referenced wherever a wet-meal claim is near. That item's own outcome bounds
what this one may borrow: on 13 independent measured series no curve enters the
two-reading band, the frozen `4.0e-10` is the honest like-for-like coefficient only
on the low-activity sunflower series, and its declared gaps include the absence of
any water branch in the uptake probe (`AttributeError: ... has no attribute
'water_retained'`, 8 refused water figures). Item 14 adds water to a hexane
trajectory; item 13 is where measured water-bearing meal actually lives, and the
two do not yet meet.

Further boundaries: the marched trajectory is the certified dry shell inside a
declared two-stage construction, not a certified moving-front solution; the
constant-rate leg is an energy balance and remains unmarched on the certified lane;
`water_departure_switch` and the sensible temperature integration are declared
constructions of this pass; the two domain extensions are typed non-qualifying
authorities; the film interval is a declared engineering interval; no plant datum
appears anywhere in this record.

`physically_qualifying: false. plant_predictive: false.`
