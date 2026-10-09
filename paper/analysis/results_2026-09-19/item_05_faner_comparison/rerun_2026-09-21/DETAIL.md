# Item 5, re-run of 2026-09-21 — the Faner march repaired on physical grounds

`physically_qualifying: false. plant_predictive: false.`
No coefficient is fitted to any Faner curve anywhere in the marched record.
**A refused comparison is a recorded result.**

This supersedes `rerun_2026-09-20/`. That folder is left byte-identical and
nothing in it is deleted; every number of it that is recomputed here is
recomputed from the same instruments, and where this record disagrees with it
the disagreement is stated with its cause rather than smoothed over.

Five changes were made, each derived from a published source or from the
source's own tables, each with its extrapolations and conventions declared.

| step | what changed | verdict |
|---|---|---|
| A | stage 1 posed on the geometry of the experiment | **the largest single improvement in the record**; the declared film interval now contains the measurement at every condition, where before it never did |
| B | a measured diffusivity law in place of one frozen value | the law is built and reported; marching it through a scalar-`D` solver is **convention-limited by an order of magnitude**, and the two families of curves demand opposite answers |
| C | a non-qualifying oil-solution arm on the residual floor | closes 24 to 37 per cent of the gap in the logarithm; **does not close it**, and carries a measured double-count exposure |
| D | the journal-family sentence withdrawn | replaced by an under-determination statement the sources themselves support, and step B turns out to strengthen the case for withdrawing it |
| E | identification kept out of the validation claim | supplement only; no coefficient of the marched record comes from a Faner curve |

Scripts, all run on this workstation under the pinned environment
(`PYTHONHASHSEED=1`, single-threaded BLAS/OMP/MKL/NumExpr,
`PYTHONDONTWRITEBYTECODE=1`, `PYTHONIOENCODING=utf-8`):

| script | wall | what |
|---|---|---|
| `build_diffusivity_law.py` | 98.3 s | step B: 34 retardation measurements, the regression, the declarations |
| `run_rerun_faner.py` | 649.7 s | steps A and C and the marched comparisons, 128 tail marches |
| `run_rerun_tail_residuals.py` | — | the falling-rate-window score |
| `build_uptake_benchmark.py` | 145.5 s | the Cardarelli benchmark at the converted coefficient |
| `run_identification_crosscheck.py` | — | step E |
| `oil_solution_arm.py` + `test_oil_solution_arm.py` | — | the arm and its 26 tests, all passing |

---

## A. Stage 1 on the geometry of the experiment

### A.1 The derivation, with the numbers

The 2026-09-20 march referred the frozen film coefficient to the particle's own
external area,

```
A/m_dry = 3/(R rho_dm,p) = 3/(0.885e-3 * 1159.65) = 2.923149664532036 m2 per kg dry meal
```

and over-predicted the measured constant-rate demand by 7.0 to 8.1 times.
Section 6 of that record measured why: the source's printed evaporation flux and
heat-transfer coefficient are referred to the **0.011 m² sample-holder
footprint**, over a layer 3 to 5 particle diameters deep, while the frozen
Coletto pair is a per-particle-area coefficient. The two are not the same
quantity.

The transfer area per unit dry mass **of the experiment** follows from two
printed numbers and nothing else:

```
A/m_dry = A_bed / m_dry ,      m_dry = W_sample / (1 + X_0)
```

| condition | `W_sample` (kg) | provenance | `X_0` (kg/kg) | `m_dry` (kg) | `A_bed/m_dry` (m²/kg dry) |
|---|---|---|---|---|---|
| 2008 soybean 136 °C | 0.02993 | thesis table 4.6 row 06, printed `W_harina`, `T_g` 137.1 °C | 0.735094 | 0.0172498 | **0.637689** |
| 2008 sunflower 136 °C | 0.03576 | thesis table 4.5 row 06, printed `W_harina`, `T_g` 138.2 °C | 0.735094 | 0.0206098 | **0.533726** |
| 2019 soybean 120 °C | 0.030 | `conditions.csv`, Sec. 2.2 p. 2, printed 30 g | 0.490628 | 0.0201257 | **0.546564** |
| 2019 sunflower 100 °C | 0.030 | same | 0.499279 | 0.0200096 | **0.549736** |

Identity of the tabulated runs with the plotted runs is **UNESTABLISHED** and is
in the refused register: the dataset record forbids merging the thesis and
journal experiments, and these are the nearest tabulated temperatures, not the
same runs.

The constant-rate balance is otherwise unchanged:
`dX/dt = -h (T_gas - T_b) (A/m_dry) / dh_vap`, with `h` the frozen Coletto pair
(131.3373155755338 W m⁻² K⁻¹) at all three ends of the declared interval,
`T_b = 68.71451253613918 °C` and `dh_vap = 334929.46668371576 J/kg`, all
reproduced bit for bit from 2026-09-20.

### A.2 The result

| condition | predicted / measured, `h×0.5 / ×1.0 / ×2.0` | 2026-09-20, same band |
|---|---|---|
| 2008 soybean | **0.765 / 1.530 / 3.061** | 3.507 / 7.015 / 14.029 |
| 2008 sunflower | **0.665 / 1.331 / 2.661** | 3.644 / 7.287 / 14.575 |
| 2019 soybean | **0.757 / 1.513 / 3.026** | 4.046 / 8.093 / 16.185 |
| 2019 sunflower | **0.676 / 1.351 / 2.703** | 3.593 / 7.186 / 14.372 |

**The declared film interval now contains the measured constant-rate demand at
every one of the four conditions. On 2026-09-20 it contained it at none.** The
constant-rate duration to the model's own `X_c = 0.1993532212219423` follows:

| condition | duration `h×0.5 / ×1.0 / ×2.0` (s) | measured `X_c` crossing (s) | duration / crossing at `h×1.0` |
|---|---|---|---|
| 2008 soybean | 63.68 / **31.84** / 15.92 | 54.930 | 0.580 |
| 2008 sunflower | 76.09 / **38.04** / 19.02 | 57.699 | 0.659 |
| 2019 soybean | 53.00 / **26.50** / 13.25 | 45.457 | 0.583 |
| 2019 sunflower | 88.94 / **44.47** / 22.24 | 65.085 | 0.683 |

Against the 2026-09-20 ratios of 0.126, 0.120, 0.109 and 0.129, this is the
measured crossing reached to within a factor 1.7 at the nominal film and
bracketed by the declared interval (the `h×0.5` ratios are 1.159, 1.319, 1.166
and 1.367).

### A.3 The disclosure that must travel with it

**This is not numerically the same construction as dividing the 2026-09-20
per-particle prediction by that record's area ratio, and the difference is
stated rather than hidden.** Section 6 of 2026-09-20 recovered a ratio of 7.82
to 9.71 and reported that dividing by it lands the corrected duration within
1.1 per cent (2008 soybean) and 5.5 per cent (2019 soybean) of the measured
crossing. The ratio marched here is 4.58 to 5.48, and the duration lands at
0.58 to 0.68 of the crossing. The cause is exact and arithmetic:

| condition | source particle area per kg dry (m²/kg) | frozen particle area per kg dry (m²/kg) | source / frozen | 2026-09-20 ratio | this record's ratio |
|---|---|---|---|---|---|
| 2008 soybean | 4.9868 | 2.9231 | **1.706** | 7.819 | 4.584 |
| 2008 sunflower | 5.1948 | 2.9231 | **1.773** | 9.711 | 5.477 |
| 2019 soybean | 4.3560 | 2.9231 | **1.490** | 7.970 | 5.348 |
| 2019 sunflower | 4.8266 | 2.9231 | **1.651** | 8.777 | 5.317 |

The 2026-09-20 ratio divides the **source's** particle area — built from the
source's own equivalent-sphere diameter and apparent particle density, applied
to the **total** (wet) sample mass — by the mesh area, and then applies that
ratio to the **frozen** particle area (radius 0.885 mm, `rho_dm,p` 1159.65, per
kg **dry** mass). Those two per-particle areas differ by 1.49 to 1.77, and that
factor is precisely the difference between the two records' stage-1 numbers.

The construction marched here uses only quantities of the experiment itself —
the printed sample mass, the printed initial loading and the printed holder
area — and never mixes a frozen geometry with a source geometry. It is the
weaker-looking number and it is the defensible one.

---

## B. A measured diffusivity law

### B.1 The conversion is an identity, and it is checked

Cardarelli fitted thesis equation (3.19), the Crank sphere solution of
`dq/dt = D grad^2 q`, which carries no storage capacity, so his coefficient is
an **apparent** one. The certified form carries the partition explicitly and at
the same numerical `D` runs slower by a retardation factor

```
r(T, a_h) = t63_certified(D) / t63_crank(3.19)(D)
```

measured here on the certified solver at each of Cardarelli's 34 tabulated
set-of-spheres coefficients. `r` does not depend on `D` — both time scales are
exactly inversely proportional to it — and that is **measured, not asserted**:
over two decades in `D` the ratio moves by `3.72e-9` relative. The
storage-consistent coefficient is then the algebraic identity

```
D*(T, a_h) = r(T, a_h) * D_table(T, a_h)
```

The isotherm behind `r` is the frozen PHY-043 native GAB of Cardarelli &
Crapiste (1996) with the PHY-048 residual-oil closure at the frozen
`w_o = 0.0195`, evaluated inside the certified solver; `r` is measured at the
saturated end state of the uptake march, and the loading it is reported against
is that state's own equilibrium loading
`X(a_h,T) = storage(a_h csat(T), T)/rho_dm,p`.

Measured across the 34 points: `r` = **1.190 to 18.983**, falling steeply with
temperature (8.6 to 19.0 at 50 °C, 1.19 to 1.42 at 95 °C) and rising with
activity. `D*` spans **1.207e-10 to 7.382e-9 m²/s**.

### B.2 The regression, and what it says

`ln D* = ln A - (Ea/Rg)(1/T) + beta ln(X/0.01)`, ordinary least squares, per
meal and pooled, with two alternatives fitted beside it. The carried form was
chosen by the lowest residual standard deviation on the **soybean** set, which
is the parameter set every Faner march of this record uses; the rule was fixed
before the numbers were seen.

| set | form | n | `Ea` (kJ/mol) | loading term | `sigma(ln D*)` | as a factor |
|---|---|---|---|---|---|---|
| sunflower | Arrhenius only | 17 | −29.06 | — | 0.7627 | ×2.144 |
| sunflower | **log loading** | 17 | **−0.455** | **0.9305** | **0.3229** | ×1.381 |
| sunflower | linear loading | 17 | −3.17 | 181.93 | 0.4275 | ×1.533 |
| soybean | Arrhenius only | 17 | −27.53 | — | 0.7695 | ×2.159 |
| soybean | **log loading (carried)** | 17 | **+0.306** | **0.9592** | **0.4104** | ×1.507 |
| soybean | linear loading | 17 | −4.37 | 168.64 | 0.5046 | ×1.656 |
| pooled | log loading | 34 | +0.54 | 0.9649 | 0.5637 | ×1.757 |

The carried soybean law is
`A = 3.513e-9 m²/s`, `Ea = 0.306 kJ/mol` (95 % interval ±16.09 kJ/mol),
`beta = 0.9592` (95 % interval ±0.3306), `sigma(ln D*) = 0.4104`, 14 degrees of
freedom.

**The physical reading, stated with its caveat.** Once the storage capacity is
carried explicitly the pore-gas diffusivity is indistinguishable from
temperature-independent and is close to proportional to loading. The apparent
activation energies Cardarelli reports (19.3 to 57.1 kJ/mol, his second table
3.7) are therefore largely the temperature dependence of the **sorption
capacity**, not of the pore transport. The caveat: the equilibrium loading at
fixed activity itself falls steeply with temperature, so the two regressors are
correlated (Pearson 0.590, variance inflation 1.535) and the activation energy
carries a ±16 kJ/mol interval. The correct statement is that `Ea` is
**indistinguishable from zero once the loading is carried**, not that the pore
transport has been measured to be athermal.

### B.3 Two extrapolations, both declared — and the loading one is the larger

*Temperature.* The measured range is 50 to 95 °C, i.e. `1/T` in
[2.7163e-3, 3.0945e-3] K⁻¹, a span of 3.7817e-4.

| target | `1/T` (K⁻¹) | distance beyond the hot end | as a fraction of the measured span |
|---|---|---|---|
| 100 °C | 2.6799e-3 | 3.640e-5 | **0.096** |
| 120 °C | 2.5436e-3 | 1.727e-4 | **0.457** |
| 136 °C | 2.4441e-3 | 2.722e-4 | **0.720** |

*Loading, which is worse.* Every tabulated coefficient is measured over a
sorption isotherm at activity 0.10 to 0.81, so the particle holds only its
sorbed inventory: the measured loadings run **7.677e-4 to 1.5616e-2 kg/kg**.
Faner's falling-rate window starts at the critical loading, where the pore still
holds condensed n-hexane, and runs from 0.19935 down to 0.1006 (2008) or
0.007415 (2019). **The two ranges do not overlap above 1.56e-2 at all.**
Carrying the law into the 2008 window extrapolates its loading regressor by a
factor 6.4 to 12.8. That is a larger extrapolation than the temperature one and
it is labelled on every number it touches.

Four diffusivity cases are therefore marched, and a declared-hold case refuses
to evaluate the law outside its measured loading envelope:

| case | what |
|---|---|
| `phy019_frozen_4.0e-10` | the 2026-09-20 comparison trace |
| `law_log_loading` | the carried law, extrapolated in `T` **and** in `X` |
| `law_log_loading_held` | the same law with its loading regressor **held** at 1.5616e-2 kg/kg above that maximum — a declared refusal to evaluate outside the measured envelope, not a clamp on a solver refusal |
| `law_arrhenius_only` | the temperature-only fit; no loading extrapolation at all, at the cost of a residual factor 2.16 |

The two ends of the 2026-09-20 Cardarelli band (7.119e-11 and 8.169e-10) are
kept so that the earlier bracket is not lost.

### B.4 The construction limit that governs every loading-dependent number

The certified dry shell is posed with **one scalar `D_eff` per step**. A
loading-dependent law has no correct place inside its flux divergence, and this
record does not pretend otherwise: the coefficient is evaluated once per stride,
at the volume mean and, separately, at the surface, and the two are reported as
a **bracket**.

| condition | case | mean-evaluated time to the last measured loading (s) | surface-evaluated (s) |
|---|---|---|---|
| 2008 soybean | `law_log_loading` | 3.77 | 298.88 |
| 2008 soybean | `law_log_loading_held` | 29.91 | 306.06 |
| 2008 sunflower | `law_log_loading` | 3.38 | 264.72 |
| 2008 sunflower | `law_log_loading_held` | 27.51 | 271.91 |
| 2019 soybean | `law_log_loading` | 185.76 | not reached in 2400 s |
| 2019 sunflower | `law_log_loading` | 400.46 | not reached in 2400 s |

**The bracket spans nearly two orders of magnitude**, because `D` varies by
about a decade across the particle. No loading-dependent number in this record
is quoted without it, the figure draws the law as a band rather than a line, and
"a space-resolved loading-dependent diffusivity `D(X(r))` inside the certified
flux divergence" is now a named element of the refused register.

### B.5 The falling-rate window, before and after

Every measured sample at or after the **measured** crossing of the source's
declared `X_c = 0.20`; the window is fixed by the measurement. Sensible
temperature case, nominal film. RMS as a fraction of the window span.

| condition | 2026-09-20, frozen `D` | **2026-09-21, frozen `D`** | **2026-09-21, measured law** | 2026-09-21, held law | 2026-09-21, Arrhenius only |
|---|---|---|---|---|---|
| 2008 soybean 136 °C | 0.1629 (8/15 in band) | **0.3221** (1/15) | 1.2576 (0/15) | 0.8131 (0/15) | 0.3775 (1/15) |
| 2008 sunflower 136 °C | 0.1205 (13/15) | **0.2534** (2/15) | 1.1421 (0/15) | 0.7280 (0/15) | 0.2892 (2/15) |
| 2019 soybean 120 °C | 0.6420 | **0.7124** | **0.0883** | 0.3043 | 0.6983 |
| 2019 sunflower 100 °C | 0.6260 | **0.6806** | **0.1216** | 0.3691 | 0.6392 |

Mean signed residuals at the frozen `D` move from +0.00845, +0.00140, +0.10580,
+0.10050 kg/kg (2026-09-20) to +0.02445, +0.02094, +0.11937, +0.11199 kg/kg.

**Three findings, reported as they came out.**

1. **Step A makes the thesis comparison worse at the frozen diffusivity.** The
   2026-09-20 march reached the crossing at 6.9 s and then took a slow tail; the
   fast stage 1 was compensating for a slow stage 2. With stage 1 posed
   honestly, the compensation is gone and the model is wetter throughout the
   window: 0.1629 becomes 0.3221 and the in-band count falls from 8 of 15 to 1
   of 15. **The 2026-09-20 thesis agreement was in part an artefact of the wrong
   transfer area.** That is worth more than the better-looking number it
   replaces.
2. **The measured diffusivity law transforms the journal comparison.** 0.642 and
   0.626 become **0.088 and 0.122** of the window span, and the systematic
   +0.10 kg/kg wet bias becomes −0.005 and −0.002 kg/kg. The journal traces do
   not need a different model; they need a faster shell.
3. **The same law ruins the thesis comparison**, to 1.26 and 1.14 of the span.
   The two published families demand diffusivities an order of magnitude apart,
   and which one the model follows is decided by the transport closure, not by
   the model's structure. This is reported as it is; no case was selected to
   make a figure look better, and all four cases are in
   `rerun_tail_residuals.csv`.

### B.6 Numerical resolution

Bounding case at 409.15 K from `X_c`, two diffusivity cases:

| case | cells | stride cap (s) | steps | t to 0.1005939 kg/kg | t to 0.0735094 kg/kg |
|---|---|---|---|---|---|
| frozen | 60 | 5.00 | 480 | 365.185 | 660.677 |
| frozen | 120 | 5.00 | 480 | 364.951 | 660.474 |
| frozen | 60 | 1.25 | 1920 | 365.012 | 660.488 |
| held law | 60 | 5.00 | 480 | 29.905 | 53.944 |
| held law | 120 | 5.00 | 480 | 29.887 | 53.928 |
| held law | 60 | 1.25 | 1920 | 29.709 | 53.712 |

Spread 0.064 % (frozen) and 0.66 % (held law). Neither result is a
discretization artefact; the loading-dependent case is convention limited, as
B.4 says, and not mesh limited.

---

## C. The oil-solution arm on the residual floor

### C.1 What it is

`oil_solution_arm.py` in this folder — analysis-side, **not** a promoted `src`
module, every object typed `physically_qualifying = False`. n-Hexane dissolved
in the residual triglyceride oil, in equilibrium with the pore gas, by Raoult's
law with an activity coefficient:

```
x_h = a_h / gamma_h ,     W_oil = w_o * (M_h / M_oil) * x_h / (1 - x_h)
```

`docs/LITERATURE_MANIFEST.md` was searched for a hexane / vegetable-oil
vapour-liquid source and carries **none**: the nearest entries are Cardarelli's
meal sorption isotherms, Khudaida's aqueous C6 systems and Gmehling's azeotropic
data, none of which is a hexane-miscella vapour pressure. The **ideal solution**
`gamma_h = 1` is therefore the central case, declared as such, with `gamma_h`
bracketed over [1.0, 1.5]; a source-pinned coefficient is **CITETODO for the
owner**. The oil molar mass is computed from atomic masses for the
triacylglycerols that bracket these oils — trilinolein C57H98O6 (879.405 g/mol)
and triolein C57H104O6 (885.453 g/mol), with tripalmitin C51H98O6 (807.339) as
the saturated end member — and the central value carried is 882.429 g/mol. The
fatty-acid composition that would fix it exactly is also CITETODO.

Refusals, not clamps: a non-finite or negative activity, a mole fraction that
reaches unity (`the liquid is no longer an oil-rich solution`), a non-positive
activity coefficient or molar mass, and an oil fraction outside `[0, 1)`.

### C.2 Continuity and the mass ledger, measured

* The arm at `w_o = 0` is **exactly `0.0`**, the same double-precision number as
  the arm-free model, and the corrected floor at zero oil is bit-identical to
  the sorbate floor.
* The arm is **exactly linear** in the oil fraction: the slope spread over
  `w_o` from 1e-1 down to 1e-12 is exactly `0.0`.
* The three-compartment ledger (sorbate, pore gas, oil solution) closes with a
  residual of **exactly `0.0`**, and removing the arm recovers the arm-free
  total exactly.
* 26 tests in `test_oil_solution_arm.py`, run from this folder, all passing.

### C.3 The corrected floor

| condition | `a_h` at `T_gas` | `w_o` | sorbate floor (kg/kg) | arm, ideal (kg/kg) | arm at `gamma = 1.5` | corrected floor, ideal | measured final | measured / sorbate | measured / corrected |
|---|---|---|---|---|---|---|---|---|---|
| 2019 soybean 120 °C | 0.2534628 | 0.0195 | 1.2169e-3 | 6.4655e-4 | 3.8721e-4 | **1.8634e-3** | 7.415e-3 | 6.09 | **3.98** (4.62 at `gamma = 1.5`) |
| 2019 sunflower 100 °C | 0.4113983 | 0.0180 | 2.6546e-3 | 1.2286e-3 | 6.6431e-4 | **3.8832e-3** | 7.415e-3 | 2.79 | **1.91** (2.23 at `gamma = 1.5`) |
| 2008 soybean 136 °C | 0.1786177 | 0.0195 | 7.6072e-4 | 4.1411e-4 | 2.5742e-4 | 1.1748e-3 | 0.1006 (where the trace stops) | 132.2 | 85.6 |
| 2008 sunflower 136 °C | 0.1786177 | 0.0180 | 7.6072e-4 | 3.8226e-4 | 2.3761e-4 | 1.1430e-3 | 0.1042 (where the trace stops) | 136.9 | 91.2 |

The two 2008 rows measure only where the experiment stopped, not a floor. On the
two 2019 runs that do reach a plateau, the arm **closes 23.6 and 37.0 per cent of the
logarithmic gap and does not close it.** The model still predicts a drier
particle than the experiment reaches, by a factor 2 to 4 instead of 3 to 6.

### C.4 The double-count exposure, measured rather than hidden

The frozen PHY-043 native GAB is a **total-meal** isotherm measured on a meal
that already contained its residual oil; PHY-048's own docstring says the native
anchor already contains its native oil. Adding this arm at the native oil
fraction therefore risks double counting whatever part of the oil solution the
measurement already contains. The size of that exposure is measured:

| temperature | `arm / native GAB` at `a_h` 0.179 / 0.253 / 0.411 |
|---|---|
| 50 °C | 0.174 / 0.196 / 0.249 |
| 65 °C | 0.234 / 0.258 / 0.315 |
| 80 °C | 0.318 / 0.342 / 0.406 |
| 95 °C | 0.430 / 0.456 / 0.527 |
| 120 °C | 0.693 / 0.725 / 0.818 |
| 136 °C | 0.923 / 0.961 / 1.076 |

At the isotherm's own measured conditions the arm is **17 to 53 per cent** of
the measured native retention, so a material part of it may already be inside
the GAB fit. At desolventizer temperature the native branch has collapsed while
the arm has not, and the arm is 69 to 108 per cent of it — which is precisely
why a floor that survives the isotherm's collapse is the physically right
candidate, and precisely why it cannot simply be added.

**The corrected floor is therefore reported as the UPPER bound of a bracket
whose lower bound is the arm-free floor.** Both bounds sit below the measured
final loading, so the finding of section 4.3 of the 2026-09-20 record stands:
*the model predicts a drier particle than the experiment reaches, and the oil
arm narrows but does not remove the discrepancy.*

---

## D. The journal-family sentence, withdrawn

The 2026-09-20 record concluded: *"The two sources cannot both be right, and the
model sides with the thesis."* **That is stronger than the evidence, and it is
withdrawn.**

Two reasons, the first from the sources and the second from this record.

*The sources.* `paper/analysis/datasets/faner2019/DATASET_RECORD.md` and
`conditions.csv` state that gas flow rate, meal amount and initial solvent
loading were **not** held constant across the journal's temperature series, and
the Figure 2 conditions are not tabulated per run. The comparison to the journal
curves is therefore **under-determined by the source's own statement**: there is
no set of per-run conditions against which the model could be scored, so a
disagreement cannot be attributed to either party.

*This record.* The claim that the model "sides with the thesis" was a property
of the transport closure, not of the model. With the frozen scalar diffusivity
the march follows the thesis (0.32 and 0.25 of the window span against 0.71 and
0.68); with Cardarelli's own measured diffusivity converted to this model's
definition, **the same march follows the journal** (0.088 and 0.122 against 1.26
and 1.14). A conclusion that reverses when the diffusivity closure is replaced
by a measured one is not a conclusion about which published family is right.

What stays, as a finding and not a validation:

* at the thesis's own tabulated conditions and the frozen diffusivity the
  formulation reproduces the thesis falling-rate leg to 25 to 32 per cent of the
  window span;
* the regime map says what would produce the journal trajectory: the journal
  traces identify a mass Biot number of **5 to 10** against the thesis's **269
  and 281** under the same closed-form law, a gap of about thirtyfold. At a
  fixed external coefficient `Bi_m` falls as `D` rises, so the journal
  trajectory is reached by a shell diffusivity about thirty times the frozen
  value — and the measured Cardarelli law, evaluated at the loadings of the
  falling-rate window, supplies 1.5e-8 m²/s at `X = 0.05`, which is 37 times
  4.0e-10. The alternative, a particle smaller by the same factor in `R²/D`,
  would need a radius 6 times smaller than either source prints and is not
  supported;
* the two published families are separated by a quantity this work can name and
  measure, and are not adjudicated here.

---

## E. The identification, kept out of the validation claim

Step B is an identification **on Cardarelli** and then a prediction **on
Faner**: the diffusivity law is regressed on Cardarelli's tabulated coefficients
and never touched afterwards. No film coefficient, no critical loading and no
diffusivity in the marched record is fitted to the curve it is scored on.

As a supplementary cross-check, `run_identification_crosscheck.py` identifies
**one** scalar diffusivity on the thesis **sunflower** 136 °C trace by
minimising its falling-rate-window residual, and then **predicts** the thesis
**soybean** 136 °C trace with it, nothing else changed. The reported direction
was fixed before the numbers were seen; the reverse direction is run and
reported beside it. The identified residual is a fit and is labelled one; the
held-out residual is the number that means something. **These numbers live in
the supplement and never in the validation section.**

---

## F. The refused register — 51 elements, ten classes

`rerun_refused_register.csv`.

| class | count |
|---|---|
| data-side, not invertible (the twelve normalized traces) | 12 |
| not data (the authors' own model lines) | 12 |
| model-side declared-domain boundary | 8 |
| instrument refusal (trace leaves the band) | 6 |
| data-side, identity unestablished (the four tabulated sample masses) | 4 |
| model-side root cause, routed by the extensions | 3 |
| data-side, metric not reached | 2 |
| data-side, no reading band | 2 |
| data-side, under-determined (the journal's per-run conditions) | 1 |
| literature not held (the hexane / vegetable-oil activity coefficient) | 1 |

The 2026-09-20 register had 43. The eight new elements are: the four
sample-mass identity refusals of step A; the journal's per-run conditions
(step D); the loading-range non-overlap of the diffusivity law (step B.3); the
space-resolved `D(X(r))` the certified divergence cannot carry (step B.4); and
the missing hexane / vegetable-oil activity coefficient (step C).

Still standing, unchanged: the constant-rate **rate** is posed from a
single-sphere energy balance, not from the certified moving-front march; the
free-liquid regime above `X_c` remains unmarched; the critical loading puts the
particle on the GAB branch above `a_h = 1`, inside the solver's own ceiling but
outside the isotherm's measured activity domain.

---

## G. The Cardarelli uptake benchmark, at the converted coefficient

`build_uptake_benchmark.py` re-marches every series item 13 posed, at
`D* = r D_table`, and scores it on the same points with the same metrics. Thirteen
independent series, 100 measured points, 50 °C.

| coefficient | rms relative span | 63 % time ratio span | points inside the declared band | inside the two-reading band |
|---|---|---|---|---|
| the source's tabulated apparent value | 0.445 – 1.193 | 6.13 – 24.43 | 3 of 100 | 0 |
| the frozen 4.0e-10 | 0.096 – 0.635 | 0.58 – 18.29 | 11 of 100 | 0 |
| **the converted `D*`** | **0.041 – 0.297** | **0.57 – 2.19** | **35 of 100** | **2** |

The conversion is the source's own fitted coefficient carried into this model's
definition by a measured identity, so a curve drawn at it is an
**IDENTIFICATION on Cardarelli**, never a prediction; it is labelled that way in
every output and in the manuscript. It is nevertheless the sharpest statement in
the record that the 2026-09-19 disagreement was a definition mismatch and not a
failure of the transport form: the time-constant ratio band collapses from a
factor 32 to a factor 3.9.

---

## H. Claim boundary

`physically_qualifying: false. plant_predictive: false.`

Both domain extensions remain typed non-qualifying authorities. The oil-solution
arm is analysis-side, non-qualifying and carries a CITETODO for its activity
coefficient. The marched trajectory is a two-stage declared construction, not a
certified moving-front solution: the constant-rate leg is an energy balance on a
transfer area read from the source's tables, the falling-rate leg starts on a
GAB branch above `a_h = 1`, the non-isothermal coupling is a declared
construction, and a loading-dependent diffusivity is evaluated once per stride
and reported as a bracket. The film interval is a declared engineering interval,
not a measurement of the rig, and step A changes only the area it is referred
to. Every inversion of measured data is reported as what the rig requires and is
never a model coefficient. Plant data is nowhere in this record.

---

## I. Addendum (2026-09-21, same day): both evaluation edges, and the physics
that ranks them

B.5 above and the paper's `fig_faner_march` / `tab:fanermarch` originally
quoted only the `law_log_loading` case — the volume-mean evaluation of the
measured law — against the frozen coefficient, even though B.4 already
established that the law is marched and read at **both** edges (mean and
surface) every stride, and reported the two as a bracket. Quoting one edge in
the window-RMS table and the figure annotation was a case selection by
omission. Both edges of `rerun_tail_residuals.csv` (sensible, nominal film)
are below, alongside B.5's frozen and volume-mean columns:

| condition | surface (`law_log_loading_surface`) | volume mean (`law_log_loading`) |
|---|---|---|
| 2008 soybean 136 °C | **0.1413** (11/15), mean +0.0058 | 1.2576 (0/15), mean −0.0992 |
| 2008 sunflower 136 °C | **0.1102** (13/15), mean +0.0009 | 1.1421 (0/15), mean −0.0994 |
| 2019 soybean 120 °C | 0.5739, no band, mean +0.0938 | **0.0883**, no band, mean −0.0046 |
| 2019 sunflower 100 °C | 0.5118, no band, mean +0.0804 | **0.1216**, no band, mean −0.0023 |

For a diffusivity that increases with loading, desorption from a sphere is
governed by the region of lowest diffusivity, which is the surface, and the
surface-evaluated coefficient is the scalar the space-resolved solution
approaches (Crank, *The Mathematics of Diffusion*, 1956); the volume-mean
evaluation over-predicts the flux by the ratio of the two coefficients, about
a decade here. The surface edge is therefore the physically indicated one: it
is the closer match on the thesis traces, and the volume-mean edge, the
closer match on the journal traces, is the unphysical end of the same
bracket. Both are now reported everywhere the law is quoted — the figure
panel notes, `tab:fanermarch`, `RESULT_2026-09-21.md`, and the supplement's
addendum at `sup:bothedges` — and neither is selected.

---

## J. Addendum (2026-09-21): the floor at desolventizer conditions

The owner, reading the draft, read the paragraph's `1.217e-3` /
`2.655e-3` kg/kg floor numbers as the model's residual-hexane floor and
concluded the model cannot represent the industry's KPI (residual hexane in
meal, hundreds of ppm and below). Those numbers are the dry-shell equilibrium
floor at the **Faner experiments' own gas** — pure superheated hexane vapour,
`a_h = P/Psat(T) =` 0.253 at 120 °C and 0.411 at 100 °C — read against the
same experiment's own measured plateau, `7.415e-3` kg/kg. They are not a
desolventizer reading, because a desolventizer's last trays run mostly steam
with a small hexane mole fraction, not pure hexane vapour.

`floor_at_desolventizer_conditions.py` reruns the same isotherm chain over a
generic grid, T = 95–120 °C and hexane mole fraction y = 0.001–1.0 (`a_h =
min(1, y·P/Psat(T))`), writing `floor_at_desolventizer_conditions.csv`. The
mole fractions are generic desolventizer-scale values, not a plant reading —
no plant gas composition is asserted anywhere in this record.

**Two constructions, reconciled.** The brief's own worked table (§ "Facts
measured this session") calls `sorption.retained_hexane(a_h, T, w_o,
GabParams(), OilIsotherm())` directly: the GAB+oil composite retention
isotherm alone, sorbed hexane per kg dry meal at equilibrium with a gas at
activity `a_h`. At the anchor condition (120 °C, pure vapour, `a_h =
0.25346`) this returns **892.10 ppm**. But the paper's printed `1.217e-3`
kg/kg (**1216.88 ppm**) comes from the rerun's `rerun_floor.csv` column
`sorbate_floor_kg_kg`, computed as `K.ppm_from_cg(a_h·cg_saturation(T), T,
p)`, which calls the certified dry-shell solver's own storage functional
`dry_shell.storage(c_g, T, p) = eps_g·c_g + rho_dm_p·W_h(a_h, T, w_o)` — the
pore **gas-phase** hexane inventory PLUS the same sorbed term, both divided
by `rho_dm_p` and reported in ppm of dry meal. Despite its name,
`sorbate_floor_kg_kg` is a TOTAL (gas + sorbed) particle inventory, not a
sorbed-only isotherm floor. The two reconcile exactly:

`892.10 ppm (sorbed) + 324.79 ppm (pore gas-phase term, eps_g·c_g/rho_dm_p) =
1216.88 ppm (total)`,

asserted to better than `1e-9` ppm in the script. **The paper's validation
number is the total (gas + sorbed) construction.** Both constructions are
carried in the output CSV at every grid point, for both this reason and
because the two do not read the same physical question: the sorbed-only
number is the isotherm's own equilibrium retention, the total is what a
particle's whole pore-plus-sorbed inventory holds at that gas condition.

**The desolventizer reading (sorbed construction, matching the brief's
worked table).** At y = 0.003 to 0.03 (100–110 °C), the same isotherm floors
at 3.5 to 55.8 mg/kg (rounds to the paper's "4 to 56"); at y = 0.1 over the
same temperature band, 118.9 to 188.2 mg/kg ("120 to 190"). These sit at or
below the low end of the industry's typical residual-hexane specification
(hundreds of ppm and below), unlike the experiment's own pure-vapour floor,
because the desolventizer's last trays never see pure hexane vapour. The
residual-hexane KPI is therefore set by the kinetics of the last trays and
the gas composition there, not by the equilibrium of this isotherm at the
Faner experiments' gas.

Claim boundary unchanged: `physically_qualifying: false. plant_predictive:
false.` No hexane mole fraction above is a plant reading; no law moved; no
pinned source touched.
