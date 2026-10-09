# Item 19 — fatal limit 3: the two typed domain extensions, against the evidence they can have

**Executes** limit 3 of `paper/01_particle_jfpe/LIMITS_AUDIT_2026-09-22.md`
(weight **fatal**, bearing on "every march number"): *the single-component gas
limit and the zero-water continuation that make the march posable lie outside
the qualified domain; the pure-hexane boundary is the moisture boundary in
another coordinate.*

**Outcome, in one paragraph.** Neither extension becomes qualified, and neither
is claimed to be. What changes is that the fatal limit is now **bounded and
measured instead of declared**, and it is bounded in two different ways for the
two extensions. The **dry-limit moisture continuation** is confirmed to be a
genuine extrapolation: of twenty-eight entries from five published sources, the
lowest activity anyone measured on soybean or sunflower material is
**0.07**, half a decade above the continuation's upper edge at
$a_w = 0.046973$, so **no held source measures anywhere inside the
continuation's domain** — but the march never integrates over that domain
either, because Faner's charge is water free and the continuation's only role is
to admit the single endpoint $W_w = 0$, where the retained-water flux is
identically zero by the absence of the component. The **single-component
pore-gas limit** turns out to exercise a strict **subset** of the binary path's
closures: at $y_h = 1$ the cross second virial, its $k_{wh}$ interval, its
below-direct-evidence flag and the entire IAPWS-95 water branch are annihilated
exactly, the vertex state being **bitwise invariant over the whole declared
$k_{wh}$ interval**, and the Maxwell-Stefan flux and its unidentified mobility
interval leave the problem altogether. The one number either module contributes
to the marched arithmetic of item 5 is an **inversion of a frozen pure-component
correlation** — the saturation temperature 68.71451 °C — and that number is
scored against measurement to $-0.137$ and $+0.093$ °C.

`physically_qualifying: false. plant_predictive: false.` Nothing is fitted,
nothing is clamped, no frozen law is moved, no module under `src/`, `tests/`,
`scripts/`, `tools/` or `docs/` is touched, and the manuscript is not edited.
Plant data appear nowhere; every placed number comes from a published source
cited by its bibliography key and its printed locator.

**Method.** `run_item19_dry_limit.py` and `run_item19_single_component.py`
(0.03 s and 0.3 s), then `make_item19_figure.py`. The published sources are not
re-transcribed here: they are imported from
`../item_18_water_branch_placement/sources.py`, which carries each number with
its locator, and from item 11's digitized sunflower JSON, so this item cannot
drift away from item 18. Ruff clean on all three files.

---

## Part A — the dry-limit moisture continuation

### A.1 The placement table

The continued branch is $W_w(a_w) = A_1 / (1 + A_2 \ln(1/a_w))$ with the frozen
$A_1 = 0.880$, $A_2 = 12.184$; its qualified band is $W_w \in [0.023,
0.2357]$, that is $a_w \in [0.046973, 0.799]$. The inverse used here round-trips
through `dry_limit_moisture_branch.dry_limit_water_activity` to machine
precision at every row, so it is the same object and not a second law.

Per source, at the **lowest activity that source actually measured**:

| source | material, branch | T (K) | reading | lowest measured $a_w$ | measured $W_w$ | continuation $W_w$ | ratio |
|---|---|---:|---|---:|---:|---:|---:|
| `pixton1975moisture` | soya meal, desorption | 288.15–308.15 | printed points | 0.183–0.224 | 0.06270 | 0.04057–0.04577 | 0.647–0.730 |
| `pixton1975moisture` | soya meal, adsorption | 288.15–308.15 | printed points | 0.299–0.363 | 0.07411 | 0.05602–0.06593 | 0.756–0.890 |
| `aviara2004sorption` | whole soya bean, desorption | 313–343 | printed points | **0.070** | 0.0528–0.0852 | 0.02635 | 0.309–0.499 |
| `zeymer2023moisture` | soybean grain, both branches | 283.15–323.15 | own fit at its stated 0.11 | 0.110 | 0.03487–0.05071 | 0.03155 | 0.622–0.905 |
| `cassini2006water` | texturized soy protein, adsorption | 283.15–313.15 | own GAB at its stated 0.07 | **0.070** | 0.02173–0.03498 | 0.02635 | 0.753–**1.213** |
| `cardarelli1998thesis` | sunflower meal, adsorption | 322.87–367.86 | digitized points | 0.269–0.329 | 0.05000 | 0.05184–0.06049 | 1.037–1.210 |

Full rows in `item19_low_activity.csv`; every measured point of every source in
`item19_low_activity_points.csv`.

### A.2 What the table says

1. **The continuation is an extrapolation, not an interpolation, and the
   interval is now named.** Zero of twenty-eight entries have a lowest measured
   activity inside the continuation's domain. The lowest activity measured
   anywhere on this material class is **0.070** (Aviara's printed Table 2 at
   four temperatures; the lower limit Cassini states for its own range), which
   is **1.49 times** the continuation's upper edge $a_w = 0.046973$. The whole
   continued branch, $a_w \in [0,\,0.046973)$ and $W_w \in [0,\,0.023)$, is
   unmeasured by every source held.

2. **In the loading coordinate the gap is narrower than that, and two entries
   cross it.** Cassini's own GAB read at the lowest activity Cassini measured
   gives **0.02173** and **0.02295** kg/kg dry at 10 and 30 °C — *below* the
   model's qualified anchor 0.023. So two of twenty-eight held entries place a
   measured-activity loading inside the continuation's loading domain, and the
   lowest measured loading held anywhere is 0.0217, a factor 1.06 below the
   anchor rather than the order of magnitude a reader might assume.

3. **Where measurement does exist, the frozen branch lies inside the spread of
   it.** Over the twenty-four printed or digitized points below $a_w = 0.35$
   (activity interval 0.070 to 0.337; Aviara 15, Pixton 5, Cardarelli 4) the
   ratio of the branch to the measurement runs **0.309 to 1.210** and straddles
   unity: the branch reads about half of the whole-bean loading, two-thirds to
   nine-tenths of the soya meal, and about 1.04 to 1.21 times the sunflower
   meal; adding the two sources read through their own fits at their own lowest
   measured activity extends the same straddle to 0.75 and 1.21 on the soy
   protein. The branch is therefore not systematically outside the low-activity
   evidence; it sits inside a spread that is itself a factor of about four wide
   because the five sources are five different materials: soya meal, whole
   soya bean, soybean grain, texturized soy protein and sunflower meal.

4. **Containment below the anchor is a weak constraint, and is reported as
   such.** Continued below $a_w = 0.046973$, the continuation lies inside the
   band formed by the three published sources that carry their own fitted form
   (Aviara's modified Halsey, Zeymer's modified Halsey, Cassini's GAB) at all 33
   grid points down to $a_w = 4.7\times10^{-6}$. But every value in that band is
   a published fit read **below the range its source measured**, and the band
   widens from a factor of **5.13** at the anchor to a factor of
   **2.6 × 10⁴** four decades below it. Containment inside a band that wide
   constrains nothing; it is recorded so that a referee sees it was checked and
   found weak, not as evidence.

5. **The march does not traverse the continued branch.** In the item-5
   marches the continuation appears through one call,
   `dlmb.qualify_dry_limit()`, whose result is recorded and not used in any
   trajectory. Faner's charge is water free, so the continuation's operative
   content is the single endpoint $W_w = 0$, where the module types the
   retained-water potential as $-\infty$, refuses to let it enter a force
   difference, and states the retained-water flux as identically zero **by the
   absence of the component**. The module's own measurement makes the endpoint
   exact rather than approximate: below $W_w = 9.6919\times10^{-5}$ kg/kg dry
   the frozen exponential underflows to exactly `0.0`, so a charge at
   $10^{-5}$ or $10^{-6}$ kg/kg and a charge at exactly zero are the **same
   double-precision number**.

### A.3 Verdict, extension 1

**Extrapolation, unmeasured over its whole domain, but exercised at one
endpoint only, where its value is exact by underflow and its flux is zero by
absence.** The referee's reading of limit 3 — that a fitted isotherm is being
evaluated far outside its data and the drying numbers depend on it — is not
what the evidence shows for this extension: no march number of item 5 depends on
a value of the continued branch at any interior point. The limit that survives
is narrower and should be stated in its narrow form.

---

## Part B — the single-component pore-gas limit

### B.1 The state limit, over the whole falling-rate window

The build note qualified this limit at three gas temperatures. The item-5 traces
carry particle temperatures from the frozen saturation temperature
**341.86451 K (68.71451 °C)** up to the gas temperature, **409.15 K** for both
2008 thesis traces (`rerun_2026-09-22_rig/rerun_traces.csv`). The continuity is
therefore re-measured here on a 33-point grid spanning that whole window at
101 325 Pa, with carrier mole fractions $10^{-2}$ down to $10^{-6}$, and
extended from the six gas-state quantities of the build note to **thirteen**,
including the storage and caloric quantities the march actually transports.

| transported quantity | difference | worst at $\varepsilon = 10^{-6}$ | order over the window |
|---|---|---:|---|
| molar density | relative | 8.66e-08 | 0.9991–0.9992 |
| hexane mass density | relative | 1.09e-06 | 0.9998–0.9999 |
| hexane fugacity | relative | 1.00e-06 | 0.9999–1.0000 |
| hexane partial pressure | relative | 1.00e-06 | 1.0000 |
| compressibility | absolute | 8.25e-08 | 0.9992–0.9993 |
| $B_{\mathrm{mix}}$ | absolute, m³/mol | 2.32e-09 | 0.9992–0.9993 |
| residual molar enthalpy | absolute, J/mol | 8.68e-04 | 0.9991–0.9992 |
| pore hexane concentration | relative | 1.09e-06 | 0.9998–0.9999 |
| **retained hexane loading** | relative | **3.94e-06** | 0.9964–0.9999 |
| **total hexane concentration** | relative | **3.88e-06** | 0.9964–0.9999 |
| $\ln\varphi_h$ | absolute | 5.88e-14 | **2.0000** |
| partial residual enthalpy, hexane | absolute, J/mol | 7.26e-10 | **2.0000** |
| water fugacity | absolute, Pa | 1.05e-01 → 0 | 0.9999 |

Every transported quantity converges to the vertex value at first order, two of
them at second order because the vertex is a stationary point of the hexane
fugacity coefficient in composition. The water fugacity goes to zero at first
order, as the component vanishing requires. The exchange force
$\varphi = \ln(f_w/f_h)$ is the one object that diverges, logarithmically, with
measured decade increments **2.30259 to 2.31112** against $\ln 10 = 2.302585$.
`item19_carrier_continuity.csv`, 165 rows.

### B.2 The pore-gas flux does **not** converge, and that is the right result

The brief asks for the continuity of the pore-gas flux as well. It is measured
and the answer is negative, in a way that supports the module rather than
contradicting it.

`coupled_pore.independent_maxwell_stefan_flux` is
$J_w = -L\,(\varphi_R - \varphi_L)/\delta$ with $L$ a caller-supplied,
composition-**independent** positive constant
(`coupled_pore.MobilityInterval`, no endpoint value approved; limit 16). The
face difference of $\varphi = \ln(f_w/f_h)$ depends on the **ratio** of the two
carrier fractions and not on their size. The scaled flux $J_w\delta/L$ therefore
needs no mobility value at all, and it has **no path-independent limit**:

* two faces approaching the vertex at a **fixed ratio of 2**: the scaled flux
  tends to $-\ln 2 = -0.693147$, worst deviation **9.58e-03** at
  $\varepsilon_L = 10^{-2}$, and it does **not** approach zero;
* two faces at a **fixed absolute difference** of $10^{-8}$: the scaled flux
  runs from $-1.01\times10^{-6}$ to $-2.3979$ as $\varepsilon_L$ falls from
  $10^{-2}$ to $10^{-9}$, and grows without bound.

The reduction's value is exactly `0.0`. So the reduction is **not** the limit of
the binary flux formula, which is exactly what the module's docstring says and
what the binary topology's refusal encodes: the zero holds because the second
component is absent, not because a limit was taken. `item19_flux_paths.csv`,
Figure panel (c). **The manuscript must not write that the pore-gas flux is
continuous into the limit.** What is continuous is the state.

### B.3 Which closures the limit exercises that the binary path does not

Measured at the simplex vertex at five temperatures spanning the window
(`item19_closure_inventory.csv`, and the `vertex_identity` block of
`item19_single_component.json`):

* $B_{\mathrm{mix}} = B_{hh}$ **bitwise** at all five, the cross-virial
  contribution to $B_{\mathrm{mix}}$ being exactly `0.0`;
* $\ln\varphi_h$ equals the pure closed form $PB_{hh}/RT$ with an absolute
  difference of exactly **0.0**;
* the vertex state is **bitwise invariant over the whole declared $k_{wh}$
  interval** $[0.477, 0.511]$, while the binary state at the coupled-pore
  carrier floor $y_h = 0.857$ moves by $2.3\times10^{-5}$ to
  $4.0\times10^{-5}$ in $\ln\varphi_h$ over the same interval;
* the authority's own `below_cross_direct_evidence` flag reads **False** at the
  vertex at every temperature, including 341.86 K where the cross correlation
  is below its 363.2 K direct-evidence temperature, because the authority sets
  that flag only when both mole fractions are positive
  (`binary_gas.py:240`). The evidence question the flag exists to raise does
  not arise at the vertex, by the binary authority's own construction.

So the inventory is:

| closure | binary path | the limit |
|---|---|---|
| pure n-hexane second virial $B_{hh}$ | yes | yes, at $y_h = 1$ |
| cross second virial $B_{wh}$, $k_{wh}\in[0.477,0.511]$ | yes, enters $B_{\mathrm{mix}}$ and both $\varphi$ | **removed**: multiplied by $y_w = 0$, contribution exactly zero |
| IAPWS-95 water branch | yes | **removed**, same way |
| PHY-031 modified-Luikov retained water | yes — and it is the clause that refuses the approach at $y_h = 0.857$ | **removed**; the endpoint is the other extension's |
| Maxwell-Stefan flux + unidentified mobility interval | yes (limit 16) | **removed**; the reduction returns exactly `0.0` |
| PHY-020/043 GAB + PHY-048 oil arm | yes | yes, and its 50–95 °C fit range is exceeded at the top of the window **by both paths alike** — limit 13, not a new limit of this extension |
| Span–Wagner saturation surface **inverted** for $T_{sat}(P)$ | forward only | **the construction is new**; the closure under it is not |
| hexane activity convention | fugacity over the liquid root at $(T,P)$ | the ideal $p_h/P_{sat}(T)$, as PHY-020 is fitted and as the certified dry shell already uses |

**No closure is new to the limit. One construction is, and one convention
differs.** The construction is a bisection inverse of a frozen correlation that
reproduces the committed harness anchor bit for bit. The convention is not
free: over the window the ideal activity and the fugacity ratio against the
saturated vapour differ by up to **11.86 %**, the ideal activity and the
convention the binary path uses (against the equation of state's liquid root,
which above the saturation temperature is not a stable phase at all) by up to
**15.57 %**, and carried into the retained-hexane loading through the frozen GAB
the three conventions spread that loading by up to **17.9 %**. That number is
new here and belongs in the supplement beside the existing convention limit 38.

### B.4 What the Faner comparison does and does not establish

Faner's experiment is a pure-hexane one, so the honest question is which of the
limit's closures it puts against a measurement.

| what the limit supplies | does the Faner comparison test it? |
|---|---|
| $T_{sat}(101325\ \mathrm{Pa}) = 68.71451$ °C, the inverted saturation surface | **yes, directly**: the constant-rate plateau, $-0.137$ and $+0.093$ °C from the two digitized plateaux against a $\pm0.7$ °C read band |
| the latent heat at that temperature, through the constant-rate duty | **yes, compositely**: duty $0.82$–$0.98$ and transition time $0.89$–$1.12$ of measured at the independent Whitaker coefficient, nothing fitted |
| the pure-hexane pore-gas **state** (density, fugacity, $\varphi_h$) | **no**: the item-5 march never evaluates it. Stage 2 runs on the certified PHY-019 dry shell, whose own declared scope is "isothermal, pure-n-hexane dry shell — exactly the qualified oracle case PHY-019 freezes", with $a_h = c_g/c_g^{sat}(T)$ |
| the classification "no composition degree of freedom, transport is the total Stefan flux set by heat arrival" | **no measurement can separate it here**: in a one-component gas every pore-gas treatment reduces to the same thing, so the comparison cannot falsify the classification and therefore cannot qualify it either |
| the falling-rate leg, $9.9$ % and $20.6$ % of the window span (intervals $5.3$–$18.9$ and $10.1$–$30.4$ %) | tests the **solid-side** diffusivity law and the certified dry shell, not the limit's pore-gas reduction |

**Said exactly.** The constant-rate and falling-rate agreements are acceptance
evidence for the *composite march* — film coefficient, single-sphere energy
balance, frozen saturation surface, certified dry shell, measured diffusivity
law — and, within that composite, they are direct acceptance evidence for the
one quantity the single-component module contributes, the saturation
temperature. They are **not** a qualification of the single-component pore-gas
reduction, because qualification means a closure exercised inside its measured
range with acceptance tests passed, and the reduction's distinctive content is
the assertion that there is no pore-gas diffusive transport to measure. A
comparison that cannot come out differently under the competing hypothesis is
not evidence for the hypothesis. **Turning this comparison into a qualification
by wording is precisely the move the paper must not make.**

### B.5 Verdict, extension 2

**Not qualified, and not qualifiable by this experiment — but it is a strict
narrowing of the binary path, not a widening of it.** Every closure it uses is
one the binary path already uses, evaluated at a composition at which the
mixing rules degenerate exactly to the pure-component ones; four closure
families, including the two carrying the largest unidentified intervals in the
model ($k_{wh}$ and the pore mobility), drop out exactly. The one number it
contributes to the marched arithmetic is scored against measurement. The open
items it does carry are the activity convention (up to 17.9 % in retained
loading) and the fact that no measurement in this comparison has any sensitivity
to its distinctive claim.

---

## The sentences the paper could carry

Offered for `sec_validation` 5.1 and S10.2; the manuscript is not edited here.

1. *"Both extensions lie outside the qualified domain. Their reach into the
   reported numbers is bounded and was measured: of the two, only the
   single-component module contributes a value to the marched arithmetic, the
   inverted saturation temperature 68.7145 °C, and that value is scored against
   the two digitized plateaux at −0.137 and +0.093 °C."*

2. *"The dry-limit continuation is an extrapolation and is used at one point.
   No held source measures water sorption on soybean or sunflower material
   below a water activity of 0.07, half a decade above the continuation's upper
   edge at 0.046973; the charge of these experiments is water free, so the
   continuation is exercised only at $W_w = 0$, where the retained-water flux
   is identically zero by the absence of the component and where the frozen
   activity has already underflowed to exactly zero for every loading below
   9.69 × 10⁻⁵ kg/kg dry."*

3. *"Where low-activity measurement does exist, the frozen branch lies inside
   its spread. Over twenty-four printed or digitized points between activities
   0.070 and 0.337, on soya meal, whole soya bean and sunflower meal, the
   branch reads 0.31 to 1.21 times the measured loading, straddling unity."*

4. *"The single-component reduction removes closures rather than adding them.
   At $y_h = 1$ the mixture second virial equals the pure n-hexane virial
   bitwise, $\ln\varphi_h$ equals its pure closed form exactly, and the vertex
   state is bitwise invariant over the whole declared interval of the cross
   interaction parameter, while the binary state at the carrier floor moves
   over that interval; the cross virial, the IAPWS-95 water branch, the
   retained-water isotherm and the Maxwell-Stefan flux with its unidentified
   mobility all leave the problem exactly."*

5. *"What converges is the state, not the flux. Thirteen transported
   quantities approach the single-component vertex at first order (two at
   second order) with worst differences below 4 × 10⁻⁶ at a carrier fraction of
   10⁻⁶, across the whole falling-rate window from 68.71 to 136 °C. The
   Maxwell-Stefan flux has no path-independent limit at a composition-
   independent mobility: it tends to −ln 2 along a fixed-ratio approach and
   grows without bound along a fixed-difference one. The reduction's exact zero
   rests on the absence of the second component, not on a limit."*

6. *"The comparison cannot qualify the reduction. In a one-component gas every
   pore-gas treatment reduces to the same statement, so no measurement in these
   experiments has sensitivity to it; the constant-rate and falling-rate
   agreements are acceptance evidence for the composite march and for the
   saturation surface, and are reported as such."*

7. For the supplement, beside limit 38: *"The hexane activity convention is
   open on this path too. The ideal ratio $p_h/P_{sat}$, the fugacity ratio
   against the saturated vapour and the ratio against the equation of state's
   liquid root differ by up to 11.9 % and 15.6 % over the window, and carry a
   spread of up to 17.9 % into the retained-hexane loading."*

---

## What would fully qualify each

Qualification, in this paper's sense, is every closure exercised inside its
measured range with its acceptance tests passed. Neither extension reaches that
today. The shortest honest route for each:

**The dry-limit moisture continuation** needs **one low-activity isotherm point
on the same meal**: a single equilibrium measurement of retained water at a
water activity below 0.047 — that is, at a loading below 0.023 kg/kg dry — on
defatted soybean meal of the class the model is built for, at a temperature
inside or above the frozen branch's 323.15–343.15 K fit range. One point
converts the continuation from an extrapolation into an interpolation over
$[a_{\mathrm{point}}, 0.046973]$ and bounds the branch where it is now free;
three points at one temperature would additionally test the branch's shape
there. A saturated-salt jar at lithium chloride (about 0.11 at 25 °C) is too
wet; the region needs a desiccant or a controlled dry-gas stream, and the
measurement is a mass determination at the 10⁻³ kg/kg level. **Failing that**,
the narrow statement in sentence 2 above is what the evidence supports, and it
should replace the broad one.

**The single-component pore-gas limit** needs **one measured pure-hexane
pore-diffusion series**: a transient uptake or release of n-hexane on the same
meal in a *pure* superheated n-hexane atmosphere at a known superheat and a
known external coefficient, resolved in time, at two or more particle sizes or
two or more film coefficients. Two sizes are what makes it decisive: if the rate
scales with the external area and not with $R^2$, the Stefan/heat-arrival
classification is confirmed against an alternative with a finite pore
diffusivity; if it scales with $R^2$, the classification is falsified and the
reduction must be replaced rather than annotated. The existing Cardarelli
chapter-3 uptake series cannot serve, because they are run in a carrier gas and
are therefore binary. **Failing that**, sentence 6 above is the boundary of what
may be claimed, and the reduction stands as a *forced* modelling statement —
forced by the absence of a second component, not selected against alternatives.

Two cheaper items would also help and need no new experiment: settling the
hexane activity convention (B.3, up to 17.9 % in retained loading), and stating
in the supplement that the flux has no path-independent limit (B.2), which is a
strength of the refusal machinery and currently goes unsaid.

---

## Files

| file | what it is |
|---|---|
| `run_item19_dry_limit.py` | part A; imports item 18's transcribed sources and item 11's digitized sunflower branch |
| `run_item19_single_component.py` | part B; the window continuity, the flux paths, the closure inventory, the convention sizes |
| `make_item19_figure.py` | the three-panel figure |
| `item19_low_activity.csv` | the 28-entry placement table |
| `item19_low_activity_points.csv` | every measured or digitized point of every source, with the branch value and ratio at its activity |
| `item19_continuation_band.csv` | the extrapolated source band over the continuation's own activity interval, with its width |
| `item19_carrier_continuity.csv` | 165 rows: 33 window temperatures × 5 carrier fractions × 13 quantities |
| `item19_flux_paths.csv` | the two approach families to the vertex |
| `item19_closure_inventory.csv` | closure by closure, binary path against the limit |
| `item19_dry_limit.json`, `item19_single_component.json` | the full records with their summaries |
| `item19_extensions.pdf`, `.png` | (a) the placement, (b) what converges, (c) what does not |

Environment: `PYTHONHASHSEED=1`, single-threaded BLAS/OMP/MKL/NUMEXPR,
`PYTHONDONTWRITEBYTECODE=1`, `.venv\Scripts\python.exe`, `PYTHONPATH=src`.
