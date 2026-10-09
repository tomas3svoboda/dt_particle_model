# Detail: the rig's own hydrodynamic state, the mass convention, and the
# re-march (2026-09-22)

This folder supersedes `rerun_2026-09-21b/` **for stage 1 only**. That folder,
`rerun_2026-09-21/` and `rerun_2026-09-20/` beside it are byte-identical and
nothing in any of them was deleted. Stage 2 — the certified dry shell, the
frozen PHY-019 value, the measured Cardarelli law at both of its evaluation
edges — the oil-solution arm and the residual floor are unchanged and were
recomputed here on one run of one code path.

`physically_qualifying: false. plant_predictive: false.`

The three mitigations executed are M7, M2 and M3 of
`paper/01_particle_jfpe/LIMITS_AUDIT_2026-09-22.md`, which close or bound
limits 44, 7, 8 and 9 of that audit's table.

---

## 1. The rig's own hydrodynamic state (M7, M2)

`rig_state.py` reconstructs it from printed numbers only. Table 4.7 of the
thesis gives, per run, the approach velocity, the gas temperature, the
coefficient, the Prandtl number and **both** Reynolds numbers, so

```text
eps = 1 - Re/Re_eps            (Eqs. 4.22 and 4.23, as printed)
g   = h/Nu = lambda_g/DEE      (Eq. 4.20)
k   = Nu/Nu_eps = (1-eps)/eps  (Eqs. 4.20 and 4.21)
```

follow from the table alone. A second, independent route to each run's voidage
— the Nusselt pair, `eps = 1/(1 + Nu/Nu_eps)` — agrees with the first to
**1.6e-4** over the twelve hexane-vapour runs and to 6.1e-4 over all
thirty-six. Interpolating linearly in the printed gas temperature across the
six hexane-vapour runs of each meal (extrapolation is refused and the code
raises) puts the four marched conditions at:

| condition | v_g (m/s) | Re | Re_eps | Pr | eps | lambda_g (W/m K) |
|---|---:|---:|---:|---:|---:|---:|
| faner2008 sunflower 136 C | 0.1285 | 74.38 | 194.71 | 0.8048 | 0.6173 | 0.02392 |
| faner2008 soybean 136 C | 0.1258 | 67.95 | 236.09 | 0.8049 | 0.7121 | 0.02364 |
| faner2019 soybean 120 C | 0.1232 | 72.17 | 229.47 | 0.8032 | 0.6854 | 0.02200 |
| faner2019 sunflower 100 C | 0.1165 | 81.71 | 230.72 | 0.8008 | 0.6457 | 0.02026 |

### The project's property authorities enter only as a check, and they close it

The project holds a certified hexane equation of state and **no** hexane-vapour
viscosity and **no** hexane-vapour thermal conductivity. Inventing either would
put an unpinned number inside a prediction, so `rig_state.py` does not. It
inverts the source's own columns instead — `mu_g = DEE rho_g v_g / Re` with
`rho_g` from the certified equation of state at the printed temperature and the
assumed atmospheric pressure, and `lambda_g = g DEE` from Eq. 4.20 — and closes
the Prandtl number against the same equation of state's heat capacity:

```text
implied mu_g      7.85e-6 to 8.92e-6 Pa s   (n-hexane vapour, 92 to 138 C)
implied lambda_g  0.0193 to 0.0241 W/(m K)
Pr(computed)/Pr(printed)   0.9996 to 1.0151   over all twelve runs
```

The printed table therefore closes against an authority outside the thesis to
within 1.5 %. The pressure remains a READING: the source prints none.

## 2. The coefficient table, and what each entry is (M2, M7)

All in W m-2 K-1, evaluated at the state above on the thesis's own definitions
(4.20) to (4.23). Nothing here is fitted to a drying curve.

| coefficient | 2008 sf 136 | 2008 sb 136 | 2019 sb 120 | 2019 sf 100 | what it is |
|---|---:|---:|---:|---:|---|
| demanded by the measurement, printed mass = charge | 97.1 | 85.5 | 86.5 | 91.6 | measurement inversion |
| demanded, printed mass = dry meal | 168.5 | 148.4 | 129.0 | 137.4 | measurement inversion |
| matched experiments, mean | 97.3 | 94.2 | 75.7 | 76.7 | same-rig measurement |
| matched experiments, band | 91.1-106.5 | 84.6-103.9 | 58.2-103.9 | 73.9-78.6 | same-rig measurement |
| **Whitaker (4.18)** | **96.9** | **75.7** | **78.6** | **80.0** | **independent literature** |
| **Bird et al. (4.17)** | **183.8** | **188.3** | **180.3** | **162.8** | **independent literature** |
| **Bradshaw-Myers (4.19)** | **247.5** | **190.3** | **198.1** | **201.5** | **independent literature** |
| Ranz-Marshall (4.14), isolated sphere | 83.6 | 86.7 | 82.4 | 73.1 | independent literature |
| thesis bed fit (4.27) | 103.9 | 81.1 | 84.3 | 85.8 | same-rig fit |
| thesis particle fit (4.26) | 92.2 | 93.9 | 90.3 | 82.1 | same-rig fit |
| frozen Coletto pair at its OWN state | 131.3 | 131.3 | 131.3 | 131.3 | borrowed |
| frozen Coletto pair at the RIG's state | 104.0 | 81.2 | 84.4 | 85.9 | borrowed |

`INDEPENDENCE` in `rig_state.py` carries these labels in the data, not only in
this table. Bird, Whitaker and Bradshaw-Myers are fitted by other authors to
other beds and to no datum of either Faner source; Eqs. 4.26 and 4.27 are the
thesis's own regressions on the same thirty-six coefficients that its Eq. 4.9
had already inferred from each run's own evaporation rate.

### M7: limit 44 closes, and the cause is measured

The frozen pair is read at `Re_eps = 45.93` and `Pr = 1.0796`; the rig sits at
`Re_eps` 194.7 to 236.1 and `Pr` 0.801 to 0.805, a Reynolds ratio of **4.24 to
5.14**. Read at the rig's own state the same correlation, with the repository's
own prefactor 0.6949 unchanged, returns **81.2 to 104.0** instead of 131.3 — a
factor 0.62 to 0.79 — and lands within 0.12 % of the thesis's own Eq. 4.27,
which is what Coletto's B.7 is a transcription of. **The 1.33 to 1.53
over-prediction the 2026-09-21 record attributed to the correlation was the
evaluation state.** The rejecting outcome written before the arithmetic (that
re-evaluation would not move the prediction toward the measurement, in which
case the area convention carried it) did not fire.

No constant of the model was changed. The prefactor deviation stands as a
documented reading: the thesis prints 0.6941, the repository binds 0.6949,
1.15e-3 relative, 0.115 % on every coefficient built on it.

### M2: the rejecting outcome did not fire either

It was written first: "if the three correlations at the rig's state do not
bracket the source's measured coefficients, the area or the group convention is
wrong and the stage-1 posing is withdrawn." They do bracket it. At all four
conditions the coefficient the measured duty demands on the printed footprint
lies between the lowest (Whitaker) and the highest (Bradshaw-Myers), and
Whitaker gives 0.998, 0.885, 0.909 and 0.873 of that demand — *slightly low at
three of four*, which is the source's own printed verdict on Whitaker, reached
here independently of it. Against the mean of the matched experiments Whitaker
gives 0.996, 0.803, 1.038 and 1.043; the tabulated mean lies inside the
three-correlation band at the two thesis conditions and 3.8 % and 4.1 % below
its low edge at the two journal conditions.

The band is wide — Bradshaw-Myers runs 2.0 to 2.6 times the demand and Bird 1.8
to 2.1 — and is reported as the spread of three published correlations, never as
an uncertainty interval.

## 3. The sample-mass convention (M3)

`mass_convention.py`. Equation 4.25 defines the holder voidage from the printed
sample mass and the printed envelope density; the voidage is already fixed per
run by the Reynolds pair, so the equation inverts for the cell volume, which is
fixed hardware and therefore a test of the whole set.

### What is settled

| group | n | V_celda (m3) | full spread | depth | in particle diameters | inside the printed 3-5 bound |
|---|---:|---:|---:|---:|---:|---|
| nitrogen holder | 18 | 6.5999e-5 | 3.6e-3 | 5.999 mm | 3.07-3.34 | yes |
| vapour holder | 18 | 9.0201e-5 | 1.3e-3 | 8.200 mm | 4.20-4.56 | yes |
| the 12 hexane-vapour runs | 12 | 9.0199e-5 | **1.3e-4** | 8.200 mm | 4.20-4.56 | yes |

Two holders, each recovered to four significant figures from eighteen
independent runs, each inside the source's own printed thin-layer bound. The
printed mass and the printed envelope density are therefore on **one common
basis**, and the voidage Table 4.7 carries is the voidage Eq. 4.25 produced.

### What is not settled, and why it cannot be settled this way

Stated in the module before the arithmetic was run. A non-swelling particle's
envelope volume does not change when solvent enters its pores, so
`rho_wet = (1+X_0) rho_dry` exactly and

```text
W / rho_env  ->  (1+X_0) W / ((1+X_0) rho_env)  =  W / rho_env
```

Equation 4.25 is **invariant** under a joint change of basis. No arithmetic on
it can separate the wet charge from the dry meal.

### The two external routes, both reported, pointing opposite ways

**Toward the dry-meal reading.** Read as dry meal the printed mass gives the
layer a solid fraction of 0.266 to 0.469 and a dry bulk density of 308 to
482 kg/m3 — a poured bed of oilseed flakes. Read as the charge it gives
0.148 to 0.260 and 171 to 268 kg/m3, looser than any poured meal bed. And at
the family's printed charge loading X_0 = 0.80 (p. 102) a solvent-bearing
envelope density of 1023 or 1160 kg/m3 would put 0.74 or 0.84 of the particle's
own envelope volume inside its pores as liquid hexane at 615 kg/m3, forcing a
solid skeleton of 2180 or 3970 kg/m3, which no oilseed-meal solid has: the
printed densities are solvent-free, and on a matched pairing so is the printed
mass.

**Toward the total-charge reading.** Equation 4.8 inverted — the 2026-09-21b
`circularity_probe` — reproduces the digitized slope of the plotted curve to
1 to 10 % at the closest matched run of every condition only on the charge
reading; on the dry reading every matched experiment implies 0.38 to 0.68 of
the plotted slope. And only on the charge reading does the coefficient the
plotted curve demands, 85.5 to 97.1, agree with the coefficient the same rig
measured, 75.7 to 97.3, and with the independent Whitaker prediction, 75.7 to
96.9; on the dry reading the demand rises to 129.0 to 168.5 and exceeds all
three.

### The verdict

Neither route is decisive, because the thesis never states that a tabulated run
is a plotted run and never prints the plotted runs' own charge: a plotted charge
of 0.58 to 0.67 of the tabulated one reconciles the dry reading exactly. The
convention is recorded as **open with its factor**, status
`CONVENTION_UNESTABLISHED_BOUNDED`, and **both transfer areas are marched**:

| condition | charge reading (m2/kg dry) | dry-meal reading | ratio |
|---|---:|---:|---:|
| faner2008 sunflower 136 C | 0.5337 | 0.3076 | 0.5764 |
| faner2008 soybean 136 C | 0.6377 | 0.3675 | 0.5764 |
| faner2019 soybean 120 C | 0.5466 | 0.3667 | 0.6708 |
| faner2019 sunflower 100 C | 0.5497 | 0.3667 | 0.6669 |

That is a move of 1/(1+X_0), 0.58 to 0.67, on every constant-rate rate and its
reciprocal on every crossing time. The march carries both for the two headline
coefficients rather than choosing.

## 4. The re-march (stage 2 unchanged)

Thirteen stage-1 cases per condition, eight diffusivity cases each, one run of
one code path, wall time 3123.7 s. Stage 2 is the same certified dry shell at
the same eight diffusivity cases read from
`../rerun_2026-09-21/diffusivity_law.json`; the oil arm and the residual floor
are imported unchanged from `../rerun_2026-09-21/oil_solution_arm.py`. The
`frozen_1.0` rows reproduce the borrowed-correlation reference exactly (1.33,
1.53, 1.51, 1.35 and crossings 0.66, 0.58, 0.58, 0.68), so the before/after is
on one code path.

### Constant rate simulated over measured, and the crossing time over the
### measured crossing

| stage-1 case | 2008 sf 136 | 2008 sb 136 | 2019 sb 120 | 2019 sf 100 |
|---|---|---|---|---|
| **lit_whitaker** (independent) | **0.98 / 0.89** | **0.88 / 1.01** | **0.91 / 0.97** | **0.82 / 1.12** |
| lit_bird (independent) | 1.86 / 0.47 | 2.19 / 0.40 | 2.08 / 0.42 | 1.67 / 0.55 |
| lit_bradshaw_myers (independent) | 2.51 / 0.35 | 2.22 / 0.40 | 2.28 / 0.39 | 2.07 / 0.45 |
| source_eq427 (same-rig fit) | 1.05 / 0.83 | 0.95 / 0.94 | 0.97 / 0.91 | 0.88 / 1.05 |
| source_eq426 (same-rig fit) | 0.93 / 0.94 | 1.09 / 0.81 | 1.04 / 0.85 | 0.84 / 1.09 |
| source_central (same-rig measurement) | 0.99 / 0.89 | 1.10 / 0.81 | 0.87 / 1.01 | 0.79 / 1.17 |
| frozen_1.0, own state (borrowed) | 1.33 / 0.66 | 1.53 / 0.58 | 1.51 / 0.58 | 1.35 / 0.68 |
| frozen_rig_state (borrowed) | 1.05 / 0.83 | 0.95 / 0.94 | 0.97 / 0.91 | 0.88 / 1.05 |
| lit_whitaker on the dry-meal area | 0.57 / 1.55 | 0.51 / 1.74 | 0.61 / 1.45 | 0.55 / 1.68 |
| source_central on the dry-meal area | 0.57 / 1.54 | 0.63 / 1.40 | 0.59 / 1.51 | 0.53 / 1.75 |

`frozen_rig_state` and `source_eq427` agree to the third digit at every
condition, as they must: Coletto's B.7 is a transcription of Faner's 4.27 and
the only difference is the prefactor, 0.6949 against 0.6941.

### The falling-rate window, at the law read at the surface loading

RMS as a fraction of the window span, band counts of the fifteen points that
carry a digitization reading band.

| condition | stage 1 | RMS/span | in band |
|---|---|---:|---:|
| 2008 sunflower 136 C | lit_whitaker | 0.0994 | 14/15 |
| 2008 sunflower 136 C | source_central | 0.0992 | 14/15 |
| 2008 sunflower 136 C | frozen_rig_state | 0.0968 | 15/15 |
| 2008 sunflower 136 C | lit_whitaker, dry area | 0.7789 | 0/15 |
| 2008 soybean 136 C | lit_whitaker | 0.2059 | 3/15 |
| 2008 soybean 136 C | source_central | 0.1629 | 7/15 |
| 2008 soybean 136 C | frozen_rig_state | 0.1874 | 3/15 |
| 2008 soybean 136 C | lit_whitaker, dry area | 1.1061 | 0/15 |

**Reported as it came out.** On the sunflower trace the independent prediction
is as good as the same-rig coefficient was, 9.9 % of the window span with 14 of
15 points in the reading band. On the soybean trace it is worse, 20.6 % against
16.3 %, because Whitaker gives that condition 75.7 W m-2 K-1 against the
tabulated 94.2 and the longer constant-rate leg grafts the tail on later; the
same mechanism the 2026-09-21b record identified, working the other way. The
best of the three at the sunflower condition is the frozen correlation read at
the rig's state, 9.7 % with 15 of 15 in band, which is the source's own bed fit
in another transcription.

### The identified scalar, re-identified on the new stage 1

`identification_crosscheck_source_h.py`, stage 1 `lit_whitaker`, the same grid
and the same standard-error construction. Identified on the thesis sunflower
trace: `D_eff = 1.1654e-9 +/- 1.4848e-10 m2/s`, 12.74 % relative standard
error, 14 degrees of freedom, 2.914 times the frozen reference. It reproduces
the trace it was fitted to at 15.84 % of the window span (7 of 15 in band) and
the held-out soybean trace at 26.45 % (5 of 15). Both are worse than the
parameter-free law at the surface loading, and the held-out residual is now
1.67 times the identified one, so the scalar transfers less well than the law.
The reverse direction is also recorded: identified on soybean,
`1.5193e-9 +/- 2.5788e-10` (16.97 %), 24.60 % fitted and 18.26 % held out.

## 5. What did not move

No constant of the model was changed: not the 0.6949 prefactor, not the
repository's `Re_eps` convention, not a tolerance, budget, gate or seed family.
Nothing under `src/`, `tests/`, `scripts/`, `tools/`, `docs/` or
`docs/evidence/` was touched. The refused register is 54 elements in twelve
classes, unchanged in count; no element was removed.

`status: NOT_QUALIFYING. physically_qualifying: false. plant_predictive: false.
production_wired: false. qsc10_complete: false.` This record does not complete
drying-rate validation, establish plant prediction, authenticate industrial
data, release production software or complete QSC-10.
