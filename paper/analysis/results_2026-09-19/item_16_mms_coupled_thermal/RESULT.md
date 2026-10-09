# Item 16 — manufactured-solution convergence for the coupled thermal operator (2026-09-23)

**Executes** mitigation **M14** of `paper/01_particle_jfpe/LIMITS_AUDIT_2026-09-22.md`,
which bears on **limit 5** ("No number exists for the fixed-radius thermal
contract"). The design, written before the campaign, is `DETAIL.md`.

Throughout, `Linf` is the maximum-norm error over cells and `L2` the
volume-weighted two-norm, both taken at the final time against the manufactured
field evaluated at the cell centres of the final grid.

## 0. The rejecting outcome, as it was written before the first run

> An observed order **below 1.8 in space** or **below 0.9 in time** on the
> finest three levels means the coupled thermal operator is **not verified to
> the constant-coefficient standard the paper prints**, and the record says so
> in those words. The audit's own weaker trigger is also carried: an observed
> order **below 1** degrades the thermal discretization claim and is reported
> as a degradation.

**It fired, on one of the three legs.** The fixed-radius wet-core energy
operator and the dry-shell transport operator pass. The **wet-core energy
operator on a moving front-fitted grid fails it in space**: the observed
spatial order is **0.719** in `Linf` and **0.747** in `L2` over 64, 128 and 256
radial cells at a 12 % front recession — below 1.8 and below 1. Its temporal
order is unaffected and reads 1.00. The cause is measured, not guessed, and is
given in Section 4.

## 1. The operator

The continuum statement is the paper's own. In the smooth wet region
Eq. `(localenergy)` with `j_dm = j_Lo = 0` and no relative n-hexane flux
(Eq. `(shellconserved)`) is Fourier conduction of the complete common-datum
internal energy,

```
    d/dt ( rho_dm,p e_w(T, X_h,w^a) ) = (1/r^2) d/dr ( r^2 k dT/dr ) + S(r,t),
    e_w = h_dm(T,w_o) + X_w h_w,l(T,P) + X_h,w^a u_h,l(T,P) - B_h(a_h=1,T,w_o),
```

every datum term of Eq. `(wetenergy)` carried: the PHY-048 reference-anchored
composite meal/oil sensible enthalpy, the retained-water branch, the
Span-Wagner saturated-liquid n-hexane internal energy, and the saturated
binding deficit. The discretization under test is the shipped one,
`wet_core.step` and `wet_core.residual_and_jacobian`: conservative cell-centred
spherical finite volumes with exact areas and volumes, a two-point face flux,
backward Euler, a damped Newton inside the caloric bracket 310 to 380 K, and,
when the front moves, an ALE mesh flux built on the exact swept spherical face
volumes with the cell energy **upwinded**. The dry side is `dry_shell.step`:
the same spatial form with the nonlinear GAB-plus-oil storage law
`C = eps_g c_g + rho_dm,p W_h(a_h,T,w_o)` and a scalar Fick pore flux.

Nothing is re-implemented here. The manufactured source enters each kernel
through its own public `source=` argument; setting it to zero recovers the
shipped operator exactly. The blocked composition screen of
`item_04_fixed_radius_thermal` is not entered, not repaired and not required.

Frozen inputs: `R0 = 0.885e-3 m`, activation at 330 K giving the uniform shell
label `X_h,w^a = 0.20315759488860166`, `k = 0.24 W/(m K)`,
`rho_dm,p = 1159.65 kg/m3`, `X_w = 0.10`, `w_o = W_O_REF`, the bare film
`h = 131.3373155755338 W/(m2 K)`, `D_eff = 4.0e-10 m2/s`, `eps_g = 0.141`.

## 2. Two prerequisite instruments

**Caloric consistency** — may the manufactured source take `de_w/dT` from
`wet_core.specific_heat_capacity`? **Yes.** Over the band the manufactured
fields visit, 334 to 348 K, the analytic capacity and a fourth-order central
difference of `specific_energy` agree to **5.70e-12** relative. Over the whole
caloric bracket the worst disagreement is 3.29e-9, at 312 K, and a step sweep
from 0.4 K down to 0.0125 K walks the difference to the analytic value
(4279.06180, 4279.11692, 4279.12032, 4279.12053, 4279.120540, 4279.1205410
against the analytic 4279.1205410), so that figure is finite-difference
truncation near the `a_h = 1` GAB wall and not a gap in the code. The same
check for `dry_shell.storage_capacity` against `storage`: **2.03e-12**.

**Spatial exactness** — does the field the temporal study uses really leave
only the time integrator's error? **Yes.** A steady `a + b r^2` with the exact
Neumann surface flux and the exact source `6 k b`, pushed through the shipped
ALE residual, gives a scaled residual of **7.60e-19, 3.50e-18, 2.60e-17,
9.35e-17** at N = 12, 24, 48, 96 — machine zero. Corroborated dynamically
(study B0): at 128 steps the error is **0.094287, 0.094519, 0.094576 K** at
N = 12, 24, 48, constant to 0.3 % over a fourfold mesh change, the residue
being only where the outermost cell centre samples the `r^2` profile.

## 3. The measured orders

Orders are the least-squares slope of `ln(error)` against `ln(resolution)` over
the **finest three levels**, with the RMS log residual and `R2` beside them.

### 3.1 Wet-core energy, fixed radius — PASSES

| study | boundary | levels | norm | order | RMS log residual | R2 |
|---|---|---|---|---|---|---|
| A1 | Robin film | N = 16, 32, 64, 128, 256 at 128 steps | `Linf` | **2.0036** | 1.280e-3 | 0.999999 |
| A1 | Robin film | " | `L2` | **2.0057** | 1.051e-3 | 0.999999 |
| A2 | Dirichlet trace | N = 16, 32, 64, 128 at 128 steps | `Linf` | **1.9963** | 1.125e-3 | 0.999999 |
| A2 | Dirichlet trace | " | `L2` | **2.0024** | 1.723e-4 | 1.000000 |

A1 levels, `Linf` (K): 8.566390e-2, 2.177909e-2, 5.472013e-3, 1.368281e-3,
3.402871e-4; pairwise orders 1.976, 1.993, 2.000, 2.008.

*Temporal contamination is measured, not assumed.* The manufactured field is
exactly linear in `t`, so the backward-Euler storage truncation is starved of
its leading term. Study A3 repeats the finest mesh at half the number of steps:
at N = 256 the error moves from 3.402871e-4 to 3.369768e-4 K, so the whole
time-discretization contribution at the reported level is **0.97 %** of the
spatial error (0.96 % in `L2`), which can shift the finest pairwise order by at
most about 0.014.

### 3.2 Wet-core energy, time refinement — PASSES

Study B1, on the radially quadratic field the spatial operator reproduces
exactly, N = 24, exact Neumann surface flux, 32 to 512 steps:

| dt (s) | 1.25e-2 | 6.25e-3 | 3.125e-3 | 1.5625e-3 | 7.8125e-4 |
|---|---|---|---|---|---|
| `Linf` (K) | 0.374031 | 0.188363 | 0.094519 | 0.047344 | 0.023693 |

Order **0.9981** in both norms (RMS log residual 2.10e-4 and 2.08e-4,
`R2 = 1.000000`); pairwise 0.990, 0.995, 0.997, 0.999; self-convergence 0.984,
0.992, 0.996. This is backward Euler, and the caloric nonlinearity does not
degrade it.

### 3.3 Wet-core energy, prescribed receding front — FAILS in space

Dirichlet interface trace at `r = s(t)`, with `s(t) = R0 (1 - sigma t/t_end)`.

| sigma | norm | N = 16 | 32 | 64 | 128 | 256 | finest-three order |
|---|---|---|---|---|---|---|---|
| 12 % | `Linf` (K) | 1.002866e-1 | 6.915456e-2 | 4.217565e-2 | 2.508890e-2 | 1.556861e-2 | **0.7189** (R2 0.9994) |
| 12 % | `L2` (K) | 3.541483e-2 | 2.171587e-2 | 1.244962e-2 | 7.183085e-3 | 4.421691e-3 | **0.7467** (R2 0.9987) |
| 1.2 % | `Linf` (K) | 6.981754e-2 | 1.332589e-2 | 1.234200e-3 | 1.420355e-3 | 1.167843e-3 | not readable |
| 1.2 % | `L2` (K) | 9.332097e-3 | 1.675846e-3 | 5.670340e-4 | 4.339662e-4 | 3.021303e-4 | not readable |

At the 12 % recession the power law is clean — `R2` of 0.9994 and 0.9987 over a
fourfold mesh range — and its exponent is **0.72 to 0.75**. It is not a
pre-asymptotic artefact: the exponent does not rise with refinement, it falls
(pairwise `Linf` 0.536, 0.713, 0.749, 0.688).

At the 1.2 % recession **no order may be quoted**, and the record says so
rather than printing one. The error is non-monotone — it dips at N = 64 and
rises again at N = 128 — because the second-order interior term and the
first-order mesh term carry opposite signs and cancel near that mesh. The
least-squares fit over the finest three levels returns 1.615 with an RMS log
residual of 0.594 and `R2 = 0.703`, which is the fit reporting its own
invalidity. The `L2` pairwise sequence 2.477, 1.563, 0.386, 0.522 shows the
crossover: second order while the mesh is coarse, first order and worse once
the mesh term dominates.

**Time on the moving grid is unaffected.** Study D1 (N = 48, 32 to 512 steps,
12 % recession) reads **0.991, 0.995, 0.998** by three-level self-convergence.
Its direct manufactured error reads an apparent order of 0.055; that number is
the spatial error floor of the moving-grid operator at N = 48, not a temporal
order, and it is recorded here only so that it is not mistaken for one later.

### 3.4 Dry-shell hexane transport — PASSES

| study | refinement | norm | order | RMS log residual | R2 |
|---|---|---|---|---|---|
| E1w | N = 8 to 256 at 2048 steps, `t_end = 16 s` | `Linf` | **1.9454** | 6.433e-3 | 0.999966 |
| E1w | " | `L2` | **1.9945** | 1.103e-3 | 0.999999 |
| E2 | 16 to 256 steps at N = 48, `t_end = 4 s` | `Linf` | **0.9983** | 1.995e-4 | 1.000000 |
| E2 | " | `L2` | **0.9724** | 3.074e-3 | 0.999969 |

E1w pairwise `L2`: 1.272, 1.722, 1.957, 1.991, 1.998. The two coarsest levels
are pre-asymptotic because the field's highest radial mode, `3 pi r / R0`, has
fewer than three cells per half-wave there; this is why the ladder was extended
to six levels rather than read at four.

E2 levels, `Linf` (kg/m3) at dt = 0.25, 0.125, 6.25e-2, 3.125e-2, 1.5625e-2 s:
2.405492e-2, 1.215307e-2, 6.106792e-3, 3.058362e-3, 1.530372e-3; pairwise
0.985, 0.993, 0.998, 0.999; self-convergence 0.977, 0.989, 0.994. Its sixth
level, dt = 7.8125e-3 s, refuses (Section 5).

The manuscript currently carries this operator's order only in the **steady,
storage-free, linear** `steady_diffusion_solve` form. What is added here is the
transient solve with the full nonlinear GAB-plus-oil storage law inside the
Newton loop.

### 3.5 Ledgers, residuals and reproducibility

The accumulated discrete ledger
`L = (E_end - E_0) + sum_n dt H_surface - sum_n dt sum_i V_i S_i`, normalized by
the largest of its three terms, over **every** wet-core level of every study:
worst **1.787e-12** (study B1, N = 24, 32 steps). Over every dry-shell level:
worst **3.926e-13** (study E1, N = 8). This is the scheme's exact telescoping
identity, so it measures conservation and not accuracy.

The worst scaled Newton residual, re-evaluated on the shipped residual at each
accepted state across the whole campaign, is **3.43e-14**, a factor 29 inside
the 1e-12 this driver requested — itself tighter than the shipped 1e-11. That
scaled residual bounds the temperature by which an accepted state can miss its
own discrete equation at about `3.43e-14 x 340 K = 1.2e-11 K`, seven orders
below the smallest error reported above, so no order here is solver-limited.

Studies C1 and C1b ran in separate concurrent processes and share the levels
N = 64 and N = 128. Not only their error norms but their whole final
temperature fields are bit-identical, element for element; the same holds for
E1 and E1w at all four of their shared levels. The figure builder asserts this
rather than assuming it.

## 4. The cause of the moving-grid failure, measured

The interior ALE mesh flux is `H_f = A_f k (T_L - T_R)/dr - Q_f U_upwind`, with
`U_upwind` the **cell** energy density of the upwind side rather than the
energy at the face. Studies Cd1 and Cd2 inject the exact manufactured field
into that term with no solve at all and accumulate the resulting per-cell
temperature defect.

| sigma | quantity | N = 16 | 32 | 64 | 128 | 256 |
|---|---|---|---|---|---|---|
| 1.2 % | max defect (K) | 2.438376e-2 | 1.226403e-2 | 6.132744e-3 | 3.064540e-3 | 1.531767e-3 |
| 1.2 % | pairwise order | — | 0.992 | 1.000 | 1.001 | 1.000 |
| 12 % | max defect (K) | 2.448087e-1 | 1.224052e-1 | 6.123913e-2 | 4.969565e-2 | 4.714226e-2 |
| 12 % | pairwise order | — | 1.000 | 0.999 | 0.301 | 0.076 |

Two things are established.

**The truncation is exactly first order in the mesh spacing and exactly
proportional to the mesh speed.** At the slow recession it halves at every
level to four figures, and the two recession rates differ by 10.0398 at N = 16
against a prescribed ratio of exactly 10.

**It stops converging at all in one cell.** At the fast recession the location
of the maximum moves, between N = 64 and N = 128, from `r/R = 0.2734` to
`r/R = 0.9961` — the cell adjacent to the receding boundary — and from there
the defect is flat under refinement, 0.049696 to 0.047142 K. The mechanism is
structural rather than incidental. For an interior cell the defects of its two
faces are differenced, so the `O(h)` face defects telescope down to `O(h^2)`,
which over an `O(h)` cell volume leaves the familiar `O(h)` error. The cell
next to the moving boundary has no partner, because the surface rate
`H = H_conduction - Q_surface U(T_interface)` is exact there, so its single
undifferenced `O(h)` face defect is divided by an `O(h)` cell volume and leaves
an `O(1)` error rate confined to that one cell. Conduction then spreads and
partly relaxes it, which is why the **solved** error still falls — at exponent
0.72 to 0.75 rather than at 1 or at 0.

This is a property of the shipped discretization and not of the manufactured
solution: the same field, on the same operator, with the same boundary kind and
the same number of steps, converges at order 2.00 the moment the front is held
still (study A2 against study C1).

## 5. Refusals, reported and not routed around

* `dry_shell.step` **refuses below a minimum time step**. Its shipped residual
  norm is an absolute density rate, whose floating-point floor is about
  `eps C / dt`, so a small enough `dt` puts the shipped 1e-12 tolerance below
  what the arithmetic can reach and the Newton line search fails inside the
  admissible storage domain. The mesh study accepts `dt = 7.8125e-3 s` at every
  level (`out/E1.json`, `out/E1w.json`) and refuses at step 0 at every level for
  `dt = 3.90625e-3 s`, `1.953125e-3 s` and `9.765625e-4 s` (`out/E1c.json`,
  `out/E1d.json`, `out/E1r.json`). The floor is **state-dependent**, as an
  absolute tolerance must be: the time study, on its own field, accepts down to
  `dt = 1.5625e-2 s` and then refuses at `dt = 7.8125e-3 s` — the same step the
  mesh study accepts — partway through its march, at step 189 of 512
  (`out/E2.json`).
* The first pass of this campaign hit the same wall from the other side. The
  driver had asked for `tol = 1e-13`, tighter than shipped, and every mesh level
  of the spatial study then refused partway through its march, at steps 61 to
  118 of 128, while the time study refused at its three finest levels. The
  repair was to use the shipped default, not to relax anything.
* No wet-core level refused anywhere in the campaign.
* No tolerance, budget, seed family, matrix or physical law was relaxed, and no
  refusal was worked around.

## 6. Verdict

1. **The fixed-radius coupled thermal operator is verified to the
   constant-coefficient standard the paper prints.** Space 2.00 under both the
   physical Robin film closure and the Dirichlet interface trace; time 1.00;
   ledgers at 1.8e-12; and the nonlinear common-datum caloric law, including
   its binding deficit, is in the loop at every level. Limit 5's "no number
   exists" no longer holds **for the discretization**.
2. **The dry-shell transport operator is verified** in its full transient
   nonlinear-storage form, at space 1.95 to 1.99 and time 0.97 to 1.00, which
   extends the manuscript's steady scalar-Fick entry.
3. **The same thermal operator on a moving front-fitted grid is not verified
   in space.** At an engineering recession rate the observed order is 0.72 to
   0.75, below the 1.8 of the rejecting outcome and below the audit's own
   degradation trigger of 1. The temporal order is unaffected at 1.00. The
   cause is the upwinded ALE mesh flux and, specifically, its undifferenced
   defect in the cell adjacent to the exactly-treated moving boundary.

## 7. What this does and does not establish

* It does **not** score the gate1h fixed-radius thermal **contract**, the
  resolved-versus-lumped comparison of `item_04`. That contract is unchanged,
  unscored, and still behind the blocked composition screen. What limit 5
  acquires here is the discretization's observed orders, which that contract
  was one of several routes to.
* It does **not** exercise the front jump conditions Eqs. `(rhcomp)`,
  `(rhenergy)`, `(frontvelocity)` or the complementarity system: the front is
  prescribed in every run. Nothing here says the coupled two-region solve
  converges.
* It does **not** reproduce the manuscript's own Table 1 row "Manufactured
  spherical conduction, wet core — order 2 — moving front-fitted grid". That
  row is not contradicted on its own terms, because the mesh speed at which it
  was measured is not stated; what this record establishes is that the
  moving-grid spatial order **depends strongly on the recession rate**, so the
  row's scope needs that rate beside it.
* It carries **no physical claim whatsoever**. An observed order certifies a
  discretization against the continuous equations that were discretized. No
  transport coefficient, no isotherm, no film coefficient and no plant quantity
  is identified, calibrated or validated by anything in this folder.
* The remedy for finding 3 — a face-interpolated or otherwise
  boundary-consistent ALE mesh flux — touches the certified core and is a
  Class-B decision. It is named here and not taken.

## 8. The sentence the paper could carry

> The coupled thermal operator, the wet-core energy equation
> \eqref{eq:wetenergy} on the common caloric datum with its binding potential,
> advanced by the production kernel, converges on manufactured solutions at
> observed order $2.00$ in space (16 to 256 radial cells, the Robin film and
> the Dirichlet interface trace, maximum and volume-weighted two-norms,
> least-squares log residual below $1.3\times10^{-3}$) and $0.998$ in time,
> with the accumulated energy ledger closing to $1.8\times10^{-12}$ of its
> largest term at every level; the dry-shell transport operator converges at
> $1.95$--$1.99$ and $0.97$--$1.00$ in its full transient nonlinear-storage
> form. Under a prescribed receding front the temporal order is unchanged at
> $0.99$, but the spatial order falls to $0.72$ and $0.75$ in the two norms at
> a \SI{12}{\percent} recession, because the arbitrary Lagrangian--Eulerian
> mesh flux upwinds the cell energy and its truncation, first order in the mesh
> spacing and proportional to the mesh speed, is undifferenced in the cell
> adjacent to the moving boundary; the discretization contribution on a moving
> front is therefore bounded by that exponent and not by the fixed-topology
> order.

A shorter form, if only one sentence fits: *manufactured-solution orders for
the coupled thermal operator are 2.00 in space and 1.00 in time at a fixed
radius, and 0.72 to 0.75 in space at a 12 % front recession, the loss traced to
the upwinded ALE mesh flux at the cell adjacent to the moving boundary.*

## 9. Files

| file | what it is |
|---|---|
| `DETAIL.md` | the design, written before the campaign, with the rejecting outcome |
| `mms_coupled_thermal.py` | the driver; `--study <tag>` runs one study and writes `out/<tag>.json` |
| `make_item16_figures.py` | reads `out/*.json`, writes the figure and the flat table |
| `out/caloric.json`, `out/exactness.json` | the two prerequisite instruments |
| `out/A1.json`, `A2`, `A3`, `B0`, `B1` | wet core, fixed radius |
| `out/C1.json`, `C1b`, `C2`, `C2b`, `D1` | wet core, prescribed receding front |
| `out/Cd1.json`, `Cd2` | the ALE mesh-flux truncation diagnostic |
| `out/E1.json`, `E1b`, `E1c`, `E1d`, `E1r`, `E1w`, `E2` | dry-shell transport, including the four refused time steps |
| `fig_mms_coupled_thermal.pdf`, `.png` | six panels: (a) wet core fixed radius, (b) wet core exact in space, (c) wet core receding front, (d) receding front self-convergence, (e) dry shell mesh, (f) dry shell time |
| `convergence_levels.csv` | every level of every study, flat, with its norms, ledger and refusal |

Environment on every call: `PYTHONHASHSEED=1`, `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`,
`PYTHONDONTWRITEBYTECODE=1`, `PYTHONIOENCODING=utf-8`,
`PYTHONPATH=<repo>/src`, interpreter `.venv\Scripts\python.exe`. Total compute
6533 s of single-thread work summed over the twenty-one studies, run as
independent concurrent processes; the longest single study, A1, is 1101 s.
Nothing under `src/`, `tests/`, `scripts/`, `tools/` or `docs/` was edited, and
nothing in the manuscript was touched.
