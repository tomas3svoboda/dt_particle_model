# Item 32 — the space-resolved law marched at each species' own sphere

`physically_qualifying: false. plant_predictive: false. nothing fitted.`

## 0. Rejecting outcomes, written 2026-09-27 before any declared-geometry march

Brief (sixth referee read, minor item 12): the space-resolved reading of the
carried diffusivity law, `D(T, X(r))` inside the flux divergence (the one
reading no evaluation choice affects; `dry_shell_variable_coefficient.py`,
driver `item_05_faner_comparison/rerun_2026-09-21c_dxr/`), has been marched
only at the 0.885 mm sphere. March it at the declared species-own spheres of
item 31, R_d,sun = 0.74 x 1.95 mm / 2 = 0.7215 mm and R_d,soy =
0.74 x 1.80 mm / 2 = 0.6660 mm, on both thesis traces (136 C, temperature
case `sensible`), on the leg the paper's printed values use (the borrowed
frozen Coletto film pair at its central multiplier 1.0) and on the Whitaker
leg (the stage-1 film coefficient and crossing of `rerun_2026-09-22_rig`,
the leg of items 24 to 31), with the film transfer area per kg dry held at
the committed sample-holder value for every radius.

* **Q0 (identity first).** The committed space-resolved record
  (`dxr_tail_residuals.csv`, 24 rows: 4 conditions x 2 laws x 3 film
  multipliers; the simulated rows of `dxr_traces.csv`; `dxr_coefficient_profiles.csv`)
  is re-marched through this item's driver, with its step recorder in place,
  and compared cell by cell by string equality. The paper's 0.770 (soybean)
  and 0.697 (sunflower) are two of these rows. One mismatch stops the item.
  Item 31's two reported-level decimated traces are re-scored with the window
  rule used here and must return item 31's 13.3 and 19.7 per cent to the
  float.
* **R1 (the space-resolved march refuses inside the scored window at the
  declared radii).** Fires if any declared-radius march (either leg, either
  trace, 60 or 120 cells) ends before the last scored sample (300.1 s after
  the start of drying). Refusals after the window are recorded verbatim and
  do not fire R1. Nothing (tolerance, retry budget, stride floor, horizon,
  law) is relaxed.
* **R2 (the space-resolved residual at the declared radii is not reproducible
  to float identity when the driver is re-run on the paper's existing
  configuration).** Fires if Q0 fails, or if the 60-cell 0.885 mm borrowed-leg
  marches of the declared phase (the paper's configuration run through the
  same code path as the declared radii) differ in their window residual from
  the committed rows by any amount.
* **Reported without a threshold:** per march the window RMS over the span
  (decimated "as printed", the committed rule, and full-step), points of 15
  in the reading band, strides, the frozen-coefficient fallback count, the
  ledgers (per-stride inventory against the discrete face outflow, between
  strides across the constant-storage remap, whole march), the within-particle
  coefficient ratio, the volume fraction of the particle above the law's
  measured loading range in the window; the 60-to-120 cell difference at
  every radius and leg.

## 0b. Declared extension, written 2026-09-27 after the declared phase and before the extension ran

The declared phase (section 2) fired R1: at 120 cells every march refuses
before the end of the window, the 0.885 mm sphere included, and at 60 cells
the soybean species-own march refuses at 88.7 s (borrowed) and 147.1 s
(Whitaker). The refusal is recorded verbatim and nothing below replaces it.
The committed driver's first stride is 1.8 to 3.2 s and its ceiling 5 s,
independent of the mesh, and its retry ladder retries a refused stride at a
LARGER stride. As a declared extension, not a repair, the same 16 marches
are run with the stride ceiling passed through the committed argument
`dt_max_s` at 5/4 s and 5/16 s (item 24's construction of the ceiling at time
levels 4 and 16; the first-stride rule, the growth, the retry ladder, the
tolerance, the budgets, the horizon and the law unchanged). The extension
cannot un-fire R1. It answers only whether the refusal is a property of the
stride the committed driver takes, and, where a 60/120 pair covers the
window at one ceiling, it gives the mesh difference at that ceiling.

---

## 1. Identity first (Q0): passed

**3,732 of 3,732 checks exact over 24 marches** (`identity.json`): every cell of
the committed `dxr_tail_residuals.csv` (4 conditions x 2 laws x 3 film
multipliers), every simulated row of `dxr_traces.csv` and every row of
`dxr_coefficient_profiles.csv`, re-marched through this item's driver with the
step recorder in place, string-equal to the committed files. The paper's
printed values are two of these rows: **0.7700807507969353** (soybean) and
**0.6967799577128457** (sunflower). Item 31's two reported-level decimated
traces, re-scored with the window rule used here, return item 31's
0.13297066864133203 (10 of 15) and 0.19708986980030843 (6 of 15) to the float.
The declared phase's 60-cell 0.885 mm borrowed-leg marches, run through the
same code path as the declared radii, return the committed 0.7700807507969353
and 0.6967799577128457 to the float.

## 2. The declared phase: 16 marches, committed stride rule

Window RMS over the window span, decimated "as printed" (the committed rule);
in band = measured points of 15 inside the reading band. A refused march is
not scored (its truncated trace would be held flat by the interpolation); its
end time is given from the start of drying, against the last scored sample at
300.1 s. Refusals are verbatim `RuntimeError: dry-shell Newton line search
failed inside the admissible storage domain`, after the committed seven
retries.

| leg | trace | sphere | 60 cells | 120 cells |
|---|---|---|---|---|
| borrowed | sunflower | 0.885 mm | **0.6968**, 0 of 15 (the paper's) | refused at 58.4 s |
| borrowed | sunflower | 0.7215 mm (own) | **0.8748**, 0 of 15 | refused at 59.6 s |
| borrowed | soybean | 0.885 mm | **0.7701**, 0 of 15 (the paper's) | refused at 52.2 s |
| borrowed | soybean | 0.6660 mm (own) | **refused at 88.7 s** | refused at 40.5 s |
| Whitaker | sunflower | 0.885 mm | 0.5642, 5 of 15 | refused at 81.9 s |
| Whitaker | sunflower | 0.7215 mm (own) | **0.7355**, 2 of 15 | refused at 82.8 s |
| Whitaker | soybean | 0.885 mm | 0.5384, 6 of 15 | refused at 85.6 s |
| Whitaker | soybean | 0.6660 mm (own) | **refused at 147.1 s** | refused at 70.0 s |

Beside them, the one-scalar boundary read (Whitaker leg, 60/1): item 31 at
the species-own spheres **0.133** (10 of 15, sunflower) and **0.197** (6 of 15,
soybean); item 27 at 0.885 mm 0.244 (2) and 0.366 (0).

Where a space-resolved march covers the window its mean signed residual is
negative (-0.37 to -0.82 of the span): the space-resolved law dries the
particle too FAST, the opposite side from the one-scalar boundary read, which
is too slow. The smaller sphere makes it faster still, so the residual rises
from 0.697 to 0.875 (borrowed) and 0.564 to 0.736 (Whitaker) on the sunflower
trace. In the window the particle volume above the law's measured loading
range (0.0156 kg/kg) falls from 0.95 to 1.0 at the window start to 0 to 0.42 at
its end; the within-particle coefficient ratio reaches 97 to 101.

**No 60/120 pair exists at any radius, leg or trace.** Every 120-cell march
refuses 40 to 86 s after the start of drying, 214 to 260 s before the last
scored sample, at the 0.885 mm sphere as well as at the declared radii. The
numerical uncertainty of the space-resolved residual therefore cannot be
stated from a mesh pair, and this includes the paper's 0.770 and 0.697.

## 2b. The declared extension: 32 marches at stride ceilings 5/4 s and 5/16 s

All 32 refuse before the end of the window, the 60-cell 0.885 mm borrowed-leg
configuration of the paper included: at the 5/4 s ceiling between 45.6 and
144.0 s after the start of drying, with the outer cell's loading 0.003 to 0.015
kg/kg and a within-particle coefficient ratio of 12 to 61 at the last accepted
stride; at the 5/16 s ceiling between 32.5 and 57.4 s, two to seven strides
after the crossing, with the coefficient ratio still 1 to 7. The refusal is
therefore not a property of the large stride: a smaller stride refuses
earlier. The 60-cell marches at the committed 5 s ceiling are the only ones of
48 that reach the horizon (6 of 48), and they do so only at the coarsest mesh
and stride of the item.

## 3. Rejecting outcomes

* **R1 FIRES.** At the declared radii the space-resolved march refuses inside
  the scored window on the soybean trace at 60 cells on both legs (88.7 s and
  147.1 s) and on both traces at 120 cells on both legs (six of eight
  declared-radius marches). The sunflower trace at its own sphere covers the
  window only at 60 cells: 0.875 (borrowed) and 0.736 (Whitaker). Nothing was
  relaxed.
* **R2 does not fire.** The committed space-resolved record is reproduced to
  the float (3,732 checks), and the paper's configuration through the
  declared phase's code path returns 0.7700807507969353 and 0.6967799577128457
  exactly.
* **Beyond the rejecting outcomes, measured:** the paper's own 0.885 mm
  space-resolved values have no mesh or stride partner that covers the window
  (16 of 16 refinements refuse, sections 2 and 2b).

## 4. Ledgers

Over all 48 marches (72 including the identity): every accepted stride closes
the particle inventory against the discrete Dirichlet-face outflow of the
implicit state to at most 9.6e-15 of the initial inventory; across the
constant-storage remap between strides to at most 7.7e-15; every whole march,
remaps included, to at most 6.7e-13. The frozen-coefficient fallback produced
one accepted stride in each of two marches; every other stride came
from the exact Newton.

## 5. The sentence the paper could carry

"At the declared volume-to-surface spheres the space-resolved solve is not a
usable comparator: on the sunflower trace it covers the window only at the
reported mesh and stride, giving 0.875 (borrowed correlation) and 0.736
(Whitaker) of the window span with the law drying the particle too fast, on
the soybean trace it refuses inside the window, and at the \SI{0.885}{\milli\metre}
sphere its 0.770 and 0.697 have no refined partner, every 120-cell or
smaller-stride march refusing before the window ends (Sec.~S10.7)."

## 6. Claim and non-claim

**Claimed:** the committed space-resolved record reproduces to the float; at
the declared radii the space-resolved march is refused inside the window on
the soybean trace and at every refinement, with the refusals and ledgers
recorded; where it covers the window it overshoots the drying (negative mean
residual) and the residual grows as the sphere shrinks; the paper's 0.770 and
0.697 exist only at 60 cells and the committed 5 s ceiling.

**Not claimed:** a converged space-resolved residual at any radius; that the
refusal is a defect of the law rather than of the variable-coefficient solver
at steep coefficient contrast (not diagnosed here; the module is read-only);
anything about the one-scalar numbers of items 27 to 31, which are unchanged;
any physical, plant, calibration or qualification statement. Nothing was
fitted, no tolerance, budget, floor, horizon, law or retry rule was changed,
and no paper file was edited.
