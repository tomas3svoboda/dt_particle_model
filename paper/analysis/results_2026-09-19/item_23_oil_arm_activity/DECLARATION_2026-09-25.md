# Declaration — the oil-arm activity coefficient, 2026-09-25

Written **before** the arm was re-run, as P-1 requires: criterion, value,
bracket and rejecting outcome first; the measurement after.  The run that
scores this declaration is `run_item23_oil_arm.py` in this folder and its
outcome is in `RESULT.md`.

This closes item **B1** of `paper/01_particle_jfpe/SUBMISSION_TODO_2026-09-24.md`
and mitigation **M6** / limit **17** of
`paper/01_particle_jfpe/LIMITS_AUDIT_2026-09-22.md`, both of which record that
no activity coefficient for n-hexane in vegetable oil was held.  Four sources
were supplied on 2026-09-25 and are digitized in
`paper/analysis/datasets/hexane_in_oil_activity/`.

## 1. The quantity being declared

The oil-solution arm of the residual floor states how much n-hexane the meal's
residual triglyceride oil holds at equilibrium with a pore gas at hexane
activity `a_h`.  Until now it was written on the **mole** fraction,

    x_h = a_h / gamma_h,   W_oil = w_o * (M_h / M_oil) * x_h / (1 - x_h),

with `gamma_h` **declared ideal** and bracketed over 1.0 to 1.5, and with
`M_oil` **computed from assumed triacylglycerols** (882.4 g/mol, bracket 807.3
to 885.5).  Two declared quantities, neither source-pinned.

It is re-declared on the **weight** fraction,

    w_h = a_h / Omega_w,   W_oil = w_o * w_h / (1 - w_h),

where `Omega_w` is the **weight-fraction activity coefficient of n-hexane in
vegetable oil at infinite dilution**, equivalently the weight-basis Henry
constant `H_w = Omega_w * p_sat(T)` of `p_hex = H_w * w_hex`.

**Why this form.** Three reasons, in order of weight.

1. **It is what the sources print.**  King & List (1990) tabulate `Omega_w`
   directly (their Table 2) and `H_w` directly (their Table 4), and their
   Eq. 2 obtains it from the retention volume, the *solute* molar mass and the
   *solute* vapour pressure — the **solvent molar mass does not enter**.  The
   same authors re-read Smith's static-cell data on the same weight basis
   ("the limiting slopes of solute vapor pressure–weight fraction
   relationship", p. 429), which again carries no oil molar mass.
2. **It retires a declaration instead of adding one.**  The oil molar mass
   leaves the arm entirely: one declared quantity becomes zero declared
   quantities and one source-pinned one.
3. **It is the conservative of the two linear readings at the arm's own
   activities.**  The arm is evaluated at `a_h` of 0.18 to 0.41, which is
   above the concentration limit below which the sources state Henry's law
   applies.  The weight-basis and mole-basis linear extrapolations of the same
   measured infinite-dilution coefficient diverge there, and the weight-basis
   one gives the smaller arm.  That is stated as a *reason for the choice*
   and is **not** a licence to report the smaller number alone: the run
   reports both, and the mole-basis reading is scored against the measured
   vapour–liquid equilibrium the sources also publish.

## 2. The criterion

**The four-method spread at the nearest temperatures to 100 to 120 °C.**

Four independent determinations of `Omega_w` for n-hexane in a refined
vegetable oil are held, each reduced to the weight basis using only quantities
its own source prints:

| # | source | method | temperatures measured |
|---|---|---|---|
| 1 | Smith & Wechter 1950 | static vapour pressure over hexane–soybean oil | 75–120 °C |
| 2 | King & List 1990 | inverse gas chromatography, soybean oil stationary phase | 58.7–123.4 °C |
| 3 | Belting et al. 2014 | dilutor (inert gas stripping), refined soybean oil | 40–80 °C |
| 4 | Belting et al. 2015 | computer-driven static cell VLE + the authors' UNIQUAC correlation of it | 75 and 100 °C |

Each is evaluated at 100, 105, 110, 115 and 120 °C, by its own source's own
temperature form where the source's range does not reach (the van't Hoff line
the source itself regresses for methods 2 and 3; the UNIQUAC correlation the
source itself fits for method 4).  The declared value is the mean of the four
over that band; the declared bracket is the envelope of the four over that
band, rounded outward to two significant figures.

## 3. The declared value and bracket

| item | value |
|---|---|
| **declared quantity** | `Omega_w`, the weight-fraction activity coefficient of n-hexane in refined vegetable oil at infinite dilution, dimensionless |
| **declared value** | **5.2** |
| **declared bracket** | **4.6 to 5.9** |
| **equivalent Henry constant** | `H_w = Omega_w * p_sat(T)`; at 105 °C this is **14.4 atm = 1.46 MPa** at the declared value |
| **equivalent mole-basis coefficient** | `gamma_x = Omega_w * M_h / M_oil`; **0.51** at the declared value with the sources' own oil molar masses (870 to 876 g/mol) |
| **declared domain** | 100 to 120 °C; refined soybean, sunflower or rapeseed oil; the Henry region, which the sources place below `a_h` of about 0.10 (120 °C) to 0.16 (100 °C) |
| **temperature dependence** | none declared. The coefficient is carried constant over 100–120 °C, because the sources place the heat of mixing at 0.14 ± 0.06 kcal/mol and 1.57 ± 0.31 kJ/mol, which moves `Omega_w` by less than the four-method spread across the band |
| **oil molar mass** | **no longer declared**; it leaves the construction |

## 4. The rejecting outcomes, written first

**R1 — the bracket.**  Any of the four methods' weight-basis coefficient, at
any temperature in 100 to 120 °C, falls outside 4.6 to 5.9.  Then the bracket
does not cover the evidence and must be widened before it is used.

**R2 — the Henry cross-check.**  The arm's own `p` versus `w` relation at the
declared value, expressed as `H_w` at 100.9 and 102.2 °C, differs from the two
printed Henry constants at those temperatures — 14.2 atm (King & List's own
IGC value, their Table 4) and 15.5 atm (their re-reading of Smith's static
data, their p. 429) — by **more than the four-method spread** (22.7 per cent
full width about the mean).  Then the pinned coefficient is not consistent
with the independent Henry constants and the pin fails.

**R3 — the arm's role in the paper.**  The arm at the declared value, added to
the arm-free floor, **exceeds the measured final loading** of
7.415e-3 kg/kg at either journal condition.  Then the arm is no longer the
upper bound of a bracket whose lower bound is the arm-free floor, the paper's
framing of it is falsified, and the framing must be withdrawn rather than the
number adjusted.  (This is the rejecting outcome
`LIMITS_AUDIT_2026-09-22.md` M6 already carries, restated for a coefficient
that is now pinned from below as well as above.)

**R4 — the overlap.**  The arm at the declared value, measured against the
frozen total-meal retention isotherm **at that isotherm's own measured
conditions**, exceeds 100 per cent of the measured retention.  Then the arm
could be duplicating the whole of the measured isotherm rather than part of
it, and the arm may not be added at all until the overlap is measured.

## 4a. One correction made to this declaration before the run

The bracket was first written **4.6 to 5.8**, which contradicts this
declaration's own rule in section 2: the envelope's upper end over the band is
5.804, and *rounded outward* to two significant figures that is **5.9**, not
5.8.  The rule is unchanged and the correction was made before
`run_item23_oil_arm.py` was scored; it is recorded here rather than silently
applied, because a bracket that fails its own rule would have fired R1 on a
rounding artifact rather than on evidence.

## 5. What this declaration does NOT claim

- It does not make the arm qualifying. The arm stays analysis-side and
  `physically_qualifying = False`; nothing here is promoted into
  `src/dtdc_simulator`.
- It does not claim the sources measured the meal's residual oil.  They
  measured refined bulk vegetable oils.  The residual oil of an extracted meal
  is not a refined oil, and no source for its activity coefficient exists.
- It does not extend the Henry region.  The two journal conditions lie above
  the limit the sources state, and the arm there is an extrapolation of a
  measured coefficient, not a measurement.
- It does not touch the tower model's own oil closure
  (`src/dtdc_simulator/core2/props/sorption.py`, PHY-048).  That comparison is
  reported for the lead and is a separate Class-B question.
- No parameter is fitted.  The coefficient is source-pinned; the one
  regression performed anywhere in this item is the reproduction of each
  source's *own* published temperature regression, checked against that
  source's own printed result.
