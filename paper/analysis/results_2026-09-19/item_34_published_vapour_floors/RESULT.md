# Item 34: the equilibrium residual floor at published vapour compositions

Answers discrepancy D1 of `paper/analysis/datasets/REDIGITIZATION_AUDIT_ARTICLES_2026-09-27.md`.
Read-only on everything outside this folder. Nothing fitted, nothing under `src/` touched.
The source's industrial plant is not named here.

## Verdict

1. **The 0.883/0.117 split is a published number with the wrong citation, the
   wrong location and the wrong temperature.** It is printed in Cardarelli (1998),
   thesis Figure 5.1 (p. 146), as the composition of the vapour leaving the top of
   the vessel at 75 °C (`y_H = 0.882905`, `y_W = 0.117095`, basis not printed). It
   is not in Cardarelli et al. (2002). The thesis records it as an external stream
   of an industrial plant. Read as a mass split, which is the reading that closes
   the printed balance, it is the vapour leaving above the pre-desolventizing
   stages. That vapour sits on the water dew line at 75 °C and is in contact with
   meal still carrying free solvent at 65 °C. item_06 (after the D4 settling
   measurement) paired it with the 105 °C sparged-section temperature. The dry meal
   never meets that gas.
2. **Published compositions for the gas the sparged meal actually meets exist.** The
   same thesis gives them, and every one floors below the 100 to 500 ppm class,
   with or without the oil arm. So do both Cardarelli models' vapour.
3. **The paper's "above the published class" scoping statement has no published
   basis.** The 1216.4 and 2085.5 ppm floors come only from the mis-paired
   composition.

## Provenance chain

| step | record | quote / content |
|---|---|---|
| source | Cardarelli 1998 thesis Fig. 5.1 p. 146; `paper/analysis/datasets/cardarelli1998_thesis_ch5/fig5_1_industrial_balance.csv` rows `vapores_salida` | "y_H_hexano, 0.882905, mol/mol or kg/kg (basis not printed), vapour outlet at the top"; temperature 75 C. DATASET_RECORD.md classes the external streams as "Plant data", taken from an industrial plant, published in the thesis |
| intake | `docs/GT_PS2_CARDARELLI1998_THESIS_CH5_INTAKE_2026-09-15.md` §2 | "vapours 30,770 kg/h at 75 °C, 88.3 percent hexane and 11.7 percent water by MASS" |
| basis | `docs/GT_PS2_FLASH_LOCATION_MEASUREMENT_RECORD_2026-09-15.md` l.16 | "that 0.117 is a mass fraction. The plant's exit vapour is therefore at 0.987 of water saturation at 75 °C" |
| first use | `docs/evidence/d4_settling_2026-09-15/s1_anchors.py` l.56-58 | "# the DT vapour composition of Figure 5.1 (88.3 % hexane, 11.7 % water by MASS)"; `D4_SETTLING.md` l.55: "At the plant's DT vapour composition the surface activity is 0.2219 (Figure 5.1's 88.3 percent hexane by mass gives mole fraction 0.612 ...)". The brief asked for "the plant's DT vapour composition"; the 75 °C exit stream was used |
| paper analysis | `item_06_cardarelli_tail/run_item06_tail.py` l.34-35 | "# The plant DT vapour composition used by the D4 settling measurement." evaluated at T_DT = 378.15 K |
| propagation | `item_23_oil_arm_activity/run_item23_oil_arm.py` l.74-82 | key `cardarelli2002_vessel_105C`, "the one vapour composition the industrial-vessel source reports: y = 0.612 n-hexane in steam at 105 C" (the citation drifts from the 1998 thesis to the 2002 article here) |
| paper | `sec_validation.tex` l.563-567, `supplementary.tex` l.1949, Table 2 row, `sec_discussion.tex` l.203, `sec_conclusions.tex` l.123 | "the one vapour composition \citet{cardarelli2002modeling} report for the industrial vessel, 0.612" |

The number is not a project assumption and not the owner's plant. It is a published
external-stream value of the thesis's plant, and it is used at a location and
temperature it does not describe.

## Table: floors at published compositions (ppm of dry meal; total storage, item_06 construction)

Arm is item_23's weight-basis oil arm at w_o = 0.0195 and Ω_w = 5.2. The bracket in
`item34_floors.csv` is 4.6 to 5.9. The class is 100 to 500 ppm.

| composition | source | T (°C) | y_hex | a_h | floor, no arm | floor, arm | vs class |
|---|---|---:|---:|---:|---:|---:|---|
| (a) diffusive-stage entry vapour, 0.03 kg/kg (text p. 149) | Cardarelli 1998 | 105 | 0.00642 | 0.00233 | **11.6** | **20.3** | below |
| same, at the gas below the PD stages | Cardarelli 1998 | 99 | 0.00642 | 0.00271 | 14.7 | 24.8 | below |
| (a) Fig. 5.3 vapour, 98.0 / 99.0 / 99.4 % flashed | Cardarelli 1998 | 105 | 0.0149 / 0.0069 / 0.0038 | 0.0054 / 0.0025 / 0.0014 | 26.8 / 12.4 / 6.8 | 47.1 / 21.7 / 12.0 | below |
| (a) sparged-section top vapour, this item's arithmetic on Fig. 5.1 + p. 149 (entry 3,500 / 2,850 ppm) | Cardarelli 1998 | 105 | 0.0091 / 0.0071 | 0.0033 / 0.0026 | 16.5 / 12.9 | 28.9 / 22.6 | below |
| (a) vessel-top exit vapour, Fig. 5.1, at its own 75 °C | Cardarelli 1998 | 75 | 0.6118 | 0.504 | 4834 | 6928* | above (not a dry-meal gas) |
| (a) same, at 105 °C (item_06's pairing; the paper's numbers) | — | 105 | 0.6118 | 0.222 | 1216 | 2085 | above (mis-paired) |
| (b) model vapour u, max of Fig. 2 (10.28 g/kg) / audit's ~11 g/kg | Cardarelli et al. 2002 | 105 | 0.00217 / 0.00232 | 0.00079 / 0.00084 | 3.9 / 4.2 | 6.9 / 7.3 | below |
| (b) thesis model vapour, max of Fig. 5.12 (10.16 g/kg) | Cardarelli 1998 | 105 | 0.00214 | 0.00078 | 3.9 | 6.8 | below |
| (c) water-saturated exit vapour, 0.068 to 0.47 kg water/kg hexane | Witte 1995 via Cardarelli 1998 Fig. 5.2 | 65 to 90 | 0.755 to 0.309 | 0.85 to 0.17 | 13,624 to 1,090 | 17,433 to 1,731* | above (PD exit gas) |
| (c) same at 60 °C (0.053 kg/kg) | Witte 1995 | 60 | 0.797 | 1.06 | refused: a_h ≥ 1, liquid hexane | — | — |
| (d) last trays, y 0.003 to 0.03 | paper's declared grid | 100 to 110 | 0.003 to 0.03 | 0.001 to 0.012 | 4.5 to 66.0 | 8.1 to 112.4 | below (arm: up to inside) |
| (d) last trays, y 0.1 | paper's declared grid | 100 to 110 | 0.1 | 0.032 to 0.041 | 152 to 222 | 273 to 378 | inside |

\* The oil-arm coefficient is outside its declared 100 to 120 °C band below 100 °C.
The 90 °C Witte reading is unresolved between its two digitizations.

## The sentence the paper can carry

In place of the vessel-gas sentence (sec_validation l.563-574 and every place that
carries 1216.4 and 2085.5 ppm):

> At the vapour composition Cardarelli (1998, p. 149) gives for the gas at the entry
> of the sparged section, 0.03 kg of n-hexane per kg of vapour (mole fraction 0.0064)
> at 105 °C, the same construction floors at 11.6 ppm without the oil-solution arm
> and 20.3 ppm with it, and across that source's Fig. 5.3 (0.004 to 0.015) at 6.8 to
> 47.1 ppm; at the model vapour of Cardarelli et al. (2002), at most about 10 g per
> kg (0.0022), it floors at 3.9 to 6.9 ppm. All lie below the published class, so
> equilibrium does not bind at any published composition of the gas the meal meets
> in the sparged section. The 88.3 per cent n-hexane printed for the vessel's exit
> vapour (Cardarelli 1998, Fig. 5.1) is the gas leaving above the
> pre-desolventizing stages at 75 °C and bears on no dry-meal floor.

The Table 2 row and the supplement (l.1946-1953) then carry "sparged-section gas:
6.8 to 47.1 ppm, below the class". The conclusions' "1216 ppm at the vapour
composition reported for the industrial vessel" is withdrawn. The supplement's
tail-sweep sentence keeps its hexane-free surface and loses the 0.612 floor. If the
lead prefers to cite only printed numbers, the minimal form uses p. 149 and Fig. 2
of 2002 only and omits this item's balance arithmetic (row a4), which is a
cross-check.

## Claim boundary

These are equilibrium statements of the frozen isotherm with the storage functional,
not predictions of discharge. The p. 149 and Fig. 5.3 compositions are the thesis
author's balance estimates. The 2002 and Fig. 5.12 compositions are model outputs.
None is a measurement of the gas. No plant number of this project enters.
