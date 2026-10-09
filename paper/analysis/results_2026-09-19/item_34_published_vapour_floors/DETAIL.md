# Item 34 detail

## Construction (unchanged reuse)

- Floor: `c_g = y P M_h / (8.31451 T)` and `ppm = d4common.ppm_from_cg(c_g, T, DryShellParams())`,
  exactly `item_06_cardarelli_tail/run_item06_tail.py` l.35-36 and section C. This is the total
  (pore gas plus sorbed) storage of the certified dry shell at w_o = 0.0195.
- Arm: `oil_activity_coefficient.corrected_floor_weight_basis(floor, min(1, yP/Psat), 0.0195)`
  from item_23, Ω_w = 5.2 (bracket 4.6 to 5.9). An extra column in the CSV repeats the arm at the
  source's own oil where the source gives one. For Cardarelli 1998 that is Fig. 5.1 product oil over
  solid, 0.005492/0.802499 = 0.00684. For Cardarelli et al. 2002 it is Fig. 5's 0.015.
- Reproductions asserted before use:
  - item_06: 1216.358141953 ppm at y = 0.6120581 and 105 °C, equal to `item06_tail.json`.
  - item_23: 2085.51 ppm with the arm.
  - item_05's last-tray grid: every y < 1 row to 5.7e-14 ppm.
- The activity from `ds.activity(c_g, T)` equals `yP/Psat` to 2.2e-16 relative in every row.
- Refusal: when `yP/Psat ≥ 1` the row is refused, because liquid hexane would condense and no
  dry-shell floor exists. This fires for Witte at 60 °C. Nothing is clamped.

## Compositions and arithmetic

**(a1) Fig. 5.1 exit vapour.**
- 0.882905/86.177 = 1.02452e-2 mol/kg and 0.117095/18.015 = 6.49986e-3 mol/kg, so y = 0.61184.
- The paper's 0.61206 uses the rounded 0.883/0.117.
- a_h at 75 °C = 0.504. The mole-basis reading (y = 0.883, a_h = 0.728) is carried as a basis
  sensitivity. The flash-location record rejects that reading because the balance closes on
  mass, and at 75 °C the mass reading lands on Witte's saturation curve (row c at 75 °C,
  y = 0.623).

**(a2) p. 149 text.** "al ingreso de la etapa difusiva ... 3500-4000 ppm para el solido y 0.03 kg
de hexano/kg vapor para la fase vapor", with 0.97 kg water per kg vapour. This is the author's
balance estimate (`text_numbers.csv`).
- (0.03/86.177)/(0.03/86.177 + 0.97/18.015) = 3.4812e-4/(3.4812e-4 + 5.38440e-2) = 0.006424.
- It pairs with 105 °C, the desolventizing temperature of Table 5.2 and the sparged-tray gas of
  Fig. 5.1. At 99 °C, the gas below the last pre-desolventizing stage, it is also reported.

**(a3) Fig. 5.3 vapour panels (p. 150).** Digitized hexane and water kg per kg vapour at 98.0,
99.0 and 99.4 % flashed. The 99 % point (0.0320, 0.9677) reproduces the p. 149 text.

**(a4) This item's balance arithmetic (Fig. 5.1 external streams + p. 149 water loading).**
- Dry solid: 95,235 × 0.61374 = 58,449.5 kg/h.
- Oil: 95,235 × 0.0042 = 400.0 kg/h.
- Water at the diffusive entry: 0.1775 × 58,449.5 = 10,374.8 kg/h.
- Product water: 72,838 × 0.191543 = 13,951.6 kg/h, so condensed in the sparged section =
  3,576.8 kg/h.
- Direct steam: 8,200 + 70 = 8,270 kg/h. Water vapour leaving the sparged section top =
  8,270 − 3,576.8 = 4,693.2 kg/h.
- Hexane in the meal at entry: 69,224.3 × c/(1 − c) at c = 3,500 or 2,850 ppm (wet-meal basis,
  as the product's 500 ppm = x_H 0.0005), which is 243.1 or 197.9 kg/h. Product hexane is
  36.4 kg/h, so 206.7 or 161.4 kg/h is stripped.
- Vapour hexane mass fraction 0.0422 or 0.0333, giving y = 0.00912 or 0.00714.
- Assumptions:
  - the side steam enters the sparged section;
  - the stage numerals are wet-meal ppm (the dataset's reading, not a printed label);
  - the p. 149 water loading holds at the section entry.
- The result brackets the author's own 0.03 kg/kg (text) and is a cross-check, not a source value.

**(b) Model vapour.**
- Cardarelli et al. 2002 Fig. 2 u, maximum 10.276 g solvent per kg vapour at z* = 1 (both Pe
  series), from `paper/analysis/datasets/cardarelli2002/curves.csv`. The audit's "about 11" is
  also run.
- Cardarelli 1998 Fig. 5.12, maximum 10.157 g/kg.
- u falls to about 0.2 g/kg at the bed bottom, so the floor there is lower still.
- Both are at 105 °C (Table 1 / Table 5.2).

**(c) Witte 1995 (Fig. 5.2, p. 147).** Digitized 0.0533 to 0.4682 kg water per kg hexane over
60 to 90 °C, giving y = 1/(1 + r·86.177/18.015). This is the water-saturated exit vapour of the
pre-desolventizing section, whose meal carries free solvent. It is not a dry-shell situation, and
it is reported only because the brief asked. At 60 °C, yP > Psat,hex, so the row is refused. The
90 °C point is flagged `b_outside_envelope` in the dataset.

**(d) Last trays.** The paper's declared generic grid, y 0.003, 0.01, 0.03 and 0.1 at 100, 105
and 110 °C, reproduces the printed 4.5 to 66.0 ppm and 152 to 222 ppm. With the arm the same
rows are 8.1 to 112.4 and 273 to 378 ppm. At y = 0.03 and 100 °C the arm lifts the floor to
112 ppm, just inside the class; the paper's "below the class" for y 0.003 to 0.03 holds for the
arm-free floor only.

## Domains

The isotherm's measured band is 50 to 95 °C, so every 99 to 110 °C row extrapolates it, as all
the paper's 105 °C numbers already do. The arm coefficient's declared band is 100 to 120 °C, so
the 60 to 99 °C rows extrapolate it. Every row at or above 100 °C is in the arm band, and every
row below 100 °C is flagged in the CSV.

## Files

- `run_item34_floors.py`: the script (about 1 s).
- `run_stdout.txt`: its output.
- `item34_floors.csv` and `item34_floors.json`: all rows, the reproductions and the (a4)
  arithmetic.
- `MANIFEST.sha256`.
