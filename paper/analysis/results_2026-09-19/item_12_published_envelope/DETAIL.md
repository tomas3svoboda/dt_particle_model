# Item 12 — industrial validation against the published envelope

Fills **declared missing result 12** (`sec_missing.tex` item 12) within the
owner's ruling of 2026-09-03: plant data are not a benchmark, the published
operating envelope is. The particle model is scored here against every band of
`docs/GT_PS2_PUBLISHED_OPERATING_ENVELOPE_2026-09-04.md` that a single-particle
model can reach, and the bands it cannot reach are named rather than left
unmentioned. No calibration was performed and none is claimed.

## Method

`run_item12_envelope.py` (wall time 2.6 s) uses the frozen isotherm
(`props/sorption.py`), the IAPWS-95 water surface and the Span n-hexane
surface. Nothing is fitted. Band endpoints are quoted with their envelope row
identifiers.

## A. The sorption-elevated boiling point, over the published moisture band

Solving a_w(W) P_sat^w(T*) = P at the published dome pressure P4 = 0.993 bar,
over the published DT-outlet moisture band M4:

| moisture (wt % wb) | W (kg/kg dry) | a_w | T* (°C) | inside T10 = 104.6–106.3 °C |
|---|---|---|---|---|
| 16.0 | 0.19048 | 0.7430 | 107.922 | no (above) |
| 17.0 | 0.20482 | 0.7630 | 107.143 | no (above) |
| 18.0 | 0.21951 | 0.7812 | 106.454 | no (above, by 0.15 K) |
| 19.0 | 0.23457 | 0.7978 | **105.840** | **yes** |
| 19.072 (the model's cap) | 0.23564 | 0.7990 | 105.799 | yes |
| 21.0 | 0.26582 | — | refused | `retained-water loading 0.26582278481012656 exceeds the evidence cap 0.235670904; apply the explicit active set, do not clamp` |

The 19.0 wt % row reproduces the paper's 105.84 °C exactly. **The prediction
enters the published T10 band at 18.240 wt % wb and stays inside it up to the
model's own evidence cap at 19.072 wt %.** Over the whole admissible part of
M4 (16.0 to 19.072 wt %) the prediction spans 105.799 to 107.922 °C, so it
brackets the measured band from above: the published moisture band is wider
than the temperature band it maps onto, and only its upper 0.83 points of
moisture are consistent with the published discharge temperatures. Pure water
at the same pressure boils at 99.410 °C, so the sorption elevation is 6.4 to
8.5 K across the admissible band. At 1 atm the same sweep gives 106.386 to
108.517 °C.

**Score: reached, and the model lands inside the band on the upper 27 % of the
admissible moisture range.**

## B. The pure n-hexane anchor

Inverting the frozen saturation curve: **68.7145 °C at 1 atm**, 68.073 °C at
P4. This sits inside the published dome design band V1 (66–78 °C, typical 71)
and above the heteroazeotrope anchor V5 (60–62 °C). The dome band is a tower
vapour-outlet temperature, not a particle observable, so it is reported for
placement only; the particle observable this anchor actually scores is the
constant-rate solid temperature of item 5(c).

## C. Admissibility of the published moisture bands

The qualified modified-Luikov branch runs from W_ref = 0.023 to
W_cap = 0.2356709 kg/kg dry, that is **2.248 to 19.072 wt % wet basis**.

| envelope row | published band (wt % wb) | admissible sub-band | fraction admissible |
|---|---|---|---|
| M2 feed water at the DT inlet | 5.0–10.0 | 5.0–10.0 | 100 % |
| M3 leaving the first countercurrent tray | 17.0–22.0 | 17.0–19.072 | 41.4 % |
| M4 leaving the DT | 16.0–21.0 | 16.0–19.072 | 61.4 % |
| M5 product after the DC | 12.0–13.5 | 12.0–13.5 | 100 % |

**Score: the two bands that matter for the condensation step are only partly
reachable.** The upper 2.93 points of M3 and the upper 1.93 points of M4 lie
above the isotherm's evidence cap and are refused, not extrapolated.

## D. Residual hexane, M6

Published specification class 100–500 ppm. At the plant DT vapour composition
the certified pore-vapour arm settles at **1216.4 ppm** for every diffusivity
in the declared interval (item 6): **2.43 times** the upper end and **12.16
times** the lower end of the band. **Score: not reachable by the pore-vapour
arm at any coefficient.**

## E. Temperature placement

Declared ranges: 323.15 K minimum, 343.15 K direct-evidence maximum, 433.0 K
sensitivity maximum. Of the nine published temperature bands, **none sits
entirely inside the direct-evidence range**: T1 (49–60 °C) and T14 (30–108 °C)
fall partly *below* 323.15 K, and T3, T5, T6, T7, T9, T10 and T11 all reach
above 343.15 K into the sensitivity extension. Every band lies inside the
sensitivity range. **Score: reachable only as declared extrapolation.**

## F. What a single-particle model cannot reach

Named, not skipped: Q1–Q17 (tower scale), S1–S15 (steam balances), J1–J6
(jacket and deck), P1–P8 (bed and sparge hydraulics), R1–R7 (mechanical; R7 is
a consumed property, not a predicted band), V3–V4 (dome outlet streams),
M1 and M7 (an input, and a tray-residence quantity), M8 (protein quality), and
T2, T4, T7, T12, T13 (deck-surface, bed-profile, pilot and compartmented-unit
readings). Full register in `item12_envelope.json`, section F.

## Claim boundary

Three of the envelope's fifty-odd rows are reachable by a particle model and
are scored here; one of those three (M6) the model fails, with the number. The
boiling-point comparison takes measured moisture as an *input* and predicts a
different observable, exactly as `sec_validation.tex` already states; widening
it over the published band does not make it a transient-model validation. The
isotherm is evaluated about 36 K above its data. No calibration was performed.
Not physically qualifying, not plant predictive.
