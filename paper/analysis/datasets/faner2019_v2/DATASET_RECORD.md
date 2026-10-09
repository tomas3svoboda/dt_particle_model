# Dataset record: faner2019_v2/Figures 2-4 and conditions (version 2, 2026-09-27)

- **Source:** Faner, S. A., Perez, E. E., & Crapiste, G. H. (2019), "Desolventizing kinetics of oilseed meals with superheated hexane," *Journal of Food Process Engineering* 42:e12987, DOI 10.1111/jfpe.12987; bibliography key `faner2019kinetics`.
- **Source file:** `literature_sources/Faner_Desolventizing_Kinetics_of_Oilseed_Meals_with_Superheated_Hexane.pdf`
- **Source sha256:** `37c57202c4d5562f08becbfb77b6822f909cf82df0ec37012aa234d2fe86fb55` (verified before extraction).
- **Extracted objects:** Figures 2, 3 and 4, printed page 5 (PDF page 5); conditions from Sections 2.1, 2.2 and 3 (pages 2-6).
- **Quantities and units:** as version 1. Figure 2: solvent loading X, kg hexane/kg dry meal; particle temperature Tp and gas temperature Tg, °C; versus time, s. Figures 3 and 4: normalized loading (X-Xe)/(Xo-Xe), dimensionless, versus time, s.
- **Conditions:** pure superheated n-hexane, no water or steam. Figure 2 is sunflower at 100 °C and soybean at 120 °C. The temperature labels of Figures 3 and 4 are kept per series. See `conditions.csv`.
- **Supersedes:** `../faner2019/` (version 1, 2026-08-22). Version 1 is not edited. Its folder carries only the pointer `SUPERSEDED_BY_V2_2026-09-27.md`.
- **Claim boundary:** as version 1. This is digitized evidence from thin-layer batch runs in pure superheated hexane. It is not binary-film, steam, water-condensation, industrial-bed, calibration or production-qualification evidence. No march was re-run to build it. No result item, paper file or `docs/` file was changed.

## 1. Why version 2

The independent audit `../REDIGITIZATION_AUDIT_ARTICLES_2026-09-27.md` (§3.1, §3.2, §4 items D2-D5, D8, D10-D12) reported these problems in version 1:

- **Figure 4, 92 °C.** The committed points from 140 to 200 s lie on the model line, not on the open squares, and read 0.012 to 0.019 high.
- **Figure 4, 109 °C and 120 °C.** The two triangle series are swapped or duplicated, with errors up to 0.028.
- **Figure 2, sunflower loading.** After about 175 s the sunflower markers merge with the soybean markers. The final sunflower loading was set equal to the soybean value, 7.415e-3, and cannot be read. Its declared ±0.003 is ±40 %.
- **Figure 2, gas temperatures.** These read about 0.5 °C high because of the triangle anchor convention.
- **Conditions.** Several are misstated or overstated: the 2 s mass-log interval, the sphericity provenance and the dry basis of the initial ratio.

This version adds a third reading that is independent of both of those readings. It adjudicates the three readings point by point, and it transcribes the conditions again word for word.

## 2. Resolution and axis calibration (third reading)

The three figures are the page-5 DCT image streams. The third reader extracted them byte for byte with pypdf 6.16.1; they are in `third_reading/source_images/`.

| Figure | Pixels | Colour | Placement on the page | Resolution |
|---|---|---|---|---|
| 2 | 993 × 672 | 8-bit grey | 238.28 × 161.18 pt | 300.0 ppi |
| 3 | 996 × 733 | 8-bit grey | 239.13 × 175.86 pt | 299.9 ppi |
| 4 | 993 × 724 | 8-bit grey | 238.23 × 173.65 pt | 300.1 ppi |

The version-1 files `figure*_source_embedded_300ppi.jpg` were re-encoded. They differ from the PDF streams by up to 6 grey levels. They were not used here.

One pixel corresponds to:

- **Figure 2:** 0.309 s, 1.236e-3 kg/kg and 0.251 °C;
- **Figure 3:** 0.279 s and 1.615e-3;
- **Figure 4:** 0.283 s and 1.637e-3.

**Calibration.** Every printed tick was located as the darkness-weighted centroid of the tick stroke outside the frame, and a straight line was fitted by least squares. Residuals are in px (tick minus fit). Source: `third_reading/calib.json` and `third_reading/scripts/calib.py`, `calfit.py`.

| Axis | Ticks used (value @ px, residual) | RMS / max residual, px |
|---|---|---|
| Fig. 2 time | 0@106.03 (+0.05); 30@203.15 (+0.00); 60@300.26 (-0.03); 90@397.35 (-0.05); 120@494.45 (-0.07); 150@591.52 (-0.07); 180@688.37 (+0.15); 210@785.45 (+0.14); 240@882.78 (-0.13) | 0.090 / 0.151 |
| Fig. 2 X | 0.6@11.32 (-0.03); 0.4@173.12 (+0.03); 0.3@254.14 (-0.06); 0.2@334.89 (+0.11); 0.1@415.91 (+0.02); 0@496.93 (-0.07) | 0.063 / 0.115 |
| Fig. 2 T | 140@11.52 (-0.03); 120@91.29 (-0.05); 100@170.99 (+0.01); 80@250.69 (+0.06); 60@330.46 (+0.04); 40@410.17 (+0.07); 0@569.84 (-0.09) | 0.056 / 0.093 |
| Fig. 3 time | 0@111.07 (+0.03); 30@218.57 (-0.02); 60@326.07 (-0.06); 90@433.58 (-0.11); 120@540.78 (+0.14); 150@648.31 (+0.07); 180@755.79 (+0.04); 210@863.31 (-0.03); 240@970.80 (-0.06) | 0.072 / 0.137 |
| Fig. 3 ordinate | 1@12.81 (-0.20); 0.8@136.25 (+0.16); 0.6@259.96 (+0.26); 0.4@384.20 (-0.17); 0.2@507.88 (-0.05); 0@631.62 (+0.01) | 0.165 / 0.255 |
| Fig. 4 time | 0@118.95 (-0.05); 30@224.93 (+0.11); 60@331.18 (-0.01); 90@437.42 (-0.11); 120@543.40 (+0.04); 150@649.62 (-0.04); 180@755.59 (+0.12); 210@861.90 (-0.06) | 0.078 / 0.117 |
| Fig. 4 ordinate | 1@12.67 (-0.10); 0.8@134.67 (+0.08); 0.6@256.98 (-0.04); 0.4@378.98 (+0.14); 0.2@501.29 (+0.02); 0@623.58 (-0.10) | 0.087 / 0.137 |

Four ticks were excluded because a data marker overlaps the tick:

- Fig. 2 X = 0.5, by the t = 0 markers;
- Fig. 2 T = 20 °C, by the 240 s circle;
- Fig. 4 t = 240 s, by the 240 s circle;
- Fig. 4 t = 210 s, by a marker. This tick was centred by hand from the pixel column (861.9), not by the centroid routine.

The value-1 ticks of Figures 3 and 4 coincide with the top frame line, so the centroid of that line was used for them.

All calibrations agree with version 1 and with the auditor to within 0.5 px.

## 3. Series identity (re-derived from the legend, the caption and the text)

**Figure 2.**

- The legend reads, from top to bottom: "X - Sunflower" (×), "X - Soybean" (open circle), "Tp - Sunflower" (filled diamond), "Tg - Sunflower" (filled triangle), "Tp - Soybean" (open diamond) and "Tg - Soybean" (open triangle).
- The caption reads: "Experimental desolventization of oilseed meals with superheated hexane al 100 °C (sunflower) and 120 °C (soybean)."
- The text on p. 5 reads: "at 100 °C for sunflower and 120 °C for soybean".
- The soybean circles are drawn over the sunflower crosses. Their white interiors hide the crosses, which is the legend order.
- The Tp and Tg markers are spaced 15 s apart: sunflower at about 0, 19, 34 ... 229 s and soybean at about 0, 10, 25 ... 235 s. The text says the temperatures were logged "every 10 s".

**Figure 3 (sunflower; "Figure 3 shows the effect of temperature on the desolventization of sunflower meals").**

- The legend reads: open square 93 °C, grey square 101 °C, open triangle 110 °C, grey triangle 121 °C, open circle "125°" (printed without C) and grey circle 136 °C, plus "Model" as a solid line.

**Figure 4 (soybean; "A similar trend was observed for soybean meals (Figure 4)").**

- The legend reads: open square 92 °C, grey square 101 °C, open triangle 109 °C, grey triangle 120 °C, open circle "125°" and grey circle 136 °C, plus "Model".
- The fill is the only attribute that tells 109 °C from 120 °C: open against grey. The source gives no other key.
- The drawing order follows the legend order. A later series covers an earlier one, and the grey circles lie on top of everything.

### 3.1 Figure 4: resolving the 109 °C and 120 °C series

**20 s.** The fill was checked pixel by pixel:

- The grey triangle (120 °C) has its apex at row 151.3 and its base at row 167-168. Its interior darkness is about 100.
- The open triangle (109 °C) has its base at row 159. The grey triangle interrupts that base line, so the open triangle lies underneath. Its interior is white (darkness 0 to 4), and that white interior covers the fill of the grey square beneath it.
- The open triangle is therefore the higher one: 109 °C = 0.769 and 120 °C = 0.755.

**30 to 60 s.** Both triangles are visible, with the open one above the grey one:

| t (s) | 109 °C | 120 °C |
|---|---|---|
| 30 | 0.654 | 0.637 |
| 40 | 0.541 | 0.530 |
| 50 | 0.441 | 0.415 |
| 60 | 0.330 | 0.319 |

**70 s.** Only the grey triangle is visible. The base of the open triangle shows 2-3 px below the grey base, which gives 109 °C = 0.249 (merged).

**80 to 110 s.** The two series cross. The grey triangle now lies above the open one:

| t (s) | 120 °C | 109 °C |
|---|---|---|
| 80 | 0.208 | 0.183 |
| 90 | 0.174 | 0.160 |
| 100 | 0.150 | 0.127 |
| 110 | 0.120 | 0.105 |

The text states that flow rate, meal amount and initial content differed between runs. The source therefore does not require the series to be ordered by temperature, and version 2 does not impose that order.

**120 s.** The base line of the open triangle shows 3 px below the grey triangle, giving 0.089 (merged).

**Other times.**

- **After 120 s** no open triangle can be seen anywhere. It is covered by the grey triangle, the grey square and the circles. The 109 °C points from 130 to 240 s are flagged `inferred`, and none of the three readers could read them. Version 1's values there track the 109 °C model line to within -0.008 to +0.012.
- **At 10 s** the grey square, the grey triangle and a hidden open triangle are stacked at 0.88-0.89. Both triangle points are `merged` or `inferred`.

**What version 1 did.** In version 1 the 109 °C and 120 °C values are:

- swapped at 20 s and at 40 s;
- duplicated at 50 s (both 0.4444, the value of the open triangle) and at 60 s;
- at 80 s, the 109 °C point carries the 120 °C value;
- at 130-150 s, the 120 °C points lie 0.011-0.019 below the grey triangles.

The auditor and the third reader agree on the identity at every time where either triangle is visible.

**What the source supports.** It supports the series identity given in the legend and resolved above for 20-60 s and 80-110 s. At 10 s, 70 s and 120 s it supports only merged values. It supports no reading of 109 °C after 120 s.

### 3.2 Figure 2: separating the sunflower and soybean loading markers

The third reader looked at the pixel columns from 150 to 240 s (`third_reading/scripts/dumpmulti.py`). The soybean circles form a chain drawn on top. The solid top of that chain lies at rows 477-484 and its solid bottom at rows 487-494. The circle radius is 5.75 px and the cross half-height is 5 px.

| t (s) | sunflower cross centre (row) | soybean circle top / centre (rows) | state |
|---|---|---|---|
| 150-170 | 456 / 464.5 / 471 | 472.5 / 482 | separate bands |
| 180 | 474 | 481 / 486 | the lower arms of the crosses touch the circle tops |
| 190 | 480 | 483 / 488 | centres visible 3 px above the circle tops; lower arms hidden |
| 200-220 | 483-486, from the upper tips only | 483 / 488 | the cross centres lie under the circle tops |
| 230-236 | about 487-488, from tips at most 1 px high | 483.5-484 / 489 | inseparable |

- **When the markers become inseparable.** From about 195-200 s the centre of each sunflower cross is hidden under the soybean circles. From then on only the upper tips of the crosses can be seen. By 230-236 s those tips stand no more than 1 px above the circle tops, and the two traces cannot be told apart. The auditor's "merge after about 175 s" refers to the touching arms, which the third reading also sees at 170-180 s.
- **Last separable sunflower value:** 0.0208 kg/kg at 190 s, uncertainty ±0.003 (flag near-overlap). The version-2 row at 190 s is the adjudicated 0.0222 ± 0.006: all three readings (v1 0.0222, auditor 0.023, third 0.0208) agree within the widened tolerance. The last value with both arms of the cross visible is 0.028 ± 0.003 at 180 s.
- **Final sunflower loading at 240 s:** 0.0105 kg/kg, flag `inferred`. The plausible range is 0.0075-0.0125, carried as ±0.003, which is ±29 % and one-sided in practice.
  - The upper bound comes from the height of the cross tips at 233-236 s.
  - The lower bound is the soybean value, because below it the crosses would be completely hidden.
  - The three readings are v1 0.007415 (set equal to soybean), the auditor "about 0.009-0.013 or less" (midpoint 0.011 used) and the third reading 0.0105. The median is taken.
- **Final soybean loading at 240 s:** 0.0080 ± 0.006 (flag merged; the median of v1 0.007415, auditor 0.008 and third 0.0092). The marker at 240 s is hidden by the right axis line. The circles at 233-238 s read 0.0088-0.0094.

## 4. Method of the third reading

**Scope.** The third reader read every experimental marker and sampled every model line. The scripts in `third_reading/scripts/` are kept exactly as they were run. They were run from the session scratch directory, so their absolute paths point there.

**Figures 3 and 4, experimental markers.**

- Eye seeds were taken from 3x-8x crops with a grid in data units and from pixel dumps (`gcrop.py`, `dump.py`). They are in `third_reading/seeds_f3.txt` and `seeds_f4.txt`.
- Each seed was refined by normalized cross-correlation in a ±3 px window, using that series' own legend symbol as the template (`tmpl.py`, `ncc.py`, `refine.py`).
- The refinement was accepted only when the score was at least 0.5 and the peak was not on the edge of the window. Otherwise the seed was kept.
- The result is in `mine_f3.json` and `mine_f4.json`, with flags `c` (clear), `n` (near-overlap), `m` (merged), `i` (inferred) and `u` (unreadable).

**Triangle anchor.** In the Figure 2 legend the connecting line crosses the triangles 1.9-2.0 px below the centre of their bounding box. That is 0.147 of the anti-aliased box height. The data lines do the same. Triangles are therefore anchored at that height:

- Fig. 3 110 °C: -0.0036, and 121 °C: -0.0038;
- Fig. 4 109 °C: -0.0043, and 120 °C: -0.0041;
- Fig. 2 Tg: -0.49 °C.

Version 1 used the box centre. Before adjudication its triangle values were converted to the line-height anchor (column `v1_value_convention_corrected`). The auditor already used the line-height anchor.

**Figure 2, loading.**

- **Up to 160 s (sunflower) and 150 s (soybean):** the midpoint of the dark band of each marker chain in the column at each grid time (`f2load.py`). This cancels the slope bias of the overlapping chain to first order.
- **Later times:** read from the pixel dumps as described in §3.2 (`mine_f2.py`).
- **t = 0:** read from the part of the marker still visible on the axis (flag merged).

**Figure 2, temperatures.** Every marker was found by cross-correlation with its own legend symbol, including the connecting line (`f2temp.py`, `f2temp2.py`). The anchor is the legend line height for triangles and the symbol centre for diamonds.

**Model lines.** Stroke centres were taken in each column at the version-1 abscissae (15, 25 ... s), between the symbol columns. Two touching strokes were split at the middle (`lines.py`). Each value was assigned to the version-1 series whose value was nearest within 0.012.

**Blindness, stated honestly.**

- Before reading, the third reader had seen the audit record: its findings, calibrations and conditions, but not its CSV values. The reader had also seen the first four rows of the version-1 `curves.csv`.
- The Figure 3 and 4 seeds were written down before the version-1 and auditor values for those figures were opened. The two were compared afterwards with `compare.py`.
- For Figure 2 the version-1 loading table and the auditor's Figure 2 CSV were opened before the loading was read, because the version-1 time grid was needed. The Figure 2 loading reading is therefore not blind.
- The Figure 2 temperatures and all model lines were measured by algorithm, without seeds from the other readers.

**Agreement of the third reading with the auditor.** On the 110 Figure 3 and 4 symbols that both read and that are flagged clear, the two agree with a median absolute difference of 0.0002 and a maximum of 0.0014. Both used the same PDF pixels and legend-symbol templates. The two readings are not fully independent in method. They are independent in their seeds and in how they judged overlaps.

## 5. Adjudication rule (`third_reading/scripts/adjudicate.py`, run by `build_v2.py`)

Each point has up to three readings: version 1 (anchor-corrected), the auditor and the third reading.

**Which auditor value is used.**

- Figures 3 and 4: the template value, or the eye value where the auditor flagged a low template score.
- Figure 2 loading: the marker centre when it lies within 0.6 s of the grid time, otherwise the eye value.
- Figure 2 temperatures: the marker within 1.5 s.

**Base tolerance.**

| Quantity | Base tolerance |
|---|---|
| Figure 2 loading | 0.003 kg/kg |
| Figure 2 temperatures | 0.7 °C |
| Figure 3 and 4 symbols | 0.005 |
| Model lines | 0.007 |

The tolerance is doubled when any reader flags the point as near-overlap, merged or inferred.

**Three readings.**

- If all three agree pairwise, the median is taken.
- If only one pair agrees, the two readings that agree prevail. The value is the mean of that pair, and the third reading is recorded as outvoted.
- If no pair agrees, the median is taken and the uncertainty is widened to max(2 × base, half the range).

**Two readings.**

- If they agree, the mean is taken.
- If they disagree and the third reader saw the marker directly (flag clear or near-overlap), the third reading is used, with uncertainty max(2 × base, half the difference).
- Otherwise the mean is taken, with uncertainty half the difference plus the base.

**One reading.** Where only version 1 has a value, it is kept as version 1 printed it, not anchor-corrected. The point is flagged `inferred` with uncertainty 3 × base. The note gives version 1 minus the model line at that time, because most such values follow the model line.

**Fixed points.** At t = 0 the normalized experimental loading is 1 by definition, with uncertainty 0. The model lines at t = 0 keep the version-1 value.

**Output.** Each point in `adjudication_per_point.csv` carries the three readings, the version-2 value, `uncertainty_abs`, `flag` (clear, near-overlap, merged or inferred), the method used and a note.

**Point counts by method.**

| Method | Points |
|---|---|
| all three agree | 391 |
| two agree, one outvoted | 34 (v1 outvoted 30 times, the auditor 0, the third reader 4) |
| three disagree | 1 (Fig. 2 soybean at 9.884 s) |
| two readings agree | 148 |
| two readings disagree | 14 |
| single reading (v1) | 49 |
| defined at t = 0 | 24 |

## 6. What changed against version 1

"Changed" means |v2 - v1| is larger than the base tolerance. Smaller shifts come from the anchor convention and from taking the median, and they are listed in `adjudication_per_point.csv`.

**Figure 2.**

| Series | Points (n) | Changed | Detail |
|---|---|---|---|
| X sunflower | 25 | 1 | 240 s: 0.007415 → 0.0105 (inferred, ±0.003). The 200-230 s points are now flagged merged, ±0.006. t = 0: 0.4993 → 0.5016 (merged, ±0.006). |
| X soybean | 25 | 2 | 9.884 s: 0.4264 → 0.4209, a marker time error in version 1. 40.154 s: 0.2212 → 0.2250. 240 s: 0.007415 → 0.0080 (merged, ±0.006). t = 0: 0.4906 → 0.4924. |
| Tp sunflower | 16 | 0 | 19-64 s plateau mean 68.62 → 68.57 °C (offset from 68.7145: -0.15 °C). |
| Tg sunflower | 16 | 0 | All points -0.25 to -0.48 °C (anchor). |
| Tp soybean | 17 | 1 | t = 0: 68.68 → 69.74 °C, on the axis, merged. 10-40 s plateau 68.85 → 68.91 °C (+0.19 °C). |
| Tg soybean | 17 | 1 | t = 0: 119.57 → 117.87 °C, on the axis, merged. All others -0.39 to -0.64 °C (anchor). |

**Figure 3 (sunflower).**

| Series | Points (n) | Changed | Detail |
|---|---|---|---|
| 93 °C | 25 | 1 | 10 s: +0.007. The 210-240 s points are inferred (unreadable) and kept at the version-1 value. |
| 101 °C | 25 | 2 | 10 s: +0.015. 160 s: -0.005 (merged). 180-240 s inferred. |
| 110 °C | 25 | 1 | 10 s: +0.018 (merged). 160-240 s inferred. Version 1 there lies +0.024 to +0.038 above the model line. |
| 121 °C | 25 | 5 | 20 s: 0.706 → 0.674; version 1 had no marker there, and the grey triangle has its apex at row 206 behind the 125 °C circle. 50 s: -0.007; version 1 duplicated the 125 °C circle. 80 s: -0.006. 160 s: +0.010. 170 s: +0.007. 190-240 s inferred. |
| 125 °C | 25 | 4 | 160-200 s: +0.005 to +0.012 (merged). |
| 136 °C | 25 | 2 | 160 s: -0.009. 170 s: -0.006. |
| Model lines, 6 series | 127 | 0 | All shifts ≤0.004. Samples at or below 0.008 are flagged merged with the axis. |

**Figure 4 (soybean).**

| Series | Points (n) | Changed | Detail |
|---|---|---|---|
| 92 °C | 25 | 11 | 70 s: -0.011; 90 s: -0.007; 140-210 s: -0.005 to -0.019, now on the open squares; 220 s: -0.010 (merged, ±0.015). |
| 101 °C | 25 | 2 | 10 s: -0.013. 100 s: +0.005. |
| 109 °C | 25 | 4 | 20 s: +0.014; 40 s: +0.011; 80 s: -0.018; 100 s: -0.008. 130-240 s inferred (unreadable). |
| 120 °C | 25 | 11 | 20 s: -0.032; 40 s: -0.015; 50 s: -0.029; 60 s: -0.017; 130-150 s: +0.011 to +0.019; 10 s and 160-230 s merged, +0.001 to +0.010. |
| 125 °C | 25 | 3 | 110 s: +0.015; 120 s: +0.012; 130 s: +0.006 (merged). 150, 170 and 180 s inferred. |
| 136 °C | 25 | 1 | 220 s: -0.012. |
| Model lines, 6 series | 118 | 0 | All shifts ≤0.005. |

## 7. Uncertainty per point

`uncertainty_abs` in `adjudication_per_point.csv` is the adjudicated half-width:

| Quantity | clear | near-overlap or merged | inferred |
|---|---|---|---|
| Figure 2 loading | ±0.003 kg/kg | ±0.006 | ±0.009 |
| Figure 2 temperatures | ±0.7 °C | ±1.4 °C | — |
| Figure 3 and 4 symbols | ±0.005 | ±0.010 | ±0.015 |
| Model lines | ±0.007 | ±0.014 | ±0.021 |

Where the readings disagree, the widened values described in §5 are used instead. The sunflower final loading is set to ±0.003, as §3.2 explains.

These bounds are larger than the spread of repeated reads, which is under 1 px for clear symbols. They include the calibration residual (≤0.26 px), the choice of anchor (±0.5 px) and the judgment of overlaps.

The 49 points flagged `inferred` because only version 1 has a value should not be scored as data. None of the three readers could read the marker there, and most of those values follow the model line.

## 8. Conditions (`conditions.csv`, each row with the article's exact wording in `notes`)

`conditions.csv` has the same columns as version 1. The `notes` column now carries the article's wording for each condition. The main ones:

- **Mean diameters.** "For sunflower, approximately 97.1% of the particles had a Deq value in the range 1 to 3 10-3 m, with an average equivalent diameter of 1.81 10-3 m. For soybean, the percentage of particles that had a Deq value between 1 and 3 10-3 m was 95.1%, with an average equivalent diameter of 1.77 10-3 m." (pp. 4-5). **1.81 mm is sunflower and 1.77 mm is soybean. They are not interchangeable.**
- **Sphericity.** "The average sphericity shape factor (ψ) was 0.74 for both sunflower and soybean meals." (p. 5). It was measured by image analysis on the 2019 study's own sunflower meal and on its own soybean meal (Sec. 2.1). It was not measured on soybean meal and then borrowed for sunflower.
- **Sample and holder.**
  - "…(~ 30 g of sample) but a thickness that assure a thin layer condition (~ 3–5Dp)." (p. 2)
  - "rectangular sample-holder (cross area: 0.011m2) made of stainless steel mesh" (p. 2)
  - Chamber "110 mm length, 100 mm inner diameter" (p. 2)
- **Initial ratio.** "embedded with hexane to a solvent-to-meal ratio of 0.45–0.5 (w : w), higher than … the industrial extractor outlet stream (0.3–0.4 w : w)" (pp. 2-3). The basis is not printed, so the unit is now `kg_hexane_per_kg_meal_basis_not_stated`. Version 1 had per kg dry meal.
- **Logging.**
  - "the weight loss was measured continuously at second intervals" (p. 2). The value is now `not_printed`. Version 1 carried 2 s, which came from a brief.
  - "recorded the temperature of the meal and the superheated vapor every 10 s" (p. 3).
  - A new row, `figure2_temperature_marker_spacing` = 15 s, has origin `digitized`.
- **Gas flow.** "The tests were performed at different flow rates and temperatures, ranging from 90 °C to 135 °C." (p. 3). No flow rate is printed. There is a new row `gas_flow_rate` = `not_printed`.
- **Not controlled across runs.** "gas flow rate, amount of meal, and initial solvent content were not exactly the same in all the experiments" (p. 5).
- **Bed porosity.** The article prints only "ranging from 0.638 to 0.709" (p. 4). The per-meal values are now labelled `derived` (1 - ρb/ρp), where version 1 had `tabulated`.
- **Rows kept as in version 1, with the wording added.** Densities, particle porosities, bulk densities (with the ×1000 unit error), oil contents, Xc = 0.20, Tb = 68.7 °C, preheating, the 0.01 g readability and the Figure 2 temperatures.
- **New rows, all tabulated.** Superheater range 80-140 °C, run range 90-135 °C, Deq range 0.5-4 mm, and the model-versus-data bands (±14/±11 % at 60 s; ±15/±10 % on t90).

## 9. Consumer inputs taken from these curves and conditions (for re-scoring; nothing was re-run)

These uses were located read-only on 2026-09-27. File and line numbers are as found then.

**item_05, `rerun_2026-09-22_rig/run_rerun_faner.py`.**

- **Final loading.** `FANER2019_FINAL_LOADING = 0.007415` (line 216) is applied to both Figure 2 traces (lines 994, 1023) as the measured plateau. RESULT_2026-09-22.md line 50 prints "measured final 0.007415".
  - Version 2: sunflower 0.0105 ± 0.003 (inferred) and soybean 0.0080 ± 0.006.
  - The two traces no longer share one value.
- **Initial loading.** X0 is taken from the first curve point (line 705).
  - Sunflower 0.4993 → 0.5016 and soybean 0.4906 → 0.4924, both merged and ±0.006.
- **Constant-rate slope.** The least-squares fit over the constant-rate window (lines 738-741) will move slightly.
  - Over X ≥ 0.20: sunflower -4.769e-3 → -4.808e-3 s-1 (+0.8 %) and soybean -6.845e-3 → -6.749e-3 s-1 (-1.4 %), mostly through the soybean point at 9.884 s.
  - The X = 0.20 crossing moves from 65.08 to 65.06 s (sunflower) and from 45.46 to 45.86 s (soybean).
- **Tp plateau and Tg mean** (lines 748-753).
  - Tp: sunflower 68.62 → 68.57 °C and soybean 68.85 → 68.91 °C.
  - Tg mean up to the crossing: sunflower 101.82 → 101.55 °C and soybean 120.26 → 119.44 °C. The soybean shift includes the t = 0 marker on the axis, 119.57 → 117.87 °C. Leaving that marker out, the shift is about -0.5 °C.
  - The sunflower and soybean diameters (1.81e-3 and 1.77e-3, lines 189 and 203) are already species-correct. They are unchanged.
- **Conditions.** The sample mass of 30 g, the 0.011 m2 area and the 3-5 Dp bed are unchanged.

**item_22 (`collect_f32.py`, `analyze_f32.py`).**

- It scores only the two Figure 2 traces, through item_05's conditions (collect_f32.py lines 59-60). It uses no Figure 3 or 4 series.
- The corrections to the Figure 4 series at 92, 109 and 120 °C therefore do not reach item_22 directly.
- item_22 inherits the Figure 2 changes: the final loadings, X0, and the soybean point at 9.884 s. It also inherits item_05's acceptance band, and the sweep range for X0 (0.40-0.55) still covers the new X0.
- The auditor names a validation script that computes t50 and t90 from Figures 3 and 4. That script is the consumer of the Figure 3 and 4 changes, and its t50 and t90 for Fig. 4 at 92, 109 and 120 °C, and Fig. 3 at 121 °C, need recomputing. Leave out the points flagged `inferred`.

**item_23, `run_item23_oil_arm.py`.**

- `MEASURED_FINAL_FANER2019 = 0.007415` (line 46) is applied to both conditions (lines 55, 60).
- Version 2 gives sunflower 0.0105 (range 0.0075-0.0125) and soybean 0.0080 ± 0.006.
- The ratios in RESULT.md line 116 are 2.79 against the floor without the oil arm (2.6546e-3) and 1.77 against the floor with it (4.2010e-3), with a log-gap closure of 44.7 %. At 0.0105 they become about 3.96 and 2.50, and the log-gap closure about 34 %; the arithmetic is only indicative and was not re-run.
- Over the plausible range 0.0075-0.0125 the ratios run from 2.83 to 4.71 and from 1.79 to 2.98.
- The oil fractions 0.018 and 0.0195 (conditions.csv) are unchanged.

**items 29-32.**

- These score the thesis traces (faner2008), not these curves.
- They take `D_EQ_MM = 1.77` (item_30 run_item30.py line 89; item_31 line 91; item_32 line 80) and `PSI = 0.74` (item_30 line 88; item_32 line 79). They use them for both species to build the 2019-diameter spheres.
- Per the article, the 2019 sunflower diameter is 1.81 mm. The 2019-diameter volume-to-surface sphere for sunflower is therefore 0.74 × 1.81/2 = 0.670 mm, not 0.655 mm, and that conventions row needs re-labelling or re-running for the sunflower trace.

**item_34.** It uses no value from this dataset beyond item_23's soybean oil fraction.

**Paper (read-only locations, for the next agent).**

- **Sphericity described as measured on soybean meal:** `sec_validation.tex` line 258 ("measured on industrial soybean meal") and line 415, and supplement lines 6123 and 6401.
  - Per the article it was measured on both 2019 meals.
  - It is still not measured on the thesis samples, so the non-claim for the thesis samples stands.
- **"1.77 mm for both traces … borrows both factors from the soybean meal":** `sec_validation.tex` line 267, and supplement lines 5915-5922 and 6364.
- **Final loading, "1.77 times above it instead of … 2.79":** supplement lines 2358-2359, the item_23 material.
- **"Sample mass is logged every 2 s":** supplement line 1794. The article says "at second intervals".

## 10. Files

- `curves.csv` has the same columns, row order and point_ids as version 1. Only the `value` column changes.
- `conditions.csv` has the same columns as version 1. Rows were added and notes were rewritten (§8).
- `adjudication_per_point.csv` has, per point: the version-1, anchor-corrected version-1, auditor, third-reader and version-2 values, the uncertainty, the flag, the method and a note.
- `qa_overlay_figure2.png`, `qa_overlay_figure3.png` and `qa_overlay_figure4.png` show every version-2 point on the PDF image stream at 2x.
  - Marker points are shown as crosses: green for clear, orange for near-overlap, magenta for merged, red for inferred.
  - Model-line samples are shown as blue circles.
- `build_v2.py` rebuilds the CSVs, the overlays and the manifest from `third_reading/` and the read-only inputs.
- `third_reading/` holds:
  - the calibration (`calib.json`, `calib_raw.json`);
  - the seeds and readings (`seeds_f3.txt`, `seeds_f4.txt`, `mine_f3.json`, `mine_f4.json`, `mine_f2_loading.json`, `f2load_mid.json`, `f2temp.json`, `f2temp_full.json`, `lines_f3.json`, `lines_f4.json`);
  - the scripts as run (`scripts/`);
  - the byte-exact PDF image streams (`source_images/`).
- The other two readings are read in place and are not copied:
  - version 1, `../faner2019/curves.csv`, sha256 `529e32f3dae5adad4760762a02f47439dd0dd09a91242a133a76f6c7b2aab783`;
  - the auditor, `../REDIGITIZATION_AUDIT_ARTICLES_2026-09-27_faner2019_{fig2_loading,fig2_temperature,fig3_fig4_experimental,fig3_fig4_model_lines}.csv`.
- `SHA256SUMS` lists the pinned source PDF and every file in this folder except itself.

## 11. Known source defects and cautions (carried over from version 1, amended)

1. The printed bulk densities carry a ×1000 unit error (0.373 and 0.338 kg/m3).
2. Equation (9)'s third branch is printed with the inequality reversed. Equation (8) has an anomalous bed-porosity factor. Both are recorded here, not repaired.
3. The series of Figures 3 and 4 are not controlled single-factor sweeps. In Figure 4 the 109 °C and 120 °C series cross between 70 and 80 s.
4. The samples were preheated, so the curves show no warming-up stage.
5. The mass-log interval is not printed (§8). The temperature markers in Figure 2 are plotted every 15 s, although the text says they were logged every 10 s.
6. The model-versus-data bands (±14/±11 %, ±15/±10 %) are comparison context, not digitization tolerances.
7. From about 195-200 s the Figure 2 sunflower loading is inferred from the tips of the crosses, and at 240 s it is not readable. Do not quote the final loading to more than two significant figures.

## Integrity

- `curves.csv`, `conditions.csv`, `adjudication_per_point.csv` and all other files: see `SHA256SUMS`.
- **Date and extractor:** 2026-09-27, third (adjudicating) reader, an Opus agent. Built with the pinned environment (PYTHONHASHSEED=1, single-thread BLAS, `.venv` interpreter).
- Any edit to the cleaned data requires re-running `build_v2.py` and adding a dated note here.
