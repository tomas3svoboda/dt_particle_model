# Faner multiplier: full-source SSE and calibration-window analysis

Loading and temperature SSE have different units and are reported separately.
All 25 loading and 17 particle-temperature points remain in these scores.

| Setup | Multiplier | Loading SSE | Temperature SSE, K2 | Core end, s |
|---|---:|---:|---:|---:|
| ref_P_F16 | 0.145785 | 0.01624729 | 1506.863 | 85.514 |
| f0p10 | 0.1 | 0.05283304 | 301.8316 | 122.29 |
| fstar4_0p1795 | 0.1795 | 0.03984012 | 3638.954 | 70.344 |
| f0p20 | 0.2 | 0.05998769 | 4979.07 | 63.578 |
| f0p30 | 0.3 | 0.1598246 | 10267.54 | 43.791 |
| R_F20_fg | 0.3843 | 0.03788581 | 3507.793 | 76.482 |
| R_F40 | 0.3843 | 0.03789672 | 3508.061 | 76.193 |
| phi_0p20000 | 0.2 | 0.07102745 | 504.196 | 137.59 |
| phi_0p25000 | 0.25 | 0.02176075 | 257.5935 | 112.24 |
| phi_0p28000 | 0.28 | 0.01398655 | 700.8702 | 101.28 |
| phi_0p28500 | 0.285 | 0.01363472 | 801.8699 | 99.702 |
| phi_0p29000 | 0.29 | 0.0135009 | 909.1337 | 98.146 |
| phi_0p30000 | 0.3 | 0.01392434 | 1139.503 | 95.191 |
| phi_0p34000 | 0.34 | 0.02188223 | 2203.617 | 85.145 |
| phi_0p36500 | 0.365 | 0.03028632 | 2935.72 | 80.009 |

## Early-time fitting choices at R = 1.5 mm

Bounds below propagate reading half-widths; they are not confidence intervals.

| Window | Slope-matched phi | Reading-box slope range | Fixed-initial-state SSE phi | Quality-scaled SSE phi | All early boxes feasible? |
|---|---:|---|---:|---:|---|
| original_four_points | 0.384285 | 0.360330 to 0.408240 | 0.349062 | 0.356426 | False |
| omit_t0 | 0.374581 | 0.342532 to 0.406631 | 0.349062 | 0.356426 | False |
| omit_29p961_s | 0.410797 | 0.386823 to 0.434771 | 0.350984 | 0.360911 | False |
| omit_both_endpoints | 0.432502 | 0.384858 to 0.480146 | 0.350984 | 0.360911 | False |

The first model loading is fixed at 0.475; the source reads 0.492364 +/- 0.006.
Dropping t=0 from prediction SSE removes a constant residual and cannot move its optimum.
The early reading boxes are also mutually incompatible with the fixed-state feed response after dropping t=0.
A lower global-SSE multiplier is a boundary calibration with a visible early/late tradeoff; reading precision alone does not identify it.

## Ranking sensitivity to the reading bounds

The listed candidates can minimize SSE among the actual simulated sample for some source readings within their existing bounds.
This discrete ranking sensitivity is not a confidence interval or a claim that every residual lies within the bounds.

- solvent_loading_unweighted: best sampled phi 0.29; reading-ambiguous sampled values [0.285, 0.29, 0.3].
- solvent_loading_reading_scaled: best sampled phi 0.28; reading-ambiguous sampled values [0.28, 0.285, 0.29].
- particle_temperature_unweighted: best sampled phi 0.25; reading-ambiguous sampled values [0.25].
- particle_temperature_reading_scaled: best sampled phi 0.25; reading-ambiguous sampled values [0.25].
- joint_equal_quantity_reading_scaled: best sampled phi 0.25; reading-ambiguous sampled values [0.25].