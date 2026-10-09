# Item 5, detail of 2026-09-24 -- the loading at which the law is evaluated, step by step

`physically_qualifying: false. plant_predictive: false.`
Read `RESULT.md` beside this file first. This document carries the method, the
verification, and all twenty-four marches in full. Every number below is
written by `make_detail.py` from `surface_loading.json` and
`surface_loading_window.csv`; none is transcribed by hand.

## 1. What was asked

The limits audit of 2026-09-22 carries limit 2 at weight **fatal**: the
measured diffusivity law is evaluated 6.4 to 12.8 times above the largest
loading at which it was measured, and the measured range and the falling-rate
windows do not overlap at all. The audit computed that factor on the window's
loading coordinate, which is the particle's volume-mean loading. The
manuscript's falling-rate line is the diffusivity case
`law_log_loading_surface`, which evaluates the law at the **surface** loading.
The two are not the same number, and nobody had measured the second.

## 2. The instrument

`run_surface_loading.py` imports `rerun_2026-09-22_rig/run_rerun_faner.py` as a
module; its `main()` is guarded and is never called, so the constants, the
frozen anchors, the property calls, the diffusivity-case construction and the
law record are all the committed objects. The committed `march_tail` is taken
by `inspect.getsource`, patched textually with recorder lines, and re-bound by
`exec` in a copy of the rig module's own globals, so the marched function sees
exactly the namespace the committed one saw.

Fourteen lines are inserted, each ending in the marker `# +REC`. The script
asserts, before it marches anything, that removing exactly the lines carrying
that marker restores the committed source string byte for byte. Two anchors are
used and each is checked to be unique in the source: the `try:` that guards the
certified step, and the `ppm = K.inventory_ppm(...)` that closes it. The
recorder captures, for every **accepted** step -- a step refused and retried at
a larger stride is written down only once it is accepted:

| recorded | meaning |
|---|---|
| `mean_pre_kg_kg` | the volume-mean loading entering the stride, `prev_ppm * 1e-6`; the audit's basis |
| `surface_pre_kg_kg` | the surface-cell loading entering the stride, `storage(c[-1], T, p) / rho_dm_p` |
| `eval_pre_kg_kg` | the loading the law was actually evaluated at for this stride |
| `d_eff_m2_s` | the coefficient the stride was marched at |
| `t_model_k`, `dt_s`, `t_before_s`, `t_after_s` | the model temperature, the stride and its endpoints |
| `mean_post_kg_kg`, `surface_post_kg_kg` | the same two loadings after the accepted step, before any temperature remap |

For `law_log_loading_surface` the evaluation loading is the surface loading;
the columns are kept separate anyway so the file reads the same for any case.

## 3. The verification, run before anything was read

Twenty-four marches: the two thesis traces at 136 C, the `sensible` temperature
case, the `law_log_loading_surface` diffusivity case, and all twelve stage-1
film closures of the rig run. Each is compared against the committed
`rerun_residuals.csv` and `rerun_tail_residuals.csv` by exact float identity
through the CSV round trip:

* the march: step count, remap count, `d_start`, `d_end`, end loading, end
  particle temperature, first stride, retarded time constant, start activity,
  status -- **10 checks**;
* the re-scored window: point count, window loading span, RMS residual, RMS
  normalized by the span, mean signed residual, maximum absolute residual,
  points inside the reading band -- **7 checks**.

**408 of 408 exact**, twenty-four marches times seventeen checks. The script
raises and writes nothing on a single mismatch. The window residual is rebuilt
the way `run_rerun_tail_residuals.py` builds it: stage 1 as the two points
`(0, X_0)` and `(t_c, X_c)`, then the tail trace sampled every five accepted
steps plus its endpoint, linearly interpolated at the measured sample times at
or after the measured crossing of the source's declared `X_c = 0.20`.

## 4. The law and its fitted region

| quantity | value |
|---|---|
| expression | `ln D* = ln A - (Ea/Rg)(1/T) + beta*ln(X/0.01)` |
| meal set the march uses | soybean, n = 17 |
| loading reference `x_ref` | 0.01 kg/kg |
| loading exponent | 0.9592058934994389 |
| activation energy | 0.3058629946934234 kJ/mol |
| residual sigma in ln D* | 0.410374372529405 |
| measured loading range | 0.0007676965718689591 to 0.015616094592894784 kg/kg |
| fitted temperature range | 50 to 95 C |

The temperature extrapolation, stated separately as the brief asks, in `1/T`:

| quantity | value |
|---|---|
| fitted range in 1/T | 2.716284e-03 to 3.094538e-03 1/K |
| fitted span in 1/T | 3.782540e-04 1/K |
| marched at 136 C | 2.444091e-03 1/K |
| beyond the hot end | 2.721927e-04 1/K |
| as a fraction of the fitted span | 0.7196 |
| kelvin above the highest fitted temperature | 41.0 |

The design region, measured rather than described. The law's two regressors are
correlated (Pearson 0.5904, VIF 1.5352 in the law record), so "inside in loading,
outside in temperature" is not the same as "inside the fitted design". The
leverage `h = x'(X'X)^-1 x` in the law's own parameterization
`ln D = b0 + b1 * (-1/(Rg T)) + b2 * ln(X/x_ref)` measures the joint position.
The recorded `cov_unscaled` was confirmed to be `(X'X)^-1` on the seventeen
soybean points to `|I - (X'X)C| = 5.966e-10` before any leverage was computed.

| quantity | value |
|---|---|
| fitted points | 17, 3 parameters |
| mean fitted leverage | 0.176471 (= 3/17) |
| smallest fitted leverage | 0.068496 |
| largest fitted leverage | **0.304594** |
| two-sided 95 per cent prediction factor on D at h = 0 | 2.4115 |
| the same at the largest fitted leverage | 2.7330 |

## 5. Table 1 -- the loading at which the law was evaluated, over the scored window

The window is the committed one: every measured sample at or after the measured
crossing of `X_c = 0.20`. "Steps" counts accepted march steps whose post-step
absolute time lies between the first and the last measured sample of that
window. "Duration inside" weights those steps by their own stride, because the
stride grows by 1.12 and they are not of equal length.

| trace | stage-1 closure | steps | evaluation loading, kg/kg | / measured max | / measured min | inside, steps | inside, duration | window RMS |
|---|---|---:|---|---|---|---:|---:|---:|
| sunflower 136 C | `frozen_1.0` | 48 | 1.062e-03 to 6.370e-03 | 0.068 to 0.408 | 1.38 to 8.30 | 48 of 48 | 1.0000 | 0.1102 |
| sunflower 136 C | `lit_whitaker` | 48 | 1.097e-03 to 2.193e-02 | 0.070 to 1.404 | 1.43 to 28.56 | 47 of 48 | 0.9831 | 0.0994 |
| sunflower 136 C | `source_central` | 48 | 1.097e-03 to 2.189e-02 | 0.070 to 1.402 | 1.43 to 28.52 | 47 of 48 | 0.9831 | 0.0992 |
| sunflower 136 C | `frozen_rig_state` | 48 | 1.085e-03 to 1.408e-02 | 0.069 to 0.902 | 1.41 to 18.34 | 48 of 48 | 1.0000 | 0.0968 |
| sunflower 136 C | `source_eq427` | 48 | 1.085e-03 to 1.409e-02 | 0.069 to 0.902 | 1.41 to 18.36 | 48 of 48 | 1.0000 | 0.0968 |
| sunflower 136 C | `source_eq426` | 49 | 1.104e-03 to 3.109e-02 | 0.071 to 1.991 | 1.44 to 40.50 | 47 of 49 | 0.9685 | 0.1052 |
| sunflower 136 C | `source_low` | 48 | 1.112e-03 to 3.109e-02 | 0.071 to 1.991 | 1.45 to 40.50 | 46 of 48 | 0.9678 | 0.1070 |
| sunflower 136 C | `source_high` | 48 | 1.083e-03 to 1.388e-02 | 0.069 to 0.889 | 1.41 to 18.08 | 48 of 48 | 1.0000 | 0.0967 |
| sunflower 136 C | `lit_bird` | 48 | 1.045e-03 to 3.012e-03 | 0.067 to 0.193 | 1.36 to 3.92 | 48 of 48 | 1.0000 | 0.1287 |
| sunflower 136 C | `lit_bradshaw_myers` | 48 | 1.037e-03 to 2.091e-03 | 0.066 to 0.134 | 1.35 to 2.72 | 48 of 48 | 1.0000 | 0.1476 |
| sunflower 136 C | `lit_whitaker_dry_area` | 43 | 1.345e-03 to 1.994e-01 | 0.086 to 12.766 | 1.75 to 259.68 | 39 of 43 | 0.9270 | 0.7789 |
| sunflower 136 C | `source_central_dry_area` | 43 | 1.343e-03 to 1.994e-01 | 0.086 to 12.766 | 1.75 to 259.68 | 39 of 43 | 0.9270 | 0.7708 |
| soybean 136 C | `frozen_1.0` | 48 | 1.053e-03 to 4.238e-03 | 0.067 to 0.271 | 1.37 to 5.52 | 48 of 48 | 1.0000 | 0.1413 |
| soybean 136 C | `lit_whitaker` | 48 | 1.113e-03 to 3.109e-02 | 0.071 to 1.991 | 1.45 to 40.50 | 46 of 48 | 0.9678 | 0.2059 |
| soybean 136 C | `source_central` | 48 | 1.074e-03 to 9.585e-03 | 0.069 to 0.614 | 1.40 to 12.49 | 48 of 48 | 1.0000 | 0.1629 |
| soybean 136 C | `frozen_rig_state` | 48 | 1.097e-03 to 2.191e-02 | 0.070 to 1.403 | 1.43 to 28.54 | 47 of 48 | 0.9831 | 0.1874 |
| soybean 136 C | `source_eq427` | 48 | 1.097e-03 to 2.192e-02 | 0.070 to 1.404 | 1.43 to 28.55 | 47 of 48 | 0.9831 | 0.1877 |
| soybean 136 C | `source_eq426` | 48 | 1.075e-03 to 9.613e-03 | 0.069 to 0.616 | 1.40 to 12.52 | 48 of 48 | 1.0000 | 0.1633 |
| soybean 136 C | `source_low` | 49 | 1.088e-03 to 2.154e-02 | 0.070 to 1.380 | 1.42 to 28.06 | 48 of 49 | 0.9834 | 0.1790 |
| soybean 136 C | `source_high` | 48 | 1.069e-03 to 8.918e-03 | 0.068 to 0.571 | 1.39 to 11.62 | 48 of 48 | 1.0000 | 0.1552 |
| soybean 136 C | `lit_bird` | 48 | 1.039e-03 to 2.223e-03 | 0.067 to 0.142 | 1.35 to 2.90 | 48 of 48 | 1.0000 | 0.1378 |
| soybean 136 C | `lit_bradshaw_myers` | 48 | 1.039e-03 to 2.208e-03 | 0.067 to 0.141 | 1.35 to 2.88 | 48 of 48 | 1.0000 | 0.1379 |
| soybean 136 C | `lit_whitaker_dry_area` | 41 | 1.437e-03 to 1.994e-01 | 0.092 to 12.766 | 1.87 to 259.68 | 36 of 41 | 0.8984 | 1.1061 |
| soybean 136 C | `source_central_dry_area` | 45 | 1.237e-03 to 1.994e-01 | 0.079 to 12.766 | 1.61 to 259.68 | 41 of 45 | 0.9303 | 0.6058 |

## 6. Table 2 -- the same on the audit's basis, the volume-mean loading

| trace | stage-1 closure | volume-mean loading, kg/kg | / measured max | inside, steps |
|---|---|---|---|---:|
| sunflower 136 C | `frozen_1.0` | 0.1118 to 0.1657 | 7.16 to 10.61 | 0 of 48 |
| sunflower 136 C | `lit_whitaker` | 0.1147 to 0.1769 | 7.35 to 11.33 | 0 of 48 |
| sunflower 136 C | `source_central` | 0.1147 to 0.1769 | 7.35 to 11.33 | 0 of 48 |
| sunflower 136 C | `frozen_rig_state` | 0.1138 to 0.1731 | 7.29 to 11.09 | 0 of 48 |
| sunflower 136 C | `source_eq427` | 0.1138 to 0.1732 | 7.29 to 11.09 | 0 of 48 |
| sunflower 136 C | `source_eq426` | 0.1149 to 0.1801 | 7.36 to 11.53 | 0 of 49 |
| sunflower 136 C | `source_low` | 0.1156 to 0.1801 | 7.40 to 11.53 | 0 of 48 |
| sunflower 136 C | `source_high` | 0.1137 to 0.1731 | 7.28 to 11.08 | 0 of 48 |
| sunflower 136 C | `lit_bird` | 0.1097 to 0.1587 | 7.02 to 10.16 | 0 of 48 |
| sunflower 136 C | `lit_bradshaw_myers` | 0.1085 to 0.1552 | 6.95 to 9.94 | 0 of 48 |
| sunflower 136 C | `lit_whitaker_dry_area` | 0.1220 to 0.1994 | 7.81 to 12.77 | 0 of 43 |
| sunflower 136 C | `source_central_dry_area` | 0.1219 to 0.1994 | 7.81 to 12.77 | 0 of 43 |
| soybean 136 C | `frozen_1.0` | 0.1107 to 0.1620 | 7.09 to 10.37 | 0 of 48 |
| soybean 136 C | `lit_whitaker` | 0.1156 to 0.1801 | 7.41 to 11.53 | 0 of 48 |
| soybean 136 C | `source_central` | 0.1129 to 0.1695 | 7.23 to 10.86 | 0 of 48 |
| soybean 136 C | `frozen_rig_state` | 0.1147 to 0.1769 | 7.35 to 11.33 | 0 of 48 |
| soybean 136 C | `source_eq427` | 0.1147 to 0.1769 | 7.35 to 11.33 | 0 of 48 |
| soybean 136 C | `source_eq426` | 0.1129 to 0.1695 | 7.23 to 10.86 | 0 of 48 |
| soybean 136 C | `source_low` | 0.1139 to 0.1768 | 7.29 to 11.32 | 0 of 49 |
| soybean 136 C | `source_high` | 0.1126 to 0.1691 | 7.21 to 10.83 | 0 of 48 |
| soybean 136 C | `lit_bird` | 0.1087 to 0.1556 | 6.96 to 9.96 | 0 of 48 |
| soybean 136 C | `lit_bradshaw_myers` | 0.1086 to 0.1555 | 6.96 to 9.96 | 0 of 48 |
| soybean 136 C | `lit_whitaker_dry_area` | 0.1238 to 0.1994 | 7.93 to 12.77 | 0 of 41 |
| soybean 136 C | `source_central_dry_area` | 0.1198 to 0.1994 | 7.67 to 12.77 | 0 of 45 |

## 7. Table 3 -- the joint position, and the temperature axis

`h` is the leverage of the evaluation point in the law's own two regressors,
against a largest fitted value of 0.3046. The last two columns count window steps
whose model temperature lies inside the fitted 50 to 95 C, and steps inside on
both axes at once.

| trace | stage-1 closure | h at the surface, min .. med .. max | 95 pc factor | h at the volume mean, med .. max | 95 pc factor | model T, C | T inside | both inside |
|---|---|---|---:|---|---:|---|---:|---:|
| sunflower 136 C | `frozen_1.0` | 0.361 .. 0.704 .. 0.727 | 3.180 | 4.37 .. 4.43 | 7.780 | 93.18 to 135.33 | 1 of 48 | 1 of 48 |
| sunflower 136 C | `lit_whitaker` | 0.330 .. 0.625 .. 0.713 | 3.165 | 4.24 .. 4.34 | 7.645 | 74.30 to 134.41 | 5 of 48 | 4 of 48 |
| sunflower 136 C | `source_central` | 0.331 .. 0.626 .. 0.714 | 3.165 | 4.24 .. 4.34 | 7.647 | 74.33 to 134.43 | 5 of 48 | 4 of 48 |
| sunflower 136 C | `frozen_rig_state` | 0.336 .. 0.651 .. 0.719 | 3.171 | 4.30 .. 4.36 | 7.678 | 79.22 to 134.76 | 4 of 48 | 4 of 48 |
| sunflower 136 C | `source_eq427` | 0.336 .. 0.651 .. 0.719 | 3.171 | 4.30 .. 4.36 | 7.678 | 79.21 to 134.75 | 4 of 48 | 4 of 48 |
| sunflower 136 C | `source_eq426` | 0.325 .. 0.612 .. 0.719 | 3.171 | 4.17 .. 4.32 | 7.621 | 70.19 to 134.18 | 7 of 49 | 5 of 49 |
| sunflower 136 C | `source_low` | 0.324 .. 0.608 .. 0.718 | 3.170 | 4.16 .. 4.32 | 7.615 | 70.15 to 134.01 | 7 of 48 | 5 of 48 |
| sunflower 136 C | `source_high` | 0.338 .. 0.656 .. 0.720 | 3.172 | 4.32 .. 4.37 | 7.689 | 79.50 to 134.83 | 4 of 48 | 4 of 48 |
| sunflower 136 C | `lit_bird` | 0.455 .. 0.732 .. 0.735 | 3.188 | 4.43 .. 4.52 | 7.910 | 110.99 to 135.62 | 0 of 48 | 0 of 48 |
| sunflower 136 C | `lit_bradshaw_myers` | 0.582 .. 0.736 .. 0.742 | 3.196 | 4.46 .. 4.59 | 8.008 | 122.07 to 135.74 | 0 of 48 | 0 of 48 |
| sunflower 136 C | `lit_whitaker_dry_area` | 0.290 .. 0.443 .. 2.291 | 4.938 | 3.50 .. 4.05 | 7.226 | 68.71 to 126.83 | 12 of 43 | 8 of 43 |
| sunflower 136 C | `source_central_dry_area` | 0.290 .. 0.444 .. 2.291 | 4.938 | 3.51 .. 4.05 | 7.230 | 68.71 to 126.90 | 12 of 43 | 8 of 43 |
| soybean 136 C | `frozen_1.0` | 0.391 .. 0.724 .. 0.731 | 3.184 | 4.40 .. 4.48 | 7.851 | 102.36 to 135.52 | 0 of 48 | 0 of 48 |
| soybean 136 C | `lit_whitaker` | 0.324 .. 0.606 .. 0.718 | 3.170 | 4.16 .. 4.32 | 7.612 | 70.12 to 133.97 | 7 of 48 | 5 of 48 |
| soybean 136 C | `source_central` | 0.343 .. 0.674 .. 0.723 | 3.175 | 4.33 .. 4.39 | 7.714 | 85.07 to 135.03 | 3 of 48 | 3 of 48 |
| soybean 136 C | `frozen_rig_state` | 0.330 .. 0.625 .. 0.713 | 3.165 | 4.24 .. 4.34 | 7.646 | 74.31 to 134.42 | 5 of 48 | 4 of 48 |
| soybean 136 C | `source_eq427` | 0.330 .. 0.625 .. 0.713 | 3.165 | 4.24 .. 4.34 | 7.645 | 74.30 to 134.42 | 5 of 48 | 4 of 48 |
| soybean 136 C | `source_eq426` | 0.343 .. 0.673 .. 0.723 | 3.175 | 4.33 .. 4.39 | 7.712 | 85.01 to 135.02 | 3 of 48 | 3 of 48 |
| soybean 136 C | `source_low` | 0.334 .. 0.636 .. 0.717 | 3.169 | 4.27 .. 4.35 | 7.665 | 74.62 to 134.66 | 5 of 49 | 4 of 49 |
| soybean 136 C | `source_high` | 0.354 .. 0.691 .. 0.726 | 3.179 | 4.35 .. 4.42 | 7.757 | 86.73 to 135.22 | 2 of 48 | 2 of 48 |
| soybean 136 C | `lit_bird` | 0.552 .. 0.735 .. 0.739 | 3.193 | 4.45 .. 4.57 | 7.979 | 119.87 to 135.71 | 0 of 48 | 0 of 48 |
| soybean 136 C | `lit_bradshaw_myers` | 0.556 .. 0.735 .. 0.740 | 3.193 | 4.45 .. 4.57 | 7.982 | 120.12 to 135.71 | 0 of 48 | 0 of 48 |
| soybean 136 C | `lit_whitaker_dry_area` | 0.287 .. 0.421 .. 2.291 | 4.938 | 3.38 .. 3.96 | 7.106 | 68.56 to 124.54 | 13 of 41 | 8 of 41 |
| soybean 136 C | `source_central_dry_area` | 0.299 .. 0.500 .. 2.291 | 4.938 | 3.72 .. 4.17 | 7.400 | 68.71 to 130.00 | 11 of 45 | 7 of 45 |

## 8. Table 4 -- the conservative window, and where the envelope is entered

The window of Tables 1 to 3 opens at the first measured sample at or after the
measured crossing. This table opens it at the **crossing itself**, which adds
the opening strides of the tail, at which the surface cell has only just left
the uniform critical loading `X_c = 0.19935` and the evaluation is therefore at
its largest. It is the harder reading of the same question and is reported
beside the other. The last column gives the tail time at which the surface
evaluation first falls at or below the measured maximum; it is the same in all
twenty-four marches, because the stride schedule and the surface boundary
condition are.

| trace | stage-1 closure | steps from the crossing | evaluation max, kg/kg | inside, steps | envelope entered at, s of tail time |
|---|---|---:|---|---:|---:|
| sunflower 136 C | `frozen_1.0` | 49 | 8.5462e-03 | 49 of 49 | 10.846 |
| sunflower 136 C | `lit_whitaker` | 49 | 3.1092e-02 | 47 of 49 | 10.846 |
| sunflower 136 C | `source_central` | 49 | 3.1092e-02 | 47 of 49 | 10.846 |
| sunflower 136 C | `frozen_rig_state` | 49 | 2.1284e-02 | 48 of 49 | 10.846 |
| sunflower 136 C | `source_eq427` | 49 | 2.1294e-02 | 48 of 49 | 10.846 |
| sunflower 136 C | `source_eq426` | 49 | 3.1092e-02 | 47 of 49 | 10.846 |
| sunflower 136 C | `source_low` | 49 | 1.9935e-01 | 46 of 49 | 10.846 |
| sunflower 136 C | `source_high` | 49 | 2.1070e-02 | 48 of 49 | 10.846 |
| sunflower 136 C | `lit_bird` | 48 | 3.0121e-03 | 48 of 48 | 10.846 |
| sunflower 136 C | `lit_bradshaw_myers` | 48 | 2.0913e-03 | 48 of 48 | 6.814 |
| sunflower 136 C | `lit_whitaker_dry_area` | 43 | 1.9935e-01 | 39 of 43 | 15.362 |
| sunflower 136 C | `source_central_dry_area` | 43 | 1.9935e-01 | 39 of 43 | 15.362 |
| soybean 136 C | `frozen_1.0` | 49 | 5.4885e-03 | 49 of 49 | 10.846 |
| soybean 136 C | `lit_whitaker` | 49 | 1.9935e-01 | 46 of 49 | 10.846 |
| soybean 136 C | `source_central` | 50 | 2.0560e-02 | 49 of 50 | 10.846 |
| soybean 136 C | `frozen_rig_state` | 49 | 3.1092e-02 | 47 of 49 | 10.846 |
| soybean 136 C | `source_eq427` | 49 | 3.1092e-02 | 47 of 49 | 10.846 |
| soybean 136 C | `source_eq426` | 50 | 2.0597e-02 | 49 of 50 | 10.846 |
| soybean 136 C | `source_low` | 50 | 3.1092e-02 | 48 of 50 | 10.846 |
| soybean 136 C | `source_high` | 49 | 1.2596e-02 | 49 of 49 | 10.846 |
| soybean 136 C | `lit_bird` | 49 | 2.5676e-03 | 49 of 49 | 6.814 |
| soybean 136 C | `lit_bradshaw_myers` | 49 | 2.5476e-03 | 49 of 49 | 6.814 |
| soybean 136 C | `lit_whitaker_dry_area` | 41 | 1.9935e-01 | 36 of 41 | 20.362 |
| soybean 136 C | `source_central_dry_area` | 45 | 1.9935e-01 | 41 of 45 | 15.362 |

## 9. What the measurement does and does not establish

**Establishes.** Over the scored falling-rate window the manuscript's reported
falling-rate line evaluates the measured diffusivity law at loadings inside the
range over which that law was measured, at 46 to 48 of the 48 or 49 accepted
steps and 96.8 to 100 per cent of the window's duration, at every stage-1
closure on the particle-mass area convention. The 6.4 to 12.8 factor of the
limits audit is a property of the volume-mean loading, which is 6.95 to 11.53
times the measured maximum over the same window and inside the measured range
at zero steps of every one of the twenty-four marches. The temperature
extrapolation is unchanged: 41 K above the hot end of the fitted range, 0.7196
of that range's own span in 1/T. Jointly, the surface evaluation sits at a
design leverage of 0.29 to 0.74 against a largest fitted value of 0.3046, and
the volume mean at 3.38 to 4.59.

**Does not establish.** Nothing here validates the march, the law, the
diffusivity value or any drying time. The surface loading is a model quantity
produced by the Dirichlet surface condition and the GAB isotherm, not a
measurement; this record measures where the law is read, not whether reading it
there is right. The surface state sits at hexane activity 1.0 while
Cardarelli's coefficients were measured at activities 0.10 to 0.81, so the
loading coordinate overlaps while the `(T, a)` state does not -- the law's
regressors are `1/T` and `ln X` only, which is why the leverage is the
operative measure. No space-resolved solve was re-run and no radial loading
profile was recorded. No number in any other folder, in the manuscript, or in
`src/`, `tests/`, `scripts/` or `docs/` was changed, and nothing was fitted.

## 10. Files

| file | what it holds |
|---|---|
| `run_surface_loading.py` | the driver: imports the committed rig module, re-binds `march_tail` with fourteen audited recorder lines, verifies, measures |
| `make_detail.py` | writes this document from the two files below |
| `surface_loading.json` | the whole report: the law, the design region, the temperature axis, the verification rows, the twenty-four window rows |
| `surface_loading_window.csv` | the twenty-four window rows, every column |
| `surface_loading_steps.csv` | 11 544 accepted steps, one row each: surface, volume-mean and evaluation loading, the coefficient, the stride, the model temperature, and whether the step is in the window |
