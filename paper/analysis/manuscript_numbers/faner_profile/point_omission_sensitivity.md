# Effect of point omission on full-loading-curve calibration

Every final reported score retains all 25 points. Omissions below change only the parameter-estimation objective.

| Fit points | Best sampled multiplier | All-25 SSE | Core end, s |
|---|---:|---:|---:|
| all_25 | 0.29 | 0.01350090 | 98.146 |
| omit_t0 | 0.29 | 0.01350090 | 98.146 |
| omit_29p961_s | 0.285 | 0.01363472 | 99.702 |
| omit_t0_and_29p961_s | 0.285 | 0.01363472 | 99.702 |
| omit_40p154_s | 0.285 | 0.01363472 | 99.702 |
| omit_29p961_and_40p154_s | 0.28 | 0.01398655 | 101.277 |

The 29.961 s and 40.154 s markers were already flagged near-overlap with half-width 0.006 kg/kg.
This supports treating their precision differently; it does not prove either observation invalid.
This whole-curve SSE analysis is distinct from omitting 29.961 s in a four-point slope fit.
