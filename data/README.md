# Numerical records

The lean Git tree contains the case catalog, source variants and derived figure data. Full histories are the eight checksum-pinned GitHub release downloads below. Each archive carries an inner `EVIDENCE_MANIFEST.sha256`; producer-relative records are kept byte for byte. No pickled checkpoints are archived.

- `final_reference.tar.xz`: 8 members; the complete `Faner_R15_N20_L1` reference.
- `final_ladder_records.tar.xz`: 196 members; JSON records and audits for all catalogued cases, plus the generated ladder report.
- `final_adopted_dt.tar.xz`: 99 members; cases `DT_N10_L1`, `DT_N20_L1`, `DT_N40_L1`, `DT_N80_L1`, `DT_N20_L1_WSAT025`, `DT_N20_L1_WSAT032`, `DT_N20_L1_WSAT075`, `DT_N20_L1_WSAT097`, `DT_N20_L1_DB190`, `DT_N20_L1_DB200`, `DT_N20_L1_DB600`.
- `final_adopted_dt_time_pair.tar.xz`: 18 members; cases `DT_N20_L2`, `DT_N40_L2`.
- `final_adopted_faner.tar.xz`: 40 members; cases `Faner_R15_N10_L1`, `Faner_R15_N20_L1`, `Faner_R15_N40_L1`, `Faner_R15_N20_L2`, `Faner_R15_N40_L2`.
- `final_faner_profile.tar.xz`: 80 members; cases `Faner_PHI_0p28000_N40_L1`, `Faner_PHI_0p28500_N40_L1`, `Faner_InitialSlope_N20_L1`, `Faner_PHI_0p20000`, `Faner_PHI_0p25000`, `Faner_PHI_0p28000`, `Faner_PHI_0p28500`, `Faner_PHI_0p30000`, `Faner_PHI_0p34000`, `Faner_PHI_0p36500`.
- `final_pd_exposure.tar.xz`: 18 members; cases `DT_N20_L1_PD_x0p3`, `DT_N20_L1_PD_x3p0`.
- `supporting_records.tar.xz`: 179 members; Supporting size-ensemble, feed-refinement, datum and calibration-generation records..

`final.json` pins all eight downloads by SHA-256 and records all 30 cases' runner arguments, switches and source variants. `ladder_report.json` is generated from these shipped records; its observable definitions are included in the report. `variants/` contains exact source replacements and their diffs. `reproduction_verification.json` records the successful fresh public-source reference run and its independent original-feed audit. Original producer audits remain inside each case archive.

To extract all archives without running a model:

```sh
python scripts/reproduce_final_reference.py --extract-only
```

The script downloads missing assets, verifies their hashes and member manifests, and refuses changing an existing extracted file. Absolute paths in producer records identify their original execution context.
