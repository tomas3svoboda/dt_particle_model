XZ-compressed tar archives are kept here. Each holds an inner
`EVIDENCE_MANIFEST.sha256` that pins every member. The members keep the
unchanged bytes of the producing runs under producer-relative, Git-ignored
paths (`paper/analysis/results_v1.4.0/<case>/`). No pickled checkpoint is archived.

- `final_reference.tar.xz` (8 members): the reference case
  `Faner_R15_N20_L1` (case, endpoint, birth and front records, the front and
  dry histories, the runner's original-feed prefix audit and the independent
  whole-history audit). `scripts/reproduce_final_reference.py` reruns it.
- `final_ladder_records.tar.xz` (196 members): the JSON records
  of every DT and Faner case, their audits, and `ladder_report.json`.
- `final_adopted_dt.tar.xz` (99 members): the shipped file set (the JSON records, both JSONL histories and the audits) of the cases `DT_N10_L1`, `DT_N20_L1`, `DT_N40_L1`, `DT_N80_L1`, `DT_N20_L1_WSAT025`, `DT_N20_L1_WSAT032`, `DT_N20_L1_WSAT075`, `DT_N20_L1_WSAT097`, `DT_N20_L1_DB190`, `DT_N20_L1_DB200`, `DT_N20_L1_DB600`.
- `final_adopted_dt_time_pair.tar.xz` (18 members): the shipped file set (the JSON records, both JSONL histories and the audits) of the cases `DT_N20_L2`, `DT_N40_L2`.
- `final_adopted_faner.tar.xz` (40 members): the shipped file set (the JSON records, both JSONL histories and the audits) of the cases `Faner_R15_N10_L1`, `Faner_R15_N20_L1`, `Faner_R15_N40_L1`, `Faner_R15_N20_L2`, `Faner_R15_N40_L2`.
- `final_faner_profile.tar.xz` (80 members): the shipped file set (the JSON records, both JSONL histories and the audits) of the cases `Faner_PHI_0p28000_N40_L1`, `Faner_PHI_0p28500_N40_L1`, `Faner_InitialSlope_N20_L1`, `Faner_PHI_0p20000`, `Faner_PHI_0p25000`, `Faner_PHI_0p28000`, `Faner_PHI_0p28500`, `Faner_PHI_0p30000`, `Faner_PHI_0p34000`, `Faner_PHI_0p36500`.
- `final_pd_exposure.tar.xz` (18 members): the shipped file set (the JSON records, both JSONL histories and the audits) of the cases `DT_N20_L1_PD_x0p3`, `DT_N20_L1_PD_x3p0`.

`final.json` pins the three archives and records each case's runner arguments,
declared switches and source variant. `ladder_report.json` is a plain copy of
the report inside the records archive. `variants/` holds the declared source
variants (exact replacements and a diff). `public_tree_audits.json` and
`reproduction_record.json` record the audits and the reproduction run from this
tree for the release. To extract without running anything:

    python scripts/reproduce_final_reference.py --extract-only

A modified extracted file is refused, never overwritten. Absolute paths inside
the records are the producer's paths and are kept as provenance.
