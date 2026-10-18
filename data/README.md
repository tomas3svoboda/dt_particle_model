Four XZ-compressed tar archives are kept here. Each holds an inner
`EVIDENCE_MANIFEST.sha256` that pins every member. The members keep the
unchanged bytes and producer-relative paths of the original runs.

- `runtime151_reference.tar.xz` holds the two reference cases of the frozen
  solver runtime 151 (Faner and DT, 8 radial cells, middle time level): feed,
  birth, front and dry-stage records, histories and checkpoints.
  `scripts/reproduce_runtime151_reference.py` extracts it and reruns both cases.
- `runtime151_adopted_cases.tar.xz` holds the complete case trees of the rows
  that the manuscript adopts (Faner and DT, 64 radial cells, middle time level).
- `runtime151_ladder_records.tar.xz` holds the case, endpoint, birth and front
  records of every radial-mesh and time-step ladder case, and of every
  saturated-loading sweep case. It also holds the runtime's provenance:
  the freeze manifest `execution_freeze_151.json` with the SHA-256 of every
  file of the runtime, the launch and audit instruments, the independent
  original-feed audit of every history, the ladder report with observed orders
  and GCI estimates, the windowed history analysis, the saturated-loading score,
  and the exported figure histories and paper numbers.
- `paper_evidence.tar.xz` holds the 2026-10-15 evidence of the earlier
  submission snapshot. It is unchanged, and `evidence.json` pins it.

`runtime151.json` pins the three runtime-151 archives. To extract them
without running anything numerical:

    python scripts/reproduce_runtime151_reference.py --extract-only
    python paper/analysis/harness/submission_evidence.py    # the 2026-10-15 archive

Extracted paths are ignored by Git. A modified extracted file is refused,
never overwritten. Absolute paths inside the records are the producer's paths
and are kept as provenance.
