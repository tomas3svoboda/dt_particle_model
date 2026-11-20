# Revised PD reference data

This folder accompanies the 20 November 2026 manuscript revision. It is a
continued, audited history from the original 49 C feed, with uniform regime A
and 64 base radial cells in B/C. It is not a new validation dataset.

- `bath.json` contains the inherited waypoint file; its pre-PD-end waypoints
  are overridden by `pd_source.json` and the constant indirect heat recipe.
- `accepted_particle_history.csv` contains the operative gas and accepted
  particle states, including external water, with no fabricated tail.
- `DT_regime_C_fluxes.csv` reconstructs the original upwind motion/diffusion
  decomposition from native recorded Stefan recovery and boundary transfers.
- `reported_values.json` records the exact values behind the rounded text.
- `independent_audit.json`, `restart_identity_audit.json` and
  `prefix_restart.json` distinguish conservation from restart provenance.
- `runtime_source_sha256.json` identifies the exact producer source generation;
  hashes alone are not a replacement for a published matching source snapshot.

The underlying native profiles and checkpoints are included in the versioned
`revised_pd_native_records.zip` release asset; see `../NATIVE_ARCHIVE.json`.
Their original research location was
`scratchpad/cold_start_2026_11_18/paper1_experiment2_complete_n64_l0_c_restart_v2_2026_11_19`.
The preceding cold-feed/PD/B prefix is in the sibling
`paper1_experiment2_complete_n64_l0_2026_11_19` directory. Its physical source
generation is unchanged; extinction and five accepted C overlap states were
reproduced bit-identically after a workstation RAM-reserve stop.

The fresh-case recipe is `run_pd_prescribed_particle_case.py` with
`--pd-heating constant --cells 64 --time-level 0 --water-film-energy body`,
a new output directory and cooperative stop file, under the repository's
one-worker resource guard. The checked-prefix extinction restart is recorded
by `continue_complete_particle_c_2026_11_19.py`. The isolated producer source
closure is distributed in `../runtime`; use
`python scripts/reproduce_revised_pd.py --verify` from the public repository
root to inspect it, or `--run` to request a fresh complete calculation.
Release v1.5.0-pd-reference adds this case separately from v1.4.0-submission.

Full-history mesh/time refinement remains open. The inlet and PD-exit external
water are zero; a small transient PD film from the moisture-bearing source is
reported rather than suppressed, and main-tray condensation supplies the large
film at birth. Pressure, homogenized jacket allocation and reference gas-film
properties are explicitly qualified engineering assumptions in the supplement.
