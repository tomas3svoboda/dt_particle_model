# Prescribed DT particle history

`runtime/inputs/publication_bath.json` is the operative prescription.
`data/cases.json` lists the reference and its mesh, time-step and gas-resistance comparisons.
`data/accepted_particle_history.csv` contains the reference's accepted particle states;
`data/qualification_comparison.json` gives the numerical comparisons and their definitions.
Native solver records are available through `NATIVE_ARCHIVE.json`.
Use `scripts/reproduce_matched_pd.py` from the repository root to verify, extract or rerun a case.
