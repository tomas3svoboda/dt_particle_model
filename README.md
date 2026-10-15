# Soybean-meal particle model

Source and reproducibility evidence for **A particle model of hexane removal
and water uptake in soybean meal desolventizing**, by Tomáš Svoboda, Marcela
Sluková, Tomáš Moucha and Phuong Hoang. Submission snapshot: **submission-2026-10-15**.

The nominal model uses `D_b = 4e-10 (T/378.15)^1.75 (101325/P)` m²/s,
`c_p,oil = 1024 + 3T` J/(kg K), total liquid-water mobility factor **100**,
and unchanged `K_l = 0.32` mol/(m² s). T is in kelvin. The prescribed DT
particle exits at **17.4618% moisture and 138.389 ppm hexane, wet basis**.
The Faner falling interval is **17.7456 s**, versus **25.3504 s** measured.
The effective pore law remains literature-informed and unvalidated; these
are prescribed-bath feasibility studies, not plant validation or a claim of
whole-history mesh/time convergence. References and derivations are in the paper.

Use Python 3.14 and the pinned environment:

```sh
python -m pip install -r requirements-dev.txt
python scripts/reproduce_submission_2026_10_15.py
python scripts/reproduce_submission_2026_10_15.py --render
python -m pytest tests/test_paper1_temperature_properties.py -q
python paper/01_particle_jfpe/regen/check_refs.py
```

Set `PYTHONPATH=src`, `PYTHONUTF8=1`, `PYTHONHASHSEED=1`, and
`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1`.
Use one numerical worker. In PowerShell, set each variable with `$env:NAME='value'`.

The first command verifies and extracts **one compressed evidence archive**
to Git-ignored paths, then independently checks all five current DT histories
and all three current Faner lineages and exterior-flow ledgers. The second
regenerates figures, numeric macros and the graphical abstract. Raw run and
checkpoint bytes retain their SHA-256 identities; modified extracted evidence
is refused rather than overwritten. PDF timestamps can change on regeneration.

For a new nominal solve:

```sh
python scripts/run_particle_dt.py --Db-ref 4e-10 --out-dir reruns/nominal_F100
```

Use `2e-10` or `6e-10` for the prescribed reference sensitivity, or add
`--constant-D` for the constant-diffusivity comparison. Oil cp(T) stays active.
The source uses the existing native birth/departure recovery controller;
fresh run records are written under `reruns/`.

| Material | Location |
|---|---|
| Particle code and required import closure | `src/` |
| Property laws, native drivers and auditors | `paper/analysis/harness/` |
| Active manuscript/supplement sources and submitted figure PDFs | `paper/01_particle_jfpe/` |
| Shared preamble and cited bibliography | `paper/shared/` |
| Compressed raw evidence and inner checksum inventory | `data/` |
| Export notes and committed-file checksums | `BUILD_RECORD.json`, `MANIFEST.sha256` |

Compile from `paper/01_particle_jfpe/` with its `build.sh` or `build.bat`
and a full LaTeX/BibTeX installation. Tables remain editable LaTeX. Historical
constant-property operator and sensitivity records are distinguished from the
current nominal study inside the archive. Use the [previous public snapshot](https://github.com/tomas3svoboda/dt_particle_model/tree/c1c5ba98b56e59ca184029b7a91ef271337512b3)
to rerun earlier configurations with their original source and feed.

Only source, submitted figure PDFs and one compressed evidence bundle are
committed. Compiled manuscripts, duplicate image formats, editorial letters,
review annotations and development logs belong to the submission workspace.
The original producer paths in raw records are provenance; public readers
resolve them locally. A preservation DOI is pending. License: MIT; digitized
readings retain their individual source and redistribution notices.
