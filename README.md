# Particle model of soybean-meal desolventizing

Reproducibility package for **A particle model of hexane removal and water
uptake in soybean meal desolventizing**, by Tomáš Svoboda, Marcela Sluková,
Tomáš Moucha and Phuong Hoang. Submission snapshot: **submission-2026-10-15**.
The paper is being prepared for submission to the *Journal of Food Process
Engineering*. This snapshot supersedes the scientific defaults in
`v1.0.0-submission`; the old tag and historical calculations remain available.

The current nominal calculation uses

- `D_b(T,P) = 4.0e-10 (T/378.15)^1.75 (101325/P)` m²/s, following Fuller's
  gas-phase temperature/pressure scaling with an assumed effective-pore
  reference. This is a literature-informed, unvalidated effective coefficient.
- `c_p,oil(T) = 1024 + 3T` J/(kg K), with T in kelvin, from García Rojas et al.
  (2013), Table 7. Energy uses its analytic integral.
- Total retained-water shell mobility factor **100** (added contribution 99),
  anchored through the Hemmingsen et al. soybean-meal uptake timescale estimate
  documented in the manuscript and Supplementary Material.
- Unchanged liquid-contact conductance `K_l = 0.32` mol/(m² s).

The prescribed DT particle exits at 1766.5355 s with **17.4618% wet-basis
moisture, 138.389 ppm wet-basis hexane, and 104.965 °C**. Faner's calculated
0.20-to-0.10 kg/kg falling interval is **17.7456 s**, versus **25.3504 s**
measured. These are feasibility comparisons in prescribed baths. They are
not plant validation, fitted transport parameters, or a claim of established
whole-history spatial/time convergence.

## Reproduce the submission checks

Use Python 3.14 (the tested environment) and one numerical worker. From the
repository root:

```sh
python -m venv .venv
# Activate .venv using your platform's standard command.
python -m pip install -r requirements-dev.txt
```

Set `PYTHONPATH=src`, `PYTHONHASHSEED=1`, `PYTHONUTF8=1`, and set
`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and
`NUMEXPR_NUM_THREADS` to `1`. In PowerShell, use `$env:NAME='value'`.

```sh
python scripts/reproduce_submission_2026_10_15.py
python scripts/reproduce_submission_2026_10_15.py --render
python -m pytest tests/test_paper1_temperature_properties.py -q
python paper/01_particle_jfpe/regen/check_refs.py
```

The first command checks every manifest entry and independently re-evaluates
the complete native DT water/hexane/energy balances for five cases. It also
re-exports all three current Faner histories, checking parent SHA-256 links,
constitutive inventories and exterior-flow ledgers against the archived
results. The second regenerates the figure arrays, number macros and
graphical abstract after those checks. Outputs from verification go under
`reruns/`; rendering intentionally updates generated manuscript assets.
PDF creation timestamps can change on rendering, so the original manifest
describes the released files before regeneration.

## Run a new nominal particle history

```sh
python paper/analysis/results_2026-10-15/dt_factor100/run_case.py --Db-ref 4e-10 --out-dir reruns/nominal_F100
```

This uses the existing native birth/departure controller and refuses outputs
outside `reruns/`. References `2e-10` and `6e-10`, or `--constant-D` with
`4e-10`, reproduce the corresponding diffusivity comparisons. The heat-capacity
law remains active in all these cases. Matching-law restart JSON from a new public run can be
passed with `--resume-tip`; a matching birth diagnostic can be passed with
`--birth-seed`. Refused attempts and native recovery choices remain recorded.

For the C-step comparison, run the property driver with the accepted nominal
B checkpoint from your new public run, native `experiment_core2_radial_fast_wetting.py`, `--factor 100`,
`--contact-conductance .32`, `--dt 2.5`, `--until 1766.5355008060146`,
`--steps 1500` and `--max-halvings 12`, writing a fresh output under `reruns/`.
The complete executed arguments are retained with the published cases.

## What is where

| Material | Location |
|---|---|
| Active particle source and required import closure | `src/dtdc_simulator/` |
| Property laws and native drivers/auditors | `paper/analysis/harness/` |
| Current portable study index and recomputed public DT audit holders | `paper/analysis/submission_2026-10-15/` |
| Current F100 DT raw histories, checkpoints and refused attempts | `paper/analysis/results_2026-10-15/dt_factor100/` |
| Current Faner property-study histories | `paper/analysis/temperature_property_runs/faner_Db*_Fuller/` |
| Manuscript, supplement, editable equation tables, figures | `paper/01_particle_jfpe/` |
| Measured readings and their source/reading uncertainty records | `paper/analysis/` and `benchmarks/` |
| Export provenance and per-file checksums | `BUILD_RECORD.json`, `MANIFEST.sha256` |

Compile the main manuscript and supplement from `paper/01_particle_jfpe/`
using its existing `build.sh` or `build.bat` with a full LaTeX/BibTeX
installation. The source package uses the `paper/shared/` files and mutual
cross-document references; retain this layout and complete the build passes.

Earlier constant-property sensitivity campaigns under the 2026-10-06 to
2026-10-08 result folders are historical, as identified in the supplement.
They are not the new nominal temperature-property/F100 calculation. To rerun
those historical campaigns, check out `v1.0.0-submission`; the current root
feed is the temperature-dependent-heat-capacity feed. In
particular, factor 1000 and 111 ppm in historical records are not the current
paper's nominal assumptions or prediction.

Raw result/checkpoint bytes and their original producer-path metadata are
preserved. `release_paths.py` resolves those paths inside this repository;
the authors' private workspace is not needed. Public index/audit/figure
holders are derived and explicitly distinguished from raw records. Existing
curated export patches, and the new path/portability changes, are documented
with hashes in `BUILD_RECORD.json`. Native acceptance tolerances and physical
laws are unchanged. Publisher PDFs, coauthor annotations, working notes and
the full equipment solver are excluded. The prescribed equipment-bath output
is included. A preservation-archive DOI remains pending; none is implied by
the GitHub snapshot.

Source-code license: MIT; see `LICENSE`. Digitized literature readings retain
their individual source and redistribution notices.
