# Soybean-meal particle model

Source and reproducibility material for **A particle model of hexane removal
and water uptake in soybean meal desolventizing**, by Tomáš Svoboda, Marcela
Sluková, Tomáš Moucha and Phuong Hoang. Snapshot: **v1.1.0-submission**.

This version contains the particle solver that produced the manuscript's
convergence study and figures. It is the frozen solver runtime 151: the native
particle model under `src/` and its drivers under `paper/analysis/harness/`.
Every solver file is byte-identical to that runtime. The only exceptions are
the six public-export patches that `BUILD_RECORD.json` lists. Those patches
change path handling and the omission of internal records. They do not change
any physical law or acceptance rule.

## What changed since v1.0.0-submission

1. The temperature-dependent oil heat capacity, cp = 1024 + 3T J/(kg K), and
   the Fuller-shaped diffusivity, D_b = 4e-10 (T/378.15)^1.75 (101325/P) m²/s,
   are now explicit native laws. Earlier snapshots applied them by rewriting
   source text at run time. The retained-water law above the sorption cap
   (the film-imbibition bridge, saturated loading 0.55 kg/kg dry) is one shared
   implementation in `src/dtdc_simulator/core2/through_bed_sorbed_water_law.py`.
2. The solver controls were revised under the hygiene principle. A search
   target is never tighter than its gate. Energy rows are gated at
   max(1e-10, 2 × the row's representability floor). The birth certificate is
   evaluated on the column-scaled matrix that the Newton search solves (cap 1e8).
   The line search measures sufficient decrease in gate units. The
   film-depletion instant is located as an event. Original-feed components are
   conserved by construction.
3. A single case runner (`run_particle_qualification_case.py`) marches either
   study from its original feed on any radial mesh and time level and audits
   it. The radial-mesh and time-step ladder of the paper can therefore be rerun.

Check out `v1.0.0-submission` or `submission-2026-10-15` to rerun the earlier
configurations with their own source. The 2026-10-15 evidence archive
(`data/paper_evidence.tar.xz`) is kept unchanged for the records that the
manuscript still cites. However, the entry points that rebuilt it are not part
of this version, because they used the source-rewriting adapters.

## Environment

Use Python 3.14 and the pinned packages:

```sh
python -m pip install -r requirements-dev.txt
```

Set `PYTHONPATH=src:paper/analysis/harness` (on Windows use `;` instead of `:`),
`PYTHONHASHSEED=1`, `PYTHONUTF8=1` and
`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1`.
Run from the repository root with one numerical worker. In PowerShell, set
each variable with `$env:NAME='value'`.

## Reproduce the reference cases

```sh
python scripts/reproduce_runtime151_reference.py
```

The script verifies and extracts `data/runtime151_reference.tar.xz` to
Git-ignored paths. It then reruns the Faner and DT cases at 8 radial cells and
the middle time level (about 10 and 20 minutes on one core) and compares every
recorded value with the archived original, without a tolerance. Wall-clock
fields, output paths and source-digest maps are the only differences it
excuses. The tree here holds a subset of the producing runtime's files, so its
digest maps differ.

## Run any ladder case

```sh
python paper/analysis/harness/run_particle_qualification_case.py --study Faner --cells 16 --time-level 1 \
    --per-area-qualification --dry-local-tolerance 1e-10 --conservative-components \
    --memoize-proof-points --out reruns/Faner_N16_L1
python paper/analysis/harness/run_particle_qualification_case.py --study DT --cells 16 --time-level 1 \
    --per-area-qualification --dry-local-tolerance 1e-10 --conservative-components \
    --event-first-film-depletion --eliminate-event-absent-components \
    --curvature-corrected-depletion-prediction --memoize-proof-points \
    --bath paper/analysis/inputs/qualification_DT_bath.json --out reruns/DT_N16_L1
```

`--cells` takes 2 to 64 and `--time-level` takes 0, 1 or 2. `--saturated-loading`
sets the declared bridge parameter (default 0.55; the paper sweeps 0.25 to
0.97). `data/runtime151_ladder_records.tar.xz` holds the case, endpoint, birth
and front records of every ladder and sweep case behind the paper. It also
holds the runtime's freeze manifest with the hash of every file, the
independent original-feed audit of every history, and the ladder report with
its observed orders and GCI estimates. `data/runtime151_adopted_cases.tar.xz`
holds the complete 64-cell, middle-time-level histories that the manuscript
reports. `python scripts/reproduce_runtime151_reference.py --extract-only`
unpacks all three runtime-151 archives. A rerun can be compared with any
archived case:

```sh
python scripts/compare_runtime151_case.py reruns/DT_N16_L1 paper/analysis/results_2026-10-17/ladder_151/DT_N16_L1
```

Only the case, endpoint, birth and front records are archived for that row, so
the comparison reports the histories and checkpoints as missing.

## Layout

| Material | Location |
|---|---|
| Native particle model and its import closure | `src/` |
| Case runner, drivers, solver controls and auditors | `paper/analysis/harness/` |
| Prescribed DT bath of the qualification cases | `paper/analysis/inputs/` |
| Runtime-151 reference cases, adopted 64-cell histories, ladder records and freeze manifest | `data/runtime151_*.tar.xz`, `data/runtime151.json` |
| Earlier (2026-10-15) evidence archive | `data/paper_evidence.tar.xz`, `data/evidence.json` |
| Manuscript and supplement sources, submitted figure PDFs | `paper/01_particle_jfpe/`, `paper/shared/` |
| Build record and checksums of every committed file | `BUILD_RECORD.json`, `MANIFEST.sha256` |

## Scope

These are prescribed-bath studies of one particle. The effective pore law is
literature-informed and has not been validated. No plant validation,
calibration, production use or whole-history convergence beyond the ladder
reported in the paper is claimed. License: MIT. Digitized readings keep their
individual source and redistribution notices.
