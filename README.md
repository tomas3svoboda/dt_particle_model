# DT particle model

Particle-level model and numerical records for *A particle model of hexane
removal and water uptake in soybean meal desolventizing*.
Release: **v1.4.0-submission**. Software license: MIT.

The formulation follows external hexane depletion, a receding liquid-hexane
core and the subsequent drying tail, with coupled water uptake and energy
conservation. Both reported histories use a 1.5-mm-radius sphere and the
same surface-graded mesh family. Oil heat capacity and binary pore
diffusivity depend on temperature. The Faner external-transfer factor
**0.29** is estimated from all 25 solvent-loading observations;
temperature is assessed separately. Thermophysical and intraparticle
transport coefficients retain their literature-derived values.

## Reproduce the results

Use Python with `requirements.txt` installed. Run from this repository's root:

```sh
python -m pip install -r requirements.txt
python scripts/reproduce_final_reference.py --verify-manifest
python scripts/reproduce_final_reference.py
python scripts/reproduce_final_reference.py --case Faner_R15_N40_L1
python scripts/reproduce_final_reference.py --case DT_N80_L1
```

The default reruns the 20-cell Faner reference. Each command downloads the
needed history archive, verifies its SHA-256 and member manifest, runs the
declared case with single-threaded BLAS, and compares its physical results
exactly with the archived producer. Paths, source identities, execution times
and explicitly declared option stamps are handled separately by the comparison.
Existing results are preserved; a reused output directory is refused.

Case inputs, time caps, parameter alternatives and option settings are in
`data/final.json`. Exact radius, exposure and diffusivity variants are
applied to a scratch source copy under `reruns/`. Numerical gates and
transport equations retain their declared settings. The DT exposure
sensitivity includes complete 0.3-, 1- and 3-times histories. The Faner
profile retains every loading observation. No manuscript text or publisher
figure images are distributed here.

The Git tree contains code, recipes and derived figure data. **Bulk histories
are [release downloads](https://github.com/tomas3svoboda/dt_particle_model/releases/tag/v1.4.0-submission)**
and are fetched only when requested. To verify and unpack the full record:

```sh
python scripts/reproduce_final_reference.py --extract-only
python scripts/regenerate_faner_results.py
```

The second command regenerates the 40-cell Faner numbers, pointwise loading
and temperature scores, history CSV and plots. Add `--case <folder>` to score
another extracted or freshly rerun Faner case. DT figure data are under
`paper/analysis/manuscript_numbers/`.

## Numerical evidence and provenance

`data/ladder_report.json`, the case records and the independent audits record
spatial and temporal refinement, material sensitivities and whole-history
water, hexane and energy closure. The original producer pinned the authors'
full source tree; this repository ships the particle import closure. Its
auditors check fresh reruns against their own source identities. Original
audits remain in the archives, separately from the public-rerun verification
in `data/reproduction_verification.json`. Source checks remain active.

`BUILD_RECORD.json` describes the full reproduction bundle; `BUNDLE_LAYOUT.json`
records the separation of bulk archives from the lean Git tree. Each download
is pinned in `data/final.json`. The complete bundle is also prepared for
versioned preservation. See `data/README.md` for the case inventory.

These are particle calculations under prescribed gas histories. Industrial
outlet ranges provide process context; the equipment model is separate.
Third-party publications and digitized readings retain their own terms.
