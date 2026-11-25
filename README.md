# Soybean-meal particle model

Code, inputs and numerical records supporting *A particle model of hexane
removal and water uptake in soybean meal desolventizing*.

The two experiments are Faner's superheated-hexane experiment and the
prescribed desolventizer-toaster (DT) particle history. The DT prescription
uses indirect PD heating followed by main-tray steam contact. Sources,
assumptions and equations are specified in the article and supplement.

## Install

Use Python 3.14 and run commands from this repository's root:

```sh
python -m pip install -r requirements.txt
```

The optional compiled equation-of-state evaluator uses
`python -m pip install -r requirements-compiled.txt`; select it with
`--compiled-eos`.

## Verify the supplied data

```sh
python scripts/reproduce_final_reference.py --verify-manifest
python scripts/reproduce_matched_pd.py --verify
```

These commands check file identities and the recorded data. The DT check
also recomputes every CSV wet-basis solvent and moisture value. The native
histories and independent conservation-audit results accompany the release.

## Reproduce the calculations

```sh
python scripts/reproduce_final_reference.py --case Faner_R15_N40_L1
python scripts/reproduce_matched_pd.py --run --case fresh_n64_a16_l1
```

Run one case at a time. Each command preserves its original inputs, numerical
controls and independent conservation checks. Fresh outputs are written to
ignored reproduction folders; existing runs are preserved. Calculations can
take hours. Windows runs use the supplied CPU/memory guard and single-threaded
numerical libraries.

The Faner case catalog is `data/final.json`. The DT reference and its spatial,
time-step and diffusion-geometry comparisons are listed in
`matched_pd/data/cases.json`; select a member with `--case`.

## Inspect histories and regenerate figures

```sh
python scripts/reproduce_final_reference.py --extract-only
python scripts/regenerate_faner_results.py
python scripts/reproduce_matched_pd.py --extract-native
python scripts/plot_matched_history.py --case matched_pd/runtime/scratchpad/archived_reference --bath matched_pd/runtime/inputs/publication_bath.json --out reruns/dt_figure
```

Downloaded archives and their members are checksum-verified. The figure
scripts use accepted particle states. The DT plots contain solvent, total
moisture, temperature, prescribed-bath overlays, core volume and external
water. Numerical comparisons are supplied in
`matched_pd/data/qualification_comparison.json`.

`CITATION.cff` identifies this release. Code is distributed under the MIT
license; experimental data retain the attribution recorded with their inputs.

Auxiliary supplementary tables: `python scripts/reproduce_supplementary_tables.py`.
Numerical inputs and stored curves are in `data/supplementary/`.
