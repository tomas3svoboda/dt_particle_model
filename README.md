# DT particle model

Particle model of hexane removal and water uptake in soybean meal. MIT license.

## Current prescribed PD/DT experiment

Use [`revised_pd/`](revised_pd/README.md) and
[v1.5.0-pd-reference](https://github.com/tomas3svoboda/dt_particle_model/releases/tag/v1.5.0-pd-reference).
Solid heated PD disks receive indirect heat; sparge bypasses PD meal.
Vapor withdrawal is represented by a prescribed reservoir; main-tray steam
contact follows PD. The recorded reference uses uniform thermal regime A
and 64 base radial cells in B/C.

```sh
python -m pip install -r requirements.txt
python scripts/reproduce_revised_pd.py --verify
python scripts/reproduce_revised_pd.py --verify-native
```

Verification checks identities and every CSV wet-basis value without a solver.
`--probe` advances 0.1 s; `--run` requests a fresh complete history and may
take hours. Read the reference's resource limits before starting one.

The recorded DT/toaster outlet is 132.756 ppm wet-basis hexane, 18.3855%
moisture and 104.8102 C. Full-history refinement and PD film applicability
remain open; this deposit is a numerical reference, with the scientific
limits specified in its accompanying documentation.

## Faner experiment

The soybean experiment in superheated hexane retains the 1.5-mm radius and
all-point external-transfer factor 0.29. Temperature was assessed separately.
Its original case values, source variants and refinement evidence are retained
in the now **Faner-only** `data/final.json`.

```sh
python scripts/reproduce_final_reference.py --verify-manifest
python scripts/reproduce_final_reference.py --case Faner_R15_N40_L1
python scripts/regenerate_faner_results.py
```

The second command runs a model. `--extract-only` now obtains only the three
Faner archive groups. Source and independent conservation checks remain active.

## Superseded PD experiment

The previous steam-rich PD bath is **withdrawn from active use**. Its DT
recipes, PD/material sensitivities, old bath file, figure CSV and DT figures
have been removed from the main checkout. Old DT case names are explicitly
refused by the public reproduction entry point.

Historical tags and archive bytes remain for traceability, with a supersession
notice on the old release. Their PD results are not current article results or
qualification evidence for this revised bath. Valid Faner evidence is retained.

Shared particle operators and inert legacy waypoints in the new pinned runtime
remain where exact reproduction requires them. The new driver always applies
its declared PD override; these are not an active old-bath recipe. Manuscripts,
coauthor artwork and publisher source figures are not distributed here.
