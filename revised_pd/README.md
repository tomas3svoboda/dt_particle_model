# Revised PD prescribed-bath particle reference

Version **v1.5.0-pd-reference**, 20 November 2026. This is a complete,
independently conservation-audited **64-cell L0 reference**, with full-history
mesh/time refinement still open. It is a prerelease numerical reference,
not a validated operating tower or a newly converged DT campaign.

The 49 C cold feed receives constant indirect heat, 190.8872106 W/kg dry,
for 490.202242 s on solid heated PD disks. Sparge steam bypasses the PD meal.
The prescribed source gas is 63.901310 C, 0.863186905 mole fraction hexane,
at 101226.9335 Pa. Its water comes from the source meal's retained moisture.
The source is equilibrium/energy-derived, not measured local gas. Jacket duty
2.30 MW is an inherited engineering seed, not an independently measured duty.
Continuous vapor withdrawal is represented by the prescribed reservoir;
header mixing and washer suction are not dynamically solved.

The particle retains the canonical uniform thermal regime A and native radial
B/C equations. Native signed exchange allows small initial hexane uptake;
external hexane must deplete before core recession. PD boiling starts near
383.9 s. At 1766.535501 s the DT/toaster outlet is 132.756052 ppm wet-basis
hexane, 18.385533% wet-basis moisture and 104.810216 C. This is not the outlet
of the downstream dryer/cooler. Temperature is 0.19 K below the cited 105 C
lower reference edge. Moisture includes external film water while present.

## What is distributed

`runtime/` is an isolated particle import closure, byte-identical to the
original producer's pinned sources. The older public runtime and recipes
remain available separately. `SOURCE_MANIFEST.json` checks every selected
source against the original full-tree pin in `data/runtime_source_sha256.json`.
Omitted unrelated development modules are not needed by the particle closure.
Fresh reruns retain all native source-identity and conservation checks against
their own distributed-source generation. Original full-tree records are not
silently re-keyed or claimed byte-identical to a fresh public rerun.

`data/` contains operative bath/source values, accepted-state and flux CSVs,
reported values, independent conservation audit and restart identity records.
The CSV corrects the producer's stale pre-override bath label at t=0 to the
operative PD gas. All particle samples are unchanged; original native records
retain their original metadata. See `data/entry_bath_correction.json`.
The release's `revised_pd_native_records.zip` contains unmodified native
profiles, checkpoints and producer metadata; its member hashes and archive
identity are in `NATIVE_ARCHIVE.json`. Absolute historical path strings are
metadata, not prerequisites for fresh runs. The recorded completion reused an
unchanged checked cold-feed/B prefix after a RAM-reserve stop, with native
extinction and five C overlap states recomputed bit-identically.

No manuscript, coauthor artwork or publisher-owned images are published here.
The source license is MIT; the new numerical outputs are original model data.

## Reproduce from the public repository root

Use Python 3.14 and the repository's pinned `requirements.txt`.

```sh
python scripts/reproduce_revised_pd.py --verify
python scripts/reproduce_revised_pd.py --verify-native
python scripts/reproduce_revised_pd.py --probe
python scripts/reproduce_revised_pd.py --run
```

The first command checks identities and reports the archived outlet without
running a solver. The second downloads and verifies native records without
unpickling checkpoints. The probe advances only 0.1 s and is not a full
reproduction. `--run` starts a fresh original-feed history, may take hours,
retains native refusals and reports outlet differences against the record.
Existing output directories are refused. Touch the printed `.STOP` path for
cooperative stopping. Windows uses the workstation resource guard (half CPU
affinity, BelowNormal priority, one numerical thread, 35% and 25 GiB free RAM
reserve). Linux uses half the available CPUs and a lower priority; its memory
reserve is not monitored by that Windows guard.

## Remaining scientific limits

The 0.005-s birth attempt refuses; the successful 0.01-s L0 birth does not
prove time convergence. Old-bath refinement/sensitivity tables cannot qualify
this revised bath. Actual-density PD Schmidt number is near 0.35, outside
the film-correlation applicability guard; reused reference film properties
remain an unqualified transfer assumption. The draft, exposure, constant
heat allocation and local gas closure remain prescribed engineering choices.
Conservation and agreement with outlet bands do not validate those choices.
Release preparation checks source/data identity and a short probe, not a new
fresh complete-history rerun or full-history convergence.
