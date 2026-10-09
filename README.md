# A particle model of hexane removal and water uptake in soybean meal desolventizing

Reproducibility repository for the paper

> Svoboda, Sluková, Henke, Moucha. *A particle model of hexane removal and water
> uptake in soybean meal desolventizing.* Journal of Food Process Engineering
> (submitted).

## What this repository is, and what it is not

It contains what is needed to

1. rebuild the paper's three computed figures and check its printed numbers against
   the records they come from;
2. re-run the two particle histories the paper presents (one meal particle through a
   desolventizer-toaster at a liquid-contact coefficient K_l = 0.32 mol m^-2 s^-1, and
   the Faner et al. (2019) laboratory sample) and the sensitivity marches tabulated in
   Supplement S12, with the particle model's own source code;
3. compile the submitted manuscript and supplement.

It is a particle-model release. The equipment (desolventizer-toaster) model is **not**
included. The particle in the desolventizer is driven by a prescribed bath (gas
temperature and hexane mole fraction along the equipment), which that equipment model
computed; the bath is shipped as a declared input file
(`paper/analysis/results_2026-10-01/dtdc_particle_destiny/destiny_v3_exit105.json`) and
cannot be regenerated from this repository. Some tray-level modules are present under
`src/dtdc_simulator/core2/` only because the particle's liquid-film kernel imports
them; the equipment solver itself is not shipped (`src/dtdc_simulator/core/` is a stub).

Every file keeps the relative path it has in the authors' working tree, so the scripts
run unchanged from the repository root. `MANIFEST.sha256` lists every file with its
SHA-256; `BUILD_RECORD.json` records the source commit the repository was cut from and
everything declared below as not reproducible here.

## Layout

| Path | Contents |
|---|---|
| `src/dtdc_simulator/core2/` | the particle model: properties (`props/`), the particle with its wet core, dry shell and external film (`particle/`), the sorption interface, and the film kernel with the modules it imports |
| `src/dtdc_simulator/core/` | stub package; only `critical_volume.py` (arithmetic re-exported by `core2.props`) |
| `benchmarks/qsc_mech01_engineering_authority_v1.yaml` | configuration read when the film kernel is imported (declared engineering assumptions for the tray shaft; no machine data; not used by the particle calculations) |
| `scripts/` | six modules of the particle model's receding-front solver, imported by the harness |
| `paper/analysis/harness/` | the scripts that produced every shipped particle history, plus their helpers |
| `scratchpad/dt_contact_sensitivity_2026_10_08/` | wrapper scripts for the regime-B segments of the desolventizer marches (numerical recovery after a rejected step: smaller first steps, seed construction, longer time budget); kept at the path the run records name |
| `scratchpad/particle_destiny_2026-10-06/faner_recovery/active_faner/` | helper modules of the Faner march, kept at the path the harness expects |
| `paper/analysis/results_2026-10-08/dt_contact_0p32/` | the presented desolventizer history: commands, segments, checkpoints, audits, outlet metrics |
| `paper/analysis/results_2026-10-08/dt_contact_sensitivity/`, `dt_mobility_f100_k0p32/`, `dt_mobility_f10000_k0p32/`, `paper/analysis/results_2026-10-07/finite_wetting/` | the sensitivity marches of Supplement S12 |
| `paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/` | the presented Faner history: step records, checkpoints, commands, the accepted history CSV |
| `paper/analysis/results_2026-10-06/`, `paper/analysis/results_2026-10-01/` | the feed march (regime A), the prescribed bath, and the earlier records the numbers check reads |
| `paper/analysis/datasets/faner2019_v2/` | the digitized Faner et al. (2019) curves (CSV/JSON with their record and checksums; no page images) |
| `docs/evidence/` | computed records that bind some printed numbers of the supplement |
| `paper/01_particle_jfpe/` | the submitted LaTeX sources, the built PDFs, `figures/`, `regen/` (figure scripts, reference check) and `check_numbers_2026-09-23.py` |
| `reruns/` | not shipped; the suggested place for your own re-run outputs (ignored by git) |

## Environment

Tested on Windows 11 with Python 3.14.7 and the pinned packages in `requirements.txt`
(numpy 2.4.6, scipy 1.18.0, matplotlib 3.11.1, PyYAML 6.0.3); the manuscript was compiled with MiKTeX. Everything runs in one Python
process with one thread.

PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\pip install -r requirements.txt
.\.venv\Scripts\Activate.ps1
$env:PYTHONPATH = (Resolve-Path src).Path
$env:PYTHONHASHSEED = '1'
$env:OPENBLAS_NUM_THREADS = '1'; $env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'; $env:NUMEXPR_NUM_THREADS = '1'
```

POSIX shell:

```sh
python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
export PYTHONPATH="$PWD/src" PYTHONHASHSEED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
```

On Windows, clone into a short folder: the longest path inside the repository has
137 characters, and Windows limits full paths to 260 characters unless long
paths are enabled.

Run the commands below from the repository root unless a `cd` is shown. The run records
keep relative paths with Windows separators; the scripts that read paths from inside a
record convert them, and all other paths come from the command line, but the
repository has been tested on Windows only.

## Rebuild the figures and check the numbers

```sh
cd paper/01_particle_jfpe/regen
python fig_particle_experiments.py --case k0p32
cd ..
```

This takes a few seconds and rewrites, in `paper/01_particle_jfpe/figures/`,
`fig_dt_particle_history_core2.pdf/.png` (the desolventizer particle),
`fig_faner_particle_history_core2.pdf/.png` (the Faner sample against the measurements)
and `fig_faner_supporting_states.pdf/.png`, together with
`fig_particle_experiments.sources.json`, which lists the SHA-256 of every input the
script read. Before drawing, the script checks that the plotted segments are the ones
the shipped audits link, that each segment starts from the checkpoint the previous one
saved, and that the plotted end point equals the recorded exit state. In the tested
run the PNG files and the sources record were byte-identical to the shipped ones and
the PDFs differed only in their creation date. `fig_particle_model.pdf` (Figure 1) is
a drawing supplied by the authors and is not generated by a script.

Then, still in `paper/01_particle_jfpe`:

```powershell
$env:DTDC_PUBLIC_PACKAGE = '1'; python check_numbers_2026-09-23.py
```

```sh
DTDC_PUBLIC_PACKAGE=1 python check_numbers_2026-09-23.py
```

The script (about half a minute) reads the submitted LaTeX sources, binds each
asserted printed quantity to the record that holds it and checks that superseded
values are no longer printed. It was kept over the whole drafting history of the paper,
so it also holds assertions whose numbers are no longer printed, and assertions whose
record is not part of this repository (the authors' working notes, and superseded
results); with `DTDC_PUBLIC_PACKAGE=1` the latter are reported as `RECORD-NOT-SHIPPED`,
not as a mismatch. Expected final counts:

```text
assertions checked against a record : 1521
superseded literals checked absent  : 233
mismatch                            : 0
awaiting landed evidence            : 0
record not shipped in the package   : 690
stale                               : 0
unrecorded                          : 0
assertions printed in submitted sources: 333
assertions only in preserved internal archive: 0
assertions no longer in any preserved source: 1165
```

Of the 333 assertions whose numbers appear in the submitted documents,
233 are checked against shipped records and 100 rest only on records
that are not shipped (listed in `BUILD_RECORD.json` under `numbers_check_in_this_cut`).
The script counts a number as printed when its digits occur anywhere in the submitted
text, so part of these are coincidences (for example a point count of 7 matching a
LaTeX setting) bound to superseded results that are deliberately left out.
Its first pass finds 1 of the 381 numbers printed in the main text without a match in the shipped records (1 declared in the script as a derived or declared value, 0 unrecorded).

`python regen/check_refs.py` checks that every cross-reference between the manuscript
and the supplement resolves.

## Re-run the two histories

Every script writes a new output file and stops if that file already exists, so write
re-runs to a new folder such as `reruns/`. The regime-B wrappers also write their
numerical-recovery records next to the output file. A re-run reads its starting
checkpoint (`*.accepted.pkl`) from the path given with `--resume`; to repeat a whole
chain, point each `--resume` at the checkpoint your previous step wrote.

### One particle through the desolventizer-toaster (K_l = 0.32 mol m^-2 s^-1)

The 1.5 mm particle (two radial cells, retained-water mobility factor 1000) enters
with the feed at t = 0. Regime A (pores full of liquid hexane) is the feed march
`paper/analysis/results_2026-10-06/shared_model/dt_feed_massfilm_gated_R1p5.json`; it
ends at 492.26 s, just after the predesolventizer exit, when the liquid-hexane core
forms. Regime B (shrinking liquid-hexane core) is the segment `b_02.json` (to
497.03 s) followed by `r3_b_01.json`, which ends with the core's extinction at
503.56 s; the attempts between them (`b_03` to `r2_b_01`) each stopped at a rejected
first step and saved the same state unchanged, and are kept because the audit links
them. Regime C (dry shell with liquid-water imbibition) is `c_dt5.json`, which reaches
the exit of the prescribed bath at 1766.54 s.

Shortest checks, all tested in this repository with outputs under `reruns/`:

```sh
python paper/analysis/harness/audit_finite_wetting_b.py --result paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --parts paper/analysis/results_2026-10-08/dt_contact_0p32/b_02.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_03.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_04.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_05.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_06.json paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_01.json paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_02.json paper/analysis/results_2026-10-08/dt_contact_0p32/r2_b_01.json paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --dry paper/analysis/results_2026-10-08/dt_contact_0p32/c_dt5.json --out reruns/dt_contact_0p32/through_exit.audit.json
python paper/analysis/harness/experiment_core2_radial_fast_wetting.py --resume paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.accepted.pkl --factor 1000 --contact-conductance 0.32 --dt 5 --out reruns/dt_contact_0p32/c_dt5.json
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032_budget.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out reruns/dt_contact_0p32/r3_b_01.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/r2_b_01.accepted.pkl
```

1. The audit (seconds) recomputes the component and energy balances of the linked
   history and writes a file byte-identical to the shipped `through_exit.audit.json`.
2. Regime C (recorded wall time 12 min) ends at 1766.54 s with particle water
   0.222554 kg/kg dry meal, hexane 135.231 mg/kg dry meal (18.2020 % wet-basis
   moisture, 110.601 ppm wet basis) and a volume-mean temperature of 104.98 °C; the
   external water film is consumed at 624.549 s. In the tested re-run (640 s) the history, the final state and the saved checkpoint state were identical to the shipped record; only provenance fields (file paths, timings, the text of the tracebacks of rejected step attempts) differed.
3. The last regime-B segment (recorded wall time 17 min) repeats `r3_b_01` from the
   checkpoint `r2_b_01.accepted.pkl`. In the tested re-run (924 s) the segment completed regime B at 503.56 s with a history, attempt records and checkpoint state identical to the shipped ones; only provenance fields (file paths, timings) differed.

The full chain, in the order it was run (the commands recorded in the
`*.command.json` files; the recorded wall time of `b_02` was 7 min). As
recorded, they write into the shipped folder; to repeat them, change that folder in
`--out`, `--resume` and `--birth-seed` to your own:

```sh
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/b_01.json
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/b_02.json --birth-seed paper/analysis/results_2026-10-08/dt_contact_0p32/b_01.birth_debug.json
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/b_03.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/b_02.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py --factor 1000 --contact-conductance 0.32 --dt 0.05 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/b_04.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/b_03.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py --factor 1000 --contact-conductance 0.32 --dt 0.025 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/b_05.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/b_04.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py --factor 1000 --contact-conductance 0.32 --dt 0.0125 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/b_06.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/b_05.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032_seed_bound.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_01.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/b_06.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032_seed_bound.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_02.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_01.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_admissible_reference_032.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/r2_b_01.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_02.accepted.pkl
python scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032_budget.py --factor 1000 --contact-conductance 0.32 --dt 0.1 --steps 300 --until 510 --out paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --resume paper/analysis/results_2026-10-08/dt_contact_0p32/r2_b_01.accepted.pkl
python paper/analysis/harness/audit_finite_wetting_b.py --result paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --parts paper/analysis/results_2026-10-08/dt_contact_0p32/b_02.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_03.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_04.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_05.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_06.json paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_01.json paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_02.json paper/analysis/results_2026-10-08/dt_contact_0p32/r2_b_01.json paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --out paper/analysis/results_2026-10-08/dt_contact_0p32/b.audit.json
python paper/analysis/harness/experiment_core2_radial_fast_wetting.py --resume paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.accepted.pkl --factor 1000 --contact-conductance 0.32 --dt 5 --out paper/analysis/results_2026-10-08/dt_contact_0p32/c_dt5.json
python paper/analysis/harness/audit_radial_fast_wetting.py --result paper/analysis/results_2026-10-08/dt_contact_0p32/c_dt5.json --out paper/analysis/results_2026-10-08/dt_contact_0p32/c.audit.json
python paper/analysis/harness/audit_finite_wetting_b.py --result paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --parts paper/analysis/results_2026-10-08/dt_contact_0p32/b_02.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_03.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_04.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_05.json paper/analysis/results_2026-10-08/dt_contact_0p32/b_06.json paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_01.json paper/analysis/results_2026-10-08/dt_contact_0p32/r1_b_02.json paper/analysis/results_2026-10-08/dt_contact_0p32/r2_b_01.json paper/analysis/results_2026-10-08/dt_contact_0p32/r3_b_01.json --dry paper/analysis/results_2026-10-08/dt_contact_0p32/c_dt5.json --out paper/analysis/results_2026-10-08/dt_contact_0p32/through_exit.audit.json
```

Regime A, the feed march (the arguments are reconstructed from the record: its time
steps of 1 s in the predesolventizer and 0.01 s after it, radius and cells; the march
stops when the liquid core forms, so any `--max-time-s` above 493 s gives the same
record; the script writes only under `paper/`):

```sh
python paper/analysis/harness/experiment_core2_dt_feed.py --out paper/analysis/reruns/dt_feed_massfilm_gated_R1p5.json --dt-s 1 --steam-dt-s 0.01 --radius-m 0.0015 --cells 2 --max-time-s 600
```

In the tested re-run (37 s) all 698 rows, the activation state and the closing step record were identical to the shipped record; only the provenance stamps differed.


### The Faner et al. (2019) sample

A 0.885 mm sphere (the reported volume-mean diameter halved; four radial cells;
initial water 0.025 kg/kg dry meal) in a 120 °C bath for 240 s. The accepted history
is the chain of step records below; `history/accepted_history.csv` and
`history/summary.json` are exported from them. The falling-rate interval (from 0.20 to
0.10 kg/kg dry meal) is 17.7 s; at 240 s the particle holds
0.00161 kg/kg dry meal of hexane and 0.0144 kg/kg dry meal of water at a
volume-mean temperature of 114.35 °C.

Shortest checks:

```sh
python paper/analysis/harness/export_core2_faner_history.py --tip paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/dry_01.json --out reruns/faner_export --ledger-tolerance 4e-10
python paper/analysis/harness/experiment_particle_dry.py --out reruns/faner/dry_01.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_24.accepted.pkl --dt 2.5 --steps 150 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out reruns/faner/birth_01.json --radius 0.000885 --water 0.025 --cells 4 --birth-dt 0.01 --coefficient-case db_max_dapp_max --continued-water --leading-seed --steps 0
```

1. The export (seconds) rebuilds `accepted_history.csv` and `summary.json` from the
   shipped step records; in the tested run both were byte-identical to the shipped
   files.
2. The last step, the dry continuation to 240 s, repeats `dry_01.json` from the
   checkpoint `event_24.accepted.pkl`. In the tested re-run (about 1 min) the record was identical to the shipped one apart from provenance fields (file paths, timings).
3. The first step, the birth of the moving front, repeats `birth_01.json`. In the tested re-run (about 15 s), made with the shipped, later revision of the harness, the record was identical to the shipped one apart from provenance fields.

The full chain in its accepted order (the commands recorded in the `*.command.json`
files, or for the steps without one, the arguments the step stored in its own record;
radial steps ran with a wall budget of 180 s each). As above, change the folder in
`--out`, `--resume` and `--tip` to repeat them:

```sh
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/birth_01.json --radius 0.000885 --water 0.025 --cells 4 --birth-dt 0.01 --coefficient-case db_max_dapp_max --continued-water --leading-seed --steps 0
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_01.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/birth_01.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.15 --horizon 54.27864018157967 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_02.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_01.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.5 --horizon 61.21355100224768 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_03.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_02.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.5 --horizon 63.79642341751772 --wall-budget 180
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_01.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_03.accepted.pkl --event arrival --dt 2.1284441491980335 --upper 10 --regional-tolerance 1e-10
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_01.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_06.json --event departure --dt 0.01 --upper 10 --budget 2400 --native-budget 6000 --wall-budget 180.0 --regional-tolerance 2e-10 --tangent-tolerance 1e-09 --event-time-tolerance 0.002
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_04.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_06.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.5 --horizon 71.66450290293075 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_05.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_04.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.5 --horizon 72.14060121152042 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_06.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_05.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.5 --horizon 72.16034788231417 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_07.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_06.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.5 --horizon 74.83336427035175 --wall-budget 180
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_07.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_07.accepted.pkl --event arrival --dt 1.7016155336531353 --upper 10 --regional-tolerance 1e-10
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_08.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_07.accepted.pkl --event departure --dt 0.01 --regional-tolerance 2e-10 --tangent-tolerance 1e-09 --native-budget 6000
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_08.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_08.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.015 --horizon 79.75538320367302 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_09.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_08.accepted.pkl --continued-water --regional-tolerance 1e-09 --difference-mode central --steps 100 --continuation-dt 0.01265625 --horizon 79.76876722206585 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_10.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_09.accepted.pkl --continued-water --regional-tolerance 1e-08 --extended-regional-root --difference-mode central --steps 100 --continuation-dt 0.018984375 --horizon 79.93590581572552 --wall-budget 180
python paper/analysis/harness/experiment_core2_faner.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_11.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_10.accepted.pkl --continued-water --regional-tolerance 1e-08 --extended-regional-root --difference-mode central --steps 100 --continuation-dt 0.028476562500000004 --horizon 80.6996945994 --wall-budget 180
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_09.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/radial_11.accepted.pkl --event arrival --dt 1.8303778553672996 --upper 10 --regional-tolerance 1e-10
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_10.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_09.accepted.pkl --event departure --dt 0.01 --regional-tolerance 2e-10 --tangent-tolerance 1e-09 --native-budget 6000
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_11.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_10.accepted.pkl --event target --target 0.001 --dt 1.291165352504791 --budget 4000 --regional-tolerance 1e-10 --condition-ceiling 1e+18
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_12.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_11.accepted.pkl --event target --target 0.0001 --dt 0.2010716996314188 --budget 4000 --regional-tolerance 1e-10 --condition-ceiling 1e+18
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_13.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_12.accepted.pkl --event target --target 1e-05 --dt 0.08480154358733259 --budget 4000 --regional-tolerance 1e-10 --condition-ceiling 1e+18
python paper/analysis/harness/experiment_particle_event.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_14.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_13.accepted.pkl --event target --target 1e-06 --dt 0.038512378376559735 --budget 4000 --regional-tolerance 1e-10 --condition-ceiling 1e+18
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_14.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_18.json --event target --target 1e-07 --dt 0.083 --upper 10 --budget 12000 --wall-budget 180.0 --regional-tolerance 1e-09 --tangent-tolerance 1e-10 --condition-ceiling 1e+18 --event-time-tolerance 0.002
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_18.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_20.json --event target --target 1e-08 --dt 0.0386 --upper 10 --budget 12000 --wall-budget 180.0 --regional-tolerance 2e-08 --tangent-tolerance 1e-10 --condition-ceiling 1e+18 --event-time-tolerance 0.002
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_20.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_21.json --event target --target 1e-09 --dt 0.0179 --upper 10 --budget 12000 --wall-budget 180.0 --regional-tolerance 2e-08 --tangent-tolerance 1e-10 --condition-ceiling 1e+18 --event-time-tolerance 0.002
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_21.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_22.json --event target --target 1e-10 --dt 0.0083 --upper 10 --budget 12000 --wall-budget 180.0 --regional-tolerance 2e-08 --tangent-tolerance 1e-10 --condition-ceiling 1e+18 --event-time-tolerance 0.002
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_22.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_23.json --event target --target 1e-11 --dt 0.00386 --upper 10 --budget 12000 --wall-budget 180.0 --regional-tolerance 2e-08 --tangent-tolerance 1e-10 --condition-ceiling 1e+18 --event-time-tolerance 0.002
python paper/analysis/harness/experiment_particle_event.py --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_23.accepted.pkl --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_24.json --event projection --dt 0.1 --upper 10 --budget 2400 --wall-budget 180.0 --regional-tolerance 2e-08 --tangent-tolerance 1e-10 --condition-ceiling 1e+18 --event-time-tolerance 0.02
python paper/analysis/harness/experiment_particle_dry.py --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/dry_01.json --resume paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/event_24.accepted.pkl --dt 2.5 --steps 150 --wall-budget 180
python paper/analysis/harness/export_core2_faner_history.py --tip paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/dry_01.json --out paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/history --ledger-tolerance 4e-10
```

## Sensitivity marches (Supplement S12)

The same particle and bath as the presented desolventizer history, with the liquid
contact coefficient K_l or the retained-water mobility factor changed. Values at the
bath exit (1766.54 s), read from the regime-C records:

| K_l (mol m^-2 s^-1) | Mobility factor | Moisture (% wet basis) | Hexane (mg/kg dry meal) | Hexane (ppm wet basis) | Film consumed (s) | Regime-C record |
|---|---|---|---|---|---|---|
| 0.05 | 1000 | 16.68 | 194.5 | 162.0 | 816.8 | `paper/analysis/results_2026-10-08/dt_contact_sensitivity/k0p05_r10_c_dt5.json` |
| 0.10 | 1000 | 17.40 | 165.3 | 136.5 | 733.6 | `paper/analysis/results_2026-10-07/finite_wetting/c_f1000_k0p1_dt5.json` |
| 0.20 | 1000 | 17.94 | 144.8 | 118.8 | 662.4 | `paper/analysis/results_2026-10-08/dt_contact_sensitivity/k0p20_r3_c_dt5.json` |
| 0.32 | 1000 | 18.20 | 135.2 | 110.6 | 624.5 | `paper/analysis/results_2026-10-08/dt_contact_0p32/c_dt5.json` |
| 0.32 | 100 | 17.62 | 158.0 | 130.1 | 702.5 | `paper/analysis/results_2026-10-08/dt_mobility_f100_k0p32/c_dt5.json` |
| 0.32 | 10000 | 18.27 | 132.5 | 108.3 | 614.0 | `paper/analysis/results_2026-10-08/dt_mobility_f10000_k0p32/c_dt5.json` |

Each case folder holds the completed regime-B segment, the segments its audit links,
their checkpoints and numerical-recovery records, the regime-C result and three
audits; the `*.command.json` files give the exact commands. The completed K_l = 0.05
segment (`k0p05_r10_b_01.json`) was started by a controller without a command file:
it resumed `k0p05_r3_b_10.accepted.pkl` through the wrapper chain named in its
`*.reference_policy.json` and `*.neighbor_seed_policy.json` records. The K_l = 0.10
case (2026-10-07) has no command files; its records hold every argument.

## Compile the manuscript

```sh
cd paper/01_particle_jfpe
build.bat        # Windows
sh build.sh      # elsewhere
```

This needs `pdflatex` and `bibtex` (TeX Live or MiKTeX) with the apacite, natbib,
siunitx and xr packages. It writes `main.pdf` (28 pages) and
`supplementary_submission.pdf` (24 pages); the two documents
cross-reference each other, which is why each is compiled after the other. In the
tested build the text of both PDFs (extracted with `pdftotext`) was identical to the
submitted PDFs shipped in the same folder.

## What is not reproducible here

* The prescribed bath of the desolventizer march
  (`paper/analysis/results_2026-10-01/dtdc_particle_destiny/destiny_v3_exit105.json`)
  was computed by the equipment model, which is not part of this repository. It is a
  declared input.
* `paper/analysis/results_2026-10-06/repairs/dt_shared_complete_trial/accepted_history.csv`
  is read by the figure script only for the rows before regime B, which are the rows of
  the feed march record. Its later rows are an earlier history of the same particle
  without liquid-water imbibition, whose producer is not included.
* The early steps of the Faner chain (birth, radial steps 1 to 9, events 1 to 8) and the
  K_l = 0.10 segments (regime B and regime C) were run with earlier revisions of harness
  scripts
  (`experiment_core2_faner.py`, `positive_water_continuation.py`,
  `experiment_particle_event.py`, and for K_l = 0.10 `experiment_core2_finite_wetting_b.py`
  and `experiment_core2_radial_fast_wetting.py`). The model source under `src/` was the
  same; the SHA-256 of every script a step used is stored in its record
  (`source_sha256_at_entry`). Re-running those steps uses the shipped revisions.
* `REPORT.md`, `outlet_metrics.json`, `outlet_comparison.csv`, `completed_case.json` and
  `provenance.json` of the K_l = 0.32 case were written by a reporting script of the
  authors that is not included; they are records.
* Three scripts differ from the authors' working copies (recorded in
  `BUILD_RECORD.json`): the figure script no longer hashes two published books whose
  outlet ranges it quotes (Witte; Kemper; not redistributable) and looks up the audit's
  Windows-style keys on any platform; the Faner export resolves the parent of each step
  inside this repository; and the loader of the tray-shaft configuration no longer
  requires the authors' note that the configuration cites, which is not distributed.
  No calculation is affected.
* The records keep the absolute paths of the authors' machine in their provenance
  fields; they are informative only.

## Citation and licence

Please cite the paper above. The code is released under the MIT licence (`LICENSE`).
The digitized Faner et al. (2019) data are readings of published figures; cite the
original article when you use them.
