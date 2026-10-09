# Exploratory DT march at K_l=0.32

Date: 2026-10-08. Completed and independently audited through the prescribed bath exit. Saved for review; no adoption into manuscript or primary figures.

| K_l, mol/(m² s) | Moisture, % wet basis | Hexane, mg/kg dry | Hexane, ppm wet basis | Mean T, °C | Film depletion, s |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 17.40 | 165.3 | 136.5 | 105.03 | 733.6 |
| 0.20 | 17.94 | 144.8 | 118.8 | 104.99 | 662.4 |
| 0.32 | 18.20 | 135.2 | 110.6 | 104.98 | 624.5 |

Common outlet time: 1766.535501 s. External water is zero at the outlet. Moisture and hexane wet-basis values exclude external liquid water.

Owner-selected literature construction: 80 C coarse solvent-extracted soybean meal, first plateau sample at 90 s, initial water 0.118/(1-0.118) kg/kg dry, assumed sphere radius 0.715 mm and endpoint W_cap-0.01. Exact native cap-to-retained activity force gives K_l=0.253170 at 80 C; assumed inverse-viscosity scaling gives 0.317809 at 100 C, rounded to 0.32. Endpoint tolerance 0.02/0.005 gives conditional scenarios 0.211/0.428. This is a one-sided construction under declared assumptions; the filter-cake method does not separately measure retained interior water or identify contact and internal transport. No benchmark outlet was used to select this value. Full arithmetic and source hashes are in nonlinear_80C_construction.json.

Relative to 0.10: moisture change +0.799 percentage points, hexane change -25.87 ppm.
Relative to 0.20: moisture change +0.262 percentage points, hexane change -8.24 ppm.

Fixed: radius 1.5 mm, two radial cells, retained-water mobility factor 1000, binary diffusivity 6e-10 m²/s, native thermodynamics and capacities, actual-activity force, feed and bath. All 263 native source identities match the completed reference and nearby sensitivity.

Numerical controls: new B ordinary dt=0.1 s (see segment commands for any recovery reductions), C maximum dt=5 s and existing adaptive/event/terminal policies. Reused reference B dt=0.2 s; nearby 0.20 B dt=0.1 s. Component/energy acceptance 1e-10 unchanged. Same-IVP birth iterates may serve as numerical guesses; accepted material from another coefficient is never imported. Native face-departure retries and local terminal projection are documented in policy and attempt records. At the first face, bounded reconstruction had residual 2.189e-10 against the earlier 2e-10 seed limit; the explicitly recorded 3e-10 tangent limit is used only to construct a numerical guess. The native corrector, accepted-step controls and 1e-10 W/H/U checks remain unchanged. No physical-law changes.

Native tangent guesses also refused a phase-safe central Jacobian probe. The existing low-conductance recovery was reused: the completed K_l=0.10 departure supplies an uncommitted numerical initial guess, the actual before-state and all residuals remain those of K_l=0.32, and numerical chart/row normalization initializes from the admissible guess while retaining tangent physical column scales. See neighbor-seed and reference initialization audits. No accepted trajectory row or material state is copied from 0.10.

The harness's 120 s per-step native proof wall budget expired during departure; the recovery declares 600 s per step, retaining the 3,000,000 proof-point budget and all native accepted-state tolerances. A computational wall allowance changes solve effort and does not change the physical equations or relaxation time.

Independent linked W/H/U maximum normalized defects: 4.2385e-11, 8.6132e-13, 6.1624e-12.

This is a conditional result at a declared exploratory coefficient. It does not identify the coefficient, the internal mobility, or whole-history time/mesh error. All manuscript source/PDF and figure pins remain identical to their pre-run hashes.

Artifacts: `completed_case.json`, `outlet_metrics.json`, `outlet_comparison.csv`, `accepted_trajectories.csv`, `comparison.pdf/png`, `through_exit.audit.json`, `provenance.json`, `source_snapshot/MANIFEST.json`.
