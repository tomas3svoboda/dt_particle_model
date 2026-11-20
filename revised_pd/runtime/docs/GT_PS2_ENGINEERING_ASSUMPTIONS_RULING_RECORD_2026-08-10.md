# GT-PS-2 engineering-assumptions ruling record (owner amendment)

**Date:** 2026-08-10.
**Authority:** explicit owner amendment, ruled in-session by the owner
("I am afraid we need to work with what we have. Fill the gaps with
defensible engineering assumptions. ... I accept your recommendation as
authorized rulings" — D0, D1, D2 initially; D2 amended and D3, D4, D5,
D6, D8 ruled later the same session; D7 ruled after a dynamic-fidelity
clarification; the owner corrected an initial D3/D2 transposition
in-session). ALL NINE ITEMS D0-D8 ARE RULED. This record is the durable
register of that amendment. It supplements — and does not
reopen — the frozen PHY set (53/53) and the PR-01..PR-11 production rulings.

**Evidence basis:** docs/DT_DTDC_evidence_register_2026-08-10.md (e9d22a9),
Addendum A (a2b92e9), Addendum B (44bf5c1), and the frozen uncertainty
register in docs/GT_PS2_PRODUCTION_RULINGS_CLAUDE_RECOMMENDATIONS_2026-08-09.md.

---

## D0 — Assumption regime and claim ceiling — **RULED 2026-08-10**

1. A new typed provenance class **`DECLARED_ENGINEERING_ASSUMPTION`** is
   created: owner-adopted, source-cited, bracketed, mandatorily
   sensitivity-swept, and permanently distinguishable in code, manifests and
   evidence from `AUTHENTICATED_MACHINE_DATUM`. Assumptions may fill
   **parameter** gaps only; they may never substitute for an unqualified
   physical law, weaken a fail-closed boundary, or silently overwrite a
   frozen PHY selection.
2. The terminal validation gate splits: **F8a "design-basis qualified"** —
   reachable with declared assumptions, every assumption swept and its KPI
   sensitivity published; **F8b "plant-validated"** — remains false until
   authenticated plant data exists. No output may describe F8a results as
   plant-validated.
3. Any later authenticated machine datum immediately supersedes the
   corresponding assumption; the swap is a narrow, recorded commit.

## D1 — Common-shaft speed (PR-08 envelope) — **RULED 2026-08-10**

- Nominal **11 r/min**; operating envelope **8-12 r/min**; design envelope
  **4-15 r/min**; low-speed alarm 8 r/min, trip 4 r/min (Zuo 2024 §5.3
  convention).
- Scaling: RPM treated as approximately diameter-independent across the
  class (0.9 m pilot at 15 r/min; 4.57 m Crown drive output 11 r/min);
  6.0 m tip speed at nominal ≈ 3.46 m/s (derived arithmetic, disclosed).
- Legacy 3 / 33 / 60 r/min and the GUI 0-10 range are **retired as
  unsupported**; the 14-42 r/min figure remains excluded (discharge-valve
  speed).
- Sources: Fu 2024 (primary, 10-12 r/min); Zuo 2024 (corroboration +
  alarms); S4M gearbox output 11 r/min; fleet record 8-12 provisional.
- Consequence: `f_rev` nominal = 11/60; PR-09 exposure fixtures and future
  production carriers use this envelope. Provenance class:
  DECLARED_ENGINEERING_ASSUMPTION (class evidence, not machine datum).

## D2 — Shaft-power magnitude (PR-11) — **RULED 2026-08-10, AMENDED 2026-08-10**

- **Amended after owner challenge and bottom-up granular-torque estimate:**
  nominal **50 kW** post-gearbox shaft input; bracket **25-165 kW**; both
  endpoints propagated through every energy result. (Original 120 kW
  nominal retired; retained inside the bracket as the conservative
  0.5x-installed scenario.)
- Amended basis (three converging routes, disclosed arithmetic):
  (1) granular-drag estimate — passive-state frontal stress K*rho*g*z ≈
  12-18 kPa on ~0.4-0.8 m² effective blade area per deep tray at ~2 m
  radius gives 15-30 kN·m/deep tray → 20-35 kW/deep tray at 11 r/min,
  shallow predesolv trays negligible (stress ∝ depth) → machine total
  ~30-60 kW; (2) the only operating record in the corpus (CN205035341U,
  50-55 A) reads ~27-50 kW under plausible LV voltages; (3) Fu 2024
  states main-shaft power is sized with margin for loaded start —
  installed 250 kW >> running load, load factor 0.2-0.4, not 0.5.
  Bracket floor = current-record reading; bracket ceiling = 164 kW
  gearbox continuous rating (hard physical cap).
- Physics note (recorded from the owner discussion): granular drag is
  Coulombic and rate-independent — torque is set by weight-loaded
  friction, not shear rate; power scales linearly with shaft speed.
- Materiality at amended nominal: ~+1.0 K equivalent outlet contribution
  (bracket +0.5 to +3.4 K) at 48.7 kW/K; in bubble-point-pinned wet trays
  the power appears as ~0.15 kg/s additional vapor (<1 % of steam duty),
  surfacing as temperature only where the latent sink is exhausted
  (toasting outlet, DC).
- Omission logic unchanged: a DECLARED_ENGINEERING_ASSUMPTION can never
  constitute the source-pinned full-envelope bound below 65 kW that PR-11
  requires, and the bracket ceiling far exceeds it — **the shaft-work
  ledger row remains LIVE**.
- Provenance class: DECLARED_ENGINEERING_ASSUMPTION; superseded instantly
  by any authenticated nameplate/torque/VFD record.

## D3 — Dissipation partition (PR-11) — **RULED 2026-08-10**
*(ruled after the owner-led physics discussion; the tray energy solve —
not the ledger row — determines the per-tray sensible/latent split, so wet
bubble-point-pinned trays convert shaft heat to vapor, not temperature)*

- Nominal partition: **all-to-meal (100 %)**. Physical basis: dissipation
  occurs at the shear planes in the solvent-wet bed (meal-on-meal,
  meal-on-floor; blade floor clearance 20-25 mm; no wall contact);
  Schumacher US4622760A locates the dominant load in the wet region;
  blade-geometry changes move motor current (CN205035341U), i.e. the load
  is bed rheology. Bearing/gearbox losses are outside the PR-11
  post-gearbox boundary by definition.
- The **all-to-wall endpoint is retained solely as the conservative bracket
  extreme** for uncertainty propagation; it is not a physical claim.
- Within-meal distribution: **wet-holdup-proportional per tray**, declared
  assumption, swept against a bottom-heavy alternative.
- Consequence if adopted: shaft heat enters `Q_other_s` (BAL-ENERGY-S) per
  tray in proportion to wet holdup. Provenance class:
  DECLARED_ENGINEERING_ASSUMPTION.

## D4 — Arm-dispersion geometry (PR-09) — **RULED 2026-08-10**

`Δz_ex` nominal 0.17 m (geometric mean of the frozen 0.09-0.33 m bracket),
swept over the full bracket; `φ_ex` nominal 0.5 (prior 0.3-0.7);
`f_state ≡ 1` declared with disclosure. First plant blade drawing replaces
`Δz_ex`. Derived nominal D_arm ≈ 1.3e-3 m²/s at 11 r/min (disclosed
arithmetic).

## D5 — Bed bulk density — **RULED 2026-08-10**

600 kg/m³ nominal wet-bed retained; bracket [480, 640] kg/m³ source-pinned
(US20210199379A1 30-40 lb/ft³); Kong 720 tapped remains the ceiling;
residence times inherit the sweep.

## D6 — Discharge configuration — **RULED 2026-08-10**

Chutes on middle trays + variable-speed rotary valve at the sparge tray
(CPM Crown, August 2024) as declared configuration; generic uncalibrated
valve law, speed as MV; PR-07 draw-footprint bracket unchanged.

## D7 — Architecture reconciliation — **RULED 2026-08-10 (Tier A)**
*(ruled after the owner's dynamic-response-fidelity question was answered:
packets are the transient state — tanks-in-series with composition memory,
exact ledgers, causal step propagation, event exactness — with stiff
gas/film physics solved implicitly per cell and C5b's switchable
dynamic-EOS gas; the retired single-DAE plan was a different
discretization, not a more dynamic one)*

Packet-carrier C5b per frozen PR-02/04/06/07; DAE/analytic-Jacobian work
scoped to per-cell stiff closures (Phase-4); tray-type object formalized
(four floor families with own free-area/pressure-loss laws, populated from
the patent corpus); stale sections of docs/CORE2_ARCHITECTURE.md marked
superseded.

## D8 — Generic campaign supervisor — **RULED 2026-08-10 (Tier B)**

One audited, manifest-parameterized supervisor replacing per-packet
specialized copies; identical check content; built as its own packet after
T3C6.

---

## Standing constraints (unchanged by this amendment)

- All 53 frozen PHY decisions and PR-01..PR-11 forms stand.
- No assumption enters a certified-particle or frozen-oracle source file.
- The six protected owner paths remain untouched.
- T3C6 manufactured-foundation evidence remains manufactured-only; the
  assumption regime feeds the SUBSEQUENT design-basis campaign tier, with
  its own packets and sweeps.
