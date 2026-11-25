"""Universal-kernel envelope hygiene for the Paper-1 radial particle kernel (2026-10-19).

Four declared knobs, each read from the environment as exactly "1" and off when unset:

* ``DTDC_KERNEL_ENVELOPE_FRONT_FILM_DEPLETION_EVENT``: in the sliding-front stage, an
  exterior-film depletion event solved by the event-first pattern of the dry stage (the
  film fixed at exactly zero, the event duration an unknown, the hexane front free),
  followed by the film-free surface of owner ruling
  GT_PS2_OWNER_RULING_PAPER1_FILM_IMBIBITION_BRIDGE_2026-10-17 point 4 for the rest of the
  history (the front stage's existing gas-only surface mode, whose surface composition is
  a gas state at the surface temperature, so its water activity is P y_w phi / f_liq(T_s)
  through the installed bridge; the exterior-liquid fast path and the film-to-shell
  transfer leave with the film, as in the dry stage).
* ``DTDC_KERNEL_ENVELOPE_FILM_NODE_CERTIFICATE_SCALED``: the external water-film node's
  condition certificate taken on the column-2-norm-scaled physical Jacobian the
  domain-aware Newton solves, against its own cap; the native-pullback value is still
  computed and reported; the shared ``controls.maximum_condition_proxy`` (1e12) is
  untouched at every other site.
* ``DTDC_KERNEL_ENVELOPE_PRESSURE_SCHEDULE``: the bath's pressure honoured (feed stage at
  the bath's first pressure, birth through the native birth assembly's pressure-only
  target configuration, front and dry stages through the native pressure-only target
  configuration per advance).
* ``DTDC_KERNEL_ENVELOPE_BIRTH_SURFACE_BOUND``: the surface-temperature upper bound of the
  birth and of the film node stated physically: the larger of the bath temperature and
  the hottest particle temperature of the step's start, plus a declared margin, instead
  of the bath temperature.

Owner ruling GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19 items 3 and 2 add
``CONDUCTION_BOUNDARY``, ``STALLED_FRONT`` and ``GAS_DEW_LINE_PROJECTION``
(particle_kernel_pd_boundary); its item 7 adds:

* ``DTDC_KERNEL_ENVELOPE_DRY_ALLOWANCE_PER_UNKNOWN``: the dry advance's residual-assembly
  allowance (the former fixed 1,400 per physical advance and per dry Newton solve) declared
  per unknown at three times the measured minimum, 3 (4 u + 6).
* ``DTDC_KERNEL_ENVELOPE_NATIVE_LEDGER_FACTOR_CEILING``: a report-only declared ceiling on
  the noise-aware energy-row gate's native-ledger factor (a typed stamp when exceeded,
  never a refusal).

The lead's record GT_PS2_UNIVERSAL_KERNEL_DRY_WALL_AND_DEWLINE_EDGE_2026-10-20 section 3
item 1 (Class A) adds:

* ``DTDC_KERNEL_ENVELOPE_FILM_FREE_SEED_STORED_SURFACE``: the dry-stage dispatcher's seed for
  the first film-free advance after the exterior film depletes (the absent-film active set
  on the advance whose previous phase was the film-bearing one) is the stored surface state
  of the last accepted row (its temperature and composition) instead of the bath gas state.
  A seed only: the path proof, its frozen temperature-derivative envelope, subdivision floor
  and point limit, every gate, tolerance, budget and law are unchanged.

The lead's numerics queue of 2026-10-20 (ruling packet 2, section 8; record
GT_PS2_UNIVERSAL_KERNEL_NUMERICS_QUEUE_BUILD_2026-10-20; Class A, numerics only) adds:

* ``DTDC_KERNEL_ENVELOPE_FILM_FREE_SEED_GAS_FALLBACK``: the companion of the film-free seed
  knob, armed only together with it (alone it is refused at launch). The seed order becomes
  stored surface first, then the bath gas state: when the absent-film active set refuses the
  stored-surface seed with the conservative absent elimination's typed refusal
  (ABSENT_FILM_CHART_REFUSAL), the same active set is tried once more on the same advance and
  the same shared allowance from the former bath gas state seed, before the contact and liquid
  sets. Any other refusal of the stored seed falls through as before.
* ``DTDC_KERNEL_ENVELOPE_DRY_HALF_STEP_FIRST``: a dry advance that starts with an exterior
  film abandons its full step once the shared per-advance assembly counter reaches the
  declared numerics threshold DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD and continues with the
  ordinary halving from the half step on the same allowance (the allowance itself unchanged).
* ``DTDC_KERNEL_ENVELOPE_DEPLETION_CARRIES_CONDENSATE``: under the dew-line projection, the
  eliminated dry-stage depletion event carries the projection's booked condensate in its
  exact-zero water closure (the film node's water row already books it), in place of the
  typed refusal "the eliminated depletion event does not carry the dew-line condensate". The
  condensate rate, the projection and its offset are unchanged.
* ``DTDC_KERNEL_ENVELOPE_INSTRUMENT_SHARED_LAW``: report-only instruments of the harness read
  the water activity on the shared equilibrium (the kernel's own chart function, PHY-053
  second-virial gas over immiscible pure liquids) in place of the ideal immiscible
  Raoult-Dalton reading (1 - y_h) P / P_sat,w(T); both readings are stamped. No decision
  depends on an instrument.

The lead's record GT_PS2_UNIVERSAL_KERNEL_PD_ADVANCE_DEMAND_2026-10-20 section 2 (ruling
packet 2, item 5; Class A as to the build, PHY-014 unchanged) adds:

* ``DTDC_KERNEL_ENVELOPE_STALL_FILM_DEPLETION_EVENT``: the companion of ``STALLED_FRONT``,
  armed only together with it (alone it is refused at launch). While the front is held in a
  stalled step, the exterior film's depletion is located as an event: the stalled
  construction (the front a literal, its slot and the water jump row out of the search, the
  water row the q = 0 consistency row at the unchanged gate, the native stall's accept chart)
  with the film unknown replaced by the event duration exactly as the front-stage event of
  ``FRONT_FILM_DEPLETION_EVENT`` does it (the film fixed at exactly zero, the drain held on its
  capacity branch, transfer share 1). After the event the front stage is film-free and the
  dispatcher tries the moving (receding) chart in its gas-only surface mode before any further
  stalled step. Every accepted stalled step and every stalled event carries a typed materiality
  row: the counterfactual advance demand of evidence pd_advance_demand_170 (the signed front
  velocity that closes the three jump rows at the held front, evaluated in a scratch overlay
  and never adopted), integrated over the step, with its cumulative value per packet and its
  ratio to the particle's hexane. Report only: the front never advances (PHY-014).

Owner ruling GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19 item 7's last part (the
surface-graded radial ladder; debt D2 of 2026-10-22, record
GT_PS2_PAPER1_SURFACE_GRADED_LADDER_2026-10-22; Class A numerics) adds:

* ``DTDC_KERNEL_ENVELOPE_SURFACE_GRADED_MESH``: the master grid of the birth is the declared
  surface-graded face set of particle_surface_graded_mesh (a uniform interior and a
  geometric band under the surface, the packet GT_PS2_PAPER1_SURFACE_LAYER_MESH_PACKET_
  2026-10-18), its three parameters given to the runner (``--graded-interior-cells``,
  ``--graded-band-cells``, ``--graded-band-ratio``) and refused without the knob, the knob
  refused without them; ``--cells`` is then the graded total. Off, the uniform grid.

The lead's Class-A seed route at the Paper-1 core extinction (record
GT_PS2_PAPER1_CORE_TARGET_SEED_ROUTE_BUILD_2026-10-22; the graded 80-cell stop of 2026-10-22, its
diagnosis and its refutation; evidence core_target_seed_177) adds one VALUED knob, armed by exactly
one of its declared values (not "1"), off when unset, any other value refused at launch:

* ``DTDC_KERNEL_ENVELOPE_CORE_TARGET_SEED=history_power_law``: the sliding-front stage's positive
  core target (the half-volume target in the core cell) seeds its duration unknown at the halving
  time of dz/dt = -c z^p fitted through the secants of the last two accepted moving front steps,
  each secant placed at the geometric mean of its end points: T = (z^(1-p) - (z/2)^(1-p)) /
  (c (1 - p)), ln 2 / c at p = 1, in place of the driver's 0.5 z / (last secant). Fewer than two
  positive accepted secants, or a fitted seed that is not finite and positive, keep the driver's
  seed. A seed only: the target, the 2,400-call search budget and its 800-call cycles, the 1e-5
  seed-reuse threshold, the 1e-11 search target, the 1e-10 gate, the 1e12 condition cap, the
  per-area chart, the four-target 0.82-contraction endpoint rule, the half-volume spacing, the
  minimum interval and every law are unchanged. Read only through a live envelope activation
  (``core_target_seed`` validates the whole environment first).

The particle's surface and path-face law (the final adoption of 2026-10-25, built on the owner's
direction "finalize it"; record docs/GT_PS2_PAPER1_SURFACE_LAW_FINAL_ADOPTION_2026-10-25.md) is
NOT a knob. Three items, each revertible by ONE line of this module (the record states each line):

* A, the film-covered surface law (experiment_core2_radial_fast_wetting.surface_contact_face).
  While an external water film covers the particle, the water taken up from the film across the
  surface face, over its three paths (the native retained path on the bridged retained activity
  of the actual states, the added liquid path and the pore-gas route), is in series with the one
  declared contact conductance K_L_MOL_M2_S, except the water vapour that leaves the pores with
  the outward gas flow: N_exempt = max(0, N_stefan), N_total = N_exempt - F / (d / L_eff +
  1 / K_l). Where the pore gas flows inward the whole gas route is in the series (built and
  measured as V4 of the 2026-10-24 study knob); where it flows outward its Stefan part is exempt
  (measured as candidate C3 of the stefan-domain study, 2026-10-25). It acts on every surface
  face under the film, in the birth, the front stage, the core-extinction projection and the dry
  stage. One-line revert: MODEL_REVERT_SURFACE_LAW = True (the former face).
* B, the path-face rules (experiment_core2_radial_fast_wetting, the block above
  liquid_connected_fast_path; measured as K1 and K2 of the cap-and-liquid-force study,
  2026-10-25): K1, no zero-force RETAINED_CAP booking on an actual-path face inside any
  liquid-connected scope (the moving wet-front interface face, the film-free surface face, the
  newborn face tangent read the bridged LUIKOV_SMOOTH force); K2, the added liquid path's
  interior-face force at the common face temperature. One-line revert: MODEL_REVERT_PATH_FORCES
  = True (the former cap rule and own-temperature force).
* C, the contact conductance K_L_MOL_M2_S = 0.062 (the immersion anchor read with the
  film-covered surface at unit activity). One-line revert: K_L_MOL_M2_S = 0.32.

The study knobs that carried the former and alternative laws are removed and refused if set
(REMOVED_ITEMS): ``FINITE_CONTACT_WHOLE_UPTAKE`` (V1 to V4), ``FINITE_CONTACT_K_L``,
``SURFACE_RETAINED_PATH``, ``SURFACE_LAW_VARIANT`` (C1 to C3), ``CAP_RULE`` (K1) and
``LIQUID_FORCE`` (K2). The former statements of A and B stay reachable, besides the one-line
reverts, only through the labelled historical study switch ``SURFACE_LAW_HISTORICAL``
(``former_unresisted_routes`` for A, ``former_cap_rule_and_own_temperature_force`` for B; NOT a
model option), which the tests the record lists use.

The surface-rules study of 2026-10-24 (record docs/GT_PS2_PAPER1_SURFACE_RULES_KNOBS_BUILD_
2026-10-24.md; a Class-B physics alternative measured for the owner, nothing adopted) keeps one
VALUED knob, read by the dry-stage dispatcher only:

* ``DTDC_KERNEL_ENVELOPE_LIQUID_MOBILITY_RULE=local_state | local_state_lagged``: the added
  (factor - 1) L_ret liquid-connected path acts on an interior face when the retained loading of
  at least one of its two cells exceeds that cell's W_cap, whether or not an exterior film exists;
  the surface contact stays film-gated. local_state reads the loadings of the trial state (an
  implicit switch); local_state_lagged reads them on the advance's accepted start state (the
  active set is fixed for the advance). Off (unset), the former rule: every interior face while
  the film exists, none after its depletion.

The feed-stage refinement study of 2026-10-24 (record docs/GT_PS2_PAPER1_FEED_STAGE_REFINEMENT_
2026-10-24.md; formulation audit finding A10; a numerical refinement study, nothing adopted) adds a
VALUED knob:

* ``DTDC_KERNEL_ENVELOPE_FEED_CELLS=1 | 2 | 4 | 8``: the DT feed stage (regime A, the high-loading
  body before the radial activation) on that many uniform radial shells
  (particle_feed_shells.MultiShellFeedMarch) in place of the one-cell thermal body of
  experiment_core2_dt_feed.FeedMarch, and the activated shell profile carried onto the birth's
  master grid by a conservative remap. "1" runs the shell code at one shell (the one-cell body to
  rounding, a code-path check). Off, the feed stage and the importer execute the former
  statements.

Every knob off reproduces the former code path. A value other than "1" is refused (for the valued
knob, a value other than its declared ones). ``physically_qualifying`` is false.

The Paper-1 front-stage memory build of 2026-10-22 (record
docs/GT_PS2_PAPER1_FRONT_STAGE_MEMORY_RELEASE_BUILD_2026-10-22.md; the lead's Class-A numerics,
memory only) adds one result-neutral knob:

* ``DTDC_KERNEL_ENVELOPE_FRONT_RELEASE_SOLVE_MEMORY``: the D9 liquid-contact limiter's face
  registry (experiment_core2_radial_fast_wetting.liquid_connected_fast_path) keeps each
  residual evaluation's entry only while the exterior face it was made for is alive, instead of
  for the whole scope. The scope is the whole front stage there (the birth's and the dry
  stage's scopes are one solve or one advance), and the former registry keeps every
  evaluation's exterior face pair and its composition-path audit, about 56 to 63 KB per
  evaluation (evidence front_memory_178). The registry is read only through a live face, so
  every lookup is answered as before; no value, gate, tolerance, budget, law or history row
  changes. Result-neutral knobs (RESULT_NEUTRAL_ITEMS) do not arm the front driver's envelope
  audit; the front record and the case record carry the knob's declared row.

Paper-1 ruling packet 6 of 2026-10-22 (GT_PS2_PAPER1_RULING_PACKET_6_PER_UNKNOWN_BUDGETS_2026-10-22;
record docs/GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md; Class B, built for the owner's
confirmation, nothing adopted; ported onto the adoption baseline of 2026-10-25 with two forms
corrected, below) adds a VALUED knob:

* ``DTDC_KERNEL_ENVELOPE_BUDGETS_PER_UNKNOWN=packet6``: five fixed numerical budgets of the Paper-1
  harness declared per unknown or per planned advance, each with its present value as the floor, so
  that a run whose every solve stayed within the present values runs unchanged (a solver consults a
  budget only to stop). With n a solve's search unknowns: (1) the front search (the water-node
  restart loop of every partial, positive target and film event, and the arrival solve): the first
  two restart cycles max(800, 4 (n + 1) + 1) residual evaluations, the total max(2,400,
  12 (n + 1)); (2) the DT birth Newton: max(2,400, 11 (n + 1)), corrected from the build's
  10 (n + 1), which its own refutation measured nine assemblies short of every nine-step birth
  (10 (n + 1) + 9); (3) the front's proof-point allowance per advance: max(3,000,000,
  15,000 (n + 1)); (4) the dry half-step-first threshold: max(1,200, floor(5 A / 12)), A the
  advance's assembly allowance in force, corrected from the build's floor(6 A / 7), which left
  A / 7 for the halving, less than its own measured half-step cost; (5) the dry-stage advance
  ceiling: max(1,500, 2 P), P the requests the declared time
  controls plan from the dry stage's first clock to the horizon, the same count dividing the
  per-advance component budgets. Unset, every site passes its former literal. Read only through a
  live envelope activation; the owner rules on arming it.

The three-knob build of 2026-10-28/29 (record THREE_KNOBS_RECORD; each knob off by default and the
former statements executed bit for bit when it is unset) adds three knobs armed by exactly "1"
and one valued companion:

* ``DTDC_KERNEL_ENVELOPE_BIRTH_SEED_BOUNDED_TO_OUTER_CELL`` (the lead's Class-A seed choice on the
  measurement of GT_PS2_PAPER1_FANER_ON_THE_GRADED_MESH_LADDER_2026-10-28 sections 18 to 23): the
  energy-based newborn-shell thickness seed of experiment_core2_dt_birth.advance_birth replaced by
  min(seed, f x outer cell width) before its fit test, f the declared
  experiment_core2_dt_birth.BIRTH_SEED_OUTER_CELL_FRACTION (0.5). A search seed only: the birth
  solve, its budget, gates, certificate and every law unchanged. Its companion
  ``DTDC_KERNEL_ENVELOPE_BIRTH_SEED_OUTER_CELL_FRACTION`` (a float in (0, 1], refused otherwise and
  refused without the knob) overrides f; the birth record stamps f, whether the bound acted and
  both seed thicknesses.
* ``DTDC_KERNEL_ENVELOPE_FANER_MEASURED_VAPOUR_TRACE`` (the owner's ruling reported with
  GT_PS2_PAPER1_FANER_MEASURED_VAPOUR_TRACE_AND_DT_BATH_CONSTANTS_2026-10-28: Faner's measured
  vapour trace is the base-case bath of the Faner march): the Faner study's constant
  --bath-temperature replaced in every stage by the measured gas-temperature trace of
  faner_measured_trace (17 points, linear in time, constant after the last), only the bulk
  temperature of the boundary and of its film boundary; the trace file declared by path and its
  sha256 stamped.
* ``DTDC_KERNEL_ENVELOPE_FRONT_STAGE_NOISE_AWARE_ENERGY_GATE`` (open decision 1 of
  GT_PS2_UNIVERSAL_KERNEL_RULING_PACKET_8_FRONT_STAGE_NUMERICS_AND_STALLED_TRANSIENT_2026-10-28; the
  owner's ruling GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17 point 2 extended to the front
  stage; Class B, built off by default, NOT adopted): the shell, jump and surface energy rows of
  the front-stage film-node search accepted at max(target, 2 x their representation floor
  sum_j |J_ij| ulp(x_j), x the physical coordinates), every water and hexane row at its unchanged
  target; its arming in any reported history waits on the owner's ruling.

The kernel's Jacobian and path-proof levers of 2026-10-31 (record JACOBIAN_LEVERS_RECORD; the
lead's Class-A numerics on the measurement JACOBIAN_LEVERS_MEASUREMENT; each armed by exactly "1",
off by default, the history bit for bit the former one with the three absent):

* ``DTDC_KERNEL_ENVELOPE_JACOBIAN_BROYDEN_UPDATES``: in the front and dry stages (the scopes the
  front driver and the runner enter), the domain-aware Newton's search Jacobian held within a
  solve, Broyden-updated after every accepted step and rebuilt fresh after a stalled contraction
  (particle_numerical_repairs.BROYDEN_REBUILD_CONTRACTION_RATIO) or a failed direction; the
  convergence test, max_nfev and the dense certificate at every returned or provisional iterate
  unchanged.
* ``DTDC_KERNEL_ENVELOPE_JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES`` (the companion of the first): in
  the dry stage only, the previous converged solve's certificate Jacobian of the same layout as the
  next solve's first search Jacobian, then Broyden; refused for the front stage.
* ``DTDC_KERNEL_ENVELOPE_PATH_PROOF_BRANCH_MEMO``: the composition-path proof's interval branch
  (native_water_path_bounds.water_path_derivative_bound) memoized on the exact bits of its six
  inputs, bounded (native_water_path_bounds.PATH_PROOF_BRANCH_MEMO_BOUND); every verdict, leaf and
  point count unchanged.

The declared n-hexane-exhausted dry regime of 2026-11-01 (record HEXANE_EXHAUSTED_REGIME_RECORD; the
owner's ruling of 2026-11-01 on the diagnosis HEXANE_EXHAUSTED_REGIME_DIAGNOSIS; armed by exactly
"1", off by default, every history bit for bit the former one with it absent):

* ``DTDC_KERNEL_ENVELOPE_HEXANE_EXHAUSTED_REGIME``: at the start of every dry advance, a film-free
  binary dry state holding less than HEXANE_EXHAUSTED_LOADING_KG_KG of n-hexane per kg dry matter
  under a presented gas at or below HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP enters the pure-water vertex
  (dtdc_simulator.core2.particle.pure_water_pore_gas) by one stamped release event booked in the
  stream integral and the ledgers; the vertex advances water and energy only with y_hexane exactly
  0; a presented gas above the clamp returns the particle to the binary stage by one stamped
  re-entry event (the pore gas dosed at the gas's own y, the sorbed n-hexane at the frozen
  isotherm's equilibrium, both taken from the stream); particle_hexane_exhausted_regime carries it.

The kernel birth-cost levers of 2026-11-02 (record BIRTH_COST_LEVERS_RECORD; the lead's Class-A
numerics on the record's own profile of a tower-born particle; each armed by exactly "1", off by
default, every history bit for bit the former one with both absent):

* ``DTDC_KERNEL_ENVELOPE_JACOBIAN_LEVERS_AT_BIRTH`` (the companion of JACOBIAN_BROYDEN_UPDATES): the
  Broyden lever's scope extended to the birth's one domain-aware Newton solve (the Paper-1 runner's
  fresh_dt_birth and the cohort tower's born job); the carried companion refused there; the feed
  stage (scipy's least_squares, which the levers do not reach) unchanged.
* ``DTDC_KERNEL_ENVELOPE_DRY_FACE_FLUX_MEMO``: the pinned dry face flux memoized on the exact bits of
  its arguments and the harness's in-memory substitutions, one frame per Newton solve, every hit
  replaying its proof accounting (particle_face_flux_memo); every value, verdict and point count
  unchanged.

The tower lane's core-extinction floor of 2026-11-02 (record TOWER_CORE_EXTINCTION_FLOOR_RECORD; a
measurement knob for the owner's decision on the lead's brief, nothing adopted; a VALUED knob armed by
exactly one of its declared values, off when unset, refused by the Paper-1 runner and by the lane
worker alone):

* ``DTDC_KERNEL_ENVELOPE_TOWER_CORE_EXTINCTION_FLOOR_M=<r_f>``: the declared front radius r_f in
  metres below which the sliding-front stage stops halving the positive core and takes the core
  endpoint (particle_core_extinction.project_core_endpoint, the native fixed-time extinction
  projection) at once: the terminal wet-volume ceiling becomes (r_f / R)^3 where that lies above
  the ruled 1e-12, and inside that band the remaining-time tolerance (0.002 / 2^L) is stamped, not
  enforced; the four linked positive targets and the 0.82 contraction stay enforced; at or below
  1e-12 the ruled endpoint rule applies unchanged. Every event stamped
  (particle_core_extinction.tower_floor_event_stamp).

The kernel Jacobian column groups of 2026-11-02 (record JACOBIAN_COLUMN_GROUPS_RECORD; the lead's
Class-A numerics; a VALUED knob, off when unset, any other value refused):

* ``DTDC_KERNEL_ENVELOPE_JACOBIAN_COLUMN_GROUPS=1``: the domain-aware Newton's finite-difference
  Jacobians at the declared sites (the birth, the front stage, the dry stage's film-bearing and
  film-free solves) taken by column groups of a declared structure, every column bit for bit the
  one-at-a-time loop's (particle_jacobian_column_groups); a row moving outside the structure falls
  back to the one-at-a-time Jacobian, the first solves of a process are bit-checked.
* ``DTDC_KERNEL_ENVELOPE_JACOBIAN_COLUMN_GROUPS=dry-rowsum``: as ``1``, and the dry stage's search
  Jacobian assembles its dense rows from the per-cell inventory terms (chain 257's Broyden class; the
  certificate exact).

The deferred dense certificate of 2026-11-03 (record DEFERRED_CERTIFICATE_RECORD; the lead's Class-A
numerics, a measurement; a VALUED knob armed only by exactly "1", off when unset, any other value
refused):

* ``DTDC_KERNEL_ENVELOPE_DEFERRED_CERTIFICATE=1``: the front driver's run (qualify_sliding_front_history)
  returns each water-node solve's iterate at gtol with its dense certificate pending and takes the
  certificate after the next solve's search, then its consumers (the film node's condition gate and
  stamps, the ledger's condition number); rows and saves wait for it; a certificate refused after the
  driver moved on replays the run with it taken at once (particle_numerical_repairs). The certificate
  is never skipped or weakened, only moved later; every value is the former one.
The compiled IAPWS-95 water residual of 2026-11-03 (record IAPWS_COMPILED_RECORD; hint C of the
external parallelization study, briefed by the lead; a measurement and a prototype, not adopted):

* ``DTDC_KERNEL_ENVELOPE_IAPWS_COMPILED=1``: the water property layer's uncached residual
  (props.water._residual_uncached, the 56-term residual Helmholtz energy with its derivatives and the
  ideal part) replaced in-process by native code generated through llvmlite from the same Python
  source, every IEEE operation in the interpreter's order, exp, log and pow the C runtime functions the
  interpreter calls (particle_iapws_compiled); any input or intermediate outside the plain-double
  domain falls back to the Python residual; the exact-bit caches and their keys unchanged; proven
  bit-identical over the logged points of chain 257's two replays; RESULT-NEUTRAL (no value, no
  history row changes); Windows only.

The steam-quench arrival regime of 2026-11-03 (record STEAM_QUENCH_ARRIVAL_RECORD; the lead's decision of
2026-11-03 under the owner's direction of the same day, a study-lane engineering knob of the tower
lane, Class B, on the diagnosis STEAM_QUENCH_ARRIVAL_DIAGNOSIS; armed by exactly "1", off by default,
every history bit for bit the former one with it absent):

* ``DTDC_KERNEL_ENVELOPE_STEAM_QUENCH_ARRIVAL``: at the start of every dry advance, a binary dry state
  whose stored surface lies below the water saturation temperature at its pressure, under a presented
  gas above that temperature, with a water activity below 1 and an n-hexane fraction below
  STEAM_QUENCH_GAS_Y_HEXANE_FILM (film-bearing) or STEAM_QUENCH_GAS_Y_HEXANE_FILM_FREE (film-free) at
  both ends of the span, takes one stamped condensation event: steam from the presented gas condenses
  until the outer cell and the exterior film stand at the saturation temperature (the outer cell's
  n-hexane held), the mass from the kernel's own energy balance and film caloric, booked as an inward
  stream in the stream integral and the ledgers; the binary liquid set then advances the event's state
  under the cell's own gas (particle_steam_quench_arrival).

The true-gas condensing-film branch of 2026-11-06 (record STEAM_CONDENSING_FILM_RECORD; the owner's
ruling STEAM_CONDENSING_FILM_RULING, the scoped true-gas steam-condensation branch, on the reviewed
candidate STEAM_CONDENSING_FILM_REVIEW; armed by exactly "1", off by default, every history bit for bit
the former one with it absent):

* ``DTDC_KERNEL_ENVELOPE_STEAM_CONDENSING_FILM``: at the start of every dry advance, a binary dry state
  under a presented gas below STEAM_CONDENSING_FILM_GAS_Y_HEXANE in n-hexane, superheated (water
  activity below 1) and hotter than the stored surface at both ends of the span, whose stored surface
  lies below the gas's water dew point on the kernel's own chart or which carries an exterior film,
  advances through the condensing-film search (particle_steam_condensing_film): a finite inward search
  box for the native Stefan inverse derived from the native coefficients and composition, Newton
  retries and a surface-energy polish charged to the original allowances, the liquid set of a film-free
  particle solved from exactly zero film, then the former dispatch; every accepted root's beta inside
  the ORIGINAL domain (|beta| <= 3) and its film profile audited, stamped; a typed refusal otherwise.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

PREFIX = "DTDC_KERNEL_ENVELOPE_"
FRONT_FILM_DEPLETION_EVENT = "FRONT_FILM_DEPLETION_EVENT"
FILM_NODE_CERTIFICATE_SCALED = "FILM_NODE_CERTIFICATE_SCALED"
PRESSURE_SCHEDULE = "PRESSURE_SCHEDULE"
BIRTH_SURFACE_BOUND = "BIRTH_SURFACE_BOUND"
# Owner ruling GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19 items 3a, 3b and 2
# (interim); particle_kernel_pd_boundary carries their behaviour.
CONDUCTION_BOUNDARY = "CONDUCTION_BOUNDARY"
STALLED_FRONT = "STALLED_FRONT"
GAS_DEW_LINE_PROJECTION = "GAS_DEW_LINE_PROJECTION"
# Owner ruling GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19 item 7 (2026-10-20
# completion build): the dry advance's residual-assembly allowance declared per unknown, and
# the report-only ceiling on the noise-aware energy-row gate's native-ledger factor.
DRY_ALLOWANCE_PER_UNKNOWN = "DRY_ALLOWANCE_PER_UNKNOWN"
NATIVE_LEDGER_FACTOR_CEILING = "NATIVE_LEDGER_FACTOR_CEILING"
# The lead's record of 2026-10-20 (the 378 K dry wall), section 3 item 1: the film-free seed.
FILM_FREE_SEED_STORED_SURFACE = "FILM_FREE_SEED_STORED_SURFACE"
# The lead's numerics queue of 2026-10-20 (ruling packet 2, section 8).
FILM_FREE_SEED_GAS_FALLBACK = "FILM_FREE_SEED_GAS_FALLBACK"
DRY_HALF_STEP_FIRST = "DRY_HALF_STEP_FIRST"
DEPLETION_CARRIES_CONDENSATE = "DEPLETION_CARRIES_CONDENSATE"
INSTRUMENT_SHARED_LAW = "INSTRUMENT_SHARED_LAW"
# The Paper-1 front-stage memory build of 2026-10-22 (result-neutral; memory only).
FRONT_RELEASE_SOLVE_MEMORY = "FRONT_RELEASE_SOLVE_MEMORY"
# The lead's record of 2026-10-20 on the pre-desolventizer advance demand, section 2.
STALL_FILM_DEPLETION_EVENT = "STALL_FILM_DEPLETION_EVENT"
# Owner ruling of 2026-10-19 item 7's graded ladder (debt D2, 2026-10-22).
SURFACE_GRADED_MESH = "SURFACE_GRADED_MESH"
# The three-knob build of 2026-10-28/29 (record THREE_KNOBS_RECORD; module docstring).
BIRTH_SEED_BOUNDED_TO_OUTER_CELL = "BIRTH_SEED_BOUNDED_TO_OUTER_CELL"
FANER_MEASURED_VAPOUR_TRACE = "FANER_MEASURED_VAPOUR_TRACE"
FRONT_STAGE_NOISE_AWARE_ENERGY_GATE = "FRONT_STAGE_NOISE_AWARE_ENERGY_GATE"
# BIRTH_SEED_BOUNDED_TO_OUTER_CELL's valued companion: a float in (0, 1], not one of ITEMS.
BIRTH_SEED_OUTER_CELL_FRACTION = "BIRTH_SEED_OUTER_CELL_FRACTION"
THREE_KNOBS_RECORD = "docs/GT_PS2_PAPER1_THREE_HARNESS_KNOBS_BUILD_2026-10-29.md"
# The kernel's Jacobian and path-proof levers of 2026-10-31 (module docstring): Class-A numerics,
# off by default; CARRIED is the companion of BROYDEN (validate_environment).
JACOBIAN_BROYDEN_UPDATES = "JACOBIAN_BROYDEN_UPDATES"
JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES = "JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES"
PATH_PROOF_BRANCH_MEMO = "PATH_PROOF_BRANCH_MEMO"
JACOBIAN_LEVERS_RECORD = "docs/GT_PS2_DTDC_KERNEL_JACOBIAN_LEVERS_BUILD_2026-10-31.md"
# The declared n-hexane-exhausted dry regime (owner ruling 2026-11-01; module docstring).
HEXANE_EXHAUSTED_REGIME = "HEXANE_EXHAUSTED_REGIME"
HEXANE_EXHAUSTED_REGIME_RECORD = (
    "docs/GT_PS2_UNIVERSAL_KERNEL_HEXANE_EXHAUSTED_REGIME_BUILD_2026-11-01.md"
)
HEXANE_EXHAUSTED_REGIME_DIAGNOSIS = (
    "docs/GT_PS2_UNIVERSAL_KERNEL_STRIPPED_PARTICLE_IN_SUPERHEATED_STEAM_DIAGNOSIS_2026-11-01.md"
)
# X_ex, kg n-hexane per kg dry matter: the diagnosis's suggested exhaustion loading (its section 7:
# five decades below the 100 to 500 ppm DT outlet band, 1.35e5 times below Paper 1's smallest
# loading 1.35e-4), ruled 2026-11-01. WAS: no regime (the binary dry stage at every loading).
HEXANE_EXHAUSTED_LOADING_KG_KG = 1e-9
# The gas condition: the kernel's bath clamp, dtdc_particle_destiny_march.interpolate_bath's
# y = min(1 - 1e-6, max(1e-6, y)) (the diagnosis's section 1, "the kernel never sees a hexane-free
# gas"); a presented gas at or below it is hexane-free for the regime. WAS: no regime.
HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP = 1e-6
# The kernel birth-cost levers of 2026-11-02 (module docstring): Class-A numerics, off by default;
# JACOBIAN_LEVERS_AT_BIRTH is the companion of JACOBIAN_BROYDEN_UPDATES (validate_environment).
JACOBIAN_LEVERS_AT_BIRTH = "JACOBIAN_LEVERS_AT_BIRTH"
DRY_FACE_FLUX_MEMO = "DRY_FACE_FLUX_MEMO"
# The compiled IAPWS-95 water residual (2026-11-03; module docstring): armed by exactly "1",
# result-neutral (bit-identical by construction and by the record's byte comparison).
# WAS (until 2026-11-03): no knob; the residual evaluated by the interpreter only.
IAPWS_COMPILED = "IAPWS_COMPILED"
IAPWS_COMPILED_RECORD = "docs/GT_PS2_DTDC_KERNEL_IAPWS_COMPILED_BIT_IDENTITY_2026-11-03.md"
# The steam-quench arrival regime (2026-11-03; module docstring): armed by exactly "1".
STEAM_QUENCH_ARRIVAL = "STEAM_QUENCH_ARRIVAL"
STEAM_QUENCH_ARRIVAL_RECORD = (
    "docs/GT_PS2_UNIVERSAL_KERNEL_STEAM_QUENCH_ARRIVAL_REGIME_BUILD_2026-11-03.md"
)
STEAM_QUENCH_ARRIVAL_DIAGNOSIS = (
    "docs/GT_PS2_DTDC_COHORT_TOWER_NEW_CELL_ARRIVAL_REFUSAL_DIAGNOSIS_2026-11-03.md"
)
# The trigger's n-hexane thresholds of the presented gas (mol/mol), read from the measured grids: a
# film-bearing particle (the diagnosis's P000397, section 3c) advances at 3e-3 and above and refuses
# at 1e-3 and below; a film-free one (the record's grid on P000004 and P000001) refuses at 0.01, 0.03
# and 0.1 (P000004) and advances at 0.2 (both), but 0.2 fires on the Paper-1 DT particle once its
# film is gone (its bath y 0.05, which it advances in), so the film-free threshold is the same 3e-3
# and the film-free gap 3e-3 to 0.2 is not covered (stated in the record).  WAS: no regime.
STEAM_QUENCH_GAS_Y_HEXANE_FILM = 3e-3
# The n-hexane-exhausted regime's re-entry hysteresis (2026-11-03, the steam-quench build's standing
# mandate; record STEAM_QUENCH_ARRIVAL_RECORD): armed by exactly "1".  Armed, a vertex state re-enters
# the binary stage only above HEXANE_EXHAUSTED_REENTRY_GAS_Y_HEXANE (the record's grid on the re-entered
# PD2 particle P000397: the binary stage refuses at 5e-5 and 1e-4 and advances at the declared value).
# WAS: re-entry above HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP.
HEXANE_EXHAUSTED_REENTRY_HYSTERESIS = "HEXANE_EXHAUSTED_REENTRY_HYSTERESIS"
HEXANE_EXHAUSTED_REENTRY_GAS_Y_HEXANE = 1e-3
# The true-gas condensing-film branch (2026-11-06; module docstring): armed by exactly "1".
STEAM_CONDENSING_FILM = "STEAM_CONDENSING_FILM"
# The condensing-film regime's two Newton levers (numerics_b, 2026-11-07; record
# CONVERGENCE_ROOT_CAUSES_RECORD; Class-A numerics, armed by exactly "1", off by default; no residual,
# gate, allowance, tolerance or beta domain moved).
CONVERGENCE_ROOT_CAUSES_RECORD = (
    "docs/GT_PS2_DTDC_COHORT_TOWER_CONVERGENCE_ROOT_CAUSES_2026-11-07.md"
)
# The seeded start (the companion of STEAM_CONDENSING_FILM and JACOBIAN_BROYDEN_UPDATES): a regime
# solve whose state unknowns (cell temperatures and compositions, the film) equal the previous
# converged regime solve's final state of the same job to 1e-12 starts from that state plus the
# previous increment scaled by the advance durations, when the residual there is admissible and
# smaller than at the former start; both evaluations charged.
STEAM_CONDENSING_FILM_SEEDED_START = "STEAM_CONDENSING_FILM_SEEDED_START"
# The floor exit (the companion of STEAM_CONDENSING_FILM): the regime's first Newton is handed to
# the regime's own surface-energy polish once it stagnates at the rounding floor (its best largest
# scaled residual within FLOOR_EXIT_FACTOR x gtol and not halved over FLOOR_EXIT_PATIENCE x the
# unknowns' evaluations); the polish, its fresh native re-certification on the full residual and
# every gate unchanged; when the polish does not close, the former retry from the best iterate.
STEAM_CONDENSING_FILM_FLOOR_EXIT = "STEAM_CONDENSING_FILM_FLOOR_EXIT"
# The sub-step request memory (the companion of STEAM_CONDENSING_FILM and JACOBIAN_BROYDEN_UPDATES):
# after a regime advance accepted below the regime's cap, the next request is at most twice the
# accepted duration (the memory lives in the job's dry lever scope); the cap itself unchanged.
STEAM_CONDENSING_FILM_REQUEST_MEMORY = "STEAM_CONDENSING_FILM_REQUEST_MEMORY"
# The film-free request cap in hexane-free steam (the companion of STEAM_CONDENSING_FILM): a
# film-free binary dry state under a superheated gas below STEAM_CONDENSING_FILM_GAS_Y_HEXANE at
# both ends, hotter than its surface, with the condensing regime's trigger not holding, requests at
# most the dry schedule's film-present band value (0.3125 s); the halving ladder unchanged.
STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP = "STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP"
# The dry budget halving (numerics_b, 2026-11-07; record CONVERGENCE_ROOT_CAUSES_RECORD; Class-A
# numerics, armed by exactly "1"): a dry sub-step attempt that exhausts its residual-assembly
# allowance enters the dry halving ladder (a refusal row, the next halving on a fresh allowance of
# the same declared size) instead of ending the job.  A declared change of the allowance's SCOPE:
# its size is unchanged but it now bounds one ATTEMPT, not one sub-step, so the assemblies spent
# on one sub-step may exceed it, by at most DRY_BUDGET_HALVING_MAX further attempts (at most
# (1 + DRY_BUDGET_HALVING_MAX) allowances per sub-step); every retry stamped with its attempt
# count and the assemblies spent so far.  Knob absent: the former re-raise.
DRY_BUDGET_HALVING = "DRY_BUDGET_HALVING"
# UNDERIVED numerical choice (numerics_b, 2026-11-07, on the lead's request for a bound): the
# budget halvings allowed per dry sub-step
DRY_BUDGET_HALVING_MAX = 3
# The film-free wet-set skip (numerics_b, 2026-11-07; record CONVERGENCE_ROOT_CAUSES_RECORD;
# Class-A numerics, armed by exactly "1"): a film-free dry advance whose surface lies more than
# FILM_FREE_WET_SET_SKIP_MARGIN_K above the presented gas's water dew point at both ends of the
# span tries the absent set only (the contact and liquid sets need a liquid supply that cannot
# form on that surface).  Knob absent: the former set order.
FILM_FREE_WET_SET_SKIP = "FILM_FREE_WET_SET_SKIP"
# UNDERIVED numerical choice (numerics_b, 2026-11-07): the surface's margin above the gas dew
# point, in kelvin, below which every set is still tried
FILM_FREE_WET_SET_SKIP_MARGIN_K = 5.0
# The absent seed above the dew point (numerics_b, 2026-11-08; record
# CONVERGENCE_ROOT_CAUSES_RECORD; Class-A numerics, armed by exactly "1"): a film-free advance
# whose absent-set seed maps its zero-flux interface outside the gas-only chart (the seed below
# the gas's water dew point) is seeded at the gas temperature plus the declared offset with the
# gas composition.  A seed only: the set order and every law unchanged.
ABSENT_SEED_ABOVE_DEW = "ABSENT_SEED_ABOVE_DEW"
# UNDERIVED numerical choice (numerics_b, 2026-11-08): the seed's offset above the gas
# temperature, in kelvin (roots at 0.001, 0.01 and 0.1 K equal to 1.3e-12 K)
ABSENT_SEED_ABOVE_DEW_OFFSET_K = 0.01
# UNDERIVED numerical choice (numerics_b, 2026-11-08): the seed is kept only when its zero-flux
# interface lies inside the chart by this relative margin (the coupling gate's width in ln a)
ABSENT_SEED_ABOVE_DEW_CHART_MARGIN = 1e-3
STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_S = 0.3125  # the dry schedule's film-present band
STEAM_CONDENSING_FILM_RECORD = "docs/GT_PS2_DTDC_COLD_START_CONDENSATION_INTEGRATION_2026-11-07.md"
STEAM_CONDENSING_FILM_RULING = (
    "docs/GT_PS2_OWNER_RULING_STEAM_CONDENSATION_AND_CONTINUATION_2026-11-04.md"
)
STEAM_CONDENSING_FILM_REVIEW = (
    "docs/GT_PS2_CONDENSATION_REVIEWED_PARTICLE_AND_VAPOUR_ASSESSMENT_2026-11-05.md"
)
# The trigger's n-hexane threshold of the presented gas (mol/mol): the diagnosis's measured film-bearing
# grid (docs/GT_PS2_DTDC_COHORT_TOWER_NEW_CELL_ARRIVAL_REFUSAL_DIAGNOSIS_2026-11-03.md, section 3c): the
# former laws advance a film-bearing arrival at 3e-3 and above and refuse at 1e-3 and below, so the
# branch takes the gases below it.  WAS: no regime.
STEAM_CONDENSING_FILM_GAS_Y_HEXANE = 3e-3
# The regime's dry-advance time cap, s: the reviewed candidate's coarse history resolution (64 of 64
# accepted 0.0625 s advances on both kept 56 s arrivals and on the 100.125 s P001896; the record's
# measurement), applied while the trigger holds over the requested advance; the kernel's own dry
# request at the 2 s tower span (0.3125 s with a film) spent the shared allowance on both kept
# arrivals.  WAS: no regime (the runner's dry schedule).
STEAM_CONDENSING_FILM_MAXIMUM_DT_S = 0.0625
STEAM_QUENCH_GAS_Y_HEXANE_FILM_FREE = 3e-3
BIRTH_COST_LEVERS_RECORD = "docs/GT_PS2_DTDC_KERNEL_BIRTH_COST_LEVERS_2026-11-02.md"
JACOBIAN_LEVERS_MEASUREMENT = "docs/GT_PS2_DTDC_KERNEL_JACOBIAN_LEVER_MEASURED_2026-10-31.md"
BIRTH_SEED_BOUND_MEASUREMENT = "docs/GT_PS2_PAPER1_FANER_ON_THE_GRADED_MESH_LADDER_2026-10-28.md"
FANER_TRACE_MEASUREMENT = (
    "docs/GT_PS2_PAPER1_FANER_MEASURED_VAPOUR_TRACE_AND_DT_BATH_CONSTANTS_2026-10-28.md"
)
FRONT_GATE_MEASUREMENT = (
    "docs/GT_PS2_UNIVERSAL_KERNEL_THIN_BIRTH_FIRST_STEP_MECHANISM_2026-10-28.md"
)
FRONT_GATE_PACKET = (
    "GT_PS2_UNIVERSAL_KERNEL_RULING_PACKET_8_FRONT_STAGE_NUMERICS_AND_STALLED_TRANSIENT_2026-10-28"
)
ITEMS = (
    FRONT_FILM_DEPLETION_EVENT,
    FILM_NODE_CERTIFICATE_SCALED,
    PRESSURE_SCHEDULE,
    BIRTH_SURFACE_BOUND,
    CONDUCTION_BOUNDARY,
    STALLED_FRONT,
    GAS_DEW_LINE_PROJECTION,
    DRY_ALLOWANCE_PER_UNKNOWN,
    NATIVE_LEDGER_FACTOR_CEILING,
    FILM_FREE_SEED_STORED_SURFACE,
    FILM_FREE_SEED_GAS_FALLBACK,
    DRY_HALF_STEP_FIRST,
    DEPLETION_CARRIES_CONDENSATE,
    INSTRUMENT_SHARED_LAW,
    FRONT_RELEASE_SOLVE_MEMORY,
    STALL_FILM_DEPLETION_EVENT,
    SURFACE_GRADED_MESH,
    BIRTH_SEED_BOUNDED_TO_OUTER_CELL,
    FANER_MEASURED_VAPOUR_TRACE,
    FRONT_STAGE_NOISE_AWARE_ENERGY_GATE,
    JACOBIAN_BROYDEN_UPDATES,
    JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES,
    PATH_PROOF_BRANCH_MEMO,
    HEXANE_EXHAUSTED_REGIME,
    JACOBIAN_LEVERS_AT_BIRTH,
    DRY_FACE_FLUX_MEMO,
    IAPWS_COMPILED,
    STEAM_QUENCH_ARRIVAL,
    HEXANE_EXHAUSTED_REENTRY_HYSTERESIS,
    STEAM_CONDENSING_FILM,
    STEAM_CONDENSING_FILM_SEEDED_START,
    STEAM_CONDENSING_FILM_FLOOR_EXIT,
    STEAM_CONDENSING_FILM_REQUEST_MEMORY,
    STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP,
    DRY_BUDGET_HALVING,
    FILM_FREE_WET_SET_SKIP,
    ABSENT_SEED_ABOVE_DEW,
)
# The lead's Class-A seed route at the Paper-1 core extinction (2026-10-22): a VALUED knob, armed
# by exactly one of its declared values; kept out of ITEMS, whose knobs are armed by exactly "1".
CORE_TARGET_SEED = "CORE_TARGET_SEED"
CORE_TARGET_SEED_HISTORY_POWER_LAW = "history_power_law"
# The surface-rules study knob of 2026-10-24 (record SURFACE_RULES_RECORD; a Class-B physics
# alternative for the owner, nothing adopted): a VALUED knob, one declared value per reading.
LIQUID_MOBILITY_RULE = "LIQUID_MOBILITY_RULE"
LIQUID_MOBILITY_LOCAL_STATE = "local_state"
LIQUID_MOBILITY_LOCAL_STATE_LAGGED = "local_state_lagged"
# Paper-1 ruling packet 6 (record BUDGETS_PER_UNKNOWN_RECORD; Class B for the owner's confirmation):
# a VALUED knob with one declared value, the five per-unknown budget forms of the section below.
BUDGETS_PER_UNKNOWN = "BUDGETS_PER_UNKNOWN"
BUDGETS_PER_UNKNOWN_PACKET6 = "packet6"
# The feed-stage refinement study knob of 2026-10-24 (record FEED_CELLS_RECORD; numerics of the
# feed stage measured for the lead, nothing adopted): a VALUED knob, the shell count.
FEED_CELLS = "FEED_CELLS"
FEED_CELLS_VALUES = ("1", "2", "4", "8")
# The labelled HISTORICAL study switch of the surface-law adoption (2026-10-25), never a model
# option, a VALUED knob off by default, kept reachable only for the tests the final record lists:
# former_unresisted_routes, the former film-covered surface face (the added liquid path alone in
# series with K_l, the native retained path and the pore-gas route unresisted; item A's former
# statements); former_cap_rule_and_own_temperature_force, the former path-face rules (the
# zero-force RETAINED_CAP booking on the actual-path faces and the added liquid path's force at
# each cell's own temperature; item B's former statements). One value at a time; the other item
# stays the model's.
SURFACE_LAW_HISTORICAL = "SURFACE_LAW_HISTORICAL"
SURFACE_LAW_HISTORICAL_FORMER = "former_unresisted_routes"
SURFACE_LAW_HISTORICAL_FORMER_PATH_FORCES = "former_cap_rule_and_own_temperature_force"
# The tower lane's front-stage time level (2026-11-02, the owner's request to know what a coarser
# front step changes; a study-lane MEASUREMENT knob, nothing adopted): a VALUED knob read only by
# the cohort tower's kernel jobs (cohort_tower_kernel_jobs.front_time_level_scope), which set the
# lane worker's front-stage step cap FRONT_MAXIMUM_DT_S to 0.1 / 2**L for the duration of a born or
# advance job. Declared levels: 1 (0.05 s, the ruled level, equal to the knob absent), 0 (0.1 s),
# -1 (0.2 s), -2 (0.4 s). Only the cap moves: the front driver's minimum interval, the core-time
# tolerance, the feed-stage and birth steps, the dry schedule, the Newton, allowances, gates,
# certificate and ledgers stay as at L = 1. The Paper-1 runner and the lane worker alone refuse it.
# WAS (until 2026-11-02): no knob; the lane worker's TIME_LEVEL 1 cap 0.1 / 2**1 = 0.05 s in every
# job ("not a knob", design memo 6.1, which stays so for the lane worker itself).
TOWER_FRONT_TIME_LEVEL = "TOWER_FRONT_TIME_LEVEL"
TOWER_FRONT_TIME_LEVEL_VALUES = ("1", "0", "-1", "-2")
TOWER_FRONT_TIME_LEVEL_RECORD = "docs/GT_PS2_DTDC_TOWER_FRONT_TIME_LEVEL_MATERIALITY_2026-11-02.md"
# The tower lane's dry-stage time level (2026-11-02, the owner's direction to find speed in the
# algorithms; a study-lane MEASUREMENT knob, nothing adopted): a VALUED knob read only by the
# cohort tower's advance jobs (cohort_tower_kernel_jobs.dry_time_level_scope), which scale every cap
# of the lane worker's dry-stage bands (the Paper-1 runner's dt-film-resolved bands at L = 1,
# particle_history_time_control) by 2**(1 - L) for the duration of one advance job, so a band reads
# base / 2**L. Declared levels: 1 (the ruled bands 0.15625 s before the first landmark, 0.3125 s to
# the second and while the external film exists, 2.5 s after it; equal to the knob absent), 0, -1
# and -2 (each band 2, 4 and 8 times wider: the film-bearing 0.3125 s to 0.625, 1.25 and 2.5 s, the
# film-free 2.5 s to 5, 10 and 20 s). Only the caps move: the two landmarks (on the birth origin),
# the film rule, the 13 halvings, the depletion advance allowance, the Newton, its allowances, every
# gate, tolerance, certificate and ledger stay. The Paper-1 runner and the lane worker alone refuse
# it.
# WAS (until 2026-11-02): no knob; the lane worker's DRY_TIME_CONTROL at its TIME_LEVEL 1 in every
# advance job (bands 0.15625 s, 0.3125 s, then 0.3125 s with the film or 2.5 s without).
TOWER_DRY_TIME_LEVEL = "TOWER_DRY_TIME_LEVEL"
# The depletion guard (2026-11-02, the lead's addition after the level -1 smoke refused in its
# oldest parcels' film-depletion window): the declared value "-1+guard" runs level -1 only in an
# advance job whose particle is far from its film's depletion: in the front stage, or in the dry
# stage with the external film present and its predicted depletion (the runner's own
# curvature_corrected_depletion_duration of the state's film points; unknown counts as far)
# more than the job's span plus TOWER_DRY_GUARD_COARSE_STEPS coarse film-bearing steps ahead.
# Otherwise (film near depletion, or the film gone) the job runs the ruled bands, which then
# hold until the second landmark, where the film-free band exceeds the tower's 2 s step at every
# level. The decision is taken per job at its start and stamped on the job.
# WAS (until 2026-11-02): the declared values 1, 0, -1, -2 only.
TOWER_DRY_TIME_LEVEL_VALUES = ("1", "0", "-1", "-2", "-1+guard")
TOWER_DRY_GUARD_SUFFIX = "+guard"
# UNDERIVED engineering margin of the guard (2026-11-02, this record): four coarse steps.
TOWER_DRY_GUARD_COARSE_STEPS = 4
TOWER_DRY_TIME_LEVEL_RECORD = "docs/GT_PS2_DTDC_TOWER_DRY_TIME_LEVEL_MATERIALITY_2026-11-02.md"
# The tower lane's core-extinction floor (2026-11-02, record TOWER_CORE_EXTINCTION_FLOOR_RECORD; a
# measurement knob for the owner's decision, nothing adopted): a VALUED knob, the declared front
# radius in metres below which the positive core is closed by the core endpoint at once.
# WAS (until 2026-11-02): no knob; the positive core halved to the ruled 1e-12 wet volume fraction
# with the remaining-time bound enforced at every endpoint (the tower lane as Paper 1).
TOWER_CORE_EXTINCTION_FLOOR_M = "TOWER_CORE_EXTINCTION_FLOOR_M"
TOWER_CORE_EXTINCTION_FLOOR_M_VALUES = ("1e-5", "3e-6", "1e-6", "1e-7", "1e-8")
TOWER_CORE_EXTINCTION_FLOOR_RECORD = (
    "docs/GT_PS2_DTDC_TOWER_CORE_EXTINCTION_FLOOR_MATERIALITY_2026-11-02.md"
)
# The tower lane's RULED time resolution (owner ruling 2026-11-02, "implement step cap 0.2 s and an
# extinction floor 1 micrometer", on docs/GT_PS2_DTDC_TOWER_FRONT_TIME_LEVEL_MATERIALITY_2026-11-02.md
# and docs/GT_PS2_DTDC_TOWER_CORE_EXTINCTION_FLOOR_MATERIALITY_2026-11-02.md; adopted in
# TOWER_RULED_TIME_RESOLUTION_RECORD): the cohort tower's "+ruled" knob sets arm the two valued knobs
# above at exactly these values, and the cohort kernel jobs refuse a "+ruled" set at any other. The
# tower lane only: the Paper-1 runner and the lane worker alone refuse both knobs; Paper 1 untouched.
# WAS (until 2026-11-02): level 1 by design memo 6.1 (the front-stage cap 0.1 / 2**1 = 0.05 s).
TOWER_RULED_FRONT_TIME_LEVEL = "-1"
# WAS (until 2026-11-02): no floor, the ruled 1e-12 wet-volume ceiling with the 1 ms remaining-time
# bound enforced at every core endpoint.
TOWER_RULED_CORE_EXTINCTION_FLOOR_M = "1e-6"
TOWER_RULED_TIME_RESOLUTION_RECORD = (
    "docs/GT_PS2_DTDC_TOWER_RULED_TIME_RESOLUTION_ADOPTED_2026-11-02.md"
)
# The kernel Jacobian column groups (2026-11-02, record JACOBIAN_COLUMN_GROUPS_RECORD; the lead's
# Class-A numerics): a VALUED knob, "1" (exact column groups at the declared sites) or "dry-rowsum"
# (as "1", the dry stage's search Jacobian's dense rows assembled from the inventory terms).
# WAS (until 2026-11-02): no knob; every finite-difference Jacobian one column at a time.
JACOBIAN_COLUMN_GROUPS = "JACOBIAN_COLUMN_GROUPS"
JACOBIAN_COLUMN_GROUPS_VALUES = ("1", "dry-rowsum")
JACOBIAN_COLUMN_GROUPS_RECORD = "docs/GT_PS2_DTDC_KERNEL_JACOBIAN_COLUMN_GROUPS_2026-11-02.md"
# The deferred dense certificate (2026-11-03, record DEFERRED_CERTIFICATE_RECORD; the lead's Class-A
# numerics, a measurement, nothing adopted): a VALUED knob with the one value "1".
# WAS (until 2026-11-03): no knob; every dense certificate taken where its iterate is returned.
DEFERRED_CERTIFICATE = "DEFERRED_CERTIFICATE"
DEFERRED_CERTIFICATE_VALUES = ("1",)
DEFERRED_CERTIFICATE_RECORD = "docs/GT_PS2_DTDC_KERNEL_DEFERRED_CERTIFICATE_PROOF_2026-11-03.md"
VALUED_ITEMS = {
    FEED_CELLS: FEED_CELLS_VALUES,
    TOWER_FRONT_TIME_LEVEL: TOWER_FRONT_TIME_LEVEL_VALUES,
    TOWER_DRY_TIME_LEVEL: TOWER_DRY_TIME_LEVEL_VALUES,
    BUDGETS_PER_UNKNOWN: (BUDGETS_PER_UNKNOWN_PACKET6,),
    CORE_TARGET_SEED: (CORE_TARGET_SEED_HISTORY_POWER_LAW,),
    LIQUID_MOBILITY_RULE: (LIQUID_MOBILITY_LOCAL_STATE, LIQUID_MOBILITY_LOCAL_STATE_LAGGED),
    SURFACE_LAW_HISTORICAL: (
        SURFACE_LAW_HISTORICAL_FORMER,
        SURFACE_LAW_HISTORICAL_FORMER_PATH_FORCES,
    ),
    TOWER_CORE_EXTINCTION_FLOOR_M: TOWER_CORE_EXTINCTION_FLOOR_M_VALUES,
    JACOBIAN_COLUMN_GROUPS: JACOBIAN_COLUMN_GROUPS_VALUES,
    DEFERRED_CERTIFICATE: DEFERRED_CERTIFICATE_VALUES,
}
SURFACE_RULES_RECORD = "docs/GT_PS2_PAPER1_SURFACE_RULES_KNOBS_BUILD_2026-10-24.md"
BUDGETS_PER_UNKNOWN_RECORD = "docs/GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md"
BUDGETS_PER_UNKNOWN_PACKET = "GT_PS2_PAPER1_RULING_PACKET_6_PER_UNKNOWN_BUDGETS_2026-10-22"
BUDGETS_PER_UNKNOWN_EVIDENCE = "paper/analysis/results_2026-10-24/budgets_per_unknown"
FEED_CELLS_RECORD = "docs/GT_PS2_PAPER1_FEED_STAGE_REFINEMENT_2026-10-24.md"

# ---- The particle's surface and path-face law (final adoption 2026-10-25; not a knob) ----------
#
# The three items of the final adoption, each revertible by ONE line (record
# SURFACE_LAW_FINAL_RECORD, section "The one-line reverts"):
#
# A. The film-covered surface law: the contact on every path with the outward Stefan water flux
#    exempt (experiment_core2_radial_fast_wetting.surface_contact_face). The one-line revert to the
#    former face (the added liquid path alone in series with K_l):
#    MODEL_REVERT_SURFACE_LAW = True
MODEL_REVERT_SURFACE_LAW = False
# B. The path-face rules K1 (no zero-force RETAINED_CAP booking on an actual-path face) and K2 (the
#    added liquid path's interior-face force at the common face temperature). The one-line revert to
#    the former cap rule and own-temperature force: MODEL_REVERT_PATH_FORCES = True
MODEL_REVERT_PATH_FORCES = False
# C. The liquid contact conductance K_l of manuscript Eq. (finitecontact), mol/(m2 s) per unit
#    ln a_w,ret: THE one point of definition. Every scope that carries the film-covered surface (the
#    birth, the front stage and the core-extinction projection through
#    experiment_core2_finite_wetting_b.finite_shell_path, the dry stage through
#    particle_postcore_solver.advance_dry_phases, the tower lane worker's finite-shell scope) reads
#    it here at call time. Provenance: the supplement's contact-scale anchor (sec. contactscale;
#    Hemmingsen et al. 2008, the 1.43-mm meal's 90-s hot-water uptake: a lumped sphere filled from
#    W0 = 0.134 to W_cap - 0.01 in 90 s, rho 1159.65 kg/m3, R 0.715 mm, isotherm A1 0.88, A2 12.184,
#    W_cap 0.235670904) read with the film-covered surface at UNIT activity, psi_s = 0, as the main
#    text and the code put it: the anchor integral 0.28889 (the printed reading, the surface at
#    W_cap, gave 1.48428), so K_l = 0.0492 at 80 C and 0.0619 at 100 C by inverse water viscosity
#    (factor 1.257), declared 0.062 (formulation audit a04_anchor.py; record
#    GT_PS2_PAPER1_FINITE_CONTACT_ALL_PATHS_BUILD_2026-10-24.md section 2; the final record).
#    WAS 0.32 (the printed anchor with the surface at W_cap; sec_formulation.tex parameter table
#    "Liquid contact K_l", Conditional), until the final adoption of 2026-10-25. One-line revert:
#    K_L_MOL_M2_S = 0.32
# sec_formulation.tex parameter table "Liquid contact K_l" (to be revised); supplement sec. contactscale:
K_L_MOL_M2_S = 0.062
# D. The total liquid-mobility factor F of manuscript Eq. (finitecontact) (the supplement's
#    "F = 100": the added liquid path L_l = (F - 1) L_wr, L_wr the native retained-water
#    mobility, L_l 8.39e-5 mol/(m s)): THE one point of definition. Read at call time by the four
#    sites that carried the literal 100 until 2026-10-31 (run_particle_qualification_case's birth
#    and core-endpoint scopes, qualify_sliding_front_history's DT front scope, the tower lane
#    worker's endpoint scope, particle_postcore_solver.advance_dry_phases's default), and mirrored
#    by the tower's port (src/dtdc_simulator/core2/study_lane/paper1_particle_port.py
#    LIQUID_MOBILITY_TOTAL_FACTOR; equality tested). The value is the paper's, unchanged: the
#    four sites read the same number they carried (record
#    GT_PS2_UNIVERSAL_KERNEL_PORT_REALIGNED_TO_FINAL_PAPER1_MODEL_2026-10-31.md, item 1).
LIQUID_MOBILITY_TOTAL_FACTOR = 100.0
# The law's name (V4's name extended by the exemption) and its records: the V4 derivation build
# of 2026-10-24, the adoption build and the final adoption of 2026-10-25, the measured candidates.
SURFACE_LAW = "whole_uptake_all_paths_outward_stefan_exempt"
SURFACE_LAW_RECORD = "docs/GT_PS2_PAPER1_FINITE_CONTACT_ALL_PATHS_BUILD_2026-10-24.md"
SURFACE_LAW_ADOPTION_RECORD = "docs/GT_PS2_PAPER1_SURFACE_LAW_ADOPTION_2026-10-25.md"
SURFACE_LAW_FINAL_RECORD = "docs/GT_PS2_PAPER1_SURFACE_LAW_FINAL_ADOPTION_2026-10-25.md"
SURFACE_LAW_STEFAN_DOMAIN_RECORD = "docs/GT_PS2_PAPER1_SURFACE_LAW_STEFAN_DOMAIN_2026-10-25.md"
SURFACE_LAW_STEFAN_C3_RECORD = "docs/GT_PS2_PAPER1_SURFACE_LAW_STEFAN_DOMAIN_C3_2026-10-25.md"
CAP_FORCE_RECORD = "docs/GT_PS2_PAPER1_CAP_AND_LIQUID_FORCE_KNOBS_2026-10-25.md"
# The path-face rules' name (K1 and K2 together, item B).
PATH_FORCES = "bridged_cap_rule_common_temperature_liquid_force"
SURFACE_LAW_RULING = (
    "owner, 2026-10-25: 'that is fine, we need to finalize it' (make the corrected law the model: "
    "flip the default, remove the switch, re-pin and record), then the direction 'finalize it' for "
    "the final law (the outward Stefan exemption, K1, K2, K_l 0.062), built on the lead's "
    "recommendation; the owner's one-word go was pending at the build's hand-back"
)
# The study knobs the adoptions removed: refused at launch if set, with the reason.
REMOVED_ITEMS = {
    "FINITE_CONTACT_WHOLE_UPTAKE": "removed by the surface-law adoption of 2026-10-25: its value "
    "whole_uptake_all_paths (V4) is the model's law where the pore gas flows inward; V1 to V3 are "
    "not reachable",
    "FINITE_CONTACT_K_L": "removed by the surface-law adoption of 2026-10-25: the contact "
    "conductance is the one declared constant particle_kernel_envelope.K_L_MOL_M2_S",
    "SURFACE_RETAINED_PATH": "removed by the surface-law adoption of 2026-10-25: its film-covered "
    "reading (the native retained path on the bridged force) is the adopted law's native path",
    "SURFACE_LAW_VARIANT": "removed by the final surface-law adoption of 2026-10-25: its value "
    "outward_carried_exempt (C3) is the model's law; C1 and C2 are not reachable",
    "CAP_RULE": "removed by the final surface-law adoption of 2026-10-25: its value "
    "bridged_everywhere (K1) is the model's cap rule",
    "LIQUID_FORCE": "removed by the final surface-law adoption of 2026-10-25: its value "
    "common_temperature (K2) is the model's liquid force",
}
RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_ENVELOPE_HYGIENE_BUILD_2026-10-19.md"
PD_RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_PD_BOUNDARY_BUILD_2026-10-19.md"
PD_RULING = "GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19"
COMPLETION_RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_ENVELOPE_COMPLETION_BUILD_2026-10-20.md"
SEED_RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_FILM_FREE_SEED_BUILD_2026-10-20.md"
SEED_DECISION = "GT_PS2_UNIVERSAL_KERNEL_DRY_WALL_AND_DEWLINE_EDGE_2026-10-20 section 3 item 1"
QUEUE_RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_NUMERICS_QUEUE_BUILD_2026-10-20.md"
QUEUE_DECISION = (
    "GT_PS2_UNIVERSAL_KERNEL_RULING_PACKET_2_2026-10-20 section 8 (the lead's Class-A "
    "numerics queue; item 7's depletion-event bookkeeping)"
)
RELEASE_RECORD = "docs/GT_PS2_PAPER1_FRONT_STAGE_MEMORY_RELEASE_BUILD_2026-10-22.md"
# Knobs that change no computed value, gate, budget, law or history row (memory only). Armed
# alone they do not arm the front driver's envelope audit (behaviour_armed); they are still
# validated, declared and stamped like every other knob.
# IAPWS_COMPILED added 2026-11-03 (WAS: FRONT_RELEASE_SOLVE_MEMORY alone).
RESULT_NEUTRAL_ITEMS = (FRONT_RELEASE_SOLVE_MEMORY, IAPWS_COMPILED)
STALL_RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_STALL_DEPLETION_BUILD_2026-10-20.md"
STALL_DECISION = (
    "GT_PS2_UNIVERSAL_KERNEL_PD_ADVANCE_DEMAND_2026-10-20 section 2 (ruling packet 2, item 5: "
    "PHY-014 unchanged; the stalled-front mode completed as numerics, Class A as to the build)"
)
GRADED_RECORD = "docs/GT_PS2_PAPER1_SURFACE_GRADED_LADDER_2026-10-22.md"
CORE_TARGET_SEED_RECORD = "docs/GT_PS2_PAPER1_CORE_TARGET_SEED_ROUTE_BUILD_2026-10-22.md"
CORE_TARGET_SEED_DECISION = (
    "the lead's Class-A seed route at the Paper-1 core extinction of 2026-10-22 (the graded 80-cell "
    "stop diagnosed and refuted: a duration seed from accepted history passes the stored target with "
    "every budget unchanged)"
)
CORE_TARGET_SEED_EVIDENCE = (
    "paper/analysis/results_2026-10-16/solver_qualification/core_target_seed_177"
)
# CORE_TARGET_SEED (armed only): the checkpoint key of the last accepted moving front points, and
# the history row kinds that are not such points (no state accepted, or the front held).
CORE_TARGET_SEED_POINTS_KEY = "kernel_envelope_core_target_seed_points"
HELD_FRONT_ROW_KINDS = ("refusal", "stalled_front", "stalled_film_depletion_event")

# The typed refusal of the conservative absent-film elimination
# (particle_conservative_surface.absent_film_surface_for_component_targets) that the gas-state
# fallback answers; any other refusal of the stored-surface seed is not answered.
ABSENT_FILM_CHART_REFUSAL = (
    "absent-film component targets need a surface outside the gas-only chart"
)

# The half-step-first threshold, a declared NUMERICS threshold (not a physics value): the
# residual assemblies of the shared per-advance counter at which a dry advance that starts
# with an exterior film abandons its full step for the half step. Chosen from the census of
# the accepted film-bearing dry advances of the complete steam-side cases (evidence
# film_free_seed_168, 27 complete cases, 147 film-bearing advances: median 132, maximum
# 1,102 at the depletion advance of T378/0.02/101.3; evidence 166's 914 is that census's
# maximum before the 378 K points completed) and from the two film-bearing stops of evidence
# 166 (the full step failing after 2,100 to 2,800 assemblies, the half step closing in 106):
# above every accepted film-bearing advance by about 9 percent, and leaving 200 of the
# 1,400 for the halving (1.9 times the measured half-step cost). Evidence
# numerics_queue_169 records/census_film_bearing_168R.json.
DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD = 1200

# Item 7, the former allowance: 1,400 residual assemblies per dry physical advance (the
# runner's shared per-advance counter) and per dry Newton solve (its max_nfev).
FORMER_DRY_ASSEMBLY_ALLOWANCE = 1400
# Item 7, the measured minimum (evidence 151/152 and the surface-layer mesh packet of
# 2026-10-18): an accepted dry advance needs at least 4 u + 6 residual assemblies, u its
# search unknowns; under the conservative elimination the ladders run, u = 3 N + 1 (cell
# temperature, composition fraction and Stefan flux per cell, plus one surface or event
# unknown): 106, 202, 394, 778 at 8 to 64 cells.
DRY_MINIMUM_CALLS_PER_UNKNOWN = 4
DRY_MINIMUM_CALLS_OFFSET = 6
# Item 7, as ruled: the allowance is three times the measured minimum.
DRY_ALLOWANCE_MULTIPLE_OF_MINIMUM = 3
# Item 7, the native-ledger factor of the noise-aware energy-row gate (gate_i / tolerance on
# the cell energy rows, >= 1) is accepted as stamped; this ceiling is for reporting only (a
# typed stamp when exceeded, never a refusal). Set at the maximum observed over the 5,269
# accepted dry advances of the 22 complete cases of envelope ladder 164 (T388_y050_P1013,
# 8 cells, middle level); evidence envelope_completion_166.
NATIVE_LEDGER_FACTOR_REPORTING_CEILING = 6.093855124424927

# Item 2: the film node's own cap on the column-2-norm-scaled physical Jacobian. Proposed
# from the measurement (evidence 161): the accepted front-stage film-node solves of the
# stored S1 and qualification states and the scaled form of the S1 refusals (2.4e5 to
# 9.3e5 in evidence 160); the same value as the ruled birth cap
# (MAXIMUM_BIRTH_CONDITION_PROXY), on the same matrix family. Open to the owner.
FILM_NODE_SCALED_CONDITION_CAP = 1.0e8

# Item 4: the declared margin above max(bath, hottest particle temperature at the step's
# start). A surface without a heat source cannot exceed the hotter of its two thermal
# reservoirs (the gas and the particle it bounds); the margin covers the sorption and
# condensation heat released at the surface, which can lift it above both by a small
# amount. Measured need (evidence 161): the accepted counterfactual birth of evidence 160
# sits 1.08 K below the activation temperature, so the bound binds at no accepted root
# with margin 0; the declared margin is one kelvin. Open to the owner.
SURFACE_TEMPERATURE_MARGIN_K = 1.0

DEFAULT_PRESSURE_PA = 101325.0


class EnvelopeKnobRefused(ValueError):
    """A DTDC_KERNEL_ENVELOPE_* variable carries a value other than exactly "1"."""


def _step_budget_exceeded():
    from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded

    return StepBudgetExceeded


class DryFullStepAbandoned(_step_budget_exceeded()):
    """DRY_HALF_STEP_FIRST: the full step of a film-bearing dry advance reached the declared
    threshold of the shared per-advance assembly counter and is abandoned for the half step.

    A StepBudgetExceeded so that every assembler and phase handler propagates it untouched (no
    active set absorbs it as a rejected trial); only the runner's halving loop catches it, and
    only at the full step. Not a refusal of the advance: the ordinary halving continues on the
    same shared allowance.
    """

    def __init__(self, message, *, used, threshold):
        super().__init__(message)
        self.used = used
        self.threshold = threshold


HALF_STEP_ABANDON_KEY = "kernel_envelope_abandon_full_step_at"


def armed(item, environ=None):
    """True when the knob is exactly "1"; False when unset; refused otherwise."""
    if item not in ITEMS:
        raise KeyError(f"unknown kernel-envelope knob {item!r}")
    environ = os.environ if environ is None else environ
    value = environ.get(PREFIX + item)
    if value is None:
        return False
    if value != "1":
        raise EnvelopeKnobRefused(
            f"{PREFIX + item}={value!r}: a kernel-envelope knob is armed by exactly '1' "
            "and off when unset; any other value is refused"
        )
    return True


def valued(item, environ=None):
    """A valued knob's declared value when set; None when unset; refused otherwise."""
    if item not in VALUED_ITEMS:
        raise KeyError(f"unknown valued kernel-envelope knob {item!r}")
    environ = os.environ if environ is None else environ
    value = environ.get(PREFIX + item)
    if value is None:
        return None
    if value not in VALUED_ITEMS[item]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + item}={value!r}: a valued kernel-envelope knob is armed by exactly one of "
            f"{', '.join(repr(v) for v in VALUED_ITEMS[item])} and off when unset; any other value "
            "is refused"
        )
    return value


def birth_seed_outer_cell_fraction_override(environ=None):
    """BIRTH_SEED_OUTER_CELL_FRACTION: the companion's float when set (None when unset).

    Refused unless it parses as a finite float f with 0 < f <= 1 written as a plain decimal
    (no sign, exponent, whitespace, "nan" or "inf"). Read by experiment_core2_dt_birth only while
    BIRTH_SEED_BOUNDED_TO_OUTER_CELL is armed; validate_environment refuses it without the knob.
    """
    environ = os.environ if environ is None else environ
    value = environ.get(PREFIX + BIRTH_SEED_OUTER_CELL_FRACTION)
    if value is None:
        return None
    plain = value != "" and all(c in "0123456789." for c in value) and value.count(".") <= 1
    fraction = float(value) if plain and value != "." else math.nan
    if not (math.isfinite(fraction) and 0.0 < fraction <= 1.0):
        raise EnvelopeKnobRefused(
            f"{PREFIX + BIRTH_SEED_OUTER_CELL_FRACTION}={value!r}: the bounded birth seed's "
            "fraction of the outer cell is a plain decimal f with 0 < f <= 1; any other value is "
            "refused"
        )
    return fraction


def validate_environment(environ=None):
    """Refuse at launch every malformed knob and every unknown DTDC_KERNEL_ENVELOPE_* name.

    Returns the states of the knobs armed by "1" (ITEMS); a valued knob (VALUED_ITEMS) is refused
    here when malformed and is read through ``valued`` or ``core_target_seed``.
    """
    environ = os.environ if environ is None else environ
    for name in environ:
        if name.startswith(PREFIX) and name[len(PREFIX) :] in REMOVED_ITEMS:
            raise EnvelopeKnobRefused(f"{name!r}: {REMOVED_ITEMS[name[len(PREFIX) :]]}")
        if name.startswith(PREFIX) and name[len(PREFIX) :] not in (
            *ITEMS,
            *VALUED_ITEMS,
            BIRTH_SEED_OUTER_CELL_FRACTION,
        ):
            raise EnvelopeKnobRefused(f"unknown kernel-envelope knob {name!r}")
    state = {item: armed(item, environ) for item in ITEMS}
    for item in VALUED_ITEMS:
        valued(item, environ)
    if PREFIX + BIRTH_SEED_OUTER_CELL_FRACTION in environ:
        if not state[BIRTH_SEED_BOUNDED_TO_OUTER_CELL]:
            raise EnvelopeKnobRefused(
                f"{PREFIX + BIRTH_SEED_OUTER_CELL_FRACTION} is the companion of "
                f"{PREFIX + BIRTH_SEED_BOUNDED_TO_OUTER_CELL} (the bounded birth seed's fraction "
                "of the outer cell) and is refused without it"
            )
        birth_seed_outer_cell_fraction_override(environ)
    if state[FILM_FREE_SEED_GAS_FALLBACK] and not state[FILM_FREE_SEED_STORED_SURFACE]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + FILM_FREE_SEED_GAS_FALLBACK} is the companion of "
            f"{PREFIX + FILM_FREE_SEED_STORED_SURFACE} (the seed order stored surface, then the "
            "gas state) and is refused without it"
        )
    if state[JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES] and not state[JACOBIAN_BROYDEN_UPDATES]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES} is the companion of "
            f"{PREFIX + JACOBIAN_BROYDEN_UPDATES} (the carried Jacobian starts a Broyden-updated "
            "sequence, as measured) and is refused without it"
        )
    if state[STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP] and not state[STEAM_CONDENSING_FILM]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP} is the companion of "
            f"{PREFIX + STEAM_CONDENSING_FILM} and is refused without it"
        )
    if state[STEAM_CONDENSING_FILM_REQUEST_MEMORY] and not (
        state[STEAM_CONDENSING_FILM] and state[JACOBIAN_BROYDEN_UPDATES]
    ):
        raise EnvelopeKnobRefused(
            f"{PREFIX + STEAM_CONDENSING_FILM_REQUEST_MEMORY} is the companion of "
            f"{PREFIX + STEAM_CONDENSING_FILM} and {PREFIX + JACOBIAN_BROYDEN_UPDATES} (its "
            "memory lives in the job's dry lever scope) and is refused without them"
        )
    if state[STEAM_CONDENSING_FILM_FLOOR_EXIT] and not state[STEAM_CONDENSING_FILM]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + STEAM_CONDENSING_FILM_FLOOR_EXIT} is the companion of "
            f"{PREFIX + STEAM_CONDENSING_FILM} and is refused without it"
        )
    if state[STEAM_CONDENSING_FILM_SEEDED_START] and not (
        state[STEAM_CONDENSING_FILM] and state[JACOBIAN_BROYDEN_UPDATES]
    ):
        raise EnvelopeKnobRefused(
            f"{PREFIX + STEAM_CONDENSING_FILM_SEEDED_START} is the companion of "
            f"{PREFIX + STEAM_CONDENSING_FILM} and {PREFIX + JACOBIAN_BROYDEN_UPDATES} (its "
            "memory lives in the job's dry lever scope) and is refused without them"
        )
    if state[JACOBIAN_LEVERS_AT_BIRTH] and not state[JACOBIAN_BROYDEN_UPDATES]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + JACOBIAN_LEVERS_AT_BIRTH} is the companion of "
            f"{PREFIX + JACOBIAN_BROYDEN_UPDATES} (it extends that lever's scope to the birth) and "
            "is refused without it"
        )
    if state[STALL_FILM_DEPLETION_EVENT] and not state[STALLED_FRONT]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + STALL_FILM_DEPLETION_EVENT} is the companion of {PREFIX + STALLED_FRONT} "
            "(the film-depletion event of a held front) and is refused without it"
        )
    return state


def armed_items(environ=None):
    """The armed knobs: those of ITEMS armed by "1", then every valued knob that is set."""
    state = validate_environment(environ)
    return tuple(item for item, on in state.items() if on) + tuple(
        item for item in VALUED_ITEMS if valued(item, environ) is not None
    )


def behaviour_armed(knobs):
    """True when an armed knob (validate_environment's mapping) can change a computed value or
    a record row, i.e. one outside RESULT_NEUTRAL_ITEMS; the front driver's envelope audit is
    armed by these only."""
    return any(on for item, on in knobs.items() if item not in RESULT_NEUTRAL_ITEMS)


def release_solve_memory_row():
    """FRONT_RELEASE_SOLVE_MEMORY's typed declared row (the front record's stamp)."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "the lead's Class-A numerics (memory only) on the mesh-cost probe of "
        "2026-10-22; no value, gate, tolerance, budget, law or record row moved",
        "record": RELEASE_RECORD,
        "physically_qualifying": False,
        "name": "d9_face_registry_entries_released_with_their_faces",
        "what": "experiment_core2_radial_fast_wetting.liquid_connected_fast_path's registry "
        "(the D9 limiter's map from an exterior face, by id, to its native face, the added "
        "liquid flux, its enthalpy and dissipation) keeps an entry only while the face it was "
        "made for is alive (weakref.finalize on the face removes it), and does not hold the "
        "face itself",
        "former": "every entry and both faces kept for the whole scope: in the front stage the "
        "whole stage, one entry per residual evaluation (about 56 to 63 KB each with the "
        "face's composition-path audit)",
        "lookups": "the registry is read only by the limiter's supply_limited_surface through "
        "a live face object, so every lookup returns the entry the former registry returns; an "
        "entry is removed before its id can be reused",
        "result_neutral": True,
        "unchanged": "every residual, Jacobian, certificate, gate, tolerance, budget, law, "
        "history row, checkpoint and restart; the knob does not arm the front driver's "
        "envelope audit",
    }


def declared_rows(environ=None):
    """Typed DECLARED_ENGINEERING_ASSUMPTION rows of the armed knobs (empty when none)."""
    common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "universal-kernel envelope hygiene (lead, Class A under the hygiene "
        "principle of GT_PS2_OWNER_RULING_PAPER1_NUMERICAL_GATE_HYGIENE_2026-10-17), "
        "items marked owner_item await the owner",
        "record": RECORD,
        "physically_qualifying": False,
    }
    rows = {
        FRONT_FILM_DEPLETION_EVENT: {
            "name": "front_stage_film_depletion_event_then_film_free_surface",
            "event": "the front-stage water-node advance with the exterior film fixed at exactly "
            "zero and the event duration an unknown in place of the film, the hexane front free; "
            "attempted first when the film's linear extrapolation from the last two accepted "
            "rows reaches zero inside the attempted interval; a refused event falls back to the "
            "former order",
            "after_event": "the front stage's gas-only surface mode (surface T, flux coordinate, "
            "gas-only composition fraction); the film-free condition a_w = P y_w phi / "
            "f_liq(T_s) of GT_PS2_OWNER_RULING_PAPER1_FILM_IMBIBITION_BRIDGE_2026-10-17 point 4 "
            "through the installed bridge equilibrium; exterior-liquid fast path and film-to-"
            "shell transfer absent (as the dry stage's exterior_liquid_present=wet); the dry "
            "stage starts in its absent-film active set",
            "unchanged": "every gate, tolerance, budget, law and the 1e-10 node rows",
        },
        FILM_NODE_CERTIFICATE_SCALED: {
            "name": "film_node_condition_certificate_on_searched_matrix",
            "gated": "np.linalg.cond of the column-2-norm-scaled physical film-node Jacobian "
            "(particle_numerical_repairs.searched_matrix_condition)",
            "cap": FILM_NODE_SCALED_CONDITION_CAP,
            "reported_not_gated": "the native-pullback condition (former gate, cap 1e12)",
            "unchanged_maximum_condition_proxy_elsewhere": 1e12,
            "owner_item": "the extension of GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17 "
            "point 3 to the film node and the cap value",
        },
        PRESSURE_SCHEDULE: {
            "name": "bath_pressure_schedule",
            "interpolation": "linear in time between waypoints, as T and y; a bath without "
            f"pressure fields is at {DEFAULT_PRESSURE_PA} Pa",
            "feed": "the bath's first pressure in place of the hard-coded 101,325 Pa",
            "birth": "cut_birth_event.assemble_birth_backward_euler(target_config=...) at the "
            "bath pressure of the birth's end",
            "front_and_dry": "the native pressure-only target configuration per advance at the "
            "bath pressure of the advance's end (event durations: re-evaluated per trial)",
            "not_threaded": "the exact-face arrival, the core-extinction projection and the "
            "endpoint (no native target argument): those advances keep their start pressure; "
            "the declared Coletto film pair state (frozen source state, not a function of "
            "the bath's T or P)",
        },
        BIRTH_SURFACE_BOUND: {
            "name": "surface_temperature_upper_bound_stated_physically",
            "bound": "max(T_bath, hottest particle temperature at the step's start) + margin",
            "declared_margin_k": SURFACE_TEMPERATURE_MARGIN_K,
            "former": "the bath temperature (the chart 'the gas heats the particle')",
            "sites": "birth (two bound sites) and the film node (two bound sites)",
            "secant_heat_capacity": "measured and reported only (the harness carries no "
            "analytic partial-enthalpy derivative)",
            "owner_item": "the bound's statement and the margin value",
        },
        # its own authority and record override the common ones
        FRONT_RELEASE_SOLVE_MEMORY: release_solve_memory_row(),
    }
    pd_common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": f"{PD_RULING} (items 3 and 2 interim; Class-A numerics and plumbing)",
        "record": PD_RECORD,
        "physically_qualifying": False,
    }
    pd_rows = {
        CONDUCTION_BOUNDARY: {
            "name": "pre_desolventizer_jacket_conduction_boundary",
            "bath_fields": "jacket_mode ('conduction' | 'gas_contact'), "
            "jacket_wall_temperature_K, jacket_contact_conductance_w_m2_k on every waypoint "
            "or none; per interval (the right waypoint's value on (t_{k-1}, t_k])",
            "law": "in a conduction interval the surface's sensible gas term h_film A (T_s - "
            "T_b) is replaced by h_c (T_s - T_wall); the vapour leaving or condensing keeps "
            "the kernel's own film mass law against the bath gas",
            "tower_mechanism": "the PD shared-wall node: UA_floor (8000 W/K) to the bottom "
            "K layer shared by dry mass, (T_wall - T_packet), the gas contact fraction 0 "
            "(pd_k_dome_bypass_tray_integration.py:6846-6857)",
            "stages": "front (film-bearing and film-free), arrival, dry (every surface active "
            "set through the film node), and the birth through the same film node",
            "bath_without_field": "ignored; the former boundary",
        },
        STALLED_FRONT: {
            "name": "stalled_front_native_interior_stall_at_the_film_node",
            "entry": "a front-stage partial advance refused (no receding root at this "
            "duration); the moving chart is tried first on every step",
            "construction": "src/dtdc_simulator/core2/particle/cut_interior_stall.py items "
            "1-5 (ratified 2026-08-21): the front a literal (zero sweep), its slot removed, "
            "the Rankine-Hugoniot water row removed from the searched rows and judged with "
            "every other row at the unchanged gate (q = 0 water-flux continuity); the "
            "hexane row (zero net hexane release) and the energy row solved",
            "unchanged": "every gate, tolerance, budget, ledger and law",
            "stamped": "front.history.jsonl kernel_envelope.stalled_front per stalled row; "
            "stall refusals with the consistency residual",
        },
        GAS_DEW_LINE_PROJECTION: {
            "name": "bath_gas_projected_onto_the_kernel_dew_line",
            "projection": "y_h below the gas-only interval's water-saturation end y_dew(T_b, "
            "P) (a_w > 1 on the kernel's equilibrium, binary gas k_wh of the film boundary) "
            "-> y_dew",
            "presented_point": "the first representable gas state strictly inside the open "
            "gas-only chart: the smallest ulp offset of y_h from the bisected dew line at which "
            "the kernel's water activity rounds below 1 (an ulp ladder at the bath, at most "
            "4,096 ulps; stamped with the offset and the activity); a representational choice, "
            "no physical offset declared (GT_PS2_UNIVERSAL_KERNEL_DRY_WALL_AND_DEWLINE_EDGE_"
            "2026-10-20 section 3 item 3; record " + SEED_RECORD + ")",
            "condensate": "removed water 1 - y_h / y_dew per mole of gas x P / (R T_b) x M_w "
            "x the film's scaled mass-transfer coefficient, booked onto the external film "
            "through the film node's water row with the vapour partial enthalpy at "
            "(T_b, P, y_dew) in its energy row",
            "film_free_surface": "refused typed (DewLineCondensateOnFilmFreeSurface)",
            "inside_interval": "nothing changes",
            "owner_item": "interim rule of item 2; the law is ruled after the two-equilibria "
            "measurement",
        },
    }
    completion_common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": f"{PD_RULING} item 7 (as ruled; the values below are the ruling's or "
        "measured, none chosen by this build)",
        "record": COMPLETION_RECORD,
        "physically_qualifying": False,
    }
    completion_rows = {
        DRY_ALLOWANCE_PER_UNKNOWN: {
            "name": "dry_advance_residual_assembly_allowance_per_unknown",
            "allowance": f"{DRY_ALLOWANCE_MULTIPLE_OF_MINIMUM} x ({DRY_MINIMUM_CALLS_PER_UNKNOWN} "
            f"u + {DRY_MINIMUM_CALLS_OFFSET}) residual assemblies",
            "unknowns_per_physical_advance": "u = 3 N + 1, N the dry grid's cells (the "
            "conservative elimination's search unknowns, the measured minimum's basis)",
            "unknowns_per_newton_solve": "u = the size of that solve's search vector",
            "former": f"{FORMER_DRY_ASSEMBLY_ALLOWANCE} at both sites",
            "measured_minimum": "4 u + 6 (evidence 151/152; "
            "GT_PS2_PAPER1_SURFACE_LAYER_MESH_PACKET_2026-10-18 section 3 item 1)",
            "unchanged": "every gate, tolerance, law, the 13 halvings and every other budget",
        },
        NATIVE_LEDGER_FACTOR_CEILING: {
            "name": "noise_aware_energy_gate_native_ledger_factor_reporting_ceiling",
            "factor": "gate_i / tolerance on the cell energy rows passed to the native ledger "
            "(experiment_core2_dt_dry native_ledger_energy_scale_factors), accepted as stamped",
            "ceiling": NATIVE_LEDGER_FACTOR_REPORTING_CEILING,
            "ceiling_origin": "the maximum observed over the accepted dry advances of the 22 "
            "complete cases of envelope ladder 164",
            "effect": "report only: a typed stamp on the advance and a case-level count when "
            "the factor exceeds the ceiling; never a refusal",
        },
    }
    seed_rows = {
        FILM_FREE_SEED_STORED_SURFACE: {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": f"{SEED_DECISION} (Class A numerics: a search seed, no value moved)",
            "record": SEED_RECORD,
            "physically_qualifying": False,
            "name": "dry_dispatcher_film_free_seed_at_the_stored_surface",
            "seed": "the absent-film active set of the first dry advance after the exterior film "
            "depletes (previous phase the film-bearing one, the film now exactly zero) is seeded "
            "at the stored surface of the last accepted row (its temperature and composition), "
            "not at the bath gas state",
            "former": "the bath gas state (T_b, y_b) of the attempted interval's end",
            "unchanged": "the path proof (its frozen temperature-derivative envelope, subdivision "
            "floor and point limit), the active-set order, every gate, tolerance, budget and law; "
            "the seed after a contact advance stays the gas state",
            "stamped": "native_report.kernel_envelope_film_free_seed on the accepted advance; "
            "the absent trial's rejected entry when refused",
        },
    }
    queue_common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": f"{QUEUE_DECISION} (Class A numerics: no law, gate, tolerance, budget, "
        "seed family, bound, band or declared physical value moved)",
        "record": QUEUE_RECORD,
        "physically_qualifying": False,
    }
    queue_rows = {
        FILM_FREE_SEED_GAS_FALLBACK: {
            "name": "dry_dispatcher_film_free_seed_order_stored_then_gas",
            "order": "the absent-film set of the first film-free advance after depletion tried "
            "at the stored surface first; on the conservative absent elimination's typed "
            f"refusal {ABSENT_FILM_CHART_REFUSAL!r} the same set tried at the former bath gas "
            "state seed, on the same advance and shared allowance, before the contact and "
            "liquid sets",
            "companion_of": FILM_FREE_SEED_STORED_SURFACE,
            "unchanged": "the path proof, the active-set order otherwise, every gate, tolerance, "
            "budget and law; any other refusal of the stored seed falls through as before",
            "stamped": "native_report.kernel_envelope_film_free_seed.accepted_seed and "
            "first_seed_refusal",
        },
        DRY_HALF_STEP_FIRST: {
            "name": "dry_half_step_first_for_film_bearing_advances",
            "policy": "a dry advance whose start carries an exterior film abandons its full step "
            "when the shared per-advance assembly counter reaches the threshold and continues "
            "the ordinary halving from the half step on the same allowance",
            "declared_numerics_threshold_assemblies": DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD,
            "threshold_origin": "above the maximum of the accepted film-bearing dry advances of "
            "the 27 complete steam-side cases of evidence film_free_seed_168 (1,102); not a "
            "physics value",
            "unchanged": "the allowance (1,400 with the per-unknown knob off), the 13 halvings, "
            "every gate, tolerance and law; film-free advances never abandon",
            "stamped": "a refusal row kernel_envelope_dry_half_step_first at the full step and "
            "the accepted row's kernel_envelope_dry_half_step_first",
        },
        DEPLETION_CARRIES_CONDENSATE: {
            "name": "eliminated_depletion_event_carries_the_dew_line_condensate",
            "bookkeeping": "the exact-zero water closure of the eliminated dry-stage depletion "
            "event L_old - water_gain - d (j_w - deposition) = 0, the deposition the "
            "projection's booked dew-line condensate (already in the film node's water row "
            "and, with the vapour enthalpy, its energy row)",
            "former": "the typed refusal 'the eliminated depletion event does not carry the "
            "dew-line condensate'",
            "unchanged": "the condensate rate, the projection, its offset, every gate and law; "
            "a film-free surface after the event still refuses booked condensate typed",
            "acts_only_under": GAS_DEW_LINE_PROJECTION,
        },
        INSTRUMENT_SHARED_LAW: {
            "name": "instrument_water_activity_on_the_shared_equilibrium",
            "reading": "a_w on the shared law (the kernel's chart function "
            "coupled_pore.evaluate_equilibrium: PHY-053 second-virial binary gas, immiscible pure "
            "liquids) in place of the ideal immiscible Raoult-Dalton reading (1 - y_h) P / "
            "P_sat,w(T)",
            "report_only": True,
            "stamped": "both readings and their difference",
            "remaining_src_sites": "the tower's I3 (through_bed_condensation_instruments), "
            "sorption_interface.hexane_activity_ceiling, core2/tray.py, the retired core/ path",
        },
    }
    stall_rows = {
        STALL_FILM_DEPLETION_EVENT: {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": f"{STALL_DECISION}; no law, gate, tolerance, budget, bound, band or "
            "declared value moved",
            "record": STALL_RECORD,
            "physically_qualifying": False,
            "name": "stalled_front_film_depletion_event_and_declined_advance_materiality",
            "companion_of": STALLED_FRONT,
            "event": "while the front is held (the last accepted step stalled, or the moving "
            "chart refused at the attempted duration), the stalled step with the film unknown "
            "replaced by the event duration (declared duration bounds of the film node), the "
            "film fixed at exactly zero, the D9 drain held on its capacity branch (transfer "
            "share 1, the front-stage event's construction); the front a literal, the hexane "
            "and energy jump rows solved, the water jump row the q = 0 consistency row at the "
            "unchanged gate; accepted on the native stall's accept chart; seeded at the linear "
            "zero of the last two accepted films",
            "after_event": "film-free; the moving (receding) chart in the gas-only surface "
            "mode is tried first on every step, a stalled step only after it is refused",
            "materiality_row": "per accepted stalled step and stalled event: the counterfactual "
            "advance demand (evidence pd_advance_demand_170's signed front velocity closing the "
            "water, hexane and energy jump rows at the held front, evaluated from the step's "
            "start state over the step in a scratch overlay, never adopted), the declined "
            "hexane q (C_h,wet - C_h,dry) dt per packet, its cumulative value and its ratio "
            "to the particle's hexane; report only",
            "phy_014": "unchanged: the front never advances",
            "unchanged": "every gate, tolerance, budget, bound, band, law and declared value",
        },
    }
    graded_rows = {
        SURFACE_GRADED_MESH: {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": f"{PD_RULING} item 7 (the surface-graded radial ladder is run) and "
            "GT_PS2_PAPER1_SURFACE_LAYER_MESH_PACKET_2026-10-18 section 4 (Class A numerics)",
            "record": GRADED_RECORD,
            "physically_qualifying": False,
            "name": "surface_graded_radial_master_grid",
            "grid": "particle_surface_graded_mesh.SurfaceGradedFaceSet(m, n, q): m uniform "
            "interior cells and n band cells under the surface at geometric ratio q, the "
            "image of the uniform computational grid under one C1 map; built by the native "
            "_from_faces and checked by the native grid validators before use",
            "parameters": "the runner's --graded-interior-cells, --graded-band-cells, "
            "--graded-band-ratio (the case record's surface_graded_mesh section carries the "
            "faces and widths); --cells equals m + n",
            "former": "the uniform grid uniform_grid(N, R)",
            "unchanged": "every gate, tolerance, budget, law and declared value; the per-cell "
            "budgets read the grid's N as before",
        },
    }
    three_knob_rows = {
        BIRTH_SEED_BOUNDED_TO_OUTER_CELL: {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": "the lead's Class-A seed choice (a search seed, no value moved) on the "
            f"measurement {BIRTH_SEED_BOUND_MEASUREMENT} sections 18 to 23 (the birth root "
            "independent of the seed wherever the birth converges)",
            "record": THREE_KNOBS_RECORD,
            "physically_qualifying": False,
            "name": "birth_seed_bounded_to_outer_cell",
            "seed": "the energy-based newborn-shell thickness seed of advance_birth replaced by "
            "min(seed, f x outer cell width) before its fit test (outer cell width R - "
            "faces[-2] of the birth's grid)",
            "declared_fraction": "experiment_core2_dt_birth.BIRTH_SEED_OUTER_CELL_FRACTION, "
            "declared 0.5 on 2026-10-28",
            "fraction_override": birth_seed_outer_cell_fraction_override(environ),
            "companion": PREFIX + BIRTH_SEED_OUTER_CELL_FRACTION,
            "former": "the seed refused when it did not fit the outer material cell (it still is "
            "when it is not positive)",
            "stamped": "the birth audit's kernel_envelope_birth_seed_bound (the fraction in "
            "force, whether the bound acted, the unbounded and bounded seed thicknesses)",
            "unchanged": "the birth solve, its budget, gates, certificate and every law",
        },
        FANER_MEASURED_VAPOUR_TRACE: {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": "the owner's ruling, reported with "
            f"{FANER_TRACE_MEASUREMENT}: Faner's measured vapour temperature trace is the "
            "base-case bath of the Faner march, not a sensitivity",
            "record": THREE_KNOBS_RECORD,
            "physically_qualifying": False,
            "name": "faner_measured_vapour_trace_as_the_bath",
            "trace": "faner_measured_trace: the 17 digitized experimental gas_temperature points "
            "of Faner et al. 2019 figure 2, soybean, 120 C (paper/analysis/datasets/faner2019_v2/"
            "curves.csv, superseded rows excluded)",
            "rule": "linear in time between the points, constant before the first (0 s) and "
            "after the last (235.058 s); only the bulk temperature of the boundary and of its "
            "film boundary replaced",
            "stages": "the feed rates at each attempted interval's end clock, the birth and its "
            "seeds at activation + birth dt, the front and dry boundary factories at their "
            "clocks, the in-run prefix audit at each feed row's clock",
            "former": "the constant --bath-temperature (393.15 K) in every stage",
            "stamped": "the birth record's and the case record's faner_measured_vapour_trace "
            "(the file path, its sha256, the points and the rule)",
            "studies": "Faner only; refused at launch for the DT study",
            "unchanged": "every film coefficient law, gate, tolerance, budget and law",
        },
        FRONT_STAGE_NOISE_AWARE_ENERGY_GATE: {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": f"{FRONT_GATE_PACKET} decision 1, open: the owner's ruling "
            "GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17 point 2 extended to the front "
            f"stage (Class B); built off by default on the measurement {FRONT_GATE_MEASUREMENT} "
            "section 12 (variant A); NOT adopted, its arming in a reported history waits on the "
            "owner's ruling",
            "record": THREE_KNOBS_RECORD,
            "physically_qualifying": False,
            "owner_item": "the extension of the noise-aware energy-row gate to the front stage",
            "name": "front_stage_noise_aware_energy_row_gate",
            "rows": "the front-stage film-node search's dry-piece (shell) energy rows, the "
            "Rankine-Hugoniot (jump) energy row and the surface node's energy row",
            "gate": "max(target, m x floor_i), floor_i = sum_j |J_ij| numpy.spacing(|x_j|), J "
            "the search's physical Jacobian, x the physical coordinates (the temperature "
            "origins added back as physical_chart_solve translates them)",
            "declared_multiple": "core2_water_film_node.FRONT_STAGE_ENERGY_GATE_MULTIPLE (2)",
            "lower_bounds": "the front search target in the search (1e-11 as run), the native "
            "nonlinear residual tolerance in the post-search check (1e-10)",
            "other_rows": "every water and hexane row at its unchanged search target and native "
            "tolerance",
            "stamped": "per accepted solve the envelope audit's front_stage_noise_aware_energy_"
            "gate (rows, residuals, floors, gates, absorbed-by-floor flags); on a refusal the "
            "error's study_debug",
            "unchanged": "the native particle ledger, the condition certificate, the front "
            "driver's conservation check, every other gate, tolerance, budget and law",
        },
    }
    lever_common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "the lead's Class-A numerics in the study lane (the delegation charter of "
        "2026-08-11), on the measurement "
        f"{JACOBIAN_LEVERS_MEASUREMENT}; no gate, tolerance, cap, allowance or law moved",
        "record": JACOBIAN_LEVERS_RECORD,
        "physically_qualifying": False,
    }
    lever_rows = {
        JACOBIAN_BROYDEN_UPDATES: {
            "name": "newton_search_jacobian_broyden_updates_with_stall_rebuild",
            "where": "every domain_aware_newton solve inside the front driver's front-stage scope "
            "and the runner's dry-stage scope (particle_numerical_repairs.JacobianLeverStage); "
            "the feed, birth, core endpoint and prefix audit unchanged",
            "search_jacobian": "held within a solve and updated after every accepted step by "
            "J + (dr - J s) s^T / (s^T s); a fresh dense difference Jacobian at the first "
            "iterate, at the next iterate after an accepted step with merit_after > r x "
            "merit_before, and at the same iterate when the held direction finds no admissible "
            "decreasing step",
            "declared_ratio": "particle_numerical_repairs.BROYDEN_REBUILD_CONTRACTION_RATIO, "
            "declared 0.5 on 2026-10-31",
            "former": "a fresh dense difference Jacobian at every iterate",
            "stamped": "each accepted history row's kernel_envelope_jacobian_levers (per solve: "
            "the calling site, the Jacobian source of every iteration, the update and rebuild "
            "counts, the outcome); the stage totals in the record",
            "unchanged": "the convergence test (gtol and the acceptance clause), max_nfev and its "
            "counting, the dense certificate Jacobian at every returned or provisional iterate, "
            "the line search and its merit, the damping ladder, the bounds handling, every gate, "
            "tolerance, cap and allowance",
        },
        JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES: {
            "name": "dry_stage_certificate_jacobian_carried_to_the_next_solve",
            "where": "the runner's dry stage only; refused for the front stage",
            "search_jacobian": "a solve's first search Jacobian is the previous converged "
            "solve's dense certificate Jacobian of the same layout (the calling site, the unknown "
            "count, the acceptance class and the difference mode), then Broyden as declared by "
            f"{PREFIX + JACOBIAN_BROYDEN_UPDATES}",
            "companion": PREFIX + JACOBIAN_BROYDEN_UPDATES,
            "former": "a fresh dense difference Jacobian at a solve's first iterate",
            "stamped": "the 'carried' source of a solve's first iteration in the row's "
            "kernel_envelope_jacobian_levers",
            "unchanged": "as for the Broyden updates",
        },
        PATH_PROOF_BRANCH_MEMO: {
            "name": "composition_path_proof_interval_branch_exact_memo",
            "where": "native_water_path_bounds.water_path_derivative_bound, every stage",
            "memo": "keyed on the exact IEEE-754 bits of the six inputs the bound reads; a hit "
            "returns a deep copy of the identical result or raises the identical domain refusal; "
            "an ArithmeticError is not memoized; bounded LRU",
            "declared_bound": "native_water_path_bounds.PATH_PROOF_BRANCH_MEMO_BOUND, declared "
            "65536 on 2026-10-31",
            "former": "the bound evaluated at every call",
            "stamped": "each accepted history row's kernel_envelope_path_proof_branch_memo (the "
            "cumulative hits, misses, evictions and entries)",
            "unchanged": "every proof verdict, leaf and point count, every value",
        },
    }
    out = []
    for item in armed_items(environ):
        if item == CORE_TARGET_SEED:
            out.append(core_target_seed_declared_row(valued(item, environ)))
        elif item == BUDGETS_PER_UNKNOWN:
            out.append(budgets_per_unknown_declared_row(valued(item, environ)))
        elif item == LIQUID_MOBILITY_RULE:
            out.append(surface_rules_declared_row(item, valued(item, environ)))
        elif item == SURFACE_LAW_HISTORICAL:
            out.append(surface_law_historical_declared_row(valued(item, environ)))
        elif item == FEED_CELLS:
            out.append(feed_cells_declared_row(int(valued(item, environ))))
        elif item == TOWER_FRONT_TIME_LEVEL:
            out.append(tower_front_time_level_declared_row(int(valued(item, environ))))
        elif item == TOWER_DRY_TIME_LEVEL:
            out.append(
                tower_dry_time_level_declared_row(
                    tower_dry_time_level(environ), tower_dry_time_level_guarded(environ)
                )
            )
        elif item == TOWER_CORE_EXTINCTION_FLOOR_M:
            out.append(tower_core_extinction_floor_declared_row(valued(item, environ)))
        elif item == JACOBIAN_COLUMN_GROUPS:
            out.append(jacobian_column_groups_declared_row(valued(item, environ)))
        elif item == DEFERRED_CERTIFICATE:
            out.append(deferred_certificate_declared_row())
        elif item in three_knob_rows:
            out.append(three_knob_rows[item])
        elif item in lever_rows:
            out.append({**lever_common, **lever_rows[item]})
        elif item == HEXANE_EXHAUSTED_REGIME:
            out.append(hexane_exhausted_regime_declared_row())
        elif item in (JACOBIAN_LEVERS_AT_BIRTH, DRY_FACE_FLUX_MEMO):
            out.append(birth_cost_lever_declared_row(item))
        elif item == IAPWS_COMPILED:
            out.append(iapws_compiled_declared_row())
        elif item == STEAM_QUENCH_ARRIVAL:
            out.append(steam_quench_arrival_declared_row())
        elif item == HEXANE_EXHAUSTED_REENTRY_HYSTERESIS:
            out.append(hexane_exhausted_reentry_hysteresis_declared_row())
        elif item == STEAM_CONDENSING_FILM:
            out.append(steam_condensing_film_declared_row())
        elif item in (
            STEAM_CONDENSING_FILM_SEEDED_START,
            STEAM_CONDENSING_FILM_FLOOR_EXIT,
            STEAM_CONDENSING_FILM_REQUEST_MEMORY,
            STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP,
            DRY_BUDGET_HALVING,
            FILM_FREE_WET_SET_SKIP,
            ABSENT_SEED_ABOVE_DEW,
        ):
            out.append(condensing_film_newton_lever_declared_row(item))
        elif item in graded_rows:
            out.append(graded_rows[item])
        elif item in stall_rows:
            out.append(stall_rows[item])
        elif item in queue_rows:
            out.append({**queue_common, **queue_rows[item]})
        elif item in seed_rows:
            out.append(seed_rows[item])
        elif item in completion_rows:
            out.append({**completion_common, **completion_rows[item]})
        elif item in pd_rows:
            out.append({**pd_common, **pd_rows[item]})
        else:
            out.append({**common, **rows[item]})
    return out


def hexane_exhausted_reentry_hysteresis_declared_row():
    """HEXANE_EXHAUSTED_REENTRY_HYSTERESIS's typed declared row (2026-11-03)."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "name": "hexane_exhausted_regime_reentry_hysteresis",
        "authority": "the steam-quench build's standing mandate of 2026-11-03 (a study-lane extension "
        "of the owner-ruled n-hexane-exhausted regime, on the lead's authority)",
        "record": STEAM_QUENCH_ARRIVAL_RECORD,
        "physically_qualifying": False,
        "reentry_gas_y_hexane": HEXANE_EXHAUSTED_REENTRY_GAS_Y_HEXANE,
        "reentry_was": "above the clamp HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP",
        "between": "the vertex continues with the presented gas taken as hexane-free (the vertex "
        "stage's own stamp); the ignored y stamped on the advance",
    }


def condensing_film_newton_lever_declared_row(item):
    """The two condensing-film Newton levers' typed declared rows (2026-11-07)."""
    rows = {
        STEAM_CONDENSING_FILM_SEEDED_START: {
            "name": "steam_condensing_film_seeded_start",
            "rule": "a regime solve continuing the previous converged regime solve of the same "
            "job starts from its state plus the previous increment x (this duration / the "
            "previous duration) when that residual is admissible and smaller than the former "
            "start's",
            "former": "every regime solve from the advance's former initial point",
        },
        STEAM_CONDENSING_FILM_FLOOR_EXIT: {
            "name": "steam_condensing_film_floor_exit",
            "rule": "the regime's first Newton handed to its surface-energy polish once its best "
            "largest scaled residual lies within the declared factor of gtol and has not "
            "halved over the declared patience",
            "former": "the first Newton run to its allowance before the retry and the polish",
        },
        STEAM_CONDENSING_FILM_REQUEST_MEMORY: {
            "name": "steam_condensing_film_request_memory",
            "rule": "after a regime advance accepted below the regime's cap, the next request "
            "is at most twice the accepted duration (doubling back to the cap)",
            "former": "every sub-step requests the cap again",
        },
        STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP: {
            "name": "steam_condensing_film_film_free_request_cap",
            "rule": "a film-free binary dry state under a superheated hexane-free gas hotter "
            "than its surface (the condensing trigger not holding) requests at most the dry "
            "schedule's film-present band value, 0.3125 s",
            "former": "the film-free bands, up to 2.5 s, as one sub-step",
        },
        DRY_BUDGET_HALVING: {
            "name": "dry_budget_halving",
            "rule": "a dry sub-step attempt that exhausts its residual-assembly allowance enters "
            "the dry halving ladder, the next halving on a fresh allowance of the same size, at "
            "most DRY_BUDGET_HALVING_MAX (3) times per sub-step",
            "former": "StepBudgetExceeded ends the job",
            "changed": "the allowance's scope: its size per attempt unchanged, but the "
            "assemblies spent on one sub-step across attempts may exceed it, up to "
            "(1 + DRY_BUDGET_HALVING_MAX) allowances; each retry stamped with its attempt "
            "count and the assemblies spent on the sub-step",
            "maximum_budget_halvings_per_substep": DRY_BUDGET_HALVING_MAX,
        },
        FILM_FREE_WET_SET_SKIP: {
            "name": "film_free_wet_set_skip",
            "rule": "a film-free dry advance whose surface lies more than "
            "FILM_FREE_WET_SET_SKIP_MARGIN_K (5 K) above the presented gas's water dew point at "
            "both ends of the span tries the absent set only",
            "former": "absent, then contact, then liquid (or contact first after a wet set)",
            "margin_k": FILM_FREE_WET_SET_SKIP_MARGIN_K,
        },
        ABSENT_SEED_ABOVE_DEW: {
            "name": "absent_seed_above_dew",
            "rule": "a film-free advance whose absent-set seed maps its zero-flux interface "
            "outside the gas-only chart is seeded at the gas temperature plus "
            "ABSENT_SEED_ABOVE_DEW_OFFSET_K (0.01 K) with the gas composition; the set order "
            "unchanged",
            "former": "the stored (or gas-state) seed, refused at its first evaluation",
            "offset_k": ABSENT_SEED_ABOVE_DEW_OFFSET_K,
        },
    }
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        **rows[item],
        "record": CONVERGENCE_ROOT_CAUSES_RECORD,
        "unchanged": "the residual, gtol, the acceptance clause, the certificate, the polish and "
        "its native re-certification, the allowances (every extra evaluation charged), the beta "
        "domain",
        "physically_qualifying": False,
    }


def steam_condensing_film_declared_row():
    """STEAM_CONDENSING_FILM's typed declared row (2026-11-06)."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "name": "dry_stage_true_gas_condensing_film_branch",
        "authority": "the owner's ruling of 2026-11-04 (the scoped true-gas steam-condensation branch) "
        f"on the reviewed candidate {STEAM_CONDENSING_FILM_REVIEW}; integrated by the tower lane's "
        "standing build on the lead's brief of 2026-11-06",
        "ruling": STEAM_CONDENSING_FILM_RULING,
        "record": STEAM_CONDENSING_FILM_RECORD,
        "physically_qualifying": False,
        "gas_y_hexane_threshold": STEAM_CONDENSING_FILM_GAS_Y_HEXANE,
        "threshold_was": "no regime",
        "trigger": "at the start of a dry advance: a binary dry state; the presented gas below the "
        "threshold in n-hexane, superheated ((1 - y) P / p_sat(T_gas) below 1) and hotter than the "
        "stored surface at both ends of the span; the stored surface below the gas's water dew point "
        "on the kernel's gas-only chart (the gas leaner than the chart's lower bound at the surface "
        "temperature) or an exterior film present",
        "search": "a finite inward search box for the native Stefan inverse, max(3, rho h_M "
        "(|ln w_bulk| + 1) max(cp) / h_Q), from the advance's native coefficients; the native inverse "
        "first at the caller's limit clipped to the original band; one ulp-neighbour Newton retry and "
        "a single-coordinate surface-energy polish re-certified by a fresh native Newton, all charged "
        "to the original search and shared allowances; a film-free particle's liquid set from exactly "
        "zero film, then the former dispatch",
        "accepted_domain": "|beta| <= the boundary's original limit (3); the external-film profile of "
        "every accepted root audited",
        "maximum_dry_advance_s": STEAM_CONDENSING_FILM_MAXIMUM_DT_S,
        "maximum_dry_advance_was": "the runner's dry schedule",
        "positive_film_stop": "lifted for every eliminated liquid advance under a gas at the bath "
        "clamp (the reviewed candidate's rule; the measured alternatives are stamped in the record)",
        "carried_jacobian": "dropped at the start of every regime advance (a fresh finite-difference "
        "Jacobian, Broyden updates inside the advance unchanged), as each of the reviewed candidate's "
        "one-advance jobs started; WAS the previous advance's certificate (inside a 2 s job the carried "
        "Jacobian left the surface-energy row stagnating for 742 to 876 assemblies per solve)",
        "unchanged": "every law, gate, tolerance, allowance, budget, presented gas and the accepted "
        "beta domain; no clamp, film dose or state projection",
    }


def steam_quench_arrival_declared_row():
    """STEAM_QUENCH_ARRIVAL's typed declared row (2026-11-03)."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "name": "dry_stage_steam_quench_condensation_to_saturation_arrival_event",
        "authority": "the lead's decision of 2026-11-03 under the owner's direction of the same day "
        "(a study-lane engineering knob of the tower lane, Class B, built on the lead's authority "
        f"as a stamped interim knob), on the diagnosis {STEAM_QUENCH_ARRIVAL_DIAGNOSIS}",
        "record": STEAM_QUENCH_ARRIVAL_RECORD,
        "physically_qualifying": False,
        "gas_y_hexane_threshold_film": STEAM_QUENCH_GAS_Y_HEXANE_FILM,
        "gas_y_hexane_threshold_film_free": STEAM_QUENCH_GAS_Y_HEXANE_FILM_FREE,
        "thresholds_was": "no regime",
        "trigger": "at the start of a dry advance: a binary dry state; its stored surface below "
        "props.water.saturation_temperature at the state's pressure; the presented gas above that "
        "temperature, with a water activity (1 - y) P / p_sat(T_gas) below 1 and an n-hexane "
        "fraction below the threshold (film-bearing or film-free) at both ends of the span",
        "event": "steam from the presented gas condenses until the declared heated depth, the outer "
        "cell and the exterior film, stands at the saturation temperature (the outer cell only ever "
        "raised, its n-hexane inventory held by its pore-gas y); the condensed mass in closed form "
        "from the kernel's own caloric (coupled_pore's equilibrium storage, props.water's liquid "
        "internal energy for the film, its vapour enthalpy at the presented gas temperature); the "
        "combined storage change booked as one inward stream in the stream integral and the "
        "original-feed ledger; the binary liquid set then advances under the cell's own gas",
        "heated_depth": "the outer cell and the exterior film",
        "heated_depth_was": "no regime",
        "unchanged": "every law, gate, tolerance, allowance and presented gas; no clamp",
    }


def iapws_compiled_declared_row():
    """IAPWS_COMPILED's typed declared row (2026-11-03)."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "name": "iapws95_residual_compiled_bit_identical",
        "authority": "the lead's brief of 2026-11-03 (hint C of the external parallelization "
        "study); a measurement and a prototype, not adopted; a core property law's evaluation, "
        "so its adoption is the owner's even in the exact case",
        "record": IAPWS_COMPILED_RECORD,
        "where": "props.water._residual_uncached replaced in-process (the module attribute; no src "
        "file edited) by particle_iapws_compiled.compiled_residual_uncached; reached through the "
        "exact-bit residual cache by every water property of the process",
        "law": "unchanged: the same Python source traced operation by operation into LLVM IR "
        "(no optimisation pass, no fast-math flag, generic x86-64 target without FMA); exp, log "
        "and pow bound by address to ucrtbase.dll, the runtime the interpreter imports; CPython's "
        "float_pow special cases mirrored",
        "fallback": "the Python residual itself for a non-float argument, a zero divisor, a "
        "transcendental or pow outside its finite non-raising domain, or a non-finite output",
        "former": "the interpreter's evaluation of the same 56-term source",
        "unchanged": "every value (bit for bit), the exact-bit caches and their keys, every gate, "
        "tolerance, cap, budget and allowance",
        "physically_qualifying": False,
    }


def birth_cost_lever_declared_row(item):
    """The typed declared rows of the kernel birth-cost levers of 2026-11-02."""
    common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "the lead's Class-A numerics in the study lane (the delegation charter of "
        "2026-08-11), on the record's own profile of a tower-born particle; no gate, tolerance, "
        "cap, allowance or law moved",
        "record": BIRTH_COST_LEVERS_RECORD,
        "physically_qualifying": False,
    }
    if item == JACOBIAN_LEVERS_AT_BIRTH:
        return {
            **common,
            "name": "birth_newton_search_jacobian_broyden_updates",
            "where": "the birth's one domain-aware Newton solve (the Paper-1 runner's "
            "fresh_dt_birth, the cohort tower's born job; particle_numerical_repairs "
            "JacobianLeverStage 'birth'); the carried companion refused there; the feed stage "
            "(scipy least_squares) unchanged",
            "search_jacobian": f"as declared by {PREFIX + JACOBIAN_BROYDEN_UPDATES} for the front "
            "stage (particle_numerical_repairs.BROYDEN_REBUILD_CONTRACTION_RATIO)",
            "companion": PREFIX + JACOBIAN_BROYDEN_UPDATES,
            "former": "a fresh dense difference Jacobian at every iterate of the birth's solve",
            "stamped": "the birth record's kernel_envelope_jacobian_levers (per solve: the calling "
            "site, the unknown count, the Jacobian source of every iteration, the update and "
            "rebuild counts, the outcome) and the stage totals",
            "unchanged": "the convergence test (gtol and the acceptance clause), max_nfev and its "
            "counting, the dense certificate Jacobian at every returned or provisional iterate, "
            "the birth condition certificate, the line search and its merit, every gate, "
            "tolerance, cap and allowance",
        }
    return {
        **common,
        "name": "dry_face_flux_exact_memo_per_solve",
        "where": "cut_transport._dry_flux as reached through core.cut._dry_flux (the pinned file "
        "unedited), wrapped at the innermost level; one frame per domain-aware Newton solve",
        "memo": "keyed on the exact bits of every argument and the identity of the harness's "
        "in-memory substitutions in the flux's call tree and the K1 suppression depth; a hit "
        "returns a deep copy of the identical result and replays the proof allowance's "
        "accounting of the call (points started, memoized points, the last path and point); a "
        "hit that would cross the allowance's point clause, or any hit under the wall-clock "
        "clause, runs the flux; a call that raises is not memoized",
        "declared_bound": "particle_face_flux_memo.FACE_FLUX_MEMO_BOUND per frame, declared 4096 "
        "on 2026-11-02",
        "former": "the flux evaluated at every call",
        "stamped": "each accepted history row's kernel_envelope_dry_face_flux_memo (the cumulative "
        "calls, hits, misses, unkeyable calls, declined hits, evictions, frames and replayed "
        "points) and the record's last",
        "unchanged": "every flux value, proof verdict and point count, every gate, tolerance, cap "
        "and allowance",
    }


def tower_dry_time_level(environ=None):
    """TOWER_DRY_TIME_LEVEL as an int when set (one of its declared values, the guard suffix
    removed), None when unset."""
    value = valued(TOWER_DRY_TIME_LEVEL, environ)
    return None if value is None else int(value.removesuffix(TOWER_DRY_GUARD_SUFFIX))


def tower_dry_time_level_guarded(environ=None):
    """TOWER_DRY_TIME_LEVEL carries the depletion guard (the declared value "-1+guard")."""
    value = valued(TOWER_DRY_TIME_LEVEL, environ)
    return value is not None and value.endswith(TOWER_DRY_GUARD_SUFFIX)


def tower_dry_band_factor(level):
    """The factor on every cap of the lane worker's dry bands (its runner bands at L = 1) at a
    declared level: 2**(1 - L), an exact power of two (1 at level 1, every cap bit for bit the
    ruled one)."""
    if str(level) not in TOWER_DRY_TIME_LEVEL_VALUES:
        raise KeyError(f"undeclared tower dry time level {level!r}")
    return 2 ** (1 - level)


def tower_dry_time_level_declared_row(level, guarded=False):
    """TOWER_DRY_TIME_LEVEL's typed declared row (with the depletion guard's rule when guarded)."""
    factor = tower_dry_band_factor(level)
    guard = {
        "depletion_guard": {
            "coarse_steps": TOWER_DRY_GUARD_COARSE_STEPS,
            "rule": "coarse bands for a job only in the front stage, or in the dry stage with the "
            "external film present and its predicted depletion (curvature_corrected_depletion_"
            "duration of the state's film points; unknown counts as far) more than the job's span "
            "plus coarse_steps coarse film-bearing steps ahead; else the ruled bands",
        }
    }
    return {
        **(guard if guarded else {}),
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "a study-lane measurement knob on the lead's brief of 2026-11-02 for the "
        "owner's decision whether the tower lane may run the dry stage at coarser declared bands; "
        "nothing adopted; Paper 1 untouched",
        "record": TOWER_DRY_TIME_LEVEL_RECORD,
        "physically_qualifying": False,
        "name": "tower_dry_stage_time_level",
        "level": level,
        "band_factor": factor,
        "band_before_first_landmark_s": 0.15625 * factor,
        "film_bearing_band_s": 0.3125 * factor,
        "film_free_band_after_second_landmark_s": 2.5 * factor,
        "former": "time level 1, the lane worker's DRY_TIME_CONTROL (the runner's "
        "dt-film-resolved bands at L = 1: 0.15625 s, 0.3125 s, then 0.3125 s with the external "
        "film or 2.5 s without)",
        "where": "the cohort tower's advance jobs only (every dry segment of the lane worker's "
        "advance); refused by the Paper-1 runner and by the lane worker alone",
        "unchanged": "the two landmarks (on the birth origin), the film-present rule, the 13 "
        "halvings, the depletion advance allowance factor, the half-step-first rule, the Newton, "
        "its allowances, every gate, tolerance, certificate and ledger; the front stage, the core "
        "endpoint, the feed stage and the birth",
    }


def tower_front_time_level(environ=None):
    """TOWER_FRONT_TIME_LEVEL as an int when set (one of its declared values), None when unset."""
    value = valued(TOWER_FRONT_TIME_LEVEL, environ)
    return None if value is None else int(value)


def tower_front_maximum_dt_s(level):
    """The front-stage step cap at a declared level: the lane worker's own expression 0.1 / 2**L
    (particle_kernel_lane_worker.FRONT_MAXIMUM_DT_S at L = 1, bit for bit)."""
    if str(level) not in TOWER_FRONT_TIME_LEVEL_VALUES:
        raise KeyError(f"undeclared tower front time level {level!r}")
    return 0.1 / 2**level


def tower_front_time_level_declared_row(level):
    """TOWER_FRONT_TIME_LEVEL's typed declared row."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "a study-lane measurement knob on the lead's brief of 2026-11-02 for the "
        "owner's decision whether the tower lane may run the front stage at a coarser declared "
        "level; nothing adopted; Paper 1 untouched",
        "record": TOWER_FRONT_TIME_LEVEL_RECORD,
        "physically_qualifying": False,
        "name": "tower_front_stage_time_level",
        "level": level,
        "front_maximum_dt_s": tower_front_maximum_dt_s(level),
        "former": "time level 1, the front-stage cap 0.05 s (the lane worker's TIME_LEVEL)",
        "where": "the cohort tower's born and advance jobs only (the newborn's front driver "
        "after the birth and every front segment of the lane worker's advance); refused by the "
        "Paper-1 runner and by the lane worker alone",
        "unchanged": "the front driver's minimum interval and its controller (next step at most "
        "1.5 times the last, never above the cap), the core-time tolerance 0.001 s, the feed-stage "
        "and birth steps, the dry-stage schedule, the Newton, its allowances, every gate, "
        "tolerance, certificate and ledger",
    }


# ---- TOWER_CORE_EXTINCTION_FLOOR_M: the tower lane's core-extinction floor (2026-11-02) ----


def tower_core_extinction_floor_m(environ=None):
    """TOWER_CORE_EXTINCTION_FLOOR_M's front radius in metres when armed, None when unset.

    Read through a live envelope activation: the whole environment is validated first.
    """
    validate_environment(environ)
    value = valued(TOWER_CORE_EXTINCTION_FLOOR_M, environ)
    return None if value is None else float(value)


def tower_core_extinction_floor_declared_row(value):
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of TOWER_CORE_EXTINCTION_FLOOR_M."""
    if value not in VALUED_ITEMS[TOWER_CORE_EXTINCTION_FLOOR_M]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + TOWER_CORE_EXTINCTION_FLOOR_M}={value!r} is not a declared value"
        )
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "measurement knob of the lead's brief of 2026-11-02 (the owner's question "
        "'so now we know that bottleneck is core extinction. What do we do about it?'); nothing "
        "adopted, the owner's decision",
        "record": TOWER_CORE_EXTINCTION_FLOOR_RECORD,
        "physically_qualifying": False,
        "name": "tower_core_extinction_floor_front_radius_m",
        "knob": f"{PREFIX + TOWER_CORE_EXTINCTION_FLOOR_M}={value}",
        "floor_front_radius_m": float(value),
        "terminal_volume_ceiling": "max(1e-12, (r_f / R)^3), R the particle radius",
        "terminal_volume_ceiling_was": "1e-12 (particle_core_extinction.estimate_core_time, the "
        "predeclared four-target 0.82-contraction endpoint rule)",
        "remaining_time_tolerance": "inside the floor's band (1e-12 < z <= (r_f / R)^3) the "
        "conditional remaining-time bound is stamped, not enforced",
        "remaining_time_tolerance_was": "the conditional remaining-time bound at or below the "
        "caller's tolerance (0.002 / 2^L; 0.001 s in the tower lane) at every endpoint",
        "event": "at the first accepted positive target in the core cell whose wet volume fraction "
        "lies at or below the ceiling (four linked targets, every ratio below 0.82): the core "
        "endpoint at once (particle_core_extinction.project_core_endpoint, the native fixed-time "
        "extinction projection over the estimated remaining time), the wet core's water and "
        "n-hexane passing to the dry pore gas and sorbed phase in that native step, its combined "
        "ledger checked at the unchanged 1e-10 gate; in the tower lane the endpoint's storage "
        "change booked in the job's stream as for the ruled endpoint",
        "stamped": "the driver record's kernel_envelope_core_extinction_floor and the worker "
        "endpoint segment's core_extinction_floor (z, radius, clock, the wet core's water and "
        "n-hexane released, the four durations, the estimate, the bound, the tolerance)",
        "unchanged": "the half-volume targets and their search, the four-target linkage, the 0.82 "
        "contraction, the native projection, its residual contract, its budget and ledgers, the "
        "clock-lead limit, every law; no state clamped or floored but the declared z",
        "refused_by": "the Paper-1 runner (run_particle_qualification_case) and the lane worker "
        "alone (verify_program_knobs)",
    }


def tower_ruled_time_resolution_declared_row():
    """The ruled tower-lane time resolution's typed declared row (the "+ruled" cohort knob sets)."""
    level = int(TOWER_RULED_FRONT_TIME_LEVEL)
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": 'owner ruling 2026-11-02 ("Ok, implement step cap 0.2 s and an extintion '
        'flooe 1 micrometer. Report back speed improvement"), on '
        f"{TOWER_FRONT_TIME_LEVEL_RECORD} and {TOWER_CORE_EXTINCTION_FLOOR_RECORD}",
        "record": TOWER_RULED_TIME_RESOLUTION_RECORD,
        "physically_qualifying": False,
        "name": "tower_ruled_time_resolution",
        "front_time_level": level,
        "front_maximum_dt_s": tower_front_maximum_dt_s(level),
        "front_time_level_was": "1 (design memo 6.1; the front-stage cap 0.05 s)",
        "core_extinction_floor_m": float(TOWER_RULED_CORE_EXTINCTION_FLOOR_M),
        "core_extinction_floor_was": "none (the ruled 1e-12 wet-volume ceiling with the 1 ms "
        "remaining-time bound)",
        "knobs": {
            PREFIX + TOWER_FRONT_TIME_LEVEL: TOWER_RULED_FRONT_TIME_LEVEL,
            PREFIX + TOWER_CORE_EXTINCTION_FLOOR_M: TOWER_RULED_CORE_EXTINCTION_FLOOR_M,
        },
        "where": "the cohort tower's knob sets ending '+ruled' only; the Paper-1 runner and the "
        "lane worker alone refuse both knobs",
    }


def hexane_exhausted_regime_declared_row():
    """HEXANE_EXHAUSTED_REGIME's typed declared row."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "the owner's Class-B ruling of 2026-11-01 (\"I think this 'ghost' regime for "
        "hexane free vapors is completely fine. Let's build it.\"), on the diagnosis "
        f"{HEXANE_EXHAUSTED_REGIME_DIAGNOSIS}",
        "record": HEXANE_EXHAUSTED_REGIME_RECORD,
        "physically_qualifying": False,
        "name": "dry_stage_hexane_exhausted_pure_water_vertex",
        "exhaustion_loading_kg_kg": HEXANE_EXHAUSTED_LOADING_KG_KG,
        "exhaustion_loading_was": "no regime",
        "gas_clamp_y_hexane": HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP,
        "gas_clamp_was": "no regime",
        "trigger": "at the start of a dry advance: no exterior film, total n-hexane per kg dry "
        "matter below the exhaustion loading, the presented gas's n-hexane at or below the clamp "
        "at both ends of the span",
        "entry": "one stamped release event: the binary state converted cell by cell to the "
        "pure-water vertex at the same temperatures and pressure, the storage difference (water, "
        "n-hexane, common-datum energy by the kernel's own caloric) booked to the gas stream and "
        "the cumulative ledger",
        "vertex": "dtdc_simulator.core2.particle.pure_water_pore_gas: water and energy rows only, "
        "y_hexane exactly 0, the presented gas taken as hexane-free (stamped), the independent "
        "diffusive flux zero by the absence of a component",
        "reentry": "a presented gas above the clamp: one stamped re-entry event, every cell's pore "
        "gas at the gas's own y_hexane and the sorbed n-hexane at the frozen isotherm's equilibrium, "
        "the dose taken from the gas stream and booked; the binary dry stage then unchanged",
        "refusal": "particle_hexane_exhausted_regime.KernelHexaneExhaustedRegimeGasHexaneReturned "
        "when the binary chart refuses the re-entry state",
        "unchanged": "every binary dry advance, the Newton, its certificate, the gate units, the "
        "allowances, the 1e-10 ledgers; no state clamped or floored",
    }


def dry_minimum_calls(unknowns):
    """Item 7's measured minimum residual assemblies of an accepted dry advance: 4 u + 6."""
    u = int(unknowns)
    if u != unknowns or u < 1:
        raise ValueError("a positive integer unknown count is required")
    return DRY_MINIMUM_CALLS_PER_UNKNOWN * u + DRY_MINIMUM_CALLS_OFFSET


def dry_assembly_allowance(unknowns, environ=None):
    """The dry advance's residual-assembly allowance.

    Knob off: the former fixed 1,400 (the argument is not read). Armed: three times the
    measured minimum, 3 (4 u + 6).
    """
    if not armed(DRY_ALLOWANCE_PER_UNKNOWN, environ):
        return FORMER_DRY_ASSEMBLY_ALLOWANCE
    return DRY_ALLOWANCE_MULTIPLE_OF_MINIMUM * dry_minimum_calls(unknowns)


def dry_physical_advance_unknowns(cells):
    """u = 3 N + 1 for a dry grid of N cells (the measured minimum's basis)."""
    n = int(cells)
    if n != cells or n < 1:
        raise ValueError("a positive integer cell count is required")
    return 3 * n + 1


def native_ledger_factor_stamp(factors, time_s=None, environ=None):
    """Report-only stamp when a native-ledger factor exceeds the declared ceiling (else None).

    Never a refusal; None when the knob is off.
    """
    if not armed(NATIVE_LEDGER_FACTOR_CEILING, environ):
        return None
    values = [float(f) for f in factors]
    if not values or max(values) <= NATIVE_LEDGER_FACTOR_REPORTING_CEILING:
        return None
    return {
        "kind": "REPORT_ONLY_DECLARED_CEILING_EXCEEDED",
        "quantity": "noise_aware_energy_gate_native_ledger_factor",
        "ceiling": NATIVE_LEDGER_FACTOR_REPORTING_CEILING,
        "maximum_factor": max(values),
        "rows_above_ceiling": [
            k for k, v in enumerate(values) if v > NATIVE_LEDGER_FACTOR_REPORTING_CEILING
        ],
        "time_s": time_s,
        "refusal": False,
        "record": COMPLETION_RECORD,
        "physically_qualifying": False,
    }


def bath_pressure_pa(bath, t_s):
    """Linear bath pressure at t_s (Pa); the bath's own convention for T and y."""
    waypoints = bath["waypoints"]
    carried = ["pressure_pa" in w for w in waypoints]
    if not any(carried):
        return DEFAULT_PRESSURE_PA
    if not all(carried):
        raise ValueError("bath pressure must be carried by every waypoint or by none")
    if t_s <= waypoints[0]["t_s"]:
        value = float(waypoints[0]["pressure_pa"])
    elif t_s >= waypoints[-1]["t_s"]:
        value = float(waypoints[-1]["pressure_pa"])
    else:
        value = None
        for left, right in zip(waypoints, waypoints[1:]):
            t0, t1 = float(left["t_s"]), float(right["t_s"])
            if t0 <= t_s <= t1:
                w = 0.0 if t1 == t0 else (t_s - t0) / (t1 - t0)
                p0, p1 = float(left["pressure_pa"]), float(right["pressure_pa"])
                value = p0 + w * (p1 - p0)
                break
        if value is None:
            raise RuntimeError(f"bath pressure interpolation missed t={t_s}")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("bath pressure must be positive and finite")
    return value


def surface_temperature_upper_k(bath_temperature_k, particle_temperatures_k):
    """Item 4's bound: max(bath, hottest particle temperature) + declared margin."""
    hottest = max(float(t) for t in particle_temperatures_k)
    return max(float(bath_temperature_k), hottest) + SURFACE_TEMPERATURE_MARGIN_K


def two_liquid_root_k(pressure_pa):
    """The native transparent two-liquid interface temperature at a pressure (K).

    The same root the feed stage takes (qsc_k2_native_source_kernel._transparent_phase_root
    with its continuation seed (334.5 K, 0.8)).
    """
    from dtdc_simulator.core2 import qsc_k2_native_source_kernel as kernel

    temperature, _ = kernel._transparent_phase_root(
        pressure_pa=pressure_pa, continuation=(334.5, 0.8)
    )
    return float(temperature)


def linear_film_zero_duration(points):
    """Duration from the last accepted film to the zero of the line through the last two.

    points: [(time_s, film_kg_m2), ...]. None when the last two do not decrease positively.
    Search seeding only.
    """
    if len(points) < 2:
        return None
    (t0, l0), (t1, l1) = points[-2], points[-1]
    if not (t1 > t0 and l0 > 0 and l1 > 0 and l1 < l0):
        return None
    duration = l1 * (t1 - t0) / (l0 - l1)
    return duration if math.isfinite(duration) and duration > 0 else None


# ---- CORE_TARGET_SEED: the positive core target's duration seed from accepted history ----


def core_target_seed(environ=None):
    """CORE_TARGET_SEED's declared value when armed, None when unset.

    Read only through a live envelope activation: the whole environment is validated first, so a
    malformed or unknown DTDC_KERNEL_ENVELOPE_* variable, or a companion without its knob, refuses
    the seed together with every other knob.
    """
    validate_environment(environ)
    return valued(CORE_TARGET_SEED, environ)


def core_target_seed_declared_row(value=CORE_TARGET_SEED_HISTORY_POWER_LAW):
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of CORE_TARGET_SEED at a declared value."""
    if value not in VALUED_ITEMS[CORE_TARGET_SEED]:
        raise EnvelopeKnobRefused(f"{PREFIX + CORE_TARGET_SEED}={value!r} is not a declared value")
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": f"{CORE_TARGET_SEED_DECISION} (Class A numerics: a search seed; no law, gate, "
        "tolerance, budget, seed-reuse threshold, bound, band or declared value moved)",
        "record": CORE_TARGET_SEED_RECORD,
        "evidence": CORE_TARGET_SEED_EVIDENCE,
        "physically_qualifying": False,
        "name": "front_stage_positive_core_target_duration_seed_from_accepted_history",
        "knob": f"{PREFIX + CORE_TARGET_SEED}={value}",
        "seed": "the duration unknown of the positive core target (z -> z/2 in the core cell) seeded "
        "at the halving time of dz/dt = -c z^p through the secants of the last two accepted moving "
        "front steps, each at the geometric mean of its end points: r_k = (z_{k-1} - z_k) / (t_k - "
        "t_{k-1}), m_k = sqrt(z_{k-1} z_k), p = ln(r_1 / r_2) / ln(m_1 / m_2), c = r_2 / m_2^p, "
        "T = (z^(1-p) - (z/2)^(1-p)) / (c (1 - p)), ln 2 / c at p = 1",
        "history": "the last three accepted moving front points (clock, front volume fraction) of "
        "the run, carried in the checkpoint while armed; on a resume from a checkpoint without them "
        "the parent's accepted history rows beside it, used only when the last equals the "
        "checkpoint's state exactly",
        "fallback": "fewer than two positive accepted secants: the driver's seed 0.5 z / (last "
        "secant)",
        "clamp": "none beyond positivity: a fitted seed that is not finite and positive takes the "
        "driver's seed",
        "former": "the driver's seed 0.5 z / (last secant)",
        "unchanged": "the target (half volume), the decision to attempt it, the 2,400-call search "
        "budget and its 800-call cycles, the 1e-5 seed-reuse threshold, the 1e-11 search target, "
        "the 1e-10 ledger gate, the 1e12 condition cap, the per-area chart, the interface "
        "temperature seed, the four-target 0.82-contraction endpoint rule, the minimum interval and "
        "every law",
        "stamped": "front.json kernel_envelope_core_target_seed (the row, the history bootstrap, every "
        "attempt's seed); front.history.jsonl kernel_envelope_core_target_seed on accepted targets",
    }


# ---- FEED_CELLS: the radially resolved DT feed stage (regime A) ----


def feed_cells(environ=None):
    """FEED_CELLS's shell count when armed, None when unset (refused if malformed). Read by the
    DT feed driver only (experiment_core2_dt_feed.run)."""
    value = valued(FEED_CELLS, environ)
    return None if value is None else int(value)


def feed_cells_declared_row(shells):
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of FEED_CELLS."""
    if str(shells) not in FEED_CELLS_VALUES or isinstance(shells, bool):
        raise EnvelopeKnobRefused(f"{PREFIX + FEED_CELLS}={shells!r} is not a declared value")
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "numerical refinement study of the feed stage (formulation audit A10); "
        "nothing adopted",
        "record": FEED_CELLS_RECORD,
        "physically_qualifying": False,
        "name": "radially_resolved_high_loading_feed_stage",
        "knob": f"{PREFIX + FEED_CELLS}={shells}",
        "shells": int(shells),
        "discretization": "the native high-loading state on uniform_grid(N, R): per-shell T_i "
        "with the live X_c(T_i) and the native wet-core energy; attached n-hexane and exterior "
        "water at the outer shell's temperature; Fourier k A_face / (c_(i+1) - c_i) between "
        "shells and A k / (R - c_N) to the interface; backward Euler",
        "declared_choices": "the PD source per kg dry (by shell dry mass); the liquid n-hexane "
        "expelled or taken up by a shell moving to or from the attached inventory at the shell's "
        "own liquid internal energy",
        "hand_over": "the activated shell profile remapped onto the birth's master grid "
        "conservatively: one historical n-hexane label for every cell, the dry-mass-weighted "
        "mean of the shells' X_c(T_i) (the kernel's sliding front needs identical adjacent "
        "labels; the whole-particle n-hexane exact), each cell's temperature from its "
        "overlapped shell pieces' summed native wet-cell energy at that label",
        "unchanged": "the interface row, its film, Ackermann and source laws, the 1e-10 residual "
        "and ledger gates, bounds, the 150-evaluation allowance, the depletion event, the birth "
        "and every later stage, every law and declared value",
        "stamped": "feed.json kernel_envelope_feed_cells, rows shell_temperatures_k; "
        "case.json kernel_envelope",
    }


# ---- The film-covered surface law (adopted 2026-10-25) and the surface-rules study knob ----------


def surface_contact_conductance():
    """The contact conductance K_l of the adopted surface law, read at call time from the one
    declared constant K_L_MOL_M2_S (so that the owner's choice of the value drops in there alone).
    """
    value = K_L_MOL_M2_S
    if isinstance(value, bool) or not isinstance(value, float):
        raise TypeError("K_L_MOL_M2_S must be a float")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("K_L_MOL_M2_S must be finite and positive")
    return value


def surface_contact_law_declared_row():
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of the model's film-covered surface law, its
    path-face rules and its one declared constant K_L_MOL_M2_S (stamped in birth.json, front.json
    and case.json). ``in_force`` states what runs: the model, or a one-line revert."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": SURFACE_LAW_RULING,
        "record": SURFACE_LAW_FINAL_RECORD,
        "adoption_record": SURFACE_LAW_ADOPTION_RECORD,
        "derivation_record": SURFACE_LAW_RECORD,
        "measurement_records": [
            SURFACE_LAW_STEFAN_DOMAIN_RECORD,
            SURFACE_LAW_STEFAN_C3_RECORD,
            CAP_FORCE_RECORD,
        ],
        "physically_qualifying": False,
        "name": "finite_contact_conductance_on_the_film_to_particle_water_uptake",
        "law_name": SURFACE_LAW,
        "law": "while an external water film covers the surface: N_exempt = max(0, N_stefan), "
        "N_stefan = N_gas - J_w; N_total = N_exempt + s R, R = N_nat,free + N_liq,free + N_gas - "
        "N_exempt, s = K_l F / (K_l F - R), i.e. N_total = N_exempt - F / (d / L_eff + 1 / K_l), "
        "L_eff / d = -R / F; where N_stefan <= 0 every path (the whole gas-frame water flux "
        "included) is multiplied by s, where N_stefan > 0 the Stefan part crosses unresisted and "
        "J_w, the native retained and the added liquid paths are multiplied by s",
        "symbols": "outward-positive fluxes; F = ln a_w,ret(surface) - ln a_w,ret(outer cell) on "
        "the bridged isotherm of the actual retained states; d the outer half cell; N_nat,free = "
        "-(L_ret / d) F the native retained path on that force; N_liq,free = -((factor - 1) "
        "L_ret / d) F; N_gas = y_w,adv N_t + J_w the gas-frame water flux, N_stefan = y_w,adv N_t "
        "its advective (Stefan) part, J_w its diffusive part; K_l the declared contact conductance",
        "continuity": "N_exempt is continuous in N_t (zero at N_t = 0), its derivative jumping "
        "there; not smoothed",
        "path_forces": PATH_FORCES,
        "path_forces_law": "K1: on every actual-path face inside a liquid-connected scope (the "
        "moving wet-front interface face, the film-free surface face, the newborn face tangent) a "
        "path segment above the Luikov cap (a_w > 0.799) carries the bridged LUIKOV_SMOOTH retained "
        "force (ideal plus second-virial quadrature), no zero-force RETAINED_CAP booking; the "
        "film-covered surface face's native evaluation keeps the native classification (its "
        "retained path is N_nat,free on F). K2: on an interior face the added liquid path's force "
        "is ln(a(W_R*) / a(W_L*)), W* each cell's equilibrium retained loading at the common face "
        "temperature (T_L + T_R) / 2 and its own composition; the donor's actual retained "
        "enthalpy, the mobility (factor - 1) L_ret and the dissipation's form unchanged",
        "K_l_mol_m2_s": surface_contact_conductance(),
        "K_l_constant": "particle_kernel_envelope.K_L_MOL_M2_S (the one point of definition)",
        "K_l_provenance": "the supplement's contact-scale anchor (sec. contactscale: Hemmingsen et "
        "al. 2008, the 1.43-mm meal's 90-s hot-water uptake, a lumped sphere filled from 0.134 to "
        "W_cap - 0.01) read with the film-covered surface at unit activity (psi_s = 0): anchor "
        "integral 0.28889, K_l 0.0492 at 80 C, 0.0619 at 100 C by inverse water viscosity, "
        "declared 0.062",
        "K_l_was": "0.32 (the same anchor printed with the surface at W_cap, integral 1.48428; "
        "the manuscript's parameter table, Conditional), until the final adoption of 2026-10-25",
        "in_force": {
            "surface_law": "former_face_reverted" if former_surface_face() else SURFACE_LAW,
            "path_forces": path_forces_in_force(),
            "one_line_reverts": {
                "A": "MODEL_REVERT_SURFACE_LAW = True",
                "B": "MODEL_REVERT_PATH_FORCES = True",
                "C": "K_L_MOL_M2_S = 0.32",
            },
        },
        "booking": "the native retained row s N_nat,free and the added liquid s N_liq,free, both "
        "at the donor actual state's retained enthalpy (dissipation R s L d^-2 F^2 for each "
        "mobility L); the gas-frame water row s N_gas with its enthalpy flux s times the native "
        "(N_stefan <= 0), or N_stefan + s J_w at the native upwind gas water partial enthalpy "
        "(N_stefan > 0); the hexane row, J_w and N_t (the Newton's unknown, pinned by the hexane "
        "rows) unchanged",
        "domain": "a trial point where K_l F - R does not have the sign of F, or F = 0, is refused "
        "(FiniteContactUndefined, a chart refusal), never clamped",
        "scope": "every surface face under a positive external water film: the birth, the front "
        "stage, the core-extinction projection and the dry stage's film-bearing active sets; K1 "
        "and K2 in every liquid-connected scope (src's classifier is not edited)",
        "unchanged": "the wet-front trace, the film node, the film-stock supply limiter's law (it "
        "acts on the added liquid path), every gate, tolerance, budget and declared value",
        "removed": "the study knobs " + ", ".join(sorted(REMOVED_ITEMS)),
        "stamped": "front.history.jsonl and dry.history.jsonl kernel_envelope_finite_contact on "
        "every accepted step (the per-step conductance K_l_mol_m2_s, the Stefan branch, the "
        "exempt flux, the path forces in force)",
    }


def surface_law_historical(environ=None):
    """SURFACE_LAW_HISTORICAL's declared value when armed, None when unset (refused if malformed)."""
    return valued(SURFACE_LAW_HISTORICAL, environ)


def former_surface_face(environ=None):
    """True when the former film-covered surface face runs in place of the model's law: the
    one-line revert MODEL_REVERT_SURFACE_LAW, or the historical switch at former_unresisted_routes.
    Read at scope entry by experiment_core2_radial_fast_wetting.liquid_connected_fast_path."""
    if MODEL_REVERT_SURFACE_LAW is not False and MODEL_REVERT_SURFACE_LAW is not True:
        raise TypeError("MODEL_REVERT_SURFACE_LAW must be True or False")
    return (
        MODEL_REVERT_SURFACE_LAW or surface_law_historical(environ) == SURFACE_LAW_HISTORICAL_FORMER
    )


def former_path_forces(environ=None):
    """True when the former path-face rules run in place of K1 and K2: the one-line revert
    MODEL_REVERT_PATH_FORCES, or the historical switch at former_cap_rule_and_own_temperature_force.
    Read at scope entry by experiment_core2_radial_fast_wetting.liquid_connected_fast_path."""
    if MODEL_REVERT_PATH_FORCES is not False and MODEL_REVERT_PATH_FORCES is not True:
        raise TypeError("MODEL_REVERT_PATH_FORCES must be True or False")
    return (
        MODEL_REVERT_PATH_FORCES
        or surface_law_historical(environ) == SURFACE_LAW_HISTORICAL_FORMER_PATH_FORCES
    )


def path_forces_in_force(environ=None):
    """The name of the path-face rules in force (report only)."""
    if former_path_forces(environ):
        return SURFACE_LAW_HISTORICAL_FORMER_PATH_FORCES
    return PATH_FORCES


def surface_law_historical_declared_row(value):
    """The typed row of the labelled historical study switch SURFACE_LAW_HISTORICAL."""
    if value not in VALUED_ITEMS[SURFACE_LAW_HISTORICAL]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + SURFACE_LAW_HISTORICAL}={value!r} is not a declared value"
        )
    common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "historical study switch of the surface-law adoptions (2026-10-25), kept "
        "reachable only for the tests the final record lists; NOT a model option",
        "record": SURFACE_LAW_FINAL_RECORD,
        "physically_qualifying": False,
        "knob": f"{PREFIX + SURFACE_LAW_HISTORICAL}={value}",
    }
    if value == SURFACE_LAW_HISTORICAL_FORMER_PATH_FORCES:
        return {
            **common,
            "name": "historical_former_path_face_rules",
            "law": "the former path-face rules: on an actual-path face every path segment above the "
            "Luikov cap booked RETAINED_CAP with zero retained force (the native classifier, not "
            "rebound); the added liquid path's interior-face force from each cell's equilibrium "
            "retained loading at the cell's own temperature",
            "scope": "every liquid-connected scope; the film-covered surface law stays the model's",
            "stamped": "the per-step kernel_envelope_finite_contact stamps carry path_forces "
            "former_cap_rule_and_own_temperature_force",
        }
    return {
        **common,
        "name": "historical_former_film_covered_surface_face",
        "law": "the former surface face: the added liquid path alone in series with the contact, "
        "N_liq = -F / (d / ((factor - 1) L_ret) + 1 / K_l); the native retained path and the "
        "pore-gas route cross the outer half cell unresisted (their uptake grows without bound "
        "as d -> 0)",
        "scope": "every film-covered surface face the model's law would carry (birth, front "
        "stage, core-extinction projection, dry stage); the contact conductance is still "
        "K_L_MOL_M2_S; the path-face rules stay the model's",
        "stamped": "the per-step kernel_envelope_finite_contact stamps carry paths None (the "
        "model's law is not evaluated)",
    }


def liquid_mobility_rule(environ=None):
    """LIQUID_MOBILITY_RULE's declared value when armed, None when unset (refused if malformed).
    Read by the dry-stage dispatcher only (particle_postcore_solver)."""
    return valued(LIQUID_MOBILITY_RULE, environ)


def surface_rules_declared_row(item, value):
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of the surface-rules study knob."""
    if item != LIQUID_MOBILITY_RULE or value not in VALUED_ITEMS[item]:
        raise EnvelopeKnobRefused(f"{PREFIX + str(item)}={value!r} is not a declared value")
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "study knob for the owner's ruling (a Class-B physics alternative; nothing "
        "adopted, nothing recommended)",
        "record": SURFACE_RULES_RECORD,
        "physically_qualifying": False,
        "knob": f"{PREFIX + item}={value}",
        "scope": "the dry stage only (particle_postcore_solver.advance_dry_phases: positive film, "
        "depletion event, wet contact and film-free active sets); the birth, the front stage and "
        "the core-extinction projection keep the former rules",
        "unchanged": "every gate, tolerance, budget, declared value, the adopted surface law and "
        "the film-stock supply limiter's law; every face and path the row does not name",
        "name": "liquid_connected_mobility_as_a_property_of_the_local_retained_state",
        "law": "the added liquid path (factor - 1) L_ret on the actual retained activities acts on "
        "an interior face when W_left > W_cap(left) or W_right > W_cap(right), W the actual cells' "
        "retained loadings on the ruled isotherm, independent of the exterior film; for one "
        "material the cell above the cap is the donor (upwind) cell of the positive-force path",
        "state_read": {
            LIQUID_MOBILITY_LOCAL_STATE: "the trial state's loadings (an implicit switch inside "
            "the Newton solve)",
            LIQUID_MOBILITY_LOCAL_STATE_LAGGED: "the loadings of the advance's accepted start "
            "state (coupled_transport._equilibrium_states of the old state); the active set is "
            "fixed for the advance, every trial, halving and active set of it; a face is "
            "identified by exact equality of its two cells' (T, y) with consecutive entries of "
            "the trial's latest _equilibrium_states call; an unidentified face falls back to the "
            "trial-state rule and is counted",
        }[value],
        "surface": "the surface face's added path is the contact with the exterior film and stays "
        "film-gated (positive film, depletion event and wet contact active sets) under the "
        "adopted surface law, as before",
        "former": "every interior face carries the added path while the exterior film exists and "
        "none after its exact-zero depletion event",
        "stamped": "dry.history.jsonl kernel_envelope_surface_rules on every accepted advance",
    }


# ---- BUDGETS_PER_UNKNOWN: the per-unknown budget forms of Paper-1 ruling packet 6 ----------------
#
# Every value below is declared by GT_PS2_PAPER1_RULING_PACKET_6_PER_UNKNOWN_BUDGETS_2026-10-22 (items
# 1 to 3) or by the lead's brief recorded in GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md
# (items 4 and 5); each floor is the site's present literal, so no form returns less than today.
# Item 1, the front search (packet 6 item 1): forward-difference Newton steps of n + 1 residual
# evaluations each; four steps per restart cycle (a positive core target needs four from the
# driver's duration seed, extinction diagnosis of 2026-10-22), twelve in the whole search (up to
# nine search steps measured for a late core-cell partial on the 80-cell graded mesh).
# GT_PS2_PAPER1_RULING_PACKET_6_PER_UNKNOWN_BUDGETS_2026-10-22 item 1:
FRONT_SEARCH_NEWTON_STEPS_PER_CYCLE = 4
FRONT_SEARCH_NEWTON_STEPS_TOTAL = 12  # GT_PS2_PAPER1_RULING_PACKET_6_PER_UNKNOWN_BUDGETS_2026-10-22
# Item 2, the DT birth Newton (packet 6 item 2: seven to nine Newton steps measured at 32 to 112
# cells, plus margin). Corrected on 2026-10-25 (GT_PS2_PAPER1_SURFACE_LAW_ADOPTION_2026-10-25.md):
# the budgets build measured the births' assembly counter exactly in n, 10 (n + 1) + 9 for every
# nine-step birth (1,389, 1,709, 2,029, 2,349, 2,669 at 64 to 128 cells: nine Newton steps and the
# certificate Jacobian at n + 1 each, plus nine single evaluations), so the packet's 10 (n + 1)
# refused the 128-cell birth by nine assemblies. 11 (n + 1) >= 10 (n + 1) + 9 whenever n >= 8 (every
# DT birth: n = 2 N + 9), with n - 8 single evaluations to spare, in the packet's per-unknown form.
BIRTH_NEWTON_STEPS = 11  # GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md section 5 (births)
# Item 3, the front's proof-point allowance per advance (packet 6 item 3: 5,600 to 6,300 points per
# unknown measured at 64 and 80 cells, times 2.5); the floor is the front driver's present literal.
# GT_PS2_PAPER1_RULING_PACKET_6_PER_UNKNOWN_BUDGETS_2026-10-22 item 3:
FRONT_PROOF_POINTS_PER_UNKNOWN = 15_000
# Item 4, the half-step-first threshold as a fraction of the advance's allowance in force. The
# threshold was declared (evidence 168/169) to sit above every accepted film-bearing full step AND to
# leave the halving at least its measured cost; the build's 6 / 7 kept the ratio 1,200 / 1,400 but
# left A / 7 for the halving (416 of 2,910 at 80 cells) against the 970 to 1,696 its own census
# measured. Corrected on 2026-10-25 (GT_PS2_PAPER1_SURFACE_LAW_ADOPTION_2026-10-25.md) on that
# build's item-4 arithmetic: the allowance in force is 3 (4 u + 6) = 12 (u + 1) + 6; the accepted
# film-bearing full steps cost 5 (u + 1) + 2 (1,212 at 80 cells, 2,412 at 160) and the largest
# half step after an abandonment 7 (u + 1) + 2 at 80 cells (1,696; 2,916 at 160), so the fraction f
# must satisfy 1,212 <= f 2,910 and (1 - f) 2,910 >= 1,696, i.e. 0.41649 <= f <= 0.41718: 5 / 12
# (floor(5 A / 12) = 1,212 and 2,412 admit those full steps; the halving keeps 1,698 and 3,378).
HALF_STEP_FIRST_FRACTION_NUMERATOR = 5  # GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md s5
HALF_STEP_FIRST_FRACTION_DENOMINATOR = 12  # records/m03_halfstep_census.json of that build
# Item 5, the dry-stage advance ceiling: the planned requests times a margin of two (every planned
# request accepted as two half steps; the complete DT histories used 1.00 to 1.27 times their plan),
# the runner's present 1,500 the floor.
DRY_ADVANCE_CEILING_MARGIN = 2  # GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md
# The planned-request walk refuses past this many requests (a walk guard, not a budget).
DRY_PLANNED_REQUEST_WALK_LIMIT = 10_000_000  # GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md


def budgets_per_unknown(environ=None):
    """BUDGETS_PER_UNKNOWN's declared value when armed, None when unset.

    Read only through a live envelope activation: the whole environment is validated first, so a
    malformed or unknown DTDC_KERNEL_ENVELOPE_* variable refuses the budgets with every other knob.
    """
    validate_environment(environ)
    return valued(BUDGETS_PER_UNKNOWN, environ)


def _positive_integer(value, what):
    n = int(value)
    if n != value or n < 1:
        raise ValueError(f"a positive integer {what} is required")
    return n


def front_search_total_budget(present, unknowns):
    """Item 1: the front search's total residual evaluations, max(present, 12 (n + 1))."""
    n = _positive_integer(unknowns, "unknown count")
    return max(
        _positive_integer(present, "present budget"), FRONT_SEARCH_NEWTON_STEPS_TOTAL * (n + 1)
    )


def front_search_cycle_budget(present, unknowns):
    """Item 1: each of the first two restart cycles, max(present, 4 (n + 1) + 1); the caller still
    caps it by what the total has left."""
    n = _positive_integer(unknowns, "unknown count")
    return max(
        _positive_integer(present, "present budget"),
        FRONT_SEARCH_NEWTON_STEPS_PER_CYCLE * (n + 1) + 1,
    )


def birth_residual_budget(present, unknowns):
    """Item 2: the DT birth Newton's residual-assembly budget, max(present, 10 (n + 1))."""
    n = _positive_integer(unknowns, "unknown count")
    return max(_positive_integer(present, "present budget"), BIRTH_NEWTON_STEPS * (n + 1))


def front_proof_point_allowance(present, unknowns):
    """Item 3: the front's proof-point allowance per advance, max(present, 15,000 (n + 1))."""
    n = _positive_integer(unknowns, "unknown count")
    return max(
        _positive_integer(present, "present budget"), FRONT_PROOF_POINTS_PER_UNKNOWN * (n + 1)
    )


def front_search_unknowns(transport):
    """n of the moving-front water-node solve on a front-stage state: wet temperature and water per
    wet piece, dry temperature and composition per dry piece, one Stefan flux per dry face, the front
    (or the target's duration), the interface temperature and the surface's three: 2 n_w + 2 n_d +
    n_f + 5, the Newton's starting-vector size on the stored G64 and N64 core-cell states (245, 197).
    """
    layout = transport.layout
    return 2 * layout.wet_piece_count + 2 * layout.dry_piece_count + layout.dry_face_count + 5


def dry_half_step_threshold(allowance):
    """Item 4: the half-step-first threshold, max(1,200, floor(5 A / 12)), A the advance's
    allowance in force. 1,200 wherever A is at most 2,881 (every allowance to 79 cells); where A is
    below 1,200 it is above the allowance and never fires (as today)."""
    a = _positive_integer(allowance, "assembly allowance")
    return max(
        DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD,
        (HALF_STEP_FIRST_FRACTION_NUMERATOR * a) // HALF_STEP_FIRST_FRACTION_DENOMINATOR,
    )


def dry_planned_requests(time_control, start_clock_s, horizon_s):
    """The requests the declared time controls plan from a clock to the horizon when every request
    is accepted whole (time_control.dry_request walked from its own landings)."""
    clock, count = float(start_clock_s), 0
    while clock < horizon_s:
        # The plan has no film loading: film_present=False, the former bands (ruling packet 7).
        clock = clock + time_control.dry_request(clock, horizon_s, film_present=False)
        count += 1
        if count > DRY_PLANNED_REQUEST_WALK_LIMIT:
            raise RuntimeError("planned dry-request walk did not reach the horizon")
    return count


def dry_advance_ceiling(time_control, start_clock_s, horizon_s, present):
    """Item 5: the dry stage's declared maximum physical advances, max(present, 2 P) (the runner's
    present 1,500), P the planned requests from the dry stage's first clock; returns the ceiling and
    its stamp."""
    floor = _positive_integer(present, "present budget")
    planned = dry_planned_requests(time_control, start_clock_s, horizon_s)
    ceiling = max(floor, DRY_ADVANCE_CEILING_MARGIN * planned)
    return ceiling, {
        "start_clock_s": start_clock_s,
        "horizon_s": horizon_s,
        "schedule": time_control.dry_schedule,
        "time_level": time_control.level,
        "planned_requests": planned,
        "margin": DRY_ADVANCE_CEILING_MARGIN,
        "floor": floor,
        "ceiling": ceiling,
    }


def budgets_per_unknown_declared_row(value=BUDGETS_PER_UNKNOWN_PACKET6):
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of BUDGETS_PER_UNKNOWN at its declared value."""
    if value not in VALUED_ITEMS[BUDGETS_PER_UNKNOWN]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + BUDGETS_PER_UNKNOWN}={value!r} is not a declared value"
        )
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": f"{BUDGETS_PER_UNKNOWN_PACKET} (Class B: declared solver budgets of the "
        "Paper-1 harness; built for the owner's confirmation, nothing adopted)",
        "record": BUDGETS_PER_UNKNOWN_RECORD,
        "evidence": BUDGETS_PER_UNKNOWN_EVIDENCE,
        "physically_qualifying": False,
        "name": "paper1_harness_budgets_per_unknown",
        "knob": f"{PREFIX + BUDGETS_PER_UNKNOWN}={value}",
        "n": "the solve's search unknowns (its Newton's starting-vector size)",
        "front_search": "the water-node restart loop (core2_water_film_node.advance_water_node) of "
        "every front partial, positive target and film event called through "
        "experiment_core2_finite_wetting_b.advance_partial: each of the first two cycles max(800, "
        "4 (n + 1) + 1) residual evaluations within the total max(2,400, 12 (n + 1)); the front "
        "driver's arrival solve (experiment_core2_dt_face.advance_arrival) max(2,400, 12 (n + 1))",
        "birth": "the DT birth Newton (run_particle_qualification_case.fresh_dt_birth into "
        "experiment_core2_dt_birth.advance_birth): max(2,400, 11 (n + 1)) (corrected 2026-10-25 "
        "from 10 (n + 1): every nine-step birth measured 10 (n + 1) + 9 assemblies)",
        "front_proof_points": "per front advance (qualify_sliding_front_history): max(3,000,000, "
        "15,000 (n + 1)), n = 2 n_w + 2 n_d + n_f + 5 of the advance's start state",
        "dry_half_step_threshold": "max(1,200, floor(5 A / 12)), A the advance's assembly allowance "
        "in force (corrected 2026-10-25 from floor(6 A / 7), which left A / 7 for the halving): "
        "the measured full step 5 (u + 1) + 2 admitted and the measured largest half step "
        "7 (u + 1) + 2 left for the halving within A = 12 (u + 1) + 6",
        "corrections": SURFACE_LAW_ADOPTION_RECORD,
        "dry_advance_ceiling": "max(1,500, 2 P), P the requests the declared time controls plan from "
        "the dry stage's first clock to the horizon; the same count divides the per-advance "
        "component budgets (particle_original_feed_balance.component_error_budgets_per_advance)",
        "floors": "every present value: 800 and 2,400 (front search), 2,400 (birth), 3,000,000 "
        "(proof points), 1,200 (half-step threshold), 1,500 (advance ceiling)",
        "unchanged": "every tolerance, gate, law and seed, the 1e-5 seed-reuse threshold, the 1e-11 "
        "search target, the 1e-10 ledger gate, the 1e12 condition cap, the dry assembly allowance "
        "and its knob, the 13 halvings, the four-target endpoint rule, the minimum interval, the "
        "per-path 1,200 proof limit, the non-physical solver path, the Faner birth and every other "
        "budget",
        "stamped": "front.json kernel_envelope_budgets_per_unknown; every front history row's proof "
        "point_budget; birth.json audit kernel_envelope_budgets_per_unknown; case.json "
        "kernel_envelope_budgets_per_unknown (the dry ceiling) and history_time_control"
        ".dry_physical_advance_limit; the half-step stamps' declared_numerics_threshold",
    }


def core_target_seed_points(checkpoint, time_s, front_z, resume=None):
    """The accepted moving front points CORE_TARGET_SEED reads at a (re)start, and their source.

    Points are (clock, front volume fraction), at most the last three, the last the start state.
    A checkpoint written with the knob armed carries them. Without them, a resume reads the
    parent's accepted history beside its checkpoint (<name>.history.jsonl for
    <name>.accepted.pkl): every row but HELD_FRONT_ROW_KINDS, one point per distinct clock, used
    only when the last equals the start state exactly. Otherwise the start state alone (the
    driver's seed until two secants are accepted). Seeding only.
    """
    state = (time_s, front_z)
    stored = checkpoint.get(CORE_TARGET_SEED_POINTS_KEY)
    if stored:
        points = [tuple(point) for point in stored][-3:]
        if points[-1] == state:
            return points, {"source": "checkpoint", "points": points}
        return [state], {
            "source": "start_state",
            "reason": "the checkpoint's seed points do not end at its state",
            "rejected_points": points,
        }
    history = None
    if resume is not None:
        resume = Path(resume)
        history = resume.with_name(resume.name.replace(".accepted.pkl", ".history.jsonl"))
    if history is None or history == resume or not history.is_file():
        return [state], {
            "source": "start_state",
            "reason": "no accepted front history beside the parent checkpoint",
        }
    points = []
    with history.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") in HELD_FRONT_ROW_KINDS or "front_z" not in row:
                continue
            point = (row["time_s"], row["front_z"])
            if points and points[-1][0] == point[0]:
                if points[-1] != point:
                    return [state], {
                        "source": "start_state",
                        "reason": "the parent history repeats a clock with another front",
                        "history": str(history),
                    }
                continue
            points.append(point)
    points = points[-3:]
    if not points or points[-1] != state:
        return [state], {
            "source": "start_state",
            "reason": "the parent history does not end at the checkpoint's state",
            "history": str(history),
            "history_last_point": points[-1] if points else None,
        }
    return points, {"source": "parent_history", "history": str(history), "points": points}


def core_target_duration_seed(rule, points, time_s, front_z, driver_seed_s):
    """CORE_TARGET_SEED: the positive core target's duration seed and its stamp.

    history_power_law: dz/dt = -c z^p through the secants of the last two accepted moving steps
    (points: the last three (clock, z), the last the current state), each secant placed at the
    geometric mean of its end points; the seed is the halving time from the current z. Fewer than
    two positive accepted secants, or a fitted seed that is not finite and positive: the driver's
    seed. No clamp beyond positivity.
    """
    if rule not in VALUED_ITEMS[CORE_TARGET_SEED]:
        raise EnvelopeKnobRefused(f"{PREFIX + CORE_TARGET_SEED}={rule!r} is not a declared value")
    stamp = {
        "rule": rule,
        "time_s": time_s,
        "front_z": front_z,
        "target_z": 0.5 * front_z,
        "driver_seed_s": driver_seed_s,
    }

    def driver_seed(reason):
        return driver_seed_s, {**stamp, "seed_s": driver_seed_s, "seed": "driver", "reason": reason}

    recent = [tuple(point) for point in points][-3:]
    if len(recent) < 3:
        return driver_seed("fewer than two accepted secants")
    if recent[-1] != (time_s, front_z):
        return driver_seed("the seed points do not end at the current state")
    (t0, z0), (t1, z1), (t2, z2) = recent
    if not (t0 < t1 < t2 and z0 > z1 > z2 > 0.0):
        return driver_seed("the last two accepted secants are not both positive")
    try:
        r1 = (z0 - z1) / (t1 - t0)
        r2 = (z1 - z2) / (t2 - t1)
        m1 = math.sqrt(z0 * z1)
        m2 = math.sqrt(z1 * z2)
        p = math.log(r1 / r2) / math.log(m1 / m2)
        c = r2 / m2**p
        if 1 - p == 0:
            seed = math.log(2.0) / c
        else:
            seed = (z2 ** (1 - p) - (z2 / 2) ** (1 - p)) / (c * (1 - p))
    except ArithmeticError:  # an overflow or a zero divisor
        return driver_seed("the fitted seed is not representable")
    except ValueError:  # a logarithm outside its domain
        return driver_seed("the fitted seed is not representable")
    if not (math.isfinite(seed) and seed > 0.0):
        return driver_seed("the fitted seed is not finite and positive")
    return seed, {
        **stamp,
        "seed_s": seed,
        "seed": rule,
        "points": recent,
        "secant_rates": [r1, r2],
        "power_p": p,
        "coefficient_c": c,
    }


# ---- INSTRUMENT_SHARED_LAW: the report-only water-activity instrument -------------------


def ideal_raoult_dalton_water_activity(temperature_k, pressure_pa, y_hexane):
    """The ideal immiscible Raoult-Dalton reading (evidence 163's L0): (1 - y_h) P / P_sat,w(T)
    with the frozen IAPWS-95 saturation pressure. An instrument reading, never a model law."""
    from dtdc_simulator.core2.props import water

    return (
        (1.0 - float(y_hexane))
        * float(pressure_pa)
        / water.saturation_pressure(float(temperature_k))
    )


def shared_law_water_activity(temperature_k, pressure_pa, y_hexane, k_wh=None):
    """The water activity on the shared equilibrium (evidence 163's L1) through the kernel's own
    chart activity function coupled_pore._gas_component_activities (PHY-053 second-virial
    binary gas at the kernel's k_wh over the pure compressed liquid's fugacity, immiscible pure
    liquids). The same function the gas-only chart and evaluate_equilibrium read, but without
    the dry-pore chart's refusal above 1, so a supersaturated gas reads its activity."""
    from dataclasses import replace

    from dtdc_simulator.core2.particle import coupled_pore as cp

    params = cp.CoupledPoreParams()
    if k_wh is not None:
        params = replace(params, k_wh=float(k_wh))
    water_activity, _ = cp._gas_component_activities(
        float(temperature_k), float(pressure_pa), float(y_hexane), params
    )
    return float(water_activity)


def instrument_water_activity(temperature_k, pressure_pa, y_hexane, environ=None, k_wh=None):
    """A report-only water-activity reading. Knob off: the former ideal reading alone. Armed:
    the shared law's reading, with the ideal one stamped beside it and their difference."""
    ideal = ideal_raoult_dalton_water_activity(temperature_k, pressure_pa, y_hexane)
    if not armed(INSTRUMENT_SHARED_LAW, environ):
        return {"water_activity": ideal, "reading": "ideal_raoult_dalton"}
    shared = shared_law_water_activity(temperature_k, pressure_pa, y_hexane, k_wh)
    return {
        "water_activity": shared,
        "reading": "shared_law",
        "water_activity_shared_law": shared,
        "water_activity_ideal_raoult_dalton": ideal,
        "shared_minus_ideal": shared - ideal,
        "report_only": True,
        "record": QUEUE_RECORD,
        "physically_qualifying": False,
    }


def bath_instrument(bath, environ=None):
    """INSTRUMENT_SHARED_LAW armed: the bath gas's water activity at every waypoint, at the
    waypoint's own pressure (101,325 Pa for a bath without pressure fields), on both readings.
    The shadow harness's bath ``water_activity`` instrument (evidence 161 bath_topology.py,
    evidence 163 section 1.2) brought into the harness. None when the knob is off or there is
    no bath. Report only: nothing reads it back."""
    if bath is None or not armed(INSTRUMENT_SHARED_LAW, environ):
        return None
    rows = []
    for waypoint in bath["waypoints"]:
        t = float(waypoint["t_s"])
        pressure = bath_pressure_pa(bath, t)
        temperature = float(waypoint["T_bath_K"])
        y = float(waypoint["y_hexane"])
        reading = instrument_water_activity(temperature, pressure, y, environ)
        rows.append(
            {
                "t_s": t,
                "temperature_k": temperature,
                "pressure_pa": pressure,
                "y_hexane": y,
                "water_activity_shared_law": reading["water_activity_shared_law"],
                "water_activity_ideal_raoult_dalton": reading["water_activity_ideal_raoult_dalton"],
            }
        )
    shared_over = [r for r in rows if r["water_activity_shared_law"] > 1.0]
    disagree = [
        r
        for r in rows
        if (r["water_activity_shared_law"] > 1.0) != (r["water_activity_ideal_raoult_dalton"] > 1.0)
    ]
    return {
        "instrument": "bath gas water activity",
        "reading": "shared_law",
        "report_only": True,
        "waypoints": rows,
        "waypoints_water_supersaturated_on_the_shared_law": len(shared_over),
        "waypoints_where_the_readings_disagree_on_saturation": len(disagree),
        "record": QUEUE_RECORD,
        "physically_qualifying": False,
    }


# ---- DEFERRED_CERTIFICATE: the dense certificate taken later (2026-11-03) ------------------------


def deferred_certificate(environ=None):
    """DEFERRED_CERTIFICATE's value ("1") when armed; None when unset; refused otherwise."""
    return valued(DEFERRED_CERTIFICATE, environ)


def deferred_certificate_declared_row():
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of DEFERRED_CERTIFICATE."""
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "the lead's Class-A numerics in the study lane (the delegation charter of "
        "2026-08-11), on the external parallelization study's hint A of 2026-11-03; a measurement, "
        "nothing adopted; no gate, tolerance, budget, cap, allowance or law moved",
        "record": DEFERRED_CERTIFICATE_RECORD,
        "physically_qualifying": False,
        "name": "dense_certificate_taken_after_the_next_search",
        "knob": f"{PREFIX + DEFERRED_CERTIFICATE}=1",
        "where": "the front driver's run (qualify_sliding_front_history): every water-node solve "
        "(core2_water_film_node, the physical chart solve) whose iterate meets gtol; not the birth, "
        "the dry stage, the endpoint surface solves or an acceptance clause's candidate",
        "certificate": "the same dense finite-difference Jacobian at the same iterate by the same "
        "statements, counted outside max_nfev as before, its proof points on the step's own "
        "allowance, inside the residual's process-global hooks as they were at the return",
        "consumers": "the film node's condition gate and its stamps, and the condition number the "
        "accepted step's ledger stores (a sentinel until the certificate is in, then filled where "
        "it was copied); rows and saves wait for it",
        "refusal": "inside its own step, the former error; after the driver moved on, the run "
        "replayed from its start with that certificate taken at once (stamped)",
        "former": "every dense certificate taken where its iterate is returned",
    }


# ---- JACOBIAN_COLUMN_GROUPS: finite-difference Jacobians by column groups (2026-11-02) ----------


def jacobian_column_groups_declared_row(value):
    """The typed DECLARED_ENGINEERING_ASSUMPTION row of JACOBIAN_COLUMN_GROUPS."""
    if value not in VALUED_ITEMS[JACOBIAN_COLUMN_GROUPS]:
        raise EnvelopeKnobRefused(
            f"{PREFIX + JACOBIAN_COLUMN_GROUPS}={value!r} is not a declared value"
        )
    row = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": "the lead's Class-A numerics in the study lane (the delegation charter of "
        "2026-08-11), on the owner's direction of 2026-11-02 to speed up the algorithms; no gate, "
        "tolerance, budget, cap, allowance or law moved",
        "record": JACOBIAN_COLUMN_GROUPS_RECORD,
        "physically_qualifying": False,
        "name": "finite_difference_jacobian_by_declared_column_groups",
        "knob": f"{PREFIX + JACOBIAN_COLUMN_GROUPS}={value}",
        "where": "the domain-aware Newton's finite-difference Jacobians at the declared sites: the "
        "birth (experiment_core2_dt_birth), the front stage (core2_water_film_node's physical "
        "chart solve) and the dry stage's film-bearing and film-free solves "
        "(experiment_core2_dt_dry.advance); every other solve one column at a time",
        "structure": "declared supersets: cut_sparsity.build_cut_jacobian_sparsity (both "
        "widenings) and cut_birth_integrator.birth_sparsity with the surface node dense; the dry "
        "cell stencil with the conservative elimination's dense rows "
        "(particle_jacobian_column_groups)",
        "exact": "each column the same step, sign, subtraction and division as the one-at-a-time "
        "loop; a row outside the declared structure the same signed zero; the last probe the "
        "former last probe; a probe moving a row outside the structure falls back to the "
        "one-at-a-time Jacobian (stamped, the layout disabled); the first "
        "particle_jacobian_column_groups.BIT_CHECK_SOLVES solves of a process bit-checked",
        "former": "every finite-difference Jacobian one residual evaluation per unknown",
        "stamped": "kernel_envelope_jacobian_column_groups (per solve: site, unknowns, groups, "
        "grouped and fallback Jacobians, group and solo evaluations, violations, bit checks; "
        "dry-rowsum: the assembled rows, the far groups and the fit residual)",
        "unchanged": "the convergence test, the acceptance clause, max_nfev and its counting rule, "
        "the line search, the damping ladder, the bounds handling, every gate, tolerance, cap and "
        "allowance",
    }
    if value == "dry-rowsum":
        row["dry_rowsum"] = (
            "the dry stage's search Jacobian (never a certificate): the outer cell's three rows and "
            "the surface rows, whose dependence on the far cells' T and q passes through the "
            "residual's two whole-particle inventory sums (water and n-hexane), assembled as "
            "a_w dW_c/dx + a_h dH_c/dx, the per-cell terms from the grouped probe and the two "
            "coefficients per row fitted by least squares over the far groups; chain 257's Broyden "
            "class, the iterate path not bit-identical; the dense certificate at every returned or "
            "provisional iterate the exact Jacobian"
        )
    return row
