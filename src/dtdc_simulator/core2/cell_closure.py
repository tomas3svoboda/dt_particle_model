"""Cell-layer implicit closure (S0) -- the per-cell DAE D7 section 2.2 rules.

Implements ``docs/GT_PS2_CELL_LAYER_CLOSURE_DESIGN_BASIS_2026-08-12.md`` under
``docs/GT_PS2_D7_ARCHITECTURE_RECONCILIATION_2026-08-10.md:66-73``: "Within each
cell and step, the fast physics is solved implicitly with analytic Jacobians:
quasi-steady or dynamic gas, the PHY-022 film / interface closure, wall-node
thermal state, per-layer pressure. This is the legitimate and retained home of
the residual/jacobian/solver work, rescoped from a global system to per-cell
systems."

EXACTLY FOUR SUBSYSTEMS, NO FIFTH, AND NO GLOBAL UNKNOWN. There is no extensive
state vector spanning cells, no globally assembled residual or Jacobian, and no
solid state smeared into continuum coordinates: the superseded single-DAE
formulation does not return here even as a convenience. The closure is a pure
function of (packet-derived aggregates, cell field state, boundary conditions,
macro step) returning sources, fluxes and one wall coupling set.

THE WALL IS THE SLOWEST MODE, SO IT IS NOT IN THE FAST BLOCK. A naive reading of
"solve the fast physics implicitly" would absorb all four subsystems into one
Newton block. FROZEN PHY-008 forbids it: the metal time constant EXCEEDS the
per-tray solid residence by 1.4-8.3x (``release/physics_decisions.yaml:953``), so
the wall is slower than the solids march.

THAT SEPARATION IS DERIVED HERE, NOT DECLARED. ``wall_time_constant_s`` is
computed from the declared capacity and the declared conductance sum, the ratio
against the solids-paced macro step is reported as
``wall_time_constant_to_macro_step_ratio``, and ``WallNodeSolve.is_slowest_mode``
is a PROPERTY that tests that ratio against
``PHY008_WALL_RESIDENCE_RATIO_BRACKET``. A parameterisation whose metal is not
slower than the solids march therefore reads False and is surfaced on the solved
row, rather than being asserted by a literal that no parameterisation could
falsify. The resulting structure:

* gas (PHY-004): ALGEBRAIC by default, DIFFERENTIAL under the dynamic-EOS
  switch -- three rows in the fast Newton block;
* per-layer pressure (PHY-017): ALGEBRAIC, index-1, SIMULTANEOUS with the gas --
  one row in the fast Newton block;
* film / interface (PHY-022 coupled two-component flux with the PHY-007
  four-branch active set): three rows in the fast Newton block;
* wall metal internal energy (PHY-008): ONE DIFFERENTIAL ROW with four coupling
  terms, advanced by ``advance_wall_node`` at the macro-step scale and never
  absorbed into the fast block.

``solve_cell_fast_block`` therefore takes the wall temperature as a FROZEN
parameter and ``advance_wall_node`` takes the converged gas and interface
temperatures as frozen parameters; ``march_cell_macro_step`` composes them in
that order. The gas and pressure rows being algebraic is what keeps the
differential block free of a fast time constant, i.e. index-1 by construction.

PRESSURE IS SIMULTANEOUS, NEVER LAGGED. Pressure and evaporation couple through
2-5 K of boiling-point shift (``release/physics_decisions.yaml:1217``), so the
layer pressure is an unknown of the same implicit block as the interface state. A
lagged pressure would break the feedback the frozen decision exists to
represent: more evaporation -> more throughput -> higher local P -> higher T_sat
-> less superheat. The wall's one-way lag is a different thing entirely and is
justified by its own time constant; the pressure carries no lag at all.

TWO LENGTH SCALES, SEPARATELY NAMED, AND THE PARTICLE SCALE NEVER REACHES THE
HYDRAULIC RESIDUAL. FROZEN PHY-017's consistency requirement
(``release/physics_decisions.yaml:1223-1224``) requires the Ergun hydraulic
diameter to be the true external granule dimension and NEVER the intraparticle
diffusion radius, which would give a spurious ~0.4-0.8 bar drop. Said exactly, so
that no enforcement is claimed that the code does not perform:

* the pressure residual reads ``LayerHydraulics`` (its series and its network
  mode), ``GasBoundary.downstream_boundary_pressure_pa`` and the
  ``loss_dispatch`` seam -- and no field of ``LayerHydraulics`` or of any element
  names a particle dimension;
* the interfacial area is taken WHOLE from the inherited tray-type dispatch, so
  no particle dimension is re-derived here at all;
* the tripwire that refuses a hydraulic diameter at or below the intraparticle
  scale lives on ``BedErgunElement.__post_init__`` -- on the element that carries
  the diameter, NOT on ``LayerHydraulics``, which has no diameter field;
* the caller-supplied ``loss_dispatch`` is a seam that sits OUTSIDE that
  tripwire. A caller handing in its own composer can evaluate any arithmetic it
  likes at any length scale; what still stands in its way by default is
  ``require_corroborated_layer_drop``, which refuses the ~0.4-0.8 bar drop the
  frozen consistency requirement names. That is a second line of defence, not the
  tripwire, and it is opt-out per ``require_corroborated_pressure_drop``.

NON-SMOOTH, HONESTLY -- SWITCHED, EXACT, AND FAIL-CLOSED.

* PHY-007's four-branch interface active set (``:902-916``) is selected by EXACT
  comparisons on two solid-side inventories the packet-carrier layer already
  owns. The cell layer READS them and never re-derives them. Branch changes take
  exact one-sided limits: there is no smoothing, no mollifier and no epsilon
  shell -- a hexane loading of 1e-300 is the wet branch, not the dry one.
  Inactive components are pinned by an exact IDENTITY ROW, which is an active-set
  row and not a penalty, a clamp or a variable bound.
* Co-located liquids take ONE area scalar, computed once from the total external
  liquid saturation and multiplied once into both components and into the heat
  path (``:916`` no-duplicate-area rule, PHY-032's shared ``a_Lg``). There is no
  per-component area field anywhere in this module.
* PHY-023's area derivative GENUINELY DIVERGES: d a_hg / d S = -(2/3) a_b
  (1-S)^(-1/3) as S -> 1 (``:1495-1510``). The analytic Jacobian carries that
  divergence rather than a finite surrogate, and refuses -- typed, at the
  singular limit -- rather than emitting a silent inf or nan.
* Saturation admissibility is FAIL-CLOSED. A trial state outside 0 <= S_L < 1 is
  a typed rejection that propagates out of the implicit solve as a step
  rejection. It is NEVER a projection back into the feasible set, and no clip,
  clamp or projection exists here as a safeguard, as a line-search fallback or as
  a variable bound. The solver has no safeguard at all: non-convergence is a
  typed refusal.

THE PRE-INTEGRATION ADMISSIBILITY GATE IS OPAQUE TO THIS MODULE, AND IT IS NOT
OPTIONAL. The ruled refusal behaviour is that an out-of-domain excursion is
refused at seed assembly and never integrated. ``solve_cell_fast_block``,
``march_cell_macro_step`` and ``relax_to_zero_accumulation_fixed_point`` therefore
take ``seed_gate`` as a REQUIRED keyword: either a callable, called ONCE before
any residual evaluation, or the explicit token ``SEED_GATE_DECLARED_ABSENT``. An
ungated solve is a named choice, never a forgotten default, and which of the two
was chosen is recorded as an artefact fact on ``CellSolve.seed_gate_witness`` and
echoed into the boundary ledger -- so a gated and an ungated converged solve are
not equal objects. A gate refusal is TERMINAL -- never retried and never softened.
This module never learns what the gate checks.
The per-discretization declared assumption the gate comes from lives in
``face_regime_declaration``, which is its single home; the residual rows, the
Jacobian sparsity, the solver and the coupling variables here encode nothing
about it, which is what lets a future Option-3 closure touch that declaration and
the particle internals without touching this structure.

OPTION-3 TOLERANT BY CONSTRUCTION. The unknown count is ``len(UNKNOWN_NAMES)``
and the row count is ``len(FAST_BLOCK_ROWS)``; nothing is dimensioned by a
literal, the differentiation carries one partial per declared unknown, and the
linear solve is dimension-generic. Adding the ruled successor's single pore-wall
condensed-phase state variable is one name, one row builder and no structural
change.

THE PRESSURE PATH IS AN ORDERED SERIES OF INDEPENDENTLY TYPED LOSS ELEMENTS
(owner ruling, 2026-08-12): "The gas pressure path shall support an ordered
series of independently typed loss elements. For the Schumacher family-2
archetype, use bed Ergun resistance in series with a tapered-bore floor-passage
relation. The printed taper ratio constrains geometry only; the floor
discharge/loss coefficient remains a declared engineering assumption with a
predeclared sensitivity range. Bed and floor coefficients shall not be jointly
inferred from total pressure drop alone."

So the family-2 series is RULED, not provisional: ``LayerLossSeries`` carries an
ordered tuple of typed elements, the inherited ``PressureLossRelation`` tags
element KINDS, and ``RULED_ELEMENT_KIND_SERIES`` maps each family to its ordered
kinds. The inherited single selector per family is read as the family's PRIMARY
element kind, and that reading is asserted at import against
``CANONICAL_LOSS_RELATION``; the inherited typed-floor field itself is neither
changed nor reinterpreted. Each element carries its OWN analytic derivatives, the
series sums them, and the composition still reaches the unknowns through exactly
two derived channels (superficial velocity and gas density) with the chain rule
applied ONCE outside the elements. Adding an element therefore moves numbers and
moves neither the residual's shape nor the Jacobian's sparsity.

THE PRINTED TAPER RATIO CONSTRAINS GEOMETRY ONLY. The floor element's loss
arithmetic uses the THROAT: the declared free-area fraction is the throat area
fraction, and the printed area ratio relates the wide end to it without entering
the drop. The floor discharge coefficient is a separate
DECLARED_ENGINEERING_ASSUMPTION carrying a PREDECLARED, SWEEPABLE sensitivity
range (``FLOOR_DISCHARGE_COEFFICIENT_SENSITIVITY_RANGE``), and no nominal value is
minted here: the caller declares one and it is checked against the range.

IDENTIFIABILITY PROHIBITION, ROUTED THROUGH ONE TYPED ENTRY POINT. Two
coefficients cannot be fitted from one aggregate, so the bed and floor
coefficients may not be jointly inferred from a total pressure drop alone. This
module offers no inference, calibration or fitting path at all, and the
prohibition is not left as a bare refusing stub that a differently-named caller
could walk past: ``fit_series_coefficients_to_observation`` IS the entry point a
caller reaching for calibration arrives at, it takes a TYPED observation, and the
observation's own shape decides the refusal. An observation carrying only a total
(``LayerDropObservation.per_element_drops_pa is None``) is refused with
``JointIdentifiabilityRefusedError`` -- the ruled refusal, selected by the target's
type and not by a name. An observation carrying per-element drops is refused with
``CoefficientInferenceNotImplementedError``, because no fitting path exists here
at all; the admitted operation is ``sweep_floor_discharge_coefficient`` over the
predeclared range. ``infer_series_coefficients_from_total_drop`` is retained as
the entry point that spells the forbidden operation out in the ruling's own words
and refuses unconditionally. Every series evaluation returns its PER-ELEMENT
contributions, so an aggregate total is never the only observable.

LIMITATIONS.

* TWO ELEMENT TYPES SHIP, AND THEY ARE BOTH FAMILY 2's. ``BedErgunElement`` and
  ``TaperedBoreFloorPassageElement`` cover the PACKED_BED_ERGUN and nozzle kinds
  OF THE FAMILY-2 ARCHETYPE, and each declares the family it is admitted on
  (``admitted_families``), checked by ``LayerLossSeries``. DOME_BYPASS and
  AIR_PLENUM have NO element type here, so constructing a series for family 1 or
  family 4 refuses: the bypass path has no printed per-layer drop inside the
  corroborated window, and family 4 has no reference-machine member and no
  admitted anchor. Refusing is deliberate; a zero would silently assert a
  negligible drop, which is a nested LIMIT of the frozen network and not its
  default.
* FAMILY 3 HAS NO ELEMENT TYPE OF ITS OWN, SO ITS SERIES REFUSES TOO. Family 3's
  ruled series is a single nozzle kind, and the only nozzle element written here
  is family 2's tapered bore, whose constructor enforces the SCHUMACHER printed
  taper band. Building family 3's series from it would import a family-2 anchor's
  printed geometry onto a Desmet deck -- the cross-provenance blending the
  per-provenance admissibility gate exists to prevent -- so the construction is
  refused instead. A dedicated family-3 element (Desmet's printed sparge holes
  and slots carry no printed taper ratio and no printed open-area band) is named
  as owed work, not written here. Families 3 and 4 remain declared UNANCHORED BY
  THE REFERENCE MACHINE.
* THE NOZZLE KIND TAG IS THE NEAREST INHERITED MEMBER. The ruled tapered-bore
  floor passage is a designed nozzle, and the inherited enum's nozzle member is
  named for its family-3 role. That member is used as the KIND tag here; minting a
  dedicated inherited member would edit the inherited typed-floor module, which
  this packet does not do. The inherited family-to-relation equality it enforces
  there constrains a typed floor's SELECTOR FIELD and is untouched.
* NO NEW PRESSURE BAND, AND THE ENFORCED GATE IS THE HULL RATHER THAN THE
  INTERSECTION. The corroboration window is the union of two printed bands that
  OVERLAP on 2.0-5.0 kPa, so its hull has printed support throughout; nothing new
  is minted. Read strictly, the design basis's "must land inside the
  family-independent corroboration window -- frozen Ergun 2.0-5.6 kPa,
  independently corroborated 1-5 kPa" is the INTERSECTION [2000, 5000] Pa, which
  is narrower than the enforced hull [1000, 5600] Pa: the hull admits
  1000-2000 Pa and 5000-5600 Pa, each with single-source printed support only.
  The gate stays on the hull, because tightening a printed band is an owner ruling
  rather than a code edit; the intersection is REPORTED instead, as
  ``PressureDropCorroboration.inside_both_printed_bands``, and the reference-scale
  family-2 total sits inside both bands. ``opening_count_density_per_m2 = None`` is
  the reference case, so this loss law cannot require an opening population and
  does not read one (``tray_type.py:648-661``).
* THE DIFFUSION-LIMITED FILM DEGENERATES AT ITS OWN BOUNDARY. The Stefan-corrected
  form assumes an inert carrier; as the sum of interfacial partial pressures
  approaches the layer pressure the Stefan factor is singular, which is physically
  the crossover to PHY-022's heat-controlled (Faner) limit. That crossover is NOT
  selected here: the degeneracy is refused as a typed step rejection instead of
  being smoothed, and selecting between the diffusion-limited and heat-controlled
  regimes is named as work this packet does not do.
* CONSTANT-PROPERTY, DECLARED-INPUT CLOSURE. Every property is a caller-declared
  input -- one molar heat capacity for the gas mixture, one latent heat per
  species, and a Clausius-Clapeyron saturation law whose reference point and
  slope the caller declares. This module mints no property number and consumes no
  frozen property module, so it authenticates no property source. All of it is
  provenance class DECLARED_ENGINEERING_ASSUMPTION.
* PHY-023's a_pg IS DECLARED AT ONE STATE THE FROZEN TEXT DOES NOT COVER. The
  ledger prints a_pg at the two hexane endpoints only; PHY-032 generalizes the
  liquid-gas area to the shared ``a_Lg`` whenever external liquid exists. This
  module therefore DECLARES the reconciliation: one gas-side area for all
  interface transport, equal to the inherited dispatch area at the dry state and
  to the dispatch area times (1-S_L)^n whenever any external liquid exists. Said
  here rather than implied.
* THE WALL / FAST-BLOCK SEAM DOES NOT CONSERVE ENERGY WITHIN ONE MACRO STEP, AND
  THAT IS THE DECLARED COST OF THE RULED ONE-WAY PHY-008 SPLIT. The fast block
  credits the two wall couplings at the OLD wall temperature (a frozen parameter),
  while ``advance_wall_node`` debits the same two paths at the NEW one, so the
  composed macro step mis-accounts exactly
  ``(wall_to_gas_ua + wall_to_interface_ua) * (T_wall_new - T_wall_old)`` watts. On
  the shipped fixture's first step that is 528.6 W on a 5.286 K move, i.e. 1.16 % of
  the 45.4 kW steam-side input; it decays geometrically over a march and is EXACTLY
  0.0 at the C5a fixed point, where the wall temperature stops moving. The row residuals
  cannot see it, because each is self-consistent inside its own subsystem, which is
  why the gap is differenced explicitly and exposed as
  ``CellBoundaryLedger.wall_coupling_split_imbalance_w`` rather than left for a
  consumer to discover. Closing it would mean iterating the wall against the fast
  block, i.e. absorbing the slowest mode into the fast block, which FROZEN PHY-008
  forbids.
* PHY-023's AREA LAW IS THE MEMORYLESS PRIMARY-DRAINAGE FORM, AND A NON-MONOTONE
  FILM TRAJECTORY IS NOT FLAGGED HERE. PHY-023's second admissibility bullet
  (``release/physics_decisions.yaml:1511``) puts rewetting and later n-hexane
  condensation on a receded core outside this closure and under FROZEN PHY-036,
  whose response-outside-envelope clause requires such a trajectory to be FLAGGED
  rather than modelled. ``require_saturation_admissible`` implements the first
  bullet only (0 <= S_L < 1, rejected not clamped). A condensing (negative)
  attached-hexane flux at X_f > 0 grows the film non-monotonically and is outside
  that qualification; this packet does NOT flag it. The strictly excluded
  sub-case -- condensation initiating on a receded/dry core -- is structurally
  unreachable, because at X_f = 0 the branch deactivates hexane and pins its flux
  to exactly 0.0 by an identity row, so no X_f is ever silently reactivated from
  zero. If a mechanical flag is wanted it belongs as a post-convergence reported
  diagnostic or a pre-integration refusal, NEVER as an in-Newton orientation
  branch.
* THE DYNAMIC-EOS GAS VOLUME AND THE TRANSIENT FILM CAPACITY BOTH CLAIM THE WHOLE
  BED VOID, AND THE PARTITION IS DECLARED RATHER THAN MODELLED. The gas storage
  moles are computed over ``bed_volume * bed_void_fraction`` treated as constant in
  the external liquid saturation, while the film inventory capacities are
  ``bed_void_fraction * rho_liquid * bed_volume``, i.e. S_L = 1 would be a void full
  of liquid. Engage both switches at a high saturation and the two rows of the same
  DAE describe overlapping volumes. DECLARED SIMPLIFICATION, with its consequence
  stated: the gas capacity carries no partial with respect to saturation, which is
  exact for the model as written and wrong for the bed it represents. A
  liquid-occupancy correction ``(1 - S_L)`` on the gas volume is a physics change,
  not a code tidy, so it is booked here rather than taken.
* THE IDENTIFIABILITY PROHIBITION IS A LIMITATION AS WELL AS A ROUTED REFUSAL.
  Stated as the owner ruling's enactment annotation requires: this module offers NO
  joint inference of the bed and floor coefficients against a total-only target,
  and no fitting or calibration path of any kind. What it offers instead is the
  predeclared sensitivity sweep and per-element contributions on every evaluation.
  The refusal is selected by the observation's own type (see
  ``fit_series_coefficients_to_observation``), but no typed API can stop a caller
  bisecting a real-valued geometry field from OUTSIDE the module; the mitigation
  that carries weight is that the aggregate is never the only observable.
* PER-PROVENANCE GEOMETRY ADMISSIBILITY LIVES IN ITS OWN MODULE AND IS NOT CALLED
  FROM HERE. ``geometry_provenance_admissibility`` checks a deck's open-area
  fraction against the printed bands of its OWN declared ``geometry_provenance``
  rather than against the inherited provenance-blind hull. It reads this module's
  types one way and is never imported here, so the closure stays provenance-blind
  by construction and the single-home source scan keeps passing. What this module
  does enforce is the reconciliation between the two consumers of that one
  quantity: the floor element's throat area fraction must equal the declared deck's
  free area fraction, so the loss arithmetic cannot silently use an open fraction
  the deck never declared. A caller that skips the provenance gate therefore still
  gets the inherited hull through ``TrayTypeSpec``, and nothing narrower.
* NO PHYSICS REOPENED, NO WIRING, NO GATE. PHY-008's well-mixed default stands,
  the Pe 0.4-6 bracket stays a bracket, PHY-009's mesh-independence duties remain
  owed, and D4 nominals feed the design-basis tier only. The facade contract is
  untouched, no config is flipped, no steam source node is created (a gas source
  is origination; countercurrent decks are ascent paths), and the
  activation-to-fully-wet handoff stays outside this DAE with
  ``complete_particle_state_transport_selected`` false. No endpoint identity is
  accepted inside the cell solve; it arrives from the event side.
* ``evidentiary_standing_rank`` IS ORDINAL. It is readable to gate or report,
  through one reporting function, and is forbidden as a weight in any residual,
  average or blend (``tray_type.py:341-346``).

physically_qualifying: false. plant_predictive: false. No F-gate advances.
"""

from __future__ import annotations

import enum
import hashlib
import math
import struct
from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import ClassVar

from .dtdc_stack import TraySpec
from .tray_type import (
    CANONICAL_LOSS_RELATION,
    FAMILY2_SCHUMACHER_TAPER_BORE_AREA_RATIO_BAND,
    DeclaredGeometryProvenance,
    PressureLossRelation,
    TrayFloorFamily,
    TrayTypeSpec,
    gas_side_contact_area_m2,
    require_type_matches_tray,
)


#: This packet mints its own identity; the inherited tray-type module carries no
#: schema id and is deliberately not retrofitted with one
#: (``tray_type.py:159-163``).
CELL_CLOSURE_MAGIC = "GT-PS-2-CELL-LAYER-IMPLICIT-CLOSURE"
CELL_CLOSURE_VERSION = 1
CELL_CLOSURE_DIGEST_DOMAIN = "GT-PS-2/cell-layer-implicit-closure/v1"

#: Provenance class of every declared value that reaches this closure (D0).
DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"

#: A physical constant, not a fitted parameter and not a band.
UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.314462618

#: The seven frozen carrier extensives, mirrored as value strings for boundary
#: inspectability rather than imported, so this module stays detached from the
#: seam modules the frozen source guards scan for.
EXTENSIVE_NAMES: tuple[str, ...] = (
    "dry_matter_kg",
    "residual_oil_label_kg",
    "attached_hexane_kg",
    "internal_hexane_kg",
    "external_water_kg",
    "retained_water_kg",
    "common_datum_energy_j",
)

#: THE TWO LENGTH SCALES, NAMED, AND BOTH NAMES REAL. The hydraulic one is a field
#: of ``BedErgunElement`` -- the element that carries it -- and of no other type
#: here. The particle-model one is ``particle_radius_m`` -- the name the INHERITED
#: particle side actually uses: a declared field of the particle-state parameters
#: (``src/dtdc_simulator/core2/sorption_interface.py:133``) and the argument name of
#: the inherited gas-side area law (``src/dtdc_simulator/core2/area.py:109``). It is
#: written here only so that the separation is machine-checkable against the REAL
#: name rather than against a string that names nothing anywhere. No value of the
#: particle-model scale is read, stored or derived in this module.
HYDRAULIC_LENGTH_SCALE_FIELD = "hydraulic_granule_diameter_m"
PARTICLE_MODEL_LENGTH_SCALE_FIELD = "particle_radius_m"
#: FROZEN PHY-034's intraparticle diffusion radius, carried ONLY as the tripwire
#: value the hydraulic residual must never be handed
#: (``release/physics_decisions.yaml:1223-1224``).
INTRAPARTICLE_DIFFUSION_RADIUS_M = 0.208e-3

#: The frozen Ergun per-layer band and the independently corroborated printed
#: band (``tray_type.py:117-123``). They OVERLAP on 2.0-5.0 kPa, so the hull below
#: has printed support throughout its length: no new band is introduced, and the
#: hull is stated as a hull rather than presented as a measurement.
FROZEN_ERGUN_LAYER_DROP_BAND_PA: tuple[float, float] = (2000.0, 5600.0)
CORROBORATED_LAYER_DROP_BAND_PA: tuple[float, float] = (1000.0, 5000.0)
CORROBORATED_LAYER_DROP_HULL_PA: tuple[float, float] = (1000.0, 5600.0)
#: The INTERSECTION of the same two printed bands -- the strict reading of the
#: design basis's "inside the family-independent corroboration window". Narrower
#: than the hull the gate enforces; reported, never enforced, because tightening a
#: printed band is an owner ruling and not a code edit.
DOUBLY_PRINTED_LAYER_DROP_INTERSECTION_PA: tuple[float, float] = (2000.0, 5000.0)

#: FROZEN PHY-008's AUDITED SEPARATION, as a bracket rather than a claim: the metal
#: time constant exceeds the per-tray solid residence by 1.4-8.3x
#: (``release/physics_decisions.yaml:953`` -- "~6.4-38.5 min for a 6 m steam-chest
#: tray versus ~27.8 min total and ~4.6 min per-tray solid residence"). The wall
#: row's ``is_slowest_mode`` is DERIVED against this bracket from the declared
#: capacity and conductances; it is not asserted.
PHY008_WALL_RESIDENCE_RATIO_BRACKET: tuple[float, float] = (1.4, 8.3)

#: The two standard Ergun coefficients. FROZEN PHY-017 calls Ergun "a standard
#: non-fitted correlation"; these are that correlation's own constants and are
#: not calibration knobs.
ERGUN_VISCOUS_COEFFICIENT = 150.0
ERGUN_INERTIAL_COEFFICIENT = 1.75

#: THE RULED SERIES (owner ruling 2026-08-12). Each family maps to an ORDERED
#: tuple of element KINDS, tagged with the inherited relation enum. Family 2 ships
#: bed Ergun resistance FOLLOWED BY the tapered-bore floor passage -- the ruled
#: pair, no longer a provisional single relation.
PRESSURE_SERIES_RULING = "GT-PS-2 owner ruling, 2026-08-12: ordered series of typed loss elements"
RULED_ELEMENT_KIND_SERIES: Mapping[TrayFloorFamily, tuple[PressureLossRelation, ...]] = (
    MappingProxyType(
        {
            TrayFloorFamily.UNPERFORATED_HEATED_DISC: (PressureLossRelation.DOME_BYPASS,),
            TrayFloorFamily.STRIPPING_DECK: (
                PressureLossRelation.PACKED_BED_ERGUN,
                PressureLossRelation.SPARGE_NOZZLE,
            ),
            TrayFloorFamily.DIRECT_STEAM_SPARGE: (PressureLossRelation.SPARGE_NOZZLE,),
            TrayFloorFamily.AIR_DISTRIBUTION: (PressureLossRelation.AIR_PLENUM,),
        }
    )
)

#: The inherited single selector per family is the family's PRIMARY element kind.
#: Asserted rather than assumed, so a change to either side is caught at import.
for _family, _kinds in RULED_ELEMENT_KIND_SERIES.items():  # pragma: no cover - invariant
    if _kinds[0] is not CANONICAL_LOSS_RELATION[_family]:
        raise AssertionError(
            "the inherited per-family selector must be the family's primary element kind"
        )
del _family, _kinds

#: Families with no member in the accepted reference machine. Constructible and
#: invariant-checked, but unexercised by plant-anchored geometry.
UNANCHORED_BY_REFERENCE_MACHINE_FAMILIES: tuple[TrayFloorFamily, ...] = (
    TrayFloorFamily.DIRECT_STEAM_SPARGE,
    TrayFloorFamily.AIR_DISTRIBUTION,
)

#: THE PREDECLARED SENSITIVITY RANGE of the floor discharge coefficient, ruled to
#: remain a DECLARED_ENGINEERING_ASSUMPTION. Both edges are textbook
#: discharge-coefficient limits, neither is fitted, and nothing between them is
#: minted as a nominal: a caller declares its own value inside the range and the
#: range is swept by ``sweep_floor_discharge_coefficient``.
#:
#: Lower edge 0.61 -- the classical sharp-edged thin-plate orifice. Upper edge
#: 0.98 -- a well-formed converging nozzle. The printed family-2 archetype bore is
#: a CONVERGING TAPER (printed area ratio f:g about 2.0-1.5), so its true
#: coefficient sits between a sharp-edged plate and an ideal converging nozzle;
#: the range brackets that interval from both sides rather than asserting a point
#: inside it. A caller may NARROW the range with its own evidence; widening it is
#: an owner ruling, not a code edit.
FLOOR_DISCHARGE_COEFFICIENT_SENSITIVITY_RANGE: tuple[float, float] = (0.61, 0.98)
#: The printed taper area ratio, imported from the inherited typed-floor module
#: and carried here for the GEOMETRY check only. It never enters the loss
#: arithmetic, which uses the throat.
TAPER_AREA_RATIO_BAND = FAMILY2_SCHUMACHER_TAPER_BORE_AREA_RATIO_BAND

#: The fast-block unknowns, in assembly order. THE COUNT IS DERIVED FROM THIS
#: TUPLE EVERYWHERE; no dimension is written as a literal, which is what makes
#: the ruled successor's one extra state variable a name and a row rather than a
#: restructuring.
UNKNOWN_NAMES: tuple[str, ...] = (
    "gas_temperature_k",
    "gas_hexane_mole_fraction",
    "gas_water_mole_fraction",
    "layer_pressure_pa",
    "interface_temperature_k",
    "external_hexane_saturation",
    "external_water_saturation",
)
_N_UNKNOWNS = len(UNKNOWN_NAMES)
UNKNOWN_INDEX: Mapping[str, int] = MappingProxyType(
    {name: index for index, name in enumerate(UNKNOWN_NAMES)}
)


class CellSubsystem(enum.Enum):
    """The four subsystems. There is no fifth."""

    GAS = "gas"
    LAYER_PRESSURE = "layer_pressure"
    FILM_INTERFACE = "film_interface"
    #: One differential row, advanced separately. NOT in the fast Newton block.
    WALL_METAL = "wall_metal"

    @property
    def is_fast_block(self) -> bool:
        return self is not CellSubsystem.WALL_METAL


class ResidualRow(enum.Enum):
    """Every row of the cell DAE, named, with its owning subsystem."""

    GAS_HEXANE_COMPONENT = "gas_hexane_component"
    GAS_WATER_COMPONENT = "gas_water_component"
    GAS_ENERGY = "gas_energy"
    LAYER_PRESSURE_LOSS = "layer_pressure_loss"
    INTERFACE_ENERGY = "interface_energy"
    FILM_HEXANE_INVENTORY = "film_hexane_inventory"
    FILM_WATER_INVENTORY = "film_water_inventory"
    #: PHY-008. The slowest mode; advanced at the macro-step scale.
    WALL_METAL_INTERNAL_ENERGY = "wall_metal_internal_energy"

    @property
    def subsystem(self) -> CellSubsystem:
        return _ROW_SUBSYSTEM[self]

    @property
    def is_differential_by_default(self) -> bool:
        """True only for the wall row; the rest are algebraic by default."""

        return self is ResidualRow.WALL_METAL_INTERNAL_ENERGY


_ROW_SUBSYSTEM: Mapping[ResidualRow, CellSubsystem] = MappingProxyType(
    {
        ResidualRow.GAS_HEXANE_COMPONENT: CellSubsystem.GAS,
        ResidualRow.GAS_WATER_COMPONENT: CellSubsystem.GAS,
        ResidualRow.GAS_ENERGY: CellSubsystem.GAS,
        ResidualRow.LAYER_PRESSURE_LOSS: CellSubsystem.LAYER_PRESSURE,
        ResidualRow.INTERFACE_ENERGY: CellSubsystem.FILM_INTERFACE,
        ResidualRow.FILM_HEXANE_INVENTORY: CellSubsystem.FILM_INTERFACE,
        ResidualRow.FILM_WATER_INVENTORY: CellSubsystem.FILM_INTERFACE,
        ResidualRow.WALL_METAL_INTERNAL_ENERGY: CellSubsystem.WALL_METAL,
    }
)

#: The fast Newton block: three subsystems, seven rows, one per unknown.
FAST_BLOCK_ROWS: tuple[ResidualRow, ...] = tuple(
    row for row in ResidualRow if row.subsystem.is_fast_block
)
#: The wall row, alone, outside the fast block.
WALL_ROWS: tuple[ResidualRow, ...] = tuple(
    row for row in ResidualRow if not row.subsystem.is_fast_block
)

if len(FAST_BLOCK_ROWS) != _N_UNKNOWNS:  # pragma: no cover - structural invariant
    raise AssertionError("the fast block must be square in its own unknowns")


class GasClosure(enum.Enum):
    """FROZEN PHY-004: algebraic nominal, dynamic-EOS as a switchable mode."""

    QUASI_STEADY = "quasi_steady"
    DYNAMIC_EOS = "dynamic_eos"


class FilmClosure(enum.Enum):
    """The film inventory rows: algebraic by default, switchable to transient."""

    QUASI_STEADY_ALGEBRAIC = "quasi_steady_algebraic"
    TRANSIENT_ACCUMULATION = "transient_accumulation"


class WallClosure(enum.Enum):
    """FROZEN PHY-008: dynamic node selected, quasi-steady retained as the limit."""

    DYNAMIC_NODE = "dynamic_node"
    QUASI_STEADY_LIMIT = "quasi_steady_limit"


class AccumulationMode(enum.Enum):
    """The C5a switch. ZERO_ACCUMULATION drives every storage term to exactly 0.

    C5a is the zero-accumulation limit of the identical representation, so it is
    reached by parameter through the SAME code path rather than by a separate
    steady-state implementation (``D7:105-108``).
    """

    TRANSIENT = "transient"
    ZERO_ACCUMULATION = "zero_accumulation"


class PressureNetworkMode(enum.Enum):
    """FROZEN PHY-017's nested limits, reachable by parameter alone (``:1217``).

    The listed alternatives of the frozen decision ARE limits of the nodal
    network: resistances to zero recovers the negligible-drop oracle exactly
    (``resistance_scale = 0.0`` under ``NODAL_HYDRAULIC``), frozen
    flows/resistances recover a prescribed axial profile
    (``FROZEN_RESISTANCE``), and a single node recovers prescribed uniform
    (``SINGLE_NODE_UNIFORM``).
    """

    NODAL_HYDRAULIC = "nodal_hydraulic"
    FROZEN_RESISTANCE = "frozen_resistance"
    SINGLE_NODE_UNIFORM = "single_node_uniform"


class InterfaceBranch(enum.Enum):
    """FROZEN PHY-007's four-branch active set (``:902-916``). Exactly four.

    Selected by EXACT comparison on the two external solid-side inventories the
    packet-carrier layer owns. Read, never re-derived, never smoothed.
    """

    DRY_SURFACE = "dry_surface"
    EXTERNAL_HEXANE = "external_hexane"
    EXTERNAL_WATER = "external_water"
    TWO_EXTERNAL_LIQUIDS = "two_external_liquids"

    @property
    def hexane_is_active(self) -> bool:
        return self in (
            InterfaceBranch.EXTERNAL_HEXANE,
            InterfaceBranch.TWO_EXTERNAL_LIQUIDS,
        )

    @property
    def water_is_active(self) -> bool:
        return self in (
            InterfaceBranch.EXTERNAL_WATER,
            InterfaceBranch.TWO_EXTERNAL_LIQUIDS,
        )

    @property
    def has_external_liquid(self) -> bool:
        return self is not InterfaceBranch.DRY_SURFACE


class FilmAreaExponent(enum.Enum):
    """FROZEN PHY-023's fixed model-form bracket. NOT a calibration knob.

    "The exponent is not a calibration knob" (``:1512``), so it is typed rather
    than numeric: an arbitrary fitted exponent is unrepresentable.
    """

    #: a_hg/a_b = 1-S -- lower, strong shielding.
    STRONG_SHIELDING = 1.0
    #: a_hg/a_b = (1-S)^(2/3) -- the frozen reference isotropic law.
    ISOTROPIC_REFERENCE = 2.0 / 3.0
    #: a_hg/a_b = sqrt(1-S) -- upper, capillary.
    CAPILLARY_UPPER = 0.5
    #: a_hg/a_b = 1 -- DIAGNOSTIC ONLY, never a modelling choice.
    DIAGNOSTIC_UNITY = 0.0

    @property
    def is_diagnostic_only(self) -> bool:
        return self is FilmAreaExponent.DIAGNOSTIC_UNITY

    @property
    def derivative_diverges_at_saturation(self) -> bool:
        """True where d/dS of (1-S)^n is unbounded as S -> 1."""

        return self.value < 1.0 and self.value > 0.0


# --- the typed refusal family ----------------------------------------------
# Two tiers, and the tiers are the point. A SEED refusal is terminal: the ruled
# behaviour is that an out-of-domain excursion is refused at seed assembly and
# never integrated. An ITERATE refusal is a step rejection that propagates out of
# the implicit solve. Neither tier ever repairs a state.


class CellClosureError(ValueError):
    """Base exception for a refused cell-layer closure operation."""


class StepRejected(CellClosureError):
    """A trial iterate is inadmissible: the step is REJECTED, never repaired.

    Every subclass propagates out of the implicit solve as a step rejection. No
    projection back into the feasible set exists anywhere in this module.
    """


class SaturationAdmissibilityError(StepRejected):
    """FROZEN PHY-023: reject S_L outside [0, 1) rather than clamping it."""


class SingularAreaDerivativeError(StepRejected):
    """The PHY-023 area derivative reached its singular limit at S -> 1.

    The derivative -(2/3) a_b (1-S)^(-1/3) genuinely diverges. The analytic
    Jacobian carries the divergence; at the point where it stops being
    representable in binary64 this refuses rather than emitting inf or nan.
    """


class InterfaceStateError(StepRejected):
    """The diffusion-limited film form reached its own degeneracy.

    Raised when the interfacial partial-pressure sum reaches the layer pressure,
    which is the crossover to the heat-controlled limit this packet does not
    select. Refused, not smoothed.
    """


class NewtonConvergenceError(StepRejected):
    """The undamped Newton iteration did not converge in the admitted count.

    There is no safeguard, no damping and no line search, because the admitted
    safeguards would all have to clip or project. Non-convergence is refused.
    """


class SingularJacobianError(StepRejected):
    """The analytic Jacobian is numerically singular at this iterate."""


class SeedAdmissibilityError(CellClosureError):
    """TERMINAL. The seed was refused before any integration took place.

    Deliberately NOT a ``StepRejected``: a rejected step invites a retry, and the
    ruled behaviour for an out-of-domain seed is that it is never integrated.
    """


class LengthScaleConfusionError(CellClosureError):
    """The hydraulic residual was handed an intraparticle length scale."""


class PressureBandError(CellClosureError):
    """A converged layer drop fell outside the corroborated printed window."""


class LossRelationNotImplementedError(CellClosureError):
    """The dispatched element kind has no element type in this packet."""


class LossSeriesError(CellClosureError):
    """The ordered loss series disagrees with the family's ruled element kinds."""


class JointIdentifiabilityRefusedError(CellClosureError):
    """Raised on any attempt to infer two series coefficients from one total.

    The owner ruling of 2026-08-12 forbids it: "Bed and floor coefficients shall
    not be jointly inferred from total pressure drop alone."
    """


class CoefficientInferenceNotImplementedError(CellClosureError):
    """Raised when any coefficient inference is asked for at all.

    Distinct from ``JointIdentifiabilityRefusedError`` on purpose. The ruled
    prohibition is specifically about two coefficients against ONE AGGREGATE;
    this is the wider and duller fact that no fitting path exists here even
    against a decomposed observation, so the two are not conflated.
    """


class WallNodeError(CellClosureError):
    """The single wall row is not solvable as posed."""


class ClosureConfigurationError(CellClosureError):
    """The declared inputs are inconsistent with each other."""


# --- exact forward differentiation -----------------------------------------
# The Jacobian is ANALYTIC: every rule below is a closed-form derivative applied
# by the chain rule, evaluated to machine precision. It is not a finite-difference
# surrogate, and the two are distinguishable -- the companion suite checks the
# assembled Jacobian AGAINST finite differences, which would be circular if the
# assembly were itself differencing. Carrying one partial per declared unknown is
# also what makes the ruled successor's extra state variable a name rather than a
# restructuring.


class _Dual:
    """A value and its exact partials with respect to the declared unknowns."""

    __slots__ = ("value", "partials")

    def __init__(self, value: float, partials: tuple[float, ...]) -> None:
        self.value = value
        self.partials = partials

    def __add__(self, other: _Dual | float) -> _Dual:
        if isinstance(other, _Dual):
            return _Dual(
                self.value + other.value,
                tuple(a + b for a, b in zip(self.partials, other.partials, strict=True)),
            )
        return _Dual(self.value + other, self.partials)

    __radd__ = __add__

    def __neg__(self) -> _Dual:
        return _Dual(-self.value, tuple(-a for a in self.partials))

    def __sub__(self, other: _Dual | float) -> _Dual:
        if isinstance(other, _Dual):
            return _Dual(
                self.value - other.value,
                tuple(a - b for a, b in zip(self.partials, other.partials, strict=True)),
            )
        return _Dual(self.value - other, self.partials)

    def __rsub__(self, other: float) -> _Dual:
        return _Dual(other - self.value, tuple(-a for a in self.partials))

    def __mul__(self, other: _Dual | float) -> _Dual:
        if isinstance(other, _Dual):
            return _Dual(
                self.value * other.value,
                tuple(
                    a * other.value + self.value * b
                    for a, b in zip(self.partials, other.partials, strict=True)
                ),
            )
        return _Dual(self.value * other, tuple(a * other for a in self.partials))

    __rmul__ = __mul__

    def __truediv__(self, other: _Dual | float) -> _Dual:
        if isinstance(other, _Dual):
            if other.value == 0.0:
                raise StepRejected("division by an exactly zero intermediate")
            inverse = 1.0 / other.value
            scale = self.value * inverse * inverse
            return _Dual(
                self.value * inverse,
                tuple(
                    a * inverse - b * scale
                    for a, b in zip(self.partials, other.partials, strict=True)
                ),
            )
        if other == 0.0:
            raise StepRejected("division by an exactly zero declared input")
        return _Dual(self.value / other, tuple(a / other for a in self.partials))

    def __rtruediv__(self, other: float) -> _Dual:
        if self.value == 0.0:
            raise StepRejected("division by an exactly zero intermediate")
        inverse = 1.0 / self.value
        scale = other * inverse * inverse
        return _Dual(other * inverse, tuple(-a * scale for a in self.partials))


def _constant(value: float) -> _Dual:
    return _Dual(value, (0.0,) * _N_UNKNOWNS)


def _variable(value: float, index: int) -> _Dual:
    partials = [0.0] * _N_UNKNOWNS
    partials[index] = 1.0
    return _Dual(value, tuple(partials))


def _exp(argument: _Dual) -> _Dual:
    # TWO DISTINCT INGRESSES, BOTH TYPED. ``math.exp`` RAISES ``OverflowError`` for a
    # finite argument above about 709.78 and never returns an infinity there, so the
    # ``isfinite`` test below cannot see that case -- it never receives a value. The
    # test is still live and still needed: ``math.exp`` returns inf for inf and nan
    # for nan without raising, and a later float multiplication can overflow
    # silently. Both paths therefore fail closed as the same typed step rejection.
    try:
        value = math.exp(argument.value)
    except OverflowError as error:
        raise StepRejected(
            "an exponential intermediate is not representable"
        ) from error
    if not math.isfinite(value):
        raise StepRejected("an exponential intermediate is not representable")
    return _Dual(value, tuple(value * a for a in argument.partials))


def _from_scalar_rule(argument: _Dual, value: float, derivative: float) -> _Dual:
    """Lift a closed-form (value, derivative) pair onto the chain rule."""

    if not math.isfinite(value) or not math.isfinite(derivative):
        raise StepRejected("a closed-form kernel is not representable at this iterate")
    return _Dual(value, tuple(derivative * a for a in argument.partials))


#: Below this magnitude the two transport kernels are evaluated by their own
#: Taylor series. THIS IS NOT A MOLLIFIER AND NOT AN EPSILON SHELL: the series IS
#: the analytic function and its analytic derivative, truncated below binary64
#: resolution (the first dropped term is under 1e-19 relative at the threshold).
#: It exists because the closed forms differ two nearly equal quantities there,
#: and a cancellation artefact would be indistinguishable from a Jacobian error.
_KERNEL_SERIES_THRESHOLD = 2.0**-10


def _stefan_scalar(b: float) -> tuple[float, float]:
    """log1p(B)/B and its derivative -- the Stefan-flow correction factor.

    Exactly 1 at B = 0 (a removable singularity, taken as the exact limit, with
    derivative exactly -1/2), so the low-flux independent-linear film is
    recovered EXACTLY rather than approached.
    """

    if b <= -1.0:
        raise InterfaceStateError(
            "the Stefan correction is undefined once the bulk carrier vanishes"
        )
    if b == 0.0:
        return (1.0, -0.5)
    if abs(b) < _KERNEL_SERIES_THRESHOLD:
        value = (
            1.0
            - b / 2.0
            + b**2 / 3.0
            - b**3 / 4.0
            + b**4 / 5.0
            - b**5 / 6.0
            + b**6 / 7.0
        )
        derivative = (
            -0.5
            + 2.0 * b / 3.0
            - 3.0 * b**2 / 4.0
            + 4.0 * b**3 / 5.0
            - 5.0 * b**4 / 6.0
            + 6.0 * b**5 / 7.0
        )
        return (value, derivative)
    log_term = math.log1p(b)
    value = log_term / b
    derivative = (b / (1.0 + b) - log_term) / (b * b)
    return (value, derivative)


def _ackermann_scalar(z: float) -> tuple[float, float]:
    """z/expm1(z) and its derivative -- the thermodynamically consistent
    heat-transfer correction PHY-022's coupled film requires.

    Exactly 1 at z = 0 (no blowing), derivative exactly -1/2.
    """

    if z == 0.0:
        return (1.0, -0.5)
    if abs(z) < _KERNEL_SERIES_THRESHOLD:
        value = 1.0 - z / 2.0 + z**2 / 12.0 - z**4 / 720.0 + z**6 / 30240.0
        derivative = -0.5 + z / 6.0 - z**3 / 180.0 + z**5 / 5040.0
        return (value, derivative)
    # Same two ingresses as ``_exp``: a raise for a large finite argument, and a
    # silent inf/nan for a non-finite one. Both are the same typed refusal.
    try:
        lower = math.expm1(z)
        upper = math.exp(z)
    except OverflowError as error:
        raise InterfaceStateError(
            "the blowing correction is not representable here"
        ) from error
    if lower == 0.0 or not math.isfinite(lower):
        raise InterfaceStateError("the blowing correction is not representable here")
    if not math.isfinite(upper):
        raise InterfaceStateError("the blowing correction is not representable here")
    return (z / lower, (lower - z * upper) / (lower * lower))


def _stefan_factor(b: _Dual) -> _Dual:
    value, derivative = _stefan_scalar(b.value)
    return _from_scalar_rule(b, value, derivative)


def _ackermann_factor(z: _Dual) -> _Dual:
    value, derivative = _ackermann_scalar(z.value)
    return _from_scalar_rule(z, value, derivative)


# --- FROZEN PHY-023 / PHY-032: ONE area scalar, and its diverging derivative --


def require_saturation_admissible(saturation: float, *, at_seed: bool) -> None:
    """FROZEN PHY-023 admissibility: 0 <= S_L < 1, REJECTED not clamped (``:1510``).

    At the seed the refusal is terminal, because the ruled behaviour is that an
    out-of-domain state is refused at seed assembly and never integrated. At an
    iterate it is a step rejection. Neither path returns a repaired value: there
    is no projection onto the feasible set here, as a safeguard or otherwise.
    """

    if type(saturation) is not float or not math.isfinite(saturation):
        message = "total external liquid saturation must be a finite binary64 value"
        raise SeedAdmissibilityError(message) if at_seed else SaturationAdmissibilityError(message)
    if saturation < 0.0 or saturation >= 1.0:
        message = (
            f"total external liquid saturation {saturation!r} is outside the frozen "
            "admissible set 0 <= S_L < 1; the frozen decision rejects an infeasible "
            "feed/bed state rather than clamping it"
        )
        raise SeedAdmissibilityError(message) if at_seed else SaturationAdmissibilityError(message)


def require_component_saturations_admissible(
    *,
    branch: InterfaceBranch,
    hexane_saturation: float,
    water_saturation: float,
    at_seed: bool,
) -> float:
    """The PER-COMPONENT feasible set, at either tier, from ONE implementation.

    The total alone would admit a negative component cancelled by a larger
    positive one, which is not a bed state: at the co-located branch
    ``(S_h, S_w) = (-0.30, 0.50)`` sums to an admissible 0.20. So each PRESENT
    component is tested on its own sign and the total is then tested against the
    frozen set.

    Both tiers call THIS function, so the seed tier and the iterate tier cannot
    drift apart -- which is exactly how they had drifted: the seed carried the
    per-component test and the iterate did not.

    An ABSENT component is deliberately not tested. Its coordinate is not an
    inventory: it is pinned by an exact identity row, the branch gives it a
    structural zero flux, and the seed tier already requires it to be exactly
    zero before any integration. Constraining it at the iterate tier as well would
    also refuse the central-difference probes the companion suite differentiates
    the identity rows with, i.e. it would bound the Jacobian gate rather than the
    physics.

    Returns the branch-relevant total so a caller does not recompute it. Refuses;
    never repairs, never projects.
    """

    if type(branch) is not InterfaceBranch:
        raise ClosureConfigurationError("branch must be an InterfaceBranch member")
    total = 0.0
    if branch.hexane_is_active:
        if hexane_saturation < 0.0:
            message = "a present external hexane film cannot carry a negative saturation"
            raise (
                SeedAdmissibilityError(message)
                if at_seed
                else SaturationAdmissibilityError(message)
            )
        total += hexane_saturation
    if branch.water_is_active:
        if water_saturation < 0.0:
            message = "a present external water film cannot carry a negative saturation"
            raise (
                SeedAdmissibilityError(message)
                if at_seed
                else SaturationAdmissibilityError(message)
            )
        total += water_saturation
    if branch.has_external_liquid:
        require_saturation_admissible(total, at_seed=at_seed)
    return total


def film_area_saturation_factor(*, saturation: float, exponent: FilmAreaExponent) -> float:
    """FROZEN PHY-023 / PHY-032 area factor (1 - S_L)^n. One scalar, no per-component form."""

    if type(exponent) is not FilmAreaExponent:
        raise CellClosureError("exponent must be a FilmAreaExponent member")
    require_saturation_admissible(saturation, at_seed=False)
    if exponent is FilmAreaExponent.DIAGNOSTIC_UNITY:
        return 1.0
    return (1.0 - saturation) ** exponent.value


def film_area_saturation_derivative(
    *, reference_area_m2: float, saturation: float, exponent: FilmAreaExponent
) -> float:
    """d/dS of A_ref (1 - S)^n = -n A_ref (1 - S)^(n-1).

    THE DIVERGENCE IS THE POINT. At the frozen reference exponent this is
    -(2/3) a_b (1-S)^(-1/3), which grows without bound as S -> 1. The analytic
    Jacobian carries exactly this expression; no finite surrogate is substituted.
    Where the magnitude stops being representable in binary64 the function raises
    ``SingularAreaDerivativeError`` rather than returning an infinity, so nothing
    downstream can consume a silent inf or nan.
    """

    if type(exponent) is not FilmAreaExponent:
        raise CellClosureError("exponent must be a FilmAreaExponent member")
    # The exported entry point validates its AREA too. Without this a negative area
    # returns a POSITIVE derivative for a shrinking-area law -- the divergence
    # pointing the wrong way with no refusal -- and an int area raises a bare
    # TypeError out of a path documented to refuse by type.
    _require_positive("reference_area_m2", reference_area_m2)
    require_saturation_admissible(saturation, at_seed=False)
    if exponent is FilmAreaExponent.DIAGNOSTIC_UNITY:
        return 0.0
    order = exponent.value
    base = 1.0 - saturation
    derivative = -order * reference_area_m2 * base ** (order - 1.0)
    if not math.isfinite(derivative):
        raise SingularAreaDerivativeError(
            f"the frozen area derivative -{order!r}*a*(1-S)^({order - 1.0!r}) is not "
            f"representable at S = {saturation!r}: the law diverges at the saturated "
            "limit and this refuses rather than emitting a non-finite Jacobian entry"
        )
    return derivative


def _saturated_area_dual(
    reference_area_m2: float, saturation: _Dual, exponent: FilmAreaExponent
) -> _Dual:
    """The ONE gas-side area scalar, differentiated by the frozen closed form."""

    if exponent is FilmAreaExponent.DIAGNOSTIC_UNITY:
        return _constant(reference_area_m2)
    value = reference_area_m2 * film_area_saturation_factor(
        saturation=saturation.value, exponent=exponent
    )
    derivative = film_area_saturation_derivative(
        reference_area_m2=reference_area_m2, saturation=saturation.value, exponent=exponent
    )
    return _from_scalar_rule(saturation, value, derivative)


# --- declared-input validation ---------------------------------------------


def _is_negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _require_binary64(name: str, value: float) -> None:
    if type(value) is not float:
        raise ClosureConfigurationError(f"{name} must be an exact binary64 float")
    if not math.isfinite(value):
        raise ClosureConfigurationError(f"{name} must be finite")
    if _is_negative_zero(value):
        raise ClosureConfigurationError(f"{name} must use canonical positive zero")


def _require_nonnegative(name: str, value: float) -> None:
    _require_binary64(name, value)
    if value < 0.0:
        raise ClosureConfigurationError(f"{name} must be non-negative")


def _require_positive(name: str, value: float) -> None:
    _require_binary64(name, value)
    if value <= 0.0:
        raise ClosureConfigurationError(f"{name} must be strictly positive")


def _require_fraction(name: str, value: float) -> None:
    _require_binary64(name, value)
    if value < 0.0 or value >= 1.0:
        raise ClosureConfigurationError(f"{name} must lie in [0.0, 1.0)")


def _require_nested_limit_switch(name: str, value: float) -> None:
    """A TWO-VALUED switch, not a continuous multiplier.

    ``resistance_scale`` exists to reach FROZEN PHY-017's negligible-drop nested
    limit by parameter alone: 0.0 is that limit and 1.0 is the frozen network. A
    real-valued multiplier on a correlation this module declares NON-FITTED would
    be a coefficient by another name -- unbounded, dimensionless, with no physical
    referent and no printed source -- so the enforced domain is exactly the two
    documented values. Refused, never clipped.
    """

    _require_nonnegative(name, value)
    if value not in (0.0, 1.0):
        raise ClosureConfigurationError(
            f"{name} is the two-valued nested-limit switch (0.0 recovers the "
            f"negligible-drop limit, 1.0 is the frozen network); {value!r} would be "
            "a continuous multiplier on a declared non-fitted correlation"
        )


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _float_part(value: float) -> bytes:
    canonical = 0.0 if value == 0.0 else value
    return struct.pack(">d", canonical)


# --- declared inputs -------------------------------------------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class GridAuthority:
    """Discretization identity, exposed for boundary inspectability.

    This is the pair the caller's pre-integration gate is keyed on, and the grid
    authority the deferred conservation proof needs. It is NOT the macro step:
    ``CellClosureInputs.macro_step_s`` is the caller-supplied solids-paced step
    the fast block solves inside.
    """

    cells: int
    dt_s: float

    def __post_init__(self) -> None:
        if type(self.cells) is not int or self.cells <= 0:
            raise ClosureConfigurationError("cells must be a strictly positive int")
        _require_positive("dt_s", self.dt_s)


@dataclass(frozen=True, slots=True, kw_only=True)
class CellGeometry:
    """The cell's geometry, CONSUMED from the inherited typed dispatch.

    No area is re-derived here. ``gas_side_reference_area_m2`` calls the inherited
    tray-type dispatch, so the six reference trays' areas stay bit-identical, and
    the gas-side name is the ruled one throughout.
    """

    tray: TraySpec
    tray_type: TrayTypeSpec

    def __post_init__(self) -> None:
        if type(self.tray) is not TraySpec:
            raise ClosureConfigurationError("tray must be an exact TraySpec")
        if type(self.tray_type) is not TrayTypeSpec:
            raise ClosureConfigurationError("tray_type must be an exact TrayTypeSpec")
        # Fail closed on a typed floor that disagrees with its tray's passage.
        require_type_matches_tray(self.tray, self.tray_type)

    @property
    def gas_side_reference_area_m2(self) -> float:
        """The inherited dispatch's area, unmodified."""

        return gas_side_contact_area_m2(self.tray, self.tray_type)

    @property
    def bed_volume_m3(self) -> float:
        return self.tray.bed_volume_m3

    @property
    def cross_section_m2(self) -> float:
        return self.tray.cross_section_m2

    def void_volume_m3(self, void_fraction: float) -> float:
        return self.bed_volume_m3 * void_fraction


@dataclass(frozen=True, slots=True, kw_only=True)
class PressureLossContribution:
    """One element's own drop and its own analytic derivatives.

    Every element computes its own closed-form partials. Nothing here is
    differenced, and the series never re-differentiates a composed total.
    """

    drop_pa: float
    d_drop_d_superficial_velocity: float
    d_drop_d_gas_density: float
    kind: PressureLossRelation
    element_label: str

    def __post_init__(self) -> None:
        _require_binary64("drop_pa", self.drop_pa)
        _require_binary64(
            "d_drop_d_superficial_velocity", self.d_drop_d_superficial_velocity
        )
        _require_binary64("d_drop_d_gas_density", self.d_drop_d_gas_density)
        if type(self.kind) is not PressureLossRelation:
            raise ClosureConfigurationError("kind must be a PressureLossRelation member")
        if type(self.element_label) is not str or not self.element_label.strip():
            raise ClosureConfigurationError("element_label must be a nonblank string")


@dataclass(frozen=True, slots=True, kw_only=True)
class BedErgunElement:
    """THE HYDRAULIC LENGTH SCALE LIVES HERE, and nothing particle-scale does.

    The grounded Ergun bed resistance, evaluated with the true external granule
    hydraulic diameter. This type carries no field naming a particle dimension,
    and its constructor refuses a hydraulic diameter at or below the frozen
    intraparticle diffusion scale -- the value FROZEN PHY-017's consistency
    requirement says would give a spurious ~0.4-0.8 bar drop.

    ``opening_count_density_per_m2`` is deliberately absent from this element and
    from the floor element: None is the reference case on the inherited typed
    floors, so neither law may require an opening population, and neither reads one.
    """

    hydraulic_granule_diameter_m: float
    bed_void_fraction: float
    flow_length_m: float
    gas_dynamic_viscosity_pa_s: float
    #: The zero-resistance nested limit, reachable by parameter alone: 0.0 gives a
    #: drop of exactly zero and recovers the negligible-drop oracle. EXACTLY 0.0 or
    #: 1.0 -- a switch, never a continuous multiplier (see
    #: ``_require_nested_limit_switch``).
    resistance_scale: float = 1.0

    kind: ClassVar[PressureLossRelation] = PressureLossRelation.PACKED_BED_ERGUN
    label: ClassVar[str] = "grounded_ergun_bed_resistance"
    carries_particle_length_scale: ClassVar[bool] = False
    #: TRUE for every constructible instance, because the only scale factor on the
    #: correlation is the two-valued nested-limit switch enforced below. A
    #: real-valued scale would make this flag unable to report its own instance,
    #: which is why the domain and not the flag is what carries the claim.
    coefficients_are_non_fitted: ClassVar[bool] = True
    resistance_scale_is_a_nested_limit_switch: ClassVar[bool] = True
    #: The family this element is admitted on. The bed Ergun resistance is the
    #: family-2 archetype's first ruled element and belongs to no other family.
    admitted_families: ClassVar[tuple[TrayFloorFamily, ...]] = (
        TrayFloorFamily.STRIPPING_DECK,
    )

    def __post_init__(self) -> None:
        _require_positive("hydraulic_granule_diameter_m", self.hydraulic_granule_diameter_m)
        _require_fraction("bed_void_fraction", self.bed_void_fraction)
        if self.bed_void_fraction <= 0.0:
            raise ClosureConfigurationError("bed_void_fraction must be strictly positive")
        _require_positive("flow_length_m", self.flow_length_m)
        _require_positive("gas_dynamic_viscosity_pa_s", self.gas_dynamic_viscosity_pa_s)
        _require_nested_limit_switch("resistance_scale", self.resistance_scale)
        if self.hydraulic_granule_diameter_m <= 2.0 * INTRAPARTICLE_DIFFUSION_RADIUS_M:
            raise LengthScaleConfusionError(
                f"the hydraulic diameter {self.hydraulic_granule_diameter_m!r} m is at or "
                f"below the frozen intraparticle diffusion scale "
                f"({INTRAPARTICLE_DIFFUSION_RADIUS_M!r} m radius); FROZEN PHY-017 requires "
                "the true external granule dimension and the two length scales are "
                "distinct by construction"
            )

    def contribution(
        self, *, superficial_velocity_m_s: float, gas_density_kg_m3: float
    ) -> PressureLossContribution:
        void = self.bed_void_fraction
        diameter = self.hydraulic_granule_diameter_m
        solid = 1.0 - void
        void_cubed = void * void * void
        viscous = (
            ERGUN_VISCOUS_COEFFICIENT
            * self.gas_dynamic_viscosity_pa_s
            * solid
            * solid
            / (void_cubed * diameter * diameter)
        )
        inertial = ERGUN_INERTIAL_COEFFICIENT * solid / (void_cubed * diameter)
        span = self.resistance_scale * self.flow_length_m
        velocity = superficial_velocity_m_s
        return PressureLossContribution(
            drop_pa=span
            * (viscous * velocity + inertial * gas_density_kg_m3 * velocity * velocity),
            d_drop_d_superficial_velocity=span
            * (viscous + 2.0 * inertial * gas_density_kg_m3 * velocity),
            d_drop_d_gas_density=span * inertial * velocity * velocity,
            kind=self.kind,
            element_label=self.label,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class DerivedWideEndNozzleInput:
    """The wide-end area fraction AS A DERIVED, PROVENANCE-BOUND NOZZLE INPUT.

    Owner ruling item 4 (2026-08-12): "Retain 25 mm as the printed bore maximum.
    Any wide-end value derived from the printed area ratio shall be stored only as
    a derived, provenance-bound nozzle input." This type is that binding: it is
    constructed on read from the throat and the printed ratio, it is never a
    stored field of anything, and it carries the printed band and the anchor the
    ratio belongs to so the number cannot travel without its source.

    It carries no loss arithmetic. The taper ratio constrains geometry only; the
    loss-relevant dimension is the throat.
    """

    wide_end_area_fraction: float
    throat_area_fraction: float
    taper_area_ratio: float
    taper_ratio_band: tuple[float, float]
    provenance: DeclaredGeometryProvenance

    is_derived_never_stored: ClassVar[bool] = True
    enters_the_loss_arithmetic: ClassVar[bool] = False
    provenance_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION

    def __post_init__(self) -> None:
        _require_fraction("wide_end_area_fraction", self.wide_end_area_fraction)
        _require_fraction("throat_area_fraction", self.throat_area_fraction)
        _require_positive("taper_area_ratio", self.taper_area_ratio)
        if type(self.provenance) is not DeclaredGeometryProvenance:
            raise ClosureConfigurationError(
                "a derived wide end must name the DeclaredGeometryProvenance whose "
                "printed area ratio it was derived from"
            )
        low, high = self.taper_ratio_band
        if not low <= self.taper_area_ratio <= high:
            raise ClosureConfigurationError(
                "a derived wide end must come from a ratio inside its printed band"
            )

    @property
    def evidentiary_standing_rank(self) -> int:
        """The ORDINAL standing of the anchor this number is bound to."""

        return self.provenance.evidentiary_standing_rank


@dataclass(frozen=True, slots=True, kw_only=True)
class TaperedBoreFloorPassageElement:
    """The ruled tapered-bore floor passage, IN SERIES with the bed.

    Owner ruling 2026-08-12. The arithmetic is a throat-referenced passage loss,
    ``dP = rho u_throat^2 / (2 Cd^2)`` with ``u_throat = u_superficial / phi``,
    where ``phi`` is the declared THROAT area fraction of the deck.

    THE PRINTED TAPER RATIO CONSTRAINS GEOMETRY ONLY. ``taper_area_ratio`` is
    checked against the printed band and relates the wide end of the bore to the
    throat; it does not appear in the drop, because the loss-relevant dimension is
    the throat. ``discharge_coefficient`` is the separate declared engineering
    assumption, checked against its predeclared sensitivity range, which is a FIELD
    so that it travels with the element and can be swept.
    """

    throat_area_fraction: float
    taper_area_ratio: float
    discharge_coefficient: float
    discharge_coefficient_sensitivity_range: tuple[float, float] = (
        FLOOR_DISCHARGE_COEFFICIENT_SENSITIVITY_RANGE
    )
    resistance_scale: float = 1.0

    kind: ClassVar[PressureLossRelation] = PressureLossRelation.SPARGE_NOZZLE
    label: ClassVar[str] = "ruled_tapered_bore_floor_passage"
    #: The coefficient is declared, never fitted, and never inferred jointly with
    #: the bed resistance from an aggregate drop.
    discharge_coefficient_provenance: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION
    taper_ratio_constrains_geometry_only: ClassVar[bool] = True
    carries_particle_length_scale: ClassVar[bool] = False
    resistance_scale_is_a_nested_limit_switch: ClassVar[bool] = True
    #: THE ANCHOR THE PRINTED TAPER RATIO COMES FROM, bound rather than implied.
    #: ``TAPER_AREA_RATIO_BAND`` is Schumacher US4619053A's printed f:g ratio, so any
    #: value derived from it is bound to that anchor and to no other -- which is also
    #: why this element is admitted on family 2 alone.
    taper_ratio_provenance: ClassVar[DeclaredGeometryProvenance] = (
        DeclaredGeometryProvenance.SCHUMACHER_US4619053A
    )
    admitted_families: ClassVar[tuple[TrayFloorFamily, ...]] = (
        TrayFloorFamily.STRIPPING_DECK,
    )

    def __post_init__(self) -> None:
        _require_fraction("throat_area_fraction", self.throat_area_fraction)
        if self.throat_area_fraction <= 0.0:
            raise ClosureConfigurationError(
                "a floor passage must have a strictly positive throat area fraction"
            )
        _require_positive("taper_area_ratio", self.taper_area_ratio)
        low, high = TAPER_AREA_RATIO_BAND
        if not low <= self.taper_area_ratio <= high:
            raise ClosureConfigurationError(
                f"the taper area ratio {self.taper_area_ratio!r} is outside the printed "
                f"band [{low}, {high}] -- refused, never clipped and never averaged"
            )
        if (
            type(self.discharge_coefficient_sensitivity_range) is not tuple
            or len(self.discharge_coefficient_sensitivity_range) != 2
        ):
            raise ClosureConfigurationError(
                "the discharge-coefficient sensitivity range must be a 2-tuple"
            )
        range_low, range_high = self.discharge_coefficient_sensitivity_range
        _require_positive("sensitivity range low", range_low)
        _require_positive("sensitivity range high", range_high)
        if not range_low <= range_high:
            raise ClosureConfigurationError("the sensitivity range edges are out of order")
        declared_low, declared_high = FLOOR_DISCHARGE_COEFFICIENT_SENSITIVITY_RANGE
        if range_low < declared_low or range_high > declared_high:
            raise ClosureConfigurationError(
                f"the sensitivity range [{range_low}, {range_high}] widens the predeclared "
                f"range [{declared_low}, {declared_high}]; narrowing it with evidence is "
                "admitted, widening it is an owner ruling"
            )
        _require_positive("discharge_coefficient", self.discharge_coefficient)
        if not range_low <= self.discharge_coefficient <= range_high:
            raise ClosureConfigurationError(
                f"the declared discharge coefficient {self.discharge_coefficient!r} is "
                f"outside its own predeclared sensitivity range [{range_low}, {range_high}]"
            )
        _require_nested_limit_switch("resistance_scale", self.resistance_scale)
        # THE DERIVED WIDE END MUST BE A POSSIBLE DECK. throat * ratio is an area
        # fraction of the same deck face, so at or above 1.0 it is a deck more than
        # fully open -- impossible geometry, refused here exactly as the inherited
        # module refuses its own derived open fraction. Refused, never clipped.
        wide_end = self.throat_area_fraction * self.taper_area_ratio
        if wide_end >= 1.0:
            raise ClosureConfigurationError(
                f"the throat fraction {self.throat_area_fraction!r} at the printed "
                f"taper ratio {self.taper_area_ratio!r} implies a wide-end area "
                f"fraction of {wide_end!r} -- impossible geometry, refused"
            )

    @property
    def wide_end_area_fraction(self) -> float:
        """GEOMETRY ONLY: the wide end implied by the printed taper ratio.

        DERIVED on read and stored nowhere, per owner ruling item 4. Bounded below
        1.0 at construction, and bound to the anchor whose printed ratio it comes
        from by ``derived_wide_end_nozzle_input``, which is the accessor a nozzle
        consumer should take.
        """

        return self.throat_area_fraction * self.taper_area_ratio

    @property
    def derived_wide_end_nozzle_input(self) -> DerivedWideEndNozzleInput:
        """The wide end AS THE RULED NOZZLE INPUT: derived, and provenance-bound.

        Owner ruling item 4: "Any wide-end value derived from the printed area
        ratio shall be stored only as a derived, provenance-bound nozzle input."
        So it is computed on read, never a field, and it travels with the anchor
        whose printed ratio produced it -- a consumer cannot receive the number
        without receiving which printed record it came from.
        """

        return DerivedWideEndNozzleInput(
            wide_end_area_fraction=self.wide_end_area_fraction,
            throat_area_fraction=self.throat_area_fraction,
            taper_area_ratio=self.taper_area_ratio,
            taper_ratio_band=TAPER_AREA_RATIO_BAND,
            provenance=self.taper_ratio_provenance,
        )

    def contribution(
        self, *, superficial_velocity_m_s: float, gas_density_kg_m3: float
    ) -> PressureLossContribution:
        fraction = self.throat_area_fraction
        coefficient = self.discharge_coefficient
        # dP = scale * rho * (u/phi)^2 / (2 Cd^2). The taper ratio is absent by
        # ruling: the loss-relevant dimension is the throat.
        factor = self.resistance_scale / (
            2.0 * coefficient * coefficient * fraction * fraction
        )
        velocity = superficial_velocity_m_s
        return PressureLossContribution(
            drop_pa=factor * gas_density_kg_m3 * velocity * velocity,
            d_drop_d_superficial_velocity=2.0 * factor * gas_density_kg_m3 * velocity,
            d_drop_d_gas_density=factor * velocity * velocity,
            kind=self.kind,
            element_label=self.label,
        )


def sweep_floor_discharge_coefficient(
    element: TaperedBoreFloorPassageElement, *, points: int = 3
) -> tuple[TaperedBoreFloorPassageElement, ...]:
    """The predeclared sensitivity sweep, as elements ready to be evaluated.

    Endpoints included, so the declared range's edges are always exercised. This
    is a SWEEP, not a fit: it produces candidate declarations to run the closure
    at, and reads no measurement.
    """

    if type(element) is not TaperedBoreFloorPassageElement:
        raise ClosureConfigurationError("element must be a TaperedBoreFloorPassageElement")
    if type(points) is not int or points < 2:
        raise ClosureConfigurationError("a sweep needs at least the two range edges")
    low, high = element.discharge_coefficient_sensitivity_range
    span = high - low
    swept: list[TaperedBoreFloorPassageElement] = []
    for index in range(points):
        value = low + span * index / (points - 1)
        swept.append(replace(element, discharge_coefficient=value))
    return tuple(swept)


_IDENTIFIABILITY_REFUSAL_MESSAGE = (
    "the bed resistance and the floor discharge coefficient may not be jointly "
    "inferred from a total pressure drop: two coefficients are not identifiable "
    "from one aggregate. The floor coefficient is a declared engineering "
    "assumption with a predeclared sensitivity range, swept by "
    "sweep_floor_discharge_coefficient, and every series evaluation returns its "
    "per-element contributions so the total is never the only observable."
)


def infer_series_coefficients_from_total_drop(*_args: object, **_kwargs: object) -> None:
    """ALWAYS REFUSES. The entry point that names the forbidden operation.

    Owner ruling 2026-08-12: "Bed and floor coefficients shall not be jointly
    inferred from total pressure drop alone." Two coefficients are not identifiable
    from one aggregate, so this module offers no inference path and this function
    exists to make the refusal executable rather than editorial.

    It is deliberately NOT the only guard: a caller reaching for calibration
    reaches for ``fit_series_coefficients_to_observation``, where the refusal is
    selected by the TYPE OF THE TARGET rather than by the name of the function.
    """

    raise JointIdentifiabilityRefusedError(_IDENTIFIABILITY_REFUSAL_MESSAGE)


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerDropObservation:
    """A measured or targeted layer drop, WITH ITS DECOMPOSITION OR WITHOUT IT.

    This type exists so that the ruled prohibition can be enforced on the SHAPE OF
    THE TARGET instead of on the name of a function. ``per_element_drops_pa`` is
    ``None`` for a total-only observation -- one aggregate number, which is exactly
    the target the owner ruling forbids fitting two coefficients against -- and a
    mapping from element label to that element's own measured drop otherwise.

    The distinction is structural, not conventional: an empty mapping is refused,
    because "I measured nothing per element" is a total-only observation wearing a
    mapping's clothes.
    """

    drop_pa: float
    per_element_drops_pa: Mapping[str, float] | None = None

    def __post_init__(self) -> None:
        _require_nonnegative("drop_pa", self.drop_pa)
        if self.per_element_drops_pa is None:
            return
        if not isinstance(self.per_element_drops_pa, Mapping):
            raise ClosureConfigurationError(
                "per_element_drops_pa must be None or a mapping of element label to drop"
            )
        if not self.per_element_drops_pa:
            raise ClosureConfigurationError(
                "an empty per-element mapping is a total-only observation; declare it "
                "as None so the identifiability refusal can see it"
            )
        for label, value in self.per_element_drops_pa.items():
            if type(label) is not str or not label.strip():
                raise ClosureConfigurationError("every element label must be a nonblank string")
            _require_nonnegative(f"per-element drop {label}", value)

    @property
    def is_total_only(self) -> bool:
        """TRUE when the aggregate is the only observable in this observation."""

        return self.per_element_drops_pa is None


def fit_series_coefficients_to_observation(
    series: LayerLossSeries, observation: LayerDropObservation
) -> None:
    """ALWAYS REFUSES, AND THE TARGET'S OWN SHAPE CHOOSES THE REFUSAL.

    This is the entry point a caller reaching for calibration arrives at, so the
    ruled prohibition is routed rather than decorative:

    * a TOTAL-ONLY observation raises ``JointIdentifiabilityRefusedError`` -- the
      ruled refusal, selected because the target carries one aggregate and two
      coefficients are not identifiable from it;
    * an observation carrying PER-ELEMENT drops raises
      ``CoefficientInferenceNotImplementedError``, because this module implements no
      fitting, calibration or inference path at all. Per-element data would not be
      the forbidden operation, but the admitted operation on the declared
      coefficient is ``sweep_floor_discharge_coefficient`` over its predeclared
      range, and no fitter exists here to reach for.

    Either way nothing is returned and no coefficient is ever determined from a
    measurement inside this module.
    """

    if type(series) is not LayerLossSeries:
        raise ClosureConfigurationError("series must be an exact LayerLossSeries")
    if type(observation) is not LayerDropObservation:
        raise ClosureConfigurationError(
            "observation must be an exact LayerDropObservation, so the target's shape "
            "is typed rather than inferred"
        )
    if observation.is_total_only:
        raise JointIdentifiabilityRefusedError(_IDENTIFIABILITY_REFUSAL_MESSAGE)
    raise CoefficientInferenceNotImplementedError(
        "this module implements no fitting, calibration or inference path, not even "
        "against a decomposed observation: the floor discharge coefficient is a "
        "declared engineering assumption and the admitted operation on it is the "
        "predeclared sweep, sweep_floor_discharge_coefficient"
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerLossSeries:
    """An ORDERED series of independently typed loss elements (ruled 2026-08-12).

    The order is the ruled order and is checked against
    ``RULED_ELEMENT_KIND_SERIES`` for the declared family, so a series labelled
    family 2 cannot ship the bed alone and cannot ship the floor first.

    THE KIND TAG IS NOT ENOUGH, so the element's own admitted family is checked
    too. The inherited nozzle member tags the family-2 tapered bore here, and a
    kind-only check would let that family-2 element -- whose constructor enforces a
    SCHUMACHER printed taper band -- satisfy family 3's ruled single-nozzle series.
    Each element therefore declares ``admitted_families`` and a series refuses a
    member that is not admitted on the family it is labelled with.

    Note the scope of the label: this type checks a series against ITS OWN declared
    family. Agreeing with the DECK's family is ``CellClosureInputs``'s check, and it
    performs it.
    """

    family: TrayFloorFamily
    elements: tuple[object, ...]

    ruling: ClassVar[str] = PRESSURE_SERIES_RULING

    def __post_init__(self) -> None:
        if type(self.family) is not TrayFloorFamily:
            raise ClosureConfigurationError("family must be a TrayFloorFamily member")
        if type(self.elements) is not tuple or not self.elements:
            raise LossSeriesError("a loss series carries at least one typed element")
        expected = RULED_ELEMENT_KIND_SERIES[self.family]
        kinds = []
        for element in self.elements:
            kind = getattr(type(element), "kind", None)
            if type(kind) is not PressureLossRelation or not callable(
                getattr(element, "contribution", None)
            ):
                raise LossSeriesError(
                    "every series member must be a typed loss element with a kind and "
                    "its own contribution"
                )
            admitted = getattr(type(element), "admitted_families", None)
            if type(admitted) is not tuple or not admitted:
                raise LossSeriesError(
                    "every series member must declare the families it is admitted on"
                )
            if self.family not in admitted:
                raise LossSeriesError(
                    f"{type(element).__name__} is admitted on "
                    f"{tuple(item.value for item in admitted)} and not on family "
                    f"{self.family.d7_family_number}; the shared kind TAG does not make "
                    "one family's printed geometry admissible on another's deck"
                )
            kinds.append(kind)
        if tuple(kinds) != expected:
            raise LossSeriesError(
                f"family {self.family.d7_family_number} is ruled to carry the ordered "
                f"element kinds {tuple(k.value for k in expected)}, not "
                f"{tuple(k.value for k in kinds)}"
            )

    @property
    def element_kinds(self) -> tuple[PressureLossRelation, ...]:
        return tuple(type(element).kind for element in self.elements)

    @property
    def bed_element(self) -> BedErgunElement | None:
        for element in self.elements:
            if type(element) is BedErgunElement:
                return element
        return None

    @property
    def is_unanchored_by_reference_machine(self) -> bool:
        return self.family in UNANCHORED_BY_REFERENCE_MACHINE_FAMILIES


def build_family_loss_series(
    family: TrayFloorFamily, *, elements: tuple[object, ...]
) -> LayerLossSeries:
    """Construct a family's ruled series, refusing families with no element type."""

    if type(family) is not TrayFloorFamily:
        raise ClosureConfigurationError("family must be a TrayFloorFamily member")
    if family is TrayFloorFamily.DIRECT_STEAM_SPARGE:
        raise LossRelationNotImplementedError(
            "family 3 has no loss element of its own here. Its ruled series is a "
            "single nozzle kind, and the only nozzle element written in this packet "
            "is family 2's tapered bore, whose constructor enforces the SCHUMACHER "
            "US4619053A printed taper band -- so building family 3's series from it "
            "would import a family-2 anchor's printed geometry onto a Desmet deck. "
            "Refused rather than blended: a dedicated family-3 element, grounded on "
            "the printed Desmet sparge passages, is owed work and is not written here"
        )
    for kind in RULED_ELEMENT_KIND_SERIES[family]:
        if kind is PressureLossRelation.DOME_BYPASS:
            raise LossRelationNotImplementedError(
                "the bypass path has no printed per-layer drop inside the corroborated "
                "window, so no element type is written for it; a zero would silently "
                "assert a negligible drop, which is a nested limit of the frozen network "
                "and not its default"
            )
        if kind is PressureLossRelation.AIR_PLENUM:
            raise LossRelationNotImplementedError(
                "the air-plenum kind has no element type here: its family has no "
                "reference-machine member and no admitted dimensioned anchor"
            )
    return LayerLossSeries(family=family, elements=elements)


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerHydraulics:
    """The layer's ordered loss series and its network mode.

    This is the only type the per-layer pressure residual reads, and every length
    scale it can reach belongs to an element. No field of this type or of any
    element names a particle dimension.
    """

    series: LayerLossSeries
    network_mode: PressureNetworkMode = PressureNetworkMode.NODAL_HYDRAULIC
    #: The prescribed-profile limit: a frozen resistance whose derivatives are
    #: exactly zero. Required under FROZEN_RESISTANCE, forbidden otherwise.
    frozen_resistance_drop_pa: float | None = None

    carries_particle_length_scale: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.series) is not LayerLossSeries:
            raise ClosureConfigurationError("series must be an exact LayerLossSeries")
        if type(self.network_mode) is not PressureNetworkMode:
            raise ClosureConfigurationError("network_mode must be a PressureNetworkMode member")
        if self.network_mode is PressureNetworkMode.FROZEN_RESISTANCE:
            if self.frozen_resistance_drop_pa is None:
                raise ClosureConfigurationError(
                    "the frozen-resistance limit needs its prescribed drop"
                )
            _require_nonnegative("frozen_resistance_drop_pa", self.frozen_resistance_drop_pa)
        elif self.frozen_resistance_drop_pa is not None:
            raise ClosureConfigurationError(
                "a prescribed drop is meaningful only under the frozen-resistance limit"
            )

    @property
    def bed_void_fraction(self) -> float | None:
        bed = self.series.bed_element
        return None if bed is None else bed.bed_void_fraction


@dataclass(frozen=True, slots=True, kw_only=True)
class FilmAreaLaw:
    """FROZEN PHY-023 / PHY-032 area and inventory law. ONE shared area scalar.

    The area itself is the inherited dispatch area times (1 - S_L)^n; the void
    fraction and the two liquid densities are here because the frozen saturation
    definitions convert an external liquid mass to a saturation through exactly
    those quantities.
    """

    exponent: FilmAreaExponent = FilmAreaExponent.ISOTROPIC_REFERENCE
    bed_void_fraction: float
    hexane_liquid_density_kg_m3: float
    water_liquid_density_kg_m3: float

    #: FROZEN PHY-023 ``:916`` / PHY-032: one common area, never a per-component
    #: pair summed together.
    per_component_area_admitted: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.exponent) is not FilmAreaExponent:
            raise ClosureConfigurationError("exponent must be a FilmAreaExponent member")
        _require_fraction("bed_void_fraction", self.bed_void_fraction)
        if self.bed_void_fraction <= 0.0:
            raise ClosureConfigurationError("bed_void_fraction must be strictly positive")
        _require_positive("hexane_liquid_density_kg_m3", self.hexane_liquid_density_kg_m3)
        _require_positive("water_liquid_density_kg_m3", self.water_liquid_density_kg_m3)

    def hexane_inventory_capacity_kg(self, bed_volume_m3: float) -> float:
        """Kilograms of external hexane per unit saturation, from the frozen
        definition S_f = X_f (1-eps) rho_dm / (eps rho_h,l)."""

        return self.bed_void_fraction * self.hexane_liquid_density_kg_m3 * bed_volume_m3

    def water_inventory_capacity_kg(self, bed_volume_m3: float) -> float:
        return self.bed_void_fraction * self.water_liquid_density_kg_m3 * bed_volume_m3


@dataclass(frozen=True, slots=True, kw_only=True)
class DeclaredGasProperties:
    """Constant-property gas closure. EVERY value is a declared caller input.

    This module mints no property number and consumes no frozen property module,
    so it authenticates no property source. One molar heat capacity carries the
    whole gas mixture and one latent heat carries each species; both are declared
    simplifications, and the datum is explicit so the energy rows are checkable.
    """

    molar_heat_capacity_j_mol_k: float
    hexane_molar_mass_kg_mol: float
    water_molar_mass_kg_mol: float
    inert_molar_mass_kg_mol: float
    hexane_latent_heat_j_mol: float
    water_latent_heat_j_mol: float
    energy_datum_temperature_k: float
    energy_datum_id: str

    provenance_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION

    def __post_init__(self) -> None:
        _require_positive("molar_heat_capacity_j_mol_k", self.molar_heat_capacity_j_mol_k)
        _require_positive("hexane_molar_mass_kg_mol", self.hexane_molar_mass_kg_mol)
        _require_positive("water_molar_mass_kg_mol", self.water_molar_mass_kg_mol)
        _require_positive("inert_molar_mass_kg_mol", self.inert_molar_mass_kg_mol)
        _require_positive("hexane_latent_heat_j_mol", self.hexane_latent_heat_j_mol)
        _require_positive("water_latent_heat_j_mol", self.water_latent_heat_j_mol)
        _require_positive("energy_datum_temperature_k", self.energy_datum_temperature_k)
        if type(self.energy_datum_id) is not str or not self.energy_datum_id.strip():
            raise ClosureConfigurationError("energy_datum_id must be a nonblank string")
        if self.molar_heat_capacity_j_mol_k <= UNIVERSAL_GAS_CONSTANT_J_MOL_K:
            raise ClosureConfigurationError(
                "a molar heat capacity at or below R leaves a non-positive c_v"
            )

    @property
    def molar_heat_capacity_constant_volume_j_mol_k(self) -> float:
        """c_v = c_p - R, used only by the dynamic-EOS internal-energy storage."""

        return self.molar_heat_capacity_j_mol_k - UNIVERSAL_GAS_CONSTANT_J_MOL_K


@dataclass(frozen=True, slots=True, kw_only=True)
class ClausiusClapeyronSaturation:
    """A declared saturation law p_sat(T) = p_ref exp(k (T - T_ref)).

    Derived form, zero fitted structure, and the caller declares all three
    numbers: this module refuses to mint a reference point or a slope, so no
    property value originates here.
    """

    reference_temperature_k: float
    reference_pressure_pa: float
    slope_per_k: float

    provenance_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION

    def __post_init__(self) -> None:
        _require_positive("reference_temperature_k", self.reference_temperature_k)
        _require_positive("reference_pressure_pa", self.reference_pressure_pa)
        _require_positive("slope_per_k", self.slope_per_k)

    def saturation_pressure_pa(self, temperature_k: float) -> float:
        """The declared law, with its own overflow typed rather than escaping.

        A public scalar accessor, so a caller can hand it any temperature. The
        exponential raises rather than returning an infinity for a large finite
        argument, and an untyped ``OverflowError`` out of this module would break
        the two-tier refusal contract; it is re-raised as a step rejection.
        """

        try:
            return self.reference_pressure_pa * math.exp(
                self.slope_per_k * (temperature_k - self.reference_temperature_k)
            )
        except OverflowError as error:
            raise StepRejected(
                f"the declared saturation law is not representable at "
                f"{temperature_k!r} K"
            ) from error

    def _dual(self, temperature: _Dual) -> _Dual:
        exponent = (temperature - self.reference_temperature_k) * self.slope_per_k
        return _exp(exponent) * self.reference_pressure_pa


@dataclass(frozen=True, slots=True, kw_only=True)
class InterfaceEquilibrium:
    """FROZEN PHY-007: two IMMISCIBLE PURE liquid phases at the layer pressure.

    Each present liquid sets its own interfacial partial pressure from its own
    pure-component saturation law; they share one gas film, one gas-side area and
    one coupled component-flux vector. That is exactly the frozen two-liquid
    branch, and it is why no mixture activity model appears here: the frozen
    liquid-phase rule keeps water and hexane distinct pure phases, and the shared
    external holdup coordinate does not imply liquid mixing.
    """

    hexane: ClausiusClapeyronSaturation
    water: ClausiusClapeyronSaturation

    def __post_init__(self) -> None:
        if type(self.hexane) is not ClausiusClapeyronSaturation:
            raise ClosureConfigurationError("hexane law must be a ClausiusClapeyronSaturation")
        if type(self.water) is not ClausiusClapeyronSaturation:
            raise ClosureConfigurationError("water law must be a ClausiusClapeyronSaturation")


@dataclass(frozen=True, slots=True, kw_only=True)
class TransferCoefficients:
    """Declared transfer conductances. The three wall couplings live here too.

    ``solid_to_interface_ua_w_k`` is the meal-side path; ``wall_to_interface_ua_w_k``
    is FROZEN PHY-008's tray-surface conduction into the meal above; and
    ``wall_to_gas_ua_w_k`` is its convection to the gas below. The steam side and
    the ambient loss belong to the wall node itself.
    """

    gas_side_molar_conductance_mol_m2_s: float
    gas_side_heat_coefficient_w_m2_k: float
    solid_to_interface_ua_w_k: float
    wall_to_interface_ua_w_k: float
    wall_to_gas_ua_w_k: float

    provenance_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION

    def __post_init__(self) -> None:
        _require_positive(
            "gas_side_molar_conductance_mol_m2_s", self.gas_side_molar_conductance_mol_m2_s
        )
        _require_positive(
            "gas_side_heat_coefficient_w_m2_k", self.gas_side_heat_coefficient_w_m2_k
        )
        _require_nonnegative("solid_to_interface_ua_w_k", self.solid_to_interface_ua_w_k)
        _require_nonnegative("wall_to_interface_ua_w_k", self.wall_to_interface_ua_w_k)
        _require_nonnegative("wall_to_gas_ua_w_k", self.wall_to_gas_ua_w_k)


@dataclass(frozen=True, slots=True, kw_only=True)
class SolidSideAggregates:
    """PACKET-DERIVED aggregates. READ, never re-derived.

    The two external inventories are the frozen active-set conditions' own
    arguments. The cell layer owns neither solid holdup, cohort identity, RTD
    stage nor exit logic, and it does not recompute these loadings from anything.
    """

    attached_hexane_loading_kg_kg: float
    external_water_loading_kg_kg: float
    dry_matter_holdup_kg: float
    solid_temperature_k: float
    hexane_supply_to_film_kg_s: float
    water_supply_to_film_kg_s: float

    inventories_are_read_never_rederived: ClassVar[bool] = True

    def __post_init__(self) -> None:
        _require_nonnegative(
            "attached_hexane_loading_kg_kg", self.attached_hexane_loading_kg_kg
        )
        _require_nonnegative(
            "external_water_loading_kg_kg", self.external_water_loading_kg_kg
        )
        _require_positive("dry_matter_holdup_kg", self.dry_matter_holdup_kg)
        _require_positive("solid_temperature_k", self.solid_temperature_k)
        _require_binary64("hexane_supply_to_film_kg_s", self.hexane_supply_to_film_kg_s)
        _require_binary64("water_supply_to_film_kg_s", self.water_supply_to_film_kg_s)

    @property
    def interface_branch(self) -> InterfaceBranch:
        """FROZEN PHY-007's active set, by EXACT comparison on the inventories.

        The comparisons are exact: no epsilon shell, no tolerance and no smoothed
        transition. A loading of 1e-300 is present, because presence is what the
        frozen condition tests and endpoint identity arrives from the event side.
        """

        hexane_present = self.attached_hexane_loading_kg_kg > 0.0
        water_present = self.external_water_loading_kg_kg > 0.0
        if hexane_present and water_present:
            return InterfaceBranch.TWO_EXTERNAL_LIQUIDS
        if hexane_present:
            return InterfaceBranch.EXTERNAL_HEXANE
        if water_present:
            return InterfaceBranch.EXTERNAL_WATER
        return InterfaceBranch.DRY_SURFACE


@dataclass(frozen=True, slots=True, kw_only=True)
class GasBoundary:
    """The cell's gas inlet and its downstream pressure boundary."""

    inlet_molar_flow_mol_s: float
    inlet_hexane_mole_fraction: float
    inlet_water_mole_fraction: float
    inlet_temperature_k: float
    downstream_boundary_pressure_pa: float

    def __post_init__(self) -> None:
        _require_positive("inlet_molar_flow_mol_s", self.inlet_molar_flow_mol_s)
        _require_fraction("inlet_hexane_mole_fraction", self.inlet_hexane_mole_fraction)
        _require_fraction("inlet_water_mole_fraction", self.inlet_water_mole_fraction)
        if self.inlet_hexane_mole_fraction + self.inlet_water_mole_fraction >= 1.0:
            raise ClosureConfigurationError(
                "the inlet gas must retain an inert carrier fraction"
            )
        _require_positive("inlet_temperature_k", self.inlet_temperature_k)
        _require_positive(
            "downstream_boundary_pressure_pa", self.downstream_boundary_pressure_pa
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class WallNodeParameters:
    """FROZEN PHY-008's dynamic metal node, plus FROZEN PHY-026's ambient UA.

    ``ambient_ua_w_k = 0.0`` is reachable EXACTLY, which is the mandatory
    adiabatic limiting oracle acting on this node (``:1616-1618``).

    ``thermal_capacity_j_k`` is a DECLARED_ENGINEERING_ASSUMPTION whose PREDECLARED
    RANGE is not on the capacity itself but on what the capacity is for: the frozen
    ratio bracket ``PHY008_WALL_RESIDENCE_RATIO_BRACKET``, which at a given
    conductance sum and macro step pins the admissible capacity. The bracket is
    stated here and tested by ``is_slowest_mode`` rather than enforced in this
    constructor, because the quasi-steady limit and the zero-conductance refusal are
    both legitimate parameterisations outside it.
    """

    thermal_capacity_j_k: float
    steam_side_ua_w_k: float
    steam_temperature_k: float
    ambient_ua_w_k: float
    ambient_temperature_k: float
    closure: WallClosure = WallClosure.DYNAMIC_NODE

    provenance_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION
    #: The predeclared range of the ratio the capacity exists to produce.
    time_constant_ratio_predeclared_range: ClassVar[tuple[float, float]] = (
        PHY008_WALL_RESIDENCE_RATIO_BRACKET
    )

    def __post_init__(self) -> None:
        _require_positive("thermal_capacity_j_k", self.thermal_capacity_j_k)
        _require_nonnegative("steam_side_ua_w_k", self.steam_side_ua_w_k)
        _require_positive("steam_temperature_k", self.steam_temperature_k)
        _require_nonnegative("ambient_ua_w_k", self.ambient_ua_w_k)
        _require_positive("ambient_temperature_k", self.ambient_temperature_k)
        if type(self.closure) is not WallClosure:
            raise ClosureConfigurationError("closure must be a WallClosure member")


@dataclass(frozen=True, slots=True, kw_only=True)
class CellFieldState:
    """The fast block's unknowns as a named state. PER CELL, never global."""

    gas_temperature_k: float
    gas_hexane_mole_fraction: float
    gas_water_mole_fraction: float
    layer_pressure_pa: float
    interface_temperature_k: float
    external_hexane_saturation: float
    external_water_saturation: float

    def __post_init__(self) -> None:
        for name in UNKNOWN_NAMES:
            _require_binary64(name, getattr(self, name))

    def as_vector(self) -> tuple[float, ...]:
        return tuple(getattr(self, name) for name in UNKNOWN_NAMES)

    @classmethod
    def from_vector(cls, vector: tuple[float, ...]) -> CellFieldState:
        if len(vector) != _N_UNKNOWNS:
            raise ClosureConfigurationError(
                f"a cell field state carries {_N_UNKNOWNS} unknowns"
            )
        return cls(**dict(zip(UNKNOWN_NAMES, vector, strict=True)))


# --- the ordered series composition ----------------------------------------
# Each element contributes its own analytic derivatives and the series SUMS them.
# A relation reaches the residual through exactly two declared channels and
# nothing else, so admitting a further element moves numbers and moves neither
# the residual's shape nor the Jacobian's sparsity.


@dataclass(frozen=True, slots=True, kw_only=True)
class PressureLossDrop:
    """The composed layer drop, with its PER-ELEMENT contributions retained.

    The decomposition is retained deliberately: the identifiability prohibition
    means an aggregate total must never be the only observable, so every series
    evaluation hands back what each element contributed.
    """

    drop_pa: float
    d_drop_d_superficial_velocity: float
    d_drop_d_gas_density: float
    contributions: tuple[PressureLossContribution, ...]
    composition: str

    #: THE ONLY CHANNELS. The chain rule to the unknowns is applied once, outside
    #: the elements, over exactly these two derived scalars.
    channels: ClassVar[tuple[str, ...]] = ("superficial_velocity_m_s", "gas_density_kg_m3")
    ordered_series_is_ruled: ClassVar[bool] = True
    aggregate_is_never_the_only_observable: ClassVar[bool] = True

    def __post_init__(self) -> None:
        _require_binary64("drop_pa", self.drop_pa)
        _require_binary64(
            "d_drop_d_superficial_velocity", self.d_drop_d_superficial_velocity
        )
        _require_binary64("d_drop_d_gas_density", self.d_drop_d_gas_density)
        if type(self.contributions) is not tuple:
            raise ClosureConfigurationError("contributions must be a tuple")
        for item in self.contributions:
            if type(item) is not PressureLossContribution:
                raise ClosureConfigurationError(
                    "every contribution must be an exact PressureLossContribution"
                )
        if type(self.composition) is not str or not self.composition.strip():
            raise ClosureConfigurationError("composition must be a nonblank string")

    @property
    def element_kinds(self) -> tuple[PressureLossRelation, ...]:
        return tuple(item.kind for item in self.contributions)

    def contribution_for(self, label: str) -> PressureLossContribution:
        for item in self.contributions:
            if item.element_label == label:
                return item
        raise LossSeriesError(f"no element labelled {label!r} contributed to this drop")


def layer_pressure_drop(
    *,
    series: LayerLossSeries,
    superficial_velocity_m_s: float,
    gas_density_kg_m3: float,
) -> PressureLossDrop:
    """Compose the layer drop over the ORDERED series of typed elements.

    Each element is asked for its own contribution in the ruled order; the drops
    and both channel derivatives are summed. The element kinds are tags: this
    function never interprets one, it only walks the series it was handed.
    """

    if type(series) is not LayerLossSeries:
        raise ClosureConfigurationError("series must be an exact LayerLossSeries")
    _require_nonnegative("superficial_velocity_m_s", superficial_velocity_m_s)
    _require_positive("gas_density_kg_m3", gas_density_kg_m3)

    contributions: list[PressureLossContribution] = []
    total = 0.0
    d_velocity = 0.0
    d_density = 0.0
    for element in series.elements:
        item = element.contribution(
            superficial_velocity_m_s=superficial_velocity_m_s,
            gas_density_kg_m3=gas_density_kg_m3,
        )
        if type(item) is not PressureLossContribution:
            raise LossSeriesError(
                "an element returned something other than a PressureLossContribution"
            )
        contributions.append(item)
        total += item.drop_pa
        d_velocity += item.d_drop_d_superficial_velocity
        d_density += item.d_drop_d_gas_density
    return PressureLossDrop(
        drop_pa=total,
        d_drop_d_superficial_velocity=d_velocity,
        d_drop_d_gas_density=d_density,
        contributions=tuple(contributions),
        composition="ordered_series:" + "+".join(item.element_label for item in contributions),
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PressureDropCorroboration:
    """Where a converged layer drop sits against the two PRINTED bands."""

    drop_pa: float
    inside_frozen_ergun_band: bool
    inside_corroborated_band: bool
    inside_printed_hull: bool
    #: The INTERSECTION of the two printed bands: doubly-printed support, which is
    #: the strict reading of the design basis's corroboration sentence. REPORTED
    #: only -- the enforced gate is the hull, because narrowing a printed window is
    #: an owner ruling. A True here means both printed bands cover the drop; a
    #: False with ``inside_printed_hull`` True means single-source support only.
    inside_both_printed_bands: bool
    frozen_ergun_band_pa: tuple[float, float] = FROZEN_ERGUN_LAYER_DROP_BAND_PA
    corroborated_band_pa: tuple[float, float] = CORROBORATED_LAYER_DROP_BAND_PA
    printed_hull_pa: tuple[float, float] = CORROBORATED_LAYER_DROP_HULL_PA
    doubly_printed_intersection_pa: tuple[float, float] = (
        DOUBLY_PRINTED_LAYER_DROP_INTERSECTION_PA
    )

    #: The two printed bands OVERLAP, so the hull introduces no unsupported
    #: interior and no new band.
    hull_has_printed_support_throughout: ClassVar[bool] = True
    #: ... but its edges have SINGLE-source support, and the gate is the hull.
    enforced_window_is_the_hull_not_the_intersection: ClassVar[bool] = True


def classify_layer_pressure_drop(drop_pa: float) -> PressureDropCorroboration:
    """Report, never refuse: where the drop sits against the printed bands."""

    _require_nonnegative("drop_pa", drop_pa)
    ergun_low, ergun_high = FROZEN_ERGUN_LAYER_DROP_BAND_PA
    corr_low, corr_high = CORROBORATED_LAYER_DROP_BAND_PA
    hull_low, hull_high = CORROBORATED_LAYER_DROP_HULL_PA
    return PressureDropCorroboration(
        drop_pa=drop_pa,
        inside_frozen_ergun_band=ergun_low <= drop_pa <= ergun_high,
        inside_corroborated_band=corr_low <= drop_pa <= corr_high,
        inside_printed_hull=hull_low <= drop_pa <= hull_high,
        inside_both_printed_bands=(
            ergun_low <= drop_pa <= ergun_high and corr_low <= drop_pa <= corr_high
        ),
    )


def require_corroborated_layer_drop(drop_pa: float) -> PressureDropCorroboration:
    """Fail closed outside the hull of the two printed bands.

    The check is on the CONVERGED drop, never on a Newton iterate: an iterate
    outside the window is not a physics claim, and refusing one would be a bound
    on the solver rather than on the arithmetic.
    """

    corroboration = classify_layer_pressure_drop(drop_pa)
    if not corroboration.inside_printed_hull:
        low, high = CORROBORATED_LAYER_DROP_HULL_PA
        raise PressureBandError(
            f"the converged layer drop {drop_pa!r} Pa falls outside the printed hull "
            f"[{low}, {high}] Pa; this packet introduces no new pressure band"
        )
    return corroboration


# --- the assembled inputs --------------------------------------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class CellClosureInputs:
    """Everything one cell's closure reads. Per cell; nothing spans cells.

    ``wall_temperature_k`` is a FROZEN PARAMETER of the fast block, not an
    unknown: the metal is the slowest mode, so it is advanced by
    ``advance_wall_node`` at the macro-step scale from the converged gas and
    interface temperatures. ``macro_step_s`` is the caller-supplied solids-paced
    step the fast block solves inside.
    """

    grid: GridAuthority
    geometry: CellGeometry
    gas_boundary: GasBoundary
    solid: SolidSideAggregates
    hydraulics: LayerHydraulics
    film_area_law: FilmAreaLaw
    transfer: TransferCoefficients
    properties: DeclaredGasProperties
    equilibrium: InterfaceEquilibrium
    wall: WallNodeParameters
    wall_temperature_k: float
    macro_step_s: float
    gas_closure: GasClosure = GasClosure.QUASI_STEADY
    film_closure: FilmClosure = FilmClosure.QUASI_STEADY_ALGEBRAIC
    accumulation: AccumulationMode = AccumulationMode.TRANSIENT
    previous: CellFieldState | None = None
    #: The loss-law seam, injectable so a future relation's arrival can be
    #: exercised against the frozen sparsity before it is an enum member.
    loss_dispatch: object = layer_pressure_drop
    require_corroborated_pressure_drop: bool = True

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    #: TRUE, on this packet's own types only. The inherited typed floor keeps its
    #: own flag false; this packet is the first permitted to write the arithmetic
    #: and flips the flag nowhere but here.
    pressure_loss_arithmetic_implemented: ClassVar[bool] = True
    geometry_values_are_declared_assumptions: ClassVar[bool] = True
    #: There is no fifth subsystem and no global unknown.
    subsystem_count: ClassVar[int] = 4
    option3_tolerant_assembly: ClassVar[bool] = True

    def __post_init__(self) -> None:
        for name, expected in (
            ("grid", GridAuthority),
            ("geometry", CellGeometry),
            ("gas_boundary", GasBoundary),
            ("solid", SolidSideAggregates),
            ("hydraulics", LayerHydraulics),
            ("film_area_law", FilmAreaLaw),
            ("transfer", TransferCoefficients),
            ("properties", DeclaredGasProperties),
            ("equilibrium", InterfaceEquilibrium),
            ("wall", WallNodeParameters),
        ):
            if type(getattr(self, name)) is not expected:
                raise ClosureConfigurationError(f"{name} must be an exact {expected.__name__}")
        _require_positive("wall_temperature_k", self.wall_temperature_k)
        _require_positive("macro_step_s", self.macro_step_s)
        for name, expected_enum in (
            ("gas_closure", GasClosure),
            ("film_closure", FilmClosure),
            ("accumulation", AccumulationMode),
        ):
            if type(getattr(self, name)) is not expected_enum:
                raise ClosureConfigurationError(
                    f"{name} must be a {expected_enum.__name__} member"
                )
        # THE DECK AND THE SERIES MUST AGREE ON FAMILY. The series validates its own
        # label against the ruled element kinds, and the deck carries the inherited
        # family lattice, but nothing compared the two: a family-2 deck could run a
        # family-3 single-nozzle path, i.e. a bed-less pressure path on a bed deck,
        # and every band gate would still pass. Refused here, in the series' own
        # error type, because it is the series that is wrong for the deck.
        if self.hydraulics.series.family is not self.geometry.tray_type.family:
            raise LossSeriesError(
                f"the declared deck is family "
                f"{self.geometry.tray_type.family.d7_family_number} and the loss series "
                f"is labelled family {self.hydraulics.series.family.d7_family_number}; "
                "one cell has one floor family and the pressure path must be that "
                "family's ruled series"
            )
        # ONE void fraction, two consumers. A silent disagreement between the
        # hydraulic and the area law would be two different beds.
        bed_void = self.hydraulics.bed_void_fraction
        if bed_void is not None and bed_void != self.film_area_law.bed_void_fraction:
            raise ClosureConfigurationError(
                "the bed element and the area law must declare the same void fraction"
            )
        if bed_void is None and (
            PressureLossRelation.PACKED_BED_ERGUN
            in RULED_ELEMENT_KIND_SERIES[self.hydraulics.series.family]
        ):  # pragma: no cover - unreachable once the series/family checks hold
            raise LossSeriesError(
                "a family whose ruled series contains the packed-bed Ergun element "
                "must carry that element, so its void fraction is reconcilable"
            )
        # ONE OPEN-AREA FRACTION, TWO CONSUMERS, exactly as for the void fraction.
        # The floor element's throat fraction IS the deck's open-area fraction (the
        # ruled convention: the printed taper ratio relates the wide end to the
        # throat and the declared free area is the throat), and the inherited typed
        # floor enforces the printed family bracket on it. Without this check the
        # loss arithmetic could run on an open fraction the deck never declared and
        # the inherited provenance lattice would have refused.
        for element in self.hydraulics.series.elements:
            throat = getattr(element, "throat_area_fraction", None)
            if throat is None:
                continue
            if throat != self.geometry.tray_type.free_area_fraction:
                raise ClosureConfigurationError(
                    f"the floor element's throat area fraction {throat!r} disagrees with "
                    f"the declared deck's free area fraction "
                    f"{self.geometry.tray_type.free_area_fraction!r}; the throat IS the "
                    "declared open area under the ruled taper convention, so a "
                    "disagreement is two different decks -- refused, never reconciled "
                    "by averaging"
                )
        if not callable(self.loss_dispatch):
            raise ClosureConfigurationError("loss_dispatch must be callable")
        if type(self.require_corroborated_pressure_drop) is not bool:
            raise ClosureConfigurationError(
                "require_corroborated_pressure_drop must be an exact bool"
            )
        needs_previous = self.gas_storage_scale != 0.0 or self.film_storage_scale != 0.0
        if needs_previous and self.previous is None:
            raise ClosureConfigurationError(
                "a storage term needs the previous cell field state"
            )
        if self.previous is not None and type(self.previous) is not CellFieldState:
            raise ClosureConfigurationError("previous must be an exact CellFieldState")

    # --- the storage scales: the C5a switch, by parameter alone -------------

    @property
    def gas_storage_scale(self) -> float:
        if self.accumulation is AccumulationMode.ZERO_ACCUMULATION:
            return 0.0
        return 1.0 if self.gas_closure is GasClosure.DYNAMIC_EOS else 0.0

    @property
    def film_storage_scale(self) -> float:
        if self.accumulation is AccumulationMode.ZERO_ACCUMULATION:
            return 0.0
        return 1.0 if self.film_closure is FilmClosure.TRANSIENT_ACCUMULATION else 0.0

    @property
    def wall_storage_scale(self) -> float:
        if self.accumulation is AccumulationMode.ZERO_ACCUMULATION:
            return 0.0
        return 1.0 if self.wall.closure is WallClosure.DYNAMIC_NODE else 0.0

    @property
    def storage_scales_are_all_exactly_zero(self) -> bool:
        """The C5a zero-accumulation limit, reached through this same code path."""

        return (
            self.gas_storage_scale == 0.0
            and self.film_storage_scale == 0.0
            and self.wall_storage_scale == 0.0
        )

    @property
    def interface_branch(self) -> InterfaceBranch:
        return self.solid.interface_branch

    # --- the PHY-008 separation, DERIVED from the declared inputs -------------

    @property
    def wall_conductance_sum_w_k(self) -> float:
        """The four declared conductances the metal node loses heat through.

        The same four the wall row's denominator carries, WITHOUT the storage
        conductance: the time constant is a property of the declared metal and its
        declared couplings, not of the step it is integrated with.
        """

        return (
            self.wall.steam_side_ua_w_k
            + self.transfer.wall_to_interface_ua_w_k
            + self.transfer.wall_to_gas_ua_w_k
            + self.wall.ambient_ua_w_k
        )

    @property
    def wall_time_constant_s(self) -> float:
        """The metal time constant C / sum(UA), COMPUTED, never declared.

        FROZEN PHY-008's separation is the ratio of this to the solids-paced step,
        so this module derives the ratio rather than asserting the conclusion. A
        zero conductance sum is refused for the same reason ``advance_wall_node``
        refuses a zero denominator: the node's response is undetermined, and a
        default would be an invention.
        """

        conductance = self.wall_conductance_sum_w_k
        if conductance == 0.0:
            raise WallNodeError(
                "the wall node has no conductance at all, so it has no time constant; "
                "refused rather than defaulted"
            )
        return self.wall.thermal_capacity_j_k / conductance

    @property
    def wall_time_constant_to_macro_step_ratio(self) -> float:
        """tau_wall / t_res -- the quantity FROZEN PHY-008 brackets at 1.4-8.3x."""

        return self.wall_time_constant_s / self.macro_step_s

    @property
    def wall_is_slowest_mode(self) -> bool:
        """TRUE only where the DERIVED ratio lands inside the frozen bracket.

        This is the structural warrant for keeping the wall row outside the fast
        block, so it is tested rather than declared. A metal that is faster than the
        solids march reads False here, and the same value is surfaced on the solved
        wall row, where a consumer of one step can see it.
        """

        low, high = PHY008_WALL_RESIDENCE_RATIO_BRACKET
        return low <= self.wall_time_constant_to_macro_step_ratio <= high

    def inputs_digest(self) -> str:
        """Framed digest in this packet's own domain."""

        parts: list[bytes] = [
            CELL_CLOSURE_MAGIC.encode("ascii"),
            CELL_CLOSURE_VERSION.to_bytes(8, "big"),
            self.grid.cells.to_bytes(8, "big"),
            _float_part(self.grid.dt_s),
            _float_part(self.macro_step_s),
            self.geometry.tray.tray_id.encode("utf-8"),
            self.geometry.tray_type.pressure_loss_relation.value.encode("ascii"),
            # THE SERIES' OWN IDENTITY, not just the deck's selector. Without these
            # the ruled bed+floor pair and a bed-less single nozzle were
            # digest-identical, so a receipt could not tell them apart.
            self.hydraulics.series.family.d7_family_number.to_bytes(8, "big"),
            b"+".join(
                kind.value.encode("ascii")
                for kind in self.hydraulics.series.element_kinds
            ),
            self.gas_closure.value.encode("ascii"),
            self.film_closure.value.encode("ascii"),
            self.accumulation.value.encode("ascii"),
            self.hydraulics.network_mode.value.encode("ascii"),
            self.film_area_law.exponent.name.encode("ascii"),
            # THE LOSS-LAW SEAM'S IDENTITY. The dispatch is a caller-settable field
            # that sits outside the bed element's length-scale tripwire, so a
            # substituted composer must at least be visible in the packet's identity
            # rather than digest-identical to the module's own.
            (
                getattr(self.loss_dispatch, "__qualname__", None)
                or type(self.loss_dispatch).__name__
            ).encode("utf-8"),
            self.interface_branch.value.encode("ascii"),
            _float_part(self.geometry.gas_side_reference_area_m2),
            _float_part(self.wall_temperature_k),
        ]
        return _framed_digest(CELL_CLOSURE_DIGEST_DOMAIN, tuple(parts))


def residual_scales(inputs: CellClosureInputs) -> tuple[float, ...]:
    """Per-row reference magnitudes for the convergence test ONLY.

    The rows are a molar balance, an energy balance and a pressure balance, so a
    single absolute tolerance would mean three different things. These scales are
    computed from declared inputs, are used only to normalize the convergence test
    and the reported norm, and appear nowhere in the residual itself.
    """

    flow = inputs.gas_boundary.inlet_molar_flow_mol_s
    area = inputs.geometry.gas_side_reference_area_m2
    # Conductance times one kelvin, so a normalized energy residual reads as an
    # equivalent temperature imbalance in kelvin.
    conductance = (
        inputs.transfer.gas_side_heat_coefficient_w_m2_k * area
        + inputs.transfer.solid_to_interface_ua_w_k
        + inputs.transfer.wall_to_interface_ua_w_k
        + inputs.transfer.wall_to_gas_ua_w_k
        + flow * inputs.properties.molar_heat_capacity_j_mol_k
    )
    pressure = inputs.gas_boundary.downstream_boundary_pressure_pa
    volume = inputs.geometry.bed_volume_m3
    if inputs.film_closure is FilmClosure.TRANSIENT_ACCUMULATION:
        # A mass rate: the film's own storage over the macro step.
        film_hexane = (
            inputs.film_area_law.hexane_inventory_capacity_kg(volume) / inputs.macro_step_s
            + abs(inputs.solid.hexane_supply_to_film_kg_s)
        )
        film_water = (
            inputs.film_area_law.water_inventory_capacity_kg(volume) / inputs.macro_step_s
            + abs(inputs.solid.water_supply_to_film_kg_s)
        )
    else:
        # A saturation: already of order one, so its own scale is exactly one.
        film_hexane = 1.0
        film_water = 1.0
    if film_hexane == 0.0:
        film_hexane = 1.0
    if film_water == 0.0:
        film_water = 1.0
    return (flow, flow, conductance, pressure, conductance, film_hexane, film_water)


def declared_external_saturations(inputs: CellClosureInputs) -> tuple[float, float]:
    """FROZEN PHY-032's saturations, evaluated from the READ loadings.

    ``S_h = C_b X_f / rho_h,l`` and ``S_w = C_b X_w,f / rho_w,l`` with
    ``C_b = (1-eps) rho_dm,p / eps`` (``release/physics_decisions.yaml:1830-1833``).
    The bed dry-matter concentration is taken from the packet-derived holdup and
    the cell's own bed volume, so the frozen definition is applied to inventories
    the packet-carrier layer owns rather than to anything reconstructed here.
    """

    volume = inputs.geometry.bed_volume_m3
    holdup = inputs.solid.dry_matter_holdup_kg
    hexane = (
        inputs.solid.attached_hexane_loading_kg_kg
        * holdup
        / inputs.film_area_law.hexane_inventory_capacity_kg(volume)
    )
    water = (
        inputs.solid.external_water_loading_kg_kg
        * holdup
        / inputs.film_area_law.water_inventory_capacity_kg(volume)
    )
    return (hexane, water)


# --- the per-cell residual and its analytic Jacobian ------------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class CellEvaluation:
    """One residual/Jacobian evaluation, with every reportable intermediate."""

    branch: InterfaceBranch
    gas_side_active_area_m2: float
    external_liquid_saturation_total: float
    stefan_factor: float
    blowing_corrected_heat_factor: float
    hexane_molar_flux_mol_s: float
    water_molar_flux_mol_s: float
    #: The film's net accumulation rate: supply in minus evaporation out. Under
    #: the default algebraic film rows this is the rate the packet-carrier layer
    #: applies to its own inventory; under the transient rows it is driven to zero
    #: against the saturation change and is reported for the same audit.
    film_hexane_net_accumulation_kg_s: float
    film_water_net_accumulation_kg_s: float
    outlet_molar_flow_mol_s: float
    gas_density_kg_m3: float
    superficial_velocity_m_s: float
    layer_pressure_drop_pa: float
    #: The ordered series' PER-ELEMENT contributions, retained so the aggregate is
    #: never the only observable (the identifiability prohibition). Empty under the
    #: two nested pressure limits, which bypass the series by construction.
    loss_contributions: tuple[PressureLossContribution, ...]
    loss_composition: str
    interface_convection_w: float
    wall_to_gas_w: float
    wall_to_interface_w: float
    solid_to_interface_w: float
    residual: tuple[float, ...]
    jacobian: tuple[tuple[float, ...], ...]
    scaled_residual_norm: float

    rows: ClassVar[tuple[ResidualRow, ...]] = FAST_BLOCK_ROWS
    unknowns: ClassVar[tuple[str, ...]] = UNKNOWN_NAMES


def _pressure_row(
    inputs: CellClosureInputs, pressure: _Dual, velocity: _Dual, density: _Dual
) -> tuple[_Dual, float, tuple[PressureLossContribution, ...], str]:
    """The per-layer pressure row. SIMULTANEOUS with the gas, never lagged.

    The ordered series composes the drop and its two channel partials, which are
    lifted onto the unknowns HERE -- once, outside the elements. The row's shape
    and this row's Jacobian sparsity are fixed by the two channels and by nothing
    any element does.
    """

    boundary = inputs.gas_boundary.downstream_boundary_pressure_pa
    mode = inputs.hydraulics.network_mode
    if mode is PressureNetworkMode.SINGLE_NODE_UNIFORM:
        return (pressure - boundary, 0.0, (), "single_node_uniform_limit")
    if mode is PressureNetworkMode.FROZEN_RESISTANCE:
        frozen = inputs.hydraulics.frozen_resistance_drop_pa
        return (
            pressure - (boundary + frozen),
            frozen,
            (),
            "frozen_resistance_prescribed_profile_limit",
        )
    drop = inputs.loss_dispatch(
        series=inputs.hydraulics.series,
        superficial_velocity_m_s=velocity.value,
        gas_density_kg_m3=density.value,
    )
    if type(drop) is not PressureLossDrop:
        raise ClosureConfigurationError("the loss dispatch must return a PressureLossDrop")
    lifted = (
        _constant(drop.drop_pa)
        + (velocity - velocity.value) * drop.d_drop_d_superficial_velocity
        + (density - density.value) * drop.d_drop_d_gas_density
    )
    return (
        pressure - boundary - lifted,
        drop.drop_pa,
        drop.contributions,
        drop.composition,
    )


def _evaluate(inputs: CellClosureInputs, vector: tuple[float, ...]) -> CellEvaluation:
    """Assemble the seven fast-block rows and their exact partials."""

    if len(vector) != _N_UNKNOWNS:
        raise ClosureConfigurationError(f"the fast block takes {_N_UNKNOWNS} unknowns")
    properties = inputs.properties
    transfer = inputs.transfer
    branch = inputs.interface_branch
    zero = _constant(0.0)

    unknowns = tuple(_variable(vector[index], index) for index in range(_N_UNKNOWNS))
    gas_temperature = unknowns[0]
    hexane_fraction = unknowns[1]
    water_fraction = unknowns[2]
    layer_pressure = unknowns[3]
    interface_temperature = unknowns[4]
    hexane_saturation = unknowns[5]
    water_saturation = unknowns[6]

    if gas_temperature.value <= 0.0 or interface_temperature.value <= 0.0:
        raise StepRejected("a trial temperature left the physical domain")
    if layer_pressure.value <= 0.0:
        raise StepRejected("a trial layer pressure left the physical domain")
    # EACH COMPONENT FRACTION, THEN THE CARRIER. The carrier test alone admits a
    # negative component cancelled by a larger positive one: at x_h = -0.3 the
    # carrier is a healthy 1.04 while the mixture molar mass goes NEGATIVE, and the
    # negative density that follows then trips a validator reserved for construction
    # faults -- escaping the solver as the wrong tier. The seed tier already refuses
    # a negative fraction; this is the iterate tier's analogue.
    if hexane_fraction.value < 0.0 or water_fraction.value < 0.0:
        raise StepRejected("a trial bulk mole fraction left the physical domain")
    carrier_fraction = 1.0 - hexane_fraction.value - water_fraction.value
    if carrier_fraction <= 0.0:
        raise StepRejected("a trial bulk composition left no inert carrier")

    # --- ONE area scalar, multiplied ONCE ---------------------------------
    # PER COMPONENT AND THEN THE TOTAL, at the ITERATE tier, through the same helper
    # the seed tier calls: the total alone admits (-0.30, 0.50) as 0.20, which is
    # not a bed state. Refused, never projected back into the feasible set.
    require_component_saturations_admissible(
        branch=branch,
        hexane_saturation=hexane_saturation.value,
        water_saturation=water_saturation.value,
        at_seed=False,
    )
    saturation_total = zero
    if branch.hexane_is_active:
        saturation_total = saturation_total + hexane_saturation
    if branch.water_is_active:
        saturation_total = saturation_total + water_saturation
    reference_area = inputs.geometry.gas_side_reference_area_m2
    if branch.has_external_liquid:
        area = _saturated_area_dual(
            reference_area, saturation_total, inputs.film_area_law.exponent
        )
    else:
        # FROZEN PHY-023 at the hexane endpoint: a_pg = a_b, a_hg = 0.
        area = _constant(reference_area)
    if area.value <= 0.0:
        raise SaturationAdmissibilityError("the gas-side area vanished at this iterate")

    # --- FROZEN PHY-022: one coupled two-component vector flux -------------
    star_hexane = (
        inputs.equilibrium.hexane._dual(interface_temperature) / layer_pressure
        if branch.hexane_is_active
        else zero
    )
    star_water = (
        inputs.equilibrium.water._dual(interface_temperature) / layer_pressure
        if branch.water_is_active
        else zero
    )
    star_total = star_hexane + star_water
    bulk_total = (hexane_fraction if branch.hexane_is_active else zero) + (
        water_fraction if branch.water_is_active else zero
    )
    remaining_carrier = 1.0 - star_total
    if remaining_carrier.value <= 0.0:
        raise InterfaceStateError(
            "the interfacial partial-pressure sum reached the layer pressure: this is "
            "the crossover to the heat-controlled limit, which this packet does not "
            "select, and it is refused rather than smoothed"
        )
    driving = (star_total - bulk_total) / remaining_carrier
    stefan = _stefan_factor(driving)

    conductance = transfer.gas_side_molar_conductance_mol_m2_s
    hexane_flux = (
        area * ((star_hexane - hexane_fraction) * conductance) * stefan
        if branch.hexane_is_active
        else zero
    )
    water_flux = (
        area * ((star_water - water_fraction) * conductance) * stefan
        if branch.water_is_active
        else zero
    )
    total_flux = hexane_flux + water_flux

    heat_conductance = area * transfer.gas_side_heat_coefficient_w_m2_k
    blowing = _ackermann_factor(
        total_flux * properties.molar_heat_capacity_j_mol_k / heat_conductance
    )
    effective_heat = heat_conductance * blowing

    inlet = inputs.gas_boundary
    inlet_flow = inlet.inlet_molar_flow_mol_s
    outlet_flow = total_flux + inlet_flow
    if outlet_flow.value <= 0.0:
        raise StepRejected("the trial outlet flow is not positive")

    mixture_molar_mass = (
        hexane_fraction * properties.hexane_molar_mass_kg_mol
        + water_fraction * properties.water_molar_mass_kg_mol
        + (1.0 - hexane_fraction - water_fraction) * properties.inert_molar_mass_kg_mol
    )
    density = (
        layer_pressure
        * mixture_molar_mass
        / (gas_temperature * UNIVERSAL_GAS_CONSTANT_J_MOL_K)
    )
    velocity = (
        outlet_flow * mixture_molar_mass / (density * inputs.geometry.cross_section_m2)
    )
    # BELT AND BRACES ON THE TWO DERIVED CHANNELS. Any route to a non-physical
    # density or velocity is an inadmissible ITERATE, so it is refused here, in the
    # iterate tier, before the series ever sees it. The composer's own validators
    # then stay what they are -- construction-fault checks on declared inputs -- and
    # can no longer be the first thing an out-of-domain iterate meets.
    if density.value <= 0.0:
        raise StepRejected("a trial gas density is not positive")
    if velocity.value < 0.0:
        raise StepRejected("a trial superficial velocity is negative")

    pressure_row, drop_pa, contributions, composition = _pressure_row(
        inputs, layer_pressure, velocity, density
    )

    # --- energy, on the declared common datum ------------------------------
    datum = properties.energy_datum_temperature_k
    heat_capacity = properties.molar_heat_capacity_j_mol_k
    hexane_latent = properties.hexane_latent_heat_j_mol
    water_latent = properties.water_latent_heat_j_mol

    inlet_enthalpy = inlet_flow * (
        inlet.inlet_hexane_mole_fraction * hexane_latent
        + inlet.inlet_water_mole_fraction * water_latent
        + heat_capacity * (inlet.inlet_temperature_k - datum)
    )
    outlet_enthalpy = outlet_flow * (
        hexane_fraction * hexane_latent
        + water_fraction * water_latent
        + (gas_temperature - datum) * heat_capacity
    )
    evaporated_enthalpy = (
        hexane_flux * hexane_latent
        + water_flux * water_latent
        + total_flux * ((interface_temperature - datum) * heat_capacity)
    )
    interface_convection = effective_heat * (interface_temperature - gas_temperature)
    wall_to_gas = transfer.wall_to_gas_ua_w_k * (inputs.wall_temperature_k - gas_temperature)
    wall_to_interface = transfer.wall_to_interface_ua_w_k * (
        inputs.wall_temperature_k - interface_temperature
    )
    solid_to_interface = transfer.solid_to_interface_ua_w_k * (
        inputs.solid.solid_temperature_k - interface_temperature
    )

    # --- FROZEN PHY-004 storage: exactly zero unless the switch is engaged ---
    # Declared simplification when it IS engaged: the storage volume is the bed
    # void only. FROZEN PHY-024 declares explicit gas storage volumes (including
    # the omitted vessel headspace) only in the dynamic-EOS mode, and the
    # headspace is not declared here.
    #
    # AND THE VOID IS NOT PARTITIONED AGAINST THE LIQUID. This volume is the FULL
    # bed void, held constant in the external liquid saturation, while the film
    # capacities treat S_L = 1 as a void full of liquid -- so with both switches
    # engaged at a high saturation the two rows describe overlapping volumes. That
    # is a declared simplification with a stated consequence (no partial of the gas
    # capacity with respect to saturation, exact for the model as written and wrong
    # for the bed it represents), booked in LIMITATIONS. Multiplying by (1 - S_L)
    # here would be a physics change, so it is not taken silently.
    gas_scale = inputs.gas_storage_scale
    if gas_scale != 0.0:
        previous = inputs.previous
        void_volume = inputs.geometry.void_volume_m3(inputs.film_area_law.bed_void_fraction)
        moles = (
            layer_pressure
            * void_volume
            / (gas_temperature * UNIVERSAL_GAS_CONSTANT_J_MOL_K)
        )
        previous_moles = (
            previous.layer_pressure_pa
            * void_volume
            / (previous.gas_temperature_k * UNIVERSAL_GAS_CONSTANT_J_MOL_K)
        )
        rate = gas_scale / inputs.macro_step_s
        internal_capacity = properties.molar_heat_capacity_constant_volume_j_mol_k
        internal = (
            hexane_fraction * hexane_latent
            + water_fraction * water_latent
            + (gas_temperature - datum) * internal_capacity
        )
        previous_internal = (
            previous.gas_hexane_mole_fraction * hexane_latent
            + previous.gas_water_mole_fraction * water_latent
            + (previous.gas_temperature_k - datum) * internal_capacity
        )
        hexane_storage = (
            moles * hexane_fraction - previous_moles * previous.gas_hexane_mole_fraction
        ) * rate
        water_storage = (
            moles * water_fraction - previous_moles * previous.gas_water_mole_fraction
        ) * rate
        energy_storage = (moles * internal - previous_moles * previous_internal) * rate
    else:
        hexane_storage = zero
        water_storage = zero
        energy_storage = zero

    gas_hexane_row = (
        inlet_flow * inlet.inlet_hexane_mole_fraction
        + hexane_flux
        - outlet_flow * hexane_fraction
        - hexane_storage
    )
    gas_water_row = (
        inlet_flow * inlet.inlet_water_mole_fraction
        + water_flux
        - outlet_flow * water_fraction
        - water_storage
    )
    gas_energy_row = (
        inlet_enthalpy
        + evaporated_enthalpy
        + interface_convection
        + wall_to_gas
        - outlet_enthalpy
        - energy_storage
    )
    interface_energy_row = (
        effective_heat * (gas_temperature - interface_temperature)
        + solid_to_interface
        + wall_to_interface
        - (hexane_flux * hexane_latent + water_flux * water_latent)
    )

    # --- the film inventory rows, and the active set's IDENTITY ROWS --------
    # An inactive component is pinned by an exact identity row. That is an
    # active-set row: it is not a penalty, not a clamp and not a variable bound,
    # and it carries no tolerance.
    #
    # A STRUCTURAL FACT THE FROZEN NO-DUPLICATE-AREA RULE FORCES, recorded here
    # because it is not obvious and it decided the default. The co-located
    # two-liquid state owns ONE gas-side area, which is a function of the TOTAL
    # external liquid saturation alone. Every quantity that reaches the flux and
    # heat rows through that area therefore has IDENTICAL partials with respect to
    # the two saturations. A "supply equals evaporation" pair of algebraic film
    # rows is consequently RANK-DEFICIENT in the two-liquid branch: the shared
    # area constrains the total and says nothing about the split. So the default
    # algebraic film rows are the frozen PHY-032 DEFINITION rows -- the saturation
    # IS the read inventory's saturation -- which is also the architecturally
    # correct reading of the cell layer as a pure function returning sources and
    # fluxes while the packet-carrier layer owns the inventory. The switchable
    # transient mode closes the inventory here instead, and is well posed because
    # each row then carries its own capacity on its own diagonal.
    film_scale = inputs.film_storage_scale
    volume = inputs.geometry.bed_volume_m3
    declared_hexane, declared_water = declared_external_saturations(inputs)
    if branch.hexane_is_active and film_scale != 0.0:
        capacity = inputs.film_area_law.hexane_inventory_capacity_kg(volume)
        film_hexane_row = (hexane_saturation - inputs.previous.external_hexane_saturation) * (
            film_scale * capacity / inputs.macro_step_s
        ) - (
            inputs.solid.hexane_supply_to_film_kg_s
            - hexane_flux * properties.hexane_molar_mass_kg_mol
        )
    else:
        film_hexane_row = hexane_saturation - declared_hexane
    if branch.water_is_active and film_scale != 0.0:
        capacity = inputs.film_area_law.water_inventory_capacity_kg(volume)
        film_water_row = (water_saturation - inputs.previous.external_water_saturation) * (
            film_scale * capacity / inputs.macro_step_s
        ) - (
            inputs.solid.water_supply_to_film_kg_s
            - water_flux * properties.water_molar_mass_kg_mol
        )
    else:
        film_water_row = water_saturation - declared_water

    rows = (
        gas_hexane_row,
        gas_water_row,
        gas_energy_row,
        pressure_row,
        interface_energy_row,
        film_hexane_row,
        film_water_row,
    )
    residual = tuple(row.value for row in rows)
    jacobian = tuple(row.partials for row in rows)
    norm = 0.0
    for value, scale in zip(residual, residual_scales(inputs), strict=True):
        ratio = abs(value) / scale
        if ratio > norm:
            norm = ratio
    return CellEvaluation(
        branch=branch,
        gas_side_active_area_m2=area.value,
        external_liquid_saturation_total=saturation_total.value,
        stefan_factor=stefan.value,
        blowing_corrected_heat_factor=blowing.value,
        hexane_molar_flux_mol_s=hexane_flux.value,
        water_molar_flux_mol_s=water_flux.value,
        film_hexane_net_accumulation_kg_s=(
            inputs.solid.hexane_supply_to_film_kg_s
            - hexane_flux.value * properties.hexane_molar_mass_kg_mol
        ),
        film_water_net_accumulation_kg_s=(
            inputs.solid.water_supply_to_film_kg_s
            - water_flux.value * properties.water_molar_mass_kg_mol
        ),
        outlet_molar_flow_mol_s=outlet_flow.value,
        gas_density_kg_m3=density.value,
        superficial_velocity_m_s=velocity.value,
        layer_pressure_drop_pa=drop_pa,
        loss_contributions=contributions,
        loss_composition=composition,
        interface_convection_w=interface_convection.value,
        wall_to_gas_w=wall_to_gas.value,
        wall_to_interface_w=wall_to_interface.value,
        solid_to_interface_w=solid_to_interface.value,
        residual=residual,
        jacobian=jacobian,
        scaled_residual_norm=norm,
    )


def assemble_fast_block(
    inputs: CellClosureInputs, state: CellFieldState
) -> CellEvaluation:
    """The public per-cell assembly. PER CELL: nothing here spans cells."""

    if type(inputs) is not CellClosureInputs:
        raise ClosureConfigurationError("inputs must be an exact CellClosureInputs")
    if type(state) is not CellFieldState:
        raise ClosureConfigurationError("state must be an exact CellFieldState")
    return _evaluate(inputs, state.as_vector())


# --- the per-cell linear solve ---------------------------------------------


def _solve_linear_system(
    matrix: tuple[tuple[float, ...], ...], rhs: tuple[float, ...]
) -> tuple[float, ...]:
    """Dimension-generic Gaussian elimination with equilibration and pivoting.

    The size comes from the arguments, so the ruled successor's extra row and
    unknown need no change here.

    A NOTE ON PIVOTING, because the word invites a misreading: the row exchange
    below is a permutation of the LINEAR ALGEBRA, chosen for numerical
    conditioning. It carries no physical orientation, it changes no sign
    convention of any flux, and it is not a projection of anything onto anything.
    The forbidden orientation reversal and the forbidden projection onto the
    contact-admissible set are physical operations on the transported quantity,
    and neither exists anywhere in this module.
    """

    size = len(rhs)
    if len(matrix) != size:
        raise ClosureConfigurationError("the linear system must be square")
    augmented: list[list[float]] = []
    for index, row in enumerate(matrix):
        if len(row) != size:
            raise ClosureConfigurationError("the linear system must be square")
        augmented.append([*row, rhs[index]])
    for row in augmented:
        largest = 0.0
        for column in range(size):
            magnitude = abs(row[column])
            if magnitude > largest:
                largest = magnitude
        if largest == 0.0:
            raise SingularJacobianError("the analytic Jacobian has an empty row")
        for column in range(size + 1):
            row[column] = row[column] / largest
    for column in range(size):
        pivot_index = column
        best = abs(augmented[column][column])
        for candidate in range(column + 1, size):
            magnitude = abs(augmented[candidate][column])
            if magnitude > best:
                best = magnitude
                pivot_index = candidate
        if best == 0.0:
            raise SingularJacobianError("the analytic Jacobian is singular at this iterate")
        if pivot_index != column:
            augmented[column], augmented[pivot_index] = (
                augmented[pivot_index],
                augmented[column],
            )
        pivot = augmented[column][column]
        for target in range(column + 1, size):
            factor = augmented[target][column] / pivot
            if factor == 0.0:
                continue
            for index in range(column, size + 1):
                augmented[target][index] -= factor * augmented[column][index]
    solution = [0.0] * size
    for row_index in reversed(range(size)):
        total = augmented[row_index][size]
        for column in range(row_index + 1, size):
            total -= augmented[row_index][column] * solution[column]
        pivot = augmented[row_index][row_index]
        if pivot == 0.0:
            raise SingularJacobianError("the analytic Jacobian is singular at this iterate")
        total = total / pivot
        if not math.isfinite(total):
            raise SingularJacobianError("the Newton step is not representable")
        solution[row_index] = total
    return tuple(solution)


# --- the pre-integration admissibility check, whose failure is TERMINAL -----

#: THE EXPLICIT WAIVER. ``seed_gate`` is a required keyword on every entry point, so
#: an ungated solve has to be spelled out with this token: a caller that has no
#: per-discretization declaration to gate against says so, and the choice is
#: recorded on the solve as an artefact fact instead of being a forgotten default.
#: It is deliberately not a callable, so it cannot be mistaken for one.
SEED_GATE_DECLARED_ABSENT = "seed_gate_declared_absent"
#: The witness string a gated solve carries when the gate is anonymous.
UNNAMED_SEED_GATE_WITNESS = "seed_gate_callable"


def seed_gate_witness(seed_gate: object) -> str:
    """The recorded identity of whatever was passed as the gate.

    A gated and an ungated converged solve must not be equal objects, so the
    outcome of the gate decision is a FACT on the artefact rather than something a
    reader infers from the call site. Refuses anything that is neither a callable
    nor the declared-absent token.
    """

    if seed_gate is SEED_GATE_DECLARED_ABSENT:
        return SEED_GATE_DECLARED_ABSENT
    if not callable(seed_gate):
        raise ClosureConfigurationError(
            "seed_gate must be a callable pre-integration gate or the explicit "
            f"{SEED_GATE_DECLARED_ABSENT!r} token; an ungated solve is a declared "
            "choice and never a default"
        )
    name = getattr(seed_gate, "__qualname__", None) or getattr(
        seed_gate, "__name__", None
    )
    if type(name) is not str or not name.strip():
        return UNNAMED_SEED_GATE_WITNESS
    module = getattr(seed_gate, "__module__", None)
    return f"{module}.{name}" if type(module) is str and module else name


def require_seed_admissible(
    inputs: CellClosureInputs,
    seed: CellFieldState,
    *,
    seed_gate: object,
) -> str:
    """Refuse an inadmissible seed BEFORE any integration takes place.

    The gate is called exactly once, first, and is entirely opaque to this
    module: it receives the grid authority and nothing else, and any refusal it
    raises is re-raised as a TERMINAL ``SeedAdmissibilityError``. It is never
    caught and retried, never softened, and never followed by a repaired seed.

    ``seed_gate`` is REQUIRED: either a callable or ``SEED_GATE_DECLARED_ABSENT``.
    Returns the witness string naming which of the two was used.
    """

    if type(inputs) is not CellClosureInputs:
        raise ClosureConfigurationError("inputs must be an exact CellClosureInputs")
    if type(seed) is not CellFieldState:
        raise ClosureConfigurationError("seed must be an exact CellFieldState")
    witness = seed_gate_witness(seed_gate)
    if seed_gate is not SEED_GATE_DECLARED_ABSENT:
        try:
            seed_gate(inputs.grid)
        except Exception as error:
            raise SeedAdmissibilityError(
                f"the pre-integration gate refused this seed: {error}"
            ) from error

    branch = inputs.interface_branch
    if branch.has_external_liquid:
        # Each present component's own saturation, and then their total -- the same
        # helper the iterate tier calls, at the terminal tier.
        require_component_saturations_admissible(
            branch=branch,
            hexane_saturation=seed.external_hexane_saturation,
            water_saturation=seed.external_water_saturation,
            at_seed=True,
        )
        # The READ loadings must themselves imply an admissible bed state. A
        # loading pair whose frozen saturations sum to one or more is an
        # infeasible feed/bed state, refused before integration and never scaled
        # back into the feasible set.
        declared_hexane, declared_water = declared_external_saturations(inputs)
        declared_total = 0.0
        if branch.hexane_is_active:
            declared_total += declared_hexane
        if branch.water_is_active:
            declared_total += declared_water
        require_saturation_admissible(declared_total, at_seed=True)
    if not branch.hexane_is_active:
        if seed.external_hexane_saturation != 0.0:
            raise SeedAdmissibilityError(
                "the active set has no external hexane, so its saturation seed must "
                "be exactly zero"
            )
        if inputs.solid.hexane_supply_to_film_kg_s != 0.0:
            raise SeedAdmissibilityError(
                "a hexane supply into an absent external film has nowhere to go; "
                "phase appearance is a branch switch driven by the packet layer's "
                "own inventory, not something the cell solve initiates"
            )
    if not branch.water_is_active:
        if seed.external_water_saturation != 0.0:
            raise SeedAdmissibilityError(
                "the active set has no external water, so its saturation seed must "
                "be exactly zero"
            )
        if inputs.solid.water_supply_to_film_kg_s != 0.0:
            raise SeedAdmissibilityError(
                "a water supply into an absent external film has nowhere to go"
            )
    if seed.gas_temperature_k <= 0.0 or seed.interface_temperature_k <= 0.0:
        raise SeedAdmissibilityError("a seed temperature is outside the physical domain")
    if seed.layer_pressure_pa <= 0.0:
        raise SeedAdmissibilityError("a seed layer pressure is outside the physical domain")
    if seed.gas_hexane_mole_fraction < 0.0 or seed.gas_water_mole_fraction < 0.0:
        raise SeedAdmissibilityError("a seed gas mole fraction is negative")
    if seed.gas_hexane_mole_fraction + seed.gas_water_mole_fraction >= 1.0:
        raise SeedAdmissibilityError("the seed gas composition leaves no inert carrier")
    return witness


# --- the fast Newton block --------------------------------------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class CellSolve:
    """One converged fast block. Three subsystems; the wall is NOT here."""

    state: CellFieldState
    evaluation: CellEvaluation
    iterations: int
    tolerance: float
    #: WHICH PRE-INTEGRATION GATE RAN, as an artefact fact. Either the gate's own
    #: identity or ``SEED_GATE_DECLARED_ABSENT``, so a gated and an ungated
    #: converged solve are not equal objects and a receipt can tell them apart.
    seed_gate_witness: str

    rows: ClassVar[tuple[ResidualRow, ...]] = FAST_BLOCK_ROWS
    unknowns: ClassVar[tuple[str, ...]] = UNKNOWN_NAMES
    wall_row_is_absent_from_the_fast_block: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.seed_gate_witness) is not str or not self.seed_gate_witness.strip():
            raise ClosureConfigurationError(
                "a converged solve must record which pre-integration gate ran"
            )

    @property
    def was_gated(self) -> bool:
        """FALSE only where the caller declared the gate absent, by name."""

        return self.seed_gate_witness != SEED_GATE_DECLARED_ABSENT

    @property
    def scaled_residual_norm(self) -> float:
        return self.evaluation.scaled_residual_norm

    @property
    def branch(self) -> InterfaceBranch:
        return self.evaluation.branch


DEFAULT_NEWTON_TOLERANCE = 1.0e-10
DEFAULT_NEWTON_MAX_ITERATIONS = 60


def solve_cell_fast_block(
    inputs: CellClosureInputs,
    *,
    seed: CellFieldState,
    seed_gate: object,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
) -> CellSolve:
    """Solve the three fast subsystems implicitly, with the analytic Jacobian.

    The wall temperature is a FROZEN PARAMETER here; the wall's own row is
    advanced separately by ``advance_wall_node``, at the pace its time constant
    demands. The pressure, by contrast, is an unknown of this very block: it is
    simultaneous, never lagged.

    UNDAMPED, UNSAFEGUARDED, AND DELIBERATELY SO. There is no line search, no
    damping and no trust region, because the safeguards that would help here all
    reduce to clipping or projecting a trial iterate back into the admissible set
    -- forbidden as Newton machinery. Non-convergence is a typed refusal, and an
    inadmissible iterate is a typed step rejection.

    ``seed_gate`` is a REQUIRED keyword: a callable, or the explicit
    ``SEED_GATE_DECLARED_ABSENT`` token. Whichever it was is recorded on the
    returned solve.

    AND THE CONVERGED STATE IS CHECKED TOO. A converged answer this module's own
    seed gate would refuse is not an answer, so the feasible set is re-tested once
    the iteration stops. That refusal is an ITERATE-tier step rejection, because the
    state is a converged iterate and not a seed: a caller must be able to reject the
    step rather than be told its seed was inadmissible.
    """

    _require_positive("tolerance", tolerance)
    if type(max_iterations) is not int or max_iterations <= 0:
        raise ClosureConfigurationError("max_iterations must be a positive int")
    witness = require_seed_admissible(inputs, seed, seed_gate=seed_gate)

    vector = seed.as_vector()
    evaluation = _evaluate(inputs, vector)
    iterations = 0
    while evaluation.scaled_residual_norm > tolerance:
        if iterations >= max_iterations:
            raise NewtonConvergenceError(
                f"the fast block did not converge in {max_iterations} iterations "
                f"(scaled residual {evaluation.scaled_residual_norm!r}); the step is "
                "refused rather than safeguarded"
            )
        step = _solve_linear_system(
            evaluation.jacobian, tuple(-value for value in evaluation.residual)
        )
        vector = tuple(
            current + increment
            for current, increment in zip(vector, step, strict=True)
        )
        evaluation = _evaluate(inputs, vector)
        iterations += 1
    state = CellFieldState.from_vector(vector)
    require_component_saturations_admissible(
        branch=evaluation.branch,
        hexane_saturation=state.external_hexane_saturation,
        water_saturation=state.external_water_saturation,
        at_seed=False,
    )
    return CellSolve(
        state=state,
        evaluation=evaluation,
        iterations=iterations,
        tolerance=tolerance,
        seed_gate_witness=witness,
    )


# --- the one differential row, outside the fast block -----------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class WallNodeSolve:
    """FROZEN PHY-008's single metal row, with its four coupling terms exposed."""

    wall_temperature_k: float
    steam_side_w: float
    to_interface_w: float
    to_gas_w: float
    to_ambient_w: float
    storage_w: float
    residual_w: float
    row_derivative_w_k: float
    storage_scale: float
    #: THE DERIVED SEPARATION, carried on the solved row. ``C / sum(UA)`` from the
    #: declared capacity and the declared conductances, and the solids-paced step it
    #: is compared against -- both exposed so a consumer of one step can see the
    #: ratio the structure rests on rather than take it on trust.
    wall_time_constant_s: float
    solids_paced_macro_step_s: float

    row: ClassVar[ResidualRow] = ResidualRow.WALL_METAL_INTERNAL_ENERGY
    subsystem: ClassVar[CellSubsystem] = CellSubsystem.WALL_METAL
    coupling_term_count: ClassVar[int] = 4
    #: FROZEN PHY-008's audited bracket, carried so the comparison is on the record.
    phy008_ratio_bracket: ClassVar[tuple[float, float]] = PHY008_WALL_RESIDENCE_RATIO_BRACKET

    @property
    def time_constant_to_macro_step_ratio(self) -> float:
        """tau_wall / t_res, computed from the two carried quantities."""

        return self.wall_time_constant_s / self.solids_paced_macro_step_s

    @property
    def is_slowest_mode(self) -> bool:
        """FROZEN PHY-008's separation, DERIVED and therefore falsifiable.

        The metal time constant must exceed the per-tray solid residence by
        1.4-8.3x for this row to be the slowest mode, which is the structural
        warrant for advancing it at the macro-step scale rather than inside the fast
        block. This tests the DERIVED ratio against that bracket, so a
        parameterisation whose metal is faster than the solids march reads False
        here instead of being told it is the slowest mode by a literal.
        """

        low, high = PHY008_WALL_RESIDENCE_RATIO_BRACKET
        return low <= self.time_constant_to_macro_step_ratio <= high

    @property
    def ambient_loss_is_exactly_zero(self) -> bool:
        """FROZEN PHY-026's adiabatic limiting oracle, reached exactly."""

        return self.to_ambient_w == 0.0


def advance_wall_node(
    inputs: CellClosureInputs,
    *,
    gas_temperature_k: float,
    interface_temperature_k: float,
    previous_wall_temperature_k: float,
) -> WallNodeSolve:
    """Advance the wall metal internal energy by ONE differential row.

    Four coupling terms, exactly as the frozen decision routes them: steam-side
    condensation in, tray-surface conduction into the meal above, convection to
    the gas below, and ambient loss. The row is linear in the wall temperature, so
    the backward-Euler solve is closed-form and its analytic derivative is the
    denominator, reported as ``row_derivative_w_k``.

    The gas and interface temperatures are FROZEN PARAMETERS here, converged by
    the fast block first. That one-way ordering is the multirate separation the
    frozen time-constant audit requires; it is not the forbidden pressure lag,
    which does not exist because pressure is an unknown of the fast block.

    The solved row carries the DERIVED metal time constant and the step it is
    compared against, so ``is_slowest_mode`` is a test rather than an assertion. A
    node with a capacity but no conductance at all has no finite time constant, so
    the separation cannot be evaluated for it: that case refuses (``WallNodeError``)
    rather than reporting a fabricated ratio.
    """

    if type(inputs) is not CellClosureInputs:
        raise ClosureConfigurationError("inputs must be an exact CellClosureInputs")
    _require_positive("gas_temperature_k", gas_temperature_k)
    _require_positive("interface_temperature_k", interface_temperature_k)
    _require_positive("previous_wall_temperature_k", previous_wall_temperature_k)

    wall = inputs.wall
    transfer = inputs.transfer
    scale = inputs.wall_storage_scale
    storage_conductance = scale * wall.thermal_capacity_j_k / inputs.macro_step_s
    denominator = (
        storage_conductance
        + wall.steam_side_ua_w_k
        + transfer.wall_to_interface_ua_w_k
        + transfer.wall_to_gas_ua_w_k
        + wall.ambient_ua_w_k
    )
    if denominator == 0.0:
        raise WallNodeError(
            "the wall row has no storage and no conductance, so its temperature is "
            "undetermined; refused rather than defaulted"
        )
    numerator = (
        storage_conductance * previous_wall_temperature_k
        + wall.steam_side_ua_w_k * wall.steam_temperature_k
        + transfer.wall_to_interface_ua_w_k * interface_temperature_k
        + transfer.wall_to_gas_ua_w_k * gas_temperature_k
        + wall.ambient_ua_w_k * wall.ambient_temperature_k
    )
    temperature = numerator / denominator
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise WallNodeError("the wall row produced a non-physical temperature")
    steam_side = wall.steam_side_ua_w_k * (wall.steam_temperature_k - temperature)
    to_interface = transfer.wall_to_interface_ua_w_k * (temperature - interface_temperature_k)
    to_gas = transfer.wall_to_gas_ua_w_k * (temperature - gas_temperature_k)
    to_ambient = wall.ambient_ua_w_k * (temperature - wall.ambient_temperature_k)
    storage = storage_conductance * (temperature - previous_wall_temperature_k)
    residual = storage - steam_side + to_interface + to_gas + to_ambient
    return WallNodeSolve(
        wall_temperature_k=temperature,
        steam_side_w=steam_side,
        to_interface_w=to_interface,
        to_gas_w=to_gas,
        to_ambient_w=to_ambient,
        storage_w=storage,
        residual_w=residual,
        row_derivative_w_k=denominator,
        storage_scale=scale,
        wall_time_constant_s=inputs.wall_time_constant_s,
        solids_paced_macro_step_s=inputs.macro_step_s,
    )


# --- boundary inspectability ------------------------------------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class CellBoundaryLedger:
    """The conservation-checkable quantities the deferred proof will need.

    Grid authority, the seven extensive names, the energy datum and every
    converged row residual, exposed as one object. Stated honestly: the row
    residuals are the converged rows of the same closure, so they witness the
    solve rather than audit it from an independent formulation. What they do give
    is a single place where a later proof can read the boundary rates, the datum
    and the grid identity without reaching inside the solver.

    AND THEY CANNOT SEE THE ONE GAP THAT MATTERS, so it is differenced separately.
    Each row residual is self-consistent inside its own subsystem, which is exactly
    why none of them notices that the fast block credited the wall couplings at the
    OLD wall temperature while the wall row debited them at the NEW one.
    ``wall_coupling_split_imbalance_w`` is that difference, computed and exposed, so
    a step that does not conserve energy across the seam cannot read as one that
    does.
    """

    grid: GridAuthority
    energy_datum_id: str
    energy_datum_temperature_k: float
    inlet_molar_flow_mol_s: float
    outlet_molar_flow_mol_s: float
    hexane_molar_flux_mol_s: float
    water_molar_flux_mol_s: float
    layer_pressure_drop_pa: float
    hexane_molar_balance_residual_mol_s: float
    water_molar_balance_residual_mol_s: float
    gas_energy_balance_residual_w: float
    layer_pressure_balance_residual_pa: float
    interface_energy_balance_residual_w: float
    film_hexane_balance_residual_kg_s: float
    film_water_balance_residual_kg_s: float
    wall_row_residual_w: float | None
    #: THE WALL / FAST-BLOCK SEAM, DIFFERENCED. What the wall row debited to the gas
    #: and the interface minus what the fast block credited from them, i.e. exactly
    #: ``(UA_wg + UA_wi) * (T_wall_new - T_wall_old)``. ``None`` when no wall row was
    #: advanced, because then there is no seam to measure. Nonzero during a transient
    #: is EXPECTED and declared: it is the cost of the ruled one-way PHY-008 split,
    #: not a solver defect, and it is exactly 0.0 at the C5a fixed point.
    wall_coupling_split_imbalance_w: float | None
    #: Which pre-integration gate ran, carried through from the solve.
    seed_gate_witness: str
    pressure_drop_corroboration: PressureDropCorroboration

    extensive_names: ClassVar[tuple[str, ...]] = EXTENSIVE_NAMES
    rows: ClassVar[tuple[ResidualRow, ...]] = FAST_BLOCK_ROWS
    wall_rows: ClassVar[tuple[ResidualRow, ...]] = WALL_ROWS


def cell_boundary_ledger(
    inputs: CellClosureInputs,
    solve: CellSolve,
    *,
    wall: WallNodeSolve | None = None,
) -> CellBoundaryLedger:
    """Expose the converged cell's boundary quantities for inspection."""

    if type(solve) is not CellSolve:
        raise ClosureConfigurationError("solve must be an exact CellSolve")
    if wall is not None and type(wall) is not WallNodeSolve:
        raise ClosureConfigurationError("wall must be an exact WallNodeSolve")
    evaluation = solve.evaluation
    residual = evaluation.residual
    return CellBoundaryLedger(
        grid=inputs.grid,
        energy_datum_id=inputs.properties.energy_datum_id,
        energy_datum_temperature_k=inputs.properties.energy_datum_temperature_k,
        inlet_molar_flow_mol_s=inputs.gas_boundary.inlet_molar_flow_mol_s,
        outlet_molar_flow_mol_s=evaluation.outlet_molar_flow_mol_s,
        hexane_molar_flux_mol_s=evaluation.hexane_molar_flux_mol_s,
        water_molar_flux_mol_s=evaluation.water_molar_flux_mol_s,
        layer_pressure_drop_pa=evaluation.layer_pressure_drop_pa,
        hexane_molar_balance_residual_mol_s=residual[0],
        water_molar_balance_residual_mol_s=residual[1],
        gas_energy_balance_residual_w=residual[2],
        layer_pressure_balance_residual_pa=residual[3],
        interface_energy_balance_residual_w=residual[4],
        film_hexane_balance_residual_kg_s=residual[5],
        film_water_balance_residual_kg_s=residual[6],
        wall_row_residual_w=None if wall is None else wall.residual_w,
        wall_coupling_split_imbalance_w=(
            None
            if wall is None
            else (wall.to_gas_w + wall.to_interface_w)
            - (evaluation.wall_to_gas_w + evaluation.wall_to_interface_w)
        ),
        seed_gate_witness=solve.seed_gate_witness,
        pressure_drop_corroboration=classify_layer_pressure_drop(
            evaluation.layer_pressure_drop_pa
        ),
    )


# --- one macro step: fast block first, then the slow row --------------------


@dataclass(frozen=True, slots=True, kw_only=True)
class CellMacroStep:
    """One solids-paced macro step of one cell: the DAE, marched."""

    fast: CellSolve
    wall: WallNodeSolve
    ledger: CellBoundaryLedger

    #: The structural claim this type exists to make: the fast block solved with
    #: the wall temperature frozen, and the wall row advanced afterwards from the
    #: converged fast state. Two paces, one step.
    wall_row_advanced_outside_the_fast_block: ClassVar[bool] = True


def march_cell_macro_step(
    inputs: CellClosureInputs,
    *,
    seed: CellFieldState,
    previous_wall_temperature_k: float,
    seed_gate: object,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
) -> CellMacroStep:
    """The multirate nesting, made explicit and ordered.

    ONE SOURCE OF TRUTH FOR THE WALL TEMPERATURE. The wall temperature is both a
    validated field of ``inputs`` and this function's own argument, and the field is
    what ``inputs_digest`` hashes -- so silently preferring the argument would make a
    digest-keyed replay reproduce a different solve. A disagreement is therefore
    REFUSED rather than resolved: a caller marching in time re-declares the field
    each step, which is what keeps the recorded identity honest.
    """

    _require_positive("previous_wall_temperature_k", previous_wall_temperature_k)
    if inputs.wall_temperature_k != previous_wall_temperature_k:
        raise ClosureConfigurationError(
            f"the declared wall temperature {inputs.wall_temperature_k!r} K disagrees "
            f"with the previous wall temperature {previous_wall_temperature_k!r} K this "
            "step is to be marched from; the field is what the inputs digest records, "
            "so the conflict is refused rather than silently overridden"
        )
    frozen_wall = replace(inputs, wall_temperature_k=previous_wall_temperature_k)
    fast = solve_cell_fast_block(
        frozen_wall,
        seed=seed,
        tolerance=tolerance,
        max_iterations=max_iterations,
        seed_gate=seed_gate,
    )
    wall = advance_wall_node(
        frozen_wall,
        gas_temperature_k=fast.state.gas_temperature_k,
        interface_temperature_k=fast.state.interface_temperature_k,
        previous_wall_temperature_k=previous_wall_temperature_k,
    )
    if frozen_wall.require_corroborated_pressure_drop:
        require_corroborated_layer_drop(fast.evaluation.layer_pressure_drop_pa)
    return CellMacroStep(
        fast=fast,
        wall=wall,
        ledger=cell_boundary_ledger(frozen_wall, fast, wall=wall),
    )


# --- the C5a zero-accumulation fixed point, through the SAME code path -------


@dataclass(frozen=True, slots=True, kw_only=True)
class ZeroAccumulationFixedPoint:
    """The C5a limit, relaxed to, not reimplemented."""

    inputs: CellClosureInputs
    fast: CellSolve
    wall: WallNodeSolve
    sweeps: int
    wall_temperature_change_k: float
    tolerance: float

    #: C5a is the zero-accumulation limit of the identical representation, so this
    #: oracle uses ``solve_cell_fast_block`` and ``advance_wall_node`` unchanged.
    #: There is no separate steady-state implementation anywhere in this module.
    uses_the_transient_code_path: ClassVar[bool] = True

    @property
    def storage_terms_are_exactly_zero(self) -> bool:
        return self.inputs.storage_scales_are_all_exactly_zero and self.wall.storage_w == 0.0


def relax_to_zero_accumulation_fixed_point(
    inputs: CellClosureInputs,
    *,
    seed: CellFieldState,
    wall_seed_temperature_k: float,
    seed_gate: object,
    tolerance: float = 1.0e-9,
    max_sweeps: int = 400,
) -> ZeroAccumulationFixedPoint:
    """Relax to the C5a fixed point by setting every storage term to zero.

    Reached by PARAMETER: ``AccumulationMode.ZERO_ACCUMULATION`` drives the gas,
    film and wall storage scales to exactly 0.0, and the same residual, the same
    Jacobian and the same wall row then define the steady limit. The sweep is a
    fixed-point iteration on the wall temperature alone, because with no storage
    anywhere the wall row is algebraic in it.

    In this limit the wall row needs at least one nonzero conductance, or its
    temperature is undetermined and ``advance_wall_node`` refuses.

    TWO TIER CORRECTIONS THIS ORACLE OWES ITS CALLER. First, only sweep 0 is a SEED:
    every later sweep hands the previous sweep's CONVERGED state back in, so an
    admissibility failure there is a RELAXATION ITERATE leaving the domain and is
    re-raised as a step rejection rather than as the terminal seed refusal the
    caller is told never to retry. Second, this oracle honours
    ``require_corroborated_pressure_drop`` on the relaxed drop exactly as
    ``march_cell_macro_step`` does, so the same declared inputs do not get two
    different gate behaviours depending on which entry point was used.
    """

    _require_positive("wall_seed_temperature_k", wall_seed_temperature_k)
    _require_positive("tolerance", tolerance)
    if type(max_sweeps) is not int or max_sweeps <= 0:
        raise ClosureConfigurationError("max_sweeps must be a positive int")
    limit = replace(inputs, accumulation=AccumulationMode.ZERO_ACCUMULATION, previous=None)
    if not limit.storage_scales_are_all_exactly_zero:  # pragma: no cover - invariant
        raise ClosureConfigurationError("the zero-accumulation limit failed to engage")

    state = seed
    temperature = wall_seed_temperature_k
    change = float("inf")
    sweeps = 0
    fast: CellSolve | None = None
    wall: WallNodeSolve | None = None
    while sweeps < max_sweeps:
        stepped = replace(limit, wall_temperature_k=temperature)
        if sweeps == 0:
            fast = solve_cell_fast_block(stepped, seed=state, seed_gate=seed_gate)
        else:
            # A LATER SWEEP IS AN ITERATE, NOT A SEED. The gate callable is correctly
            # called once, at sweep 0; what was wrong is the TIER of the intrinsic
            # checks on the sweeps after it, which blamed a caller's seed for a state
            # this oracle itself produced.
            try:
                fast = solve_cell_fast_block(
                    stepped, seed=state, seed_gate=SEED_GATE_DECLARED_ABSENT
                )
            except SeedAdmissibilityError as error:
                raise StepRejected(
                    f"a relaxation iterate left the admissible domain on sweep "
                    f"{sweeps}: {error}"
                ) from error
        wall = advance_wall_node(
            stepped,
            gas_temperature_k=fast.state.gas_temperature_k,
            interface_temperature_k=fast.state.interface_temperature_k,
            previous_wall_temperature_k=temperature,
        )
        change = abs(wall.wall_temperature_k - temperature)
        state = fast.state
        temperature = wall.wall_temperature_k
        sweeps += 1
        if change <= tolerance:
            if limit.require_corroborated_pressure_drop:
                require_corroborated_layer_drop(fast.evaluation.layer_pressure_drop_pa)
            return ZeroAccumulationFixedPoint(
                inputs=replace(limit, wall_temperature_k=temperature),
                fast=fast,
                wall=wall,
                sweeps=sweeps,
                wall_temperature_change_k=change,
                tolerance=tolerance,
            )
    raise NewtonConvergenceError(
        f"the zero-accumulation relaxation did not reach its fixed point in "
        f"{max_sweeps} sweeps (last wall change {change!r} K)"
    )


# --- reporting --------------------------------------------------------------


def report_geometry_evidentiary_standing(geometry: CellGeometry) -> int:
    """The typed floor's ORDINAL standing rank, for gating or reporting only.

    Forbidden as a weight in any residual, average or blend. This is the single
    place the rank is read in this module, it is returned unmodified, and no
    arithmetic is performed on it anywhere.
    """

    if type(geometry) is not CellGeometry:
        raise ClosureConfigurationError("geometry must be an exact CellGeometry")
    return geometry.tray_type.geometry_provenance.evidentiary_standing_rank


@dataclass(frozen=True, slots=True)
class ClaimBoundary:
    """This packet's claim boundary, as machine-readable flags."""

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    #: Flipped on this packet's OWN types only; the inherited typed floor keeps
    #: its own flag false and is not modified.
    pressure_loss_arithmetic_implemented: ClassVar[bool] = True
    geometry_values_are_declared_assumptions: ClassVar[bool] = True
    production_wiring_added: ClassVar[bool] = False
    new_gas_source_node_added: ClassVar[bool] = False
    global_state_vector_present: ClassVar[bool] = False
    f_gate_advanced: ClassVar[bool] = False


CLAIM_BOUNDARY = ClaimBoundary()


__all__ = (
    "AccumulationMode",
    "BedErgunElement",
    "CELL_CLOSURE_DIGEST_DOMAIN",
    "CELL_CLOSURE_MAGIC",
    "CELL_CLOSURE_VERSION",
    "CLAIM_BOUNDARY",
    "CORROBORATED_LAYER_DROP_BAND_PA",
    "CORROBORATED_LAYER_DROP_HULL_PA",
    "CellBoundaryLedger",
    "CellClosureError",
    "CellClosureInputs",
    "CellEvaluation",
    "CellFieldState",
    "CellGeometry",
    "CellMacroStep",
    "CellSolve",
    "CellSubsystem",
    "ClaimBoundary",
    "ClausiusClapeyronSaturation",
    "ClosureConfigurationError",
    "CoefficientInferenceNotImplementedError",
    "DECLARED_ENGINEERING_ASSUMPTION",
    "DEFAULT_NEWTON_MAX_ITERATIONS",
    "DEFAULT_NEWTON_TOLERANCE",
    "DOUBLY_PRINTED_LAYER_DROP_INTERSECTION_PA",
    "DeclaredGasProperties",
    "DerivedWideEndNozzleInput",
    "ERGUN_INERTIAL_COEFFICIENT",
    "ERGUN_VISCOUS_COEFFICIENT",
    "EXTENSIVE_NAMES",
    "FAST_BLOCK_ROWS",
    "FLOOR_DISCHARGE_COEFFICIENT_SENSITIVITY_RANGE",
    "FROZEN_ERGUN_LAYER_DROP_BAND_PA",
    "FilmAreaExponent",
    "FilmAreaLaw",
    "FilmClosure",
    "GasBoundary",
    "GasClosure",
    "GridAuthority",
    "HYDRAULIC_LENGTH_SCALE_FIELD",
    "INTRAPARTICLE_DIFFUSION_RADIUS_M",
    "InterfaceBranch",
    "InterfaceEquilibrium",
    "InterfaceStateError",
    "JointIdentifiabilityRefusedError",
    "LayerDropObservation",
    "LayerHydraulics",
    "LayerLossSeries",
    "LengthScaleConfusionError",
    "LossRelationNotImplementedError",
    "LossSeriesError",
    "NewtonConvergenceError",
    "PARTICLE_MODEL_LENGTH_SCALE_FIELD",
    "PHY008_WALL_RESIDENCE_RATIO_BRACKET",
    "PRESSURE_SERIES_RULING",
    "PressureBandError",
    "PressureDropCorroboration",
    "PressureLossContribution",
    "PressureLossDrop",
    "PressureNetworkMode",
    "RULED_ELEMENT_KIND_SERIES",
    "ResidualRow",
    "SEED_GATE_DECLARED_ABSENT",
    "SaturationAdmissibilityError",
    "SeedAdmissibilityError",
    "SingularAreaDerivativeError",
    "SingularJacobianError",
    "SolidSideAggregates",
    "StepRejected",
    "TAPER_AREA_RATIO_BAND",
    "TaperedBoreFloorPassageElement",
    "TransferCoefficients",
    "UNANCHORED_BY_REFERENCE_MACHINE_FAMILIES",
    "UNIVERSAL_GAS_CONSTANT_J_MOL_K",
    "UNKNOWN_INDEX",
    "UNKNOWN_NAMES",
    "UNNAMED_SEED_GATE_WITNESS",
    "WALL_ROWS",
    "WallClosure",
    "WallNodeError",
    "WallNodeParameters",
    "WallNodeSolve",
    "ZeroAccumulationFixedPoint",
    "advance_wall_node",
    "assemble_fast_block",
    "build_family_loss_series",
    "cell_boundary_ledger",
    "classify_layer_pressure_drop",
    "declared_external_saturations",
    "film_area_saturation_derivative",
    "film_area_saturation_factor",
    "fit_series_coefficients_to_observation",
    "infer_series_coefficients_from_total_drop",
    "layer_pressure_drop",
    "march_cell_macro_step",
    "relax_to_zero_accumulation_fixed_point",
    "report_geometry_evidentiary_standing",
    "require_component_saturations_admissible",
    "require_corroborated_layer_drop",
    "require_saturation_admissible",
    "require_seed_admissible",
    "residual_scales",
    "seed_gate_witness",
    "solve_cell_fast_block",
    "sweep_floor_discharge_coefficient",
)
