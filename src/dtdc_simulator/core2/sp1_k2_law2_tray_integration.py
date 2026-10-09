"""Executable, nonqualifying SP1 K=2 / RTD=8 engineering macro-step.

This module joins two deliberately separate engineering slices without
retrofitting either one.  The zero-inert Law-2 cell solve supplies component
fluxes and the shared-old-wall update for each of two vertical SP1 cells.  A
deterministic adapter converts those accepted fluxes into weighted packet
inventory changes, while :mod:`k_cell_tray_host` owns RTD clocking, whole-packet
motion, cell cancellation ledgers, and its sealed evaluate/commit transaction.

The two gas cells are coupled serially.  The lower-cell outlet is the exact
upper-cell inlet, and a fixed-point iteration closes the one internal face
pressure.  Both cells must carry one exact gas-property/enthalpy convention,
including molar heat capacity and the numeric and named energy datum.  The two
bed resistances span their own K depths, while exactly one active floor passage
lives at SP1's lower physical boundary; the upper series carries only the typed
zero-resistance nested-limit placeholder required by the inherited series
schema.  Packet mass changes are allocated without clipping: evaporation is
proportional to the available external phase and condensation is proportional
to dry matter.  Any attempted depletion beyond the available external phase is
rejected.  For the context-pinned tag-1 geometry, static layer UA is partitioned
by representative-particle count.  The cell therefore sees the exact
packet-weighted UA temperature moment without collapsing the packet states.
Energy is not assigned from an invented particle heat-capacity law: the
accepted wall storage minus the accepted gas, jacket, and ambient boundary
transfer is the unique packet-energy transfer required by the enclosing cell
control volume.  Its conductive part is allocated from each packet's own
temperature and UA share; the remaining nonconductive energy is allocated by
representative-particle share.  A decoded-candidate ledger checks the realized
mass and energy changes.  A whole-tray ledger proves that the internal K-face
enthalpy cancels and that only bottom inlet, top outlet, jacket, and ambient
energy remain.

This is an ``ENGINEERING_FEASIBILITY_ONLY`` integration.  It is not F4,
physical qualification, plant calibration, production wiring, a production
particle codec, durable atomic persistence, or a complete inter-tray model.
The detached tag-1 high-loading adapter now reconstructs each layer's one exact
solid temperature closure from the accepted packet energy/state.  That adapter supports
only the pre-activation high-loading variant and evaluates free-water calorics
at one fixed tray pressure; it refuses the activation handoff and admits only
the replay-validated detached engineering K-event planner.  The mixed-temperature
closure remains bounded to identical pinned
representative geometry, static UA, and one shared K-cell film/interface
environment.  All qualification and release claims remain false.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import ClassVar, Protocol

from . import cell_closure as cc
from . import cell_engineering_feasibility as ef
from . import cell_native_caloric as native_caloric
from . import component_energy_datum_adapter as component_datum
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import qsc_layer_falling_rate_law as falling_rate
from . import sorption_interface
from . import through_bed_sorbed_water_law as sorbed_water
from .dtdc_stack import REFERENCE_TRAYS
from .tray_type import reference_tray_type

INTEGRATION_DIGEST_DOMAIN = "GT-PS-2/sp1-k2-law2-tray-integration/v1"
ENGINEERING_PACKET_MAGIC = b"GTPS2-EF-PACKET-V1\x00"
ENGINEERING_PACKET_CODEC_REGISTRY_DIGEST = (
    "sha256:" + hashlib.sha256(ENGINEERING_PACKET_MAGIC + b"/codec").hexdigest()
)
ENGINEERING_PACKET_AUDITOR_IDENTITY_DIGEST = (
    "sha256:" + hashlib.sha256(ENGINEERING_PACKET_MAGIC + b"/auditor").hexdigest()
)

SP1_LAYER_COUNT = 2
SP1_RTD_STAGE_COUNT = 8
DEFAULT_FACE_TOLERANCE = 1.0e-11
DEFAULT_FACE_MAX_ITERATIONS = 20
DEFAULT_RELATIVE_LIMIT = 1.0e-10
#: W3-1 (owner ruling recorded in docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md,
#: 2026-09-03 morning session; executed at HEAD f5a70bb): the fixed row scale
#: the Newton actually converges the cell hydraulic-pressure row against
#: (``cell_engineering_feasibility`` ``row_scales`` index 3 = 1.0e4 Pa, at
#: :1823 and :2276, and the ``/ 10_000.0`` scaling of the same row in every
#: native-source kernel).  One converged pressure row therefore guarantees no
#: better than ``newton_tolerance * NEWTON_PRESSURE_ROW_SCALE_PA`` Pa of
#: absolute closure, whatever the tray drop happens to be.
NEWTON_PRESSURE_ROW_SCALE_PA = 1.0e4
AGGREGATE_PRESSURE_BINARY64_RESOLUTION_LAW_ID = (
    "telescoped_newton_pressure_row_scale_plus_internal_face_slip_plus_fsum_recombination_v1"
)
GAS_BOUNDARY_ENTHALPY_CONVENTION = (
    "declared_molar_cp_temperature_offset_plus_composition_weighted_species_offsets_"
    "bound_to_native_packet_datum_v1"
)
UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID = (
    "static_uniform_representative_particle_ua_v1"
)
PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID = (
    "packet_local_conduction_plus_representative_particle_weight_share_remainder_v1"
)
#: B8 (owner-ruled GO 2026-08-31): declared engineering-lane low-steam drop
#: band variant - certified printed hull untouched, downward extension only.
DECLARED_LANE_DROP_BAND_LAW_ID = "DECLARED_LANE_LOW_STEAM_DROP_BAND_V1"
#: B1 stage 3c (owner-ruled GO 2026-09-01): the declared mass Biot number of
#: the receding-front series resistance, Bi_m = k_c R / D_eff at DT scale -
#: a first-order DECLARED_ENGINEERING_ASSUMPTION from the signed design
#: packet (~1e5 with k_c ~ 3e-3 mol/m2/s film conductance, R ~ 0.9 mm,
#: D_eff ~ 1e-9 m2/s).  It multiplies only the CORE_SUPPLIED arm's film
#: factor; the two-film law never reads it.
DECLARED_FALLING_RATE_MASS_BIOT = 1.0e5
#: B1 stage 3c crossing (measured wall1, 2026-09-01): a boundary state
#: (total = X_c*m + eps) is film-active BY CONSTRUCTION, so a
#: constant-state re-encode can never cross the variant seam and per-packet
#: proportional draws Zeno against the nil film.  The declared crossing:
#: when a layer's weighted attached-film aggregate is at or below this
#: threshold (engineering-nil: <= 1e-4 kg on tonne-scale trays, ~1e-8
#: relative) AND the bound codec speaks the subcritical variant, the CELL
#: VIEW reclassifies the remnant into the core, binds the falling-rate arm,
#: and the packets then cross one by one THROUGH the hindered draw (the
#: dual codec's overlap-band routing lands each side lawfully).  Packet
#: payloads and conservation are untouched - only the layer-aggregate view
#: reclassifies, and the reclassified mass is recorded first-class.
NIL_FILM_RECLASSIFICATION_THRESHOLD_KG = 1.0e-4
#: B1 STAGE 4 (owner ruling 2026-09-19, "I rule towards B1 stage 4 - proceed
#: with audit and gated build"; audit
#: ``docs/evidence/ruling_audit_2026-09-20/B1S4/AUDIT_B1S4.md``): the SAME
#: engineering-nil band, read for the EXTERNAL WATER film.  NO new number is
#: declared - the constant above is reused, and the declared row that owns
#: the criterion lives in ``through_bed_sorbed_water_law`` as
#: ``external_water_nil_band_reclassification``.  The criterion is GEOMETRIC,
#: not mass: from the shipped 3/R area law at R = 0.885 mm and 600 kg/m3 the
#: band is at most 0.035 monolayers of liquid water on any layer of this
#: tower (0.01039 on MN2 K1, 0.03462 on the wettest layer PD1 K1), two orders
#: below a continuous film.
#:
#: WHY IT IS NOT A DEADBAND ON THE FROZEN PREDICATE (the D7-c-pattern
#: sentence this build owes ``area.py``): the frozen PHY-007/023 event
#: convention requires the EXACT inventory predicate, taken from the exact
#: zero side with no deadband.  The band does NOT loosen that predicate.  It
#: reclassifies the layer INVENTORY VIEW - the remnant is booked as sorbed
#: (retained) water before any predicate is evaluated - so every consumer
#: downstream, the exhaustion gate and ``_core_supplied_water_closure``
#: included, still receives an EXACT zero and still compares it exactly.
#: Loosening the predicate itself would move frozen physics and is owner-only.
#:
#: The real mass is not moved here and is never discarded: the same step's
#: D9-a drain imbibes the sub-band film into the hexane-free shell on every
#: packet (capacity on the measured layer is 169.10 kg against a film of
#: 6.5e-10 kg, so ``drain_split_kg`` lands the film at exactly ``+0.0`` and is
#: exactly conservative), which is the mass-conserving handover.  The F-W1
#: deposit operator is NOT used: its landed residual leaves the packet
#: inventory.
#:
#: GATED ON A DECLARED ARM.  Exactly as the hexane band is gated on the bound
#: subcritical codec, the water band fires only where the layer's model
#: declares a C8 Tier 1 sorbed-water arm.  ``declared_sorbed_water_arm`` is
#: ``None`` on every certified run and every shipped model, so on those runs
#: this branch is unreachable and the kernel is bit-identical by construction.
NIL_WATER_FILM_RECLASSIFICATION_LAW_ID = "b1-stage4-external-water-nil-band-v1"
#: TAG-3 / C8 Tier 2b (owner-ruled 2026-09-06, "I rule all four as
#: recommended, proceed."): provenance class of the regime-C binding's one
#: declared value, and the law id of the binding as a whole.
DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"
#: D9-e FIX (2026-09-14): THIS kernel is what carries a layer's accepted
#: cell field state into the next interval as that layer's Newton seed
#: (``_solve_serial_k_cells`` below, and the same line in the single-wall
#: and generic integrations), and owner ruling D9-e made the ten-unknown
#: dry-shell state one that may be accepted and carried.  The rule that
#: settles what a carried seed means when the layer has LEFT the branch is
#: declared by the cell, which owns both widths, and is READ here so the
#: assumption is named where the carry happens.  Never re-declared.
D9E_CARRIED_SEED_WIDTH_ROW = ef.D9E_CARRIED_SEED_DECLARED_ASSUMPTIONS[
    "d9e_carried_seed_width_from_record"
]
REGIME_C_HEXANE_BINDING_LAW_ID = "tag3-regime-c-sorbed-hexane-activity-binding-v1"
#: The bed pressure the regime-C activity cap is evaluated at.  NOT a new
#: number: it is the value every authority already used inside this closure
#: evaluates at by default (``sorption_interface.critical_loading_kg_kg``,
#: ``sorption_interface.equilibrium_loading_kg_kg`` and
#: ``qsc_layer_falling_rate_law.equilibrium_floor_loading`` all default to it)
#: and the fixed caloric pressure the packet codecs are pinned to
#: (``particle.wet_core.WetCoreParams.pressure_pa``).  Written here only so
#: the cap's basis is explicit; every caller may override it.
REGIME_C_LAYER_PRESSURE_PA = 101_325.0
PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID = (
    "weighted_before_after_common_energy_nextafter_spacing_v1"
)
WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID = (
    "weighted_packet_cell_state_and_external_operation_nextafter_budget_v1"
)
COUPLED_WALL_ENERGY_BINARY64_RESOLUTION_LAW_ID = (
    "packet_cell_state_wall_product_and_row_operation_nextafter_budget_v1"
)


class SP1K2Law2IntegrationError(ValueError):
    """Base refusal for the bounded SP1 integration."""


class SP1K2Law2ConfigurationError(SP1K2Law2IntegrationError):
    """The declared host, geometry, or closure inputs are incompatible."""


class SP1K2Law2StepError(SP1K2Law2IntegrationError):
    """A coupled trial macro-step is inadmissible."""


class SP1K2LaneDropBandError(SP1K2Law2StepError):
    """A converged layer drop left the DECLARED engineering-lane band (B8)."""


class SP1K2SorptionFloorError(SP1K2Law2StepError):
    """A core-supplied layer reached the frozen a_h=1 sorption floor X_e.

    B1 stage 3c: below X_e the surface activity falls off the wet-core
    branch and the regime-C dry continuation (PART-02 machinery) owns the
    state - refused typed, never clamped or extrapolated.
    """


class SP1K2SorbedWaterNoShellError(SP1K2Law2StepError):
    """D9-b-2: a film-free layer whose hexane-free shell is closed.

    Raised by the water binder, never by the cell: under the D9-a routing the
    state does not occur, and this refusal is the guard that says so.
    """


class SP1K2SorbedWaterFloorError(SP1K2Law2StepError):
    """A sorbed-water layer left the accepted Luikov law's qualified domain.

    C8 Tier 1: below ``W_ref`` the state potential is singular and packet A5
    authorizes no continuation, so the arm refuses typed.  Continuing there
    is the Tier 2 DRY_SURFACE branch, which no code implements.
    """


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerPhaseExhaustionSignal:
    """Typed transition data riding a Law-2 phase-exhaustion refusal (B1).

    Carried unchanged through the through-bed and mechanism-interval
    rejection wrappers so a caller may derive the analytic exhaustion
    time and re-request the shortened interval under the endpoint
    layer's SHORTENED_REISSUE contract.  Pure data; no behavior.
    """

    vertical_layer_id: int
    phase_name: str
    layer_stock_kg: float
    demanded_kg: float
    layer_macro_step_s: float
    interval_end_time_s: float | None
    #: Option A (owner ruling 2026-09-28): the hexane class whose OWN stock the
    #: demand exceeded, or None on every aggregate (single-class) refusal -
    #: which is every refusal a layer without a composed per-class transfer can
    #: raise, so the payload of every pre-A refusal is unchanged field for field.
    hexane_class_key: str | None = None
    #: Level 0b (owner ruling 2026-09-28 on the film-class water-routing packet,
    #: option (a)): the executable WATER class whose own stock the demand
    #: exceeded on a layer whose transfer was composed per class, or None on
    #: every aggregate refusal - so every pre-0b payload is unchanged.
    water_class_key: str | None = None


class SP1K2Law2PhaseExhaustionError(SP1K2Law2StepError):
    """Law-2 layer phase-exhaustion refusal carrying the transition data.

    B1 stage 1 (owner-signed design, GT_PS2_B1_DESIGN_SIGNOFF_RECORD_
    2026-08-31): the layer-aggregate depletion condition — which is
    packet-independent — is additionally checked once per layer before
    the packet loop and raised TYPED, with the refusing layer, phase,
    stock, demand, and layer macro stride attached so a caller may
    derive the analytic exhaustion time and re-request the shortened
    interval under the endpoint layer's SHORTENED_REISSUE contract.
    The message is byte-identical to the historical refusal and the
    class subtypes SP1K2Law2StepError, so every existing handler,
    classifier, and test matches unchanged; the per-particle helper
    keeps its own check as defense in depth.  The step is rejected,
    never clamped — this type only ADDS evidence to the refusal.
    """

    def __init__(
        self,
        message: str,
        *,
        vertical_layer_id: int,
        phase_name: str,
        layer_stock_kg: float,
        demanded_kg: float,
        layer_macro_step_s: float,
        interval_end_time_s: float | None,
        hexane_class_key: str | None = None,
        water_class_key: str | None = None,
    ) -> None:
        super().__init__(message)
        self.vertical_layer_id = vertical_layer_id
        self.phase_name = phase_name
        self.layer_stock_kg = layer_stock_kg
        self.demanded_kg = demanded_kg
        self.layer_macro_step_s = layer_macro_step_s
        self.interval_end_time_s = interval_end_time_s
        self.hexane_class_key = hexane_class_key
        self.water_class_key = water_class_key
        self.signal = LayerPhaseExhaustionSignal(
            vertical_layer_id=vertical_layer_id,
            phase_name=phase_name,
            layer_stock_kg=layer_stock_kg,
            demanded_kg=demanded_kg,
            layer_macro_step_s=layer_macro_step_s,
            interval_end_time_s=interval_end_time_s,
            hexane_class_key=hexane_class_key,
            water_class_key=water_class_key,
        )


class StaleSP1K2Law2StateError(SP1K2Law2IntegrationError):
    """Commit was attempted against a stale or foreign coupled state."""


def _require_finite(name: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise SP1K2Law2ConfigurationError(f"{name} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise SP1K2Law2ConfigurationError(f"{name} must be strictly positive")


def _relative_residual(residual: float, *terms: float) -> float:
    return abs(residual) / max(1.0, *(abs(value) for value in terms))


def _binary64_resolution(value: float) -> float:
    """Conservative one-ULP radius around one finite accepted scalar."""

    lower = math.nextafter(value, -math.inf)
    upper = math.nextafter(value, math.inf)
    return max(value - lower, upper - value)


def _residual_beyond_resolution(residual: float, resolution_bound: float) -> float:
    if not math.isfinite(resolution_bound) or resolution_bound < 0.0:
        raise SP1K2Law2StepError("binary64 resolution bound is not finite and nonnegative")
    if abs(residual) <= resolution_bound:
        return 0.0
    return math.copysign(abs(residual) - resolution_bound, residual)


def _product_binary64_resolution(left: float, right: float, product: float) -> float:
    return math.fsum(
        (
            abs(right) * _binary64_resolution(left),
            abs(left) * _binary64_resolution(right),
            _binary64_resolution(product),
        )
    )


def aggregate_pressure_binary64_resolution_bound(
    *,
    newton_tolerance: float,
    pressure_row_count: int,
    face_tolerance: float = 0.0,
    face_pressures_pa: tuple[float, ...] = (),
    summed_terms_pa: tuple[float, ...] = (),
) -> float:
    """W3-1: the guaranteed absolute resolution (Pa) of a telescoped pressure ledger.

    Ruling W3-1 (docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md, 2026-09-03
    morning session; executed at HEAD ``f5a70bb``) replaces the aggregate
    pressure ledgers' divide-by-tray-drop acceptance normalization with the
    resolution the pressure rows themselves promise.  The construction is the
    module's OWN energy-row precedent (``_wall_row_binary64_resolution_bound``,
    ``_layer_wall_energy_binary64_resolution_bound``, compared through
    ``_residual_beyond_resolution``): a residual lying inside the row's derived
    resolution contributes exactly 0.0, and no ``relative_limit`` moves.

    Three additive mechanisms, summed with ``math.fsum`` exactly as
    ``_product_binary64_resolution`` sums one radius per operand:

    1. ``pressure_row_count * newton_tolerance * NEWTON_PRESSURE_ROW_SCALE_PA``
       - the telescoped aggregate ``endpoint_drop - ordered_element_sum`` is
       algebraically the sum of the per-cell Newton pressure-row residuals, and
       each of those rows is converged only against the fixed 1.0e4 Pa row
       scale.  K rows therefore guarantee K times 1e-6 Pa at the frozen 1e-10
       Newton gate.
    2. ``face_tolerance * max(1.0, |face pressure|)`` per internal face, for the
       arrangements that close the tray-internal faces with a Picard rather than
       as joint unknowns: the accepted face carries up to one face-tolerance of
       slip into the same telescoped sum (the bound gap named in
       docs/GT_PS2_WAVE3_SOLVER_STRUGGLE_AUDIT_2026-09-03.md section 6.1 note 2).
       Simultaneous kernels solve the faces exactly and pass no faces here.
    3. One ``_binary64_resolution`` radius per term entering the two orderings
       of the element/cell sums - the double-rounding difference that finding
       3.11 of the same audit measures at +-1 ULP (4.547e-13 Pa) and that no
       equality can be required to close.

    Mechanism 1 dominates by ~6 orders of magnitude everywhere measured.
    """

    _require_finite("aggregate pressure newton tolerance", newton_tolerance, positive=True)
    _require_finite("aggregate pressure face tolerance", face_tolerance)
    if type(pressure_row_count) is not int or pressure_row_count < 1:
        raise SP1K2Law2ConfigurationError(
            "telescoped pressure row count must be a positive exact int"
        )
    if face_tolerance < 0.0:
        raise SP1K2Law2ConfigurationError("aggregate pressure face tolerance must be nonnegative")
    for value in (*face_pressures_pa, *summed_terms_pa):
        _require_finite("aggregate pressure resolution term", value)
    return math.fsum(
        (
            pressure_row_count * newton_tolerance * NEWTON_PRESSURE_ROW_SCALE_PA,
            *(face_tolerance * max(1.0, abs(value)) for value in face_pressures_pa),
            *(_binary64_resolution(value) for value in summed_terms_pa),
        )
    )


def _bulk_gas_common_datum_enthalpy_j_mol(
    properties: cc.DeclaredGasProperties,
    *,
    temperature_k: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
    datum_adapter: component_datum.NumericComponentEnergyDatumAdapter,
) -> float:
    """Mirror the accepted cell residual's full bulk-gas datum convention."""

    datum_adapter.require_gas_properties(properties)
    hexane_shift, water_shift = datum_adapter.common_datum_species_shifts_j_mol
    return (
        properties.molar_heat_capacity_j_mol_k
        * (temperature_k - properties.energy_datum_temperature_k)
        + hexane_mole_fraction * (properties.hexane_latent_heat_j_mol + hexane_shift)
        + water_mole_fraction * (properties.water_latent_heat_j_mol + water_shift)
    )


#: CELL-02d W7: the native counterpart of GAS_BOUNDARY_ENTHALPY_CONVENTION.
#: The declared constant and its frozen digests are untouched; the native mode
#: declares its own convention identity beside them.
NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION = (
    "native_property_law_molar_mixture_enthalpy_pressure_explicit_virial_"
    "bound_to_native_packet_datum_v1"
)


def _native_datum_stream_enthalpy_j_mol(
    properties: native_caloric.NativeGasProperties,
    *,
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
) -> float:
    """CELL-02d W3: the native mirror of the bulk-gas stream enthalpy.

    Every stream is booked at the pressure of the cell state it represents
    (the W4 per-cell pressure convention); the water fraction must be the
    exact binary64 complement so the binary law receives one composition.
    """

    if water_mole_fraction != 1.0 - hexane_mole_fraction:
        raise SP1K2Law2ConfigurationError(
            "the native binary stream requires the exact complementary water fraction"
        )
    return native_caloric.mixture_enthalpy(
        temperature_k, pressure_pa, hexane_mole_fraction
    ).value_j_mol


def _stream_enthalpy_j_mol(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
    datum_adapter: object,
) -> float:
    """Dispatch one stream enthalpy across the closed caloric-mode union.

    The declared arm is the original common-datum arithmetic, verbatim, and
    ignores the pressure (its enthalpy is pressure-free by model form); the
    native arm evaluates the frozen mixture law at the W4 pressure.
    """

    if type(properties) is cc.DeclaredGasProperties:
        return _bulk_gas_common_datum_enthalpy_j_mol(
            properties,
            temperature_k=temperature_k,
            hexane_mole_fraction=hexane_mole_fraction,
            water_mole_fraction=water_mole_fraction,
            datum_adapter=datum_adapter,
        )
    return _native_datum_stream_enthalpy_j_mol(
        properties,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        hexane_mole_fraction=hexane_mole_fraction,
        water_mole_fraction=water_mole_fraction,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SP1KLayerEngineeringModel:
    """Static engineering declarations for one geometric SP1 K layer.

    Solid temperature is intentionally absent: every step decodes it from the
    prior accepted tag-1 packets in this K layer.
    """

    grid: cc.GridAuthority
    geometry: cc.CellGeometry
    hydraulics: cc.LayerHydraulics
    film_area_law: cc.FilmAreaLaw
    transfer: cc.TransferCoefficients
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties
    transport: ef.Law2TransportInputs
    wall: cc.WallNodeParameters
    hexane_supply_to_film_kg_s: float = 0.0
    water_supply_to_film_kg_s: float = 0.0
    packet_thermal_aggregation_law_id: str = UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID
    packet_energy_allocation_law_id: str = (
        PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID
    )
    #: B8 (owner-ruled GO 2026-08-31): DECLARED engineering-lane low-steam
    #: band variant.  ``None`` keeps the certified printed-hull gate inside
    #: the cell solve and every sealed model digest byte-identical.  A
    #: declared band may only EXTEND the printed hull downward (drying-bed
    #: turndown); the printed ceiling stays authoritative.  repr=False so
    #: the canonical model digest binds it explicitly and conditionally.
    declared_lane_drop_band_pa: tuple[float, float] | None = field(default=None, repr=False)
    #: C8 Tier 1 (owner-ruled D-C8-1..5, 2026-09-05): the DECLARED sorbed-water
    #: arm.  ``None`` - the default and every shipped model - keeps the sealed
    #: two-film law byte-identical AND keeps a film-free layer refusing with
    #: the frozen "Law 2 is executable only for TWO_EXTERNAL_LIQUIDS", so the
    #: capability cannot appear on any run that did not ask for it.  repr=False
    #: for the same reason as the lane band: the canonical model digest binds
    #: it explicitly and only when declared.
    declared_sorbed_water_arm: sorbed_water.SorbedWaterArm | None = field(default=None, repr=False)

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    solid_temperature_is_packet_derived: ClassVar[bool] = True
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False
    fixed_tray_pressure_calorics: ClassVar[bool] = True
    heterogeneous_packet_temperature_closure_implemented: ClassVar[bool] = True
    state_dependent_packet_ua_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for name, expected in (
            ("grid", cc.GridAuthority),
            ("geometry", cc.CellGeometry),
            ("hydraulics", cc.LayerHydraulics),
            ("film_area_law", cc.FilmAreaLaw),
            ("transfer", cc.TransferCoefficients),
            ("transport", ef.Law2TransportInputs),
            ("wall", cc.WallNodeParameters),
        ):
            if type(getattr(self, name)) is not expected:
                raise SP1K2Law2ConfigurationError(f"{name} must be an exact {expected.__name__}")
        # CELL-02d W1: the same closed caloric-mode union as the cell kernel.
        if type(self.properties) is not cc.DeclaredGasProperties and (
            type(self.properties) is not native_caloric.NativeGasProperties
        ):
            raise SP1K2Law2ConfigurationError(
                "properties must be an exact DeclaredGasProperties or NativeGasProperties"
            )
        # CELL-02d W2 forgery gate: the native datum identity may only ride on
        # the native closure type, never on declared constants.
        if type(self.properties) is cc.DeclaredGasProperties and (
            self.properties.energy_datum_id == native_caloric.NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID
        ):
            raise SP1K2Law2ConfigurationError(
                "the native zero-gauge datum identity cannot be carried by a "
                "declared constant-property closure"
            )
        _require_finite("hexane_supply_to_film_kg_s", self.hexane_supply_to_film_kg_s)
        _require_finite("water_supply_to_film_kg_s", self.water_supply_to_film_kg_s)
        if (
            self.packet_thermal_aggregation_law_id
            != UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID
        ):
            raise SP1K2Law2ConfigurationError(
                "the bounded slice requires static uniform representative-particle UA"
            )
        if (
            self.packet_energy_allocation_law_id
            != PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID
        ):
            raise SP1K2Law2ConfigurationError(
                "the bounded slice requires packet-local conduction plus a "
                "representative-particle weight-share remainder"
            )
        if self.declared_lane_drop_band_pa is not None:
            band = self.declared_lane_drop_band_pa
            if (
                type(band) is not tuple
                or len(band) != 2
                or any(type(value) is not float for value in band)
            ):
                raise SP1K2Law2ConfigurationError(
                    "the declared lane drop band must be an exact (float, float) tuple"
                )
            hull_low, hull_high = cc.CORROBORATED_LAYER_DROP_HULL_PA
            floor_pa, ceiling_pa = band
            if not 0.0 < floor_pa < hull_low or ceiling_pa != hull_high:
                raise SP1K2Law2ConfigurationError(
                    "the declared lane band may only extend the printed hull "
                    "downward: 0 < floor < 1000 Pa with the printed ceiling "
                    "unchanged (B8 scope)"
                )
        if self.declared_sorbed_water_arm is not None and (
            type(self.declared_sorbed_water_arm) is not sorbed_water.SorbedWaterArm
        ):
            raise SP1K2Law2ConfigurationError(
                "declared_sorbed_water_arm must be an exact SorbedWaterArm"
            )


def canonical_sp1_k2_layer_model_digest(
    models: tuple[SP1KLayerEngineeringModel, SP1KLayerEngineeringModel],
) -> str:
    """Bind every immutable layer-model declaration in canonical K order."""

    if (
        type(models) is not tuple
        or len(models) != 2
        or any(type(item) is not SP1KLayerEngineeringModel for item in models)
    ):
        raise SP1K2Law2ConfigurationError(
            "the layer-model digest requires two exact models in canonical K order"
        )
    digest = hashlib.sha256()
    for part in (
        (INTEGRATION_DIGEST_DOMAIN + "/layer-models").encode("ascii"),
        *(repr(model).encode("utf-8") for model in models),
    ):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    # B8: the lane band variant is repr-excluded, so it is bound explicitly
    # and only when declared - every sealed (None) digest stays byte-identical
    # and a band-variant model can never share a digest with a certified one.
    for index, model in enumerate(models):
        if model.declared_lane_drop_band_pa is not None:
            part = (f"lane-drop-band/{index}/{model.declared_lane_drop_band_pa!r}").encode("ascii")
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
    # C8 Tier 1: the sorbed-water arm is repr-excluded for the same reason and
    # is bound the same way - every sealed (None) digest stays byte-identical
    # and an armed model can never share a digest with a certified one.
    for index, model in enumerate(models):
        if model.declared_sorbed_water_arm is not None:
            part = (f"sorbed-water-arm/{index}/{model.declared_sorbed_water_arm!r}").encode("ascii")
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedSP1K2Law2TrayState:
    """One immutable in-memory state spanning host, fast seeds, and two walls."""

    tray: tray_host.AcceptedKCellTrayState
    layer_model_configuration_digest: str
    cell_field_states: tuple[ef.BinaryNoInertCellFieldState, ef.BinaryNoInertCellFieldState]
    wall_temperatures_k: tuple[float, float]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    complete_intertray_model: ClassVar[bool] = False
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.tray) is not tray_host.AcceptedKCellTrayState:
            raise SP1K2Law2ConfigurationError("tray must be an accepted K-cell host state")
        if self.tray.physical_tray_id != "SP1":
            raise SP1K2Law2ConfigurationError("the bounded integration accepts only SP1")
        if (
            type(self.layer_model_configuration_digest) is not str
            or len(self.layer_model_configuration_digest) != 71
            or not self.layer_model_configuration_digest.startswith("sha256:")
            or any(
                character not in "0123456789abcdef"
                for character in self.layer_model_configuration_digest[7:]
            )
        ):
            raise SP1K2Law2ConfigurationError(
                "accepted layer-model configuration requires a sha256 digest"
            )
        physical = next(
            item
            for item in self.tray.topology.trays
            if item.physical_tray_id == self.tray.physical_tray_id
        )
        if (
            physical.vertical_layer_count != SP1_LAYER_COUNT
            or physical.rtd_stage_count != SP1_RTD_STAGE_COUNT
        ):
            raise SP1K2Law2ConfigurationError("SP1 integration requires exactly K=2 and RTD=8")
        if type(self.cell_field_states) is not tuple or len(self.cell_field_states) != 2:
            raise SP1K2Law2ConfigurationError("accepted state requires two cell fast states")
        # D9-e: the accepted fast states may now be the TEN-unknown dry-shell
        # state as well as the shipped binary one, so the gate tests exact
        # membership in the cell's own exported tuple rather than one name.
        # Exactly as before it is an EXACT-type test: an ``isinstance`` here
        # would admit any future subclass silently, and the dry-shell state is
        # not a subclass of the binary one at all, so only naming it admits it.
        # A genuinely foreign type still refuses typed, on this same line.
        if any(
            type(item) not in ef.ACCEPTED_FAST_BLOCK_STATE_TYPES for item in self.cell_field_states
        ):
            raise SP1K2Law2ConfigurationError("accepted fast states have a foreign type")
        if type(self.wall_temperatures_k) is not tuple or len(self.wall_temperatures_k) != 2:
            raise SP1K2Law2ConfigurationError("accepted state requires two wall temperatures")
        for value in self.wall_temperatures_k:
            _require_finite("accepted wall temperature", value, positive=True)

    @property
    def state_digest(self) -> str:
        digest = hashlib.sha256()
        for part in (
            INTEGRATION_DIGEST_DOMAIN.encode("ascii"),
            self.tray.state_digest.encode("ascii"),
            self.layer_model_configuration_digest.encode("ascii"),
            repr(self.cell_field_states).encode("utf-8"),
            *(struct.pack(">d", value) for value in self.wall_temperatures_k),
        ):
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
        return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class SP1K2Law2StepRequest:
    """Declared dynamic boundaries for one coupled macro-step.

    The bottom boundary's downstream pressure is the initial internal-face
    pressure iterate.  Its inlet flow, composition, and temperature remain the
    lower-cell inlet.  ``top_boundary_pressure_pa`` is the upper-cell outlet
    pressure.
    """

    attempt_id: str
    prior: AcceptedSP1K2Law2TrayState
    end_time_s: float
    layer_models: tuple[SP1KLayerEngineeringModel, SP1KLayerEngineeringModel]
    bottom_gas_boundary: ef.EngineeringGasBoundary
    top_boundary_pressure_pa: float
    flow_segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...]
    newton_tolerance: float = ef.DEFAULT_NEWTON_TOLERANCE
    face_tolerance: float = DEFAULT_FACE_TOLERANCE
    face_max_iterations: int = DEFAULT_FACE_MAX_ITERATIONS
    relative_limit: float = DEFAULT_RELATIVE_LIMIT

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise SP1K2Law2ConfigurationError("attempt_id must be nonblank")
        if type(self.prior) is not AcceptedSP1K2Law2TrayState:
            raise SP1K2Law2ConfigurationError("prior has the wrong coupled-state type")
        _require_finite("end_time_s", self.end_time_s)
        if self.end_time_s <= self.prior.tray.time_s:
            raise SP1K2Law2ConfigurationError("the coupled macro-step must advance time")
        if type(self.layer_models) is not tuple or len(self.layer_models) != 2:
            raise SP1K2Law2ConfigurationError("the SP1 integration requires two layer models")
        if any(type(item) is not SP1KLayerEngineeringModel for item in self.layer_models):
            raise SP1K2Law2ConfigurationError("layer_models contains a foreign model type")
        if type(self.bottom_gas_boundary) is not ef.EngineeringGasBoundary:
            raise SP1K2Law2ConfigurationError("bottom_gas_boundary has the wrong type")
        if self.bottom_gas_boundary.carrier_topology is not ef.CarrierTopology.BINARY_NO_INERT:
            raise SP1K2Law2ConfigurationError("the bounded tray integration requires Law 2")
        _require_finite("top_boundary_pressure_pa", self.top_boundary_pressure_pa, positive=True)
        _require_finite("newton_tolerance", self.newton_tolerance, positive=True)
        _require_finite("face_tolerance", self.face_tolerance, positive=True)
        _require_finite("relative_limit", self.relative_limit, positive=True)
        if type(self.face_max_iterations) is not int or self.face_max_iterations <= 0:
            raise SP1K2Law2ConfigurationError(
                "face_max_iterations must be a strictly positive exact integer"
            )
        if type(self.flow_segments) is not tuple or not self.flow_segments:
            raise SP1K2Law2ConfigurationError("a nonempty accepted-flow history is required")
        if any(
            type(item) is not tray_host.AcceptedDryMatterFlowSegment for item in self.flow_segments
        ):
            raise SP1K2Law2ConfigurationError("flow history contains a foreign segment")
        if self.flow_segments[0].start_time_s != self.prior.tray.time_s:
            raise SP1K2Law2ConfigurationError("flow history does not start at accepted time")
        if self.flow_segments[-1].end_time_s != self.end_time_s:
            raise SP1K2Law2ConfigurationError("flow history does not end at requested time")
        _validate_sp1_layer_geometry(self.layer_models)
        _validate_one_gas_enthalpy_convention(self.layer_models)
        if canonical_sp1_k2_layer_model_digest(self.layer_models) != (
            self.prior.layer_model_configuration_digest
        ):
            raise SP1K2Law2ConfigurationError(
                "layer-model configuration drifted from the accepted-state digest"
            )
        for model in self.layer_models:
            if model.properties.energy_datum_id != self.prior.tray.energy_datum_id:
                raise SP1K2Law2ConfigurationError(
                    "cell model and packet host use different energy datums"
                )


def _validate_one_gas_enthalpy_convention(
    models: tuple[SP1KLayerEngineeringModel, SP1KLayerEngineeringModel],
) -> None:
    lower, upper = models
    if lower.properties != upper.properties:
        raise SP1K2Law2ConfigurationError(
            "both K cells must bind one exact gas-properties enthalpy convention, "
            "including its caloric mode and datum ID"
        )


def _validate_sp1_layer_geometry(
    models: tuple[SP1KLayerEngineeringModel, SP1KLayerEngineeringModel],
) -> None:
    reference = next(item for item in REFERENCE_TRAYS if item.tray_id == "SP1")
    tray_type = reference_tray_type(reference)
    expected_ids = ("SP1:K1", "SP1:K2")
    floors: list[cc.TaperedBoreFloorPassageElement] = []
    for expected_id, model in zip(expected_ids, models, strict=True):
        tray = model.geometry.tray
        if tray.tray_id != expected_id:
            raise SP1K2Law2ConfigurationError(
                f"SP1 layer geometry must be explicitly identified as {expected_id}"
            )
        if model.geometry.tray_type != tray_type:
            raise SP1K2Law2ConfigurationError("SP1 K layers require the inherited SP1 tray type")
        if tray.role != reference.role or tray.contact is not reference.contact:
            raise SP1K2Law2ConfigurationError("SP1 K layer changed its role or contact mode")
        if tray.diameter_m != reference.diameter_m:
            raise SP1K2Law2ConfigurationError("SP1 K layer changed the physical tray diameter")
        bed = model.hydraulics.series.bed_element
        if bed is None or bed.flow_length_m != tray.loaded_depth_m:
            raise SP1K2Law2ConfigurationError(
                "each K-layer Ergun length must equal that layer's loaded depth"
            )
        if bed.resistance_scale != 1.0:
            raise SP1K2Law2ConfigurationError(
                "each K layer must retain its active physical bed resistance"
            )
        floor_rows = tuple(
            item
            for item in model.hydraulics.series.elements
            if type(item) is cc.TaperedBoreFloorPassageElement
        )
        if len(floor_rows) != 1:
            raise SP1K2Law2ConfigurationError(
                "each typed K series must expose exactly one floor placeholder"
            )
        floors.append(floor_rows[0])
    if tuple(item.resistance_scale for item in floors) != (1.0, 0.0):
        raise SP1K2Law2ConfigurationError(
            "SP1 K=2 requires exactly one active physical floor at the lower boundary; "
            "the upper typed floor placeholder must be the exact zero-resistance limit"
        )
    if replace(floors[0], resistance_scale=0.0) != floors[1]:
        raise SP1K2Law2ConfigurationError(
            "the inactive upper floor placeholder must carry the same physical-floor "
            "declaration as the single active lower boundary"
        )
    if math.fsum(model.geometry.tray.loaded_depth_m for model in models) != (
        reference.loaded_depth_m
    ):
        raise SP1K2Law2ConfigurationError("the two K layers must partition SP1 bed depth exactly")
    if math.fsum(model.geometry.tray.indirect_duty_w for model in models) != (
        reference.indirect_duty_w
    ):
        raise SP1K2Law2ConfigurationError(
            "the two K layers must partition SP1 indirect duty exactly"
        )
    reference_geometry = cc.CellGeometry(tray=reference, tray_type=tray_type)
    if math.fsum(model.geometry.gas_side_reference_area_m2 for model in models) != (
        reference_geometry.gas_side_reference_area_m2
    ):
        raise SP1K2Law2ConfigurationError(
            "the two K layers must partition SP1 gas-side area exactly"
        )


@dataclass(frozen=True, slots=True)
class EngineeringFeasibilityPacketCodec:
    """Deterministic engineering-only payload codec and inventory auditor."""

    codec_registry_digest: str = ENGINEERING_PACKET_CODEC_REGISTRY_DIGEST
    auditor_identity_digest: str = ENGINEERING_PACKET_AUDITOR_IDENTITY_DIGEST
    source_identity: str = "engineering-feasibility-packet-codec-v1"

    production_codec_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def encode(self, inventory: tray_host.PerParticleInventory) -> bytes:
        if type(inventory) is not tray_host.PerParticleInventory:
            raise SP1K2Law2ConfigurationError("codec requires exact per-particle inventory")
        datum = inventory.energy_datum_id.encode("utf-8")
        if len(datum) > 65_535:
            raise SP1K2Law2ConfigurationError("energy datum is too long for engineering codec")
        return b"".join(
            (
                ENGINEERING_PACKET_MAGIC,
                struct.pack(">7dH", *inventory.as_tuple(), len(datum)),
                datum,
            )
        )

    def audit(self, payload_bytes: bytes) -> tray_host.PerParticleInventory:
        if type(payload_bytes) is not bytes or not payload_bytes.startswith(
            ENGINEERING_PACKET_MAGIC
        ):
            raise SP1K2Law2StepError("payload is outside the engineering codec registry")
        body = payload_bytes[len(ENGINEERING_PACKET_MAGIC) :]
        header_size = struct.calcsize(">7dH")
        if len(body) < header_size:
            raise SP1K2Law2StepError("engineering payload is truncated")
        *values, datum_length = struct.unpack(">7dH", body[:header_size])
        datum_bytes = body[header_size:]
        if len(datum_bytes) != datum_length:
            raise SP1K2Law2StepError("engineering payload datum length is inconsistent")
        try:
            datum = datum_bytes.decode("utf-8")
        except UnicodeDecodeError as error:
            raise SP1K2Law2StepError("engineering payload datum is not UTF-8") from error
        return tray_host.PerParticleInventory(
            dry_matter_kg=values[0],
            residual_oil_label_kg=values[1],
            attached_hexane_kg=values[2],
            internal_hexane_kg=values[3],
            external_water_kg=values[4],
            retained_water_kg=values[5],
            common_datum_energy_j=values[6],
            energy_datum_id=datum,
        )


class EngineeringPacketCodecAuditor(Protocol):
    codec_registry_digest: str
    auditor_identity_digest: str
    source_identity: str

    def encode(self, inventory: tray_host.PerParticleInventory) -> bytes: ...

    def audit(self, payload_bytes: bytes) -> tray_host.PerParticleInventory: ...


@dataclass(frozen=True, slots=True, kw_only=True)
class PacketLayerContribution:
    """One audited packet's weighted contribution to a static-UA K closure."""

    packet_id: str
    representative_particles: float
    temperature_k: float
    ua_share_w_k: float
    weighted_dry_matter_kg: float
    weighted_attached_hexane_kg: float
    weighted_external_water_kg: float
    weighted_common_datum_energy_j: float
    #: D19 (owner ruling 2026-09-23, option 1 of
    #: ``docs/GT_PS2_RULING_PACKET_D19_MIXED_LAYER_FRONT_2026-09-21.md``): the
    #: packet's own CORE hexane and OIL-LABEL aggregates.  The D9-a front is
    #: now read per packet, and the three authorities it reads need the
    #: packet's own internal loading and oil fraction; this row carried
    #: neither.  Both are the packet's weighted totals, in the same units and
    #: from the same ``weighted_totals.totals`` tuple as the four above
    #: (indices 3 and 1), so no new quantity is introduced and the layer
    #: aggregates remain their exact ``math.fsum``.  Required, not defaulted:
    #: a default would let a builder hand the front an oil-free, hexane-free
    #: packet that the regime authority would classify EXHAUSTED, which is a
    #: wrong answer rather than a refusal.
    weighted_internal_hexane_kg: float
    weighted_residual_oil_label_kg: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    state_dependent_ua: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class PacketDerivedLayerInventory:
    layer: int
    dry_matter_kg: float
    attached_hexane_kg: float
    #: B1 stage 3c: the sorbed/pore (core) hexane aggregate - the arm's
    #: supply once the attached film is exactly exhausted.
    internal_hexane_kg: float
    #: Layer oil-label aggregate; binds w_o for the frozen sorption floor.
    residual_oil_label_kg: float
    external_water_kg: float
    #: C8 Tier 1: the SORBED (retained) water aggregate - the sorbed arm's own
    #: stock, and the loading its activity is read at.  Additive: every
    #: pre-existing consumer of this inventory is unchanged, and the default
    #: keeps every pre-existing DIRECT construction (diagnostic fixtures that
    #: probe the closures) valid.  ``_layer_inventory`` always supplies it, so
    #: no production path can reach the default.
    retained_water_kg: float = 0.0
    common_datum_energy_j: float
    solid_temperature_k: float
    #: B1 stage 3c crossing: the engineering-nil film mass the CELL VIEW
    #: reclassified into the core this step (0.0 outside the crossing).
    nil_film_reclassified_kg: float
    #: B1 stage 4: the engineering-nil EXTERNAL WATER film the CELL VIEW
    #: reclassified into the retained (sorbed) water this step, 0.0 outside
    #: the crossing and on every run whose model declares no sorbed-water
    #: arm.  Additive with a default, for the reason the retained aggregate
    #: above carries one: every pre-existing DIRECT construction (the
    #: diagnostic fixtures that probe the closures) stays valid, and
    #: ``_layer_inventory`` always supplies it, so no production path can
    #: reach the default.
    nil_water_film_reclassified_kg: float = 0.0
    solid_to_interface_ua_w_k: float
    ua_share_sum_w_k: float
    ua_partition_residual_w_k: float
    ua_temperature_moment_residual_w: float
    contributions: tuple[PacketLayerContribution, ...]
    aggregation_law_id: str
    energy_allocation_law_id: str
    packet_count: int

    physically_qualifying: ClassVar[bool] = False
    solid_temperature_is_packet_derived: ClassVar[bool] = True
    heterogeneous_packet_temperatures_preserved: ClassVar[bool] = True


def _layer_inventory(
    tray: tray_host.AcceptedKCellTrayState,
    layer: int,
    *,
    model: SP1KLayerEngineeringModel,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> PacketDerivedLayerInventory:
    packets = tuple(
        sorted(
            (
                packet
                for packet in tray.packets
                if packet.owner_key.vertical_layer_id.value == layer
            ),
            key=lambda packet: packet.packet_id,
        )
    )
    if not packets:
        raise SP1K2Law2StepError(f"SP1:K{layer} has no weighted packets")
    decoded: list[high_loading.CompleteHighLoadingPacket] = []
    for packet in packets:
        codec.require_host_binding(
            configuration_digest=packet.configuration_digest,
            energy_datum_id=packet.entry.inventory.energy_datum_id,
            pressure_pa=packet.pressure_pa,
        )
        complete = codec.decode(packet.entry.payload_bytes)
        if complete.inventory != packet.entry.inventory:
            raise SP1K2Law2StepError(
                f"SP1:K{layer} packet envelope disagrees with its decoded tag-1 state"
            )
        if packet.entry.surface is None or packet.entry.surface != complete.surface:
            raise SP1K2Law2StepError(
                f"SP1:K{layer} packet lacks the exact decoded B2 surface envelope"
            )
        if struct.pack(">d", complete.state.dry_meal_mass_kg) != struct.pack(
            ">d", codec.sphere_params.dry_meal_mass_kg
        ):
            raise SP1K2Law2StepError(
                f"SP1:K{layer} packet changed the pinned representative geometry"
            )
        decoded.append(complete)
    representative_particles = tuple(
        packet.entry.weight.representative_particles for packet in packets
    )
    total_representative_particles = math.fsum(representative_particles)
    if total_representative_particles <= 0.0:  # host validates each weight; defensive aggregate
        raise SP1K2Law2StepError(f"SP1:K{layer} has no representative-particle weight")
    temperatures = tuple(item.state.temperature_k for item in decoded)
    temperature_anchor = min(temperatures)
    solid_temperature = (
        temperature_anchor
        + math.fsum(
            count * (temperature - temperature_anchor)
            for count, temperature in zip(representative_particles, temperatures, strict=True)
        )
        / total_representative_particles
    )
    declared_ua = model.transfer.solid_to_interface_ua_w_k
    ua_shares = tuple(
        declared_ua * count / total_representative_particles for count in representative_particles
    )
    rows = tuple(packet.weighted_totals.totals for packet in packets)
    contributions = tuple(
        PacketLayerContribution(
            packet_id=packet.packet_id,
            representative_particles=count,
            temperature_k=complete.state.temperature_k,
            ua_share_w_k=ua_share,
            weighted_dry_matter_kg=row[0],
            weighted_attached_hexane_kg=row[2],
            weighted_external_water_kg=row[4],
            weighted_common_datum_energy_j=row[6],
            # D19: the same weighted-totals tuple, indices 3 and 1.
            weighted_internal_hexane_kg=row[3],
            weighted_residual_oil_label_kg=row[1],
        )
        for packet, complete, count, ua_share, row in zip(
            packets,
            decoded,
            representative_particles,
            ua_shares,
            rows,
            strict=True,
        )
    )
    ua_share_sum = math.fsum(ua_shares)
    ua_temperature_moment = math.fsum(
        contribution.ua_share_w_k * contribution.temperature_k for contribution in contributions
    )
    attached_aggregate_kg = math.fsum(row[2] for row in rows)
    internal_aggregate_kg = math.fsum(row[3] for row in rows)
    nil_film_reclassified_kg = 0.0
    if (
        0.0 < attached_aggregate_kg <= NIL_FILM_RECLASSIFICATION_THRESHOLD_KG
        and getattr(codec, "tag2_codec", None) is not None
    ):
        # B1 stage 3c crossing (see NIL_FILM_RECLASSIFICATION_THRESHOLD_KG):
        # the cell view treats the nil film remnant as core mass so the
        # falling-rate arm binds; packet payloads and totals are untouched
        # and the packets cross the variant seam through the hindered draw.
        nil_film_reclassified_kg = attached_aggregate_kg
        internal_aggregate_kg = math.fsum((internal_aggregate_kg, attached_aggregate_kg))
        attached_aggregate_kg = 0.0
    external_water_aggregate_kg = math.fsum(row[4] for row in rows)
    retained_water_aggregate_kg = math.fsum(row[5] for row in rows)
    nil_water_film_reclassified_kg = 0.0
    if (
        0.0 < external_water_aggregate_kg <= NIL_FILM_RECLASSIFICATION_THRESHOLD_KG
        and getattr(model, "declared_sorbed_water_arm", None) is not None
    ):
        # B1 STAGE 4 (see NIL_WATER_FILM_RECLASSIFICATION_LAW_ID above): the
        # CELL VIEW books the engineering-nil water film as sorbed water, so
        # the layer-aggregate exhaustion gate and the C8 Tier 1 binder both
        # receive an EXACT zero and the sorbed arm binds for this step.  The
        # frozen exact-zero predicates are untouched; the packets' own
        # payloads and totals are untouched; the real mass is moved, in this
        # same step, by the D9-a drain, and the reclassified mass is recorded
        # first-class here.  Unreachable wherever no sorbed-water arm is
        # declared - which is every certified run.
        nil_water_film_reclassified_kg = external_water_aggregate_kg
        retained_water_aggregate_kg = math.fsum(
            (retained_water_aggregate_kg, external_water_aggregate_kg)
        )
        external_water_aggregate_kg = 0.0
    return PacketDerivedLayerInventory(
        layer=layer,
        dry_matter_kg=math.fsum(row[0] for row in rows),
        attached_hexane_kg=attached_aggregate_kg,
        internal_hexane_kg=internal_aggregate_kg,
        residual_oil_label_kg=math.fsum(row[1] for row in rows),
        external_water_kg=external_water_aggregate_kg,
        retained_water_kg=retained_water_aggregate_kg,
        nil_film_reclassified_kg=nil_film_reclassified_kg,
        nil_water_film_reclassified_kg=nil_water_film_reclassified_kg,
        common_datum_energy_j=math.fsum(row[6] for row in rows),
        solid_temperature_k=solid_temperature,
        solid_to_interface_ua_w_k=declared_ua,
        ua_share_sum_w_k=ua_share_sum,
        ua_partition_residual_w_k=ua_share_sum - declared_ua,
        ua_temperature_moment_residual_w=(ua_temperature_moment - declared_ua * solid_temperature),
        contributions=contributions,
        aggregation_law_id=model.packet_thermal_aggregation_law_id,
        energy_allocation_law_id=model.packet_energy_allocation_law_id,
        packet_count=len(packets),
    )


REGIME_C_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "regime_c_activity_cap_at_solid_temperature": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the pressure-feasible activity ceiling that keeps the interface "
                    "hexane fraction at or below one is evaluated at the layer's solid "
                    "temperature because the interface temperature is a solve unknown; "
                    "the cell's own composition domain guard remains the backstop"
                ),
                "value": "min(1, P_layer / Psat_h(T_solid))",
                "layer_pressure_pa": REGIME_C_LAYER_PRESSURE_PA,
                "sources": (
                    "sorption_interface.py F6 (the `equilibrium_loading_kg_kg` "
                    "docstring, and `hexane_activity_ceiling` which is the same "
                    "expression this binding calls rather than re-implements); "
                    "docs/GT_PS2_C8_TIER2A_DOUBLY_SORBED_INTERFACE_BUILD_RECORD_"
                    "2026-09-05.md section 8"
                ),
                "rejecting_outcome": (
                    "a tag-3 layer refused by the cell with 'trial intrinsic-binary "
                    "composition left [0, 1]' means the cap at the solid temperature "
                    "was insufficient: the cap must then move inside the cell to the "
                    "interface temperature (a design change), never a wider tolerance"
                ),
            }
        ),
    }
)


@dataclass(frozen=True, slots=True, kw_only=True)
class RegimeCHexaneBinding:
    """What the regime-C (tag-3) hexane arm bound, and whether the cap fired.

    The record the specification calls the step's declared-assumptions block.
    It is NOT a field of ``PacketDerivedLayerInventory``: that dataclass's repr
    is hashed into ``through_bed_k2_single_wall_tray_integration``'s prepared
    inventory digest, which is outside this packet's write set, so a field
    there would move a digest for a diagnostic.  This record is instead
    reproducible from the layer inventory alone through
    :func:`regime_c_hexane_binding`, which the closure itself calls.
    """

    layer: int
    mean_loading_kg_kg: float
    solid_temperature_k: float
    layer_pressure_pa: float
    lagged_hexane_flux_kg_m2_s: float
    #: The equilibrium/lagged-flux surface activity before any cap.
    surface_activity: float
    #: min(1, P / Psat_h(T_solid)) - the SAME expression
    #: ``sorption_interface.equilibrium_loading_kg_kg`` uses, called not copied.
    feasible_activity_cap: float
    #: What the cell is handed: min(surface_activity, feasible_activity_cap).
    hexane_activity: float
    hexane_activity_capped: bool
    #: Exactly 1.0 in regime C (the B1 double-count guard).
    hexane_conductance_factor: float
    declared_assumption_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION
    law_id: ClassVar[str] = REGIME_C_HEXANE_BINDING_LAW_ID

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def regime_c_hexane_binding(
    inventory: PacketDerivedLayerInventory,
    *,
    lagged_hexane_flux_kg_m2_s: float = 0.0,
    pressure_pa: float = REGIME_C_LAYER_PRESSURE_PA,
) -> RegimeCHexaneBinding:
    """The tag-3 (regime C) sorbed-hexane interface binding for one layer.

    The activity is the sorption authority's own surface activity at the
    layer's mean sorbed loading and the LAGGED hexane flux, capped at the
    pressure-feasible ceiling (the declared assumption
    ``regime_c_activity_cap_at_solid_temperature``).  The conductance factor
    is exactly 1.0: in regime C the hindrance IS the activity, and applying
    the receding-front factor as well would double-count it (the B1 guard;
    the Tier 2a diagnostic measured that every below-unit factor is refused
    by the cell anyway).

    ``lagged_hexane_flux_kg_m2_s`` defaults to 0.0, which ``surface_state``
    treats as the equilibrium surface.  See the build record: the accepted
    per-area hexane flux is not carried across macro-steps in
    ``AcceptedSP1K2Law2TrayState`` (its ``cell_field_states`` are the eight
    Newton unknowns, none of which is a species mass flux), so every binding
    in this build uses the equilibrium surface.  The parameter exists so a
    later packet can wire a carried flux without touching this law.
    """

    dry_mass = inventory.dry_matter_kg
    mean_loading = inventory.internal_hexane_kg / dry_mass
    oil_fraction = inventory.residual_oil_label_kg / dry_mass
    temperature_k = inventory.solid_temperature_k
    try:
        params = sorption_interface.ParticleSorptionParams(residual_oil_mass_fraction=oil_fraction)
        surface = sorption_interface.surface_state(
            mean_loading_kg_kg=mean_loading,
            temperature_k=temperature_k,
            hexane_flux_kg_m2_s=lagged_hexane_flux_kg_m2_s,
            params=params,
            pressure_pa=pressure_pa,
        )
        cap = sorption_interface.hexane_activity_ceiling(temperature_k, pressure_pa)
    except (sorption_interface.SorptionInterfaceError, ValueError) as error:
        raise SP1K2Law2StepError(
            f"K{inventory.layer} regime-C hexane binding authorities refused: {error}"
        ) from error
    activity = min(surface.surface_activity, cap)
    return RegimeCHexaneBinding(
        layer=inventory.layer,
        mean_loading_kg_kg=mean_loading,
        solid_temperature_k=temperature_k,
        layer_pressure_pa=pressure_pa,
        lagged_hexane_flux_kg_m2_s=lagged_hexane_flux_kg_m2_s,
        surface_activity=surface.surface_activity,
        feasible_activity_cap=cap,
        hexane_activity=activity,
        hexane_activity_capped=activity != surface.surface_activity,
        hexane_conductance_factor=1.0,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class _LayerHexaneAuthorities:
    """D9-2 (risk R-K1): the three hexane authorities of a layer, read ONCE.

    ``_core_supplied_hexane_closure`` and ``_layer_front_fraction`` both need
    the regime, the critical loading and the equilibrium floor loading at the
    same layer state.  The D9-1 proposal read them TWICE, in two places that
    would drift; they are read here instead, in one order, from one
    ``ParticleSorptionParams`` record, and both callers consume this.

    What is NOT shared: each caller keeps its own refusal wording and its own
    treatment of the regimes, because those differ by design - the closure
    refuses EXHAUSTED and binds DRY_SORBATE through the regime-C path, while
    the front reports zero for both.
    """

    internal_loading_kg_kg: float
    residual_oil_mass_fraction: float
    solid_temperature_k: float
    regime: sorption_interface.HexaneRegime
    critical_loading_kg_kg: float
    equilibrium_floor_loading_kg_kg: float

    @property
    def front_is_at_the_surface(self) -> bool:
        """At or above the critical loading: no dry shell, film control."""

        return self.internal_loading_kg_kg >= self.critical_loading_kg_kg

    @property
    def core_is_fully_receded(self) -> bool:
        """No liquid hexane anywhere: the whole particle is hexane-free shell."""

        return self.regime in (
            sorption_interface.HexaneRegime.EXHAUSTED,
            sorption_interface.HexaneRegime.DRY_SORBATE,
        )


def _layer_hexane_authorities(
    inventory: PacketDerivedLayerInventory,
    *,
    pressure_pa: float,
    context: str,
) -> _LayerHexaneAuthorities:
    """Read the regime, the critical loading and the floor, in that order.

    ``context`` names the caller in the refusal, so the two call sites keep
    the messages they had before the extraction.

    D19 (2026-09-23): the body moved WHOLE into :func:`_hexane_authorities_at`
    so that the same three reads, in the same order, can be taken at one
    PACKET's state as well as at the layer aggregate.  Nothing in the reads
    changed - this function forwards the layer's four scalars and is
    bit-identical to what it was.
    """

    return _hexane_authorities_at(
        dry_matter_kg=inventory.dry_matter_kg,
        internal_hexane_kg=inventory.internal_hexane_kg,
        residual_oil_label_kg=inventory.residual_oil_label_kg,
        temperature_k=inventory.solid_temperature_k,
        pressure_pa=pressure_pa,
        context=context,
    )


def _hexane_authorities_at(
    *,
    dry_matter_kg: float,
    internal_hexane_kg: float,
    residual_oil_label_kg: float,
    temperature_k: float,
    pressure_pa: float,
    context: str,
) -> _LayerHexaneAuthorities:
    """The three authorities at one dry-matter-bearing state, layer or packet.

    D19: the caller supplies the four scalars instead of an inventory, because
    the D9-a front now reads them once per PACKET.  The order of the three
    reads, their refusal class and the refusal wording are exactly what
    :func:`_layer_hexane_authorities` carried before the extraction.
    """

    dry_mass = dry_matter_kg
    internal_loading = internal_hexane_kg / dry_mass
    oil_fraction = residual_oil_label_kg / dry_mass
    try:
        params = sorption_interface.ParticleSorptionParams(residual_oil_mass_fraction=oil_fraction)
        regime = sorption_interface.classify(
            internal_loading, temperature_k, params, pressure_pa=pressure_pa
        )
        critical = sorption_interface.critical_loading_kg_kg(temperature_k, params)
        floor = falling_rate.equilibrium_floor_loading(temperature_k, oil_fraction)
    except (
        sorption_interface.SorptionInterfaceError,
        falling_rate.FallingRateLawError,
        ValueError,
    ) as error:
        raise SP1K2Law2StepError(f"{context} authorities refused: {error}") from error
    return _LayerHexaneAuthorities(
        internal_loading_kg_kg=internal_loading,
        residual_oil_mass_fraction=oil_fraction,
        solid_temperature_k=temperature_k,
        regime=regime,
        critical_loading_kg_kg=critical,
        equilibrium_floor_loading_kg_kg=floor,
    )


def _core_supplied_hexane_closure(
    inventory: PacketDerivedLayerInventory,
    *,
    lagged_hexane_flux_kg_m2_s: float = 0.0,
    pressure_pa: float = REGIME_C_LAYER_PRESSURE_PA,
) -> ef.CoreSuppliedHexaneClosure | None:
    """Bind the core-supplied hexane closure for a film-exhausted layer.

    ``None`` (any positive attached film) keeps the sealed two-film law
    byte-identical.  The REGIME decides the rest, and the regime comes from
    ``sorption_interface.classify`` - the authority, so the Q-TAG2-FLOOR label
    change above the n-hexane boiling point moves this switch with it:

    * FILM_ACTIVE (at or above pore saturation, the instant after the crossing
      before any front forms): a_h = 1 with factor exactly 1.0, the
      constant-rate continuation.  Unchanged.
    * RECEDING_FRONT (B1 stage 3c, the stage-3b binding): a_h = 1 exactly - a
      wet core persists, so the equilibrium surface activity is unity and ALL
      hindrance is transport - with the stage-3a receding-front series factor
      at the declared mass Biot number.  Unchanged.
    * DRY_SORBATE (TAG-3, owner ruling 2026-09-06): the sorbed-phase activity
      from the sorption authority, capped at the pressure-feasible ceiling,
      with factor exactly 1.0.  Until 2026-09-06 this raised
      ``SP1K2SorptionFloorError``.
    * EXHAUSTED (no positive hexane at all): refused typed, never clamped -
      there is no surface state to report.
    """

    if inventory.attached_hexane_kg != 0.0:
        return None
    # D9-2 (R-K1): the three authorities, read through the shared record.
    authorities = _layer_hexane_authorities(
        inventory, pressure_pa=pressure_pa, context="core-supplied closure"
    )
    internal_loading = authorities.internal_loading_kg_kg
    temperature_k = authorities.solid_temperature_k
    regime = authorities.regime
    critical = authorities.critical_loading_kg_kg
    floor = authorities.equilibrium_floor_loading_kg_kg
    if regime is sorption_interface.HexaneRegime.EXHAUSTED:
        raise SP1K2SorptionFloorError(
            f"K{inventory.layer} core hexane loading {internal_loading:.6e} kg/kg "
            f"is not positive at {temperature_k:.2f} K; there is no sorbed "
            "surface state to bind - reject, do not clamp"
        )
    if regime is sorption_interface.HexaneRegime.DRY_SORBATE:
        binding = regime_c_hexane_binding(
            inventory,
            lagged_hexane_flux_kg_m2_s=lagged_hexane_flux_kg_m2_s,
            pressure_pa=pressure_pa,
        )
        return ef.CoreSuppliedHexaneClosure(
            hexane_activity=binding.hexane_activity,
            hexane_conductance_factor=binding.hexane_conductance_factor,
            # D9-b (d): the DRY_SORBATE regime is the fully receded core - no
            # liquid hexane anywhere, so the whole particle is hexane-free
            # shell and the front fraction is exactly zero.
            front_fraction=0.0,
        )
    if authorities.front_is_at_the_surface:
        # The front sits AT the surface: zero shell resistance, film control.
        front = 1.0
        factor = 1.0
    else:
        try:
            front = falling_rate.receding_front_fraction(internal_loading, critical, floor)
            factor = falling_rate.falling_rate_flux_factor(front, DECLARED_FALLING_RATE_MASS_BIOT)
        except falling_rate.FallingRateLawError as error:
            raise SP1K2Law2StepError(f"falling-rate closure refused: {error}") from error
    # D9-b (d): the front travels to the cell on the closure, so the dry-shell
    # interface and the D9-a water routing read the kernel's OWN value and can
    # never drift from the series factor computed one line above.
    return ef.CoreSuppliedHexaneClosure(
        hexane_activity=1.0, hexane_conductance_factor=factor, front_fraction=front
    )


def _core_supplied_water_closure(
    inventory: PacketDerivedLayerInventory,
    arm: sorbed_water.SorbedWaterArm | None,
) -> ef.CoreSuppliedWaterClosure | None:
    """C8 Tier 1: bind the water arm for a layer with no external water film.

    ``None`` on three distinct grounds, and they are different findings:
    no arm is DECLARED for this model (the default: the frozen refusal for a
    film-free layer is kept exactly as it stands today), or the layer carries
    a positive external water film (the sealed two-film law, byte-identical).

    On the arm itself the activity is the accepted Luikov law with the
    high-loading bridge above its fitted ceiling
    (``through_bed_sorbed_water_law``), read at the layer's own
    dry-matter-based retained loading, and the conductance factor is the
    declared series factor.

    C8 Tier 2a (owner ruling 2026-09-05 item 3): a layer whose hexane film is
    ALSO exhausted is the frozen PHY-007 DRY_SURFACE branch, and that branch
    is now executable at the cell.  The water arm therefore binds here exactly
    as it does beside a hexane film - same law, same activity, same declared
    factor - and the hexane arm binds beside it, from its own binder, on its
    own B1 stage-3c terms.  Until 2026-09-05 this raised a typed refusal.
    """

    if arm is None:
        return None
    if inventory.external_water_kg != 0.0:
        return None
    # D9-a-3 (owner ruling 2026-09-13): the arm exchanges with the gas through
    # the hexane-free shell ONLY.  A film-free layer at a CLOSED shell has no
    # exchangeable water and no interface state to report, so the arm refuses
    # typed rather than binding on a shell of zero volume.  Under the D9-a
    # routing that state does not occur, because such a layer keeps its water
    # film; this refusal is the guard that says so, and its firing is kill
    # criterion K8.
    #
    # OWNER RULING D9-e (2026-09-13) MOVED THE BOUNDARY, and this comparison
    # with it.  The declared ceiling is now 1.0 exactly - the front AT the
    # surface, a shell of zero volume - so the test is "at or above", not
    # "above": at 0.99 <= f < 1 the shell has real volume, the arm binds, and
    # the cell's dry-shell interface carries it (closure (vi), which has a
    # root at every front from the crossing to half radius).  The withdrawn
    # 0.99 boundary refused those fronts, which is what killed the D9-3
    # honest lane at f = 0.991391.  Nothing else in this binder moved.
    front_fraction = _layer_front_fraction(inventory)
    if front_fraction >= sorbed_water.NO_SHELL_FRONT_FRACTION_CEILING:
        raise SP1K2SorbedWaterNoShellError(
            f"K{inventory.layer} has no shell to carry the hexane duty "
            f"(front fraction {front_fraction:.6f} at the declared "
            f"{sorbed_water.NO_SHELL_FRONT_FRACTION_CEILING}); a film-free layer "
            "at a closed shell is outside D9-a's routing - refuse, do not clamp"
        )
    dry_mass = inventory.dry_matter_kg
    loading = inventory.retained_water_kg / dry_mass
    try:
        activity = arm.activity(loading)
    except sorbed_water.SorbedWaterLawError as error:
        raise SP1K2SorbedWaterFloorError(
            f"K{inventory.layer} sorbed-water arm refused at retained loading "
            f"{loading:.6e} kg/kg dry: {error}"
        ) from error
    return ef.CoreSuppliedWaterClosure(
        water_activity=activity,
        water_conductance_factor=arm.water_conductance_factor,
    )


#: D19 (owner ruling 2026-09-23 on
#: ``docs/GT_PS2_RULING_PACKET_D19_MIXED_LAYER_FRONT_2026-09-21.md``,
#: option 1): the law id of the front a layer reports to the D9-a water
#: routing.  It declares NO constant - it is a weighting rule, and the row
#: below owns its criterion, its basis, its bracket and the four outcomes
#: that reject it.
D19_DRY_MATTER_WEIGHTED_FRONT_LAW_ID = "d19-dry-matter-weighted-packet-front-v1"

D9A_FRONT_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "dry_matter_weighted_packet_front": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D19_DRY_MATTER_WEIGHTED_FRONT_LAW_ID,
                "criterion": (
                    "the hexane front a MIXED layer reports to the D9-a water "
                    "routing: the DRY-MATTER-WEIGHTED mean of its packets' own "
                    "fronts, each read by the same three authorities, in the same "
                    "order, that the layer aggregate was read by.  The ground is "
                    "additivity of exchange area, not a fit: the sorbed-water arm "
                    "exchanges with the gas through each particle's hexane-free "
                    "shell, so the layer's exchangeable shell is the SUM over its "
                    "packets, and at the pinned monodisperse representative "
                    "geometry (one radius and one dry-meal mass per particle, "
                    "enforced packet by packet in _layer_inventory) a packet's "
                    "shell area is proportional to its dry matter.  Dry matter is "
                    "therefore the area weight, not a proxy for one.  What it "
                    "replaces is not a coarser version of the same rule: the "
                    "superseded read reported the WHOLE layer as film-bearing "
                    "whenever ANY packet carried attached film, so one 3.27 g "
                    "arrival sliver closed the shell of 6,229.9 kg of dry bulk"
                ),
                "value": (
                    "f_layer = sum_i m_i f_i / sum_i m_i over the layer's packet "
                    "rows, with m_i the packet's weighted dry matter and f_i the "
                    "existing three-authority read on that packet's own state "
                    "(attached film: 1; fully receded core, i.e. DRY_SORBATE or "
                    "EXHAUSTED: 0; otherwise the falling-rate law's own "
                    "receding_front_fraction).  Both sums are math.fsum"
                ),
                "no_new_number": (
                    "the rule declares no constant, no threshold and no "
                    "tolerance.  Every number it reads - the critical loading, "
                    "the equilibrium floor, the regime boundaries, the no-shell "
                    "ceiling of exactly 1.0 - is the one the superseded read "
                    "already read, from the same authorities"
                ),
                "bracket": [1.9707288582193954e-03, 1.0],
                "bracket_ends": (
                    "both ends are MEASUREMENTS at one state, MN1 layer 1 of lane "
                    "op200s4pc4 at its last committed step 23: "
                    "1.9707288582193954e-03 is what this rule returns there and "
                    "1.0 is what the superseded maximum-over-packets read "
                    "returned at the same state.  The control lane "
                    "op200s4pc2ctl's pair at its step 30 is "
                    "2.288013544749551e-03 and 1.0.  No intermediate rule was "
                    "selected; the span is printed so the size of the change is "
                    "on the record and not inferred"
                ),
                "superseded_prediction": (
                    "the ruling packet's section 3 predicted 5.244e-07 and "
                    "9.213e-07 at these two states, reading the adversarial "
                    "pass's FILM-BEARING DRY-MATTER SHARE as the answer and "
                    "calling it an upper bound.  This build measured the rule "
                    "instead, and the share is neither the front nor a bound on "
                    "it: at the pc4 state nineteen sub-kg PD3 transfer packets "
                    "(13.44 kg of 6229.92 kg dry matter) are wet cores reading "
                    "fronts of 0.8327 to 1.0, not the zero the share assumed, "
                    "and they carry the whole weighted mean.  The prediction is "
                    "SUPERSEDED BY MEASUREMENT.  The conclusion it supported is "
                    "unchanged and is what matters: the front sits three orders "
                    "below the declared no-shell ceiling of exactly 1.0, so the "
                    "arm binds on the dry bulk's open shell instead of refusing"
                ),
                "uniform_layer_identity": (
                    "a UNIFORM layer - every packet reading the same f_i - "
                    "returns that f_i bit for bit, because the weighted mean of "
                    "one repeated value IS that value and the implementation "
                    "returns it rather than recomputing m*f/m, which is not "
                    "exact in binary64.  All film-bearing therefore returns "
                    "exactly 1.0, all fully receded exactly 0.0, and a layer with "
                    "one packet row its own front.  A layer with NO packet rows "
                    "(the diagnostic inventories built directly, which carry "
                    "contributions=()) is its own single row and reads exactly as "
                    "it did before this rule.  The extension is therefore live "
                    "only where the packets actually disagree"
                ),
                "sources": (
                    "owner ruling 2026-09-23 (the D19 line of 'I rule D19 and D20 "
                    "as recommended'), on the packet's recommendation 'D19 option "
                    "1, build and refreeze'; packet "
                    "docs/GT_PS2_RULING_PACKET_D19_MIXED_LAYER_FRONT_2026-09-21.md "
                    "section 3 option 1 and section 4; the measurement "
                    "docs/evidence/d9a_no_shell_mixed_layer_2026-09-21/ "
                    "(MEASUREMENT.md, ADVERSARIAL.md, adversarial_states.csv, "
                    "2,784 committed layer-states over six lanes); the D9-a design "
                    "brief docs/GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md and the "
                    "owner rulings docs/GT_PS2_OWNER_RULINGS_D9_2026-09-13.md; the "
                    "build record docs/GT_PS2_D19_LEVEL0_BUILD_RECORD_2026-09-24.md"
                ),
                "scope_amended_2026_09_28": (
                    "level 0b (owner ruling 2026-09-28 on the film-class water-routing "
                    "packet, option (a), component (i); law "
                    "level0b-class-own-water-routing-v1): this weighted front stays the "
                    "LAYER's report - the aggregate water binder, the layer water gate, "
                    "and the D9-a-2 drain and D9-a-3 floor of every layer whose accepted "
                    "transfer is not composed per class, which is every layer of every "
                    "lane with the seam off.  On a layer whose accepted transfer carries "
                    "water_by_class, the drain and the floor read each packet's OWN "
                    "three-authority front instead (_packet_front_fraction on the class "
                    "sub-inventory's rows): a packet with an attached hexane film has no "
                    "open shell, its front is exactly 1, its capacity exactly 0.0, and its "
                    "water film stays.  This is D9-a-2 as its design brief wrote it "
                    "('f = 1 while the layer carries an attached hexane film'), applied at "
                    "the object that carries the film, and this row's own criterion taken "
                    "literally: the exchangeable shell is the SUM over the packets, and a "
                    "sum of (1 - f_i^3) is not (1 - f^3) at the mean"
                ),
                "rejecting_outcome": (
                    "FOUR outcomes reject this row, each written before the build. "
                    "(1) ACCEPTANCE: the acceptance lane - the production point "
                    "with Stage B, stage 4 and the per-class composition on - "
                    "refuses at or before 1.3408 plant s for the same D9-a reason; "
                    "the rule would then not be describing the state that stopped "
                    "it. (2) CONSERVATION: the water ledger of any committed mixed "
                    "state fails the kernel's pinned 1e-10 relative clause; a "
                    "re-routed drain would then be leaking mass. (3) IDENTITY: a "
                    "uniform-layer read differs from the superseded aggregate read "
                    "by one ulp on any existing fixture; the rule is an EXTENSION, "
                    "and an extension that moves a uniform layer is a replacement. "
                    "(4) ATTRIBUTION: a leaf of the before/after twin moves on a "
                    "lane state that is not one of the 410 committed mixed states "
                    "the adversarial pass predicted, nor downstream of one"
                ),
            }
        ),
    }
)


@dataclass(frozen=True, slots=True, kw_only=True)
class _PacketFrontRow:
    """D19: one row the dry-matter-weighted front is read over.

    ``is_layer_aggregate`` is True for the single row a layer with no packet
    rows stands in with, and it exists only so that such a layer keeps the
    refusal wording it had before D19 - the row is otherwise read exactly as
    a packet row is.
    """

    packet_id: str
    is_layer_aggregate: bool
    dry_matter_kg: float
    attached_hexane_kg: float
    internal_hexane_kg: float
    residual_oil_label_kg: float
    temperature_k: float


def _packet_front_rows(
    inventory: PacketDerivedLayerInventory,
) -> tuple[_PacketFrontRow, ...]:
    """The rows the D19 weighted front reads, in the CELL VIEW's own terms.

    Two things are read through the view rather than off the payloads, and
    both keep a pre-existing law intact:

    * the B1 stage 3c nil-film band.  When the layer view has reclassified
      its engineering-nil attached film into the core, every packet's own
      attached film is reclassified the same way, packet by packet, so the
      per-packet sums still reproduce the layer aggregates and a layer in the
      band still reads as the film-free layer the band declares it to be;
    * a layer with NO packet rows is its own single row.  The diagnostic
      inventories built directly carry ``contributions=()`` with a declared
      ``packet_count``, and the layer aggregate IS the only row that exists
      for them, so the weighted mean over it is the aggregate read.
    """

    contributions = inventory.contributions
    if not contributions:
        return (
            _PacketFrontRow(
                packet_id="",
                is_layer_aggregate=True,
                dry_matter_kg=inventory.dry_matter_kg,
                attached_hexane_kg=inventory.attached_hexane_kg,
                internal_hexane_kg=inventory.internal_hexane_kg,
                residual_oil_label_kg=inventory.residual_oil_label_kg,
                temperature_k=inventory.solid_temperature_k,
            ),
        )
    reclassified = inventory.nil_film_reclassified_kg
    rows: list[_PacketFrontRow] = []
    for item in contributions:
        attached = item.weighted_attached_hexane_kg
        internal = item.weighted_internal_hexane_kg
        if reclassified != 0.0:
            internal = math.fsum((internal, attached))
            attached = 0.0
        rows.append(
            _PacketFrontRow(
                packet_id=item.packet_id,
                is_layer_aggregate=False,
                dry_matter_kg=item.weighted_dry_matter_kg,
                attached_hexane_kg=attached,
                internal_hexane_kg=internal,
                residual_oil_label_kg=item.weighted_residual_oil_label_kg,
                temperature_k=item.temperature_k,
            )
        )
    return tuple(rows)


def _packet_front_fraction(
    row: _PacketFrontRow,
    *,
    pressure_pa: float,
    layer: int,
) -> float:
    """D19: ONE row's front, by the three authorities D9-a-4 already named.

    The cases the D9-a design brief names, unchanged in every respect except
    the state they are read at:

    * a row carrying a positive ATTACHED hexane film has no dry shell at all
      - the front is at the surface and ``f = 1``;
    * a row at or above the critical loading likewise has ``f = 1``, which is
      the factor-1.0 branch of the hexane closure;
    * the DRY_SORBATE and EXHAUSTED regimes are the fully receded core,
      ``f = 0``: the whole particle is shell;
    * otherwise the falling-rate law's own ``receding_front_fraction``.
    """

    if row.attached_hexane_kg != 0.0:
        return 1.0
    context = "D9-a front" if row.is_layer_aggregate else f"D9-a front (packet {row.packet_id})"
    # D9-2 (R-K1): the SAME three authorities the hexane closure reads, from
    # the same record, so the water side can never gate on a front the hexane
    # side did not compute.
    authorities = _hexane_authorities_at(
        dry_matter_kg=row.dry_matter_kg,
        internal_hexane_kg=row.internal_hexane_kg,
        residual_oil_label_kg=row.residual_oil_label_kg,
        temperature_k=row.temperature_k,
        pressure_pa=pressure_pa,
        context=context,
    )
    if authorities.core_is_fully_receded:
        return 0.0
    if authorities.front_is_at_the_surface:
        return 1.0
    try:
        return falling_rate.receding_front_fraction(
            authorities.internal_loading_kg_kg,
            authorities.critical_loading_kg_kg,
            authorities.equilibrium_floor_loading_kg_kg,
        )
    except falling_rate.FallingRateLawError as error:
        raise SP1K2Law2StepError(f"D9-a front closure refused: {error}") from error


def _dry_matter_weighted_front(
    rows: tuple[_PacketFrontRow, ...],
    *,
    pressure_pa: float,
    layer: int,
) -> float:
    """D19: ``sum_i m_i f_i / sum_i m_i``, exactly, over the rows given.

    A row with exactly zero dry matter carries zero weight in a
    dry-matter-weighted mean and has no loading to read, so it is left out of
    both sums rather than divided by; a row with NEGATIVE dry matter is
    inadmissible and refuses typed, and so does a set of rows carrying no
    positive dry matter at all.  Nothing is clamped and nothing is defaulted.
    """

    weights: list[float] = []
    fronts: list[float] = []
    for row in rows:
        mass = row.dry_matter_kg
        if mass < 0.0:
            raise SP1K2Law2StepError(
                f"K{layer} D9-a front: packet {row.packet_id!r} carries "
                f"{mass!r} kg of dry matter, which is not admissible"
            )
        if mass == 0.0:
            continue
        weights.append(mass)
        fronts.append(_packet_front_fraction(row, pressure_pa=pressure_pa, layer=layer))
    if not weights:
        raise SP1K2Law2StepError(
            f"K{layer} D9-a front: the layer's packet rows carry no positive dry matter"
        )
    first = fronts[0]
    if all(front == first for front in fronts):
        # A UNIFORM layer.  The weighted mean of one repeated value IS that
        # value, so it is returned rather than recomputed: m*f/m is not exact
        # in binary64, and the uniform-layer identity this rule is an
        # extension of has to hold to the bit, not to a tolerance.
        return first
    numerator = math.fsum(weight * front for weight, front in zip(weights, fronts, strict=True))
    return numerator / math.fsum(weights)


def _layer_front_fraction(
    inventory: PacketDerivedLayerInventory,
    *,
    pressure_pa: float = REGIME_C_LAYER_PRESSURE_PA,
) -> float:
    """D9-a-4: the hexane front fraction the WATER routing gates on.

    D19 (owner ruling 2026-09-23, option 1 of
    ``docs/GT_PS2_RULING_PACKET_D19_MIXED_LAYER_FRONT_2026-09-21.md``): the
    layer reports the DRY-MATTER-WEIGHTED front of its packets,
    ``f_layer = sum_i m_i f_i / sum_i m_i``, with ``f_i`` per packet by the
    three-authority read in :func:`_packet_front_fraction`.  The declared row
    ``D9A_FRONT_DECLARED_ASSUMPTIONS["dry_matter_weighted_packet_front"]``
    above owns the criterion, the basis, the measured bracket and the four
    rejecting outcomes.

    THE RULE UNTIL 2026-09-23, which this one extends and does not replace,
    was the layer AGGREGATE read: a layer whose aggregate
    ``inventory.attached_hexane_kg != 0.0`` returned 1.0 outright, and every
    other layer was read by the same three authorities at the layer's own
    AGGREGATE loading.  A uniform layer returns exactly what that read
    returned, bit for bit, and the D9-2 (R-K1) property that the water side
    and the hexane side consume one authority record is kept: both still read
    :func:`_hexane_authorities_at`, now once per packet.

    Why it had to move: measured on the production point (the packet's
    section 2), a layer of 26 fully-receded dry-sorbate packets and ONE
    3.27 g arrival sliver carrying 2.612e-4 kg of attached film reported a
    front of exactly 1.0 - a closed shell - at the instant its external water
    film drained to exactly zero, and the water arm then refused typed on a
    layer whose 6,229.9 kg of dry bulk had its shell fully open.  The sliver
    was read as the layer.
    """

    if inventory.dry_matter_kg <= 0.0:
        raise SP1K2Law2StepError("a K layer has no positive dry-matter holdup")
    return _dry_matter_weighted_front(
        _packet_front_rows(inventory),
        pressure_pa=pressure_pa,
        layer=inventory.layer,
    )


def _cell_inputs(
    *,
    model: SP1KLayerEngineeringModel,
    inventory: PacketDerivedLayerInventory,
    gas_boundary: ef.EngineeringGasBoundary,
    wall_temperature_k: float,
    macro_step_s: float,
    lagged_hexane_flux_kg_m2_s: float = 0.0,
    composed_class_channel: bool = False,
) -> ef.EngineeringCellInputs:
    if inventory.dry_matter_kg <= 0.0:
        raise SP1K2Law2StepError("a K layer has no positive dry-matter holdup")
    # B1 stage 3c / C8 Tier 1 / C8 Tier 2a: each arm is decided by its OWN
    # binder on its own exact condition, and both may bind at once.  The
    # hexane arm binds when the attached film is exactly zero (regime B:
    # a_h = 1 with the receding-front factor); the water arm binds when a
    # sorbed-water arm is declared AND the external water is exactly zero
    # (the Luikov activity with the declared series factor).  Both together
    # are the frozen PHY-007 DRY_SURFACE active set, which the cell admits
    # since 2026-09-05.  An exhausted hexane film beside a positive water
    # film is still EXTERNAL_WATER and the water binder still returns None
    # there, so every pre-existing pairing is unchanged.
    core_supplied_hexane = _core_supplied_hexane_closure(
        inventory, lagged_hexane_flux_kg_m2_s=lagged_hexane_flux_kg_m2_s
    )
    # Level 0b (iii): a COMPOSED class channel reads the free-film pin's binder;
    # every aggregate layer reads the unchanged binder and its K8 refusal.
    core_supplied_water = (
        _composed_class_water_closure(inventory, model.declared_sorbed_water_arm)
        if composed_class_channel
        else _core_supplied_water_closure(inventory, model.declared_sorbed_water_arm)
    )
    solid = cc.SolidSideAggregates(
        attached_hexane_loading_kg_kg=(inventory.attached_hexane_kg / inventory.dry_matter_kg),
        external_water_loading_kg_kg=(inventory.external_water_kg / inventory.dry_matter_kg),
        dry_matter_holdup_kg=inventory.dry_matter_kg,
        solid_temperature_k=inventory.solid_temperature_k,
        hexane_supply_to_film_kg_s=model.hexane_supply_to_film_kg_s,
        water_supply_to_film_kg_s=model.water_supply_to_film_kg_s,
    )
    return ef.EngineeringCellInputs(
        grid=model.grid,
        geometry=model.geometry,
        gas_boundary=gas_boundary,
        solid=solid,
        hydraulics=model.hydraulics,
        film_area_law=model.film_area_law,
        transfer=model.transfer,
        properties=model.properties,
        transport=model.transport,
        wall=model.wall,
        wall_temperature_k=wall_temperature_k,
        macro_step_s=macro_step_s,
        # B8: a declared lane band replaces the in-solve printed-hull gate
        # with the band gate applied at the same converged-drop site below.
        require_corroborated_pressure_drop=(model.declared_lane_drop_band_pa is None),
        # B1 stage 3c: None for every film-active layer (byte-identical);
        # the falling-rate closure for an exactly film-exhausted layer.
        core_supplied_hexane=core_supplied_hexane,
        # C8 Tier 1: None unless a sorbed-water arm is DECLARED on this model
        # AND the layer's external water is exactly zero.
        core_supplied_water=core_supplied_water,
    )


def _lane_band_checked(
    model: SP1KLayerEngineeringModel, solve: ef.BinaryNoInertCellSolve
) -> ef.BinaryNoInertCellSolve:
    """B8 declared lane gate on the converged layer drop (never an iterate).

    Applied exactly where the certified printed-hull gate would have fired;
    the printed hull stays authoritative for every corroboration CLAIM -
    the band only widens lawful ENGINEERING operation downward for drying
    beds under steam turndown, and the model digest carries the variant.
    """

    band = model.declared_lane_drop_band_pa
    if band is None:
        return solve
    drop = solve.evaluation.layer_pressure_drop_pa
    if not band[0] <= drop <= band[1]:
        raise SP1K2LaneDropBandError(
            f"the converged layer drop {drop!r} Pa left the declared "
            f"engineering-lane band [{band[0]}, {band[1]}] Pa "
            f"({DECLARED_LANE_DROP_BAND_LAW_ID})"
        )
    return solve


def _macro_step_from_solve(
    inputs: ef.EngineeringCellInputs, solve: ef.BinaryNoInertCellSolve
) -> ef.EngineeringCellMacroStep:
    accepted = ef.accept_engineering_fast_step(inputs, solve)
    wall = ef.advance_engineering_wall_node(accepted, wall=inputs.wall)
    ledger = ef.engineering_boundary_ledger(accepted, wall)
    return ef.EngineeringCellMacroStep(
        fast=solve,
        accepted=accepted,
        wall=wall,
        ledger=ledger,
        claims=ef.ENGINEERING_FEASIBILITY_CLAIMS,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class InterKGasFaceLedger:
    lower_outlet_molar_flow_mol_s: float
    upper_inlet_molar_flow_mol_s: float
    lower_outlet_hexane_mole_fraction: float
    upper_inlet_hexane_mole_fraction: float
    lower_outlet_temperature_k: float
    upper_inlet_temperature_k: float
    lower_downstream_pressure_pa: float
    upper_layer_pressure_pa: float
    pressure_residual_pa: float
    pressure_relative_residual: float
    pressure_relative_limit: float
    iterations: int

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.lower_outlet_molar_flow_mol_s == self.upper_inlet_molar_flow_mol_s
            and self.lower_outlet_hexane_mole_fraction == self.upper_inlet_hexane_mole_fraction
            and self.lower_outlet_temperature_k == self.upper_inlet_temperature_k
            and self.pressure_relative_residual <= self.pressure_relative_limit
        )


def _solve_serial_k_cells(
    request: SP1K2Law2StepRequest,
    inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory],
) -> tuple[
    tuple[ef.EngineeringCellInputs, ef.EngineeringCellInputs],
    tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep],
    InterKGasFaceLedger,
]:
    macro_step_s = request.end_time_s - request.prior.tray.time_s
    face_pressure = request.bottom_gas_boundary.downstream_boundary_pressure_pa
    # D9-e FIX (2026-09-14), declared row ``d9e_carried_seed_width_from_record``
    # (:data:`D9E_CARRIED_SEED_WIDTH_ROW`): these seeds are the PREVIOUS
    # interval's accepted cell field states, and since D9-e one of them may be
    # the ten-unknown dry-shell state.  The layer it belongs to may have left
    # the branch in between - regaining an attached hexane film, an external
    # water film, or both - so the width of the block that is actually solved
    # comes from the RECORD the cell inputs present, and the cell entry point
    # projects a carried ten-unknown seed onto its eight shipped components
    # when the record no longer presents the branch.  Nothing is carried
    # differently here; the rule lives at the one place that can see both.
    seeds = list(request.prior.cell_field_states)
    for iteration in range(1, request.face_max_iterations + 1):
        lower_boundary = replace(
            request.bottom_gas_boundary,
            downstream_boundary_pressure_pa=face_pressure,
        )
        lower_inputs = _cell_inputs(
            model=request.layer_models[0],
            inventory=inventories[0],
            gas_boundary=lower_boundary,
            wall_temperature_k=request.prior.wall_temperatures_k[0],
            macro_step_s=macro_step_s,
        )
        lower_solve = _lane_band_checked(
            request.layer_models[0],
            ef.solve_binary_no_inert_fast_block(
                lower_inputs,
                seed=seeds[0],
                tolerance=request.newton_tolerance,
            ),
        )
        seeds[0] = lower_solve.state
        upper_boundary = ef.EngineeringGasBoundary(
            inlet_molar_flow_mol_s=lower_solve.evaluation.outlet_molar_flow_mol_s,
            inlet_hexane_mole_fraction=lower_solve.state.gas_hexane_mole_fraction,
            inlet_water_mole_fraction=lower_solve.state.gas_water_mole_fraction,
            inlet_temperature_k=lower_solve.state.gas_temperature_k,
            downstream_boundary_pressure_pa=request.top_boundary_pressure_pa,
            carrier_topology=ef.CarrierTopology.BINARY_NO_INERT,
        )
        upper_inputs = _cell_inputs(
            model=request.layer_models[1],
            inventory=inventories[1],
            gas_boundary=upper_boundary,
            wall_temperature_k=request.prior.wall_temperatures_k[1],
            macro_step_s=macro_step_s,
        )
        upper_solve = _lane_band_checked(
            request.layer_models[1],
            ef.solve_binary_no_inert_fast_block(
                upper_inputs,
                seed=seeds[1],
                tolerance=request.newton_tolerance,
            ),
        )
        seeds[1] = upper_solve.state
        next_face_pressure = upper_solve.state.layer_pressure_pa
        pressure_residual = next_face_pressure - face_pressure
        pressure_relative = _relative_residual(pressure_residual, next_face_pressure, face_pressure)
        if pressure_relative <= request.face_tolerance:
            lower_step = _macro_step_from_solve(lower_inputs, lower_solve)
            upper_step = _macro_step_from_solve(upper_inputs, upper_solve)
            face = InterKGasFaceLedger(
                lower_outlet_molar_flow_mol_s=(lower_solve.evaluation.outlet_molar_flow_mol_s),
                upper_inlet_molar_flow_mol_s=upper_boundary.inlet_molar_flow_mol_s,
                lower_outlet_hexane_mole_fraction=(lower_solve.state.gas_hexane_mole_fraction),
                upper_inlet_hexane_mole_fraction=(upper_boundary.inlet_hexane_mole_fraction),
                lower_outlet_temperature_k=lower_solve.state.gas_temperature_k,
                upper_inlet_temperature_k=upper_boundary.inlet_temperature_k,
                lower_downstream_pressure_pa=face_pressure,
                upper_layer_pressure_pa=next_face_pressure,
                pressure_residual_pa=pressure_residual,
                pressure_relative_residual=pressure_relative,
                pressure_relative_limit=request.face_tolerance,
                iterations=iteration,
            )
            return (lower_inputs, upper_inputs), (lower_step, upper_step), face
        face_pressure = next_face_pressure
    raise SP1K2Law2StepError(
        "the serial K-cell internal pressure face did not converge without clipping"
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class AggregateSP1PressureLedger:
    lower_bed_drop_pa: float
    upper_bed_drop_pa: float
    physical_floor_drop_pa: float
    inactive_upper_floor_drop_pa: float
    active_physical_floor_count: int
    ordered_element_sum_pa: float
    layer_drop_sum_pa: float
    bottom_to_top_pressure_drop_pa: float
    lower_series_residual_pa: float
    upper_series_residual_pa: float
    aggregate_pressure_residual_pa: float
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.active_physical_floor_count == 1
            and self.physical_floor_drop_pa > 0.0
            and self.inactive_upper_floor_drop_pa == 0.0
            and self.maximum_relative_residual <= self.relative_limit
        )


def _aggregate_pressure_ledger(
    request: SP1K2Law2StepRequest,
    inputs: tuple[ef.EngineeringCellInputs, ef.EngineeringCellInputs],
    steps: tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep],
) -> AggregateSP1PressureLedger:
    drops = tuple(
        cc.layer_pressure_drop(
            series=cell_inputs.hydraulics.series,
            superficial_velocity_m_s=step.fast.evaluation.superficial_velocity_m_s,
            gas_density_kg_m3=step.fast.evaluation.gas_density_kg_m3,
        )
        for cell_inputs, step in zip(inputs, steps, strict=True)
    )
    lower_bed = drops[0].contribution_for(cc.BedErgunElement.label).drop_pa
    upper_bed = drops[1].contribution_for(cc.BedErgunElement.label).drop_pa
    lower_floor = drops[0].contribution_for(cc.TaperedBoreFloorPassageElement.label).drop_pa
    upper_floor = drops[1].contribution_for(cc.TaperedBoreFloorPassageElement.label).drop_pa
    ordered_sum = math.fsum((lower_bed, lower_floor, upper_bed, upper_floor))
    layer_sum = math.fsum(step.fast.evaluation.layer_pressure_drop_pa for step in steps)
    endpoint_drop = steps[0].fast.state.layer_pressure_pa - request.top_boundary_pressure_pa
    lower_residual = drops[0].drop_pa - steps[0].fast.evaluation.layer_pressure_drop_pa
    upper_residual = drops[1].drop_pa - steps[1].fast.evaluation.layer_pressure_drop_pa
    aggregate_residual = endpoint_drop - ordered_sum
    # W3-1 (owner ruling 2026-09-03, executed at HEAD f5a70bb).  PREVIOUSLY the
    # four folded residuals were divided by the tray pressure drop (~3.25e3 Pa)
    # and compared to 1.0e-10, an absolute allowance of ~3.25e-7 Pa - 9.302x
    # tighter than the 3.03e-6 Pa the two Newton pressure rows plus the internal
    # face Picard can deliver, so the gate refused solves that met the solver's
    # own contract at 73.4 % of the Newton gate.  The INTENT is unchanged:
    # pressure closure across the tray must still hold, now to the resolution
    # the solver actually promises.  `relative_limit` is untouched.
    resolution_bound = aggregate_pressure_binary64_resolution_bound(
        newton_tolerance=request.newton_tolerance,
        pressure_row_count=SP1_LAYER_COUNT,
        face_tolerance=request.face_tolerance,
        face_pressures_pa=(steps[1].fast.state.layer_pressure_pa,),
        summed_terms_pa=(
            lower_bed,
            lower_floor,
            upper_bed,
            upper_floor,
            ordered_sum,
            layer_sum,
            endpoint_drop,
        ),
    )
    maximum = max(
        _relative_residual(
            _residual_beyond_resolution(lower_residual, resolution_bound),
            drops[0].drop_pa,
            steps[0].fast.evaluation.layer_pressure_drop_pa,
        ),
        _relative_residual(
            _residual_beyond_resolution(upper_residual, resolution_bound),
            drops[1].drop_pa,
            steps[1].fast.evaluation.layer_pressure_drop_pa,
        ),
        _relative_residual(
            _residual_beyond_resolution(ordered_sum - layer_sum, resolution_bound),
            ordered_sum,
            layer_sum,
        ),
        _relative_residual(
            _residual_beyond_resolution(aggregate_residual, resolution_bound),
            endpoint_drop,
            ordered_sum,
        ),
    )
    return AggregateSP1PressureLedger(
        lower_bed_drop_pa=lower_bed,
        upper_bed_drop_pa=upper_bed,
        physical_floor_drop_pa=lower_floor,
        inactive_upper_floor_drop_pa=upper_floor,
        active_physical_floor_count=1,
        ordered_element_sum_pa=ordered_sum,
        layer_drop_sum_pa=layer_sum,
        bottom_to_top_pressure_drop_pa=endpoint_drop,
        lower_series_residual_pa=lower_residual,
        upper_series_residual_pa=upper_residual,
        aggregate_pressure_residual_pa=aggregate_residual,
        maximum_relative_residual=maximum,
        relative_limit=request.relative_limit,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class ClassHexaneTransfer:
    """Option A: one executable class's own accepted hexane transfer.

    ``class_key`` is the composing solve's executable class (a hexane class
    paired with the water state); ``hexane_class_key`` is the declared row's
    hexane class that class belongs to.  ``hexane_from_packets_kg`` is
    ``class_hexane_demand_kg`` of the class's own converged block,
    ``area * flux * M * dt`` in the kernel's own operation order, and the
    layer's ``hexane_from_packets_kg`` is the exact ``math.fsum`` of these rows
    (``_hexane_gate_rows`` refuses otherwise).
    """

    class_key: str
    hexane_class_key: str
    packet_ids: tuple[str, ...]
    gas_side_active_area_m2: float
    hexane_molar_flux_mol_m2_s: float
    macro_step_s: float
    hexane_from_packets_kg: float

    physically_qualifying: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class ClassWaterTransfer:
    """Level 0b: one executable class's own accepted WATER transfer.

    The water mirror of :class:`ClassHexaneTransfer` (owner ruling 2026-09-28,
    film-class water-routing packet option (a), component (ii)).
    ``water_from_packets_kg`` is :func:`class_water_demand_kg` of the class's own
    converged block, ``area * flux * M_w * dt`` in the kernel's own operation
    order, and the layer's ``water_from_packets_kg`` is the exact ``math.fsum``
    of these rows (``_water_gate_rows`` refuses otherwise).
    ``retained_water_kg_by_packet`` is each member packet's WEIGHTED retained
    water, aligned with ``packet_ids`` - the same weighted-totals entry
    ``_layer_inventory`` sums - carried because the class's sorbed stock needs
    it and the contributions carry no retained-water row; the kernel checks
    every entry against its own packets, bit for bit, and refuses otherwise.
    """

    class_key: str
    packet_ids: tuple[str, ...]
    retained_water_kg_by_packet: tuple[float, ...]
    gas_side_active_area_m2: float
    water_molar_flux_mol_m2_s: float
    macro_step_s: float
    water_from_packets_kg: float

    physically_qualifying: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerAcceptedTransfers:
    layer: int
    hexane_from_packets_kg: float
    water_from_packets_kg: float
    energy_from_packets_j: float
    gas_boundary_hexane_to_cell_kg: float
    gas_boundary_water_to_cell_kg: float
    gas_and_wall_boundary_energy_to_cell_j: float
    #: Option A: the per-class hexane transfers the layer total is the exact
    #: fsum of, supplied ONLY by a composing per-class solve of a layer holding
    #: more than one hexane class.  ``None`` - the default, and what
    #: ``_accepted_layer_transfers`` always returns - keeps the aggregate
    #: depletion gate and allocation byte for byte.
    hexane_by_class: tuple[ClassHexaneTransfer, ...] | None = None
    #: Level 0b: the per-class WATER transfers the layer total is the exact
    #: fsum of, supplied ONLY by a composing per-class solve.  ``None`` - the
    #: default, and what ``_accepted_layer_transfers`` always returns - keeps
    #: the layer water gate, allocation, drain and floor byte for byte.
    water_by_class: tuple[ClassWaterTransfer, ...] | None = None

    physically_qualifying: ClassVar[bool] = False


def _accepted_layer_transfers(
    inputs: ef.EngineeringCellInputs,
    step: ef.EngineeringCellMacroStep,
    *,
    layer: int,
    datum_adapter: component_datum.NumericComponentEnergyDatumAdapter,
) -> LayerAcceptedTransfers:
    evaluation = step.fast.evaluation
    state = step.fast.state
    dt = inputs.macro_step_s
    area = evaluation.gas_side_active_area_m2
    hexane_from_packets = (
        area * evaluation.hexane_molar_flux_mol_m2_s * ef.CANONICAL_HEXANE_MOLAR_MASS_KG_MOL * dt
    )
    water_from_packets = (
        area * evaluation.water_molar_flux_mol_m2_s * ef.CANONICAL_WATER_MOLAR_MASS_KG_MOL * dt
    )
    boundary = inputs.gas_boundary
    inlet_hexane = (
        boundary.inlet_molar_flow_mol_s
        * boundary.inlet_hexane_mole_fraction
        * ef.CANONICAL_HEXANE_MOLAR_MASS_KG_MOL
    )
    outlet_hexane = (
        evaluation.outlet_molar_flow_mol_s
        * state.gas_hexane_mole_fraction
        * ef.CANONICAL_HEXANE_MOLAR_MASS_KG_MOL
    )
    inlet_water = (
        boundary.inlet_molar_flow_mol_s
        * boundary.inlet_water_mole_fraction
        * ef.CANONICAL_WATER_MOLAR_MASS_KG_MOL
    )
    outlet_water = (
        evaluation.outlet_molar_flow_mol_s
        * state.gas_water_mole_fraction
        * ef.CANONICAL_WATER_MOLAR_MASS_KG_MOL
    )
    inlet_enthalpy = _stream_enthalpy_j_mol(
        inputs.properties,
        temperature_k=boundary.inlet_temperature_k,
        pressure_pa=state.layer_pressure_pa,
        hexane_mole_fraction=boundary.inlet_hexane_mole_fraction,
        water_mole_fraction=boundary.inlet_water_mole_fraction,
        datum_adapter=datum_adapter,
    )
    outlet_enthalpy = _stream_enthalpy_j_mol(
        inputs.properties,
        temperature_k=state.gas_temperature_k,
        pressure_pa=state.layer_pressure_pa,
        hexane_mole_fraction=state.gas_hexane_mole_fraction,
        water_mole_fraction=state.gas_water_mole_fraction,
        datum_adapter=datum_adapter,
    )
    gas_boundary_energy_w = (
        boundary.inlet_molar_flow_mol_s * inlet_enthalpy
        - evaluation.outlet_molar_flow_mol_s * outlet_enthalpy
    )
    external_energy = dt * (gas_boundary_energy_w + step.wall.steam_side_w - step.wall.to_ambient_w)
    # The host cell is the complementary gas-plus-wall control volume and its
    # quasi-steady gas has no storage.  Therefore its accepted energy change is
    # exactly the dynamic wall storage.  Solve the enclosing control-volume
    # identity for the packet contribution; do not invent a particle Cp law or
    # re-evaluate either shared-old-wall transfer.
    wall_storage = step.wall.storage_w * dt
    energy_from_packets = wall_storage - external_energy
    return LayerAcceptedTransfers(
        layer=layer,
        hexane_from_packets_kg=hexane_from_packets,
        water_from_packets_kg=water_from_packets,
        energy_from_packets_j=energy_from_packets,
        gas_boundary_hexane_to_cell_kg=dt * (inlet_hexane - outlet_hexane),
        gas_boundary_water_to_cell_kg=dt * (inlet_water - outlet_water),
        gas_and_wall_boundary_energy_to_cell_j=external_energy,
    )


def _phase_decrement_per_particle(
    *,
    total_to_cell: float,
    inventory: tray_host.PerParticleInventory,
    weighted_phase_total: float,
    weighted_dry_total: float,
    phase_value: float,
) -> float:
    if total_to_cell >= 0.0:
        if total_to_cell > weighted_phase_total:
            raise SP1K2Law2StepError(
                "Law-2 flux would deplete an external packet phase; step rejected, not clamped"
            )
        if weighted_phase_total == 0.0:
            if total_to_cell == 0.0:
                return 0.0
            raise SP1K2Law2StepError("positive transfer has no external packet inventory")
        return total_to_cell * phase_value / weighted_phase_total
    return total_to_cell * inventory.dry_matter_kg / weighted_dry_total


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerPacketExchangeLedger:
    """Decoded-candidate proof for one layer's packet mass/energy allocation."""

    layer: int
    aggregation_law_id: str
    energy_allocation_law_id: str
    declared_solid_to_interface_ua_w_k: float
    ua_share_sum_w_k: float
    ua_partition_residual_w_k: float
    packet_solid_to_interface_w: float
    cell_solid_to_interface_w: float
    thermal_moment_residual_w: float
    nonconductive_remainder_energy_j: float
    requested_hexane_transfer_kg: float
    realized_hexane_transfer_kg: float
    requested_water_transfer_kg: float
    realized_water_transfer_kg: float
    requested_energy_transfer_j: float
    realized_energy_transfer_j: float
    energy_transfer_residual_j: float
    energy_binary64_resolution_bound_j: float
    energy_transfer_residual_beyond_resolution_j: float
    energy_resolution_law_id: str
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    candidate_payloads_redecoded: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        return (
            self.aggregation_law_id == UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID
            and self.energy_allocation_law_id
            == PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID
            and self.energy_resolution_law_id == PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID
            and self.maximum_relative_residual <= self.relative_limit
        )


def _hexane_depletion_stock(layer_total) -> tuple[str, float]:
    """The drawable hexane stock of a populated layer, by its phase.

    * a film-bearing layer draws its attached film;
    * a receding-front layer (B1 stage 3c) draws the sorbed/pore inventory
      ABOVE the regime-C floor ``X_e(T)`` of the falling-rate law, the
      boundary that was refused rather than crossed until the tag-3 ruling;
    * a dry-sorbate layer (TAG-3, owner ruling 2026-09-06, Q-C8-2B) - mean
      loading at or below that floor - draws the WHOLE sorbed inventory: the
      GAB activity falls with the loading all the way to the EXHAUSTED regime,
      which the core-supplied binder refuses typed at exactly zero.  Before
      this branch the stock of such a layer read NEGATIVE (inventory minus the
      floor), and the first profile-start march was refused at step 1 on
      exactly that reading.

    The floor is the same call the falling-rate law and ``sorption_interface``
    share (the F6 instrument pins them equal), evaluated here as before.
    """

    if layer_total.attached_hexane_kg != 0.0:
        return "attached_hexane", layer_total.attached_hexane_kg
    floor_loading = falling_rate.equilibrium_floor_loading(
        layer_total.solid_temperature_k,
        layer_total.residual_oil_label_kg / layer_total.dry_matter_kg,
    )
    mean_loading = layer_total.internal_hexane_kg / layer_total.dry_matter_kg
    if mean_loading <= floor_loading:
        return "internal_hexane", layer_total.internal_hexane_kg
    return "internal_hexane", layer_total.internal_hexane_kg - (
        floor_loading * layer_total.dry_matter_kg
    )



#: OPTION A, level 0.  Owner ruling 2026-09-28 on option A of
#: docs/GT_PS2_RULING_PACKET_TAG2_FILM_EXHAUSTED_TRANSITION_2026-09-25.md as
#: corrected by docs/GT_PS2_OPTION_A_LEVEL0_DESIGN_2026-09-28.md (build record
#: docs/GT_PS2_OPTION_A_LEVEL0_BUILD_RECORD_2026-09-28.md).  The two
#: layer-AGGREGATE hexane reads - the arm binder's exact-zero film test and the
#: depletion stock above - become CLASS-aware as a per-class flux sum with a
#: per-class stock, never a weighted closure.  Neither body changes: the
#: binder and the stock law are read once per hexane CLASS.  Everything below
#: acts only on a layer whose accepted transfer carries ``hexane_by_class``,
#: which only a composing per-class solve supplies
#: (``study_lane/per_class_commit.py``); ``_accepted_layer_transfers`` never
#: sets it, so every other path is the pre-A code byte for byte.
#:
#: The one declared consequence the design names (its section 1): on this path
#: the D11 row ``per_class_film_floor_kg = 0.0`` is superseded - the B1 stage
#: 3c nil-film band applies to each hexane class's own attached aggregate, with
#: the unchanged NIL_FILM_RECLASSIFICATION_THRESHOLD_KG criterion, because a
#: class-level film has the same Zeno problem a layer-level one had.
OPTION_A_PER_CLASS_HEXANE_LAW_ID = "option-a-per-class-hexane-flux-and-stock-v1"
HEXANE_CLASS_TWO_FILM = "two_film"
HEXANE_CLASS_CORE_SUPPLIED_FRONT = "core_supplied_front"
HEXANE_CLASS_CORE_SUPPLIED_DRY_SORBATE = "core_supplied_dry_sorbate"
HEXANE_CLASS_ORDER = (
    HEXANE_CLASS_TWO_FILM,
    HEXANE_CLASS_CORE_SUPPLIED_FRONT,
    HEXANE_CLASS_CORE_SUPPLIED_DRY_SORBATE,
)
PER_CLASS_HEXANE_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "per_class_hexane_flux_and_stock": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": OPTION_A_PER_CLASS_HEXANE_LAW_ID,
                "declares": (
                    "the CLASS-RESOLVED hexane film flux and depletion stock of a layer, "
                    "replacing the two layer-AGGREGATE reads (the arm binder's "
                    "attached_hexane_kg != 0.0 on the layer, and the depletion stock's "
                    "attached_hexane_kg != 0.0 on the layer): a layer's hexane classes are "
                    "the sets of its packets on which the UNCHANGED binder takes the same "
                    "branch (two_film, core_supplied_front, core_supplied_dry_sorbate); "
                    "each class is solved on its own branch at its own state and draws on "
                    "its own stock"
                ),
                "value": (
                    "hexane_from_packets = sum_c A_c J_c M dt, an exact math.fsum over the "
                    "classes, each term A_c J_c M dt of the class's own converged block in "
                    "the kernel's own operation order; the depletion gate compares, per "
                    "class, the class's own demand against _hexane_depletion_stock(class "
                    "inventory); never a weighted closure - no activity, loading, front or "
                    "flux is averaged across classes"
                ),
                "criterion": (
                    "the sidecar measurement docs/evidence/sidecar_class_flux_2026-09-25/: on "
                    "every mixed layer of both F31 production-point lanes (MN1 layers 1 and "
                    "2; 66 and 80 scored rows) the flux the film-bearing class can supply "
                    "from its own film over its own area is 2.877e-07 to 9.003e-04 of the "
                    "aggregate flux the kernel charged (F31-001) and 9.615e-07 to 5.047e-04 "
                    "(F31-002), from the first scorable interval (step 2, t = 0.263 s, "
                    "ratio 3.762e-04); at the wall 5.054e-04 and 5.048e-04, the film's life "
                    "24.79 s against 12.5 ms; dry-matter and particle apportionments agree "
                    "to 2.7e-16, so the class share is not a choice"
                ),
                "no_new_number": (
                    "declares no constant: the classes are the binder's own branches, the "
                    "stock is the unchanged stock law, the nil-film band is the unchanged "
                    "NIL_FILM_RECLASSIFICATION_THRESHOLD_KG read on the CLASS aggregate, and "
                    "the class channel share is the D11 row per_class_channel_share "
                    "(representative-particle fraction)"
                ),
                "bracket": ["uniform-layer identity", "per-class sum"],
                "bracket_ends": (
                    "IDENTITY END: a layer of one class returns the layer inventory object "
                    "itself with share exactly 1.0, the composing solve bypasses, "
                    "hexane_by_class is None, and the gate, the stock, the allocation and "
                    "the demand are the pre-A code path byte for byte (tests: 20, 40 and "
                    "230 uniform committed states of the D19 lanes bit-identical; the "
                    "one-term demand reproduces F31-002's recorded transfer to the bit on 4 "
                    "of 4 stamped fixture layers).  PER-CLASS END: a mixed layer's demand is "
                    "sum_c A_c J_c M dt and its gate one row per class; measured at the wall "
                    "the per-class stocks are 0.2520 (film), 0.566 (front, above its floor) "
                    "and 3.152 kg (dry sorbate) on F31-001 L1 against the aggregate's 0.2520 "
                    "alone"
                ),
                "rejecting_outcome": (
                    "IDENTITY: any all-one-class inventory (all tag-1, all tag-2 or all "
                    "tag-3) whose gate row, stock, allocation or demand differs from the "
                    "pre-A kernel by one ulp.  ACCEPTANCE: the production point "
                    "op400s4pcd20 (seam ON, F31-001 argv) refuses at or before 3.161 plant s "
                    "for any reason.  CONSERVATION: any exact-fsum I2 or packet-water "
                    "closure of a committed interval of either acceptance lane above the "
                    "current 8e-12 kg (today 2.34e-12, 7.96e-12 and 7.53e-13 kg on MN1, MN2 "
                    "and SP1 of F31-001; 1.59e-13, 1.39e-13 and 1.10e-13 kg on the "
                    "packet-water ledgers), or any requested-versus-realized hexane residual "
                    "of a LayerPacketExchangeLedger leaving its pinned 1e-10 clause.  "
                    "FIXTURE: the fixture lane pcc360d20 refusing at any interval it "
                    "completed on the D20 tree (360 of 360 to 27.207 s).  ATTRIBUTION (the "
                    "packet's own): the per-class flux sum and the aggregate flux agreeing "
                    "on the mixed layers to within the 1e-10 ledger clause - already "
                    "measured NOT to fire (1.1e3 to 3.5e6 apart)"
                ),
                "sources": (
                    "docs/GT_PS2_RULING_PACKET_TAG2_FILM_EXHAUSTED_TRANSITION_2026-09-25.md "
                    "sections 2, 3, 8 (option A), 9, 10 and 13; "
                    "docs/evidence/sidecar_class_flux_2026-09-25/MEASUREMENT.md; "
                    "docs/evidence/farm_f31_2026-09-25/; the D19 pattern "
                    "docs/GT_PS2_D19_LEVEL0_BUILD_RECORD_2026-09-24.md; the D11 class "
                    "channel study_lane/per_class_interface.py "
                    "PER_CLASS_INTERFACE_DECLARATIONS"
                ),
            }
        ),
    }
)


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerHexaneClass:
    """Option A: one hexane class of a layer - the packets on which the
    UNCHANGED binder takes the same branch - with its own inventory."""

    class_key: str
    packet_ids: tuple[str, ...]
    inventory: PacketDerivedLayerInventory
    representative_particle_share: float
    is_whole_layer: bool

    physically_qualifying: ClassVar[bool] = False


def _packet_hexane_class_key(
    item: PacketLayerContribution, *, layer: int, pressure_pa: float
) -> str:
    """The branch the unchanged binder takes on a layer of this packet alone.

    Attached film on the packet: two-film (the binder returns ``None``).  A
    film-free packet is read by the same three authorities, in the same order,
    at its own state: a fully receded core (DRY_SORBATE or EXHAUSTED) is the
    dry-sorbate class, anything else the receding-front class.  Absent or
    negative dry matter has no loading to read and refuses typed.
    """

    dry = item.weighted_dry_matter_kg
    if not dry > 0.0:
        raise SP1K2Law2StepError(
            f"K{layer} option A: packet {item.packet_id!r} carries {dry!r} kg of dry "
            "matter; a hexane class needs positive dry matter - refuse, do not skip"
        )
    if item.weighted_attached_hexane_kg != 0.0:
        return HEXANE_CLASS_TWO_FILM
    authorities = _hexane_authorities_at(
        dry_matter_kg=dry,
        internal_hexane_kg=item.weighted_internal_hexane_kg,
        residual_oil_label_kg=item.weighted_residual_oil_label_kg,
        temperature_k=item.temperature_k,
        pressure_pa=pressure_pa,
        context=f"option A class (packet {item.packet_id})",
    )
    if authorities.core_is_fully_receded:
        return HEXANE_CLASS_CORE_SUPPLIED_DRY_SORBATE
    return HEXANE_CLASS_CORE_SUPPLIED_FRONT


def _class_sub_inventory(
    inventory: PacketDerivedLayerInventory,
    members: tuple[PacketLayerContribution, ...],
    *,
    retained_water_kg_by_packet=None,
) -> PacketDerivedLayerInventory:
    """A class's own inventory, built as ``_layer_inventory`` builds a layer's.

    Packet-id order (the members arrive in the layer's contribution order),
    ``math.fsum`` aggregates, the representative-particle temperature anchor,
    the layer's declared UA partitioned over the class's particles (the D11
    channel share, ``layer UA * class particles / layer particles``), and the
    B1 stage 3c nil-film band applied to the CLASS's own attached aggregate -
    the criterion unchanged, the object it reads the class.  The band's codec
    gate is met by construction wherever this runs on a mixed layer: a layer
    holding a film-free packet next to a film-bearing one carries a tag-2 or
    tag-3 payload, which only a codec speaking the subcritical variant can
    decode.  The external-water nil band is NOT re-applied per class (the
    design names only the hexane band).

    The contributions carry no retained-water row.  Where the caller supplies
    ``retained_water_kg_by_packet`` (the composing solve does, from the same
    weighted totals ``_layer_inventory`` reads) the class carries its exact
    ``math.fsum``; otherwise the field is NaN - a value every finiteness check
    refuses - rather than a number nobody computed.  The gate and the
    allocation below read only the hexane, dry-matter, oil and temperature
    fields.
    """

    counts = tuple(item.representative_particles for item in members)
    total = math.fsum(counts)
    share = total / math.fsum(item.representative_particles for item in inventory.contributions)
    temperatures = tuple(item.temperature_k for item in members)
    anchor = min(temperatures)
    solid_temperature = (
        anchor
        + math.fsum(
            count * (temperature - anchor)
            for count, temperature in zip(counts, temperatures, strict=True)
        )
        / total
    )
    declared_ua = inventory.solid_to_interface_ua_w_k * share
    ua_shares = tuple(declared_ua * count / total for count in counts)
    contributions = tuple(
        replace(item, ua_share_w_k=ua_share)
        for item, ua_share in zip(members, ua_shares, strict=True)
    )
    attached = math.fsum(item.weighted_attached_hexane_kg for item in members)
    internal = math.fsum(item.weighted_internal_hexane_kg for item in members)
    reclassified = 0.0
    if 0.0 < attached <= NIL_FILM_RECLASSIFICATION_THRESHOLD_KG:
        reclassified = attached
        internal = math.fsum((internal, attached))
        attached = 0.0
    if retained_water_kg_by_packet is None:
        retained = math.nan
    else:
        retained = math.fsum(retained_water_kg_by_packet[item.packet_id] for item in members)
    ua_share_sum = math.fsum(ua_shares)
    ua_moment = math.fsum(item.ua_share_w_k * item.temperature_k for item in contributions)
    return PacketDerivedLayerInventory(
        layer=inventory.layer,
        dry_matter_kg=math.fsum(item.weighted_dry_matter_kg for item in members),
        attached_hexane_kg=attached,
        internal_hexane_kg=internal,
        residual_oil_label_kg=math.fsum(item.weighted_residual_oil_label_kg for item in members),
        external_water_kg=math.fsum(item.weighted_external_water_kg for item in members),
        retained_water_kg=retained,
        nil_film_reclassified_kg=reclassified,
        common_datum_energy_j=math.fsum(item.weighted_common_datum_energy_j for item in members),
        solid_temperature_k=solid_temperature,
        solid_to_interface_ua_w_k=declared_ua,
        ua_share_sum_w_k=ua_share_sum,
        ua_partition_residual_w_k=ua_share_sum - declared_ua,
        ua_temperature_moment_residual_w=(ua_moment - declared_ua * solid_temperature),
        contributions=contributions,
        aggregation_law_id=inventory.aggregation_law_id,
        energy_allocation_law_id=inventory.energy_allocation_law_id,
        packet_count=len(members),
    )


def layer_hexane_classes(
    inventory: PacketDerivedLayerInventory,
    *,
    pressure_pa: float = REGIME_C_LAYER_PRESSURE_PA,
    retained_water_kg_by_packet=None,
) -> tuple[LayerHexaneClass, ...]:
    """Option A: the layer's hexane classes in HEXANE_CLASS_ORDER.

    A layer with no packet rows, or whose packets all fall in ONE class, is
    returned as ONE class whose inventory IS the layer inventory object, share
    exactly 1.0: the uniform-layer identity is by construction (the D19
    pattern), so every consumer of such a layer reads exactly what it read
    before.  Otherwise one sub-inventory per class.  Absent or negative dry
    matter on any packet refuses typed - every packet is checked before any is
    classified - and nothing is skipped or clamped.
    """

    if not inventory.dry_matter_kg > 0.0:
        raise SP1K2Law2StepError("a K layer has no positive dry-matter holdup")
    contributions = inventory.contributions
    for item in contributions:
        if not item.weighted_dry_matter_kg > 0.0:
            raise SP1K2Law2StepError(
                f"K{inventory.layer} option A: packet {item.packet_id!r} carries "
                f"{item.weighted_dry_matter_kg!r} kg of dry matter; a hexane class needs "
                "positive dry matter - refuse, do not skip"
            )
    if not contributions:
        key = _packet_hexane_class_key(
            PacketLayerContribution(
                packet_id="",
                representative_particles=1.0,
                temperature_k=inventory.solid_temperature_k,
                ua_share_w_k=inventory.solid_to_interface_ua_w_k,
                weighted_dry_matter_kg=inventory.dry_matter_kg,
                weighted_attached_hexane_kg=inventory.attached_hexane_kg,
                weighted_external_water_kg=inventory.external_water_kg,
                weighted_common_datum_energy_j=inventory.common_datum_energy_j,
                weighted_internal_hexane_kg=inventory.internal_hexane_kg,
                weighted_residual_oil_label_kg=inventory.residual_oil_label_kg,
            ),
            layer=inventory.layer,
            pressure_pa=pressure_pa,
        )
        return (
            LayerHexaneClass(
                class_key=key,
                packet_ids=(),
                inventory=inventory,
                representative_particle_share=1.0,
                is_whole_layer=True,
            ),
        )
    keyed: dict[str, list[PacketLayerContribution]] = {}
    for item in contributions:
        keyed.setdefault(
            _packet_hexane_class_key(item, layer=inventory.layer, pressure_pa=pressure_pa), []
        ).append(item)
    if len(keyed) == 1:
        (key,) = keyed
        return (
            LayerHexaneClass(
                class_key=key,
                packet_ids=tuple(item.packet_id for item in contributions),
                inventory=inventory,
                representative_particle_share=1.0,
                is_whole_layer=True,
            ),
        )
    total = math.fsum(item.representative_particles for item in contributions)
    return tuple(
        LayerHexaneClass(
            class_key=key,
            packet_ids=tuple(item.packet_id for item in keyed[key]),
            inventory=_class_sub_inventory(
                inventory,
                tuple(keyed[key]),
                retained_water_kg_by_packet=retained_water_kg_by_packet,
            ),
            representative_particle_share=(
                math.fsum(item.representative_particles for item in keyed[key]) / total
            ),
            is_whole_layer=False,
        )
        for key in HEXANE_CLASS_ORDER
        if key in keyed
    )


def class_hexane_closures(
    classes: tuple[LayerHexaneClass, ...],
    *,
    lagged_hexane_flux_kg_m2_s: float = 0.0,
    pressure_pa: float = REGIME_C_LAYER_PRESSURE_PA,
) -> tuple[ef.CoreSuppliedHexaneClosure | None, ...]:
    """The UNCHANGED binder, read once per class - never on a weighted mean."""

    return tuple(
        _core_supplied_hexane_closure(
            item.inventory,
            lagged_hexane_flux_kg_m2_s=lagged_hexane_flux_kg_m2_s,
            pressure_pa=pressure_pa,
        )
        for item in classes
    )


def hexane_depletion_stocks_by_class(
    classes: tuple[LayerHexaneClass, ...],
) -> tuple[tuple[str, str, float], ...]:
    """``(class_key, phase_name, stock)`` - the UNCHANGED stock law per class."""

    return tuple((item.class_key, *_hexane_depletion_stock(item.inventory)) for item in classes)


def class_hexane_demand_kg(
    gas_side_active_area_m2: float, hexane_molar_flux_mol_m2_s: float, macro_step_s: float
) -> float:
    """``A_c J_c M dt`` in ``_accepted_layer_transfers``'s own operation order."""

    return (
        gas_side_active_area_m2
        * hexane_molar_flux_mol_m2_s
        * ef.CANONICAL_HEXANE_MOLAR_MASS_KG_MOL
        * macro_step_s
    )


def _option_a_partition_refusal(layer: int, what: str) -> SP1K2Law2StepError:
    return SP1K2Law2StepError(
        f"K{layer} option A: the composing solve's classes do not partition the layer's "
        f"hexane classes exactly ({what}) - refuse, do not remap"
    )


def _hexane_gate_rows(
    layer_total: PacketDerivedLayerInventory,
    transfer: LayerAcceptedTransfers,
) -> tuple[tuple[str | None, str, float, float], ...]:
    """``(hexane_class_key, phase_name, demanded, stock)`` rows of the gate.

    ``transfer.hexane_by_class is None``: ONE row, the aggregate read,
    ``(None, *_hexane_depletion_stock(layer_total))`` against the layer demand -
    byte for byte the pre-A gate.  Otherwise one row per hexane class present,
    in HEXANE_CLASS_ORDER: its demand is the ``math.fsum`` of its executable
    classes' own transfers, its stock the unchanged stock law on its own
    sub-inventory.  The layer demand must be the exact ``math.fsum`` of the
    class transfers, and the composing solve's classes must partition the
    layer's packets exactly as ``layer_hexane_classes`` does; either miss
    refuses typed.
    """

    if transfer.hexane_by_class is None:
        phase, stock = _hexane_depletion_stock(layer_total)
        return ((None, phase, transfer.hexane_from_packets_kg, stock),)
    rows = transfer.hexane_by_class
    layer = layer_total.layer
    if transfer.hexane_from_packets_kg != math.fsum(r.hexane_from_packets_kg for r in rows):
        raise _option_a_partition_refusal(
            layer, "the layer's hexane transfer is not the exact fsum of its class transfers"
        )
    classes = layer_hexane_classes(layer_total)
    known = {item.class_key for item in classes}
    if any(r.hexane_class_key not in known for r in rows):
        raise _option_a_partition_refusal(layer, "a class names a hexane class the layer lacks")
    if sorted(pid for r in rows for pid in r.packet_ids) != sorted(
        pid for item in classes for pid in item.packet_ids
    ):
        raise _option_a_partition_refusal(layer, "the class packets are not the layer's packets")
    gate: list[tuple[str | None, str, float, float]] = []
    for item in classes:
        members = [r for r in rows if r.hexane_class_key == item.class_key]
        if sorted(pid for r in members for pid in r.packet_ids) != sorted(item.packet_ids):
            raise _option_a_partition_refusal(layer, f"hexane class {item.class_key!r}")
        phase, stock = _hexane_depletion_stock(item.inventory)
        gate.append(
            (
                item.class_key,
                phase,
                math.fsum(r.hexane_from_packets_kg for r in members),
                stock,
            )
        )
    return tuple(gate)


#: LEVEL 0b, component (ii) of option (a), with component (i) on the same
#: layers.  Owner ruling 2026-09-28 on
#: docs/GT_PS2_RULING_PACKET_FILM_CLASS_WATER_ROUTING_2026-09-28.md section 5
#: ("(a) as recommended").  The water mirror of option A's class-resolved
#: hexane flux: on a layer whose accepted transfer carries ``water_by_class``
#: (supplied only by a composing per-class solve), each executable class books
#: its OWN water demand to its OWN packets, gates it against its OWN stock, and
#: the D9-a drain and floor read each packet's OWN front.  Everything below acts
#: only there; ``_accepted_layer_transfers`` never sets ``water_by_class``, so
#: every other layer - every layer of every lane with the seam off - takes the
#: pre-0b code byte for byte.
LEVEL0B_CLASS_OWN_WATER_ROUTING_LAW_ID = "level0b-class-own-water-routing-v1"
PER_CLASS_WATER_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "class_own_water_routing": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": LEVEL0B_CLASS_OWN_WATER_ROUTING_LAW_ID,
                "declares": (
                    "on a layer composed per class, the CLASS-RESOLVED water booking and "
                    "the PACKET-OWN D9-a front: (ii) the layer's water transfer is the "
                    "exact math.fsum over its executable classes of each class's own "
                    "A_c J_w,c M_w dt; the depletion gate reads one row per class, the "
                    "class's own demand against the class's own stock (its external "
                    "film where the class carries one, else its exchangeable sorbed "
                    "water summed over its packets at each packet's own front); each "
                    "packet draws on its class's demand on its class's basis (the "
                    "external-film basis if the class carries a film, the dry-matter "
                    "basis if it does not); (i) the D9-a-2 drain and the D9-a-3 floor "
                    "read each packet's own three-authority front, read on the class "
                    "sub-inventory's rows, instead of the layer's dry-matter-weighted "
                    "mean"
                ),
                "value": (
                    "water_from_packets = sum_c A_c J_w,c M_w dt, an exact math.fsum, "
                    "each term in _accepted_layer_transfers' own operation order; "
                    "stock_c = external_water_c if external_water_c != 0.0 or no arm is "
                    "declared, else fsum_i exchangeable_water_kg(m_i, W_i, f_i) over the "
                    "class's packets; f_i = _packet_front_fraction on the class "
                    "sub-inventory's row of packet i; never a weighted closure - no "
                    "activity, loading, front or flux is averaged across classes"
                ),
                "criterion": (
                    "measured on the option A level-0 smokes (routing packet sections "
                    "1.3 and 1.4, docs/evidence/option_a_level0_2026-09-28/): at the "
                    "layer's weighted front (about 3.8e-4) the D9-a drain moved a "
                    "film-bearing arrival sliver's whole 0.02340675 kg water film into "
                    "the sliver's own sorbed store at step 2 - the capacity it read, "
                    "0.3511 kg, is fifteen times the film, where the sliver's own front "
                    "of exactly 1 gives a capacity of exactly 0.0 - and the layer-level "
                    "dry-matter basis booked each film-bearing sliver its 3.76e-4 share "
                    "of the layer's condensate (8.122e-5 kg at step 3 of the pre-edit "
                    "twin) into the sorbed store behind its hexane film; D9-a-1 says a "
                    "closed shell cannot take that water in.  The layer-level booking "
                    "is the water side of the misattribution option A removed for "
                    "hexane"
                ),
                "no_new_number": (
                    "declares no constant: the classes are the composing solve's "
                    "executable classes (hexane class by water state), the stock is the "
                    "unchanged D9-a-3 exchangeable-water law, the front is the unchanged "
                    "three-authority read, the bases are the kernel's existing two"
                ),
                "bracket": ["uniform-layer identity", "per-class sum"],
                "bracket_ends": (
                    "IDENTITY END: a layer not composed per class carries "
                    "water_by_class None, and its gate, allocation, drain, floor and "
                    "binder are the pre-0b code path byte for byte (the weighted front "
                    "of a uniform layer IS every packet's own front, bit for bit, by the "
                    "D19 uniform-layer identity).  PER-CLASS END: a composed layer's "
                    "water is sum_c A_c J_w,c M_w dt booked class by class"
                ),
                "scope": (
                    "composed layers only.  Component (i) is NOT applied to a layer the "
                    "seam did not compose, because there the layer is still solved and "
                    "booked as one aggregate: a film-bearing sliver keeping its water "
                    "film would make the layer's aggregate external water positive and "
                    "the aggregate binder would then solve the whole layer on the "
                    "two-film law, which is the harm the routing packet's option (d1) "
                    "predicts for (i) without (ii)"
                ),
                "rejecting_outcome": (
                    "ACCEPTANCE: op400s4pcd20 (F31-001 argv, sorbed token) refuses at "
                    "or before 3.161 plant s.  CONSERVATION: a committed interval whose "
                    "I2 or packet-water closure exceeds its pre-edit level, or a "
                    "composed interval whose tray gas-water, tray gas-hexane or tower "
                    "two-path residuals exceed the 1e-10 relative ledger clause.  "
                    "IDENTITY: a uniform layer whose drain, floor, water allocation, "
                    "binder result or outlet differs by one ulp from the D19 path; "
                    "pcc360's first three intervals differing from the level-0 twin.  "
                    "ATTRIBUTION: water booked to a class's packets differing from that "
                    "class channel's demand by more than the ledger clause.  FIXTURE: a "
                    "mixed-layer fixture of tests/fixtures/option_a_mixed_layers/ whose "
                    "per-class hexane result moves.  EVAPORATING ZERO FILM: a "
                    "film-bearing class at zero water film whose class solve evaporates "
                    "water - its stock is exactly 0.0 and the gate refuses typed.  D23 "
                    "CLOSURE: a composed interval whose tray gas-side water or hexane "
                    "exchange, or outlet energy, differs from the sum over class "
                    "channels by more than the ledger clause"
                ),
                "sources": (
                    "docs/GT_PS2_RULING_PACKET_FILM_CLASS_WATER_ROUTING_2026-09-28.md "
                    "sections 1 to 5 (option (a), components (i) to (iii)); "
                    "docs/GT_PS2_D23_SEAM_UPPER_OUTLET_GAP_RECORD_2026-09-28.md; the D9 "
                    "design brief row D9-a-2 ('f = 1 while the layer carries an attached "
                    "hexane film (no shell; capacity 0)'); the option A pattern "
                    "PER_CLASS_HEXANE_DECLARED_ASSUMPTIONS above"
                ),
            }
        ),
    }
)


#: LEVEL 0b, component (iii): the free-film pin on a composed, film-bearing,
#: water-free class channel.  Owner ruling 2026-09-28, "(a) as recommended".
LEVEL0B_FREE_FILM_PIN_LAW_ID = "level0b-free-film-pin-v1"
#: The pin's water activity and series factor.  NOT a new number: the arm's
#: own ``activity()`` returns exactly 1.0 at and above x_bf, and 1.0 is the
#: factor of an unhindered film - the two values the cell already solves a
#: hexane-film layer with when the arm is bound.
FREE_FILM_PIN_WATER_ACTIVITY = 1.0
FREE_FILM_PIN_WATER_CONDUCTANCE_FACTOR = 1.0
FREE_FILM_PIN_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "free_film_pin_on_a_composed_film_class": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": LEVEL0B_FREE_FILM_PIN_LAW_ID,
                "declares": (
                    "on a COMPOSED class channel only (a class sub-inventory the per-class "
                    "commit seam solves as its own channel), a film-bearing class - "
                    "attached hexane film present, external water exactly 0.0, a declared "
                    "sorbed-water arm, front at the declared no-shell ceiling - binds the "
                    "water arm at a_w = 1.0 with factor 1.0 instead of raising "
                    "SP1K2SorbedWaterNoShellError; condensate landing on it is booked to "
                    "its own external film (D9-a-1) and stays there (component (i): its "
                    "own front is exactly 1, its drain capacity exactly 0.0)"
                ),
                "value": "water_activity = 1.0, water_conductance_factor = 1.0",
                "criterion": (
                    "D9-a-1: liquid water at the surface sets the surface activity to one "
                    "and a closed shell cannot take it in; the arm's own activity() is "
                    "exactly 1.0 at and above x_bf; the certified retained-upper closure "
                    "(GT_PS2_RETAINED_UPPER_GATE_CLOSURE_RECORD_2026-08-21.md) resolved the "
                    "same exact-zero-water condensing state by lawful water appearance; the "
                    "option A level-0 smoke refused at step 3 (0.5 plant s) on exactly this "
                    "state, a 2.34 kg arrival sliver whose water film the layer-front drain "
                    "had emptied (routing packet sections 1.1 to 1.3)"
                ),
                "no_new_number": (
                    "declares no constant: 1.0 is the arm's own activity at the film-onset "
                    "loading and the unhindered factor; the ceiling it reads is the unchanged "
                    "NO_SHELL_FRONT_FRACTION_CEILING"
                ),
                "bracket": ["refusal (K8)", "free-film pin"],
                "scope": (
                    "composed class channels only: _cell_inputs(composed_class_channel=True), "
                    "which only the per-class commit seam passes; the aggregate-layer binder "
                    "and its K8 refusal are untouched on every aggregate path"
                ),
                "rejecting_outcome": (
                    "EVAPORATING ZERO FILM: a pinned class whose class solve evaporates water "
                    "has no stock (its exchangeable sorbed water is exactly 0.0 at front 1 and "
                    "its film is 0.0), the unchanged class gate refuses typed, and the "
                    "one-sided law (a film if condensing, blocked if not) goes to the owner; "
                    "ACCEPTANCE, CONSERVATION, IDENTITY and ATTRIBUTION as the class-own "
                    "water routing row words them"
                ),
                "sources": (
                    "docs/GT_PS2_RULING_PACKET_FILM_CLASS_WATER_ROUTING_2026-09-28.md "
                    "section 2(a) component (iii) and section 3; the D9 owner rulings "
                    "D9-a-1 and D9-b-2 (K8)"
                ),
            }
        ),
    }
)


def _composed_class_water_closure(
    inventory: PacketDerivedLayerInventory,
    arm: sorbed_water.SorbedWaterArm | None,
) -> ef.CoreSuppliedWaterClosure | None:
    """Level 0b (iii): the water binder of a COMPOSED class channel.

    Exactly ``_core_supplied_water_closure`` except on the one state K8 would
    refuse on a film-bearing class: a declared arm, external water exactly 0.0,
    an attached hexane film, and a front at or above the declared no-shell
    ceiling.  There the free-film pin binds a_w = 1.0 with factor 1.0.  A
    film-FREE class at a closed shell (a wet core at the critical loading, or
    a class whose film the class-level nil band reclassified) still refuses.
    """

    if (
        arm is not None
        and inventory.external_water_kg == 0.0
        and inventory.attached_hexane_kg != 0.0
        and _layer_front_fraction(inventory) >= sorbed_water.NO_SHELL_FRONT_FRACTION_CEILING
    ):
        return ef.CoreSuppliedWaterClosure(
            water_activity=FREE_FILM_PIN_WATER_ACTIVITY,
            water_conductance_factor=FREE_FILM_PIN_WATER_CONDUCTANCE_FACTOR,
        )
    return _core_supplied_water_closure(inventory, arm)


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerWaterClass:
    """Level 0b: one executable water class of a composed layer, read once.

    ``inventory`` is the class sub-inventory (``_class_sub_inventory``, retained
    water carried), ``own_front_by_packet`` each member's own three-authority
    front read on that inventory's rows, and ``(phase_name, stock_kg)`` the
    class's gate row against ``demand_kg``, the fsum-exact class transfer.
    """

    class_key: str
    packet_ids: tuple[str, ...]
    inventory: PacketDerivedLayerInventory
    own_front_by_packet: MappingProxyType
    sorbed_arm_active: bool
    phase_name: str
    stock_kg: float
    demand_kg: float
    demand_rate_kg_s: float

    physically_qualifying: ClassVar[bool] = False


def class_water_demand_kg(
    gas_side_active_area_m2: float, water_molar_flux_mol_m2_s: float, macro_step_s: float
) -> float:
    """``A_c J_w,c M_w dt`` in ``_accepted_layer_transfers``'s own operation order."""

    return (
        gas_side_active_area_m2
        * water_molar_flux_mol_m2_s
        * ef.CANONICAL_WATER_MOLAR_MASS_KG_MOL
        * macro_step_s
    )


def _level0b_partition_refusal(layer: int, what: str) -> SP1K2Law2StepError:
    return SP1K2Law2StepError(
        f"K{layer} level 0b: the composing solve's water classes do not partition the "
        f"layer's packets exactly ({what}) - refuse, do not remap"
    )


def layer_water_classes(
    layer_total: PacketDerivedLayerInventory,
    transfer: LayerAcceptedTransfers,
    arm: sorbed_water.SorbedWaterArm | None,
    *,
    pressure_pa: float,
    retained_water_kg_by_packet_check=None,
) -> tuple[LayerWaterClass, ...]:
    """Level 0b: the composed layer's water classes, in the transfer's row order.

    Refuses typed - never remaps - when the rows do not partition the layer's
    packets exactly, when the layer's water transfer is not the exact fsum of
    the class rows, when a row carries fewer retained-water entries than
    packets, when a member has no positive dry matter, or (kernel only, where
    ``retained_water_kg_by_packet_check`` is the packets' own weighted retained
    water) when a carried retained-water entry is not the packet's own, bit for
    bit.
    """

    rows = transfer.water_by_class
    layer = layer_total.layer
    if rows is None:
        raise _level0b_partition_refusal(layer, "the transfer carries no class water rows")
    if len(rows) < 2:
        raise _level0b_partition_refusal(layer, "a composed water booking needs two classes")
    if transfer.water_from_packets_kg != math.fsum(r.water_from_packets_kg for r in rows):
        raise _level0b_partition_refusal(
            layer, "the layer's water transfer is not the exact fsum of its class transfers"
        )
    contribution_ids = [item.packet_id for item in layer_total.contributions]
    if sorted(pid for r in rows for pid in r.packet_ids) != sorted(contribution_ids):
        raise _level0b_partition_refusal(layer, "the class packets are not the layer's packets")
    classes: list[LayerWaterClass] = []
    for r in rows:
        if len(r.retained_water_kg_by_packet) != len(r.packet_ids):
            raise _level0b_partition_refusal(
                layer, f"class {r.class_key!r} carries a retained-water row per packet"
            )
        retained = dict(zip(r.packet_ids, r.retained_water_kg_by_packet, strict=True))
        if retained_water_kg_by_packet_check is not None:
            for packet_id, value in retained.items():
                if struct.pack(">d", value) != struct.pack(
                    ">d", retained_water_kg_by_packet_check[packet_id]
                ):
                    raise _level0b_partition_refusal(
                        layer, f"packet {packet_id!r}'s carried retained water is not its own"
                    )
        members = tuple(item for item in layer_total.contributions if item.packet_id in retained)
        for item in members:
            if not item.weighted_dry_matter_kg > 0.0:
                raise SP1K2Law2StepError(
                    f"K{layer} level 0b: packet {item.packet_id!r} carries "
                    f"{item.weighted_dry_matter_kg!r} kg of dry matter; a water class needs "
                    "positive dry matter - refuse, do not skip"
                )
        inventory = _class_sub_inventory(layer_total, members, retained_water_kg_by_packet=retained)
        fronts = {
            row.packet_id: _packet_front_fraction(row, pressure_pa=pressure_pa, layer=layer)
            for row in _packet_front_rows(inventory)
        }
        external = inventory.external_water_kg
        if external != 0.0 or arm is None:
            phase_name, stock = "external_water", external
        else:
            phase_name = "retained_water"
            stock = math.fsum(
                sorbed_water.exchangeable_water_kg(
                    dry_matter_kg=item.weighted_dry_matter_kg,
                    retained_water_kg=retained[item.packet_id],
                    front_fraction=fronts[item.packet_id],
                )
                for item in members
            )
        classes.append(
            LayerWaterClass(
                class_key=r.class_key,
                packet_ids=tuple(item.packet_id for item in members),
                inventory=inventory,
                own_front_by_packet=MappingProxyType(fronts),
                sorbed_arm_active=external == 0.0 and arm is not None,
                phase_name=phase_name,
                stock_kg=stock,
                demand_kg=r.water_from_packets_kg,
                demand_rate_kg_s=(
                    r.gas_side_active_area_m2
                    * r.water_molar_flux_mol_m2_s
                    * ef.CANONICAL_WATER_MOLAR_MASS_KG_MOL
                ),
            )
        )
    return tuple(classes)


def _layer_sorbed_arm(request, layer: int) -> sorbed_water.SorbedWaterArm | None:
    """The DECLARED sorbed-water arm of one layer's model, or ``None``.

    Read from the request's own layer models, which the through-bed host
    shares with this kernel; ``None`` is every certified run.
    """

    models = getattr(request, "layer_models", None)
    if models is None or not 1 <= layer <= len(models):
        return None
    return models[layer - 1].declared_sorbed_water_arm


def _packet_entries_after(
    request: SP1K2Law2StepRequest,
    transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers],
    inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory],
    macro_steps: tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep],
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[
    dict[str, tray_host.PacketPlanEntry],
    tuple[LayerPacketExchangeLedger, LayerPacketExchangeLedger],
]:
    by_layer = {item.layer: item for item in transfers}
    totals = {item.layer: item for item in inventories}
    steps = {layer: macro_steps[layer - 1] for layer in (1, 2)}
    contribution_by_id = {
        contribution.packet_id: contribution
        for inventory in inventories
        for contribution in inventory.contributions
    }
    remainder_by_layer: dict[int, float] = {}
    # D9-a-4: ONE front fraction per layer per interval, at the accepted layer
    # pressure the cell solved with, read by the drain and by the arm's floor.
    front_by_layer: dict[int, float] = {}
    for layer in (1, 2):
        layer_total = totals[layer]
        step = steps[layer]
        front_by_layer[layer] = _layer_front_fraction(
            layer_total, pressure_pa=step.fast.state.layer_pressure_pa
        )
        conductive_energy = (
            step.accepted.macro_step_s
            * layer_total.solid_to_interface_ua_w_k
            * (layer_total.solid_temperature_k - step.fast.state.interface_temperature_k)
        )
        remainder_by_layer[layer] = by_layer[layer].energy_from_packets_j - conductive_energy
    # B1 stage 1 (owner-signed): the depletion condition below is
    # layer-aggregate and packet-independent, so it is checked once per
    # populated layer here — in the exact order the packet loop would
    # encounter it (canonical packet order: layer ascending, hexane
    # before water) — and raised TYPED with the transition data.  The
    # per-particle helper keeps the identical check as defense in depth;
    # for every reachable input the refusal outcome and message are
    # unchanged.
    populated_layers = sorted(
        {packet.owner_key.vertical_layer_id.value for packet in request.prior.tray.packets}
    )
    # Level 0b: each composed layer's water classes, read ONCE, with each
    # member's own front at the accepted layer pressure the cell solved with and
    # the carried retained water checked against the packets' own, bit for bit.
    # Empty for every layer whose transfer carries no ``water_by_class``.
    water_classes_by_layer: dict[int, tuple[LayerWaterClass, ...]] = {}
    if any(by_layer[layer].water_by_class is not None for layer in populated_layers):
        own_retained = {
            packet.packet_id: packet.weighted_totals.totals[5]
            for packet in request.prior.tray.packets
        }
        for layer in populated_layers:
            if by_layer[layer].water_by_class is not None:
                water_classes_by_layer[layer] = layer_water_classes(
                    totals[layer],
                    by_layer[layer],
                    _layer_sorbed_arm(request, layer),
                    pressure_pa=steps[layer].fast.state.layer_pressure_pa,
                    retained_water_kg_by_packet_check=own_retained,
                )
    for layer in populated_layers:
        transfer = by_layer[layer]
        layer_total = totals[layer]
        # Option A: one aggregate row (byte for byte the pre-A gate) unless the
        # layer's transfer was composed per hexane class; then one row per class.
        hexane_rows = _hexane_gate_rows(layer_total, transfer)
        if layer in water_classes_by_layer:
            # Level 0b: one water row per class - its own demand against its
            # own stock - in the composing solve's class order.
            water_rows = tuple(
                (item.class_key, item.phase_name, item.demand_kg, item.stock_kg)
                for item in water_classes_by_layer[layer]
            )
        elif layer_total.external_water_kg != 0.0 or _layer_sorbed_arm(request, layer) is None:
            water_phase_name = "external_water"
            water_stock = layer_total.external_water_kg
        else:
            # C8 Tier 1: on the sorbed arm the drawable stock is the retained
            # water ABOVE the accepted law's qualified floor.  Sorbed water may
            # evaporate (the interface can sit above the co-boiling lock) down
            # to that floor, and below it the step is refused typed - the dry
            # continuation is Tier 2, not a clamp.
            water_phase_name = "retained_water"
            water_stock = layer_total.retained_water_kg - (
                sorbed_water.LUIKOV_PARAMS.W_ref * layer_total.dry_matter_kg
            )
        if layer not in water_classes_by_layer:
            water_rows = ((None, water_phase_name, transfer.water_from_packets_kg, water_stock),)
        for is_water, (class_key, phase_name, demanded, stock) in (
            *((False, row) for row in hexane_rows),
            *((True, row) for row in water_rows),
        ):
            if demanded >= 0.0 and demanded > stock:
                raise SP1K2Law2PhaseExhaustionError(
                    "Law-2 flux would deplete an external packet phase; "
                    "step rejected, not clamped",
                    vertical_layer_id=layer,
                    phase_name=phase_name,
                    layer_stock_kg=stock,
                    demanded_kg=demanded,
                    layer_macro_step_s=steps[layer].accepted.macro_step_s,
                    interval_end_time_s=getattr(request, "end_time_s", None),
                    hexane_class_key=None if is_water else class_key,
                    water_class_key=class_key if is_water else None,
                )
    # Option A: each packet's hexane class, read ONCE per composed layer by the
    # same partition the gate above used; empty for every aggregate layer.
    class_of_packet: dict[str, LayerHexaneClass] = {
        packet_id: item
        for layer in populated_layers
        if by_layer[layer].hexane_by_class is not None
        for item in layer_hexane_classes(totals[layer])
        for packet_id in item.packet_ids
    }
    # Level 0b: each packet's water class on a composed layer; empty elsewhere.
    water_class_of: dict[str, LayerWaterClass] = {
        packet_id: item
        for classes in water_classes_by_layer.values()
        for item in classes
        for packet_id in item.packet_ids
    }
    entries: dict[str, tray_host.PacketPlanEntry] = {}
    realized_rows: dict[int, list[tuple[float, float, float]]] = {1: [], 2: []}
    energy_resolution_rows: dict[int, list[float]] = {1: [], 2: []}
    for packet in request.prior.tray.packets:
        layer = packet.owner_key.vertical_layer_id.value
        transfer = by_layer[layer]
        layer_total = totals[layer]
        inventory = packet.entry.inventory
        contribution = contribution_by_id[packet.packet_id]
        hexane_total_to_cell = transfer.hexane_from_packets_kg
        if transfer.hexane_by_class is not None:
            # Option A: the packet draws on ITS class's own demand and stock, by
            # the same proportional law the layer used, on the class basis.
            own = class_of_packet[packet.packet_id]
            hexane_total_to_cell = math.fsum(
                r.hexane_from_packets_kg
                for r in transfer.hexane_by_class
                if r.hexane_class_key == own.class_key
            )
            basis = own.inventory
        else:
            basis = layer_total
        if basis.attached_hexane_kg != 0.0:
            hexane_basis_total = basis.attached_hexane_kg
            hexane_basis_value = inventory.attached_hexane_kg
        else:
            # B1 stage 3c: film-exhausted layer - allocate the core draw in
            # proportion to each packet's TOTAL hexane (the same proportional
            # law the film draw uses, on the arm's supply).  The per-packet
            # total matches the reclassified layer basis whether the packet
            # still carries a nil film remnant or crossed already, keeping
            # the exchange-ledger residual at machine resolution.
            hexane_basis_total = basis.internal_hexane_kg
            hexane_basis_value = math.fsum(
                (inventory.attached_hexane_kg, inventory.internal_hexane_kg)
            )
        hexane_decrement = _phase_decrement_per_particle(
            total_to_cell=hexane_total_to_cell,
            inventory=inventory,
            weighted_phase_total=hexane_basis_total,
            weighted_dry_total=basis.dry_matter_kg,
            phase_value=hexane_basis_value,
        )
        arm = _layer_sorbed_arm(request, layer)
        own_water = water_class_of.get(packet.packet_id)
        if own_water is None:
            sorbed_arm_active = layer_total.external_water_kg == 0.0 and arm is not None
            water_total_to_cell = transfer.water_from_packets_kg
            water_basis = layer_total
            # D9-a-4 / D19: the layer's weighted front - on a layer not composed
            # per class, the layer's own report.
            packet_front = front_by_layer[layer]
        else:
            # Level 0b (ii): the packet draws on ITS class's own water demand,
            # on its class's basis, by the same proportional law; (i) the drain
            # and the floor below read the packet's OWN front.
            sorbed_arm_active = own_water.sorbed_arm_active
            water_total_to_cell = own_water.demand_kg
            water_basis = own_water.inventory
            packet_front = own_water.own_front_by_packet[packet.packet_id]
        if sorbed_arm_active:
            # C8 Tier 1: allocation basis DRY MATTER, so every packet's
            # retained loading moves by the same increment - which IS the
            # instant-sorption case of D-C8-3.  Both branches of the helper
            # then reduce to the same dry-matter proportional split, and the
            # drawable-stock guard is the layer-aggregate one above plus the
            # per-packet floor check below.
            water_basis_total = water_basis.dry_matter_kg
            water_basis_value = inventory.dry_matter_kg
        else:
            water_basis_total = water_basis.external_water_kg
            water_basis_value = inventory.external_water_kg
        water_decrement = _phase_decrement_per_particle(
            total_to_cell=water_total_to_cell,
            inventory=inventory,
            weighted_phase_total=water_basis_total,
            weighted_dry_total=water_basis.dry_matter_kg,
            phase_value=water_basis_value,
        )
        representative_share = contribution.representative_particles / math.fsum(
            item.representative_particles for item in layer_total.contributions
        )
        weighted_packet_energy_decrement = math.fsum(
            (
                steps[layer].accepted.macro_step_s
                * contribution.ua_share_w_k
                * (contribution.temperature_k - steps[layer].fast.state.interface_temperature_k),
                representative_share * remainder_by_layer[layer],
            )
        )
        energy_decrement = weighted_packet_energy_decrement / contribution.representative_particles
        total_hexane_after = math.fsum(
            (
                inventory.attached_hexane_kg,
                inventory.internal_hexane_kg,
                -hexane_decrement,
            )
        )
        retained_water_after = inventory.retained_water_kg
        external_water_after = inventory.external_water_kg
        if sorbed_arm_active:
            # D9-a-1 (owner ruling 2026-09-13): the routing is no longer a
            # switch.  CONDENSATE is booked to the external film ALWAYS -
            # liquid water at the surface sets the surface activity to one and
            # a closed shell cannot take it in (the closure record, section 4).
            # An evaporative DRAW comes out of the open shell, and only out of
            # it: behind the front the water is inert (D9-a-3).
            if water_decrement > 0.0:
                retained_water_after = inventory.retained_water_kg - water_decrement
                # D9-a-3 (brief item (c)): the floor is retained MINUS the
                # exchangeable stock, not the bare Luikov floor.
                floor_kg = arm.shell_exhaustion_floor_kg(
                    dry_matter_kg=inventory.dry_matter_kg,
                    retained_water_kg=inventory.retained_water_kg,
                    front_fraction=packet_front,
                )
                if retained_water_after < floor_kg:
                    if own_water is None:
                        raise SP1K2Law2PhaseExhaustionError(
                            "Law-2 flux would deplete an external packet phase; "
                            "step rejected, not clamped",
                            vertical_layer_id=layer,
                            phase_name="retained_water",
                            layer_stock_kg=sorbed_water.exchangeable_water_kg(
                                dry_matter_kg=layer_total.dry_matter_kg,
                                retained_water_kg=layer_total.retained_water_kg,
                                front_fraction=front_by_layer[layer],
                            ),
                            demanded_kg=transfer.water_from_packets_kg,
                            layer_macro_step_s=steps[layer].accepted.macro_step_s,
                            interval_end_time_s=getattr(request, "end_time_s", None),
                        )
                    raise SP1K2Law2PhaseExhaustionError(
                        "Law-2 flux would deplete an external packet phase; "
                        "step rejected, not clamped",
                        vertical_layer_id=layer,
                        phase_name="retained_water",
                        layer_stock_kg=own_water.stock_kg,
                        demanded_kg=own_water.demand_kg,
                        layer_macro_step_s=steps[layer].accepted.macro_step_s,
                        interval_end_time_s=getattr(request, "end_time_s", None),
                        water_class_key=own_water.class_key,
                    )
            else:
                external_water_after = math.fsum((inventory.external_water_kg, -water_decrement))
        else:
            external_water_after = inventory.external_water_kg - water_decrement
        if arm is not None:
            # D9-a-2 (owner ruling 2026-09-13): THE DRAIN.  After the Law-2
            # advance and before the film-onset check, the surface film
            # imbibes into the hexane-free shell - the film-onset deposit's
            # pattern run in the REVERSE direction, on every packet of the
            # layer in proportion to its dry matter (the capacity is
            # proportional to it, row D9-a-4), at the nominal or at zero as
            # the arm's token selects (row D9-a-1).  The moved water carries
            # its enthalpy at the layer temperature and no sorption heat: the
            # shipped net excess binding enthalpy is 0.0 (props/sorption.py
            # :525-532), so no caloric term accompanies the move.  Not a
            # clamp: the mass is moved, never discarded, and the split is the
            # law module's, which closes the pair total under fsum.
            # Level 0b (i): on a composed layer ``packet_front`` is the packet's
            # OWN front - exactly 1 under an attached hexane film, where the
            # capacity is exactly 0.0 and the water film stays; elsewhere it is
            # the layer's weighted front, as before.
            drain_kg = arm.layer_film_drain_kg(
                film_kg=external_water_after,
                dry_matter_kg=inventory.dry_matter_kg,
                retained_water_kg=retained_water_after,
                front_fraction=packet_front,
            )
            external_water_after, retained_water_after = sorbed_water.drain_split_kg(
                film_kg=external_water_after,
                retained_water_kg=retained_water_after,
                drain_kg=drain_kg,
            )
        if sorbed_arm_active:
            # C8 Tier 1 section 1.5: film onset as a DEPOSIT inside the same
            # advance call, kept verbatim as defense in depth.  Under D9-a it
            # is UNREACHABLE by construction: the retained water can only rise
            # through the drain, whose capacity is exactly the distance to
            # this ceiling (row D9-a-2).  D9-2 must assert that, not assume it.
            ceiling_kg = arm.film_onset_ceiling_kg(inventory.dry_matter_kg)
            if retained_water_after > ceiling_kg:
                surplus_kg = retained_water_after - ceiling_kg
                retained_water_after = ceiling_kg
                external_water_after = math.fsum((external_water_after, surplus_kg))
        if total_hexane_after < 0.0 or external_water_after < 0.0 or retained_water_after < 0.0:
            raise SP1K2Law2StepError(
                "packet allocation left the mass domain; step rejected, not clamped"
            )
        # C8 Tier 1: the retained-water keyword is passed ONLY when the sorbed
        # arm actually moved it.  Two reasons, both load-bearing: every film-arm
        # call reaches the codec with the exact argument list it had before
        # 2026-09-05 (byte-identity by construction, not by argument-value
        # coincidence), and the tag-2 / dual-variant codecs, which do not carry
        # the keyword and are outside this packet's write set, keep working
        # unchanged on the film arm they already serve.
        retained_water_keyword = (
            {}
            if retained_water_after == inventory.retained_water_kg
            else {"retained_water_mass_after_kg": retained_water_after}
        )
        complete_after = codec.advance(
            packet.entry.payload_bytes,
            total_hexane_mass_after_kg=total_hexane_after,
            external_free_water_mass_after_kg=external_water_after,
            envelope_energy_after_j=inventory.common_datum_energy_j - energy_decrement,
            **retained_water_keyword,
        )
        energy_resolution_rows[layer].append(
            contribution.representative_particles
            * math.fsum(
                (
                    _binary64_resolution(inventory.common_datum_energy_j),
                    _binary64_resolution(complete_after.inventory.common_datum_energy_j),
                )
            )
        )
        realized_rows[layer].append(
            (
                contribution.representative_particles
                * math.fsum(
                    (
                        inventory.attached_hexane_kg,
                        inventory.internal_hexane_kg,
                        -complete_after.inventory.attached_hexane_kg,
                        -complete_after.inventory.internal_hexane_kg,
                    )
                ),
                # C8 Tier 1: the realized WATER draw is the packet's TOTAL
                # water change - external plus sorbed - so the exchange ledger
                # closes on the sorbed arm and on the film-onset deposit alike.
                # On every film-arm packet the sorbed term is exactly zero and
                # ``math.fsum`` is correctly rounded, so this is bit-for-bit
                # the single subtraction it replaces.
                contribution.representative_particles
                * math.fsum(
                    (
                        inventory.external_water_kg,
                        inventory.retained_water_kg,
                        -complete_after.inventory.external_water_kg,
                        -complete_after.inventory.retained_water_kg,
                    )
                ),
                contribution.representative_particles
                * (
                    inventory.common_datum_energy_j - complete_after.inventory.common_datum_energy_j
                ),
            )
        )
        entries[packet.packet_id] = tray_host.PacketPlanEntry(
            packet_id=packet.packet_id,
            payload_bytes=complete_after.payload_bytes,
            inventory=complete_after.inventory,
            weight=packet.entry.weight,
            surface=complete_after.surface,
        )
    ledgers: list[LayerPacketExchangeLedger] = []
    for layer in (1, 2):
        layer_total = totals[layer]
        transfer = by_layer[layer]
        step = steps[layer]
        packet_conductive_power = math.fsum(
            contribution.ua_share_w_k
            * (contribution.temperature_k - step.fast.state.interface_temperature_k)
            for contribution in layer_total.contributions
        )
        cell_conductive_power = step.fast.evaluation.solid_to_interface_w
        realized_hexane = math.fsum(row[0] for row in realized_rows[layer])
        realized_water = math.fsum(row[1] for row in realized_rows[layer])
        realized_energy = math.fsum(row[2] for row in realized_rows[layer])
        energy_residual = realized_energy - transfer.energy_from_packets_j
        energy_resolution_bound = math.fsum(energy_resolution_rows[layer])
        energy_residual_beyond_resolution = _residual_beyond_resolution(
            energy_residual,
            energy_resolution_bound,
        )
        thermal_residual = packet_conductive_power - cell_conductive_power
        maximum = max(
            _relative_residual(
                layer_total.ua_partition_residual_w_k,
                layer_total.ua_share_sum_w_k,
                layer_total.solid_to_interface_ua_w_k,
            ),
            _relative_residual(
                layer_total.ua_temperature_moment_residual_w,
                layer_total.solid_to_interface_ua_w_k * layer_total.solid_temperature_k,
            ),
            _relative_residual(
                thermal_residual,
                packet_conductive_power,
                cell_conductive_power,
            ),
            _relative_residual(
                realized_hexane - transfer.hexane_from_packets_kg,
                realized_hexane,
                transfer.hexane_from_packets_kg,
            ),
            _relative_residual(
                realized_water - transfer.water_from_packets_kg,
                realized_water,
                transfer.water_from_packets_kg,
            ),
            _relative_residual(
                energy_residual_beyond_resolution,
                realized_energy,
                transfer.energy_from_packets_j,
            ),
        )
        ledgers.append(
            LayerPacketExchangeLedger(
                layer=layer,
                aggregation_law_id=layer_total.aggregation_law_id,
                energy_allocation_law_id=layer_total.energy_allocation_law_id,
                declared_solid_to_interface_ua_w_k=(layer_total.solid_to_interface_ua_w_k),
                ua_share_sum_w_k=layer_total.ua_share_sum_w_k,
                ua_partition_residual_w_k=layer_total.ua_partition_residual_w_k,
                packet_solid_to_interface_w=packet_conductive_power,
                cell_solid_to_interface_w=cell_conductive_power,
                thermal_moment_residual_w=thermal_residual,
                nonconductive_remainder_energy_j=remainder_by_layer[layer],
                requested_hexane_transfer_kg=transfer.hexane_from_packets_kg,
                realized_hexane_transfer_kg=realized_hexane,
                requested_water_transfer_kg=transfer.water_from_packets_kg,
                realized_water_transfer_kg=realized_water,
                requested_energy_transfer_j=transfer.energy_from_packets_j,
                realized_energy_transfer_j=realized_energy,
                energy_transfer_residual_j=energy_residual,
                energy_binary64_resolution_bound_j=energy_resolution_bound,
                energy_transfer_residual_beyond_resolution_j=(energy_residual_beyond_resolution),
                energy_resolution_law_id=(PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID),
                maximum_relative_residual=maximum,
                relative_limit=request.relative_limit,
            )
        )
    return entries, (ledgers[0], ledgers[1])


@dataclass(frozen=True, slots=True, kw_only=True)
class _AcceptedLaw2PacketAdvance:
    entries: dict[str, tray_host.PacketPlanEntry]
    residual_by_layer: tuple[float, float]
    maximum_scaled_residual: float
    callback_id: str = "sp1-k2-law2-accepted-packet-advance-v1"
    source_identity: str = "cell-engineering-feasibility-accepted-flux"

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def propose(
        self, request: tray_host.ManufacturedPacketAdvanceRequest
    ) -> tray_host.ManufacturedPacketAdvanceProposal:
        layer = request.packet.owner_key.vertical_layer_id.value
        return tray_host.ManufacturedPacketAdvanceProposal(
            packet_id=request.packet.packet_id,
            entry_after=self.entries[request.packet.packet_id],
            callback_id=self.callback_id,
            source_identity=self.source_identity,
            residual_contract_passed=True,
            maximum_scaled_residual=self.residual_by_layer[layer - 1],
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class WholeTrayExternalEnergyLedger:
    boundary_enthalpy_convention: str
    component_datum_adapter_digest: str
    component_datum_reference_temperature_k: float
    component_datum_reference_pressure_pa: float
    component_numeric_offsets: tuple[float, float, float, float]
    component_common_datum_species_shifts_j_mol: tuple[float, float]
    molar_heat_capacity_j_mol_k: float
    energy_datum_temperature_k: float
    energy_datum_id: str
    hexane_latent_heat_j_mol: float
    water_latent_heat_j_mol: float
    internal_face_outlet_energy_w: float
    internal_face_inlet_energy_w: float
    internal_face_cancellation_w: float
    bottom_inlet_energy_w: float
    top_outlet_energy_w: float
    jacket_and_ambient_energy_w: float
    per_cell_external_energy_j: float
    external_only_tray_energy_j: float
    internal_face_energy_residual_j: float
    whole_tray_system_energy_change_j: float
    whole_tray_external_energy_residual_j: float
    whole_tray_energy_binary64_resolution_bound_j: float
    whole_tray_external_energy_residual_beyond_resolution_j: float
    energy_resolution_law_id: str
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False
    numeric_component_datum_bound: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        return (
            self.internal_face_cancellation_w == 0.0
            and self.energy_resolution_law_id == WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID
            and self.whole_tray_external_energy_residual_beyond_resolution_j
            == _residual_beyond_resolution(
                self.whole_tray_external_energy_residual_j,
                self.whole_tray_energy_binary64_resolution_bound_j,
            )
            and self.maximum_relative_residual <= self.relative_limit
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class NativeWholeTrayExternalEnergyLedger:
    """CELL-02d W5/W6: the native-mode whole-tray external energy ledger.

    The internal gas face between the K cells is not enthalpy-free in the
    native mode: each cell books its state at its own layer pressure, so the
    verbatim hand-off stream carries the declared throttling defect
    ``flow * h(T, p_upper, y) - flow * h(T, p_lower, y)`` - the unmodeled
    Joule-Thomson relaxation of the well-mixed-cell form. The defect is
    booked as its own field, every closure residual is corrected by exactly
    ``dt * defect``, and the pass gate reconstructs the defect bit-exactly
    from the stored face facts - never against a tolerance. The declared
    ledger and its exact-zero cancellation gate are untouched.
    """

    boundary_enthalpy_convention: str
    energy_datum_id: str
    face_temperature_k: float
    face_hexane_mole_fraction: float
    face_molar_flow_mol_s: float
    face_lower_layer_pressure_pa: float
    face_upper_layer_pressure_pa: float
    face_lower_pressure_enthalpy_j_mol: float
    face_upper_pressure_enthalpy_j_mol: float
    internal_face_outlet_energy_w: float
    internal_face_inlet_energy_w: float
    internal_face_cancellation_w: float
    declared_face_throttling_defect_w: float
    bottom_inlet_energy_w: float
    top_outlet_energy_w: float
    jacket_and_ambient_energy_w: float
    per_cell_external_energy_j: float
    external_only_tray_energy_j: float
    defect_corrected_internal_face_energy_residual_j: float
    whole_tray_system_energy_change_j: float
    defect_corrected_whole_tray_external_energy_residual_j: float
    whole_tray_energy_binary64_resolution_bound_j: float
    whole_tray_external_energy_residual_beyond_resolution_j: float
    energy_resolution_law_id: str
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False
    numeric_component_datum_bound: ClassVar[bool] = False
    face_throttling_defect_is_booked_not_smoothed: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        reconstructed_defect = (
            self.face_molar_flow_mol_s * self.face_upper_pressure_enthalpy_j_mol
            - self.face_molar_flow_mol_s * self.face_lower_pressure_enthalpy_j_mol
        )
        return (
            self.boundary_enthalpy_convention == NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION
            and self.declared_face_throttling_defect_w == reconstructed_defect
            and self.internal_face_cancellation_w == self.declared_face_throttling_defect_w
            and self.energy_resolution_law_id == WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID
            and self.whole_tray_external_energy_residual_beyond_resolution_j
            == _residual_beyond_resolution(
                self.defect_corrected_whole_tray_external_energy_residual_j,
                self.whole_tray_energy_binary64_resolution_bound_j,
            )
            and self.maximum_relative_residual <= self.relative_limit
        )


def _weighted_packet_common_energy_resolution(packet: tray_host.AcceptedTrayPacket) -> float:
    weight = packet.entry.weight.representative_particles
    per_particle = packet.entry.inventory.common_datum_energy_j
    weighted = packet.weighted_totals.totals[6]
    return math.fsum(
        (
            abs(per_particle) * _binary64_resolution(weight),
            abs(weight) * _binary64_resolution(per_particle),
            _binary64_resolution(weighted),
        )
    )


def _whole_tray_external_energy_ledger(
    request: SP1K2Law2StepRequest,
    inputs: tuple[ef.EngineeringCellInputs, ef.EngineeringCellInputs],
    steps: tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep],
    transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers],
    host_step: tray_host.ValidatedKCellTrayStep,
    *,
    request_codec: high_loading.HighLoadingV1PacketAdapter,
) -> WholeTrayExternalEnergyLedger | NativeWholeTrayExternalEnergyLedger:
    properties = inputs[0].properties
    dt = inputs[0].macro_step_s
    lower_evaluation = steps[0].fast.evaluation
    lower_state = steps[0].fast.state
    upper_boundary = inputs[1].gas_boundary
    upper_evaluation = steps[1].fast.evaluation
    upper_state = steps[1].fast.state
    internal_outlet = lower_evaluation.outlet_molar_flow_mol_s * (
        _stream_enthalpy_j_mol(
            properties,
            temperature_k=lower_state.gas_temperature_k,
            pressure_pa=lower_state.layer_pressure_pa,
            hexane_mole_fraction=lower_state.gas_hexane_mole_fraction,
            water_mole_fraction=lower_state.gas_water_mole_fraction,
            datum_adapter=request_codec.component_datum_adapter,
        )
    )
    internal_inlet = upper_boundary.inlet_molar_flow_mol_s * (
        _stream_enthalpy_j_mol(
            properties,
            temperature_k=upper_boundary.inlet_temperature_k,
            pressure_pa=upper_state.layer_pressure_pa,
            hexane_mole_fraction=upper_boundary.inlet_hexane_mole_fraction,
            water_mole_fraction=upper_boundary.inlet_water_mole_fraction,
            datum_adapter=request_codec.component_datum_adapter,
        )
    )
    internal_cancellation = internal_inlet - internal_outlet
    bottom = inputs[0].gas_boundary
    bottom_inlet = bottom.inlet_molar_flow_mol_s * (
        _stream_enthalpy_j_mol(
            properties,
            temperature_k=bottom.inlet_temperature_k,
            pressure_pa=lower_state.layer_pressure_pa,
            hexane_mole_fraction=bottom.inlet_hexane_mole_fraction,
            water_mole_fraction=bottom.inlet_water_mole_fraction,
            datum_adapter=request_codec.component_datum_adapter,
        )
    )
    top_outlet = upper_evaluation.outlet_molar_flow_mol_s * (
        _stream_enthalpy_j_mol(
            properties,
            temperature_k=upper_state.gas_temperature_k,
            pressure_pa=upper_state.layer_pressure_pa,
            hexane_mole_fraction=upper_state.gas_hexane_mole_fraction,
            water_mole_fraction=upper_state.gas_water_mole_fraction,
            datum_adapter=request_codec.component_datum_adapter,
        )
    )
    jacket_ambient = math.fsum(step.wall.steam_side_w - step.wall.to_ambient_w for step in steps)
    per_cell_external = math.fsum(
        transfer.gas_and_wall_boundary_energy_to_cell_j for transfer in transfers
    )
    external_power = math.fsum((bottom_inlet, -top_outlet, jacket_ambient))
    external_only = dt * external_power
    internal_residual = per_cell_external - external_only
    packet_before = math.fsum(
        packet.weighted_totals.totals[6] for packet in request.prior.tray.packets
    )
    packet_after = math.fsum(
        packet.weighted_totals.totals[6] for packet in host_step.candidate.packets
    )
    cell_before = math.fsum(
        cell.conserved_state.totals.common_datum_energy_j for cell in request.prior.tray.cells
    )
    cell_after = math.fsum(
        cell.conserved_state.totals.common_datum_energy_j for cell in host_step.candidate.cells
    )
    system_change = math.fsum((packet_after, cell_after, -packet_before, -cell_before))
    system_residual = system_change - external_only
    state_resolution = math.fsum(
        (
            *(
                _weighted_packet_common_energy_resolution(packet)
                for packet in request.prior.tray.packets
            ),
            *(
                _weighted_packet_common_energy_resolution(packet)
                for packet in host_step.candidate.packets
            ),
            *(
                _binary64_resolution(cell.conserved_state.totals.common_datum_energy_j)
                for cell in request.prior.tray.cells
            ),
            *(
                _binary64_resolution(cell.conserved_state.totals.common_datum_energy_j)
                for cell in host_step.candidate.cells
            ),
        )
    )
    aggregate_operation_resolution = math.fsum(
        _binary64_resolution(value)
        for value in (
            packet_before,
            packet_after,
            cell_before,
            cell_after,
            system_change,
        )
    )
    external_operation_resolution = math.fsum(
        (
            abs(dt)
            * math.fsum(
                _binary64_resolution(value) for value in (bottom_inlet, top_outlet, jacket_ambient)
            ),
            _product_binary64_resolution(dt, external_power, external_only),
        )
    )
    system_resolution_bound = math.fsum(
        (state_resolution, aggregate_operation_resolution, external_operation_resolution)
    )
    system_residual_beyond_resolution = _residual_beyond_resolution(
        system_residual,
        system_resolution_bound,
    )
    if type(properties) is not cc.DeclaredGasProperties:
        # CELL-02d W5 native arm: the verbatim face hand-off carries the
        # booked throttling defect; every closure residual is corrected by
        # exactly dt * defect and the correction is reconstructible.
        face_flow = upper_boundary.inlet_molar_flow_mol_s
        face_temperature = upper_boundary.inlet_temperature_k
        face_hexane = upper_boundary.inlet_hexane_mole_fraction
        face_lower_enthalpy = native_caloric.mixture_enthalpy(
            face_temperature, lower_state.layer_pressure_pa, face_hexane
        ).value_j_mol
        face_upper_enthalpy = native_caloric.mixture_enthalpy(
            face_temperature, upper_state.layer_pressure_pa, face_hexane
        ).value_j_mol
        face_defect = face_flow * face_upper_enthalpy - face_flow * face_lower_enthalpy
        if internal_cancellation != face_defect:
            raise SP1K2Law2StepError(
                "the native internal gas face is not a verbatim lower-outlet "
                "hand-off; the booked throttling defect cannot be reconstructed"
            )
        defect_energy = dt * face_defect
        native_internal_residual = per_cell_external - external_only - defect_energy
        native_system_residual = system_change - external_only - defect_energy
        native_resolution_bound = math.fsum(
            (
                system_resolution_bound,
                _product_binary64_resolution(dt, face_defect, defect_energy),
            )
        )
        native_beyond = _residual_beyond_resolution(
            native_system_residual,
            native_resolution_bound,
        )
        native_maximum = max(
            host_step.tray_ledger.maximum_relative_residual,
            _relative_residual(native_internal_residual, per_cell_external, external_only),
            _relative_residual(native_beyond, system_change, external_only),
        )
        return NativeWholeTrayExternalEnergyLedger(
            boundary_enthalpy_convention=NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION,
            energy_datum_id=properties.energy_datum_id,
            face_temperature_k=face_temperature,
            face_hexane_mole_fraction=face_hexane,
            face_molar_flow_mol_s=face_flow,
            face_lower_layer_pressure_pa=lower_state.layer_pressure_pa,
            face_upper_layer_pressure_pa=upper_state.layer_pressure_pa,
            face_lower_pressure_enthalpy_j_mol=face_lower_enthalpy,
            face_upper_pressure_enthalpy_j_mol=face_upper_enthalpy,
            internal_face_outlet_energy_w=internal_outlet,
            internal_face_inlet_energy_w=internal_inlet,
            internal_face_cancellation_w=internal_cancellation,
            declared_face_throttling_defect_w=face_defect,
            bottom_inlet_energy_w=bottom_inlet,
            top_outlet_energy_w=top_outlet,
            jacket_and_ambient_energy_w=jacket_ambient,
            per_cell_external_energy_j=per_cell_external,
            external_only_tray_energy_j=external_only,
            defect_corrected_internal_face_energy_residual_j=native_internal_residual,
            whole_tray_system_energy_change_j=system_change,
            defect_corrected_whole_tray_external_energy_residual_j=native_system_residual,
            whole_tray_energy_binary64_resolution_bound_j=native_resolution_bound,
            whole_tray_external_energy_residual_beyond_resolution_j=native_beyond,
            energy_resolution_law_id=WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID,
            maximum_relative_residual=native_maximum,
            relative_limit=request.relative_limit,
        )
    cp = properties.molar_heat_capacity_j_mol_k
    datum = properties.energy_datum_temperature_k
    maximum = max(
        host_step.tray_ledger.maximum_relative_residual,
        _relative_residual(internal_cancellation, internal_inlet, internal_outlet),
        _relative_residual(internal_residual, per_cell_external, external_only),
        _relative_residual(
            system_residual_beyond_resolution,
            system_change,
            external_only,
        ),
    )
    return WholeTrayExternalEnergyLedger(
        boundary_enthalpy_convention=GAS_BOUNDARY_ENTHALPY_CONVENTION,
        component_datum_adapter_digest=(request_codec.component_datum_adapter.definition_digest),
        component_datum_reference_temperature_k=(
            request_codec.component_datum_adapter.reference_temperature_k
        ),
        component_datum_reference_pressure_pa=(
            request_codec.component_datum_adapter.reference_pressure_pa
        ),
        component_numeric_offsets=request_codec.component_datum_adapter.numeric_offsets,
        component_common_datum_species_shifts_j_mol=(
            request_codec.component_datum_adapter.common_datum_species_shifts_j_mol
        ),
        molar_heat_capacity_j_mol_k=cp,
        energy_datum_temperature_k=datum,
        energy_datum_id=properties.energy_datum_id,
        hexane_latent_heat_j_mol=properties.hexane_latent_heat_j_mol,
        water_latent_heat_j_mol=properties.water_latent_heat_j_mol,
        internal_face_outlet_energy_w=internal_outlet,
        internal_face_inlet_energy_w=internal_inlet,
        internal_face_cancellation_w=internal_cancellation,
        bottom_inlet_energy_w=bottom_inlet,
        top_outlet_energy_w=top_outlet,
        jacket_and_ambient_energy_w=jacket_ambient,
        per_cell_external_energy_j=per_cell_external,
        external_only_tray_energy_j=external_only,
        internal_face_energy_residual_j=internal_residual,
        whole_tray_system_energy_change_j=system_change,
        whole_tray_external_energy_residual_j=system_residual,
        whole_tray_energy_binary64_resolution_bound_j=system_resolution_bound,
        whole_tray_external_energy_residual_beyond_resolution_j=(system_residual_beyond_resolution),
        energy_resolution_law_id=WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID,
        maximum_relative_residual=maximum,
        relative_limit=request.relative_limit,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class CoupledTrayLedger:
    layer_hexane_mass_residuals_kg: tuple[float, float]
    layer_water_mass_residuals_kg: tuple[float, float]
    layer_wall_energy_residuals_j: tuple[float, float]
    layer_wall_energy_binary64_resolution_bounds_j: tuple[float, float]
    layer_wall_energy_residuals_beyond_resolution_j: tuple[float, float]
    layer_wall_relative_residuals: tuple[float, float]
    layer_wall_row_residuals_w: tuple[float, float]
    layer_wall_row_binary64_resolution_bounds_w: tuple[float, float]
    layer_wall_row_residuals_beyond_resolution_w: tuple[float, float]
    layer_wall_row_relative_residuals: tuple[float, float]
    energy_resolution_law_id: str
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        expected_energy = tuple(
            _residual_beyond_resolution(residual, bound)
            for residual, bound in zip(
                self.layer_wall_energy_residuals_j,
                self.layer_wall_energy_binary64_resolution_bounds_j,
                strict=True,
            )
        )
        expected_rows = tuple(
            _residual_beyond_resolution(residual, bound)
            for residual, bound in zip(
                self.layer_wall_row_residuals_w,
                self.layer_wall_row_binary64_resolution_bounds_w,
                strict=True,
            )
        )
        return (
            self.energy_resolution_law_id == COUPLED_WALL_ENERGY_BINARY64_RESOLUTION_LAW_ID
            and self.layer_wall_energy_residuals_beyond_resolution_j == expected_energy
            and self.layer_wall_row_residuals_beyond_resolution_w == expected_rows
            and self.maximum_relative_residual <= self.relative_limit
        )


def _wall_row_binary64_resolution_bound(macro: ef.EngineeringCellMacroStep) -> float:
    accepted = macro.accepted
    wall = accepted.wall_parameters
    dt = accepted.macro_step_s
    scale = accepted.wall_storage_scale
    capacity = wall.thermal_capacity_j_k
    storage_conductance = scale * capacity / dt
    conductance_resolution = math.fsum(
        (
            abs(capacity / dt) * _binary64_resolution(scale),
            abs(scale / dt) * _binary64_resolution(capacity),
            abs(storage_conductance / dt) * _binary64_resolution(dt),
            _binary64_resolution(storage_conductance),
        )
    )
    temperature_delta = macro.wall.wall_temperature_k - macro.wall.old_wall_temperature_k
    temperature_delta_resolution = math.fsum(
        (
            _binary64_resolution(macro.wall.wall_temperature_k),
            _binary64_resolution(macro.wall.old_wall_temperature_k),
            _binary64_resolution(temperature_delta),
        )
    )
    storage_resolution = math.fsum(
        (
            abs(storage_conductance) * temperature_delta_resolution,
            abs(temperature_delta) * conductance_resolution,
            _binary64_resolution(macro.wall.storage_w),
        )
    )
    steam_delta = wall.steam_temperature_k - macro.wall.wall_temperature_k
    steam_resolution = math.fsum(
        (
            abs(steam_delta) * _binary64_resolution(wall.steam_side_ua_w_k),
            abs(wall.steam_side_ua_w_k)
            * math.fsum(
                (
                    _binary64_resolution(wall.steam_temperature_k),
                    _binary64_resolution(macro.wall.wall_temperature_k),
                )
            ),
            _binary64_resolution(macro.wall.steam_side_w),
        )
    )
    ambient_delta = macro.wall.wall_temperature_k - wall.ambient_temperature_k
    ambient_resolution = math.fsum(
        (
            abs(ambient_delta) * _binary64_resolution(wall.ambient_ua_w_k),
            abs(wall.ambient_ua_w_k)
            * math.fsum(
                (
                    _binary64_resolution(macro.wall.wall_temperature_k),
                    _binary64_resolution(wall.ambient_temperature_k),
                )
            ),
            _binary64_resolution(macro.wall.to_ambient_w),
        )
    )
    row_operation_resolution = math.fsum(
        _binary64_resolution(value)
        for value in (
            macro.wall.storage_w,
            macro.wall.steam_side_w,
            macro.wall.to_ambient_w,
            macro.wall.shared_to_gas_w,
            macro.wall.shared_to_interface_w,
            macro.wall.residual_w,
        )
    )
    return math.fsum(
        (
            storage_resolution,
            steam_resolution,
            ambient_resolution,
            row_operation_resolution,
        )
    )


def _layer_wall_energy_binary64_resolution_bound(
    *,
    before: tray_host.AcceptedKCell,
    after: tray_host.AcceptedKCell,
    macro: ef.EngineeringCellMacroStep,
    transfer: LayerAcceptedTransfers,
    packet_ledger: LayerPacketExchangeLedger,
) -> float:
    before_energy = before.conserved_state.totals.common_datum_energy_j
    after_energy = after.conserved_state.totals.common_datum_energy_j
    delta_energy = after_energy - before_energy
    wall_storage = macro.wall.storage_w * macro.accepted.macro_step_s
    return math.fsum(
        (
            packet_ledger.energy_binary64_resolution_bound_j,
            _binary64_resolution(before_energy),
            _binary64_resolution(after_energy),
            _binary64_resolution(delta_energy),
            _product_binary64_resolution(
                macro.wall.storage_w,
                macro.accepted.macro_step_s,
                wall_storage,
            ),
            _binary64_resolution(transfer.energy_from_packets_j),
            _binary64_resolution(transfer.gas_and_wall_boundary_energy_to_cell_j),
        )
    )


def _coupled_ledger(
    request: SP1K2Law2StepRequest,
    host_step: tray_host.ValidatedKCellTrayStep,
    macro_steps: tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep],
    transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers],
    whole_tray_energy: WholeTrayExternalEnergyLedger | NativeWholeTrayExternalEnergyLedger,
    tray_pressure: AggregateSP1PressureLedger,
    packet_exchange_ledgers: tuple[LayerPacketExchangeLedger, LayerPacketExchangeLedger],
) -> CoupledTrayLedger:
    hexane_residuals: list[float] = []
    water_residuals: list[float] = []
    energy_residuals: list[float] = []
    energy_resolution_bounds: list[float] = []
    energy_residuals_beyond_resolution: list[float] = []
    wall_row_residuals: list[float] = []
    wall_row_resolution_bounds: list[float] = []
    wall_row_residuals_beyond_resolution: list[float] = []
    wall_row_relative_residuals: list[float] = []
    relative: list[float] = [
        host_step.tray_ledger.maximum_relative_residual,
        whole_tray_energy.maximum_relative_residual,
        tray_pressure.maximum_relative_residual,
        *(ledger.maximum_relative_residual for ledger in packet_exchange_ledgers),
    ]
    for index, (before, after, macro, transfer) in enumerate(
        zip(
            request.prior.tray.cells,
            host_step.candidate.cells,
            macro_steps,
            transfers,
            strict=True,
        )
    ):
        delta_hexane = (
            after.conserved_state.totals.hexane_kg - before.conserved_state.totals.hexane_kg
        )
        delta_water = after.conserved_state.totals.water_kg - before.conserved_state.totals.water_kg
        delta_energy = (
            after.conserved_state.totals.common_datum_energy_j
            - before.conserved_state.totals.common_datum_energy_j
        )
        wall_storage = macro.wall.storage_w * macro.accepted.macro_step_s
        energy_residual = delta_energy - wall_storage
        energy_resolution_bound = _layer_wall_energy_binary64_resolution_bound(
            before=before,
            after=after,
            macro=macro,
            transfer=transfer,
            packet_ledger=packet_exchange_ledgers[index],
        )
        energy_residual_beyond_resolution = _residual_beyond_resolution(
            energy_residual,
            energy_resolution_bound,
        )
        wall_row_residual = macro.ledger.wall_row_residual_w
        wall_row_resolution_bound = _wall_row_binary64_resolution_bound(macro)
        wall_row_residual_beyond_resolution = _residual_beyond_resolution(
            wall_row_residual,
            wall_row_resolution_bound,
        )
        wall_row_relative = _relative_residual(
            wall_row_residual_beyond_resolution,
            macro.wall.storage_w,
            macro.wall.steam_side_w,
            macro.wall.to_ambient_w,
            macro.wall.shared_to_gas_w,
            macro.wall.shared_to_interface_w,
        )
        hexane_residuals.append(delta_hexane)
        water_residuals.append(delta_water)
        energy_residuals.append(energy_residual)
        energy_resolution_bounds.append(energy_resolution_bound)
        energy_residuals_beyond_resolution.append(energy_residual_beyond_resolution)
        wall_row_residuals.append(wall_row_residual)
        wall_row_resolution_bounds.append(wall_row_resolution_bound)
        wall_row_residuals_beyond_resolution.append(wall_row_residual_beyond_resolution)
        wall_row_relative_residuals.append(wall_row_relative)
        relative.extend(
            (
                _relative_residual(
                    delta_hexane,
                    transfer.hexane_from_packets_kg,
                    transfer.gas_boundary_hexane_to_cell_kg,
                ),
                _relative_residual(
                    delta_water,
                    transfer.water_from_packets_kg,
                    transfer.gas_boundary_water_to_cell_kg,
                ),
                _relative_residual(
                    energy_residual_beyond_resolution,
                    delta_energy,
                    wall_storage,
                ),
                wall_row_relative,
                macro.fast.evaluation.scaled_residual_norm,
            )
        )
        if index not in (0, 1):  # pragma: no cover - exact K=2 zip guard
            raise AssertionError("unexpected SP1 K-layer index")
    wall_relative = tuple(
        _relative_residual(
            residual,
            raw_residual + macro.wall.storage_w * macro.accepted.macro_step_s,
        )
        for residual, raw_residual, macro in zip(
            energy_residuals_beyond_resolution,
            energy_residuals,
            macro_steps,
            strict=True,
        )
    )
    relative.extend(wall_relative)
    return CoupledTrayLedger(
        layer_hexane_mass_residuals_kg=tuple(hexane_residuals),
        layer_water_mass_residuals_kg=tuple(water_residuals),
        layer_wall_energy_residuals_j=tuple(energy_residuals),
        layer_wall_energy_binary64_resolution_bounds_j=tuple(energy_resolution_bounds),
        layer_wall_energy_residuals_beyond_resolution_j=tuple(energy_residuals_beyond_resolution),
        layer_wall_relative_residuals=wall_relative,
        layer_wall_row_residuals_w=tuple(wall_row_residuals),
        layer_wall_row_binary64_resolution_bounds_w=tuple(wall_row_resolution_bounds),
        layer_wall_row_residuals_beyond_resolution_w=tuple(wall_row_residuals_beyond_resolution),
        layer_wall_row_relative_residuals=tuple(wall_row_relative_residuals),
        energy_resolution_law_id=COUPLED_WALL_ENERGY_BINARY64_RESOLUTION_LAW_ID,
        maximum_relative_residual=max(relative),
        relative_limit=request.relative_limit,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SP1K2Law2StepRejection:
    attempt_id: str
    reason: str
    rollback_state: AcceptedSP1K2Law2TrayState

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


_VALIDATED_COUPLED_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedSP1K2Law2TrayStep:
    request: SP1K2Law2StepRequest
    candidate: AcceptedSP1K2Law2TrayState
    host_transaction: tray_host.ValidatedKCellTrayStep
    layer_inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory]
    cell_inputs: tuple[ef.EngineeringCellInputs, ef.EngineeringCellInputs]
    cell_macro_steps: tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep]
    layer_transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers]
    inter_k_face: InterKGasFaceLedger
    tray_pressure: AggregateSP1PressureLedger
    whole_tray_energy: WholeTrayExternalEnergyLedger | NativeWholeTrayExternalEnergyLedger
    packet_exchange_ledgers: tuple[LayerPacketExchangeLedger, LayerPacketExchangeLedger]
    coupled_ledger: CoupledTrayLedger
    codec_registry_digest: str
    auditor_identity_digest: str
    auditor_source_identity: str
    _seal: object = field(repr=False, compare=False)

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_COUPLED_SEAL:
            raise TypeError("validated coupled steps are issued only by the evaluator")


def evaluate_sp1_k2_law2_tray_step(
    request: SP1K2Law2StepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
) -> ValidatedSP1K2Law2TrayStep | SP1K2Law2StepRejection:
    """Evaluate one coupled SP1 macro-step through the shared K2 kernel."""

    if type(request) is not SP1K2Law2StepRequest:
        raise TypeError("SP1 integration requires an exact step request")

    # Local import avoids a module-initialization cycle: the shared kernel uses
    # the legacy ledger and private numerical functions, while this wrapper
    # retains the public SP1 types and digest domain byte for byte.
    from . import through_bed_k2_law2_kernel as through_bed_kernel

    contract = through_bed_kernel.reference_through_bed_k2_contract("SP1")
    normalized_prior = through_bed_kernel.KernelAcceptedThroughBedK2State(
        contract=contract,
        tray=request.prior.tray,
        layer_model_configuration_digest=request.prior.layer_model_configuration_digest,
        cell_field_states=request.prior.cell_field_states,
        wall_temperatures_k=request.prior.wall_temperatures_k,
    )
    normalized_request = through_bed_kernel.ThroughBedK2KernelRequest(
        attempt_id=request.attempt_id,
        prior=normalized_prior,
        end_time_s=request.end_time_s,
        layer_models=request.layer_models,
        layer_model_configuration_digest=canonical_sp1_k2_layer_model_digest(request.layer_models),
        bottom_gas_boundary=request.bottom_gas_boundary,
        bottom_gas_binding=through_bed_kernel.ThroughBedGasBoundaryBinding(
            physical_tray_id="SP1",
            boundary_kind=through_bed_kernel.ThroughBedGasBoundaryKind.DIRECT_STEAM_SOURCE,
            source_physical_tray_id=None,
            authority_id="legacy-sp1-direct-steam-boundary-v1",
        ),
        top_boundary_pressure_pa=request.top_boundary_pressure_pa,
        flow_segments=request.flow_segments,
        newton_tolerance=request.newton_tolerance,
        face_tolerance=request.face_tolerance,
        face_max_iterations=request.face_max_iterations,
        relative_limit=request.relative_limit,
    )
    result = through_bed_kernel.evaluate_through_bed_k2_kernel(
        normalized_request,
        codec=codec,
        k_planner=k_planner,
    )
    if type(result) is through_bed_kernel.ThroughBedK2KernelRejection:
        return SP1K2Law2StepRejection(
            attempt_id=request.attempt_id,
            reason=result.reason,
            rollback_state=request.prior,
        )
    if type(result) is not through_bed_kernel.ValidatedThroughBedK2KernelStep:
        return SP1K2Law2StepRejection(
            attempt_id=request.attempt_id,
            reason="shared through-bed kernel returned a foreign result",
            rollback_state=request.prior,
        )
    candidate = AcceptedSP1K2Law2TrayState(
        tray=result.candidate.tray,
        layer_model_configuration_digest=request.prior.layer_model_configuration_digest,
        cell_field_states=result.candidate.cell_field_states,
        wall_temperatures_k=result.candidate.wall_temperatures_k,
    )
    return ValidatedSP1K2Law2TrayStep(
        request=request,
        candidate=candidate,
        host_transaction=result.host_transaction,
        layer_inventories=result.layer_inventories,
        cell_inputs=result.cell_inputs,
        cell_macro_steps=result.cell_macro_steps,
        layer_transfers=result.layer_transfers,
        inter_k_face=result.inter_k_face,
        tray_pressure=result.tray_pressure,
        whole_tray_energy=result.whole_tray_energy,
        packet_exchange_ledgers=result.packet_exchange_ledgers,
        coupled_ledger=result.coupled_ledger,
        codec_registry_digest=result.codec_registry_digest,
        auditor_identity_digest=result.auditor_identity_digest,
        auditor_source_identity=result.auditor_source_identity,
        _seal=_VALIDATED_COUPLED_SEAL,
    )


def commit_sp1_k2_law2_tray_step(
    current: AcceptedSP1K2Law2TrayState,
    transaction: ValidatedSP1K2Law2TrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner | None = None,
) -> AcceptedSP1K2Law2TrayState:
    """Return one host-plus-walls candidate after all reversible checks.

    This is an in-memory all-or-nothing acceptance seam.  The caller still owns
    its lock and any durable transaction; persistent atomicity remains false.
    """

    if type(current) is not AcceptedSP1K2Law2TrayState:
        raise TypeError("coupled commit requires an exact accepted state")
    if (
        type(transaction) is not ValidatedSP1K2Law2TrayStep
        or transaction._seal is not _VALIDATED_COUPLED_SEAL
    ):
        raise SP1K2Law2IntegrationError("coupled commit requires a validated transaction")
    if current is not transaction.request.prior:
        raise StaleSP1K2Law2StateError("coupled transaction is stale or foreign")
    if current.state_digest != transaction.request.prior.state_digest:
        raise StaleSP1K2Law2StateError("coupled prior digest changed before commit")
    if (
        codec.codec_registry_digest != transaction.codec_registry_digest
        or codec.auditor_identity_digest != transaction.auditor_identity_digest
        or codec.source_identity != transaction.auditor_source_identity
    ):
        raise SP1K2Law2IntegrationError("coupled commit codec violates transaction pins")
    if k_planner is None:
        k_planner = tray_host.NoKRekeyPlanner()
    reevaluated = evaluate_sp1_k2_law2_tray_step(
        transaction.request,
        codec=codec,
        k_planner=k_planner,
    )
    if type(reevaluated) is not ValidatedSP1K2Law2TrayStep or reevaluated != transaction:
        raise SP1K2Law2IntegrationError(
            "coupled transaction differs from a complete deterministic commit-time reevaluation"
        )
    committed_tray = tray_host.commit_k_cell_tray_step(
        current.tray,
        transaction.host_transaction,
        auditor=codec,
    )
    if committed_tray is not transaction.candidate.tray:
        raise SP1K2Law2IntegrationError("host and coupled candidates diverged at commit")
    if (
        not transaction.inter_k_face.passed
        or not transaction.tray_pressure.passed
        or not transaction.whole_tray_energy.passed
        or not all(ledger.passed for ledger in transaction.packet_exchange_ledgers)
        or not transaction.coupled_ledger.passed
    ):
        raise SP1K2Law2IntegrationError("coupled acceptance ledgers no longer pass")
    return transaction.candidate


__all__ = (
    "AcceptedSP1K2Law2TrayState",
    "AggregateSP1PressureLedger",
    "CoupledTrayLedger",
    "COUPLED_WALL_ENERGY_BINARY64_RESOLUTION_LAW_ID",
    "DEFAULT_FACE_MAX_ITERATIONS",
    "DEFAULT_FACE_TOLERANCE",
    "DEFAULT_RELATIVE_LIMIT",
    "ENGINEERING_PACKET_AUDITOR_IDENTITY_DIGEST",
    "ENGINEERING_PACKET_CODEC_REGISTRY_DIGEST",
    "GAS_BOUNDARY_ENTHALPY_CONVENTION",
    "EngineeringFeasibilityPacketCodec",
    "EngineeringPacketCodecAuditor",
    "INTEGRATION_DIGEST_DOMAIN",
    "InterKGasFaceLedger",
    "LayerAcceptedTransfers",
    "ClassHexaneTransfer",
    "HEXANE_CLASS_ORDER",
    "LayerHexaneClass",
    "OPTION_A_PER_CLASS_HEXANE_LAW_ID",
    "PER_CLASS_HEXANE_DECLARED_ASSUMPTIONS",
    "ClassWaterTransfer",
    "LEVEL0B_CLASS_OWN_WATER_ROUTING_LAW_ID",
    "LayerWaterClass",
    "PER_CLASS_WATER_DECLARED_ASSUMPTIONS",
    "class_water_demand_kg",
    "layer_water_classes",
    "FREE_FILM_PIN_DECLARED_ASSUMPTIONS",
    "FREE_FILM_PIN_WATER_ACTIVITY",
    "FREE_FILM_PIN_WATER_CONDUCTANCE_FACTOR",
    "LEVEL0B_FREE_FILM_PIN_LAW_ID",
    "class_hexane_closures",
    "class_hexane_demand_kg",
    "hexane_depletion_stocks_by_class",
    "layer_hexane_classes",
    "LayerPacketExchangeLedger",
    "PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID",
    "PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID",
    "D19_DRY_MATTER_WEIGHTED_FRONT_LAW_ID",
    "D9A_FRONT_DECLARED_ASSUMPTIONS",
    "DECLARED_FALLING_RATE_MASS_BIOT",
    "DECLARED_LANE_DROP_BAND_LAW_ID",
    "PacketDerivedLayerInventory",
    "PacketLayerContribution",
    "SP1K2LaneDropBandError",
    "REGIME_C_DECLARED_ASSUMPTIONS",
    "REGIME_C_HEXANE_BINDING_LAW_ID",
    "REGIME_C_LAYER_PRESSURE_PA",
    "RegimeCHexaneBinding",
    "SP1K2SorptionFloorError",
    "regime_c_hexane_binding",
    "SP1K2Law2ConfigurationError",
    "SP1K2Law2IntegrationError",
    "SP1K2Law2StepError",
    "SP1K2Law2StepRejection",
    "SP1K2Law2StepRequest",
    "SP1KLayerEngineeringModel",
    "StaleSP1K2Law2StateError",
    "ValidatedSP1K2Law2TrayStep",
    "UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID",
    "WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID",
    "NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION",
    "NativeWholeTrayExternalEnergyLedger",
    "WholeTrayExternalEnergyLedger",
    "canonical_sp1_k2_layer_model_digest",
    "commit_sp1_k2_law2_tray_step",
    "evaluate_sp1_k2_law2_tray_step",
)
