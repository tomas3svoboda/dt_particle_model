"""The hexane appearance/re-evaporation law of the pre-birth condensed-hexane epoch.

Design basis: ``docs/GT_PS2_RR2_BIRTH_CONTINUATION_DESIGN_2026-08-21.md``
(section 1 E5 - the measured property anchors; section 2 D1-D3 - the epoch
law, the deposit path and the epoch dynamics; section 3 BC-1 - this build
step).  Authority: TEMPORARY RULING TR-2 under the overnight autonomy grant,
itself standing on TR-1 (``docs/GT_PS2_O9C_BIRTH_DEW_BARRIER_FINDING_2026-08-21.md``
section 3), which un-parked RR2 as the birth-continuation packet after the
birth dew barrier was measured on three independent routes.  BOTH RULINGS
AWAIT MORNING RATIFICATION, so this module is a LEAF: nothing in production
imports it, and deleting it plus its test file removes it without trace.

Lineage.  This is the A22 wedge Option 3 - "a condensed phase at the pore
wall, the virtual supersaturated state made real", recorded 2026-08-11 as
the long-term physically complete answer in
``docs/GT_PS2_A22_WEDGE_OPTION1_RULING_RECORD_2026-08-11.md`` (item 5) while
Option 1 was ruled as the bounded near-term wedge - arriving at the front
door: under the corrected donors the newborn shell is held ON the hexane
dew line, so the corrected physics does not merely permit the mobile-liquid
topology, it REQUIRES it at birth.

THE LAW.  The exact role-swapped mirror of the ruled water-appearance law
(RW-2, :mod:`dtdc_simulator.core2.qsc_water_appearance_transition_law`) with
the roles AS PRINTED in Tutkun - HEXANE condensing, water/steam
noncondensing:

    N_h   = k_y * ln((1 - y_h,I) / (1 - y_h,b))        [Tutkun Eq. (2),
            Colburn-Hougen film form; == Perre Eq. (A.7)]
    y_h,I = a_h * p_sat_h(T_I) / P,  a_h = 1           [equilibrium over its
            OWN liquid - there is NO isotherm anywhere in this law; the
            sorbed-phase machinery (PHY-031) is untouched]
    phi   = |N_h| * C_p,h / h_y                        [Tutkun Eq. (5)]
    q     = h_y * (phi / (1 - exp(-phi))) * (T_b - T_I)
            + N_h * lambda_h                           [Tutkun Eq. (4),
            Ackermann-corrected sensible + latent]

The interface temperature is the PARTICLE surface temperature under exactly
the lumped-flake convention RW-2 documents (Sipos 1961): on a lumped meal
flake the interface IS the solid surface, so Tutkun's coolant-side interface
root-find collapses to an explicit evaluation - no inner solver, no
termination floors.  Water is the noncondensing species of the immiscible
pair [Webb; Coletto Sec. 2.3.4 VLLE topology], and its interface enrichment
y_w,I = y_w,b * exp(N_h / k_y) is recorded as evidence [Webb, unnumbered].

PROPERTY AUTHORITY.  Saturation pressure and the latent enthalpy come from
the frozen pure n-hexane authority PHY-049 (Span-Wagner 2003 technical
Helmholtz EOS, :mod:`dtdc_simulator.core2.props.hexane`) and from nothing
else.  Per the design's E5 anchors the latent enthalpy is the vapor-minus-
liquid ``state_Tp`` enthalpy difference at (T_I, P); at the measured dew-line
temperature T* = 338.42486498350115 K and P = 101325 Pa that is
29025.454658883784 J/mol == 336818.4903304585 J/kg.  Those are declared CHECK
VALUES with construction-appropriate relative bounds - never bit-pins.

DOMAIN (D1).  The FULLY-WET state only.  Under frozen PHY-036 (X_f > 0 <=>
s = R) the epoch law never coexists with a receded core, so the request
carries an explicit front-position/topology assertion and refuses anything
whose front has left the particle surface.  The irreversible primary-drainage
latch is honored BY CONSTRUCTION: the request refuses any drainage history
other than PRE_CORE_RECESSION, which is precisely the predicate
:func:`dtdc_simulator.core2.particle.surface_active_set.validate_primary_drainage_admission`
enforces before attached n-hexane may exist at all.

DEPOSIT PATH (D2).  Condensate enters the EXTERNAL attached inventory - the
intensive twin ``ExternalSurfaceInventories.attached_hexane_kg_per_kg_dry`` -
through a per-step accumulator mirroring CC-2a's water-film accumulator
(:func:`dtdc_simulator.core2.particle.cut_cap_active_continuation.free_water_film_accumulator`):
exact ``math.fsum`` identity, non-negative cone, typed refusals.  The radial
activation endpoint requires ``attached_hexane == 0`` (``sphere.py:152-153``),
so the film lives on the external-inventories channel and NEVER on the feed
channel.

TWO DELIBERATE DIVERGENCES FROM RW-2 (both required by the design):

1.  SIGNED, BIDIRECTIONAL MODE (D3).  RW-2 admits the condensing direction
    only and refuses a desorbing driving force outright, because its
    desorbing direction belongs to the certified binary film.  Here the SAME
    law serves both halves of the epoch - the film grows while the flake is
    cold and SHRINKS as it heats, and D3 states plainly that re-evaporation
    is "the SAME law both ways".  The Colburn-Hougen flux therefore carries
    its sign as typed evidence (:class:`QSCHexaneEpochMode`), the Ackermann
    correction is formed on |N_h| so it is symmetric, and every fail-closed
    gate (heteroazeotrope, receded core, drainage latch, boiling limit,
    non-negative film) applies identically in both directions.  A zero
    driving force is still a typed refusal, never a silent no-op, and a step
    whose layers disagree in sign - or whose sign contradicts a caller's
    declared mode - is a typed refusal too.
2.  NO KERNEL PRE-EVENT EVIDENCE.  RW-2 consumes the certified K2 kernel's
    own ABSENT_APPEARANCE_DEMANDED NO-GO.  No such hexane-side kernel
    evidence exists (E3: there is no composition channel), so the epoch law's
    admissibility evidence is the explicit fully-wet assertion instead.  This
    also keeps the module a genuine leaf: its import graph is the frozen
    property authorities, the frozen surface active-set foundation and the
    certified heteroazeotrope solver (the G5-derived gate), and it does not
    drag the K2 kernel behind it.

MODE RESTRICTION.  Case (i) only: the interface must sit strictly above the
hexane-water heteroazeotrope temperature T_het(P), derived at the layer's
OWN pressure from the certified production solver
(:func:`dtdc_simulator.core2.heteroazeotrope.grounded_heteroazeotrope_temperature_k`;
G5 closed 2026-08-22, owner ruling G5-R1,
``docs/GT_PS2_G5_HETEROAZEOTROPE_PRESSURE_GROUNDING_RECORD_2026-08-22.md``),
else a typed refusal.  A layer pressure outside the G5-grounded band is its
own typed refusal - the gate is never clamped onto an ungrounded pressure.
The gate is SHARED with RW-2 by construction: both laws call the SAME
derived-gate function, so the sharing is machine-checked by object identity
rather than by a restated constant.  The retired provisional constant
335.15 K (Kemper, Edible Oil Processing 2nd ed., Fig. 4.6: 62 C at
760 mmHg) stays as corroboration only, in the T2b Kemper oracle band.

Every claim flag is false.  This law produces one epoch step; it advances no
accepted state, implements no birth gate (D3's typed gate is BC-2) and no
eventual-birth seam (D4), and is not physically qualifying.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import InitVar, dataclass
from enum import Enum
from typing import ClassVar

from .heteroazeotrope import (
    HeteroazeotropeError,
    grounded_heteroazeotrope_temperature_k,
)
from .particle import surface_active_set as sas
from .props import hexane as hexane_props

QSC_HEXANE_APPEARANCE_LAW_ID = (
    "colburn-hougen-film-condensation-of-hexane-onto-fully-wet-flake-epoch-v1"
)
QSC_HEXANE_APPEARANCE_SCHEMA_ID = "qsc-hexane-appearance-transition-step-v1"
QSC_HEXANE_APPEARANCE_SCHEMA_REVISION = 1
QSC_HEXANE_APPEARANCE_SOURCE_IDENTITY = "dtdc_simulator.core2.qsc_hexane_appearance_transition_law"
#: The G5 anchor pressure: the pressure at which the certified solver's
#: heteroazeotrope check value 334.4804 K was verified (G5 record section
#: 3.1).  G5 CLOSED 2026-08-22 (owner ruling G5-R1): the mode-(ii) gate is
#: no longer a constant but T_het derived at the layer's own pressure from
#: the certified heteroazeotrope solver - the design basis (D1) requires ONE
#: gate for both halves of the immiscible pair, and both laws now call the
#: SAME :func:`grounded_heteroazeotrope_temperature_k`.  See
#: ``docs/GT_PS2_G5_HETEROAZEOTROPE_PRESSURE_GROUNDING_RECORD_2026-08-22.md``.
QSC_HEXANE_APPEARANCE_REFERENCE_PRESSURE_PA = 101_325.0
#: A2 (RATIFIED 2026-08-22, ``docs/GT_PS2_BC4_COMMITTED_BIRTH_AND_BURNOFF_
#: SEAM_2026-08-22.md`` section 3): the SECOND admission axis of the SAME
#: law - the partial (post-birth) state whose newborn shell carries the
#: committed attached film of BC-4a.  The kernel is unchanged, both signs;
#: only the admission differs, so the axis carries its own schema identity.
QSC_HEXANE_POST_BIRTH_ADMISSION_AXIS_ID = (
    "post-birth-committed-film-newborn-shell-surface-admission-axis-v1"
)
QSC_HEXANE_POST_BIRTH_SCHEMA_ID = "qsc-hexane-post-birth-film-transition-step-v1"
QSC_HEXANE_POST_BIRTH_SCHEMA_REVISION = 1
#: PHY-049 molar mass of n-hexane (Span-Wagner 2003, kg/mol).
HEXANE_MOLAR_MASS_KG_MOL = hexane_props.M


class QSCHexaneEpochMode(Enum):
    """The signed direction of one epoch step (design D3)."""

    #: N_h > 0: flux toward the surface, the attached film GROWS.
    CONDENSING = "condensing"
    #: N_h < 0: flux away from the surface, the attached film SHRINKS.
    RE_EVAPORATING = "re_evaporating"


class QSCHexaneAppearanceRefusalCode(Enum):
    RECEDED_CORE_DOMAIN = "receded_core_domain"
    DRAINAGE_LATCH_VIOLATED = "drainage_latch_violated"
    MODE_SWITCH_EUTECTIC_REACHED = "mode_switch_eutectic_reached"
    INTERFACE_ABOVE_BOILING = "interface_above_boiling"
    PROPERTY_AUTHORITY_NOT_CONVERGED = "property_authority_not_converged"
    NULL_DRIVING_FORCE = "null_driving_force"
    MODE_SIGN_CONTRADICTION = "mode_sign_contradiction"
    NEGATIVE_FILM_ACCUMULATOR = "negative_film_accumulator"
    NONFINITE_INPUT = "nonfinite_input"
    NONPOSITIVE_STEP = "nonpositive_step"
    # -- A2 (RATIFIED 2026-08-22): the post-birth committed-film admission
    #    axis's own typed refusals.  Additive members only; nothing above moves.
    POST_BIRTH_COMPANION_MISSING = "post_birth_film_companion_missing"
    POST_BIRTH_COMPANION_UNBOUND = "post_birth_film_companion_unbound"
    POST_BIRTH_FRONT_NOT_RECEDED = "post_birth_front_not_receded"
    POST_BIRTH_FILM_ABSENT = "post_birth_film_absent"
    POST_BIRTH_DRAINAGE_HISTORY_MISMATCH = "post_birth_drainage_history_mismatch"
    # -- G5 (owner ruling G5-R1, 2026-08-22): the derived heteroazeotrope
    #    gate's own refusal for a layer pressure outside the grounded band.
    #    Additive member only; nothing above moves.
    HETEROAZEOTROPE_GATE_UNGROUNDED = "heteroazeotrope_gate_ungrounded_pressure"


class QSCHexaneAppearanceTransitionError(ValueError):
    """Typed refusal for the hexane appearance/re-evaporation epoch law."""

    def __init__(self, code: QSCHexaneAppearanceRefusalCode, message: str) -> None:
        super().__init__(f"{code.value}: {message}")
        self.code = code


_STEP_SEAL = object()


def _frame(value: object) -> bytes:
    if type(value) is str:
        raw = value.encode("utf-8")
        return b"s" + str(len(raw)).encode("ascii") + b":" + raw
    if type(value) is int:
        raw = str(value).encode("ascii")
        return b"i" + str(len(raw)).encode("ascii") + b":" + raw
    if type(value) is float:
        raw = value.hex().encode("ascii")
        return b"f" + str(len(raw)).encode("ascii") + b":" + raw
    raise QSCHexaneAppearanceTransitionError(
        QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
        "digest folding admits only str, int, and float leaves",
    )


def _digest(label: str, values: tuple[object, ...]) -> str:
    hasher = hashlib.sha256()
    hasher.update(label.encode("utf-8"))
    for value in values:
        hasher.update(_frame(value))
    return "sha256:" + hasher.hexdigest()


def _require_finite(label: str, value: object) -> float:
    if type(value) is not float or not math.isfinite(value):
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
            f"{label} must be a finite binary64 value",
        )
    return value


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneFullyWetAssertion:
    """The caller's explicit fully-wet front/topology assertion (design D1).

    Frozen PHY-036 reads ``X_f > 0 <=> s = R``: while free liquid is present
    the drainage front SITS AT the particle surface.  The epoch law is
    admissible on that state and no other, so the caller must state the front
    position, the free-liquid loading, the irreversible drainage history and
    the external-surface inventories the film accumulates onto.  Every
    departure is a typed refusal - the front is never nudged, the history is
    never assumed, and the inventories are never clamped.
    """

    particle_radius_m: float
    front_position_m: float
    free_liquid_loading_kg_per_kg_dry: float
    drainage_history: sas.HexanePrimaryDrainageHistory
    surface_inventories: sas.ExternalSurfaceInventories
    basis_kg_dry: float

    def __post_init__(self) -> None:
        if type(self.drainage_history) is not sas.HexanePrimaryDrainageHistory:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the primary-drainage history must be the frozen typed enum",
            )
        if type(self.surface_inventories) is not sas.ExternalSurfaceInventories:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the film basis must be typed ExternalSurfaceInventories",
            )
        for label, value in (
            ("particle radius", self.particle_radius_m),
            ("front position", self.front_position_m),
            ("free-liquid loading", self.free_liquid_loading_kg_per_kg_dry),
            ("dry-mass basis", self.basis_kg_dry),
        ):
            _require_finite(label, value)
        if self.particle_radius_m <= 0.0 or self.basis_kg_dry <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the particle radius and the dry-mass basis must be strictly positive",
            )
        if self.front_position_m > self.particle_radius_m:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the drainage front may not sit outside the particle surface",
            )
        if self.free_liquid_loading_kg_per_kg_dry <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.RECEDED_CORE_DOMAIN,
                "the epoch law is admissible on the fully-wet state only: PHY-036 "
                "reads X_f > 0 <=> s = R, and a vanished free-liquid loading is "
                "not that state; fail closed, never extrapolate",
            )
        if self.front_position_m != self.particle_radius_m:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.RECEDED_CORE_DOMAIN,
                f"the drainage front sits at {self.front_position_m} m inside a "
                f"{self.particle_radius_m} m particle: frozen PHY-036 forbids the "
                "epoch law from coexisting with a receded core",
            )
        if self.drainage_history is not sas.HexanePrimaryDrainageHistory.PRE_CORE_RECESSION:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.DRAINAGE_LATCH_VIOLATED,
                "attached n-hexane after core recession/dry-out is outside the "
                "qualified primary-drainage envelope; this is exactly the predicate "
                "surface_active_set.validate_primary_drainage_admission enforces, "
                "honored here BY CONSTRUCTION before any state is built",
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneAppearanceLayerInput:
    """One layer's gas/surface state and declared film coefficients.

    There is no sorbed loading here and no isotherm: the epoch law's interface
    is a condensed hexane liquid over which the hexane activity is exactly one.
    """

    group_id: str
    bulk_temperature_k: float
    bulk_hexane_mole_fraction: float
    pressure_pa: float
    particle_temperature_k: float
    molar_transfer_coefficient_mol_m2_s: float
    heat_transfer_coefficient_w_m2_k: float
    active_area_m2: float
    coefficient_provenance: str

    def __post_init__(self) -> None:
        if type(self.group_id) is not str or not self.group_id.strip():
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "layer group id must be nonblank",
            )
        if type(self.coefficient_provenance) is not str or not self.coefficient_provenance.strip():
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "layer coefficient provenance must be nonblank",
            )
        for label, value in (
            ("bulk temperature", self.bulk_temperature_k),
            ("bulk hexane mole fraction", self.bulk_hexane_mole_fraction),
            ("pressure", self.pressure_pa),
            ("particle temperature", self.particle_temperature_k),
            ("molar transfer coefficient", self.molar_transfer_coefficient_mol_m2_s),
            ("heat transfer coefficient", self.heat_transfer_coefficient_w_m2_k),
            ("active area", self.active_area_m2),
        ):
            _require_finite(label, value)
        if not 0.0 < self.bulk_hexane_mole_fraction < 1.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "bulk hexane mole fraction must lie strictly inside the simplex",
            )
        if self.pressure_pa <= 0.0 or self.active_area_m2 <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "pressure and active area must be strictly positive",
            )
        if (
            self.molar_transfer_coefficient_mol_m2_s <= 0.0
            or self.heat_transfer_coefficient_w_m2_k <= 0.0
        ):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "film coefficients must be strictly positive",
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneAppearanceTransitionRequest:
    """One epoch step over an explicitly asserted fully-wet flake state."""

    attempt_id: str
    fully_wet_state: QSCHexaneFullyWetAssertion
    layers: tuple[QSCHexaneAppearanceLayerInput, ...]
    step_duration_s: float
    #: Optional caller expectation.  When supplied, a computed mode that
    #: disagrees is a typed MODE_SIGN_CONTRADICTION refusal - the epoch's
    #: direction is evidence, never a knob.
    declared_mode: QSCHexaneEpochMode | None = None

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "attempt id must be nonblank",
            )
        if type(self.fully_wet_state) is not QSCHexaneFullyWetAssertion:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.RECEDED_CORE_DOMAIN,
                "the epoch law consumes the exact fully-wet assertion type",
            )
        if type(self.layers) is not tuple or not self.layers:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "at least one layer input is required",
            )
        if any(type(item) is not QSCHexaneAppearanceLayerInput for item in self.layers):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "layer inputs must be the exact layer-input type",
            )
        if self.declared_mode is not None and type(self.declared_mode) is not QSCHexaneEpochMode:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
                "a declared epoch mode must be the typed enum",
            )
        _require_finite("step duration", self.step_duration_s)
        if self.step_duration_s <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONPOSITIVE_STEP,
                "the epoch step duration must be strictly positive",
            )


class _FalseTransitionClaims:
    pre_birth_condensed_hexane_epoch_only: ClassVar[bool] = True
    fully_wet_domain_only_phy_036: ClassVar[bool] = True
    bidirectional_condensing_and_re_evaporating: ClassVar[bool] = True
    deposit_enters_external_attached_inventory_per_d2: ClassVar[bool] = True
    mode_i_hexane_condenses_water_noncondensing: ClassVar[bool] = True
    owner_ratified: ClassVar[bool] = False
    imported_by_production: ClassVar[bool] = False
    birth_gate_implemented: ClassVar[bool] = False
    eventual_birth_seam_implemented: ClassVar[bool] = False
    advances_accepted_tray_state: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneAppearanceLayerResult:
    """One layer's signed epoch fluxes and Ackermann-corrected energy split."""

    group_id: str
    hexane_activity: float
    interface_hexane_mole_fraction: float
    bulk_hexane_mole_fraction: float
    interface_saturation_pressure_pa: float
    mode: QSCHexaneEpochMode
    #: Tutkun Eq. (2), POSITIVE toward the surface (condensing).
    signed_molar_flux_mol_m2_s: float
    ackermann_rate_factor: float
    ackermann_correction: float
    sensible_heat_flux_w_m2: float
    #: Signed with the flux: released to the solid when condensing, drawn from
    #: it when the film re-evaporates.
    latent_heat_flux_w_m2: float
    latent_molar_enthalpy_j_mol: float
    latent_mass_enthalpy_j_kg: float
    signed_hexane_mass_kg: float
    signed_film_rate_kg_per_kg_dry_s: float
    energy_to_solid_j: float
    noncondensing_water_interface_enrichment: float


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneFilmAccumulator:
    """The attached-hexane film gained (or lost) across one epoch step.

    The exact mirror of CC-2a's ``FreeWaterFilmAccumulator``: the intensive
    twin is carried on the caller's declared dry-mass basis, the accumulator
    identity holds construction-exactly through ``math.fsum``, and leaving the
    non-negative cone is a typed refusal rather than a clamp.  A hexane film
    step may not move the free water, exactly as a water film step may not
    move the attached hexane.
    """

    before: sas.ExternalSurfaceInventories
    after: sas.ExternalSurfaceInventories
    mode: QSCHexaneEpochMode
    basis_kg_dry: float
    duration_s: float
    rate_kg_per_kg_dry_s: float
    deposit_kg_per_kg_dry: float
    deposit_kg: float
    deposit_mol: float
    deposit_energy_j: float
    molar_enthalpy_j_mol: float
    exact_accumulator_identity: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for payload in (self.before, self.after):
            if type(payload) is not sas.ExternalSurfaceInventories:
                raise QSCHexaneAppearanceTransitionError(
                    QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                    "the film accumulator needs typed ExternalSurfaceInventories",
                )
        if not math.isfinite(self.basis_kg_dry) or self.basis_kg_dry <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the film accumulator basis must be positive and finite",
            )
        if self.after.free_water_kg_per_kg_dry != self.before.free_water_kg_per_kg_dry:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "a hexane film step may not move the free water",
            )
        if self.after.attached_hexane_kg_per_kg_dry != math.fsum(
            (self.before.attached_hexane_kg_per_kg_dry, self.deposit_kg_per_kg_dry)
        ):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the film accumulator identity must hold construction-exactly",
            )
        if self.deposit_kg != self.deposit_kg_per_kg_dry * self.basis_kg_dry:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the extensive film deposit must be the intensive gain times basis",
            )
        if self.after.attached_hexane_kg_per_kg_dry < 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR,
                "the external attached-hexane accumulator may never go negative",
            )
        if (self.deposit_kg_per_kg_dry > 0.0) is not (self.mode is QSCHexaneEpochMode.CONDENSING):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
                "the accumulated film sign must agree with the epoch mode",
            )

    @property
    def topology_after(self) -> sas.SurfaceTopology:
        return self.after.topology


def hexane_film_accumulator(
    film_before: sas.ExternalSurfaceInventories,
    *,
    mode: QSCHexaneEpochMode,
    basis_kg_dry: float,
    duration_s: float,
    rate_kg_per_kg_dry_s: float,
    molar_enthalpy_j_mol: float,
) -> QSCHexaneFilmAccumulator:
    """Return the attached-hexane film gained over one epoch step, exactly."""

    deposit_intensive = rate_kg_per_kg_dry_s * duration_s
    deposit_kg = deposit_intensive * basis_kg_dry
    try:
        film_after = sas.ExternalSurfaceInventories(
            math.fsum((film_before.attached_hexane_kg_per_kg_dry, deposit_intensive)),
            film_before.free_water_kg_per_kg_dry,
        )
    except ValueError as exc:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR,
            f"the external attached-hexane accumulator left its non-negative cone: {exc}; "
            "the film cannot re-evaporate more hexane than it holds - refuse, never clamp",
        ) from exc
    deposit_mol = deposit_kg / HEXANE_MOLAR_MASS_KG_MOL
    return QSCHexaneFilmAccumulator(
        before=film_before,
        after=film_after,
        mode=mode,
        basis_kg_dry=basis_kg_dry,
        duration_s=duration_s,
        rate_kg_per_kg_dry_s=rate_kg_per_kg_dry_s,
        deposit_kg_per_kg_dry=deposit_intensive,
        deposit_kg=deposit_kg,
        deposit_mol=deposit_mol,
        deposit_energy_j=deposit_mol * molar_enthalpy_j_mol,
        molar_enthalpy_j_mol=molar_enthalpy_j_mol,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedQSCHexaneAppearanceTransitionStep(_FalseTransitionClaims):
    """The sealed lawful epoch step on an asserted fully-wet flake state."""

    request: QSCHexaneAppearanceTransitionRequest
    layer_results: tuple[QSCHexaneAppearanceLayerResult, ...]
    mode: QSCHexaneEpochMode
    film: QSCHexaneFilmAccumulator
    net_hexane_mass_change_kg: float
    total_energy_to_solid_j: float
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_APPEARANCE_SCHEMA_ID
    schema_revision: ClassVar[int] = QSC_HEXANE_APPEARANCE_SCHEMA_REVISION
    law_id: ClassVar[str] = QSC_HEXANE_APPEARANCE_LAW_ID

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _STEP_SEAL:
            raise TypeError("hexane epoch transition steps are law-sealed")

    @property
    def definition_digest(self) -> str:
        state = self.request.fully_wet_state
        leaves: list[object] = [
            QSC_HEXANE_APPEARANCE_SCHEMA_ID,
            QSC_HEXANE_APPEARANCE_SCHEMA_REVISION,
            QSC_HEXANE_APPEARANCE_LAW_ID,
            self.request.attempt_id,
            self.request.step_duration_s,
            self.mode.value,
            state.particle_radius_m,
            state.front_position_m,
            state.free_liquid_loading_kg_per_kg_dry,
            state.drainage_history.value,
            state.surface_inventories.attached_hexane_kg_per_kg_dry,
            state.surface_inventories.free_water_kg_per_kg_dry,
            state.basis_kg_dry,
            self.film.deposit_kg_per_kg_dry,
            self.film.after.attached_hexane_kg_per_kg_dry,
            self.net_hexane_mass_change_kg,
            self.total_energy_to_solid_j,
        ]
        for result in self.layer_results:
            leaves.extend(
                (
                    result.group_id,
                    result.hexane_activity,
                    result.interface_hexane_mole_fraction,
                    result.signed_molar_flux_mol_m2_s,
                    result.signed_hexane_mass_kg,
                    result.energy_to_solid_j,
                )
            )
        return _digest("qsc-hexane-appearance-transition-step", tuple(leaves))


def _interface_hexane_fraction(layer: QSCHexaneAppearanceLayerInput) -> tuple[float, float, float]:
    """Return (a_h, p_sat, y_h,I) at the layer's own interface temperature."""

    interface_temperature = layer.particle_temperature_k
    # G5 closed 2026-08-22 (owner ruling G5-R1): the gate is T_het derived at
    # the layer's own pressure from the certified heteroazeotrope solver; see
    # docs/GT_PS2_G5_HETEROAZEOTROPE_PRESSURE_GROUNDING_RECORD_2026-08-22.md.
    try:
        heteroazeotrope_temperature = grounded_heteroazeotrope_temperature_k(
            layer.pressure_pa
        )
    except HeteroazeotropeError as exc:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.HETEROAZEOTROPE_GATE_UNGROUNDED,
            "the mode-(ii) gate cannot be derived at the layer pressure "
            f"{layer.pressure_pa} Pa: {exc}",
        ) from exc
    if interface_temperature <= heteroazeotrope_temperature:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.MODE_SWITCH_EUTECTIC_REACHED,
            "the interface reached the derived hexane-water heteroazeotrope "
            f"({heteroazeotrope_temperature} K at {layer.pressure_pa} Pa, certified "
            "solver, G5 closed 2026-08-22): mode (ii) is not implemented; typed "
            "refusal in BOTH epoch directions, never a silent branch",
        )
    if interface_temperature >= hexane_props.TC:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.INTERFACE_ABOVE_BOILING,
            f"no hexane saturation locus above Tc = {hexane_props.TC} K",
        )
    saturation = hexane_props.saturation(interface_temperature)
    if not saturation.converged or not math.isfinite(saturation.p):
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.PROPERTY_AUTHORITY_NOT_CONVERGED,
            f"the frozen PHY-049 saturation solve did not converge at {interface_temperature} K",
        )
    # a_h = 1 EXACTLY: the interface is the law's own condensed hexane liquid.
    # There is no isotherm on this axis and none is invented here.
    activity = 1.0
    interface_fraction = activity * saturation.p / layer.pressure_pa
    if interface_fraction >= 1.0:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.INTERFACE_ABOVE_BOILING,
            f"p_sat_h({interface_temperature} K) = {saturation.p} Pa reaches the layer "
            f"pressure {layer.pressure_pa} Pa: the y_dew locus does not exist at or "
            "above the boiling point (the T_max bracketing limit measured in the "
            "march era); fail closed rather than leave the simplex",
        )
    if interface_fraction <= 0.0:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
            "interface hexane fraction left the simplex",
        )
    return activity, saturation.p, interface_fraction


def evaluate_qsc_hexane_appearance_transition_step(
    request: QSCHexaneAppearanceTransitionRequest,
) -> ValidatedQSCHexaneAppearanceTransitionStep:
    """Evaluate one pre-birth epoch step; refuse anything off the design basis."""

    state = request.fully_wet_state
    results: list[QSCHexaneAppearanceLayerResult] = []
    modes: set[QSCHexaneEpochMode] = set()
    for layer in request.layers:
        interface_temperature = layer.particle_temperature_k
        activity, saturation_pressure, interface_fraction = _interface_hexane_fraction(layer)
        flux = layer.molar_transfer_coefficient_mol_m2_s * math.log(
            (1.0 - interface_fraction) / (1.0 - layer.bulk_hexane_mole_fraction)
        )
        # Sign convention: Tutkun's N_1 is POSITIVE toward the surface.  With
        # the bulk above the interface the log ratio exceeds one and hexane
        # condenses; with the bulk below it the same expression turns negative
        # and the SAME law re-evaporates the film (design D3).
        if flux == 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NULL_DRIVING_FORCE,
                "the bulk hexane fraction sits exactly on the interface dew value: "
                "no epoch step exists; refuse rather than emit a silent no-op",
            )
        mode = QSCHexaneEpochMode.CONDENSING if flux > 0.0 else QSCHexaneEpochMode.RE_EVAPORATING
        modes.add(mode)
        magnitude = abs(flux)
        # Ackermann on |N_h|: the rate factor is a magnitude, so the correction
        # is symmetric under the epoch's sign reversal (deliberate divergence 1).
        cp_molar = hexane_props.cp0(interface_temperature)
        rate_factor = magnitude * cp_molar / layer.heat_transfer_coefficient_w_m2_k
        # expm1 keeps the Ackermann factor exact at small rate factors, where
        # 1 - exp(-phi) cancels catastrophically.
        correction = rate_factor / (-math.expm1(-rate_factor))
        sensible = (
            layer.heat_transfer_coefficient_w_m2_k
            * correction
            * (layer.bulk_temperature_k - interface_temperature)
        )
        # E5: the latent enthalpy is the frozen authority's vapor-minus-liquid
        # state_Tp enthalpy difference at the interface state.
        vapor = hexane_props.state_Tp(interface_temperature, layer.pressure_pa, "vapor")
        liquid = hexane_props.state_Tp(interface_temperature, layer.pressure_pa, "liquid")
        latent_mass = vapor.h_mass - liquid.h_mass
        latent_molar = latent_mass * HEXANE_MOLAR_MASS_KG_MOL
        latent = flux * latent_molar
        signed_mass = (
            flux * layer.active_area_m2 * request.step_duration_s * HEXANE_MOLAR_MASS_KG_MOL
        )
        energy = (sensible + latent) * layer.active_area_m2 * request.step_duration_s
        # Webb: the noncondensing partner is enriched at a condensing interface
        # and depleted at a re-evaporating one - the signed flux carries both.
        enrichment = (1.0 - layer.bulk_hexane_mole_fraction) * math.exp(
            flux / layer.molar_transfer_coefficient_mol_m2_s
        )
        results.append(
            QSCHexaneAppearanceLayerResult(
                group_id=layer.group_id,
                hexane_activity=activity,
                interface_hexane_mole_fraction=interface_fraction,
                bulk_hexane_mole_fraction=layer.bulk_hexane_mole_fraction,
                interface_saturation_pressure_pa=saturation_pressure,
                mode=mode,
                signed_molar_flux_mol_m2_s=flux,
                ackermann_rate_factor=rate_factor,
                ackermann_correction=correction,
                sensible_heat_flux_w_m2=sensible,
                latent_heat_flux_w_m2=latent,
                latent_molar_enthalpy_j_mol=latent_molar,
                latent_mass_enthalpy_j_kg=latent_mass,
                signed_hexane_mass_kg=signed_mass,
                signed_film_rate_kg_per_kg_dry_s=(
                    signed_mass / state.basis_kg_dry / request.step_duration_s
                ),
                energy_to_solid_j=energy,
                noncondensing_water_interface_enrichment=enrichment,
            )
        )
    if len(modes) != 1:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
            "the epoch step's layers disagree in flux sign: one step is one mode, "
            "so split the step rather than net two directions together",
        )
    (mode,) = modes
    if request.declared_mode is not None and request.declared_mode is not mode:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
            f"the caller declared {request.declared_mode.value} but the Colburn-Hougen "
            f"driving force is {mode.value}; the epoch's direction is evidence, not a knob",
        )
    net_mass = math.fsum(item.signed_hexane_mass_kg for item in results)
    total_energy = math.fsum(item.energy_to_solid_j for item in results)
    weight = math.fsum(abs(item.signed_hexane_mass_kg) for item in results)
    # A reported aggregate only: the |mass|-weighted mean latent enthalpy of
    # the layers that fed the film.  No closure rests on it.
    molar_enthalpy = (
        math.fsum(
            abs(item.signed_hexane_mass_kg) * item.latent_molar_enthalpy_j_mol for item in results
        )
        / weight
    )
    film = hexane_film_accumulator(
        state.surface_inventories,
        mode=mode,
        basis_kg_dry=state.basis_kg_dry,
        duration_s=request.step_duration_s,
        rate_kg_per_kg_dry_s=math.fsum(item.signed_film_rate_kg_per_kg_dry_s for item in results),
        molar_enthalpy_j_mol=molar_enthalpy,
    )
    return ValidatedQSCHexaneAppearanceTransitionStep(
        request=request,
        layer_results=tuple(results),
        mode=mode,
        film=film,
        net_hexane_mass_change_kg=net_mass,
        total_energy_to_solid_j=total_energy,
        _seal=_STEP_SEAL,
    )


# ---------------------------------------------------------------------------
# A2 - THE SECOND ADMISSION AXIS (RATIFIED 2026-08-22)
#
# Authority: ``docs/GT_PS2_BC4_COMMITTED_BIRTH_AND_BURNOFF_SEAM_2026-08-22.md``
# section 3 A2, RULED block ("I ratify A2 and A3").  BC-4b measured, on the
# committed corrected-physics birth (BC-4a, commit a1ff942), that the newborn
# state carries an attached hexane film with the front honestly RECEDED - and
# that every admission axis of this law refused it typed (GAP 1).  A2 extends
# the admission with a typed SECOND axis: the partial (post-birth) state whose
# newborn shell carries the COMMITTED attached film, free of the fully-wet
# front requirement on this axis and gated instead on the film companion's
# presence and identity binding.
#
# THE KERNEL IS UNCHANGED - BOTH SIGNS.  The measured post-birth surface is
# CONDENSING (y_I = p_sat(T_s)/P = 0.882 < y_b = 0.90 at the committed
# references; the sign flips at 338.5473743733437 K): the film GROWS first and
# burns off later, through exactly the bidirectional Colburn-Hougen kernel the
# fully-wet axis evaluates.  The evaluation below repeats that kernel verbatim
# - same expressions, same operand order - so the two axes are one law, and
# the fully-wet axis's own entry point is untouched (additive change only).
#
# THE DRAINAGE-LATCH INTERPLAY, STATED.  The film on the newborn state is the
# PRE-RECESSION epoch's continuation (A1's amendment scope, enforced at the
# BC-4a booking: the D3 birth gate required the epoch film at EXACTLY zero, so
# the newborn film is entirely the birth step's own condensate).  The state's
# drainage history is CORE_RECESSION_OCCURRED - the birth happened - and
# FURTHER recession stays refused until the film is EXACTLY zero: this axis
# evolves the film on a shell whose front is pinned, it never moves a front,
# and at film == 0.0 the axis closes (a zero film is a typed refusal here, not
# an admission - the latch clearance belongs to the ordinary machinery).
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexanePostBirthFilmAssertion:
    """The typed newborn-shell-surface assertion (A2's second admission axis).

    Admissible exactly on the committed post-birth state: the drainage front
    strictly INSIDE the particle, the attached-hexane film strictly POSITIVE,
    the history CORE_RECESSION_OCCURRED, and the committed-film companion
    PRESENT and IDENTITY-BOUND to the declared committed state (BC-4a's
    ``is_bound_to`` discipline) with the companion's film and the asserted
    inventories one number, bit for bit.  Every departure is a typed refusal -
    the fully-wet state belongs to the first axis, a dried-out core is outside
    A1's scope, and a vanished film means the latch has cleared.
    """

    particle_radius_m: float
    front_position_m: float
    drainage_history: sas.HexanePrimaryDrainageHistory
    surface_inventories: sas.ExternalSurfaceInventories
    basis_kg_dry: float
    committed_state: object
    film_companion: object

    #: A1's amendment scope: the newborn film is the pre-recession epoch's
    #: continuation - the D3 gate certified the epoch film at EXACTLY zero, so
    #: this film is entirely the birth step's own committed condensate.
    film_is_the_pre_recession_epochs_continuation: ClassVar[bool] = True
    #: The latch interplay: this axis is PRE-RECESSION with respect to any
    #: FURTHER front motion - recession stays refused until film == 0.0, and
    #: the law on this axis never moves a front.
    recession_stays_refused_until_film_exactly_zero: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.drainage_history) is not sas.HexanePrimaryDrainageHistory:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the primary-drainage history must be the frozen typed enum",
            )
        if type(self.surface_inventories) is not sas.ExternalSurfaceInventories:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the film basis must be typed ExternalSurfaceInventories",
            )
        for label, value in (
            ("particle radius", self.particle_radius_m),
            ("front position", self.front_position_m),
            ("dry-mass basis", self.basis_kg_dry),
        ):
            _require_finite(label, value)
        if self.particle_radius_m <= 0.0 or self.basis_kg_dry <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the particle radius and the dry-mass basis must be strictly positive",
            )
        if self.front_position_m > self.particle_radius_m:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "the drainage front may not sit outside the particle surface",
            )
        if self.front_position_m == self.particle_radius_m:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_FRONT_NOT_RECEDED,
                "the front sits AT the particle surface: that is the fully-wet "
                "state and the FIRST admission axis owns it; the post-birth "
                "axis admits only an honestly receded front",
            )
        if self.front_position_m <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_FRONT_NOT_RECEDED,
                "the post-birth axis requires a strictly interior drainage "
                "front: without a wet core behind the newborn shell there is "
                "no strict-partial state to evolve a film on",
            )
        if self.surface_inventories.attached_hexane_kg_per_kg_dry <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_FILM_ABSENT,
                "the post-birth axis exists only while a committed attached "
                "film covers the newborn surface; at EXACTLY zero the drainage "
                "latch clears and the ordinary machinery owns the state - "
                "refuse here, never re-open the axis on a vanished film",
            )
        if self.drainage_history is not (
            sas.HexanePrimaryDrainageHistory.CORE_RECESSION_OCCURRED
        ):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_DRAINAGE_HISTORY_MISMATCH,
                f"the post-birth axis admits history "
                f"{sas.HexanePrimaryDrainageHistory.CORE_RECESSION_OCCURRED.value!r} "
                f"only, got {self.drainage_history.value!r}: PRE_CORE_RECESSION "
                "belongs to the fully-wet axis and a dried-out core is outside "
                "A1's committed-film scope",
            )
        if self.committed_state is None:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_COMPANION_UNBOUND,
                "the post-birth axis needs the committed state the film "
                "companion is identity-bound to",
            )
        if self.film_companion is None:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_COMPANION_MISSING,
                "the post-birth axis is scoped to the A1 film's continuation: "
                "the committed-film companion must be present",
            )
        bound = getattr(self.film_companion, "is_bound_to", None)
        if not callable(bound):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_COMPANION_MISSING,
                "the film companion carries no is_bound_to identity predicate; "
                "an unbindable companion is no companion",
            )
        if bound(self.committed_state) is not True:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_COMPANION_UNBOUND,
                "the film companion is not identity-bound to the declared "
                "committed state (BC-4a's is_bound_to discipline): the record "
                "can neither be re-attached to another state nor survive a "
                "state it does not describe",
            )
        companion_film = getattr(
            self.film_companion, "attached_hexane_kg_per_kg_dry", None
        )
        if companion_film != self.surface_inventories.attached_hexane_kg_per_kg_dry:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_COMPANION_UNBOUND,
                "the companion's attached film and the asserted inventories "
                f"must be one number, bit for bit: companion {companion_film!r} "
                "vs inventories "
                f"{self.surface_inventories.attached_hexane_kg_per_kg_dry!r}",
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexanePostBirthTransitionRequest:
    """One epoch-law step over the asserted post-birth committed-film state."""

    attempt_id: str
    post_birth_state: QSCHexanePostBirthFilmAssertion
    layers: tuple[QSCHexaneAppearanceLayerInput, ...]
    step_duration_s: float
    #: Optional caller expectation, exactly the fully-wet axis's discipline:
    #: a computed mode that disagrees is a typed MODE_SIGN_CONTRADICTION.
    declared_mode: QSCHexaneEpochMode | None = None

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "attempt id must be nonblank",
            )
        if type(self.post_birth_state) is not QSCHexanePostBirthFilmAssertion:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.POST_BIRTH_FRONT_NOT_RECEDED,
                "the post-birth axis consumes the exact post-birth film "
                "assertion type; the fully-wet assertion belongs to the first "
                "axis and no other payload is admitted",
            )
        if type(self.layers) is not tuple or not self.layers:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "at least one layer input is required",
            )
        if any(type(item) is not QSCHexaneAppearanceLayerInput for item in self.layers):
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONFINITE_INPUT,
                "layer inputs must be the exact layer-input type",
            )
        if self.declared_mode is not None and type(self.declared_mode) is not QSCHexaneEpochMode:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
                "a declared epoch mode must be the typed enum",
            )
        _require_finite("step duration", self.step_duration_s)
        if self.step_duration_s <= 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NONPOSITIVE_STEP,
                "the epoch step duration must be strictly positive",
            )


class _FalsePostBirthClaims:
    post_birth_committed_film_axis_only: ClassVar[bool] = True
    same_kernel_as_the_fully_wet_axis_both_signs: ClassVar[bool] = True
    film_is_the_pre_recession_epochs_continuation: ClassVar[bool] = True
    recession_stays_refused_until_film_exactly_zero: ClassVar[bool] = True
    a2_owner_ratified_2026_08_22: ClassVar[bool] = True
    imported_by_production: ClassVar[bool] = False
    birth_gate_implemented: ClassVar[bool] = False
    advances_accepted_tray_state: ClassVar[bool] = False
    moves_the_drainage_front: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedQSCHexanePostBirthTransitionStep(_FalsePostBirthClaims):
    """The sealed lawful epoch-law step on the post-birth committed-film state."""

    request: QSCHexanePostBirthTransitionRequest
    layer_results: tuple[QSCHexaneAppearanceLayerResult, ...]
    mode: QSCHexaneEpochMode
    film: QSCHexaneFilmAccumulator
    net_hexane_mass_change_kg: float
    total_energy_to_solid_j: float
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_POST_BIRTH_SCHEMA_ID
    schema_revision: ClassVar[int] = QSC_HEXANE_POST_BIRTH_SCHEMA_REVISION
    #: The SAME law: the kernel is BC-1's, unchanged, both signs.
    law_id: ClassVar[str] = QSC_HEXANE_APPEARANCE_LAW_ID
    admission_axis_id: ClassVar[str] = QSC_HEXANE_POST_BIRTH_ADMISSION_AXIS_ID

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _STEP_SEAL:
            raise TypeError("hexane epoch transition steps are law-sealed")

    @property
    def definition_digest(self) -> str:
        state = self.request.post_birth_state
        leaves: list[object] = [
            QSC_HEXANE_POST_BIRTH_SCHEMA_ID,
            QSC_HEXANE_POST_BIRTH_SCHEMA_REVISION,
            QSC_HEXANE_POST_BIRTH_ADMISSION_AXIS_ID,
            QSC_HEXANE_APPEARANCE_LAW_ID,
            self.request.attempt_id,
            self.request.step_duration_s,
            self.mode.value,
            state.particle_radius_m,
            state.front_position_m,
            state.drainage_history.value,
            state.surface_inventories.attached_hexane_kg_per_kg_dry,
            state.surface_inventories.free_water_kg_per_kg_dry,
            state.basis_kg_dry,
            self.film.deposit_kg_per_kg_dry,
            self.film.after.attached_hexane_kg_per_kg_dry,
            self.net_hexane_mass_change_kg,
            self.total_energy_to_solid_j,
        ]
        for result in self.layer_results:
            leaves.extend(
                (
                    result.group_id,
                    result.hexane_activity,
                    result.interface_hexane_mole_fraction,
                    result.signed_molar_flux_mol_m2_s,
                    result.signed_hexane_mass_kg,
                    result.energy_to_solid_j,
                )
            )
        return _digest("qsc-hexane-post-birth-film-transition-step", tuple(leaves))


def evaluate_qsc_hexane_post_birth_transition_step(
    request: QSCHexanePostBirthTransitionRequest,
) -> ValidatedQSCHexanePostBirthTransitionStep:
    """Evaluate one epoch-law step on the post-birth axis (A2).

    THE KERNEL BELOW IS BC-1'S, VERBATIM - the same expressions in the same
    operand order as ``evaluate_qsc_hexane_appearance_transition_step``, both
    signs, so the two admission axes evaluate one law bit for bit.  Only the
    admission (already enforced by the typed request/assertion constructors)
    differs.  The fully-wet entry point is untouched: this is an additive
    second door into the same room, and the suite asserts the bit-identity of
    the two doors' layer results on identical layer inputs.
    """

    state = request.post_birth_state
    results: list[QSCHexaneAppearanceLayerResult] = []
    modes: set[QSCHexaneEpochMode] = set()
    for layer in request.layers:
        interface_temperature = layer.particle_temperature_k
        activity, saturation_pressure, interface_fraction = _interface_hexane_fraction(layer)
        flux = layer.molar_transfer_coefficient_mol_m2_s * math.log(
            (1.0 - interface_fraction) / (1.0 - layer.bulk_hexane_mole_fraction)
        )
        if flux == 0.0:
            raise QSCHexaneAppearanceTransitionError(
                QSCHexaneAppearanceRefusalCode.NULL_DRIVING_FORCE,
                "the bulk hexane fraction sits exactly on the interface dew value: "
                "no epoch step exists; refuse rather than emit a silent no-op",
            )
        mode = QSCHexaneEpochMode.CONDENSING if flux > 0.0 else QSCHexaneEpochMode.RE_EVAPORATING
        modes.add(mode)
        magnitude = abs(flux)
        cp_molar = hexane_props.cp0(interface_temperature)
        rate_factor = magnitude * cp_molar / layer.heat_transfer_coefficient_w_m2_k
        correction = rate_factor / (-math.expm1(-rate_factor))
        sensible = (
            layer.heat_transfer_coefficient_w_m2_k
            * correction
            * (layer.bulk_temperature_k - interface_temperature)
        )
        vapor = hexane_props.state_Tp(interface_temperature, layer.pressure_pa, "vapor")
        liquid = hexane_props.state_Tp(interface_temperature, layer.pressure_pa, "liquid")
        latent_mass = vapor.h_mass - liquid.h_mass
        latent_molar = latent_mass * HEXANE_MOLAR_MASS_KG_MOL
        latent = flux * latent_molar
        signed_mass = (
            flux * layer.active_area_m2 * request.step_duration_s * HEXANE_MOLAR_MASS_KG_MOL
        )
        energy = (sensible + latent) * layer.active_area_m2 * request.step_duration_s
        enrichment = (1.0 - layer.bulk_hexane_mole_fraction) * math.exp(
            flux / layer.molar_transfer_coefficient_mol_m2_s
        )
        results.append(
            QSCHexaneAppearanceLayerResult(
                group_id=layer.group_id,
                hexane_activity=activity,
                interface_hexane_mole_fraction=interface_fraction,
                bulk_hexane_mole_fraction=layer.bulk_hexane_mole_fraction,
                interface_saturation_pressure_pa=saturation_pressure,
                mode=mode,
                signed_molar_flux_mol_m2_s=flux,
                ackermann_rate_factor=rate_factor,
                ackermann_correction=correction,
                sensible_heat_flux_w_m2=sensible,
                latent_heat_flux_w_m2=latent,
                latent_molar_enthalpy_j_mol=latent_molar,
                latent_mass_enthalpy_j_kg=latent_mass,
                signed_hexane_mass_kg=signed_mass,
                signed_film_rate_kg_per_kg_dry_s=(
                    signed_mass / state.basis_kg_dry / request.step_duration_s
                ),
                energy_to_solid_j=energy,
                noncondensing_water_interface_enrichment=enrichment,
            )
        )
    if len(modes) != 1:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
            "the epoch step's layers disagree in flux sign: one step is one mode, "
            "so split the step rather than net two directions together",
        )
    (mode,) = modes
    if request.declared_mode is not None and request.declared_mode is not mode:
        raise QSCHexaneAppearanceTransitionError(
            QSCHexaneAppearanceRefusalCode.MODE_SIGN_CONTRADICTION,
            f"the caller declared {request.declared_mode.value} but the Colburn-Hougen "
            f"driving force is {mode.value}; the epoch's direction is evidence, not a knob",
        )
    net_mass = math.fsum(item.signed_hexane_mass_kg for item in results)
    total_energy = math.fsum(item.energy_to_solid_j for item in results)
    weight = math.fsum(abs(item.signed_hexane_mass_kg) for item in results)
    molar_enthalpy = (
        math.fsum(
            abs(item.signed_hexane_mass_kg) * item.latent_molar_enthalpy_j_mol for item in results
        )
        / weight
    )
    film = hexane_film_accumulator(
        state.surface_inventories,
        mode=mode,
        basis_kg_dry=state.basis_kg_dry,
        duration_s=request.step_duration_s,
        rate_kg_per_kg_dry_s=math.fsum(item.signed_film_rate_kg_per_kg_dry_s for item in results),
        molar_enthalpy_j_mol=molar_enthalpy,
    )
    return ValidatedQSCHexanePostBirthTransitionStep(
        request=request,
        layer_results=tuple(results),
        mode=mode,
        film=film,
        net_hexane_mass_change_kg=net_mass,
        total_energy_to_solid_j=total_energy,
        _seal=_STEP_SEAL,
    )


__all__ = [
    "HEXANE_MOLAR_MASS_KG_MOL",
    "QSC_HEXANE_APPEARANCE_LAW_ID",
    "QSC_HEXANE_APPEARANCE_REFERENCE_PRESSURE_PA",
    "QSC_HEXANE_APPEARANCE_SCHEMA_ID",
    "QSC_HEXANE_APPEARANCE_SCHEMA_REVISION",
    "QSC_HEXANE_APPEARANCE_SOURCE_IDENTITY",
    "QSC_HEXANE_POST_BIRTH_ADMISSION_AXIS_ID",
    "QSC_HEXANE_POST_BIRTH_SCHEMA_ID",
    "QSC_HEXANE_POST_BIRTH_SCHEMA_REVISION",
    "QSCHexaneAppearanceLayerInput",
    "QSCHexaneAppearanceLayerResult",
    "QSCHexaneAppearanceRefusalCode",
    "QSCHexaneAppearanceTransitionError",
    "QSCHexaneAppearanceTransitionRequest",
    "QSCHexaneEpochMode",
    "QSCHexaneFilmAccumulator",
    "QSCHexaneFullyWetAssertion",
    "QSCHexanePostBirthFilmAssertion",
    "QSCHexanePostBirthTransitionRequest",
    "ValidatedQSCHexaneAppearanceTransitionStep",
    "ValidatedQSCHexanePostBirthTransitionStep",
    "evaluate_qsc_hexane_appearance_transition_step",
    "evaluate_qsc_hexane_post_birth_transition_step",
    "hexane_film_accumulator",
]
