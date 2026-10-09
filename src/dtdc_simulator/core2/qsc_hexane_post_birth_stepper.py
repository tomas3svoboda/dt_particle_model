"""BC-5b: the composed post-birth stepper - film growth, burn-off, latch clear.

AUTHORITY.  ``docs/GT_PS2_BC4_COMMITTED_BIRTH_AND_BURNOFF_SEAM_2026-08-22.md``,
RULED block 2026-08-22 ("I ratify A2 and A3"): with the two ratified
admissions built (BC-5a, commit 5409d4c), this module composes them into the
seam BC-4b measured missing - the post-birth film evolution on the committed
corrected-physics birth state (BC-4a, commit a1ff942), mirroring the BC-2
epoch pattern (:mod:`dtdc_simulator.core2.qsc_hexane_epoch_prefix`) one state
family later.

THE COMPOSED INTERVAL.  One post-birth interval over one duration is:

  (i)   the A2 film-law step at the newborn surface - the BC-1 kernel on the
        SECOND admission axis (``evaluate_qsc_hexane_post_birth_transition_
        step``), evaluated at the START-of-interval surface state exactly as
        the BC-2 epoch's declared explicit-surface split.  The interface
        temperature is the accepted state's own reduced-film surface root
        (the surface film audit's ``surface_temperature_k`` - the same
        reference the BC-4b sign probe measured at T* - 0.508 K), carried on
        the film companion;
  (ii)  the A3 film-aware interior stall advance over the SAME interval - the
        pinned-front stall on the condensed-solvent boundary variant with a
        FRESH sweep journal and a FRESH single-admission channel per solve,
        under the typed :class:`cut_interior_stall.FilmCoveredStallAdmission`.
        The front pin is exactly the drainage latch's requirement while
        ``film > 0``;
  (iii) the film companion carried forward with the accumulator discipline:
        the A2 accumulator's own ``after`` inventories, identity-bound to the
        stall's committed state;
  (iv)  the two lawful bookings of the surface hexane throughput REPORTED
        beside one another (the A2 Colburn-Hougen deposit vs the stall
        ledger's boundary outflow), never reconciled by force - exactly
        BC-4a's "this record says where that outflow LIVES" discipline.

THE TYPED FILM-EXTINCTION EVENT.  The drainage latch clears at film EXACTLY
zero (float equality), so the burn-off's last interval is a LOCATED event in
the BC-2 discipline: the film is linear in time over the interval, the
analytic root is bracketed to one binary64 duration, the endpoint inventory
is CONSTRUCTED at exactly zero and any sub-resolution remainder is REPORTED
as a balance residual - never a clamp, never an epsilon film.  The ULP window
and probe-halving budgets are the BC-2 epoch's own, imported so this module
declares no bracket budget of its own.

THE RULED RECESSION BRANCH (A5, RULED 2026-08-22; BC-6b).  The frozen PHY-014
text carries its own carve-out sentence - "existing admissible full-core or
attached-liquid states may evolve normally" - and the A5 ruling
(``docs/GT_PS2_A5_LATCH_CARVEOUT_RULING_PACKET_2026-08-22.md`` section 3)
completed the drainage latch to it: recession is ADMISSIBLE while the
attached film is a lawfully committed, identity-bound film of the A1/A2
lineage.  A composed interval whose SA-1 demand classifies RECESSION
therefore no longer refuses: the lineage is verified through the latch's own
:func:`surface_active_set.require_committed_film_lineage`, and the interval
runs the ordinary moving-front march entry
(``cut_front_demand.advance_with_typed_front_demand`` - the recession arm of
``cut_interior_stall.advance_with_stall_admission``) WITH the film companion
carried and evolving: the A2 law step composes per interval in both signs,
the front recedes, and the carrier re-binds to the committed state.  The
interior water balance drying the flake while its surface is wet with
solvent is ordinary simultaneous heat/mass transfer, and every quantity in
it stays typed, committed, and mass-closed.

THE LATCH CLEARANCE AND THE HAND-OFF.  Once the carrier's film is exactly
zero, :func:`require_post_birth_latch_clear` seals the typed clearance record
and the state belongs to the ordinary stall-aware march entry
(``cut_interior_stall.advance_with_stall_admission``) - this module never
takes a post-clearance step.  Under A5 the clearance semantics are: film ==
0.0 ENDS THE FILM PHASE - it no longer gates recession, which the ruled
branch above admits under the committed film.  An ADVANCE demand raises
SA-1's own :class:`cut_front_demand.PrimaryDrainageAdvanceDemanded`
untouched; a film without the typed lineage keeps the frozen refusal.  All
first-class findings, not defects.

Every claim flag reachable from this module is False.  Nothing here makes the
particle model physically qualifying, plant predictive, or production wired.
"""

from __future__ import annotations

import math
from dataclasses import InitVar, dataclass, field
from enum import Enum
from typing import ClassVar

from . import qsc_hexane_appearance_transition_law as bc1
from . import qsc_hexane_epoch_prefix as bc2
from .particle import condensed_solvent_reduced_film as csrf
from .particle import coupled_transport as ct
from .particle import cut_birth_condensed_film as bc4
from .particle import cut_event_orchestrator as ceo
from .particle import cut_front_demand as cfd
from .particle import cut_integrator as ci
from .particle import cut_interior_stall as cis
from .particle import surface_active_set as sas

QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_ID = "qsc-hexane-post-birth-composed-interval-v1"
QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_REVISION = 1
QSC_HEXANE_POST_BIRTH_STEPPER_ORCHESTRATION_ID = (
    "post-birth-film-growth-burnoff-latch-clear-composed-stepper-v1"
)

#: SHARED ORIGIN: the BC-2 epoch's own located-event budgets, imported rather
#: than restated so this module declares no bracket budget of its own.
FILM_EXTINCTION_ULP_WINDOW = bc2.QSC_HEXANE_EPOCH_EXHAUSTION_ULP_WINDOW
FILM_EXTINCTION_PROBE_HALVINGS = bc2.QSC_HEXANE_EPOCH_EXHAUSTION_PROBE_HALVINGS


class QSCHexanePostBirthStepperRefusalCode(Enum):
    NONFINITE_INPUT = "nonfinite_input"
    NONPOSITIVE_STEP = "nonpositive_step"
    NOT_A_COMMITTED_BIRTH_FILM = "not_a_committed_birth_film"
    CARRIER_NOT_TYPED = "carrier_not_typed"
    CARRIER_UNBOUND = "carrier_unbound"
    BOUNDARY_NOT_PLAIN = "boundary_not_plain_reduced_film"
    LAW_REFUSED = "post_birth_law_refused"
    FILM_ALREADY_EXACTLY_ZERO = "film_already_exactly_zero"
    FILM_EXTINCTION_OVERSHOOT = "film_extinction_overshoot"
    FILM_EXTINCTION_NOT_EXACT = "film_extinction_not_exact"
    FILM_EXTINCTION_BEYOND_BOUND = "film_extinction_beyond_bound"
    NO_FILM_TO_EXTINGUISH = "no_film_to_extinguish"
    FILM_NOT_RE_EVAPORATING = "film_not_re_evaporating"
    #: Retired as a raise site by A5/BC-6b (RULED 2026-08-22): the ruled
    #: recession branch proceeds under a lineage-verified committed film.
    #: The member stays for record readers; nothing raises it any more.
    RECESSION_DEMANDED_UNDER_FILM = "recession_demanded_under_film"
    SURFACE_AUDIT_MISSING = "surface_audit_missing"
    LATCH_NOT_CLEAR = "latch_not_clear"
    #: A5/BC-6b: the ruled recession branch's own typed refusals.
    FILM_LINEAGE_UNVERIFIED = "film_lineage_unverified_for_recession"
    RECESSION_SEED_MISSING = "recession_under_film_same_cell_seed_missing"
    FACE_ROUTE_UNDER_FILM = "face_event_under_film_not_composed"


class QSCHexanePostBirthStepperError(ValueError):
    """Typed refusal for the composed post-birth stepper; exact rollback.

    No entry point in this module ever mutates its inputs, so a caller that
    catches this has lost nothing.
    """

    def __init__(
        self,
        code: QSCHexanePostBirthStepperRefusalCode,
        message: str,
        *,
        state: object | None = None,
        evidence: dict | None = None,
    ) -> None:
        super().__init__(f"{code.value}: {message}")
        self.code = code
        self.state = state
        self.evidence = dict(evidence) if evidence else {}


_STEPPER_SEAL = object()


class _FalseStepperClaims:
    post_birth_film_seam_only: ClassVar[bool] = True
    composes_the_two_ratified_admissions_a2_a3: ClassVar[bool] = True
    explicit_surface_implicit_interior_split: ClassVar[bool] = True
    front_pinned_under_film_per_the_latch: ClassVar[bool] = True
    bc2_located_event_discipline: ClassVar[bool] = True
    moves_the_drainage_front: ClassVar[bool] = False
    advances_accepted_tray_state: ClassVar[bool] = False
    imported_by_production: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False


def _require_finite(label: str, value: object) -> float:
    if type(value) is not float or not math.isfinite(value):
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NONFINITE_INPUT,
            f"{label} must be a finite binary64 value",
        )
    return value


def _require_seed_typed(same_cell_seed: object) -> None:
    """A5/BC-6b: a supplied recession seed must be the exact typed seed."""

    if same_cell_seed is not None and not isinstance(
        same_cell_seed, ci.CutStepSeed
    ):
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NONFINITE_INPUT,
            "the ruled recession branch consumes a typed cut_integrator."
            f"CutStepSeed; got {type(same_cell_seed).__name__}",
        )


# ---------------------------------------------------------------------------
# The carried film companion
# ---------------------------------------------------------------------------
@dataclass(frozen=True, kw_only=True)
class PostBirthFilmCarrier(_FalseStepperClaims):
    """The committed film companion carried across the post-birth intervals.

    The BC-4a sealed record binds the birth condensate to the birth-committed
    state BY IDENTITY; each accepted interval re-binds the running film to the
    stall's committed state through this carrier, with the accumulator
    discipline (the A2 accumulator's own ``after`` inventories, cumulative
    counters via ``math.fsum``).  ``is_bound_to`` and
    ``attached_hexane_kg_per_kg_dry`` are exactly the two facts the A2
    assertion and the A3 admission gate on, so one carrier serves both doors.
    """

    committed_state: object
    inventories: sas.ExternalSurfaceInventories
    basis_kg_dry: float
    #: The accepted reduced-film surface root the NEXT law step evaluates at.
    surface_temperature_k: float
    #: Lineage: the BC-4a sealed birth record this film continues.
    birth_film_record: object
    clock_s: float
    interval_index: int
    cumulative_condensed_kg: float
    cumulative_burned_off_kg: float
    cumulative_film_energy_j: float

    def __post_init__(self) -> None:
        if type(self.inventories) is not sas.ExternalSurfaceInventories:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.CARRIER_NOT_TYPED,
                "the carrier needs typed ExternalSurfaceInventories",
            )
        for label, value in (
            ("carrier basis", self.basis_kg_dry),
            ("carrier surface temperature", self.surface_temperature_k),
            ("carrier clock", self.clock_s),
            ("carrier cumulative condensed", self.cumulative_condensed_kg),
            ("carrier cumulative burned off", self.cumulative_burned_off_kg),
            ("carrier cumulative film energy", self.cumulative_film_energy_j),
        ):
            _require_finite(label, value)
        if self.basis_kg_dry <= 0.0 or self.clock_s < 0.0:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.NONFINITE_INPUT,
                "the carrier basis must be positive and the clock non-negative",
            )
        if (
            isinstance(self.interval_index, bool)
            or type(self.interval_index) is not int
            or self.interval_index < 0
        ):
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.NONFINITE_INPUT,
                "the carrier interval index must be a non-negative exact int",
            )
        if self.committed_state is None:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.CARRIER_UNBOUND,
                "the carrier must be bound to a committed state",
            )
        if self.birth_film_record is None:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.NOT_A_COMMITTED_BIRTH_FILM,
                "the carrier must carry its BC-4a birth-record lineage",
            )

    def is_bound_to(self, state: object) -> bool:
        return state is self.committed_state

    @property
    def attached_hexane_kg_per_kg_dry(self) -> float:
        return self.inventories.attached_hexane_kg_per_kg_dry

    @property
    def film_kg(self) -> float:
        return self.attached_hexane_kg_per_kg_dry * self.basis_kg_dry

    def journal_record(self) -> dict:
        return {
            "schema_id": QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_ID,
            "record": "post-birth-film-carrier",
            "interval_index": self.interval_index,
            "clock_s": self.clock_s,
            "film_kg_per_kg_dry": self.attached_hexane_kg_per_kg_dry,
            "film_kg": self.film_kg,
            "free_water_kg_per_kg_dry": self.inventories.free_water_kg_per_kg_dry,
            "surface_topology": self.inventories.topology.value,
            "basis_kg_dry": self.basis_kg_dry,
            "surface_temperature_k": self.surface_temperature_k,
            "cumulative_condensed_kg": self.cumulative_condensed_kg,
            "cumulative_burned_off_kg": self.cumulative_burned_off_kg,
            "cumulative_film_energy_j": self.cumulative_film_energy_j,
            "physically_qualifying": False,
        }


def initialize_post_birth_film_state(
    birth_record: bc4.CommittedBirthAttachedFilm,
) -> PostBirthFilmCarrier:
    """Enter the post-birth seam on the BC-4a sealed film record."""

    if type(birth_record) is not bc4.CommittedBirthAttachedFilm:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NOT_A_COMMITTED_BIRTH_FILM,
            "the post-birth seam opens on exactly one sealed "
            f"CommittedBirthAttachedFilm; got {type(birth_record).__name__}",
        )
    committed = birth_record.committed_state
    if not isinstance(committed, ci.CutIntegratorState):
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NOT_A_COMMITTED_BIRTH_FILM,
            "the birth record's committed state is not a cut integrator state",
        )
    if not birth_record.is_bound_to(committed):
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_UNBOUND,
            "the birth record does not bind its own committed state",
        )
    film = birth_record.film.after
    if film.attached_hexane_kg_per_kg_dry <= 0.0:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.FILM_ALREADY_EXACTLY_ZERO,
            "the committed birth booked no positive film; there is no "
            "post-birth film seam to step",
        )
    return PostBirthFilmCarrier(
        committed_state=committed,
        inventories=film,
        basis_kg_dry=birth_record.film.basis_kg_dry,
        surface_temperature_k=birth_record.admission.surface_temperature_k,
        birth_film_record=birth_record,
        clock_s=0.0,
        interval_index=0,
        cumulative_condensed_kg=0.0,
        cumulative_burned_off_kg=0.0,
        cumulative_film_energy_j=0.0,
    )


# ---------------------------------------------------------------------------
# The composed interval
# ---------------------------------------------------------------------------
@dataclass(frozen=True, kw_only=True)
class QSCHexanePostBirthIntervalStep(_FalseStepperClaims):
    """One accepted composed interval: the A2 law step + the interior advance.

    The interior advance is EXACTLY ONE of the two lawful branches (the
    stall-aware march's own dichotomy, one layer earlier): the A3 film-covered
    pinned-front stall (``stall_step``), or - A5/BC-6b, RULED 2026-08-22 - the
    ordinary moving-front macrostep under the lineage-verified committed film
    (``macrostep``).  ``branch`` names which; the carrier is re-bound to the
    committed state either way.
    """

    before_state: ci.CutIntegratorState
    after_state: ci.CutIntegratorState
    before_carrier: PostBirthFilmCarrier
    after_carrier: PostBirthFilmCarrier
    law_step: bc1.ValidatedQSCHexanePostBirthTransitionStep
    stall_step: cis.InteriorStallStep | None = None
    #: A5/BC-6b: the ordinary moving-front macrostep taken under the film.
    macrostep: ceo.CutEventMacrostep | None = None
    branch: cis.StallAwareBranch = cis.StallAwareBranch.INTERIOR_STALL
    #: The SA-1 verdict that selected the branch, carried as evidence.
    entry_front_demand: cfd.FrontDemandEvidence | None = None
    duration_s: float
    mode: bc1.QSCHexaneEpochMode
    film_before_kg_per_kg_dry: float
    film_after_kg_per_kg_dry: float
    film_extinction_event: bool
    #: Strictly zero except on the located FILM EXTINCTION EVENT, where it is
    #: the sub-resolution remainder the endpoint construction retired -
    #: reported exactly as the BC-2 epoch reports its own.
    film_extinction_residual_kg_per_kg_dry: float
    surface_temperature_before_k: float
    surface_temperature_after_k: float
    condensed_root_admitted: bool
    #: The two lawful bookings of the surface hexane throughput, REPORTED
    #: beside one another and never reconciled by force (BC-4a's discipline).
    law_deposit_mol: float
    stall_boundary_hexane_out_mol: float
    film_vs_boundary_association_mol: float
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_ID
    schema_revision: ClassVar[int] = QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_REVISION
    orchestration_id: ClassVar[str] = QSC_HEXANE_POST_BIRTH_STEPPER_ORCHESTRATION_ID
    explicit_surface_coupling: ClassVar[bool] = True

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _STEPPER_SEAL:
            raise TypeError("post-birth composed intervals are stepper-sealed")
        if (self.stall_step is None) == (self.macrostep is None):
            raise ValueError(
                "a composed interval carries exactly one committed interior "
                "advance: the A3 stall or the A5-ruled recession macrostep"
            )
        if self.branch is cis.StallAwareBranch.INTERIOR_STALL:
            if self.stall_step is None:
                raise ValueError("the stall branch needs its sealed stall step")
        elif self.macrostep is None:
            raise ValueError("the recession branch needs its committed macrostep")

    def journal_record(self) -> dict:
        record = self.after_carrier.journal_record()
        record.update(
            {
                "record": "post-birth-composed-interval",
                "branch": self.branch.value,
                "dt_s": self.duration_s,
                "mode": self.mode.value,
                "film_before_kg_per_kg_dry": self.film_before_kg_per_kg_dry,
                "film_after_kg_per_kg_dry": self.film_after_kg_per_kg_dry,
                "film_extinction_event": self.film_extinction_event,
                "film_extinction_residual_kg_per_kg_dry": (
                    self.film_extinction_residual_kg_per_kg_dry
                ),
                "surface_temperature_before_k": self.surface_temperature_before_k,
                "surface_temperature_after_k": self.surface_temperature_after_k,
                "condensed_root_admitted": self.condensed_root_admitted,
                "signed_molar_flux_mol_m2_s": (
                    self.law_step.layer_results[0].signed_molar_flux_mol_m2_s
                ),
                "interface_hexane_mole_fraction": (
                    self.law_step.layer_results[0].interface_hexane_mole_fraction
                ),
                "law_energy_to_solid_j": self.law_step.total_energy_to_solid_j,
                "law_deposit_mol": self.law_deposit_mol,
                "stall_boundary_hexane_out_mol": self.stall_boundary_hexane_out_mol,
                "film_vs_boundary_association_mol": (
                    self.film_vs_boundary_association_mol
                ),
                "front_z": self.after_state.transport.geometry.front.z,
                "explicit_surface_coupling": True,
                "physically_qualifying": False,
            }
        )
        if self.entry_front_demand is not None:
            record.update(
                {
                    "entry_q_demand_m3_s": (
                        self.entry_front_demand.demanded_front_volume_rate_m3_s
                    ),
                    "entry_resolution_m3_s": self.entry_front_demand.resolution_m3_s,
                    "entry_normalized_demand": (
                        self.entry_front_demand.normalized_demand
                    ),
                }
            )
        if self.stall_step is not None:
            ledger = self.stall_step.ledger
            record.update(
                {
                    "front_pinned_bit_exactly": ledger.front_pinned_bit_exactly,
                    "max_independent_scaled_residual": (
                        ledger.maximum_independent_scaled_residual
                    ),
                    "normalized_stall_consistency_row": (
                        ledger.normalized_stall_consistency_row
                    ),
                    "max_step_ledger_residual": ledger.maximum_step_ledger_residual,
                    "max_cumulative_ledger_residual": (
                        ledger.maximum_cumulative_ledger_residual
                    ),
                    "condition_proxy": ledger.condition_proxy,
                    "entry_q_demand_m3_s": (
                        ledger.entry_demanded_front_volume_rate_m3_s
                    ),
                    "entry_resolution_m3_s": ledger.entry_resolution_m3_s,
                }
            )
        else:
            macro = self.macrostep.ledger
            record.update(
                {
                    "front_pinned_bit_exactly": False,
                    "front_z_before": (
                        self.before_state.transport.geometry.front.z
                    ),
                    "macro_branch": getattr(
                        macro.branch, "value", str(macro.branch)
                    ),
                    "max_stage_scaled_residual": (
                        macro.maximum_stage_scaled_residual
                    ),
                    "max_step_ledger_residual": (
                        macro.maximum_stage_step_ledger_residual
                    ),
                    "max_cumulative_ledger_residual": (
                        macro.maximum_final_cumulative_ledger_residual
                    ),
                }
            )
        return record


def _law_step(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    duration_s: float,
    declared_mode: bc1.QSCHexaneEpochMode | None,
) -> bc1.ValidatedQSCHexanePostBirthTransitionStep:
    """Evaluate the A2 axis on the current newborn-shell surface state."""

    geometry = state.transport.geometry
    assertion = bc1.QSCHexanePostBirthFilmAssertion(
        particle_radius_m=geometry.master_grid.R,
        front_position_m=geometry.front.radius_m,
        drainage_history=sas.HexanePrimaryDrainageHistory.CORE_RECESSION_OCCURRED,
        surface_inventories=carrier.inventories,
        basis_kg_dry=carrier.basis_kg_dry,
        committed_state=state,
        film_companion=carrier,
    )
    layer = bc1.QSCHexaneAppearanceLayerInput(
        group_id=f"post-birth-surface-interval-{carrier.interval_index}",
        bulk_temperature_k=gas_boundary.bulk_temperature_k,
        bulk_hexane_mole_fraction=gas_boundary.bulk_hexane_mole_fraction,
        pressure_pa=gas_boundary.pressure_pa,
        particle_temperature_k=carrier.surface_temperature_k,
        molar_transfer_coefficient_mol_m2_s=(
            gas_boundary.molar_transfer_coefficient_mol_m2_s
        ),
        heat_transfer_coefficient_w_m2_k=gas_boundary.heat_transfer_coefficient_w_m2_k,
        active_area_m2=geometry.master_grid.areas[-1],
        coefficient_provenance=gas_boundary.coefficient_provenance,
    )
    request = bc1.QSCHexanePostBirthTransitionRequest(
        attempt_id=f"post-birth-interval-{carrier.interval_index}",
        post_birth_state=assertion,
        layers=(layer,),
        step_duration_s=duration_s,
        declared_mode=declared_mode,
    )
    return bc1.evaluate_qsc_hexane_post_birth_transition_step(request)


def _gate_entry(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    plain_boundary: ct.ReducedFilmPoreBoundary,
    duration_s: float,
) -> None:
    if not isinstance(state, ci.CutIntegratorState):
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_UNBOUND,
            "the composed interval advances a cut integrator state",
        )
    if type(carrier) is not PostBirthFilmCarrier:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_NOT_TYPED,
            "the composed interval consumes the exact typed film carrier",
        )
    if carrier.is_bound_to(state) is not True:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_UNBOUND,
            "the film carrier is not identity-bound to the state it would "
            "advance; a film record can neither be re-attached to another "
            "state nor survive a state it does not describe",
            state=state,
        )
    if type(gas_boundary) is not bc2.QSCHexaneEpochBoundary:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NONFINITE_INPUT,
            "the composed interval consumes the exact typed epoch boundary "
            "(the same declared bulk-gas pair BC-2 stepped the epoch on)",
        )
    if type(plain_boundary) is not ct.ReducedFilmPoreBoundary:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.BOUNDARY_NOT_PLAIN,
            "the composed interval promotes an exact plain reduced-film pore "
            "boundary per solve (fresh journal + fresh channel); got "
            f"{type(plain_boundary).__name__}",
        )
    _require_finite("interval duration", duration_s)
    if duration_s <= 0.0:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NONPOSITIVE_STEP,
            "the interval duration must be strictly positive",
            state=state,
        )
    if carrier.attached_hexane_kg_per_kg_dry == 0.0:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.FILM_ALREADY_EXACTLY_ZERO,
            "the film is EXACTLY zero: the drainage latch has cleared and the "
            "state belongs to the ordinary stall-aware march entry, not to "
            "this stepper",
            state=state,
        )


def _composed_interval(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    law: bc1.ValidatedQSCHexanePostBirthTransitionStep,
    plain_boundary: ct.ReducedFilmPoreBoundary,
    duration_s: float,
    *,
    stall_controls: cis.InteriorStallControls | None,
    journal_record_capacity: int,
    film_extinction_event: bool,
    film_extinction_residual_kg_per_kg_dry: float,
    same_cell_seed: ci.CutStepSeed | None = None,
) -> QSCHexanePostBirthIntervalStep:
    """Fold one accepted A2 law step into one lawful interior advance.

    The SA-1 verdict selects the branch: STALL runs the A3 film-covered
    pinned-front stall; RECESSION (A5, RULED 2026-08-22 - the ruled latch
    carve-out) runs the ordinary moving-front march entry under the
    lineage-verified committed film; ADVANCE raises SA-1's typed refusal
    untouched.
    """

    journal = csrf.CondensedSolventSweepJournal(
        record_capacity=journal_record_capacity
    )
    channel = csrf.CondensedRootCommitChannel()
    promoted = csrf.promote_reduced_film_boundary(
        plain_boundary,
        sweep_journal=journal,
        condensed_root_commit_channel=channel,
    )
    evidence = cfd.evaluate_front_demand(state, promoted, duration_s)
    cfd.refuse_on_advance_demand(
        evidence,
        state,
        surfaced_at=(
            "qsc_hexane_post_birth_stepper._composed_interval "
            "(film-covered entry; A2/A3 RATIFIED 2026-08-22, A5 RULED "
            "2026-08-22)"
        ),
    )
    stall: cis.InteriorStallStep | None = None
    macrostep: ceo.CutEventMacrostep | None = None
    if evidence.outcome is cfd.FrontDemandOutcome.RECESSION_DEMAND:
        # A5/BC-6b: THE RULED RECESSION BRANCH.  The latch carve-out admits
        # recession exactly when the attached film is the lawfully committed,
        # identity-bound A1/A2-lineage film - verified through the latch's own
        # public lineage predicate, never re-implemented here.
        try:
            sas.require_committed_film_lineage(
                sas.CommittedFilmLineageEvidence(
                    film_companion=carrier, committed_state=state
                ),
                carrier.attached_hexane_kg_per_kg_dry,
            )
        except sas.PrimaryDrainageEnvelopeError as exc:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.FILM_LINEAGE_UNVERIFIED,
                "the assembled balances demand recession but the attached "
                f"film's A1/A2 lineage did not verify: {exc}; a film without "
                "the typed committed lineage keeps the frozen latch refusal",
                state=state,
                evidence={
                    "q_demand_m3_s": evidence.demanded_front_volume_rate_m3_s,
                    "resolution_m3_s": evidence.resolution_m3_s,
                    "film_kg_per_kg_dry": carrier.attached_hexane_kg_per_kg_dry,
                },
            ) from exc
        if same_cell_seed is None:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.RECESSION_SEED_MISSING,
                "the ruled recession branch runs the ordinary moving-front "
                "macrostep, which needs the caller's typed same-cell seed; "
                "pass same_cell_seed (the march controller owns seeds)",
                state=state,
                evidence={
                    "q_demand_m3_s": evidence.demanded_front_volume_rate_m3_s,
                    "resolution_m3_s": evidence.resolution_m3_s,
                    "normalized_demand": evidence.normalized_demand,
                    "probe_duration_s": evidence.probe_duration_s,
                },
            )
        macrostep = cfd.advance_with_typed_front_demand(
            state,
            duration_s,
            promoted,
            same_cell_seed,
        )
        if macrostep.same_cell_step is None or not isinstance(
            macrostep.after, ci.CutIntegratorState
        ):
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.FACE_ROUTE_UNDER_FILM,
                "the recession macrostep left the same-cell branch "
                f"({macrostep.ledger.branch.value!r}): a face-arrival event "
                "under the committed film has no composed seam yet; this is "
                "a first-class typed halt, never worked around",
                state=state,
                evidence={"macro_branch": macrostep.ledger.branch.value},
            )
        after_state = macrostep.after
        underlying = macrostep.same_cell_step
        audit = getattr(underlying, "surface_film_audit", None)
        if audit is None:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.SURFACE_AUDIT_MISSING,
                "the accepted recession macrostep carries no surface film "
                "audit; the next law step's interface reference is missing",
                state=state,
            )
        boundary_out_mol = underlying.ledger.boundary_hexane_out_mol
        condensed_admitted = channel.admission is not None
    else:
        admission = cis.FilmCoveredStallAdmission(
            film_companion=carrier,
            attached_hexane_kg_per_kg_dry=carrier.attached_hexane_kg_per_kg_dry,
        )
        stall = cis.solve_interior_stall(
            state,
            duration_s,
            promoted,
            controls=stall_controls,
            entry_evidence=evidence,
            film_covered_admission=admission,
        )
        after_state = stall.after
        audit = stall.step.surface_film_audit
        if audit is None:
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.SURFACE_AUDIT_MISSING,
                "the accepted film-covered stall carries no surface film "
                "audit; the next law step's interface reference is missing",
                state=state,
            )
        boundary_out_mol = stall.step.ledger.boundary_hexane_out_mol
        condensed_admitted = stall.condensed_root_admission is not None
    if film_extinction_residual_kg_per_kg_dry == 0.0:
        inventories_after = law.film.after
    else:
        # THE LOCATED FILM EXTINCTION EVENT (BC-2's discipline): the endpoint
        # is CONSTRUCTED at exactly zero and the remainder is REPORTED.
        inventories_after = sas.ExternalSurfaceInventories(
            0.0, law.film.after.free_water_kg_per_kg_dry
        )
    realized_kg = (
        inventories_after.attached_hexane_kg_per_kg_dry
        - carrier.attached_hexane_kg_per_kg_dry
    ) * carrier.basis_kg_dry
    condensed = carrier.cumulative_condensed_kg
    burned_off = carrier.cumulative_burned_off_kg
    if law.mode is bc1.QSCHexaneEpochMode.CONDENSING:
        condensed = math.fsum((condensed, realized_kg))
    else:
        burned_off = math.fsum((burned_off, -realized_kg))
    after_carrier = PostBirthFilmCarrier(
        committed_state=after_state,
        inventories=inventories_after,
        basis_kg_dry=carrier.basis_kg_dry,
        surface_temperature_k=audit.surface_temperature_k,
        birth_film_record=carrier.birth_film_record,
        clock_s=math.fsum((carrier.clock_s, duration_s)),
        interval_index=carrier.interval_index + 1,
        cumulative_condensed_kg=condensed,
        cumulative_burned_off_kg=burned_off,
        cumulative_film_energy_j=math.fsum(
            (carrier.cumulative_film_energy_j, law.film.deposit_energy_j)
        ),
    )
    law_deposit_mol = law.film.deposit_mol
    return QSCHexanePostBirthIntervalStep(
        before_state=state,
        after_state=after_state,
        before_carrier=carrier,
        after_carrier=after_carrier,
        law_step=law,
        stall_step=stall,
        macrostep=macrostep,
        branch=(
            cis.StallAwareBranch.INTERIOR_STALL
            if stall is not None
            else cis.StallAwareBranch.ORDINARY_RECESSION
        ),
        entry_front_demand=evidence,
        duration_s=duration_s,
        mode=law.mode,
        film_before_kg_per_kg_dry=carrier.attached_hexane_kg_per_kg_dry,
        film_after_kg_per_kg_dry=(
            inventories_after.attached_hexane_kg_per_kg_dry
        ),
        film_extinction_event=film_extinction_event,
        film_extinction_residual_kg_per_kg_dry=(
            film_extinction_residual_kg_per_kg_dry
        ),
        surface_temperature_before_k=carrier.surface_temperature_k,
        surface_temperature_after_k=audit.surface_temperature_k,
        condensed_root_admitted=condensed_admitted,
        law_deposit_mol=law_deposit_mol,
        stall_boundary_hexane_out_mol=boundary_out_mol,
        # The two lawful bookings of one surface throughput, associated by
        # fsum and REPORTED - a positive law deposit against a negative
        # boundary outflow cancels exactly when the two books agree.
        film_vs_boundary_association_mol=math.fsum(
            (law_deposit_mol, boundary_out_mol)
        ),
        _seal=_STEPPER_SEAL,
    )


def probe_post_birth_film_law(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    duration_s: float,
    *,
    declared_mode: bc1.QSCHexaneEpochMode | None = None,
) -> bc1.ValidatedQSCHexanePostBirthTransitionStep:
    """Evaluate ONE A2 law step on the current state without taking any step.

    A diagnostic-only public entry: the same typed assertion and the same
    kernel the composed interval evaluates, with no stall advance, no film
    commit and no state change - the probe's film accumulator is evidence,
    never adopted.  Every gate of the A2 axis applies unchanged.
    """

    if not isinstance(state, ci.CutIntegratorState):
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_UNBOUND,
            "the film-law probe evaluates on a cut integrator state",
        )
    if type(carrier) is not PostBirthFilmCarrier:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_NOT_TYPED,
            "the film-law probe consumes the exact typed film carrier",
        )
    if carrier.is_bound_to(state) is not True:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_UNBOUND,
            "the film carrier is not identity-bound to the probed state",
            state=state,
        )
    if type(gas_boundary) is not bc2.QSCHexaneEpochBoundary:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NONFINITE_INPUT,
            "the film-law probe consumes the exact typed epoch boundary",
        )
    return _law_step(state, carrier, gas_boundary, duration_s, declared_mode)


def advance_post_birth_interval(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    plain_boundary: ct.ReducedFilmPoreBoundary,
    duration_s: float,
    *,
    stall_controls: cis.InteriorStallControls | None = None,
    journal_record_capacity: int = 64,
    declared_mode: bc1.QSCHexaneEpochMode | None = None,
    same_cell_seed: ci.CutStepSeed | None = None,
) -> QSCHexanePostBirthIntervalStep:
    """Advance one composed post-birth interval; refuse anything off-basis.

    A law step that would overdraw the film is a typed
    ``FILM_EXTINCTION_OVERSHOOT`` refusal carrying the analytic extinction
    duration - take :func:`advance_post_birth_to_film_extinction` instead.
    ``same_cell_seed`` (A5/BC-6b) is the caller's typed seed for the ruled
    recession branch's ordinary macrostep; with ``None`` a RECESSION demand
    is a typed ``RECESSION_SEED_MISSING`` refusal - the branch never invents
    a seed of its own.
    """

    _gate_entry(state, carrier, gas_boundary, plain_boundary, duration_s)
    _require_seed_typed(same_cell_seed)
    try:
        law = _law_step(state, carrier, gas_boundary, duration_s, declared_mode)
    except bc1.QSCHexaneAppearanceTransitionError as exc:
        if exc.code is bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR:
            evidence = _extinction_evidence(
                state, carrier, gas_boundary, duration_s
            )
            raise QSCHexanePostBirthStepperError(
                QSCHexanePostBirthStepperRefusalCode.FILM_EXTINCTION_OVERSHOOT,
                "this interval would burn off more hexane than the film holds "
                f"({exc}); refuse rather than clamp - take the located "
                "extinction step instead (advance_post_birth_to_film_"
                "extinction), whose analytic duration at this state is "
                f"{evidence.get('extinction_duration_s')!r} s",
                state=state,
                evidence=evidence,
            ) from exc
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.LAW_REFUSED,
            f"the A2 post-birth law refused this interval: {exc}",
            state=state,
            evidence={"law_refusal_code": exc.code.value, "duration_s": duration_s},
        ) from exc
    return _composed_interval(
        state,
        carrier,
        law,
        plain_boundary,
        duration_s,
        stall_controls=stall_controls,
        journal_record_capacity=journal_record_capacity,
        film_extinction_event=False,
        film_extinction_residual_kg_per_kg_dry=0.0,
        same_cell_seed=same_cell_seed,
    )


# ---------------------------------------------------------------------------
# The located film-extinction event (BC-2's discipline)
# ---------------------------------------------------------------------------
def _probe_film_rate(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    maximum_duration_s: float,
) -> tuple[bc1.ValidatedQSCHexanePostBirthTransitionStep, float]:
    """Read the film rate off the A2 law through a bounded duration bracket."""

    duration = maximum_duration_s
    for _ in range(FILM_EXTINCTION_PROBE_HALVINGS):
        if not duration > 0.0:
            break
        try:
            return _law_step(state, carrier, gas_boundary, duration, None), duration
        except bc1.QSCHexaneAppearanceTransitionError as exc:
            if exc.code is not (
                bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR
            ):
                raise
            duration = duration / 2.0
    raise QSCHexanePostBirthStepperError(
        QSCHexanePostBirthStepperRefusalCode.FILM_EXTINCTION_NOT_EXACT,
        "the film rate could not be read within "
        f"{FILM_EXTINCTION_PROBE_HALVINGS} halvings of {maximum_duration_s} s; "
        "refuse rather than guess a rate",
        state=state,
    )


def _extinction_evidence(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    probe_duration_s: float,
) -> dict:
    """Report the analytic extinction duration without taking any step."""

    film = carrier.attached_hexane_kg_per_kg_dry
    try:
        probe, _ = _probe_film_rate(state, carrier, gas_boundary, probe_duration_s)
    except (QSCHexanePostBirthStepperError, bc1.QSCHexaneAppearanceTransitionError):
        return {"film_kg_per_kg_dry": film, "extinction_duration_s": None}
    rate = probe.film.rate_kg_per_kg_dry_s
    if rate >= 0.0:
        return {
            "film_kg_per_kg_dry": film,
            "film_rate_kg_per_kg_dry_s": rate,
            "extinction_duration_s": None,
        }
    return {
        "film_kg_per_kg_dry": film,
        "film_rate_kg_per_kg_dry_s": rate,
        "extinction_duration_s": film / -rate,
    }


def _bracket_extinction_duration(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    analytic: float,
) -> tuple[float, bc1.ValidatedQSCHexanePostBirthTransitionStep] | None:
    """Return the largest duration whose accumulator stays non-negative."""

    duration = analytic
    accepted: (
        tuple[float, bc1.ValidatedQSCHexanePostBirthTransitionStep] | None
    ) = None
    for _ in range(FILM_EXTINCTION_ULP_WINDOW):
        if not duration > 0.0:
            return None
        try:
            accepted = (
                duration,
                _law_step(state, carrier, gas_boundary, duration, None),
            )
            break
        except bc1.QSCHexaneAppearanceTransitionError as exc:
            if exc.code is not (
                bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR
            ):
                raise
            duration = math.nextafter(duration, 0.0)
    if accepted is None:
        return None
    for _ in range(FILM_EXTINCTION_ULP_WINDOW):
        if accepted[1].film.after.attached_hexane_kg_per_kg_dry == 0.0:
            return accepted
        candidate = math.nextafter(accepted[0], math.inf)
        try:
            law = _law_step(state, carrier, gas_boundary, candidate, None)
        except bc1.QSCHexaneAppearanceTransitionError as exc:
            if exc.code is not (
                bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR
            ):
                raise
            return accepted
        accepted = (candidate, law)
    return None


def advance_post_birth_to_film_extinction(
    state: ci.CutIntegratorState,
    carrier: PostBirthFilmCarrier,
    gas_boundary: bc2.QSCHexaneEpochBoundary,
    plain_boundary: ct.ReducedFilmPoreBoundary,
    *,
    maximum_duration_s: float,
    stall_controls: cis.InteriorStallControls | None = None,
    journal_record_capacity: int = 64,
    same_cell_seed: ci.CutStepSeed | None = None,
) -> QSCHexanePostBirthIntervalStep:
    """Take the one composed interval that lands the film on EXACTLY zero.

    The film-extinction event stays a LOCATED typed event under A5: analytic
    root, one-ULP binary64 bracket, endpoint CONSTRUCTED at exactly zero with
    the sub-resolution remainder reported.  At film == 0.0 the FILM PHASE
    ends (the A5 semantics - it no longer gates recession).  The interior
    advance runs over the SAME located duration on whichever lawful branch
    the SA-1 demand selects - the A3 stall or, with ``same_cell_seed``
    supplied, the ruled recession macrostep - so the state and the film
    share one clock.
    """

    _gate_entry(state, carrier, gas_boundary, plain_boundary, maximum_duration_s)
    _require_seed_typed(same_cell_seed)
    film = carrier.attached_hexane_kg_per_kg_dry
    if film == 0.0:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.NO_FILM_TO_EXTINGUISH,
            "the attached film is already exactly zero; there is no "
            "extinction interval to take",
            state=state,
        )
    probe, _probe_duration = _probe_film_rate(
        state, carrier, gas_boundary, maximum_duration_s
    )
    if probe.mode is not bc1.QSCHexaneEpochMode.RE_EVAPORATING:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.FILM_NOT_RE_EVAPORATING,
            "the Colburn-Hougen driving force is CONDENSING at this surface: "
            "the film is growing, not burning off; the direction is evidence, "
            "not a knob",
            state=state,
            evidence={"mode": probe.mode.value},
        )
    rate = probe.film.rate_kg_per_kg_dry_s
    analytic = film / -rate
    if analytic > maximum_duration_s:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.FILM_EXTINCTION_BEYOND_BOUND,
            f"the analytic extinction duration {analytic} s exceeds the "
            f"declared bound {maximum_duration_s} s; step further rather than "
            "widen the bound implicitly",
            state=state,
            evidence={
                "extinction_duration_s": analytic,
                "maximum_duration_s": maximum_duration_s,
                "film_kg_per_kg_dry": film,
                "film_rate_kg_per_kg_dry_s": rate,
            },
        )
    bracket = _bracket_extinction_duration(state, carrier, gas_boundary, analytic)
    if bracket is None:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.FILM_EXTINCTION_NOT_EXACT,
            "the film-extinction root could not be bracketed in binary64 "
            f"within {FILM_EXTINCTION_ULP_WINDOW} durations of {analytic} s; "
            "refuse rather than construct an endpoint on an unlocated root",
            state=state,
            evidence={
                "extinction_duration_s": analytic,
                "film_kg_per_kg_dry": film,
                "film_rate_kg_per_kg_dry_s": rate,
            },
        )
    duration, law = bracket
    residue = law.film.after.attached_hexane_kg_per_kg_dry
    return _composed_interval(
        state,
        carrier,
        law,
        plain_boundary,
        duration,
        stall_controls=stall_controls,
        journal_record_capacity=journal_record_capacity,
        film_extinction_event=True,
        film_extinction_residual_kg_per_kg_dry=residue,
        same_cell_seed=same_cell_seed,
    )


# ---------------------------------------------------------------------------
# The latch clearance and the hand-off gate
# ---------------------------------------------------------------------------
@dataclass(frozen=True, kw_only=True)
class QSCHexanePostBirthLatchClearance(_FalseStepperClaims):
    """The typed clearance record: the film is EXACTLY zero, recession may
    lawfully be demanded again, and the state belongs to the ordinary
    stall-aware march entry."""

    cleared_at_clock_s: float
    interval_index: int
    surface_temperature_k: float
    cumulative_condensed_kg: float
    cumulative_burned_off_kg: float
    cumulative_film_energy_j: float
    film_exactly_zero: ClassVar[bool] = True

    def journal_record(self) -> dict:
        return {
            "schema_id": QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_ID,
            "record": "post-birth-latch-clearance",
            "cleared_at_clock_s": self.cleared_at_clock_s,
            "interval_index": self.interval_index,
            "surface_temperature_k": self.surface_temperature_k,
            "cumulative_condensed_kg": self.cumulative_condensed_kg,
            "cumulative_burned_off_kg": self.cumulative_burned_off_kg,
            "cumulative_film_energy_j": self.cumulative_film_energy_j,
            "film_exactly_zero": True,
            "physically_qualifying": False,
        }


def require_post_birth_latch_clear(
    carrier: PostBirthFilmCarrier,
) -> QSCHexanePostBirthLatchClearance:
    """Seal the clearance when the film is EXACTLY zero; refuse otherwise."""

    if type(carrier) is not PostBirthFilmCarrier:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.CARRIER_NOT_TYPED,
            "the latch-clearance gate consumes the exact typed film carrier",
        )
    if carrier.attached_hexane_kg_per_kg_dry != 0.0:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.LATCH_NOT_CLEAR,
            "the attached film still holds "
            f"{carrier.attached_hexane_kg_per_kg_dry!r} kg/kg dry: recession "
            "stays refused until the film is EXACTLY zero; keep stepping the "
            "burn-off, never clamp",
        )
    return QSCHexanePostBirthLatchClearance(
        cleared_at_clock_s=carrier.clock_s,
        interval_index=carrier.interval_index,
        surface_temperature_k=carrier.surface_temperature_k,
        cumulative_condensed_kg=carrier.cumulative_condensed_kg,
        cumulative_burned_off_kg=carrier.cumulative_burned_off_kg,
        cumulative_film_energy_j=carrier.cumulative_film_energy_j,
    )


def demote_condensed_boundary(
    boundary: csrf.CondensedSolventReducedFilmPoreBoundary,
) -> ct.ReducedFilmPoreBoundary:
    """Return the exact plain base of one condensed-solvent boundary variant.

    The per-solve promotion discipline needs the PLAIN boundary (a fresh
    journal and a fresh channel are attached per stall solve); a committed
    birth carries the promoted variant, so the runner demotes it once here.
    Every base field is copied unchanged - the demotion IS the same film
    boundary, exactly as the promotion was.
    """

    if type(boundary) is not csrf.CondensedSolventReducedFilmPoreBoundary:
        raise QSCHexanePostBirthStepperError(
            QSCHexanePostBirthStepperRefusalCode.BOUNDARY_NOT_PLAIN,
            "demotion consumes the exact condensed-solvent boundary variant; "
            f"got {type(boundary).__name__}",
        )
    from dataclasses import fields as dataclass_fields

    return ct.ReducedFilmPoreBoundary(
        **{
            base_field.name: getattr(boundary, base_field.name)
            for base_field in dataclass_fields(ct.ReducedFilmPoreBoundary)
        }
    )


__all__ = [
    "FILM_EXTINCTION_PROBE_HALVINGS",
    "FILM_EXTINCTION_ULP_WINDOW",
    "QSC_HEXANE_POST_BIRTH_STEPPER_ORCHESTRATION_ID",
    "QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_ID",
    "QSC_HEXANE_POST_BIRTH_STEPPER_SCHEMA_REVISION",
    "PostBirthFilmCarrier",
    "QSCHexanePostBirthIntervalStep",
    "QSCHexanePostBirthLatchClearance",
    "QSCHexanePostBirthStepperError",
    "QSCHexanePostBirthStepperRefusalCode",
    "advance_post_birth_interval",
    "advance_post_birth_to_film_extinction",
    "demote_condensed_boundary",
    "initialize_post_birth_film_state",
    "probe_post_birth_film_law",
    "require_post_birth_latch_clear",
]
