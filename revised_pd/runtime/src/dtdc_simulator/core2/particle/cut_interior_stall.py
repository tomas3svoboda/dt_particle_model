r"""SA-2/A1 - the INTERIOR STALL step for a stationary primary-drainage front.

Design basis: ``docs/GT_PS2_O10B_FRONT_CHART_DESIGN_REVISION_B_2026-08-21.md``
section 2 A1 ("THE INTERIOR STALL BRANCH"), ratified 2026-08-21 for Scope A.
Every file:line claim repeated below is anchored in the read-only survey
``docs/evidence/GT_PS2_O10B_ENDPOINT_MONOTONICITY_SURVEY_2026-08-21.md``.
The construction is MIRRORED from the on-face stall chart
``cut_face_tangent`` (whose own gate ``s_dot <= 0`` ADMITS zero) with the
CC-1a rank-exchange discipline of ``cut_cap_contact_integrator``; the entry
condition is SA-1's ``cut_front_demand`` STALL classification.

WHAT THIS MODULE IS FOR.  Survey F8: a front that lawfully STOPS is refused at
six independent sites, and ``cut_geometry.FrontStepClassification.NO_MOTION``
(``cut_geometry.py:54``, ``:473-488``) exists as a classification with no step
behind it.  The corrected physics needs fronts that stop.  This module supplies
the missing step: one backward-Euler step over a caller-declared duration with
the front PINNED at its own position, bit-exactly.

THE TANGENT TRANSFER, ITEM BY ITEM (the design's own list, answered).

1. HOW THE FRONT IS PINNED.  ``cut_face_tangent`` pins the front ON the master
   face and carries no front-position unknown at all; the arrival chart does
   the same and puts the EVENT TIME in the vacated slot (survey D2.1).  Here
   the front is pinned at its own INTERIOR position: the candidate carries
   ``before.transport.geometry.front.z`` as a LITERAL - never encoded, never
   decoded, never perturbed - so ``cg.partition_master_grid`` reproduces the
   source geometry and ``cg.swept_cut_geometry`` returns an exactly zero swept
   volume.  Both facts are certified on every accepted step, not assumed.

2. THE UNKNOWN SET (rank bookkeeping).  The ordinary same-cell rank is
   ``2*nw + 3*nd + 3`` (``cut_transport.CutTransportLayout``).  Pinning the
   front removes exactly one unknown, so the stall rank is
   ``2*nw + 3*nd + 2``.  Unlike CC-1a, NOTHING replaces the vacated slot: a
   finite-duration stall step cannot carry a free front speed the way the
   zero-volume tangent does, because a nonzero ``s_dot`` in the jump rows with
   a frozen geometry would convect material no inventory change books - the
   conservation ledger would not close.  The exchange is therefore
   UNKNOWN-FOR-ROW, not unknown-for-unknown, and that is the honest difference
   from both precedents, stated rather than hidden.

3. HOW THE RH ROWS REDUCE AT ``q = 0``.  With the sweep exactly zero the three
   jump rows (``cut_transport.py:3185-3217``) lose every ``q``-term and become
   pure FLUX CONTINUITY statements:

       ``rh_water  -> A * (J_dry,w - J_wet,w)``
       ``rh_hexane -> A * J_dry,h``
       ``rh_energy -> A * (E_dry - E_wet)``

   Three conditions on two remaining interface unknowns: the stalled system is
   OVER-DETERMINED BY EXACTLY ONE ROW.  That is not a defect - it is what a
   stall IS - and the design names its resolution: "flux continuity as the
   stall CONSISTENCY residual".  One row leaves the solved set and becomes the
   certificate, exactly as ``cut_face_tangent``'s summed-component recovery
   drops ``dry_water_flux_continuity`` from its independent vector while
   keeping it in the diagnostic one (``FaceTangentLayout``).

4. WHICH ROW.  ``rh_water``, for three independent reasons and no tuning:
   (i) it is the row SA-1 already inverts into a signed front demand
   (``q = A*(J_dry - J_wet) / (C_dry - C_donor)``), so the ENTRY condition and
   the EXIT certificate are the SAME instrument at the same contract;
   (ii) the datum-covariant basis (``cut_integrator._datum_covariant_residual_
   vector``) leaves ``rh_water`` untouched while it mixes hexane into every
   energy row, so deleting ``rh_water`` - and only ``rh_water`` - keeps the
   remaining rows an exact invertible row basis and the water caloric datum an
   exact gauge freedom;
   (iii) an outward drift then surfaces through SA-1's typed refusal by
   construction rather than by a second, differently-scaled criterion.

5. ACCEPTANCE GATES.  Unchanged and non-widenable: residual ``2e-11``, ledger
   ``1e-10``, condition proxy ``1e12``, plus ``solution.success``, a zero
   ``active_mask``, strictly interior bounded primitives and a production
   reassembly with the caller's reduced-film thresholds enforced - the union of
   ``cut_face_tangent``'s and CC-1a's ladders.  The CONSISTENCY row is judged
   at the SAME ``2e-11`` against its own ``cut_integrator._residual_scales``
   row scale; no new tolerance is declared and none is relaxed.

6. HAND-OFF / SEAL.  An accepted stall returns a sealed
   :class:`InteriorStallStep` carrying one ordinary
   :class:`cut_integrator.CutIntegratorStep` (the unchanged conservation ledger
   at its ``1e-10`` contract), the SA-1 :class:`cut_front_demand.
   FrontDemandEvidence` that admitted the entry, and a second one re-classifying
   the COMMITTED state.  The committed ``after`` state is an ordinary
   accepted cut state whose front is bit-identical to the source front, so the
   ordinary recession machinery consumes it with no adapter.

WHAT DID **NOT** TRANSFER - THE FINDING THE DESIGN ASKED FOR.  The tangent's
outward-flux gate ``q_h >= 0`` (``cut_face_tangent.py:520``, ``:1026-1029``)
does NOT transfer to the interior, and it is not needed here.  On the face it
is the CONTINUATION gate: the tangent exists to hand a zero-volume seam to a
finite positive sweep, and a negative conserved hexane flux there would drive
the newborn dry cell backwards.  At an interior stall there is no newborn
volume, no continuation, and no positive sweep to protect: the hexane RH row
is a SOLVED row of the stall system and its own reduction, ``A * J_dry,h = 0``,
already pins the interface hexane flux to zero within contract.  Adding a
one-sided ``q_h >= 0`` gate on top would constrain a quantity the residual
already determines and could refuse an otherwise-certified root for a sign that
carries no interior meaning.  This is reported as a structural non-transfer,
not silently dropped; the strictness it protects on the face is carried here by
the SA-1 typed advance refusal, which is a statement about the FRONT and is
retained everywhere.

STRICTNESS AGAINST ADVANCE IS RETAINED, PHY-014 IS UNTOUCHED.  No frozen veto
is edited, no chart is widened, no tolerance is relaxed.  The front cannot move
in this step - it is a literal, not a coordinate - so the drainage vetoes
(``cut_geometry.py:418-419``, ``front.py:279-280``, ``dry_region.py:455-456``,
``sphere_coupling.py:161-162``) are all evaluated at exact equality, which each
of them already admits.  An entry demand that is ADVANCE is refused with SA-1's
:class:`cut_front_demand.PrimaryDrainageAdvanceDemanded` before any solve, and
a converged stall whose own consistency row demands an outward front is refused
the same way afterwards.

THE SIX F8 SITES ARE ROUTED AROUND, NOT WEAKENED (per-site verdicts).

* ``cut_event_orchestrator.py:770-773`` (ordering evidence ``>= 0.0``) and
  ``:305-307`` (bracket velocity) guard a FACE-ARRIVAL event.  A stalled front
  never arrives: at zero velocity the predicted event duration is not finite,
  so admitting zero there would certify an arrival that does not happen.  The
  stall step raises no arrival, builds no ordering evidence, and never enters
  the face route.  SITES STAY.
* ``accepted_history_same_cell_controller.py:174-176`` and ``:202-203`` require
  a strictly negative accepted front rate because they EXTRAPOLATE a front
  seed, ``seed_front_z = current + rate*dt``.  At zero rate that seed is the
  current front, which cannot encode in the same-cell chart at all (the chart's
  upper front bound IS the current front, survey D2.3) - the admission would
  produce an unusable seed.  The stall step constructs NO front seed.  SITES
  STAY.
* ``trajectory.py:540-543`` belongs to the legacy partial-trajectory solver
  whose unknown is a strictly positive recession INCREMENT on a one-sided logit
  (``:767-776``, survey F3's fifth encoding).  Admitting zero there is a chart
  redesign of a module this step never calls.  SITE STAYS.
* ``water_topology_oracle.py:89-90`` is a MANUFACTURED weak-solution oracle for
  the moving-front benchmark; at zero front speed its jump term vanishes and it
  verifies nothing about the stall.  This module never calls it.  SITE STAYS.
* ``cut_face_departure.py:603-606`` ("not an interior stall") is CONFIRMED to
  stay: departure into cell ``k-1`` genuinely requires strict inward recession,
  and a stalled on-face state is the tangent chart's business, not departure's.

DECLARED LIMITS, STATED RATHER THAN PAPERED OVER.

* A stall is a statement AT A DURATION.  Both the entry resolution and the
  consistency row scale with ``dt`` (the row scale carries ``inventory/dt``),
  so the same physical state is a stall over a short enough step and a resolved
  motion over a long one.  That is SA-1's declared property, inherited here
  verbatim and recorded on every ledger; it is not hidden inside a gate.
* This step SOLVES; it does not decide physics.  It cannot make a front stop -
  it represents a front the assembled balances cannot move at this duration.
* Only the DENSE bounded trust-region route is certified, for CC-1a's reason:
  this solve optimises the original datum-covariant basis with one row removed,
  not the face-limit regularized basis the structural-sparse route uses.
* No macrostep composition, no predictor, no re-march driver.  The stall-aware
  step entry here is the composition point the O9c re-march adopts.

Every claim flag on every payload in this module is ``False``.  Nothing here
makes the particle model physically qualifying, plant predictive, or production
wired.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Sequence

import numpy as np
from scipy import optimize

from dtdc_simulator.core2.particle import condensed_solvent_reduced_film as csrf
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_event_orchestrator as ceo
from dtdc_simulator.core2.particle import cut_front_demand as cfd
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut


SA2_DESIGN_BASIS = "docs/GT_PS2_O10B_FRONT_CHART_DESIGN_REVISION_B_2026-08-21.md"
SA2_SURVEY_EVIDENCE = (
    "docs/evidence/GT_PS2_O10B_ENDPOINT_MONOTONICITY_SURVEY_2026-08-21.md"
)
#: The on-face stall chart whose construction this module transfers inward.
SA2_MIRRORED_CONSTRUCTION = "src/dtdc_simulator/core2/particle/cut_face_tangent.py"
#: The pinned-coordinate rank-exchange precedent (CC-1a).
SA2_RANK_EXCHANGE_PRECEDENT = (
    "src/dtdc_simulator/core2/particle/cut_cap_contact_integrator.py"
)

#: The UNCHANGED acceptance contract.  These are ceilings, never targets: a
#: caller may declare something stricter, and a declaration looser than either
#: these values or the source state's own solver controls is refused rather
#: than adopted.
INTERIOR_STALL_RESIDUAL_TOLERANCE = 2.0e-11
INTERIOR_STALL_LEDGER_TOLERANCE = 1.0e-10
INTERIOR_STALL_MAXIMUM_CONDITION_PROXY = 1.0e12

_NONLINEAR_METHOD = "trf"
_NONLINEAR_ALGORITHM = (
    "bounded_dense_trf_finite_difference_interior_stall_pinned_front"
)
_NONLINEAR_RESIDUAL_BASIS = (
    "original_datum_covariant_minus_the_rh_water_stall_consistency_row"
)
_JACOBIAN_ROUTE_REASON = "interior_stall_pinned_front_small_rank_dense_oracle"
#: The literal written into the vacated front slot before it is DISCARDED.
#: Nothing physical depends on it: :func:`decode_interior_stall_candidate`
#: overwrites the decoded front with the pinned source literal.  Zero is used
#: because it decodes to the chart midpoint under either front policy and can
#: therefore never itself trip a decode guard.
_DISCARDED_FRONT_COORDINATE = 0.0


class InteriorStallRefusalCode(str, Enum):
    """Named refusal reasons for the interior-stall ceremony."""

    NOT_A_STALL_ENTRY = "entry_demand_is_not_a_stall"
    FRONT_NOT_PINNED = "front_was_not_pinned_bit_exactly"
    RANK_EXCHANGE_LOST = "interior_stall_rank_exchange_lost"
    TOLERANCE_WIDENED = "tolerance_widened"
    NONLINEAR_CONTRACT_FAILED = "nonlinear_residual_contract_failed"
    CONDITION_UNCERTIFIED = "jacobian_condition_proxy_uncertified"
    LEDGER_CONTRACT_FAILED = "conservation_ledger_contract_failed"
    STALL_CONSISTENCY_FAILED = "stall_consistency_row_exceeded_its_contract"
    OUTWARD_DRIFT = "stall_solve_drifted_outward"
    #: A3 (RATIFIED 2026-08-22): the film-covered admission's own typed
    #: refusal reason.  Additive member; nothing above moves.
    FILM_ADMISSION_INVALID = "film_covered_admission_invalid"


class InteriorStallError(RuntimeError):
    """Rejected interior stall; base of the typed refusal family.

    ``rollback_state`` is the exact immutable source state.  A caller that
    catches this has lost nothing: no field of the source was mutated, no
    candidate was adopted, and no tolerance was relaxed.
    """

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState | None = None,
        *,
        refusal_code: InteriorStallRefusalCode | None = None,
        nonlinear_evaluations: int = 0,
        rejected_trial_evaluations: int = 0,
        nonlinear_optimizer_function_evaluations: int = 0,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        stall_consistency_residual: float = math.inf,
        last_candidate: cut.CutTransportUnknowns | None = None,
        scaled_residuals: tuple[float, ...] = (),
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.refusal_code = refusal_code
        self.nonlinear_evaluations = nonlinear_evaluations
        self.rejected_trial_evaluations = rejected_trial_evaluations
        self.nonlinear_optimizer_function_evaluations = (
            nonlinear_optimizer_function_evaluations
        )
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.stall_consistency_residual = stall_consistency_residual
        self.last_candidate = last_candidate
        self.scaled_residuals = scaled_residuals


class InteriorStallEntryError(InteriorStallError):
    """The state handed to the stall step is not a stall at this duration.

    A ``RECESSION_DEMAND`` entry belongs to the ordinary same-cell solve
    (design section 2 A3, the only branch that still moves a front); an
    ``ADVANCE_DEMAND`` entry never reaches this error - it is refused earlier
    with SA-1's :class:`cut_front_demand.PrimaryDrainageAdvanceDemanded`, which
    names the frozen boundary.
    """

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState | None = None,
        *,
        evidence: cfd.FrontDemandEvidence | None = None,
        refusal_code: InteriorStallRefusalCode = (
            InteriorStallRefusalCode.NOT_A_STALL_ENTRY
        ),
    ) -> None:
        super().__init__(message, rollback_state, refusal_code=refusal_code)
        self.evidence = evidence


class InteriorStallConsistencyError(InteriorStallError):
    """The converged step is not consistent with a stationary front.

    The stall consistency residual is the ``q = 0`` reduction of the
    Rankine-Hugoniot WATER row - the one row the rank exchange removed from the
    solved set.  Exceeding its contract means the assembled interface balance
    demands a front motion the step did not perform, so the step is refused
    with exact rollback instead of being committed as if it had stalled.
    """


@dataclass(frozen=True)
class InteriorStallLayout:
    """The rank of one pinned-front stall, in the tangent's own shape.

    ``independent_residual_blocks`` is what the solver drives to zero;
    ``diagnostic_residual_blocks`` is every row the assembly reports, including
    the ``rh_water`` consistency row.  This mirrors
    :class:`cut_face_tangent.FaceTangentLayout`, whose dependent
    ``dry_water_flux_continuity`` row lives in exactly the same relationship to
    its own independent vector.
    """

    wet_piece_count: int
    dry_piece_count: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        for label, value in (
            ("wet", self.wet_piece_count),
            ("dry", self.dry_piece_count),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"an interior stall needs a positive {label} rank")

    @property
    def dry_face_count(self) -> int:
        return self.dry_piece_count + 1

    @property
    def unknown_blocks(self) -> tuple[tuple[str, int], ...]:
        """Return the ordinary blocks with ``front_z`` REMOVED (pinned)."""

        return (
            ("wet_temperature", self.wet_piece_count),
            ("wet_retained_water", self.wet_piece_count),
            ("dry_temperature", self.dry_piece_count),
            ("dry_y_hexane", self.dry_piece_count),
            ("dry_total_stefan_flux", self.dry_face_count),
            ("interface_temperature", 1),
        )

    @property
    def independent_residual_blocks(self) -> tuple[tuple[str, int], ...]:
        """Return the solved rows: the ordinary set minus ``rh_water``."""

        return (
            ("wet_water", self.wet_piece_count),
            ("wet_energy", self.wet_piece_count),
            ("dry_water", self.dry_piece_count),
            ("dry_hexane", self.dry_piece_count),
            ("dry_energy", self.dry_piece_count),
            ("rh_hexane", 1),
            ("rh_energy", 1),
        )

    @property
    def diagnostic_residual_blocks(self) -> tuple[tuple[str, int], ...]:
        """Return every reported row, consistency row included."""

        return (
            ("wet_water", self.wet_piece_count),
            ("wet_energy", self.wet_piece_count),
            ("dry_water", self.dry_piece_count),
            ("dry_hexane", self.dry_piece_count),
            ("dry_energy", self.dry_piece_count),
            ("rh_water", 1),
            ("rh_hexane", 1),
            ("rh_energy", 1),
        )

    @property
    def unknown_count(self) -> int:
        return sum(size for _, size in self.unknown_blocks)

    @property
    def independent_residual_count(self) -> int:
        return sum(size for _, size in self.independent_residual_blocks)

    @property
    def diagnostic_residual_count(self) -> int:
        return sum(size for _, size in self.diagnostic_residual_blocks)

    @property
    def ordinary_rank(self) -> int:
        return 2 * self.wet_piece_count + 3 * self.dry_piece_count + 3

    @property
    def rank_formula(self) -> int:
        """Return the documented stall rank ``2*nw + 3*nd + 2``."""

        return 2 * self.wet_piece_count + 3 * self.dry_piece_count + 2

    @property
    def is_square(self) -> bool:
        return (
            self.unknown_count
            == self.independent_residual_count
            == self.rank_formula
        )


@dataclass(frozen=True)
class InteriorStallControls:
    """The UNCHANGED acceptance contract for one interior-stall step.

    Mirrors :class:`cut_cap_contact_integrator.CapContactControls`: the three
    numerical contracts are validated against their frozen ceilings here and
    against the source state's own controls in :func:`solve_interior_stall`.
    Neither check can ever widen anything.
    """

    nonlinear_residual_tolerance: float = INTERIOR_STALL_RESIDUAL_TOLERANCE
    ledger_tolerance: float = INTERIOR_STALL_LEDGER_TOLERANCE
    maximum_condition_proxy: float = INTERIOR_STALL_MAXIMUM_CONDITION_PROXY
    maximum_function_evaluations: int = 800
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        for name, value, ceiling in (
            (
                "residual",
                self.nonlinear_residual_tolerance,
                INTERIOR_STALL_RESIDUAL_TOLERANCE,
            ),
            ("ledger", self.ledger_tolerance, INTERIOR_STALL_LEDGER_TOLERANCE),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"interior-stall {name} tolerance must be positive")
            if value > ceiling:
                raise ValueError(
                    f"interior-stall {name} tolerance cannot exceed its frozen "
                    f"{ceiling:.1e} contract"
                )
        if (
            not math.isfinite(self.maximum_condition_proxy)
            or self.maximum_condition_proxy <= 1.0
            or self.maximum_condition_proxy
            > INTERIOR_STALL_MAXIMUM_CONDITION_PROXY
        ):
            raise ValueError(
                "interior-stall condition-proxy ceiling must lie in "
                f"(1, {INTERIOR_STALL_MAXIMUM_CONDITION_PROXY:.1e}]"
            )
        try:
            budget = operator.index(self.maximum_function_evaluations)
        except TypeError as exc:
            raise ValueError(
                "interior-stall evaluation budget must be an integer"
            ) from exc
        if isinstance(self.maximum_function_evaluations, bool) or budget < 1:
            raise ValueError("interior-stall evaluation budget must be positive")


@dataclass(frozen=True)
class InteriorStallSeed:
    """One explicit seed for the non-front unknowns plus the Stefan scales.

    The seed candidate is an ORDINARY same-cell candidate.  Its front
    coordinate is IGNORED and replaced by the source front literal - this
    module owns the pin - so a seed carrying any other front is refused rather
    than quietly corrected.
    """

    candidate: cut.CutTransportUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, cut.CutTransportUnknowns):
            raise TypeError("interior-stall seed requires CutTransportUnknowns")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0
            for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every interior-stall Stefan face needs a positive scale")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("interior-stall seed needs a provenance label")


@dataclass(frozen=True)
class FilmCoveredStallAdmission:
    """A3 (RATIFIED 2026-08-22): the typed film-covered stall declaration.

    Authority: ``docs/GT_PS2_BC4_COMMITTED_BIRTH_AND_BURNOFF_SEAM_2026-08-22.md``
    section 3 A3, RULED block.  BC-4b GAP 3 measured that no operator could
    lawfully advance the partial-state interior while a liquid film covers the
    surface: this step's accepted reassembly owned the GAS-ONLY reduced-film
    surface, which the mobile-liquid branch forbids under an attached film.

    Under this admission - DEFAULT-OFF, the ``None`` path is byte-identical to
    today's - the pinned-front stall runs on the condensed-solvent boundary
    variant (BC-3's own dispatch at the D4 seam): sweep samples strictly above
    the hexane dew locus evaluate byte-identically to the legacy closure, and
    a sample or root at/below the locus takes the exact-chart ``y_dew``
    construction with typed evidence.  A converged surface root on the
    condensed branch is then the LAWFUL film-covered surface - recorded
    through the boundary's fresh single-admission channel - instead of the
    pre-A3 typed refusal.

    THE LATCH IS THE ENTRY CONDITION, NOT A CASUALTY.  The stall pins the
    front bit-exactly, which is precisely what the drainage latch requires
    while ``film > 0``: recession stays refused until the film is EXACTLY
    zero, and this admission never relaxes that - it only lets the interior
    heat through a surface the film lawfully covers.  The companion carries
    the film's identity: it must be bound (``is_bound_to``) to the exact
    source state this stall advances, and the declared film must be strictly
    positive - a vanished film means the latch has cleared and the ordinary
    gas-only stall owns the state.
    """

    film_companion: object
    attached_hexane_kg_per_kg_dry: float
    film_covered_state_is_pre_recession: bool = field(default=True, init=False)
    recession_stays_refused_until_film_exactly_zero: bool = field(
        default=True, init=False
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.film_companion is None:
            raise InteriorStallError(
                "the film-covered stall admission needs the committed film "
                "companion; an absent companion is the gas-only stall's case",
                refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
            )
        if not callable(getattr(self.film_companion, "is_bound_to", None)):
            raise InteriorStallError(
                "the film companion carries no is_bound_to identity predicate; "
                "an unbindable companion is no companion",
                refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
            )
        value = self.attached_hexane_kg_per_kg_dry
        if type(value) is not float or not math.isfinite(value) or value <= 0.0:
            raise InteriorStallError(
                "the film-covered stall admission requires a strictly positive "
                f"finite attached film, got {value!r}: at EXACTLY zero the "
                "drainage latch clears and the ordinary gas-only stall owns "
                "the state",
                refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
            )


@dataclass(frozen=True)
class InteriorStallRankAudit:
    """Structural proof of the unknown-for-row exchange at a pinned front."""

    ordinary_unknown_count: int
    ordinary_residual_count: int
    stall_unknown_count: int
    independent_residual_count: int
    rank_formula: int
    wet_piece_count: int
    dry_piece_count: int
    pinned_front_coordinate_slot: int
    stall_consistency_row_index: int
    front_unknown_count: int = field(default=0, init=False)
    replacement_unknown_count: int = field(default=0, init=False)
    removed_residual_row_count: int = field(default=1, init=False)
    new_phase_or_source_was_added: bool = field(default=False, init=False)
    uses_unchanged_backward_euler_rows: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.ordinary_unknown_count != self.ordinary_residual_count:
            raise ValueError("the underlying same-cell residual must stay square")
        if self.stall_unknown_count + 1 != self.ordinary_unknown_count:
            raise ValueError("pinning the front must remove exactly one unknown")
        if self.independent_residual_count + 1 != self.ordinary_residual_count:
            raise ValueError("the stall must remove exactly one residual row")
        if self.stall_unknown_count != self.independent_residual_count:
            raise ValueError("the pinned-front stall system is not square")
        if self.rank_formula != self.stall_unknown_count:
            raise ValueError(
                "interior-stall rank must remain the documented 2*nw + 3*nd + 2"
            )
        if self.pinned_front_coordinate_slot != self.ordinary_unknown_count - 2:
            raise ValueError(
                "the pinned slot is not the ordinary chart's front_z coordinate"
            )
        if self.stall_consistency_row_index != (
            2 * self.wet_piece_count + 3 * self.dry_piece_count
        ):
            raise ValueError(
                "the stall consistency row is not the datum-covariant rh_water row"
            )

    @property
    def is_square(self) -> bool:
        return self.stall_unknown_count == self.independent_residual_count


@dataclass(frozen=True)
class InteriorStallLedger:
    """Pin, rank, nonlinear, consistency, and conservation certificate."""

    dt_s: float
    pinned_front_z: float
    pinned_front_radius_m: float
    committed_front_z: float
    front_pinned_bit_exactly: bool
    swept_wet_to_dry_volume_m3: float
    interface_swept_volume_rate_m3_s: float
    wet_hexane_sweep_identity_mol_s: float
    layout: InteriorStallLayout
    rank_audit: InteriorStallRankAudit
    #: The ORDINARY same-cell chart's front interval.  Its upper endpoint IS
    #: the current front (``cut_integrator.py:1624``), which is exactly why the
    #: pinned front cannot be a solver coordinate here.
    same_cell_chart_front_z: tuple[float, float]
    pinned_front_is_on_the_same_cell_chart_bound: bool
    #: The host CELL's own mesh-geometric span, used ONLY for the accepted-step
    #: bounded-primitive diagnostic.  Nothing is ever encoded into it - the
    #: stall coordinate vector has no front slot at all - so this is not the D1
    #: widened solve chart Revision B declined (design section 2 A3).
    front_diagnostic_interval_z: tuple[float, float]
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    nonlinear_optimizer_function_evaluations: int
    nonlinear_message: str
    nonlinear_algorithm: str
    nonlinear_residual_basis: str
    jacobian_route_reason: str
    used_sparse_jacobian: bool
    independent_scaled_residuals: tuple[float, ...]
    maximum_independent_scaled_residual: float
    diagnostic_scaled_residuals: tuple[float, ...]
    stall_consistency_row_mol_s: float
    stall_consistency_row_scale_mol_s: float
    normalized_stall_consistency_row: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    declared_nonlinear_residual_tolerance: float
    declared_ledger_tolerance: float
    declared_maximum_condition_proxy: float
    maximum_step_ledger_residual: float
    maximum_cumulative_ledger_residual: float
    entry_demanded_front_volume_rate_m3_s: float
    entry_resolution_m3_s: float
    #: The consistency row of the CONVERGED assembly inverted through its own
    #: jump; the exit fields below re-read the same row on the committed state,
    #: where the O9a donor sits at the new time level.
    converged_consistency_demanded_rate_m3_s: float
    exit_demanded_front_volume_rate_m3_s: float
    exit_resolution_m3_s: float
    exit_outcome: cfd.FrontDemandOutcome
    accepted: bool
    front_encoded_as_a_solver_coordinate: bool = field(default=False, init=False)
    front_z_chart_widened: bool = field(default=False, init=False)
    tolerance_widened: bool = field(default=False, init=False)
    componentwise_clipping_used: bool = field(default=False, init=False)
    outward_front_admitted: bool = field(default=False, init=False)
    #: Reported, not hidden: both the entry resolution and this consistency row
    #: carry the assembly's own ``inventory/dt`` row scale, so a stall is a
    #: statement AT A DURATION.  SA-1 declares the same property.
    stall_verdict_depends_on_duration: bool = field(default=True, init=False)
    macrostep_ledger_composed: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square or not self.rank_audit.is_square:
            raise ValueError("an interior-stall ledger must record a square rank")
        if not self.front_pinned_bit_exactly:
            raise InteriorStallConsistencyError(
                "an interior-stall ledger cannot record a front that moved",
                refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            )
        if self.committed_front_z != self.pinned_front_z:
            raise InteriorStallConsistencyError(
                "the committed front is not the pinned source front",
                refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            )
        if self.swept_wet_to_dry_volume_m3 != 0.0 or (
            self.interface_swept_volume_rate_m3_s != 0.0
        ):
            raise InteriorStallConsistencyError(
                "an interior stall must sweep exactly zero volume; got "
                f"{self.swept_wet_to_dry_volume_m3!r} m3 at "
                f"{self.interface_swept_volume_rate_m3_s!r} m3/s",
                refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            )
        if self.wet_hexane_sweep_identity_mol_s != 0.0:
            raise InteriorStallConsistencyError(
                "a pinned front consumes no wet n-hexane, so the sweep identity "
                f"must be exactly zero; got {self.wet_hexane_sweep_identity_mol_s!r}",
                refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            )
        if self.exit_outcome is cfd.FrontDemandOutcome.ADVANCE_DEMAND:
            raise InteriorStallConsistencyError(
                "an interior-stall ledger may never record an outward exit demand",
                refusal_code=InteriorStallRefusalCode.OUTWARD_DRIFT,
            )
        if len(self.independent_scaled_residuals) != (
            self.layout.independent_residual_count
        ):
            raise ValueError("the independent residual record lost its rank")
        if len(self.diagnostic_scaled_residuals) != (
            self.layout.diagnostic_residual_count
        ):
            raise ValueError("the diagnostic residual record lost a row")


@dataclass(frozen=True)
class InteriorStallStep:
    """Sealed committed interior-stall step: the front stayed, the clock moved.

    ``step`` is an ordinary :class:`cut_integrator.CutIntegratorStep`; its
    ``after`` state is an accepted cut state at ``time + dt`` whose front is
    bit-identical to the source front and whose unchanged conservation ledger
    closes at the declared ``1e-10`` contract.  ``entry_evidence`` is the SA-1
    :class:`cut_front_demand.FrontDemandEvidence` that admitted the entry.
    """

    before: ci.CutIntegratorState
    boundary: ct.PoreBoundary
    seed: InteriorStallSeed
    candidate: cut.CutTransportUnknowns
    step: ci.CutIntegratorStep
    entry_evidence: cfd.FrontDemandEvidence
    exit_evidence: cfd.FrontDemandEvidence
    ledger: InteriorStallLedger
    front_pinned_not_charted: bool = field(default=True, init=False)
    entered_after_same_cell_refusal: bool = False
    #: A3 (RATIFIED 2026-08-22): the typed film-covered declaration this stall
    #: ran under, or ``None`` for the gas-only stall (the byte-identical
    #: legacy path).
    film_covered_admission: FilmCoveredStallAdmission | None = None
    #: The condensed-root admission the boundary's channel recorded when the
    #: accepted surface root took the condensed branch, or ``None`` when the
    #: root evaluated strictly above the locus (byte-identical legacy surface).
    #: Never non-``None`` without the film-covered admission.
    condensed_root_admission: csrf.CondensedRootCommitAdmission | None = None
    macrostep_ledger_composed: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)
    plant_predictive: bool = field(default=False, init=False)
    production_wired: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.step, ci.CutIntegratorStep):
            raise TypeError("an interior stall needs a committed same-cell step")
        if self.step.before is not self.before:
            raise ValueError("the interior stall lost its exact source identity")
        if self.step.boundary is not self.boundary:
            raise ValueError("the interior stall changed the current pore boundary")
        if self.step.ledger.dt_s != self.ledger.dt_s:
            raise ValueError("the committed step duration is not the stall duration")
        for label, payload in (
            ("entry", self.entry_evidence),
            ("exit", self.exit_evidence),
        ):
            if not isinstance(payload, cfd.FrontDemandEvidence):
                raise TypeError(
                    f"an interior stall needs typed SA-1 {label} evidence"
                )
            if payload.outcome is not cfd.FrontDemandOutcome.STALL:
                raise ValueError(
                    "an interior stall may only seal a STALL "
                    f"{label} classification; got {payload.outcome.value}"
                )
        if self.after.transport.geometry.front.z != (
            self.before.transport.geometry.front.z
        ):
            raise InteriorStallConsistencyError(
                "the committed stall state moved its front",
                self.before,
                refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            )
        if self.after.transport.geometry.front.radius_m != (
            self.before.transport.geometry.front.radius_m
        ):
            raise InteriorStallConsistencyError(
                "the committed stall state moved its front radius",
                self.before,
                refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            )
        if not isinstance(self.entered_after_same_cell_refusal, bool):
            raise TypeError("the stall entry provenance flag must be boolean")
        if self.film_covered_admission is not None and (
            type(self.film_covered_admission) is not FilmCoveredStallAdmission
        ):
            raise TypeError(
                "a film-covered stall seals the exact FilmCoveredStallAdmission"
            )
        if self.condensed_root_admission is not None:
            if self.film_covered_admission is None:
                raise ValueError(
                    "a condensed surface root may seal only under the typed "
                    "film-covered admission (A3); a gas-only stall has no "
                    "lawful condensed root"
                )
            if type(self.condensed_root_admission) is not (
                csrf.CondensedRootCommitAdmission
            ):
                raise TypeError(
                    "the condensed-root evidence must be the exact typed "
                    "CondensedRootCommitAdmission"
                )

    @property
    def after(self) -> ci.CutIntegratorState:
        return self.step.after

    @property
    def dt_s(self) -> float:
        return self.ledger.dt_s

    @property
    def front_z(self) -> float:
        return self.ledger.pinned_front_z


class StallAwareBranch(str, Enum):
    """Which lawful branch of design section 2 one step actually took."""

    #: The front was pinned and the step committed (design section 2 A1).
    INTERIOR_STALL = "interior_stall"
    #: The ordinary moving-front solve ran (design section 2 A3).
    ORDINARY_RECESSION = "ordinary_recession"


@dataclass(frozen=True)
class StallAwareStep:
    """One step taken under the typed three-way front verdict.

    This is the composition the O9c re-march adopts: exactly one of
    ``stall_step`` / ``macrostep`` is present, the SA-1 verdict that selected it
    is carried, and the forbidden third outcome never returns a payload at all
    - it raises :class:`cut_front_demand.PrimaryDrainageAdvanceDemanded`.
    """

    before: ci.CutIntegratorState
    boundary: ct.PoreBoundary
    dt_s: float
    entry_evidence: cfd.FrontDemandEvidence
    branch: StallAwareBranch
    stall_step: InteriorStallStep | None
    macrostep: ceo.CutEventMacrostep | None
    physically_qualifying: bool = field(default=False, init=False)
    plant_predictive: bool = field(default=False, init=False)
    production_wired: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if type(self.branch) is not StallAwareBranch:
            raise TypeError("the stall-aware branch must use the typed enum")
        if (self.stall_step is None) == (self.macrostep is None):
            raise ValueError("a stall-aware step carries exactly one committed payload")
        if self.branch is StallAwareBranch.INTERIOR_STALL:
            if not isinstance(self.stall_step, InteriorStallStep):
                raise TypeError("the stall branch needs a sealed InteriorStallStep")
            if self.entry_evidence.outcome is not cfd.FrontDemandOutcome.STALL:
                raise ValueError("the stall branch requires a STALL verdict")
        elif not isinstance(self.macrostep, ceo.CutEventMacrostep):
            raise TypeError("the recession branch needs a committed macrostep")
        elif self.entry_evidence.outcome is not (
            cfd.FrontDemandOutcome.RECESSION_DEMAND
        ):
            raise ValueError("the recession branch requires a RECESSION verdict")

    @property
    def after(self) -> ci.CutIntegratorState:
        if self.stall_step is not None:
            return self.stall_step.after
        assert self.macrostep is not None  # narrowed by __post_init__
        return self.macrostep.after


# ---------------------------------------------------------------------------
# Rank bookkeeping
# ---------------------------------------------------------------------------
def interior_stall_layout(layout: cut.CutTransportLayout) -> InteriorStallLayout:
    """Return the pinned-front rank of one ordinary same-cell layout."""

    if not isinstance(layout, cut.CutTransportLayout):
        raise TypeError("the interior-stall layout needs a CutTransportLayout")
    return InteriorStallLayout(
        wet_piece_count=layout.wet_piece_count,
        dry_piece_count=layout.dry_piece_count,
    )


def pinned_front_coordinate_slot(layout: cut.CutTransportLayout) -> int:
    """Return the coordinate slot the pinned front vacates.

    ``front_z`` is the second-to-last ordinary unknown block
    (``cut_transport.CutTransportLayout.unknown_blocks``), which is also the
    slot ``cut_integrator._coordinate_solver_bounds_for_chart`` addresses as
    ``[-2]``.  The slot is DELETED from the stall coordinate vector; it is not
    re-used.
    """

    if not isinstance(layout, cut.CutTransportLayout):
        raise TypeError("the pinned front slot needs a CutTransportLayout")
    blocks = layout.unknown_blocks
    if blocks[-2][0] != "front_z":
        raise RuntimeError(
            "the ordinary same-cell chart no longer carries front_z where the "
            "interior stall pins it"
        )
    return layout.unknown_count - 2


def stall_consistency_row_index(layout: cut.CutTransportLayout) -> int:
    """Return the ``rh_water`` row index in the datum-covariant basis.

    ``cut_integrator._datum_covariant_residual_vector`` preserves the block
    ordering and leaves ``rh_water`` UNTRANSFORMED, so this index addresses the
    same physical row in both bases.
    """

    if not isinstance(layout, cut.CutTransportLayout):
        raise TypeError("the stall consistency row needs a CutTransportLayout")
    blocks = layout.residual_blocks
    if blocks[-3][0] != "rh_water":
        raise RuntimeError(
            "the ordinary same-cell residual no longer carries rh_water where "
            "the interior stall reads its consistency row"
        )
    return 2 * layout.wet_piece_count + 3 * layout.dry_piece_count


def audit_interior_stall_rank(
    before: ci.CutIntegratorState,
) -> InteriorStallRankAudit:
    """Return the exact unknown-for-row exchange at one pinned front."""

    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("the interior-stall rank audit needs a cut integrator state")
    layout = before.transport.layout
    stall = interior_stall_layout(layout)
    return InteriorStallRankAudit(
        ordinary_unknown_count=layout.unknown_count,
        ordinary_residual_count=layout.residual_count,
        stall_unknown_count=stall.unknown_count,
        independent_residual_count=stall.independent_residual_count,
        rank_formula=stall.rank_formula,
        wet_piece_count=layout.wet_piece_count,
        dry_piece_count=layout.dry_piece_count,
        pinned_front_coordinate_slot=pinned_front_coordinate_slot(layout),
        stall_consistency_row_index=stall_consistency_row_index(layout),
    )


# ---------------------------------------------------------------------------
# The pin, the seed, and the coordinate map
# ---------------------------------------------------------------------------
def stall_seed_from_state(
    state: ci.CutIntegratorState,
    *,
    stefan_flux_scale_floor_mol_m2_s: float = 1.0e-3,
    label: str = "accepted-state primitives at a pinned front",
) -> InteriorStallSeed:
    """Return the accepted-state seed for one stall, front already pinned.

    The primitives are the state's own - the accepted-history controller's
    discipline (``accepted_history_same_cell_controller.py:191-222``) minus its
    front extrapolation, which a stall has nothing to extrapolate.  No rejected
    candidate, clipping, or field tangent is involved.
    """

    if not isinstance(state, ci.CutIntegratorState):
        raise TypeError("the interior-stall seed needs a cut integrator state")
    if (
        not math.isfinite(stefan_flux_scale_floor_mol_m2_s)
        or stefan_flux_scale_floor_mol_m2_s <= 0.0
    ):
        raise ValueError("the Stefan scale floor must be positive and finite")
    candidate = cut.candidate_from_state(
        state.transport,
        dry_total_stefan_fluxes_mol_m2_s=state.last_total_stefan_fluxes_mol_m2_s,
        interface_temperature_k=state.last_interface_temperature_k,
    )
    scales = tuple(
        max(abs(value), stefan_flux_scale_floor_mol_m2_s)
        for value in state.last_total_stefan_fluxes_mol_m2_s
    )
    return InteriorStallSeed(candidate, scales, label)


def certify_front_is_pinned(
    before: ci.CutIntegratorState,
    candidate: cut.CutTransportUnknowns,
) -> None:
    """Refuse any candidate whose front is not the source front, bit for bit.

    Exact comparison, never a tolerance: the pin is a LITERAL copy of
    ``before.transport.geometry.front.z``, so anything else is a different
    front and a different chart.
    """

    if not isinstance(candidate, cut.CutTransportUnknowns):
        raise TypeError("the pin certification needs CutTransportUnknowns")
    pinned = before.transport.geometry.front.z
    if candidate.front_z != pinned:
        raise InteriorStallError(
            "an interior-stall candidate must carry the source front exactly; "
            f"got {candidate.front_z!r} against {pinned!r}",
            before,
            refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
            last_candidate=candidate,
        )


def encode_interior_stall_coordinates(
    candidate: cut.CutTransportUnknowns,
    chart: "ci._CoordinateChart",
    stefan_flux_scales_mol_m2_s: Sequence[float],
    *,
    layout: cut.CutTransportLayout,
) -> np.ndarray:
    """Encode the non-front coordinates; the front slot is DELETED.

    The candidate's own front cannot be encoded at all - the ordinary chart's
    upper front bound IS the current front (``cut_integrator.py:1624``) and
    ``_encode_open`` refuses an endpoint.  The encode therefore runs on a
    throwaway front at the chart midpoint and the slot is removed immediately;
    that value reaches no residual, no bound, and no committed state.
    """

    slot = pinned_front_coordinate_slot(layout)
    midpoint = 0.5 * math.fsum(chart.front_z)
    encodable = replace(candidate, front_z=midpoint)
    full = np.asarray(
        ci._encode_candidate_for_chart(  # noqa: SLF001
            encodable,
            chart,
            stefan_flux_scales_mol_m2_s,
        ),
        dtype=float,
    )
    coordinates = np.delete(full, slot)
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("encoded interior-stall coordinates must be finite")
    return coordinates


def decode_interior_stall_candidate(
    coordinates: Sequence[float],
    layout: cut.CutTransportLayout,
    chart: "ci._CoordinateChart",
    stefan_flux_scales_mol_m2_s: Sequence[float],
    *,
    pinned_front_z: float,
) -> cut.CutTransportUnknowns:
    """Decode the reduced vector and write the pinned front back in exactly."""

    slot = pinned_front_coordinate_slot(layout)
    values = np.asarray(tuple(float(value) for value in coordinates), dtype=float)
    if values.shape != (layout.unknown_count - 1,) or not np.all(np.isfinite(values)):
        raise ValueError(
            "interior-stall coordinate vector is non-finite or has the wrong rank"
        )
    restored = np.insert(values, slot, _DISCARDED_FRONT_COORDINATE)
    candidate = ci._decode_candidate_for_chart(  # noqa: SLF001
        restored,
        layout,
        chart,
        stefan_flux_scales_mol_m2_s,
    )
    pinned = replace(candidate, front_z=pinned_front_z)
    if pinned.front_z != pinned_front_z:
        raise RuntimeError("the interior-stall pin did not survive its own decode")
    return pinned


# ---------------------------------------------------------------------------
# The step
# ---------------------------------------------------------------------------
def solve_interior_stall(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    seed: InteriorStallSeed | None = None,
    *,
    controls: InteriorStallControls | None = None,
    entry_evidence: cfd.FrontDemandEvidence | None = None,
    entered_after_same_cell_refusal: bool = False,
    film_covered_admission: FilmCoveredStallAdmission | None = None,
) -> InteriorStallStep:
    """Commit one interior-stall step or refuse with exact rollback.

    The front is pinned at ``before``'s own position, the vacated unknown is
    removed from the coordinate vector, the ``rh_water`` jump row is removed
    from the solved set and becomes the STALL CONSISTENCY residual, and every
    remaining row is the unchanged
    :func:`cut_transport.assemble_backward_euler` assembly.

    ``entry_evidence`` is SA-1's verdict at ``(before, dt_s)``; when omitted it
    is evaluated here.  A ``STALL`` verdict is the entry condition, an
    ``ADVANCE_DEMAND`` verdict raises
    :class:`cut_front_demand.PrimaryDrainageAdvanceDemanded`, and a
    ``RECESSION_DEMAND`` verdict raises :class:`InteriorStallEntryError`
    pointing at the ordinary solve.

    ``film_covered_admission`` (A3, RATIFIED 2026-08-22; DEFAULT-OFF) declares
    that this stall advances the strict-partial state UNDER a committed
    attached film: the boundary must then be the exact condensed-solvent
    reduced-film variant carrying a FRESH single-admission channel, sweep
    samples strictly above the hexane dew locus evaluate byte-identically to
    the legacy closure, and an accepted surface root at/below the locus is the
    lawful film-covered surface, recorded through the channel and sealed on
    the step.  With ``None`` - the default - every path below is
    byte-identical to the pre-A3 step, and a channel-armed condensed boundary
    is refused typed so the new lawful path can never be taken silently.
    """

    evaluations = 0
    rejected_trials = 0
    optimizer_function_evaluations = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    consistency = math.inf
    last_candidate: cut.CutTransportUnknowns | None = None
    last_scaled: tuple[float, ...] = ()
    try:
        if not isinstance(before, ci.CutIntegratorState):
            raise TypeError("an interior stall needs a cut integrator state")
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("an interior stall requires a supported pore boundary")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("the interior-stall duration must be positive and finite")
        declared = InteriorStallControls() if controls is None else controls
        if not isinstance(declared, InteriorStallControls):
            raise TypeError("an interior stall needs typed InteriorStallControls")
        _certify_declared_contract_is_not_wider(declared, before.controls)
        if not isinstance(entered_after_same_cell_refusal, bool):
            raise TypeError("the stall entry provenance flag must be boolean")

        # A3 (RATIFIED 2026-08-22): the film-covered admission's typed gates.
        # Everything here fires BEFORE any solve, with exact rollback.
        if film_covered_admission is not None:
            if type(film_covered_admission) is not FilmCoveredStallAdmission:
                raise InteriorStallError(
                    "the film-covered stall consumes the exact typed "
                    "FilmCoveredStallAdmission",
                    before,
                    refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
                )
            if film_covered_admission.film_companion.is_bound_to(before) is not True:
                raise InteriorStallError(
                    "the film companion is not identity-bound to the exact "
                    "source state this stall advances; a film record can "
                    "neither be re-attached to another state nor survive a "
                    "state it does not describe",
                    before,
                    refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
                )
            if type(boundary) is not csrf.CondensedSolventReducedFilmPoreBoundary:
                raise InteriorStallError(
                    "the film-covered stall requires the exact condensed-"
                    "solvent reduced-film boundary variant (BC-3's D4-seam "
                    f"dispatch); got {type(boundary).__name__}: a gas-only "
                    "boundary cannot lawfully evaluate the film-covered "
                    "surface",
                    before,
                    refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
                )
            film_channel = boundary.condensed_root_commit_channel
            if film_channel is None:
                raise InteriorStallError(
                    "the film-covered stall requires an armed condensed-root "
                    "commit channel on its boundary: the accepted condensed "
                    "surface root is recorded there as the step's typed "
                    "evidence; construct a fresh channel per solve",
                    before,
                    refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
                )
            if film_channel.admitted:
                raise InteriorStallError(
                    "the boundary's condensed-root commit channel already "
                    "carries an admission: one channel books exactly one "
                    "accepted root, so construct a fresh channel per solve, "
                    "never reuse one across solves",
                    before,
                    refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
                )
        elif (
            type(boundary) is csrf.CondensedSolventReducedFilmPoreBoundary
            and boundary.condensed_root_commit_channel is not None
        ):
            raise InteriorStallError(
                "a condensed-root commit channel on the stall boundary "
                "requires the typed film-covered admission (A3 is DEFAULT-"
                "OFF): without the declaration the lawful film-covered path "
                "may never be taken silently",
                before,
                refusal_code=InteriorStallRefusalCode.FILM_ADMISSION_INVALID,
            )

        layout = before.transport.layout
        if layout.cut_cell_index is None:
            raise InteriorStallError(
                "an interior stall requires a strict partial-front cut chart",
                before,
                refusal_code=InteriorStallRefusalCode.RANK_EXCHANGE_LOST,
            )
        rank = audit_interior_stall_rank(before)
        stall_layout = interior_stall_layout(layout)
        if not stall_layout.is_square or not rank.is_square:
            raise InteriorStallError(
                "the interior-stall rank exchange is not square",
                before,
                refusal_code=InteriorStallRefusalCode.RANK_EXCHANGE_LOST,
            )
        if layout.unknown_count > before.controls.maximum_dense_unknowns:
            raise InteriorStallError(
                "SA-2 certifies only the dense bounded trust-region stall route; "
                f"unknown rank {layout.unknown_count} exceeds the declared dense "
                f"ceiling {before.controls.maximum_dense_unknowns}",
                before,
                refusal_code=InteriorStallRefusalCode.RANK_EXCHANGE_LOST,
            )

        # ENTRY: SA-1's typed verdict, at this state and this duration.
        evidence = (
            cfd.evaluate_front_demand(
                before,
                boundary,
                dt_s,
                residual_tolerance=declared.nonlinear_residual_tolerance,
            )
            if entry_evidence is None
            else entry_evidence
        )
        if not isinstance(evidence, cfd.FrontDemandEvidence):
            raise TypeError("an interior stall needs typed SA-1 entry evidence")
        if entry_evidence is not None and evidence.probe_duration_s != dt_s:
            raise InteriorStallEntryError(
                "the supplied SA-1 evidence was judged over a different duration: "
                f"{evidence.probe_duration_s!r} against {dt_s!r}",
                before,
                evidence=evidence,
            )
        cfd.refuse_on_advance_demand(
            evidence,
            before,
            surfaced_at=(
                "cut_interior_stall.solve_interior_stall (entry; design 2 A1/A2)"
            ),
        )
        if evidence.outcome is not cfd.FrontDemandOutcome.STALL:
            raise InteriorStallEntryError(
                "the interior-stall step admits only a STALL entry; this state "
                f"demands {evidence.outcome.value} at {dt_s!r} s "
                f"({evidence.normalized_demand:.3e} times its own resolution). "
                "A resolved recession belongs to the ordinary same-cell solve.",
                before,
                evidence=evidence,
            )

        stall_seed = stall_seed_from_state(before) if seed is None else seed
        if not isinstance(stall_seed, InteriorStallSeed):
            raise TypeError("an interior stall needs a typed InteriorStallSeed")
        if len(stall_seed.candidate.equation_rank_vector()) != layout.unknown_count:
            raise ValueError("the interior-stall seed does not match the current rank")
        if len(stall_seed.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
            raise ValueError("one interior-stall Stefan scale is required per dry face")
        certify_front_is_pinned(before, stall_seed.candidate)

        pinned_front_z = before.transport.geometry.front.z
        chart = ci._coordinate_chart(before, boundary)  # noqa: SLF001
        if chart.front_z[1] != pinned_front_z:
            raise RuntimeError(
                "the same-cell chart's upper front bound is no longer the current "
                "front; the interior-stall rank exchange assumes it"
            )
        diagnostic_chart, diagnostic_interval = _stall_diagnostic_chart(before, chart)
        stefan_scales = stall_seed.stefan_flux_scales_mol_m2_s
        slot = rank.pinned_front_coordinate_slot
        consistency_row = rank.stall_consistency_row_index

        coordinates = encode_interior_stall_coordinates(
            stall_seed.candidate,
            chart,
            stefan_scales,
            layout=layout,
        )

        def decode(values: Sequence[float]) -> cut.CutTransportUnknowns:
            return decode_interior_stall_candidate(
                values,
                layout,
                chart,
                stefan_scales,
                pinned_front_z=pinned_front_z,
            )

        pinned_seed_candidate = decode(coordinates)
        last_candidate = pinned_seed_candidate
        initial_assembly = cut.assemble_backward_euler(
            before.transport,
            pinned_seed_candidate,
            dt_s,
            boundary,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        _certify_frozen_sweep(before, initial_assembly)
        scales = ci._residual_scales(  # noqa: SLF001
            before.transport,
            initial_assembly,
            dt_s,
        )
        scale_vector = np.asarray(scales.vector, dtype=float)
        independent_scales = np.delete(scale_vector, consistency_row)

        def residual(values: np.ndarray) -> np.ndarray:
            nonlocal evaluations, rejected_trials, last_candidate
            candidate = decode(values)
            last_candidate = candidate
            try:
                assembly = cut.assemble_backward_euler(
                    before.transport,
                    candidate,
                    dt_s,
                    boundary,
                    enforce_reduced_film_thresholds=False,
                )
            except cut.CutTransportStepError:
                rejected_trials += 1
                raise
            evaluations += 1
            rows = np.asarray(
                ci._datum_covariant_residual_vector(assembly.residuals),  # noqa: SLF001
                dtype=float,
            )
            return np.delete(rows, consistency_row) / independent_scales

        lower, upper = ci._coordinate_solver_bounds_for_chart(  # noqa: SLF001
            layout,
            before.controls.maximum_logit_magnitude,
            chart,
        )
        # The front slot is DELETED, not widened: the stall solve owns no front
        # box at all.  This is the mechanical statement of design section 2 A3.
        lower = np.delete(lower, slot)
        upper = np.delete(upper, slot)

        solution = optimize.least_squares(
            residual,
            coordinates,
            method=_NONLINEAR_METHOD,
            bounds=(lower, upper),
            ftol=before.controls.nonlinear_step_tolerance,
            xtol=before.controls.nonlinear_step_tolerance,
            # CC-1a's measured route: a zero-residual root is certified by the
            # residual contract itself, never by an early optimality stall.
            gtol=None,
            max_nfev=declared.maximum_function_evaluations,
            x_scale="jac",
        )
        optimizer_function_evaluations = int(solution.nfev)
        solution_coordinates = np.asarray(solution.x, dtype=float)
        candidate = decode(solution_coordinates)
        last_candidate = candidate
        signed = residual(solution_coordinates)
        last_scaled = tuple(float(value) for value in signed)
        maximum_scaled = max(abs(value) for value in last_scaled)
        condition_proxy = ci._jacobian_condition_proxy(solution.jac)  # noqa: SLF001

        if not solution.success:
            raise RuntimeError(
                f"the interior-stall root did not certify success: {solution.message}"
            )
        if np.any(solution.active_mask != 0):
            raise RuntimeError(
                "the interior-stall root touched a numerical coordinate guard; the "
                "open-chart root is not certified"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > declared.maximum_condition_proxy
        ):
            raise InteriorStallError(
                "the interior-stall Jacobian condition proxy is uncertified: "
                f"{condition_proxy:.3e}",
                before,
                refusal_code=InteriorStallRefusalCode.CONDITION_UNCERTIFIED,
            )
        if maximum_scaled > declared.nonlinear_residual_tolerance:
            raise InteriorStallError(
                "the interior-stall residual contract failed: "
                f"{maximum_scaled:.3e} > "
                f"{declared.nonlinear_residual_tolerance:.3e}",
                before,
                refusal_code=InteriorStallRefusalCode.NONLINEAR_CONTRACT_FAILED,
            )

        # Only a numerically certified root is entitled to exercise the
        # caller-owned reduced-film thresholds; the reassembly is
        # constitutively identical and merely enforces its audit.
        converged = cut.assemble_backward_euler(
            before.transport,
            candidate,
            dt_s,
            boundary,
            enforce_reduced_film_thresholds=True,
        )
        evaluations += 1
        certify_front_is_pinned(before, converged.candidate)
        _certify_frozen_sweep(before, converged)
        diagnostic_signed = tuple(
            value / scale
            for value, scale in zip(
                ci._datum_covariant_residual_vector(converged.residuals),  # noqa: SLF001
                scales.vector,
            )
        )
        production_independent = tuple(
            value
            for index, value in enumerate(diagnostic_signed)
            if index != consistency_row
        )
        last_scaled = production_independent
        maximum_scaled = max(abs(value) for value in production_independent)
        if maximum_scaled > declared.nonlinear_residual_tolerance:
            raise InteriorStallError(
                "the interior-stall production reassembly failed the unchanged "
                f"residual contract: {maximum_scaled:.3e} > "
                f"{declared.nonlinear_residual_tolerance:.3e}",
                before,
                refusal_code=InteriorStallRefusalCode.NONLINEAR_CONTRACT_FAILED,
            )
        if isinstance(boundary, ct.ReducedFilmPoreBoundary):
            film = converged.surface_film_audit
            if film is None or not (
                film.reduction_thresholds_enforced
                and film.fast_mass_guard_passed
                and film.validity_band_checked
                and film.monotone_unique_root_diagnostic_passed
            ):
                raise RuntimeError(
                    "the interior-stall candidate lacks a fully enforced film audit"
                )

        # THE STALL CONSISTENCY RESIDUAL.  The row the rank exchange removed,
        # read back at the converged point through SA-1's own frozen-front
        # certification: ``_demanded_front_volume_rate`` re-proves, bit-exactly,
        # that the sweep is zero and that the row IS the interface flux
        # difference before it inverts the jump.
        converged_rate, terms = cfd._demanded_front_volume_rate(  # noqa: SLF001
            before,
            converged,
        )
        consistency = abs(terms.row_mol_s) / scales.rh_water_mol_s

        # ``_accept_step`` builds objects and mutates nothing, so the gates
        # below still refuse with exact rollback: the caller gets ``before``
        # back and this step exists nowhere.
        step = ci._accept_step(  # noqa: SLF001
            before,
            converged,
            boundary,
            ci.CutStepSeed(
                stall_seed.candidate,
                stefan_scales,
                f"interior-stall seed: {stall_seed.label}",
            ),
            scales,
            scales,
            tuple(abs(value) for value in diagnostic_signed),
            dt_s,
            evaluations,
            rejected_trials,
            str(solution.message),
            condition_proxy,
            diagnostic_chart,
            None,
            _JACOBIAN_ROUTE_REASON,
            _NONLINEAR_METHOD,
        )
        if (
            step.ledger.maximum_step_ledger_residual > declared.ledger_tolerance
            or step.ledger.maximum_cumulative_ledger_residual
            > declared.ledger_tolerance
        ):
            raise InteriorStallError(
                "the interior-stall accepted/cumulative conservation ledger "
                f"exceeded its {declared.ledger_tolerance:.1e} contract",
                before,
                refusal_code=InteriorStallRefusalCode.LEDGER_CONTRACT_FAILED,
            )

        # THE EXIT VERDICT, in SA-1's own instrument.  Re-classifying the
        # COMMITTED state closes the loop the entry opened: the same criterion,
        # the same duration, the same non-widened tolerance.  Its RH water row
        # is the converged row above bit-for-bit - only the O9a donor is read at
        # the new time level - so this is the entry criterion applied to the
        # step's own result, not a second, differently scaled test.
        exit_evidence = cfd.evaluate_front_demand(
            step.after,
            boundary,
            dt_s,
            residual_tolerance=declared.nonlinear_residual_tolerance,
        )
        # STRICTNESS AGAINST ADVANCE, retained at the exit exactly as at the
        # entry: an outward drift is SA-1's typed refusal, never a stall.
        cfd.refuse_on_advance_demand(
            exit_evidence,
            before,
            surfaced_at=(
                "cut_interior_stall.solve_interior_stall (converged stall "
                "consistency row; design 2 A1/A2)"
            ),
        )
        if consistency > declared.nonlinear_residual_tolerance:
            raise InteriorStallConsistencyError(
                "the interior-stall consistency row exceeded the unchanged "
                f"residual contract: {consistency:.3e} > "
                f"{declared.nonlinear_residual_tolerance:.3e}",
                before,
                refusal_code=InteriorStallRefusalCode.STALL_CONSISTENCY_FAILED,
                last_candidate=candidate,
                stall_consistency_residual=consistency,
            )
        if exit_evidence.outcome is not cfd.FrontDemandOutcome.STALL:
            raise InteriorStallConsistencyError(
                "the committed pinned-front state is not stalled: its own "
                "Rankine-Hugoniot water row demands "
                f"{exit_evidence.demanded_front_volume_rate_m3_s!r} m3/s against a "
                f"{exit_evidence.resolution_m3_s!r} m3/s resolution "
                f"({exit_evidence.outcome.value})",
                before,
                refusal_code=InteriorStallRefusalCode.STALL_CONSISTENCY_FAILED,
                last_candidate=candidate,
                stall_consistency_residual=consistency,
            )

        ledger = InteriorStallLedger(
            dt_s=dt_s,
            pinned_front_z=pinned_front_z,
            pinned_front_radius_m=before.transport.geometry.front.radius_m,
            committed_front_z=step.after.transport.geometry.front.z,
            front_pinned_bit_exactly=(
                step.after.transport.geometry.front.z == pinned_front_z
                and converged.candidate.front_z == pinned_front_z
            ),
            swept_wet_to_dry_volume_m3=(
                converged.swept_geometry.wet_to_dry_volume_m3
            ),
            interface_swept_volume_rate_m3_s=(
                converged.swept_geometry.interface_swept_volume_rate_m3_s
            ),
            wet_hexane_sweep_identity_mol_s=(
                converged.ledger.wet_hexane_sweep_identity_mol_s
            ),
            layout=stall_layout,
            rank_audit=rank,
            same_cell_chart_front_z=chart.front_z,
            pinned_front_is_on_the_same_cell_chart_bound=(
                chart.front_z[1] == pinned_front_z
            ),
            front_diagnostic_interval_z=diagnostic_interval,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            nonlinear_optimizer_function_evaluations=(
                optimizer_function_evaluations
            ),
            nonlinear_message=str(solution.message),
            nonlinear_algorithm=_NONLINEAR_ALGORITHM,
            nonlinear_residual_basis=_NONLINEAR_RESIDUAL_BASIS,
            jacobian_route_reason=_JACOBIAN_ROUTE_REASON,
            used_sparse_jacobian=False,
            independent_scaled_residuals=production_independent,
            maximum_independent_scaled_residual=maximum_scaled,
            diagnostic_scaled_residuals=diagnostic_signed,
            stall_consistency_row_mol_s=terms.row_mol_s,
            stall_consistency_row_scale_mol_s=scales.rh_water_mol_s,
            normalized_stall_consistency_row=consistency,
            condition_proxy=condition_proxy,
            minimum_fractional_distance_to_bound=(
                step.ledger.minimum_fractional_distance_to_bound
            ),
            declared_nonlinear_residual_tolerance=(
                declared.nonlinear_residual_tolerance
            ),
            declared_ledger_tolerance=declared.ledger_tolerance,
            declared_maximum_condition_proxy=declared.maximum_condition_proxy,
            maximum_step_ledger_residual=(
                step.ledger.maximum_step_ledger_residual
            ),
            maximum_cumulative_ledger_residual=(
                step.ledger.maximum_cumulative_ledger_residual
            ),
            entry_demanded_front_volume_rate_m3_s=(
                evidence.demanded_front_volume_rate_m3_s
            ),
            entry_resolution_m3_s=evidence.resolution_m3_s,
            converged_consistency_demanded_rate_m3_s=converged_rate,
            exit_demanded_front_volume_rate_m3_s=(
                exit_evidence.demanded_front_volume_rate_m3_s
            ),
            exit_resolution_m3_s=exit_evidence.resolution_m3_s,
            exit_outcome=exit_evidence.outcome,
            accepted=True,
        )
        # A3: the accepted condensed surface root's typed evidence, when the
        # film-covered reassembly took the condensed branch.  ``None`` when
        # the root evaluated strictly above the locus (the byte-identical
        # legacy surface) - late in burn-off the film-covered surface lawfully
        # clears the locus before the film clears the latch.
        condensed_root_admission = None
        if film_covered_admission is not None:
            condensed_root_admission = (
                boundary.condensed_root_commit_channel.admission
            )
        return InteriorStallStep(
            before=before,
            boundary=boundary,
            seed=stall_seed,
            candidate=converged.candidate,
            step=step,
            entry_evidence=evidence,
            exit_evidence=exit_evidence,
            ledger=ledger,
            entered_after_same_cell_refusal=entered_after_same_cell_refusal,
            film_covered_admission=film_covered_admission,
            condensed_root_admission=condensed_root_admission,
        )
    except cfd.PrimaryDrainageAdvanceDemanded:
        raise
    except InteriorStallError as error:
        if error.rollback_state is None and isinstance(
            before,
            ci.CutIntegratorState,
        ):
            error.rollback_state = before
        raise
    except Exception as exc:
        if not isinstance(before, ci.CutIntegratorState):
            raise
        raise InteriorStallError(
            f"interior-stall step rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            nonlinear_optimizer_function_evaluations=(
                optimizer_function_evaluations
            ),
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            stall_consistency_residual=consistency,
            last_candidate=last_candidate,
            scaled_residuals=last_scaled,
        ) from exc


# ---------------------------------------------------------------------------
# CONTROLLER INTEGRATION - composed beside, in the SA-1 wrapper pattern
# ---------------------------------------------------------------------------
def advance_with_stall_admission(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    same_cell_seed: ci.CutStepSeed,
    *,
    stall_seed: InteriorStallSeed | None = None,
    stall_controls: InteriorStallControls | None = None,
    **macrostep_kwargs: object,
) -> StallAwareStep:
    """Take the one lawful step design section 2 allows at this state.

    The SA-1 verdict at ``(before, dt_s)`` selects the branch:

    * ``ADVANCE_DEMAND`` - :class:`cut_front_demand.
      PrimaryDrainageAdvanceDemanded`, the typed refusal naming PHY-014.  No
      step is taken and no envelope is changed; the refusal IS the Scope-B
      evidence record (design section 1).
    * ``STALL`` - :func:`solve_interior_stall` runs for the whole duration and
      commits with the front pinned.
    * ``RECESSION_DEMAND`` - the ordinary moving-front macrostep runs, wrapped
      in SA-1's own typed composition point so a refusal that turns out to be
      an advance is still NAMED.

    This is the entry the O9c re-march adopts.  It composes BESIDE the
    orchestrator: no F8 refusal site is edited, and a caller that does not
    adopt this entry sees exactly today's behaviour.
    """

    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("the stall-aware entry needs a cut integrator state")
    if not isinstance(same_cell_seed, ci.CutStepSeed):
        raise TypeError("the stall-aware entry needs a typed same-cell seed")
    if not isinstance(
        boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("the stall-aware entry requires a supported pore boundary")
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("the stall-aware duration must be positive and finite")

    evidence = cfd.evaluate_front_demand(before, boundary, dt_s)
    cfd.refuse_on_advance_demand(
        evidence,
        before,
        surfaced_at=(
            "cut_interior_stall.advance_with_stall_admission "
            "(pre-checked entry; design 2 A2)"
        ),
    )
    if evidence.outcome is cfd.FrontDemandOutcome.STALL:
        return StallAwareStep(
            before=before,
            boundary=boundary,
            dt_s=dt_s,
            entry_evidence=evidence,
            branch=StallAwareBranch.INTERIOR_STALL,
            stall_step=solve_interior_stall(
                before,
                dt_s,
                boundary,
                stall_seed,
                controls=stall_controls,
                entry_evidence=evidence,
            ),
            macrostep=None,
        )
    macrostep = cfd.advance_with_typed_front_demand(
        before,
        dt_s,
        boundary,
        same_cell_seed,
        **macrostep_kwargs,
    )
    return StallAwareStep(
        before=before,
        boundary=boundary,
        dt_s=dt_s,
        entry_evidence=evidence,
        branch=StallAwareBranch.ORDINARY_RECESSION,
        stall_step=None,
        macrostep=macrostep,
    )


def stall_step_for_refused_same_cell(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    refusal: BaseException,
    *,
    stall_seed: InteriorStallSeed | None = None,
    stall_controls: InteriorStallControls | None = None,
) -> InteriorStallStep:
    """Re-examine one refused same-cell step and stall it if it is a stall.

    This is the design's second entry ("after a same-cell refusal").  The
    refused step's own source state is re-classified by SA-1 at the SAME
    duration; a ``STALL`` verdict runs the pinned-front step, and ANY other
    verdict re-raises the caller's ORIGINAL refusal untouched, chained as the
    cause.  A refusal never acquires a stall it did not earn.
    """

    if not isinstance(refusal, BaseException):
        raise TypeError("the refused-step entry needs the original refusal")
    last_candidate = getattr(refusal, "last_candidate", None)
    if last_candidate is None:
        same_cell_error = getattr(refusal, "same_cell_error", None)
        last_candidate = getattr(same_cell_error, "last_candidate", None)
    evidence = cfd.evaluate_front_demand(
        before,
        boundary,
        dt_s,
        last_candidate=last_candidate,
    )
    cfd.refuse_on_advance_demand(
        evidence,
        before,
        surfaced_at=(
            "cut_interior_stall.stall_step_for_refused_same_cell "
            "(post-refusal entry; design 2 A2)"
        ),
        original_refusal=refusal,
    )
    if evidence.outcome is not cfd.FrontDemandOutcome.STALL:
        raise refusal
    return solve_interior_stall(
        before,
        dt_s,
        boundary,
        stall_seed,
        controls=stall_controls,
        entry_evidence=evidence,
        entered_after_same_cell_refusal=True,
    )


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------
def _certify_frozen_sweep(
    before: ci.CutIntegratorState,
    assembly: cut.CutTransportAssembly,
) -> None:
    """Refuse any assembly whose pinned front nonetheless swept volume.

    Three exact identities, never tolerances: the candidate front is the source
    front, the swept wet-to-dry volume is exactly zero, and the wet n-hexane
    sweep identity - which at ``q = 0`` reduces to "the wet region lost no
    n-hexane" - is exactly zero.
    """

    if assembly.candidate.front_z != before.transport.geometry.front.z:
        raise InteriorStallError(
            "the interior-stall assembly lost its pinned front",
            before,
            refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
        )
    swept = assembly.swept_geometry
    if swept.wet_to_dry_volume_m3 != 0.0 or (
        swept.interface_swept_volume_rate_m3_s != 0.0
    ):
        raise InteriorStallError(
            "an interior stall must sweep exactly zero volume; the assembly swept "
            f"{swept.wet_to_dry_volume_m3!r} m3",
            before,
            refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
        )
    if assembly.ledger.wet_hexane_sweep_identity_mol_s != 0.0:
        raise InteriorStallError(
            "a pinned front consumes no wet n-hexane; the sweep identity is "
            f"{assembly.ledger.wet_hexane_sweep_identity_mol_s!r}",
            before,
            refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
        )


def _stall_diagnostic_chart(
    before: ci.CutIntegratorState,
    chart: "ci._CoordinateChart",
) -> tuple["ci._CoordinateChart", tuple[float, float]]:
    """Return the accepted-step diagnostic chart and its front interval.

    WHY THIS EXISTS, stated plainly.  ``cut_integrator._bounded_fractions``
    (``:4185-4261``) raises when an accepted primitive sits ON a chart bound,
    and the ordinary same-cell chart's UPPER front bound IS the current front
    (``:1624``).  A pinned front therefore cannot be reported against that
    chart at all - the diagnostic would refuse the very state it is meant to
    describe.

    WHAT IT IS NOT.  It is not the D1 widened SOLVE chart that Revision B
    declined (design section 2 A3).  The stall coordinate vector has no front
    slot, its solver bounds have no front entry, and nothing is ever encoded
    into this interval.  The interval is the host CELL's own mesh-geometric
    span, so the reported fraction is a GEOMETRIC fact - how far the pinned
    front sits from its bounding master faces, which is also the ULP face-snap
    hazard the survey's H1 flags - and never a solver-boundary approach.
    """

    geometry = before.transport.geometry
    grid = geometry.master_grid
    cut_index = before.transport.layout.cut_cell_index
    inner_z = (grid.faces[cut_index] / grid.R) ** 3
    outer_z = (grid.faces[cut_index + 1] / grid.R) ** 3
    front_z = geometry.front.z
    if not inner_z < front_z < outer_z:
        raise InteriorStallError(
            "the pinned front is not strictly inside its own host cell",
            before,
            refusal_code=InteriorStallRefusalCode.FRONT_NOT_PINNED,
        )
    interval = (inner_z, outer_z)
    return replace(chart, front_z=interval), interval


def _certify_declared_contract_is_not_wider(
    controls: InteriorStallControls,
    solver_controls: ci.CutSolverControls,
) -> None:
    """Refuse any declared contract looser than the source state's own.

    The same non-widening enforcement CC-1a applies
    (``cut_cap_contact_integrator._certify_declared_contract_is_not_wider``),
    with the stall's own wording.
    """

    if not isinstance(solver_controls, ci.CutSolverControls):
        raise TypeError("an interior stall needs the state's typed solver controls")
    if (
        controls.nonlinear_residual_tolerance
        > solver_controls.nonlinear_residual_tolerance
    ):
        raise InteriorStallError(
            "an interior stall may not widen the residual contract: "
            f"{controls.nonlinear_residual_tolerance:.3e} > "
            f"{solver_controls.nonlinear_residual_tolerance:.3e}",
            refusal_code=InteriorStallRefusalCode.TOLERANCE_WIDENED,
        )
    if controls.ledger_tolerance > solver_controls.ledger_tolerance:
        raise InteriorStallError(
            "an interior stall may not widen the ledger contract: "
            f"{controls.ledger_tolerance:.3e} > "
            f"{solver_controls.ledger_tolerance:.3e}",
            refusal_code=InteriorStallRefusalCode.TOLERANCE_WIDENED,
        )
    if controls.maximum_condition_proxy > solver_controls.maximum_condition_proxy:
        raise InteriorStallError(
            "an interior stall may not widen the condition ceiling: "
            f"{controls.maximum_condition_proxy:.3e} > "
            f"{solver_controls.maximum_condition_proxy:.3e}",
            refusal_code=InteriorStallRefusalCode.TOLERANCE_WIDENED,
        )
    if (
        controls.maximum_function_evaluations
        > solver_controls.maximum_function_evaluations
    ):
        raise InteriorStallError(
            "an interior stall may not exceed the state's evaluation budget: "
            f"{controls.maximum_function_evaluations} > "
            f"{solver_controls.maximum_function_evaluations}",
            refusal_code=InteriorStallRefusalCode.TOLERANCE_WIDENED,
        )


__all__ = [
    "INTERIOR_STALL_LEDGER_TOLERANCE",
    "INTERIOR_STALL_MAXIMUM_CONDITION_PROXY",
    "INTERIOR_STALL_RESIDUAL_TOLERANCE",
    "SA2_DESIGN_BASIS",
    "SA2_MIRRORED_CONSTRUCTION",
    "SA2_RANK_EXCHANGE_PRECEDENT",
    "SA2_SURVEY_EVIDENCE",
    "FilmCoveredStallAdmission",
    "InteriorStallConsistencyError",
    "InteriorStallControls",
    "InteriorStallEntryError",
    "InteriorStallError",
    "InteriorStallLayout",
    "InteriorStallLedger",
    "InteriorStallRankAudit",
    "InteriorStallRefusalCode",
    "InteriorStallSeed",
    "InteriorStallStep",
    "StallAwareBranch",
    "StallAwareStep",
    "advance_with_stall_admission",
    "audit_interior_stall_rank",
    "certify_front_is_pinned",
    "decode_interior_stall_candidate",
    "encode_interior_stall_coordinates",
    "interior_stall_layout",
    "pinned_front_coordinate_slot",
    "solve_interior_stall",
    "stall_consistency_row_index",
    "stall_seed_from_state",
    "stall_step_for_refused_same_cell",
]
