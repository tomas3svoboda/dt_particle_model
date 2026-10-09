r"""SA-1/A2 - the typed advance-demand outcome for the primary-drainage front.

Design basis: ``docs/GT_PS2_O10B_FRONT_CHART_DESIGN_REVISION_B_2026-08-21.md``
section 2 A2 ("THE TYPED ADVANCE REFUSAL"), ratified 2026-08-21 for Scope A.
Every file:line claim repeated below is anchored in the read-only survey
``docs/evidence/GT_PS2_O10B_ENDPOINT_MONOTONICITY_SURVEY_2026-08-21.md``.

WHAT THIS MODULE IS FOR.  The survey found two places where a front that the
corrected physics wants to move OUTWARD dies untyped:

* D2.2 - an advancing front reaching the outer master face raises
  ``cut_transport.py:1644-1648``
  ``CutTransportTopologyError("front reached or crossed a master face; event
  restart required...")``, which the atomic orchestrator turns into an
  ``_AtomicRejection(TOPOLOGY_HANDOFF_FAILED)``.  An untyped rejection, not an
  event.
* F1/F8 - inside the residual itself, ``cut_transport.assemble_backward_euler``
  calls ``cg.swept_cut_geometry`` (``cut_transport.py:1649``), whose
  ``cut_geometry.py:418-419`` drainage veto aborts the first outward trial, so
  the same-cell solve dies as a generic ``CutIntegratorStepError`` - the
  cond-inf saturation family.

Neither death names the boundary it hit.  This module converts both into the
house's typed refusal, :class:`PrimaryDrainageAdvanceDemanded`, carrying the
demanded front volumetric rate as evidence - the exact analog of the K2
kernel's typed pre-event NO-GO whose phase condition is
``surface_active_set.PhaseResidualClassification.ABSENT_APPEARANCE_DEMANDED``
(``qsc_k2_native_source_kernel.K2NativeWaterPreEventNoGo``), the evidence that
triggered the ruled re-wetting law.  The envelope speaks; nothing clamps.

A BESIDE MODULE, NOT AN EDIT - AND THE FROZEN VETOES ARE NOT TOUCHED.  PHY-014
/ PHY-023 / PHY-036 are FROZEN physics rulings (survey F6), enforced verbatim
at ``cut_geometry.py:418-419``, ``front.py:279-280``,
``surface_active_set.py:704-765`` and ``coupled_pore.py:2506-2516``.  This
module edits none of them, widens no chart, and relaxes no tolerance.  The
drainage-latch site was COMPLETED to PHY-014's own carve-out sentence by the
A5 ruling (2026-08-22, ``docs/GT_PS2_A5_LATCH_CARVEOUT_RULING_PACKET_
2026-08-22.md``): a lawfully committed, identity-bound A1/A2-lineage film
admits recession while it evolves; the no-lineage refusal and the three named
prohibitions (s=R reset, X_f recreation, internal-core nucleation) are
byte-identical to the pre-A5 latch.  It
composes at the CONTROLLER/ORCHESTRATOR level exactly as design section 2 A2
prescribes: :func:`assemble_with_typed_front_demand` and
:func:`advance_with_typed_front_demand` wrap call sites that already exist and
already fail closed, and they re-raise a NAMED refusal in place of an unnamed
one.  A caller that does not adopt the wrappers sees exactly today's behaviour.

WHY THE COMPOSITION IS A WRAPPER AND NOT AN EDIT INSIDE THE ORCHESTRATOR.
Three reasons, all structural: (i) the demand evaluation costs two extra
assemblies, and the orchestrator's failure path is walked inside bisection
loops where that cost is not free; (ii) importing this module from
``cut_event_orchestrator`` would close an import cycle, since the evaluation
needs the orchestrator's own typed failure payload; (iii) the refusal is a
CALLER-OWNED routing decision - the Scope-B evidence stream (design section 1)
is collected by the march controller, not manufactured inside the solver.

THE CRITERION IS DERIVATIVE-FREE, ON ASSEMBLED FLUXES.  It is the
``cut_cap_active_continuation.evaluate_cap_contact_demand`` pattern: freeze the
refused step's SOURCE state, assemble one unchanged
:func:`cut_transport.assemble_backward_euler` at it, and read the demanded
front rate off the Rankine-Hugoniot water row under the O9a-corrected donor.
With the front pinned at its own position the swept-volume rate is exactly
zero, so the RH water row IS ``A_Gamma * (J_dry - J_wet)`` bit-exactly
(certified, not assumed), and the row's own jump condition

    ``q_Gamma * (C_dry - C_donor) = A_Gamma * (J_dry - J_wet)``

inverts to the demanded rate

    ``q_demand = A_Gamma * (J_dry - J_wet) / (C_dry - C_donor)``.

``interface_swept_volume_rate_m3_s`` is the signed ``dV_wet/dt``
(``cut_geometry.py:243-246``): NEGATIVE is primary-drainage recession, POSITIVE
is the forbidden advance.

THE RESOLUTION IS NOT AN INVENTED EPSILON.  It is the SAME acceptance band the
refused step itself was judged against, mapped through the SAME jump:

    ``resolution = nonlinear_residual_tolerance * rh_water_row_scale / |C_dry - C_donor|``

with ``nonlinear_residual_tolerance`` the state's own declared contract
(``cut_integrator.CutSolverControls``) and ``rh_water_row_scale`` the
assembly's own ``cut_integrator._residual_scales(...).rh_water_mol_s`` at the
refused step's duration.  A demand inside that band is a front motion the
assembly could not have resolved, which is precisely a STALL.  Note the
consequence, reported rather than hidden: the row scale carries ``inventory /
dt`` terms, so the RESOLUTION depends on the refused step's duration while the
demanded RATE does not - both facts are certified on every evaluation.

DECLARED LIMITS, STATED RATHER THAN PAPERED OVER.

* This module CLASSIFIES; it solves nothing.  ``ADVANCE_DEMAND`` is a refusal
  with evidence, never a step.  Whether the primary-drainage envelope re-opens
  for a sustained resolved advance is the owner's Scope-B question (design
  section 1), and this module is the instrument that collects the evidence for
  it - not an answer to it.
* The criterion freezes storage, so it excludes BY CONSTRUCTION the ALE storage
  term a genuinely moving front contributes; it answers "what front rate does
  this state's assembled interface balance demand", not "where will the front
  be at the end of the step".  This is the identical, deliberate scope limit
  CC-2a declares for its own demand criterion.
* The evaluation is defined only where the RH water jump ``C_dry - C_donor`` is
  finite and non-zero.  A vanishing jump is reported as a typed refusal, never
  as a zero rate.
* The interior STALL BRANCH itself (design section 2 A1 - the six F8 refusal
  sites gaining a typed stall admission) is SA-2 and is not built here.
  :attr:`FrontDemandOutcome.STALL` names the classification; it performs no
  step.

Every claim flag on every payload in this module is ``False``.  Nothing here
makes the particle model physically qualifying, plant predictive, or
production wired.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_event_orchestrator as ceo
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut


SA1_DESIGN_BASIS = "docs/GT_PS2_O10B_FRONT_CHART_DESIGN_REVISION_B_2026-08-21.md"
SA1_SURVEY_EVIDENCE = (
    "docs/evidence/GT_PS2_O10B_ENDPOINT_MONOTONICITY_SURVEY_2026-08-21.md"
)
FROZEN_PHYSICS_AUTHORITY = "release/physics_decisions.yaml"


@dataclass(frozen=True)
class FrozenPhysicsCitation:
    """One frozen decision named and quoted verbatim from its own authority.

    ``statement`` is the decision's ``selection.statement`` copied character for
    character out of ``release/physics_decisions.yaml``; ``anchor`` addresses
    that decision BY ID in the form ``<authority>#<decision id>``, never by line
    range, so inserting a decision into the register cannot silently repoint it.
    A refusal that names a frozen boundary must be able to show the boundary's
    own words, not a paraphrase of them.
    """

    decision_id: str
    authority: str
    anchor: str
    statement: str
    enforcement_sites: tuple[str, ...]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        for label, value in (
            ("decision id", self.decision_id),
            ("authority", self.authority),
            ("anchor", self.anchor),
            ("statement", self.statement),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"a frozen-physics citation needs a {label}")
        if not self.enforcement_sites or not all(
            isinstance(site, str) and site.strip() for site in self.enforcement_sites
        ):
            raise ValueError("a frozen-physics citation needs its enforcement sites")


#: PHY-014, verbatim from ``release/physics_decisions.yaml#PHY-014``
#: (``decisions.PHY-014.selection.statement``).
PHY_014 = FrozenPhysicsCitation(
    decision_id="PHY-014",
    authority=FROZEN_PHYSICS_AUTHORITY,
    anchor="release/physics_decisions.yaml#PHY-014",
    statement=(
        "The qualified model is a primary-drainage topology. Once the central "
        "mobile-n-hexane core has receded, no supersaturation or condensation "
        "event may set s=R, recreate X_f, or nucleate a new internal core. "
        "Existing admissible full-core or attached-liquid states may evolve "
        "normally."
    ),
    enforcement_sites=(
        "cut_geometry.py:418-419",
        "front.py:279-280",
        # The latch, completed to the statement's own "existing admissible
        # ... attached-liquid states may evolve normally" sentence by the A5
        # ruling (2026-08-22): the committed A1/A2-lineage film admits
        # recession; every other refusal is byte-identical.
        "surface_active_set.py:704-765",
        "coupled_pore.py:2506-2516",
    ),
)

#: PHY-023, verbatim from ``release/physics_decisions.yaml#PHY-023`` - the
#: second ``decisions.PHY-023.admissibility`` clause, which binds the memoryless
#: primary-drainage area closure to the same monotone envelope.
PHY_023 = FrozenPhysicsCitation(
    decision_id="PHY-023",
    authority=FROZEN_PHYSICS_AUTHORITY,
    anchor="release/physics_decisions.yaml#PHY-023",
    statement=(
        "Normal qualified material trajectories are monotonic primary "
        "drainage. Rewetting or later n-hexane condensation on a receded core "
        "remains outside this closure and under PHY-036."
    ),
    enforcement_sites=("cut_geometry.py:462-463", "dry_region.py:455-456"),
)

#: PHY-036, verbatim from ``release/physics_decisions.yaml#PHY-036``
#: (``decisions.PHY-036.selection.statement``).
PHY_036 = FrozenPhysicsCitation(
    decision_id="PHY-036",
    authority=FROZEN_PHYSICS_AUTHORITY,
    anchor="release/physics_decisions.yaml#PHY-036",
    statement=(
        "For qualified material trajectories, X_f>0 implies s=R and the "
        "particle-gas pathway is inactive. Drain X_f conservatively to exactly "
        "zero before activating front recession. A trajectory requiring later "
        "n-hexane condensation or rewetting on a receded/dry core is outside "
        "the qualified morphology envelope and must be flagged rather than "
        "modeled with an ungrounded history closure."
    ),
    enforcement_sites=("sphere_coupling.py:161-162", "sphere_coupling.py:211-212"),
)

#: The three frozen rulings a forbidden-advance refusal must name.  PHY-036's
#: own text supplies this module's required response verbatim: such a
#: trajectory "must be flagged rather than modeled with an ungrounded history
#: closure".  Flagging it is exactly what this module does.
PRIMARY_DRAINAGE_CITATIONS = (PHY_014, PHY_023, PHY_036)


class FrontDemandOutcome(str, Enum):
    """The typed three-way verdict of design section 2 A2/A3.

    The verdict is a PHYSICAL statement about the assembled Rankine-Hugoniot
    water row and never a conditioning proxy: per the Amendment-18 principle a
    representability failure alone may not flip topology, so nothing about the
    Jacobian, its condition, or any solver diagnostic enters here.
    """

    #: ``q_demand < -resolution`` - strict inward sweep, the only branch that
    #: still solves a moving front (design section 2 A3).
    RECESSION_DEMAND = "recession_demand"
    #: ``|q_demand| <= resolution`` - the demanded motion is below the refused
    #: step's own acceptance band.  SA-2 owns the stall BRANCH; this is the
    #: classification only.
    STALL = "stall"
    #: ``q_demand > +resolution`` - the forbidden outward demand.  PHY-014.
    ADVANCE_DEMAND = "advance_demand"


class FrontDemandRefusalCode(str, Enum):
    """Named refusal reasons for the front-demand ceremony."""

    FRONT_NOT_FROZEN = "front_not_frozen_on_the_demand_probe"
    JUMP_NOT_USABLE = "rankine_hugoniot_water_jump_not_usable"
    DURATION_DEPENDENT_RATE = "duration_dependent_demanded_rate"
    TOLERANCE_WIDENED = "tolerance_widened"
    ADVANCE_DEMANDED = "primary_drainage_advance_demanded"


class FrontDemandError(RuntimeError):
    """Base of the typed front-demand family; carries the exact source state.

    ``rollback_state`` is the immutable state the evaluation was asked about.
    A caller that catches this has lost nothing: no field was mutated, no
    candidate adopted, no tolerance relaxed.
    """

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState | None = None,
        *,
        refusal_code: FrontDemandRefusalCode | None = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.refusal_code = refusal_code


class FrontDemandEvaluationError(FrontDemandError):
    """The demand criterion could not be evaluated on the given state."""


class PrimaryDrainageAdvanceDemanded(FrontDemandError):
    """TYPED REFUSAL: the assembled interface balance demands an outward front.

    This is design section 2 A2's refusal and the analog of the K2 kernel's
    ``ABSENT_APPEARANCE_DEMANDED`` pre-event NO-GO: the state is not repaired,
    not clamped and not solved - it is NAMED, together with the frozen boundary
    it reached and the rate it demanded.  PHY-036's own text requires exactly
    this response, that such a trajectory "must be flagged rather than modeled
    with an ungrounded history closure".

    A caught instance is the Scope-B evidence record.  It changes no envelope
    by itself: the owner's PHY-014 question opens only on a SUSTAINED, resolved
    advance across the O9c re-march (design section 1), never on one step.
    """

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState | None = None,
        *,
        evidence: "FrontDemandEvidence",
        surfaced_at: str,
        original_refusal: BaseException | None = None,
    ) -> None:
        if not isinstance(evidence, FrontDemandEvidence):
            raise TypeError("the typed advance refusal requires typed demand evidence")
        if evidence.outcome is not FrontDemandOutcome.ADVANCE_DEMAND:
            raise ValueError(
                "the typed advance refusal may only carry ADVANCE_DEMAND evidence; "
                f"got {evidence.outcome.value}"
            )
        if not isinstance(surfaced_at, str) or not surfaced_at.strip():
            raise ValueError("the typed advance refusal needs its surfacing site")
        super().__init__(
            message,
            rollback_state,
            refusal_code=FrontDemandRefusalCode.ADVANCE_DEMANDED,
        )
        self.evidence = evidence
        self.surfaced_at = surfaced_at
        self.original_refusal = original_refusal
        self.citations = PRIMARY_DRAINAGE_CITATIONS

    @property
    def demanded_front_volume_rate_m3_s(self) -> float:
        """Return the demanded signed ``dV_wet/dt``; positive is the advance."""

        return self.evidence.demanded_front_volume_rate_m3_s

    @property
    def resolution_m3_s(self) -> float:
        return self.evidence.resolution_m3_s

    @property
    def cited_decision_ids(self) -> tuple[str, ...]:
        return tuple(citation.decision_id for citation in self.citations)


@dataclass(frozen=True)
class FrontDemandEvidence:
    """Assembled-flux evidence and the typed verdict of design section 2 A2.

    Every number here is read off ONE unchanged
    :func:`cut_transport.assemble_backward_euler` assembly at the refused
    step's source state with the front FROZEN at its own position, so the swept
    volume rate is exactly zero, the RH water row is exactly the interface flux
    difference, and the demanded rate is derivative-free and duration
    invariant.  All three of those facts are CERTIFIED on construction rather
    than asserted in prose.
    """

    #: The signed ``dV_wet/dt`` the RH water row demands, m3/s.  Negative is
    #: primary-drainage recession, positive is the forbidden advance.
    demanded_front_volume_rate_m3_s: float
    #: ``residual_tolerance * rh_water_row_scale / |jump|``; see the module
    #: docstring.  Never an invented epsilon.
    resolution_m3_s: float
    outcome: FrontDemandOutcome
    probe_duration_s: float
    invariance_probe_duration_s: float
    invariance_demanded_rate_m3_s: float
    front_z: float
    front_radius_m: float
    interface_area_m2: float
    frozen_swept_volume_rate_m3_s: float
    dry_interface_water_flux_mol_m2_s: float
    wet_interface_water_flux_mol_m2_s: float
    interface_water_flux_difference_mol_s: float
    dry_total_water_concentration_mol_m3: float
    corrected_wet_donor_concentration_mol_m3: float
    rankine_hugoniot_water_jump_mol_m3: float
    rankine_hugoniot_water_row_mol_s: float
    rankine_hugoniot_row_scale_mol_s: float
    declared_residual_tolerance: float
    normalized_rankine_hugoniot_row: float
    #: The refused solve's own last candidate front position when the failure
    #: carried one, else ``None``.  EVIDENCE ONLY - it never enters the verdict,
    #: which is decided on assembled fluxes alone.
    last_candidate_front_z: float | None
    last_candidate_front_z_change: float | None
    front_frozen_exactly: bool = field(default=True, init=False)
    storage_frozen_exactly: bool = field(default=True, init=False)
    decided_on_assembled_fluxes_only: bool = field(default=True, init=False)
    used_conditioning_proxy: bool = field(default=False, init=False)
    uses_corrected_o9a_front_donor: bool = field(default=True, init=False)
    rate_duration_invariance_certified: bool = field(default=True, init=False)
    #: Reported, not hidden: the row SCALE carries ``inventory / dt`` terms, so
    #: the resolution is a property of the refused step's duration.
    resolution_depends_on_probe_duration: bool = field(default=True, init=False)
    stall_branch_implemented: bool = field(default=False, init=False)
    advance_branch_implemented: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)
    plant_predictive: bool = field(default=False, init=False)
    production_wired: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if type(self.outcome) is not FrontDemandOutcome:
            raise TypeError("the front-demand verdict must use the typed trichotomy")
        if self.frozen_swept_volume_rate_m3_s != 0.0:
            raise FrontDemandEvaluationError(
                "the front-demand criterion is only defined with the front frozen; "
                f"the probe swept {self.frozen_swept_volume_rate_m3_s!r} m3/s",
                refusal_code=FrontDemandRefusalCode.FRONT_NOT_FROZEN,
            )
        if not math.isfinite(self.resolution_m3_s) or self.resolution_m3_s < 0.0:
            raise ValueError("the stall resolution must be finite and non-negative")
        if self.demanded_front_volume_rate_m3_s != (
            self.invariance_demanded_rate_m3_s
        ):
            raise FrontDemandEvaluationError(
                "the frozen-front demanded rate is not duration invariant; it is "
                "therefore not the derivative-free criterion design section 2 names",
                refusal_code=FrontDemandRefusalCode.DURATION_DEPENDENT_RATE,
            )
        expected = classify_front_demand(
            self.demanded_front_volume_rate_m3_s,
            self.resolution_m3_s,
        )
        if expected is not self.outcome:
            raise ValueError("the front-demand verdict disagrees with its own evidence")

    @property
    def demands_forbidden_advance(self) -> bool:
        return self.outcome is FrontDemandOutcome.ADVANCE_DEMAND

    @property
    def normalized_demand(self) -> float:
        """Return the demanded rate in units of its own stall resolution."""

        if self.resolution_m3_s == 0.0:
            return math.inf if self.demanded_front_volume_rate_m3_s != 0.0 else 0.0
        return self.demanded_front_volume_rate_m3_s / self.resolution_m3_s

    @property
    def citations(self) -> tuple[FrozenPhysicsCitation, ...]:
        return PRIMARY_DRAINAGE_CITATIONS


def classify_front_demand(
    demanded_rate_m3_s: float,
    resolution_m3_s: float,
) -> FrontDemandOutcome:
    """Return the typed three-way verdict for one signed demanded front rate.

    The band is CLOSED on the stall side, exactly as
    ``cut_cap_active_continuation._classify_demand`` closes its own: a demand
    sitting exactly on the resolution is a motion the assembly could not have
    resolved, so it is a stall and not a refusal.  ``cut_face_tangent``'s frozen
    on-face gate makes the same choice (``s_dot <= 0`` admits zero).
    """

    if not math.isfinite(demanded_rate_m3_s):
        raise FrontDemandEvaluationError("the demanded front rate must be finite")
    if not math.isfinite(resolution_m3_s) or resolution_m3_s < 0.0:
        raise FrontDemandEvaluationError(
            "the stall resolution must be finite and non-negative"
        )
    if demanded_rate_m3_s > resolution_m3_s:
        return FrontDemandOutcome.ADVANCE_DEMAND
    if demanded_rate_m3_s < -resolution_m3_s:
        return FrontDemandOutcome.RECESSION_DEMAND
    return FrontDemandOutcome.STALL


def evaluate_front_demand(
    state: ci.CutIntegratorState,
    boundary: ct.PoreBoundary,
    dt_s: float,
    *,
    last_candidate: cut.CutTransportUnknowns | None = None,
    residual_tolerance: float | None = None,
) -> FrontDemandEvidence:
    """Return the typed front-demand verdict at one refused step's source state.

    ``dt_s`` is the duration the refused step was judged against; it sets the
    STALL RESOLUTION through the assembly's own residual scales and nothing
    else.  ``last_candidate`` is the refused solve's own last trial
    (``cut_integrator.CutIntegratorStepError.last_candidate``) when the failure
    carried one; it is recorded as evidence and never enters the verdict.

    The evaluation freezes the state: the probe candidate carries the state's
    own primitives and its own front position, so the swept volume is exactly
    zero, every inventory change is exactly zero, and the RH water row IS the
    interface flux difference.  Duration invariance of the DEMANDED RATE is
    certified against a half-duration re-assembly, not assumed.
    """

    if not isinstance(state, ci.CutIntegratorState):
        raise TypeError("the front-demand criterion needs a cut integrator state")
    if not isinstance(
        boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("the front-demand criterion requires a supported pore boundary")
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("the front-demand probe duration must be positive and finite")
    declared_tolerance = (
        state.controls.nonlinear_residual_tolerance
        if residual_tolerance is None
        else float(residual_tolerance)
    )
    if not math.isfinite(declared_tolerance) or not (
        0.0 < declared_tolerance <= state.controls.nonlinear_residual_tolerance
    ):
        raise FrontDemandError(
            "the stall resolution may not be declared looser than the state's own "
            f"{state.controls.nonlinear_residual_tolerance:.1e} residual contract",
            state,
            refusal_code=FrontDemandRefusalCode.TOLERANCE_WIDENED,
        )
    if last_candidate is not None and not isinstance(
        last_candidate,
        cut.CutTransportUnknowns,
    ):
        raise TypeError("the last-candidate evidence must be CutTransportUnknowns")

    transport = state.transport
    candidate = cut.candidate_from_state(
        transport,
        dry_total_stefan_fluxes_mol_m2_s=state.last_total_stefan_fluxes_mol_m2_s,
        interface_temperature_k=state.last_interface_temperature_k,
    )
    assembly = cut.assemble_backward_euler(
        transport,
        candidate,
        dt_s,
        boundary,
        enforce_reduced_film_thresholds=False,
    )
    rate, terms = _demanded_front_volume_rate(state, assembly)

    invariance_duration = 0.5 * dt_s
    invariance_assembly = cut.assemble_backward_euler(
        transport,
        candidate,
        invariance_duration,
        boundary,
        enforce_reduced_film_thresholds=False,
    )
    invariance_rate, _ = _demanded_front_volume_rate(state, invariance_assembly)

    scales = ci._residual_scales(transport, assembly, dt_s)  # noqa: SLF001
    row_scale = scales.rh_water_mol_s
    # THE RESOLUTION, construction-stated.  ``|rh_water| <= tolerance *
    # row_scale`` is the SAME acceptance band ``cut_integrator`` judges the row
    # by; dividing it by the SAME jump that maps the row into a volumetric rate
    # carries that band, unchanged, onto the front.  No new tolerance is
    # declared and none is widened.
    resolution = declared_tolerance * row_scale / abs(terms.jump_mol_m3)
    front = transport.geometry.front
    return FrontDemandEvidence(
        demanded_front_volume_rate_m3_s=rate,
        resolution_m3_s=resolution,
        outcome=classify_front_demand(rate, resolution),
        probe_duration_s=dt_s,
        invariance_probe_duration_s=invariance_duration,
        invariance_demanded_rate_m3_s=invariance_rate,
        front_z=front.z,
        front_radius_m=front.radius_m,
        interface_area_m2=terms.area_m2,
        frozen_swept_volume_rate_m3_s=terms.frozen_swept_rate_m3_s,
        dry_interface_water_flux_mol_m2_s=terms.dry_flux_mol_m2_s,
        wet_interface_water_flux_mol_m2_s=terms.wet_flux_mol_m2_s,
        interface_water_flux_difference_mol_s=terms.flux_difference_mol_s,
        dry_total_water_concentration_mol_m3=terms.dry_concentration_mol_m3,
        corrected_wet_donor_concentration_mol_m3=terms.donor_concentration_mol_m3,
        rankine_hugoniot_water_jump_mol_m3=terms.jump_mol_m3,
        rankine_hugoniot_water_row_mol_s=terms.row_mol_s,
        rankine_hugoniot_row_scale_mol_s=row_scale,
        declared_residual_tolerance=declared_tolerance,
        normalized_rankine_hugoniot_row=abs(terms.row_mol_s) / row_scale,
        last_candidate_front_z=(
            None if last_candidate is None else last_candidate.front_z
        ),
        last_candidate_front_z_change=(
            None
            if last_candidate is None
            else math.fsum((last_candidate.front_z, -front.z))
        ),
    )


@dataclass(frozen=True)
class _RankineHugoniotTerms:
    """The assembled RH water-row terms one frozen-front probe exposes."""

    area_m2: float
    frozen_swept_rate_m3_s: float
    dry_flux_mol_m2_s: float
    wet_flux_mol_m2_s: float
    flux_difference_mol_s: float
    dry_concentration_mol_m3: float
    donor_concentration_mol_m3: float
    jump_mol_m3: float
    row_mol_s: float


def _demanded_front_volume_rate(
    state: ci.CutIntegratorState,
    assembly: cut.CutTransportAssembly,
) -> tuple[float, _RankineHugoniotTerms]:
    """Return ``(q_demand, terms)`` from one frozen-front assembly.

    The row is formed in ``cut_transport._residual_blocks`` (``:3185-3195``) as

        ``fsum(A*J_dry, -A*J_wet, -q*C_dry, +q*C_donor)``

    with ``C_donor`` the O9a-corrected consumed-piece BULK concentration.  On a
    frozen front ``q`` is exactly zero, so the last two terms are signed zeros
    and cannot change a correctly rounded ``fsum``: the row must therefore equal
    the pure flux difference BIT-EXACTLY.  Testing that equality is the exact
    test - back-subtracting the sweep out of an already rounded row would only
    measure that rounding.  The same reasoning is
    ``cut_cap_active_continuation.evaluate_cap_contact_demand``'s storage-frozen
    certification, transplanted from the wet-water row onto the jump row.
    """

    area = 4.0 * math.pi * assembly.candidate_geometry.front.radius_m**2
    swept_rate = assembly.swept_geometry.interface_swept_volume_rate_m3_s
    dry_flux = assembly.dry_face_fluxes[0].component.conserved_water_flux_mol_m2_s
    wet_flux = assembly.wet_face_fluxes[-1].retained_water_flux_mol_m2_s
    row = assembly.residuals.rh_water_mol_s
    difference = math.fsum((area * dry_flux, -area * wet_flux))
    if swept_rate != 0.0 or row != difference:
        raise FrontDemandEvaluationError(
            "the front-demand criterion is only defined on a frozen front; the "
            f"probe swept {swept_rate!r} m3/s and its RH water row {row!r} "
            f"differs from the pure interface flux difference {difference!r}",
            state,
            refusal_code=FrontDemandRefusalCode.FRONT_NOT_FROZEN,
        )
    # The O9a-corrected wet-side donor, read from ``cut_integrator``'s own
    # reconstruction so this criterion cannot silently drift away from the
    # donor the solver books (``cut_integrator._assembly_sweep_donor_concentration``,
    # ruled in docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md).
    donor = ci._assembly_sweep_donor_concentration(assembly.before)  # noqa: SLF001
    dry_concentration = assembly.interface.dry.total_water_concentration_mol_m3
    jump = math.fsum((dry_concentration, -donor))
    if jump == 0.0 or not math.isfinite(jump) or not math.isfinite(difference):
        raise FrontDemandEvaluationError(
            "the Rankine-Hugoniot water jump is unusable, so no front rate is "
            f"demanded by this state: jump {jump!r}, flux difference "
            f"{difference!r}",
            state,
            refusal_code=FrontDemandRefusalCode.JUMP_NOT_USABLE,
        )
    return difference / jump, _RankineHugoniotTerms(
        area_m2=area,
        frozen_swept_rate_m3_s=swept_rate,
        dry_flux_mol_m2_s=dry_flux,
        wet_flux_mol_m2_s=wet_flux,
        flux_difference_mol_s=difference,
        dry_concentration_mol_m3=dry_concentration,
        donor_concentration_mol_m3=donor,
        jump_mol_m3=jump,
        row_mol_s=row,
    )


def refuse_on_advance_demand(
    evidence: FrontDemandEvidence,
    state: ci.CutIntegratorState,
    *,
    surfaced_at: str,
    original_refusal: BaseException | None = None,
) -> None:
    """Raise the typed refusal when, and only when, the demand is an advance.

    A ``RECESSION_DEMAND`` or ``STALL`` verdict returns silently: PHY-014 is
    strict against ADVANCE and against nothing else, so a receding or stalled
    front keeps exactly the refusal its own solver produced.
    """

    if not isinstance(evidence, FrontDemandEvidence):
        raise TypeError("the advance refusal needs typed demand evidence")
    if evidence.outcome is not FrontDemandOutcome.ADVANCE_DEMAND:
        return
    names = ", ".join(citation.decision_id for citation in PRIMARY_DRAINAGE_CITATIONS)
    # ``from`` keeps the original refusal as the CAUSE: the typed outcome
    # renames the death, it never hides the site that produced it.
    raise PrimaryDrainageAdvanceDemanded(
        "primary-drainage front advance demanded at "
        f"{surfaced_at}: the assembled Rankine-Hugoniot water row demands "
        f"dV_wet/dt = {evidence.demanded_front_volume_rate_m3_s!r} m3/s, "
        f"{evidence.normalized_demand:.3e} times its own stall resolution "
        f"{evidence.resolution_m3_s!r} m3/s. The frozen monotone "
        f"primary-drainage envelope ({names}) admits no outward front: "
        f"{PHY_014.statement} This trajectory is FLAGGED, not modelled - "
        f"{PHY_036.decision_id} requires exactly that response.",
        state,
        evidence=evidence,
        surfaced_at=surfaced_at,
        original_refusal=original_refusal,
    ) from original_refusal


# ---------------------------------------------------------------------------
# COMPOSITION POINT (a) - the survey's D2.2 untyped outer-face death
# ---------------------------------------------------------------------------
def assemble_with_typed_front_demand(
    state: ci.CutIntegratorState,
    candidate: cut.CutTransportUnknowns,
    dt_s: float,
    boundary: ct.PoreBoundary,
    **assembly_kwargs: object,
) -> cut.CutTransportAssembly:
    """Assemble one step, naming a forbidden advance instead of dying untyped.

    This wraps the call that survey D2.2 identifies as the outer-face gap:
    ``cut_transport.assemble_backward_euler`` raises
    ``CutTransportTopologyError`` at ``cut_transport.py:1644-1648`` when the
    candidate front reaches or crosses a master face, and
    ``cut_geometry.py:418-419`` raises ``PrimaryDrainageError`` from inside the
    same assembly (``cut_transport.py:1649``) when the candidate front sits
    ABOVE the current one.  Both are FROZEN and are re-raised unchanged unless
    the source state's own assembled fluxes demand an advance, in which case
    the typed refusal is raised FROM them - the original stays as the cause.

    Nothing about the frozen sites changes: this is a catch-site wrapper, which
    design section 2 A2 names as acceptable composition.

    THE SCOPE IS EXACTLY THOSE TWO SITES.  ``assemble_backward_euler``'s own
    outer handler (``cut_transport.py:1592-1600``) re-wraps everything as
    ``CutTransportStepError``, flagging the face-crossing case with
    ``event_restart_required`` and chaining the drainage veto as ``__cause__``.
    Any OTHER rejection - a band violation, a non-finite candidate, a failed
    interface root - is re-raised untouched without paying for a demand
    evaluation and without acquiring a physics label it did not earn.
    """

    try:
        return cut.assemble_backward_euler(
            state.transport,
            candidate,
            dt_s,
            boundary,
            **assembly_kwargs,  # type: ignore[arg-type]
        )
    except cut.CutTransportStepError as exc:
        if not exc.event_restart_required and not isinstance(
            exc.__cause__,
            cg.PrimaryDrainageError,
        ):
            raise
        evidence = evaluate_front_demand(
            state,
            boundary,
            dt_s,
            last_candidate=candidate,
        )
        refuse_on_advance_demand(
            evidence,
            state,
            surfaced_at=(
                "cut_transport.assemble_backward_euler "
                "(cut_transport.py:1644-1649; survey D2.2)"
            ),
            original_refusal=exc,
        )
        raise


# ---------------------------------------------------------------------------
# COMPOSITION POINT (b) - the same-cell / macrostep refusal family
# ---------------------------------------------------------------------------
def advance_with_typed_front_demand(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    same_cell_seed: ci.CutStepSeed,
    **macrostep_kwargs: object,
) -> ceo.CutEventMacrostep:
    """Advance one atomic macrostep, naming a forbidden advance on refusal.

    This wraps ``cut_event_orchestrator.advance_one_interior_face_macrostep``,
    whose typed rejection is the funnel for BOTH families the survey names: the
    ``TOPOLOGY_HANDOFF_FAILED`` path (``cut_event_orchestrator.py:1085-1132``,
    ``:986-993``) and the same-cell / cond-inf saturation family that arrives as
    ``CutIntegratorStepError`` (``cut_event_orchestrator.py:1069-1082``, carried
    on ``CutEventMacrostepError.same_cell_error``).

    An accepted macrostep is returned untouched.  A refused one is re-examined
    on ITS OWN source state's assembled fluxes: an advance demand becomes
    :class:`PrimaryDrainageAdvanceDemanded` chained from the original refusal;
    a recession or stall demand re-raises the original refusal EXACTLY, so no
    existing failure code, message or red test changes meaning.
    """

    try:
        return ceo.advance_one_interior_face_macrostep(
            before,
            dt_s,
            boundary,
            same_cell_seed,
            **macrostep_kwargs,  # type: ignore[arg-type]
        )
    except ceo.CutEventMacrostepError as exc:
        last_candidate = (
            exc.same_cell_error.last_candidate
            if exc.same_cell_error is not None
            else None
        )
        evidence = evaluate_front_demand(
            before,
            boundary,
            dt_s,
            last_candidate=last_candidate,
        )
        refuse_on_advance_demand(
            evidence,
            before,
            surfaced_at=(
                "cut_event_orchestrator.advance_one_interior_face_macrostep "
                f"[{exc.code.value}] (survey D2.2/F1/F8)"
            ),
            original_refusal=exc,
        )
        raise


__all__ = [
    "FROZEN_PHYSICS_AUTHORITY",
    "PHY_014",
    "PHY_023",
    "PHY_036",
    "PRIMARY_DRAINAGE_CITATIONS",
    "SA1_DESIGN_BASIS",
    "SA1_SURVEY_EVIDENCE",
    "FrontDemandError",
    "FrontDemandEvaluationError",
    "FrontDemandEvidence",
    "FrontDemandOutcome",
    "FrontDemandRefusalCode",
    "FrozenPhysicsCitation",
    "PrimaryDrainageAdvanceDemanded",
    "advance_with_typed_front_demand",
    "assemble_with_typed_front_demand",
    "classify_front_demand",
    "evaluate_front_demand",
    "refuse_on_advance_demand",
]
