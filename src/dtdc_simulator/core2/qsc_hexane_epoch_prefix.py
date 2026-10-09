"""BC-2: the pre-birth condensed-hexane epoch on the fully-wet flake.

Design basis: ``docs/GT_PS2_RR2_BIRTH_CONTINUATION_DESIGN_2026-08-21.md``
(D1-D4, and section 3 BC-2 - "prefix integration: activation -> epoch
stepping (chunked, with the film trajectory journaled) -> the typed birth
gate -> the first corrected-physics birth").  Authority: the owner RATIFIED
M1 and M2 on 2026-08-22 (``docs/GT_PS2_MORNING_REVIEW_PACKET_2026-08-22.md``,
Status REVIEWED AND RULED), so TR-1 and TR-2 stand ratified and the design
basis is the ruled basis for this module.

This module orchestrates the epoch.  It owns no constitutive law of its own:
the surface law is the BC-1 leaf
(:mod:`dtdc_simulator.core2.qsc_hexane_appearance_transition_law`), the
interior caloric/conduction advance is the FROZEN wet-core operator
(:mod:`dtdc_simulator.core2.particle.wet_core`), and the dew-clearance
criterion is the FROZEN gas-only chart's own strict admission
(:mod:`dtdc_simulator.core2.particle.coupled_pore`).  Nothing here invents a
criterion, widens a tolerance, or clamps a state.

THE SURFACE-HEAT SEAM (measured at build time, per the packet's instruction)
===========================================================================
The state handed over by stage 1 -
``radial_solvent_tail_engineering_slice.locate_activation_and_lift_wet_state``
- is a :class:`~dtdc_simulator.core2.particle.wet_water.WetWaterState`, and
ITS stepper ``wet_water.step(state, dt_s, boundary)`` takes exactly one
boundary object, ``WetWaterDirichletBoundary(temperature_k,
retained_water_loading, ...)``.  That is a DIRICHLET TRACE and nothing else:
the outer-face energy flux is derived inside ``wet_water._one_face_flux`` from
the prescribed boundary temperature, and there is no additive surface-heat
term, no Neumann/Robin branch and no source vector anywhere in that entry.
The 2N predictor therefore has NO hook for the epoch law's ``energy_to_solid``.

The lawful hook exists one layer down, in the frozen wet-core operator that
``wet_water`` is itself built on::

    wet_core.step(grid_prev, grid, T_prev, dt, inventory, p,
                  surface=("flux", q_surface), source=None, ...)

``surface=("flux", q)`` is a genuine Neumann surface condition
(``wet_core.residual_and_jacobian`` :``489-494``: ``H_conduction =
areas[-1] * q_surface``, OUTWARD positive, with ``dH_conduction = 0``), and
``source`` is a per-cell volumetric energy source.  The composition is EXACT
rather than approximate because ``wet_water.evaluate_cell`` delegates its own
energy density to ``wet_core.energy_density(T, X_h, replace(model.wet,
X_water=W_i, w_o=wo_i))``: under a retained-water field and an oil-label field
that are each UNIFORM across cells, one ``WetCoreParams`` reproduces the 2N
predictor's caloric identity for every cell BIT FOR BIT.  That uniformity is
not assumed - it is a typed admission gate here, and it is exactly the
condition under which ``wet_water``'s own retained-water flux
(``-mobility * d ln a_w(W) / dr``, a function of W alone) is IDENTICALLY zero
at every face.  So during the epoch the water channel does not merely stay
put by fiat; it is stationary by construction, which is also precisely what
mode (i) means - HEXANE condensing, water NONCONDENSING [Tutkun, as printed;
Webb immiscibility].

Under ``surface=("flux", q)`` the frozen discretisation introduces no separate
surface temperature at all: the flux is imposed at the outer face and applied
to the outermost cell.  The epoch law's interface temperature is therefore
unambiguously the outermost wet cell's temperature, which is exactly the
LUMPED-FLAKE convention BC-1 declares in its own docstring (Sipos 1961: on a
lumped meal flake the interface IS the solid surface, so Tutkun's interface
root-find collapses to an explicit evaluation - no inner solver, no
termination floors).  The half-cell conduction trace that a finite-film
boundary would reconstruct is journaled as EVIDENCE ONLY
(``reconstructed_surface_trace_k``), never fed back into the law.

THE EPOCH STEP (D3)
===================
One epoch step over one interval is:

  (i)   evaluate the BC-1 law on the current wet surface state - condensing or
        re-evaporating per its own bidirectional Colburn-Hougen sign, with the
        interface at the outermost wet cell temperature;
  (ii)  advance the interior thermally with ``wet_core.step`` under
        ``surface=("flux", -energy_to_solid_j / (dt * A))`` - the law's own
        energy folded into the surface heat input over the interval;
  (iii) the attached-hexane film accumulator is BC-1's own
        (``hexane_film_accumulator``): exact ``math.fsum`` identity,
        non-negative cone, typed refusals, CC-2a's pattern;
  (iv)  journal T profile, film mass, flux, mode and the gate evidence.

The surface law is evaluated at the START of the interval and the interior is
advanced implicitly over it - a declared first-order explicit-surface /
implicit-interior split.  It is stated, journaled
(``explicit_surface_coupling`` on every record) and never hidden.

THE TYPED BIRTH GATE (D3)
=========================
A birth attempt is admitted only when BOTH hold:

  (i)  DEW CLEARANCE across the would-be newborn piece AND its surface - the
       gas-only chart's OWN strict admission ``coupled_pore.encode_gas_only_y``
       at the boundary composition, which refuses at ``a_h >= 1`` with the
       chart's own message.  No criterion is invented here and no threshold is
       declared: the chart is the authority, evaluated with the existing
       machinery, and the reported activities come from the frozen public
       diagnostic ``coupled_pore.binary_gas_component_activities``.
  (ii) the attached-hexane film has re-evaporated to EXACTLY zero - float
       equality with ``0.0``, which is also the radial activation endpoint's
       own requirement (``sphere.py:152-153``) and the drainage latch's.

Until both hold the epoch simply steps.  A birth attempt below either
condition is a typed refusal, never a clamp and never a nudge.

HOW THE FILM REACHES EXACTLY ZERO
=================================
Condition (ii) is float equality, so the epoch's last step is a LOCATED FILM
EXTINCTION EVENT, in the same discipline ``sphere.locate_activation_event``
uses for the attached film it extinguishes: the inventory is linear in time
over the interval, so the root is analytic; the root is then bracketed to one
binary64 duration; and at the located event the endpoint inventory is
constructed at exactly zero with any sub-resolution remainder REPORTED as
``film_extinction_residual_kg_per_kg_dry``, exactly as ``ActivationEvent``
reports its own balance residuals.  A step that would overdraw the film is a
typed refusal carrying the located duration - the inventory never goes
negative and is never clamped up.  Once the film is gone the epoch is over by
construction: the accumulator's own cone forbids another re-evaporating step.

Every claim flag is false.  This module advances no accepted tray state,
implements no birth solve, and is not physically qualifying.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import InitVar, dataclass, replace
from enum import Enum
from typing import ClassVar

from . import qsc_hexane_appearance_transition_law as bc1
from .particle import coupled_pore as cp
from .particle import surface_active_set as sas
from .particle import wet_core
from .particle import wet_water

QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID = "qsc-hexane-epoch-prefix-step-v1"
#: Revision 2 (BC-6, 2026-08-22): A6 tightened the declared wet-core Newton
#: contract to 1e-13 and R3 added the cumulative energy-closure ledger + gate
#: (two new state fields, journaled and digest-folded).
QSC_HEXANE_EPOCH_PREFIX_SCHEMA_REVISION = 2
QSC_HEXANE_EPOCH_PREFIX_ORCHESTRATION_ID = (
    "pre-birth-condensed-hexane-epoch-on-fully-wet-flake-with-typed-birth-gate-v1"
)
QSC_HEXANE_EPOCH_PREFIX_SOURCE_IDENTITY = "dtdc_simulator.core2.qsc_hexane_epoch_prefix"

#: The epoch's OWN DECLARED wet-core Newton contract for the interior advance.
#: A6 (RULED 2026-08-22, R1 of ``docs/GT_PS2_A4_CUMULATIVE_LEDGER_ANALYSIS_
#: 2026-08-22.md``): the A4 analysis measured that ``wet_core.step``'s frozen
#: default 1e-11 lets the one-Newton-correction early exit leak a same-signed
#: per-step energy defect (drift linear in the one-correction step count, NOT
#: a random walk), which surfaced two layers downstream as the birth's
#: cumulative-ledger refusal.  The owner ratified tightening the epoch's
#: declared contract to 1e-13 - below the measured one-correction plateau,
#: two orders above the measured two-correction floor (~2e-16..2e-15
#: normalized) - which removes the leak at source (measured 1.77e5x reduction
#: at N=4 dt 5e-4).  THE NON-WIDENING PRINCIPLE, stated: declaring TIGHTER
#: than the frozen ``wet_core.step`` default is lawful (the same direction
#: ``coupled_transport.py:203`` permits at construction); declaring WIDER
#: refuses.  ``wet_core.step`` itself is untouched - its other call sites
#: keep the frozen default bit-identically.  This deliberately supersedes the
#: earlier shared-origin declaration; the test suite machine-checks the
#: tightening direction instead of float equality with the default.
QSC_HEXANE_EPOCH_WET_CORE_NEWTON_TOLERANCE = 1.0e-13
#: The frozen wet-core Newton iteration budget, same shared-origin discipline.
QSC_HEXANE_EPOCH_WET_CORE_MAX_ITERATIONS = 30
#: Bounded ULP window for the exact film-exhaustion duration search.  This is
#: not a tolerance: the acceptance test is float EQUALITY with ``0.0``, and the
#: window only bounds how many candidate durations are tried before the search
#: refuses.
QSC_HEXANE_EPOCH_EXHAUSTION_ULP_WINDOW = 4096
#: Bounded halving budget for READING the film rate off the BC-1 law when the
#: caller's own duration already overshoots the film.  The Colburn-Hougen flux
#: does not depend on the film at all - only the accumulator's non-negative
#: cone does - so halving the probe duration recovers the same rate.  This is a
#: bracket budget, not a tolerance: nothing is accepted on it, and exhausting
#: it is a typed refusal.
QSC_HEXANE_EPOCH_EXHAUSTION_PROBE_HALVINGS = 60


class QSCHexaneEpochPrefixRefusalCode(Enum):
    NONFINITE_INPUT = "nonfinite_input"
    NONPOSITIVE_STEP = "nonpositive_step"
    NOT_A_WET_WATER_STATE = "not_a_wet_water_state"
    NONUNIFORM_RETAINED_WATER = "nonuniform_retained_water"
    NONUNIFORM_OIL_LABEL = "nonuniform_oil_label"
    NO_FREE_LIQUID = "no_free_liquid"
    EPOCH_LAW_REFUSED = "epoch_law_refused"
    WET_INTERIOR_ADVANCE_REFUSED = "wet_interior_advance_refused"
    ENERGY_CLOSURE_CONTRACT_VIOLATED = "energy_closure_contract_violated"
    RETAINED_WATER_MOVED = "retained_water_moved"
    FILM_EXHAUSTION_OVERSHOOT = "film_exhaustion_overshoot"
    FILM_EXHAUSTION_NOT_EXACT = "film_exhaustion_not_exact"
    FILM_EXHAUSTION_BEYOND_BOUND = "film_exhaustion_beyond_bound"
    NO_FILM_TO_EXHAUST = "no_film_to_exhaust"
    FILM_NOT_RE_EVAPORATING = "film_not_re_evaporating"
    BIRTH_GATE_FILM_PRESENT = "birth_gate_film_present"
    BIRTH_GATE_DEW_CLEARANCE_ABSENT = "birth_gate_dew_clearance_absent"
    CHART_EVALUATION_REFUSED = "chart_evaluation_refused"
    #: R3 (RULED 2026-08-22, folded into BC-6): the epoch's own cumulative
    #: energy-closure ledger exceeded its construction-derived accumulated
    #: bound - the drift is refused HERE, never handed across the seam.
    CUMULATIVE_ENERGY_CLOSURE_BEYOND_BOUND = (
        "cumulative_energy_closure_beyond_construction_bound"
    )


class QSCHexaneEpochPrefixError(ValueError):
    """Typed refusal for the pre-birth condensed-hexane epoch orchestration.

    The refusal carries the epoch state it refused ON, so the caller's rollback
    is exact: no epoch entry point ever mutates its input.
    """

    def __init__(
        self,
        code: QSCHexaneEpochPrefixRefusalCode,
        message: str,
        *,
        state: object | None = None,
        evidence: dict | None = None,
    ) -> None:
        super().__init__(f"{code.value}: {message}")
        self.code = code
        self.state = state
        self.evidence = dict(evidence) if evidence else {}


_EPOCH_SEAL = object()


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
    if type(value) is bool:
        return b"b1:" + (b"1" if value else b"0")
    raise QSCHexaneEpochPrefixError(
        QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
        "digest folding admits only str, int, bool, and float leaves",
    )


def _digest(label: str, values: tuple[object, ...]) -> str:
    hasher = hashlib.sha256()
    hasher.update(label.encode("utf-8"))
    for value in values:
        hasher.update(_frame(value))
    return "sha256:" + hasher.hexdigest()


def _require_finite(label: str, value: object) -> float:
    if type(value) is not float or not math.isfinite(value):
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
            f"{label} must be a finite binary64 value",
        )
    return value


class _FalseEpochClaims:
    pre_birth_condensed_hexane_epoch_only: ClassVar[bool] = True
    fully_wet_domain_only_phy_036: ClassVar[bool] = True
    explicit_surface_implicit_interior_split: ClassVar[bool] = True
    lumped_flake_interface_convention: ClassVar[bool] = True
    surface_heat_enters_through_frozen_wet_core_flux_bc: ClassVar[bool] = True
    owner_ratified_design_basis: ClassVar[bool] = True
    birth_solve_implemented: ClassVar[bool] = False
    eventual_birth_seam_d4_implemented: ClassVar[bool] = False
    advances_accepted_tray_state: ClassVar[bool] = False
    imported_by_production: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False


# ---------------------------------------------------------------------------
# The declared epoch boundary
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneEpochBoundary:
    """The caller's declared bulk gas state and film coefficients.

    The coefficients are DECLARED with provenance and are never fitted here.
    They are the same pair of quantities BC-1's layer input consumes; this
    object simply carries them once for a whole epoch instead of once per step.
    """

    bulk_temperature_k: float
    bulk_hexane_mole_fraction: float
    pressure_pa: float
    molar_transfer_coefficient_mol_m2_s: float
    heat_transfer_coefficient_w_m2_k: float
    coefficient_provenance: str
    label: str = "declared pre-birth epoch boundary"
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for name in ("coefficient_provenance", "label"):
            value = getattr(self, name)
            if type(value) is not str or not value.strip():
                raise QSCHexaneEpochPrefixError(
                    QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
                    f"epoch boundary {name} must be nonblank",
                )
        for label, value in (
            ("bulk temperature", self.bulk_temperature_k),
            ("bulk hexane mole fraction", self.bulk_hexane_mole_fraction),
            ("pressure", self.pressure_pa),
            ("molar transfer coefficient", self.molar_transfer_coefficient_mol_m2_s),
            ("heat transfer coefficient", self.heat_transfer_coefficient_w_m2_k),
        ):
            _require_finite(label, value)
        if not 0.0 < self.bulk_hexane_mole_fraction < 1.0:
            raise QSCHexaneEpochPrefixError(
                QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
                "the bulk hexane mole fraction must lie strictly inside the simplex",
            )
        if self.pressure_pa <= 0.0:
            raise QSCHexaneEpochPrefixError(
                QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
                "the epoch boundary pressure must be strictly positive",
            )
        if (
            self.molar_transfer_coefficient_mol_m2_s <= 0.0
            or self.heat_transfer_coefficient_w_m2_k <= 0.0
        ):
            raise QSCHexaneEpochPrefixError(
                QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
                "declared film coefficients must be strictly positive",
            )


# ---------------------------------------------------------------------------
# The epoch state
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneEpochState(_FalseEpochClaims):
    """One sealed pre-birth epoch state: a fully-wet flake plus its film."""

    wet_state: wet_water.WetWaterState
    surface_inventories: sas.ExternalSurfaceInventories
    epoch_time_s: float
    step_index: int
    cumulative_energy_to_solid_j: float
    cumulative_condensed_hexane_kg: float
    cumulative_re_evaporated_hexane_kg: float
    #: R3 (RULED 2026-08-22): the epoch's OWN cumulative energy-closure
    #: ledger - the exactly-rounded running sum of every accepted step's
    #: ``energy_closure_residual_j`` (== the wet-core residual sum by the
    #: proved telescoping identity) - beside its construction-derived bound,
    #: the accumulated ``n * tol * scale`` of the declared Newton contract
    #: (the epoch suite's own accumulated-bound construction).  The commit
    #: gate refuses when the ledger leaves its bound, so a long prefix can
    #: never hand a silently-drifted state across the birth seam.
    cumulative_energy_closure_residual_j: float
    cumulative_energy_closure_bound_j: float
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID
    schema_revision: ClassVar[int] = QSC_HEXANE_EPOCH_PREFIX_SCHEMA_REVISION
    orchestration_id: ClassVar[str] = QSC_HEXANE_EPOCH_PREFIX_ORCHESTRATION_ID

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _EPOCH_SEAL:
            raise TypeError("pre-birth epoch states are orchestration-sealed")

    # -- derived, all from the wet state's own frozen primitives -----------
    @property
    def grid(self):
        return self.wet_state.grid

    @property
    def particle_radius_m(self) -> float:
        return self.wet_state.grid.R

    @property
    def surface_area_m2(self) -> float:
        return self.wet_state.grid.areas[-1]

    @property
    def basis_kg_dry(self) -> float:
        return math.fsum(
            volume * self.wet_state.model.wet.rho_dm_p
            for volume in self.wet_state.grid.volumes
        )

    @property
    def interface_temperature_k(self) -> float:
        """The outermost wet cell temperature - BC-1's lumped-flake interface."""

        return self.wet_state.temperatures_k[-1]

    @property
    def film_kg_per_kg_dry(self) -> float:
        return self.surface_inventories.attached_hexane_kg_per_kg_dry

    @property
    def film_kg(self) -> float:
        return self.film_kg_per_kg_dry * self.basis_kg_dry

    @property
    def cell_params(self) -> wet_core.WetCoreParams:
        """The one ``WetCoreParams`` that reproduces the 2N caloric identity.

        Admissible only on the uniform retained-water / oil-label field the
        epoch gate enforces; see the module docstring's seam section.
        """

        return replace(
            self.wet_state.model.wet,
            X_water=self.wet_state.retained_water_loadings[0],
            w_o=self.wet_state.oil_fraction_labels[0],
        )

    @property
    def free_liquid_loading_kg_per_kg_dry(self) -> float:
        """The outermost cell's mobile (free) n-hexane loading, PHY-036's X_f."""

        return wet_core.partition(
            self.interface_temperature_k,
            self.wet_state.historical_hexane_loadings[-1],
            self.cell_params,
        ).mobile_liquid

    @property
    def wet_internal_energy_j(self) -> float:
        params = self.cell_params
        return math.fsum(
            volume * wet_core.energy_density(temperature, hexane, params)
            for volume, temperature, hexane in zip(
                self.wet_state.grid.volumes,
                self.wet_state.temperatures_k,
                self.wet_state.historical_hexane_loadings,
                strict=True,
            )
        )

    @property
    def definition_digest(self) -> str:
        leaves: list[object] = [
            QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID,
            QSC_HEXANE_EPOCH_PREFIX_SCHEMA_REVISION,
            QSC_HEXANE_EPOCH_PREFIX_ORCHESTRATION_ID,
            self.epoch_time_s,
            self.step_index,
            self.surface_inventories.attached_hexane_kg_per_kg_dry,
            self.surface_inventories.free_water_kg_per_kg_dry,
            self.cumulative_energy_to_solid_j,
            self.cumulative_condensed_hexane_kg,
            self.cumulative_re_evaporated_hexane_kg,
            self.cumulative_energy_closure_residual_j,
            self.cumulative_energy_closure_bound_j,
            self.wet_state.grid.R,
            self.wet_state.time_s,
        ]
        leaves.extend(self.wet_state.temperatures_k)
        leaves.extend(self.wet_state.retained_water_loadings)
        leaves.extend(self.wet_state.historical_hexane_loadings)
        return _digest("qsc-hexane-epoch-prefix-state", tuple(leaves))

    def journal_record(self) -> dict:
        """A JSON-ready line for the epoch journal (design BC-2 (iv))."""

        return {
            "schema_id": QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID,
            "step_index": self.step_index,
            "epoch_time_s": self.epoch_time_s,
            "wet_time_s": self.wet_state.time_s,
            "wet_temperatures_k": list(self.wet_state.temperatures_k),
            "interface_temperature_k": self.interface_temperature_k,
            "retained_water_loadings": list(self.wet_state.retained_water_loadings),
            "historical_hexane_loadings": list(self.wet_state.historical_hexane_loadings),
            "film_kg_per_kg_dry": self.film_kg_per_kg_dry,
            "film_kg": self.film_kg,
            "free_water_kg_per_kg_dry": (
                self.surface_inventories.free_water_kg_per_kg_dry
            ),
            "surface_topology": self.surface_inventories.topology.value,
            "cumulative_energy_to_solid_j": self.cumulative_energy_to_solid_j,
            "cumulative_condensed_hexane_kg": self.cumulative_condensed_hexane_kg,
            "cumulative_re_evaporated_hexane_kg": (
                self.cumulative_re_evaporated_hexane_kg
            ),
            "cumulative_energy_closure_residual_j": (
                self.cumulative_energy_closure_residual_j
            ),
            "cumulative_energy_closure_bound_j": (
                self.cumulative_energy_closure_bound_j
            ),
            "wet_internal_energy_j": self.wet_internal_energy_j,
            "state_digest": self.definition_digest,
            "physically_qualifying": False,
        }


def _require_epoch_admissible_wet_state(state: wet_water.WetWaterState) -> None:
    """Fail closed unless the epoch's exact composition preconditions hold."""

    if type(state) is not wet_water.WetWaterState:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NOT_A_WET_WATER_STATE,
            "the epoch consumes the exact frozen WetWaterState handed over by "
            "the activation/lift stage; no other payload is admitted",
        )
    water = state.retained_water_loadings
    oil = state.oil_fraction_labels
    if len(set(water)) != 1:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONUNIFORM_RETAINED_WATER,
            "the epoch's exact composition with the frozen wet-core operator "
            f"needs ONE retained-water loading across cells, got {water!r}: a "
            "nonuniform field would both break the bit-exact caloric identity "
            "and drive a nonzero wet-water face flux, which mode (i) "
            "(water NONCONDENSING) forbids; fail closed, never average",
        )
    if len(set(oil)) != 1:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONUNIFORM_OIL_LABEL,
            f"the epoch needs ONE oil-fraction label across cells, got {oil!r}; "
            "fail closed, never average",
        )


def initialize_qsc_hexane_epoch(
    wet_state: wet_water.WetWaterState,
    *,
    attached_hexane_kg_per_kg_dry: float = 0.0,
    free_water_kg_per_kg_dry: float = 0.0,
    epoch_time_s: float = 0.0,
) -> QSCHexaneEpochState:
    """Enter the pre-birth epoch on the activation/lift hand-off state."""

    _require_epoch_admissible_wet_state(wet_state)
    for label, value in (
        ("attached hexane film", attached_hexane_kg_per_kg_dry),
        ("external free water", free_water_kg_per_kg_dry),
        ("epoch entry time", epoch_time_s),
    ):
        _require_finite(label, value)
    if epoch_time_s < 0.0:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
            "the epoch entry time must be non-negative",
        )
    try:
        inventories = sas.ExternalSurfaceInventories(
            attached_hexane_kg_per_kg_dry, free_water_kg_per_kg_dry
        )
    except ValueError as exc:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
            f"the external surface inventories refused their own cone: {exc}",
        ) from exc
    state = QSCHexaneEpochState(
        wet_state=wet_state,
        surface_inventories=inventories,
        epoch_time_s=epoch_time_s,
        step_index=0,
        cumulative_energy_to_solid_j=0.0,
        cumulative_condensed_hexane_kg=0.0,
        cumulative_re_evaporated_hexane_kg=0.0,
        cumulative_energy_closure_residual_j=0.0,
        cumulative_energy_closure_bound_j=0.0,
        _seal=_EPOCH_SEAL,
    )
    # PHY-036: X_f > 0 <=> s = R.  A vanished free-liquid loading is not the
    # fully-wet state and the epoch law refuses it; say so at entry rather than
    # one step later.
    if state.free_liquid_loading_kg_per_kg_dry <= 0.0:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NO_FREE_LIQUID,
            "PHY-036 reads X_f > 0 <=> s = R: the hand-off state carries no "
            "mobile n-hexane at its surface, so it is not the fully-wet state "
            "the epoch law is admissible on",
            state=state,
        )
    return state


# ---------------------------------------------------------------------------
# One epoch step
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedQSCHexaneEpochStep(_FalseEpochClaims):
    """One accepted epoch step: the BC-1 law folded into the frozen interior."""

    before: QSCHexaneEpochState
    after: QSCHexaneEpochState
    law_step: bc1.ValidatedQSCHexaneAppearanceTransitionStep
    mode: bc1.QSCHexaneEpochMode
    duration_s: float
    #: Outward-positive, exactly the sign convention ``wet_core`` uses.
    surface_heat_flux_w_m2: float
    energy_to_solid_j: float
    wet_energy_gain_j: float
    energy_closure_residual_j: float
    normalized_wet_core_residual: float
    reconstructed_surface_trace_k: float
    film_before_kg_per_kg_dry: float
    film_after_kg_per_kg_dry: float
    film_deposit_kg: float
    #: Strictly zero except on a located FILM EXTINCTION EVENT, where it is the
    #: sub-resolution remainder the endpoint construction retired, reported as a
    #: balance residual exactly as ``sphere.ActivationEvent`` reports its own.
    film_extinction_residual_kg_per_kg_dry: float
    signed_molar_flux_mol_m2_s: float
    interface_hexane_mole_fraction: float
    exhaustion_step: bool
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID
    schema_revision: ClassVar[int] = QSC_HEXANE_EPOCH_PREFIX_SCHEMA_REVISION
    orchestration_id: ClassVar[str] = QSC_HEXANE_EPOCH_PREFIX_ORCHESTRATION_ID
    explicit_surface_coupling: ClassVar[bool] = True

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _EPOCH_SEAL:
            raise TypeError("pre-birth epoch steps are orchestration-sealed")

    def journal_record(self) -> dict:
        record = self.after.journal_record()
        record.update(
            {
                "mode": self.mode.value,
                "dt_s": self.duration_s,
                "surface_heat_flux_w_m2": self.surface_heat_flux_w_m2,
                "energy_to_solid_j": self.energy_to_solid_j,
                "wet_energy_gain_j": self.wet_energy_gain_j,
                "energy_closure_residual_j": self.energy_closure_residual_j,
                "normalized_wet_core_residual": self.normalized_wet_core_residual,
                "reconstructed_surface_trace_k": self.reconstructed_surface_trace_k,
                "film_before_kg_per_kg_dry": self.film_before_kg_per_kg_dry,
                "film_after_kg_per_kg_dry": self.film_after_kg_per_kg_dry,
                "film_deposit_kg": self.film_deposit_kg,
                "film_extinction_residual_kg_per_kg_dry": (
                    self.film_extinction_residual_kg_per_kg_dry
                ),
                "signed_molar_flux_mol_m2_s": self.signed_molar_flux_mol_m2_s,
                "interface_hexane_mole_fraction": self.interface_hexane_mole_fraction,
                "interface_temperature_before_k": self.before.interface_temperature_k,
                "wet_temperatures_before_k": list(self.before.wet_state.temperatures_k),
                "exhaustion_step": self.exhaustion_step,
                "explicit_surface_coupling": True,
                "lumped_flake_interface_convention": True,
            }
        )
        return record


def _law_step(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    duration_s: float,
    declared_mode: bc1.QSCHexaneEpochMode | None,
) -> bc1.ValidatedQSCHexaneAppearanceTransitionStep:
    """Evaluate the BC-1 leaf on the current wet surface state."""

    radius = state.particle_radius_m
    assertion = bc1.QSCHexaneFullyWetAssertion(
        particle_radius_m=radius,
        # PHY-036: while free liquid is present the drainage front SITS AT the
        # particle surface.  This is not a nudge - the fully-wet state has no
        # other front position.
        front_position_m=radius,
        free_liquid_loading_kg_per_kg_dry=state.free_liquid_loading_kg_per_kg_dry,
        drainage_history=sas.HexanePrimaryDrainageHistory.PRE_CORE_RECESSION,
        surface_inventories=state.surface_inventories,
        basis_kg_dry=state.basis_kg_dry,
    )
    layer = bc1.QSCHexaneAppearanceLayerInput(
        group_id=f"epoch-surface-step-{state.step_index}",
        bulk_temperature_k=boundary.bulk_temperature_k,
        bulk_hexane_mole_fraction=boundary.bulk_hexane_mole_fraction,
        pressure_pa=boundary.pressure_pa,
        particle_temperature_k=state.interface_temperature_k,
        molar_transfer_coefficient_mol_m2_s=(
            boundary.molar_transfer_coefficient_mol_m2_s
        ),
        heat_transfer_coefficient_w_m2_k=boundary.heat_transfer_coefficient_w_m2_k,
        active_area_m2=state.surface_area_m2,
        coefficient_provenance=boundary.coefficient_provenance,
    )
    request = bc1.QSCHexaneAppearanceTransitionRequest(
        attempt_id=f"epoch-step-{state.step_index}",
        fully_wet_state=assertion,
        layers=(layer,),
        step_duration_s=duration_s,
        declared_mode=declared_mode,
    )
    return bc1.evaluate_qsc_hexane_appearance_transition_step(request)


def _newton_scale(state: QSCHexaneEpochState) -> float:
    """``wet_core.step``'s own capacity-based Newton scale, restated exactly."""

    params = state.cell_params
    grid = state.wet_state.grid
    return max(
        max(
            grid.volumes[index]
            * wet_core.energy_capacity(
                state.wet_state.temperatures_k[index],
                state.wet_state.historical_hexane_loadings[index],
                params,
            )
            * max(abs(state.wet_state.temperatures_k[index]), 1.0),
            1.0e-12,
        )
        for index in range(grid.n)
    )


def _advance_interior(
    state: QSCHexaneEpochState,
    duration_s: float,
    surface_heat_flux_w_m2: float,
) -> tuple[tuple[float, ...], float, float]:
    """Advance the interior implicitly under the frozen wet-core flux BC.

    Returns ``(temperatures, surface_rate_w, normalized_residual)``.
    """

    grid = state.wet_state.grid
    params = state.cell_params
    hexane = list(state.wet_state.historical_hexane_loadings)
    previous = list(state.wet_state.temperatures_k)
    surface = ("flux", surface_heat_flux_w_m2)
    try:
        temperatures = wet_core.step(
            grid,
            grid,
            previous,
            duration_s,
            hexane,
            params,
            surface,
            None,
            None,
            tol=QSC_HEXANE_EPOCH_WET_CORE_NEWTON_TOLERANCE,
            max_iter=QSC_HEXANE_EPOCH_WET_CORE_MAX_ITERATIONS,
        )
    except Exception as exc:  # noqa: BLE001 - every frozen refusal is the datum
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.WET_INTERIOR_ADVANCE_REFUSED,
            f"the frozen wet-core interior advance refused with exact rollback: {exc}",
            state=state,
            evidence={
                "duration_s": duration_s,
                "surface_heat_flux_w_m2": surface_heat_flux_w_m2,
                "temperatures_before_k": list(state.wet_state.temperatures_k),
            },
        ) from exc
    residual, _sub, _diag, _sup, surface_rate = wet_core.residual_and_jacobian(
        grid,
        grid,
        previous,
        temperatures,
        duration_s,
        hexane,
        params,
        surface,
        None,
        None,
    )
    normalized = max(abs(value) for value in residual) / _newton_scale(state)
    return tuple(temperatures), surface_rate, normalized


def _commit(
    state: QSCHexaneEpochState,
    law: bc1.ValidatedQSCHexaneAppearanceTransitionStep,
    duration_s: float,
    *,
    exhaustion_step: bool,
    film_extinction_residual_kg_per_kg_dry: float = 0.0,
) -> ValidatedQSCHexaneEpochStep:
    """Fold one law step into the interior and seal the accepted epoch step.

    ``film_extinction_residual_kg_per_kg_dry`` is nonzero only on a located
    FILM EXTINCTION EVENT (see :func:`advance_to_film_exhaustion`), where the
    endpoint inventory is constructed at exactly zero and the sub-resolution
    remainder is REPORTED rather than hidden - the same discipline
    ``sphere.locate_activation_event`` applies to the attached film it
    extinguishes ("epsilon films are forbidden"), with its own balance residual
    carried as evidence.
    """

    area = state.surface_area_m2
    energy_to_solid = law.total_energy_to_solid_j
    # wet_core's surface flux is OUTWARD positive; the law's energy_to_solid is
    # INWARD positive.  One negation, stated once, is the whole coupling.
    surface_heat_flux = -energy_to_solid / duration_s / area
    temperatures, surface_rate, normalized = _advance_interior(
        state, duration_s, surface_heat_flux
    )
    if normalized > QSC_HEXANE_EPOCH_WET_CORE_NEWTON_TOLERANCE:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.ENERGY_CLOSURE_CONTRACT_VIOLATED,
            "the interior advance did not meet the FROZEN wet-core Newton "
            f"contract ({QSC_HEXANE_EPOCH_WET_CORE_NEWTON_TOLERANCE}): measured "
            f"{normalized}; the epoch declares no tolerance of its own and "
            "refuses rather than accept a looser closure",
            state=state,
            evidence={"normalized_wet_core_residual": normalized},
        )
    params = state.cell_params
    grid = state.wet_state.grid
    gain = math.fsum(
        volume
        * (
            wet_core.energy_density(after_t, hexane, params)
            - wet_core.energy_density(before_t, hexane, params)
        )
        for volume, before_t, after_t, hexane in zip(
            grid.volumes,
            state.wet_state.temperatures_k,
            temperatures,
            state.wet_state.historical_hexane_loadings,
            strict=True,
        )
    )
    # R3 (RULED 2026-08-22): the epoch's own cumulative energy-closure ledger
    # and gate.  The per-step closure is the wet-core residual SUM (proved
    # telescoping identity); the declared per-step contract
    # ``max_i |r_i| / scale <= tol`` bounds it by ``n * tol * scale`` with the
    # operator's own capacity scale, so the accumulated sum of those bounds is
    # the contract's OWN bound on the cumulative closure - the epoch suite's
    # accumulated-bound construction, now enforced at commit rather than only
    # asserted in tests.  A ledger outside its bound refuses HERE with exact
    # rollback; downstream layers never inherit an ungated drift.
    closure_residual = math.fsum((gain, -energy_to_solid))
    step_closure_bound = (
        grid.n * QSC_HEXANE_EPOCH_WET_CORE_NEWTON_TOLERANCE * _newton_scale(state)
    )
    cumulative_closure = math.fsum(
        (state.cumulative_energy_closure_residual_j, closure_residual)
    )
    cumulative_closure_bound = math.fsum(
        (state.cumulative_energy_closure_bound_j, step_closure_bound)
    )
    if abs(cumulative_closure) > cumulative_closure_bound:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.CUMULATIVE_ENERGY_CLOSURE_BEYOND_BOUND,
            "the epoch's cumulative energy-closure ledger left its own "
            f"construction-derived bound: |{cumulative_closure}| J > "
            f"{cumulative_closure_bound} J accumulated over "
            f"{state.step_index + 1} steps; the epoch refuses rather than hand "
            "a drifted state across the birth seam",
            state=state,
            evidence={
                "cumulative_energy_closure_residual_j": cumulative_closure,
                "cumulative_energy_closure_bound_j": cumulative_closure_bound,
                "step_energy_closure_residual_j": closure_residual,
                "step_energy_closure_bound_j": step_closure_bound,
            },
        )
    boundary_energy_out = duration_s * surface_rate
    advanced_wet = replace(
        state.wet_state,
        time_s=state.wet_state.time_s + duration_s,
        temperatures_k=temperatures,
        cumulative_boundary_energy_out_j=math.fsum(
            (state.wet_state.cumulative_boundary_energy_out_j, boundary_energy_out)
        ),
        cumulative_absolute_energy_exchange_j=math.fsum(
            (
                state.wet_state.cumulative_absolute_energy_exchange_j,
                abs(boundary_energy_out),
            )
        ),
    )
    # Mode (i): water is the NONCONDENSING partner, so the retained-water field
    # must be untouched.  Assert it rather than trust it.
    if advanced_wet.retained_water_loadings != state.wet_state.retained_water_loadings:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.RETAINED_WATER_MOVED,
            "an epoch step moved the retained water; mode (i) carries water as "
            "the noncondensing species and the epoch never transports it",
            state=state,
        )
    if film_extinction_residual_kg_per_kg_dry == 0.0:
        inventories_after = law.film.after
    else:
        # THE LOCATED FILM EXTINCTION EVENT.  The accumulator's own endpoint
        # still holds a strictly positive sub-resolution remainder, and the
        # located root lies inside the next binary64 duration, which overdraws.
        # The house forbids epsilon films at an extinction endpoint
        # (``sphere.activate_radially``: "radial activation requires exactly
        # zero attached n-hexane"), so the endpoint is CONSTRUCTED at exactly
        # zero and the remainder is carried as a reported balance residual -
        # never silently absorbed, and never a negative inventory clamped up.
        inventories_after = sas.ExternalSurfaceInventories(
            0.0, law.film.after.free_water_kg_per_kg_dry
        )
    realized_deposit_kg = (
        inventories_after.attached_hexane_kg_per_kg_dry - state.film_kg_per_kg_dry
    ) * state.basis_kg_dry
    condensed = state.cumulative_condensed_hexane_kg
    re_evaporated = state.cumulative_re_evaporated_hexane_kg
    if law.mode is bc1.QSCHexaneEpochMode.CONDENSING:
        condensed = math.fsum((condensed, realized_deposit_kg))
    else:
        re_evaporated = math.fsum((re_evaporated, -realized_deposit_kg))
    after = QSCHexaneEpochState(
        wet_state=advanced_wet,
        surface_inventories=inventories_after,
        epoch_time_s=math.fsum((state.epoch_time_s, duration_s)),
        step_index=state.step_index + 1,
        cumulative_energy_to_solid_j=math.fsum(
            (state.cumulative_energy_to_solid_j, energy_to_solid)
        ),
        cumulative_condensed_hexane_kg=condensed,
        cumulative_re_evaporated_hexane_kg=re_evaporated,
        cumulative_energy_closure_residual_j=cumulative_closure,
        cumulative_energy_closure_bound_j=cumulative_closure_bound,
        _seal=_EPOCH_SEAL,
    )
    # Evidence only: the half-cell conduction trace a finite-film boundary
    # would reconstruct.  It is never fed back into the law - the epoch runs on
    # BC-1's declared lumped-flake interface (the outer cell temperature).
    half_cell = grid.R - grid.centers[-1]
    trace = state.interface_temperature_k - (
        surface_heat_flux * half_cell / params.conductivity
    )
    result = law.layer_results[0]
    return ValidatedQSCHexaneEpochStep(
        before=state,
        after=after,
        law_step=law,
        mode=law.mode,
        duration_s=duration_s,
        surface_heat_flux_w_m2=surface_heat_flux,
        energy_to_solid_j=energy_to_solid,
        wet_energy_gain_j=gain,
        energy_closure_residual_j=closure_residual,
        normalized_wet_core_residual=normalized,
        reconstructed_surface_trace_k=trace,
        film_before_kg_per_kg_dry=state.film_kg_per_kg_dry,
        film_after_kg_per_kg_dry=after.film_kg_per_kg_dry,
        film_deposit_kg=realized_deposit_kg,
        film_extinction_residual_kg_per_kg_dry=film_extinction_residual_kg_per_kg_dry,
        signed_molar_flux_mol_m2_s=result.signed_molar_flux_mol_m2_s,
        interface_hexane_mole_fraction=result.interface_hexane_mole_fraction,
        exhaustion_step=exhaustion_step,
        _seal=_EPOCH_SEAL,
    )


def advance_qsc_hexane_epoch_step(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    duration_s: float,
    *,
    declared_mode: bc1.QSCHexaneEpochMode | None = None,
) -> ValidatedQSCHexaneEpochStep:
    """Advance one pre-birth epoch step; refuse anything off the design basis."""

    if type(state) is not QSCHexaneEpochState:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NOT_A_WET_WATER_STATE,
            "the epoch advance consumes the exact sealed epoch state",
        )
    if type(boundary) is not QSCHexaneEpochBoundary:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
            "the epoch advance consumes the exact typed epoch boundary",
        )
    _require_finite("epoch step duration", duration_s)
    if duration_s <= 0.0:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONPOSITIVE_STEP,
            "the epoch step duration must be strictly positive",
            state=state,
        )
    _require_epoch_admissible_wet_state(state.wet_state)
    try:
        law = _law_step(state, boundary, duration_s, declared_mode)
    except bc1.QSCHexaneAppearanceTransitionError as exc:
        if exc.code is bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR:
            evidence = _exhaustion_evidence(state, boundary, duration_s)
            raise QSCHexaneEpochPrefixError(
                QSCHexaneEpochPrefixRefusalCode.FILM_EXHAUSTION_OVERSHOOT,
                "this interval would re-evaporate more hexane than the film "
                f"holds ({exc}); the epoch refuses rather than clamp - take the "
                "exact exhaustion step instead (advance_to_film_exhaustion), "
                f"whose duration at this state is {evidence.get('exhaustion_duration_s')!r} s",
                state=state,
                evidence=evidence,
            ) from exc
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.EPOCH_LAW_REFUSED,
            f"the BC-1 epoch law refused this step: {exc}",
            state=state,
            evidence={"law_refusal_code": exc.code.value, "duration_s": duration_s},
        ) from exc
    return _commit(state, law, duration_s, exhaustion_step=False)


def _probe_film_rate(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    maximum_duration_s: float,
) -> tuple[bc1.ValidatedQSCHexaneAppearanceTransitionStep, float]:
    """Read the film rate off BC-1 through a bounded duration bracket.

    The Colburn-Hougen flux is a function of the surface/bulk state alone, so
    the intensive film rate does not depend on the step duration; only the
    accumulator's non-negative cone does.  When the caller's own duration
    already overshoots the film, halving it recovers the SAME rate without
    clamping anything.  Exhausting the bracket budget is a typed refusal.
    """

    duration = maximum_duration_s
    for _ in range(QSC_HEXANE_EPOCH_EXHAUSTION_PROBE_HALVINGS):
        if not duration > 0.0:
            break
        try:
            return _law_step(state, boundary, duration, None), duration
        except bc1.QSCHexaneAppearanceTransitionError as exc:
            if exc.code is not bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR:
                raise
            # The probe overshot, so the exhaustion time is shorter than this.
            duration = duration / 2.0
    raise QSCHexaneEpochPrefixError(
        QSCHexaneEpochPrefixRefusalCode.FILM_EXHAUSTION_NOT_EXACT,
        "the film rate could not be read within "
        f"{QSC_HEXANE_EPOCH_EXHAUSTION_PROBE_HALVINGS} halvings of "
        f"{maximum_duration_s} s; refuse rather than guess a rate",
        state=state,
    )


def _exhaustion_evidence(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    probe_duration_s: float,
) -> dict:
    """Report the analytic exhaustion duration without taking any step."""

    film = state.film_kg_per_kg_dry
    try:
        probe, _ = _probe_film_rate(state, boundary, probe_duration_s)
    except (QSCHexaneEpochPrefixError, bc1.QSCHexaneAppearanceTransitionError):
        return {"film_kg_per_kg_dry": film, "exhaustion_duration_s": None}
    rate = probe.film.rate_kg_per_kg_dry_s
    if rate >= 0.0:
        return {
            "film_kg_per_kg_dry": film,
            "film_rate_kg_per_kg_dry_s": rate,
            "exhaustion_duration_s": None,
        }
    return {
        "film_kg_per_kg_dry": film,
        "film_rate_kg_per_kg_dry_s": rate,
        "exhaustion_duration_s": film / -rate,
    }


def advance_to_film_exhaustion(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    *,
    maximum_duration_s: float,
) -> ValidatedQSCHexaneEpochStep:
    """Take the one epoch step that lands the film on EXACTLY zero.

    The birth gate's second condition is float equality with ``0.0`` (D3, and
    the radial activation endpoint's own requirement), so the epoch needs a step
    that reaches it exactly rather than approximately.  The event is located the
    way ``sphere.locate_activation_event`` locates its own attached-film
    extinction: the inventory is linear in time over the interval, so the root
    is analytic, ``t_e = X_film / |dX_film/dt|``.

    The located root is then BRACKETED IN BINARY64: the routine finds the
    LARGEST duration whose BC-1 accumulator still leaves a non-negative film,
    so the very next representable duration overdraws and the true root lies
    strictly inside that one-ULP interval.  Two outcomes, both exact:

    * the bracketing duration's own accumulator lands on ``0.0`` - an ordinary
      accumulator step, nothing is constructed;
    * it leaves a strictly positive SUB-RESOLUTION remainder - the located
      FILM EXTINCTION EVENT.  The endpoint inventory is constructed at exactly
      zero (the house forbids epsilon films at an extinction endpoint) and the
      remainder is reported as
      ``film_extinction_residual_kg_per_kg_dry``, exactly as
      ``sphere.ActivationEvent`` reports its own balance residuals.

    A negative inventory is never produced and never clamped up; the search
    only ever walks the non-negative side of the bracket.
    """

    if type(state) is not QSCHexaneEpochState:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NOT_A_WET_WATER_STATE,
            "the exhaustion step consumes the exact sealed epoch state",
        )
    _require_finite("maximum exhaustion duration", maximum_duration_s)
    if maximum_duration_s <= 0.0:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONPOSITIVE_STEP,
            "the exhaustion duration bound must be strictly positive",
            state=state,
        )
    film = state.film_kg_per_kg_dry
    if film == 0.0:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NO_FILM_TO_EXHAUST,
            "the attached-hexane film is already exactly zero; there is no "
            "exhaustion step to take",
            state=state,
        )
    probe, _probe_duration = _probe_film_rate(state, boundary, maximum_duration_s)
    if probe.mode is not bc1.QSCHexaneEpochMode.RE_EVAPORATING:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.FILM_NOT_RE_EVAPORATING,
            "the Colburn-Hougen driving force is CONDENSING at this state: the "
            "film is growing, not exhausting; the epoch's direction is evidence, "
            "not a knob",
            state=state,
            evidence={"mode": probe.mode.value},
        )
    rate = probe.film.rate_kg_per_kg_dry_s
    analytic = film / -rate
    if analytic > maximum_duration_s:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.FILM_EXHAUSTION_BEYOND_BOUND,
            f"the analytic exhaustion duration {analytic} s exceeds the declared "
            f"bound {maximum_duration_s} s; step the epoch further rather than "
            "widen the bound implicitly",
            state=state,
            evidence={
                "exhaustion_duration_s": analytic,
                "maximum_duration_s": maximum_duration_s,
                "film_kg_per_kg_dry": film,
                "film_rate_kg_per_kg_dry_s": rate,
            },
        )
    bracket = _bracket_extinction_duration(state, boundary, analytic)
    if bracket is None:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.FILM_EXHAUSTION_NOT_EXACT,
            "the film-extinction root could not be bracketed in binary64 within "
            f"{QSC_HEXANE_EPOCH_EXHAUSTION_ULP_WINDOW} durations of {analytic} s; "
            "refuse rather than construct an endpoint on an unlocated root",
            state=state,
            evidence={
                "exhaustion_duration_s": analytic,
                "film_kg_per_kg_dry": film,
                "film_rate_kg_per_kg_dry_s": rate,
            },
        )
    duration, law = bracket
    residue = law.film.after.attached_hexane_kg_per_kg_dry
    return _commit(
        state,
        law,
        duration,
        exhaustion_step=True,
        film_extinction_residual_kg_per_kg_dry=residue,
    )


def _bracket_extinction_duration(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    analytic: float,
) -> tuple[float, bc1.ValidatedQSCHexaneAppearanceTransitionStep] | None:
    """Return the largest duration whose accumulator stays non-negative.

    The next representable duration above the returned one overdraws the film,
    so the extinction root is located to one ULP.  ``None`` means the bracket
    could not be closed inside the declared window.
    """

    duration = analytic
    accepted: tuple[float, bc1.ValidatedQSCHexaneAppearanceTransitionStep] | None = None
    for _ in range(QSC_HEXANE_EPOCH_EXHAUSTION_ULP_WINDOW):
        if not duration > 0.0:
            return None
        try:
            accepted = (duration, _law_step(state, boundary, duration, None))
            break
        except bc1.QSCHexaneAppearanceTransitionError as exc:
            if exc.code is not bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR:
                raise
            duration = math.nextafter(duration, 0.0)
    if accepted is None:
        return None
    # Walk UP to the last duration that still keeps the film non-negative.
    for _ in range(QSC_HEXANE_EPOCH_EXHAUSTION_ULP_WINDOW):
        if accepted[1].film.after.attached_hexane_kg_per_kg_dry == 0.0:
            return accepted
        candidate = math.nextafter(accepted[0], math.inf)
        try:
            law = _law_step(state, boundary, candidate, None)
        except bc1.QSCHexaneAppearanceTransitionError as exc:
            if exc.code is not bc1.QSCHexaneAppearanceRefusalCode.NEGATIVE_FILM_ACCUMULATOR:
                raise
            # The next representable duration overdraws: the root is bracketed.
            return accepted
        accepted = (candidate, law)
    return None


# ---------------------------------------------------------------------------
# The typed birth gate (design D3)
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneDewClearanceSample(_FalseEpochClaims):
    """One chart sample: is the gas-only dry topology lawful at ``(T, y_b)``?"""

    site: str
    temperature_k: float
    pressure_pa: float
    hexane_mole_fraction: float
    hexane_activity: float
    water_activity: float
    admitted: bool
    chart_refusal: str | None

    def journal_record(self) -> dict:
        return {
            "site": self.site,
            "temperature_k": self.temperature_k,
            "hexane_mole_fraction": self.hexane_mole_fraction,
            "hexane_activity": self.hexane_activity,
            "water_activity": self.water_activity,
            "admitted": self.admitted,
            "chart_refusal": self.chart_refusal,
        }


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCHexaneBirthGateVerdict(_FalseEpochClaims):
    """The typed pre-birth admission verdict (design D3).

    ``admitted`` is the conjunction of the chart's OWN strict clearance across
    every sampled site and float equality of the attached-hexane film with
    ``0.0``.  Nothing else enters, and nothing here is a threshold.
    """

    admitted: bool
    film_kg_per_kg_dry: float
    film_exactly_zero: bool
    dew_clearance: bool
    newborn_temperature_k: float
    surface_temperature_k: float
    lumped_flake_surface_identity: bool
    boundary_hexane_mole_fraction: float
    pressure_pa: float
    samples: tuple[QSCHexaneDewClearanceSample, ...]
    epoch_time_s: float
    step_index: int

    def journal_record(self) -> dict:
        return {
            "schema_id": QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID,
            "gate": "qsc-hexane-birth-gate-d3",
            "admitted": self.admitted,
            "film_kg_per_kg_dry": self.film_kg_per_kg_dry,
            "film_exactly_zero": self.film_exactly_zero,
            "dew_clearance": self.dew_clearance,
            "newborn_temperature_k": self.newborn_temperature_k,
            "surface_temperature_k": self.surface_temperature_k,
            "lumped_flake_surface_identity": self.lumped_flake_surface_identity,
            "boundary_hexane_mole_fraction": self.boundary_hexane_mole_fraction,
            "pressure_pa": self.pressure_pa,
            "epoch_time_s": self.epoch_time_s,
            "step_index": self.step_index,
            "samples": [sample.journal_record() for sample in self.samples],
            "physically_qualifying": False,
        }


def _clearance_sample(
    site: str,
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
) -> QSCHexaneDewClearanceSample:
    """Ask the FROZEN gas-only chart itself whether this state is lawful."""

    refusal: str | None = None
    admitted = True
    try:
        cp.encode_gas_only_y(temperature_k, pressure_pa, y_hexane)
    except cp.CoupledPoreTopologyError as exc:
        admitted = False
        refusal = str(exc)
    except Exception as exc:  # noqa: BLE001 - anything else is not a verdict
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.CHART_EVALUATION_REFUSED,
            f"the gas-only chart could not be evaluated at {site}: {exc}",
        ) from exc
    try:
        water_activity, hexane_activity = cp.binary_gas_component_activities(
            temperature_k, pressure_pa, y_hexane
        )
    except Exception:  # noqa: BLE001 - diagnostics never override the verdict
        water_activity, hexane_activity = math.nan, math.nan
    return QSCHexaneDewClearanceSample(
        site=site,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        hexane_mole_fraction=y_hexane,
        hexane_activity=hexane_activity,
        water_activity=water_activity,
        admitted=admitted,
        chart_refusal=refusal,
    )


def evaluate_qsc_hexane_birth_gate(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    *,
    newborn_temperature_k: float | None = None,
) -> QSCHexaneBirthGateVerdict:
    """Evaluate D3's typed birth gate without taking any action.

    The would-be newborn piece is the outermost wet piece - the one the
    receding birth turns into the newborn dry shell - so its temperature
    defaults to that piece's own temperature, which under BC-1's declared
    lumped-flake convention is also the flake's surface temperature.  A caller
    holding an independent newborn-temperature candidate may declare it.
    """

    if type(state) is not QSCHexaneEpochState:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NOT_A_WET_WATER_STATE,
            "the birth gate consumes the exact sealed epoch state",
        )
    if type(boundary) is not QSCHexaneEpochBoundary:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.NONFINITE_INPUT,
            "the birth gate consumes the exact typed epoch boundary",
        )
    surface_temperature = state.interface_temperature_k
    if newborn_temperature_k is None:
        newborn = surface_temperature
        lumped_identity = True
    else:
        newborn = _require_finite("declared newborn temperature", newborn_temperature_k)
        lumped_identity = newborn == surface_temperature
    y_b = boundary.bulk_hexane_mole_fraction
    pressure = boundary.pressure_pa
    # The refusal locus the O9c campaign measured samples the gas at the
    # BOUNDARY composition across the film temperature span whose cold end is
    # the newborn dry piece.  Sample the whole span: the newborn piece, its
    # surface, and the bulk end.
    samples = (
        _clearance_sample("newborn_piece", newborn, pressure, y_b),
        _clearance_sample("newborn_surface", surface_temperature, pressure, y_b),
        _clearance_sample("film_bulk_end", boundary.bulk_temperature_k, pressure, y_b),
    )
    clearance = all(sample.admitted for sample in samples)
    film = state.film_kg_per_kg_dry
    film_zero = film == 0.0
    return QSCHexaneBirthGateVerdict(
        admitted=clearance and film_zero,
        film_kg_per_kg_dry=film,
        film_exactly_zero=film_zero,
        dew_clearance=clearance,
        newborn_temperature_k=newborn,
        surface_temperature_k=surface_temperature,
        lumped_flake_surface_identity=lumped_identity,
        boundary_hexane_mole_fraction=y_b,
        pressure_pa=pressure,
        samples=samples,
        epoch_time_s=state.epoch_time_s,
        step_index=state.step_index,
    )


def require_qsc_hexane_birth_admission(
    state: QSCHexaneEpochState,
    boundary: QSCHexaneEpochBoundary,
    *,
    newborn_temperature_k: float | None = None,
) -> QSCHexaneBirthGateVerdict:
    """Return the verdict when D3 admits a birth attempt; refuse otherwise."""

    verdict = evaluate_qsc_hexane_birth_gate(
        state, boundary, newborn_temperature_k=newborn_temperature_k
    )
    if not verdict.dew_clearance:
        blocking = [sample for sample in verdict.samples if not sample.admitted]
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.BIRTH_GATE_DEW_CLEARANCE_ABSENT,
            "the would-be newborn piece has not cleared the hexane dew line: "
            + "; ".join(
                f"{sample.site} at {sample.temperature_k} K -> {sample.chart_refusal}"
                for sample in blocking
            )
            + "; the epoch is lawfully DELAYED, birth is not refused",
            state=state,
            evidence=verdict.journal_record(),
        )
    if not verdict.film_exactly_zero:
        raise QSCHexaneEpochPrefixError(
            QSCHexaneEpochPrefixRefusalCode.BIRTH_GATE_FILM_PRESENT,
            f"the attached-hexane film still holds {verdict.film_kg_per_kg_dry} "
            "kg/kg dry: the radial activation endpoint requires attached hexane "
            "to be EXACTLY zero, so the epoch continues re-evaporating; refuse, "
            "never clamp",
            state=state,
            evidence=verdict.journal_record(),
        )
    return verdict


__all__ = [
    "QSC_HEXANE_EPOCH_EXHAUSTION_PROBE_HALVINGS",
    "QSC_HEXANE_EPOCH_EXHAUSTION_ULP_WINDOW",
    "QSC_HEXANE_EPOCH_PREFIX_ORCHESTRATION_ID",
    "QSC_HEXANE_EPOCH_PREFIX_SCHEMA_ID",
    "QSC_HEXANE_EPOCH_PREFIX_SCHEMA_REVISION",
    "QSC_HEXANE_EPOCH_PREFIX_SOURCE_IDENTITY",
    "QSC_HEXANE_EPOCH_WET_CORE_MAX_ITERATIONS",
    "QSC_HEXANE_EPOCH_WET_CORE_NEWTON_TOLERANCE",
    "QSCHexaneBirthGateVerdict",
    "QSCHexaneDewClearanceSample",
    "QSCHexaneEpochBoundary",
    "QSCHexaneEpochPrefixError",
    "QSCHexaneEpochPrefixRefusalCode",
    "QSCHexaneEpochState",
    "ValidatedQSCHexaneEpochStep",
    "advance_qsc_hexane_epoch_step",
    "advance_to_film_exhaustion",
    "evaluate_qsc_hexane_birth_gate",
    "initialize_qsc_hexane_epoch",
    "require_qsc_hexane_birth_admission",
]
