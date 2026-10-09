r"""Safeguarded same-cut-cell nonlinear integrator for Gate 1g.

The integrator wraps :mod:`cut_transport` without changing its equations.  It
solves the simultaneous backward-Euler wet/front/dry residual in transformed,
dimensionless coordinates:

* every temperature, retained-water loading, dry composition, and interface
  temperature lies strictly inside a caller-declared open interval;
* ``z`` lies strictly between the current cut cell's inner-face coordinate and
  the previous ``z``, so every accepted step is receding and remains on one
  fixed-rank chart; and
* every dry-face total Stefan flux is an unbounded algebraic coordinate divided
  by its caller-supplied positive physical scale.

Residual scaling is fixed before iteration.  The sparse production route uses
an exact unimodular row operation that exposes the conservative full-cut-cell
limit as the wet cut piece vanishes; the small dense route retains the original
rows as an independent numerical oracle.  Both routes are accepted only after
the original rows pass their unchanged scales and tolerance.  Energy scales use
datum-independent capacity rates and the dry and interface energy equations are
solved in the exactly equivalent basis ``R_E - H_REF*R_hexane``.  No penalty
residual, clipping, tolerance relaxation, contact fallback, or physics-changing
predictor is used.  A topology-invalid trial aborts the attempt rather than
being converted to a fake smooth residual.

A successful sparse reduced-film root that stops on ``xtol`` within the
frozen ``(1, 2]`` multiple of the unchanged residual tolerance may exercise
one deterministic scalar post-root safeguard.  It holds every coordinate but
the strongest structurally paired dry composition fixed, requires a
representable sign bracket, and then repeats the full residual, film, exact
smooth-cap, open-bound, rank, conditioning, and conservation certificates.
It is unavailable at cap contact, after any rejected topology trial, or under
a provisional work budget; failure commits nothing.

This stage deliberately excludes master-face arrival/crossing.  Such an event
requires a separate rank-changing event solve.  All failures carry the exact
immutable input :class:`CutIntegratorState` as ``rollback_state``.  The outer
Dirichlet trace and all coefficients remain numerical/non-qualifying.
"""

from __future__ import annotations

import enum
import math
import operator
from dataclasses import dataclass, field, replace
from typing import Callable, Sequence

import numpy as np
from scipy import optimize
from scipy.special import expit

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_sparsity
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.particle import wet_retained_cap as wrc
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wt

SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER = 2.0
# GT-PS-2-P1E-21 section 3.4.  The near-miss repair band is pinned at the
# absolute ceiling it was qualified against so that relaxing the convergence
# gate cannot silently widen a repair path.  At the A21 production gate the
# admission test requires ``maximum_scaled > 1e-10`` and ``<= 4e-11`` at once,
# so the band is empty and the safeguard is unreachable; for any caller at
# ``tolerance <= 2e-11`` the band is bit-for-bit unchanged.
SCALAR_POST_ROOT_NEAR_MISS_ABSOLUTE_CEILING = 4.0e-11
SCALAR_POST_ROOT_CENTRAL_RELATIVE_STEP = 1.0e-6
SCALAR_POST_ROOT_MAXIMUM_BRACKET_EXPANSIONS = 4
SCALAR_POST_ROOT_MAXIMUM_ITERATIONS = 64
SCALAR_POST_ROOT_REQUIRED_INTERFACE_TOLERANCE = 2.0e-11
P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING = 0.23354519496
NONLINEAR_RESIDUAL_SCALE_FLOOR = 1.0e-300
SPARSE_LSMR_NUMERICAL_AMENDMENT_ID = "GT-PS-2-P1E-08"
SPARSE_LSMR_ABSOLUTE_TOLERANCE = 1.0e-15
SPARSE_LSMR_RELATIVE_TOLERANCE = 1.0e-15
SPARSE_LSMR_MAXIMUM_ITERATIONS = 1000


class SameCellNonlinearSolverPolicy(enum.Enum):
    """Explicit numerical route for the unchanged same-cell residual.

    ``ESTABLISHED_ROUTE`` preserves the historical selection: sparse dogbox
    for a structurally coloured reduced-film Jacobian and dense TRF for the
    small-rank oracle.  ``SPARSE_TRF_LSMR`` is an opt-in alternative for the
    same bounded coordinates, row basis, structural Jacobian, LSMR accuracy,
    and original-residual certification.  It disables transformed-objective
    ``gtol`` termination because that gradient can converge before the
    separately certified original rows; ``ftol`` and ``xtol`` retain the
    caller's unchanged nonlinear-step tolerance.  It is numerical policy only:
    it changes no equation, chart, scale, acceptance tolerance, or seed.
    """

    ESTABLISHED_ROUTE = "established_route"
    SPARSE_TRF_LSMR = "sparse_trf_lsmr"


class SameCellFrontCoordinatePolicy(enum.Enum):
    """Numerical coordinate for the unchanged physical front primitive.

    ``WHOLE_CELL_OPEN_LOGIT`` is the established chart and remains the
    default.  ``DIRECT_PHYSICAL_FRONT_Z`` is an opt-in same-cell chart whose
    coordinate is the decoded binary64 physical ``front_z`` itself.  Its open
    bounds remain ``(z_inner_face, z_before)``.  It therefore changes only
    the nonlinear coordinate: the decoded candidate, residual equations,
    residual scales, physical gates, conservation ledgers, and tolerances are
    identical.  The optimizer guards are the first representable strict
    physical fronts at each end, so every coordinate it may evaluate decodes
    without clipping, subtraction, or endpoint projection.
    """

    WHOLE_CELL_OPEN_LOGIT = "whole_cell_open_logit"
    DIRECT_PHYSICAL_FRONT_Z = "direct_physical_front_z"


class NoUsableDirectPhysicalFrontInteriorError(ValueError):
    """The physical front interval has fewer than two ordered interior floats."""


def scalar_post_root_near_miss_ceiling(nonlinear_residual_tolerance: float) -> float:
    """Return the frozen GT-PS-2-P1E-21 near-miss admission ceiling.

    Runtime admission and every verifier that reconstructs the band must call
    this one expression, so evidence and solver cannot disagree about which
    candidates were eligible for the deterministic scalar post-root polish.
    """

    return min(
        SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER * nonlinear_residual_tolerance,
        SCALAR_POST_ROOT_NEAR_MISS_ABSOLUTE_CEILING,
    )


@dataclass(frozen=True)
class WetRetainedWaterApplicabilityAudit:
    """Exact bounded-P1EF audit for explicit wet retained-water states.

    The guard is an engineering applicability margin, not a capacity or phase
    boundary.  Dry local-equilibrium retained water is deliberately outside
    this audit.
    """

    policy_enabled: bool
    guard_loading: float
    maximum_loading: float
    minimum_eta_cap: float
    active_piece_count: int
    passed: bool
    physically_qualifying: bool = field(default=False, init=False)


class WetRetainedWaterApplicabilityError(RuntimeError):
    """A fully certified root lies outside the selected smooth P1EF scope."""

    def __init__(
        self,
        message: str,
        audit: WetRetainedWaterApplicabilityAudit,
    ) -> None:
        super().__init__(message)
        self.audit = audit


class ConditionProxyUncertifiedError(RuntimeError):
    """Typed cause: the certified route's condition gate refused the solve.

    Added by GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29 so the
    condition-cap cause is type-discriminable from a solver not-success
    refusal (whose wrapper proxy may also exceed the cap).  Additive: the
    message, trigger condition, and gate value are unchanged, and every
    existing ``RuntimeError`` handler behaves exactly as before.  The type
    is meaningful in-process only via the wrapper's ``__cause__`` chain;
    serialized evidence retains the proxy as an exact binary64 value and
    never reconstructs this object.
    """

    def __init__(self, *, condition_proxy: float, maximum_condition_proxy: float) -> None:
        super().__init__(
            f"nonlinear Jacobian condition proxy is uncertified: {condition_proxy:.3e}"
        )
        self.condition_proxy = condition_proxy
        self.maximum_condition_proxy = maximum_condition_proxy


class CutIntegratorStepError(RuntimeError):
    """Rejected nonlinear step with exact state identity for rollback."""

    def __init__(
        self,
        message: str,
        rollback_state: "CutIntegratorState",
        *,
        nonlinear_evaluations: int,
        rejected_trial_evaluations: int,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        event_restart_required: bool = False,
        last_candidate: cut.CutTransportUnknowns | None = None,
        scaled_residuals: tuple[float, ...] = (),
        residual_blocks: cut.CutResidualBlocks | None = None,
        minimum_fractional_distance_to_bound: float = 0.0,
        nonlinear_optimizer_function_evaluations: int = 0,
        provisional_work_budget: int | None = None,
        provisional_work_budget_exhausted: bool = False,
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: WetRetainedWaterApplicabilityAudit | None = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.rejected_trial_evaluations = rejected_trial_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.event_restart_required = event_restart_required
        self.last_candidate = last_candidate
        self.scaled_residuals = scaled_residuals
        self.residual_blocks = residual_blocks
        self.minimum_fractional_distance_to_bound = minimum_fractional_distance_to_bound
        self.nonlinear_optimizer_function_evaluations = nonlinear_optimizer_function_evaluations
        self.provisional_work_budget = provisional_work_budget
        self.provisional_work_budget_exhausted = provisional_work_budget_exhausted
        self.retained_water_applicability_guard_failed = retained_water_applicability_guard_failed
        self.retained_water_applicability_audit = retained_water_applicability_audit


@dataclass(frozen=True)
class CutSolverControls:
    """Caller-declared open charts and immutable numerical acceptance contract."""

    wet_temperature_bounds_k: tuple[float, float]
    wet_water_bounds: tuple[float, float]
    interface_temperature_bounds_k: tuple[float, float]
    # GT-PS-2-P1E-21.  The scaled residual is dimensionless, so this gate is a
    # significant-digit requirement on the balance closure.  It was relaxed
    # from 2.0e-11 to the 1.0e-10 raw conservation authority it exists to
    # serve: the old value sat below the ~2.27e-11 binary64 attainable floor
    # with no engineering margin.  Raw ledgers, condition ceiling, leaf floor,
    # work budgets, and every equation are unchanged.
    nonlinear_residual_tolerance: float = 1.0e-10
    ledger_tolerance: float = 1.0e-10
    nonlinear_step_tolerance: float = 1.0e-10
    maximum_function_evaluations: int = 800
    maximum_condition_proxy: float = 1.0e14
    maximum_logit_magnitude: float = 30.0
    maximum_dense_unknowns: int = 20
    reduced_film_structural_sparse_jacobian: bool = True
    enforce_p1ef_smooth_wet_retained_water_guard: bool = False
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        for name, bounds in (
            ("wet-temperature", self.wet_temperature_bounds_k),
            ("wet-water", self.wet_water_bounds),
            ("interface-temperature", self.interface_temperature_bounds_k),
        ):
            lower, upper = bounds
            if not all(math.isfinite(value) for value in bounds) or not lower < upper:
                raise ValueError(f"{name} bounds must be finite and ordered")
        tolerances = (
            self.nonlinear_residual_tolerance,
            self.ledger_tolerance,
            self.nonlinear_step_tolerance,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in tolerances):
            raise ValueError("solver tolerances must be positive and finite")
        if self.nonlinear_residual_tolerance > 1.0e-10:
            raise ValueError("regional/RH tolerance cannot exceed 1e-10")
        if self.ledger_tolerance > 1.0e-10:
            raise ValueError("accepted/cumulative ledger tolerance cannot exceed 1e-10")
        try:
            maximum = operator.index(self.maximum_function_evaluations)
        except TypeError as exc:
            raise ValueError("maximum function evaluations must be an integer") from exc
        if isinstance(self.maximum_function_evaluations, bool) or maximum < 1:
            raise ValueError("maximum function evaluations must be positive")
        if not math.isfinite(self.maximum_condition_proxy) or self.maximum_condition_proxy <= 1.0:
            raise ValueError("condition-proxy ceiling must be finite and exceed one")
        if (
            not math.isfinite(self.maximum_logit_magnitude)
            or not 10.0 <= self.maximum_logit_magnitude <= 35.0
        ):
            raise ValueError("maximum logit magnitude must lie in [10, 35]")
        try:
            dense_unknowns = operator.index(self.maximum_dense_unknowns)
        except TypeError as exc:
            raise ValueError("maximum dense unknowns must be an integer") from exc
        if isinstance(self.maximum_dense_unknowns, bool) or dense_unknowns < 1:
            raise ValueError("maximum dense unknown count must be positive")
        if not isinstance(self.reduced_film_structural_sparse_jacobian, bool):
            raise ValueError("reduced-film sparse-Jacobian policy must be boolean")
        if not isinstance(
            self.enforce_p1ef_smooth_wet_retained_water_guard,
            bool,
        ):
            raise ValueError("P1EF smooth retained-water guard policy must be boolean")


def certify_explicit_wet_retained_water_applicability(
    controls: CutSolverControls,
    retained_water_loadings: Sequence[float],
    capacity_duals_over_rt: Sequence[float],
    luikov: sp.LuikovParams,
) -> WetRetainedWaterApplicabilityAudit:
    """Audit, and when selected enforce, the exact smooth P1EF wet guard.

    ``W <= W_guard`` and ``lambda == 0`` are exact comparisons.  There is no
    tolerance, clipping, projection, phase selection, or dry-pore check here.
    """

    if not isinstance(controls, CutSolverControls):
        raise TypeError("wet retained-water applicability needs CutSolverControls")
    loadings = tuple(retained_water_loadings)
    duals = tuple(capacity_duals_over_rt)
    if not loadings or len(loadings) != len(duals):
        raise ValueError("wet retained-water applicability needs equally ranked non-empty states")
    if not all(math.isfinite(value) for value in (*loadings, *duals)):
        raise ValueError("wet retained-water applicability states must be finite")
    lower = controls.wet_water_bounds[0]
    denominator = luikov.W_cap - lower
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("wet retained-water applicability chart has no cap span")
    maximum_loading = max(loadings)
    minimum_eta_cap = min((luikov.W_cap - loading) / denominator for loading in loadings)
    active_piece_count = sum(dual != 0.0 for dual in duals)
    enabled = controls.enforce_p1ef_smooth_wet_retained_water_guard
    passed = bool(
        not enabled
        or (
            maximum_loading <= P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING
            and active_piece_count == 0
        )
    )
    audit = WetRetainedWaterApplicabilityAudit(
        policy_enabled=enabled,
        guard_loading=P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING,
        maximum_loading=maximum_loading,
        minimum_eta_cap=minimum_eta_cap,
        active_piece_count=active_piece_count,
        passed=passed,
    )
    if not passed:
        reason = (
            "an explicit wet retained-water loading exceeded the frozen "
            f"P1EF guard ({maximum_loading!r} > "
            f"{P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING!r})"
            if maximum_loading > P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING
            else (
                "the selected smooth P1EF retained-water route encountered "
                f"{active_piece_count} capacity-active wet piece(s)"
            )
        )
        raise WetRetainedWaterApplicabilityError(reason, audit)
    return audit


def certify_enabled_wet_retained_water_state_applicability(
    controls: object,
    state: object,
    *,
    capacity_duals_known_zero: bool = False,
    luikov: sp.LuikovParams | None = None,
) -> WetRetainedWaterApplicabilityAudit | None:
    """Enforce the selected guard without burdening disabled protocol states.

    Several orchestration tests deliberately use narrow structural state
    doubles.  More importantly, the default-disabled policy is intended to be
    a genuine no-op.  Attribute access is therefore deferred until the typed
    controls explicitly enable the guard.  Enabled production states must
    expose the complete explicit wet-water payload and are then audited by the
    exact public helper above.
    """

    enabled = bool(
        getattr(
            controls,
            "enforce_p1ef_smooth_wet_retained_water_guard",
            False,
        )
    )
    if not enabled:
        return None
    if not isinstance(controls, CutSolverControls):
        raise TypeError("enabled wet retained-water guard needs CutSolverControls")
    try:
        loadings = tuple(getattr(state, "wet_retained_water_loadings"))
        if capacity_duals_known_zero:
            duals = (0.0,) * len(loadings)
        else:
            duals = tuple(
                getattr(
                    state,
                    "effective_wet_retained_water_capacity_duals_over_rt",
                )
            )
        selected_luikov = (
            getattr(getattr(getattr(state, "config"), "wet"), "luikov")
            if luikov is None
            else luikov
        )
    except AttributeError as exc:
        raise TypeError(
            "enabled wet retained-water guard needs a complete explicit wet state"
        ) from exc
    return certify_explicit_wet_retained_water_applicability(
        controls,
        loadings,
        duals,
        selected_luikov,
    )


@dataclass(frozen=True)
class CutStepSeed:
    """One explicit nonlinear seed and physical scaling for algebraic ``N_t``."""

    candidate: cut.CutTransportUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, cut.CutTransportUnknowns):
            raise TypeError("step seed requires CutTransportUnknowns")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0 for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every algebraic Stefan face needs a positive physical scale")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("step seed needs a non-empty provenance label")


@dataclass(frozen=True)
class MaterialInventorySnapshot:
    """Component and common-energy inventories on stationary material cells."""

    water_cell_mol: tuple[float, ...]
    hexane_cell_mol: tuple[float, ...]
    energy_cell_j: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def total_water_mol(self) -> float:
        return math.fsum(self.water_cell_mol)

    @property
    def total_hexane_mol(self) -> float:
        return math.fsum(self.hexane_cell_mol)

    @property
    def total_energy_j(self) -> float:
        return math.fsum(self.energy_cell_j)


@dataclass(frozen=True)
class WaterPhaseInventorySnapshot:
    """Whole-particle water partition on the frozen Gate-1g topology.

    ``retained_water_mol`` includes the wet core and dry-matrix retained
    water.  ``pore_vapor_water_mol`` is the dry gas-connected pore inventory.
    The current no-external-liquid topology has exactly zero free-liquid
    inventory; reporting that zero explicitly prevents it being confused with
    an untracked reservoir.
    """

    retained_water_mol: float
    pore_vapor_water_mol: float
    free_liquid_water_mol: float = 0.0
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def total_water_mol(self) -> float:
        return math.fsum(
            (
                self.retained_water_mol,
                self.pore_vapor_water_mol,
                self.free_liquid_water_mol,
            )
        )


@dataclass(frozen=True)
class CutIntegratorState:
    """Accepted differential state plus immutable cumulative references.

    Signed boundary and per-cell material changes use persisted Neumaier
    ``high + compensation`` pairs.  The ``corrected_cumulative_boundary_*``
    properties expose the signed boundary totals used by the ledger; monotone
    absolute-transfer fields are independent scale histories.
    """

    transport: cut.CutTransportState
    controls: CutSolverControls
    reference_inventory: MaterialInventorySnapshot
    reference_capacity_energy_scale_j: float
    last_interface_temperature_k: float
    last_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    cumulative_boundary_water_out_mol: float = 0.0
    cumulative_boundary_hexane_out_mol: float = 0.0
    cumulative_boundary_energy_out_j: float = 0.0
    boundary_water_compensation_mol: float = 0.0
    boundary_hexane_compensation_mol: float = 0.0
    boundary_energy_compensation_j: float = 0.0
    cumulative_absolute_water_transfer_mol: float = 0.0
    cumulative_absolute_hexane_transfer_mol: float = 0.0
    cumulative_absolute_energy_transfer_j: float = 0.0
    cumulative_material_water_change_cell_mol: tuple[float, ...] = ()
    cumulative_material_hexane_change_cell_mol: tuple[float, ...] = ()
    cumulative_material_energy_change_cell_j: tuple[float, ...] = ()
    water_change_compensation_cell_mol: tuple[float, ...] = ()
    hexane_change_compensation_cell_mol: tuple[float, ...] = ()
    energy_change_compensation_cell_j: tuple[float, ...] = ()
    accepted_steps: int = 0
    cumulative_nonlinear_evaluations: int = 0
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.transport, cut.CutTransportState):
            raise TypeError("integrator state requires a cut-transport state")
        if not isinstance(self.controls, CutSolverControls):
            raise TypeError("integrator state requires solver controls")
        _validate_controls_against_state(self.transport, self.controls)
        if not math.isfinite(self.reference_capacity_energy_scale_j) or (
            self.reference_capacity_energy_scale_j <= 0.0
        ):
            raise ValueError("reference capacity-energy scale must be positive")
        if not _inside(
            self.last_interface_temperature_k,
            self.controls.interface_temperature_bounds_k,
        ):
            raise ValueError("last interface temperature left its open chart")
        if len(self.last_total_stefan_fluxes_mol_m2_s) != (self.transport.layout.dry_face_count):
            raise ValueError("last Stefan fluxes must align with dry faces")
        if not all(math.isfinite(value) for value in self.last_total_stefan_fluxes_mol_m2_s):
            raise ValueError("last Stefan fluxes must be finite")
        values = (
            self.cumulative_boundary_water_out_mol,
            self.cumulative_boundary_hexane_out_mol,
            self.cumulative_boundary_energy_out_j,
            self.boundary_water_compensation_mol,
            self.boundary_hexane_compensation_mol,
            self.boundary_energy_compensation_j,
            self.cumulative_absolute_water_transfer_mol,
            self.cumulative_absolute_hexane_transfer_mol,
            self.cumulative_absolute_energy_transfer_j,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("cumulative ledger values must remain finite")
        if any(value < 0.0 for value in values[6:]):
            raise ValueError("cumulative absolute transfers cannot be negative")
        for name, vector in (
            ("cumulative material-water changes", self.cumulative_material_water_change_cell_mol),
            ("cumulative material-hexane changes", self.cumulative_material_hexane_change_cell_mol),
            ("cumulative material-energy changes", self.cumulative_material_energy_change_cell_j),
            ("water summation compensation", self.water_change_compensation_cell_mol),
            ("hexane summation compensation", self.hexane_change_compensation_cell_mol),
            ("energy summation compensation", self.energy_change_compensation_cell_j),
        ):
            if len(vector) != self.transport.geometry.master_grid.n:
                raise ValueError(f"{name} must align with stationary material cells")
            if not all(math.isfinite(value) for value in vector):
                raise ValueError(f"{name} must remain finite")
        if self.accepted_steps < 0 or self.cumulative_nonlinear_evaluations < 0:
            raise ValueError("cumulative counters cannot be negative")

    @property
    def corrected_cumulative_boundary_water_out_mol(self) -> float:
        return math.fsum(
            (
                self.cumulative_boundary_water_out_mol,
                self.boundary_water_compensation_mol,
            )
        )

    @property
    def corrected_cumulative_boundary_hexane_out_mol(self) -> float:
        return math.fsum(
            (
                self.cumulative_boundary_hexane_out_mol,
                self.boundary_hexane_compensation_mol,
            )
        )

    @property
    def corrected_cumulative_boundary_energy_out_j(self) -> float:
        return math.fsum(
            (
                self.cumulative_boundary_energy_out_j,
                self.boundary_energy_compensation_j,
            )
        )


@dataclass(frozen=True)
class DryEnergyCapacityStencilAudit:
    """Certificate for one Amendment-19 dry-energy normalization derivative."""

    temperature_k: float
    pressure_pa: float
    y_hexane: float
    actual_oil_fraction_label: float
    gas_accessible_fraction: float
    primitive_domain_kind: str
    stencil_kind: str
    delta_temperature_k: float
    sample_temperatures_k: tuple[float, ...]
    sample_full_coupled_pore_certified: tuple[bool, ...]
    center_gas_only_admissible: bool | None
    central_lower_gas_only_admissible: bool | None
    central_upper_gas_only_admissible: bool | None
    forward_second_gas_only_admissible: bool | None
    backward_second_gas_only_admissible: bool | None
    capacity_j_m3_k: float
    numerical_amendment_id: str = field(default="GT-PS-2-P1E-19", init=False)
    normalization_only: bool = field(default=True, init=False)
    accepted_residual_equations_changed: bool = field(default=False, init=False)
    composition_or_activity_clipping_used: bool = field(default=False, init=False)
    property_extrapolation_used: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        scalar_values = (
            self.temperature_k,
            self.pressure_pa,
            self.y_hexane,
            self.actual_oil_fraction_label,
            self.gas_accessible_fraction,
            self.delta_temperature_k,
            self.capacity_j_m3_k,
        )
        if not all(math.isfinite(value) for value in scalar_values):
            raise ValueError("dry-capacity stencil audit values must be finite")
        if self.pressure_pa <= 0.0 or not 0.0 < self.y_hexane < 1.0:
            raise ValueError("dry-capacity stencil primitive is invalid")
        if self.gas_accessible_fraction != 1.0:
            raise ValueError("dry-capacity normalization requires the full dry volume")
        if self.delta_temperature_k <= 0.0 or self.capacity_j_m3_k <= 0.0:
            raise ValueError("dry-capacity stencil increment and result must be positive")
        if not self.sample_temperatures_k or not all(
            math.isfinite(value) for value in self.sample_temperatures_k
        ):
            raise ValueError("dry-capacity stencil sample temperatures are invalid")
        if len(self.sample_full_coupled_pore_certified) != len(
            self.sample_temperatures_k
        ) or not all(self.sample_full_coupled_pore_certified):
            raise ValueError("every dry-capacity sample needs full pore certification")
        allowed = {
            "gas_only_second_order_centered",
            "gas_only_second_order_forward",
            "gas_only_second_order_backward",
            "joint_primitive_legacy_centered_unchanged",
        }
        if self.stencil_kind not in allowed:
            raise ValueError("dry-capacity stencil kind is unknown")
        if self.primitive_domain_kind == "GasOnlyPrimitiveDomain":
            flags = (
                self.center_gas_only_admissible,
                self.central_lower_gas_only_admissible,
                self.central_upper_gas_only_admissible,
            )
            if not all(isinstance(value, bool) for value in flags):
                raise ValueError("gas-only capacity audit lacks topology decisions")
            if self.center_gas_only_admissible is not True:
                raise ValueError("dry-capacity center must remain in the gas-only chart")
            if self.stencil_kind == "gas_only_second_order_centered" and not (
                self.central_lower_gas_only_admissible is True
                and self.central_upper_gas_only_admissible is True
                and len(self.sample_temperatures_k) == 2
            ):
                raise ValueError("centered dry-capacity audit lacks both admissible samples")
            if self.stencil_kind == "gas_only_second_order_forward" and not (
                self.central_upper_gas_only_admissible is True
                and self.forward_second_gas_only_admissible is True
                and len(self.sample_temperatures_k) == 3
            ):
                raise ValueError("forward dry-capacity audit lacks both admissible samples")
            if self.stencil_kind == "gas_only_second_order_backward" and not (
                self.central_lower_gas_only_admissible is True
                and self.backward_second_gas_only_admissible is True
                and len(self.sample_temperatures_k) == 3
            ):
                raise ValueError("backward dry-capacity audit lacks both admissible samples")
        elif self.primitive_domain_kind == "JointPrimitiveBand":
            if self.stencil_kind != "joint_primitive_legacy_centered_unchanged":
                raise ValueError("legacy primitive band must retain its centered stencil")
            if any(
                value is not None
                for value in (
                    self.center_gas_only_admissible,
                    self.central_lower_gas_only_admissible,
                    self.central_upper_gas_only_admissible,
                    self.forward_second_gas_only_admissible,
                    self.backward_second_gas_only_admissible,
                )
            ):
                raise ValueError("legacy primitive band cannot claim gas-only certification")
            if len(self.sample_temperatures_k) != 2:
                raise ValueError("legacy centered dry-capacity audit needs two samples")
        else:
            raise ValueError("dry-capacity primitive-domain kind is unknown")


@dataclass(frozen=True)
class CutResidualScales:
    """Fixed row-local physical scales in the residual block ordering."""

    wet_water_mol_s: tuple[float, ...]
    wet_energy_w: tuple[float, ...]
    dry_water_mol_s: tuple[float, ...]
    dry_hexane_mol_s: tuple[float, ...]
    dry_energy_w: tuple[float, ...]
    rh_water_mol_s: float
    rh_hexane_mol_s: float
    rh_energy_w: float
    dry_energy_capacity_stencil_audits: tuple[
        DryEnergyCapacityStencilAudit,
        ...,
    ]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if len(self.dry_energy_capacity_stencil_audits) < len(self.dry_energy_w):
            raise ValueError("every dry-energy residual scale needs an Amendment-19 stencil audit")
        if not all(
            isinstance(audit, DryEnergyCapacityStencilAudit)
            for audit in self.dry_energy_capacity_stencil_audits
        ):
            raise TypeError("dry-energy capacity audit payload has the wrong type")

    @property
    def vector(self) -> tuple[float, ...]:
        values = (
            *self.wet_water_mol_s,
            *self.wet_energy_w,
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_energy_w,
            self.rh_water_mol_s,
            self.rh_hexane_mol_s,
            self.rh_energy_w,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise RuntimeError("every nonlinear row scale must be positive and finite")
        return values


@dataclass(frozen=True)
class ScalarPostRootPolishAudit:
    """Deterministic certificate for the narrowly gated scalar safeguard.

    The safeguard is an outer nonlinear-algorithm operation.  It changes no
    residual equation, physical bound, interface solve, acceptance tolerance,
    or conservation limit.  A payload exists only on a completely accepted
    step; every failed polish is discarded with the caller's exact immutable
    input state.
    """

    near_miss_multiplier_ceiling: float
    maximum_scaled_residual_before: float
    dominant_residual_block: str
    dominant_residual_local_index: int
    dominant_residual_global_index: int
    selected_coordinate_block: str
    selected_coordinate_local_index: int
    selected_coordinate_global_index: int
    selected_coordinate_derivative: float
    selection_stencil_evaluations: int
    predicted_scalar_correction: float
    predicted_scalar_correction_ulps: float
    bracket_probe_count: int
    bracket_lower_coordinate: float
    bracket_upper_coordinate: float
    bracket_lower_row_residual: float
    bracket_upper_row_residual: float
    scalar_iterations: int
    scalar_function_calls: int
    scalar_root_coordinate: float
    scalar_root_delta: float
    maximum_scaled_residual_after: float
    final_condition_proxy: float
    final_jacobian_rank: int
    final_jacobian_dimension: int
    final_condition_stencil_evaluations: int
    used_generic_interface_tolerance_unchanged: bool
    no_clipping_projection_or_tolerance_relaxation: bool = field(
        default=True,
        init=False,
    )
    same_topology_and_rank: bool = field(default=True, init=False)
    partial_candidate_committed: bool = field(default=False, init=False)
    exact_rollback_on_failure: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        finite_values = (
            self.near_miss_multiplier_ceiling,
            self.maximum_scaled_residual_before,
            self.selected_coordinate_derivative,
            self.predicted_scalar_correction,
            self.predicted_scalar_correction_ulps,
            self.bracket_lower_coordinate,
            self.bracket_upper_coordinate,
            self.bracket_lower_row_residual,
            self.bracket_upper_row_residual,
            self.scalar_root_coordinate,
            self.scalar_root_delta,
            self.maximum_scaled_residual_after,
            self.final_condition_proxy,
        )
        if not all(math.isfinite(value) for value in finite_values):
            raise ValueError("scalar post-root audit values must be finite")
        if self.near_miss_multiplier_ceiling != SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER:
            raise ValueError("scalar post-root near-miss ceiling is not the frozen value")
        if not self.bracket_lower_coordinate < self.bracket_upper_coordinate:
            raise ValueError("scalar post-root bracket must be strictly ordered")
        if self.bracket_lower_row_residual * self.bracket_upper_row_residual > 0.0:
            raise ValueError("scalar post-root bracket does not contain a row sign change")
        integer_values = (
            self.dominant_residual_local_index,
            self.dominant_residual_global_index,
            self.selected_coordinate_local_index,
            self.selected_coordinate_global_index,
            self.selection_stencil_evaluations,
            self.bracket_probe_count,
            self.scalar_iterations,
            self.scalar_function_calls,
            self.final_jacobian_rank,
            self.final_jacobian_dimension,
            self.final_condition_stencil_evaluations,
        )
        if any(value < 0 for value in integer_values):
            raise ValueError("scalar post-root audit counters and indices must be nonnegative")
        if (
            self.selection_stencil_evaluations < 2
            or self.bracket_probe_count < 1
            or self.scalar_function_calls < 1
        ):
            raise ValueError("scalar post-root audit must disclose nonzero work")
        if self.final_jacobian_rank != self.final_jacobian_dimension:
            raise ValueError("scalar post-root final Jacobian is not full rank")
        if not self.used_generic_interface_tolerance_unchanged:
            raise ValueError("scalar safeguard cannot change the interface tolerance")


@dataclass(frozen=True)
class CutIntegratorLedger:
    """Nonlinear, step, cumulative, conditioning, and bound diagnostics."""

    dt_s: float
    pressure_before_pa: float
    pressure_after_pa: float
    finite_pressure_transition: bool
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    nonlinear_message: str
    residual_scales: CutResidualScales
    optimization_residual_scales: CutResidualScales
    scaled_residuals: tuple[float, ...]
    maximum_scaled_residual: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    minimum_front_z_distance_to_chart_boundary: float
    used_sparse_jacobian: bool
    jacobian_route_reason: str
    front_coordinate_policy: str
    nonlinear_algorithm: str
    nonlinear_residual_basis: str
    scalability_class: str
    jacobian_structural_color_count: int
    jacobian_structural_nnz: int
    water_step_residual_mol: float
    hexane_step_residual_mol: float
    energy_step_residual_j: float
    normalized_water_step_residual: float
    normalized_hexane_step_residual: float
    normalized_energy_step_residual: float
    water_cumulative_residual_mol: float
    hexane_cumulative_residual_mol: float
    energy_cumulative_residual_j: float
    normalized_water_cumulative_residual: float
    normalized_hexane_cumulative_residual: float
    normalized_energy_cumulative_residual: float
    boundary_water_out_mol: float
    boundary_hexane_out_mol: float
    boundary_energy_out_j: float
    surface_film_audit: ct.SurfaceFilmAudit | None
    vanishing_wet_cut_treatment_audit: cut.VanishingWetCutTreatmentAudit
    wet_retained_cap_active_piece_indices: tuple[int, ...]
    wet_retained_cap_active_piece_count: int
    wet_retained_cap_maximum_loading: float
    wet_retained_cap_maximum_capacity_dual_over_rt: float
    wet_retained_cap_maximum_complementarity_product: float
    wet_retained_cap_exact_graph_without_tolerance: bool
    wet_retained_cap_piece_transitions: tuple[
        wrc.WetRetainedCapTransition,
        ...,
    ]
    wet_retained_cap_entry_piece_indices: tuple[int, ...]
    wet_retained_cap_entry_piece_count: int
    wet_retained_cap_continuation_piece_indices: tuple[int, ...]
    wet_retained_cap_continuation_piece_count: int
    wet_retained_cap_exit_piece_indices: tuple[int, ...]
    wet_retained_cap_exit_piece_count: int
    accepted: bool
    scalar_post_root_polish_audit: ScalarPostRootPolishAudit | None = None
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def maximum_step_ledger_residual(self) -> float:
        return max(
            self.normalized_water_step_residual,
            self.normalized_hexane_step_residual,
            self.normalized_energy_step_residual,
        )

    @property
    def maximum_cumulative_ledger_residual(self) -> float:
        return max(
            self.normalized_water_cumulative_residual,
            self.normalized_hexane_cumulative_residual,
            self.normalized_energy_cumulative_residual,
        )


@dataclass(frozen=True)
class CutIntegratorStep:
    """One certified accepted same-cut-cell primary-drainage step."""

    before: CutIntegratorState
    after: CutIntegratorState
    assembly: cut.CutTransportAssembly
    boundary: ct.PoreBoundary
    seed: CutStepSeed
    ledger: CutIntegratorLedger
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def surface_film_audit(self) -> ct.SurfaceFilmAudit | None:
        return self.assembly.surface_film_audit


@dataclass(frozen=True)
class _PieceData:
    wet_cells: tuple[ww.WetWaterCellState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_water_mol: tuple[float, ...]
    wet_hexane_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float, ...]
    dry_hexane_mol: tuple[float, ...]
    dry_energy_j: tuple[float, ...]
    wet_capacity_energy_j: tuple[float, ...]
    dry_capacity_energy_j: tuple[float, ...]
    dry_capacity_stencil_audits: tuple[DryEnergyCapacityStencilAudit, ...]


def initialize_integrator_state(
    transport: cut.CutTransportState,
    controls: CutSolverControls,
    *,
    interface_temperature_k: float,
    total_stefan_fluxes_mol_m2_s: Sequence[float] | None = None,
) -> CutIntegratorState:
    """Create cumulative references from one validated strict partial state."""

    _validate_controls_against_state(transport, controls)
    cut_index = transport.layout.cut_cell_index
    cut.evaluate_interface_state(
        interface_temperature_k,
        transport.config,
        transport.historical_hexane_loadings[cut_index],
        transport.oil_fraction_labels[cut_index],
    )
    snapshot = inventory_snapshot(transport)
    capacity_scale = capacity_energy_scale(transport)
    if total_stefan_fluxes_mol_m2_s is None:
        stefan = (0.0,) * transport.layout.dry_face_count
    else:
        stefan = tuple(total_stefan_fluxes_mol_m2_s)
    return CutIntegratorState(
        transport=transport,
        controls=controls,
        reference_inventory=snapshot,
        reference_capacity_energy_scale_j=capacity_scale,
        last_interface_temperature_k=interface_temperature_k,
        last_total_stefan_fluxes_mol_m2_s=stefan,
        cumulative_material_water_change_cell_mol=(0.0,) * transport.geometry.master_grid.n,
        cumulative_material_hexane_change_cell_mol=(0.0,) * transport.geometry.master_grid.n,
        cumulative_material_energy_change_cell_j=(0.0,) * transport.geometry.master_grid.n,
        water_change_compensation_cell_mol=(0.0,) * transport.geometry.master_grid.n,
        hexane_change_compensation_cell_mol=(0.0,) * transport.geometry.master_grid.n,
        energy_change_compensation_cell_j=(0.0,) * transport.geometry.master_grid.n,
    )


def inventory_snapshot(state: cut.CutTransportState) -> MaterialInventorySnapshot:
    """Evaluate exact wet+dry inventories on each stationary material cell."""

    # Inventory is a state-function evaluation.  It must not depend on the
    # auxiliary temperature stencil used only to scale energy residuals: at a
    # valid gas-only state that fixed-composition stencil may leave the local
    # phase chart even though the state and all three inventories are valid.
    pieces = _piece_data(state, evaluate_capacity=False)
    n = state.geometry.master_grid.n
    water = [0.0] * n
    hexane = [0.0] * n
    energy = [0.0] * n
    for position, cell in enumerate(state.layout.wet_cell_indices):
        water[cell] += pieces.wet_water_mol[position]
        hexane[cell] += pieces.wet_hexane_mol[position]
        energy[cell] += pieces.wet_energy_j[position]
    for position, cell in enumerate(state.layout.dry_cell_indices):
        water[cell] += pieces.dry_water_mol[position]
        hexane[cell] += pieces.dry_hexane_mol[position]
        energy[cell] += pieces.dry_energy_j[position]
    return MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def water_phase_inventory_snapshot(
    state: cut.CutTransportState,
) -> WaterPhaseInventorySnapshot:
    """Evaluate retained, dry-pore-vapor, and free-liquid water separately."""

    pieces = _piece_data(state, evaluate_capacity=False)
    dry_volumes = tuple(
        state.geometry.cells[cell].dry_volume_m3 for cell in state.layout.dry_cell_indices
    )
    dry_retained = math.fsum(
        volume * cell.retained_water_concentration_kg_m3 / wt.M
        for volume, cell in zip(dry_volumes, pieces.dry_cells)
    )
    dry_pore_vapor = math.fsum(
        volume * cell.pore_water_concentration_kg_m3 / wt.M
        for volume, cell in zip(dry_volumes, pieces.dry_cells)
    )
    snapshot = WaterPhaseInventorySnapshot(
        retained_water_mol=math.fsum((*pieces.wet_water_mol, dry_retained)),
        pore_vapor_water_mol=dry_pore_vapor,
    )
    total = inventory_snapshot(state).total_water_mol
    scale = max(abs(total), abs(snapshot.total_water_mol), 1.0e-300)
    if abs(math.fsum((snapshot.total_water_mol, -total))) / scale > 2.0e-14:
        raise RuntimeError("water phase partition does not close total inventory")
    return snapshot


def capacity_energy_scale(state: cut.CutTransportState) -> float:
    """Return the datum-independent whole-particle ``sum(V C_T T)`` scale."""

    pieces = _piece_data(state)
    scale = math.fsum((*pieces.wet_capacity_energy_j, *pieces.dry_capacity_energy_j))
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("capacity-energy scale must be positive and finite")
    return scale


def _jacobian_route_reason(
    before: CutIntegratorState,
    boundary: ct.PoreBoundary,
) -> str:
    """Select a numerical Jacobian route without changing the residual.

    An exact-A reduced-film residual evaluates a safeguarded scalar surface
    closure on every coupled residual call.  Repeating that closure once per
    dense finite-difference column is avoidable work even at small rank, so
    the explicitly enabled reduced-film policy uses the already verified
    bordered-banded structural coloring at every mesh.  Legacy Dirichlet
    small-rank roots retain the dense numerical oracle unless their declared
    unknown-count threshold is exceeded.
    """

    if (
        isinstance(boundary, ct.ReducedFilmPoreBoundary)
        and before.controls.reduced_film_structural_sparse_jacobian
    ):
        return "exact_ackermann_reduced_film_structural_coloring"
    if before.transport.layout.unknown_count > before.controls.maximum_dense_unknowns:
        return "declared_unknown_count_threshold"
    return "small_rank_dense_oracle"


def advance_same_cell_backward_euler(
    before: CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    seed: CutStepSeed,
    *,
    target_config: cut.CutTransportConfig | None = None,
    provisional_maximum_function_evaluations: int | None = None,
    nonlinear_solver_policy: SameCellNonlinearSolverPolicy = (
        SameCellNonlinearSolverPolicy.ESTABLISHED_ROUTE
    ),
    front_coordinate_policy: SameCellFrontCoordinatePolicy = (
        SameCellFrontCoordinatePolicy.WHOLE_CELL_OPEN_LOGIT
    ),
) -> CutIntegratorStep:
    """Solve and certify one strict same-cell receding step or roll back exactly.

    ``provisional_maximum_function_evaluations`` is an optional controller work
    cap for a trial leaf that may subsequently be bisected.  It cannot exceed
    the state's normal solver budget and changes no residual, chart, physical
    guard, or acceptance tolerance.  Calls that omit it retain the established
    production budget.

    ``nonlinear_solver_policy`` makes the sparse nonlinear method explicit
    while retaining the identical residual and acceptance contract.  The
    established route remains the default; sparse TRF/LSMR must be requested
    deliberately by a trajectory controller.  ``front_coordinate_policy``
    independently selects either the established whole-cell logit or the
    direct physical front coordinate.  Both decode to the same strict
    same-cell front primitive and share every physical acceptance gate.
    """

    evaluations = 0
    rejected_trials = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    last_candidate = None
    last_scaled: tuple[float, ...] = ()
    last_blocks = None
    minimum_bound_distance = 0.0
    optimizer_function_evaluations = 0
    provisional_work_budget_exhausted = False
    scalar_post_root_polish_audit: ScalarPostRootPolishAudit | None = None
    try:
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("time step must be positive and finite")
        if type(nonlinear_solver_policy) is not SameCellNonlinearSolverPolicy:
            raise TypeError("same-cell nonlinear solver policy must use the typed enum")
        if type(front_coordinate_policy) is not SameCellFrontCoordinatePolicy:
            raise TypeError("same-cell front-coordinate policy must use the typed enum")
        if provisional_maximum_function_evaluations is None:
            maximum_function_evaluations = before.controls.maximum_function_evaluations
        else:
            try:
                maximum_function_evaluations = operator.index(
                    provisional_maximum_function_evaluations
                )
            except TypeError as exc:
                raise ValueError("provisional nonlinear work budget must be an integer") from exc
            if (
                isinstance(provisional_maximum_function_evaluations, bool)
                or maximum_function_evaluations < 1
                or maximum_function_evaluations > before.controls.maximum_function_evaluations
            ):
                raise ValueError(
                    "provisional nonlinear work budget must lie in [1, "
                    f"{before.controls.maximum_function_evaluations}]"
                )
        config = before.transport.config if target_config is None else target_config
        cut.validate_pressure_only_target_config(before.transport.config, config)
        certify_explicit_wet_retained_water_applicability(
            before.controls,
            before.transport.wet_retained_water_loadings,
            before.transport.effective_wet_retained_water_capacity_duals_over_rt,
            before.transport.config.wet.luikov,
        )
        chart = _coordinate_chart(
            before,
            boundary,
            target_config=config,
            front_coordinate_policy=front_coordinate_policy,
        )
        _validate_seed(before, seed, boundary, chart=chart)
        jacobian_route_reason = _jacobian_route_reason(before, boundary)
        use_sparse = jacobian_route_reason != "small_rank_dense_oracle"
        if nonlinear_solver_policy is SameCellNonlinearSolverPolicy.SPARSE_TRF_LSMR:
            if not use_sparse:
                raise ValueError(
                    "sparse TRF/LSMR policy requires the declared structural sparse Jacobian"
                )
            nonlinear_method = "trf"
        else:
            nonlinear_method = "dogbox" if use_sparse else "trf"
        sparsity = (
            cut_sparsity.build_cut_jacobian_sparsity(
                before.transport.layout,
                face_conditioned_composition=(chart.dry_y_hexane is None),
                face_limit_regularized=True,
            )
            if use_sparse
            else None
        )
        initial_coordinates = _encode_candidate_for_chart(
            seed.candidate,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        initial_assembly = cut.assemble_backward_euler(
            before.transport,
            seed.candidate,
            dt_s,
            boundary,
            target_config=config,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        certification_scales = _residual_scales(
            before.transport,
            initial_assembly,
            dt_s,
        )
        optimization_scales = (
            _face_limit_regularized_residual_scales(
                before.transport,
                initial_assembly,
                dt_s,
            )
            if use_sparse
            else certification_scales
        )
        certification_scale_vector = np.asarray(
            certification_scales.vector,
            dtype=float,
        )
        optimization_scale_vector = np.asarray(
            optimization_scales.vector,
            dtype=float,
        )

        def assemble_scaled_residual_bases(
            coordinates: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray, cut.CutTransportAssembly]:
            nonlocal evaluations, rejected_trials
            candidate = _decode_candidate_for_chart(
                coordinates,
                before.transport.layout,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
            )
            try:
                assembly = cut.assemble_backward_euler(
                    before.transport,
                    candidate,
                    dt_s,
                    boundary,
                    target_config=config,
                    enforce_reduced_film_thresholds=False,
                )
            except cut.CutTransportStepError:
                rejected_trials += 1
                raise
            evaluations += 1
            original = _datum_covariant_residual_vector(assembly.residuals)
            optimization = (
                _face_limit_regularized_residual_vector(
                    original,
                    before.transport.layout,
                )
                if use_sparse
                else original
            )
            return (
                np.asarray(optimization, dtype=float) / optimization_scale_vector,
                np.asarray(original, dtype=float) / certification_scale_vector,
                assembly,
            )

        def assemble_scaled_residual(
            coordinates: np.ndarray,
        ) -> tuple[np.ndarray, cut.CutTransportAssembly]:
            optimization, _, assembly = assemble_scaled_residual_bases(coordinates)
            return optimization, assembly

        def assemble_original_scaled_residual(
            coordinates: np.ndarray,
        ) -> tuple[np.ndarray, cut.CutTransportAssembly]:
            _, original, assembly = assemble_scaled_residual_bases(coordinates)
            return original, assembly

        def scaled_residual(coordinates: np.ndarray) -> np.ndarray:
            return assemble_scaled_residual(coordinates)[0]

        lower_coordinates, upper_coordinates = _coordinate_solver_bounds_for_chart(
            before.transport.layout,
            before.controls.maximum_logit_magnitude,
            chart,
        )
        nonlinear_solver_options = (
            {"tr_solver": "lsmr"}
            if nonlinear_solver_policy is SameCellNonlinearSolverPolicy.SPARSE_TRF_LSMR
            else {}
        )
        solution = optimize.least_squares(
            scaled_residual,
            initial_coordinates,
            # The established sparse route uses rectangular dogleg trust
            # regions.  The explicit TRF/LSMR alternative uses the identical
            # bounded transformed coordinates, exact unit-determinant row
            # basis, and original-row acceptance check.  These are nonlinear
            # algorithm choices, not model fallbacks.  The smaller dense route
            # retains the original rows as an independent oracle.
            method=nonlinear_method,
            bounds=(lower_coordinates, upper_coordinates),
            ftol=before.controls.nonlinear_step_tolerance,
            xtol=before.controls.nonlinear_step_tolerance,
            gtol=(
                None
                if nonlinear_solver_policy is SameCellNonlinearSolverPolicy.SPARSE_TRF_LSMR
                else before.controls.nonlinear_step_tolerance
            ),
            max_nfev=maximum_function_evaluations,
            x_scale="jac",
            jac_sparsity=None if sparsity is None else sparsity.jac_sparsity,
            # With a structural sparse Jacobian SciPy solves each trust-region
            # subproblem by LSMR.  Its generic defaults are too loose for this
            # conservation problem: the outer dogleg iteration can otherwise
            # stall above the declared residual contract even though a
            # well-conditioned root exists.  These are linear-algebra
            # tolerances only; they do not alter the nonlinear equations,
            # physical charts, or acceptance thresholds.
            tr_options=(
                {
                    "atol": SPARSE_LSMR_ABSOLUTE_TOLERANCE,
                    "btol": SPARSE_LSMR_RELATIVE_TOLERANCE,
                    "maxiter": SPARSE_LSMR_MAXIMUM_ITERATIONS,
                }
                if sparsity is not None
                else None
            ),
            **nonlinear_solver_options,
        )
        optimizer_function_evaluations = int(solution.nfev)
        provisional_work_budget_exhausted = bool(
            provisional_maximum_function_evaluations is not None
            and not solution.success
            and solution.status == 0
            and solution.nfev >= maximum_function_evaluations
        )
        solution_coordinates = np.asarray(solution.x, dtype=float)
        converged_candidate = _decode_candidate_for_chart(
            solution_coordinates,
            before.transport.layout,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        last_candidate = converged_candidate
        (
            optimization_signed_vector,
            original_signed_vector,
            converged,
        ) = assemble_scaled_residual_bases(solution_coordinates)
        optimization_scaled = tuple(float(abs(value)) for value in optimization_signed_vector)
        scaled = tuple(float(abs(value)) for value in original_signed_vector)
        last_scaled = tuple(float(value) for value in original_signed_vector)
        last_blocks = converged.residuals
        candidate_fractions = _bounded_fractions(converged_candidate, chart)
        minimum_bound_distance = min(min(value, 1.0 - value) for value in candidate_fractions)
        maximum_scaled = max(scaled)
        condition_proxy = _jacobian_condition_proxy(solution.jac)
        if not solution.success:
            raise RuntimeError(f"nonlinear root did not certify success: {solution.message}")
        if np.any(solution.active_mask != 0):
            raise RuntimeError(
                "nonlinear solution touched a numerical coordinate guard; open-chart "
                "root is not certified"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > before.controls.maximum_condition_proxy
        ):
            raise ConditionProxyUncertifiedError(
                condition_proxy=condition_proxy,
                maximum_condition_proxy=before.controls.maximum_condition_proxy,
            )
        if maximum_scaled > before.controls.nonlinear_residual_tolerance:
            near_miss_ceiling = scalar_post_root_near_miss_ceiling(
                before.controls.nonlinear_residual_tolerance
            )
            eligible_for_scalar_polish = (
                use_sparse
                and sparsity is not None
                and isinstance(boundary, ct.ReducedFilmPoreBoundary)
                and solution.status == 3
                and rejected_trials == 0
                and provisional_maximum_function_evaluations is None
                and maximum_scaled <= near_miss_ceiling
            )
            if not eligible_for_scalar_polish:
                raise RuntimeError(
                    "regional/RH residual contract failed outside the deterministic "
                    "scalar near-miss safeguard: "
                    f"{maximum_scaled:.3e} > "
                    f"{before.controls.nonlinear_residual_tolerance:.3e}"
                )
            (
                solution_coordinates,
                converged,
                last_scaled,
                condition_proxy,
                scalar_post_root_polish_audit,
            ) = _attempt_scalar_post_root_polish(
                before,
                chart,
                sparsity,
                solution_coordinates,
                last_scaled,
                lower_coordinates,
                upper_coordinates,
                assemble_original_scaled_residual,
            )
            converged_candidate = converged.candidate
            last_candidate = converged_candidate
            last_blocks = converged.residuals
            scaled = tuple(abs(value) for value in last_scaled)
            maximum_scaled = max(scaled)
            candidate_fractions = _bounded_fractions(converged_candidate, chart)
            minimum_bound_distance = min(min(value, 1.0 - value) for value in candidate_fractions)
            polished_original = _datum_covariant_residual_vector(converged.residuals)
            polished_optimization = (
                _face_limit_regularized_residual_vector(
                    polished_original,
                    before.transport.layout,
                )
                if use_sparse
                else polished_original
            )
            optimization_scaled = tuple(
                abs(residual / scale)
                for residual, scale in zip(
                    polished_optimization,
                    optimization_scales.vector,
                )
            )
        maximum_optimization_scaled = max(optimization_scaled)
        if maximum_optimization_scaled > before.controls.nonlinear_residual_tolerance:
            raise RuntimeError(
                "face-limit optimization residual contract failed after original-basis "
                "certification: "
                f"{maximum_optimization_scaled:.3e} > "
                f"{before.controls.nonlinear_residual_tolerance:.3e}"
            )
        # Only a numerically certified root is entitled to exercise the
        # caller-owned fast-mass reduction threshold.  Reassembly is constitutively
        # identical but returns an explicitly enforced audit in the accepted
        # assembly; arbitrary trust-region iterates remain non-enforcing.
        converged = cut.assemble_backward_euler(
            before.transport,
            converged_candidate,
            dt_s,
            boundary,
            target_config=config,
            enforce_reduced_film_thresholds=True,
        )
        evaluations += 1
        production_covariant_residuals = _datum_covariant_residual_vector(converged.residuals)
        last_scaled = tuple(
            residual / scale
            for residual, scale in zip(
                production_covariant_residuals,
                certification_scales.vector,
            )
        )
        scaled = tuple(abs(value) for value in last_scaled)
        maximum_scaled = max(scaled)
        last_blocks = converged.residuals
        if maximum_scaled > before.controls.nonlinear_residual_tolerance:
            raise RuntimeError(
                "production reassembly failed the unchanged full residual contract: "
                f"{maximum_scaled:.3e} > "
                f"{before.controls.nonlinear_residual_tolerance:.3e}"
            )
        production_optimization_residuals = (
            _face_limit_regularized_residual_vector(
                production_covariant_residuals,
                before.transport.layout,
            )
            if use_sparse
            else production_covariant_residuals
        )
        maximum_production_optimization_scaled = max(
            abs(residual / scale)
            for residual, scale in zip(
                production_optimization_residuals,
                optimization_scales.vector,
            )
        )
        if maximum_production_optimization_scaled > before.controls.nonlinear_residual_tolerance:
            raise RuntimeError(
                "production reassembly failed the face-limit optimization residual "
                "contract: "
                f"{maximum_production_optimization_scaled:.3e} > "
                f"{before.controls.nonlinear_residual_tolerance:.3e}"
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
                    "production scalar-polish candidate lacks a fully enforced film audit"
                )
        if scalar_post_root_polish_audit is not None:
            cap = converged.wet_retained_cap_certificate
            if not cap.exact_graph_without_tolerance or cap.active_piece_count != 0:
                raise RuntimeError(
                    "production scalar-polish candidate left the exact smooth cap graph"
                )
            candidate_fractions = _bounded_fractions(converged_candidate, chart)
            minimum_bound_distance = min(min(value, 1.0 - value) for value in candidate_fractions)
        step = _accept_step(
            before,
            converged,
            boundary,
            seed,
            certification_scales,
            optimization_scales,
            scaled,
            dt_s,
            evaluations,
            rejected_trials,
            str(solution.message),
            condition_proxy,
            chart,
            sparsity,
            jacobian_route_reason,
            nonlinear_method,
            scalar_post_root_polish_audit,
        )
        if (
            step.ledger.maximum_step_ledger_residual > before.controls.ledger_tolerance
            or step.ledger.maximum_cumulative_ledger_residual > before.controls.ledger_tolerance
        ):
            raise RuntimeError("accepted/cumulative conservation ledger exceeded 1e-10 contract")
        certify_explicit_wet_retained_water_applicability(
            before.controls,
            before.transport.wet_retained_water_loadings,
            before.transport.effective_wet_retained_water_capacity_duals_over_rt,
            before.transport.config.wet.luikov,
        )
        certify_explicit_wet_retained_water_applicability(
            before.controls,
            step.after.transport.wet_retained_water_loadings,
            step.after.transport.effective_wet_retained_water_capacity_duals_over_rt,
            step.after.transport.config.wet.luikov,
        )
        return step
    except CutIntegratorStepError:
        raise
    except Exception as exc:
        event_restart = isinstance(exc, cut.CutTransportStepError) and (exc.event_restart_required)
        applicability_error = exc if isinstance(exc, WetRetainedWaterApplicabilityError) else None
        raise CutIntegratorStepError(
            f"same-cell cut solve rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            event_restart_required=event_restart,
            last_candidate=last_candidate,
            scaled_residuals=last_scaled,
            residual_blocks=last_blocks,
            minimum_fractional_distance_to_bound=minimum_bound_distance,
            nonlinear_optimizer_function_evaluations=(optimizer_function_evaluations),
            provisional_work_budget=provisional_maximum_function_evaluations,
            provisional_work_budget_exhausted=provisional_work_budget_exhausted,
            retained_water_applicability_guard_failed=(applicability_error is not None),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


@dataclass(frozen=True)
class _CoordinateChart:
    wet_temperature: tuple[float, float]
    wet_water: tuple[float, float]
    wet_retained_cap: wrc.WetRetainedCapSemismoothChart | None
    dry_temperature: tuple[float, float]
    dry_y_hexane: tuple[float, float] | None
    dry_config: ct.FullyDryTransportConfig
    dry_pores: tuple[cp.CoupledPoreParams, ...]
    surface_boundary: ct.PoreBoundary | None
    front_z: tuple[float, float]
    interface_temperature: tuple[float, float]
    moving_interface_composition_force_authority: cut.MovingInterfaceCompositionForceAuthority
    front_coordinate_policy: SameCellFrontCoordinatePolicy = (
        SameCellFrontCoordinatePolicy.WHOLE_CELL_OPEN_LOGIT
    )


def _coordinate_chart(
    before: CutIntegratorState,
    surface_boundary: ct.PoreBoundary | None = None,
    *,
    target_config: cut.CutTransportConfig | None = None,
    front_coordinate_policy: SameCellFrontCoordinatePolicy = (
        SameCellFrontCoordinatePolicy.WHOLE_CELL_OPEN_LOGIT
    ),
) -> _CoordinateChart:
    config = before.transport.config if target_config is None else target_config
    cut.validate_pressure_only_target_config(before.transport.config, config)
    geometry = before.transport.geometry
    cut_index = before.transport.layout.cut_cell_index
    inner_radius = geometry.master_grid.faces[cut_index]
    inner_z = (inner_radius / geometry.master_grid.R) ** 3
    if not inner_z < geometry.front.z:
        raise cut.CutTransportTopologyError("current cut chart has no positive z width")
    if type(front_coordinate_policy) is not SameCellFrontCoordinatePolicy:
        raise TypeError("same-cell front-coordinate policy must use the typed enum")
    if surface_boundary is not None and not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("surface boundary must be Dirichlet or bounded reduced film")
    return _CoordinateChart(
        wet_temperature=before.controls.wet_temperature_bounds_k,
        wet_water=before.controls.wet_water_bounds,
        wet_retained_cap=(
            wrc.WetRetainedCapSemismoothChart(
                before.controls.wet_water_bounds[0],
                config.wet.luikov,
            )
            if before.controls.wet_water_bounds[1] == config.wet.luikov.W_cap
            else None
        ),
        dry_temperature=(config.dry.conditioned_temperature_domain.solver_bounds_k),
        dry_y_hexane=(
            config.dry.primitive_band.y_hexane_bounds
            if isinstance(
                config.dry.primitive_band,
                ct.JointPrimitiveBand,
            )
            else None
        ),
        dry_config=config.dry,
        dry_pores=tuple(
            replace(
                config.dry.pore,
                w_o=before.transport.oil_fraction_labels[cell],
            )
            for cell in before.transport.layout.dry_cell_indices
        ),
        surface_boundary=surface_boundary,
        front_z=(inner_z, geometry.front.z),
        interface_temperature=before.controls.interface_temperature_bounds_k,
        moving_interface_composition_force_authority=(
            config.moving_interface_composition_force_authority
        ),
        front_coordinate_policy=front_coordinate_policy,
    )


def _encode_candidate(
    candidate: cut.CutTransportUnknowns,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> np.ndarray:
    wet_duals = candidate.effective_wet_retained_water_capacity_duals_over_rt
    if chart.wet_retained_cap is None:
        if any(value != 0.0 for value in wet_duals):
            raise ValueError(
                "a positive wet retained-cap dual requires the exact-cap semismooth chart"
            )
        encoded_wet_water = tuple(
            _encode_open(value, chart.wet_water) for value in candidate.wet_retained_water_loadings
        )
    else:
        encoded_wet_water = tuple(
            chart.wet_retained_cap.encode(loading, dual)
            for loading, dual in zip(
                candidate.wet_retained_water_loadings,
                wet_duals,
            )
        )
    values = (
        *(_encode_open(value, chart.wet_temperature) for value in candidate.wet_temperatures_k),
        *encoded_wet_water,
        *(_encode_open(value, chart.dry_temperature) for value in candidate.dry_temperatures_k),
        *(
            _encode_dry_y(
                piece,
                candidate.dry_temperatures_k,
                candidate.interface_temperature_k,
                value,
                chart,
            )
            for piece, value in enumerate(candidate.dry_y_hexane)
        ),
        *(
            value / scale
            for value, scale in zip(
                candidate.dry_total_stefan_fluxes_mol_m2_s,
                stefan_scales,
            )
        ),
        _encode_open(candidate.front_z, chart.front_z),
        _encode_open(candidate.interface_temperature_k, chart.interface_temperature),
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("encoded nonlinear coordinates must be finite")
    return np.asarray(values, dtype=float)


def _decode_candidate(
    coordinates: Sequence[float],
    layout: cut.CutTransportLayout,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> cut.CutTransportUnknowns:
    values = tuple(float(value) for value in coordinates)
    if len(values) != layout.unknown_count or not all(math.isfinite(value) for value in values):
        raise ValueError("nonlinear coordinate vector is non-finite or has wrong rank")
    nw = layout.wet_piece_count
    nd = layout.dry_piece_count
    nf = layout.dry_face_count
    cursor = 0

    def take(count: int) -> tuple[float, ...]:
        nonlocal cursor
        result = values[cursor : cursor + count]
        cursor += count
        return result

    wet_t = tuple(_decode_open(value, chart.wet_temperature) for value in take(nw))
    wet_water_coordinates = take(nw)
    if chart.wet_retained_cap is None:
        wet_w = tuple(_decode_open(value, chart.wet_water) for value in wet_water_coordinates)
        wet_duals = (0.0,) * nw
    else:
        wet_graph = tuple(chart.wet_retained_cap.decode(value) for value in wet_water_coordinates)
        wet_w = tuple(point.retained_water_loading for point in wet_graph)
        wet_duals = tuple(point.capacity_dual_over_rt for point in wet_graph)
    dry_t = tuple(_decode_open(value, chart.dry_temperature) for value in take(nd))
    dry_y_coordinates = take(nd)
    nt_coordinates = take(nf)
    front_coordinate = take(1)[0]
    interface_coordinate = take(1)[0]
    interface_t = _decode_open(interface_coordinate, chart.interface_temperature)
    dry_y = tuple(
        _decode_dry_y(
            piece,
            dry_t,
            interface_t,
            value,
            chart,
        )
        for piece, value in enumerate(dry_y_coordinates)
    )
    nt = tuple(value * scale for value, scale in zip(nt_coordinates, stefan_scales))
    front_z = _decode_open(front_coordinate, chart.front_z)
    return cut.CutTransportUnknowns(
        wet_temperatures_k=wet_t,
        wet_retained_water_loadings=wet_w,
        dry_temperatures_k=dry_t,
        dry_y_hexane=dry_y,
        dry_total_stefan_fluxes_mol_m2_s=nt,
        front_z=front_z,
        interface_temperature_k=interface_t,
        wet_retained_water_capacity_duals_over_rt=wet_duals,
    )


def _encode_direct_physical_front_z(
    front_z: float,
    bounds: tuple[float, float],
) -> float:
    """Return the strict physical front coordinate without arithmetic."""

    if not _inside(front_z, bounds):
        raise ValueError("cannot encode a direct physical front outside its open interval")
    return front_z


def _decode_direct_physical_front_z(
    front_z: float,
    bounds: tuple[float, float],
) -> float:
    """Return the binary64 coordinate as physical ``front_z`` exactly."""

    if not _inside(front_z, bounds):
        raise ValueError("direct physical front is outside its open interval")
    return front_z


def _encode_candidate_for_chart(
    candidate: cut.CutTransportUnknowns,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> np.ndarray:
    """Encode with the typed front policy; legacy mapping stays bitwise intact."""

    if chart.front_coordinate_policy is SameCellFrontCoordinatePolicy.WHOLE_CELL_OPEN_LOGIT:
        return _encode_candidate(candidate, chart, stefan_scales)
    if chart.front_coordinate_policy is not SameCellFrontCoordinatePolicy.DIRECT_PHYSICAL_FRONT_Z:
        raise ValueError("unknown same-cell front-coordinate policy")
    lower, upper = chart.front_z
    midpoint = lower + 0.5 * (upper - lower)
    coordinates = _encode_candidate(
        replace(candidate, front_z=midpoint),
        chart,
        stefan_scales,
    )
    coordinates[-2] = _encode_direct_physical_front_z(
        candidate.front_z,
        chart.front_z,
    )
    return coordinates


def _decode_candidate_for_chart(
    coordinates: Sequence[float],
    layout: cut.CutTransportLayout,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> cut.CutTransportUnknowns:
    """Decode with the typed front policy and no physical-state fallback."""

    if chart.front_coordinate_policy is SameCellFrontCoordinatePolicy.WHOLE_CELL_OPEN_LOGIT:
        return _decode_candidate(coordinates, layout, chart, stefan_scales)
    if chart.front_coordinate_policy is not SameCellFrontCoordinatePolicy.DIRECT_PHYSICAL_FRONT_Z:
        raise ValueError("unknown same-cell front-coordinate policy")
    values = np.asarray(tuple(float(value) for value in coordinates), dtype=float)
    if values.shape != (layout.unknown_count,) or not np.all(np.isfinite(values)):
        raise ValueError("nonlinear coordinate vector is non-finite or has wrong rank")
    front_z = _decode_direct_physical_front_z(
        float(values[-2]),
        chart.front_z,
    )
    legacy_coordinates = values.copy()
    legacy_coordinates[-2] = 0.0
    decoded = _decode_candidate(
        legacy_coordinates,
        layout,
        chart,
        stefan_scales,
    )
    return replace(decoded, front_z=front_z)


def _encode_open(value: float, bounds: tuple[float, float]) -> float:
    lower, upper = bounds
    if not _inside(value, bounds):
        raise ValueError("cannot encode a bounded primitive on/outside its open interval")
    fraction = (value - lower) / (upper - lower)
    return math.log(fraction) - math.log1p(-fraction)


def _decode_open(value: float, bounds: tuple[float, float]) -> float:
    lower, upper = bounds
    decoded = lower + (upper - lower) * float(expit(value))
    if not _inside(decoded, bounds):
        raise ValueError("transformed coordinate collapsed onto an open-interval boundary")
    return decoded


def _inside(value: float, bounds: tuple[float, float]) -> bool:
    return math.isfinite(value) and bounds[0] < value < bounds[1]


def _encode_dry_y(
    piece: int,
    dry_temperatures_k: Sequence[float],
    interface_temperature_k: float,
    y_hexane: float,
    chart: _CoordinateChart,
) -> float:
    if chart.dry_y_hexane is not None:
        return _encode_open(y_hexane, chart.dry_y_hexane)
    temperatures = _dry_composition_context_temperatures(
        piece,
        dry_temperatures_k,
        interface_temperature_k,
        chart,
    )
    return _encode_face_conditioned_dry_y(
        temperatures,
        chart.dry_config.pressure_pa,
        y_hexane,
        chart.dry_pores[piece],
    )


def _decode_dry_y(
    piece: int,
    dry_temperatures_k: Sequence[float],
    interface_temperature_k: float,
    coordinate: float,
    chart: _CoordinateChart,
) -> float:
    if chart.dry_y_hexane is not None:
        return _decode_open(coordinate, chart.dry_y_hexane)
    temperatures = _dry_composition_context_temperatures(
        piece,
        dry_temperatures_k,
        interface_temperature_k,
        chart,
    )
    return _decode_face_conditioned_dry_y(
        temperatures,
        chart.dry_config.pressure_pa,
        coordinate,
        chart.dry_pores[piece],
    )


def _dry_composition_context_temperatures(
    piece: int,
    dry_temperatures_k: Sequence[float],
    interface_temperature_k: float,
    chart: _CoordinateChart,
) -> tuple[float, ...]:
    """Return the cell and every adjacent face temperature for one ``y_i``.

    ``cut_transport._dry_flux`` evaluates each cell composition at its own
    storage temperature and at every adjacent common-temperature face.  The
    legacy moving RH face also contributes ``T_Gamma``.  Under Amendment 20's
    explicit actual-path authority it does not: the moving path is certified
    by its own residual operator, while internal and surface faces retain
    their symmetric endpoint temperatures.  The outer face is included when
    the actual boundary is available.
    """

    temperatures = tuple(float(value) for value in dry_temperatures_k)
    if len(temperatures) != len(chart.dry_pores):
        raise ValueError("dry temperatures do not align with exact per-cell labels")
    if not 0 <= piece < len(temperatures):
        raise IndexError("dry piece is outside the face-conditioned chart")
    cell_temperature = temperatures[piece]
    result = [cell_temperature]
    if piece == 0:
        if chart.moving_interface_composition_force_authority is (
            cut.MovingInterfaceCompositionForceAuthority.LEGACY_COMMON_T_GAMMA
        ):
            result.append(interface_temperature_k)
    else:
        result.append(0.5 * (temperatures[piece - 1] + cell_temperature))
    if piece + 1 < len(temperatures):
        right_temperature = temperatures[piece + 1]
        result.append(0.5 * (cell_temperature + right_temperature))
    elif chart.surface_boundary is not None:
        result.append(0.5 * (cell_temperature + chart.surface_boundary.temperature_k))
    if not all(math.isfinite(value) and value > 0.0 for value in result):
        raise cp.CoupledPoreTopologyError(
            "face-conditioned dry temperatures must be positive and finite"
        )
    return tuple(result)


def _face_conditioned_composition_interval(
    temperatures_k: Sequence[float],
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
) -> tuple[float, float]:
    """Intersect the exact open gas-only intervals used by one cell/its faces."""

    temperatures = tuple(float(value) for value in temperatures_k)
    if not temperatures:
        raise ValueError("face-conditioned composition chart needs a temperature")
    intervals = tuple(
        cp.gas_only_composition_interval(temperature, pressure_pa, pore)
        for temperature in temperatures
    )
    lower = max(interval.lower_y_hexane for interval in intervals)
    upper = min(interval.upper_y_hexane for interval in intervals)
    if not lower < upper:
        raise cp.CoupledPoreTopologyError(
            "cell and adjacent symmetric faces have no common open gas-only "
            "composition interval; another phase/topology is required"
        )
    return lower, upper


def _face_conditioned_coordinate_and_derivative(
    temperatures_k: Sequence[float],
    pressure_pa: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
) -> tuple[float, float]:
    r"""Return the average local barrier and its ``d/dy_h``.

    Combining all barriers (rather than selecting the currently tightest
    interval endpoint)
    preserves a smooth, strictly increasing map when the active phase boundary
    changes with temperature.  At either endpoint of the common interval at
    least one constituent barrier diverges, so the composite map still covers
    all of :math:`(-\infty,\infty)` one-to-one in real arithmetic.  Dividing
    by the exact stencil count keeps the raw chart independent of whether a
    piece has two or three adjacent evaluation temperatures.

    A context whose representable interior cannot bracket a requested
    coordinate fails closed in :func:`_decode_face_conditioned_dry_y`.
    """

    temperatures = tuple(float(value) for value in temperatures_k)
    if not temperatures:
        raise ValueError("face-conditioned composition chart needs a temperature")
    activity_at_ref = sp.water_activity(pore.luikov.W_ref, pore.luikov)
    coordinates = []
    derivatives = []
    for temperature in temperatures:
        values = cp._gas_component_activities_with_derivatives(  # noqa: SLF001
            temperature,
            pressure_pa,
            y_hexane,
            pore,
        )
        coordinates.append(
            cp._gas_only_barrier_coordinate(  # noqa: SLF001
                y_hexane,
                values[0],
                values[1],
                activity_at_ref,
            )
        )
        derivatives.append(
            cp._gas_only_barrier_coordinate_derivative(  # noqa: SLF001
                y_hexane,
                values[0],
                values[1],
                values[2],
                values[3],
                activity_at_ref,
            )
        )
    count = len(temperatures)
    coordinate = math.fsum(coordinates) / count
    derivative = math.fsum(derivatives) / count
    if not math.isfinite(coordinate) or not math.isfinite(derivative) or derivative <= 0.0:
        raise cp.CoupledPoreTopologyError(
            "face-conditioned gas-only barrier lost finite strict monotonicity"
        )
    return coordinate, derivative


def _encode_face_conditioned_dry_y(
    temperatures_k: Sequence[float],
    pressure_pa: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
) -> float:
    """Encode one composition that is strictly feasible in its whole stencil."""

    lower, upper = _face_conditioned_composition_interval(
        temperatures_k,
        pressure_pa,
        pore,
    )
    if not math.isfinite(y_hexane) or not lower < y_hexane < upper:
        raise cp.CoupledPoreTopologyError(
            "dry composition is not strictly feasible at its cell and both "
            "adjacent symmetric face temperatures; do not clip"
        )
    return _face_conditioned_coordinate_and_derivative(
        temperatures_k,
        pressure_pa,
        y_hexane,
        pore,
    )[0]


def _decode_face_conditioned_dry_y(
    temperatures_k: Sequence[float],
    pressure_pa: float,
    coordinate: float,
    pore: cp.CoupledPoreParams,
) -> float:
    """Safeguarded inverse of the exact all-cell/all-face barrier."""

    if not math.isfinite(coordinate):
        raise ValueError("face-conditioned composition coordinate must be finite")
    lower, upper = _face_conditioned_composition_interval(
        temperatures_k,
        pressure_pa,
        pore,
    )
    composition = lower + 0.5 * (upper - lower)
    best_composition = composition
    best_absolute_residual = math.inf
    saw_negative_residual = False
    saw_positive_residual = False
    last_composition = math.nan
    for _ in range(128):
        if not lower < composition < upper or composition == last_composition:
            break
        last_composition = composition
        encoded, derivative = _face_conditioned_coordinate_and_derivative(
            temperatures_k,
            pressure_pa,
            composition,
            pore,
        )
        residual = encoded - coordinate
        if abs(residual) < best_absolute_residual:
            best_composition = composition
            best_absolute_residual = abs(residual)
        coordinate_tolerance = 64.0 * math.ulp(max(1.0, abs(coordinate)))
        if abs(residual) <= coordinate_tolerance:
            return composition
        if residual < 0.0:
            saw_negative_residual = True
            lower = composition
        else:
            saw_positive_residual = True
            upper = composition
        newton_composition = composition - residual / derivative
        if lower < newton_composition < upper:
            next_composition = newton_composition
        else:
            next_composition = lower + 0.5 * (upper - lower)
        if next_composition == composition:
            next_composition = lower + 0.5 * (upper - lower)
        if not lower < next_composition < upper:
            break
        composition = next_composition
    if saw_negative_residual and saw_positive_residual:
        # The target lies between adjacent representable interior values.
        # Nearest directed rounding is not a phase-boundary projection.
        return best_composition
    raise cp.CoupledPoreTopologyError(
        "face-conditioned gas-only coordinate has no converged representable "
        "interior value; do not clip"
    )


def _validate_controls_against_state(
    state: cut.CutTransportState, controls: CutSolverControls
) -> None:
    wet_t_lo, wet_t_hi = controls.wet_temperature_bounds_k
    if not state.config.wet.wet.T_min <= wet_t_lo < wet_t_hi <= state.config.wet.wet.T_max:
        raise ValueError("wet solver chart must lie inside the wet caloric authority")
    water_lo, water_hi = controls.wet_water_bounds
    if not state.config.wet.luikov.W_ref <= water_lo < water_hi <= state.config.wet.luikov.W_cap:
        raise ValueError("wet-water solver chart must lie inside the Luikov authority")
    gamma_lo, gamma_hi = controls.interface_temperature_bounds_k
    pore_lo, pore_hi = state.config.dry.pore.temperature_bounds_k
    if (
        not max(state.config.wet.wet.T_min, pore_lo)
        <= gamma_lo
        < gamma_hi
        <= min(
            state.config.wet.wet.T_max,
            pore_hi,
        )
    ):
        raise ValueError("interface-temperature chart left its shared property authority")
    for value in state.wet_temperatures_k:
        if not _inside(value, controls.wet_temperature_bounds_k):
            raise ValueError("wet state temperature is not strictly inside its solver chart")
    wet_duals = state.effective_wet_retained_water_capacity_duals_over_rt
    wrc.certify_exact_graph(
        state.wet_retained_water_loadings,
        wet_duals,
        state.config.wet.luikov,
    )
    cap_chart_enabled = water_hi == state.config.wet.luikov.W_cap
    for value, dual in zip(state.wet_retained_water_loadings, wet_duals):
        if cap_chart_enabled:
            if not water_lo < value <= water_hi:
                raise ValueError(
                    "wet state water loading is outside its semiclosed cap solver chart"
                )
        elif not _inside(value, controls.wet_water_bounds) or dual != 0.0:
            raise ValueError(
                "wet state water loading is not strictly inside its smooth solver chart"
            )
    certify_explicit_wet_retained_water_applicability(
        controls,
        state.wet_retained_water_loadings,
        wet_duals,
        state.config.wet.luikov,
    )
    dry_t_bounds = state.config.dry.conditioned_temperature_domain.solver_bounds_k
    for value in state.dry_temperatures_k:
        if not _inside(value, dry_t_bounds):
            raise ValueError("dry state temperature is not strictly inside its solver chart")
    for cell, temperature, value in zip(
        state.layout.dry_cell_indices,
        state.dry_temperatures_k,
        state.dry_y_hexane,
    ):
        cut._validate_open_dry_band(  # noqa: SLF001
            temperature,
            value,
            state.config,
            pore=replace(
                state.config.dry.pore,
                w_o=state.oil_fraction_labels[cell],
            ),
        )
    cut_index = state.layout.cut_cell_index
    for temperature in (
        gamma_lo,
        0.5 * (gamma_lo + gamma_hi),
        gamma_hi,
    ):
        cut.evaluate_interface_state(
            temperature,
            state.config,
            state.historical_hexane_loadings[cut_index],
            state.oil_fraction_labels[cut_index],
        )


def _validate_seed(
    before: CutIntegratorState,
    seed: CutStepSeed,
    surface_boundary: ct.PoreBoundary | None = None,
    *,
    chart: _CoordinateChart | None = None,
) -> None:
    layout = before.transport.layout
    if len(seed.candidate.equation_rank_vector()) != layout.unknown_count:
        raise ValueError("seed candidate does not match the current nonlinear rank")
    if len(seed.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
        raise ValueError("one explicit Stefan scale is required per dry face")
    active_chart = _coordinate_chart(before, surface_boundary) if chart is None else chart
    coordinates = _encode_candidate_for_chart(
        seed.candidate,
        active_chart,
        seed.stefan_flux_scales_mol_m2_s,
    )
    lower, upper = _coordinate_solver_bounds_for_chart(
        layout,
        before.controls.maximum_logit_magnitude,
        active_chart,
    )
    outside = np.flatnonzero((coordinates < lower) | (coordinates > upper))
    if outside.size:
        indices = ", ".join(str(int(index)) for index in outside)
        raise ValueError(
            "encoded seed is outside the existing solver-coordinate bounds "
            f"at indices {indices}; do not clip or project"
        )


def _piece_data(
    state: cut.CutTransportState,
    *,
    evaluate_capacity: bool = True,
) -> _PieceData:
    geometry = state.geometry
    layout = state.layout
    wet_cells = tuple(
        ww.evaluate_cell(
            temperature,
            water,
            state.historical_hexane_loadings[cell],
            state.oil_fraction_labels[cell],
            state.config.wet,
        )
        for cell, temperature, water in zip(
            layout.wet_cell_indices,
            state.wet_temperatures_k,
            state.wet_retained_water_loadings,
        )
    )
    dry_cells = tuple(
        cp.evaluate_equilibrium(
            temperature,
            state.config.dry.pressure_pa,
            y_hexane,
            replace(
                state.config.dry.pore,
                w_o=state.oil_fraction_labels[cell],
            ),
        )
        for cell, temperature, y_hexane in zip(
            layout.dry_cell_indices,
            state.dry_temperatures_k,
            state.dry_y_hexane,
        )
    )
    wet_volumes = tuple(geometry.cells[cell].wet_volume_m3 for cell in layout.wet_cell_indices)
    dry_volumes = tuple(geometry.cells[cell].dry_volume_m3 for cell in layout.dry_cell_indices)
    wet_water = tuple(
        volume * cell.retained_water_concentration_mol_m3
        for volume, cell in zip(wet_volumes, wet_cells)
    )
    wet_hexane = tuple(
        volume * state.config.wet.wet.rho_dm_p * state.historical_hexane_loadings[index] / hx.M
        for volume, index in zip(wet_volumes, layout.wet_cell_indices)
    )
    wet_energy = tuple(
        volume * cell.energy_density_j_m3 for volume, cell in zip(wet_volumes, wet_cells)
    )
    dry_water = tuple(
        volume * cell.total_water_concentration_mol_m3
        for volume, cell in zip(dry_volumes, dry_cells)
    )
    dry_hexane = tuple(
        volume * cell.total_hexane_concentration_mol_m3
        for volume, cell in zip(dry_volumes, dry_cells)
    )
    dry_energy = tuple(
        volume * cell.energy_density_j_m3 for volume, cell in zip(dry_volumes, dry_cells)
    )
    if evaluate_capacity:
        wet_capacity = tuple(
            volume
            * wet_core.energy_capacity(
                cell.temperature_k,
                cell.historical_hexane_loading,
                replace(
                    state.config.wet.wet,
                    X_water=cell.retained_water_loading,
                    w_o=cell.oil_fraction_label,
                ),
            )
            * max(abs(cell.temperature_k), 1.0)
            for volume, cell in zip(wet_volumes, wet_cells)
        )
        dry_capacity_results = tuple(
            _dry_energy_capacity_with_audit(
                cell.temperature_k,
                cell.y_hexane,
                replace(
                    state.config.dry.pore,
                    w_o=state.oil_fraction_labels[index],
                ),
                state.config,
            )
            for cell, index in zip(
                dry_cells,
                layout.dry_cell_indices,
            )
        )
        dry_capacity = tuple(
            volume * result[0] * max(abs(cell.temperature_k), 1.0)
            for volume, cell, result in zip(
                dry_volumes,
                dry_cells,
                dry_capacity_results,
            )
        )
        dry_capacity_audits = tuple(result[1] for result in dry_capacity_results)
    else:
        wet_capacity = ()
        dry_capacity = ()
        dry_capacity_audits = ()
    return _PieceData(
        wet_cells,
        dry_cells,
        wet_water,
        wet_hexane,
        wet_energy,
        dry_water,
        dry_hexane,
        dry_energy,
        wet_capacity,
        dry_capacity,
        dry_capacity_audits,
    )


def _dry_energy_capacity(
    temperature: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
    config: cut.CutTransportConfig,
) -> float:
    """Compatibility scalar view of the audited Amendment-19 stencil."""

    return _dry_energy_capacity_with_audit(
        temperature,
        y_hexane,
        pore,
        config,
    )[0]


def _dry_energy_capacity_with_audit(
    temperature: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
    config: cut.CutTransportConfig,
) -> tuple[float, DryEnergyCapacityStencilAudit]:
    """Return one datum-reduced dry capacity and its topology certificate.

    ``JointPrimitiveBand`` retains the historical centered arithmetic exactly.
    Only the exact ``GasOnlyPrimitiveDomain`` may choose a second-order
    one-sided stencil when a centered artificial sample leaves that topology.
    The helper changes normalization only; no returned sample enters physical
    storage, residual assembly, or a conservation ledger.
    """

    t_lo, t_hi = config.dry.conditioned_temperature_domain.solver_bounds_k
    delta = min(1.0e-3, 0.25 * (temperature - t_lo), 0.25 * (t_hi - temperature))
    machine_step = math.sqrt(np.finfo(float).eps) * max(abs(temperature), 1.0)
    if delta <= machine_step:
        raise ValueError("dry capacity stencil is too near its open temperature bound")

    def reduced_energy(state: cp.EquilibriumPoreState) -> float:
        return math.fsum(
            (
                state.energy_density_j_m3,
                -hx.H_REF * state.total_hexane_concentration_mol_m3,
            )
        )

    primitive_domain = config.dry.primitive_band
    if isinstance(primitive_domain, ct.JointPrimitiveBand):
        # Preserve the historical evaluation order and arithmetic bit for bit.
        hotter = cp.evaluate_equilibrium(
            temperature + delta,
            config.dry.pressure_pa,
            y_hexane,
            pore,
        )
        colder = cp.evaluate_equilibrium(
            temperature - delta,
            config.dry.pressure_pa,
            y_hexane,
            pore,
        )
        hotter_reduced_energy = reduced_energy(hotter)
        colder_reduced_energy = reduced_energy(colder)
        capacity = (hotter_reduced_energy - colder_reduced_energy) / (2.0 * delta)
        if not math.isfinite(capacity) or capacity <= 0.0:
            raise ValueError("dry energy capacity must be positive and finite")
        return capacity, DryEnergyCapacityStencilAudit(
            temperature_k=temperature,
            pressure_pa=config.dry.pressure_pa,
            y_hexane=y_hexane,
            actual_oil_fraction_label=pore.w_o,
            gas_accessible_fraction=1.0,
            primitive_domain_kind="JointPrimitiveBand",
            stencil_kind="joint_primitive_legacy_centered_unchanged",
            delta_temperature_k=delta,
            sample_temperatures_k=(temperature + delta, temperature - delta),
            sample_full_coupled_pore_certified=(True, True),
            center_gas_only_admissible=None,
            central_lower_gas_only_admissible=None,
            central_upper_gas_only_admissible=None,
            forward_second_gas_only_admissible=None,
            backward_second_gas_only_admissible=None,
            capacity_j_m3_k=capacity,
        )
    if not isinstance(primitive_domain, ct.GasOnlyPrimitiveDomain):
        raise TypeError("dry capacity stencil requires a known primitive domain")

    sample_cache: dict[
        float,
        tuple[bool, cp.EquilibriumPoreState | None, float | None],
    ] = {}

    def gas_only_sample(
        sample_temperature: float,
    ) -> tuple[bool, cp.EquilibriumPoreState | None, float | None]:
        cached = sample_cache.get(sample_temperature)
        if cached is not None:
            return cached
        if not t_lo < sample_temperature < t_hi:
            result = (False, None, None)
        else:
            try:
                cp.encode_gas_only_y(
                    sample_temperature,
                    config.dry.pressure_pa,
                    y_hexane,
                    pore,
                )
                state = cp.evaluate_equilibrium(
                    sample_temperature,
                    config.dry.pressure_pa,
                    y_hexane,
                    pore,
                )
            except cp.CoupledPoreTopologyError:
                result = (False, None, None)
            else:
                result = (True, state, reduced_energy(state))
        sample_cache[sample_temperature] = result
        return result

    center_ok, center_state, center_energy = gas_only_sample(temperature)
    if not center_ok or center_state is None or center_energy is None:
        raise ValueError("dry capacity center is outside the exact gas-only topology")

    colder_temperature = temperature - delta
    hotter_temperature = temperature + delta
    hotter_second_temperature = temperature + 2.0 * delta
    colder_second_temperature = temperature - 2.0 * delta
    colder_ok, colder, colder_energy = gas_only_sample(colder_temperature)
    hotter_ok, hotter, hotter_energy = gas_only_sample(hotter_temperature)
    forward_second_ok: bool | None = None
    backward_second_ok: bool | None = None

    if (
        colder_ok
        and hotter_ok
        and colder is not None
        and hotter is not None
        and colder_energy is not None
        and hotter_energy is not None
    ):
        capacity = (hotter_energy - colder_energy) / (2.0 * delta)
        stencil_kind = "gas_only_second_order_centered"
        sample_temperatures = (hotter_temperature, colder_temperature)
        sample_certified = (True, True)
    else:
        forward_second_ok, forward_second, forward_second_energy = gas_only_sample(
            hotter_second_temperature
        )
        if (
            hotter_ok
            and forward_second_ok
            and hotter is not None
            and hotter_energy is not None
            and forward_second is not None
            and forward_second_energy is not None
        ):
            capacity = math.fsum(
                (-3.0 * center_energy, 4.0 * hotter_energy, -forward_second_energy)
            ) / (2.0 * delta)
            stencil_kind = "gas_only_second_order_forward"
            sample_temperatures = (
                temperature,
                hotter_temperature,
                hotter_second_temperature,
            )
            sample_certified = (True, True, True)
        else:
            backward_second_ok, backward_second, backward_second_energy = gas_only_sample(
                colder_second_temperature
            )
            if not (
                colder_ok
                and backward_second_ok
                and colder is not None
                and colder_energy is not None
                and backward_second is not None
                and backward_second_energy is not None
            ):
                raise ValueError("dry capacity stencil has no resolved open direction")
            capacity = math.fsum(
                (3.0 * center_energy, -4.0 * colder_energy, backward_second_energy)
            ) / (2.0 * delta)
            stencil_kind = "gas_only_second_order_backward"
            sample_temperatures = (
                temperature,
                colder_temperature,
                colder_second_temperature,
            )
            sample_certified = (True, True, True)

    if not math.isfinite(capacity) or capacity <= 0.0:
        raise ValueError("dry energy capacity must be positive and finite")
    return capacity, DryEnergyCapacityStencilAudit(
        temperature_k=temperature,
        pressure_pa=config.dry.pressure_pa,
        y_hexane=y_hexane,
        actual_oil_fraction_label=pore.w_o,
        gas_accessible_fraction=1.0,
        primitive_domain_kind="GasOnlyPrimitiveDomain",
        stencil_kind=stencil_kind,
        delta_temperature_k=delta,
        sample_temperatures_k=sample_temperatures,
        sample_full_coupled_pore_certified=sample_certified,
        center_gas_only_admissible=True,
        central_lower_gas_only_admissible=colder_ok,
        central_upper_gas_only_admissible=hotter_ok,
        forward_second_gas_only_admissible=forward_second_ok,
        backward_second_gas_only_admissible=backward_second_ok,
        capacity_j_m3_k=capacity,
    )


def _datum_covariant_residual_vector(
    residuals: cut.CutResidualBlocks,
) -> tuple[float, ...]:
    """Return an equation-equivalent residual basis independent of ``H_REF``.

    Shifting every hexane molar enthalpy by a constant ``a`` transforms each
    dry or interface energy balance as ``R_E -> R_E + a*R_h``.  Subtracting
    ``H_REF*R_h`` is therefore an invertible row operation on the simultaneous
    component/energy system, not a new closure or a deleted energy term.
    Wet-shell hexane is a frozen material label, so its datum contribution
    cancels identically inside each wet energy balance and has no separate row.
    """

    dry_energy = tuple(
        math.fsum((energy, -hx.H_REF * hexane))
        for energy, hexane in zip(
            residuals.dry_energy_w,
            residuals.dry_hexane_mol_s,
        )
    )
    rh_energy = math.fsum(
        (
            residuals.rh_energy_w,
            -hx.H_REF * residuals.rh_hexane_mol_s,
        )
    )
    return (
        *residuals.wet_water_mol_s,
        *residuals.wet_energy_w,
        *residuals.dry_water_mol_s,
        *residuals.dry_hexane_mol_s,
        *dry_energy,
        residuals.rh_water_mol_s,
        residuals.rh_hexane_mol_s,
        rh_energy,
    )


def _face_limit_residual_row_indices(
    layout: cut.CutTransportLayout,
) -> tuple[int, int, int, int, int, int, int, int]:
    """Return the cut-wet, first-dry, and RH rows used by the seam basis."""

    if not isinstance(layout, cut.CutTransportLayout) or not layout.is_square:
        raise TypeError("face-limit row basis requires a square CutTransportLayout")
    nw = layout.wet_piece_count
    nd = layout.dry_piece_count
    return (
        nw - 1,
        2 * nw - 1,
        2 * nw,
        2 * nw + nd,
        2 * nw + 2 * nd,
        2 * nw + 3 * nd,
        2 * nw + 3 * nd + 1,
        2 * nw + 3 * nd + 2,
    )


def _face_limit_regularized_residual_vector(
    original_covariant: Sequence[float],
    layout: cut.CutTransportLayout,
) -> tuple[float, ...]:
    """Apply the exact same-cell row basis with a regular inner-face limit.

    Only the first dry water, hexane, and common-datum energy rows are
    replaced.  The retained wet-cut and RH rows make this a unit-determinant
    triangular row operation, so it changes neither rank nor the exact root
    set for any strictly positive cut-piece volume.
    """

    values = tuple(float(value) for value in original_covariant)
    if len(values) != layout.residual_count:
        raise ValueError("face-limit row basis received the wrong residual rank")
    ww, we, dw, dh, de, rhw, rhh, rhe = _face_limit_residual_row_indices(layout)
    transformed = list(values)
    transformed[dw] = math.fsum((values[ww], values[dw], values[rhw]))
    transformed[dh] = math.fsum((values[dh], values[rhh]))
    transformed[de] = math.fsum((values[we], values[de], values[rhe]))
    return tuple(transformed)


def _inverse_face_limit_regularized_residual_vector(
    regularized_covariant: Sequence[float],
    layout: cut.CutTransportLayout,
) -> tuple[float, ...]:
    """Invert :func:`_face_limit_regularized_residual_vector` algebraically."""

    values = tuple(float(value) for value in regularized_covariant)
    if len(values) != layout.residual_count:
        raise ValueError("inverse face-limit row basis received the wrong residual rank")
    ww, we, dw, dh, de, rhw, rhh, rhe = _face_limit_residual_row_indices(layout)
    original = list(values)
    original[dw] = math.fsum((values[dw], -values[ww], -values[rhw]))
    original[dh] = math.fsum((values[dh], -values[rhh]))
    original[de] = math.fsum((values[de], -values[we], -values[rhe]))
    return tuple(original)


def _stable_material_changes(
    before: cut.CutTransportState,
    after: cut.CutTransportState,
    swept,
) -> MaterialInventorySnapshot:
    """Return cancellation-safe ``V_new*dC + dV*C_old`` cell changes.

    The formula and the signed volume changes are the same ones used by the
    cut residual foundation.  Accumulating these increments avoids recovering
    a tiny accepted ALE transfer by subtracting two nearly equal inventories.
    """

    if before.layout != after.layout:
        raise ValueError("stable material changes require one fixed cut chart")
    old = _piece_data(before)
    new = _piece_data(after)
    layout = before.layout
    n = before.geometry.master_grid.n
    water = [0.0] * n
    hexane = [0.0] * n
    energy = [0.0] * n

    def delta(
        new_volume: float,
        volume_change: float,
        old_concentration: float,
        new_concentration: float,
    ) -> float:
        return math.fsum(
            (
                new_volume * (new_concentration - old_concentration),
                volume_change * old_concentration,
            )
        )

    for position, cell_index in enumerate(layout.wet_cell_indices):
        old_volume = before.geometry.cells[cell_index].wet_volume_m3
        new_volume = after.geometry.cells[cell_index].wet_volume_m3
        volume_change = swept.cell_wet_volume_changes_m3[cell_index]
        if math.fsum((new_volume, -old_volume, -volume_change)) != 0.0:
            tolerance = 64.0 * math.ulp(before.geometry.master_grid.total_volume)
            if abs(math.fsum((new_volume, -old_volume, -volume_change))) > tolerance:
                raise RuntimeError("wet volume increment disagrees with shared ALE sweep")
        old_cell = old.wet_cells[position]
        new_cell = new.wet_cells[position]
        changes = (
            delta(
                new_volume,
                volume_change,
                old_cell.retained_water_concentration_mol_m3,
                new_cell.retained_water_concentration_mol_m3,
            ),
            delta(
                new_volume,
                volume_change,
                before.config.wet.wet.rho_dm_p
                * before.historical_hexane_loadings[cell_index]
                / hx.M,
                after.config.wet.wet.rho_dm_p * after.historical_hexane_loadings[cell_index] / hx.M,
            ),
            delta(
                new_volume,
                volume_change,
                old_cell.energy_density_j_m3,
                new_cell.energy_density_j_m3,
            ),
        )
        water[cell_index] = math.fsum((water[cell_index], changes[0]))
        hexane[cell_index] = math.fsum((hexane[cell_index], changes[1]))
        energy[cell_index] = math.fsum((energy[cell_index], changes[2]))

    for position, cell_index in enumerate(layout.dry_cell_indices):
        old_volume = before.geometry.cells[cell_index].dry_volume_m3
        new_volume = after.geometry.cells[cell_index].dry_volume_m3
        volume_change = swept.cell_dry_volume_changes_m3[cell_index]
        tolerance = 64.0 * math.ulp(before.geometry.master_grid.total_volume)
        if abs(math.fsum((new_volume, -old_volume, -volume_change))) > tolerance:
            raise RuntimeError("dry volume increment disagrees with shared ALE sweep")
        old_cell = old.dry_cells[position]
        new_cell = new.dry_cells[position]
        changes = (
            delta(
                new_volume,
                volume_change,
                old_cell.total_water_concentration_mol_m3,
                new_cell.total_water_concentration_mol_m3,
            ),
            delta(
                new_volume,
                volume_change,
                old_cell.total_hexane_concentration_mol_m3,
                new_cell.total_hexane_concentration_mol_m3,
            ),
            delta(
                new_volume,
                volume_change,
                old_cell.energy_density_j_m3,
                new_cell.energy_density_j_m3,
            ),
        )
        water[cell_index] = math.fsum((water[cell_index], changes[0]))
        hexane[cell_index] = math.fsum((hexane[cell_index], changes[1]))
        energy[cell_index] = math.fsum((energy[cell_index], changes[2]))
    return MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _compensated_add_vectors(
    sums: Sequence[float],
    compensations: Sequence[float],
    increments: Sequence[float],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Neumaier accumulation of signed per-cell material changes."""

    if not len(sums) == len(compensations) == len(increments):
        raise ValueError("compensated material-change vectors must align")
    new_sums: list[float] = []
    new_compensations: list[float] = []
    for current, compensation, increment in zip(sums, compensations, increments):
        total, new_compensation = _compensated_add_scalar(
            current,
            compensation,
            increment,
        )
        new_sums.append(total)
        new_compensations.append(new_compensation)
    return tuple(new_sums), tuple(new_compensations)


def _compensated_add_scalar(
    current: float,
    compensation: float,
    increment: float,
) -> tuple[float, float]:
    """One Neumaier update retaining the signed low-order correction."""

    total = current + increment
    if abs(current) >= abs(increment):
        correction = (current - total) + increment
    else:
        correction = (increment - total) + current
    return total, math.fsum((compensation, correction))


def _assembly_sweep_donor_concentration(
    before: cut.CutTransportState,
) -> float:
    """Return the wet-side ALE mass donor the assembly actually books.

    Since the O9a front-donor ruling
    (docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md, ruled
    2026-08-21) that donor is the CONSUMED piece's own bulk retained-water
    concentration - the outermost wet piece, read at the OLD time level, the
    same level ``cut_transport._piece_inventory_changes`` books the paired
    geometric storage term at - and no longer the zero-volume equilibrium
    interface trace.  This is the identical expression
    ``wet_water.evaluate_cell`` uses for
    ``retained_water_concentration_mol_m3``, in the same operation order, so it
    reproduces the assembled donor bit for bit.

    Every reconstruction of the assembled front-face water rate in this module
    reads the donor from here, so a residual SCALE can never silently drift
    away from the row it is scaling.
    """

    return before.config.wet.wet.rho_dm_p * before.wet_retained_water_loadings[-1] / wt.M


def _assembly_sweep_donor_energy_density(
    before: cut.CutTransportState,
) -> float:
    """Return the wet-side ALE energy donor the assembly actually books.

    Since the O10a energy-donor completion
    (docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md, ruled
    2026-08-21) that donor is the CONSUMED piece's OWN bulk energy density -
    the outermost wet piece, read at the same OLD time level as its mass donor
    above - and no longer the zero-volume equilibrium interface trace.  The
    PHYSICS GROUND recorded in that ruling is that two-phase moving-boundary
    (Rankine-Hugoniot) jump conditions must convect mass and enthalpy OF THE
    SAME MATERIAL STATE, an internal-consistency requirement of the balance
    laws; the material state itself - bulk, undrained - is the O9-measured fact
    (Peclet 3.4e3 forbids pre-drainage).

    This is the same old outermost wet cell ``_piece_data`` evaluates, built
    from the same primitives in the same order, so it reproduces the assembled
    donor ``old_inventory.wet_energy_density_j_m3[-1]`` bit for bit.  The
    retained-water capacity dual is deliberately left at its default: it enters
    the cell's potential and active set, never
    :attr:`wet_water.WetWaterCellState.energy_density_j_m3`.

    Every reconstruction of the assembled front-face ENERGY rate in this module
    reads the donor from here, so a residual SCALE can never silently drift
    away from the row it is scaling.
    """

    cut_cell = before.layout.wet_cell_indices[-1]
    return ww.evaluate_cell(
        before.wet_temperatures_k[-1],
        before.wet_retained_water_loadings[-1],
        before.historical_hexane_loadings[cut_cell],
        before.oil_fraction_labels[cut_cell],
        before.config.wet,
    ).energy_density_j_m3


def _residual_scales(
    before: cut.CutTransportState,
    initial: cut.CutTransportAssembly,
    dt_s: float,
) -> CutResidualScales:
    pieces = _piece_data(before)
    layout = before.layout
    grid = before.geometry.master_grid
    area_gamma = 4.0 * math.pi * initial.candidate_geometry.front.radius_m**2
    q_gamma = initial.swept_geometry.interface_swept_volume_rate_m3_s
    # O9a, Stage 1b.  ``wet_h_water`` is a VERBATIM RECONSTRUCTION of the
    # assembly's own ``wet_water_h`` face rate, so it must read the same
    # corrected donor the assembly books; otherwise this scale would bound a
    # rate no module computes.  This is the same "copied verbatim from the
    # assembly" repair Stage 1 applied to
    # ``cut_cap_active_continuation._wet_water_face_rates``.  It cannot loosen
    # a gate: every ``max`` this value enters below also carries the piece's
    # own ``inventory / dt`` term, and ``|q_gamma * C| = |dV| * C / dt`` is
    # bounded by ``V * C / dt`` for any admissible sweep.
    wet_h_water = math.fsum(
        (
            area_gamma * initial.wet_face_fluxes[-1].retained_water_flux_mol_m2_s,
            -q_gamma * _assembly_sweep_donor_concentration(before),
        )
    )
    wet_h_hexane = -q_gamma * (
        before.config.wet.wet.rho_dm_p * initial.interface.historical_hexane_loading / hx.M
    )
    dry_h_water = math.fsum(
        (
            area_gamma * initial.dry_face_fluxes[0].component.conserved_water_flux_mol_m2_s,
            -q_gamma * initial.interface.dry.total_water_concentration_mol_m3,
        )
    )
    dry_h_hexane = math.fsum(
        (
            area_gamma * initial.dry_face_fluxes[0].component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * initial.interface.dry.total_hexane_concentration_mol_m3,
        )
    )
    wet_water_rates = [0.0]
    for position, cell in enumerate(layout.wet_cell_indices[1:], start=1):
        wet_water_rates.append(
            grid.areas[cell] * initial.wet_face_fluxes[position].retained_water_flux_mol_m2_s
        )
    wet_water_rates.append(wet_h_water)
    dry_water_rates = [dry_h_water]
    dry_hexane_rates = [dry_h_hexane]
    for position in range(1, layout.dry_face_count):
        area = (
            grid.areas[-1]
            if position == layout.dry_piece_count
            else grid.areas[layout.dry_cell_indices[position]]
        )
        dry_water_rates.append(
            area * initial.dry_face_fluxes[position].component.conserved_water_flux_mol_m2_s
        )
        dry_hexane_rates.append(
            area * initial.dry_face_fluxes[position].component.conserved_hexane_flux_mol_m2_s
        )
    wet_water_scales = tuple(
        max(
            pieces.wet_water_mol[index] / dt_s,
            abs(wet_water_rates[index]),
            abs(wet_water_rates[index + 1]),
        )
        for index in range(layout.wet_piece_count)
    )
    wet_energy_scales = tuple(value / dt_s for value in pieces.wet_capacity_energy_j)
    dry_water_scales = tuple(
        max(
            pieces.dry_water_mol[index] / dt_s,
            abs(dry_water_rates[index]),
            abs(dry_water_rates[index + 1]),
        )
        for index in range(layout.dry_piece_count)
    )
    dry_hexane_scales = tuple(
        max(
            pieces.dry_hexane_mol[index] / dt_s,
            abs(dry_hexane_rates[index]),
            abs(dry_hexane_rates[index + 1]),
        )
        for index in range(layout.dry_piece_count)
    )
    dry_energy_scales = tuple(value / dt_s for value in pieces.dry_capacity_energy_j)
    cut_wet_position = layout.wet_piece_count - 1
    cut_dry_position = 0
    return CutResidualScales(
        wet_water_mol_s=wet_water_scales,
        wet_energy_w=wet_energy_scales,
        dry_water_mol_s=dry_water_scales,
        dry_hexane_mol_s=dry_hexane_scales,
        dry_energy_w=dry_energy_scales,
        rh_water_mol_s=max(
            pieces.wet_water_mol[cut_wet_position] / dt_s,
            pieces.dry_water_mol[cut_dry_position] / dt_s,
            abs(wet_h_water),
            abs(dry_h_water),
        ),
        rh_hexane_mol_s=max(
            pieces.wet_hexane_mol[cut_wet_position] / dt_s,
            pieces.dry_hexane_mol[cut_dry_position] / dt_s,
            abs(wet_h_hexane),
            abs(dry_h_hexane),
        ),
        rh_energy_w=(
            pieces.wet_capacity_energy_j[cut_wet_position]
            + pieces.dry_capacity_energy_j[cut_dry_position]
        )
        / dt_s,
        dry_energy_capacity_stencil_audits=pieces.dry_capacity_stencil_audits,
    )


def _candidate_piece_data_for_scales(
    before: cut.CutTransportState,
    initial: cut.CutTransportAssembly,
) -> _PieceData:
    """Evaluate candidate inventories/capacities without constructing a state."""

    layout = initial.layout
    geometry = initial.candidate_geometry
    config = initial.target_config
    wet_volumes = tuple(geometry.cells[cell].wet_volume_m3 for cell in layout.wet_cell_indices)
    dry_volumes = tuple(geometry.cells[cell].dry_volume_m3 for cell in layout.dry_cell_indices)
    wet_water = tuple(
        volume * cell.retained_water_concentration_mol_m3
        for volume, cell in zip(wet_volumes, initial.wet_cells)
    )
    wet_hexane = tuple(
        volume * config.wet.wet.rho_dm_p * before.historical_hexane_loadings[index] / hx.M
        for volume, index in zip(wet_volumes, layout.wet_cell_indices)
    )
    wet_energy = tuple(
        volume * cell.energy_density_j_m3 for volume, cell in zip(wet_volumes, initial.wet_cells)
    )
    dry_water = tuple(
        volume * cell.total_water_concentration_mol_m3
        for volume, cell in zip(dry_volumes, initial.dry_cells)
    )
    dry_hexane = tuple(
        volume * cell.total_hexane_concentration_mol_m3
        for volume, cell in zip(dry_volumes, initial.dry_cells)
    )
    dry_energy = tuple(
        volume * cell.energy_density_j_m3 for volume, cell in zip(dry_volumes, initial.dry_cells)
    )
    wet_capacity = tuple(
        volume
        * wet_core.energy_capacity(
            cell.temperature_k,
            cell.historical_hexane_loading,
            replace(
                config.wet.wet,
                X_water=cell.retained_water_loading,
                w_o=cell.oil_fraction_label,
            ),
        )
        * max(abs(cell.temperature_k), 1.0)
        for volume, cell in zip(wet_volumes, initial.wet_cells)
    )
    dry_capacity_results = tuple(
        _dry_energy_capacity_with_audit(
            cell.temperature_k,
            cell.y_hexane,
            replace(
                config.dry.pore,
                w_o=before.oil_fraction_labels[index],
            ),
            config,
        )
        for cell, index in zip(
            initial.dry_cells,
            layout.dry_cell_indices,
        )
    )
    dry_capacity = tuple(
        volume * result[0] * max(abs(cell.temperature_k), 1.0)
        for volume, cell, result in zip(
            dry_volumes,
            initial.dry_cells,
            dry_capacity_results,
        )
    )
    return _PieceData(
        initial.wet_cells,
        initial.dry_cells,
        wet_water,
        wet_hexane,
        wet_energy,
        dry_water,
        dry_hexane,
        dry_energy,
        wet_capacity,
        dry_capacity,
        tuple(result[1] for result in dry_capacity_results),
    )


def _regularized_extensive_rates(
    before: cut.CutTransportState,
    initial: cut.CutTransportAssembly,
) -> tuple[
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[float, ...],
]:
    """Return water, hexane, and datum-reduced energy face rates."""

    layout = initial.layout
    grid = before.geometry.master_grid
    area_gamma = 4.0 * math.pi * initial.candidate_geometry.front.radius_m**2
    q_gamma = initial.swept_geometry.interface_swept_volume_rate_m3_s
    # O9a, Stage 1b.  Same "copied verbatim from the assembly" contract as
    # ``_residual_scales``: this reconstruction of the assembled front-face
    # water rate reads the corrected donor so the face-limit regularized scales
    # bound the rows the solver actually assembles.
    #
    # O10a, same family pass: ``wet_energy_h`` below is the matching verbatim
    # reconstruction of the assembly's own front-face ENERGY rate, so it moves
    # with the mass rate onto the consumed piece's own bulk energy density.
    # The ruling's physics ground is that two-phase moving-boundary
    # (Rankine-Hugoniot) jump conditions must convect mass and enthalpy OF THE
    # SAME MATERIAL STATE - an internal-consistency requirement of the balance
    # laws - at the bulk, undrained state the O9 traverse measured (Peclet
    # 3.4e3).  Leaving this one on the trace would make the energy scale bound
    # a rate no module computes.  The hexane sweep is untouched; only the water
    # mass and the energy donors moved.
    wet_water_h = math.fsum(
        (
            area_gamma * initial.wet_face_fluxes[-1].retained_water_flux_mol_m2_s,
            -q_gamma * _assembly_sweep_donor_concentration(before),
        )
    )
    wet_hexane_h = -q_gamma * (
        before.config.wet.wet.rho_dm_p * initial.interface.historical_hexane_loading / hx.M
    )
    wet_energy_h = math.fsum(
        (
            area_gamma * initial.wet_face_fluxes[-1].total_energy_flux_w_m2,
            # O10a: the corrected energy donor, read from the same consumed
            # piece and the same time level as the water donor above.
            -q_gamma * _assembly_sweep_donor_energy_density(before),
            -hx.H_REF * wet_hexane_h,
        )
    )
    dry_water_h = math.fsum(
        (
            area_gamma * initial.dry_face_fluxes[0].component.conserved_water_flux_mol_m2_s,
            -q_gamma * initial.interface.dry.total_water_concentration_mol_m3,
        )
    )
    dry_hexane_h = math.fsum(
        (
            area_gamma * initial.dry_face_fluxes[0].component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * initial.interface.dry.total_hexane_concentration_mol_m3,
        )
    )
    dry_energy_h = math.fsum(
        (
            area_gamma * initial.dry_face_fluxes[0].energy.total_energy_flux_w_m2,
            -q_gamma * initial.interface.dry.energy_density_j_m3,
            -hx.H_REF * dry_hexane_h,
        )
    )

    wet_water_rates = [0.0]
    wet_energy_rates = [0.0]
    for position, cell in enumerate(layout.wet_cell_indices[1:], start=1):
        area = grid.areas[cell]
        wet_water_rates.append(
            area * initial.wet_face_fluxes[position].retained_water_flux_mol_m2_s
        )
        wet_energy_rates.append(area * initial.wet_face_fluxes[position].total_energy_flux_w_m2)
    wet_water_rates.append(wet_water_h)
    wet_energy_rates.append(wet_energy_h)

    dry_water_rates = [dry_water_h]
    dry_hexane_rates = [dry_hexane_h]
    dry_energy_rates = [dry_energy_h]
    for position in range(1, layout.dry_face_count):
        area = (
            grid.areas[-1]
            if position == layout.dry_piece_count
            else grid.areas[layout.dry_cell_indices[position]]
        )
        face_flux = initial.dry_face_fluxes[position]
        dry_water_rates.append(area * face_flux.component.conserved_water_flux_mol_m2_s)
        dry_hexane_rate = area * face_flux.component.conserved_hexane_flux_mol_m2_s
        dry_hexane_rates.append(dry_hexane_rate)
        dry_energy_rates.append(
            math.fsum(
                (
                    area * face_flux.energy.total_energy_flux_w_m2,
                    -hx.H_REF * dry_hexane_rate,
                )
            )
        )
    return (
        tuple(wet_water_rates),
        tuple(wet_energy_rates),
        tuple(dry_water_rates),
        tuple(dry_hexane_rates),
        tuple(dry_energy_rates),
    )


def _face_limit_regularized_residual_scales(
    before: cut.CutTransportState,
    initial: cut.CutTransportAssembly,
    dt_s: float,
) -> CutResidualScales:
    """Return fixed seed-local scales for the sparse face-limit row basis."""

    old = _piece_data(before)
    new = _candidate_piece_data_for_scales(before, initial)
    layout = before.layout
    floor = NONLINEAR_RESIDUAL_SCALE_FLOOR
    original = _datum_covariant_residual_vector(initial.residuals)
    regularized = _face_limit_regularized_residual_vector(original, layout)
    (
        wet_water_rates,
        wet_energy_rates,
        dry_water_rates,
        dry_hexane_rates,
        dry_energy_rates,
    ) = _regularized_extensive_rates(before, initial)
    nw = layout.wet_piece_count
    nd = layout.dry_piece_count
    _, _, dw, dh, de, _, _, _ = _face_limit_residual_row_indices(layout)

    wet_water_scales = tuple(
        max(
            old.wet_water_mol[index] / dt_s,
            new.wet_water_mol[index] / dt_s,
            abs(wet_water_rates[index]),
            abs(wet_water_rates[index + 1]),
            abs(initial.residuals.wet_water_mol_s[index]),
            floor,
        )
        for index in range(nw)
    )
    wet_energy_scales = tuple(
        max(
            old.wet_capacity_energy_j[index] / dt_s,
            new.wet_capacity_energy_j[index] / dt_s,
            abs(wet_energy_rates[index]),
            abs(wet_energy_rates[index + 1]),
            abs(initial.residuals.wet_energy_w[index]),
            floor,
        )
        for index in range(nw)
    )
    dry_water_scales = [
        max(
            old.dry_water_mol[index] / dt_s,
            new.dry_water_mol[index] / dt_s,
            abs(dry_water_rates[index]),
            abs(dry_water_rates[index + 1]),
            abs(initial.residuals.dry_water_mol_s[index]),
            floor,
        )
        for index in range(nd)
    ]
    dry_hexane_scales = [
        max(
            old.dry_hexane_mol[index] / dt_s,
            new.dry_hexane_mol[index] / dt_s,
            abs(dry_hexane_rates[index]),
            abs(dry_hexane_rates[index + 1]),
            abs(initial.residuals.dry_hexane_mol_s[index]),
            floor,
        )
        for index in range(nd)
    ]
    dry_energy_scales = [
        max(
            old.dry_capacity_energy_j[index] / dt_s,
            new.dry_capacity_energy_j[index] / dt_s,
            abs(dry_energy_rates[index]),
            abs(dry_energy_rates[index + 1]),
            abs(original[2 * nw + 2 * nd + index]),
            floor,
        )
        for index in range(nd)
    ]

    cut_wet = nw - 1
    dry_water_scales[0] = max(
        math.fsum((old.wet_water_mol[cut_wet], old.dry_water_mol[0])) / dt_s,
        math.fsum((new.wet_water_mol[cut_wet], new.dry_water_mol[0])) / dt_s,
        abs(wet_water_rates[cut_wet]),
        abs(dry_water_rates[1]),
        abs(regularized[dw]),
        floor,
    )
    dry_hexane_scales[0] = max(
        math.fsum((old.wet_hexane_mol[cut_wet], old.dry_hexane_mol[0])) / dt_s,
        math.fsum((new.wet_hexane_mol[cut_wet], new.dry_hexane_mol[0])) / dt_s,
        abs(dry_hexane_rates[1]),
        abs(regularized[dh]),
        floor,
    )
    dry_energy_scales[0] = max(
        math.fsum(
            (
                old.wet_capacity_energy_j[cut_wet],
                old.dry_capacity_energy_j[0],
            )
        )
        / dt_s,
        math.fsum(
            (
                new.wet_capacity_energy_j[cut_wet],
                new.dry_capacity_energy_j[0],
            )
        )
        / dt_s,
        abs(wet_energy_rates[cut_wet]),
        abs(dry_energy_rates[1]),
        abs(regularized[de]),
        floor,
    )
    return CutResidualScales(
        wet_water_mol_s=wet_water_scales,
        wet_energy_w=wet_energy_scales,
        dry_water_mol_s=tuple(dry_water_scales),
        dry_hexane_mol_s=tuple(dry_hexane_scales),
        dry_energy_w=tuple(dry_energy_scales),
        rh_water_mol_s=max(
            old.wet_water_mol[cut_wet] / dt_s,
            old.dry_water_mol[0] / dt_s,
            new.wet_water_mol[cut_wet] / dt_s,
            new.dry_water_mol[0] / dt_s,
            abs(wet_water_rates[-1]),
            abs(dry_water_rates[0]),
            abs(initial.residuals.rh_water_mol_s),
            floor,
        ),
        rh_hexane_mol_s=max(
            old.wet_hexane_mol[cut_wet] / dt_s,
            old.dry_hexane_mol[0] / dt_s,
            new.wet_hexane_mol[cut_wet] / dt_s,
            new.dry_hexane_mol[0] / dt_s,
            abs(dry_hexane_rates[0]),
            abs(initial.residuals.rh_hexane_mol_s),
            floor,
        ),
        rh_energy_w=max(
            math.fsum(
                (
                    old.wet_capacity_energy_j[cut_wet],
                    old.dry_capacity_energy_j[0],
                )
            )
            / dt_s,
            math.fsum(
                (
                    new.wet_capacity_energy_j[cut_wet],
                    new.dry_capacity_energy_j[0],
                )
            )
            / dt_s,
            abs(wet_energy_rates[-1]),
            abs(dry_energy_rates[0]),
            abs(original[-1]),
            floor,
        ),
        dry_energy_capacity_stencil_audits=(
            *old.dry_capacity_stencil_audits,
            *new.dry_capacity_stencil_audits,
        ),
    )


def _coordinate_solver_bounds(
    layout: cut.CutTransportLayout,
    maximum_logit: float,
    *,
    wet_retained_cap_enabled: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(wet_retained_cap_enabled, bool):
        raise TypeError("wet retained-cap coordinate policy must be boolean")
    wet_count = layout.wet_piece_count
    dry_bounded_count = 2 * layout.dry_piece_count
    nt_count = layout.dry_face_count
    bounded_after_nt = 2
    lower = np.asarray(
        (-maximum_logit,) * wet_count
        + (-maximum_logit,) * wet_count
        + (-maximum_logit,) * dry_bounded_count
        + (-math.inf,) * nt_count
        + (-maximum_logit,) * bounded_after_nt,
        dtype=float,
    )
    upper = np.asarray(
        (maximum_logit,) * wet_count
        + ((math.inf,) * wet_count if wet_retained_cap_enabled else (maximum_logit,) * wet_count)
        + (maximum_logit,) * dry_bounded_count
        + (math.inf,) * nt_count
        + (maximum_logit,) * bounded_after_nt,
        dtype=float,
    )
    return lower, upper


def _coordinate_solver_bounds_for_chart(
    layout: cut.CutTransportLayout,
    maximum_logit: float,
    chart: _CoordinateChart,
) -> tuple[np.ndarray, np.ndarray]:
    """Return unchanged bounds except for an opted-in physical ``front_z``."""

    lower, upper = _coordinate_solver_bounds(
        layout,
        maximum_logit,
        wet_retained_cap_enabled=(chart.wet_retained_cap is not None),
    )
    if chart.front_coordinate_policy is SameCellFrontCoordinatePolicy.WHOLE_CELL_OPEN_LOGIT:
        return lower, upper
    if chart.front_coordinate_policy is not SameCellFrontCoordinatePolicy.DIRECT_PHYSICAL_FRONT_Z:
        raise ValueError("unknown same-cell front-coordinate policy")
    first_interior_front = math.nextafter(chart.front_z[0], chart.front_z[1])
    last_interior_front = math.nextafter(chart.front_z[1], chart.front_z[0])
    if not first_interior_front < last_interior_front:
        raise NoUsableDirectPhysicalFrontInteriorError(
            "direct physical front chart has no strictly ordered binary64 interior guards"
        )
    lower = lower.copy()
    upper = upper.copy()
    lower[-2] = first_interior_front
    upper[-2] = last_interior_front
    return lower, upper


def _structural_block(
    blocks: Sequence[cut_sparsity.StructuralBlock],
    name: str,
) -> cut_sparsity.StructuralBlock:
    matches = tuple(block for block in blocks if block.name == name)
    if len(matches) != 1:
        raise RuntimeError(f"scalar post-root certificate needs one {name!r} block")
    return matches[0]


def _structural_location(
    blocks: Sequence[cut_sparsity.StructuralBlock],
    global_index: int,
) -> tuple[cut_sparsity.StructuralBlock, int]:
    for block in blocks:
        if block.start <= global_index < block.stop:
            return block, global_index - block.start
    raise RuntimeError("scalar post-root index is outside the declared square rank")


def _dominant_scaled_residual_location(
    sparsity: cut_sparsity.CutJacobianSparsity,
    signed_scaled_residual: Sequence[float],
) -> tuple[cut_sparsity.StructuralBlock, int, int, float]:
    """Locate the largest absolute row, breaking exact ties by global index."""

    residual = np.asarray(signed_scaled_residual, dtype=float)
    if residual.shape != (sparsity.layout.residual_count,) or not np.all(np.isfinite(residual)):
        raise RuntimeError("scalar post-root residual has the wrong finite square rank")
    absolute = np.abs(residual)
    global_index = int(np.argmax(absolute))
    block, local_index = _structural_location(
        sparsity.residual_blocks,
        global_index,
    )
    return block, local_index, global_index, float(absolute[global_index])


def _validate_scalar_post_root_point(
    before: CutIntegratorState,
    chart: _CoordinateChart,
    sparsity: cut_sparsity.CutJacobianSparsity,
    coordinates: Sequence[float],
    lower_coordinates: np.ndarray,
    upper_coordinates: np.ndarray,
    assembly: cut.CutTransportAssembly | None = None,
) -> np.ndarray:
    """Require one strict same-rank, same-topology, smooth-cap trial point."""

    point = np.asarray(coordinates, dtype=float)
    expected_shape = (sparsity.layout.unknown_count,)
    if point.shape != expected_shape or not np.all(np.isfinite(point)):
        raise RuntimeError("scalar post-root point has the wrong finite square rank")
    if not (
        lower_coordinates.shape == expected_shape
        and upper_coordinates.shape == expected_shape
        and np.all(point > lower_coordinates)
        and np.all(point < upper_coordinates)
    ):
        raise RuntimeError("scalar post-root point left the strict solver chart")

    if chart.wet_retained_cap is not None:
        wet_water = _structural_block(
            sparsity.unknown_blocks,
            "wet_retained_water",
        )
        if np.any(point[wet_water.slice] >= 0.0):
            raise RuntimeError(
                "scalar post-root safeguard excludes retained-cap contact and active branches"
            )

    before_cap = wrc.certify_exact_graph(
        before.transport.wet_retained_water_loadings,
        before.transport.effective_wet_retained_water_capacity_duals_over_rt,
        before.transport.config.wet.luikov,
    )
    if before_cap.active_piece_count != 0:
        raise RuntimeError(
            "scalar post-root safeguard requires a wholly smooth retained-cap source state"
        )

    if assembly is not None:
        if assembly.layout != sparsity.layout or assembly.layout != before.transport.layout:
            raise RuntimeError("scalar post-root evaluation changed topology or equation rank")
        if len(assembly.residuals.vector) != sparsity.layout.residual_count:
            raise RuntimeError("scalar post-root evaluation changed residual rank")
        cap = assembly.wet_retained_cap_certificate
        if not cap.exact_graph_without_tolerance or cap.active_piece_count != 0:
            raise RuntimeError(
                "scalar post-root evaluation left the exact smooth retained-cap graph"
            )
        _bounded_fractions(assembly.candidate, chart)
    return point


def _central_coordinate_points(
    coordinates: np.ndarray,
    column: int,
) -> tuple[np.ndarray, np.ndarray]:
    base = float(coordinates[column])
    half_width = SCALAR_POST_ROOT_CENTRAL_RELATIVE_STEP * max(1.0, abs(base))
    lower_value = base - half_width
    upper_value = base + half_width
    if not (
        math.isfinite(lower_value)
        and math.isfinite(upper_value)
        and lower_value < base < upper_value
    ):
        raise RuntimeError("scalar post-root central stencil is not representable")
    lower = coordinates.copy()
    upper = coordinates.copy()
    lower[column] = lower_value
    upper[column] = upper_value
    return lower, upper


def _select_strongest_paired_dry_composition_coordinate(
    before: CutIntegratorState,
    chart: _CoordinateChart,
    sparsity: cut_sparsity.CutJacobianSparsity,
    coordinates: np.ndarray,
    residual_row: int,
    lower_coordinates: np.ndarray,
    upper_coordinates: np.ndarray,
    evaluate_scaled: Callable[
        [np.ndarray],
        tuple[np.ndarray, cut.CutTransportAssembly],
    ],
) -> tuple[cut_sparsity.StructuralBlock, int, int, float, int]:
    """Select the strongest structurally paired dry-composition coordinate."""

    dry_y = _structural_block(sparsity.unknown_blocks, "dry_y_hexane")
    supported = tuple(
        column
        for column in sparsity.row_columns[residual_row]
        if dry_y.start <= column < dry_y.stop
    )
    if not supported:
        raise RuntimeError(
            "dominant dry-water row has no physically paired dry-composition coordinate"
        )

    best_column: int | None = None
    best_derivative = 0.0
    evaluations = 0
    for column in supported:
        lower_point, upper_point = _central_coordinate_points(coordinates, column)
        lower_vector, lower_assembly = evaluate_scaled(lower_point)
        upper_vector, upper_assembly = evaluate_scaled(upper_point)
        evaluations += 2
        _validate_scalar_post_root_point(
            before,
            chart,
            sparsity,
            lower_point,
            lower_coordinates,
            upper_coordinates,
            lower_assembly,
        )
        _validate_scalar_post_root_point(
            before,
            chart,
            sparsity,
            upper_point,
            lower_coordinates,
            upper_coordinates,
            upper_assembly,
        )
        denominator = float(upper_point[column] - lower_point[column])
        derivative = float((upper_vector[residual_row] - lower_vector[residual_row]) / denominator)
        if not math.isfinite(derivative):
            raise RuntimeError("scalar post-root coordinate derivative is not finite")
        if best_column is None or abs(derivative) > abs(best_derivative):
            best_column = column
            best_derivative = derivative

    if best_column is None or best_derivative == 0.0:
        raise RuntimeError("scalar post-root dry-composition sensitivity is exactly zero")
    return (
        dry_y,
        best_column - dry_y.start,
        best_column,
        best_derivative,
        evaluations,
    )


def _central_condition_certificate(
    before: CutIntegratorState,
    chart: _CoordinateChart,
    sparsity: cut_sparsity.CutJacobianSparsity,
    coordinates: np.ndarray,
    lower_coordinates: np.ndarray,
    upper_coordinates: np.ndarray,
    evaluate_scaled: Callable[
        [np.ndarray],
        tuple[np.ndarray, cut.CutTransportAssembly],
    ],
) -> tuple[float, int, int]:
    """Independently certify full rank and conditioning at the polished root."""

    dimension = sparsity.layout.unknown_count
    jacobian = np.empty((dimension, dimension), dtype=float)
    evaluations = 0
    for column in range(dimension):
        lower_point, upper_point = _central_coordinate_points(coordinates, column)
        lower_vector, lower_assembly = evaluate_scaled(lower_point)
        upper_vector, upper_assembly = evaluate_scaled(upper_point)
        evaluations += 2
        _validate_scalar_post_root_point(
            before,
            chart,
            sparsity,
            lower_point,
            lower_coordinates,
            upper_coordinates,
            lower_assembly,
        )
        _validate_scalar_post_root_point(
            before,
            chart,
            sparsity,
            upper_point,
            lower_coordinates,
            upper_coordinates,
            upper_assembly,
        )
        denominator = float(upper_point[column] - lower_point[column])
        jacobian[:, column] = (upper_vector - lower_vector) / denominator
    if not np.all(np.isfinite(jacobian)):
        raise RuntimeError("scalar post-root final central Jacobian is not finite")
    rank = int(np.linalg.matrix_rank(jacobian))
    condition = _jacobian_condition_proxy(jacobian)
    if rank != dimension:
        raise RuntimeError(f"scalar post-root final Jacobian lost rank: {rank} < {dimension}")
    if not math.isfinite(condition) or condition > before.controls.maximum_condition_proxy:
        raise RuntimeError(
            f"scalar post-root final Jacobian condition proxy is uncertified: {condition:.3e}"
        )
    return condition, rank, evaluations


def _attempt_scalar_post_root_polish(
    before: CutIntegratorState,
    chart: _CoordinateChart,
    sparsity: cut_sparsity.CutJacobianSparsity,
    coordinates: Sequence[float],
    signed_scaled_residual: Sequence[float],
    lower_coordinates: np.ndarray,
    upper_coordinates: np.ndarray,
    evaluate_scaled: Callable[
        [np.ndarray],
        tuple[np.ndarray, cut.CutTransportAssembly],
    ],
) -> tuple[
    np.ndarray,
    cut.CutTransportAssembly,
    tuple[float, ...],
    float,
    ScalarPostRootPolishAudit,
]:
    """Apply the frozen deterministic scalar safeguard or fail without a state."""

    tolerance = before.controls.nonlinear_residual_tolerance
    ceiling = scalar_post_root_near_miss_ceiling(tolerance)
    point = _validate_scalar_post_root_point(
        before,
        chart,
        sparsity,
        coordinates,
        lower_coordinates,
        upper_coordinates,
    ).copy()
    residual = np.asarray(signed_scaled_residual, dtype=float)
    residual_block, residual_local, residual_global, maximum_before = (
        _dominant_scaled_residual_location(sparsity, residual)
    )
    if not tolerance < maximum_before <= ceiling:
        raise RuntimeError(
            "scalar post-root safeguard is restricted to the frozen near-miss band "
            "bounded by the GT-PS-2-P1E-21 pinned absolute ceiling"
        )
    if residual_block.name != "dry_water":
        raise RuntimeError(
            "scalar post-root safeguard is restricted to a dominant dry-water balance"
        )
    interface_tolerance = before.transport.config.interface_composition.log_fugacity_tolerance
    if interface_tolerance != SCALAR_POST_ROOT_REQUIRED_INTERFACE_TOLERANCE:
        raise RuntimeError(
            "scalar post-root safeguard requires the unchanged generic interface tolerance"
        )

    (
        coordinate_block,
        coordinate_local,
        coordinate_global,
        derivative,
        selection_evaluations,
    ) = _select_strongest_paired_dry_composition_coordinate(
        before,
        chart,
        sparsity,
        point,
        residual_global,
        lower_coordinates,
        upper_coordinates,
        evaluate_scaled,
    )
    predicted_correction = -float(residual[residual_global]) / derivative
    if not math.isfinite(predicted_correction) or predicted_correction == 0.0:
        raise RuntimeError("scalar post-root predicted correction is not finite and nonzero")
    direction = math.inf if predicted_correction > 0.0 else -math.inf
    adjacent = math.nextafter(float(point[coordinate_global]), direction)
    ulp = abs(adjacent - float(point[coordinate_global]))
    ulp_ratio = abs(predicted_correction) / ulp
    if not math.isfinite(ulp) or ulp <= 0.0 or not math.isfinite(ulp_ratio):
        raise RuntimeError("scalar post-root predicted correction has no finite ULP count")
    predicted_ulps = max(1, math.ceil(ulp_ratio))

    base_coordinate = float(point[coordinate_global])
    base_row = float(residual[residual_global])
    bracket_probe_count = 0
    probe_coordinate: float | None = None
    probe_row: float | None = None
    for expansion in range(SCALAR_POST_ROOT_MAXIMUM_BRACKET_EXPANSIONS + 1):
        bracket_probe_count += 1
        count = predicted_ulps * (2**expansion)
        trial_coordinate = base_coordinate + math.copysign(float(count) * ulp, predicted_correction)
        if trial_coordinate == base_coordinate:
            trial_coordinate = adjacent
        if not math.isfinite(trial_coordinate):
            raise RuntimeError("scalar post-root bracket probe is not finite")
        trial = point.copy()
        trial[coordinate_global] = trial_coordinate
        trial_vector, trial_assembly = evaluate_scaled(trial)
        _validate_scalar_post_root_point(
            before,
            chart,
            sparsity,
            trial,
            lower_coordinates,
            upper_coordinates,
            trial_assembly,
        )
        trial_row = float(trial_vector[residual_global])
        if not math.isfinite(trial_row):
            raise RuntimeError("scalar post-root bracket row is not finite")
        if base_row * trial_row <= 0.0:
            probe_coordinate = trial_coordinate
            probe_row = trial_row
            break
    if probe_coordinate is None or probe_row is None:
        raise RuntimeError("scalar post-root safeguard found no representable row sign bracket")

    if base_coordinate < probe_coordinate:
        bracket_lower = base_coordinate
        bracket_upper = probe_coordinate
        row_lower = base_row
        row_upper = probe_row
    else:
        bracket_lower = probe_coordinate
        bracket_upper = base_coordinate
        row_lower = probe_row
        row_upper = base_row

    def scalar_row(value: float) -> float:
        trial = point.copy()
        trial[coordinate_global] = value
        trial_vector, trial_assembly = evaluate_scaled(trial)
        _validate_scalar_post_root_point(
            before,
            chart,
            sparsity,
            trial,
            lower_coordinates,
            upper_coordinates,
            trial_assembly,
        )
        return float(trial_vector[residual_global])

    scalar_root, scalar_result = optimize.brentq(
        scalar_row,
        bracket_lower,
        bracket_upper,
        xtol=np.nextafter(0.0, 1.0),
        rtol=4.0 * np.finfo(float).eps,
        maxiter=SCALAR_POST_ROOT_MAXIMUM_ITERATIONS,
        full_output=True,
        disp=False,
    )
    if not scalar_result.converged:
        raise RuntimeError("scalar post-root safeguarded solve did not converge")
    polished = point.copy()
    polished[coordinate_global] = float(scalar_root)
    polished_vector, polished_assembly = evaluate_scaled(polished)
    _validate_scalar_post_root_point(
        before,
        chart,
        sparsity,
        polished,
        lower_coordinates,
        upper_coordinates,
        polished_assembly,
    )
    maximum_after = float(np.max(np.abs(polished_vector)))
    if maximum_after > tolerance:
        raise RuntimeError(
            "scalar post-root full residual contract failed after the row solve: "
            f"{maximum_after:.3e} > {tolerance:.3e}"
        )

    condition, rank, condition_evaluations = _central_condition_certificate(
        before,
        chart,
        sparsity,
        polished,
        lower_coordinates,
        upper_coordinates,
        evaluate_scaled,
    )
    audit = ScalarPostRootPolishAudit(
        near_miss_multiplier_ceiling=SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER,
        maximum_scaled_residual_before=maximum_before,
        dominant_residual_block=residual_block.name,
        dominant_residual_local_index=residual_local,
        dominant_residual_global_index=residual_global,
        selected_coordinate_block=coordinate_block.name,
        selected_coordinate_local_index=coordinate_local,
        selected_coordinate_global_index=coordinate_global,
        selected_coordinate_derivative=derivative,
        selection_stencil_evaluations=selection_evaluations,
        predicted_scalar_correction=predicted_correction,
        predicted_scalar_correction_ulps=float(predicted_ulps),
        bracket_probe_count=bracket_probe_count,
        bracket_lower_coordinate=bracket_lower,
        bracket_upper_coordinate=bracket_upper,
        bracket_lower_row_residual=row_lower,
        bracket_upper_row_residual=row_upper,
        scalar_iterations=int(scalar_result.iterations),
        scalar_function_calls=int(scalar_result.function_calls),
        scalar_root_coordinate=float(scalar_root),
        scalar_root_delta=float(scalar_root - base_coordinate),
        maximum_scaled_residual_after=maximum_after,
        final_condition_proxy=condition,
        final_jacobian_rank=rank,
        final_jacobian_dimension=sparsity.layout.unknown_count,
        final_condition_stencil_evaluations=condition_evaluations,
        used_generic_interface_tolerance_unchanged=True,
    )
    return (
        polished,
        polished_assembly,
        tuple(float(value) for value in polished_vector),
        condition,
        audit,
    )


def _jacobian_condition_proxy(jacobian: np.ndarray) -> float:
    matrix = (
        np.asarray(jacobian.toarray(), dtype=float)
        if hasattr(jacobian, "toarray")
        else np.asarray(jacobian, dtype=float)
    )
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return math.inf
    condition = float(np.linalg.cond(matrix))
    return condition


def _accept_step(
    before: CutIntegratorState,
    assembly: cut.CutTransportAssembly,
    boundary: ct.PoreBoundary,
    seed: CutStepSeed,
    scales: CutResidualScales,
    optimization_scales: CutResidualScales,
    scaled_residuals: tuple[float, ...],
    dt_s: float,
    evaluations: int,
    rejected_trials: int,
    nonlinear_message: str,
    condition_proxy: float,
    chart: _CoordinateChart,
    sparsity,
    jacobian_route_reason: str,
    nonlinear_method: str,
    scalar_post_root_polish_audit: ScalarPostRootPolishAudit | None = None,
) -> CutIntegratorStep:
    candidate = assembly.candidate
    after_transport = cut.CutTransportState(
        geometry=assembly.candidate_geometry,
        config=assembly.target_config,
        time_s=before.transport.time_s + dt_s,
        wet_temperatures_k=candidate.wet_temperatures_k,
        wet_retained_water_loadings=candidate.wet_retained_water_loadings,
        dry_temperatures_k=candidate.dry_temperatures_k,
        dry_y_hexane=candidate.dry_y_hexane,
        historical_hexane_loadings=before.transport.historical_hexane_loadings,
        oil_fraction_labels=before.transport.oil_fraction_labels,
        wet_retained_water_capacity_duals_over_rt=(
            candidate.effective_wet_retained_water_capacity_duals_over_rt
        ),
    )
    old_inventory = inventory_snapshot(before.transport)
    new_inventory = inventory_snapshot(after_transport)
    area = after_transport.geometry.master_grid.areas[-1]
    surface = assembly.dry_face_fluxes[-1]
    water_out = dt_s * area * surface.component.conserved_water_flux_mol_m2_s
    hexane_out = dt_s * area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_out = dt_s * area * surface.energy.total_energy_flux_w_m2
    stable_changes = _stable_material_changes(
        before.transport,
        after_transport,
        assembly.swept_geometry,
    )
    water_step_residual = math.fsum((*stable_changes.water_cell_mol, water_out))
    hexane_step_residual = math.fsum((*stable_changes.hexane_cell_mol, hexane_out))
    energy_step_residual = math.fsum((*stable_changes.energy_cell_j, energy_out))
    step_water_scale = max(
        old_inventory.total_water_mol,
        new_inventory.total_water_mol,
        abs(water_out),
    )
    step_hexane_scale = max(
        old_inventory.total_hexane_mol,
        new_inventory.total_hexane_mol,
        abs(hexane_out),
    )
    step_energy_scale = max(
        capacity_energy_scale(before.transport),
        capacity_energy_scale(after_transport),
    )
    cumulative_water_out_high, boundary_water_compensation = _compensated_add_scalar(
        before.cumulative_boundary_water_out_mol,
        before.boundary_water_compensation_mol,
        water_out,
    )
    cumulative_hexane_out_high, boundary_hexane_compensation = _compensated_add_scalar(
        before.cumulative_boundary_hexane_out_mol,
        before.boundary_hexane_compensation_mol,
        hexane_out,
    )
    cumulative_energy_out_high, boundary_energy_compensation = _compensated_add_scalar(
        before.cumulative_boundary_energy_out_j,
        before.boundary_energy_compensation_j,
        energy_out,
    )
    cumulative_water_out = math.fsum((cumulative_water_out_high, boundary_water_compensation))
    cumulative_hexane_out = math.fsum((cumulative_hexane_out_high, boundary_hexane_compensation))
    cumulative_energy_out = math.fsum((cumulative_energy_out_high, boundary_energy_compensation))
    cumulative_water_changes, water_compensation = _compensated_add_vectors(
        before.cumulative_material_water_change_cell_mol,
        before.water_change_compensation_cell_mol,
        stable_changes.water_cell_mol,
    )
    cumulative_hexane_changes, hexane_compensation = _compensated_add_vectors(
        before.cumulative_material_hexane_change_cell_mol,
        before.hexane_change_compensation_cell_mol,
        stable_changes.hexane_cell_mol,
    )
    cumulative_energy_changes, energy_compensation = _compensated_add_vectors(
        before.cumulative_material_energy_change_cell_j,
        before.energy_change_compensation_cell_j,
        stable_changes.energy_cell_j,
    )
    corrected_cumulative_water_changes = tuple(
        math.fsum((value, compensation))
        for value, compensation in zip(
            cumulative_water_changes,
            water_compensation,
        )
    )
    corrected_cumulative_hexane_changes = tuple(
        math.fsum((value, compensation))
        for value, compensation in zip(
            cumulative_hexane_changes,
            hexane_compensation,
        )
    )
    corrected_cumulative_energy_changes = tuple(
        math.fsum((value, compensation))
        for value, compensation in zip(
            cumulative_energy_changes,
            energy_compensation,
        )
    )
    cumulative_water_residual = math.fsum(
        (*corrected_cumulative_water_changes, cumulative_water_out)
    )
    cumulative_hexane_residual = math.fsum(
        (*corrected_cumulative_hexane_changes, cumulative_hexane_out)
    )
    cumulative_energy_residual = math.fsum(
        (*corrected_cumulative_energy_changes, cumulative_energy_out)
    )
    cumulative_water_scale = max(
        before.reference_inventory.total_water_mol,
        new_inventory.total_water_mol,
        abs(cumulative_water_out),
        before.cumulative_absolute_water_transfer_mol + abs(water_out),
    )
    cumulative_hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        new_inventory.total_hexane_mol,
        abs(cumulative_hexane_out),
        before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out),
    )
    cumulative_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        step_energy_scale,
        before.cumulative_absolute_energy_transfer_j + abs(energy_out),
    )
    after = CutIntegratorState(
        transport=after_transport,
        controls=before.controls,
        reference_inventory=before.reference_inventory,
        reference_capacity_energy_scale_j=before.reference_capacity_energy_scale_j,
        last_interface_temperature_k=candidate.interface_temperature_k,
        last_total_stefan_fluxes_mol_m2_s=(candidate.dry_total_stefan_fluxes_mol_m2_s),
        cumulative_boundary_water_out_mol=cumulative_water_out_high,
        cumulative_boundary_hexane_out_mol=cumulative_hexane_out_high,
        cumulative_boundary_energy_out_j=cumulative_energy_out_high,
        boundary_water_compensation_mol=boundary_water_compensation,
        boundary_hexane_compensation_mol=boundary_hexane_compensation,
        boundary_energy_compensation_j=boundary_energy_compensation,
        cumulative_absolute_water_transfer_mol=(
            before.cumulative_absolute_water_transfer_mol + abs(water_out)
        ),
        cumulative_absolute_hexane_transfer_mol=(
            before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out)
        ),
        cumulative_absolute_energy_transfer_j=(
            before.cumulative_absolute_energy_transfer_j + abs(energy_out)
        ),
        cumulative_material_water_change_cell_mol=cumulative_water_changes,
        cumulative_material_hexane_change_cell_mol=cumulative_hexane_changes,
        cumulative_material_energy_change_cell_j=cumulative_energy_changes,
        water_change_compensation_cell_mol=water_compensation,
        hexane_change_compensation_cell_mol=hexane_compensation,
        energy_change_compensation_cell_j=energy_compensation,
        accepted_steps=before.accepted_steps + 1,
        cumulative_nonlinear_evaluations=(before.cumulative_nonlinear_evaluations + evaluations),
    )
    fractions = _bounded_fractions(candidate, chart)
    front_distances = (
        candidate.front_z - chart.front_z[0],
        chart.front_z[1] - candidate.front_z,
    )
    wet_retained_cap_before = wrc.certify_exact_graph(
        before.transport.wet_retained_water_loadings,
        before.transport.effective_wet_retained_water_capacity_duals_over_rt,
        before.transport.config.wet.luikov,
    )
    wet_retained_cap_after = assembly.wet_retained_cap_certificate
    wet_retained_cap_transition = wrc.certify_exact_transition(
        wet_retained_cap_before,
        wet_retained_cap_after,
    )
    ledger = CutIntegratorLedger(
        dt_s=dt_s,
        pressure_before_pa=before.transport.config.dry.pressure_pa,
        pressure_after_pa=assembly.target_config.dry.pressure_pa,
        finite_pressure_transition=(
            assembly.target_config.dry.pressure_pa != before.transport.config.dry.pressure_pa
        ),
        nonlinear_evaluations=evaluations,
        rejected_trial_evaluations=rejected_trials,
        nonlinear_message=nonlinear_message,
        residual_scales=scales,
        optimization_residual_scales=optimization_scales,
        scaled_residuals=scaled_residuals,
        maximum_scaled_residual=max(scaled_residuals),
        condition_proxy=condition_proxy,
        minimum_fractional_distance_to_bound=min(min(value, 1.0 - value) for value in fractions),
        minimum_front_z_distance_to_chart_boundary=min(front_distances),
        used_sparse_jacobian=sparsity is not None,
        jacobian_route_reason=jacobian_route_reason,
        front_coordinate_policy=chart.front_coordinate_policy.value,
        nonlinear_algorithm=(
            (
                "bounded_sparse_dogbox_face_limit_basis_plus_original_scalar_post_root"
                if nonlinear_method == "dogbox"
                else "bounded_sparse_trf_lsmr_face_limit_basis_plus_original_scalar_post_root"
            )
            if scalar_post_root_polish_audit is not None
            else (
                (
                    "bounded_sparse_dogbox_face_limit_basis_colored_finite_difference"
                    if nonlinear_method == "dogbox"
                    else "bounded_sparse_trf_lsmr_face_limit_basis_colored_finite_difference"
                )
                if sparsity is not None
                else "bounded_dense_trf_finite_difference"
            )
        ),
        nonlinear_residual_basis=(
            "same_cell_face_limit_regularized_with_original_certification"
            if sparsity is not None
            else "original_datum_covariant"
        ),
        scalability_class=(
            "rare near-miss scalar safeguard plus a full final central-condition "
            "certificate; runtime qualification is disclosed separately"
            if scalar_post_root_polish_audit is not None
            else (
                "bordered-banded grouped finite-difference sparse Jacobian; runtime "
                "qualification remains open"
                if sparsity is not None
                else "certified small-N dense finite-difference Jacobian oracle"
            )
        ),
        # This is the deterministic coloring exposed by our declared structural
        # pattern.  SciPy may regroup the columns internally, so nonlinear_evaluations
        # remains the authoritative work counter.
        jacobian_structural_color_count=(
            before.transport.layout.unknown_count if sparsity is None else sparsity.color_count
        ),
        jacobian_structural_nnz=(
            before.transport.layout.unknown_count**2 if sparsity is None else sparsity.nnz
        ),
        water_step_residual_mol=water_step_residual,
        hexane_step_residual_mol=hexane_step_residual,
        energy_step_residual_j=energy_step_residual,
        normalized_water_step_residual=abs(water_step_residual) / step_water_scale,
        normalized_hexane_step_residual=(abs(hexane_step_residual) / step_hexane_scale),
        normalized_energy_step_residual=abs(energy_step_residual) / step_energy_scale,
        water_cumulative_residual_mol=cumulative_water_residual,
        hexane_cumulative_residual_mol=cumulative_hexane_residual,
        energy_cumulative_residual_j=cumulative_energy_residual,
        normalized_water_cumulative_residual=(
            abs(cumulative_water_residual) / cumulative_water_scale
        ),
        normalized_hexane_cumulative_residual=(
            abs(cumulative_hexane_residual) / cumulative_hexane_scale
        ),
        normalized_energy_cumulative_residual=(
            abs(cumulative_energy_residual) / cumulative_energy_scale
        ),
        boundary_water_out_mol=water_out,
        boundary_hexane_out_mol=hexane_out,
        boundary_energy_out_j=energy_out,
        surface_film_audit=assembly.surface_film_audit,
        vanishing_wet_cut_treatment_audit=(assembly.vanishing_wet_cut_treatment_audit),
        wet_retained_cap_active_piece_indices=(wet_retained_cap_after.active_piece_indices),
        wet_retained_cap_active_piece_count=(wet_retained_cap_after.active_piece_count),
        wet_retained_cap_maximum_loading=max(
            wet_retained_cap_before.maximum_loading,
            wet_retained_cap_after.maximum_loading,
        ),
        wet_retained_cap_maximum_capacity_dual_over_rt=max(
            wet_retained_cap_before.maximum_capacity_dual_over_rt,
            wet_retained_cap_after.maximum_capacity_dual_over_rt,
        ),
        wet_retained_cap_maximum_complementarity_product=max(
            wet_retained_cap_before.maximum_complementarity_product,
            wet_retained_cap_after.maximum_complementarity_product,
        ),
        wet_retained_cap_exact_graph_without_tolerance=(
            wet_retained_cap_before.exact_graph_without_tolerance
            and wet_retained_cap_after.exact_graph_without_tolerance
            and wet_retained_cap_transition.exact_branch_comparisons_without_tolerance
        ),
        wet_retained_cap_piece_transitions=(wet_retained_cap_transition.piece_transitions),
        wet_retained_cap_entry_piece_indices=(wet_retained_cap_transition.entry_piece_indices),
        wet_retained_cap_entry_piece_count=(wet_retained_cap_transition.entry_piece_count),
        wet_retained_cap_continuation_piece_indices=(
            wet_retained_cap_transition.continuation_piece_indices
        ),
        wet_retained_cap_continuation_piece_count=(
            wet_retained_cap_transition.continuation_piece_count
        ),
        wet_retained_cap_exit_piece_indices=(wet_retained_cap_transition.exit_piece_indices),
        wet_retained_cap_exit_piece_count=(wet_retained_cap_transition.exit_piece_count),
        accepted=True,
        scalar_post_root_polish_audit=scalar_post_root_polish_audit,
    )
    return CutIntegratorStep(before, after, assembly, boundary, seed, ledger)


def _bounded_fractions(
    candidate: cut.CutTransportUnknowns, chart: _CoordinateChart
) -> tuple[float, ...]:
    pairs = [
        *((value, chart.wet_temperature) for value in candidate.wet_temperatures_k),
        *((value, chart.dry_temperature) for value in candidate.dry_temperatures_k),
        (candidate.front_z, chart.front_z),
        (candidate.interface_temperature_k, chart.interface_temperature),
    ]
    fractions = [(value - bounds[0]) / (bounds[1] - bounds[0]) for value, bounds in pairs]
    if chart.wet_retained_cap is None:
        fractions.extend(
            (value - chart.wet_water[0]) / (chart.wet_water[1] - chart.wet_water[0])
            for value in candidate.wet_retained_water_loadings
        )
    else:
        certificate = wrc.certify_exact_graph(
            candidate.wet_retained_water_loadings,
            candidate.effective_wet_retained_water_capacity_duals_over_rt,
            chart.wet_retained_cap.luikov,
        )
        # The cap is an admitted active manifold, not a solver boundary.
        # Map only the unauthorized lower-chart approach into a symmetric
        # diagnostic fraction; exact cap contact maps to 1/2.
        fractions.extend(
            (point.retained_water_loading - chart.wet_retained_cap.lower_loading)
            / (
                chart.wet_retained_cap.loading_span
                + point.retained_water_loading
                - chart.wet_retained_cap.lower_loading
            )
            for point in certificate.points
        )
    if chart.dry_y_hexane is not None:
        fractions.extend(
            (value - chart.dry_y_hexane[0]) / (chart.dry_y_hexane[1] - chart.dry_y_hexane[0])
            for value in candidate.dry_y_hexane
        )
    else:
        activity_at_ref = sp.water_activity(
            chart.dry_config.pore.luikov.W_ref,
            chart.dry_config.pore.luikov,
        )
        for piece, (composition, pore) in enumerate(zip(candidate.dry_y_hexane, chart.dry_pores)):
            temperatures = _dry_composition_context_temperatures(
                piece,
                candidate.dry_temperatures_k,
                candidate.interface_temperature_k,
                chart,
            )
            lower, upper = _face_conditioned_composition_interval(
                temperatures,
                chart.dry_config.pressure_pa,
                pore,
            )
            fractions.append((composition - lower) / (upper - lower))
            # These are one-sided topology margins.  ``g/(1+g)`` maps the
            # exact positive dimensionless activity gap into (0,1/2), so the
            # existing symmetric distance calculation below treats only
            # collapse toward zero as a boundary approach.
            for temperature in temperatures:
                state = cp.evaluate_equilibrium(
                    temperature,
                    chart.dry_config.pressure_pa,
                    composition,
                    pore,
                )
                for gap in (
                    1.0 - state.water_activity,
                    1.0 - state.hexane_activity,
                    state.water_activity - activity_at_ref,
                ):
                    fractions.append(gap / (1.0 + gap))
    fractions_tuple = tuple(fractions)
    if not all(0.0 < value < 1.0 for value in fractions_tuple):
        raise RuntimeError("accepted bounded primitive touched its chart boundary")
    return fractions_tuple


__all__ = [
    "CutIntegratorLedger",
    "CutIntegratorState",
    "CutIntegratorStep",
    "CutIntegratorStepError",
    "CutResidualScales",
    "CutSolverControls",
    "CutStepSeed",
    "MaterialInventorySnapshot",
    "NoUsableDirectPhysicalFrontInteriorError",
    "P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING",
    "SameCellNonlinearSolverPolicy",
    "SameCellFrontCoordinatePolicy",
    "ScalarPostRootPolishAudit",
    "WaterPhaseInventorySnapshot",
    "WetRetainedWaterApplicabilityAudit",
    "WetRetainedWaterApplicabilityError",
    "advance_same_cell_backward_euler",
    "capacity_energy_scale",
    "certify_enabled_wet_retained_water_state_applicability",
    "certify_explicit_wet_retained_water_applicability",
    "initialize_integrator_state",
    "inventory_snapshot",
    "water_phase_inventory_snapshot",
]
