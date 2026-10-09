r"""Positive-core extinction-time targets and fail-closed limit certification.

The exact ``s=0`` endpoint has no interface equations and therefore cannot
determine its own event time.  This module localises time only on the strict
positive-core chart.  For a caller-supplied ``0 < z_* < z_old`` it holds
``z=z_*`` and replaces the usual front-coordinate unknown by the positive
step duration.  The residual is otherwise exactly
:func:`cut_transport.assemble_backward_euler`, including every wet, dry,
ALE, and Rankine--Hugoniot equation.

A finite extinction time is not inferred from nonlinear success.  A sequence
of positive-target roots must additionally have decreasing targets,
monotone time, contracting tail increments, a sufficiently small last target,
and a geometric tail bound below the declared tolerance.  Failure produces
no certificate and therefore cannot be handed to the exact-dry fixed-time
projection.  In particular, an asymptotically non-extinguishing front is a
valid target sequence but not a finite extinction event.  This is not a
solver-selected contact or pinned state: every positive-target root still
satisfies the moving-front equations.

The coefficient and Dirichlet cases remain numerical oracles and report
``physically_qualifying=False``.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy import optimize

from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_extinction_projection as projection
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_sparsity
from dtdc_simulator.core2.particle import cut_transport as cut


class PositiveCoreTargetStepError(RuntimeError):
    """Rejected target solve carrying the exact immutable rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState,
        *,
        nonlinear_evaluations: int = 0,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        last_candidate: cut.CutTransportUnknowns | None = None,
        last_event_duration_s: float | None = None,
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (
            ci.WetRetainedWaterApplicabilityAudit | None
        ) = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.last_candidate = last_candidate
        self.last_event_duration_s = last_event_duration_s
        self.retained_water_applicability_guard_failed = (
            retained_water_applicability_guard_failed
        )
        self.retained_water_applicability_audit = (
            retained_water_applicability_audit
        )


class ExtinctionTimeLimitError(RuntimeError):
    """The positive-target sequence did not certify a finite time limit."""

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState,
        *,
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (
            ci.WetRetainedWaterApplicabilityAudit | None
        ) = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.retained_water_applicability_guard_failed = (
            retained_water_applicability_guard_failed
        )
        self.retained_water_applicability_audit = (
            retained_water_applicability_audit
        )


@dataclass(frozen=True)
class EnvironmentConditionedResidualResolution:
    """Owner-ruled declared residual RESOLUTION for one certificate path.

    P-1 RULING (2026-08-22, option (a), verbatim "P-1 I approve your
    recommended option (a)"): the PART-02 extinction certificate carries a
    DECLARED environment-conditioned residual bound sized from fresh
    measurement of this workstation's positive-target residual floor (the
    BLAS/ULP-fragile cancellation band BC-7 section 4.4 documented); the
    frozen ``2e-11`` positive-target residual contract stays UNTOUCHED for
    every other consumer.

    This object is that declaration, typed.  It is NOT a tolerance change:
    ``frozen_contract_value`` must equal the state's own frozen contract
    exactly (a declaration written against another contract refuses), the
    declared bound must strictly exceed it (a declaration at or below the
    frozen contract is not an environment-conditioned resolution), and every
    ledger, certificate, and journal that accepts under it carries the
    declaration verbatim so no consumer can mistake the conditioned tier for
    the frozen one.  Only the residual-magnitude clause is conditioned; the
    condition-proxy, bounded-chart, ledger, and structural gates are
    untouched.  Refuse-not-clamp everywhere: residuals above the declared
    bound still refuse typed with exact rollback.

    P-2 RULING (2026-08-23, verbatim "P-2 approved"; floor record section 6
    recommendation): the same declaration MAY additionally carry a DECLARED
    EVALUATION BUDGET for the same single certificate path - the measured
    solver-termination floor remedy (the S7c escalated-retry pattern "the
    identical pipeline ... under the request's declared evaluation
    budget").  ``declared_maximum_function_evaluations`` and
    ``frozen_budget_value`` come together or not at all; the frozen value
    must equal the state's own ``maximum_function_evaluations`` exactly and
    the declared budget must strictly exceed it.  The frozen budget stays
    untouched for every other consumer; a declaration without the budget
    pair conditions the residual clause alone.
    """

    declared_bound: float
    frozen_contract_value: float
    ruling: str
    measurement: str
    declared_maximum_function_evaluations: int | None = None
    frozen_budget_value: int | None = None
    physically_qualifying: bool = field(default=False, init=False)
    widens_only_the_declared_certificate_path: bool = field(
        default=True, init=False
    )

    def __post_init__(self) -> None:
        for name, value in (
            ("declared bound", self.declared_bound),
            ("frozen contract value", self.frozen_contract_value),
        ):
            if not isinstance(value, float) or not math.isfinite(value) or (
                value <= 0.0
            ):
                raise ValueError(
                    f"environment-conditioned {name} must be a positive "
                    "finite binary64"
                )
        if self.declared_bound <= self.frozen_contract_value:
            raise ValueError(
                "an environment-conditioned resolution must strictly exceed "
                "the frozen contract it conditions; at or below it the frozen "
                "contract already governs and no declaration exists"
            )
        for name, value in (
            ("ruling", self.ruling),
            ("measurement", self.measurement),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    "an environment-conditioned resolution needs a non-empty "
                    f"{name} provenance string"
                )
        budget = self.declared_maximum_function_evaluations
        frozen_budget = self.frozen_budget_value
        if (budget is None) != (frozen_budget is None):
            raise ValueError(
                "a declared evaluation budget and the frozen budget it "
                "conditions come together or not at all"
            )
        if budget is not None:
            for name, value in (
                ("declared evaluation budget", budget),
                ("frozen budget value", frozen_budget),
            ):
                if isinstance(value, bool) or not isinstance(value, int) or (
                    value < 1
                ):
                    raise ValueError(
                        f"environment-conditioned {name} must be a positive "
                        "integer"
                    )
            if budget <= frozen_budget:
                raise ValueError(
                    "a declared evaluation budget must strictly exceed the "
                    "frozen budget it conditions; at or below it the frozen "
                    "budget already governs and no declaration exists"
                )


def _residual_limit_in_force(
    before: ci.CutIntegratorState,
    resolution: EnvironmentConditionedResidualResolution | None,
) -> float:
    """Return the residual limit governing one declared solve or validation.

    ``None`` is the frozen contract, bit-identical to the pre-declaration
    path.  A declaration must be the typed object above and must condition
    EXACTLY the state's own frozen contract; anything else refuses.
    """

    frozen_value = before.controls.nonlinear_residual_tolerance
    if resolution is None:
        return frozen_value
    if not isinstance(resolution, EnvironmentConditionedResidualResolution):
        raise TypeError(
            "the environment-conditioned resolution must be the typed "
            "EnvironmentConditionedResidualResolution declaration; got "
            f"{type(resolution).__name__}"
        )
    if resolution.frozen_contract_value != frozen_value:
        raise ValueError(
            "the environment-conditioned resolution declares the frozen "
            f"contract {resolution.frozen_contract_value!r} but this state's "
            f"frozen contract is {frozen_value!r}; a declaration written "
            "against another contract cannot condition this one"
        )
    return resolution.declared_bound


def _evaluation_budget_in_force(
    before: ci.CutIntegratorState,
    resolution: EnvironmentConditionedResidualResolution | None,
) -> int:
    """Return the nonlinear evaluation budget governing one declared solve.

    ``None`` - or a declaration without the P-2 budget pair - is the state's
    own frozen budget, bit-identical to the pre-declaration path.  A
    declared budget must condition EXACTLY the state's frozen budget;
    anything else refuses.  Callers validate the declaration's TYPE through
    :func:`_residual_limit_in_force` first.
    """

    frozen_budget = before.controls.maximum_function_evaluations
    if resolution is None or (
        resolution.declared_maximum_function_evaluations is None
    ):
        return frozen_budget
    if resolution.frozen_budget_value != frozen_budget:
        raise ValueError(
            "the environment-conditioned resolution declares the frozen "
            f"evaluation budget {resolution.frozen_budget_value!r} but this "
            f"state's frozen budget is {frozen_budget!r}; a declaration "
            "written against another budget cannot condition this one"
        )
    return resolution.declared_maximum_function_evaluations


@dataclass(frozen=True)
class PositiveCoreTargetSeed:
    """One target candidate, positive duration, and explicit Stefan scales."""

    candidate: cut.CutTransportUnknowns
    event_duration_s: float
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, cut.CutTransportUnknowns):
            raise TypeError("positive-core target seed needs cut unknowns")
        if not math.isfinite(self.event_duration_s) or self.event_duration_s <= 0.0:
            raise ValueError("positive-core target seed duration must be positive")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0
            for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every target Stefan scale must be positive and finite")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("positive-core target seed needs a provenance label")


@dataclass(frozen=True)
class PositiveCoreTargetRankAudit:
    """Structural proof that target ``z`` and unknown duration exchange rank."""

    residual_count: int
    ordinary_unknown_count: int
    fixed_target_state_unknown_count: int
    duration_unknown_count: int
    target_z_equation_count: int = field(default=0, init=False)
    endpoint_interface_equation_was_added: bool = field(default=False, init=False)
    uses_unchanged_positive_core_rh_rows: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.duration_unknown_count != 1:
            raise ValueError("positive-core target solve needs one duration unknown")
        if self.fixed_target_state_unknown_count + 1 != self.residual_count:
            raise ValueError("fixed-target state plus duration must be square")
        if self.ordinary_unknown_count != self.residual_count:
            raise ValueError("underlying positive-core residual must remain square")

    @property
    def target_unknown_count(self) -> int:
        return self.fixed_target_state_unknown_count + self.duration_unknown_count

    @property
    def is_square(self) -> bool:
        return self.target_unknown_count == self.residual_count


@dataclass(frozen=True)
class PositiveCoreTargetLedger:
    """One accepted positive-target nonlinear and structural certificate."""

    target_z: float
    target_fraction_of_initial_z: float
    event_duration_s: float
    nonlinear_evaluations: int
    nonlinear_message: str
    scaled_residuals: tuple[float, ...]
    maximum_scaled_residual: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    used_sparse_jacobian: bool
    jacobian_structural_nnz: int
    jacobian_structural_color_count: int
    rank_audit: PositiveCoreTargetRankAudit
    accepted: bool
    environment_conditioned_resolution: (
        EnvironmentConditionedResidualResolution | None
    ) = None
    event_time_is_finite_limit_certified: bool = field(default=False, init=False)
    endpoint_rh_equation_was_used: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class PositiveCoreTargetStep:
    """Uncommitted root from the same input state to one positive target."""

    before: ci.CutIntegratorState
    boundary: ct.PoreBoundary
    seed: PositiveCoreTargetSeed
    assembly: cut.CutTransportAssembly
    ledger: PositiveCoreTargetLedger
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class ExtinctionTimeRefinementControls:
    """Numerical contract for a finite positive-target time limit."""

    duration_bounds_s: tuple[float, float]
    minimum_levels: int = 4
    tail_increment_count: int = 3
    maximum_target_contraction_ratio: float = 0.75
    maximum_tail_increment_ratio: float = 0.5
    maximum_final_target_fraction: float = 1.0e-6
    absolute_time_tolerance_s: float = 1.0e-9
    relative_time_tolerance: float = 1.0e-8
    environment_conditioned_resolution: (
        EnvironmentConditionedResidualResolution | None
    ) = None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.environment_conditioned_resolution is not None and not (
            isinstance(
                self.environment_conditioned_resolution,
                EnvironmentConditionedResidualResolution,
            )
        ):
            raise TypeError(
                "refinement controls accept only the typed "
                "EnvironmentConditionedResidualResolution declaration or None"
            )
        lower, upper = self.duration_bounds_s
        if not all(math.isfinite(value) and value > 0.0 for value in (lower, upper)):
            raise ValueError("duration bounds must be positive and finite")
        if not lower < upper:
            raise ValueError("duration bounds must be ordered")
        for name, value, minimum in (
            ("minimum levels", self.minimum_levels, 3),
            ("tail increment count", self.tail_increment_count, 2),
        ):
            try:
                integer = operator.index(value)
            except TypeError as exc:
                raise ValueError(f"{name} must be an integer") from exc
            if isinstance(value, bool) or integer < minimum:
                raise ValueError(f"{name} must be at least {minimum}")
        if self.minimum_levels < self.tail_increment_count + 1:
            raise ValueError("minimum levels must cover every tail increment")
        for name, value in (
            ("target contraction", self.maximum_target_contraction_ratio),
            ("tail increment contraction", self.maximum_tail_increment_ratio),
        ):
            if not math.isfinite(value) or not 0.0 < value < 1.0:
                raise ValueError(f"maximum {name} ratio must lie in (0,1)")
        if not math.isfinite(self.maximum_final_target_fraction) or not (
            0.0 < self.maximum_final_target_fraction < 1.0
        ):
            raise ValueError("final target fraction must lie in (0,1)")
        if not math.isfinite(self.absolute_time_tolerance_s) or (
            self.absolute_time_tolerance_s <= 0.0
        ):
            raise ValueError("absolute time tolerance must be positive")
        if not math.isfinite(self.relative_time_tolerance) or not (
            0.0 < self.relative_time_tolerance <= 1.0e-4
        ):
            raise ValueError("relative time tolerance must lie in (0,1e-4]")


@dataclass(frozen=True)
class ExtinctionTimeLimitCertificate:
    """Numerical finite-limit certificate; never an endpoint closure."""

    before: ci.CutIntegratorState
    boundary: ct.PoreBoundary
    steps: tuple[PositiveCoreTargetStep, ...]
    extrapolated_event_duration_s: float
    last_positive_target_duration_s: float
    estimated_tail_bound_s: float
    maximum_observed_tail_increment_ratio: float
    controls: ExtinctionTimeRefinementControls
    uses_only_positive_core_equations: bool = field(default=True, init=False)
    endpoint_rh_equation_was_used: bool = field(default=False, init=False)
    exact_extinction_was_committed: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.steps:
            raise ValueError("finite-limit certificate needs positive-target roots")
        if self.estimated_tail_bound_s < 0.0:
            raise ValueError("finite-limit tail bound cannot be negative")
        if self.extrapolated_event_duration_s < self.last_positive_target_duration_s:
            raise ValueError("finite-limit extrapolation cannot precede its last target")


def audit_positive_core_target_rank(
    before: ci.CutIntegratorState,
) -> PositiveCoreTargetRankAudit:
    """Return the exact rank exchange on a cell-zero positive-core chart."""

    _require_cell_zero_positive_core(before)
    layout = before.transport.layout
    return PositiveCoreTargetRankAudit(
        residual_count=layout.residual_count,
        ordinary_unknown_count=layout.unknown_count,
        fixed_target_state_unknown_count=layout.unknown_count - 1,
        duration_unknown_count=1,
    )


def solve_positive_core_target(
    before: ci.CutIntegratorState,
    target_z: float,
    boundary: ct.PoreBoundary,
    seed: PositiveCoreTargetSeed,
    *,
    duration_bounds_s: tuple[float, float],
    environment_conditioned_resolution: (
        EnvironmentConditionedResidualResolution | None
    ) = None,
) -> PositiveCoreTargetStep:
    """Solve unchanged positive-core ALE/RH equations at fixed ``z_*``.

    ``environment_conditioned_resolution`` defaults to ``None`` - the frozen
    residual contract, bit-identical to the pre-declaration path.  When the
    P-1-ruled typed declaration is passed, ONLY the residual-magnitude
    acceptance clause reads the declared bound; the accepted ledger then
    carries the declaration verbatim, and downstream certification refuses
    any mixing of declared and undeclared evidence.
    """

    evaluations = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    last_candidate = None
    last_duration = None
    try:
        _require_cell_zero_positive_core(before)
        residual_limit = _residual_limit_in_force(
            before, environment_conditioned_resolution
        )
        evaluation_budget = _evaluation_budget_in_force(
            before, environment_conditioned_resolution
        )
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("target solve requires a supported pore boundary")
        _validate_duration_bounds(duration_bounds_s)
        old_z = before.transport.geometry.front.z
        if not math.isfinite(target_z) or not 0.0 < target_z < old_z:
            raise ValueError("positive target z must lie strictly inside (0,z_old)")
        layout = before.transport.layout
        rank = audit_positive_core_target_rank(before)
        if seed.candidate.front_z != target_z:
            raise ValueError("target seed front z must equal the supplied target exactly")
        if len(seed.candidate.vector()) != layout.unknown_count:
            raise ValueError("target seed does not match the positive-core rank")
        if len(seed.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
            raise ValueError("one target Stefan scale is required per dry face")
        if not _inside(seed.event_duration_s, duration_bounds_s):
            raise ValueError("target seed duration left its declared open interval")

        chart = ci._coordinate_chart(before)  # noqa: SLF001
        coordinates = ci._encode_candidate(  # noqa: SLF001
            seed.candidate,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        coordinates[-2] = ci._encode_open(  # noqa: SLF001
            seed.event_duration_s,
            duration_bounds_s,
        )
        initial = cut.assemble_backward_euler(
            before.transport,
            seed.candidate,
            seed.event_duration_s,
            boundary,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        scales = ci._residual_scales(  # noqa: SLF001
            before.transport,
            initial,
            seed.event_duration_s,
        )
        scale_vector = np.asarray(scales.vector, dtype=float)
        structural = (
            cut_sparsity.build_cut_jacobian_sparsity(layout)
            if layout.unknown_count > before.controls.maximum_dense_unknowns
            else None
        )

        def decode(values: Sequence[float]):
            duration = ci._decode_open(  # noqa: SLF001
                float(values[-2]),
                duration_bounds_s,
            )
            candidate_coordinates = np.asarray(values, dtype=float).copy()
            candidate_coordinates[-2] = ci._encode_open(  # noqa: SLF001
                target_z,
                chart.front_z,
            )
            candidate = ci._decode_candidate(  # noqa: SLF001
                candidate_coordinates,
                layout,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
            )
            if candidate.front_z != target_z:
                candidate = cut.CutTransportUnknowns(
                    wet_temperatures_k=candidate.wet_temperatures_k,
                    wet_retained_water_loadings=(
                        candidate.wet_retained_water_loadings
                    ),
                    dry_temperatures_k=candidate.dry_temperatures_k,
                    dry_y_hexane=candidate.dry_y_hexane,
                    dry_total_stefan_fluxes_mol_m2_s=(
                        candidate.dry_total_stefan_fluxes_mol_m2_s
                    ),
                    front_z=target_z,
                    interface_temperature_k=candidate.interface_temperature_k,
                )
            return candidate, duration

        def residual(values: np.ndarray) -> np.ndarray:
            nonlocal evaluations
            candidate, duration = decode(values)
            assembly = cut.assemble_backward_euler(
                before.transport,
                candidate,
                duration,
                boundary,
                enforce_reduced_film_thresholds=False,
            )
            evaluations += 1
            return np.asarray(
                ci._datum_covariant_residual_vector(assembly.residuals),  # noqa: SLF001
                dtype=float,
            ) / scale_vector

        lower, upper = ci._coordinate_solver_bounds(  # noqa: SLF001
            layout,
            before.controls.maximum_logit_magnitude,
        )
        solution = optimize.least_squares(
            residual,
            coordinates,
            method="dogbox" if structural is not None else "trf",
            bounds=(lower, upper),
            ftol=before.controls.nonlinear_step_tolerance,
            xtol=before.controls.nonlinear_step_tolerance,
            gtol=before.controls.nonlinear_step_tolerance,
            max_nfev=evaluation_budget,
            x_scale="jac",
            jac_sparsity=None if structural is None else structural.jac_sparsity,
            tr_options=(
                {"atol": 1.0e-14, "btol": 1.0e-14, "maxiter": 1000}
                if structural is not None
                else None
            ),
        )
        candidate, duration = decode(solution.x)
        last_candidate = candidate
        last_duration = duration
        converged = cut.assemble_backward_euler(
            before.transport,
            candidate,
            duration,
            boundary,
        )
        evaluations += 1
        signed_scaled = tuple(
            value / scale
            for value, scale in zip(
                ci._datum_covariant_residual_vector(converged.residuals),  # noqa: SLF001
                scales.vector,
            )
        )
        maximum_scaled = max(abs(value) for value in signed_scaled)
        condition_proxy = ci._jacobian_condition_proxy(solution.jac)  # noqa: SLF001
        fractions = (*ci._bounded_fractions(candidate, chart), _fraction(duration, duration_bounds_s))  # noqa: SLF001
        minimum_fraction = min(min(value, 1.0 - value) for value in fractions)
        if not solution.success:
            raise RuntimeError(f"positive-target root did not succeed: {solution.message}")
        if np.any(solution.active_mask != 0):
            raise RuntimeError("positive-target root touched the numerical logit guard")
        if maximum_scaled > residual_limit:
            if environment_conditioned_resolution is None:
                raise RuntimeError(
                    "positive-target residual contract failed: "
                    f"{maximum_scaled:.3e} > "
                    f"{before.controls.nonlinear_residual_tolerance:.3e}"
                )
            raise RuntimeError(
                "positive-target residual exceeds even the declared "
                "environment-conditioned resolution: "
                f"{maximum_scaled:.3e} > declared bound "
                f"{environment_conditioned_resolution.declared_bound:.3e} "
                "(frozen contract "
                f"{before.controls.nonlinear_residual_tolerance:.3e})"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > before.controls.maximum_condition_proxy
        ):
            raise RuntimeError(
                "positive-target Jacobian condition proxy is uncertified: "
                f"{condition_proxy:.3e}"
            )
        ledger = PositiveCoreTargetLedger(
            target_z=target_z,
            target_fraction_of_initial_z=target_z / old_z,
            event_duration_s=duration,
            nonlinear_evaluations=evaluations,
            nonlinear_message=str(solution.message),
            scaled_residuals=signed_scaled,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            minimum_fractional_distance_to_bound=minimum_fraction,
            used_sparse_jacobian=structural is not None,
            jacobian_structural_nnz=(
                layout.unknown_count**2 if structural is None else structural.nnz
            ),
            jacobian_structural_color_count=(
                layout.unknown_count if structural is None else structural.color_count
            ),
            rank_audit=rank,
            accepted=True,
            environment_conditioned_resolution=(
                environment_conditioned_resolution
            ),
        )
        ci.certify_enabled_wet_retained_water_state_applicability(
            before.controls,
            before.transport,
        )
        if before.controls.enforce_p1ef_smooth_wet_retained_water_guard:
            ci.certify_explicit_wet_retained_water_applicability(
                before.controls,
                converged.candidate.wet_retained_water_loadings,
                converged.candidate.effective_wet_retained_water_capacity_duals_over_rt,
                before.transport.config.wet.luikov,
            )
        return PositiveCoreTargetStep(before, boundary, seed, converged, ledger)
    except PositiveCoreTargetStepError:
        raise
    except Exception as exc:
        if not isinstance(before, ci.CutIntegratorState):
            raise
        applicability_error = (
            exc
            if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
            else None
        )
        raise PositiveCoreTargetStepError(
            f"positive-core target solve rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            last_candidate=last_candidate,
            last_event_duration_s=last_duration,
            retained_water_applicability_guard_failed=(
                applicability_error is not None
            ),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


def continuation_seed(
    previous: PositiveCoreTargetStep,
    target_z: float,
    *,
    label: str,
) -> PositiveCoreTargetSeed:
    """Create a geometry-scaled next seed without committing the trial root."""

    old_z = previous.before.transport.geometry.front.z
    if not 0.0 < target_z < previous.ledger.target_z:
        raise ValueError("continuation target must be positive and strictly smaller")
    old_sweep = old_z - previous.ledger.target_z
    new_sweep = old_z - target_z
    duration = previous.ledger.event_duration_s * new_sweep / old_sweep
    return PositiveCoreTargetSeed(
        candidate=cut.CutTransportUnknowns(
            wet_temperatures_k=previous.assembly.candidate.wet_temperatures_k,
            wet_retained_water_loadings=(
                previous.assembly.candidate.wet_retained_water_loadings
            ),
            dry_temperatures_k=previous.assembly.candidate.dry_temperatures_k,
            dry_y_hexane=previous.assembly.candidate.dry_y_hexane,
            dry_total_stefan_fluxes_mol_m2_s=(
                previous.assembly.candidate.dry_total_stefan_fluxes_mol_m2_s
            ),
            front_z=target_z,
            interface_temperature_k=(
                previous.assembly.candidate.interface_temperature_k
            ),
        ),
        event_duration_s=duration,
        stefan_flux_scales_mol_m2_s=(
            previous.seed.stefan_flux_scales_mol_m2_s
        ),
        label=label,
    )


def certify_finite_extinction_time(
    steps: Sequence[PositiveCoreTargetStep],
    controls: ExtinctionTimeRefinementControls,
) -> ExtinctionTimeLimitCertificate:
    """Certify a numerically resolved finite ``z_* -> 0+`` time limit."""

    roots = tuple(steps)
    if not roots:
        raise ValueError("finite-limit audit needs positive-target roots")
    candidate_before = getattr(roots[0], "before", None)
    if not isinstance(candidate_before, ci.CutIntegratorState):
        raise TypeError(
            "finite-limit audit needs genuine PositiveCoreTargetStep roots"
        )
    before = candidate_before
    try:
        if not isinstance(controls, ExtinctionTimeRefinementControls):
            raise TypeError("finite-limit audit needs typed refinement controls")
        if len(roots) < controls.minimum_levels:
            raise ValueError("too few positive-target refinement levels")
        boundary = getattr(roots[0], "boundary", None)
        for step in roots:
            validate_positive_core_target_step(
                step,
                before,
                boundary,
                controls,
            )
        targets = tuple(step.ledger.target_z for step in roots)
        durations = tuple(step.ledger.event_duration_s for step in roots)
        if any(new >= old for old, new in zip(targets, targets[1:])):
            raise ValueError("positive targets must decrease strictly toward zero")
        if any(new <= old for old, new in zip(durations, durations[1:])):
            raise ValueError("receding target times must increase strictly")
        contractions = tuple(new / old for old, new in zip(targets, targets[1:]))
        tail_contractions = contractions[-controls.tail_increment_count :]
        if any(
            value > controls.maximum_target_contraction_ratio
            for value in tail_contractions
        ):
            raise ValueError("tail targets were not refined by the declared contraction")
        initial_z = before.transport.geometry.front.z
        final_fraction = targets[-1] / initial_z
        if final_fraction > controls.maximum_final_target_fraction:
            raise ValueError("last positive target is not close enough to zero")

        increments = tuple(
            new - old for old, new in zip(durations, durations[1:])
        )
        tail = increments[-controls.tail_increment_count :]
        ratios = tuple(new / old for old, new in zip(tail, tail[1:]))
        maximum_ratio = max(ratios)
        if maximum_ratio > controls.maximum_tail_increment_ratio:
            raise ValueError(
                "tail duration increments do not contract; finite extinction "
                "is not certified (the branch may be asymptotically "
                "non-extinguishing)"
            )
        remainder = tail[-1] * maximum_ratio / (1.0 - maximum_ratio)
        extrapolated = durations[-1] + remainder
        tolerance = max(
            controls.absolute_time_tolerance_s,
            controls.relative_time_tolerance * extrapolated,
        )
        if remainder > tolerance:
            raise ValueError(
                "positive-target time tail exceeds the declared finite-limit tolerance"
            )
        return ExtinctionTimeLimitCertificate(
            before=before,
            boundary=boundary,
            steps=roots,
            extrapolated_event_duration_s=extrapolated,
            last_positive_target_duration_s=durations[-1],
            estimated_tail_bound_s=remainder,
            maximum_observed_tail_increment_ratio=maximum_ratio,
            controls=controls,
        )
    except ExtinctionTimeLimitError:
        raise
    except Exception as exc:
        applicability_error = (
            exc
            if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
            else None
        )
        raise ExtinctionTimeLimitError(
            f"finite extinction-time limit rejected with exact rollback: {exc}",
            before,
            retained_water_applicability_guard_failed=(
                applicability_error is not None
            ),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


def validate_positive_core_target_step(
    step: PositiveCoreTargetStep,
    before: ci.CutIntegratorState,
    boundary: ct.PoreBoundary,
    controls: ExtinctionTimeRefinementControls,
) -> None:
    """Validate that one limit root is a coherent live target solve.

    This is deliberately stricter than duck typing.  The limit certificate is
    allowed to carry only the immutable payload returned by
    :func:`solve_positive_core_target`; manufactured metadata and edited
    ledgers are not event evidence.  The checks duplicate inexpensive
    identities from the nonlinear assembly so the extinction orchestrator can
    reject a stale or tampered certificate before projection work begins.
    """

    if not isinstance(step, PositiveCoreTargetStep):
        raise TypeError(
            "finite-limit audit needs genuine PositiveCoreTargetStep roots"
        )
    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("positive-target validation needs CutIntegratorState")
    if not isinstance(
        boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("positive-target validation needs one supported boundary")
    if not isinstance(controls, ExtinctionTimeRefinementControls):
        raise TypeError("positive-target validation needs typed refinement controls")
    if not isinstance(step.seed, PositiveCoreTargetSeed):
        raise TypeError("positive-target root lost its typed nonlinear seed")
    if not isinstance(step.assembly, cut.CutTransportAssembly):
        raise TypeError("positive-target root lost its live transport assembly")
    if not isinstance(step.ledger, PositiveCoreTargetLedger):
        raise TypeError("positive-target root lost its typed solver ledger")
    if not isinstance(step.ledger.rank_audit, PositiveCoreTargetRankAudit):
        raise TypeError("positive-target root lost its typed rank audit")
    if step.before is not before:
        raise ValueError("every target root must share the exact rollback state")
    if step.boundary != boundary:
        raise ValueError("every target root must share one boundary oracle")

    assembly = step.assembly
    ledger = step.ledger
    if assembly.before is not before.transport:
        raise ValueError("positive-target assembly changed its immutable input state")
    if assembly.surface_boundary != boundary:
        raise ValueError("positive-target assembly changed its boundary oracle")
    if assembly.target_config != before.transport.config:
        raise ValueError("positive-target assembly changed its transport authority")
    if assembly.layout != before.transport.layout:
        raise ValueError("positive-target assembly changed the cell-zero chart rank")
    if assembly.swept_geometry.previous is not before.transport.geometry:
        raise ValueError("positive-target sweep changed its immutable old geometry")
    if assembly.swept_geometry.current is not assembly.candidate_geometry:
        raise ValueError("positive-target sweep and candidate geometry disagree")

    target = ledger.target_z
    initial_z = before.transport.geometry.front.z
    target_values = (
        step.seed.candidate.front_z,
        assembly.candidate.front_z,
        assembly.candidate_geometry.front.z,
    )
    if not math.isfinite(target) or not 0.0 < target < initial_z:
        raise ValueError("positive-target root target left (0,z_old)")
    if any(value != target for value in target_values):
        raise ValueError("positive-target seed, assembly, and ledger targets disagree")
    if assembly.candidate_geometry.front.regime != "partial":
        raise ValueError("positive-target root left the strict partial-front chart")
    if assembly.candidate_geometry.cut_cell_index != 0:
        raise ValueError("positive-target root left master cell zero")
    if ledger.target_fraction_of_initial_z != target / initial_z:
        raise ValueError("positive-target fraction is inconsistent with its state")
    if assembly.swept_geometry.wet_to_dry_volume_m3 <= 0.0:
        raise ValueError("positive-target root did not perform a receding wet-core sweep")

    duration = ledger.event_duration_s
    if not _inside(duration, controls.duration_bounds_s):
        raise ValueError("solved positive-target duration left its declared interval")
    if assembly.swept_geometry.dt_s != duration:
        raise ValueError("positive-target assembly and ledger durations disagree")
    if not _inside(step.seed.event_duration_s, controls.duration_bounds_s):
        raise ValueError("positive-target seed duration left its declared interval")

    expected_rank = audit_positive_core_target_rank(before)
    if ledger.rank_audit != expected_rank:
        raise ValueError("positive-target rank audit is stale or inconsistent")
    if len(assembly.candidate.vector()) != assembly.layout.unknown_count:
        raise ValueError("positive-target candidate does not match its rank")
    raw_residuals = tuple(assembly.residuals.vector)
    scaled_residuals = tuple(ledger.scaled_residuals)
    if len(raw_residuals) != assembly.layout.residual_count or len(
        scaled_residuals
    ) != assembly.layout.residual_count:
        raise ValueError("positive-target residual vector does not match its rank")
    if not all(math.isfinite(value) for value in (*raw_residuals, *scaled_residuals)):
        raise ValueError("positive-target residual evidence must be finite")
    maximum_scaled = max(abs(value) for value in scaled_residuals)
    if maximum_scaled != ledger.maximum_scaled_residual:
        raise ValueError("positive-target maximum scaled residual is inconsistent")
    declared = controls.environment_conditioned_resolution
    if ledger.environment_conditioned_resolution != declared:
        raise ValueError(
            "positive-target root and refinement controls disagree about the "
            "environment-conditioned residual declaration; declared and "
            "undeclared evidence never mix in one certificate"
        )
    if maximum_scaled > _residual_limit_in_force(before, declared):
        if declared is None:
            raise ValueError(
                "positive-target root exceeds the nonlinear tolerance"
            )
        raise ValueError(
            "positive-target root exceeds even the declared "
            "environment-conditioned resolution"
        )

    evaluations = ledger.nonlinear_evaluations
    if (
        isinstance(evaluations, bool)
        or not isinstance(evaluations, int)
        or evaluations <= 0
    ):
        raise ValueError("positive-target root must disclose nonzero nonlinear work")
    if not ledger.accepted:
        raise ValueError("a rejected target cannot enter the limit audit")
    if not math.isfinite(ledger.condition_proxy) or not (
        0.0 <= ledger.condition_proxy <= before.controls.maximum_condition_proxy
    ):
        raise ValueError("positive-target condition proxy is not certified")
    if not math.isfinite(ledger.minimum_fractional_distance_to_bound) or not (
        0.0 < ledger.minimum_fractional_distance_to_bound <= 0.5
    ):
        raise ValueError("positive-target bounded-chart distance is not interior")
    ci.certify_enabled_wet_retained_water_state_applicability(
        before.controls,
        before.transport,
    )
    if getattr(
        before.controls,
        "enforce_p1ef_smooth_wet_retained_water_guard",
        False,
    ):
        ci.certify_explicit_wet_retained_water_applicability(
            before.controls,
            assembly.candidate.wet_retained_water_loadings,
            assembly.candidate.effective_wet_retained_water_capacity_duals_over_rt,
            before.transport.config.wet.luikov,
        )


def project_certified_extinction(
    certificate: ExtinctionTimeLimitCertificate,
    seed: projection.ExtinctionProjectionSeed,
) -> projection.ExtinctionProjectionStep:
    """Hand a certified duration to the separate exact-dry projection."""

    if not isinstance(certificate, ExtinctionTimeLimitCertificate):
        raise TypeError("exact-dry handoff requires a finite-limit certificate")
    duration = certificate.extrapolated_event_duration_s
    if seed.candidate.event_time_s != duration:
        raise ExtinctionTimeLimitError(
            "projection seed duration does not exactly match the certified limit",
            certificate.before,
        )
    return projection.solve_fixed_time_extinction_projection(
        certificate.before,
        duration,
        certificate.boundary,
        seed,
        environment_conditioned_resolution=(
            certificate.controls.environment_conditioned_resolution
        ),
    )


def _require_cell_zero_positive_core(before: ci.CutIntegratorState) -> None:
    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("positive-core target solve requires CutIntegratorState")
    geometry = before.transport.geometry
    if geometry.front.regime != "partial" or not 0.0 < geometry.front.z < 1.0:
        raise ValueError("extinction localisation requires a strict positive core")
    if before.transport.layout.cut_cell_index != 0:
        raise ValueError("extinction localisation begins only in master cell zero")
    ci.certify_enabled_wet_retained_water_state_applicability(
        before.controls,
        before.transport,
    )


def _validate_duration_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not all(math.isfinite(value) and value > 0.0 for value in bounds) or not (
        lower < upper
    ):
        raise ValueError("target duration bounds must be positive, finite, and ordered")


def _inside(value: float, bounds: tuple[float, float]) -> bool:
    return math.isfinite(value) and bounds[0] < value < bounds[1]


def _fraction(value: float, bounds: tuple[float, float]) -> float:
    return (value - bounds[0]) / (bounds[1] - bounds[0])


__all__ = [
    "EnvironmentConditionedResidualResolution",
    "ExtinctionTimeLimitCertificate",
    "ExtinctionTimeLimitError",
    "ExtinctionTimeRefinementControls",
    "PositiveCoreTargetLedger",
    "PositiveCoreTargetRankAudit",
    "PositiveCoreTargetSeed",
    "PositiveCoreTargetStep",
    "PositiveCoreTargetStepError",
    "audit_positive_core_target_rank",
    "certify_finite_extinction_time",
    "continuation_seed",
    "project_certified_extinction",
    "solve_positive_core_target",
    "validate_positive_core_target_step",
]
