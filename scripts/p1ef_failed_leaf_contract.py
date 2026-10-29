"""Exact release verifier for P1EF failed-leaf and physical-event partitions.

This module is deliberately reporting-only.  It does not call a solver or
construct a trajectory.  Both the process-isolated case runner and the final
consolidator use it to recompute Amendment-09 compliance from disclosed raw
attempt and accepted-leaf records instead of trusting a producer boolean.
Physical event intervals are verified independently and never inherit the
failed-leaf controller's binary-grid or absolute-floor restrictions.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence


POLICY_ID = "frozen_amendment09_absolute_floor_failed_leaf_only"
ABSOLUTE_MINIMUM_LEAF_DT_S = 0.00029296875
MAXIMUM_BINARY_DEPTH_CAP = 8
PROVISIONAL_FUNCTION_EVALUATIONS_PER_SEED = 16
REVERSAL_MAXIMUM_SEED_ATTEMPTS = 2
NONREVERSAL_MAXIMUM_SEED_ATTEMPTS = 3
PREDICTOR_BACKTRACK_MAXIMUM_BINARY_EXPONENT = 8
LEDGER_LIMIT = 1.0e-10
TIME_ABS_TOLERANCE_S = 1.0e-15
FACE_CONDITIONED_PREDICTOR_AMENDMENT_ID = "GT-PS-2-P1E-17"
FACE_CONDITIONED_PREDICTOR_COORDINATE = (
    "exact_face_conditioned_gas_only_barrier_coordinate"
)
FACE_CONDITIONED_PREDICTOR_CONTRACT_KEY = (
    "face_conditioned_dry_composition_predictor"
)
INTERFACE_ROOT_AMENDMENT_ID = "GT-PS-2-P1E-18"
INTERFACE_ROOT_CONTRACT_KEY = "interface_saturation_endpoint_selection"
DRY_ENERGY_CAPACITY_AMENDMENT_ID = "GT-PS-2-P1E-19"
DRY_ENERGY_CAPACITY_CONTRACT_KEY = "dry_energy_capacity_normalization_stencil"

_REQUESTED_ROLE = "requested_backward_euler"
_ANCHOR_ROLE = "current_state_anchored_front_seed_fallback_backward_euler"
_PREDICTOR_CHART_VALIDATION_ROLE = "predictor_chart_validation"
_ANCHOR_CHART_VALIDATION_ROLE = (
    "current_state_anchored_front_seed_fallback_chart_validation"
)
_SEED_BUDGET_EXHAUSTED_ROLE = "provisional_nonlinear_seed_attempt_budget_exhausted"
_NONLINEAR_ROLES = {_REQUESTED_ROLE, _ANCHOR_ROLE}
_NONLINEAR_EXEMPT_ROLES = {
    _PREDICTOR_CHART_VALIDATION_ROLE,
    _ANCHOR_CHART_VALIDATION_ROLE,
    _SEED_BUDGET_EXHAUSTED_ROLE,
}
_NORMALIZED_LEDGER_FIELDS = (
    "normalized_water_step_residual",
    "normalized_hexane_step_residual",
    "normalized_energy_step_residual",
    "normalized_water_cumulative_residual",
    "normalized_hexane_cumulative_residual",
    "normalized_energy_cumulative_residual",
)
_RAW_MACRO_LEDGER_FIELDS = (
    "raw_macro_water_residual_mol",
    "raw_macro_hexane_residual_mol",
    "raw_macro_energy_residual_j",
)
_EVENT_PARTITION_ROLES = (
    "exact_face_arrival",
    "exact_face_departure_remainder",
)
_EVENT_ORDERING_MODES = {
    "rigorous_interval_certificate",
    "open_interval_numerical_ordering",
}
_EVENT_CHAIN_FIELDS = (
    "macro_input_identity_preserved",
    "event_input_identity_preserved",
    "event_to_remainder_identity_preserved",
    "macro_output_identity_preserved",
    "boundary_identity_preserved",
)
_EVENT_RAW_LEDGER_FIELDS = (
    "water_step_residual_mol",
    "hexane_step_residual_mol",
    "energy_step_residual_j",
    "water_cumulative_residual_mol",
    "hexane_cumulative_residual_mol",
    "energy_cumulative_residual_j",
)


def _finite_number(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _close(left: Any, right: float) -> bool:
    return _finite_number(left) and math.isclose(
        float(left),
        right,
        rel_tol=0.0,
        abs_tol=TIME_ABS_TOLERANCE_S,
    )


def _integer(value: Any, *, lower: int | None = None) -> bool:
    if isinstance(value, bool) or not isinstance(value, int):
        return False
    return lower is None or value >= lower


def _path_interval(path: str) -> tuple[float, float]:
    numerator = 0
    for symbol in path:
        numerator = 2 * numerator + (symbol == "R")
    denominator = 2 ** len(path)
    return numerator / denominator, (numerator + 1) / denominator


def effective_maximum_binary_depth(requested_dt_s: float) -> int:
    """Return the exact Amendment-09 depth for one requested macro width."""

    if not _finite_number(requested_dt_s) or float(requested_dt_s) <= 0.0:
        raise ValueError("requested macro width must be positive and finite")
    eligible = tuple(
        depth
        for depth in range(1, MAXIMUM_BINARY_DEPTH_CAP + 1)
        if math.ldexp(float(requested_dt_s), -depth) >= ABSOLUTE_MINIMUM_LEAF_DT_S
    )
    return max(eligible, default=0)


def floor_contract(requested_dt_s: float) -> dict[str, Any]:
    """Reconstruct the immutable producer declaration for one macro width."""

    requested = float(requested_dt_s)
    depth = effective_maximum_binary_depth(requested)
    terminal = math.ldexp(requested, -depth) if depth else requested
    next_leaf = math.ldexp(requested, -(depth + 1))
    return {
        "requested_dt_s": requested,
        "effective_maximum_binary_depth": depth,
        "effective_terminal_leaf_dt_s": terminal,
        "maximum_leaf_count": 2**depth if depth else 1,
        "adaptive_recovery_available": depth >= 1,
        "no_leaf_below_absolute_floor": (
            depth == 0 or terminal >= ABSOLUTE_MINIMUM_LEAF_DT_S
        ),
        "floor_was_binding": (
            depth >= 1 and next_leaf < ABSOLUTE_MINIMUM_LEAF_DT_S
        ),
        "depth_cap_was_binding": (
            depth == MAXIMUM_BINARY_DEPTH_CAP
            and next_leaf >= ABSOLUTE_MINIMUM_LEAF_DT_S
        ),
    }


def face_conditioned_predictor_contract() -> dict[str, Any]:
    """Return the exact Amendment-17 history-predictor declaration."""

    return {
        "numerical_amendment_id": FACE_CONDITIONED_PREDICTOR_AMENDMENT_ID,
        "gas_only_predictor_coordinate": FACE_CONDITIONED_PREDICTOR_COORDINATE,
        "face_conditioned_coordinate_applies_only_to_gas_only_primitive_domain": True,
        "source_surface_boundary_included_in_source_chart": True,
        "target_surface_boundary_included_in_target_chart": True,
        "reduced_film_conditioning_temperature_source": (
            "declared_bulk_boundary_temperature_never_solved_surface_temperature"
        ),
        "predicted_temperatures_constructed_before_target_composition_decode": True,
        "joint_primitive_band_raw_mole_fraction_predictor_unchanged": True,
        "accepted_face_event_predictor_is_front_only": True,
        "composition_clipping_or_projection_allowed": False,
        "accepted_equations_changed": False,
        "acceptance_tolerances_changed": False,
        "material_inventory_changed": False,
        "solver_work_budget_changed": False,
        "face_event_solver_or_partition_changed": False,
        "component_energy_ledger_or_acceptance_limit_changed": False,
        "authoritative_target_seed_validation_required": True,
    }


def interface_root_contract() -> dict[str, Any]:
    """Return the exact Amendment-18 interface-root declaration."""

    return {
        "numerical_amendment_id": INTERFACE_ROOT_AMENDMENT_ID,
        "phy039_log_fugacity_residual_unchanged": True,
        "log_fugacity_tolerance_unchanged": True,
        "caller_declared_composition_bracket_unchanged": True,
        "actual_consumed_shell_oil_label_used": True,
        "product_form_hexane_activity_is_gas_side_authority": True,
        "full_coupled_pore_equilibrium_certified_inside_bisection": True,
        "representability_certified_inside_log_tolerance_plateau": True,
        "maximum_representable_endpoint_selected_after_roundoff_plateau": True,
        "ordinary_nonplateau_root_returns_first_fully_certified_in_tolerance_candidate": True,
        "hexane_activity_above_one_classified_as_upper_endpoint": True,
        "nonhexane_coupled_pore_failures_are_fatal": True,
        "nextafter_retreat_after_root_allowed": False,
        "post_root_composition_clipping_or_projection_allowed": False,
        "closed_saturated_interface_distinct_from_open_dry_cell_chart": True,
        "open_dry_cell_chart_used_to_validate_closed_interface_endpoint": False,
        "accepted_equations_changed": False,
        "acceptance_tolerances_changed": False,
        "solver_work_budget_changed": False,
        "face_event_solver_or_partition_changed": False,
        "material_inventory_changed": False,
        "component_energy_ledger_or_acceptance_limit_changed": False,
    }


def dry_energy_capacity_stencil_contract() -> dict[str, Any]:
    """Return the exact Amendment-19 normalization-only declaration."""

    return {
        "numerical_amendment_id": DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        "applies_to_exact_gas_only_primitive_domain": True,
        "stencil_selection_order": [
            "second_order_centered",
            "second_order_forward",
            "second_order_backward",
            "reject_without_nonlinear_iteration",
        ],
        "joint_primitive_band_centered_route_bit_identical": True,
        "full_coupled_pore_certification_required_for_every_sample": True,
        "finite_positive_capacity_required": True,
        "same_helper_used_by_same_cell_birth_face_arrival_and_face_departure": True,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_allowed": False,
        "property_extrapolation_allowed": False,
        "acceptance_tolerances_changed": False,
        "every_accepted_residual_scale_evaluation_must_be_serialized": True,
        "physically_qualifying": False,
    }


def solver_contract_errors(
    value: Any,
    *,
    expected_timesteps_s: Sequence[float],
) -> list[str]:
    """Verify the immutable controller declaration in a raw backend payload."""

    errors: list[str] = []
    if not isinstance(value, Mapping):
        return ["solver acceptance contract is missing"]
    expected_grids = list(expected_timesteps_s)
    if value.get("backward_euler_grids_s") != expected_grids:
        errors.append("solver requested macro grids differ from the formal campaign")
    if value.get("accepted_subdivision_policy") != POLICY_ID:
        errors.append("solver does not declare the frozen failed-leaf policy")
    if value.get("arbitrary_accepted_subdivision_allowed") is not False:
        errors.append("solver permits arbitrary accepted subdivision")
    if "accepted_subdivision_allowed" in value:
        errors.append("solver retains the superseded no-subdivision contract")
    predictor = value.get(FACE_CONDITIONED_PREDICTOR_CONTRACT_KEY)
    if not isinstance(predictor, Mapping):
        errors.append("solver face-conditioned predictor declaration is missing")
    elif dict(predictor) != face_conditioned_predictor_contract():
        errors.append("solver face-conditioned predictor declaration is not exact")
    interface_root = value.get(INTERFACE_ROOT_CONTRACT_KEY)
    if not isinstance(interface_root, Mapping):
        errors.append("solver interface-root declaration is missing")
    elif dict(interface_root) != interface_root_contract():
        errors.append("solver interface-root declaration is not exact")
    dry_energy_capacity = value.get(DRY_ENERGY_CAPACITY_CONTRACT_KEY)
    if not isinstance(dry_energy_capacity, Mapping):
        errors.append("solver dry-energy capacity declaration is missing")
    elif dict(dry_energy_capacity) != dry_energy_capacity_stencil_contract():
        errors.append("solver dry-energy capacity declaration is not exact")
    controller = value.get("failed_leaf_controller")
    expected = {
        "policy_id": POLICY_ID,
        "licensed_only_for_rejected_same_cell_first_material_step_after_imposed_boundary_discontinuity": True,
        "absolute_minimum_leaf_dt_s": ABSOLUTE_MINIMUM_LEAF_DT_S,
        "maximum_binary_depth_cap": MAXIMUM_BINARY_DEPTH_CAP,
        "formal_grid_effective_depths": [
            floor_contract(value) for value in expected_timesteps_s
        ],
        "failed_leaf_only_left_to_right_binary_recursion": True,
        "provisional_optimizer_function_evaluations_per_ordered_seed": (
            PROVISIONAL_FUNCTION_EVALUATIONS_PER_SEED
        ),
        "reversal_maximum_nonlinear_seed_attempts": (
            REVERSAL_MAXIMUM_SEED_ATTEMPTS
        ),
        "nonreversal_maximum_nonlinear_seed_attempts": (
            NONREVERSAL_MAXIMUM_SEED_ATTEMPTS
        ),
        "predictor_backtrack_maximum_binary_exponent": (
            PREDICTOR_BACKTRACK_MAXIMUM_BINARY_EXPONENT
        ),
        "requested_macro_nodes_remain_exact": True,
        "accepted_internal_leaf_sizes_and_nodes_disclosed": True,
        "atomic_macrostep_commit": True,
        "exact_rollback_on_rejection": True,
        "raw_macro_component_energy_ledger_limit": LEDGER_LIMIT,
        "wall_clock_timeout_controls_acceptance": False,
        "per_seed_and_per_leaf_work_budget_unchanged_from_amendment03": True,
        "total_tree_work_envelope_changed_by_amendment09": True,
    }
    if not isinstance(controller, Mapping):
        errors.append("solver failed-leaf controller declaration is missing")
    elif dict(controller) != expected:
        errors.append("solver failed-leaf controller declaration is not exact")
    return errors


def _controller_evidence_errors(
    value: Any,
    *,
    adaptive: bool,
    requested_dt_s: float,
    leaf_durations: Sequence[Any],
) -> list[str]:
    errors: list[str] = []
    if not isinstance(value, Mapping):
        return ["step lacks its failed-leaf controller evidence"]
    contract = floor_contract(requested_dt_s)
    finite_leaf_durations = tuple(
        float(item) for item in leaf_durations if _finite_number(item)
    )
    minimum_accepted_leaf = (
        min(finite_leaf_durations) if finite_leaf_durations else None
    )
    expected_scalars = {
        "passed": True,
        "policy_id": POLICY_ID,
        "licensed_only_after_imposed_boundary_discontinuity": True,
        "absolute_minimum_leaf_dt_s": ABSOLUTE_MINIMUM_LEAF_DT_S,
        "maximum_binary_depth_cap": MAXIMUM_BINARY_DEPTH_CAP,
        "effective_maximum_binary_depth": contract[
            "effective_maximum_binary_depth"
        ],
        "effective_terminal_leaf_dt_s": contract[
            "effective_terminal_leaf_dt_s"
        ],
        "accepted_binary_depth_bounded": True,
        "maximum_leaf_count": contract["maximum_leaf_count"],
        "leaf_grid_disclosed_and_bounded": True,
        "minimum_accepted_leaf_dt_s": minimum_accepted_leaf,
        "no_leaf_below_absolute_floor": True,
        "floor_was_binding": contract["floor_was_binding"],
        "depth_cap_was_binding": contract["depth_cap_was_binding"],
        "per_seed_and_per_leaf_work_budget_unchanged_from_amendment03": True,
        "total_tree_work_envelope_changed_by_amendment09": True,
        "predictor_backtrack_maximum_binary_exponent": (
            PREDICTOR_BACKTRACK_MAXIMUM_BINARY_EXPONENT
        ),
        "provisional_optimizer_function_evaluations_per_seed": (
            PROVISIONAL_FUNCTION_EVALUATIONS_PER_SEED
        ),
        "reversal_maximum_nonlinear_seed_attempts": (
            REVERSAL_MAXIMUM_SEED_ATTEMPTS
        ),
        "nonreversal_maximum_nonlinear_seed_attempts": (
            NONREVERSAL_MAXIMUM_SEED_ATTEMPTS
        ),
        "fixed_work_budget_observed": True,
        "attempt_paths_form_failed_leaf_binary_tree": True,
        "initial_requested_macro_failure_disclosed": True,
        "one_returned_attempt_per_accepted_leaf": True,
        "per_leaf_seed_attempt_budget_observed": True,
        "atomic_macro_input_identity_preserved": True,
        "immutable_leaf_chain_identity_preserved": True,
        "exact_post_jump_boundary_identity_preserved": True,
        "every_attempt_rollback_identity_preserved": True,
        "raw_macro_ledger_recomputed_from_original_and_final_inventory": adaptive,
        "raw_macro_ledger_le_1e_10": True,
        "arbitrary_subdivision_allowed": False,
    }
    for key, expected in expected_scalars.items():
        observed = value.get(key)
        if isinstance(expected, float) and math.isfinite(expected):
            matches = _close(observed, expected)
        else:
            matches = observed == expected
        if not matches:
            errors.append(f"step controller evidence changed {key}")
    return errors


def _attempt_suffix(role: Any, *, depth: int, path: str) -> str | None:
    if not isinstance(role, str):
        return None
    prefix = f"adaptive_d{depth}_{path}::"
    return role[len(prefix) :] if role.startswith(prefix) else None


def _adaptive_seed_candidate_sequence(
    *,
    is_reversal: bool,
) -> tuple[tuple[frozenset[str], float], ...]:
    """Return the exact ordered candidate lattice used by adaptive leaves.

    An exact open-chart refusal consumes no nonlinear-seed work, but it does
    consume its frozen seed position.  Amendment 03 licenses only the
    half-history predictor followed by the fixed anchor on a reversal leaf,
    and only full-history, anchor, then half-history on a non-reversal leaf.
    The broader binary predictor scan remains an ordinary-solve facility; it
    is not a license for a fourth provisional failed-leaf seed.
    """

    predictor_roles = frozenset(
        {_REQUESTED_ROLE, _PREDICTOR_CHART_VALIDATION_ROLE}
    )
    first_predictor = (predictor_roles, 0.5 if is_reversal else 1.0)
    anchor = (
        frozenset({_ANCHOR_ROLE, _ANCHOR_CHART_VALIDATION_ROLE}),
        1.0,
    )
    if is_reversal:
        return first_predictor, anchor
    return (
        first_predictor,
        anchor,
        (
            predictor_roles,
            0.5,
        ),
    )


def _adaptive_seed_sequence_errors(
    entries: Sequence[tuple[Mapping[str, Any], str]],
    *,
    is_reversal: bool,
    status: Any,
) -> list[str]:
    """Verify candidate order, zero-work chart refusals, and budget closure."""

    errors: list[str] = []
    maximum_attempts = (
        REVERSAL_MAXIMUM_SEED_ATTEMPTS
        if is_reversal
        else NONREVERSAL_MAXIMUM_SEED_ATTEMPTS
    )
    budget_markers = [
        index
        for index, (_, suffix) in enumerate(entries)
        if suffix == _SEED_BUDGET_EXHAUSTED_ROLE
    ]
    marker_is_exact = budget_markers == [len(entries) - 1]
    if budget_markers and not marker_is_exact:
        errors.append("adaptive nonlinear seed-budget marker is misplaced")
    candidate_entries = entries[:-1] if marker_is_exact else entries
    expected_candidates = _adaptive_seed_candidate_sequence(
        is_reversal=is_reversal
    )
    if len(candidate_entries) > len(expected_candidates):
        errors.append("adaptive ordered seed candidate sequence changed")
    for index, (attempt, suffix) in enumerate(candidate_entries):
        if index >= len(expected_candidates):
            break
        expected_roles, expected_fraction = expected_candidates[index]
        if suffix not in expected_roles or not _close(
            attempt.get("predictor_fraction"),
            expected_fraction,
        ):
            errors.append("adaptive ordered seed candidate sequence changed")
            break
        if suffix in {
            _PREDICTOR_CHART_VALIDATION_ROLE,
            _ANCHOR_CHART_VALIDATION_ROLE,
        } and (
            attempt.get("accepted") is not False
            or attempt.get("nonlinear_evaluations") != 0
            or attempt.get("provisional_work_budget") is not None
            or attempt.get("provisional_seed_attempt_budget_exhausted") is not False
        ):
            errors.append("adaptive chart refusal consumed nonlinear work or was accepted")

    nonlinear = [
        (attempt, suffix)
        for attempt, suffix in candidate_entries
        if suffix in _NONLINEAR_ROLES
    ]
    accepted_candidate_indices = [
        index
        for index, (attempt, _) in enumerate(candidate_entries)
        if attempt.get("accepted") is True
    ]
    if accepted_candidate_indices and accepted_candidate_indices != [
        len(candidate_entries) - 1
    ]:
        errors.append("adaptive seed candidate sequence continued after acceptance")
    if not 0 <= len(nonlinear) <= maximum_attempts:
        errors.append("adaptive per-leaf nonlinear seed budget changed")
    if any(
        attempt.get("provisional_work_budget")
        != PROVISIONAL_FUNCTION_EVALUATIONS_PER_SEED
        for attempt, _ in nonlinear
    ):
        errors.append("adaptive per-seed nonlinear work budget changed")

    marker = entries[-1][0] if marker_is_exact else None
    if marker is not None and (
        marker.get("accepted") is not False
        or not _close(marker.get("predictor_fraction"), 0.0)
        or marker.get("nonlinear_evaluations") != 0
        or marker.get("provisional_work_budget")
        != PROVISIONAL_FUNCTION_EVALUATIONS_PER_SEED
        or marker.get("provisional_seed_attempt_budget_exhausted") is not True
    ):
        errors.append("adaptive nonlinear seed-budget marker changed")
    if status == "failed_leaf_bisected":
        if len(candidate_entries) != len(expected_candidates):
            errors.append("bisected adaptive leaf did not exhaust its ordered seed slots")
        if marker_is_exact:
            if len(nonlinear) != maximum_attempts:
                errors.append(
                    "adaptive nonlinear seed-budget marker precedes budget exhaustion"
                )
        elif len(nonlinear) == maximum_attempts:
            errors.append("bisected adaptive leaf omitted its seed-budget exhaustion marker")
    if status == "returned":
        if not nonlinear:
            errors.append("returned adaptive leaf never evaluated a nonlinear seed")
        if budget_markers:
            errors.append("returned adaptive leaf carries a seed-budget exhaustion marker")
    return errors


def _partition_rows(
    value: Any,
    *,
    label: str,
) -> tuple[list[Mapping[str, Any]], list[Any], list[Any], list[str]]:
    """Return disclosed partition rows without assigning them policy meaning."""

    errors: list[str] = []
    if not isinstance(value, list):
        return [], [], [], [f"step does not disclose {label}"]
    rows: list[Mapping[str, Any]] = []
    durations: list[Any] = []
    nodes: list[Any] = []
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            errors.append(f"{label} partition {index} is not an object")
            continue
        rows.append(item)
        duration = item.get("duration_s")
        node = item.get("time_after_s")
        durations.append(duration)
        nodes.append(node)
        if not _finite_number(duration) or float(duration) <= 0.0:
            errors.append(f"{label} partition {index} has invalid duration")
        if not _finite_number(node):
            errors.append(f"{label} partition {index} has invalid endpoint")
    return rows, durations, nodes, errors


def _partition_telescope_errors(
    durations: Sequence[Any],
    nodes: Sequence[Any],
    *,
    label: str,
    requested_dt_s: float,
    expected_time_before_s: float,
    expected_time_after_s: float,
) -> list[str]:
    """Verify one authoritative time partition against the macro endpoints."""

    errors: list[str] = []
    if not durations:
        return [f"{label} partition is empty"]
    if len(durations) != len(nodes):
        return [f"{label} durations and endpoints have different lengths"]
    if any(
        not _finite_number(duration) or float(duration) <= 0.0
        for duration in durations
    ) or any(not _finite_number(node) for node in nodes):
        return errors
    if not _close(math.fsum(float(item) for item in durations), requested_dt_s):
        errors.append(f"{label} partitions do not sum to the macro width")
    running = expected_time_before_s
    for index, (duration, node) in enumerate(zip(durations, nodes, strict=True)):
        running = math.fsum((running, float(duration)))
        if not _close(node, running):
            errors.append(f"{label} partition {index} does not telescope")
    if not _close(nodes[-1], expected_time_after_s):
        errors.append(f"{label} partition misses the macro endpoint")
    return errors


def _legacy_partition_alias_errors(
    record: Mapping[str, Any],
    *,
    authoritative_alias: str,
    authoritative_durations: Sequence[Any],
    authoritative_nodes: Sequence[Any],
) -> list[str]:
    """Check the explicitly non-authoritative backward-compatibility alias."""

    errors: list[str] = []
    marker = record.get("legacy_accepted_internal_substep_reporting")
    expected_marker = {
        "authoritative": False,
        "alias_of": authoritative_alias,
        "semantic_role": "backward_compatibility_balance_interval_alias_only",
        "must_not_be_used_to_apply_amendment09_floor": True,
    }
    if not isinstance(marker, Mapping) or dict(marker) != expected_marker:
        errors.append("legacy accepted-substep reporting is not explicitly non-authoritative")

    count = record.get("accepted_internal_substep_count")
    durations = record.get("accepted_internal_substep_dt_s")
    nodes = record.get("accepted_internal_time_nodes_s")
    if not _integer(count, lower=0):
        errors.append("legacy accepted-substep count is invalid")
        count = -1
    if not isinstance(durations, list) or len(durations) != count:
        errors.append("legacy accepted-substep durations are not fully disclosed")
        durations = []
    if not isinstance(nodes, list) or len(nodes) != count:
        errors.append("legacy accepted-substep endpoints are not fully disclosed")
        nodes = []
    if count != len(authoritative_durations):
        errors.append("legacy accepted-substep count differs from its authoritative alias")
    if len(durations) != len(authoritative_durations) or any(
        not _close(observed, float(expected))
        for observed, expected in zip(durations, authoritative_durations)
        if _finite_number(expected)
    ):
        errors.append("legacy accepted-substep durations differ from their authoritative alias")
    if len(nodes) != len(authoritative_nodes) or any(
        not _close(observed, float(expected))
        for observed, expected in zip(nodes, authoritative_nodes)
        if _finite_number(expected)
    ):
        errors.append("legacy accepted-substep endpoints differ from their authoritative alias")
    return errors


def _physical_event_contract_errors(
    record: Mapping[str, Any],
    partitions: Sequence[Mapping[str, Any]],
    durations: Sequence[Any],
    nodes: Sequence[Any],
    attempts: Sequence[Any],
    *,
    requested_dt_s: float,
    expected_time_before_s: float,
    expected_time_after_s: float,
) -> list[str]:
    """Recompute the accepted face-event transaction without A09 leaf rules."""

    errors = _partition_telescope_errors(
        durations,
        nodes,
        label="physical event",
        requested_dt_s=requested_dt_s,
        expected_time_before_s=expected_time_before_s,
        expected_time_after_s=expected_time_after_s,
    )
    roles = tuple(partition.get("role") for partition in partitions)
    if roles != _EVENT_PARTITION_ROLES:
        errors.append("physical event partitions are not arrival then departure")
    if not _integer(record.get("interior_face_event_index"), lower=0):
        errors.append("physical event does not disclose its interior face index")

    value = record.get("physical_event_contract")
    if not isinstance(value, Mapping):
        return [*errors, "physical event contract is missing"]
    if value.get("event_kind") != "interior_master_face_arrival_and_departure":
        errors.append("physical event kind changed")
    if value.get("accepted_balance_interval_count") != 2:
        errors.append("physical event does not disclose exactly two balance intervals")
    if len(durations) == 2:
        if not _close(value.get("event_duration_s"), float(durations[0])):
            errors.append("physical event duration differs from its arrival partition")
        if not _close(value.get("remainder_duration_s"), float(durations[1])):
            errors.append("physical event remainder differs from its departure partition")
        if not _close(record.get("interior_face_event_duration_s"), float(durations[0])):
            errors.append("top-level physical event duration is inconsistent")
        if not _close(
            record.get("interior_face_departure_remainder_s"),
            float(durations[1]),
        ):
            errors.append("top-level physical event remainder is inconsistent")
    if not _close(value.get("duration_closure_error_s"), 0.0):
        errors.append("physical event duration closure is not exact")
    if not _close(value.get("absolute_time_closure_error_s"), 0.0):
        errors.append("physical event absolute-time closure is not exact")

    mode = value.get("ordering_mode")
    if mode not in _EVENT_ORDERING_MODES:
        errors.append("physical event ordering mode is unlicensed")
    if record.get("interior_face_event_ordering_mode") != mode:
        errors.append("top-level physical event ordering mode is inconsistent")
    expected_ordering = [
        "same_cell_rejected",
        mode,
        *_EVENT_PARTITION_ROLES,
    ]
    if value.get("ordering") != expected_ordering:
        errors.append("physical event ordering evidence changed")
    error_bound = value.get("event_time_error_bound_s")
    if mode == "rigorous_interval_certificate":
        if not _finite_number(error_bound) or float(error_bound) < 0.0:
            errors.append("rigorous physical event lacks a finite time-error bound")
    elif error_bound is not None:
        errors.append("open-interval physical event reports a false rigorous error bound")

    chain = value.get("immutable_chain")
    expected_chain = {field: True for field in _EVENT_CHAIN_FIELDS}
    if not isinstance(chain, Mapping) or dict(chain) != expected_chain:
        errors.append("physical event immutable object chain is not exact")
    if value.get("every_attempt_rollback_identity_preserved") is not True:
        errors.append("physical event contract does not preserve attempt rollback")
    if any(
        not isinstance(attempt, Mapping)
        or attempt.get("rollback_identity_preserved") is not True
        for attempt in attempts
    ):
        errors.append("physical event attempt rollback is not exact")

    ledger = value.get("raw_macro_component_energy_ledger")
    if not isinstance(ledger, Mapping):
        return [*errors, "physical event raw component/energy ledger is missing"]
    for field in _EVENT_RAW_LEDGER_FIELDS:
        if not _finite_number(ledger.get(field)):
            errors.append(f"physical event raw ledger does not disclose {field}")
    for field, top_field in zip(_EVENT_RAW_LEDGER_FIELDS[:3], _RAW_MACRO_LEDGER_FIELDS):
        observed = ledger.get(field)
        top = record.get(top_field)
        if _finite_number(observed) and (
            not _finite_number(top) or not _close(top, float(observed))
        ):
            errors.append(f"physical event raw ledger differs from {top_field}")

    normalized_values: list[float] = []
    for field in _NORMALIZED_LEDGER_FIELDS:
        observed = ledger.get(field)
        top = record.get(field)
        if not _finite_number(observed):
            errors.append(f"physical event normalized ledger does not disclose {field}")
            continue
        normalized_values.append(float(observed))
        if not _finite_number(top) or not _close(top, float(observed)):
            errors.append(f"physical event normalized ledger differs from {field}")
        if abs(float(observed)) > LEDGER_LIMIT:
            errors.append(f"physical event normalized ledger violates {field}")
    if len(normalized_values) == len(_NORMALIZED_LEDGER_FIELDS):
        maximum_step = max(abs(value) for value in normalized_values[:3])
        maximum_cumulative = max(abs(value) for value in normalized_values[3:])
        if not _close(ledger.get("maximum_normalized_step_residual"), maximum_step):
            errors.append("physical event maximum step ledger residual is misreported")
        if not _close(
            ledger.get("maximum_normalized_cumulative_residual"),
            maximum_cumulative,
        ):
            errors.append("physical event maximum cumulative ledger residual is misreported")
        if not _close(record.get("interior_face_raw_macro_ledger_residual"), maximum_step):
            errors.append("top-level physical event ledger residual is inconsistent")
    if not _close(ledger.get("limit"), LEDGER_LIMIT):
        errors.append("physical event component/energy ledger limit changed")
    return errors


def record_contract_errors(
    record: Any,
    *,
    requested_dt_s: float,
    expected_time_before_s: float,
    expected_time_after_s: float,
    starts_after_boundary_discontinuity: bool,
) -> list[str]:
    """Recompute ordinary, adaptive, or physical-event macro acceptance."""

    errors: list[str] = []
    if not isinstance(record, Mapping):
        return ["accepted step record is not an object"]
    adaptive_flag = record.get("adaptive_subdivision_used")
    adaptive = adaptive_flag is True
    event_flag = record.get("interior_face_event_used")
    event_used = event_flag is True
    if adaptive_flag not in (True, False):
        errors.append("step adaptive flag is not boolean")
    if event_flag not in (True, False):
        errors.append("step physical-event flag is not boolean")
    if adaptive and event_used:
        errors.append("step cannot be both an A09 failed-leaf recovery and a physical event")
    if not _close(record.get("dt_s"), requested_dt_s):
        errors.append("accepted macro width differs from its requested grid")
    if not _close(record.get("time_before_s"), expected_time_before_s):
        errors.append("accepted macro start is off the requested node")
    if not _close(record.get("time_after_s"), expected_time_after_s):
        errors.append("accepted macro end is off the requested node")
    if (
        record.get("starts_after_imposed_boundary_discontinuity")
        is not starts_after_boundary_discontinuity
    ):
        errors.append("step discontinuity position is misreported")

    failed_rows, leaf_durations, leaf_nodes, partition_errors = _partition_rows(
        record.get("failed_leaf_partitions"),
        label="failed-leaf",
    )
    errors.extend(partition_errors)
    event_rows, event_durations, event_nodes, partition_errors = _partition_rows(
        record.get("physical_event_partitions"),
        label="physical-event",
    )
    errors.extend(partition_errors)

    authoritative_alias = (
        "physical_event_partitions" if event_used else "failed_leaf_partitions"
    )
    authoritative_durations = event_durations if event_used else leaf_durations
    authoritative_nodes = event_nodes if event_used else leaf_nodes
    errors.extend(
        _legacy_partition_alias_errors(
            record,
            authoritative_alias=authoritative_alias,
            authoritative_durations=authoritative_durations,
            authoritative_nodes=authoritative_nodes,
        )
    )

    floor = floor_contract(requested_dt_s)
    expected_floor_fields = {
        "adaptive_absolute_minimum_leaf_dt_s": ABSOLUTE_MINIMUM_LEAF_DT_S,
        "adaptive_effective_maximum_binary_depth": floor[
            "effective_maximum_binary_depth"
        ],
        "adaptive_effective_terminal_leaf_dt_s": floor[
            "effective_terminal_leaf_dt_s"
        ],
        "adaptive_floor_was_binding": floor["floor_was_binding"],
        "adaptive_depth_cap_was_binding": floor["depth_cap_was_binding"],
        "adaptive_no_leaf_below_absolute_floor": True,
    }
    for field, expected in expected_floor_fields.items():
        observed = record.get(field)
        matches = _close(observed, expected) if isinstance(expected, float) else observed == expected
        if not matches:
            errors.append(f"step absolute-floor evidence changed {field}")
    finite_leaf_durations = tuple(float(item) for item in leaf_durations if _finite_number(item))
    minimum_accepted_leaf = (
        min(finite_leaf_durations) if finite_leaf_durations else None
    )
    observed_minimum = record.get("adaptive_minimum_accepted_leaf_dt_s")
    minimum_matches = (
        observed_minimum is None
        if minimum_accepted_leaf is None
        else _close(observed_minimum, minimum_accepted_leaf)
    )
    if not minimum_matches:
        errors.append("step minimum accepted leaf is misreported")

    errors.extend(
        _controller_evidence_errors(
            record.get("frozen_failed_leaf_controller_contract"),
            adaptive=adaptive,
            requested_dt_s=requested_dt_s,
            leaf_durations=leaf_durations,
        )
    )
    attempts = record.get("attempts")
    if not isinstance(attempts, list) or not attempts:
        errors.append("step does not disclose its solver attempts")
        attempts = []
    if any(
        not isinstance(attempt, Mapping)
        or attempt.get("rollback_identity_preserved") is not True
        for attempt in attempts
    ):
        errors.append("a disclosed solver attempt lacks exact rollback identity")

    if event_used:
        if failed_rows:
            errors.append("physical event is incorrectly reported as failed-leaf subdivision")
        if adaptive:
            errors.append("physical event incorrectly invokes the failed-leaf controller")
        if record.get("recovery_mode", "").startswith("adaptive_"):
            errors.append("physical event carries an adaptive recovery label")
        if record.get("adaptive_binary_depth") != 0:
            errors.append("physical event has nonzero adaptive depth")
        if record.get("adaptive_attempted_binary_depths") != []:
            errors.append("physical event reports adaptive depth attempts")
        if record.get("adaptive_predictor_backtrack_maximum_binary_exponent") is not None:
            errors.append("physical event reports an adaptive backtrack budget")
        errors.extend(
            _physical_event_contract_errors(
                record,
                event_rows,
                event_durations,
                event_nodes,
                attempts,
                requested_dt_s=requested_dt_s,
                expected_time_before_s=expected_time_before_s,
                expected_time_after_s=expected_time_after_s,
            )
        )
        return errors

    if event_rows:
        errors.append("non-event macrostep discloses physical event partitions")
    if record.get("physical_event_contract") is not None:
        errors.append("non-event macrostep discloses a physical event contract")
    errors.extend(
        _partition_telescope_errors(
            leaf_durations,
            leaf_nodes,
            label="failed-leaf",
            requested_dt_s=requested_dt_s,
            expected_time_before_s=expected_time_before_s,
            expected_time_after_s=expected_time_after_s,
        )
    )
    count = len(failed_rows)

    if not adaptive:
        if record.get("recovery_mode", "").startswith("adaptive_"):
            errors.append("ordinary macrostep carries an adaptive recovery label")
        if count != 1 or len(leaf_durations) != 1 or not _close(
            leaf_durations[0] if leaf_durations else None,
            requested_dt_s,
        ):
            errors.append("ordinary macrostep is not one requested-width BE leaf")
        elif (
            failed_rows[0].get("role") != "direct_requested_backward_euler"
            or failed_rows[0].get("binary_depth") != 0
            or failed_rows[0].get("binary_path") is not None
        ):
            errors.append("ordinary macrostep failed-leaf partition metadata changed")
        if record.get("adaptive_binary_depth") != 0:
            errors.append("ordinary macrostep has nonzero adaptive depth")
        if record.get("adaptive_attempted_binary_depths") != []:
            errors.append("ordinary macrostep reports adaptive depth attempts")
        if record.get("adaptive_predictor_backtrack_maximum_binary_exponent") is not None:
            errors.append("ordinary macrostep reports an adaptive backtrack budget")
        if any(
            isinstance(attempt, Mapping)
            and (
                attempt.get("adaptive_depth") is not None
                or attempt.get("adaptive_path") is not None
                or (
                    attempt.get("macro_branch_status") != "not_applicable"
                    # GT-PS-2-P1E-R12: the strict same-cell chart's documented
                    # FIRST REFUSAL leaves its mandatory audit record on the
                    # accepted ordinary step (A20 replay depends on it).  A
                    # refusal carries no adaptive depth and no adaptive path;
                    # it is not adaptive execution and does not breach the
                    # ordinary-macrostep cleanliness this sweep protects.
                    and not (
                        attempt.get("macro_branch_status") == "same_cell"
                        and attempt.get("adaptive_depth") is None
                        and attempt.get("adaptive_path") is None
                        and attempt.get("role")
                        == "accepted_history_bracket_same_cell_first_refusal"
                    )
                )
            )
            for attempt in attempts
        ):
            errors.append("ordinary macrostep contains adaptive attempt metadata")
        return errors

    if not starts_after_boundary_discontinuity:
        errors.append("adaptive recovery occurred away from a boundary discontinuity")
    effective_depth = int(floor["effective_maximum_binary_depth"])
    if effective_depth < 1:
        errors.append("adaptive recovery was used when the absolute floor forbids a child")
        return errors
    depth = record.get("adaptive_binary_depth")
    if not _integer(depth, lower=1) or depth > effective_depth:
        errors.append(
            "adaptive binary depth exceeds the requested grid's absolute-floor limit"
        )
        return errors
    if record.get("recovery_mode") != f"adaptive_same_cell_failed_leaf_depth_{depth}":
        errors.append("adaptive recovery mode does not disclose its exact depth")
    if record.get("adaptive_attempted_binary_depths") != list(range(1, depth + 1)):
        errors.append("adaptive attempted depths are not the exact consecutive sequence")
    if (
        record.get("adaptive_predictor_backtrack_maximum_binary_exponent")
        != PREDICTOR_BACKTRACK_MAXIMUM_BINARY_EXPONENT
    ):
        errors.append("adaptive predictor backtrack budget changed")
    if not 2 <= count <= 2**effective_depth:
        errors.append("adaptive accepted leaf count exceeds the effective tree bound")

    for index, leaf_dt in enumerate(leaf_durations):
        if not _finite_number(leaf_dt):
            continue
        matches_binary_grid = any(
            _close(leaf_dt, math.ldexp(requested_dt_s, -candidate_depth))
            for candidate_depth in range(1, depth + 1)
        )
        if not matches_binary_grid or float(leaf_dt) < ABSOLUTE_MINIMUM_LEAF_DT_S:
            errors.append(f"adaptive leaf {index} is outside the licensed binary grid")

    for field in _NORMALIZED_LEDGER_FIELDS:
        value = record.get(field)
        if not _finite_number(value) or abs(float(value)) > LEDGER_LIMIT:
            errors.append(f"adaptive raw macro ledger violates {field}")
    for field in _RAW_MACRO_LEDGER_FIELDS:
        if not _finite_number(record.get(field)):
            errors.append(f"adaptive raw macro ledger does not disclose {field}")

    root_attempts: list[Mapping[str, Any]] = []
    grouped: dict[tuple[int, str], list[tuple[Mapping[str, Any], str]]] = {}
    for attempt in attempts:
        if not isinstance(attempt, Mapping):
            continue
        attempt_depth = attempt.get("adaptive_depth")
        path = attempt.get("adaptive_path")
        if attempt_depth == 0 and path == "root":
            root_attempts.append(attempt)
            continue
        if (
            not _integer(attempt_depth, lower=1)
            or attempt_depth > depth
            or not isinstance(path, str)
            or len(path) != attempt_depth
            or not path
            or set(path) - {"L", "R"}
        ):
            errors.append("adaptive attempt has an invalid depth/path")
            continue
        suffix = _attempt_suffix(attempt.get("role"), depth=attempt_depth, path=path)
        if suffix not in _NONLINEAR_ROLES | _NONLINEAR_EXEMPT_ROLES:
            errors.append("adaptive attempt has an unlicensed seed role")
            continue
        if not _close(attempt.get("attempted_dt_s"), math.ldexp(requested_dt_s, -attempt_depth)):
            errors.append("adaptive attempt width differs from its binary path")
        if attempt.get("macro_branch_status") not in {
            "failed_leaf_bisected",
            "returned",
        }:
            errors.append("adaptive attempt has an unlicensed branch status")
        grouped.setdefault((attempt_depth, path), []).append((attempt, suffix))

    if not root_attempts:
        errors.append("adaptive record lacks the rejected requested macro attempt")
    for attempt in root_attempts:
        if (
            attempt.get("accepted") is not False
            or attempt.get("macro_branch_status") != "initial_requested_failure"
            or not _close(attempt.get("attempted_dt_s"), requested_dt_s)
        ):
            errors.append("adaptive root attempt is not the exact rejected macro solve")

    paths = {path for _, path in grouped}
    if not {"L", "R"} <= paths:
        errors.append("adaptive binary tree does not start with both requested halves")
    returned_paths: list[str] = []
    failed_paths: set[str] = set()
    for (attempt_depth, path), entries in grouped.items():
        statuses = {attempt.get("macro_branch_status") for attempt, _ in entries}
        if len(statuses) != 1:
            errors.append("one adaptive leaf mixes returned and bisected statuses")
            continue
        status = next(iter(statuses))
        is_reversal = path == "L" * attempt_depth
        errors.extend(
            _adaptive_seed_sequence_errors(
                entries,
                is_reversal=is_reversal,
                status=status,
            )
        )
        nonlinear = [
            (attempt, suffix)
            for attempt, suffix in entries
            if suffix in _NONLINEAR_ROLES
        ]
        accepted = [attempt for attempt, _ in nonlinear if attempt.get("accepted") is True]
        if status == "returned":
            returned_paths.append(path)
            if len(accepted) != 1 or nonlinear[-1][0].get("accepted") is not True:
                errors.append("returned adaptive leaf lacks one final accepted root")
        else:
            failed_paths.add(path)
            if accepted or attempt_depth >= effective_depth:
                errors.append("failed adaptive leaf was accepted or bisected too deeply")

    for path in failed_paths:
        if f"{path}L" not in paths or f"{path}R" not in paths:
            errors.append("adaptive controller did not bisect only the failed leaf")
    for path in returned_paths:
        if any(other.startswith(path) and other != path for other in paths):
            errors.append("adaptive controller subdivided an already returned sibling")

    returned_intervals = sorted(
        (*_path_interval(path), path) for path in returned_paths
    )
    cursor = 0.0
    ordered_paths: list[str] = []
    for start, end, path in returned_intervals:
        if not math.isclose(start, cursor, rel_tol=0.0, abs_tol=TIME_ABS_TOLERANCE_S):
            errors.append("returned adaptive leaves do not tile the macro interval")
            break
        cursor = end
        ordered_paths.append(path)
    if not math.isclose(cursor, 1.0, rel_tol=0.0, abs_tol=TIME_ABS_TOLERANCE_S):
        errors.append("returned adaptive leaves do not cover the macro endpoint")
    expected_leaf_durations = [requested_dt_s / (2 ** len(path)) for path in ordered_paths]
    if len(expected_leaf_durations) != len(leaf_durations) or any(
        not _close(observed, expected)
        for observed, expected in zip(leaf_durations, expected_leaf_durations)
    ):
        errors.append("disclosed accepted leaves differ from returned attempt paths")
    if len(failed_rows) != len(ordered_paths) or any(
        partition.get("role") != "amendment09_failed_leaf_binary_partition"
        or partition.get("binary_depth") != len(path)
        or partition.get("binary_path") != path
        for partition, path in zip(failed_rows, ordered_paths)
    ):
        errors.append("failed-leaf partition metadata differs from returned attempt paths")
    if sorted({len(path) for path in paths}) != list(range(1, depth + 1)):
        errors.append("adaptive attempt tree differs from disclosed attempted depths")
    return errors


def campaign_contract_errors(
    raw: Any,
    *,
    expected_timesteps_s: Sequence[float],
    segment_duration_s: float,
) -> list[str]:
    """Verify every job and record in one raw dynamic campaign."""

    if not isinstance(raw, Mapping):
        return ["raw backend payload is missing"]
    errors = solver_contract_errors(
        raw.get("solver_acceptance_contract"),
        expected_timesteps_s=expected_timesteps_s,
    )
    jobs = raw.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        return [*errors, "raw backend contains no jobs"]
    for job_index, job in enumerate(jobs):
        label = f"job[{job_index}]"
        if not isinstance(job, Mapping):
            errors.append(f"{label} is not an object")
            continue
        requested_dt = job.get("dt_s")
        if requested_dt not in expected_timesteps_s:
            errors.append(f"{label} uses an unlicensed requested macro width")
            continue
        if not _close(job.get("segment_duration_s"), segment_duration_s):
            errors.append(f"{label} segment duration differs from the formal history")
        requested_steps = round(segment_duration_s / float(requested_dt))
        if not math.isclose(
            requested_steps * float(requested_dt),
            segment_duration_s,
            rel_tol=0.0,
            abs_tol=TIME_ABS_TOLERANCE_S,
        ):
            errors.append(f"{label} requested grid does not align the segment")
            continue
        records = job.get("step_records")
        expected_count = 3 * requested_steps
        if not isinstance(records, list) or len(records) != expected_count:
            errors.append(f"{label} does not disclose every requested macro interval")
            continue
        for index, record in enumerate(records):
            time_before = index * float(requested_dt)
            time_after = (index + 1) * float(requested_dt)
            starts_after_jump = index in (requested_steps, 2 * requested_steps)
            for error in record_contract_errors(
                record,
                requested_dt_s=float(requested_dt),
                expected_time_before_s=time_before,
                expected_time_after_s=time_after,
                starts_after_boundary_discontinuity=starts_after_jump,
            ):
                errors.append(f"{label} record[{index}]: {error}")
        gates = job.get("gates")
        if not isinstance(gates, Mapping):
            errors.append(f"{label} lacks numerical gates")
            continue
        if "no_step_subdivision" in gates:
            errors.append(f"{label} retains the superseded no-subdivision gate")
        if gates.get("exact_requested_time_grid") is not True:
            errors.append(f"{label} failed the exact requested macro-node gate")
        if gates.get("frozen_failed_leaf_controller_contract_passed") is not True:
            errors.append(f"{label} failed the frozen failed-leaf controller gate")
        scenario = job.get("scenario")
        expected_jh_gate = scenario in {"composition_only", "combined"}
        if gates.get("surface_hexane_component_N_h_direction_is_gate") is not False:
            errors.append(f"{label} incorrectly makes conserved N_h a direction gate")
        if gates.get("surface_total_stefan_flux_is_direction_gate") is not False:
            errors.append(f"{label} incorrectly makes total Stefan flux a direction gate")
        if gates.get("surface_independent_hexane_J_h_direction_is_gate") is not expected_jh_gate:
            errors.append(f"{label} has the wrong independent J_h direction authority")
        instantaneous = job.get("instantaneous_flux_sign_diagnostics")
        if not isinstance(instantaneous, Mapping) or instantaneous.get(
            "excluded_from_mesh_time_acceptance"
        ) is not True:
            errors.append(f"{label} gates instantaneous discontinuity flux signs")
    return errors
