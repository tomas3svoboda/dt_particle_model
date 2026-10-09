"""Qualify the Gate-1g same-cell solver under separated boundary pulses.

This campaign is numerical evidence, not soybean-DT parameter identification.
The transport coefficients and outer Dirichlet trace are deliberately labelled
nonqualifying.  The script changes neither a residual equation nor an accepted
state: mesh/time continuation is used only to construct nonlinear seeds.

The dry primitive chart requests 330--346 K at 1 atm.  The exact coupled
gas-only phase calculation raises its lower endpoint to about 334.4804 K.  The
wider request is intentional: 339.15 K was an artificial campaign floor that
became an active numerical bound as a refined dry-cell centre approached the
roughly 337 K RH interface.  Phase feasibility, rather than a rectangular
temperature/composition box, remains the actual lower authority.

Every reported trajectory advances to the same 0.3 s endpoint: a hot/lean
baseline for 0.1 s, one finite 0.1 s square pulse, and an exact return to the
same baseline object for 0.1 s.  The matrix separates composition-only,
temperature-only, and combined forcing so a response is never credited to a
mechanism that the experiment did not isolate.  Both discontinuities are
aligned exactly with every requested time grid.  Failed solves never enter a
trajectory.  Optional boundary-homotopy roots share the same old state and
are only nonlinear seeds for the final, unchanged backward-Euler problem.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field as dataclass_field, is_dataclass, replace
from enum import Enum
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import struct
from time import perf_counter
from typing import Iterable, Mapping, Sequence

import numpy as np

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import actual_composition_force_path as acfp
from dtdc_simulator.core2.particle import cut_continuation as cc
from dtdc_simulator.core2.particle import cut_event_orchestrator as event_orchestrator
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import front
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.particle import wet_retained_cap as wrc
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.particle.grid import uniform_grid
from dtdc_simulator.core2.props import sorption as sp
from scripts import p1ef_material_volume_profiles as material_profiles
from scripts import p1ef_surface_constitutive_oracle as surface_oracle


RADIUS_M = 0.885e-3
MAXIMUM_RADIUS_SEED_CONTINUATION_INCREMENT_M = 0.025e-3
UPPER_RADIUS_SEED_ANCHOR_M = 1.25e-3
CROSS_RADIUS_SEED_MAXIMUM_ATTEMPTS_PER_STEP = 1
# A rejected requested macrostep may be retried by recursively bisecting only
# the leaf that failed.  A successful sibling remains provisional and is not
# recomputed.  The accepted boundary is the unchanged post-jump object on every
# leaf, and the complete macro interval commits atomically only after its raw
# endpoint/flux-sum ledger passes.  Amendment-09 freezes one absolute lower
# leaf bound across the formal time grids.  The effective maximum depth is
# derived from the requested macro width and capped at eight; it is a numerical
# work allocation, not a physical time scale or a relaxed acceptance criterion.
ADAPTIVE_FAILED_LEAF_POLICY_ID = "frozen_amendment09_absolute_floor_failed_leaf_only"
ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S = 0.00029296875
ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH = 8
# Adaptive leaf solves use a bounded seed-globalization budget before the
# controller abandons that provisional branch and tries the next exact binary
# partition.  Ordinary fixed-step calls retain the exhaustive floating-point
# budget below.  This changes work only, never equations, bounds, tolerances,
# the requested macro interval, or an accepted state.
ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT = 8
# A discontinuity macro and every provisional child keep the normal nonlinear
# equations and acceptance contract, but a single doomed coarse solve may not
# monopolize the campaign.  This cap is passed only to provisional jump-branch
# attempts and is below the immutable state's general 800-evaluation ceiling.
# A converged candidate still has to pass every ordinary residual, ledger,
# conditioning, open-chart, and reduced-film guard before it can be retained.
DISCONTINUITY_LEAF_MAXIMUM_OPTIMIZER_FUNCTION_EVALUATIONS = 16
# Preserve the approved deterministic seed order: history-continuous tangent
# first, current-state/front-only anchor second.  If both bounded attempts
# reject, the leaf is bisected rather than trying dozens of same-dt seeds.
DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS = 2
# After a discontinuity has already been crossed by a certified leaf, a failed
# non-reversal sibling gets one final established continuation seed: the
# history predictor contracted to fraction 0.5.  It is ordered strictly after
# history-1.0 and the current-state anchor.  No fourth seed is licensed.
DISCONTINUITY_NONREVERSAL_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS = 3
PRESSURE_PA = 101_325.0
INITIAL_FRONT_RADIUS_FRACTION = 0.60
INITIAL_FRONT_Z = INITIAL_FRONT_RADIUS_FRACTION**3
SEGMENT_DURATION_S = 0.1
DEFAULT_MESHES = (12, 24, 48, 96)
DEFAULT_TIMESTEPS_S = (0.1, 0.05, 0.025)
BOOTSTRAP_MESHES = (2, 3, 4, 6, 7, 8, *DEFAULT_MESHES)
CONSERVATION_LIMIT = 1.0e-10
# GT-PS-2-P1E-21 section 3.2.  The nonlinear convergence gate has exactly one
# definition, the ``CutSolverControls`` authority; this module derives it and
# never restates the literal.  The former single ``RESIDUAL_LIMIT`` constant
# was overloaded: it also stood for the Amendment-18 saturated-interface
# log-fugacity tolerance, which is a different quantity and stays 2.0e-11
# (A21 section 3.3).  Reusing one name for both would have silently relaxed
# the A18 root tolerance and broken its unchanged-tolerance equality check.
NONLINEAR_RESIDUAL_LIMIT = ci.CutSolverControls.nonlinear_residual_tolerance
INTERFACE_ROOT_LOG_FUGACITY_LIMIT = 2.0e-11
STEFAN_SCALE_FLOOR = 1.0e-4
TRAJECTORY_SEGMENT_COUNT = 3
AUDIT_RECORD_ID = "G1G-COEFFICIENT-EVIDENCE-AUDIT-2026-08-01"
AUDIT_RECORD_PATH = "docs/GATE1G_COEFFICIENT_EVIDENCE_AUDIT_2026-08-01.md"
RETAINED_WATER_DIRECT_TEMPERATURE_MAX_K = 343.15
# Seed-only predictor backtracking follows exact binary fractions down to the
# floating-point resolution of a physical increment.  It never changes an
# accepted dt or boundary trace.  A candidate that rounds back onto the old
# open front bound is rejected by the existing chart validator.
PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT = 53
SCENARIO_COMPOSITION_ONLY = "composition_only"
SCENARIO_TEMPERATURE_ONLY = "temperature_only"
SCENARIO_COMBINED = "combined"
TEMPERATURE_ONLY_THERMAL_ORACLE_AMENDMENT_ID = "GT-PS-2-P1E-14"
FACE_CONDITIONED_PREDICTOR_AMENDMENT_ID = "GT-PS-2-P1E-17"
INTERFACE_ROOT_AMENDMENT_ID = "GT-PS-2-P1E-18"
INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION = 1
INTERFACE_ROOT_CERTIFICATION_SCHEMA_VERSION = 1
DRY_ENERGY_CAPACITY_AMENDMENT_ID = "GT-PS-2-P1E-19"
DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION = 1
ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID = "GT-PS-2-P1E-20"
# Version 2 (2026-08-11, owner batch ruling item 1, the R5S6-deferred "A20
# schema fix"): the committed face-0 total Stefan flux joins raw_replay_inputs,
# so a history-dependent successor operator's Pe input is part of the
# replayable schema. Version-1 artifacts are era-bound: valid against the
# qualifier that produced them, refused by this one (equality check), never
# retro-edited.
# Version 3 (2026-09-01, GT_PS2_A20_REPLAY_SCHEMA_V3_AMENDMENT_2026-09-01.md):
# the closed-hexane derivative authority joins raw_replay_inputs, recorded
# from the accepted production path's own audit. Commit 7a5d3bf (PART-01
# native derivative certificate) switched the production call site to
# NATIVE_LOCAL_ANALYTIC_ENVELOPE; a replay that omits the authority takes
# the GLOBAL default and the exact-equality gate correctly refuses. Version-2
# artifacts are era-bound under the same discipline as v1.
ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION = 3
TEMPERATURE_ONLY_THERMAL_DRIVING_FORCE_DEFINITION = (
    "delta_temperature_k=boundary_temperature_k-surface_film_audit.surface_temperature_k"
)
SCHEMA_VERSION = 5
FIXED_RADIUS_FORCING_HISTORY_ID = "gate1g_three_segment_hot_lean_cool_rich_hot_lean_v1"
ACCEPTED_TEMPERATURE_EXTREMA_SCOPE = (
    "exact_initial_and_accepted_BE_time_nodes_all_wet_dry_interface_pieces"
)


def _effective_adaptive_maximum_binary_depth(
    requested_dt_s: float,
    *,
    absolute_minimum_leaf_dt_s: float = (ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S),
    maximum_binary_depth_cap: int = ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH,
) -> int:
    """Return the deepest licensed binary leaf without crossing the floor.

    Exact ``ldexp`` comparisons make the policy deterministic and avoid a
    tolerance-adjusted logarithm.  Zero means that even the first bisection
    would cross the absolute floor; adaptive recovery must then fail closed
    before solving a child.
    """

    if not math.isfinite(requested_dt_s) or requested_dt_s <= 0.0:
        raise ValueError("adaptive requested dt must be positive and finite")
    if not math.isfinite(absolute_minimum_leaf_dt_s) or absolute_minimum_leaf_dt_s <= 0.0:
        raise ValueError("adaptive absolute leaf floor must be positive and finite")
    if (
        isinstance(maximum_binary_depth_cap, bool)
        or not isinstance(maximum_binary_depth_cap, int)
        or not 1 <= maximum_binary_depth_cap <= ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH
    ):
        raise ValueError(
            "adaptive binary-depth cap must be an integer in "
            f"[1, {ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH}]"
        )
    eligible = tuple(
        depth
        for depth in range(1, maximum_binary_depth_cap + 1)
        if math.ldexp(requested_dt_s, -depth) >= absolute_minimum_leaf_dt_s
    )
    return max(eligible, default=0)


def _adaptive_floor_contract(requested_dt_s: float) -> dict:
    """Report the immutable Amendment-09 work envelope for one macro width."""

    depth = _effective_adaptive_maximum_binary_depth(requested_dt_s)
    terminal = math.ldexp(requested_dt_s, -depth) if depth else requested_dt_s
    next_leaf = math.ldexp(requested_dt_s, -(depth + 1))
    return {
        "requested_dt_s": requested_dt_s,
        "effective_maximum_binary_depth": depth,
        "effective_terminal_leaf_dt_s": terminal,
        "maximum_leaf_count": 2**depth if depth else 1,
        "adaptive_recovery_available": depth >= 1,
        "no_leaf_below_absolute_floor": (
            depth == 0 or terminal >= ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S
        ),
        "floor_was_binding": (
            depth >= 1 and next_leaf < ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S
        ),
        "depth_cap_was_binding": (
            depth == ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH
            and next_leaf >= ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S
        ),
    }


def _face_conditioned_predictor_contract() -> dict[str, object]:
    """Report the exact Amendment-17 history-predictor policy."""

    return {
        "numerical_amendment_id": FACE_CONDITIONED_PREDICTOR_AMENDMENT_ID,
        "gas_only_predictor_coordinate": (
            cc.FACE_CONDITIONED_DRY_COMPOSITION_PREDICTOR
        ),
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


def _interface_root_contract() -> dict[str, object]:
    """Report the exact Amendment-18 interface-root policy."""

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


def _dry_energy_capacity_stencil_contract() -> dict[str, object]:
    """Report the exact Amendment-19 normalization-only policy."""

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


RETAINED_MOBILITY = cp.MobilityInterval(
    2.0e-9,
    2.0e-9,
    "Gate-1g numerical stress point; not soybean-DT identification",
).select(0.0)
DRY_CONFIG = ct.FullyDryTransportConfig(
    pressure_pa=PRESSURE_PA,
    binary_diffusivity=tc.EffectiveBinaryDiffusivityInterval(
        2.0e-9,
        2.0e-9,
        "Gate-1g numerical stress point; not soybean-DT identification",
    ).select(0.0),
    retained_water_mobility=RETAINED_MOBILITY,
    thermal_conductivity=ct.ThermalConductivitySelection(
        0.12,
        "Gate-1g numerical stress point; not soybean-DT identification",
    ),
    primitive_band=ct.GasOnlyPrimitiveDomain(
        (330.0, 346.0),
        (
            "nonqualifying stress request; exact pressure-conditioned gas-only "
            "phase domain remains authoritative"
        ),
    ),
)
WET_MODEL = ww.WetWaterModel(RETAINED_MOBILITY)
CONFIG = cut.CutTransportConfig(
    DRY_CONFIG,
    WET_MODEL,
    cut.InterfaceCompositionBracket(
        (0.50, 0.999),
        "interface-only saturated-root bracket; not a dry-cell band",
    ),
    moving_interface_composition_force_authority=(
        cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
    ),
)
CONTROLS = ci.CutSolverControls(
    wet_temperature_bounds_k=(330.0, 345.0),
    wet_water_bounds=(0.0231, 0.23),
    interface_temperature_bounds_k=(334.5, 340.0),
    nonlinear_residual_tolerance=NONLINEAR_RESIDUAL_LIMIT,
    ledger_tolerance=CONSERVATION_LIMIT,
    nonlinear_step_tolerance=1.0e-12,
    maximum_function_evaluations=800,
    maximum_condition_proxy=1.0e12,
)
HOT_LEAN_BOUNDARY = ct.DirichletPoreBoundary(
    345.9,
    0.7901,
    "nonqualifying hot/lean dynamic-step oracle",
)
REVERSAL_TEMPERATURE_K = 342.0
REVERSAL_EXACT_CHART_COORDINATE = 8.0
COLDEST_REVERSAL_FACE_TEMPERATURE_K = 0.5 * (
    DRY_CONFIG.conditioned_temperature_domain.solver_bounds_k[0] + REVERSAL_TEMPERATURE_K
)
REVERSAL_Y_HEXANE = cp.decode_gas_only_y(
    COLDEST_REVERSAL_FACE_TEMPERATURE_K,
    PRESSURE_PA,
    REVERSAL_EXACT_CHART_COORDINATE,
    DRY_CONFIG.pore,
)
# Composition-only pulse: temperature is exactly fixed at the baseline value.
# Coordinate q=8 is decoded at the coldest possible symmetric boundary face,
# the conservative endpoint for the pressure-conditioned gas-only chart.
COLDEST_COMPOSITION_FACE_TEMPERATURE_K = 0.5 * (
    DRY_CONFIG.conditioned_temperature_domain.solver_bounds_k[0] + HOT_LEAN_BOUNDARY.temperature_k
)
COMPOSITION_REVERSAL_Y_HEXANE = cp.decode_gas_only_y(
    COLDEST_COMPOSITION_FACE_TEMPERATURE_K,
    PRESSURE_PA,
    REVERSAL_EXACT_CHART_COORDINATE,
    DRY_CONFIG.pore,
)
HOT_RICH_BOUNDARY = ct.DirichletPoreBoundary(
    HOT_LEAN_BOUNDARY.temperature_k,
    COMPOSITION_REVERSAL_Y_HEXANE,
    (
        "nonqualifying composition-only q=8 oracle; y decoded at the coldest "
        "possible symmetric surface-face temperature"
    ),
)
COOL_LEAN_BOUNDARY = ct.DirichletPoreBoundary(
    REVERSAL_TEMPERATURE_K,
    HOT_LEAN_BOUNDARY.y_hexane,
    "nonqualifying temperature-only dynamic-step oracle",
)
# Correlated finite-pulse numerical stress trace.  The 3.9 K temperature drop
# and exact gas-only chart coordinate 8 give a phase-interior rich boundary at
# the coldest possible symmetric face.  It is a nonqualifying numerical stress
# input, not an asserted industrial feasibility-envelope corner.
COOL_RICH_BOUNDARY = ct.DirichletPoreBoundary(
    REVERSAL_TEMPERATURE_K,
    REVERSAL_Y_HEXANE,
    (
        "nonqualifying correlated reversal oracle; y decoded at the coldest "
        "possible symmetric surface-face temperature"
    ),
)


@dataclass(frozen=True)
class StressScenario:
    """One mechanism-labelled three-segment square-pulse experiment."""

    identifier: str
    forcing_class: str
    pulse_stage: str
    pulse_boundary: ct.DirichletPoreBoundary
    gate_authority: str

    @property
    def segments(self) -> tuple[tuple[str, ct.DirichletPoreBoundary], ...]:
        return (
            ("hot_lean_baseline_pre_pulse", HOT_LEAN_BOUNDARY),
            (self.pulse_stage, self.pulse_boundary),
            # Exact object identity makes the third segment a return, not an
            # independently rounded third boundary condition.
            ("hot_lean_baseline_return", HOT_LEAN_BOUNDARY),
        )


SCENARIOS = {
    SCENARIO_COMPOSITION_ONLY: StressScenario(
        SCENARIO_COMPOSITION_ONLY,
        "composition_only_at_fixed_boundary_temperature",
        "hot_rich_composition_only_finite_pulse",
        HOT_RICH_BOUNDARY,
        (
            "gate independent Maxwell-Stefan hexane flux J_h direction, "
            "relative step response, and gas entropy production; conserved "
            "component flux N_h and total Stefan flux N_t are diagnostic only"
        ),
    ),
    SCENARIO_TEMPERATURE_ONLY: StressScenario(
        SCENARIO_TEMPERATURE_ONLY,
        "temperature_only_at_fixed_boundary_composition",
        "cool_lean_temperature_only_finite_pulse",
        COOL_LEAN_BOUNDARY,
        (
            "gate signed boundary-to-film-surface thermal-driving-force response, "
            "independently reconstructed exact production conductive heat-flux "
            "response/sign, nonzero local-equilibrium storage response, and gas "
            "entropy; do not prescribe outer dry-cell-centre temperature or "
            "storage direction and do not attribute a component-flux response "
            "to Maxwell-Stefan diffusion"
        ),
    ),
    SCENARIO_COMBINED: StressScenario(
        SCENARIO_COMBINED,
        "correlated_temperature_and_composition",
        "cool_rich_finite_pulse",
        COOL_RICH_BOUNDARY,
        (
            "gate the raw-primitive Amendment-15 independent n-hexane "
            "counterflux J_h history response; conserved component flux N_h "
            "and total Stefan flux N_t remain ledger/transport diagnostics"
        ),
    ),
}
COMBINED_SCENARIO = SCENARIOS[SCENARIO_COMBINED]
# Backward-compatible alias for existing combined-only artifacts and imports.
DYNAMIC_BOUNDARY_SEGMENTS = COMBINED_SCENARIO.segments

# This N=2 root is a regression oracle for the unchanged dt=0.1 s equations.
N2_ROOT = cut.CutTransportUnknowns(
    wet_temperatures_k=(335.16876612065596, 336.2174581674662),
    wet_retained_water_loadings=(0.09999999999927532, 0.09999873535816925),
    dry_temperatures_k=(343.04902583944363,),
    dry_y_hexane=(0.8196136171414955,),
    dry_total_stefan_fluxes_mol_m2_s=(
        9.252425167926228e-05,
        0.0019907859101169433,
    ),
    front_z=0.21599507825461278,
    interface_temperature_k=336.97604077678653,
)


@dataclass(frozen=True)
class PhaseActiveSetCertificate:
    """Independent authority required before a rejected solve selects a phase."""

    component: str
    accepted_limit_state_certified: bool
    complementarity_condition_certified: bool
    source: str

    def __post_init__(self) -> None:
        if self.component not in {"water", "hexane"}:
            raise ValueError("phase certificate component must be water or hexane")
        if not self.accepted_limit_state_certified:
            raise ValueError("phase authority needs an accepted-limit state")
        if not self.complementarity_condition_certified:
            raise ValueError("phase authority needs a complementarity certificate")
        if not self.source.strip():
            raise ValueError("phase authority needs an auditable source")


@dataclass(frozen=True)
class AttemptRecord:
    role: str
    predictor_fraction: float
    wall_time_s: float
    accepted: bool
    nonlinear_evaluations: int
    maximum_scaled_residual: float
    condition_proxy: float
    rollback_identity_preserved: bool
    error: str | None = None
    # A same-cell solver can observe a trial candidate beyond the current cut
    # chart.  That is not an event-location result and carries no authority to
    # change topology.  The established face-event orchestrator owns any
    # accepted endpoint/bracket certificate and subsequent rank change.
    event_restart_required: bool = False
    face_event_trial_topology_excursion_observed: bool = False
    external_mobile_liquid_active_set_required: bool = False
    external_free_water_active_set_required: bool = False
    external_mobile_liquid_trial_excursion_observed: bool = False
    external_free_water_trial_excursion_observed: bool = False
    retained_water_applicability_guard_failed: bool = False
    phase_active_set_certificates: tuple[PhaseActiveSetCertificate, ...] = ()
    scaled_residual_block_maxima: tuple[tuple[str, float], ...] = ()
    macro_branch_status: str = "not_applicable"
    provisional_work_budget: int | None = None
    provisional_work_budget_exhausted: bool = False
    provisional_seed_attempt_budget_exhausted: bool = False
    attempted_dt_s: float | None = None
    adaptive_depth: int | None = None
    adaptive_path: str | None = None
    pre_adaptive_macro_branch_status: str | None = None
    event_failure_code: str | None = None
    event_failure_stage_type: str | None = None
    event_failure_nonlinear_evaluations: int | None = None
    event_failure_rejected_trial_evaluations: int | None = None
    event_failure_maximum_scaled_residual: float | None = None
    event_failure_condition_proxy: float | None = None

    def __post_init__(self) -> None:
        if self.attempted_dt_s is not None and (
            not math.isfinite(self.attempted_dt_s) or self.attempted_dt_s <= 0.0
        ):
            raise ValueError("attempted dt must be positive and finite")
        if self.adaptive_depth is not None and (
            isinstance(self.adaptive_depth, bool)
            or not isinstance(self.adaptive_depth, int)
            or self.adaptive_depth < 0
        ):
            raise ValueError("adaptive attempt depth must be a nonnegative integer")
        if self.adaptive_path is not None and not self.adaptive_path:
            raise ValueError("adaptive attempt path cannot be empty")
        if (
            self.event_failure_nonlinear_evaluations is not None
            and self.event_failure_nonlinear_evaluations < 0
        ) or (
            self.event_failure_rejected_trial_evaluations is not None
            and self.event_failure_rejected_trial_evaluations < 0
        ):
            raise ValueError("event failure evaluation counts cannot be negative")
        if self.provisional_work_budget is not None and (
            isinstance(self.provisional_work_budget, bool)
            or not isinstance(self.provisional_work_budget, int)
            or self.provisional_work_budget < 1
        ):
            raise ValueError("provisional nonlinear work budget must be positive")
        if self.provisional_work_budget_exhausted and self.provisional_work_budget is None:
            raise ValueError("provisional work-budget exhaustion needs its declared budget")
        if self.accepted and (
            self.provisional_work_budget_exhausted or self.provisional_seed_attempt_budget_exhausted
        ):
            raise ValueError("an accepted attempt cannot exhaust a recovery budget")
        if self.event_restart_required:
            raise ValueError(
                "a rejected same-cell attempt cannot certify a face event; "
                "route an accepted endpoint through the face-event orchestrator"
            )
        required_components = {
            *(("hexane",) if self.external_mobile_liquid_active_set_required else ()),
            *(("water",) if self.external_free_water_active_set_required else ()),
        }
        certificate_components = {
            certificate.component for certificate in self.phase_active_set_certificates
        }
        if len(certificate_components) != len(self.phase_active_set_certificates):
            raise ValueError("phase active-set certificates must be component-unique")
        if required_components != certificate_components:
            raise ValueError(
                "phase-active-set requirement needs a matching independent "
                "accepted-limit/complementarity certificate"
            )


@dataclass(frozen=True)
class AdaptiveSameCellMacroLedger:
    """Independent raw ledger for one all-or-nothing subdivided macrostep."""

    dt_s: float
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    nonlinear_message: str
    maximum_scaled_residual: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    minimum_front_z_distance_to_chart_boundary: float
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
    wet_retained_cap_active_piece_indices: tuple[int, ...]
    wet_retained_cap_active_piece_count: int
    wet_retained_cap_maximum_loading: float
    wet_retained_cap_maximum_capacity_dual_over_rt: float
    wet_retained_cap_maximum_complementarity_product: float
    wet_retained_cap_exact_graph_without_tolerance: bool
    wet_retained_cap_piece_transition_history: tuple[
        tuple[wrc.WetRetainedCapTransition, ...],
        ...,
    ]
    wet_retained_cap_entry_piece_count: int
    wet_retained_cap_continuation_piece_count: int
    wet_retained_cap_exit_piece_count: int

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
class AdaptiveSameCellMacrostep:
    """Certified requested interval composed from provisional binary leaves."""

    before: ci.CutIntegratorState
    after: ci.CutIntegratorState
    assembly: cut.CutTransportAssembly
    boundary: ct.DirichletPoreBoundary
    seed: ci.CutStepSeed
    ledger: AdaptiveSameCellMacroLedger
    substeps: tuple[ci.CutIntegratorStep, ...]
    binary_depth: int
    attempted_binary_depths: tuple[int, ...]
    predictor_backtrack_maximum_binary_exponent: int
    initial_failure_audit: dict | None = dataclass_field(
        default=None,
        compare=False,
        repr=False,
    )
    physically_qualifying: bool = dataclass_field(default=False, init=False)


@dataclass(frozen=True)
class AcceptedFaceEventMacrostep:
    """Trajectory-facing view of one atomically accepted interior-face step.

    The underlying orchestrator payload remains authoritative.  This adapter
    only exposes the common diagnostics consumed by the existing dynamic
    campaign and cannot manufacture an event, state, or balance interval.
    """

    event_macrostep: event_orchestrator.CutEventMacrostep
    bracket: (
        event_orchestrator.AcceptedHistoryFaceEventBracket
        | cc.MasterFaceEventContinuation
    )
    ledger: AdaptiveSameCellMacroLedger
    physically_qualifying: bool = dataclass_field(default=False, init=False)

    def __post_init__(self) -> None:
        macro = self.event_macrostep
        if macro.before is not self.bracket.before:
            raise ValueError("trajectory face adapter lost the bracketed input")
        if macro.ledger.branch is not (
            event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE
        ):
            raise ValueError("trajectory requires arrival plus strict departure")
        if not isinstance(macro.after, ci.CutIntegratorState):
            raise TypeError("trajectory face route did not return a strict cut state")
        if macro.face_departure_step is None:
            raise ValueError("trajectory face route lost its departure stage")
        if self.ledger.dt_s != macro.ledger.requested_dt_s:
            raise ValueError("trajectory face adapter changed the requested duration")
        if (
            max(
                self.ledger.maximum_step_ledger_residual,
                self.ledger.maximum_cumulative_ledger_residual,
            )
            > CONSERVATION_LIMIT
        ):
            raise ValueError("trajectory face adapter exceeds conservation contract")

    @property
    def before(self) -> ci.CutIntegratorState:
        return self.event_macrostep.before

    @property
    def after(self) -> ci.CutIntegratorState:
        after = self.event_macrostep.after
        if not isinstance(after, ci.CutIntegratorState):
            raise RuntimeError("accepted event macrostep lost strict departure")
        return after

    @property
    def assembly(self):
        departure = self.event_macrostep.face_departure_step
        if departure is None:
            raise RuntimeError("accepted event macrostep lost departure assembly")
        return departure.assembly

    @property
    def boundary(self) -> ct.PoreBoundary:
        return self.event_macrostep.boundary

    @property
    def seed(self) -> ci.CutStepSeed:
        return self.bracket.same_cell_seed


AcceptedSameCellStep = ci.CutIntegratorStep | AdaptiveSameCellMacrostep | AcceptedFaceEventMacrostep


@dataclass(frozen=True)
class ProvisionalAcceptedLeafSummary:
    """Immutable scalar evidence for a leaf discarded by macro rollback.

    No state, assembly, or mutable solver object is retained here.  The
    summary therefore cannot be mistaken for a committed trajectory node.
    """

    adaptive_depth: int
    adaptive_path: str
    starts_at_reversal: bool
    time_before_s: float
    time_after_s: float
    dt_s: float
    front_z: float
    cut_cell_index: int
    inner_face_z: float
    outer_face_z: float
    nearest_master_face_distance_z: float
    maximum_wet_retained_water_loading: float
    maximum_loading_wet_piece_position: int
    maximum_loading_master_cell_index: int
    retained_water_chart_upper_loading: float
    distance_to_retained_water_chart_upper: float
    retained_water_guard_loading: float | None
    retained_water_guard_margin: float | None
    wet_retained_cap_active_piece_indices: tuple[int, ...]
    wet_retained_cap_active_piece_count: int
    wet_retained_cap_maximum_capacity_dual_over_rt: float
    wet_retained_cap_maximum_complementarity_product: float
    wet_retained_cap_exact_graph_without_tolerance: bool
    wet_retained_cap_piece_transitions: tuple[
        wrc.WetRetainedCapTransition,
        ...,
    ]
    wet_retained_cap_entry_piece_count: int
    wet_retained_cap_continuation_piece_count: int
    wet_retained_cap_exit_piece_count: int
    surface_conserved_water_flux_mol_m2_s: float
    surface_conserved_hexane_flux_mol_m2_s: float
    surface_independent_hexane_flux_mol_m2_s: float
    surface_total_stefan_flux_mol_m2_s: float
    maximum_scaled_residual: float
    scaled_residual_block_maxima: tuple[tuple[str, float], ...]
    condition_proxy: float
    nonlinear_evaluations: int
    water_step_residual_mol: float
    hexane_step_residual_mol: float
    energy_step_residual_j: float
    maximum_step_ledger_residual: float
    maximum_cumulative_ledger_residual: float
    surface_fast_mass_departure: float | None
    surface_fast_mass_guard_passed: bool | None
    amendment_18_19_diagnostic_execution_evidence: dict[str, object] | None = None
    amendment_18_19_diagnostic_execution_evidence_contract_passed: bool = False
    provisional_only: bool = dataclass_field(default=True, init=False)
    contains_committable_state_object: bool = dataclass_field(
        default=False,
        init=False,
    )


class SameCellRecoveryError(RuntimeError):
    """A requested same-cell solve failed without committing its old state."""

    def __init__(
        self,
        message: str,
        *,
        rollback_state: ci.CutIntegratorState,
        attempts: tuple[AttemptRecord, ...],
        last_step_error: ci.CutIntegratorStepError | None = None,
        failure_jacobian_audit: dict | None = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.attempts = attempts
        self.last_step_error = last_step_error
        self.failure_jacobian_audit = failure_jacobian_audit


class AcceptedHistoryFaceRouteRecoveryError(SameCellRecoveryError):
    """Uncertified event work converted to same-cell-only recovery input."""

    def __init__(
        self,
        message: str,
        *,
        rollback_state: ci.CutIntegratorState,
        attempts: tuple[AttemptRecord, AttemptRecord],
        last_step_error: ci.CutIntegratorStepError,
        failure_jacobian_audit: dict,
        bracket: event_orchestrator.AcceptedHistoryFaceEventBracket,
        macrostep_error: event_orchestrator.CutEventMacrostepError,
        original_same_cell_seed: ci.CutStepSeed,
    ) -> None:
        if rollback_state is not bracket.before:
            raise ValueError("face-route recovery lost the bracket input state")
        if macrostep_error.rollback_state is not rollback_state:
            raise ValueError("face-route recovery lost atomic rollback identity")
        if macrostep_error.same_cell_error is not last_step_error or (
            last_step_error.rollback_state is not rollback_state
        ):
            raise ValueError("face-route recovery lost the original same-cell failure")
        if original_same_cell_seed is not bracket.same_cell_seed:
            raise ValueError("face-route recovery changed the original same-cell seed")
        if len(attempts) != 2 or any(attempt.accepted for attempt in attempts):
            raise ValueError("face-route recovery needs exactly two rejected records")
        if not all(attempt.rollback_identity_preserved for attempt in attempts):
            raise ValueError("face-route recovery requires exact rollback records")
        super().__init__(
            message,
            rollback_state=rollback_state,
            attempts=attempts,
            last_step_error=last_step_error,
            failure_jacobian_audit=failure_jacobian_audit,
        )
        self.bracket = bracket
        self.macrostep_error = macrostep_error
        self.original_same_cell_seed = original_same_cell_seed
        self.event_candidate_entered_recovery = False
        self.topology_selected = False


class RollbackIdentityViolation(RuntimeError):
    """A failed solver returned anything except its exact immutable input."""

    def __init__(self, context: str, *, expected_input_state, rollback_state) -> None:
        super().__init__(
            f"{context}: failed solver did not return the exact input state; recovery is forbidden"
        )
        self.expected_input_state = expected_input_state
        self.rollback_state = rollback_state


def _require_exact_rollback_identity(
    expected_input_state,
    rollback_state,
    *,
    context: str,
) -> None:
    if rollback_state is not expected_input_state:
        raise RollbackIdentityViolation(
            context,
            expected_input_state=expected_input_state,
            rollback_state=rollback_state,
        )


def _require_failed_attempt_rollbacks(
    attempts: Sequence[AttemptRecord],
    *,
    context: str,
) -> None:
    if any(
        not attempt.accepted and not attempt.rollback_identity_preserved for attempt in attempts
    ):
        raise RollbackIdentityViolation(
            context,
            expected_input_state=None,
            rollback_state=None,
        )


class AdaptiveSameCellMacrostepError(RuntimeError):
    """Binary recovery exhausted or met a licensed non-same-cell transition."""

    def __init__(
        self,
        message: str,
        *,
        requested_input_state: ci.CutIntegratorState,
        rollback_state: ci.CutIntegratorState,
        requested_dt_s: float,
        minimum_leaf_dt_s: float,
        attempted_binary_depths: tuple[int, ...],
        attempts: tuple[AttemptRecord, ...],
        classification: str,
        last_step_error: ci.CutIntegratorStepError | None = None,
        initial_failure_audit: dict | None = None,
        last_auditable_failure_jacobian_audit: dict | None = None,
        provisional_accepted_leaf_summaries: tuple[ProvisionalAcceptedLeafSummary, ...] = (),
        predictor_backtrack_maximum_binary_exponent: int = (
            ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
        ),
    ) -> None:
        _require_exact_rollback_identity(
            requested_input_state,
            rollback_state,
            context="adaptive macrostep failure",
        )
        super().__init__(message)
        self.requested_input_state = requested_input_state
        self.rollback_state = rollback_state
        self.requested_dt_s = requested_dt_s
        self.minimum_leaf_dt_s = minimum_leaf_dt_s
        self.attempted_binary_depths = attempted_binary_depths
        self.attempts = attempts
        self.classification = classification
        self.last_step_error = last_step_error
        self.initial_failure_audit = initial_failure_audit
        self.last_auditable_failure_jacobian_audit = last_auditable_failure_jacobian_audit
        if not all(
            isinstance(item, ProvisionalAcceptedLeafSummary)
            for item in provisional_accepted_leaf_summaries
        ):
            raise TypeError("provisional leaf evidence must use immutable summaries")
        self.provisional_accepted_leaf_summaries = provisional_accepted_leaf_summaries
        self.predictor_backtrack_maximum_binary_exponent = (
            predictor_backtrack_maximum_binary_exponent
        )

    @property
    def rollback_identity_preserved(self) -> bool:
        return self.rollback_state is self.requested_input_state

    @property
    def failed_state_committed(self) -> bool:
        return self.rollback_state is not self.requested_input_state


@dataclass(frozen=True)
class BaselinePrefix:
    """One immutable accepted hot/lean prefix shared before branch forcing."""

    last_step: AcceptedSameCellStep
    accepted_steps: tuple[AcceptedSameCellStep, ...]
    step_records: tuple[dict, ...]
    incremental_wall_time_s: float
    incremental_nonlinear_evaluations: int
    # Report-only: the seed route of a face-route time-refined first step
    # (DTDC_P1EF_FACE_ROUTE_TIME_REFINED_PREFIX); None on the on-disk path.
    face_route_time_refined_seed_route: dict | None = None


@dataclass(frozen=True)
class FixedRadiusFirstStepAudit:
    """Seed-only route to the unchanged fixed-radius first-step equations."""

    radius_m: float
    direct_representative_seed_accepted: bool
    direct_rejection_rollback_identity_preserved: bool | None
    seed_route: str
    seed_only_radius_stages_m: tuple[float, ...]
    cross_radius_seed_source_m: float | None = None
    cross_radius_attempts: tuple[AttemptRecord, ...] = ()
    amendment_18_19_diagnostic_execution_occurrences: tuple[
        dict[str, object], ...
    ] = ()
    accepted_state_carried_between_radii: bool = dataclass_field(
        default=False,
        init=False,
    )
    conserved_inventory_remapped_between_radii: bool = dataclass_field(
        default=False,
        init=False,
    )


@dataclass(frozen=True)
class CoupledFixedRadiusTrajectory:
    """Stable coupled-output view of one unchanged Gate-1g trajectory.

    The dimensional radius is an input to the actual grid and moving-front
    equations.  Radius continuation, when needed, supplies nonlinear seeds
    only: every accepted state starts from its own fixed-radius initial state.
    The raw diagnostic payload remains available for audit, but the named
    fields below are the backend contract consumed by radius studies.
    """

    radius_m: float
    cells: int
    time_step_s: float
    duration_s: float
    scenario_id: str
    forcing_history_id: str
    initial_front_position_over_radius: float
    final_front_position_over_radius: float
    minimum_front_position_over_radius: float
    maximum_front_position_over_radius: float
    initial_water_inventory_mol: float
    initial_hexane_inventory_mol: float
    initial_energy_inventory_j: float
    final_water_inventory_mol: float
    final_retained_water_inventory_mol: float
    final_pore_vapor_water_inventory_mol: float
    final_free_liquid_water_inventory_mol: float
    final_hexane_inventory_mol: float
    final_energy_inventory_j: float
    cumulative_boundary_water_out_mol: float
    cumulative_boundary_hexane_out_mol: float
    cumulative_boundary_energy_out_j: float
    minimum_material_temperature_k: float
    peak_material_temperature_k: float
    maximum_instantaneous_radial_temperature_spread_k: float
    temperature_extrema_scope: str
    unrepresented_continuous_intra_step_extrema_claimed: bool
    maximum_normalized_water_residual: float
    maximum_normalized_hexane_residual: float
    maximum_normalized_energy_residual: float
    maximum_scaled_nonlinear_residual: float
    first_step_audit: FixedRadiusFirstStepAudit
    nonqualifying_coefficient_fixture: bool
    dirichlet_outer_trace_without_qualified_film: bool
    physically_qualifying: bool
    diagnostics: dict = dataclass_field(compare=False, repr=False)
    accepted_step_seed_candidates: tuple[
        cut.CutTransportUnknowns,
        ...,
    ] = dataclass_field(compare=False, repr=False)

    @property
    def maximum_normalized_component_residual(self) -> float:
        return max(
            self.maximum_normalized_water_residual,
            self.maximum_normalized_hexane_residual,
        )

    @property
    def energy_uptake_j(self) -> float:
        return math.fsum((self.final_energy_inventory_j, -self.initial_energy_inventory_j))

    @property
    def cumulative_heat_into_particle_j(self) -> float:
        return -self.cumulative_boundary_energy_out_j


class CrossRadiusSeedTrajectoryError(RuntimeError):
    """One capped cross-radius seed attempt failed with exact rollback."""

    def __init__(
        self,
        message: str,
        *,
        target_radius_m: float,
        source_radius_m: float,
        attempts: tuple[AttemptRecord, ...],
        requested_input_state: ci.CutIntegratorState,
        rollback_state: ci.CutIntegratorState,
        last_step_error: ci.CutIntegratorStepError | None = None,
        dt_s: float | None = None,
        boundary: ct.DirichletPoreBoundary | None = None,
        seed: ci.CutStepSeed | None = None,
    ) -> None:
        _require_exact_rollback_identity(
            requested_input_state,
            rollback_state,
            context="cross-radius seed failure",
        )
        super().__init__(message)
        self.target_radius_m = target_radius_m
        self.source_radius_m = source_radius_m
        self.attempts = attempts
        self.requested_input_state = requested_input_state
        self.last_step_error = last_step_error
        self.rollback_state = rollback_state
        self.dt_s = dt_s
        self.boundary = boundary
        self.seed = seed

    @property
    def rollback_identity_preserved(self) -> bool:
        return self.rollback_state is self.requested_input_state

    @property
    def failed_state_committed(self) -> bool:
        return self.rollback_state is not self.requested_input_state


def _validate_radius(radius_m: float) -> float:
    if not math.isfinite(radius_m) or radius_m <= 0.0:
        raise ValueError("particle radius must be positive and finite")
    return radius_m


def initial_integrator(
    cells: int,
    radius_m: float | None = None,
) -> ci.CutIntegratorState:
    """Construct the unchanged Gate-1g initial state at one fixed radius."""

    if radius_m is None:
        radius_m = RADIUS_M
    radius_m = _validate_radius(radius_m)
    grid = uniform_grid(cells, radius_m)
    geometry = cg.partition_master_grid(
        grid,
        radius_m=INITIAL_FRONT_RADIUS_FRACTION * radius_m,
    )
    layout = cut.layout_for_geometry(geometry)
    history = wet_core.activate(grid, 335.0, WET_MODEL.wet).loadings
    transport = cut.initialize_cut_state(
        geometry,
        CONFIG,
        [335.0] * layout.wet_piece_count,
        [0.10] * layout.wet_piece_count,
        [343.0] * layout.dry_piece_count,
        [0.82] * layout.dry_piece_count,
        history,
    )
    return ci.initialize_integrator_state(
        transport,
        CONTROLS,
        interface_temperature_k=337.0,
    )


# Backward-compatible private spelling retained for checked-in campaign tests
# and artifacts created before the radius-aware backend was exposed.
_initial_integrator = initial_integrator


def _seed(candidate: cut.CutTransportUnknowns, label: str) -> ci.CutStepSeed:
    return ci.CutStepSeed(
        candidate,
        tuple(
            max(abs(value), STEFAN_SCALE_FLOOR)
            for value in candidate.dry_total_stefan_fluxes_mol_m2_s
        ),
        label,
    )


def _current_state_seed(
    previous: AcceptedSameCellStep,
    *,
    front_z: float | None,
    label: str,
) -> ci.CutStepSeed:
    """Keep every differential primitive current; optionally predict only front z."""

    state = previous.after
    candidate = cut.candidate_from_state(
        state.transport,
        dry_total_stefan_fluxes_mol_m2_s=(state.last_total_stefan_fluxes_mol_m2_s),
        interface_temperature_k=state.last_interface_temperature_k,
        front_z=front_z,
    )
    return _seed(candidate, label)


def _accepted_face_history_source(
    previous: AcceptedSameCellStep,
) -> event_orchestrator.AcceptedFaceHistoryStep:
    """Return the last concrete accepted front-motion interval."""

    if isinstance(previous, ci.CutIntegratorStep):
        return previous
    if isinstance(previous, AdaptiveSameCellMacrostep):
        if not previous.substeps:
            raise RuntimeError("adaptive macrostep lost its accepted leaf history")
        source = previous.substeps[-1]
        if source.after is not previous.after:
            raise RuntimeError("adaptive macrostep last leaf lost endpoint identity")
        return source
    if isinstance(previous, AcceptedFaceEventMacrostep):
        source = previous.event_macrostep.face_departure_step
        if source is None or source.after is not previous.after:
            raise RuntimeError("face macrostep lost accepted departure history")
        return source
    raise TypeError("unsupported accepted trajectory history")


def _accepted_history_front_motion(
    previous: AcceptedSameCellStep,
) -> tuple[float, float, float]:
    """Return accepted ``(start_z, end_z, dt)`` without rank interpolation."""

    source = _accepted_face_history_source(previous)
    if isinstance(source, ci.CutIntegratorStep):
        start_z = source.before.transport.geometry.front.z
    else:
        start_z = source.before.geometry.front.z
    end_z = source.after.transport.geometry.front.z
    return start_z, end_z, source.ledger.dt_s


def _current_state_anchored_front_seed(
    previous: AcceptedSameCellStep,
    *,
    next_dt_s: float,
) -> ci.CutStepSeed:
    """Use current stored fields and extrapolate only accepted front motion."""

    if not math.isfinite(next_dt_s) or next_dt_s <= 0.0:
        raise ValueError("current-state anchored seed needs a positive finite dt")
    old_front_z, current_front_z, history_dt_s = _accepted_history_front_motion(previous)
    ratio = next_dt_s / history_dt_s
    predicted_front_z = current_front_z + ratio * (current_front_z - old_front_z)
    return _current_state_seed(
        previous,
        front_z=predicted_front_z,
        label=(
            "exact current differential state with front-only accepted-motion "
            "predictor; no clipping or state interpolation"
        ),
    )


def _predict_from_step(
    step: AcceptedSameCellStep,
    *,
    target_boundary: ct.PoreBoundary,
    next_dt_s: float,
    predictor_fraction: float,
    label: str,
) -> ci.CutStepSeed:
    """Extrapolate a differential increment in its exact solver coordinates."""

    if not math.isfinite(predictor_fraction) or predictor_fraction <= 0.0:
        raise ValueError("predictor fraction must be positive and finite")
    if not isinstance(
        target_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("history prediction requires the actual target boundary")
    if isinstance(step, AcceptedFaceEventMacrostep):
        start_z, current_z, history_dt_s = _accepted_history_front_motion(step)
        candidate = cut.candidate_from_state(
            step.after.transport,
            dry_total_stefan_fluxes_mol_m2_s=(step.after.last_total_stefan_fluxes_mol_m2_s),
            interface_temperature_k=step.after.last_interface_temperature_k,
            front_z=math.fsum(
                (
                    current_z,
                    predictor_fraction * next_dt_s / history_dt_s * (current_z - start_z),
                )
            ),
        )
        result = _seed(
            candidate,
            (
                f"{label}; accepted face-departure front tangent fraction "
                f"{predictor_fraction:g}; current rank-safe primitives; no clipping"
            ),
        )
        ci._validate_seed(step.after, result, target_boundary)  # noqa: SLF001
        return result
    ratio = predictor_fraction * next_dt_s / step.ledger.dt_s
    old = step.before.transport
    new = step.after.transport

    def extrapolate(old_values: Sequence[float], new_values: Sequence[float]) -> tuple[float, ...]:
        return tuple(
            new_value + ratio * (new_value - old_value)
            for old_value, new_value in zip(old_values, new_values)
        )

    if step.after.controls.wet_water_bounds[1] == new.config.wet.luikov.W_cap:
        retained_cap_chart = wrc.WetRetainedCapSemismoothChart(
            step.after.controls.wet_water_bounds[0],
            new.config.wet.luikov,
        )
        old_retained_coordinates = tuple(
            retained_cap_chart.encode(loading, dual)
            for loading, dual in zip(
                old.wet_retained_water_loadings,
                old.effective_wet_retained_water_capacity_duals_over_rt,
            )
        )
        new_retained_coordinates = tuple(
            retained_cap_chart.encode(loading, dual)
            for loading, dual in zip(
                new.wet_retained_water_loadings,
                new.effective_wet_retained_water_capacity_duals_over_rt,
            )
        )
        retained_points = tuple(
            retained_cap_chart.decode(coordinate)
            for coordinate in extrapolate(
                old_retained_coordinates,
                new_retained_coordinates,
            )
        )
        predicted_retained_loadings = tuple(
            point.retained_water_loading for point in retained_points
        )
        predicted_retained_duals = tuple(point.capacity_dual_over_rt for point in retained_points)
        retained_predictor_label = "exact retained-cap graph-coordinate"
    else:
        old_retained_duals = old.effective_wet_retained_water_capacity_duals_over_rt
        new_retained_duals = new.effective_wet_retained_water_capacity_duals_over_rt
        if any((*old_retained_duals, *new_retained_duals)):
            raise RuntimeError("legacy smooth wet-water predictor cannot carry a retained-cap dual")
        predicted_retained_loadings = extrapolate(
            old.wet_retained_water_loadings,
            new.wet_retained_water_loadings,
        )
        predicted_retained_duals = (0.0,) * len(predicted_retained_loadings)
        retained_predictor_label = "legacy smooth physical-loading"

    predicted_dry_temperatures = extrapolate(
        old.dry_temperatures_k,
        new.dry_temperatures_k,
    )
    predicted_interface_temperature = step.after.last_interface_temperature_k + ratio * (
        step.after.last_interface_temperature_k
        - step.before.last_interface_temperature_k
    )
    predicted_dry_y, dry_composition_route = cc.affine_dry_composition_seed_values(
        step.after,
        source_surface_boundary=step.boundary,
        target_surface_boundary=target_boundary,
        old_dry_temperatures_k=old.dry_temperatures_k,
        old_interface_temperature_k=step.before.last_interface_temperature_k,
        old_dry_y_hexane=old.dry_y_hexane,
        new_dry_temperatures_k=new.dry_temperatures_k,
        new_interface_temperature_k=step.after.last_interface_temperature_k,
        new_dry_y_hexane=new.dry_y_hexane,
        predicted_dry_temperatures_k=predicted_dry_temperatures,
        predicted_interface_temperature_k=predicted_interface_temperature,
        affine_multiplier=1.0 + ratio,
    )

    candidate = cut.CutTransportUnknowns(
        wet_temperatures_k=extrapolate(
            old.wet_temperatures_k,
            new.wet_temperatures_k,
        ),
        wet_retained_water_loadings=predicted_retained_loadings,
        dry_temperatures_k=predicted_dry_temperatures,
        dry_y_hexane=predicted_dry_y,
        # Algebraic face rates are carried, not multiplied by a time ratio.
        dry_total_stefan_fluxes_mol_m2_s=(step.assembly.candidate.dry_total_stefan_fluxes_mol_m2_s),
        front_z=(new.geometry.front.z + ratio * (new.geometry.front.z - old.geometry.front.z)),
        interface_temperature_k=predicted_interface_temperature,
        wet_retained_water_capacity_duals_over_rt=predicted_retained_duals,
    )
    result = _seed(
        candidate,
        (
            f"{label}; unbounded physical-increment extrapolation fraction "
            f"{predictor_fraction:g}; {retained_predictor_label} predictor; "
            f"dry composition in {dry_composition_route}; no clipping"
        ),
    )
    # Fail explicitly if extrapolation left the exact chart.
    ci._validate_seed(step.after, result, target_boundary)  # noqa: SLF001
    return result


def _interpolate_boundary(
    left: ct.PoreBoundary,
    right: ct.PoreBoundary,
    fraction: float,
) -> ct.PoreBoundary:
    if not 0.0 < fraction <= 1.0:
        raise ValueError("boundary-homotopy fraction must lie in (0,1]")
    if type(left) is not type(right):
        raise TypeError("boundary homotopy cannot change the boundary-law class")
    if (
        replace(
            left,
            temperature_k=right.temperature_k,
            y_hexane=right.y_hexane,
            label=right.label,
        )
        != right
    ):
        raise ValueError(
            "boundary homotopy may vary only temperature and composition; "
            "film coefficients, guards, and active-set authority must be identical"
        )
    temperature = left.temperature_k + fraction * (right.temperature_k - left.temperature_k)
    phase_lower_temperature = DRY_CONFIG.conditioned_temperature_domain.solver_bounds_k[0]

    def coldest_face_temperature(boundary_temperature: float) -> float:
        return 0.5 * (phase_lower_temperature + boundary_temperature)

    left_coordinate = cp.encode_gas_only_y(
        coldest_face_temperature(left.temperature_k),
        PRESSURE_PA,
        left.y_hexane,
        DRY_CONFIG.pore,
    )
    right_coordinate = cp.encode_gas_only_y(
        coldest_face_temperature(right.temperature_k),
        PRESSURE_PA,
        right.y_hexane,
        DRY_CONFIG.pore,
    )
    coordinate = left_coordinate + fraction * (right_coordinate - left_coordinate)
    composition = cp.decode_gas_only_y(
        coldest_face_temperature(temperature),
        PRESSURE_PA,
        coordinate,
        DRY_CONFIG.pore,
    )
    return replace(
        right,
        temperature_k=temperature,
        y_hexane=composition,
        label=(
            f"nonqualifying exact-chart boundary seed homotopy fraction {fraction:g}; "
            "composition decoded at coldest possible symmetric face temperature"
        ),
    )


def _scaled_residual_block_maxima(
    layout: cut.CutTransportLayout,
    scaled_residuals: Sequence[float],
) -> tuple[tuple[str, float], ...]:
    """Group the signed solver residual vector by its declared equation blocks."""

    values = tuple(float(value) for value in scaled_residuals)
    if not values:
        return ()
    if len(values) != layout.residual_count:
        raise ValueError("scaled residual vector does not match the declared rank")
    cursor = 0
    maxima = []
    for name, size in layout.residual_blocks:
        block = values[cursor : cursor + size]
        cursor += size
        maxima.append((name, max(abs(value) for value in block)))
    return tuple(maxima)


def _reports_external_mobile_liquid_trial_excursion(message: str) -> bool:
    """Recognize an inadmissible trial without promoting it to phase authority."""

    return "activate the mobile-liquid topology, do not clamp" in message


def _reports_external_free_water_trial_excursion(message: str) -> bool:
    """Recognize an inadmissible water trial, not a certified phase transition."""

    return "activate the external/free-water topology, do not clamp" in message


def _candidate_physical_vector(
    candidate: cut.CutTransportUnknowns,
) -> np.ndarray:
    return np.asarray(
        (
            *candidate.wet_temperatures_k,
            *candidate.wet_retained_water_loadings,
            *candidate.dry_temperatures_k,
            *candidate.dry_y_hexane,
            *candidate.dry_total_stefan_fluxes_mol_m2_s,
            candidate.front_z,
            candidate.interface_temperature_k,
        ),
        dtype=float,
    )


def _candidate_from_physical_vector(
    values: Sequence[float],
    layout: cut.CutTransportLayout,
) -> cut.CutTransportUnknowns:
    physical = tuple(float(value) for value in values)
    if len(physical) != layout.unknown_count:
        raise ValueError("physical candidate vector does not match declared rank")
    nw = layout.wet_piece_count
    nd = layout.dry_piece_count
    nf = layout.dry_face_count
    cursor = 0

    def take(size: int) -> tuple[float, ...]:
        nonlocal cursor
        result = physical[cursor : cursor + size]
        cursor += size
        return result

    candidate = cut.CutTransportUnknowns(
        wet_temperatures_k=take(nw),
        wet_retained_water_loadings=take(nw),
        dry_temperatures_k=take(nd),
        dry_y_hexane=take(nd),
        dry_total_stefan_fluxes_mol_m2_s=take(nf),
        front_z=take(1)[0],
        interface_temperature_k=take(1)[0],
    )
    if cursor != len(physical):
        raise RuntimeError("physical candidate unpack did not consume its vector")
    return candidate


def _central_difference_jacobian(
    function,
    point: np.ndarray,
    column_scales: np.ndarray,
) -> np.ndarray:
    """Evaluate a bounded read-only Jacobian, shrinking invalid probes."""

    base = np.asarray(function(point), dtype=float)
    if base.ndim != 1 or not np.all(np.isfinite(base)):
        raise ValueError("Jacobian audit needs one finite residual vector")
    matrix = np.empty((base.size, point.size), dtype=float)
    for column in range(point.size):
        step = 1.0e-6 * float(column_scales[column])
        if not math.isfinite(step) or step <= 0.0:
            raise ValueError("Jacobian audit column scale must be positive and finite")
        plus = minus = None
        for _ in range(16):
            plus_point = point.copy()
            minus_point = point.copy()
            plus_point[column] += step
            minus_point[column] -= step
            try:
                plus = np.asarray(function(plus_point), dtype=float)
                minus = np.asarray(function(minus_point), dtype=float)
            except Exception:
                plus = minus = None
                step *= 0.5
                continue
            if np.all(np.isfinite(plus)) and np.all(np.isfinite(minus)):
                break
            plus = minus = None
            step *= 0.5
        if plus is None or minus is None:
            raise RuntimeError(
                f"Jacobian audit could not form an interior probe for column {column}"
            )
        matrix[:, column] = (plus - minus) / (2.0 * step)
    return matrix


def _matrix_rank_condition(jacobian: np.ndarray) -> dict:
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    rank = int(np.linalg.matrix_rank(jacobian))
    full_rank = min(jacobian.shape)
    condition = float(np.linalg.cond(jacobian))
    return {
        "rows": int(jacobian.shape[0]),
        "columns": int(jacobian.shape[1]),
        "rank": rank,
        "full_rank": full_rank,
        "is_full_rank": rank == full_rank,
        "condition_2": condition if math.isfinite(condition) else None,
        "condition_2_is_finite": math.isfinite(condition),
        "largest_singular_value": float(singular_values[0]),
        "smallest_singular_value": float(singular_values[-1]),
    }


def audit_failed_same_cell_candidate(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.DirichletPoreBoundary,
    seed: ci.CutStepSeed,
    error: ci.CutIntegratorStepError,
) -> dict:
    """Audit the last auditable rejected candidate without retrying it."""

    if error.rollback_state is not before:
        raise ValueError("failure audit needs the exact rollback state")
    candidate = error.last_candidate
    if candidate is None or not error.scaled_residuals:
        raise ValueError("failure audit needs a last candidate and scaled residuals")
    chart = ci._coordinate_chart(before, boundary)  # noqa: SLF001
    initial_assembly = cut.assemble_backward_euler(
        before.transport,
        seed.candidate,
        dt_s,
        boundary,
    )
    residual_scales = ci._residual_scales(  # noqa: SLF001
        before.transport,
        initial_assembly,
        dt_s,
    )
    row_scales = np.asarray(residual_scales.vector, dtype=float)

    def residual_for_candidate(value: cut.CutTransportUnknowns) -> np.ndarray:
        assembly = cut.assemble_backward_euler(
            before.transport,
            value,
            dt_s,
            boundary,
        )
        return (
            np.asarray(
                ci._datum_covariant_residual_vector(assembly.residuals),  # noqa: SLF001
                dtype=float,
            )
            / row_scales
        )

    layout = before.transport.layout
    retained_cap_certificate = wrc.certify_exact_graph(
        candidate.wet_retained_water_loadings,
        candidate.effective_wet_retained_water_capacity_duals_over_rt,
        before.transport.config.wet.luikov,
    )
    active_retained_cap = retained_cap_certificate.active_piece_count > 0
    if active_retained_cap:
        # W alone is not a lossless coordinate on the vertical normal-cone
        # branch.  A physical-W central difference would silently discard
        # lambda, so only the production graph-coordinate audit is valid here.
        physical_column_scales = None
        physical_jacobian = None
        normalized_physical_jacobian = None
    else:
        physical_point = _candidate_physical_vector(candidate)
        physical_column_scales = np.asarray(
            (
                *((chart.wet_temperature[1] - chart.wet_temperature[0],) * layout.wet_piece_count),
                *((chart.wet_water[1] - chart.wet_water[0],) * layout.wet_piece_count),
                *((chart.dry_temperature[1] - chart.dry_temperature[0],) * layout.dry_piece_count),
                *((1.0,) * layout.dry_piece_count),
                *seed.stefan_flux_scales_mol_m2_s,
                chart.front_z[1] - chart.front_z[0],
                chart.interface_temperature[1] - chart.interface_temperature[0],
            ),
            dtype=float,
        )
        physical_jacobian = _central_difference_jacobian(
            lambda values: residual_for_candidate(_candidate_from_physical_vector(values, layout)),
            physical_point,
            physical_column_scales,
        )
        normalized_physical_jacobian = physical_jacobian * physical_column_scales[np.newaxis, :]

    coordinate_point = ci._encode_candidate(  # noqa: SLF001
        candidate,
        chart,
        seed.stefan_flux_scales_mol_m2_s,
    )
    coordinate_jacobian = _central_difference_jacobian(
        lambda values: residual_for_candidate(
            ci._decode_candidate(  # noqa: SLF001
                values,
                layout,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
            )
        ),
        coordinate_point,
        np.maximum(np.abs(coordinate_point), 1.0),
    )
    reevaluated = residual_for_candidate(candidate)
    solver_residual = np.asarray(error.scaled_residuals, dtype=float)
    return {
        "audit_completed": True,
        "equation_rank": layout.residual_count,
        "unknown_rank": layout.unknown_count,
        "square_system": layout.is_square,
        "maximum_scaled_residual": float(np.max(np.abs(reevaluated))),
        "reevaluation_matches_solver_scaled_residual_max_abs": float(
            np.max(np.abs(reevaluated - solver_residual))
        ),
        "scaled_residual_block_maxima": dict(_scaled_residual_block_maxima(layout, reevaluated)),
        "solver_reported_condition_proxy": (
            error.condition_proxy if math.isfinite(error.condition_proxy) else None
        ),
        "solver_reported_condition_proxy_is_finite": math.isfinite(error.condition_proxy),
        "minimum_fractional_distance_to_bound": (error.minimum_fractional_distance_to_bound),
        "face_event_trial_topology_excursion_observed": (error.event_restart_required),
        "certified_face_event_restart_required": False,
        "physical_primitive_jacobian_raw_units": (
            None
            if physical_jacobian is None
            else {
                **_matrix_rank_condition(physical_jacobian),
                "condition_is_unit_dependent": True,
            }
        ),
        "physical_primitive_jacobian_column_scaled": (
            None
            if normalized_physical_jacobian is None or physical_column_scales is None
            else {
                **_matrix_rank_condition(normalized_physical_jacobian),
                "column_scales": physical_column_scales.tolist(),
            }
        ),
        "physical_primitive_jacobian_unavailable_reason": (
            "W-only coordinates discard the positive normal-cone dual"
            if active_retained_cap
            else None
        ),
        "solver_transformed_coordinate_jacobian": _matrix_rank_condition(coordinate_jacobian),
        "solver_transformed_coordinate_is_exact_retained_cap_graph": (active_retained_cap),
        "finite_difference_relative_step": 1.0e-6,
        "physics_or_acceptance_contract_changed": False,
    }


def _unavailable_failure_jacobian_audit(
    before: ci.CutIntegratorState,
    reason: str,
) -> dict:
    layout = before.transport.layout
    return {
        "audit_completed": False,
        "reason": reason,
        "equation_rank": layout.residual_count,
        "unknown_rank": layout.unknown_count,
        "square_system": layout.is_square,
        "physical_primitive_jacobian_raw_units": None,
        "physical_primitive_jacobian_column_scaled": None,
        "solver_transformed_coordinate_jacobian": None,
        "physics_or_acceptance_contract_changed": False,
    }


def _failure_jacobian_audit(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.DirichletPoreBoundary,
    seed: ci.CutStepSeed,
    error: ci.CutIntegratorStepError | None,
) -> dict:
    if error is None:
        return _unavailable_failure_jacobian_audit(
            before,
            "no evaluated nonlinear candidate was returned",
        )
    if error.last_candidate is None or not error.scaled_residuals:
        return _unavailable_failure_jacobian_audit(
            before,
            (
                "failure occurred during seed/trial evaluation before a valid "
                "last candidate and scaled residual vector existed"
            ),
        )
    try:
        return audit_failed_same_cell_candidate(
            before,
            dt_s,
            boundary,
            seed,
            error,
        )
    except Exception as exc:
        payload = _unavailable_failure_jacobian_audit(
            before,
            "finite-difference failure audit did not complete",
        )
        payload["error_type"] = type(exc).__name__
        payload["error"] = str(exc)
        return payload


def _try_step_detailed(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.DirichletPoreBoundary,
    seed: ci.CutStepSeed,
    *,
    role: str,
    predictor_fraction: float,
    provisional_maximum_function_evaluations: int | None = None,
) -> tuple[
    ci.CutIntegratorStep | None,
    AttemptRecord,
    ci.CutIntegratorStepError | None,
]:
    started = perf_counter()
    try:
        if provisional_maximum_function_evaluations is None:
            step = ci.advance_same_cell_backward_euler(
                before,
                dt_s,
                boundary,
                seed,
            )
        else:
            step = ci.advance_same_cell_backward_euler(
                before,
                dt_s,
                boundary,
                seed,
                provisional_maximum_function_evaluations=(provisional_maximum_function_evaluations),
            )
    except ci.CutIntegratorStepError as exc:
        _require_exact_rollback_identity(
            before,
            exc.rollback_state,
            context=f"{role} same-cell rejection",
        )
        message = str(exc)
        return (
            None,
            AttemptRecord(
                role=role,
                predictor_fraction=predictor_fraction,
                wall_time_s=perf_counter() - started,
                accepted=False,
                nonlinear_evaluations=exc.nonlinear_evaluations,
                maximum_scaled_residual=exc.maximum_scaled_residual,
                condition_proxy=exc.condition_proxy,
                rollback_identity_preserved=exc.rollback_state is before,
                error=message,
                face_event_trial_topology_excursion_observed=(exc.event_restart_required),
                external_mobile_liquid_trial_excursion_observed=(
                    _reports_external_mobile_liquid_trial_excursion(message)
                ),
                external_free_water_trial_excursion_observed=(
                    _reports_external_free_water_trial_excursion(message)
                ),
                scaled_residual_block_maxima=_scaled_residual_block_maxima(
                    before.transport.layout,
                    exc.scaled_residuals,
                ),
                provisional_work_budget=exc.provisional_work_budget,
                provisional_work_budget_exhausted=(exc.provisional_work_budget_exhausted),
                retained_water_applicability_guard_failed=(
                    exc.retained_water_applicability_guard_failed
                ),
            ),
            exc,
        )
    if step.before is not before:
        raise RuntimeError("accepted same-cell step does not reference its exact input")
    return (
        step,
        AttemptRecord(
            role=role,
            predictor_fraction=predictor_fraction,
            wall_time_s=perf_counter() - started,
            accepted=True,
            nonlinear_evaluations=step.ledger.nonlinear_evaluations,
            maximum_scaled_residual=step.ledger.maximum_scaled_residual,
            condition_proxy=step.ledger.condition_proxy,
            rollback_identity_preserved=step.before is before,
            scaled_residual_block_maxima=_scaled_residual_block_maxima(
                before.transport.layout,
                step.ledger.scaled_residuals,
            ),
            provisional_work_budget=provisional_maximum_function_evaluations,
        ),
        None,
    )


def _try_step(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.DirichletPoreBoundary,
    seed: ci.CutStepSeed,
    *,
    role: str,
    predictor_fraction: float,
) -> tuple[ci.CutIntegratorStep | None, AttemptRecord]:
    """Compatibility wrapper retaining the established two-result API."""

    step, record, _ = _try_step_detailed(
        before,
        dt_s,
        boundary,
        seed,
        role=role,
        predictor_fraction=predictor_fraction,
    )
    return step, record


def _seed_validation_attempt(
    *,
    role: str,
    predictor_fraction: float,
    error: Exception,
) -> AttemptRecord:
    message = str(error)
    return AttemptRecord(
        role=role,
        predictor_fraction=predictor_fraction,
        wall_time_s=0.0,
        accepted=False,
        nonlinear_evaluations=0,
        maximum_scaled_residual=math.inf,
        condition_proxy=math.inf,
        rollback_identity_preserved=True,
        error=message,
        external_mobile_liquid_trial_excursion_observed=(
            _reports_external_mobile_liquid_trial_excursion(message)
        ),
        external_free_water_trial_excursion_observed=(
            _reports_external_free_water_trial_excursion(message)
        ),
    )


def _advance_with_recovery(
    previous: AcceptedSameCellStep,
    dt_s: float,
    previous_boundary: ct.DirichletPoreBoundary,
    boundary: ct.DirichletPoreBoundary,
    *,
    is_reversal: bool,
    try_current_state_anchor: bool = False,
    predictor_backtrack_maximum_binary_exponent: int = (PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT),
    provisional_maximum_function_evaluations: int | None = None,
    maximum_nonlinear_seed_attempts: int | None = None,
) -> tuple[ci.CutIntegratorStep, tuple[AttemptRecord, ...], str]:
    """Solve one requested BE step; all recovery roots are seed-only."""

    if (
        isinstance(predictor_backtrack_maximum_binary_exponent, bool)
        or not isinstance(predictor_backtrack_maximum_binary_exponent, int)
        or predictor_backtrack_maximum_binary_exponent < 0
        or predictor_backtrack_maximum_binary_exponent > PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
    ):
        raise ValueError(
            "predictor backtrack maximum binary exponent must be an integer "
            f"in [0, {PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT}]"
        )
    if is_reversal and predictor_backtrack_maximum_binary_exponent < 1:
        raise ValueError("a reversal predictor budget must include exponent 1")
    if provisional_maximum_function_evaluations is not None and (
        isinstance(provisional_maximum_function_evaluations, bool)
        or not isinstance(provisional_maximum_function_evaluations, int)
        or provisional_maximum_function_evaluations < 1
        or provisional_maximum_function_evaluations
        > previous.after.controls.maximum_function_evaluations
    ):
        raise ValueError(
            "provisional nonlinear work budget must be an integer in [1, "
            f"{previous.after.controls.maximum_function_evaluations}]"
        )
    if maximum_nonlinear_seed_attempts is not None and (
        isinstance(maximum_nonlinear_seed_attempts, bool)
        or not isinstance(maximum_nonlinear_seed_attempts, int)
        or maximum_nonlinear_seed_attempts < 1
    ):
        raise ValueError("nonlinear seed-attempt budget must be a positive integer")
    if maximum_nonlinear_seed_attempts is not None:
        frozen_attempt_count = (
            DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
            if is_reversal
            else DISCONTINUITY_NONREVERSAL_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
        )
        if maximum_nonlinear_seed_attempts != frozen_attempt_count:
            raise ValueError(
                "bounded discontinuity recovery requires the frozen ordered "
                f"seed count {frozen_attempt_count}"
            )

    attempts: list[AttemptRecord] = []
    nonlinear_seed_attempts = 0
    last_valid_seed: ci.CutStepSeed | None = None
    last_step_error: ci.CutIntegratorStepError | None = None
    last_failure_context: (
        tuple[
            ci.CutStepSeed,
            ct.DirichletPoreBoundary,
            ci.CutIntegratorStepError,
        ]
        | None
    ) = None
    last_auditable_context: (
        tuple[
            ci.CutStepSeed,
            ct.DirichletPoreBoundary,
            ci.CutIntegratorStepError,
        ]
        | None
    ) = None

    def remember_failure(
        failed_seed: ci.CutStepSeed,
        failed_boundary: ct.DirichletPoreBoundary,
        failed_error: ci.CutIntegratorStepError | None,
    ) -> None:
        nonlocal last_failure_context, last_auditable_context
        if failed_error is None:
            return
        last_failure_context = (failed_seed, failed_boundary, failed_error)
        if failed_error.last_candidate is not None and failed_error.scaled_residuals:
            last_auditable_context = last_failure_context

    def failure_audit() -> dict:
        context = last_auditable_context or last_failure_context
        if context is None:
            return _unavailable_failure_jacobian_audit(
                previous.after,
                "no evaluated nonlinear failure candidate existed",
            )
        failed_seed, failed_boundary, failed_error = context
        return _failure_jacobian_audit(
            previous.after,
            dt_s,
            failed_boundary,
            failed_seed,
            failed_error,
        )

    def try_seed(
        failed_seed: ci.CutStepSeed,
        failed_boundary: ct.DirichletPoreBoundary,
        *,
        role: str,
        predictor_fraction: float,
    ) -> tuple[
        ci.CutIntegratorStep | None,
        AttemptRecord,
        ci.CutIntegratorStepError | None,
    ]:
        nonlocal nonlinear_seed_attempts, last_step_error
        if (
            maximum_nonlinear_seed_attempts is not None
            and nonlinear_seed_attempts >= maximum_nonlinear_seed_attempts
        ):
            raise RuntimeError("internal nonlinear seed budget check was bypassed")
        nonlinear_seed_attempts += 1
        if provisional_maximum_function_evaluations is None:
            result = _try_step_detailed(
                previous.after,
                dt_s,
                failed_boundary,
                failed_seed,
                role=role,
                predictor_fraction=predictor_fraction,
            )
        else:
            result = _try_step_detailed(
                previous.after,
                dt_s,
                failed_boundary,
                failed_seed,
                role=role,
                predictor_fraction=predictor_fraction,
                provisional_maximum_function_evaluations=(provisional_maximum_function_evaluations),
            )
        _, record, step_error = result
        if step_error is not None and step_error.retained_water_applicability_guard_failed:
            attempts.append(record)
            last_step_error = step_error
            remember_failure(failed_seed, failed_boundary, step_error)
            raise SameCellRecoveryError(
                "accepted root exceeded the selected smooth wet retained-water "
                "applicability guard; alternate seeds and homotopy are forbidden",
                rollback_state=previous.after,
                attempts=tuple(attempts),
                last_step_error=step_error,
                failure_jacobian_audit=failure_audit(),
            )
        return result

    def seed_attempt_budget_is_exhausted() -> bool:
        return bool(
            maximum_nonlinear_seed_attempts is not None
            and nonlinear_seed_attempts >= maximum_nonlinear_seed_attempts
        )

    def raise_seed_attempt_budget_exhausted() -> None:
        message = (
            "provisional nonlinear seed-attempt budget exhausted after "
            f"{nonlinear_seed_attempts} unchanged-equation attempt(s)"
        )
        attempts.append(
            AttemptRecord(
                role="provisional_nonlinear_seed_attempt_budget_exhausted",
                predictor_fraction=0.0,
                wall_time_s=0.0,
                accepted=False,
                nonlinear_evaluations=0,
                maximum_scaled_residual=(
                    math.inf if last_step_error is None else last_step_error.maximum_scaled_residual
                ),
                condition_proxy=(
                    math.inf if last_step_error is None else last_step_error.condition_proxy
                ),
                rollback_identity_preserved=True,
                error=message,
                provisional_work_budget=(provisional_maximum_function_evaluations),
                provisional_seed_attempt_budget_exhausted=True,
            )
        )
        raise SameCellRecoveryError(
            message,
            rollback_state=previous.after,
            attempts=tuple(attempts),
            last_step_error=last_step_error,
            failure_jacobian_audit=failure_audit(),
        )

    def try_current_state_anchor_fallback() -> ci.CutIntegratorStep | None:
        """Try the current-state/front-only seed after the primary tangent seed.

        A backward-Euler residual can have more than one admissible nonlinear
        root.  The temporal tangent predictor is therefore always tried first:
        it is the discrete continuation of the already accepted trajectory.
        This auxiliary seed is a globalization fallback only and must never
        replace a root already certified from that history-continuous seed.
        """

        nonlocal last_valid_seed, last_step_error
        try:
            anchored_seed = _current_state_anchored_front_seed(
                previous,
                next_dt_s=dt_s,
            )
            ci._validate_seed(  # noqa: SLF001
                previous.after,
                anchored_seed,
                boundary,
            )
        except ValueError as exc:
            attempts.append(
                _seed_validation_attempt(
                    role=("current_state_anchored_front_seed_fallback_chart_validation"),
                    predictor_fraction=1.0,
                    error=exc,
                )
            )
            return None
        else:
            last_valid_seed = anchored_seed
            anchored_step, record, step_error = try_seed(
                anchored_seed,
                boundary,
                role=("current_state_anchored_front_seed_fallback_backward_euler"),
                predictor_fraction=1.0,
            )
            attempts.append(record)
            _require_failed_attempt_rollbacks(
                (record,),
                context="current-state anchored recovery attempt",
            )
            last_step_error = step_error
            remember_failure(anchored_seed, boundary, step_error)
            if anchored_step is not None:
                return anchored_step
            if record.external_mobile_liquid_active_set_required:
                raise SameCellRecoveryError(
                    "same-cell recovery reached the certified external mobile-liquid active set",
                    rollback_state=previous.after,
                    attempts=tuple(attempts),
                    last_step_error=step_error,
                    failure_jacobian_audit=failure_audit(),
                )
            if record.external_free_water_active_set_required:
                raise SameCellRecoveryError(
                    "same-cell recovery reached the certified external free-water active set",
                    rollback_state=previous.after,
                    attempts=tuple(attempts),
                    last_step_error=step_error,
                    failure_jacobian_audit=failure_audit(),
                )
            return None

    first_exponent = 1 if is_reversal else 0
    last_exponent = (
        1
        if maximum_nonlinear_seed_attempts is not None
        else predictor_backtrack_maximum_binary_exponent
    )
    fractions = tuple(
        math.ldexp(1.0, -exponent)
        for exponent in range(
            first_exponent,
            last_exponent + 1,
        )
    )
    for predictor_index, fraction in enumerate(fractions):
        try:
            seed = _predict_from_step(
                previous,
                target_boundary=boundary,
                next_dt_s=dt_s,
                predictor_fraction=fraction,
                label="requested-step predictor",
            )
        except ValueError as exc:
            attempts.append(
                _seed_validation_attempt(
                    role="predictor_chart_validation",
                    predictor_fraction=fraction,
                    error=exc,
                )
            )
        else:
            last_valid_seed = seed
            if seed_attempt_budget_is_exhausted():
                raise_seed_attempt_budget_exhausted()
            step, record, step_error = try_seed(
                seed,
                boundary,
                role="requested_backward_euler",
                predictor_fraction=fraction,
            )
            attempts.append(record)
            _require_failed_attempt_rollbacks(
                (record,),
                context="requested-step predictor recovery attempt",
            )
            last_step_error = step_error
            remember_failure(seed, boundary, step_error)
            if step is not None:
                return step, tuple(attempts), "direct_exact"
            if record.external_mobile_liquid_active_set_required:
                raise SameCellRecoveryError(
                    "same-cell recovery reached the external mobile-liquid active set",
                    rollback_state=previous.after,
                    attempts=tuple(attempts),
                    last_step_error=step_error,
                    failure_jacobian_audit=failure_audit(),
                )
            if record.external_free_water_active_set_required:
                raise SameCellRecoveryError(
                    "same-cell recovery reached the external free-water active set",
                    rollback_state=previous.after,
                    attempts=tuple(attempts),
                    last_step_error=step_error,
                    failure_jacobian_audit=failure_audit(),
                )

        # The first predictor is the primary temporal-continuation seed.  Only
        # after it fails may the current-state/front-only fallback compete for
        # this step; otherwise a secondary seed could silently select another
        # admissible implicit root and make the trajectory seed-order dependent.
        if predictor_index == 0 and try_current_state_anchor:
            if seed_attempt_budget_is_exhausted():
                raise_seed_attempt_budget_exhausted()
            anchored_step = try_current_state_anchor_fallback()
            if anchored_step is not None:
                return (
                    anchored_step,
                    tuple(attempts),
                    "current_state_anchored_fallback_exact",
                )
        if seed_attempt_budget_is_exhausted():
            raise_seed_attempt_budget_exhausted()

    if maximum_nonlinear_seed_attempts is not None:
        raise SameCellRecoveryError(
            "frozen provisional seed sequence exhausted without an accepted root",
            rollback_state=previous.after,
            attempts=tuple(attempts),
            last_step_error=last_step_error,
            failure_jacobian_audit=failure_audit(),
        )

    # Boundary continuation changes only auxiliary seed equations.  Every root
    # has the exact same immutable old state; only the fraction-one root is
    # returned and entered into the accepted trajectory.
    if last_valid_seed is None:
        raise SameCellRecoveryError(
            "no strict-open physical-increment predictor exists for boundary homotopy",
            rollback_state=previous.after,
            attempts=tuple(attempts),
            failure_jacobian_audit=failure_audit(),
        )
    seed = last_valid_seed
    final_step = None
    for fraction in (0.25, 0.5, 0.75, 1.0):
        if seed_attempt_budget_is_exhausted():
            raise_seed_attempt_budget_exhausted()
        trial_boundary = (
            boundary
            if fraction == 1.0
            else _interpolate_boundary(
                previous_boundary,
                boundary,
                fraction,
            )
        )
        role = (
            "requested_backward_euler_after_boundary_homotopy"
            if fraction == 1.0
            else "auxiliary_same-old-state_boundary_homotopy_seed"
        )
        final_step, record, step_error = try_seed(
            seed,
            trial_boundary,
            role=role,
            predictor_fraction=fraction,
        )
        attempts.append(record)
        _require_failed_attempt_rollbacks(
            (record,),
            context="boundary-homotopy recovery attempt",
        )
        last_step_error = step_error
        remember_failure(seed, trial_boundary, step_error)
        if final_step is None:
            classification = (
                "external mobile-liquid active set"
                if record.external_mobile_liquid_active_set_required
                else (
                    "external free-water active set"
                    if record.external_free_water_active_set_required
                    else "nonlinear residual failure"
                )
            )
            raise SameCellRecoveryError(
                "boundary seed homotopy failed without changing the accepted state: "
                f"{classification}: {record.error}",
                rollback_state=previous.after,
                attempts=tuple(attempts),
                last_step_error=last_step_error,
                failure_jacobian_audit=failure_audit(),
            )
        seed = _seed(
            final_step.assembly.candidate,
            f"same-old-state boundary homotopy root at fraction {fraction:g}",
        )
    assert final_step is not None
    return final_step, tuple(attempts), "boundary_homotopy_exact"


def _advance_with_cross_radius_seed(
    previous: AcceptedSameCellStep,
    dt_s: float,
    boundary: ct.DirichletPoreBoundary,
    candidate: cut.CutTransportUnknowns,
    *,
    source_radius_m: float,
) -> tuple[ci.CutIntegratorStep, tuple[AttemptRecord, ...], str]:
    """Try exactly one neighboring-radius candidate on unchanged equations."""

    target_radius_m = previous.after.transport.geometry.master_grid.R
    seed = _seed(
        candidate,
        (
            f"certified {source_radius_m:g} m trajectory root reused only as "
            f"a nonlinear seed for independent {target_radius_m:g} m equations"
        ),
    )
    try:
        ci._validate_seed(previous.after, seed)  # noqa: SLF001
    except ValueError as exc:
        record = AttemptRecord(
            role="cross_radius_seed_chart_validation",
            predictor_fraction=1.0,
            wall_time_s=0.0,
            accepted=False,
            nonlinear_evaluations=0,
            maximum_scaled_residual=math.inf,
            condition_proxy=math.inf,
            rollback_identity_preserved=True,
            error=str(exc),
        )
        raise CrossRadiusSeedTrajectoryError(
            "cross-radius seed is outside the target's exact chart",
            target_radius_m=target_radius_m,
            source_radius_m=source_radius_m,
            attempts=(record,),
            requested_input_state=previous.after,
            rollback_state=previous.after,
            dt_s=dt_s,
            boundary=boundary,
            seed=seed,
        ) from exc
    step, record, step_error = _try_step_detailed(
        previous.after,
        dt_s,
        boundary,
        seed,
        role="requested_backward_euler_from_cross_radius_seed_only_candidate",
        predictor_fraction=1.0,
    )
    if step is None:
        if not record.rollback_identity_preserved:
            raise RuntimeError("cross-radius seed rejection did not roll back exactly")
        raise CrossRadiusSeedTrajectoryError(
            (
                "single capped cross-radius seed attempt failed without state "
                f"commit: {record.error}"
            ),
            target_radius_m=target_radius_m,
            source_radius_m=source_radius_m,
            attempts=(record,),
            requested_input_state=previous.after,
            last_step_error=step_error,
            rollback_state=step_error.rollback_state,
            dt_s=dt_s,
            boundary=boundary,
            seed=seed,
        )
    return step, (record,), "cross_radius_seed_only_exact_root"


def _attempt_classification(attempts: Sequence[AttemptRecord]) -> str:
    if any(attempt.retained_water_applicability_guard_failed for attempt in attempts):
        return "wet_retained_water_applicability_guard_exceeded"
    if any(attempt.external_mobile_liquid_active_set_required for attempt in attempts):
        return "external_mobile_liquid_active_set_required"
    if any(attempt.external_free_water_active_set_required for attempt in attempts):
        return "external_free_water_active_set_required"
    return "same_cell_nonlinear_failure"


def _adapt_accepted_face_event_macrostep(
    bracket: (
        event_orchestrator.AcceptedHistoryFaceEventBracket
        | cc.MasterFaceEventContinuation
    ),
    macro: event_orchestrator.CutEventMacrostep,
) -> AcceptedFaceEventMacrostep:
    """Expose an accepted event macro through the campaign ledger protocol."""

    event = macro.face_event_step
    departure = macro.face_departure_step
    if event is None or departure is None:
        raise ValueError("dynamic trajectory needs both arrival and departure")
    before_cap = wrc.certify_exact_graph(
        macro.before.transport.wet_retained_water_loadings,
        macro.before.transport.effective_wet_retained_water_capacity_duals_over_rt,
        macro.before.transport.config.wet.luikov,
    )
    after_cap = wrc.certify_exact_graph(
        macro.after.transport.wet_retained_water_loadings,
        macro.after.transport.effective_wet_retained_water_capacity_duals_over_rt,
        macro.after.transport.config.wet.luikov,
    )
    if before_cap.active_piece_count or after_cap.active_piece_count:
        raise RuntimeError(
            "rank-changing face trajectory is not yet licensed to carry an active "
            "wet retained-cap normal-cone piece"
        )
    ledger = AdaptiveSameCellMacroLedger(
        dt_s=macro.ledger.requested_dt_s,
        nonlinear_evaluations=(macro.ledger.total_attempted_nonlinear_evaluations),
        rejected_trial_evaluations=sum(
            (
                event.ledger.rejected_trial_evaluations,
                departure.ledger.rejected_trial_evaluations,
            )
        ),
        nonlinear_message=(
            "accepted-history bracket -> exact face arrival -> exact adjacent-cell "
            "departure; atomic raw macro ledger"
        ),
        maximum_scaled_residual=macro.ledger.maximum_stage_scaled_residual,
        condition_proxy=max(
            event.ledger.condition_proxy,
            departure.ledger.condition_proxy,
        ),
        minimum_fractional_distance_to_bound=min(
            event.ledger.minimum_fractional_distance_to_bound,
            departure.ledger.minimum_fractional_distance_to_bound,
        ),
        minimum_front_z_distance_to_chart_boundary=(
            departure.ledger.minimum_front_z_distance_to_chart_boundary
        ),
        water_step_residual_mol=macro.ledger.water_step_residual_mol,
        hexane_step_residual_mol=macro.ledger.hexane_step_residual_mol,
        energy_step_residual_j=macro.ledger.energy_step_residual_j,
        normalized_water_step_residual=(macro.ledger.normalized_water_step_residual),
        normalized_hexane_step_residual=(macro.ledger.normalized_hexane_step_residual),
        normalized_energy_step_residual=(macro.ledger.normalized_energy_step_residual),
        water_cumulative_residual_mol=macro.ledger.water_cumulative_residual_mol,
        hexane_cumulative_residual_mol=(macro.ledger.hexane_cumulative_residual_mol),
        energy_cumulative_residual_j=macro.ledger.energy_cumulative_residual_j,
        normalized_water_cumulative_residual=(macro.ledger.normalized_water_cumulative_residual),
        normalized_hexane_cumulative_residual=(macro.ledger.normalized_hexane_cumulative_residual),
        normalized_energy_cumulative_residual=(macro.ledger.normalized_energy_cumulative_residual),
        boundary_water_out_mol=macro.ledger.boundary_water_out_mol,
        boundary_hexane_out_mol=macro.ledger.boundary_hexane_out_mol,
        boundary_energy_out_j=macro.ledger.boundary_energy_out_j,
        wet_retained_cap_active_piece_indices=after_cap.active_piece_indices,
        wet_retained_cap_active_piece_count=after_cap.active_piece_count,
        wet_retained_cap_maximum_loading=max(
            before_cap.maximum_loading,
            after_cap.maximum_loading,
        ),
        wet_retained_cap_maximum_capacity_dual_over_rt=max(
            before_cap.maximum_capacity_dual_over_rt,
            after_cap.maximum_capacity_dual_over_rt,
        ),
        wet_retained_cap_maximum_complementarity_product=max(
            before_cap.maximum_complementarity_product,
            after_cap.maximum_complementarity_product,
        ),
        wet_retained_cap_exact_graph_without_tolerance=(
            before_cap.exact_graph_without_tolerance and after_cap.exact_graph_without_tolerance
        ),
        # The face route changes rank, so same-piece transition labels are not
        # manufactured.  Active pieces are rejected above until that event
        # graph has its own conservative primal/dual handoff certificate.
        wet_retained_cap_piece_transition_history=(),
        wet_retained_cap_entry_piece_count=0,
        wet_retained_cap_continuation_piece_count=0,
        wet_retained_cap_exit_piece_count=0,
    )
    return AcceptedFaceEventMacrostep(macro, bracket, ledger)


def _event_failure_stage_from_exception_chain(
    error: event_orchestrator.CutEventMacrostepError,
    same_cell_error: ci.CutIntegratorStepError,
) -> object | None:
    """Find diagnostics for the rejected event stage without using its candidate."""

    cause = error.__cause__
    visited: set[int] = set()
    while cause is not None and id(cause) not in visited:
        visited.add(id(cause))
        if cause is not same_cell_error and hasattr(cause, "nonlinear_evaluations"):
            return cause
        cause = cause.__cause__
    return None


def _uncertified_face_route_recovery_error(
    previous: AcceptedSameCellStep,
    bracket: event_orchestrator.AcceptedHistoryFaceEventBracket,
    macrostep_error: event_orchestrator.CutEventMacrostepError,
    boundary: ct.PoreBoundary,
    *,
    wall_time_s: float,
) -> AcceptedHistoryFaceRouteRecoveryError:
    """Quarantine rejected event work and retain same-cell recovery authority."""

    original = previous.after
    _require_exact_rollback_identity(
        original,
        macrostep_error.rollback_state,
        context="accepted-history face macrostep",
    )
    same_cell_error = macrostep_error.same_cell_error
    if same_cell_error is None:
        raise macrostep_error
    _require_exact_rollback_identity(
        original,
        same_cell_error.rollback_state,
        context="accepted-history same-cell first refusal",
    )
    stage = _event_failure_stage_from_exception_chain(
        macrostep_error,
        same_cell_error,
    )
    stage_evaluations = int(getattr(stage, "nonlinear_evaluations", 0))
    stage_rejected = int(getattr(stage, "rejected_trial_evaluations", 0))
    stage_maximum_residual = float(getattr(stage, "maximum_scaled_residual", math.inf))
    stage_condition = float(getattr(stage, "condition_proxy", math.inf))
    same_cell_record = AttemptRecord(
        role="accepted_history_bracket_same_cell_first_refusal",
        predictor_fraction=1.0,
        wall_time_s=0.0,
        accepted=False,
        nonlinear_evaluations=same_cell_error.nonlinear_evaluations,
        maximum_scaled_residual=same_cell_error.maximum_scaled_residual,
        condition_proxy=same_cell_error.condition_proxy,
        rollback_identity_preserved=True,
        error=str(same_cell_error),
        face_event_trial_topology_excursion_observed=(same_cell_error.event_restart_required),
        scaled_residual_block_maxima=_scaled_residual_block_maxima(
            original.transport.layout,
            same_cell_error.scaled_residuals,
        ),
        macro_branch_status="same_cell_first_refusal_rejected",
        provisional_work_budget=same_cell_error.provisional_work_budget,
        provisional_work_budget_exhausted=(same_cell_error.provisional_work_budget_exhausted),
        attempted_dt_s=bracket.requested_dt_s,
        pre_adaptive_macro_branch_status="same_cell_first_refusal_rejected",
    )
    face_record = AttemptRecord(
        role="accepted_history_uncertified_face_dispatch",
        predictor_fraction=1.0,
        wall_time_s=wall_time_s,
        accepted=False,
        nonlinear_evaluations=stage_evaluations,
        maximum_scaled_residual=stage_maximum_residual,
        condition_proxy=stage_condition,
        rollback_identity_preserved=True,
        error=str(macrostep_error),
        macro_branch_status="uncertified_face_dispatch_rejected",
        attempted_dt_s=bracket.requested_dt_s,
        pre_adaptive_macro_branch_status="uncertified_face_dispatch_rejected",
        event_failure_code=macrostep_error.code.value,
        event_failure_stage_type=(None if stage is None else type(stage).__name__),
        event_failure_nonlinear_evaluations=stage_evaluations,
        event_failure_rejected_trial_evaluations=stage_rejected,
        event_failure_maximum_scaled_residual=stage_maximum_residual,
        event_failure_condition_proxy=stage_condition,
    )
    failure_audit = _failure_jacobian_audit(
        original,
        bracket.requested_dt_s,
        boundary,
        bracket.same_cell_seed,
        same_cell_error,
    )
    return AcceptedHistoryFaceRouteRecoveryError(
        (
            "accepted-history tangent route remained uncertified after exact "
            "atomic rollback; only the original same-cell failure may enter "
            "independent adaptive recovery"
        ),
        rollback_state=original,
        attempts=(same_cell_record, face_record),
        last_step_error=same_cell_error,
        failure_jacobian_audit=failure_audit,
        bracket=bracket,
        macrostep_error=macrostep_error,
        original_same_cell_seed=bracket.same_cell_seed,
    )


def _advance_bracketed_face_from_accepted_history(
    previous: AcceptedSameCellStep,
    dt_s: float,
    boundary: ct.PoreBoundary,
) -> (
    tuple[
        ci.CutIntegratorStep | AcceptedFaceEventMacrostep,
        tuple[AttemptRecord, ...],
        str,
    ]
    | None
):
    """Use accepted history only to license an atomic first-refusal route."""

    source = _accepted_face_history_source(previous)
    bracket = event_orchestrator.bracket_inner_face_from_accepted_history(
        source,
        dt_s,
        provenance=(
            "Gate-1g accepted trajectory/departure history; deterministic "
            "constant-front-tangent routing input; no rejected-root authority"
        ),
    )
    if bracket is None:
        return None
    started = perf_counter()
    try:
        macro = event_orchestrator.advance_bracketed_interior_face_macrostep(
            bracket,
            boundary,
        )
    except event_orchestrator.CutEventMacrostepError as exc:
        if exc.code is (
            event_orchestrator.CutMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED
        ):
            raise
        raise _uncertified_face_route_recovery_error(
            previous,
            bracket,
            exc,
            boundary,
            wall_time_s=perf_counter() - started,
        ) from exc
    elapsed = perf_counter() - started

    if macro.before is not previous.after:
        raise RuntimeError("accepted-history macrostep did not preserve the exact input object")
    if macro.ledger.branch is event_orchestrator.CutMacrostepBranch.SAME_CELL:
        same_cell = macro.same_cell_step
        if (
            same_cell is None
            or macro.face_event_step is not None
            or macro.face_departure_step is not None
            or same_cell.before is not previous.after
            or same_cell.after is not macro.after
        ):
            raise RuntimeError("accepted-history first-refusal branch lost its same-cell payload")
        record = AttemptRecord(
            role="accepted_history_bracket_same_cell_first_refusal",
            predictor_fraction=1.0,
            wall_time_s=elapsed,
            accepted=True,
            nonlinear_evaluations=(macro.ledger.total_attempted_nonlinear_evaluations),
            maximum_scaled_residual=macro.ledger.maximum_stage_scaled_residual,
            condition_proxy=same_cell.ledger.condition_proxy,
            rollback_identity_preserved=True,
            macro_branch_status=(event_orchestrator.CutMacrostepBranch.SAME_CELL.value),
            attempted_dt_s=dt_s,
        )
        return (
            same_cell,
            (record,),
            "accepted_history_bracket_same_cell_first_refusal",
        )

    if macro.ledger.branch is not (
        event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE
    ):
        raise RuntimeError("accepted-history dynamic route cannot commit an exact-face endpoint")
    adapted = _adapt_accepted_face_event_macrostep(bracket, macro)
    record = AttemptRecord(
        role="accepted_history_atomic_interior_face_macrostep",
        predictor_fraction=1.0,
        wall_time_s=elapsed,
        accepted=True,
        nonlinear_evaluations=(macro.ledger.total_attempted_nonlinear_evaluations),
        maximum_scaled_residual=macro.ledger.maximum_stage_scaled_residual,
        condition_proxy=adapted.ledger.condition_proxy,
        rollback_identity_preserved=macro.before is previous.after,
        macro_branch_status=(
            event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE.value
        ),
        attempted_dt_s=dt_s,
    )
    return adapted, (record,), "accepted_history_atomic_face_event"


def _compose_adaptive_same_cell_macrostep(
    before_step: AcceptedSameCellStep,
    boundary: ct.DirichletPoreBoundary,
    requested_dt_s: float,
    substeps: Sequence[ci.CutIntegratorStep],
    *,
    binary_depth: int,
    attempted_binary_depths: tuple[int, ...],
    predictor_backtrack_maximum_binary_exponent: int = (
        ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
    ),
    initial_failure_audit: dict | None = None,
) -> AdaptiveSameCellMacrostep:
    """Compose raw endpoint/transfer ledgers after every provisional leaf passes."""

    leaves = tuple(substeps)
    if not leaves:
        raise ValueError("adaptive macrostep needs at least one accepted leaf")
    original = before_step.after
    if leaves[0].before is not original:
        raise RuntimeError("adaptive branch did not start from the exact old state")
    for left, right in zip(leaves, leaves[1:]):
        if right.before is not left.after:
            raise RuntimeError("adaptive leaves do not form one immutable state chain")
    if any(leaf.boundary is not boundary for leaf in leaves):
        raise RuntimeError("adaptive leaves did not use the exact post-jump boundary")
    final = leaves[-1].after
    expected_time = original.transport.time_s + requested_dt_s
    if not math.isclose(
        final.transport.time_s,
        expected_time,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise RuntimeError("adaptive leaves do not close the requested macro interval")
    if any(
        leaf.after.transport.geometry.master_grid.R != original.transport.geometry.master_grid.R
        for leaf in leaves
    ):
        raise RuntimeError("adaptive branch changed the target particle radius")

    leaf_cap_certificates = []
    for leaf in leaves:
        before_cap = wrc.certify_exact_graph(
            leaf.before.transport.wet_retained_water_loadings,
            leaf.before.transport.effective_wet_retained_water_capacity_duals_over_rt,
            leaf.before.transport.config.wet.luikov,
        )
        after_cap = wrc.certify_exact_graph(
            leaf.after.transport.wet_retained_water_loadings,
            leaf.after.transport.effective_wet_retained_water_capacity_duals_over_rt,
            leaf.after.transport.config.wet.luikov,
        )
        transition = wrc.certify_exact_transition(before_cap, after_cap)
        if isinstance(leaf, ci.CutIntegratorStep):
            recorded = leaf.ledger
            exact_match = all(
                (
                    recorded.wet_retained_cap_active_piece_indices
                    == after_cap.active_piece_indices,
                    recorded.wet_retained_cap_active_piece_count == after_cap.active_piece_count,
                    recorded.wet_retained_cap_maximum_loading
                    == max(before_cap.maximum_loading, after_cap.maximum_loading),
                    recorded.wet_retained_cap_maximum_capacity_dual_over_rt
                    == max(
                        before_cap.maximum_capacity_dual_over_rt,
                        after_cap.maximum_capacity_dual_over_rt,
                    ),
                    recorded.wet_retained_cap_maximum_complementarity_product == 0.0,
                    recorded.wet_retained_cap_piece_transitions == transition.piece_transitions,
                )
            )
            if not exact_match:
                raise RuntimeError(
                    "accepted leaf retained-cap ledger disagrees with its exact state"
                )
        leaf_cap_certificates.append((before_cap, after_cap, transition))

    old_inventory = ci.inventory_snapshot(original.transport)
    new_inventory = ci.inventory_snapshot(final.transport)
    boundary_water_out = math.fsum(leaf.ledger.boundary_water_out_mol for leaf in leaves)
    boundary_hexane_out = math.fsum(leaf.ledger.boundary_hexane_out_mol for leaf in leaves)
    boundary_energy_out = math.fsum(leaf.ledger.boundary_energy_out_j for leaf in leaves)
    water_step_residual = math.fsum(
        (
            new_inventory.total_water_mol,
            -old_inventory.total_water_mol,
            boundary_water_out,
        )
    )
    hexane_step_residual = math.fsum(
        (
            new_inventory.total_hexane_mol,
            -old_inventory.total_hexane_mol,
            boundary_hexane_out,
        )
    )
    energy_step_residual = math.fsum(
        (
            new_inventory.total_energy_j,
            -old_inventory.total_energy_j,
            boundary_energy_out,
        )
    )
    water_step_scale = max(
        abs(old_inventory.total_water_mol),
        abs(new_inventory.total_water_mol),
        abs(boundary_water_out),
        math.fsum(abs(leaf.ledger.boundary_water_out_mol) for leaf in leaves),
        1.0e-300,
    )
    hexane_step_scale = max(
        abs(old_inventory.total_hexane_mol),
        abs(new_inventory.total_hexane_mol),
        abs(boundary_hexane_out),
        math.fsum(abs(leaf.ledger.boundary_hexane_out_mol) for leaf in leaves),
        1.0e-300,
    )
    energy_step_scale = max(
        ci.capacity_energy_scale(original.transport),
        ci.capacity_energy_scale(final.transport),
        abs(boundary_energy_out),
        math.fsum(abs(leaf.ledger.boundary_energy_out_j) for leaf in leaves),
        1.0e-300,
    )

    reference = original.reference_inventory
    cumulative_water_out = final.corrected_cumulative_boundary_water_out_mol
    cumulative_hexane_out = final.corrected_cumulative_boundary_hexane_out_mol
    cumulative_energy_out = final.corrected_cumulative_boundary_energy_out_j
    water_cumulative_residual = math.fsum(
        (
            new_inventory.total_water_mol,
            -reference.total_water_mol,
            cumulative_water_out,
        )
    )
    hexane_cumulative_residual = math.fsum(
        (
            new_inventory.total_hexane_mol,
            -reference.total_hexane_mol,
            cumulative_hexane_out,
        )
    )
    energy_cumulative_residual = math.fsum(
        (
            new_inventory.total_energy_j,
            -reference.total_energy_j,
            cumulative_energy_out,
        )
    )
    water_cumulative_scale = max(
        abs(reference.total_water_mol),
        abs(new_inventory.total_water_mol),
        abs(cumulative_water_out),
        final.cumulative_absolute_water_transfer_mol,
        1.0e-300,
    )
    hexane_cumulative_scale = max(
        abs(reference.total_hexane_mol),
        abs(new_inventory.total_hexane_mol),
        abs(cumulative_hexane_out),
        final.cumulative_absolute_hexane_transfer_mol,
        1.0e-300,
    )
    energy_cumulative_scale = max(
        original.reference_capacity_energy_scale_j,
        ci.capacity_energy_scale(final.transport),
        abs(cumulative_energy_out),
        final.cumulative_absolute_energy_transfer_j,
        1.0e-300,
    )
    ledger = AdaptiveSameCellMacroLedger(
        dt_s=requested_dt_s,
        nonlinear_evaluations=sum(leaf.ledger.nonlinear_evaluations for leaf in leaves),
        rejected_trial_evaluations=sum(leaf.ledger.rejected_trial_evaluations for leaf in leaves),
        nonlinear_message=(
            f"certified exact binary subdivision depth {binary_depth}; "
            f"{len(leaves)} provisional leaves committed atomically"
        ),
        maximum_scaled_residual=max(leaf.ledger.maximum_scaled_residual for leaf in leaves),
        condition_proxy=max(leaf.ledger.condition_proxy for leaf in leaves),
        minimum_fractional_distance_to_bound=min(
            leaf.ledger.minimum_fractional_distance_to_bound for leaf in leaves
        ),
        minimum_front_z_distance_to_chart_boundary=min(
            leaf.ledger.minimum_front_z_distance_to_chart_boundary for leaf in leaves
        ),
        water_step_residual_mol=water_step_residual,
        hexane_step_residual_mol=hexane_step_residual,
        energy_step_residual_j=energy_step_residual,
        normalized_water_step_residual=(abs(water_step_residual) / water_step_scale),
        normalized_hexane_step_residual=(abs(hexane_step_residual) / hexane_step_scale),
        normalized_energy_step_residual=(abs(energy_step_residual) / energy_step_scale),
        water_cumulative_residual_mol=water_cumulative_residual,
        hexane_cumulative_residual_mol=hexane_cumulative_residual,
        energy_cumulative_residual_j=energy_cumulative_residual,
        normalized_water_cumulative_residual=(
            abs(water_cumulative_residual) / water_cumulative_scale
        ),
        normalized_hexane_cumulative_residual=(
            abs(hexane_cumulative_residual) / hexane_cumulative_scale
        ),
        normalized_energy_cumulative_residual=(
            abs(energy_cumulative_residual) / energy_cumulative_scale
        ),
        boundary_water_out_mol=boundary_water_out,
        boundary_hexane_out_mol=boundary_hexane_out,
        boundary_energy_out_j=boundary_energy_out,
        wet_retained_cap_active_piece_indices=(leaf_cap_certificates[-1][1].active_piece_indices),
        wet_retained_cap_active_piece_count=(leaf_cap_certificates[-1][1].active_piece_count),
        wet_retained_cap_maximum_loading=max(
            certificate.maximum_loading
            for before_cap, after_cap, _ in leaf_cap_certificates
            for certificate in (before_cap, after_cap)
        ),
        wet_retained_cap_maximum_capacity_dual_over_rt=max(
            certificate.maximum_capacity_dual_over_rt
            for before_cap, after_cap, _ in leaf_cap_certificates
            for certificate in (before_cap, after_cap)
        ),
        wet_retained_cap_maximum_complementarity_product=max(
            certificate.maximum_complementarity_product
            for before_cap, after_cap, _ in leaf_cap_certificates
            for certificate in (before_cap, after_cap)
        ),
        wet_retained_cap_exact_graph_without_tolerance=True,
        wet_retained_cap_piece_transition_history=tuple(
            transition.piece_transitions for _, _, transition in leaf_cap_certificates
        ),
        wet_retained_cap_entry_piece_count=sum(
            transition.entry_piece_count for _, _, transition in leaf_cap_certificates
        ),
        wet_retained_cap_continuation_piece_count=sum(
            transition.continuation_piece_count for _, _, transition in leaf_cap_certificates
        ),
        wet_retained_cap_exit_piece_count=sum(
            transition.exit_piece_count for _, _, transition in leaf_cap_certificates
        ),
    )
    if (
        max(
            ledger.maximum_step_ledger_residual,
            ledger.maximum_cumulative_ledger_residual,
        )
        > CONSERVATION_LIMIT
    ):
        raise RuntimeError("adaptive macrostep raw component/energy ledger exceeded 1e-10")
    for state in (original, *(leaf.after for leaf in leaves)):
        ci.certify_explicit_wet_retained_water_applicability(
            state.controls,
            state.transport.wet_retained_water_loadings,
            state.transport.effective_wet_retained_water_capacity_duals_over_rt,
            state.transport.config.wet.luikov,
        )
    return AdaptiveSameCellMacrostep(
        before=original,
        after=final,
        assembly=leaves[-1].assembly,
        boundary=boundary,
        seed=leaves[0].seed,
        ledger=ledger,
        substeps=leaves,
        binary_depth=binary_depth,
        attempted_binary_depths=attempted_binary_depths,
        predictor_backtrack_maximum_binary_exponent=(predictor_backtrack_maximum_binary_exponent),
        initial_failure_audit=initial_failure_audit,
    )


def _provisional_accepted_leaf_summary(
    leaf,
    *,
    depth: int,
    path: str,
    starts_at_reversal: bool,
    retained_water_guard_loading: float | None,
) -> ProvisionalAcceptedLeafSummary | None:
    """Copy release-audit scalars without retaining a provisional state."""

    # Unit tests use deliberately tiny synthetic leaf sentinels to exercise
    # recursion.  Production leaves are always the concrete immutable step.
    if not isinstance(leaf, ci.CutIntegratorStep):
        return None
    transport = leaf.after.transport
    geometry = transport.geometry
    layout = transport.layout
    cut_index = geometry.cut_cell_index
    inner_face_z = (geometry.master_grid.faces[cut_index] / geometry.master_grid.R) ** 3
    outer_face_z = (geometry.master_grid.faces[cut_index + 1] / geometry.master_grid.R) ** 3
    front_z = geometry.front.z
    loadings = transport.wet_retained_water_loadings
    maximum_position = max(range(len(loadings)), key=loadings.__getitem__)
    maximum_loading = loadings[maximum_position]
    chart_upper = leaf.after.controls.wet_water_bounds[1]
    surface = leaf.assembly.dry_face_fluxes[-1]
    film = leaf.assembly.surface_film_audit
    amendment_18_19_evidence = _amendment_18_19_execution_occurrence(
        leaf,
        execution_classification=(
            "noncommitted_provisional_atomic_macro_rollback"
        ),
        execution_context={
            "execution_partition_id": "provisional_atomic_macro_rollback",
            "execution_role": "adaptive_same_cell_provisional_leaf",
            "adaptive_depth": depth,
            "adaptive_path": path,
            "starts_at_reversal": starts_at_reversal,
            "time_before_s": leaf.before.transport.time_s,
            "time_after_s": transport.time_s,
            "dt_s": leaf.ledger.dt_s,
        },
    )
    return ProvisionalAcceptedLeafSummary(
        adaptive_depth=depth,
        adaptive_path=path,
        starts_at_reversal=starts_at_reversal,
        time_before_s=leaf.before.transport.time_s,
        time_after_s=transport.time_s,
        dt_s=leaf.ledger.dt_s,
        front_z=front_z,
        cut_cell_index=cut_index,
        inner_face_z=inner_face_z,
        outer_face_z=outer_face_z,
        nearest_master_face_distance_z=min(
            front_z - inner_face_z,
            outer_face_z - front_z,
        ),
        maximum_wet_retained_water_loading=maximum_loading,
        maximum_loading_wet_piece_position=maximum_position,
        maximum_loading_master_cell_index=(layout.wet_cell_indices[maximum_position]),
        retained_water_chart_upper_loading=chart_upper,
        distance_to_retained_water_chart_upper=chart_upper - maximum_loading,
        retained_water_guard_loading=retained_water_guard_loading,
        retained_water_guard_margin=(
            None
            if retained_water_guard_loading is None
            else retained_water_guard_loading - maximum_loading
        ),
        wet_retained_cap_active_piece_indices=(leaf.ledger.wet_retained_cap_active_piece_indices),
        wet_retained_cap_active_piece_count=(leaf.ledger.wet_retained_cap_active_piece_count),
        wet_retained_cap_maximum_capacity_dual_over_rt=(
            leaf.ledger.wet_retained_cap_maximum_capacity_dual_over_rt
        ),
        wet_retained_cap_maximum_complementarity_product=(
            leaf.ledger.wet_retained_cap_maximum_complementarity_product
        ),
        wet_retained_cap_exact_graph_without_tolerance=(
            leaf.ledger.wet_retained_cap_exact_graph_without_tolerance
        ),
        wet_retained_cap_piece_transitions=(leaf.ledger.wet_retained_cap_piece_transitions),
        wet_retained_cap_entry_piece_count=(leaf.ledger.wet_retained_cap_entry_piece_count),
        wet_retained_cap_continuation_piece_count=(
            leaf.ledger.wet_retained_cap_continuation_piece_count
        ),
        wet_retained_cap_exit_piece_count=(leaf.ledger.wet_retained_cap_exit_piece_count),
        surface_conserved_water_flux_mol_m2_s=(surface.component.conserved_water_flux_mol_m2_s),
        surface_conserved_hexane_flux_mol_m2_s=(surface.component.conserved_hexane_flux_mol_m2_s),
        surface_independent_hexane_flux_mol_m2_s=(-surface.independent_water_flux_mol_m2_s),
        surface_total_stefan_flux_mol_m2_s=(surface.component.total_stefan_flux_mol_m2_s),
        maximum_scaled_residual=leaf.ledger.maximum_scaled_residual,
        scaled_residual_block_maxima=_scaled_residual_block_maxima(
            layout,
            leaf.ledger.scaled_residuals,
        ),
        condition_proxy=leaf.ledger.condition_proxy,
        nonlinear_evaluations=leaf.ledger.nonlinear_evaluations,
        water_step_residual_mol=leaf.ledger.water_step_residual_mol,
        hexane_step_residual_mol=leaf.ledger.hexane_step_residual_mol,
        energy_step_residual_j=leaf.ledger.energy_step_residual_j,
        maximum_step_ledger_residual=leaf.ledger.maximum_step_ledger_residual,
        maximum_cumulative_ledger_residual=(leaf.ledger.maximum_cumulative_ledger_residual),
        surface_fast_mass_departure=(
            None if film is None else film.mass_departure.absolute_surface_departure
        ),
        surface_fast_mass_guard_passed=(None if film is None else film.fast_mass_guard_passed),
        amendment_18_19_diagnostic_execution_evidence=(amendment_18_19_evidence),
        amendment_18_19_diagnostic_execution_evidence_contract_passed=True,
    )


def _advance_with_adaptive_same_cell_macrostep(
    previous: AcceptedSameCellStep,
    requested_dt_s: float,
    previous_boundary: ct.DirichletPoreBoundary,
    boundary: ct.DirichletPoreBoundary,
    *,
    is_reversal: bool,
    initial_attempts: Sequence[AttemptRecord] = (),
    initial_step_error: ci.CutIntegratorStepError | None = None,
    initial_failure_audit: dict | None = None,
    maximum_binary_depth: int = ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH,
    retained_water_guard_loading: float | None = None,
    predictor_backtrack_maximum_binary_exponent: int = (
        ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
    ),
) -> tuple[AdaptiveSameCellMacrostep, tuple[AttemptRecord, ...], str]:
    """Retry one rejected interval by bisecting only each failed BE leaf.

    The caller must first have tried the requested macro interval and supplied
    its rejected-attempt audit.  Children are then visited deterministically
    from left to right.  An accepted child is retained as an immutable,
    provisional state while only a failed sibling is subdivided.  Nothing is
    returned (and therefore nothing can enter the trajectory) until all leaves
    form the requested interval and the independently recomposed raw ledger
    passes the unchanged ``1e-10`` contract.
    """

    if (
        isinstance(maximum_binary_depth, bool)
        or not isinstance(maximum_binary_depth, int)
        or not 1 <= maximum_binary_depth <= ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH
    ):
        raise ValueError(
            "adaptive binary-depth cap must be an integer in "
            f"[1, {ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH}]"
        )
    if not math.isfinite(requested_dt_s) or requested_dt_s <= 0.0:
        raise ValueError("adaptive requested dt must be positive and finite")
    if (
        isinstance(predictor_backtrack_maximum_binary_exponent, bool)
        or not isinstance(predictor_backtrack_maximum_binary_exponent, int)
        or predictor_backtrack_maximum_binary_exponent < 1
        or predictor_backtrack_maximum_binary_exponent > PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
    ):
        raise ValueError(
            "adaptive predictor budget must be an integer in "
            f"[1, {PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT}]"
        )
    original = previous.after
    if retained_water_guard_loading is not None:
        if (
            not math.isfinite(retained_water_guard_loading)
            or not original.controls.wet_water_bounds[0]
            < retained_water_guard_loading
            < original.controls.wet_water_bounds[1]
        ):
            raise ValueError("retained-water evidence guard must lie strictly inside the chart")
    recorded_attempts = [
        replace(
            attempt,
            macro_branch_status="initial_requested_failure",
            attempted_dt_s=requested_dt_s,
            adaptive_depth=0,
            adaptive_path="root",
        )
        for attempt in initial_attempts
    ]
    if not recorded_attempts:
        raise ValueError("adaptive subdivision requires an audited rejected requested leaf")
    if any(attempt.accepted for attempt in recorded_attempts):
        raise ValueError("adaptive subdivision cannot replace an accepted requested leaf")
    _require_failed_attempt_rollbacks(
        recorded_attempts,
        context="initial adaptive same-cell failure",
    )
    if initial_step_error is not None:
        _require_exact_rollback_identity(
            original,
            initial_step_error.rollback_state,
            context="initial adaptive same-cell step error",
        )
    initial_classification = _attempt_classification(recorded_attempts)
    effective_maximum_binary_depth = _effective_adaptive_maximum_binary_depth(
        requested_dt_s,
        maximum_binary_depth_cap=maximum_binary_depth,
    )
    minimum_leaf_dt_s = (
        math.ldexp(requested_dt_s, -effective_maximum_binary_depth)
        if effective_maximum_binary_depth
        else requested_dt_s
    )
    if initial_classification != "same_cell_nonlinear_failure":
        raise AdaptiveSameCellMacrostepError(
            "adaptive same-cell subdivision is forbidden for a non-same-cell "
            f"transition: {initial_classification}",
            requested_input_state=original,
            rollback_state=original,
            requested_dt_s=requested_dt_s,
            minimum_leaf_dt_s=minimum_leaf_dt_s,
            attempted_binary_depths=(),
            attempts=tuple(recorded_attempts),
            classification=initial_classification,
            last_step_error=initial_step_error,
            initial_failure_audit=initial_failure_audit,
            provisional_accepted_leaf_summaries=(),
            predictor_backtrack_maximum_binary_exponent=(
                predictor_backtrack_maximum_binary_exponent
            ),
        )
    if effective_maximum_binary_depth == 0:
        raise AdaptiveSameCellMacrostepError(
            "same-cell failed-leaf recovery has no licensed child at or above "
            f"the absolute minimum dt "
            f"{ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S:g} s",
            requested_input_state=original,
            rollback_state=original,
            requested_dt_s=requested_dt_s,
            minimum_leaf_dt_s=minimum_leaf_dt_s,
            attempted_binary_depths=(),
            attempts=tuple(recorded_attempts),
            classification="same_cell_absolute_minimum_leaf_floor_prevents_subdivision",
            last_step_error=initial_step_error,
            initial_failure_audit=initial_failure_audit,
            provisional_accepted_leaf_summaries=(),
            predictor_backtrack_maximum_binary_exponent=(
                predictor_backtrack_maximum_binary_exponent
            ),
        )

    attempted_depths: list[int] = []
    provisional_accepted_leaf_summaries: list[ProvisionalAcceptedLeafSummary] = []
    last_step_error = initial_step_error
    last_auditable_failure_jacobian_audit: dict | None = None
    adaptive_attempt_start = len(recorded_attempts)

    def note_depth(depth: int) -> None:
        if depth not in attempted_depths:
            attempted_depths.append(depth)

    def solve_leaf(
        current: AcceptedSameCellStep,
        leaf_dt_s: float,
        *,
        depth: int,
        path: str,
        starts_at_reversal: bool,
    ) -> list[ci.CutIntegratorStep]:
        """Return one certified leaf chain or propagate its exact failure."""

        nonlocal last_step_error, last_auditable_failure_jacobian_audit
        note_depth(depth)
        try:
            leaf, attempts, _ = _advance_with_recovery(
                current,
                leaf_dt_s,
                previous_boundary if starts_at_reversal else boundary,
                boundary,
                is_reversal=starts_at_reversal,
                try_current_state_anchor=True,
                predictor_backtrack_maximum_binary_exponent=(
                    predictor_backtrack_maximum_binary_exponent
                ),
                provisional_maximum_function_evaluations=(
                    DISCONTINUITY_LEAF_MAXIMUM_OPTIMIZER_FUNCTION_EVALUATIONS
                ),
                maximum_nonlinear_seed_attempts=(
                    DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
                    if starts_at_reversal
                    else DISCONTINUITY_NONREVERSAL_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
                ),
            )
        except SameCellRecoveryError as exc:
            context = f"adaptive depth {depth} path {path}"
            _require_exact_rollback_identity(
                current.after,
                exc.rollback_state,
                context=f"{context} failure",
            )
            _require_failed_attempt_rollbacks(
                exc.attempts,
                context=f"{context} attempts",
            )
            recorded_attempts.extend(
                replace(
                    attempt,
                    role=f"adaptive_d{depth}_{path}::{attempt.role}",
                    macro_branch_status="failed_leaf_bisected",
                    attempted_dt_s=leaf_dt_s,
                    adaptive_depth=depth,
                    adaptive_path=path,
                )
                for attempt in exc.attempts
            )
            last_step_error = exc.last_step_error
            if exc.failure_jacobian_audit is not None:
                last_auditable_failure_jacobian_audit = exc.failure_jacobian_audit
            classification = _attempt_classification(exc.attempts)
            if classification != "same_cell_nonlinear_failure":
                raise
            if depth >= effective_maximum_binary_depth:
                raise

            half_dt_s = 0.5 * leaf_dt_s
            left = solve_leaf(
                current,
                half_dt_s,
                depth=depth + 1,
                path=f"{path}L",
                starts_at_reversal=starts_at_reversal,
            )
            left_endpoint: AcceptedSameCellStep = left[-1]
            right = solve_leaf(
                left_endpoint,
                half_dt_s,
                depth=depth + 1,
                path=f"{path}R",
                starts_at_reversal=False,
            )
            return [*left, *right]

        _require_failed_attempt_rollbacks(
            attempts,
            context=f"adaptive depth {depth} path {path} attempts",
        )
        if leaf.before is not current.after:
            raise RuntimeError("adaptive leaf acceptance does not reference its exact input")
        summary = _provisional_accepted_leaf_summary(
            leaf,
            depth=depth,
            path=path,
            starts_at_reversal=starts_at_reversal,
            retained_water_guard_loading=retained_water_guard_loading,
        )
        if summary is not None:
            provisional_accepted_leaf_summaries.append(summary)
        recorded_attempts.extend(
            replace(
                attempt,
                role=f"adaptive_d{depth}_{path}::{attempt.role}",
                macro_branch_status="returned",
                attempted_dt_s=leaf_dt_s,
                adaptive_depth=depth,
                adaptive_path=path,
            )
            for attempt in attempts
        )
        return [leaf]

    try:
        half_dt_s = 0.5 * requested_dt_s
        leaves = solve_leaf(
            previous,
            half_dt_s,
            depth=1,
            path="L",
            starts_at_reversal=is_reversal,
        )
        right_start: AcceptedSameCellStep = leaves[-1]
        leaves.extend(
            solve_leaf(
                right_start,
                half_dt_s,
                depth=1,
                path="R",
                starts_at_reversal=False,
            )
        )
    except SameCellRecoveryError as exc:
        classification = _attempt_classification(exc.attempts)
        if classification == "same_cell_nonlinear_failure":
            classification = "same_cell_nonlinear_minimum_dt_exhausted"
            message = (
                "same-cell failed-leaf subdivision exhausted its documented "
                f"minimum dt {minimum_leaf_dt_s:g} s with exact macro rollback"
            )
        else:
            message = (
                "adaptive failed-leaf branch met a licensed non-same-cell "
                f"transition: {classification}"
            )
        recorded_attempts[adaptive_attempt_start:] = [
            replace(attempt, macro_branch_status="discarded")
            for attempt in recorded_attempts[adaptive_attempt_start:]
        ]
        raise AdaptiveSameCellMacrostepError(
            message,
            requested_input_state=original,
            rollback_state=original,
            requested_dt_s=requested_dt_s,
            minimum_leaf_dt_s=minimum_leaf_dt_s,
            attempted_binary_depths=tuple(attempted_depths),
            attempts=tuple(recorded_attempts),
            classification=classification,
            last_step_error=last_step_error,
            initial_failure_audit=initial_failure_audit,
            last_auditable_failure_jacobian_audit=(last_auditable_failure_jacobian_audit),
            provisional_accepted_leaf_summaries=tuple(provisional_accepted_leaf_summaries),
            predictor_backtrack_maximum_binary_exponent=(
                predictor_backtrack_maximum_binary_exponent
            ),
        ) from exc

    maximum_accepted_depth = max(
        int(attempt.role.split("_", 2)[1][1:])
        for attempt in recorded_attempts[adaptive_attempt_start:]
        if attempt.accepted and attempt.macro_branch_status == "returned"
    )
    try:
        macrostep = _compose_adaptive_same_cell_macrostep(
            previous,
            boundary,
            requested_dt_s,
            leaves,
            binary_depth=maximum_accepted_depth,
            attempted_binary_depths=tuple(attempted_depths),
            predictor_backtrack_maximum_binary_exponent=(
                predictor_backtrack_maximum_binary_exponent
            ),
            initial_failure_audit=initial_failure_audit,
        )
    except Exception as exc:
        recorded_attempts[adaptive_attempt_start:] = [
            replace(attempt, macro_branch_status="discarded")
            for attempt in recorded_attempts[adaptive_attempt_start:]
        ]
        raise AdaptiveSameCellMacrostepError(
            f"adaptive macrostep composition rejected: {exc}",
            requested_input_state=original,
            rollback_state=original,
            requested_dt_s=requested_dt_s,
            minimum_leaf_dt_s=minimum_leaf_dt_s,
            attempted_binary_depths=tuple(attempted_depths),
            attempts=tuple(recorded_attempts),
            classification=(
                "wet_retained_water_applicability_guard_exceeded"
                if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
                else "macro_ledger_or_chain_failure"
            ),
            last_step_error=last_step_error,
            initial_failure_audit=initial_failure_audit,
            last_auditable_failure_jacobian_audit=(last_auditable_failure_jacobian_audit),
            provisional_accepted_leaf_summaries=tuple(provisional_accepted_leaf_summaries),
            predictor_backtrack_maximum_binary_exponent=(
                predictor_backtrack_maximum_binary_exponent
            ),
        ) from exc
    return (
        macrostep,
        tuple(recorded_attempts),
        f"adaptive_same_cell_failed_leaf_depth_{maximum_accepted_depth}",
    )


def _recover_after_uncertified_face_route(
    previous: AcceptedSameCellStep,
    requested_dt_s: float,
    previous_boundary: ct.DirichletPoreBoundary,
    boundary: ct.DirichletPoreBoundary,
    *,
    is_reversal: bool,
    failure: AcceptedHistoryFaceRouteRecoveryError,
) -> tuple[AdaptiveSameCellMacrostep, tuple[AttemptRecord, ...], str]:
    """Restart adaptive work from immutable same-cell authority only."""

    if failure.rollback_state is not previous.after:
        raise RuntimeError("uncertified face recovery changed the macro input")
    if failure.last_step_error is not failure.macrostep_error.same_cell_error:
        raise RuntimeError("uncertified face recovery changed the first refusal")
    if failure.original_same_cell_seed is not failure.bracket.same_cell_seed:
        raise RuntimeError("uncertified face recovery changed the original seed")
    if failure.event_candidate_entered_recovery or failure.topology_selected:
        raise RuntimeError("uncertified face work cannot authorize adaptive state")
    return _advance_with_adaptive_same_cell_macrostep(
        previous,
        requested_dt_s,
        previous_boundary,
        boundary,
        is_reversal=is_reversal,
        initial_attempts=failure.attempts,
        initial_step_error=failure.last_step_error,
        initial_failure_audit=failure.failure_jacobian_audit,
    )


@dataclass(frozen=True)
class _AcceptedBalanceIntervalRoot:
    """One accepted BE balance interval and the root closing that interval."""

    role: str
    interval_dt_s: float
    time_s: float
    assembly: object
    temperatures_k: tuple[float, ...]
    front_z: float
    scalar_post_root_polish_audit: ci.ScalarPostRootPolishAudit | None
    residual_scale_evaluations: tuple[tuple[str, object], ...]

    def __post_init__(self) -> None:
        if not self.role.strip():
            raise ValueError("accepted balance root needs a role")
        if not math.isfinite(self.interval_dt_s) or self.interval_dt_s <= 0.0:
            raise ValueError("accepted balance interval must be positive and finite")
        if not math.isfinite(self.time_s) or self.time_s < 0.0:
            raise ValueError("accepted balance root time must be finite and nonnegative")
        if not self.temperatures_k or not all(
            math.isfinite(value) for value in self.temperatures_k
        ):
            raise ValueError("accepted balance root temperatures must be finite")
        if not math.isfinite(self.front_z):
            raise ValueError("accepted balance root front coordinate must be finite")
        if not self.residual_scale_evaluations or not all(
            isinstance(kind, str) and kind.strip()
            for kind, _scales in self.residual_scale_evaluations
        ):
            raise ValueError("accepted balance root needs residual-scale identities")
        kinds = tuple(kind for kind, _scales in self.residual_scale_evaluations)
        if len(set(kinds)) != len(kinds):
            raise ValueError("accepted balance root residual-scale identities must be unique")


SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID = "GT-PS-2-P1E-07"


def _scalar_post_root_polish_residual_evaluations(
    audit: ci.ScalarPostRootPolishAudit,
) -> int:
    """Count residual assemblies uniquely attributable to one scalar polish."""

    if not isinstance(audit, ci.ScalarPostRootPolishAudit):
        raise TypeError("scalar post-root work accounting requires its frozen audit")
    return (
        audit.selection_stencil_evaluations
        + audit.bracket_probe_count
        + audit.scalar_function_calls
        + 1  # explicit full polished-root residual
        + audit.final_condition_stencil_evaluations
    )


def _scalar_post_root_polish_summary(
    roots: Sequence[
        tuple[
            str,
            float,
            float,
            ci.ScalarPostRootPolishAudit | None,
        ]
    ],
) -> dict:
    """Serialize every accepted-root use; an empty audit list means zero uses."""

    serialized = [
        {
            "root_role": role,
            "balance_interval_dt_s": interval_dt_s,
            "time_s": time_s,
            "audit": asdict(audit),
        }
        for role, interval_dt_s, time_s, audit in roots
        if audit is not None
    ]
    audits = [audit for *_, audit in roots if audit is not None]
    return {
        "numerical_amendment_id": SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID,
        "accepted_root_count_considered": len(roots),
        "used_accepted_root_count": len(audits),
        "used_step": bool(audits),
        "maximum_scaled_residual_before": (
            max(audit.maximum_scaled_residual_before for audit in audits) if audits else None
        ),
        "maximum_scaled_residual_after": (
            max(audit.maximum_scaled_residual_after for audit in audits) if audits else None
        ),
        "total_scalar_polish_residual_evaluations": sum(
            _scalar_post_root_polish_residual_evaluations(audit) for audit in audits
        ),
        "audits": serialized,
        "absence_means_zero_uses": True,
        "physics_equations_or_acceptance_tolerances_changed": False,
    }


def _aggregate_scalar_post_root_polish_summaries(
    summaries: Sequence[Mapping[str, object]],
) -> dict:
    """Aggregate explicit step summaries without silently dropping zero-use steps."""

    summary_tuple = tuple(summaries)
    required = {
        "accepted_root_count_considered",
        "used_accepted_root_count",
        "used_step",
        "maximum_scaled_residual_before",
        "maximum_scaled_residual_after",
        "total_scalar_polish_residual_evaluations",
        "absence_means_zero_uses",
    }
    if any(not required <= set(summary) for summary in summary_tuple):
        raise ValueError("scalar post-root aggregate encountered an incomplete summary")
    before = [
        float(summary["maximum_scaled_residual_before"])
        for summary in summary_tuple
        if summary["maximum_scaled_residual_before"] is not None
    ]
    after = [
        float(summary["maximum_scaled_residual_after"])
        for summary in summary_tuple
        if summary["maximum_scaled_residual_after"] is not None
    ]
    return {
        "numerical_amendment_id": SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID,
        "accepted_step_record_count_considered": len(summary_tuple),
        "used_step_record_count": sum(bool(summary["used_step"]) for summary in summary_tuple),
        "accepted_root_count_considered": sum(
            int(summary["accepted_root_count_considered"]) for summary in summary_tuple
        ),
        "used_accepted_root_count": sum(
            int(summary["used_accepted_root_count"]) for summary in summary_tuple
        ),
        "maximum_scaled_residual_before": max(before) if before else None,
        "maximum_scaled_residual_after": max(after) if after else None,
        "total_scalar_polish_residual_evaluations": sum(
            int(summary["total_scalar_polish_residual_evaluations"]) for summary in summary_tuple
        ),
        "all_steps_explicitly_serialized": all(
            summary["absence_means_zero_uses"] is True for summary in summary_tuple
        ),
        "absence_means_zero_uses": True,
        "physics_equations_or_acceptance_tolerances_changed": False,
    }


def _cut_step_scalar_post_root_polish_summary(
    step: ci.CutIntegratorStep,
    role: str,
) -> dict:
    if not isinstance(step, ci.CutIntegratorStep):
        raise TypeError("scalar post-root step summary requires a cut-integrator step")
    return _scalar_post_root_polish_summary(
        (
            (
                role,
                step.ledger.dt_s,
                step.after.transport.time_s,
                step.ledger.scalar_post_root_polish_audit,
            ),
        )
    )


def _combine_scalar_post_root_polish_accounting(
    accounts: Sequence[Mapping[str, object]],
) -> dict:
    """Combine already aggregated, non-overlapping accepted-work partitions."""

    account_tuple = tuple(accounts)
    required = {
        "accepted_step_record_count_considered",
        "used_step_record_count",
        "accepted_root_count_considered",
        "used_accepted_root_count",
        "maximum_scaled_residual_before",
        "maximum_scaled_residual_after",
        "total_scalar_polish_residual_evaluations",
        "all_steps_explicitly_serialized",
    }
    if any(not required <= set(account) for account in account_tuple):
        raise ValueError("scalar post-root work partition is incomplete")
    before = [
        float(account["maximum_scaled_residual_before"])
        for account in account_tuple
        if account["maximum_scaled_residual_before"] is not None
    ]
    after = [
        float(account["maximum_scaled_residual_after"])
        for account in account_tuple
        if account["maximum_scaled_residual_after"] is not None
    ]
    return {
        "numerical_amendment_id": SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID,
        "accepted_step_record_count_considered": sum(
            int(account["accepted_step_record_count_considered"]) for account in account_tuple
        ),
        "used_step_record_count": sum(
            int(account["used_step_record_count"]) for account in account_tuple
        ),
        "accepted_root_count_considered": sum(
            int(account["accepted_root_count_considered"]) for account in account_tuple
        ),
        "used_accepted_root_count": sum(
            int(account["used_accepted_root_count"]) for account in account_tuple
        ),
        "maximum_scaled_residual_before": max(before) if before else None,
        "maximum_scaled_residual_after": max(after) if after else None,
        "total_scalar_polish_residual_evaluations": sum(
            int(account["total_scalar_polish_residual_evaluations"]) for account in account_tuple
        ),
        "all_steps_explicitly_serialized": all(
            account["all_steps_explicitly_serialized"] is True for account in account_tuple
        ),
        "absence_means_zero_uses": True,
        "physics_equations_or_acceptance_tolerances_changed": False,
    }


def _face_event_balance_interval_roots(
    macro: event_orchestrator.CutEventMacrostep,
) -> tuple[_AcceptedBalanceIntervalRoot, _AcceptedBalanceIntervalRoot]:
    """Expose arrival and departure roots with their own balance durations."""

    if macro.ledger.branch is not (
        event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE
    ):
        raise ValueError("balance-root reporting needs arrival plus departure")
    event = macro.face_event_step
    departure = macro.face_departure_step
    if event is None or departure is None:
        raise ValueError("face macrostep lost an accepted balance root")
    if (
        event.before is not macro.before
        or departure.before is not event.after
        or departure.after is not macro.after
        or event.boundary is not macro.boundary
        or departure.boundary is not macro.boundary
    ):
        raise RuntimeError("face balance roots do not form one immutable chain")

    tau = macro.ledger.event_duration_s
    remainder = macro.ledger.departure_remainder_s
    if tau != event.ledger.event_time_s or remainder != departure.ledger.dt_s:
        raise RuntimeError("face macro ledger changed an accepted root duration")
    if not math.isclose(
        math.fsum((tau, remainder)),
        macro.ledger.requested_dt_s,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise RuntimeError("face balance-root durations do not close the macrostep")
    if not math.isclose(
        event.after.time_s,
        math.fsum((macro.before.transport.time_s, tau)),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ) or not math.isclose(
        departure.after.transport.time_s,
        math.fsum((event.after.time_s, remainder)),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise RuntimeError("face balance-root absolute times do not telescope")

    return (
        _AcceptedBalanceIntervalRoot(
            role="exact_face_arrival",
            interval_dt_s=tau,
            time_s=event.after.time_s,
            assembly=event.assembly,
            temperatures_k=(
                *event.after.wet_temperatures_k,
                *event.after.dry_temperatures_k,
                event.after.interface_temperature_k,
            ),
            front_z=event.after.geometry.front.z,
            scalar_post_root_polish_audit=None,
            residual_scale_evaluations=(
                ("exact_face_arrival_certification", event.ledger.residual_scales),
            ),
        ),
        _AcceptedBalanceIntervalRoot(
            role="exact_face_departure_remainder",
            interval_dt_s=remainder,
            time_s=departure.after.transport.time_s,
            assembly=departure.assembly,
            temperatures_k=(
                *departure.after.transport.wet_temperatures_k,
                *departure.after.transport.dry_temperatures_k,
                departure.after.last_interface_temperature_k,
            ),
            front_z=departure.after.transport.geometry.front.z,
            # The face-departure integrator has no scalar post-root polish
            # route.  Report zero uses explicitly instead of reading the
            # optional same-cell-ledger field from a different solver type.
            scalar_post_root_polish_audit=None,
            residual_scale_evaluations=(
                (
                    "strict_face_departure_certification",
                    departure.ledger.residual_scales,
                ),
            ),
        ),
    )


def _accepted_balance_interval_roots(
    step: AcceptedSameCellStep,
) -> tuple[_AcceptedBalanceIntervalRoot, ...]:
    """Return every committed balance root in physical time order."""

    if isinstance(step, AcceptedFaceEventMacrostep):
        return _face_event_balance_interval_roots(step.event_macrostep)
    leaves = step.substeps if isinstance(step, AdaptiveSameCellMacrostep) else (step,)
    return tuple(
        _AcceptedBalanceIntervalRoot(
            role=(
                "adaptive_same_cell_leaf"
                if isinstance(step, AdaptiveSameCellMacrostep)
                else "same_cell_root"
            ),
            interval_dt_s=leaf.ledger.dt_s,
            time_s=leaf.after.transport.time_s,
            assembly=leaf.assembly,
            temperatures_k=(
                *leaf.after.transport.wet_temperatures_k,
                *leaf.after.transport.dry_temperatures_k,
                leaf.after.last_interface_temperature_k,
            ),
            front_z=leaf.after.transport.geometry.front.z,
            scalar_post_root_polish_audit=(leaf.ledger.scalar_post_root_polish_audit),
            residual_scale_evaluations=(
                (("same_cell_original_certification", leaf.ledger.residual_scales),)
                if leaf.ledger.optimization_residual_scales
                is leaf.ledger.residual_scales
                else (
                    (
                        "same_cell_original_certification",
                        leaf.ledger.residual_scales,
                    ),
                    (
                        "same_cell_face_limit_optimization",
                        leaf.ledger.optimization_residual_scales,
                    ),
                )
            ),
        )
        for leaf in leaves
    )


def _json_native_finite_evidence(value: object, *, label: str) -> object:
    """Return a strict JSON-native copy of finite scientific evidence."""

    if is_dataclass(value):
        return _json_native_finite_evidence(asdict(value), label=label)
    if isinstance(value, Enum):
        return _json_native_finite_evidence(value.value, label=label)
    if isinstance(value, Mapping):
        return {
            str(key): _json_native_finite_evidence(
                item,
                label=f"{label}.{key}",
            )
            for key, item in value.items()
        }
    if isinstance(value, (tuple, list, np.ndarray)):
        return [
            _json_native_finite_evidence(item, label=f"{label}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, (float, np.floating)):
        scalar = float(value)
        if not math.isfinite(scalar):
            raise ValueError(f"{label} contains a nonfinite number")
        return scalar
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError(f"{label} contains non-JSON evidence type {type(value)!r}")


def _exact_float(value: object, label: str) -> float:
    """Require a JSON floating-point scalar without coercing booleans/strings."""

    if isinstance(value, bool) or not isinstance(value, float) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite JSON floating-point number")
    return value


def _exact_float_dataclass_payload(
    value: object,
    template: object,
    label: str,
) -> dict[str, float]:
    """Validate one all-float dataclass payload against its exact field set."""

    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be an object")
    expected = set(asdict(template))
    if set(value) != expected:
        raise ValueError(f"{label} has an incomplete or expanded field set")
    return {name: _exact_float(value[name], f"{label}.{name}") for name in expected}


def _coupled_pore_params_from_evidence(value: object) -> cp.CoupledPoreParams:
    """Rebuild the exact pore authority serialized with one A20 replay."""

    if not isinstance(value, Mapping):
        raise TypeError("A20 pore parameters must be an object")
    expected = set(asdict(cp.CoupledPoreParams()))
    if set(value) != expected:
        raise ValueError("A20 pore parameters have an incomplete or expanded field set")
    gab = sp.GabParams(
        **_exact_float_dataclass_payload(
            value["gab"],
            sp.GabParams(),
            "A20 pore.gab",
        )
    )
    oil = sp.OilIsotherm(
        **_exact_float_dataclass_payload(
            value["oil"],
            sp.OilIsotherm(),
            "A20 pore.oil",
        )
    )
    luikov = sp.LuikovParams(
        **_exact_float_dataclass_payload(
            value["luikov"],
            sp.LuikovParams(),
            "A20 pore.luikov",
        )
    )
    temperature_bounds = value["temperature_bounds_k"]
    if not isinstance(temperature_bounds, list) or len(temperature_bounds) != 2:
        raise ValueError("A20 pore temperature bounds must be a two-value JSON array")
    scalar_names = expected - {"gab", "oil", "luikov", "temperature_bounds_k"}
    rebuilt = cp.CoupledPoreParams(
        **{
            name: _exact_float(value[name], f"A20 pore.{name}")
            for name in scalar_names
        },
        gab=gab,
        oil=oil,
        luikov=luikov,
        temperature_bounds_k=tuple(
            _exact_float(item, f"A20 pore.temperature_bounds_k[{index}]")
            for index, item in enumerate(temperature_bounds)
        ),
    )
    if _json_native_finite_evidence(rebuilt, label="rebuilt A20 pore") != dict(value):
        raise RuntimeError("A20 pore parameter reconstruction changed exact values")
    return rebuilt


def _binary_diffusivity_from_evidence(
    value: object,
) -> tc.EffectiveBinaryDiffusivitySelection:
    """Rebuild the exact diffusivity coordinate used by the face mobility."""

    if not isinstance(value, Mapping) or set(value) != {
        "lower_m2_s",
        "upper_m2_s",
        "fraction",
        "provenance",
    }:
        raise ValueError("A20 binary diffusivity evidence has the wrong field set")
    provenance = value["provenance"]
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("A20 binary diffusivity provenance must be nonempty")
    interval = tc.EffectiveBinaryDiffusivityInterval(
        _exact_float(value["lower_m2_s"], "A20 diffusivity.lower_m2_s"),
        _exact_float(value["upper_m2_s"], "A20 diffusivity.upper_m2_s"),
        provenance,
    )
    return interval.select(
        _exact_float(value["fraction"], "A20 diffusivity.fraction")
    )


def _endpoint_enthalpy_selection(
    path: acfp.ActualCompositionForcePath,
    flux_mol_m2_s: float,
    state_field: str,
) -> dict[str, object]:
    """Serialize the unchanged upwind endpoint used by component energy."""

    if not math.isfinite(flux_mol_m2_s):
        raise ValueError("A20 endpoint enthalpy selection needs a finite flux")
    if flux_mol_m2_s < 0.0:
        state = path.right_state
        donor = "right_actual_endpoint"
    else:
        state = path.left_state
        donor = "left_actual_endpoint"
    value = getattr(state, state_field)
    if not math.isfinite(value):
        raise RuntimeError("A20 selected endpoint enthalpy is nonfinite")
    return {
        "donor": donor,
        "selected_enthalpy_j_mol": value,
    }


def _a20_route_geometry(
    root: _AcceptedBalanceIntervalRoot,
    geometry: cg.CutGeometry,
) -> dict[str, object]:
    """Serialize the exact topology facts controlling A20 applicability."""

    dry_indices = tuple(getattr(root.assembly.layout, "dry_cell_indices", ()))
    return {
        "root_role": root.role,
        "front_regime": geometry.front.regime,
        "front_z": geometry.front.z,
        "front_radius_m": geometry.front.radius_m,
        "particle_radius_m": geometry.master_grid.R,
        "front_at_master_face": geometry.front.at_master_face,
        "front_master_face_index": geometry.front.master_face_index,
        "cut_cell_index": geometry.cut_cell_index,
        "dry_piece_count": len(dry_indices),
        "dry_face_flux_count": len(getattr(root.assembly, "dry_face_fluxes", ())),
        "has_finite_positive_cut_subpiece": geometry.cut_cell_index is not None,
    }


def _moving_interface_actual_path_inapplicable_evidence(
    root: _AcceptedBalanceIntervalRoot,
    config: cut.CutTransportConfig,
    geometry: cg.CutGeometry,
    *,
    reason: str,
) -> dict[str, object]:
    """Serialize a topology-owned route on which no finite A20 face exists."""

    route_geometry = _a20_route_geometry(root, geometry)
    if reason == "exact_zero_volume_face_arrival_event_route":
        if not (
            root.role == "exact_face_arrival"
            and geometry.front.at_master_face
            and geometry.cut_cell_index is None
            and not route_geometry["has_finite_positive_cut_subpiece"]
        ):
            raise RuntimeError("A20 face-arrival inapplicability geometry is inconsistent")
    elif reason == "fully_dry_no_moving_interface":
        if geometry.front.regime != "fully_dry":
            raise RuntimeError("A20 fully-dry inapplicability geometry is inconsistent")
    elif reason == "no_moving_interface_dry_face":
        if route_geometry["dry_face_flux_count"] != 0:
            raise RuntimeError("A20 no-face inapplicability geometry is inconsistent")
    else:
        raise ValueError("unknown A20 inapplicability reason")
    binding = {
        "schema_version": ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION,
        "engineering_amendment_id": ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID,
        "accepted_balance_root": {
            "root_role": root.role,
            "balance_interval_dt_s": root.interval_dt_s,
            "time_s": root.time_s,
        },
        "applicability": "topology_route_inapplicable",
        "moving_interface_composition_force_authority": (
            config.moving_interface_composition_force_authority.value
        ),
        "exact_inapplicability_reason": reason,
        "route_geometry": route_geometry,
        "path_replay_required": False,
        "production_path_consumed_as_evidence": False,
        "missing_path_is_not_an_error_on_this_exact_route": True,
        "physics_equations_or_acceptance_tolerances_changed": False,
    }
    bundle = {
        **binding,
        "canonical_binding_sha256": _canonical_evidence_sha256(binding),
        "contract_passed": True,
        "physically_qualifying": False,
    }
    _validate_moving_interface_actual_path_root_evidence(bundle)
    return bundle


def _moving_interface_actual_path_root_evidence(
    root: _AcceptedBalanceIntervalRoot,
    config: cut.CutTransportConfig,
) -> dict[str, object]:
    """Serialize and directly replay one accepted-root Amendment-20 face."""

    if config.moving_interface_composition_force_authority is not (
        cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
    ):
        raise RuntimeError("accepted A20 root lost the actual-(T,y) path authority")
    assembly = root.assembly
    geometry = getattr(assembly, "candidate_geometry", None)
    if geometry is None:
        geometry = getattr(assembly, "event_geometry", None)
    if not isinstance(geometry, cg.CutGeometry):
        raise TypeError("accepted A20 root lost its exact geometry")
    if root.role == "exact_face_arrival":
        return _moving_interface_actual_path_inapplicable_evidence(
            root,
            config,
            geometry,
            reason="exact_zero_volume_face_arrival_event_route",
        )
    if geometry.front.regime == "fully_dry":
        return _moving_interface_actual_path_inapplicable_evidence(
            root,
            config,
            geometry,
            reason="fully_dry_no_moving_interface",
        )
    if not getattr(assembly, "dry_face_fluxes", ()):
        return _moving_interface_actual_path_inapplicable_evidence(
            root,
            config,
            geometry,
            reason="no_moving_interface_dry_face",
        )
    face = assembly.dry_face_fluxes[0]
    if not isinstance(face, cut.DryFaceFlux):
        raise TypeError("accepted A20 root first dry face has the wrong type")
    if face.thermodynamic_force_reconstruction is not (
        cut.DryThermodynamicForceReconstruction.MOVING_INTERFACE_ACTUAL_T_Y_PATH
    ):
        raise RuntimeError("accepted A20 root used the wrong moving-face reconstruction")
    if face.thermodynamic_force_temperature_k is not None:
        raise RuntimeError("accepted A20 root constructed a fictitious force temperature")
    production_path = face.composition_force_path
    if not isinstance(production_path, acfp.ActualCompositionForcePath):
        raise RuntimeError("accepted A20 root omitted its actual composition-force path")

    dry_indices = tuple(getattr(assembly.layout, "dry_cell_indices", ()))
    dry_cells = tuple(getattr(assembly, "dry_cells", ()))
    if not dry_indices or not dry_cells:
        raise RuntimeError("accepted A20 root has no first dry endpoint")
    first_master_cell_index = dry_indices[0]
    first_center_m = cut._dry_piece_center(  # noqa: SLF001
        geometry.cells[first_master_cell_index],
        geometry.front.radius_m,
    )
    distance_m = math.fsum((first_center_m, -geometry.front.radius_m))
    interface = assembly.interface
    first_dry = dry_cells[0]
    oil_labels = tuple(getattr(assembly.before, "oil_fraction_labels", ()))
    if len(oil_labels) != geometry.master_grid.n:
        raise RuntimeError("accepted A20 root lost its material-label vector")
    effective_pore = replace(
        config.dry.pore,
        w_o=oil_labels[first_master_cell_index],
    )
    # Schema v2: the committed face-0 total Stefan flux, as a history-dependent
    # successor operator would consume it at call entry -- the SEED-carried
    # algebraic face rate ("carried, not multiplied"), i.e. the value committed
    # BEFORE this root's own evaluation, not the root's output flux. A20 itself
    # never reads it; recording it makes the replay schema complete for the
    # swap era (R5S6 roadmap item, owner batch ruling 2026-08-11 item 1).
    seed_candidate = getattr(root.assembly, "candidate", None)
    carried_rates = getattr(seed_candidate, "dry_total_stefan_fluxes_mol_m2_s", None)
    if not carried_rates or not math.isfinite(carried_rates[0]):
        raise RuntimeError(
            "accepted A20 root has no finite seed-carried face-0 Stefan flux; "
            "the v2 replay schema fails closed rather than omit it"
        )
    raw_inputs = {
        "left_temperature_k": interface.temperature_k,
        "left_y_hexane": interface.y_hexane,
        "right_temperature_k": first_dry.temperature_k,
        "right_y_hexane": first_dry.y_hexane,
        "pressure_pa": config.dry.pressure_pa,
        "distance_m": distance_m,
        "committed_face0_total_stefan_flux_mol_m2_s": carried_rates[0],
        "first_dry_master_cell_index": first_master_cell_index,
        # Schema v3: the derivative authority the production call actually
        # used, from the path's own audit -- never assumed from a default.
        "closed_hexane_derivative_authority": (
            production_path.phase_audit.closed_hexane_derivative_authority.value
        ),
        "pore_parameters": _json_native_finite_evidence(
            effective_pore,
            label="A20 effective pore parameters",
        ),
        "binary_diffusivity": {
            "lower_m2_s": config.dry.binary_diffusivity.interval.lower_m2_s,
            "upper_m2_s": config.dry.binary_diffusivity.interval.upper_m2_s,
            "fraction": config.dry.binary_diffusivity.fraction,
            "provenance": config.dry.binary_diffusivity.interval.provenance,
        },
    }
    if not (
        production_path.left_state == interface.dry
        and production_path.right_state == first_dry
        and production_path.left_state.temperature_k == raw_inputs["left_temperature_k"]
        and production_path.left_state.y_hexane == raw_inputs["left_y_hexane"]
        and production_path.right_state.temperature_k == raw_inputs["right_temperature_k"]
        and production_path.right_state.y_hexane == raw_inputs["right_y_hexane"]
        and production_path.left_state.pressure_pa == raw_inputs["pressure_pa"]
        and production_path.right_state.pressure_pa == raw_inputs["pressure_pa"]
        and production_path.distance_m == distance_m
    ):
        raise RuntimeError("accepted A20 production path detached from exact face endpoints")

    replay = acfp.evaluate_actual_composition_force_path(
        raw_inputs["left_temperature_k"],
        raw_inputs["left_y_hexane"],
        raw_inputs["right_temperature_k"],
        raw_inputs["right_y_hexane"],
        raw_inputs["pressure_pa"],
        raw_inputs["distance_m"],
        effective_pore,
        closed_hexane_derivative_authority=acfp.ClosedHexaneDerivativeAuthority(
            raw_inputs["closed_hexane_derivative_authority"]
        ),
    )
    if replay != production_path:
        raise RuntimeError("accepted A20 path did not exactly equal its independent replay")
    if not (
        face.gas_chemical_force_gradient_m_inv
        == production_path.gas_force_gradient_m_inv
        == replay.gas_force_gradient_m_inv
        and face.retained_water_chemical_force_gradient_m_inv
        == production_path.retained_force_gradient_m_inv
        == replay.retained_force_gradient_m_inv
    ):
        raise RuntimeError("accepted A20 production composition gradients changed path identity")

    symmetric_face = tc.symmetric_binary_face_state(
        replay.left_state.binary_gas.molar_density_mol_m3,
        replay.left_state.y_water,
        replay.left_state.y_hexane,
        replay.right_state.binary_gas.molar_density_mol_m3,
        replay.right_state.y_water,
        replay.right_state.y_hexane,
    )
    expected_mobility = tc.evaluate_binary_mobility_interval(
        config.dry.binary_diffusivity.interval,
        symmetric_face,
    ).select(config.dry.binary_diffusivity.fraction)
    # GT-PS-2-P1E-20 section 3.4 / revision record R1.
    #
    # The previous contract compared the production mobility bit-for-bit against
    # one rebuilt from the serialized endpoint primitives.  That is unsatisfiable:
    # the rebuild reconstructs the intermediate symmetric face state and follows a
    # different floating-point operation order, landing exactly one ULP away
    # (measured 1.488116110436873e-09 against 1.4881161104368728e-09, relative
    # 1.389644e-16).  The two agree physically; the contract, not the physics, was
    # wrong, and the finite-positive moving-face route could never satisfy it.
    #
    # The check is therefore split at the point where the divergence actually
    # arises, rather than relaxed at the end where a tolerance would silently
    # cover every upstream step at once:
    #
    #   (a) mobility recomputed from the PRODUCTION face state  -> exact;
    #   (b) production face state against the state rebuilt from the serialized
    #       primitives                                          -> bounded.
    #
    # (a) is the mobility law itself and stays an exact contract.  (b) isolates
    # the reconstruction round-trip, which is the only step that can differ, and
    # is the sole place a floating-point budget is admitted.
    production_face_state = face.selected_binary_mobility_face_state
    if production_face_state is None:
        raise RuntimeError(
            "A20 accepted moving face did not serialize its mobility face state; "
            "section 3.4 recomputability primitive is missing"
        )
    # The adaptation to the coupled-pore selection is part of the production
    # arithmetic, not an implementation detail: it evaluates the mobility as
    # ``lower + fraction*(upper - lower)`` over the evaluated interval bounds,
    # whereas ``EvaluatedBinaryMobilitySelection`` evaluates the algebraically
    # identical ``c_g * D * y_w*y_h`` directly.  The two group the same three
    # factors differently and therefore differ in the last bit.  Recomputing
    # through the identical adaptation is what makes this contract exact rather
    # than merely close.
    exact_mobility = (
        tc.evaluate_binary_mobility_interval(
            config.dry.binary_diffusivity.interval,
            production_face_state,
        )
        .select(config.dry.binary_diffusivity.fraction)
        .as_coupled_pore_selection()
    )
    if face.selected_binary_mobility_mol_m_s != exact_mobility.value_mol_m_s:
        raise RuntimeError(
            "accepted A20 binary mobility is not the exact mobility of its own "
            "serialized face state"
        )
    _require_a20_face_state_round_trip(production_face_state, symmetric_face)

    component_fluxes = {
        "gas_water_flux_mol_m2_s": face.component.gas_water_flux_mol_m2_s,
        "gas_hexane_flux_mol_m2_s": face.component.gas_hexane_flux_mol_m2_s,
        "retained_water_flux_mol_m2_s": (face.component.retained_water_flux_mol_m2_s),
    }
    enthalpy_selection = {
        "gas_water": _endpoint_enthalpy_selection(
            replay,
            component_fluxes["gas_water_flux_mol_m2_s"],
            "water_gas_partial_enthalpy_j_mol",
        ),
        "gas_hexane": _endpoint_enthalpy_selection(
            replay,
            component_fluxes["gas_hexane_flux_mol_m2_s"],
            "hexane_gas_partial_enthalpy_j_mol",
        ),
        "retained_water": _endpoint_enthalpy_selection(
            replay,
            component_fluxes["retained_water_flux_mol_m2_s"],
            "retained_water_enthalpy_j_mol",
        ),
    }
    enthalpy_fluxes = {
        "gas_water_enthalpy_flux_w_m2": face.energy.gas_water_enthalpy_flux_w_m2,
        "gas_hexane_enthalpy_flux_w_m2": (face.energy.gas_hexane_enthalpy_flux_w_m2),
        "retained_water_enthalpy_flux_w_m2": (
            face.energy.retained_water_enthalpy_flux_w_m2
        ),
    }
    expected_enthalpy_fluxes = {
        "gas_water_enthalpy_flux_w_m2": (
            component_fluxes["gas_water_flux_mol_m2_s"]
            * enthalpy_selection["gas_water"]["selected_enthalpy_j_mol"]
        ),
        "gas_hexane_enthalpy_flux_w_m2": (
            component_fluxes["gas_hexane_flux_mol_m2_s"]
            * enthalpy_selection["gas_hexane"]["selected_enthalpy_j_mol"]
        ),
        "retained_water_enthalpy_flux_w_m2": (
            component_fluxes["retained_water_flux_mol_m2_s"]
            * enthalpy_selection["retained_water"]["selected_enthalpy_j_mol"]
        ),
    }
    if enthalpy_fluxes != expected_enthalpy_fluxes:
        raise RuntimeError("accepted A20 energy flux lost actual-endpoint upwinding")

    total_stefan_flux = face.component.total_stefan_flux_mol_m2_s
    if total_stefan_flux > 0.0:
        advective_donor = "left_actual_endpoint"
        advective_y_hexane = replay.left_state.y_hexane
    elif total_stefan_flux < 0.0:
        advective_donor = "right_actual_endpoint"
        advective_y_hexane = replay.right_state.y_hexane
    else:
        advective_donor = "symmetric_actual_endpoint_average_at_zero_stefan_flux"
        advective_y_hexane = 0.5 * (
            replay.left_state.y_hexane + replay.right_state.y_hexane
        )
    independent_hexane_flux = -face.independent_water_flux_mol_m2_s
    if not (
        face.component.gas_hexane_flux_mol_m2_s
        == math.fsum((advective_y_hexane * total_stefan_flux, independent_hexane_flux))
    ):
        raise RuntimeError("accepted A20 advective hexane endpoint selection changed")
    advective_selection = {
        "total_stefan_flux_mol_m2_s": total_stefan_flux,
        "donor": advective_donor,
        "selected_y_hexane": advective_y_hexane,
        "independent_hexane_flux_mol_m2_s": independent_hexane_flux,
        "gas_hexane_flux_mol_m2_s": face.component.gas_hexane_flux_mol_m2_s,
    }
    entropy = {
        "gas_entropy_production_w_m3_k": face.gas_entropy_production_w_m3_k,
        "retained_water_entropy_production_w_m3_k": (
            face.retained_water_entropy_production_w_m3_k
        ),
        "conduction_entropy_production_w_m3_k": (
            face.conduction_entropy_production_w_m3_k
        ),
        "total_entropy_production_w_m3_k": face.total_entropy_production_w_m3_k,
    }
    if not (
        all(math.isfinite(value) and value >= 0.0 for value in entropy.values())
        and entropy["total_entropy_production_w_m3_k"]
        == math.fsum(
            (
                entropy["gas_entropy_production_w_m3_k"],
                entropy["retained_water_entropy_production_w_m3_k"],
                entropy["conduction_entropy_production_w_m3_k"],
            )
        )
    ):
        raise RuntimeError("accepted A20 local entropy separation is invalid")
    heat = {
        "fourier_heat_flux_w_m2": face.fourier_heat_flux_w_m2,
        "reciprocal_heat_of_transport_flux_w_m2": (
            face.reciprocal_heat_of_transport_flux_w_m2
        ),
        "reduced_heat_flux_w_m2": face.reduced_heat_flux_w_m2,
        "dry_face_conductive_heat_flux_w_m2": face.conductive_heat_flux_w_m2,
        "component_energy_conductive_heat_flux_w_m2": (
            face.energy.conductive_heat_flux_w_m2
        ),
    }
    expected_reduced_heat = (
        heat["fourier_heat_flux_w_m2"]
        if heat["reciprocal_heat_of_transport_flux_w_m2"] == 0.0
        else math.fsum(
            (
                heat["fourier_heat_flux_w_m2"],
                heat["reciprocal_heat_of_transport_flux_w_m2"],
            )
        )
    )
    if not (
        heat["reduced_heat_flux_w_m2"] == expected_reduced_heat
        and heat["dry_face_conductive_heat_flux_w_m2"]
        == heat["fourier_heat_flux_w_m2"]
        and heat["component_energy_conductive_heat_flux_w_m2"]
        == heat["reduced_heat_flux_w_m2"]
    ):
        raise RuntimeError("accepted A20 heat-flux separation is invalid")

    phase_audit = production_path.phase_audit
    if not (
        phase_audit.deterministic_guard_passed
        and phase_audit.continuous_segment_certificate
        and phase_audit.conditional_on_frozen_engineering_envelopes
        and not phase_audit.unconditional_property_interval_arithmetic_proof
        and phase_audit.no_clipping_or_projection_used
        and phase_audit.no_metastable_storage_used
    ):
        raise RuntimeError("accepted A20 path lost its conditional continuous certificate")
    production_path_json = _json_native_finite_evidence(
        production_path,
        label="A20 production path dataclass",
    )
    production_values = {
        "composition_force_path_mode": production_path.mode.value,
        "gas_composition_force_gradient_m_inv": (
            face.gas_chemical_force_gradient_m_inv
        ),
        "retained_water_composition_force_gradient_m_inv": (
            face.retained_water_chemical_force_gradient_m_inv
        ),
        "binary_thermal_force_gradient_m_inv": (
            face.binary_thermal_force_gradient_m_inv
        ),
        "retained_water_thermal_force_gradient_m_inv": (
            face.retained_water_thermal_force_gradient_m_inv
        ),
        "selected_binary_mobility_mol_m_s": (
            face.selected_binary_mobility_mol_m_s
        ),
        # A20 section 3.4 recomputability primitive (revision record R1): the
        # exact symmetric face state the mobility was evaluated from, so a
        # verifier can reproduce it without rebuilding it from primitives.
        "selected_binary_mobility_face_state": {
            name: getattr(production_face_state, name)
            for name in A20_FACE_STATE_ROUND_TRIP_FIELDS
        },
        "thermodynamic_force_temperature_k": face.thermodynamic_force_temperature_k,
        "component_molar_fluxes": component_fluxes,
        "actual_endpoint_advective_hexane_selection": advective_selection,
        "selected_actual_endpoint_enthalpies": enthalpy_selection,
        "component_enthalpy_fluxes": enthalpy_fluxes,
        "local_entropy_production": entropy,
        "heat_flux_separation": heat,
    }
    binding = {
        "schema_version": ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION,
        "engineering_amendment_id": ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID,
        "accepted_balance_root": {
            "root_role": root.role,
            "balance_interval_dt_s": root.interval_dt_s,
            "time_s": root.time_s,
        },
        "applicability": "finite_positive_moving_interface_face",
        "moving_interface_composition_force_authority": (
            config.moving_interface_composition_force_authority.value
        ),
        "thermodynamic_force_reconstruction": (
            face.thermodynamic_force_reconstruction.value
        ),
        "raw_replay_inputs": raw_inputs,
        "production_values": production_values,
        "conditional_continuous_phase_audit": _json_native_finite_evidence(
            phase_audit,
            label="A20 conditional continuous phase audit",
        ),
        "production_path_dataclass_identity_sha256": (
            _canonical_evidence_sha256(production_path_json)
        ),
        "separation_and_replay_checks": {
            "replay_dataclass_exactly_equal_to_production": True,
            "production_gradients_exactly_equal_to_path": True,
            "actual_endpoint_states_exactly_equal_to_assembly_endpoints": True,
            "binary_mobility_recomputed_from_actual_endpoint_states": True,
            "component_enthalpies_selected_only_from_actual_endpoints": True,
            "composition_and_thermal_force_terms_serialized_separately": True,
            "nonnegative_local_entropy_terms_serialized_separately": True,
            "fourier_reciprocal_and_reduced_heat_serialized_separately": True,
            "actual_endpoint_advective_hexane_selection_serialized": True,
            "fictitious_thermodynamic_force_temperature_constructed": False,
            "full_conditional_continuous_phase_audit_serialized": True,
            "balance_row_consumption_identity_verified_here": False,
            "balance_row_consumption_identity_deferred_to_dynamic_ledgers": True,
            "balance_row_consumption_identity_defer_reason": (
                "this root contract replays the constitutive face; conservative "
                "row consumption remains owned by the accepted component/energy ledgers"
            ),
        },
    }
    bundle = {
        **binding,
        "canonical_binding_sha256": _canonical_evidence_sha256(binding),
        "contract_passed": True,
        "physically_qualifying": False,
    }
    _validate_moving_interface_actual_path_root_evidence(bundle)
    return bundle


def _validate_a20_inapplicable_root_evidence(
    value: Mapping[str, object],
) -> dict[str, object]:
    """Validate an exact topology route with no finite-positive A20 face."""

    binding_fields = {
        "schema_version",
        "engineering_amendment_id",
        "accepted_balance_root",
        "applicability",
        "moving_interface_composition_force_authority",
        "exact_inapplicability_reason",
        "route_geometry",
        "path_replay_required",
        "production_path_consumed_as_evidence",
        "missing_path_is_not_an_error_on_this_exact_route",
        "physics_equations_or_acceptance_tolerances_changed",
    }
    if set(value) != binding_fields | {
        "canonical_binding_sha256",
        "contract_passed",
        "physically_qualifying",
    }:
        raise ValueError("A20 inapplicable-root evidence has the wrong field set")
    if not (
        value["schema_version"] == ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION
        and value["engineering_amendment_id"]
        == ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID
        and value["applicability"] == "topology_route_inapplicable"
        and value["moving_interface_composition_force_authority"]
        == cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH.value
        and value["path_replay_required"] is False
        and value["production_path_consumed_as_evidence"] is False
        and value["missing_path_is_not_an_error_on_this_exact_route"] is True
        and value["physics_equations_or_acceptance_tolerances_changed"] is False
        and value["contract_passed"] is True
        and value["physically_qualifying"] is False
    ):
        raise RuntimeError("A20 inapplicable-root authority or flags changed")
    binding = {name: value[name] for name in binding_fields}
    if value["canonical_binding_sha256"] != _canonical_evidence_sha256(binding):
        raise RuntimeError("A20 inapplicable-root canonical binding hash mismatch")
    root = value["accepted_balance_root"]
    if not isinstance(root, Mapping) or set(root) != {
        "root_role",
        "balance_interval_dt_s",
        "time_s",
    }:
        raise ValueError("A20 inapplicable accepted-root envelope is incomplete")
    if not (
        isinstance(root["root_role"], str)
        and root["root_role"].strip()
        and _exact_float(root["balance_interval_dt_s"], "A20 inapplicable root dt")
        > 0.0
        and _exact_float(root["time_s"], "A20 inapplicable root time") >= 0.0
    ):
        raise ValueError("A20 inapplicable accepted-root envelope is invalid")
    geometry = value["route_geometry"]
    required_geometry = {
        "root_role",
        "front_regime",
        "front_z",
        "front_radius_m",
        "particle_radius_m",
        "front_at_master_face",
        "front_master_face_index",
        "cut_cell_index",
        "dry_piece_count",
        "dry_face_flux_count",
        "has_finite_positive_cut_subpiece",
    }
    if not isinstance(geometry, Mapping) or set(geometry) != required_geometry:
        raise ValueError("A20 inapplicable route geometry is incomplete")
    if geometry["root_role"] != root["root_role"]:
        raise RuntimeError("A20 inapplicable route geometry names the wrong root")
    for name in ("front_z", "front_radius_m", "particle_radius_m"):
        _exact_float(geometry[name], f"A20 route geometry {name}")
    for name in ("dry_piece_count", "dry_face_flux_count"):
        count = geometry[name]
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f"A20 route geometry {name} is invalid")
    reason = value["exact_inapplicability_reason"]
    if reason == "exact_zero_volume_face_arrival_event_route":
        if not (
            root["root_role"] == "exact_face_arrival"
            and geometry["front_regime"] == "partial"
            and geometry["front_at_master_face"] is True
            and isinstance(geometry["front_master_face_index"], int)
            and not isinstance(geometry["front_master_face_index"], bool)
            and geometry["cut_cell_index"] is None
            and geometry["has_finite_positive_cut_subpiece"] is False
        ):
            raise RuntimeError("A20 face-arrival inapplicability facts changed")
    elif reason == "fully_dry_no_moving_interface":
        if not (
            geometry["front_regime"] == "fully_dry"
            and geometry["front_z"] == 0.0
            and geometry["front_radius_m"] == 0.0
            and geometry["has_finite_positive_cut_subpiece"] is False
        ):
            raise RuntimeError("A20 fully-dry inapplicability facts changed")
    elif reason == "no_moving_interface_dry_face":
        if geometry["dry_face_flux_count"] != 0:
            raise RuntimeError("A20 no-face inapplicability facts changed")
    else:
        raise ValueError("A20 inapplicability reason is unknown")
    return {
        "accepted_balance_root": dict(root),
        "applicability": "topology_route_inapplicable",
        "exact_inapplicability_reason": reason,
        "canonical_binding_sha256": value["canonical_binding_sha256"],
        "path_replayed": False,
        "contract_passed": True,
    }


A20_FACE_STATE_ROUND_TRIP_ULP_BUDGET = 4
A20_FACE_STATE_ROUND_TRIP_FIELDS = (
    "left_molar_density_mol_m3",
    "right_molar_density_mol_m3",
    "molar_density_mol_m3",
    "left_y_water",
    "left_y_hexane",
    "right_y_water",
    "right_y_hexane",
    "y_water_y_hexane",
)


def _binary64_ulp_distance(left: float, right: float) -> int:
    """Return the exact representable-step distance between two binary64 values."""

    if left == right:
        return 0
    if not (math.isfinite(left) and math.isfinite(right)):
        raise ValueError("ULP distance requires finite values")
    ordered = []
    for value in (left, right):
        bits = struct.unpack("<q", struct.pack("<d", value))[0]
        ordered.append(-0x8000000000000000 - bits if bits < 0 else bits)
    return abs(ordered[0] - ordered[1])


def _require_a20_face_state_round_trip(
    production: object,
    replayed: object,
) -> None:
    """Bound the ONLY step that can legitimately differ (revision record R1).

    The production symmetric face state and the state rebuilt from the
    serialized endpoint primitives describe the same physical face.  They can
    differ in their last bits because the rebuild recomputes the endpoint
    equilibrium states rather than reusing them.  That reconstruction round-trip
    is isolated here and bounded by a frozen ULP budget.

    This is the sole floating-point budget in the accepted-root contract.  The
    mobility law itself remains an exact comparison against the production face
    state, and the A20 path-dataclass replay remains bit-exact.  Sign, finiteness
    and positivity are rejected before the budget is consulted, so the budget can
    never admit a physically different state.
    """

    for name in A20_FACE_STATE_ROUND_TRIP_FIELDS:
        produced = getattr(production, name)
        rebuilt = getattr(replayed, name)
        if not (isinstance(produced, float) and isinstance(rebuilt, float)):
            raise ValueError(f"A20 face-state field {name} is not a JSON float")
        if not (math.isfinite(produced) and math.isfinite(rebuilt)):
            raise RuntimeError(f"A20 face-state field {name} is not finite")
        if (produced > 0.0) != (rebuilt > 0.0) or (produced < 0.0) != (rebuilt < 0.0):
            raise RuntimeError(f"A20 face-state field {name} changed sign on replay")
        distance = _binary64_ulp_distance(produced, rebuilt)
        if distance > A20_FACE_STATE_ROUND_TRIP_ULP_BUDGET:
            raise RuntimeError(
                "A20 production face state and its primitive replay differ beyond "
                f"the frozen round-trip budget: {name} is {distance} ULP > "
                f"{A20_FACE_STATE_ROUND_TRIP_ULP_BUDGET} "
                f"({produced!r} vs {rebuilt!r})"
            )


def _validate_a20_advective_entropy_and_heat_evidence(
    production: Mapping[str, object],
    replay: object,
) -> None:
    """Recompute the Amendment-20 section 3.4 advective, entropy, and heat terms.

    These three bundles are required by A20 section 3.4 but were absent from
    the accepted-root expected field set.  Every identity below is recomputed
    from the serialized evidence itself, and the advective donor is recomputed
    against the independently replayed endpoint states, so the verifier does
    not inherit the producer's own arithmetic.

    Scope limit: the entropy and heat magnitudes are checked for finiteness,
    sign, and exact internal closure only.  They are NOT independently replayed
    from the constitutive laws, because they originate in the face-flux object
    rather than the actual-path replay.  That replay remains open A20
    qualification work and must not be claimed as complete.
    """

    advective = production["actual_endpoint_advective_hexane_selection"]
    entropy = production["local_entropy_production"]
    heat = production["heat_flux_separation"]
    if not isinstance(advective, Mapping) or set(advective) != {
        "total_stefan_flux_mol_m2_s",
        "donor",
        "selected_y_hexane",
        "independent_hexane_flux_mol_m2_s",
        "gas_hexane_flux_mol_m2_s",
    }:
        raise ValueError("A20 actual-endpoint advective selection is incomplete")
    if not isinstance(entropy, Mapping) or set(entropy) != {
        "gas_entropy_production_w_m3_k",
        "retained_water_entropy_production_w_m3_k",
        "conduction_entropy_production_w_m3_k",
        "total_entropy_production_w_m3_k",
    }:
        raise ValueError("A20 local entropy production evidence is incomplete")
    if not isinstance(heat, Mapping) or set(heat) != {
        "fourier_heat_flux_w_m2",
        "reciprocal_heat_of_transport_flux_w_m2",
        "reduced_heat_flux_w_m2",
        "dry_face_conductive_heat_flux_w_m2",
        "component_energy_conductive_heat_flux_w_m2",
    }:
        raise ValueError("A20 heat-flux separation evidence is incomplete")

    stefan = _exact_float(
        advective["total_stefan_flux_mol_m2_s"],
        "A20 advective total Stefan flux",
    )
    selected_y = _exact_float(
        advective["selected_y_hexane"],
        "A20 advective selected composition",
    )
    independent_hexane = _exact_float(
        advective["independent_hexane_flux_mol_m2_s"],
        "A20 advective independent n-hexane flux",
    )
    gas_hexane = _exact_float(
        advective["gas_hexane_flux_mol_m2_s"],
        "A20 advective gas n-hexane flux",
    )
    # The upwind donor is a deterministic function of the Stefan sign at the
    # two independently replayed endpoints; no fabricated midpoint is allowed.
    if stefan > 0.0:
        expected_donor = "left_actual_endpoint"
        expected_y = replay.left_state.y_hexane
    elif stefan < 0.0:
        expected_donor = "right_actual_endpoint"
        expected_y = replay.right_state.y_hexane
    else:
        expected_donor = "symmetric_actual_endpoint_average_at_zero_stefan_flux"
        expected_y = 0.5 * (replay.left_state.y_hexane + replay.right_state.y_hexane)
    if advective["donor"] != expected_donor or selected_y != expected_y:
        raise RuntimeError("A20 advective endpoint selection did not replay")
    if gas_hexane != math.fsum((selected_y * stefan, independent_hexane)):
        raise RuntimeError("A20 advective n-hexane flux identity did not recompute")

    entropy_terms = {
        name: _exact_float(entropy[name], f"A20 {name}") for name in entropy
    }
    if any(value < 0.0 for value in entropy_terms.values()):
        raise RuntimeError("A20 serialized local entropy production is negative")
    if entropy_terms["total_entropy_production_w_m3_k"] != math.fsum(
        (
            entropy_terms["gas_entropy_production_w_m3_k"],
            entropy_terms["retained_water_entropy_production_w_m3_k"],
            entropy_terms["conduction_entropy_production_w_m3_k"],
        )
    ):
        raise RuntimeError("A20 serialized entropy separation does not close exactly")

    heat_terms = {name: _exact_float(heat[name], f"A20 {name}") for name in heat}
    fourier = heat_terms["fourier_heat_flux_w_m2"]
    reciprocal = heat_terms["reciprocal_heat_of_transport_flux_w_m2"]
    expected_reduced = (
        fourier if reciprocal == 0.0 else math.fsum((fourier, reciprocal))
    )
    if not (
        heat_terms["reduced_heat_flux_w_m2"] == expected_reduced
        and heat_terms["dry_face_conductive_heat_flux_w_m2"] == fourier
        and heat_terms["component_energy_conductive_heat_flux_w_m2"]
        == heat_terms["reduced_heat_flux_w_m2"]
    ):
        raise RuntimeError("A20 serialized heat-flux separation does not recompute")


def _validate_moving_interface_actual_path_root_evidence(
    value: object,
) -> dict[str, object]:
    """Independently replay or route-audit one accepted-root A20 bundle."""

    if not isinstance(value, Mapping):
        raise TypeError("A20 accepted-root evidence must be an object")
    if value.get("applicability") == "topology_route_inapplicable":
        return _validate_a20_inapplicable_root_evidence(value)
    if value.get("applicability") != "finite_positive_moving_interface_face":
        raise RuntimeError("A20 accepted-root applicability is missing or unknown")
    binding_fields = {
        "schema_version",
        "engineering_amendment_id",
        "accepted_balance_root",
        "applicability",
        "moving_interface_composition_force_authority",
        "thermodynamic_force_reconstruction",
        "raw_replay_inputs",
        "production_values",
        "conditional_continuous_phase_audit",
        "production_path_dataclass_identity_sha256",
        "separation_and_replay_checks",
    }
    expected_fields = binding_fields | {
        "canonical_binding_sha256",
        "contract_passed",
        "physically_qualifying",
    }
    if set(value) != expected_fields:
        raise ValueError("A20 accepted-root evidence has an incomplete or expanded field set")
    if not (
        value["schema_version"] == ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION
        and value["engineering_amendment_id"]
        == ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID
        and value["applicability"] == "finite_positive_moving_interface_face"
        and value["moving_interface_composition_force_authority"]
        == cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH.value
        and value["thermodynamic_force_reconstruction"]
        == cut.DryThermodynamicForceReconstruction.MOVING_INTERFACE_ACTUAL_T_Y_PATH.value
        and value["contract_passed"] is True
        and value["physically_qualifying"] is False
    ):
        raise RuntimeError("A20 accepted-root authority or qualification flags changed")
    binding = {name: value[name] for name in binding_fields}
    if value["canonical_binding_sha256"] != _canonical_evidence_sha256(binding):
        raise RuntimeError("A20 accepted-root canonical binding hash mismatch")

    root = value["accepted_balance_root"]
    if not isinstance(root, Mapping) or set(root) != {
        "root_role",
        "balance_interval_dt_s",
        "time_s",
    }:
        raise ValueError("A20 accepted balance root envelope is incomplete")
    if not (
        isinstance(root["root_role"], str)
        and root["root_role"].strip()
        and _exact_float(root["balance_interval_dt_s"], "A20 root dt") > 0.0
        and _exact_float(root["time_s"], "A20 root time") >= 0.0
    ):
        raise ValueError("A20 accepted balance root envelope is invalid")

    raw = value["raw_replay_inputs"]
    if not isinstance(raw, Mapping) or set(raw) != {
        "left_temperature_k",
        "left_y_hexane",
        "right_temperature_k",
        "right_y_hexane",
        "pressure_pa",
        "distance_m",
        "committed_face0_total_stefan_flux_mol_m2_s",
        "first_dry_master_cell_index",
        "closed_hexane_derivative_authority",
        "pore_parameters",
        "binary_diffusivity",
    }:
        raise ValueError("A20 raw replay inputs are incomplete")
    # v3: the recorded authority must name a real enum member; the replay
    # below fails closed on an unknown value rather than defaulting.
    if raw["closed_hexane_derivative_authority"] not in {
        member.value for member in acfp.ClosedHexaneDerivativeAuthority
    }:
        raise ValueError("A20 closed-hexane derivative authority is invalid")
    # v2: committed flux is sign-free (evaporative positive, condensing
    # negative) but must be a finite float -- the Pe input a successor
    # operator consumes cannot be reconstructed if this is absent or NaN.
    _exact_float(
        raw["committed_face0_total_stefan_flux_mol_m2_s"],
        "A20 committed face-0 Stefan flux",
    )
    first_index = raw["first_dry_master_cell_index"]
    if isinstance(first_index, bool) or not isinstance(first_index, int) or first_index < 0:
        raise ValueError("A20 first dry master-cell index is invalid")
    pore = _coupled_pore_params_from_evidence(raw["pore_parameters"])
    diffusivity = _binary_diffusivity_from_evidence(raw["binary_diffusivity"])
    replay = acfp.evaluate_actual_composition_force_path(
        _exact_float(raw["left_temperature_k"], "A20 left temperature"),
        _exact_float(raw["left_y_hexane"], "A20 left composition"),
        _exact_float(raw["right_temperature_k"], "A20 right temperature"),
        _exact_float(raw["right_y_hexane"], "A20 right composition"),
        _exact_float(raw["pressure_pa"], "A20 pressure"),
        _exact_float(raw["distance_m"], "A20 distance"),
        pore,
        closed_hexane_derivative_authority=acfp.ClosedHexaneDerivativeAuthority(
            raw["closed_hexane_derivative_authority"]
        ),
    )
    replay_json = _json_native_finite_evidence(
        replay,
        label="independently replayed A20 path",
    )
    if value["production_path_dataclass_identity_sha256"] != (
        _canonical_evidence_sha256(replay_json)
    ):
        raise RuntimeError("A20 replay path dataclass identity mismatch")
    replay_phase = _json_native_finite_evidence(
        replay.phase_audit,
        label="independently replayed A20 phase audit",
    )
    if value["conditional_continuous_phase_audit"] != replay_phase:
        raise RuntimeError("A20 conditional continuous phase audit mismatch")
    phase = replay.phase_audit
    if not (
        phase.deterministic_guard_passed
        and phase.continuous_segment_certificate
        and phase.conditional_on_frozen_engineering_envelopes
        and not phase.unconditional_property_interval_arithmetic_proof
        and phase.no_clipping_or_projection_used
        and phase.no_metastable_storage_used
    ):
        raise RuntimeError("A20 replay did not retain the conditional phase certificate")

    production = value["production_values"]
    if not isinstance(production, Mapping) or set(production) != {
        "composition_force_path_mode",
        "gas_composition_force_gradient_m_inv",
        "retained_water_composition_force_gradient_m_inv",
        "binary_thermal_force_gradient_m_inv",
        "retained_water_thermal_force_gradient_m_inv",
        "selected_binary_mobility_mol_m_s",
        "selected_binary_mobility_face_state",
        "thermodynamic_force_temperature_k",
        "component_molar_fluxes",
        "selected_actual_endpoint_enthalpies",
        "component_enthalpy_fluxes",
        # Amendment 20 section 3.4 requires every accepted moving-face flux to
        # serialize its advective composition, reciprocal heat, and entropy
        # terms.  The producer emits all three; this expected set previously
        # omitted them, so the finite-positive moving-face route could never
        # validate its own bundle.  The route was unreachable behind the
        # superseded convergence gate, which is why the contradiction was not
        # observed before GT-PS-2-P1E-21.
        "actual_endpoint_advective_hexane_selection",
        "local_entropy_production",
        "heat_flux_separation",
    }:
        raise ValueError("A20 production values are incomplete")
    if not (
        production["composition_force_path_mode"] == replay.mode.value
        and production["gas_composition_force_gradient_m_inv"]
        == replay.gas_force_gradient_m_inv
        and production["retained_water_composition_force_gradient_m_inv"]
        == replay.retained_force_gradient_m_inv
        and production["thermodynamic_force_temperature_k"] is None
    ):
        raise RuntimeError("A20 production gradients or path mode changed")
    _exact_float(
        production["binary_thermal_force_gradient_m_inv"],
        "A20 binary thermal-force gradient",
    )
    _exact_float(
        production["retained_water_thermal_force_gradient_m_inv"],
        "A20 retained thermal-force gradient",
    )
    # Revision record R1, validator side.  The same split applied at the producer
    # is applied here, for the same reason: the mobility is recomputed EXACTLY
    # from the serialized production face state through the identical
    # coupled-pore adaptation, and the reconstruction of that face state from the
    # serialized endpoint primitives is bounded separately.  Comparing the
    # production mobility directly against a face rebuilt from primitives mixes
    # both steps and cannot be exact.
    replayed_face = tc.symmetric_binary_face_state(
        replay.left_state.binary_gas.molar_density_mol_m3,
        replay.left_state.y_water,
        replay.left_state.y_hexane,
        replay.right_state.binary_gas.molar_density_mol_m3,
        replay.right_state.y_water,
        replay.right_state.y_hexane,
    )
    serialized_face = production["selected_binary_mobility_face_state"]
    if not isinstance(serialized_face, Mapping) or set(serialized_face) != set(
        A20_FACE_STATE_ROUND_TRIP_FIELDS
    ):
        raise ValueError("A20 serialized mobility face state is incomplete")
    production_face = tc.symmetric_binary_face_state(
        _exact_float(serialized_face["left_molar_density_mol_m3"], "A20 left face density"),
        _exact_float(serialized_face["left_y_water"], "A20 left face y_water"),
        _exact_float(serialized_face["left_y_hexane"], "A20 left face y_hexane"),
        _exact_float(serialized_face["right_molar_density_mol_m3"], "A20 right face density"),
        _exact_float(serialized_face["right_y_water"], "A20 right face y_water"),
        _exact_float(serialized_face["right_y_hexane"], "A20 right face y_hexane"),
    )
    # The reconstructed state must reproduce the serialized derived quantities
    # exactly; only the primitive round-trip below is permitted to differ.
    for name in ("molar_density_mol_m3", "y_water_y_hexane"):
        if getattr(production_face, name) != _exact_float(
            serialized_face[name], f"A20 serialized face {name}"
        ):
            raise RuntimeError(
                f"A20 serialized face state is internally inconsistent in {name}"
            )
    expected_mobility = (
        tc.evaluate_binary_mobility_interval(diffusivity.interval, production_face)
        .select(diffusivity.fraction)
        .as_coupled_pore_selection()
    )
    if production["selected_binary_mobility_mol_m_s"] != (
        expected_mobility.value_mol_m_s
    ):
        raise RuntimeError(
            "A20 selected mobility is not the exact mobility of its own serialized "
            "face state"
        )
    _require_a20_face_state_round_trip(production_face, replayed_face)

    fluxes = production["component_molar_fluxes"]
    enthalpies = production["selected_actual_endpoint_enthalpies"]
    energy_fluxes = production["component_enthalpy_fluxes"]
    if not isinstance(fluxes, Mapping) or set(fluxes) != {
        "gas_water_flux_mol_m2_s",
        "gas_hexane_flux_mol_m2_s",
        "retained_water_flux_mol_m2_s",
    }:
        raise ValueError("A20 component flux evidence is incomplete")
    if not isinstance(enthalpies, Mapping) or set(enthalpies) != {
        "gas_water",
        "gas_hexane",
        "retained_water",
    }:
        raise ValueError("A20 endpoint enthalpy selections are incomplete")
    if not isinstance(energy_fluxes, Mapping) or set(energy_fluxes) != {
        "gas_water_enthalpy_flux_w_m2",
        "gas_hexane_enthalpy_flux_w_m2",
        "retained_water_enthalpy_flux_w_m2",
    }:
        raise ValueError("A20 component enthalpy flux evidence is incomplete")
    selection_spec = (
        (
            "gas_water",
            "gas_water_flux_mol_m2_s",
            "water_gas_partial_enthalpy_j_mol",
            "gas_water_enthalpy_flux_w_m2",
        ),
        (
            "gas_hexane",
            "gas_hexane_flux_mol_m2_s",
            "hexane_gas_partial_enthalpy_j_mol",
            "gas_hexane_enthalpy_flux_w_m2",
        ),
        (
            "retained_water",
            "retained_water_flux_mol_m2_s",
            "retained_water_enthalpy_j_mol",
            "retained_water_enthalpy_flux_w_m2",
        ),
    )
    for component, flux_name, state_field, energy_name in selection_spec:
        flux = _exact_float(fluxes[flux_name], f"A20 {flux_name}")
        expected_selection = _endpoint_enthalpy_selection(replay, flux, state_field)
        if enthalpies[component] != expected_selection:
            raise RuntimeError(f"A20 {component} endpoint enthalpy selection changed")
        if energy_fluxes[energy_name] != (
            flux * expected_selection["selected_enthalpy_j_mol"]
        ):
            raise RuntimeError(f"A20 {component} enthalpy flux changed")
    checks = value["separation_and_replay_checks"]
    expected_checks = {
        "replay_dataclass_exactly_equal_to_production": True,
        "production_gradients_exactly_equal_to_path": True,
        "actual_endpoint_states_exactly_equal_to_assembly_endpoints": True,
        "binary_mobility_recomputed_from_actual_endpoint_states": True,
        "component_enthalpies_selected_only_from_actual_endpoints": True,
        "composition_and_thermal_force_terms_serialized_separately": True,
        "fictitious_thermodynamic_force_temperature_constructed": False,
        "full_conditional_continuous_phase_audit_serialized": True,
        # Amendment 20 section 3.4 disclosures that the producer emits and this
        # expected set previously omitted.  The first three cover the advective,
        # entropy, and heat bundles admitted above.  The last three are the
        # producer's honest negative disclosure that the conservative balance-row
        # consumption identity is NOT proven by this constitutive root contract
        # and remains owned by the accepted component/energy ledgers; that
        # deferral is preserved verbatim rather than upgraded to a claim.
        "nonnegative_local_entropy_terms_serialized_separately": True,
        "fourier_reciprocal_and_reduced_heat_serialized_separately": True,
        "actual_endpoint_advective_hexane_selection_serialized": True,
        "balance_row_consumption_identity_verified_here": False,
        "balance_row_consumption_identity_deferred_to_dynamic_ledgers": True,
        "balance_row_consumption_identity_defer_reason": (
            "this root contract replays the constitutive face; conservative "
            "row consumption remains owned by the accepted component/energy ledgers"
        ),
    }
    if checks != expected_checks:
        raise RuntimeError("A20 separation/replay declarations changed")

    _validate_a20_advective_entropy_and_heat_evidence(production, replay)

    return {
        "accepted_balance_root": dict(root),
        "production_path_dataclass_identity_sha256": (
            value["production_path_dataclass_identity_sha256"]
        ),
        "canonical_binding_sha256": value["canonical_binding_sha256"],
        "contract_passed": True,
    }


def _moving_interface_actual_path_audit_summary(
    samples: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Revalidate every A20 root sample in one committed step."""

    validations: list[dict[str, object]] = []
    for index, sample in enumerate(samples):
        if not isinstance(sample, Mapping):
            raise TypeError(f"A20 accepted-root sample {index} is not an object")
        evidence = sample.get("moving_interface_actual_composition_force_path_audit")
        validation = _validate_moving_interface_actual_path_root_evidence(evidence)
        root = validation["accepted_balance_root"]
        if not (
            root["root_role"] == sample.get("root_role")
            and root["balance_interval_dt_s"] == sample.get("balance_interval_dt_s")
            and root["time_s"] == sample.get("time_s")
        ):
            raise RuntimeError("A20 accepted-root envelope detached from its sample")
        validations.append(validation)
    identities = [item["canonical_binding_sha256"] for item in validations]
    if len(set(identities)) != len(identities):
        raise RuntimeError("A20 accepted-root evidence identities are not one-to-one")
    # A20 evidence-contract defect 5 (audit doc, 2026-08-03): a committed step
    # may legitimately mix finite-positive moving-face roots with exact
    # topology-route-inapplicable roots (e.g. an exact zero-volume face
    # arrival).  The original summary assumed every root was finite-positive
    # and crashed on the first mixed step ever executed.  The partition below
    # keeps replay claims scoped to the roots that require replay, and counts
    # the inapplicable routes explicitly instead of misdescribing them.
    finite = [
        item
        for item in validations
        if item.get("applicability") != "topology_route_inapplicable"
    ]
    inapplicable = [item for item in validations if item not in finite]
    reasons: dict[str, int] = {}
    for item in inapplicable:
        reason = str(item.get("exact_inapplicability_reason"))
        reasons[reason] = reasons.get(reason, 0) + 1
    return {
        "schema_version": ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION,
        "engineering_amendment_id": ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID,
        "accepted_root_count": len(validations),
        "accepted_root_identity_count": len(identities),
        "accepted_root_identities": validations,
        "finite_positive_root_count": len(finite),
        "inapplicable_route_root_count": len(inapplicable),
        "inapplicable_route_reasons": dict(sorted(reasons.items())),
        # Scoped claim: every root that REQUIRES an independent path replay
        # (the finite-positive ones) was replayed; inapplicable routes carry
        # path_replayed=False by frozen design and are audited, not replayed.
        "all_roots_independently_replayed": bool(validations),
        "all_finite_positive_roots_independently_replayed": bool(validations),
        "all_roots_use_actual_linear_T_y_path_authority": bool(validations),
        "all_roots_serialize_full_conditional_continuous_phase_audit": (
            bool(validations)
        ),
        "fictitious_thermodynamic_force_temperature_count": 0,
        "contract_passed": bool(validations),
        "physically_qualifying": False,
    }


def _surface_flux_sample(
    root: _AcceptedBalanceIntervalRoot,
    *,
    include_constitutive_evidence: bool = True,
) -> dict:
    """Serialize the surface root that closes one accepted balance interval."""

    assembly = root.assembly
    surface = assembly.dry_face_fluxes[-1]
    config = getattr(assembly, "target_config", None)
    if config is None:
        config = assembly.before.config
    if not isinstance(config, cut.CutTransportConfig):
        raise TypeError("interface-root evidence lost its cut configuration")
    interface_root_evidence = _interface_composition_root_evidence(root, config)
    dry_capacity_evidence = _dry_energy_capacity_stencil_evidence(root, config)
    if not include_constitutive_evidence:
        return {
            "root_role": root.role,
            "balance_interval_dt_s": root.interval_dt_s,
            "time_s": root.time_s,
            "scalar_post_root_polish_audit": (
                None
                if root.scalar_post_root_polish_audit is None
                else asdict(root.scalar_post_root_polish_audit)
            ),
            "total_stefan_flux_mol_m2_s": (
                root.assembly.candidate.dry_total_stefan_fluxes_mol_m2_s[-1]
            ),
            "conserved_water_flux_mol_m2_s": (surface.component.conserved_water_flux_mol_m2_s),
            "conserved_hexane_flux_mol_m2_s": (surface.component.conserved_hexane_flux_mol_m2_s),
            "independent_hexane_flux_mol_m2_s": (-surface.independent_water_flux_mol_m2_s),
            **interface_root_evidence,
            **dry_capacity_evidence,
            **(
                {"surface_film_audit": asdict(surface.surface_film_audit)}
                if surface.surface_film_audit is not None
                else {}
            ),
        }
    actual_path_evidence = (
        _moving_interface_actual_path_root_evidence(root, config)
        if config.moving_interface_composition_force_authority
        is cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
        else None
    )
    geometry = getattr(assembly, "candidate_geometry", None)
    if geometry is None:
        geometry = getattr(assembly, "event_geometry", None)
    if not isinstance(geometry, cg.CutGeometry):
        raise TypeError("surface constitutive evidence lost its accepted geometry")
    dry_indices = tuple(assembly.layout.dry_cell_indices)
    if not dry_indices:
        raise RuntimeError("surface constitutive evidence needs a dry outer piece")
    outer_cell = geometry.cells[dry_indices[-1]]
    outer_center_m = cut._dry_piece_center(  # noqa: SLF001
        outer_cell,
        geometry.front.radius_m,
    )
    boundary = assembly.surface_boundary
    if surface.surface_film_audit is None:
        if not isinstance(boundary, ct.DirichletPoreBoundary):
            raise RuntimeError("resolved surface state lost its finite-film audit")
        surface_temperature_k = boundary.temperature_k
        surface_y_hexane = boundary.y_hexane
    else:
        film = surface.surface_film_audit
        surface_temperature_k = film.surface_temperature_k
        surface_y_hexane = film.assumed_surface_y_hexane
        if (
            film.pressure_pa != config.dry.pressure_pa
            or film.binary_gas_interaction_k_wh != config.dry.pore.k_wh
        ):
            raise RuntimeError("surface film and pore constitutive authorities diverged")
    diffusivity = config.dry.binary_diffusivity
    effective_force = math.fsum(
        (
            surface.gas_chemical_force_gradient_m_inv,
            surface.binary_thermal_force_gradient_m_inv,
        )
    )
    if surface.thermodynamic_force_temperature_k is None:
        raise RuntimeError("a physical surface face lost its force temperature")
    primitives = {
        "schema_version": surface_oracle.RAW_SAMPLE_SCHEMA_VERSION,
        "sample_kind": surface_oracle.RAW_SAMPLE_KIND,
        "outer_dry_state": {
            "temperature_k": assembly.candidate.dry_temperatures_k[-1],
            "y_hexane": assembly.candidate.dry_y_hexane[-1],
            "center_radius_m": outer_center_m,
        },
        "resolved_surface_state": {
            "temperature_k": surface_temperature_k,
            "y_hexane": surface_y_hexane,
            "radius_m": geometry.master_grid.R,
        },
        "center_to_surface_distance_m": geometry.master_grid.R - outer_center_m,
        "pressure_pa": config.dry.pressure_pa,
        "binary_gas_interaction_k_wh": config.dry.pore.k_wh,
        "binary_diffusivity": {
            "lower_m2_s": diffusivity.interval.lower_m2_s,
            "upper_m2_s": diffusivity.interval.upper_m2_s,
            "fraction": diffusivity.fraction,
            "selected_m2_s": diffusivity.value_m2_s,
        },
        "mass_force": {
            "mode": surface.mass_force_mode.value,
            "binary_thermal_diffusion_factor": (surface.binary_thermal_diffusion_factor),
        },
        "production_values": {
            "thermodynamic_force_temperature_k": (surface.thermodynamic_force_temperature_k),
            "selected_binary_mobility_mol_m_s": (surface.selected_binary_mobility_mol_m_s),
            "gas_chemical_force_gradient_m_inv": (surface.gas_chemical_force_gradient_m_inv),
            "binary_thermal_force_gradient_m_inv": (surface.binary_thermal_force_gradient_m_inv),
            "effective_binary_force_gradient_m_inv": effective_force,
            "independent_hexane_flux_mol_m2_s": (-surface.independent_water_flux_mol_m2_s),
            "binary_gas_entropy_production_w_m3_k": (surface.gas_entropy_production_w_m3_k),
        },
    }
    constitutive_audit = surface_oracle.audit_surface_constitutive_sample(primitives)
    return {
        "root_role": root.role,
        "balance_interval_dt_s": root.interval_dt_s,
        "time_s": root.time_s,
        "scalar_post_root_polish_audit": (
            None
            if root.scalar_post_root_polish_audit is None
            else asdict(root.scalar_post_root_polish_audit)
        ),
        "total_stefan_flux_mol_m2_s": (
            root.assembly.candidate.dry_total_stefan_fluxes_mol_m2_s[-1]
        ),
        "conserved_water_flux_mol_m2_s": (surface.component.conserved_water_flux_mol_m2_s),
        "conserved_hexane_flux_mol_m2_s": (surface.component.conserved_hexane_flux_mol_m2_s),
        "independent_hexane_flux_mol_m2_s": (-surface.independent_water_flux_mol_m2_s),
        "surface_constitutive_primitives": primitives,
        "surface_constitutive_oracle": constitutive_audit,
        **interface_root_evidence,
        **dry_capacity_evidence,
        **(
            {
                "moving_interface_actual_composition_force_path_audit": (
                    actual_path_evidence
                )
            }
            if actual_path_evidence is not None
            else {}
        ),
        **(
            {"surface_film_audit": asdict(surface.surface_film_audit)}
            if surface.surface_film_audit is not None
            else {}
        ),
    }


def _canonical_evidence_sha256(value: Mapping[str, object]) -> str:
    """Hash finite, JSON-native evidence with one deterministic encoding."""

    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _interface_composition_root_evidence(
    root: _AcceptedBalanceIntervalRoot,
    config: cut.CutTransportConfig,
) -> dict:
    """Serialize one canonical A18 input/audit/certificate binding."""

    interface = root.assembly.interface
    audit = interface.interface_composition_root_audit
    if not isinstance(audit, cut.InterfaceCompositionRootAudit):
        raise TypeError("accepted balance root lost its frozen A18 root audit")
    interface_oil_fraction_label = getattr(
        interface,
        "oil_fraction_label",
        getattr(interface, "dry_oil_fraction_label", None),
    )
    if (
        interface.temperature_k != interface.dry.temperature_k
        or interface.y_hexane != interface.dry.y_hexane
        or interface.y_hexane != audit.returned_y_hexane
        or interface.log_hexane_fugacity_residual
        != audit.returned_log_fugacity_residual
        or interface.dry.hexane_activity
        != audit.returned_product_hexane_activity
        or interface_oil_fraction_label != audit.actual_oil_fraction_label
        or interface.dry.pressure_pa != config.dry.pressure_pa
        or config.interface_composition.y_hexane_bounds
        != (audit.caller_lower_y_hexane, audit.caller_upper_y_hexane)
        or config.interface_composition.log_fugacity_tolerance
        != audit.log_fugacity_tolerance
    ):
        raise RuntimeError("accepted balance root and its A18 certificate diverged")
    effective_pore = replace(
        config.dry.pore,
        w_o=audit.actual_oil_fraction_label,
    )
    if effective_pore.k_wh != config.dry.pore.k_wh:
        raise RuntimeError("A18 evidence changed the binary-gas interaction authority")
    recomputed_dry = cp.evaluate_equilibrium(
        interface.temperature_k,
        config.dry.pressure_pa,
        interface.y_hexane,
        effective_pore,
    )
    if recomputed_dry != interface.dry:
        raise RuntimeError("A18 evidence could not reproduce the exact coupled-pore state")
    recomputed_log_fugacity_residual = front.log_fugacity_residual(
        interface.temperature_k,
        config.dry.pressure_pa,
        1.0 - interface.y_hexane,
        interface.y_hexane,
        k_wh=effective_pore.k_wh,
    )[0]
    _, recomputed_product_hexane_activity = cp._gas_component_activities(  # noqa: SLF001
        interface.temperature_k,
        config.dry.pressure_pa,
        interface.y_hexane,
        effective_pore,
    )
    if (
        recomputed_log_fugacity_residual
        != audit.returned_log_fugacity_residual
        or recomputed_product_hexane_activity
        != audit.returned_product_hexane_activity
    ):
        raise RuntimeError("A18 stored-bit log/product re-evaluation is not bit exact")
    serialized_state = asdict(interface.dry)
    serialized_pore = asdict(effective_pore)
    audit_record = asdict(audit)
    context = {
        "schema_version": INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION,
        "interface_temperature_t_gamma_k": interface.temperature_k,
        "interface_temperature_t_gamma_k_binary64_hex": interface.temperature_k.hex(),
        "pressure_pa": interface.dry.pressure_pa,
        "pressure_pa_binary64_hex": interface.dry.pressure_pa.hex(),
        "binary_gas_interaction_k_wh": effective_pore.k_wh,
        "binary_gas_interaction_k_wh_binary64_hex": effective_pore.k_wh.hex(),
        "actual_oil_fraction_label": audit.actual_oil_fraction_label,
        "actual_oil_fraction_label_binary64_hex": (
            audit.actual_oil_fraction_label.hex()
        ),
        "caller_lower_y_hexane": audit.caller_lower_y_hexane,
        "caller_lower_y_hexane_binary64_hex": audit.caller_lower_y_hexane.hex(),
        "caller_upper_y_hexane": audit.caller_upper_y_hexane,
        "caller_upper_y_hexane_binary64_hex": audit.caller_upper_y_hexane.hex(),
        "log_fugacity_tolerance": audit.log_fugacity_tolerance,
        "log_fugacity_tolerance_binary64_hex": audit.log_fugacity_tolerance.hex(),
        "returned_y_hexane": audit.returned_y_hexane,
        "returned_y_hexane_binary64_hex": audit.returned_y_hexane.hex(),
        "final_lower_y_hexane_binary64_hex": audit.final_lower_y_hexane.hex(),
        "final_upper_y_hexane_binary64_hex": audit.final_upper_y_hexane.hex(),
    }
    certification = {
        "schema_version": INTERFACE_ROOT_CERTIFICATION_SCHEMA_VERSION,
        "coupled_pore_params": serialized_pore,
        "coupled_pore_params_canonical_sha256": (
            _canonical_evidence_sha256(serialized_pore)
        ),
        "equilibrium_pore_state": serialized_state,
        "equilibrium_pore_state_canonical_sha256": (
            _canonical_evidence_sha256(serialized_state)
        ),
    }
    stored_bit_re_evaluation = {
        "schema_version": 1,
        "interface_temperature_t_gamma_k": interface.temperature_k,
        "interface_temperature_t_gamma_k_binary64_hex": interface.temperature_k.hex(),
        "pressure_pa": config.dry.pressure_pa,
        "pressure_pa_binary64_hex": config.dry.pressure_pa.hex(),
        "returned_y_hexane": interface.y_hexane,
        "returned_y_hexane_binary64_hex": interface.y_hexane.hex(),
        "binary_gas_interaction_k_wh": effective_pore.k_wh,
        "binary_gas_interaction_k_wh_binary64_hex": effective_pore.k_wh.hex(),
        "actual_oil_fraction_label": effective_pore.w_o,
        "actual_oil_fraction_label_binary64_hex": effective_pore.w_o.hex(),
        "recomputed_log_fugacity_residual": recomputed_log_fugacity_residual,
        "recomputed_log_fugacity_residual_binary64_hex": (
            recomputed_log_fugacity_residual.hex()
        ),
        "recomputed_product_hexane_activity": recomputed_product_hexane_activity,
        "recomputed_product_hexane_activity_binary64_hex": (
            recomputed_product_hexane_activity.hex()
        ),
        "log_fugacity_residual_bit_exact": True,
        "product_hexane_activity_bit_exact": True,
        "full_coupled_pore_state_bit_exact": True,
        "tolerance_or_approximate_comparison_used": False,
    }
    accepted_balance_root = {
        "root_role": root.role,
        "balance_interval_dt_s": root.interval_dt_s,
        "time_s": root.time_s,
    }
    binding_payload = {
        "accepted_balance_root": accepted_balance_root,
        "root_audit": audit_record,
        "exact_input_context": context,
        "full_coupled_pore_certification": certification,
        "stored_bit_log_product_re_evaluation": stored_bit_re_evaluation,
    }
    bundle = {
        "schema_version": INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": INTERFACE_ROOT_AMENDMENT_ID,
        **binding_payload,
        "canonical_binding_sha256": _canonical_evidence_sha256(binding_payload),
        "physically_qualifying": False,
    }
    _validate_interface_saturation_endpoint_audit(bundle)
    return {"interface_saturation_endpoint_audit": bundle}


def _evidence_mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    return value


def _require_binary64_hex(
    context: Mapping[str, object],
    field: str,
    expected: float,
) -> None:
    if context.get(field) != expected or context.get(f"{field}_binary64_hex") != expected.hex():
        raise RuntimeError(f"A18 exact input context changed {field}")


def _validate_interface_saturation_endpoint_audit(
    value: Mapping[str, object],
) -> dict[str, object]:
    """Fail closed on any incomplete or internally inconsistent A18 bundle."""

    bundle = _evidence_mapping(value, "interface saturation endpoint audit")
    required = {
        "schema_version",
        "numerical_amendment_id",
        "root_audit",
        "accepted_balance_root",
        "exact_input_context",
        "full_coupled_pore_certification",
        "stored_bit_log_product_re_evaluation",
        "canonical_binding_sha256",
        "physically_qualifying",
    }
    if not required <= set(bundle):
        raise ValueError("interface saturation endpoint audit is incomplete")
    if (
        bundle["schema_version"] != INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION
        or bundle["numerical_amendment_id"] != INTERFACE_ROOT_AMENDMENT_ID
        or bundle["physically_qualifying"] is not False
    ):
        raise ValueError("interface saturation endpoint audit changed its authority")
    audit = _evidence_mapping(bundle["root_audit"], "A18 root audit")
    accepted_balance_root = _evidence_mapping(
        bundle["accepted_balance_root"],
        "A18 accepted balance root",
    )
    if not {
        "root_role",
        "balance_interval_dt_s",
        "time_s",
    } <= set(accepted_balance_root) or not (
        isinstance(accepted_balance_root["root_role"], str)
        and accepted_balance_root["root_role"].strip()
        and math.isfinite(float(accepted_balance_root["balance_interval_dt_s"]))
        and float(accepted_balance_root["balance_interval_dt_s"]) > 0.0
        and math.isfinite(float(accepted_balance_root["time_s"]))
        and float(accepted_balance_root["time_s"]) >= 0.0
    ):
        raise ValueError("A18 accepted balance-root identity is invalid")
    context = _evidence_mapping(bundle["exact_input_context"], "A18 exact input context")
    if context.get("schema_version") != INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION:
        raise ValueError("A18 exact input context schema changed")
    certification = _evidence_mapping(
        bundle["full_coupled_pore_certification"],
        "A18 coupled-pore certification",
    )
    stored_bit_re_evaluation = _evidence_mapping(
        bundle["stored_bit_log_product_re_evaluation"],
        "A18 stored-bit log/product re-evaluation",
    )
    required_audit = {
        "caller_lower_y_hexane",
        "caller_upper_y_hexane",
        "final_lower_y_hexane",
        "final_upper_y_hexane",
        "returned_y_hexane",
        "returned_log_fugacity_residual",
        "returned_product_hexane_activity",
        "actual_oil_fraction_label",
        "log_fugacity_tolerance",
        "rounding_plateau_encountered",
        "no_distinct_midpoint_encountered",
        "bisection_iteration_count",
        "final_upper_classification",
        "final_upper_log_fugacity_residual",
        "final_upper_product_hexane_activity",
        "full_coupled_pore_certification_passed",
        "numerical_amendment_id",
        "postprojection_or_nextafter_used",
        "log_fugacity_equation_changed",
        "log_fugacity_tolerance_changed",
        "physically_qualifying",
    }
    if not required_audit <= set(audit):
        raise ValueError("A18 root audit is incomplete")
    if not (
        audit["numerical_amendment_id"] == INTERFACE_ROOT_AMENDMENT_ID
        and audit["full_coupled_pore_certification_passed"] is True
        and audit["postprojection_or_nextafter_used"] is False
        and audit["log_fugacity_equation_changed"] is False
        and audit["log_fugacity_tolerance_changed"] is False
        and audit["physically_qualifying"] is False
        and audit["log_fugacity_tolerance"] == INTERFACE_ROOT_LOG_FUGACITY_LIMIT
    ):
        raise RuntimeError("A18 root changed an equation, tolerance, or endpoint")
    lower = float(audit["caller_lower_y_hexane"])
    upper = float(audit["caller_upper_y_hexane"])
    final_lower = float(audit["final_lower_y_hexane"])
    returned = float(audit["returned_y_hexane"])
    final_upper = float(audit["final_upper_y_hexane"])
    residual = float(audit["returned_log_fugacity_residual"])
    activity = float(audit["returned_product_hexane_activity"])
    bisection_iteration_count = audit["bisection_iteration_count"]
    if (
        isinstance(bisection_iteration_count, bool)
        or not isinstance(bisection_iteration_count, int)
        or not 0 <= bisection_iteration_count <= 120
    ):
        raise RuntimeError("A18 bisection contraction count is invalid")
    if not (
        0.0 < lower < upper < 1.0
        and lower <= final_lower <= returned <= final_upper <= upper
        and 0.0 <= residual <= INTERFACE_ROOT_LOG_FUGACITY_LIMIT
        and 0.0 < activity <= 1.0
    ):
        raise RuntimeError("A18 returned root is outside its frozen certificate")
    plateau = audit["rounding_plateau_encountered"] is True
    no_midpoint = audit["no_distinct_midpoint_encountered"] is True
    if plateau != no_midpoint or (
        plateau and math.nextafter(final_lower, math.inf) != final_upper
    ):
        raise RuntimeError("A18 rounding-plateau certificate is inconsistent")
    upper_classification = audit["final_upper_classification"]
    allowed_upper_classifications = {
        "negative_log_fugacity_residual_liquid_side",
        "product_activity_supersaturated_on_log_tolerance_plateau",
        "returned_fully_certified_upper_endpoint",
    }
    if upper_classification not in allowed_upper_classifications:
        raise RuntimeError("A18 final upper-endpoint classification is not live")
    upper_residual = float(audit["final_upper_log_fugacity_residual"])
    upper_activity = audit["final_upper_product_hexane_activity"]
    if upper_classification == "negative_log_fugacity_residual_liquid_side":
        upper_consistent = upper_residual < 0.0 and upper_activity is None
    elif upper_classification == (
        "product_activity_supersaturated_on_log_tolerance_plateau"
    ):
        upper_consistent = (
            0.0 <= upper_residual <= INTERFACE_ROOT_LOG_FUGACITY_LIMIT
            and isinstance(upper_activity, (int, float))
            and math.isfinite(float(upper_activity))
            and float(upper_activity) > 1.0
        )
    else:
        upper_consistent = (
            returned == final_upper
            and residual == upper_residual
            and activity == upper_activity
        )
    if not upper_consistent:
        raise RuntimeError("A18 final upper-endpoint classification is inconsistent")
    for field, expected in (
        ("interface_temperature_t_gamma_k", float(context["interface_temperature_t_gamma_k"])),
        ("pressure_pa", float(context["pressure_pa"])),
        ("binary_gas_interaction_k_wh", float(context["binary_gas_interaction_k_wh"])),
        ("actual_oil_fraction_label", float(audit["actual_oil_fraction_label"])),
        ("caller_lower_y_hexane", lower),
        ("caller_upper_y_hexane", upper),
        ("log_fugacity_tolerance", INTERFACE_ROOT_LOG_FUGACITY_LIMIT),
        ("returned_y_hexane", returned),
    ):
        _require_binary64_hex(context, field, expected)
    if (
        context.get("final_lower_y_hexane_binary64_hex") != final_lower.hex()
        or context.get("final_upper_y_hexane_binary64_hex") != final_upper.hex()
    ):
        raise RuntimeError("A18 final bracket lost binary64 identity")
    pore = _evidence_mapping(certification.get("coupled_pore_params"), "coupled-pore params")
    state = _evidence_mapping(
        certification.get("equilibrium_pore_state"),
        "equilibrium pore state",
    )
    if certification.get("schema_version") != INTERFACE_ROOT_CERTIFICATION_SCHEMA_VERSION:
        raise ValueError("A18 coupled-pore certificate schema changed")
    if certification.get("coupled_pore_params_canonical_sha256") != (
        _canonical_evidence_sha256(pore)
    ) or certification.get("equilibrium_pore_state_canonical_sha256") != (
        _canonical_evidence_sha256(state)
    ):
        raise RuntimeError("A18 coupled-pore certificate hash mismatch")
    if not (
        pore.get("w_o") == audit["actual_oil_fraction_label"]
        and pore.get("k_wh") == context["binary_gas_interaction_k_wh"]
        and state.get("temperature_k") == context["interface_temperature_t_gamma_k"]
        and state.get("pressure_pa") == context["pressure_pa"]
        and state.get("y_hexane") == returned
        and state.get("hexane_activity") == activity
    ):
        raise RuntimeError("A18 exact input and coupled-pore certificate diverged")
    if stored_bit_re_evaluation.get("schema_version") != 1:
        raise ValueError("A18 stored-bit re-evaluation schema changed")
    for field, expected in (
        ("interface_temperature_t_gamma_k", float(context["interface_temperature_t_gamma_k"])),
        ("pressure_pa", float(context["pressure_pa"])),
        ("returned_y_hexane", returned),
        ("binary_gas_interaction_k_wh", float(context["binary_gas_interaction_k_wh"])),
        ("actual_oil_fraction_label", float(audit["actual_oil_fraction_label"])),
    ):
        _require_binary64_hex(stored_bit_re_evaluation, field, expected)
    re_evaluation_pore = cp.CoupledPoreParams(
        w_o=float(stored_bit_re_evaluation["actual_oil_fraction_label"]),
        k_wh=float(stored_bit_re_evaluation["binary_gas_interaction_k_wh"]),
    )
    re_evaluated_log = front.log_fugacity_residual(
        float(stored_bit_re_evaluation["interface_temperature_t_gamma_k"]),
        float(stored_bit_re_evaluation["pressure_pa"]),
        1.0 - float(stored_bit_re_evaluation["returned_y_hexane"]),
        float(stored_bit_re_evaluation["returned_y_hexane"]),
        k_wh=re_evaluation_pore.k_wh,
    )[0]
    _, re_evaluated_product_activity = cp._gas_component_activities(  # noqa: SLF001
        float(stored_bit_re_evaluation["interface_temperature_t_gamma_k"]),
        float(stored_bit_re_evaluation["pressure_pa"]),
        float(stored_bit_re_evaluation["returned_y_hexane"]),
        re_evaluation_pore,
    )
    if not (
        stored_bit_re_evaluation.get("recomputed_log_fugacity_residual")
        == re_evaluated_log
        == residual
        and stored_bit_re_evaluation.get(
            "recomputed_log_fugacity_residual_binary64_hex"
        )
        == residual.hex()
        and stored_bit_re_evaluation.get("recomputed_product_hexane_activity")
        == re_evaluated_product_activity
        == activity
        and stored_bit_re_evaluation.get(
            "recomputed_product_hexane_activity_binary64_hex"
        )
        == activity.hex()
        and stored_bit_re_evaluation.get("log_fugacity_residual_bit_exact") is True
        and stored_bit_re_evaluation.get("product_hexane_activity_bit_exact") is True
        and stored_bit_re_evaluation.get("full_coupled_pore_state_bit_exact") is True
        and stored_bit_re_evaluation.get("tolerance_or_approximate_comparison_used")
        is False
    ):
        raise RuntimeError("A18 stored-bit log/product re-evaluation is not bit exact")
    binding_payload = {
        "accepted_balance_root": accepted_balance_root,
        "root_audit": audit,
        "exact_input_context": context,
        "full_coupled_pore_certification": certification,
        "stored_bit_log_product_re_evaluation": stored_bit_re_evaluation,
    }
    binding = _canonical_evidence_sha256(binding_payload)
    if bundle["canonical_binding_sha256"] != binding:
        raise RuntimeError("A18 canonical root binding hash mismatch")
    return {
        "canonical_binding_sha256": binding,
        "accepted_balance_root": accepted_balance_root,
        "rounding_plateau_encountered": plateau,
        "bisection_iteration_count": bisection_iteration_count,
        "full_coupled_pore_certification_passed": True,
        "tolerance_changed": False,
        "clipping_projection_or_nextafter_used": False,
    }


def _interface_composition_root_audit_summary(
    samples: Sequence[Mapping[str, object]],
) -> dict:
    """Validate and identify every canonical accepted-root A18 bundle."""

    if not samples:
        raise ValueError("A18 step summary needs at least one accepted root")
    roots = []
    plateau_count = 0
    total_bisection_iterations = 0
    for sample in samples:
        bundle = _evidence_mapping(
            sample.get("interface_saturation_endpoint_audit"),
            "accepted-root A18 bundle",
        )
        validation = _validate_interface_saturation_endpoint_audit(bundle)
        accepted_balance_root = validation["accepted_balance_root"]
        if not (
            accepted_balance_root["root_role"] == sample["root_role"]
            and accepted_balance_root["balance_interval_dt_s"]
            == sample["balance_interval_dt_s"]
            and accepted_balance_root["time_s"] == sample["time_s"]
        ):
            raise RuntimeError("A18 bundle and accepted-root envelope diverged")
        plateau_count += int(validation["rounding_plateau_encountered"] is True)
        total_bisection_iterations += int(validation["bisection_iteration_count"])
        identity_payload = {
            "root_role": sample["root_role"],
            "balance_interval_dt_s": sample["balance_interval_dt_s"],
            "time_s": sample["time_s"],
            "canonical_binding_sha256": validation["canonical_binding_sha256"],
            "bisection_iteration_count": validation[
                "bisection_iteration_count"
            ],
        }
        roots.append(
            {
                **identity_payload,
                "accepted_root_identity_sha256": _canonical_evidence_sha256(
                    identity_payload
                ),
            }
        )
    root_identities = [root["accepted_root_identity_sha256"] for root in roots]
    if len(set(root_identities)) != len(root_identities):
        raise RuntimeError("A18 accepted-root identities are not one-to-one")
    return {
        "schema_version": INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": INTERFACE_ROOT_AMENDMENT_ID,
        "accepted_root_count": len(roots),
        "accepted_root_identity_count": len(root_identities),
        "rounding_plateau_root_count": plateau_count,
        "rounding_plateau_encountered": plateau_count > 0,
        "total_bisection_iteration_count": total_bisection_iterations,
        "maximum_bisection_iteration_count": max(
            int(root["bisection_iteration_count"]) for root in roots
        ),
        "full_coupled_pore_certified_root_count": len(roots),
        "tolerance_changed_root_count": 0,
        "clipping_projection_or_nextafter_root_count": 0,
        "all_roots_full_coupled_pore_certification_passed": True,
        "all_roots_use_unchanged_log_fugacity_tolerance": True,
        "all_roots_avoid_clipping_projection_and_nextafter": True,
        "accepted_root_identities": roots,
        "contract_passed": True,
        "physically_qualifying": False,
    }


def _expected_dry_capacity_audit_count(
    evaluation_kind: str,
    dry_energy_row_count: int,
) -> int:
    """Return the exact number of A19 primitive evaluations for one scaler."""

    if dry_energy_row_count < 1:
        raise ValueError("A19 residual scales need at least one dry-energy row")
    if evaluation_kind == "same_cell_original_certification":
        return dry_energy_row_count
    if evaluation_kind == "same_cell_face_limit_optimization":
        return 2 * dry_energy_row_count
    if evaluation_kind == "birth_original_certification":
        return dry_energy_row_count
    if evaluation_kind == "strict_face_departure_certification":
        return 2 * dry_energy_row_count - 1
    if evaluation_kind == "exact_face_arrival_certification":
        return 2 * dry_energy_row_count
    raise ValueError(f"unknown A19 residual-scale evaluation kind: {evaluation_kind}")


def _validate_dry_energy_capacity_stencil_audit(
    value: Mapping[str, object],
    *,
    expected_pressure_pa: float,
    expected_primitive_domain_kind: str,
    conditioned_temperature_bounds_k: tuple[float, float],
) -> None:
    """Validate one serialized A19 derivative certificate without tolerance."""

    audit = _evidence_mapping(value, "A19 dry-energy stencil audit")
    required = {
        "temperature_k",
        "pressure_pa",
        "y_hexane",
        "actual_oil_fraction_label",
        "gas_accessible_fraction",
        "primitive_domain_kind",
        "stencil_kind",
        "delta_temperature_k",
        "sample_temperatures_k",
        "sample_full_coupled_pore_certified",
        "center_gas_only_admissible",
        "central_lower_gas_only_admissible",
        "central_upper_gas_only_admissible",
        "forward_second_gas_only_admissible",
        "backward_second_gas_only_admissible",
        "capacity_j_m3_k",
        "numerical_amendment_id",
        "normalization_only",
        "accepted_residual_equations_changed",
        "composition_or_activity_clipping_used",
        "property_extrapolation_used",
        "physically_qualifying",
    }
    if not required <= set(audit):
        raise ValueError("A19 dry-energy stencil audit is incomplete")
    if not (
        audit["numerical_amendment_id"] == DRY_ENERGY_CAPACITY_AMENDMENT_ID
        and audit["normalization_only"] is True
        and audit["accepted_residual_equations_changed"] is False
        and audit["composition_or_activity_clipping_used"] is False
        and audit["property_extrapolation_used"] is False
        and audit["physically_qualifying"] is False
    ):
        raise RuntimeError("A19 stencil changed equations or used an unauthorized repair")
    temperature = float(audit["temperature_k"])
    pressure = float(audit["pressure_pa"])
    y_hexane = float(audit["y_hexane"])
    oil_label = float(audit["actual_oil_fraction_label"])
    delta = float(audit["delta_temperature_k"])
    capacity = float(audit["capacity_j_m3_k"])
    if not all(
        math.isfinite(item)
        for item in (temperature, pressure, y_hexane, oil_label, delta, capacity)
    ) or not (
        pressure == expected_pressure_pa
        and 0.0 < y_hexane < 1.0
        and audit["gas_accessible_fraction"] == 1.0
        and delta > 0.0
        and capacity > 0.0
    ):
        raise RuntimeError("A19 stencil primitive or capacity is invalid")
    if audit["primitive_domain_kind"] != expected_primitive_domain_kind:
        raise RuntimeError("A19 stencil changed the declared primitive domain")
    lower, upper = conditioned_temperature_bounds_k
    expected_delta = min(
        1.0e-3,
        0.25 * (temperature - lower),
        0.25 * (upper - temperature),
    )
    if delta != expected_delta:
        raise RuntimeError("A19 stencil changed the deterministic temperature increment")
    sample_temperatures = tuple(audit["sample_temperatures_k"])
    sample_certifications = tuple(audit["sample_full_coupled_pore_certified"])
    if not sample_temperatures or not all(
        isinstance(item, (int, float)) and math.isfinite(float(item))
        for item in sample_temperatures
    ) or len(sample_certifications) != len(sample_temperatures) or not all(
        item is True for item in sample_certifications
    ):
        raise RuntimeError("A19 stencil samples lack full coupled-pore certification")
    kind = audit["stencil_kind"]
    centered_samples = (temperature + delta, temperature - delta)
    if kind == "joint_primitive_legacy_centered_unchanged":
        if not (
            expected_primitive_domain_kind == "JointPrimitiveBand"
            and sample_temperatures == centered_samples
            and all(
                audit[field] is None
                for field in (
                    "center_gas_only_admissible",
                    "central_lower_gas_only_admissible",
                    "central_upper_gas_only_admissible",
                    "forward_second_gas_only_admissible",
                    "backward_second_gas_only_admissible",
                )
            )
        ):
            raise RuntimeError("A19 legacy centered route changed")
        return
    if expected_primitive_domain_kind != "GasOnlyPrimitiveDomain" or (
        audit["center_gas_only_admissible"] is not True
    ):
        raise RuntimeError("A19 gas-only stencil lacks its exact center")
    colder_ok = audit["central_lower_gas_only_admissible"]
    hotter_ok = audit["central_upper_gas_only_admissible"]
    forward_second_ok = audit["forward_second_gas_only_admissible"]
    backward_second_ok = audit["backward_second_gas_only_admissible"]
    if not isinstance(colder_ok, bool) or not isinstance(hotter_ok, bool):
        raise RuntimeError("A19 gas-only stencil lacks centered topology decisions")
    if kind == "gas_only_second_order_centered":
        valid = (
            colder_ok is True
            and hotter_ok is True
            and forward_second_ok is None
            and backward_second_ok is None
            and sample_temperatures == centered_samples
        )
    elif kind == "gas_only_second_order_forward":
        valid = (
            colder_ok is False
            and hotter_ok is True
            and forward_second_ok is True
            and backward_second_ok is None
            and sample_temperatures
            == (temperature, temperature + delta, temperature + 2.0 * delta)
        )
    elif kind == "gas_only_second_order_backward":
        valid = (
            colder_ok is True
            and not (hotter_ok is True and forward_second_ok is True)
            and isinstance(forward_second_ok, bool)
            and backward_second_ok is True
            and sample_temperatures
            == (temperature, temperature - delta, temperature - 2.0 * delta)
        )
    else:
        raise RuntimeError("A19 stencil kind is unknown")
    if not valid:
        raise RuntimeError("A19 stencil selection order or sample identity changed")


def _dry_energy_capacity_scale_evaluation_evidence(
    root: _AcceptedBalanceIntervalRoot,
    config: cut.CutTransportConfig,
    evaluation_kind: str,
    scales: object,
) -> dict:
    """Serialize one retained A19 residual-scale evaluation."""

    audits = getattr(scales, "dry_energy_capacity_stencil_audits", None)
    if not isinstance(audits, tuple) or not all(
        isinstance(audit, ci.DryEnergyCapacityStencilAudit) for audit in audits
    ):
        raise TypeError("accepted residual scales lost their A19 stencil audits")
    dry_energy_scales = tuple(getattr(scales, "dry_energy_w", ()))
    residual_scale_vector = tuple(getattr(scales, "vector", ()))
    if not dry_energy_scales or not residual_scale_vector or not all(
        math.isfinite(float(value)) and float(value) > 0.0
        for value in residual_scale_vector
    ):
        raise RuntimeError("A19 residual-scale evidence lost its positive scale vector")
    expected_count = _expected_dry_capacity_audit_count(
        evaluation_kind,
        len(dry_energy_scales),
    )
    if len(audits) != expected_count:
        raise RuntimeError(
            "accepted residual-scale evaluation has an unexpected A19 audit count"
        )
    primitive_domain_kind = type(config.dry.primitive_band).__name__
    temperature_bounds = config.dry.conditioned_temperature_domain.solver_bounds_k
    audit_records = []
    for index, audit in enumerate(audits):
        serialized = asdict(audit)
        _validate_dry_energy_capacity_stencil_audit(
            serialized,
            expected_pressure_pa=config.dry.pressure_pa,
            expected_primitive_domain_kind=primitive_domain_kind,
            conditioned_temperature_bounds_k=temperature_bounds,
        )
        audit_records.append(
            {
                "capacity_evaluation_index": index,
                "audit": serialized,
                "audit_canonical_sha256": _canonical_evidence_sha256(serialized),
            }
        )
    vector_payload = {"residual_scale_vector": list(residual_scale_vector)}
    vector_hash = _canonical_evidence_sha256(vector_payload)
    identity_payload = {
        "root_role": root.role,
        "balance_interval_dt_s": root.interval_dt_s,
        "time_s": root.time_s,
        "residual_scale_evaluation_kind": evaluation_kind,
        "residual_scale_vector_canonical_sha256": vector_hash,
        "audit_canonical_sha256": [
            record["audit_canonical_sha256"] for record in audit_records
        ],
    }
    bundle = {
        "schema_version": DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        "root_role": root.role,
        "balance_interval_dt_s": root.interval_dt_s,
        "time_s": root.time_s,
        "residual_scale_evaluation_kind": evaluation_kind,
        "residual_scale_evaluation_count": 1,
        "dry_energy_row_count": len(dry_energy_scales),
        "expected_capacity_evaluation_count": expected_count,
        "capacity_evaluation_count": len(audit_records),
        "primitive_domain_kind": primitive_domain_kind,
        "pressure_pa": config.dry.pressure_pa,
        "conditioned_temperature_solver_bounds_k": list(temperature_bounds),
        "dry_energy_scales_w": list(dry_energy_scales),
        "rh_energy_scale_w": float(getattr(scales, "rh_energy_w")),
        **vector_payload,
        "residual_scale_vector_canonical_sha256": vector_hash,
        "capacity_evaluation_audits": audit_records,
        "residual_scale_evaluation_identity_sha256": (
            _canonical_evidence_sha256(identity_payload)
        ),
        "all_capacity_evaluations_finite_positive_and_fully_certified": True,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "physically_qualifying": False,
    }
    _validate_dry_energy_capacity_evaluation_evidence(bundle)
    return bundle


def _validate_dry_energy_capacity_evaluation_evidence(
    value: Mapping[str, object],
) -> dict[str, object]:
    """Fail closed on one serialized accepted-root A19 scale evaluation."""

    bundle = _evidence_mapping(value, "A19 residual-scale evaluation evidence")
    required = {
        "schema_version",
        "numerical_amendment_id",
        "root_role",
        "balance_interval_dt_s",
        "time_s",
        "residual_scale_evaluation_kind",
        "residual_scale_evaluation_count",
        "dry_energy_row_count",
        "expected_capacity_evaluation_count",
        "capacity_evaluation_count",
        "primitive_domain_kind",
        "pressure_pa",
        "conditioned_temperature_solver_bounds_k",
        "dry_energy_scales_w",
        "rh_energy_scale_w",
        "residual_scale_vector",
        "residual_scale_vector_canonical_sha256",
        "capacity_evaluation_audits",
        "residual_scale_evaluation_identity_sha256",
        "all_capacity_evaluations_finite_positive_and_fully_certified",
        "accepted_residual_equations_changed",
        "composition_or_activity_clipping_used",
        "property_extrapolation_used",
        "acceptance_tolerances_changed",
        "physically_qualifying",
    }
    if not required <= set(bundle):
        raise ValueError("A19 residual-scale evaluation evidence is incomplete")
    if not (
        bundle["schema_version"] == DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION
        and bundle["numerical_amendment_id"] == DRY_ENERGY_CAPACITY_AMENDMENT_ID
        and bundle["residual_scale_evaluation_count"] == 1
        and bundle["all_capacity_evaluations_finite_positive_and_fully_certified"] is True
        and bundle["accepted_residual_equations_changed"] is False
        and bundle["composition_or_activity_clipping_used"] is False
        and bundle["property_extrapolation_used"] is False
        and bundle["acceptance_tolerances_changed"] is False
        and bundle["physically_qualifying"] is False
    ):
        raise RuntimeError("A19 residual-scale evidence changed its frozen authority")
    row_count = int(bundle["dry_energy_row_count"])
    expected_count = _expected_dry_capacity_audit_count(
        str(bundle["residual_scale_evaluation_kind"]),
        row_count,
    )
    records = bundle["capacity_evaluation_audits"]
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        raise TypeError("A19 capacity evaluations must be a sequence")
    if not (
        bundle["expected_capacity_evaluation_count"] == expected_count
        and bundle["capacity_evaluation_count"] == expected_count
        and len(records) == expected_count
        and len(bundle["dry_energy_scales_w"]) == row_count
    ):
        raise RuntimeError("A19 residual-scale row/audit counts do not match")
    pressure = float(bundle["pressure_pa"])
    bounds = tuple(bundle["conditioned_temperature_solver_bounds_k"])
    if len(bounds) != 2:
        raise ValueError("A19 conditioned temperature bounds are incomplete")
    audit_hashes = []
    stencil_kinds: dict[str, int] = {}
    for index, raw_record in enumerate(records):
        record = _evidence_mapping(raw_record, "A19 capacity evaluation record")
        audit = _evidence_mapping(record.get("audit"), "A19 capacity stencil audit")
        if record.get("capacity_evaluation_index") != index:
            raise RuntimeError("A19 capacity evaluation indices are not contiguous")
        audit_hash = _canonical_evidence_sha256(audit)
        if record.get("audit_canonical_sha256") != audit_hash:
            raise RuntimeError("A19 capacity evaluation hash mismatch")
        _validate_dry_energy_capacity_stencil_audit(
            audit,
            expected_pressure_pa=pressure,
            expected_primitive_domain_kind=str(bundle["primitive_domain_kind"]),
            conditioned_temperature_bounds_k=(float(bounds[0]), float(bounds[1])),
        )
        audit_hashes.append(audit_hash)
        kind = str(audit["stencil_kind"])
        stencil_kinds[kind] = stencil_kinds.get(kind, 0) + 1
    vector = tuple(bundle["residual_scale_vector"])
    if not vector or not all(
        isinstance(item, (int, float)) and math.isfinite(float(item)) and float(item) > 0.0
        for item in vector
    ) or not all(
        isinstance(item, (int, float)) and math.isfinite(float(item)) and float(item) > 0.0
        for item in (*tuple(bundle["dry_energy_scales_w"]), bundle["rh_energy_scale_w"])
    ):
        raise RuntimeError("A19 serialized residual scales are not positive and finite")
    vector_hash = _canonical_evidence_sha256({"residual_scale_vector": list(vector)})
    if bundle["residual_scale_vector_canonical_sha256"] != vector_hash:
        raise RuntimeError("A19 residual-scale vector hash mismatch")
    identity_payload = {
        "root_role": bundle["root_role"],
        "balance_interval_dt_s": bundle["balance_interval_dt_s"],
        "time_s": bundle["time_s"],
        "residual_scale_evaluation_kind": bundle["residual_scale_evaluation_kind"],
        "residual_scale_vector_canonical_sha256": vector_hash,
        "audit_canonical_sha256": audit_hashes,
    }
    evaluation_identity = _canonical_evidence_sha256(identity_payload)
    if bundle["residual_scale_evaluation_identity_sha256"] != evaluation_identity:
        raise RuntimeError("A19 residual-scale evaluation identity hash mismatch")
    return {
        "residual_scale_evaluation_identity_sha256": evaluation_identity,
        "capacity_evaluation_count": expected_count,
        "dry_energy_row_count": row_count,
        "stencil_kind_counts": stencil_kinds,
        "repairs_or_tolerance_changes_used": False,
        "scale_authority_sha256": _canonical_evidence_sha256(identity_payload),
    }


def _dry_energy_capacity_stencil_evidence(
    root: _AcceptedBalanceIntervalRoot,
    config: cut.CutTransportConfig,
) -> dict:
    """Serialize every retained A19 scale evaluation for one accepted root."""

    evaluations = [
        _dry_energy_capacity_scale_evaluation_evidence(
            root,
            config,
            evaluation_kind,
            scales,
        )
        for evaluation_kind, scales in root.residual_scale_evaluations
    ]
    evaluation_identities = [
        evaluation["residual_scale_evaluation_identity_sha256"]
        for evaluation in evaluations
    ]
    identity_payload = {
        "root_role": root.role,
        "balance_interval_dt_s": root.interval_dt_s,
        "time_s": root.time_s,
        "residual_scale_evaluation_identity_sha256": evaluation_identities,
    }
    bundle = {
        "schema_version": DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        "root_role": root.role,
        "balance_interval_dt_s": root.interval_dt_s,
        "time_s": root.time_s,
        "residual_scale_evaluation_count": len(evaluations),
        "residual_scale_evaluations": evaluations,
        "dry_energy_row_count": sum(
            int(evaluation["dry_energy_row_count"]) for evaluation in evaluations
        ),
        "capacity_evaluation_audit_count": sum(
            int(evaluation["capacity_evaluation_count"])
            for evaluation in evaluations
        ),
        "residual_scale_evaluation_kind_counts": {
            kind: sum(
                evaluation["residual_scale_evaluation_kind"] == kind
                for evaluation in evaluations
            )
            for kind in sorted(
                {
                    str(evaluation["residual_scale_evaluation_kind"])
                    for evaluation in evaluations
                }
            )
        },
        "accepted_root_scale_identity_sha256": _canonical_evidence_sha256(
            identity_payload
        ),
        "all_capacity_evaluations_finite_positive_and_fully_certified": True,
        "all_used_residual_scale_evaluations_retained": True,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "physically_qualifying": False,
    }
    _validate_dry_energy_capacity_root_evidence(bundle)
    return {"dry_energy_capacity_stencil_audit": bundle}


def _validate_dry_energy_capacity_root_evidence(
    value: Mapping[str, object],
) -> dict[str, object]:
    """Fail closed unless a root retains every scale basis that accepted it."""

    bundle = _evidence_mapping(value, "accepted-root A19 evidence")
    required = {
        "schema_version",
        "numerical_amendment_id",
        "root_role",
        "balance_interval_dt_s",
        "time_s",
        "residual_scale_evaluation_count",
        "residual_scale_evaluations",
        "dry_energy_row_count",
        "capacity_evaluation_audit_count",
        "residual_scale_evaluation_kind_counts",
        "accepted_root_scale_identity_sha256",
        "all_capacity_evaluations_finite_positive_and_fully_certified",
        "all_used_residual_scale_evaluations_retained",
        "accepted_residual_equations_changed",
        "composition_or_activity_clipping_used",
        "property_extrapolation_used",
        "acceptance_tolerances_changed",
        "physically_qualifying",
    }
    if not required <= set(bundle):
        raise ValueError("accepted-root A19 evidence is incomplete")
    if not (
        bundle["schema_version"] == DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION
        and bundle["numerical_amendment_id"] == DRY_ENERGY_CAPACITY_AMENDMENT_ID
        and bundle["all_capacity_evaluations_finite_positive_and_fully_certified"] is True
        and bundle["all_used_residual_scale_evaluations_retained"] is True
        and bundle["accepted_residual_equations_changed"] is False
        and bundle["composition_or_activity_clipping_used"] is False
        and bundle["property_extrapolation_used"] is False
        and bundle["acceptance_tolerances_changed"] is False
        and bundle["physically_qualifying"] is False
    ):
        raise RuntimeError("accepted-root A19 evidence changed its frozen authority")
    evaluations = bundle["residual_scale_evaluations"]
    if not isinstance(evaluations, Sequence) or isinstance(evaluations, (str, bytes)):
        raise TypeError("accepted-root A19 scale evaluations must be a sequence")
    evaluation_count = int(bundle["residual_scale_evaluation_count"])
    if evaluation_count not in (1, 2) or len(evaluations) != evaluation_count:
        raise RuntimeError("accepted-root A19 scale-evaluation count is invalid")
    validations = []
    kinds = []
    for evaluation in evaluations:
        evaluation_map = _evidence_mapping(evaluation, "A19 scale evaluation")
        validation = _validate_dry_energy_capacity_evaluation_evidence(evaluation_map)
        if not (
            evaluation_map["root_role"] == bundle["root_role"]
            and evaluation_map["balance_interval_dt_s"]
            == bundle["balance_interval_dt_s"]
            and evaluation_map["time_s"] == bundle["time_s"]
        ):
            raise RuntimeError("A19 scale evaluation changed its accepted-root envelope")
        validations.append(validation)
        kinds.append(str(evaluation_map["residual_scale_evaluation_kind"]))
    if len(set(kinds)) != len(kinds):
        raise RuntimeError("accepted-root A19 scale-evaluation kinds are duplicated")
    role = bundle["root_role"]
    if role in {"same_cell_root", "adaptive_same_cell_leaf"}:
        expected_kinds = ["same_cell_original_certification"]
        if evaluation_count == 2:
            expected_kinds.append("same_cell_face_limit_optimization")
    elif role == "exact_face_arrival":
        expected_kinds = ["exact_face_arrival_certification"]
    elif role == "exact_face_departure_remainder":
        expected_kinds = ["strict_face_departure_certification"]
    elif role == "receding_birth":
        expected_kinds = ["birth_original_certification"]
    else:
        raise RuntimeError("accepted-root A19 role is unknown")
    if kinds != expected_kinds:
        raise RuntimeError("accepted-root A19 scale-basis identities are incomplete")
    kind_counts = {kind: kinds.count(kind) for kind in sorted(set(kinds))}
    if not (
        bundle["dry_energy_row_count"]
        == sum(int(item["dry_energy_row_count"]) for item in validations)
        and bundle["capacity_evaluation_audit_count"]
        == sum(int(item["capacity_evaluation_count"]) for item in validations)
        and bundle["residual_scale_evaluation_kind_counts"] == kind_counts
    ):
        raise RuntimeError("accepted-root A19 aggregate counts do not match evaluations")
    evaluation_identities = [
        item["residual_scale_evaluation_identity_sha256"] for item in validations
    ]
    identity_payload = {
        "root_role": bundle["root_role"],
        "balance_interval_dt_s": bundle["balance_interval_dt_s"],
        "time_s": bundle["time_s"],
        "residual_scale_evaluation_identity_sha256": evaluation_identities,
    }
    root_identity = _canonical_evidence_sha256(identity_payload)
    if bundle["accepted_root_scale_identity_sha256"] != root_identity:
        raise RuntimeError("accepted-root A19 identity hash mismatch")
    stencil_kind_counts: dict[str, int] = {}
    for item in validations:
        for kind, count in item["stencil_kind_counts"].items():
            stencil_kind_counts[kind] = stencil_kind_counts.get(kind, 0) + int(count)
    return {
        "accepted_root_scale_identity_sha256": root_identity,
        "residual_scale_evaluation_count": evaluation_count,
        "capacity_evaluation_count": sum(
            int(item["capacity_evaluation_count"]) for item in validations
        ),
        "dry_energy_row_count": sum(
            int(item["dry_energy_row_count"]) for item in validations
        ),
        "stencil_kind_counts": stencil_kind_counts,
        "residual_scale_evaluation_kind_counts": kind_counts,
        "repairs_or_tolerance_changes_used": False,
    }


def _dry_energy_capacity_stencil_audit_summary(
    samples: Sequence[Mapping[str, object]],
) -> dict:
    """Validate and count every accepted-root A19 scale evaluation."""

    if not samples:
        raise ValueError("A19 step summary needs at least one accepted root")
    identities = []
    residual_scale_evaluation_count = 0
    capacity_count = 0
    dry_energy_row_count = 0
    stencil_kind_counts: dict[str, int] = {}
    evaluation_kind_counts: dict[str, int] = {}
    for sample in samples:
        bundle = _evidence_mapping(
            sample.get("dry_energy_capacity_stencil_audit"),
            "accepted-root A19 evidence",
        )
        validation = _validate_dry_energy_capacity_root_evidence(bundle)
        if not (
            bundle["root_role"] == sample["root_role"]
            and bundle["balance_interval_dt_s"] == sample["balance_interval_dt_s"]
            and bundle["time_s"] == sample["time_s"]
        ):
            raise RuntimeError("A19 bundle and accepted-root envelope diverged")
        identities.append(
            {
                "root_role": sample["root_role"],
                "balance_interval_dt_s": sample["balance_interval_dt_s"],
                "time_s": sample["time_s"],
                "accepted_root_scale_identity_sha256": validation[
                    "accepted_root_scale_identity_sha256"
                ],
            }
        )
        residual_scale_evaluation_count += int(
            validation["residual_scale_evaluation_count"]
        )
        capacity_count += int(validation["capacity_evaluation_count"])
        dry_energy_row_count += int(validation["dry_energy_row_count"])
        for kind, count in validation["stencil_kind_counts"].items():
            stencil_kind_counts[kind] = stencil_kind_counts.get(kind, 0) + int(count)
        for kind, count in validation["residual_scale_evaluation_kind_counts"].items():
            evaluation_kind_counts[kind] = evaluation_kind_counts.get(kind, 0) + int(count)
    identity_hashes = [item["accepted_root_scale_identity_sha256"] for item in identities]
    if len(set(identity_hashes)) != len(identity_hashes):
        raise RuntimeError("A19 accepted-root scale identities are not one-to-one")
    return {
        "schema_version": DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        "accepted_root_count": len(samples),
        "residual_scale_evaluation_count": residual_scale_evaluation_count,
        "dry_energy_row_count": dry_energy_row_count,
        "capacity_evaluation_audit_count": capacity_count,
        "stencil_kind_counts": dict(sorted(stencil_kind_counts.items())),
        "residual_scale_evaluation_kind_counts": dict(
            sorted(evaluation_kind_counts.items())
        ),
        "accepted_root_scale_identities": identities,
        "all_capacity_evaluations_finite_positive_and_fully_certified": True,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": True,
        "physically_qualifying": False,
    }


_DIAGNOSTIC_EXECUTION_CLASSIFICATIONS = {
    "diagnostic_bootstrap_or_radius_seed_execution",
    "noncommitted_determinism_or_restart_replay",
    "noncommitted_provisional_atomic_macro_rollback",
    "noncommitted_supporting_qualification_execution",
}


def _accepted_root_a18_a19_pair_identities(
    samples: Sequence[Mapping[str, object]],
    a18: Mapping[str, object],
    a19: Mapping[str, object],
) -> list[dict[str, object]]:
    """Bind corresponding A18 and A19 identities without requiring global uniqueness."""

    a18_roots = a18.get("accepted_root_identities")
    a19_roots = a19.get("accepted_root_scale_identities")
    if not (
        isinstance(a18_roots, Sequence)
        and not isinstance(a18_roots, (str, bytes))
        and isinstance(a19_roots, Sequence)
        and not isinstance(a19_roots, (str, bytes))
        and len(samples) == len(a18_roots) == len(a19_roots)
    ):
        raise RuntimeError("A18/A19 accepted-root identity counts diverged")
    pairs: list[dict[str, object]] = []
    for root_index, (sample, a18_root, a19_root) in enumerate(
        zip(samples, a18_roots, a19_roots, strict=True)
    ):
        a18_map = _evidence_mapping(a18_root, "A18 accepted-root identity")
        a19_map = _evidence_mapping(a19_root, "A19 accepted-root identity")
        envelope = {
            "root_role": sample.get("root_role"),
            "balance_interval_dt_s": sample.get("balance_interval_dt_s"),
            "time_s": sample.get("time_s"),
        }
        if not all(
            identity.get(name) == value
            for identity in (a18_map, a19_map)
            for name, value in envelope.items()
        ):
            raise RuntimeError("paired A18/A19 roots changed their accepted envelope")
        payload = {
            "accepted_root_index": root_index,
            **envelope,
            "a18_accepted_root_identity_sha256": a18_map.get(
                "accepted_root_identity_sha256"
            ),
            "a19_accepted_root_scale_identity_sha256": a19_map.get(
                "accepted_root_scale_identity_sha256"
            ),
        }
        pairs.append(
            {
                **payload,
                "physical_root_evidence_pair_identity_sha256": (
                    _canonical_evidence_sha256(payload)
                ),
            }
        )
    return pairs


def _validate_amendment_18_19_execution_occurrence(
    value: Mapping[str, object],
) -> dict[str, object]:
    """Revalidate one accepted diagnostic execution and its exact context binding."""

    occurrence = _evidence_mapping(value, "A18/A19 diagnostic execution occurrence")
    required = {
        "schema_version",
        "numerical_amendment_ids",
        "execution_classification",
        "execution_context",
        "execution_context_identity_sha256",
        "accepted_root_count",
        "accepted_root_samples",
        "accepted_root_evidence_pair_count",
        "accepted_root_evidence_pairs",
        "execution_occurrence_identity_sha256",
        "committed_trajectory_root_count",
        "adds_committed_trajectory_root_count",
        "diagnostic_execution_partition",
        "physical_root_identity_uniqueness_required_across_executions",
        "interface_saturation_endpoint_audit_summary",
        "dry_energy_capacity_stencil_audit_summary",
        "contract_passed",
        "physically_qualifying",
    }
    if set(occurrence) != required:
        raise ValueError("A18/A19 diagnostic execution occurrence fields are not exact")
    classification = occurrence["execution_classification"]
    if not (
        occurrence["schema_version"] == 1
        and occurrence["numerical_amendment_ids"]
        == [INTERFACE_ROOT_AMENDMENT_ID, DRY_ENERGY_CAPACITY_AMENDMENT_ID]
        and classification in _DIAGNOSTIC_EXECUTION_CLASSIFICATIONS
        and occurrence["committed_trajectory_root_count"] == 0
        and occurrence["adds_committed_trajectory_root_count"] is False
        and occurrence["diagnostic_execution_partition"] is True
        and occurrence["physical_root_identity_uniqueness_required_across_executions"]
        is False
        and occurrence["contract_passed"] is True
        and occurrence["physically_qualifying"] is False
    ):
        raise RuntimeError("A18/A19 diagnostic execution changed its frozen boundary")
    context = _evidence_mapping(occurrence["execution_context"], "execution context")
    if not (
        isinstance(context.get("execution_partition_id"), str)
        and context["execution_partition_id"].strip()
        and isinstance(context.get("execution_role"), str)
        and context["execution_role"].strip()
        and context.get("diagnostic_execution_partition") is True
        and context.get("adds_committed_trajectory_root_count") is False
        and isinstance(
            context.get("may_reference_committed_trajectory_root_elsewhere"),
            bool,
        )
    ):
        raise ValueError("A18/A19 diagnostic execution context is incomplete")
    context_hash = _canonical_evidence_sha256(context)
    if occurrence["execution_context_identity_sha256"] != context_hash:
        raise RuntimeError("A18/A19 diagnostic execution context hash mismatch")
    samples = occurrence["accepted_root_samples"]
    if not isinstance(samples, Sequence) or isinstance(samples, (str, bytes)) or not samples:
        raise ValueError("A18/A19 diagnostic execution has no accepted-root samples")
    sample_maps = [
        _evidence_mapping(sample, "diagnostic accepted-root sample") for sample in samples
    ]
    expected_sample_fields = {
        "root_role",
        "balance_interval_dt_s",
        "time_s",
        "interface_saturation_endpoint_audit",
        "dry_energy_capacity_stencil_audit",
    }
    if any(set(sample) != expected_sample_fields for sample in sample_maps):
        raise ValueError("diagnostic accepted-root sample fields are not exact")
    a18 = _interface_composition_root_audit_summary(sample_maps)
    a19 = _dry_energy_capacity_stencil_audit_summary(sample_maps)
    pairs = _accepted_root_a18_a19_pair_identities(sample_maps, a18, a19)
    accepted_root_count = len(sample_maps)
    if not (
        occurrence["accepted_root_count"] == accepted_root_count
        and occurrence["accepted_root_evidence_pair_count"] == accepted_root_count
        and occurrence["accepted_root_evidence_pairs"] == pairs
        and occurrence["interface_saturation_endpoint_audit_summary"] == a18
        and occurrence["dry_energy_capacity_stencil_audit_summary"] == a19
    ):
        raise RuntimeError("A18/A19 diagnostic execution child binding is not exact")
    pair_hashes = [
        item["physical_root_evidence_pair_identity_sha256"] for item in pairs
    ]
    occurrence_payload = {
        "execution_classification": classification,
        "execution_context_identity_sha256": context_hash,
        "accepted_root_count": accepted_root_count,
        "physical_root_evidence_pair_identity_sha256": pair_hashes,
    }
    occurrence_hash = _canonical_evidence_sha256(occurrence_payload)
    if occurrence["execution_occurrence_identity_sha256"] != occurrence_hash:
        raise RuntimeError("A18/A19 diagnostic execution occurrence hash mismatch")
    return {
        "execution_context_identity_sha256": context_hash,
        "execution_occurrence_identity_sha256": occurrence_hash,
        "physical_root_evidence_pair_identity_sha256": pair_hashes,
        "accepted_root_count": accepted_root_count,
        "interface_saturation_endpoint_audit_summary": a18,
        "dry_energy_capacity_stencil_audit_summary": a19,
    }


def _amendment_18_19_execution_occurrence(
    step: AcceptedSameCellStep,
    *,
    execution_classification: str,
    execution_context: Mapping[str, object],
) -> dict[str, object]:
    """Serialize every retained A18/A19 scale basis for one accepted execution."""

    if execution_classification not in _DIAGNOSTIC_EXECUTION_CLASSIFICATIONS:
        raise ValueError("unknown A18/A19 diagnostic execution classification")
    context = dict(execution_context)
    context.update(
        diagnostic_execution_partition=True,
        adds_committed_trajectory_root_count=False,
        may_reference_committed_trajectory_root_elsewhere=(
            execution_classification
            == "diagnostic_bootstrap_or_radius_seed_execution"
        ),
    )
    roots = _accepted_balance_interval_roots(step)
    samples = []
    for root in roots:
        serialized = _surface_flux_sample(root, include_constitutive_evidence=False)
        samples.append(
            {
                "root_role": serialized["root_role"],
                "balance_interval_dt_s": serialized["balance_interval_dt_s"],
                "time_s": serialized["time_s"],
                "interface_saturation_endpoint_audit": serialized[
                    "interface_saturation_endpoint_audit"
                ],
                "dry_energy_capacity_stencil_audit": serialized[
                    "dry_energy_capacity_stencil_audit"
                ],
            }
        )
    a18 = _interface_composition_root_audit_summary(samples)
    a19 = _dry_energy_capacity_stencil_audit_summary(samples)
    pairs = _accepted_root_a18_a19_pair_identities(samples, a18, a19)
    context_hash = _canonical_evidence_sha256(context)
    occurrence_payload = {
        "execution_classification": execution_classification,
        "execution_context_identity_sha256": context_hash,
        "accepted_root_count": len(samples),
        "physical_root_evidence_pair_identity_sha256": [
            item["physical_root_evidence_pair_identity_sha256"] for item in pairs
        ],
    }
    occurrence = {
        "schema_version": 1,
        "numerical_amendment_ids": [
            INTERFACE_ROOT_AMENDMENT_ID,
            DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "execution_classification": execution_classification,
        "execution_context": context,
        "execution_context_identity_sha256": context_hash,
        "accepted_root_count": len(samples),
        "accepted_root_samples": samples,
        "accepted_root_evidence_pair_count": len(pairs),
        "accepted_root_evidence_pairs": pairs,
        "execution_occurrence_identity_sha256": _canonical_evidence_sha256(
            occurrence_payload
        ),
        "committed_trajectory_root_count": 0,
        "adds_committed_trajectory_root_count": False,
        "diagnostic_execution_partition": True,
        "physical_root_identity_uniqueness_required_across_executions": False,
        "interface_saturation_endpoint_audit_summary": a18,
        "dry_energy_capacity_stencil_audit_summary": a19,
        "contract_passed": True,
        "physically_qualifying": False,
    }
    _validate_amendment_18_19_execution_occurrence(occurrence)
    return occurrence


def _merge_nonnegative_count_mapping(
    target: dict[str, int],
    source: Mapping[str, object],
) -> None:
    for name, value in source.items():
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"A18/A19 diagnostic count {name!r} is invalid")
        target[name] = target.get(name, 0) + value


def _aggregate_amendment_18_19_diagnostic_partition(
    partition_id: str,
    execution_occurrences: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Aggregate unique execution contexts while allowing repeated physical root bits."""

    if not isinstance(partition_id, str) or not partition_id.strip():
        raise ValueError("A18/A19 diagnostic partition needs an identifier")
    if not execution_occurrences:
        raise ValueError("A18/A19 diagnostic partition needs accepted executions")
    occurrences = [dict(item) for item in execution_occurrences]
    validations = [
        _validate_amendment_18_19_execution_occurrence(item) for item in occurrences
    ]
    if any(
        item["execution_context"]["execution_partition_id"] != partition_id
        for item in occurrences
    ):
        raise RuntimeError("A18/A19 occurrence is bound to the wrong partition")
    context_hashes = [
        item["execution_context_identity_sha256"] for item in validations
    ]
    occurrence_hashes = [
        item["execution_occurrence_identity_sha256"] for item in validations
    ]
    if not (
        len(set(context_hashes)) == len(context_hashes)
        and len(set(occurrence_hashes)) == len(occurrence_hashes)
    ):
        raise RuntimeError("A18/A19 diagnostic execution contexts are not one-to-one")
    physical_hashes = [
        root_hash
        for item in validations
        for root_hash in item["physical_root_evidence_pair_identity_sha256"]
    ]
    accepted_root_count = sum(int(item["accepted_root_count"]) for item in validations)
    plateau_count = 0
    total_bisection_iterations = 0
    maximum_bisection_iterations = 0
    residual_scale_evaluation_count = 0
    dry_energy_row_count = 0
    capacity_evaluation_count = 0
    stencil_kind_counts: dict[str, int] = {}
    evaluation_kind_counts: dict[str, int] = {}
    for item in validations:
        a18 = item["interface_saturation_endpoint_audit_summary"]
        a19 = item["dry_energy_capacity_stencil_audit_summary"]
        plateau_count += int(a18["rounding_plateau_root_count"])
        total_bisection_iterations += int(a18["total_bisection_iteration_count"])
        maximum_bisection_iterations = max(
            maximum_bisection_iterations,
            int(a18["maximum_bisection_iteration_count"]),
        )
        residual_scale_evaluation_count += int(a19["residual_scale_evaluation_count"])
        dry_energy_row_count += int(a19["dry_energy_row_count"])
        capacity_evaluation_count += int(a19["capacity_evaluation_audit_count"])
        _merge_nonnegative_count_mapping(stencil_kind_counts, a19["stencil_kind_counts"])
        _merge_nonnegative_count_mapping(
            evaluation_kind_counts,
            a19["residual_scale_evaluation_kind_counts"],
        )
    a18_accounting = {
        "accepted_diagnostic_root_execution_count": accepted_root_count,
        "full_coupled_pore_certified_root_execution_count": accepted_root_count,
        "rounding_plateau_root_execution_count": plateau_count,
        "total_bisection_iteration_count": total_bisection_iterations,
        "maximum_bisection_iteration_count": maximum_bisection_iterations,
        "tolerance_changed_root_execution_count": 0,
        "clipping_projection_or_nextafter_root_execution_count": 0,
        "contract_passed": True,
    }
    a19_accounting = {
        "accepted_diagnostic_root_execution_count": accepted_root_count,
        "residual_scale_evaluation_count": residual_scale_evaluation_count,
        "dry_energy_row_count": dry_energy_row_count,
        "capacity_evaluation_audit_count": capacity_evaluation_count,
        "stencil_kind_counts": dict(sorted(stencil_kind_counts.items())),
        "residual_scale_evaluation_kind_counts": dict(
            sorted(evaluation_kind_counts.items())
        ),
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": True,
    }
    partition_identity_payload = {
        "partition_id": partition_id,
        "execution_context_identity_sha256": context_hashes,
        "execution_occurrence_identity_sha256": occurrence_hashes,
        "interface_saturation_endpoint_audit_accounting": a18_accounting,
        "dry_energy_capacity_stencil_audit_accounting": a19_accounting,
    }
    return {
        "schema_version": 1,
        "numerical_amendment_ids": [
            INTERFACE_ROOT_AMENDMENT_ID,
            DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "partition_id": partition_id,
        "execution_occurrence_count": len(occurrences),
        "execution_context_identity_count": len(context_hashes),
        "execution_occurrence_identity_count": len(occurrence_hashes),
        "accepted_diagnostic_root_execution_count": accepted_root_count,
        "committed_trajectory_root_count": 0,
        "execution_occurrences": occurrences,
        "execution_context_identity_sha256": context_hashes,
        "execution_occurrence_identity_sha256": occurrence_hashes,
        "physical_root_evidence_pair_identity_occurrence_count": len(physical_hashes),
        "unique_physical_root_evidence_pair_identity_count": len(set(physical_hashes)),
        "physical_root_identity_uniqueness_required_across_executions": False,
        "execution_context_identity_uniqueness_required": True,
        "may_reference_committed_trajectory_root_elsewhere": any(
            occurrence["execution_context"][
                "may_reference_committed_trajectory_root_elsewhere"
            ]
            for occurrence in occurrences
        ),
        "interface_saturation_endpoint_audit_accounting": a18_accounting,
        "dry_energy_capacity_stencil_audit_accounting": a19_accounting,
        "partition_identity_sha256": _canonical_evidence_sha256(
            partition_identity_payload
        ),
        "contract_passed": True,
        "physically_qualifying": False,
    }


def _aggregate_amendment_18_19_diagnostic_partitions(
    partitions: Mapping[str, Mapping[str, object]],
    *,
    expected_partition_ids: Sequence[str],
) -> dict[str, object]:
    """Revalidate exact diagnostic partitions and bind their execution identities."""

    expected_ids = list(expected_partition_ids)
    if list(partitions) != expected_ids or len(set(expected_ids)) != len(expected_ids):
        raise ValueError("A18/A19 diagnostic partition set or order changed")
    validated_partitions: dict[str, dict[str, object]] = {}
    for partition_id in expected_ids:
        value = _evidence_mapping(partitions[partition_id], "A18/A19 diagnostic partition")
        occurrences = value.get("execution_occurrences")
        if not isinstance(occurrences, Sequence) or isinstance(occurrences, (str, bytes)):
            raise TypeError("A18/A19 diagnostic partition occurrences are invalid")
        recomputed = _aggregate_amendment_18_19_diagnostic_partition(
            partition_id,
            occurrences,
        )
        if dict(value) != recomputed:
            raise RuntimeError("A18/A19 diagnostic partition binding is not exact")
        validated_partitions[partition_id] = recomputed
    partition_hashes = [
        validated_partitions[name]["partition_identity_sha256"] for name in expected_ids
    ]
    context_hashes = [
        context_hash
        for name in expected_ids
        for context_hash in validated_partitions[name][
            "execution_context_identity_sha256"
        ]
    ]
    occurrence_hashes = [
        occurrence_hash
        for name in expected_ids
        for occurrence_hash in validated_partitions[name][
            "execution_occurrence_identity_sha256"
        ]
    ]
    if not (
        len(set(partition_hashes)) == len(partition_hashes)
        and len(set(context_hashes)) == len(context_hashes)
        and len(set(occurrence_hashes)) == len(occurrence_hashes)
    ):
        raise RuntimeError("A18/A19 diagnostic partition execution identities collide")
    physical_hashes = [
        pair["physical_root_evidence_pair_identity_sha256"]
        for name in expected_ids
        for occurrence in validated_partitions[name]["execution_occurrences"]
        for pair in occurrence["accepted_root_evidence_pairs"]
    ]
    execution_count = sum(
        int(value["execution_occurrence_count"])
        for value in validated_partitions.values()
    )
    root_count = sum(
        int(value["accepted_diagnostic_root_execution_count"])
        for value in validated_partitions.values()
    )
    a18_accounting = {
        "accepted_diagnostic_root_execution_count": root_count,
        "full_coupled_pore_certified_root_execution_count": sum(
            int(
                value["interface_saturation_endpoint_audit_accounting"][
                    "full_coupled_pore_certified_root_execution_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "rounding_plateau_root_execution_count": sum(
            int(
                value["interface_saturation_endpoint_audit_accounting"][
                    "rounding_plateau_root_execution_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "total_bisection_iteration_count": sum(
            int(
                value["interface_saturation_endpoint_audit_accounting"][
                    "total_bisection_iteration_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "maximum_bisection_iteration_count": max(
            int(
                value["interface_saturation_endpoint_audit_accounting"][
                    "maximum_bisection_iteration_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "tolerance_changed_root_execution_count": 0,
        "clipping_projection_or_nextafter_root_execution_count": 0,
        "contract_passed": True,
    }
    stencil_kind_counts: dict[str, int] = {}
    evaluation_kind_counts: dict[str, int] = {}
    for value in validated_partitions.values():
        a19_partition = value["dry_energy_capacity_stencil_audit_accounting"]
        _merge_nonnegative_count_mapping(
            stencil_kind_counts,
            a19_partition["stencil_kind_counts"],
        )
        _merge_nonnegative_count_mapping(
            evaluation_kind_counts,
            a19_partition["residual_scale_evaluation_kind_counts"],
        )
    a19_accounting = {
        "accepted_diagnostic_root_execution_count": root_count,
        "residual_scale_evaluation_count": sum(
            int(
                value["dry_energy_capacity_stencil_audit_accounting"][
                    "residual_scale_evaluation_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "dry_energy_row_count": sum(
            int(
                value["dry_energy_capacity_stencil_audit_accounting"][
                    "dry_energy_row_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "capacity_evaluation_audit_count": sum(
            int(
                value["dry_energy_capacity_stencil_audit_accounting"][
                    "capacity_evaluation_audit_count"
                ]
            )
            for value in validated_partitions.values()
        ),
        "stencil_kind_counts": dict(sorted(stencil_kind_counts.items())),
        "residual_scale_evaluation_kind_counts": dict(
            sorted(evaluation_kind_counts.items())
        ),
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": True,
    }
    identity_payload = {
        "partition_id": expected_ids,
        "partition_identity_sha256": partition_hashes,
        "execution_context_identity_sha256": context_hashes,
        "execution_occurrence_identity_sha256": occurrence_hashes,
        "accepted_diagnostic_root_execution_count": root_count,
        "interface_saturation_endpoint_audit_accounting": a18_accounting,
        "dry_energy_capacity_stencil_audit_accounting": a19_accounting,
    }
    return {
        "schema_version": 1,
        "numerical_amendment_ids": [
            INTERFACE_ROOT_AMENDMENT_ID,
            DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "partition_count": len(expected_ids),
        "partition_identity_count": len(partition_hashes),
        "partition_ids": expected_ids,
        "partition_identity_sha256": partition_hashes,
        "execution_occurrence_count": execution_count,
        "execution_context_identity_count": len(context_hashes),
        "execution_occurrence_identity_count": len(occurrence_hashes),
        "accepted_diagnostic_root_execution_count": root_count,
        "committed_trajectory_root_count": 0,
        "physical_root_evidence_pair_identity_occurrence_count": len(physical_hashes),
        "unique_physical_root_evidence_pair_identity_count": len(set(physical_hashes)),
        "physical_root_identity_uniqueness_required_across_executions": False,
        "execution_context_identity_uniqueness_required": True,
        "interface_saturation_endpoint_audit_accounting": a18_accounting,
        "dry_energy_capacity_stencil_audit_accounting": a19_accounting,
        "combined_diagnostic_partition_identity_sha256": (
            _canonical_evidence_sha256(identity_payload)
        ),
        "contract_passed": True,
        "physically_qualifying": False,
    }


def _is_r7_wall_job(job) -> bool:
    """Context-free wall-job test for evidence accounting (R7).

    A wall job is a completed=False record whose (scenario, cells, dt) and
    failure signature match a committed R3/R4 wall.  The N12 corner wall is
    recognized by grid+signature: sibling corner jobs at the same grid PASS,
    so only the R4-documented failure matches; any NEW failure elsewhere
    still fails the campaign through the selection gate's exact excluded-set
    accounting, so this predicate cannot hide undocumented failures.
    """

    if not isinstance(job, Mapping) or job.get("completed") is not False:
        return False
    key = (job.get("scenario"), job.get("cells"), job.get("dt_s"))
    error = str(job.get("error") or "")
    signature = R7_WALL_SIGNATURES.get(key)
    if signature is None and key == (SCENARIO_COMPOSITION_ONLY, 12, 0.075):
        signature = "subdivision exhausted its documented minimum dt"
    return bool(signature) and signature in error


def _amendment_18_19_raw_structural_gate(
    jobs: Sequence[Mapping[str, object]],
    diagnostic_partitions: Mapping[str, Mapping[str, object]],
    diagnostic_contract: Mapping[str, object],
) -> dict[str, object]:
    """Bind committed job references separately from diagnostic executions."""

    errors: list[str] = []
    expected_partition_ids = (
        "bootstrap_radius_seed_executions",
        "determinism_restart_replay_executions",
    )
    try:
        expected_diagnostic = _aggregate_amendment_18_19_diagnostic_partitions(
            diagnostic_partitions,
            expected_partition_ids=expected_partition_ids,
        )
        if dict(diagnostic_contract) != expected_diagnostic:
            raise RuntimeError("diagnostic partition aggregate is detached")
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        errors.append(f"noncommitted diagnostic A18/A19 evidence: {exc}")
        expected_diagnostic = None
    job_references: list[dict[str, object]] = []
    committed_root_reference_count = 0
    # R7: verified wall jobs carry no committed evidence by construction and
    # are excluded from evidence accounting on BOTH the producer and the
    # validator side of this shared function.
    wall_job_count = sum(1 for job in jobs if _is_r7_wall_job(job))
    jobs = [job for job in jobs if not _is_r7_wall_job(job)]
    # R4/R7: a campaign whose EVERY job is a verified fail-closed wall (the
    # excluded corner invocation) has no committed evidence children by
    # construction; its verification is the excluded-cell fail-closed
    # accounting, so the committed-evidence obligation is delegated there
    # rather than failing an obligation that cannot exist.
    all_jobs_are_verified_wall_exclusions = wall_job_count > 0 and not jobs
    for job_index, job in enumerate(jobs):
        if not isinstance(job, Mapping) or job.get("completed") is not True:
            errors.append(f"campaign job {job_index} is incomplete")
            continue
        records = job.get("step_records")
        try:
            if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
                raise TypeError("step records are not a sequence")
            expected = _aggregate_amendment_18_19_evidence(records)
            gates = _evidence_mapping(job.get("gates"), "campaign job gates")
            if not (
                job.get("amendment_18_19_evidence_contract") == expected
                and gates.get("interface_saturation_endpoint_audit_contract_passed")
                is True
                and gates.get("dry_energy_capacity_stencil_audit_contract_passed")
                is True
                and gates.get("amendment_18_19_evidence_contract_passed") is True
            ):
                raise RuntimeError("job A18/A19 evidence is not exact")
            step_hashes = [
                item["accepted_step_evidence_identity_sha256"]
                for item in expected["accepted_step_evidence_identities"]
            ]
            payload = {
                "job_reference_index": job_index,
                "job_id": job.get("job_id"),
                "scenario": job.get("scenario"),
                "cells": job.get("cells"),
                "dt_s": job.get("dt_s"),
                "accepted_step_record_count": expected["accepted_step_record_count"],
                "accepted_root_reference_count": expected["accepted_root_count"],
                "accepted_step_evidence_identity_sha256": step_hashes,
            }
            job_references.append(
                {
                    **payload,
                    "committed_job_reference_identity_sha256": (
                        _canonical_evidence_sha256(payload)
                    ),
                }
            )
            committed_root_reference_count += int(expected["accepted_root_count"])
        except (KeyError, TypeError, ValueError, RuntimeError) as exc:
            errors.append(f"campaign job {job.get('job_id', job_index)!r}: {exc}")
    reference_hashes = [
        item["committed_job_reference_identity_sha256"] for item in job_references
    ]
    if len(set(reference_hashes)) != len(reference_hashes):
        errors.append("committed job reference execution contexts are not one-to-one")
    all_jobs_bound = bool(jobs) and len(job_references) == len(jobs)
    if all_jobs_are_verified_wall_exclusions:
        all_jobs_bound = True
    elif not all_jobs_bound:
        errors.append("not every campaign job has exact committed A18/A19 evidence")
    diagnostic_passed = (
        expected_diagnostic is not None
        and expected_diagnostic["contract_passed"] is True
    )
    identity_payload = {
        "committed_job_reference_identity_sha256": reference_hashes,
        "diagnostic_partition_identity_sha256": (
            []
            if expected_diagnostic is None
            else expected_diagnostic["partition_identity_sha256"]
        ),
        "combined_diagnostic_partition_identity_sha256": (
            None
            if expected_diagnostic is None
            else expected_diagnostic[
                "combined_diagnostic_partition_identity_sha256"
            ]
        ),
    }
    return {
        "schema_version": 1,
        "numerical_amendment_ids": [
            INTERFACE_ROOT_AMENDMENT_ID,
            DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "committed_job_reference_count": len(job_references),
        "committed_job_reference_identity_count": len(reference_hashes),
        "committed_trajectory_root_reference_count": committed_root_reference_count,
        "committed_job_references": job_references,
        "committed_job_reference_identity_sha256": reference_hashes,
        "diagnostic_partition_count": (
            0 if expected_diagnostic is None else expected_diagnostic["partition_count"]
        ),
        "diagnostic_partition_identity_sha256": (
            []
            if expected_diagnostic is None
            else expected_diagnostic["partition_identity_sha256"]
        ),
        "diagnostic_execution_occurrence_count": (
            0
            if expected_diagnostic is None
            else expected_diagnostic["execution_occurrence_count"]
        ),
        "accepted_diagnostic_root_execution_count": (
            0
            if expected_diagnostic is None
            else expected_diagnostic["accepted_diagnostic_root_execution_count"]
        ),
        "diagnostic_partition_committed_trajectory_root_count": 0,
        "physical_root_identity_uniqueness_required_across_execution_and_reference_partitions": False,
        "execution_context_identity_uniqueness_required_within_each_partition": True,
        "raw_structural_evidence_identity_sha256": _canonical_evidence_sha256(
            identity_payload
        ),
        "all_committed_job_references_bound": all_jobs_bound,
        "all_jobs_are_verified_wall_exclusions": (
            all_jobs_are_verified_wall_exclusions
        ),
        "all_noncommitted_diagnostic_executions_bound": diagnostic_passed,
        "contract_passed": all_jobs_bound and diagnostic_passed and not errors,
        "physically_qualifying": False,
        "errors": errors,
    }


def _surface_flux_time_integrals(
    roots: Sequence[_AcceptedBalanceIntervalRoot],
    ledger: object,
    area_m2: float,
) -> dict[str, float]:
    """Integrate each accepted surface root over its own BE interval."""

    if not roots:
        raise ValueError("surface-flux integration needs an accepted root")
    if not math.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("surface-flux integration needs a positive finite area")
    return {
        "surface_conserved_water_flux_mol_m2_s": (ledger.boundary_water_out_mol / area_m2),
        "surface_conserved_hexane_flux_mol_m2_s": (ledger.boundary_hexane_out_mol / area_m2),
        "surface_independent_hexane_flux_mol_m2_s": math.fsum(
            root.interval_dt_s
            * (-root.assembly.dry_face_fluxes[-1].independent_water_flux_mol_m2_s)
            for root in roots
        ),
    }


def _step_record(
    step: AcceptedSameCellStep,
    attempts: Iterable[AttemptRecord],
    stage: str,
    recovery_mode: str,
    *,
    starts_after_boundary_discontinuity: bool = False,
) -> dict:
    face_macro = step.event_macrostep if isinstance(step, AcceptedFaceEventMacrostep) else None
    inventory = ci.inventory_snapshot(step.after.transport)
    water_phases = ci.water_phase_inventory_snapshot(step.after.transport)
    ledger = step.ledger
    surface_nt = step.assembly.candidate.dry_total_stefan_fluxes_mol_m2_s[-1]
    surface = step.assembly.dry_face_fluxes[-1]
    independent_hexane = -surface.independent_water_flux_mol_m2_s
    hexane_stefan_carrier = math.fsum(
        (
            surface.component.conserved_hexane_flux_mol_m2_s,
            -independent_hexane,
        )
    )
    accepted_leaves = step.substeps if isinstance(step, AdaptiveSameCellMacrostep) else (step,)
    if face_macro is not None:
        face_departure = face_macro.face_departure_step
        if face_departure is None:
            raise RuntimeError("accepted face event lost its departure remainder")
        wet_cut_treatment_applicability = (
            "exact_face_arrival_nonapplicable_strict_departure_remainder_applicable"
        )
        wet_cut_treatment_history = (face_departure.assembly.vanishing_wet_cut_treatment_audit,)
        wet_cut_treatment_route_history = ("strict_positive_departure_remainder",)
    else:
        wet_cut_treatment_applicability = "applicable_ordinary_same_cell_operator"
        wet_cut_treatment_history = tuple(
            leaf.assembly.vanishing_wet_cut_treatment_audit for leaf in accepted_leaves
        )
        wet_cut_treatment_route_history = tuple(
            "ordinary_same_cell_operator" for _leaf in accepted_leaves
        )
    accepted_balance_roots = _accepted_balance_interval_roots(step)
    scalar_post_root_polish = _scalar_post_root_polish_summary(
        tuple(
            (
                root.role,
                root.interval_dt_s,
                root.time_s,
                root.scalar_post_root_polish_audit,
            )
            for root in accepted_balance_roots
        )
    )
    accepted_node_temperatures = tuple(root.temperatures_k for root in accepted_balance_roots)
    maximum_particle_temperature = max(
        value for values in accepted_node_temperatures for value in values
    )
    minimum_particle_temperature = min(
        value for values in accepted_node_temperatures for value in values
    )
    maximum_node_radial_spread = max(
        max(values) - min(values) for values in accepted_node_temperatures
    )
    accepted_node_front_z = tuple(root.front_z for root in accepted_balance_roots)
    wet_retained_cap_transition_history = tuple(
        tuple(transition.value for transition in transitions)
        for leaf in accepted_leaves
        for transitions in (
            (leaf.ledger.wet_retained_cap_piece_transitions,)
            if isinstance(leaf, ci.CutIntegratorStep)
            else leaf.ledger.wet_retained_cap_piece_transition_history
        )
    )
    wet_retained_cap_active_master_cell_indices = tuple(
        step.after.transport.layout.wet_cell_indices[index]
        for index in ledger.wet_retained_cap_active_piece_indices
    )
    wet_retained_water_guard_loading = (
        ci.P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING
        if step.after.controls.enforce_p1ef_smooth_wet_retained_water_guard
        else None
    )
    wet_retained_water_guard_margin = (
        None
        if wet_retained_water_guard_loading is None
        else math.fsum(
            (
                wet_retained_water_guard_loading,
                -ledger.wet_retained_cap_maximum_loading,
            )
        )
    )
    wet_retained_water_capacity_loading = step.after.transport.config.wet.luikov.W_cap
    wet_retained_water_capacity_margin = math.fsum(
        (
            wet_retained_water_capacity_loading,
            -ledger.wet_retained_cap_maximum_loading,
        )
    )
    surface_flux_samples = tuple(_surface_flux_sample(root) for root in accepted_balance_roots)
    interface_saturation_endpoint_audit_summary = (
        _interface_composition_root_audit_summary(surface_flux_samples)
    )
    dry_energy_capacity_stencil_audit_summary = (
        _dry_energy_capacity_stencil_audit_summary(surface_flux_samples)
    )
    a20_evidence_presence = tuple(
        "moving_interface_actual_composition_force_path_audit" in sample
        for sample in surface_flux_samples
    )
    if all(a20_evidence_presence):
        moving_interface_actual_path_audit_summary = (
            _moving_interface_actual_path_audit_summary(surface_flux_samples)
        )
    elif any(a20_evidence_presence):
        raise RuntimeError("accepted step has only a partial A20 root partition")
    else:
        moving_interface_actual_path_audit_summary = {
            "schema_version": ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION,
            "engineering_amendment_id": ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID,
            "applicable": False,
            "accepted_root_count_expected": len(accepted_balance_roots),
            "accepted_root_count": 0,
            "accepted_root_identity_count": 0,
            "accepted_root_identities": [],
            "all_roots_independently_replayed": False,
            "all_roots_use_actual_linear_T_y_path_authority": False,
            "all_roots_serialize_full_conditional_continuous_phase_audit": False,
            "fictitious_thermodynamic_force_temperature_count": 0,
            "contract_passed": False,
            "physically_qualifying": False,
            "errors": [
                "moving-interface actual-(T,y) composition-force authority is not installed"
            ],
        }
    if not (
        interface_saturation_endpoint_audit_summary["accepted_root_count"]
        == len(accepted_balance_roots)
        == dry_energy_capacity_stencil_audit_summary["accepted_root_count"]
    ):
        raise RuntimeError("accepted-root A18/A19 evidence counts diverged")
    if all(a20_evidence_presence) and (
        moving_interface_actual_path_audit_summary["accepted_root_count"]
        != len(accepted_balance_roots)
    ):
        raise RuntimeError("accepted-root A20 evidence count diverged")
    endpoint_surface_sample = surface_flux_samples[-1]
    accepted_root_surface_oracle_passed = all(
        sample["surface_constitutive_oracle"]["passed"] is True for sample in surface_flux_samples
    )
    area = step.after.transport.geometry.master_grid.areas[-1]
    surface_flux_time_integrals = _surface_flux_time_integrals(
        accepted_balance_roots,
        ledger,
        area,
    )
    records = tuple(attempts)
    adaptive = isinstance(step, AdaptiveSameCellMacrostep)
    floor_contract = _adaptive_floor_contract(ledger.dt_s)
    effective_maximum_binary_depth = floor_contract["effective_maximum_binary_depth"]
    effective_terminal_leaf_dt_s = floor_contract["effective_terminal_leaf_dt_s"]
    balance_interval_durations = tuple(root.interval_dt_s for root in accepted_balance_roots)
    balance_interval_time_nodes = tuple(root.time_s for root in accepted_balance_roots)
    # Failed-leaf controller partitions and physical-event balance intervals are
    # different time structures.  In particular, an independently localized
    # event and its conservative remainder may be nonbinary and may be shorter
    # than the Amendment-09 failed-leaf floor.  Never feed those event widths
    # into the failed-leaf policy audit.
    failed_leaf_durations = (
        () if face_macro is not None else tuple(leaf.ledger.dt_s for leaf in accepted_leaves)
    )
    failed_leaf_time_nodes = (
        ()
        if face_macro is not None
        else tuple(leaf.after.transport.time_s for leaf in accepted_leaves)
    )
    input_identity = accepted_leaves[0].before is step.before
    leaf_chain_identity = all(
        right.before is left.after for left, right in zip(accepted_leaves, accepted_leaves[1:])
    )
    boundary_identity = all(leaf.boundary is step.boundary for leaf in accepted_leaves)
    attempts_rollback = all(record.rollback_identity_preserved for record in records)
    adaptive_attempts = tuple(
        record
        for record in records
        if record.adaptive_depth is not None and record.adaptive_depth > 0
    )
    adaptive_nonlinear_attempts = tuple(
        record
        for record in adaptive_attempts
        if not record.role.endswith("chart_validation")
        and not record.provisional_seed_attempt_budget_exhausted
    )
    adaptive_work_budget_exact = all(
        record.provisional_work_budget == DISCONTINUITY_LEAF_MAXIMUM_OPTIMIZER_FUNCTION_EVALUATIONS
        for record in adaptive_nonlinear_attempts
    )
    adaptive_attempt_paths = {
        (record.adaptive_depth, record.adaptive_path) for record in adaptive_attempts
    }
    adaptive_path_shape_contract = all(
        depth is not None
        and path is not None
        and len(path) == depth
        and set(path) <= {"L", "R"}
        and record.attempted_dt_s is not None
        and math.isclose(
            record.attempted_dt_s,
            math.ldexp(ledger.dt_s, -depth),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        and record.macro_branch_status in {"failed_leaf_bisected", "returned"}
        for record in adaptive_attempts
        for depth, path in ((record.adaptive_depth, record.adaptive_path),)
    )
    root_attempt_contract = not adaptive or (
        bool(records)
        and all(
            record.adaptive_depth == 0
            and record.adaptive_path == "root"
            and record.macro_branch_status == "initial_requested_failure"
            and record.accepted is False
            and record.attempted_dt_s is not None
            and math.isclose(
                record.attempted_dt_s,
                ledger.dt_s,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
            for record in records
            if record.adaptive_depth == 0
        )
        and any(record.adaptive_depth == 0 for record in records)
    )
    accepted_attempts = tuple(record for record in adaptive_attempts if record.accepted)
    accepted_leaf_attempt_contract = not adaptive or (
        len(accepted_attempts) == len(accepted_leaves)
        and all(record.macro_branch_status == "returned" for record in accepted_attempts)
        and all(
            math.isclose(
                record.attempted_dt_s or math.nan,
                leaf_dt,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
            for record, leaf_dt in zip(accepted_attempts, failed_leaf_durations)
        )
    )
    failed_paths = {
        record.adaptive_path
        for record in adaptive_attempts
        if record.macro_branch_status == "failed_leaf_bisected"
    }
    adaptive_binary_tree_contract = not adaptive or (
        {(1, "L"), (1, "R")} <= adaptive_attempt_paths
        and all(
            path is not None
            and len(path) < effective_maximum_binary_depth
            and (len(path) + 1, f"{path}L") in adaptive_attempt_paths
            and (len(path) + 1, f"{path}R") in adaptive_attempt_paths
            for path in failed_paths
        )
    )
    seed_attempt_budget_contract = True
    for depth, path in sorted(adaptive_attempt_paths):
        nonlinear_count = sum(
            record in adaptive_nonlinear_attempts
            for record in adaptive_attempts
            if record.adaptive_depth == depth and record.adaptive_path == path
        )
        maximum_seed_attempts = (
            DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
            if path == "L" * depth
            else DISCONTINUITY_NONREVERSAL_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
        )
        seed_attempt_budget_contract = (
            seed_attempt_budget_contract and 1 <= nonlinear_count <= maximum_seed_attempts
        )
    adaptive_depth_contract = not adaptive or (
        1 <= step.binary_depth <= effective_maximum_binary_depth
        and step.attempted_binary_depths == tuple(range(1, max(step.attempted_binary_depths) + 1))
        and max(step.attempted_binary_depths) == step.binary_depth
    )
    adaptive_leaf_contract = not adaptive or (
        1 < len(accepted_leaves) <= 2**effective_maximum_binary_depth
        and all(
            ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S <= leaf_dt <= ledger.dt_s
            for leaf_dt in failed_leaf_durations
        )
        and math.isclose(
            math.fsum(failed_leaf_durations),
            ledger.dt_s,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        and all(
            math.isclose(
                time_node,
                math.fsum(
                    (
                        step.before.transport.time_s,
                        *failed_leaf_durations[:index],
                    )
                ),
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
            for index, time_node in enumerate(failed_leaf_time_nodes, start=1)
        )
    )
    adaptive_macro_ledger_contract = (
        not adaptive
        or max(
            ledger.maximum_step_ledger_residual,
            ledger.maximum_cumulative_ledger_residual,
        )
        <= CONSERVATION_LIMIT
    )
    frozen_failed_leaf_controller_passed = all(
        (
            not adaptive or starts_after_boundary_discontinuity,
            adaptive_depth_contract,
            adaptive_leaf_contract,
            not adaptive
            or step.predictor_backtrack_maximum_binary_exponent
            == ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT,
            adaptive_work_budget_exact,
            adaptive_path_shape_contract,
            root_attempt_contract,
            accepted_leaf_attempt_contract,
            adaptive_binary_tree_contract,
            seed_attempt_budget_contract,
            adaptive_macro_ledger_contract,
            input_identity,
            leaf_chain_identity,
            boundary_identity,
            attempts_rollback,
        )
    )
    if face_macro is not None:
        failed_leaf_partitions: list[dict[str, object]] = []
    elif adaptive:
        failed_leaf_partitions = [
            {
                "role": "amendment09_failed_leaf_binary_partition",
                "duration_s": leaf.ledger.dt_s,
                "time_after_s": leaf.after.transport.time_s,
                "binary_depth": attempt.adaptive_depth,
                "binary_path": attempt.adaptive_path,
            }
            for leaf, attempt in zip(accepted_leaves, accepted_attempts, strict=True)
        ]
    else:
        failed_leaf_partitions = [
            {
                "role": "direct_requested_backward_euler",
                "duration_s": ledger.dt_s,
                "time_after_s": step.after.transport.time_s,
                "binary_depth": 0,
                "binary_path": None,
            }
        ]

    physical_event_partitions = (
        [
            {
                "role": root.role,
                "duration_s": root.interval_dt_s,
                "time_after_s": root.time_s,
            }
            for root in accepted_balance_roots
        ]
        if face_macro is not None
        else []
    )
    physical_event_contract: dict[str, object] | None = None
    if face_macro is not None:
        macro = face_macro
        event = macro.face_event_step
        departure = macro.face_departure_step
        evidence = macro.event_ordering_evidence
        if event is None or departure is None or evidence is None:
            raise RuntimeError("accepted face event lost its physical transaction evidence")
        threshold_authority = macro.before.transport.config.vanishing_wet_cut_threshold_authority
        threshold_authority_preserved = all(
            candidate is threshold_authority
            for candidate in (
                event.before.transport.config.vanishing_wet_cut_threshold_authority,
                event.after.config.vanishing_wet_cut_threshold_authority,
                departure.before.config.vanishing_wet_cut_threshold_authority,
                departure.after.transport.config.vanishing_wet_cut_threshold_authority,
                macro.after.transport.config.vanishing_wet_cut_threshold_authority,
                step.after.transport.config.vanishing_wet_cut_threshold_authority,
            )
        )
        physical_event_contract = {
            "event_kind": "interior_master_face_arrival_and_departure",
            "event_duration_s": macro.ledger.event_duration_s,
            "remainder_duration_s": macro.ledger.departure_remainder_s,
            "event_time_error_bound_s": macro.ledger.event_time_error_bound_s,
            "duration_closure_error_s": macro.ledger.duration_closure_error_s,
            "absolute_time_closure_error_s": (macro.ledger.absolute_time_closure_error_s),
            "accepted_balance_interval_count": (macro.ledger.accepted_balance_interval_count),
            "ordering_mode": macro.ledger.event_ordering_mode.value,
            "ordering": list(macro.ledger.event_ordering),
            "immutable_chain": {
                "macro_input_identity_preserved": macro.before is step.before,
                "event_input_identity_preserved": event.before is macro.before,
                "event_to_remainder_identity_preserved": (departure.before is event.after),
                "macro_output_identity_preserved": (
                    departure.after is macro.after and step.after is macro.after
                ),
                "boundary_identity_preserved": (
                    event.boundary is macro.boundary
                    and departure.boundary is macro.boundary
                    and step.boundary is macro.boundary
                ),
            },
            "vanishing_wet_cut_threshold_authority_transaction": {
                "threshold_authority": threshold_authority.value,
                "target_wet_volume_fraction": (threshold_authority.target_wet_volume_fraction),
                "ordinary_same_cell_operator_applicable": False,
                "exact_face_arrival_operator_applicable": False,
                "strict_positive_departure_remainder_operator_applicable": True,
                "departure_remainder_treatment_audit": asdict(
                    departure.assembly.vanishing_wet_cut_treatment_audit
                ),
                "amendment_16_applicability_passed": True,
                "preserved_across_arrival_and_departure": (threshold_authority_preserved),
            },
            "every_attempt_rollback_identity_preserved": attempts_rollback,
            "raw_macro_component_energy_ledger": {
                "water_step_residual_mol": ledger.water_step_residual_mol,
                "hexane_step_residual_mol": ledger.hexane_step_residual_mol,
                "energy_step_residual_j": ledger.energy_step_residual_j,
                "water_cumulative_residual_mol": (ledger.water_cumulative_residual_mol),
                "hexane_cumulative_residual_mol": (ledger.hexane_cumulative_residual_mol),
                "energy_cumulative_residual_j": ledger.energy_cumulative_residual_j,
                "normalized_water_step_residual": (ledger.normalized_water_step_residual),
                "normalized_hexane_step_residual": (ledger.normalized_hexane_step_residual),
                "normalized_energy_step_residual": (ledger.normalized_energy_step_residual),
                "normalized_water_cumulative_residual": (
                    ledger.normalized_water_cumulative_residual
                ),
                "normalized_hexane_cumulative_residual": (
                    ledger.normalized_hexane_cumulative_residual
                ),
                "normalized_energy_cumulative_residual": (
                    ledger.normalized_energy_cumulative_residual
                ),
                "maximum_normalized_step_residual": (ledger.maximum_step_ledger_residual),
                "maximum_normalized_cumulative_residual": (
                    ledger.maximum_cumulative_ledger_residual
                ),
                "limit": CONSERVATION_LIMIT,
            },
        }
    legacy_partition_alias = (
        "physical_event_partitions" if face_macro is not None else "failed_leaf_partitions"
    )
    return {
        "stage": stage,
        "material_volume_profile": material_profiles.serialize_integrator_state(step.after),
        "time_before_s": step.before.transport.time_s,
        "time_after_s": step.after.transport.time_s,
        "dt_s": ledger.dt_s,
        "recovery_mode": recovery_mode,
        "adaptive_subdivision_used": adaptive,
        "starts_after_imposed_boundary_discontinuity": (starts_after_boundary_discontinuity),
        "interior_face_event_used": face_macro is not None,
        "interior_face_event_index": (
            None
            if face_macro is None or face_macro.face_event_step is None
            else face_macro.face_event_step.after.arrival_face_index
        ),
        "interior_face_event_duration_s": (
            None if face_macro is None else face_macro.ledger.event_duration_s
        ),
        "interior_face_departure_remainder_s": (
            None if face_macro is None else face_macro.ledger.departure_remainder_s
        ),
        "interior_face_event_ordering_mode": (
            None
            if face_macro is None or face_macro.ledger.event_ordering_mode is None
            else face_macro.ledger.event_ordering_mode.value
        ),
        "interior_face_raw_macro_ledger_residual": (
            None if face_macro is None else face_macro.ledger.maximum_step_ledger_residual
        ),
        "adaptive_binary_depth": (step.binary_depth if adaptive else 0),
        "adaptive_attempted_binary_depths": (
            list(step.attempted_binary_depths) if adaptive else []
        ),
        "adaptive_predictor_backtrack_maximum_binary_exponent": (
            step.predictor_backtrack_maximum_binary_exponent if adaptive else None
        ),
        "adaptive_absolute_minimum_leaf_dt_s": (ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S),
        "adaptive_effective_maximum_binary_depth": (effective_maximum_binary_depth),
        "adaptive_effective_terminal_leaf_dt_s": effective_terminal_leaf_dt_s,
        "adaptive_minimum_accepted_leaf_dt_s": (
            min(failed_leaf_durations) if failed_leaf_durations else None
        ),
        "adaptive_no_leaf_below_absolute_floor": all(
            leaf_dt >= ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S
            for leaf_dt in failed_leaf_durations
        ),
        "adaptive_floor_was_binding": floor_contract["floor_was_binding"],
        "adaptive_depth_cap_was_binding": floor_contract["depth_cap_was_binding"],
        "initial_requested_failure_jacobian_audit": (
            step.initial_failure_audit if adaptive else None
        ),
        "frozen_failed_leaf_controller_contract": {
            "passed": frozen_failed_leaf_controller_passed,
            "policy_id": ADAPTIVE_FAILED_LEAF_POLICY_ID,
            "licensed_only_after_imposed_boundary_discontinuity": (
                not adaptive or starts_after_boundary_discontinuity
            ),
            "absolute_minimum_leaf_dt_s": (ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S),
            "maximum_binary_depth_cap": ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH,
            "effective_maximum_binary_depth": effective_maximum_binary_depth,
            "effective_terminal_leaf_dt_s": effective_terminal_leaf_dt_s,
            "accepted_binary_depth_bounded": adaptive_depth_contract,
            "maximum_leaf_count": 2**effective_maximum_binary_depth,
            "leaf_grid_disclosed_and_bounded": adaptive_leaf_contract,
            "minimum_accepted_leaf_dt_s": (
                min(failed_leaf_durations) if failed_leaf_durations else None
            ),
            "no_leaf_below_absolute_floor": all(
                leaf_dt >= ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S
                for leaf_dt in failed_leaf_durations
            ),
            "floor_was_binding": floor_contract["floor_was_binding"],
            "depth_cap_was_binding": floor_contract["depth_cap_was_binding"],
            "per_seed_and_per_leaf_work_budget_unchanged_from_amendment03": True,
            "total_tree_work_envelope_changed_by_amendment09": True,
            "predictor_backtrack_maximum_binary_exponent": (
                ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
            ),
            "provisional_optimizer_function_evaluations_per_seed": (
                DISCONTINUITY_LEAF_MAXIMUM_OPTIMIZER_FUNCTION_EVALUATIONS
            ),
            "reversal_maximum_nonlinear_seed_attempts": (
                DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
            ),
            "nonreversal_maximum_nonlinear_seed_attempts": (
                DISCONTINUITY_NONREVERSAL_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
            ),
            "fixed_work_budget_observed": adaptive_work_budget_exact,
            "attempt_paths_form_failed_leaf_binary_tree": (
                adaptive_path_shape_contract and adaptive_binary_tree_contract
            ),
            "initial_requested_macro_failure_disclosed": root_attempt_contract,
            "one_returned_attempt_per_accepted_leaf": (accepted_leaf_attempt_contract),
            "per_leaf_seed_attempt_budget_observed": seed_attempt_budget_contract,
            "atomic_macro_input_identity_preserved": input_identity,
            "immutable_leaf_chain_identity_preserved": leaf_chain_identity,
            "exact_post_jump_boundary_identity_preserved": boundary_identity,
            "every_attempt_rollback_identity_preserved": attempts_rollback,
            "raw_macro_ledger_recomputed_from_original_and_final_inventory": (adaptive),
            "raw_macro_ledger_le_1e_10": adaptive_macro_ledger_contract,
            "arbitrary_subdivision_allowed": False,
        },
        "failed_leaf_partitions": failed_leaf_partitions,
        "physical_event_partitions": physical_event_partitions,
        "physical_event_contract": physical_event_contract,
        "legacy_accepted_internal_substep_reporting": {
            "authoritative": False,
            "alias_of": legacy_partition_alias,
            "semantic_role": "backward_compatibility_balance_interval_alias_only",
            "must_not_be_used_to_apply_amendment09_floor": True,
        },
        "accepted_internal_substep_count": len(accepted_balance_roots),
        "accepted_internal_substep_dt_s": list(balance_interval_durations),
        "accepted_internal_time_nodes_s": list(balance_interval_time_nodes),
        "scalar_post_root_polish": scalar_post_root_polish,
        "interface_saturation_endpoint_audit_summary": (
            interface_saturation_endpoint_audit_summary
        ),
        "dry_energy_capacity_stencil_audit_summary": (
            dry_energy_capacity_stencil_audit_summary
        ),
        "moving_interface_actual_composition_force_path_audit_summary": (
            moving_interface_actual_path_audit_summary
        ),
        "amendment_18_19_evidence_contract_passed": True,
        "amendment_20_evidence_contract_passed": (
            moving_interface_actual_path_audit_summary["contract_passed"]
        ),
        "attempts": [asdict(record) for record in records],
        "wall_time_s": math.fsum(record.wall_time_s for record in records),
        "nonlinear_evaluations": ledger.nonlinear_evaluations,
        "condition_proxy": ledger.condition_proxy,
        "maximum_scaled_residual": ledger.maximum_scaled_residual,
        "maximum_step_ledger_residual": ledger.maximum_step_ledger_residual,
        "maximum_cumulative_ledger_residual": (ledger.maximum_cumulative_ledger_residual),
        "normalized_water_step_residual": (ledger.normalized_water_step_residual),
        "normalized_hexane_step_residual": (ledger.normalized_hexane_step_residual),
        "normalized_energy_step_residual": (ledger.normalized_energy_step_residual),
        "normalized_water_cumulative_residual": (ledger.normalized_water_cumulative_residual),
        "normalized_hexane_cumulative_residual": (ledger.normalized_hexane_cumulative_residual),
        "normalized_energy_cumulative_residual": (ledger.normalized_energy_cumulative_residual),
        "minimum_fractional_distance_to_bound": (ledger.minimum_fractional_distance_to_bound),
        "vanishing_wet_cut_treatment_applicability": (wet_cut_treatment_applicability),
        "vanishing_wet_cut_treatment": (
            None if not wet_cut_treatment_history else asdict(wet_cut_treatment_history[-1])
        ),
        "vanishing_wet_cut_treatment_history": [
            asdict(audit) for audit in wet_cut_treatment_history
        ],
        "vanishing_wet_cut_treatment_route_history": list(wet_cut_treatment_route_history),
        "vanishing_wet_cut_activation_count": sum(
            audit.activation_count for audit in wet_cut_treatment_history
        ),
        "maximum_vanishing_wet_cut_water_pair_sum_difference_mol_s": max(
            (abs(audit.water_pair_sum_difference_mol_s) for audit in wet_cut_treatment_history),
            default=0.0,
        ),
        "maximum_vanishing_wet_cut_common_energy_pair_sum_difference_w": max(
            (abs(audit.common_energy_pair_sum_difference_w) for audit in wet_cut_treatment_history),
            default=0.0,
        ),
        "wet_retained_cap_active_piece_indices": list(ledger.wet_retained_cap_active_piece_indices),
        "wet_retained_cap_active_master_cell_indices": list(
            wet_retained_cap_active_master_cell_indices
        ),
        "wet_retained_cap_active_piece_count": (ledger.wet_retained_cap_active_piece_count),
        "wet_retained_cap_maximum_loading": (ledger.wet_retained_cap_maximum_loading),
        "wet_retained_water_guard_loading": wet_retained_water_guard_loading,
        "wet_retained_water_guard_margin": wet_retained_water_guard_margin,
        "wet_retained_water_capacity_loading": (wet_retained_water_capacity_loading),
        "wet_retained_water_capacity_margin": (wet_retained_water_capacity_margin),
        "wet_retained_cap_maximum_capacity_dual_over_rt": (
            ledger.wet_retained_cap_maximum_capacity_dual_over_rt
        ),
        "wet_retained_cap_maximum_complementarity_product": (
            ledger.wet_retained_cap_maximum_complementarity_product
        ),
        "wet_retained_cap_complementarity_exact": (
            ledger.wet_retained_cap_maximum_complementarity_product == 0.0
        ),
        "wet_retained_cap_exact_graph_without_tolerance": (
            ledger.wet_retained_cap_exact_graph_without_tolerance
        ),
        "wet_retained_cap_piece_transition_history": list(wet_retained_cap_transition_history),
        "wet_retained_cap_entry_piece_count": (ledger.wet_retained_cap_entry_piece_count),
        "wet_retained_cap_continuation_piece_count": (
            ledger.wet_retained_cap_continuation_piece_count
        ),
        "wet_retained_cap_exit_piece_count": (ledger.wet_retained_cap_exit_piece_count),
        "solver_controls_match_campaign_contract": all(
            leaf.before.controls == CONTROLS and leaf.after.controls == CONTROLS
            for leaf in accepted_leaves
        ),
        "front_z": step.after.transport.geometry.front.z,
        "minimum_front_z": min(accepted_node_front_z),
        "maximum_front_z": max(accepted_node_front_z),
        "particle_radius_m": step.after.transport.geometry.master_grid.R,
        "interface_temperature_k": step.after.last_interface_temperature_k,
        "maximum_particle_temperature_k": maximum_particle_temperature,
        "minimum_particle_temperature_k": minimum_particle_temperature,
        "instantaneous_radial_temperature_spread_k": (maximum_node_radial_spread),
        "minimum_dry_temperature_k": min(step.after.transport.dry_temperatures_k),
        "surface_dry_temperature_k": step.after.transport.dry_temperatures_k[-1],
        "boundary_temperature_k": step.boundary.temperature_k,
        "retained_water_source_temperature_extrapolation": (
            max(maximum_particle_temperature, step.boundary.temperature_k)
            > RETAINED_WATER_DIRECT_TEMPERATURE_MAX_K
        ),
        "surface_total_stefan_flux_mol_m2_s": surface_nt,
        "surface_conserved_water_flux_mol_m2_s": (surface.component.conserved_water_flux_mol_m2_s),
        "surface_conserved_hexane_flux_mol_m2_s": (
            surface.component.conserved_hexane_flux_mol_m2_s
        ),
        "surface_hexane_stefan_carrier_flux_mol_m2_s": (hexane_stefan_carrier),
        "surface_independent_hexane_flux_mol_m2_s": independent_hexane,
        "surface_independent_water_flux_mol_m2_s": (surface.independent_water_flux_mol_m2_s),
        "surface_retained_water_flux_mol_m2_s": (surface.retained_water_flux_mol_m2_s),
        "surface_gas_entropy_production_w_m3_k": (surface.gas_entropy_production_w_m3_k),
        "surface_binary_gas_entropy_production_w_m3_k": (surface.gas_entropy_production_w_m3_k),
        "surface_total_entropy_production_w_m3_k": (surface.total_entropy_production_w_m3_k),
        "surface_constitutive_primitives": endpoint_surface_sample[
            "surface_constitutive_primitives"
        ],
        "surface_constitutive_oracle": endpoint_surface_sample["surface_constitutive_oracle"],
        "accepted_root_surface_constitutive_oracle_passed": (accepted_root_surface_oracle_passed),
        **(
            {"surface_film_audit": asdict(surface.surface_film_audit)}
            if surface.surface_film_audit is not None
            else {}
        ),
        "accepted_node_surface_flux_samples": list(surface_flux_samples),
        "surface_flux_time_integrals_mol_s_m2": surface_flux_time_integrals,
        "raw_macro_water_residual_mol": ledger.water_step_residual_mol,
        "raw_macro_hexane_residual_mol": ledger.hexane_step_residual_mol,
        "raw_macro_energy_residual_j": ledger.energy_step_residual_j,
        "water_inventory_mol": inventory.total_water_mol,
        "retained_water_inventory_mol": water_phases.retained_water_mol,
        "pore_vapor_water_inventory_mol": water_phases.pore_vapor_water_mol,
        "free_liquid_water_inventory_mol": water_phases.free_liquid_water_mol,
        "hexane_inventory_mol": inventory.total_hexane_mol,
        "energy_inventory_j": inventory.total_energy_j,
        "cumulative_boundary_water_out_mol": (
            step.after.corrected_cumulative_boundary_water_out_mol
        ),
        "cumulative_boundary_hexane_out_mol": (
            step.after.corrected_cumulative_boundary_hexane_out_mol
        ),
        "cumulative_boundary_energy_out_j": (step.after.corrected_cumulative_boundary_energy_out_j),
    }


def _radius_seed_stages(
    start_radius_m: float,
    target_radius_m: float,
) -> tuple[float, ...]:
    """Return deterministic seed-only stages ending exactly at ``target``."""

    start = _validate_radius(start_radius_m)
    target = _validate_radius(target_radius_m)
    if start == target:
        return ()
    stage_count = math.ceil(abs(target - start) / MAXIMUM_RADIUS_SEED_CONTINUATION_INCREMENT_M)
    return tuple(
        target if index == stage_count else start + (target - start) * index / stage_count
        for index in range(1, stage_count + 1)
    )


def _solve_n2_first_step_at_radius(
    radius_m: float,
    *,
    cross_radius_seed_candidate: cut.CutTransportUnknowns | None = None,
    cross_radius_seed_source_m: float | None = None,
) -> tuple[ci.CutIntegratorStep, FixedRadiusFirstStepAudit]:
    """Solve the unchanged N=2 first step with seed-only radius recovery."""

    radius_m = _validate_radius(radius_m)
    direct_before = _initial_integrator(2, radius_m)
    if cross_radius_seed_candidate is not None:
        if cross_radius_seed_source_m is None:
            raise ValueError("cross-radius first-step seed needs its source radius")
        source_radius_m = _validate_radius(cross_radius_seed_source_m)
        cross_seed = _seed(
            cross_radius_seed_candidate,
            (
                f"certified {source_radius_m:g} m first-step root reused only "
                f"as nonlinear seed for independent {radius_m:g} m equations"
            ),
        )
        try:
            ci._validate_seed(direct_before, cross_seed)  # noqa: SLF001
        except ValueError as exc:
            record = AttemptRecord(
                role="cross_radius_first_step_seed_chart_validation",
                predictor_fraction=1.0,
                wall_time_s=0.0,
                accepted=False,
                nonlinear_evaluations=0,
                maximum_scaled_residual=math.inf,
                condition_proxy=math.inf,
                rollback_identity_preserved=True,
                error=str(exc),
            )
            raise CrossRadiusSeedTrajectoryError(
                "cross-radius first-step seed is outside the target chart",
                target_radius_m=radius_m,
                source_radius_m=source_radius_m,
                attempts=(record,),
                requested_input_state=direct_before,
                rollback_state=direct_before,
            ) from exc
        cross_step, cross_record, cross_step_error = _try_step_detailed(
            direct_before,
            SEGMENT_DURATION_S,
            HOT_LEAN_BOUNDARY,
            cross_seed,
            role="fixed_radius_first_step_from_cross_radius_seed_only_candidate",
            predictor_fraction=1.0,
        )
        if cross_step is None:
            if not cross_record.rollback_identity_preserved:
                raise RuntimeError("cross-radius first-step rejection did not roll back")
            raise CrossRadiusSeedTrajectoryError(
                (
                    "single capped cross-radius first-step seed failed without "
                    f"state commit: {cross_record.error}"
                ),
                target_radius_m=radius_m,
                source_radius_m=source_radius_m,
                attempts=(cross_record,),
                requested_input_state=direct_before,
                rollback_state=cross_step_error.rollback_state,
                last_step_error=cross_step_error,
            )
        cross_occurrence = _amendment_18_19_execution_occurrence(
            cross_step,
            execution_classification="diagnostic_bootstrap_or_radius_seed_execution",
            execution_context={
                "execution_partition_id": "bootstrap_radius_seed_executions",
                "execution_role": (
                    "fixed_radius_first_step_from_cross_radius_seed_only_candidate"
                ),
                "radius_seed_execution_index": 0,
                "mesh_cells": 2,
                "target_radius_m": radius_m,
                "cross_radius_seed_source_m": source_radius_m,
            },
        )
        return cross_step, FixedRadiusFirstStepAudit(
            radius_m=radius_m,
            direct_representative_seed_accepted=False,
            direct_rejection_rollback_identity_preserved=None,
            seed_route="cross_radius_seed_only_exact_root",
            seed_only_radius_stages_m=(),
            cross_radius_seed_source_m=source_radius_m,
            cross_radius_attempts=(cross_record,),
            amendment_18_19_diagnostic_execution_occurrences=(cross_occurrence,),
        )
    direct_seed = _seed(
        N2_ROOT,
        "representative-radius exact root reused only as a nonlinear seed",
    )
    direct, direct_record = _try_step(
        direct_before,
        SEGMENT_DURATION_S,
        HOT_LEAN_BOUNDARY,
        direct_seed,
        role="fixed_radius_direct_representative_root_seed",
        predictor_fraction=1.0,
    )
    if direct is not None:
        direct_occurrence = _amendment_18_19_execution_occurrence(
            direct,
            execution_classification="diagnostic_bootstrap_or_radius_seed_execution",
            execution_context={
                "execution_partition_id": "bootstrap_radius_seed_executions",
                "execution_role": "fixed_radius_direct_representative_root_seed",
                "radius_seed_execution_index": 0,
                "mesh_cells": 2,
                "target_radius_m": radius_m,
            },
        )
        return direct, FixedRadiusFirstStepAudit(
            radius_m=radius_m,
            direct_representative_seed_accepted=True,
            direct_rejection_rollback_identity_preserved=None,
            seed_route="direct_representative_root_seed",
            seed_only_radius_stages_m=(),
            amendment_18_19_diagnostic_execution_occurrences=(direct_occurrence,),
        )
    if not direct_record.rollback_identity_preserved:
        raise RuntimeError("rejected fixed-radius first step did not roll back exactly")

    # The 1.25 mm equation root is directly seedable by the representative
    # candidate and shortens the known upper-envelope recovery.  No accepted
    # state, conserved inventory, or time history is carried across radii.
    anchor_radius = (
        UPPER_RADIUS_SEED_ANCHOR_M if radius_m > UPPER_RADIUS_SEED_ANCHOR_M else RADIUS_M
    )
    anchor_before = _initial_integrator(2, anchor_radius)
    anchor, anchor_record = _try_step(
        anchor_before,
        SEGMENT_DURATION_S,
        HOT_LEAN_BOUNDARY,
        direct_seed,
        role="fixed_radius_seed_only_anchor",
        predictor_fraction=1.0,
    )
    if anchor is None:
        raise RuntimeError(
            f"fixed-radius seed anchor failed without state commit: {anchor_record.error}"
        )
    diagnostic_occurrences = [
        _amendment_18_19_execution_occurrence(
            anchor,
            execution_classification="diagnostic_bootstrap_or_radius_seed_execution",
            execution_context={
                "execution_partition_id": "bootstrap_radius_seed_executions",
                "execution_role": "fixed_radius_seed_only_anchor",
                "radius_seed_execution_index": 0,
                "mesh_cells": 2,
                "target_radius_m": anchor_radius,
                "final_requested_radius_m": radius_m,
            },
        )
    ]
    stages = _radius_seed_stages(anchor_radius, radius_m)
    previous = anchor
    for stage_index, stage_radius in enumerate(stages, start=1):
        stage_before = _initial_integrator(2, stage_radius)
        stage_seed = _seed(
            previous.assembly.candidate,
            (
                f"accepted {previous.after.transport.geometry.master_grid.R:g} m "
                "root reused only as next-radius nonlinear seed"
            ),
        )
        current, record = _try_step(
            stage_before,
            SEGMENT_DURATION_S,
            HOT_LEAN_BOUNDARY,
            stage_seed,
            role="fixed_radius_seed_only_continuation",
            predictor_fraction=1.0,
        )
        if current is None:
            raise RuntimeError(
                f"seed-only radius continuation failed without state commit: {record.error}"
            )
        diagnostic_occurrences.append(
            _amendment_18_19_execution_occurrence(
                current,
                execution_classification=(
                    "diagnostic_bootstrap_or_radius_seed_execution"
                ),
                execution_context={
                    "execution_partition_id": (
                        "bootstrap_radius_seed_executions"
                    ),
                    "execution_role": "fixed_radius_seed_only_continuation",
                    "radius_seed_execution_index": stage_index,
                    "mesh_cells": 2,
                    "target_radius_m": stage_radius,
                    "final_requested_radius_m": radius_m,
                },
            )
        )
        previous = current
    if previous.after.transport.geometry.master_grid.R != radius_m:
        raise RuntimeError("radius continuation did not reach the requested geometry")
    return previous, FixedRadiusFirstStepAudit(
        radius_m=radius_m,
        direct_representative_seed_accepted=False,
        direct_rejection_rollback_identity_preserved=True,
        seed_route="seed_only_radius_continuation",
        seed_only_radius_stages_m=stages,
        amendment_18_19_diagnostic_execution_occurrences=tuple(
            diagnostic_occurrences
        ),
    )


def _bootstrap_first_steps(
    maximum_mesh: int,
    radius_m: float | None = None,
    *,
    cross_radius_n2_seed_candidate: cut.CutTransportUnknowns | None = None,
    cross_radius_seed_source_m: float | None = None,
) -> tuple[
    dict[int, ci.CutIntegratorStep],
    dict[int, dict],
    FixedRadiusFirstStepAudit,
]:
    if maximum_mesh not in BOOTSTRAP_MESHES:
        raise ValueError("maximum mesh must be one of the declared campaign meshes")
    if radius_m is None:
        radius_m = RADIUS_M
    radius_m = _validate_radius(radius_m)
    started = perf_counter()
    step, first_step_audit = _solve_n2_first_step_at_radius(
        radius_m,
        cross_radius_seed_candidate=cross_radius_n2_seed_candidate,
        cross_radius_seed_source_m=cross_radius_seed_source_m,
    )
    source = step.before
    diagnostics = {
        2: {
            "wall_time_s": perf_counter() - started,
            "nonlinear_evaluations": step.ledger.nonlinear_evaluations,
            "maximum_scaled_residual": step.ledger.maximum_scaled_residual,
            "condition_proxy": step.ledger.condition_proxy,
            "minimum_dry_temperature_k": min(step.after.transport.dry_temperatures_k),
            "minimum_fractional_distance_to_bound": (
                step.ledger.minimum_fractional_distance_to_bound
            ),
            "scalar_post_root_polish": (
                _cut_step_scalar_post_root_polish_summary(
                    step,
                    "mesh_bootstrap_exact_root",
                )
            ),
            "radius_m": radius_m,
            "fixed_radius_first_step_audit": asdict(first_step_audit),
        }
    }
    requested: dict[int, ci.CutIntegratorStep] = {2: step} if maximum_mesh == 2 else {}
    candidate = step.assembly.candidate
    for bootstrap_index, cells in enumerate(BOOTSTRAP_MESHES[1:], start=1):
        if cells > maximum_mesh:
            break
        target = _initial_integrator(cells, radius_m)
        continuation = cc.prolongate_same_front_candidate(
            source,
            candidate,
            target,
            stefan_flux_scale_floor_mol_m2_s=STEFAN_SCALE_FLOOR,
            label=f"N={source.transport.geometry.master_grid.n} to N={cells}",
            master_face_event_chart=cc.MasterFaceEventChart(
                requested_dt_s=SEGMENT_DURATION_S,
                label=(
                    "mesh bootstrap master-face event route; the exact event "
                    "chart the mesh-continuation refusal names"
                ),
            ),
        )
        started = perf_counter()
        if isinstance(continuation, cc.MasterFaceEventContinuation):
            macro = event_orchestrator.advance_one_interior_face_macrostep(
                target,
                SEGMENT_DURATION_S,
                HOT_LEAN_BOUNDARY,
                continuation.same_cell_seed,
                face_event_controls=continuation.face_event_controls,
                face_event_seed=continuation.face_event_seed,
            )
            if macro.ledger.branch is not (
                event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE
            ):
                raise RuntimeError(
                    "mesh bootstrap master-face route did not commit one exact "
                    "arrival plus its strict departure"
                )
            step = _adapt_accepted_face_event_macrostep(continuation, macro)
        else:
            step = ci.advance_same_cell_backward_euler(
                target,
                SEGMENT_DURATION_S,
                HOT_LEAN_BOUNDARY,
                continuation.seed,
            )
        amendment_occurrence = _amendment_18_19_execution_occurrence(
            step,
            execution_classification="diagnostic_bootstrap_or_radius_seed_execution",
            execution_context={
                "execution_partition_id": "bootstrap_radius_seed_executions",
                "execution_role": "mesh_bootstrap_continuation",
                "mesh_bootstrap_execution_index": bootstrap_index,
                "source_mesh_cells": source.transport.geometry.master_grid.n,
                "target_mesh_cells": cells,
                "radius_m": radius_m,
            },
        )
        diagnostics[cells] = {
            "wall_time_s": perf_counter() - started,
            "nonlinear_evaluations": step.ledger.nonlinear_evaluations,
            "maximum_scaled_residual": step.ledger.maximum_scaled_residual,
            "maximum_cumulative_ledger_residual": (step.ledger.maximum_cumulative_ledger_residual),
            "condition_proxy": step.ledger.condition_proxy,
            "minimum_dry_temperature_k": min(step.after.transport.dry_temperatures_k),
            "minimum_fractional_distance_to_bound": (
                step.ledger.minimum_fractional_distance_to_bound
            ),
            "scalar_post_root_polish": (
                _cut_step_scalar_post_root_polish_summary(
                    step,
                    "mesh_bootstrap_exact_root",
                )
                if isinstance(step, ci.CutIntegratorStep)
                else _scalar_post_root_polish_summary(
                    tuple(
                        (
                            role,
                            interval_dt_s,
                            step.after.transport.time_s,
                            None,
                        )
                        for role, interval_dt_s in (
                            (
                                "mesh_bootstrap_exact_face_arrival",
                                step.event_macrostep.ledger.event_duration_s,
                            ),
                            (
                                "mesh_bootstrap_strict_face_departure",
                                step.event_macrostep.ledger.departure_remainder_s,
                            ),
                        )
                    )
                )
            ),
            "master_face_event_route": (
                None
                if isinstance(step, ci.CutIntegratorStep)
                else asdict(continuation.record)
            ),
            "amendment_18_19_diagnostic_execution_occurrence": (
                amendment_occurrence
            ),
        }
        if cells in DEFAULT_MESHES or cells == maximum_mesh:
            requested[cells] = step
        source = target
        candidate = step.assembly.candidate
    radius_occurrences = list(
        first_step_audit.amendment_18_19_diagnostic_execution_occurrences
    )
    mesh_occurrences = [
        diagnostics[cells]["amendment_18_19_diagnostic_execution_occurrence"]
        for cells in sorted(diagnostics)
        if cells != 2
    ]
    diagnostics[2]["amendment_18_19_diagnostic_execution_partition"] = (
        _aggregate_amendment_18_19_diagnostic_partition(
            "bootstrap_radius_seed_executions",
            (*radius_occurrences, *mesh_occurrences),
        )
    )
    return requested, diagnostics, first_step_audit


def _signed_extrema(values: Sequence[float]) -> dict:
    """Summarize one signed flux sequence without imposing a dead band."""

    finite = bool(values) and all(math.isfinite(value) for value in values)
    if not finite:
        return {
            "sample_count": len(values),
            "all_finite": False,
            "minimum_mol_m2_s": None,
            "maximum_mol_m2_s": None,
            "first_mol_m2_s": None,
            "last_mol_m2_s": None,
            "contains_strictly_outward": False,
            "contains_strictly_inward": False,
        }
    return {
        "sample_count": len(values),
        "all_finite": True,
        "minimum_mol_m2_s": min(values),
        "maximum_mol_m2_s": max(values),
        "first_mol_m2_s": values[0],
        "last_mol_m2_s": values[-1],
        "contains_strictly_outward": max(values) > 0.0,
        "contains_strictly_inward": min(values) < 0.0,
    }


def _finite_pulse_component_response(
    baseline_before_fluxes: Sequence[float],
    pulse_fluxes: Sequence[float],
    baseline_return_fluxes: Sequence[float],
) -> dict:
    """Audit outward/inward/outward response of one conserved component.

    This function is appropriate for a component whose expected response is
    outward under the declared baseline and inward under the declared pulse.
    It must not be called with algebraic total Stefan flux as a surrogate.
    """

    sequences = (
        baseline_before_fluxes,
        pulse_fluxes,
        baseline_return_fluxes,
    )
    finite = all(sequence for sequence in sequences) and all(
        math.isfinite(value) for sequence in sequences for value in sequence
    )
    outward_before = finite and max(baseline_before_fluxes) > 0.0
    inward_during = finite and min(pulse_fluxes) < 0.0
    outward_after = finite and max(baseline_return_fluxes) > 0.0
    return {
        "baseline_before": _signed_extrema(baseline_before_fluxes),
        "finite_pulse": _signed_extrema(pulse_fluxes),
        "baseline_return": _signed_extrema(baseline_return_fluxes),
        "strict_sign_only_no_dead_band": True,
        "outward_to_inward_observed": outward_before and inward_during,
        "inward_to_outward_return_observed": inward_during and outward_after,
        "outward_inward_outward_observed": (outward_before and inward_during and outward_after),
    }


def _sign_transition_bracket(
    step_records: Sequence[dict],
    before_stage: str,
    after_stage: str,
    field: str,
    *,
    before_sign: int,
) -> dict:
    """Bracket a sampled sign transition without interpolating a zero time."""

    before = [record for record in step_records if record["stage"] == before_stage]
    after = [record for record in step_records if record["stage"] == after_stage]
    if not before or not after or before_sign not in (-1, 1):
        return {"observed": False, "lower_s": None, "upper_s": None, "width_s": None}
    previous = before[-1]
    previous_value = previous[field]
    expected_before = before_sign * previous_value > 0.0
    if not expected_before or not math.isfinite(previous_value):
        return {"observed": False, "lower_s": None, "upper_s": None, "width_s": None}
    for current in after:
        value = current[field]
        if not math.isfinite(value):
            break
        if before_sign * value < 0.0:
            lower = current["time_before_s"]
            upper = current["time_after_s"]
            return {
                "observed": True,
                "lower_s": lower,
                "upper_s": upper,
                "midpoint_s": 0.5 * (lower + upper),
                "width_s": upper - lower,
                "interpolated_zero_time": False,
            }
        previous_value = value
    return {"observed": False, "lower_s": None, "upper_s": None, "width_s": None}


def _stage_flux_time_integral(
    step_records: Sequence[dict],
    stage: str,
    field: str,
) -> float:
    """Return the accepted BE right-endpoint flux integral in mol s/m2."""

    return math.fsum(
        record.get("surface_flux_time_integrals_mol_s_m2", {}).get(
            field,
            record["dt_s"] * record[field],
        )
        for record in step_records
        if record["stage"] == stage
    )


def _pulse_layer_resolution_diagnostic(
    cells: int,
    radius_m: float | None = None,
) -> dict:
    """Disclose the optimistic unretarded pore-diffusion length scale."""

    if radius_m is None:
        radius_m = RADIUS_M
    radius_m = _validate_radius(radius_m)
    diffusivity = DRY_CONFIG.binary_diffusivity.value_m2_s
    length = math.sqrt(diffusivity * SEGMENT_DURATION_S)
    spacing = radius_m / cells
    return {
        "fixed_particle_radius_m": radius_m,
        "configured_effective_binary_diffusivity_m2_s": diffusivity,
        "pulse_duration_s": SEGMENT_DURATION_S,
        "unretarded_pore_diffusion_length_m": length,
        "uniform_master_grid_spacing_m": spacing,
        "master_cells_per_unretarded_diffusion_length": length / spacing,
        "includes_local_equilibrium_storage_retardation": False,
        "is_a_resolution_certificate": False,
        "interpretation": (
            "sqrt(D_eff*tau)/dr is an optimistic pore scale; local-equilibrium "
            "storage can shorten the apparent layer, so empirical Cauchy "
            "contraction or a finer/boundary-refined reference remains required"
        ),
        "ideal_Dirichlet_jump_surface_flux_samples_are_diagnostic": True,
        "instantaneous_jump_flux_is_mesh_acceptance_quantity": False,
        "finite_window_integrated_fluxes_and_inventories_are_convergence_quantities": (True),
        "finite_external_film_installed": False,
    }


def _accepted_grid_contract(
    step_records: Sequence[dict],
    dt_s: float,
    scenario: StressScenario | None = None,
) -> dict:
    """Derive the exact requested-grid claim from accepted trajectory records."""

    if scenario is None:
        scenario = COMBINED_SCENARIO
    requested_steps = round(SEGMENT_DURATION_S / dt_s)
    expected_count = TRAJECTORY_SEGMENT_COUNT * requested_steps
    expected_stages = tuple(stage for stage, _ in scenario.segments for _ in range(requested_steps))
    actual_stages = tuple(record["stage"] for record in step_records)
    count_matches = len(step_records) == expected_count
    durations_match = count_matches and all(record["dt_s"] == dt_s for record in step_records)
    times_match = count_matches and all(
        math.isclose(
            record["time_before_s"],
            index * dt_s,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        and math.isclose(
            record["time_after_s"],
            (index + 1) * dt_s,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        for index, record in enumerate(step_records)
    )
    stages_match = actual_stages == expected_stages
    transition_indices = (requested_steps, 2 * requested_steps)
    transitions_aligned = count_matches and all(
        math.isclose(
            step_records[index]["time_before_s"],
            segment * SEGMENT_DURATION_S,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        for segment, index in enumerate(transition_indices, start=1)
    )
    passed = all((count_matches, durations_match, times_match, stages_match, transitions_aligned))
    adaptive_subdivision_used = any(
        record.get("adaptive_subdivision_used", False) for record in step_records
    )
    frozen_failed_leaf_controller_contract_passed = all(
        isinstance(record.get("frozen_failed_leaf_controller_contract"), dict)
        and record["frozen_failed_leaf_controller_contract"].get("passed") is True
        for record in step_records
    )
    internal_nodes_respect_macro_intervals = all(
        (
            not record.get("accepted_internal_time_nodes_s")
            or (
                all(
                    record["time_before_s"] < time_s <= record["time_after_s"]
                    for time_s in record["accepted_internal_time_nodes_s"]
                )
                and all(
                    right > left
                    for left, right in zip(
                        record["accepted_internal_time_nodes_s"],
                        record["accepted_internal_time_nodes_s"][1:],
                    )
                )
                and math.isclose(
                    record["accepted_internal_time_nodes_s"][-1],
                    record["time_after_s"],
                    rel_tol=0.0,
                    abs_tol=1.0e-15,
                )
            )
        )
        for record in step_records
    )
    return {
        "requested_dt_s": dt_s,
        "expected_accepted_step_count": expected_count,
        "actual_accepted_step_count": len(step_records),
        "every_accepted_step_has_requested_dt": durations_match,
        "accepted_times_match_requested_grid": times_match,
        "stage_sequence_matches_schedule": stages_match,
        "square_transitions_align_with_grid": transitions_aligned,
        "adaptive_subdivision_used": adaptive_subdivision_used,
        # Retain the historical diagnostic key for consumers that distinguish
        # a direct requested-grid solve from the now-licensed atomic
        # failed-leaf subdivision route.
        "no_accepted_step_subdivision": (durations_match and not adaptive_subdivision_used),
        "internal_nodes_respect_macro_intervals": (internal_nodes_respect_macro_intervals),
        "reporting_time_comparison_abs_tolerance_s": 1.0e-15,
        "exact_requested_time_grid": passed and internal_nodes_respect_macro_intervals,
        "frozen_failed_leaf_controller_contract_passed": (
            passed and frozen_failed_leaf_controller_contract_passed
        ),
        "accepted_subdivision_policy": ADAPTIVE_FAILED_LEAF_POLICY_ID,
        "arbitrary_accepted_subdivision_allowed": False,
    }


def _relative_transition_response(
    step_records: Sequence[dict],
    scenario: StressScenario,
    field: str,
    *,
    pulse_direction: int,
) -> dict:
    """Compare accepted endpoint samples without interpolation or a dead band."""

    if pulse_direction not in (-1, 1):
        raise ValueError("pulse direction must be -1 or +1")
    baseline_stage, pulse_stage, return_stage = (item[0] for item in scenario.segments)
    groups = {
        stage: [record for record in step_records if record["stage"] == stage]
        for stage in (baseline_stage, pulse_stage, return_stage)
    }
    if not all(groups.values()):
        return {
            "all_finite": False,
            "pulse_minus_baseline": None,
            "return_minus_pulse": None,
            "pulse_response_in_declared_direction": False,
            "return_response_in_opposite_direction": False,
        }
    baseline_value = groups[baseline_stage][-1][field]
    pulse_value = groups[pulse_stage][0][field]
    pulse_end_value = groups[pulse_stage][-1][field]
    return_value = groups[return_stage][0][field]
    finite = all(
        math.isfinite(value)
        for value in (baseline_value, pulse_value, pulse_end_value, return_value)
    )
    pulse_delta = pulse_value - baseline_value
    return_delta = return_value - pulse_end_value
    return {
        "field": field,
        "baseline_end": baseline_value,
        "pulse_first": pulse_value,
        "pulse_end": pulse_end_value,
        "return_first": return_value,
        "pulse_minus_baseline": pulse_delta,
        "return_minus_pulse": return_delta,
        "all_finite": finite,
        "strict_sign_only_no_dead_band": True,
        "pulse_response_in_declared_direction": (finite and pulse_direction * pulse_delta > 0.0),
        "return_response_in_opposite_direction": (finite and pulse_direction * return_delta < 0.0),
    }


def _temperature_only_thermal_driving_force_response(
    step_records: Sequence[dict],
    scenario: StressScenario,
    *,
    expected_segment_duration_s: float | None = None,
) -> dict:
    """Recompute the Amendment-14 thermal oracle from serialized film primitives.

    A particle may continue warming after the boundary is made cooler when the
    boundary remains hotter than the particle.  The separated temperature-only
    experiment therefore gates the signed boundary-to-surface driving force,
    not the sign of the absolute surface-temperature response.  No serialized
    derived difference, heat response, or pass flag is trusted here.
    """

    segment_duration_s = (
        SEGMENT_DURATION_S if expected_segment_duration_s is None else expected_segment_duration_s
    )
    if (
        isinstance(segment_duration_s, bool)
        or not isinstance(segment_duration_s, (int, float))
        or not math.isfinite(float(segment_duration_s))
        or float(segment_duration_s) <= 0.0
    ):
        raise ValueError("Amendment-14 expected segment duration must be positive")
    segment_duration_s = float(segment_duration_s)

    applicability = {
        "amendment_id": TEMPERATURE_ONLY_THERMAL_ORACLE_AMENDMENT_ID,
        "scope": "separated_temperature_only_square_pulse",
        "definition": TEMPERATURE_ONLY_THERMAL_DRIVING_FORCE_DEFINITION,
        "pulse_requirement": "cooler_pulse_strictly_reduces_delta_temperature",
        "return_requirement": ("exact_baseline_return_strictly_increases_delta_temperature"),
        "production_heat_flux_definition": ("q_out=h*Ackermann_factor*(T_film_surface-T_boundary)"),
        "production_heat_flux_sign_convention": "positive_outward_from_particle",
        "production_heat_response_requirement": (
            "cooler_pulse_strictly_reduces_abs_q_out_and_exact_baseline_return_"
            "strictly_increases_abs_q_out"
        ),
        "outer_dry_cell_center_temperature_direction_prescribed": False,
        "water_storage_direction_prescribed": False,
        "strict_sign_only_no_dead_band": True,
        "physics_or_numerical_equations_changed": False,
        "physically_qualifying": False,
    }
    failure: dict = {
        "applicability": applicability,
        "serialized_input_fields": [
            "stage",
            "time_before_s",
            "time_after_s",
            "dt_s",
            "boundary_temperature_k",
            "surface_film_audit.bulk_temperature_k",
            "surface_film_audit.surface_temperature_k",
            "surface_film_audit.heat_transfer_w_m2_k",
            "surface_film_audit.ackermann.factor",
            "surface_film_audit.production_conductive_heat_flux_w_m2",
            "surface_film_audit.heat_closure_roundoff_bound_w_m2",
        ],
        "derived_values_recomputed_not_trusted": True,
        "all_serialized_inputs_finite": False,
        "every_boundary_matches_declared_stage_exactly": False,
        "accepted_stage_sequence_matches_declared_contiguous_order": False,
        "expected_contiguous_stage_sequence": None,
        "observed_contiguous_stage_sequence": None,
        "every_record_has_strict_chronological_interval": False,
        "every_record_time_increment_matches_dt_within_roundoff": False,
        "maximum_time_increment_residual_s": None,
        "maximum_time_increment_roundoff_bound_s": None,
        "accepted_record_times_are_exactly_contiguous": False,
        "every_stage_covers_its_exact_declared_interval": False,
        "expected_stage_duration_s": segment_duration_s,
        "observed_stage_durations_s": None,
        "observed_stage_time_intervals_s": None,
        "every_film_bulk_temperature_matches_boundary_exactly": False,
        "every_production_heat_flux_reconstruction_within_roundoff_bound": False,
        "every_production_heat_flux_has_constitutive_sign": False,
        "maximum_heat_flux_reconstruction_residual_w_m2": None,
        "maximum_heat_flux_reconstruction_roundoff_bound_w_m2": None,
        "exact_baseline_boundary_return": False,
        "pulse_boundary_is_strictly_cooler": False,
        "samples": None,
        "pulse_minus_baseline_delta_temperature_k": None,
        "return_minus_pulse_delta_temperature_k": None,
        "pulse_strictly_reduced_driving_force": False,
        "return_strictly_increased_driving_force": False,
        "pulse_strictly_reduced_absolute_inward_heat_flux": False,
        "return_strictly_increased_absolute_inward_heat_flux": False,
        "oracle_passed": False,
        "validation_errors": [],
    }
    errors: list[str] = []
    if scenario.identifier != SCENARIO_TEMPERATURE_ONLY:
        errors.append("Amendment-14 oracle is applicable only to temperature_only")
        failure["validation_errors"] = errors
        return failure

    segments = scenario.segments
    expected_boundaries = {stage: boundary.temperature_k for stage, boundary in segments}
    expected_stages = tuple(expected_boundaries)
    groups: dict[str, list[dict]] = {stage: [] for stage in expected_stages}
    all_inputs_finite = bool(step_records)
    boundaries_match = bool(step_records)
    film_bulk_matches = bool(step_records)
    heat_flux_closure_passed = bool(step_records)
    heat_flux_sign_passed = bool(step_records)
    heat_flux_residuals: list[float] = []
    heat_flux_roundoff_bounds: list[float] = []
    observed_stage_sequence: list[str] = []
    chronological_intervals = bool(step_records)
    exact_time_increments = bool(step_records)
    contiguous_times = bool(step_records)
    previous_time_after: float | None = None
    time_increment_residuals: list[float] = []
    time_increment_roundoff_bounds: list[float] = []

    def finite_number(value: object) -> bool:
        return (
            not isinstance(value, bool)
            and isinstance(value, (int, float))
            and math.isfinite(float(value))
        )

    for index, record in enumerate(step_records):
        if not isinstance(record, Mapping):
            errors.append(f"accepted record {index} is not an object")
            all_inputs_finite = False
            boundaries_match = False
            film_bulk_matches = False
            heat_flux_closure_passed = False
            heat_flux_sign_passed = False
            continue
        stage = record.get("stage")
        if not isinstance(stage, str) or stage not in groups:
            errors.append(f"accepted record {index} has an unexpected stage")
            all_inputs_finite = False
            boundaries_match = False
            film_bulk_matches = False
            heat_flux_closure_passed = False
            heat_flux_sign_passed = False
            continue
        if not observed_stage_sequence or observed_stage_sequence[-1] != stage:
            observed_stage_sequence.append(stage)
        time_before = record.get("time_before_s")
        time_after = record.get("time_after_s")
        dt_s = record.get("dt_s")
        finite_times = all(finite_number(value) for value in (time_before, time_after, dt_s))
        chronological = (
            finite_times
            and float(time_before) >= 0.0
            and float(time_after) > float(time_before)
            and float(dt_s) > 0.0
        )
        target_time = math.fsum((float(time_before), float(dt_s))) if chronological else None
        time_increment_residual = (
            math.fsum((float(time_after), -target_time)) if target_time is not None else None
        )
        time_increment_roundoff = (
            8.0
            * max(
                math.ulp(max(abs(float(time_before)), 1.0)),
                math.ulp(max(abs(float(time_after)), 1.0)),
                math.ulp(max(abs(float(dt_s)), 1.0)),
            )
            if chronological
            else None
        )
        increment_matches = (
            chronological and abs(time_increment_residual) <= time_increment_roundoff
        )
        if time_increment_residual is not None:
            time_increment_residuals.append(abs(time_increment_residual))
        if time_increment_roundoff is not None:
            time_increment_roundoff_bounds.append(time_increment_roundoff)
        contiguous = chronological and (
            previous_time_after is None or float(time_before) == previous_time_after
        )
        if not chronological:
            errors.append(f"accepted record {index} has a nonchronological time interval")
            chronological_intervals = False
        if not finite_times:
            all_inputs_finite = False
        if not increment_matches:
            errors.append(
                f"accepted record {index} time increment contradicts its dt beyond roundoff"
            )
            exact_time_increments = False
        if not contiguous:
            errors.append(f"accepted record {index} is not exactly contiguous in time")
            contiguous_times = False
        if chronological:
            previous_time_after = float(time_after)
        boundary = record.get("boundary_temperature_k")
        film = record.get("surface_film_audit")
        ackermann = film.get("ackermann") if isinstance(film, Mapping) else None
        film_bulk = film.get("bulk_temperature_k") if isinstance(film, Mapping) else None
        film_surface = film.get("surface_temperature_k") if isinstance(film, Mapping) else None
        heat_transfer = film.get("heat_transfer_w_m2_k") if isinstance(film, Mapping) else None
        ackermann_factor = ackermann.get("factor") if isinstance(ackermann, Mapping) else None
        production_heat_flux = (
            film.get("production_conductive_heat_flux_w_m2") if isinstance(film, Mapping) else None
        )
        closure_bound = (
            film.get("heat_closure_roundoff_bound_w_m2") if isinstance(film, Mapping) else None
        )
        finite = all(
            finite_number(value)
            for value in (
                boundary,
                film_bulk,
                film_surface,
                heat_transfer,
                ackermann_factor,
                production_heat_flux,
                closure_bound,
            )
        )
        positive_coefficients = (
            finite
            and float(heat_transfer) > 0.0
            and float(ackermann_factor) > 0.0
            and float(closure_bound) >= 0.0
        )
        if not positive_coefficients:
            errors.append(f"accepted record {index} has incomplete or inadmissible film primitives")
            all_inputs_finite = False
        if boundary != expected_boundaries[stage]:
            errors.append(f"accepted record {index} boundary contradicts its declared stage")
            boundaries_match = False
        if film_bulk != boundary:
            errors.append(f"accepted record {index} film bulk temperature differs from boundary")
            film_bulk_matches = False
        if not positive_coefficients or not finite_times:
            heat_flux_closure_passed = False
            heat_flux_sign_passed = False
            continue
        delta_temperature = math.fsum((float(boundary), -float(film_surface)))
        expected_heat_flux = (
            float(heat_transfer)
            * float(ackermann_factor)
            * math.fsum((float(film_surface), -float(film_bulk)))
        )
        heat_flux_residual = math.fsum((float(production_heat_flux), -expected_heat_flux))
        heat_flux_residuals.append(abs(heat_flux_residual))
        heat_flux_roundoff_bounds.append(float(closure_bound))
        closure_passed = abs(heat_flux_residual) <= float(closure_bound)
        sign_passed = (
            (delta_temperature > 0.0 and float(production_heat_flux) < 0.0)
            or (delta_temperature < 0.0 and float(production_heat_flux) > 0.0)
            or (delta_temperature == 0.0 and float(production_heat_flux) == 0.0)
        )
        if not closure_passed:
            errors.append(f"accepted record {index} production heat flux fails reconstruction")
            heat_flux_closure_passed = False
        if not sign_passed:
            errors.append(f"accepted record {index} production heat flux has the wrong sign")
            heat_flux_sign_passed = False
        groups[stage].append(
            {
                "stage": stage,
                "time_before_s": float(time_before),
                "time_after_s": float(time_after),
                "dt_s": float(dt_s),
                "boundary_temperature_k": float(boundary),
                "film_surface_temperature_k": float(film_surface),
                "delta_temperature_k": delta_temperature,
                "heat_transfer_w_m2_k": float(heat_transfer),
                "ackermann_factor": float(ackermann_factor),
                "production_conductive_heat_flux_outward_w_m2": (float(production_heat_flux)),
                "reconstructed_production_conductive_heat_flux_outward_w_m2": (expected_heat_flux),
                "heat_flux_reconstruction_residual_w_m2": heat_flux_residual,
                "heat_flux_reconstruction_roundoff_bound_w_m2": (float(closure_bound)),
                "heat_flux_reconstruction_within_roundoff_bound": closure_passed,
                "heat_flux_sign_matches_driving_force": sign_passed,
            }
        )

    for stage in expected_stages:
        if not groups[stage]:
            errors.append(f"accepted records omit stage {stage}")
    baseline_stage, pulse_stage, return_stage = expected_stages
    baseline_boundary = expected_boundaries[baseline_stage]
    pulse_boundary = expected_boundaries[pulse_stage]
    return_boundary = expected_boundaries[return_stage]
    exact_return = baseline_boundary == return_boundary
    pulse_is_cooler = pulse_boundary < baseline_boundary
    stage_sequence_matches = tuple(observed_stage_sequence) == expected_stages
    observed_stage_intervals = {
        stage: (
            [
                group[0]["time_before_s"],
                group[-1]["time_after_s"],
            ]
            if group
            else None
        )
        for stage, group in groups.items()
    }
    observed_stage_durations = {
        stage: math.fsum(item["dt_s"] for item in group) for stage, group in groups.items()
    }
    exact_stage_intervals = bool(groups[baseline_stage]) and (
        groups[baseline_stage][0]["time_before_s"] == 0.0
        and all(observed_stage_durations[stage] == segment_duration_s for stage in expected_stages)
    )

    failure.update(
        all_serialized_inputs_finite=all_inputs_finite,
        every_boundary_matches_declared_stage_exactly=boundaries_match,
        accepted_stage_sequence_matches_declared_contiguous_order=(stage_sequence_matches),
        expected_contiguous_stage_sequence=list(expected_stages),
        observed_contiguous_stage_sequence=list(observed_stage_sequence),
        every_record_has_strict_chronological_interval=chronological_intervals,
        every_record_time_increment_matches_dt_within_roundoff=(exact_time_increments),
        maximum_time_increment_residual_s=(
            max(time_increment_residuals) if time_increment_residuals else None
        ),
        maximum_time_increment_roundoff_bound_s=(
            max(time_increment_roundoff_bounds) if time_increment_roundoff_bounds else None
        ),
        accepted_record_times_are_exactly_contiguous=contiguous_times,
        every_stage_covers_its_exact_declared_interval=exact_stage_intervals,
        expected_stage_duration_s=segment_duration_s,
        observed_stage_durations_s=observed_stage_durations,
        observed_stage_time_intervals_s=observed_stage_intervals,
        every_film_bulk_temperature_matches_boundary_exactly=film_bulk_matches,
        every_production_heat_flux_reconstruction_within_roundoff_bound=(heat_flux_closure_passed),
        every_production_heat_flux_has_constitutive_sign=heat_flux_sign_passed,
        maximum_heat_flux_reconstruction_residual_w_m2=(
            max(heat_flux_residuals) if heat_flux_residuals else None
        ),
        maximum_heat_flux_reconstruction_roundoff_bound_w_m2=(
            max(heat_flux_roundoff_bounds) if heat_flux_roundoff_bounds else None
        ),
        exact_baseline_boundary_return=exact_return,
        pulse_boundary_is_strictly_cooler=pulse_is_cooler,
    )
    if not exact_return:
        errors.append("third-stage boundary is not the exact baseline return")
    if not pulse_is_cooler:
        errors.append("temperature-only pulse boundary is not strictly cooler")
    if not stage_sequence_matches:
        errors.append("accepted records do not follow the contiguous declared stage order")
    if not exact_stage_intervals:
        errors.append("accepted records do not cover the exact declared stage intervals")
    if errors:
        failure["validation_errors"] = errors
        return failure

    selected = {
        "baseline_end": groups[baseline_stage][-1],
        "pulse_first": groups[pulse_stage][0],
        "pulse_end": groups[pulse_stage][-1],
        "return_first": groups[return_stage][0],
    }
    samples = selected
    pulse_change = math.fsum(
        (
            samples["pulse_first"]["delta_temperature_k"],
            -samples["baseline_end"]["delta_temperature_k"],
        )
    )
    return_change = math.fsum(
        (
            samples["return_first"]["delta_temperature_k"],
            -samples["pulse_end"]["delta_temperature_k"],
        )
    )
    pulse_reduced = pulse_change < 0.0
    return_increased = return_change > 0.0
    pulse_heat_reduced = abs(
        samples["pulse_first"]["production_conductive_heat_flux_outward_w_m2"]
    ) < abs(samples["baseline_end"]["production_conductive_heat_flux_outward_w_m2"])
    return_heat_increased = abs(
        samples["return_first"]["production_conductive_heat_flux_outward_w_m2"]
    ) > abs(samples["pulse_end"]["production_conductive_heat_flux_outward_w_m2"])
    failure.update(
        samples=samples,
        pulse_minus_baseline_delta_temperature_k=pulse_change,
        return_minus_pulse_delta_temperature_k=return_change,
        pulse_strictly_reduced_driving_force=pulse_reduced,
        return_strictly_increased_driving_force=return_increased,
        pulse_strictly_reduced_absolute_inward_heat_flux=pulse_heat_reduced,
        return_strictly_increased_absolute_inward_heat_flux=(return_heat_increased),
        oracle_passed=(
            pulse_reduced and return_increased and pulse_heat_reduced and return_heat_increased
        ),
        validation_errors=[],
    )
    return failure


def _nonzero_storage_response(
    step_records: Sequence[dict],
    scenario: StressScenario,
) -> dict:
    """Audit aggregate equilibrium-storage changes without assigning a cause."""

    fields = (
        "retained_water_inventory_mol",
        "pore_vapor_water_inventory_mol",
    )
    diagnostics = {
        field: _relative_transition_response(
            step_records,
            scenario,
            field,
            pulse_direction=-1,
        )
        for field in fields
    }
    finite = all(item["all_finite"] for item in diagnostics.values())
    pulse_changed = finite and any(
        item["pulse_minus_baseline"] != 0.0 for item in diagnostics.values()
    )
    return_changed = finite and any(
        item["return_minus_pulse"] != 0.0 for item in diagnostics.values()
    )
    return {
        "inventories": diagnostics,
        "all_finite": finite,
        "pulse_changed_at_machine_resolution": pulse_changed,
        "return_changed_at_machine_resolution": return_changed,
        "direction_is_not_prescribed": True,
        "interpretation": (
            "retained and dry-pore-vapor inventories are aggregate state-function "
            "diagnostics; nonzero response is gated, but no flux mechanism or "
            "monotone sign is inferred from them"
        ),
    }


def _surface_constitutive_oracle_contract(
    step_records: Sequence[dict],
    scenario: StressScenario,
) -> dict:
    """Recompute A15 from every raw root and four chronological endpoints."""

    errors: list[str] = []
    accepted_root_count = 0
    every_root_passed = True
    every_serialized_audit_matched = True
    for record_index, record in enumerate(step_records):
        samples = record.get("accepted_node_surface_flux_samples")
        if not isinstance(samples, list) or not samples:
            errors.append(f"step record {record_index} has no accepted-root surface samples")
            every_root_passed = False
            continue
        accepted_root_count += len(samples)
        for sample_index, sample in enumerate(samples):
            if not isinstance(sample, Mapping):
                errors.append(f"step record {record_index} root {sample_index} is not an object")
                every_root_passed = False
                continue
            raw = sample.get("surface_constitutive_primitives")
            recomputed = surface_oracle.audit_surface_constitutive_sample(raw)
            if recomputed["passed"] is not True:
                errors.append(f"step record {record_index} root {sample_index} failed A15")
                every_root_passed = False
            if sample.get("surface_constitutive_oracle") != recomputed:
                errors.append(
                    f"step record {record_index} root {sample_index} changed its A15 audit"
                )
                every_serialized_audit_matched = False
            if sample_index == len(samples) - 1 and (
                record.get("surface_constitutive_primitives") != raw
                or record.get("surface_constitutive_oracle") != recomputed
            ):
                errors.append(f"step record {record_index} endpoint A15 sample is detached")
                every_serialized_audit_matched = False

    baseline_stage, pulse_stage, return_stage = (item[0] for item in scenario.segments)
    grouped = {
        stage: [record for record in step_records if record.get("stage") == stage]
        for stage in (baseline_stage, pulse_stage, return_stage)
    }
    history_samples: dict[str, object] | None = None
    history_audit: dict[str, object] | None = None
    history_applicable = scenario.identifier in (
        SCENARIO_COMPOSITION_ONLY,
        SCENARIO_COMBINED,
    )
    if not all(grouped.values()):
        errors.append("A15 history roles do not cover all three forcing stages")
    else:
        history_samples = {
            "baseline_end": grouped[baseline_stage][-1].get("surface_constitutive_primitives"),
            "pulse_first": grouped[pulse_stage][0].get("surface_constitutive_primitives"),
            "pulse_end": grouped[pulse_stage][-1].get("surface_constitutive_primitives"),
            "return_first": grouped[return_stage][0].get("surface_constitutive_primitives"),
        }
        history_audit = surface_oracle.audit_history_aware_flux_response(history_samples)

        def _r15_pulse_end_flux() -> float | None:
            samples_by_role = history_audit.get("sample_audits") or {}
            reconstruction = (
                (samples_by_role.get("pulse_end") or {}).get("reconstruction")
            )
            if not isinstance(reconstruction, Mapping):
                return None
            flux = reconstruction.get("independent_hexane_flux_mol_m2_s")
            return flux if isinstance(flux, float) else None

        def _r15_acceptable() -> bool:
            if not (
                history_audit["baseline_end_outward"] is True
                and history_audit["pulse_strictly_decreased_J_h"] is True
                and history_audit["return_strictly_increased_J_h"] is True
                and history_audit["pulse_first_inward"] is False
            ):
                return False
            # GT-PS-2-P1E-R16 (revised): a single pulse sample sits AT the
            # pulse end, where the crossing IS observable — the R15
            # pulse-end criterion therefore applies without any fallback.
            # (The first R16 draft accepted single-sample cells on the
            # invariant alone; the joint selector then picked a candidate
            # whose pulse-end flux is measurably outward, and gate1h's
            # finer grids caught it.  Temporal-resolution disclosure for
            # coarse auxiliary cells lives in the gate1h assessment, not
            # here.)
            pulse_end_flux = _r15_pulse_end_flux()
            return pulse_end_flux is not None and pulse_end_flux < 0.0

        if history_applicable and history_audit["passed"] is not True:
            if _r15_acceptable():
                pass
            elif A15_FULL_REVERSAL_GATING:
                errors.append("the history-aware A15 J_h response failed")
            elif not (
                history_audit["pulse_strictly_decreased_J_h"] is True
                and history_audit["baseline_end_outward"] is True
            ):
                errors.append(
                    "the R11 stress-corner A15 invariant failed: the pulse did "
                    "not strictly decrease J_h from a strictly outward baseline"
                )

    def _pulse_end_crossing() -> bool:
        # GT-PS-2-P1E-R15: at fine dt the first accepted node after the
        # pulse switch samples the film transient BEFORE the inward
        # crossing completes; the crossing itself is demonstrated by the
        # pulse-end sample.  GT-PS-2-P1E-R16: with exactly ONE pulse
        # sample (pulse_first IS pulse_end) the crossing is temporally
        # unobservable, and the measured invariant legs gate instead.
        if history_audit is None:
            return False
        if not (
            history_audit["baseline_end_outward"] is True
            and history_audit["pulse_strictly_decreased_J_h"] is True
            and history_audit["return_strictly_increased_J_h"] is True
            and history_audit["pulse_first_inward"] is False
        ):
            return False
        samples = history_audit.get("sample_audits") or {}
        end_rec = (samples.get("pulse_end") or {}).get("reconstruction")
        if not isinstance(end_rec, Mapping):
            return False
        end_flux = end_rec.get("independent_hexane_flux_mol_m2_s")
        return isinstance(end_flux, float) and end_flux < 0.0

    def _history_acceptable() -> bool:
        if not history_applicable:
            return True
        if history_audit is None:
            return False
        if history_audit["passed"] is True:
            return True
        if _pulse_end_crossing():
            return True
        # GT-PS-2-P1E-R11: at stress-varied film corners the full reversal
        # demonstration is measured infeasible for every chart-admissible
        # composition (floor above ceiling at lambda=0.5, return film-lag at
        # lambda=2); the gate there is the measured universal invariant.
        return (
            not A15_FULL_REVERSAL_GATING
            and history_audit["pulse_strictly_decreased_J_h"] is True
            and history_audit["baseline_end_outward"] is True
        )

    passed = (
        accepted_root_count > 0
        and every_root_passed
        and every_serialized_audit_matched
        and history_samples is not None
        and _history_acceptable()
        and not errors
    )
    return {
        "amendment_id": surface_oracle.AMENDMENT_ID,
        "raw_sample_schema": surface_oracle.raw_sample_schema(),
        "accepted_root_count": accepted_root_count,
        "every_accepted_root_reconstructed_from_raw_primitives": (every_root_passed),
        "every_serialized_root_audit_matches_reconstruction": (every_serialized_audit_matched),
        "history_response_applicable": history_applicable,
        "history_samples": history_samples,
        "history_audit": history_audit,
        "r11_full_reversal_gating": A15_FULL_REVERSAL_GATING,
        "r11_history_acceptance": (
            "full_reversal_demonstration"
            if A15_FULL_REVERSAL_GATING
            else "stress_corner_pulse_decrease_invariant"
        ),
        # R15 disclosure is added ONLY when the pulse-end-crossing branch
        # fires, so contracts produced before R15 recompute byte-identically.
        **(
            {"r15_pulse_end_crossing_accepted": True}
            if (
                history_applicable
                and history_audit is not None
                and history_audit["passed"] is not True
                and _pulse_end_crossing()
            )
            else {}
        ),
        "absolute_return_sign_prescribed": False,
        "N_h_or_N_t_substitution_allowed": False,
        "physics_or_numerical_equations_changed": False,
        "engineering_foundation_only": True,
        "physically_qualifying": False,
        "plant_predictive": False,
        "normative_F2_solvent_particle": False,
        "normative_F3_water_interface": False,
        "contract_passed": passed,
        "errors": errors,
    }


def _mechanism_specific_contract(
    scenario: StressScenario,
    step_records: Sequence[dict],
    hexane_component_response: dict,
    independent_hexane_response: dict,
    *,
    temperature_only_segment_duration_s: float | None = None,
) -> dict:
    """Apply only the response gates justified by the scenario forcing."""

    surface_gas_entropy = [
        record["surface_binary_gas_entropy_production_w_m3_k"] for record in step_records
    ]
    entropy_pass = bool(surface_gas_entropy) and all(
        math.isfinite(value) and value >= 0.0 for value in surface_gas_entropy
    )
    constitutive = _surface_constitutive_oracle_contract(step_records, scenario)
    if scenario.identifier == SCENARIO_COMPOSITION_ONLY:
        relative = _relative_transition_response(
            step_records,
            scenario,
            "surface_independent_hexane_flux_mol_m2_s",
            pulse_direction=-1,
        )
        passed = constitutive["contract_passed"] and entropy_pass
        return {
            "scenario": scenario.identifier,
            "gated_mechanism": "independent_Maxwell_Stefan_hexane_flux_J_h",
            "j_hexane_direction_response": independent_hexane_response,
            "j_hexane_relative_transition_response": relative,
            "legacy_absolute_return_sign_diagnostic_is_gate": False,
            "surface_constitutive_oracle_contract": constitutive,
            "surface_gas_entropy_production": {
                "minimum_w_m3_k": min(surface_gas_entropy),
                "nonnegative_without_tolerance": entropy_pass,
                "dimensionless_exchange_force_gradient_sign_follows_J_h_for_positive_mobility": True,
                "force_sign_relation": (
                    "J_h = L*d[(mu_w-mu_h)/(R*T)]/dr for the positive scalar "
                    "mobility L used here, so J_h follows that exchange-force "
                    "gradient and R*L*gradient^2 is nonnegative"
                ),
            },
            "conserved_hexane_component_flux_N_h_is_direction_gate": False,
            "total_stefan_flux_N_t_is_direction_gate": False,
            "mechanism_contract_passed": passed,
        }
    if scenario.identifier == SCENARIO_TEMPERATURE_ONLY:
        absolute_temperature = _relative_transition_response(
            step_records,
            scenario,
            "surface_dry_temperature_k",
            pulse_direction=-1,
        )
        absolute_temperature.update(
            direction_not_prescribed=True,
            directional_response_flags_are_diagnostic_only=True,
            diagnostic_location=("outer_dry_cell_center_not_resolved_film_surface"),
            interpretation=(
                "thermal history can make the outer dry-cell centre continue "
                "warming under a cooler boundary that remains hotter than the "
                "resolved film surface; Amendment 14 therefore does not gate the "
                "sign of absolute cell-centre temperature"
            ),
        )
        thermal_driving_force = _temperature_only_thermal_driving_force_response(
            step_records,
            scenario,
            expected_segment_duration_s=temperature_only_segment_duration_s,
        )
        storage = _nonzero_storage_response(step_records, scenario)
        passed = (
            thermal_driving_force["oracle_passed"]
            and storage["pulse_changed_at_machine_resolution"]
            and storage["return_changed_at_machine_resolution"]
            and constitutive["contract_passed"]
            and entropy_pass
        )
        return {
            "scenario": scenario.identifier,
            "gated_mechanism": (
                "signed_boundary_to_film_surface_thermal_driving_force_exact_"
                "production_heat_flux_and_local_equilibrium_storage_response"
            ),
            "evidence_oracle_applicability": thermal_driving_force["applicability"],
            "boundary_to_film_surface_thermal_driving_force_response": (thermal_driving_force),
            "surface_dry_temperature_response": absolute_temperature,
            "water_storage_response": storage,
            "surface_constitutive_oracle_contract": constitutive,
            "Maxwell_Stefan_attribution_made": False,
            "hexane_J_h_direction_is_gate": False,
            "conserved_hexane_component_flux_N_h_is_direction_gate": False,
            "surface_binary_gas_entropy_production_nonnegative_without_tolerance": (entropy_pass),
            "mechanism_contract_passed": passed,
        }
    if scenario.identifier == SCENARIO_COMBINED:
        relative_j = _relative_transition_response(
            step_records,
            scenario,
            "surface_independent_hexane_flux_mol_m2_s",
            pulse_direction=-1,
        )
        passed = constitutive["contract_passed"] and entropy_pass
        return {
            "scenario": scenario.identifier,
            "gated_mechanism": "independent_Maxwell_Stefan_hexane_flux_J_h",
            "conserved_hexane_component_N_h_response": hexane_component_response,
            "conserved_hexane_component_flux_N_h_is_direction_gate": False,
            "j_hexane_direction_response": independent_hexane_response,
            "j_hexane_relative_transition_response": relative_j,
            "legacy_absolute_return_sign_diagnostic_is_gate": False,
            "surface_constitutive_oracle_contract": constitutive,
            "surface_gas_entropy_production_nonnegative_without_tolerance": (entropy_pass),
            "total_stefan_flux_N_t_is_direction_gate": False,
            "mechanism_contract_passed": passed,
        }
    raise ValueError(f"unknown Gate-1g stress scenario: {scenario.identifier}")


def _aggregate_amendment_18_19_evidence(
    step_records: Sequence[Mapping[str, object]],
) -> dict:
    """Revalidate every accepted step/root and expose exact A18/A19 counts."""

    if not step_records:
        raise ValueError("A18/A19 job accounting needs accepted step records")
    step_identities = []
    accepted_root_count = 0
    plateau_root_count = 0
    total_bisection_iterations = 0
    maximum_bisection_iterations = 0
    residual_scale_evaluation_count = 0
    capacity_evaluation_count = 0
    dry_energy_row_count = 0
    stencil_kind_counts: dict[str, int] = {}
    evaluation_kind_counts: dict[str, int] = {}
    for step_index, record in enumerate(step_records):
        samples = record.get("accepted_node_surface_flux_samples")
        if not isinstance(samples, Sequence) or isinstance(samples, (str, bytes)) or not samples:
            raise ValueError("accepted step lost its accepted-root evidence samples")
        a18 = _interface_composition_root_audit_summary(samples)
        a19 = _dry_energy_capacity_stencil_audit_summary(samples)
        if not (
            record.get("interface_saturation_endpoint_audit_summary") == a18
            and record.get("dry_energy_capacity_stencil_audit_summary") == a19
            and record.get("amendment_18_19_evidence_contract_passed") is True
            and record.get("accepted_internal_substep_count") == len(samples)
            and a18["accepted_root_count"] == a19["accepted_root_count"] == len(samples)
            and a18["contract_passed"] is True
            and a19["contract_passed"] is True
        ):
            raise RuntimeError("accepted step A18/A19 summary does not match its roots")
        step_identity_payload = {
            "accepted_step_record_index": step_index,
            "stage": record.get("stage"),
            "time_before_s": record.get("time_before_s"),
            "time_after_s": record.get("time_after_s"),
            "dt_s": record.get("dt_s"),
            "a18_accepted_root_identity_sha256": [
                root["accepted_root_identity_sha256"]
                for root in a18["accepted_root_identities"]
            ],
            "a19_accepted_root_scale_identity_sha256": [
                root["accepted_root_scale_identity_sha256"]
                for root in a19["accepted_root_scale_identities"]
            ],
        }
        step_identities.append(
            {
                **step_identity_payload,
                "accepted_step_evidence_identity_sha256": _canonical_evidence_sha256(
                    step_identity_payload
                ),
            }
        )
        accepted_root_count += len(samples)
        plateau_root_count += int(a18["rounding_plateau_root_count"])
        total_bisection_iterations += int(a18["total_bisection_iteration_count"])
        maximum_bisection_iterations = max(
            maximum_bisection_iterations,
            int(a18["maximum_bisection_iteration_count"]),
        )
        capacity_evaluation_count += int(a19["capacity_evaluation_audit_count"])
        residual_scale_evaluation_count += int(
            a19["residual_scale_evaluation_count"]
        )
        dry_energy_row_count += int(a19["dry_energy_row_count"])
        for kind, count in a19["stencil_kind_counts"].items():
            stencil_kind_counts[kind] = stencil_kind_counts.get(kind, 0) + int(count)
        for kind, count in a19["residual_scale_evaluation_kind_counts"].items():
            evaluation_kind_counts[kind] = evaluation_kind_counts.get(kind, 0) + int(count)
    step_hashes = [item["accepted_step_evidence_identity_sha256"] for item in step_identities]
    if len(set(step_hashes)) != len(step_hashes):
        raise RuntimeError("accepted-step A18/A19 evidence identities are not one-to-one")
    return {
        "schema_version": 1,
        "numerical_amendment_ids": [
            INTERFACE_ROOT_AMENDMENT_ID,
            DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "accepted_step_record_count": len(step_records),
        "accepted_step_evidence_identity_count": len(step_identities),
        "accepted_root_count": accepted_root_count,
        "accepted_step_evidence_identities": step_identities,
        "interface_saturation_endpoint_audit_accounting": {
            "accepted_root_count": accepted_root_count,
            "full_coupled_pore_certified_root_count": accepted_root_count,
            "rounding_plateau_root_count": plateau_root_count,
            "total_bisection_iteration_count": total_bisection_iterations,
            "maximum_bisection_iteration_count": maximum_bisection_iterations,
            "tolerance_changed_root_count": 0,
            "clipping_projection_or_nextafter_root_count": 0,
            "contract_passed": True,
        },
        "dry_energy_capacity_stencil_audit_accounting": {
            "residual_scale_evaluation_count": residual_scale_evaluation_count,
            "dry_energy_row_count": dry_energy_row_count,
            "capacity_evaluation_audit_count": capacity_evaluation_count,
            "stencil_kind_counts": dict(sorted(stencil_kind_counts.items())),
            "residual_scale_evaluation_kind_counts": dict(
                sorted(evaluation_kind_counts.items())
            ),
            "accepted_residual_equations_changed": False,
            "composition_or_activity_clipping_used": False,
            "property_extrapolation_used": False,
            "acceptance_tolerances_changed": False,
            "contract_passed": True,
        },
        "all_accepted_steps_and_roots_have_complete_consistent_evidence": True,
        "accepted_state_clipping_or_projection_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": True,
        "physically_qualifying": False,
    }


def _aggregate_amendment_20_evidence(
    step_records: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Revalidate every committed step/root actual-path evidence bundle."""

    if not step_records:
        raise ValueError("A20 job accounting needs accepted step records")
    step_identities: list[dict[str, object]] = []
    root_count = 0
    root_path_identities: list[str] = []
    for step_index, record in enumerate(step_records):
        samples = record.get("accepted_node_surface_flux_samples")
        if (
            not isinstance(samples, Sequence)
            or isinstance(samples, (str, bytes))
            or not samples
        ):
            raise ValueError("accepted step lost its A20 accepted-root samples")
        summary = _moving_interface_actual_path_audit_summary(samples)
        if not (
            record.get("moving_interface_actual_composition_force_path_audit_summary")
            == summary
            and record.get("amendment_20_evidence_contract_passed") is True
            and record.get("accepted_internal_substep_count") == len(samples)
            and summary["accepted_root_count"] == len(samples)
            and summary["accepted_root_identity_count"] == len(samples)
            and summary["all_roots_independently_replayed"] is True
            and summary["contract_passed"] is True
        ):
            raise RuntimeError("accepted step A20 summary does not match its roots")
        identities = summary["accepted_root_identities"]
        canonical_hashes = [item["canonical_binding_sha256"] for item in identities]
        # Defect 5: only finite-positive roots carry a production-path identity;
        # inapplicable topology routes are counted, not path-hashed.
        path_hashes = [
            item["production_path_dataclass_identity_sha256"]
            for item in identities
            if item.get("applicability") != "topology_route_inapplicable"
        ]
        payload = {
            "accepted_step_record_index": step_index,
            "stage": record.get("stage"),
            "time_before_s": record.get("time_before_s"),
            "time_after_s": record.get("time_after_s"),
            "dt_s": record.get("dt_s"),
            "accepted_root_count": len(samples),
            "accepted_root_canonical_binding_sha256": canonical_hashes,
            "accepted_root_path_dataclass_identity_sha256": path_hashes,
        }
        step_identities.append(
            {
                **payload,
                "accepted_step_a20_evidence_identity_sha256": (
                    _canonical_evidence_sha256(payload)
                ),
            }
        )
        root_count += len(samples)
        root_path_identities.extend(path_hashes)
    step_hashes = [
        item["accepted_step_a20_evidence_identity_sha256"] for item in step_identities
    ]
    if len(set(step_hashes)) != len(step_hashes):
        raise RuntimeError("accepted-step A20 evidence identities are not one-to-one")
    return {
        "schema_version": ACTUAL_COMPOSITION_FORCE_PATH_EVIDENCE_SCHEMA_VERSION,
        "engineering_amendment_id": ACTUAL_COMPOSITION_FORCE_PATH_AMENDMENT_ID,
        "accepted_step_record_count": len(step_records),
        "accepted_step_evidence_identity_count": len(step_identities),
        "accepted_root_count": root_count,
        "accepted_root_path_identity_count": len(root_path_identities),
        "accepted_step_evidence_identities": step_identities,
        "every_accepted_root_independently_replayed": True,
        "every_accepted_root_uses_actual_linear_T_y_path_authority": True,
        "every_accepted_root_serializes_full_conditional_continuous_phase_audit": True,
        "fictitious_thermodynamic_force_temperature_count": 0,
        "accepted_state_clipping_or_projection_used": False,
        "physics_equations_or_acceptance_tolerances_changed": False,
        "contract_passed": True,
        "physically_qualifying": False,
    }


def _attempt_evaluation_count(step_records: Sequence[dict]) -> int:
    return sum(
        attempt["nonlinear_evaluations"]
        for record in step_records
        for attempt in record["attempts"]
        if attempt["role"] != "mesh_bootstrap_exact_root"
    )


FACE_ROUTE_TIME_REFINED_PREFIX_ENVIRONMENT_VARIABLE = "DTDC_P1EF_FACE_ROUTE_TIME_REFINED_PREFIX"


class FaceRouteTimeRefinedPrefixKnobError(ValueError):
    """The time-refined prefix route knob holds a value other than unset, 0 or 1."""


class FaceRouteTimeRefinedPrefixError(RuntimeError):
    """Typed refusal: the face-route time-refined prefix route cannot be taken."""


def face_route_time_refined_prefix_armed() -> bool:
    """Read the knob at call time: unset or "0" is the on-disk path, "1" arms it.

    Any other text refuses (before any solve when read by ``run_campaign``).
    """

    value = os.environ.get(FACE_ROUTE_TIME_REFINED_PREFIX_ENVIRONMENT_VARIABLE)
    if value is None or value == "0":
        return False
    if value == "1":
        return True
    raise FaceRouteTimeRefinedPrefixKnobError(
        f"{FACE_ROUTE_TIME_REFINED_PREFIX_ENVIRONMENT_VARIABLE} must be unset, '0' or "
        f"'1', not {value!r}"
    )


def _face_route_time_refined_first_step(
    cells: int,
    dt_s: float,
    coarse_first_step: "AcceptedFaceEventMacrostep",
    bootstrap_first_steps: Mapping[int, object] | None,
) -> tuple[AcceptedSameCellStep, ci.CutIntegratorState, tuple[AttemptRecord, ...], dict]:
    """Seed the first time-refined step of a mesh whose bootstrap step took the face route.

    Knob ``DTDC_P1EF_FACE_ROUTE_TIME_REFINED_PREFIX`` only.  The coarse
    (``SEGMENT_DURATION_S``) first step of this mesh is the mesh bootstrap's
    atomic arrival-plus-departure macrostep, so its candidate is a departure
    candidate in the adjacent cell's chart and cannot be rescaled in time from
    the initial state.  The bootstrap's source mesh still has a strict
    same-cell root of the same initial state.  That root is continued in time
    on its own mesh by the unchanged ``rescale_same_before_increment`` and
    solved by the unchanged same-cell solve at the refined step; the refined
    root is carried to this mesh by the unchanged
    ``prolongate_same_front_candidate`` with the exact master-face event chart
    at the refined step, exactly as the mesh bootstrap carries the coarse root;
    the target is then solved by the unchanged same-cell solve, or by the
    unchanged arrival-plus-departure macrostep when the proposal crosses the
    face.  Only the seed path changes: every accepted root is solved and gated
    by the unchanged solvers at the unchanged tolerances, and any failure is a
    typed refusal with nothing committed.
    """

    bracket = coarse_first_step.bracket
    if not isinstance(bracket, cc.MasterFaceEventContinuation):
        raise FaceRouteTimeRefinedPrefixError(
            "face-route time-refined prefix: the coarse first step is not a "
            "mesh-bootstrap master-face macrostep"
        )
    source_cells = bracket.record.source_cell_count
    if bootstrap_first_steps is None or source_cells not in bootstrap_first_steps:
        raise FaceRouteTimeRefinedPrefixError(
            "face-route time-refined prefix: the bootstrap step of the source "
            f"mesh N={source_cells} was not handed to the prefix"
        )
    source_first = bootstrap_first_steps[source_cells]
    if not isinstance(source_first, ci.CutIntegratorStep):
        raise FaceRouteTimeRefinedPrefixError(
            "face-route time-refined prefix: the source mesh's bootstrap step is "
            "not a strict same-cell root"
        )
    radius_m = coarse_first_step.before.transport.geometry.master_grid.R
    if source_first.before.transport.geometry.master_grid.R != radius_m:
        raise FaceRouteTimeRefinedPrefixError(
            "face-route time-refined prefix: the source mesh's bootstrap step has "
            "another particle radius"
        )
    source_before = _initial_integrator(source_cells, radius_m)
    source_started = perf_counter()
    source_continuation = cc.rescale_same_before_increment(
        source_before,
        source_first.assembly.candidate,
        surface_boundary=HOT_LEAN_BOUNDARY,
        source_dt_s=SEGMENT_DURATION_S,
        target_dt_s=dt_s,
        stefan_flux_scale_floor_mol_m2_s=STEFAN_SCALE_FLOOR,
        label=f"N={source_cells} dt={SEGMENT_DURATION_S:g} to dt={dt_s:g}",
    )
    source_step, source_record = _try_step(
        source_before,
        dt_s,
        HOT_LEAN_BOUNDARY,
        source_continuation.seed,
        role="face_route_prefix_source_mesh_time_continuation_root",
        predictor_fraction=dt_s / SEGMENT_DURATION_S,
    )
    source_wall_s = perf_counter() - source_started
    route: dict = {
        "environment_variable": FACE_ROUTE_TIME_REFINED_PREFIX_ENVIRONMENT_VARIABLE,
        "armed": True,
        "report_only": True,
        "source_cells": source_cells,
        "target_cells": cells,
        "dt_s": dt_s,
        "coarse_first_step_route": "mesh_bootstrap_master_face_event",
        "coarse_first_step_event_duration_s": (
            coarse_first_step.event_macrostep.ledger.event_duration_s
        ),
        "source_root_accepted": source_step is not None,
        "source_root_nonlinear_evaluations": source_record.nonlinear_evaluations,
        "source_root_maximum_scaled_residual": source_record.maximum_scaled_residual,
        "source_root_condition_proxy": source_record.condition_proxy,
        "source_root_wall_time_s": source_wall_s,
        "source_root_front_z": (
            None if source_step is None else source_step.assembly.candidate.front_z
        ),
        "source_root_counted_in_prefix_work": True,
        "tolerances_budgets_gates_changed": False,
    }
    if source_step is None:
        raise FaceRouteTimeRefinedPrefixError(
            "face-route time-refined prefix: the source-mesh time-continuation "
            f"root failed: {source_record.error}"
        )
    if source_step.ledger.scalar_post_root_polish_audit is not None:
        raise FaceRouteTimeRefinedPrefixError(
            "face-route time-refined prefix: the source-mesh root used the scalar "
            "post-root polish, whose work this seed route cannot enter into the "
            "campaign's scalar accounting"
        )
    target_before = _initial_integrator(cells, radius_m)
    continuation = cc.prolongate_same_front_candidate(
        source_before,
        source_step.assembly.candidate,
        target_before,
        stefan_flux_scale_floor_mol_m2_s=STEFAN_SCALE_FLOOR,
        label=f"N={source_cells} to N={cells} at dt={dt_s:g}",
        master_face_event_chart=cc.MasterFaceEventChart(
            requested_dt_s=dt_s,
            label=(
                "face-route time-refined prefix master-face event route; the "
                "exact event chart the mesh-continuation refusal names"
            ),
        ),
    )
    if isinstance(continuation, cc.MasterFaceEventContinuation):
        route["target_route"] = "master_face_event"
        started = perf_counter()
        try:
            macro = event_orchestrator.advance_one_interior_face_macrostep(
                target_before,
                dt_s,
                HOT_LEAN_BOUNDARY,
                continuation.same_cell_seed,
                face_event_controls=continuation.face_event_controls,
                face_event_seed=continuation.face_event_seed,
            )
        except event_orchestrator.CutEventMacrostepError as exc:
            raise FaceRouteTimeRefinedPrefixError(
                "face-route time-refined prefix: the master-face macrostep at the "
                f"refined step refused: {exc}"
            ) from exc
        if macro.ledger.branch is not (
            event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE
        ):
            raise FaceRouteTimeRefinedPrefixError(
                "face-route time-refined prefix: the master-face route did not "
                "commit one exact arrival plus its strict departure"
            )
        first = _adapt_accepted_face_event_macrostep(continuation, macro)
        first_record = AttemptRecord(
            role="time_continuation_master_face_event_macrostep",
            predictor_fraction=dt_s / SEGMENT_DURATION_S,
            wall_time_s=perf_counter() - started,
            accepted=True,
            nonlinear_evaluations=macro.ledger.total_attempted_nonlinear_evaluations,
            maximum_scaled_residual=macro.ledger.maximum_stage_scaled_residual,
            condition_proxy=first.ledger.condition_proxy,
            rollback_identity_preserved=macro.before is target_before,
            macro_branch_status=(
                event_orchestrator.CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE.value
            ),
            attempted_dt_s=dt_s,
        )
    else:
        route["target_route"] = "same_cell"
        first, first_record = _try_step(
            target_before,
            dt_s,
            HOT_LEAN_BOUNDARY,
            continuation.seed,
            role="time_continuation_exact_root",
            predictor_fraction=dt_s / SEGMENT_DURATION_S,
        )
        if first is None:
            raise FaceRouteTimeRefinedPrefixError(
                "face-route time-refined prefix: first time-refined step failed: "
                f"{first_record.error}"
            )
    route["target_root_nonlinear_evaluations"] = first_record.nonlinear_evaluations
    route["target_root_maximum_scaled_residual"] = first_record.maximum_scaled_residual
    route["target_root_front_z"] = first.after.transport.geometry.front.z
    return first, target_before, (first_record,), route


def _prepare_baseline_prefix(
    cells: int,
    dt_s: float,
    coarse_first_step: ci.CutIntegratorStep,
    *,
    bootstrap_first_steps: Mapping[int, object] | None = None,
) -> BaselinePrefix:
    """Solve the common baseline once; no post-disturbance state enters it.

    ``bootstrap_first_steps`` is read only when the knob
    ``DTDC_P1EF_FACE_ROUTE_TIME_REFINED_PREFIX`` is armed and the coarse
    first step is a mesh-bootstrap master-face macrostep; otherwise the
    on-disk path runs unchanged.
    """

    face_route_seed_route = None
    requested_steps = round(SEGMENT_DURATION_S / dt_s)
    if not math.isclose(requested_steps * dt_s, SEGMENT_DURATION_S, abs_tol=1.0e-15):
        raise ValueError("each timestep must divide the fixed segment duration exactly")
    if dt_s == SEGMENT_DURATION_S:
        first = coarse_first_step
        before = first.before
        first_attempts = (
            AttemptRecord(
                "mesh_bootstrap_exact_root",
                1.0,
                0.0,
                True,
                first.ledger.nonlinear_evaluations,
                first.ledger.maximum_scaled_residual,
                first.ledger.condition_proxy,
                first.before is before,
            ),
        )
    elif (
        isinstance(coarse_first_step, AcceptedFaceEventMacrostep)
        and face_route_time_refined_prefix_armed()
    ):
        (
            first,
            before,
            first_attempts,
            face_route_seed_route,
        ) = _face_route_time_refined_first_step(
            cells,
            dt_s,
            coarse_first_step,
            bootstrap_first_steps,
        )
    else:
        before = _initial_integrator(
            cells,
            coarse_first_step.before.transport.geometry.master_grid.R,
        )
        continuation = cc.rescale_same_before_increment(
            before,
            coarse_first_step.assembly.candidate,
            surface_boundary=HOT_LEAN_BOUNDARY,
            source_dt_s=SEGMENT_DURATION_S,
            target_dt_s=dt_s,
            stefan_flux_scale_floor_mol_m2_s=STEFAN_SCALE_FLOOR,
            label=f"dt={SEGMENT_DURATION_S:g} to dt={dt_s:g}",
        )
        first, first_record = _try_step(
            before,
            dt_s,
            HOT_LEAN_BOUNDARY,
            continuation.seed,
            role="time_continuation_exact_root",
            predictor_fraction=dt_s / SEGMENT_DURATION_S,
        )
        if first is None:
            raise RuntimeError(f"first time-refined step failed: {first_record.error}")
        first_attempts = (first_record,)

    baseline_stage = COMBINED_SCENARIO.segments[0][0]
    step_records = [
        _step_record(
            first,
            first_attempts,
            baseline_stage,
            "time_or_mesh_continuation",
        )
    ]
    accepted_steps = [first]
    previous = first
    for _ in range(1, requested_steps):
        routed = _advance_bracketed_face_from_accepted_history(
            previous,
            dt_s,
            HOT_LEAN_BOUNDARY,
        )
        if routed is None:
            current, attempts, mode = _advance_with_recovery(
                previous,
                dt_s,
                HOT_LEAN_BOUNDARY,
                HOT_LEAN_BOUNDARY,
                is_reversal=False,
            )
        else:
            current, attempts, mode = routed
        step_records.append(_step_record(current, attempts, baseline_stage, mode))
        accepted_steps.append(current)
        previous = current
    if face_route_seed_route is not None:
        # The source-mesh root is prefix work outside the step records.
        return BaselinePrefix(
            last_step=previous,
            accepted_steps=tuple(accepted_steps),
            step_records=tuple(step_records),
            incremental_wall_time_s=math.fsum(
                (
                    *(record["wall_time_s"] for record in step_records),
                    face_route_seed_route["source_root_wall_time_s"],
                )
            ),
            incremental_nonlinear_evaluations=(
                _attempt_evaluation_count(step_records)
                + face_route_seed_route["source_root_nonlinear_evaluations"]
            ),
            face_route_time_refined_seed_route=face_route_seed_route,
        )
    return BaselinePrefix(
        last_step=previous,
        accepted_steps=tuple(accepted_steps),
        step_records=tuple(step_records),
        incremental_wall_time_s=math.fsum(record["wall_time_s"] for record in step_records),
        incremental_nonlinear_evaluations=_attempt_evaluation_count(step_records),
    )


def _run_trajectory(
    cells: int,
    dt_s: float,
    coarse_first_step: ci.CutIntegratorStep,
    scenario: StressScenario | None = None,
    baseline_prefix: BaselinePrefix | None = None,
    audit_final_state_sink: list[ci.CutIntegratorState] | None = None,
    cross_radius_suffix_seed_candidates: Sequence[cut.CutTransportUnknowns] | None = None,
    cross_radius_seed_source_m: float | None = None,
    accepted_steps_out: list[AcceptedSameCellStep] | None = None,
) -> dict:
    if scenario is None:
        scenario = COMBINED_SCENARIO
    radius_m = coarse_first_step.before.transport.geometry.master_grid.R
    _validate_radius(radius_m)
    requested_steps = round(SEGMENT_DURATION_S / dt_s)
    if not math.isclose(requested_steps * dt_s, SEGMENT_DURATION_S, abs_tol=1.0e-15):
        raise ValueError("each timestep must divide the fixed segment duration exactly")
    shared_prefix = baseline_prefix is not None
    prefix = baseline_prefix or _prepare_baseline_prefix(
        cells,
        dt_s,
        coarse_first_step,
    )
    prefix_radius = prefix.last_step.after.transport.geometry.master_grid.R
    if prefix_radius != radius_m:
        raise RuntimeError("baseline continuation changed the fixed particle radius")
    baseline_stage, _ = scenario.segments[0]
    pulse_stage, _ = scenario.segments[1]
    return_stage, _ = scenario.segments[2]
    step_records = list(prefix.step_records)
    previous = prefix.last_step
    if accepted_steps_out is not None:
        accepted_steps_out.extend(prefix.accepted_steps)
    previous_boundary = HOT_LEAN_BOUNDARY
    suffix_records: list[dict] = []
    expected_suffix_steps = 2 * requested_steps
    if cross_radius_suffix_seed_candidates is not None:
        if cross_radius_seed_source_m is None:
            raise ValueError("cross-radius suffix seeds need their source radius")
        if len(cross_radius_suffix_seed_candidates) != expected_suffix_steps:
            raise ValueError("cross-radius suffix seed count must match requested steps")
    suffix_step_index = 0
    for segment_index, (stage, boundary) in enumerate(scenario.segments[1:], start=1):
        for index in range(requested_steps):
            try:
                if cross_radius_suffix_seed_candidates is None:
                    routed = _advance_bracketed_face_from_accepted_history(
                        previous,
                        dt_s,
                        boundary,
                    )
                    if routed is None:
                        bounded_discontinuity_leaf = index == 0 and isinstance(
                            boundary,
                            ct.ReducedFilmPoreBoundary,
                        )
                        current, attempts, mode = _advance_with_recovery(
                            previous,
                            dt_s,
                            previous_boundary,
                            boundary,
                            is_reversal=index == 0,
                            # The history-tangent seed remains first.  If that
                            # exact requested-boundary root rejects, try the
                            # documented same-old-state/front-only anchor before
                            # spending the bounded predictor-backtrack sequence.
                            # Every attempt retains the identical BE equations,
                            # dt, post-jump boundary, and immutable rollback state.
                            try_current_state_anchor=isinstance(
                                boundary,
                                ct.ReducedFilmPoreBoundary,
                            ),
                            provisional_maximum_function_evaluations=(
                                DISCONTINUITY_LEAF_MAXIMUM_OPTIMIZER_FUNCTION_EVALUATIONS
                                if bounded_discontinuity_leaf
                                else None
                            ),
                            maximum_nonlinear_seed_attempts=(
                                DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
                                if bounded_discontinuity_leaf
                                else None
                            ),
                        )
                    else:
                        current, attempts, mode = routed
                else:
                    current, attempts, mode = _advance_with_cross_radius_seed(
                        previous,
                        dt_s,
                        boundary,
                        cross_radius_suffix_seed_candidates[suffix_step_index],
                        source_radius_m=float(cross_radius_seed_source_m),
                    )
            except (
                CrossRadiusSeedTrajectoryError,
                SameCellRecoveryError,
            ) as exc:
                if isinstance(exc, AcceptedHistoryFaceRouteRecoveryError):
                    current, attempts, mode = _recover_after_uncertified_face_route(
                        previous,
                        dt_s,
                        previous_boundary,
                        boundary,
                        is_reversal=index == 0,
                        failure=exc,
                    )
                    record = _step_record(
                        current,
                        attempts,
                        stage,
                        mode,
                        starts_after_boundary_discontinuity=index == 0,
                    )
                    step_records.append(record)
                    suffix_records.append(record)
                    previous = current
                    previous_boundary = boundary
                    suffix_step_index += 1
                    if accepted_steps_out is not None:
                        accepted_steps_out.append(current)
                    continue
                if isinstance(exc, SameCellRecoveryError):
                    failure_audit = exc.failure_jacobian_audit
                elif (
                    isinstance(exc, CrossRadiusSeedTrajectoryError)
                    and exc.rollback_state is not None
                    and exc.dt_s is not None
                    and exc.boundary is not None
                    and exc.seed is not None
                ):
                    failure_audit = _failure_jacobian_audit(
                        exc.rollback_state,
                        exc.dt_s,
                        exc.boundary,
                        exc.seed,
                        exc.last_step_error,
                    )
                else:
                    failure_audit = None
                current, attempts, mode = _advance_with_adaptive_same_cell_macrostep(
                    previous,
                    dt_s,
                    previous_boundary,
                    boundary,
                    is_reversal=index == 0,
                    initial_attempts=exc.attempts,
                    initial_step_error=exc.last_step_error,
                    initial_failure_audit=failure_audit,
                )
            record = _step_record(
                current,
                attempts,
                stage,
                mode,
                starts_after_boundary_discontinuity=index == 0,
            )
            step_records.append(record)
            suffix_records.append(record)
            previous = current
            previous_boundary = boundary
            suffix_step_index += 1
            if accepted_steps_out is not None:
                accepted_steps_out.append(current)

    final_inventory = ci.inventory_snapshot(previous.after.transport)
    final_water_phases = ci.water_phase_inventory_snapshot(previous.after.transport)
    maximum_step = max(record["maximum_step_ledger_residual"] for record in step_records)
    maximum_cumulative = max(
        record["maximum_cumulative_ledger_residual"] for record in step_records
    )
    maximum_residual = max(record["maximum_scaled_residual"] for record in step_records)
    maximum_water_residual = max(
        max(
            record["normalized_water_step_residual"],
            record["normalized_water_cumulative_residual"],
        )
        for record in step_records
    )
    maximum_hexane_residual = max(
        max(
            record["normalized_hexane_step_residual"],
            record["normalized_hexane_cumulative_residual"],
        )
        for record in step_records
    )
    maximum_energy_residual = max(
        max(
            record["normalized_energy_step_residual"],
            record["normalized_energy_cumulative_residual"],
        )
        for record in step_records
    )
    strictly_inside_open_bounds = all(
        record["minimum_fractional_distance_to_bound"] > 0.0 for record in step_records
    )
    controls_unchanged = all(
        record["solver_controls_match_campaign_contract"] for record in step_records
    )
    scalar_post_root_polish_accounting = _aggregate_scalar_post_root_polish_summaries(
        tuple(record["scalar_post_root_polish"] for record in step_records)
    )
    incremental_scalar_post_root_polish_accounting = _aggregate_scalar_post_root_polish_summaries(
        tuple(record["scalar_post_root_polish"] for record in suffix_records)
    )
    amendment_18_19_evidence = _aggregate_amendment_18_19_evidence(step_records)
    amendment_20_evidence = _aggregate_amendment_20_evidence(step_records)

    def stage_fluxes(stage: str, field: str) -> list[float]:
        return [record[field] for record in step_records if record["stage"] == stage]

    component_fields = {
        "water": "surface_conserved_water_flux_mol_m2_s",
        "hexane": "surface_conserved_hexane_flux_mol_m2_s",
        "hexane_independent_diffusion": ("surface_independent_hexane_flux_mol_m2_s"),
    }
    component_flux_diagnostics = {
        component: {
            stage: _signed_extrema(stage_fluxes(stage, field)) for stage, _ in scenario.segments
        }
        for component, field in component_fields.items()
    }
    hexane_pulse_response = _finite_pulse_component_response(
        stage_fluxes(baseline_stage, component_fields["hexane"]),
        stage_fluxes(pulse_stage, component_fields["hexane"]),
        stage_fluxes(return_stage, component_fields["hexane"]),
    )
    independent_hexane_pulse_response = _finite_pulse_component_response(
        stage_fluxes(
            baseline_stage,
            component_fields["hexane_independent_diffusion"],
        ),
        stage_fluxes(
            pulse_stage,
            component_fields["hexane_independent_diffusion"],
        ),
        stage_fluxes(
            return_stage,
            component_fields["hexane_independent_diffusion"],
        ),
    )
    total_stefan_flux_diagnostic = {
        stage: _signed_extrema(stage_fluxes(stage, "surface_total_stefan_flux_mol_m2_s"))
        for stage, _ in scenario.segments
    }
    accepted_grid = _accepted_grid_contract(step_records, dt_s, scenario)
    schedule_contract = _dynamic_schedule_contract(scenario)
    mechanism_contract = _mechanism_specific_contract(
        scenario,
        step_records,
        hexane_pulse_response,
        independent_hexane_pulse_response,
    )
    surface_constitutive_contract = mechanism_contract["surface_constitutive_oracle_contract"]
    binary_gas_entropy_contract_passed = all(
        math.isfinite(record["surface_binary_gas_entropy_production_w_m3_k"])
        and record["surface_binary_gas_entropy_production_w_m3_k"] >= 0.0
        for record in step_records
    )
    final = previous.after
    initial_temperatures = (
        *coarse_first_step.before.transport.wet_temperatures_k,
        *coarse_first_step.before.transport.dry_temperatures_k,
        coarse_first_step.before.last_interface_temperature_k,
    )
    initial_temperature_minimum = min(initial_temperatures)
    initial_temperature_maximum = max(initial_temperatures)
    trajectory_minimum_temperature = min(
        initial_temperature_minimum,
        *(record["minimum_particle_temperature_k"] for record in step_records),
    )
    trajectory_maximum_temperature = max(
        initial_temperature_maximum,
        *(record["maximum_particle_temperature_k"] for record in step_records),
    )
    maximum_instantaneous_radial_spread = max(
        initial_temperature_maximum - initial_temperature_minimum,
        *(record["instantaneous_radial_temperature_spread_k"] for record in step_records),
    )
    minimum_front_z = min(
        coarse_first_step.before.transport.geometry.front.z,
        *(record.get("minimum_front_z", record["front_z"]) for record in step_records),
    )
    maximum_front_z = max(
        coarse_first_step.before.transport.geometry.front.z,
        *(record.get("maximum_front_z", record["front_z"]) for record in step_records),
    )
    conserved_inward_bracket = _sign_transition_bracket(
        step_records,
        baseline_stage,
        pulse_stage,
        component_fields["hexane"],
        before_sign=1,
    )
    conserved_return_bracket = _sign_transition_bracket(
        step_records,
        pulse_stage,
        return_stage,
        component_fields["hexane"],
        before_sign=-1,
    )
    independent_inward_bracket = _sign_transition_bracket(
        step_records,
        baseline_stage,
        pulse_stage,
        component_fields["hexane_independent_diffusion"],
        before_sign=1,
    )
    independent_return_bracket = _sign_transition_bracket(
        step_records,
        pulse_stage,
        return_stage,
        component_fields["hexane_independent_diffusion"],
        before_sign=-1,
    )
    stage_records = {
        stage: [record for record in step_records if record["stage"] == stage]
        for stage, _ in scenario.segments
    }
    convergence_observables = {
        "final_front_z": final.transport.geometry.front.z,
        "final_front_position_over_radius": (final.transport.geometry.front.z ** (1.0 / 3.0)),
        "minimum_front_position_over_radius": minimum_front_z ** (1.0 / 3.0),
        "maximum_front_position_over_radius": maximum_front_z ** (1.0 / 3.0),
        "final_interface_temperature_k": final.last_interface_temperature_k,
        "minimum_material_temperature_k": trajectory_minimum_temperature,
        "peak_particle_temperature_k": trajectory_maximum_temperature,
        "maximum_instantaneous_radial_temperature_spread_k": (maximum_instantaneous_radial_spread),
        "final_retained_water_mol": final_water_phases.retained_water_mol,
        "final_pore_vapor_water_mol": final_water_phases.pore_vapor_water_mol,
        "final_free_liquid_water_mol": final_water_phases.free_liquid_water_mol,
        "final_total_water_mol": final_inventory.total_water_mol,
        "final_total_hexane_mol": final_inventory.total_hexane_mol,
        "final_total_energy_j": final_inventory.total_energy_j,
        "cumulative_boundary_hexane_out_mol": (final.corrected_cumulative_boundary_hexane_out_mol),
        "cumulative_boundary_energy_out_j": (final.corrected_cumulative_boundary_energy_out_j),
    }
    endpoint_roles = (
        ("baseline_end", stage_records[baseline_stage][-1]),
        ("pulse_first", stage_records[pulse_stage][0]),
        ("pulse_end", stage_records[pulse_stage][-1]),
        ("return_first", stage_records[return_stage][0]),
        ("return_end", stage_records[return_stage][-1]),
    )
    for role, record in endpoint_roles:
        for component, field in component_fields.items():
            convergence_observables[f"{role}_{component}_flux_mol_m2_s"] = record[field]
    for role, (stage, _) in zip(
        ("baseline", "pulse", "return"),
        scenario.segments,
    ):
        for component in ("water", "hexane"):
            convergence_observables[f"{role}_{component}_flux_time_integral_mol_s_m2"] = (
                _stage_flux_time_integral(
                    step_records,
                    stage,
                    component_fields[component],
                )
            )
    job_id = (
        f"conditioned-cut-n{cells}-dt{dt_s:g}"
        if scenario.identifier == SCENARIO_COMBINED
        else f"conditioned-cut-{scenario.identifier}-n{cells}-dt{dt_s:g}"
    )
    suffix_wall_time = math.fsum(record["wall_time_s"] for record in suffix_records)
    if audit_final_state_sink is not None:
        if audit_final_state_sink:
            raise ValueError("audit final-state sink must be empty on entry")
        # Diagnostic plumbing only: exposing the already committed immutable
        # state changes no equation, seed, tolerance, or accepted trajectory.
        audit_final_state_sink.append(final)
    vanishing_wet_cut_treatment_history = [
        audit for record in step_records for audit in record["vanishing_wet_cut_treatment_history"]
    ]
    vanishing_wet_cut_treatment_route_history = [
        route
        for record in step_records
        for route in record["vanishing_wet_cut_treatment_route_history"]
    ]
    if len(vanishing_wet_cut_treatment_route_history) != len(vanishing_wet_cut_treatment_history):
        raise RuntimeError("wet-cut treatment route/audit histories lost alignment")
    job_maximum_wet_retained_water_loading = max(
        float(record["wet_retained_cap_maximum_loading"]) for record in step_records
    )
    job_guard_loadings = {record["wet_retained_water_guard_loading"] for record in step_records}
    job_capacity_loadings = {
        record["wet_retained_water_capacity_loading"] for record in step_records
    }
    if len(job_guard_loadings) != 1 or len(job_capacity_loadings) != 1:
        raise RuntimeError("accepted trajectory changed its retained-water guard/cap")
    job_guard_loading = next(iter(job_guard_loadings))
    job_capacity_loading = next(iter(job_capacity_loadings))
    return {
        "job_id": job_id,
        "completed": True,
        "physically_qualifying": False,
        "scenario": scenario.identifier,
        "forcing_class": scenario.forcing_class,
        "mechanism_gate_authority": scenario.gate_authority,
        "cells": cells,
        "radius_m": radius_m,
        "dt_s": dt_s,
        "segment_duration_s": SEGMENT_DURATION_S,
        "accepted_steps": len(step_records),
        "accepted_internal_substeps": sum(
            record.get("accepted_internal_substep_count", 1) for record in step_records
        ),
        "wall_time_s": math.fsum(record["wall_time_s"] for record in step_records),
        "incremental_campaign_wall_time_s": (
            suffix_wall_time
            if shared_prefix
            else math.fsum(record["wall_time_s"] for record in step_records)
        ),
        "shared_hot_lean_baseline_prefix_reused": shared_prefix,
        "post_disturbance_state_shared_between_scenarios": False,
        "total_nonlinear_evaluations": sum(
            record["nonlinear_evaluations"] for record in step_records
        ),
        "incremental_campaign_nonlinear_evaluations": (
            _attempt_evaluation_count(suffix_records)
            if shared_prefix
            else _attempt_evaluation_count(step_records)
        ),
        "scalar_post_root_polish_accounting": (scalar_post_root_polish_accounting),
        "incremental_campaign_scalar_post_root_polish_accounting": (
            incremental_scalar_post_root_polish_accounting
        ),
        "amendment_18_19_evidence_contract": amendment_18_19_evidence,
        "interface_saturation_endpoint_audit_accounting": (
            amendment_18_19_evidence[
                "interface_saturation_endpoint_audit_accounting"
            ]
        ),
        "dry_energy_capacity_stencil_audit_accounting": (
            amendment_18_19_evidence[
                "dry_energy_capacity_stencil_audit_accounting"
            ]
        ),
        "amendment_20_actual_composition_force_path_evidence_contract": (
            amendment_20_evidence
        ),
        "step_records": step_records,
        "vanishing_wet_cut_treatment_applicability": (
            "aggregate_of_applicable_strict_cut_records"
            if vanishing_wet_cut_treatment_history
            else "not_applicable_only_exact_zero_transactions"
        ),
        "vanishing_wet_cut_treatment": (
            None
            if not vanishing_wet_cut_treatment_history
            else vanishing_wet_cut_treatment_history[-1]
        ),
        "vanishing_wet_cut_treatment_history": (vanishing_wet_cut_treatment_history),
        "vanishing_wet_cut_treatment_route_history": (vanishing_wet_cut_treatment_route_history),
        "vanishing_wet_cut_activation_count": sum(
            int(audit["activation_count"]) for audit in vanishing_wet_cut_treatment_history
        ),
        "maximum_wet_retained_water_loading": (job_maximum_wet_retained_water_loading),
        "wet_retained_water_guard_loading": job_guard_loading,
        "minimum_wet_retained_water_guard_margin": (
            None
            if job_guard_loading is None
            else math.fsum((job_guard_loading, -job_maximum_wet_retained_water_loading))
        ),
        "wet_retained_water_capacity_loading": job_capacity_loading,
        "minimum_wet_retained_water_capacity_margin": math.fsum(
            (job_capacity_loading, -job_maximum_wet_retained_water_loading)
        ),
        "maximum_vanishing_wet_cut_water_pair_sum_difference_mol_s": max(
            (
                abs(float(audit["water_pair_sum_difference_mol_s"]))
                for audit in vanishing_wet_cut_treatment_history
            ),
            default=0.0,
        ),
        "maximum_vanishing_wet_cut_common_energy_pair_sum_difference_w": max(
            (
                abs(float(audit["common_energy_pair_sum_difference_w"]))
                for audit in vanishing_wet_cut_treatment_history
            ),
            default=0.0,
        ),
        "initial_material_volume_profile": (
            material_profiles.serialize_integrator_state(coarse_first_step.before)
        ),
        "component_surface_flux_diagnostics": component_flux_diagnostics,
        "hexane_finite_pulse_response": hexane_pulse_response,
        "independent_hexane_finite_pulse_response": (independent_hexane_pulse_response),
        "mechanism_specific_contract": mechanism_contract,
        "surface_constitutive_oracle_contract": surface_constitutive_contract,
        "instantaneous_flux_sign_diagnostics": {
            "excluded_from_mesh_time_acceptance": True,
            "independent_hexane_J_h": {
                "outward_to_inward": independent_inward_bracket,
                "inward_to_outward_return": independent_return_bracket,
            },
            "conserved_hexane_N_h": {
                "outward_to_inward": conserved_inward_bracket,
                "inward_to_outward_return": conserved_return_bracket,
            },
        },
        "pulse_layer_resolution_diagnostic": (_pulse_layer_resolution_diagnostic(cells, radius_m)),
        "retained_water_source_domain": {
            "direct_temperature_max_k": (RETAINED_WATER_DIRECT_TEMPERATURE_MAX_K),
            "any_temperature_extrapolation": any(
                record["retained_water_source_temperature_extrapolation"] for record in step_records
            ),
            "flag_is_explicit_on_every_step": True,
            "interpretation": (
                "use above 343.15 K is model-form extrapolation, not direct "
                "retained-water experimental validation"
            ),
        },
        "total_stefan_flux_diagnostic_non_gating": total_stefan_flux_diagnostic,
        "temperature_extrema_contract": {
            "scope": ACCEPTED_TEMPERATURE_EXTREMA_SCOPE,
            "includes_initial_state": True,
            "includes_every_accepted_BE_time_node": True,
            "includes_all_wet_dry_and_interface_material_pieces": True,
            "unrepresented_continuous_intra_step_extrema_claimed": False,
        },
        "convergence_observables": convergence_observables,
        "final": {
            "time_s": final.transport.time_s,
            "front_z": final.transport.geometry.front.z,
            "interface_temperature_k": final.last_interface_temperature_k,
            "water_inventory_mol": final_inventory.total_water_mol,
            "retained_water_inventory_mol": final_water_phases.retained_water_mol,
            "pore_vapor_water_inventory_mol": (final_water_phases.pore_vapor_water_mol),
            "free_liquid_water_inventory_mol": (final_water_phases.free_liquid_water_mol),
            "hexane_inventory_mol": final_inventory.total_hexane_mol,
            "energy_inventory_j": final_inventory.total_energy_j,
            "cumulative_boundary_water_out_mol": (
                final.corrected_cumulative_boundary_water_out_mol
            ),
            "cumulative_boundary_hexane_out_mol": (
                final.corrected_cumulative_boundary_hexane_out_mol
            ),
            "cumulative_boundary_energy_out_j": (final.corrected_cumulative_boundary_energy_out_j),
            "wet_temperatures_k": list(final.transport.wet_temperatures_k),
            "dry_temperatures_k": list(final.transport.dry_temperatures_k),
            "dry_y_hexane": list(final.transport.dry_y_hexane),
        },
        "gates": {
            "accepted_grid_contract": accepted_grid,
            "exact_requested_time_grid": accepted_grid["exact_requested_time_grid"],
            "frozen_failed_leaf_controller_contract_passed": accepted_grid[
                "frozen_failed_leaf_controller_contract_passed"
            ],
            "two_square_boundary_steps_on_requested_grid": (
                schedule_contract["true_finite_pulse"]
                and accepted_grid["square_transitions_align_with_grid"]
            ),
            "true_finite_pulse_boundary_return": schedule_contract["true_finite_pulse"],
            "strictly_inside_all_open_solver_bounds": strictly_inside_open_bounds,
            "no_clipping_or_bound_projection": strictly_inside_open_bounds,
            "no_tolerance_relaxation": controls_unchanged,
            "interface_saturation_endpoint_audit_contract_passed": (
                amendment_18_19_evidence[
                    "interface_saturation_endpoint_audit_accounting"
                ]["contract_passed"]
            ),
            "dry_energy_capacity_stencil_audit_contract_passed": (
                amendment_18_19_evidence[
                    "dry_energy_capacity_stencil_audit_accounting"
                ]["contract_passed"]
            ),
            "amendment_18_19_evidence_contract_passed": (
                amendment_18_19_evidence["contract_passed"]
            ),
            "amendment_20_actual_composition_force_path_contract_passed": (
                amendment_20_evidence["contract_passed"]
            ),
            "maximum_scaled_residual": maximum_residual,
            "maximum_step_ledger_residual": maximum_step,
            "maximum_cumulative_ledger_residual": maximum_cumulative,
            "maximum_normalized_water_residual": maximum_water_residual,
            "maximum_normalized_hexane_residual": maximum_hexane_residual,
            "maximum_normalized_energy_residual": maximum_energy_residual,
            "residual_contract_passed": maximum_residual <= NONLINEAR_RESIDUAL_LIMIT,
            "conservation_contract_passed": (
                max(maximum_step, maximum_cumulative) <= CONSERVATION_LIMIT
            ),
            # N_h remains the exact conserved component flux.  Its sign is not
            # the sign of the independent Maxwell--Stefan mechanism because
            # Stefan carriage can dominate it.  Only J_h supplies the declared
            # composition-response direction gate; N_h remains a mandatory
            # residual, ledger, and finite-window transfer diagnostic.
            "surface_hexane_component_flux_reversal_observed": (
                hexane_pulse_response["outward_to_inward_observed"]
            ),
            "surface_hexane_component_return_reversal_observed": (
                hexane_pulse_response["inward_to_outward_return_observed"]
            ),
            "surface_hexane_component_finite_pulse_response_observed": (
                hexane_pulse_response["outward_inward_outward_observed"]
            ),
            "surface_hexane_component_N_h_direction_is_gate": False,
            "surface_independent_hexane_counterdiffusion_pulse_observed": (
                independent_hexane_pulse_response["outward_inward_outward_observed"]
            ),
            "surface_independent_hexane_J_h_direction_is_gate": (
                scenario.identifier in (SCENARIO_COMPOSITION_ONLY, SCENARIO_COMBINED)
            ),
            "surface_total_stefan_flux_reversal_observed_diagnostic": (
                _outward_then_inward(
                    stage_fluxes(
                        baseline_stage,
                        "surface_total_stefan_flux_mol_m2_s",
                    ),
                    stage_fluxes(
                        pulse_stage,
                        "surface_total_stefan_flux_mol_m2_s",
                    ),
                )
            ),
            "surface_total_stefan_flux_is_direction_gate": False,
            "surface_binary_gas_entropy_contract_passed": (binary_gas_entropy_contract_passed),
            "surface_constitutive_oracle_contract_passed": (
                surface_constitutive_contract["contract_passed"]
            ),
            "mechanism_specific_contract_passed": mechanism_contract["mechanism_contract_passed"],
            "physically_qualifying": False,
        },
    }


def run_fixed_radius_trajectory(
    radius_m: float | None = None,
    *,
    cells: int = 2,
    time_step_s: float | None = None,
    scenario: StressScenario | None = None,
    seed_only_reference_trajectory: CoupledFixedRadiusTrajectory | None = None,
) -> CoupledFixedRadiusTrajectory:
    """Run one coupled three-segment pulse at an independently fixed radius.

    This is the radius-study backend, not a new solver.  It calls the same
    bootstrap, baseline, reversal, recovery, residual, and ledger paths used
    by :func:`run_campaign`.  The only varied physical input is the radius
    supplied to the master grid and front geometry.
    """

    if radius_m is None:
        radius_m = RADIUS_M
    if time_step_s is None:
        time_step_s = SEGMENT_DURATION_S
    if scenario is None:
        scenario = COMBINED_SCENARIO
    radius_m = _validate_radius(radius_m)
    if cells not in BOOTSTRAP_MESHES:
        raise ValueError("fixed-radius backend needs a declared bootstrap mesh")
    if time_step_s not in DEFAULT_TIMESTEPS_S:
        raise ValueError("fixed-radius backend needs a declared campaign timestep")
    if scenario.identifier not in SCENARIOS or SCENARIOS[scenario.identifier] != scenario:
        raise ValueError("fixed-radius backend needs an unchanged declared scenario")
    if seed_only_reference_trajectory is not None:
        reference = seed_only_reference_trajectory
        if time_step_s != SEGMENT_DURATION_S:
            raise ValueError(
                "cross-radius trajectory continuation requires the installed "
                "segment-duration timestep"
            )
        if (
            reference.cells != cells
            or reference.time_step_s != time_step_s
            or reference.scenario_id != scenario.identifier
            or reference.forcing_history_id != FIXED_RADIUS_FORCING_HISTORY_ID
        ):
            raise ValueError(
                "cross-radius seed reference must use the identical mesh, time "
                "grid, scenario, and forcing identity"
            )
        if len(reference.accepted_step_seed_candidates) != TRAJECTORY_SEGMENT_COUNT:
            raise ValueError("cross-radius seed reference is missing accepted roots")
        source_radius_m = reference.radius_m
        reference_candidates = reference.accepted_step_seed_candidates
    else:
        source_radius_m = None
        reference_candidates = None

    first_steps, bootstrap, first_step_audit = _bootstrap_first_steps(
        cells,
        radius_m,
        cross_radius_n2_seed_candidate=(
            None if reference_candidates is None else reference_candidates[0]
        ),
        cross_radius_seed_source_m=source_radius_m,
    )
    coarse_first_step = first_steps[cells]
    prefix = _prepare_baseline_prefix(
        cells,
        time_step_s,
        coarse_first_step,
        bootstrap_first_steps=first_steps,
    )
    accepted_steps: list[AcceptedSameCellStep] = []
    result = _run_trajectory(
        cells,
        time_step_s,
        coarse_first_step,
        scenario,
        prefix,
        cross_radius_suffix_seed_candidates=(
            None if reference_candidates is None else reference_candidates[1:]
        ),
        cross_radius_seed_source_m=source_radius_m,
        accepted_steps_out=accepted_steps,
    )
    if not all(record["particle_radius_m"] == radius_m for record in result["step_records"]):
        raise RuntimeError("accepted trajectory changed its fixed particle radius")

    initial = ci.inventory_snapshot(coarse_first_step.before.transport)
    final = result["final"]
    observables = result["convergence_observables"]
    gates = result["gates"]
    water_partition_total = math.fsum(
        (
            final["retained_water_inventory_mol"],
            final["pore_vapor_water_inventory_mol"],
            final["free_liquid_water_inventory_mol"],
        )
    )
    water_scale = max(abs(final["water_inventory_mol"]), 1.0e-300)
    if abs(water_partition_total - final["water_inventory_mol"]) / water_scale > 2.0e-14:
        raise RuntimeError("fixed-radius backend water partitions do not close")
    if not gates["conservation_contract_passed"]:
        raise RuntimeError("fixed-radius backend violated the frozen conservation limit")

    return CoupledFixedRadiusTrajectory(
        radius_m=radius_m,
        cells=cells,
        time_step_s=time_step_s,
        duration_s=TRAJECTORY_SEGMENT_COUNT * SEGMENT_DURATION_S,
        scenario_id=scenario.identifier,
        forcing_history_id=FIXED_RADIUS_FORCING_HISTORY_ID,
        initial_front_position_over_radius=INITIAL_FRONT_RADIUS_FRACTION,
        final_front_position_over_radius=observables["final_front_position_over_radius"],
        minimum_front_position_over_radius=observables["minimum_front_position_over_radius"],
        maximum_front_position_over_radius=observables["maximum_front_position_over_radius"],
        initial_water_inventory_mol=initial.total_water_mol,
        initial_hexane_inventory_mol=initial.total_hexane_mol,
        initial_energy_inventory_j=initial.total_energy_j,
        final_water_inventory_mol=final["water_inventory_mol"],
        final_retained_water_inventory_mol=final["retained_water_inventory_mol"],
        final_pore_vapor_water_inventory_mol=final["pore_vapor_water_inventory_mol"],
        final_free_liquid_water_inventory_mol=final["free_liquid_water_inventory_mol"],
        final_hexane_inventory_mol=final["hexane_inventory_mol"],
        final_energy_inventory_j=final["energy_inventory_j"],
        cumulative_boundary_water_out_mol=final["cumulative_boundary_water_out_mol"],
        cumulative_boundary_hexane_out_mol=final["cumulative_boundary_hexane_out_mol"],
        cumulative_boundary_energy_out_j=final["cumulative_boundary_energy_out_j"],
        minimum_material_temperature_k=observables["minimum_material_temperature_k"],
        peak_material_temperature_k=observables["peak_particle_temperature_k"],
        maximum_instantaneous_radial_temperature_spread_k=observables[
            "maximum_instantaneous_radial_temperature_spread_k"
        ],
        temperature_extrema_scope=ACCEPTED_TEMPERATURE_EXTREMA_SCOPE,
        unrepresented_continuous_intra_step_extrema_claimed=False,
        maximum_normalized_water_residual=gates["maximum_normalized_water_residual"],
        maximum_normalized_hexane_residual=gates["maximum_normalized_hexane_residual"],
        maximum_normalized_energy_residual=gates["maximum_normalized_energy_residual"],
        maximum_scaled_nonlinear_residual=gates["maximum_scaled_residual"],
        first_step_audit=first_step_audit,
        nonqualifying_coefficient_fixture=True,
        dirichlet_outer_trace_without_qualified_film=True,
        physically_qualifying=False,
        diagnostics={
            "bootstrap": bootstrap,
            "trajectory": result,
            "instantaneous_surface_flux_claim_boundary": {
                "ideal_Dirichlet_square_jump": True,
                "surface_flux_samples_are_diagnostic": True,
                "surface_flux_samples_are_mesh_acceptance_quantities": False,
                "finite_window_integrated_component_fluxes_are_convergence_quantities": (True),
                "component_inventories_are_convergence_quantities": True,
                "finite_external_film_installed": False,
            },
            "cross_radius_seed_only_contract": {
                "used": seed_only_reference_trajectory is not None,
                "source_radius_m": source_radius_m,
                "maximum_attempts_per_accepted_step": (CROSS_RADIUS_SEED_MAXIMUM_ATTEMPTS_PER_STEP),
                "target_initial_state_constructed_independently": True,
                "accepted_state_carried_between_radii": False,
                "conserved_inventory_interpolated_or_remapped": False,
                "equations_forcing_dt_and_tolerances_unchanged": True,
                "all_rejections_require_exact_rollback": True,
            },
        },
        accepted_step_seed_candidates=tuple(step.assembly.candidate for step in accepted_steps),
    )


def _outward_then_inward(
    first_stage_fluxes: Sequence[float],
    reversed_stage_fluxes: Sequence[float],
) -> bool:
    """Return the signed component-flux reversal contract.

    Flux is positive outward.  Both stages must be represented, and the
    declared reversal must contain an outward value before an inward value.
    No tolerance, clipping, or total-Stefan surrogate is used here.
    """

    if not first_stage_fluxes or not reversed_stage_fluxes:
        return False
    values = (*first_stage_fluxes, *reversed_stage_fluxes)
    if not all(math.isfinite(value) for value in values):
        return False
    return max(first_stage_fluxes) > 0.0 and min(reversed_stage_fluxes) < 0.0


def _maximum_vector_difference(
    left: cut.CutTransportUnknowns,
    right: cut.CutTransportUnknowns,
) -> float:
    return float(
        np.max(
            np.abs(
                np.asarray(left.lossless_state_vector()) - np.asarray(right.lossless_state_vector())
            )
        )
    )


def _determinism_audit(reference: ci.CutIntegratorStep) -> dict:
    seed = reference.seed
    started = perf_counter()
    serial = ci.advance_same_cell_backward_euler(
        reference.before,
        reference.ledger.dt_s,
        reference.boundary,
        seed,
    )
    serial_wall = perf_counter() - started
    restarted_before, restarted_boundary, restarted_seed = pickle.loads(
        pickle.dumps((reference.before, reference.boundary, seed))
    )
    started = perf_counter()
    restarted = ci.advance_same_cell_backward_euler(
        restarted_before,
        reference.ledger.dt_s,
        restarted_boundary,
        restarted_seed,
    )
    restart_wall = perf_counter() - started
    serial_difference = _maximum_vector_difference(
        reference.assembly.candidate,
        serial.assembly.candidate,
    )
    restart_difference = _maximum_vector_difference(
        reference.assembly.candidate,
        restarted.assembly.candidate,
    )
    replay_occurrences = (
        _amendment_18_19_execution_occurrence(
            serial,
            execution_classification=(
                "noncommitted_determinism_or_restart_replay"
            ),
            execution_context={
                "execution_partition_id": (
                    "determinism_restart_replay_executions"
                ),
                "execution_role": "determinism_serial_replay",
                "replay_execution_index": 0,
                "mesh_cells": reference.before.transport.geometry.master_grid.n,
                "dt_s": reference.ledger.dt_s,
                "pickle_restart_used": False,
            },
        ),
        _amendment_18_19_execution_occurrence(
            restarted,
            execution_classification=(
                "noncommitted_determinism_or_restart_replay"
            ),
            execution_context={
                "execution_partition_id": (
                    "determinism_restart_replay_executions"
                ),
                "execution_role": "determinism_restart_replay",
                "replay_execution_index": 1,
                "mesh_cells": reference.before.transport.geometry.master_grid.n,
                "dt_s": reference.ledger.dt_s,
                "pickle_restart_used": True,
            },
        ),
    )
    amendment_partition = _aggregate_amendment_18_19_diagnostic_partition(
        "determinism_restart_replay_executions",
        replay_occurrences,
    )
    return {
        "mesh_cells": reference.before.transport.geometry.master_grid.n,
        "dt_s": reference.ledger.dt_s,
        "serial_replay_wall_time_s": serial_wall,
        "restart_replay_wall_time_s": restart_wall,
        "serial_candidate_max_abs_difference": serial_difference,
        "restart_candidate_max_abs_difference": restart_difference,
        "serial_scaled_residuals_bit_exact": (
            serial.ledger.scaled_residuals == reference.ledger.scaled_residuals
        ),
        "restart_scaled_residuals_bit_exact": (
            restarted.ledger.scaled_residuals == reference.ledger.scaled_residuals
        ),
        "serial_evaluations_equal": (
            serial.ledger.nonlinear_evaluations == reference.ledger.nonlinear_evaluations
        ),
        "restart_evaluations_equal": (
            restarted.ledger.nonlinear_evaluations == reference.ledger.nonlinear_evaluations
        ),
        "scalar_post_root_polish_accounting": (
            _aggregate_scalar_post_root_polish_summaries(
                (
                    _cut_step_scalar_post_root_polish_summary(
                        serial,
                        "determinism_serial_replay",
                    ),
                    _cut_step_scalar_post_root_polish_summary(
                        restarted,
                        "determinism_restart_replay",
                    ),
                )
            )
        ),
        "amendment_18_19_diagnostic_execution_partition": amendment_partition,
        "bit_exact_passed": (
            serial_difference == 0.0
            and restart_difference == 0.0
            and serial.ledger.scaled_residuals == reference.ledger.scaled_residuals
            and restarted.ledger.scaled_residuals == reference.ledger.scaled_residuals
        ),
    }


INTEGRATED_CONVERGENCE_OBSERVABLES = (
    "final_front_z",
    "final_front_position_over_radius",
    "minimum_front_position_over_radius",
    "maximum_front_position_over_radius",
    "final_interface_temperature_k",
    "minimum_material_temperature_k",
    "peak_particle_temperature_k",
    "maximum_instantaneous_radial_temperature_spread_k",
    "final_retained_water_mol",
    "final_pore_vapor_water_mol",
    "final_free_liquid_water_mol",
    "final_total_water_mol",
    "final_total_hexane_mol",
    "final_total_energy_j",
    "cumulative_boundary_hexane_out_mol",
    "cumulative_boundary_energy_out_j",
    *(
        f"{role}_{component}_flux_time_integral_mol_s_m2"
        for role in ("baseline", "pulse", "return")
        for component in ("water", "hexane")
    ),
)

# A mathematical Dirichlet square jump has no resolved external-film time or
# length scale.  Right-endpoint surface-flux samples are retained for sign and
# mechanism diagnostics, but are not Gate-1g mesh-acceptance observables.
INSTANTANEOUS_JUMP_FLUX_DIAGNOSTICS = (
    *(
        f"{role}_{component}_flux_mol_m2_s"
        for role in (
            "baseline_end",
            "pulse_first",
            "pulse_end",
            "return_first",
            "return_end",
        )
        for component in ("water", "hexane", "hexane_independent_diffusion")
    ),
)
PRIMARY_OBSERVABLES = (
    *INTEGRATED_CONVERGENCE_OBSERVABLES,
    *INSTANTANEOUS_JUMP_FLUX_DIAGNOSTICS,
)


def _difference(left: dict, right: dict) -> dict:
    left_observables = left["convergence_observables"]
    right_observables = right["convergence_observables"]
    differences = {}
    for name in PRIMARY_OBSERVABLES:
        left_value = left_observables[name]
        right_value = right_observables[name]
        differences[name] = (
            abs(left_value - right_value)
            if isinstance(left_value, (int, float)) and isinstance(right_value, (int, float))
            else None
        )
    return differences


def _observed_order(
    coarse_delta: float | None,
    fine_delta: float | None,
) -> float | None:
    if coarse_delta is None or fine_delta is None or coarse_delta <= 0.0 or fine_delta <= 0.0:
        return None
    return math.log2(coarse_delta / fine_delta)


def _refinement_report(results: Sequence[dict]) -> dict:
    completed = [item for item in results if item.get("completed", False)]
    by_key = {(item["cells"], item["dt_s"]): item for item in completed}
    meshes = sorted({item["cells"] for item in completed})
    timesteps = sorted({item["dt_s"] for item in completed}, reverse=True)
    mesh_differences = []
    for dt_s in timesteps:
        available_meshes = [cells for cells in meshes if (cells, dt_s) in by_key]
        sequence = [by_key[(cells, dt_s)] for cells in available_meshes]
        differences = [_difference(left, right) for left, right in zip(sequence, sequence[1:])]
        orders = []
        for left, right in zip(differences, differences[1:]):
            orders.append(
                {name: _observed_order(left[name], right[name]) for name in PRIMARY_OBSERVABLES}
            )
        mesh_differences.append(
            {
                "dt_s": dt_s,
                "adjacent_pairs": [
                    {
                        "coarse_cells": coarse,
                        "fine_cells": fine,
                        "absolute_differences": difference,
                    }
                    for coarse, fine, difference in zip(
                        available_meshes,
                        available_meshes[1:],
                        differences,
                    )
                ],
                "triplet_observed_orders": orders,
            }
        )
    time_differences = []
    for cells in meshes:
        available_timesteps = [dt_s for dt_s in timesteps if (cells, dt_s) in by_key]
        sequence = [by_key[(cells, dt_s)] for dt_s in available_timesteps]
        differences = [_difference(left, right) for left, right in zip(sequence, sequence[1:])]
        orders = []
        for left, right in zip(differences, differences[1:]):
            orders.append(
                {name: _observed_order(left[name], right[name]) for name in PRIMARY_OBSERVABLES}
            )
        time_differences.append(
            {
                "cells": cells,
                "adjacent_pairs": [
                    {
                        "coarse_dt_s": coarse,
                        "fine_dt_s": fine,
                        "absolute_differences": difference,
                    }
                    for coarse, fine, difference in zip(
                        available_timesteps,
                        available_timesteps[1:],
                        differences,
                    )
                ],
                "triplet_observed_orders": orders,
            }
        )
    return {
        "mesh": mesh_differences,
        "time": time_differences,
        "mesh_acceptance_observables": list(INTEGRATED_CONVERGENCE_OBSERVABLES),
        "diagnostic_only_observables": list(INSTANTANEOUS_JUMP_FLUX_DIAGNOSTICS),
        "instantaneous_jump_flux_mesh_acceptance_required": False,
        "finite_external_film_installed": False,
        "excluded_failed_jobs": [
            item["job_id"] for item in results if not item.get("completed", False)
        ],
        "interpretation": (
            "absolute adjacent differences and diagnostic observed orders only; "
            "finite-window flux integrals, inventories, front, and temperature "
            "metrics are the convergence quantities. Instantaneous surface-flux "
            "samples adjacent to the ideal Dirichlet jumps are sign/mechanism "
            "diagnostics until a finite external film is installed. No physical "
            "relevance tolerance is inferred from this numerical oracle"
        ),
    }


def _pressure_conditioned_contracts() -> list[dict]:
    records = []
    for pressure_pa in (101_325.0, 135_000.0, 150_000.0, 170_000.0):
        config = replace(
            DRY_CONFIG,
            pressure_pa=pressure_pa,
            primitive_band=ct.GasOnlyPrimitiveDomain(
                (330.0, 408.15),
                "multi-pressure topology contract; no dynamic solve",
            ),
        )
        diagnostic = config.conditioned_temperature_domain
        lower, upper = diagnostic.solver_bounds_k
        temperature = lower + 0.35 * (upper - lower)
        composition = cp.decode_gas_only_y(
            temperature,
            pressure_pa,
            0.0,
            config.pore,
        )
        state = cp.evaluate_equilibrium(
            temperature,
            pressure_pa,
            composition,
            config.pore,
        )
        interval = cp.gas_only_composition_interval(
            temperature,
            pressure_pa,
            config.pore,
        )
        records.append(
            {
                "pressure_pa": pressure_pa,
                "requested_temperature_bounds_k": list(diagnostic.requested_bounds_k),
                "solver_temperature_bounds_k": list(diagnostic.solver_bounds_k),
                "lower_is_phase_boundary": diagnostic.lower_is_phase_boundary,
                "sample_temperature_k": temperature,
                "sample_y_hexane": composition,
                "sample_interval_y_hexane": [
                    interval.lower_y_hexane,
                    interval.upper_y_hexane,
                ],
                "sample_hexane_activity": state.hexane_activity,
                "sample_water_activity": state.water_activity,
                "strictly_phase_feasible": (
                    interval.lower_y_hexane < composition < interval.upper_y_hexane
                    and state.hexane_activity < 1.0
                    and state.water_activity < 1.0
                ),
                "dynamic_cut_solve_was_run": pressure_pa == PRESSURE_PA,
            }
        )
    return records


def _external_binary_gas_phase_context(
    temperature_k: float,
    boundary: ct.PoreBoundary,
) -> dict:
    water_activity, hexane_activity = cp.binary_gas_component_activities(
        temperature_k,
        PRESSURE_PA,
        boundary.y_hexane,
        DRY_CONFIG.pore,
    )
    mole_fraction_interior = 0.0 < boundary.y_hexane < 1.0
    water_below_unity = water_activity < 1.0
    hexane_below_unity = hexane_activity < 1.0
    return {
        "temperature_k": temperature_k,
        "y_hexane": boundary.y_hexane,
        "water_activity": water_activity,
        "hexane_activity": hexane_activity,
        "mole_fraction_strictly_interior": mole_fraction_interior,
        "water_activity_below_unity": water_below_unity,
        "hexane_activity_below_unity": hexane_below_unity,
        "strictly_phase_feasible": (
            mole_fraction_interior and water_below_unity and hexane_below_unity
        ),
        "retained_sorption_authority_used": False,
    }


def _boundary_face_thermodynamic_context(
    face_temperature_k: float,
    boundary: ct.PoreBoundary,
) -> dict:
    phase = _external_binary_gas_phase_context(face_temperature_k, boundary)
    interval = cp.gas_only_composition_interval(
        face_temperature_k,
        PRESSURE_PA,
        DRY_CONFIG.pore,
    )
    chart_compatible = interval.lower_y_hexane < boundary.y_hexane < interval.upper_y_hexane
    dry_chart = {
        "open_interval_y_hexane": [
            interval.lower_y_hexane,
            interval.upper_y_hexane,
        ],
        "lower_constraint": interval.lower_constraint.value,
        "upper_constraint": interval.upper_constraint.value,
        "composition_inside_open_chart": chart_compatible,
        "includes_modified_luikov_lower_storage_authority": True,
        "is_external_bulk_gas_phase_envelope": False,
    }
    accepted = phase["strictly_phase_feasible"] and chart_compatible
    return {
        "temperature_k": face_temperature_k,
        "external_binary_gas_phase": phase,
        "particle_dry_storage_chart": dry_chart,
        "accepted_numerical_boundary_context": accepted,
        # Backward-compatible fields.  They are now explicitly labelled as a
        # combined phase-plus-storage-chart result, not a bulk phase boundary.
        "gas_only_interval_y_hexane": dry_chart["open_interval_y_hexane"],
        "hexane_activity": phase["hexane_activity"],
        "water_activity": phase["water_activity"],
        "strictly_interior": accepted,
        "strictly_interior_semantics": ("combined_external_phase_and_particle_dry_storage_chart"),
    }


def _boundary_phase_safety(boundary: ct.PoreBoundary) -> dict:
    domain = DRY_CONFIG.conditioned_temperature_domain
    face_temperatures = {
        "coldest_possible_symmetric_face": 0.5
        * (domain.solver_bounds_k[0] + boundary.temperature_k),
        "hottest_possible_symmetric_face": 0.5
        * (domain.solver_bounds_k[1] + boundary.temperature_k),
    }
    contexts = {
        name: _boundary_face_thermodynamic_context(face_temperature, boundary)
        for name, face_temperature in face_temperatures.items()
    }
    external_bulk = _external_binary_gas_phase_context(
        boundary.temperature_k,
        boundary,
    )
    all_face_phases_feasible = all(
        context["external_binary_gas_phase"]["strictly_phase_feasible"]
        for context in contexts.values()
    )
    all_face_charts_compatible = all(
        context["particle_dry_storage_chart"]["composition_inside_open_chart"]
        for context in contexts.values()
    )
    external_phase_feasible = external_bulk["strictly_phase_feasible"] and all_face_phases_feasible
    accepted = external_phase_feasible and all_face_charts_compatible
    coldest = contexts["coldest_possible_symmetric_face"]
    return {
        "boundary": asdict(boundary),
        "external_bulk_binary_gas_phase": external_bulk,
        "conservative_symmetric_face_contexts": contexts,
        "external_binary_gas_phase_feasible": external_phase_feasible,
        "particle_dry_storage_chart_compatible": all_face_charts_compatible,
        "accepted_numerical_boundary_context": accepted,
        # Backward-compatible cold-face fields.
        "coldest_possible_symmetric_face_temperature_k": coldest["temperature_k"],
        "gas_only_interval_y_hexane": coldest["gas_only_interval_y_hexane"],
        "hexane_activity": coldest["hexane_activity"],
        "water_activity": coldest["water_activity"],
        "strictly_interior": accepted,
        "strictly_interior_semantics": ("combined_external_phase_and_particle_dry_storage_chart"),
        "physically_qualifying": False,
        "physical_feasibility_envelope_identified": False,
        "treated_as_source_extrapolation_for_physical_claims": True,
        "retained_water_direct_source_temperature_max_k": (RETAINED_WATER_DIRECT_TEMPERATURE_MAX_K),
        "boundary_temperature_above_retained_water_direct_source_range": (
            boundary.temperature_k > RETAINED_WATER_DIRECT_TEMPERATURE_MAX_K
        ),
        "interpretation": (
            "the external bulk and both limiting symmetric-face gas states are "
            "phase checked independently of retained sorption; particle-face "
            "compatibility with the separate dry retained-storage chart is also "
            "required. The endpoint remains a nonqualifying numerical source "
            "extrapolation because no soybean-DT feasibility envelope is identified"
        ),
    }


def _dynamic_schedule_contract(
    scenario: StressScenario | None = None,
) -> dict:
    if scenario is None:
        scenario = COMBINED_SCENARIO
    segments = [
        {
            "stage": stage,
            "start_s": index * SEGMENT_DURATION_S,
            "end_s": (index + 1) * SEGMENT_DURATION_S,
            "boundary": asdict(boundary),
        }
        for index, (stage, boundary) in enumerate(scenario.segments)
    ]
    baseline_return_is_exact = scenario.segments[0][1] is (scenario.segments[2][1])
    pulse_differs_from_baseline = scenario.segments[1][1] != (scenario.segments[0][1])
    return {
        "scenario": scenario.identifier,
        "forcing_class": scenario.forcing_class,
        "mechanism_gate_authority": scenario.gate_authority,
        "segments": segments,
        "square_step_transition_times_s": [
            SEGMENT_DURATION_S,
            2.0 * SEGMENT_DURATION_S,
        ],
        "accepted_boundary_interpolation": False,
        "finite_pulse_duration_s": SEGMENT_DURATION_S,
        "pulse_differs_from_baseline": pulse_differs_from_baseline,
        "return_uses_exact_same_baseline_object": baseline_return_is_exact,
        "true_finite_pulse": (pulse_differs_from_baseline and baseline_return_is_exact),
    }


def _numerical_case_authority_metadata(
    selected_scenarios: Sequence[StressScenario] | None = None,
) -> dict:
    """Expose the configured case without promoting it to a physical band."""

    binary = DRY_CONFIG.binary_diffusivity
    retained = DRY_CONFIG.retained_water_mobility
    audited_binary_bracket = (2.0e-10, 6.0e-10)
    scenarios = (
        tuple(SCENARIOS.values()) if selected_scenarios is None else tuple(selected_scenarios)
    )
    boundaries = tuple(boundary for scenario in scenarios for _, boundary in scenario.segments)
    all_reduced_film = bool(boundaries) and all(
        isinstance(boundary, ct.ReducedFilmPoreBoundary) for boundary in boundaries
    )
    all_dirichlet = bool(boundaries) and all(
        isinstance(boundary, ct.DirichletPoreBoundary) for boundary in boundaries
    )
    if all_reduced_film:
        outer_boundary_type = "ReducedFilmPoreBoundary_fast_mass_finite_heat"
    elif all_dirichlet:
        outer_boundary_type = "DirichletPoreBoundary"
    else:
        outer_boundary_type = "mixed_explicit_pore_boundary_types"
    return {
        "local_evidence_audit": {
            "record_id": AUDIT_RECORD_ID,
            "path": AUDIT_RECORD_PATH,
            "status": "evidence_audit_only_no_authority_or_pass_flag_change",
        },
        "physical_feasibility_envelope": None,
        "physical_coefficient_corner_set": None,
        "configured_single_numerical_case": {
            "pore_effective_binary_diffusivity": {
                "interval_m2_s": [
                    binary.interval.lower_m2_s,
                    binary.interval.upper_m2_s,
                ],
                "selected_fraction": binary.fraction,
                "selected_value_m2_s": binary.value_m2_s,
                "provenance": binary.interval.provenance,
                "physically_qualifying": binary.physically_qualifying,
            },
            "retained_matrix_water_mobility": {
                "interval_mol_m_s": [
                    retained.interval.lower_mol_m_s,
                    retained.interval.upper_mol_m_s,
                ],
                "selected_fraction": retained.fraction,
                "selected_value_mol_m_s": retained.value_mol_m_s,
                "provenance": retained.interval.label,
                "physically_qualifying": retained.physically_qualifying,
            },
            "dry_thermal_conductivity": {
                "selected_value_w_m_k": (DRY_CONFIG.thermal_conductivity.value_w_m_k),
                "provenance": DRY_CONFIG.thermal_conductivity.label,
                "physically_qualifying": (DRY_CONFIG.thermal_conductivity.physically_qualifying),
            },
            # The P1EF reduced boundary keeps a finite exact-A heat film while
            # taking only the gas-side mass resistance to its audited fast-film
            # limit.  It is neither a Dirichlet trace nor the general finite
            # binary mass/heat/contact model that remains outside P1EF.
            "finite_external_binary_film_installed": False,
            "reduced_fast_mass_finite_heat_film_installed": all_reduced_film,
            "outer_boundary_type": outer_boundary_type,
            "physically_qualifying": False,
        },
        "audited_nonqualifying_brackets_not_claimed_as_physical": {
            "pore_effective_binary_diffusivity_m2_s": {
                "interval": list(audited_binary_bracket),
                "audit_location": "section 6, table row 1",
                "authorized_use": (
                    "manufactured numerical interval centered on the restricted "
                    "Cardarelli scalar-limit oracle"
                ),
                "configured_case_is_inside_interval": (
                    audited_binary_bracket[0] <= binary.value_m2_s <= audited_binary_bracket[1]
                ),
                "physically_qualifying": False,
            },
            "retained_matrix_water_mobility_mol_m_s": None,
            "dry_thermal_conductivity_w_m_k": None,
        },
        "coefficient_corner_coverage_performed": False,
        "interpretation": (
            "The current root-preserving stress fixture is one manufactured case. "
            "Its selected binary diffusivity lies outside the audit's documented "
            "nonqualifying numerical bracket, and no audited numerical interval "
            "maps to the retained-matrix mobility or dry conductivity used here. "
            "Therefore this run is not coefficient-corner evidence."
        ),
        "physically_qualifying": False,
    }


def _resolve_scenarios(requested: Sequence[str]) -> tuple[StressScenario, ...]:
    names = tuple(requested)
    if not names:
        raise ValueError("campaign needs at least one stress scenario")
    if "all" in names:
        if len(names) != 1:
            raise ValueError("'all' cannot be combined with named scenarios")
        names = tuple(SCENARIOS)
    unknown = set(names) - set(SCENARIOS)
    if unknown:
        raise ValueError(f"unknown Gate-1g stress scenarios: {sorted(unknown)}")
    selected = set(names)
    return tuple(scenario for name, scenario in SCENARIOS.items() if name in selected)


def _scenario_job_id(scenario: StressScenario, cells: int, dt_s: float) -> str:
    if scenario.identifier == SCENARIO_COMBINED:
        return f"conditioned-cut-n{cells}-dt{dt_s:g}"
    return f"conditioned-cut-{scenario.identifier}-n{cells}-dt{dt_s:g}"


def _strict_json_failure_evidence_copy(
    value,
    *,
    location: str,
    nonfinite_paths: list[str],
):
    """Copy failure evidence while losslessly tagging nonfinite diagnostics."""

    if is_dataclass(value):
        return _strict_json_failure_evidence_copy(
            asdict(value),
            location=location,
            nonfinite_paths=nonfinite_paths,
        )
    if isinstance(value, Enum):
        return _strict_json_failure_evidence_copy(
            value.value,
            location=location,
            nonfinite_paths=nonfinite_paths,
        )
    if isinstance(value, (float, np.floating)):
        scalar = float(value)
        if not math.isfinite(scalar):
            if math.isnan(scalar):
                label = "NaN"
            elif scalar > 0.0:
                label = "+Infinity"
            else:
                label = "-Infinity"
            nonfinite_paths.append(location)
            return {"__nonfinite_float__": label}
        return scalar
    if isinstance(value, Mapping):
        return {
            str(key): _strict_json_failure_evidence_copy(
                item,
                location=f"{location}.{key}",
                nonfinite_paths=nonfinite_paths,
            )
            for key, item in value.items()
        }
    if isinstance(value, (tuple, list, np.ndarray)):
        return [
            _strict_json_failure_evidence_copy(
                item,
                location=f"{location}[{index}]",
                nonfinite_paths=nonfinite_paths,
            )
            for index, item in enumerate(value)
        ]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return repr(value)


def _strict_json_failure_evidence(payload: Mapping[str, object]) -> dict[str, object]:
    """Return strict JSON without making a nonfinite diagnostic look finite."""

    nonfinite_paths: list[str] = []
    encoded = _strict_json_failure_evidence_copy(
        payload,
        location="$",
        nonfinite_paths=nonfinite_paths,
    )
    if not isinstance(encoded, dict):
        raise TypeError("strict-JSON failure evidence root must be an object")
    encoded["strict_json_nonfinite_encoding"] = {
        "schema_version": 1,
        "output_contains_only_standard_json_numbers": True,
        "source_nonfinite_float_count": len(nonfinite_paths),
        "source_nonfinite_float_paths": nonfinite_paths,
        "tag_format": {"__nonfinite_float__": "+Infinity|-Infinity|NaN"},
        "finite_values_changed": False,
        "acceptance_fields_recomputed_or_relaxed": False,
        "tagged_values_are_nonnumeric_and_cannot_pass_finite_checks": True,
        "interpretation": (
            "nonfinite diagnostic sentinels retain their exact class and sign; "
            "no solver result or qualification decision is changed"
        ),
    }
    return encoded


def _requested_macro_input_evidence(state: object) -> dict[str, object]:
    """Serialize the immutable adaptive input without retaining a state object."""

    if not isinstance(state, ci.CutIntegratorState):
        return {
            "serialization_available": False,
            "python_type": f"{type(state).__module__}.{type(state).__qualname__}",
            "reason": "requested input is not a CutIntegratorState",
        }
    try:
        return _requested_macro_input_evidence_unchecked(state)
    except Exception as error:  # noqa: BLE001 - failure reporting must not mask solver error
        transport = getattr(state, "transport", None)
        geometry = getattr(transport, "geometry", None)
        front_state = getattr(geometry, "front", None)
        grid = getattr(geometry, "master_grid", None)
        return {
            "serialization_available": False,
            "python_type": f"{type(state).__module__}.{type(state).__qualname__}",
            "time_s": getattr(transport, "time_s", None),
            "front_z": getattr(front_state, "z", None),
            "mesh_cells": getattr(grid, "n", None),
            "serialization_error_type": type(error).__name__,
            "serialization_error": str(error),
            "reason": "requested input evidence serialization failed closed",
            "contains_committable_state_object": False,
            "accepted_state_changed_by_serialization": False,
            "physically_qualifying": False,
        }


def _requested_macro_input_evidence_unchecked(
    state: ci.CutIntegratorState,
) -> dict[str, object]:
    """Build the complete scalar/array view for a typed immutable input."""

    transport = state.transport
    return {
        "serialization_available": True,
        "material_volume_profile": material_profiles.serialize_integrator_state(state),
        "differential_state": {
            "time_s": transport.time_s,
            "front_z": transport.geometry.front.z,
            "mesh_cells": transport.geometry.master_grid.n,
            "cut_cell_index": transport.geometry.cut_cell_index,
            "wet_temperatures_k": transport.wet_temperatures_k,
            "wet_retained_water_loadings": transport.wet_retained_water_loadings,
            "wet_retained_water_capacity_duals_over_rt": (
                transport.effective_wet_retained_water_capacity_duals_over_rt
            ),
            "dry_temperatures_k": transport.dry_temperatures_k,
            "dry_y_hexane": transport.dry_y_hexane,
            "historical_hexane_loadings": transport.historical_hexane_loadings,
            "oil_fraction_labels": transport.oil_fraction_labels,
            "last_interface_temperature_k": state.last_interface_temperature_k,
            "last_total_stefan_fluxes_mol_m2_s": (state.last_total_stefan_fluxes_mol_m2_s),
        },
        "direct_material_inventory": ci.inventory_snapshot(transport),
        "direct_water_phase_inventory": ci.water_phase_inventory_snapshot(transport),
        "integrator_reference_and_compensated_history": {
            "reference_inventory": state.reference_inventory,
            "reference_capacity_energy_scale_j": state.reference_capacity_energy_scale_j,
            "cumulative_boundary_water_out_mol": (state.cumulative_boundary_water_out_mol),
            "cumulative_boundary_hexane_out_mol": (state.cumulative_boundary_hexane_out_mol),
            "cumulative_boundary_energy_out_j": state.cumulative_boundary_energy_out_j,
            "boundary_water_compensation_mol": state.boundary_water_compensation_mol,
            "boundary_hexane_compensation_mol": state.boundary_hexane_compensation_mol,
            "boundary_energy_compensation_j": state.boundary_energy_compensation_j,
            "corrected_cumulative_boundary_water_out_mol": (
                state.corrected_cumulative_boundary_water_out_mol
            ),
            "corrected_cumulative_boundary_hexane_out_mol": (
                state.corrected_cumulative_boundary_hexane_out_mol
            ),
            "corrected_cumulative_boundary_energy_out_j": (
                state.corrected_cumulative_boundary_energy_out_j
            ),
            "cumulative_absolute_water_transfer_mol": (
                state.cumulative_absolute_water_transfer_mol
            ),
            "cumulative_absolute_hexane_transfer_mol": (
                state.cumulative_absolute_hexane_transfer_mol
            ),
            "cumulative_absolute_energy_transfer_j": (state.cumulative_absolute_energy_transfer_j),
            "cumulative_material_water_change_cell_mol": (
                state.cumulative_material_water_change_cell_mol
            ),
            "cumulative_material_hexane_change_cell_mol": (
                state.cumulative_material_hexane_change_cell_mol
            ),
            "cumulative_material_energy_change_cell_j": (
                state.cumulative_material_energy_change_cell_j
            ),
            "water_change_compensation_cell_mol": (state.water_change_compensation_cell_mol),
            "hexane_change_compensation_cell_mol": (state.hexane_change_compensation_cell_mol),
            "energy_change_compensation_cell_j": state.energy_change_compensation_cell_j,
            "accepted_steps": state.accepted_steps,
            "cumulative_nonlinear_evaluations": state.cumulative_nonlinear_evaluations,
        },
        "contains_committable_state_object": False,
        "accepted_state_changed_by_serialization": False,
        "physically_qualifying": False,
    }


def _adaptive_failure_in_exception_chain(
    error: BaseException,
) -> AdaptiveSameCellMacrostepError | None:
    """Find adaptive evidence even when a shared-prefix wrapper is outermost."""

    current: BaseException | None = error
    visited: set[int] = set()
    while current is not None and id(current) not in visited:
        visited.add(id(current))
        if isinstance(current, AdaptiveSameCellMacrostepError):
            return current
        current = current.__cause__ if current.__cause__ is not None else current.__context__
    return None


def _last_step_error_evidence(
    error: ci.CutIntegratorStepError | None,
    *,
    requested_input_state: object,
) -> dict[str, object] | None:
    """Copy the terminal rejected candidate diagnostics without its rollback state."""

    if error is None:
        return None
    return {
        "error_type": type(error).__name__,
        "error": str(error),
        "rollback_identity_is_requested_macro_input": (
            error.rollback_state is requested_input_state
        ),
        "nonlinear_evaluations": error.nonlinear_evaluations,
        "rejected_trial_evaluations": error.rejected_trial_evaluations,
        "maximum_scaled_residual": error.maximum_scaled_residual,
        "condition_proxy": error.condition_proxy,
        "event_restart_required": error.event_restart_required,
        "last_candidate": error.last_candidate,
        "scaled_residuals": error.scaled_residuals,
        "residual_blocks": error.residual_blocks,
        "minimum_fractional_distance_to_bound": (error.minimum_fractional_distance_to_bound),
        "nonlinear_optimizer_function_evaluations": (
            error.nonlinear_optimizer_function_evaluations
        ),
        "provisional_work_budget": error.provisional_work_budget,
        "provisional_work_budget_exhausted": error.provisional_work_budget_exhausted,
        "retained_water_applicability_guard_failed": (
            error.retained_water_applicability_guard_failed
        ),
        "retained_water_applicability_audit": error.retained_water_applicability_audit,
    }


def _adaptive_same_cell_failure_evidence(
    error: AdaptiveSameCellMacrostepError,
) -> dict[str, object]:
    """Serialize every immutable Amendment-09 failure field at the catch boundary."""

    requested_transport = getattr(error.requested_input_state, "transport", None)
    raw: dict[str, object] = {
        "schema_version": 1,
        "serialization_completed": True,
        "error_type": type(error).__name__,
        "error": str(error),
        "classification": error.classification,
        "requested_macro_input_time_s": getattr(requested_transport, "time_s", None),
        "requested_macro_input": _requested_macro_input_evidence(error.requested_input_state),
        "requested_dt_s": error.requested_dt_s,
        "minimum_leaf_dt_s": error.minimum_leaf_dt_s,
        "attempted_binary_depths": error.attempted_binary_depths,
        "attempts": error.attempts,
        "last_step_error": _last_step_error_evidence(
            error.last_step_error,
            requested_input_state=error.requested_input_state,
        ),
        "initial_failure_audit": error.initial_failure_audit,
        "last_auditable_failure_jacobian_audit": (error.last_auditable_failure_jacobian_audit),
        "provisional_accepted_leaf_count": len(error.provisional_accepted_leaf_summaries),
        "provisional_accepted_leaf_summaries": (error.provisional_accepted_leaf_summaries),
        "predictor_backtrack_maximum_binary_exponent": (
            error.predictor_backtrack_maximum_binary_exponent
        ),
        "rollback_identity_preserved": error.rollback_identity_preserved,
        "every_rejected_attempt_rollback_identity_preserved": all(
            attempt.accepted or attempt.rollback_identity_preserved for attempt in error.attempts
        ),
        "failed_state_committed": error.failed_state_committed,
        "trajectory_state_committed_from_provisional_leaves": False,
        "physically_qualifying": False,
    }
    return _strict_json_failure_evidence(raw)


def _adaptive_failure_serialization_fallback(
    error: AdaptiveSameCellMacrostepError,
    serialization_error: Exception,
) -> dict[str, object]:
    """Preserve fail-closed identity if optional rich reporting itself fails."""

    requested_transport = getattr(error.requested_input_state, "transport", None)
    return _strict_json_failure_evidence(
        {
            "schema_version": 1,
            "serialization_completed": False,
            "serialization_error_type": type(serialization_error).__name__,
            "serialization_error": str(serialization_error),
            "rich_fields_omitted_due_to_serialization_failure": True,
            "error_type": type(error).__name__,
            "error": str(error),
            "classification": error.classification,
            "requested_macro_input_time_s": getattr(
                requested_transport,
                "time_s",
                None,
            ),
            "requested_dt_s": error.requested_dt_s,
            "minimum_leaf_dt_s": error.minimum_leaf_dt_s,
            "attempted_binary_depths": error.attempted_binary_depths,
            "rollback_identity_preserved": error.rollback_identity_preserved,
            "failed_state_committed": error.failed_state_committed,
            "trajectory_state_committed_from_provisional_leaves": False,
            "physically_qualifying": False,
        }
    )


def _campaign_failure_result(
    scenario: StressScenario,
    cells: int,
    dt_s: float,
    error: Exception,
) -> dict[str, object]:
    """Build one fail-closed job record, retaining adaptive failure evidence."""

    result: dict[str, object] = {
        "job_id": _scenario_job_id(scenario, cells, dt_s),
        "completed": False,
        "scenario": scenario.identifier,
        "forcing_class": scenario.forcing_class,
        "cells": cells,
        "dt_s": dt_s,
        "error_type": type(error).__name__,
        "error": str(error),
        "failed_state_was_committed": False,
        "physically_qualifying": False,
    }
    adaptive_error = _adaptive_failure_in_exception_chain(error)
    if adaptive_error is not None:
        result["failed_state_was_committed"] = adaptive_error.failed_state_committed
        try:
            adaptive_evidence = _adaptive_same_cell_failure_evidence(adaptive_error)
        except Exception as serialization_error:  # noqa: BLE001 - preserve solver failure
            adaptive_evidence = _adaptive_failure_serialization_fallback(
                adaptive_error,
                serialization_error,
            )
        result["adaptive_same_cell_failure_evidence"] = adaptive_evidence
    return result


# --- GT-PS-2-P1E-R7: pinned wall evidence -----------------------------------
#
# The four Phase-1 walls are deterministic: for one frozen source state and
# one complete physical job identity, the failing job record is bit-stable.
# Re-solving a proven wall costs ~5-10 minutes of subdivision per campaign
# while producing an already-committed outcome.  R7 pins the wall job record
# once per (source, physical-identity) key and replays the pin as evidence,
# exactly as every other pinned artifact in this project: the pin stores the
# FULL original job record, the complete key, and the documented failure
# signature; any change to sources or job identity invalidates the key and
# forces a genuine re-solve.  Pins are evidence caching, not gate relaxation:
# the R3/R4 verified fail-closed contracts evaluate the replayed record
# identically.  The registry is evidence-based: only walls with committed
# mechanism records (R2/R3/R4) are pinnable, only under the scenario they
# were measured in.

# GT-PS-2-P1E-R11: the A15 full direction-reversal demonstration gates only
# at nominal film stress (lambda=1); stress-varied corners gate on the
# measured universal invariant (pulse strictly decreases J_h from a strictly
# outward baseline).  Installed per case by the qualifier; the selector's
# offline recompute installs the same value per invocation, keeping the
# producer and validator sides symmetric.
A15_FULL_REVERSAL_GATING = True

WALL_PIN_DIRECTORY = Path("numerical_convergence_results/p1ef_wall_pins")
WALL_PIN_CONTEXT: dict | None = None  # set by qualifier.run_case per case
R7_WALL_SIGNATURES = {
    (SCENARIO_COMPOSITION_ONLY, 48, 0.01875): (
        "subdivision exhausted its documented minimum dt"
    ),
    (SCENARIO_COMPOSITION_ONLY, 96, 0.0375): (
        "subdivision exhausted its documented minimum dt"
    ),
    (SCENARIO_COMPOSITION_ONLY, 96, 0.01875): (
        "cannot be enclosed inside the caller-oriented phase domain"
    ),
    # R14 (2026-08-05): the combined (cool/rich) scenario at the pinned
    # joint forcing hits the SAME three walls at the SAME cells with the
    # SAME measured signatures (pocket at 48/0.01875 and 96/0.0375; water
    # decidability barrier at 96/0.01875) — measured in the first gate1h
    # campaign at the selected forcing.  Same mechanisms, same fail-closed
    # accounting, scenario-keyed registration.
    (SCENARIO_COMBINED, 48, 0.01875): (
        "subdivision exhausted its documented minimum dt"
    ),
    (SCENARIO_COMBINED, 96, 0.0375): (
        "subdivision exhausted its documented minimum dt"
    ),
    (SCENARIO_COMBINED, 96, 0.01875): (
        "cannot be enclosed inside the caller-oriented phase domain"
    ),
}
R7_CORNER_WALL_SIGNATURES = {
    ("db_min_dapp_max", 0.5, SCENARIO_COMPOSITION_ONLY, 12, 0.075): (
        "subdivision exhausted its documented minimum dt"
    ),
    # R10: composition-band corner wall, measured for w95p1 only.
    ("db_max_dapp_max", 0.5, SCENARIO_COMPOSITION_ONLY, 12, 0.075): (
        "subdivision exhausted its documented minimum dt"
    ),
}


def _wall_pin_key(scenario_id: str, cells: int, dt_s: float) -> dict | None:
    """Complete physical-identity key, or None when pinning is inapplicable."""

    context = WALL_PIN_CONTEXT
    if not isinstance(context, dict):
        return None
    signature = R7_WALL_SIGNATURES.get((scenario_id, cells, dt_s))
    if signature is None:
        corner = (
            context.get("coefficient_case_id"),
            context.get("stress_multiplier"),
            scenario_id,
            cells,
            dt_s,
        )
        signature = R7_CORNER_WALL_SIGNATURES.get(corner)
    if signature is None:
        return None
    key = {
        "record": "GT-PS-2-P1E-R7",
        "source_combined_sha256": context.get("source_combined_sha256"),
        "coefficient_case_id": context.get("coefficient_case_id"),
        "stress_multiplier": context.get("stress_multiplier"),
        "structural_conductivity_w_m_k": context.get(
            "structural_conductivity_w_m_k"
        ),
        "pulse_y_hexane": context.get("pulse_y_hexane"),
        "scenario": scenario_id,
        "cells": cells,
        "dt_s": dt_s,
        "documented_failure_signature": signature,
    }
    # B-F1 (2026-09-29): a face-tangent knob that is not at its on-disk
    # default is part of the physical identity; a knob-on wall is never
    # replayed into a knob-off run or the reverse.  Absent at the defaults, so
    # a knob-off key keeps its pre-B-F1 form.
    for knob in ("face_tangent_inward_seed_armed", "face_tangent_residual_bound"):
        if knob in context:
            key[knob] = context[knob]
    # The face-route time-refined prefix knob likewise, only when armed.
    if face_route_time_refined_prefix_armed():
        key["face_route_time_refined_prefix_armed"] = True
    return key


def _wall_pin_path(key: dict) -> Path:
    digest = hashlib.sha256(
        json.dumps(key, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:24]
    return WALL_PIN_DIRECTORY / f"wall_pin_{digest}.json"


def _load_wall_pin(scenario_id: str, cells: int, dt_s: float) -> dict | None:
    key = _wall_pin_key(scenario_id, cells, dt_s)
    if key is None:
        return None
    path = _wall_pin_path(key)
    if not path.is_file():
        return None
    try:
        stored = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if stored.get("key") != key:
        return None
    record = stored.get("job_record")
    if not isinstance(record, dict) or record.get("completed") is not False:
        return None
    if key["documented_failure_signature"] not in str(record.get("error") or ""):
        return None
    replay = dict(record)
    replay["replayed_from_wall_pin"] = True
    replay["wall_pin_key"] = key
    replay["wall_pin_path"] = str(path)
    return replay


def _store_wall_pin(
    scenario_id: str, cells: int, dt_s: float, record: dict
) -> None:
    key = _wall_pin_key(scenario_id, cells, dt_s)
    if key is None or record.get("completed") is not False:
        return
    if key["documented_failure_signature"] not in str(record.get("error") or ""):
        return
    WALL_PIN_DIRECTORY.mkdir(parents=True, exist_ok=True)
    payload = {"key": key, "job_record": record}
    _wall_pin_path(key).write_text(
        json.dumps(payload, indent=1, sort_keys=True, allow_nan=True),
        encoding="utf-8",
    )


def run_campaign(
    meshes: Sequence[int],
    timesteps_s: Sequence[float],
    scenarios: Sequence[str] = (SCENARIO_COMBINED,),
) -> dict:
    meshes = tuple(sorted(set(meshes)))
    timesteps = tuple(sorted(set(timesteps_s), reverse=True))
    selected_scenarios = _resolve_scenarios(scenarios)
    if not meshes or not timesteps:
        raise ValueError("campaign needs at least one mesh and timestep")
    if not set(meshes) <= set(DEFAULT_MESHES):
        raise ValueError("meshes must be selected from 12, 24, 48, and 96")
    if not set(timesteps) <= set(DEFAULT_TIMESTEPS_S):
        choices = ", ".join(f"{value:g}" for value in DEFAULT_TIMESTEPS_S)
        raise ValueError(f"timesteps must be selected from {choices} s")
    # The face-route time-refined prefix knob: a malformed value refuses here.
    face_route_prefix_armed = face_route_time_refined_prefix_armed()
    campaign_started = perf_counter()
    first_steps, bootstrap, _ = _bootstrap_first_steps(max(meshes))
    results = []
    baseline_prefixes: dict[str, dict] = {}
    for cells in meshes:
        for dt_s in timesteps:
            if selected_scenarios and all(
                _load_wall_pin(scenario.identifier, cells, dt_s) is not None
                for scenario in selected_scenarios
            ):
                for scenario in selected_scenarios:
                    results.append(
                        _load_wall_pin(scenario.identifier, cells, dt_s)
                    )
                continue
            try:
                prefix = _prepare_baseline_prefix(
                    cells,
                    dt_s,
                    first_steps[cells],
                    bootstrap_first_steps=first_steps,
                )
                baseline_prefixes[f"n{cells}_dt{dt_s:g}"] = {
                    "cells": cells,
                    "dt_s": dt_s,
                    "accepted_steps": len(prefix.step_records),
                    "last_time_s": prefix.last_step.after.transport.time_s,
                    "last_front_z": prefix.last_step.after.transport.geometry.front.z,
                    "incremental_wall_time_s": prefix.incremental_wall_time_s,
                    "incremental_nonlinear_evaluations": (prefix.incremental_nonlinear_evaluations),
                    "scalar_post_root_polish_accounting": (
                        _aggregate_scalar_post_root_polish_summaries(
                            tuple(
                                record["scalar_post_root_polish"]
                                for record in prefix.step_records
                                if any(
                                    attempt["role"] != "mesh_bootstrap_exact_root"
                                    for attempt in record["attempts"]
                                )
                            )
                        )
                    ),
                    "consumer_scenarios": [scenario.identifier for scenario in selected_scenarios],
                    "same_immutable_solver_state_reused": True,
                    "post_disturbance_state_shared": False,
                    "work_counted_once_in_campaign_accounting": True,
                }
                if prefix.face_route_time_refined_seed_route is not None:
                    baseline_prefixes[f"n{cells}_dt{dt_s:g}"][
                        "face_route_time_refined_seed_route"
                    ] = prefix.face_route_time_refined_seed_route
                prefix_error = None
            except Exception as exc:  # noqa: BLE001 - preserve prefix failure
                prefix = None
                prefix_error = exc
            for scenario in selected_scenarios:
                pinned = _load_wall_pin(scenario.identifier, cells, dt_s)
                if pinned is not None:
                    results.append(pinned)
                    continue
                try:
                    if prefix_error is not None:
                        raise RuntimeError(
                            "shared hot/lean baseline prefix failed without state "
                            f"commit: {prefix_error}"
                        ) from prefix_error
                    assert prefix is not None
                    result = _run_trajectory(
                        cells,
                        dt_s,
                        first_steps[cells],
                        scenario,
                        prefix,
                    )
                except Exception as exc:  # noqa: BLE001 - preserve campaign failures
                    result = _campaign_failure_result(scenario, cells, dt_s, exc)
                    _store_wall_pin(scenario.identifier, cells, dt_s, result)
                results.append(result)
                latest = results[-1]
                progress = {
                    "job_id": latest["job_id"],
                    "scenario": scenario.identifier,
                    "completed": latest["completed"],
                }
                if latest["completed"]:
                    progress.update(
                        wall_time_s=latest["wall_time_s"],
                        evaluations=latest["total_nonlinear_evaluations"],
                        gates=latest["gates"],
                    )
                else:
                    progress.update(
                        error_type=latest["error_type"],
                        error=latest["error"],
                    )
                print(json.dumps(progress), flush=True)
    determinism = _determinism_audit(first_steps[min(meshes)])
    contracts = _pressure_conditioned_contracts()
    schedule_contracts = {
        scenario.identifier: _dynamic_schedule_contract(scenario) for scenario in selected_scenarios
    }
    boundary_phase_safety_by_scenario = {
        scenario.identifier: {
            stage: _boundary_phase_safety(boundary) for stage, boundary in scenario.segments
        }
        for scenario in selected_scenarios
    }
    diagnostic_partitions = {
        "bootstrap_radius_seed_executions": bootstrap[2][
            "amendment_18_19_diagnostic_execution_partition"
        ],
        "determinism_restart_replay_executions": determinism[
            "amendment_18_19_diagnostic_execution_partition"
        ],
    }
    diagnostic_evidence_contract = (
        _aggregate_amendment_18_19_diagnostic_partitions(
            diagnostic_partitions,
            expected_partition_ids=tuple(diagnostic_partitions),
        )
    )
    amendment_18_19_structural_gate = _amendment_18_19_raw_structural_gate(
        results,
        diagnostic_partitions,
        diagnostic_evidence_contract,
    )
    # R7: verified wall jobs are fail-closed by design and are exempt from
    # the completed/gates aggregation; their verification lives in the
    # selection contract's excluded-cell reports.
    structural_results = [
        item for item in results if not _is_r7_wall_job(item)
    ]
    structural_pass = (
        all(item.get("completed", False) for item in structural_results)
        and all(
            item["gates"]["residual_contract_passed"]
            and item["gates"]["conservation_contract_passed"]
            and item["gates"]["exact_requested_time_grid"]
            and item["gates"]["frozen_failed_leaf_controller_contract_passed"]
            and item["gates"]["no_clipping_or_bound_projection"]
            and item["gates"]["no_tolerance_relaxation"]
            and item["gates"]["interface_saturation_endpoint_audit_contract_passed"]
            and item["gates"]["dry_energy_capacity_stencil_audit_contract_passed"]
            and item["gates"]["amendment_18_19_evidence_contract_passed"]
            and item["gates"][
                "amendment_20_actual_composition_force_path_contract_passed"
            ]
            and item["gates"]["two_square_boundary_steps_on_requested_grid"]
            and item["gates"]["true_finite_pulse_boundary_return"]
            and item["gates"]["mechanism_specific_contract_passed"]
            and item["gates"]["surface_constitutive_oracle_contract_passed"]
            and item["gates"]["surface_binary_gas_entropy_contract_passed"]
            and not item["gates"]["surface_total_stefan_flux_is_direction_gate"]
            for item in structural_results
        )
        and determinism["bit_exact_passed"]
        and amendment_18_19_structural_gate["contract_passed"]
        and all(item["true_finite_pulse"] for item in schedule_contracts.values())
        and all(
            endpoint["external_binary_gas_phase_feasible"]
            and endpoint["particle_dry_storage_chart_compatible"]
            and endpoint["accepted_numerical_boundary_context"]
            and not endpoint["physically_qualifying"]
            and endpoint["treated_as_source_extrapolation_for_physical_claims"]
            for scenario_safety in boundary_phase_safety_by_scenario.values()
            for endpoint in scenario_safety.values()
        )
    )
    domain = DRY_CONFIG.conditioned_temperature_domain
    reversal_context = _boundary_face_thermodynamic_context(
        COLDEST_REVERSAL_FACE_TEMPERATURE_K,
        COOL_RICH_BOUNDARY,
    )
    requested_scenario_names = [scenario.identifier for scenario in selected_scenarios]
    omitted_scenario_names = [name for name in SCENARIOS if name not in requested_scenario_names]
    full_matrix_requested = not omitted_scenario_names
    refinements = {
        scenario.identifier: _refinement_report(
            [item for item in results if item.get("scenario") == scenario.identifier]
        )
        for scenario in selected_scenarios
    }
    bootstrap_wall = math.fsum(item["wall_time_s"] for item in bootstrap.values())
    bootstrap_evaluations = sum(item["nonlinear_evaluations"] for item in bootstrap.values())
    prefix_wall = math.fsum(item["incremental_wall_time_s"] for item in baseline_prefixes.values())
    prefix_evaluations = sum(
        item["incremental_nonlinear_evaluations"] for item in baseline_prefixes.values()
    )
    suffix_wall = math.fsum(item.get("incremental_campaign_wall_time_s", 0.0) for item in results)
    suffix_evaluations = sum(
        item.get("incremental_campaign_nonlinear_evaluations", 0) for item in results
    )
    determinism_wall = math.fsum(
        (
            determinism["serial_replay_wall_time_s"],
            determinism["restart_replay_wall_time_s"],
        )
    )
    determinism_evaluations = 2 * first_steps[min(meshes)].ledger.nonlinear_evaluations
    bootstrap_scalar_post_root_accounting = _aggregate_scalar_post_root_polish_summaries(
        tuple(item["scalar_post_root_polish"] for item in bootstrap.values())
    )
    scalar_post_root_polish_accounting = _combine_scalar_post_root_polish_accounting(
        (
            bootstrap_scalar_post_root_accounting,
            *(item["scalar_post_root_polish_accounting"] for item in baseline_prefixes.values()),
            *(
                item["incremental_campaign_scalar_post_root_polish_accounting"]
                for item in results
                if item.get("completed") is True
            ),
            determinism["scalar_post_root_polish_accounting"],
        )
    )
    combined_safety = boundary_phase_safety_by_scenario.get(
        SCENARIO_COMBINED,
        {},
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "scope": "Gate_1g_conditioned_same_cut_cell_dynamic_stress",
        "generated_nonqualifying_numerical_evidence": True,
        "physically_qualifying": False,
        "physical_qualification_blockers": [
            "soybean water/hexane pore, matrix, and film coefficient intervals are not identified",
            (
                "the bounded reduced fast-mass/exact-A film is an engineering "
                "reduction; a target-identified general finite film/contact law is "
                "not installed"
                if all(
                    isinstance(boundary, ct.ReducedFilmPoreBoundary)
                    for scenario in selected_scenarios
                    for _, boundary in scenario.segments
                )
                else (
                    "outer trace is Dirichlet; the complete finite film/contact law "
                    "is not installed"
                )
            ),
        ],
        "gate_1g_coverage_not_provided_by_this_harness": [
            (
                "approved industry feasibility-envelope corners: no such envelope "
                "is identified in the local evidence audit"
            ),
            (
                "coefficient corners: this is one manufactured fixture and its "
                "binary diffusivity is outside the audit's documented numerical bracket"
            ),
            "dynamic topology-event ordering and event-time convergence",
            "parallel determinism (this harness checks serial replay and restart only)",
            "predictive validation against held-out soybean DT pulse data",
        ],
        "equations_changed_by_harness": False,
        "accepted_state_clipping_or_projection": False,
        "pressure_pa": PRESSURE_PA,
        "radius_m": RADIUS_M,
        "initial_front_z": INITIAL_FRONT_Z,
        "stress_matrix_scenario_selection": {
            "requested": requested_scenario_names,
            "omitted": omitted_scenario_names,
            "full_matrix_requested": full_matrix_requested,
            "full_matrix_completed": (
                full_matrix_requested and all(item.get("completed", False) for item in results)
            ),
            "section_8_3_numerical_matrix_coverage_passed": (
                full_matrix_requested and structural_pass
            ),
            "combined_only_is_not_full_matrix_coverage": True,
            "physically_qualifying": False,
        },
        "dynamic_schedules_by_scenario": schedule_contracts,
        "boundary_phase_safety_by_scenario": (boundary_phase_safety_by_scenario),
        # Backward-compatible combined-only views are present only when that
        # scenario was actually requested.
        "dynamic_schedule": (
            schedule_contracts[SCENARIO_COMBINED]["segments"]
            if SCENARIO_COMBINED in schedule_contracts
            else None
        ),
        "dynamic_schedule_contract": schedule_contracts.get(SCENARIO_COMBINED),
        "boundary_phase_safety": combined_safety or None,
        "reversal_boundary_phase_safety": {
            "coldest_possible_symmetric_face_temperature_k": (COLDEST_REVERSAL_FACE_TEMPERATURE_K),
            "exact_chart_coordinate": REVERSAL_EXACT_CHART_COORDINATE,
            "decoded_y_hexane": COOL_RICH_BOUNDARY.y_hexane,
            **reversal_context,
            "interpretation": (
                "the outer trace is binary-gas phase feasible at the coldest "
                "possible common face and separately lies inside the particle dry "
                "retained-storage chart; neither result identifies the cause of a "
                "failed nonlinear trial"
            ),
            "executed_in_this_campaign": SCENARIO_COMBINED in schedule_contracts,
            "physically_qualifying": False,
            "treated_as_source_extrapolation_for_physical_claims": True,
        },
        "primitive_chart": {
            "requested_temperature_bounds_k": list(domain.requested_bounds_k),
            "exact_pressure_conditioned_solver_bounds_k": list(domain.solver_bounds_k),
            "lower_is_phase_boundary": domain.lower_is_phase_boundary,
            "old_artificial_campaign_floor_k": 339.15,
            "band_correction_classification": (
                "nonqualifying numerical chart correction; exact phase inequalities retained"
            ),
        },
        "solver_acceptance_contract": {
            "backward_euler_grids_s": list(timesteps),
            "accepted_subdivision_policy": ADAPTIVE_FAILED_LEAF_POLICY_ID,
            "arbitrary_accepted_subdivision_allowed": False,
            "face_conditioned_dry_composition_predictor": (
                _face_conditioned_predictor_contract()
            ),
            "interface_saturation_endpoint_selection": (
                _interface_root_contract()
            ),
            "dry_energy_capacity_normalization_stencil": (
                _dry_energy_capacity_stencil_contract()
            ),
            "failed_leaf_controller": {
                "policy_id": ADAPTIVE_FAILED_LEAF_POLICY_ID,
                "licensed_only_for_rejected_same_cell_first_material_step_after_imposed_boundary_discontinuity": True,
                "absolute_minimum_leaf_dt_s": (ADAPTIVE_MACROSTEP_ABSOLUTE_MINIMUM_LEAF_DT_S),
                "maximum_binary_depth_cap": ADAPTIVE_MACROSTEP_MAXIMUM_BINARY_DEPTH,
                "formal_grid_effective_depths": [
                    _adaptive_floor_contract(value) for value in timesteps
                ],
                "failed_leaf_only_left_to_right_binary_recursion": True,
                "provisional_optimizer_function_evaluations_per_ordered_seed": (
                    DISCONTINUITY_LEAF_MAXIMUM_OPTIMIZER_FUNCTION_EVALUATIONS
                ),
                "reversal_maximum_nonlinear_seed_attempts": (
                    DISCONTINUITY_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
                ),
                "nonreversal_maximum_nonlinear_seed_attempts": (
                    DISCONTINUITY_NONREVERSAL_LEAF_MAXIMUM_NONLINEAR_SEED_ATTEMPTS
                ),
                "predictor_backtrack_maximum_binary_exponent": (
                    ADAPTIVE_PREDICTOR_BACKTRACK_MAX_BINARY_EXPONENT
                ),
                "requested_macro_nodes_remain_exact": True,
                "accepted_internal_leaf_sizes_and_nodes_disclosed": True,
                "atomic_macrostep_commit": True,
                "exact_rollback_on_rejection": True,
                "raw_macro_component_energy_ledger_limit": CONSERVATION_LIMIT,
                "wall_clock_timeout_controls_acceptance": False,
                "per_seed_and_per_leaf_work_budget_unchanged_from_amendment03": True,
                "total_tree_work_envelope_changed_by_amendment09": True,
            },
            "accepted_state_clipping_or_projection_allowed": False,
            "residual_tolerance": CONTROLS.nonlinear_residual_tolerance,
            "ledger_tolerance": CONTROLS.ledger_tolerance,
            "nonlinear_step_tolerance": CONTROLS.nonlinear_step_tolerance,
            "maximum_function_evaluations": (CONTROLS.maximum_function_evaluations),
            "maximum_condition_proxy": CONTROLS.maximum_condition_proxy,
            "tolerance_relaxation_allowed": False,
            "seed_only_boundary_homotopy_may_change_accepted_equations": False,
            "scalar_post_root_polish": {
                "numerical_amendment_id": (SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID),
                "near_miss_multiplier_ceiling": (ci.SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER),
                "residual_tolerance_unchanged": True,
                "interface_tolerance_unchanged": True,
                "physics_equations_changed": False,
                "active_or_contact_retained_cap_branch_eligible": False,
                "provisional_work_budget_step_eligible": False,
                "failed_polish_commits_partial_candidate": False,
                "exact_rollback_on_failure": True,
                "every_use_must_be_serialized": True,
                "absence_means_zero_uses": True,
            },
        },
        "numerical_case_authority": _numerical_case_authority_metadata(selected_scenarios),
        "bootstrap": bootstrap,
        "shared_hot_lean_baseline_prefixes": baseline_prefixes,
        **(
            {
                "face_route_time_refined_prefix": {
                    "environment_variable": FACE_ROUTE_TIME_REFINED_PREFIX_ENVIRONMENT_VARIABLE,
                    "armed": True,
                    "routes_taken": sorted(
                        key
                        for key, item in baseline_prefixes.items()
                        if "face_route_time_refined_seed_route" in item
                    ),
                    "tolerances_budgets_gates_changed": False,
                }
            }
            if face_route_prefix_armed
            else {}
        ),
        "jobs": results,
        "refinement": refinements.get(SCENARIO_COMBINED),
        "refinement_by_scenario": refinements,
        "determinism": determinism,
        "amendment_18_19_diagnostic_execution_partitions": (
            diagnostic_partitions
        ),
        "amendment_18_19_diagnostic_execution_contract": (
            diagnostic_evidence_contract
        ),
        "amendment_18_19_structural_gate": amendment_18_19_structural_gate,
        "pressure_conditioned_topology_contracts": contracts,
        "campaign_work_accounting": {
            "bootstrap_wall_time_s": bootstrap_wall,
            "shared_baseline_prefix_wall_time_s": prefix_wall,
            "scenario_suffix_wall_time_s": suffix_wall,
            "determinism_audit_wall_time_s": determinism_wall,
            "bootstrap_nonlinear_evaluations": bootstrap_evaluations,
            "shared_baseline_prefix_nonlinear_evaluations": prefix_evaluations,
            "scenario_suffix_nonlinear_evaluations": suffix_evaluations,
            "determinism_audit_nonlinear_evaluations": determinism_evaluations,
            "total_nonlinear_evaluations_counting_shared_work_once": (
                bootstrap_evaluations
                + prefix_evaluations
                + suffix_evaluations
                + determinism_evaluations
            ),
            "shared_work_is_not_multiplied_by_scenario_count": True,
            "wall_components_exclude_python_orchestration_overhead": True,
        },
        "scalar_post_root_polish_accounting": (scalar_post_root_polish_accounting),
        "campaign_wall_time_s": perf_counter() - campaign_started,
        "all_structural_numerical_contracts_passed": structural_pass,
        "approved_feasibility_envelope_coverage_passed": False,
        "coefficient_corner_coverage_passed": False,
        "phase_1_gate_1g_physical_pass": False,
    }


def _parse_csv(values: str, cast) -> tuple:
    return tuple(cast(value.strip()) for value in values.split(",") if value.strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshes", default="12,24,48,96")
    parser.add_argument("--timesteps", default="0.1,0.05,0.025")
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=(*SCENARIOS, "all"),
        default=(SCENARIO_COMBINED,),
        help=(
            "stress branches to execute; use 'all' for the full separated "
            "composition/temperature/combined matrix (default: combined only)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("numerical_convergence_results/gate1g_conditioned_cut_stress.json"),
    )
    arguments = parser.parse_args()
    payload = run_campaign(
        _parse_csv(arguments.meshes, int),
        _parse_csv(arguments.timesteps, float),
        arguments.scenarios,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(arguments.output),
                "job_count": len(payload["jobs"]),
                "campaign_wall_time_s": payload["campaign_wall_time_s"],
                "all_structural_numerical_contracts_passed": payload[
                    "all_structural_numerical_contracts_passed"
                ],
                "phase_1_gate_1g_physical_pass": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
