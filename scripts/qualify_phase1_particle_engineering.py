"""Run one bounded P1EF dynamic-stress qualification case.

The mature Gate-1g numerical orchestration is intentionally reused rather
than copied.  This script installs one owner-authorized engineering
configuration into that backend inside a dedicated process, obtains an exact
N=2 root for the unchanged target equations, and then runs the requested
mesh/time/scenario matrix.  Module-global installation makes this runner
process-isolated and not thread-safe; independent coefficient/film cases may
be parallelized only as separate processes.

The resulting artifact is engineering evidence only.  It never changes the
normative F2/F3 status and never calls itself physically qualified or plant
predictive.
"""

from __future__ import annotations

import argparse
import csv
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from functools import wraps
import hashlib
import importlib.metadata
import inspect
import json
import math
import operator
import os
import platform
from pathlib import Path
from typing import Any, Mapping, ParamSpec, Sequence, TypeVar

import numpy as np

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import cut_continuation as cc
from dtdc_simulator.core2.particle import cut_face_tangent as face_tangent
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import engineering_foundation as p1ef
from dtdc_simulator.core2.particle import external_film_reduced as efr
from scripts import p1ef_failed_leaf_contract as failed_leaf_contract
from scripts import p1ef_initial_snapshot_provenance as initial_snapshot
from scripts import p1ef_surface_constitutive_oracle as surface_oracle
from scripts import qualify_gate1g_conditioned_cut_stress as backend


SCHEMA_VERSION = 1
ENGINEERING_SCOPE = "P1EF_fixed_radius_dynamic_stress_case"
BASELINE_TEMPERATURE_K = 408.15
PULSE_TEMPERATURE_K = 363.15
BASELINE_Y_HEXANE = 0.79
# GT-PS-2-P1E-R16: the A15 full-reversal demonstration binds at this
# nominal structural conductivity; coefficient corners gate on the
# measured invariant (see install_engineering_case).
NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K = 0.24
PULSE_Y_HEXANE = 0.85
INTERFACE_TEMPERATURE_BOUNDS_K = (334.5, 341.0)
WET_TEMPERATURE_BOUNDS_K = (323.15, 360.0)
REFERENCE_N2_SEED_TIMESTEP_S = 0.1
FORCING_TIMESCALE_S = 0.075
FORMAL_TIMESTEPS_S = (0.075, 0.0375, 0.01875)
FORMAL_HISTORY_ID = "p1ef_symmetric_0p075_hot_lean__cool_rich__exact_hot_lean_return_v1"
FORMAL_HISTORY_SELECTION_STATUS = "prospective_qualification_pending"
INITIAL_DRY_TEMPERATURE_K = 343.0
FAST_MASS_DEPARTURE_LIMIT = 0.01
NORMALIZED_SURFACE_HEAT_CLOSURE_LIMIT = 2.0e-12
P1EF_MAXIMUM_CONDITION_PROXY = 1.0e12
LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT = (
    efr.LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
)

# `install_engineering_case` intentionally installs a case process-locally for
# the standalone diagnostics.  Public `run_case` calls are instead
# transactional: every backend binding below is restored to the exact object
# that the caller had installed before entry.  Keep this inventory explicit so
# a newly introduced backend mutation cannot silently escape the transaction.
_BACKEND_CASE_GLOBAL_NAMES = (
    "PRESSURE_PA",
    "RADIUS_M",
    "DRY_CONFIG",
    "WET_MODEL",
    "CONFIG",
    "RETAINED_MOBILITY",
    "SEGMENT_DURATION_S",
    "DEFAULT_TIMESTEPS_S",
    "CONTROLS",
    "HOT_LEAN_BOUNDARY",
    "HOT_RICH_BOUNDARY",
    "COOL_LEAN_BOUNDARY",
    "COOL_RICH_BOUNDARY",
    "REVERSAL_TEMPERATURE_K",
    "COLDEST_REVERSAL_FACE_TEMPERATURE_K",
    "COLDEST_COMPOSITION_FACE_TEMPERATURE_K",
    "REVERSAL_Y_HEXANE",
    "COMPOSITION_REVERSAL_Y_HEXANE",
    "REVERSAL_EXACT_CHART_COORDINATE",
    "SCENARIOS",
    "COMBINED_SCENARIO",
    "DYNAMIC_BOUNDARY_SEGMENTS",
    "FIXED_RADIUS_FORCING_HISTORY_ID",
    "A15_FULL_REVERSAL_GATING",
    "N2_ROOT",
)
_P = ParamSpec("_P")
_R = TypeVar("_R")


def _snapshot_backend_case_state() -> dict[str, object]:
    """Capture exact backend bindings; installed objects are not copied."""

    return {name: getattr(backend, name) for name in _BACKEND_CASE_GLOBAL_NAMES}


def _restore_backend_case_state(snapshot: Mapping[str, object]) -> None:
    """Restore a previously captured backend installation by identity."""

    for name in reversed(_BACKEND_CASE_GLOBAL_NAMES):
        setattr(backend, name, snapshot[name])


@contextmanager
def _preserve_backend_case_state() -> Iterator[None]:
    """Restore the complete caller installation on every context exit."""

    snapshot = _snapshot_backend_case_state()
    try:
        yield
    finally:
        _restore_backend_case_state(snapshot)


def _rollback_backend_case_state_on_failure(
    function: Callable[_P, _R],
) -> Callable[_P, _R]:
    """Keep successful direct installs persistent but make failure atomic."""

    @wraps(function)
    def wrapped(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        snapshot = _snapshot_backend_case_state()
        committed = False
        try:
            result = function(*args, **kwargs)
            committed = True
            return result
        finally:
            if not committed:
                _restore_backend_case_state(snapshot)

    return wrapped


def _restore_backend_case_state_after_call(
    function: Callable[_P, _R],
) -> Callable[_P, _R]:
    """Make a public runner transactional without changing its signature."""

    @wraps(function)
    def wrapped(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        with _preserve_backend_case_state():
            return function(*args, **kwargs)

    return wrapped
SOURCE_PROVENANCE_PATHS = (
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_DECISION.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_02.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_03.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_04.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_05.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_07.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_08.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_09.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_10.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_11.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_12.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_13.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_14.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_15.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_16.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_17.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_18.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_19.md",
    "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_20.md",
    "scripts/qualify_phase1_particle_engineering.py",
    "scripts/select_p1ef_pulse_composition.py",
    "scripts/qualify_gate1g_conditioned_cut_stress.py",
    "scripts/p1ef_failed_leaf_contract.py",
    "scripts/p1ef_initial_snapshot_provenance.py",
    "scripts/p1ef_material_volume_profiles.py",
    "scripts/p1ef_surface_constitutive_oracle.py",
    "src/dtdc_simulator/core2/particle/actual_composition_force_path.py",
    "src/dtdc_simulator/core2/particle/bed_film.py",
    "src/dtdc_simulator/core2/particle/coefficient_manifold.py",
    "src/dtdc_simulator/core2/particle/conditioning.py",
    "src/dtdc_simulator/core2/particle/coupled_pore.py",
    "src/dtdc_simulator/core2/particle/coupled_transport.py",
    "src/dtdc_simulator/core2/particle/cut_birth_event.py",
    "src/dtdc_simulator/core2/particle/cut_continuation.py",
    "src/dtdc_simulator/core2/particle/cut_event_orchestrator.py",
    "src/dtdc_simulator/core2/particle/cut_face_departure.py",
    "src/dtdc_simulator/core2/particle/cut_face_departure_integrator.py",
    "src/dtdc_simulator/core2/particle/cut_face_event.py",
    "src/dtdc_simulator/core2/particle/cut_face_event_integrator.py",
    "src/dtdc_simulator/core2/particle/cut_face_tangent.py",
    "src/dtdc_simulator/core2/particle/cut_face_to_face.py",
    "src/dtdc_simulator/core2/particle/cut_face_to_face_integrator.py",
    "src/dtdc_simulator/core2/particle/cut_geometry.py",
    "src/dtdc_simulator/core2/particle/cut_integrator.py",
    "src/dtdc_simulator/core2/particle/cut_sparsity.py",
    "src/dtdc_simulator/core2/particle/cut_transport.py",
    "src/dtdc_simulator/core2/particle/dry_thermo.py",
    "src/dtdc_simulator/core2/particle/engineering_foundation.py",
    "src/dtdc_simulator/core2/particle/external_film_reduced.py",
    "src/dtdc_simulator/core2/particle/front.py",
    "src/dtdc_simulator/core2/particle/grid.py",
    "src/dtdc_simulator/core2/particle/nonisothermal_potential.py",
    "src/dtdc_simulator/core2/particle/surface_active_set.py",
    "src/dtdc_simulator/core2/particle/transport_coefficients.py",
    "src/dtdc_simulator/core2/particle/wet_core.py",
    "src/dtdc_simulator/core2/particle/wet_retained_cap.py",
    "src/dtdc_simulator/core2/particle/wet_water.py",
    "src/dtdc_simulator/core/critical_volume.py",
    "src/dtdc_simulator/core2/exact_cache.py",
    "src/dtdc_simulator/core2/props/binary_gas.py",
    "src/dtdc_simulator/core2/props/critical_volume.py",
    "src/dtdc_simulator/core2/props/hexane.py",
    "src/dtdc_simulator/core2/props/sorption.py",
    "src/dtdc_simulator/core2/props/state.py",
    "src/dtdc_simulator/core2/props/water.py",
)
AMENDMENT_06_PROVENANCE_HOOK = "docs/GT_PS2_P1E_PHASE1_ENGINEERING_FOUNDATION_AMENDMENT_06.md"
AMENDMENT_06_PROSPECTIVE_ACTIVE_PROVENANCE_PATHS = (
    AMENDMENT_06_PROVENANCE_HOOK,
    "src/dtdc_simulator/core2/particle/wet_retained_cap.py",
    "src/dtdc_simulator/core2/particle/wet_water.py",
    "src/dtdc_simulator/core2/particle/cut_transport.py",
    "src/dtdc_simulator/core2/particle/cut_integrator.py",
    "src/dtdc_simulator/core2/particle/cut_continuation.py",
    "tests/test_core2_wet_retained_cap.py",
    "tests/test_core2_cut_integrator.py",
)

# Exact-A nominal target root captured after a strict accepted N=2, 0.1 s
# solve.  It is reused only as a rescaled nonlinear seed for the independent
# 0.075 s root: every campaign reconstructs its own initial conserved state
# and re-solves the complete selected equations.
# Re-earned 2026-09-01 under the corrected O9a ALE donor (owner ruling:
# "measurement-first re-earn of reference fixture"; record
# GT_PS2_V1_N2_ROOT_STRUCTURAL_DIAGNOSIS_2026-09-01.md).  The pre-O9a seed
# landed the solver in a saturating front-logit corner (exactly singular FD
# Jacobian, FAIL_CLOSED_TARGET_N2_ROOT on both campaign arms); this root is
# the accepted 0.1 s reference solve on the corrected chain (max scaled
# residual 1.13e-11, condition proxy 6.8e3; gate re-solve 1.45e-12 / 1144).
# Pre-amendment artifacts are era-bound to the prior literal, never
# retro-edited.
P1EF_REFERENCE_N2_SEED = cut.CutTransportUnknowns(
    wet_temperatures_k=(334.32323795054174, 337.16275993022816),
    wet_retained_water_loadings=(0.10001358068175062, 0.09997754874353512),
    dry_temperatures_k=(343.445024525004,),
    dry_y_hexane=(0.8166715619835563,),
    dry_total_stefan_fluxes_mol_m2_s=(
        0.13016954488527765,
        0.07035076148852501,
    ),
    front_z=0.21027104349936449,
    interface_temperature_k=337.715349918534,
)


@dataclass(frozen=True)
class InstalledEngineeringCase:
    """Complete immutable description of one process-isolated backend case."""

    coefficient_case: p1ef.EngineeringCoefficientCase
    configuration: p1ef.EngineeringParticleConfiguration
    stress_multiplier: float
    pulse_temperature_k: float
    pulse_y_hexane: float
    p1ef_smooth_retained_water_guard_enabled: bool
    vanishing_wet_cut_threshold_authority: cut.VanishingWetCutThresholdAuthority
    hot_lean: object
    hot_rich: object
    cool_lean: object
    cool_rich: object


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_source_hash(files: Mapping[str, str]) -> str:
    canonical = json.dumps(
        dict(files),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def source_provenance(repository_root: Path | None = None) -> dict:
    """Hash every source that defines or orchestrates a P1EF case."""

    root = (
        Path(__file__).resolve().parents[1]
        if repository_root is None
        else repository_root.resolve()
    )
    files: dict[str, str] = {}
    for relative in SOURCE_PROVENANCE_PATHS:
        path = (root / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise ValueError("P1EF provenance path escapes the repository") from exc
        if not path.is_file():
            raise FileNotFoundError(f"P1EF provenance source is missing: {relative}")
        files[relative] = _sha256(path)
    files = dict(sorted(files.items()))
    hook_path = (root / AMENDMENT_06_PROVENANCE_HOOK).resolve()
    try:
        hook_path.relative_to(root)
    except ValueError as exc:
        raise ValueError("P1EF Amendment-06 provenance hook escapes the repository") from exc
    hook_present = hook_path.is_file()
    return {
        "files_sha256": files,
        "combined_sha256": _canonical_source_hash(files),
        "amendment_hooks": {
            AMENDMENT_06_PROVENANCE_HOOK: {
                "present": hook_present,
                "sha256": _sha256(hook_path) if hook_present else None,
                "activation_status": ("pending_qualification_not_case_selection_authority"),
                "prospective_active_provenance_paths": list(
                    AMENDMENT_06_PROSPECTIVE_ACTIVE_PROVENANCE_PATHS
                ),
                "activation_requires": (
                    "live_N12_0p075_entry_active_continuation_exit_exact_graph_"
                    "deterministic_and_component_energy_ledgers_le_1e_10"
                ),
                "smooth_guard_supersession_scope_after_activation": (
                    "whole_trajectory_selected_before_first_step_not_per_root"
                ),
            }
        },
    }


def runtime_provenance() -> dict[str, str]:
    """Record the numerical runtime without making it an acceptance input."""

    return {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "numpy": importlib.metadata.version("numpy"),
        "scipy": importlib.metadata.version("scipy"),
    }


def _normalized_meshes(meshes: Sequence[int]) -> tuple[int, ...]:
    normalized: list[int] = []
    for value in meshes:
        try:
            cells = operator.index(value)
        except TypeError as exc:
            raise ValueError("P1EF mesh counts must be integers") from exc
        if isinstance(value, bool) or cells <= 0:
            raise ValueError("P1EF mesh counts must be positive integers")
        normalized.append(cells)
    if not normalized:
        raise ValueError("P1EF case needs at least one mesh")
    return tuple(sorted(set(normalized)))


def build_initial_condition_contract(
    installed: InstalledEngineeringCase,
    meshes: Sequence[int],
) -> dict[str, Any]:
    """Serialize the independently reconstructed manufactured ``t=0`` states."""

    normalized_meshes = _normalized_meshes(meshes)
    states_by_mesh = {
        cells: backend.initial_integrator(cells, installed.configuration.radius_m)
        for cells in normalized_meshes
    }
    snapshots, canonical_digest = initial_snapshot.serialize_mesh_snapshots(
        states_by_mesh
    )
    contract: dict[str, Any] = {
        "state_kind": initial_snapshot.STATE_KIND,
        "front_radius_fraction": backend.INITIAL_FRONT_RADIUS_FRACTION,
        "front_z": backend.INITIAL_FRONT_Z,
        "wet_temperature_k": 335.0,
        "wet_retained_water_loading": 0.10,
        "dry_temperature_k": INITIAL_DRY_TEMPERATURE_K,
        "dry_y_hexane": 0.82,
        "interface_temperature_k": 337.0,
        "historical_wet_hexane_trace_temperature_k": 335.0,
        **initial_snapshot.PROVENANCE_CLAIMS,
        "t0_snapshots_by_mesh": snapshots,
        "canonical_inventory_sha256": canonical_digest,
        "canonical_significant_decimal_digits": (
            initial_snapshot.CANONICAL_SIGNIFICANT_DECIMAL_DIGITS
        ),
        "all_mesh_snapshots_share_canonical_inventory_sha256": True,
    }
    validate_initial_condition_contract(contract)
    return contract


def validate_initial_condition_contract(value: Mapping[str, Any]) -> None:
    """Fail closed on a malformed or internally inconsistent initial contract."""

    expected_keys = {
        "state_kind",
        "front_radius_fraction",
        "front_z",
        "wet_temperature_k",
        "wet_retained_water_loading",
        "dry_temperature_k",
        "dry_y_hexane",
        "interface_temperature_k",
        "historical_wet_hexane_trace_temperature_k",
        *initial_snapshot.PROVENANCE_CLAIMS,
        "t0_snapshots_by_mesh",
        "canonical_inventory_sha256",
        "canonical_significant_decimal_digits",
        "all_mesh_snapshots_share_canonical_inventory_sha256",
    }
    if not isinstance(value, Mapping) or set(value) != expected_keys:
        raise ValueError("P1EF manufactured initial-condition fields changed")
    if value.get("state_kind") != initial_snapshot.STATE_KIND:
        raise ValueError("P1EF manufactured initial-condition kind changed")
    expected_primitives = {
        "front_radius_fraction": backend.INITIAL_FRONT_RADIUS_FRACTION,
        "front_z": backend.INITIAL_FRONT_Z,
        "wet_temperature_k": 335.0,
        "wet_retained_water_loading": 0.10,
        "dry_temperature_k": INITIAL_DRY_TEMPERATURE_K,
        "dry_y_hexane": 0.82,
        "interface_temperature_k": 337.0,
        "historical_wet_hexane_trace_temperature_k": 335.0,
    }
    if any(value.get(name) != expected for name, expected in expected_primitives.items()):
        raise ValueError("P1EF manufactured initial primitives changed")
    if any(
        value.get(name) is not expected
        for name, expected in initial_snapshot.PROVENANCE_CLAIMS.items()
    ):
        raise ValueError("P1EF manufactured initial-state claims changed")
    if (
        value.get("canonical_significant_decimal_digits")
        != initial_snapshot.CANONICAL_SIGNIFICANT_DECIMAL_DIGITS
    ):
        raise ValueError("P1EF initial-inventory canonical precision changed")
    if value.get("all_mesh_snapshots_share_canonical_inventory_sha256") is not True:
        raise ValueError("P1EF initial-inventory mesh identity is not certified")

    snapshots = value.get("t0_snapshots_by_mesh")
    if not isinstance(snapshots, Mapping) or not snapshots:
        raise ValueError("P1EF initial condition lacks mesh snapshots")
    digests: set[str] = set()
    for mesh_key, snapshot in snapshots.items():
        if not isinstance(mesh_key, str) or not mesh_key.isdecimal():
            raise ValueError("P1EF initial snapshot mesh key is invalid")
        initial_snapshot.validate_t0_snapshot(snapshot)
        if snapshot["mesh_cells"] != int(mesh_key):
            raise ValueError("P1EF initial snapshot mesh key changed")
        digests.add(snapshot["canonical_inventory_sha256"])
    if len(digests) != 1 or value.get("canonical_inventory_sha256") not in digests:
        raise ValueError("P1EF initial integrated inventory is not mesh invariant")


def attach_initial_snapshots_to_raw_jobs(
    raw: dict[str, Any],
    initial_condition: Mapping[str, Any],
) -> None:
    """Attach the exact per-mesh ``t=0`` provenance to every raw job in place."""

    validate_initial_condition_contract(initial_condition)
    jobs = raw.get("jobs") if isinstance(raw, Mapping) else None
    if not isinstance(jobs, list) or not jobs:
        raise ValueError("P1EF raw campaign lacks jobs for t=0 provenance")
    snapshots = initial_condition["t0_snapshots_by_mesh"]
    for job in jobs:
        if not isinstance(job, dict):
            raise ValueError("P1EF raw campaign job is malformed")
        try:
            cells = operator.index(job.get("cells"))
        except TypeError as exc:
            raise ValueError("P1EF raw job lacks an integer mesh count") from exc
        snapshot = snapshots.get(str(cells))
        if snapshot is None:
            raise ValueError("P1EF raw job mesh lacks a t=0 snapshot")
        initial_snapshot.validate_t0_snapshot(snapshot)
        if "initial_t0_snapshot" in job:
            raise ValueError("P1EF raw job attempted to overwrite t=0 provenance")
        profile = job.get("initial_material_volume_profile")
        if profile is not None:
            if (
                not isinstance(profile, Mapping)
                or profile.get("mesh_cells") != cells
                or profile.get("time_s") != 0.0
                or profile.get("front_z") != backend.INITIAL_FRONT_Z
            ):
                raise ValueError("P1EF raw job initial material profile is inconsistent")
        job["initial_t0_snapshot"] = deepcopy(snapshot)

    raw["initial_t0_snapshot_contract"] = {
        "state_kind": initial_snapshot.STATE_KIND,
        **initial_snapshot.PROVENANCE_CLAIMS,
        "mesh_cells": sorted(int(key) for key in snapshots),
        "canonical_inventory_sha256": initial_condition[
            "canonical_inventory_sha256"
        ],
        "all_jobs_carry_valid_per_mesh_snapshot": True,
    }


def _finite(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _film_audit_errors(
    audit: Any,
    *,
    multiplier: float,
    label: str,
) -> tuple[list[str], float | None, float | None, float | None]:
    errors: list[str] = []
    if not isinstance(audit, Mapping):
        return [f"{label} has no accepted-root surface film audit"], None, None, None
    for key in (
        "coefficient_covariance_preserved",
        "fast_mass_guard_passed",
        "validity_band_checked",
        "reduction_thresholds_enforced",
        "production_used_exact_ackermann_factor",
        "production_used_exact_phy053_partial_enthalpy_secant_heat_capacities",
        "legacy_engineering_regression_heat_capacities_are_counterfactual_only",
        "unity_ackermann_is_counterfactual_only",
        "production_used_fast_mass_surface_equals_bulk",
        "conserved_water_crosses_external_film_as_gas_once",
    ):
        if audit.get(key) is not True:
            errors.append(f"{label} film claim {key} is not true")
    if audit.get("production_used_exact_unity_ackermann_factor") is not False:
        errors.append(f"{label} promoted the rejected unity-A shortcut")
    if (
        not _finite(audit.get("stress_multiplier"))
        or float(audit["stress_multiplier"]) != multiplier
    ):
        errors.append(f"{label} changed the common film multiplier")
    mass = audit.get("mass_departure")
    departure: float | None = None
    if not isinstance(mass, Mapping):
        errors.append(f"{label} lacks the fast-mass departure audit")
    else:
        if mass.get("within_caller_threshold") is not True:
            errors.append(f"{label} failed the accepted-root fast-mass guard")
        if mass.get("supplied_fluxes_determine_departure_not_rate") is not True:
            errors.append(f"{label} let the film determine the particle flux")
        if mass.get("profile_admissibility_audited") is not True:
            errors.append(f"{label} lacks the conditioned film-profile audit")
        threshold = mass.get("caller_acceptance_threshold_mass_fraction")
        if not _finite(threshold) or float(threshold) != FAST_MASS_DEPARTURE_LIMIT:
            errors.append(f"{label} changed the 0.01 fast-mass guard")
        value = mass.get("absolute_surface_departure")
        if not _finite(value) or not 0.0 <= float(value) <= FAST_MASS_DEPARTURE_LIMIT:
            errors.append(f"{label} fast-mass departure exceeds 0.01")
        else:
            departure = float(value)
    normalized = audit.get("normalized_heat_closure_residual")
    dimensional = audit.get("heat_closure_residual_w_m2")
    roundoff = audit.get("heat_closure_roundoff_bound_w_m2")
    closure: float | None = None
    if not all(_finite(value) for value in (normalized, dimensional, roundoff)):
        errors.append(f"{label} lacks finite exact-A closure diagnostics")
    else:
        closure = abs(float(normalized))
        tight = closure <= NORMALIZED_SURFACE_HEAT_CLOSURE_LIMIT
        roundoff_branch = abs(float(dimensional)) <= 2.0 * abs(float(roundoff))
        if not (tight or roundoff_branch):
            errors.append(f"{label} exact-A scalar surface closure did not converge")
    if not (
        audit.get("monotone_unique_root_diagnostic_passed") is True
        or audit.get("exact_equal_temperature_branch") is True
    ):
        errors.append(f"{label} lacks a unique scalar surface closure")
    evaluations = audit.get("scalar_closure_function_evaluations")
    if isinstance(evaluations, bool) or not isinstance(evaluations, int) or evaluations <= 0:
        errors.append(f"{label} has an invalid scalar-closure evaluation count")
    if audit.get("full_ackermann_conductive_heat_flux_w_m2") != audit.get(
        "production_conductive_heat_flux_w_m2"
    ):
        errors.append(f"{label} production heat is not the exact-A heat")
    if audit.get("full_ackermann_surface_temperature_k") != audit.get("surface_temperature_k"):
        errors.append(f"{label} surface temperature lost exact-A identity")
    if audit.get("closed_pair_energy_residual_w_m2") != 0.0:
        errors.append(f"{label} external film energy pair does not close")
    exact_water_cp = audit.get("exact_water_partial_enthalpy_secant_heat_capacity_j_kg_k")
    exact_hexane_cp = audit.get("exact_hexane_partial_enthalpy_secant_heat_capacity_j_kg_k")
    if not all(
        _finite(value) and float(value) > 0.0 for value in (exact_water_cp, exact_hexane_cp)
    ):
        errors.append(f"{label} lacks positive exact PHY-053 secant capacities")
    else:
        try:
            recomputed_capacities = efr.phy053_partial_enthalpy_secant_heat_capacities(
                surface_temperature_k=float(audit["surface_temperature_k"]),
                bulk_temperature_k=float(audit["bulk_temperature_k"]),
                pressure_pa=float(audit["pressure_pa"]),
                y_hexane=float(audit["bulk_y_hexane"]),
                k_wh=float(audit["binary_gas_interaction_k_wh"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            errors.append(f"{label} cannot reproduce exact PHY-053 secant cp: {exc}")
        else:
            if float(exact_water_cp) != recomputed_capacities.water_j_kg_k:
                errors.append(f"{label} water cp is not the exact PHY-053 endpoint secant")
            if float(exact_hexane_cp) != recomputed_capacities.hexane_j_kg_k:
                errors.append(f"{label} hexane cp is not the exact PHY-053 endpoint secant")
            if (
                audit.get("exact_secant_equal_temperature_derivative_limit_used")
                is not recomputed_capacities.equal_temperature_derivative_limit_used
            ):
                errors.append(f"{label} exact-cp equal-temperature branch is inconsistent")
            if audit.get("exact_secant_derivative_stencil_half_width_k") != (
                recomputed_capacities.derivative_stencil_half_width_k
            ):
                errors.append(f"{label} exact-cp derivative stencil is inconsistent")
    legacy_water_cp = audit.get("legacy_engineering_regression_water_heat_capacity_j_kg_k")
    legacy_hexane_cp = audit.get("legacy_engineering_regression_hexane_heat_capacity_j_kg_k")
    if legacy_water_cp != (
        p1ef.LEGACY_ENGINEERING_REGRESSION_WATER_HEAT_CAPACITY_J_KG_K
    ) or legacy_hexane_cp != (p1ef.LEGACY_ENGINEERING_REGRESSION_HEXANE_HEAT_CAPACITY_J_KG_K):
        errors.append(f"{label} changed the declared legacy-cp counterfactual")
    counterfactual_values = (
        audit.get("legacy_engineering_regression_ackermann_factor"),
        audit.get("legacy_engineering_regression_ackermann_relative_error"),
        audit.get("legacy_engineering_regression_series_heat_transfer_w_m2_k"),
        audit.get("legacy_engineering_regression_series_conductance_relative_error"),
    )
    series_error: float | None = None
    if not all(_finite(value) and float(value) >= 0.0 for value in counterfactual_values):
        errors.append(f"{label} lacks finite legacy-cp counterfactual diagnostics")
    else:
        series_error = float(counterfactual_values[-1])
        ackermann = audit.get("ackermann")
        exact_factor = ackermann.get("factor") if isinstance(ackermann, Mapping) else None
        exact_series = audit.get("series_heat_transfer_w_m2_k")
        if not all(_finite(value) and float(value) > 0.0 for value in (exact_factor, exact_series)):
            errors.append(f"{label} lacks exact Ackermann comparison denominators")
        else:
            expected_factor_error = abs(
                float(counterfactual_values[0]) - float(exact_factor)
            ) / abs(float(exact_factor))
            expected_series_error = abs(
                float(counterfactual_values[2]) - float(exact_series)
            ) / abs(float(exact_series))
            if float(counterfactual_values[1]) != expected_factor_error:
                errors.append(f"{label} legacy-cp Ackermann relative error is inconsistent")
            if series_error != expected_series_error:
                errors.append(f"{label} legacy-cp series relative error is inconsistent")
    threshold = audit.get(
        "legacy_engineering_regression_series_conductance_relative_error_threshold"
    )
    if threshold != LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT:
        errors.append(f"{label} changed the 0.005 legacy-cp comparison threshold")
    counterfactual_pass = audit.get("legacy_engineering_regression_series_conductance_guard_passed")
    # B-F1 (2026-09-29): on post-departure roots the producer's comparison
    # yields a numpy bool (a comparison of numpy floats); it is the same
    # boolean and is normalized, not reinterpreted.  Every other type (int,
    # float, None, str) is still refused below, and the value is still checked
    # against ``series_error <= limit``.
    if isinstance(counterfactual_pass, np.bool_):
        counterfactual_pass = bool(counterfactual_pass)
    if not isinstance(counterfactual_pass, bool):
        errors.append(f"{label} lacks the legacy-cp counterfactual pass/fail result")
    elif series_error is not None and counterfactual_pass is not (
        series_error <= LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
    ):
        errors.append(f"{label} legacy-cp counterfactual pass/fail is inconsistent")
    return errors, departure, closure, series_error


def _empty_interface_saturation_endpoint_summary() -> dict[str, Any]:
    """Return explicit zero A18 accounting for a rejected evidence branch."""

    return {
        "schema_version": backend.INTERFACE_ROOT_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": backend.INTERFACE_ROOT_AMENDMENT_ID,
        "accepted_root_count": 0,
        "accepted_root_identity_count": 0,
        "rounding_plateau_root_count": 0,
        "rounding_plateau_encountered": False,
        "total_bisection_iteration_count": 0,
        "maximum_bisection_iteration_count": 0,
        "full_coupled_pore_certified_root_count": 0,
        "tolerance_changed_root_count": 0,
        "clipping_projection_or_nextafter_root_count": 0,
        "all_roots_full_coupled_pore_certification_passed": False,
        "all_roots_use_unchanged_log_fugacity_tolerance": False,
        "all_roots_avoid_clipping_projection_and_nextafter": False,
        "accepted_root_identities": [],
        "contract_passed": False,
        "physically_qualifying": False,
    }


def _empty_dry_energy_capacity_stencil_summary() -> dict[str, Any]:
    """Return explicit zero A19 accounting for a rejected evidence branch."""

    return {
        "schema_version": backend.DRY_ENERGY_CAPACITY_EVIDENCE_SCHEMA_VERSION,
        "numerical_amendment_id": backend.DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        "accepted_root_count": 0,
        "residual_scale_evaluation_count": 0,
        "dry_energy_row_count": 0,
        "capacity_evaluation_audit_count": 0,
        "stencil_kind_counts": {},
        "residual_scale_evaluation_kind_counts": {},
        "accepted_root_scale_identities": [],
        "all_capacity_evaluations_finite_positive_and_fully_certified": False,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": False,
        "physically_qualifying": False,
    }


def _failed_target_amendment_18_19_evidence(errors: Sequence[str]) -> dict[str, Any]:
    """Serialize a complete nonpassing target-root evidence contract."""

    return {
        "schema_version": 1,
        "numerical_amendment_ids": [
            backend.INTERFACE_ROOT_AMENDMENT_ID,
            backend.DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "accepted_root_count": 0,
        "accepted_root_evidence_identity_count": 0,
        "accepted_root_evidence_identities": [],
        "interface_saturation_endpoint_audit_summary": (
            _empty_interface_saturation_endpoint_summary()
        ),
        "dry_energy_capacity_stencil_audit_summary": (
            _empty_dry_energy_capacity_stencil_summary()
        ),
        "all_accepted_roots_have_complete_consistent_evidence": False,
        "contract_passed": False,
        "physically_qualifying": False,
        "errors": list(errors),
    }


def _accepted_root_amendment_18_19_evidence(
    sample: Mapping[str, Any],
) -> dict[str, Any]:
    """Revalidate and bind the exact A18/A19 children of one accepted root."""

    a18 = backend._interface_composition_root_audit_summary((sample,))  # noqa: SLF001
    a19 = backend._dry_energy_capacity_stencil_audit_summary((sample,))  # noqa: SLF001
    if not (
        a18["accepted_root_count"]
        == a18["accepted_root_identity_count"]
        == a19["accepted_root_count"]
        == 1
        and len(a18["accepted_root_identities"]) == 1
        and len(a19["accepted_root_scale_identities"]) == 1
        and a18["contract_passed"] is True
        and a19["contract_passed"] is True
    ):
        raise RuntimeError("target root A18/A19 evidence is not exactly one-to-one")
    a18_identity = a18["accepted_root_identities"][0]
    a19_identity = a19["accepted_root_scale_identities"][0]
    identity_payload = {
        "root_role": sample["root_role"],
        "balance_interval_dt_s": sample["balance_interval_dt_s"],
        "time_s": sample["time_s"],
        "a18_accepted_root_identity_sha256": a18_identity[
            "accepted_root_identity_sha256"
        ],
        "a19_accepted_root_scale_identity_sha256": a19_identity[
            "accepted_root_scale_identity_sha256"
        ],
    }
    root_identity = {
        **identity_payload,
        "accepted_root_evidence_identity_sha256": (
            backend._canonical_evidence_sha256(identity_payload)  # noqa: SLF001
        ),
    }
    return {
        "schema_version": 1,
        "numerical_amendment_ids": [
            backend.INTERFACE_ROOT_AMENDMENT_ID,
            backend.DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "accepted_root_count": 1,
        "accepted_root_evidence_identity_count": 1,
        "accepted_root_evidence_identities": [root_identity],
        "interface_saturation_endpoint_audit_summary": a18,
        "dry_energy_capacity_stencil_audit_summary": a19,
        "all_accepted_roots_have_complete_consistent_evidence": True,
        "contract_passed": True,
        "physically_qualifying": False,
        "errors": [],
    }


def _empty_release_interface_saturation_endpoint_accounting() -> dict[str, Any]:
    return {
        "target_n2_root_count": 0,
        "accepted_job_root_count": 0,
        "accepted_root_count": 0,
        "accepted_diagnostic_root_execution_count": 0,
        "evidence_root_reference_and_execution_count": 0,
        "full_coupled_pore_certified_root_count": 0,
        "rounding_plateau_root_count": 0,
        "total_bisection_iteration_count": 0,
        "maximum_bisection_iteration_count": 0,
        "tolerance_changed_root_count": 0,
        "clipping_projection_or_nextafter_root_count": 0,
        "diagnostic_full_coupled_pore_certified_root_execution_count": 0,
        "diagnostic_rounding_plateau_root_execution_count": 0,
        "diagnostic_total_bisection_iteration_count": 0,
        "diagnostic_maximum_bisection_iteration_count": 0,
        "total_full_coupled_pore_certified_root_reference_and_execution_count": 0,
        "total_rounding_plateau_root_reference_and_execution_count": 0,
        "total_bisection_iteration_count_across_references_and_executions": 0,
        "maximum_bisection_iteration_count_across_references_and_executions": 0,
        "all_roots_full_coupled_pore_certification_passed": False,
        "all_roots_use_unchanged_log_fugacity_tolerance": False,
        "all_roots_avoid_clipping_projection_and_nextafter": False,
        "contract_passed": False,
    }


def _empty_release_dry_energy_capacity_stencil_accounting() -> dict[str, Any]:
    return {
        "target_n2_root_count": 0,
        "accepted_job_root_count": 0,
        "accepted_root_count": 0,
        "accepted_diagnostic_root_execution_count": 0,
        "evidence_root_reference_and_execution_count": 0,
        "residual_scale_evaluation_count": 0,
        "dry_energy_row_count": 0,
        "capacity_evaluation_audit_count": 0,
        "diagnostic_residual_scale_evaluation_count": 0,
        "diagnostic_dry_energy_row_count": 0,
        "diagnostic_capacity_evaluation_audit_count": 0,
        "total_residual_scale_evaluation_count": 0,
        "total_dry_energy_row_count": 0,
        "total_capacity_evaluation_audit_count": 0,
        "stencil_kind_counts": {},
        "residual_scale_evaluation_kind_counts": {},
        "diagnostic_stencil_kind_counts": {},
        "diagnostic_residual_scale_evaluation_kind_counts": {},
        "total_stencil_kind_counts": {},
        "total_residual_scale_evaluation_kind_counts": {},
        "all_capacity_evaluations_finite_positive_and_fully_certified": False,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": False,
    }


def _failed_release_amendment_18_19_evidence(
    errors: Sequence[str],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Return fail-closed release evidence with no implied accepted child."""

    a18 = _empty_release_interface_saturation_endpoint_accounting()
    a19 = _empty_release_dry_energy_capacity_stencil_accounting()
    combined = {
        "schema_version": 1,
        "numerical_amendment_ids": [
            backend.INTERFACE_ROOT_AMENDMENT_ID,
            backend.DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "target_n2_root_count": 0,
        "accepted_job_count": 0,
        "accepted_step_record_count": 0,
        "accepted_root_count": 0,
        "committed_or_target_root_reference_count": 0,
        "accepted_diagnostic_root_execution_count": 0,
        "evidence_root_reference_and_execution_count": 0,
        "target_n2_root_evidence_identity_sha256": None,
        "accepted_job_evidence_identity_count": 0,
        "accepted_job_evidence_identities": [],
        "diagnostic_partition_count": 0,
        "diagnostic_partition_identity_count": 0,
        "diagnostic_partition_ids": [],
        "diagnostic_partition_identity_sha256": [],
        "diagnostic_execution_occurrence_count": 0,
        "diagnostic_execution_context_identity_count": 0,
        "diagnostic_execution_occurrence_identity_count": 0,
        "diagnostic_execution_partition_bindings": {},
        "combined_diagnostic_partition_identity_sha256": None,
        "physical_root_evidence_pair_identity_occurrence_count": 0,
        "unique_physical_root_evidence_pair_identity_count": 0,
        "physical_root_identity_uniqueness_required_across_execution_and_reference_partitions": False,
        "child_evidence_identity_count": 0,
        "child_evidence_identity_sha256": [],
        "release_evidence_identity_sha256": None,
        "interface_saturation_endpoint_audit_accounting": a18,
        "dry_energy_capacity_stencil_audit_accounting": a19,
        "all_target_job_step_and_root_evidence_bound": False,
        "accepted_state_clipping_or_projection_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": False,
        "physically_qualifying": False,
        "errors": list(errors),
    }
    return combined, a18, a19


def _merge_count_mapping(
    target: dict[str, int],
    source: Mapping[str, Any],
) -> None:
    """Accumulate a validated integer-count mapping without silent coercion."""

    for name, value in source.items():
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"evidence count {name!r} is not a nonnegative integer")
        target[name] = target.get(name, 0) + value


def _release_amendment_18_19_evidence(
    root_audit: Mapping[str, Any],
    jobs: Sequence[Any],
    raw: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], list[str]]:
    """Bind target, committed references, and noncommitted execution evidence."""

    evidence_errors: list[str] = []
    target_contract: dict[str, Any] | None = None
    try:
        a18_bundle = root_audit.get("interface_saturation_endpoint_audit")
        a19_bundle = root_audit.get("dry_energy_capacity_stencil_audit")
        if not isinstance(a18_bundle, Mapping) or not isinstance(a19_bundle, Mapping):
            raise TypeError("target N=2 root omits its A18 or A19 evidence bundle")
        accepted_root = a18_bundle.get("accepted_balance_root")
        if not isinstance(accepted_root, Mapping):
            raise TypeError("target N=2 A18 evidence omits its accepted-root envelope")
        target_sample = {
            "root_role": accepted_root.get("root_role"),
            "balance_interval_dt_s": accepted_root.get("balance_interval_dt_s"),
            "time_s": accepted_root.get("time_s"),
            "interface_saturation_endpoint_audit": a18_bundle,
            "dry_energy_capacity_stencil_audit": a19_bundle,
        }
        target_contract = _accepted_root_amendment_18_19_evidence(target_sample)
        if (
            root_audit.get("amendment_18_19_evidence_contract") != target_contract
            or root_audit.get("amendment_18_19_evidence_contract_passed") is not True
        ):
            raise RuntimeError("target N=2 A18/A19 evidence binding is not exact")
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        evidence_errors.append(str(exc))

    job_identities: list[dict[str, Any]] = []
    job_aggregates: list[dict[str, Any]] = []
    # R7: verified wall jobs are not evidence children (see backend predicate).
    # The UNFILTERED list is retained for the shared structural gate, which
    # applies the same wall filter internally — passing the pre-filtered list
    # would hide the walls from its delegation predicate and desynchronize
    # the producer and validator sides.
    unfiltered_jobs = list(jobs)
    wall_job_count = sum(
        1 for job in jobs if backend._is_r7_wall_job(job)  # noqa: SLF001
    )
    jobs = [job for job in jobs if not backend._is_r7_wall_job(job)]  # noqa: SLF001
    # R4/R7: a campaign whose every job is a verified fail-closed wall (the
    # excluded corner invocation) has no evidence children by construction;
    # the obligation is delegated to the excluded-cell fail-closed accounting.
    all_jobs_are_verified_wall_exclusions = wall_job_count > 0 and not jobs
    for job_index, job in enumerate(jobs):
        if not isinstance(job, Mapping) or job.get("completed") is not True:
            evidence_errors.append(
                f"campaign job {job_index} is not a completed A18/A19 evidence child"
            )
            continue
        records = job.get("step_records")
        try:
            if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
                raise TypeError("accepted step records are not a sequence")
            expected = backend._aggregate_amendment_18_19_evidence(records)  # noqa: SLF001
            gates = job.get("gates")
            if not isinstance(gates, Mapping):
                raise TypeError("job gates are missing")
            if not (
                job.get("amendment_18_19_evidence_contract") == expected
                and job.get("interface_saturation_endpoint_audit_accounting")
                == expected["interface_saturation_endpoint_audit_accounting"]
                and job.get("dry_energy_capacity_stencil_audit_accounting")
                == expected["dry_energy_capacity_stencil_audit_accounting"]
                and gates.get("interface_saturation_endpoint_audit_contract_passed")
                is True
                and gates.get("dry_energy_capacity_stencil_audit_contract_passed")
                is True
                and gates.get("amendment_18_19_evidence_contract_passed") is True
                and job.get("accepted_steps") == expected["accepted_step_record_count"]
                and job.get("accepted_internal_substeps")
                == expected["accepted_root_count"]
            ):
                raise RuntimeError("serialized job A18/A19 accounting is not exact")
            job_id = job.get("job_id")
            cells = job.get("cells")
            dt_s = job.get("dt_s")
            scenario = job.get("scenario")
            if not (
                isinstance(job_id, str)
                and job_id.strip()
                and isinstance(scenario, str)
                and scenario.strip()
                and isinstance(cells, int)
                and not isinstance(cells, bool)
                and cells > 0
                and _finite(dt_s)
                and float(dt_s) > 0.0
            ):
                raise ValueError("job A18/A19 identity coordinates are invalid")
            step_hashes = [
                item["accepted_step_evidence_identity_sha256"]
                for item in expected["accepted_step_evidence_identities"]
            ]
            if not (
                expected["accepted_step_evidence_identity_count"]
                == expected["accepted_step_record_count"]
                == len(step_hashes)
                and len(set(step_hashes)) == len(step_hashes)
            ):
                raise RuntimeError("job A18/A19 step identities are not one-to-one")
            identity_payload = {
                "job_id": job_id,
                "scenario": scenario,
                "cells": cells,
                "dt_s": float(dt_s),
                "accepted_step_record_count": expected["accepted_step_record_count"],
                "accepted_root_count": expected["accepted_root_count"],
                "accepted_step_evidence_identity_sha256": step_hashes,
                "interface_saturation_endpoint_audit_accounting": expected[
                    "interface_saturation_endpoint_audit_accounting"
                ],
                "dry_energy_capacity_stencil_audit_accounting": expected[
                    "dry_energy_capacity_stencil_audit_accounting"
                ],
            }
            job_identities.append(
                {
                    **identity_payload,
                    "accepted_job_evidence_identity_sha256": (
                        backend._canonical_evidence_sha256(identity_payload)  # noqa: SLF001
                    ),
                }
            )
            job_aggregates.append(expected)
        except (KeyError, TypeError, ValueError, RuntimeError) as exc:
            evidence_errors.append(
                f"campaign job {job.get('job_id', job_index)!r} A18/A19: {exc}"
            )

    if not jobs and not all_jobs_are_verified_wall_exclusions:
        evidence_errors.append("raw campaign contains no A18/A19 job evidence")
    job_ids = [item["job_id"] for item in job_identities]
    job_hashes = [
        item["accepted_job_evidence_identity_sha256"] for item in job_identities
    ]
    if len(set(job_ids)) != len(job_ids) or len(set(job_hashes)) != len(job_hashes):
        evidence_errors.append("accepted-job A18/A19 evidence identities are not one-to-one")
    if len(job_identities) != len(jobs):
        evidence_errors.append("not every accepted job has a bound A18/A19 identity")

    diagnostic_partitions: dict[str, Mapping[str, Any]] = {}
    diagnostic_contract: dict[str, Any] | None = None
    expected_diagnostic_partition_ids = (
        "bootstrap_radius_seed_executions",
        "determinism_restart_replay_executions",
    )
    try:
        raw_partitions = raw.get("amendment_18_19_diagnostic_execution_partitions")
        if not isinstance(raw_partitions, Mapping):
            raise TypeError("raw diagnostic execution partitions are missing")
        if list(raw_partitions) != list(expected_diagnostic_partition_ids):
            raise ValueError("raw diagnostic execution partition set or order changed")
        diagnostic_partitions = {
            name: raw_partitions[name] for name in expected_diagnostic_partition_ids
        }
        diagnostic_contract = (
            backend._aggregate_amendment_18_19_diagnostic_partitions(  # noqa: SLF001
                diagnostic_partitions,
                expected_partition_ids=expected_diagnostic_partition_ids,
            )
        )
        if raw.get("amendment_18_19_diagnostic_execution_contract") != (
            diagnostic_contract
        ):
            raise RuntimeError("raw diagnostic execution aggregate is detached")
        structural_gate = backend._amendment_18_19_raw_structural_gate(  # noqa: SLF001
            unfiltered_jobs,
            diagnostic_partitions,
            diagnostic_contract,
        )
        if (
            raw.get("amendment_18_19_structural_gate") != structural_gate
            or structural_gate["contract_passed"] is not True
        ):
            raise RuntimeError("raw A18/A19 structural gate is not exact")
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        evidence_errors.append(f"raw noncommitted A18/A19 evidence: {exc}")

    if evidence_errors or target_contract is None or diagnostic_contract is None:
        combined, a18, a19 = _failed_release_amendment_18_19_evidence(
            evidence_errors
        )
        return combined, a18, a19, evidence_errors

    target_a18 = target_contract["interface_saturation_endpoint_audit_summary"]
    target_a19 = target_contract["dry_energy_capacity_stencil_audit_summary"]
    target_identity = target_contract["accepted_root_evidence_identities"][0][
        "accepted_root_evidence_identity_sha256"
    ]
    accepted_job_root_count = sum(
        int(item["accepted_root_count"]) for item in job_aggregates
    )
    accepted_step_record_count = sum(
        int(item["accepted_step_record_count"]) for item in job_aggregates
    )
    total_root_count = 1 + accepted_job_root_count
    diagnostic_root_count = int(
        diagnostic_contract["accepted_diagnostic_root_execution_count"]
    )
    diagnostic_a18 = diagnostic_contract[
        "interface_saturation_endpoint_audit_accounting"
    ]
    diagnostic_a19 = diagnostic_contract[
        "dry_energy_capacity_stencil_audit_accounting"
    ]
    a18 = {
        "target_n2_root_count": 1,
        "accepted_job_root_count": accepted_job_root_count,
        "accepted_root_count": total_root_count,
        "accepted_diagnostic_root_execution_count": diagnostic_root_count,
        "evidence_root_reference_and_execution_count": (
            total_root_count + diagnostic_root_count
        ),
        "full_coupled_pore_certified_root_count": (
            int(target_a18["full_coupled_pore_certified_root_count"])
            + sum(
                int(item["interface_saturation_endpoint_audit_accounting"]
                         ["full_coupled_pore_certified_root_count"])
                for item in job_aggregates
            )
        ),
        "rounding_plateau_root_count": (
            int(target_a18["rounding_plateau_root_count"])
            + sum(
                int(item["interface_saturation_endpoint_audit_accounting"]
                         ["rounding_plateau_root_count"])
                for item in job_aggregates
            )
        ),
        "total_bisection_iteration_count": (
            int(target_a18["total_bisection_iteration_count"])
            + sum(
                int(item["interface_saturation_endpoint_audit_accounting"]
                         ["total_bisection_iteration_count"])
                for item in job_aggregates
            )
        ),
        "maximum_bisection_iteration_count": max(
            [int(target_a18["maximum_bisection_iteration_count"])]
            + [
                int(item["interface_saturation_endpoint_audit_accounting"]
                         ["maximum_bisection_iteration_count"])
                for item in job_aggregates
            ]
        ),
        "tolerance_changed_root_count": 0,
        "clipping_projection_or_nextafter_root_count": 0,
        "diagnostic_full_coupled_pore_certified_root_execution_count": (
            diagnostic_a18[
                "full_coupled_pore_certified_root_execution_count"
            ]
        ),
        "diagnostic_rounding_plateau_root_execution_count": diagnostic_a18[
            "rounding_plateau_root_execution_count"
        ],
        "diagnostic_total_bisection_iteration_count": diagnostic_a18[
            "total_bisection_iteration_count"
        ],
        "diagnostic_maximum_bisection_iteration_count": diagnostic_a18[
            "maximum_bisection_iteration_count"
        ],
        "total_full_coupled_pore_certified_root_reference_and_execution_count": (
            total_root_count
            + int(
                diagnostic_a18[
                    "full_coupled_pore_certified_root_execution_count"
                ]
            )
        ),
        "total_rounding_plateau_root_reference_and_execution_count": (
            int(target_a18["rounding_plateau_root_count"])
            + sum(
                int(
                    item["interface_saturation_endpoint_audit_accounting"][
                        "rounding_plateau_root_count"
                    ]
                )
                for item in job_aggregates
            )
            + int(diagnostic_a18["rounding_plateau_root_execution_count"])
        ),
        "total_bisection_iteration_count_across_references_and_executions": (
            int(target_a18["total_bisection_iteration_count"])
            + sum(
                int(
                    item["interface_saturation_endpoint_audit_accounting"][
                        "total_bisection_iteration_count"
                    ]
                )
                for item in job_aggregates
            )
            + int(diagnostic_a18["total_bisection_iteration_count"])
        ),
        "maximum_bisection_iteration_count_across_references_and_executions": max(
            int(target_a18["maximum_bisection_iteration_count"]),
            *(
                int(
                    item["interface_saturation_endpoint_audit_accounting"][
                        "maximum_bisection_iteration_count"
                    ]
                )
                for item in job_aggregates
            ),
            int(diagnostic_a18["maximum_bisection_iteration_count"]),
        ),
        "all_roots_full_coupled_pore_certification_passed": True,
        "all_roots_use_unchanged_log_fugacity_tolerance": True,
        "all_roots_avoid_clipping_projection_and_nextafter": True,
        "contract_passed": True,
    }
    stencil_kind_counts: dict[str, int] = {}
    evaluation_kind_counts: dict[str, int] = {}
    _merge_count_mapping(stencil_kind_counts, target_a19["stencil_kind_counts"])
    _merge_count_mapping(
        evaluation_kind_counts,
        target_a19["residual_scale_evaluation_kind_counts"],
    )
    for aggregate in job_aggregates:
        job_a19 = aggregate["dry_energy_capacity_stencil_audit_accounting"]
        _merge_count_mapping(stencil_kind_counts, job_a19["stencil_kind_counts"])
        _merge_count_mapping(
            evaluation_kind_counts,
            job_a19["residual_scale_evaluation_kind_counts"],
        )
    total_stencil_kind_counts = dict(stencil_kind_counts)
    total_evaluation_kind_counts = dict(evaluation_kind_counts)
    _merge_count_mapping(
        total_stencil_kind_counts,
        diagnostic_a19["stencil_kind_counts"],
    )
    _merge_count_mapping(
        total_evaluation_kind_counts,
        diagnostic_a19["residual_scale_evaluation_kind_counts"],
    )
    a19 = {
        "target_n2_root_count": 1,
        "accepted_job_root_count": accepted_job_root_count,
        "accepted_root_count": total_root_count,
        "accepted_diagnostic_root_execution_count": diagnostic_root_count,
        "evidence_root_reference_and_execution_count": (
            total_root_count + diagnostic_root_count
        ),
        "residual_scale_evaluation_count": (
            int(target_a19["residual_scale_evaluation_count"])
            + sum(
                int(item["dry_energy_capacity_stencil_audit_accounting"]
                         ["residual_scale_evaluation_count"])
                for item in job_aggregates
            )
        ),
        "dry_energy_row_count": (
            int(target_a19["dry_energy_row_count"])
            + sum(
                int(item["dry_energy_capacity_stencil_audit_accounting"]
                         ["dry_energy_row_count"])
                for item in job_aggregates
            )
        ),
        "capacity_evaluation_audit_count": (
            int(target_a19["capacity_evaluation_audit_count"])
            + sum(
                int(item["dry_energy_capacity_stencil_audit_accounting"]
                         ["capacity_evaluation_audit_count"])
                for item in job_aggregates
            )
        ),
        "diagnostic_residual_scale_evaluation_count": diagnostic_a19[
            "residual_scale_evaluation_count"
        ],
        "diagnostic_dry_energy_row_count": diagnostic_a19[
            "dry_energy_row_count"
        ],
        "diagnostic_capacity_evaluation_audit_count": diagnostic_a19[
            "capacity_evaluation_audit_count"
        ],
        "total_residual_scale_evaluation_count": (
            int(target_a19["residual_scale_evaluation_count"])
            + sum(
                int(
                    item["dry_energy_capacity_stencil_audit_accounting"][
                        "residual_scale_evaluation_count"
                    ]
                )
                for item in job_aggregates
            )
            + int(diagnostic_a19["residual_scale_evaluation_count"])
        ),
        "total_dry_energy_row_count": (
            int(target_a19["dry_energy_row_count"])
            + sum(
                int(
                    item["dry_energy_capacity_stencil_audit_accounting"][
                        "dry_energy_row_count"
                    ]
                )
                for item in job_aggregates
            )
            + int(diagnostic_a19["dry_energy_row_count"])
        ),
        "total_capacity_evaluation_audit_count": (
            int(target_a19["capacity_evaluation_audit_count"])
            + sum(
                int(
                    item["dry_energy_capacity_stencil_audit_accounting"][
                        "capacity_evaluation_audit_count"
                    ]
                )
                for item in job_aggregates
            )
            + int(diagnostic_a19["capacity_evaluation_audit_count"])
        ),
        "stencil_kind_counts": dict(sorted(stencil_kind_counts.items())),
        "residual_scale_evaluation_kind_counts": dict(
            sorted(evaluation_kind_counts.items())
        ),
        "diagnostic_stencil_kind_counts": diagnostic_a19[
            "stencil_kind_counts"
        ],
        "diagnostic_residual_scale_evaluation_kind_counts": diagnostic_a19[
            "residual_scale_evaluation_kind_counts"
        ],
        "total_stencil_kind_counts": dict(
            sorted(total_stencil_kind_counts.items())
        ),
        "total_residual_scale_evaluation_kind_counts": dict(
            sorted(total_evaluation_kind_counts.items())
        ),
        "all_capacity_evaluations_finite_positive_and_fully_certified": True,
        "accepted_residual_equations_changed": False,
        "composition_or_activity_clipping_used": False,
        "property_extrapolation_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": True,
    }
    diagnostic_partition_bindings = {
        partition_id: {
            "partition_identity_sha256": partition[
                "partition_identity_sha256"
            ],
            "execution_occurrence_count": partition[
                "execution_occurrence_count"
            ],
            "execution_context_identity_count": partition[
                "execution_context_identity_count"
            ],
            "execution_occurrence_identity_count": partition[
                "execution_occurrence_identity_count"
            ],
            "accepted_diagnostic_root_execution_count": partition[
                "accepted_diagnostic_root_execution_count"
            ],
            "committed_trajectory_root_count": 0,
            "execution_context_identity_sha256": partition[
                "execution_context_identity_sha256"
            ],
            "execution_occurrence_identity_sha256": partition[
                "execution_occurrence_identity_sha256"
            ],
            "may_reference_committed_trajectory_root_elsewhere": partition[
                "may_reference_committed_trajectory_root_elsewhere"
            ],
            "interface_saturation_endpoint_audit_accounting": partition[
                "interface_saturation_endpoint_audit_accounting"
            ],
            "dry_energy_capacity_stencil_audit_accounting": partition[
                "dry_energy_capacity_stencil_audit_accounting"
            ],
        }
        for partition_id, partition in diagnostic_partitions.items()
    }
    diagnostic_partition_hashes = diagnostic_contract[
        "partition_identity_sha256"
    ]
    child_hashes = [
        target_identity,
        *job_hashes,
        *diagnostic_partition_hashes,
    ]
    release_identity_payload = {
        "target_n2_root_evidence_identity_sha256": target_identity,
        "accepted_job_evidence_identity_sha256": job_hashes,
        "accepted_step_record_count": accepted_step_record_count,
        "accepted_root_count": total_root_count,
        "accepted_diagnostic_root_execution_count": diagnostic_root_count,
        "diagnostic_partition_identity_sha256": diagnostic_partition_hashes,
        "combined_diagnostic_partition_identity_sha256": diagnostic_contract[
            "combined_diagnostic_partition_identity_sha256"
        ],
        "interface_saturation_endpoint_audit_accounting": a18,
        "dry_energy_capacity_stencil_audit_accounting": a19,
    }
    combined = {
        "schema_version": 1,
        "numerical_amendment_ids": [
            backend.INTERFACE_ROOT_AMENDMENT_ID,
            backend.DRY_ENERGY_CAPACITY_AMENDMENT_ID,
        ],
        "target_n2_root_count": 1,
        "accepted_job_count": len(job_identities),
        "accepted_step_record_count": accepted_step_record_count,
        "accepted_root_count": total_root_count,
        "committed_or_target_root_reference_count": total_root_count,
        "accepted_diagnostic_root_execution_count": diagnostic_root_count,
        "evidence_root_reference_and_execution_count": (
            total_root_count + diagnostic_root_count
        ),
        "target_n2_root_evidence_identity_sha256": target_identity,
        "accepted_job_evidence_identity_count": len(job_identities),
        "accepted_job_evidence_identities": job_identities,
        "diagnostic_partition_count": diagnostic_contract["partition_count"],
        "diagnostic_partition_identity_count": diagnostic_contract[
            "partition_identity_count"
        ],
        "diagnostic_partition_ids": diagnostic_contract["partition_ids"],
        "diagnostic_partition_identity_sha256": diagnostic_partition_hashes,
        "diagnostic_execution_occurrence_count": diagnostic_contract[
            "execution_occurrence_count"
        ],
        "diagnostic_execution_context_identity_count": diagnostic_contract[
            "execution_context_identity_count"
        ],
        "diagnostic_execution_occurrence_identity_count": diagnostic_contract[
            "execution_occurrence_identity_count"
        ],
        "diagnostic_execution_partition_bindings": (
            diagnostic_partition_bindings
        ),
        "combined_diagnostic_partition_identity_sha256": diagnostic_contract[
            "combined_diagnostic_partition_identity_sha256"
        ],
        "physical_root_evidence_pair_identity_occurrence_count": (
            diagnostic_contract[
                "physical_root_evidence_pair_identity_occurrence_count"
            ]
        ),
        "unique_physical_root_evidence_pair_identity_count": (
            diagnostic_contract[
                "unique_physical_root_evidence_pair_identity_count"
            ]
        ),
        "physical_root_identity_uniqueness_required_across_execution_and_reference_partitions": False,
        "child_evidence_identity_count": len(child_hashes),
        "child_evidence_identity_sha256": child_hashes,
        "release_evidence_identity_sha256": (
            backend._canonical_evidence_sha256(release_identity_payload)  # noqa: SLF001
        ),
        "interface_saturation_endpoint_audit_accounting": a18,
        "dry_energy_capacity_stencil_audit_accounting": a19,
        "all_target_job_step_and_root_evidence_bound": True,
        "accepted_state_clipping_or_projection_used": False,
        "acceptance_tolerances_changed": False,
        "contract_passed": True,
        "physically_qualifying": False,
        "errors": [],
    }
    return combined, a18, a19, []


def amendment_17_18_19_case_contract_errors(
    payload: Mapping[str, Any],
) -> list[str]:
    """Independently revalidate one serialized case's A17--A19 evidence."""

    if not isinstance(payload, Mapping):
        return ["case artifact is not a mapping"]

    errors: list[str] = []
    root_audit = payload.get("target_n2_root")
    raw = payload.get("raw_backend")
    release = payload.get("release_evidence_contract")
    if not isinstance(root_audit, Mapping):
        errors.append("case artifact omits its target N=2 root audit")
    if not isinstance(raw, Mapping):
        errors.append("case artifact omits its raw backend evidence")
    if not isinstance(release, Mapping):
        errors.append("case artifact omits its release evidence contract")
    if errors:
        return errors

    solver_contract = raw.get("solver_acceptance_contract")
    if not isinstance(solver_contract, Mapping):
        errors.append("raw solver acceptance contract is missing")
    else:
        if (
            solver_contract.get("maximum_condition_proxy")
            != P1EF_MAXIMUM_CONDITION_PROXY
        ):
            errors.append(
                "raw solver maximum condition proxy is not the exact frozen P1EF ceiling"
            )
        exact_solver_declarations = (
            (
                failed_leaf_contract.FACE_CONDITIONED_PREDICTOR_CONTRACT_KEY,
                failed_leaf_contract.face_conditioned_predictor_contract(),
                "A17 face-conditioned predictor",
            ),
            (
                failed_leaf_contract.INTERFACE_ROOT_CONTRACT_KEY,
                failed_leaf_contract.interface_root_contract(),
                "A18 interface-root selection",
            ),
            (
                failed_leaf_contract.DRY_ENERGY_CAPACITY_CONTRACT_KEY,
                failed_leaf_contract.dry_energy_capacity_stencil_contract(),
                "A19 dry-energy normalization",
            ),
        )
        for key, expected, label in exact_solver_declarations:
            value = solver_contract.get(key)
            if not isinstance(value, Mapping) or dict(value) != expected:
                errors.append(f"raw solver {label} declaration is not exact")

    expected_predictor = failed_leaf_contract.face_conditioned_predictor_contract()
    exposed_predictor = release.get("face_conditioned_dry_composition_predictor")
    if (
        not isinstance(exposed_predictor, Mapping)
        or dict(exposed_predictor) != expected_predictor
    ):
        errors.append("release A17 face-conditioned predictor binding is not exact")
    release_solver = release.get("solver")
    if (
        not isinstance(release_solver, Mapping)
        or release_solver.get(
            "face_conditioned_dry_composition_predictor_contract_passed"
        )
        is not True
    ):
        errors.append("release A17 pass gate is absent or nonpassing")

    jobs = raw.get("jobs")
    if not isinstance(jobs, list):
        errors.append("raw campaign jobs are not a serialized list")
        return errors
    try:
        expected_combined, expected_a18, expected_a19, evidence_errors = (
            _release_amendment_18_19_evidence(root_audit, jobs, raw)
        )
    except (KeyError, TypeError, ValueError, RuntimeError, OverflowError) as exc:
        errors.append(f"A18/A19 evidence recomputation failed closed: {exc}")
        return errors
    errors.extend(f"A18/A19 evidence: {error}" for error in evidence_errors)
    if release.get("amendment_18_19_evidence") != expected_combined:
        errors.append("release combined A18/A19 evidence is not the exact recomputation")
    if release.get("interface_saturation_endpoint_audit") != expected_a18:
        errors.append("release A18 accounting is not the exact recomputation")
    if release.get("dry_energy_capacity_stencil_audit") != expected_a19:
        errors.append("release A19 accounting is not the exact recomputation")
    if expected_combined.get("contract_passed") is not True:
        errors.append("recomputed combined A18/A19 evidence is nonpassing")
    if expected_a18.get("contract_passed") is not True:
        errors.append("recomputed A18 accounting is nonpassing")
    if expected_a19.get("contract_passed") is not True:
        errors.append("recomputed A19 accounting is nonpassing")
    return errors


def _failed_release_evidence_contract(errors: Sequence[str]) -> dict:
    """Return a complete nonpassing contract instead of omitting evidence."""

    amendment, interface_endpoint, dry_capacity = (
        _failed_release_amendment_18_19_evidence(errors)
    )
    return {
        "contract_passed": False,
        "errors": list(errors),
        "accepted_residuals": {
            "all_square_and_full_rank": False,
            "basis": "not available because the case did not complete",
        },
        "rollback": {
            "every_failed_attempt_preserved_input_identity": False,
            "no_failed_state_committed": False,
        },
        "solver": {
            "no_clipping_or_bound_projection": False,
            "no_coefficient_floor": False,
            "no_penalty_source": False,
            "no_tolerance_relaxation": False,
            "no_equation_changing_fallback": False,
            "no_topology_selection_from_rejected_trial": False,
            "frozen_failed_leaf_controller_contract_passed": False,
            "face_conditioned_dry_composition_predictor_contract_passed": False,
            "accepted_subdivision_policy": failed_leaf_contract.POLICY_ID,
            "arbitrary_accepted_subdivision_allowed": False,
        },
        "scalar_post_root_polish": {
            "numerical_amendment_id": (backend.SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID),
            "used_step_count": 0,
            "used_accepted_root_count": 0,
            "serialized_job_occurrence_step_use_count": 0,
            "serialized_job_occurrence_root_use_count": 0,
            "serialized_job_occurrence_residual_evaluations": 0,
            "maximum_scaled_residual_before": None,
            "maximum_scaled_residual_after": None,
            "total_scalar_polish_residual_evaluations": 0,
            "all_uses_fully_serialized": False,
            "absence_means_zero_uses": True,
            "physics_equations_or_acceptance_tolerances_changed": False,
        },
        "film": {
            "all_accepted_roots_audited": False,
            "total_accepted_root_count": 0,
            "audited_accepted_root_count": 0,
            "fast_mass_surface_equals_bulk": True,
            "fast_mass_guard_enforced_at_accepted_root": False,
            "fast_mass_departure_threshold": FAST_MASS_DEPARTURE_LIMIT,
            "maximum_absolute_mass_fraction_departure": None,
            "production_used_exact_ackermann_factor": True,
            "production_used_exact_phy053_partial_enthalpy_secant_heat_capacities": False,
            "legacy_engineering_regression_heat_capacities_counterfactual_only": True,
            "legacy_cp_counterfactual_comparison_audited_at_all_accepted_roots": False,
            "legacy_cp_series_conductance_relative_error_threshold": (
                LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
            ),
            "maximum_legacy_cp_series_conductance_relative_error": None,
            "all_legacy_cp_counterfactuals_within_0_005": False,
            "legacy_cp_counterfactual_is_release_acceptance_guard": False,
            "unity_ackermann_is_counterfactual_only": True,
            "unity_ackermann_shortcut_used": False,
            "a1_counterfactual_promoted_to_production": False,
            "exact_scalar_surface_closure_converged": False,
            "normalized_surface_heat_closure_limit": (NORMALIZED_SURFACE_HEAT_CLOSURE_LIMIT),
            "maximum_normalized_surface_heat_closure_residual": None,
            "exact_ackermann_series_identity_verified": False,
            "monotone_unique_surface_root_diagnostic_passed": False,
            "common_heat_mass_multiplier_preserved": False,
            "external_liquid_saturation_zero": True,
        },
        "amendment_18_19_evidence": amendment,
        "interface_saturation_endpoint_audit": interface_endpoint,
        "dry_energy_capacity_stencil_audit": dry_capacity,
        "face_conditioned_dry_composition_predictor": None,
    }


def _scalar_post_root_audit_metrics(
    value: Any,
    *,
    label: str,
    errors: list[str],
) -> tuple[float, float, int] | None:
    """Validate one serialized use and return its before/after/work metrics."""

    if not isinstance(value, Mapping):
        errors.append(f"{label} scalar post-root audit is not an object")
        return None
    required_exact = {
        "near_miss_multiplier_ceiling": (ci.SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER),
        "dominant_residual_block": "dry_water",
        "selected_coordinate_block": "dry_y_hexane",
        "used_generic_interface_tolerance_unchanged": True,
        "no_clipping_projection_or_tolerance_relaxation": True,
        "same_topology_and_rank": True,
        "partial_candidate_committed": False,
        "exact_rollback_on_failure": True,
    }
    for name, expected in required_exact.items():
        if value.get(name) != expected:
            errors.append(f"{label} scalar post-root field {name} is contradictory")
    before = value.get("maximum_scaled_residual_before")
    after = value.get("maximum_scaled_residual_after")
    condition = value.get("final_condition_proxy")
    if not (
        _finite(before)
        and backend.CONTROLS.nonlinear_residual_tolerance
        < float(before)
        <= ci.scalar_post_root_near_miss_ceiling(
            backend.CONTROLS.nonlinear_residual_tolerance
        )
    ):
        errors.append(f"{label} scalar post-root near-miss band is invalid")
    if not (
        _finite(after) and 0.0 <= float(after) <= backend.CONTROLS.nonlinear_residual_tolerance
    ):
        errors.append(f"{label} scalar post-root final full residual is invalid")
    if not (
        _finite(condition) and 0.0 < float(condition) <= backend.CONTROLS.maximum_condition_proxy
    ):
        errors.append(f"{label} scalar post-root final condition is invalid")
    rank = value.get("final_jacobian_rank")
    dimension = value.get("final_jacobian_dimension")
    if (
        isinstance(rank, bool)
        or not isinstance(rank, int)
        or isinstance(dimension, bool)
        or not isinstance(dimension, int)
        or rank <= 0
        or rank != dimension
    ):
        errors.append(f"{label} scalar post-root final Jacobian is not full rank")
    counters = (
        "selection_stencil_evaluations",
        "bracket_probe_count",
        "scalar_function_calls",
        "final_condition_stencil_evaluations",
    )
    work_values: list[int] = []
    for name in counters:
        item = value.get(name)
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            errors.append(f"{label} scalar post-root work field {name} is invalid")
        else:
            work_values.append(item)
    if len(work_values) != len(counters) or not (_finite(before) and _finite(after)):
        return None
    return float(before), float(after), sum(work_values) + 1


def _scalar_post_root_step_summary_metrics(
    value: Any,
    *,
    label: str,
    errors: list[str],
) -> tuple[int, int, float | None, float | None, int] | None:
    """Validate an explicit accepted-step summary; zero-use is represented, not omitted."""

    if not isinstance(value, Mapping):
        errors.append(f"{label} lacks explicit scalar post-root serialization")
        return None
    if value.get("numerical_amendment_id") != (backend.SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID):
        errors.append(f"{label} scalar post-root amendment provenance is wrong")
    if value.get("absence_means_zero_uses") is not True:
        errors.append(f"{label} scalar post-root absence semantics are not explicit")
    if value.get("physics_equations_or_acceptance_tolerances_changed") is not False:
        errors.append(f"{label} scalar post-root serialization claims a model change")
    accepted = value.get("accepted_root_count_considered")
    used = value.get("used_accepted_root_count")
    work = value.get("total_scalar_polish_residual_evaluations")
    if any(
        isinstance(item, bool) or not isinstance(item, int) or item < 0
        for item in (accepted, used, work)
    ) or not (isinstance(accepted, int) and isinstance(used, int) and used <= accepted):
        errors.append(f"{label} scalar post-root counts are invalid")
        return None
    if value.get("used_step") is not (used > 0):
        errors.append(f"{label} scalar post-root used-step flag is inconsistent")
    serialized = value.get("audits")
    if not isinstance(serialized, list) or len(serialized) != used:
        errors.append(f"{label} scalar post-root audit count is inconsistent")
        return None
    metrics: list[tuple[float, float, int]] = []
    for index, entry in enumerate(serialized):
        if not isinstance(entry, Mapping):
            errors.append(f"{label} scalar post-root audit entry {index} is invalid")
            continue
        metric = _scalar_post_root_audit_metrics(
            entry.get("audit"),
            label=f"{label} root {index}",
            errors=errors,
        )
        if metric is not None:
            metrics.append(metric)
    before = value.get("maximum_scaled_residual_before")
    after = value.get("maximum_scaled_residual_after")
    expected_before = max((item[0] for item in metrics), default=None)
    expected_after = max((item[1] for item in metrics), default=None)
    expected_work = sum(item[2] for item in metrics)
    if before != expected_before or after != expected_after or work != expected_work:
        errors.append(f"{label} scalar post-root maxima/work do not match its audits")
    return accepted, used, expected_before, expected_after, expected_work


def _scalar_post_root_accounting_metrics(
    value: Any,
    *,
    label: str,
    errors: list[str],
) -> tuple[int, int, int, int, float | None, float | None, int] | None:
    """Validate one aggregate accounting payload."""

    if not isinstance(value, Mapping):
        errors.append(f"{label} lacks scalar post-root aggregate accounting")
        return None
    exact = {
        "numerical_amendment_id": backend.SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID,
        "all_steps_explicitly_serialized": True,
        "absence_means_zero_uses": True,
        "physics_equations_or_acceptance_tolerances_changed": False,
    }
    for name, expected in exact.items():
        if value.get(name) != expected:
            errors.append(f"{label} scalar post-root aggregate field {name} is wrong")
    names = (
        "accepted_step_record_count_considered",
        "used_step_record_count",
        "accepted_root_count_considered",
        "used_accepted_root_count",
        "total_scalar_polish_residual_evaluations",
    )
    counts = [value.get(name) for name in names]
    if any(isinstance(item, bool) or not isinstance(item, int) or item < 0 for item in counts):
        errors.append(f"{label} scalar post-root aggregate counts are invalid")
        return None
    accepted_steps, used_steps, accepted_roots, used_roots, work = counts
    if used_steps > accepted_steps or used_roots > accepted_roots or used_steps > used_roots:
        errors.append(f"{label} scalar post-root aggregate counts are inconsistent")
    before = value.get("maximum_scaled_residual_before")
    after = value.get("maximum_scaled_residual_after")
    if used_roots == 0:
        if before is not None or after is not None or work != 0:
            errors.append(f"{label} zero scalar uses carry nonzero metrics")
    elif not (
        _finite(before)
        and backend.CONTROLS.nonlinear_residual_tolerance
        < float(before)
        <= ci.scalar_post_root_near_miss_ceiling(
            backend.CONTROLS.nonlinear_residual_tolerance
        )
        and _finite(after)
        and 0.0 <= float(after) <= backend.CONTROLS.nonlinear_residual_tolerance
        and work > 0
    ):
        errors.append(f"{label} scalar post-root aggregate metrics are invalid")
    return (
        accepted_steps,
        used_steps,
        accepted_roots,
        used_roots,
        None if before is None else float(before),
        None if after is None else float(after),
        work,
    )


def build_release_evidence_contract(
    installed: InstalledEngineeringCase,
    root_audit: Mapping[str, Any],
    raw: Mapping[str, Any],
) -> dict:
    """Derive the release summary from accepted records and solver audits."""

    errors: list[str] = []
    film_errors: list[str] = []
    audits: list[Mapping[str, Any]] = []
    departures: list[float] = []
    closures: list[float] = []
    legacy_cp_series_errors: list[float] = []
    expected_audit_count = 0
    scalar_errors: list[str] = []
    scalar_before_values: list[float] = []
    scalar_after_values: list[float] = []
    scalar_root_use_count = 0
    scalar_step_use_count = 0
    scalar_work = 0
    serialized_job_scalar_root_use_count = 0
    serialized_job_scalar_step_use_count = 0
    serialized_job_scalar_work = 0
    surface_constitutive_errors: list[str] = []
    surface_constitutive_root_count = 0
    surface_constitutive_job_count = 0

    if "scalar_post_root_polish_audit" not in root_audit:
        scalar_errors.append(
            "target N=2 root omits explicit scalar post-root zero/use serialization"
        )
    root_scalar_audit = root_audit.get("scalar_post_root_polish_audit")
    if root_scalar_audit is not None:
        root_metric = _scalar_post_root_audit_metrics(
            root_scalar_audit,
            label="target N=2 root",
            errors=scalar_errors,
        )
        if root_metric is not None:
            scalar_before_values.append(root_metric[0])
            scalar_after_values.append(root_metric[1])
            scalar_work += root_metric[2]
            scalar_root_use_count += 1
            scalar_step_use_count += 1

    def append_audit(value: Any, label: str) -> None:
        nonlocal expected_audit_count
        expected_audit_count += 1
        audit_errors, departure, closure, legacy_cp_series_error = _film_audit_errors(
            value,
            multiplier=installed.stress_multiplier,
            label=label,
        )
        film_errors.extend(audit_errors)
        errors.extend(audit_errors)
        if isinstance(value, Mapping):
            audits.append(value)
        if departure is not None:
            departures.append(departure)
        if closure is not None:
            closures.append(closure)
        if legacy_cp_series_error is not None:
            legacy_cp_series_errors.append(legacy_cp_series_error)

    append_audit(root_audit.get("surface_film_audit"), "target N=2 root")
    jobs = raw.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        errors.append("raw campaign contains no accepted jobs")
        jobs = []
    (
        amendment_18_19_evidence,
        interface_saturation_endpoint_audit,
        dry_energy_capacity_stencil_audit,
        amendment_18_19_errors,
    ) = _release_amendment_18_19_evidence(root_audit, jobs, raw)
    errors.extend(amendment_18_19_errors)
    condition_proxies: list[float] = []
    root_attempt = root_audit.get("attempt")
    attempts: list[Mapping[str, Any]] = []
    if isinstance(root_attempt, Mapping):
        attempts.append(root_attempt)
        if _finite(root_attempt.get("condition_proxy")):
            condition_proxies.append(float(root_attempt["condition_proxy"]))
    else:
        errors.append("target root lacks its attempt audit")
    step_records: list[Mapping[str, Any]] = []
    for job in jobs:
        # R7/R12: verified wall jobs are excluded-cell evidence with their
        # own fail-closed accounting; they are not release children.
        if backend._is_r7_wall_job(job):  # noqa: SLF001
            continue
        if not isinstance(job, Mapping) or job.get("completed") is not True:
            errors.append("release contract encountered an incomplete campaign job")
            continue
        records = job.get("step_records")
        if not isinstance(records, list) or not records:
            errors.append("completed campaign job lacks accepted step records")
            continue
        scenario = backend.SCENARIOS.get(job.get("scenario"))
        if scenario is None:
            surface_constitutive_errors.append(
                "completed campaign job has no declared A15 scenario"
            )
        else:
            recomputed_surface_contract = (
                backend._surface_constitutive_oracle_contract(  # noqa: SLF001
                    records,
                    scenario,
                )
            )
            if (
                recomputed_surface_contract.get("contract_passed") is not True
                or job.get("surface_constitutive_oracle_contract")
                != recomputed_surface_contract
                or not isinstance(job.get("mechanism_specific_contract"), Mapping)
                or job["mechanism_specific_contract"].get(
                    "surface_constitutive_oracle_contract"
                )
                != recomputed_surface_contract
                or not isinstance(job.get("gates"), Mapping)
                or job["gates"].get(
                    "surface_constitutive_oracle_contract_passed"
                )
                is not True
            ):
                surface_constitutive_errors.append(
                    f"campaign job {job.get('job_id')} failed A15 recomputation/binding"
                )
            else:
                surface_constitutive_job_count += 1
                surface_constitutive_root_count += int(
                    recomputed_surface_contract["accepted_root_count"]
                )
        for record_index, record in enumerate(records):
            if not isinstance(record, Mapping):
                errors.append("accepted step record is not an object")
                continue
            step_records.append(record)
            if _finite(record.get("condition_proxy")):
                condition_proxies.append(float(record["condition_proxy"]))
            raw_attempts = record.get("attempts")
            if isinstance(raw_attempts, list):
                attempts.extend(attempt for attempt in raw_attempts if isinstance(attempt, Mapping))
            else:
                errors.append("accepted step lacks attempt/rollback records")
            scalar_summary = record.get("scalar_post_root_polish")
            scalar_metric = _scalar_post_root_step_summary_metrics(
                scalar_summary,
                label=f"accepted job step {record_index}",
                errors=scalar_errors,
            )
            if scalar_metric is not None:
                _, used_roots, _, _, work = scalar_metric
                serialized_job_scalar_root_use_count += used_roots
                serialized_job_scalar_step_use_count += int(used_roots > 0)
                serialized_job_scalar_work += work
            samples = record.get("accepted_node_surface_flux_samples")
            count = record.get("accepted_internal_substep_count")
            if (
                not isinstance(samples, list)
                or isinstance(count, bool)
                or not isinstance(count, int)
                or count <= 0
                or len(samples) != count
            ):
                errors.append("accepted step does not audit every internal root")
                continue
            sample_audits: list[Mapping[str, Any]] = []
            sample_scalar_audits: list[Mapping[str, Any]] = []
            for sample_index, sample in enumerate(samples):
                value = sample.get("surface_film_audit") if isinstance(sample, Mapping) else None
                append_audit(
                    value,
                    f"accepted job root {record_index}/{sample_index}",
                )
                if isinstance(value, Mapping):
                    sample_audits.append(value)
                if not isinstance(sample, Mapping) or (
                    "scalar_post_root_polish_audit" not in sample
                ):
                    scalar_errors.append(
                        "accepted root omits explicit scalar post-root zero/use field"
                    )
                elif isinstance(
                    sample.get("scalar_post_root_polish_audit"),
                    Mapping,
                ):
                    sample_scalar_audits.append(sample["scalar_post_root_polish_audit"])
            if not sample_audits or record.get("surface_film_audit") != sample_audits[-1]:
                errors.append("macrostep film audit differs from its final accepted leaf")
            if isinstance(scalar_summary, Mapping):
                serialized_entries = scalar_summary.get("audits")
                serialized_audits = (
                    [
                        entry.get("audit")
                        for entry in serialized_entries
                        if isinstance(entry, Mapping)
                    ]
                    if isinstance(serialized_entries, list)
                    else []
                )
                if serialized_audits != sample_scalar_audits:
                    scalar_errors.append(
                        "accepted step scalar summary differs from its accepted-root samples"
                    )

    raw_scalar_metric = _scalar_post_root_accounting_metrics(
        raw.get("scalar_post_root_polish_accounting"),
        label="raw campaign",
        errors=scalar_errors,
    )
    if raw_scalar_metric is not None:
        (
            _,
            raw_used_steps,
            _,
            raw_used_roots,
            raw_before,
            raw_after,
            raw_work,
        ) = raw_scalar_metric
        scalar_step_use_count += raw_used_steps
        scalar_root_use_count += raw_used_roots
        scalar_work += raw_work
        if raw_before is not None:
            scalar_before_values.append(raw_before)
        if raw_after is not None:
            scalar_after_values.append(raw_after)

    accounting_partitions_present = any(
        name in raw
        for name in (
            "bootstrap",
            "shared_hot_lean_baseline_prefixes",
            "determinism",
        )
    )
    if accounting_partitions_present:
        try:
            bootstrap_mapping = raw["bootstrap"]
            prefix_mapping = raw["shared_hot_lean_baseline_prefixes"]
            determinism_mapping = raw["determinism"]
            if not (
                isinstance(bootstrap_mapping, Mapping)
                and isinstance(prefix_mapping, Mapping)
                and isinstance(determinism_mapping, Mapping)
            ):
                raise TypeError("scalar work partitions are not mappings")
            bootstrap_account = backend._aggregate_scalar_post_root_polish_summaries(  # noqa: SLF001
                tuple(item["scalar_post_root_polish"] for item in bootstrap_mapping.values())
            )
            expected_raw_account = backend._combine_scalar_post_root_polish_accounting(  # noqa: SLF001
                (
                    bootstrap_account,
                    *(
                        item["scalar_post_root_polish_accounting"]
                        for item in prefix_mapping.values()
                    ),
                    *(
                        job["incremental_campaign_scalar_post_root_polish_accounting"]
                        for job in jobs
                        # R7/R12: wall jobs carry no scalar post-root work.
                        if not backend._is_r7_wall_job(job)  # noqa: SLF001
                    ),
                    determinism_mapping["scalar_post_root_polish_accounting"],
                )
            )
        except (KeyError, TypeError, ValueError) as exc:
            scalar_errors.append(f"raw scalar post-root work partitions are incomplete: {exc}")
        else:
            if raw.get("scalar_post_root_polish_accounting") != expected_raw_account:
                scalar_errors.append(
                    "raw scalar post-root aggregate does not equal its work partitions"
                )

    solver_contract = raw.get("solver_acceptance_contract")
    solver_mapping = solver_contract if isinstance(solver_contract, Mapping) else {}
    face_conditioned_predictor = solver_mapping.get(
        failed_leaf_contract.FACE_CONDITIONED_PREDICTOR_CONTRACT_KEY
    )
    face_conditioned_predictor_exact = (
        isinstance(face_conditioned_predictor, Mapping)
        and dict(face_conditioned_predictor)
        == failed_leaf_contract.face_conditioned_predictor_contract()
    )
    scalar_policy = solver_mapping.get("scalar_post_root_polish")
    expected_scalar_policy = {
        "numerical_amendment_id": backend.SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID,
        "near_miss_multiplier_ceiling": ci.SCALAR_POST_ROOT_NEAR_MISS_MULTIPLIER,
        "residual_tolerance_unchanged": True,
        "interface_tolerance_unchanged": True,
        "physics_equations_changed": False,
        "active_or_contact_retained_cap_branch_eligible": False,
        "provisional_work_budget_step_eligible": False,
        "failed_polish_commits_partial_candidate": False,
        "exact_rollback_on_failure": True,
        "every_use_must_be_serialized": True,
        "absence_means_zero_uses": True,
    }
    if scalar_policy != expected_scalar_policy:
        scalar_errors.append("raw campaign changed or omitted scalar post-root policy")
    errors.extend(scalar_errors)
    errors.extend(surface_constitutive_errors)
    scalar_serialization_passed = not scalar_errors
    maximum_condition = solver_mapping.get("maximum_condition_proxy")
    all_full_rank = (
        bool(condition_proxies)
        and _finite(maximum_condition)
        and all(0.0 < value <= float(maximum_condition) for value in condition_proxies)
    )
    if not all_full_rank:
        errors.append("an accepted square residual lacks a certified finite Jacobian")
    rollback_preserved = bool(attempts) and all(
        attempt.get("rollback_identity_preserved") is True for attempt in attempts
    )
    if not rollback_preserved:
        errors.append("an attempted solve did not preserve exact input identity")
    # R7/R12: verified wall jobs carry no gates by construction; their
    # verification is the excluded-cell fail-closed accounting.
    gated_jobs = [
        job
        for job in jobs
        if not backend._is_r7_wall_job(job)  # noqa: SLF001
    ]
    no_clipping = raw.get("accepted_state_clipping_or_projection") is False and all(
        isinstance(job, Mapping)
        and isinstance(job.get("gates"), Mapping)
        and job["gates"].get("no_clipping_or_bound_projection") is True
        for job in gated_jobs
    )
    no_relaxation = solver_mapping.get("tolerance_relaxation_allowed") is False and all(
        isinstance(job, Mapping)
        and isinstance(job.get("gates"), Mapping)
        and job["gates"].get("no_tolerance_relaxation") is True
        for job in gated_jobs
    )
    no_equation_fallback = (
        raw.get("equations_changed_by_harness") is False
        and solver_mapping.get("seed_only_boundary_homotopy_may_change_accepted_equations") is False
    )
    declared_grids = solver_mapping.get("backward_euler_grids_s")
    if (
        not isinstance(declared_grids, list)
        or not declared_grids
        or any(value not in FORMAL_TIMESTEPS_S for value in declared_grids)
        or declared_grids != sorted(set(declared_grids), reverse=True)
    ):
        failed_leaf_errors = ["raw campaign requested macro grids are not an ordered formal subset"]
    else:
        # R7/R12: verified wall jobs are excluded-cell evidence; a campaign
        # whose EVERY job is a wall delegates entirely to that accounting.
        raw_job_list = raw.get("jobs") if isinstance(raw.get("jobs"), list) else []
        nonwall_job_list = [
            job
            for job in raw_job_list
            if not backend._is_r7_wall_job(job)  # noqa: SLF001
        ]
        if raw_job_list and not nonwall_job_list:
            failed_leaf_errors = []
        else:
            failed_leaf_errors = failed_leaf_contract.campaign_contract_errors(
                {**raw, "jobs": nonwall_job_list},
                expected_timesteps_s=tuple(declared_grids),
                segment_duration_s=FORCING_TIMESCALE_S,
            )
    errors.extend(failed_leaf_errors)
    frozen_failed_leaf_controller_passed = not failed_leaf_errors
    common_multiplier = bool(audits) and all(
        _finite(audit.get("stress_multiplier"))
        and float(audit["stress_multiplier"]) == installed.stress_multiplier
        for audit in audits
    )
    exact_a = bool(audits) and all(
        audit.get("production_used_exact_ackermann_factor") is True
        and audit.get("production_used_exact_unity_ackermann_factor") is False
        and audit.get("unity_ackermann_is_counterfactual_only") is True
        for audit in audits
    )
    exact_phy053_secant_cp = bool(audits) and all(
        audit.get("production_used_exact_phy053_partial_enthalpy_secant_heat_capacities") is True
        and audit.get("legacy_engineering_regression_heat_capacities_are_counterfactual_only")
        is True
        for audit in audits
    )
    legacy_cp_comparison_audited = bool(audits) and (len(legacy_cp_series_errors) == len(audits))
    legacy_cp_all_within_threshold = legacy_cp_comparison_audited and all(
        value <= LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
        for value in legacy_cp_series_errors
    )
    mass_guard = bool(audits) and len(departures) == len(audits)
    heat_closure = bool(audits) and len(closures) == len(audits)
    unique_closure = bool(audits) and all(
        audit.get("monotone_unique_root_diagnostic_passed") is True
        or audit.get("exact_equal_temperature_branch") is True
        for audit in audits
    )
    boundaries = (
        installed.hot_lean,
        installed.hot_rich,
        installed.cool_lean,
        installed.cool_rich,
    )
    no_external_liquid = all(boundary.liquid_saturation == 0.0 for boundary in boundaries)
    all_audited = (
        expected_audit_count > 0 and len(audits) == expected_audit_count and not film_errors
    )
    contract = {
        "accepted_residuals": {
            "all_square_and_full_rank": all_full_rank,
            "basis": (
                "unchanged square cut residual plus finite accepted solver-Jacobian "
                "condition proxy below the predeclared maximum"
            ),
            "accepted_condition_proxy_count": len(condition_proxies),
            "maximum_observed_condition_proxy": (
                max(condition_proxies) if condition_proxies else None
            ),
        },
        "rollback": {
            "every_failed_attempt_preserved_input_identity": rollback_preserved,
            "no_failed_state_committed": rollback_preserved,
        },
        "solver": {
            "no_clipping_or_bound_projection": no_clipping,
            "no_coefficient_floor": True,
            "no_penalty_source": True,
            "no_tolerance_relaxation": no_relaxation,
            "no_equation_changing_fallback": no_equation_fallback,
            "no_topology_selection_from_rejected_trial": rollback_preserved,
            "seed_only_globalization_may_change_accepted_equations": False,
            "frozen_failed_leaf_controller_contract_passed": (frozen_failed_leaf_controller_passed),
            "face_conditioned_dry_composition_predictor_contract_passed": (
                face_conditioned_predictor_exact
            ),
            "accepted_subdivision_policy": failed_leaf_contract.POLICY_ID,
            "arbitrary_accepted_subdivision_allowed": False,
        },
        "scalar_post_root_polish": {
            "numerical_amendment_id": (backend.SCALAR_POST_ROOT_NUMERICAL_AMENDMENT_ID),
            "used_step_count": scalar_step_use_count,
            "used_accepted_root_count": scalar_root_use_count,
            "serialized_job_occurrence_step_use_count": (serialized_job_scalar_step_use_count),
            "serialized_job_occurrence_root_use_count": (serialized_job_scalar_root_use_count),
            "serialized_job_occurrence_residual_evaluations": (serialized_job_scalar_work),
            "maximum_scaled_residual_before": (
                max(scalar_before_values) if scalar_before_values else None
            ),
            "maximum_scaled_residual_after": (
                max(scalar_after_values) if scalar_after_values else None
            ),
            "total_scalar_polish_residual_evaluations": scalar_work,
            "all_uses_fully_serialized": scalar_serialization_passed,
            "absence_means_zero_uses": True,
            "physics_equations_or_acceptance_tolerances_changed": False,
        },
        "film": {
            "all_accepted_roots_audited": all_audited,
            "total_accepted_root_count": expected_audit_count,
            "audited_accepted_root_count": len(audits),
            "fast_mass_surface_equals_bulk": mass_guard,
            "fast_mass_guard_enforced_at_accepted_root": mass_guard,
            "fast_mass_departure_threshold": FAST_MASS_DEPARTURE_LIMIT,
            "maximum_absolute_mass_fraction_departure": (max(departures) if departures else None),
            "production_used_exact_ackermann_factor": exact_a,
            "production_used_exact_phy053_partial_enthalpy_secant_heat_capacities": (
                exact_phy053_secant_cp
            ),
            "legacy_engineering_regression_heat_capacities_counterfactual_only": (
                exact_phy053_secant_cp
            ),
            "legacy_cp_counterfactual_comparison_audited_at_all_accepted_roots": (
                legacy_cp_comparison_audited
            ),
            "legacy_cp_series_conductance_relative_error_threshold": (
                LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
            ),
            "maximum_legacy_cp_series_conductance_relative_error": (
                max(legacy_cp_series_errors) if legacy_cp_series_errors else None
            ),
            "all_legacy_cp_counterfactuals_within_0_005": (legacy_cp_all_within_threshold),
            "legacy_cp_counterfactual_is_release_acceptance_guard": False,
            "unity_ackermann_is_counterfactual_only": exact_a,
            "unity_ackermann_shortcut_used": False,
            "a1_counterfactual_promoted_to_production": False,
            "exact_scalar_surface_closure_converged": heat_closure,
            "normalized_surface_heat_closure_limit": (NORMALIZED_SURFACE_HEAT_CLOSURE_LIMIT),
            "maximum_normalized_surface_heat_closure_residual": (
                max(closures) if closures else None
            ),
            "exact_ackermann_series_identity_verified": (
                heat_closure and exact_a and not film_errors
            ),
            "monotone_unique_surface_root_diagnostic_passed": unique_closure,
            "common_heat_mass_multiplier_preserved": common_multiplier,
            "external_liquid_saturation_zero": no_external_liquid,
        },
        "surface_constitutive_oracle": {
            "amendment_id": surface_oracle.AMENDMENT_ID,
            # R14: registered walls carry no A15 evidence; completeness is
            # judged over the non-wall job population.
            "every_accepted_root_reconstructed_from_raw_primitives": (
                bool(jobs)
                and surface_constitutive_job_count == len(gated_jobs)
                and not surface_constitutive_errors
            ),
            "accepted_job_count": surface_constitutive_job_count,
            "accepted_root_count": surface_constitutive_root_count,
            "N_h_or_N_t_substitution_allowed": False,
            "absolute_return_sign_prescribed": False,
            "physics_or_numerical_equations_changed": False,
            "errors": surface_constitutive_errors,
        },
        "amendment_18_19_evidence": amendment_18_19_evidence,
        "interface_saturation_endpoint_audit": (
            interface_saturation_endpoint_audit
        ),
        "dry_energy_capacity_stencil_audit": dry_energy_capacity_stencil_audit,
        "face_conditioned_dry_composition_predictor": (
            dict(face_conditioned_predictor)
            if isinstance(face_conditioned_predictor, Mapping)
            else None
        ),
    }
    required = (
        all_full_rank,
        rollback_preserved,
        no_clipping,
        no_relaxation,
        no_equation_fallback,
        frozen_failed_leaf_controller_passed,
        face_conditioned_predictor_exact,
        scalar_serialization_passed,
        all_audited,
        mass_guard,
        exact_a,
        exact_phy053_secant_cp,
        legacy_cp_comparison_audited,
        heat_closure,
        unique_closure,
        common_multiplier,
        no_external_liquid,
        # R7/R12: every NON-WALL job must carry exact A15 evidence; wall
        # jobs are excluded-cell evidence (a wall-only campaign delegates).
        bool(jobs)
        and surface_constitutive_job_count == len(gated_jobs)
        and not surface_constitutive_errors,
        amendment_18_19_evidence["contract_passed"] is True,
    )
    contract["contract_passed"] = all(required) and not errors
    contract["errors"] = errors
    return contract


def _coefficient_case(identifier: str) -> p1ef.EngineeringCoefficientCase:
    cases = {case.identifier: case for case in p1ef.engineering_coefficient_cases()}
    try:
        return cases[identifier]
    except KeyError as exc:
        raise ValueError(
            f"unknown coefficient case {identifier!r}; choose one of {tuple(cases)}"
        ) from exc


def _film_boundary(
    temperature_k: float,
    y_hexane: float,
    stress_multiplier: float,
    label: str,
):
    return p1ef.make_reduced_film_boundary(
        bulk_temperature_k=temperature_k,
        bulk_y_hexane=y_hexane,
        stress_multiplier=stress_multiplier,
        forcing_timescale_s=FORCING_TIMESCALE_S,
        label=label,
    )


@_rollback_backend_case_state_on_failure
def install_engineering_case(
    coefficient_case_id: str,
    *,
    stress_multiplier: float,
    structural_conductivity_w_m_k: float,
    pulse_temperature_k: float | None = None,
    pulse_y_hexane: float | None = None,
    enforce_p1ef_smooth_wet_retained_water_guard: bool = False,
    vanishing_wet_cut_threshold_authority: (
        cut.VanishingWetCutThresholdAuthority
    ) = cut.P1EF_PRODUCTION_VANISHING_WET_CUT_AUTHORITY,
) -> InstalledEngineeringCase:
    """Persist one valid P1EF case process-locally; roll back failed installs."""

    selected_pulse_temperature_k = (
        PULSE_TEMPERATURE_K
        if pulse_temperature_k is None
        else pulse_temperature_k
    )
    if (
        isinstance(selected_pulse_temperature_k, bool)
        or not isinstance(selected_pulse_temperature_k, (int, float))
        or not math.isfinite(float(selected_pulse_temperature_k))
        or float(selected_pulse_temperature_k) <= 0.0
    ):
        raise ValueError("pulse temperature must be positive and finite")
    selected_pulse_temperature_k = float(selected_pulse_temperature_k)
    selected_pulse_y_hexane = (
        PULSE_Y_HEXANE if pulse_y_hexane is None else pulse_y_hexane
    )
    if (
        isinstance(selected_pulse_y_hexane, bool)
        or not isinstance(selected_pulse_y_hexane, (int, float))
        or not math.isfinite(float(selected_pulse_y_hexane))
        or not BASELINE_Y_HEXANE < float(selected_pulse_y_hexane) < 1.0
    ):
        raise ValueError(
            "pulse n-hexane mole fraction must be finite and lie strictly "
            f"between the baseline {BASELINE_Y_HEXANE} and one"
        )
    selected_pulse_y_hexane = float(selected_pulse_y_hexane)
    if not isinstance(enforce_p1ef_smooth_wet_retained_water_guard, bool):
        raise ValueError("P1EF smooth retained-water guard policy must be boolean")
    if vanishing_wet_cut_threshold_authority not in {
        cut.P1EF_PRODUCTION_VANISHING_WET_CUT_AUTHORITY,
        cut.P1EF_SENSITIVITY_VANISHING_WET_CUT_AUTHORITY,
    }:
        raise ValueError(
            "P1EF wet-cut treatment needs the frozen production or sensitivity authority"
        )

    coefficient_case = _coefficient_case(coefficient_case_id)
    configuration = p1ef.build_engineering_particle_configuration(
        coefficient_case,
        structural_conductivity_w_m_k=structural_conductivity_w_m_k,
    )
    configuration = replace(
        configuration,
        cut=replace(
            configuration.cut,
            vanishing_wet_cut_threshold_authority=(
                vanishing_wet_cut_threshold_authority
            ),
        ),
    )
    hot_lean = _film_boundary(
        BASELINE_TEMPERATURE_K,
        BASELINE_Y_HEXANE,
        stress_multiplier,
        "P1EF hot/lean baseline and exact return",
    )
    hot_rich = _film_boundary(
        BASELINE_TEMPERATURE_K,
        selected_pulse_y_hexane,
        stress_multiplier,
        "P1EF composition-only hot/rich finite pulse",
    )
    cool_lean = _film_boundary(
        selected_pulse_temperature_k,
        BASELINE_Y_HEXANE,
        stress_multiplier,
        "P1EF temperature-only cool/lean finite pulse",
    )
    cool_rich = _film_boundary(
        selected_pulse_temperature_k,
        selected_pulse_y_hexane,
        stress_multiplier,
        "P1EF combined cool/rich finite pulse",
    )

    backend.PRESSURE_PA = configuration.process_pressure_pa
    backend.RADIUS_M = configuration.radius_m
    # GT-PS-2-P1E-R11/R16: the full A15 reversal demonstration gates at the
    # NOMINAL operating point (nominal film stress AND nominal structural
    # conductivity); coefficient/stress corners gate on the measured
    # invariant.  Measured basis for the conductivity clause: at k=0.29 the
    # combined pulse-end flux stays outward (+1.44e-7 after an 87 percent
    # decay) — the crossing is a property of the nominal point.
    backend.A15_FULL_REVERSAL_GATING = (
        stress_multiplier == 1.0
        and structural_conductivity_w_m_k
        == NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K
    )
    backend.DRY_CONFIG = configuration.dry
    backend.WET_MODEL = configuration.wet
    backend.CONFIG = configuration.cut
    backend.RETAINED_MOBILITY = configuration.dry.retained_water_mobility
    # Amendment 04 supersedes the former 0.1 s passing campaign.  The mature
    # backend is installed process-locally with the prospective symmetric
    # segment and requested macro grids; the former trace remains a separate
    # adverse challenge and is never inserted into this job matrix.
    backend.SEGMENT_DURATION_S = FORCING_TIMESCALE_S
    backend.DEFAULT_TIMESTEPS_S = FORMAL_TIMESTEPS_S
    backend.CONTROLS = replace(
        backend.CONTROLS,
        wet_temperature_bounds_k=WET_TEMPERATURE_BOUNDS_K,
        # Use the already frozen Luikov/GAB capacity authority as the open
        # numerical chart endpoint.  The endpoint itself remains excluded by
        # the transformed solver coordinate; an earlier rounded 0.235 ceiling
        # sat 6.70904e-4 below this authority and became an artificial active
        # bound during the exact-A cool/rich pulse tail.
        wet_water_bounds=(0.0231, configuration.wet.luikov.W_cap),
        interface_temperature_bounds_k=INTERFACE_TEMPERATURE_BOUNDS_K,
        # GT-PS-2-P1E-21 section 3.2: derive the gate, never restate its
        # literal.  The formal campaign runs the same authority the solver
        # declares, so evidence and runtime cannot disagree about the gate.
        nonlinear_residual_tolerance=ci.CutSolverControls.nonlinear_residual_tolerance,
        ledger_tolerance=1.0e-10,
        nonlinear_step_tolerance=1.0e-12,
        maximum_function_evaluations=1200,
        maximum_condition_proxy=1.0e12,
        enforce_p1ef_smooth_wet_retained_water_guard=(
            enforce_p1ef_smooth_wet_retained_water_guard
        ),
    )
    backend.HOT_LEAN_BOUNDARY = hot_lean
    backend.HOT_RICH_BOUNDARY = hot_rich
    backend.COOL_LEAN_BOUNDARY = cool_lean
    backend.COOL_RICH_BOUNDARY = cool_rich
    backend.REVERSAL_TEMPERATURE_K = selected_pulse_temperature_k
    backend.COLDEST_REVERSAL_FACE_TEMPERATURE_K = 0.5 * (
        configuration.dry.conditioned_temperature_domain.solver_bounds_k[0]
        + selected_pulse_temperature_k
    )
    backend.COLDEST_COMPOSITION_FACE_TEMPERATURE_K = 0.5 * (
        configuration.dry.conditioned_temperature_domain.solver_bounds_k[0] + BASELINE_TEMPERATURE_K
    )
    backend.REVERSAL_Y_HEXANE = selected_pulse_y_hexane
    backend.COMPOSITION_REVERSAL_Y_HEXANE = selected_pulse_y_hexane
    backend.REVERSAL_EXACT_CHART_COORDINATE = cp.encode_gas_only_y(
        backend.COLDEST_REVERSAL_FACE_TEMPERATURE_K,
        configuration.process_pressure_pa,
        selected_pulse_y_hexane,
        configuration.dry.pore,
    )

    scenarios = {
        backend.SCENARIO_COMPOSITION_ONLY: backend.StressScenario(
            backend.SCENARIO_COMPOSITION_ONLY,
            "composition_only_at_fixed_boundary_temperature",
            "hot_rich_composition_only_finite_pulse",
            hot_rich,
            "P1EF independent binary counter-diffusion and conserved-flux response",
        ),
        backend.SCENARIO_TEMPERATURE_ONLY: backend.StressScenario(
            backend.SCENARIO_TEMPERATURE_ONLY,
            "temperature_only_at_fixed_boundary_composition",
            "cool_lean_temperature_only_finite_pulse",
            cool_lean,
            "P1EF resolved conduction and local-equilibrium storage response",
        ),
        backend.SCENARIO_COMBINED: backend.StressScenario(
            backend.SCENARIO_COMBINED,
            "correlated_temperature_and_composition",
            "cool_rich_finite_pulse",
            cool_rich,
            "P1EF simultaneous bounded step-change stress",
        ),
    }
    backend.SCENARIOS = scenarios
    backend.COMBINED_SCENARIO = scenarios[backend.SCENARIO_COMBINED]
    backend.DYNAMIC_BOUNDARY_SEGMENTS = backend.COMBINED_SCENARIO.segments
    backend.FIXED_RADIUS_FORCING_HISTORY_ID = FORMAL_HISTORY_ID

    _validate_installed_trace(configuration, (hot_lean, hot_rich, cool_lean, cool_rich))
    return InstalledEngineeringCase(
        coefficient_case=coefficient_case,
        configuration=configuration,
        stress_multiplier=stress_multiplier,
        pulse_temperature_k=selected_pulse_temperature_k,
        pulse_y_hexane=selected_pulse_y_hexane,
        p1ef_smooth_retained_water_guard_enabled=(
            enforce_p1ef_smooth_wet_retained_water_guard
        ),
        vanishing_wet_cut_threshold_authority=(
            vanishing_wet_cut_threshold_authority
        ),
        hot_lean=hot_lean,
        hot_rich=hot_rich,
        cool_lean=cool_lean,
        cool_rich=cool_rich,
    )


def _validate_installed_trace(
    configuration: p1ef.EngineeringParticleConfiguration,
    boundaries: Sequence[object],
) -> None:
    """Fail before solving if a predeclared trace leaves the dry topology."""

    for boundary in boundaries:
        temperatures = (
            boundary.temperature_k,
            0.5 * (INITIAL_DRY_TEMPERATURE_K + boundary.temperature_k),
        )
        for temperature_k in temperatures:
            interval = cp.gas_only_composition_interval(
                temperature_k,
                configuration.process_pressure_pa,
                configuration.dry.pore,
            )
            if not interval.lower_y_hexane < boundary.y_hexane < interval.upper_y_hexane:
                raise ValueError(
                    "P1EF boundary trace is outside the exact gas-only topology; "
                    "do not clip the engineering stress"
                )


def solve_target_n2_root(
    installed: InstalledEngineeringCase,
    seed_candidate: cut.CutTransportUnknowns,
) -> tuple[ci.CutIntegratorStep | None, dict]:
    """Solve the exact target first step or return a fail-closed audit."""

    before = backend.initial_integrator(2, installed.configuration.radius_m)
    continuation = cc.rescale_same_before_increment(
        before,
        seed_candidate,
        surface_boundary=installed.hot_lean,
        source_dt_s=REFERENCE_N2_SEED_TIMESTEP_S,
        target_dt_s=FORCING_TIMESCALE_S,
        stefan_flux_scale_floor_mol_m2_s=backend.STEFAN_SCALE_FLOOR,
        label=(
            "P1EF 0.1 s accepted reference root rescaled only as a nonlinear "
            "seed for the independent 0.075 s target equations"
        ),
    )
    step, record, error = backend._try_step_detailed(  # noqa: SLF001
        before,
        FORCING_TIMESCALE_S,
        installed.hot_lean,
        continuation.seed,
        role="P1EF exact target N=2 root",
        predictor_fraction=1.0,
    )
    audit = {
        "accepted": step is not None,
        "attempt": asdict(record),
        "rollback_identity_preserved": (True if error is None else error.rollback_state is before),
        "failed_state_committed": False,
        "error_type": None if error is None else type(error).__name__,
        "error": None if error is None else str(error),
        "seed_only_time_continuation": asdict(continuation.record),
        "reference_seed_accepted_state_or_inventory_reused": False,
        "target_timestep_s": FORCING_TIMESCALE_S,
        "scalar_post_root_polish_audit": None,
        "interface_saturation_endpoint_audit": None,
        "dry_energy_capacity_stencil_audit": None,
        "amendment_18_19_evidence_contract": (
            _failed_target_amendment_18_19_evidence(
                ("target N=2 root was not accepted",)
            )
        ),
        "amendment_18_19_evidence_contract_passed": False,
    }
    if step is not None:
        try:
            balance_roots = backend._accepted_balance_interval_roots(step)  # noqa: SLF001
            if len(balance_roots) != 1:
                raise RuntimeError(
                    "exact target N=2 solve did not produce exactly one balance root"
                )
            evidence_sample = backend._surface_flux_sample(  # noqa: SLF001
                balance_roots[0],
                include_constitutive_evidence=False,
            )
            target_amendment_evidence = _accepted_root_amendment_18_19_evidence(
                evidence_sample
            )
        except (KeyError, TypeError, ValueError, RuntimeError) as exc:
            evidence_error = f"target N=2 A18/A19 evidence failed: {exc}"
            audit["accepted"] = False
            audit["error_type"] = type(exc).__name__
            audit["error"] = evidence_error
            audit["amendment_18_19_evidence_contract"] = (
                _failed_target_amendment_18_19_evidence((evidence_error,))
            )
            return None, audit
        audit["interface_saturation_endpoint_audit"] = evidence_sample[
            "interface_saturation_endpoint_audit"
        ]
        audit["dry_energy_capacity_stencil_audit"] = evidence_sample[
            "dry_energy_capacity_stencil_audit"
        ]
        audit["amendment_18_19_evidence_contract"] = target_amendment_evidence
        audit["amendment_18_19_evidence_contract_passed"] = True
        audit["accepted_root"] = asdict(step.assembly.candidate)
        audit["maximum_scaled_residual"] = step.ledger.maximum_scaled_residual
        audit["condition_proxy"] = step.ledger.condition_proxy
        audit["maximum_step_ledger_residual"] = step.ledger.maximum_step_ledger_residual
        audit["maximum_cumulative_ledger_residual"] = step.ledger.maximum_cumulative_ledger_residual
        audit["surface_film_audit"] = (
            None if step.surface_film_audit is None else asdict(step.surface_film_audit)
        )
        audit["scalar_post_root_polish_audit"] = (
            None
            if step.ledger.scalar_post_root_polish_audit is None
            else asdict(step.ledger.scalar_post_root_polish_audit)
        )
        audit["wet_retained_cap_active_piece_count"] = (
            step.ledger.wet_retained_cap_active_piece_count
        )
        audit["wet_retained_cap_maximum_loading"] = (
            step.ledger.wet_retained_cap_maximum_loading
        )
        audit["wet_retained_cap_maximum_capacity_dual_over_rt"] = (
            step.ledger.wet_retained_cap_maximum_capacity_dual_over_rt
        )
        audit["wet_retained_cap_maximum_complementarity_product"] = (
            step.ledger.wet_retained_cap_maximum_complementarity_product
        )
        audit["wet_retained_cap_exact_graph_without_tolerance"] = (
            step.ledger.wet_retained_cap_exact_graph_without_tolerance
        )
        audit["vanishing_wet_cut_treatment"] = asdict(
            step.assembly.vanishing_wet_cut_treatment_audit
        )
        audit["vanishing_wet_cut_activation_count"] = (
            step.assembly.vanishing_wet_cut_treatment_audit.activation_count
        )
        audit["maximum_vanishing_wet_cut_water_pair_sum_difference_mol_s"] = abs(
            step.assembly.vanishing_wet_cut_treatment_audit
            .water_pair_sum_difference_mol_s
        )
        audit[
            "maximum_vanishing_wet_cut_common_energy_pair_sum_difference_w"
        ] = abs(
            step.assembly.vanishing_wet_cut_treatment_audit
            .common_energy_pair_sum_difference_w
        )
        wet_retained_water_applicability = (
            ci.certify_explicit_wet_retained_water_applicability(
                step.after.controls,
                step.after.transport.wet_retained_water_loadings,
                (
                    step.after.transport
                    .effective_wet_retained_water_capacity_duals_over_rt
                ),
                step.after.transport.config.wet.luikov,
            )
        )
        audit["wet_retained_water_applicability"] = asdict(
            wet_retained_water_applicability
        )
        audit["wet_retained_water_guard_loading"] = (
            wet_retained_water_applicability.guard_loading
        )
        audit["wet_retained_water_guard_margin"] = math.fsum(
            (
                wet_retained_water_applicability.guard_loading,
                -step.ledger.wet_retained_cap_maximum_loading,
            )
        )
        audit["wet_retained_water_capacity_loading"] = (
            installed.configuration.wet.luikov.W_cap
        )
        audit["wet_retained_water_capacity_margin"] = math.fsum(
            (
                installed.configuration.wet.luikov.W_cap,
                -step.ledger.wet_retained_cap_maximum_loading,
            )
        )
    return step, audit


# ---------------------------------------------------------------------------
# Brief B-F2 (2026-09-30; paper_gaps PLAN.md item 5): the return-hold knob for
# a longer ablation window.  The three 0.075 s segments stay exactly as they
# are.  With ``DTDC_ABLATION_RETURN_HOLD_S`` set to a positive number of
# seconds, every completed three-segment trajectory of ``run_case`` is
# continued past its third segment at the same stride (its own dt), holding the
# return forcing (the scenario's exact return boundary object, the hot/lean
# baseline) for that many seconds, through the backend's own step functions
# (the history-routed face route, the same-cell recovery, the uncertified-face
# recovery and the frozen failed-leaf macrostep, exactly as the return
# segment's steps after its first), its ledgers and its step record.  The hold
# is reported in a separate ``return_hold`` block of the job; no three-segment
# gate, contract, observable or per-step row reads it, and the caller's
# accepted-step list and audit sink stay the three-segment ones.  A refused hold
# step ends the hold, typed and uncommitted, and leaves the three-segment job as
# it was.  A held step is committed only when its frozen failed-leaf controller
# contract passes, the per-step test of the three-segment gate; a held step
# whose contract fails (for example a failed-leaf subdivision on a step that does
# not start after an imposed boundary discontinuity) refuses with
# ``ReturnHoldFailedLeafContractError`` and ends the hold there.  Unset (or 0) is
# today's behaviour: nothing is installed and nothing is added to any record.
# ---------------------------------------------------------------------------

ABLATION_RETURN_HOLD_ENVIRONMENT_VARIABLE = "DTDC_ABLATION_RETURN_HOLD_S"
RETURN_HOLD_DEFINITION = {
    "record": "B-F2 2026-09-30 return-hold knob for the longer ablation window",
    "environment_variable": ABLATION_RETURN_HOLD_ENVIRONMENT_VARIABLE,
    "boundary": "the scenario's exact return boundary object (segments[2][1])",
    "stride": "the job's own requested dt; the hold must be a whole number of strides",
    "step_functions": (
        "backend return-segment body for a step after the first: history-routed "
        "face route, same-cell recovery with the current-state anchor, "
        "uncertified-face recovery, frozen failed-leaf macrostep (is_reversal false)"
    ),
    "held_step_commit_requires": (
        "frozen_failed_leaf_controller_contract.passed is True, the per-step test of "
        "the three-segment gate frozen_failed_leaf_controller_contract_passed; "
        "otherwise ReturnHoldFailedLeafContractError: the step uncommitted, the hold ended"
    ),
    "read_by_three_segment_gates": False,
    "counts_toward_release_or_campaign_gates": False,
    "physically_qualifying": False,
}


class ReturnHoldKnobError(ValueError):
    """``DTDC_ABLATION_RETURN_HOLD_S`` cannot be honoured as declared."""


def ablation_return_hold_s() -> float | None:
    """Return the return-hold duration in seconds, or None when the knob is off.

    Unset, ``"0"`` or any spelling of zero is off.  A value that is not a finite
    non-negative number refuses, so a mistyped knob never silently selects a path.
    """

    raw = os.environ.get(ABLATION_RETURN_HOLD_ENVIRONMENT_VARIABLE)
    if raw is None:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise ReturnHoldKnobError(
            f"{ABLATION_RETURN_HOLD_ENVIRONMENT_VARIABLE} is not a number: {raw!r}"
        ) from exc
    if not math.isfinite(value) or value < 0.0:
        raise ReturnHoldKnobError(
            f"{ABLATION_RETURN_HOLD_ENVIRONMENT_VARIABLE} must be finite and non-negative, "
            f"got {raw!r}"
        )
    return None if value == 0.0 else value


class ReturnHoldFailedLeafContractError(RuntimeError):
    """A held step failed the frozen failed-leaf controller contract; it is not committed.

    ``details`` names the step (hold index, times, recovery mode), the contract
    terms that are not true and the contract block as the backend recorded it.
    """

    def __init__(self, message: str, details: dict[str, Any]) -> None:
        super().__init__(message)
        self.details = details


# The boolean terms of the backend's frozen failed-leaf controller predicate as
# ``_step_record`` writes them (``frozen_failed_leaf_controller_contract``);
# reported by name when a held step's contract fails.  The decision itself is
# the recorded ``passed`` flag, exactly as the three-segment gate reads it.
FROZEN_FAILED_LEAF_CONTRACT_TERMS = (
    "licensed_only_after_imposed_boundary_discontinuity",
    "accepted_binary_depth_bounded",
    "leaf_grid_disclosed_and_bounded",
    "fixed_work_budget_observed",
    "attempt_paths_form_failed_leaf_binary_tree",
    "initial_requested_macro_failure_disclosed",
    "one_returned_attempt_per_accepted_leaf",
    "per_leaf_seed_attempt_budget_observed",
    "raw_macro_ledger_le_1e_10",
    "atomic_macro_input_identity_preserved",
    "immutable_leaf_chain_identity_preserved",
    "exact_post_jump_boundary_identity_preserved",
    "every_attempt_rollback_identity_preserved",
)


def require_held_step_failed_leaf_contract(record: Mapping[str, Any], hold_index: int) -> None:
    """Refuse a held step whose frozen failed-leaf controller contract did not pass.

    The predicate is the three-segment gate's own per-step test
    (``_accepted_grid_contract``: the contract is a dict and its ``passed`` is
    ``True``), applied before the held step is committed instead of after the
    trajectory, so no held step past a failing one is ever taken.
    """

    contract = record.get("frozen_failed_leaf_controller_contract")
    if isinstance(contract, dict) and contract.get("passed") is True:
        return
    failing_terms = (
        [term for term in FROZEN_FAILED_LEAF_CONTRACT_TERMS if contract.get(term) is not True]
        if isinstance(contract, dict)
        else ["frozen_failed_leaf_controller_contract missing"]
    )
    details = {
        "hold_step_index": hold_index,
        "stage": record.get("stage"),
        "time_before_s": record.get("time_before_s"),
        "time_after_s": record.get("time_after_s"),
        "recovery_mode": record.get("recovery_mode"),
        "adaptive_subdivision_used": record.get("adaptive_subdivision_used"),
        "adaptive_binary_depth": record.get("adaptive_binary_depth"),
        "starts_after_imposed_boundary_discontinuity": record.get(
            "starts_after_imposed_boundary_discontinuity"
        ),
        "maximum_scaled_residual": record.get("maximum_scaled_residual"),
        "maximum_step_ledger_residual": record.get("maximum_step_ledger_residual"),
        "maximum_cumulative_ledger_residual": record.get("maximum_cumulative_ledger_residual"),
        "wall_time_s": record.get("wall_time_s"),
        "failing_contract_terms": failing_terms,
        "frozen_failed_leaf_controller_contract": (
            dict(contract) if isinstance(contract, dict) else None
        ),
    }
    raise ReturnHoldFailedLeafContractError(
        f"held step {hold_index} ({details['time_before_s']!r} to "
        f"{details['time_after_s']!r} s, recovery mode {details['recovery_mode']!r}) "
        "failed the frozen failed-leaf controller contract (not true: "
        f"{', '.join(failing_terms) or 'passed'}); the step is not committed and the "
        "hold ends here, as the three-segment gate "
        "frozen_failed_leaf_controller_contract_passed refuses an ordinary step",
        details,
    )


def return_hold_step_count(hold_s: float, dt_s: float) -> int:
    """Whole strides in the hold; anything else refuses and names the neighbours."""

    ratio = hold_s / dt_s
    steps = round(ratio)
    if steps >= 1 and math.isclose(steps * dt_s, hold_s, rel_tol=1.0e-12, abs_tol=0.0):
        return steps
    lower = max(1, math.floor(ratio))
    raise ReturnHoldKnobError(
        f"{ABLATION_RETURN_HOLD_ENVIRONMENT_VARIABLE}={hold_s!r} s is not a whole number "
        f"of dt={dt_s!r} s strides (nearest admissible: {lower * dt_s!r} s or "
        f"{(lower + 1) * dt_s!r} s)"
    )


def _advance_return_hold_step(previous: Any, dt_s: float, boundary: Any) -> tuple:
    """One held step: the backend's return-segment body for a step after the first."""

    try:
        routed = backend._advance_bracketed_face_from_accepted_history(  # noqa: SLF001
            previous,
            dt_s,
            boundary,
        )
        if routed is not None:
            return routed
        return backend._advance_with_recovery(  # noqa: SLF001
            previous,
            dt_s,
            boundary,
            boundary,
            is_reversal=False,
            try_current_state_anchor=isinstance(boundary, backend.ct.ReducedFilmPoreBoundary),
            provisional_maximum_function_evaluations=None,
            maximum_nonlinear_seed_attempts=None,
        )
    except backend.SameCellRecoveryError as exc:
        if isinstance(exc, backend.AcceptedHistoryFaceRouteRecoveryError):
            return backend._recover_after_uncertified_face_route(  # noqa: SLF001
                previous,
                dt_s,
                boundary,
                boundary,
                is_reversal=False,
                failure=exc,
            )
        return backend._advance_with_adaptive_same_cell_macrostep(  # noqa: SLF001
            previous,
            dt_s,
            boundary,
            boundary,
            is_reversal=False,
            initial_attempts=exc.attempts,
            initial_step_error=exc.last_step_error,
            initial_failure_audit=exc.failure_jacobian_audit,
        )


def _return_hold_state_summary(step: Any) -> dict[str, Any]:
    state = step.after
    inventory = ci.inventory_snapshot(state.transport)
    water = ci.water_phase_inventory_snapshot(state.transport)
    front_z = state.transport.geometry.front.z
    return {
        "time_s": state.transport.time_s,
        "front_z": front_z,
        "front_position_over_radius": front_z ** (1.0 / 3.0),
        "interface_temperature_k": state.last_interface_temperature_k,
        "water_inventory_mol": inventory.total_water_mol,
        "retained_water_inventory_mol": water.retained_water_mol,
        "hexane_inventory_mol": inventory.total_hexane_mol,
        "energy_inventory_j": inventory.total_energy_j,
        "cumulative_boundary_hexane_out_mol": state.corrected_cumulative_boundary_hexane_out_mol,
        "cumulative_boundary_energy_out_j": state.corrected_cumulative_boundary_energy_out_j,
    }


def run_return_hold(
    last_step: Any,
    scenario: Any,
    cells: int,
    dt_s: float,
    hold_s: float,
) -> dict[str, Any]:
    """Continue a committed three-segment trajectory at its exact return boundary."""

    steps = return_hold_step_count(hold_s, dt_s)
    return_stage, boundary = scenario.segments[2]
    start = _return_hold_state_summary(last_step)
    records: list[dict[str, Any]] = []
    refusal: dict[str, Any] | None = None
    previous = last_step
    for hold_index in range(steps):
        try:
            current, attempts, mode = _advance_return_hold_step(previous, dt_s, boundary)
            record = backend._step_record(  # noqa: SLF001
                current,
                attempts,
                return_stage,
                mode,
                starts_after_boundary_discontinuity=False,
            )
            require_held_step_failed_leaf_contract(record, hold_index)
        except Exception as exc:  # noqa: BLE001 - a refused held step ends the hold uncommitted
            try:
                evidence: Any = backend._campaign_failure_result(  # noqa: SLF001
                    scenario,
                    cells,
                    dt_s,
                    exc,
                )
            except Exception as evidence_exc:  # noqa: BLE001 - keep the refusal itself
                evidence = {"error_type": type(evidence_exc).__name__, "error": str(evidence_exc)}
            refusal = {
                "hold_step_index": hold_index,
                "time_before_s": previous.after.transport.time_s,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "failed_state_was_committed": False,
                "failure_evidence": evidence,
            }
            if isinstance(exc, ReturnHoldFailedLeafContractError):
                refusal["refused_step"] = exc.details
            break
        record["return_hold_step_index"] = hold_index
        records.append(record)
        previous = current
    residuals = [record["maximum_scaled_residual"] for record in records]
    step_ledgers = [record["maximum_step_ledger_residual"] for record in records]
    cumulative_ledgers = [record["maximum_cumulative_ledger_residual"] for record in records]
    return {
        "definition": RETURN_HOLD_DEFINITION,
        "ablation_return_hold_s": hold_s,
        "dt_s": dt_s,
        "requested_steps": steps,
        "accepted_steps": len(records),
        "completed": refusal is None,
        "refusal": refusal,
        "stage": return_stage,
        "boundary_is_scenario_return_object": boundary is scenario.segments[2][1],
        "boundary_is_scenario_baseline_object": boundary is scenario.segments[0][1],
        "start": start,
        "final": _return_hold_state_summary(previous),
        "maximum_scaled_residual": max(residuals) if residuals else None,
        "maximum_step_ledger_residual": max(step_ledgers) if step_ledgers else None,
        "maximum_cumulative_ledger_residual": (
            max(cumulative_ledgers) if cumulative_ledgers else None
        ),
        "residuals_within_backend_limit": all(
            value <= backend.NONLINEAR_RESIDUAL_LIMIT for value in residuals
        ),
        "ledgers_within_conservation_limit": all(
            value <= backend.CONSERVATION_LIMIT for value in (*step_ledgers, *cumulative_ledgers)
        ),
        "strictly_inside_all_open_solver_bounds": all(
            record["minimum_fractional_distance_to_bound"] > 0.0 for record in records
        ),
        "solver_controls_unchanged": all(
            record["solver_controls_match_campaign_contract"] for record in records
        ),
        "frozen_failed_leaf_controller_contract_passed_every_step": all(
            record["frozen_failed_leaf_controller_contract"]["passed"] is True for record in records
        ),
        "interior_face_events": sum(
            1 for record in records if record.get("interior_face_event_used") is True
        ),
        "wall_time_s": math.fsum(record["wall_time_s"] for record in records),
        "step_records": records,
        "read_by_three_segment_gates": False,
        "physically_qualifying": False,
    }


def _return_hold_trajectory(original: Callable[..., dict], hold_s: float) -> Callable[..., dict]:
    """Wrap the backend's trajectory so a completed job is continued by the hold."""

    signature = inspect.signature(original)

    @wraps(original)
    def wrapped(*args: Any, **kwargs: Any) -> dict:
        bound = signature.bind(*args, **kwargs)
        sink = bound.arguments.get("accepted_steps_out")
        if sink is None:
            sink = []
            bound.arguments["accepted_steps_out"] = sink
        job = original(*bound.args, **bound.kwargs)
        scenario = bound.arguments.get("scenario") or backend.COMBINED_SCENARIO
        last_step = sink[-1]
        if last_step.after.transport.time_s != job["final"]["time_s"]:
            raise RuntimeError("return hold lost the committed three-segment final state")
        job["return_hold"] = run_return_hold(
            last_step,
            scenario,
            bound.arguments["cells"],
            bound.arguments["dt_s"],
            hold_s,
        )
        return job

    return wrapped


@contextmanager
def _return_hold_installed(hold_s: float | None) -> Iterator[None]:
    """Install the hold around one campaign; with the knob off, install nothing."""

    if hold_s is None:
        yield
        return
    original = backend._run_trajectory  # noqa: SLF001
    backend._run_trajectory = _return_hold_trajectory(original, hold_s)  # noqa: SLF001
    try:
        yield
    finally:
        backend._run_trajectory = original  # noqa: SLF001


@_restore_backend_case_state_after_call
def run_case(
    *,
    coefficient_case_id: str,
    stress_multiplier: float,
    structural_conductivity_w_m_k: float,
    meshes: Sequence[int],
    timesteps_s: Sequence[float],
    scenarios: Sequence[str],
    pulse_temperature_k: float | None = None,
    pulse_y_hexane: float | None = None,
    enforce_p1ef_smooth_wet_retained_water_guard: bool = False,
    vanishing_wet_cut_threshold_authority: (
        cut.VanishingWetCutThresholdAuthority
    ) = cut.P1EF_PRODUCTION_VANISHING_WET_CUT_AUTHORITY,
) -> dict:
    """Run one process-isolated case and restore the caller's backend state."""

    # B-F1: the face-tangent knobs in force for this whole case (a malformed
    # value refuses here, before any solve).
    face_tangent_seed_armed = face_tangent.face_tangent_inward_seed_armed()
    face_tangent_bound = face_tangent.face_tangent_residual_bound()
    face_tangent_bound_declared = face_tangent.face_tangent_residual_bound_declared()
    # B-F2: the return-hold knob refuses here, before any solve, when it is
    # malformed or not a whole number of any requested stride.
    return_hold_s = ablation_return_hold_s()
    if return_hold_s is not None:
        for dt_s in timesteps_s:
            return_hold_step_count(return_hold_s, float(dt_s))
    provenance_before = source_provenance()
    installed = install_engineering_case(
        coefficient_case_id,
        stress_multiplier=stress_multiplier,
        structural_conductivity_w_m_k=structural_conductivity_w_m_k,
        pulse_temperature_k=pulse_temperature_k,
        pulse_y_hexane=pulse_y_hexane,
        enforce_p1ef_smooth_wet_retained_water_guard=(
            enforce_p1ef_smooth_wet_retained_water_guard
        ),
        vanishing_wet_cut_threshold_authority=(
            vanishing_wet_cut_threshold_authority
        ),
    )
    initial_condition = build_initial_condition_contract(installed, meshes)
    first, root_audit = solve_target_n2_root(installed, P1EF_REFERENCE_N2_SEED)
    common = {
        "schema_version": SCHEMA_VERSION,
        "scope": ENGINEERING_SCOPE,
        "decision_id": p1ef.P1EF_DECISION_ID,
        "physically_qualifying": False,
        "plant_predictive": False,
        "normative_gate_status_unchanged": {
            "F2_solvent_particle": False,
            "F3_water_interface": False,
        },
        "process_isolation_required": True,
        "thread_safe": False,
        "coefficient_case": {
            "identifier": installed.coefficient_case.identifier,
            "binary_diffusivity_m2_s": (installed.coefficient_case.binary_diffusivity_m2_s),
            "apparent_water_diffusivity_m2_s": (
                installed.coefficient_case.apparent_water_diffusivity_m2_s
            ),
            "retained_water_mobility_mol_m_s": (
                installed.coefficient_case.partition.retained_water_mobility_mol_m_s
            ),
        },
        "structural_conductivity_w_m_k": structural_conductivity_w_m_k,
        "film_common_stress_multiplier": stress_multiplier,
        "radius_m": installed.configuration.radius_m,
        "pressure_pa": installed.configuration.process_pressure_pa,
        "forcing": {
            "history_id": FORMAL_HISTORY_ID,
            "selection_status": FORMAL_HISTORY_SELECTION_STATUS,
            "baseline": [BASELINE_TEMPERATURE_K, BASELINE_Y_HEXANE],
            "pulse": [installed.pulse_temperature_k, installed.pulse_y_hexane],
            "composition_only_pulse": [
                BASELINE_TEMPERATURE_K,
                installed.pulse_y_hexane,
            ],
            "temperature_only_pulse": [
                installed.pulse_temperature_k,
                BASELINE_Y_HEXANE,
            ],
            "combined_pulse": [
                installed.pulse_temperature_k,
                installed.pulse_y_hexane,
            ],
            "return_reuses_exact_baseline_object": all(
                scenario.segments[0][1] is installed.hot_lean
                and scenario.segments[2][1] is installed.hot_lean
                for scenario in backend.SCENARIOS.values()
            ),
            "segment_duration_s": FORCING_TIMESCALE_S,
            "requested_macro_timesteps_s": list(FORMAL_TIMESTEPS_S),
            "former_0p1_history_counts_as_passing_job": False,
        },
        "initial_condition": initial_condition,
        "retained_water_applicability_guard": {
            "policy_enabled": (
                installed.p1ef_smooth_retained_water_guard_enabled
            ),
            "guard_loading": (
                ci.P1EF_SMOOTH_WET_RETAINED_WATER_GUARD_LOADING
            ),
            "exact_loading_and_zero_dual_comparisons": True,
            "terminal_on_failure_without_retry_or_bisection": True,
            "capacity_or_phase_boundary_claimed": False,
        },
        "vanishing_wet_cut_treatment": {
            "amendment_id": "GT-PS-2-P1E-12",
            "threshold_authority": (
                installed.vanishing_wet_cut_threshold_authority.value
            ),
            "target_wet_volume_fraction": (
                installed.vanishing_wet_cut_threshold_authority
                .target_wet_volume_fraction
            ),
            "role": (
                "production"
                if installed.vanishing_wet_cut_threshold_authority
                is cut.P1EF_PRODUCTION_VANISHING_WET_CUT_AUTHORITY
                else "nonproduction_threshold_sensitivity"
            ),
            "scope": "bounded_p1ef_numerical_regularization_only",
            "same_wet_flux_authority_used_for_reconstruction": True,
            "transformed_rows": [
                "last_full_wet_cell_water",
                "last_full_wet_cell_energy",
                "wet_cut_cell_water",
                "wet_cut_cell_energy",
            ],
            "pair_sum_conservation_required": True,
            "runtime_fallback_allowed": False,
            "physically_qualifying": False,
            "plant_predictive": False,
            "normative_gate_status_unchanged": {
                "F2_solvent_particle": False,
                "F3_water_interface": False,
            },
        },
        "target_n2_root": root_audit,
        "source_provenance": provenance_before,
        "model_source_provenance_before": provenance_before,
        "runtime_provenance": runtime_provenance(),
    }
    if first is None:
        provenance_after = source_provenance()
        return {
            **common,
            "status": "FAIL_CLOSED_TARGET_N2_ROOT",
            "passed": False,
            "release_evidence_contract": _failed_release_evidence_contract(
                ("target N=2 root did not pass",)
            ),
            "raw_backend": None,
            "model_source_provenance_after": provenance_after,
            "model_sources_unchanged_during_case": (
                provenance_before == provenance_after
            ),
        }

    backend.N2_ROOT = first.assembly.candidate
    # GT-PS-2-P1E-R7: complete physical-identity context for wall-evidence
    # pinning; cleared afterward so no other caller inherits it implicitly.
    backend.WALL_PIN_CONTEXT = {
        "source_combined_sha256": provenance_before["combined_sha256"],
        "coefficient_case_id": coefficient_case_id,
        "stress_multiplier": stress_multiplier,
        "structural_conductivity_w_m_k": structural_conductivity_w_m_k,
        "pulse_y_hexane": pulse_y_hexane,
    }
    # B-F1: a knob that changes the tangent path is part of the physical
    # identity of a wall pin; entered only when not the on-disk default so a
    # knob-off key keeps its pre-B-F1 form.
    if face_tangent_seed_armed:
        backend.WALL_PIN_CONTEXT["face_tangent_inward_seed_armed"] = True
    if face_tangent_bound_declared:
        backend.WALL_PIN_CONTEXT["face_tangent_residual_bound"] = face_tangent_bound
    face_tangent.drain_face_tangent_call_reports()
    try:
        with _return_hold_installed(return_hold_s):
            raw = backend.run_campaign(meshes, timesteps_s, scenarios)
    finally:
        backend.WALL_PIN_CONTEXT = None
    attach_initial_snapshots_to_raw_jobs(raw, initial_condition)
    release_contract = build_release_evidence_contract(
        installed,
        root_audit,
        raw,
    )
    provenance_after = source_provenance()
    sources_unchanged = provenance_before == provenance_after
    release_contract["model_source_stability"] = {
        "before": provenance_before,
        "after": provenance_after,
        "unchanged_during_case": sources_unchanged,
        "required_for_pass": True,
    }
    if not sources_unchanged:
        release_contract["errors"].append(
            "model sources changed during the public case"
        )
        release_contract["contract_passed"] = False
    numerical_pass = bool(
        raw["all_structural_numerical_contracts_passed"]
        and release_contract["contract_passed"]
        and sources_unchanged
    )
    # B-F1: report-only, attached after every gate and contract was decided.
    face_tangent_block = attach_face_tangent_instrument(
        raw,
        face_tangent.drain_face_tangent_call_reports(),
        armed=face_tangent_seed_armed,
        bound=face_tangent_bound,
        declared=face_tangent_bound_declared,
    )
    payload = {
        **common,
        "status": (
            "PASS_ENGINEERING_NUMERICAL_CASE"
            if numerical_pass
            else "FAIL_CLOSED_ENGINEERING_NUMERICAL_CASE"
        ),
        "passed": numerical_pass,
        "release_evidence_contract": release_contract,
        "raw_backend": raw,
        "model_source_provenance_after": provenance_after,
        "model_sources_unchanged_during_case": sources_unchanged,
        "face_tangent_inward_seed": face_tangent_block,
    }
    if return_hold_s is not None:
        # B-F2: report only, after every gate and contract was decided.
        payload["ablation_return_hold"] = {
            **RETURN_HOLD_DEFINITION,
            "ablation_return_hold_s": return_hold_s,
            "job_holds": [
                {
                    "job_id": job.get("job_id"),
                    "completed": (job.get("return_hold") or {}).get("completed"),
                    "requested_steps": (job.get("return_hold") or {}).get("requested_steps"),
                    "accepted_steps": (job.get("return_hold") or {}).get("accepted_steps"),
                }
                for job in raw.get("jobs") or ()
                if isinstance(job, Mapping) and job.get("completed") is True
            ],
        }
    return payload


FACE_TANGENT_INSTRUMENT_DEFINITION = {
    "record": "B-F1 2026-09-29 face-tangent inward seed build",
    "report_only": True,
    "committed_state_front_speed_m_s": (
        "front speed recovered by the stock flux-matched reconstruction at the committed "
        "arrival state's own (T_Gamma, q_seed): the speed the committed arrival fluxes imply "
        "(M-F1's seed speed)"
    ),
    "tangent_root_front_speed_m_s": "front speed of the accepted exact-face tangent root",
    "front_speed_jump_ratio": "tangent_root_front_speed_m_s / committed_state_front_speed_m_s",
    "mesh_cells": "the job's master-grid cell count",
    "matching": (
        "the tangent-call report of the same cells and arrival face, with an accepted path, "
        "whose time lies in (time_before_s, time_after_s] nearest to "
        "time_before_s + interior_face_event_duration_s"
    ),
    "gates_or_contracts_read_this": False,
}
_ACCEPTED_TANGENT_PATHS = ("stock_certified", "inward_seed_admitted")


def _face_event_instrument_row(
    job: Mapping[str, Any],
    step_index: int,
    step: Mapping[str, Any],
    reports: Sequence[Mapping[str, Any]],
    *,
    bound: float,
    declared: bool,
) -> dict[str, Any]:
    cells = job.get("cells")
    face = step.get("interior_face_event_index")
    start = step.get("time_before_s")
    end = step.get("time_after_s")
    duration = step.get("interior_face_event_duration_s")
    expected = (
        start + duration
        if isinstance(start, (int, float)) and isinstance(duration, (int, float))
        else None
    )
    candidates = [
        report
        for report in reports
        if report.get("cells") == cells
        and report.get("arrival_face_index") == face
        and report.get("path") in _ACCEPTED_TANGENT_PATHS
        and isinstance(report.get("time_s"), float)
        and isinstance(start, (int, float))
        and isinstance(end, (int, float))
        and start < report["time_s"] <= end
    ]
    match = (
        min(
            candidates,
            key=lambda report: (
                abs(report["time_s"] - expected) if expected is not None else 0.0,
                report.get("call_index", 0),
            ),
        )
        if candidates
        else None
    )
    committed = match.get("seed_front_speed_m_s") if match else None
    root = match.get("root_front_speed_m_s") if match else None
    return {
        "step_index": step_index,
        "mesh_cells": cells,
        "dt_s": job.get("dt_s"),
        "arrival_face_index": face,
        "matched": match is not None,
        "matched_report_call_index": match.get("call_index") if match else None,
        "arrival_time_s": match.get("time_s") if match else None,
        "arrival_time_match_error_s": (
            match["time_s"] - expected if match and expected is not None else None
        ),
        "tangent_path": match.get("path") if match else None,
        "committed_state_front_speed_m_s": committed,
        "tangent_root_front_speed_m_s": root,
        "front_speed_jump_ratio": (
            root / committed if root is not None and committed not in (None, 0.0) else None
        ),
        "root_flux_ratio_over_seed": match.get("root_flux_ratio_over_seed") if match else None,
        "root_five_row_scaled_residual": (
            match.get("root_five_row_scaled_residual") if match else None
        ),
        "residual_bound": bound,
        "residual_bound_declared": declared,
    }


def attach_face_tangent_instrument(
    raw: Mapping[str, Any] | None,
    reports: Sequence[Mapping[str, Any]],
    *,
    armed: bool,
    bound: float,
    declared: bool,
) -> dict[str, Any]:
    """Attach the B-F1 knob record and the face-speed instrument (report only).

    Every job record receives the seed knob, the residual bound in force and
    ``residual_bound_declared``; every accepted interior-face step receives a
    ``face_tangent_instrument`` row, collected per job in
    ``face_event_front_speed_instrument``.  Called after every gate and
    contract of the case was decided; nothing reads these fields back.
    """

    report_rows = [dict(report) for report in reports]
    path_counts: dict[str, int] = {}
    for report in report_rows:
        path = str(report.get("path"))
        path_counts[path] = path_counts.get(path, 0) + 1
    jobs = raw.get("jobs") if isinstance(raw, Mapping) else None
    for job in jobs if isinstance(jobs, list) else ():
        if not isinstance(job, dict):
            continue
        job["face_tangent_inward_seed_armed"] = armed
        job["face_tangent_residual_bound"] = bound
        job["residual_bound_declared"] = declared
        rows = []
        for step_index, step in enumerate(job.get("step_records") or ()):
            if not isinstance(step, dict) or step.get("interior_face_event_used") is not True:
                continue
            row = _face_event_instrument_row(
                job, step_index, step, report_rows, bound=bound, declared=declared
            )
            step["face_tangent_instrument"] = row
            rows.append(row)
        job["face_event_front_speed_instrument"] = rows
        hold = job.get("return_hold")
        if isinstance(hold, dict):
            # B-F2: the held steps' face events, matched the same way.
            hold_rows = []
            for hold_index, step in enumerate(hold.get("step_records") or ()):
                if not isinstance(step, dict) or step.get("interior_face_event_used") is not True:
                    continue
                row = _face_event_instrument_row(
                    job, hold_index, step, report_rows, bound=bound, declared=declared
                )
                step["face_tangent_instrument"] = row
                hold_rows.append(row)
            hold["face_event_front_speed_instrument"] = hold_rows
    return {
        "inward_seed_armed": armed,
        "inward_seed_environment_variable": (
            face_tangent.FACE_TANGENT_INWARD_SEED_ENVIRONMENT_VARIABLE
        ),
        "residual_bound": bound,
        "residual_bound_literal": face_tangent._ROOT_RESIDUAL_TOLERANCE,  # noqa: SLF001
        "residual_bound_declared": declared,
        "residual_bound_environment_variable": (
            face_tangent.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE
        ),
        "tangent_call_count": len(report_rows),
        "tangent_call_path_counts": dict(sorted(path_counts.items())),
        "tangent_call_reports": report_rows,
        "instrument_definition": FACE_TANGENT_INSTRUMENT_DEFINITION,
        "physically_qualifying": False,
    }


PER_STEP_CSV_IDENTITY_COLUMNS = ("job_id", "scenario", "cells", "dt_s", "step_index")
PER_STEP_CSV_EXCLUDED_KEYS = frozenset({"wall_time_s"})


def _csv_scalar(value: Any) -> Any:
    if isinstance(value, Mapping) and set(value) == {"__nonfinite_float__"}:
        return value["__nonfinite_float__"]
    if isinstance(value, float):
        return repr(value)
    return value


def per_step_csv_rows(payload: Mapping[str, Any]) -> tuple[list[str], list[dict[str, Any]]]:
    """One row per accepted step record: every scalar field plus the instrument.

    Nested fields are left to the JSON; ``wall_time_s`` is excluded so that two
    runs of the same trajectory give identical rows.
    """

    raw = payload.get("raw_backend")
    jobs = raw.get("jobs") if isinstance(raw, Mapping) else None
    rows: list[dict[str, Any]] = []
    keys: set[str] = set()
    for job in jobs if isinstance(jobs, list) else ():
        if not isinstance(job, Mapping):
            continue
        for step_index, step in enumerate(job.get("step_records") or ()):
            if not isinstance(step, Mapping):
                continue
            row: dict[str, Any] = {
                "job_id": job.get("job_id"),
                "scenario": job.get("scenario"),
                "cells": job.get("cells"),
                "dt_s": _csv_scalar(job.get("dt_s")),
                "step_index": step_index,
            }
            for key, value in step.items():
                if key in PER_STEP_CSV_EXCLUDED_KEYS:
                    continue
                if isinstance(value, Mapping) and set(value) != {"__nonfinite_float__"}:
                    continue
                if isinstance(value, (list, tuple)):
                    continue
                row[key] = _csv_scalar(value)
                keys.add(key)
            instrument = step.get("face_tangent_instrument")
            if isinstance(instrument, Mapping):
                for key, value in instrument.items():
                    if isinstance(value, (Mapping, list, tuple)):
                        continue
                    row[f"ft_{key}"] = _csv_scalar(value)
                    keys.add(f"ft_{key}")
            rows.append(row)
    columns = [*PER_STEP_CSV_IDENTITY_COLUMNS]
    columns.extend(sorted(keys - set(PER_STEP_CSV_IDENTITY_COLUMNS)))
    return columns, rows


def return_hold_csv_rows(
    payload: Mapping[str, Any],
) -> tuple[list[str], list[dict[str, Any]]]:
    """The held steps' rows, in the per-step CSV's own format (B-F2)."""

    raw = payload.get("raw_backend")
    jobs = raw.get("jobs") if isinstance(raw, Mapping) else None
    held_jobs = [
        {**job, "step_records": job["return_hold"].get("step_records") or []}
        for job in (jobs if isinstance(jobs, list) else ())
        if isinstance(job, Mapping) and isinstance(job.get("return_hold"), Mapping)
    ]
    return per_step_csv_rows({"raw_backend": {"jobs": held_jobs}})


def write_per_step_csv(payload: Mapping[str, Any], path: Path) -> None:
    columns, rows = per_step_csv_rows(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, restval="")
        writer.writeheader()
        writer.writerows(rows)


def _parse_csv(values: str, cast) -> tuple:
    return tuple(cast(value.strip()) for value in values.split(",") if value.strip())


def _strict_json_evidence_copy(
    value: Any,
    *,
    location: str,
    nonfinite_paths: list[str],
) -> Any:
    """Copy evidence while losslessly tagging nonfinite diagnostic floats.

    B-F1 (2026-09-29): numpy scalars (a numpy bool from a comparison of numpy
    floats on post-departure roots broke ``json.dumps``) are converted to the
    Python scalar of the same value by ``.item()``; the value is unchanged and
    a nonfinite numpy float is then tagged like any other.
    """

    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            label = "NaN"
        elif value > 0.0:
            label = "+Infinity"
        else:
            label = "-Infinity"
        nonfinite_paths.append(location)
        return {"__nonfinite_float__": label}
    if isinstance(value, Mapping):
        return {
            key: _strict_json_evidence_copy(
                item,
                location=f"{location}.{key}",
                nonfinite_paths=nonfinite_paths,
            )
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [
            _strict_json_evidence_copy(
                item,
                location=f"{location}[{index}]",
                nonfinite_paths=nonfinite_paths,
            )
            for index, item in enumerate(value)
        ]
    return deepcopy(value)


def strict_json_evidence_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Return strict-JSON evidence without turning sentinels into finite data.

    Qualification and solver acceptance must be evaluated before this reporting
    boundary.  A nonfinite diagnostic keeps its exact class and sign in a tagged
    object, which remains nonnumeric (and therefore cannot satisfy a finite-number
    acceptance check).  Finite values and already-computed contract flags are not
    changed.
    """

    nonfinite_paths: list[str] = []
    encoded = _strict_json_evidence_copy(
        payload,
        location="$",
        nonfinite_paths=nonfinite_paths,
    )
    if not isinstance(encoded, dict):
        raise TypeError("strict-JSON evidence payload root must be an object")
    if "strict_json_nonfinite_encoding" in encoded:
        raise ValueError("strict-JSON evidence encoding metadata already exists")
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coefficient-case",
        choices=tuple(case.identifier for case in p1ef.engineering_coefficient_cases()),
        default="nominal_log_center",
    )
    parser.add_argument(
        "--film-multiplier",
        type=float,
        choices=p1ef.COMPACT_FILM_STRESS_MULTIPLIERS,
        default=1.0,
    )
    parser.add_argument("--conductivity", type=float, choices=(0.24, 0.29), default=0.24)
    parser.add_argument("--pulse-temperature-k", type=float, default=PULSE_TEMPERATURE_K)
    parser.add_argument("--pulse-y-hexane", type=float, default=PULSE_Y_HEXANE)
    parser.add_argument(
        "--enforce-p1ef-smooth-wet-retained-water-guard",
        action="store_true",
        help="enable the formal post-root smooth-wet applicability guard",
    )
    parser.add_argument("--meshes", default="12,24,48,96")
    parser.add_argument("--timesteps", default="0.075")
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=(*backend.SCENARIOS, "all"),
        default=("all",),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("numerical_convergence_results/core2_p1ef_dynamic_case_20260801.json"),
    )
    arguments = parser.parse_args()
    payload = run_case(
        coefficient_case_id=arguments.coefficient_case,
        stress_multiplier=arguments.film_multiplier,
        structural_conductivity_w_m_k=arguments.conductivity,
        meshes=_parse_csv(arguments.meshes, int),
        timesteps_s=_parse_csv(arguments.timesteps, float),
        scenarios=arguments.scenarios,
        pulse_temperature_k=arguments.pulse_temperature_k,
        pulse_y_hexane=arguments.pulse_y_hexane,
        enforce_p1ef_smooth_wet_retained_water_guard=(
            arguments.enforce_p1ef_smooth_wet_retained_water_guard
        ),
    )
    strict_payload = strict_json_evidence_payload(payload)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(strict_payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    # B-F1: the per-step CSV beside the record (same stem, ``.steps.csv``).
    write_per_step_csv(strict_payload, arguments.output.with_suffix(".steps.csv"))
    if "ablation_return_hold" in strict_payload:
        # B-F2: the held steps beside it (``.return_hold.steps.csv``).
        columns, rows = return_hold_csv_rows(strict_payload)
        with arguments.output.with_suffix(".return_hold.steps.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, restval="")
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps({"status": payload["status"], "passed": payload["passed"]}))
    raise SystemExit(0 if payload["passed"] else 2)


if __name__ == "__main__":
    main()
