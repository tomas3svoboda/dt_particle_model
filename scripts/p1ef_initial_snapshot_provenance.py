"""Fail-closed provenance for the manufactured P1EF ``t=0`` snapshot.

This module evaluates inventories from an already constructed immutable
``CutIntegratorState``.  It does not initialize, project, or otherwise change
that state.  Raw binary64 values remain in the artifact.  A separate
15-significant-digit decimal projection gives the integrated inventory a
mesh-stable canonical identity while an exact raw digest detects any payload
tampering, including changes below that canonical precision.
"""

from __future__ import annotations

import hashlib
import json
import math
import operator
from typing import Any, Mapping

from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wt


SCHEMA_VERSION = 1
STATE_KIND = "manufactured_partially_wet_t0_snapshot"
CANONICAL_SIGNIFICANT_DECIMAL_DIGITS = 15
WATER_PARTITION_RELATIVE_CLOSURE_LIMIT = 2.0e-14

FALSE_CLAIMS = {
    "feed_prehistory_claimed": False,
    "activation_prehistory_claimed": False,
    "process_prehistory_claimed": False,
    "activation_to_snapshot_continuity_claimed": False,
}
TRUE_CLAIMS = {
    "ledger_origin_is_t0": True,
    "first_relaxation_is_part_of_solved_trajectory": True,
    "historical_wet_hexane_label_is_constitutive_trace_only": True,
    "formal_a1_birth_and_face_event_evidence_remains_separate": True,
}
PROVENANCE_CLAIMS = {**FALSE_CLAIMS, **TRUE_CLAIMS}

_WATER_KEYS = {
    "retained_water_mol",
    "pore_vapor_water_mol",
    "free_liquid_water_mol",
    "partition_total_water_mol",
    "partition_minus_direct_total_water_mol",
    "relative_partition_closure",
    "relative_partition_closure_limit",
}
_SNAPSHOT_KEYS = {
    "schema_version",
    "state_kind",
    "mesh_cells",
    "time_s",
    "radius_m",
    "dry_meal_mass_kg",
    "water_partition_mol",
    "total_water_mol",
    "total_hexane_mol",
    "common_datum_energy_j",
    "whole_particle_equivalent_water_kg_per_kg_dry_meal",
    "whole_particle_equivalent_hexane_kg_per_kg_dry_meal",
    "direct_inventory_equals_reference_inventory",
    "direct_water_cell_inventory_equals_reference_inventory",
    "direct_hexane_cell_inventory_equals_reference_inventory",
    "direct_energy_cell_inventory_equals_reference_inventory",
    "provenance_claims",
    "canonical_significant_decimal_digits",
    "canonical_inventory_projection",
    "canonical_inventory_sha256",
    "raw_snapshot_sha256",
    "accepted_state_changed_by_serialization",
    "physically_qualifying",
}


def _finite(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _canonical_decimal(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError("canonical initial inventory values must be finite")
    return format(value, f".{CANONICAL_SIGNIFICANT_DECIMAL_DIGITS - 1}e")


def _inventory_projection(snapshot: Mapping[str, Any]) -> dict[str, Any]:
    water = snapshot["water_partition_mol"]
    return {
        "schema_version": SCHEMA_VERSION,
        "state_kind": STATE_KIND,
        "dry_meal_mass_kg": _canonical_decimal(float(snapshot["dry_meal_mass_kg"])),
        "water_partition_mol": {
            "retained_water_mol": _canonical_decimal(
                float(water["retained_water_mol"])
            ),
            "pore_vapor_water_mol": _canonical_decimal(
                float(water["pore_vapor_water_mol"])
            ),
            "free_liquid_water_mol": _canonical_decimal(
                float(water["free_liquid_water_mol"])
            ),
        },
        "total_water_mol": _canonical_decimal(float(snapshot["total_water_mol"])),
        "total_hexane_mol": _canonical_decimal(float(snapshot["total_hexane_mol"])),
        "common_datum_energy_j": _canonical_decimal(
            float(snapshot["common_datum_energy_j"])
        ),
        "whole_particle_equivalent_water_kg_per_kg_dry_meal": (
            _canonical_decimal(
                float(
                    snapshot[
                        "whole_particle_equivalent_water_kg_per_kg_dry_meal"
                    ]
                )
            )
        ),
        "whole_particle_equivalent_hexane_kg_per_kg_dry_meal": (
            _canonical_decimal(
                float(
                    snapshot[
                        "whole_particle_equivalent_hexane_kg_per_kg_dry_meal"
                    ]
                )
            )
        ),
    }


def _raw_digest_projection(snapshot: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in snapshot.items()
        if key != "raw_snapshot_sha256"
    }


def serialize_t0_snapshot(state: ci.CutIntegratorState) -> dict[str, Any]:
    """Evaluate and hash one exact manufactured initial state.

    The direct component/energy inventory is deliberately recomputed from the
    transport state.  The initializer's stored reference may certify that
    result only by exact per-cell equality; it is never copied into the output
    as a substitute for the direct evaluation.
    """

    if not isinstance(state, ci.CutIntegratorState):
        raise TypeError("t=0 provenance requires CutIntegratorState")
    transport = state.transport
    if transport.time_s != 0.0:
        raise ValueError("manufactured initial snapshot must have exact time_s=0")

    direct = ci.inventory_snapshot(transport)
    reference = state.reference_inventory
    exact_water = direct.water_cell_mol == reference.water_cell_mol
    exact_hexane = direct.hexane_cell_mol == reference.hexane_cell_mol
    exact_energy = direct.energy_cell_j == reference.energy_cell_j
    if not (exact_water and exact_hexane and exact_energy and direct == reference):
        raise ValueError(
            "direct t=0 component/energy inventory differs from reference_inventory"
        )

    phases = ci.water_phase_inventory_snapshot(transport)
    partition_total = phases.total_water_mol
    partition_residual = math.fsum((partition_total, -direct.total_water_mol))
    partition_scale = max(
        abs(partition_total),
        abs(direct.total_water_mol),
        1.0e-300,
    )
    relative_partition_closure = abs(partition_residual) / partition_scale
    if relative_partition_closure > WATER_PARTITION_RELATIVE_CLOSURE_LIMIT:
        raise ValueError("direct t=0 water partition does not close total water")

    grid = transport.geometry.master_grid
    dry_meal_mass_kg = transport.config.wet.wet.rho_dm_p * grid.total_volume
    if not math.isfinite(dry_meal_mass_kg) or dry_meal_mass_kg <= 0.0:
        raise ValueError("manufactured initial snapshot has invalid dry-meal mass")

    snapshot: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "state_kind": STATE_KIND,
        "mesh_cells": grid.n,
        "time_s": transport.time_s,
        "radius_m": grid.R,
        "dry_meal_mass_kg": dry_meal_mass_kg,
        "water_partition_mol": {
            "retained_water_mol": phases.retained_water_mol,
            "pore_vapor_water_mol": phases.pore_vapor_water_mol,
            "free_liquid_water_mol": phases.free_liquid_water_mol,
            "partition_total_water_mol": partition_total,
            "partition_minus_direct_total_water_mol": partition_residual,
            "relative_partition_closure": relative_partition_closure,
            "relative_partition_closure_limit": (
                WATER_PARTITION_RELATIVE_CLOSURE_LIMIT
            ),
        },
        "total_water_mol": direct.total_water_mol,
        "total_hexane_mol": direct.total_hexane_mol,
        "common_datum_energy_j": direct.total_energy_j,
        "whole_particle_equivalent_water_kg_per_kg_dry_meal": (
            direct.total_water_mol * wt.M / dry_meal_mass_kg
        ),
        "whole_particle_equivalent_hexane_kg_per_kg_dry_meal": (
            direct.total_hexane_mol * hx.M / dry_meal_mass_kg
        ),
        "direct_inventory_equals_reference_inventory": True,
        "direct_water_cell_inventory_equals_reference_inventory": exact_water,
        "direct_hexane_cell_inventory_equals_reference_inventory": exact_hexane,
        "direct_energy_cell_inventory_equals_reference_inventory": exact_energy,
        "provenance_claims": dict(PROVENANCE_CLAIMS),
        "canonical_significant_decimal_digits": (
            CANONICAL_SIGNIFICANT_DECIMAL_DIGITS
        ),
        "canonical_inventory_projection": None,
        "canonical_inventory_sha256": None,
        "raw_snapshot_sha256": None,
        "accepted_state_changed_by_serialization": False,
        "physically_qualifying": False,
    }
    projection = _inventory_projection(snapshot)
    snapshot["canonical_inventory_projection"] = projection
    snapshot["canonical_inventory_sha256"] = _sha256(projection)
    snapshot["raw_snapshot_sha256"] = _sha256(_raw_digest_projection(snapshot))
    validate_t0_snapshot(snapshot)
    return snapshot


def validate_t0_snapshot(value: Mapping[str, Any]) -> None:
    """Fail closed on malformed, mislabelled, or tampered snapshot evidence."""

    if not isinstance(value, Mapping) or set(value) != _SNAPSHOT_KEYS:
        raise ValueError("t=0 snapshot fields changed")
    if value.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("t=0 snapshot schema changed")
    if value.get("state_kind") != STATE_KIND:
        raise ValueError("t=0 snapshot state kind changed")
    try:
        mesh_cells = operator.index(value.get("mesh_cells"))
    except TypeError as exc:
        raise ValueError("t=0 snapshot mesh count is invalid") from exc
    if isinstance(value.get("mesh_cells"), bool) or mesh_cells <= 0:
        raise ValueError("t=0 snapshot mesh count is invalid")
    if value.get("time_s") != 0.0:
        raise ValueError("t=0 snapshot ledger origin changed")
    if not _finite(value.get("radius_m")) or float(value["radius_m"]) <= 0.0:
        raise ValueError("t=0 snapshot radius is invalid")

    scalar_names = (
        "dry_meal_mass_kg",
        "total_water_mol",
        "total_hexane_mol",
        "common_datum_energy_j",
        "whole_particle_equivalent_water_kg_per_kg_dry_meal",
        "whole_particle_equivalent_hexane_kg_per_kg_dry_meal",
    )
    if not all(_finite(value.get(name)) for name in scalar_names):
        raise ValueError("t=0 snapshot has a nonfinite integrated inventory")
    if float(value["dry_meal_mass_kg"]) <= 0.0:
        raise ValueError("t=0 snapshot dry-meal mass is not positive")
    if any(float(value[name]) < 0.0 for name in scalar_names[1:3]):
        raise ValueError("t=0 snapshot component inventory is negative")

    water = value.get("water_partition_mol")
    if not isinstance(water, Mapping) or set(water) != _WATER_KEYS:
        raise ValueError("t=0 snapshot water-partition fields changed")
    if not all(_finite(water.get(name)) for name in _WATER_KEYS):
        raise ValueError("t=0 snapshot water partition is nonfinite")
    component_names = (
        "retained_water_mol",
        "pore_vapor_water_mol",
        "free_liquid_water_mol",
    )
    if any(float(water[name]) < 0.0 for name in component_names):
        raise ValueError("t=0 snapshot water phase inventory is negative")
    recomputed_partition_total = math.fsum(float(water[name]) for name in component_names)
    if float(water["partition_total_water_mol"]) != recomputed_partition_total:
        raise ValueError("t=0 snapshot water partition total was tampered")
    recomputed_residual = math.fsum(
        (recomputed_partition_total, -float(value["total_water_mol"]))
    )
    if float(water["partition_minus_direct_total_water_mol"]) != recomputed_residual:
        raise ValueError("t=0 snapshot water partition residual was tampered")
    scale = max(
        abs(recomputed_partition_total),
        abs(float(value["total_water_mol"])),
        1.0e-300,
    )
    relative = abs(recomputed_residual) / scale
    if float(water["relative_partition_closure"]) != relative:
        raise ValueError("t=0 snapshot water partition closure was tampered")
    if (
        float(water["relative_partition_closure_limit"])
        != WATER_PARTITION_RELATIVE_CLOSURE_LIMIT
        or relative > WATER_PARTITION_RELATIVE_CLOSURE_LIMIT
    ):
        raise ValueError("t=0 snapshot water partition closure failed")

    dry_meal_mass_kg = float(value["dry_meal_mass_kg"])
    expected_water_loading = float(value["total_water_mol"]) * wt.M / dry_meal_mass_kg
    expected_hexane_loading = (
        float(value["total_hexane_mol"]) * hx.M / dry_meal_mass_kg
    )
    if (
        float(value["whole_particle_equivalent_water_kg_per_kg_dry_meal"])
        != expected_water_loading
        or float(value["whole_particle_equivalent_hexane_kg_per_kg_dry_meal"])
        != expected_hexane_loading
    ):
        raise ValueError("t=0 snapshot equivalent loading was tampered")

    exact_reference_flags = (
        "direct_inventory_equals_reference_inventory",
        "direct_water_cell_inventory_equals_reference_inventory",
        "direct_hexane_cell_inventory_equals_reference_inventory",
        "direct_energy_cell_inventory_equals_reference_inventory",
    )
    if any(value.get(name) is not True for name in exact_reference_flags):
        raise ValueError("t=0 snapshot lacks exact reference-inventory equality")
    if value.get("provenance_claims") != PROVENANCE_CLAIMS:
        raise ValueError("t=0 snapshot provenance claims changed")
    if value.get("accepted_state_changed_by_serialization") is not False:
        raise ValueError("t=0 snapshot claims state mutation")
    if value.get("physically_qualifying") is not False:
        raise ValueError("t=0 snapshot was promoted to physical evidence")
    if (
        value.get("canonical_significant_decimal_digits")
        != CANONICAL_SIGNIFICANT_DECIMAL_DIGITS
    ):
        raise ValueError("t=0 snapshot canonical precision changed")

    projection = _inventory_projection(value)
    if value.get("canonical_inventory_projection") != projection:
        raise ValueError("t=0 canonical inventory projection was tampered")
    if value.get("canonical_inventory_sha256") != _sha256(projection):
        raise ValueError("t=0 canonical inventory digest was tampered")
    if value.get("raw_snapshot_sha256") != _sha256(_raw_digest_projection(value)):
        raise ValueError("t=0 raw snapshot digest was tampered")


def serialize_mesh_snapshots(
    states_by_mesh: Mapping[int, ci.CutIntegratorState],
) -> tuple[dict[str, dict[str, Any]], str]:
    """Serialize independently constructed mesh states and require one identity."""

    if not isinstance(states_by_mesh, Mapping) or not states_by_mesh:
        raise ValueError("t=0 snapshot contract needs at least one mesh")
    snapshots: dict[str, dict[str, Any]] = {}
    for raw_cells in sorted(states_by_mesh):
        try:
            cells = operator.index(raw_cells)
        except TypeError as exc:
            raise ValueError("t=0 snapshot mesh key is invalid") from exc
        if isinstance(raw_cells, bool) or cells <= 0 or cells != raw_cells:
            raise ValueError("t=0 snapshot mesh key is invalid")
        snapshot = serialize_t0_snapshot(states_by_mesh[raw_cells])
        if snapshot["mesh_cells"] != cells:
            raise ValueError("t=0 snapshot state does not match its mesh key")
        snapshots[str(cells)] = snapshot

    digests = {snapshot["canonical_inventory_sha256"] for snapshot in snapshots.values()}
    if len(digests) != 1:
        raise ValueError("manufactured t=0 integrated inventory is not mesh invariant")
    return snapshots, next(iter(digests))


__all__ = [
    "CANONICAL_SIGNIFICANT_DECIMAL_DIGITS",
    "FALSE_CLAIMS",
    "PROVENANCE_CLAIMS",
    "SCHEMA_VERSION",
    "STATE_KIND",
    "TRUE_CLAIMS",
    "WATER_PARTITION_RELATIVE_CLOSURE_LIMIT",
    "serialize_mesh_snapshots",
    "serialize_t0_snapshot",
    "validate_t0_snapshot",
]
