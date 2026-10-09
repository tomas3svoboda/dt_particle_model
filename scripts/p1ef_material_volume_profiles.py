"""Strict material-volume profiles for P1EF mesh/time convergence.

This module only serializes an already accepted immutable cut-integrator
state.  It changes no equation, state, grid, or acceptance decision.  The
coordinate is the normalized spherical material volume ``z=(r/R)^3`` so
profiles from different radial meshes can be compared without interpolating
the conserved state itself.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = 1
COORDINATE = "normalized_material_volume_z=(r/R)^3"


def _finite(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _cell_z_interval(
    cell: object,
    *,
    particle_radius_m: float,
) -> tuple[float, float]:
    lower = (float(cell.inner_radius_m) / particle_radius_m) ** 3
    upper = (float(cell.outer_radius_m) / particle_radius_m) ** 3
    if cell.index == 0:
        lower = 0.0
    return lower, upper


def _validate_partition(
    intervals: Sequence[Sequence[float]],
    temperatures_k: Sequence[float],
    *,
    lower: float,
    upper: float,
    label: str,
) -> None:
    if len(intervals) != len(temperatures_k):
        raise ValueError(f"{label} interval and temperature counts differ")
    if lower == upper:
        if intervals:
            raise ValueError(f"zero-volume {label} region has profile pieces")
        return
    if not intervals:
        raise ValueError(f"positive-volume {label} region has no profile pieces")
    if intervals[0][0] != lower or intervals[-1][1] != upper:
        raise ValueError(f"{label} profile does not meet its exact phase bounds")
    for index, (interval, temperature_k) in enumerate(
        zip(intervals, temperatures_k, strict=True)
    ):
        if (
            len(interval) != 2
            or not all(_finite(value) for value in interval)
            or not 0.0 <= float(interval[0]) < float(interval[1]) <= 1.0
        ):
            raise ValueError(f"{label} profile interval {index} is invalid")
        if not _finite(temperature_k) or float(temperature_k) <= 0.0:
            raise ValueError(f"{label} profile temperature {index} is invalid")
        if index and intervals[index - 1][1] != interval[0]:
            raise ValueError(f"{label} profile intervals are not contiguous")


def serialize_integrator_state(state: object) -> dict[str, Any]:
    """Serialize one accepted partial-front state without modifying it."""

    transport = state.transport
    geometry = transport.geometry
    layout = transport.layout
    radius_m = float(geometry.master_grid.R)
    front_z = float(geometry.front.z)
    if not _finite(radius_m) or radius_m <= 0.0:
        raise ValueError("material-volume profile requires a positive radius")
    if not _finite(front_z) or not 0.0 <= front_z <= 1.0:
        raise ValueError("material-volume profile front lies outside [0,1]")

    wet_intervals: list[list[float]] = []
    for cell_index in layout.wet_cell_indices:
        cell = geometry.cells[cell_index]
        lower, upper = _cell_z_interval(cell, particle_radius_m=radius_m)
        if cell_index == layout.cut_cell_index:
            upper = front_z
        wet_intervals.append([lower, upper])

    dry_intervals: list[list[float]] = []
    for cell_index in layout.dry_cell_indices:
        cell = geometry.cells[cell_index]
        lower, upper = _cell_z_interval(cell, particle_radius_m=radius_m)
        if cell_index == layout.cut_cell_index:
            lower = front_z
        if cell_index == geometry.master_grid.n - 1:
            upper = 1.0
        dry_intervals.append([lower, upper])

    wet_temperatures = [float(value) for value in transport.wet_temperatures_k]
    dry_temperatures = [float(value) for value in transport.dry_temperatures_k]
    _validate_partition(
        wet_intervals,
        wet_temperatures,
        lower=0.0,
        upper=front_z,
        label="wet",
    )
    _validate_partition(
        dry_intervals,
        dry_temperatures,
        lower=front_z,
        upper=1.0,
        label="dry",
    )
    interface_temperature_k = float(state.last_interface_temperature_k)
    time_s = float(transport.time_s)
    if not _finite(time_s) or time_s < 0.0:
        raise ValueError("material-volume profile time is invalid")
    if not _finite(interface_temperature_k) or interface_temperature_k <= 0.0:
        raise ValueError("material-volume interface temperature is invalid")
    return {
        "schema_version": SCHEMA_VERSION,
        "coordinate": COORDINATE,
        "time_s": time_s,
        "front_z": front_z,
        "interface_temperature_k": interface_temperature_k,
        "wet_z_intervals": wet_intervals,
        "wet_temperatures_k": wet_temperatures,
        "dry_z_intervals": dry_intervals,
        "dry_temperatures_k": dry_temperatures,
        "mesh_cells": geometry.master_grid.n,
        "cut_cell_index": geometry.cut_cell_index,
        "complete_exact_partition": True,
        "accepted_state_changed_by_serialization": False,
        "physically_qualifying": False,
    }


def validate_serialized_profile(value: Mapping[str, Any]) -> None:
    """Fail closed on a malformed serialized profile."""

    if value.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("material-volume profile schema changed")
    if value.get("coordinate") != COORDINATE:
        raise ValueError("material-volume profile coordinate changed")
    front_z = value.get("front_z")
    if not _finite(front_z) or not 0.0 <= float(front_z) <= 1.0:
        raise ValueError("material-volume profile front is invalid")
    _validate_partition(
        value.get("wet_z_intervals", ()),
        value.get("wet_temperatures_k", ()),
        lower=0.0,
        upper=float(front_z),
        label="wet",
    )
    _validate_partition(
        value.get("dry_z_intervals", ()),
        value.get("dry_temperatures_k", ()),
        lower=float(front_z),
        upper=1.0,
        label="dry",
    )


__all__ = [
    "COORDINATE",
    "SCHEMA_VERSION",
    "serialize_integrator_state",
    "validate_serialized_profile",
]
