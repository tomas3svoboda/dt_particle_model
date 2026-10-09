"""Deterministic mesh-continuation seeds for the Gate-1g cut solver.

This module interpolates an already converged *candidate* onto another master
mesh solely as a nonlinear initial guess.  It never remaps or commits a
conserved state, so it cannot create material inventory or qualify a mesh.
Every target mesh still assembles and closes its own full ALE/RH residual.

Profiles are interpolated in physical radius within their own wet or dry
region.  The front volume coordinate and interface temperature are retained,
and algebraic Stefan fluxes are interpolated on their actual face radii.  No
epsilon phase piece, bound clipping, coefficient interpolation, or residual
relaxation is introduced.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass, field
from typing import Sequence

from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_face_event as face
from dtdc_simulator.core2.particle import cut_face_event_integrator as fi
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_retained_cap as wrc


FACE_CONDITIONED_DRY_COMPOSITION_PREDICTOR = (
    "exact_face_conditioned_gas_only_barrier_coordinate"
)
LEGACY_DRY_COMPOSITION_PREDICTOR = "legacy_raw_physical_mole_fraction"


@dataclass(frozen=True)
class MeshContinuationRecord:
    """Auditable source/target coordinates for one non-conservative seed."""

    source_cell_count: int
    target_cell_count: int
    front_z: float
    source_wet_nodes_m: tuple[float, ...]
    target_wet_nodes_m: tuple[float, ...]
    source_dry_nodes_m: tuple[float, ...]
    target_dry_nodes_m: tuple[float, ...]
    source_dry_face_nodes_m: tuple[float, ...]
    target_dry_face_nodes_m: tuple[float, ...]
    label: str
    conservative_state_remap: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class MeshContinuationSeed:
    """A solver seed plus proof that no accepted state was remapped."""

    seed: ci.CutStepSeed
    record: MeshContinuationRecord
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class MasterFaceEventChart:
    """The exact event chart a crossing continuation candidate must take.

    The mesh-continuation path cannot represent a candidate whose front has
    left the before-state's master cell: the same-cell chart encodes the front
    as a logit inside that one cell, so the face itself is not a representable
    seed coordinate.  Supplying this object does not relax that fact.  It only
    declares the requested macrostep duration the crossing candidate was
    solved over, which is the one datum needed to land the proposal on its
    first master face and to hand the remainder to the exact face chart.

    Nothing here is a tolerance, a clamp, or a new physical number.  Every
    quantity of the route is re-derived from the target's own master grid by
    :func:`cut_geometry.plan_primary_drainage_step`.
    """

    requested_dt_s: float
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.requested_dt_s) or self.requested_dt_s <= 0.0:
            raise ValueError(
                "master-face event chart needs a positive finite requested duration"
            )
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("master-face event chart provenance label is required")


@dataclass(frozen=True)
class MasterFaceEventRecord:
    """Auditable coordinates of one crossing candidate landed on its face."""

    source_cell_count: int
    target_cell_count: int
    before_front_z: float
    requested_front_z: float
    accepted_front_z: float
    first_event_face_index: int
    encountered_face_indices: tuple[int, ...]
    unresolved_face_indices: tuple[int, ...]
    requested_dt_s: float
    accepted_dt_s: float
    accepted_proposal_fraction: float
    remaining_dt_s: float
    front_velocity_z_s: float
    classification: str
    requires_refinement: bool
    label: str
    conservative_state_remap: bool = field(default=False, init=False)
    must_resolve_accepted_event: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class MasterFaceEventContinuation:
    """Exact-event-chart inputs for a candidate that crosses a master face.

    ``accepted_dt_s`` on the record is a volume-coordinate interpolation of the
    crossing proposal, exactly as :class:`cut_geometry.FrontStepPlan` defines
    it.  It is a seed for the event-time unknown and never a converged event
    time: the arrival root and its separate ordering evidence stay mandatory,
    and the unresolved remainder must be continued from the face.
    """

    before: ci.CutIntegratorState
    same_cell_seed: ci.CutStepSeed
    face_event_controls: fi.FaceEventControls
    face_event_seed: fi.FaceEventSeed
    record: MasterFaceEventRecord
    conservative_state_remap: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class TimeContinuationRecord:
    """Audit record for linear-in-step predictor coordinates only."""

    source_dt_s: float
    target_dt_s: float
    increment_ratio: float
    dry_composition_predictor_coordinate: str
    label: str
    source_surface_boundary_included: bool = field(default=True, init=False)
    target_surface_boundary_included: bool = field(default=True, init=False)
    composition_clipped_or_projected: bool = field(default=False, init=False)
    algebraic_fluxes_rescaled: bool = field(default=False, init=False)
    accepted_state_changed: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class TimeContinuationSeed:
    seed: ci.CutStepSeed
    record: TimeContinuationRecord
    physically_qualifying: bool = field(default=False, init=False)


def prolongate_same_front_candidate(
    source_before: ci.CutIntegratorState,
    source_candidate: cut.CutTransportUnknowns,
    target_before: ci.CutIntegratorState,
    *,
    stefan_flux_scale_floor_mol_m2_s: float,
    label: str,
    master_face_event_chart: MasterFaceEventChart | None = None,
) -> MeshContinuationSeed | MasterFaceEventContinuation:
    """Interpolate a strict same-cell root into a target-mesh initial guess.

    A candidate whose front has left the before-state's master cell has no
    same-cell representation on the target mesh.  Without
    ``master_face_event_chart`` that case is refused with its established
    text.  With one, the proposal is landed on its first master face and the
    exact face-arrival chart's own inputs are returned instead: the step is
    accepted at the face and its remainder is continued from the face.
    """

    if not isinstance(source_before, ci.CutIntegratorState) or not isinstance(
        target_before, ci.CutIntegratorState
    ):
        raise TypeError("mesh continuation requires source and target integrator states")
    if not isinstance(source_candidate, cut.CutTransportUnknowns):
        raise TypeError("mesh continuation requires one converged cut candidate")
    if (
        not math.isfinite(stefan_flux_scale_floor_mol_m2_s)
        or stefan_flux_scale_floor_mol_m2_s <= 0.0
    ):
        raise ValueError("Stefan-flux scale floor must be positive and finite")
    if not isinstance(label, str) or not label.strip():
        raise ValueError("mesh-continuation provenance label is required")

    source = source_before.transport
    target = target_before.transport
    if source.config != target.config:
        raise ValueError("mesh continuation cannot change the transport configuration")
    if source.geometry.master_grid.R != target.geometry.master_grid.R:
        raise ValueError("mesh continuation requires one physical particle radius")
    if source.geometry.front.z != target.geometry.front.z:
        raise ValueError("mesh continuation requires the exact same before-front coordinate")
    if len(source_candidate.equation_rank_vector()) != source.layout.unknown_count:
        raise ValueError("source candidate does not match its source cut rank")
    if not 0.0 < source_candidate.front_z < source.geometry.front.z:
        raise ValueError("source candidate must be a strict receding partial-front step")

    source_geometry = cg.partition_master_grid(
        source.geometry.master_grid,
        z=source_candidate.front_z,
    )
    target_geometry = cg.partition_master_grid(
        target.geometry.master_grid,
        z=source_candidate.front_z,
    )
    if (
        source_geometry.cut_cell_index != source.layout.cut_cell_index
        or target_geometry.cut_cell_index != target.layout.cut_cell_index
    ):
        if master_face_event_chart is None:
            raise ValueError(
                "candidate crosses a master face; use the exact event chart, not mesh continuation"
            )
        return _route_master_face_event(
            source_before,
            source_candidate,
            target_before,
            source_geometry=source_geometry,
            target_geometry=target_geometry,
            chart=master_face_event_chart,
            stefan_flux_scale_floor_mol_m2_s=stefan_flux_scale_floor_mol_m2_s,
            label=label,
        )

    source_wet_nodes = _wet_nodes(source_geometry, source.layout)
    target_wet_nodes = _wet_nodes(target_geometry, target.layout)
    source_dry_nodes = _dry_nodes(source_geometry, source.layout)
    target_dry_nodes = _dry_nodes(target_geometry, target.layout)
    source_face_nodes = _dry_face_nodes(source_geometry, source.layout)
    target_face_nodes = _dry_face_nodes(target_geometry, target.layout)

    wet_t = _interpolate_profile(
        source_wet_nodes,
        source_candidate.wet_temperatures_k,
        target_wet_nodes,
    )
    source_cap_chart = _wet_cap_chart(source_before)
    target_cap_chart = _wet_cap_chart(target_before)
    if (source_cap_chart is None) != (target_cap_chart is None):
        raise ValueError("mesh continuation cannot change the wet retained-cap chart")
    if source_cap_chart is None:
        if any(
            value != 0.0
            for value in source_candidate.effective_wet_retained_water_capacity_duals_over_rt
        ):
            raise ValueError("smooth mesh continuation cannot discard a wet retained-cap dual")
        wet_w = _interpolate_profile(
            source_wet_nodes,
            source_candidate.wet_retained_water_loadings,
            target_wet_nodes,
        )
        wet_duals: tuple[float, ...] = (0.0,) * len(wet_w)
    else:
        source_graph = tuple(
            source_cap_chart.encode(loading, dual)
            for loading, dual in zip(
                source_candidate.wet_retained_water_loadings,
                source_candidate.effective_wet_retained_water_capacity_duals_over_rt,
            )
        )
        target_graph = _interpolate_profile(
            source_wet_nodes,
            source_graph,
            target_wet_nodes,
        )
        target_points = tuple(target_cap_chart.decode(value) for value in target_graph)
        wet_w = tuple(point.retained_water_loading for point in target_points)
        wet_duals = tuple(point.capacity_dual_over_rt for point in target_points)
    dry_t = _interpolate_profile(
        source_dry_nodes,
        source_candidate.dry_temperatures_k,
        target_dry_nodes,
    )
    dry_y = _interpolate_profile(
        source_dry_nodes,
        source_candidate.dry_y_hexane,
        target_dry_nodes,
    )
    nt = _interpolate_profile(
        source_face_nodes,
        source_candidate.dry_total_stefan_fluxes_mol_m2_s,
        target_face_nodes,
    )
    candidate = cut.CutTransportUnknowns(
        wet_temperatures_k=wet_t,
        wet_retained_water_loadings=wet_w,
        dry_temperatures_k=dry_t,
        dry_y_hexane=dry_y,
        dry_total_stefan_fluxes_mol_m2_s=nt,
        front_z=source_candidate.front_z,
        interface_temperature_k=source_candidate.interface_temperature_k,
        wet_retained_water_capacity_duals_over_rt=wet_duals,
    )
    scales = tuple(
        max(abs(value), stefan_flux_scale_floor_mol_m2_s) for value in nt
    )
    seed = ci.CutStepSeed(
        candidate=candidate,
        stefan_flux_scales_mol_m2_s=scales,
        label=(
            f"{label}; radial primitive/face interpolation for nonlinear seeding only; "
            "no accepted-state or inventory remap"
        ),
    )
    # Exercise the target's exact strict-open chart before returning the seed.
    ci._validate_seed(target_before, seed)  # noqa: SLF001
    record = MeshContinuationRecord(
        source_cell_count=source.geometry.master_grid.n,
        target_cell_count=target.geometry.master_grid.n,
        front_z=source_candidate.front_z,
        source_wet_nodes_m=source_wet_nodes,
        target_wet_nodes_m=target_wet_nodes,
        source_dry_nodes_m=source_dry_nodes,
        target_dry_nodes_m=target_dry_nodes,
        source_dry_face_nodes_m=source_face_nodes,
        target_dry_face_nodes_m=target_face_nodes,
        label=label,
    )
    return MeshContinuationSeed(seed=seed, record=record)


def _route_master_face_event(
    source_before: ci.CutIntegratorState,
    source_candidate: cut.CutTransportUnknowns,
    target_before: ci.CutIntegratorState,
    *,
    source_geometry: cg.CutGeometry,
    target_geometry: cg.CutGeometry,
    chart: MasterFaceEventChart,
    stefan_flux_scale_floor_mol_m2_s: float,
    label: str,
) -> MasterFaceEventContinuation:
    """Land one crossing proposal on its first master face and seed the event.

    The source candidate must still be a strict same-cell root of its own
    mesh: the exact event chart of a mesh is taken on that mesh, never on a
    neighbour's behalf.  On the target mesh the proposal is canonicalised by
    the unchanged :func:`cut_geometry.plan_primary_drainage_step`, which
    accepts it only to the first encountered face.  The source candidate's
    converged profiles are interpolated in physical radius onto the arrival
    configuration exactly as the same-cell path interpolates them; nothing is
    clipped, remapped, or rescaled, and no epsilon phase piece is introduced.
    """

    source = source_before.transport
    target = target_before.transport
    if source_geometry.cut_cell_index != source.layout.cut_cell_index:
        raise ValueError(
            "the source candidate crosses a master face of its own mesh; the "
            "exact event chart must be taken on the mesh that produced it"
        )
    cut_index = target.layout.cut_cell_index
    if cut_index is None or cut_index <= 0:
        raise ValueError(
            "master-face event routing requires one strict interior target cut cell"
        )
    if any(
        value != 0.0
        for value in source_candidate.effective_wet_retained_water_capacity_duals_over_rt
    ) or any(
        value != 0.0
        for value in target.effective_wet_retained_water_capacity_duals_over_rt
    ):
        raise ValueError(
            "the exact face-arrival chart does not carry the wet retained-cap "
            "graph dual; an active cap fails closed"
        )

    plan = cg.plan_primary_drainage_step(
        target.geometry,
        target_geometry,
        chart.requested_dt_s,
    )
    if not plan.must_resolve_accepted_event:
        raise ValueError(
            "master-face event routing requires a proposal that crosses a face"
        )
    first_face = plan.first_event_face_index
    if first_face != cut_index:
        raise ValueError(
            "the first master-face event is not the target cut cell's inner face"
        )
    arrival_layout = face.layout_for_arrival(target)
    if arrival_layout.arrival_face_index != first_face:
        raise ValueError("exact face-arrival rank does not match the planned event")

    accepted = plan.accepted
    arrival_radius_m = accepted.front.radius_m
    arrival_wet_nodes = tuple(
        cut._wet_piece_center(accepted.cells[cell], arrival_radius_m)  # noqa: SLF001
        for cell in arrival_layout.wet_cell_indices
    )
    arrival_dry_nodes = tuple(
        cut._dry_piece_center(accepted.cells[cell], arrival_radius_m)  # noqa: SLF001
        for cell in arrival_layout.dry_cell_indices
    )
    arrival_face_nodes = (
        arrival_radius_m,
        *(
            accepted.master_grid.faces[cell + 1]
            for cell in arrival_layout.dry_cell_indices
        ),
    )

    source_wet_nodes = _wet_nodes(source_geometry, source.layout)
    source_dry_nodes = _dry_nodes(source_geometry, source.layout)
    source_face_nodes = _dry_face_nodes(source_geometry, source.layout)
    arrival = face.FaceArrivalUnknowns(
        wet_temperatures_k=_interpolate_profile(
            source_wet_nodes,
            source_candidate.wet_temperatures_k,
            arrival_wet_nodes,
        ),
        wet_retained_water_loadings=_interpolate_profile(
            source_wet_nodes,
            source_candidate.wet_retained_water_loadings,
            arrival_wet_nodes,
        ),
        dry_temperatures_k=_interpolate_profile(
            source_dry_nodes,
            source_candidate.dry_temperatures_k,
            arrival_dry_nodes,
        ),
        dry_y_hexane=_interpolate_profile(
            source_dry_nodes,
            source_candidate.dry_y_hexane,
            arrival_dry_nodes,
        ),
        dry_total_stefan_fluxes_mol_m2_s=_interpolate_profile(
            source_face_nodes,
            source_candidate.dry_total_stefan_fluxes_mol_m2_s,
            arrival_face_nodes,
        ),
        event_time_s=plan.accepted_dt_s,
        interface_temperature_k=source_candidate.interface_temperature_k,
    )
    if len(arrival.vector()) != arrival_layout.unknown_count:
        raise ValueError("master-face event seed lost the exact arrival rank")
    flux_scales = tuple(
        max(abs(value), stefan_flux_scale_floor_mol_m2_s)
        for value in arrival.dry_total_stefan_fluxes_mol_m2_s
    )
    face_event_seed = fi.FaceEventSeed(
        arrival,
        flux_scales,
        (
            f"{label}; radial primitive/face interpolation of the crossing "
            "source root onto exact-face rank for nonlinear seeding only; "
            "planned first-face volume-coordinate arrival time; no accepted "
            "state, inventory, or event time is committed here"
        ),
    )

    before_front_z = target.geometry.front.z
    front_velocity_z_s = (
        math.fsum((source_candidate.front_z, -before_front_z)) / chart.requested_dt_s
    )
    same_cell_seed = _master_face_same_cell_seed(
        target_before,
        front_velocity_z_s=front_velocity_z_s,
        predicted_event_duration_s=plan.accepted_dt_s,
        inner_face_z=accepted.front.z,
        flux_scales_mol_m2_s=flux_scales,
        label=label,
    )
    record = MasterFaceEventRecord(
        source_cell_count=source.geometry.master_grid.n,
        target_cell_count=target.geometry.master_grid.n,
        before_front_z=before_front_z,
        requested_front_z=source_candidate.front_z,
        accepted_front_z=accepted.front.z,
        first_event_face_index=first_face,
        encountered_face_indices=plan.encountered_face_indices,
        unresolved_face_indices=(
            () if plan.remaining is None else plan.remaining.unresolved_face_indices
        ),
        requested_dt_s=plan.requested_dt_s,
        accepted_dt_s=plan.accepted_dt_s,
        accepted_proposal_fraction=plan.accepted_proposal_fraction,
        remaining_dt_s=plan.remaining_dt_s,
        front_velocity_z_s=front_velocity_z_s,
        classification=plan.classification.value,
        requires_refinement=plan.requires_refinement,
        label=label,
    )
    return MasterFaceEventContinuation(
        before=target_before,
        same_cell_seed=same_cell_seed,
        face_event_controls=fi.FaceEventControls((0.0, chart.requested_dt_s)),
        face_event_seed=face_event_seed,
        record=record,
    )


def _master_face_same_cell_seed(
    target_before: ci.CutIntegratorState,
    *,
    front_velocity_z_s: float,
    predicted_event_duration_s: float,
    inner_face_z: float,
    flux_scales_mol_m2_s: tuple[float, ...],
    label: str,
) -> ci.CutStepSeed:
    """Return the strict-interior same-cell seed the event route must refuse first.

    The exact event chart is entered only after the unchanged same-cell
    equations reject with exact rollback, so this routine must still hand them
    a seed inside their own chart.  The planned face-arrival horizon is
    backtracked dyadically along the source root's own front tangent; no other
    primitive is extrapolated, clipped, or projected.
    """

    state = target_before.transport
    before_front_z = state.geometry.front.z
    last_error: Exception | None = None
    for backtrack_exponent in range(48):
        seed_duration_s = math.ldexp(
            predicted_event_duration_s,
            -(backtrack_exponent + 1),
        )
        front_z = math.fsum(
            (before_front_z, front_velocity_z_s * seed_duration_s)
        )
        if not inner_face_z < front_z < before_front_z:
            last_error = ValueError(
                "master-face same-cell seed is not strictly inside the current "
                "material-cell chart"
            )
            continue
        seed = ci.CutStepSeed(
            cut.CutTransportUnknowns(
                wet_temperatures_k=state.wet_temperatures_k,
                wet_retained_water_loadings=state.wet_retained_water_loadings,
                dry_temperatures_k=state.dry_temperatures_k,
                dry_y_hexane=state.dry_y_hexane,
                dry_total_stefan_fluxes_mol_m2_s=(
                    target_before.last_total_stefan_fluxes_mol_m2_s
                ),
                front_z=front_z,
                interface_temperature_k=(
                    target_before.last_interface_temperature_k
                ),
                wet_retained_water_capacity_duals_over_rt=(
                    state.effective_wet_retained_water_capacity_duals_over_rt
                ),
            ),
            flux_scales_mol_m2_s,
            (
                f"{label}; master-face event route same-cell seed; "
                "half-face-horizon dyadic backtrack exponent "
                f"{backtrack_exponent}; source-root front tangent only, no "
                "field extrapolation and no clipping"
            ),
        )
        try:
            ci._validate_seed(target_before, seed)  # noqa: SLF001
        except ValueError as exc:
            last_error = exc
            continue
        return seed
    raise ValueError(
        "master-face event route could not enter the strict same-cell chart "
        "after deterministic dyadic backtracking"
    ) from last_error


def rescale_same_before_increment(
    before: ci.CutIntegratorState,
    source_candidate: cut.CutTransportUnknowns,
    *,
    surface_boundary: ct.PoreBoundary,
    source_dt_s: float,
    target_dt_s: float,
    stefan_flux_scale_floor_mol_m2_s: float,
    label: str,
) -> TimeContinuationSeed:
    """Predict another timestep root from a root with the exact same before state.

    Differential primitive and front increments are scaled linearly in ``dt``.
    Algebraic Stefan fluxes are rates and are therefore retained, not scaled.
    Interface temperature is an algebraic trace, but its departure from the
    last accepted trace is used as a smooth continuation predictor.  The
    exact gas-only dry composition is scaled in the solver's face-conditioned
    barrier coordinate, including the actual surface boundary.  The returned
    object remains only a seed and is validated against the unchanged
    strict-open target chart.
    """

    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("time continuation requires a cut integrator state")
    if not isinstance(source_candidate, cut.CutTransportUnknowns):
        raise TypeError("time continuation requires a cut candidate")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("time continuation requires the actual surface boundary")
    durations = (source_dt_s, target_dt_s)
    if not all(math.isfinite(value) and value > 0.0 for value in durations):
        raise ValueError("continuation timesteps must be positive and finite")
    if (
        not math.isfinite(stefan_flux_scale_floor_mol_m2_s)
        or stefan_flux_scale_floor_mol_m2_s <= 0.0
    ):
        raise ValueError("Stefan-flux scale floor must be positive and finite")
    if not isinstance(label, str) or not label.strip():
        raise ValueError("time-continuation provenance label is required")
    state = before.transport
    layout = state.layout
    if len(source_candidate.equation_rank_vector()) != layout.unknown_count:
        raise ValueError("source candidate does not match the before-state rank")
    ratio = target_dt_s / source_dt_s

    def scaled(old: Sequence[float], new: Sequence[float]) -> tuple[float, ...]:
        if len(old) != len(new):
            raise ValueError("time-continuation primitive blocks are misaligned")
        return tuple(
            old_value + ratio * (new_value - old_value)
            for old_value, new_value in zip(old, new)
        )

    cap_chart = _wet_cap_chart(before)
    if cap_chart is None:
        if any(
            value != 0.0
            for value in (
                *state.effective_wet_retained_water_capacity_duals_over_rt,
                *source_candidate.effective_wet_retained_water_capacity_duals_over_rt,
            )
        ):
            raise ValueError("smooth time continuation cannot discard a wet retained-cap dual")
        wet_w = scaled(
            state.wet_retained_water_loadings,
            source_candidate.wet_retained_water_loadings,
        )
        wet_duals = (0.0,) * len(wet_w)
    else:
        old_graph = tuple(
            cap_chart.encode(loading, dual)
            for loading, dual in zip(
                state.wet_retained_water_loadings,
                state.effective_wet_retained_water_capacity_duals_over_rt,
            )
        )
        new_graph = tuple(
            cap_chart.encode(loading, dual)
            for loading, dual in zip(
                source_candidate.wet_retained_water_loadings,
                source_candidate.effective_wet_retained_water_capacity_duals_over_rt,
            )
        )
        predicted_graph = scaled(old_graph, new_graph)
        predicted_points = tuple(cap_chart.decode(value) for value in predicted_graph)
        wet_w = tuple(point.retained_water_loading for point in predicted_points)
        wet_duals = tuple(point.capacity_dual_over_rt for point in predicted_points)

    predicted_wet_temperatures = scaled(
        state.wet_temperatures_k,
        source_candidate.wet_temperatures_k,
    )
    predicted_dry_temperatures = scaled(
        state.dry_temperatures_k,
        source_candidate.dry_temperatures_k,
    )
    predicted_interface_temperature = before.last_interface_temperature_k + ratio * (
        source_candidate.interface_temperature_k
        - before.last_interface_temperature_k
    )
    predicted_dry_y, dry_composition_route = affine_dry_composition_seed_values(
        before,
        source_surface_boundary=surface_boundary,
        target_surface_boundary=surface_boundary,
        old_dry_temperatures_k=state.dry_temperatures_k,
        old_interface_temperature_k=before.last_interface_temperature_k,
        old_dry_y_hexane=state.dry_y_hexane,
        new_dry_temperatures_k=source_candidate.dry_temperatures_k,
        new_interface_temperature_k=source_candidate.interface_temperature_k,
        new_dry_y_hexane=source_candidate.dry_y_hexane,
        predicted_dry_temperatures_k=predicted_dry_temperatures,
        predicted_interface_temperature_k=predicted_interface_temperature,
        affine_multiplier=ratio,
    )

    candidate = cut.CutTransportUnknowns(
        wet_temperatures_k=predicted_wet_temperatures,
        wet_retained_water_loadings=wet_w,
        dry_temperatures_k=predicted_dry_temperatures,
        dry_y_hexane=predicted_dry_y,
        dry_total_stefan_fluxes_mol_m2_s=(
            source_candidate.dry_total_stefan_fluxes_mol_m2_s
        ),
        front_z=state.geometry.front.z
        + ratio * (source_candidate.front_z - state.geometry.front.z),
        interface_temperature_k=predicted_interface_temperature,
        wet_retained_water_capacity_duals_over_rt=wet_duals,
    )
    scales = tuple(
        max(abs(value), stefan_flux_scale_floor_mol_m2_s)
        for value in candidate.dry_total_stefan_fluxes_mol_m2_s
    )
    seed = ci.CutStepSeed(
        candidate=candidate,
        stefan_flux_scales_mol_m2_s=scales,
        label=(
            f"{label}; linear differential-increment timestep predictor only; "
            f"dry composition in {dry_composition_route}; "
            "algebraic flux rates unchanged; no accepted-state change"
        ),
    )
    ci._validate_seed(before, seed, surface_boundary)  # noqa: SLF001
    return TimeContinuationSeed(
        seed=seed,
        record=TimeContinuationRecord(
            source_dt_s=source_dt_s,
            target_dt_s=target_dt_s,
            increment_ratio=ratio,
            dry_composition_predictor_coordinate=dry_composition_route,
            label=label,
        ),
    )


def affine_dry_composition_seed_values(
    chart_state: ci.CutIntegratorState,
    *,
    source_surface_boundary: ct.PoreBoundary,
    target_surface_boundary: ct.PoreBoundary,
    old_dry_temperatures_k: Sequence[float],
    old_interface_temperature_k: float,
    old_dry_y_hexane: Sequence[float],
    new_dry_temperatures_k: Sequence[float],
    new_interface_temperature_k: float,
    new_dry_y_hexane: Sequence[float],
    predicted_dry_temperatures_k: Sequence[float],
    predicted_interface_temperature_k: float,
    affine_multiplier: float,
) -> tuple[tuple[float, ...], str]:
    """Return one affine dry-composition seed in the authoritative chart.

    ``old + affine_multiplier * (new - old)`` is applied to the exact
    face-conditioned gas-only barrier coordinate.  Old and new accepted
    values are encoded with their source boundary and their own temperature
    contexts; the predicted coordinate is decoded with the target boundary
    and predicted temperatures.  This transports a seed across a boundary
    jump without clipping a physical mole fraction or changing an accepted
    state.  The legacy rectangular primitive band deliberately retains its
    established raw-``y_h`` arithmetic bit for bit.
    """

    if not isinstance(chart_state, ci.CutIntegratorState):
        raise TypeError("dry-composition prediction requires a cut integrator state")
    for name, boundary in (
        ("source", source_surface_boundary),
        ("target", target_surface_boundary),
    ):
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError(f"{name} dry-composition predictor boundary is invalid")
    if not math.isfinite(affine_multiplier):
        raise ValueError("dry-composition affine multiplier must be finite")

    old_t = tuple(old_dry_temperatures_k)
    old_y = tuple(old_dry_y_hexane)
    new_t = tuple(new_dry_temperatures_k)
    new_y = tuple(new_dry_y_hexane)
    predicted_t = tuple(predicted_dry_temperatures_k)
    source_chart = ci._coordinate_chart(  # noqa: SLF001
        chart_state,
        source_surface_boundary,
    )
    target_chart = ci._coordinate_chart(  # noqa: SLF001
        chart_state,
        target_surface_boundary,
    )
    expected_rank = len(source_chart.dry_pores)
    if len(target_chart.dry_pores) != expected_rank or any(
        len(values) != expected_rank
        for values in (old_t, old_y, new_t, new_y, predicted_t)
    ):
        raise ValueError("dry-composition predictor blocks do not match the cut rank")
    if (source_chart.dry_y_hexane is None) != (target_chart.dry_y_hexane is None):
        raise ValueError("source and target dry-composition charts changed representation")

    if source_chart.dry_y_hexane is not None:
        values = tuple(
            old_value + affine_multiplier * (new_value - old_value)
            for old_value, new_value in zip(old_y, new_y)
        )
        return values, LEGACY_DRY_COMPOSITION_PREDICTOR

    old_coordinates = tuple(
        ci._encode_dry_y(  # noqa: SLF001
            piece,
            old_t,
            old_interface_temperature_k,
            value,
            source_chart,
        )
        for piece, value in enumerate(old_y)
    )
    new_coordinates = tuple(
        ci._encode_dry_y(  # noqa: SLF001
            piece,
            new_t,
            new_interface_temperature_k,
            value,
            source_chart,
        )
        for piece, value in enumerate(new_y)
    )
    predicted_coordinates = tuple(
        old_value + affine_multiplier * (new_value - old_value)
        for old_value, new_value in zip(old_coordinates, new_coordinates)
    )
    values = tuple(
        ci._decode_dry_y(  # noqa: SLF001
            piece,
            predicted_t,
            predicted_interface_temperature_k,
            coordinate,
            target_chart,
        )
        for piece, coordinate in enumerate(predicted_coordinates)
    )
    return values, FACE_CONDITIONED_DRY_COMPOSITION_PREDICTOR


def _wet_cap_chart(
    state: ci.CutIntegratorState,
) -> wrc.WetRetainedCapSemismoothChart | None:
    """Return the prospectively enabled exact cap chart, if any."""

    lower, upper = state.controls.wet_water_bounds
    if upper != state.transport.config.wet.luikov.W_cap:
        return None
    return wrc.WetRetainedCapSemismoothChart(
        lower,
        state.transport.config.wet.luikov,
    )


def _wet_nodes(
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
) -> tuple[float, ...]:
    return tuple(
        cut._wet_piece_center(geometry.cells[cell], geometry.front.radius_m)  # noqa: SLF001
        for cell in layout.wet_cell_indices
    )


def _dry_nodes(
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
) -> tuple[float, ...]:
    return tuple(
        cut._dry_piece_center(geometry.cells[cell], geometry.front.radius_m)  # noqa: SLF001
        for cell in layout.dry_cell_indices
    )


def _dry_face_nodes(
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
) -> tuple[float, ...]:
    return (
        geometry.front.radius_m,
        *(geometry.master_grid.faces[cell + 1] for cell in layout.dry_cell_indices),
    )


def _interpolate_profile(
    source_nodes: Sequence[float],
    source_values: Sequence[float],
    target_nodes: Sequence[float],
) -> tuple[float, ...]:
    nodes = tuple(source_nodes)
    values = tuple(source_values)
    targets = tuple(target_nodes)
    if len(nodes) != len(values) or not nodes:
        raise ValueError("source continuation profile must be non-empty and aligned")
    if not all(math.isfinite(value) for value in (*nodes, *values, *targets)):
        raise ValueError("continuation profiles must be finite")
    if any(right <= left for left, right in zip(nodes, nodes[1:])):
        raise ValueError("source continuation nodes must be strictly increasing")
    result: list[float] = []
    for target in targets:
        if target <= nodes[0]:
            result.append(values[0])
            continue
        if target >= nodes[-1]:
            result.append(values[-1])
            continue
        right = bisect.bisect_right(nodes, target)
        left = right - 1
        fraction = (target - nodes[left]) / (nodes[right] - nodes[left])
        result.append(values[left] + fraction * (values[right] - values[left]))
    return tuple(result)


__all__ = [
    "FACE_CONDITIONED_DRY_COMPOSITION_PREDICTOR",
    "LEGACY_DRY_COMPOSITION_PREDICTOR",
    "MasterFaceEventChart",
    "MasterFaceEventContinuation",
    "MasterFaceEventRecord",
    "MeshContinuationRecord",
    "MeshContinuationSeed",
    "TimeContinuationRecord",
    "TimeContinuationSeed",
    "affine_dry_composition_seed_values",
    "prolongate_same_front_candidate",
    "rescale_same_before_increment",
]
