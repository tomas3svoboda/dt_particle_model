r"""Safeguarded exact-face to adjacent-face event solve.

This is the finite-cell continuation of :mod:`cut_face_event_integrator`.
The immutable start is exact master face ``K > 1`` and the solved endpoint is
exact face ``K - 1``.  Event duration is one of the simultaneous unknowns.
There is no overlap interval and no positive-volume surrogate at either end.

Only the frozen wet-water/conduction, binary-pore, retained-water, ALE, and
Rankine--Hugoniot equations assembled by :mod:`cut_face_to_face` occur here.
The coordinate maps, scales, nonlinear method, and rank tests are numerical
oracles and are explicitly nonqualifying physics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np
from scipy import optimize

from dtdc_simulator.core2.particle import conditioning
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_face_departure_integrator as di
from dtdc_simulator.core2.particle import cut_face_event_integrator as fi
from dtdc_simulator.core2.particle import cut_face_to_face as ff
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.props import hexane as hx


class FaceToFaceIntegratorStepError(RuntimeError):
    """Rejected solve carrying the identical immutable exact-face state."""

    def __init__(
        self,
        message: str,
        rollback_state: fi.FaceArrivalCommittedState,
        *,
        nonlinear_evaluations: int,
        rejected_trial_evaluations: int,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        last_candidate: ff.FaceToFaceUnknowns | None = None,
        scaled_residuals: tuple[float, ...] = (),
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (
            ci.WetRetainedWaterApplicabilityAudit | None
        ) = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.rejected_trial_evaluations = rejected_trial_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.last_candidate = last_candidate
        self.scaled_residuals = scaled_residuals
        self.retained_water_applicability_guard_failed = (
            retained_water_applicability_guard_failed
        )
        self.retained_water_applicability_audit = (
            retained_water_applicability_audit
        )


@dataclass(frozen=True)
class FaceToFaceSeed:
    """One explicit physical seed and one scale per dry Stefan face."""

    candidate: ff.FaceToFaceUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, ff.FaceToFaceUnknowns):
            raise TypeError("face-to-face seed requires FaceToFaceUnknowns")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0
            for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every face-to-face Stefan face needs a positive scale")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("face-to-face seed needs a provenance label")


@dataclass(frozen=True)
class FaceToFaceStep:
    """One accepted adjacent exact-face event and independent rank audit."""

    before: fi.FaceArrivalCommittedState
    after: fi.FaceArrivalCommittedState
    assembly: ff.FaceToFaceAssembly
    boundary: ct.PoreBoundary
    solver_controls: ci.CutSolverControls
    event_controls: fi.FaceEventControls
    seed: FaceToFaceSeed
    ledger: fi.FaceEventLedger
    rank_audit: ff.FaceToFaceRankAudit
    physically_qualifying: bool = field(default=False, init=False)


_STEFAN_SCALE_FLOOR_MOL_M2_S = 1.0e-3


def solve_adjacent_face_event(
    before: fi.FaceArrivalCommittedState,
    solver_controls: ci.CutSolverControls,
    event_controls: fi.FaceEventControls,
    boundary: ct.PoreBoundary,
    seed: FaceToFaceSeed,
) -> FaceToFaceStep:
    """Solve and certify one whole-cell primary-drainage event or roll back."""

    evaluations = 0
    rejected_trials = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    last_candidate = None
    last_scaled: tuple[float, ...] = ()
    try:
        if not isinstance(before, fi.FaceArrivalCommittedState):
            raise TypeError("face-to-face before state must be exact-face committed")
        if not isinstance(solver_controls, ci.CutSolverControls):
            raise TypeError("face-to-face solve requires CutSolverControls")
        if solver_controls != before.solver_controls:
            raise ValueError(
                "face-to-face solve must preserve the exact arrival solver controls"
            )
        ci.certify_explicit_wet_retained_water_applicability(
            solver_controls,
            before.wet_retained_water_loadings,
            (0.0,) * len(before.wet_retained_water_loadings),
            before.config.wet.luikov,
        )
        if not isinstance(event_controls, fi.FaceEventControls):
            raise TypeError("face-to-face solve requires FaceEventControls")
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("face-to-face solve requires a supported pore boundary")
        layout = ff.layout_for_face_to_face(before)
        _validate_chart(before, solver_controls, event_controls)
        chart = _coordinate_chart(before, solver_controls, event_controls)
        _validate_seed(layout, chart, seed)
        stefan_scales = seed.stefan_flux_scales_mol_m2_s
        initial_coordinates = _encode_candidate(seed.candidate, chart, stefan_scales)
        initial = ff.assemble_face_to_face(
            before,
            seed.candidate,
            boundary,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        scales = _residual_scales(before, initial)
        scale_vector = np.asarray(scales.vector, dtype=float)
        structural = (
            fi._event_sparsity(layout)  # noqa: SLF001
            if layout.unknown_count > solver_controls.maximum_dense_unknowns
            else None
        )

        def scaled_residual(coordinates: np.ndarray) -> np.ndarray:
            nonlocal evaluations, rejected_trials
            candidate = _decode_candidate(coordinates, layout, chart, stefan_scales)
            try:
                assembly = ff.assemble_face_to_face(
                    before,
                    candidate,
                    boundary,
                    enforce_reduced_film_thresholds=False,
                )
            except ff.FaceToFaceStepError:
                rejected_trials += 1
                raise
            evaluations += 1
            return np.asarray(
                ff.datum_covariant_residual_vector(assembly.residuals),
                dtype=float,
            ) / scale_vector

        lower, upper = ci._coordinate_solver_bounds(  # noqa: SLF001
            layout,
            solver_controls.maximum_logit_magnitude,
        )
        solution = optimize.least_squares(
            scaled_residual,
            initial_coordinates,
            method="dogbox" if structural is not None else "trf",
            bounds=(lower, upper),
            ftol=solver_controls.nonlinear_step_tolerance,
            xtol=solver_controls.nonlinear_step_tolerance,
            gtol=solver_controls.nonlinear_step_tolerance,
            max_nfev=solver_controls.maximum_function_evaluations,
            x_scale="jac",
            jac_sparsity=None if structural is None else structural.matrix,
            tr_options=(
                {"atol": 1.0e-14, "btol": 1.0e-14, "maxiter": 1000}
                if structural is not None
                else None
            ),
        )
        candidate = _decode_candidate(solution.x, layout, chart, stefan_scales)
        converged = ff.assemble_face_to_face(before, candidate, boundary)
        evaluations += 1
        signed_scaled = tuple(
            residual / scale
            for residual, scale in zip(
                ff.datum_covariant_residual_vector(converged.residuals),
                scales.vector,
            )
        )
        last_candidate = candidate
        last_scaled = signed_scaled
        maximum_scaled = max(abs(value) for value in signed_scaled)
        condition_proxy = ci._jacobian_condition_proxy(solution.jac)  # noqa: SLF001
        if not solution.success:
            raise RuntimeError(
                "face-to-face nonlinear root did not certify success: "
                f"{solution.message}"
            )
        if np.any(solution.active_mask != 0):
            raise RuntimeError(
                "face-to-face solution touched a coordinate guard; open root failed"
            )
        if maximum_scaled > solver_controls.nonlinear_residual_tolerance:
            raise RuntimeError(
                "face-to-face regional/RH residual contract failed: "
                f"{maximum_scaled:.3e} > "
                f"{solver_controls.nonlinear_residual_tolerance:.3e}"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > solver_controls.maximum_condition_proxy
        ):
            raise RuntimeError(
                "face-to-face coordinate Jacobian condition is uncertified: "
                f"{condition_proxy:.3e}"
            )
        fractions = fi._bounded_fractions(candidate, chart)  # noqa: SLF001
        if not fractions:
            raise RuntimeError("face-to-face root lost its bounded coordinates")
        labels, variable_scales = _physical_variable_contract(
            before,
            layout,
            chart,
        )
        rank_audit = ff.audit_face_to_face_rank(
            converged,
            scales.vector,
            labels,
            variable_scales,
        )
        if not rank_audit.structurally_square or not rank_audit.numerically_full_rank:
            raise conditioning.JacobianCertificateError(
                "face-to-face physical Jacobian is not full rank"
            )
        physical_condition = rank_audit.jacobian_certificate.condition_number_2
        if physical_condition > solver_controls.maximum_condition_proxy:
            raise conditioning.JacobianCertificateError(
                "face-to-face physical Jacobian exceeds the declared condition "
                f"ceiling: {physical_condition:.3e} > "
                f"{solver_controls.maximum_condition_proxy:.3e}"
            )
        step = _accept_event(
            before,
            solver_controls,
            event_controls,
            boundary,
            seed,
            converged,
            scales,
            signed_scaled,
            evaluations,
            rejected_trials,
            str(solution.message),
            condition_proxy,
            chart,
            structural,
            rank_audit,
        )
        ci.certify_explicit_wet_retained_water_applicability(
            solver_controls,
            before.wet_retained_water_loadings,
            (0.0,) * len(before.wet_retained_water_loadings),
            before.config.wet.luikov,
        )
        ci.certify_explicit_wet_retained_water_applicability(
            solver_controls,
            step.after.wet_retained_water_loadings,
            (0.0,) * len(step.after.wet_retained_water_loadings),
            step.after.config.wet.luikov,
        )
        return step
    except FaceToFaceIntegratorStepError:
        raise
    except Exception as exc:
        if not isinstance(before, fi.FaceArrivalCommittedState):
            raise
        applicability_error = (
            exc
            if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
            else None
        )
        raise FaceToFaceIntegratorStepError(
            f"exact face-to-face solve rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            last_candidate=last_candidate,
            scaled_residuals=last_scaled,
            retained_water_applicability_guard_failed=(
                applicability_error is not None
            ),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


def _coordinate_chart(
    before: fi.FaceArrivalCommittedState,
    controls: ci.CutSolverControls,
    event_controls: fi.FaceEventControls,
) -> fi._CoordinateChart:
    layout = ff.layout_for_face_to_face(before)
    dry_config = before.config.dry
    return fi._CoordinateChart(  # noqa: SLF001
        wet_temperature=controls.wet_temperature_bounds_k,
        wet_water=controls.wet_water_bounds,
        dry_temperature=dry_config.conditioned_temperature_domain.solver_bounds_k,
        dry_y_hexane=(
            dry_config.primitive_band.y_hexane_bounds
            if isinstance(dry_config.primitive_band, ct.JointPrimitiveBand)
            else None
        ),
        dry_config=dry_config,
        dry_pores=tuple(
            replace(dry_config.pore, w_o=before.oil_fraction_labels[index])
            for index in layout.dry_cell_indices
        ),
        event_time=event_controls.event_time_bounds_s,
        interface_temperature=controls.interface_temperature_bounds_k,
    )


def _encode_candidate(
    candidate: ff.FaceToFaceUnknowns,
    chart: fi._CoordinateChart,
    stefan_scales: Sequence[float],
) -> np.ndarray:
    return fi._encode_candidate(candidate, chart, stefan_scales)  # noqa: SLF001


def _decode_candidate(
    coordinates: Sequence[float],
    layout: ff.FaceToFaceLayout,
    chart: fi._CoordinateChart,
    stefan_scales: Sequence[float],
) -> ff.FaceToFaceUnknowns:
    decoded = fi._decode_candidate(  # noqa: SLF001
        coordinates,
        layout,
        chart,
        stefan_scales,
    )
    return ff.FaceToFaceUnknowns.from_vector(layout, decoded.vector())


def _validate_seed(
    layout: ff.FaceToFaceLayout,
    chart: fi._CoordinateChart,
    seed: FaceToFaceSeed,
) -> None:
    if len(seed.candidate.vector()) != layout.unknown_count:
        raise ValueError("face-to-face seed does not match the event rank")
    if len(seed.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
        raise ValueError("one explicit Stefan scale is required per endpoint dry face")
    _encode_candidate(seed.candidate, chart, seed.stefan_flux_scales_mol_m2_s)


def _validate_chart(
    before: fi.FaceArrivalCommittedState,
    controls: ci.CutSolverControls,
    event_controls: fi.FaceEventControls,
) -> None:
    layout = ff.layout_for_face_to_face(before)
    wet_t_lo, wet_t_hi = controls.wet_temperature_bounds_k
    wet = before.config.wet.wet
    if not wet.T_min <= wet_t_lo < wet_t_hi <= wet.T_max:
        raise ValueError("wet solver chart left the wet caloric authority")
    water_lo, water_hi = controls.wet_water_bounds
    if not before.config.wet.luikov.W_ref <= water_lo < water_hi <= (
        before.config.wet.luikov.W_cap
    ):
        raise ValueError("wet-water chart left the Luikov authority")
    dry_t_bounds = before.config.dry.conditioned_temperature_domain.solver_bounds_k
    for value in before.wet_temperatures_k:
        if not wet_t_lo < value < wet_t_hi:
            raise ValueError("start wet temperature is outside the open chart")
    for value in before.wet_retained_water_loadings:
        if not water_lo < value < water_hi:
            raise ValueError("start retained water is outside the open chart")
    for value in before.dry_temperatures_k:
        if not dry_t_bounds[0] < value < dry_t_bounds[1]:
            raise ValueError("start dry temperature is outside the open chart")
    for index, temperature, value in zip(
        range(before.arrival_face_index, before.geometry.master_grid.n),
        before.dry_temperatures_k,
        before.dry_y_hexane,
    ):
        cp.evaluate_equilibrium(
            temperature,
            before.config.dry.pressure_pa,
            value,
            replace(before.config.dry.pore, w_o=before.oil_fraction_labels[index]),
        )
    gamma_lo, gamma_hi = controls.interface_temperature_bounds_k
    pore_lo, pore_hi = before.config.dry.pore.temperature_bounds_k
    if not max(wet.T_min, pore_lo) <= gamma_lo < gamma_hi <= min(
        wet.T_max,
        pore_hi,
    ):
        raise ValueError("interface-temperature chart left shared authority")
    k = layout.arrival_face_index
    for temperature in (gamma_lo, 0.5 * math.fsum((gamma_lo, gamma_hi)), gamma_hi):
        # This is a domain audit only.  It adds no interface equation.
        from dtdc_simulator.core2.particle import cut_face_event as face

        face.evaluate_face_interface_state(
            temperature,
            before.config,
            before.historical_hexane_loadings[k],
            before.oil_fraction_labels[k],
            before.oil_fraction_labels[k],
        )
    if event_controls.event_time_bounds_s[0] < 0.0:
        raise ValueError("face-to-face time chart cannot start before zero")


def _residual_scales(
    before: fi.FaceArrivalCommittedState,
    initial: ff.FaceToFaceAssembly,
) -> fi.FaceEventResidualScales:
    """Freeze component-inventory and datum-covariant capacity row scales."""

    tau = initial.candidate.event_time_s
    grid = before.geometry.master_grid
    layout = initial.layout
    floor = 1.0e-300
    covariant = ff.datum_covariant_residual_vector(initial.residuals)
    wet_water: list[float] = []
    wet_energy: list[float] = []
    for position, (index, old_cell, new_cell) in enumerate(
        zip(layout.wet_cell_indices, initial.old_wet_cells, initial.wet_cells)
    ):
        volume = grid.volumes[index]
        wet_water.append(
            max(
                volume * old_cell.retained_water_concentration_mol_m3 / tau,
                volume * new_cell.retained_water_concentration_mol_m3 / tau,
                abs(initial.residuals.wet_water_mol_s[position]),
                floor,
            )
        )
        wet_energy.append(
            max(
                fi._wet_capacity_energy(volume, old_cell, before) / tau,  # noqa: SLF001
                fi._wet_capacity_energy(volume, new_cell, before) / tau,  # noqa: SLF001
                abs(initial.residuals.wet_energy_w[position]),
                floor,
            )
        )

    transformed = layout.transformed_cell_index
    transformed_volume = grid.volumes[transformed]
    old_transformed = initial.old_wet_cells[-1]
    new_transformed = initial.dry_cells[0]
    old_hexane = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[transformed]
        / hx.M
    )
    transformed_old_capacity = fi._wet_capacity_energy(  # noqa: SLF001
        transformed_volume,
        old_transformed,
        before,
    )
    # The first endpoint dry cell is the finite continuation of the newborn
    # gas-only trace.  Reuse the departure chart's phase-safe one-sided
    # capacity stencil; this changes scaling only, never the energy residual.
    new_capacity, new_capacity_audit = _dry_capacity_energy_with_audit(
        transformed_volume,
        new_transformed,
        before,
        transformed,
    )
    dry_capacity_audits = [new_capacity_audit]
    dry_water = [
        max(
            transformed_volume
            * old_transformed.retained_water_concentration_mol_m3
            / tau,
            transformed_volume
            * new_transformed.total_water_concentration_mol_m3
            / tau,
            abs(initial.residuals.dry_water_mol_s[0]),
            floor,
        )
    ]
    dry_hexane = [
        max(
            transformed_volume * old_hexane / tau,
            transformed_volume
            * new_transformed.total_hexane_concentration_mol_m3
            / tau,
            abs(initial.residuals.dry_hexane_mol_s[0]),
            floor,
        )
    ]
    dry_energy = [
        max(
            transformed_old_capacity / tau,
            new_capacity / tau,
            abs(
                math.fsum(
                    (
                        initial.residuals.dry_energy_w[0],
                        -hx.H_REF * initial.residuals.dry_hexane_mol_s[0],
                    )
                )
            ),
            floor,
        )
    ]
    for position, index in enumerate(layout.dry_cell_indices[1:], start=1):
        old_cell = initial.old_dry_cells[position - 1]
        new_cell = initial.dry_cells[position]
        volume = grid.volumes[index]
        dry_water.append(
            max(
                volume * old_cell.total_water_concentration_mol_m3 / tau,
                volume * new_cell.total_water_concentration_mol_m3 / tau,
                abs(initial.residuals.dry_water_mol_s[position]),
                floor,
            )
        )
        dry_hexane.append(
            max(
                volume * old_cell.total_hexane_concentration_mol_m3 / tau,
                volume * new_cell.total_hexane_concentration_mol_m3 / tau,
                abs(initial.residuals.dry_hexane_mol_s[position]),
                floor,
            )
        )
        row_old_capacity, old_capacity_audit = _dry_capacity_energy_with_audit(
            volume,
            old_cell,
            before,
            index,
        )
        row_new_capacity, new_capacity_audit = _dry_capacity_energy_with_audit(
            volume,
            new_cell,
            before,
            index,
        )
        dry_capacity_audits.extend((old_capacity_audit, new_capacity_audit))
        dry_energy.append(
            max(
                row_old_capacity / tau,
                row_new_capacity / tau,
                abs(
                    math.fsum(
                        (
                            initial.residuals.dry_energy_w[position],
                            -hx.H_REF
                            * initial.residuals.dry_hexane_mol_s[position],
                        )
                    )
                ),
                floor,
            )
        )
    return fi.FaceEventResidualScales(
        wet_water_mol_s=tuple(wet_water),
        wet_energy_w=tuple(wet_energy),
        dry_water_mol_s=tuple(dry_water),
        dry_hexane_mol_s=tuple(dry_hexane),
        dry_energy_w=tuple(dry_energy),
        rh_water_mol_s=max(
            transformed_volume
            * old_transformed.retained_water_concentration_mol_m3
            / tau,
            abs(initial.residuals.rh_water_mol_s),
            floor,
        ),
        rh_hexane_mol_s=max(
            transformed_volume * old_hexane / tau,
            abs(initial.residuals.rh_hexane_mol_s),
            floor,
        ),
        rh_energy_w=max(
            transformed_old_capacity / tau,
            abs(covariant[-1]),
            floor,
        ),
        dry_energy_capacity_stencil_audits=tuple(dry_capacity_audits),
    )


def _dry_capacity_energy_with_audit(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    state: fi.FaceArrivalCommittedState,
    cell_index: int,
) -> tuple[float, ci.DryEnergyCapacityStencilAudit]:
    capacity, audit = ci._dry_energy_capacity_with_audit(  # noqa: SLF001
        cell.temperature_k,
        cell.y_hexane,
        replace(state.config.dry.pore, w_o=state.oil_fraction_labels[cell_index]),
        state.config,
    )
    result = volume_m3 * capacity * max(abs(cell.temperature_k), 1.0)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("dry face-to-face capacity scale must be positive")
    return result, audit


def _physical_variable_contract(
    before: fi.FaceArrivalCommittedState,
    layout: ff.FaceToFaceLayout,
    chart: fi._CoordinateChart,
) -> tuple[tuple[str, ...], tuple[float, ...]]:
    labels = (
        *(f"wet_temperature[{i}]" for i in range(layout.wet_piece_count)),
        *(f"wet_retained_water[{i}]" for i in range(layout.wet_piece_count)),
        *(f"dry_temperature[{i}]" for i in range(layout.dry_piece_count)),
        *(f"dry_y_hexane[{i}]" for i in range(layout.dry_piece_count)),
        *(f"total_stefan_flux[{i}]" for i in range(layout.dry_face_count)),
        "event_time",
        "interface_temperature",
    )
    start_fluxes = before.dry_total_stefan_fluxes_mol_m2_s
    stefan = tuple(
        max(
            abs(start_fluxes[min(index, len(start_fluxes) - 1)]),
            _STEFAN_SCALE_FLOOR_MOL_M2_S,
        )
        for index in range(layout.dry_face_count)
    )
    scales = (
        *((chart.wet_temperature[1] - chart.wet_temperature[0],) * layout.wet_piece_count),
        *((chart.wet_water[1] - chart.wet_water[0],) * layout.wet_piece_count),
        *((chart.dry_temperature[1] - chart.dry_temperature[0],) * layout.dry_piece_count),
        *((1.0,) * layout.dry_piece_count),
        *stefan,
        chart.event_time[1] - chart.event_time[0],
        chart.interface_temperature[1] - chart.interface_temperature[0],
    )
    if len(labels) != layout.unknown_count or len(scales) != layout.unknown_count:
        raise RuntimeError("face-to-face physical scale contract lost rank")
    return tuple(labels), tuple(float(value) for value in scales)


def _event_inventory(assembly: ff.FaceToFaceAssembly) -> ci.MaterialInventorySnapshot:
    grid = assembly.event_geometry.master_grid
    water = [0.0] * grid.n
    hexane = [0.0] * grid.n
    energy = [0.0] * grid.n
    for index, cell in zip(assembly.layout.wet_cell_indices, assembly.wet_cells):
        volume = grid.volumes[index]
        water[index] = volume * cell.retained_water_concentration_mol_m3
        hexane[index] = (
            volume
            * assembly.before.config.wet.wet.rho_dm_p
            * assembly.before.historical_hexane_loadings[index]
            / hx.M
        )
        energy[index] = volume * cell.energy_density_j_m3
    for index, cell in zip(assembly.layout.dry_cell_indices, assembly.dry_cells):
        volume = grid.volumes[index]
        water[index] = volume * cell.total_water_concentration_mol_m3
        hexane[index] = volume * cell.total_hexane_concentration_mol_m3
        energy[index] = volume * cell.energy_density_j_m3
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _material_changes(assembly: ff.FaceToFaceAssembly) -> ci.MaterialInventorySnapshot:
    n = assembly.event_geometry.master_grid.n
    water = [0.0] * n
    hexane = [0.0] * n
    energy = [0.0] * n
    for position, index in enumerate(assembly.layout.wet_cell_indices):
        water[index] = assembly.inventory_changes.wet_water_mol[position]
        energy[index] = assembly.inventory_changes.wet_energy_j[position]
    for position, index in enumerate(assembly.layout.dry_cell_indices):
        water[index] = assembly.inventory_changes.dry_water_mol[position]
        hexane[index] = assembly.inventory_changes.dry_hexane_mol[position]
        energy[index] = assembly.inventory_changes.dry_energy_j[position]
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _capacity_scale(assembly: ff.FaceToFaceAssembly) -> float:
    grid = assembly.event_geometry.master_grid
    values = [
        fi._wet_capacity_energy(grid.volumes[index], cell, assembly.before)  # noqa: SLF001
        for index, cell in zip(assembly.layout.wet_cell_indices, assembly.wet_cells)
    ]
    values.extend(
        di._dry_capacity_energy(  # noqa: SLF001
            grid.volumes[index],
            cell,
            assembly.before,
            index,
        )
        for index, cell in zip(assembly.layout.dry_cell_indices, assembly.dry_cells)
    )
    result = math.fsum(values)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("face-to-face capacity scale must be positive")
    return result


def _accept_event(
    before: fi.FaceArrivalCommittedState,
    solver_controls: ci.CutSolverControls,
    event_controls: fi.FaceEventControls,
    boundary: ct.PoreBoundary,
    seed: FaceToFaceSeed,
    assembly: ff.FaceToFaceAssembly,
    scales: fi.FaceEventResidualScales,
    scaled_residuals: tuple[float, ...],
    evaluations: int,
    rejected_trials: int,
    nonlinear_message: str,
    condition_proxy: float,
    chart: fi._CoordinateChart,
    structural,
    rank_audit: ff.FaceToFaceRankAudit,
) -> FaceToFaceStep:
    tau = assembly.candidate.event_time_s
    current = _event_inventory(assembly)
    changes = _material_changes(assembly)
    area = assembly.event_geometry.master_grid.areas[-1]
    surface = assembly.dry_face_fluxes[-1]
    water_out = tau * area * surface.component.conserved_water_flux_mol_m2_s
    hexane_out = tau * area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_out = tau * area * surface.energy.total_energy_flux_w_m2
    water_step = math.fsum((*changes.water_cell_mol, water_out))
    hexane_step = math.fsum((*changes.hexane_cell_mol, hexane_out))
    energy_step = math.fsum((*changes.energy_cell_j, energy_out))
    after_capacity = _capacity_scale(assembly)
    step_water_scale = max(
        before.current_inventory.total_water_mol,
        current.total_water_mol,
        abs(water_out),
        1.0e-300,
    )
    step_hexane_scale = max(
        before.current_inventory.total_hexane_mol,
        current.total_hexane_mol,
        abs(hexane_out),
        1.0e-300,
    )
    step_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        after_capacity,
        1.0e-300,
    )
    boundary_water, boundary_water_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_water_out_mol,
        before.boundary_water_compensation_mol,
        water_out,
    )
    boundary_hexane, boundary_hexane_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_hexane_out_mol,
        before.boundary_hexane_compensation_mol,
        hexane_out,
    )
    boundary_energy, boundary_energy_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_energy_out_j,
        before.boundary_energy_compensation_j,
        energy_out,
    )
    cumulative_water, water_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_water_change_cell_mol,
        before.water_change_compensation_cell_mol,
        changes.water_cell_mol,
    )
    cumulative_hexane, hexane_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_hexane_change_cell_mol,
        before.hexane_change_compensation_cell_mol,
        changes.hexane_cell_mol,
    )
    cumulative_energy, energy_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_energy_change_cell_j,
        before.energy_change_compensation_cell_j,
        changes.energy_cell_j,
    )
    cumulative_water_residual = math.fsum(
        (*cumulative_water, *water_comp, boundary_water, boundary_water_comp)
    )
    cumulative_hexane_residual = math.fsum(
        (*cumulative_hexane, *hexane_comp, boundary_hexane, boundary_hexane_comp)
    )
    cumulative_energy_residual = math.fsum(
        (*cumulative_energy, *energy_comp, boundary_energy, boundary_energy_comp)
    )
    corrected_boundary_water = math.fsum((boundary_water, boundary_water_comp))
    corrected_boundary_hexane = math.fsum((boundary_hexane, boundary_hexane_comp))
    cumulative_water_scale = max(
        before.reference_inventory.total_water_mol,
        current.total_water_mol,
        abs(corrected_boundary_water),
        before.cumulative_absolute_water_transfer_mol + abs(water_out),
        1.0e-300,
    )
    cumulative_hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        current.total_hexane_mol,
        abs(corrected_boundary_hexane),
        before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out),
        1.0e-300,
    )
    cumulative_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        after_capacity,
        before.cumulative_absolute_energy_transfer_j + abs(energy_out),
        1.0e-300,
    )
    normalized = (
        abs(water_step) / step_water_scale,
        abs(hexane_step) / step_hexane_scale,
        abs(energy_step) / step_energy_scale,
        abs(cumulative_water_residual) / cumulative_water_scale,
        abs(cumulative_hexane_residual) / cumulative_hexane_scale,
        abs(cumulative_energy_residual) / cumulative_energy_scale,
    )
    if max(normalized) > solver_controls.ledger_tolerance:
        raise RuntimeError(
            "face-to-face signed ledger failed before commit: "
            f"{max(normalized):.3e} > {solver_controls.ledger_tolerance:.3e}"
        )
    candidate = assembly.candidate
    after = fi.FaceArrivalCommittedState(
        geometry=assembly.event_geometry,
        config=before.config,
        solver_controls=solver_controls,
        time_s=before.time_s + tau,
        arrival_face_index=assembly.layout.arrival_face_index,
        wet_temperatures_k=candidate.wet_temperatures_k,
        wet_retained_water_loadings=candidate.wet_retained_water_loadings,
        dry_temperatures_k=candidate.dry_temperatures_k,
        dry_y_hexane=candidate.dry_y_hexane,
        historical_hexane_loadings=before.historical_hexane_loadings,
        oil_fraction_labels=before.oil_fraction_labels,
        interface_temperature_k=candidate.interface_temperature_k,
        dry_total_stefan_fluxes_mol_m2_s=(
            candidate.dry_total_stefan_fluxes_mol_m2_s
        ),
        current_inventory=current,
        reference_inventory=before.reference_inventory,
        reference_capacity_energy_scale_j=before.reference_capacity_energy_scale_j,
        cumulative_boundary_water_out_mol=boundary_water,
        cumulative_boundary_hexane_out_mol=boundary_hexane,
        cumulative_boundary_energy_out_j=boundary_energy,
        boundary_water_compensation_mol=boundary_water_comp,
        boundary_hexane_compensation_mol=boundary_hexane_comp,
        boundary_energy_compensation_j=boundary_energy_comp,
        cumulative_absolute_water_transfer_mol=(
            before.cumulative_absolute_water_transfer_mol + abs(water_out)
        ),
        cumulative_absolute_hexane_transfer_mol=(
            before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out)
        ),
        cumulative_absolute_energy_transfer_j=(
            before.cumulative_absolute_energy_transfer_j + abs(energy_out)
        ),
        cumulative_material_water_change_cell_mol=cumulative_water,
        cumulative_material_hexane_change_cell_mol=cumulative_hexane,
        cumulative_material_energy_change_cell_j=cumulative_energy,
        water_change_compensation_cell_mol=water_comp,
        hexane_change_compensation_cell_mol=hexane_comp,
        energy_change_compensation_cell_j=energy_comp,
        accepted_steps=before.accepted_steps + 1,
        cumulative_nonlinear_evaluations=(
            before.cumulative_nonlinear_evaluations + evaluations
        ),
    )
    fractions = fi._bounded_fractions(candidate, chart)  # noqa: SLF001
    ledger = fi.FaceEventLedger(
        event_time_s=tau,
        nonlinear_evaluations=evaluations,
        rejected_trial_evaluations=rejected_trials,
        nonlinear_message=nonlinear_message,
        residual_scales=scales,
        scaled_residuals=scaled_residuals,
        maximum_scaled_residual=max(abs(value) for value in scaled_residuals),
        condition_proxy=condition_proxy,
        minimum_fractional_distance_to_bound=min(
            min(value, 1.0 - value) for value in fractions
        ),
        used_sparse_jacobian=structural is not None,
        scalability_class=(
            "bordered-banded sparse adjacent-face event Jacobian"
            if structural is not None
            else "certified small-N dense adjacent-face event oracle"
        ),
        jacobian_structural_color_count=(
            assembly.layout.unknown_count
            if structural is None
            else structural.color_count
        ),
        jacobian_structural_nnz=(
            assembly.layout.unknown_count**2 if structural is None else structural.nnz
        ),
        water_step_residual_mol=water_step,
        hexane_step_residual_mol=hexane_step,
        energy_step_residual_j=energy_step,
        normalized_water_step_residual=normalized[0],
        normalized_hexane_step_residual=normalized[1],
        normalized_energy_step_residual=normalized[2],
        water_cumulative_residual_mol=cumulative_water_residual,
        hexane_cumulative_residual_mol=cumulative_hexane_residual,
        energy_cumulative_residual_j=cumulative_energy_residual,
        normalized_water_cumulative_residual=normalized[3],
        normalized_hexane_cumulative_residual=normalized[4],
        normalized_energy_cumulative_residual=normalized[5],
        boundary_water_out_mol=water_out,
        boundary_hexane_out_mol=hexane_out,
        boundary_energy_out_j=energy_out,
        accepted=True,
    )
    return FaceToFaceStep(
        before,
        after,
        assembly,
        boundary,
        solver_controls,
        event_controls,
        seed,
        ledger,
        rank_audit,
    )


__all__ = [
    "FaceToFaceIntegratorStepError",
    "FaceToFaceSeed",
    "FaceToFaceStep",
    "solve_adjacent_face_event",
]
