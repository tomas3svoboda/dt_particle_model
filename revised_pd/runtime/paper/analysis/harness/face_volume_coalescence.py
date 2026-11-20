"""Conservative local mesh coalescence at an exact particle-front face.

This candidate spatial discretization removes the master face coincident with
the physical front. The old wet and dry phase pieces remain separately stored
inside one ordinary native cut cell, each with finite unchanged volume. No
zero-volume newborn dry state is introduced. The native ALE/RH operator and
the original constitutive/phase laws then own the next physical advance.

This is an engineering experiment, not a qualified replacement mesh policy.
No time advances, no ledger resets, and no parameter fitting occur here.
"""

from dataclasses import dataclass, fields, replace
import math

from dtdc_simulator.core2.particle import cut_face_event_integrator as face
from dtdc_simulator.core2.particle import cut_geometry as geometry
from dtdc_simulator.core2.particle import cut_integrator as integrator
from dtdc_simulator.core2.particle import cut_transport as transport
from dtdc_simulator.core2.particle import grid
from dtdc_simulator.core2.particle import coupled_transport as dry_transport


class FaceCoalescenceError(RuntimeError):
    """A refused mesh operation carrying the exact original physical state."""

    def __init__(self, message, before):
        super().__init__(message)
        self.rollback_state = before


def restore_fully_dry_face(before, missing_radius_m):
    """Restore the innermost base face after the native dry projection.

    Both child cells keep the identical native T/y/material state. Scalar
    reference inventories and all compensated cumulative histories survive.
    """
    try:
        if not isinstance(before, dry_transport.FullyDryTransportState):
            raise TypeError("face restoration requires the native fully dry state")
        old_grid = before.grid
        if (
            not math.isfinite(missing_radius_m)
            or not old_grid.faces[0] < missing_radius_m < old_grid.faces[1]
        ):
            raise ValueError("missing face must split the first finite dry cell")
        new_grid = grid._from_faces(
            (old_grid.faces[0], missing_radius_m, *old_grid.faces[1:]), old_grid.R
        )
        after = replace(
            before,
            grid=new_grid,
            temperatures_k=(before.temperatures_k[0], *before.temperatures_k),
            y_hexane=(before.y_hexane[0], *before.y_hexane),
        )

        def inventory(state):
            return dry_transport._inventories(
                state.grid,
                dry_transport._equilibrium_states(
                    state.temperatures_k, state.y_hexane, state.config
                ),
            )

        old, new = inventory(before), inventory(after)
        old_values = _sum_inventory(old)
        new_values = _sum_inventory(new)
        energy_scale = math.fsum(
            dry_transport._equilibrium_energy_capacity(t, y, before.config) * v * max(abs(t), 1.0)
            for t, y, v in zip(before.temperatures_k, before.y_hexane, old_grid.volumes)
        )
        scales = (
            max(abs(old_values[0]), abs(new_values[0])),
            max(abs(old_values[1]), abs(new_values[1])),
            energy_scale,
        )
        normalized = tuple(
            abs(math.fsum((a, -b))) / s for a, b, s in zip(new_values, old_values, scales)
        )
        if max(normalized) > before.config.ledger_tolerance:
            raise ValueError("native conservation refused dry base-face restoration")
        # Every other field is copied by dataclasses.replace; check the full
        # declared scalar ledger/config/reference set, not selected totals.
        preserved = all(
            getattr(after, f.name) == getattr(before, f.name)
            for f in fields(before)
            if f.init and f.name not in ("grid", "temperatures_k", "y_hexane")
        )
        if not preserved:
            raise ValueError("dry mesh operation changed its native reference/history")
        return after, {
            "restored_face_radius_m": missing_radius_m,
            "old_master_cells": old_grid.n,
            "new_master_cells": new_grid.n,
            "old_W_H_U": old_values,
            "new_W_H_U": new_values,
            "normalized_inventory_differences": normalized,
            "all_reference_and_history_fields_preserved": preserved,
            "physically_qualifying": False,
        }
    except FaceCoalescenceError:
        raise
    except Exception as error:
        raise FaceCoalescenceError(str(error), before) from error


@dataclass(frozen=True)
class FaceCoalescenceAudit:
    removed_face_radius_m: float
    old_master_cells: int
    new_master_cells: int
    old_water_hexane_energy: tuple[float, float, float]
    new_water_hexane_energy: tuple[float, float, float]
    inventory_differences: tuple[float, float, float]
    normalized_inventory_differences: tuple[float, float, float]
    cumulative_ledger_differences: tuple[float, float, float]
    native_reference_inventory_identity_preserved: bool
    native_controls_identity_preserved: bool
    physical_time_and_step_count_preserved: bool
    old_phase_piece_volumes_m3: tuple[float, float]
    new_phase_piece_volumes_m3: tuple[float, float]
    physically_qualifying: bool = False
    restored_face_radius_m: float | None = None
    first_dry_volume_split_fraction: float | None = None
    inserted_stefan_flux_is_uncommitted_guess: bool = False


def _sum_inventory(snapshot):
    return snapshot.total_water_mol, snapshot.total_hexane_mol, snapshot.total_energy_j


def _merge_compensated_pair(high, compensation, index):
    value, remainder = integrator._compensated_add_scalar(
        high[index], compensation[index], high[index + 1]
    )
    value, remainder = integrator._compensated_add_scalar(value, remainder, compensation[index + 1])
    return (
        (*high[:index], value, *high[index + 2 :]),
        (*compensation[:index], remainder, *compensation[index + 2 :]),
    )


def coalesce_arrival_face(before):
    """Keep the two original phase states inside a finite native cut cell.

    Merging is permitted only for identical exact material/history labels.
    Cumulative material changes are conservatively summed with their native
    compensation pairs. The original reference object, source/boundary totals,
    capacities, counters and physical time are retained, including any existing
    conservation defect. The original ledger tolerance is not changed.
    """
    try:
        if not isinstance(before, face.FaceArrivalCommittedState):
            raise TypeError("coalescence requires an exact native face-arrival state")
        old_grid = before.geometry.master_grid
        k = before.arrival_face_index
        for labels in (before.historical_hexane_loadings, before.oil_fraction_labels):
            if labels[k - 1] != labels[k]:
                raise ValueError("coalescence requires identical adjacent material/history labels")
        new_grid = grid._from_faces((*old_grid.faces[:k], *old_grid.faces[k + 1 :]), old_grid.R)
        partition = geometry.partition_master_grid(
            new_grid, radius_m=before.geometry.front.radius_m
        )
        if partition.front.z != before.geometry.front.z:
            raise ValueError("mesh operation changed the persisted front coordinate")
        cut = transport.initialize_cut_state(
            partition,
            before.config,
            before.wet_temperatures_k,
            before.wet_retained_water_loadings,
            before.dry_temperatures_k,
            before.dry_y_hexane,
            (*before.historical_hexane_loadings[:k], *before.historical_hexane_loadings[k + 1 :]),
            oil_fraction_labels=(
                *before.oil_fraction_labels[:k],
                *before.oil_fraction_labels[k + 1 :],
            ),
            time_s=before.time_s,
        )
        kwargs = {
            f.name: getattr(before, f.name)
            for f in fields(integrator.CutIntegratorState)
            if f.init and hasattr(before, f.name)
        }
        kwargs.update(
            transport=cut,
            controls=before.solver_controls,
            last_interface_temperature_k=before.interface_temperature_k,
            last_total_stefan_fluxes_mol_m2_s=before.dry_total_stefan_fluxes_mol_m2_s,
        )
        cumulative_differences = []
        for quantity, unit in (("water", "mol"), ("hexane", "mol"), ("energy", "j")):
            high_name = f"cumulative_material_{quantity}_change_cell_{unit}"
            comp_name = f"{quantity}_change_compensation_cell_{unit}"
            old_high, old_comp = getattr(before, high_name), getattr(before, comp_name)
            new_high, new_comp = _merge_compensated_pair(old_high, old_comp, k - 1)
            kwargs.update({high_name: new_high, comp_name: new_comp})
            cumulative_differences.append(
                math.fsum((*new_high, *new_comp, *(-v for v in old_high), *(-v for v in old_comp)))
            )
        after = integrator.CutIntegratorState(**kwargs)
        old_values = _sum_inventory(before.current_inventory)
        new_values = _sum_inventory(integrator.inventory_snapshot(after.transport))
        differences = tuple(math.fsum((n, -o)) for o, n in zip(old_values, new_values))
        scales = (
            max(abs(old_values[0]), abs(new_values[0])),
            max(abs(old_values[1]), abs(new_values[1])),
            before.reference_capacity_energy_scale_j,
        )
        normalized = tuple(abs(d) / s for d, s in zip(differences, scales))
        if max(normalized) > before.solver_controls.ledger_tolerance:
            raise ValueError("native component/energy conservation refused mesh coalescence")
        if (
            max(abs(d) / s for d, s in zip(cumulative_differences, scales))
            > before.solver_controls.ledger_tolerance
        ):
            raise ValueError("native cumulative ledger refused mesh coalescence")
        original_volumes = (
            before.geometry.cells[k - 1].wet_volume_m3,
            before.geometry.cells[k].dry_volume_m3,
        )
        new_volumes = (partition.cells[k - 1].wet_volume_m3, partition.cells[k - 1].dry_volume_m3)
        if (
            max(abs(a - b) for a, b in zip(original_volumes, new_volumes)) / old_grid.total_volume
            > before.solver_controls.ledger_tolerance
        ):
            raise ValueError("native geometric conservation refused mesh coalescence")
        audit = FaceCoalescenceAudit(
            old_grid.faces[k],
            old_grid.n,
            new_grid.n,
            old_values,
            new_values,
            differences,
            normalized,
            tuple(cumulative_differences),
            after.reference_inventory is before.reference_inventory,
            after.controls is before.solver_controls,
            after.transport.time_s == before.time_s
            and after.accepted_steps == before.accepted_steps,
            original_volumes,
            new_volumes,
        )
        return after, audit
    except FaceCoalescenceError:
        raise
    except Exception as error:
        raise FaceCoalescenceError(str(error), before) from error


def _split_compensated_pair(high, compensation, fraction):
    """Split an extensive history without discarding its low-order part."""
    part = high * fraction
    correction = math.fsum((math.fma(high, fraction, -part), compensation * fraction))
    remaining, remainder = integrator._compensated_add_scalar(high, compensation, -part)
    remaining, remainder = integrator._compensated_add_scalar(remaining, remainder, -correction)
    return (part, correction), (remaining, remainder)


def slide_arrival_face(before, previously_removed_radius_m):
    """Restore the prior dry face while coalescing the current front face.

    Exactly one base face remains absent. The fully dry first cell is split
    with unchanged intensive states; its inner part and the adjacent full wet
    cell become the next native cut cell. This avoids both a zero-volume dry
    newborn and cumulative loss of radial resolution. The new internal Stefan
    flux is only a conservative interpolated algebraic seed for the next solve.
    """
    try:
        if not isinstance(before, face.FaceArrivalCommittedState):
            raise TypeError("sliding requires an exact native face-arrival state")
        old_grid = before.geometry.master_grid
        k = before.arrival_face_index
        restored = previously_removed_radius_m
        if not math.isfinite(restored) or not old_grid.faces[k] < restored < old_grid.faces[k + 1]:
            raise ValueError("previously removed face must lie inside the first full dry cell")
        for labels in (before.historical_hexane_loadings, before.oil_fraction_labels):
            if labels[k - 1] != labels[k]:
                raise ValueError("sliding requires identical adjacent material/history labels")
        new_grid = grid._from_faces(
            (*old_grid.faces[:k], restored, *old_grid.faces[k + 1 :]), old_grid.R
        )
        partition = geometry.partition_master_grid(new_grid, radius_m=old_grid.faces[k])
        if partition.front.z != before.geometry.front.z:
            raise ValueError("mesh operation changed the persisted front coordinate")
        fraction = partition.cells[k - 1].dry_volume_m3 / before.geometry.cells[k].dry_volume_m3
        if not 0.0 < fraction < 1.0:
            raise ValueError("sliding requires two finite dry pieces")
        cut = transport.initialize_cut_state(
            partition,
            before.config,
            before.wet_temperatures_k,
            before.wet_retained_water_loadings,
            (before.dry_temperatures_k[0], *before.dry_temperatures_k),
            (before.dry_y_hexane[0], *before.dry_y_hexane),
            before.historical_hexane_loadings,
            oil_fraction_labels=before.oil_fraction_labels,
            time_s=before.time_s,
        )
        # Interpolate the extensive area*flux in spherical volume. This field
        # supplies an initial algebraic guess, never a committed flux integral.
        flux = before.dry_total_stefan_fluxes_mol_m2_s
        inserted = (
            math.fsum(
                (
                    (1.0 - fraction) * old_grid.areas[k] * flux[0],
                    fraction * old_grid.areas[k + 1] * flux[1],
                )
            )
            / new_grid.areas[k]
        )
        kwargs = {
            f.name: getattr(before, f.name)
            for f in fields(integrator.CutIntegratorState)
            if f.init and hasattr(before, f.name)
        }
        kwargs.update(
            transport=cut,
            controls=before.solver_controls,
            last_interface_temperature_k=before.interface_temperature_k,
            last_total_stefan_fluxes_mol_m2_s=(flux[0], inserted, *flux[1:]),
        )
        cumulative_differences = []
        for quantity, unit in (("water", "mol"), ("hexane", "mol"), ("energy", "j")):
            high_name = f"cumulative_material_{quantity}_change_cell_{unit}"
            comp_name = f"{quantity}_change_compensation_cell_{unit}"
            old_high, old_comp = getattr(before, high_name), getattr(before, comp_name)
            part, remainder = _split_compensated_pair(old_high[k], old_comp[k], fraction)
            merged = integrator._compensated_add_scalar(old_high[k - 1], old_comp[k - 1], part[0])
            merged = integrator._compensated_add_scalar(*merged, part[1])
            new_high = (*old_high[: k - 1], merged[0], remainder[0], *old_high[k + 1 :])
            new_comp = (*old_comp[: k - 1], merged[1], remainder[1], *old_comp[k + 1 :])
            kwargs.update({high_name: new_high, comp_name: new_comp})
            cumulative_differences.append(
                math.fsum((*new_high, *new_comp, *(-v for v in old_high), *(-v for v in old_comp)))
            )
        after = integrator.CutIntegratorState(**kwargs)
        old_values = _sum_inventory(before.current_inventory)
        new_values = _sum_inventory(integrator.inventory_snapshot(cut))
        differences = tuple(math.fsum((n, -o)) for o, n in zip(old_values, new_values))
        scales = (
            max(abs(old_values[0]), abs(new_values[0])),
            max(abs(old_values[1]), abs(new_values[1])),
            before.reference_capacity_energy_scale_j,
        )
        normalized = tuple(abs(d) / s for d, s in zip(differences, scales))
        if max(normalized) > before.solver_controls.ledger_tolerance:
            raise ValueError("native component/energy conservation refused mesh sliding")
        if (
            max(abs(d) / s for d, s in zip(cumulative_differences, scales))
            > before.solver_controls.ledger_tolerance
        ):
            raise ValueError("native cumulative ledger refused mesh sliding")
        old_volumes = tuple(
            math.fsum(getattr(cell, phase + "_volume_m3") for cell in before.geometry.cells)
            for phase in ("wet", "dry")
        )
        new_volumes = tuple(
            math.fsum(getattr(cell, phase + "_volume_m3") for cell in partition.cells)
            for phase in ("wet", "dry")
        )
        if (
            max(abs(a - b) for a, b in zip(old_volumes, new_volumes)) / old_grid.total_volume
            > before.solver_controls.ledger_tolerance
        ):
            raise ValueError("native geometric conservation refused mesh sliding")
        return after, FaceCoalescenceAudit(
            old_grid.faces[k],
            old_grid.n,
            new_grid.n,
            old_values,
            new_values,
            differences,
            normalized,
            tuple(cumulative_differences),
            after.reference_inventory is before.reference_inventory,
            after.controls is before.solver_controls,
            cut.time_s == before.time_s and after.accepted_steps == before.accepted_steps,
            old_volumes,
            new_volumes,
            restored_face_radius_m=restored,
            first_dry_volume_split_fraction=fraction,
            inserted_stefan_flux_is_uncommitted_guess=True,
        )
    except FaceCoalescenceError:
        raise
    except Exception as error:
        raise FaceCoalescenceError(str(error), before) from error
