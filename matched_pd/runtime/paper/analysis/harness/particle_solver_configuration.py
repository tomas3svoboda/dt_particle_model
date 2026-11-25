"""Explicit Paper 1 property configuration and face-search policy.

No source inspection, rewriting or installation is needed for these numerical
and material selections. Build the initial feed with ``wet_parameters`` and
use ``cut_configuration`` when constructing its resolved wet/dry model.
"""

from dataclasses import dataclass, replace
import math
import numpy as np

from dtdc_simulator.core2.particle import cut_face_departure_integrator as departure
from dtdc_simulator.core2.particle import cut_geometry as geometry
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.props import sorption as sp

PAPER_CHECKPOINT_CODEC = "paper1-explicit-properties-stable-controls-v2"


def radius_normalized_birth_jacobian(physical_jacobian, wet_piece_count, radius):
    """Jacobian for (g_T*R,g_y*R), retaining every other physical coordinate."""
    matrix = np.asarray(physical_jacobian, dtype=float)
    if not math.isfinite(radius) or radius <= 0 or wet_piece_count <= 0:
        raise ValueError("positive native radius and wet-piece count required")
    if (
        matrix.ndim != 2
        or matrix.shape[0] != matrix.shape[1]
        or matrix.shape[1] <= 2 * wet_piece_count + 1
    ):
        raise ValueError("physical birth Jacobian does not contain both gradient columns")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("physical birth Jacobian must be finite")
    normalized = matrix.copy()
    normalized[:, [2 * wet_piece_count, 2 * wet_piece_count + 1]] /= radius
    return normalized


@dataclass(frozen=True)
class ParticleModelParameters:
    binary_reference_m2_s: float = 4e-10
    oil_law: sp.OilHeatCapacityLaw = sp.OilHeatCapacityLaw.GARCIA_ROJAS_2013
    binary_law: tc.BinaryDiffusivityLaw = tc.BinaryDiffusivityLaw.FULLER_EFFECTIVE_REFERENCE

    def __post_init__(self):
        if not math.isfinite(self.binary_reference_m2_s) or self.binary_reference_m2_s <= 0:
            raise ValueError("binary pore reference must be positive and finite")
        if not isinstance(self.oil_law, sp.OilHeatCapacityLaw) or not isinstance(
            self.binary_law, tc.BinaryDiffusivityLaw
        ):
            raise TypeError("particle property laws require explicit supported identities")

    def wet_parameters(self, before):
        return replace(before, oil_heat_capacity_law=self.oil_law)

    def cut_configuration(self, before):
        interval = tc.EffectiveBinaryDiffusivityInterval(
            self.binary_reference_m2_s,
            self.binary_reference_m2_s,
            "Literature-informed effective-pore assumption at 378.15 K, 101325 Pa; not identified binary measurements",
        )
        binary = replace(interval.select(0), temperature_law=self.binary_law)
        return replace(
            before,
            wet=replace(before.wet, wet=self.wet_parameters(before.wet.wet)),
            dry=replace(
                before.dry,
                pore=replace(before.dry.pore, oil_heat_capacity_law=self.oil_law),
                binary_diffusivity=binary,
            ),
            interface_composition=replace(
                before.interface_composition, finish_at_representable_bracket=True
            ),
        )


PAPER_FACE_SEARCH = departure.FaceDepartureNumericalPolicy(
    feasible_seed_reference=True,
    root_search_tolerance=1e-11,
    cauchy_globalization=True,
    powell_dogleg_rescue=True,
    require_inward_phase_direction=True,
)


def warming_face_seed(tangent, dt_s):
    """Existing independently admissible uncommitted seed, made explicit.

    From the October 16 recovery's warming constant-composition seed. This is
    a finite-root guess, not a certified asymptotic tangent or accepted step.
    Native seed, phase, root, conditioning and ledger checks own acceptance.
    """
    if not math.isfinite(dt_s) or dt_s <= 0:
        raise ValueError("face-seed duration must be positive and finite")
    before, t = tangent.before, tangent.candidate
    radius = math.fsum((before.geometry.front.radius_m, dt_s * t.front_speed_m_s))
    grid = before.geometry.master_grid
    partition = geometry.partition_master_grid(grid, radius_m=radius)
    partition = geometry.partition_master_grid(grid, z=partition.front.z)
    distance = departure._newborn_center_distance(partition, before.arrival_face_index - 1)
    gradients = departure.departure.NewbornDryTraceGradients(t.dry_temperature_gradient_k_m, 0.0)
    fluxes = (
        tangent.direct_interface_flux.recovered_total_stefan_flux_mol_m2_s,
        t.outer_total_stefan_flux_mol_m2_s,
        *before.dry_total_stefan_fluxes_mol_m2_s[1:],
    )
    candidate = departure.departure.FaceDepartureUnknowns(
        before.wet_temperatures_k,
        before.wet_retained_water_loadings,
        (
            math.fsum((t.interface_temperature_k, gradients.temperature_gradient_k_m * distance)),
            *before.dry_temperatures_k,
        ),
        (tangent.interface.y_hexane, *before.dry_y_hexane),
        fluxes,
        partition.front.z,
        t.interface_temperature_k,
    )
    return departure.FaceDepartureSeed(
        candidate,
        tuple(max(abs(v), 0.01) for v in fluxes),
        "Uncommitted warming constant-composition guess; no physical state projection",
        gradients,
    )
