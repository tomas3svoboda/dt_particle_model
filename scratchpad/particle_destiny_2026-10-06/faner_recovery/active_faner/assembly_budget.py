"""Scoped hard accounting for the existing native radial/film study.

Every attempted native cut assembly counts, including initial assembly,
finite-difference probes, and calls that raise on invalid trial states. This
wrapper does not replace a residual, change a physical law, or accept a state.
"""

from contextlib import contextmanager
import time

import numpy as np
import finite_film_core2_study as study

core = study.core


class AssemblyBudgetExceeded(RuntimeError):
    classification = "ASSEMBLY_BUDGET_EXHAUSTED"


class AssemblyWallBudgetExceeded(RuntimeError):
    classification = "ASSEMBLY_WALL_BUDGET_EXHAUSTED"


@contextmanager
def assembly_budget(boundary, limit=400, remaining_wall_seconds=None):
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("assembly budget must be a positive integer")
    original_assemble = core.cut.assemble_backward_euler
    original_scales = core.ci._residual_scales
    started = time.perf_counter()
    audit = {
        "declared_limit": limit,
        "assemblies_started": 0,
        "assemblies_returned": 0,
        "invalid_trial_assemblies": 0,
        "budget_exhausted": False,
        "stage": "initial_surface_and_radial_assembly",
        "last_complete_scaled_residuals": None,
        "last_complete_maximum_scaled_residual": None,
        "last_native_raw_residuals": None,
        "last_assembly_error": None,
    }
    captured = {"scales": None, "last_assembly": None}
    mass_scale = max(boundary.film_boundary.correlated_pair.binary_mass_transfer_m_s, 1e-8)

    def record_complete(assembly):
        raw = np.asarray(core.ci._datum_covariant_residual_vector(assembly.residuals))
        audit["last_native_raw_residuals"] = raw.tolist()
        scales = captured["scales"]
        if scales is None:
            return
        native = raw / np.asarray(scales.vector)
        film = assembly.surface_film_audit
        # Exactly the existing joint study's two film-row scales, retaining
        # its seed outer-cell temperature as the fixed heat-scale reference.
        tc = captured["seed_outer_cell_temperature"]
        heat_scale = max(
            boundary.film_boundary.correlated_pair.heat_transfer_w_m2_k
            * abs(boundary.temperature_k - tc),
            1.0,
        )
        if hasattr(film, "water_node_residual_kg_m2_s"):
            complete = np.r_[
                native,
                film.water_node_residual_kg_m2_s / 0.05,
                film.hexane_node_residual_kg_m2_s / 0.05,
                film.energy_node_residual_w_m2 / 5000.0,
            ]
            audit["surface_row_scale_authority"] = (
                "absent-film conservative node; no liquid storage"
            )
            surface_temperature = film.temperature_k
        else:
            complete = np.r_[
                native,
                film.binary_mass_residual_kg_m2_s / mass_scale,
                film.heat_residual_w_m2 / heat_scale,
            ]
            surface_temperature = film.surface_temperature_k
        audit["last_complete_scaled_residuals"] = complete.tolist()
        audit["last_complete_maximum_scaled_residual"] = float(np.max(np.abs(complete)))
        audit["last_surface_unknowns"] = [surface_temperature, film.surface_y_hexane]

    def budgeted_assemble(*args, **kwargs):
        if audit["assemblies_started"] >= limit:
            audit["budget_exhausted"] = True
            raise AssemblyBudgetExceeded(
                f"native assembly budget exhausted: {limit} calls; "
                f"stage={audit['stage']}; "
                f"last scaled residual={audit['last_complete_maximum_scaled_residual']}"
            )
        if (
            remaining_wall_seconds is not None
            and time.perf_counter() - started >= remaining_wall_seconds
        ):
            audit["wall_budget_exhausted"] = True
            raise AssemblyWallBudgetExceeded(
                "bounded active-Faner wall budget exhausted between assemblies"
            )
        audit["assemblies_started"] += 1
        if "seed_outer_cell_temperature" not in captured:
            candidate = args[1] if len(args) > 1 else kwargs["candidate"]
            captured["seed_outer_cell_temperature"] = candidate.dry_temperatures_k[-1]
        try:
            assembly = original_assemble(*args, **kwargs)
        except Exception as exc:
            audit["invalid_trial_assemblies"] += 1
            audit["last_assembly_error"] = str(exc)
            audit["last_assembly_error_type"] = type(exc).__name__
            raise
        audit["assemblies_returned"] += 1
        captured["last_assembly"] = assembly
        record_complete(assembly)
        return assembly

    def capturing_scales(*args, **kwargs):
        result = original_scales(*args, **kwargs)
        captured["scales"] = result
        audit["stage"] = "joint_radial_film_nonlinear_solve"
        if captured["last_assembly"] is not None:
            record_complete(captured["last_assembly"])
        return result

    core.cut.assemble_backward_euler = budgeted_assemble
    core.ci._residual_scales = capturing_scales
    try:
        yield audit
    finally:
        core.cut.assemble_backward_euler = original_assemble
        core.ci._residual_scales = original_scales
        audit["wall_seconds"] = time.perf_counter() - started
        audit["native_callables_restored"] = (
            core.cut.assemble_backward_euler is original_assemble
            and core.ci._residual_scales is original_scales
        )
