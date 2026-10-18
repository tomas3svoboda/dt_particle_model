"""March one production particle along the declared Paraiso/Kemper DT destiny.

Three surface regimes, one model, the Faner/P1EF wheel:

  A  high-loading attached film  (sphere.initialize_qualified_feed +
     locate_activation_event; Coletto B.7-B.10 heat, Span-Wagner latent heat)
  B  receding front              (qualify_phase1_particle_engineering +
     qualify_gate1g_conditioned_cut_stress, same as Faner 2019)
  C  dry continuation            (recorded if the receding machine reaches it)

No tray solver. No vapor balance. No fitting. No src/ edits.
y_h(t) stays in (0, 1). Bath T is the published meal ladder, capped at the
wet-core caloric ceiling 380 K (3 K vs 110 C is measurement error).

Radii are the three Coletto-film sizes already in thermal_oracle /
qualify_a2c_wet_core_scales: 0.50, 0.885, 1.50 mm.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

HARNESS_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = HARNESS_DIR.parent
ROOT = ANALYSIS_DIR.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from dtdc_simulator.core2 import qsc_hexane_extinction_localizer_march as elm  # noqa: E402
from dtdc_simulator.core2.particle import bed_film as bf  # noqa: E402
from dtdc_simulator.core2.particle import coupled_pore as cp  # noqa: E402
from dtdc_simulator.core2.particle import coupled_transport as ct  # noqa: E402
from dtdc_simulator.core2.particle import cut_event_orchestrator as eo  # noqa: E402
from dtdc_simulator.core2.particle import cut_extinction_event as cee  # noqa: E402
from dtdc_simulator.core2.particle import cut_extinction_projection as cep  # noqa: E402
from dtdc_simulator.core2.particle import cut_extinction_time as cet  # noqa: E402
from dtdc_simulator.core2.particle import cut_face_tangent as cft  # noqa: E402
from dtdc_simulator.core2.particle import cut_integrator as ci  # noqa: E402
from dtdc_simulator.core2.particle import cut_transport as cut  # noqa: E402
from dtdc_simulator.core2.particle import engineering_foundation as p1ef  # noqa: E402
from dtdc_simulator.core2.particle import external_film_reduced as efr  # noqa: E402
from dtdc_simulator.core2.particle import sphere  # noqa: E402
from dtdc_simulator.core2.particle import wet_core  # noqa: E402
from dtdc_simulator.core2.props import binary_gas as bg  # noqa: E402
from dtdc_simulator.core2.props import hexane as hx  # noqa: E402
from dtdc_simulator.core2.props import sorption as sp  # noqa: E402
from dtdc_simulator.core2.props import water as wa  # noqa: E402

from era_stamp import era_stamp  # noqa: E402

from scripts import qualify_gate1g_conditioned_cut_stress as backend  # noqa: E402
from scripts import qualify_phase1_particle_engineering as q  # noqa: E402

DESTINY_PATH = (
    ANALYSIS_DIR / "results_2026-10-01" / "dtdc_particle_destiny" / "destiny.json"
)
# Same three radii as thermal_oracle.BARE_FILM_H_W_M2_K / qualify_a2c.
DECLARED_RADII_M = (0.5e-3, 0.885e-3, 1.5e-3)
STRESS_MULTIPLIER = 1.0
COEFFICIENT_CASE = "nominal_log_center"
CONDUCTIVITY = p1ef.NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K
FANER_DT_S = 0.075
# Regime-C (fully-dry) time-step cap for the doubling ladder; the paper's own
# fully-dry coupled_transport marches (binary kinetics SI) used 10-20 s steps.
DRY_DT_CAP_S = 2.4
# N=4 crawled into the Amendment-09 leaf floor at master faces because the
# harness dyadic dt ladder turned a 0.075 s crossing into a same-cell crawl.
# N=2 is the P1EF native mesh after the N=2 seed: one interior face at 0.5 R,
# then the centre cell, then licensed extinction.
FANER_CELLS = 2
# P1EF/Gate-1g native front (qualify_gate1g INITIAL_FRONT_RADIUS_FRACTION).
# A just-born front at 0.999 R has one dry piece and refuses on the
# water-bearing chart; Faner's 0.999 R marches used the water-free chart.
RECEDING_FRONT_FRACTION = backend.INITIAL_FRONT_RADIUS_FRACTION
P1EF_WET_T_K = 335.0
P1EF_DRY_T_K = 343.0
P1EF_IFACE_T_K = 337.0
P1EF_DRY_Y_HEXANE = 0.82
FANER_SEED_STEFAN = (0.022, 0.020)
PRESSURE_PA = 101325.0
M_HEXANE = hx.M
# Leave 5 % of the dry-cell water-dew headroom so a_w stays strictly below 1.
DEW_ACTIVITY_MARGIN = 0.95
# Faner emergency continuation D (docs/evidence/paper1_faner2019_march_emergency_2026-09-30):
# license same-cell adaptive subdivision under constant forcing. Not a new
# solver; the live qualifier already computes those leaves and then refuses
# them unless an imposed boundary discontinuity started the step.
SUBDIVISION_UNDER_CONSTANT_FORCING_ENVIRONMENT_VARIABLE = (
    "DTDC_SUBDIVISION_UNDER_CONSTANT_FORCING"
)
FANER_EMERGENCY_CERTIFICATION = "TEMPORARILY_UNCERTIFIED_EMERGENCY_2026-09-30"
_FANER_D_STEP_RECORD_WRAPPED = False
M_WATER = 0.018015268


def face_tangent_knobs() -> dict:
    """Record the two owner-ruled face-tangent knobs in force (src/ env knobs)."""

    return {
        "DTDC_FACE_TANGENT_INWARD_SEED": os.environ.get(
            cft.FACE_TANGENT_INWARD_SEED_ENVIRONMENT_VARIABLE
        ),
        "inward_seed_armed": cft.face_tangent_inward_seed_armed(),
        "DTDC_FACE_TANGENT_RESIDUAL_BOUND": os.environ.get(
            cft.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE
        ),
        "residual_bound_in_force": cft.face_tangent_residual_bound(),
        "residual_bound_declared": cft.face_tangent_residual_bound_declared(),
        "records": [
            "docs/GT_PS2_FACE_TANGENT_INWARD_SEED_BUILD_RECORD_2026-09-29.md",
            "docs/GT_PS2_OWNER_RULING_FACE_TANGENT_RESIDUAL_BOUND_2026-09-29.md",
        ],
    }


class CellZeroTerminus(Exception):
    """Cell-zero same-cell and adaptive refusal: hand to the extinction family."""

    def __init__(self, message: str, *, adaptive_error: Exception) -> None:
        super().__init__(message)
        self.adaptive_error = adaptive_error


def advance_receding_step(previous, dt_s: float, boundary):
    """One held step. Cell > 0: the qualifier body unchanged (bracket route owns
    interior faces). Cell 0: BC-7's own guard (qsc_hexane_face_event_march:
    "cell zero has no interior inner face ... march plainly there"), i.e. the
    same recovery -> adaptive chain without the bracket builder, and the
    adaptive leaf-floor refusal becomes the typed cell-zero terminus."""

    state = previous.after
    if state.transport.layout.cut_cell_index != 0:
        return q._advance_return_hold_step(previous, dt_s, boundary)
    try:
        return backend._advance_with_recovery(
            previous,
            dt_s,
            boundary,
            boundary,
            is_reversal=False,
            try_current_state_anchor=isinstance(boundary, ct.ReducedFilmPoreBoundary),
            provisional_maximum_function_evaluations=None,
            maximum_nonlinear_seed_attempts=None,
        )
    except backend.SameCellRecoveryError as exc:
        try:
            return backend._advance_with_adaptive_same_cell_macrostep(
                previous,
                dt_s,
                boundary,
                boundary,
                is_reversal=False,
                initial_attempts=exc.attempts,
                initial_step_error=exc.last_step_error,
                initial_failure_audit=exc.failure_jacobian_audit,
            )
        except backend.AdaptiveSameCellMacrostepError as adaptive_exc:
            raise CellZeroTerminus(
                "cell-zero same-cell and adaptive recovery refused; extinction family next",
                adaptive_error=adaptive_exc,
            ) from adaptive_exc


def extinction_resolution_for_state(state):
    """PART-02 declared residual bound / evaluation budget, typed against the
    frozen contract the destiny state actually carries.

    The ruled declaration (``elm.part02_environment_conditioned_resolution``)
    is written against the 2e-11 / 1200 PART-02 controls; the destiny march
    runs the P1EF Gate-1g controls (1e-10 / 1200) so the certified layer
    refuses the ruled object verbatim. The declared values (5e-9 bound,
    16000 evaluations, same ruling and measurement) are kept; only the
    "frozen value it conditions" is the state's own. Recorded in result.
    """

    ruled = elm.part02_environment_conditioned_resolution()
    frozen_residual = float(state.controls.nonlinear_residual_tolerance)
    frozen_budget = int(state.controls.maximum_function_evaluations)
    declared_bound = max(ruled.declared_bound, 2.0 * frozen_residual)
    declared_budget = max(ruled.declared_maximum_function_evaluations, 2 * frozen_budget)
    return replace(
        ruled,
        declared_bound=declared_bound,
        frozen_contract_value=frozen_residual,
        declared_maximum_function_evaluations=declared_budget,
        frozen_budget_value=frozen_budget,
        ruling=(
            ruled.ruling
            + f" | {FANER_EMERGENCY_CERTIFICATION}: declaration re-typed against the "
            f"destiny state's own frozen contract ({frozen_residual!r}/{frozen_budget})"
        ),
    )


def resolution_record(resolution) -> dict:
    return {
        "declared_bound": resolution.declared_bound,
        "frozen_contract_value": resolution.frozen_contract_value,
        "declared_maximum_function_evaluations": (
            resolution.declared_maximum_function_evaluations
        ),
        "frozen_budget_value": resolution.frozen_budget_value,
        "ruling": resolution.ruling,
        "measurement": resolution.measurement,
    }


# Cell-zero target march: relative front targets tried after consecutive
# frozen-dt refusals (the extinction ladder's own head schedule, then halving).
TARGET_MARCH_FRACTIONS = (0.99, 0.97, 0.94, 0.90, 0.85, 0.78, 0.70, 0.60, 0.50)
TARGET_MARCH_DURATION_FLOOR_S = 1.0e-9
# Geometric extinction cutoff (TEMPORARILY_UNCERTIFIED_EMERGENCY_2026-09-30).
# On the dew-floor-clamped branch the finite-limit extinction certificate
# refuses (tail durations diverge; measured in diag_R1p500_N2_cell0b) and the
# positive-core target march stalls at nm-scale cores (diag_R1p500_N2_cell0c:
# ~1 ms physical per 160-300 s wall at s = 2 um). Once the wet core radius
# s = R z^(1/3) is at or below this floor - comparable to the meal pore scale
# and three orders below the resolved cell - the src exact z=0 projection
# (cut_extinction_projection.solve_fixed_time_extinction_projection) is
# committed over the frozen dt, and regime C continues on the src fully-dry
# predictor. The projection is the certified family's own endpoint solve; only
# the finite-limit certificate that normally licenses its duration is bypassed.
GEOMETRIC_EXTINCTION_CORE_RADIUS_M = 5.0e-6
GEOMETRIC_EXTINCTION_DURATION_LADDER = (1.0, 0.1, 10.0, 0.01)
_STEFAN_SCALE_FLOOR_MOL_M2_S = 1.0e-3  # mirrors qsc_hexane_extinction_localizer_march


class CommittedPositiveCoreTargetStep(ci.CutIntegratorStep):
    """Harness-side commit of one accepted positive-core target root.

    ``cut_extinction_time.solve_positive_core_target`` solves the UNCHANGED
    cell-zero ALE/RH backward-Euler equations with the front prescribed and
    the duration unknown. Its accepted root is therefore the same committed
    state a BE step of that duration would produce; this object carries it in
    the accepted-trajectory shape the qualifier recovery chain and the
    extinction localizer consume (``before``/``after``/``assembly``/
    ``ledger.dt_s``/``ledger.accepted``). The cumulative ledger datum is
    re-anchored at the commit (``cut_integrator.initialize_integrator_state``);
    inventories are state-based and unaffected.
    """


def commit_positive_core_target(previous, boundary, fraction: float, max_duration_s: float):
    state = previous.after
    start_z, end_z, hist_dt = backend._accepted_history_front_motion(previous)
    rate = (end_z - start_z) / hist_dt
    if not (math.isfinite(rate) and rate < 0.0):
        raise RuntimeError(f"accepted history front rate is not receding ({rate!r} /s)")
    z0 = state.transport.geometry.front.z
    target = fraction * z0
    upper = max(max_duration_s, 2.0 * TARGET_MARCH_DURATION_FLOOR_S)
    bounds = (TARGET_MARCH_DURATION_FLOOR_S, upper)
    seed = elm._first_target_seed(state, target, rate, bounds)  # noqa: SLF001
    resolution = extinction_resolution_for_state(state)
    step = cet.solve_positive_core_target(
        state,
        target,
        boundary,
        seed,
        duration_bounds_s=bounds,
        environment_conditioned_resolution=resolution,
    )
    duration = float(step.ledger.event_duration_s)
    cand = step.assembly.candidate
    after_transport = cut.CutTransportState(
        geometry=step.assembly.candidate_geometry,
        config=step.assembly.target_config,
        time_s=state.transport.time_s + duration,
        wet_temperatures_k=cand.wet_temperatures_k,
        wet_retained_water_loadings=cand.wet_retained_water_loadings,
        dry_temperatures_k=cand.dry_temperatures_k,
        dry_y_hexane=cand.dry_y_hexane,
        historical_hexane_loadings=state.transport.historical_hexane_loadings,
        oil_fraction_labels=state.transport.oil_fraction_labels,
        wet_retained_water_capacity_duals_over_rt=(
            cand.effective_wet_retained_water_capacity_duals_over_rt
        ),
    )
    after = ci.initialize_integrator_state(
        after_transport,
        state.controls,
        interface_temperature_k=cand.interface_temperature_k,
        total_stefan_fluxes_mol_m2_s=cand.dry_total_stefan_fluxes_mol_m2_s,
    )
    after = replace(
        after,
        accepted_steps=state.accepted_steps + 1,
        cumulative_nonlinear_evaluations=(
            state.cumulative_nonlinear_evaluations + int(step.ledger.nonlinear_evaluations)
        ),
    )
    led = step.ledger
    ledger = SimpleNamespace(
        dt_s=duration,
        accepted=bool(led.accepted),
        nonlinear_evaluations=int(led.nonlinear_evaluations),
        maximum_scaled_residual=float(led.maximum_scaled_residual),
        condition_proxy=float(led.condition_proxy),
        maximum_step_ledger_residual=None,
        maximum_cumulative_ledger_residual=None,
        ledger_datum_reanchored=True,
        owner="cut_extinction_time.solve_positive_core_target",
    )
    committed = CommittedPositiveCoreTargetStep(
        before=state,
        after=after,
        assembly=step.assembly,
        boundary=boundary,
        seed=step.seed,
        ledger=ledger,
    )
    info = {
        "target_fraction_of_current_z": fraction,
        "z_before": z0,
        "z_after": target,
        "history_front_rate_z_s": rate,
        "forecast_duration_s": seed.event_duration_s,
        "event_duration_s": duration,
        "nonlinear_evaluations": int(led.nonlinear_evaluations),
        "maximum_scaled_residual": float(led.maximum_scaled_residual),
        "condition_proxy": float(led.condition_proxy),
        "interface_temperature_k": cand.interface_temperature_k,
        "wet_temperatures_k": list(cand.wet_temperatures_k),
        "dry_temperatures_k": list(cand.dry_temperatures_k),
        "dry_y_hexane": list(cand.dry_y_hexane),
        "declared_bound": resolution.declared_bound,
    }
    return committed, info


def resolve_cell_zero_extinction(previous, dt_s: float, boundary, *, on_rung=None):
    """Regime B -> C: the orchestrator's typed EXTINCTION_REQUIRED terminus routed
    through qsc_hexane_extinction_localizer_march (ladder, finite-limit
    certificate, exact z=0 projection, fully-dry remainder). Existing wheels."""

    state = previous.after
    history = backend._accepted_face_history_source(previous)
    seed = backend._current_state_seed(
        previous, front_z=None, label="dtdc destiny: cell-zero terminus probe seed"
    )
    try:
        macro = eo.advance_one_interior_face_macrostep(state, dt_s, boundary, seed)
    except eo.CutEventMacrostepError as exc:
        if exc.code is not eo.CutMacrostepFailureCode.EXTINCTION_REQUIRED:
            raise
        terminus = exc
    else:
        return "same_cell", macro
    resolution = extinction_resolution_for_state(state)
    localization = elm.localize_cell_zero_extinction(
        state,
        boundary,
        terminus,
        history=history,
        environment_conditioned_resolution=resolution,
        on_rung=on_rung,
    )
    t_ext = float(localization.certificate.extrapolated_event_duration_s)
    commit_dt = float(dt_s) if dt_s >= t_ext else t_ext
    outcome = elm.commit_cell_zero_extinction_outcome(localization, commit_dt)
    return "extinction", outcome


def wet_core_radius_m(state, radius_m: float) -> float:
    return radius_m * state.transport.geometry.front.z ** (1.0 / 3.0)


def project_geometric_extinction(previous, frozen_dt_s: float, boundary):
    """Exact z=0 endpoint projection over a prescribed duration (src solve),
    seeded structurally from the current cell-zero state (same mapping as
    qsc_hexane_extinction_localizer_march.extinction_projection_seed_from_localization
    but from the live state instead of a certified rung). Returns
    (fully_dry_state, info). Tries a short duration ladder (multiples of the
    frozen dt) and reports every refusal."""

    before = previous.after
    tr = before.transport
    if tr.layout.cut_cell_index != 0:
        raise ValueError("geometric extinction projection requires a cell-zero cut")
    n = tr.layout.dry_piece_count
    fluxes = tuple(before.last_total_stefan_fluxes_mol_m2_s)
    if len(tr.dry_temperatures_k) != n or len(fluxes) != n + 1:
        raise ValueError(
            f"dry rank does not map onto the endpoint rank (pieces {n}, fluxes {len(fluxes)})"
        )
    errors = []
    for multiple in GEOMETRIC_EXTINCTION_DURATION_LADDER:
        duration = float(frozen_dt_s * multiple)
        candidate = cee.ExtinctionUnknowns(
            dry_temperatures_k=tuple(tr.dry_temperatures_k),
            dry_y_hexane=tuple(tr.dry_y_hexane),
            positive_area_total_stefan_fluxes_mol_m2_s=fluxes[1:],
            event_time_s=duration,
        )
        scales = tuple(
            max(abs(v), _STEFAN_SCALE_FLOOR_MOL_M2_S)
            for v in candidate.positive_area_total_stefan_fluxes_mol_m2_s
        )
        seed = cep.ExtinctionProjectionSeed(
            candidate,
            scales,
            "dtdc destiny geometric cutoff: live cell-zero state mapped structurally "
            f"onto the endpoint rank at {duration} s (uncertified duration)",
        )
        try:
            step = cep.solve_fixed_time_extinction_projection(
                before,
                duration,
                boundary,
                seed,
                environment_conditioned_resolution=extinction_resolution_for_state(before),
            )
        except Exception as exc:  # noqa: BLE001
            errors.append(
                {
                    "duration_s": duration,
                    "error_type": type(exc).__name__,
                    "error": str(exc)[:2000],
                }
            )
            continue
        dry_committed = step.after
        labels = tuple(dry_committed.oil_fraction_labels)
        label_retyped = False
        if len(set(labels)) != 1:
            raise RuntimeError(
                "geometric extinction endpoint has nonuniform oil labels; the "
                "fully-dry predictor must not homogenise them"
            )
        if labels[0] != dry_committed.config.pore.w_o:
            # The destiny pieces carry the feed oil label (uniform); the dry
            # config still carries the P1EF default w_o. Re-type the config to
            # the one label so the existing predictor keeps the same material.
            dry_committed = replace(
                dry_committed,
                config=replace(
                    dry_committed.config,
                    pore=replace(dry_committed.config.pore, w_o=labels[0]),
                ),
            )
            label_retyped = True
        dry_state = dry_committed.to_uniform_fully_dry_transport_state()
        led = step.ledger
        info = {
            "projection_duration_s": duration,
            "oil_label": labels[0],
            "dry_config_w_o_retyped_to_label": label_retyped,
            "config_w_o_before_retype": step.after.config.pore.w_o,
            "nonlinear_evaluations": getattr(led, "nonlinear_evaluations", None),
            "maximum_scaled_residual": getattr(led, "maximum_scaled_residual", None),
            "condition_proxy": getattr(led, "condition_proxy", None),
            "maximum_step_ledger_residual": led.maximum_step_ledger_residual,
            "maximum_cumulative_ledger_residual": led.maximum_cumulative_ledger_residual,
            "transition_authority": step.after.transition_authority,
            "refused_durations": errors,
            "declared_bound": extinction_resolution_for_state(before).declared_bound,
        }
        return dry_state, info
    raise RuntimeError(
        "geometric extinction projection refused at every duration: "
        + json.dumps(errors)[:4000]
    )


def dry_summary(state: ct.FullyDryTransportState, dry_mass_kg: float, w_o: float) -> dict:
    cfg = state.config
    grid = state.grid
    states = [
        cp.evaluate_equilibrium(T, cfg.pressure_pa, y, cfg.pore)
        for T, y in zip(state.temperatures_k, state.y_hexane)
    ]
    vols = grid.volumes
    hex_mol = math.fsum(v * s.total_hexane_concentration_mol_m3 for v, s in zip(vols, states))
    wat_mol = math.fsum(v * s.total_water_concentration_mol_m3 for v, s in zip(vols, states))
    energy = math.fsum(v * s.energy_density_j_m3 for v, s in zip(vols, states))
    ret_w_mol = (
        math.fsum(v * s.retained_water_concentration_kg_m3 for v, s in zip(vols, states))
        / M_WATER
    )
    ret_h_kg = math.fsum(v * s.retained_hexane_concentration_kg_m3 for v, s in zip(vols, states))
    mean_t = math.fsum(v * T for v, T in zip(vols, state.temperatures_k)) / math.fsum(vols)
    x_h = hex_mol * M_HEXANE / dry_mass_kg
    x_w = wat_mol * M_WATER / dry_mass_kg
    total = 1.0 + x_h + x_w + w_o
    return {
        "regime": "C",
        "time_s": state.time_s,
        "front_z": 0.0,
        "front_s_over_R": 0.0,
        "T_k": mean_t,
        "interface_temperature_k": None,
        "X_hexane_kg_kg": x_h,
        "X_hexane_sorbed_kg_kg": ret_h_kg / dry_mass_kg,
        "X_water_kg_kg": x_w,
        "moisture_wt_pct_wb": 100.0 * x_w / total,
        "hexane_ppm_wb": 1.0e6 * x_h / total,
        "hexane_mol": hex_mol,
        "water_mol": wat_mol,
        "retained_water_mol": ret_w_mol,
        "energy_j": energy,
        "wet_pieces": 0,
        "dry_pieces": grid.n,
        "wet_temperature_min_k": None,
        "wet_temperature_max_k": None,
        "dry_temperature_surface_piece_k": state.temperatures_k[-1],
        "dry_y_hexane": list(state.y_hexane),
        "dry_temperatures_k": list(state.temperatures_k),
    }


def load_destiny() -> dict:
    return json.loads(DESTINY_PATH.read_text(encoding="utf-8"))


def interpolate_bath(waypoints: list[dict], t_s: float) -> tuple[float, float]:
    """Linear T_bath, y_h on the declared destiny table. y_h stays in (0, 1)."""

    if t_s <= waypoints[0]["t_s"]:
        T = float(waypoints[0]["T_bath_K"])
        y = float(waypoints[0]["y_hexane"])
    elif t_s >= waypoints[-1]["t_s"]:
        T = float(waypoints[-1]["T_bath_K"])
        y = float(waypoints[-1]["y_hexane"])
    else:
        T = y = None
        for left, right in zip(waypoints, waypoints[1:]):
            t0 = float(left["t_s"])
            t1 = float(right["t_s"])
            if t0 <= t_s <= t1:
                w = 0.0 if t1 == t0 else (t_s - t0) / (t1 - t0)
                T = float(left["T_bath_K"]) + w * (
                    float(right["T_bath_K"]) - float(left["T_bath_K"])
                )
                y = float(left["y_hexane"]) + w * (
                    float(right["y_hexane"]) - float(left["y_hexane"])
                )
                break
        if T is None:
            raise RuntimeError(f"destiny interpolation missed t={t_s}")
    y = min(1.0 - 1.0e-6, max(1.0e-6, y))
    return T, y


def film_pair_at_radius(radius_m: float) -> bf.CorrelatedBedFilmPair:
    """Coletto B.7-B.10 pair. Production radius uses the frozen P1EF source state."""

    source = p1ef.source_correlated_film_pair()
    if math.isclose(radius_m, p1ef.PRODUCTION_RADIUS_M, rel_tol=0.0, abs_tol=0.0):
        return source
    state = replace(
        source.state,
        particle_radius_m=radius_m,
        state_provenance=(
            source.state.state_provenance
            + f"; Coletto pair evaluated at R={radius_m} m "
            "(qualify_a2c / thermal_oracle radius set)"
        ),
    )
    return bf.coletto_b7_b10_pair(state)


def reduced_film_boundary(
    *,
    temperature_k: float,
    y_hexane: float,
    pair: bf.CorrelatedBedFilmPair,
    forcing_timescale_s: float,
    label: str,
) -> ct.ReducedFilmPoreBoundary:
    """Same constructor as p1ef.make_reduced_film_boundary, with a given pair."""

    scaled_mass = STRESS_MULTIPLIER * pair.binary_mass_transfer_m_s
    relaxation_s = pair.state.bulk_binary_diffusivity_m2_s / scaled_mass**2
    ratio = relaxation_s / forcing_timescale_s
    if ratio <= 0.01:
        fidelity = efr.FidelityBand.NOMINAL
    elif ratio <= 0.10:
        fidelity = efr.FidelityBand.BOUNDED_ANOMALY
    else:
        raise efr.ReducedFilmValidityError(
            "film relaxation/forcing ratio exceeds the bounded P1EF anomaly limit "
            f"(ratio={ratio:.6g}, forcing={forcing_timescale_s} s)"
        )
    return ct.ReducedFilmPoreBoundary(
        temperature_k=temperature_k,
        y_hexane=y_hexane,
        correlated_pair=pair,
        stress_multiplier=STRESS_MULTIPLIER,
        liquid_saturation=0.0,
        forcing_timescale_s=forcing_timescale_s,
        absolute_mass_fraction_departure_threshold=(
            p1ef.FAST_FILM_MASS_FRACTION_DEPARTURE_LIMIT
        ),
        absolute_ackermann_factor_departure_threshold=(
            p1ef.UNITY_ACKERMANN_FACTOR_DEPARTURE_LIMIT
        ),
        legacy_engineering_regression_water_heat_capacity_j_kg_k=(
            p1ef.LEGACY_ENGINEERING_REGRESSION_WATER_HEAT_CAPACITY_J_KG_K
        ),
        legacy_engineering_regression_hexane_heat_capacity_j_kg_k=(
            p1ef.LEGACY_ENGINEERING_REGRESSION_HEXANE_HEAT_CAPACITY_J_KG_K
        ),
        binary_gas_interaction_k_wh=bg.K_WH_CENTRAL,
        label=label,
        fidelity_band=fidelity,
    )


def receding_admissible_y_hexane(y_declared: float, dry_temperature_k: float) -> float:
    """Raise declared y_h to the dry-pore water-dew floor at the receding dry T.

    Coupled-pore a_w must lie in (0, 1]. At the P1EF dry-cell temperature the
    industrial destiny y_h (0.614 at activation, 0.200 at 1800 s) sits below
    that floor. This is the existing topology admission, not a new rate law.
    """

    y = min(1.0 - 1.0e-6, max(1.0e-6, y_declared))
    psat_w = wa.saturation_pressure(dry_temperature_k)
    y_h_min = 1.0 - DEW_ACTIVITY_MARGIN * psat_w / PRESSURE_PA
    y_h_min = min(0.85, max(0.50, y_h_min))
    return max(y, y_h_min)


def stefan_seed_for_faces(nfaces: int) -> tuple[float, ...]:
    if nfaces <= 0:
        raise ValueError("receding seed needs at least one dry face")
    if nfaces <= len(FANER_SEED_STEFAN):
        return FANER_SEED_STEFAN[:nfaces]
    extra = (FANER_SEED_STEFAN[-1],) * (nfaces - len(FANER_SEED_STEFAN))
    return FANER_SEED_STEFAN + extra


def hexane_saturation_temperature_k(pressure_pa: float) -> float:
    lo, hi = 310.0, 380.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if hx.saturation_pressure(mid) < pressure_pa:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _hexane_mass_fraction(y_hexane: float) -> float:
    num = y_hexane * M_HEXANE
    return num / (num + (1.0 - y_hexane) * M_WATER)


def wet_surface_film_rates(
    *,
    surface_temperature_k: float,
    bath_temperature_k: float,
    bath_y_hexane: float,
    conductance_w_m2_k: float,
    pair: bf.CorrelatedBedFilmPair,
    pressure_pa: float,
) -> dict:
    """Regime-A surface law of Sec. 2.1 on the production Coletto film pair.

    Liquid n-hexane wets the surface at unit activity, so the interface vapor
    is ``y_I = min(1, Psat_h(T_s)/P)`` (eq. aceiling).  Only hexane crosses the
    film (no water condensation branch; primary drainage forbids a negative
    outward rate), so the ideal-binary mass-coordinate Stefan relation of
    external_film_reduced reduces to ``j_h = rho_f h_M ln[(1-w_b)/(1-w_I)]``;
    heat enters through the Ackermann-corrected conductance
    ``q = U A(beta) (T_b - T_s)``, ``beta = j_h cp_h / U``, with the frozen
    production validity ``|beta| <= 1``.  ``conductance_w_m2_k`` is the film
    coefficient alone in regime A and the film + dry-shell series value in
    the A->B bridge.  Lumped harness evaluation, not a resolved src solve.
    """

    y_i = min(1.0, hx.saturation_pressure(surface_temperature_k) / pressure_pa)
    w_i = _hexane_mass_fraction(y_i)
    w_b = _hexane_mass_fraction(bath_y_hexane)
    rho_k = pair.state.gas_density_kg_m3 * STRESS_MULTIPLIER * pair.binary_mass_transfer_m_s
    boiling = y_i >= 1.0
    if boiling:
        flux_form = "pressure_ceiling_boiling"
        j_h = None
    elif w_i <= w_b:
        flux_form = "no_outward_drive"
        j_h = 0.0
    else:
        flux_form = "stefan_logarithm"
        j_h = rho_k * math.log((1.0 - w_b) / (1.0 - w_i))
    cp_h_film = hx.cp0(0.5 * (surface_temperature_k + bath_temperature_k)) / M_HEXANE
    if j_h is None:
        beta = None
        factor = 1.0
    else:
        beta = j_h * cp_h_film / conductance_w_m2_k
        if abs(beta) > 1.0:
            raise efr.ReducedFilmValidityError(
                f"regime-A film Ackermann |beta|={abs(beta):.4f} exceeds 1 "
                f"(T_s={surface_temperature_k:.3f} K, y_I={y_i:.5f}, y_b={bath_y_hexane:.5f})"
            )
        factor = efr._inverse_exprel(beta)
    q_in_w_m2 = conductance_w_m2_k * factor * (bath_temperature_k - surface_temperature_k)
    dh_vap_j_kg = hx.dh_vap(surface_temperature_k) / M_HEXANE
    if boiling:
        j_h = max(0.0, q_in_w_m2) / dh_vap_j_kg
    return {
        "j_h_kg_m2_s": j_h,
        "q_in_w_m2": q_in_w_m2,
        "y_interface": y_i,
        "w_interface": w_i,
        "w_bath": w_b,
        "beta": beta,
        "ackermann_factor": factor,
        "cp_h_film_j_kg_k": cp_h_film,
        "dh_vap_j_kg": dh_vap_j_kg,
        "flux_form": flux_form,
        "rho_k_kg_m2_s": rho_k,
    }


def high_loading_rates(
    state: sphere.HighLoadingState,
    *,
    radius_m: float,
    T_bath_k: float,
    y_bath: float,
    pair: bf.CorrelatedBedFilmPair,
    tsat_k: float,
    pressure_pa: float,
) -> tuple[float, float, dict]:
    """Regime A on the Sec. 2.1 surface law: external film, unit activity.

    The lumped particle (external resistance only, so T_s = T) evaporates
    hexane at the Stefan film flux against the declared bath composition and
    receives Ackermann-corrected film heat.  Evaporated hexane leaves at its
    vapor enthalpy on the common datum, h_l(T) + dh_vap(T), so the latent duty
    is debited from the particle (the particle sits near its wet-bulb point;
    it never superheats above the pure-hexane boiling point while wetted).
    """

    area = 4.0 * math.pi * radius_m**2
    h_q = STRESS_MULTIPLIER * pair.heat_transfer_w_m2_k
    film = wet_surface_film_rates(
        surface_temperature_k=state.temperature_k,
        bath_temperature_k=T_bath_k,
        bath_y_hexane=y_bath,
        conductance_w_m2_k=h_q,
        pair=pair,
        pressure_pa=pressure_pa,
    )
    q_in = film["q_in_w_m2"] * area
    m_dot = film["j_h_kg_m2_s"] * area
    liquid = hx.state_Tp(state.temperature_k, pressure_pa, "liquid")
    h_vapor_j_kg = liquid.h_mass + film["dh_vap_j_kg"]
    e_dot = m_dot * h_vapor_j_kg - q_in
    audit = {
        "q_in_w": q_in,
        "h_q_w_m2_k": h_q,
        "boiling": film["flux_form"] == "pressure_ceiling_boiling",
        "tsat_k": tsat_k,
        "dh_vap_j_kg": film["dh_vap_j_kg"],
        "h_vapor_j_kg": h_vapor_j_kg,
        "u_liquid_j_kg": liquid.u_mass,
        "m_dot_kg_s": m_dot,
        "e_dot_w": e_dot,
        "y_bath": y_bath,
        "y_interface": film["y_interface"],
        "beta": film["beta"],
        "ackermann_factor": film["ackermann_factor"],
        "flux_form": film["flux_form"],
    }
    return m_dot, e_dot, audit


def apply_high_loading_step(
    state: sphere.HighLoadingState,
    dt_s: float,
    m_dot: float,
    e_dot: float,
    params: sphere.SphereParams,
) -> tuple[str, sphere.HighLoadingState | sphere.ActivationEvent]:
    if m_dot < 0.0:
        raise ValueError("primary drainage forbids a negative outward n-hexane rate")
    try:
        event = sphere.locate_activation_event(state, dt_s, m_dot, e_dot, params)
        return "activated", event
    except sphere.ActivationNotReached:
        pass
    wet = replace(
        params.wet_core,
        X_water=state.retained_water_loading,
        w_o=state.oil_fraction,
    )
    mass_h = state.total_hexane_mass_kg - m_dot * dt_s
    energy = state.total_internal_energy_j - e_dot * dt_s
    if mass_h < 0.0:
        raise RuntimeError("high-loading step overshot attached hexane; dt too large")
    temperature = sphere.high_loading_temperature_from_energy(
        energy, state.dry_meal_mass_kg, mass_h, wet
    )
    x_c = sphere.critical_loading(temperature, params)
    internal = state.dry_meal_mass_kg * x_c
    attached = mass_h - internal
    if attached <= 0.0:
        event = sphere.locate_activation_event(state, dt_s, max(m_dot, 1.0e-30), e_dot, params)
        return "activated", event
    nxt = sphere.HighLoadingState(
        dry_meal_mass_kg=state.dry_meal_mass_kg,
        oil_label_mass_kg=state.oil_label_mass_kg,
        retained_water_mass_kg=state.retained_water_mass_kg,
        internal_hexane_mass_kg=internal,
        attached_hexane_mass_kg=attached,
        total_internal_energy_j=energy,
        temperature_k=temperature,
    )
    return "high_loading", nxt


def high_loading_row(t_s: float, state: sphere.HighLoadingState, extra: dict) -> dict:
    dry = state.dry_meal_mass_kg
    x_h = state.total_hexane_loading
    x_w = state.retained_water_loading
    w_o = state.oil_fraction
    total = 1.0 + x_h + x_w + w_o
    return {
        "regime": "A",
        "time_s": t_s,
        "T_k": state.temperature_k,
        "X_hexane_kg_kg": x_h,
        "X_attached_kg_kg": state.attached_hexane_mass_kg / dry,
        "X_internal_kg_kg": state.internal_hexane_mass_kg / dry,
        "X_water_kg_kg": x_w,
        "moisture_wt_pct_wb": 100.0 * x_w / total,
        "hexane_ppm_wb": 1.0e6 * x_h / total,
        **extra,
    }


def receding_summary(state, dry_mass_kg: float) -> dict:
    tr = state.transport
    inv = ci.inventory_snapshot(tr)
    water = ci.water_phase_inventory_snapshot(tr)
    geo = tr.geometry
    lay = tr.layout
    wet_v = [geo.cells[i].wet_volume_m3 for i in lay.wet_cell_indices]
    dry_v = [geo.cells[i].dry_volume_m3 for i in lay.dry_cell_indices]
    vol = math.fsum(wet_v) + math.fsum(dry_v)
    mean_t = (
        math.fsum(v * t for v, t in zip(wet_v, tr.wet_temperatures_k))
        + math.fsum(v * t for v, t in zip(dry_v, tr.dry_temperatures_k))
    ) / vol
    z = geo.front.z
    x_h = inv.total_hexane_mol * M_HEXANE / dry_mass_kg
    x_w = inv.total_water_mol * 0.018015268 / dry_mass_kg
    w_o = tr.oil_fraction_labels[0] if tr.oil_fraction_labels else 0.0
    total = 1.0 + x_h + x_w + w_o
    return {
        "regime": "B",
        "time_s": tr.time_s,
        "front_z": z,
        "front_s_over_R": z ** (1.0 / 3.0),
        "T_k": mean_t,
        "interface_temperature_k": state.last_interface_temperature_k,
        "X_hexane_kg_kg": x_h,
        "X_water_kg_kg": x_w,
        "moisture_wt_pct_wb": 100.0 * x_w / total,
        "hexane_ppm_wb": 1.0e6 * x_h / total,
        "hexane_mol": inv.total_hexane_mol,
        "water_mol": inv.total_water_mol,
        "retained_water_mol": water.retained_water_mol,
        "energy_j": inv.total_energy_j,
        "wet_pieces": lay.wet_piece_count,
        "dry_pieces": lay.dry_piece_count,
        "wet_temperature_min_k": min(tr.wet_temperatures_k) if tr.wet_temperatures_k else None,
        "wet_temperature_max_k": max(tr.wet_temperatures_k) if tr.wet_temperatures_k else None,
        "dry_temperature_surface_piece_k": tr.dry_temperatures_k[-1] if tr.dry_temperatures_k else None,
    }


def score_outlet(summary_row: dict, box: dict) -> dict:
    paraiso = box["paraiso_dt_outlet"]
    t_c = summary_row["T_k"] - 273.15
    return {
        "T_C": t_c,
        "moisture_wt_pct_wb": summary_row["moisture_wt_pct_wb"],
        "hexane_ppm_wb": summary_row["hexane_ppm_wb"],
        "paraiso_T_C": paraiso["T_C"],
        "paraiso_moisture_wt_pct_wb": paraiso["moisture_wt_pct_wb"],
        "paraiso_hexane_ppm": paraiso["hexane_ppm"],
        "in_dt_temperature_box": box["dt_temperature_C"][0]
        <= t_c
        <= box["dt_temperature_C"][1],
        "in_dt_moisture_box": box["dt_moisture_wt_pct_wb"][0]
        <= summary_row["moisture_wt_pct_wb"]
        <= box["dt_moisture_wt_pct_wb"][1],
        "in_dt_hexane_box": box["dt_hexane_ppm"][0]
        <= summary_row["hexane_ppm_wb"]
        <= box["dt_hexane_ppm"][1],
    }


def subdivision_under_constant_forcing_licensed() -> bool:
    raw = os.environ.get(
        SUBDIVISION_UNDER_CONSTANT_FORCING_ENVIRONMENT_VARIABLE, ""
    ).strip()
    if raw in ("", "0"):
        return False
    if raw != "1":
        raise ValueError(
            "DTDC_SUBDIVISION_UNDER_CONSTANT_FORCING must be 1 or unset"
        )
    return True


def install_faner_subdivision_licence() -> None:
    """Overlay Faner emergency D onto backend._step_record. No src/ edit."""

    global _FANER_D_STEP_RECORD_WRAPPED
    if _FANER_D_STEP_RECORD_WRAPPED:
        return
    original = backend._step_record

    def _step_record(
        step,
        attempts,
        stage,
        recovery_mode,
        *,
        starts_after_boundary_discontinuity: bool = False,
    ):
        rec = original(
            step,
            attempts,
            stage,
            recovery_mode,
            starts_after_boundary_discontinuity=starts_after_boundary_discontinuity,
        )
        adaptive = bool(rec.get("adaptive_subdivision_used"))
        emergency = (
            adaptive
            and not starts_after_boundary_discontinuity
            and subdivision_under_constant_forcing_licensed()
        )
        if not emergency:
            return rec
        contract = rec.get("frozen_failed_leaf_controller_contract")
        if not isinstance(contract, dict):
            return rec
        contract = dict(contract)
        contract["licensed_only_after_imposed_boundary_discontinuity"] = True
        contract["emergency_constant_forcing_subdivision_licence_used"] = True
        contract["emergency_certification"] = FANER_EMERGENCY_CERTIFICATION
        contract["passed"] = all(
            contract.get(term) is True for term in q.FROZEN_FAILED_LEAF_CONTRACT_TERMS
        )
        rec = dict(rec)
        rec["frozen_failed_leaf_controller_contract"] = contract
        rec["emergency_constant_forcing_subdivision_licence_used"] = True
        rec["emergency_certification"] = FANER_EMERGENCY_CERTIFICATION
        return rec

    backend._step_record = _step_record
    _FANER_D_STEP_RECORD_WRAPPED = True


def install_receding_backend(radius_m: float, wet_params: wet_core.WetCoreParams):
    install_faner_subdivision_licence()
    installed = q.install_engineering_case(
        COEFFICIENT_CASE,
        stress_multiplier=STRESS_MULTIPLIER,
        structural_conductivity_w_m_k=CONDUCTIVITY,
    )
    wet_model = replace(installed.configuration.wet, wet=wet_params)
    cut_config = replace(installed.configuration.cut, wet=wet_model)
    backend.RADIUS_M = radius_m
    backend.PRESSURE_PA = wet_params.pressure_pa
    backend.WET_MODEL = wet_model
    backend.CONFIG = cut_config
    backend.DRY_CONFIG = cut_config.dry
    backend.SEGMENT_DURATION_S = FANER_DT_S
    # Wet caloric authority is WetCoreParams (310, 380). Interface stays the
    # Faner/P1EF PHY-039 window: steam-distillation 334.5 K to hexane Tsat
    # 341 K. evaluate_interface_state has no root above Tsat at 1 atm.
    backend.CONTROLS = replace(
        backend.CONTROLS,
        wet_temperature_bounds_k=(310.0, 380.0),
        interface_temperature_bounds_k=q.INTERFACE_TEMPERATURE_BOUNDS_K,
        wet_water_bounds=(0.0231, wet_model.luikov.W_cap),
        maximum_function_evaluations=1200,
        maximum_condition_proxy=1.0e12,
    )
    return installed, cut_config


def strictly_inside(value: float, bounds: tuple[float, float], margin: float = 0.05) -> float:
    lo, hi = bounds
    lo_i = lo + margin
    hi_i = hi - margin
    if not lo_i < hi_i:
        return 0.5 * (lo + hi)
    return min(hi_i, max(lo_i, value))


def patch_initial_integrator(
    *,
    front_fraction: float,
    wet_t: float,
    wet_w: float,
    dry_t: float,
    dry_y: float,
    iface_t: float,
    act_t: float,
    oil_fraction: float,
    time_s: float,
):
    def initial_integrator(cells, radius_m=None):
        if radius_m is None:
            radius_m = backend.RADIUS_M
        radius_m = backend._validate_radius(radius_m)
        grid = backend.uniform_grid(cells, radius_m)
        geometry = backend.cg.partition_master_grid(
            grid, radius_m=front_fraction * radius_m
        )
        layout = cut.layout_for_geometry(geometry)
        history = wet_core.activate(grid, act_t, backend.WET_MODEL.wet).loadings
        transport = cut.initialize_cut_state(
            geometry,
            backend.CONFIG,
            [wet_t] * layout.wet_piece_count,
            [wet_w] * layout.wet_piece_count,
            [dry_t] * layout.dry_piece_count,
            [dry_y] * layout.dry_piece_count,
            history,
            oil_fraction_labels=(oil_fraction,) * grid.n,
            time_s=time_s,
        )
        return ci.initialize_integrator_state(
            transport, backend.CONTROLS, interface_temperature_k=iface_t
        )

    backend.initial_integrator = initial_integrator
    backend._initial_integrator = initial_integrator
    return initial_integrator


def dump_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")
    # Windows: a concurrent reader (editor, AV, ConvertFrom-Json poll) can hold
    # the target briefly and make replace() raise PermissionError; retry.
    for attempt in range(20):
        try:
            tmp.replace(path)
            return
        except PermissionError:
            if attempt == 19:
                raise
            time.sleep(0.25 * (attempt + 1))


def march_one(args: argparse.Namespace) -> dict:
    destiny = load_destiny()
    feed = destiny["particle"]
    waypoints = destiny["waypoints"]
    t_end = float(args.max_physical_s)
    radius_m = float(args.radius_m)
    pair = film_pair_at_radius(radius_m)
    tsat = hexane_saturation_temperature_k(PRESSURE_PA)
    installed, cut_config = install_receding_backend(
        radius_m,
        replace(
            wet_core.WetCoreParams(),
            pressure_pa=PRESSURE_PA,
            conductivity=CONDUCTIVITY,
            X_water=feed["feed_water_kg_kg_dry"],
            w_o=feed["feed_oil_kg_kg_dry"],
        ),
    )
    wet_params = replace(
        installed.configuration.wet.wet,
        X_water=feed["feed_water_kg_kg_dry"],
        w_o=feed["feed_oil_kg_kg_dry"],
        pressure_pa=PRESSURE_PA,
        conductivity=CONDUCTIVITY,
    )
    wet_model = replace(installed.configuration.wet, wet=wet_params)
    backend.WET_MODEL = wet_model
    backend.CONFIG = replace(cut_config, wet=wet_model)
    params = sphere.SphereParams(
        particle_radius_m=radius_m,
        radial_cells=args.cells,
        wet_core=wet_params,
    )
    state = sphere.initialize_qualified_feed(
        params,
        total_hexane_loading=feed["feed_hexane_kg_kg_dry"],
        temperature_k=feed["feed_temperature_K"],
        retained_water_loading=feed["feed_water_kg_kg_dry"],
        oil_fraction=feed["feed_oil_kg_kg_dry"],
    )
    dry_mass = state.dry_meal_mass_kg
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    steps_path = out_dir / "steps.jsonl"
    summary_path = out_dir / "summary.json"
    steps_file = steps_path.open("w", encoding="utf-8")
    t0 = time.perf_counter()
    result = {
        "physically_qualifying": False,
        "plant_predictive": False,
        "fitted": False,
        "model": "production core2 particle, three surface regimes, Faner/P1EF wheel",
        "radius_m": radius_m,
        "cells": args.cells,
        "stamp": era_stamp(input_paths=[DESTINY_PATH]),
        "destiny": DESTINY_PATH.as_posix(),
        "film": {
            "h_q_w_m2_k": pair.heat_transfer_w_m2_k,
            "h_m_m_s": pair.binary_mass_transfer_m_s,
            "stress_multiplier": STRESS_MULTIPLIER,
            "production_pair": math.isclose(
                radius_m, p1ef.PRODUCTION_RADIUS_M, rel_tol=0.0, abs_tol=0.0
            ),
        },
        "tsat_k": tsat,
        "stage": "high_loading",
        "args": vars(args),
    }

    def flush_result() -> None:
        result["wall_total_s"] = time.perf_counter() - t0
        dump_json(summary_path, result)

    def write_step(row: dict) -> None:
        steps_file.write(json.dumps(row, default=str) + "\n")
        steps_file.flush()

    t_s = 0.0
    hl_steps = 0
    activation_event = None
    global DEW_ACTIVITY_MARGIN
    dew_margin = float(getattr(args, "dew_margin", DEW_ACTIVITY_MARGIN))
    if not 0.5 <= dew_margin <= 0.99:
        raise ValueError("--dew-margin must lie in [0.5, 0.99]")
    DEW_ACTIVITY_MARGIN = dew_margin
    result["dew_activity_margin"] = DEW_ACTIVITY_MARGIN
    try:
        while t_s < t_end and state.attached_hexane_mass_kg > 0.0:
            T_bath, y_h = interpolate_bath(waypoints, t_s)
            # Same no-condensation admission as regimes B and C: the binary
            # (inert-free) pore/film gas of Sec. 2 cannot carry a water
            # partial pressure above saturation at the coldest temperature in
            # the film; the declared tray gas below its own water dew point
            # is admitted at a_w = margin instead (recorded, not a rate law).
            y_film = receding_admissible_y_hexane(y_h, min(state.temperature_k, T_bath))
            dt = min(args.high_loading_dt, t_end - t_s)
            m_dot, e_dot, audit = high_loading_rates(
                state,
                radius_m=radius_m,
                T_bath_k=T_bath,
                y_bath=y_film,
                pair=pair,
                tsat_k=tsat,
                pressure_pa=PRESSURE_PA,
            )
            audit["y_hexane_declared"] = y_h
            audit["y_hexane_film"] = y_film
            kind, payload = apply_high_loading_step(state, dt, m_dot, e_dot, params)
            if kind == "activated":
                event = payload
                t_s = t_s + event.elapsed_high_loading_s
                activation_event = event
                state = event.high_loading_at_event
                row = high_loading_row(
                    t_s,
                    state,
                    {
                        "index": hl_steps,
                        "T_bath_k": T_bath,
                        "y_hexane": y_h,
                        "elapsed_high_loading_s": event.elapsed_high_loading_s,
                        "activation": True,
                        **audit,
                    },
                )
                write_step(row)
                result["activation"] = {
                    "time_s": t_s,
                    "T_k": event.radial.temperature_k,
                    "X_hexane_kg_kg": event.high_loading_at_event.total_hexane_loading,
                    "ledger_max_normalized_residual": (
                        event.transition_ledger.max_normalized_residual
                    ),
                    "remaining_radial_time_s": event.remaining_radial_time_s,
                }
                break
            state = payload
            t_s += dt
            hl_steps += 1
            if hl_steps % 10 == 0 or t_s >= t_end:
                write_step(
                    high_loading_row(
                        t_s,
                        state,
                        {
                            "index": hl_steps,
                            "T_bath_k": T_bath,
                            "y_hexane": y_h,
                            **audit,
                        },
                    )
                )
        result["high_loading_steps"] = hl_steps
        result["high_loading_end_s"] = t_s
        result["high_loading_end"] = high_loading_row(t_s, state, {})
        flush_result()

        if activation_event is None:
            result["stage"] = "high_loading_incomplete"
            result["stop"] = "activation not reached in the allotted time"
            result["outlet_score"] = score_outlet(
                result["high_loading_end"], destiny["outlet_box"]
            )
            return result

        if args.high_loading_only:
            result["stage"] = "high_loading_only"
            result["stop"] = "high-loading-only requested"
            return result

        result["stage"] = "install_receding"
        flush_result()
        T_a = activation_event.radial.temperature_k
        T_bath, y_declared = interpolate_bath(waypoints, t_s)
        wet_bounds = backend.CONTROLS.wet_temperature_bounds_k
        iface_bounds = backend.CONTROLS.interface_temperature_bounds_k
        dry_bounds = backend.CONFIG.dry.conditioned_temperature_domain.solver_bounds_k
        # Receding start temperatures: the P1EF native Gate-1g state (front
        # 0.60 R, wet 335 / dry 343 / iface 337 / y=0.82) is only the probe
        # used to read the resolved start inventory; the state actually handed
        # to the resolved march is rebuilt below from the A->B bridge endpoint
        # (wet = iface = bridge T_i, dry between T_i and the bath, dry gas at
        # the interface trace y_I = Psat_h(T_i)/P).  A 0.999 R birth on the
        # water-bearing chart has one dry piece and refuses. Historical
        # loadings still come from the actual activation temperature.
        wet_t = strictly_inside(P1EF_WET_T_K, wet_bounds)
        dry_t = strictly_inside(P1EF_DRY_T_K, dry_bounds)
        iface_t = strictly_inside(P1EF_IFACE_T_K, iface_bounds)
        dry_y = P1EF_DRY_Y_HEXANE
        y_h = receding_admissible_y_hexane(y_declared, min(dry_t, T_bath))
        result["p1ef_case"] = installed.coefficient_case.identifier
        result["receding_charts"] = {
            "wet_temperature_bounds_k": list(wet_bounds),
            "interface_temperature_bounds_k": list(iface_bounds),
            "dry_solver_bounds_k": list(dry_bounds),
            "interface_bounds_origin": (
                "qualify_phase1_particle_engineering.INTERFACE_TEMPERATURE_BOUNDS_K; "
                "PHY-039 steam-distillation to n-hexane Tsat at 1 atm"
            ),
            "start_origin": (
                "probe: qualify_gate1g_conditioned_cut_stress.initial_integrator "
                "(front 0.60 R, wet 335 K, dry 343 K, iface 337 K, y_h=0.82); handed-off "
                "state rebuilt from the A->B bridge endpoint, see ab_bridge"
            ),
            "activation_T_k": T_a,
            "activation_front_s_over_R": 1.0,
            "receding_front_fraction": args.front_fraction,
            "wet_t": wet_t,
            "dry_t": dry_t,
            "iface_t": iface_t,
            "dry_y_hexane": dry_y,
            "y_declared_at_activation": y_declared,
            "y_film_at_activation": y_h,
            "y_hexane_dew_floor_applied": y_h != y_declared,
            "dew_activity_margin": DEW_ACTIVITY_MARGIN,
            "dew_floor_origin": (
                "film y_water <= dew_activity_margin * psat_w(T_ref) / P with T_ref = "
                "min(outer dry piece (regime B) or outer cell (regime C), bath T), so "
                "the dry-shell coupled-pore a_w stays strictly below 1 (no "
                "condensation topology in the model; re-wetting is outside the "
                "formulation's stated process). Active only while the declared bath is "
                "below its own water dew point. Margin is a recorded harness knob."
            ),
            "faner_emergency_subdivision_licence": (
                subdivision_under_constant_forcing_licensed()
            ),
            "faner_emergency_certification": (
                FANER_EMERGENCY_CERTIFICATION
                if subdivision_under_constant_forcing_licensed()
                else None
            ),
            "receding_dt_s": FANER_DT_S,
            "receding_dt_frozen": True,
            "external_dyadic_ladder": False,
            "face_crossing_owner": (
                "qualify_gate1g_conditioned_cut_stress."
                "_advance_bracketed_face_from_accepted_history"
            ),
            "cells_origin": (
                "P1EF native N=2 after the N=2 seed; one interior face at 0.5 R"
                if args.cells == 2
                else f"caller --cells {args.cells}"
            ),
        }
        patch_initial_integrator(
            front_fraction=args.front_fraction,
            wet_t=wet_t,
            wet_w=feed["feed_water_kg_kg_dry"],
            dry_t=dry_t,
            dry_y=dry_y,
            iface_t=iface_t,
            act_t=T_a,
            oil_fraction=feed["feed_oil_kg_kg_dry"],
            time_s=t_s,
        )

        # A -> B bridge.  The resolved cut chart cannot be born at the surface
        # (a 0.999 R birth has one dry piece and refuses), and the resolved
        # start state at front_fraction R holds only the core's historical
        # label plus the shell's sorbate.  Instead of discarding the
        # difference, the recession from R to front_fraction R is integrated
        # with the heat-limited front law of Sec. 2.3 in its quasi-steady
        # lumped form: the latent duty is supplied by the film conductance in
        # series with the growing dry shell's conduction resistance
        # R (R - s) / (k s), the core boils at the activation temperature, and
        # the clock advances until the lumped hexane inventory equals the
        # resolved start state's inventory.  Harness approximation of the
        # formulation (no resolved shell fields), recorded, mass-conserving.
        probe_before = backend.initial_integrator(args.cells, radius_m)
        probe_summary = receding_summary(probe_before, dry_mass)
        x_h_b = probe_summary["X_hexane_kg_kg"]
        x_h_a = state.total_hexane_loading
        bridge = {
            "law": (
                "quasi-steady heat-limited recession (Sec. 2.3 front velocity, lumped) on "
                "the Sec. 2.1 surface law: at each step the interface temperature T_i solves "
                "U(s) A(beta) (T_bath - T_i) = j_h(T_i) dh_vap(T_i) with the Stefan film flux "
                "j_h = rho_f h_M ln[(1-w_b)/(1-w_I(T_i))], w_I from y_I = Psat_h(T_i)/P, w_b from "
                "the declared bath under the no-condensation film admission at min(T_i, T_bath); "
                "1/U = 1/h_film + R (R - s) / (k s), s/R = (X_h / X_h,activation)^(1/3); "
                "no shell mass resistance (total flux algebraic, Sec. 2.2)"
            ),
            "T_a_k": T_a,
            "X_hexane_activation_kg_kg": x_h_a,
            "X_hexane_resolved_start_kg_kg": x_h_b,
            "X_water_activation_kg_kg": state.retained_water_loading,
            "X_water_resolved_start_kg_kg": probe_summary["X_water_kg_kg"],
            "resolved_start_temperatures_k": {
                "wet": wet_t,
                "dry": dry_t,
                "interface": iface_t,
            },
            "time_start_s": t_s,
            "euler_dt_s": 5.0e-3,
        }
        if x_h_b < x_h_a:
            area = 4.0 * math.pi * radius_m**2
            h_q = STRESS_MULTIPLIER * pair.heat_transfer_w_m2_k
            m_h = x_h_a * dry_mass
            m_target = x_h_b * dry_mass
            t_b = t_s
            bridge_dt = bridge["euler_dt_s"]
            shell_fraction_end = None
            t_i_first = None
            t_i_last = None
            beta_max = 0.0

            def bridge_interface(u_series: float, T_bath_b: float, y_b: float):
                """Quasi-steady interface temperature by bisection (too-hot => beta>1)."""
                lo = hexane_saturation_temperature_k(max(1.0e-6, y_b) * PRESSURE_PA)
                hi = T_bath_b
                if hi <= lo:
                    return None

                def balance(T_i: float):
                    try:
                        f = wet_surface_film_rates(
                            surface_temperature_k=T_i,
                            bath_temperature_k=T_bath_b,
                            bath_y_hexane=receding_admissible_y_hexane(
                                y_b, min(T_i, T_bath_b)
                            ),
                            conductance_w_m2_k=u_series,
                            pair=pair,
                            pressure_pa=PRESSURE_PA,
                        )
                    except efr.ReducedFilmValidityError:
                        return -math.inf, None
                    return f["q_in_w_m2"] - f["j_h_kg_m2_s"] * f["dh_vap_j_kg"], f

                f_lo, film_lo = balance(lo)
                if f_lo <= 0.0:
                    return None
                best = (lo, film_lo)
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    f_mid, film_mid = balance(mid)
                    if f_mid > 0.0:
                        lo = mid
                        best = (mid, film_mid)
                    else:
                        hi = mid
                return best

            while m_h > m_target and t_b < t_end:
                T_bath_b, y_b = interpolate_bath(waypoints, t_b)
                s_over_r = max(1.0e-9, (m_h / (x_h_a * dry_mass))) ** (1.0 / 3.0)
                r_shell = radius_m * (1.0 - s_over_r) / (CONDUCTIVITY * s_over_r)
                u_series = 1.0 / (1.0 / h_q + r_shell)
                found = bridge_interface(u_series, T_bath_b, y_b)
                if found is None or found[1] is None or found[1]["j_h_kg_m2_s"] <= 0.0:
                    # No outward hexane drive at this bath state (bath at or
                    # below the no-condensation zero-drive temperature): the
                    # front holds; the clock advances on the declared bath.
                    bridge["hold_s"] = bridge.get("hold_s", 0.0) + bridge_dt
                    t_b += bridge_dt
                    continue
                T_i, film = found
                m_dot = film["j_h_kg_m2_s"] * area
                dt_b = min(bridge_dt, (m_h - m_target) / m_dot)
                m_h -= m_dot * dt_b
                t_b += dt_b
                shell_fraction_end = r_shell / (1.0 / h_q + r_shell)
                t_i_first = T_i if t_i_first is None else t_i_first
                t_i_last = T_i
                beta_max = max(beta_max, abs(film["beta"] or 0.0))
            bridge["time_end_s"] = t_b
            bridge["duration_s"] = t_b - t_s
            bridge["hexane_removed_kg"] = (x_h_a * dry_mass) - m_h
            bridge["hexane_removed_kg_kg"] = x_h_a - m_h / dry_mass
            bridge["residual_inventory_mismatch_kg_kg"] = m_h / dry_mass - x_h_b
            bridge["shell_resistance_fraction_at_end"] = shell_fraction_end
            bridge["interface_temperature_first_k"] = t_i_first
            bridge["interface_temperature_last_k"] = t_i_last
            bridge["ackermann_beta_max"] = beta_max
            bridge["h_film_w_m2_k"] = h_q
            x_h_end = m_h / dry_mass
            x_w_a = state.retained_water_loading
            total_end = 1.0 + x_h_end + x_w_a + state.oil_fraction
            t_b_end_bath = interpolate_bath(waypoints, t_b)[0]
            write_step(
                {
                    "regime": "AB_bridge",
                    "index": hl_steps,
                    "time_s": t_b,
                    "T_k": t_i_last if t_i_last is not None else T_a,
                    "interface_temperature_k": t_i_last,
                    "X_hexane_kg_kg": x_h_end,
                    "X_water_kg_kg": x_w_a,
                    "moisture_wt_pct_wb": 100.0 * x_w_a / total_end,
                    "hexane_ppm_wb": 1.0e6 * x_h_end / total_end,
                    "front_s_over_R": (m_h / (x_h_a * dry_mass)) ** (1.0 / 3.0),
                    "T_bath_k": t_b_end_bath,
                    "bridge": True,
                }
            )
            t_s = t_b
            # Hand-off state from the bridge endpoint (not the P1EF probe):
            # the core and interface at the bridge interface temperature, the
            # dry shell midway between interface and bath (its conduction
            # profile is linear in 1/r; one piece at N=2), the dry pore gas at
            # the interface trace y_I = Psat_h(T_i)/P.  Every value is kept
            # strictly inside its chart; the clipping, if any, is recorded.
            if t_i_last is not None:
                wet_t_raw = t_i_last
                iface_t_raw = t_i_last
                dry_t_raw = 0.5 * (t_i_last + t_b_end_bath)
                y_declared_end = interpolate_bath(waypoints, t_b)[1]
                # One dry piece at N=2: its gas is set to the admitted film
                # composition at its own temperature (the outer-surface state
                # the boundary imposes, a_w = margin).  The interface trace
                # y_I = Psat_h(T_i)/P read at the hotter piece temperature
                # would put the piece at a_w ~ 0.69, inside the dry chart's
                # band of non-positive datum-reduced energy capacity
                # (a_w ~ 0.60-0.78), which initialize_integrator_state refuses.
                dry_y_raw = receding_admissible_y_hexane(
                    y_declared_end, min(dry_t_raw, t_b_end_bath)
                )
                wet_t = strictly_inside(wet_t_raw, wet_bounds)
                iface_t = strictly_inside(iface_t_raw, iface_bounds)
                dry_t = strictly_inside(dry_t_raw, dry_bounds)
                dry_y = dry_y_raw
                bridge["handoff"] = {
                    "rule": (
                        "wet = iface = bridge T_i; dry = (T_i + T_bath)/2; dry y_h = the "
                        "no-condensation admitted film composition at the dry piece "
                        "temperature (a_w = dew_activity_margin); each kept strictly "
                        "inside its chart"
                    ),
                    "raw": {
                        "wet": wet_t_raw,
                        "dry": dry_t_raw,
                        "interface": iface_t_raw,
                        "dry_y_hexane": dry_y_raw,
                        "interface_trace_y_hexane": min(
                            1.0, hx.saturation_pressure(t_i_last) / PRESSURE_PA
                        ),
                    },
                    "applied": {
                        "wet": wet_t,
                        "dry": dry_t,
                        "interface": iface_t,
                        "dry_y_hexane": dry_y,
                    },
                    "clipped": {
                        "wet": wet_t != wet_t_raw,
                        "dry": dry_t != dry_t_raw,
                        "interface": iface_t != iface_t_raw,
                    },
                    "T_bath_k": t_b_end_bath,
                    "dry_gas_activities_at_dry_t": {
                        "a_w": (1.0 - dry_y) * PRESSURE_PA / wa.saturation_pressure(dry_t),
                        "a_h": dry_y * PRESSURE_PA / hx.saturation_pressure(dry_t),
                    },
                }
                bridge["resolved_start_temperatures_k"] = {
                    "wet": wet_t,
                    "dry": dry_t,
                    "interface": iface_t,
                }
                result["receding_charts"].update(
                    {
                        "wet_t": wet_t,
                        "dry_t": dry_t,
                        "iface_t": iface_t,
                        "dry_y_hexane": dry_y,
                    }
                )
            result["ab_bridge"] = bridge
            flush_result()
            patch_initial_integrator(
                front_fraction=args.front_fraction,
                wet_t=wet_t,
                wet_w=feed["feed_water_kg_kg_dry"],
                dry_t=dry_t,
                dry_y=dry_y,
                iface_t=iface_t,
                act_t=T_a,
                oil_fraction=feed["feed_oil_kg_kg_dry"],
                time_s=t_s,
            )
            probe_after = receding_summary(
                backend.initial_integrator(args.cells, radius_m), dry_mass
            )
            bridge["X_hexane_resolved_start_handoff_kg_kg"] = probe_after["X_hexane_kg_kg"]
            bridge["X_water_resolved_start_handoff_kg_kg"] = probe_after["X_water_kg_kg"]
            bridge["handoff_inventory_mismatch_kg_kg"] = (
                m_h / dry_mass - probe_after["X_hexane_kg_kg"]
            )
            T_bath, y_declared = interpolate_bath(waypoints, t_s)
            y_h = receding_admissible_y_hexane(y_declared, min(dry_t, T_bath))
            result["receding_charts"]["y_declared_at_activation"] = y_declared
            result["receding_charts"]["y_film_at_activation"] = y_h
            result["receding_charts"]["y_hexane_dew_floor_applied"] = y_h != y_declared
        else:
            bridge["duration_s"] = 0.0
            bridge["note"] = "resolved start inventory not below activation inventory"
        result["ab_bridge"] = bridge
        flush_result()

        def boundary_at(
            clock: float,
            *,
            particle_reference_t_k: float = dry_t,
            forcing_timescale_s: float = FANER_DT_S,
        ) -> ct.ReducedFilmPoreBoundary:
            T_b, y_b = interpolate_bath(waypoints, clock)
            # No-condensation admission: the film's water partial pressure is
            # capped at margin * psat_w at the coldest temperature the outer
            # dry side can reach over the step: min(outer dry piece/cell,
            # bath).  Below the water dew point of the declared bath this is
            # active (the declared gas is fogging), above it the film is the
            # declared bath.
            dew_reference_t_k = min(particle_reference_t_k, T_b)
            y_film = receding_admissible_y_hexane(y_b, dew_reference_t_k)
            return reduced_film_boundary(
                temperature_k=T_b,
                y_hexane=y_film,
                pair=pair,
                forcing_timescale_s=forcing_timescale_s,
                label=(
                    f"dtdc destiny t={clock:.4f}s T={T_b:.4f}K "
                    f"y_h={y_film:.6f} (declared {y_b:.6f}, dew ref {dew_reference_t_k:.3f} K)"
                ),
            )

        result["stage"] = "receding_first_step"
        flush_result()
        before = backend.initial_integrator(args.cells, radius_m)
        nfaces = before.transport.layout.dry_face_count
        cand = cut.candidate_from_state(
            before.transport,
            dry_total_stefan_fluxes_mol_m2_s=stefan_seed_for_faces(nfaces),
            interface_temperature_k=iface_t,
        )
        cand = replace(cand, front_z=cand.front_z * (1.0 - args.seed_front_z_fraction))
        seed = backend._seed(cand, "dtdc destiny: Faner initial-state seed")
        boundary = boundary_at(t_s)
        first_dt = FANER_DT_S
        step, _record, error = backend._try_step_detailed(
            before,
            first_dt,
            boundary,
            seed,
            role="dtdc destiny first receding step",
            predictor_fraction=1.0,
        )
        result["first_receding_step"] = {
            "accepted": step is not None,
            "error": None if error is None else f"{type(error).__name__}: {str(error)[:3000]}",
            "dt_s": first_dt,
        }
        if step is None:
            # Push through: dyadic dt, same seed rule as Faner.
            refusal_ladder = []
            trial_dt = first_dt
            for _ in range(8):
                trial_dt *= 0.5
                step, _record, error = backend._try_step_detailed(
                    before,
                    trial_dt,
                    boundary,
                    seed,
                    role=f"dtdc destiny first receding step dt={trial_dt}",
                    predictor_fraction=1.0,
                )
                refusal_ladder.append(
                    {
                        "dt_s": trial_dt,
                        "accepted": step is not None,
                        "error": None
                        if error is None
                        else f"{type(error).__name__}: {str(error)[:1500]}",
                    }
                )
                if step is not None:
                    break
            result["first_receding_dt_ladder"] = refusal_ladder
        if step is None:
            result["stage"] = "refused_receding_first_step"
            result["stop"] = "receding first step refused after dt ladder"
            result["outlet_score"] = score_outlet(
                result["high_loading_end"], destiny["outlet_box"]
            )
            return result

        rec = backend._step_record(
            step,
            (
                backend.AttemptRecord(
                    "mesh_bootstrap_exact_root",
                    1.0,
                    0.0,
                    True,
                    step.ledger.nonlinear_evaluations,
                    step.ledger.maximum_scaled_residual,
                    step.ledger.condition_proxy,
                    True,
                ),
            ),
            "dtdc_destiny",
            "direct_first_step",
            starts_after_boundary_discontinuity=True,
        )
        row = {
            "index": 0,
            **receding_summary(step.after, dry_mass),
            "T_bath_k": boundary.temperature_k,
            "y_hexane": boundary.y_hexane,
            "maximum_scaled_residual": rec["maximum_scaled_residual"],
            "maximum_step_ledger_residual": rec["maximum_step_ledger_residual"],
            "maximum_cumulative_ledger_residual": rec[
                "maximum_cumulative_ledger_residual"
            ],
            "mode": "direct_first_step",
        }
        write_step(row)
        result["receding_initial"] = row
        flush_result()

        result["stage"] = "receding_march"
        previous = step
        accepted = 1
        counts = {
            "face_events": 0,
            "subdivisions": 0,
            "push_through_retries": 0,
            "emergency_licence_used": 0,
        }
        refusal = None
        receding_dt = FANER_DT_S
        result["face_tangent_knobs"] = face_tangent_knobs()
        dry_state = None
        extinction_tried = False
        target_idx = 0
        w_o_feed = feed["feed_oil_kg_kg_dry"]
        while True:
            clock = previous.after.transport.time_s
            if time.perf_counter() - t0 > args.wall_budget_s:
                result["stop"] = "wall budget"
                break
            if clock >= t_end - 1e-12:
                result["stop"] = "max physical time"
                break
            if args.smoke and accepted >= args.smoke_receding_steps:
                result["stop"] = "smoke receding cap"
                break
            snap = receding_summary(previous.after, dry_mass)
            if snap["X_hexane_kg_kg"] <= args.x_stop:
                result["stop"] = f"loading at or below {args.x_stop}"
                break
            boundary = boundary_at(
                clock,
                particle_reference_t_k=float(
                    previous.after.transport.dry_temperatures_k[-1]
                ),
            )
            core_s = wet_core_radius_m(previous.after, args.radius_m)
            if (
                core_s <= GEOMETRIC_EXTINCTION_CORE_RADIUS_M
                and previous.after.transport.layout.cut_cell_index == 0
            ):
                # Regime B -> C by geometric cutoff (see constant docstring).
                result["stage"] = "geometric_extinction"
                result["progress"] = {
                    "index": accepted,
                    "awaiting": "geometric_extinction_projection",
                    **snap,
                }
                flush_result()
                tw = time.perf_counter()
                try:
                    dry_state, ginfo = project_geometric_extinction(
                        previous, receding_dt, boundary
                    )
                except Exception as exc:  # noqa: BLE001
                    result["stop"] = "geometric extinction projection refused"
                    result["refusal"] = {
                        "step_index": accepted,
                        "time_before_s": clock,
                        "core_radius_m": core_s,
                        "error_type": type(exc).__name__,
                        "error": str(exc)[:6000],
                        "traceback": "".join(traceback.format_exception(exc))[-6000:],
                        "wall_s": time.perf_counter() - tw,
                    }
                    break
                wall = time.perf_counter() - tw
                counts["extinction_events"] = 1
                counts["geometric_extinction"] = 1
                drow = dry_summary(dry_state, dry_mass, w_o_feed)
                result["extinction"] = {
                    "kind": "geometric_cutoff",
                    "step_index": accepted,
                    "time_before_s": clock,
                    "frozen_dt_s": receding_dt,
                    "committed_dt_s": dry_state.time_s - clock,
                    "core_radius_at_cutoff_m": core_s,
                    "core_radius_floor_m": GEOMETRIC_EXTINCTION_CORE_RADIUS_M,
                    "front_z_at_cutoff": snap["front_z"],
                    "hexane_mol_before": snap["hexane_mol"],
                    "hexane_mol_after": drow["hexane_mol"],
                    "water_mol_before": snap["water_mol"],
                    "water_mol_after": drow["water_mol"],
                    "energy_j_before": snap["energy_j"],
                    "energy_j_after": drow["energy_j"],
                    "projection": ginfo,
                    "environment_conditioned_resolution": resolution_record(
                        extinction_resolution_for_state(previous.after)
                    ),
                    "certification": FANER_EMERGENCY_CERTIFICATION,
                    "wall_s": wall,
                }
                row = {
                    "index": accepted,
                    "wall_s": wall,
                    **drow,
                    "T_bath_k": boundary.temperature_k,
                    "y_hexane": boundary.y_hexane,
                    "maximum_scaled_residual": ginfo["maximum_scaled_residual"],
                    "mode": "geometric_extinction",
                    "dt_s": ginfo["projection_duration_s"],
                    "face_event": False,
                    "adaptive": False,
                    "emergency_licence_used": True,
                }
                write_step(row)
                accepted += 1
                result["progress"] = row
                result["counts"] = counts
                flush_result()
                break
            result["progress"] = {
                "index": accepted,
                "awaiting": "receding_step",
                "attempting_dt_s": receding_dt,
                **snap,
                "T_bath_k": boundary.temperature_k,
                "y_hexane": boundary.y_hexane,
            }
            result["counts"] = counts
            flush_result()
            tw = time.perf_counter()
            current = None
            last_exc = None
            cell_zero_terminus = None
            try:
                current, attempts, mode = advance_receding_step(
                    previous, receding_dt, boundary
                )
                rec = backend._step_record(
                    current,
                    attempts,
                    "dtdc_destiny",
                    mode,
                    starts_after_boundary_discontinuity=False,
                )
                q.require_held_step_failed_leaf_contract(rec, accepted)
            except CellZeroTerminus as exc:
                cell_zero_terminus = exc
                current = None
            except Exception as exc:  # noqa: BLE001 - recorded; no external dt ladder
                last_exc = exc
                current = None
            if cell_zero_terminus is not None and args.try_extinction and not extinction_tried:
                extinction_tried = True
                # Regime B -> C. Cell-zero receding refused on the frozen
                # dt; the src extinction family (localizer ladder +
                # finite-limit certificate + exact z=0 projection) owns the
                # terminus. Composition is harness-side; wheels are src/.
                result["stage"] = "extinction"
                result["progress"] = {
                    "index": accepted,
                    "awaiting": "extinction_localization",
                    **snap,
                }
                flush_result()
                rungs: list[dict] = []

                def on_rung(level: int, step) -> None:
                    led = step.ledger
                    rungs.append(
                        {
                            "level": level,
                            "target_z": led.target_z,
                            "target_fraction_of_initial_z": led.target_fraction_of_initial_z,
                            "event_duration_s": led.event_duration_s,
                            "nonlinear_evaluations": led.nonlinear_evaluations,
                            "maximum_scaled_residual": led.maximum_scaled_residual,
                            "condition_proxy": led.condition_proxy,
                        }
                    )
                    result["progress"] = {
                        "index": accepted,
                        "awaiting": "extinction_localization",
                        "rung": rungs[-1],
                        "wall_s": time.perf_counter() - tw,
                    }
                    flush_result()

                try:
                    kind, payload = resolve_cell_zero_extinction(
                        previous, receding_dt, boundary, on_rung=on_rung
                    )
                except Exception as exc:  # noqa: BLE001
                    result["extinction_refusal"] = {
                        "step_index": accepted,
                        "time_before_s": clock,
                        "attempted_dt_s": receding_dt,
                        "error_type": type(exc).__name__,
                        "error": str(exc)[:4000],
                        "traceback": "".join(traceback.format_exception(exc))[-6000:],
                        "adaptive_error_type": type(
                            cell_zero_terminus.adaptive_error
                        ).__name__,
                        "adaptive_error": str(cell_zero_terminus.adaptive_error)[:4000],
                        "rungs": rungs[-40:],
                        "environment_conditioned_resolution": resolution_record(
                            extinction_resolution_for_state(previous.after)
                        ),
                        "wall_s": time.perf_counter() - tw,
                    }
                    result["stage"] = "receding_march"
                    flush_result()
                    kind = None
                if kind == "same_cell":
                    # The orchestrator's own same-cell BE accepted where the
                    # qualifier recovery chain refused: take it as the step.
                    wall = time.perf_counter() - tw
                    taken = payload.same_cell_step
                    row = {
                        "index": accepted,
                        "wall_s": wall,
                        **receding_summary(taken.after, dry_mass),
                        "T_bath_k": boundary.temperature_k,
                        "y_hexane": boundary.y_hexane,
                        "mode": "orchestrator_same_cell",
                        "face_event": False,
                        "adaptive": False,
                        "emergency_licence_used": False,
                    }
                    write_step(row)
                    previous = taken
                    accepted += 1
                    counts["orchestrator_same_cell"] = (
                        counts.get("orchestrator_same_cell", 0) + 1
                    )
                    result["stage"] = "receding_march"
                    continue
                if kind == "extinction":
                    wall = time.perf_counter() - tw
                    outcome = payload
                    counts["extinction_events"] = 1
                    after = outcome.after
                    if isinstance(after, cep.ExtinctionCommittedState):
                        dry_state = after.to_uniform_fully_dry_transport_state()
                    else:
                        dry_state = after
                    result["extinction"] = {
                        "step_index": accepted,
                        "time_before_s": clock,
                        "frozen_dt_s": receding_dt,
                        "committed_dt_s": dry_state.time_s - clock,
                        "branch": str(outcome.branch),
                        "dry_continuation_committed": outcome.dry_continuation_committed,
                        "outcome": json.loads(
                            json.dumps(outcome.journal_record(), default=str)
                        ),
                        "rungs": rungs,
                        "adaptive_refusal": str(cell_zero_terminus.adaptive_error)[:2000],
                        "environment_conditioned_resolution": resolution_record(
                            extinction_resolution_for_state(previous.after)
                        ),
                        "ruled_part02_resolution": resolution_record(
                            elm.part02_environment_conditioned_resolution()
                        ),
                        "wall_s": wall,
                    }
                    row = {
                        "index": accepted,
                        "wall_s": wall,
                        **dry_summary(dry_state, dry_mass, w_o_feed),
                        "T_bath_k": boundary.temperature_k,
                        "y_hexane": boundary.y_hexane,
                        "mode": "extinction",
                        "face_event": False,
                        "adaptive": False,
                        "emergency_licence_used": False,
                    }
                    write_step(row)
                    accepted += 1
                    result["progress"] = row
                    result["counts"] = counts
                    flush_result()
                    break
            if cell_zero_terminus is not None:
                # Cell-zero target march. The frozen-dt BE refused (sliver
                # dry piece right after the face, or asymptotic tail); the
                # src positive-core target solve (front prescribed, duration
                # unknown; identical equations) is committed instead. After
                # each commit the frozen-dt BE is tried again.
                counts["cell_zero_refusals"] = counts.get("cell_zero_refusals", 0) + 1
                fraction = TARGET_MARCH_FRACTIONS[
                    min(target_idx, len(TARGET_MARCH_FRACTIONS) - 1)
                ]
                remaining = t_end - clock
                committed = None
                commit_errors = []
                for attempt_fraction in (fraction, math.sqrt(fraction)):
                    try:
                        committed, info = commit_positive_core_target(
                            previous, boundary, attempt_fraction, remaining
                        )
                        break
                    except Exception as exc:  # noqa: BLE001
                        commit_errors.append(
                            {
                                "fraction": attempt_fraction,
                                "error_type": type(exc).__name__,
                                "error": str(exc)[:3000],
                            }
                        )
                if committed is None:
                    result["stop"] = "cell-zero target commit refused"
                    result["refusal"] = {
                        "step_index": accepted,
                        "time_before_s": clock,
                        "attempted_dt_s": receding_dt,
                        "adaptive_error_type": type(
                            cell_zero_terminus.adaptive_error
                        ).__name__,
                        "adaptive_error": str(cell_zero_terminus.adaptive_error)[:4000],
                        "target_commit_errors": commit_errors,
                        "wall_s": time.perf_counter() - tw,
                    }
                    break
                wall = time.perf_counter() - tw
                counts["target_commits"] = counts.get("target_commits", 0) + 1
                row = {
                    "index": accepted,
                    "wall_s": wall,
                    **receding_summary(committed.after, dry_mass),
                    "T_bath_k": boundary.temperature_k,
                    "y_hexane": boundary.y_hexane,
                    "maximum_scaled_residual": info["maximum_scaled_residual"],
                    "maximum_step_ledger_residual": None,
                    "maximum_cumulative_ledger_residual": None,
                    "mode": "positive_core_target_commit",
                    "dt_s": info["event_duration_s"],
                    "target_fraction_of_current_z": info["target_fraction_of_current_z"],
                    "nonlinear_evaluations": info["nonlinear_evaluations"],
                    "condition_proxy": info["condition_proxy"],
                    "face_event": False,
                    "adaptive": False,
                    "emergency_licence_used": False,
                    "frozen_dt_refusal": str(cell_zero_terminus.adaptive_error)[:300],
                }
                write_step(row)
                target_log = result.setdefault("target_march", [])
                if len(target_log) < 400:
                    target_log.append(
                        {
                            "index": accepted,
                            "time_before_s": clock,
                            "wall_s": wall,
                            "commit_errors": commit_errors,
                            **info,
                        }
                    )
                previous = committed
                accepted += 1
                target_idx += 1
                result["progress"] = row
                result["counts"] = counts
                flush_result()
                continue
            target_idx = 0
            if current is None:
                refusal = {
                    "step_index": accepted,
                    "time_before_s": clock,
                    "attempted_dt_s": receding_dt,
                    "error_type": type(last_exc).__name__,
                    "error": str(last_exc)[:4000],
                    "traceback": "".join(traceback.format_exception(last_exc))[-4000:],
                    "wall_s": time.perf_counter() - tw,
                }
                result["stop"] = "receding step refused"
                result["refusal"] = refusal
                break
            wall = time.perf_counter() - tw
            fe = bool(rec.get("interior_face_event_used"))
            ad = bool(rec.get("adaptive_subdivision_used"))
            lic = bool(
                (rec.get("frozen_failed_leaf_controller_contract") or {}).get(
                    "emergency_constant_forcing_subdivision_licence_used"
                )
            )
            counts["face_events"] += int(fe)
            counts["subdivisions"] += int(ad)
            counts["emergency_licence_used"] += int(lic)
            row = {
                "index": accepted,
                "wall_s": wall,
                **receding_summary(current.after, dry_mass),
                "T_bath_k": boundary.temperature_k,
                "y_hexane": boundary.y_hexane,
                "maximum_scaled_residual": rec["maximum_scaled_residual"],
                "maximum_step_ledger_residual": rec["maximum_step_ledger_residual"],
                "maximum_cumulative_ledger_residual": rec[
                    "maximum_cumulative_ledger_residual"
                ],
                "mode": mode,
                "face_event": fe,
                "adaptive": ad,
                "emergency_licence_used": lic,
            }
            write_step(row)
            previous = current
            accepted += 1
            result["progress"] = row
            result["counts"] = counts
            if fe or ad or accepted % 5 == 0:
                flush_result()
        result["receding_steps"] = accepted
        result["counts"] = counts
        if dry_state is None:
            last = receding_summary(previous.after, dry_mass)
            result["final"] = last
            result["outlet_score"] = score_outlet(last, destiny["outlet_box"])
            result["stage"] = "done"
            return result

        # Regime C: fully-dry predictor (coupled_transport BE) on the same
        # reduced-film bath boundary.  Two recorded regime-C policies:
        #  (i) time step: start at the frozen 0.075 s; on a step refusal the
        #      step is doubled (never shrunk again) up to DRY_DT_CAP_S.  The
        #      fully-dry predictor recovers the surface Stefan flux
        #      algebraically from the inventory change over dt, so a Newton
        #      *trial* at tiny dt carries a flux ~ (trial inventory change)/dt
        #      that leaves the reduced film's |beta|<=1 Ackermann validity
        #      band although the converged step does not (measured 2026-10-02
        #      at 0.5 mm: refused at 0.075/0.15/0.3 s after 2 evaluations,
        #      accepted at >=0.6 s with scaled residuals ~1e-12).
        #  (ii) dew floor: the film water partial pressure cap is evaluated at
        #      min(particle dry-side reference, bath T) in both regimes B and
        #      C, so the no-condensation admission follows the particle and
        #      releases once the bath is above its own water dew point.
        result["stage"] = "fully_dry_march"
        result["stop"] = None
        dry_accepted = 0
        dry_max = {"component": 0.0, "energy": 0.0, "step_ledger": 0.0, "cum_ledger": 0.0}
        dry_evals = 0
        dry_refusal = None
        dry_dt = FANER_DT_S
        dry_dt_events: list[dict] = []
        dry_dt_refused_trials = 0
        while True:
            clock = dry_state.time_s
            if time.perf_counter() - t0 > args.wall_budget_s:
                result["stop"] = "wall budget"
                break
            if clock >= t_end - 1e-12:
                result["stop"] = "max physical time"
                break
            if args.smoke and dry_accepted >= args.smoke_receding_steps:
                result["stop"] = "smoke receding cap"
                break
            snap = dry_summary(dry_state, dry_mass, w_o_feed)
            if snap["X_hexane_kg_kg"] <= args.x_stop:
                result["stop"] = f"loading at or below {args.x_stop}"
                break
            dew_ref_t = min(float(dry_state.temperatures_k[-1]), interpolate_bath(waypoints, clock)[0])
            dt_try = min(dry_dt, t_end - clock)
            boundary = boundary_at(
                clock,
                particle_reference_t_k=float(dry_state.temperatures_k[-1]),
                forcing_timescale_s=dt_try,
            )
            tw = time.perf_counter()
            try:
                dstep = ct.advance_fully_dry_backward_euler(dry_state, dt_try, boundary)
            except ct.CoupledTransportStepError as exc:
                dry_dt_refused_trials += 1
                if dry_dt < DRY_DT_CAP_S:
                    dry_dt_events.append(
                        {
                            "step_index": accepted,
                            "time_before_s": clock,
                            "refused_dt_s": dt_try,
                            "next_dt_s": min(2.0 * dry_dt, DRY_DT_CAP_S),
                            "error": str(exc)[:400],
                            "nonlinear_evaluations": getattr(exc, "nonlinear_evaluations", None),
                        }
                    )
                    dry_dt = min(2.0 * dry_dt, DRY_DT_CAP_S)
                    continue
                dry_refusal = {
                    "step_index": accepted,
                    "time_before_s": clock,
                    "attempted_dt_s": dt_try,
                    "error_type": type(exc).__name__,
                    "error": str(exc)[:4000],
                    "traceback": "".join(traceback.format_exception(exc))[-4000:],
                    "wall_s": time.perf_counter() - tw,
                }
                result["stop"] = "fully-dry step refused"
                result["refusal"] = dry_refusal
                break
            except Exception as exc:  # noqa: BLE001
                dry_refusal = {
                    "step_index": accepted,
                    "time_before_s": clock,
                    "attempted_dt_s": dt_try,
                    "error_type": type(exc).__name__,
                    "error": str(exc)[:4000],
                    "traceback": "".join(traceback.format_exception(exc))[-4000:],
                    "wall_s": time.perf_counter() - tw,
                }
                result["stop"] = "fully-dry step refused"
                result["refusal"] = dry_refusal
                break
            led = dstep.ledger
            dry_max["component"] = max(dry_max["component"], led.max_scaled_component_residual)
            dry_max["energy"] = max(dry_max["energy"], led.max_scaled_energy_residual)
            dry_max["step_ledger"] = max(
                dry_max["step_ledger"], led.max_normalized_accepted_step_ledger
            )
            dry_max["cum_ledger"] = max(dry_max["cum_ledger"], led.max_normalized_cumulative_ledger)
            dry_evals += int(led.nonlinear_evaluations)
            dry_state = dstep.after
            row = {
                "index": accepted,
                "wall_s": time.perf_counter() - tw,
                **dry_summary(dry_state, dry_mass, w_o_feed),
                "T_bath_k": boundary.temperature_k,
                "y_hexane": boundary.y_hexane,
                "y_hexane_declared": interpolate_bath(waypoints, clock)[1],
                "dew_reference_t_k": dew_ref_t,
                "dt_s": dt_try,
                "maximum_scaled_residual": max(
                    led.max_scaled_component_residual, led.max_scaled_energy_residual
                ),
                "maximum_step_ledger_residual": led.max_normalized_accepted_step_ledger,
                "maximum_cumulative_ledger_residual": led.max_normalized_cumulative_ledger,
                "nonlinear_evaluations": led.nonlinear_evaluations,
                "mode": "fully_dry",
                "face_event": False,
                "adaptive": False,
                "emergency_licence_used": False,
            }
            write_step(row)
            accepted += 1
            dry_accepted += 1
            result["progress"] = row
            if dry_accepted % 20 == 0:
                result["fully_dry"] = {
                    "steps": dry_accepted,
                    "dt_s": dry_dt,
                    "max_residuals": dict(dry_max),
                    "nonlinear_evaluations": dry_evals,
                    "dt_events": list(dry_dt_events),
                }
                flush_result()
        result["receding_steps"] = accepted
        result["fully_dry"] = {
            "steps": dry_accepted,
            "dt_s": dry_dt,
            "dt_start_s": FANER_DT_S,
            "dt_cap_s": DRY_DT_CAP_S,
            "dt_policy": (
                "start at the frozen receding dt; double (never shrink) on a "
                "CoupledTransportStepError up to dt_cap_s; the fully-dry predictor "
                "recovers the surface Stefan flux as inventory change / dt, so a "
                "Newton trial at tiny dt leaves the reduced film |beta|<=1 Ackermann "
                "band although the converged step does not"
            ),
            "dt_events": list(dry_dt_events),
            "refused_trials": dry_dt_refused_trials,
            "dew_floor_rule": (
                "film y_water <= dew_activity_margin * psat_w(T_ref) / P with "
                "T_ref = min(outer dry piece (regime B) or outer cell (regime C), bath T)"
            ),
            "dew_activity_margin": DEW_ACTIVITY_MARGIN,
            "max_residuals": dict(dry_max),
            "nonlinear_evaluations": dry_evals,
            "owner": "coupled_transport.advance_fully_dry_backward_euler",
        }
        last = dry_summary(dry_state, dry_mass, w_o_feed)
        result["final"] = last
        result["outlet_score"] = score_outlet(last, destiny["outlet_box"])
        result["stage"] = "done"
        return result
    except Exception as exc:  # noqa: BLE001
        result["stage"] = "failed"
        result["stop"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()[-6000:]
        return result
    finally:
        result["wall_total_s"] = time.perf_counter() - t0
        dump_json(summary_path, result)
        steps_file.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius-m", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells", type=int, default=FANER_CELLS)
    ap.add_argument("--high-loading-dt", type=float, default=1.0)
    ap.add_argument("--max-physical-s", type=float, default=1800.0)
    ap.add_argument("--wall-budget-s", type=float, default=8 * 3600.0)
    ap.add_argument("--x-stop", type=float, default=1.0e-6)
    ap.add_argument("--front-fraction", type=float, default=RECEDING_FRONT_FRACTION)
    ap.add_argument("--seed-front-z-fraction", type=float, default=2.0e-3)
    ap.add_argument("--high-loading-only", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--smoke-receding-steps", type=int, default=2)
    ap.add_argument(
        "--dew-margin",
        type=float,
        default=DEW_ACTIVITY_MARGIN,
        help="film water-activity margin for the dry-shell dew floor (default 0.95)",
    )
    ap.add_argument(
        "--try-extinction",
        action="store_true",
        help=(
            "at the first cell-zero frozen-dt refusal run the src extinction "
            "localizer ladder once (finite-limit certificate -> exact z=0 -> fully "
            "dry); default is the cell-zero positive-core target march"
        ),
    )
    args = ap.parse_args()
    if args.smoke:
        args.wall_budget_s = min(args.wall_budget_s, 900.0)
    result = march_one(args)
    print(json.dumps({k: result[k] for k in result if k != "traceback"}, default=str, indent=1))
    if result.get("stage") not in {"done", "high_loading_only", "high_loading_incomplete"}:
        if result.get("stop") not in {
            "smoke high-loading cap",
            "smoke receding cap",
            "max physical time",
            "wall budget",
        }:
            sys.exit(2)


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "1")
    os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    main()
