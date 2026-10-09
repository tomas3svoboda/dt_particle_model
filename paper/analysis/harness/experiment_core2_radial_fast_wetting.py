"""Conditional radial fast-wetting comparison, starting at accepted core extinction.

D9 supplies fast/no-drain limits, not a measured liquid uptake coefficient.
This experiment approaches fast chemical-potential relaxation in the opened
radial shell while exterior liquid is present. It is NOT the aggregate D9
instantaneous allocation, a Darcy/Washburn law, or an identified rate.
The base native retained-water flux remains; a parallel dissipative path
uses the actual cells' retained-water activities, with (factor - 1) times
the base mobility. Native component and energy assembly,
capacity active sets, film storage and accepted-step ledgers are reused.
No additional water storage or sorption capacity is introduced.
"""

import argparse
from contextlib import contextmanager
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys
import time
import traceback

import experiment_core2_dt_dry as dry
from dtdc_simulator.core2 import through_bed_sorbed_water_law as d9


def inventory_bounded_transfer_fraction(available_flux, demand_flux, dt_s):
    """Normalize the D9 drain in mass units, preserving its exact full branch.

    Converting a full drain back to a rate before normalizing can produce
    1 + one ulp. Dividing the drained mass by the requested mass returns
    exactly one when D9 supplies the complete request, without a clamp.
    """
    if available_flux <= 0:
        return 0.0
    requested_mass = demand_flux * dt_s
    drained_mass = d9.film_drain_kg(
        film_kg=float(available_flux * dt_s),
        capacity_kg=float(requested_mass),
        bracket_end=d9.DRAIN_BRACKET_END_NOMINAL,
    )
    return drained_mass / requested_mass


def add_retained_transfer(native, addition, enthalpy, dissipation):
    """Book one additional liquid flux in component, energy and entropy rows."""
    component = replace(
        native.component,
        retained_water_flux_mol_m2_s=native.component.retained_water_flux_mol_m2_s + addition,
        conserved_water_flux_mol_m2_s=native.component.conserved_water_flux_mol_m2_s + addition,
    )
    energy = replace(
        native.energy,
        retained_water_enthalpy_flux_w_m2=native.energy.retained_water_enthalpy_flux_w_m2
        + addition * enthalpy,
        total_energy_flux_w_m2=native.energy.total_energy_flux_w_m2 + addition * enthalpy,
    )
    return replace(
        native,
        component=component,
        energy=energy,
        retained_water_flux_mol_m2_s=native.retained_water_flux_mol_m2_s + addition,
        retained_water_entropy_production_w_m3_k=native.retained_water_entropy_production_w_m3_k
        + dissipation,
    )


@contextmanager
def liquid_connected_fast_path(
    factor, *, exterior_liquid_present, contact_conductance=None, exclude_wet_front=False
):
    """Add a named parallel path only on an exterior-liquid-connected dry grid."""
    if not math.isfinite(factor) or factor < 1:
        raise ValueError("fast-wetting multiplier must be finite and >= 1")
    if contact_conductance is not None and (
        not math.isfinite(contact_conductance) or contact_conductance <= 0
    ):
        raise ValueError("liquid contact conductance must be finite and positive")
    original = dry.core.cut._dry_flux
    surface_pairs = {}

    def flux(*args, **kwargs):
        is_wet_front = (
            "temperature_gradient_k_m" in kwargs or "composition_gradient_m_inv" in kwargs
        )
        if factor == 1 or not exterior_liquid_present or (exclude_wet_front and is_wet_front):
            return original(*args, **kwargs)
        config = args[8]
        native = original(*args, **kwargs)
        left = dry.cp.evaluate_equilibrium(args[0], config.dry.pressure_pa, args[1], args[6])
        right = dry.cp.evaluate_equilibrium(args[2], config.dry.pressure_pa, args[3], args[7])
        # D9's available capacity is that of the ACTUAL retained inventory.
        # Re-evaluating it at a fictitious common face T can reverse which
        # cell is capped. The added liquid path therefore reads actual W,
        # and removes pure-liquid T dependence by using water activity.
        aw_left = dry.core.sp.water_activity(left.retained_water_loading, args[6].luikov)
        aw_right = dry.core.sp.water_activity(right.retained_water_loading, args[7].luikov)
        gradient = math.log(aw_right / aw_left) / args[4]
        mobility = (factor - 1) * config.dry.retained_water_mobility.value_mol_m_s
        # Only the exterior liquid contact gets a series resistance. Cut
        # assembly marks its wet-core trace with the gradient keywords.
        if contact_conductance is not None and "thermodynamic_force_temperature_k" in kwargs:
            mobility = args[4] / (args[4] / mobility + 1.0 / contact_conductance)
        addition = -mobility * gradient
        donor = left if addition >= 0 else right
        dissipation = dry.core.bg.R * mobility * gradient * gradient
        augmented = add_retained_transfer(
            native, addition, donor.retained_water_enthalpy_j_mol, dissipation
        )
        if "thermodynamic_force_temperature_k" in kwargs:
            surface_pairs[id(augmented)] = (
                augmented,
                native,
                addition,
                donor.retained_water_enthalpy_j_mol,
                dissipation,
            )
        return augmented

    def supply_limited_surface(
        face,
        *,
        old_liquid_loading_kg_m2,
        dt_s,
        external_water_flux_kg_m2_s,
        requested_share=None,
    ):
        if factor == 1 or not exterior_liquid_present:
            return face
        _, native, addition, enthalpy, dissipation = surface_pairs[id(face)]
        available = (
            old_liquid_loading_kg_m2 / dt_s
            + native.component.conserved_water_flux_mol_m2_s * dry.core.wa.M
            - external_water_flux_kg_m2_s
        )
        demand = -addition * dry.core.wa.M
        if demand <= 0:
            return native  # D9's boundary transfer is film -> shell only.
        bounded_share = inventory_bounded_transfer_fraction(available, demand, dt_s)
        share = bounded_share if requested_share is None else requested_share
        if not 0 <= share <= 1:
            raise ValueError("wet-contact transfer fraction outside [0,1]")
        transferred = demand * share
        supply_limited_surface.last_audit = {
            "available_water_flux_kg_m2_s": float(available),
            "requested_liquid_uptake_kg_m2_s": float(demand),
            "transferred_liquid_uptake_kg_m2_s": float(transferred),
            "fraction_of_requested_path": float(share),
            "D9_supply_bounded_fraction": float(bounded_share),
            "old_liquid_loading_kg_m2": float(old_liquid_loading_kg_m2),
            "dt_s": float(dt_s),
            "D9_inventory_bounded_transfer_used": True,
        }
        return add_retained_transfer(native, addition * share, enthalpy, dissipation * share)

    supply_limited_surface.last_audit = None

    dry.core.cut._dry_flux = flux
    try:
        yield supply_limited_surface
    finally:
        dry.core.cut._dry_flux = original


def inventory(state):
    values = dry.ct._inventories(
        state.grid,
        dry.ct._equilibrium_states(state.temperatures_k, state.y_hexane, state.config),
    )
    return (
        values.total_water_mol * dry.core.wa.M,
        values.total_hexane_mol * dry.core.hx.M,
        values.total_energy_j,
    )


def row(state, surface, carry):
    w, h, u = inventory(state)
    volume = math.fsum(state.grid.volumes)
    mass = state.config.pore.rho_dm_p * volume
    film = state.grid.areas[-1] * carry["external_water_loading_kg_m2"]
    cells = dry.ct._equilibrium_states(state.temperatures_k, state.y_hexane, state.config)
    return {
        "time_s": state.time_s,
        "temperature_k": math.fsum(t * v for t, v in zip(state.temperatures_k, state.grid.volumes))
        / volume,
        "particle_water_kg_kg_dry": w / mass,
        "external_water_kg_kg_dry": film / mass,
        "hexane_kg_kg_dry": h / mass,
        "moisture_wet_basis": (w + film) / (mass + w + film + h),
        "retained_water_cells_kg_kg_dry": [c.retained_water_loading for c in cells],
        "y_hexane_cells": list(state.y_hexane),
        "combined_W_H_U": [w + film, h, u + carry["external_water_energy_j"]],
        "surface": asdict(surface),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--factor", type=float, required=True)
    parser.add_argument("--contact-conductance", type=float)
    parser.add_argument("--dt", type=float, default=2.5)
    parser.add_argument("--until", type=float, default=1766.5355008060146)
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--max-halvings", type=int, default=12)
    args = parser.parse_args()
    if args.dt <= 0 or args.steps <= 0 or args.max_halvings < 0:
        raise ValueError("positive timestep and step budget required")
    if args.out.exists() or args.out.with_suffix(".accepted.pkl").exists():
        raise ValueError("refusing to overwrite evidence")
    checkpoint = pickle.loads(args.resume.read_bytes())
    current = checkpoint["current"]
    if hasattr(current, "to_uniform_fully_dry_transport_state"):
        current = current.to_uniform_fully_dry_transport_state()
    if not isinstance(current, dry.ct.FullyDryTransportState):
        raise TypeError("this first comparison requires an accepted fully opened radial grid")
    current = replace(current, config=replace(current.config, nonlinear_residual_tolerance=1e-10))
    config, surface, carry = checkpoint["cut_config"], checkpoint["surface"], checkpoint["carry"]
    bath_path = dry.core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    bath = json.loads(bath_path.read_text())
    root = dry.driver.REPO_ROOT.resolve()
    sources = {Path(__file__).resolve(), bath_path.resolve()}
    sources.update(root.joinpath("src/dtdc_simulator/core2").rglob("*.py"))
    for module in list(sys.modules.values()):
        name = getattr(module, "__file__", None)
        if name and Path(name).suffix == ".py" and Path(name).resolve().is_relative_to(root):
            sources.add(Path(name).resolve())
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    for name, expected in checkpoint["source_sha256"].items():
        path = Path(name)
        if path.is_relative_to(root / "src/dtdc_simulator/core2"):
            if hashes.get(str(path.resolve())) != expected:
                raise ValueError(f"parent native source changed: {path}")
    initial = row(current, surface, carry)
    mass = current.config.pore.rho_dm_p * math.fsum(current.grid.volumes)
    area = current.grid.areas[-1]
    integrated = [0.0, 0.0, 0.0]
    maximum = [0.0, 0.0, 0.0]
    scales = [max(abs(v), 1e-300) for v in initial["combined_W_H_U"]]
    # Energy uses an inventory/capacity scale rather than an arbitrary datum cancellation.
    scales[2] = max(scales[2], mass * current.config.pore.cp_dry_meal)
    record = {
        "experiment": "post-front radial liquid-connected fast-wetting limit approach",
        "physically_qualifying": False,
        "full_feed_to_exit_fast_imbibition_march": False,
        "aggregate_D9_instantaneous_drain_implemented": False,
        "base_transport_coefficient_changed": False,
        "additional_parallel_path_total_to_base_mobility_factor": args.factor,
        "finite_liquid_contact_conductance_mol_m2_s": args.contact_conductance,
        "parallel_path_force": "gradient of log retained-water activity from actual local retained inventory; no explicit thermal-diffusion term",
        "parallel_path_caloric": "native retained-water partial enthalpy at actual donor temperature; zero native excess binding heat",
        "activation": "positive exterior water or supply-limited saturated contact with resolved condensate; fully opened radial shell only",
        "native_capacity_kg_kg_dry": config.dry.pore.luikov.W_cap,
        "capacity_extension_to_D9_x_bf_used": False,
        "parent_checkpoint": str(args.resume),
        "parent_sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
        "source_sha256_at_entry": hashes,
        "initial": initial,
        "attempts": [],
        "history": [initial],
        "dt_requested_s": args.dt,
        "step_halving": f"on rejected step, exact rollback; at most {args.max_halvings} halvings",
        "absent_surface_initial_guess": "bulk bath T and y after depletion or wet contact; otherwise old surface; physical particle and old film storage unchanged",
        "radial_fast_wetting_force_version": "actual_retained_activity_v2",
        "zero_inventory_wet_contact": "water activity unity; extra liquid uptake bounded by resolved condensate/base-pore supply using D9 film_drain_kg; may form new positive exterior film",
        "native_and_combined_ledger_tolerance": 1e-10,
    }
    started = time.perf_counter()
    fatal = None
    for index in range(args.steps):
        if current.time_s >= args.until:
            break
        requested = min(args.dt, args.until - current.time_s)
        accepted = False
        for halving in range(args.max_halvings + 1):
            duration = requested / (2**halving)
            present = bool(carry["external_water_loading_kg_m2"] > 0)
            event = False
            try:
                boundary = dry.boundary_driver.bath_boundary(
                    bath, current.grid.R, current.time_s + duration
                )
                surface_guess = surface
                just_depleted = (
                    len(record["history"]) >= 2
                    and record["history"][-1]["external_water_kg_kg_dry"] == 0
                    and record["history"][-2]["external_water_kg_kg_dry"] > 0
                )
                previous_contact = any(
                    a.get("accepted") and a.get("zero_inventory_wet_contact")
                    for a in record["attempts"][-1:]
                )
                if not present and (just_depleted or previous_contact):
                    # At film depletion the old saturated trace can put every
                    # finite-difference trial on the wrong side of the phase
                    # boundary. The new gas-only node has independent T/y:
                    # initialize those unknowns from the admissible bath.
                    # This surface argument is ONLY a nonlinear initial guess;
                    # old film energy/T remain in carry, and current is exact.
                    surface_guess = replace(
                        surface,
                        temperature_k=boundary.temperature_k,
                        surface_y_hexane=boundary.y_hexane,
                    )
                phases = (
                    ("liquid",)
                    if present
                    else (
                        ("contact", "absent", "liquid")
                        if surface.external_water_flux_kg_m2_s < 0
                        and (just_depleted or previous_contact)
                        else ("absent", "contact", "liquid")
                    )
                )
                result = None
                last_error = None
                limiter_audit = None
                for phase in phases:
                    wet_contact = phase in ("liquid", "contact")
                    guess = surface if wet_contact else surface_guess
                    with (
                        liquid_connected_fast_path(
                            args.factor,
                            exterior_liquid_present=wet_contact,
                            contact_conductance=args.contact_conductance,
                        ) as limiter,
                        dry.driver.bounded_path_proof(3000000, 120.0) as path_audit,
                        dry.driver.path_bounds.monotone_water_path_enclosure([]),
                        dry.study.magnitude_storage_scale(),
                    ):
                        try:
                            result = dry.advance(
                                current,
                                config,
                                guess,
                                carry,
                                duration,
                                boundary,
                                surface_phase_override=None if present else phase,
                                surface_liquid_transfer_adapter=(
                                    limiter if phase == "contact" else None
                                ),
                            )
                            if phase == "contact":
                                contact_audit = limiter.last_audit
                                if result[1].liquid_water_kg_m2 != 0.0:
                                    raise RuntimeError(
                                        "wet contact must resolve exactly zero film storage"
                                    )
                                if contact_audit is None:
                                    raise RuntimeError(
                                        "wet contact has no active inward liquid transfer"
                                    )
                                deficit = abs(
                                    contact_audit["transferred_liquid_uptake_kg_m2_s"]
                                    - contact_audit["available_water_flux_kg_m2_s"]
                                )
                                if deficit * duration * area / scales[0] > 1e-10:
                                    raise RuntimeError("wet-contact resolved supply mismatch")
                        except Exception as root_error:
                            result = None
                            last_error = root_error
                            record["attempts"].append(
                                {
                                    "accepted": False,
                                    "time_s": current.time_s,
                                    "attempt_dt_s": duration,
                                    "surface_phase_candidate": phase,
                                    "error": str(root_error),
                                    "traceback": traceback.format_exc(),
                                    "optimizer_debug": getattr(root_error, "optimizer_debug", None),
                                    "path_audit": dict(path_audit),
                                }
                            )
                            if present:
                                result = dry.advance(
                                    current,
                                    config,
                                    surface,
                                    carry,
                                    duration,
                                    boundary,
                                    zero_film_event=True,
                                    boundary_factory=lambda d: dry.boundary_driver.bath_boundary(
                                        bath, current.grid.R, current.time_s + d
                                    ),
                                )
                                event = True
                                if not 0 < result[0].time_s - current.time_s <= duration:
                                    raise ValueError(
                                        "depletion root lies outside attempted interval"
                                    )
                        if result is not None:
                            limiter_audit = limiter.last_audit
                            break
                if result is None:
                    raise last_error
                after, audit, next_carry, report = result
                report["fast_path_enabled"] = wet_contact and args.factor > 1
                report["fast_path_factor"] = args.factor if wet_contact else 1.0
                report["zero_inventory_wet_contact"] = wet_contact and not present
                report["surface_liquid_supply_limiter"] = limiter_audit
                if wet_contact and not present:
                    report["external_water_phase"] = "wet contact with resolved supply"
                report["native_config_branch_is_coefficient_metadata_only"] = False
                candidate = row(after, audit, next_carry)
                dt_actual = after.time_s - current.time_s
                flows = [
                    audit.external_water_flux_kg_m2_s,
                    audit.external_hexane_flux_kg_m2_s,
                    audit.external_energy_flux_w_m2,
                ]
                next_integrated = [s + dt_actual * area * f for s, f in zip(integrated, flows)]
                defects = [
                    abs(v - i + s) / scale
                    for v, i, s, scale in zip(
                        candidate["combined_W_H_U"],
                        initial["combined_W_H_U"],
                        next_integrated,
                        scales,
                    )
                ]
                if max(defects) > 1e-10:
                    raise RuntimeError("independent cumulative exterior-flux ledger refused")
                integrated = next_integrated
                maximum = [max(a, b) for a, b in zip(maximum, defects)]
                report["independent_cumulative_W_H_U"] = defects
                record["attempts"].append(report)
                record["history"].append(candidate)
                current, surface, carry = after, audit, next_carry
                accepted = True
                args.out.parent.mkdir(parents=True, exist_ok=True)
                args.out.with_suffix(".progress.pkl").write_bytes(
                    pickle.dumps(
                        {
                            **checkpoint,
                            "current": current,
                            "surface": surface,
                            "carry": carry,
                            "source_sha256": hashes,
                            "result_path": str(args.out),
                            "fast_wetting_study_factor": args.factor,
                            "radial_fast_wetting_force_version": "actual_retained_activity_v2",
                        },
                        protocol=5,
                    )
                )
                if index % 5 == 0 or event:
                    args.out.with_suffix(".progress.json").write_text(
                        json.dumps(
                            {
                                **record,
                                "final": candidate,
                                "status": "RUNNING; not final evidence",
                            },
                            indent=2,
                        ),
                        encoding="utf-8",
                    )
                if index % 5 == 0 or event:
                    print(
                        json.dumps(
                            {
                                "step": index + 1,
                                "time_s": current.time_s,
                                "W": candidate["particle_water_kg_kg_dry"],
                                "L": candidate["external_water_kg_kg_dry"],
                                "event": event,
                                "phase": phase,
                            }
                        ),
                        flush=True,
                    )
                break
            except Exception as exc:
                record["attempts"].append(
                    {
                        "accepted": False,
                        "time_s": current.time_s,
                        "attempt_dt_s": duration,
                        "error": str(exc),
                        "traceback": traceback.format_exc(),
                    }
                )
                print(
                    json.dumps(
                        {
                            "rejected_at_s": current.time_s,
                            "attempt_dt_s": duration,
                            "error": str(exc),
                        }
                    ),
                    flush=True,
                )
        if not accepted:
            fatal = record["attempts"][-1]["error"]
            break
    record.update(
        final=row(current, surface, carry),
        requested_window_reached=bool(current.time_s >= args.until),
        bath_exit_reached=bool(current.time_s >= bath["waypoints"][-1]["t_s"]),
        error=fatal,
        wall_s=time.perf_counter() - started,
        independent_maximum_normalized_W_H_U=maximum,
        independent_physical_error_bounds_per_kg_dry_W_H_U=[
            v * s / mass for v, s in zip(maximum, scales)
        ],
        source_sha256_at_exit={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)
        },
    )
    record["sources_unchanged_during_run"] = record["source_sha256_at_exit"] == hashes
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.with_suffix(".accepted.pkl").write_bytes(
        pickle.dumps(
            {
                **checkpoint,
                "current": current,
                "surface": surface,
                "carry": carry,
                "result_path": str(args.out),
                "fast_wetting_study_factor": args.factor,
                "radial_fast_wetting_force_version": "actual_retained_activity_v2",
                "source_sha256": hashes,
            },
            protocol=5,
        )
    )
    args.out.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "bath_exit_reached": record["bath_exit_reached"],
                "final": record["final"],
                "error": fatal,
                "wall_s": record["wall_s"],
            }
        ),
        flush=True,
    )
    if fatal:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
