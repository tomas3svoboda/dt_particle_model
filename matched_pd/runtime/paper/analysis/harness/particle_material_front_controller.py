"""Material-interface continuation using existing native material-face events.

This preserves every historical label on a fixed master grid when history is
nonuniform. It reuses the existing controller's ordinary/corner/core logic;
its native exact-face departure replaces coalescence in that case only.
The emitted source and diff are reviewable. No particle equation is changed.
"""

import difflib
import hashlib
from pathlib import Path
import types

DEPARTURE = """                if getattr(current, "requires_face_departure_chart", False):
                    duration = min(args.maximum_dt, args.until - transport.time_s)
                    committed = current
                    while duration >= args.minimum_dt:
                        path_work = begin_advance()
                        try:
                            step, new_surface, departure_audit = engine.face.advance_departure(
                                current, carry, surface, bath, duration,
                                tangent_residual_bound=2e-10,
                            )
                            combined = engine.driver.combined_ledger(
                                current, config, carry, step, new_surface
                            )
                            if combined["maximum_normalized_residual"] > current.solver_controls.ledger_tolerance:
                                raise RuntimeError("native material-face combined ledger refused")
                            new_carry = {**carry, **combined,
                                "external_water_loading_kg_m2": new_surface.liquid_water_kg_m2}
                            new_current = step.after
                            row = engine.primitive_row(new_current, config, new_carry,
                                new_surface, "native_material_face_departure")
                            elapsed = row["time_s"] - transport.time_s
                            if elapsed <= 0 or new_current.transport.geometry.front.z >= transport.geometry.front.z:
                                raise RuntimeError("native material-face departure did not recede in positive time")
                            area = transport.geometry.master_grid.areas[-1]
                            streams = [s.copy() for s in checkpoint["audit_streams"]]
                            for stream, flux in zip(streams, (
                                new_surface.external_water_flux_kg_m2_s,
                                new_surface.external_hexane_flux_kg_m2_s,
                                new_surface.external_energy_flux_w_m2)):
                                stream.append(elapsed * area * flux)
                            origin = checkpoint["audit_origin"]
                            defects = [math.fsum((a, -b, *s)) for a,b,s in
                                zip(row["combined_W_H_U"], origin["W_H_U"], streams)]
                            normalized = [abs(v)/s for v,s in zip(defects,origin["scales"])]
                            if max(normalized) > current.solver_controls.ledger_tolerance:
                                raise RuntimeError("native material-face independent history ledger refused")
                            row.update(ledger=asdict(step.ledger),combined_conservation=combined,
                                independent_normalized_residuals=normalized,
                                native_material_face=departure_audit,path_work=dict(path_work),
                                material_history_labels_preserved=True)
                            checkpoint.update(
                                front_rate=(transport.geometry.front.z-new_current.transport.geometry.front.z)/elapsed,
                                next_dt=min(args.maximum_dt,elapsed*1.5),audit_streams=streams,
                                positive_core_targets=[],missing_face=None)
                            current,carry,surface=new_current,new_carry,new_surface
                            record["accepted_steps"] += 1
                            record["native_material_departures"] = record.get("native_material_departures",0)+1
                            append(row)
                            record["last_accepted"]={k:v for k,v in row.items() if k not in ("after_transport",)}
                            save()
                            print(json.dumps(dict(kind="native_material_face_departure",
                                time_s=row["time_s"],front_z=current.transport.geometry.front.z)),flush=True)
                            break
                        except Exception as error:
                            if current is not committed:
                                raise RuntimeError("native material-face persistence failed after commit") from error
                            append(dict(kind="refusal",time_s=transport.time_s,dt_s=duration,
                                attempted="native_material_face_departure",error=str(error),exact_rollback=True))
                            record["refusals"]+=1
                            if isinstance(error,StepBudgetExceeded):
                                raise
                            duration*=.5
                    else:
                        raise RuntimeError("native material-face departure exhausted the original minimum dt")
                    continue
"""


def build(out, *, stop_file=None):
    import qualify_sliding_front_history as original

    out = out.resolve()
    path = Path(original.__file__)
    source = path.read_text(encoding="utf-8")
    variant = source.replace(
        'CODEC = "paper1-conservative-sliding-front-v3-approved-chart-control"',
        'CODEC = "paper1-native-material-front-v5"',
    )
    # Exact-face payloads retain their own config/time directly; they are not
    # strict-cut states and must not be wrapped or coerced into one.
    variant = variant.replace("state.transport.", 'getattr(state, "transport", state).')
    needle = "                transport = current.transport\n"
    assert variant.count(needle) == 1
    variant = variant.replace(
        needle, '                transport = getattr(current, "transport", current)\n'
    )
    needle = "                cut_index = transport.geometry.cut_cell_index\n"
    assert variant.count(needle) == 1
    variant = variant.replace(needle, DEPARTURE + needle)
    needle = (
        "                        if event_attempt and not positive_target and not film_event:\n"
    )
    assert variant.count(needle) == 1
    variant = variant.replace(
        needle,
        '                        if event_attempt and not positive_target and not film_event and not record["preserve_material_interfaces"]:\n',
    )
    needle = "                        new_z = new_current.transport.geometry.front.z\n"
    assert variant.count(needle) == 1
    variant = variant.replace(
        needle,
        '                        new_z = getattr(new_current, "transport", new_current).geometry.front.z\n',
    )
    needle = '        "convergence_established": False,\n'
    assert variant.count(needle) == 1
    variant = variant.replace(
        needle,
        needle
        + '        "preserve_material_interfaces": len(set(current.transport.historical_hexane_loadings)) > 1,\n',
    )
    if stop_file is not None:
        needle = "                if transport.time_s >= args.until:\n"
        assert variant.count(needle) == 1
        variant = variant.replace(
            needle,
            f'                if Path({str(stop_file)!r}).exists():\n                    record["stop"]="COOPERATIVE_RESOURCE_STOP"\n                    save()\n                    break\n'
            + needle,
        )
    target = out / "material_front_runtime.py"
    target.write_text(variant, encoding="utf-8")
    (out / "material_front_runtime.diff").write_text(
        "".join(
            difflib.unified_diff(
                source.splitlines(True),
                variant.splitlines(True),
                fromfile=str(path),
                tofile=str(target),
            )
        ),
        encoding="utf-8",
    )
    module = types.ModuleType("particle_material_front_v5")
    module.__file__ = str(target)
    exec(compile(variant, str(target), "exec"), vars(module))
    old_hashes = module.source_hashes

    def hashes():
        return old_hashes() | {
            str(target.relative_to(original.ROOT)): hashlib.sha256(target.read_bytes()).hexdigest(),
            str(Path(__file__).relative_to(original.ROOT)): hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
        }

    module.source_hashes = hashes
    return module
