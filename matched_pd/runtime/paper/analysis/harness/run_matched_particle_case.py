"""Fresh corrected prescribed-bath particle case with frozen source closure.

Native birth, exact material-face motion, core projection and C operators.
One cold feed, one runtime, no reused trajectory and no adopted convergence.
"""

# Harness imports follow the explicit repository search path.
# ruff: noqa: E402
import argparse
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import time
import traceback
import types

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "paper/analysis/harness")]
import particle_matched_history_audit as audit
from particle_matched_feed import FreshMatchedFeed
import particle_material_front_controller as material_front_controller
import particle_feed_shells as shells
import particle_original_feed_balance as prefix
import run_particle_qualification_case as producer
from run_pd_prescribed_particle_case import CooperativeStop


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--stop-file", type=Path, required=True)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--main-shells", type=int, required=True)
    parser.add_argument("--time-level", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument(
        "--geometry", choices=("unconfined_sphere", "equal_volume_cell"), required=True
    )
    parser.add_argument("--bath", type=Path, required=True)
    parser.add_argument(
        "--compiled-eos",
        action="store_true",
        help="optional result-neutral evaluator; 13 host bit-identity/replay checks passed",
    )
    args = parser.parse_args()
    args.out, args.stop_file, args.bath = (
        args.out.resolve(),
        args.stop_file.resolve(),
        args.bath.resolve(),
    )
    if args.out.exists() or args.stop_file.exists():
        raise ValueError("use a new output folder and an unrequested stop marker")
    if args.main_shells < 2 or args.cells % args.main_shells:
        raise ValueError("main-A shell count must divide the B/C master-grid count")
    args.out.mkdir(parents=True)
    runtime = material_front_controller.build(args.out, stop_file=args.stop_file)
    feed = FreshMatchedFeed(
        args.out, args.stop_file, geometry=args.geometry, main_shells=args.main_shells
    )
    os.environ["DTDC_KERNEL_ENVELOPE_PRESSURE_SCHEDULE"] = "1"
    os.environ["DTDC_KERNEL_ENVELOPE_FRONT_RELEASE_SOLVE_MEMORY"] = "1"
    os.environ["DTDC_KERNEL_ENVELOPE_CORE_TARGET_SEED"] = "history_power_law"
    if args.compiled_eos:
        os.environ["DTDC_KERNEL_ENVELOPE_IAPWS_COMPILED"] = "1"
    elif producer.envelope.armed(producer.envelope.IAPWS_COMPILED):
        raise ValueError("compiled evaluator requires the explicit --compiled-eos option")

    extension_files = [
        Path(__file__),
        Path(audit.__file__),
        Path(material_front_controller.__file__),
        args.out / "material_front_runtime.py",
        ROOT / "paper/analysis/harness/particle_matched_feed.py",
        ROOT / "paper/analysis/harness/particle_pd_diffusion_boundary.py",
        ROOT / "paper/analysis/harness/particle_a_component_elimination.py",
    ]

    def extension_hashes():
        return {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in extension_files
        }

    entry = extension_hashes()
    old_hashes, old_front_hashes = producer.hashes, runtime.source_hashes
    old_feed, old_front = producer.feed_driver.run, producer.front.run
    old_remap, old_prefix = shells.remap_to_master_grid, prefix.prepare_prefix
    producer.hashes = lambda: old_hashes() | extension_hashes()
    runtime.source_hashes = lambda: old_front_hashes() | extension_hashes()
    manifest = dict(
        status="running",
        fresh_original_feed=True,
        single_runtime_source_generation=True,
        geometry=args.geometry,
        B_C_cells=args.cells,
        main_A_shells=args.main_shells,
        PD_thermal_cells=1,
        time_level=args.time_level,
        extension_sources_at_entry=entry,
        bit_identical_evaluator=args.compiled_eos,
        bath_sha256=hashlib.sha256(args.bath.read_bytes()).hexdigest(),
        scope="candidate evidence; independent audit and qualification required",
    )
    producer.write_json(args.out / "production.json", manifest)
    started, last_check = time.monotonic(), 0.0

    def check_stop(*unused):
        nonlocal last_check
        now = time.monotonic()
        if now - last_check < 0.1:
            return
        last_check = now
        if args.stop_file.exists() or now - started > 8 * 3600:
            raise CooperativeStop("cooperative resource or eight-hour wall stop")

    def prepare(folder, study, bath_path, live):
        if study != "DT" or live != producer.hashes():
            raise ValueError("fresh producer source identity changed before original-feed audit")
        result = json.loads((folder / "case.json").read_text())
        if not result["fresh_original_feed"] or not result["single_runtime_source_generation"]:
            raise ValueError("fresh production requires one cold feed and runtime source set")
        return audit.audit_prefix(
            folder,
            [folder / "front.history.jsonl"],
            endpoint=json.loads((folder / "endpoint.json").read_text()),
        )

    tool = next(i for i in (3, 4) if sys.monitoring.get_tool(i) is None)
    sys.monitoring.use_tool_id(tool, "fresh-matched-particle-resource-stop")
    sys.monitoring.register_callback(tool, sys.monitoring.events.PY_START, check_stop)
    for module in list(sys.modules.values()):
        if isinstance(module, types.ModuleType):
            for value in list(vars(module).values()):
                if inspect.isfunction(value) and "advance" in value.__name__:
                    file = Path(value.__code__.co_filename)
                    if file.is_relative_to(ROOT / "src") or file.is_relative_to(
                        ROOT / "paper/analysis/harness"
                    ):
                        sys.monitoring.set_local_events(
                            tool, value.__code__, sys.monitoring.events.PY_START
                        )
    for function in (producer.write_json, producer.HistoryTimeControl.dry_request):
        sys.monitoring.set_local_events(tool, function.__code__, sys.monitoring.events.PY_START)
    try:
        producer.feed_driver.run, producer.front.run = feed.run, runtime.run
        shells.remap_to_master_grid = lambda *values, **options: old_remap(
            *values, **(options | {"labels": "per_shell"})
        )
        prefix.prepare_prefix = prepare
        sys.argv = [
            "run_particle_qualification_case.py",
            "--study",
            "DT",
            "--cells",
            str(args.cells),
            "--time-level",
            str(args.time_level),
            "--bath",
            str(args.bath),
            "--out",
            str(args.out / "case"),
            "--per-area-qualification",
            "--dry-local-tolerance",
            "1e-10",
            "--conservative-components",
            "--event-first-film-depletion",
            "--eliminate-event-absent-components",
            "--curvature-corrected-depletion-prediction",
            "--memoize-proof-points",
            "--dry-time-schedule",
            "uniform",
            "--conservative-core-endpoint",
        ]
        producer.main()
        result = json.loads((args.out / "case/case.json").read_text())
        manifest.update(
            {
                name: result.get(name)
                for name in ("status", "complete_history", "error", "C_accepted_steps")
            }
        )
        if result["complete_history"]:
            _, report = audit.audit_prefix(
                args.out / "case",
                [args.out / "case/front.history.jsonl"],
                endpoint=json.loads((args.out / "case/endpoint.json").read_text()),
                dry_histories=[args.out / "case/dry.history.jsonl"],
            )
            producer.write_json(args.out / "independent_whole_history_audit.json", report)
            manifest["independent_whole_history_passed"] = True
    except CooperativeStop as error:
        manifest.update(status="COOPERATIVE_RESOURCE_STOP", error=str(error))
    except Exception as error:
        manifest.update(status="REFUSED", error=str(error), traceback=traceback.format_exc())
    finally:
        sys.monitoring.free_tool_id(tool)
        producer.hashes, runtime.source_hashes = old_hashes, old_front_hashes
        producer.feed_driver.run, producer.front.run = old_feed, old_front
        shells.remap_to_master_grid, prefix.prepare_prefix = old_remap, old_prefix
        manifest.update(
            wall_seconds=time.monotonic() - started,
            extension_sources_at_exit=extension_hashes(),
            extension_sources_unchanged=entry == extension_hashes(),
        )
        producer.write_json(args.out / "production.json", manifest)
    print(
        json.dumps(
            {
                key: value
                for key, value in manifest.items()
                if not key.startswith("extension_sources")
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
