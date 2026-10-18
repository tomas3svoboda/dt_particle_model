"""Audit fresh DT storage and exterior-rate histories from the original feed.

Reuses native state functions and the existing primitive-inventory evaluator;
does not consume the producer's cumulative ledgers as conservation evidence.
An accepted prefix audit cannot qualify an incomplete physical history.
"""

# Legacy study imports change cwd; resolve command paths before importing them.
# ruff: noqa: E402
from pathlib import Path

INVOCATION_CWD = Path.cwd()

import argparse
from dataclasses import replace
import hashlib
import json
import math
import pickle

import audit_core2_dt_external_ledger as independent
import experiment_core2_dt_feed as feed_driver
import experiment_core2_dt_birth as radial
import experiment_core2_dt_face as boundary_driver
from particle_solver_configuration import ParticleModelParameters
from dtdc_simulator.core2.particle import grid as spherical_grid

core = radial.core


class OriginalFeedBalance:
    """Keep the original datum and every signed exterior contribution."""

    def __init__(self, initial, scales, tolerance):
        if len(initial) != 3 or len(scales) != 3:
            raise ValueError("three W/H/U inventories and scales are required")
        if not all(math.isfinite(x) for x in (*initial, *scales, tolerance)):
            raise ValueError("balance inputs must be finite")
        if min(scales) <= 0 or tolerance <= 0:
            raise ValueError("strictly positive original scales/tolerance are required")
        self.initial = tuple(initial)
        self.scales = tuple(scales)
        self.tolerance = tolerance
        self.clock = 0.0
        self.streams = [[], [], []]
        self.maximum = [0.0, 0.0, 0.0]
        self.rows = []
        self.last_actual = self.initial
        self.last_signed_defects = (0.0, 0.0, 0.0)

    def accept(self, clock, actual, outward_rates, source):
        if len(actual) != 3 or len(outward_rates) != 3:
            raise ValueError("three W/H/U values are required")
        if not all(math.isfinite(x) for x in (clock, *actual, *outward_rates)):
            raise ValueError("nonfinite primitive storage or exterior rate")
        elapsed = clock - self.clock
        if elapsed <= 0:
            raise ValueError("disconnected or nonpositive physical clock in " + source)
        terms = [elapsed * rate for rate in outward_rates]
        residuals = [
            math.fsum((value, -origin, *history, term))
            for value, origin, history, term in zip(actual, self.initial, self.streams, terms)
        ]
        normalized = [abs(r) / scale for r, scale in zip(residuals, self.scales)]
        if max(normalized) > self.tolerance:
            raise ValueError("original-feed W/H/U conservation refused: " + str(normalized))
        # A failed check cannot move the independent audit origin or its clock.
        for history, term in zip(self.streams, terms):
            history.append(term)
        self.clock = clock
        self.last_actual = tuple(actual)
        self.last_signed_defects = tuple(residuals)
        self.maximum = [max(a, b) for a, b in zip(self.maximum, normalized)]
        self.rows.append({"time_s": clock, "source": source, "normalized_residuals": normalized})


def validate_producer_sources(record, live_source_hashes=None):
    """A live phase-boundary audit requires the exact current entry source set."""
    verdict = record.get("sources_unchanged_during_run")
    if verdict is False:
        raise ValueError("producer sources changed during execution")
    if verdict is None:
        if not live_source_hashes or live_source_hashes != record["source_sha256_at_entry"]:
            raise ValueError("live prefix requires its complete unchanged runtime source set")
        if any(
            hashlib.sha256((radial.REPO_ROOT / p).read_bytes()).hexdigest() != h
            for p, h in live_source_hashes.items()
        ):
            raise ValueError("live source stamp does not match the actual files")


def radial_storage(transport, surface, config, area):
    """Native W/H/U storage from archived cut/dry primitives and surface mass."""
    if "geometry" in transport:
        inv, _ = independent.inventory(transport, config)
    else:
        grid = spherical_grid._from_faces(tuple(transport["grid"]["faces"]), transport["grid"]["R"])
        if tuple(transport["grid"]["volumes"]) != grid.volumes:
            raise ValueError("archived geometry differs from native spherical geometry")
        cells = core.ct._equilibrium_states(
            transport["temperatures_k"], transport["y_hexane"], config.dry
        )
        native = core.ct._inventories(grid, cells)
        inv = [native.total_water_mol, native.total_hexane_mol, native.total_energy_j]
    film = area * surface["liquid_water_kg_m2"]
    return [
        inv[0] * core.wa.M + film,
        inv[1] * core.hx.M,
        inv[2]
        + film * core.wa.state_Tp(surface["temperature_k"], core.PRESSURE_PA, "liquid").u_mass,
    ]


def audit(
    folder,
    bath_path,
    *,
    accepted_prefix=False,
    live_source_hashes=None,
    return_balance=False,
    until_core_endpoint=False,
):
    record = json.loads((folder / "case.json").read_text())
    if record["study"] != "DT":
        raise ValueError("this audit currently supports the DT original-feed law only")
    validate_producer_sources(record, live_source_hashes)
    if not record["fresh_original_feed"] or not record["single_runtime_source_generation"]:
        raise ValueError("fresh original-feed single-generation case required")
    if not record["complete_history"] and not accepted_prefix:
        raise ValueError("incomplete history; explicitly request an accepted-prefix audit")
    if until_core_endpoint and not accepted_prefix:
        raise ValueError("endpoint-only audit must explicitly declare its prefix scope")
    if record["bath_sha256"] != hashlib.sha256(bath_path.read_bytes()).hexdigest():
        raise ValueError("the prescribed bath differs from the producer's exact input")
    native_pins = {
        p: digest for p, digest in record["source_sha256_at_entry"].items() if p.startswith("src/")
    }
    if not native_pins or any(
        hashlib.sha256((radial.REPO_ROOT / p).read_bytes()).hexdigest() != digest
        for p, digest in native_pins.items()
    ):
        raise ValueError("audit native source differs from the original-feed producer")
    birth = json.loads((folder / "birth.json").read_text())
    birth_checkpoint = pickle.loads((folder / "birth.accepted.pkl").read_bytes())
    config = birth_checkpoint["config"]
    feed_path = folder / "feed.json"
    if birth["feed_sha256"] != hashlib.sha256(feed_path.read_bytes()).hexdigest():
        raise ValueError("feed identity differs from the native birth parent")
    feed, bath = json.loads(feed_path.read_text()), json.loads(bath_path.read_text())
    march = feed_driver.FeedMarch(
        bath,
        feed["radius_m"],
        record["base_cells"],
        property_parameters=ParticleModelParameters(),
    )
    if (
        config.wet.wet.oil_heat_capacity_law is not march.wet.oil_heat_capacity_law
        or config.dry.pore.oil_heat_capacity_law is not march.wet.oil_heat_capacity_law
        or config.dry.binary_diffusivity.temperature_law is not ParticleModelParameters().binary_law
        or config.dry.binary_diffusivity.value_m2_s
        != ParticleModelParameters().binary_reference_m2_s
    ):
        raise ValueError("audit/producer property identities differ")
    mass, area = march.dry, march.area

    def feed_storage(row):
        state = core.sphere.initialize_qualified_feed(
            march.params,
            total_hexane_loading=row["hexane_kg_kg_dry"],
            temperature_k=row["temperature_k"],
            retained_water_loading=row["retained_water_kg_kg_dry"],
            oil_fraction=march.wet.w_o,
        )
        film = mass * row["external_water_kg_kg_dry"]
        return [
            state.retained_water_mass_kg + film,
            state.total_hexane_mass_kg,
            state.total_internal_energy_j
            + film * core.wa.state_Tp(row["temperature_k"], core.PRESSURE_PA, "liquid").u_mass,
        ]

    balance = OriginalFeedBalance(
        feed_storage(feed["rows"][0]),
        [
            mass * feed["rows"][0]["retained_water_kg_kg_dry"],
            mass * feed["rows"][0]["hexane_kg_kg_dry"],
            birth_checkpoint["carry"]["import_capacity_scale_j"],
        ],
        config.dry.ledger_tolerance,
    )
    for old, row in zip(feed["rows"], feed["rows"][1:]):
        march.time = old["t_s"]
        boundary = boundary_driver.bath_boundary(bath, march.radius, row["t_s"])
        if march.time < march.pd_end:
            native = core.reduced_film_boundary(
                temperature_k=345.15,
                y_hexane=0.6644429974143522,
                pair=march.pair,
                forcing_timescale_s=0.075,
                label="unchanged prescribed PD dome",
            )
            boundary = replace(
                boundary,
                temperature_k=345.15,
                y_hexane=0.6644429974143522,
                film_boundary=native,
            )
        total, *_ = radial.film_at_beta(
            boundary,
            config,
            march.radius,
            march.ti,
            row["beta"],
            surface_y_hexane=march.yi,
        )
        native, _, _, energy, _, _, active_area = march.interface(
            row["temperature_k"], total, row["t_s"]
        )
        balance.accept(
            row["t_s"],
            feed_storage(row),
            [
                active_area * native.water_molar_flux_mol_m2_s * core.wa.M,
                active_area * native.hexane_molar_flux_mol_m2_s * core.hx.M,
                energy,
            ],
            "native original-feed EOS and independently reconstructed exterior/source rates",
        )

    def accept_radial(transport, surface, kind):
        balance.accept(
            transport["time_s"],
            radial_storage(transport, surface, config, area),
            [
                area * surface[k]
                for k in (
                    "external_water_flux_kg_m2_s",
                    "external_hexane_flux_kg_m2_s",
                    "external_energy_flux_w_m2",
                )
            ],
            kind,
        )

    accept_radial(birth["audit"]["after_transport"], birth["audit"]["surface_water_node"], "birth")
    for line in (folder / "front.history.jsonl").read_text().splitlines():
        row = json.loads(line)
        if row["kind"] not in ("initial", "refusal"):
            accept_radial(row["after_transport"], row["surface"], "front " + row["kind"])
    if (folder / "endpoint.json").exists():
        row = json.loads((folder / "endpoint.json").read_text())
        accept_radial(
            row["after_transport"], row["surface"], "native endpoint and base-face restoration"
        )
    if not until_core_endpoint and (folder / "dry.history.jsonl").exists():
        for line in (folder / "dry.history.jsonl").read_text().splitlines():
            row = json.loads(line)
            if row["kind"] not in ("initial", "refusal"):
                accept_radial(row["after_transport"], row["surface"], "dry " + row["kind"])
    complete = record["complete_history"] and not until_core_endpoint
    if complete and balance.clock != 1766.5355008060146:
        raise ValueError("complete-history flag without the prescribed outlet time")
    report = {
        "passed": True,
        "complete_history": complete,
        "scope": (
            "accepted original-feed prefix through C activation"
            if until_core_endpoint
            else (
                "complete original-feed history"
                if complete
                else (
                    "accepted original-feed prefix at a live phase boundary"
                    if live_source_hashes is not None
                    else "accepted original-feed prefix only; producer stopped"
                )
            )
        ),
        "case_json_sha256": hashlib.sha256((folder / "case.json").read_bytes()).hexdigest(),
        "native_source_pins_checked": len(native_pins),
        "physical_rows_audited": len(balance.rows),
        "last_time_s": balance.clock,
        "maximum_independent_normalized_W_H_U": balance.maximum,
        "original_feed_W_H_U": balance.initial,
        "original_scales": balance.scales,
        "producer_conditioning_counterfactual": record["experimental_per_area_conditioning"],
        "physically_qualifying": False,
        "convergence_established": False,
        "submission_released": False,
        "rows": balance.rows,
        "last_actual_W_H_U": balance.last_actual,
        "last_signed_defects_W_H_U": balance.last_signed_defects,
    }
    return (balance, report) if return_balance else report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--bath", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--accepted-prefix", action="store_true")
    args = parser.parse_args()
    case, bath, out = [(INVOCATION_CWD / p).resolve() for p in (args.case, args.bath, args.out)]
    if out.exists():
        parser.error("refusing to overwrite independent evidence")
    report = audit(case, bath, accepted_prefix=args.accepted_prefix)
    with out.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}), flush=True)


if __name__ == "__main__":
    main()
