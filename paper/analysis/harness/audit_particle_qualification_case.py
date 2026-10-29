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
from contextlib import contextmanager
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


def radial_storage(transport, surface, config, area, pressure_pa=None):
    """Native W/H/U storage from archived cut/dry primitives and surface mass.

    pressure_pa: the film's liquid energy pressure (default None: core.PRESSURE_PA, the former
    constant). A fully dry row's film area is the outer face area of the archived grid this
    function rebuilds (the producer's state geometry), not the caller's 4 pi R^2; at R = 0.885
    mm the two are 2 ulp apart, equal at the qualification radius (2026-10-20 completion build).
    """
    pressure_pa = core.PRESSURE_PA if pressure_pa is None else pressure_pa
    if "geometry" in transport:
        inv, _ = independent.inventory(transport, config)
    else:
        grid = spherical_grid._from_faces(tuple(transport["grid"]["faces"]), transport["grid"]["R"])
        if tuple(transport["grid"]["volumes"]) != grid.volumes:
            raise ValueError("archived geometry differs from native spherical geometry")
        area = grid.areas[-1]
        cells = core.ct._equilibrium_states(
            transport["temperatures_k"], transport["y_hexane"], config.dry
        )
        native = core.ct._inventories(grid, cells)
        inv = [native.total_water_mol, native.total_hexane_mol, native.total_energy_j]
    film = area * surface["liquid_water_kg_m2"]
    return [
        inv[0] * core.wa.M + film,
        inv[1] * core.hx.M,
        inv[2] + film * core.wa.state_Tp(surface["temperature_k"], pressure_pa, "liquid").u_mass,
    ]


def bath_outlet_time_s(bath):
    """The bath's own horizon: its last waypoint's time (1766.5355008060146 s for the
    qualification bath, the former hard-coded outlet)."""
    return float(bath["waypoints"][-1]["t_s"])


def feed_pressure_schedule(feed, bath):
    """The producer's feed pressure and the row pressure law.

    A feed without the producer's kernel-envelope stamp: (None, constant core.PRESSURE_PA),
    the former path. With it: the bath's first pressure, required equal to the stamp, and the
    bath pressure at each row's time.
    """
    stamp = feed.get("kernel_envelope_feed_pressure_pa")
    if stamp is None:
        return None, lambda t_s: core.PRESSURE_PA
    import particle_kernel_envelope as envelope

    first = envelope.bath_pressure_pa(bath, bath["waypoints"][0]["t_s"])
    if float(stamp) != first:
        raise ValueError("producer feed pressure differs from the bath's first pressure")
    return first, lambda t_s: envelope.bath_pressure_pa(bath, t_s)


def feed_row_total_hexane(row, critical_loading):
    """The total hexane loading the producer passed for a feed row.

    At a zero-attached row (the activation) the producer passed exactly critical(T) + 0.0;
    its stored total is the state's derived loading, which may sit an ulp away (below the
    critical loading, so the native initializer refuses it as subcritical). The row is then
    rebuilt from the exact stamped quantities (T and the exact zero), after a codec check
    that the stored total is within 4 ulp of it. Every other row keeps its stored total.
    """
    total = row["hexane_kg_kg_dry"]
    if row["attached_hexane_kg_kg_dry"] != 0.0:
        return total
    critical = critical_loading(row["temperature_k"])
    if abs(critical - total) > 4 * max(math.ulp(critical), math.ulp(total)):
        raise ValueError("primary/derived feed inventory codec discrepancy")
    return critical + row["attached_hexane_kg_kg_dry"]


def presented_feed_boundary(boundary, row):
    """The feed row's boundary as the producer presented it.

    A boundary on which no dew-line projection acted is returned unchanged (the former
    object: the identity on every unprojected case). When the auditor's rebuild projected the
    bath gas (DTDC_KERNEL_ENVELOPE_GAS_DEW_LINE_PROJECTION armed and the bath supersaturated),
    the producer's stamp of the bath it presented at that row (``bath_temperature_k``,
    ``bath_y_hexane``) decides: the projected gas is kept when the stamp is the projected
    composition, the raw bath gas is rebuilt when the stamp is the raw composition (the feed
    stage presents the raw bath), anything else is refused. Returns (boundary, unprojected).
    """
    raw = getattr(boundary, "dew_bath_y_hexane", None)
    if raw is None:
        return boundary, False
    stamped_y = row.get("bath_y_hexane")
    if row.get("bath_temperature_k") != boundary.temperature_k:
        raise ValueError("producer's presented feed bath temperature differs from the rebuild")
    if stamped_y == boundary.y_hexane:
        return boundary, False
    if stamped_y != raw:
        raise ValueError(
            "producer's presented feed bath composition is neither the raw nor the projected gas"
        )
    fields = dict(
        temperature_k=boundary.temperature_k,
        y_hexane=raw,
        label=boundary.label,
        film_boundary=replace(boundary.film_boundary, y_hexane=raw),
        accepted_beta_limit=boundary.accepted_beta_limit,
        actual_surface_path=boundary.actual_surface_path,
        numerical_cap_seam=boundary.numerical_cap_seam,
    )
    if boundary.jacket_wall_temperature_k is None:
        return boundary_driver.study.FiniteFilmStudyBoundary(**fields), True
    return (
        boundary_driver.pd_boundary.KernelEnvelopeBoundary(
            **fields,
            jacket_wall_temperature_k=boundary.jacket_wall_temperature_k,
            jacket_contact_conductance_w_m2_k=boundary.jacket_contact_conductance_w_m2_k,
        ),
        True,
    )


@contextmanager
def declared_retained_water_law(record):
    """Use the producer's stamped law; preserve an already matching installation."""
    import particle_imbibition_bridge as bridge

    try:
        declaration = record["numerical_gate_ruling_settings"]["imbibition_bridge"]
        expected = (declaration["mode"], declaration["saturated_loading_kg_kg_dry"])
    except KeyError as error:
        raise ValueError("producer retained-water law declaration is required") from error
    mode, loading = expected
    if mode not in bridge.MODES:
        raise ValueError("unknown producer retained-water law")
    if mode == "plain-cap":
        if loading is not None:
            raise ValueError("plain-cap producer must declare no saturated loading")
    elif mode == "former-bridge":
        if loading != bridge.FORMER_BRIDGE_W_SAT:
            raise ValueError("former-bridge producer saturated loading differs from its law")
    elif type(loading) is not float or not math.isfinite(loading):
        raise ValueError("producer saturated loading must be a finite binary64")
    elif not bridge.W_SAT_SWEEP[0] <= loading <= bridge.W_SAT_SWEEP[1]:
        raise ValueError("producer saturated loading is outside the ruled range")
    current = bridge.active()
    if current is not None:
        if current != expected:
            raise ValueError("audit/producer retained-water laws differ")
        yield declaration
    else:
        with bridge.installed(mode, loading):
            yield declaration


def audit(
    folder,
    bath_path,
    *,
    accepted_prefix=False,
    live_source_hashes=None,
    return_balance=False,
    until_core_endpoint=False,
):
    """Independently rebuild storage with the exact declared retained-water law."""
    record = json.loads((folder / "case.json").read_text())
    with declared_retained_water_law(record) as declaration:
        result = _audit_with_active_retained_water_law(
            folder,
            bath_path,
            accepted_prefix=accepted_prefix,
            live_source_hashes=live_source_hashes,
            return_balance=return_balance,
            until_core_endpoint=until_core_endpoint,
        )
        report = result[1] if return_balance else result
        report["audit_runtime_bridge"] = {
            "mode": declaration["mode"],
            "saturated_loading_kg_kg_dry": declaration["saturated_loading_kg_kg_dry"],
            "source": "case.json numerical_gate_ruling_settings.imbibition_bridge",
        }
        return result


def _audit_with_active_retained_water_law(
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
    feed_pressure, row_pressure = feed_pressure_schedule(feed, bath)
    march = feed_driver.FeedMarch(
        bath,
        feed["radius_m"],
        record["base_cells"],
        property_parameters=ParticleModelParameters(),
        pressure_pa=feed_pressure,
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
    feed_shells = feed.get("kernel_envelope_feed_cells")
    if feed_shells is not None:
        # FEED_CELLS: the shell-resolved feed's storage rebuilt independently from each row's
        # stamped primitives (shell temperatures, attached n-hexane, exterior water) with the
        # native high-loading caloric per shell; the outer shell's temperature carries the
        # attached liquid and the exterior water.
        shell_grid = core.q.backend.uniform_grid(int(feed_shells), feed["radius_m"])
        shell_wet = core.sphere._material_wet_params(
            march.params, march.state.retained_water_loading, march.state.oil_fraction
        )
        shell_dry = [volume * shell_wet.rho_dm_p for volume in shell_grid.volumes]

    def shell_feed_storage(row):
        temperatures = row["shell_temperatures_k"]
        if len(temperatures) != len(shell_dry) or temperatures[-1] != row["temperature_k"]:
            raise ValueError("shell-resolved feed row does not match its declared shells")
        attached = row["attached_hexane_kg"]
        film = mass * row["external_water_kg_kg_dry"]
        surface = temperatures[-1]
        hexane = math.fsum(
            (
                *(
                    m * core.sphere.critical_loading(t, march.params)
                    for m, t in zip(shell_dry, temperatures)
                ),
                attached,
            )
        )
        energy = math.fsum(
            (
                *(
                    m
                    * core.sphere.combined_high_loading_specific_energy(
                        t, core.sphere.critical_loading(t, march.params), shell_wet
                    )
                    for m, t in zip(shell_dry, temperatures)
                ),
                attached * core.hx.state_Tp(surface, march.pressure, "liquid").u_mass,
                film * core.wa.state_Tp(surface, march.pressure, "liquid").u_mass,
            )
        )
        return [march.state.retained_water_mass_kg + film, hexane, energy]

    def feed_storage(row):
        if feed_shells is not None:
            return shell_feed_storage(row)
        state = core.sphere.initialize_qualified_feed(
            march.params,
            total_hexane_loading=feed_row_total_hexane(
                row, lambda t: core.sphere.critical_loading(t, march.params)
            ),
            temperature_k=row["temperature_k"],
            retained_water_loading=row["retained_water_kg_kg_dry"],
            oil_fraction=march.wet.w_o,
        )
        film = mass * row["external_water_kg_kg_dry"]
        return [
            state.retained_water_mass_kg + film,
            state.total_hexane_mass_kg,
            state.total_internal_energy_j
            + film * core.wa.state_Tp(row["temperature_k"], march.pressure, "liquid").u_mass,
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
    presented_raw_rows = 0
    for old, row in zip(feed["rows"], feed["rows"][1:]):
        march.time = old["t_s"]
        boundary = boundary_driver.bath_boundary(bath, march.radius, row["t_s"])
        if march.time >= march.pd_end:
            # Outside the inherited PD dome the producer presented the bath; its stamp decides.
            boundary, unprojected = presented_feed_boundary(boundary, row)
            presented_raw_rows += unprojected
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
            radial_storage(transport, surface, config, area, row_pressure(transport["time_s"])),
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
    if complete and balance.clock != bath_outlet_time_s(bath):
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
    if presented_raw_rows:
        # Only when a projection acted in the rebuild (absent otherwise: the former report).
        report["feed_rows_rebuilt_with_the_producers_presented_raw_bath"] = presented_raw_rows
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
