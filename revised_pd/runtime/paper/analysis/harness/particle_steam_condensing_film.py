"""Kernel-envelope STEAM_CONDENSING_FILM: the true-gas condensing-film branch of the dry stage (armed only).

The owner's ruling of 2026-11-04 (docs/GT_PS2_OWNER_RULING_STEAM_CONDENSATION_AND_CONTINUATION_2026-11-04.md,
the scoped true-gas steam-condensation branch: the actual gas, conservative water/hexane/energy exchange,
the positive-beta domain retained, the existing band on every former path, Paper 1 unchanged, no global
beta extension) and the reviewed candidate of 2026-11-05
(docs/GT_PS2_CONDENSATION_REVIEWED_PARTICLE_AND_VAPOUR_ASSESSMENT_2026-11-05.md; its isolated source
tools/qualify_steam_condensation.py, sha256 6030b210...).  This module ports that candidate's process-
local monkey-patches into a declared regime of the dry dispatcher; record
docs/GT_PS2_DTDC_COLD_START_CONDENSATION_INTEGRATION_2026-11-07.md.  Imported only by
particle_postcore_solver.advance_dry_phases when DTDC_KERNEL_ENVELOPE_STEAM_CONDENSING_FILM is exactly
"1"; nothing on the knob-absent path imports it.

The physics.  A particle colder than the dew point of a hexane-lean superheated steam (the sparge-swept
trays' quench) condenses steam on its surface: its water-saturated interface vapour carries n-hexane
that the gas cannot, so the native Stefan film relation's root for the particle's own n-hexane supply
lies at an inward (condensing) total flux.  The native law's inverse refuses such a root when the
TRIAL beta of its bracket leaves the harness band |beta| <= 3 (an exploratory point, not an accepted
state: the review measured accepted beta -0.198 to +0.040 on both kept arrivals).

What it does, at the start of every dry advance while armed (``trigger``):

* the trigger: a binary dry state (not the pure-water vertex) under a presented gas whose n-hexane lies
  below STEAM_CONDENSING_FILM_GAS_Y_HEXANE and whose water activity on the kernel's saturation law lies
  below 1 (superheated) at both ends of the span, hotter than the stored surface, and EITHER the stored
  surface lies below the gas's water dew point on the kernel's own chart (the gas's n-hexane below the
  lower bound of the gas-only composition interval at the surface temperature: water condenses) OR the
  particle carries an exterior water film (the branch carried to the film's native depletion);
* the search, inside the advance only (every patch restored on exit): the boundary's beta enclosure is
  a finite inward search box derived from the native initial coefficients and composition
  (``search_bound``, at least the original 3); the native inverse is tried first at the caller's limit
  clipped to the original band, and only on its typed beta refusal is the inward root searched on
  [lower, 0] where the native zero-hexane-flux relation bounds it; Newton failures are retried once from
  a one-ulp neighbour seed and, if the last surface-energy row alone fails, polished along a single
  coordinate and re-certified by a fresh native Newton on the FULL original residual, every evaluation
  charged to the original search's and the shared physical-advance allowances; the liquid set of a
  film-free particle is solved from EXACTLY ZERO incoming film (the film created through the native
  mass and energy balances), the former dispatch (native liquid, depletion event, contact and absent
  sets) taken when that refuses;
* the acceptance: every accepted root's beta inside the ORIGINAL harness domain (|beta| <= 3, the
  boundary's original limit) and its reconstructed external-film profile audited
  (particle_numerical_repairs.audit_film_profile), else the set is refused; every accepted beta, the
  search bound and the trigger readings stamped on the advance (native_report key REPORT_KEY); a typed
  ``SteamCondensingFilmRefused`` (a DryPhaseAdvanceError) when no set has a root inside the domain.

No gate, tolerance, allowance, budget, presented gas or accepted beta domain moves; no clamp, no film
dose, no state projection.  physically_qualifying: false.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import fields, make_dataclass
import math

import numpy as np
from scipy.optimize import brentq

import experiment_core2_dt_dry as dry
import experiment_core2_dt_face as face
import particle_conservative_surface as conservative
import particle_kernel_envelope as envelope

core = conservative.core
REGIME = "steam_condensing_film"
REPORT_KEY = "kernel_envelope_steam_condensing_film"
RECORD = envelope.STEAM_CONDENSING_FILM_RECORD
#: the candidate this module ports (tools/qualify_steam_condensation.py at the review)
CANDIDATE_SHA256 = "6030b210caacf10e98fa9a134b5934394fb11c71a9d8b975d737efb5db4fa66e"
#: The positive-film stop rule (labelled measurement variants; the declared rule is the default).
#: "candidate" (declared, the reviewed candidate's rule; see "candidate" below).
#: "depletion_horizon": the native stop at the first nonpositive-film trial of an eliminated
#: liquid advance is kept only when the runner's own curvature-corrected depletion prediction lies
#: within STOP_HORIZON_INTERVALS attempted intervals (a film close to its exact-zero event), and lifted
#: otherwise: a growing or still-thick film's exploratory nonpositive trials are rejected as ordinary
#: inadmissible trials, a film about to deplete keeps the native stop and its event.  Measured
#: alternatives: "candidate" (the reviewed candidate: lifted for every eliminated liquid advance
#: under the bath clamp) spent the shared allowance in a nearly depleted film's Newton (fine P000004,
#: kernel 51.125 s, 1,200 of 1,400 assemblies, with the carried Jacobian) and passed every history
#: once the carried Jacobian is dropped (CARRIED_JACOBIAN); "creation_only"
#: (lifted for the creation from zero only) spent it in a thick film's spurious depletion events (a
#: 2 s job on P000004, kernel 50.3 to 50.5 s, film 1.1e-3 kg/m2, prediction tens of intervals ahead);
#: "unless_film_decreasing" (lifted while the film grows) the same thick decreasing film.
STOP_RULE = "candidate"
STOP_RULES = ("depletion_horizon", "unless_film_decreasing", "creation_only", "candidate")
#: The declared horizon, attempted intervals: between the two measured failures (4.5 and tens).
STOP_HORIZON_INTERVALS = 8.0
#: The half-step-first abandonment threshold inside the regime (a labelled measurement variant; None:
#: the lane's declared envelope.DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD).
HALF_STEP_THRESHOLD = None
#: The first Newton attempt's share of the solve's own evaluation allowance inside the regime (a
#: labelled measurement variant; None: the whole allowance, the candidate's order).  The retry and
#: the surface-energy polish then run on the remainder of the same allowance.
FIRST_NEWTON_NFEV = None
#: The carried dry Jacobian (kernel-envelope JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES) at a regime
#: advance's start (a labelled measurement variant): "drop" starts every regime advance from a fresh
#: finite-difference Jacobian, Broyden updates inside the advance unchanged, as each of the reviewed
#: candidate's one-advance jobs did; "keep" (WAS) seeds it with the previous advance's certificate.
CARRIED_JACOBIAN = "drop"
CARRIED_JACOBIAN_RULES = ("drop", "keep")


def _presented_y(boundary):
    bath = getattr(boundary, "dew_bath_y_hexane", None)
    return float(boundary.y_hexane if bath is None else bath)


def _water_activity(boundary, pressure_pa):
    t = float(boundary.temperature_k)
    if t >= core.wa.TC:
        return 0.0
    return (1.0 - _presented_y(boundary)) * pressure_pa / core.wa.saturation_pressure(t)


def interface_y_lower_bound(surface_temperature_k, pressure_pa, pore):
    """The lowest gas n-hexane fraction the kernel's gas-only chart admits at the surface temperature
    (a water-saturated interface): a gas leaner than this condenses water at that surface."""
    return core.cp.gas_only_composition_interval(
        surface_temperature_k, pressure_pa, pore
    ).y_hexane_bounds[0]


def trigger(before, config, carry, surface, duration, boundary_factory):
    """The trigger readings when the regime engages on this advance, else None."""
    if not isinstance(before, core.ct.FullyDryTransportState):
        return None
    if getattr(before, "regime", None) is not None:
        return None
    start = boundary_factory(before.time_s)
    end = boundary_factory(before.time_s + duration)
    pressure = float(before.config.pressure_pa)
    ys = (_presented_y(start), _presented_y(end))
    if not all(y < envelope.STEAM_CONDENSING_FILM_GAS_Y_HEXANE for y in ys):
        return None
    activities = tuple(_water_activity(b, pressure) for b in (start, end))
    if not all(a < 1.0 for a in activities):
        return None
    t_s = float(surface.temperature_k)
    gas_t = (float(start.temperature_k), float(end.temperature_k))
    if not all(t > t_s for t in gas_t):
        return None
    film = float(carry["external_water_loading_kg_m2"])
    interface_y = interface_y_lower_bound(t_s, pressure, before.config.pore)
    below_dew = all(y < interface_y for y in ys)
    if not (below_dew or film > 0.0):
        return None
    return {
        "regime": REGIME,
        "record": RECORD,
        "surface_temperature_k": t_s,
        "pressure_pa": pressure,
        "film_kg_m2": film,
        "surface_below_gas_dew_point": below_dew,
        "interface_y_lower_bound": interface_y,
        "gas_y_hexane_start_end": list(ys),
        "gas_temperature_k_start_end": list(gas_t),
        "gas_water_activity_start_end": list(activities),
        "threshold_y_hexane": envelope.STEAM_CONDENSING_FILM_GAS_Y_HEXANE,
        "physically_qualifying": False,
    }


def search_bound(before, config, boundary, surface_temperature_k):
    """The finite inward search box (the candidate's, from the native initial coefficients): a search
    device, never an accepted domain."""
    cell = core.ct._equilibrium_states(before.temperatures_k, before.y_hexane, before.config)[-1]
    film = conservative._reduced_film(
        cell,
        config,
        before.grid.R,
        before.grid.R - before.grid.centers[-1],
        boundary,
        surface_temperature_k,
    )
    capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=surface_temperature_k,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=before.config.pressure_pa,
        y_hexane=cell.y_hexane,
        k_wh=before.config.pore.k_wh,
    )
    bulk = core._hexane_mass_fraction(boundary.y_hexane)
    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    return max(
        float(boundary.accepted_beta_limit),
        conductance
        * (abs(math.log(bulk)) + 1.0)
        * max(capacities.water_j_kg_k, capacities.hexane_j_kg_k)
        / film.heat_transfer_w_m2_k,
    )


_BOUNDARY_TYPES: dict = {}


def _validate_condensing_boundary(self):
    core.ct.DirichletPoreBoundary.__post_init__(self)
    if self.film_boundary is None or not math.isfinite(self.accepted_beta_limit):
        raise ValueError("finite condensing-film search envelope required")
    if (self.temperature_k, self.y_hexane) != (
        self.film_boundary.temperature_k,
        self.film_boundary.y_hexane,
    ):
        raise ValueError("actual gas identity must be retained")


def _widened(existing, bound):
    parent = type(existing)
    if parent not in _BOUNDARY_TYPES:
        _BOUNDARY_TYPES[parent] = make_dataclass(
            "Condensing" + parent.__name__,
            [],
            bases=(parent,),
            frozen=True,
            namespace={"__post_init__": _validate_condensing_boundary},
        )
    parameters = {item.name: getattr(existing, item.name) for item in fields(existing)}
    parameters["accepted_beta_limit"] = bound
    return _BOUNDARY_TYPES[parent](**parameters)


def _closure_cells(function):
    return dict(
        zip(
            function.__code__.co_freevars,
            (cell.cell_contents for cell in function.__closure__ or ()),
        )
    )


#: STEAM_CONDENSING_FILM_FLOOR_EXIT: the floor band (a multiple of gtol) and the patience (a
#: multiple of the unknowns: evaluations without a halving of the best largest residual)
FLOOR_EXIT_FACTOR = 10.0
FLOOR_EXIT_PATIENCE = 2


class _FloorStagnation(Exception):
    """The regime's first Newton stagnated at the rounding floor (not a refusal; handled)."""


def _forward_jacobian(function, x, residual, lower, upper):
    """A forward-difference Jacobian at x (the polish's search slopes; its certificate is the
    native Newton's own), each probe inside the open bounds, else the backward probe."""
    columns = []
    for column in range(x.size):
        h = math.sqrt(np.finfo(float).eps) * max(1.0, abs(x[column]))
        trial = x.copy()
        trial[column] += h
        if not lower[column] < trial[column] < upper[column]:
            trial[column] = x[column] - h
        step = trial[column] - x[column]
        columns.append((np.asarray(function(trial), dtype=float) - residual) / step)
    return np.column_stack(columns)


def _floor_exit_first(function, initial, options, original_newton, polish, stamp):
    """STEAM_CONDENSING_FILM_FLOOR_EXIT: the first Newton, handed to the polish at the floor.

    Returns the native result when it converges or fails before the floor; at the floor the
    polished result when the polish closes, else an unconverged result at the best iterate (its
    Jacobian the forward one) for the former retry.  Every evaluation is counted in nfev.
    """
    from scipy.optimize import OptimizeResult

    gtol = float(options.get("gtol", 1e-11))
    x0 = np.asarray(initial, dtype=float)
    patience = FLOOR_EXIT_PATIENCE * x0.size
    track = {"calls": 0, "best": math.inf, "x": None, "r": None, "since": 0, "mark": math.inf}
    counts = stamp.setdefault("floor_exit", {"exits": 0, "polished": 0, "handed_back": 0})

    def tracked(values):
        residual = np.asarray(function(values), dtype=float)
        track["calls"] += 1
        merit = float(np.max(np.abs(residual)))
        if math.isfinite(merit) and merit < track["best"]:
            track["best"], track["x"], track["r"] = merit, np.array(values, dtype=float), residual
        if track["best"] <= FLOOR_EXIT_FACTOR * gtol:
            if track["best"] <= 0.5 * track["mark"]:
                track["mark"], track["since"] = track["best"], 0
            else:
                track["since"] += 1
            if track["since"] >= patience:
                raise _FloorStagnation()
        return residual

    try:
        result = original_newton(tracked, initial, **options)
        return result
    except _FloorStagnation:
        pass
    counts["exits"] += 1
    lower, upper = map(np.asarray, options["bounds"])
    x, residual = track["x"], track["r"]
    jac = _forward_jacobian(function, x, residual, lower, upper)
    used = track["calls"] + x.size
    failed = OptimizeResult(
        x=x,
        fun=residual,
        jac=jac,
        success=False,
        nfev=used,
        message="steam condensing film: the first Newton stagnated at the rounding floor",
    )
    polished, used = polish(function, failed, options, used)
    polished.nfev = used
    if polished.success:
        counts["polished"] += 1
    else:
        counts["handed_back"] += 1
    return polished


#: STEAM_CONDENSING_FILM_SEEDED_START: the state-continuity tolerance (absolute, in the solve's
#: shifted coordinates: K for the cell temperatures, the composition fraction, kg/m2 for the film)
SEEDED_START_CONTINUITY = 1e-12


def _seeded_newton(function, initial, options, duration, stamp, solve):
    """STEAM_CONDENSING_FILM_SEEDED_START: the regime solve from the extrapolated previous solve.

    The memory is the job's dry lever scope (one per kernel job, so a job's values never depend on
    the jobs a worker ran before); without a scope the former start.  The state unknowns are the
    cell temperatures and compositions and the last unknown; the third block holds per-advance
    increments and is extrapolated, not compared.  Both extra residual calls are charged to the
    solve's allowance (the shared physical-advance allowance charges itself in the residual).
    """
    import particle_numerical_repairs as repairs

    levers = repairs._ACTIVE_JACOBIAN_LEVERS
    x0 = np.asarray(initial, dtype=float)
    origins = _closure_cells(function).get("origins")
    if levers is None or origins is None or (x0.size - 1) % 3 or duration is None:
        return solve(function, initial, **options)
    origins = np.asarray(origins, dtype=float)
    memory = levers.__dict__.setdefault("steam_condensing_film_seed", {})
    counts = stamp.setdefault(
        "seeded_start", {"continued": 0, "seeded": 0, "kept_former": 0, "refused": 0}
    )
    n = (x0.size - 1) // 3
    state = np.r_[np.arange(2 * n), x0.size - 1]
    previous = memory.get(x0.size)
    start, spent = x0, 0
    if previous is not None and np.all(
        np.abs((x0 + origins - previous[0])[state]) <= SEEDED_START_CONTINUITY
    ):
        counts["continued"] += 1
        seed = x0 + previous[1] * (duration / previous[2])
        lower, upper = map(np.asarray, options["bounds"])
        if np.all(lower < seed) and np.all(seed < upper) and options["max_nfev"] > 8:
            try:
                spent += 1
                at_seed = float(np.max(np.abs(function(seed))))
                spent += 1
                at_former = float(np.max(np.abs(function(x0))))
            except (ValueError, RuntimeError) as error:
                if "budget" in str(error).lower() or "source" in str(error).lower():
                    raise
                counts["refused"] += 1
            else:
                if math.isfinite(at_seed) and at_seed < at_former:
                    start = seed
                    counts["seeded"] += 1
                else:
                    counts["kept_former"] += 1
    trial = dict(options)
    trial["max_nfev"] = options["max_nfev"] - spent
    result = solve(function, start, **trial)
    result.nfev += spent
    if result.success:
        memory[x0.size] = (np.asarray(result.x) + origins, np.asarray(result.x) - x0, duration)
    return result


@contextmanager
def scope(bound, original_limit, stamp, lift_stop=False):
    """The regime's search patches for one dry advance; every patched callable restored on exit."""
    numerics = dry.study.numerics
    if face.study.numerics is not numerics:
        raise RuntimeError("the face and dry modules must share one numerics module")
    saved = {
        "bath_boundary": face.bath_boundary,
        "inverse": conservative.total_flux_for_hexane_target,
        "newton": numerics.domain_aware_newton,
        "advance": dry.advance,
    }
    original_boundary = saved["bath_boundary"]
    original_inverse = saved["inverse"]
    original_newton = saved["newton"]
    original_advance = saved["advance"]
    inverse_calls = stamp.setdefault("inverse_fallbacks", [])
    # Kernel-envelope STEAM_CONDENSING_FILM_SEEDED_START (armed only): the advance's duration
    # for the seeded start (set by advance below, read by newton)
    seeded = envelope.armed(envelope.STEAM_CONDENSING_FILM_SEEDED_START)
    floor_exit = envelope.armed(envelope.STEAM_CONDENSING_FILM_FLOOR_EXIT)
    request_memory = envelope.armed(envelope.STEAM_CONDENSING_FILM_REQUEST_MEMORY)
    current = {"duration": None}
    polish_trials = stamp.setdefault("polish_trials", [])
    accepted = stamp.setdefault("accepted_roots", [])

    def boundary(*values, **options):
        return _widened(original_boundary(*values, **options), bound)

    def inverse(film, target, interface, bulk, cp_water, cp_hexane, limit):
        try:
            return original_inverse(
                film, target, interface, bulk, cp_water, cp_hexane, min(limit, original_limit)
            )
        except ValueError as error:
            if (
                limit <= original_limit
                or "beta" not in str(error)
                or not (0.0 < bulk < interface < 1.0)
            ):
                raise
        coefficient = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
        zero_hexane_flux = coefficient * math.log(bulk / interface)
        count = 0

        def evaluate(total):
            nonlocal count
            count += 1
            return core.efr.component_fluxes_for_supplied_total_mass_flux(
                film,
                total_mass_flux_kg_m2_s=total,
                interface_hexane_mass_fraction=interface,
                bulk_hexane_mass_fraction=bulk,
            )

        if target == 0.0:
            total = zero_hexane_flux
        else:
            zero = evaluate(0.0)
            beta_zero = core.efr.ackermann_heat_transfer(
                film,
                fluxes=zero,
                water_film_heat_capacity_j_kg_k=cp_water,
                hexane_film_heat_capacity_j_kg_k=cp_hexane,
                enforce_production_domain=False,
            ).beta
            lower = min(
                (-limit - beta_zero) * film.heat_transfer_w_m2_k / capacity
                for capacity in (cp_water, cp_hexane)
            )
            total = brentq(
                lambda value: evaluate(value).hexane_kg_m2_s - target,
                lower,
                0.0,
                xtol=math.ulp(0.0),
                rtol=4 * np.finfo(float).eps,
            )
        fluxes = evaluate(total)
        ackermann = core.efr.ackermann_heat_transfer(
            film,
            fluxes=fluxes,
            water_film_heat_capacity_j_kg_k=cp_water,
            hexane_film_heat_capacity_j_kg_k=cp_hexane,
            enforce_production_domain=False,
        )
        if ackermann.beta < -bound:
            raise ValueError("candidate left its finite condensation search envelope")
        if len(inverse_calls) < 64:
            inverse_calls.append({"total_flux": total, "beta": ackermann.beta})
        stamp["inverse_fallback_count"] = stamp.get("inverse_fallback_count", 0) + 1
        stamp["inverse_fallback_min_beta"] = min(
            stamp.get("inverse_fallback_min_beta", math.inf), ackermann.beta
        )
        return total, fluxes, ackermann.beta, count

    def polish_surface_energy(function, failed, options, used):
        closure = _closure_cells(function)
        cells = (len(failed.x) - 1) // 3
        if (
            "origins" not in closure
            or "acceptance" not in options
            or len(failed.x) != 3 * cells + 1
        ):
            return failed, used
        origins = closure["origins"]
        gates = options["acceptance"].gates(failed.x, failed.jac)[0]
        ratios = np.abs(failed.fun) / gates
        if not (np.all(ratios[:-1] <= 1.0) and ratios[-1] > 1.0):
            return failed, used
        lower, upper = map(np.asarray, options["bounds"])
        driver = _closure_cells(closure["evaluate"]) if "evaluate" in closure else {}
        counter = driver.get("physical_advance_assembly_budget")
        reserve = len(failed.x) + 2
        for column in list(range(cells, 3 * cells)) + [3 * cells] + list(range(cells)):
            slope = failed.jac[-1, column]
            if not math.isfinite(slope) or slope == 0.0:
                continue
            target = failed.x[column] - failed.fun[-1] / slope
            physical_target = target + origins[column]
            for value in (
                target,
                np.nextafter(physical_target, -math.inf) - origins[column],
                np.nextafter(physical_target, math.inf) - origins[column],
            ):
                if options["max_nfev"] - used <= reserve or (
                    counter is not None and counter["limit"] - counter["used"] <= reserve
                ):
                    return failed, used
                if value == failed.x[column] or not lower[column] < value < upper[column]:
                    continue
                point = failed.x.copy()
                point[column] = value
                used += 1
                try:
                    residual = function(point)
                except (ValueError, RuntimeError) as error:
                    if "budget" in str(error).lower() or "source" in str(error).lower():
                        raise
                    continue
                if np.all(np.isfinite(residual)) and np.all(np.abs(residual) <= gates):
                    retry = dict(options)
                    retry["max_nfev"] = options["max_nfev"] - used
                    certified = original_newton(function, point, **retry)
                    used += certified.nfev
                    return certified, used
        return failed, used

    def newton(function, initial, **options):
        if seeded:
            return _seeded_newton(function, initial, options, current["duration"], stamp, _newton)
        return _newton(function, initial, **options)

    def _newton(function, initial, **options):
        first = dict(options)
        capped = (
            FIRST_NEWTON_NFEV is not None
            and "max_nfev" in options
            and FIRST_NEWTON_NFEV < options["max_nfev"]
        )
        if not capped and floor_exit:
            # STEAM_CONDENSING_FILM_FLOOR_EXIT: the first Newton ends at the rounding floor
            result = _floor_exit_first(
                function, initial, first, original_newton, polish_surface_energy, stamp
            )
            if result.success:
                return result
        elif not capped:
            result = original_newton(function, initial, **first)
        else:
            # the first attempt on its declared share; the best evaluated iterate (smallest largest
            # scaled residual) seeds the retry on the remainder of the same allowance
            first["max_nfev"] = FIRST_NEWTON_NFEV
            best = [math.inf, None]

            def tracked(values):
                residual = function(values)
                merit = float(np.max(np.abs(residual)))
                if math.isfinite(merit) and merit < best[0]:
                    best[0], best[1] = merit, np.array(values, dtype=float)
                return residual

            # restarted shares: each attempt on at most FIRST_NEWTON_NFEV evaluations, restarted from
            # the best iterate while that improves, on the one solve allowance; a stagnation return
            # goes to the surface-energy polish
            used, seed, restarts = 0, initial, []
            while True:
                share = min(FIRST_NEWTON_NFEV, options["max_nfev"] - used)
                if share <= 2:
                    raise RuntimeError("domain-aware Newton residual-call budget exhausted")
                before_best = best[0]
                attempt = dict(options)
                attempt["max_nfev"] = share
                try:
                    result = original_newton(tracked, seed, **attempt)
                except RuntimeError as error:
                    if "residual-call budget exhausted" not in str(error):
                        raise
                    used += share
                    restarts.append({"share": share, "best_merit": best[0]})
                    if best[1] is None or not best[0] < before_best:
                        raise
                    seed = best[1]
                    continue
                used += result.nfev
                trial = {
                    "restarts": restarts,
                    "last_used": result.nfev,
                    "success": bool(result.success),
                }
                if not result.success:
                    result, used = polish_surface_energy(function, result, options, used)
                    trial["surface_energy_polish_success"] = bool(result.success)
                result.nfev = used
                if restarts and len(polish_trials) < 32:
                    polish_trials.append(trial)
                return result
        remaining = options.get("max_nfev", 0) - result.nfev
        if result.success or remaining <= 2:
            return result
        closure = _closure_cells(function)
        retry_seed = result.x
        probe_count = 0
        if "origins" in closure and "acceptance" in options:
            origins = closure["origins"]
            physical = result.x + origins
            column = (len(initial) - 1) // 3 - 1
            gates = options["acceptance"].gates(result.x, result.jac)[0]
            best = float(np.max(np.abs(result.fun) / gates))
            for direction in (math.inf, -math.inf):
                trial = result.x.copy()
                trial[column] = np.nextafter(physical[column], direction) - origins[column]
                probe_count += 1
                try:
                    residual = function(trial)
                    merit = float(np.max(np.abs(residual) / gates))
                    if merit < best:
                        retry_seed, best = trial, merit
                    if merit <= 1.0:
                        break
                except ValueError, RuntimeError:
                    pass
        next_options = dict(options)
        next_options["max_nfev"] = remaining - probe_count
        polished = original_newton(function, retry_seed, **next_options)
        trial = {
            "first_used": result.nfev,
            "remaining": remaining,
            "probe_count": probe_count,
            "polish_used": polished.nfev,
            "success": bool(polished.success),
        }
        used = polished.nfev + result.nfev + probe_count
        if not polished.success:
            polished, used = polish_surface_energy(function, polished, options, used)
            trial["surface_energy_polish_success"] = bool(polished.success)
        polished.nfev = used
        if len(polish_trials) < 32:
            polish_trials.append(trial)
        return polished

    def advance(*values, **options):
        # The positive-film stop at the first nonpositive-film trial is lifted only for the film's
        # creation from zero (the route's own liquid call passes it); a depleting film keeps the
        # native stop and its exact-zero event.  WAS (the candidate): lifted for every eliminated
        # liquid advance under a gas at the bath clamp, which on the fine P000004 history spent the
        # shared allowance in the depleting film's Newton (1,200 of 1,400 assemblies at 51.125 s).
        if STOP_RULE not in STOP_RULES:
            raise ValueError(f"undeclared stop rule {STOP_RULE!r}")
        if STOP_RULE == "candidate":
            condensing_trial = values[5].y_hexane <= envelope.HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP
            if condensing_trial and options.get("conservative_liquid_elimination"):
                options["stop_at_nonpositive_liquid"] = False
        elif (
            STOP_RULE in ("depletion_horizon", "unless_film_decreasing")
            and lift_stop
            and options.get("conservative_liquid_elimination")
            and options.get("stop_at_nonpositive_liquid")
        ):
            options["stop_at_nonpositive_liquid"] = False
            stamp["stop_lifted"] = stamp.get("stop_lifted", 0) + 1
        current["duration"] = float(values[4])
        returned = original_advance(*values, **options)
        accepted_duration = float(values[4])
        after, surface = returned[:2]
        if abs(surface.beta) > original_limit:
            raise ValueError(
                "steam condensing film: accepted root outside the original physical beta domain"
            )
        outer = core.ct._equilibrium_states(after.temperatures_k, after.y_hexane, after.config)[-1]
        resolved_film = conservative._reduced_film(
            outer,
            values[1],
            after.grid.R,
            after.grid.R - after.grid.centers[-1],
            values[5],
            surface.temperature_k,
        )
        numerics.audit_film_profile(
            resolved_film,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(surface.surface_y_hexane),
            bulk_hexane_mass_fraction=core._hexane_mass_fraction(values[5].y_hexane),
            fluxes=core.efr.ComponentMassFluxes(
                water_kg_m2_s=surface.external_water_flux_kg_m2_s,
                hexane_kg_m2_s=surface.external_hexane_flux_kg_m2_s,
            ),
            closure_residual_budget_kg_m2_s=1e-12
            * resolved_film.film_density_kg_m3
            * resolved_film.binary_mass_transfer_m_s,
        )
        accepted.append(
            {
                "clock_s": after.time_s,
                "beta": surface.beta,
                "surface_temperature_k": surface.temperature_k,
                "film_kg_m2": surface.liquid_water_kg_m2,
                "water_flux_kg_m2_s": surface.external_water_flux_kg_m2_s,
                "hexane_flux_kg_m2_s": surface.external_hexane_flux_kg_m2_s,
                "profile_audited": True,
            }
        )
        if request_memory:
            # STEAM_CONDENSING_FILM_REQUEST_MEMORY: the accepted duration for the next request
            _remember_accepted_duration(accepted_duration)
        return returned

    if CARRIED_JACOBIAN not in CARRIED_JACOBIAN_RULES:
        raise ValueError(f"undeclared carried-Jacobian rule {CARRIED_JACOBIAN!r}")
    import particle_numerical_repairs as repairs

    levers = repairs._ACTIVE_JACOBIAN_LEVERS
    if CARRIED_JACOBIAN == "drop" and levers is not None and levers.carried_across_advances:
        stamp["carried_jacobians_dropped"] = len(levers.carried)
        levers.carried = {}
    face.bath_boundary = boundary
    conservative.total_flux_for_hexane_target = inverse
    numerics.domain_aware_newton = newton
    dry.advance = advance
    try:
        yield
    finally:
        face.bath_boundary = saved["bath_boundary"]
        conservative.total_flux_for_hexane_target = saved["inverse"]
        numerics.domain_aware_newton = saved["newton"]
        dry.advance = saved["advance"]


def route(
    before,
    config,
    carry,
    surface,
    duration,
    boundary_factory,
    *,
    previous_phase,
    dispatch,
    advance_type,
    error_type,
    options,
):
    """The regime's advance when the trigger holds, else None.

    ``dispatch(**overrides)`` is the former dispatcher with the regime marked inner; ``options`` the
    caller's dispatcher options (assembly_budget, search_tolerance, ...)."""
    readings = trigger(before, config, carry, surface, duration, boundary_factory)
    if readings is None:
        return None
    import experiment_core2_radial_fast_wetting as liquid
    from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded

    end_boundary = boundary_factory(before.time_s + duration)
    original_limit = float(end_boundary.accepted_beta_limit)
    bound = search_bound(before, config, end_boundary, float(surface.temperature_k))
    stamp = {
        **readings,
        "original_beta_limit": original_limit,
        "finite_negative_search_bound": -bound,
        "search_bound_is_not_a_domain": True,
        "candidate_sha256": CANDIDATE_SHA256,
    }
    rejected = []
    result = None
    predicted = options.get("predicted_depletion_duration_s")
    decreasing = predicted is not None and math.isfinite(predicted) and predicted > 0.0
    if STOP_RULE == "depletion_horizon":
        keep_stop = decreasing and predicted <= STOP_HORIZON_INTERVALS * duration
    else:
        keep_stop = decreasing
    stamp["stop_rule"] = STOP_RULE
    stamp["film_decreasing"] = decreasing
    stamp["predicted_depletion_duration_s"] = predicted
    stamp["native_stop_kept"] = keep_stop
    with scope(bound, original_limit, stamp, lift_stop=not keep_stop):
        # the lane worker's factory reads face.bath_boundary at call time: inside the scope it
        # returns the widened search boundary (checked, never assumed)
        probe = boundary_factory(before.time_s + duration)
        if probe.accepted_beta_limit != bound or (probe.temperature_k, probe.y_hexane) != (
            end_boundary.temperature_k,
            end_boundary.y_hexane,
        ):
            raise RuntimeError("the condensing-film scope did not reach the advance's boundary")
        film_free = carry["external_water_loading_kg_m2"] == 0.0 and previous_phase != "liquid"
        creation = None
        if film_free and readings["surface_below_gas_dew_point"]:
            boundary = probe
            with liquid.liquid_connected_fast_path(
                options.get("liquid_mobility_factor", envelope.LIQUID_MOBILITY_TOTAL_FACTOR),
                exterior_liquid_present=True,
                contact_conductance=options.get("contact_conductance")
                or envelope.surface_contact_conductance(),
            ) as limiter:
                try:
                    after, audit, next_carry, report = dry.advance(
                        before,
                        config,
                        surface,
                        carry,
                        duration,
                        boundary,
                        surface_phase_override="liquid",
                        conservative_liquid_elimination=True,
                        stop_at_nonpositive_liquid=False,
                        search_tolerance=options.get("search_tolerance"),
                        physical_advance_assembly_budget=options["assembly_budget"],
                        global_component_mass_error_budgets_kg=options.get(
                            "global_component_mass_error_budgets_kg"
                        ),
                        **(
                            {}
                            if options.get("target_pressure_pa") is None
                            else {"target_pressure_pa": options["target_pressure_pa"]}
                        ),
                    )
                    native_duration = report["dt_s"]
                    if not (
                        math.isfinite(native_duration)
                        and 0 < native_duration <= duration
                        and after.time_s == before.time_s + native_duration
                    ):
                        raise error_type(
                            "dry/event duration is outside its declared attempted interval",
                            before,
                            rejected,
                        )
                    next_carry = {
                        **next_carry,
                        "external_water_mass_kg": after.grid.areas[-1] * audit.liquid_water_kg_m2,
                    }
                    creation = {
                        "outcome": "ok",
                        "old_film_kg_m2": 0.0,
                        "new_film_kg_m2": audit.liquid_water_kg_m2,
                    }
                    result = advance_type(
                        after, audit, next_carry, report, "liquid", False, limiter.last_audit, ()
                    )
                except StepBudgetExceeded:
                    raise
                except error_type:
                    raise
                except Exception as error:
                    creation = {"outcome": "refused", "error": str(error)[:300]}
                    rejected.append(
                        {
                            "phase": "liquid from zero film (steam condensing film)",
                            "error": str(error),
                            "optimizer_debug": getattr(error, "optimizer_debug", None),
                        }
                    )
        stamp["film_creation_from_zero"] = creation
        if result is None:
            try:
                result = dispatch()
            except error_type as error:
                stamp["refused"] = True
                refused = error_type(
                    "steam condensing film: no accepted root inside the original beta domain ("
                    + str(error)
                    + ")",
                    before,
                    tuple(rejected) + tuple(error.rejected_trials),
                )
                refused.kernel_envelope_steam_condensing_film = stamp
                raise refused from error
    stamp["accepted_phase"] = result.phase
    stamp["accepted_beta"] = getattr(result.surface, "beta", None)
    stamp["rejected_before_acceptance"] = len(rejected)
    report = {**result.native_report, REPORT_KEY: stamp}
    return advance_type(
        result.after,
        result.surface,
        result.carry,
        report,
        result.phase,
        result.depletion_event,
        result.liquid_supply_audit,
        tuple(rejected) + tuple(result.rejected_trials),
        result.regime_events,
    )


def _lever_memory():
    import particle_numerical_repairs as repairs

    levers = repairs._ACTIVE_JACOBIAN_LEVERS
    return None if levers is None else levers.__dict__.setdefault("steam_condensing_film_dt", {})


def _remember_accepted_duration(duration):
    memory = _lever_memory()
    if memory is not None:
        memory["last"] = float(duration)


def _last_accepted_duration():
    memory = _lever_memory()
    return None if memory is None else memory.get("last")


def _film_free_steam_cap(before, carry, surface, requested, boundary_factory):
    """STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP: (request, stamp) when it applies, else None."""
    if not isinstance(before, core.ct.FullyDryTransportState):
        return None
    if getattr(before, "regime", None) is not None:
        return None
    if float(carry["external_water_loading_kg_m2"]) > 0.0:
        return None
    start = boundary_factory(before.time_s)
    end = boundary_factory(before.time_s + requested)
    pressure = float(before.config.pressure_pa)
    ys = (_presented_y(start), _presented_y(end))
    if not all(y < envelope.STEAM_CONDENSING_FILM_GAS_Y_HEXANE for y in ys):
        return None
    if not all(_water_activity(b, pressure) < 1.0 for b in (start, end)):
        return None
    t_s = float(surface.temperature_k)
    if not all(float(b.temperature_k) > t_s for b in (start, end)):
        return None
    cap = envelope.STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_S
    if requested <= cap:
        return None
    return cap, {
        "requested_s": requested,
        "capped_s": cap,
        "capped": True,
        "film_free_request_cap": True,
        "half_step_first_threshold": half_step_threshold(),
        "record": envelope.CONVERGENCE_ROOT_CAUSES_RECORD,
        "physically_qualifying": False,
    }


def request_cap(before, config, carry, surface, requested, boundary_factory):
    """The regime's declared dry-advance time cap (envelope.STEAM_CONDENSING_FILM_MAXIMUM_DT_S) when
    the trigger holds over the requested advance: (request, stamp), the stamp None when it does not.
    """
    if trigger(before, config, carry, surface, requested, boundary_factory) is None:
        if envelope.armed(envelope.STEAM_CONDENSING_FILM_FILM_FREE_REQUEST_CAP):
            capped = _film_free_steam_cap(before, carry, surface, requested, boundary_factory)
            if capped is not None:
                return capped
        return requested, None
    cap = envelope.STEAM_CONDENSING_FILM_MAXIMUM_DT_S
    remembered = None
    if envelope.armed(envelope.STEAM_CONDENSING_FILM_REQUEST_MEMORY):
        # STEAM_CONDENSING_FILM_REQUEST_MEMORY: at most twice the last accepted duration
        last = _last_accepted_duration()
        if last is not None and 2.0 * last < cap:
            remembered = 2.0 * last
    limit = cap if remembered is None else remembered
    return min(requested, limit), {
        "requested_s": requested,
        "capped_s": min(requested, limit),
        "capped": requested > limit,
        "request_memory_s": remembered,
        "half_step_first_threshold": half_step_threshold(),
        "record": RECORD,
        "physically_qualifying": False,
    }


def half_step_threshold():
    """The half-step-first abandonment threshold of an advance the regime takes (the labelled
    measurement variant HALF_STEP_THRESHOLD; None: the lane's declared threshold)."""
    if HALF_STEP_THRESHOLD is None:
        return envelope.DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD
    return HALF_STEP_THRESHOLD
