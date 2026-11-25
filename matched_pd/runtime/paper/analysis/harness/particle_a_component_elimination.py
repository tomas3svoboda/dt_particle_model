"""Equation-equivalent elimination of A's two linear component stock rows."""

import math
import numpy as np
from scipy.optimize import least_squares
import particle_feed_shells as shells


def solve(march, requested_dt, event=False, *, topology="two_liquid"):
    n = march.shells
    has_water = topology == "two_liquid"
    old_t = list(march.temperatures)
    old_storage = march.storage(old_t, march.attached, march.free_water)
    old_shell = [m * march.shell_energy(t) for m, t in zip(march.shell_dry, old_t)]
    old_critical = [march.critical(t) for t in old_t]
    old_surface = march.attached * march.liquid_hexane_u(
        old_t[-1]
    ) + march.free_water * march.liquid_water_u(old_t[-1])
    latest = {}
    calls = 0

    def residual(x):
        nonlocal calls
        calls += 1
        temperatures = list(map(float, x[:n]))
        dt = float(x[n + 1]) if event else requested_dt
        total_flux = float(x[n])
        ti, yi = march.two_liquid_trace if has_water else (float(x[-2]), float(x[-1]))
        saved_trace = march.ti, march.yi
        try:
            march.ti, march.yi = ti, yi
            native, ack, interface_e, out, tb, yb, area = march.interface(
                temperatures[-1], total_flux, march.time + dt
            )
            water_u = march.liquid_water_u(ti)
        finally:
            march.ti, march.yi = saved_trace
        critical = [march.critical(t) for t in temperatures]
        hexane = old_storage[0] - dt * area * native.hexane_molar_flux_mol_m2_s * shells.core.hx.M
        water = (
            march.free_water - dt * area * native.water_molar_flux_mol_m2_s * shells.core.wa.M
            if has_water
            else 0.0
        )
        internal_hexane = math.fsum(m * c for m, c in zip(march.shell_dry, critical))
        attached = 0.0 if event else hexane - internal_hexane
        shell_e = [m * march.shell_energy(t) for m, t in zip(march.shell_dry, temperatures)]
        surface_e = attached * march.liquid_hexane_u(temperatures[-1]) + water * water_u
        qface = [
            g * (temperatures[i] - temperatures[i + 1])
            for i, g in enumerate(march.face_conductance_w_k)
        ]
        moved = [
            m * (c - c0) * march.liquid_hexane_u(t)
            for m, c, c0, t in zip(march.shell_dry, critical, old_critical, temperatures)
        ]
        rows = []
        for i in range(n):
            qin = (qface[i - 1] if i else 0.0) - (qface[i] if i < n - 1 else 0.0)
            terms = [shell_e[i], -old_shell[i], -dt * qin]
            if i < n - 1:
                terms.append(-moved[i])
            else:
                terms.extend([surface_e, -old_surface, dt * out, *moved[:-1]])
            rows.append(math.fsum(terms) / (march.shell_dry[i] * 2e5))
        rows.append(interface_e / (march.area * 1e4))
        if event:
            rows.append((internal_hexane - hexane) / march.dry)
        gas = shells.core.bg.state(ti, march.pressure, 1 - yi, yi)
        h_ratio = (
            gas.fugacity_hexane_pa / shells.core.hx.state_Tp(ti, march.pressure, "liquid").fugacity
        )
        w_ratio = (
            gas.fugacity_water_pa / shells.core.wa.state_Tp(ti, march.pressure, "liquid").fugacity
        )
        if not has_water:
            conductance = march.pair.state.gas_density_kg_m3 * march.pair.binary_mass_transfer_m_s
            rows.extend(
                [
                    native.water_molar_flux_mol_m2_s * shells.core.wa.M / conductance,
                    math.log(h_ratio),
                ]
            )
        latest.update(
            temperatures=temperatures,
            attached=attached,
            water=water,
            energy=math.fsum((*shell_e, surface_e)),
            hexane=internal_hexane + attached,
            dt=dt,
            native=native,
            ack=ack,
            energy_out=out,
            tb=tb,
            yb=yb,
            active_area=area,
            residual=np.asarray(rows),
            ti=ti,
            yi=yi,
            h_ratio=h_ratio,
            w_ratio=w_ratio,
        )
        return latest["residual"]

    x0 = [*old_t, 0.2] + ([requested_dt * 0.5] if event else [])
    bounds = (
        [march.wet.T_min + 1e-4] * n + [-100.0] + ([1e-9] if event else []),
        [march.wet.T_max - 1e-4] * n + [100.0] + ([requested_dt] if event else []),
    )
    if not has_water:
        x0.extend([(march.two_liquid_trace[0] + march.pure_hexane_temperature_k) * 0.5, 0.95])
        bounds[0].extend([march.two_liquid_trace[0], 0.5])
        bounds[1].extend([march.pure_hexane_temperature_k, 1 - 1e-10])
    result = least_squares(
        residual,
        x0,
        bounds=bounds,
        x_scale="jac",
        jac="3-point",
        ftol=1e-14,
        xtol=1e-14,
        gtol=None,
        max_nfev=150,
    )
    residual(result.x)
    if not result.success or np.max(np.abs(latest["residual"])) > 1e-10:
        error = shells.FeedShellsRefused("eliminated radial A residual did not close")
        error.details = dict(
            maximum=float(np.max(np.abs(latest["residual"]))),
            coordinates=result.x.tolist(),
            calls=calls,
        )
        raise error
    if (
        latest["water"] < 0
        or (not event and latest["attached"] < 0)
        or abs(latest["ack"].beta) > 3
        or (not has_water and latest["w_ratio"] > 1 + 1e-10)
    ):
        error = shells.FeedShellsRefused(
            "closed eliminated radial A trial violates a stock, phase or Ackermann bound"
        )
        error.details = {key: latest[key] for key in ("water", "attached", "w_ratio", "ti", "yi")}
        error.details["beta"] = latest["ack"].beta
        raise error
    if event:
        if latest["attached"] != 0.0 or latest["native"].hexane_molar_flux_mol_m2_s <= 0:
            raise shells.FeedShellsRefused("eliminated event does not close exact absent-H balance")
    increment = np.array(
        [
            latest["dt"]
            * latest["active_area"]
            * latest["native"].hexane_molar_flux_mol_m2_s
            * shells.core.hx.M,
            latest["dt"]
            * latest["active_area"]
            * latest["native"].water_molar_flux_mol_m2_s
            * shells.core.wa.M,
            latest["dt"] * latest["energy_out"],
        ]
    )
    cumulative = march.cumulative + increment
    old_trace = march.ti, march.yi
    march.ti, march.yi = latest["ti"], latest["yi"]
    actual = np.asarray(march.storage(latest["temperatures"], latest["attached"], latest["water"]))
    ledger = actual - np.asarray(march.initial) + cumulative
    scales = np.array(
        [
            max(march.initial[0], actual[0], abs(cumulative[0])),
            march.dry * march.state.retained_water_loading + max(actual[1], abs(cumulative[1])),
            march.dry * 2e5,
        ]
    )
    normalized = np.abs(ledger) / scales
    if max(normalized) > 1e-10:
        march.ti, march.yi = old_trace
        raise shells.FeedShellsRefused("eliminated radial A cumulative ledger failed")
    latest["normalized_cumulative_ledger"] = normalized.tolist()
    latest["eliminated_stock_diagnostic"] = dict(residual_calls=calls, nfev=result.nfev)
    march.cumulative = cumulative
    march.temperatures = latest["temperatures"]
    march.attached = latest["attached"]
    march.free_water = latest["water"]
    march.total_energy = float(actual[2])
    march.time += latest["dt"]
    march.rows.append(march.row(latest, ledger))
    return latest
