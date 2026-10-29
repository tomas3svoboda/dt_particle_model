"""Radially resolved high-loading feed stage (regime A) behind DTDC_KERNEL_ENVELOPE_FEED_CELLS.

A numerical refinement study of the DT feed stage (record
docs/GT_PS2_PAPER1_FEED_STAGE_REFINEMENT_2026-10-24.md; formulation audit finding A10). With the
knob absent nothing here is called and the feed stage is the native ONE-cell thermal body of
experiment_core2_dt_feed.FeedMarch, bit for bit.

What the one-cell feed is. FeedMarch carries the particle before its radial activation as the
native high-loading state of the gate-1f seam (sphere.HighLoadingState, PHY-018): one material
temperature T, the live internal n-hexane X_c(T) and the excess as attached liquid outside,
with the exterior water at the body temperature. Conduction enters only through the interface
energy row, as the native one-cell Fourier conductance A k / (R - r_c) of uniform_grid(1, R)
between the body and the two-liquid interface temperature.

What this module changes, and nothing else. The same state on N uniform shells of the native
grid uniform_grid(N, R) (N = 1 is the one-cell body again, to rounding):

* every shell i carries its own temperature T_i, the live internal loading X_c(T_i) and the
  native wet-core energy m_i e_w(T_i, X_c(T_i)) (the same caloric function the one-cell body's
  combined_high_loading_specific_energy evaluates with zero attached excess);
* the attached n-hexane and the exterior water are surface inventories at the outer shell's
  temperature T_N (the one-cell body's "exterior water shares body temperature");
* neighbouring shells exchange heat by the native Fourier stencil k A_face / (c_{i+1} - c_i) on
  the face-midpoint centres; the outer shell meets the interface through the same one-cell form
  A k / (R - c_N) at T_N (the interface row and its film and source law unchanged);
* the declared PD source (85 W per kg dry) is deposited per kg dry, i.e. in proportion to each
  shell's dry mass;
* the liquid n-hexane that a shell expels or takes up as X_c(T_i) moves with T_i passes to or
  from the attached inventory at the shell's own liquid internal energy u_l(T_i): the same
  Span-Wagner liquid state on both sides, so the transfer conserves energy exactly and carries
  no latent impulse, as the one-cell caloric declares.

Backward Euler in every row, the same least-squares closure, residual gate (1e-10), bounds,
Ackermann band, depletion event and whole-particle cumulative ledger (1e-10) as FeedMarch.

Hand-over. The radial activation (zero attached n-hexane) stores X_c(T_i) per shell
(wet_core.activate's loading, per shell) with the shell temperatures; remap_to_master_grid then
carries that profile onto the birth's master grid conservatively. The kernel's front stage slides
its mesh only across identical adjacent history labels, so the hand-over gives every master cell
one historical loading, the dry-mass-weighted mean of the shells' X_c(T_i) (the whole-particle
n-hexane exact, its radial placement moved by at most the shells' label spread, 2e-4 kg per kg
at four shells), and each master cell the temperature whose native wet-cell energy at that label
equals the summed energy of the shell pieces it overlaps (Brent). The per-shell labels are kept
as labels="per_shell" for the record's measurement (the front refused them). The importer's own
native energy check verifies the whole hand-over. physically_qualifying is false.
"""

from __future__ import annotations

from dataclasses import asdict
import math

import numpy as np
from scipy.optimize import brentq, least_squares

import experiment_core2_dt_feed as feed_driver
from dtdc_simulator.core2.particle import front as native_front
from dtdc_simulator.core2.particle import grid as spherical_grid
from dtdc_simulator.core2.particle import sphere as native_sphere
from dtdc_simulator.core2.particle import wet_core as native_wet_core

core = feed_driver.core

RECORD = "docs/GT_PS2_PAPER1_FEED_STAGE_REFINEMENT_2026-10-24.md"
SOURCE_W_PER_KG_DRY = 85.0  # the FeedMarch.interface declared PD source, read back below


class FeedShellsRefused(RuntimeError):
    """A multi-shell feed step, activation or hand-over that does not close (never repaired)."""


class MultiShellFeedMarch(feed_driver.FeedMarch):
    """FeedMarch on N uniform shells; see the module docstring."""

    def __init__(self, bath, radius, cells, *, shells, property_parameters=None, pressure_pa=None):
        if isinstance(shells, bool) or not isinstance(shells, int) or shells < 1:
            raise FeedShellsRefused(f"shell count must be a positive integer, not {shells!r}")
        super().__init__(
            bath,
            radius,
            cells,
            property_parameters=property_parameters,
            **({} if pressure_pa is None else {"pressure_pa": pressure_pa}),
        )
        lumped = self.state
        self.shells = shells
        self.wet_material = native_sphere._material_wet_params(
            self.params, lumped.retained_water_loading, lumped.oil_fraction
        )
        grid = core.q.backend.uniform_grid(shells, radius)
        self.shell_grid = grid
        rho = self.wet_material.rho_dm_p
        self.shell_dry = [volume * rho for volume in grid.volumes]
        self.shell_dry_total = math.fsum(self.shell_dry)
        self.internal_conductance_w_m2_k = core.CONDUCTIVITY / (radius - grid.centers[-1])
        self.model_view.transfer.solid_to_interface_ua_w_k = (
            self.area * self.internal_conductance_w_m2_k
        )
        self.face_conductance_w_k = [
            core.CONDUCTIVITY * grid.areas[i + 1] / (grid.centers[i + 1] - grid.centers[i])
            for i in range(shells - 1)
        ]
        temperature = lumped.temperature_k
        self.temperatures = [temperature] * shells
        total_hexane = lumped.total_hexane_mass_kg
        internal = math.fsum(m * self.critical(temperature) for m in self.shell_dry)
        self.attached = total_hexane - internal
        if self.attached < 0.0:
            raise FeedShellsRefused("the feed is subcritical on the shell grid")
        self.free_water = 0.0
        self.total_energy = self.storage(self.temperatures, self.attached, 0.0)[2]
        self.initial = (total_hexane, 0.0, self.total_energy)
        self.cumulative = np.zeros(3)
        self.rows = [self.row(None, None)]

    # -- native state functions on the shells ---------------------------------------------------
    def critical(self, temperature):
        return native_sphere.critical_loading(temperature, self.params)

    def shell_energy(self, temperature):
        """J per kg dry of a shell at the live X_c(T) with zero attached excess."""
        return native_sphere.combined_high_loading_specific_energy(
            temperature, self.critical(temperature), self.wet_material
        )

    def liquid_hexane_u(self, temperature):
        return core.hx.state_Tp(temperature, self.pressure, "liquid").u_mass

    def liquid_water_u(self, temperature):
        return core.wa.state_Tp(temperature, self.pressure, "liquid").u_mass

    def storage(self, temperatures, attached_kg, free_water_kg):
        """(hexane kg, exterior water kg, energy J) of the shells plus the surface inventories."""
        hexane = math.fsum(
            (*(m * self.critical(t) for m, t in zip(self.shell_dry, temperatures)), attached_kg)
        )
        surface = temperatures[-1]
        energy = math.fsum(
            (
                *(m * self.shell_energy(t) for m, t in zip(self.shell_dry, temperatures)),
                attached_kg * self.liquid_hexane_u(surface),
                free_water_kg * self.liquid_water_u(surface),
            )
        )
        return hexane, free_water_kg, energy

    @property
    def mass_mean_temperature_k(self):
        return math.fsum(m * t for m, t in zip(self.shell_dry, self.temperatures)) / (
            self.shell_dry_total
        )

    # -- one backward-Euler step -----------------------------------------------------------------
    def solve(self, requested_dt, event=False):
        n = self.shells
        old_t = list(self.temperatures)
        old_attached = self.attached
        old_w = self.free_water / self.dry
        old_hexane = math.fsum(
            (*(m * self.critical(t) for m, t in zip(self.shell_dry, old_t)), old_attached)
        )
        old_shell_e = [m * self.shell_energy(t) for m, t in zip(self.shell_dry, old_t)]
        old_surface_e = math.fsum(
            (
                old_attached * self.liquid_hexane_u(old_t[-1]),
                self.free_water * self.liquid_water_u(old_t[-1]),
            )
        )
        old_critical = [self.critical(t) for t in old_t]
        in_pd = self.time < self.pd_end
        history = {}

        def evaluate(x):
            temperatures = [float(v) for v in x[:n]]
            attached = 0.0 if event else float(x[n]) * self.dry
            free_w = float(x[n] if event else x[n + 1])
            total_flux = float(x[n + 1] if event else x[n + 2])
            dt = float(x[n + 2]) if event else requested_dt
            surface = temperatures[-1]
            critical = [self.critical(t) for t in temperatures]
            shell_e = [m * self.shell_energy(t) for m, t in zip(self.shell_dry, temperatures)]
            water_kg = self.dry * free_w
            surface_e = math.fsum(
                (
                    attached * self.liquid_hexane_u(surface),
                    water_kg * self.liquid_water_u(surface),
                )
            )
            native, ack, interface_e, net_out, tb, yb, active_area = self.interface(
                surface, total_flux, self.time + dt
            )
            source = SOURCE_W_PER_KG_DRY * self.dry if in_pd else 0.0
            boundary_out = net_out + source  # the enthalpy outflow minus the gas heat
            # Conduction between shells (implicit), positive outward.
            q_face = [
                g * (temperatures[i] - temperatures[i + 1])
                for i, g in enumerate(self.face_conductance_w_k)
            ]
            # Liquid n-hexane moved into shell i from the attached inventory (negative: out).
            moved = [m * (c - c0) for m, c, c0 in zip(self.shell_dry, critical, old_critical)]
            moved_e = [dm * self.liquid_hexane_u(t) for dm, t in zip(moved, temperatures)]
            energy_rows = []
            for i in range(n):
                q_in = (q_face[i - 1] if i > 0 else 0.0) - (q_face[i] if i < n - 1 else 0.0)
                share = source * self.shell_dry[i] / self.shell_dry_total
                terms = [shell_e[i], -old_shell_e[i], -dt * q_in, -dt * share]
                if i < n - 1:
                    terms.append(-moved_e[i])
                else:
                    terms.extend((surface_e, -old_surface_e, dt * boundary_out, *moved_e[:-1]))
                energy_rows.append(math.fsum(terms) / (self.shell_dry[i] * 2e5))
            hexane = math.fsum((*(m * c for m, c in zip(self.shell_dry, critical)), attached))
            residual = np.array(
                [
                    (
                        hexane
                        - old_hexane
                        + dt * active_area * native.hexane_molar_flux_mol_m2_s * core.hx.M
                    )
                    / self.dry,
                    free_w
                    - old_w
                    + dt * active_area / self.dry * native.water_molar_flux_mol_m2_s * core.wa.M,
                    *energy_rows,
                    interface_e / (self.area * 1e4),
                ]
            )
            history.update(
                temperatures=temperatures,
                attached=attached,
                water=water_kg,
                hexane=hexane,
                energy=math.fsum((*shell_e, surface_e)),
                dt=dt,
                native=native,
                ack=ack,
                energy_out=net_out,
                residual=residual,
                tb=tb,
                yb=yb,
                active_area=active_area,
            )
            return residual

        lo_t, hi_t = self.wet.T_min + 1e-4, self.wet.T_max - 1e-4
        if event:
            x0 = [*old_t, old_w + 1e-4, 0.2, requested_dt * 0.5]
            lower = [lo_t] * n + [0.0, -100.0, 1e-9]
            upper = [hi_t] * n + [1.0, 100.0, requested_dt]
        else:
            x0 = [*old_t, old_attached / self.dry, old_w + 1e-4, 0.2]
            lower = [lo_t] * n + [0.0, 0.0, -100.0]
            upper = [hi_t] * n + [1.0, 1.0, 100.0]
        solution = least_squares(
            evaluate,
            x0,
            bounds=(lower, upper),
            x_scale="jac",
            xtol=1e-14,
            ftol=1e-14,
            gtol=None,
            max_nfev=150,
            jac="3-point",
        )
        residual = evaluate(solution.x)
        if not solution.success or max(abs(residual)) > 1e-10:
            error = RuntimeError(
                f"native high-loading {'event' if event else 'step'} not closed: "
                f"{max(abs(residual)):.6g}"
            )
            error.details = {
                "residual": residual.tolist(),
                "coordinates": solution.x.tolist(),
                "nfev": solution.nfev,
            }
            raise error
        if abs(history["ack"].beta) > 3.0:
            raise RuntimeError(
                f"closed trial outside declared signed Ackermann study band: {history['ack'].beta}"
            )
        if event and history["native"].hexane_molar_flux_mol_m2_s <= 0:
            raise RuntimeError("attached-film depletion event lacks outward hexane flux")
        if (
            self.free_water == 0
            and history["water"] > 0
            and history["native"].water_molar_flux_mol_m2_s >= 0
        ):
            raise RuntimeError("water appearance lacks inward component flux")
        dt = history["dt"]
        increments = np.array(
            [
                dt
                * history["active_area"]
                * history["native"].hexane_molar_flux_mol_m2_s
                * core.hx.M,
                dt
                * history["active_area"]
                * history["native"].water_molar_flux_mol_m2_s
                * core.wa.M,
                dt * history["energy_out"],
            ]
        )
        trial_cumulative = self.cumulative + increments
        actual = np.array(
            self.storage(history["temperatures"], history["attached"], history["water"])
        )
        ledger = actual - np.array(self.initial) + trial_cumulative
        ledger_scales = np.array(
            [
                max(self.initial[0], actual[0], abs(trial_cumulative[0])),
                self.dry * self.state.retained_water_loading
                + max(actual[1], abs(trial_cumulative[1])),
                self.dry * 2e5,
            ]
        )
        normalized = np.abs(ledger) / ledger_scales
        if np.max(normalized) > 1e-10:
            raise RuntimeError(
                f"whole particle/film cumulative ledger refuses before commit: {normalized.tolist()}"
            )
        history["normalized_cumulative_ledger"] = normalized.tolist()
        self.cumulative = trial_cumulative
        self.temperatures = history["temperatures"]
        self.attached = history["attached"]
        self.free_water = history["water"]
        self.total_energy = float(actual[2])
        self.time += dt
        self.rows.append(self.row(history, ledger))
        return history

    def row(self, step, ledger):
        if not hasattr(self, "temperatures"):
            # FeedMarch.__init__'s own initial row, replaced once the shells exist.
            return super().row(step, ledger)
        surface = self.temperatures[-1]
        hexane, _, _ = self.storage(self.temperatures, self.attached, self.free_water)
        return {
            "t_s": self.time,
            # The outer shell's temperature: the interface's particle temperature and the
            # exterior water's (the one-cell body's single temperature at N = 1).
            "temperature_k": surface,
            "hexane_kg_kg_dry": hexane / self.dry,
            "attached_hexane_kg_kg_dry": self.attached / self.dry,
            "retained_water_kg_kg_dry": self.state.retained_water_loading,
            "external_water_kg_kg_dry": self.free_water / self.dry,
            "total_particle_plus_film_energy_j": self.total_energy,
            "beta": step["ack"].beta if step else None,
            "cumulative_component_energy_residuals": (
                ledger.tolist() if ledger is not None else None
            ),
            "normalized_cumulative_ledger": (
                step.get("normalized_cumulative_ledger") if step else None
            ),
            "bath_temperature_k": step["tb"] if step else None,
            "bath_y_hexane": step["yb"] if step else None,
            "shell_temperatures_k": list(self.temperatures),
            "mass_mean_temperature_k": self.mass_mean_temperature_k,
            "attached_hexane_kg": self.attached,
        }

    # -- radial activation on the shells ---------------------------------------------------------
    def activate(self):
        """The z = 1 radial state on the shell grid (X_c(T_i) stored once per shell)."""
        if self.attached != 0.0:
            raise FeedShellsRefused("radial activation requires exactly zero attached n-hexane")
        grid = self.shell_grid
        wet = self.wet_material
        loadings = [native_wet_core.activate(grid, t, wet).loadings[0] for t in self.temperatures]
        temperatures = list(self.temperatures)
        n = grid.n
        M_dm = sum(volume * wet.rho_dm_p for volume in grid.volumes)
        M_h = native_wet_core.total_hexane(grid, loadings, wet)
        U = native_wet_core.total_energy(grid, temperatures, loadings, wet)
        before_h, _, before_total_e = self.storage(self.temperatures, 0.0, self.free_water)
        film_e = self.free_water * self.liquid_water_u(temperatures[-1])
        before_e = before_total_e - film_e
        scale_e = native_sphere._energy_scale(self.shell_dry_total, wet, temperatures[-1])
        residuals = {
            "dry_meal_residual_kg": M_dm - self.shell_dry_total,
            "hexane_residual_kg": M_h - before_h,
            "energy_residual_j": U - before_e,
        }
        ledger = {
            **residuals,
            "max_normalized_residual": max(
                abs(residuals["dry_meal_residual_kg"]) / self.shell_dry_total,
                abs(residuals["hexane_residual_kg"]) / before_h,
                abs(residuals["energy_residual_j"]) / scale_e,
            ),
        }
        radial = {
            "grid": asdict(grid),
            "coordinate": asdict(native_front.FrontCoordinate(grid.R, 1.0)),
            "wet_inventory": {"loadings": loadings},
            "temperatures_k": temperatures,
            "oil_fraction_labels": [wet.w_o] * n,
            "retained_water_loadings": [wet.X_water] * n,
            "dry_meal_mass_kg": M_dm,
            "retained_water_mass_kg": M_dm * wet.X_water,
            "internal_hexane_mass_kg": M_h,
            "attached_hexane_mass_kg": 0.0,
            "total_internal_energy_j": U,
            # The outer shell's temperature (the exterior water's), stamped as such.
            "temperature_k": temperatures[-1],
            "temperature_meaning": "outer shell (exterior water and interface); "
            "temperatures_k carries the shell profile",
            "mass_mean_temperature_k": self.mass_mean_temperature_k,
        }
        return {
            "radial": radial,
            "ledger": ledger,
            "external_water_mass_kg_carried_separately": self.free_water,
            "external_water_energy_j_carried_separately": film_e,
        }


# -- the hand-over onto the birth's master grid ---------------------------------------------------
def _overlap_volume(a, b, c, d):
    lo, hi = max(a, c), min(b, d)
    if hi <= lo:
        return 0.0
    return spherical_grid.spherical_volume_difference(lo, hi)


def remap_to_master_grid(radial, grid, model, water, oil, evaluate_cell, *, labels="uniform"):
    """Per-cell temperatures and historical loadings of the master grid from the shell profile.

    Returns (temperatures, historical loadings, report).

    labels="uniform" (the hand-over): every master cell carries ONE historical n-hexane loading,
    the dry-mass-weighted mean of the shells' X_c(T_i), so the wet core keeps the uniform
    activation label the kernel's sliding front requires ("sliding requires identical adjacent
    material/history labels"); the whole-particle n-hexane is conserved exactly, its radial
    placement moves by at most max |X_c(T_i) - mean|. Each master cell's temperature solves its
    own energy: the native wet-cell energy at the mean label equals the summed native energy of
    the shell pieces it overlaps (Brent), so the energy is conserved cell by cell.

    labels="per_shell" (measured and refused by the front stage, record section 5): a cell
    inside one shell takes that shell's T and X_c(T_i) exactly; a straddling cell the
    volume-weighted loading and the Brent temperature of its pieces' summed energy.
    """
    if labels not in ("uniform", "per_shell"):
        raise FeedShellsRefused(f"unknown hand-over label rule {labels!r}")
    faces = radial["grid"]["faces"]
    shell_t = radial["temperatures_k"]
    shell_x = radial["wet_inventory"]["loadings"]
    shell_v = radial["grid"]["volumes"]
    shell_density = [
        evaluate_cell(t, water, x, oil, model).energy_density_j_m3 for t, x in zip(shell_t, shell_x)
    ]
    mean_label = math.fsum(v * x for v, x in zip(shell_v, shell_x)) / math.fsum(shell_v)
    temperatures, loadings, inverted = [], [], []
    for j in range(grid.n):
        a, b = grid.faces[j], grid.faces[j + 1]
        pieces = [(i, _overlap_volume(a, b, faces[i], faces[i + 1])) for i in range(len(shell_t))]
        pieces = [(i, v) for i, v in pieces if v > 0.0]
        if labels == "per_shell" and len(pieces) == 1:
            i = pieces[0][0]
            temperatures.append(shell_t[i])
            loadings.append(shell_x[i])
            continue
        covered = math.fsum(v for _, v in pieces)
        if labels == "uniform":
            loading = mean_label
        else:
            loading = math.fsum(v * shell_x[i] for i, v in pieces) / covered
        target = math.fsum(v * shell_density[i] for i, v in pieces) / covered
        lo = min(shell_t[i] for i, _ in pieces)
        hi = max(shell_t[i] for i, _ in pieces)

        def difference(t):
            return evaluate_cell(t, water, loading, oil, model).energy_density_j_m3 - target

        pad = 0.05 * (hi - lo) + 0.05
        f_lo, f_hi = difference(lo - pad), difference(hi + pad)
        if f_lo * f_hi > 0.0:
            raise FeedShellsRefused(f"master cell {j}: the energy temperature is not bracketed")
        temperature = brentq(difference, lo - pad, hi + pad, xtol=1e-13, rtol=8.881784197001252e-16)
        temperatures.append(temperature)
        loadings.append(loading)
        inverted.append(
            {
                "cell": j,
                "faces_m": [a, b],
                "shells": [i for i, _ in pieces],
                "overlap_fractions": [v / grid.volumes[j] for _, v in pieces],
                "temperature_k": temperature,
                "historical_loading": loading,
                "energy_density_residual_j_m3": difference(temperature),
            }
        )
    report = {
        "map": "conservative remap of the shell profile onto the master grid",
        "labels": labels,
        "mean_historical_loading": mean_label,
        "shell_historical_loadings": list(shell_x),
        "largest_label_shift": max(abs(x - mean_label) for x in shell_x),
        "shell_faces_m": list(faces),
        "inverted_cells": inverted,
        "physically_qualifying": False,
    }
    return tuple(temperatures), tuple(loadings), report
