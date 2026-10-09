r"""Gate 1d: resolved wet-core conduction and shell-conserved n-hexane.

This module implements the two deliberately separate wet-side closures:

* PHY-038: at radial activation, each material shell stores the historical
  loading ``X_h,w^a = X_c(theta_a)``.  The wet core has no pore-gas or liquid
  diffusion operator.  Retained and mobile liquid locally repartition at
  ``a_h=1`` while their sum remains the immutable activation label.
* PHY-013: the complete common-datum internal-energy density is advanced by
  spherical conduction on a front-fitted moving grid.

The ALE discretisation uses exact swept spherical face volumes.  Thus its mesh
flux satisfies the discrete geometric conservation law exactly, including for
finite front motion; replacing it by ``A_new * face_speed`` would not.

The wet energy per kg dry meal is

    e_w = h_dm(T,w_o) + X_w h_w,l(T,P)
          + X_h,w^a u_h,l(T,P) - B_h(a_h=1,T,w_o),

because ``U_ret + U_mobile = X_h,w^a u_h,l - B_h`` under the PHY-051 shared
liquid-like retained/mobile volume convention.  This form includes binding
once and makes the retained/mobile repartition explicitly source-free.
``h_dm`` uses PHY-048's reference-anchored native-meal/oil weights; the native
meal heat capacity already contains its reference oil and is counted once.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from dtdc_simulator.core2.particle._exact_memo import exact_memo
from dtdc_simulator.core2.particle.grid import SphericalGrid
from dtdc_simulator.core2.props import critical_volume as cv
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa


@dataclass(frozen=True)
class WetCoreParams:
    """Introspectable wet-core properties and frozen soybean priors."""

    rho_dm_p: float = 1159.65
    epsilon_p: float = 0.141
    pressure_pa: float = 101325.0
    conductivity: float = 0.24
    cp_dry_meal: float = sp.CP_NATIVE_DRY_MEAL
    cp_oil: float = sp.CP_OIL
    T_ref_solid: float = 298.15
    X_water: float = 0.10
    w_o: float = sp.W_O_REF
    gab: sp.GabParams = sp.GabParams()
    oil: sp.OilIsotherm = sp.OilIsotherm()
    # a_h=1 requires K(T)<1 for the frozen GAB branch; its singular threshold
    # is about 305 K.  That is a real mathematical wall and it anchors T_min.
    # The CEILING has no matching wall: K(T) falls monotonically with T, so the
    # a_h=1 singularity is only ever approached from BELOW.  T_max is therefore
    # a chosen VALIDATION SCOPE, not a physical limit.
    #
    # RAISED 360.0 -> 380.0 on 2026-09-03 under owner ruling Q-PD-1 (RATIFIED,
    # docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md, "RULED 2026-09-03 evening"),
    # WS-B of the PD phase-change program
    # (docs/GT_PS2_PD_PHASE_CHANGE_PROGRAM_DESIGN_2026-09-03.md, section 2.2).
    # Ceremony record, with the measured validity sweep and the pin cascade:
    # docs/GT_PS2_PD_WSB_BRACKET_EXTENSION_CEREMONY_2026-09-03.md.
    #
    # Measured over 360-385 K at the frozen soybean priors: X_c(T) stays
    # strictly positive and strictly decreasing with no new extremum,
    # K(T)*a_h reaches only 0.6960 at 360 K and 0.6256 at 380 K (strictly
    # below 1 with a wide margin), and the Span-Wagner liquid n-hexane branch
    # still resolves at 101 325 Pa throughout.  The extension is what lets
    # `pd_k_dome_bypass_tray_integration._cached_hexane_saturation_temperature_k`
    # find the 371.5 K design-duty root, which the 360.0 ceiling refused.
    #
    # NOT bit-neutral for bracket-SEEDED inversions: the ceiling also seeds
    # every bisection bracket in this module, in `sphere.py`, in
    # `high_loading_packet_adapter.py` and in
    # `engineering_subcritical_packet_state.py`, so converged temperatures
    # inside the OLD range shift within each solver's own declared resolution.
    # The ceremony record measures all 33 of them.
    T_min: float = 310.0
    T_max: float = 380.0


@dataclass(frozen=True)
class ContinuedWaterWetCoreParams(WetCoreParams):
    """Native wet storage with the explicitly continued positive-water domain."""

    luikov: sp.ContinuedPositiveLuikovParams = sp.ContinuedPositiveLuikovParams()


@dataclass(frozen=True)
class WetPartition:
    """Algebraic wet-side n-hexane partition, kg/kg dry meal."""

    total_historical: float
    retained: float
    mobile_liquid: float
    pore_gas_specific_volume: float
    condensed_specific_volume: float


@dataclass(frozen=True)
class WetCoreInventory:
    """Immutable material-shell labels stored once at the activation event."""

    activation_temperature: float
    loadings: tuple[float, ...]
    activation_specific_volumes: tuple[float, ...]


def activate(grid: SphericalGrid, T_activation: float, p: WetCoreParams) -> WetCoreInventory:
    """Store ``X_c(theta_a)`` once in every shell (PHY-035/038/050).

    No later API recalculates these loadings from the current temperature.
    """
    _validate_params(p)
    if not math.isfinite(T_activation) or not p.T_min <= T_activation <= p.T_max:
        raise ValueError(
            f"activation temperature {T_activation} outside wet-core caloric bracket "
            f"[{p.T_min}, {p.T_max}]"
        )
    X_a = cv.critical_hexane_loading(
        T_activation,
        p.epsilon_p,
        p.rho_dm_p,
        pressure_pa=p.pressure_pa,
    )
    # Fail at the event itself if the selected capacity cannot support the
    # required active mobile-liquid branch; do not create a nominally wet state
    # that only fails on its first energy evaluation.
    activated_partition = partition(T_activation, X_a, p)
    if activated_partition.mobile_liquid <= 1.0e-13 * max(abs(X_a), 1.0):
        raise ValueError("radial activation requires strictly positive mobile liquid")
    nu_a = cv.effective_critical_specific_volume(p.epsilon_p, p.rho_dm_p)
    return WetCoreInventory(
        activation_temperature=T_activation,
        loadings=(X_a,) * grid.n,
        activation_specific_volumes=(nu_a,) * grid.n,
    )


def partition(T: float, X_h_a: float, p: WetCoreParams) -> WetPartition:
    """Repartition a conserved wet-shell loading at saturated activity.

    The pore-gas volume is exactly zero.  If retained capacity alone exceeds
    the historical total, the one-front liquid-filled topology has failed and
    the state is rejected rather than clipped.
    """
    W = sp.retained_hexane(1.0, T, p.w_o, p.gab, p.oil)
    X_liq = X_h_a - W
    scale = max(abs(X_h_a), abs(W), 1.0)
    if X_liq < -1.0e-13 * scale:
        raise ValueError(
            "wet-shell mobile liquid disappeared away from the front: "
            f"X_h,w^a={X_h_a}, W_eq(a_h=1)={W}; reject, do not clamp"
        )
    if X_liq < 0.0:
        X_liq = 0.0
    rho_l = hx.state_Tp(T, p.pressure_pa, "liquid").rho_mass
    return WetPartition(
        total_historical=X_h_a,
        retained=W,
        mobile_liquid=X_liq,
        pore_gas_specific_volume=0.0,
        condensed_specific_volume=X_h_a / rho_l,
    )


def volume_departure(
    T: float, X_h_a: float, activation_specific_volume: float, p: WetCoreParams
) -> float:
    """Post-activation condensed-volume departure, m3/kg dry meal.

    This is a gate diagnostic, never a mass-reset equation.  Gate 1h applies
    the fixed-geometry tolerance over the phase-feasible active-core range.
    """
    return partition(T, X_h_a, p).condensed_specific_volume - activation_specific_volume


def total_hexane(
    grid: SphericalGrid, inventory: WetCoreInventory | Sequence[float], p: WetCoreParams
) -> float:
    """Wet-domain historical n-hexane inventory, kg."""
    X = _loadings(inventory, grid.n)
    return sum(V * p.rho_dm_p * x for V, x in zip(grid.volumes, X))


# ---------------------------------------------------------------------------
# Common-datum wet caloric state and analytic temperature derivative
# ---------------------------------------------------------------------------
_GL5_NODES = (
    -0.906179845938664,
    -0.538469310105683,
    0.0,
    0.538469310105683,
    0.906179845938664,
)
_GL5_WEIGHTS = (
    0.236926885056189,
    0.478628670499366,
    0.568888888888889,
    0.478628670499366,
    0.236926885056189,
)


def _integrate_activity_0_1(function, panels: int = 24) -> float:
    total = 0.0
    width = 1.0 / panels
    for panel in range(panels):
        left = panel * width
        for node, weight in zip(_GL5_NODES, _GL5_WEIGHTS):
            a = left + 0.5 * width * (node + 1.0)
            total += 0.5 * width * weight * function(a)
    return total


@exact_memo(ambient=lambda: (hx.T_REF, hx.H_REF, hx.S_REF))
def saturated_binding_deficit(T: float, p: WetCoreParams) -> float:
    """B_h at a_h=1, J/kg dry meal, without inverse-isotherm quadrature.

    Changing variables from loading to activity gives the exact identity
    ``B = -R*T^2/M * integral_0^1 (dW/dT)_a / a da``.
    Exact-argument memoized (GT_PS2_EXACT_MEMOIZATION_RULING 2026-08-30).
    """
    integral = _integrate_activity_0_1(
        lambda a: sp.retained_hexane_dT(a, T, p.w_o, p.gab, p.oil) / a
    )
    return -sp.R * T * T * integral / hx.M


def saturated_binding_deficit_dT(T: float, p: WetCoreParams) -> float:
    """Analytic total derivative dB_h(a_h=1,T)/dT, J/(kg dry K)."""
    first = _integrate_activity_0_1(lambda a: sp.retained_hexane_dT(a, T, p.w_o, p.gab, p.oil) / a)
    second = _integrate_activity_0_1(
        lambda a: sp.retained_hexane_dTT(a, T, p.w_o, p.gab, p.oil) / a
    )
    return -sp.R * (2.0 * T * first + T * T * second) / hx.M


def _hexane_liquid_u_du(T: float, pressure_pa: float) -> tuple[float, float]:
    """Return (u, du/dT|P) in J/kg and J/(kg K), analytically from Helmholtz."""
    state = hx.state_Tp(T, pressure_pa, "liquid")
    rho = state.rho
    res = hx.residual(rho / hx.RHOC, hx.TC / T)
    dp_dT_rho = rho * hx.R * (res.compressibility - res.delta * res.tau * res.phir_dt)
    dp_drho_T = hx.R * T * (1.0 + 2.0 * res.delta * res.phir_d + res.delta**2 * res.phir_dd)
    dv_mass_dT_P = dp_dT_rho / (hx.M * rho * rho * dp_drho_T)
    du_dT_P = state.cp / hx.M - pressure_pa * dv_mass_dT_P
    return state.u_mass, du_dT_P


def hexane_liquid_specific_energy(T: float, p: WetCoreParams) -> float:
    """Common-datum liquid n-hexane internal energy, J/kg.

    This public scalar is also the energy authority for an attached mobile-
    liquid node.  Retained-solvent binding belongs to the meal-side caloric
    state and must not be added to this free-liquid energy.
    """
    if not math.isfinite(T) or not p.T_min <= T <= p.T_max:
        raise ValueError("liquid n-hexane temperature is outside the caloric bracket")
    return _hexane_liquid_u_du(T, p.pressure_pa)[0]


def hexane_liquid_specific_heat_capacity(T: float, p: WetCoreParams) -> float:
    """``du_h,l/dT|P`` for an attached mobile-liquid node, J/(kg K)."""
    if not math.isfinite(T) or not p.T_min <= T <= p.T_max:
        raise ValueError("liquid n-hexane temperature is outside the caloric bracket")
    capacity = _hexane_liquid_u_du(T, p.pressure_pa)[1]
    if capacity <= 0.0:
        raise ValueError("liquid n-hexane caloric capacity must be positive")
    return capacity


def hexane_liquid_density_and_derivative(T: float, p: WetCoreParams) -> tuple[float, float]:
    """Return ``rho_h,l`` and ``d rho_h,l/dT|P`` in mass-specific units."""
    if not math.isfinite(T) or not p.T_min <= T <= p.T_max:
        raise ValueError("liquid n-hexane temperature is outside the caloric bracket")
    state = hx.state_Tp(T, p.pressure_pa, "liquid")
    rho = state.rho
    res = hx.residual(rho / hx.RHOC, hx.TC / T)
    dp_dT_rho = rho * hx.R * (res.compressibility - res.delta * res.tau * res.phir_dt)
    dp_drho_T = hx.R * T * (1.0 + 2.0 * res.delta * res.phir_d + res.delta**2 * res.phir_dd)
    specific_volume_derivative = dp_dT_rho / (hx.M * rho * rho * dp_drho_T)
    density_derivative = -state.rho_mass**2 * specific_volume_derivative
    return state.rho_mass, density_derivative


def _water_liquid_u_du(T: float, pressure_pa: float) -> tuple[float, float]:
    """Pure-liquid-partial-volume water sensitivity, not the nominal ledger.

    G1G-01 selects zero retained-water caloric partial volume nominally.  This
    helper is retained only for the mandatory upper-endpoint diagnostic where
    retained water is assigned the pure-liquid specific volume.
    """
    state = wa.state_Tp(T, pressure_pa, "liquid")
    rho = state.rho_mass
    res = wa.residual(rho / wa.RHOC, wa.TC / T)
    dp_dT_rho = rho * wa.R * (res.compressibility - res.delta * res.tau * res.phir_dt)
    dp_drho_T = wa.R * T * (1.0 + 2.0 * res.delta * res.phir_d + res.delta**2 * res.phir_dd)
    dv_mass_dT_P = dp_dT_rho / (rho * rho * dp_drho_T)
    du_dT_P = state.cp_mass - pressure_pa * dv_mass_dT_P
    return state.u_mass, du_dT_P


def _water_retained_h_dh(T: float, pressure_pa: float) -> tuple[float, float]:
    """Nominal retained-water ``(u, du/dT)=(h_l, cp_l)``, mass-specific.

    This is the G1G-01 zero-caloric-partial-volume mixture convention.  The
    ordinary vaporisation enthalpy is then present exactly once in the phase
    state difference; no latent source is added to the energy equation.
    """
    state = wa.state_Tp(T, pressure_pa, "liquid")
    return state.h_mass, state.cp_mass


@exact_memo(ambient=lambda: (hx.T_REF, hx.H_REF, hx.S_REF))
def specific_energy(T: float, X_h_a: float, p: WetCoreParams) -> float:
    """Complete wet-core internal energy, J/kg dry meal.

    Exact-argument memoized (GT_PS2_EXACT_MEMOIZATION_RULING 2026-08-30).
    """
    partition(T, X_h_a, p)  # fail closed if the saturated wet topology is infeasible
    u_h, _ = _hexane_liquid_u_du(T, p.pressure_pa)
    h_w_retained, _ = _water_retained_h_dh(T, p.pressure_pa)
    return (
        sp.composite_dry_meal_sensible_energy(T, p.T_ref_solid, p.w_o, p.cp_dry_meal, p.cp_oil)
        + p.X_water * h_w_retained
        + X_h_a * u_h
        - saturated_binding_deficit(T, p)
    )


def specific_heat_capacity(T: float, X_h_a: float, p: WetCoreParams) -> float:
    """Analytic de_w/dT at fixed historical shell inventories, J/(kg dry K)."""
    partition(T, X_h_a, p)
    _, du_h = _hexane_liquid_u_du(T, p.pressure_pa)
    _, dh_w_retained = _water_retained_h_dh(T, p.pressure_pa)
    capacity = (
        sp.composite_dry_meal_heat_capacity(p.w_o, p.cp_dry_meal, p.cp_oil)
        + p.X_water * dh_w_retained
        + X_h_a * du_h
        - saturated_binding_deficit_dT(T, p)
    )
    if capacity <= 0.0:
        raise ValueError(f"non-positive wet-core caloric capacity {capacity} J/(kg dry K)")
    return capacity


def energy_density(T: float, X_h_a: float, p: WetCoreParams) -> float:
    """Wet-core internal-energy density, J/m3 particle."""
    return p.rho_dm_p * specific_energy(T, X_h_a, p)


def energy_capacity(T: float, X_h_a: float, p: WetCoreParams) -> float:
    """Analytic wet-mixture energy derivative, J/(m3 particle K)."""
    return p.rho_dm_p * specific_heat_capacity(T, X_h_a, p)


def temperature_from_energy_density(U: float, X_h_a: float, p: WetCoreParams) -> float:
    """Positive-capacity caloric inverse U -> T (bracketed, no clamp)."""
    lo, hi = p.T_min, p.T_max
    U_lo = energy_density(lo, X_h_a, p)
    U_hi = energy_density(hi, X_h_a, p)
    if not U_lo <= U <= U_hi:
        raise ValueError(
            f"wet-core energy density {U} outside caloric bracket [{U_lo}, {U_hi}]; "
            "reject, do not clamp"
        )
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if energy_density(mid, X_h_a, p) < U:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-10:
            break
    return 0.5 * (lo + hi)


def total_energy(
    grid: SphericalGrid,
    temperatures: Sequence[float],
    inventory: WetCoreInventory | Sequence[float],
    p: WetCoreParams,
) -> float:
    """Integrated wet-domain internal energy, J."""
    _require_length("temperatures", temperatures, grid.n)
    X = _loadings(inventory, grid.n)
    return sum(V * energy_density(T, x, p) for V, T, x in zip(grid.volumes, temperatures, X))


def isothermal_temperature(
    grid: SphericalGrid,
    energy: float,
    inventory: WetCoreInventory | Sequence[float],
    p: WetCoreParams,
) -> float:
    """Faner limiting oracle: uniform wet-core T at a specified total energy."""
    X = _loadings(inventory, grid.n)

    def energy_at(T: float) -> float:
        return sum(V * energy_density(T, x, p) for V, x in zip(grid.volumes, X))

    lo, hi = p.T_min, p.T_max
    E_lo, E_hi = energy_at(lo), energy_at(hi)
    if not E_lo <= energy <= E_hi:
        raise ValueError(
            f"wet-core total energy {energy} outside isothermal caloric bracket "
            f"[{E_lo}, {E_hi}]; reject, do not clamp"
        )
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if energy_at(mid) < energy:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-10:
            break
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Conservative front-fitted ALE conduction
# ---------------------------------------------------------------------------
def _swept_volume_rates(old: SphericalGrid, new: SphericalGrid, dt: float) -> list[float]:
    return [
        4.0 / 3.0 * math.pi * (r_new**3 - r_old**3) / dt
        for r_old, r_new in zip(old.faces, new.faces)
    ]


def residual_and_jacobian(
    grid_prev: SphericalGrid,
    grid: SphericalGrid,
    T_prev: Sequence[float],
    T: Sequence[float],
    dt: float,
    inventory: WetCoreInventory | Sequence[float],
    p: WetCoreParams,
    surface=("closed",),
    source: Sequence[float] | None = None,
    interface_temperature: float | None = None,
) -> tuple[list[float], list[float], list[float], list[float], float]:
    """Return backward-Euler ALE residual, analytic tridiagonal Jacobian,
    and the outward surface energy rate relative to the moving boundary (W).
    """
    _validate_step_inputs(grid_prev, grid, T_prev, T, dt, p)
    X = _loadings(inventory, grid.n)
    S = [0.0] * grid.n if source is None else list(source)
    _require_length("source", S, grid.n)
    if not all(math.isfinite(value) for value in S):
        raise ValueError("wet-core energy sources must be finite")
    U_prev = [energy_density(t, x, p) for t, x in zip(T_prev, X)]
    U = [energy_density(t, x, p) for t, x in zip(T, X)]
    C = [energy_capacity(t, x, p) for t, x in zip(T, X)]
    residual = [
        grid.volumes[i] * U[i] - grid_prev.volumes[i] * U_prev[i] - dt * grid.volumes[i] * S[i]
        for i in range(grid.n)
    ]
    sub = [0.0] * grid.n
    diag = [grid.volumes[i] * C[i] for i in range(grid.n)]
    sup = [0.0] * grid.n
    swept = _swept_volume_rates(grid_prev, grid, dt)

    # Each interior extensive face rate H is outward from cell left to right:
    # H = A*q - Q_mesh*u_upwind.  Primary drainage has Q_mesh <= 0, so the
    # relative flow is outward and the left cell is the upwind wet material.
    # The sign branch also makes the operator correct under local ALE rezoning.
    for face in range(1, grid.n):
        left, right = face - 1, face
        dr = grid.centers[right] - grid.centers[left]
        conductance = grid.areas[face] * p.conductivity / dr
        if swept[face] <= 0.0:
            mesh_energy = -swept[face] * U[left]
            dmesh_left, dmesh_right = -swept[face] * C[left], 0.0
        else:
            mesh_energy = -swept[face] * U[right]
            dmesh_left, dmesh_right = 0.0, -swept[face] * C[right]
        H = conductance * (T[left] - T[right]) + mesh_energy
        dH_left = conductance + dmesh_left
        dH_right = -conductance + dmesh_right
        residual[left] += dt * H
        residual[right] -= dt * H
        diag[left] += dt * dH_left
        sup[left] += dt * dH_right
        sub[right] -= dt * dH_left
        diag[right] -= dt * dH_right

    last = grid.n - 1
    bc = surface[0]
    if bc == "closed":
        H_conduction = 0.0
        dH_conduction = 0.0
    elif bc == "dirichlet":
        T_surface = surface[1]
        if not math.isfinite(T_surface) or not p.T_min <= T_surface <= p.T_max:
            raise ValueError("wet-core Dirichlet temperature is outside the caloric bracket")
        dr = grid.R - grid.centers[last]
        conductance = grid.areas[-1] * p.conductivity / dr
        H_conduction = conductance * (T[last] - T_surface)
        dH_conduction = conductance
    elif bc == "flux":
        q_surface = surface[1]
        if not math.isfinite(q_surface):
            raise ValueError("wet-core surface heat flux must be finite")
        H_conduction = grid.areas[-1] * q_surface
        dH_conduction = 0.0
    elif bc == "robin":
        if len(surface) != 3:
            raise ValueError("wet-core Robin boundary requires ('robin', h, T_infinity)")
        heat_transfer_coefficient, ambient_temperature = surface[1:]
        conductance = finite_film_conductance(grid, p, heat_transfer_coefficient)
        if not math.isfinite(ambient_temperature):
            raise ValueError("wet-core Robin ambient temperature must be finite")
        # The external film and the outer half cell are resistances in series.
        # H is positive outward, consistently with every other surface BC.
        H_conduction = conductance * (T[last] - ambient_temperature)
        dH_conduction = conductance
    else:
        raise ValueError(f"unknown wet-core surface BC {bc!r}")
    if interface_temperature is None and bc == "dirichlet":
        interface_temperature = surface[1]
    if (
        interface_temperature is not None
        and bc == "dirichlet"
        and interface_temperature != surface[1]
    ):
        raise ValueError("Dirichlet and explicit wet-core interface temperatures must be identical")
    if interface_temperature is None:
        U_interface = U[last]
        dU_interface_dcell = C[last]
    else:
        if (
            not math.isfinite(interface_temperature)
            or not p.T_min <= interface_temperature <= p.T_max
        ):
            raise ValueError("wet-core interface temperature is outside the caloric bracket")
        U_interface = energy_density(interface_temperature, X[last], p)
        # The interface trace is a separate coupled unknown/value; this local
        # tridiagonal Jacobian differentiates only with respect to cell T.
        dU_interface_dcell = 0.0
    surface_rate = H_conduction - swept[-1] * U_interface
    residual[last] += dt * surface_rate
    diag[last] += dt * (dH_conduction - swept[-1] * dU_interface_dcell)
    return residual, sub, diag, sup, surface_rate


def finite_film_conductance(
    grid: SphericalGrid, p: WetCoreParams, heat_transfer_coefficient: float
) -> float:
    """Outer-face conductance for a finite film in series with a half cell.

    The coefficient is based on the physical particle surface area.  Keeping
    ``(R-r_N)/k`` in the resistance is required for radial mesh convergence;
    applying ``h`` directly to the outer cell-center temperature would silently
    remove the final half-cell conduction resistance.
    """
    if not math.isfinite(heat_transfer_coefficient) or heat_transfer_coefficient <= 0.0:
        raise ValueError("wet-core Robin heat-transfer coefficient must be positive")
    radial_resistance = (grid.R - grid.centers[-1]) / p.conductivity
    return grid.areas[-1] / (radial_resistance + 1.0 / heat_transfer_coefficient)


def finite_film_surface_temperature(
    grid: SphericalGrid,
    cell_temperature: float,
    p: WetCoreParams,
    heat_transfer_coefficient: float,
    ambient_temperature: float,
) -> float:
    """Reconstruct the physical surface trace for the series-film boundary."""
    if not math.isfinite(cell_temperature) or not math.isfinite(ambient_temperature):
        raise ValueError("finite-film temperatures must be finite")
    conductance_per_area = (
        finite_film_conductance(grid, p, heat_transfer_coefficient) / grid.areas[-1]
    )
    outward_heat_flux = conductance_per_area * (cell_temperature - ambient_temperature)
    return ambient_temperature + outward_heat_flux / heat_transfer_coefficient


def step(
    grid_prev: SphericalGrid,
    grid: SphericalGrid,
    T_prev: Sequence[float],
    dt: float,
    inventory: WetCoreInventory | Sequence[float],
    p: WetCoreParams,
    surface=("closed",),
    source: Sequence[float] | None = None,
    interface_temperature: float | None = None,
    *,
    tol: float = 1.0e-11,
    max_iter: int = 30,
) -> list[float]:
    """Advance resolved wet-core conduction by one implicit ALE step."""
    T = list(T_prev)
    X = _loadings(inventory, grid.n)
    # A capacity-based scale is independent of the arbitrary common caloric
    # datum.  Scaling on |U| would make Newton convergence datum-dependent.
    scale = max(
        max(
            grid_prev.volumes[i] * energy_capacity(T_prev[i], x, p) * max(abs(T_prev[i]), 1.0),
            1.0e-12,
        )
        for i, x in enumerate(X)
    )
    for _ in range(max_iter):
        residual, sub, diag, sup, _ = residual_and_jacobian(
            grid_prev,
            grid,
            T_prev,
            T,
            dt,
            inventory,
            p,
            surface,
            source,
            interface_temperature,
        )
        norm = max(abs(value) for value in residual) / scale
        if norm < tol:
            return T
        delta = _thomas(sub, diag, sup, residual)
        damping = 1.0
        while damping > 1.0e-8:
            trial = [T[i] - damping * delta[i] for i in range(grid.n)]
            if all(p.T_min <= value <= p.T_max for value in trial):
                break
            damping *= 0.5
        if damping <= 1.0e-8:
            raise RuntimeError("wet-core energy Newton left the caloric temperature bracket")
        T = trial
    raise RuntimeError("wet-core energy Newton did not converge")


def steady_conduction_solve(
    grid: SphericalGrid, conductivity: float, source: Sequence[float], surface
) -> list[float]:
    """Linear steady conduction solve for spatial-order qualification."""
    _require_length("source", source, grid.n)
    sub = [0.0] * grid.n
    diag = [0.0] * grid.n
    sup = [0.0] * grid.n
    rhs = [-grid.volumes[i] * source[i] for i in range(grid.n)]
    for face in range(1, grid.n):
        left, right = face - 1, face
        conductance = grid.areas[face] * conductivity / (grid.centers[right] - grid.centers[left])
        diag[left] -= conductance
        sup[left] += conductance
        sub[right] += conductance
        diag[right] -= conductance
    if surface[0] == "dirichlet":
        conductance = grid.areas[-1] * conductivity / (grid.R - grid.centers[-1])
        diag[-1] -= conductance
        rhs[-1] -= conductance * surface[1]
    elif surface[0] != "closed":
        raise ValueError("steady conduction qualifier supports closed or dirichlet surface")
    return _thomas(sub, diag, sup, rhs)


def _loadings(inventory: WetCoreInventory | Sequence[float], expected: int) -> tuple[float, ...]:
    values = inventory.loadings if isinstance(inventory, WetCoreInventory) else tuple(inventory)
    _require_length("wet-shell loadings", values, expected)
    return values


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


def _validate_params(p: WetCoreParams) -> None:
    scalar_values = (
        p.rho_dm_p,
        p.epsilon_p,
        p.pressure_pa,
        p.conductivity,
        p.cp_dry_meal,
        p.cp_oil,
        p.T_ref_solid,
        p.X_water,
        p.w_o,
        p.T_min,
        p.T_max,
    )
    if not all(math.isfinite(value) for value in scalar_values):
        raise ValueError("wet-core parameters must be finite")
    if p.rho_dm_p <= 0.0 or p.conductivity <= 0.0 or p.cp_dry_meal <= 0.0 or p.cp_oil <= 0.0:
        raise ValueError(
            "wet-core density, conductivity, native-meal cp, and oil cp must be positive"
        )
    if not 0.0 < p.epsilon_p < 1.0:
        raise ValueError("wet-core porosity must lie strictly between zero and one")
    if p.pressure_pa <= 0.0:
        raise ValueError("wet-core pressure must be positive")
    water_bounds = p.luikov if isinstance(p, ContinuedWaterWetCoreParams) else sp.LuikovParams()
    if isinstance(p, ContinuedWaterWetCoreParams) and p.X_water <= 0.0:
        raise ValueError(
            "continued wet material requires positive water; zero needs absent topology"
        )
    if not water_bounds.W_ref <= p.X_water <= water_bounds.W_cap:
        raise ValueError(
            "wet-core retained-water loading must lie on the qualified "
            f"positive-moisture Luikov branch [{water_bounds.W_ref}, {water_bounds.W_cap}] "
            "inside the evidence cap"
        )
    if not sp.W_O_MIN <= p.w_o <= sp.W_O_MAX:
        raise ValueError("wet-core residual-oil fraction is outside the PHY-048 envelope")
    if not p.T_min < p.T_max:
        raise ValueError("wet-core caloric temperature bracket is empty")
    if p.gab.K(p.T_min) >= 1.0:
        raise ValueError("wet-core lower temperature reaches the a_h=1 GAB singularity")


def _validate_step_inputs(
    old: SphericalGrid,
    new: SphericalGrid,
    T_prev: Sequence[float],
    T: Sequence[float],
    dt: float,
    p: WetCoreParams,
) -> None:
    _validate_params(p)
    if old.n != new.n:
        raise ValueError("front-fitted wet grids must retain their cell count")
    if new.R > old.R * (1.0 + 1.0e-14):
        raise ValueError("primary-drainage wet front cannot expand")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("wet-core time step must be positive and finite")
    _require_length("previous temperatures", T_prev, old.n)
    _require_length("trial temperatures", T, new.n)
    temperatures = (*T_prev, *T)
    if not all(math.isfinite(value) and p.T_min <= value <= p.T_max for value in temperatures):
        raise ValueError("wet-core cell temperature is outside the caloric bracket")


def _thomas(sub, diag, sup, rhs):
    """Solve a tridiagonal system without mutating its inputs."""
    n = len(rhs)
    upper = [0.0] * n
    transformed = [0.0] * n
    upper[0] = sup[0] / diag[0]
    transformed[0] = rhs[0] / diag[0]
    for i in range(1, n):
        pivot = diag[i] - sub[i] * upper[i - 1]
        upper[i] = sup[i] / pivot if i < n - 1 else 0.0
        transformed[i] = (rhs[i] - sub[i] * transformed[i - 1]) / pivot
    solution = [0.0] * n
    solution[-1] = transformed[-1]
    for i in range(n - 2, -1, -1):
        solution[i] = transformed[i] - upper[i] * solution[i + 1]
    return solution
