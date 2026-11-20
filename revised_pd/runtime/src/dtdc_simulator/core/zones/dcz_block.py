"""Phase-active block-implicit DCZ solver.

The semismooth residual and its analytic sparse Jacobian remain adjacent and
are verified as a pair. The independent affine phase-fixed solve lives in
``dcz_fixed_block`` and is re-exported here for historical compatibility.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu, spsolve

from dtdc_simulator.core import thermo
from dtdc_simulator.core.zones import particle as pt
from dtdc_simulator.core.zones.dcz_fixed_block import (
    _ThermalBlockResult as _ThermalBlockResult,
    _solve_phase_fixed_thermal_block as _solve_phase_fixed_thermal_block,
)
from dtdc_simulator.core.zones.dcz_types import DCZConstants

# Qualified interval used as a nonlinear-solver safeguard. It is deliberately
# much wider than the physical operating envelope and is not an active process
# constraint at accepted operating points.
_TP_MIN_K = 250.0
_TP_MAX_K = 480.0

# ---------------------------------------------------------------------------
# Coupled phase/component block (semismooth Newton).
#
# This is the phase-active generalisation of ``_solve_phase_fixed_thermal_block``.
# It solves particle radial temperatures, vapor temperature, vapor water flow,
# bound and free solid moisture, and one signed net water phase rate as one
# residual system with a scaled minimum-map complementarity.  Positive phase
# rate is condensation and negative phase rate is free-liquid evaporation.
# With sorption disabled and no active phase change it reduces to the affine
# phase-free block.  See
# ``docs/DCZ_BLOCK_IMPLICIT_THERMAL_COUPLING_DESIGN_2026-07-26.md`` for the
# residual and Jacobian derivation this code implements.
# ---------------------------------------------------------------------------

_DEW_BOUNDS_K = (230.0, 450.0)


def _dew_point_temperature_and_dwater(
    hexane_kg_s: float,
    water_kg_s: float,
    antoine_water: thermo.AntoineParams,
    pressure_pa: float,
) -> tuple[float, float]:
    """``T_dew(hexane/water)`` and ``d T_dew / d water`` at fixed hexane.

    The value is identical to :func:`thermo.dew_point_temperature`; the
    derivative is the exact analytic chain rule through the closed-form
    Antoine inverse.  Returns a zero derivative in the clamped regions, matching
    the value function's flat behaviour there.
    """
    lo, hi = _DEW_BOUNDS_K
    water = max(water_kg_s, 1.0e-12)
    y_v2 = hexane_kg_s / water
    rmw = thermo.M_WATER / thermo.M_HEXANE
    denom = 1.0 + rmw * y_v2
    y_water = 1.0 / denom
    target = y_water * pressure_pa
    if target <= 0.0:
        return lo, 0.0
    ln10 = math.log(10.0)
    log10_p_bar = math.log10(target / 1.0e5)
    base = antoine_water.A - log10_p_bar
    t_raw = antoine_water.B / base - antoine_water.C
    if t_raw <= lo or t_raw >= hi:
        return min(max(t_raw, lo), hi), 0.0
    d_y_v2_d_w = -hexane_kg_s / water**2
    d_y_water_d_y_v2 = -rmw * y_water**2
    d_target_d_y_water = pressure_pa
    d_log_d_target = 1.0 / (target * ln10)
    d_t_d_log = antoine_water.B / base**2
    d_t_d_w = (
        d_t_d_log
        * d_log_d_target
        * d_target_d_y_water
        * d_y_water_d_y_v2
        * d_y_v2_d_w
    )
    return t_raw, d_t_d_w


def _water_activity_and_grads(
    hexane_kg_s: float,
    water_kg_s: float,
    surface_temperature_k: float,
    antoine_water: thermo.AntoineParams,
    pressure_pa: float,
) -> tuple[float, float, float, bool]:
    """``a_w`` and its derivatives w.r.t. water flow and surface temperature.

    Returns ``(a_w_clamped, da_dw, da_dT, clamped)``.  ``a_w`` is clamped to the
    Luikov-validated range exactly as the partitioned solver clamps it; the
    ``clamped`` flag lets the Jacobian zero the derivative there so it stays
    consistent with the flat residual.
    """
    water = max(water_kg_s, 1.0e-12)
    y_v2 = hexane_kg_s / water
    rmw = thermo.M_WATER / thermo.M_HEXANE
    denom = 1.0 + rmw * y_v2
    y_water = 1.0 / denom
    p_sat = thermo.antoine_pressure_pa(surface_temperature_k, antoine_water)
    a_raw = y_water * pressure_pa / p_sat
    lo = 1.0e-9
    hi = thermo.LUIKOV_MAX_VALIDATED_UR
    if a_raw <= lo or a_raw >= hi:
        return min(max(a_raw, lo), hi), 0.0, 0.0, True
    ln10 = math.log(10.0)
    d_y_v2_d_w = -hexane_kg_s / water**2
    d_y_water_d_y_v2 = -rmw * y_water**2
    da_d_y_water = pressure_pa / p_sat
    da_d_w = da_d_y_water * d_y_water_d_y_v2 * d_y_v2_d_w
    # d p_sat/dT = p_sat * ln10 * B/(C+T)^2 ; a_w = num/p_sat.
    d_psat_d_t = p_sat * ln10 * antoine_water.B / (antoine_water.C + surface_temperature_k) ** 2
    da_d_t = -a_raw / p_sat * d_psat_d_t
    return a_raw, da_d_w, da_d_t, False


def _luikov_xe_and_slope(a_w: float, luikov: thermo.LuikovParams) -> tuple[float, float]:
    """``Xe(a_w)`` and ``dXe/da_w`` for the modified-Luikov isotherm."""
    denom = 1.0 - luikov.A2 * math.log(a_w)
    xe = luikov.A1 / denom
    d_xe = luikov.A1 * luikov.A2 / (a_w * denom**2)
    return xe, d_xe


@dataclass(frozen=True)
class _CoupledBlockContext:
    """Fixed inputs to the coupled phase/component residual and Jacobian."""

    nz: int
    layers: int
    particle_offsets: tuple[int, ...]
    particle_unknowns: int
    scalar_offset: int
    unknowns: int
    time_substeps: tuple[int, ...]
    cell_residence_s: tuple[float, ...]
    dt_sub: tuple[float, ...]
    # geometry
    volumes: tuple[float, ...]
    face_areas: tuple[float, ...]
    dr: float
    total_volume: float
    # material / transport
    k_mix: float
    base_cv: float
    water_cv: float
    surface_conductance_w_k: float
    reference_temperature_k: float
    m_dry_kg_s: float
    kappa_w: float
    bound_capacity: float
    lambda_scale: float
    dew_scale_k: float
    # fixed per-cell inputs
    vapor_hexane_kg_s: tuple[float, ...]
    q_indirect_w_m3: tuple[float, ...]
    q_axial_w_m3: tuple[float, ...]
    cell_volume_m3: tuple[float, ...]
    representative_particles: tuple[float, ...]
    sorption_source: tuple[tuple[float, ...], ...]
    # boundary conditions
    bottom_water_kg_s: float
    bottom_hexane_kg_s: float
    bottom_temperature_k: float
    moisture_in: float
    bound_in: float
    free_in: float
    initial_particle_tp: tuple[float, ...]
    # enthalpy datum
    ref: thermo.VaporEnthalpyRef
    cp_water_liquid: float
    antoine_water: thermo.AntoineParams
    luikov: thermo.LuikovParams
    pressure_pa: float

    def p_index(self, cell: int, substep: int, shell: int) -> int:
        return self.particle_offsets[cell] + substep * self.layers + shell

    def tv_index(self, cell: int) -> int:
        return self.scalar_offset + 5 * cell

    def w_index(self, cell: int) -> int:
        return self.scalar_offset + 5 * cell + 1

    def xb_index(self, cell: int) -> int:
        return self.scalar_offset + 5 * cell + 2

    def xf_index(self, cell: int) -> int:
        return self.scalar_offset + 5 * cell + 3

    def lam_index(self, cell: int) -> int:
        """Index of the signed net phase rate (legacy accessor name)."""
        return self.scalar_offset + 5 * cell + 4


def _coupled_cell_derived(u: np.ndarray, ctx: _CoupledBlockContext, j: int) -> dict:
    """Shared per-cell derived quantities used by both residual and Jacobian."""
    tv = u[ctx.tv_index(j)]
    w = u[ctx.w_index(j)]
    xb = u[ctx.xb_index(j)]
    xf = u[ctx.xf_index(j)]
    lam = u[ctx.lam_index(j)]
    last = ctx.time_substeps[j] - 1
    surf_index = ctx.p_index(j, last, ctx.layers - 1)
    t_surf = u[surf_index]

    if j == 0:
        xb_prev = ctx.bound_in
        xf_prev = ctx.free_in
        x_prev = ctx.moisture_in
    else:
        xb_prev = u[ctx.xb_index(j - 1)]
        xf_prev = u[ctx.xf_index(j - 1)]
        x_prev = u[ctx.xb_index(j - 1)] + u[ctx.xf_index(j - 1)]

    hex_j = ctx.vapor_hexane_kg_s[j]
    # Bound-sorption LDF toward the local equilibrium moisture.
    a_w, da_dw, da_dt, a_clamped = _water_activity_and_grads(
        hex_j, w, t_surf, ctx.antoine_water, ctx.pressure_pa
    )
    xe, dxe_da = _luikov_xe_and_slope(a_w, ctx.luikov)
    decay = math.exp(-ctx.kappa_w * ctx.cell_residence_s[j])
    xb_ldf_raw = xe * (1.0 - decay) + xb_prev * decay
    xb_clamped = xb_ldf_raw <= 0.0 or xb_ldf_raw >= ctx.bound_capacity
    xb_ldf = min(max(xb_ldf_raw, 0.0), ctx.bound_capacity)
    # d(xb_ldf)/d(...) through Xe(a_w(w, t_surf)); zero in clamp regions.
    if xb_clamped or a_clamped:
        dxbldf_dw = 0.0
        dxbldf_dt = 0.0
    else:
        dxbldf_dxe = 1.0 - decay
        dxbldf_dw = dxbldf_dxe * dxe_da * da_dw
        dxbldf_dt = dxbldf_dxe * dxe_da * da_dt
    dxbldf_dxbprev = 0.0 if xb_clamped else decay

    t_dew, dtdew_dw = _dew_point_temperature_and_dwater(
        hex_j, w, ctx.antoine_water, ctx.pressure_pa
    )
    return {
        "tv": tv,
        "w": w,
        "xb": xb,
        "xf": xf,
        "lam": lam,
        "t_surf": t_surf,
        "surf_index": surf_index,
        "xb_prev": xb_prev,
        "xf_prev": xf_prev,
        "x_prev": x_prev,
        "hex_j": hex_j,
        "xb_ldf": xb_ldf,
        "dxbldf_dw": dxbldf_dw,
        "dxbldf_dt": dxbldf_dt,
        "dxbldf_dxbprev": dxbldf_dxbprev,
        "t_dew": t_dew,
        "dtdew_dw": dtdew_dw,
    }


def _coupled_block_residual(u: np.ndarray, ctx: _CoupledBlockContext) -> np.ndarray:
    """Residual vector of the coupled phase/component thermal block."""
    ref = ctx.ref
    tbw = ref.T_boil_water
    tbh = ref.T_boil_hexane
    cpwl = ctx.cp_water_liquid
    f = np.zeros(ctx.unknowns, dtype=float)

    def h_water(t: float) -> float:
        return ref.dH_vap_water + ref.cp_water_vapor * (t - tbw)

    def h_hexane(t: float) -> float:
        return ref.dH_vap_hexane + ref.cp_hexane_vapor * (t - tbh)

    # -- particle radial energy rows --
    for j in range(ctx.nz):
        substeps = ctx.time_substeps[j]
        dt_sub = ctx.dt_sub[j]
        cell = _coupled_cell_derived(u, ctx, j)
        x_prev = cell["x_prev"]
        x_here = u[ctx.xb_index(j)] + u[ctx.xf_index(j)]
        tv = cell["tv"]
        for step in range(substeps):
            a_old = step / substeps
            a_new = (step + 1) / substeps
            moisture_old = x_prev + (x_here - x_prev) * a_old
            moisture_new = x_prev + (x_here - x_prev) * a_new
            cv_old = ctx.base_cv + ctx.water_cv * moisture_old
            cv_new = ctx.base_cv + ctx.water_cv * moisture_new
            for i in range(ctx.layers):
                row = ctx.p_index(j, step, i)
                volume = ctx.volumes[i]
                storage_new = cv_new * volume / dt_sub
                storage_old = cv_old * volume / dt_sub
                t_here = u[row]
                res = storage_new * (t_here - ctx.reference_temperature_k)
                # radial conduction
                if i > 0:
                    cond = ctx.k_mix * ctx.face_areas[i - 1] / ctx.dr
                    res += cond * (t_here - u[ctx.p_index(j, step, i - 1)])
                if i < ctx.layers - 1:
                    cond = ctx.k_mix * ctx.face_areas[i] / ctx.dr
                    res += cond * (t_here - u[ctx.p_index(j, step, i + 1)])
                else:
                    res += ctx.surface_conductance_w_k * (t_here - tv)
                # residence advection / accumulation of previous state
                if step > 0:
                    t_prev = u[ctx.p_index(j, step - 1, i)]
                elif j > 0:
                    t_prev = u[ctx.p_index(j - 1, ctx.time_substeps[j - 1] - 1, i)]
                else:
                    t_prev = ctx.initial_particle_tp[i]
                res -= storage_old * (t_prev - ctx.reference_temperature_k)
                res -= (
                    ctx.sorption_source[j][i] + ctx.q_axial_w_m3[j]
                ) * volume
                f[row] = res

    # -- cell scalar rows: vapor enthalpy, water conservation, split, NCP --
    for j in range(ctx.nz):
        cell = _coupled_cell_derived(u, ctx, j)
        tv = cell["tv"]
        w = cell["w"]
        xb = u[ctx.xb_index(j)]
        xf = u[ctx.xf_index(j)]
        lam = cell["lam"]
        t_surf = cell["t_surf"]
        hex_j = cell["hex_j"]
        xb_prev = cell["xb_prev"]
        xf_prev = cell["xf_prev"]
        x_here = xb + xf
        x_prev = cell["x_prev"]

        if j < ctx.nz - 1:
            w_in = u[ctx.w_index(j + 1)]
            tv_in = u[ctx.tv_index(j + 1)]
            hex_in = ctx.vapor_hexane_kg_s[j + 1]
        else:
            w_in = ctx.bottom_water_kg_s
            tv_in = ctx.bottom_temperature_k
            hex_in = ctx.bottom_hexane_kg_s

        bound_rate = ctx.m_dry_kg_s * (xb - xb_prev)
        enthalpy_out = w * h_water(tv) + hex_j * h_hexane(tv)
        enthalpy_in = w_in * h_water(tv_in) + hex_in * h_hexane(tv_in)
        transferred_hex = (hex_j - hex_in) * (
            ref.dH_vap_hexane + ref.cp_hexane_vapor * (t_surf - tbh)
        )
        bound_liquid = bound_rate * cpwl * (t_surf - tbw)
        # One signed net phase rate prevents simultaneous condensation and
        # evaporation.  The liquid crosses the interface at the local binary
        # dew temperature; on a wet branch T_V=T_dew, while a dry superheated
        # branch can only evaporate the finite entering inventory.
        phase_liquid = lam * cpwl * (cell["t_dew"] - tbw)
        external_duty = ctx.q_indirect_w_m3[j] * ctx.cell_volume_m3[j]
        r_noninterface = (
            enthalpy_out
            - enthalpy_in
            - external_duty
            - transferred_hex
            + bound_liquid
            + phase_liquid
        )
        # paired interface flux (equal and opposite to the particle surface rows)
        rep = ctx.representative_particles[j]
        substeps = ctx.time_substeps[j]
        surf_mean = sum(
            u[ctx.p_index(j, s, ctx.layers - 1)] for s in range(substeps)
        ) / substeps
        f[ctx.tv_index(j)] = r_noninterface + rep * ctx.surface_conductance_w_k * (
            tv - surf_mean
        )

        # water conservation
        f[ctx.w_index(j)] = w - w_in + ctx.m_dry_kg_s * (x_here - x_prev)

        # bound-sorption LDF closure
        f[ctx.xb_index(j)] = xb - cell["xb_ldf"]

        # mechanism split (free water)
        f[ctx.xf_index(j)] = (
            ctx.m_dry_kg_s * (xf - xf_prev) - lam
        )

        # Wet/dry complementarity (scaled minimum map):
        #
        #   0 <= X_free,j  ⟂  T_V,j - T_dew,j >= 0.
        #
        # Together with m_dry*(X_free,j-X_free,j-1)=lambda_j, this is a
        # bounded signed phase law.  A positive lambda creates condensate; a
        # negative lambda evaporates at most the entering free inventory.  If
        # liquid remains the vapor is saturated, and if it is exhausted the
        # vapor may be superheated.  There is no separate condensation and
        # evaporation pair that could become algebraically degenerate.
        free_scale = ctx.lambda_scale / ctx.m_dry_kg_s
        a_hat = xf / free_scale
        b_hat = (tv - cell["t_dew"]) / ctx.dew_scale_k
        f[ctx.lam_index(j)] = min(a_hat, b_hat)

    return f


def _coupled_block_jacobian(u: np.ndarray, ctx: _CoupledBlockContext):
    """Analytic sparse Jacobian of :func:`_coupled_block_residual`."""
    ref = ctx.ref
    tbw = ref.T_boil_water
    cpwl = ctx.cp_water_liquid
    cpwv = ref.cp_water_vapor
    cphv = ref.cp_hexane_vapor
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    def add(r: int, col: int, value: float) -> None:
        if value != 0.0:
            rows.append(r)
            cols.append(col)
            data.append(value)

    def h_water(t: float) -> float:
        return ref.dH_vap_water + cpwv * (t - tbw)

    # -- particle rows --
    for j in range(ctx.nz):
        substeps = ctx.time_substeps[j]
        dt_sub = ctx.dt_sub[j]
        cell = _coupled_cell_derived(u, ctx, j)
        x_prev = cell["x_prev"]
        xb_here_index = ctx.xb_index(j)
        xf_here_index = ctx.xf_index(j)
        x_here = u[xb_here_index] + u[xf_here_index]
        for step in range(substeps):
            a_old = step / substeps
            a_new = (step + 1) / substeps
            moisture_old = x_prev + (x_here - x_prev) * a_old
            moisture_new = x_prev + (x_here - x_prev) * a_new
            cv_old = ctx.base_cv + ctx.water_cv * moisture_old
            cv_new = ctx.base_cv + ctx.water_cv * moisture_new
            for i in range(ctx.layers):
                row = ctx.p_index(j, step, i)
                volume = ctx.volumes[i]
                storage_new = cv_new * volume / dt_sub
                storage_old = cv_old * volume / dt_sub
                t_here = u[row]
                diagonal = storage_new
                if i > 0:
                    cond = ctx.k_mix * ctx.face_areas[i - 1] / ctx.dr
                    diagonal += cond
                    add(row, ctx.p_index(j, step, i - 1), -cond)
                if i < ctx.layers - 1:
                    cond = ctx.k_mix * ctx.face_areas[i] / ctx.dr
                    diagonal += cond
                    add(row, ctx.p_index(j, step, i + 1), -cond)
                else:
                    diagonal += ctx.surface_conductance_w_k
                    add(row, ctx.tv_index(j), -ctx.surface_conductance_w_k)
                add(row, row, diagonal)
                # previous-state advection term
                if step > 0:
                    prev_index = ctx.p_index(j, step - 1, i)
                    t_prev = u[prev_index]
                    add(row, prev_index, -storage_old)
                elif j > 0:
                    prev_index = ctx.p_index(j - 1, ctx.time_substeps[j - 1] - 1, i)
                    t_prev = u[prev_index]
                    add(row, prev_index, -storage_old)
                else:
                    t_prev = ctx.initial_particle_tp[i]
                # moisture dependence of storage (bilinear coupling)
                coeff_here = ctx.water_cv * volume / dt_sub * (
                    a_new * (t_here - ctx.reference_temperature_k)
                    - a_old * (t_prev - ctx.reference_temperature_k)
                )
                add(row, xb_here_index, coeff_here)
                add(row, xf_here_index, coeff_here)
                if j > 0:
                    coeff_prev = ctx.water_cv * volume / dt_sub * (
                        (1.0 - a_new) * (t_here - ctx.reference_temperature_k)
                        - (1.0 - a_old) * (t_prev - ctx.reference_temperature_k)
                    )
                    add(row, ctx.xb_index(j - 1), coeff_prev)
                    add(row, ctx.xf_index(j - 1), coeff_prev)

    # -- cell scalar rows --
    for j in range(ctx.nz):
        cell = _coupled_cell_derived(u, ctx, j)
        tv = cell["tv"]
        w = cell["w"]
        xb = u[ctx.xb_index(j)]
        xf = u[ctx.xf_index(j)]
        lam = cell["lam"]
        t_surf = cell["t_surf"]
        hex_j = cell["hex_j"]
        xb_prev = cell["xb_prev"]
        surf_index = cell["surf_index"]
        substeps = ctx.time_substeps[j]
        rep = ctx.representative_particles[j]
        sc = ctx.surface_conductance_w_k

        if j < ctx.nz - 1:
            w_in = u[ctx.w_index(j + 1)]
            tv_in = u[ctx.tv_index(j + 1)]
            hex_in = ctx.vapor_hexane_kg_s[j + 1]
            interior_in = True
        else:
            w_in = ctx.bottom_water_kg_s
            tv_in = ctx.bottom_temperature_k
            hex_in = ctx.bottom_hexane_kg_s
            interior_in = False

        # ---- vapor enthalpy row ----
        row = ctx.tv_index(j)
        add(row, row, w * cpwv + hex_j * cphv + rep * sc)
        add(
            row,
            ctx.w_index(j),
            h_water(tv) + lam * cpwl * cell["dtdew_dw"],
        )
        add(row, ctx.lam_index(j), cpwl * (cell["t_dew"] - tbw))
        bound_coeff = ctx.m_dry_kg_s * cpwl * (t_surf - tbw)
        add(row, ctx.xb_index(j), bound_coeff)
        if j > 0:
            add(row, ctx.xb_index(j - 1), -bound_coeff)
        # particle-surface dependence: transferred hexane enthalpy + bound liquid,
        # plus the paired interface flux (last substep gets both).
        surf_enth_coeff = (
            -(hex_j - hex_in) * cphv
            + ctx.m_dry_kg_s * (xb - xb_prev) * cpwl
        )
        for s in range(substeps):
            idx = ctx.p_index(j, s, ctx.layers - 1)
            coeff = -rep * sc / substeps
            if idx == surf_index:
                coeff += surf_enth_coeff
            add(row, idx, coeff)
        if interior_in:
            add(row, ctx.tv_index(j + 1), -(w_in * cpwv + hex_in * cphv))
            add(row, ctx.w_index(j + 1), -h_water(tv_in))

        # ---- water conservation row ----
        row = ctx.w_index(j)
        add(row, ctx.w_index(j), 1.0)
        add(row, ctx.xb_index(j), ctx.m_dry_kg_s)
        add(row, ctx.xf_index(j), ctx.m_dry_kg_s)
        if interior_in:
            add(row, ctx.w_index(j + 1), -1.0)
        if j > 0:
            add(row, ctx.xb_index(j - 1), -ctx.m_dry_kg_s)
            add(row, ctx.xf_index(j - 1), -ctx.m_dry_kg_s)

        # ---- bound-sorption LDF row ----
        row = ctx.xb_index(j)
        add(row, ctx.xb_index(j), 1.0)
        add(row, ctx.w_index(j), -cell["dxbldf_dw"])
        add(row, surf_index, -cell["dxbldf_dt"])
        if j > 0:
            add(row, ctx.xb_index(j - 1), -cell["dxbldf_dxbprev"])

        # ---- mechanism-split (free water) row ----
        row = ctx.xf_index(j)
        add(row, ctx.xf_index(j), ctx.m_dry_kg_s)
        add(row, ctx.lam_index(j), -1.0)
        if j > 0:
            add(row, ctx.xf_index(j - 1), -ctx.m_dry_kg_s)

        # ---- wet/dry complementarity (minimum map) row ----
        # The active branch is the smaller of the two scaled arguments; its
        # semismooth Jacobian is that branch's gradient.
        row = ctx.lam_index(j)
        free_scale = ctx.lambda_scale / ctx.m_dry_kg_s
        a_hat = xf / free_scale
        b_hat = (tv - cell["t_dew"]) / ctx.dew_scale_k
        if a_hat <= b_hat:
            add(row, ctx.xf_index(j), 1.0 / free_scale)
        else:
            add(row, ctx.tv_index(j), 1.0 / ctx.dew_scale_k)
            add(row, ctx.w_index(j), -cell["dtdew_dw"] / ctx.dew_scale_k)

    matrix = coo_matrix(
        (data, (rows, cols)), shape=(ctx.unknowns, ctx.unknowns)
    ).tocsr()
    return matrix


def _solve_coupled_newton_step_schur(
    matrix,
    residual: np.ndarray,
    particle_unknowns: int,
) -> np.ndarray:
    """Solve one exact Newton system by eliminating particle temperatures.

    The coupled Jacobian has the block form

    ``[A B; C D] [dTp; dy] = -[Fp; Fy]``,

    where ``A`` is the structured particle time/radial march and ``y`` holds
    only five scalar unknowns per axial cell.  Factoring the full matrix with
    SuperLU creates avoidable fill across the conjugate interface.  This
    routine factors ``A``, forms the exact Schur complement
    ``D - C A^-1 B``, solves the reduced scalar system, and back-substitutes.
    No equation, tolerance, or Newton linearisation is approximated.
    """
    n = matrix.shape[0]
    p = particle_unknowns
    if matrix.shape[1] != n or residual.shape != (n,):
        raise ValueError("coupled Newton system dimensions are inconsistent")
    if p <= 0 or p >= n:
        return np.asarray(spsolve(matrix.tocsc(), -residual), dtype=float)

    a = matrix[:p, :p].tocsc()
    b = matrix[:p, p:].tocsc()
    c = matrix[p:, :p].tocsr()
    d = matrix[p:, p:].toarray()
    rhs_p = -residual[:p]
    rhs_y = -residual[p:]

    particle_lu = splu(a)
    a_inv_rhs = np.asarray(particle_lu.solve(rhs_p), dtype=float)
    # SuperLU handles multiple dense right-hand sides in one call.  Fortran
    # ordering avoids an internal transpose/copy in the numeric solve.
    a_inv_b = np.asarray(
        particle_lu.solve(np.asfortranarray(b.toarray())),
        dtype=float,
    )
    schur = d - np.asarray(c @ a_inv_b, dtype=float)
    reduced_rhs = rhs_y - np.asarray(c @ a_inv_rhs, dtype=float)
    delta_y = np.linalg.solve(schur, reduced_rhs)
    delta_p = a_inv_rhs - a_inv_b @ delta_y
    return np.concatenate((delta_p, delta_y))


@dataclass(frozen=True)
class _CoupledBlockResult:
    particles: tuple[pt.ParticleState, ...]
    vapor_temperature_k: tuple[float, ...]
    vapor_water_kg_s: tuple[float, ...]
    bound_moisture: tuple[float, ...]
    free_moisture: tuple[float, ...]
    condensed_water_kg_s: tuple[float, ...]
    evaporated_water_kg_s: tuple[float, ...]
    particle_energy_terms_w: tuple[tuple[float, float, float, float], ...]
    vapor_row_residual_w: tuple[float, ...]
    residual_max_w: float
    interface_mismatch_max_w: float
    iterations: int
    converged: bool


def _build_coupled_block_context(
    *,
    initial_particle: pt.ParticleState,
    particles: list[pt.ParticleState],
    vapor_hexane_kg_s: list[float],
    bottom_water_kg_s: float,
    bottom_hexane_kg_s: float,
    bottom_temperature_k: float,
    moisture_in: float,
    free_in: float,
    q_indirect_w_m3: tuple[float, ...],
    q_axial_w_m3: tuple[float, ...],
    cell_widths_m: tuple[float, ...],
    cell_residence_s: tuple[float, ...],
    time_substeps: tuple[int, ...],
    bed_area_m2: float,
    geometry: pt.ShellGeometry,
    c: DCZConstants,
    kappa_w: float,
    bound_capacity: float,
    m_dry_kg_s: float,
    lambda_scale: float,
    dew_scale_k: float,
) -> _CoupledBlockContext:
    nz = len(particles)
    layers = c.particle.Np
    particle_offsets: list[int] = []
    particle_unknowns = 0
    for substeps in time_substeps:
        particle_offsets.append(particle_unknowns)
        particle_unknowns += substeps * layers
    scalar_offset = particle_unknowns
    unknowns = particle_unknowns + 5 * nz

    pc = c.particle
    k_mix = pc.alpha_pg * pc.k_pg + pc.alpha_ps * pc.k_ps
    base_cv = pc.alpha_pg * pc.rho_pg * pc.cp_pg + pc.alpha_ps * pc.rho_ps * pc.cp_ps
    water_cv = pc.alpha_ps * pc.rho_ps * pc.cp_water_liquid
    surface_area_m2 = geometry.face_areas[-1]
    surface_conductance_w_k = c.hQ * surface_area_m2
    reference_temperature_k = c.vapor_enthalpy_ref.T_boil_water

    # Hexane heat of sorption is a lagged input (fixed wpg2, lagged Tp), exactly
    # as in the phase-free block; the outer iteration converges the lag.
    sorption_source: list[tuple[float, ...]] = []
    for j in range(nz):
        entering_w = initial_particle.wpg2 if j == 0 else particles[j - 1].wpg2
        sorption_entering = pt.ParticleState(
            wpg2=entering_w,
            Tp=initial_particle.Tp if j == 0 else particles[j - 1].Tp,
        )
        sources = pt.finite_sorption_heat_source_per_layer_w_m3(
            sorption_entering, particles[j], cell_residence_s[j], pc
        )
        sorption_source.append(tuple(sources))

    cell_volume_m3 = tuple(bed_area_m2 * width for width in cell_widths_m)
    representative_particles = tuple(
        c.alpha_L * volume / geometry.total_volume for volume in cell_volume_m3
    )
    dt_sub = tuple(
        cell_residence_s[j] / time_substeps[j] for j in range(nz)
    )
    return _CoupledBlockContext(
        nz=nz,
        layers=layers,
        particle_offsets=tuple(particle_offsets),
        particle_unknowns=particle_unknowns,
        scalar_offset=scalar_offset,
        unknowns=unknowns,
        time_substeps=tuple(time_substeps),
        cell_residence_s=tuple(cell_residence_s),
        dt_sub=dt_sub,
        volumes=tuple(geometry.volumes),
        face_areas=tuple(geometry.face_areas),
        dr=geometry.dr,
        total_volume=geometry.total_volume,
        k_mix=k_mix,
        base_cv=base_cv,
        water_cv=water_cv,
        surface_conductance_w_k=surface_conductance_w_k,
        reference_temperature_k=reference_temperature_k,
        m_dry_kg_s=m_dry_kg_s,
        kappa_w=kappa_w,
        bound_capacity=bound_capacity,
        lambda_scale=lambda_scale,
        dew_scale_k=dew_scale_k,
        vapor_hexane_kg_s=tuple(vapor_hexane_kg_s),
        q_indirect_w_m3=tuple(q_indirect_w_m3),
        q_axial_w_m3=tuple(q_axial_w_m3),
        cell_volume_m3=cell_volume_m3,
        representative_particles=representative_particles,
        sorption_source=tuple(sorption_source),
        bottom_water_kg_s=bottom_water_kg_s,
        bottom_hexane_kg_s=bottom_hexane_kg_s,
        bottom_temperature_k=bottom_temperature_k,
        moisture_in=moisture_in,
        bound_in=moisture_in - free_in,
        free_in=free_in,
        initial_particle_tp=tuple(initial_particle.Tp),
        ref=c.vapor_enthalpy_ref,
        cp_water_liquid=pc.cp_water_liquid,
        antoine_water=c.antoine_water,
        luikov=c.luikov,
        pressure_pa=c.pressure_pa,
    )


def _coupled_block_energy_terms(
    u: np.ndarray, ctx: _CoupledBlockContext
) -> tuple[tuple[tuple[float, float, float, float], ...], float]:
    """Per-cell (accumulation, sorption, axial, surface) heat rates in W and the
    maximum paired interface mismatch, all evaluated at the accepted state."""
    terms: list[tuple[float, float, float, float]] = []
    interface_mismatch_max = 0.0
    for j in range(ctx.nz):
        substeps = ctx.time_substeps[j]
        last = substeps - 1
        x_here = u[ctx.xb_index(j)] + u[ctx.xf_index(j)]
        if j == 0:
            x_prev = ctx.moisture_in
            entering_tp = ctx.initial_particle_tp
        else:
            x_prev = u[ctx.xb_index(j - 1)] + u[ctx.xf_index(j - 1)]
            entering_tp = tuple(
                u[ctx.p_index(j - 1, ctx.time_substeps[j - 1] - 1, i)]
                for i in range(ctx.layers)
            )
        cv_new = ctx.base_cv + ctx.water_cv * x_here
        cv_old = ctx.base_cv + ctx.water_cv * x_prev
        accumulation = sum(
            (
                cv_new * (u[ctx.p_index(j, last, i)] - ctx.reference_temperature_k)
                - cv_old * (entering_tp[i] - ctx.reference_temperature_k)
            )
            * ctx.volumes[i]
            / ctx.cell_residence_s[j]
            for i in range(ctx.layers)
        )
        sorption = sum(
            ctx.sorption_source[j][i] * ctx.volumes[i] for i in range(ctx.layers)
        )
        axial = ctx.q_axial_w_m3[j] * ctx.total_volume
        surface = sum(
            ctx.surface_conductance_w_k
            * (u[ctx.p_index(j, s, ctx.layers - 1)] - u[ctx.tv_index(j)])
            for s in range(substeps)
        ) / substeps
        terms.append((accumulation, sorption, axial, float(surface)))
        rep = ctx.representative_particles[j]
        vapor_side = rep * ctx.surface_conductance_w_k * (
            sum(u[ctx.p_index(j, s, ctx.layers - 1)] for s in range(substeps))
            / substeps
            - u[ctx.tv_index(j)]
        )
        interface_mismatch_max = max(
            interface_mismatch_max, abs(vapor_side - rep * surface)
        )
    return tuple(terms), interface_mismatch_max


def _solve_phase_coupled_thermal_block(
    *,
    initial_particle: pt.ParticleState,
    particles: list[pt.ParticleState],
    vapor_temperature_guess_k: list[float],
    vapor_water_guess_kg_s: list[float],
    vapor_hexane_kg_s: list[float],
    bottom_water_kg_s: float,
    bottom_hexane_kg_s: float,
    bottom_temperature_k: float,
    bound_moisture_guess: list[float],
    free_moisture_guess: list[float],
    condensation_guess_kg_s: list[float],
    evaporation_guess_kg_s: list[float] | None = None,
    moisture_in: float,
    free_in: float,
    q_indirect_w_m3: tuple[float, ...],
    q_axial_w_m3: tuple[float, ...],
    cell_widths_m: tuple[float, ...],
    cell_residence_s: tuple[float, ...],
    time_substeps: tuple[int, ...],
    bed_area_m2: float,
    geometry: pt.ShellGeometry,
    c: DCZConstants,
    kappa_w: float,
    bound_capacity: float,
    m_dry_kg_s: float,
    max_iter: int = 150,
) -> _CoupledBlockResult:
    """Semismooth-Newton solve of the coupled phase/component thermal block.

    The unknowns, residuals, complementarity, and analytic Jacobian are
    documented in ``docs/DCZ_BLOCK_IMPLICIT_THERMAL_COUPLING_DESIGN_2026-07-26.md``.
    """
    nz = len(particles)
    if evaporation_guess_kg_s is None:
        evaporation_guess_kg_s = [0.0 for _ in range(nz)]
    if len(evaporation_guess_kg_s) != nz:
        raise ValueError("evaporation_guess_kg_s must have one value per cell")
    lambda_scale = max(bottom_water_kg_s, 1.0e-6)
    dew_scale_k = 10.0
    ctx = _build_coupled_block_context(
        initial_particle=initial_particle,
        particles=particles,
        vapor_hexane_kg_s=vapor_hexane_kg_s,
        bottom_water_kg_s=bottom_water_kg_s,
        bottom_hexane_kg_s=bottom_hexane_kg_s,
        bottom_temperature_k=bottom_temperature_k,
        moisture_in=moisture_in,
        free_in=free_in,
        q_indirect_w_m3=q_indirect_w_m3,
        q_axial_w_m3=q_axial_w_m3,
        cell_widths_m=cell_widths_m,
        cell_residence_s=cell_residence_s,
        time_substeps=time_substeps,
        bed_area_m2=bed_area_m2,
        geometry=geometry,
        c=c,
        kappa_w=kappa_w,
        bound_capacity=bound_capacity,
        m_dry_kg_s=m_dry_kg_s,
        lambda_scale=lambda_scale,
        dew_scale_k=dew_scale_k,
    )

    u = np.zeros(ctx.unknowns, dtype=float)
    for j in range(nz):
        for step in range(time_substeps[j]):
            for i in range(ctx.layers):
                u[ctx.p_index(j, step, i)] = particles[j].Tp[i]
        u[ctx.tv_index(j)] = vapor_temperature_guess_k[j]
        u[ctx.w_index(j)] = max(vapor_water_guess_kg_s[j], 1.0e-9)
        u[ctx.xb_index(j)] = bound_moisture_guess[j]
        u[ctx.xf_index(j)] = free_moisture_guess[j]
        u[ctx.lam_index(j)] = (
            max(condensation_guess_kg_s[j], 0.0)
            - max(evaporation_guess_kg_s[j], 0.0)
        )

    def scaled_merit(state: np.ndarray, weight: np.ndarray) -> float:
        residual = _coupled_block_residual(state, ctx) * weight
        return 0.5 * float(residual @ residual)

    converged = False
    iterations = 0
    for iterations in range(1, max_iter + 1):
        f = _coupled_block_residual(u, ctx)
        energy_max = float(np.max(np.abs(f[: ctx.particle_unknowns]))) if ctx.particle_unknowns else 0.0
        water_max = xb_max = xf_max = ncp_max = 0.0
        for j in range(nz):
            energy_max = max(energy_max, abs(f[ctx.tv_index(j)]))
            water_max = max(water_max, abs(f[ctx.w_index(j)]))
            xb_max = max(xb_max, abs(f[ctx.xb_index(j)]))
            xf_max = max(xf_max, abs(f[ctx.xf_index(j)]))
            ncp_max = max(ncp_max, abs(f[ctx.lam_index(j)]))
        if (
            energy_max < 1.0e-4
            and water_max < 1.0e-8
            and xb_max < 1.0e-9
            and xf_max < 1.0e-8
            and ncp_max < 1.0e-7
        ):
            converged = True
            break

        matrix = _coupled_block_jacobian(u, ctx)
        try:
            # SuperLU is faster than explicit Schur formation for the present
            # narrow-band ordering (0.060 s versus 0.231 s on the frozen
            # 20-cell/96-shell benchmark).  Keep the exact Schur routine as an
            # algebraic regression/reference, not as the production path.
            step = np.asarray(spsolve(matrix.tocsc(), -f), dtype=float)
        except Exception:  # pragma: no cover - singular Jacobian guard
            break
        if not np.all(np.isfinite(step)):
            break

        # Row-equilibrate the merit by each equation's Jacobian sensitivity.  The
        # raw first-law residuals are many orders larger than the dimensionless
        # complementarity residual, so an unscaled ``||F||^2`` line search
        # over-protects energy and starves the phase/moisture progress -- the
        # observed 1/128 step-starvation on stiff nested grids.  Scaling each
        # residual by its largest Jacobian entry makes every equation comparable
        # so a near-Newton step is admissible.
        row_scale = np.asarray(
            np.abs(matrix).max(axis=1).todense(), dtype=float
        ).ravel()
        weight = 1.0 / np.maximum(row_scale, 1.0e-30)

        # Fraction-to-boundary: keep every vapor-water flow strictly positive
        # and free-liquid inventory nonnegative so the thermodynamic closures
        # stay inside their domains.  The signed phase rate itself is unbounded
        # above and is bounded below implicitly by the free-water equation.
        alpha = 1.0
        for j in range(nz):
            widx = ctx.w_index(j)
            if step[widx] < 0.0:
                alpha = min(alpha, -0.9 * u[widx] / step[widx])
            xfidx = ctx.xf_index(j)
            if step[xfidx] < 0.0 and u[xfidx] > 0.0:
                alpha = min(alpha, -0.9 * u[xfidx] / step[xfidx])
        alpha = max(alpha, 1.0e-6)

        current = scaled_merit(u, weight)
        while alpha > 1.0e-8:
            trial = u + alpha * step
            if scaled_merit(trial, weight) <= (1.0 - 1.0e-4 * alpha) * current:
                break
            alpha *= 0.5
        u = u + alpha * step

    particles_out: list[pt.ParticleState] = []
    for j in range(nz):
        last = time_substeps[j] - 1
        temps = tuple(
            float(u[ctx.p_index(j, last, i)]) for i in range(ctx.layers)
        )
        if any(t < _TP_MIN_K or t > _TP_MAX_K for t in temps):
            raise ValueError(
                "DCZ coupled block left the qualified particle-temperature interval"
            )
        particles_out.append(pt.ParticleState(wpg2=particles[j].wpg2, Tp=temps))

    energy_terms, interface_mismatch = _coupled_block_energy_terms(u, ctx)
    final = _coupled_block_residual(u, ctx)
    residual_max_w = float(
        np.max(np.abs(final[: ctx.particle_unknowns]))
    ) if ctx.particle_unknowns else 0.0
    for j in range(nz):
        residual_max_w = max(residual_max_w, abs(final[ctx.tv_index(j)]))

    return _CoupledBlockResult(
        particles=tuple(particles_out),
        vapor_temperature_k=tuple(float(u[ctx.tv_index(j)]) for j in range(nz)),
        vapor_water_kg_s=tuple(float(u[ctx.w_index(j)]) for j in range(nz)),
        bound_moisture=tuple(float(u[ctx.xb_index(j)]) for j in range(nz)),
        free_moisture=tuple(float(u[ctx.xf_index(j)]) for j in range(nz)),
        condensed_water_kg_s=tuple(
            max(float(u[ctx.lam_index(j)]), 0.0) for j in range(nz)
        ),
        evaporated_water_kg_s=tuple(
            max(-float(u[ctx.lam_index(j)]), 0.0) for j in range(nz)
        ),
        particle_energy_terms_w=energy_terms,
        vapor_row_residual_w=tuple(
            float(final[ctx.tv_index(j)]) for j in range(nz)
        ),
        residual_max_w=residual_max_w,
        interface_mismatch_max_w=interface_mismatch,
        iterations=iterations,
        converged=converged,
    )
