"""Immutable contracts and acceptance validation for the integrated DT solver.

Keeping these structures independent of nonlinear orchestration makes the
scientific interface reviewable and lets runtime/UI code depend on stable
results without importing solver implementation details. The dt_solver module
re-exports every historical name for compatibility.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from dtdc_simulator.core import thermo
from dtdc_simulator.core.zones import dcz
from dtdc_simulator.core.zones import ftrz
from dtdc_simulator.core.zones import particle as pt
from dtdc_simulator.core.zones import phz as phz_mod


@dataclass(frozen=True)
class DTTray:
    """One real DT tray (PREDESOLV/MAIN/SPARGE), top-to-bottom solid order."""

    id: str
    role: str  # "PREDESOLV" | "MAIN" | "SPARGE"
    diameter_m: float
    bed_height_m: float
    Q_indirect_w: float
    direct_steam_kg_s: float = 0.0  # SPARGE only; must be the LAST tray if nonzero


@dataclass(frozen=True)
class SolidFeed:
    T: float  # K
    X1: float  # moisture, kg/kg dry solid
    X2: float  # hexane, kg/kg dry solid
    X3: float  # oil, kg/kg dry solid
    m_dry_kg_s: float


@dataclass(frozen=True)
class VaporFeed:
    """ "Clean" vapor arriving at the very bottom of the DT (from below the
    SPARGE tray, e.g. the rotary lock / DC side), BEFORE direct steam mixes
    in -- see module docstring on the sparge BC."""

    m_water_kg_s: float
    m_hex_kg_s: float
    T: float  # K


@dataclass(frozen=True)
class DTSolverConstants:
    phz: phz_mod.PHZConstants
    ftrz: ftrz.FTRZConstants  # carries T_boil_hexane/T_boil_water/bed_porosity/gab/oil/rho_ps/...
    particle: pt.ParticleConstants  # carries D_eff/r_P/Np/alpha_ps/alpha_pg/rho_ps/k_ps/k_pg/...
    D_ax: float  # m2/s, axial dispersion (DCZ hexane balance, eq. A.36)
    k_mixL: float  # W/(m K), bed-scale solid/gas mixture conductivity (eq. A.32)
    rho_V: float  # kg/m3, vapor reference density
    cp_V: float  # J/(kg K), vapor specific heat
    mu_V: float  # Pa.s, vapor viscosity
    D_HW: float  # m2/s, hexane-water diffusivity (Schmidt number, eq. B.10)
    direct_steam_pressure_pa: float  # Pa absolute, measured post-valve steam boundary
    T_direct_steam: float  # K, measured post-valve steam temperature/enthalpy coordinate
    luikov: thermo.LuikovParams  # water sorption/desorption isotherm, DCZ's subsaturated regime
    water_diffusivity: float  # m2/s, water's own intraparticle diffusivity (Touffet et al. 2026)
    # -- NOT hM (hexane-vapor-tuned, ~25-100x too fast once coupled to an energy balance, see
    # zones/dcz.py's own docstring for the confirmed inversion and why diffusivity, not a
    # convective coefficient, is the right basis)
    antoine_hexane: thermo.AntoineParams  # hexane saturation pressure -- for the DT dome/vapor
    # temperature (binary hexane-water dew point of the exiting vapor, see the axial-profile build)
    pressure_pa: float = thermo.ATM_PRESSURE_PA  # lower-DT internal operating pressure for the
    # water dew-point / activity calc (FTRZ + DCZ). Above atmospheric per the sparge-tray
    # pressure drop (Kemper 2019, 0.35-0.70 kg/cm2). The DOME stays at atmospheric (its own
    # binary dew point, `_binary_dew_T`, is unchanged). Default atmospheric.

    @property
    def k_V(self) -> float:
        """Bulk interstitial vapor thermal conductivity -- reused from the
        particle pore-gas conductivity (same gas, different scale); see
        module docstring."""
        return self.particle.k_pg

@dataclass(frozen=True)
class PHZPassResult:
    tray_results: tuple[phz_mod.PHZTrayResult, ...]  # fully-PHZ trays, top to bottom
    boundary_tray_index: int  # index into the caller's `trays` list
    z_star_m: float  # PHZ's own share of the boundary tray's height
    boundary_tray_result: phz_mod.PHZTrayResult  # the boundary tray's own PHZ portion
    # (height z_star_m), kept for the axial profile (DTAxialProfile) -- everything else here
    # already only needed the boundary tray's own EXIT state (`exit_state` below), but the
    # profile needs its per-cell interior too.
    # At T_boil_hexane; X2 may exceed X2,cr at the PREDESOLV hardware boundary.
    exit_state: phz_mod.SolidState
    L_PHZ_m: float

@dataclass(frozen=True)
class TraySummary:
    id: str
    T: float  # K, exit (bottom-face) solid temperature
    X1: float  # kg/kg dry solid, exit moisture
    X2: float  # kg/kg dry solid, exit hexane


@dataclass(frozen=True)
class DTAxialProfile:
    """Per-cell axial profile spanning the WHOLE DT (PHZ -> FTRZ -> DCZ, top to
    bottom), for visualization (the HMI's own "profile along the tower", not
    consumed by `Model.step` -- that only needs `tray_summaries`' per-tray exit
    values). Parallel tuples, one entry per cell, ordered top-to-bottom by
    `z_m` (cumulative distance from the DT's own top face).

    `vapor_flow_kg_s` is a real, cell-varying quantity throughout: PHZ/FTRZ
    evaporation/condensation and DCZ particle transfer/water sorption all
    change the explicit component flows as vapor travels upward.
    """

    z_m: tuple[float, ...]
    zone: tuple[str, ...]  # "PHZ" | "FTRZ" | "DCZ"
    stage_id: tuple[str, ...]  # which real tray this cell belongs to
    solid_T: tuple[float, ...]  # K
    solid_X1: tuple[float, ...]  # kg/kg dry solid, moisture
    solid_X2: tuple[float, ...]  # kg/kg dry solid, hexane
    vapor_T: tuple[float, ...]
    vapor_flow_kg_s: tuple[float, ...]  # total (water + hexane)
    vapor_hexane_frac: tuple[float, ...]  # mass fraction, 0-1
    vapor_water_frac: tuple[float, ...]  # mass fraction, 0-1
    mechanism: tuple[str, ...]  # local transfer closure, finer than hardware zone


@dataclass(frozen=True)
class DTCouplingResiduals:
    solid_interface_T: float = math.inf
    vapor_interface_T: float = math.inf
    vapor_interface_hexane_fraction: float = math.inf
    vapor_interface_water_flow: float = math.inf
    vapor_interface_hexane_flow: float = math.inf
    ftrz_length: float = math.inf

    def maximum_scaled(self, outer_tol: float) -> float:
        """Dimensionless shadow metric using the current physical tolerances."""
        length_tol = min(outer_tol, 1.0e-3)
        return max(
            self.solid_interface_T / outer_tol,
            self.vapor_interface_T / outer_tol,
            self.vapor_interface_hexane_fraction / outer_tol,
            self.vapor_interface_water_flow / outer_tol,
            self.vapor_interface_hexane_flow / outer_tol,
            self.ftrz_length / length_tol,
        )


@dataclass(frozen=True)
class DCZTransientSeed:
    """Resolved steady data required for conservative transient initialization.

    This is deliberately richer than the visualization profile: it preserves
    the actual radial particle states, physical cell widths, phase split,
    component flows, and per-cell indirect duty used by the accepted DCZ
    solve.
    """

    constants: dcz.DCZConstants
    dry_solid_flow_kg_s: float
    solid_feed_particle: pt.ParticleState
    solid_feed_bound_water_kg_kg_dry: float
    solid_feed_free_water_kg_kg_dry: float
    bottom_vapor_water_kg_s: float
    bottom_vapor_hexane_kg_s: float
    bottom_vapor_temperature_k: float
    indirect_duty_by_cell_w: tuple[float, ...]


@dataclass(frozen=True)
class DTResult:
    phz: PHZPassResult
    ftrz: ftrz.FTRZZoneResult
    dcz: dcz.DCZZoneResult
    tray_summaries: tuple[TraySummary, ...]
    axial_profile: DTAxialProfile
    L_PHZ_m: float
    L_FTRZ_m: float
    L_DCZ_m: float
    hQ: float
    hM: float
    aV: float
    outer_iterations: int
    converged: bool
    coupling_residuals: DTCouplingResiduals = DTCouplingResiduals()
    transient_seed: DCZTransientSeed | None = None

    @property
    def solid_exit_X2(self) -> float:
        """Hexane leaving the whole DT (kg/kg dry solid) -- the DCZ's own
        volumetric-mean exit value, Coletto's own KPI (Fig. 9(a))."""
        return self.dcz.solid_out_X2


def validate_dt_result(result: DTResult, solid_feed: SolidFeed, c: DTSolverConstants) -> None:
    """Raise ``ValueError`` unless a result is safe to publish.

    ``solve_dt`` intentionally returns its best iterate when it reaches the
    outer-iteration cap.  That is useful diagnostically, but a digital-twin
    acceptance boundary must not replace its last trusted state with either a
    nonconverged iterate or a numerically finite yet physically impossible
    profile.

    The lower temperature bound expresses the no-refrigeration character of
    the DT: jackets, condensing steam, and hot vapor cannot cool the bulk meal
    below the colder of the incoming meal and hexane's boiling temperature.
    In particular, a water dew-point closure may describe an interface but
    cannot overwrite the FTRZ bulk temperature governed by Coletto A.17.
    """
    if not result.converged:
        raise ValueError("DT solve did not converge")
    if not result.dcz.converged:
        raise ValueError("DT solve contains a nonconverged DCZ inner solution")
    if result.dcz.coupling_tolerance is None and result.dcz.residuals.maximum_scaled > 1.0:
        raise ValueError("DT solve contains unresolved DCZ fixed-point residuals")

    profile = result.axial_profile
    parallel_fields = (
        profile.z_m,
        profile.zone,
        profile.stage_id,
        profile.solid_T,
        profile.solid_X1,
        profile.solid_X2,
        profile.vapor_T,
        profile.vapor_flow_kg_s,
        profile.vapor_hexane_frac,
        profile.vapor_water_frac,
        profile.mechanism,
    )
    n = len(profile.z_m)
    if n == 0 or any(len(values) != n for values in parallel_fields):
        raise ValueError("DT axial profile fields must be nonempty and parallel")
    if any(not math.isfinite(value) for value in profile.z_m):
        raise ValueError("DT axial coordinates must be finite")
    if any(b < a for a, b in zip(profile.z_m, profile.z_m[1:])):
        raise ValueError("DT axial coordinates must be ordered top-to-bottom")

    numeric_profile_fields = (
        profile.solid_T,
        profile.solid_X1,
        profile.solid_X2,
        profile.vapor_T,
        profile.vapor_flow_kg_s,
        profile.vapor_hexane_frac,
        profile.vapor_water_frac,
    )
    if any(not math.isfinite(value) for values in numeric_profile_fields for value in values):
        raise ValueError("DT axial profile contains a non-finite state")

    min_bulk_T = min(solid_feed.T, c.ftrz.T_boil_hexane) - 1.0
    if any(value < min_bulk_T for value in profile.solid_T):
        raise ValueError("DT bulk-meal temperature fell below its physical lower bound")
    if any(value <= 0.0 for value in profile.vapor_T):
        raise ValueError("DT vapor temperature must remain above absolute zero")
    if any(value <= 0.0 for value in profile.vapor_flow_kg_s):
        raise ValueError("DT vapor flow must remain positive")
    for values, name in (
        (profile.solid_X1, "solid moisture"),
        (profile.solid_X2, "solid hexane"),
        (profile.vapor_hexane_frac, "vapor hexane fraction"),
        (profile.vapor_water_frac, "vapor water fraction"),
    ):
        if any(value < -1.0e-9 or value > 1.0 + 1.0e-9 for value in values):
            raise ValueError(f"DT {name} left the physical [0, 1] range")
    if any(
        abs(hexane + water - 1.0) > 1.0e-6
        for hexane, water in zip(profile.vapor_hexane_frac, profile.vapor_water_frac)
    ):
        raise ValueError("DT binary-vapor mass fractions do not sum to one")

    if not result.tray_summaries:
        raise ValueError("DT result must contain tray summaries")
    for summary in result.tray_summaries:
        if not all(math.isfinite(value) for value in (summary.T, summary.X1, summary.X2)):
            raise ValueError(f"DT tray {summary.id} contains a non-finite state")
        if summary.T < min_bulk_T:
            raise ValueError(f"DT tray {summary.id} temperature is below its physical lower bound")
        if not (0.0 <= summary.X1 <= 1.0 and 0.0 <= summary.X2 <= 1.0):
            raise ValueError(f"DT tray {summary.id} composition left the physical [0, 1] range")

    lengths = (result.L_PHZ_m, result.L_FTRZ_m, result.L_DCZ_m)
    if any(not math.isfinite(value) or value < 0.0 for value in lengths):
        raise ValueError("DT zone lengths must be finite and non-negative")
    if not all(math.isfinite(value) and value > 0.0 for value in (result.hQ, result.hM, result.aV)):
        raise ValueError("DT transport coefficients must be finite and positive")
    if abs(result.L_FTRZ_m - result.ftrz.L_FTRZ_m) > 1.0e-3:
        raise ValueError("DT FTRZ result and moving-boundary geometry are inconsistent")
    if (
        result.dcz.warm_start is None
        or abs(result.L_DCZ_m - result.dcz.warm_start.bed_height_m) > 1.0e-9
    ):
        raise ValueError("DT DCZ result and remeshed geometry are inconsistent")
