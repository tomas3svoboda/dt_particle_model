"""Bed scale of the Diffusion-Controlled Zone (DCZ) -- Coletto, Bandoni &
Blanco (2022), Table A.4 / §2.4.2 / §A.3.2-A.3.3, coupled to the particle
scale (`core/zones/particle.py`) via the REAL algorithm published in the
paper's supplementary material (Fig. 3, "Primary Internal Loop") -- not a
reconstruction. M2 Phase 3 (BuildSpec §14): standalone, pure, unit-tested --
not yet wired into `core/model.py`'s tray-by-tray sweep (M2 Phase 4).

THE PRIMARY INTERNAL LOOP (exactly Fig. 3 of the supplementary material),
iteration index `s`, across all `nz` axial cells at once each step:
  1. Energy balance at particle scale -- march every cell's particle
     temperature field top-to-bottom, using the *current* vapor temperature
     profile and the *previous* iteration's mass-march rate (sorption heat).
  2. Energy balance at bed scale -- integrate bottom-to-top using the
     just-updated particle surface temperatures, producing a new TV profile.
  3. Mass balance at particle scale -- march every cell's particle hexane
     field top-to-bottom, using the just-updated TV profile (isotherm is
     T-dependent) and the *previous* iteration's wV2 profile.
  4. Mass balance at bed scale -- integrate bottom-to-top using the
     particle's conserved hexane loss and the same-pass water transfer,
     producing explicit component-flow and derived-composition profiles.
  5. Convergence check on every exported boundary plus maximum full vapor
     temperature/component-flow profile changes. Not converged -> repeat from
     the complete lagged state.

DOCUMENTED RESOLUTION -- bed-scale marching direction and sign (not a
literal transcription of eq. A.35): the source material states `alphaV *
||uV|| * rho * dphi/dz = -kappa_phi*aV*alphaL*(phiV-phiL) + S*_phiV`, but
doesn't fully spell out (in what's available to us) the vector orientation
of `r-hat` or how the vapor's own upward (decreasing-z) velocity enters the
z-derivative's sign. Rather than risk a silent sign error from an
ambiguous convention, this module marches bottom-to-top cell-by-cell adding
the PHYSICALLY-motivated transfer (vapor gains hexane/heat from the
particle as it flows from a cell's bottom face to its top face) -- the same
category of documented gap-closure as FTRZ's eq. A.18 units resolution.

Axial conduction (energy, eq. A.32) remains a lagged conservative
finite-volume Laplacian. Axial dispersion (mass, eq. A.36) is instead solved
implicitly with the conservative counter-current component balance. The
following historical source note applies to the conduction correction:
§3.3 remark that `q_condL` "is calculated using the temperature difference
between two cells") -- lagged by one outer pass, the same technique used for
`L_FTRZ`'s free boundary.

BED-SCALE MARCHING IS IMPLICIT PER CELL, not the naive explicit (forward)
step eq. A.35 might otherwise suggest: the particle<->vapor transfer
coefficient (`hM`/`hQ * aV * alphaL`) is typically stiff relative to a
practical cell size -- its own relaxation length can be much shorter than
`dz` -- so an explicit step diverges (confirmed by hitting exactly this
during development: cell-to-cell values blew up by many orders of
magnitude). A per-cell backward-Euler-style relaxation step is
unconditionally stable for this linear form, mirroring why the particle
-scale march is implicit for the same underlying stiffness reason.

ENERGY-BALANCE FIX (M2 Phase 4, `core/dt_solver.py`'s integration work):
M2 Phase 3 shipped this module WITHOUT eq. A.34/A.37's `SVm2*Ĥ2`/`ṁ'ax,net*
Ĥ2` enthalpy-transport terms, on the documented grounds that implementing
them literally (as an absolute mass-flux-times-enthalpy quantity) caused
runaway heating. Integrating this module into the full DT solve (realistic
`Q_indirect` magnitudes, not the small illustrative ones this module's own
tests use) surfaced that dropping them isn't the bounded "runs a bit cooler"
approximation it was believed to be: confirmed by direct instrumentation,
the particle<->vapor system has no energy floor without them and drifts to
UNBOUNDED cooling over enough outer iterations, eventually crashing the GAB
isotherm's own validity range. Fixed by restoring the missing term, but as a
different (and correct) quantity: not the bed-scale *surface* mass flux
`SVm2` (tried first, negligible effect -- a particle's own diffusive
relaxation is far slower than one axial cell's residence time, so the
surface flux badly lags the true internal desorption rate), but the
particle-VOLUME-integrated sorption/desorption sink (exactly what step 1
already subtracts from the particle via eq. A.30's first two terms),
credited back to the vapor at the SAME cell with the opposite sign -- see
`march_particle_energy`'s step-1 call site below and `zones/particle.py::
sorption_heat_sink_volumetric_mean_w_m3`'s own docstring for the full
derivation. Being an EXACT transfer between the two already-computed energy
balances (not an independently-derived absolute term), it cannot manufacture
a runaway the way the first attempt did.

MOISTURE (H2O) BALANCE (found this session, not part of M2 Phase 3/4's own
scope): this module originally carried NO water balance at all -- the
`direct_steam` MV (SP1's sparge injection, mixed into the DT's bottom vapor
BC in `core/dt_solver.py`) flowed through DCZ's vapor stream without ever
being able to condense onto the solid, so changing it had literally no
effect on any tray's reported moisture. Coletto's own DCZ sub-model is
hexane-only (§7.6). A first version (M2-follow-up "H2O balance") added
bulk-vapor dew-point condensation only (mirroring `zones/ftrz.py`'s own
V-SCAL/V-SAT switch) -- correct as far as it went, but a binary
supersaturated-or-not threshold can't respond to the DEGREE of humidity,
only to crossing a cutoff, which produced a real, confirmed bug: `direct_
steam` (injected hotter than the surrounding vapor) moved the boundary
FURTHER from that cutoff as it increased, so more steam meant LESS
condensation -- backwards from real plants. Fixed by adding a genuine
sorption/desorption isotherm for the (much more common) subsaturated
regime, once a real water-sorption isotherm for soybean meal was found
(`literature_sources/Gianini_Study_of_the_equilibrium_isotherms_of_soybean_
meal.pdf` -- measured on meal sampled directly from a desolventizer/
toaster's own outlet; `thermo.LuikovParams`/`thermo.luikov_equilibrium_
moisture`) -- structurally the SAME architecture hexane's own GAB isotherm
already uses (`wpg2` -> `W2(a_h,T)`), just lumped/lagged rather than
radially resolved (no diffusivity data exists for water in this matrix to
justify a 12-layer FVM the way hexane gets one).

Two additive mechanisms use the same water-activity state. Supersaturated
vapor is pinned to its dew point and the required bulk condensate is
back-solved from the cell energy balance. The resulting wet solid then makes
a finite-rate bidirectional adjustment toward `Xe(a_w)` over its local
residence. The condensation active-set mass is relaxed between iterations;
an actively condensing saturated cell cannot evaporate its new free
condensate against the clamped bound-water isotherm.

Both regimes accumulate top-to-bottom (solid flow order) into each cell's
own `X1_bulk`, starting from what FTRZ handed off (`X1_in`) -- replacing
what `core/dt_solver.py` used to carry forward from FTRZ completely
unchanged through every DCZ-spanned tray.

EXTRAPOLATION CAVEAT, stated not hidden: the Gianini isotherm's own tested
range is 15-70 C; DCZ's own operating temperatures currently reach higher.
The cited paper's own finding that temperature barely affects the isotherm
in ITS tested range is reassuring but doesn't cover that gap.

VARIABLE VAPOR FLOW (2026-07-23): water and hexane are marched as explicit
component mass flows. Condensation/sorption therefore reduce water kg/s and
particle transfer increases hexane kg/s cell by cell; local vapor heat
capacity/velocity use the resulting total flow. `hQ`/`hM`/`aV` remain the
once-computed bed coefficients supplied by the caller, but the former
fixed-total-flow composition approximation is removed.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from dtdc_simulator.core import thermo
from dtdc_simulator.core.zones import particle as pt
from dtdc_simulator.core.zones import particle_jit
from dtdc_simulator.core.zones.dcz_block import (
    _CoupledBlockContext as _CoupledBlockContext,
    _CoupledBlockResult as _CoupledBlockResult,
    _ThermalBlockResult as _ThermalBlockResult,
    _build_coupled_block_context as _build_coupled_block_context,
    _coupled_block_energy_terms as _coupled_block_energy_terms,
    _coupled_block_jacobian as _coupled_block_jacobian,
    _coupled_block_residual as _coupled_block_residual,
    _coupled_cell_derived as _coupled_cell_derived,
    _dew_point_temperature_and_dwater as _dew_point_temperature_and_dwater,
    _luikov_xe_and_slope as _luikov_xe_and_slope,
    _solve_coupled_newton_step_schur as _solve_coupled_newton_step_schur,
    _solve_phase_coupled_thermal_block,
    _solve_phase_fixed_thermal_block,
    _water_activity_and_grads as _water_activity_and_grads,
)
from dtdc_simulator.core.zones.dcz_enthalpy import (
    evaporation_required_for_wet_temperature_constraint,
    solve_vapor_cell_enthalpy,
    solve_vapor_cell_temperature_at_state,
    vapor_cell_enthalpy_residual_at_state,
)
from dtdc_simulator.core.zones.dcz_grid import (
    axial_laplacian,
    resample_warm_start as _resample_warm_start,
    solve_vapor_hexane_dispersion as _solve_vapor_hexane_dispersion,
)
from dtdc_simulator.core.zones.dcz_types import (
    BALANCE_TOL_FLOW,
    ENERGY_BALANCE_TOL_W,
    EXIT_TOL_FLOW,
    EXIT_TOL_T,
    EXIT_TOL_X1,
    EXIT_TOL_X2,
    PHASE_EQUILIBRIUM_TOL_K,
    PROFILE_TOL_FLOW,
    PROFILE_TOL_T,
    DCZCellResult,
    DCZConstants,
    DCZResiduals,
    DCZWarmStart,
    DCZZoneResult,
    VaporState,
)
from dtdc_simulator.core.zones.dcz_water import (
    countercurrent_vapor_water_profile as _countercurrent_vapor_water_profile,
    exact_linear_driving_force_step,
    solve_bound_water_countercurrent_bvp,
)

# Physical bounds for the particle temperature (K), a robustness clamp on the coupled
# water-sorption <-> energy march (see solve_dcz_zone step 1). The DT operates entirely
# within ~ -20..205 C; at every realistic operating point the particle sits far inside
# this band, so the clamp only engages for a deeply off-design iterate (e.g. an under-set
# sparge far below the ~110 kg/t design rate) where the lagged latent-heat feedback would
# otherwise run the temperature away toward 0 K / overflow. Bounding keeps the iterate
# finite so it degrades to a converged-but-imperfect result instead of crashing.
_TP_MIN_K = 250.0
_TP_MAX_K = 480.0


def _dominant_mode_extrapolation(
    x_k: "np.ndarray", f_k: "np.ndarray", f_prev: "np.ndarray", boost_max: float
) -> "np.ndarray | None":
    """Geometric (single-dominant-mode) extrapolation of a fixed-point iterate.

    After a Picard iteration's fast modes decay, its residual ``f = H(x) − x`` is
    dominated by ONE slow eigenmode with eigenvalue λ near 1 (the DCZ energy
    coupling, ρ≈0.9998), so ``f_k ≈ λ f_{k-1}``. We estimate λ from the residual
    alignment (a Rayleigh quotient) and take the geometric-series step toward the
    fixed point, ``x* = x_k + f_k / (1 − λ)``, with the boost ``1/(1−λ)`` capped
    at `boost_max`. Returns the extrapolated iterate, or ``None`` if the residual
    does not look like a clean slow mode (then the caller takes a plain step).

    This deliberately replaces multi-vector Anderson: the DCZ map is piecewise
    (water condensation/isotherm branches switch as cells cross their dew point),
    which makes Anderson's least-squares over residual DIFFERENCES ill-behaved. A
    single scalar extrapolation along the current residual direction is far more
    robust on such a map while still collapsing the O(10^4)-iteration crawl."""
    denom = float(f_prev @ f_prev)
    if denom <= 0.0:
        return None
    lam = float(f_k @ f_prev) / denom
    if not (0.5 < lam < 0.9999):
        return None  # not a clean, contracting slow mode -> plain step
    boost = min(1.0 / (1.0 - lam), boost_max)
    x_next = x_k + f_k * boost
    if not np.all(np.isfinite(x_next)):
        return None
    return x_next


@dataclass
class _AdaptiveRelaxation:
    """Residual-monotonic scalar damping for one coupled variable family."""

    value: float
    minimum: float
    maximum: float
    previous_norm: float = math.inf

    def observe(self, norm: float, active_set_changed: bool = False) -> None:
        if not math.isfinite(norm):
            self.value = max(self.minimum, 0.5 * self.value)
        elif active_set_changed or norm > 1.05 * self.previous_norm:
            self.value = max(self.minimum, 0.5 * self.value)
        elif norm < 0.85 * self.previous_norm:
            self.value = min(self.maximum, 1.08 * self.value)
        self.previous_norm = norm


def bulk_temperature(cell: DCZCellResult, geometry: pt.ShellGeometry) -> float:
    """Volumetric-mean particle temperature (eq. 8) for a cell -- the bulk
    solid temperature a caller would report/compare against Fig. 9."""
    return pt.shell_volumetric_mean(cell.particle.Tp, geometry)


def solve_dcz_zone(
    nz: int,
    m_dry_kg_s: float,
    m_vapor_kg_s: float,
    T_L_sup: float,
    vapor_inf: VaporState,
    q_Iv_w_m3: float | tuple[float, ...],
    c: DCZConstants,
    X1_in: float = 0.0,
    outer_max_iter: int = 100,
    outer_tol: float = 1.0e-5,
    outer_relaxation: float = 0.5,
    residual_log: list[tuple[int, float, float, float, float, float]] | None = None,
    warm_start: DCZWarmStart | None = None,
    adaptive_damping: bool = True,
    particle_max_step_s: float | None = None,
    coupling_tolerance: float | None = None,
    axial_cell_fractions: tuple[float, ...] | None = None,
    thermal_coupling: str = "partitioned",
    progress_callback: Callable[[dict[str, object]], None] | None = None,
) -> DCZZoneResult:
    """Solve the DCZ, `nz` axial cells top (FTRZ handoff) to bottom, via the
    Primary Internal Loop above. The particle's own initial condition at
    zone entry is `wpg2=1.0` uniformly (pores saturated with hexane vapor,
    `a_h=1` -- exactly FTRZ's own termination condition, `X2=X2,eq(TV,inf)`,
    which by eq. 5's definition IS the `a_h=1` state) and `Tp=T_L_sup`
    uniformly (FTRZ's exit solid temperature). `X1_in` is the solid moisture
    (kg/kg dry solid) entering DCZ from FTRZ above -- see module docstring's
    moisture-balance section for how it's accumulated forward through DCZ.

    `outer_relaxation` under-relaxes the vapor profile updates (steps 2 and
    4) between passes -- reusing `ModelParams.outer_relaxation`'s existing
    convention from the tray-sweep design (§7.9). Needed in practice: without
    it, this Gauss-Seidel coupling was observed during development to drift
    for hundreds of iterations before settling (the vapor and particle
    profiles overshooting and correcting each other pass-over-pass) rather
    than converging cleanly.

    `q_Iv_w_m3` may be a single scalar (uniform volumetric indirect heat
    across the whole zone, e.g. for standalone/illustrative use -- the
    original M2 Phase 3 convention, unchanged) or a length-`nz` per-cell
    profile (top-to-bottom, matching this function's own `j` cell indexing
    below): DCZ commonly spans several real trays with materially different
    `Q_indirect` (M2 Phase 4's `core/dt_solver.py` builds this profile from
    each tray's own duty divided by its own volume, rather than smearing all
    spanned trays' duties into one artificial average).

    ``particle_max_step_s`` is an independent accuracy control for the
    particle energy and hexane marches.  When set, each axial residence step
    is divided into ``ceil(dt / particle_max_step_s)`` equal backward-Euler
    substeps.  The default ``None`` retains one particle step per axial cell.

    ``thermal_coupling="block_implicit"`` solves particle temperature, vapor
    temperature and water flow, bound/free solid moisture, and a signed net
    condensation/evaporation rate as one sparse semismooth system.  The
    default ``"partitioned"`` retains the published primary internal loop.
    """
    if thermal_coupling not in {"partitioned", "block_implicit"}:
        raise ValueError(
            "thermal_coupling must be 'partitioned' or 'block_implicit'"
        )
    A_bed = math.pi / 4.0 * c.diameter_m**2
    if axial_cell_fractions is None:
        cell_fractions = tuple(1.0 / nz for _ in range(nz))
    else:
        if len(axial_cell_fractions) != nz or any(
            not math.isfinite(value) or value <= 0.0 for value in axial_cell_fractions
        ):
            raise ValueError("axial_cell_fractions must be finite, positive, and length nz")
        fraction_sum = sum(axial_cell_fractions)
        if not math.isclose(fraction_sum, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("axial_cell_fractions must sum to one")
        cell_fractions = tuple(value / fraction_sum for value in axial_cell_fractions)
    cell_widths_m = tuple(c.bed_height_m * value for value in cell_fractions)
    q_Iv_profile = (
        q_Iv_w_m3 if isinstance(q_Iv_w_m3, tuple) else tuple(q_Iv_w_m3 for _ in range(nz))
    )
    if len(q_Iv_profile) != nz:
        raise ValueError(f"q_Iv_w_m3 profile length ({len(q_Iv_profile)}) must equal nz ({nz})")
    u_L = m_dry_kg_s / (c.particle.alpha_ps * c.alpha_L * c.particle.rho_ps * A_bed)
    cell_residence_s = tuple(width / u_L for width in cell_widths_m)
    if particle_max_step_s is not None and (
        not math.isfinite(particle_max_step_s) or particle_max_step_s <= 0.0
    ):
        raise ValueError("particle_max_step_s must be finite and positive when provided")
    particle_time_substeps = tuple(
        1 if particle_max_step_s is None else max(1, math.ceil(dt / particle_max_step_s))
        for dt in cell_residence_s
    )

    geometry = pt.build_shell_geometry(c.particle.r_P, c.particle.Np)
    zero_rates = tuple(0.0 for _ in range(c.particle.Np))

    initial_particle = pt.ParticleState(
        wpg2=tuple(1.0 for _ in range(c.particle.Np)),
        Tp=tuple(T_L_sup for _ in range(c.particle.Np)),
    )
    cascade_invariants = particle_jit.build_invariants(initial_particle, geometry)
    m_water_bottom_kg_s = (1.0 - vapor_inf.wV2) * m_vapor_kg_s
    m_hex_bottom_kg_s = vapor_inf.wV2 * m_vapor_kg_s
    saturated_bound_water_capacity = thermo.luikov_equilibrium_moisture(
        thermo.LUIKOV_MAX_VALIDATED_UR,
        c.luikov,
    )
    if warm_start is None:
        particles: list[pt.ParticleState] = [initial_particle for _ in range(nz)]
        vapor_m_water = [m_water_bottom_kg_s for _ in range(nz)]
        vapor_m_hex = [m_hex_bottom_kg_s for _ in range(nz)]
        vapor_T = [vapor_inf.T for _ in range(nz)]
        dwpg2_dt_prev: list[tuple[float, ...]] = [zero_rates for _ in range(nz)]
        X1_profile = [X1_in for _ in range(nz)]
        initial_free_water = max(X1_in - saturated_bound_water_capacity, 0.0)
        free_water_profile = [initial_free_water for _ in range(nz)]
        bound_water_mass_rate_w_m3 = [0.0 for _ in range(nz)]
        condensed_profile_kg_s = [0.0 for _ in range(nz)]
        evaporated_profile_kg_s = [0.0 for _ in range(nz)]
    else:
        warm = _resample_warm_start(
            warm_start,
            nz,
            c.bed_height_m,
            c.particle.Np,
            cell_widths_m,
        )
        particles = list(warm.particles)
        vapor_m_water = [max(value, 0.0) for value in warm.vapor_water_kg_s]
        vapor_m_hex = [max(value, 0.0) for value in warm.vapor_hexane_kg_s]
        vapor_T = list(warm.vapor_T)
        dwpg2_dt_prev = list(warm.dwpg2_dt)
        X1_profile = list(warm.X1)
        free_water_profile = [
            min(max(value, 0.0), max(X1_profile[j], 0.0))
            for j, value in enumerate(warm.free_water)
        ]
        bound_water_mass_rate_w_m3 = list(warm.bound_water_mass_rate_w_m3)
        condensed_profile_kg_s = [max(value, 0.0) for value in warm.condensed_water_kg_s]
        evaporated_profile_kg_s = [max(value, 0.0) for value in warm.evaporated_water_kg_s]
    vapor_wV2 = [
        vapor_m_hex[j] / max(vapor_m_water[j] + vapor_m_hex[j], 1.0e-12) for j in range(nz)
    ]
    kappa_w = 15.0 * c.water_diffusivity / c.particle.r_P**2  # 1/s, water's own bed-scale
    # equilibration rate -- the Glueckauf linear-driving-force (LDF) approximation for a
    # diffusing sphere (a standard, well-established result, not invented for this project).
    # Found this session (DECISIONS.md): reusing hexane's own hM*aV here (~0.05-0.17 m/s * ~1800
    # m2/m3, a Faner-correlation value tuned for hexane VAPOR transport) was ~25-100x too fast
    # once this relaxation was genuinely coupled to the vapor enthalpy balance
    # below -- confirmed directly: it inverted the basic "more indirect duty -> hotter DT exit"
    # relationship (doubled duty gave a COOLER converged profile than halved duty). The water
    # -specific CONVECTIVE coefficient tried next (literature_sources/Touffet_Moisture_sorption_
    # and_diffusion_in_pellet_animal_feed.pdf, 2026, dynamic vapor sorption on pelleted animal
    # feed) turned out to be insufficient too -- their own Biot-number analysis explains why:
    # INTERNAL particle diffusion, not external convection, is the actual rate-limiting step
    # (confirmed directly: even their measured k, combined with this project's own bed geometry,
    # still implied near-instant equilibration, ~500x faster than DCZ's own per-cell residence).
    # Using their own diffusion coefficient instead (D, not k) via the LDF form above lands in a
    # genuinely meaningful, neither-instant-nor-negligible band relative to DCZ's own residence
    # time -- confirmed this session, see DECISIONS.md for the measured comparison. NOT combined
    # with the convective coefficient in a full two-resistance series model (internal diffusion
    # dominates enough here that the convective term's own contribution is negligible by
    # comparison) -- a leaner, still real-data-grounded fix.
    # Free-liquid saturation temperature at local pressure. This constraint is
    # never applied to bound moisture; bound-water activity uses the actual
    # particle-surface temperature.
    pure_water_saturation_temperature = thermo.dew_point_temperature(
        0.0, c.antoine_water, P=c.pressure_pa
    )
    if thermo.water_phase_change_enthalpy(
        pure_water_saturation_temperature,
        c.vapor_enthalpy_ref,
        c.particle.cp_water_liquid,
    ) <= 0.0:
        raise ValueError("water vapor enthalpy must exceed liquid-water enthalpy in the DCZ")
    # Bound-water mass transfer remains one Picard pass lagged in the vapor
    # enthalpy equation and is marched explicitly after step 4.5.

    iterations = 0
    # Previous-pass exit KPIs (bottom cell) for the physical convergence test;
    # seeded at inf so the first pass can never trip.
    prev_x2_exit = math.inf
    prev_T_exit = math.inf
    prev_x1_exit = math.inf
    prev_water_out = math.inf
    prev_hex_out = math.inf
    prev_vapor_top_T = math.inf
    residuals = DCZResiduals()
    converged = False

    initial_alpha = min(max(outer_relaxation, 1.0e-3), 1.0)
    mean_cell_width_m = c.bed_height_m / nz
    minimum_width_ratio = min(cell_widths_m) / mean_cell_width_m
    # Axial solid conduction remains a Picard-lagged conservative Laplacian.
    # Its explicit fixed-point stability limit scales with the square of the
    # smallest cell width. Hexane axial dispersion is solved implicitly with
    # the component balance below and therefore must not inherit this cap.
    diffusive_damping_scale = min(
        1.0,
        max(1.0e-3, minimum_width_ratio**2),
    )
    transport_alpha_cap = initial_alpha * diffusive_damping_scale
    transport_alpha_floor = min(0.05, max(1.0e-4, 0.1 * transport_alpha_cap))
    temperature_alpha = min(
        warm.temperature_relaxation if warm_start is not None else initial_alpha,
        transport_alpha_cap,
    )
    hexane_alpha = min(
        warm.hexane_relaxation if warm_start is not None else initial_alpha,
        initial_alpha,
    )
    water_alpha = warm.water_relaxation if warm_start is not None else min(initial_alpha, 0.10)
    temperature_damping = _AdaptiveRelaxation(
        temperature_alpha,
        transport_alpha_floor,
        transport_alpha_cap,
    )
    hexane_damping = _AdaptiveRelaxation(
        hexane_alpha,
        min(0.05, initial_alpha),
        initial_alpha,
    )
    # Water is the active-set-limited family.  Start conservatively even if
    # temperature/hexane can safely use the historical 0.5 Picard step.
    water_damping = _AdaptiveRelaxation(
        water_alpha,
        0.025,
        min(initial_alpha, 0.25),
    )
    previous_water_active_set: tuple[tuple[bool, bool, bool], ...] | None = (
        warm.water_active_set if warm_start is not None else None
    )

    # Convergence acceleration of the vapor TEMPERATURE profile -- the
    # near-neutral mode (ρ≈0.9998) lives entirely in the energy coupling
    # (vapor_T <-> particle Tp); the hexane mass profile (vapor_wV2) already
    # converges fast, so it is left on plain relaxed Picard. See
    # `_dominant_mode_extrapolation`.
    #  - `acc_start`: plain Picard first, so the fast modes decay and the
    #    residual becomes a clean single slow mode before we extrapolate (an
    #    early extrapolation can overshoot to a cold T where the GAB isotherm is
    #    invalid, a_h ≥ 1/K).
    #  - `acc_boost_max`: cap on the geometric boost 1/(1−λ).
    #  - residual-decrease safeguard + `acc_cooldown`: if an accelerated step
    #    grew the residual, take several plain passes before re-engaging. So the
    #    accelerator can never do worse than base Picard.
    #  - `[_T_lo, _T_hi]`: physical guard band; an out-of-band extrapolation is
    #    discarded. `_reset_prev`: force a plain pass after each extrapolation so
    #    the next λ estimate uses clean consecutive residuals.
    acc_start = 8 if minimum_width_ratio >= 1.0 - 1.0e-12 else outer_max_iter + 1
    acc_boost_max = 40.0
    _T_lo = min(vapor_inf.T, T_L_sup) - 40.0
    _T_hi = max(vapor_inf.T, T_L_sup) + 80.0
    acc_f_prev: np.ndarray | None = None
    acc_prev_res = math.inf
    acc_cooldown = 0

    for iterations in range(1, outer_max_iter + 1):
        old_vapor_T = tuple(vapor_T)
        old_vapor_m_water = tuple(vapor_m_water)
        old_vapor_m_hex = tuple(vapor_m_hex)
        old_X1_profile = tuple(X1_profile)
        old_free_water_profile = tuple(free_water_profile)
        old_condensed_profile = tuple(condensed_profile_kg_s)
        old_evaporated_profile = tuple(evaporated_profile_kg_s)
        # Accelerator iterate x_k = the vapor TEMPERATURE profile ENTERING this
        # pass. The relaxed iteration body below is the fixed-point map H; its
        # output h_k is captured after step 4.5 and fed to the extrapolation.
        acc_x_k = np.array(vapor_T)
        lagged_particles_for_block = tuple(particles)
        # -- axial correction sources, lagged from the previous iteration's profile --
        Tp_bulk_profile = tuple(
            pt.shell_volumetric_mean(particles[j].Tp, geometry) for j in range(nz)
        )
        q_condL = axial_laplacian(Tp_bulk_profile, cell_widths_m)
        q_condL = tuple(c.k_mixL * q for q in q_condL)
        # 1. energy balance at particle scale (top -> bottom). Each outer
        # pass re-marches a FRESH Tp cascade from the zone's own entry
        # condition through all nz cells in sequence -- cell j's particle
        # must reflect (j+1)*dt of accumulated residence time, not "however
        # many outer iterations have run so far" (those are a convergence
        # device for the vapor<->particle coupling, not a proxy for
        # residence time; conflating the two was caught and fixed during
        # development -- an earlier draft let each cell's state persist and
        # advance by one dt per OUTER iteration, so every cell ended up with
        # identical residence time regardless of axial position). Each
        # cell's own wpg2 (from the previous iteration's mass cascade, step
        # 3 below) is carried through unchanged -- only Tp cascades here.
        compiled_energy = particle_jit.energy_cascade_with_diagnostics(
            initial_particle,
            particles,
            dwpg2_dt_prev,
            vapor_T,
            q_condL,
            X1_profile,
            c.hQ,
            cell_residence_s,
            geometry,
            c.particle,
            _TP_MIN_K,
            _TP_MAX_K,
            cascade_invariants,
            time_substeps=particle_time_substeps,
            moisture_in=X1_in,
            water_enthalpy_reference_k=c.vapor_enthalpy_ref.T_boil_water,
        )
        if compiled_energy is not None:
            compiled_temperatures, compiled_particle_energy_terms = compiled_energy
            particles = [
                pt.ParticleState(
                    wpg2=particles[j].wpg2,
                    Tp=tuple(float(value) for value in compiled_temperatures[j]),
                )
                for j in range(nz)
            ]
            particle_energy_terms_per_particle = [
                tuple(float(value) for value in compiled_particle_energy_terms[j])
                for j in range(nz)
            ]
        else:
            new_particles_energy = []
            particle_energy_terms_per_particle = []
            running_Tp = initial_particle.Tp
            for j in range(nz):
                seed = pt.ParticleState(wpg2=particles[j].wpg2, Tp=running_Tp)
                # eq. A.30 sorption/desorption heat source for the PARTICLE energy
                # march (step 1). No longer also fed to the vapor balance -- the
                # vapor now uses Coletto's own SVm2·Ĥ2 mass-enthalpy term (step 2),
                # not this particle-volume sorption sink.
                updated = seed
                dt = cell_residence_s[j]
                substeps = particle_time_substeps[j]
                particle_dt = dt / substeps
                entering_w = (
                    initial_particle.wpg2
                    if j == 0
                    else particles[j - 1].wpg2
                )
                sorption_entering = pt.ParticleState(
                    wpg2=entering_w,
                    Tp=(
                        initial_particle.Tp
                        if j == 0
                        else particles[j - 1].Tp
                    ),
                )
                sorption_sources = pt.finite_sorption_heat_source_per_layer_w_m3(
                    sorption_entering,
                    particles[j],
                    dt,
                    c.particle,
                )
                sorption_source_w = sum(
                    source * volume
                    for source, volume in zip(sorption_sources, geometry.volumes)
                )
                surface_heat_w = 0.0
                moisture_start = X1_in if j == 0 else X1_profile[j - 1]
                for substep in range(substeps):
                    moisture_old = moisture_start + (
                        X1_profile[j] - moisture_start
                    ) * substep / substeps
                    moisture_new = moisture_start + (
                        X1_profile[j] - moisture_start
                    ) * (substep + 1) / substeps
                    full_sources = tuple(s + q_condL[j] for s in sorption_sources)
                    updated, energy_diagnostics = pt.march_particle_energy(
                        updated,
                        vapor_T[j],
                        c.hQ,
                        full_sources,
                        particle_dt,
                        geometry,
                        c.particle,
                        X1=moisture_new,
                        X1_previous=moisture_old,
                        water_enthalpy_reference_k=c.vapor_enthalpy_ref.T_boil_water,
                    )
                    surface_heat_w += (
                        energy_diagnostics.heat_flux_to_vapor_w_m2
                        * geometry.face_areas[-1]
                        / substeps
                    )
                # Robustness clamp to a physical DT temperature band (see
                # _TP_MIN_K/_TP_MAX_K). It never engages at realistic points.
                if (
                    updated.Tp[0] < _TP_MIN_K
                    or updated.Tp[-1] > _TP_MAX_K
                    or any(t < _TP_MIN_K or t > _TP_MAX_K for t in updated.Tp)
                ):
                    updated = pt.ParticleState(
                        wpg2=updated.wpg2,
                        Tp=tuple(min(max(t, _TP_MIN_K), _TP_MAX_K) for t in updated.Tp),
                    )
                new_particles_energy.append(updated)
                base_cv = (
                    c.particle.alpha_pg * c.particle.rho_pg * c.particle.cp_pg
                    + c.particle.alpha_ps * c.particle.rho_ps * c.particle.cp_ps
                )
                water_cv = (
                    c.particle.alpha_ps
                    * c.particle.rho_ps
                    * c.particle.cp_water_liquid
                )
                T_ref_particle = c.vapor_enthalpy_ref.T_boil_water
                accumulation_w = sum(
                    (
                        (base_cv + water_cv * X1_profile[j])
                        * (temperature_out - T_ref_particle)
                        - (base_cv + water_cv * moisture_start)
                        * (temperature_in - T_ref_particle)
                    )
                    * volume
                    / dt
                    for temperature_out, temperature_in, volume in zip(
                        updated.Tp,
                        seed.Tp,
                        geometry.volumes,
                    )
                )
                axial_source_w = q_condL[j] * geometry.total_volume
                particle_energy_terms_per_particle.append(
                    (
                        accumulation_w,
                        sorption_source_w,
                        axial_source_w,
                        surface_heat_w,
                    )
                )
                running_Tp = updated.Tp
            particles = new_particles_energy
        if thermal_coupling == "block_implicit":
            # The component/active-set outer iteration still needs a
            # provisional temperature before the fixed-state block is
            # assembled later in this pass.  Retain the previous accepted
            # thermal state here instead of applying the noncontractive
            # partitioned particle update that block mode is designed to
            # replace.  The mass state is unchanged at this point, and the
            # block solve below supersedes these placeholder energy terms
            # before any residual can be accepted.
            particles = list(lagged_particles_for_block)
            particle_energy_terms_per_particle = [
                (0.0, 0.0, 0.0, 0.0) for _ in range(nz)
            ]
        particle_temperature_clamped = [
            any(t <= _TP_MIN_K or t >= _TP_MAX_K for t in particle.Tp)
            for particle in particles
        ]
        particle_interphase_temperature = [
            (
                old_vapor_T[j]
                + particle_energy_terms_per_particle[j][3]
                / (c.hQ * geometry.face_areas[-1])
                if c.hQ * geometry.face_areas[-1] > 0.0
                else pt.outer_layer_value(particles[j].Tp)
            )
            for j in range(nz)
        ]

        # 2. energy balance at bed scale (bottom -> top). A per-cell IMPLICIT
        # (backward) relaxation step -- not the naive explicit forward step
        # eq. A.35 might suggest -- because the particle<->vapor transfer
        # coefficient is typically stiff relative to a practical cell size
        # (its own relaxation length can be much shorter than dz); backward
        # Euler is unconditionally stable for this linear relaxation form,
        # mirroring why the particle-scale march is implicit for the same
        # stiffness reason.
        new_vapor_T = (
            list(vapor_T)
            if thermal_coupling == "block_implicit"
            else [0.0] * nz
        )
        raw_condensed_kg_s = [0.0] * nz
        required_evaporation_kg_s = [0.0] * nz
        T_running = vapor_inf.T
        water_running_energy_kg_s = m_water_bottom_kg_s
        hexane_running_energy_kg_s = m_hex_bottom_kg_s
        for j in range(nz - 1, -1, -1):
            if thermal_coupling == "block_implicit":
                # The monolithic projection after the conservative component
                # update supplies the accepted temperature.  Do not evaluate
                # the noncontractive partitioned map merely to obtain a
                # provisional mass-transfer temperature.
                continue
            Tp12 = pt.outer_layer_value(particles[j].Tp)
            cell_volume_m3 = A_bed * cell_widths_m[j]
            water_in_energy_kg_s = water_running_energy_kg_s
            hexane_in_energy_kg_s = hexane_running_energy_kg_s
            # The conservative particle/bed component march is authoritative.
            # Use its lagged accepted cell increment in the enthalpy equation
            # instead of independently re-evaluating hM at the interface. The
            # latter differs from the finite-volume particle inventory at
            # finite radial resolution and previously created an unpaired
            # hexane-enthalpy source even though component mass was closed.
            accepted_hexane_in_kg_s = (
                m_hex_bottom_kg_s if j == nz - 1 else vapor_m_hex[j + 1]
            )
            hexane_transfer_energy_kg_s = (
                vapor_m_hex[j] - accepted_hexane_in_kg_s
            )
            wet_equilibrium_temperature = thermo.dew_point_temperature(
                vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
                c.antoine_water,
                P=c.pressure_pa,
            )
            wet_phase_change_enthalpy_j_kg = thermo.water_phase_change_enthalpy(
                wet_equilibrium_temperature,
                c.vapor_enthalpy_ref,
                c.particle.cp_water_liquid,
            )
            if free_water_profile[j] > 1.0e-12:
                bound_surface_temperature = min(Tp12, wet_equilibrium_temperature)
            else:
                bound_surface_temperature = Tp12
            enthalpy_result = solve_vapor_cell_enthalpy(
                water_in_kg_s=water_in_energy_kg_s,
                hexane_in_kg_s=hexane_in_energy_kg_s,
                temperature_in_k=T_running,
                hexane_transfer_kg_s=hexane_transfer_energy_kg_s,
                bound_water_to_solid_kg_s=(
                    bound_water_mass_rate_w_m3[j] * cell_volume_m3
                ),
                free_water_evaporation_kg_s=evaporated_profile_kg_s[j],
                external_duty_w=q_Iv_profile[j] * cell_volume_m3,
                interphase_ua_w_k=(
                    c.hQ * c.aV * c.alpha_L * cell_volume_m3
                ),
                particle_surface_temperature_k=Tp12,
                bound_water_surface_temperature_k=bound_surface_temperature,
                free_water_temperature_k=wet_equilibrium_temperature,
                c=c,
                interphase_particle_temperature_k=particle_interphase_temperature[j],
            )
            T_running = enthalpy_result.temperature_k
            water_running_energy_kg_s = enthalpy_result.water_out_kg_s
            hexane_running_energy_kg_s = enthalpy_result.hexane_out_kg_s
            raw_condensed_kg_s[j] = enthalpy_result.condensed_water_kg_s
            cell_thermal_coefficient_w_k = (
                enthalpy_result.water_out_kg_s
                * c.vapor_enthalpy_ref.cp_water_vapor
                + enthalpy_result.hexane_out_kg_s
                * c.vapor_enthalpy_ref.cp_hexane_vapor
                + c.hQ * c.aV * c.alpha_L * cell_volume_m3
            )
            # One interface has one signed phase-change rate. Condensation
            # and evaporation are complementary, never simultaneous. The old
            # additive update could condense and immediately re-evaporate the
            # same water in one cell, producing a large artificial internal
            # reflux while both component totals still appeared conservative.
            if raw_condensed_kg_s[j] > 1.0e-12:
                required_evaporation_kg_s[j] = 0.0
            elif (
                free_water_profile[j] > 1.0e-12
                or evaporated_profile_kg_s[j] > 1.0e-12
            ):
                required_evaporation_kg_s[j] = (
                    evaporation_required_for_wet_temperature_constraint(
                        T_running,
                        wet_equilibrium_temperature,
                        evaporated_profile_kg_s[j],
                        cell_thermal_coefficient_w_k,
                        wet_phase_change_enthalpy_j_kg,
                    )
                )
            new_vapor_T[j] = T_running

            # SUPERSEDED IMPLEMENTATION NOTE: the commentary below documents
            # the former independent SVm2 evaluation. The executable
            # conservative scheme above now uses the authoritative discrete
            # hexane component increment, so mass and transferred enthalpy
            # cannot diverge at finite radial resolution.
            # eq. A.34/A.37 vapor energy source, EXACTLY as Coletto prints it
            # (D1, GROUNDING_MATRIX.md):
            #   S_VQ = −a_v α_L J_QR·ř  +  SVm2·Ĥ2  +  q̇_Iv  +  ṁ'_ax·Ĥ2
            # - `−a_v α_L J_QR` (particle<->vapor convection) is the
            #   `kappa_e*(Tp12 − TV)` relaxation just below.
            # - `SVm2·Ĥ2 + ṁ'_ax·Ĥ2` is the enthalpy carried by the hexane MASS
            #   transferred into the vapor. `SVm2` is the SAME bed<->particle
            #   hexane transfer the mass balance (step 4) uses, `hM ρV aV αL
            #   (wpg2_12 − wV2)`; `ṁ'_ax` is the axial hexane flux `m_ax_net`;
            #   `Ĥ2` is the hexane specific vapor enthalpy (same datum machinery
            #   zones/ftrz.py uses, evaluated at the local vapor temperature).
            # This REPLACES the previous non-paper `−sorption_sink·α_L`
            # substitution, which injected the particle's FULL heat of sorption
            # (mostly latent) into the vapor's SENSIBLE temperature balance
            # (A.25 LHS is ρV·CPV·TV) -> the over-heating fixed point. Paired
            # with particle.py's restored literal A.30 sign; validated as a PAIR.
            # Ĥ2 is the hexane's SENSIBLE vapor enthalpy relative to the local
            # vapor temperature: the bed-scale transfer moves hexane that is
            # ALREADY vapor (pore gas -> bulk vapor), so no phase change occurs
            # here -- the desorption phase change lives entirely at the particle
            # scale (A.30). It arrives at the particle-surface temperature Tp12
            # and equilibrates into the vapor at TV, contributing only
            # `cp_hex·(Tp12 − TV)` of sensible heat. (A.25's LHS is the SENSIBLE
            # ρV·CPV·TV; putting the hexane's LATENT `dH_vap_hexane` here instead
            # would inject phase-change energy into a sensible balance -- tested
            # directly this session: it drove the DT meal to ~187 C. The latent
            # heat is carried by the composition wV2, tracked by the mass
            # balance A.24, not by this temperature balance.)

            # MOISTURE (H2O) BALANCE, see module docstring. Two distinct
            # supersaturation checks are needed, not one:

            # (a) the INFLOW itself may already be below its own dew point
            # (e.g. SP1's own steam+upstream-vapor mix, `dt_solver.py`'s
            # `T_bottom`) -- a coarse cell's own strong indirect duty can
            # then re-superheat it within that SAME cell, so checking only
            # the cell's OUTPUT (b, below) would silently miss condensation
            # that physically must happen right at entry, before any
            # convective heating. Flash-condense against T_in alone (no
            # kappa_e/source terms -- those belong to the cell's own march,
            # not this pre-conditioning step): released latent heat reheats
            # the remaining (approximately unchanged mass) water stream by
            # `T_dew_j - T_in`.


            # (b) the cell's OWN convective/duty/sorption balance may cool it
            # back below the (possibly just-corrected) dew point -- cap and
            # back-solve the condensed mass the SAME relaxation equation
            # implies, rather than let the vapor cool below what its own
            # composition allows. Closed-form (the equation is linear in
            # `source`), not a root-find -- unlike `zones/ftrz.py`'s own
            # V-SAT branch, which needs `brentq` because its energy balance
            # is nonlinear in the condensed mass (composition-dependent
            # enthalpy); this one isn't.

            # Conservative wet-cell phase constraint. The evaporation rate
            # from the previous Picard pass enters as a real enthalpy sink.
            # The exact inverse of this cell's discrete energy equation then
            # proposes the TOTAL evaporation needed to reach T_sat. Step 4.5
            # applies the solid-water inventory bound and marches that same
            # mass into the vapor. There is deliberately no temperature clip:
            # if liquid water is unavailable, the temperature remains above
            # T_sat and the constraint is reported water-limited.
        # Under-relax the coupled particle/vapor fixed point. The accepted
        # state is still audited below, so relaxation cannot be mistaken for
        # first-law closure.
        alpha_temperature = temperature_damping.value
        vapor_T = [
            vapor_T[j] + alpha_temperature * (new_vapor_T[j] - vapor_T[j])
            for j in range(nz)
        ]
        alpha_water = water_damping.value
        condensed_profile_kg_s = [
            condensed_profile_kg_s[j]
            + alpha_water * (
                raw_condensed_kg_s[j] - condensed_profile_kg_s[j]
            )
            for j in range(nz)
        ]
        # A transient Picard proposal must not assign the same entering vapor
        # inventory to condensation in several cells.  The sequential raw
        # march normally guarantees this; retain an explicit global projection
        # after relaxation so remeshed/warm-started profiles are conservative
        # as well.  Solid evaporation is credited on the following coupled
        # pass, so this bound is deliberately conservative during iteration.
        total_condensation_kg_s = sum(condensed_profile_kg_s)
        recycled_water_kg_s = sum(evaporated_profile_kg_s) + sum(
            max(-rate * A_bed * cell_widths_m[j], 0.0)
            for j, rate in enumerate(bound_water_mass_rate_w_m3)
        )
        condensation_inventory_kg_s = (
            m_water_bottom_kg_s + recycled_water_kg_s
        )
        if total_condensation_kg_s > condensation_inventory_kg_s:
            condensation_scale = (
                condensation_inventory_kg_s / total_condensation_kg_s
            )
            condensed_profile_kg_s = [
                value * condensation_scale for value in condensed_profile_kg_s
            ]

        # 3. mass balance at particle scale (top -> bottom). Same cascade
        # logic as step 1: a fresh wpg2 cascade from the zone entry
        # condition through all nz cells, each seeded with that cell's own
        # (just-updated in step 1) Tp.
        compiled_mass = particle_jit.mass_cascade(
            initial_particle,
            particles,
            vapor_wV2,
            c.hM,
            c.rho_V,
            cell_residence_s,
            geometry,
            c.particle,
            cascade_invariants,
            time_substeps=particle_time_substeps,
        )
        if compiled_mass is None:
            new_particles_mass = []
            new_rates = []
            running_wpg2 = initial_particle.wpg2
            for j in range(nz):
                seed = pt.ParticleState(wpg2=running_wpg2, Tp=particles[j].Tp)
                updated = seed
                dt = cell_residence_s[j]
                substeps = particle_time_substeps[j]
                particle_dt = dt / substeps
                for _ in range(substeps):
                    updated, _ = pt.march_particle_mass(
                        updated,
                        vapor_wV2[j],
                        c.hM,
                        c.rho_V,
                        particle_dt,
                        geometry,
                        c.particle,
                    )
                diag = pt.MassDiagnostics(
                    dwpg2_dt=tuple(
                        (updated.wpg2[i] - seed.wpg2[i]) / dt
                        for i in range(c.particle.Np)
                    ),
                    hexane_flux_to_vapor_kg_m2_s=(
                        c.hM * c.rho_V * (updated.wpg2[-1] - vapor_wV2[j])
                    ),
                )
                new_particles_mass.append(updated)
                new_rates.append(diag.dwpg2_dt)
                running_wpg2 = updated.wpg2
            particles = new_particles_mass
            dwpg2_dt_prev = new_rates
            x2_particle = [
                pt.shell_volumetric_mean(
                    tuple(_x2_so(w, t, c.particle) for w, t in zip(cell.wpg2, cell.Tp)),
                    geometry,
                )
                for cell in particles
            ]
        else:
            compiled_w, compiled_rates, compiled_x2 = compiled_mass
            particles = [
                pt.ParticleState(
                    wpg2=tuple(float(value) for value in compiled_w[j]),
                    Tp=particles[j].Tp,
                )
                for j in range(nz)
            ]
            dwpg2_dt_prev = [tuple(float(value) for value in compiled_rates[j]) for j in range(nz)]
            x2_particle = [float(value) for value in compiled_x2]

        # 4. Hexane component balance at bed scale (bottom -> top). Use the
        # exact A.22 inventory carried by the particle march: sorbed/oil
        # hexane plus connected pore vapor on ``rho_pg``. The former code
        # injected eq. 5's pure-saturated-vapor term only at the zone inlet,
        # then used an X2 diagnostic that dropped pore vapor in later cells.
        # That mixed two density conventions and created a cross-scale source.
        # The reported meal KPI remains ``x2_particle`` (solid/oil retained);
        # only the component ledger uses the full conserved inventory.
        x2_particle_total = [
            _x2_total_inventory_loading(particle, geometry, c.particle)
            for particle in particles
        ]
        x2_zone_in_total = _x2_total_inventory_loading(
            initial_particle,
            geometry,
            c.particle,
        )
        hexane_release_kg_s = [
            m_dry_kg_s
            * (
                (
                    x2_zone_in_total
                    if j == 0
                    else x2_particle_total[j - 1]
                )
                - x2_particle_total[j]
            )
            for j in range(nz)
        ]
        new_vapor_m_hex, _ = _solve_vapor_hexane_dispersion(
            vapor_water_kg_s=vapor_m_water,
            hexane_release_kg_s=hexane_release_kg_s,
            bottom_hexane_kg_s=m_hex_bottom_kg_s,
            cell_widths_m=cell_widths_m,
            bed_area_m2=A_bed,
            alpha_v=c.alpha_V,
            axial_diffusivity_m2_s=c.D_ax,
            vapor_density_kg_m3=c.rho_V,
            initial_hexane_kg_s=vapor_m_hex,
        )

        alpha_hexane = hexane_damping.value
        vapor_m_hex = [
            max(
                vapor_m_hex[j] + alpha_hexane * (new_vapor_m_hex[j] - vapor_m_hex[j]),
                0.0,
            )
            for j in range(nz)
        ]
        vapor_wV2 = [
            vapor_m_hex[j] / max(vapor_m_water[j] + vapor_m_hex[j], 1.0e-12) for j in range(nz)
        ]

        # 4.5. solid moisture at bed scale (top -> bottom, matching solid
        # flow direction -- see module docstring's "MOISTURE (H2O) BALANCE").
        # Uses THIS iteration's own just-updated `vapor_wV2`/`vapor_T` (steps
        # 2/4 above) and `condensed_profile_kg_s` (step 2's relaxed supersaturation
        # branches); produces an explicit bound-water mass-transfer profile
        # consumed by step 2 of the next outer iteration. The vapor enthalpy
        # cell converts that same mass to a phase duty on the common datum.
        new_X1_profile = [0.0] * nz
        new_free_water_profile = [0.0] * nz
        raw_evaporated_kg_s = [0.0] * nz
        X1_running = X1_in
        free_water_running = max(X1_in - saturated_bound_water_capacity, 0.0)
        # Water-availability budget (mass conservation): the solid can adsorb at
        # most the water the vapor carries into the zone. Without this cap a
        # high-a_w regime (near-saturated steam, e.g. under an elevated DT
        # pressure) drives the subsaturated adsorption below to pull MORE water
        # onto the meal than exists in the vapor, sending the vapor's own water
        # flow negative (found via the FTRZ handoff crash under a raised
        # pressure). Condensation (step 2) is already limited by its own
        # `water_remaining_kg_s`; this budget adds the same discipline to the
        # subsaturated branch and shares the running total across both.
        #
        # The whole zone's condensation is PRE-COUNTED into the shared budget here,
        # before the isotherm adsorption draws on it: the vapor rises bottom->top, so
        # water condensed low (at the sparge) is gone before it reaches the meal higher
        # up. Counting condensation only as each condensed cell is REACHED walking the
        # solid top->bottom instead let a top isotherm cell adsorb against the FULL
        # budget before a bottom condensed cell debited it -- both branches drawing the
        # same water, so the meal could gain more than the steam actually supplied
        # (found at an under-set sparge: DT exit ~29%wb, ~2.8 kg/s water conjured from
        # nothing, dome vapor water driven to 0). Pre-counting caps the combined
        # condensation+adsorption at the actual inflowing vapor water. At the calibrated
        # operating point the sparge is strong enough that the meal reaches ~19% via the
        # isotherm alone and this is slack; it binds only when the sparge is weak.
        water_available_kg_s = (1.0 - vapor_inf.wV2) * m_vapor_kg_s
        cum_water_to_solid_kg_s = sum(condensed_profile_kg_s)
        water_limited_cells: list[bool] = []
        evaporation_limited_cells: list[bool] = []
        bound_water_bvp_profile: list[float] | None = None
        phase_free_bound_water = (
            max(X1_in - saturated_bound_water_capacity, 0.0) <= 1.0e-12
            and max(free_water_profile, default=0.0) <= 1.0e-12
            and max(condensed_profile_kg_s, default=0.0) <= 1.0e-12
            and max(raw_condensed_kg_s, default=0.0) <= 1.0e-12
            and max(evaporated_profile_kg_s, default=0.0) <= 1.0e-12
            and max(required_evaporation_kg_s, default=0.0) <= 1.0e-12
        )
        if phase_free_bound_water:
            bound_water_bvp_profile, _ = solve_bound_water_countercurrent_bvp(
                moisture_in=X1_in,
                bottom_vapor_water_kg_s=m_water_bottom_kg_s,
                vapor_hexane_top_kg_s=vapor_m_hex,
                particle_surface_temperature_k=[
                    pt.outer_layer_value(particle.Tp) for particle in particles
                ],
                dry_solid_flow_kg_s=m_dry_kg_s,
                rate_s_inv=kappa_w,
                cell_residence_s=cell_residence_s,
                bound_capacity=saturated_bound_water_capacity,
                luikov=c.luikov,
                antoine_water=c.antoine_water,
                pressure_pa=c.pressure_pa,
            )
        for j in range(nz):
            if bound_water_bvp_profile is not None:
                new_X1_profile[j] = bound_water_bvp_profile[j]
                new_free_water_profile[j] = 0.0
                raw_evaporated_kg_s[j] = 0.0
                water_limited_cells.append(False)
                evaporation_limited_cells.append(False)
                X1_running = new_X1_profile[j]
                free_water_running = 0.0
                continue
            # Bulk condensation is an immediate vapor-energy event. Credit
            # its relaxed active-set mass first, then let the wet solid make
            # a finite-rate bidirectional adjustment toward Xe. Keeping the
            # mechanisms additive removes the discontinuous either/or branch
            # that produced the low-duty three-cycle.
            bound_water_before = max(X1_running - free_water_running, 0.0)
            free_water_after_condensation = (
                free_water_running + condensed_profile_kg_s[j] / m_dry_kg_s
            )
            water_available_from_solid_kg_s = m_dry_kg_s * free_water_after_condensation
            evaporation_kg_s = (
                0.0
                if condensed_profile_kg_s[j] > 1.0e-12
                else min(
                    required_evaporation_kg_s[j],
                    water_available_from_solid_kg_s,
                )
            )
            evaporation_limited = (
                required_evaporation_kg_s[j] - evaporation_kg_s > PROFILE_TOL_FLOW
            )
            evaporation_limited_cells.append(evaporation_limited)
            raw_evaporated_kg_s[j] = evaporation_kg_s
            free_water_after_evaporation = (
                free_water_after_condensation - evaporation_kg_s / m_dry_kg_s
            )
            # Evaporation returns water to the counter-current vapor stream,
            # increasing the amount that later cells may sorb.
            cum_water_to_solid_kg_s -= evaporation_kg_s
            Y_V2_j = vapor_wV2[j] / max(1.0 - vapor_wV2[j], 1.0e-12)
            particle_surface_temperature = pt.outer_layer_value(particles[j].Tp)
            # Only an actual free-liquid interface can impose the saturation
            # temperature. Bound moisture sees the particle surface's real
            # temperature; retaining the old saturation cap here would keep
            # superheated meal artificially humid and suppress its physical
            # desorption/latent-heat sink.
            if free_water_after_evaporation > 1.0e-12:
                T_surface = min(
                    particle_surface_temperature,
                    thermo.dew_point_temperature(
                        vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
                        c.antoine_water,
                        P=c.pressure_pa,
                    ),
                )
            else:
                T_surface = particle_surface_temperature
            a_w = thermo.water_activity(Y_V2_j, T_surface, c.antoine_water, P=c.pressure_pa)
            a_w = min(max(a_w, 1.0e-9), thermo.LUIKOV_MAX_VALIDATED_UR)
            Xe = thermo.luikov_equilibrium_moisture(a_w, c.luikov)
            bound_water_new = exact_linear_driving_force_step(
                bound_water_before,
                Xe,
                kappa_w,
                cell_residence_s[j],
            )
            mass_rate_kg_s = m_dry_kg_s * (bound_water_new - bound_water_before)
            # A cell actively condensing bulk water is saturated; do not let
            # the clamped bound-water isotherm immediately evaporate that
            # same free condensate. Positive sorption remains additive.
            if condensed_profile_kg_s[j] > 1.0e-12 and mass_rate_kg_s < 0.0:
                mass_rate_kg_s = 0.0
                bound_water_new = bound_water_before
            # The active free-water evaporation constraint takes precedence
            # over an opposing bound-water adsorption step in the same
            # lumped cell. This is the complementarity analogue of avoiding
            # simultaneous evaporation and condensation at one interface.
            if evaporation_kg_s > 1.0e-12 and mass_rate_kg_s > 0.0:
                mass_rate_kg_s = 0.0
                bound_water_new = bound_water_before
            if mass_rate_kg_s > 0.0:
                headroom_kg_s = max(water_available_kg_s - cum_water_to_solid_kg_s, 0.0)
                water_limited = mass_rate_kg_s > headroom_kg_s
                if mass_rate_kg_s > headroom_kg_s:
                    mass_rate_kg_s = headroom_kg_s
                    bound_water_new = bound_water_before + mass_rate_kg_s / m_dry_kg_s
            else:
                water_limited = False
            water_limited_cells.append(water_limited)
            cum_water_to_solid_kg_s += mass_rate_kg_s
            new_free_water_profile[j] = free_water_after_evaporation
            new_X1_profile[j] = bound_water_new + free_water_after_evaporation
            X1_running = new_X1_profile[j]
            free_water_running = free_water_after_evaporation
        # A nonlinear trial profile can be globally water-bounded yet locally
        # impossible for counter-current vapor: water released by a higher
        # solid cell cannot travel downward to supply a lower cell.  Backtrack
        # only the Picard step until every bottom-to-top vapor control volume
        # has non-negative water.  This line search does not alter a feasible
        # fixed point; it prevents an intermediate iterate from being
        # misclassified as a thermodynamic impossibility.
        accepted_alpha_water = alpha_water
        while True:
            candidate_X1_profile = [
                X1_profile[j]
                + accepted_alpha_water * (new_X1_profile[j] - X1_profile[j])
                for j in range(nz)
            ]
            _, minimum_trial_water_kg_s = _countercurrent_vapor_water_profile(
                candidate_X1_profile,
                X1_in,
                m_dry_kg_s,
                m_water_bottom_kg_s,
            )
            if minimum_trial_water_kg_s >= -BALANCE_TOL_FLOW:
                break
            accepted_alpha_water *= 0.5
            if accepted_alpha_water < 1.0e-12:
                raise ValueError(
                    "DCZ water-profile line search could not recover a "
                    "locally feasible counter-current component state"
                )
        if accepted_alpha_water < alpha_water:
            water_damping.value = accepted_alpha_water
        alpha_water = accepted_alpha_water
        X1_profile = candidate_X1_profile
        free_water_profile = [
            free_water_profile[j]
            + alpha_water * (new_free_water_profile[j] - free_water_profile[j])
            for j in range(nz)
        ]
        evaporated_profile_kg_s = [
            evaporated_profile_kg_s[j]
            + alpha_water * (
                raw_evaporated_kg_s[j] - evaporated_profile_kg_s[j]
            )
            for j in range(nz)
        ]
        for j in range(nz):
            net_condensation_kg_s = (
                condensed_profile_kg_s[j] - evaporated_profile_kg_s[j]
            )
            condensed_profile_kg_s[j] = (
                max(net_condensation_kg_s, 0.0)
                if net_condensation_kg_s > 1.0e-12
                else 0.0
            )
            evaporated_profile_kg_s[j] = (
                max(-net_condensation_kg_s, 0.0)
                if net_condensation_kg_s < -1.0e-12
                else 0.0
            )
        # The accepted X1 profile is the component ledger.  Recover the bound
        # transfer as the remainder after the accepted condensation and free
        # evaporation rates, instead of relaxing it as a fourth independent
        # water state.  This makes every cell satisfy
        #   m_dry*dX1 = m_bound + m_condensed - m_evaporated
        # identically while retaining the physical phase split needed by the
        # enthalpy balance.
        bound_water_mass_rate_w_m3 = []
        for j in range(nz):
            solid_water_in = X1_in if j == 0 else X1_profile[j - 1]
            total_water_to_solid_kg_s = m_dry_kg_s * (
                X1_profile[j] - solid_water_in
            )
            bound_water_mass_rate_w_m3.append(
                (
                    total_water_to_solid_kg_s
                    - condensed_profile_kg_s[j]
                    + evaporated_profile_kg_s[j]
                )
                / (A_bed * cell_widths_m[j])
            )

        # March the accepted solid-moisture increments into the counter-current
        # vapor as an exact discrete component ledger.  The solid profile is
        # authoritative, just as the particle X2 cascade is authoritative for
        # hexane above.  Relaxing vapor water independently used to preserve a
        # warm-start inconsistency by (1-alpha_water) on every pass; at the
        # minimum adaptive alpha=0.001 this produced an apparently converged
        # profile with a technologically material whole-zone water residual.
        # Algebraic reconstruction removes that numerical storage term:
        #
        #   m_w,V,out,j = m_w,V,in,j - m_dry (X1_out,j - X1_in,j).
        #
        # Individual condensation, free evaporation, and bound-water rates
        # remain in the Picard state because the enthalpy equation needs their
        # distinct phase enthalpies.  Their fixed-point consistency is checked
        # separately through the profile residuals.
        new_vapor_m_water, minimum_vapor_water_kg_s = (
            _countercurrent_vapor_water_profile(
                X1_profile,
                X1_in,
                m_dry_kg_s,
                m_water_bottom_kg_s,
            )
        )
        if minimum_vapor_water_kg_s < -BALANCE_TOL_FLOW:
            raise ValueError(
                "DCZ accepted water profile has a negative local vapor "
                f"inventory ({minimum_vapor_water_kg_s:.6g} kg/s)"
            )
        new_vapor_m_water = [max(value, 0.0) for value in new_vapor_m_water]
        vapor_m_water = new_vapor_m_water
        vapor_wV2 = [
            vapor_m_hex[j] / max(vapor_m_water[j] + vapor_m_hex[j], 1.0e-12) for j in range(nz)
        ]

        # Project the accepted component and phase-transfer state onto its
        # exact enthalpy solution.  In block mode the phase-free affine block is
        # the fast path whenever no cell condenses and every cell stays above
        # its dew curve; otherwise the coupled phase/component semismooth-Newton
        # block solves particle temperature, vapor temperature, vapor water
        # flow, solid moisture, and the condensation rate simultaneously.
        if thermal_coupling == "block_implicit":
            phase_free_guess = (
                max(condensed_profile_kg_s, default=0.0) <= 1.0e-12
                and max(evaporated_profile_kg_s, default=0.0) <= 1.0e-12
                and max(free_water_profile, default=0.0) <= 1.0e-12
            )
            used_coupled_block = True
            if phase_free_guess:
                thermal_block = _solve_phase_fixed_thermal_block(
                    initial_particle=initial_particle,
                    particles=particles,
                    vapor_temperature_guess_k=vapor_T,
                    vapor_water_kg_s=vapor_m_water,
                    vapor_hexane_kg_s=vapor_m_hex,
                    bottom_water_kg_s=m_water_bottom_kg_s,
                    bottom_hexane_kg_s=m_hex_bottom_kg_s,
                    bottom_temperature_k=vapor_inf.T,
                    bound_water_mass_rate_w_m3=bound_water_mass_rate_w_m3,
                    free_water_profile=free_water_profile,
                    condensed_water_kg_s=condensed_profile_kg_s,
                    evaporated_water_kg_s=evaporated_profile_kg_s,
                    moisture_profile=X1_profile,
                    moisture_in=X1_in,
                    q_indirect_w_m3=q_Iv_profile,
                    q_axial_w_m3=q_condL,
                    cell_widths_m=cell_widths_m,
                    cell_residence_s=cell_residence_s,
                    time_substeps=particle_time_substeps,
                    bed_area_m2=A_bed,
                    geometry=geometry,
                    c=c,
                )
                subdew_cells = [
                    j
                    for j in range(nz)
                    if thermal_block.vapor_temperature_k[j]
                    < thermo.dew_point_temperature(
                        vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
                        c.antoine_water,
                        P=c.pressure_pa,
                    )
                    - PHASE_EQUILIBRIUM_TOL_K
                ]
                if not subdew_cells:
                    used_coupled_block = False
                    particles = list(thermal_block.particles)
                    vapor_T = list(thermal_block.vapor_temperature_k)
                    particle_energy_terms_per_particle = list(
                        thermal_block.particle_energy_terms_w
                    )
            if used_coupled_block:
                free_in = max(X1_in - saturated_bound_water_capacity, 0.0)
                coupled_block = _solve_phase_coupled_thermal_block(
                    initial_particle=initial_particle,
                    particles=particles,
                    vapor_temperature_guess_k=vapor_T,
                    vapor_water_guess_kg_s=vapor_m_water,
                    vapor_hexane_kg_s=vapor_m_hex,
                    bottom_water_kg_s=m_water_bottom_kg_s,
                    bottom_hexane_kg_s=m_hex_bottom_kg_s,
                    bottom_temperature_k=vapor_inf.T,
                    bound_moisture_guess=[
                        X1_profile[j] - free_water_profile[j] for j in range(nz)
                    ],
                    free_moisture_guess=list(free_water_profile),
                    condensation_guess_kg_s=list(condensed_profile_kg_s),
                    evaporation_guess_kg_s=list(evaporated_profile_kg_s),
                    moisture_in=X1_in,
                    free_in=free_in,
                    q_indirect_w_m3=q_Iv_profile,
                    q_axial_w_m3=q_condL,
                    cell_widths_m=cell_widths_m,
                    cell_residence_s=cell_residence_s,
                    time_substeps=particle_time_substeps,
                    bed_area_m2=A_bed,
                    geometry=geometry,
                    c=c,
                    kappa_w=kappa_w,
                    bound_capacity=saturated_bound_water_capacity,
                    m_dry_kg_s=m_dry_kg_s,
                )
                if not coupled_block.converged:
                    raise ValueError(
                        "DCZ coupled phase/component block did not converge "
                        f"(residual {coupled_block.residual_max_w:.3g} W after "
                        f"{coupled_block.iterations} iterations)"
                    )
                particles = list(coupled_block.particles)
                vapor_T = list(coupled_block.vapor_temperature_k)
                vapor_m_water = list(coupled_block.vapor_water_kg_s)
                free_water_profile = list(coupled_block.free_moisture)
                X1_profile = [
                    coupled_block.bound_moisture[j] + coupled_block.free_moisture[j]
                    for j in range(nz)
                ]
                condensed_profile_kg_s = list(coupled_block.condensed_water_kg_s)
                evaporated_profile_kg_s = list(coupled_block.evaporated_water_kg_s)
                evaporation_limited_cells = [
                    coupled_block.evaporated_water_kg_s[j] > PROFILE_TOL_FLOW
                    and coupled_block.free_moisture[j] <= EXIT_TOL_X1
                    for j in range(nz)
                ]
                bound_water_mass_rate_w_m3 = []
                for j in range(nz):
                    bound_prev = (
                        (X1_in - free_in)
                        if j == 0
                        else coupled_block.bound_moisture[j - 1]
                    )
                    bound_water_mass_rate_w_m3.append(
                        m_dry_kg_s
                        * (coupled_block.bound_moisture[j] - bound_prev)
                        / (A_bed * cell_widths_m[j])
                    )
                vapor_wV2 = [
                    vapor_m_hex[j]
                    / max(vapor_m_water[j] + vapor_m_hex[j], 1.0e-12)
                    for j in range(nz)
                ]
                particle_energy_terms_per_particle = list(
                    coupled_block.particle_energy_terms_w
                )
            particle_interphase_temperature = [
                (
                    vapor_T[j]
                    + particle_energy_terms_per_particle[j][3]
                    / (c.hQ * geometry.face_areas[-1])
                    if c.hQ * geometry.face_areas[-1] > 0.0
                    else pt.outer_layer_value(particles[j].Tp)
                )
                for j in range(nz)
            ]
            particle_temperature_clamped = [
                any(t <= _TP_MIN_K or t >= _TP_MAX_K for t in particle.Tp)
                for particle in particles
            ]
        else:
            projected_vapor_T = [0.0] * nz
            projected_temperature_in = vapor_inf.T
            for j in range(nz - 1, -1, -1):
                cell_volume_m3 = A_bed * cell_widths_m[j]
                water_in_cell = (
                    m_water_bottom_kg_s if j == nz - 1 else vapor_m_water[j + 1]
                )
                hexane_in_cell = (
                    m_hex_bottom_kg_s if j == nz - 1 else vapor_m_hex[j + 1]
                )
                Tp_surface = pt.outer_layer_value(particles[j].Tp)
                wet_equilibrium_temperature = thermo.dew_point_temperature(
                    vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
                    c.antoine_water,
                    P=c.pressure_pa,
                )
                bound_surface_temperature = (
                    min(Tp_surface, wet_equilibrium_temperature)
                    if free_water_profile[j] > 1.0e-12
                    else Tp_surface
                )
                energy_temperature = solve_vapor_cell_temperature_at_state(
                    water_in_kg_s=water_in_cell,
                    hexane_in_kg_s=hexane_in_cell,
                    temperature_in_k=projected_temperature_in,
                    water_out_kg_s=vapor_m_water[j],
                    hexane_out_kg_s=vapor_m_hex[j],
                    bound_water_to_solid_kg_s=(
                        bound_water_mass_rate_w_m3[j] * cell_volume_m3
                    ),
                    free_water_evaporation_kg_s=evaporated_profile_kg_s[j],
                    condensed_water_kg_s=condensed_profile_kg_s[j],
                    external_duty_w=q_Iv_profile[j] * cell_volume_m3,
                    interphase_ua_w_k=(
                        c.hQ * c.aV * c.alpha_L * cell_volume_m3
                    ),
                    particle_surface_temperature_k=Tp_surface,
                    bound_water_surface_temperature_k=bound_surface_temperature,
                    free_water_temperature_k=wet_equilibrium_temperature,
                    c=c,
                    interphase_particle_temperature_k=particle_interphase_temperature[j],
                )
                # Complementarity projection: an accepted gas cell cannot lie
                # below its dew curve.  If condensate is active, retain the dew
                # temperature and let the nonzero energy residual drive the
                # damped condensate rate; if inactive, use the superheated
                # enthalpy solution whenever it is feasible.
                projected_vapor_T[j] = (
                    wet_equilibrium_temperature
                    if condensed_profile_kg_s[j] > 1.0e-12
                    else max(energy_temperature, wet_equilibrium_temperature)
                )
                projected_temperature_in = projected_vapor_T[j]
            vapor_T = projected_vapor_T

        # Thermal coupling changes the temperature used by the GAB
        # equilibrium.  Recompute the reported retained loading from the
        # accepted particle state; keeping the value produced before the
        # thermal block made ``X2_bulk`` stale and broke exact
        # steady-to-transient projection in block mode.
        x2_particle = [
            pt.shell_volumetric_mean(
                tuple(
                    _x2_so(activity, temperature, c.particle)
                    for activity, temperature in zip(
                        particle.wpg2,
                        particle.Tp,
                    )
                ),
                geometry,
            )
            for particle in particles
        ]

        # Accepted-state energy audit.  The cell root solve above is exact for
        # its Picard inputs, but temperature, component flows, and phase rates
        # are subsequently updated.  Auditing those accepted states prevents a
        # small fixed-point increment from hiding a large first-law residual.
        accepted_vapor_energy_residual_w = [0.0] * nz
        particle_accumulation_w_profile = [0.0] * nz
        particle_sorption_w_profile = [0.0] * nz
        particle_axial_w_profile = [0.0] * nz
        particle_surface_w_profile = [0.0] * nz
        particle_energy_residual_w_profile = [0.0] * nz
        interphase_mismatch_w_profile = [0.0] * nz
        combined_energy_residual_w_profile = [0.0] * nz
        phase_equilibrium_temperature_residual = [0.0] * nz
        for j in range(nz):
            cell_volume_m3 = A_bed * cell_widths_m[j]
            representative_particles_per_cell = (
                c.alpha_L * cell_volume_m3 / geometry.total_volume
            )
            particle_terms = particle_energy_terms_per_particle[j]
            particle_accumulation_w_profile[j] = (
                representative_particles_per_cell * particle_terms[0]
            )
            particle_sorption_w_profile[j] = (
                representative_particles_per_cell * particle_terms[1]
            )
            particle_axial_w_profile[j] = (
                representative_particles_per_cell * particle_terms[2]
            )
            particle_surface_w_profile[j] = (
                representative_particles_per_cell * particle_terms[3]
            )
            particle_energy_residual_w_profile[j] = (
                particle_accumulation_w_profile[j]
                - particle_sorption_w_profile[j]
                - particle_axial_w_profile[j]
                + particle_surface_w_profile[j]
            )
            Tp_surface = pt.outer_layer_value(particles[j].Tp)
            vapor_interphase_w = (
                c.hQ
                * c.aV
                * c.alpha_L
                * cell_volume_m3
                * (particle_interphase_temperature[j] - vapor_T[j])
            )
            interphase_mismatch_w_profile[j] = (
                vapor_interphase_w - particle_surface_w_profile[j]
            )
            water_in_cell = (
                m_water_bottom_kg_s if j == nz - 1 else vapor_m_water[j + 1]
            )
            hexane_in_cell = (
                m_hex_bottom_kg_s if j == nz - 1 else vapor_m_hex[j + 1]
            )
            temperature_in_cell = (
                vapor_inf.T if j == nz - 1 else vapor_T[j + 1]
            )
            # The coupled block carries bound-water enthalpy at the particle
            # surface temperature (the solid it enters), so the audit must use
            # the same datum in block mode; the partitioned route keeps its
            # dew-pinned wet-surface convention.
            bound_surface_temperature = (
                Tp_surface
                if thermal_coupling == "block_implicit"
                else (
                    min(
                        Tp_surface,
                        thermo.dew_point_temperature(
                            vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
                            c.antoine_water,
                            P=c.pressure_pa,
                        ),
                    )
                    if free_water_profile[j] > 1.0e-12
                    else Tp_surface
                )
            )
            wet_equilibrium_temperature = thermo.dew_point_temperature(
                vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
                c.antoine_water,
                P=c.pressure_pa,
            )
            if free_water_profile[j] > 1.0e-10:
                phase_equilibrium_temperature_residual[j] = abs(
                    vapor_T[j] - wet_equilibrium_temperature
                )
            else:
                phase_equilibrium_temperature_residual[j] = max(
                    wet_equilibrium_temperature - vapor_T[j],
                    0.0,
                )
            accepted_vapor_energy_residual_w[j] = (
                vapor_cell_enthalpy_residual_at_state(
                    water_in_kg_s=water_in_cell,
                    hexane_in_kg_s=hexane_in_cell,
                    temperature_in_k=temperature_in_cell,
                    water_out_kg_s=vapor_m_water[j],
                    hexane_out_kg_s=vapor_m_hex[j],
                    temperature_out_k=vapor_T[j],
                    bound_water_to_solid_kg_s=(
                        bound_water_mass_rate_w_m3[j] * cell_volume_m3
                    ),
                    free_water_evaporation_kg_s=evaporated_profile_kg_s[j],
                    condensed_water_kg_s=condensed_profile_kg_s[j],
                    external_duty_w=q_Iv_profile[j] * cell_volume_m3,
                    interphase_heat_to_vapor_w=vapor_interphase_w,
                    particle_surface_temperature_k=Tp_surface,
                    bound_water_surface_temperature_k=bound_surface_temperature,
                    free_water_temperature_k=wet_equilibrium_temperature,
                    c=c,
                )
            )
            combined_energy_residual_w_profile[j] = (
                accepted_vapor_energy_residual_w[j]
                + particle_energy_residual_w_profile[j]
                + interphase_mismatch_w_profile[j]
            )

        # -- convergence acceleration of the vapor TEMPERATURE profile --
        # h_k = H(x_k): the vapor_T the relaxed iteration body just produced.
        # f_k = h_k − x_k is the fixed-point residual; ||f_k|| → 0 at convergence.
        acc_h_k = np.array(vapor_T)  # h_k = H(x_k)
        f_k = acc_h_k - acc_x_k
        res_vap = float(np.max(np.abs(f_k)))
        apply_acc = False
        if apply_acc and res_vap > 1.05 * acc_prev_res:
            # Residual GREW -> the previous extrapolation overshot. Take several
            # plain (relaxed-Picard) passes to re-settle before re-engaging.
            acc_cooldown = acc_start
            acc_f_prev = None
            apply_acc = False
        if acc_cooldown > 0:
            acc_cooldown -= 1
        if apply_acc:
            x_next = _dominant_mode_extrapolation(acc_x_k, f_k, acc_f_prev, acc_boost_max)
            if x_next is not None and all(math.isfinite(t) and _T_lo <= t <= _T_hi for t in x_next):
                vapor_T = [float(t) for t in x_next]
                acc_f_prev = None  # force a plain pass next, for a clean λ estimate
            else:
                acc_f_prev = f_k
        else:
            acc_f_prev = f_k
        # vapor_wV2 is ALWAYS the plain relaxed step (not accelerated -- init note).
        acc_prev_res = res_vap

        # 5. Honest fixed-point convergence check.  The previous implementation
        # checked the bottom meal and top component flows but omitted the top
        # vapor temperature consumed by solve_dt.  A very slow thermal mode
        # could therefore be labelled converged while that boundary was still
        # moving by several kelvin.  Check both exported boundaries and maximum
        # full-profile residuals of all vapor states.
        last = particles[-1]
        x2_exit = x2_particle[-1]
        T_exit = pt.shell_volumetric_mean(last.Tp, geometry)
        x1_exit = X1_profile[-1]
        d_x2 = abs(x2_exit - prev_x2_exit)
        d_T = abs(T_exit - prev_T_exit)
        d_x1 = abs(x1_exit - prev_x1_exit)
        d_water_out = abs(vapor_m_water[0] - prev_water_out)
        d_hex_out = abs(vapor_m_hex[0] - prev_hex_out)
        d_vapor_top_T = abs(vapor_T[0] - prev_vapor_top_T)
        d_vapor_T_profile = max(abs(vapor_T[j] - old_vapor_T[j]) for j in range(nz))
        d_vapor_water_profile = max(abs(vapor_m_water[j] - old_vapor_m_water[j]) for j in range(nz))
        d_vapor_hex_profile = max(abs(vapor_m_hex[j] - old_vapor_m_hex[j]) for j in range(nz))
        d_evaporated_profile = max(
            abs(evaporated_profile_kg_s[j] - old_evaporated_profile[j]) for j in range(nz)
        )
        d_free_water_profile = max(
            abs(free_water_profile[j] - old_free_water_profile[j]) for j in range(nz)
        )
        water_boundary_balance = abs(
            (m_water_bottom_kg_s - vapor_m_water[0])
            - m_dry_kg_s * (X1_profile[-1] - X1_in)
        )
        residuals = DCZResiduals(
            solid_out_X2=d_x2,
            solid_out_T=d_T,
            solid_out_X1=d_x1,
            vapor_top_T=d_vapor_top_T,
            vapor_top_water_flow=d_water_out,
            vapor_top_hexane_flow=d_hex_out,
            vapor_profile_T_max=d_vapor_T_profile,
            vapor_profile_water_flow_max=d_vapor_water_profile,
            vapor_profile_hexane_flow_max=d_vapor_hex_profile,
            evaporated_water_flow_max=d_evaporated_profile,
            free_water_max=d_free_water_profile,
            water_boundary_balance=water_boundary_balance,
            vapor_energy_balance_max=max(
                abs(value) for value in accepted_vapor_energy_residual_w
            ),
            interphase_energy_mismatch_max=max(
                abs(value) for value in interphase_mismatch_w_profile
            ),
            phase_equilibrium_temperature_max=max(
                phase_equilibrium_temperature_residual
            ),
        )
        prev_x2_exit, prev_T_exit, prev_x1_exit = x2_exit, T_exit, x1_exit
        prev_water_out, prev_hex_out = vapor_m_water[0], vapor_m_hex[0]
        prev_vapor_top_T = vapor_T[0]
        if residual_log is not None:
            residual_log.append((iterations, d_x2, d_T, d_x1, d_water_out, d_hex_out))
        if progress_callback is not None:
            progress_callback(
                {
                    "level": "dcz_iteration",
                    "iteration": iterations,
                    "maximum_scaled_residual": residuals.maximum_scaled,
                    "solid_out_X2_residual": d_x2,
                    "solid_out_temperature_residual_k": d_T,
                    "solid_out_X1_residual": d_x1,
                    "vapor_top_temperature_residual_k": d_vapor_top_T,
                    "vapor_top_water_residual_kg_s": d_water_out,
                    "vapor_top_hexane_residual_kg_s": d_hex_out,
                    "vapor_profile_temperature_residual_k": d_vapor_T_profile,
                    "vapor_profile_water_residual_kg_s": d_vapor_water_profile,
                    "vapor_profile_hexane_residual_kg_s": d_vapor_hex_profile,
                    "water_boundary_residual_kg_s": water_boundary_balance,
                    "vapor_energy_residual_max_w": residuals.vapor_energy_balance_max,
                    "interphase_mismatch_max_w": (
                        residuals.interphase_energy_mismatch_max
                    ),
                    "phase_equilibrium_residual_max_k": (
                        residuals.phase_equilibrium_temperature_max
                    ),
                }
            )

        water_active_set = tuple(
            (
                raw_condensed_kg_s[j] > 1.0e-12,
                raw_evaporated_kg_s[j] > 1.0e-12,
                water_limited_cells[j] or evaporation_limited_cells[j],
            )
            for j in range(nz)
        )
        active_set_changed = (
            previous_water_active_set is not None and water_active_set != previous_water_active_set
        )
        previous_water_active_set = water_active_set
        if adaptive_damping:
            temperature_damping.observe(d_vapor_T_profile)
            hexane_damping.observe(d_vapor_hex_profile)
            water_norm = max(
                d_vapor_water_profile / PROFILE_TOL_FLOW,
                max(abs(new_X1_profile[j] - old_X1_profile[j]) / EXIT_TOL_X1 for j in range(nz)),
                max(
                    abs(raw_condensed_kg_s[j] - old_condensed_profile[j]) / PROFILE_TOL_FLOW
                    for j in range(nz)
                ),
                d_evaporated_profile / PROFILE_TOL_FLOW,
                d_free_water_profile / EXIT_TOL_X1,
            )
            water_damping.observe(water_norm, active_set_changed)

        coupled_converged = False
        if coupling_tolerance is not None:
            if not math.isfinite(coupling_tolerance) or coupling_tolerance <= 0.0:
                raise ValueError("coupling_tolerance must be finite and positive")
            coupled_converged = (
                residuals.solid_out_X2 <= EXIT_TOL_X2
                and residuals.solid_out_T <= max(EXIT_TOL_T, coupling_tolerance)
                and residuals.solid_out_X1 <= max(EXIT_TOL_X1, coupling_tolerance)
                and residuals.vapor_top_T <= max(PROFILE_TOL_T, coupling_tolerance)
                and residuals.vapor_top_water_flow
                <= max(EXIT_TOL_FLOW, coupling_tolerance)
                and residuals.vapor_top_hexane_flow
                <= max(EXIT_TOL_FLOW, coupling_tolerance)
                and residuals.vapor_profile_T_max
                <= max(PROFILE_TOL_T, coupling_tolerance)
                and residuals.vapor_profile_water_flow_max
                <= max(PROFILE_TOL_FLOW, coupling_tolerance)
                and residuals.vapor_profile_hexane_flow_max
                <= max(PROFILE_TOL_FLOW, coupling_tolerance)
                and residuals.evaporated_water_flow_max
                <= max(PROFILE_TOL_FLOW, coupling_tolerance)
                and residuals.free_water_max <= max(EXIT_TOL_X1, coupling_tolerance)
                # Conservation is an algebraic invariant, not an operating
                # fixed-point tolerance, and must never be weakened by the
                # outer solve's requested accuracy.
                and residuals.water_boundary_balance <= BALANCE_TOL_FLOW
                and residuals.vapor_energy_balance_max <= ENERGY_BALANCE_TOL_W
                and residuals.interphase_energy_mismatch_max <= ENERGY_BALANCE_TOL_W
                and residuals.phase_equilibrium_temperature_max
                <= PHASE_EQUILIBRIUM_TOL_K
            )
        if residuals.maximum_scaled <= 1.0 or coupled_converged:
            converged = True
            break

    cells = []
    for j in range(nz):
        cell_volume_m3 = A_bed * cell_widths_m[j]
        representative_particles_per_cell = (
            c.alpha_L * cell_volume_m3 / geometry.total_volume
        )
        X2_bulk = x2_particle[j]
        # X1_profile/condensed_profile_kg_s are already the FINAL outer iteration's
        # own converged values (step 4.5 above, top-j=0-to-bottom-j=nz-1,
        # same order as this loop) -- no separate recomputation needed here.
        (
            particle_accumulation_per_particle_w,
            particle_sorption_per_particle_w,
            particle_axial_per_particle_w,
            particle_surface_per_particle_w,
        ) = particle_energy_terms_per_particle[j]
        particle_accumulation_w = (
            representative_particles_per_cell
            * particle_accumulation_per_particle_w
        )
        particle_sorption_w = (
            representative_particles_per_cell * particle_sorption_per_particle_w
        )
        particle_axial_w = (
            representative_particles_per_cell * particle_axial_per_particle_w
        )
        particle_surface_w = (
            representative_particles_per_cell * particle_surface_per_particle_w
        )
        particle_energy_residual_w = (
            particle_accumulation_w
            - particle_sorption_w
            - particle_axial_w
            + particle_surface_w
        )
        Tp_surface = pt.outer_layer_value(particles[j].Tp)
        wet_equilibrium_temperature = thermo.dew_point_temperature(
            vapor_m_hex[j] / max(vapor_m_water[j], 1.0e-12),
            c.antoine_water,
            P=c.pressure_pa,
        )
        vapor_interphase_w = (
            c.hQ
            * c.aV
            * c.alpha_L
            * cell_volume_m3
            * (particle_interphase_temperature[j] - vapor_T[j])
        )
        interphase_mismatch_w = vapor_interphase_w - particle_surface_w
        water_in_cell = (
            m_water_bottom_kg_s if j == nz - 1 else vapor_m_water[j + 1]
        )
        hexane_in_cell = (
            m_hex_bottom_kg_s if j == nz - 1 else vapor_m_hex[j + 1]
        )
        temperature_in_cell = vapor_inf.T if j == nz - 1 else vapor_T[j + 1]
        # Match the accepted-state audit: the coupled block carries bound-water
        # enthalpy at the particle surface temperature.
        bound_surface_temperature = (
            Tp_surface
            if thermal_coupling == "block_implicit"
            else (
                min(Tp_surface, wet_equilibrium_temperature)
                if free_water_profile[j] > 1.0e-12
                else Tp_surface
            )
        )
        accepted_vapor_energy_residual_w = vapor_cell_enthalpy_residual_at_state(
            water_in_kg_s=water_in_cell,
            hexane_in_kg_s=hexane_in_cell,
            temperature_in_k=temperature_in_cell,
            water_out_kg_s=vapor_m_water[j],
            hexane_out_kg_s=vapor_m_hex[j],
            temperature_out_k=vapor_T[j],
            bound_water_to_solid_kg_s=(
                bound_water_mass_rate_w_m3[j] * cell_volume_m3
            ),
            free_water_evaporation_kg_s=evaporated_profile_kg_s[j],
            condensed_water_kg_s=condensed_profile_kg_s[j],
            external_duty_w=q_Iv_profile[j] * cell_volume_m3,
            interphase_heat_to_vapor_w=vapor_interphase_w,
            particle_surface_temperature_k=Tp_surface,
            bound_water_surface_temperature_k=bound_surface_temperature,
            free_water_temperature_k=wet_equilibrium_temperature,
            c=c,
        )
        combined_energy_residual_w = (
            accepted_vapor_energy_residual_w
            + particle_energy_residual_w
            + interphase_mismatch_w
        )
        cells.append(
            DCZCellResult(
                vapor_top=VaporState(wV2=vapor_wV2[j], T=vapor_T[j]),
                particle=particles[j],
                X2_bulk=X2_bulk,
                X1_bulk=X1_profile[j],
                free_water_bulk=free_water_profile[j],
                condensed_water_kg_s=condensed_profile_kg_s[j],
                vapor_water_kg_s=vapor_m_water[j],
                vapor_hexane_kg_s=vapor_m_hex[j],
                evaporated_water_kg_s=evaporated_profile_kg_s[j],
                water_evaporation_energy_w=(
                    evaporated_profile_kg_s[j]
                    * thermo.water_phase_change_enthalpy(
                        wet_equilibrium_temperature,
                        c.vapor_enthalpy_ref,
                        c.particle.cp_water_liquid,
                    )
                ),
                vapor_energy_balance_residual_w=accepted_vapor_energy_residual_w,
                particle_energy_accumulation_w=particle_accumulation_w,
                particle_sorption_source_w=particle_sorption_w,
                particle_axial_conduction_source_w=particle_axial_w,
                particle_surface_heat_to_vapor_w=particle_surface_w,
                particle_energy_balance_residual_w=particle_energy_residual_w,
                interphase_energy_mismatch_w=interphase_mismatch_w,
                combined_energy_balance_residual_w=combined_energy_residual_w,
                wet_phase_constraint_active=evaporated_profile_kg_s[j] > PROFILE_TOL_FLOW,
                wet_phase_water_limited=evaporation_limited_cells[j],
                evaporative_pin_energy_removed_w=0.0,
                particle_temperature_clamped=particle_temperature_clamped[j],
            )
        )
    warm = DCZWarmStart(
        bed_height_m=c.bed_height_m,
        particles=tuple(particles),
        vapor_T=tuple(vapor_T),
        vapor_water_kg_s=tuple(vapor_m_water),
        vapor_hexane_kg_s=tuple(vapor_m_hex),
        X1=tuple(X1_profile),
        free_water=tuple(free_water_profile),
        condensed_water_kg_s=tuple(condensed_profile_kg_s),
        evaporated_water_kg_s=tuple(evaporated_profile_kg_s),
        bound_water_mass_rate_w_m3=tuple(bound_water_mass_rate_w_m3),
        dwpg2_dt=tuple(dwpg2_dt_prev),
        temperature_relaxation=temperature_damping.value,
        hexane_relaxation=hexane_damping.value,
        water_relaxation=water_damping.value,
        water_active_set=previous_water_active_set,
        cell_widths_m=cell_widths_m,
    )
    return DCZZoneResult(
        cells=tuple(cells),
        iterations=iterations,
        converged=converged,
        residuals=residuals,
        warm_start=warm,
        coupling_tolerance=coupling_tolerance,
        cell_widths_m=cell_widths_m,
    )


def _x2_so(wpg2: float, T: float, pc: pt.ParticleConstants) -> float:
    """X2,so(wpg2,T) (eq. A.26) -- used only for the zone's own `X2_bulk`
    diagnostic (Fig. 9(a)-style reporting), dropping the free pore-gas
    hexane term (eq. 5's `Ypg2` -- the smallest of the three contributions
    at equilibrium) as a documented reporting simplification."""
    return thermo.gab_hexane_content(wpg2, T, pc.gab) + pc.X3 * thermo.oil_hexane_content(
        wpg2, pc.oil
    )


def _x2_total_inventory_loading(
    particle: pt.ParticleState,
    geometry: pt.ShellGeometry,
    pc: pt.ParticleConstants,
) -> float:
    """Total particle hexane per kg dry matrix on the A.22 inventory basis."""

    pore_per_dry = pc.alpha_pg * pc.rho_pg / (pc.alpha_ps * pc.rho_ps)
    return pt.shell_volumetric_mean(
        tuple(
            _x2_so(activity, temperature, pc) + pore_per_dry * activity
            for activity, temperature in zip(particle.wpg2, particle.Tp)
        ),
        geometry,
    )
