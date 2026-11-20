"""Integrated DT (PHZ + FTRZ + DCZ) solve — Coletto, Bandoni & Blanco (2022)
§7.8/Fig. 5 (BuildSpec), the last slice of M2 (BuildSpec §14): "the
tray-by-tray fixed-point sweep with under-relaxation and warm-start."
Standalone, pure, unit-tested — not yet wired into `core/model.py`'s
real-time `step()` (that is M3, BuildSpec §14: transport lag, per-tick
warm-starting from engine state).

KEY DESIGN INSIGHT (found while tracing the actual data dependencies between
`zones/phz.py`, `zones/ftrz.py`, `zones/dcz.py`'s own signatures, not assumed
upfront): PHZ's solid solve depends only on `Q_indirect` and the top feed BC
-- NOT on vapor state (already true of the standalone `phz.py`: vapor is
informational-only there). So PHZ is solved ONCE, up front and is capped by
the physical PREDESOLV hardware boundary. Crossing the receding-front critical
content changes the jacket-driven mechanism *inside* that section; it never
moves direct-steam FTRZ into a sealed PREDESOLV tray. Cold, boiling-point, or
superheated PHZ outlet states are handed to FTRZ at the first countercurrent
tray; no outer-loop dependency is introduced. The *only* two coupling
scalars across the remaining FTRZ/DCZ interface are `DCZZoneResult.vapor_out`
(top face) feeding FTRZ's `vapor_in` (bottom BC), and
`FTRZZoneResult.solid_out.T` (bottom-most cell) feeding DCZ's `T_L_sup`. This
collapses Fig. 5's "tray-by-tray fixed-point sweep" to a two-variable
Gauss-Seidel loop, spanning trays rather than iterating within each one --
simpler than a literal reading might suggest, and still faithful to it
("solve FTRZ before DCZ... repeat until profiles converge").

BOUNDARIES:
- PHZ/FTRZ interface (`L_PHZ`): fixed at the bottom of the last contiguous
  PREDESOLV tray. `X2,cr` is an internal PHZ mechanism switch (constant-rate
  to jacket-driven falling-rate evaporation), not a steam-contact boundary.
- FTRZ/DCZ interface (`L_FTRZ`): solved endogenously by
  `solve_ftrz_zone`. Every resolved outer pass rebuilds the DCZ mesh at that
  current physical boundary and interpolates the complete lagged DCZ state
  onto the new cell centres. Extreme operation therefore cannot retain a
  first-iterate geometry.

DUTY APPORTIONMENT (bed height as a genuine spatial domain across trays, not
one lumped quantity): each real tray's own `Q_indirect_w` is a *total*
wattage (BuildSpec §5.1). Every zone that occupies part of a tray's height
draws the *same* uniform volumetric density `q_Iv = Q_indirect_w/(A_bed*
tray_height)` as the rest of that tray (eq. A.2a's own convention, applied
consistently instead of inventing a second quantity) -- this is what RESOLVES
two prior gaps left open at M2 Phases 2-3:
- FTRZ's own `Q_cond_w` (no formula in the paper, only a qualitative
  description) is now literally the host tray's `q_Iv`.
- DCZ's own `q_Iv_w_m3` (a single scalar in Phase 3, smearing every tray DCZ
  spans into one artificial average) now accepts a per-cell profile
  (`zones/dcz.py`, this phase), built here by mapping each DCZ cell's
  z-position to whichever real tray contains it.
Neither zone's duty is double-counted: a tray split between FTRZ and DCZ
contributes disjoint sub-heights to each, both scaled from the same
uniform density.

BED TRANSPORT COEFFICIENTS (`hQ`, `hM`, `aV`) -- LITERATURE GAP, confirmed
absent (not assumed) by a targeted re-search of every PDF in
`literature_sources/`, including the supplementary material: Coletto's own
Nuε-Reε correlation (eq. B.7) is cited to Faner's unpublished 2008 PhD
thesis, and `aV` (specific interfacial area, eq. A.35) to Rhodes (2008), a
textbook -- neither source, nor a defining formula for either quantity,
appears anywhere in what we have. Closure, confirmed with the user: standard
packed-bed correlations, clearly tagged `[DERIVED]` like every
other genuinely-unrecoverable constant already in this codebase --
```
aV    = 3 / rP                                                [GEOMETRY, spheres]
Reeps = rho_V * uV,superficial * (2*rP) / (mu_V * eps_b)     [DERIVED, voidage-corrected Re]
PrV   = cp_V * mu_V / kV                                     [STD]
```
feeding the EXISTING, unchanged `thermo.nu_from_reynolds`/`hq_from_nu`/
`schmidt_number`/`hm_from_hq` chain (eqs. B.7-B.10). `kV` (bulk interstitial
vapor thermal conductivity) is itself not a field anywhere in the existing
config -- reused here from `ParticleConstants.k_pg` (the pore-gas thermal
conductivity is the same hexane/steam vapor mixture, just resolved at the
particle-pore scale rather than the bed-interstitial scale; a documented,
cheap approximation rather than inventing an unrelated number).

CROSS-CHECKED, NOT UPGRADED (2026-07-14): two further Faner publications were
tracked down and read after this closure shipped (a 2019 J. Food Process Eng.
journal article and a 2006 conference paper) -- neither is the actual 2008
thesis Coletto cites, and neither contains eq. B.7's 0.6949/0.579
coefficients or a Rhodes-type `aV` derivation. A separate, independent
reconstruction attempt (reasoning from Rhodes' own superficial-vs-interstitial
velocity convention and standard packed-bed pressure-drop correlations, not a
transcription of either source) converged on the SAME two formulas already
implemented here (`Reε = ρ·us·dp/(μ·ε)`, its own top-ranked candidate; `aV =
6(1-ε)/dp`, algebraically identical to `3(1-ε)/rP`) -- a reassuring second
line of reasoning, not a primary-source verification. Still tagged
`[DERIVED]`, not upgraded to `[PAPER]`: that reconstruction's own
stated confidence tops out at "moderate" for `Reε`, and it does not resolve
`sorption_C0`/`sorption_C1`.
(2026-09-23, D20 -- A COMMENT ONLY; no arithmetic on this retired path was
changed.) The 2008 thesis IS now held and digitized
(`paper/analysis/datasets/faner2008_thesis_4_5_3/`), and it REFUTES the
`[DERIVED]` closure recorded above. Its printed page 106, Equation 4.23,
reads `Reε = D_EE ρ_g v_g / [μ_g (1 - ε)]` -- divided by the SOLID
fraction, not by `eps_b` -- and its printed page 109, Equation 4.27, prints
the prefactor 0.6941, not 0.6949. The owner ruled packet D20 option 1 on
2026-09-23 and the three LIVE core2 sites were corrected accordingly; see
`docs/GT_PS2_D20_LEVEL0_BUILD_RECORD_2026-09-23.md`. This retired path is
deliberately NOT corrected, so every bed film coefficient it produces
carries the old convention and is 1/0.7898 times the corrected one at the
frozen ε_b = 0.40. A core2 number and a core number of this quantity are
not the same quantity and must not be compared as one. (2026-09-01: the Cardarelli 1998 thesis
sorption-energy chapter is now held and prints the measured q_st data
rather than a fitted pair -- the coefficients stay [DERIVED]; see
docs/GT_PS2_CARDARELLI1998_INTAKE_RECORD_2026-09-01.md.)

AREA-BASIS CORRECTION (2026-07-26): the earlier note above conflated area
per bed volume, `6*(1-eps_b)/dp`, with Coletto's `aV`.  Because Coletto's
source equations multiply `aV` by `alpha_L` explicitly, `aV` must be the
area per particle volume, `6/dp=3/rP`; their product recovers the standard
bed-specific area.  This distinction is required for the particle and vapor
interphase heat terms to cancel in a combined energy ledger.

INTRAPARTICLE VOLUME -- deliberately DCZ-only. PHZ (sensible-heat/surface
-flash limited) and FTRZ (Receding Front explicitly neglects dry-shell mass
-transfer resistance, §2.3) are lumped-particle by Coletto's own zone
philosophy (§7.0); DCZ alone resolves the 12-layer spherical FVM
(`zones/particle.py`) because it is the one zone where intraparticle
diffusion is rate-limiting. No change needed to make this true -- it already
falls out of composing the three zone modules as published.

SPARGE (direct steam) BC -- scoped to the reference geometry: the SPARGE tray
is the bottommost DT tray, so direct steam sets DCZ's own bottom boundary
condition (`vapor_inf`) via a mass/energy mix with whatever "clean" vapor
arrives from below the DT (documented boundary assumption: direct steam
enters as saturated steam at `T_boil_water`, mixed by a mass-weighted
temperature average against the incoming clean vapor -- both streams are
overwhelmingly water vapor at this point, so a shared `cp_V` is a reasonable
simplification, not derived from the paper). `solve_dt` raises if direct
steam is configured on any tray other than the last -- a genuine point-source
term mid-domain is future work, not implemented here.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from collections.abc import Callable
from dataclasses import replace

import numpy as np
from scipy.optimize import brentq

from dtdc_simulator.core.dt_coupling import (
    _AndersonInterfaceMixer,
    _adapt_relaxation,
)
from dtdc_simulator.core.dt_domain import (
    _average_q_Iv_over_depth,
    _build_dcz_domain,
    _tray_q_Iv as _tray_q_Iv,
)
from dtdc_simulator.core.dt_types import (
    DCZTransientSeed,
    DTAxialProfile,
    DTCouplingResiduals,
    DTResult,
    DTSolverConstants,
    DTTray,
    PHZPassResult,
    SolidFeed,
    TraySummary,
    VaporFeed,
    validate_dt_result as validate_dt_result,
)
from dtdc_simulator.core import thermo
from dtdc_simulator.core.critical_volume import (
    clean_wet_critical_hexane_loading_from_skeletal_basis,
)
from dtdc_simulator.core.zones import dcz
from dtdc_simulator.core.zones import dcz_enthalpy
from dtdc_simulator.core.zones import ftrz
from dtdc_simulator.core.zones import particle as pt
from dtdc_simulator.core.zones import phz as phz_mod

# ---------------------------------------------------------------------------
# Bed transport coefficients (hQ, hM, aV) -- see module docstring
# ---------------------------------------------------------------------------


def bed_transport_coefficients(
    u_V_superficial_m_s: float, c: DTSolverConstants
) -> tuple[float, float, float]:
    """(hQ, hM, aV) from mean bed/vapor conditions -- Faner Nu-Re correlation
    (eq. B.7-B.9), closed with the standard packed-bed `aV`/`Reeps` formulas
    documented in this module's docstring (a confirmed literature gap, not an
    assumed one). The release model uses Coletto B.7-B.9 directly, with no
    empirical sweep-arm agitation enhancement.
    """
    r_P = c.particle.r_P
    eps_b = c.ftrz.bed_porosity
    alpha_L = 1.0 - eps_b
    # Coletto's bed sources use ``alpha_L * aV``. Therefore ``aV`` is the
    # external area per particle volume (sphere geometry, 3/r), while
    # ``alpha_L*aV = 3*(1-eps_b)/r`` is the standard area per bed volume.
    # Including alpha_L here as well would double-count bed solids.
    aV = 3.0 / r_P
    Re_eps = c.rho_V * u_V_superficial_m_s * (2.0 * r_P) / (c.mu_V * eps_b)
    Pr_V = c.cp_V * c.mu_V / c.k_V
    Nu_eps = thermo.nu_from_reynolds(Re_eps, Pr_V)
    hQ = thermo.hq_from_nu(Nu_eps, r_P, c.k_V, alpha_V=eps_b, alpha_L=alpha_L)
    Sc = thermo.schmidt_number(c.mu_V, c.rho_V, c.D_HW)
    hM = thermo.hm_from_hq(hQ, c.rho_V, c.cp_V, Pr_V, Sc)
    return hQ, hM, aV


# ---------------------------------------------------------------------------
# PHZ pass (solved once, up front -- see module docstring)
# ---------------------------------------------------------------------------



def _phz_pass(
    trays: list[DTTray],
    solid_feed: SolidFeed,
    vapor_hint: phz_mod.VaporState,
    m_vapor_water_kg_s: float,
    nz_phz: int,
    c: DTSolverConstants,
) -> PHZPassResult:
    predesolv_indices = [i for i, tray in enumerate(trays) if tray.role == "PREDESOLV"]
    if not predesolv_indices:
        raise ValueError("PHZ requires at least one PREDESOLV tray")
    expected_prefix = list(range(predesolv_indices[-1] + 1))
    if predesolv_indices != expected_prefix:
        raise ValueError("PREDESOLV trays must form a contiguous prefix at the top of the DT")

    solid = phz_mod.SolidState(T=solid_feed.T, X2=solid_feed.X2)
    completed: list[phz_mod.PHZTrayResult] = []

    def _x2_cr(T: float) -> float:
        return clean_wet_critical_hexane_loading_from_skeletal_basis(
            particle_porosity=c.particle.alpha_pg,
            dry_solid_skeletal_density_kg_m3=c.particle.rho_ps,
            critical_volume_multiplier=(
                c.particle.critical_volume_multiplier
            ),
            liquid_hexane_density_kg_m3=thermo.rho_hexane_liquid(T),
        )

    x2_cr = _x2_cr(c.phz.T_boil_hexane)
    x2_eq = thermo.x2_equilibrium(
        max(vapor_hint.T, c.phz.T_boil_hexane),
        solid_feed.X3,
        c.ftrz.gab,
        c.ftrz.oil,
        c.particle.alpha_pg,
        c.particle.alpha_ps,
        c.particle.rho_ps,
    )

    for idx in predesolv_indices:
        tray = trays[idx]
        result = phz_mod.solve_phz_tray(
            nz_phz,
            tray.bed_height_m,
            tray.diameter_m,
            tray.Q_indirect_w,
            solid,
            vapor_hint,
            solid_feed.m_dry_kg_s,
            m_vapor_water_kg_s,
            solid_feed.X1,
            solid_feed.X3,
            c.phz,
            X2_critical=x2_cr,
            X2_equilibrium=x2_eq,
        )
        completed.append(result)
        solid = result.solid_out

    # Hardware boundary: direct steam cannot contact the sealed PREDESOLV
    # beds, so FTRZ always begins at the first countercurrent tray.  Reaching
    # X2_cr inside PREDESOLV now changes the local jacket-driven mechanism
    # (constant -> falling rate); it does not move the hardware handoff.
    last_idx = predesolv_indices[-1]
    last_tray = trays[last_idx]
    last_result = completed.pop()
    return PHZPassResult(
        tray_results=tuple(completed),
        boundary_tray_index=last_idx,
        z_star_m=last_tray.bed_height_m,
        boundary_tray_result=last_result,
        exit_state=last_result.solid_out,
        L_PHZ_m=sum(trays[i].bed_height_m for i in predesolv_indices),
    )


# ---------------------------------------------------------------------------
# DCZ domain assembly (per-tray duty apportionment -- see module docstring)
# ---------------------------------------------------------------------------


def solve_dt(
    trays: list[DTTray],
    solid_feed: SolidFeed,
    vapor_feed_below: VaporFeed,
    c: DTSolverConstants,
    nz_phz: int = 20,
    nz_ftrz: int = 20,
    nz_dcz: int = 20,
    outer_relaxation: float = 0.5,
    outer_tol: float = 1.0e-5,
    outer_max_iter: int = 100,
    dcz_inner_max_iter: int = 100,
    dcz_continuation_max_blocks: int = 100,
    warm_start_vapor_in: ftrz.VaporState | None = None,
    warm_start_T_L_sup: float | None = None,
    warm_start_dcz: dcz.DCZWarmStart | None = None,
    residual_log: list[tuple[int, DTCouplingResiduals, float, float, float]] | None = None,
    dcz_particle_max_step_s: float | None = None,
    dcz_axial_cell_fractions: tuple[float, ...] | None = None,
    dcz_thermal_coupling: str = "partitioned",
    outer_acceleration: bool = False,
    progress_callback: Callable[[dict[str, object]], None] | None = None,
) -> DTResult:
    """Solve the integrated DT: PHZ once, then FTRZ<->DCZ Gauss-Seidel (see
    module docstring for the full design). `warm_start_vapor_in`/
    `warm_start_T_L_sup` seed the outer loop's coupling variables and
    `warm_start_dcz` carries the complete converged inner lagged state (M2's own
    stated "with... warm-start" acceptance criterion) -- the hook a future
    real-time wrapper (M3) will drive tick-to-tick; unused here beyond
    speeding up this one steady solve's own convergence.

    ``dcz_particle_max_step_s`` independently limits the particle-scale
    residence-time step inside each DCZ axial cell.  ``None`` preserves the
    historical one-particle-step-per-axial-cell scheme.  A positive limit
    activates equal substeps, decoupling radial transient accuracy from the
    axial mesh without changing material parameters or total residence time.
    ``dcz_thermal_coupling`` selects the legacy partitioned thermal iteration
    or the phase-free sparse block-implicit solve.
    """
    if not trays:
        raise ValueError("solve_dt requires at least one DT tray")
    if solid_feed.X3 < 0.0:
        raise ValueError("feed oil content X3 must be non-negative")

    # Oil is a live feed disturbance, not a frozen material constant. PHZ
    # already consumes ``solid_feed.X3`` directly; give FTRZ's equilibrium
    # relation and DCZ's particle inventory the same current value. ``replace``
    # keeps the caller-owned frozen constants immutable across solves.
    c = replace(
        c,
        ftrz=replace(c.ftrz, X3=solid_feed.X3),
        particle=replace(c.particle, X3=solid_feed.X3),
    )
    diameters = {t.diameter_m for t in trays}
    if len(diameters) != 1:
        raise ValueError("dt_solver currently requires a uniform tray diameter across the DT")
    diameter_m = next(iter(diameters))
    A_bed_m2 = math.pi / 4.0 * diameter_m**2

    sparge_indices = [i for i, t in enumerate(trays) if t.direct_steam_kg_s > 0.0]
    if sparge_indices and sparge_indices != [len(trays) - 1]:
        raise ValueError(
            "direct (sparge) steam is only supported on the bottommost DT tray in "
            "this phase -- see dt_solver.py module docstring"
        )
    m_dir_kg_s = trays[-1].direct_steam_kg_s

    # --- sparge / bottom BC (mass+energy mix, see module docstring) ---
    m_clean_kg_s = vapor_feed_below.m_water_kg_s + vapor_feed_below.m_hex_kg_s
    m_water_bottom = vapor_feed_below.m_water_kg_s + m_dir_kg_s
    m_hex_bottom = vapor_feed_below.m_hex_kg_s
    m_vapor_total_kg_s = m_water_bottom + m_hex_bottom
    if m_vapor_total_kg_s <= 0.0:
        raise ValueError("the DT bottom vapor boundary must contain positive mass flow")
    if m_dir_kg_s > 0.0:
        T_bottom = dcz_enthalpy.adiabatic_vapor_mix_temperature_k(
            (
                (
                    vapor_feed_below.m_water_kg_s,
                    vapor_feed_below.m_hex_kg_s,
                    vapor_feed_below.T,
                ),
                (m_dir_kg_s, 0.0, c.T_direct_steam),
            ),
            c.ftrz.vapor_enthalpy_ref,
        )
    else:
        T_bottom = vapor_feed_below.T
    vapor_inf = dcz.VaporState(wV2=m_hex_bottom / m_vapor_total_kg_s, T=T_bottom)

    # --- bed transport coefficients (shared by FTRZ + DCZ) ---
    u_V_superficial = m_vapor_total_kg_s / (c.rho_V * A_bed_m2)
    hQ, hM, aV = bed_transport_coefficients(u_V_superficial, c)

    # --- PHZ pass (once) ---
    if m_clean_kg_s > 0.0:
        vapor_hint = phz_mod.VaporState(
            wV1=vapor_feed_below.m_water_kg_s / m_clean_kg_s,
            wV2=vapor_feed_below.m_hex_kg_s / m_clean_kg_s,
            T=vapor_feed_below.T,
        )
        phz_hint_water_kg_s = vapor_feed_below.m_water_kg_s
    else:
        # A vapor-sealed DT->DC rotary airlock implies no external return
        # stream below the sparge tray. The real sparge steam still supplies a
        # finite lower vapor boundary, so use that physical mixture as the
        # PHZ's provisional vapor hint rather than requiring a fictitious
        # epsilon-flow solely to avoid division by zero.
        vapor_hint = phz_mod.VaporState(
            wV1=1.0 - vapor_inf.wV2,
            wV2=vapor_inf.wV2,
            T=vapor_inf.T,
        )
        phz_hint_water_kg_s = m_water_bottom
    phz_result = _phz_pass(
        trays,
        solid_feed,
        vapor_hint,
        phz_hint_water_kg_s,
        nz_phz,
        c,
    )

    boundary_tray = trays[phz_result.boundary_tray_index]
    remaining: list[DTTray] = []
    host_remainder_m = boundary_tray.bed_height_m - phz_result.z_star_m
    if host_remainder_m > 1.0e-9:
        frac = (
            host_remainder_m / boundary_tray.bed_height_m
            if boundary_tray.bed_height_m > 0.0
            else 0.0
        )
        remaining.append(
            replace(
                boundary_tray,
                bed_height_m=host_remainder_m,
                Q_indirect_w=boundary_tray.Q_indirect_w * frac,
            )
        )
    remaining.extend(trays[phz_result.boundary_tray_index + 1 :])
    if not remaining:
        raise ValueError(
            "PHZ boundary lands exactly at the last tray's bottom -- no FTRZ/DCZ "
            "domain remains for the given trays/duties"
        )

    def ftrz_q_Iv(length_m: float) -> float:
        return _average_q_Iv_over_depth(remaining, length_m, A_bed_m2)

    X2_sup = phz_result.exit_state.X2

    # --- FTRZ<->DCZ Gauss-Seidel outer loop ---
    vapor_in = warm_start_vapor_in or ftrz.VaporState(
        m_water_kg_s=(1.0 - vapor_inf.wV2) * m_vapor_total_kg_s,
        m_hex_kg_s=vapor_inf.wV2 * m_vapor_total_kg_s,
        T=vapor_inf.T,
    )
    T_L_sup = warm_start_T_L_sup if warm_start_T_L_sup is not None else c.ftrz.T_boil_hexane

    L_FTRZ_m: float | None = None
    dcz_c: dcz.DCZConstants | None = None
    q_Iv_profile: tuple[float, ...] | None = None
    ftrz_result: ftrz.FTRZZoneResult | None = None
    dcz_result: dcz.DCZZoneResult | None = None
    dcz_warm_start = warm_start_dcz
    converged = False
    coupling_residuals = DTCouplingResiduals()
    iterations = 0
    dcz_feed_temperature_k = T_L_sup
    base_relaxation = min(max(outer_relaxation, 1.0e-3), 1.0)
    temperature_relaxation = base_relaxation
    hexane_relaxation = base_relaxation
    water_relaxation = min(base_relaxation, 0.15)
    previous_temperature_residual = math.inf
    previous_hexane_residual = math.inf
    previous_water_residual = math.inf
    previous_raw_L_FTRZ_m: float | None = None
    interface_mixer = _AndersonInterfaceMixer() if outer_acceleration else None
    previous_interface_map_norm = math.inf
    previous_damped_candidate: np.ndarray | None = None
    previous_update_was_accelerated = False

    for iterations in range(1, outer_max_iter + 1):
        ftrz_result = ftrz.solve_ftrz_zone(
            nz=nz_ftrz,
            X2_sup=X2_sup,
            m_dry_kg_s=solid_feed.m_dry_kg_s,
            vapor_in=vapor_in,
            q_Iv_w_m3=ftrz_q_Iv,
            hQ=hQ,
            hM=hM,
            aV_m2_per_m3=aV,
            diameter_m=diameter_m,
            c=c.ftrz,
            X1_sup=solid_feed.X1,
            T_solid_sup=phz_result.exit_state.T,
        )

        raw_L_FTRZ_m = ftrz_result.L_FTRZ_m
        L_FTRZ_m = raw_L_FTRZ_m
        d_length = (
            math.inf if previous_raw_L_FTRZ_m is None else abs(raw_L_FTRZ_m - previous_raw_L_FTRZ_m)
        )
        previous_raw_L_FTRZ_m = raw_L_FTRZ_m
        dcz_c, q_Iv_profile = _build_dcz_domain(
            remaining,
            L_FTRZ_m,
            A_bed_m2,
            hQ,
            hM,
            aV,
            c,
            nz_dcz,
            dcz_axial_cell_fractions,
        )
        dcz_feed_temperature_k = T_L_sup

        # A cap-limited inner result is not a valid outer-map evaluation.
        # Continue it in bounded blocks from the complete lagged state while
        # keeping the FTRZ boundary and DCZ mesh fixed.  Outer iteration
        # accounting therefore means what it says: one count per fully
        # resolved FTRZ<->DCZ map, not one count per arbitrary inner block.
        for dcz_block_index in range(dcz_continuation_max_blocks):
            dcz_progress = None
            if progress_callback is not None:

                def dcz_progress(
                    payload: dict[str, object],
                    *,
                    dt_iteration: int = iterations,
                    block_index: int = dcz_block_index + 1,
                ) -> None:
                    progress_callback(
                        {
                            "dt_outer_iteration": dt_iteration,
                            "dcz_continuation_block": block_index,
                            **payload,
                        }
                    )

            dcz_result = dcz.solve_dcz_zone(
                nz=nz_dcz,
                m_dry_kg_s=solid_feed.m_dry_kg_s,
                m_vapor_kg_s=m_vapor_total_kg_s,
                T_L_sup=T_L_sup,
                vapor_inf=vapor_inf,
                q_Iv_w_m3=q_Iv_profile,
                c=dcz_c,
                X1_in=solid_feed.X1 + ftrz_result.solid_out.X1,
                outer_max_iter=dcz_inner_max_iter,
                outer_relaxation=outer_relaxation,
                warm_start=dcz_warm_start,
                particle_max_step_s=dcz_particle_max_step_s,
                coupling_tolerance=0.1 * outer_tol,
                axial_cell_fractions=dcz_axial_cell_fractions,
                thermal_coupling=dcz_thermal_coupling,
                progress_callback=dcz_progress,
            )
            dcz_warm_start = dcz_result.warm_start
            if progress_callback is not None:
                progress_callback(
                    {
                        "level": "dcz_block_end",
                        "dt_outer_iteration": iterations,
                        "dcz_continuation_block": dcz_block_index + 1,
                        "converged": dcz_result.converged,
                        "iterations": dcz_result.iterations,
                        "maximum_scaled_residual": (
                            dcz_result.residuals.maximum_scaled
                        ),
                    }
                )
            if dcz_result.converged:
                break
        if not dcz_result.converged:
            break

        new_T_L_sup = ftrz_result.solid_out.T
        dcz_top = dcz_result.vapor_out
        # DCZ now marches explicit water and hexane component flows. Use
        # those solved boundary values directly instead of reconstructing a
        # total flow from solid moisture and then splitting it by wV2.
        new_vapor_in = ftrz.VaporState(
            m_water_kg_s=dcz_result.vapor_water_out_kg_s,
            m_hex_kg_s=dcz_result.vapor_hexane_out_kg_s,
            T=dcz_top.T,
        )

        d_T = abs(new_T_L_sup - T_L_sup)
        # `vapor_in`'s own total mass no longer necessarily equals the fixed
        # whole-DT `m_vapor_total_kg_s` once DCZ condensation is subtracted
        # (see `new_vapor_in` above) -- use its own actual total here.
        d_wV2 = abs(
            dcz_top.wV2 - vapor_in.m_hex_kg_s / (vapor_in.m_water_kg_s + vapor_in.m_hex_kg_s)
        )
        d_TV = abs(dcz_top.T - vapor_in.T)
        d_m_water = abs(new_vapor_in.m_water_kg_s - vapor_in.m_water_kg_s)
        d_m_hex = abs(new_vapor_in.m_hex_kg_s - vapor_in.m_hex_kg_s)
        coupling_residuals = DTCouplingResiduals(
            solid_interface_T=d_T,
            vapor_interface_T=d_TV,
            vapor_interface_hexane_fraction=d_wV2,
            vapor_interface_water_flow=d_m_water,
            vapor_interface_hexane_flow=d_m_hex,
            ftrz_length=d_length,
        )
        maximum_scaled_coupling_residual = coupling_residuals.maximum_scaled(outer_tol)

        coupling_state = np.array(
            [
                T_L_sup,
                vapor_in.m_water_kg_s,
                vapor_in.m_hex_kg_s,
                vapor_in.T,
            ],
            dtype=float,
        )
        coupling_map_residual = np.array(
            [
                new_T_L_sup - T_L_sup,
                new_vapor_in.m_water_kg_s - vapor_in.m_water_kg_s,
                new_vapor_in.m_hex_kg_s - vapor_in.m_hex_kg_s,
                new_vapor_in.T - vapor_in.T,
            ],
            dtype=float,
        )
        damped_candidate = coupling_state + np.array(
            [
                temperature_relaxation,
                water_relaxation,
                hexane_relaxation,
                temperature_relaxation,
            ],
            dtype=float,
        ) * coupling_map_residual
        interface_map_norm = float(
            np.max(np.abs(coupling_map_residual))
        )
        rollback_acceleration = (
            previous_update_was_accelerated
            and previous_damped_candidate is not None
            and interface_map_norm > 1.2 * previous_interface_map_norm
        )
        accelerated = False
        accepted_candidate = damped_candidate
        if rollback_acceleration:
            # The full nonlinear map, evaluated on this iteration, overrules
            # the preceding multisecant prediction. Return to the damped point
            # saved at that same preceding state and discard the invalid local
            # secant model. One map evaluation is spent on the safeguard, but
            # the iteration cannot walk away along a rejected accelerated path.
            accepted_candidate = previous_damped_candidate
            if interface_mixer is not None:
                interface_mixer.clear(cooldown=4)
        elif interface_mixer is not None:
            proposed, accelerated = interface_mixer.propose(
                coupling_state,
                coupling_map_residual,
                damped_candidate,
                tolerance=outer_tol,
            )
            if accelerated:
                # Backtrack along the multisecant direction to remain strictly
                # inside the physical component/temperature domain. Scaling
                # the common direction preserves its coupled character.
                proposal_step = proposed - coupling_state
                lower = np.array([250.0, 0.0, 0.0, 250.0])
                upper = np.array([500.0, math.inf, math.inf, 500.0])
                physical_scale = 1.0
                for index in range(len(proposal_step)):
                    if proposal_step[index] < 0.0:
                        physical_scale = min(
                            physical_scale,
                            0.9
                            * (coupling_state[index] - lower[index])
                            / -proposal_step[index],
                        )
                    elif (
                        proposal_step[index] > 0.0
                        and math.isfinite(upper[index])
                    ):
                        physical_scale = min(
                            physical_scale,
                            0.9
                            * (upper[index] - coupling_state[index])
                            / proposal_step[index],
                        )
                proposed = coupling_state + max(
                    min(physical_scale, 1.0),
                    0.0,
                ) * proposal_step
            if (
                250.0 <= proposed[0] <= 500.0
                and proposed[1] >= 0.0
                and proposed[2] >= 0.0
                and 250.0 <= proposed[3] <= 500.0
            ):
                accepted_candidate = proposed
            else:
                accelerated = False
        previous_interface_map_norm = interface_map_norm
        previous_damped_candidate = (
            damped_candidate.copy() if accelerated else None
        )
        previous_update_was_accelerated = accelerated
        T_L_sup = float(accepted_candidate[0])
        vapor_in = ftrz.VaporState(
            m_water_kg_s=float(accepted_candidate[1]),
            m_hex_kg_s=float(accepted_candidate[2]),
            T=float(accepted_candidate[3]),
        )
        temperature_residual = max(d_T, d_TV)
        hexane_residual = max(d_wV2, d_m_hex)
        water_residual = d_m_water
        temperature_relaxation = _adapt_relaxation(
            temperature_relaxation,
            temperature_residual,
            previous_temperature_residual,
            minimum=0.05,
            maximum=base_relaxation,
        )
        hexane_relaxation = _adapt_relaxation(
            hexane_relaxation,
            hexane_residual,
            previous_hexane_residual,
            minimum=0.05,
            maximum=base_relaxation,
        )
        water_relaxation = _adapt_relaxation(
            water_relaxation,
            water_residual,
            previous_water_residual,
            minimum=0.025,
            maximum=min(base_relaxation, 0.25),
        )
        previous_temperature_residual = temperature_residual
        previous_hexane_residual = hexane_residual
        previous_water_residual = water_residual
        if residual_log is not None:
            residual_log.append(
                (
                    iterations,
                    coupling_residuals,
                    temperature_relaxation,
                    hexane_relaxation,
                    water_relaxation,
                )
            )
        if progress_callback is not None:
            progress_callback(
                {
                    "level": "dt_outer_iteration",
                    "iteration": iterations,
                    "maximum_scaled_residual": maximum_scaled_coupling_residual,
                    "solid_interface_temperature_residual_k": d_T,
                    "vapor_interface_temperature_residual_k": d_TV,
                    "vapor_interface_hexane_fraction_residual": d_wV2,
                    "vapor_interface_water_residual_kg_s": d_m_water,
                    "vapor_interface_hexane_residual_kg_s": d_m_hex,
                    "ftrz_length_residual_m": d_length,
                    "temperature_relaxation": temperature_relaxation,
                    "hexane_relaxation": hexane_relaxation,
                    "water_relaxation": water_relaxation,
                    "anderson_accelerated": accelerated,
                    "anderson_rollback": rollback_acceleration,
                }
            )

        # The dimensionless metric makes every physical tolerance explicit.
        # Its scales intentionally reproduce the established gate exactly:
        # outer_tol for temperatures, fraction and component flows, and the
        # historical 1 mm ceiling for FTRZ length.
        if maximum_scaled_coupling_residual <= 1.0 and dcz_result.converged:
            converged = True
            break

    assert ftrz_result is not None and dcz_result is not None and L_FTRZ_m is not None

    # --- per-real-tray exit (bottom-face) summaries ---
    tray_summaries: list[TraySummary] = []
    for r, t in zip(phz_result.tray_results, trays[: phz_result.boundary_tray_index]):
        tray_summaries.append(
            TraySummary(id=t.id, T=r.solid_out.T, X1=solid_feed.X1, X2=r.solid_out.X2)
        )
    boundary_tray = trays[phz_result.boundary_tray_index]
    boundary_is_full_tray = abs(phz_result.z_star_m - boundary_tray.bed_height_m) <= 1.0e-9
    if boundary_is_full_tray:
        tray_summaries.append(
            TraySummary(
                id=boundary_tray.id,
                T=phz_result.exit_state.T,
                X1=solid_feed.X1,
                X2=phz_result.exit_state.X2,
            )
        )

    # X1 is carried unchanged through PHZ (its own scope excludes a moisture
    # balance) and FTRZ's own solid_out.X1 is a DELTA accumulated from zero
    # (see module docstring) -- the true absolute moisture entering DCZ is
    # their sum (this is what `X1_in` above already fed into the final
    # `dcz_result`). DCZ now DOES carry its own moisture state (see
    # zones/dcz.py's module docstring, "MOISTURE (H2O) BALANCE") -- each
    # DCZ-spanned tray below uses that cell's own accumulated `X1_bulk`
    # rather than this entry value carried forward unchanged.

    geometry = pt.build_shell_geometry(c.particle.r_P, c.particle.Np)
    dcz_cell_widths = (
        dcz_result.cell_widths_m
        if dcz_result.cell_widths_m
        else tuple(dcz_c.bed_height_m / nz_dcz for _ in range(nz_dcz))
    )
    dcz_faces = [0.0]
    for width in dcz_cell_widths:
        dcz_faces.append(dcz_faces[-1] + width)

    def _dcz_cell_at_bottom_face(z_local_m: float) -> dcz.DCZCellResult:
        # `z_local_m` is a tray's own bottom-face position in DCZ's local
        # z-coordinate (0 at the FTRZ/DCZ interface); the small epsilon keeps
        # a boundary landing exactly on a cell edge in the tray whose bottom
        # face it actually is, not the next tray's first cell.
        j = bisect_right(dcz_faces, z_local_m - 1.0e-9) - 1
        j = min(max(j, 0), nz_dcz - 1)
        return dcz_result.cells[j]

    def _ftrz_cell_at_bottom_face(depth_m: float) -> ftrz.FTRZCellResult:
        z = 0.0
        for cell in ftrz_result.cells:
            z += cell.dz_m
            if depth_m <= z + 1.0e-9:
                return cell
        return ftrz_result.cells[-1]

    # A thin free-boundary zone may cross a real-tray boundary.  Sample every
    # tray's physical bottom face from whichever zone actually contains it,
    # rather than assuming FTRZ fits wholly inside ``remaining[0]``.
    depth_from_ftrz_top_m = 0.0
    for tray in remaining:
        depth_from_ftrz_top_m += tray.bed_height_m
        if depth_from_ftrz_top_m <= L_FTRZ_m + 1.0e-9:
            cell = _ftrz_cell_at_bottom_face(depth_from_ftrz_top_m)
            tray_summaries.append(
                TraySummary(
                    id=tray.id,
                    T=cell.solid.T,
                    X1=solid_feed.X1 + cell.solid.X1,
                    X2=cell.solid.X2,
                )
            )
            continue
        cell = _dcz_cell_at_bottom_face(depth_from_ftrz_top_m - L_FTRZ_m)
        tray_summaries.append(
            TraySummary(
                id=tray.id,
                T=dcz.bulk_temperature(cell, geometry),
                X1=cell.X1_bulk,
                X2=cell.X2_bulk,
            )
        )

    # --- whole-DT axial profile (visualization only, see DTAxialProfile docstring) ---
    def _remaining_stage_id_at(depth_from_ftrz_top_m: float) -> str:
        z = 0.0
        for tray in remaining:
            if z - 1.0e-9 <= depth_from_ftrz_top_m < z + tray.bed_height_m + 1.0e-9:
                return tray.id
            z += tray.bed_height_m
        return remaining[-1].id

    def _dcz_stage_id_at(z_local_m: float) -> str:
        """Map a DCZ-local position back to its real tray."""
        return _remaining_stage_id_at(L_FTRZ_m + z_local_m)

    prof_z: list[float] = []
    prof_zone: list[str] = []
    prof_stage: list[str] = []
    prof_solid_T: list[float] = []
    prof_solid_X1: list[float] = []
    prof_solid_X2: list[float] = []
    prof_vapor_T: list[float] = []
    prof_vapor_flow: list[float] = []
    prof_vapor_hex: list[float] = []
    prof_vapor_water: list[float] = []
    prof_mechanism: list[str] = []

    def _append_profile_point(
        z: float,
        zone_name: str,
        stage: str,
        s_T: float,
        s_X1: float,
        s_X2: float,
        v_T: float,
        v_flow: float,
        v_hex: float,
        v_water: float,
        mechanism: str,
    ) -> None:
        prof_z.append(z)
        prof_zone.append(zone_name)
        prof_stage.append(stage)
        prof_solid_T.append(s_T)
        prof_solid_X1.append(s_X1)
        prof_solid_X2.append(s_X2)
        prof_vapor_T.append(v_T)
        prof_vapor_flow.append(v_flow)
        prof_vapor_hex.append(v_hex)
        prof_vapor_water.append(v_water)
        prof_mechanism.append(mechanism)

    # PHZ vapor -- PHYSICAL (not the old per-cell "informational" placeholder).
    # PHZ's SOLID solve is legitimately vapor-decoupled (the bed is jacket-heated,
    # so its temperature/hexane profile doesn't depend on the gas). But the vapor
    # RISING through PHZ physically carries hexane: it arrives hexane-rich from the
    # FTRZ flash front below (`ftrz_result.vapor_out`) and picks up MORE from each
    # PHZ cell's own surface evaporation (`hexane_evaporated_kg_s`, already computed
    # by `solve_phz_tray`) on the way up to the DT vapor outlet. The old trace used
    # each PHZ cell's decoupled `vapor_out` (which read ~pure water and did NOT
    # carry the solvent leaving the top) -- a display artifact, fixed here.
    #
    # Construction: water-vapor flow is CONSTANT across PHZ (X1 is carried unchanged
    # -- no water exchange with the solid here), equal to the FTRZ outlet water flow.
    # Hexane flow leaving a cell's top face = FTRZ hexane + every PHZ cell at or below
    # it (all cells the rising vapor has already passed), so we accumulate bottom-to-top.
    # The vapor TEMPERATURE keeps each cell's own local estimate (`vapor_out.T`, which
    # tracks the solid temperature) -- a rigorous PHZ vapor energy balance is out of
    # scope and unneeded, since the solid solve doesn't consume it.
    phz_vapor_in = ftrz_result.vapor_out  # hexane-rich vapor entering PHZ from the flash front
    m_water_phz = phz_vapor_in.m_water_kg_s  # constant across PHZ
    m_hex_ftrz = phz_vapor_in.m_hex_kg_s

    phz_points: list[tuple[float, str, phz_mod.PHZCellResult]] = []
    z_cum = 0.0
    for tray, result in zip(trays[: phz_result.boundary_tray_index], phz_result.tray_results):
        for cell in result.cells:
            phz_points.append((z_cum + cell.z_from_top_m, tray.id, cell))
        z_cum += tray.bed_height_m
    boundary_tray = trays[phz_result.boundary_tray_index]
    for cell in phz_result.boundary_tray_result.cells:
        phz_points.append((z_cum + cell.z_from_top_m, boundary_tray.id, cell))

    # Suffix-sum the hexane evaporated: cell i's top-face vapor carries the FTRZ
    # hexane plus cells i..last (bottom-to-top accumulation in the vapor's direction).
    hex_flow_at = [0.0] * len(phz_points)
    running_hex = m_hex_ftrz
    for i in range(len(phz_points) - 1, -1, -1):
        running_hex += phz_points[i][2].hexane_evaporated_kg_s
        hex_flow_at[i] = running_hex

    antoine_water = c.ftrz.antoine_water
    antoine_hexane = c.antoine_hexane
    P_atm = thermo.ATM_PRESSURE_PA

    def _binary_dew_T(m_hex: float, m_water: float, fallback: float) -> float:
        """DT vapor temperature = the BINARY (hexane + water) dew point of the
        exiting vapor: the highest temperature at which either component starts
        condensing (its partial pressure reaches saturation) -- the physical
        temperature of the vapor as it leaves the bed. For water-rich vapor the
        water dew point governs (the AOCS/Kemper 2019 dome relation: more water ->
        higher dome); for hexane-rich vapor the hexane dew point governs, so the
        result naturally floors near the ~61 C hexane-water heteroazeotrope
        instead of dropping unphysically toward the water dew point of a nearly
        dry-of-water vapor."""
        n_hex = m_hex / thermo.M_HEXANE
        n_water = m_water / thermo.M_WATER
        n_tot = n_hex + n_water
        if n_tot <= 0.0:
            return fallback
        dews = []
        if m_water > 1.0e-9:
            y_w = n_water / n_tot
            dews.append(
                brentq(
                    lambda T: thermo.antoine_pressure_pa(T, antoine_water) - y_w * P_atm,
                    230.0,
                    470.0,
                )
            )
        if m_hex > 1.0e-9:
            y_h = n_hex / n_tot
            dews.append(
                brentq(
                    lambda T: thermo.antoine_pressure_pa(T, antoine_hexane) - y_h * P_atm,
                    230.0,
                    470.0,
                )
            )
        return max(dews) if dews else fallback

    for i, (z, stage, cell) in enumerate(phz_points):
        m_hex = hex_flow_at[i]
        v_flow = m_water_phz + m_hex
        # PHZ vapor TEMPERATURE is a vapor-liquid-equilibrium quantity (the binary
        # dew point above), NOT the bed solid temperature. The PHZ solid solve is
        # legitimately vapor-decoupled (jacket-heated), but the DT DOME temperature
        # (this trace's top point) is exactly the AOCS/Kemper (2019) dome relation.
        # Reporting the solid temp instead read the top-of-bed feed temperature
        # (~59 C, BELOW the ~61 C heteroazeotrope floor -- physically impossible)
        # as the "dome". The dew point connects continuously to the FTRZ vapor
        # below (whose own T equals its dew point), so no per-tray "sawtooth".
        v_T = _binary_dew_T(m_hex, m_water_phz, cell.solid_out.T)
        _append_profile_point(
            z,
            "PHZ",
            stage,
            cell.solid_out.T,
            solid_feed.X1,
            cell.solid_out.X2,
            v_T,
            v_flow,
            m_hex / v_flow if v_flow > 0.0 else 0.0,
            m_water_phz / v_flow if v_flow > 0.0 else 1.0,
            cell.regime,
        )
    z_cum = phz_result.L_PHZ_m  # authoritative, avoids float drift from the sum above

    # FTRZ may cross a real-tray boundary. NOTE `cell.solid.X1` is the net
    # water transferred within FTRZ only (bulk condensation plus finite-rate
    # sorption, initialized at 0 and accumulated top-to-bottom), so the
    # displayed TOTAL adds the feed moisture carried through PHZ.
    ftrz_depth_m = 0.0
    for cell in ftrz_result.cells:
        z_cum += cell.dz_m
        ftrz_depth_m += cell.dz_m
        v_flow = cell.vapor_out.m_water_kg_s + cell.vapor_out.m_hex_kg_s
        _append_profile_point(
            z_cum,
            "FTRZ",
            _remaining_stage_id_at(ftrz_depth_m),
            cell.solid.T,
            solid_feed.X1 + cell.solid.X1,
            cell.solid.X2,
            cell.vapor_out.T,
            v_flow,
            cell.vapor_out.m_hex_kg_s / v_flow,
            cell.vapor_out.m_water_kg_s / v_flow,
            "STEAM_FLASH",
        )
    z_cum = phz_result.L_PHZ_m + L_FTRZ_m

    # DCZ: explicit component flows vary cell by cell.
    for j, cell in enumerate(dcz_result.cells):
        z_local_end = dcz_faces[j + 1]
        z_local_center = 0.5 * (dcz_faces[j] + dcz_faces[j + 1])
        stage = _dcz_stage_id_at(z_local_center)
        _append_profile_point(
            z_cum + z_local_end,
            "DCZ",
            stage,
            dcz.bulk_temperature(cell, geometry),
            cell.X1_bulk,
            cell.X2_bulk,
            cell.vapor_top.T,
            cell.vapor_flow_kg_s,
            cell.vapor_top.wV2,
            1.0 - cell.vapor_top.wV2,
            "DIFFUSION_CONTROLLED",
        )

    axial_profile = DTAxialProfile(
        z_m=tuple(prof_z),
        zone=tuple(prof_zone),
        stage_id=tuple(prof_stage),
        solid_T=tuple(prof_solid_T),
        solid_X1=tuple(prof_solid_X1),
        solid_X2=tuple(prof_solid_X2),
        vapor_T=tuple(prof_vapor_T),
        vapor_flow_kg_s=tuple(prof_vapor_flow),
        vapor_hexane_frac=tuple(prof_vapor_hex),
        vapor_water_frac=tuple(prof_vapor_water),
        mechanism=tuple(prof_mechanism),
    )
    dcz_feed_moisture = solid_feed.X1 + ftrz_result.solid_out.X1
    saturated_bound_water_capacity = thermo.luikov_equilibrium_moisture(
        thermo.LUIKOV_MAX_VALIDATED_UR,
        dcz_c.luikov,
    )
    dcz_feed_free_water = max(
        dcz_feed_moisture - saturated_bound_water_capacity,
        0.0,
    )
    transient_seed = DCZTransientSeed(
        constants=dcz_c,
        dry_solid_flow_kg_s=solid_feed.m_dry_kg_s,
        solid_feed_particle=pt.ParticleState(
            wpg2=tuple(1.0 for _ in range(c.particle.Np)),
            Tp=tuple(dcz_feed_temperature_k for _ in range(c.particle.Np)),
        ),
        solid_feed_bound_water_kg_kg_dry=(
            dcz_feed_moisture - dcz_feed_free_water
        ),
        solid_feed_free_water_kg_kg_dry=dcz_feed_free_water,
        bottom_vapor_water_kg_s=m_water_bottom,
        bottom_vapor_hexane_kg_s=m_hex_bottom,
        bottom_vapor_temperature_k=T_bottom,
        indirect_duty_by_cell_w=tuple(
            density
            * (math.pi / 4.0 * dcz_c.diameter_m**2)
            * width
            for density, width in zip(q_Iv_profile, dcz_cell_widths)
        ),
    )

    return DTResult(
        phz=phz_result,
        ftrz=ftrz_result,
        dcz=dcz_result,
        tray_summaries=tuple(tray_summaries),
        axial_profile=axial_profile,
        L_PHZ_m=phz_result.L_PHZ_m,
        L_FTRZ_m=L_FTRZ_m,
        L_DCZ_m=dcz_c.bed_height_m,
        hQ=hQ,
        hM=hM,
        aV=aV,
        outer_iterations=iterations,
        converged=converged,
        coupling_residuals=coupling_residuals,
        transient_seed=transient_seed,
    )
