r"""PHY-050 / PHY-051 critical-volume phase ledger for core2.

Reuse decision (documented seam): the porosity-dependent effective-critical-
volume ledger is a frozen-decision *leaf* module — `core.critical_volume` — that
depends only on `math`/`dataclasses`, has no coupling to the legacy transient
solver, and is already qualified (`tests/test_critical_volume.py`). Re-writing it
in core2 would duplicate a large, exact ledger and risk silent divergence, so
core2 reuses it through this single explicit seam. If the legacy package is ever
retired, only this file moves.

Legacy-evidence policy (ground_truth_v2_gates.yaml: automatic_carryover=false):
the reuse is justified by an EXPLICIT equation-equivalence map + a numerical-
appropriateness note (see `docs/CORE2_LEGACY_REUSE_AUDIT.md`). Summary of the
audit:
  - equivalence: nu_c,eff=kappa_v*eps_p/rho_dm,p (PHY-050); X_c=W_h+X_h,l,c+X_h,g,c;
    retained partial volume = W/rho_h,l (PHY-051); X_f split kept separate (PHY-050).
  - numerics OK: exact algebraic closure, reject-not-clamp on overfill (PHY-050),
    machine-precision X_c,ref=0.20 recovery.
  - CAVEAT enforced here: the leaf module PERMITS a nonzero internal free-water
    volume nu_w; the frozen NOMINAL (Packet C, PHY-005/032) requires nu_w=0.
    core2 therefore uses `nominal_fixed_occupancy` / `solve_nominal_ledger`,
    which forbid internal free water. The raw leaf classes remain available for
    oracles/diagnostics only.
  - CAVEAT for the caller (sphere solver): PHY-035/038/050 amendment (A2c) makes
    X_c(theta) exact ONLY at the activation event; after activation the front
    must use the conserved historical wet trace X_h,w^a, never a re-solved
    X_c(theta). This passive ledger cannot enforce that — the sphere solver must.

What core2 adds: the frozen ledger calibrates kappa_v once at the clean soybean
reference density (rho_h,l,ref = 615 kg/m3). The *state-dependent* X_c(theta),
however, must use the LIVE liquid n-hexane density from the Span-Wagner authority
(PHY-049) — `X_c(theta) = kappa_v * epsilon_p * rho_h,l(T,P) / rho_dm,p` in the
clean liquid-saturated limit. This module supplies that state-dependent value on
core2's own property authority, keeping the 615 kg/m3 constant strictly as the
kappa_v anchor (never as a live property).
"""

from __future__ import annotations

# Frozen ledger (re-exported so core2 code imports from within core2).
from dtdc_simulator.core.critical_volume import (  # noqa: F401
    FANER_SOYBEAN_REFERENCE_CRITICAL_HEXANE_LOADING,
    FANER_SOYBEAN_REFERENCE_DRY_MEAL_DENSITY_KG_M3,
    FANER_SOYBEAN_REFERENCE_POROSITY,
    CriticalFixedPhaseOccupancy,
    CriticalPhaseLedger,
    CriticalPhasePartition,
    CriticalVolumeBasis,
    DrySideFrontInventory,
    FrontMobileJumpLedger,
    HighLoadingHexanePartition,
    clean_reference_critical_volume_multiplier,
    clean_wet_critical_hexane_loading_from_skeletal_basis,
    solve_critical_phase_ledger,
    solve_front_mobile_jump_ledger,
    split_attached_interparticle_hexane,
)
from dtdc_simulator.core2.props import hexane as hx

# Frozen soybean reference constants (properties/soybean.yaml).
KAPPA_V_SOYBEAN = 2.6746237675142717
RHO_HEXANE_LIQUID_REF = 615.0          # kg/m3, kappa_v anchor ONLY (PHY-050)


def critical_hexane_loading(
    T: float,
    epsilon_p: float,
    rho_dm_p: float,
    *,
    pressure_pa: float = 101325.0,
    kappa_v: float = KAPPA_V_SOYBEAN,
) -> float:
    """X_c(theta) in the clean liquid-saturated limit, kg n-hexane/kg dry meal.

    Uses the LIVE Span-Wagner liquid n-hexane density rho_h,l(T,P) (PHY-049),
    not the fixed 615 kg/m3 kappa_v anchor:

        X_c(theta) = kappa_v * epsilon_p * rho_h,l(T,P) / rho_dm,p.
    """
    rho_h_l = hx.state_Tp(T, pressure_pa, "liquid").rho_mass
    return kappa_v * epsilon_p * rho_h_l / rho_dm_p


def effective_critical_specific_volume(
    epsilon_p: float, rho_dm_p: float, *, kappa_v: float = KAPPA_V_SOYBEAN
) -> float:
    """nu_c,eff = kappa_v * epsilon_p / rho_dm,p, m3/kg dry meal (PHY-050)."""
    return kappa_v * epsilon_p / rho_dm_p


def nominal_fixed_occupancy(
    *,
    retained_hexane_loading_kg_per_kg_dry_meal: float,
    pore_gas_hexane_loading_kg_per_kg_dry_meal: float,
    pore_gas_specific_volume_m3_per_kg_dry_meal: float,
    liquid_hexane_density_kg_m3: float,
) -> CriticalFixedPhaseOccupancy:
    """Build the FROZEN-NOMINAL fixed occupancy: internal free water nu_w = 0.

    Enforces Packet C / PHY-005 / PHY-032: free water is external, retained
    water carries no invented internal liquid partial volume, and nu_g is one
    shared pore-gas volume. Reaching a nonzero internal free-water volume would
    require a future separately approved topology, so this seam forbids it.
    """
    return CriticalFixedPhaseOccupancy(
        retained_hexane_loading_kg_per_kg_dry_meal=(
            retained_hexane_loading_kg_per_kg_dry_meal
        ),
        pore_gas_hexane_loading_kg_per_kg_dry_meal=(
            pore_gas_hexane_loading_kg_per_kg_dry_meal
        ),
        free_water_loading_kg_per_kg_dry_meal=0.0,          # nu_w = 0 (nominal)
        pore_gas_specific_volume_m3_per_kg_dry_meal=(
            pore_gas_specific_volume_m3_per_kg_dry_meal
        ),
        liquid_hexane_density_kg_m3=liquid_hexane_density_kg_m3,
        liquid_water_density_kg_m3=1.0,                     # unused (nu_w = 0)
    )


def solve_nominal_ledger(
    basis: CriticalVolumeBasis, occupancy: CriticalFixedPhaseOccupancy
) -> CriticalPhaseLedger:
    """Solve the critical ledger, rejecting any nonzero internal free water.

    The frozen nominal (PHY-005/032) has nu_w = 0. Use `nominal_fixed_occupancy`
    to build the input; this guard also catches a raw leaf occupancy that
    violates the nominal.
    """
    if occupancy.free_water_loading_kg_per_kg_dry_meal != 0.0:
        raise ValueError(
            "nominal critical ledger forbids internal free water (PHY-005/032, "
            "nu_w=0); a nonzero internal free-water topology is not frozen"
        )
    return solve_critical_phase_ledger(basis, occupancy)
