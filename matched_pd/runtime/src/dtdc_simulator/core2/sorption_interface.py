r"""T3b — the sub-critical n-hexane closure that couples the tray march to
the certified particle model.

WHAT THIS FIXES.  The T2 film instrument took the surface n-hexane
activity as unity at every loading, because `solve_heteroazeotrope` was
only ever called on its two-liquid branch.  The consequence, measured:
the interface state was BIT-IDENTICAL at loadings 0.3, 1e-2, 2.8e-4 and
2e-9 kg/kg, so n-hexane behaved as free liquid down to exhaustion, the
residual tail was structurally unreachable, and — because `a_h = 1` was
also the stated basis of the old `_BRACKET_HIGH_K = 340.5` — the
interface temperature was capped 38 K below the benchmark meal outlet.
Below the critical loading there is no mobile n-hexane phase at all, so
neither the activity nor the transport resistance of the certified
particle model was reachable from the tray.

THE REGIME MAP, taken from the certified particle model rather than
invented here.  `particle/front.py`'s Faner oracle is defined on
``X_e <= X <= X_c`` and `props/critical_volume` supplies ``X_c(T)``, so
the frozen physics already has exactly three n-hexane regimes:

  A. ``X >= X_c(T)``  FILM_ACTIVE.  Pore volume liquid-saturated plus
     attached external liquid (PHY-046/PHY-050).  A bulk liquid phase
     touches the surface, so ``a_h = 1`` exactly.  Unchanged behaviour.

  B. ``X_e(T) < X < X_c(T)``  RECEDING_FRONT.  A wet core of radius
     ``s = R((X-X_e)/(X_c-X_e))^(1/3)`` (Faner) still pins ``a_h = 1``
     AT THE FRONT, but n-hexane must now diffuse out through the dry
     shell, so the SURFACE activity falls below 1.  This is the
     resistance the dead `film.py` hook ``dry_shell_resistance`` was
     built for and never received.

  C. ``0 < X <= X_e(T)``  DRY_SORBATE.  No liquid n-hexane anywhere.
     Surface escaping tendency is the PHY-043/048 GAB sorbate activity
     ``a_h(W, T, w_o)``, and transport is the PHY-019/020 pore-diffusion
     problem of `particle/dry_shell.py`.  The benchmark's 280 ppm
     residual lives deep inside this regime.

MEASURED, so the reader can see why regime C is where the KPI lives:
``X_e = retained_hexane(a_h=1)`` is 0.0254 kg/kg at 334.5 K and 0.0105
at 373 K, while ``X_c`` is 0.2017 and 0.1888 — i.e. the sorbate holds
only ~5 % of what the saturated pore holds, so the activity stays pinned
at 1 for the first ~20x of the depletion and only then starts to fall.
The benchmark residual 3.4568e-4 kg/kg sits at ``a_h = 0.031`` (334.5 K)
to ``0.074`` (373 K).

A CORRECTION TO THE AUDIT'S REMEDIATION ORDER.  That audit called for
moving the PHY-046 refusal threshold "from 1e-9 to X_c".  Taken
literally that is wrong, and would have converted a large band of
PHYSICALLY VALID states into refusals: between ``X_e`` and ``X_c`` a
liquid n-hexane phase still exists, so ``a_h = 1`` remains correct there
and only the transport term is missing.  The threshold is not a refusal
boundary at all — it is a REGIME SWITCH, and the switch that turns the
activity off sits at ``X_e``, not ``X_c``.  What the 1e-9 constant
really did was let regimes B and C run under regime A's constitutive
law; that is the defect, and it is fixed by dispatching on the regime,
not by refusing earlier.

TRANSPORT.  Two n-hexane sources feed the surface and they are in
PARALLEL, so their conductances add:

  * the wet core through the dry shell, quasi-steady spherical
    conduction, ``k_shell = D_eff * s / (R (R - s))`` [m/s] — Faner's
    exact shell resistance, zero at ``s = 0`` and singular at ``s = R``
    (i.e. film-controlled, which is regime A);
  * the dry shell's own desorbing sorbate, the Glueckauf linear-driving-
    force reduction of the `dry_shell` PDE, ``k_ldf = 5 D_eff / R`` on a
    PORE-GAS DENSITY driving force — the bare pore diffusivity, because
    the flux is ``D_eff grad c`` and only the ACCUMULATION side carries
    the sorption retardation (`ldf_conductance_m_s` derives this).

Both terms therefore use ``D_eff`` on a density driving force, and the
parallel sum is non-singular at both ends and hands over continuously
from B to C.

DISCLOSED, and MEASURED rather than asserted: the Glueckauf 15/R^2
constant assumes a parabolic radial profile, which a strongly nonlinear
GAB isotherm does not produce.  `scripts/qualify_t3b_ldf_against_dry_
shell.py` runs the reduction against the certified radial PDE and
reports the error; the record
`docs/GT_PS2_T3B_SUBCRITICAL_HEXANE_CLOSURE_2026-08-07.md` carries the
numbers.  That qualification earned its keep immediately — it caught
this closure's own first version using the retarded ``D_app`` in the
flux coefficient, a 452 % terminal error.  The regime-B shell term has
no reachable transient oracle and is disclosed unqualified.

NOT REPRESENTED, named rather than hidden:
  * capillary (Kelvin) depression of ``a_h`` in a partly drained pore.
    Order of magnitude at a 10 nm pore and 373 K is ``a_h ~ 0.86``, so
    regime B's ``a_h = 1`` at the front is optimistic by up to ~15 %.
  * intraparticle temperature gradients: the shell is treated as
    isothermal at the interface temperature, which is the scope
    `dry_shell.py` itself freezes.
  * a particle-size distribution: one radius, monodisperse.

physically_qualifying is False throughout.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass

from .particle import dry_shell as _dry_shell
from .props import critical_volume as _critical_volume
from .props import hexane as _hexane
from .props import sorption as _sorption


class SorptionInterfaceError(RuntimeError):
    """A sub-critical closure contract was violated; the state is unusable."""


class HexaneRegime(enum.Enum):
    """Which constitutive law owns the surface n-hexane."""

    FILM_ACTIVE = "FILM_ACTIVE"          # X >= X_c : bulk liquid, a_h = 1
    RECEDING_FRONT = "RECEDING_FRONT"    # X_e < X < X_c : wet core + dry shell
    DRY_SORBATE = "DRY_SORBATE"          # 0 < X <= X_e : GAB sorbate only
    EXHAUSTED = "EXHAUSTED"              # X <= 0 : nothing left to release


@dataclass(frozen=True)
class ParticleSorptionParams:
    """Frozen soybean particle constants used by the sub-critical closure.

    Defaults are the same values the certified modules carry, so this
    closure cannot silently drift from them:
    `particle/dry_shell.DryShellParams` (porosity, dry-meal density,
    pore diffusivity, residual oil) and
    `props/critical_volume.FANER_SOYBEAN_REFERENCE_*`.
    """

    particle_radius_m: float = 0.885e-3
    particle_porosity: float = 0.141           # = DryShellParams.eps_g
    dry_meal_density_kg_m3: float = 1159.65    # = DryShellParams.rho_dm_p
    pore_diffusivity_m2_s: float = 4.0e-10     # = DryShellParams.D_eff
    residual_oil_mass_fraction: float = 0.0195
    gab: _sorption.GabParams = _sorption.GabParams()
    oil: _sorption.OilIsotherm = _sorption.OilIsotherm()

    def dry_shell_params(self) -> _dry_shell.DryShellParams:
        """The certified `dry_shell` parameter object for these constants."""
        return _dry_shell.DryShellParams(
            eps_g=self.particle_porosity,
            rho_dm_p=self.dry_meal_density_kg_m3,
            D_eff=self.pore_diffusivity_m2_s,
            w_o=self.residual_oil_mass_fraction,
            gab=self.gab,
            oil=self.oil,
        )

    @property
    def external_area_per_dry_mass_m2_kg(self) -> float:
        """3/(r_p rho_dm,p) — external particle area per kg dry meal.

        Derived from the PARTICLE, not from the bed void fraction, so a
        bed-scale voidage/bulk-density inconsistency cannot leak in.
        (One such inconsistency exists and is recorded in the T3b record:
        `dtdc_stack`'s void 0.40 with bulk density 600 implies a solid
        fraction of 0.517, not 0.60, against rho_dm,p = 1159.65.)
        """
        return 3.0 / (self.particle_radius_m * self.dry_meal_density_kg_m3)


@dataclass(frozen=True)
class SubcriticalHexaneState:
    """The surface n-hexane state of one layer's representative particle."""

    regime: HexaneRegime
    mean_loading_kg_kg: float
    critical_loading_kg_kg: float       # X_c(T)
    equilibrium_loading_kg_kg: float    # X_e(T)
    front_radius_m: float               # s; R in regime A, 0 in regime C
    surface_activity: float             # a_h at the outer surface
    surface_pore_density_kg_m3: float   # c_g(R)
    saturation_pore_density_kg_m3: float  # c_g^sat(T)
    apparent_diffusivity_m2_s: float    # D_eff / (dC/dc_g)
    shell_conductance_m_s: float        # k_shell (regime B core source)
    ldf_conductance_m_s: float          # k_ldf  (dry-shell sorbate source)
    intraparticle_conductance_m_s: float  # k_shell + k_ldf
    supply_limited_flux_kg_m2_s: float  # max j the particle can deliver

    @property
    def physically_qualifying(self) -> bool:
        return False


# ---------------------------------------------------------------------------
# Regime boundaries — both taken from frozen authorities, neither invented
# ---------------------------------------------------------------------------
def critical_loading_kg_kg(
    temperature_k: float,
    params: ParticleSorptionParams | None = None,
    *,
    pressure_pa: float = 101325.0,
) -> float:
    """X_c(T) — PHY-050 pore-saturated critical loading, on the LIVE
    Span-Wagner liquid density (`props/critical_volume`)."""
    p = params or ParticleSorptionParams()
    return _critical_volume.critical_hexane_loading(
        temperature_k,
        p.particle_porosity,
        p.dry_meal_density_kg_m3,
        pressure_pa=pressure_pa,
    )


def hexane_activity_ceiling(temperature_k: float, pressure_pa: float) -> float:
    """The largest n-hexane activity a bed at ``pressure_pa`` can support.

    The pore n-hexane is a component of the bed gas, so its partial pressure
    can never exceed the TOTAL pressure the bed is declared at.  Since
    ``a_h = p_h / Psat_h(T)`` on the frozen pressure-ratio convention the
    sorption authorities use, the pressure-feasible ceiling is

        a_ceiling(T, P) = min(1, P / Psat_h(T)).

    Below the n-hexane boiling point at ``P`` this is exactly ``1.0`` and
    nothing changes; above it the ceiling falls, and it is the same bound
    `props/sorption.activity_from_retained` already accepts as its optional
    ``a_ceiling`` argument ("bounds the search to a pressure-feasible
    activity (a_h*Psat_h(T) <= P)").  ``Psat_h`` is read from the frozen
    n-hexane Helmholtz authority; no constant is copied here.
    """
    if not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise SorptionInterfaceError(
            f"activity-ceiling temperature must be positive and finite, got {temperature_k}"
        )
    if not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise SorptionInterfaceError(
            f"activity-ceiling pressure must be positive and finite, got {pressure_pa}"
        )
    saturation_pressure_pa = _hexane.saturation_pressure(temperature_k)
    if not math.isfinite(saturation_pressure_pa) or saturation_pressure_pa <= 0.0:
        raise SorptionInterfaceError(
            f"the frozen n-hexane authority returned a non-usable saturation "
            f"pressure {saturation_pressure_pa} at T = {temperature_k} K"
        )
    return min(1.0, pressure_pa / saturation_pressure_pa)


def equilibrium_loading_kg_kg(
    temperature_k: float,
    params: ParticleSorptionParams | None = None,
    *,
    pressure_pa: float = 101325.0,
) -> float:
    """X_e(T, P) — the TOTAL particle loading at the pressure-feasible
    n-hexane activity ceiling with no liquid.

    This is the loading at which the Faner front reaches the centre and
    the last liquid n-hexane disappears; below it the surface activity
    starts to fall.  It must be the same inventory the certified
    `dry_shell.storage` counts, i.e. saturated pore GAS as well as the
    PHY-043/048 GAB sorbate:

        X_e = [eps_g c_g(a*, T) + rho_dm,p W_h(a*, T, w_o)] / rho_dm,p,
        a* = a_ceiling(T, P) = min(1, P / Psat_h(T))

    Dropping the pore-gas term (7.7 % of X_e at 373 K) would put a
    0.03 step in the surface activity at the regime B/C handover.

    F6 (2026-09-02, wave-2 defect #6).  Until this fix the activity was the
    IMPLICIT ``a* = 1``, i.e. a pore n-hexane partial pressure of Psat_h(T)
    in a bed declared at 1.01-1.70 bar.  Above the n-hexane boiling point
    that is not a reachable state, and X_e was overstated accordingly:
    measured 4.24x at 373 K, 11.78x at 423 K and 13.39x at 433 K at 1 atm
    (2.06x / 6.80x / 7.81x at 1.70 bar).  At and below the boiling point the
    ceiling is exactly 1.0 and every number is BIT-IDENTICAL - in particular
    the T3b 334.5 K anchor is unchanged at 1.000x, which
    ``tests/test_core2_codec_caloric_seam_instrument.py`` pins as a
    regression check.  X_e is both the regime switch (`classify`) and the
    Faner front-law denominator (`front_radius_m`,
    `qsc_layer_falling_rate_law.receding_front_fraction`), so the exposure is
    the 373-433 K band the desolventizer trays run in.
    """
    p = params or ParticleSorptionParams()
    ceiling = hexane_activity_ceiling(temperature_k, pressure_pa)
    return (
        _dry_shell.storage(
            ceiling * _dry_shell.cg_saturation(temperature_k),
            temperature_k,
            p.dry_shell_params(),
        )
        / p.dry_meal_density_kg_m3
    )


def classify(
    mean_loading_kg_kg: float,
    temperature_k: float,
    params: ParticleSorptionParams | None = None,
    *,
    pressure_pa: float = 101325.0,
) -> HexaneRegime:
    """Which of the three certified n-hexane regimes this loading is in."""
    if not math.isfinite(mean_loading_kg_kg):
        raise SorptionInterfaceError("n-hexane loading must be finite")
    if mean_loading_kg_kg <= 0.0:
        return HexaneRegime.EXHAUSTED
    p = params or ParticleSorptionParams()
    if mean_loading_kg_kg >= critical_loading_kg_kg(
        temperature_k, p, pressure_pa=pressure_pa
    ):
        return HexaneRegime.FILM_ACTIVE
    if mean_loading_kg_kg > equilibrium_loading_kg_kg(
        temperature_k, p, pressure_pa=pressure_pa
    ):
        return HexaneRegime.RECEDING_FRONT
    return HexaneRegime.DRY_SORBATE


def front_radius_m(
    mean_loading_kg_kg: float,
    temperature_k: float,
    params: ParticleSorptionParams | None = None,
    *,
    pressure_pa: float = 101325.0,
) -> float:
    """Faner front position s(X) — the same cube-root relation
    `particle/front.faner_radius_oracle` certifies, evaluated on the live
    X_c(T)/X_e(T) pair.  R in regime A, 0 at or below X_e."""
    p = params or ParticleSorptionParams()
    X_c = critical_loading_kg_kg(temperature_k, p, pressure_pa=pressure_pa)
    X_e = equilibrium_loading_kg_kg(temperature_k, p, pressure_pa=pressure_pa)
    if not X_c > X_e:
        raise SorptionInterfaceError(
            f"degenerate loading interval at T = {temperature_k} K: "
            f"X_e = {X_e} is not below X_c = {X_c}"
        )
    if mean_loading_kg_kg >= X_c:
        return p.particle_radius_m
    if mean_loading_kg_kg <= X_e:
        return 0.0
    return p.particle_radius_m * (
        (mean_loading_kg_kg - X_e) / (X_c - X_e)
    ) ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# Transport — both conductances on the certified dry-shell authority
# ---------------------------------------------------------------------------
def apparent_diffusivity_m2_s(
    pore_gas_density_kg_m3: float,
    temperature_k: float,
    params: ParticleSorptionParams | None = None,
) -> float:
    """D_app = D_eff / (dC/dc_g) — the SORPTION-RETARDED pore diffusivity.

    ``dC/dc_g`` is `dry_shell.storage_capacity`, the analytic derivative
    of the PHY-020 local-equilibrium storage.  Writing the certified PDE
    ``dC/dt = div(D_eff grad c_g)`` in terms of C gives
    ``dC/dt = div(D_app grad C)``, so D_app is the diffusivity that
    governs how fast the LOADING equilibrates — the quantity a tray
    march needs.  It is not a new constitutive parameter.
    """
    p = params or ParticleSorptionParams()
    capacity = _dry_shell.storage_capacity(
        pore_gas_density_kg_m3, temperature_k, p.dry_shell_params()
    )
    if not math.isfinite(capacity) or capacity <= 0.0:
        raise SorptionInterfaceError(
            f"non-positive dry-shell storage capacity {capacity} at "
            f"c_g = {pore_gas_density_kg_m3}, T = {temperature_k} K"
        )
    return p.pore_diffusivity_m2_s / capacity


def shell_conductance_m_s(
    front_radius: float, params: ParticleSorptionParams | None = None
) -> float:
    """k_shell = D_eff s / (R (R - s)) — Faner's quasi-steady spherical
    dry-shell conductance, per unit EXTERNAL area, on a pore-gas density
    driving force.  Zero at s = 0, +inf at s = R (film-controlled)."""
    p = params or ParticleSorptionParams()
    R = p.particle_radius_m
    if front_radius <= 0.0:
        return 0.0
    if front_radius >= R:
        return math.inf
    return p.pore_diffusivity_m2_s * front_radius / (R * (R - front_radius))


def ldf_conductance_m_s(
    params: ParticleSorptionParams | None = None,
) -> float:
    r"""k_ldf = 5 D_eff / R — the Glueckauf linear-driving-force
    conductance of a sphere, per unit external area, on a PORE-GAS
    DENSITY driving force.

    THE BARE D_eff, NOT THE RETARDED D_app — this is the Kirchhoff form
    and getting it wrong is a 450 % error, measured.  Derivation: the
    certified PDE is ``dC/dt = div(D_eff grad c)``.  Writing it on the
    storage variable gives ``dC/dt = div(D_app grad C)`` with the
    strongly varying ``D_app = D_eff/(dC/dc)``, so a naive LDF on C must
    pick SOME value of D_app.  The correct linearization is the integral
    mean over the driving range, and that integral collapses exactly:

        D_app_bar = 1/(C1-C0) * int_{C0}^{C1} D_app dC
                  = D_eff * (c1 - c0) / (C1 - C0)

    because ``D_app dC = D_eff dc``.  Substituting it back into
    ``dCbar/dt = -(15 D_app_bar/R^2)(Cbar - C_s)`` cancels the storage
    difference and leaves ``dCbar/dt = -(15 D_eff/R^2)(c_bar - c_s)`` —
    linear in the DENSITY difference with the bare pore diffusivity.
    Equivalently: the flux is D_eff grad c, so the flux-side coefficient
    was never allowed to carry the sorption retardation; only the
    accumulation side is retarded, and that is already in ``dC/dt``.

    Measured against the certified radial PDE, this reduces the worst
    terminal error over the qualification cases from 452 % to a few per
    cent — see `scripts/qualify_t3b_ldf_against_dry_shell.py`.  What
    remains is the genuine Glueckauf profile-shape approximation (the
    "15" assumes a parabolic radial profile), which the same script
    reports and which is DISCLOSED, not eliminated.

    `apparent_diffusivity_m2_s` is retained as the correct object for
    the TIME-CONSTANT reading tau = R^2/(15 D_app); it is not the object
    for a flux coefficient.
    """
    p = params or ParticleSorptionParams()
    return 5.0 * p.pore_diffusivity_m2_s / p.particle_radius_m


def _intraparticle_supply(
    mean_loading_kg_kg: float,
    temperature_k: float,
    params: ParticleSorptionParams,
    pressure_pa: float,
) -> tuple[HexaneRegime, float, float, float, float, float, float, float]:
    """(regime, s, c_mean, c_sat, D_app, k_shell, k_ldf, k_total).

    Everything the sub-critical closure needs that does NOT depend on the
    flux — so a caller can ask "how much can this particle deliver?"
    before deciding what flux to impose.
    """
    dsp = params.dry_shell_params()
    c_sat = _dry_shell.cg_saturation(temperature_k)
    regime = classify(
        mean_loading_kg_kg, temperature_k, params, pressure_pa=pressure_pa
    )
    s = front_radius_m(
        mean_loading_kg_kg, temperature_k, params, pressure_pa=pressure_pa
    )
    if regime is HexaneRegime.DRY_SORBATE:
        storage_mean = mean_loading_kg_kg * params.dry_meal_density_kg_m3
        try:
            c_mean = _dry_shell.cg_from_storage(storage_mean, temperature_k, dsp)
        except ValueError as exc:
            raise SorptionInterfaceError(
                f"mean loading {mean_loading_kg_kg} kg/kg exceeds the "
                f"admissible dry-shell storage at T = {temperature_k} K: {exc}"
            ) from exc
    else:
        c_mean = c_sat
    d_app = apparent_diffusivity_m2_s(c_mean, temperature_k, params)
    k_shell = shell_conductance_m_s(s, params)
    k_ldf = ldf_conductance_m_s(params)
    return regime, s, c_mean, c_sat, d_app, k_shell, k_ldf, k_shell + k_ldf


def supply_limited_flux_kg_m2_s(
    mean_loading_kg_kg: float,
    temperature_k: float,
    params: ParticleSorptionParams | None = None,
    *,
    pressure_pa: float = 101325.0,
) -> float:
    """The largest n-hexane flux this particle can deliver to its own
    surface, per unit external area: ``k_intra * c_mean``.

    Reached when the surface pore density is driven to zero.  In regime A
    it is infinite (a bulk liquid phase touches the surface).  This is the
    quantity that makes the residual tail a PLATEAU rather than a strip to
    zero, and callers use it to bound their flux estimator explicitly
    rather than discovering the bound as a refusal.
    """
    p = params or ParticleSorptionParams()
    if classify(
        mean_loading_kg_kg, temperature_k, p, pressure_pa=pressure_pa
    ) is HexaneRegime.FILM_ACTIVE:
        return math.inf
    _, _, c_mean, _, _, _, _, k_total = _intraparticle_supply(
        mean_loading_kg_kg, temperature_k, p, pressure_pa
    )
    return k_total * c_mean


def supply_limited_loading_drop(
    mean_loading_kg_kg: float,
    temperature_k: float,
    step_duration_s: float,
    params: ParticleSorptionParams | None = None,
    *,
    pressure_pa: float = 101325.0,
    substeps: int = 8,
) -> float:
    """The largest LOADING DROP the particle's own transport permits over
    a macro-step of `step_duration_s`, kg/kg dry.

    WHY THIS EXISTS AND THE INSTANTANEOUS FLUX BOUND IS NOT ENOUGH.  The
    intraparticle conductance changes by orders of magnitude as the front
    recedes, so a start-of-step value is the wrong bound in both
    directions.  In regime A it is INFINITE (liquid at the surface), and
    a tray march using it lets a plant-scale area strip a layer from
    X = 0.335 straight to zero inside one macro-step — skipping regimes B
    and C entirely.  That is a time-resolution failure, not a physics
    one, and the cure is to integrate the supply across the step while
    the regime traverses.

    REGIME-CROSSING LIMITER.  Above X_e the particle imposes no bound
    that can be integrated: regime A has bulk liquid at the surface, and
    regime B's shell conductance is correctly SINGULAR as s -> R (a
    zero-thickness shell has zero resistance), so the front leaves the
    surface under any finite driving force and the release there is
    film-limited, not particle-limited.  The bound above X_e is therefore
    the crossing itself — one macro-step may take the loading down to the
    entry of the transport-limited regime and no further, with the
    remainder resolved on the next step.  That is ordinary operator
    splitting, and it is what stops a plant-scale area from stripping a
    layer from feed loading to zero inside one step.

    At and below X_e the rate is smooth (s = 0, k = 5 D_eff/R, no
    singularity) and is integrated with RK4.

    DISCLOSED, and deliberately in the conservative direction: the bound
    is evaluated at a bulk pore-gas density of ZERO, i.e. the most
    favourable gas side the layer could possibly present, and it ignores
    fresh feed arriving mid-step.  It is a CEILING on release, never a
    rate law — the rate law is the film.
    """
    p = params or ParticleSorptionParams()
    if not math.isfinite(step_duration_s) or step_duration_s <= 0.0:
        raise SorptionInterfaceError("step duration must be positive")
    if mean_loading_kg_kg <= 0.0:
        return 0.0
    area_per_mass = p.external_area_per_dry_mass_m2_kg
    X_e = equilibrium_loading_kg_kg(temperature_k, p, pressure_pa=pressure_pa)
    if mean_loading_kg_kg > X_e:
        return mean_loading_kg_kg - X_e

    def rate(value: float) -> float:
        """d(loading)/dt ceiling, kg/(kg dry . s); negative = leaving."""
        if value <= 0.0:
            return 0.0
        _, _, c_mean, _, _, _, _, k_total = _intraparticle_supply(
            value, temperature_k, p, pressure_pa
        )
        return -area_per_mass * k_total * c_mean

    x = mean_loading_kg_kg
    h = step_duration_s / substeps
    for _ in range(substeps):
        k1 = rate(x)
        k2 = rate(max(x + 0.5 * h * k1, 0.0))
        k3 = rate(max(x + 0.5 * h * k2, 0.0))
        k4 = rate(max(x + h * k3, 0.0))
        x = max(x + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0, 0.0)
    return mean_loading_kg_kg - x


def surface_state_from_gas(
    *,
    mean_loading_kg_kg: float,
    temperature_k: float,
    gas_hexane_density_kg_m3: float,
    film_conductance_m_s: float,
    params: ParticleSorptionParams | None = None,
    pressure_pa: float = 101325.0,
) -> SubcriticalHexaneState:
    r"""Surface state from the SERIES resistance, solved in closed form.

    PREFER THIS over `surface_state` in a marched layer.  Equating the
    two transport paths at the surface,

        j = k_gas (c_s - c_bulk)          [external film]
        j = k_intra (c_mean - c_s)        [intraparticle]

    gives directly

        c_s = (k_gas c_bulk + k_intra c_mean) / (k_gas + k_intra),

    with no iteration and no lagged flux.

    WHY THE LAGGED-FLUX FORM OSCILLATES AND THIS DOES NOT.  Feeding the
    previous step's flux into the activity closes a loop whose gain is
    ``d(ln j)/d(ln a_h) * d(ln a_h)/d(ln j) = (c_mean - c_s)/c_s``.  At
    the residual plateau the surface is depleted to a few per cent of the
    mean, so that gain is order 50 and no reasonable under-relaxation
    stabilises it — the tray march simply fails to converge.  Solving the
    balance instead of lagging it removes the loop entirely.

    The physical content is worth stating plainly: with a mass Biot
    number ``k_gas R / D_eff`` of order 1e5, ``k_intra << k_gas`` and the
    formula collapses to ``c_s -> c_bulk``.  The particle surface sits in
    equilibrium with the local gas, and the rate is set wholly by the
    solid side.  That is the regime the benchmark's own QG-5 mechanism
    row demands ("local-equilibrium storage plus intraparticle
    transport", Cardarelli refs 9/10) and is precisely what the
    film-only model could not represent.
    """
    p = params or ParticleSorptionParams()
    T = temperature_k
    if not math.isfinite(T) or T <= 0.0:
        raise SorptionInterfaceError(f"temperature must be positive, got {T}")
    if not math.isfinite(film_conductance_m_s) or film_conductance_m_s <= 0.0:
        raise SorptionInterfaceError(
            f"film conductance must be positive and finite, got "
            f"{film_conductance_m_s}"
        )
    if not math.isfinite(gas_hexane_density_kg_m3) or gas_hexane_density_kg_m3 < 0.0:
        raise SorptionInterfaceError(
            f"bulk n-hexane density must be finite and non-negative, got "
            f"{gas_hexane_density_kg_m3}"
        )
    c_sat = _dry_shell.cg_saturation(T)
    X_c = critical_loading_kg_kg(T, p, pressure_pa=pressure_pa)
    X_e = equilibrium_loading_kg_kg(T, p, pressure_pa=pressure_pa)
    regime = classify(mean_loading_kg_kg, T, p, pressure_pa=pressure_pa)

    if regime is HexaneRegime.EXHAUSTED:
        raise SorptionInterfaceError(
            f"n-hexane loading {mean_loading_kg_kg} is not positive; there is "
            "no surface state to report"
        )
    if regime is HexaneRegime.FILM_ACTIVE:
        d_app = apparent_diffusivity_m2_s(c_sat, T, p)
        return SubcriticalHexaneState(
            regime=regime,
            mean_loading_kg_kg=mean_loading_kg_kg,
            critical_loading_kg_kg=X_c,
            equilibrium_loading_kg_kg=X_e,
            front_radius_m=p.particle_radius_m,
            surface_activity=1.0,
            surface_pore_density_kg_m3=c_sat,
            saturation_pore_density_kg_m3=c_sat,
            apparent_diffusivity_m2_s=d_app,
            shell_conductance_m_s=math.inf,
            ldf_conductance_m_s=ldf_conductance_m_s(p),
            intraparticle_conductance_m_s=math.inf,
            supply_limited_flux_kg_m2_s=math.inf,
        )

    regime, s, c_mean, c_sat, d_app, k_shell, k_ldf, k_total = (
        _intraparticle_supply(mean_loading_kg_kg, T, p, pressure_pa)
    )
    if not math.isfinite(k_total):
        # s -> R: a zero-thickness shell has no resistance, so the
        # surface is at the front's saturated state.
        c_surface = c_sat
    else:
        c_surface = (
            film_conductance_m_s * gas_hexane_density_kg_m3
            + k_total * c_mean
        ) / (film_conductance_m_s + k_total)
    if regime is HexaneRegime.RECEDING_FRONT and c_surface > c_sat:
        # A WET FRONT CANNOT SUPPORT A SUPERSATURATED SURFACE.  Liquid
        # n-hexane is present at s, so any excess above c_g^sat condenses
        # and the surface is pinned at saturation.  That is a PHASE
        # CONSTRAINT, not a clamp on the activity: the liquid phase which
        # does the pinning exists in this regime and is exactly what
        # `RECEDING_FRONT` means.
        #
        # It binds in a narrow band just below X_c, where the shell is
        # thin and the series formula's k_intra >> k_gas limit rounds the
        # surface a hair above saturation.  Measured before this fix: a
        # plant-scale fixture refused at X = 1.998444e-1 against
        # X_c(337.23 K) = 2.004e-1 with an activity of 1.001257 — a 0.13 %
        # overshoot at the regime boundary, not a physical excursion.
        c_surface = c_sat
    activity = c_surface / c_sat
    if activity > 1.0:
        # DRY_SORBATE: there is no liquid phase anywhere to absorb the
        # excess, so a supersaturated surface would have to create one.
        # That is outside this regime — reject, do not clamp.
        raise SorptionInterfaceError(
            f"surface n-hexane activity {activity:.6f} exceeds 1 at "
            f"X = {mean_loading_kg_kg:.6e} kg/kg, T = {T:.2f} K with a bulk "
            f"n-hexane density of {gas_hexane_density_kg_m3:.6e} kg/m3 — the "
            "gas is supersaturated against a surface with no liquid phase "
            "to condense into; reject, do not clamp"
        )
    return SubcriticalHexaneState(
        regime=regime,
        mean_loading_kg_kg=mean_loading_kg_kg,
        critical_loading_kg_kg=X_c,
        equilibrium_loading_kg_kg=X_e,
        front_radius_m=s,
        surface_activity=activity,
        surface_pore_density_kg_m3=c_surface,
        saturation_pore_density_kg_m3=c_sat,
        apparent_diffusivity_m2_s=d_app,
        shell_conductance_m_s=k_shell,
        ldf_conductance_m_s=k_ldf,
        intraparticle_conductance_m_s=k_total,
        supply_limited_flux_kg_m2_s=k_total * max(c_mean - c_surface, 0.0),
    )


# ---------------------------------------------------------------------------
# The closure
# ---------------------------------------------------------------------------
def surface_state(
    *,
    mean_loading_kg_kg: float,
    temperature_k: float,
    hexane_flux_kg_m2_s: float,
    params: ParticleSorptionParams | None = None,
    pressure_pa: float = 101325.0,
) -> SubcriticalHexaneState:
    """Surface n-hexane state for a representative particle carrying
    `mean_loading_kg_kg` and currently releasing `hexane_flux_kg_m2_s`
    per unit external area (positive solid -> gas).

    The flux is the LAGGED value from the previous macro-step, matching
    the depletion sequence's existing explicit-state discipline (PASS 1
    solves at the current lagged layer states).  A negative flux
    (condensation onto the particle) raises the surface above the mean
    and is handled by the same expression.
    """
    p = params or ParticleSorptionParams()
    T = temperature_k
    if not math.isfinite(T) or T <= 0.0:
        raise SorptionInterfaceError(f"temperature must be positive, got {T}")
    c_sat = _dry_shell.cg_saturation(T)
    X_c = critical_loading_kg_kg(T, p, pressure_pa=pressure_pa)
    X_e = equilibrium_loading_kg_kg(T, p, pressure_pa=pressure_pa)
    regime = classify(mean_loading_kg_kg, T, p, pressure_pa=pressure_pa)

    if regime is HexaneRegime.EXHAUSTED:
        raise SorptionInterfaceError(
            f"n-hexane loading {mean_loading_kg_kg} is not positive; there is "
            "no surface state to report"
        )

    if regime is HexaneRegime.FILM_ACTIVE:
        # Bulk liquid touches the surface: a_h = 1 exactly, no
        # intraparticle resistance in the PHY-046 film-active law.
        d_app = apparent_diffusivity_m2_s(c_sat, T, p)
        return SubcriticalHexaneState(
            regime=regime,
            mean_loading_kg_kg=mean_loading_kg_kg,
            critical_loading_kg_kg=X_c,
            equilibrium_loading_kg_kg=X_e,
            front_radius_m=p.particle_radius_m,
            surface_activity=1.0,
            surface_pore_density_kg_m3=c_sat,
            saturation_pore_density_kg_m3=c_sat,
            apparent_diffusivity_m2_s=d_app,
            shell_conductance_m_s=math.inf,
            ldf_conductance_m_s=ldf_conductance_m_s(p),
            intraparticle_conductance_m_s=math.inf,
            supply_limited_flux_kg_m2_s=math.inf,
        )

    # --- regimes B and C share one expression -------------------------
    # `c_mean` is the density driving the intraparticle transport: the
    # dry region's own sorbate inventory in regime C, and the saturated
    # wet front in regime B (whose liquid inventory is carried by the
    # shell conductance term instead).
    regime, s, c_mean, c_sat, d_app, k_shell, k_ldf, k_total = (
        _intraparticle_supply(mean_loading_kg_kg, T, p, pressure_pa)
    )

    if not math.isfinite(k_total) or k_total <= 0.0:
        raise SorptionInterfaceError(
            f"non-positive intraparticle conductance {k_total} at "
            f"X = {mean_loading_kg_kg}, T = {T} K"
        )
    c_surface = c_mean - hexane_flux_kg_m2_s / k_total
    if c_surface < 0.0:
        raise SorptionInterfaceError(
            f"the particle cannot supply {hexane_flux_kg_m2_s:.6e} kg/m2/s at "
            f"X = {mean_loading_kg_kg:.6e} kg/kg, T = {T:.2f} K "
            f"(intraparticle conductance {k_total:.6e} m/s gives a negative "
            "surface pore density) — the caller must cap the flux at "
            "`supply_limited_flux_kg_m2_s`, not extrapolate"
        )

    if regime is HexaneRegime.RECEDING_FRONT and c_surface > c_sat:
        # Same phase constraint as `surface_state_from_gas`: a wet front
        # is present, so excess above c_g^sat condenses and pins the
        # surface at saturation.  Kept identical between the two entry
        # points so they cannot disagree about physics.
        c_surface = c_sat
    activity = c_surface / c_sat
    # The GAB closure is singular at K(T) a_h = 1 and REJECTS beyond it
    # (PHY-021: no activity clamp).  In DRY_SORBATE there is no sorbate
    # branch above unity and no liquid phase to condense into, so a
    # supersaturated surface means the caller handed us a condensing flux
    # this regime cannot represent.
    if activity > 1.0:
        raise SorptionInterfaceError(
            f"surface n-hexane activity {activity:.6f} exceeds 1 at "
            f"X = {mean_loading_kg_kg:.6e} kg/kg, T = {T:.2f} K: a condensing "
            "flux would create a liquid phase this sub-critical regime does "
            "not carry — reject, do not clamp"
        )

    return SubcriticalHexaneState(
        regime=regime,
        mean_loading_kg_kg=mean_loading_kg_kg,
        critical_loading_kg_kg=X_c,
        equilibrium_loading_kg_kg=X_e,
        front_radius_m=s,
        surface_activity=activity,
        surface_pore_density_kg_m3=c_surface,
        saturation_pore_density_kg_m3=c_sat,
        apparent_diffusivity_m2_s=d_app,
        shell_conductance_m_s=k_shell,
        ldf_conductance_m_s=k_ldf,
        intraparticle_conductance_m_s=k_total,
        supply_limited_flux_kg_m2_s=k_total * c_mean,
    )

# ---------------------------------------------------------------------------
# T3b track (a): the WETTED FRACTION of a well-mixed layer
# (design record: docs/GT_PS2_T3B_WETTED_FRACTION_DESIGN_2026-08-08.md)
# ---------------------------------------------------------------------------
def wetted_fraction(
    *,
    upstream_loading_kg_kg: float,
    mean_loading_kg_kg: float,
    critical_loading_kg_kg: float,
) -> float:
    r"""Fraction of a CSTR layer's particles still carrying free liquid.

    THE DEFECT THIS REMOVES.  A lumped layer cannot be PARTLY wetted:
    its surface is film-active (a_h = 1, ~335 K azeotrope) or
    sub-critical (a_h << 1, ~383 K sorbed-water boiling point), and it
    switches binarily at X_c.  Both the rate and a ~48 K interface
    temperature jump there, so any layer whose balance point lands near
    X_c limit-cycles (T3B-N1/N2).  No rate law can smooth a binary
    state — at the PARTICLE scale the a_h transition really is a switch,
    spanning ~4 nm of shell (D_eff/k_gas).  What is wrong is not the
    switch but pretending every particle in a layer switches at once.

    THE CLOSURE, AND WHY IT IS PARAMETER-FREE.  PHY-032 (frozen,
    `occupancy_only`) prohibits "a fitted pore-blocking coefficient
    inferred from macro KPIs"; a tunable wetted fraction would violate
    that ruling, and none is needed.  The layer is already a CSTR, so
    particle age is exponentially distributed with the layer's own
    ``tau = M/F``; above X_c the release is HEAT-limited (T3B-N2
    ruling) and a heat-limited rate is CONSTANT in loading, so a
    particle entering at ``X_in`` reaches X_c at the age
    ``a_c = (X_in - X_c)/r`` and ``phi = 1 - exp(-a_c/tau)``.

    THE RATE THEN ELIMINATES ITSELF, and it must.  The CSTR steady-state
    mass balance is ``F(X_in - Xbar) = M r``, i.e. ``r tau = X_in -
    Xbar`` EXACTLY, so

        phi = 1 - exp( -(X_in - X_c) / (X_in - Xbar) )

    — three loadings the march already carries, no rate, no timescale,
    no parameter.

    A CLAIM I MADE HERE AND WITHDRAW (2026-08-08).  This docstring
    previously said the first version — taking ``r`` and ``tau`` as
    independent inputs — admitted a positive rate feedback that
    "MEASURABLY inverted the gas-capacity monotonicity", and that
    eliminating ``r`` through the mass balance removed it.  Both halves
    are false.  The two forms are ALGEBRAICALLY IDENTICAL at steady
    state, because ``r tau = X_in - Xbar`` IS the CSTR mass balance;
    measured, they produce BIT-IDENTICAL output.  So there was no
    feedback to remove, and the numbers the old text cited as the
    "before" case (starved 0.168 vs rich 0.152 kg/s, phi 0.385/0.780/
    0.951 against 0.000/0.530/0.823) are in fact the numbers this
    version still produces.  The gas-capacity inversion is OPEN, not
    fixed — it is recorded as T3B-N5.

    The mass-balance form is kept anyway, for the one reason that
    survives: away from steady state ``r tau`` and ``X_in - Xbar``
    genuinely differ, and only this form stays tied to the layer's own
    conserved inventory.  That is a transient-robustness argument, not
    the defect fix it was written up as.

    LIMITS (gated in tests/test_core2_sorption_interface.py):
      X_in <= X_c                  -> 0.0 exactly (nothing enters wetted)
      Xbar -> X_in (no depletion)  -> 1.0 (nothing dries)
      X_in >> X_c                  -> -> 1.0, film-active behaviour
      excess << depleted           -> excess/depleted, linear, no jump
      continuous in every argument; monotone increasing in X_in and
      decreasing in depletion.

    NOT CLAIMED: this is a mixing/population correction, not new
    particle physics.  It inherits the CSTR idealisation (a real tray
    is neither perfectly mixed nor plug flow) and monodispersity (a
    size distribution would broaden the transition further).  It also
    presumes the heat-limited rate binds above X_c — measured true at
    DT conditions; where transport binds there instead, the
    constant-rate premise weakens and phi should be taken from the
    actual trajectory rather than this closed form.
    """
    for name, value in (
        ("upstream loading", upstream_loading_kg_kg),
        ("mean loading", mean_loading_kg_kg),
        ("critical loading", critical_loading_kg_kg),
    ):
        if not math.isfinite(value):
            raise SorptionInterfaceError(f"{name} must be finite, got {value}")
    if critical_loading_kg_kg <= 0.0:
        raise SorptionInterfaceError("critical loading must be positive")

    excess = upstream_loading_kg_kg - critical_loading_kg_kg
    if excess <= 0.0:
        # Nothing enters the layer carrying free liquid.
        return 0.0
    depleted = upstream_loading_kg_kg - mean_loading_kg_kg
    if depleted <= 0.0:
        # The layer is not depleting (or is gaining): nothing dries, so
        # every particle that entered wetted still is.
        return 1.0
    return -math.expm1(-excess / depleted)
