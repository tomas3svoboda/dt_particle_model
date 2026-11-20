"""Film-based rate layer — the T2 fidelity keystone (grounding freeze
docs/T2_FILM_LAYER_GROUNDING_FREEZE_2026-08-06.md).

Replaces the equilibrium VLE rows with the finite-rate PHY-022 film,
closed per the frozen determinacy contract:

1. **Interface state (T_I, y_I)** from the PRODUCTION two-liquid
   authority stack (compressed-liquid Helmholtz fugacities + PHY-053
   φ), with the water equality scaled by the Luikov sorption activity
   CAPPED at 1 (F1 free water).  This continuously generalizes
   Tutkun's two-mode switch — a_w = 1 IS the certified heteroazeotrope
   pin T_E(P), y_E(P) bit-exactly; a_w < 1 lifts T_I above it.
2. **Total mass flux** from the interfacial energy row used exactly
   once (Krishna eq 42/43; Webb's closure): Ackermann-corrected
   gas-side heat = signed per-component latent load (PHY-053-exact
   capacities in β) + bed-side penetration sink.  Second bisection.
3. **Flux split** from the whitelisted rank-1 Stefan partition — the
   film never invents total flux (all four pinned sources concur).

The co/counter-flux regime (steam condensing while hexane evaporates)
is admitted with SIGNED per-component latent terms — the disclosed
model-form extrapolation of the grounding freeze §2.  Noncondensing
behavior emerges from the rates, never from a clamp.

DISCLOSED EXTENSION (A5-style): the Ackermann correction is evaluated
BEYOND the P1EF |β| ≤ 1 production pin (via the primitive's exposed
`enforce_production_domain=False`).  The condensation-driven flash
(FTRZ) regime lives at β ≫ 1 — superheated steam condensing onto the
wet bed pumps hexane evaporation through the latent row — and the
β/(e^β−1) form at arbitrary β is exactly the Krishna/Webb high-flux
film correction the grounding freeze pinned (Webb's q_s = e^ε·q_g
identity).  A numeric-sanity ceiling |β| ≤ 50 fails closed.  Fluxes are per
interfacial area; the caller applies the PHY-023/032 area exactly
once.  X_h below the PHY-046 threshold refuses (PHY-019 hand-off).
physically_qualifying is False.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import sorption_interface as _si
from .film import TrayLayerFilm, hexane_mass_fraction_from_mole_fraction
from .heteroazeotrope import (
    HeteroazeotropeError,
    HeteroazeotropeState,
    solve_heteroazeotrope,
)
from .props import hexane as _hexane
from .props import water as _water
from .tray import TrayParams, luikov_activity

# RETIRED (T3b).  Was 1.0e-9 and used as the PHY-046 film-active refusal
# threshold, i.e. 2.0e8 BELOW the actual critical loading X_c ~ 0.2017
# kg/kg.  Kept only so external references fail loudly rather than
# silently reading a stale constant; nothing in core2 uses it.
MOBILE_HEXANE_THRESHOLD_KG_KG = 1.0e-9
BETA_SANITY_CEILING = 50.0


class FilmBedLayerError(RuntimeError):
    """A film-layer contract was violated; the state is unusable."""


@dataclass(frozen=True)
class FilmBedLayerResult:
    """One rate-resolved interface solve (per unit interfacial area)."""

    interface_temperature_k: float
    interface_y_hexane: float
    water_activity: float
    heteroazeotrope_pinned: bool  # a_w capped at 1 (F1 free water)
    hexane_flux_kg_m2_s: float    # positive solid -> gas
    water_flux_kg_m2_s: float
    total_flux_kg_m2_s: float
    gas_side_heat_w_m2: float     # Ackermann-corrected, gas -> interface
    bed_side_heat_w_m2: float     # interface -> bed
    latent_load_w_m2: float       # signed per-component sum
    ackermann_factor: float
    energy_residual_w_m2: float
    # Film conductances at the solve state — for layer-scale NTU
    # (exchanger) integration by the caller.
    mass_conductance_kg_m2_s: float   # K = rho_film * k_m
    heat_conductance_w_m2_k: float    # h_Q
    # Rank-1 split components + interface/latent state for the
    # layer-scale energy row (F3 mode-split integration): the bulk
    # Stefan generation component scales with FULL area, the diffusive
    # component with the NTU-effective area.
    interface_w_hexane: float         # hexane MASS fraction at y_I
    diffusive_hexane_flux_kg_m2_s: float  # K·φ·(w_I − w_b); water = −this
    lambda_hexane_j_kg: float         # dh_vap/M at T_I
    lambda_water_j_kg: float          # IAPWS specific at T_I

    @property
    def physically_qualifying(self) -> bool:
        return False


def interface_state(
    *,
    water_activity: float,
    pressure_pa: float,
    hexane_activity: float = 1.0,
) -> HeteroazeotropeState:
    """Interface (T_I, y_I) from the PRODUCTION two-liquid authority
    stack (compressed-liquid fugacities + PHY-053 φ), generalized with
    BOTH sorption activities: y_w·φ_w·P = a_w·f_w^L and
    y_h·φ_h·P = a_h·f_h^L.  a_w = a_h = 1 is exactly the certified
    heteroazeotrope pin T_E(P), y_E(P)."""
    if not 0.0 < water_activity <= 1.0:
        raise FilmBedLayerError("water activity must lie in (0, 1]")
    if not 0.0 < hexane_activity <= 1.0:
        raise FilmBedLayerError("hexane activity must lie in (0, 1]")
    try:
        return solve_heteroazeotrope(
            pressure_pa,
            water_activity=water_activity,
            hexane_activity=hexane_activity,
        )
    except HeteroazeotropeError as exc:
        raise FilmBedLayerError(
            f"interface equilibrium unsolvable at a_w = {water_activity}, "
            f"a_h = {hexane_activity}, P = {pressure_pa} Pa: {exc}"
        ) from exc


def interface_temperature_k(
    *,
    water_activity: float,
    pressure_pa: float,
    hexane_activity: float = 1.0,
) -> float:
    """Interface temperature alone (see `interface_state`)."""
    return interface_state(
        water_activity=water_activity,
        pressure_pa=pressure_pa,
        hexane_activity=hexane_activity,
    ).temperature_k


@dataclass(frozen=True)
class InterfaceClosure:
    """The heat-side interface state alone — shared by the full film
    solve and the generation-only (G = 0) node of the sequence."""

    water_activity: float
    heteroazeotrope_pinned: bool
    temperature_k: float
    y_hexane: float
    w_hexane: float
    lambda_hexane_j_kg: float
    lambda_water_j_kg: float
    hexane_activity: float = 1.0
    hexane_regime: str = "FILM_ACTIVE"

    @property
    def physically_qualifying(self) -> bool:
        return False


def interface_closure(
    *,
    hexane_loading_kg_kg: float,
    retained_water_kg_kg: float,
    free_water_present: bool,
    pressure_pa: float,
    params: TrayParams | None = None,
    hexane_activity: float | None = None,
    hexane_regime: str = "FILM_ACTIVE",
    sorption_temperature_k: float | None = None,
    particle: _si.ParticleSorptionParams | None = None,
) -> InterfaceClosure:
    """Interface state + mass-basis latents with the film refusals.

    T3b: `hexane_activity` is the SURFACE n-hexane activity of the
    representative particle.  Callers that march a layer state supply it
    (with its regime label) from `sorption_interface.surface_state`,
    which carries the intraparticle gradient.  Callers that only have a
    loading may pass `sorption_temperature_k` instead and get the
    zero-gradient (thermodynamic-only) activity computed here.

    The old `MOBILE_HEXANE_THRESHOLD_KG_KG = 1e-9` refusal is GONE — it
    was placed 2.0e8 below the PHY-046 critical loading X_c ~ 0.2017, so
    every state between the two ran under the film-active a_h = 1 law
    without a flag, which is how six S2 campaign cells were reported as
    passes while discharging below the frozen floor.  What replaces it is
    a regime DISPATCH, not an earlier refusal: see `sorption_interface`
    for why moving the threshold to X_c (the audit's literal
    recommendation) would have been wrong.
    """
    p = params or TrayParams()
    if not math.isfinite(hexane_loading_kg_kg) or hexane_loading_kg_kg <= 0.0:
        raise FilmBedLayerError(
            f"n-hexane loading {hexane_loading_kg_kg} is not positive — there "
            "is no n-hexane interface to close (fail-closed)"
        )
    if hexane_activity is None:
        # Zero-gradient thermodynamic activity: the sorption half without
        # the transport half.  Honest for a standalone equilibrium query,
        # OPTIMISTIC for a marched layer (it ignores the intraparticle
        # gradient, which at the residual plateau is worth another ~40x).
        reference_t = (
            sorption_temperature_k
            if sorption_temperature_k is not None
            else None
        )
        if reference_t is None:
            raise FilmBedLayerError(
                "interface_closure needs either an explicit hexane_activity "
                "or a sorption_temperature_k at which to evaluate it"
            )
        surface = _si.surface_state(
            mean_loading_kg_kg=hexane_loading_kg_kg,
            temperature_k=reference_t,
            hexane_flux_kg_m2_s=0.0,
            params=particle,
            pressure_pa=pressure_pa,
        )
        hexane_activity = surface.surface_activity
        hexane_regime = surface.regime.name
    activity = 1.0 if free_water_present else min(
        luikov_activity(retained_water_kg_kg, p), 1.0
    )
    interface = interface_state(
        water_activity=activity,
        pressure_pa=pressure_pa,
        hexane_activity=hexane_activity,
    )
    T_I = interface.temperature_k
    # Latent heats on a MASS basis: the hexane authority reports dh_vap
    # in J/mol (molar Helmholtz form), the water authority in J/kg
    # (IAPWS-95 specific form) — divide only hexane by M.
    return InterfaceClosure(
        water_activity=activity,
        heteroazeotrope_pinned=activity >= 1.0 and hexane_activity >= 1.0,
        temperature_k=T_I,
        y_hexane=interface.y_hexane,
        w_hexane=hexane_mass_fraction_from_mole_fraction(interface.y_hexane),
        lambda_hexane_j_kg=_hexane.saturation(T_I).dh_vap / _hexane.M,
        lambda_water_j_kg=_water.saturation(T_I).dh_vap,
        hexane_activity=hexane_activity,
        hexane_regime=hexane_regime,
    )


def solve_film_bed_layer(
    *,
    film: TrayLayerFilm,
    bed_temperature_k: float,
    hexane_loading_kg_kg: float,
    retained_water_kg_kg: float,
    free_water_present: bool,
    gas_temperature_k: float,
    gas_y_hexane: float,
    pressure_pa: float,
    bed_contact_h_w_m2_k: float,
    params: TrayParams | None = None,
    max_flux_kg_m2_s: float = 1.0,
    hexane_activity: float | None = None,
    hexane_regime: str = "FILM_ACTIVE",
    sorption_temperature_k: float | None = None,
    particle: _si.ParticleSorptionParams | None = None,
) -> FilmBedLayerResult:
    """Solve one rate-resolved interface per the frozen contract."""
    if bed_contact_h_w_m2_k <= 0.0:
        raise FilmBedLayerError("bed contact conductance must be positive")
    closure = interface_closure(
        hexane_loading_kg_kg=hexane_loading_kg_kg,
        retained_water_kg_kg=retained_water_kg_kg,
        free_water_present=free_water_present,
        pressure_pa=pressure_pa,
        params=params,
        hexane_activity=hexane_activity,
        hexane_regime=hexane_regime,
        sorption_temperature_k=(
            sorption_temperature_k
            if sorption_temperature_k is not None
            else bed_temperature_k
        ),
        particle=particle,
    )
    activity = closure.water_activity
    pinned = closure.heteroazeotrope_pinned
    T_I = closure.temperature_k
    w_interface = closure.w_hexane
    w_bulk = hexane_mass_fraction_from_mole_fraction(gas_y_hexane)

    coefficients = film.coefficients(
        bulk_temperature_k=gas_temperature_k, interface_temperature_k=T_I
    )
    K = coefficients.film_density_kg_m3 * coefficients.binary_mass_transfer_m_s
    lambda_h = closure.lambda_hexane_j_kg
    lambda_w = closure.lambda_water_j_kg

    # Flux-independent quantities hoisted out of the bisection: the
    # film coefficients and the PHY-053 secant capacities depend only
    # on (T_g, T_I, P, y) — rebuilt per solve, never per iterate.
    from .particle import external_film_reduced as efr

    # Secant degeneracy guard: below the authority's own derivative
    # stencil scale (ldexp(T, −12) ≈ 0.08 K) the raw secant cancels
    # catastrophically (an NTU→∞ cascade delivers gas within ULPs of
    # T_I), so snap to the certified equal-temperature limit branch —
    # a ≤0.1 J/kgK perturbation against ~2000.  Heats keep the true
    # ΔT (linear, no cancellation).
    bulk_for_secant = gas_temperature_k
    if abs(gas_temperature_k - T_I) < math.ldexp(T_I, -12):
        bulk_for_secant = T_I
    capacities = efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=T_I,
        bulk_temperature_k=bulk_for_secant,
        pressure_pa=pressure_pa,
        y_hexane=gas_y_hexane,
    )

    def ackermann(j_h: float, j_w: float):
        # DISCLOSED EXTENSION: production-domain pin lifted — the FTRZ
        # flash regime needs beta >> 1 (module docstring); the numeric
        # sanity ceiling replaces the P1EF envelope pin.
        heat = efr.ackermann_heat_transfer(
            coefficients,
            fluxes=efr.ComponentMassFluxes(
                water_kg_m2_s=j_w, hexane_kg_m2_s=j_h
            ),
            water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
            enforce_production_domain=False,
        )
        if abs(heat.beta) > BETA_SANITY_CEILING:
            raise FilmBedLayerError(
                f"Ackermann beta {heat.beta:.1f} beyond the numeric "
                f"sanity ceiling {BETA_SANITY_CEILING}"
            )
        return heat

    flux_bound = max_flux_kg_m2_s

    def split(j_t: float) -> tuple[float, float]:
        """Rank-1 conditioned Stefan partition (φ(Pe) form)."""
        pe = j_t / K
        if abs(pe) < 1.0e-12:
            phi = 1.0
        else:
            phi = pe / math.expm1(pe)
        j_h = w_interface * j_t + K * phi * (w_interface - w_bulk)
        return j_h, j_t - j_h

    def energy_residual(j_t: float) -> float:
        j_h, j_w = split(j_t)
        heat = ackermann(j_h, j_w)
        q_gas = heat.heat_into_particle_w_m2
        q_bed = bed_contact_h_w_m2_k * (T_I - bed_temperature_k)
        latent = j_h * lambda_h + j_w * lambda_w
        return q_gas - q_bed - latent

    low, high = -flux_bound, flux_bound
    e_low, e_high = energy_residual(low), energy_residual(high)
    if e_low * e_high > 0.0:
        raise FilmBedLayerError(
            "interfacial energy row has no root inside the flux bracket "
            f"[{low:.4f}, {high:.4f}] kg/m2/s "
            f"(E(low) = {e_low:.3e}, E(high) = {e_high:.3e}) — raise "
            "max_flux_kg_m2_s only with a defended envelope"
        )
    for _ in range(200):
        mid = 0.5 * (low + high)
        e_mid = energy_residual(mid)
        if e_low * e_mid <= 0.0:
            high, e_high = mid, e_mid
        else:
            low, e_low = mid, e_mid
        if high - low < 1.0e-15:
            break
    j_t = 0.5 * (low + high)
    j_h, j_w = split(j_t)
    heat = ackermann(j_h, j_w)
    q_gas = heat.heat_into_particle_w_m2
    q_bed = bed_contact_h_w_m2_k * (T_I - bed_temperature_k)
    latent = j_h * lambda_h + j_w * lambda_w

    return FilmBedLayerResult(
        interface_temperature_k=T_I,
        interface_y_hexane=closure.y_hexane,
        water_activity=activity,
        heteroazeotrope_pinned=pinned,
        hexane_flux_kg_m2_s=j_h,
        water_flux_kg_m2_s=j_w,
        total_flux_kg_m2_s=j_t,
        gas_side_heat_w_m2=q_gas,
        bed_side_heat_w_m2=q_bed,
        latent_load_w_m2=latent,
        ackermann_factor=heat.factor,
        energy_residual_w_m2=q_gas - q_bed - latent,
        mass_conductance_kg_m2_s=K,
        heat_conductance_w_m2_k=coefficients.heat_transfer_w_m2_k,
        interface_w_hexane=w_interface,
        diffusive_hexane_flux_kg_m2_s=j_h - w_interface * j_t,
        lambda_hexane_j_kg=lambda_h,
        lambda_water_j_kg=lambda_w,
    )
