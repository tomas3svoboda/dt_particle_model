"""Film-consistent depletion sequence â€” the (b) refinement's rate-true
instrument (grounding freeze docs/T2_FILM_LAYER_GROUNDING_FREEZE_2026-08-06.md Â§5).

Replaces the equilibrium layer solve of `depletion_sequence` with the
finite-rate film closure of `film_bed_layer`: per macro-step and layer,
the interface state (T_I, y_I) and the component fluxes come from the
two-mode film solve at the layer's CURRENT inventories and gas state,
scaled by the PHY-023/032 interfacial area supplied per layer.  The
certified ensemble rates and their supply caps disappear â€” the film IS
the rate law, and gas carrying capacity emerges from it instead of
being imposed as a cap.

State per layer (bottom-up, L1 = index 0): sorbed hexane X_h and water
X_w (kg/kg dry), free-water inventory (kg), bed temperature T_bed (K)
â€” the bed temperature is now a marched state closed by the film's
bed-side heat, no longer an output of an equilibrium tray solve.

Tranche accounting rules (documented, mirroring the sibling modules):

- Gas node per layer is INLET-CLOSED with the F3 MODE-SPLIT layer
  energy row (findings record F3, adversarially verified): the
  DIFFUSIVE component of the rank-1 split uses the NTU-effective area
  AÂ·(1â’e^{â’NTU})/NTU (its driving force depletes along the gas path);
  the bulk GENERATION component is powered by heat at the T_I-pinned
  interface â€” full-area bed contact + NTU-exact gas sensible â€” and is
  the boiling / evaporation-generated-flow channel.  Generated vapor
  joins the stream at T_I (mass-weighted merge).  NTUâ†’0 recovers
  fluxĂ—A; NTUâ†’âž gives contact equilibrium PLUS wall-powered boiling.
  The layer energy row closes by construction (residual reported;
  nonzero only under asymmetric truncation).  The enthalpy-exact node
  arrives with the dynamic wall node.
- Bed temperature and inventories update by EXACT LINEARIZED
  RELAXATION over the step (sources frozen over dt) â€” the CSTR
  discipline of the ruled (b) design map.
- Water condensed beyond the Luikov cap spills to the free-water
  inventory (F1 mechanism); free water present â‡’ the interface is
  heteroazeotrope-pinned.
- Exhaustion truncation is exact bookkeeping: a layer cannot release
  more of a species than inventory/dt + advective inflow, and cannot
  condense more than the gas inflow carries.  Film heats scale with
  the same factor (proportional rule; exact in the untruncated
  converged regime).
- Fail-closed: any FilmBedLayerError (PHY-019 dry bed, unsolvable
  interface, energy row outside the Ackermann production bracket)
  aborts the sequence with step/layer context â€” the toasting boundary
  stays observable, never faked.
- W-T2-2 LIFTED (this instrument): a zero gas stream â€” predesolv
  inlet, or mid-cascade after a cold bed condenses the entire stream
  (steam-up) â€” runs as a GENERATION-ONLY node: the layer boils its
  own vapor at (T_I, y_I) powered by the bed contact; a bed below its
  interface state with no vapor present is inert.  The
  evaporation-generated-flow treatment of findings record F3.

physically_qualifying is False throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import sorption_interface as _si
from .film import TrayLayerFilm
from .film_bed_layer import (
    FilmBedLayerError,
    interface_closure,
    solve_film_bed_layer,
)
from .tray import GasStream, SolidStream, TrayParams

DISCHARGE_STABILITY_REL = 1.0e-6
STABILITY_LOOKBACK_STEPS = 10
_R_GAS = 8.31451
_HEXANE_M = 0.08618      # kg/mol, frozen Span-Wagner molar mass
# T3b night tranche: below this loading a layer is EXHAUSTED — no
# n-hexane interface exists; the interface is the water-only a_h -> 0
# limit of the same closure and the hexane release is identically zero.
EXHAUSTED_LOADING_KG_KG = 1.0e-12
# The a_h value used to evaluate the water-only limit.  Measured: the
# generalized solve converges smoothly for a_h down to at least 1e-6
# with T_I within 0.01 K of the a_h -> 0 sorbed-water boiling point;
# 1e-15 keeps the hexane partial pressure physically negligible without
# leaving the solver's (0, 1] domain.
WATER_ONLY_HEXANE_ACTIVITY = 1.0e-15
# NOTE (T3b, why there is no relaxation parameter here): the surface
# activity was first estimated from the PREVIOUS step's n-hexane flux
# under an under-relaxation.  That closes a feedback loop of gain
# (c_mean - c_s)/c_s, order 50 at the residual plateau, and the march
# measurably failed to converge at every relaxation tried (0.5 and 0.2
# gave bit-identical limit cycles).  The series balance is now solved in
# closed form by `sorption_interface.surface_state_from_gas`, so there
# is nothing to relax.


class FilmDepletionError(RuntimeError):
    """A film-sequence contract was violated; the trajectory is unusable."""


def _phi_rate_latent_candidate(
    *,
    wetted_fraction: float,
    film_hexane_rate_kg_s: float,
    series_hexane_rate_kg_s: float,
    interface_heat_w: float,
    lambda_hexane_j_kg: float,
    lambda_water_j_kg: float,
) -> tuple[float, float, float]:
    r"""Candidate reduced-host blend of wet-film and dry-series populations.

    CANDIDATE, NOT AN OWNER-APPROVED PHYSICS CLOSURE.  The latent-duty
    identity below is exact, but it is not the complete common-datum first
    law (which also carries conductive and component-enthalpy terms).  The
    reduced host carries only one
    averaged interface state and therefore does not uniquely determine
    whether unused film-hexane heat boils water, heats the bed, or is
    partitioned between separately resolved wet/dry interfaces.

    ``interface_heat_w`` is the heat already admitted by the F3 gas/bed
    interface row.  Only the wetted population's film-hexane release is
    charged to that row; the dry population's series release draws its
    latent heat directly from the bed.  Water therefore receives the
    remainder of the *same* interface heat::

        r_h = phi*r_film + (1-phi)*r_series
        L_series = (1-phi)*r_series*lambda_h
        r_w = (Q_interface - phi*r_film*lambda_h)/lambda_w

    Hence, before inventory truncation,

        r_h*lambda_h + r_w*lambda_w
            = Q_interface + L_series

    for every phi in [0, 1].  The phi=1 endpoint recovers the existing
    pure-film row.  At phi=0 the hexane rate recovers the pure-series limit,
    but assigning the remaining interface duty to water is a new closure,
    not a consequence of conservation alone.  Values outside [0, 1] are
    refused, not clamped: phi is an upstream population result, not a tuning
    knob.  This helper belongs to the reduced-host/oracle model and is not a
    substitute for the resolved particle energy equation.  Its heat
    allocation and march-level qualification remain open.  It is retained as
    an isolated algebra experiment and is intentionally not called by
    ``FilmDepletionSequence``.
    """
    values = {
        "wetted fraction": wetted_fraction,
        "film hexane rate": film_hexane_rate_kg_s,
        "series hexane rate": series_hexane_rate_kg_s,
        "interface heat": interface_heat_w,
        "hexane latent heat": lambda_hexane_j_kg,
        "water latent heat": lambda_water_j_kg,
    }
    for name, value in values.items():
        if not math.isfinite(value):
            raise FilmDepletionError(f"{name} must be finite, got {value}")
    if not 0.0 <= wetted_fraction <= 1.0:
        raise FilmDepletionError(
            "wetted fraction must lie in [0, 1]; "
            f"got {wetted_fraction} (refuse, do not clamp)"
        )
    if lambda_hexane_j_kg <= 0.0 or lambda_water_j_kg <= 0.0:
        raise FilmDepletionError("latent heats must be positive")

    dry_fraction = 1.0 - wetted_fraction
    hexane_rate = (
        wetted_fraction * film_hexane_rate_kg_s
        + dry_fraction * series_hexane_rate_kg_s
    )
    series_latent_w = (
        dry_fraction * series_hexane_rate_kg_s * lambda_hexane_j_kg
    )
    water_rate = (
        interface_heat_w
        - wetted_fraction * film_hexane_rate_kg_s * lambda_hexane_j_kg
    ) / lambda_water_j_kg
    return hexane_rate, water_rate, series_latent_w


@dataclass(frozen=True)
class FilmLayerStepRecord:
    """Per-layer film diagnostics at one macro-step."""

    interface_temperature_k: float
    heteroazeotrope_pinned: bool
    hexane_rate_kg_s: float   # positive bed -> gas, post-truncation
    water_rate_kg_s: float
    truncation_factor: float  # 1 when the film rate fit the inventory
    gas_side_heat_w: float
    bed_side_heat_w: float    # interface -> bed (negative = bed powers boil)
    energy_residual_w: float = 0.0  # latent - Q_gas - Q_bed->I after truncation
    # T3b sub-critical n-hexane diagnostics
    hexane_regime: str = "FILM_ACTIVE"
    hexane_surface_activity: float = 1.0
    front_radius_fraction: float = 1.0       # s/R
    supply_limited_flux_kg_m2_s: float = math.inf
    supply_limited: bool = False             # release hit the particle's cap
    wetted_fraction: float = 0.0             # T3b (a) phi: wetted population


@dataclass(frozen=True)
class FilmDepletionStepRecord:
    step: int
    layer_loadings_x_h: tuple[float, ...]
    layer_loadings_x_w: tuple[float, ...]
    layer_free_water_kg: tuple[float, ...]
    layer_bed_temperatures_k: tuple[float, ...]
    hexane_release_kg_s: float
    water_release_kg_s: float
    films: tuple[FilmLayerStepRecord, ...]


@dataclass(frozen=True)
class FilmDepletionResult:
    steps: tuple[FilmDepletionStepRecord, ...]
    converged: bool
    residence_avg_hexane_release_kg_s: float
    residence_avg_water_release_kg_s: float
    discharge_x_h: float
    discharge_temperature_k: float
    ledger_residual_rel: float
    # Max per-layer |latent â’ Q_gas â’ Q_bedâ†’I| over the final
    # residence window (W): the layer energy row closes by
    # construction; nonzero only where truncation binds asymmetrically.
    max_energy_residual_w: float = 0.0
    # True when the run deferred the N2 heat cap (screening tier).
    # Screening results are for architecture/geometry screening ONLY;
    # their energetics carry the disclosed fictitious-latent imbalance.
    energy_screening_tier: bool = False

    @property
    def physically_qualifying(self) -> bool:
        return False


class FilmDepletionSequence:
    """Time-marched K-layer depletion closed by the film rate law."""

    physically_qualifying = False

    def __init__(
        self,
        *,
        films: tuple[TrayLayerFilm, ...],
        interfacial_area_m2: tuple[float, ...],
        bed_contact_h_w_m2_k: tuple[float, ...],
        gas_heat_capacity_j_kg_k: float,
        # T3a: the GAS-side exchange area may differ from the bed-side
        # one.  On an unperforated predesolventizer disc the vapour can
        # only leave through the free surface, but the bed still boils
        # through its whole volume against its jacket — so restricting
        # both paths would starve the bed of its own duty.  Defaults to
        # the interfacial area, i.e. the previous behaviour exactly.
        gas_side_area_m2: tuple[float, ...] | None = None,
        layer_holdup_kg: float | tuple[float, ...],
        layer_params: tuple[TrayParams, ...] | None = None,
        wall_duty_w: float = 0.0,
        layer_duty_w: tuple[float, ...] | None = None,
        macro_steps_per_residence: int = 40,
        max_residences: float = 10.0,
        # T3b: the certified particle constants the sub-critical n-hexane
        # closure runs on.  Defaults are the frozen soybean values that
        # `particle/dry_shell` and `props/critical_volume` already carry.
        particle: _si.ParticleSorptionParams | None = None,
        # T3B SCREENING TIER (owner-accepted 2026-08-08, ruling doc
        # bottom line, option (c)): defer the N2 heat cap so the
        # instrument converges on the softer transport-only dynamics.
        # The energy imbalance this admits is NOT hidden — it lands in
        # the machine energy ledger exactly as before (measured floor
        # 0.201 relative) and every result carries this flag, so no
        # artifact can pass itself off as honest-energetics.  Default
        # OFF: the ruled physics is the default.
        energy_screening_tier: bool = False,
    ) -> None:
        K = len(films)
        if K < 1:
            raise FilmDepletionError("at least one layer is required")
        if not (K == len(interfacial_area_m2) == len(bed_contact_h_w_m2_k)):
            raise FilmDepletionError("layer configuration lengths must agree")
        if any(a <= 0.0 for a in interfacial_area_m2):
            raise FilmDepletionError("interfacial areas must be positive")
        # T3a generalization: per-cell holdup (trays carry different bed
        # depths) and per-cell wall duty (each tray has its own jacket).
        # Both default to the previous scalar behaviour EXACTLY, so the
        # single-tray gates stay bit-identical.
        holdups = (
            tuple(layer_holdup_kg)
            if isinstance(layer_holdup_kg, (tuple, list))
            else (float(layer_holdup_kg),) * K
        )
        if len(holdups) != K:
            raise FilmDepletionError("layer holdup lengths must agree")
        if any(h <= 0.0 for h in holdups):
            raise FilmDepletionError("layer holdup must be positive")
        if layer_duty_w is not None and len(layer_duty_w) != K:
            raise FilmDepletionError("layer duty lengths must agree")
        if macro_steps_per_residence < 5:
            raise FilmDepletionError("macro-step grid too coarse")
        if gas_heat_capacity_j_kg_k <= 0.0:
            raise FilmDepletionError("gas heat capacity must be positive")
        if gas_side_area_m2 is not None:
            if len(gas_side_area_m2) != K:
                raise FilmDepletionError("gas-side area lengths must agree")
            if any(a <= 0.0 for a in gas_side_area_m2):
                raise FilmDepletionError("gas-side areas must be positive")
        self._films = films
        self._areas = interfacial_area_m2
        self._gas_areas = tuple(gas_side_area_m2 or interfacial_area_m2)
        self._contact_h = bed_contact_h_w_m2_k
        self._cp_gas = gas_heat_capacity_j_kg_k
        self._K = K
        self._holdups = holdups
        self._params = layer_params or tuple(TrayParams() for _ in range(K))
        self._wall_duty_w = wall_duty_w
        self._layer_duty_w = (
            tuple(layer_duty_w)
            if layer_duty_w is not None
            else tuple(wall_duty_w if k == 0 else 0.0 for k in range(K))
        )
        self._steps_per_residence = macro_steps_per_residence
        self._max_residences = max_residences
        self._particle = particle or _si.ParticleSorptionParams()
        self._screening_tier = bool(energy_screening_tier)

    def run(
        self,
        *,
        feed: SolidStream,
        gas_inlet: GasStream,
    ) -> FilmDepletionResult:
        K = self._K
        F = feed.F_dm
        if F <= 0.0:
            raise FilmDepletionError("feed dry-mass rate must be positive")
        # Per-cell residence; the macro step is set by the TOTAL
        # residence so a machine with unequal trays still resolves its
        # own turnover.  Reduces exactly to tau_layer*K/steps when every
        # holdup is equal (the single-tray case).
        tau = [h / F for h in self._holdups]
        dt = sum(tau) / self._steps_per_residence
        max_steps = int(self._steps_per_residence * self._max_residences)

        x_h = [feed.X_h] * K
        x_w = [feed.X_w] * K
        t_bed = [feed.T] * K
        free_water = [0.0] * K
        cap = self._params[0].w_water_cap
        # T3b lagged states for the sub-critical n-hexane closure.  The
        # interface temperature is where the GAB isotherm is evaluated
        # (it is the surface's own temperature, and at cyclic steady
        # state the film pins the bed to it); the flux carries the
        # intraparticle gradient.  Both are lagged one macro-step,
        # exactly the discipline PASS 1 already uses for layer states.
        # T3b: seed the lagged INTERFACE temperature from a real
        # interface solve at the feed state.  Seeding it with feed.T
        # would hand the reduced-film envelope a temperature below its
        # 323.15 K floor on the very first step -- the reference machine
        # feeds at 322.15 K, one kelvin under it.  An interface
        # temperature is a boiling point and is always inside the band.
        # T3b track (a): per-layer specific release (kg hexane per kg dry
        # meal per second), lagged one macro-step — the wetted-fraction
        # closure needs the drying rate a particle actually experiences.
        j_h_specific = [0.0] * K
        t_interface = []
        for k in range(K):
            try:
                seed = interface_closure(
                    hexane_loading_kg_kg=feed.X_h,
                    retained_water_kg_kg=feed.X_w,
                    free_water_present=False,
                    pressure_pa=self._params[k].pressure_pa,
                    params=self._params[k],
                    sorption_temperature_k=feed.T,
                    particle=self._particle,
                )
            except (FilmBedLayerError, _si.SorptionInterfaceError) as exc:
                raise FilmDepletionError(
                    f"cannot seed the interface temperature of layer {k} "
                    f"from the feed state (X_h = {feed.X_h}, X_w = "
                    f"{feed.X_w}, T = {feed.T} K): {exc}"
                ) from exc
            t_interface.append(seed.temperature_k)

        records: list[FilmDepletionStepRecord] = []
        history: list[float] = []
        history_w: list[float] = []
        total_released_h = 0.0
        fed_h = 0.0
        discharged_h = 0.0
        initial_inventory_h = sum(
            x * h for x, h in zip(x_h, self._holdups)
        )
        converged = False

        for step in range(1, max_steps + 1):
            gas_in = gas_inlet
            step_release_h = 0.0
            step_release_w = 0.0
            film_records: list[FilmLayerStepRecord] = []
            rates: list[tuple[float, float]] = []
            supply_drops: list[float] = []

            # PASS 1 (bottom-up): film solves along the gas cascade at
            # the current (lagged) layer states.
            for k in range(K):
                gas_total = gas_in.G_w + gas_in.G_h
                # W-T2-2 LIFTED (film instrument): a zero gas stream is
                # a GENERATION-ONLY node â€” the layer boils its own
                # vapor at (T_I, y_I), powered by the bed contact.
                # Reached by predesolv inlets AND mid-cascade when a
                # cold bed condenses the entire stream (steam-up).
                generation_only = gas_total <= 0.0
                y_h = math.nan
                film = None
                try:
                    # --- T3b: the sub-critical n-hexane surface state ---
                    # Which constitutive law owns the surface, and what
                    # the intraparticle gradient does to its activity.
                    # The estimator's flux is bounded by the particle's
                    # own supply so the closure stays inside its domain;
                    # the PHYSICAL cap is applied to the release below.
                    pressure_k = self._params[k].pressure_pa
                    # SERIES resistance solved in closed form.  Two
                    # things must NOT be lagged, for different reasons.
                    #
                    # (1) the FLUX: lagging it closes a loop of gain
                    #     (c_mean - c_s)/c_s ~ 50 at the residual
                    #     plateau, which no under-relaxation stabilises.
                    #     Solved in closed form instead.
                    # (2) the temperature at which the regime boundary
                    #     X_c(T) is evaluated must not be the INTERFACE
                    #     temperature, because T_I is itself a function
                    #     of the regime: a wetted surface sits at the
                    #     ~337 K azeotrope and a dry one at the ~365 K
                    #     sorption-elevated point, and X_c moves 0.2007
                    #     -> 0.1917 between them.  A layer inside that
                    #     window then flips regime every step and the
                    #     march runs a clean limit cycle (measured:
                    #     period 15, X_h swinging 0.157 <-> 0.193).
                    #     The BED temperature is a marched thermal state
                    #     with real inertia and no dependence on the
                    #     activity, so it breaks the circularity — and
                    #     it is the physically right argument anyway,
                    #     since X_c is the particle's pore-saturation
                    #     property at the particle's own temperature.
                    #
                    # TWO DISTINCT REFERENCE TEMPERATURES, and collapsing
                    # them into one was a defect.  The SORPTION side must
                    # use the bed temperature, per (2).  The FILM
                    # coefficients must use the INTERFACE temperature, as
                    # they always did: the frozen reduced-film envelope
                    # requires both its temperatures in [323.15, 433] K,
                    # and an interface temperature is a boiling point so
                    # it never leaves that band — whereas a bed can, and
                    # the reference machine feeds at 322.15 K, one kelvin
                    # under the floor.  Driving the film off t_bed
                    # refused the whole machine at layer 0, step 1.
                    t_sorb = t_bed[k]
                    t_film = t_interface[k]
                    # lagged input for the wetted-fraction mixture: the
                    # upstream loading a particle ENTERS with.  The
                    # DEPLETION (X_in - Xbar) supplies the rate through
                    # the CSTR mass balance, so no rate is passed.
                    up_x_h_lag = x_h[k + 1] if k + 1 < K else feed.X_h
                    coefficients = self._films[k].coefficients(
                        bulk_temperature_k=gas_in.T,
                        interface_temperature_k=t_film,
                    )
                    k_gas = coefficients.binary_mass_transfer_m_s
                    y_bulk_h = (
                        _mole_fraction_hexane(gas_in)
                        if gas_total > 0.0
                        else 0.0
                    )
                    # Bulk n-hexane density on the SORPTION temperature
                    # basis, so that c_gas/c_sat is a true activity ratio
                    # against the c_sat the closure computes internally.
                    c_gas_h = (
                        y_bulk_h * pressure_k * _HEXANE_M
                        / (_R_GAS * t_sorb)
                    )
                    if x_h[k] <= EXHAUSTED_LOADING_KG_KG:
                        # EXHAUSTED: no n-hexane interface exists.  The
                        # interface is the water-only a_h -> 0 limit of
                        # the same closure (continuous — measured T_I
                        # within 0.01 K of the limit at a_h = 1e-6) and
                        # the hexane release is identically zero.
                        # Re-adsorption of gas-phase hexane onto
                        # exhausted meal is NOT represented; refuse if
                        # the gas could drive it rather than ignore it.
                        surface = None
                        a_h_surface = WATER_ONLY_HEXANE_ACTIVITY
                        surface_regime = "EXHAUSTED"
                        j_supply = 0.0
                        wetted_phi = 0.0
                        if c_gas_h > 0.5 * _si._dry_shell.cg_saturation(
                            t_sorb
                        ):
                            raise FilmDepletionError(
                                f"layer {k} is hexane-exhausted under a "
                                f"hexane-rich gas (c_gas = {c_gas_h:.3e} "
                                "kg/m3) — re-adsorption is not "
                                "represented; fail closed"
                            )
                    else:
                        surface = _si.surface_state_from_gas(
                            mean_loading_kg_kg=x_h[k],
                            temperature_k=t_sorb,
                            gas_hexane_density_kg_m3=c_gas_h,
                            film_conductance_m_s=k_gas,
                            params=self._particle,
                            pressure_pa=pressure_k,
                        )
                        a_h_surface = surface.surface_activity
                        surface_regime = surface.regime.name
                        j_supply = surface.supply_limited_flux_kg_m2_s
                        # --- T3b track (a): WETTED-FRACTION MIXTURE ---
                        # A lumped layer cannot be partly wetted, so its
                        # surface flips binarily at X_c and takes a ~48 K
                        # interface jump with it — the T3B-N1/N2 limit
                        # cycles.  A real layer holds a POPULATION: with
                        # the CSTR age distribution the module already
                        # assumes and the heat-limited (hence constant)
                        # depletion above X_c, the wetted fraction is
                        # closed-form and parameter-free (PHY-032
                        # forbids a fitted blocking coefficient; none is
                        # used).  The layer's activity is the phi-
                        # weighted mixture of the wetted surface (a_h=1)
                        # and the dry one.  phi = 0 and phi = 1 recover
                        # today's branches exactly.
                        phi = _si.wetted_fraction(
                            upstream_loading_kg_kg=up_x_h_lag,
                            mean_loading_kg_kg=x_h[k],
                            critical_loading_kg_kg=(
                                surface.critical_loading_kg_kg
                            ),
                        )
                        if phi > 0.0:
                            a_h_surface = phi * 1.0 + (1.0 - phi) * (
                                a_h_surface
                            )
                            if phi >= 1.0:
                                surface_regime = "FILM_ACTIVE"
                            elif surface_regime != "FILM_ACTIVE":
                                surface_regime = "MIXED_WETTED"
                        wetted_phi = phi
                    if generation_only:
                        closure = interface_closure(
                            hexane_loading_kg_kg=x_h[k],
                            retained_water_kg_kg=x_w[k],
                            free_water_present=free_water[k] > 0.0,
                            pressure_pa=pressure_k,
                            params=self._params[k],
                            hexane_activity=a_h_surface,
                            hexane_regime=surface_regime,
                        )
                        T_I = closure.temperature_k
                        pinned = closure.heteroazeotrope_pinned
                        w_i = closure.w_hexane
                        lam_h = closure.lambda_hexane_j_kg
                        lam_w = closure.lambda_water_j_kg
                    else:
                        y_h = _mole_fraction_hexane(gas_in)
                        film = solve_film_bed_layer(
                            film=self._films[k],
                            bed_temperature_k=t_bed[k],
                            hexane_loading_kg_kg=x_h[k],
                            retained_water_kg_kg=x_w[k],
                            free_water_present=free_water[k] > 0.0,
                            gas_temperature_k=gas_in.T,
                            gas_y_hexane=y_h,
                            pressure_pa=pressure_k,
                            bed_contact_h_w_m2_k=self._contact_h[k],
                            params=self._params[k],
                            hexane_activity=a_h_surface,
                            hexane_regime=surface_regime,
                        )
                        T_I = film.interface_temperature_k
                        pinned = film.heteroazeotrope_pinned
                        w_i = film.interface_w_hexane
                        lam_h = film.lambda_hexane_j_kg
                        lam_w = film.lambda_water_j_kg
                    t_interface[k] = T_I
                except _si.SorptionInterfaceError as exc:
                    raise FilmDepletionError(
                        f"sub-critical n-hexane closure refused at layer {k}, "
                        f"step {step} (X_h = {x_h[k]:.6e}, T_I = "
                        f"{t_interface[k]:.2f} K): {exc}"
                    ) from exc
                except FilmBedLayerError as exc:
                    raise FilmDepletionError(
                        f"film solve refused at layer {k}, step {step} "
                        f"(X_h = {x_h[k]:.3e}, T_bed = {t_bed[k]:.2f} K): "
                        f"{exc}"
                    ) from exc
                except Exception as exc:  # whitelisted fail-closed authorities
                    if type(exc).__name__ not in (
                        "ReducedFilmValidityError",
                        "ValueError",
                    ):
                        raise
                    raise FilmDepletionError(
                        f"frozen authority refused at layer {k}, step "
                        f"{step} (T_gas = {gas_in.T:.6f} K, y_h = "
                        f"{y_h:.6f}): {type(exc).__name__}: {exc}"
                    ) from exc

                # F3 MODE-SPLIT layer integration (verified fix,
                # findings record F3): the layer-scale interfacial
                # energy row is solved ONCE at layer scale â€”
                #   Q_gas + Q_bedâ†’I = R_tÂ·Î»_mix + DÂ·(Î»_h â’ Î»_w)
                # with the rank-1 split R_h = w_IÂ·R_t + D.  The
                # DIFFUSIVE component D depletes its driving force
                # along the gas path (NTU-effective area); the bulk
                # GENERATION component R_t is powered by heat at the
                # T_I-pinned interface and scales with FULL area
                # (bed side) / the NTU-exact gas sensible (gas side).
                # NTUâ†’0 recovers fluxĂ—A; NTUâ†’âž gives contact
                # equilibrium PLUS wall-powered boiling â€” the
                # evaporation-generated-flow channel.
                area = self._areas[k]          # bed-side contact area
                gas_area = self._gas_areas[k]  # gas-side exchange area

                # Bed relaxation quantities FIRST: the energy row uses
                # the STEP-MEAN bed temperature of the exact linearized
                # relaxation (time-average of the exponential), not the
                # stale start-of-step value â€” at tau_T << dt the
                # unresolved warm-up transient would otherwise book a
                # spurious megawatt condensation against a cold bed.
                # Identical to the frozen-state convention at
                # stationarity and at tau_T >> dt.
                duty = self._layer_duty_w[k]
                # T3b series-flux energy closure: intraparticle
                # evaporation draws its latent from the PARTICLE'S OWN
                # heat — the bed — not from the gas film.  Charging it
                # here as a bed sink is what re-closes the layer energy
                # row after the sub-critical release override below
                # (unclosed it measured 968 kW on the six-tray machine).
                # RAW series flux, pre-truncation, disclosed: when the
                # truncation binds, the difference lands in the REPORTED
                # energy residual, never silently.
                # Bed-side scalars hoisted above the series block: the
                # T3B-N2 heat budget needs them before the relaxation
                # quantities are formed (plain reads, no state).
                up_t = t_bed[k + 1] if k + 1 < K else feed.T
                c_s = self._params[k].cp_dry_meal
                h_a = self._contact_h[k] * area
                series_override = surface is not None and surface_regime in (
                    "RECEDING_FRONT",
                    "DRY_SORBATE",
                )
                if series_override:
                    # GAS-SIDE DEPLETION on the series flux (exchanger
                    # NTU form).  j_supply is the per-area series flux
                    # at the INLET gas composition; along the bed the
                    # gas enriches and the driving force falls.  The
                    # NTU is built on the SERIES conductance
                    # 1/(1/k_intra + 1/k_gas):
                    #   tail (k_intra << k_gas): NTU ~ 0.02 measured at
                    #     SP1 -> area factor 0.99, behaviour unchanged;
                    #   near X_c (k_intra -> k_gas): NTU large -> the
                    #     release saturates at the gas carrying
                    #     capacity instead of j*A_full, which would
                    #     otherwise exceed it by orders AND crash the
                    #     bed via the latent sink below.
                    # generation_only nodes keep the full area: the
                    # bleed generates its own flow, nothing depletes.
                    k_intra_s = surface.intraparticle_conductance_m_s
                    k_ser = 1.0 / (
                        1.0 / max(k_intra_s, 1.0e-300) + 1.0 / k_gas
                    )
                    if generation_only:
                        area_series = area
                    else:
                        g_vol = gas_total / coefficients.film_density_kg_m3
                        ntu_s = k_ser * area / max(g_vol, 1.0e-12)
                        area_series = (
                            area * (-math.expm1(-ntu_s)) / ntu_s
                            if ntu_s > 1.0e-12
                            else area
                        )
                    r_series = j_supply * area_series
                    # --- T3B-N2 RATE SELECTOR (owner-RULED 2026-08-08,
                    # GT_PS2_T3B_N2_RATE_SELECTOR_RULING) -------------
                    # r_h = min(transport, heat): Cardarelli's diffusion
                    # flux vs Faner's heat-controlled capacity, selected
                    # by whichever BINDS.  min SELECTS a restriction; it
                    # never sums them, so the July mutual-exclusivity
                    # clause's double-count concern does not arise.  The
                    # heat budget is the bed's stationary sources at the
                    # start-of-step bed state (PASS-1 lag discipline,
                    # exact at stationarity): duty + advective sensible
                    # + SIGNED interface contact — an interface hotter
                    # than the bed (condensing steam) contributes; a bed
                    # hotter than the interface is spending that heat on
                    # the F3 surface row, so it is subtracted here.  One
                    # joule never serves both channels.
                    heat_budget_w = (
                        duty
                        + F * c_s * (up_t - t_bed[k])
                        + h_a * (T_I - t_bed[k])
                    )
                    r_heat = max(heat_budget_w, 0.0) / lam_h
                    if not self._screening_tier:
                        r_series = min(r_series, r_heat)
                    # CAP BEFORE CHARGING THE BED.  The raw series flux
                    # can exceed what the step bounds will let the layer
                    # release; charging the bed with the un-capped
                    # latent left a measured 307 kW unclosed on the
                    # six-tray machine (968 kW before the NTU
                    # correction).  Both bounds are computable here:
                    # the step-integrated particle drop and the exact
                    # relaxation inventory bound (lagged upstream, as
                    # everywhere in PASS 1).  PASS 2's re-truncation
                    # against the UPDATED upstream can still shave the
                    # release; that shaving is a step transient, zero
                    # at steady state, and lands in the REPORTED
                    # energy residual.
                    _relax_e = 1.0 - math.exp(-dt / tau[k])
                    _credit_e = (1.0 - _relax_e) / _relax_e
                    _up_e = x_h[k + 1] if k + 1 < K else feed.X_h
                    _drop_e = _si.supply_limited_loading_drop(
                        x_h[k],
                        t_sorb,
                        dt,
                        self._particle,
                        pressure_pa=pressure_k,
                    )
                    r_series = min(
                        r_series,
                        F * _up_e + F * x_h[k] * _credit_e,
                        self._holdups[k] * _drop_e / dt,
                    )
                    latent_series_w = r_series * lam_h
                    # capture the fully step-bounded transport value for
                    # the second (step-mean) heat pass: min chains
                    # commute, so re-minning against the t_mean budget
                    # below preserves every first-block cap.
                    r_transport = r_series
                else:
                    r_series = 0.0
                    latent_series_w = 0.0
                    r_transport = 0.0
                gain = F * c_s + h_a
                t_eq = (
                    F * c_s * up_t + h_a * T_I + duty - latent_series_w
                ) / gain
                tau_t = self._holdups[k] * c_s / gain
                decay = math.exp(-dt / tau_t)
                t_mean = t_eq + (t_bed[k] - t_eq) * (tau_t / dt) * (
                    1.0 - decay
                )
                if series_override and not self._screening_tier:
                    # T3B-N2, SECOND (step-mean) PASS of the heat cap —
                    # the owner's blessed form says the budget is taken
                    # at the STEP-MEAN bed temperature.  The first pass
                    # above used start-of-step t_bed (causality: t_mean
                    # needs the latent, which needs the cap); one Picard
                    # iteration closes the loop.  Identical at
                    # stationarity (t_mean == t_bed); off stationarity
                    # it damps the one-step lag that flipped the cap on
                    # and off — MEASURED on the quarter-duty wet
                    # machine, which converged pre-selector, then ran a
                    # fast 8-9-flips-per-23-steps oscillation around
                    # X ~ 0.072 (far from X_c, so not the physical
                    # bistability) under the start-of-step cap.
                    heat_budget_w = (
                        duty
                        + F * c_s * (up_t - t_mean)
                        + h_a * (T_I - t_mean)
                    )
                    r_heat = max(heat_budget_w, 0.0) / lam_h
                    r_series = min(r_transport, r_heat)
                    latent_series_w = r_series * lam_h
                    t_eq = (
                        F * c_s * up_t + h_a * T_I + duty
                        - latent_series_w
                    ) / gain
                    t_mean = t_eq + (t_bed[k] - t_eq) * (tau_t / dt) * (
                        1.0 - decay
                    )

                if generation_only:
                    D = 0.0
                    q_gas_raw = 0.0
                else:
                    ntu_m = film.mass_conductance_kg_m2_s * gas_area / gas_total
                    ntu_h = (
                        film.heat_conductance_w_m2_k
                        * gas_area
                        / (gas_total * self._cp_gas)
                    )
                    eff_m = (
                        gas_area * (-math.expm1(-ntu_m)) / ntu_m
                        if ntu_m > 1e-12
                        else gas_area
                    )
                    D = film.diffusive_hexane_flux_kg_m2_s * eff_m
                    t_out_raw = T_I + (gas_in.T - T_I) * math.exp(
                        -ntu_h * film.ackermann_factor
                    )
                    q_gas_raw = gas_total * self._cp_gas * (
                        gas_in.T - t_out_raw
                    )
                q_bed_in = h_a * (t_mean - T_I)
                lam_mix = w_i * lam_h + (1.0 - w_i) * lam_w
                R_t = (q_gas_raw + q_bed_in - D * (lam_h - lam_w)) / lam_mix
                if generation_only and R_t < 0.0:
                    # INERT NODE: no vapor present and a bed below its
                    # own interface state has nothing to condense.  With
                    # no vapor phase there is NO INTERFACE, so the
                    # interfacial heat path must close too â€” zeroing the
                    # mass flux alone would leave q_bed flowing into a
                    # surface that does not exist and break the layer
                    # energy row by exactly that heat (measured 125 kW
                    # at the zero-duty/zero-inlet S2 cell before this
                    # fix; found by the energy-residual diagnostic).
                    # The bed then exchanges only by advection and wall
                    # duty, so its relaxation drops the hÂ·A coupling.
                    R_t = 0.0
                    q_bed_in = 0.0
                    gain = F * c_s
                    t_eq = (F * c_s * up_t + duty - latent_series_w) / gain
                    tau_t = self._holdups[k] * c_s / gain
                    decay = math.exp(-dt / tau_t)
                raw_h = w_i * R_t + D
                raw_w = R_t - raw_h
                # --- T3b SERIES-FLUX hexane release (sub-critical) ----
                # Below X_c the hexane leaves by INTRAPARTICLE transport,
                # and the F3 machinery structurally cannot express it:
                # its hexane share is w_I*R_t + D, both evaluated at the
                # co-boiling interface composition, which at a sorbed
                # surface is nearly hexane-free — measured on the SP1
                # fixture: release 1e-3 kg/s against an intraparticle
                # supply of 2.3 kg/s, a 2000x under-release whose refill
                # ran the layer into X_c and fired a period-10
                # relaxation oscillation (entire inventory dumped in one
                # step, 47 K bed crash, repeat).  The series flux
                # k_intra*(c_mean - c_surface) IS the closure's
                # supply_limited_flux — the through-flux of the
                # two-resistance circuit, continuous at X_c (k_shell ->
                # inf recovers the gas-film-limited FILM_ACTIVE value)
                # and it is Cardarelli's diffusion-controlled residual
                # mechanism, the benchmark's own Table-10 requirement.
                # Water and heat stay with the F3 interface machinery;
                # the hexane latent of the override lands in the
                # REPORTED energy residual, never silently.
                #
                # The phi rate/heat candidate is deliberately NOT wired.
                # Its failed attempt left an unselected heat/water law and
                # a FILM_ACTIVE regime dispatch.  Retain this pre-candidate
                # behavior until a complete population closure is qualified.
                if series_override:
                    raw_h = r_series
                elif surface is None:
                    raw_h = 0.0          # EXHAUSTED: nothing to release

                # Exhaustion truncation (exact bookkeeping): release
                # bounded by what the EXACT relaxation update can drain
                # without a negative inventory (r <= FÂ·x_up +
                # FÂ·xÂ·(1â’relax)/relax); condensation bounded by the gas
                # inflow of that species.  PASS 2 re-truncates against
                # the updated upstream â€” the gas/solid mismatch this
                # leaves is a step transient, zero at cyclic steady
                # state (two-pass discipline, as in the sibling module).
                relax = 1.0 - math.exp(-dt / tau[k])
                credit = (1.0 - relax) / relax
                up_x_h = x_h[k + 1] if k + 1 < K else feed.X_h
                up_x_w = x_w[k + 1] if k + 1 < K else feed.X_w
                bound_h_out = F * up_x_h + F * x_h[k] * credit
                bound_w_out = (
                    F * up_x_w
                    + F * x_w[k] * credit
                    + free_water[k] / dt
                )
                # T3b PHYSICAL supply cap: a particle cannot deliver
                # n-hexane to its surface faster than its own
                # intraparticle transport allows.  This — not the
                # inventory bound — is what turns the residual into a
                # PLATEAU instead of a strip to zero, and it is the
                # resistance the film-only law omitted entirely.
                # STEP-INTEGRATED, not evaluated at the start-of-step
                # conductance: the regime can traverse within one macro
                # step, and at plant-scale area the start-of-step value
                # (infinite in regime A) would let a layer strip from
                # feed loading to zero in a single step.
                # MUST use the same sorption reference temperature the
                # closure itself used (t_sorb = t_bed), or the bound is
                # computed against a DIFFERENT X_e than the regime
                # dispatch, and the layer can be drained past the
                # transport-limited entry the bound exists to protect.
                drop_max = _si.supply_limited_loading_drop(
                    x_h[k],
                    t_sorb,
                    dt,
                    self._particle,
                    pressure_pa=pressure_k,
                )
                # The layer's particles cannot deliver more than their
                # own transport allows over the step.
                #
                # This is a bound on the RELEASE ONLY.  An earlier
                # version tried to route it through the exact relaxation
                # update as `r <= F(x_up - x + drop/relax)` so that the
                # loading itself could not fall below the transport-
                # limited entry.  That was wrong twice over: it can go
                # NEGATIVE when the upstream is lean, which would force
                # condensation, and it conflates two different
                # mechanisms — the particle's transport limit governs
                # how fast n-hexane leaves a PARTICLE, not how fast lean
                # particles replace rich ones in a well-mixed layer.  A
                # layer fed by a leaner upstream genuinely becomes
                # leaner; that is advection, and it is not bounded by
                # D_eff.
                #
                # The simple form does not over-drain either, which was
                # the stated reason for changing it.  With the release at
                # this bound and a steady upstream the exact update gives
                #   x_new = x - drop*(tau/dt)*(1 - e^(-dt/tau)) >= x - drop
                # because (tau/dt)(1 - e^(-dt/tau)) <= 1 for every dt.
                bound_h_supply = self._holdups[k] * drop_max / dt
                r_h_unlimited = min(max(raw_h, -gas_in.G_h), bound_h_out)
                r_h = min(r_h_unlimited, bound_h_supply)
                supply_limited = r_h < r_h_unlimited
                r_w = min(max(raw_w, -gas_in.G_w), bound_w_out)
                raw_mag = abs(raw_h) + abs(raw_w)
                factor = (
                    (abs(r_h) + abs(r_w)) / raw_mag if raw_mag > 0.0 else 1.0
                )
                j_h_specific[k] = max(r_h, 0.0) / self._holdups[k]
                rates.append((r_h, r_w))
                supply_drops.append(drop_max)

                # Heats carry the truncation factor (proportional
                # rule); the layer energy residual measures how far
                # asymmetric truncation broke the row â€” reported,
                # ~0 in the untruncated converged regime.
                q_gas_w = q_gas_raw * factor
                q_bed_w = -q_bed_in * factor  # interface -> bed sign
                if generation_only:
                    t_gas_out = T_I
                else:
                    t_gas_out = gas_in.T - q_gas_w / (
                        gas_total * self._cp_gas
                    )
                energy_residual_w = (
                    r_h * lam_h + r_w * lam_w
                ) - q_gas_w + q_bed_w - latent_series_w
                film_records.append(
                    FilmLayerStepRecord(
                        interface_temperature_k=T_I,
                        heteroazeotrope_pinned=pinned,
                        hexane_rate_kg_s=r_h,
                        water_rate_kg_s=r_w,
                        truncation_factor=factor,
                        gas_side_heat_w=q_gas_w,
                        bed_side_heat_w=q_bed_w,
                        energy_residual_w=energy_residual_w,
                        hexane_regime=surface_regime,
                        hexane_surface_activity=a_h_surface,
                        front_radius_fraction=(
                            surface.front_radius_m
                            / self._particle.particle_radius_m
                            if surface is not None
                            else 0.0
                        ),
                        supply_limited_flux_kg_m2_s=j_supply,
                        supply_limited=supply_limited,
                        wetted_fraction=wetted_phi,
                    )
                )

                # Bed temperature: exact linearized relaxation
                # (quantities precomputed above the energy row).
                t_bed[k] = t_eq + (t_bed[k] - t_eq) * decay

                # Gas node: exact mass update; generated vapor joins
                # the stream AT T_I (mass-weighted merge, the tranche
                # rule of merge_gas_streams).  FULL condensation is a
                # legal outcome (cold bed swallows the whole stream â€”
                # steam-up); the next layer then runs generation-only.
                g_h_out = max(gas_in.G_h + r_h, 0.0)
                g_w_out = max(gas_in.G_w + r_w, 0.0)
                g_total_out = g_h_out + g_w_out
                if g_total_out <= 0.0:
                    gas_in = GasStream(G_w=0.0, G_h=0.0, T=T_I)
                else:
                    delta_g = r_h + r_w
                    if delta_g > 0.0:
                        t_gas_out = (
                            gas_total * t_gas_out + delta_g * T_I
                        ) / g_total_out
                    gas_in = GasStream(G_w=g_w_out, G_h=g_h_out, T=t_gas_out)

            # PASS 2 (top-down): inventory updates with the implicit
            # outflow convention â€” each layer's inflow IS its upstream's
            # UPDATED state, so the solid-side ledger telescopes exactly.
            for k in range(K - 1, -1, -1):
                # PER-LAYER relaxation.  This was previously computed
                # ONCE above the loop, where `k` still held the leftover
                # value K-1 from PASS 1, so every layer was updated on
                # the LAST tray's time constant.  Harmless while all
                # holdups are equal (the single-tray and 3-layer
                # fixtures); wrong for the reference machine, whose bed
                # depths run 0.30 m to 1.10 m.
                relax = 1.0 - math.exp(-dt / tau[k])
                credit = (1.0 - relax) / relax
                r_h, r_w = rates[k]
                up_x_h = x_h[k + 1] if k + 1 < K else feed.X_h
                up_x_w = x_w[k + 1] if k + 1 < K else feed.X_w
                fw_inflow = free_water[k + 1] / tau[k + 1] if k + 1 < K else 0.0

                # Re-truncate against the UPDATED upstream (see PASS 1
                # note): keeps the exact update non-negative.
                r_h = min(r_h, F * up_x_h + F * x_h[k] * credit)
                # NOTE: the T3b particle-supply bound is NOT re-truncated
                # here.  It does not depend on the upstream state — it
                # bounds what this layer's own particles can deliver —
                # so the PASS 1 value already in `rates[k]` still holds.
                x_h_eq = up_x_h - r_h / F
                x_h_new = x_h_eq + (x_h[k] - x_h_eq) * (1.0 - relax)
                if x_h_new < 0.0:
                    raise FilmDepletionError(
                        f"negative hexane inventory at layer {k} step {step}"
                    )

                # Water: evaporation draws FREE water first (external
                # liquid before sorbed â€” F1 rule); condensate fills the
                # sorbed inventory to the Luikov cap, then spills to
                # free water.  Free water advects DOWN the solid path at
                # the holdup turnover rate; discharge leaves at L1.
                draw_free = 0.0
                if r_w > 0.0:
                    draw_free = min(r_w, free_water[k] / dt + fw_inflow)
                r_w_sorbed = r_w - draw_free
                r_w_sorbed = min(
                    r_w_sorbed, F * up_x_w + F * x_w[k] * credit
                )
                x_w_eq = up_x_w - r_w_sorbed / F
                x_w_new = x_w_eq + (x_w[k] - x_w_eq) * (1.0 - relax)
                spill = 0.0
                if x_w_new > cap:
                    spill = (x_w_new - cap) * self._holdups[k]
                    x_w_new = cap
                if x_w_new < 0.0:
                    raise FilmDepletionError(
                        f"negative water inventory at layer {k} step {step}"
                    )
                free_water[k] += dt * (
                    fw_inflow - free_water[k] / tau[k] - draw_free
                ) + spill
                if free_water[k] < 0.0:
                    free_water[k] = 0.0  # exact-zero handoff, no deadband

                step_release_h += (
                    F * (up_x_h - x_h_new)
                    - self._holdups[k] * (x_h_new - x_h[k]) / dt
                )
                step_release_w += (
                    F * (up_x_w - x_w_new)
                    - self._holdups[k] * (x_w_new - x_w[k]) / dt
                    + draw_free
                )
                x_h[k], x_w[k] = x_h_new, x_w_new

            total_released_h += step_release_h * dt
            fed_h += F * feed.X_h * dt
            discharged_h += F * x_h[0] * dt
            history.append(x_h[0])
            history_w.append(x_w[0])
            records.append(
                FilmDepletionStepRecord(
                    step=step,
                    layer_loadings_x_h=tuple(x_h),
                    layer_loadings_x_w=tuple(x_w),
                    layer_free_water_kg=tuple(free_water),
                    layer_bed_temperatures_k=tuple(t_bed),
                    hexane_release_kg_s=step_release_h,
                    water_release_kg_s=step_release_w,
                    films=tuple(film_records),
                )
            )

            if len(history) > STABILITY_LOOKBACK_STEPS:
                # WINDOW SPREAD, not a single lookback compare.  The old
                # criterion compared history[-1] against the value exactly
                # STABILITY_LOOKBACK_STEPS earlier — which a limit cycle
                # of that same period satisfies EXACTLY.  Measured on the
                # SP1 machine tray: a period-10 relaxation oscillation
                # (discharge X_h swinging 0.028 <-> 0.188, the whole
                # inventory dumped every 10th step) was reported
                # "converged" because 10 == the lookback.  The spread of
                # the whole window cannot be aliased by any period.
                tail = history[-1 - STABILITY_LOOKBACK_STEPS:]
                scale = max(abs(tail[-1]), 1.0e-12)
                tail_w = history_w[-1 - STABILITY_LOOKBACK_STEPS:]
                scale_w = max(abs(tail_w[-1]), 1.0e-12)
                # BOTH species must be stationary.  Watching x_h alone
                # declared convergence while discharge moisture was
                # still drifting at ~1e-5 per step — which is exactly
                # the mismatch the machine's steady-state water ledger
                # then reports (measured 0.3057 relative).
                if (
                    (max(tail) - min(tail)) / scale
                    < DISCHARGE_STABILITY_REL
                    and (max(tail_w) - min(tail_w)) / scale_w
                    < DISCHARGE_STABILITY_REL
                ):
                    converged = True
                    break

        if not converged:
            # Report WHY, not just that.  A monotone tail means the
            # budget was too short; an alternating tail means the
            # lagged-state feedback is oscillating and the cure is
            # relaxation or a smaller step — two different faults that
            # the bare message could not tell apart.
            tail = history[-min(len(history), 24):]
            deltas = [b - a for a, b in zip(tail, tail[1:])]
            sign_changes = sum(
                1 for a, b in zip(deltas, deltas[1:]) if a * b < 0.0
            )
            spread = (max(tail) - min(tail)) / max(abs(tail[-1]), 1e-30)
            character = (
                "OSCILLATING" if sign_changes >= max(2, len(deltas) // 3)
                else "drifting monotonically"
            )
            failure = FilmDepletionError(
                f"discharge did not stabilize within {max_steps} macro-steps "
                f"({character}: {sign_changes} sign changes in the last "
                f"{len(deltas)} increments, relative spread {spread:.3e} "
                f"against a {DISCHARGE_STABILITY_REL:.0e} criterion; "
                f"last discharge X_h = {tail[-1]:.6e}, "
                f"min {min(tail):.6e}, max {max(tail):.6e})"
            )
            # Carry the trajectory so a caller can see the SHAPE of the
            # failure, not just its magnitude.
            failure.discharge_history = tuple(history)
            failure.steps = tuple(records)
            raise failure

        window = records[-self._steps_per_residence :]
        avg_h = sum(r.hexane_release_kg_s for r in window) / len(window)
        avg_w = sum(r.water_release_kg_s for r in window) / len(window)
        max_energy_residual = max(
            abs(f.energy_residual_w) for r in window for f in r.films
        )

        final_inventory_h = sum(x * h for x, h in zip(x_h, self._holdups))
        ledger_residual = (
            fed_h
            - discharged_h
            - total_released_h
            - (final_inventory_h - initial_inventory_h)
        )
        ledger_rel = abs(ledger_residual) / max(fed_h, 1.0e-12)

        return FilmDepletionResult(
            steps=tuple(records),
            converged=True,
            residence_avg_hexane_release_kg_s=avg_h,
            residence_avg_water_release_kg_s=avg_w,
            discharge_x_h=x_h[0],
            discharge_temperature_k=t_bed[0],
            ledger_residual_rel=ledger_rel,
            max_energy_residual_w=max_energy_residual,
            energy_screening_tier=self._screening_tier,
        )


def _mole_fraction_hexane(gas: GasStream) -> float:
    """Stream mass rates -> hexane mole fraction (frozen molar masses)."""
    from .props import hexane as _hexane
    from .props import water as _water

    n_h = gas.G_h / _hexane.M
    n_w = gas.G_w / _water.M
    return n_h / (n_h + n_w)
