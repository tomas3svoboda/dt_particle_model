"""The ONE film-free surface n-hexane activity, shared by every caller.

Owner ruling 2026-09-06 night ("Ok, I trust the audit and your revised
recommendations. I rule as you ask.",
``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md``), on the peer-review audit
``docs/GT_PS2_PEER_REVIEW_AUDIT_THREE_RULINGS_2026-09-06.md`` section 3.3,
condition (i): the film-free surface activity is read through
``sorption_interface.surface_state`` via ONE shared helper that the sparge
kernel and the PD cold-surface condensation term both use, not through a
regime guard wrapped around ``props.sorption.activity_from_retained``.

WHAT THIS MODULE IS.  A thin restatement of the rule the sparge kernel
already runs (``sp1_k2_law2_tray_integration.regime_c_hexane_binding`` and
``_core_supplied_hexane_closure``), expressed on a LOADING rather than on a
layer inventory or a packet, so a caller that holds neither object can ask
the same question and get the same answer.  It computes nothing of its own:
every number comes from ``sorption_interface``, which is the authority.

    FILM_ACTIVE      a_h = 1.0 exactly (bulk liquid touches the surface)
    RECEDING_FRONT   a_h = 1.0 exactly (a wet core persists, so the
                     equilibrium surface activity is unity and ALL the
                     hindrance is transport, carried by a separate
                     conductance factor that is NOT this helper's business)
    DRY_SORBATE      the authority's own storage-based surface activity
                     ``surface_state(...).surface_activity``, capped at the
                     pressure-feasible ceiling
                     ``hexane_activity_ceiling(T, P) = min(1, P/Psat_h(T))``
    EXHAUSTED        refused typed - there is no surface state to report,
                     and nothing is clamped into existence

THE CONVENTION IS THE OWNER'S, NOT THIS MODULE'S.  In DRY_SORBATE the
authority reports ``a = c_surface / c_sat`` on the STORAGE the certified
``particle/dry_shell`` counts, i.e. saturated pore gas PLUS the PHY-043/048
GAB sorbate.  The sorbate-only inverse (the isotherm inversion the PD
cold-surface term used until today) reads higher, by a measured 0.4-3.6 % at
disc-tray temperatures and 7-23 % at kernel temperatures, because it assigns
the pore-gas share to the sorbate.  The audit put that choice to the owner
and it was ruled on the science lane's recommendation: the storage-based
surface concentration is the physics of record, because it is the kernel's
ruled choice and the one consistent with F6's definition of X_e and with
continuity at the regime B/C handover.

WHY NOT INSIDE ``sorption_interface`` ITSELF.  The authority is byte-pinned
for this packet, and the panel measured that it needs no new edge to serve
this rule: every caller already imports it.  This module therefore sits one
level above it and adds only the dispatch, so the kernel can adopt it in its
own later ruled edit without any authority moving.

physically_qualifying: false.  plant_predictive: false.  No F-gate advances.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from . import sorption_interface

#: Provenance class of every declared value in this module (D0).  There are
#: none: the only literals below are the exact 0.0 and 1.0 the regime map
#: itself states, which are not engineering choices.
DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"

#: Law id of the rule as a whole.  Bump the revision if a NUMBER or a rule
#: below changes meaning; adding a reported field is additive and keeps it.
FILM_FREE_HEXANE_ACTIVITY_LAW_ID = (
    "sorption-interface-regime-map-surface-activity-capped-at-feasible-ceiling-v1"
)

#: The kernel expression this helper mirrors, named so a reader can check it.
MIRRORED_KERNEL_LAW_SOURCE = (
    "src/dtdc_simulator/core2/sp1_k2_law2_tray_integration.py "
    "regime_c_hexane_binding + _core_supplied_hexane_closure"
)


class FilmFreeHexaneActivityError(RuntimeError):
    """The film-free surface n-hexane activity is not defined for this state."""


@dataclass(frozen=True, slots=True, kw_only=True)
class FilmFreeHexaneActivity:
    """The surface n-hexane activity of a film-free particle, and its basis."""

    #: Which of the certified regimes the loading is in.
    regime: sorption_interface.HexaneRegime
    #: What a caller binds at the interface: min(surface activity, ceiling).
    hexane_activity: float
    #: min(1, P / Psat_h(T)) - the pressure-feasible ceiling, reported on
    #: every regime so a caller can see the bound whether or not it binds.
    feasible_activity_ceiling: float
    #: True only when the ceiling actually moved the number.
    hexane_activity_capped: bool
    #: The unbounded authority value the cap was applied to (exactly 1.0 in
    #: FILM_ACTIVE and RECEDING_FRONT, where no cap is applied at all).
    surface_activity: float
    source_law_id: ClassVar[str] = FILM_FREE_HEXANE_ACTIVITY_LAW_ID
    declared_assumption_class: ClassVar[str] = DECLARED_ENGINEERING_ASSUMPTION

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def film_free_surface_hexane_activity(
    *,
    hexane_loading_kg_kg_dry: float,
    temperature_k: float,
    pressure_pa: float,
    residual_oil_mass_fraction: float,
    lagged_flux_kg_m2_s: float = 0.0,
) -> FilmFreeHexaneActivity:
    """The surface n-hexane activity of a particle carrying NO free film.

    ``hexane_loading_kg_kg_dry`` is the particle's whole non-film (pore plus
    sorbate) n-hexane inventory per kg dry matter - the conserved envelope
    loading ``inventory.internal_hexane_kg / inventory.dry_matter_kg``, which
    is the same quantity the kernel binds on.  A caller holding a free film
    must not call this at all: its surface is bulk liquid and its ratio is
    exactly 1.0 by the FILM_ACTIVE law, decided before any loading is read.

    ``lagged_flux_kg_m2_s`` defaults to 0.0, which the authority treats as the
    EQUILIBRIUM surface; positive is solid to gas, and a negative (condensing)
    flux raises the surface above the mean, which is where the ceiling can
    bind.  The default is the value both current callers use, for the reason
    the kernel records: no accepted per-area hexane flux is carried across
    macro-steps, so every binding in this build is the equilibrium surface.

    Raises :class:`FilmFreeHexaneActivityError` and never clamps: on an
    EXHAUSTED loading (nothing left to release, so there is no surface state),
    and on any refusal of the sorption authority itself.
    """

    try:
        params = sorption_interface.ParticleSorptionParams(
            residual_oil_mass_fraction=residual_oil_mass_fraction
        )
        regime = sorption_interface.classify(
            hexane_loading_kg_kg_dry,
            temperature_k,
            params,
            pressure_pa=pressure_pa,
        )
    except (sorption_interface.SorptionInterfaceError, ValueError) as error:
        raise FilmFreeHexaneActivityError(
            "the n-hexane regime of a film-free surface at "
            f"X = {hexane_loading_kg_kg_dry!r} kg/kg, T = {temperature_k!r} K, "
            f"P = {pressure_pa!r} Pa is refused by the sorption authority: {error}"
        ) from error

    if regime is sorption_interface.HexaneRegime.EXHAUSTED:
        raise FilmFreeHexaneActivityError(
            f"film-free n-hexane loading {hexane_loading_kg_kg_dry!r} kg/kg is "
            f"not positive at {temperature_k!r} K; there is no sorbed surface "
            "state to bind - reject, do not clamp"
        )

    try:
        ceiling = sorption_interface.hexane_activity_ceiling(
            temperature_k, pressure_pa
        )
        if regime is sorption_interface.HexaneRegime.DRY_SORBATE:
            surface = sorption_interface.surface_state(
                mean_loading_kg_kg=hexane_loading_kg_kg_dry,
                temperature_k=temperature_k,
                hexane_flux_kg_m2_s=lagged_flux_kg_m2_s,
                params=params,
                pressure_pa=pressure_pa,
            )
            surface_activity = surface.surface_activity
        else:
            # FILM_ACTIVE and RECEDING_FRONT: a liquid phase touches the
            # surface (at the outer face, or at the front behind a dry
            # shell), so the equilibrium activity is unity exactly and the
            # hindrance is transport, not activity.  No cap is applied -
            # the ceiling bounds the pore GAS a SORBATE is in equilibrium
            # with, and there is no sorbate branch here.
            surface_activity = 1.0
    except (sorption_interface.SorptionInterfaceError, ValueError) as error:
        raise FilmFreeHexaneActivityError(
            "the film-free surface n-hexane state at "
            f"X = {hexane_loading_kg_kg_dry!r} kg/kg, T = {temperature_k!r} K, "
            f"P = {pressure_pa!r} Pa is refused by the sorption authority: {error}"
        ) from error

    if regime is sorption_interface.HexaneRegime.DRY_SORBATE:
        activity = min(surface_activity, ceiling)
    else:
        activity = surface_activity
    return FilmFreeHexaneActivity(
        regime=regime,
        hexane_activity=activity,
        feasible_activity_ceiling=ceiling,
        hexane_activity_capped=activity != surface_activity,
        surface_activity=surface_activity,
    )
