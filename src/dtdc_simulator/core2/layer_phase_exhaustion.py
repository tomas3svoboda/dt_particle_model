"""B2 (speed program 2026-09-02): the closed-form Law-2 exhaustion constraint set.

The Law-2 layer depletion condition is a set of closed-form predicates — for
the three through-bed trays, two layers x {hexane, external water} each, so
twelve rows per tower.  ``demanded = fl(K * dt)`` with

    K = fl(fl(gas_side_active_area_m2 * molar_flux_mol_m2_s) * molar_mass)

taken from the CONVERGED fast block, which never reads the stride.  ``K`` is
therefore dt-free and ``dt* = stock / K`` is the analytic admissible stride.

WHY THIS MODULE EXISTS SEPARATELY, AND WHAT THAT COSTS
------------------------------------------------------
``sp1_k2_law2_tray_integration.py`` is a BYTE-PINNED certified source: the
RS-2 positive-solver control authority
(``qsc10_k4_rs2_positive_solver_control_authority.py``) pins its canonical-LF
length and SHA-256.  The obvious single-source-of-truth refactor — lifting the
kernel's own loop body into a shared helper — would re-key that pin, so it was
REDESIGNED AWAY.  The kernel gate at
``sp1_k2_law2_tray_integration._packet_entries_after`` is left byte-for-byte
frozen and keeps running as defense in depth; the branch below is a VERBATIM
transcription of it (same order, same hexane stock branch including the B1
stage-3c nil-film reclassification and the regime-C sorption floor, the same
C8 Tier-1 water branch - retained water above the qualified Luikov floor on a
film-free layer with a DECLARED sorbed-water arm, external film otherwise -
the same ``demanded >= 0.0 and demanded > stock`` gate, same message).

F-B2-1 (2026-09-06) is the measured cost of that arrangement: the C8 Tier-1
water branch landed in the kernel on 2026-09-05 and this transcription kept
naming ``external_water`` unconditionally, so a film-free armed layer with a
positive water flux was refused here for an "external_water stock 0.0" the
kernel would never have raised.  The arm is now threaded in (``sorbed_arms``)
and the branch mirrors the kernel again.

Because this is a transcription and not a shared call, byte-identity is held
by TEST rather than by construction:
``tests/test_core2_layer_phase_exhaustion_predicates.py`` constructs an
exhaustion at every through-bed tray and bit-compares the payload this module
raises against the payload the frozen kernel gate raises with the hoisted call
site disabled.  Any future edit to either branch that breaks the agreement
turns that test red.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, ClassVar

from . import cell_engineering_feasibility as ef
from . import k_cell_tray_host as tray_host
from . import qsc_layer_falling_rate_law as falling_rate
from . import sp1_k2_law2_tray_integration as legacy
from . import through_bed_sorbed_water_law as sorbed_water

LAYER_PHASE_EXHAUSTION_PREDICATE_LAW_ID = "sp1-k2-law2-layer-phase-exhaustion-closed-form-v1"
#: Byte-identical to the frozen kernel refusal (implicit concatenation there).
LAYER_PHASE_EXHAUSTION_MESSAGE = (
    "Law-2 flux would deplete an external packet phase; step rejected, not clamped"
)


@dataclass(frozen=True, slots=True, kw_only=True)
class LayerPhaseExhaustionPredicate:
    """One (layer, phase) row of the closed-form exhaustion constraint set."""

    vertical_layer_id: int
    phase_name: str
    layer_stock_kg: float
    demanded_kg: float
    layer_macro_step_s: float
    #: ``K`` in kg/s — the dt-free demand rate the converged fast block fixed.
    demand_rate_kg_s: float
    #: ``stock / K`` where the constraint is active (``K > 0``), else ``None``.
    admissible_macro_step_s: float | None
    exhausted: bool
    #: Option A (owner ruling 2026-09-28): the hexane class this row's stock and
    #: demand belong to on a layer whose transfer was composed per hexane class,
    #: else ``None`` - every aggregate row, and every water row.
    hexane_class_key: str | None = None
    #: Level 0b (owner ruling 2026-09-28, film-class water routing option (a)):
    #: the executable water class this row's stock and demand belong to on a
    #: layer whose transfer was composed per class, else ``None`` - every
    #: aggregate row, and every hexane row.
    water_class_key: str | None = None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def populated_layer_ids(tray: tray_host.AcceptedKCellTrayState) -> tuple[int, ...]:
    """Populated layers in the kernel's canonical (ascending) order."""

    return tuple(sorted({packet.owner_key.vertical_layer_id.value for packet in tray.packets}))


def hexane_depletion_stock(layer_total) -> tuple[str, float]:
    """The drawable hexane stock of a populated layer, by its phase.

    Mirror of ``sp1_k2_law2_tray_integration._hexane_depletion_stock`` (the
    kernel cannot import this module without a cycle), kept identical by
    test:

    * a film-bearing layer draws its attached film;
    * a receding-front layer (B1 stage 3c) draws the sorbed/pore inventory
      ABOVE the regime-C floor ``X_e(T)`` of the falling-rate law;
    * a dry-sorbate layer (TAG-3, owner ruling 2026-09-06, Q-C8-2B), mean
      loading at or below that floor, draws the WHOLE sorbed inventory: the
      GAB activity falls with the loading to the EXHAUSTED regime, which the
      core-supplied binder refuses typed at exactly zero.  Before this branch
      the stock of such a layer read NEGATIVE (inventory minus the floor), and
      the first profile-start march was refused at step 1 on that reading.
    """

    if layer_total.attached_hexane_kg != 0.0:
        return "attached_hexane", layer_total.attached_hexane_kg
    floor_loading = falling_rate.equilibrium_floor_loading(
        layer_total.solid_temperature_k,
        layer_total.residual_oil_label_kg / layer_total.dry_matter_kg,
    )
    mean_loading = layer_total.internal_hexane_kg / layer_total.dry_matter_kg
    if mean_loading <= floor_loading:
        return "internal_hexane", layer_total.internal_hexane_kg
    return "internal_hexane", layer_total.internal_hexane_kg - (
        floor_loading * layer_total.dry_matter_kg
    )


def water_depletion_stock(layer_total, arm) -> tuple[str, float]:
    """The drawable water stock of a populated layer, by its phase.

    Mirror of the C8 Tier-1 branch inside
    ``sp1_k2_law2_tray_integration._packet_entries_after`` (the kernel cannot
    import this module without a cycle), kept identical by test:

    * a film-bearing layer, and any layer whose model declares NO sorbed-water
      arm, draws its external film - every certified run;
    * a film-free layer WITH a declared arm draws the retained (sorbed) water
      ABOVE the accepted law's qualified floor ``W_ref``.  Sorbed water may
      evaporate down to that floor; below it the step is refused typed - the
      dry continuation is Tier 2, not a clamp.

    ``arm`` is the object the kernel reads through
    ``sp1_k2_law2_tray_integration._layer_sorbed_arm``; only its
    declared-or-not identity is used here, exactly as in the kernel.  The floor
    is the same symbol the kernel multiplies, not a retyped number.
    """

    if layer_total.external_water_kg != 0.0 or arm is None:
        return "external_water", layer_total.external_water_kg
    return "retained_water", layer_total.retained_water_kg - (
        sorbed_water.LUIKOV_PARAMS.W_ref * layer_total.dry_matter_kg
    )


def layer_phase_exhaustion_predicates(
    *,
    populated_layers: tuple[int, ...],
    inventories: dict[int, legacy.PacketDerivedLayerInventory],
    transfers: dict[int, legacy.LayerAcceptedTransfers],
    macro_steps: dict[int, ef.EngineeringCellMacroStep],
    sorbed_arms: dict[int, sorbed_water.SorbedWaterArm | None] | None = None,
) -> tuple[LayerPhaseExhaustionPredicate, ...]:
    """Evaluate the constraint set in the kernel's canonical order.

    Order is exactly the order the packet loop would encounter the refusals:
    populated layers ascending, hexane before the layer's water phase.

    ``sorbed_arms`` carries, per layer, the DECLARED sorbed-water arm the
    kernel reads through ``legacy._layer_sorbed_arm`` (or ``None`` where none
    is declared).  Its own default, ``None``, reads as "no layer declares an
    arm" - the pre-F-B2-1 behaviour, so every caller that does not pass it is
    bit-for-bit unmoved.
    """

    rows: list[LayerPhaseExhaustionPredicate] = []
    for layer in populated_layers:
        transfer = transfers[layer]
        layer_total = inventories[layer]
        step = macro_steps[layer]
        evaluation = step.fast.evaluation
        area = evaluation.gas_side_active_area_m2
        macro_step_s = step.accepted.macro_step_s
        hexane_by_class = getattr(transfer, "hexane_by_class", None)
        if hexane_by_class is None:
            hexane_phase_name, hexane_stock = hexane_depletion_stock(layer_total)
            hexane_rows = (
                (
                    None,
                    hexane_phase_name,
                    transfer.hexane_from_packets_kg,
                    hexane_stock,
                    area
                    * evaluation.hexane_molar_flux_mol_m2_s
                    * ef.CANONICAL_HEXANE_MOLAR_MASS_KG_MOL,
                ),
            )
        else:
            # Option A (owner ruling 2026-09-28): the kernel's gate reads one row
            # per hexane class, in HEXANE_CLASS_ORDER, each the unchanged stock
            # law on the class's own sub-inventory against the fsum of the
            # class's own transfers; the dt-free rate is the fsum of its
            # executable classes' own A_c J_c M.  The partition is the kernel's
            # (``legacy.layer_hexane_classes``); its exact-partition refusal
            # stays in the kernel gate, which runs after this hoist.
            hexane_rows = []
            for item in legacy.layer_hexane_classes(layer_total):
                members = [r for r in hexane_by_class if r.hexane_class_key == item.class_key]
                class_phase_name, class_stock = hexane_depletion_stock(item.inventory)
                hexane_rows.append(
                    (
                        item.class_key,
                        class_phase_name,
                        math.fsum(r.hexane_from_packets_kg for r in members),
                        class_stock,
                        math.fsum(
                            r.gas_side_active_area_m2
                            * r.hexane_molar_flux_mol_m2_s
                            * ef.CANONICAL_HEXANE_MOLAR_MASS_KG_MOL
                            for r in members
                        ),
                    )
                )
        arm = None if sorbed_arms is None else sorbed_arms.get(layer)
        if getattr(transfer, "water_by_class", None) is None:
            water_phase_name, water_stock = water_depletion_stock(layer_total, arm)
            water_rows = (
                (
                    None,
                    water_phase_name,
                    transfer.water_from_packets_kg,
                    water_stock,
                    area
                    * evaluation.water_molar_flux_mol_m2_s
                    * ef.CANONICAL_WATER_MOLAR_MASS_KG_MOL,
                ),
            )
        else:
            # Level 0b (owner ruling 2026-09-28): the kernel's gate reads one
            # water row per executable class, in the composing solve's order,
            # each the class's own demand against its own stock (its film, or
            # its exchangeable sorbed water summed over its packets at each
            # packet's own front); the dt-free rate is the class's own
            # A_c J_w,c M_w.  The classes are the kernel's own reader
            # (``legacy.layer_water_classes``) at the accepted layer pressure;
            # its retained-water check against the packets stays in the kernel
            # gate, which runs after this hoist.
            water_rows = tuple(
                (
                    item.class_key,
                    item.phase_name,
                    item.demand_kg,
                    item.stock_kg,
                    item.demand_rate_kg_s,
                )
                for item in legacy.layer_water_classes(
                    layer_total,
                    transfer,
                    arm,
                    pressure_pa=step.fast.state.layer_pressure_pa,
                )
            )
        for is_water, (class_key, phase_name, demanded, stock, rate) in (
            *((False, row) for row in hexane_rows),
            *((True, row) for row in water_rows),
        ):
            rows.append(
                LayerPhaseExhaustionPredicate(
                    vertical_layer_id=layer,
                    phase_name=phase_name,
                    layer_stock_kg=stock,
                    demanded_kg=demanded,
                    layer_macro_step_s=macro_step_s,
                    demand_rate_kg_s=rate,
                    admissible_macro_step_s=(stock / rate if rate > 0.0 else None),
                    exhausted=demanded >= 0.0 and demanded > stock,
                    hexane_class_key=None if is_water else class_key,
                    water_class_key=class_key if is_water else None,
                )
            )
    return tuple(rows)


#: S1 (speed program 2026-09-02): an OPT-IN, diagnostic-only channel that
#: publishes the complete constraint set at the moment of a refusal.
#:
#: The typed refusal carries only the FIRST failing row, because that is what
#: the frozen kernel gate raises and the payload must stay byte-identical.  A
#: caller that wants to REPLAY the refusal law in closed form (the march's
#: first-shot stride predictor) needs the whole set, including the rows that
#: did not fail and the dt-free ``demand_rate_kg_s`` of each.
#:
#: The default is ``None`` and the hook is read only inside the ``exhausted``
#: branch below, so on every path that does not refuse -- i.e. every path a
#: certified run takes -- this module is byte-for-byte the module it was
#: before: no call, no attribute read, no cost.  Nothing here can change a
#: refusal, its payload, or its order; an observer that raises would mask the
#: refusal, so observers must not raise (the march's does not).
PREDICATE_SET_OBSERVER: Callable[[tuple[LayerPhaseExhaustionPredicate, ...]], None] | None = None


def raise_first_layer_phase_exhaustion(
    predicates: tuple[LayerPhaseExhaustionPredicate, ...],
    *,
    interval_end_time_s: float | None,
) -> None:
    """Raise the FIRST failing predicate with the frozen kernel payload."""

    for row in predicates:
        if row.exhausted:
            observer = PREDICATE_SET_OBSERVER
            if observer is not None:
                observer(predicates)
            raise legacy.SP1K2Law2PhaseExhaustionError(
                LAYER_PHASE_EXHAUSTION_MESSAGE,
                vertical_layer_id=row.vertical_layer_id,
                phase_name=row.phase_name,
                layer_stock_kg=row.layer_stock_kg,
                demanded_kg=row.demanded_kg,
                layer_macro_step_s=row.layer_macro_step_s,
                interval_end_time_s=interval_end_time_s,
                hexane_class_key=row.hexane_class_key,
                water_class_key=row.water_class_key,
            )


__all__ = (
    "LAYER_PHASE_EXHAUSTION_MESSAGE",
    "LAYER_PHASE_EXHAUSTION_PREDICATE_LAW_ID",
    "PREDICATE_SET_OBSERVER",
    "LayerPhaseExhaustionPredicate",
    "layer_phase_exhaustion_predicates",
    "populated_layer_ids",
    "raise_first_layer_phase_exhaustion",
    "water_depletion_stock",
)
