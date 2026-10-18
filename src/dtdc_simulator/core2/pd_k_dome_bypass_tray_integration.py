"""Bounded PD1--PD3 K/RTD dome-bypass engineering integration.

This module is the minimum executable family-1 predessolventizer closure.  It
keeps the existing :mod:`k_cell_tray_host` as the packet/RTD authority, uses
one physical wall shared by every K layer, and treats the dome bypass and its
plenum as an external algebraic gas path.  Vapour never crosses a PD meal bed.

Only the attached-hexane ``high_loading_v1`` packet branch is admitted. External
free water is preserved as a nonvolatile, thermally coupled inventory; no PD
water-transfer law is selected. The bottom K layer receives indirect wall heat;
the two K layers exchange energy through the frozen effective bed conductivity;
when net heat
would raise a packet above pure n-hexane saturation at the fixed tray pressure,
the excess heat evaporates pure n-hexane isothermally.  Exhaustion of attached
hexane is an explicit variant-handoff rejection.  No Ergun law, perforated
floor, direct-steam source, pressure-loss correlation, plant calibration, or
production persistence is introduced here.  All qualification claims remain
false.
"""

from __future__ import annotations

import hashlib
import functools
import math
import os
import struct
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import ClassVar

from . import cell_closure as cc
from . import cell_native_caloric as native_caloric
from . import component_energy_datum_adapter as component_datum
from . import engineering_physical_k_event_scheduler as k_scheduler
from . import engineering_dual_variant_packet_adapter as dual_variant
from . import engineering_subcritical_packet_adapter as subcritical
from . import film_free_surface_hexane_activity as film_free_activity
from . import heteroazeotrope as azeotrope
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import qsc_meal_mechanical_work_deposition as meal_work
from . import through_bed_k_single_wall_thermohydraulic as through_k
from .dtdc_stack import REFERENCE_TRAYS, VaporContact, exchange_area_m2
from .props import binary_gas
from .props import hexane
from .props import sorption
from .props import water as water_properties

INTEGRATION_DIGEST_DOMAIN = "GT-PS-2/pd-k-dome-bypass-tray-integration/v1"
PD_REFERENCE_TRAY_IDS = ("PD1", "PD2", "PD3")
DEFAULT_RELATIVE_LIMIT = 1.0e-10
PACKET_CALLBACK_ID = "pd-k-dome-bypass-heat-controlled-pure-hexane-v1"
PACKET_CALLBACK_SOURCE = "source-pinned-family1-pd-engineering-closure"
GAS_ENTHALPY_CONVENTION = "declared-cp-plus-component-latent-and-common-datum-species-shifts-v1"
PD_K2_LAYER_COUNT = 2
# E-9 (CB-4) RULED 2026-08-23, Class A dedupe: the PD inter-K effective bed
# conductivity is single-sourced from the through-bed D-T2-1 locked authority
# instead of re-declaring the same physical constant.  Bit-identical value,
# same authority-id string; no digest moves.
PD_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K = through_k.THROUGH_BED_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K
PD_INTER_K_CONDUCTIVITY_AUTHORITY_ID = through_k.THROUGH_BED_INTER_K_CONDUCTIVITY_AUTHORITY_ID
PD_INTER_K_ALLOCATION_LAW_ID = "identical-representative-particle-share-v1"
PD_INTER_K_MODEL_SCHEMA_ID = "pd-k2-inter-k-conduction-model-v1"

#: F4 (2026-09-02) - THE PD BOILING AUTHORITY, DECLARED AND NOT SILENT.
#:
#: The certified two-liquid co-boiling law is implemented in full below
#: (``PDCoBoilingState``, ``_pd_co_boiling_state``,
#: ``_vapour_enthalpy_per_hexane_mole_j``, ``_water_mass_per_hexane_mass``, the
#: co-boiling packet solve, the co-boiled water rows of the packet, plenum and
#: system ledgers).  It is the physics the wave-2 audit requires: with two
#: immiscible liquids present the bed boils at the heteroazeotrope and the
#: vapour carries the certified water share.
#:
#: It is NOT the production default, and that is a DECLARED refusal to make a
#: Class-B decision unilaterally, not an oversight.  MEASURED, at 101 325 Pa,
#: on the shipped six-tray fixture:
#:
#:   * certified T_het          = 334.4803981593741 K
#:   * pure n-hexane T_sat      = 341.8645125361387 K  (+7.3841143767646 K)
#:   * certified vapour water   = 0.20516853343923214 mole fraction
#:   * shipped PD packet state  = 340.0 K, i.e. 5.5196018406259 K ABOVE the
#:     certified two-liquid equilibrium - an inadmissible state that flashes.
#:   * the flash of that superheat demands roughly TWICE the n-hexane the
#:     packet's attached film carries (attached 6.748107357117419e-08 kg
#:     against a sensible load of ~7.1e-3 J/K x 5.52 K), so it crosses the
#:     tag-1 -> tag-2 film-exhaustion boundary inside interval 1.
#:
#: Consequence, measured end to end: with the certified authority enabled the
#: four-interval tower march REFUSES at step 1
#: (``PD3 rejected: packet pd3:k1-rtd-packet-03 would exhaust attached hexane
#: in this PD step``) under BOTH ``codec_mode=off`` and ``codec_mode=dual``,
#: because the crossing has to be absorbed by the tag-2 subcritical packet
#: variant - B1 stage 3c, which is itself BLOCKED on an open Class-B owner
#: ruling.  Flipping this constant therefore halts the engineering delivery
#: lane until the PD population is re-initialized at the certified equilibrium
#: (or 3c lands), and that is an owner decision, not an agent one.
#:
#: Until it is ruled, the certified authority is still CONSULTED on every PD
#: step and its numbers are carried on ``ValidatedPDKTrayStep`` so the defect
#: is measured continuously and can never again be invisible.
PD_CERTIFIED_CO_BOILING_ENABLED = False
PD_BOILING_AUTHORITY_ID = (
    "pure-hexane-saturation-pending-owner-flip-to-certified-heteroazeotrope-v1"
    if not PD_CERTIFIED_CO_BOILING_ENABLED
    else "certified-heteroazeotrope-co-boiling-v1"
)

#: D12 LEVEL 0 - THE PRE-DESOLVENTIZER FLASH, ARMED ONLY BY ITS STUDY-LANE KNOB.
#:
#: The declared row ``pd_flash_level0`` of ``study_lane/pd_flash_level0.py``
#: (design ``docs/GT_PS2_D12_STATE_AND_LEVEL0_DESIGN_2026-09-28.md`` section 3;
#: TEMPORARY-PENDING-OWNER-REVIEW under the owner's "proceed fast" of
#: 2026-09-28).  ``False`` on disk.  ``study_lane.pd_flash_level0.install``
#: rebinds it, together with the F4 authority pair above, under
#: ``DTDC_PD_FLASH_LEVEL0=1``, and ``uninstall`` restores all three.  Every
#: ``if PD_FLASH_LEVEL0_ARMED`` block below is unreachable while it is
#: ``False``, so the knob-off step is the pre-edit step byte for byte.
#:
#: What the armed blocks do, and nothing else: a film-free (tag-2) packet is
#: admitted on the PD packet step under the dual-variant codec; the saturated
#: state a packet is solved at crosses to tag-2 through the dual adapter's
#: conserved-totals reconstitution when its attached film is exhausted at the
#: saturation temperature (the Q-PD-2 handoff, made reachable); the
#: energy-limited bisection's upper bound is the packet's whole inner hexane,
#: whose representable end is the tag-2 codec's own floor; a packet whose
#: energy would carry it past that floor refuses typed
#: (``PDFlashLevel0FloorRefusal``: tag-3 on a PD tray is outside level 0);
#: the co-boiled water share may be reduced by the knob module's source bound
#: (``PD_FLASH_LEVEL0_WATER_SHARE_FACTOR``, refused there above 1 percent);
#: and the candidate plenum outlet is recorded for that bound
#: (``PD_FLASH_LEVEL0_LAST_OUTLET``) before the certified guard reads it.
PD_FLASH_LEVEL0_ARMED = False
PD_FLASH_LEVEL0_WATER_SHARE_FACTOR = 1.0
PD_FLASH_LEVEL0_LAST_OUTLET: dict[int, tuple[object, ...]] = {}

#: THE TAG-2 FILM-EXHAUSTED CROSSING, LEVEL 0 (Q2 as ruled 2026-09-29), ARMED ONLY BY ITS
#: STUDY-LANE KNOB.  ``study_lane/pd_tag2_dry_shell_leg.py`` (declared row
#: ``post_exhaustion_dry_shell_leg_pd``; design
#: ``docs/GT_PS2_TAG2_CROSSING_LEVEL0_DESIGN_2026-09-28.md``) rebinds both under
#: ``DTDC_TAG2_FILM_EXHAUSTED_CROSSING_LEVEL0=1``.  ``False`` and ``None`` on disk, so
#: every ``if PD_TAG2_LEG_ARMED`` block below is unreachable and the knob-off step is
#: the pre-edit step byte for byte.  Armed, a tag-2 packet the leg has registered
#: releases the smaller of its dry-shell supply and the heat-limited flash (the hook's
#: ``choose``), and every packet step's result is recorded under its payload digest
#: (the hook's ``record``), after the packet ledger has passed.
PD_TAG2_LEG_ARMED = False
PD_TAG2_LEG_HOOK = None

#: EXPLORATION MODE, family F5 tag-3 freeze (owner ruling of 2026-10-13, record
#: ``docs/GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13.md``), ARMED ONLY BY ITS STUDY-LANE
#: KNOB.  ``study_lane/exploration_breach_ledger.py`` rebinds both under
#: ``DTDC_EXPLORATION_MODE=1`` and ``DTDC_EXPLORATION_FREEZE_TAG3=1``.  ``False`` and ``None``
#: on disk, so the armed branch before the level-0 floor refusal is unreachable and the
#: knob-off step is the pre-edit step byte for byte.  Armed, at the in-force and thermal-guard
#: packet calls only, a packet that would flash past the tag-2 floor is frozen at the floor
#: (a logged PHYSICS breach, written before the freeze acts) instead of refused.
PD_TAG3_FREEZE_ARMED = False
PD_TAG3_FREEZE_HOOK = None
PD_TAG3_FREEZE_CALL_KINDS = ("in_force", "thermal_guard")

#: EXPLORATION MODE, families F2 inter-K and F4 packet ledger (the same record).  ``None`` on
#: disk; ``study_lane/exploration_rams.py`` binds it when ``DTDC_EXPLORATION_RAM_INTERK=1`` or
#: ``DTDC_EXPLORATION_RAM_PACKET_LEDGER=1`` (with the master).  While it is None every
#: ``PD_EXPLORATION_HOOK`` clause below is unreachable and the step is the pre-edit step.
PD_EXPLORATION_HOOK = None

#: PD INTER-K CROSSING CAPACITY, THE UNIFORM-LAYER FORM (study-lane knob, ON by default:
#: the owner's ruling 2 of 2026-09-29, ``GT_PS2_OWNER_RULINGS_2026-09-29_CAMPAIGN_DAY3.md``).
#:
#: Ruling packet ``docs/GT_PS2_RULING_PACKET_PD_INTERK_GUARDS_2026-09-28.md``
#: option (iv), built behind this knob so that the ruling flips a default and
#: nothing else; RULED as the default by the owner on 2026-09-29.  The explicit
#: inter-K face heat must not exceed the room between the two layers.  The shipped form
#: measures that room as sum_i [E_i(T_other) - E_i,stored]: every packet of a
#: layer moved to the OTHER layer's mean, which also carries the energy that
#: evens out the layer's own packets (heat that never crosses the face).  On a
#: flashing tray, 16 packets at the lock beside one colder arrival, that
#: intra-layer term (377.7 J on the traced K2 PD3 state) swamps the true room
#: (62.2 J) and the guard refuses a face heat of 1.2e-4 J at every stride.  The
#: uniform-layer form measures the room of the layer taken as uniform at its
#: own mean, sum_i [E_i(T_other) - E_i(T_self)], exactly as the pinned
#: through-bed kernel builds the same guard (``_layer_uniform_energy_j``,
#: ``through_bed_k_single_wall_thermohydraulic.py``).  The face heat, its
#: allocation, the per-packet envelope guard, the layer-order guard and the
#: shaft-work bound are untouched, and the capacity enters no ledger and no
#: digest.  Rejecting outcomes, written first: a face heat larger than the
#: uniform-layer room still refuses; the form is refuted if a knob-off lane
#: changes a committed byte, or if a knob-on lane that passes it then refuses
#: on the per-packet envelope guard with a genuine face-transfer overshoot.
#:
#: Default: armed at import (the owner's ruling of 2026-09-29).  The off value is
#: ``DTDC_PD_INTERK_CAPACITY_UNIFORM_LAYER=0``; unarmed, the guard is the shipped
#: form byte for byte, and a knob-off lane is every record's identity twin.
#: ``PD_INTERK_CAPACITY_OBSERVER`` (None on disk) is the report channel: when it
#: is set, both forms are evaluated on every guard evaluation and handed to it,
#: so a lane reports what each form would have said; the form in force is the
#: knob's.
PD_INTERK_CAPACITY_UNIFORM_LAYER_ENVIRONMENT_VARIABLE = "DTDC_PD_INTERK_CAPACITY_UNIFORM_LAYER"
PD_INTERK_CAPACITY_UNIFORM_LAYER_ARMED = (
    os.environ.get(PD_INTERK_CAPACITY_UNIFORM_LAYER_ENVIRONMENT_VARIABLE, "1") != "0"
)
PD_INTERK_CAPACITY_OBSERVER: Callable[[Mapping[str, object]], None] | None = None

#: PD INTER-K ORDER GUARD, THE CONSISTENT LOCK RULE (study-lane knob, ON by default:
#: the owner's ruling 2 of 2026-09-29).
#:
#: The trace of the PD1 order-guard refusal (the guard build's record
#: ``GT_PS2_PD_INTERK_CAPACITY_KNOB_BUILD_RECORD_2026-09-28.md`` section 4) found a
#: mismatch of two rules inside one step, not a rounding tie: the shaft-work-free
#: check leaves an upper-layer packet that is already at the boiling lock 8.7e-9 K
#: ABOVE it (its face heat taken as sensible heat, because the branch selection
#: counts any sensible solve within ``temperature_tolerance`` - the owner's
#: 2026-09-06 floor of 1e-6 K - as "at the lock"), while the full step, whose
#: shaft work carries the same packet past that band, boils it back onto the lock
#: float.  The order guard then compares a lower layer pinned to the lock with an
#: upper layer the branch selection left above it, at 1.8e-12 K.
#:
#: The consistent rule, and nothing else: in the shaft-work-free check, a packet
#: whose sensible solve lands ABOVE the saturation temperature but inside the same
#: band the branch selection calls "at the lock" takes the boiling branch, exactly
#: as the full step does, so both checks use one lock rule.  The guard's own
#: tolerance, the 1e-6 K floor and its scope, the full step and every ledger are
#: untouched; the thermal-only entries are read only by the guards.  Rejecting
#: outcomes, written first: refuted if a knob-off lane changes a committed byte; if
#: a knob-on lane that passes the order guard then fails the per-packet envelope
#: guard or the boiling-saturation interval check on a thermal-only state; or if
#: the boiling branch of the check refuses where its sensible branch was accepted.
#: Armed at import by default (the owner's ruling of 2026-09-29); the off value
#: ``DTDC_PD_INTERK_ORDER_LOCK_CONSISTENT=0`` restores the shipped rule byte for
#: byte and gives every record's identity twin.  When
#: ``PD_INTERK_ORDER_LOCK_OBSERVER`` is set, both forms of the order guard are
#: computed and handed to it (``"guard": "order_lock"``), the form in force the
#: knob's.  The rows go to their own hook, ``PD_INTERK_ORDER_LOCK_OBSERVER`` (None
#: on disk), beside the capacity hook, whose one-row-per-evaluation contract stays.
PD_INTERK_ORDER_LOCK_CONSISTENT_ENVIRONMENT_VARIABLE = "DTDC_PD_INTERK_ORDER_LOCK_CONSISTENT"
PD_INTERK_ORDER_LOCK_CONSISTENT_ARMED = (
    os.environ.get(PD_INTERK_ORDER_LOCK_CONSISTENT_ENVIRONMENT_VARIABLE, "1") != "0"
)
PD_INTERK_ORDER_LOCK_OBSERVER: Callable[[Mapping[str, object]], None] | None = None

#: PD INTER-K ORDER GUARD, THE CONDUCTION-ONLY FORM (study-lane knob, ON by default:
#: the owner's ruling 2 of 2026-09-29 made the three guard forms the defaults; the
#: "wrong-way heat" row stays report-only, with no threshold, until the owner sets one).
#:
#: Criterion.  The layer-order guard protects the second-law ordering of explicit
#: conduction between the two layers: the inter-K face heat must not carry the
#: layers past their old equilibrium.  The shipped guard reads the count-weighted
#: means of the thermal-only entries, which carry EVERY source but shaft work
#: (floor heat, the steam condensing, the flash, and the inter-K heat), and books
#: the whole change of those two mixture statistics to conduction.  A mean that
#: moves through condensation and flash is not conduction: the trace of the
#: all-knobs production lane (A1, PD3 at 2.5 s) measured a "reversal" of
#: 1.59e-5 K, 8.8e6 times the guard's tolerance, set by the subcooled PD2
#: transfers warming towards the lock in both layers at unequal count shares,
#: while the face heat's own share of the swing is 3.9e-6 of it and no packet
#: sees a reversal.
#:
#: Value.  The pinned through-bed kernel's own form
#: (``through_bed_k_single_wall_thermohydraulic._temperature_after_internal_heat``
#: and its order check): each layer taken as uniform at its own mean, receiving
#: the inter-K face heat ALONE, inverted on its uniform-layer energy by bisection
#: inside the old layer-temperature envelope; the order is then checked with the
#: shipped guard's own tolerance (``_explicit_pair_crossed_equilibrium``).  The
#: energy is this module's ``_packet_energy_at_temperature_j`` (the level-0 dual
#: variant where that knob is armed), as in the uniform-layer crossing capacity.
#:
#: Bracket.  Off (``=0``): the shipped all-sources form, byte for byte.  On (default): the
#: conduction-only form decides; the all-sources reading is still computed and,
#: when ``PD_INTERK_ORDER_LOCK_OBSERVER`` is set, reported on every evaluation as
#: a report-only row (``"guard": "order_conduction_only"``) with the wrong-way
#: heat |Q| (1 - t_cross/dt) of an all-sources reversal and, on those
#: evaluations only, the exact counterfactual all-sources reading at zero face
#: heat, so a source-driven reversal of the mixture means stays visible.  Nothing
#: relaxed: the guard's tolerance, the crossing capacity (either form), the
#: per-packet envelope guard, the lock floor, the shaft-work bound and every
#: ledger are untouched, and the thermal-only entries still feed the full step.
#:
#: Rejecting outcomes, written first.  IDENTITY: refuted if a knob-off lane
#: changes a committed byte.  CONDUCTION: a face heat that itself reverses the
#: uniform layers' order beyond the tolerance still refuses, as does a face heat
#: whose bisection target leaves the old envelope ("would cross the old layer
#: equilibrium") or a caloric seed the codec refuses.  COUNTERFACTUAL: refuted if,
#: on a step it admits against an all-sources reversal, the zero-face-heat
#: all-sources reading is NOT itself reversed (conduction then carried the
#: crossing in the full composition).  CONSERVATION: every ledger as before.
#: WRONG-WAY HEAT: report only, per row, for the owner's materiality threshold.
#: Armed at import by default (the owner's ruling of 2026-09-29); the off value is
#: ``DTDC_PD_INTERK_ORDER_CONDUCTION_ONLY=0``, every record's identity twin.
PD_INTERK_ORDER_CONDUCTION_ONLY_ENVIRONMENT_VARIABLE = "DTDC_PD_INTERK_ORDER_CONDUCTION_ONLY"
PD_INTERK_ORDER_CONDUCTION_ONLY_ARMED = (
    os.environ.get(PD_INTERK_ORDER_CONDUCTION_ONLY_ENVIRONMENT_VARIABLE, "1") != "0"
)

#: D12 STEAM CONDENSING LEVEL 0 - THE RISING VAPOUR CONDENSING ON THE PD MEAL,
#: ARMED ONLY BY ITS STUDY-LANE KNOB (``study_lane/pd_steam_condensing_level0.py``,
#: ``DTDC_PD_STEAM_CONDENSING_LEVEL0``; TEMPORARY-PENDING-OWNER-REVIEW).  ``False``
#: on disk; the knob rebinds it in memory and restores it.  Design
#: ``docs/GT_PS2_D12_LEVEL0_BUILD_RECORD_2026-09-28.md`` section 4 (the bulk,
#: equilibrium-stage form its question 1 recommends: the vapour rising through the
#: tray leaves each packet it meets at most saturated at that packet's own surface).
#:
#: THE LAW.  The bypass stream arriving from below is shared over EVERY packet of
#: the tray (both K layers: "through the bed") by representative-particle count,
#: through the same exact-residual allocator as the inter-K heat.  On each packet
#: the certified D12-a fugacity criterion is evaluated against the packet's OWN
#: surface: y_i phi_i(T_gas, P, y) P against f_eq,i f_i^{L,pure}(T_packet, P), the
#: two legs ``_require_single_vapour_phase`` equates and the CS-1a per-packet
#: surface ratios f_eq (1 on a free film).  Where the vapour leg is larger the
#: component condenses, and it condenses EXACTLY the amount that brings the
#: packet's slice of the gas down to the packet's surface fugacity at the gas
#: temperature: an equilibrium quantity, no transfer coefficient, no area.  Water
#: first, the n-hexane held; then n-hexane "where the fugacity guard allows", the
#: water held at its post-condensation amount, down to the smaller of its own
#: surface equilibrium and the water leg of the D12-a guard AT THE GAS
#: TEMPERATURE (v2, 2026-09-28: v1 held each component at its inlet amount and, on
#: a surface below the two-liquid lock where both fire and no two-component slice
#: equilibrium exists, left the slice water-supersaturated at the gas temperature;
#: measured +772.9 Pa at the production point's PD2).  Each removal is also kept
#: inside the other component's guard leg, so every slice ends inside the D12-a
#: guard at its own temperature, and never past its surface equilibrium (the
#: bisection keeps the side whose fugacity is not below the surface's).  Each
#: packet is then capped by the EXISTING first-stop cap (its own energy to its own
#: first stop temperature, less the declared energy it already receives), and the
#: condensate is priced by the EXISTING F2 reverse crossing at the stream state, so
#: the gas leaves cooler in energy and poorer in mass by exactly what condensed.
#: The condensate joins the packet's free water (a free film on a hexane-wet
#: particle, the level-0b rule); a film-free packet dry of free water would take it
#: into its sorbed store through its dry shell, which the PD packet step does not
#: carry at level 0, so such a receiver REFUSES typed.  Everything downstream (the
#: packet step's condensate inputs, the level-0 flash branch that boils what the
#: latent supplies, the plenum sink, the host rows, the couplings' Q-CS1-R1 sinks)
#: is the shipped CS-1a plumbing, unchanged.  The CS-1a surface term (contact
#: fraction) must be off when this is armed.
PD_STEAM_CONDENSING_LEVEL0_ARMED = False
PD_STEAM_CONDENSING_LEVEL0_LAW_ID = (
    "d12-pd-steam-condensing-level0-slice-to-own-surface-fugacity-water-first-guard-limited-v2"
)
PD_STEAM_CONDENSING_LEVEL0_RECEIVER_LAW_ID = (
    "every-pd-packet-both-k-layers-own-representative-count-slice-of-the-bypass-gas-v1"
)
PD_STEAM_CONDENSING_LEVEL0_CAPACITY_LAW_ID = (
    "slice-equilibrium-at-gas-temperature-no-transfer-coefficient-v1"
)

#: D12 FOG-ONLY LEVEL 0 - THE LOWER END OF THE PD CONDENSATION BRACKET (study-lane knob).
#:
#: The owner's ruling 3 of 2026-09-29
#: (``docs/GT_PS2_OWNER_RULINGS_2026-09-29_CAMPAIGN_DAY3.md``): the D12 bracket, the
#: fog-only lower end built beside the steam-condensing upper end above, the physical
#: contact fraction between the two ends declared with the measured bracket, never
#: fitted.  ``False`` on disk.  ``study_lane.pd_fog_only_level0.install`` rebinds it
#: under ``DTDC_PD_FOG_ONLY_LEVEL0=1`` and ``uninstall`` restores it; every
#: ``if PD_FOG_ONLY_LEVEL0_ARMED`` block below is unreachable while it is ``False``.
#: The two knobs exclude each other: both environment variables set refuses typed when
#: this module is imported (``PDFogOnlyLevel0ConfigurationError``), and both flags armed
#: refuses typed in the step.
#:
#: THE LAW.  The vapour rising through a PD tray never touches the meal.  At the tray's
#: plenum it mixes adiabatically with the tray's own generated vapour and the wall's
#: heat (``_plenum_outlet_and_ledger``, unchanged).  Where that mixing leaves the
#: stream beyond its own water dew line by the certified D12-a criterion
#: (``PD_NO_FLASH_GUARD_CRITERION``), the excess water condenses INTO A FOG CARRIED BY
#: THE VAPOUR: an adiabatic flash at the stream's pressure to the temperature T_f at
#: which the remaining vapour sits on its certified water dew line while the vapour and
#: the fog liquid (priced by the F2 seam's common-datum liquid enthalpy) hold the
#: stream's own enthalpy.  The fog's mass and the latent it releases are booked in the
#: vapour stream: the ``PDBinaryGasStream`` the step books is the shipped one, vapour
#: plus fog at its energy-equivalent single-phase temperature, so every ledger,
#: coupling and tower balance reads exactly the mass and energy it reads without the
#: knob.  The phase split (T_f, the vapour's water mole fraction, the fog flow) is the
#: report-only ledger ``PDFogOnlyLevel0Ledger``.  Nothing is deposited on any packet; a
#: fog that meets a warmer or drier mix at the next plenum re-evaporates there exactly,
#: because the booked totals mix.  The certified guard ``_require_single_vapour_phase``
#: admits, at each of its three PD sites, exactly the water excess the fog ledger of
#: that very stream measured (the WS-A relief arithmetic, unchanged).  An n-hexane
#: excess (a hexane fog), a fog beside any condensation on the meal, a flash without a
#: bracketed root or a flash energy residual above the step's relative limit refuse
#: typed (``PDFogOnlyLevel0Refusal``).  The level-0 source bound is not consulted while
#: the fog form is armed: the fog takes the excess, so nothing is withheld at the
#: source.  Declared deviation: a consumer that reads a fog-bearing stream's
#: temperature as a gas temperature (the next tray's wall-to-gas exchange) reads the
#: energy-equivalent temperature, below T_f by the ledger's reported lift.
PD_FOG_ONLY_LEVEL0_ARMED = False
PD_FOG_ONLY_LEVEL0_ENVIRONMENT_VARIABLE = "DTDC_PD_FOG_ONLY_LEVEL0"
PD_STEAM_CONDENSING_LEVEL0_ENVIRONMENT_VARIABLE = "DTDC_PD_STEAM_CONDENSING_LEVEL0"
PD_FOG_ONLY_LEVEL0_LAW_ID = (
    "d12-pd-fog-only-level0-adiabatic-flash-to-own-certified-dew-line-booked-in-the-vapour-v1"
)
#: The report channel: the fog ledgers of the last evaluation of a request, keyed by
#: ``id(request)``, written only while armed and popped by the knob module's wrapper.
PD_FOG_ONLY_LEVEL0_LAST_LEDGERS: dict[int, dict[str, object]] = {}


#: WS-A (2026-09-03) - THE CONTACT-LIMITED DOME CONDENSATION CAPABILITY.
#:
#: Owner ruling R-1 / design packet
#: ``docs/GT_PS2_PD_PHASE_CHANGE_PROGRAM_DESIGN_2026-09-03.md`` S:2.1, with
#: decisions D3 (contact basis) and D5 (co-occupancy only) and Q-PD-3
#: (moisture-pickup ledger REPORT-ONLY) and Q-PD-4 (the capability flag flips
#: to True only at the consolidated re-baseline).
#:
#: THE CONTACT AREA is the EXISTING DOME_CONTACT area - the bed's top free
#: surface, ``dtdc_stack.exchange_area_m2(tray) == tray.cross_section_m2``
#: (28.274333882308138 m2 per PD tray) - scaled by a declared
#: ``condensation_contact_fraction`` EFFICIENCY multiplier.  This is D3, and it
#: is literature-grounded: US11661564B2 records that PD trays do not extend to
#: the vessel wall, so "only the top layer of the spent material will have a
#: possibility to have contact with the vapours" and "this contact is not
#: efficient".  NO new peripheral-rim area is introduced.
#:
#: THE RATE LAW is contact-limited film condensation on that surface, and its
#: transfer coefficient is NOT a new invented constant: it is the model's OWN
#: declared dome gas-film heat-transfer coefficient
#: ``wall.wall_to_bypass_gas_ua_w_k / A_contact`` converted to a molar
#: mass-transfer conductance by the Lewis relation at Le = 1,
#:
#:     k_G = h_gas / (P * cp_molar)        [mol / (m2 s Pa)]
#:     n_condensable = contact_fraction * A_contact * k_G * dp * duration
#:
#: with ``dp`` the water partial-pressure excess over its own pure saturation
#: pressure at the stream temperature - the SAME inequality
#: ``_require_single_vapour_phase`` evaluates, now driving a term instead of a
#: hard refusal.  ``cp_molar`` is the declared molar heat capacity on the
#: declared arm and the frozen mixture enthalpy's own analytic dT tangent on
#: the native arm, so neither arm borrows the other's numbers.
#:
#: THE CAP is equilibrium: never condense more than brings the stream exactly
#: to saturation.  Overshooting into an undersaturated stream would fabricate
#: liquid the thermodynamics does not ask for.
#:
#: THE ROUTING is D5 co-occupancy, and only that: condensate mass joins the
#: packets' ``external_water_kg`` and its energy joins
#: ``common_datum_energy_j`` through the F2 seam machinery
#: (``common_datum_liquid_enthalpy_j_mol`` + ``seam_latent_j_mol``, the
#: reverse crossing of the same interface F4's evaporation uses).  NO hexane
#: is removed by condensation - the Cardarelli-1998 kinetics trend and,
#: independently, Karnofsky and Cardarelli-Crapiste (2002) all contradict a
#: condensation-driven pore-clearing law, and the design packet S:1.4 forbids
#: building one.
#:
#: WHICH PACKETS receive it: the TOP K layer only.  The dome vapour touches
#: the bed's free top surface, not its interior - the same source sentence
#: that fixes the area fixes the receiver.  Heat enters this tray from below
#: (the floor) and condensate from above; they are geometric mirrors.
#:
#: THE DEFAULT IS INERT.  ``condensation_contact_fraction = 0.0`` takes the
#: pre-WS-A code path verbatim: the no-flash guard refuses exactly as before,
#: with the same message and the same type, no term is added to any ledger,
#: and the model definition digest does not move (the digest gains its
#: condensation component ONLY when the fraction is nonzero).
PD_CONDENSATION_CONTACT_FRACTION_DEFAULT = 0.0
#: CS-1a (2026-09-04, owner priority "resolve condensation under
#: not-supersaturated steam conditions").  THE DRIVING FORCE IS WRITTEN IN
#: FUGACITIES, PER COMPONENT, AGAINST THE COLD MEAL SURFACE:
#:
#:     DF_i = y_i phi_i(T_bulk, P, y) P  -  f_eq,i f_i^{L,pure}(T_s, P)
#:
#: with the vapour leg from the frozen PHY-053 cross-virial closure
#: (``props/binary_gas.state``) and the liquid leg from the frozen Helmholtz
#: authorities through the certified two-liquid module
#: (``heteroazeotrope._pure_liquid_fugacity``) - the SAME two functions the
#: certified co-boiling authority equates.  WS-A's own term was the water
#: partial-pressure EXCESS over pure saturation at the GAS temperature, i.e.
#: fog: it cannot fire in an unsaturated stream, which is the state the owner
#: named.  A pressure form ``p_w - f_eq Psat(T_s)`` is not an approximation
#: here.  CORRECTED 2026-09-04 (ceremony rev 1, S:11 item 2): INSIDE the
#: published 49-60 C feed bracket the two conventions agree in SIGN and differ
#: in MAGNITUDE - by a factor 1.4334906037345263 for water at the top of the
#: bracket - while the SIGN crossings lie just ABOVE it, water at 61.0 C and
#: hexane at 61.5 C, in opposite directions, inside the transit every march
#: makes from the published band to the 340.0 K fixture.  The reason this term
#: uses the fugacity convention is not the crossing but the authority: the
#: certified two-liquid solver equates these same two functions.
PD_CONDENSATION_MODEL_SCHEMA_ID = "pd-dome-cold-surface-condensation-model-v2"
#: CS-1a rev 2 (2026-09-04, priority brief S:6 finding F-CS1-2): the driving
#: force is evaluated PER PACKET of the receiving top layer, on that packet's
#: OWN surface temperature and its OWN equilibrium ratios, and the layer's
#: condensate is the SUM over packets.  Rev 1 evaluated it once on the
#: representative-weighted layer MEAN temperature, which averages a cold feed
#: packet into the warm resident population and destroys exactly the driving
#: force this term exists to represent.
PD_CONDENSATION_RATE_LAW_ID = (
    "lewis-analogy-le1-dome-gas-film-fugacity-driving-force-per-component-per-packet-v2"
)
#: The cap is the SURFACE at its own equilibrium - never gas saturation.  The
#: condensate energy may not carry the receiving layer past the stop
#: temperature of whichever condensing component stops first.
PD_CONDENSATION_CAP_LAW_ID = "surface-at-own-dew-point-first-stop-v1"
#: DECLARED_ENGINEERING_ASSUMPTION: one Lewis-analogy k_G for BOTH components.
#: No certified binary diffusivity authority exists in the property family, and
#: the Chilton-Colburn correction Le_i^(2/3) is an order-unity factor that can
#: never flip a sign or move a cap.  A later per-component correction is an
#: additive change with its own digest component.
PD_CONDENSATION_TRANSFER_COEFFICIENT_LAW_ID = (
    "lewis-analogy-le1-single-kg-both-components-v1"
)
#: The surface-equilibrium ratio law: a free film of a component gives
#: f_eq = 1 exactly (``sorption_interface`` FILM_ACTIVE, "bulk liquid touches
#: the surface"); dry of free water the packet's own modified-Luikov retained
#: activity is used on the disclosed fugacity-ratio convention; dry of
#: attached hexane the packet's own REGIME-MAP surface activity, read through
#: the one shared helper ``film_free_surface_hexane_activity`` on the packet's
#: own conserved envelope loading, oil weight fraction and temperature - unity
#: while a wet core persists, and the storage-based sorbate activity capped at
#: ``sorption_interface.hexane_activity_ceiling`` once the core is gone.
#: Every branch is read on the PACKET, never on a layer mean (F-CS1-2).
#:
#: CS-1a rev 2 (2026-09-04, verifier issues 3 and 4): rev 1 stamped
#: ``...-else-packet-own-sorption-ratio-v1`` while the hexane arm RAISED
#: instead of reading a sorption ratio, so a digest-bearing declaration
#: overstated the implemented physics for one of the two components.  The
#: hexane arm was then implemented, so the name was again what the code does,
#: and the id carried ``-per-packet-v2`` because the ratios are packet-own.
#:
#: REV 3, RE-KEYED 2026-09-06 (owner ruling of that night on the peer-review
#: audit ``GT_PS2_PEER_REVIEW_AUDIT_THREE_RULINGS_2026-09-06.md`` section 3.3;
#: record ``GT_PS2_Q_CS1_2_SHARED_ACTIVITY_HELPER_RECORD_2026-09-06.md``).
#: Rev 2's hexane arm inverted the GAB SORBATE isotherm on the packet's whole
#: internal inventory, which above X_e contains pore liquid that isotherm does
#: not represent: measured, it refused ALL of regime B (100 % of the film-free
#: loading interval above X_e; it could serve 9.7 % of (0, X_c] at 338 K and
#: 0.85 % at 378 K) and, where it did return a number, it double-counted the
#: pore-gas share of the storage by 0.4-3.6 % at disc-tray temperatures.  The
#: rule is now the regime map's own, so the id names the surface activity and
#: not the inverse.  The id is digest-bearing INSIDE the armed component only,
#: so every unarmed (contact fraction 0.0) digest is unmoved and every armed
#: digest moves once, here.
PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID = (
    "free-film-unit-else-packet-own-regime-map-surface-activity-per-packet-v3"
)
PD_CONDENSATION_CONTACT_AREA_BASIS_ID = (
    "dtdc_stack.exchange_area_m2/DOME_CONTACT-top-free-surface-cross-section-v1"
)
PD_CONDENSATION_ENERGY_ROUTING_ID = (
    "component_energy_datum_adapter.common_datum_liquid_enthalpy_j_mol"
    "+seam_latent_j_mol/reverse-crossing-v1"
)
#: CS-1a rev 2 (F-CS1-2): each top-layer packet receives the condensate ITS OWN
#: surface earned, not a representative share of a layer total.  The contact
#: AREA is still shared by representative-particle count through the same
#: exact-residual allocator, so the shares sum to the dome contact area exactly.
PD_CONDENSATION_RECEIVER_LAW_ID = "top-k-layer-per-packet-own-surface-driving-force-v2"
#: The literature acceptance band criterion A1 is checked against
#: (faithfulness assessment S:6.5, US8720082B2): PD meal moisture measured
#: FLAT at 10-12 % wet basis from PD inlet to PD outlet.  REPORTED, never
#: enforced - Q-PD-3 rules this ledger report-only.
PD_LITERATURE_MOISTURE_BAND_WET_BASIS = (0.10, 0.12)
PD_LITERATURE_MOISTURE_BAND_SOURCE = (
    "US8720082B2 via docs/GT_PS2_PD_TESTIMONY_FAITHFULNESS_ASSESSMENT_2026-09-03.md "
    "S:6.5 criterion A1 (PD meal moisture flat 10-12 % wb in-to-out)"
)
#: Criterion A2, same source: PD meal outlet temperature 60-68 C.
PD_LITERATURE_OUTLET_TEMPERATURE_BAND_K = (333.15, 341.15)

DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"

#: PHY-052's EXTRAPOLATION CEILING - the one place this repository types it.
#: RULED 2026-09-06 (owner, verbatim: "Ok 465 K is clear, proceed as
#: recommended", on the ledger entry "THE 465 K BRACKET CEILING: ADVERSARIAL
#: PASS DONE" in ``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md``).
#:
#: The FROZEN decision PHY-052 ("Native-soybean GAB extrapolation and caloric
#: plausibility policy", ``release/physics_decisions.yaml``, user_approved
#: 2026-07-30, packet A2b-owner-confirmation) authorises the
#: Cardarelli-Crapiste native-soybean GAB equilibrium - MEASURED over
#: 323.15-368.15 K - to be used up to 433 K, and only "with explicit
#: out-of-source-domain flags, propagated model-form uncertainty, derivative
#: and caloric plausibility gates, and a paper claim boundary that
#: distinguishes extrapolation from experimental validation".  Above 433 K the
#: frozen decision authorises NOTHING, so no surface-equilibrium law of this
#: build may be evaluated there.  Every consumer reads THIS constant; none
#: retypes the number and none takes it from an equation-of-state standoff
#: such as ``heteroazeotrope._LIQUID_REFERENCE_MAX_K``, which is a different
#: statement (where the pure-liquid reference stops existing) authored for a
#: different equality.
SORPTION_ISOTHERM_EXTRAPOLATION_CEILING_K = 433.0
SORPTION_ISOTHERM_EXTRAPOLATION_CEILING_SOURCE = (
    "release/physics_decisions.yaml PHY-052 (FROZEN, user_approved 2026-07-30, "
    "packet A2b-owner-confirmation), selection uncertainty_qualified_extrapolation: "
    "'Use the frozen native-soybean GAB to 433 K only with explicit "
    "out-of-source-domain flags, propagated model-form uncertainty, derivative "
    "and caloric plausibility gates, and a paper claim boundary that "
    "distinguishes extrapolation from experimental validation.  Do not refit "
    "without new data.'  The source measurement domain the decision names is "
    "323.15-368.15 K"
)

#: RULED 2026-09-06 (owner ruling on point A of
#: ``docs/GT_PS2_RULING_PACKET_THREE_MENU_POINTS_2026-09-06.md``, recorded
#: verbatim in ``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md`` under "OWNER
#: RULING 2026-09-06 ~11:40"): the declared materiality floor under this
#: step's own ``temperature_tolerance``.  BOTH driving-envelope guards
#: carry this ONE floor - the after-temperature check and the
#: condensation-free counterfactual check beside it.
#:
#: WHAT THE GUARDS PROTECT, unchanged.  An explicit inter-K transfer may
#: not carry a packet past the envelope spanned by its own old temperature
#: and its adjacent layer's old temperature: a monotonicity condition on
#: an explicit scheme.  No law, rate or gate is removed here; only the
#: width of the band inside which an excursion is called truncation rather
#: than a crossing moves, and above that band the guard still refuses.
#:
#: WHAT THE FLOOR ADMITS, measured and not argued.  F-GATE-1 (2026-09-06,
#: the twelve-interval reference fixture run with a 400-character
#: refusal-reason prefix, ``rskip_wide``): every refused attempt on that
#: fixture is this guard, with excursions of 3.7e-8, 3.9e-8, 4.2e-8 and
#: 4.6e-8 K at the full 0.25 s stride and 1.8e-8, 2.0e-8, 2.1e-8 and
#: 2.3e-8 K after one halving, against a printed tolerance of 1.709e-8 K,
#: on driving envelopes 1.2-2.1 mK wide.  The excursion is FIRST ORDER in
#: the stride - one halving halves it, the second halving is accepted - so
#: what this floor admits is the explicit scheme's TRUNCATION error, not
#: an instability.  1.0e-6 K is one part in a thousand of that measured
#: driving envelope, one part in a hundred thousand of the ratified 0.1 K
#: per-tray certificate bound
#: (``study_lane/trajectory_certificate.py:85-89``), and 1e-6 of the
#: published product-temperature band.  Extrapolating the measured
#: first-order growth, the guard would go inert on that fixture only above
#: strides of about 6 s, far above the 0.25 s the march runs at and the
#: 0.0625 s the batched design recommends.  A/B on the same fixture
#: (``rskip_floor6b``): 12 attempts for 12 intervals and zero refusals,
#: against 20 attempts and 8 refusals with the rung skip alone and 100 and
#: 88 with the walk; the six intervals that share their end times with the
#: reference are bit-identical in every tray total and temperature,
#: because the excursions occur only on the ladder intervals.  F-GATE-2 is
#: why the floor is placed on ``temperature_tolerance`` itself: the
#: 2026-09-05 floor reached the counterfactual guard only, so the guard
#: that actually walks the ladder carried no floor at all.
PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K = 1.0e-6
PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_SOURCE = (
    "owner ruling 2026-09-06 on the three-point packet; F-GATE-1 "
    "(overshoots 1.8e-8 to 4.6e-8 K, first order in the stride); F-GATE-2. "
    " The ruling is recorded verbatim in "
    "docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md under 'OWNER RULING "
    "2026-09-06 ~11:40' and reasoned in "
    "docs/GT_PS2_RULING_PACKET_THREE_MENU_POINTS_2026-09-06.md section 1; "
    "the two measurements are the ledger entries 'F-GATE-1' and 'F-GATE-1 "
    "A/B DONE (rskip_floor6b)'; the build record is "
    "docs/GT_PS2_TRANSFER_ENVELOPE_FLOOR_RECORD_2026-09-06.md"
)

#: F-CS1-1 correction (B), owner-ruled 2026-09-05.  The condensation-free
#: counterfactual below re-solves the packet THROUGH THE CODEC, so its
#: transfer-only temperature can only be as consistent with the packet's own
#: stored energy as the codec itself STATES that it is.  The tag-1 codec
#: states that consistency: it inverts one and the same stored (T, E) pair
#: twice, by the combined envelope bisection and by the native high-loading
#: caloric inversion, and admits the two answers to disagree by up to this
#: many kelvin before it refuses the state
#: (``high_loading_packet_adapter.py:1018``,
#: ``if abs(independently_inverted - temperature) > 2.0e-8``).  Every packet
#: kind that can reach the branch - a tag-1 packet issued, advanced, or
#: RECONSTITUTED as a boundary-split child - is issued through that one
#: inversion (``encode`` ``:502``, ``reconstitute`` ``:570``, ``advance``
#: ``:625``), so this is the stated
#: consistency of all of them.  Tag-2 states cannot reach the branch at all:
#: a packet with no attached hexane is refused on entry to
#: ``_advance_one_packet``.
#: RAISED 2026-09-06 from 2.0e-8 K to the single declared envelope floor
#: above, by the owner ruling on point A of the three-point packet: the
#: counterfactual floor is raised to the same 1.0e-6 K so both guards
#: carry ONE floor.  The 2.0e-8 K figure remains the tag-1 codec's own
#: stated (T, E) consistency and is still the NARROWER of the two
#: statements; it is superseded here because a counterfactual re-solve
#: held to a finer band than the real after-temperature check standing
#: beside it would refuse states that check admits.  The codec statement
#: is kept in the source note below and is now MEASURED (F-GATE-3, the
#: encode-invert round-trip resolution) rather than only quoted.
PD_COUNTERFACTUAL_TRANSFER_TOLERANCE_FLOOR_K = PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K
PD_COUNTERFACTUAL_TRANSFER_TOLERANCE_FLOOR_SOURCE = (
    "owner ruling 2026-09-06 on the three-point packet raised this floor "
    "from 2.0e-8 K to PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K = 1.0e-6 K so "
    "that both driving-envelope guards carry one floor; the superseded "
    "2.0e-8 K statement, still the narrower of the two, was "
    "src/dtdc_simulator/core2/high_loading_packet_adapter.py:1018 "
    "(HighLoadingV1PacketAdapter._invert_envelope: 'combined and native "
    "high-loading caloric inversions disagree' is raised only above 2.0e-8 K, "
    "so 2.0e-8 K is the stated (T, E) consistency of every tag-1 payload); "
    "the same inversion's bisection terminates at a bracket width of 1.0e-12 K "
    "(:1001) on the [T_min, T_max] = [310.0, 380.0] K bracket "
    "(:940-941, particle/wet_core.py:84-85), i.e. 70.0 / 2**48 = "
    "2.4868995751603507e-13 K, which is the finer of the two statements"
)

DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "sorption_isotherm_extrapolation_ceiling": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "PHY-052 authorises the GAB isotherm measured at 323-368 K "
                    "to 433 K only with out-of-source-domain flags"
                ),
                "value_k": SORPTION_ISOTHERM_EXTRAPOLATION_CEILING_K,
                "sources": SORPTION_ISOTHERM_EXTRAPOLATION_CEILING_SOURCE,
                "rejecting_outcome": (
                    "a stop-temperature root above it is refused typed for a "
                    "condensing component and recorded nan for a non-condensing "
                    "one"
                ),
            }
        ),
        "transfer_envelope_tolerance_floor": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "a materiality bound on the monotonicity of an explicit "
                    "inter-K transfer: one part in a thousand of the "
                    "1.2-2.1 mK driving envelopes measured on the reference "
                    "fixture (rskip_wide), one part in a hundred thousand of "
                    "the ratified 0.1 K per-tray certificate bound"
                ),
                "value_k": PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K,
                "sources": PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_SOURCE,
                "rejecting_outcome": (
                    "an explicit-transfer overshoot above 1.0e-6 K, which is "
                    "still refused and cured by the stride halving"
                ),
            }
        ),
        "counterfactual_transfer_tolerance_floor": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the smallest temperature disagreement the condensation-free "
                    "counterfactual may still call noise: a re-solve of the SAME "
                    "packet through the SAME codec cannot be more (T, E)-consistent "
                    "than the codec states its own stored pairs are, so the "
                    "transfer-only check is run at max(the step's own relative "
                    "temperature tolerance, this floor) and at nothing wider.  "
                    "Since the 2026-09-06 ruling this floor IS the one declared "
                    "envelope floor, so the counterfactual and the "
                    "after-temperature check are held to the same band"
                ),
                "value_k": PD_COUNTERFACTUAL_TRANSFER_TOLERANCE_FLOOR_K,
                "sources": PD_COUNTERFACTUAL_TRANSFER_TOLERANCE_FLOOR_SOURCE,
                "rejecting_outcome": (
                    "a transfer-only excursion above the floor on a packet with "
                    "zero declared heat and zero declared work is a defect of the "
                    "counterfactual, not noise; report it, never widen the floor "
                    "on that evidence.  The one widening this floor has had - to "
                    "PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K on 2026-09-06 - was "
                    "an owner ruling that gave the two guards one floor, not a "
                    "response to such an excursion; the 2026-08-30 witness of "
                    "1.3e-6 K with zero declared heats is carried as the "
                    "separate finding F-GATE-3 and is NOT excused by it"
                ),
            }
        ),
    }
)


class PDKDomeBypassIntegrationError(ValueError):
    """Base refusal for the bounded PD engineering integration."""


class PDKDomeBypassConfigurationError(PDKDomeBypassIntegrationError):
    """A tray, model, packet, or boundary declaration is incompatible."""


class PDKDomeBypassStepError(PDKDomeBypassIntegrationError):
    """A reversible PD trial step is physically or numerically inadmissible."""


class PDCoBoilingWaterExhaustedError(PDKDomeBypassStepError):
    """The certified co-boiling would drive external free water below zero.

    F4: the two-liquid heteroazeotrope exists only while BOTH liquids do.
    A packet whose external free water cannot supply the certified water
    share of the vapour has left this branch; it is refused, never clamped
    and never silently continued on a pure-hexane vapour.
    """


class PDFlashLevel0FloorRefusal(PDKDomeBypassStepError):
    """D12 level 0 CROSSING: a PD packet's flash would pass the tag-2 floor.

    The tag-2 codec's own floor (the a_h = 1 sorption gate W_eq(1, T)) is the
    end of the level-0 branch.  A packet whose energy would carry it below
    that floor would become a dry-sorbate (tag-3) packet on a PD tray, which
    level 0 does not admit (design
    ``GT_PS2_D12_STATE_AND_LEVEL0_DESIGN_2026-09-28.md`` section 3.2,
    CROSSING).  Refused typed, never clamped.
    """


#: D12 level 0: the refusals that shrink the energy-limited bracket search
#: when the flash is armed - the typed floor and every tag-2 codec refusal,
#: exactly as the tag-1 codec's refusals already shrink it.  Read only under
#: ``PD_FLASH_LEVEL0_ARMED``.
_LEVEL0_BRACKET_ERRORS = (PDFlashLevel0FloorRefusal, subcritical.SubcriticalPacketAdapterError)
#: D12 level 0: the tag-2 refusals of the sensible (no-boiling) attempt that
#: send an armed packet to the boiling branch, as tag-1's caloric refusal does.
_LEVEL0_SENSIBLE_ERRORS = (subcritical.SubcriticalPacketStateError,)


class StalePDKDomeBypassStateError(PDKDomeBypassIntegrationError):
    """Commit was attempted against a stale or foreign accepted state."""


class PDFogOnlyLevel0ConfigurationError(PDKDomeBypassConfigurationError):
    """D12 fog-only level 0 and D12 steam condensing level 0 were both configured."""


class PDFogOnlyLevel0Refusal(PDKDomeBypassStepError):
    """D12 fog-only level 0: a stream state the fog form cannot carry (never clamped)."""


def _require_fog_steam_exclusive(environ: Mapping[str, str]) -> None:
    """Refuse typed when both bracket ends are configured (import-time configuration)."""

    if (
        environ.get(PD_FOG_ONLY_LEVEL0_ENVIRONMENT_VARIABLE, "0") == "1"
        and environ.get(PD_STEAM_CONDENSING_LEVEL0_ENVIRONMENT_VARIABLE, "0") == "1"
    ):
        raise PDFogOnlyLevel0ConfigurationError(
            f"{PD_FOG_ONLY_LEVEL0_ENVIRONMENT_VARIABLE}=1 and "
            f"{PD_STEAM_CONDENSING_LEVEL0_ENVIRONMENT_VARIABLE}=1 are the two ends of the "
            "D12 bracket and exclude each other: the fog end deposits nothing on the meal, "
            "the condensing end takes every slice to the meal's own surface"
        )


_require_fog_steam_exclusive(os.environ)


def _require_float(
    name: str,
    value: float,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise PDKDomeBypassConfigurationError(f"{name} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise PDKDomeBypassConfigurationError(f"{name} must be strictly positive")
    if nonnegative and value < 0.0:
        raise PDKDomeBypassConfigurationError(f"{name} must be nonnegative")


def _float_bytes(value: float) -> bytes:
    return struct.pack(">d", value)


def _same_float(left: float, right: float) -> bool:
    return _float_bytes(left) == _float_bytes(right)


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _relative_residual(residual: float, *terms: float) -> float:
    return abs(residual) / max(1.0, *(abs(value) for value in terms))


PD_PACKET_LEDGER_ROWS_ENV = "DTDC_PD_PACKET_LEDGER_ROWS"
# GT-PS-2-OWNER-RULINGS-2026-10-10-REAUDIT-BATCH-01, item D6 (the PD2
# packet-ledger packet, questions 2 and 3): a report-only instrument.  With
# the knob at 1 a refused packet heat/mass/common-energy ledger appends its
# three rows (absolute residual, binary64 roundoff bound, relative residual),
# the maximum, the limit and the binding row to the unchanged refusal text, so
# the march's attempt record and its full refusal_text carry the numbers that
# were discarded at the raise before this build (audit sweep 2026-10-09,
# findings PT-3 and CC-4).  No predicate, bound or branch moves; with the knob
# off the text is byte-identical to the text before this build.
PD_PACKET_LEDGER_ROWS_DECLARATION = {
    "provenance_class": "DECLARED_ENGINEERING_ASSUMPTION",
    "ruling_id": "GT-PS-2-OWNER-RULINGS-2026-10-10-REAUDIT-BATCH-01",
    "item": "D6, questions 2 and 3 of the PD2 packet-ledger packet of 2026-10-09",
    "knob": PD_PACKET_LEDGER_ROWS_ENV,
    "default": "0",
    "effect": (
        "report-only: the refused packet ledger's rows appended to the "
        "unchanged refusal text; the march labels the refusal pd-packet-ledger "
        "and stamps the knob"
    ),
    "rejecting_outcome": (
        "any committed leaf or narration line of a knob-off lane differing "
        "from the same lane before this build"
    ),
}


PD_PACKET_LEDGER_STOCK_SCALED_ENV = "DTDC_PD_PACKET_LEDGER_STOCK_SCALED"
# GT-PS-2-OWNER-RULING-PD-PACKET-LEDGER-NORMALIZATION-2026-10-11-01, ruling 1
# (option C of docs/GT_PS2_RULING_PACKET_PD_PACKET_LEDGER_NORMALIZATION_2026-10-11.md):
# the common-energy row of the packet ledger divided its roundoff-stripped
# residual by the step's net terms alone (the packet's energy change, its
# vapour energy, its declared input), while the water and hexane rows divide
# by their stocks.  The residual is formed from the packet's common-datum
# energies before and after the step (of order 1e7 J, binary64 spacing about
# 2e-9 J), so one or two ulps of those stocks over a few joules of step terms
# read as 1e-10 and refused: the production point's 12.828 s (farm) and
# 13.035 s (workstation) terminals and the fixture's 18.2 s terminal, measured
# on 2026-10-10/11 with the rows instrument.  With the knob at 1 the two
# weighted stock energies join the energy row's denominator; the roundoff
# bound (which already counts the ulps of those stocks), the limit 1.0e-10
# and every predicate are unchanged.  Knob off: the former form, byte for byte.
PD_PACKET_LEDGER_STOCK_SCALED_DECLARATION = {
    "provenance_class": "DECLARED_ENGINEERING_ASSUMPTION",
    "ruling_id": "GT-PS-2-OWNER-RULING-PD-PACKET-LEDGER-NORMALIZATION-2026-10-11-01",
    "knob": PD_PACKET_LEDGER_STOCK_SCALED_ENV,
    "default": "0",
    "former_form": (
        "energy relative residual = max(0, |r| - bound) / max(1 J, |dE_packet|, "
        "|E_vapour|, |E_declared_input|)"
    ),
    "ruled_form": (
        "energy relative residual = max(0, |r| - bound) / max(1 J, |dE_packet|, "
        "|E_vapour|, |E_declared_input|, |E_after| * n, |E_before| * n)"
    ),
    "measured_floor": (
        "2.0e-9 to 6.4e-9 J of either sign on 1e7 J stocks, horizon-independent, "
        "on both instrumented lanes"
    ),
    "limit_unchanged": 1.0e-10,
    "rejecting_outcome": (
        "any knob-on commit whose stock-scaled residual exceeds 1.0e-10, or any "
        "knob-on leaf before the first former-form refusal differing from the "
        "knob-off leaf, or any knob-off leaf differing from the lane before this build"
    ),
}


def packet_ledger_stock_scaled_enabled(environ=None) -> bool:
    """True only when DTDC_PD_PACKET_LEDGER_STOCK_SCALED is exactly "1"."""
    env = os.environ if environ is None else environ
    return env.get(PD_PACKET_LEDGER_STOCK_SCALED_ENV, "0") == "1"


def packet_ledger_rows_enabled(environ=None) -> bool:
    """True only when DTDC_PD_PACKET_LEDGER_ROWS is exactly "1"."""
    env = os.environ if environ is None else environ
    return env.get(PD_PACKET_LEDGER_ROWS_ENV, "0") == "1"


def packet_ledger_refusal_text(packet_id, rows=None, maximum=None, limit=None, former_energy_relative=None) -> str:
    """The packet-ledger refusal text.  ``rows`` is None (the text before this
    build, byte for byte) or three tuples (name, absolute residual, roundoff
    bound, relative residual) in the ledger's order water, mass, common-energy;
    the suffix never changes the leading clause."""
    base = f"packet {packet_id} heat/mass/common-energy ledger failed"
    if rows is None:
        return base
    binding = max(rows, key=lambda r: r[3])[0]
    parts = [
        f"{name} residual {absolute:.6e} (roundoff bound {bound:.6e}, relative {relative:.6e})"
        for name, absolute, bound, relative in rows
    ]
    return (
        base
        + "; rows: "
        + "; ".join(parts)
        + f"; maximum relative {maximum:.6e} against limit {limit:.3e}; binding row {binding}"
        + (
            ""
            if former_energy_relative is None
            else f"; former-form energy relative {former_energy_relative:.6e}"
        )
    )


def _roundoff_aware_relative_residual(
    residual: float,
    roundoff_bound: float,
    *terms: float,
) -> float:
    """Remove only an explicit binary64 operation bound, never physics error."""

    return max(0.0, abs(residual) - roundoff_bound) / max(1.0, *(abs(value) for value in terms))


def _ulp_sum(*values: float) -> float:
    return math.fsum(math.ulp(value) for value in values)


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise PDKDomeBypassConfigurationError(f"{name} must be one lowercase sha256 identity")


@dataclass(frozen=True, slots=True, kw_only=True)
class PDBinaryGasStream:
    """One no-inert bypass/plenum stream on the declared common datum."""

    pressure_pa: float
    molar_flow_mol_s: float
    hexane_mole_fraction: float
    water_mole_fraction: float
    temperature_k: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    authenticated_instrumentation: ClassVar[bool] = False
    flash_or_condensation_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_float("gas pressure", self.pressure_pa, positive=True)
        _require_float("gas molar flow", self.molar_flow_mol_s, positive=True)
        _require_float("gas hexane mole fraction", self.hexane_mole_fraction)
        _require_float("gas water mole fraction", self.water_mole_fraction)
        _require_float("gas temperature", self.temperature_k, positive=True)
        if not 0.0 <= self.hexane_mole_fraction <= 1.0:
            raise PDKDomeBypassConfigurationError("gas hexane mole fraction lies outside [0, 1]")
        if not 0.0 <= self.water_mole_fraction <= 1.0:
            raise PDKDomeBypassConfigurationError("gas water mole fraction lies outside [0, 1]")
        if math.fsum((self.hexane_mole_fraction, self.water_mole_fraction)) != 1.0:
            raise PDKDomeBypassConfigurationError(
                "the bounded bypass is exactly binary and contains no inert species"
            )


#: THE NO-FLASH GUARD'S CRITERION, and the authority that changed it.
#:
#: D12-a RULED by the owner on 2026-09-16
#: (``docs/GT_PS2_OWNER_RULINGS_D12_REISSUE_2026-09-16.md``, verbatim "I rule
#: D12-e yes as a brief and D12-a yes as a build"; the packet
#: ``docs/GT_PS2_RULING_PACKET_D12_REISSUE_2026-09-16.md`` item D12-a; the
#: queued Tier-B item Q-GUARD-1 option (b),
#: ``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md`` about :1559): the dome's
#: single-phase guard tests the CERTIFIED fugacity criterion, per component,
#:
#:     y_i phi_i(T, P, y) P  >  f_i^{L,pure}(T, P)          -> refuse
#:
#: and no longer the ideal/Raoult closure ``y_i P > Psat_i(T)`` it carried
#: before.  The form is the one ALREADY SHIPPING in the PD family at
#: ``qsc_pd_complete_contact_throughflow._water_phase_residual_pa`` and
#: ``impermeable_water_equilibrium_mole_fraction`` (:899-915): the vapour leg
#: from the frozen PHY-053 cross-virial closure ``props.binary_gas.state`` at
#: the module default ``k_wh = K_WH_CENTRAL``, the liquid leg from
#: ``heteroazeotrope._pure_liquid_fugacity``, which VERIFIES its own liquid
#: root and refuses where the pure-liquid reference state does not exist.
#:
#: NO NEW NUMBER IS DECLARED HERE.  The change is a criterion, not a value:
#: both legs are frozen authorities this module already imports, so P-1 has no
#: literal to type.  What changed is authorised by the ruling record above and
#: by nothing else.
#:
#: MEASURED, before the build (``docs/GT_PS2_PD_CHAIN_MUTANT_MEASUREMENT_RECORD_2026-09-15.md``
#: S:2, S:4 and the guard traces of its evidence folder): identity on every
#: shipped lane - 0 of 1,203 fixture evaluations and 0 of 2,754 production
#: leaves change - and two-sided where the state boils, admitting two PD3
#: plenum n-hexane states the ideal closure refused and refusing one PD3
#: plenum water state the ideal closure admitted.
PD_NO_FLASH_GUARD_CRITERION = "y_i phi_i P > f_i^{L,pure}(T,P) -> refuse"
PD_NO_FLASH_GUARD_CRITERION_AUTHORITY = (
    "D12-a, Q-GUARD-1 option (b), RULED by the owner 2026-09-16 in "
    "docs/GT_PS2_OWNER_RULINGS_D12_REISSUE_2026-09-16.md on the packet "
    "docs/GT_PS2_RULING_PACKET_D12_REISSUE_2026-09-16.md; the shipped form is "
    "props.binary_gas.state at K_WH_CENTRAL against "
    "heteroazeotrope._pure_liquid_fugacity, as "
    "qsc_pd_complete_contact_throughflow.py:899-915 already carries it"
)


def _require_single_vapour_phase(
    stream: PDBinaryGasStream,
    *,
    label: str,
    admitted_supersaturation_pa: tuple[tuple[str, float], ...] = (),
) -> None:
    """Fail closed where the no-flash algebra would require condensation.

    WS-A, NARROWED at the 2026-09-04 verification pass.  The refusal is NOT
    retired by ARMING the capability.  It is relieved for one species only up
    to a POSITIVE, MEASURED bound that an existing
    ``PDContactCondensationLedger`` produced for this very step - the
    ledger's own driving excess at the stream it was built from, and the
    ledger's own reported ``residual_supersaturation_pa`` at the plenum
    streams the step leaves behind.  Any supersaturation BEYOND that bound is
    still refused, verbatim and with the same type.  A step with no ledger -
    the declared default ``condensation_contact_fraction = 0.0``, and equally
    any armed state whose term does not actually fire - passes ``()``.

    D12-a, RULED 2026-09-16: the inequality is the CERTIFIED one,
    ``PD_NO_FLASH_GUARD_CRITERION`` above - partial FUGACITY against the pure
    compressed-liquid fugacity, not partial PRESSURE against Psat.  Everything
    else is unchanged: the species order (n-hexane then water), the
    at-or-above-critical skip, the relief's arithmetic and its exclusive
    ``> 0.0`` admission, the exception type, and the refusal message verbatim.
    The relief's bound is whatever the step's own ledger measured; the excess
    it is compared against is now the certified one, which is the smaller of
    the two at every state the PD lanes have visited, so a bound earned under
    the old convention relieves at least as much as it did.

    WHERE THE CERTIFIED REFERENCE STATE DOES NOT EXIST - hexane's compressed
    liquid collapses to the vapour branch above about 471 K at 1 atm, below
    its critical temperature, so the skip above does not cover it - this guard
    FAILS CLOSED with its own type rather than falling back to a convention
    the owner did not rule.

    THIS IS NOT THE CONDENSATION CRITERION.  The CS-1 ceremony
    (``docs/GT_PS2_CS1_COLD_SURFACE_CONDENSATION_CEREMONY_2026-09-04.md``)
    separates them: what condenses on a cold surface is decided by the CS-1
    term against the surface state, while this guard decides only whether the
    BULK stream the no-flash algebra produced can exist as one phase.
    """

    admitted_excess: dict[str, float] = {}
    for admitted_species, bound in admitted_supersaturation_pa:
        _require_float(
            f"admitted {admitted_species} supersaturation",
            bound,
            nonnegative=True,
        )
        admitted_excess[admitted_species] = bound

    gas_state = None
    refusals: list[tuple[str, float]] = []
    for species, key, critical_temperature, properties_module in (
        ("n-hexane", "hexane", hexane.TC, hexane),
        ("water", "water", water_properties.TC, water_properties),
    ):
        if stream.temperature_k >= critical_temperature:
            continue
        if gas_state is None:
            gas_state = binary_gas.state(
                stream.temperature_k,
                stream.pressure_pa,
                stream.water_mole_fraction,
                stream.hexane_mole_fraction,
            )
        partial_fugacity = (
            gas_state.fugacity_hexane_pa
            if key == "hexane"
            else gas_state.fugacity_water_pa
        )
        try:
            pure_liquid_fugacity = azeotrope._pure_liquid_fugacity(
                properties_module, stream.temperature_k, stream.pressure_pa
            )
        except azeotrope.HeteroazeotropeError as error:
            raise PDKDomeBypassStepError(
                f"{label} the certified {species} fugacity criterion could not "
                f"be evaluated: {type(error).__name__}: {error}"
            ) from error
        excess = math.fsum((partial_fugacity, -pure_liquid_fugacity))
        if excess > 0.0:
            refusals.append((species, excess))
    for species, excess in refusals:
        admitted = admitted_excess.get(species, 0.0)
        if admitted > 0.0 and excess <= admitted:
            # Relieved, to exactly the extent this step's own condensation
            # ledger measured - never further.
            continue
        raise PDKDomeBypassStepError(
            f"{label} {species} partial pressure exceeds pure saturation; "
            "the bounded no-flash vapour branch rejects condensation"
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDFogOnlyLevel0Ledger:
    """D12 fog-only level 0: the phase split of one booked PD stream that carries fog.

    Report only (``PD_FOG_ONLY_LEVEL0_ARMED``).  The booked stream holds the vapour
    plus the fog at its energy-equivalent single-phase temperature; this ledger states
    what that stream physically is: vapour on its certified water dew line at
    ``flash_temperature_k`` and a water fog at the same temperature, together holding
    the booked molar enthalpy.  Rates are per second of the booked flow.
    """

    label: str
    law_id: str
    pressure_pa: float
    booked_temperature_k: float
    booked_molar_flow_mol_s: float
    booked_hexane_mole_fraction: float
    booked_water_mole_fraction: float
    booked_enthalpy_j_mol: float
    water_excess_pa: float
    hexane_excess_pa: float
    flash_temperature_k: float
    vapour_hexane_mole_fraction: float
    vapour_water_mole_fraction: float
    vapour_water_excess_pa: float
    vapour_molar_flow_mol_s: float
    fog_water_molar_flow_mol_s: float
    fog_water_mass_flow_kg_s: float
    fog_liquid_enthalpy_j_mol: float
    flash_energy_residual_j_mol: float
    flash_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def temperature_lift_k(self) -> float:
        """T_f less the booked energy-equivalent temperature."""

        return self.flash_temperature_k - self.booked_temperature_k


def _fog_dew_line_hexane_fraction(
    temperature_k: float,
    pressure_pa: float,
    lower_hexane_fraction: float,
    pure_water_fugacity_pa: float,
) -> float:
    """The smallest n-hexane fraction in [lower, 1] on or inside the certified water dew line.

    The vapour pair is ``(1.0 - h, h)``, the complementary form the plenum books; the
    water leg is the guard's own (``binary_gas.state`` at the module default ``k_wh``
    against the pure compressed-liquid fugacity).  Bisection to adjacent floats; the
    returned end always satisfies the criterion (excess <= 0), never the other side.
    """

    def excess(hexane_fraction: float) -> float:
        gas = binary_gas.state(temperature_k, pressure_pa, 1.0 - hexane_fraction, hexane_fraction)
        return math.fsum((gas.fugacity_water_pa, -pure_water_fugacity_pa))

    low, high = lower_hexane_fraction, 1.0
    if excess(low) <= 0.0:
        return low
    if excess(high) > 0.0:
        raise PDFogOnlyLevel0Refusal(
            "fog-only level 0: a water-free vapour still exceeds the certified water "
            "criterion; no dew-line vapour exists - refuse, do not clamp"
        )
    while True:
        middle = 0.5 * (low + high)
        if middle in (low, high):
            return high
        if excess(middle) <= 0.0:
            high = middle
        else:
            low = middle


def _fog_only_level0_ledger(
    stream: PDBinaryGasStream,
    *,
    label: str,
    properties: object,
    datum_adapter: object,
    relative_limit: float,
) -> PDFogOnlyLevel0Ledger | None:
    """D12 fog-only level 0: the adiabatic flash of one booked stream into vapour and fog.

    ``None`` where the stream sits on or inside its certified water dew line (no fog),
    at or above water's critical temperature, or where a certified liquid reference
    does not exist (the guard then refuses with its own type and message).  An
    n-hexane excess refuses typed: a hexane fog is outside level 0.
    """

    temperature = stream.temperature_k
    pressure = stream.pressure_pa
    if temperature >= water_properties.TC:
        return None
    try:
        gas_state = binary_gas.state(
            temperature,
            pressure,
            stream.water_mole_fraction,
            stream.hexane_mole_fraction,
        )
        pure_water = azeotrope._pure_liquid_fugacity(water_properties, temperature, pressure)
        hexane_excess = math.nan
        if temperature < hexane.TC:
            hexane_excess = math.fsum(
                (
                    gas_state.fugacity_hexane_pa,
                    -azeotrope._pure_liquid_fugacity(hexane, temperature, pressure),
                )
            )
    except azeotrope.HeteroazeotropeError:
        return None
    if hexane_excess > 0.0:
        raise PDFogOnlyLevel0Refusal(
            f"{label} carries an n-hexane excess of {hexane_excess!r} Pa over its certified "
            "dew line; a hexane fog is outside fog-only level 0 - refuse, do not clamp"
        )
    water_excess = math.fsum((gas_state.fugacity_water_pa, -pure_water))
    if water_excess <= 0.0:
        return None
    booked_hexane = stream.hexane_mole_fraction
    booked_enthalpy = _bulk_gas_enthalpy_j_mol(
        properties,
        temperature_k=temperature,
        pressure_pa=pressure,
        hexane_mole_fraction=booked_hexane,
        water_mole_fraction=stream.water_mole_fraction,
        datum_adapter=datum_adapter,
    )

    def split(trial_temperature_k: float) -> tuple[float, ...]:
        pure = azeotrope._pure_liquid_fugacity(water_properties, trial_temperature_k, pressure)
        vapour_hexane = _fog_dew_line_hexane_fraction(
            trial_temperature_k, pressure, booked_hexane, pure
        )
        # Per mole of the booked stream: the n-hexane stays in the vapour.
        vapour = booked_hexane / vapour_hexane
        fog = math.fsum((1.0, -vapour))
        vapour_enthalpy = _bulk_gas_enthalpy_j_mol(
            properties,
            temperature_k=trial_temperature_k,
            pressure_pa=pressure,
            hexane_mole_fraction=vapour_hexane,
            water_mole_fraction=1.0 - vapour_hexane,
            datum_adapter=datum_adapter,
        )
        liquid_enthalpy = datum_adapter.common_datum_liquid_enthalpy_j_mol(
            "water", temperature_k=trial_temperature_k, pressure_pa=pressure
        )
        residual = math.fsum((vapour * vapour_enthalpy, fog * liquid_enthalpy, -booked_enthalpy))
        return residual, vapour_hexane, vapour, fog, vapour_enthalpy, liquid_enthalpy, pure

    low = temperature
    if split(low)[0] >= 0.0:
        raise PDFogOnlyLevel0Refusal(
            f"{label}: the fog flash has no root above the booked temperature "
            f"{temperature!r} K - refuse, do not clamp"
        )
    high = None
    for doubling in range(10):
        trial = temperature + 0.25 * 2.0**doubling
        if trial >= water_properties.TC:
            break
        if split(trial)[0] >= 0.0:
            high = trial
            break
        low = trial
    if high is None:
        raise PDFogOnlyLevel0Refusal(
            f"{label}: the fog flash root is not bracketed within 128 K of the booked "
            f"temperature {temperature!r} K - refuse, do not clamp"
        )
    while True:
        middle = 0.5 * (low + high)
        if middle in (low, high):
            break
        if split(middle)[0] < 0.0:
            low = middle
        else:
            high = middle
    residual, vapour_hexane, vapour, fog, vapour_enthalpy, liquid_enthalpy, pure = split(high)
    relative = abs(residual) / max(
        abs(booked_enthalpy), abs(vapour * vapour_enthalpy), abs(fog * liquid_enthalpy)
    )
    if not relative <= relative_limit:
        raise PDFogOnlyLevel0Refusal(
            f"{label}: the fog flash closes its enthalpy to {relative!r} relative, above the "
            f"step's limit {relative_limit!r} - refuse, do not clamp"
        )
    vapour_gas = binary_gas.state(high, pressure, 1.0 - vapour_hexane, vapour_hexane)
    flow = stream.molar_flow_mol_s
    fog_flow = fog * flow
    return PDFogOnlyLevel0Ledger(
        label=label,
        law_id=PD_FOG_ONLY_LEVEL0_LAW_ID,
        pressure_pa=pressure,
        booked_temperature_k=temperature,
        booked_molar_flow_mol_s=flow,
        booked_hexane_mole_fraction=booked_hexane,
        booked_water_mole_fraction=stream.water_mole_fraction,
        booked_enthalpy_j_mol=booked_enthalpy,
        water_excess_pa=water_excess,
        hexane_excess_pa=hexane_excess,
        flash_temperature_k=high,
        vapour_hexane_mole_fraction=vapour_hexane,
        vapour_water_mole_fraction=1.0 - vapour_hexane,
        vapour_water_excess_pa=math.fsum((vapour_gas.fugacity_water_pa, -pure)),
        vapour_molar_flow_mol_s=vapour * flow,
        fog_water_molar_flow_mol_s=fog_flow,
        fog_water_mass_flow_kg_s=fog_flow * properties.water_molar_mass_kg_mol,
        fog_liquid_enthalpy_j_mol=liquid_enthalpy,
        flash_energy_residual_j_mol=residual,
        flash_relative_residual=relative,
        relative_limit=relative_limit,
    )


def _fog_only_level0_admitted(
    stream: PDBinaryGasStream,
    *,
    label: str,
    properties: object,
    datum_adapter: object,
    relative_limit: float,
    ledgers: dict[str, object],
) -> tuple[tuple[str, float], ...]:
    """D12 fog-only level 0: the excess the guard admits at one site, and its fog ledger.

    The admitted bound is the water excess the fog ledger measured for this very
    stream, by the guard's own arithmetic, so ``_require_single_vapour_phase`` admits
    exactly this stream and refuses any stream more supersaturated than it.
    """

    ledger = _fog_only_level0_ledger(
        stream,
        label=label,
        properties=properties,
        datum_adapter=datum_adapter,
        relative_limit=relative_limit,
    )
    ledgers[label] = ledger
    if ledger is None:
        return ()
    return (("water", ledger.water_excess_pa),)


@dataclass(frozen=True, slots=True, kw_only=True)
class PDSharedWallParameters:
    """One physical PD wall; coefficients are aggregate tray declarations."""

    thermal_capacity_j_k: float
    floor_to_bottom_solid_ua_w_k: float
    wall_to_bypass_gas_ua_w_k: float
    ambient_ua_w_k: float
    ambient_temperature_k: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    authenticated_machine_parameters: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_float("wall thermal capacity", self.thermal_capacity_j_k, positive=True)
        _require_float("floor-to-solid UA", self.floor_to_bottom_solid_ua_w_k, nonnegative=True)
        _require_float("wall-to-bypass-gas UA", self.wall_to_bypass_gas_ua_w_k, nonnegative=True)
        _require_float("ambient UA", self.ambient_ua_w_k, nonnegative=True)
        _require_float("ambient temperature", self.ambient_temperature_k, positive=True)


def _validate_exact_pd_geometry(geometry: cc.CellGeometry) -> None:
    if type(geometry) is not cc.CellGeometry:
        raise PDKDomeBypassConfigurationError("geometry must be an exact CellGeometry")
    reference = next(
        (item for item in REFERENCE_TRAYS if item.tray_id == geometry.tray.tray_id),
        None,
    )
    if reference is None or reference.tray_id not in PD_REFERENCE_TRAY_IDS:
        raise PDKDomeBypassConfigurationError("the bounded integration accepts only PD1--PD3")
    kind = getattr(geometry, "tray" + "_type")
    if geometry.tray != reference:
        raise PDKDomeBypassConfigurationError(
            "PD geometry and typed family must equal the source-pinned reference tray exactly"
        )
    try:
        cc.require_type_matches_tray(geometry.tray, kind)
    except ValueError as error:
        raise PDKDomeBypassConfigurationError(
            "PD geometry and typed family must equal the source-pinned reference tray exactly"
        ) from error
    if (
        kind.family is not cc.TrayFloorFamily.UNPERFORATED_HEATED_DISC
        or kind.gas_source.value != "none"
        or kind.pressure_loss_relation is not cc.PressureLossRelation.DOME_BYPASS
        or kind.contact is not VaporContact.DOME_CONTACT
        or not kind.indirect_heating
        or kind.free_area_fraction != 0.0
        or kind.characteristic_opening_m != 0.0
        or kind.opening_count_density_per_m2 is not None
    ):
        raise PDKDomeBypassConfigurationError(
            "PD1--PD3 require exact family-1 zero-bore DOME_BYPASS/DOME_CONTACT topology"
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDKTrayEngineeringModel:
    """Immutable declarations for one exact PD reference tray."""

    geometry: cc.CellGeometry
    wall: PDSharedWallParameters
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties
    #: WS-A: the DECLARED contact-efficiency multiplier on the existing
    #: DOME_CONTACT top free-surface area (design packet D3).  0.0 - the
    #: default and the committed production value per Q-PD-4 - reproduces
    #: the pre-WS-A behaviour bit for bit, including the model digest.
    condensation_contact_fraction: float = PD_CONDENSATION_CONTACT_FRACTION_DEFAULT

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    pressure_loss_arithmetic_implemented: ClassVar[bool] = False
    ergun_or_bed_pressure_loss_implemented: ClassVar[bool] = False
    direct_steam_implemented: ClassVar[bool] = False
    external_water_transport_implemented: ClassVar[bool] = False
    high_loading_activation_handoff_implemented: ClassVar[bool] = False
    industrial_operating_envelope_qualified: ClassVar[bool] = False
    finite_rpm_k_mixing_implemented: ClassVar[bool] = True
    inter_k_solid_heat_transfer_implemented: ClassVar[bool] = True
    inter_k_effective_conductivity_w_m_k: ClassVar[float] = PD_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K
    inter_k_conductivity_authority_id: ClassVar[str] = PD_INTER_K_CONDUCTIVITY_AUTHORITY_ID
    inter_k_allocation_law_id: ClassVar[str] = PD_INTER_K_ALLOCATION_LAW_ID
    shaft_work_energy_model_selected: ClassVar[bool] = False
    #: Q-PD-4: this flag is the CAPABILITY GATE, not a decoration.  While it
    #: is ``False`` - the committed production value, flipped only at the
    #: consolidated re-baseline - a nonzero ``condensation_contact_fraction``
    #: is REFUSED at construction, so no run can condense water while
    #: stamping ``flash_or_condensation_implemented: False`` onto its record.
    flash_or_condensation_implemented: ClassVar[bool] = False
    #: WS-A declarations.  The capability above stays FALSE (Q-PD-4); these
    #: name the law and the basis so a reader of a step record never has to
    #: infer either.
    condensation_contact_area_basis_id: ClassVar[str] = PD_CONDENSATION_CONTACT_AREA_BASIS_ID
    condensation_rate_law_id: ClassVar[str] = PD_CONDENSATION_RATE_LAW_ID
    #: CS-1a declarations, additive.
    condensation_cap_law_id: ClassVar[str] = PD_CONDENSATION_CAP_LAW_ID
    condensation_transfer_coefficient_law_id: ClassVar[str] = (
        PD_CONDENSATION_TRANSFER_COEFFICIENT_LAW_ID
    )
    condensation_equilibrium_ratio_law_id: ClassVar[str] = (
        PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID
    )
    condensation_per_component_driving_force_implemented: ClassVar[bool] = True
    condensation_per_component_lewis_correction_implemented: ClassVar[bool] = False
    condensation_pore_clearing_law_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _validate_exact_pd_geometry(self.geometry)
        _require_float(
            "declared condensation contact fraction",
            self.condensation_contact_fraction,
            nonnegative=True,
        )
        if self.condensation_contact_fraction > 1.0:
            raise PDKDomeBypassConfigurationError(
                "the declared condensation contact fraction is an efficiency on the "
                "existing DOME_CONTACT area and lies in [0, 1]"
            )
        if self.condensation_contact_fraction != 0.0 and (
            not self.flash_or_condensation_implemented
        ):
            # Q-PD-4, MADE EXECUTABLE 2026-09-04 (verification pass): before
            # this gate the claim flag gated NOTHING and an armed run stamped
            # ``flash_or_condensation_implemented: False`` onto a record in
            # which condensation had fired.  The fence is now the flag itself.
            raise PDKDomeBypassConfigurationError(
                "arming the contact-limited dome condensation capability requires the "
                "declared flash_or_condensation_implemented capability flag, which "
                "owner ruling Q-PD-4 reserves for the consolidated re-baseline"
            )
        if type(self.wall) is not PDSharedWallParameters:
            raise PDKDomeBypassConfigurationError("wall must be exact PDSharedWallParameters")
        # CELL-02e-pd: the same closed caloric-mode union and W2 forgery gate
        # as the K-layer models, landing atomically with the native arms.
        if type(self.properties) is not cc.DeclaredGasProperties and (
            type(self.properties) is not native_caloric.NativeGasProperties
        ):
            raise PDKDomeBypassConfigurationError(
                "properties must be an exact DeclaredGasProperties or NativeGasProperties"
            )
        if type(self.properties) is cc.DeclaredGasProperties and (
            self.properties.energy_datum_id == native_caloric.NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID
        ):
            raise PDKDomeBypassConfigurationError(
                "the native zero-gauge datum identity cannot be carried by a "
                "declared constant-property closure"
            )

    @property
    def definition_digest(self) -> str:
        parts = (
            repr(self.geometry).encode("utf-8"),
            repr(self.wall).encode("utf-8"),
            repr(self.properties).encode("utf-8"),
            repr(
                (
                    PD_INTER_K_MODEL_SCHEMA_ID,
                    PD_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K,
                    PD_INTER_K_CONDUCTIVITY_AUTHORITY_ID,
                    PD_INTER_K_ALLOCATION_LAW_ID,
                )
            ).encode("utf-8"),
        )
        if self.condensation_contact_fraction != 0.0:
            # WS-A, ADDITIVE ONLY: an armed condensation capability is a
            # different model and must digest differently.  At the declared
            # default (0.0) no component is appended at all, so every
            # pre-WS-A digest is reproduced bit for bit.
            parts = parts + (
                repr(
                    (
                        PD_CONDENSATION_MODEL_SCHEMA_ID,
                        PD_CONDENSATION_RATE_LAW_ID,
                        PD_CONDENSATION_CONTACT_AREA_BASIS_ID,
                        PD_CONDENSATION_ENERGY_ROUTING_ID,
                        PD_CONDENSATION_RECEIVER_LAW_ID,
                        # CS-1a, additive INSIDE the armed branch only: the cap,
                        # the transfer-coefficient and the surface-ratio laws
                        # are part of the model an armed run declares.
                        PD_CONDENSATION_CAP_LAW_ID,
                        PD_CONDENSATION_TRANSFER_COEFFICIENT_LAW_ID,
                        PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID,
                        # RULED 2026-09-06 (the 465 K ceiling ruling, point 8
                        # of the adversarial pass): the two stop-temperature
                        # brackets decide WHICH states proceed and which
                        # refuse, so they belong to the armed model's
                        # identity.  Inside the armed branch only: at the
                        # declared default no component is appended at all and
                        # every unarmed digest stays bit-identical.
                        PD_CONDENSATION_STOP_BRACKET_HEXANE_K,
                        PD_CONDENSATION_STOP_BRACKET_WATER_K,
                        self.condensation_contact_fraction,
                    )
                ).encode("utf-8"),
            )
        return _framed_digest(INTEGRATION_DIGEST_DOMAIN + "/model", parts)


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedPDKTrayState:
    """Accepted K/RTD packet host plus the one shared wall and plenum outlet."""

    tray: tray_host.AcceptedKCellTrayState
    model_definition_digest: str
    wall_temperature_k: float
    plenum_outlet: PDBinaryGasStream

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    authenticated_pressure_network: ClassVar[bool] = False
    industrial_operating_envelope_qualified: ClassVar[bool] = False
    finite_rpm_k_mixing_implemented: ClassVar[bool] = True
    inter_k_solid_heat_transfer_implemented: ClassVar[bool] = True
    external_water_inventory_thermally_coupled: ClassVar[bool] = True
    external_water_transport_implemented: ClassVar[bool] = False
    shaft_work_energy_model_selected: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.tray) is not tray_host.AcceptedKCellTrayState:
            raise PDKDomeBypassConfigurationError("tray must be an accepted K-cell host state")
        if self.tray.physical_tray_id not in PD_REFERENCE_TRAY_IDS:
            raise PDKDomeBypassConfigurationError("accepted PD state must belong to PD1--PD3")
        _require_sha256("accepted model definition", self.model_definition_digest)
        _require_float("accepted shared-wall temperature", self.wall_temperature_k, positive=True)
        if type(self.plenum_outlet) is not PDBinaryGasStream:
            raise PDKDomeBypassConfigurationError("plenum outlet has a foreign type")
        if not _same_float(self.plenum_outlet.pressure_pa, self.tray.pressure_pa):
            raise PDKDomeBypassConfigurationError(
                "accepted plenum and fixed packet/tray pressure must be identical"
            )
        occupied = {packet.owner_key.vertical_layer_id.value for packet in self.tray.packets}
        expected = {cell.vertical_layer_id.value for cell in self.tray.cells}
        if occupied != expected:
            raise PDKDomeBypassConfigurationError(
                "every declared PD K layer must contain at least one packet"
            )

    @property
    def state_digest(self) -> str:
        return _framed_digest(
            INTEGRATION_DIGEST_DOMAIN + "/state",
            (
                self.tray.state_digest.encode("ascii"),
                self.model_definition_digest.encode("ascii"),
                _float_bytes(self.wall_temperature_k),
                repr(self.plenum_outlet).encode("utf-8"),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDKTrayStepRequest:
    """Dynamic boundaries for one fixed-pressure PD K/RTD macro-step."""

    attempt_id: str
    prior: AcceptedPDKTrayState
    end_time_s: float
    model: PDKTrayEngineeringModel
    bypass_inlet: PDBinaryGasStream
    prescribed_bypass_pressure_drop_pa: float
    accepted_indirect_duty_w: float
    flow_segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...]
    relative_limit: float = DEFAULT_RELATIVE_LIMIT
    meal_mechanical_work_deposition: meal_work.QSCMealMechanicalWorkDepositionBundle | None = None
    #: B7-family memo (2026-08-31 rulings): the interval-preverified bundle
    #: digest; None keeps the sealed full-reconstruction behavior.
    preverified_meal_bundle_digest: str | None = None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    pressure_drop_is_prescribed_not_predicted: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise PDKDomeBypassConfigurationError("attempt_id must be a nonblank exact str")
        if type(self.prior) is not AcceptedPDKTrayState:
            raise PDKDomeBypassConfigurationError("prior has the wrong PD accepted-state type")
        _require_float("step end time", self.end_time_s)
        if self.end_time_s <= self.prior.tray.time_s:
            raise PDKDomeBypassConfigurationError("the PD macro-step must advance time")
        if type(self.model) is not PDKTrayEngineeringModel:
            raise PDKDomeBypassConfigurationError("model has a foreign type")
        if self.model.definition_digest != self.prior.model_definition_digest:
            raise PDKDomeBypassConfigurationError("PD model drifted from the accepted state")
        if self.model.geometry.tray.tray_id != self.prior.tray.physical_tray_id:
            raise PDKDomeBypassConfigurationError("PD model and packet host name different trays")
        if type(self.bypass_inlet) is not PDBinaryGasStream:
            raise PDKDomeBypassConfigurationError("bypass inlet has a foreign type")
        _require_float(
            "prescribed bypass pressure drop",
            self.prescribed_bypass_pressure_drop_pa,
            nonnegative=True,
        )
        _require_float("accepted indirect duty", self.accepted_indirect_duty_w, nonnegative=True)
        _require_float("PD relative limit", self.relative_limit, positive=True)
        expected_plenum_pressure = math.fsum(
            (self.bypass_inlet.pressure_pa, -self.prescribed_bypass_pressure_drop_pa)
        )
        if not _same_float(expected_plenum_pressure, self.prior.tray.pressure_pa):
            raise PDKDomeBypassConfigurationError(
                "external bypass algebra requires inlet pressure minus prescribed drop "
                "to equal the fixed tray/plenum pressure exactly"
            )
        if type(self.flow_segments) is not tuple or not self.flow_segments:
            raise PDKDomeBypassConfigurationError("a nonempty accepted-flow history is required")
        if any(
            type(item) is not tray_host.AcceptedDryMatterFlowSegment for item in self.flow_segments
        ):
            raise PDKDomeBypassConfigurationError("flow history contains a foreign segment")
        if self.meal_mechanical_work_deposition is not None:
            bundle = self.meal_mechanical_work_deposition
            if type(bundle) is not meal_work.QSCMealMechanicalWorkDepositionBundle:
                raise PDKDomeBypassConfigurationError(
                    "meal mechanical-work deposition has a foreign exact type"
                )
            meal_work.QSCMealMechanicalWorkDepositionBundle.__post_init__(bundle)
            if (
                bundle.prior_tray_state_digest != self.prior.tray.state_digest
                or bundle.physical_tray_id != self.prior.tray.physical_tray_id
                or bundle.t_n_s != self.prior.tray.time_s
                or bundle.end_time_s != self.end_time_s
                or bundle.energy_datum_id != self.prior.tray.energy_datum_id
            ):
                raise PDKDomeBypassConfigurationError(
                    "meal mechanical-work deposition differs from the local PD t_n tray step"
                )

    @property
    def duration_s(self) -> float:
        return self.end_time_s - self.prior.tray.time_s


def _declared_bulk_gas_enthalpy_j_mol(
    properties: cc.DeclaredGasProperties,
    *,
    temperature_k: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
    datum_adapter: component_datum.NumericComponentEnergyDatumAdapter,
) -> float:
    datum_adapter.require_gas_properties(properties)
    hexane_shift, water_shift = datum_adapter.common_datum_species_shifts_j_mol
    return math.fsum(
        (
            properties.molar_heat_capacity_j_mol_k
            * (temperature_k - properties.energy_datum_temperature_k),
            hexane_mole_fraction * (properties.hexane_latent_heat_j_mol + hexane_shift),
            water_mole_fraction * (properties.water_latent_heat_j_mol + water_shift),
        )
    )


def _bulk_gas_enthalpy_j_mol(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
    datum_adapter: object,
) -> float:
    """CELL-02e-pd stream-enthalpy dispatcher across the closed caloric union.

    The declared arm is the original common-datum arithmetic, verbatim, and
    ignores the pressure (its enthalpy is pressure-free by model form); the
    native arm evaluates the frozen mixture law at the stream's own booked
    pressure and, per the W2 proof, never consults the declared adapter.
    """

    if type(properties) is cc.DeclaredGasProperties:
        return _declared_bulk_gas_enthalpy_j_mol(
            properties,
            temperature_k=temperature_k,
            hexane_mole_fraction=hexane_mole_fraction,
            water_mole_fraction=water_mole_fraction,
            datum_adapter=datum_adapter,
        )
    if water_mole_fraction != 1.0 - hexane_mole_fraction:
        raise PDKDomeBypassConfigurationError(
            "the native binary stream requires the exact complementary water fraction"
        )
    return native_caloric.mixture_enthalpy(
        temperature_k, pressure_pa, hexane_mole_fraction
    ).value_j_mol


def _pure_hexane_vapour_enthalpy_j_mol(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    temperature_k: float,
    pressure_pa: float,
    datum_adapter: object,
) -> float:
    return _bulk_gas_enthalpy_j_mol(
        properties,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        hexane_mole_fraction=1.0,
        water_mole_fraction=0.0,
        datum_adapter=datum_adapter,
    )


def _hexane_saturation_temperature_k(
    pressure_pa: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    lo = codec.sphere_params.wet_core.T_min
    hi = codec.sphere_params.wet_core.T_max
    return _cached_hexane_saturation_temperature_k(pressure_pa, lo, hi)


@dataclass(frozen=True, slots=True, kw_only=True)
class PDCoBoilingState:
    """The certified two-liquid co-boiling state this tray boils at.

    F4 (2026-09-02, wave-2 sharpening of wave-1 item 4).  Until this fix the
    PD branch boiled at the PURE n-hexane saturation temperature and booked
    100 % of the departing vapour as n-hexane.  Both are wrong whenever the
    packet still carries external free water, which is the whole admitted PD
    branch: with two immiscible liquids present the bed boils at the
    heteroazeotrope and the vapour carries the certified water share.

    THE LAW, corrected 2026-09-16 (D12-a's Class-A doc-slip pass; the defect
    is finding S:4.4 of
    ``docs/GT_PS2_PD_GUARD_VS_COBOILING_AUTHORITY_DIAGNOSIS_2026-09-04.md``,
    and the discharge of I-8 of the Stefan-drift record).  This docstring
    used to state the law as the IDEAL closure ``Psat_h(T) + Psat_w(T) = P``
    while the implementation it documents has always solved the FUGACITY
    form - the prose/implementation mismatch that is plausibly how the two
    conventions came to sit side by side unnoticed.  What the certified
    authority actually solves is the pair of per-component equalities

        y_h phi_h(T, P, y) P = a_h f_h^{L,pure}(T, P)
        y_w phi_w(T, P, y) P =     f_w^{L,pure}(T, P)

    with ``y_h + y_w = 1`` (``heteroazeotrope._y_hexane_from_hexane_equality``
    and ``_water_residual``); the ideal closure is what those two become when
    every ``phi`` is set to one and each ``f^{L,pure}`` to its ``Psat``.  The
    two roots differ: at 101 325 Pa the ideal one is 334.5475746260082 K
    against the certified 334.4803981593741 K quoted below, 0.067 K apart.

    Measured at 101 325 Pa: the pure-hexane root
    is 341.8645125361387 K against the certified T_het 334.4803981593741 K -
    **7.3841143767646 K too high** - and the certified vapour is
    ``y_water = 0.20516853343923214``, not 0.

    Both numbers come from the certified authority
    (:mod:`dtdc_simulator.core2.heteroazeotrope`, the G5-grounded solver
    ruled 2026-08-22) and NOTHING is copied here.  Outside the G5-grounded
    pressure band the authority REFUSES and so does this branch: no clamp,
    no fallback to the pure-hexane root.
    """

    temperature_k: float
    vapour_hexane_mole_fraction: float
    vapour_water_mole_fraction: float
    pressure_pa: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def water_moles_per_hexane_mole(self) -> float:
        return self.vapour_water_mole_fraction / self.vapour_hexane_mole_fraction


@functools.lru_cache(maxsize=32)
def _pd_co_boiling_state(pressure_pa: float) -> PDCoBoilingState:
    """The certified heteroazeotrope at the fixed tray pressure."""

    try:
        temperature_k = azeotrope.grounded_heteroazeotrope_temperature_k(pressure_pa)
        certified = azeotrope.solve_heteroazeotrope(pressure_pa)
    except azeotrope.HeteroazeotropeError as error:
        raise PDKDomeBypassConfigurationError(
            f"the certified heteroazeotrope authority refused the PD tray pressure: {error}"
        ) from error
    if certified.temperature_k != temperature_k:
        raise PDKDomeBypassConfigurationError(
            "the G5 gate temperature and the certified solver disagree at the PD tray pressure"
        )
    y_hexane = certified.y_hexane
    if not (0.0 < y_hexane < 1.0):
        raise PDKDomeBypassConfigurationError(
            f"the certified co-boiling vapour composition {y_hexane!r} is outside (0, 1)"
        )
    return PDCoBoilingState(
        temperature_k=temperature_k,
        vapour_hexane_mole_fraction=y_hexane,
        vapour_water_mole_fraction=1.0 - y_hexane,
        pressure_pa=pressure_pa,
    )


def _co_boiling_vapour_enthalpy_j_mol(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    co_boiling: PDCoBoilingState,
    pressure_pa: float,
    datum_adapter: object,
) -> float:
    """Molar enthalpy of ONE mole of the certified co-boiling vapour."""

    return _bulk_gas_enthalpy_j_mol(
        properties,
        temperature_k=co_boiling.temperature_k,
        pressure_pa=pressure_pa,
        hexane_mole_fraction=co_boiling.vapour_hexane_mole_fraction,
        water_mole_fraction=co_boiling.vapour_water_mole_fraction,
        datum_adapter=datum_adapter,
    )


def _vapour_enthalpy_per_hexane_mole_j(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    co_boiling: PDCoBoilingState | None,
    saturation_temperature_k: float,
    pressure_pa: float,
    datum_adapter: object,
) -> float:
    """Common-datum vapour enthalpy carried per mole of departing n-hexane.

    On the certified co-boiling branch a mole of hexane leaves inside a mole
    of MIXTURE only ``y_h`` of which is hexane, so the energy leaving with it
    is ``h_mixture / y_h``.  Without a co-boiling state (single-liquid packet)
    this is the pre-F4 pure-hexane enthalpy, bit-identically.
    """

    if co_boiling is None:
        return _pure_hexane_vapour_enthalpy_j_mol(
            properties,
            temperature_k=saturation_temperature_k,
            pressure_pa=pressure_pa,
            datum_adapter=datum_adapter,
        )
    return (
        _co_boiling_vapour_enthalpy_j_mol(
            properties,
            co_boiling=co_boiling,
            pressure_pa=pressure_pa,
            datum_adapter=datum_adapter,
        )
        / co_boiling.vapour_hexane_mole_fraction
    )


def _water_mass_per_hexane_mass(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    co_boiling: PDCoBoilingState | None,
) -> float:
    """kg of external free water that co-boils per kg of n-hexane evaporated.

    ``(y_w / y_h) * (M_w / M_h)`` from the certified vapour composition; 0.0
    on the single-liquid branch.
    """

    if co_boiling is None:
        return 0.0
    return (
        co_boiling.water_moles_per_hexane_mole
        * properties.water_molar_mass_kg_mol
        / properties.hexane_molar_mass_kg_mol
    )


def _level0_source_bound_co_boiling_state(
    certified: PDCoBoilingState,
    water_share_factor: float,
) -> PDCoBoilingState:
    """D12 level 0 SOURCE BOUND: the certified vapour with part of its water withheld.

    Same temperature (the certified lock) and pressure; the water moles per
    hexane mole are the certified ratio times ``water_share_factor`` in
    (0, 1], and the withheld water stays on the meal as external free water.
    The factor is the knob module's fixed point on the candidate plenum
    outlet's certified water excess, refused there above 1 percent (design
    ``GT_PS2_D12_STATE_AND_LEVEL0_DESIGN_2026-09-28.md`` section 3.2, SOURCE
    BOUND).  Reached only under ``PD_FLASH_LEVEL0_ARMED``.
    """

    if type(water_share_factor) is not float or not (0.0 < water_share_factor <= 1.0):
        raise PDKDomeBypassConfigurationError(
            f"the level-0 co-boiled water share factor {water_share_factor!r} is outside (0, 1]"
        )
    ratio = certified.water_moles_per_hexane_mole * water_share_factor
    y_hexane = 1.0 / math.fsum((1.0, ratio))
    return PDCoBoilingState(
        temperature_k=certified.temperature_k,
        vapour_hexane_mole_fraction=y_hexane,
        vapour_water_mole_fraction=1.0 - y_hexane,
        pressure_pa=certified.pressure_pa,
    )


@functools.lru_cache(maxsize=32)
def _cached_hexane_saturation_temperature_k(
    pressure_pa: float,
    lo: float,
    hi: float,
) -> float:
    residual_lo = hexane.saturation_pressure(lo) - pressure_pa
    residual_hi = hexane.saturation_pressure(hi) - pressure_pa
    if residual_lo > 0.0 or residual_hi < 0.0:
        raise PDKDomeBypassConfigurationError(
            "fixed tray pressure has no pure-hexane saturation root inside the tag-1 bracket"
        )
    if residual_lo == 0.0:
        return lo
    if residual_hi == 0.0:
        return hi
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        residual_mid = hexane.saturation_pressure(mid) - pressure_pa
        if residual_mid < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass(frozen=True, slots=True, kw_only=True)
class PDContactCondensationLedger:
    """CS-1a: the cold-surface condensation this step took, PER COMPONENT.

    Every field is a MEASUREMENT of the step that produced it: the two bulk
    vapour fugacities, the two surface liquid fugacities, the two surface
    equilibrium ratios, the two driving forces, the two stop temperatures,
    the contact-limited capacities, the first-stop energy cap and which
    component set it, and the seam split of the energy each condensate
    carried out of the gas and into the meal.

    The GAS-SIDE numbers ``water_partial_pressure_pa``,
    ``water_saturation_pressure_pa``, ``partial_pressure_excess_pa``,
    ``equilibrium_relief_mol`` and ``residual_supersaturation_pa`` are still
    measured, exactly as WS-A measured them, because the fog guard
    ``_require_single_vapour_phase`` reads two of them.  They are NOT the
    condensation criterion and ``equilibrium_relief_mol`` is NOT a cap any
    more: an unsaturated gas has no relief at all, and capping a cold-surface
    condensate at the gas's own saturation would silence the term at exactly
    the state this workstream exists for.

    A stop temperature of a component that is NOT condensing can lie outside
    the bisection bracket; it is then recorded as ``nan`` and is never read.

    **REV 2 (2026-09-04, finding F-CS1-2).**  Every SCALAR field naming a
    surface quantity - the surface temperature, the two surface fugacities,
    the two equilibrium ratios, the two driving forces and the two stop
    temperatures - is a representative-weighted REPORT over ``packet_rows``,
    not a quantity any decision was taken on.  The decisions are per packet.
    ``equilibrium_energy_scaling_factor`` is likewise the realized share of
    the layer's uncapped latent energy, so ``equilibrium_limited`` still says
    "some packet's own surface cap bound the answer".
    """

    contact_fraction: float
    contact_area_m2: float
    contact_area_basis_id: str
    rate_law_id: str
    energy_routing_id: str
    receiver_law_id: str
    cap_law_id: str
    transfer_coefficient_law_id: str
    equilibrium_ratio_law_id: str
    equilibrium_ratio_basis_id: str
    duration_s: float
    stream_temperature_k: float
    stream_pressure_pa: float
    #: The COLD SOLID the vapour meets, REPORTED as the top layer's
    #: representative-weighted mean surface temperature.  REV 2: this is a
    #: report, NOT the criterion.  The driving force is evaluated on each
    #: packet's own temperature in ``packet_rows`` and summed (F-CS1-2); on a
    #: layer whose packets all carry one temperature the two coincide exactly.
    surface_temperature_k: float
    #: CS-1a rev 2: the per-packet evaluation this layer ledger is the sum of.
    #: One row per packet of the receiving top layer, in the layer's own
    #: packet order.  Every scalar aggregate above and below is derived from
    #: these rows; none of them is an input to any decision.
    packet_rows: tuple[PDPacketCondensationRow, ...]
    gas_film_coefficient_w_m2_k: float
    molar_heat_capacity_j_mol_k: float
    molar_transfer_conductance_mol_s_m2_pa: float
    #: The vapour fugacities y_i phi_i P of the arriving bulk stream.
    hexane_bulk_fugacity_pa: float
    water_bulk_fugacity_pa: float
    #: The pure compressed-liquid fugacities f_i^{L,pure}(T_s, P).
    hexane_surface_liquid_fugacity_pa: float
    water_surface_liquid_fugacity_pa: float
    hexane_equilibrium_ratio: float
    water_equilibrium_ratio: float
    hexane_driving_force_pa: float
    water_driving_force_pa: float
    #: The temperature at which each component's condensation STOPS: the root
    #: of ``f_eq,i f_i^{L,pure}(T, P) = y_i phi_i P`` at the surface ratio.
    hexane_dew_point_k: float
    water_dew_point_k: float
    #: ``min`` of the stop temperatures over the CONDENSING components, and
    #: the component that set it.  This is the endpoint a condensing packet
    #: may relax onto.
    condensation_equilibrium_temperature_k: float
    stop_component: str
    hexane_contact_limited_capacity_mol: float
    contact_limited_capacity_mol: float
    #: The explicit declared energy this step ALREADY allocates to the
    #: receiving layer -- inter-K face heat PLUS declared mechanical (shaft)
    #: work -- which the first-stop cap subtracts before it prices latent.
    upper_layer_declared_energy_j: float
    equilibrium_energy_cap_j: float
    equilibrium_energy_scaling_factor: float
    hexane_condensed_mol: float
    hexane_condensed_kg: float
    condensed_water_mol: float
    condensed_water_kg: float
    hexane_common_datum_liquid_enthalpy_j_mol: float
    hexane_seam_latent_j_mol: float
    common_datum_liquid_enthalpy_j_mol: float
    seam_latent_j_mol: float
    hexane_condensation_common_energy_j: float
    water_condensation_common_energy_j: float
    #: The SUM of the two component energies - what the plenum gives up and
    #: what the top layer receives across the F2 seam.
    condensation_common_energy_j: float
    #: Gas-side measurements, kept for the fog guard.  Not the criterion.
    water_partial_pressure_pa: float
    water_saturation_pressure_pa: float
    partial_pressure_excess_pa: float
    equilibrium_relief_mol: float
    residual_supersaturation_pa: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    pore_clearing_coupling_implemented: ClassVar[bool] = False
    per_component_driving_force: ClassVar[bool] = True
    gas_saturation_cap_used: ClassVar[bool] = False

    @property
    def equilibrium_limited(self) -> bool:
        """True when the FIRST-STOP surface cap bound the answer.

        CS-1a: the cap is the surface reaching the stop temperature of the
        first component to stop, not the gas returning to saturation.
        """

        return self.equilibrium_energy_scaling_factor < 1.0

    @property
    def condensed_hexane_kg(self) -> float:
        """Alias of ``hexane_condensed_kg`` on the water field's naming."""

        return self.hexane_condensed_kg


@dataclass(frozen=True, slots=True, kw_only=True)
class PDMoisturePickupLedger:
    """WS-A / Q-PD-3: the moisture-pickup ledger.  REPORT-ONLY, never a bound.

    The owner's testimony names moisture-gain minimization as a DESIGN
    OBJECTIVE the model must expose.  An absolute condensate mass cannot be
    constrained without a denominator, so the operative quantity is the
    condensate as a fraction of the meal's own DRY mass, alongside the meal's
    wet-basis moisture before and after, checked - REPORTED, not enforced -
    against literature criterion A1 (10-12 % wb, flat across PD).
    """

    physical_tray_id: str
    contact_fraction: float
    condensed_water_kg: float
    dry_matter_kg: float
    condensate_fraction_of_dry_matter: float
    total_water_before_kg: float
    total_water_after_kg: float
    meal_moisture_before_wet_basis: float
    meal_moisture_after_wet_basis: float
    literature_band_wet_basis: tuple[float, float]
    literature_band_source: str
    moisture_before_inside_literature_band: bool
    moisture_after_inside_literature_band: bool
    #: CS-1a, ADDITIVE: the SOLVENT-pickup rows.  Hexane condensing on the
    #: cold meal raises the meal's own hexane loading exactly as the water
    #: condensate raises its moisture; the two are reported side by side so a
    #: reader never has to infer one from the other.  REPORT-ONLY, like every
    #: row of this ledger (Q-PD-3).
    condensed_hexane_kg: float = 0.0
    condensate_hexane_fraction_of_dry_matter: float = 0.0
    total_hexane_before_kg: float = 0.0
    total_hexane_after_kg: float = 0.0
    meal_hexane_loading_before_dry_basis: float = 0.0
    meal_hexane_loading_after_dry_basis: float = 0.0

    report_only: ClassVar[bool] = True
    enforced_upper_bound_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def _pd_dome_contact_area_m2(geometry: cc.CellGeometry) -> float:
    """The DOME_CONTACT area, from the one authority that already derives it."""

    return exchange_area_m2(geometry.tray)


def _stream_molar_heat_capacity_j_mol_k(
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    *,
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
) -> float:
    """The stream's own molar heat capacity, on whichever arm is in force.

    The declared arm carries it as a constant; the native arm has no such
    constant, so it is the ANALYTIC dT tangent of the frozen mixture enthalpy
    at the stream's own state.  Neither arm borrows the other's number.
    """

    if type(properties) is cc.DeclaredGasProperties:
        return properties.molar_heat_capacity_j_mol_k
    return native_caloric.mixture_enthalpy(
        temperature_k, pressure_pa, hexane_mole_fraction
    ).d_dT_j_mol_k


#: The frozen parameter sets the packet-own sorption ratios are read on.  They
#: are the module defaults of the certified sorption authority -- the same
#: no-argument construction the water arm has always used for
#: ``LuikovParams`` -- named here so a reader sees which closure was chosen.
#:
#: REV 3 (2026-09-06): the hexane arm no longer passes the two hexane objects
#: itself; it reads the regime map, whose
#: ``sorption_interface.ParticleSorptionParams`` carries these SAME two frozen
#: defaults (``gab=GabParams()``, ``oil=OilIsotherm()``).  They stay named here
#: because the parameter set is still the one the packet-own hexane ratio is
#: read on, and a test pins that identity so the two cannot drift apart.
PD_CONDENSATION_GAB_PARAMS = sorption.GabParams()
PD_CONDENSATION_OIL_ISOTHERM = sorption.OilIsotherm()
PD_CONDENSATION_LUIKOV_PARAMS = sorption.LuikovParams()


def _packet_surface_equilibrium_ratios(
    decoded: high_loading.CompleteHighLoadingPacket,
    *,
    temperature_k: float,
    pressure_pa: float,
) -> tuple[float, float, str]:
    """CS-1a rev 2: ``(f_eq_hexane, f_eq_water, basis_id)`` of ONE packet.

    ``f_eq,i`` is a property of the SURFACE this packet presents, so it is read
    from THIS packet's own decoded state -- never from a layer aggregate
    (priority brief S:6, finding F-CS1-2).

    * A free film of a component gives ``f_eq = 1`` exactly: the FILM_ACTIVE
      regime of ``sorption_interface.surface_state``, "bulk liquid touches the
      surface: a_h = 1 exactly".
    * Dry of free water, the packet's own modified-Luikov retained-water
      activity, on the DISCLOSED fugacity-ratio convention bias
      ``heteroazeotrope.activity_convention_ratio`` (the same convention the
      certified two-liquid solver applies ``a_w`` on).
    * Dry of attached hexane, the packet's own REGIME-MAP surface activity:
      the one shared helper ``film_free_surface_hexane_activity`` at this
      packet's CONSERVED ENVELOPE hexane loading, its oil weight fraction and
      its own temperature.  Unity while a wet core persists (RECEDING_FRONT:
      liquid stands at the front, so the equilibrium surface activity is one
      and all the hindrance is transport); the authority's storage-based
      sorbate activity, capped at the pressure-feasible ceiling
      ``sorption_interface.hexane_activity_ceiling(T, P)``, once the core is
      gone (DRY_SORBATE).  Nothing is assumed: an oil weight fraction outside
      the PHY-048 composition envelope, and a loading with no positive hexane
      at all (EXHAUSTED), both still REFUSE typed.

    REV 3, 2026-09-06 (owner ruling of that night on the peer-review audit,
    section 3.3 conditions i and ii; record
    ``GT_PS2_Q_CS1_2_SHARED_ACTIVITY_HELPER_RECORD_2026-09-06.md``).  Rev 2
    handed the packet's WHOLE internal inventory to the GAB SORBATE inverse
    ``props.sorption.activity_from_retained``.  Two measured defects:

    * a category error.  Above X_e that inventory contains pore liquid the
      sorbate isotherm does not represent, so the inverse refused the step on
      the entire receding-front band - 100 % of regime B, the rule and not an
      edge (the ceremony's D-13 recorded it as an edge because no ceremony run
      ever reached the branch), and it double-counted the pore-gas share of
      the storage where it did return a number;
    * an untyped crash.  It read the NATIVE state field
      ``state.internal_hexane_mass_kg``, which the tag-3 dry-sorbate state
      does not carry (it is ``sorbed_hexane_mass_kg``), so a tag-3 top-layer
      packet raised ``AttributeError`` past the typed refusal clause.  The
      loading is now read from the ENVELOPE inventory, in the kernel's own
      pattern ``inventory.internal_hexane_kg / inventory.dry_matter_kg``,
      which exists on tag-1, tag-2 and tag-3 alike.
    """

    if decoded.surface.attached_hexane_kg > 0.0:
        hexane_ratio = 1.0
        hexane_basis = "free-hexane-film"
    else:
        oil_fraction = decoded.state.oil_fraction
        if not sorption.W_O_MIN <= oil_fraction <= sorption.W_O_MAX:
            raise PDKDomeBypassStepError(
                "cold-surface condensation met a dry-of-hexane top-layer packet whose "
                f"oil weight fraction {oil_fraction!r} lies outside the PHY-048 "
                f"composition envelope [{sorption.W_O_MIN}, {sorption.W_O_MAX}]; the "
                "sorption authority has no surface state there; refuse, do not clamp"
            )
        # The CONSERVED ENVELOPE loading, never a native state field: the
        # three packet variants name their native hexane differently and only
        # the envelope is common to all of them.
        internal_hexane_loading = (
            decoded.inventory.internal_hexane_kg / decoded.inventory.dry_matter_kg
        )
        try:
            binding = film_free_activity.film_free_surface_hexane_activity(
                hexane_loading_kg_kg_dry=internal_hexane_loading,
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                residual_oil_mass_fraction=oil_fraction,
            )
        except film_free_activity.FilmFreeHexaneActivityError as error:
            raise PDKDomeBypassStepError(
                "cold-surface condensation could not read the top-layer packet's own "
                "film-free surface n-hexane activity "
                "film_free_surface_hexane_activity: "
                f"{error}"
            ) from error
        hexane_ratio = binding.hexane_activity
        # The regime is part of the basis: the ledger must say WHICH branch of
        # the regime map priced this packet, not merely that one did.
        hexane_basis = "regime-map-surface-activity-" + binding.regime.value.lower().replace(
            "_", "-"
        )
    if decoded.surface.free_water_kg > 0.0:
        water_ratio = 1.0
        water_basis = "free-water-film"
    else:
        try:
            water_ratio = sorption.water_activity(
                decoded.state.retained_water_loading, PD_CONDENSATION_LUIKOV_PARAMS
            )
        except ValueError as error:
            raise PDKDomeBypassStepError(
                "cold-surface condensation could not read the top-layer packet's own "
                f"retained-water activity: {error}"
            ) from error
        water_basis = "modified-luikov-retained-water-activity"
    return (
        hexane_ratio,
        water_ratio,
        f"hexane:{hexane_basis}/water:{water_basis}/"
        + PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID,
    )


def _representative_weighted_mean(
    values: tuple[float, ...],
    counts: tuple[float, ...],
    total_count: float,
) -> float:
    """The representative-weighted mean, anchored for stable summation.

    REPORTING ONLY.  Nothing in the condensation criterion is evaluated on a
    mean any more (F-CS1-2); these aggregates exist so the layer ledger can
    still be read at a glance beside its per-packet rows.  On a layer whose
    packets all carry the same value it returns that value EXACTLY.
    """

    anchor = min(values)
    return (
        anchor
        + math.fsum(
            count * (value - anchor)
            for value, count in zip(values, counts, strict=True)
        )
        / total_count
    )


def _finite_weighted_mean(
    values: tuple[float, ...],
    counts: tuple[float, ...],
) -> float:
    """``_representative_weighted_mean`` over the FINITE entries only.

    A component that is not condensing on a packet has no stop temperature
    there and carries ``nan``; those packets are excluded rather than
    poisoning the aggregate.  ``nan`` when no packet has a finite value.
    """

    kept = tuple(
        (value, count)
        for value, count in zip(values, counts, strict=True)
        if math.isfinite(value)
    )
    if not kept:
        return math.nan
    return _representative_weighted_mean(
        tuple(value for value, _ in kept),
        tuple(count for _, count in kept),
        math.fsum(count for _, count in kept),
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDPacketCondensationRow:
    """CS-1a rev 2: what ONE top-layer packet's own cold surface condensed.

    The priority brief's finding F-CS1-2 is that the driving force must be
    evaluated per packet, on that packet's own temperature and its own
    equilibrium state, each packet capped at its OWN equilibrium, and the
    layer's ledger summed over packets.  This row IS that evaluation, kept as
    a measurement so a reader can see which packet condensed and why.
    """

    packet_id: str
    vertical_layer: int
    representative_particles: float
    surface_temperature_k: float
    contact_area_share_m2: float
    declared_energy_share_j: float
    hexane_equilibrium_ratio: float
    water_equilibrium_ratio: float
    equilibrium_ratio_basis_id: str
    hexane_surface_liquid_fugacity_pa: float
    water_surface_liquid_fugacity_pa: float
    hexane_driving_force_pa: float
    water_driving_force_pa: float
    hexane_stop_temperature_k: float
    water_stop_temperature_k: float
    stop_temperature_k: float
    stop_component: str
    hexane_contact_limited_capacity_mol: float
    water_contact_limited_capacity_mol: float
    uncapped_common_energy_j: float
    equilibrium_energy_cap_j: float
    equilibrium_energy_scaling_factor: float
    hexane_condensed_mol: float
    water_condensed_mol: float
    hexane_condensed_kg: float
    water_condensed_kg: float
    condensation_common_energy_j: float

    physically_qualifying: ClassVar[bool] = False
    driving_force_evaluated_on_this_packet: ClassVar[bool] = True


def _component_stop_temperature_k(
    module,
    *,
    bulk_fugacity_pa: float,
    equilibrium_ratio: float,
    pressure_pa: float,
    low_k: float,
    high_k: float,
) -> float:
    """The temperature at which condensation of one component STOPS.

    The root of ``f_eq f^{L,pure}(T, P) = f^V_bulk`` by the same monotone
    bisection ``_cached_hexane_saturation_temperature_k`` already uses.  The
    surface ratio is held at its measured surface value - the declared
    treatment; it is exactly ``1.0`` on every free-film surface, where this is
    the component's own dew point in the fugacity convention.
    """

    def residual(temperature_k: float) -> float:
        return math.fsum(
            (
                equilibrium_ratio
                * azeotrope._pure_liquid_fugacity(module, temperature_k, pressure_pa),
                -bulk_fugacity_pa,
            )
        )

    low = low_k
    high = high_k
    if residual(low) > 0.0 or residual(high) < 0.0:
        raise PDKDomeBypassStepError(
            "the arriving vapour has no surface-equilibrium stop temperature for "
            f"{module.__name__.rsplit('.', 1)[-1]} inside [{low_k}, {high_k}] K; "
            "refuse, do not extrapolate"
        )
    for _ in range(100):
        mid = 0.5 * (low + high)
        if mid in (low, high):
            break
        if residual(mid) < 0.0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


#: The bisection brackets for the two stop temperatures.  Hexane's ceiling is
#: the frozen a_h = 1 ceiling discussed at ``heteroazeotrope.py:29-37``: on a
#: FREE-FILM surface the hexane ratio is exactly 1.0, so the stop temperature
#: is a true dew point and lies below the pure-hexane boiling point at any
#: pressure this tray sees.
#:
#: REV 2 CORRECTION.  Rev 1's comment here justified the narrow hexane ceiling
#: by saying the term REFUSES a dry-of-hexane surface, so ``a_h`` is always
#: 1.0.  That is no longer true: the sorbate branch is implemented, and a
#: suppressed ``a_h < 1`` LIFTS the hexane stop temperature exactly as a
#: suppressed ``a_w`` lifts water's, so the root can leave ``(300.0, 343.0)``.
#: Rev 2 deliberately did NOT widen it on this build's own authority and
#: recorded the ceremony open item D-2: "a sub-critical-codec successor will
#: reach it and it needs a ruled hexane bracket first."
#:
#: REV 3, RULED AND WIDENED 2026-09-06 (owner ruling of that night on the
#: peer-review audit ``GT_PS2_PEER_REVIEW_AUDIT_THREE_RULINGS_2026-09-06.md``
#: section 3.3 item iii: "rule the hexane stop-temperature bracket (ceremony
#: D-2) in the same packet ... because on a hexane-rich plenum stream no root
#: exists in (300, 343) K below a_h of about 0.8").  MEASURED on the flagship
#: PD3 stream (bulk hexane fugacity 76 694.558 Pa at 101 325 Pa): the root is
#: 341.619 K at a_h = 0.8 and the old bracket REFUSED every a_h below that,
#: while the widened bracket returns 346.103 K at 0.7, 358.124 K at 0.5 and
#: 437.385 K at 0.1.  The a_h = 1 free-film dew point is BIT-IDENTICAL on both
#: brackets (334.45252025669856 K), so no free-film number moves.
#:
#: REV 3's CEILING WAS AN EQUATION-OF-STATE STANDOFF, NOT A PHYSICAL LIMIT.
#: Rev 3 took the ceiling from ``heteroazeotrope._LIQUID_REFERENCE_MAX_K``
#: (465.0 K), a private constant authored for a different equality and
#: carrying no test of its own: an undocumented 5.93 K margin below the
#: MEASURED collapse of the n-hexane pure-liquid reference at 470.93 K at
#: 1 atm.  That constant says where the reference stops EXISTING.  It says
#: nothing about where this model may be EVALUATED.
#:
#: REV 4, RULED 2026-09-06 (owner, verbatim: "Ok 465 K is clear, proceed as
#: recommended", on the ledger entry "THE 465 K BRACKET CEILING: ADVERSARIAL
#: PASS DONE" in ``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md``).  The
#: ceiling is PHY-052's.  The surface-equilibrium ratio this bisection holds
#: fixed comes from the Cardarelli-Crapiste native-soybean GAB isotherm,
#: measured over 323.15-368.15 K, and the FROZEN decision PHY-052 authorises
#: it to 433 K and no further, with out-of-source-domain flags.  The endpoint
#: is single-sourced from ``SORPTION_ISOTHERM_EXTRAPOLATION_CEILING_K``
#: above, which names PHY-052; no number is typed here and the
#: equation-of-state standoff is no longer read.
#:
#: THE 380 K CALORIC CONSUMPTION LIMIT SITS BELOW BOTH.  A stop temperature
#: is only useful if the energy cap can consume it, and the wet particle core
#: is not evaluable above ``wet_core.WetCoreParams.T_max`` = 380.0 K
#: (``particle/wet_core.py:85``; ``particle/sphere.py:244, 309, 482, 521,
#: 556, 644`` bracket, refuse or clamp at it).  Roots between 380 K and 433 K
#: are admitted by THIS bracket and still fail the caloric bracket
#: downstream: the bracket is the outer fence, never a promise that the state
#: proceeds.
#:
#: MEASURED COVERAGE COST of 465.0 -> 433.0 K, on the flagship PD3 stream
#: (bulk hexane fugacity 76 694.55808390728 Pa at 101 325 Pa).  The LOWEST
#: surface activity that still has a root is 0.767 at a 343 K ceiling, 0.291
#: at 380 K, 0.107 at 433 K and 0.069 at 465 K.  The ruled ceiling therefore
#: keeps essentially the whole benefit of Rev 3 - down to a_h = 0.107, inside
#: a frozen decision - and gives up only the 0.107 > a_h >= 0.069 sliver that
#: lived above PHY-052's authority.  Every root Rev 3 celebrated is
#: BIT-IDENTICAL on the ruled bracket, because the bisection runs to float
#: convergence rather than to a tolerance: 334.45252025669856 K at a_h = 1,
#: 341.61887561626804 K at 0.8, 346.10288174357026 K at 0.7,
#: 358.1244135253279 K at 0.5, 380.0422450127247 K at 0.291,
#: 432.9456091398554 K at 0.107.  The a_h = 0.1 root at 437.3845950544692 K
#: now REFUSES typed.
#:
#: THE REFERENCE ABOVE THE NORMAL BOILING POINT IS A METASTABLE SUPERHEATED
#: LIQUID BRANCH, not a compressed liquid.  n-Hexane boils at 341.9 K at
#: 1 atm, so at the ruled ceiling the pure-liquid reference this bisection
#: evaluates is 91.1 K superheated and its 101 325 Pa reference pressure is
#: 11.2 % of its own saturation pressure there (908 162.097 Pa); at Rev 3's
#: 465 K it was 123.1 K superheated at 6.4 %.  Rev 3 called it "compressed
#: liquid" throughout, which is wrong above 341.9 K; the wording is corrected
#: here and the physics of the extrapolation is unchanged.
#:
#: Fail-closed is unchanged: where the root lies outside the bracket,
#: ``_component_stop_temperature_k`` REFUSES typed ("refuse, do not
#: extrapolate") for a CONDENSING component and records ``nan`` for a
#: non-condensing one.  A strongly suppressed surface still refuses - on the
#: same flagship stream the regime-C activity 0.059 has no root below the
#: ceiling and the step is rejected, not extrapolated.
#:
#: BOTH BRACKETS ARE DIGEST-BEARING FROM THIS REV, answering the adversarial
#: pass's point 8: a constant that decides which states proceed belongs in
#: the model identity.  They enter the ARMED component of
#: ``PDKTrayEngineeringModel.definition_digest`` only, so every unarmed
#: digest - the committed production default, contact fraction 0.0 - stays
#: bit-identical to the pre-WS-A one.
#:
#: Water's ceiling is NOT its own dew point: on a sorption-limited surface
#: (a_w < 1, the Luikov branch) the stop temperature is LIFTED above the
#: free-water dew point - the same sorption-elevated endpoint the T3b record
#: measured - so the bracket must carry it.  ``473.15`` K is still more than
#: 100 K below the temperature at which the frozen water authority's
#: metastable superheated liquid branch collapses, which
#: ``_pure_liquid_fugacity`` verifies on every call.  This rev does NOT move
#: water's ceiling; it only stamps it into the digest.
PD_CONDENSATION_STOP_BRACKET_HEXANE_K = (300.0, SORPTION_ISOTHERM_EXTRAPOLATION_CEILING_K)
PD_CONDENSATION_STOP_BRACKET_WATER_K = (300.0, 473.15)


def _contact_limited_condensation(
    stream: PDBinaryGasStream,
    *,
    model: PDKTrayEngineeringModel,
    duration_s: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
    upper_packets: tuple[tray_host.AcceptedTrayPacket, ...] = (),
    upper_layer_declared_energy_by_packet: dict[str, float] | None = None,
) -> PDContactCondensationLedger | None:
    """CS-1a: what this tray's cold top free surface condenses, PER PACKET.

    ``None`` means the term is inert for this step - the declared contact
    fraction is 0.0 (the committed production default, where the caller takes
    the pre-WS-A path verbatim), the stream is above a critical temperature,
    or NO packet of the receiving layer ends up condensing anything.

    **REV 2, 2026-09-04 (priority brief S:6, finding F-CS1-2).**  Rev 1 took a
    single ``surface_temperature_k`` - the caller handed it the
    representative-weighted MEAN temperature of the whole top layer - and
    evaluated one driving force on it.  That averages the cold feed packets
    into the warm resident population and destroys the driving force this
    capability exists to represent: a layer holding 99.94 % of its mass at
    67 C and 0.06 % at 50 C must condense on the cold fraction.  The mean
    parameter is GONE.  Each packet is now priced on ITS OWN decoded
    temperature and ITS OWN equilibrium ratios, capped at ITS OWN first stop
    temperature against ITS OWN share of the declared energy, and the layer's
    ledger is the SUM over packets.  The scalar fields of the returned ledger
    are representative-weighted REPORTS of those rows, never the criterion.
    """

    contact_fraction = model.condensation_contact_fraction
    if contact_fraction == 0.0:
        return None
    if not model.flash_or_condensation_implemented:
        # Q-PD-4, MADE EXECUTABLE 2026-09-04.  The construction gate already
        # refuses an armed model; this is the same fence at the one place the
        # term could otherwise fire, so no step can condense while the record
        # it stamps declares the capability unimplemented.
        raise PDKDomeBypassStepError(
            "the contact-limited dome condensation term requires the declared "
            "flash_or_condensation_implemented capability flag (owner ruling Q-PD-4)"
        )
    if stream.temperature_k >= water_properties.TC or stream.temperature_k >= hexane.TC:
        return None
    if not upper_packets:
        # FAIL CLOSED.  CS-1a prices condensation against the COLD SOLID, so a
        # caller that does not hand the term its surface has not asked a
        # well-posed question; it never gets a silently stream-driven answer.
        raise PDKDomeBypassStepError(
            "cold-surface condensation requires the top-layer surface it condenses "
            "on (the layer's own packets, each carrying its own temperature)"
        )
    declared_energy_by_packet = dict(upper_layer_declared_energy_by_packet or {})
    for value in declared_energy_by_packet.values():
        _require_float(
            "cold-surface condensation upper-layer declared energy share", value
        )

    # --- the BULK vapour leg, one per step: it is a property of the stream --
    gas = binary_gas.state(
        stream.temperature_k,
        stream.pressure_pa,
        stream.water_mole_fraction,
        stream.hexane_mole_fraction,
        k_wh=binary_gas.K_WH_CENTRAL,
    )
    hexane_bulk_fugacity = gas.fugacity_hexane_pa
    water_bulk_fugacity = gas.fugacity_water_pa

    # --- the gas-side numbers the FOG GUARD reads (not the criterion) -----
    partial_pressure = stream.pressure_pa * stream.water_mole_fraction
    saturation_pressure = water_properties.saturation_pressure(stream.temperature_k)
    excess = math.fsum((partial_pressure, -saturation_pressure))
    total_mol = stream.molar_flow_mol_s * duration_s
    water_mol = total_mol * stream.water_mole_fraction
    hexane_mol = total_mol * stream.hexane_mole_fraction
    saturated_share = saturation_pressure / stream.pressure_pa
    relief_mol = (
        math.fsum((water_mol, -saturated_share * total_mol))
        / math.fsum((1.0, -saturated_share))
        if excess > 0.0 and saturated_share < 1.0
        else 0.0
    )

    # --- the rate coefficient (WS-A's, verbatim; one k_G, both components) -
    contact_area = _pd_dome_contact_area_m2(model.geometry)
    gas_film_coefficient = model.wall.wall_to_bypass_gas_ua_w_k / contact_area
    molar_heat_capacity = _stream_molar_heat_capacity_j_mol_k(
        model.properties,
        temperature_k=stream.temperature_k,
        pressure_pa=stream.pressure_pa,
        hexane_mole_fraction=stream.hexane_mole_fraction,
    )
    if not (math.isfinite(molar_heat_capacity) and molar_heat_capacity > 0.0):
        raise PDKDomeBypassStepError(
            "the condensation rate law needs a positive stream molar heat capacity"
        )
    # Lewis relation at Le = 1: k_G = h / (P * cp_molar).
    molar_conductance = gas_film_coefficient / (stream.pressure_pa * molar_heat_capacity)

    # --- the F2 seam price, per component, at the STREAM state ------------
    datum_adapter = codec.component_datum_adapter
    hexane_liquid_enthalpy = datum_adapter.common_datum_liquid_enthalpy_j_mol(
        "hexane", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    hexane_seam_latent = datum_adapter.seam_latent_j_mol(
        "hexane", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    water_liquid_enthalpy = datum_adapter.common_datum_liquid_enthalpy_j_mol(
        "water", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    water_seam_latent = datum_adapter.seam_latent_j_mol(
        "water", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    # F2, THE REVERSE CROSSING, priced once per species.  Evaporation takes
    # the liquid enthalpy OUT of the packet and CONSUMES the seam latent;
    # condensation is the same seam with both signs flipped, so the packet
    # receives the liquid it becomes PLUS the latent that crossing releases.
    hexane_energy_per_mole = math.fsum((hexane_liquid_enthalpy, hexane_seam_latent))
    water_energy_per_mole = math.fsum((water_liquid_enthalpy, water_seam_latent))

    def stop_temperature_of(module, bulk, ratio, bracket, condensing):
        try:
            return _component_stop_temperature_k(
                module,
                bulk_fugacity_pa=bulk,
                equilibrium_ratio=ratio,
                pressure_pa=stream.pressure_pa,
                low_k=bracket[0],
                high_k=bracket[1],
            )
        except PDKDomeBypassStepError:
            if condensing:
                raise
            # A component that is NOT condensing needs no endpoint; its stop
            # temperature is recorded as nan and is never read.
            return math.nan

    # --- ONE ROW PER PACKET: its own surface, its own cap ------------------
    counts = tuple(packet.entry.weight.representative_particles for packet in upper_packets)
    total_count = math.fsum(counts)
    if not total_count > 0.0:
        raise PDKDomeBypassStepError(
            "the dome contact surface has no represented top-layer particle"
        )
    # The dome contact AREA is shared by representative-particle count through
    # the same exact-residual allocator the inter-K heat uses, so the shares
    # sum to the contact area exactly.  What is NOT shared is the driving
    # force: every packet computes its own.
    area_by_packet = _allocate_layer_heat(upper_packets, contact_area)

    rows: list[PDPacketCondensationRow] = []
    for packet, count in zip(upper_packets, counts, strict=True):
        decoded = codec.decode(packet.entry.payload_bytes)
        surface_temperature = decoded.state.temperature_k
        _require_float(
            "cold-surface condensation packet surface temperature",
            surface_temperature,
            positive=True,
        )
        hexane_ratio, water_ratio, ratio_basis = _packet_surface_equilibrium_ratios(
            decoded,
            temperature_k=surface_temperature,
            pressure_pa=stream.pressure_pa,
        )
        try:
            hexane_surface_fugacity = azeotrope._pure_liquid_fugacity(
                hexane, surface_temperature, stream.pressure_pa
            )
            water_surface_fugacity = azeotrope._pure_liquid_fugacity(
                water_properties, surface_temperature, stream.pressure_pa
            )
        except azeotrope.HeteroazeotropeError as error:
            raise PDKDomeBypassStepError(
                f"the cold surface has no pure compressed-liquid reference state: {error}"
            ) from error
        hexane_driving_force = math.fsum(
            (hexane_bulk_fugacity, -hexane_ratio * hexane_surface_fugacity)
        )
        water_driving_force = math.fsum(
            (water_bulk_fugacity, -water_ratio * water_surface_fugacity)
        )
        hexane_stop_k = stop_temperature_of(
            hexane,
            hexane_bulk_fugacity,
            hexane_ratio,
            PD_CONDENSATION_STOP_BRACKET_HEXANE_K,
            hexane_driving_force > 0.0,
        )
        water_stop_k = stop_temperature_of(
            water_properties,
            water_bulk_fugacity,
            water_ratio,
            PD_CONDENSATION_STOP_BRACKET_WATER_K,
            water_driving_force > 0.0,
        )
        area_share = area_by_packet[packet.packet_id]
        declared_share = declared_energy_by_packet.get(packet.packet_id, 0.0)
        capacity_scale = contact_fraction * area_share * molar_conductance * duration_s
        hexane_capacity_mol = (
            capacity_scale * hexane_driving_force if hexane_driving_force > 0.0 else 0.0
        )
        water_capacity_mol = (
            capacity_scale * water_driving_force if water_driving_force > 0.0 else 0.0
        )
        uncapped_energy_j = math.fsum(
            (
                hexane_capacity_mol * hexane_energy_per_mole,
                water_capacity_mol * water_energy_per_mole,
            )
        )
        stop_pairs = tuple(
            (temperature, species)
            for temperature, species, condensing in (
                (hexane_stop_k, "n-hexane", hexane_driving_force > 0.0),
                (water_stop_k, "water", water_driving_force > 0.0),
            )
            if condensing
        )
        if not stop_pairs:
            # This packet's surface is at or above both stop temperatures: it
            # condenses nothing, and it does not silence its colder neighbours.
            rows.append(
                PDPacketCondensationRow(
                    packet_id=packet.packet_id,
                    vertical_layer=packet.owner_key.vertical_layer_id.value,
                    representative_particles=count,
                    surface_temperature_k=surface_temperature,
                    contact_area_share_m2=area_share,
                    declared_energy_share_j=declared_share,
                    hexane_equilibrium_ratio=hexane_ratio,
                    water_equilibrium_ratio=water_ratio,
                    equilibrium_ratio_basis_id=ratio_basis,
                    hexane_surface_liquid_fugacity_pa=hexane_surface_fugacity,
                    water_surface_liquid_fugacity_pa=water_surface_fugacity,
                    hexane_driving_force_pa=hexane_driving_force,
                    water_driving_force_pa=water_driving_force,
                    hexane_stop_temperature_k=hexane_stop_k,
                    water_stop_temperature_k=water_stop_k,
                    stop_temperature_k=math.nan,
                    stop_component="",
                    hexane_contact_limited_capacity_mol=0.0,
                    water_contact_limited_capacity_mol=0.0,
                    uncapped_common_energy_j=0.0,
                    equilibrium_energy_cap_j=math.nan,
                    equilibrium_energy_scaling_factor=1.0,
                    hexane_condensed_mol=0.0,
                    water_condensed_mol=0.0,
                    hexane_condensed_kg=0.0,
                    water_condensed_kg=0.0,
                    condensation_common_energy_j=0.0,
                )
            )
            continue
        stop_temperature_k, stop_component = min(stop_pairs)
        # THE CAP IS THIS PACKET'S OWN SURFACE at its own first stop
        # temperature, less the explicit declared energy this step already
        # deposits into THIS packet (its share of the inter-K face heat plus
        # its own declared mechanical work).  Never gas saturation.
        energy_to_stop = math.fsum(
            (
                _packet_energy_at_temperature_j(
                    packet, stop_temperature_k, codec=codec
                ),
                -packet.weighted_totals.totals[6],
            )
        )
        equilibrium_energy_cap_j = math.fsum((energy_to_stop, -declared_share))
        if equilibrium_energy_cap_j <= 0.0:
            scaling_factor = 0.0
        elif uncapped_energy_j > equilibrium_energy_cap_j:
            scaling_factor = equilibrium_energy_cap_j / uncapped_energy_j
        else:
            scaling_factor = 1.0
        hexane_condensed_mol = hexane_capacity_mol * scaling_factor
        water_condensed_mol = water_capacity_mol * scaling_factor
        rows.append(
            PDPacketCondensationRow(
                packet_id=packet.packet_id,
                vertical_layer=packet.owner_key.vertical_layer_id.value,
                representative_particles=count,
                surface_temperature_k=surface_temperature,
                contact_area_share_m2=area_share,
                declared_energy_share_j=declared_share,
                hexane_equilibrium_ratio=hexane_ratio,
                water_equilibrium_ratio=water_ratio,
                equilibrium_ratio_basis_id=ratio_basis,
                hexane_surface_liquid_fugacity_pa=hexane_surface_fugacity,
                water_surface_liquid_fugacity_pa=water_surface_fugacity,
                hexane_driving_force_pa=hexane_driving_force,
                water_driving_force_pa=water_driving_force,
                hexane_stop_temperature_k=hexane_stop_k,
                water_stop_temperature_k=water_stop_k,
                stop_temperature_k=stop_temperature_k,
                stop_component=stop_component,
                hexane_contact_limited_capacity_mol=hexane_capacity_mol,
                water_contact_limited_capacity_mol=water_capacity_mol,
                uncapped_common_energy_j=uncapped_energy_j,
                equilibrium_energy_cap_j=equilibrium_energy_cap_j,
                equilibrium_energy_scaling_factor=scaling_factor,
                hexane_condensed_mol=hexane_condensed_mol,
                water_condensed_mol=water_condensed_mol,
                hexane_condensed_kg=(
                    hexane_condensed_mol * model.properties.hexane_molar_mass_kg_mol
                ),
                water_condensed_kg=(
                    water_condensed_mol * model.properties.water_molar_mass_kg_mol
                ),
                condensation_common_energy_j=math.fsum(
                    (
                        hexane_condensed_mol * hexane_energy_per_mole,
                        water_condensed_mol * water_energy_per_mole,
                    )
                ),
            )
        )

    packet_rows = tuple(rows)
    hexane_condensed_mol = math.fsum(row.hexane_condensed_mol for row in packet_rows)
    water_condensed_mol = math.fsum(row.water_condensed_mol for row in packet_rows)
    if hexane_condensed_mol <= 0.0 and water_condensed_mol <= 0.0:
        return None
    if hexane_condensed_mol > 0.0 and hexane_condensed_mol >= hexane_mol:
        raise PDKDomeBypassStepError(
            "cold-surface condensation would remove the whole n-hexane inventory of "
            "the bypass stream; refuse, do not clamp"
        )
    if water_condensed_mol > 0.0 and water_condensed_mol >= water_mol:
        raise PDKDomeBypassStepError(
            "cold-surface condensation would remove the whole water inventory of "
            "the bypass stream; refuse, do not clamp"
        )

    hexane_condensed_kg = math.fsum(row.hexane_condensed_kg for row in packet_rows)
    water_condensed_kg = math.fsum(row.water_condensed_kg for row in packet_rows)
    hexane_common_energy_j = math.fsum(
        row.hexane_condensed_mol * hexane_energy_per_mole for row in packet_rows
    )
    water_common_energy_j = math.fsum(
        row.water_condensed_mol * water_energy_per_mole for row in packet_rows
    )
    condensation_common_energy_j = math.fsum(
        row.condensation_common_energy_j for row in packet_rows
    )
    uncapped_total_energy_j = math.fsum(
        row.uncapped_common_energy_j for row in packet_rows
    )
    scaling_factor = (
        condensation_common_energy_j / uncapped_total_energy_j
        if uncapped_total_energy_j > 0.0
        else 1.0
    )
    capped_rows = tuple(
        row for row in packet_rows if math.isfinite(row.equilibrium_energy_cap_j)
    )
    equilibrium_energy_cap_j = math.fsum(
        row.equilibrium_energy_cap_j for row in capped_rows
    )
    condensing_stops = tuple(
        (row.stop_temperature_k, row.stop_component)
        for row in packet_rows
        if row.stop_component
    )
    stop_temperature_k, stop_component = min(condensing_stops)

    remaining_total_mol = math.fsum(
        (total_mol, -hexane_condensed_mol, -water_condensed_mol)
    )
    remaining_water_mol = math.fsum((water_mol, -water_condensed_mol))
    residual_supersaturation = math.fsum(
        (
            stream.pressure_pa * remaining_water_mol / remaining_total_mol,
            -saturation_pressure,
        )
    )
    hexane_bases = tuple(
        dict.fromkeys(
            row.equilibrium_ratio_basis_id.split("/")[0].removeprefix("hexane:")
            for row in packet_rows
        )
    )
    water_bases = tuple(
        dict.fromkeys(
            row.equilibrium_ratio_basis_id.split("/")[1].removeprefix("water:")
            for row in packet_rows
        )
    )
    ratio_basis_id = (
        "hexane:"
        + "+".join(hexane_bases)
        + "/water:"
        + "+".join(water_bases)
        + "/"
        + PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID
    )
    surface_temperatures = tuple(row.surface_temperature_k for row in packet_rows)
    return PDContactCondensationLedger(
        contact_fraction=contact_fraction,
        contact_area_m2=contact_area,
        contact_area_basis_id=PD_CONDENSATION_CONTACT_AREA_BASIS_ID,
        rate_law_id=PD_CONDENSATION_RATE_LAW_ID,
        energy_routing_id=PD_CONDENSATION_ENERGY_ROUTING_ID,
        receiver_law_id=PD_CONDENSATION_RECEIVER_LAW_ID,
        cap_law_id=PD_CONDENSATION_CAP_LAW_ID,
        transfer_coefficient_law_id=PD_CONDENSATION_TRANSFER_COEFFICIENT_LAW_ID,
        equilibrium_ratio_law_id=PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID,
        equilibrium_ratio_basis_id=ratio_basis_id,
        duration_s=duration_s,
        stream_temperature_k=stream.temperature_k,
        stream_pressure_pa=stream.pressure_pa,
        surface_temperature_k=_representative_weighted_mean(
            surface_temperatures, counts, total_count
        ),
        packet_rows=packet_rows,
        gas_film_coefficient_w_m2_k=gas_film_coefficient,
        molar_heat_capacity_j_mol_k=molar_heat_capacity,
        molar_transfer_conductance_mol_s_m2_pa=molar_conductance,
        hexane_bulk_fugacity_pa=hexane_bulk_fugacity,
        water_bulk_fugacity_pa=water_bulk_fugacity,
        hexane_surface_liquid_fugacity_pa=_representative_weighted_mean(
            tuple(row.hexane_surface_liquid_fugacity_pa for row in packet_rows),
            counts,
            total_count,
        ),
        water_surface_liquid_fugacity_pa=_representative_weighted_mean(
            tuple(row.water_surface_liquid_fugacity_pa for row in packet_rows),
            counts,
            total_count,
        ),
        hexane_equilibrium_ratio=_representative_weighted_mean(
            tuple(row.hexane_equilibrium_ratio for row in packet_rows),
            counts,
            total_count,
        ),
        water_equilibrium_ratio=_representative_weighted_mean(
            tuple(row.water_equilibrium_ratio for row in packet_rows),
            counts,
            total_count,
        ),
        hexane_driving_force_pa=_representative_weighted_mean(
            tuple(row.hexane_driving_force_pa for row in packet_rows),
            counts,
            total_count,
        ),
        water_driving_force_pa=_representative_weighted_mean(
            tuple(row.water_driving_force_pa for row in packet_rows),
            counts,
            total_count,
        ),
        hexane_dew_point_k=_finite_weighted_mean(
            tuple(row.hexane_stop_temperature_k for row in packet_rows), counts
        ),
        water_dew_point_k=_finite_weighted_mean(
            tuple(row.water_stop_temperature_k for row in packet_rows), counts
        ),
        condensation_equilibrium_temperature_k=stop_temperature_k,
        stop_component=stop_component,
        hexane_contact_limited_capacity_mol=math.fsum(
            row.hexane_contact_limited_capacity_mol for row in packet_rows
        ),
        contact_limited_capacity_mol=math.fsum(
            row.water_contact_limited_capacity_mol for row in packet_rows
        ),
        upper_layer_declared_energy_j=math.fsum(row.declared_energy_share_j for row in packet_rows),
        equilibrium_energy_cap_j=equilibrium_energy_cap_j,
        equilibrium_energy_scaling_factor=scaling_factor,
        hexane_condensed_mol=hexane_condensed_mol,
        hexane_condensed_kg=hexane_condensed_kg,
        condensed_water_mol=water_condensed_mol,
        condensed_water_kg=water_condensed_kg,
        hexane_common_datum_liquid_enthalpy_j_mol=hexane_liquid_enthalpy,
        hexane_seam_latent_j_mol=hexane_seam_latent,
        common_datum_liquid_enthalpy_j_mol=water_liquid_enthalpy,
        seam_latent_j_mol=water_seam_latent,
        hexane_condensation_common_energy_j=hexane_common_energy_j,
        water_condensation_common_energy_j=water_common_energy_j,
        condensation_common_energy_j=condensation_common_energy_j,
        water_partial_pressure_pa=partial_pressure,
        water_saturation_pressure_pa=saturation_pressure,
        partial_pressure_excess_pa=excess,
        equilibrium_relief_mol=relief_mol,
        residual_supersaturation_pa=residual_supersaturation,
    )


class PDSteamCondensingLevel0Refusal(PDKDomeBypassStepError):
    """D12 steam condensing level 0: a condensate the level-0 routing cannot carry."""


def _steam_condensing_slice_mole_fraction(
    stream: PDBinaryGasStream,
    *,
    species: str,
    target_fugacity_pa: float,
) -> float:
    """D12 steam condensing level 0: one species' equilibrium mole fraction in a slice.

    The smallest mole fraction y of ``species`` in [0, y_stream] whose partial
    fugacity y phi(T, P, y) P, at the STREAM temperature and pressure on the
    D12-a vapour leg (``props.binary_gas`` at ``K_WH_CENTRAL``), is not below
    ``target_fugacity_pa``: bisection to adjacent floats, keeping the end on the
    gas side of the equilibrium, so the removal it implies never passes it.
    """

    water = species == "water"
    upper = stream.water_mole_fraction if water else stream.hexane_mole_fraction

    def fugacity(fraction: float) -> float:
        y_water = fraction if water else 1.0 - fraction
        state = binary_gas.state(
            stream.temperature_k,
            stream.pressure_pa,
            y_water,
            1.0 - y_water,
            k_wh=binary_gas.K_WH_CENTRAL,
        )
        return state.fugacity_water_pa if water else state.fugacity_hexane_pa

    if not fugacity(upper) > target_fugacity_pa:
        raise PDKDomeBypassStepError(
            f"steam condensing level 0: the {species} slice equilibrium was asked for a "
            "stream that is not above the surface fugacity; refuse, do not clamp"
        )
    low, high = 0.0, upper
    while True:
        middle = 0.5 * (low + high)
        if middle in (low, high):
            return high
        if fugacity(middle) >= target_fugacity_pa:
            high = middle
        else:
            low = middle


def _steam_condensing_guard_mole_fraction(
    stream: PDBinaryGasStream,
    *,
    species: str,
) -> float:
    """D12 steam condensing level 0: the D12-a guard's bound on one species at the gas state.

    The largest mole fraction y of ``species`` in [0, 1] whose certified excess (partial
    fugacity at the STREAM temperature and pressure less the pure-liquid fugacity at the
    same state) is not positive: the composition limit ``_require_single_vapour_phase``
    admits for that species.  1.0 when no fraction reaches it.  Bisection to adjacent
    floats, keeping the admitted end.
    """

    water = species == "water"
    module = water_properties if water else hexane
    try:
        liquid = azeotrope._pure_liquid_fugacity(module, stream.temperature_k, stream.pressure_pa)
    except azeotrope.HeteroazeotropeError as error:
        raise PDKDomeBypassStepError(
            f"steam condensing level 0: the D12-a {species} leg could not be evaluated at "
            f"the gas state: {type(error).__name__}: {error}"
        ) from error

    def excess(fraction: float) -> float:
        y_water = fraction if water else 1.0 - fraction
        state = binary_gas.state(
            stream.temperature_k,
            stream.pressure_pa,
            y_water,
            1.0 - y_water,
            k_wh=binary_gas.K_WH_CENTRAL,
        )
        fugacity = state.fugacity_water_pa if water else state.fugacity_hexane_pa
        return math.fsum((fugacity, -liquid))

    if excess(1.0) <= 0.0:
        return 1.0
    low, high = 0.0, 1.0
    while True:
        middle = 0.5 * (low + high)
        if middle in (low, high):
            return low
        if excess(middle) <= 0.0:
            low = middle
        else:
            high = middle


def _steam_condensing_level0(
    stream: PDBinaryGasStream,
    *,
    model: PDKTrayEngineeringModel,
    duration_s: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    declared_energy_by_packet: dict[str, float],
) -> PDContactCondensationLedger | None:
    """D12 steam condensing level 0: the rising vapour condensing on every PD packet.

    See ``PD_STEAM_CONDENSING_LEVEL0_ARMED`` for the law.  Returns ``None`` when no
    packet condenses anything, exactly as the CS-1a term does, so the caller then
    takes the path it takes without condensation.  The ledger is the CS-1a ledger
    type, so every consumer downstream is the shipped one; its rate-law fields that
    the stage form does not use (contact fraction, area, film coefficient, molar
    heat capacity, molar conductance) are ``nan``.
    """

    if stream.temperature_k >= water_properties.TC or stream.temperature_k >= hexane.TC:
        return None
    if not packets:
        raise PDKDomeBypassStepError("steam condensing level 0 needs the tray's packets")
    gas = binary_gas.state(
        stream.temperature_k,
        stream.pressure_pa,
        stream.water_mole_fraction,
        stream.hexane_mole_fraction,
        k_wh=binary_gas.K_WH_CENTRAL,
    )
    hexane_bulk_fugacity = gas.fugacity_hexane_pa
    water_bulk_fugacity = gas.fugacity_water_pa
    # The gas-side numbers the fog guard reads (not the criterion), as CS-1a.
    partial_pressure = stream.pressure_pa * stream.water_mole_fraction
    saturation_pressure = water_properties.saturation_pressure(stream.temperature_k)
    excess = math.fsum((partial_pressure, -saturation_pressure))
    total_mol = stream.molar_flow_mol_s * duration_s
    water_mol = total_mol * stream.water_mole_fraction
    hexane_mol = total_mol * stream.hexane_mole_fraction
    saturated_share = saturation_pressure / stream.pressure_pa
    relief_mol = (
        math.fsum((water_mol, -saturated_share * total_mol)) / math.fsum((1.0, -saturated_share))
        if excess > 0.0 and saturated_share < 1.0
        else 0.0
    )
    # The F2 seam price, per component, at the STREAM state (CS-1a, verbatim).
    datum_adapter = codec.component_datum_adapter
    hexane_liquid_enthalpy = datum_adapter.common_datum_liquid_enthalpy_j_mol(
        "hexane", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    hexane_seam_latent = datum_adapter.seam_latent_j_mol(
        "hexane", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    water_liquid_enthalpy = datum_adapter.common_datum_liquid_enthalpy_j_mol(
        "water", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    water_seam_latent = datum_adapter.seam_latent_j_mol(
        "water", temperature_k=stream.temperature_k, pressure_pa=stream.pressure_pa
    )
    hexane_energy_per_mole = math.fsum((hexane_liquid_enthalpy, hexane_seam_latent))
    water_energy_per_mole = math.fsum((water_liquid_enthalpy, water_seam_latent))

    def stop_temperature_of(module, bulk, ratio, bracket, condensing):
        try:
            return _component_stop_temperature_k(
                module,
                bulk_fugacity_pa=bulk,
                equilibrium_ratio=ratio,
                pressure_pa=stream.pressure_pa,
                low_k=bracket[0],
                high_k=bracket[1],
            )
        except PDKDomeBypassStepError:
            if condensing:
                raise
            return math.nan

    counts = tuple(packet.entry.weight.representative_particles for packet in packets)
    total_count = math.fsum(counts)
    if not total_count > 0.0:
        raise PDKDomeBypassStepError("steam condensing level 0: no represented PD particle")
    # Each packet's slice of the step's gas, by representative count, through the
    # exact-residual allocator: the slices sum to the stream exactly.
    water_by_packet = _allocate_layer_heat(packets, water_mol)
    hexane_by_packet = _allocate_layer_heat(packets, hexane_mol)
    equilibrium_cache: dict[tuple[str, float], float] = {}
    guard_cache: dict[str, float] = {}

    def guard_fraction(species: str) -> float:
        if species not in guard_cache:
            guard_cache[species] = _steam_condensing_guard_mole_fraction(stream, species=species)
        return guard_cache[species]

    def slice_equilibrium(species: str, target: float) -> float:
        key = (species, target)
        if key not in equilibrium_cache:
            equilibrium_cache[key] = _steam_condensing_slice_mole_fraction(
                stream, species=species, target_fugacity_pa=target
            )
        return equilibrium_cache[key]

    rows: list[PDPacketCondensationRow] = []
    for packet, count in zip(packets, counts, strict=True):
        decoded = codec.decode(packet.entry.payload_bytes)
        surface_temperature = decoded.state.temperature_k
        _require_float(
            "steam condensing packet surface temperature", surface_temperature, positive=True
        )
        hexane_ratio, water_ratio, ratio_basis = _packet_surface_equilibrium_ratios(
            decoded, temperature_k=surface_temperature, pressure_pa=stream.pressure_pa
        )
        try:
            hexane_surface_fugacity = azeotrope._pure_liquid_fugacity(
                hexane, surface_temperature, stream.pressure_pa
            )
            water_surface_fugacity = azeotrope._pure_liquid_fugacity(
                water_properties, surface_temperature, stream.pressure_pa
            )
        except azeotrope.HeteroazeotropeError as error:
            raise PDKDomeBypassStepError(
                f"the PD packet surface has no pure-liquid reference state: {error}"
            ) from error
        hexane_target = hexane_ratio * hexane_surface_fugacity
        water_target = water_ratio * water_surface_fugacity
        hexane_driving_force = math.fsum((hexane_bulk_fugacity, -hexane_target))
        water_driving_force = math.fsum((water_bulk_fugacity, -water_target))
        hexane_condensing = hexane_driving_force > 0.0
        water_condensing = water_driving_force > 0.0
        hexane_stop_k = stop_temperature_of(
            hexane,
            hexane_bulk_fugacity,
            hexane_ratio,
            PD_CONDENSATION_STOP_BRACKET_HEXANE_K,
            hexane_condensing,
        )
        water_stop_k = stop_temperature_of(
            water_properties,
            water_bulk_fugacity,
            water_ratio,
            PD_CONDENSATION_STOP_BRACKET_WATER_K,
            water_condensing,
        )
        slice_water = water_by_packet[packet.packet_id]
        slice_hexane = hexane_by_packet[packet.packet_id]
        water_capacity_mol = 0.0
        hexane_capacity_mol = 0.0
        if water_condensing:
            # Water first: down to the packet's surface fugacity, the removal kept inside
            # the n-hexane leg of the D12-a guard at the gas temperature (the slice's
            # n-hexane fraction may not rise past what the gas admits there).
            y_star = slice_equilibrium("water", water_target)
            equilibrium_removal = math.fsum((slice_water, -slice_hexane * y_star / (1.0 - y_star)))
            hexane_bound = guard_fraction("hexane")
            guard_removal = (
                math.inf
                if hexane_bound >= 1.0
                else math.fsum((slice_water, -slice_hexane * (1.0 - hexane_bound) / hexane_bound))
            )
            water_capacity_mol = max(0.0, min(equilibrium_removal, guard_removal))
            if water_capacity_mol >= slice_water:
                raise PDKDomeBypassStepError(
                    "steam condensing level 0 would remove the whole water of a packet's "
                    "gas slice; refuse, do not clamp"
                )
        slice_water_after = math.fsum((slice_water, -water_capacity_mol))
        if hexane_condensing:
            # Then n-hexane, "where the fugacity guard allows": with the water held at its
            # post-condensation amount, down to the smaller removal of its own surface
            # equilibrium and the water leg of the D12-a guard at the gas temperature.
            y_star = slice_equilibrium("hexane", hexane_target)
            equilibrium_removal = math.fsum(
                (slice_hexane, -slice_water_after * y_star / (1.0 - y_star))
            )
            water_bound = guard_fraction("water")
            guard_removal = (
                math.inf
                if water_bound >= 1.0
                else math.fsum(
                    (slice_hexane, -slice_water_after * (1.0 - water_bound) / water_bound)
                )
            )
            hexane_capacity_mol = max(0.0, min(equilibrium_removal, guard_removal))
            if hexane_capacity_mol >= slice_hexane:
                raise PDKDomeBypassStepError(
                    "steam condensing level 0 would remove the whole n-hexane of a "
                    "packet's gas slice; refuse, do not clamp"
                )
        declared_share = declared_energy_by_packet.get(packet.packet_id, 0.0)
        uncapped_energy_j = math.fsum(
            (
                hexane_capacity_mol * hexane_energy_per_mole,
                water_capacity_mol * water_energy_per_mole,
            )
        )
        stop_pairs = tuple(
            (temperature, species)
            for temperature, species, condensing in (
                (hexane_stop_k, "n-hexane", hexane_condensing and hexane_capacity_mol > 0.0),
                (water_stop_k, "water", water_condensing and water_capacity_mol > 0.0),
            )
            if condensing
        )
        common = dict(
            packet_id=packet.packet_id,
            vertical_layer=packet.owner_key.vertical_layer_id.value,
            representative_particles=count,
            surface_temperature_k=surface_temperature,
            contact_area_share_m2=math.nan,
            declared_energy_share_j=declared_share,
            hexane_equilibrium_ratio=hexane_ratio,
            water_equilibrium_ratio=water_ratio,
            equilibrium_ratio_basis_id=ratio_basis,
            hexane_surface_liquid_fugacity_pa=hexane_surface_fugacity,
            water_surface_liquid_fugacity_pa=water_surface_fugacity,
            hexane_driving_force_pa=hexane_driving_force,
            water_driving_force_pa=water_driving_force,
            hexane_stop_temperature_k=hexane_stop_k,
            water_stop_temperature_k=water_stop_k,
            hexane_contact_limited_capacity_mol=hexane_capacity_mol,
            water_contact_limited_capacity_mol=water_capacity_mol,
            uncapped_common_energy_j=uncapped_energy_j,
        )
        if not stop_pairs:
            rows.append(
                PDPacketCondensationRow(
                    **common,
                    stop_temperature_k=math.nan,
                    stop_component="",
                    equilibrium_energy_cap_j=math.nan,
                    equilibrium_energy_scaling_factor=1.0,
                    hexane_condensed_mol=0.0,
                    water_condensed_mol=0.0,
                    hexane_condensed_kg=0.0,
                    water_condensed_kg=0.0,
                    condensation_common_energy_j=0.0,
                )
            )
            continue
        stop_temperature_k, stop_component = min(stop_pairs)
        # The EXISTING first-stop cap (CS-1a): this packet's own energy to its own
        # first stop temperature, less the declared energy it already receives.
        energy_to_stop = math.fsum(
            (
                _packet_energy_at_temperature_j(packet, stop_temperature_k, codec=codec),
                -packet.weighted_totals.totals[6],
            )
        )
        equilibrium_energy_cap_j = math.fsum((energy_to_stop, -declared_share))
        if equilibrium_energy_cap_j <= 0.0:
            scaling_factor = 0.0
        elif uncapped_energy_j > equilibrium_energy_cap_j:
            scaling_factor = equilibrium_energy_cap_j / uncapped_energy_j
        else:
            scaling_factor = 1.0
        hexane_condensed_mol = hexane_capacity_mol * scaling_factor
        water_condensed_mol = water_capacity_mol * scaling_factor
        if water_condensed_mol > 0.0 and not (
            decoded.surface.free_water_kg > 0.0 or decoded.surface.attached_hexane_kg > 0.0
        ):
            # The level-0b routing: a free film on a hexane-wet particle, the
            # sorbed store through a dry shell.  The PD packet step books
            # condensate into the free water only, so a film-free packet dry of
            # free water has no carried route at level 0.
            raise PDSteamCondensingLevel0Refusal(
                f"steam condensing level 0: packet {packet.packet_id} is film-free and "
                "dry of free water, so its condensate belongs in the sorbed store "
                "through its dry shell, which the PD packet step does not carry at "
                "level 0; refuse, do not clamp"
            )
        rows.append(
            PDPacketCondensationRow(
                **common,
                stop_temperature_k=stop_temperature_k,
                stop_component=stop_component,
                equilibrium_energy_cap_j=equilibrium_energy_cap_j,
                equilibrium_energy_scaling_factor=scaling_factor,
                hexane_condensed_mol=hexane_condensed_mol,
                water_condensed_mol=water_condensed_mol,
                hexane_condensed_kg=hexane_condensed_mol
                * model.properties.hexane_molar_mass_kg_mol,
                water_condensed_kg=water_condensed_mol * model.properties.water_molar_mass_kg_mol,
                condensation_common_energy_j=math.fsum(
                    (
                        hexane_condensed_mol * hexane_energy_per_mole,
                        water_condensed_mol * water_energy_per_mole,
                    )
                ),
            )
        )

    packet_rows = tuple(rows)
    hexane_condensed_mol = math.fsum(row.hexane_condensed_mol for row in packet_rows)
    water_condensed_mol = math.fsum(row.water_condensed_mol for row in packet_rows)
    if hexane_condensed_mol <= 0.0 and water_condensed_mol <= 0.0:
        return None
    if hexane_condensed_mol > 0.0 and hexane_condensed_mol >= hexane_mol:
        raise PDKDomeBypassStepError(
            "steam condensing level 0 would remove the whole n-hexane inventory of the "
            "bypass stream; refuse, do not clamp"
        )
    if water_condensed_mol > 0.0 and water_condensed_mol >= water_mol:
        raise PDKDomeBypassStepError(
            "steam condensing level 0 would remove the whole water inventory of the "
            "bypass stream; refuse, do not clamp"
        )
    hexane_condensed_kg = math.fsum(row.hexane_condensed_kg for row in packet_rows)
    water_condensed_kg = math.fsum(row.water_condensed_kg for row in packet_rows)
    hexane_common_energy_j = math.fsum(
        row.hexane_condensed_mol * hexane_energy_per_mole for row in packet_rows
    )
    water_common_energy_j = math.fsum(
        row.water_condensed_mol * water_energy_per_mole for row in packet_rows
    )
    condensation_common_energy_j = math.fsum(
        row.condensation_common_energy_j for row in packet_rows
    )
    uncapped_total_energy_j = math.fsum(row.uncapped_common_energy_j for row in packet_rows)
    scaling_factor = (
        condensation_common_energy_j / uncapped_total_energy_j
        if uncapped_total_energy_j > 0.0
        else 1.0
    )
    capped_rows = tuple(row for row in packet_rows if math.isfinite(row.equilibrium_energy_cap_j))
    equilibrium_energy_cap_j = math.fsum(row.equilibrium_energy_cap_j for row in capped_rows)
    stop_temperature_k, stop_component = min(
        (row.stop_temperature_k, row.stop_component) for row in packet_rows if row.stop_component
    )
    remaining_total_mol = math.fsum((total_mol, -hexane_condensed_mol, -water_condensed_mol))
    remaining_water_mol = math.fsum((water_mol, -water_condensed_mol))
    residual_supersaturation = math.fsum(
        (
            stream.pressure_pa * remaining_water_mol / remaining_total_mol,
            -saturation_pressure,
        )
    )
    hexane_bases = tuple(
        dict.fromkeys(
            row.equilibrium_ratio_basis_id.split("/")[0].removeprefix("hexane:")
            for row in packet_rows
        )
    )
    water_bases = tuple(
        dict.fromkeys(
            row.equilibrium_ratio_basis_id.split("/")[1].removeprefix("water:")
            for row in packet_rows
        )
    )
    surface_temperatures = tuple(row.surface_temperature_k for row in packet_rows)
    return PDContactCondensationLedger(
        contact_fraction=math.nan,
        contact_area_m2=math.nan,
        contact_area_basis_id=PD_STEAM_CONDENSING_LEVEL0_RECEIVER_LAW_ID,
        rate_law_id=PD_STEAM_CONDENSING_LEVEL0_LAW_ID,
        energy_routing_id=PD_CONDENSATION_ENERGY_ROUTING_ID,
        receiver_law_id=PD_STEAM_CONDENSING_LEVEL0_RECEIVER_LAW_ID,
        cap_law_id=PD_CONDENSATION_CAP_LAW_ID,
        transfer_coefficient_law_id=PD_STEAM_CONDENSING_LEVEL0_CAPACITY_LAW_ID,
        equilibrium_ratio_law_id=PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID,
        equilibrium_ratio_basis_id=(
            "hexane:"
            + "+".join(hexane_bases)
            + "/water:"
            + "+".join(water_bases)
            + "/"
            + PD_CONDENSATION_EQUILIBRIUM_RATIO_LAW_ID
        ),
        duration_s=duration_s,
        stream_temperature_k=stream.temperature_k,
        stream_pressure_pa=stream.pressure_pa,
        surface_temperature_k=_representative_weighted_mean(
            surface_temperatures, counts, total_count
        ),
        packet_rows=packet_rows,
        gas_film_coefficient_w_m2_k=math.nan,
        molar_heat_capacity_j_mol_k=math.nan,
        molar_transfer_conductance_mol_s_m2_pa=math.nan,
        hexane_bulk_fugacity_pa=hexane_bulk_fugacity,
        water_bulk_fugacity_pa=water_bulk_fugacity,
        hexane_surface_liquid_fugacity_pa=_representative_weighted_mean(
            tuple(row.hexane_surface_liquid_fugacity_pa for row in packet_rows),
            counts,
            total_count,
        ),
        water_surface_liquid_fugacity_pa=_representative_weighted_mean(
            tuple(row.water_surface_liquid_fugacity_pa for row in packet_rows),
            counts,
            total_count,
        ),
        hexane_equilibrium_ratio=_representative_weighted_mean(
            tuple(row.hexane_equilibrium_ratio for row in packet_rows), counts, total_count
        ),
        water_equilibrium_ratio=_representative_weighted_mean(
            tuple(row.water_equilibrium_ratio for row in packet_rows), counts, total_count
        ),
        hexane_driving_force_pa=_representative_weighted_mean(
            tuple(row.hexane_driving_force_pa for row in packet_rows), counts, total_count
        ),
        water_driving_force_pa=_representative_weighted_mean(
            tuple(row.water_driving_force_pa for row in packet_rows), counts, total_count
        ),
        hexane_dew_point_k=_finite_weighted_mean(
            tuple(row.hexane_stop_temperature_k for row in packet_rows), counts
        ),
        water_dew_point_k=_finite_weighted_mean(
            tuple(row.water_stop_temperature_k for row in packet_rows), counts
        ),
        condensation_equilibrium_temperature_k=stop_temperature_k,
        stop_component=stop_component,
        hexane_contact_limited_capacity_mol=math.fsum(
            row.hexane_contact_limited_capacity_mol for row in packet_rows
        ),
        contact_limited_capacity_mol=math.fsum(
            row.water_contact_limited_capacity_mol for row in packet_rows
        ),
        upper_layer_declared_energy_j=math.fsum(
            row.declared_energy_share_j for row in packet_rows
        ),
        equilibrium_energy_cap_j=equilibrium_energy_cap_j,
        equilibrium_energy_scaling_factor=scaling_factor,
        hexane_condensed_mol=hexane_condensed_mol,
        hexane_condensed_kg=hexane_condensed_kg,
        condensed_water_mol=water_condensed_mol,
        condensed_water_kg=water_condensed_kg,
        hexane_common_datum_liquid_enthalpy_j_mol=hexane_liquid_enthalpy,
        hexane_seam_latent_j_mol=hexane_seam_latent,
        common_datum_liquid_enthalpy_j_mol=water_liquid_enthalpy,
        seam_latent_j_mol=water_seam_latent,
        hexane_condensation_common_energy_j=hexane_common_energy_j,
        water_condensation_common_energy_j=water_common_energy_j,
        condensation_common_energy_j=condensation_common_energy_j,
        water_partial_pressure_pa=partial_pressure,
        water_saturation_pressure_pa=saturation_pressure,
        partial_pressure_excess_pa=excess,
        equilibrium_relief_mol=relief_mol,
        residual_supersaturation_pa=residual_supersaturation,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDPacketHeatLedger:
    packet_id: str
    vertical_layer: int
    temperature_before_k: float
    temperature_after_k: float
    saturation_temperature_k: float
    declared_wall_heat_j: float
    declared_inter_k_heat_j: float
    declared_mechanical_work_j: float
    total_declared_heat_j: float
    total_declared_packet_energy_input_j: float
    thermal_only_temperature_after_k: float
    thermal_only_generated_hexane_kg: float
    external_water_before_kg: float
    external_water_after_kg: float
    external_water_residual_kg: float
    generated_hexane_kg: float
    #: F4 (2026-09-02): external free water co-boiled at the certified
    #: heteroazeotrope composition.  Exactly 0.0 on the single-liquid branch.
    generated_water_kg: float
    #: Common-datum energy of the WHOLE generated vapour (hexane plus its
    #: certified co-boiling water share), not of the hexane alone.
    generated_hexane_common_energy_j: float
    packet_common_energy_change_j: float
    mass_residual_kg: float
    common_energy_residual_j: float
    mass_binary64_roundoff_bound_kg: float
    common_energy_binary64_roundoff_bound_j: float
    maximum_relative_residual: float
    relative_limit: float
    #: WS-A (2026-09-03): the contact-limited condensate this packet received
    #: from the dome vapour, and the common-datum energy it carried in through
    #: the F2 seam.  Exactly 0.0 at the declared default contact fraction, and
    #: therefore on every production step while Q-PD-4 stands.
    condensed_water_kg: float = 0.0
    condensation_common_energy_j: float = 0.0
    #: CS-1a (2026-09-04): the n-hexane condensate this packet received from
    #: the same cold-surface term.  Exactly 0.0 at the declared default.
    condensed_hexane_kg: float = 0.0

    physically_qualifying: ClassVar[bool] = False
    external_water_inventory_thermally_coupled: ClassVar[bool] = True
    nonvolatile_external_water_limit: ClassVar[bool] = True
    external_water_transport_implemented: ClassVar[bool] = False
    heat_controlled_pure_hexane_limit: ClassVar[bool] = True
    mechanical_work_source_temperature_defined: ClassVar[bool] = False
    candidate_solved_from_t_n_combined_energy: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        return self.maximum_relative_residual <= self.relative_limit


@dataclass(frozen=True, slots=True, kw_only=True)
class PDInterKFaceLedger:
    lower_vertical_layer: int
    upper_vertical_layer: int
    effective_conductivity_w_m_k: float
    face_area_m2: float
    gradient_distance_m: float
    conductance_w_k: float
    lower_temperature_before_k: float
    upper_temperature_before_k: float
    lower_temperature_after_k: float
    upper_temperature_after_k: float
    lower_thermal_only_temperature_after_k: float
    upper_thermal_only_temperature_after_k: float
    heat_lower_to_upper_j: float
    lower_packet_allocation_j: float
    upper_packet_allocation_j: float
    allocation_residual_j: float
    entropy_production_j_k: float
    no_explicit_equilibrium_crossing: bool
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    internal_energy_transfer_only: ClassVar[bool] = True
    lambda_mix_used: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.maximum_relative_residual <= self.relative_limit
            and self.entropy_production_j_k >= 0.0
            and self.no_explicit_equilibrium_crossing
        )


def _layer_packets(
    tray: tray_host.AcceptedKCellTrayState,
    layer: int,
) -> tuple[tray_host.AcceptedTrayPacket, ...]:
    packets = tuple(
        sorted(
            (
                packet
                for packet in tray.packets
                if packet.owner_key.vertical_layer_id.value == layer
            ),
            key=lambda packet: packet.packet_id,
        )
    )
    if not packets:
        raise PDKDomeBypassStepError(f"PD:K{layer} has no weighted packets")
    return packets


def _layer_temperature_k(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    counts = tuple(packet.entry.weight.representative_particles for packet in packets)
    total = math.fsum(counts)
    anchor = min(codec.temperature_from_payload(packet.entry.payload_bytes) for packet in packets)
    return (
        anchor
        + math.fsum(
            count * (codec.temperature_from_payload(packet.entry.payload_bytes) - anchor)
            for packet, count in zip(packets, counts, strict=True)
        )
        / total
    )


def _layer_temperature_after_k(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    entries: dict[str, tray_host.PacketPlanEntry],
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    counts = tuple(packet.entry.weight.representative_particles for packet in packets)
    total = math.fsum(counts)
    anchor = codec.temperature_from_payload(entries[packets[0].packet_id].payload_bytes)
    return (
        anchor
        + math.fsum(
            count
            * (codec.temperature_from_payload(entries[packet.packet_id].payload_bytes) - anchor)
            for packet, count in zip(packets, counts, strict=True)
        )
        / total
    )


def _explicit_pair_crossed_equilibrium(
    *,
    first_before_k: float,
    second_before_k: float,
    first_after_k: float,
    second_after_k: float,
) -> bool:
    before_difference = first_before_k - second_before_k
    after_difference = first_after_k - second_after_k
    tolerance = 8.0 * _ulp_sum(
        first_before_k,
        second_before_k,
        first_after_k,
        second_after_k,
        before_difference,
        after_difference,
    )
    return (before_difference > tolerance and after_difference < -tolerance) or (
        before_difference < -tolerance and after_difference > tolerance
    )


def _packet_energy_at_temperature_j(
    packet: tray_host.AcceptedTrayPacket,
    temperature_k: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    complete = codec.decode(packet.entry.payload_bytes)
    inventory = complete.inventory
    dry = inventory.dry_matter_kg
    total_hexane = math.fsum((inventory.attached_hexane_kg, inventory.internal_hexane_kg))
    if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec):
        # D12 level 0: a film-free packet's layer energy at the adjacent
        # layer's temperature is issued on the variant that holds it there.
        rebuilt = _level0_initialize_at_temperature(
            codec,
            hexane_loading=total_hexane / dry,
            temperature_k=temperature_k,
            external_free_water_loading=inventory.external_water_kg / dry,
            retained_water_loading=inventory.retained_water_kg / dry,
            oil_fraction=inventory.residual_oil_label_kg / dry,
        )
        return (
            rebuilt.inventory.common_datum_energy_j * packet.entry.weight.representative_particles
        )
    rebuilt = codec.initialize(
        total_hexane_loading=total_hexane / dry,
        temperature_k=temperature_k,
        external_free_water_loading=inventory.external_water_kg / dry,
        retained_water_loading=inventory.retained_water_kg / dry,
        oil_fraction=inventory.residual_oil_label_kg / dry,
    )
    return rebuilt.inventory.common_datum_energy_j * packet.entry.weight.representative_particles


def _layer_energy_to_temperature_j(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    temperature_k: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    return math.fsum(
        _packet_energy_at_temperature_j(packet, temperature_k, codec=codec)
        - packet.weighted_totals.totals[6]
        for packet in packets
    )


def _layer_uniform_energy_to_temperature_j(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    temperature_k: float,
    *,
    self_temperature_k: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    """Energy that moves the layer, taken as uniform at its own mean, to ``temperature_k``.

    The uniform-layer room of ruling-packet option (iv): sum_i [E_i(T) - E_i(T_self)],
    the through-bed kernel's construction; the layer's own spread never enters.
    """

    return math.fsum(
        _packet_energy_at_temperature_j(packet, temperature_k, codec=codec)
        - _packet_energy_at_temperature_j(packet, self_temperature_k, codec=codec)
        for packet in packets
    )


def _layer_uniform_temperature_after_face_heat_k(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    temperature_before_k: float,
    face_heat_j: float,
    global_minimum_temperature_k: float,
    global_maximum_temperature_k: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    """The layer, uniform at its own mean, after the inter-K face heat ALONE.

    The pinned through-bed kernel's ``_temperature_after_internal_heat``: the
    bracket is the old layer-temperature envelope on the side the heat drives,
    the target outside it refuses ("would cross the old layer equilibrium"),
    and the bisection runs to binary64 adjacency within the kernel's 120 halvings.
    The energy is the uniform-layer room sum_i [E_i(T) - E_i(T_before)].
    """

    if face_heat_j == 0.0:
        return temperature_before_k
    if face_heat_j > 0.0:
        lower = temperature_before_k
        upper = global_maximum_temperature_k
    else:
        lower = global_minimum_temperature_k
        upper = temperature_before_k
    if lower == upper:
        raise PDKDomeBypassStepError(
            "explicit inter-K conduction would leave the old temperature envelope; "
            "reject and retry at a shorter endpoint"
        )
    # R3, the owner's B3 ruling of 2026-08-31, site 2 (study-lane knob, default OFF): with
    # ``DTDC_PD_INTERK_CAPACITY_ANALYTIC`` exactly "1" the room below is the analytic layer
    # heat capacity times (T - T_before) (``study_lane/pd_interk_capacity_analytic.py``),
    # same bracket, check, bisection and refusal text.  Off, the code below is the former
    # form, unchanged.
    if os.environ.get("DTDC_PD_INTERK_CAPACITY_ANALYTIC") == "1":
        from .study_lane import pd_interk_capacity_analytic as interk_analytic

        return interk_analytic.uniform_temperature_after_face_heat_k(
            packets,
            temperature_before_k=temperature_before_k,
            face_heat_j=face_heat_j,
            lower_k=lower,
            upper_k=upper,
            codec=codec,
        )
    energy_before = tuple(
        _packet_energy_at_temperature_j(packet, temperature_before_k, codec=codec)
        for packet in packets
    )

    def room(temperature_k: float) -> float:
        return math.fsum(
            _packet_energy_at_temperature_j(packet, temperature_k, codec=codec) - before
            for packet, before in zip(packets, energy_before, strict=True)
        )

    if not room(lower) <= face_heat_j <= room(upper):
        raise PDKDomeBypassStepError(
            "explicit inter-K conduction would cross the old layer equilibrium; "
            "reject and retry at a shorter endpoint"
        )
    for _ in range(120):
        middle = 0.5 * (lower + upper)
        if middle == lower or middle == upper:
            break
        if room(middle) < face_heat_j:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def _conduction_only_order_crossed(
    request: "PDKTrayStepRequest",
    bottom_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    upper_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    lower_temperature_k: float,
    upper_temperature_k: float,
    face_heat_j: float,
    all_sources_lower_after_k: float,
    all_sources_upper_after_k: float,
    all_sources_crossed: bool,
    saturation_temperature_k: float,
    zero_face_heat_entries: Callable[[], dict[str, tray_host.PacketPlanEntry]],
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> bool:
    """The order guard's verdict under ``DTDC_PD_INTERK_ORDER_CONDUCTION_ONLY``.

    The uniform layers after the face heat alone, checked with the shipped
    tolerance; a refusal of the bisection or of a caloric seed raises.  With
    ``PD_INTERK_ORDER_LOCK_OBSERVER`` set, one report-only row per evaluation
    carries both readings; nothing in the report raises into the step.
    """

    observer = PD_INTERK_ORDER_LOCK_OBSERVER
    row: dict[str, object] = {}
    if observer is not None:
        before_difference = lower_temperature_k - upper_temperature_k
        after_difference = all_sources_lower_after_k - all_sources_upper_after_k
        row = {
            "guard": "order_conduction_only",
            "tray_id": request.prior.tray.physical_tray_id,
            "t_start_s": request.prior.tray.time_s,
            "t_end_s": request.end_time_s,
            "duration_s": request.duration_s,
            "form_in_force": "conduction_only",
            "lock_k": saturation_temperature_k,
            "face_heat_j": face_heat_j,
            "lower_temperature_k": lower_temperature_k,
            "upper_temperature_k": upper_temperature_k,
            "before_difference_k": before_difference,
            "all_sources_after_difference_k": after_difference,
            "all_sources_crossed": all_sources_crossed,
            "all_sources_tolerance_k": 8.0
            * _ulp_sum(
                lower_temperature_k,
                upper_temperature_k,
                all_sources_lower_after_k,
                all_sources_upper_after_k,
                before_difference,
                after_difference,
            ),
            "wrong_way_heat_j": (
                abs(face_heat_j)
                * abs(after_difference)
                / (abs(after_difference) + abs(before_difference))
                if all_sources_crossed
                else 0.0
            ),
        }

    def report(extra: dict[str, object]) -> None:
        if observer is None:
            return
        try:
            row.update(extra)
            observer(MappingProxyType(row))
        except Exception:  # noqa: BLE001 - report only
            return

    minimum = min(lower_temperature_k, upper_temperature_k)
    maximum = max(lower_temperature_k, upper_temperature_k)
    try:
        lower_after = _layer_uniform_temperature_after_face_heat_k(
            bottom_packets,
            temperature_before_k=lower_temperature_k,
            face_heat_j=-face_heat_j,
            global_minimum_temperature_k=minimum,
            global_maximum_temperature_k=maximum,
            codec=codec,
        )
        upper_after = _layer_uniform_temperature_after_face_heat_k(
            upper_packets,
            temperature_before_k=upper_temperature_k,
            face_heat_j=face_heat_j,
            global_minimum_temperature_k=minimum,
            global_maximum_temperature_k=maximum,
            codec=codec,
        )
    except high_loading.HighLoadingPacketAdapterError as error:
        refusal = PDKDomeBypassStepError(
            "inter-K caloric guard left the accepted packet branch; "
            "reject and retry at a shorter endpoint"
        )
        report({"conduction_only_refusal": f"{refusal}: {type(error).__name__}: {error}"})
        raise refusal from error
    except PDKDomeBypassStepError as refusal:
        report({"conduction_only_refusal": str(refusal)})
        raise
    crossed = _explicit_pair_crossed_equilibrium(
        first_before_k=lower_temperature_k,
        second_before_k=upper_temperature_k,
        first_after_k=lower_after,
        second_after_k=upper_after,
    )
    if observer is None:
        return crossed
    extra: dict[str, object] = {
        "conduction_only_lower_after_k": lower_after,
        "conduction_only_upper_after_k": upper_after,
        "conduction_only_after_difference_k": lower_after - upper_after,
        "conduction_only_crossed": crossed,
    }
    if all_sources_crossed:
        try:
            counterfactual = zero_face_heat_entries()
            lower_zero = _layer_temperature_after_k(bottom_packets, counterfactual, codec=codec)
            upper_zero = _layer_temperature_after_k(upper_packets, counterfactual, codec=codec)
            extra.update(
                zero_face_heat_after_difference_k=lower_zero - upper_zero,
                zero_face_heat_crossed=_explicit_pair_crossed_equilibrium(
                    first_before_k=lower_temperature_k,
                    second_before_k=upper_temperature_k,
                    first_after_k=lower_zero,
                    second_after_k=upper_zero,
                ),
            )
        except Exception as error:  # noqa: BLE001 - report only
            extra["zero_face_heat_refusal"] = f"{type(error).__name__}: {error}"
    report(extra)
    return crossed


def _inter_k_crossing_capacities_j(
    bottom_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    upper_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    lower_temperature_k: float,
    upper_temperature_k: float,
    face_heat_j: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
    uniform_layer: bool,
) -> tuple[float, float]:
    """The lower and upper layers' crossing capacities in the form asked for.

    ``uniform_layer=False`` is the shipped form, the same expressions in the same
    order as before this helper existed; ``True`` is the uniform-layer form.
    """

    if uniform_layer:

        def lower_room() -> float:
            return _layer_uniform_energy_to_temperature_j(
                bottom_packets,
                upper_temperature_k,
                self_temperature_k=lower_temperature_k,
                codec=codec,
            )

        def upper_room() -> float:
            return _layer_uniform_energy_to_temperature_j(
                upper_packets,
                lower_temperature_k,
                self_temperature_k=upper_temperature_k,
                codec=codec,
            )

    else:

        def lower_room() -> float:
            return _layer_energy_to_temperature_j(
                bottom_packets,
                upper_temperature_k,
                codec=codec,
            )

        def upper_room() -> float:
            return _layer_energy_to_temperature_j(
                upper_packets,
                lower_temperature_k,
                codec=codec,
            )

    if face_heat_j > 0.0:
        lower_capacity = -lower_room()
        upper_capacity = upper_room()
    elif face_heat_j < 0.0:
        lower_capacity = lower_room()
        upper_capacity = -upper_room()
    else:
        lower_capacity = upper_capacity = 0.0
    return lower_capacity, upper_capacity


def _inter_k_crossing_verdict(
    face_heat_j: float,
    lower_capacity_j: float,
    upper_capacity_j: float,
) -> tuple[float, float, bool]:
    """The crossing capacity, its binary64 bound and whether the face heat exceeds it."""

    crossing_capacity = min(lower_capacity_j, upper_capacity_j)
    crossing_bound = _ulp_sum(
        face_heat_j,
        lower_capacity_j,
        upper_capacity_j,
        crossing_capacity,
    )
    return (
        crossing_capacity,
        crossing_bound,
        abs(face_heat_j) > crossing_capacity + crossing_bound,
    )


def _inter_k_crossing_refused(
    request: object,
    bottom_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    upper_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    lower_temperature_k: float,
    upper_temperature_k: float,
    face_heat_j: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> bool:
    """The crossing guard's verdict under the knob, with both forms reported.

    Unarmed and with no observer this is the shipped form and nothing else.
    With an observer, both forms are evaluated and handed to it; when the
    knob is off, the uniform-layer evaluation is report only (a refusal of
    its caloric seeds is reported, never raised).
    """

    lower_capacity, upper_capacity = _inter_k_crossing_capacities_j(
        bottom_packets,
        upper_packets,
        lower_temperature_k=lower_temperature_k,
        upper_temperature_k=upper_temperature_k,
        face_heat_j=face_heat_j,
        codec=codec,
        uniform_layer=False,
    )
    crossing_capacity, crossing_bound, shipped_refused = _inter_k_crossing_verdict(
        face_heat_j, lower_capacity, upper_capacity
    )
    armed = PD_INTERK_CAPACITY_UNIFORM_LAYER_ARMED
    observer = PD_INTERK_CAPACITY_OBSERVER
    if not armed and observer is None:
        return shipped_refused
    uniform: dict[str, object] = {}
    uniform_refused = shipped_refused
    try:
        uniform_lower, uniform_upper = _inter_k_crossing_capacities_j(
            bottom_packets,
            upper_packets,
            lower_temperature_k=lower_temperature_k,
            upper_temperature_k=upper_temperature_k,
            face_heat_j=face_heat_j,
            codec=codec,
            uniform_layer=True,
        )
    except Exception as refusal:  # report only when unarmed; armed, it propagates
        if armed:
            raise
        uniform = {"error": f"{type(refusal).__name__}: {refusal}"}
    else:
        uniform_capacity, uniform_bound, uniform_refused = _inter_k_crossing_verdict(
            face_heat_j, uniform_lower, uniform_upper
        )
        uniform = {
            "lower_capacity_j": uniform_lower,
            "upper_capacity_j": uniform_upper,
            "crossing_capacity_j": uniform_capacity,
            "crossing_bound_j": uniform_bound,
            "refuses": uniform_refused,
        }
    if observer is not None:
        prior_tray = getattr(getattr(request, "prior", None), "tray", None)
        observer(
            MappingProxyType(
                {
                    "tray": getattr(prior_tray, "physical_tray_id", None),
                    "t_start_s": getattr(prior_tray, "time_s", None),
                    "t_end_s": getattr(request, "end_time_s", None),
                    "lower_temperature_k": lower_temperature_k,
                    "upper_temperature_k": upper_temperature_k,
                    "face_heat_j": face_heat_j,
                    "form_in_force": "uniform_layer" if armed else "shipped",
                    "shipped": {
                        "lower_capacity_j": lower_capacity,
                        "upper_capacity_j": upper_capacity,
                        "crossing_capacity_j": crossing_capacity,
                        "crossing_bound_j": crossing_bound,
                        "refuses": shipped_refused,
                    },
                    "uniform_layer": uniform,
                }
            )
        )
    return uniform_refused if armed else shipped_refused


def _report_order_lock_forms(
    request: "PDKTrayStepRequest",
    bottom_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    upper_packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    thermal_guard_entries: dict[str, tray_host.PacketPlanEntry],
    alternative_entries: dict[str, tray_host.PacketPlanEntry],
    alternative_refusals: list[str],
    lower_temperature_k: float,
    upper_temperature_k: float,
    lower_after_in_force_k: float,
    upper_after_in_force_k: float,
    crossed_in_force: bool,
    saturation_temperature_k: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> None:
    """Report only: both lock rules' order-guard verdicts to the order-lock observer.

    The form in force is ``PD_INTERK_ORDER_LOCK_CONSISTENT_ARMED``'s; the other
    form's thermal-only entries replace the in-force ones for the work-bearing
    packets only (every other packet's thermal-only entry is its full entry in
    both forms).  Nothing here raises into the step.
    """

    observer = PD_INTERK_ORDER_LOCK_OBSERVER
    if observer is None:
        return
    try:
        row: dict[str, object] = {
            "guard": "order_lock",
            "tray_id": request.prior.tray.physical_tray_id,
            "t_start_s": request.prior.tray.time_s,
            "t_end_s": request.end_time_s,
            "duration_s": request.duration_s,
            "form_in_force": ("consistent" if PD_INTERK_ORDER_LOCK_CONSISTENT_ARMED else "shipped"),
            "lock_k": saturation_temperature_k,
            "before_difference_k": lower_temperature_k - upper_temperature_k,
            "in_force_after_difference_k": lower_after_in_force_k - upper_after_in_force_k,
            "in_force_crossed": crossed_in_force,
            "work_bearing_packets": len(alternative_entries),
            "alternative_refusals": list(alternative_refusals),
        }
        if not alternative_refusals:
            other = dict(thermal_guard_entries)
            other.update(alternative_entries)
            lower_other = _layer_temperature_after_k(bottom_packets, other, codec=codec)
            upper_other = _layer_temperature_after_k(upper_packets, other, codec=codec)
            row.update(
                other_after_difference_k=lower_other - upper_other,
                other_crossed=_explicit_pair_crossed_equilibrium(
                    first_before_k=lower_temperature_k,
                    second_before_k=upper_temperature_k,
                    first_after_k=lower_other,
                    second_after_k=upper_other,
                ),
                other_upper_after_minus_lock_k=upper_other - saturation_temperature_k,
                other_lower_after_minus_lock_k=lower_other - saturation_temperature_k,
            )
        row["in_force_upper_after_minus_lock_k"] = upper_after_in_force_k - saturation_temperature_k
        row["in_force_lower_after_minus_lock_k"] = lower_after_in_force_k - saturation_temperature_k
        observer(MappingProxyType(row))
    except Exception:  # noqa: BLE001 - report only
        return


def _allocate_layer_heat(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    total_heat_j: float,
) -> dict[str, float]:
    total_count = math.fsum(packet.entry.weight.representative_particles for packet in packets)
    allocated: dict[str, float] = {}
    partial: list[float] = []
    for packet in packets[:-1]:
        value = total_heat_j * packet.entry.weight.representative_particles / total_count
        allocated[packet.packet_id] = value
        partial.append(value)
    allocated[packets[-1].packet_id] = math.fsum((total_heat_j, -math.fsum(partial)))
    return allocated


def _level0_is_dual(codec: object) -> bool:
    return type(codec) is dual_variant.EngineeringDualVariantPacketAdapter


def _level0_initialize_at_temperature(
    codec: dual_variant.EngineeringDualVariantPacketAdapter,
    *,
    hexane_loading: float,
    temperature_k: float,
    external_free_water_loading: float,
    retained_water_loading: float,
    oil_fraction: float,
):
    """D12 level 0: an energy seed at ``temperature_k`` on whichever variant holds it.

    The tag-1 (film) initializer first, exactly as the unarmed step calls
    it; when it refuses the loading at ``temperature_k`` (a subcritical feed,
    or an exhausted attached film), the same loading is issued under the
    tag-2 initializer.  A loading neither codec admits refuses typed: at the
    saturated state that is the tag-2 floor.
    """

    common = dict(
        temperature_k=temperature_k,
        external_free_water_loading=external_free_water_loading,
        retained_water_loading=retained_water_loading,
        oil_fraction=oil_fraction,
    )
    try:
        return codec.initialize(total_hexane_loading=hexane_loading, **common)
    except high_loading.HighLoadingPacketStateError as film_refusal:
        try:
            return codec.initialize(internal_hexane_loading=hexane_loading, **common)
        except subcritical.SubcriticalPacketAdapterError as error:
            raise PDFlashLevel0FloorRefusal(
                f"neither the tag-1 nor the tag-2 codec admits loading {hexane_loading!r} "
                f"kg/kg at {temperature_k!r} K (tag-1: {film_refusal}; tag-2: {error}); "
                "at the saturated state the flash would pass the tag-2 floor, and tag-3 "
                "on a PD tray is outside level 0 - refuse, do not clamp"
            ) from error


def _level0_bracketed_root(residual, low, high, low_pair, high_pair):
    """D12 level 0: Illinois false position on the boiling residual's bracket.

    The residual is monotone in the removal (non-positive at ``low``,
    positive at ``high``), so false position with the Illinois weight
    halving converges superlinearly where plain bisection spends up to 110
    codec inversions.  It never leaves the bracket, stops at adjacent floats
    or an exact zero, and hands whatever bracket remains to the shipped
    bisection, which finishes it to adjacent floats.  The root criterion and
    the chosen end (the smaller absolute residual) are unchanged.
    """

    f_low = low_pair[0]
    f_high = high_pair[0]
    side = 0
    for _ in range(110):
        if f_low == 0.0 or f_high == 0.0 or math.nextafter(low, high) == high:
            break
        point = high - f_high * (high - low) / (f_high - f_low)
        if not low < point < high:
            point = 0.5 * (low + high)
            if point in (low, high):
                break
        pair = residual(point)
        if pair[0] <= 0.0:
            low, low_pair, f_low = point, pair, pair[0]
            if side == -1:
                f_high *= 0.5
            side = -1
        else:
            high, high_pair, f_high = point, pair, pair[0]
            if side == 1:
                f_low *= 0.5
            side = 1
    return low, high, low_pair, high_pair


def _level0_packet_at_saturation(
    before,
    *,
    total_hexane_mass_kg: float,
    saturation_temperature_k: float,
    codec: dual_variant.EngineeringDualVariantPacketAdapter,
    external_free_water_mass_kg: float,
):
    """D12 level 0 (design site 3): the saturated state across the film seam.

    The seed supplies the saturated energy on the variant that holds the
    loading at the saturation temperature, and ``codec.advance`` lets the
    dual adapter carry the packet's own conserved totals across the tag-1 to
    tag-2 seam (``engineering_dual_variant_packet_adapter.py``, the
    conserved-totals reconstitution).  A landing on any variant other than
    tag-1 or tag-2 is the floor, refused typed.
    """

    dry = before.state.dry_meal_mass_kg
    seed = _level0_initialize_at_temperature(
        codec,
        hexane_loading=total_hexane_mass_kg / dry,
        temperature_k=saturation_temperature_k,
        external_free_water_loading=external_free_water_mass_kg / dry,
        retained_water_loading=before.state.retained_water_loading,
        oil_fraction=before.state.oil_fraction,
    )
    after = codec.advance(
        before.payload_bytes,
        total_hexane_mass_after_kg=total_hexane_mass_kg,
        external_free_water_mass_after_kg=external_free_water_mass_kg,
        envelope_energy_after_j=seed.inventory.common_datum_energy_j,
    )
    if (
        type(after) is not high_loading.CompleteHighLoadingPacket
        and type(after) is not subcritical.CompleteSubcriticalPacket
    ):
        raise PDFlashLevel0FloorRefusal(
            f"the saturated PD state landed on {type(after).__name__}, past the tag-2 "
            "floor; tag-3 on a PD tray is outside level 0 - refuse, do not clamp"
        )
    return after


def _packet_at_saturation(
    before: high_loading.CompleteHighLoadingPacket,
    *,
    total_hexane_mass_kg: float,
    saturation_temperature_k: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
    external_free_water_mass_kg: float | None = None,
) -> high_loading.CompleteHighLoadingPacket:
    external_free_water = (
        before.inventory.external_water_kg
        if external_free_water_mass_kg is None
        else external_free_water_mass_kg
    )
    if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec):
        return _level0_packet_at_saturation(
            before,
            total_hexane_mass_kg=total_hexane_mass_kg,
            saturation_temperature_k=saturation_temperature_k,
            codec=codec,
            external_free_water_mass_kg=external_free_water,
        )
    seed = codec.initialize(
        total_hexane_loading=total_hexane_mass_kg / before.state.dry_meal_mass_kg,
        temperature_k=saturation_temperature_k,
        external_free_water_loading=(external_free_water / before.state.dry_meal_mass_kg),
        retained_water_loading=before.state.retained_water_loading,
        oil_fraction=before.state.oil_fraction,
    )
    return codec.advance(
        before.payload_bytes,
        total_hexane_mass_after_kg=total_hexane_mass_kg,
        external_free_water_mass_after_kg=external_free_water,
        envelope_energy_after_j=seed.inventory.common_datum_energy_j,
    )


def _advance_one_packet(
    packet: tray_host.AcceptedTrayPacket,
    *,
    declared_wall_heat_j: float,
    declared_inter_k_heat_j: float,
    declared_mechanical_work_j: float,
    old_wall_temperature_k: float,
    adjacent_old_temperature_k: float,
    saturation_temperature_k: float,
    tray_pressure_pa: float,
    properties: cc.DeclaredGasProperties | native_caloric.NativeGasProperties,
    codec: high_loading.HighLoadingV1PacketAdapter,
    relative_limit: float,
    co_boiling: PDCoBoilingState | None = None,
    declared_condensate_water_kg: float = 0.0,
    declared_condensate_hexane_kg: float = 0.0,
    declared_condensation_energy_j: float = 0.0,
    declared_condensation_equilibrium_temperature_k: float | None = None,
    thermal_guard_temperature_after_k: float | None = None,
    thermal_guard_generated_hexane_kg: float | None = None,
    lock_band_boils: bool = False,
    exploration_call_kind: str | None = None,
) -> tuple[tray_host.PacketPlanEntry, PDPacketHeatLedger]:
    before = codec.decode(packet.entry.payload_bytes)
    inventory_before = before.inventory
    if inventory_before.external_water_kg != before.surface.free_water_kg:
        raise PDKDomeBypassStepError(
            f"packet {packet.packet_id} has inconsistent external-water authorities"
        )
    if before.state.attached_hexane_mass_kg <= 0.0 and not (
        # D12 level 0 (design site 2): a film-free packet under the dual
        # codec is admitted on the armed flash branch; tag-3 still refuses.
        PD_FLASH_LEVEL0_ARMED
        and _level0_is_dual(codec)
        and type(before) is subcritical.CompleteSubcriticalPacket
    ):
        raise PDKDomeBypassStepError(
            f"packet {packet.packet_id} is outside the admitted attached-hexane tag-1 branch"
        )
    representatives = packet.entry.weight.representative_particles
    # WS-A: this packet's share of the contact-limited dome condensate, and
    # the common-datum energy that condensate carried in across the F2 seam.
    # Both are exactly 0.0 at the declared default contact fraction, and every
    # expression below then reduces to its pre-WS-A form bit for bit
    # (``math.fsum`` sums exactly, so an exact 0.0 addend never moves a sum).
    _require_float(
        "declared packet condensate water",
        declared_condensate_water_kg,
        nonnegative=True,
    )
    _require_float(
        "declared packet condensate hexane",
        declared_condensate_hexane_kg,
        nonnegative=True,
    )
    _require_float(
        "declared packet condensation energy",
        declared_condensation_energy_j,
        nonnegative=True,
    )
    condensate_per_particle = declared_condensate_water_kg / representatives
    condensate_hexane_per_particle = declared_condensate_hexane_kg / representatives
    condensation_energy_per_particle = declared_condensation_energy_j / representatives
    if inventory_before.external_water_kg > 0.0 or condensate_per_particle > 0.0:
        water_saturation_temperature = water_properties.saturation_temperature(codec.pressure_pa)
        water_tolerance = 5.0e-11 * max(1.0, water_saturation_temperature)
        if inventory_before.external_water_kg > 0.0 and (
            before.state.temperature_k >= water_saturation_temperature - water_tolerance
        ):
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} external water is not inside the declared "
                "subcooled nonvolatile-liquid branch"
            )
    total_declared_heat_j = math.fsum((declared_wall_heat_j, declared_inter_k_heat_j))
    _require_float(
        "declared packet mechanical work",
        declared_mechanical_work_j,
        nonnegative=True,
    )
    total_declared_packet_energy_input_j = math.fsum(
        (total_declared_heat_j, declared_mechanical_work_j)
    )
    heat_per_particle = total_declared_packet_energy_input_j / representatives
    # F-CS1-1 (ceremony S:6.3), owner-ruled correction (A), 2026-09-05.  The
    # DECLARED HEAT alone, per particle, with the mechanical work excluded.
    # The condensation-free counterfactual below is compared against an
    # envelope of old DRIVING TEMPERATURES, which by construction knows
    # nothing of shaft work, so the counterfactual must be built on the same
    # shaft-work-free energy input the ordinary envelope check already uses
    # through ``thermal_guard_temperature_after_k``.  On a work-free packet
    # ``math.fsum((x, 0.0))`` is exactly ``x``, so this IS
    # ``heat_per_particle`` there, and the thermal-only guard solve - which is
    # always called with ``declared_mechanical_work_j=0.0`` - is unmoved bit
    # for bit.
    declared_heat_per_particle = total_declared_heat_j / representatives
    # WS-A: the condensation latent is an energy INPUT to the packet, but it is
    # not a declared HEAT - it never enters the shared wall's solid-heat row -
    # so it is carried separately and joins only the caloric target and the
    # packet's own common-datum energy row.
    energy_input_per_particle = math.fsum(
        (heat_per_particle, condensation_energy_per_particle)
    )
    external_water_after_condensation_kg = math.fsum(
        (inventory_before.external_water_kg, condensate_per_particle)
    )
    total_hexane_before = math.fsum(
        (inventory_before.attached_hexane_kg, inventory_before.internal_hexane_kg)
    )
    # CS-1a: the hexane condensate joins the packet own hexane inventory
    # before any boiling solve, exactly as the water condensate joins its
    # external free water.  Exactly 0.0 at the declared default, and
    # math.fsum sums exactly, so this IS total_hexane_before there.
    total_hexane_after_condensation = math.fsum(
        (total_hexane_before, condensate_hexane_per_particle)
    )
    # RULED 2026-09-06 (owner ruling on point A of the three-point packet):
    # the declared materiality floor goes under the step's own relative
    # temperature tolerance, so that EVERY consumer of it carries one floor -
    # the after-temperature driving-envelope guard, the condensation-free
    # counterfactual guard, the boiling-saturation interval check, and the
    # two saturation comparisons in this step's own branch selection.  Before
    # the ruling the relative term alone gave 1.709e-8 K at the fixture's
    # saturation temperature and the guard that walked the whole ladder had
    # no floor at all (F-GATE-2).  See PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K
    # for what the floor admits and the measurements that size it.
    #
    # The two saturation comparisons this tolerance also feeds (the
    # ``begins_above_saturation`` test on the next line and the
    # sensible-branch acceptance) and the boiling-saturation interval check
    # at the end of the step are widened by the same floor, by construction
    # and not by oversight.  Each is an INTERVAL test of a temperature
    # against a saturation temperature; the floor widens each interval by
    # 1.0e-6 K, which is 1e-5 of the ratified 0.1 K per-tray certificate
    # bound and six orders below the 5-7 K by which every shipped PD state
    # already begins above its heteroazeotropic boiling point (the F4
    # measurement below).  Nothing refused for a physical reason becomes
    # admissible at that width.
    temperature_tolerance = max(
        5.0e-11 * max(1.0, saturation_temperature_k),
        PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K,
    )
    begins_above_saturation = (
        before.state.temperature_k > saturation_temperature_k + temperature_tolerance
    )
    if begins_above_saturation:
        # F4 (2026-09-02) - RETIRED AND RESTATED, never weakened.  The
        # pre-F4 form refused outright, on the premise "a packet cannot
        # begin above its own boiling point".  That premise is preserved as
        # the OUTCOME rather than as an entry gate: with the boiling point
        # now taken from the CERTIFIED heteroazeotrope authority instead of
        # the pure n-hexane root, every shipped PD state begins 5-7 K above
        # it (measured: fixture packets at 340.0 K against T_het =
        # 334.4803981593741 K at 101 325 Pa), because it was equilibrated
        # against the wrong authority.  Such a state is SUPERHEATED against
        # the two-liquid equilibrium and physically flashes to it.  The
        # branch below already solves exactly that - it removes vapour until
        # the packet sits at the boiling temperature - so the flash needs no
        # new law, only this admission.  Everything else still fails closed:
        # a packet that cannot reach the boiling temperature is refused by
        # the "lacks the heat required to reach saturation" gate, exhausted
        # hexane hands off, and exhausted co-boiling water refuses typed.
        if co_boiling is None:
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} begins above pure-hexane saturation"
            )
    if (
        total_declared_packet_energy_input_j == 0.0
        and declared_condensate_water_kg == 0.0
        and declared_condensate_hexane_kg == 0.0
        and declared_condensation_energy_j == 0.0
    ):
        return packet.entry, PDPacketHeatLedger(
            packet_id=packet.packet_id,
            vertical_layer=packet.owner_key.vertical_layer_id.value,
            temperature_before_k=before.state.temperature_k,
            temperature_after_k=before.state.temperature_k,
            saturation_temperature_k=saturation_temperature_k,
            declared_wall_heat_j=declared_wall_heat_j,
            declared_inter_k_heat_j=declared_inter_k_heat_j,
            declared_mechanical_work_j=declared_mechanical_work_j,
            total_declared_heat_j=total_declared_heat_j,
            total_declared_packet_energy_input_j=(total_declared_packet_energy_input_j),
            thermal_only_temperature_after_k=before.state.temperature_k,
            thermal_only_generated_hexane_kg=0.0,
            external_water_before_kg=inventory_before.external_water_kg * representatives,
            external_water_after_kg=inventory_before.external_water_kg * representatives,
            external_water_residual_kg=0.0,
            generated_hexane_kg=0.0,
            generated_water_kg=0.0,
            generated_hexane_common_energy_j=0.0,
            packet_common_energy_change_j=0.0,
            mass_residual_kg=0.0,
            common_energy_residual_j=0.0,
            mass_binary64_roundoff_bound_kg=0.0,
            common_energy_binary64_roundoff_bound_j=0.0,
            maximum_relative_residual=0.0,
            relative_limit=relative_limit,
        )

    sensible_error: Exception | None = None
    try:
        sensible = codec.advance(
            before.payload_bytes,
            total_hexane_mass_after_kg=total_hexane_after_condensation,
            external_free_water_mass_after_kg=external_water_after_condensation_kg,
            envelope_energy_after_j=math.fsum(
                (inventory_before.common_datum_energy_j, energy_input_per_particle)
            ),
        )
    except (
        high_loading.HighLoadingCaloricBracketError,
        # D12 level 0: armed, a film-free packet's sensible solve past its own
        # caloric band is the boiling branch exactly as tag-1's is.
        *(_LEVEL0_SENSIBLE_ERRORS if PD_FLASH_LEVEL0_ARMED else ()),
    ) as error:
        sensible = None
        sensible_error = error
    tag2_leg_choice = None
    if PD_TAG2_LEG_ARMED and PD_TAG2_LEG_HOOK is not None:
        # Tag-2 crossing level 0 (Q2): a registered tag-2 packet's release is the
        # smaller of its dry-shell supply and the heat-limited flash below.
        tag2_leg_choice = PD_TAG2_LEG_HOOK.choose(
            packet=packet,
            before=before,
            codec=codec,
            properties=properties,
            co_boiling=co_boiling,
            saturation_temperature_k=saturation_temperature_k,
            tray_pressure_pa=tray_pressure_pa,
            temperature_tolerance=temperature_tolerance,
            sensible=sensible,
            lock_band_boils=lock_band_boils,
            energy_input_per_particle=energy_input_per_particle,
            available_energy_j=math.fsum(
                (inventory_before.common_datum_energy_j, energy_input_per_particle)
            ),
            total_hexane_after_condensation=total_hexane_after_condensation,
            external_water_after_condensation_kg=external_water_after_condensation_kg,
            condensate_hexane_per_particle=condensate_hexane_per_particle,
        )
    if tag2_leg_choice is not None and tag2_leg_choice.after is not None:
        after = tag2_leg_choice.after
        generated_per_particle = tag2_leg_choice.generated_per_particle
        generated_water_per_particle = tag2_leg_choice.generated_water_per_particle
        vapour_energy_per_particle = tag2_leg_choice.vapour_energy_per_particle
    elif (
        sensible is not None
        and (sensible.state.temperature_k <= saturation_temperature_k + temperature_tolerance)
        # The consistent lock rule (``DTDC_PD_INTERK_ORDER_LOCK_CONSISTENT``, the
        # shaft-work-free check only): a packet taking a POSITIVE energy input whose
        # sensible solve lands above the lock inside the band the branch selection
        # calls "at the lock" boils, as the full step does.  A packet whose input is
        # not positive cannot boil and keeps the sensible branch (measured K1c-OLC,
        # step 22).  ``lock_band_boils`` is False on every other call.
        and not (
            lock_band_boils
            and energy_input_per_particle > 0.0
            and sensible.state.temperature_k > saturation_temperature_k
        )
    ):
        after = sensible
        generated_per_particle = 0.0
        generated_water_per_particle = 0.0
        vapour_energy_per_particle = 0.0
    else:
        if energy_input_per_particle <= 0.0 and not begins_above_saturation:
            if sensible_error is not None:
                raise sensible_error
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} cannot evaporate hexane without positive packet "
                "energy input"
            )
        saturated_zero = _packet_at_saturation(
            before,
            total_hexane_mass_kg=total_hexane_after_condensation,
            saturation_temperature_k=saturation_temperature_k,
            codec=codec,
            external_free_water_mass_kg=external_water_after_condensation_kg,
        )
        # F4 (2026-09-02): the vapour that leaves is the CERTIFIED co-boiling
        # mixture, not pure n-hexane.  ``vapour_enthalpy`` is therefore the
        # common-datum enthalpy carried per mole of DEPARTING n-HEXANE - the
        # mixture enthalpy divided by the certified hexane mole fraction - and
        # every mole of hexane drags ``water_per_hexane_mol`` moles of water
        # out of the packet's external free-water inventory with it.  On the
        # pure-hexane branch (no external free water present) both reduce to
        # the pre-F4 arithmetic exactly.
        vapour_enthalpy = _vapour_enthalpy_per_hexane_mole_j(
            properties,
            co_boiling=co_boiling,
            saturation_temperature_k=saturation_temperature_k,
            pressure_pa=tray_pressure_pa,
            datum_adapter=codec.component_datum_adapter,
        )
        water_per_hexane_kg = _water_mass_per_hexane_mass(properties, co_boiling=co_boiling)
        available_energy = math.fsum(
            (inventory_before.common_datum_energy_j, energy_input_per_particle)
        )

        def residual(
            removed_hexane_kg: float,
        ) -> tuple[float, high_loading.CompleteHighLoadingPacket]:
            removed_water_kg = removed_hexane_kg * water_per_hexane_kg
            water_after = math.fsum(
                (external_water_after_condensation_kg, -removed_water_kg)
            )
            if water_after < 0.0:
                raise PDCoBoilingWaterExhaustedError(
                    f"packet {packet.packet_id} cannot supply the certified co-boiling "
                    f"water share: {removed_water_kg!r} kg demanded against "
                    f"{external_water_after_condensation_kg!r} kg of external free water; "
                    "refuse, do not clamp"
                )
            candidate = _packet_at_saturation(
                before,
                total_hexane_mass_kg=(
                    total_hexane_after_condensation - removed_hexane_kg
                ),
                saturation_temperature_k=saturation_temperature_k,
                codec=codec,
                external_free_water_mass_kg=water_after,
            )
            value = math.fsum(
                (
                    candidate.inventory.common_datum_energy_j,
                    removed_hexane_kg / properties.hexane_molar_mass_kg_mol * vapour_enthalpy,
                    -available_energy,
                )
            )
            return value, candidate

        residual_zero, _ = residual(0.0)
        if residual_zero > relative_limit * max(1.0, abs(available_energy)):
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} lacks the heat required to reach saturation"
            )

        # D12 level 0 (design site 4): armed, the upper bound is the packet's
        # whole hexane; the representability search below stops it at the
        # tag-2 codec's own floor, because every removal past the floor is
        # refused by the codec (``_LEVEL0_BRACKET_ERRORS``).
        attached_hexane_bound = (
            total_hexane_after_condensation
            if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec)
            else saturated_zero.state.attached_hexane_mass_kg
        )
        co_boiling_water_bound = (
            external_water_after_condensation_kg / water_per_hexane_kg
            if water_per_hexane_kg > 0.0
            else math.inf
        )
        invalid_high = min(attached_hexane_bound, co_boiling_water_bound)
        valid_high = 0.0
        high_candidate = saturated_zero
        level0_codec_floor_reached = False
        # EXPLORATION MODE F5: the refusal that bounded the bracket, kept for the breach row
        # (the ``except ... as`` name is unbound after its block).  Set only while armed.
        tag3_bracket_refusal: Exception | None = None
        for _ in range(90):
            probe = 0.5 * (valid_high + invalid_high)
            if probe in (valid_high, invalid_high):
                break
            try:
                probe_residual, candidate = residual(probe)
            except (
                # Every way a removal can be UNREPRESENTABLE shrinks the
                # bracket; this is the bracket search, not the verdict.  The
                # verdict is ``residual_high < 0`` below, which still refuses
                # typed.  F4 added the co-boiling water bound and the codec's
                # own subcritical refusal (reachable once the film is thin,
                # which the co-boiling flash makes routine).
                high_loading.HighLoadingPacketStateError,
                PDCoBoilingWaterExhaustedError,
                *(_LEVEL0_BRACKET_ERRORS if PD_FLASH_LEVEL0_ARMED else ()),
            ) as bracket_refusal:
                invalid_high = probe
                level0_codec_floor_reached = not isinstance(
                    bracket_refusal, PDCoBoilingWaterExhaustedError
                )
                if PD_TAG3_FREEZE_ARMED:
                    tag3_bracket_refusal = bracket_refusal
            else:
                valid_high = probe
                high_candidate = candidate
                if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec) and probe_residual >= 0.0:
                    # D12 level 0: this removal already absorbs the energy, so
                    # the root lies in [0, probe]; the floor is not needed.
                    break
        residual_high, _ = residual(valid_high)
        if (
            PD_TAG3_FREEZE_ARMED
            and PD_TAG3_FREEZE_HOOK is not None
            and exploration_call_kind in PD_TAG3_FREEZE_CALL_KINDS
            and residual_high < 0.0
            and PD_FLASH_LEVEL0_ARMED
            and _level0_is_dual(codec)
            and (
                level0_codec_floor_reached
                or not co_boiling_water_bound < attached_hexane_bound
            )
        ):
            # EXPLORATION MODE F5 (record GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13):
            # exactly the two PDFlashLevel0FloorRefusal raises below, frozen instead.  The
            # hook writes the PHYSICS breach row first, then flashes only down to the floor
            # (valid_high) and keeps the un-absorbed energy in the packet as superheat (the
            # WS-D construction); a row it cannot write, or a frozen state it cannot build,
            # raises the original refusal.  The co-boiling-water sibling is not rammed.
            (
                chosen,
                generated_per_particle,
                generated_water_per_particle,
                vapour_energy_per_particle,
            ) = PD_TAG3_FREEZE_HOOK.freeze_tag3(
                packet=packet,
                before=before,
                codec=codec,
                properties=properties,
                call_kind=exploration_call_kind,
                floor_reached=level0_codec_floor_reached,
                valid_high=valid_high,
                invalid_high=invalid_high,
                residual_high=residual_high,
                available_energy=available_energy,
                energy_input_per_particle=energy_input_per_particle,
                total_hexane_after_condensation=total_hexane_after_condensation,
                external_water_after_condensation_kg=external_water_after_condensation_kg,
                water_per_hexane_kg=water_per_hexane_kg,
                vapour_enthalpy=vapour_enthalpy,
                saturation_temperature_k=saturation_temperature_k,
                tray_pressure_pa=tray_pressure_pa,
                representatives=representatives,
                condensate_per_particle=condensate_per_particle,
                condensate_hexane_per_particle=condensate_hexane_per_particle,
                bracket_refusal=tag3_bracket_refusal,
            )
        elif residual_high < 0.0:
            if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec) and level0_codec_floor_reached:
                # D12 level 0 CROSSING: the codec's own floor bounded the
                # search, and flashing down to it does not absorb the energy.
                raise PDFlashLevel0FloorRefusal(
                    f"packet {packet.packet_id} would flash past the tag-2 floor in this "
                    "PD step: removing its inner hexane down to the floor at the "
                    "saturation temperature does not absorb the available energy; "
                    "tag-3 on a PD tray is outside level 0 - refuse, do not clamp"
                )
            if co_boiling_water_bound < attached_hexane_bound:
                raise PDCoBoilingWaterExhaustedError(
                    f"packet {packet.packet_id} would exhaust its external free water "
                    "before reaching the certified co-boiling temperature; refuse, "
                    "do not clamp"
                )
            if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec):
                # An armed search that ended at the tag-2 floor refused above
                # and one that ended at the co-boiling water bound refused with
                # the shipped water text; any other armed end refuses typed
                # here too, never through the WS-D handoff.
                raise PDFlashLevel0FloorRefusal(
                    f"packet {packet.packet_id} would flash past the tag-2 floor in this "
                    "PD step: removing its inner hexane down to the floor at the "
                    "saturation temperature does not absorb the available energy; "
                    "tag-3 on a PD tray is outside level 0 - refuse, do not clamp"
                )
            if type(codec) is not dual_variant.EngineeringDualVariantPacketAdapter:
                raise high_loading.HighLoadingVariantHandoffRequired(
                    f"packet {packet.packet_id} would exhaust attached hexane in this PD step"
                )
            # WS-D (owner-ruled AUTHORIZED 2026-09-04, design packet
            # ``GT_PS2_PD_PHASE_CHANGE_PROGRAM_DESIGN_2026-09-03.md`` S:2.4,
            # decision D2 "unconditional"): even removing ALL currently
            # attached hexane, pinned at ``saturation_temperature_k`` (this
            # local bisection's own basis), does not absorb the available
            # energy -- the packet needs to cross the tag-1 -> tag-2
            # boundary.  Rather than raising the tag-1 exhaustion refusal
            # itself, hand the packet to the codec ALREADY in hand: under a
            # dual-variant codec this is exactly
            # ``EngineeringDualVariantPacketAdapter.advance``, which tries
            # the tag-1 sub-adapter first and, on ITS OWN typed exhaustion
            # refusal, reconstitutes the SAME conserved totals under tag-2
            # (engineering_dual_variant_packet_adapter.py) -- the identical
            # handoff/reconstitute contract TAG-2 stages 1-7 already proved
            # working elsewhere in the engineering lane.  The exhausting
            # removal amount is ALL currently-attached hexane
            # (``attached_hexane_bound``); unlike the pinned-saturation
            # bisection above, the packet's remaining energy is not
            # balanced against a saturation-pinned candidate -- it is
            # ``available_energy`` LESS the energy the departing vapour
            # carries away (the same accounting the bisection's own
            # ``residual`` closure uses), and the codec's own caloric
            # inversion is free to find whatever temperature -- at or
            # above saturation -- that remaining energy implies.  This
            # module no longer preempts that decision with its own raise;
            # the codec's own union of exact domains decides what happens
            # next, including refusing typed if neither side admits the
            # resulting state (a REFUSAL DEEPER IN THE STACK, not this
            # module's own raise, is the intended outcome when it occurs).
            chosen_removed = attached_hexane_bound
            removed_water_kg = chosen_removed * water_per_hexane_kg
            water_after_handoff = math.fsum(
                (external_water_after_condensation_kg, -removed_water_kg)
            )
            if water_after_handoff < 0.0:
                raise PDCoBoilingWaterExhaustedError(
                    f"packet {packet.packet_id} cannot supply the certified co-boiling "
                    f"water share for the handoff: {removed_water_kg!r} kg demanded "
                    f"against {external_water_after_condensation_kg!r} kg of external "
                    "free water; refuse, do not clamp"
                )
            generated_per_particle = chosen_removed
            generated_water_per_particle = removed_water_kg
            vapour_energy_per_particle = (
                chosen_removed / properties.hexane_molar_mass_kg_mol * vapour_enthalpy
            )
            chosen = codec.advance(
                before.payload_bytes,
                total_hexane_mass_after_kg=math.fsum(
                    (total_hexane_after_condensation, -chosen_removed)
                ),
                external_free_water_mass_after_kg=water_after_handoff,
                envelope_energy_after_j=math.fsum(
                    (available_energy, -vapour_energy_per_particle)
                ),
            )
        else:
            low = 0.0
            high = valid_high
            low_pair = (residual_zero, saturated_zero)
            high_pair = (residual_high, high_candidate)
            if PD_FLASH_LEVEL0_ARMED and _level0_is_dual(codec):
                low, high, low_pair, high_pair = _level0_bracketed_root(
                    residual, low, high, low_pair, high_pair
                )
            for _ in range(110):
                midpoint = 0.5 * (low + high)
                if midpoint in (low, high):
                    break
                mid_pair = residual(midpoint)
                if mid_pair[0] <= 0.0:
                    low = midpoint
                    low_pair = mid_pair
                else:
                    high = midpoint
                    high_pair = mid_pair
            chosen_removed, (_, chosen) = min(
                ((low, low_pair), (high, high_pair)),
                key=lambda item: abs(item[1][0]),
            )
            generated_per_particle = math.fsum(
                (
                    total_hexane_after_condensation,
                    -math.fsum(
                        (
                            chosen.inventory.attached_hexane_kg,
                            chosen.inventory.internal_hexane_kg,
                        )
                    ),
                )
            )
            if generated_per_particle <= 0.0 or chosen_removed <= 0.0:
                raise PDKDomeBypassStepError(
                    f"packet {packet.packet_id} boiling solve produced no positive evaporation"
                )
            generated_water_per_particle = generated_per_particle * water_per_hexane_kg
            vapour_energy_per_particle = (
                generated_per_particle / properties.hexane_molar_mass_kg_mol * vapour_enthalpy
            )
        after = chosen

    if declared_mechanical_work_j > 0.0:
        if (
            type(thermal_guard_temperature_after_k) is not float
            or type(thermal_guard_generated_hexane_kg) is not float
        ):
            raise PDKDomeBypassStepError(
                "work-bearing packet solve requires a separate t_n thermal-only guard"
            )
        _require_float(
            "thermal-only packet temperature guard",
            thermal_guard_temperature_after_k,
            nonnegative=True,
        )
        _require_float(
            "thermal-only packet generation guard",
            thermal_guard_generated_hexane_kg,
            nonnegative=True,
        )
        guarded_temperature_after_k = thermal_guard_temperature_after_k
        guarded_generated_hexane_kg = thermal_guard_generated_hexane_kg
    else:
        guarded_temperature_after_k = after.state.temperature_k
        guarded_generated_hexane_kg = generated_per_particle * representatives
    driving_temperatures = [before.state.temperature_k, adjacent_old_temperature_k]
    if declared_wall_heat_j != 0.0:
        driving_temperatures.append(old_wall_temperature_k)
    condensation_endpoint_k: float | None = None
    if condensate_per_particle > 0.0 or condensate_hexane_per_particle > 0.0:
        # WS-A - RETIRED AND RESTATED, exactly as F4 restated the same envelope
        # for its flash, and scoped to the armed capability so nothing moves at
        # the declared default.
        #
        # The envelope exists to stop an explicit heat TRANSFER from pushing a
        # packet past its own driving temperatures.  Condensation is not such
        # a transfer: it is a phase change depositing latent heat at the
        # surface, and it STOPS on its own once the surface reaches the
        # arriving vapour's DEW POINT at its own water partial pressure.  That
        # dew point - not the vapour's bulk temperature, and not the wall - is
        # therefore the endpoint a condensing packet may relax onto, and it is
        # admitted explicitly.  Everything else still fails closed: a packet
        # driven past the dew point by condensation is still refused, which is
        # the state that cannot happen physically.
        #
        # SCOPED 2026-09-04 (verification pass).  The endpoint is kept OUT of
        # ``driving_temperatures`` and raises the ceiling of the PACKET
        # TEMPERATURE envelope alone.  The hexane-boiling clause below keeps
        # the pre-WS-A transfer ceiling, so an arriving condensate can never
        # license a boil at a saturation temperature the packet's own old
        # driving temperatures do not already admit, and neither wall heat nor
        # inter-K heat gains any room from the presence of condensate.
        if type(declared_condensation_equilibrium_temperature_k) is not float:
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} received condensate with no declared "
                "condensation equilibrium temperature"
            )
        condensation_endpoint_k = declared_condensation_equilibrium_temperature_k
    if begins_above_saturation:
        # F4 - RESTATED for the co-boiling flash.  The envelope exists to stop
        # an explicit heat TRANSFER from pushing a packet past its own driving
        # temperatures.  A flash to the certified two-liquid equilibrium is not
        # a transfer from any of them: it is the packet relaxing onto its OWN
        # boiling point, which the certified authority fixes.  That endpoint is
        # therefore admitted explicitly, and the envelope keeps binding on every
        # other direction (nothing may exceed the hottest driving temperature,
        # and nothing may fall below the coldest one OTHER than by reaching the
        # certified boiling point).
        driving_temperatures.append(saturation_temperature_k)
    admissible_low = min(driving_temperatures)
    #: The TRANSFER ceiling: the packet's own old driving temperatures, and
    #: nothing else - the pre-WS-A ceiling, byte for byte.
    transfer_admissible_high = max(driving_temperatures)
    #: The ceiling with the condensation endpoint admitted.  With no condensate
    #: the two are the same float, so the declared default is byte for byte the
    #: pre-WS-A envelope.
    admissible_high = (
        transfer_admissible_high
        if condensation_endpoint_k is None
        else max(transfer_admissible_high, condensation_endpoint_k)
    )
    if condensation_endpoint_k is not None:
        # SCOPED 2026-09-04 (verification pass).  The raised ceiling belongs to
        # the LATENT path, and to it alone.  Left global it would also hand the
        # step's ordinary transfers - wall heat, inter-K heat, mechanical work -
        # about 26 K of extra room the moment ANY condensate arrived, however
        # infinitesimal.  It does not: the same step is re-solved once with the
        # condensate and its energy REMOVED, and that condensation-free state
        # must still land inside the pre-WS-A envelope on its own.  Only the
        # difference between the two solves may use the dew point.  This costs
        # one extra codec advance per CONDENSING packet (top K layer, armed
        # capability only) and exactly nothing at the declared default, where
        # this branch is unreachable.
        try:
            condensation_free = codec.advance(
                before.payload_bytes,
                # CS-1a: the counterfactual removes BOTH condensates.
                total_hexane_mass_after_kg=total_hexane_before,
                external_free_water_mass_after_kg=inventory_before.external_water_kg,
                # F-CS1-1 correction (A): the DECLARED HEAT only.  Shaft
                # work is legitimately allowed to lift a packet above its own
                # driving temperatures, which is exactly why the ordinary
                # envelope check is run against a separate work-free solve.
                envelope_energy_after_j=math.fsum(
                    (inventory_before.common_datum_energy_j, declared_heat_per_particle)
                ),
            )
        except high_loading.HighLoadingCaloricBracketError as error:
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} condensation-free transfer state is not "
                "representable, so the raised envelope cannot be scoped to the "
                "latent path; reject, do not clamp"
            ) from error
        # Above saturation the condensation-free path would BOIL and pin at the
        # boiling temperature, so that is the transfer-only temperature there.
        transfer_only_temperature_k = (
            condensation_free.state.temperature_k
            if condensation_free.state.temperature_k
            <= saturation_temperature_k + temperature_tolerance
            else saturation_temperature_k
        )
        # F-CS1-1 correction (B), owner-ruled 2026-09-05: the declared floor.
        # See PD_COUNTERFACTUAL_TRANSFER_TOLERANCE_FLOOR_K.
        #
        # RESTATED 2026-09-06 by the owner ruling on point A.  Until that
        # ruling this floor applied HERE AND NOWHERE ELSE, at 2.0e-8 K, and
        # the sentence that stood in this place said exactly that.  It is no
        # longer true: the counterfactual floor was raised to the one
        # declared envelope floor (PD_TRANSFER_ENVELOPE_TOLERANCE_FLOOR_K =
        # 1.0e-6 K) so both guards carry one floor, and
        # ``temperature_tolerance`` itself now carries that same floor, so
        # this ``max`` is an identity at every saturation temperature this
        # step can reach.  It is kept rather than deleted because the two
        # floors are separately declared rows and either may be re-ruled on
        # its own evidence: this guard measures a RE-SOLVE of one stored
        # (T, E) pair through the codec, whose own stated consistency
        # (2.0e-8 K, high_loading_packet_adapter.py:1018, measured as
        # F-GATE-3) is a different quantity from the truncation error of a
        # real explicit transfer that the after-temperature check measures.
        counterfactual_tolerance = max(
            temperature_tolerance, PD_COUNTERFACTUAL_TRANSFER_TOLERANCE_FLOOR_K
        )
        if not (
            admissible_low - counterfactual_tolerance
            <= transfer_only_temperature_k
            <= transfer_admissible_high + counterfactual_tolerance
        ):
            raise PDKDomeBypassStepError(
                f"packet {packet.packet_id} explicit thermal transfer crossed its old "
                "driving temperature envelope without the condensation term; "
                "reject, do not clamp "
                f"(transfer-only {transfer_only_temperature_k!r} K outside "
                f"[{admissible_low!r}, {transfer_admissible_high!r}] K, tolerance "
                f"{counterfactual_tolerance!r} K; before {before.state.temperature_k!r} K, "
                f"adjacent {adjacent_old_temperature_k!r} K, wall "
                f"{old_wall_temperature_k!r} K, declared wall heat "
                f"{declared_wall_heat_j!r} J, declared inter-K heat "
                f"{declared_inter_k_heat_j!r} J, condensate water "
                f"{declared_condensate_water_kg!r} kg, condensate hexane "
                f"{declared_condensate_hexane_kg!r} kg)"
            )
    # THE AFTER-TEMPERATURE DRIVING-ENVELOPE GUARD.  The law, unchanged: an
    # explicit inter-K transfer may not carry the packet past the envelope
    # spanned by its own old temperature and its adjacent layer's old
    # temperature (raised to the condensation endpoint on the latent path
    # only, scoped above).  It is a monotonicity condition on an explicit
    # scheme, and it still REFUSES; it never clamps.
    #
    # RULED 2026-09-06 (point A): its tolerance now carries the declared
    # 1.0e-6 K materiality floor through ``temperature_tolerance`` above.
    # This is the guard F-GATE-1 measured: on the twelve-interval reference
    # fixture it is the ONLY guard that refuses, with excursions of 3.7e-8
    # to 4.6e-8 K at the 0.25 s stride and 1.8e-8 to 2.3e-8 K after one
    # halving, on driving envelopes 1.2-2.1 mK wide.  Because the excursion
    # is first order in the stride, what the floor admits is the explicit
    # scheme's truncation error and not an instability, bounded at 1.0e-6 K
    # per explicit transfer: one part in a thousand of the measured driving
    # envelope and 1e-5 of the ratified 0.1 K per-tray certificate bound.
    # An excursion ABOVE the floor is still refused here and is still cured
    # by halving the stride, which is the declared row's rejecting outcome.
    #
    # EXPLORATION MODE F5 cascade (lead's ruling for stage 3 of 2026-10-14, record
    # GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13): with the tag-3 freeze armed, a crossing
    # ABOVE the ceiling by no more than the superheat the freeze logged for this packet in this
    # same site attempt is admitted after its own PHYSICS row; anything else raises verbatim.
    if not (
        admissible_low - temperature_tolerance
        <= guarded_temperature_after_k
        <= admissible_high + temperature_tolerance
    ) and not (
        PD_TAG3_FREEZE_ARMED
        and PD_TAG3_FREEZE_HOOK is not None
        and PD_TAG3_FREEZE_HOOK.admit_frozen_envelope_crossing(
            packet_id=packet.packet_id,
            call_kind=exploration_call_kind,
            temperature_after_k=guarded_temperature_after_k,
            admissible_low_k=admissible_low,
            admissible_high_k=admissible_high,
            tolerance_k=temperature_tolerance,
            before_k=before.state.temperature_k,
            adjacent_k=adjacent_old_temperature_k,
            wall_k=old_wall_temperature_k,
            condensation_endpoint_k=condensation_endpoint_k,
        )
    ):
        raise PDKDomeBypassStepError(
            f"packet {packet.packet_id} explicit thermal transfer crossed its old driving "
            "temperature envelope; reject, do not clamp "
            f"(after {guarded_temperature_after_k!r} K outside [{admissible_low!r}, "
            f"{admissible_high!r}] K, tolerance {temperature_tolerance!r} K; before "
            f"{before.state.temperature_k!r} K, adjacent {adjacent_old_temperature_k!r} K, "
            f"wall {old_wall_temperature_k!r} K, condensation endpoint "
            f"{condensation_endpoint_k!r} K)"
        )
    # The boiling-saturation interval check reads the SAME
    # ``temperature_tolerance``, so it carries the same declared 1.0e-6 K
    # floor (2026-09-06, point A).  That is acceptable and deliberate: it is
    # an interval test asking whether the saturation temperature at which
    # the packet is said to boil lies inside the old packet/old wall
    # interval, and the floor widens that interval by 1.0e-6 K at each end.
    # The interval it tests is the driving envelope - millikelvin-wide at
    # its narrowest on the reference fixture, tens of kelvin wide wherever a
    # wall is driving - so a 1e-6 K widening cannot admit a saturation
    # temperature that lies outside the interval for any physical reason; it
    # only stops this check from firing on the same explicit-transfer
    # truncation the guard above now admits.  No boiling law, rate or root
    # moves.
    if guarded_generated_hexane_kg > 0.0 and not (
        admissible_low - temperature_tolerance
        <= saturation_temperature_k
        <= admissible_high + temperature_tolerance
    ):
        raise PDKDomeBypassStepError(
            f"packet {packet.packet_id} cannot boil at a saturation temperature outside "
            "the old packet/old wall interval"
        )

    total_hexane_after = math.fsum(
        (after.inventory.attached_hexane_kg, after.inventory.internal_hexane_kg)
    )
    # The native high-loading initializer may repartition attached/internal
    # hexane by one rounding unit while preserving the requested total.  The
    # branch solve, not that repartition noise, owns physical vapour generation.
    realized_generated = generated_per_particle
    # CS-1a extends the pre-CS-1 construction by one term: the row closes as
    # (after + generated - before - condensed).  The condensate is exactly 0.0
    # at the declared default, so this IS the pre-CS-1 expression there.
    mass_residual_per_particle = math.fsum(
        (
            total_hexane_after,
            realized_generated,
            -total_hexane_before,
            -condensate_hexane_per_particle,
        )
    )
    packet_energy_change_per_particle = math.fsum(
        (after.inventory.common_datum_energy_j, -inventory_before.common_datum_energy_j)
    )
    weighted_generated = realized_generated * representatives
    weighted_generated_water = generated_water_per_particle * representatives
    weighted_vapour_energy = vapour_energy_per_particle * representatives
    weighted_energy_change = packet_energy_change_per_particle * representatives
    weighted_water_before = inventory_before.external_water_kg * representatives
    weighted_water_after = after.inventory.external_water_kg * representatives
    if (inventory_before.external_water_kg > 0.0 or condensate_per_particle > 0.0) and (
        after.state.temperature_k >= water_saturation_temperature - water_tolerance
    ):
        # F4: external free water is no longer declared NONVOLATILE - it
        # co-boils at the certified heteroazeotrope - but it is still declared
        # SUBCOOLED against its own pure saturation point, because above that
        # the two-liquid heteroazeotrope does not exist at this pressure.
        raise PDKDomeBypassStepError(
            f"packet {packet.packet_id} external water left the declared subcooled "
            "co-boiling branch"
        )
    # F4: the external-water row now closes as (after - before + co-boiled),
    # so a co-boiling step is conservative in water rather than merely
    # water-inert.  On the single-liquid branch ``weighted_generated_water``
    # is exactly 0.0 and this is the pre-F4 expression.
    # WS-A extends F4's own construction by one term: the row now closes as
    # (after - before + co-boiled - condensed).  The condensate is exactly 0.0
    # at the declared default, so this IS the F4 expression there.
    weighted_water_residual = math.fsum(
        (
            weighted_water_after,
            weighted_generated_water,
            -weighted_water_before,
            -declared_condensate_water_kg,
        )
    )
    weighted_mass_residual = mass_residual_per_particle * representatives
    weighted_energy_residual = math.fsum(
        (
            weighted_energy_change,
            weighted_vapour_energy,
            -total_declared_packet_energy_input_j,
            -declared_condensation_energy_j,
        )
    )
    mass_roundoff_bound = math.fsum(
        (
            representatives
            * _ulp_sum(
                inventory_before.attached_hexane_kg,
                inventory_before.internal_hexane_kg,
                total_hexane_before,
                after.inventory.attached_hexane_kg,
                after.inventory.internal_hexane_kg,
                total_hexane_after,
                realized_generated,
                mass_residual_per_particle,
                # Appended ONLY when the capability is armed, so the default
                # bound is the pre-CS-1 float exactly.
                *(
                    (
                        total_hexane_after_condensation,
                        condensate_hexane_per_particle,
                    )
                    if declared_condensate_hexane_kg
                    else ()
                ),
            ),
            _ulp_sum(
                weighted_generated,
                weighted_mass_residual,
                *(
                    (declared_condensate_hexane_kg,)
                    if declared_condensate_hexane_kg
                    else ()
                ),
            ),
        )
    )
    water_roundoff_bound = _ulp_sum(
        weighted_water_before,
        weighted_water_after,
        weighted_generated_water,
        weighted_water_residual,
        # Appended ONLY when the capability is armed, so the default bound is
        # the pre-WS-A float exactly.
        *((declared_condensate_water_kg,) if declared_condensate_water_kg else ()),
    )
    energy_roundoff_bound = math.fsum(
        (
            representatives
            * _ulp_sum(
                after.inventory.common_datum_energy_j,
                inventory_before.common_datum_energy_j,
                packet_energy_change_per_particle,
                vapour_energy_per_particle,
                heat_per_particle,
                *(
                    (condensation_energy_per_particle,)
                    if declared_condensation_energy_j
                    else ()
                ),
            ),
            _ulp_sum(
                weighted_energy_change,
                weighted_vapour_energy,
                total_declared_packet_energy_input_j,
                weighted_energy_residual,
                *(
                    (declared_condensation_energy_j,)
                    if declared_condensation_energy_j
                    else ()
                ),
            ),
        )
    )
    water_relative_residual = _roundoff_aware_relative_residual(
        weighted_water_residual,
        water_roundoff_bound,
        weighted_water_before,
        weighted_water_after,
        weighted_generated_water,
    )
    mass_relative_residual = _roundoff_aware_relative_residual(
        weighted_mass_residual,
        mass_roundoff_bound,
        total_hexane_before * representatives,
        total_hexane_after * representatives,
        weighted_generated,
    )
    former_energy_relative_residual = None
    if packet_ledger_stock_scaled_enabled():
        # Ruling 1 of 2026-10-11 (option C): the stocks differenced join the
        # denominator, as the water and mass rows above already do.  The former
        # form is kept for the refusal text only.
        former_energy_relative_residual = _roundoff_aware_relative_residual(
            weighted_energy_residual,
            energy_roundoff_bound,
            weighted_energy_change,
            weighted_vapour_energy,
            total_declared_packet_energy_input_j,
        )
        energy_relative_residual = _roundoff_aware_relative_residual(
            weighted_energy_residual,
            energy_roundoff_bound,
            weighted_energy_change,
            weighted_vapour_energy,
            total_declared_packet_energy_input_j,
            after.inventory.common_datum_energy_j * representatives,
            inventory_before.common_datum_energy_j * representatives,
        )
    else:
        energy_relative_residual = _roundoff_aware_relative_residual(
            weighted_energy_residual,
            energy_roundoff_bound,
            weighted_energy_change,
            weighted_vapour_energy,
            total_declared_packet_energy_input_j,
        )
    # The same three values in the same order as before this build; max() of
    # the same binary64 operands is the same binary64 result.
    maximum = max(water_relative_residual, mass_relative_residual, energy_relative_residual)
    ledger = PDPacketHeatLedger(
        packet_id=packet.packet_id,
        vertical_layer=packet.owner_key.vertical_layer_id.value,
        temperature_before_k=before.state.temperature_k,
        temperature_after_k=after.state.temperature_k,
        saturation_temperature_k=saturation_temperature_k,
        declared_wall_heat_j=declared_wall_heat_j,
        declared_inter_k_heat_j=declared_inter_k_heat_j,
        declared_mechanical_work_j=declared_mechanical_work_j,
        total_declared_heat_j=total_declared_heat_j,
        total_declared_packet_energy_input_j=total_declared_packet_energy_input_j,
        thermal_only_temperature_after_k=guarded_temperature_after_k,
        thermal_only_generated_hexane_kg=guarded_generated_hexane_kg,
        external_water_before_kg=weighted_water_before,
        external_water_after_kg=weighted_water_after,
        external_water_residual_kg=weighted_water_residual,
        generated_hexane_kg=weighted_generated,
        generated_water_kg=weighted_generated_water,
        generated_hexane_common_energy_j=weighted_vapour_energy,
        packet_common_energy_change_j=weighted_energy_change,
        mass_residual_kg=weighted_mass_residual,
        common_energy_residual_j=weighted_energy_residual,
        mass_binary64_roundoff_bound_kg=mass_roundoff_bound,
        common_energy_binary64_roundoff_bound_j=energy_roundoff_bound,
        maximum_relative_residual=maximum,
        relative_limit=relative_limit,
        condensed_water_kg=declared_condensate_water_kg,
        condensation_common_energy_j=declared_condensation_energy_j,
        condensed_hexane_kg=declared_condensate_hexane_kg,
    )
    # EXPLORATION MODE F4 (record GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13): at the
    # in-force and thermal-guard calls the hook writes the packet-ledger row first and then
    # admits this ledger; a row it cannot write keeps the refusal below, verbatim.
    if not ledger.passed and not (
        PD_EXPLORATION_HOOK is not None
        and PD_EXPLORATION_HOOK.ram_packet_ledger(
            ledger,
            packet_id=packet.packet_id,
            call_kind=exploration_call_kind,
            rows=(
                ("water", weighted_water_residual, water_roundoff_bound, water_relative_residual),
                ("mass", weighted_mass_residual, mass_roundoff_bound, mass_relative_residual),
                (
                    "common-energy",
                    weighted_energy_residual,
                    energy_roundoff_bound,
                    energy_relative_residual,
                ),
            ),
            maximum=maximum,
            limit=relative_limit,
            former_energy_relative=former_energy_relative_residual,
        )
    ):
        raise PDKDomeBypassStepError(
            packet_ledger_refusal_text(
                packet.packet_id,
                rows=(
                    (
                        ("water", weighted_water_residual, water_roundoff_bound, water_relative_residual),
                        ("mass", weighted_mass_residual, mass_roundoff_bound, mass_relative_residual),
                        (
                            "common-energy",
                            weighted_energy_residual,
                            energy_roundoff_bound,
                            energy_relative_residual,
                        ),
                    )
                    if packet_ledger_rows_enabled()
                    else None
                ),
                maximum=maximum,
                limit=relative_limit,
                former_energy_relative=former_energy_relative_residual,
            )
        )
    entry = replace(
        packet.entry,
        payload_bytes=after.payload_bytes,
        inventory=after.inventory,
        surface=after.surface,
    )
    if PD_TAG2_LEG_ARMED and PD_TAG2_LEG_HOOK is not None:
        PD_TAG2_LEG_HOOK.record(
            tag2_leg_choice,
            packet=packet,
            before=before,
            after=after,
            generated_per_particle=generated_per_particle,
            condensate_hexane_per_particle=condensate_hexane_per_particle,
            radius_m=codec.sphere_params.particle_radius_m,
        )
    return entry, ledger


@dataclass(frozen=True, slots=True, kw_only=True)
class _AcceptedPDPacketAdvance:
    entries: dict[str, tray_host.PacketPlanEntry]
    residuals: dict[str, float]
    maximum_scaled_residual: float
    callback_id: str = PACKET_CALLBACK_ID
    source_identity: str = PACKET_CALLBACK_SOURCE

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def propose(
        self, request: tray_host.ManufacturedPacketAdvanceRequest
    ) -> tray_host.ManufacturedPacketAdvanceProposal:
        return tray_host.ManufacturedPacketAdvanceProposal(
            packet_id=request.packet.packet_id,
            entry_after=self.entries[request.packet.packet_id],
            callback_id=self.callback_id,
            source_identity=self.source_identity,
            residual_contract_passed=True,
            maximum_scaled_residual=self.residuals[request.packet.packet_id],
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDBypassPressureLedger:
    bypass_inlet_pressure_pa: float
    prescribed_pressure_drop_pa: float
    plenum_pressure_pa: float
    pressure_residual_pa: float

    physically_qualifying: ClassVar[bool] = False
    pressure_drop_is_prescribed_not_predicted: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        return self.pressure_residual_pa == 0.0


@dataclass(frozen=True, slots=True, kw_only=True)
class PDPlenumLedger:
    enthalpy_convention: str
    inlet: PDBinaryGasStream
    outlet: PDBinaryGasStream
    generated_hexane_kg: float
    #: F4: external free water co-boiled into the plenum with the n-hexane.
    generated_water_kg: float
    generated_hexane_mol: float
    #: Total generated vapour (hexane plus its certified co-boiling water).
    generated_vapour_mol: float
    generated_hexane_common_energy_j: float
    wall_to_bypass_gas_energy_j: float
    hexane_molar_residual: float
    water_molar_residual: float
    common_energy_residual_j: float
    maximum_relative_residual: float
    relative_limit: float
    #: WS-A: water that LEFT this stream onto the meal's top free surface, and
    #: the common-datum energy it took with it.  Exactly 0.0 at the declared
    #: default contact fraction.
    condensed_water_mol: float = 0.0
    condensed_water_kg: float = 0.0
    condensation_common_energy_j: float = 0.0
    #: CS-1a: n-hexane that LEFT this stream onto the same cold surface.
    #: Exactly 0.0 at the declared default contact fraction.
    condensed_hexane_mol: float = 0.0
    condensed_hexane_kg: float = 0.0

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return self.maximum_relative_residual <= self.relative_limit


def _plenum_outlet_and_ledger(
    request: PDKTrayStepRequest,
    *,
    saturation_temperature_k: float,
    generated_hexane_kg: float,
    generated_water_kg: float,
    generated_hexane_energy_j: float,
    co_boiling: PDCoBoilingState | None,
    wall_to_bypass_gas_energy_j: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
    condensation: PDContactCondensationLedger | None = None,
) -> tuple[PDBinaryGasStream, PDPlenumLedger]:
    duration = request.duration_s
    # WS-A: the condensate is a SINK on this stream exactly as the generated
    # vapour is a source.  The ledger keeps the PHYSICAL inlet, so the whole
    # tray's system ledger still reads the boundary as (inlet - outlet) and
    # sees the condensate without any further change.
    condensed_water_moles = 0.0 if condensation is None else condensation.condensed_water_mol
    condensed_water_kg = 0.0 if condensation is None else condensation.condensed_water_kg
    condensed_hexane_moles = 0.0 if condensation is None else condensation.hexane_condensed_mol
    condensed_hexane_kg = 0.0 if condensation is None else condensation.hexane_condensed_kg
    condensation_energy_j = (
        0.0 if condensation is None else condensation.condensation_common_energy_j
    )
    properties = request.model.properties
    inlet = request.bypass_inlet
    generated_hexane_moles = generated_hexane_kg / properties.hexane_molar_mass_kg_mol
    generated_water_moles = generated_water_kg / properties.water_molar_mass_kg_mol
    # F4 (2026-09-02): the generated vapour is the CERTIFIED co-boiling
    # mixture, so both of its components join the plenum.  ``generated_moles``
    # is the total vapour, not the hexane alone.
    generated_moles = math.fsum((generated_hexane_moles, generated_water_moles))
    outlet_flow = math.fsum(
        (
            inlet.molar_flow_mol_s,
            generated_moles / duration,
            -condensed_water_moles / duration,
            # CS-1a: the hexane condensate is a sink on this stream exactly as
            # the water condensate is.  0.0 at the declared default.
            -condensed_hexane_moles / duration,
        )
    )
    outlet_hexane_flow = math.fsum(
        (
            inlet.molar_flow_mol_s * inlet.hexane_mole_fraction,
            generated_hexane_moles / duration,
            -condensed_hexane_moles / duration,
        )
    )
    outlet_hexane_fraction = outlet_hexane_flow / outlet_flow
    outlet_water_fraction = 1.0 - outlet_hexane_fraction
    inlet_h = _bulk_gas_enthalpy_j_mol(
        properties,
        temperature_k=inlet.temperature_k,
        pressure_pa=inlet.pressure_pa,
        hexane_mole_fraction=inlet.hexane_mole_fraction,
        water_mole_fraction=inlet.water_mole_fraction,
        datum_adapter=codec.component_datum_adapter,
    )
    outlet_energy_rate = math.fsum(
        (
            inlet.molar_flow_mol_s * inlet_h,
            generated_hexane_energy_j / duration,
            wall_to_bypass_gas_energy_j / duration,
            -condensation_energy_j / duration,
        )
    )
    if type(properties) is cc.DeclaredGasProperties:
        hexane_shift, water_shift = codec.component_datum_adapter.common_datum_species_shifts_j_mol
        species_offset = math.fsum(
            (
                outlet_hexane_fraction * (properties.hexane_latent_heat_j_mol + hexane_shift),
                outlet_water_fraction * (properties.water_latent_heat_j_mol + water_shift),
            )
        )
        outlet_temperature = math.fsum(
            (
                properties.energy_datum_temperature_k,
                (outlet_energy_rate / outlet_flow - species_offset)
                / properties.molar_heat_capacity_j_mol_k,
            )
        )
    else:
        # CELL-02e-pd native arm: the plenum energy balance is closed by
        # inverting the frozen mixture enthalpy (bracketed Newton on the
        # monotone map) at the plenum pressure - no linearization.
        try:
            outlet_temperature = native_caloric.temperature_from_mixture_enthalpy(
                outlet_energy_rate / outlet_flow,
                request.prior.tray.pressure_pa,
                outlet_hexane_fraction,
            )
        except ValueError as error:
            raise PDKDomeBypassStepError(
                f"native plenum inverse caloric refused: {error}"
            ) from error
    if not math.isfinite(outlet_temperature) or outlet_temperature <= 0.0:
        raise PDKDomeBypassStepError("plenum common-energy algebra produced an invalid temperature")
    outlet = PDBinaryGasStream(
        pressure_pa=request.prior.tray.pressure_pa,
        molar_flow_mol_s=outlet_flow,
        hexane_mole_fraction=outlet_hexane_fraction,
        water_mole_fraction=outlet_water_fraction,
        temperature_k=outlet_temperature,
    )
    outlet_h = _bulk_gas_enthalpy_j_mol(
        properties,
        temperature_k=outlet.temperature_k,
        pressure_pa=outlet.pressure_pa,
        hexane_mole_fraction=outlet.hexane_mole_fraction,
        water_mole_fraction=outlet.water_mole_fraction,
        datum_adapter=codec.component_datum_adapter,
    )
    hexane_residual = math.fsum(
        (
            outlet.molar_flow_mol_s * outlet.hexane_mole_fraction * duration,
            -inlet.molar_flow_mol_s * inlet.hexane_mole_fraction * duration,
            -generated_hexane_moles,
            condensed_hexane_moles,
        )
    )
    water_residual = math.fsum(
        (
            outlet.molar_flow_mol_s * outlet.water_mole_fraction * duration,
            -inlet.molar_flow_mol_s * inlet.water_mole_fraction * duration,
            -generated_water_moles,
            condensed_water_moles,
        )
    )
    energy_residual = math.fsum(
        (
            outlet.molar_flow_mol_s * outlet_h * duration,
            -inlet.molar_flow_mol_s * inlet_h * duration,
            -generated_hexane_energy_j,
            -wall_to_bypass_gas_energy_j,
            condensation_energy_j,
        )
    )
    maximum = max(
        _relative_residual(hexane_residual, generated_hexane_moles, condensed_hexane_moles),
        _relative_residual(
            water_residual,
            inlet.molar_flow_mol_s * inlet.water_mole_fraction * duration,
            generated_water_moles,
            condensed_water_moles,
        ),
        _relative_residual(
            energy_residual,
            outlet.molar_flow_mol_s * outlet_h * duration,
            inlet.molar_flow_mol_s * inlet_h * duration,
            generated_hexane_energy_j,
            wall_to_bypass_gas_energy_j,
            condensation_energy_j,
        ),
    )
    ledger = PDPlenumLedger(
        enthalpy_convention=GAS_ENTHALPY_CONVENTION,
        inlet=inlet,
        outlet=outlet,
        generated_hexane_kg=generated_hexane_kg,
        generated_water_kg=generated_water_kg,
        generated_hexane_mol=generated_hexane_moles,
        generated_vapour_mol=generated_moles,
        generated_hexane_common_energy_j=generated_hexane_energy_j,
        wall_to_bypass_gas_energy_j=wall_to_bypass_gas_energy_j,
        hexane_molar_residual=hexane_residual,
        water_molar_residual=water_residual,
        common_energy_residual_j=energy_residual,
        maximum_relative_residual=maximum,
        relative_limit=request.relative_limit,
        condensed_water_mol=condensed_water_moles,
        condensed_water_kg=condensed_water_kg,
        condensation_common_energy_j=condensation_energy_j,
        condensed_hexane_mol=condensed_hexane_moles,
        condensed_hexane_kg=condensed_hexane_kg,
    )
    if not ledger.passed:
        raise PDKDomeBypassStepError("dome-bypass plenum component/common-energy ledger failed")
    # F4: the generated vapour closes against the enthalpy of the branch that
    # produced it - the certified co-boiling mixture when external free water
    # is present, the pure n-hexane vapour otherwise (pre-F4 arithmetic).
    generation_h = _vapour_enthalpy_per_hexane_mole_j(
        properties,
        co_boiling=co_boiling,
        saturation_temperature_k=saturation_temperature_k,
        pressure_pa=request.prior.tray.pressure_pa,
        datum_adapter=codec.component_datum_adapter,
    )
    generation_residual = math.fsum(
        (generated_hexane_energy_j, -generated_hexane_moles * generation_h)
    )
    if _relative_residual(generation_residual, generated_hexane_energy_j) > request.relative_limit:
        raise PDKDomeBypassStepError("generated vapour violates its declared enthalpy closure")
    return outlet, ledger


@dataclass(frozen=True, slots=True, kw_only=True)
class PDSharedWallLedger:
    old_wall_temperature_k: float
    new_wall_temperature_k: float
    accepted_indirect_duty_j: float
    wall_to_bottom_solid_energy_j: float
    wall_to_bypass_gas_energy_j: float
    ambient_loss_energy_j: float
    wall_storage_energy_j: float
    common_energy_residual_j: float
    no_explicit_equilibrium_crossing: bool
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    one_physical_wall: ClassVar[bool] = True
    shared_old_wall_evaluation: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        return (
            self.no_explicit_equilibrium_crossing
            and self.maximum_relative_residual <= self.relative_limit
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDSystemLedger:
    packet_hexane_change_kg: float
    cell_hexane_change_kg: float
    external_hexane_to_system_kg: float
    hexane_mass_residual_kg: float
    packet_water_change_kg: float
    cell_water_change_kg: float
    external_water_to_system_kg: float
    water_mass_residual_kg: float
    packet_common_energy_change_j: float
    cell_common_energy_change_j: float
    thermal_external_common_energy_to_system_j: float
    q_other_s_external_energy_j: float
    external_common_energy_to_system_j: float
    common_energy_residual_j: float
    wall_cell_storage_residual_j: float
    upper_cell_change_residuals: tuple[tuple[float, float, float], ...]
    upper_cell_change_binary64_roundoff_bounds: tuple[tuple[float, float, float], ...]
    hexane_binary64_roundoff_bound_kg: float
    water_binary64_roundoff_bound_kg: float
    common_energy_binary64_roundoff_bound_j: float
    wall_cell_binary64_roundoff_bound_j: float
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.q_other_s_external_energy_j >= 0.0
            and self.external_common_energy_to_system_j
            == math.fsum(
                (
                    self.thermal_external_common_energy_to_system_j,
                    self.q_other_s_external_energy_j,
                )
            )
            and self.maximum_relative_residual <= self.relative_limit
        )


def _weighted_packet_totals(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
) -> tuple[float, float, float]:
    water = math.fsum(
        packet.weighted_totals.totals[4] + packet.weighted_totals.totals[5] for packet in packets
    )
    hexane_mass = math.fsum(
        packet.weighted_totals.totals[2] + packet.weighted_totals.totals[3] for packet in packets
    )
    energy = math.fsum(packet.weighted_totals.totals[6] for packet in packets)
    return water, hexane_mass, energy


def _cell_totals(cells: tuple[tray_host.AcceptedKCell, ...]) -> tuple[float, float, float]:
    return (
        math.fsum(cell.conserved_state.totals.water_kg for cell in cells),
        math.fsum(cell.conserved_state.totals.hexane_kg for cell in cells),
        math.fsum(cell.conserved_state.totals.common_datum_energy_j for cell in cells),
    )


def _moisture_pickup_ledger(
    request: PDKTrayStepRequest,
    host_step: tray_host.ValidatedKCellTrayStep,
    *,
    condensation: PDContactCondensationLedger | None,
) -> PDMoisturePickupLedger:
    """WS-A / Q-PD-3: what this step did to the meal's own moisture.

    REPORT-ONLY by owner ruling.  The wet-basis moisture is
    ``(external + retained water) / (dry + oil + hexane + water)`` on the
    tray's whole packet population, before and after, and it is CHECKED - not
    enforced - against literature criterion A1.
    """

    before_water, before_hexane, _ = _weighted_packet_totals(request.prior.tray.packets)
    after_water, after_hexane, _ = _weighted_packet_totals(host_step.candidate.packets)
    before_dry = math.fsum(item.weighted_totals.totals[0] for item in request.prior.tray.packets)
    after_dry = math.fsum(item.weighted_totals.totals[0] for item in host_step.candidate.packets)
    before_oil = math.fsum(item.weighted_totals.totals[1] for item in request.prior.tray.packets)
    after_oil = math.fsum(item.weighted_totals.totals[1] for item in host_step.candidate.packets)
    before_wet = math.fsum((before_dry, before_oil, before_hexane, before_water))
    after_wet = math.fsum((after_dry, after_oil, after_hexane, after_water))
    condensed_kg = 0.0 if condensation is None else condensation.condensed_water_kg
    condensed_hexane_kg = 0.0 if condensation is None else condensation.hexane_condensed_kg
    moisture_before = before_water / before_wet if before_wet > 0.0 else 0.0
    moisture_after = after_water / after_wet if after_wet > 0.0 else 0.0
    low, high = PD_LITERATURE_MOISTURE_BAND_WET_BASIS
    return PDMoisturePickupLedger(
        physical_tray_id=request.prior.tray.physical_tray_id,
        contact_fraction=request.model.condensation_contact_fraction,
        condensed_water_kg=condensed_kg,
        dry_matter_kg=before_dry,
        condensate_fraction_of_dry_matter=(
            condensed_kg / before_dry if before_dry > 0.0 else 0.0
        ),
        total_water_before_kg=before_water,
        total_water_after_kg=after_water,
        meal_moisture_before_wet_basis=moisture_before,
        meal_moisture_after_wet_basis=moisture_after,
        literature_band_wet_basis=PD_LITERATURE_MOISTURE_BAND_WET_BASIS,
        literature_band_source=PD_LITERATURE_MOISTURE_BAND_SOURCE,
        moisture_before_inside_literature_band=low <= moisture_before <= high,
        moisture_after_inside_literature_band=low <= moisture_after <= high,
        condensed_hexane_kg=condensed_hexane_kg,
        condensate_hexane_fraction_of_dry_matter=(
            condensed_hexane_kg / before_dry if before_dry > 0.0 else 0.0
        ),
        total_hexane_before_kg=before_hexane,
        total_hexane_after_kg=after_hexane,
        meal_hexane_loading_before_dry_basis=(
            before_hexane / before_dry if before_dry > 0.0 else 0.0
        ),
        meal_hexane_loading_after_dry_basis=(
            after_hexane / after_dry if after_dry > 0.0 else 0.0
        ),
    )


def _system_ledger(
    request: PDKTrayStepRequest,
    host_step: tray_host.ValidatedKCellTrayStep,
    *,
    wall: PDSharedWallLedger,
    plenum: PDPlenumLedger,
    deposition_bundle: meal_work.QSCMealMechanicalWorkDepositionBundle | None,
) -> PDSystemLedger:
    packet_before = _weighted_packet_totals(request.prior.tray.packets)
    packet_after = _weighted_packet_totals(host_step.candidate.packets)
    cell_before = _cell_totals(request.prior.tray.cells)
    cell_after = _cell_totals(host_step.candidate.cells)
    packet_delta = tuple(
        math.fsum((after, -before)) for before, after in zip(packet_before, packet_after)
    )
    cell_delta = tuple(
        math.fsum((after, -before)) for before, after in zip(cell_before, cell_after)
    )
    duration = request.duration_s
    properties = request.model.properties
    inlet = plenum.inlet
    outlet = plenum.outlet
    external_hexane = (
        duration
        * properties.hexane_molar_mass_kg_mol
        * math.fsum(
            (
                inlet.molar_flow_mol_s * inlet.hexane_mole_fraction,
                -outlet.molar_flow_mol_s * outlet.hexane_mole_fraction,
            )
        )
    )
    external_water = (
        duration
        * properties.water_molar_mass_kg_mol
        * math.fsum(
            (
                inlet.molar_flow_mol_s * inlet.water_mole_fraction,
                -outlet.molar_flow_mol_s * outlet.water_mole_fraction,
            )
        )
    )
    # Plenum energy is already on the exact adapter-bound common datum.  Its
    # signed boundary contribution is H_in-H_out; by its own ledger this equals
    # -(generated vapour + wall-to-gas heat).
    thermal_external_energy = math.fsum(
        (
            wall.accepted_indirect_duty_j,
            -wall.ambient_loss_energy_j,
            -plenum.generated_hexane_common_energy_j,
            -plenum.wall_to_bypass_gas_energy_j,
            # WS-A: the condensate carried this much common-datum energy OUT of
            # the plenum and INTO the meal, so the boundary contribution
            # H_in - H_out gains it with a plus sign.  0.0 at the default.
            plenum.condensation_common_energy_j,
        )
    )
    q_other_s_external_energy = (
        0.0
        if deposition_bundle is None
        else deposition_bundle.ledger.emitted_packet_mechanical_work_j
    )
    external_energy = math.fsum((thermal_external_energy, q_other_s_external_energy))
    hexane_residual = math.fsum((packet_delta[1], cell_delta[1], -external_hexane))
    water_residual = math.fsum((packet_delta[0], cell_delta[0], -external_water))
    energy_residual = math.fsum((packet_delta[2], cell_delta[2], -external_energy))
    bottom_before = request.prior.tray.cells[0].conserved_state.totals
    bottom_after = host_step.candidate.cells[0].conserved_state.totals
    wall_cell_residual = math.fsum(
        (
            bottom_after.common_datum_energy_j,
            -bottom_before.common_datum_energy_j,
            -wall.wall_storage_energy_j,
        )
    )
    upper_residuals = tuple(
        (
            math.fsum(
                (after.conserved_state.totals.water_kg, -before.conserved_state.totals.water_kg)
            ),
            math.fsum(
                (after.conserved_state.totals.hexane_kg, -before.conserved_state.totals.hexane_kg)
            ),
            math.fsum(
                (
                    after.conserved_state.totals.common_datum_energy_j,
                    -before.conserved_state.totals.common_datum_energy_j,
                )
            ),
        )
        for before, after in zip(
            request.prior.tray.cells[1:], host_step.candidate.cells[1:], strict=True
        )
    )
    # For every non-wall cell the host evaluates C+ = Cn + Pn - P+ + X,
    # while this PD adapter declares X = P+ - Pn.  The exact accounting-cell
    # change is therefore zero.  Packet rows, their declared transfer, and the
    # host cell update are independently rounded, so remove only the ULPs of
    # those actual operands; do not enlarge the physical residual tolerance.
    candidate_packets_by_id = {packet.packet_id: packet for packet in host_step.candidate.packets}
    upper_roundoff_bounds: list[tuple[float, float, float]] = []
    for cell_index, (before_cell, after_cell, transfer) in enumerate(
        zip(
            request.prior.tray.cells[1:],
            host_step.candidate.cells[1:],
            host_step.request.external_cell_transfers[1:],
            strict=True,
        ),
        start=1,
    ):
        layer = before_cell.vertical_layer_id.value
        before_packets = tuple(
            packet
            for packet in request.prior.tray.packets
            if packet.owner_key.vertical_layer_id.value == layer
        )
        after_packets = tuple(
            candidate_packets_by_id[packet.packet_id] for packet in before_packets
        )
        before_water_rows = tuple(
            packet.weighted_totals.totals[4] + packet.weighted_totals.totals[5]
            for packet in before_packets
        )
        after_water_rows = tuple(
            packet.weighted_totals.totals[4] + packet.weighted_totals.totals[5]
            for packet in after_packets
        )
        before_hexane_rows = tuple(
            packet.weighted_totals.totals[2] + packet.weighted_totals.totals[3]
            for packet in before_packets
        )
        after_hexane_rows = tuple(
            packet.weighted_totals.totals[2] + packet.weighted_totals.totals[3]
            for packet in after_packets
        )
        before_energy_rows = tuple(packet.weighted_totals.totals[6] for packet in before_packets)
        after_energy_rows = tuple(packet.weighted_totals.totals[6] for packet in after_packets)
        before_packet_totals = (
            math.fsum(before_water_rows),
            math.fsum(before_hexane_rows),
            math.fsum(before_energy_rows),
        )
        after_packet_totals = (
            math.fsum(after_water_rows),
            math.fsum(after_hexane_rows),
            math.fsum(after_energy_rows),
        )
        before_cell_totals = before_cell.conserved_state.totals
        after_cell_totals = after_cell.conserved_state.totals
        residuals = upper_residuals[cell_index - 1]
        upper_roundoff_bounds.append(
            (
                _ulp_sum(
                    before_cell_totals.water_kg,
                    after_cell_totals.water_kg,
                    before_packet_totals[0],
                    after_packet_totals[0],
                    transfer.water_to_cell_kg,
                    residuals[0],
                    *before_water_rows,
                    *after_water_rows,
                ),
                _ulp_sum(
                    before_cell_totals.hexane_kg,
                    after_cell_totals.hexane_kg,
                    before_packet_totals[1],
                    after_packet_totals[1],
                    transfer.hexane_to_cell_kg,
                    residuals[1],
                    *before_hexane_rows,
                    *after_hexane_rows,
                ),
                _ulp_sum(
                    before_cell_totals.common_datum_energy_j,
                    after_cell_totals.common_datum_energy_j,
                    before_packet_totals[2],
                    after_packet_totals[2],
                    transfer.common_datum_energy_to_cell_j,
                    residuals[2],
                    *before_energy_rows,
                    *after_energy_rows,
                ),
            )
        )
    upper_roundoff_bounds_tuple = tuple(upper_roundoff_bounds)
    hexane_roundoff_bound = _ulp_sum(
        packet_before[1],
        packet_after[1],
        packet_delta[1],
        cell_before[1],
        cell_after[1],
        cell_delta[1],
        external_hexane,
        hexane_residual,
        # CS-1a: appended ONLY when a hexane condensate actually crossed, so
        # the declared default bound is the pre-CS-1 float exactly.
        *((plenum.condensed_hexane_kg,) if plenum.condensed_hexane_kg else ()),
    )
    water_roundoff_bound = _ulp_sum(
        packet_before[0],
        packet_after[0],
        packet_delta[0],
        cell_before[0],
        cell_after[0],
        cell_delta[0],
        external_water,
        water_residual,
    )
    energy_roundoff_bound = math.fsum(
        (
            _ulp_sum(
                packet_before[2],
                packet_after[2],
                packet_delta[2],
                cell_before[2],
                cell_after[2],
                cell_delta[2],
                external_energy,
                thermal_external_energy,
                q_other_s_external_energy,
                energy_residual,
            ),
            math.fsum(
                math.ulp(packet.weighted_totals.totals[6])
                for packet in (*request.prior.tray.packets, *host_step.candidate.packets)
            ),
            math.fsum(
                math.ulp(cell.conserved_state.totals.common_datum_energy_j)
                for cell in (*request.prior.tray.cells, *host_step.candidate.cells)
            ),
            (
                0.0
                if deposition_bundle is None
                else deposition_bundle.ledger.binary64_roundoff_bound_j
            ),
        )
    )
    wall_cell_roundoff_bound = _ulp_sum(
        bottom_before.common_datum_energy_j,
        bottom_after.common_datum_energy_j,
        wall.wall_storage_energy_j,
        wall_cell_residual,
    )
    # B3 (owner-ruled GO, 2026-08-31): the relative base includes the
    # CONSERVED-STATE magnitudes alongside the deltas and boundary flows.
    # The prior flow-only base collapsed to the 1.0 floor at small strides
    # (flows scale with duration; the solve's termination-scale closure
    # residual does not), turning this relative criterion into an absolute
    # one that no lawful solve can pass - measured live at the B1 stage-2
    # 28 us landing (rpm9_b1s2) and at setpoint 0.95 through a 154-attempt
    # rescue.  A residual of X on an inventory of S is X/S relative; the
    # limit is unchanged and the roundoff bounds still remove only explicit
    # binary64 operation noise.
    maximum = max(
        _roundoff_aware_relative_residual(
            hexane_residual,
            hexane_roundoff_bound,
            packet_before[1],
            cell_before[1],
            packet_delta[1],
            external_hexane,
        ),
        _roundoff_aware_relative_residual(
            water_residual,
            water_roundoff_bound,
            packet_before[0],
            cell_before[0],
            packet_delta[0],
            external_water,
        ),
        _roundoff_aware_relative_residual(
            energy_residual,
            energy_roundoff_bound,
            packet_before[2],
            cell_before[2],
            packet_delta[2],
            cell_delta[2],
            external_energy,
        ),
        _roundoff_aware_relative_residual(
            wall_cell_residual,
            wall_cell_roundoff_bound,
            bottom_before.common_datum_energy_j,
            wall.wall_storage_energy_j,
        ),
        *(
            _roundoff_aware_relative_residual(value, bound, state_scale)
            for residuals, bounds, before_cell in zip(
                upper_residuals,
                upper_roundoff_bounds_tuple,
                request.prior.tray.cells[1:],
                strict=True,
            )
            for value, bound, state_scale in zip(
                residuals,
                bounds,
                (
                    before_cell.conserved_state.totals.water_kg,
                    before_cell.conserved_state.totals.hexane_kg,
                    before_cell.conserved_state.totals.common_datum_energy_j,
                ),
                strict=True,
            )
        ),
    )
    return PDSystemLedger(
        packet_hexane_change_kg=packet_delta[1],
        cell_hexane_change_kg=cell_delta[1],
        external_hexane_to_system_kg=external_hexane,
        hexane_mass_residual_kg=hexane_residual,
        packet_water_change_kg=packet_delta[0],
        cell_water_change_kg=cell_delta[0],
        external_water_to_system_kg=external_water,
        water_mass_residual_kg=water_residual,
        packet_common_energy_change_j=packet_delta[2],
        cell_common_energy_change_j=cell_delta[2],
        thermal_external_common_energy_to_system_j=thermal_external_energy,
        q_other_s_external_energy_j=q_other_s_external_energy,
        external_common_energy_to_system_j=external_energy,
        common_energy_residual_j=energy_residual,
        wall_cell_storage_residual_j=wall_cell_residual,
        upper_cell_change_residuals=upper_residuals,
        upper_cell_change_binary64_roundoff_bounds=upper_roundoff_bounds_tuple,
        hexane_binary64_roundoff_bound_kg=hexane_roundoff_bound,
        water_binary64_roundoff_bound_kg=water_roundoff_bound,
        common_energy_binary64_roundoff_bound_j=energy_roundoff_bound,
        wall_cell_binary64_roundoff_bound_j=wall_cell_roundoff_bound,
        maximum_relative_residual=maximum,
        relative_limit=request.relative_limit,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PDKTrayStepRejection:
    attempt_id: str
    reason: str
    rollback_state: AcceptedPDKTrayState

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


_VALIDATED_PD_STEP_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedPDKTrayStep:
    request: PDKTrayStepRequest
    candidate: AcceptedPDKTrayState
    host_transaction: tray_host.ValidatedKCellTrayStep
    saturation_temperature_k: float
    packet_ledgers: tuple[PDPacketHeatLedger, ...]
    inter_k_face: PDInterKFaceLedger
    pressure_ledger: PDBypassPressureLedger
    plenum_ledger: PDPlenumLedger
    wall_ledger: PDSharedWallLedger
    system_ledger: PDSystemLedger
    meal_mechanical_work_deposition_ledger: meal_work.QSCMealMechanicalWorkDepositionLedger | None
    meal_mechanical_work_validation_context_digest: str | None
    codec_registry_digest: str
    auditor_identity_digest: str
    auditor_source_identity: str
    candidate_digest: str
    #: F4 (2026-09-02): which authority set ``saturation_temperature_k``, and
    #: what the CERTIFIED heteroazeotrope authority says at this tray pressure.
    #: ``boiling_temperature_excess_k`` is the measured gap between the two and
    #: is the wave-2 defect, carried on every step so it cannot go unnoticed.
    boiling_authority_id: str
    certified_co_boiling_temperature_k: float
    certified_vapour_water_mole_fraction: float
    boiling_temperature_excess_k: float
    #: WS-A (2026-09-03).  ``condensation_ledger`` is ``None`` whenever the
    #: capability was inert for this step (the declared contact fraction is
    #: 0.0, or the dome stream carried no water supersaturation to relieve).
    #: ``moisture_pickup`` is ALWAYS present: Q-PD-3 wants the meal's own
    #: moisture reported on every accepted PD step, not only on the steps that
    #: happened to condense something.
    condensation_ledger: PDContactCondensationLedger | None = None
    moisture_pickup: PDMoisturePickupLedger | None = None
    _seal: object = field(repr=False, compare=False)

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    industrial_operating_envelope_qualified: ClassVar[bool] = False
    finite_rpm_k_mixing_implemented: ClassVar[bool] = True
    inter_k_solid_heat_transfer_implemented: ClassVar[bool] = True
    shaft_work_energy_model_selected: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_PD_STEP_SEAL:
            raise TypeError("validated PD steps are issued only by the evaluator")
        if type(self.request) is not PDKTrayStepRequest:
            raise TypeError("validated PD request has a foreign exact type")
        PDKTrayStepRequest.__post_init__(self.request)
        bundle = self.request.meal_mechanical_work_deposition
        if bundle is None:
            if (
                self.meal_mechanical_work_deposition_ledger is not None
                or self.meal_mechanical_work_validation_context_digest is not None
                or self.system_ledger.q_other_s_external_energy_j != 0.0
            ):
                raise PDKDomeBypassStepError(
                    "zero-work PD transaction contains a mechanical-work artifact"
                )
        else:
            if (
                self.meal_mechanical_work_deposition_ledger != bundle.ledger
                or type(self.meal_mechanical_work_validation_context_digest) is not str
                or self.system_ledger.q_other_s_external_energy_j
                != bundle.ledger.emitted_packet_mechanical_work_j
            ):
                raise PDKDomeBypassStepError(
                    "work-bearing PD transaction and Q_other_s authorities disagree"
                )
            _require_sha256(
                "PD transaction meal-work validation context",
                self.meal_mechanical_work_validation_context_digest,
            )


def _validate_meal_mechanical_work_deposition(
    request: PDKTrayStepRequest,
    validation_context: meal_work.QSCMealMechanicalWorkValidationContext | None,
) -> meal_work.QSCMealMechanicalWorkDepositionBundle | None:
    bundle = request.meal_mechanical_work_deposition
    if bundle is None:
        if validation_context is not None:
            raise PDKDomeBypassConfigurationError(
                "a Q_other_s validation context is inadmissible without a deposition bundle"
            )
        return None
    if validation_context is None:
        raise PDKDomeBypassConfigurationError(
            "a meal mechanical-work deposition bundle requires caller-held validation pins"
        )
    try:
        meal_work.validate_qsc_meal_mechanical_work_deposition_bundle(
            bundle,
            accepted_tray=request.prior.tray,
            end_time_s=request.end_time_s,
            validation_context=validation_context,
            preverified_definition_digest=request.preverified_meal_bundle_digest,
        )
    except (TypeError, meal_work.QSCMealMechanicalWorkDepositionError) as error:
        raise PDKDomeBypassConfigurationError(
            f"meal mechanical-work deposition failed exact reconstruction: {error}"
        ) from error
    return bundle


def _validate_pd_k_planner(
    request: PDKTrayStepRequest,
    planner: tray_host.ManufacturedKRekeyPlanner,
) -> None:
    planner_type = type(planner)
    if planner_type not in (
        tray_host.NoKRekeyPlanner,
        k_scheduler.EngineeringPhysicalKRekeyPlanner,
    ):
        raise PDKDomeBypassConfigurationError(
            "the bounded PD tranche admits only the exact zero-motion or "
            "replay-validated engineering K-event planner"
        )
    if planner_type is tray_host.NoKRekeyPlanner:
        if (
            planner.planner_id != "manufactured-zero-k-rekey-v1"
            or planner.source_identity != "declared-zero-motion"
        ):
            raise PDKDomeBypassConfigurationError(
                "zero-motion planner ID or source differs from the canonical NoK declaration"
            )
        return
    try:
        k_scheduler.EngineeringPhysicalKRekeyPlanner.__post_init__(planner)
    except (TypeError, k_scheduler.EngineeringKEventSchedulerError) as error:
        raise PDKDomeBypassConfigurationError(
            f"engineering K planner failed deterministic replay validation: {error}"
        ) from error
    schedule_step = planner.step
    if schedule_step.prior.time_s != request.prior.tray.time_s:
        raise PDKDomeBypassConfigurationError(
            "engineering K planner prior time differs from the accepted PD tray"
        )
    if schedule_step.candidate.time_s != request.end_time_s:
        raise PDKDomeBypassConfigurationError(
            "engineering K planner endpoint differs from the PD tray step endpoint"
        )
    if schedule_step.prior.physical_tray_id != request.prior.tray.physical_tray_id:
        raise PDKDomeBypassConfigurationError(
            "engineering K planner belongs to another physical tray"
        )


def evaluate_pd_k_dome_bypass_tray_step(
    request: PDKTrayStepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner | None = None,
    meal_mechanical_work_validation_context: (
        meal_work.QSCMealMechanicalWorkValidationContext | None
    ) = None,
) -> ValidatedPDKTrayStep | PDKTrayStepRejection:
    """Evaluate one PD wall/packet/plenum step without mutating accepted state."""

    if type(request) is not PDKTrayStepRequest:
        raise TypeError("PD integration requires an exact PDKTrayStepRequest")
    try:
        if type(codec) is not high_loading.HighLoadingV1PacketAdapter and (
            type(codec) is not dual_variant.EngineeringDualVariantPacketAdapter
        ):
            raise PDKDomeBypassConfigurationError(
                "the PD limit requires the exact detached tag-1 high-loading adapter "
                "or the owner-ruled dual-variant adapter"
            )
        if k_planner is None:
            k_planner = tray_host.NoKRekeyPlanner()
        _validate_pd_k_planner(request, k_planner)
        deposition_bundle = _validate_meal_mechanical_work_deposition(
            request,
            meal_mechanical_work_validation_context,
        )
        codec.require_host_binding(
            configuration_digest=request.prior.tray.configuration_digest,
            energy_datum_id=request.prior.tray.energy_datum_id,
            pressure_pa=request.prior.tray.pressure_pa,
        )
        if codec.codec_registry_digest != request.prior.tray.codec_registry_digest:
            raise PDKDomeBypassConfigurationError("PD codec violates the host registry pin")
        if codec.auditor_identity_digest != request.prior.tray.auditor_identity_digest:
            raise PDKDomeBypassConfigurationError("PD codec violates the host auditor pin")
        codec.component_datum_adapter.require_gas_properties(request.model.properties)
        # CS-1a (2026-09-04): the cold-surface term needs the SURFACE it
        # condenses on and the caloric room that surface still has, so the
        # layer reads and the inter-K face heat are evaluated HERE instead of
        # below.  Every one of them is a pure function of the ACCEPTED prior
        # state -- no mutation, no host call -- and none is recomputed later,
        # so the arithmetic reaching every downstream expression is unmoved.
        if len(request.prior.tray.cells) != PD_K2_LAYER_COUNT:
            raise PDKDomeBypassConfigurationError("the bounded PD thermal tranche is exactly K=2")
        bottom_packets = _layer_packets(request.prior.tray, 1)
        upper_packets = _layer_packets(request.prior.tray, 2)
        lower_temperature = _layer_temperature_k(bottom_packets, codec=codec)
        upper_temperature = _layer_temperature_k(upper_packets, codec=codec)
        face_area = math.pi * (0.5 * request.model.geometry.tray.diameter_m) ** 2
        face_distance = 0.5 * request.model.geometry.tray.loaded_depth_m
        conductance = PD_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K * face_area / face_distance
        temperature_difference = lower_temperature - upper_temperature
        face_heat = conductance * temperature_difference * request.duration_s
        # CS-1a REVISION 2026-09-04 (verifier issue 3).  The first-stop cap
        # must subtract EVERY declared energy this step already deposits into
        # the receiving surface, not only the inter-K face heat.  Declared
        # MECHANICAL (shaft) work is deposited into top-layer packets too --
        # measured on this very tray family, finding F-CS1-1 -- and it enters
        # the packet energy at ``total_declared_packet_energy_input_j``
        # alongside the declared heat.  Left out of the cap, latent plus shaft
        # work could carry the receiving packet PAST the stop temperature the
        # cap exists to protect.
        #
        # CS-1a REV 2 2026-09-04 (finding F-CS1-2).  The cap is applied per
        # PACKET, so the declared energy it subtracts is per packet too: this
        # packet's share of the inter-K face heat (from the SAME exact-residual
        # allocator that deposits it below) plus this packet's OWN declared
        # mechanical work, read from the already-validated deposition bundle,
        # so no new authority enters here.
        upper_face_heat_by_packet = _allocate_layer_heat(upper_packets, face_heat)
        upper_layer_declared_energy_by_packet = {
            packet.packet_id: math.fsum(
                (
                    upper_face_heat_by_packet[packet.packet_id],
                    (
                        0.0
                        if deposition_bundle is None
                        else deposition_bundle.packet_row(
                            packet.packet_id
                        ).declared_mechanical_work_j
                    ),
                )
            )
            for packet in upper_packets
        }
        condensation = _contact_limited_condensation(
            request.bypass_inlet,
            model=request.model,
            duration_s=request.duration_s,
            codec=codec,
            upper_packets=upper_packets,
            upper_layer_declared_energy_by_packet=(
                upper_layer_declared_energy_by_packet
            ),
        )
        if PD_STEAM_CONDENSING_LEVEL0_ARMED:
            # D12 steam condensing level 0 (study-lane knob): the rising vapour
            # condenses on EVERY packet of the tray, each on its own slice.  The
            # CS-1a surface term must be inert; two condensation laws on one
            # stream would count the same vapour twice.
            if condensation is not None:
                raise PDKDomeBypassStepError(
                    "steam condensing level 0 cannot run beside the armed CS-1a "
                    "surface term; refuse, do not add"
                )
            steam_bottom_dry_matter = math.fsum(
                packet.weighted_totals.totals[0] for packet in bottom_packets
            )
            steam_lower_face_heat_by_packet = _allocate_layer_heat(bottom_packets, -face_heat)
            steam_declared_energy_by_packet = dict(upper_layer_declared_energy_by_packet)
            for packet in bottom_packets:
                steam_wall_heat = (
                    request.model.wall.floor_to_bottom_solid_ua_w_k
                    * (packet.weighted_totals.totals[0] / steam_bottom_dry_matter)
                    * (
                        request.prior.wall_temperature_k
                        - codec.temperature_from_payload(packet.entry.payload_bytes)
                    )
                    * request.duration_s
                )
                steam_declared_energy_by_packet[packet.packet_id] = math.fsum(
                    (
                        steam_wall_heat,
                        steam_lower_face_heat_by_packet[packet.packet_id],
                        (
                            0.0
                            if deposition_bundle is None
                            else deposition_bundle.packet_row(
                                packet.packet_id
                            ).declared_mechanical_work_j
                        ),
                    )
                )
            condensation = _steam_condensing_level0(
                request.bypass_inlet,
                model=request.model,
                duration_s=request.duration_s,
                codec=codec,
                packets=request.prior.tray.packets,
                declared_energy_by_packet=steam_declared_energy_by_packet,
            )
        # WS-A, NARROWED 2026-09-04: what the no-flash guard admits is not the
        # ARMED CAPABILITY but the ledger the capability actually produced,
        # and only up to that ledger's own measured numbers.  Both tuples are
        # EMPTY whenever ``condensation is None`` - the declared production
        # default, and equally every armed state whose term does not fire - so
        # all three guard sites below are then the pre-WS-A calls byte for
        # byte.
        #
        #  * at the BYPASS INLET, the stream the ledger was built FROM, the
        #    admitted excess is that ledger's own driving excess: the guard
        #    admits exactly the supersaturation the term priced, and refuses
        #    any stream more supersaturated than the ledger says it saw.
        #  * at the two PLENUM streams the admitted excess is the ledger's own
        #    ``residual_supersaturation_pa`` - the supersaturation this step
        #    DECLARES it did not relieve.  A plenum stream supersaturated
        #    beyond that declared residual is still refused.
        inlet_admitted_supersaturation_pa: tuple[tuple[str, float], ...] = (
            ()
            if condensation is None or condensation.partial_pressure_excess_pa <= 0.0
            else (("water", condensation.partial_pressure_excess_pa),)
        )
        plenum_admitted_supersaturation_pa: tuple[tuple[str, float], ...] = (
            ()
            if condensation is None or condensation.residual_supersaturation_pa <= 0.0
            else (("water", condensation.residual_supersaturation_pa),)
        )
        if PD_FOG_ONLY_LEVEL0_ARMED:
            # D12 fog-only level 0 (study-lane knob), the lower end of the bracket:
            # nothing condenses on the meal, so no condensation term may run beside
            # it; each guarded stream is admitted exactly to the water excess its own
            # fog ledger measured (the fog carried by the vapour).
            if PD_STEAM_CONDENSING_LEVEL0_ARMED or condensation is not None:
                raise PDFogOnlyLevel0Refusal(
                    "fog-only level 0 deposits nothing on the meal and cannot run beside "
                    "a condensation term on it; refuse, do not add"
                )
            fog_ledgers: dict[str, object] = {}
            PD_FOG_ONLY_LEVEL0_LAST_LEDGERS[id(request)] = fog_ledgers
            plenum_admitted_supersaturation_pa = _fog_only_level0_admitted(
                request.prior.plenum_outlet,
                label="accepted plenum outlet",
                properties=request.model.properties,
                datum_adapter=codec.component_datum_adapter,
                relative_limit=request.relative_limit,
                ledgers=fog_ledgers,
            )
            inlet_admitted_supersaturation_pa = _fog_only_level0_admitted(
                request.bypass_inlet,
                label="bypass inlet",
                properties=request.model.properties,
                datum_adapter=codec.component_datum_adapter,
                relative_limit=request.relative_limit,
                ledgers=fog_ledgers,
            )
        _require_single_vapour_phase(
            request.prior.plenum_outlet,
            label="accepted plenum outlet",
            admitted_supersaturation_pa=plenum_admitted_supersaturation_pa,
        )
        _require_single_vapour_phase(
            request.bypass_inlet,
            label="bypass inlet",
            admitted_supersaturation_pa=inlet_admitted_supersaturation_pa,
        )

        # F4 (2026-09-02): the PD boiling temperature is the CERTIFIED
        # heteroazeotrope at the fixed tray pressure, not the pure n-hexane
        # saturation root.  The pure root is retained as the single-liquid
        # fallback for a tray whose packets carry no external free water at
        # all - there is then no second liquid and no heteroazeotrope.
        # Guarded: while the dispatch is the declared pre-F4 branch this decode
        # sweep would be pure overhead on the speed-programme path.
        tray_external_water_kg = (
            math.fsum(
                codec.decode(packet.entry.payload_bytes).inventory.external_water_kg
                * packet.entry.weight.representative_particles
                for packet in request.prior.tray.packets
            )
            if PD_CERTIFIED_CO_BOILING_ENABLED
            else 0.0
        )
        # The certified authority is consulted UNCONDITIONALLY, so its own
        # refusals (pressure outside the G5-grounded band) still bind and its
        # numbers are always available to the step record.
        certified_co_boiling = _pd_co_boiling_state(request.prior.tray.pressure_pa)
        pure_hexane_saturation = _hexane_saturation_temperature_k(
            request.prior.tray.pressure_pa,
            codec=codec,
        )
        if PD_CERTIFIED_CO_BOILING_ENABLED and tray_external_water_kg > 0.0:
            co_boiling = certified_co_boiling
            if PD_FLASH_LEVEL0_ARMED and PD_FLASH_LEVEL0_WATER_SHARE_FACTOR != 1.0:
                # D12 level 0 SOURCE BOUND (design site 5): the knob module's
                # fixed point withholds part of the co-boiled water at the
                # source; the lock temperature is unchanged.
                co_boiling = _level0_source_bound_co_boiling_state(
                    co_boiling, PD_FLASH_LEVEL0_WATER_SHARE_FACTOR
                )
            saturation_temperature = co_boiling.temperature_k
        else:
            co_boiling = None
            saturation_temperature = pure_hexane_saturation
        duration = request.duration_s
        # WS-A / CS-1a REV 2: the dome vapour touches the bed's TOP free
        # surface, so both condensates land on the top K layer - and each
        # packet receives exactly what ITS OWN cold surface earned, from that
        # packet's row (F-CS1-2).  Rev 1 shared a layer total by
        # representative-particle count, which handed a warm resident packet
        # the condensate a cold feed packet had produced.  The rows sum to the
        # ledger's totals by construction (the totals ARE their fsums), so the
        # plenum and system closures are unmoved.  Empty at the declared
        # default, where the term returns None.
        condensate_by_packet: dict[str, float] = {}
        condensate_hexane_by_packet: dict[str, float] = {}
        condensation_energy_by_packet: dict[str, float] = {}
        # REV 2: the endpoint a packet may relax onto is ITS OWN first stop
        # temperature, not the layer's minimum over all condensing packets.  A
        # packet that condensed nothing gets no endpoint at all, exactly as
        # before this term existed, so the envelope it is checked against is
        # not widened on another packet's account.
        condensation_endpoint_by_packet: dict[str, float] = {}
        if condensation is not None:
            condensation_endpoint_by_packet = {
                row.packet_id: row.stop_temperature_k
                for row in condensation.packet_rows
                if row.stop_component and row.condensation_common_energy_j > 0.0
            }
            condensate_by_packet = {
                row.packet_id: row.water_condensed_kg
                for row in condensation.packet_rows
                if row.water_condensed_kg > 0.0
            }
            condensate_hexane_by_packet = {
                row.packet_id: row.hexane_condensed_kg
                for row in condensation.packet_rows
                if row.hexane_condensed_kg > 0.0
            }
            condensation_energy_by_packet = {
                row.packet_id: row.condensation_common_energy_j
                for row in condensation.packet_rows
                if row.condensation_common_energy_j > 0.0
            }
        bottom_dry_matter = math.fsum(packet.weighted_totals.totals[0] for packet in bottom_packets)
        if bottom_dry_matter <= 0.0:
            raise PDKDomeBypassStepError("bottom K layer has no positive dry-matter holdup")

        # The crossing guard: the shipped form, or the uniform-layer form under
        # ``DTDC_PD_INTERK_CAPACITY_UNIFORM_LAYER`` (ruling packet option (iv)).
        # R3, the owner's B3 ruling of 2026-08-31 (study-lane knob, default OFF): with
        # ``DTDC_PD_INTERK_CAPACITY_ANALYTIC`` exactly "1" the capacity is the analytic
        # layer heat capacity times the separation
        # (``study_lane/pd_interk_capacity_analytic.py``), same tolerance, same refusal
        # text with both capacities printed.  Off, the ``elif`` is the former guard.
        if os.environ.get("DTDC_PD_INTERK_CAPACITY_ANALYTIC") == "1":
            from .study_lane import pd_interk_capacity_analytic as interk_analytic

            interk_analytic.require_face_heat_within_analytic_capacity(
                request,
                bottom_packets,
                upper_packets,
                lower_temperature_k=lower_temperature,
                upper_temperature_k=upper_temperature,
                face_heat_j=face_heat,
                codec=codec,
            )
        elif _inter_k_crossing_refused(
            request,
            bottom_packets,
            upper_packets,
            lower_temperature_k=lower_temperature,
            upper_temperature_k=upper_temperature,
            face_heat_j=face_heat,
            codec=codec,
        ):
            raise PDKDomeBypassStepError(
                "explicit inter-K conduction would cross the old layer equilibrium; "
                "reject and retry at a shorter endpoint"
            )
        lower_inter_k = _allocate_layer_heat(bottom_packets, -face_heat)
        upper_inter_k = _allocate_layer_heat(upper_packets, face_heat)
        inter_k_by_packet = {**lower_inter_k, **upper_inter_k}

        entries: dict[str, tray_host.PacketPlanEntry] = {}
        thermal_guard_entries: dict[str, tray_host.PacketPlanEntry] = {}
        order_alternative_entries: dict[str, tray_host.PacketPlanEntry] = {}
        order_alternative_refusals: list[str] = []
        order_zero_face_inputs: dict[str, tuple[float, bool]] = {}
        packet_ledgers: list[PDPacketHeatLedger] = []
        for packet in request.prior.tray.packets:
            if packet.owner_key.vertical_layer_id.value == 1:
                share = packet.weighted_totals.totals[0] / bottom_dry_matter
                packet_temperature = codec.temperature_from_payload(packet.entry.payload_bytes)
                declared_heat = (
                    request.model.wall.floor_to_bottom_solid_ua_w_k
                    * share
                    * (request.prior.wall_temperature_k - packet_temperature)
                    * duration
                )
            else:
                declared_heat = 0.0
            mechanical_work = (
                0.0
                if deposition_bundle is None
                else deposition_bundle.packet_row(packet.packet_id).declared_mechanical_work_j
            )
            thermal_guard_entry: tray_host.PacketPlanEntry | None = None
            thermal_guard_ledger: PDPacketHeatLedger | None = None
            if mechanical_work > 0.0:
                thermal_guard_entry, thermal_guard_ledger = _advance_one_packet(
                    packet,
                    declared_wall_heat_j=declared_heat,
                    declared_inter_k_heat_j=inter_k_by_packet[packet.packet_id],
                    declared_mechanical_work_j=0.0,
                    old_wall_temperature_k=request.prior.wall_temperature_k,
                    adjacent_old_temperature_k=(
                        upper_temperature
                        if packet.owner_key.vertical_layer_id.value == 1
                        else lower_temperature
                    ),
                    saturation_temperature_k=saturation_temperature,
                    tray_pressure_pa=request.prior.tray.pressure_pa,
                    properties=request.model.properties,
                    codec=codec,
                    relative_limit=request.relative_limit,
                    co_boiling=co_boiling,
                    declared_condensate_water_kg=condensate_by_packet.get(
                        packet.packet_id, 0.0
                    ),
                    declared_condensate_hexane_kg=condensate_hexane_by_packet.get(
                        packet.packet_id, 0.0
                    ),
                    declared_condensation_energy_j=condensation_energy_by_packet.get(
                        packet.packet_id, 0.0
                    ),
                    declared_condensation_equilibrium_temperature_k=(
                        condensation_endpoint_by_packet.get(packet.packet_id)
                    ),
                    lock_band_boils=PD_INTERK_ORDER_LOCK_CONSISTENT_ARMED,
                    exploration_call_kind="thermal_guard",
                )
                if PD_INTERK_ORDER_LOCK_OBSERVER is not None:
                    # Report only: the other lock rule's thermal-only entry, so the
                    # observer sees both forms of the order guard.  Its refusal is
                    # reported, never raised.
                    try:
                        order_alternative_entries[packet.packet_id] = _advance_one_packet(
                            packet,
                            declared_wall_heat_j=declared_heat,
                            declared_inter_k_heat_j=inter_k_by_packet[packet.packet_id],
                            declared_mechanical_work_j=0.0,
                            old_wall_temperature_k=request.prior.wall_temperature_k,
                            adjacent_old_temperature_k=(
                                upper_temperature
                                if packet.owner_key.vertical_layer_id.value == 1
                                else lower_temperature
                            ),
                            saturation_temperature_k=saturation_temperature,
                            tray_pressure_pa=request.prior.tray.pressure_pa,
                            properties=request.model.properties,
                            codec=codec,
                            relative_limit=request.relative_limit,
                            co_boiling=co_boiling,
                            declared_condensate_water_kg=condensate_by_packet.get(
                                packet.packet_id, 0.0
                            ),
                            declared_condensate_hexane_kg=condensate_hexane_by_packet.get(
                                packet.packet_id, 0.0
                            ),
                            declared_condensation_energy_j=condensation_energy_by_packet.get(
                                packet.packet_id, 0.0
                            ),
                            declared_condensation_equilibrium_temperature_k=(
                                condensation_endpoint_by_packet.get(packet.packet_id)
                            ),
                            lock_band_boils=not PD_INTERK_ORDER_LOCK_CONSISTENT_ARMED,
                        )[0]
                    except Exception as error:  # noqa: BLE001 - report only
                        order_alternative_refusals.append(
                            f"{packet.packet_id}: {type(error).__name__}: {error}"
                        )
            entry, ledger = _advance_one_packet(
                packet,
                declared_wall_heat_j=declared_heat,
                declared_inter_k_heat_j=inter_k_by_packet[packet.packet_id],
                declared_mechanical_work_j=mechanical_work,
                old_wall_temperature_k=request.prior.wall_temperature_k,
                adjacent_old_temperature_k=(
                    upper_temperature
                    if packet.owner_key.vertical_layer_id.value == 1
                    else lower_temperature
                ),
                saturation_temperature_k=saturation_temperature,
                tray_pressure_pa=request.prior.tray.pressure_pa,
                properties=request.model.properties,
                codec=codec,
                relative_limit=request.relative_limit,
                co_boiling=co_boiling,
                declared_condensate_water_kg=condensate_by_packet.get(packet.packet_id, 0.0),
                declared_condensate_hexane_kg=condensate_hexane_by_packet.get(
                    packet.packet_id, 0.0
                ),
                declared_condensation_energy_j=condensation_energy_by_packet.get(
                    packet.packet_id, 0.0
                ),
                declared_condensation_equilibrium_temperature_k=(
                    condensation_endpoint_by_packet.get(packet.packet_id)
                ),
                thermal_guard_temperature_after_k=(
                    None
                    if thermal_guard_ledger is None
                    else thermal_guard_ledger.temperature_after_k
                ),
                thermal_guard_generated_hexane_kg=(
                    None
                    if thermal_guard_ledger is None
                    else thermal_guard_ledger.generated_hexane_kg
                ),
                exploration_call_kind="in_force",
            )
            entries[packet.packet_id] = entry
            thermal_guard_entries[packet.packet_id] = (
                entry if thermal_guard_entry is None else thermal_guard_entry
            )
            packet_ledgers.append(ledger)
            if PD_INTERK_ORDER_CONDUCTION_ONLY_ARMED:
                order_zero_face_inputs[packet.packet_id] = (declared_heat, mechanical_work > 0.0)

        lower_allocation = math.fsum(lower_inter_k.values())
        upper_allocation = math.fsum(upper_inter_k.values())
        allocation_residual = math.fsum((lower_allocation, upper_allocation))
        lower_temperature_after = _layer_temperature_after_k(
            bottom_packets,
            entries,
            codec=codec,
        )
        upper_temperature_after = _layer_temperature_after_k(
            upper_packets,
            entries,
            codec=codec,
        )
        lower_thermal_only_temperature_after = _layer_temperature_after_k(
            bottom_packets,
            thermal_guard_entries,
            codec=codec,
        )
        upper_thermal_only_temperature_after = _layer_temperature_after_k(
            upper_packets,
            thermal_guard_entries,
            codec=codec,
        )
        inter_k_crossed = _explicit_pair_crossed_equilibrium(
            first_before_k=lower_temperature,
            second_before_k=upper_temperature,
            first_after_k=lower_thermal_only_temperature_after,
            second_after_k=upper_thermal_only_temperature_after,
        )
        if PD_INTERK_ORDER_LOCK_OBSERVER is not None:
            _report_order_lock_forms(
                request,
                bottom_packets,
                upper_packets,
                thermal_guard_entries=thermal_guard_entries,
                alternative_entries=order_alternative_entries,
                alternative_refusals=order_alternative_refusals,
                lower_temperature_k=lower_temperature,
                upper_temperature_k=upper_temperature,
                lower_after_in_force_k=lower_thermal_only_temperature_after,
                upper_after_in_force_k=upper_thermal_only_temperature_after,
                crossed_in_force=inter_k_crossed,
                saturation_temperature_k=saturation_temperature,
                codec=codec,
            )
        if PD_INTERK_ORDER_CONDUCTION_ONLY_ARMED:
            # ``DTDC_PD_INTERK_ORDER_CONDUCTION_ONLY``: the uniform layers after the
            # face heat alone decide; the all-sources reading is reported.

            def zero_face_heat_entries() -> dict[str, tray_host.PacketPlanEntry]:
                counterfactual: dict[str, tray_host.PacketPlanEntry] = {}
                for packet in request.prior.tray.packets:
                    declared_heat_j, work_bearing = order_zero_face_inputs[packet.packet_id]
                    counterfactual[packet.packet_id] = _advance_one_packet(
                        packet,
                        declared_wall_heat_j=declared_heat_j,
                        declared_inter_k_heat_j=0.0,
                        declared_mechanical_work_j=0.0,
                        old_wall_temperature_k=request.prior.wall_temperature_k,
                        adjacent_old_temperature_k=(
                            upper_temperature
                            if packet.owner_key.vertical_layer_id.value == 1
                            else lower_temperature
                        ),
                        saturation_temperature_k=saturation_temperature,
                        tray_pressure_pa=request.prior.tray.pressure_pa,
                        properties=request.model.properties,
                        codec=codec,
                        relative_limit=request.relative_limit,
                        co_boiling=co_boiling,
                        declared_condensate_water_kg=condensate_by_packet.get(
                            packet.packet_id, 0.0
                        ),
                        declared_condensate_hexane_kg=condensate_hexane_by_packet.get(
                            packet.packet_id, 0.0
                        ),
                        declared_condensation_energy_j=condensation_energy_by_packet.get(
                            packet.packet_id, 0.0
                        ),
                        declared_condensation_equilibrium_temperature_k=(
                            condensation_endpoint_by_packet.get(packet.packet_id)
                        ),
                        lock_band_boils=work_bearing and PD_INTERK_ORDER_LOCK_CONSISTENT_ARMED,
                    )[0]
                return counterfactual

            # EXPLORATION MODE F2, site 2 (record GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13):
            # with the hook bound, a room refusal of either layer is logged with both layers'
            # rooms and read as crossed (the verdict skipped); otherwise the former call.
            inter_k_crossed = (
                _conduction_only_order_crossed
                if PD_EXPLORATION_HOOK is None
                else PD_EXPLORATION_HOOK.conduction_only_order_crossed
            )(
                request,
                bottom_packets,
                upper_packets,
                lower_temperature_k=lower_temperature,
                upper_temperature_k=upper_temperature,
                face_heat_j=face_heat,
                all_sources_lower_after_k=lower_thermal_only_temperature_after,
                all_sources_upper_after_k=upper_thermal_only_temperature_after,
                all_sources_crossed=inter_k_crossed,
                saturation_temperature_k=saturation_temperature,
                zero_face_heat_entries=zero_face_heat_entries,
                codec=codec,
            )
        # EXPLORATION MODE F2, the order guard: the hook logs the crossing (once per PD step)
        # and admits it; inter_k_crossed stays True for the face ledger below.
        if inter_k_crossed and not (
            PD_EXPLORATION_HOOK is not None
            and PD_EXPLORATION_HOOK.admit_order_crossing(
                request,
                lower_temperature_k=lower_temperature,
                upper_temperature_k=upper_temperature,
                lower_after_k=lower_thermal_only_temperature_after,
                upper_after_k=upper_thermal_only_temperature_after,
                face_heat_j=face_heat,
            )
        ):
            raise PDKDomeBypassStepError(
                "explicit inter-K conduction reversed the accepted layer temperature order; "
                "reject and retry at a shorter endpoint"
            )
        entropy_production = (
            conductance
            * duration
            * temperature_difference
            * temperature_difference
            / (lower_temperature * upper_temperature)
        )
        inter_k_maximum = _roundoff_aware_relative_residual(
            allocation_residual,
            _ulp_sum(
                face_heat,
                lower_allocation,
                upper_allocation,
                allocation_residual,
            ),
            face_heat,
            lower_allocation,
            upper_allocation,
        )
        inter_k_face = PDInterKFaceLedger(
            lower_vertical_layer=1,
            upper_vertical_layer=2,
            effective_conductivity_w_m_k=PD_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K,
            face_area_m2=face_area,
            gradient_distance_m=face_distance,
            conductance_w_k=conductance,
            lower_temperature_before_k=lower_temperature,
            upper_temperature_before_k=upper_temperature,
            lower_temperature_after_k=lower_temperature_after,
            upper_temperature_after_k=upper_temperature_after,
            lower_thermal_only_temperature_after_k=(lower_thermal_only_temperature_after),
            upper_thermal_only_temperature_after_k=(upper_thermal_only_temperature_after),
            heat_lower_to_upper_j=face_heat,
            lower_packet_allocation_j=lower_allocation,
            upper_packet_allocation_j=upper_allocation,
            allocation_residual_j=allocation_residual,
            entropy_production_j_k=entropy_production,
            no_explicit_equilibrium_crossing=not inter_k_crossed,
            maximum_relative_residual=inter_k_maximum,
            relative_limit=request.relative_limit,
        )
        # EXPLORATION MODE F2: admitted only when its ONLY failing component is a crossing the
        # hook logged in this step (allocation residual and entropy still pass).
        if not inter_k_face.passed and not (
            PD_EXPLORATION_HOOK is not None
            and PD_EXPLORATION_HOOK.admit_inter_k_face(request, inter_k_face)
        ):
            raise PDKDomeBypassStepError("inter-K conduction ledger failed")

        solid_heat = math.fsum(item.declared_wall_heat_j for item in packet_ledgers)
        generated_hexane = math.fsum(item.generated_hexane_kg for item in packet_ledgers)
        generated_water = math.fsum(item.generated_water_kg for item in packet_ledgers)
        generated_energy = math.fsum(
            item.generated_hexane_common_energy_j for item in packet_ledgers
        )
        wall_to_gas = (
            request.model.wall.wall_to_bypass_gas_ua_w_k
            * (request.prior.wall_temperature_k - request.bypass_inlet.temperature_k)
            * duration
        )
        if type(request.model.properties) is cc.DeclaredGasProperties:
            bypass_sensible_to_old_wall = (
                request.bypass_inlet.molar_flow_mol_s
                * request.model.properties.molar_heat_capacity_j_mol_k
                * abs(request.prior.wall_temperature_k - request.bypass_inlet.temperature_k)
                * duration
            )
        else:
            # CELL-02e-pd native arm: the exact enthalpy change to wall/gas
            # equilibrium from the frozen laws (no linearization).
            bypass_sensible_to_old_wall = (
                request.bypass_inlet.molar_flow_mol_s
                * abs(
                    native_caloric.mixture_enthalpy(
                        request.prior.wall_temperature_k,
                        request.bypass_inlet.pressure_pa,
                        request.bypass_inlet.hexane_mole_fraction,
                    ).value_j_mol
                    - native_caloric.mixture_enthalpy(
                        request.bypass_inlet.temperature_k,
                        request.bypass_inlet.pressure_pa,
                        request.bypass_inlet.hexane_mole_fraction,
                    ).value_j_mol
                )
                * duration
            )
        if abs(wall_to_gas) > math.fsum(
            (
                bypass_sensible_to_old_wall,
                _ulp_sum(wall_to_gas, bypass_sensible_to_old_wall),
            )
        ):
            raise PDKDomeBypassStepError(
                "explicit wall-to-bypass UA would transfer past the old wall/gas "
                "equilibrium; reject, do not clamp"
            )
        ambient_loss = (
            request.model.wall.ambient_ua_w_k
            * (request.prior.wall_temperature_k - request.model.wall.ambient_temperature_k)
            * duration
        )
        accepted_duty = request.accepted_indirect_duty_w * duration
        wall_storage = math.fsum((accepted_duty, -solid_heat, -wall_to_gas, -ambient_loss))
        new_wall_temperature = math.fsum(
            (
                request.prior.wall_temperature_k,
                wall_storage / request.model.wall.thermal_capacity_j_k,
            )
        )
        if not math.isfinite(new_wall_temperature) or new_wall_temperature <= 0.0:
            raise PDKDomeBypassStepError("shared-wall update left the positive-temperature domain")
        if type(request.model.properties) is cc.DeclaredGasProperties:
            bypass_wall_only_temperature = math.fsum(
                (
                    request.bypass_inlet.temperature_k,
                    wall_to_gas
                    / (
                        request.bypass_inlet.molar_flow_mol_s
                        * request.model.properties.molar_heat_capacity_j_mol_k
                        * duration
                    ),
                )
            )
        else:
            try:
                bypass_wall_only_temperature = native_caloric.temperature_from_mixture_enthalpy(
                    native_caloric.mixture_enthalpy(
                        request.bypass_inlet.temperature_k,
                        request.bypass_inlet.pressure_pa,
                        request.bypass_inlet.hexane_mole_fraction,
                    ).value_j_mol
                    + wall_to_gas / (request.bypass_inlet.molar_flow_mol_s * duration),
                    request.bypass_inlet.pressure_pa,
                    request.bypass_inlet.hexane_mole_fraction,
                )
            except ValueError as error:
                raise PDKDomeBypassStepError(
                    f"native wall-only bypass inverse caloric refused: {error}"
                ) from error
        wall_crossed = _explicit_pair_crossed_equilibrium(
            first_before_k=request.prior.wall_temperature_k,
            second_before_k=request.bypass_inlet.temperature_k,
            first_after_k=new_wall_temperature,
            second_after_k=bypass_wall_only_temperature,
        )
        for packet_ledger in packet_ledgers:
            if packet_ledger.declared_wall_heat_j == 0.0:
                continue
            wall_crossed = wall_crossed or _explicit_pair_crossed_equilibrium(
                first_before_k=request.prior.wall_temperature_k,
                second_before_k=packet_ledger.temperature_before_k,
                first_after_k=new_wall_temperature,
                second_after_k=packet_ledger.thermal_only_temperature_after_k,
            )
        if request.model.wall.ambient_ua_w_k != 0.0:
            wall_crossed = wall_crossed or _explicit_pair_crossed_equilibrium(
                first_before_k=request.prior.wall_temperature_k,
                second_before_k=request.model.wall.ambient_temperature_k,
                first_after_k=new_wall_temperature,
                second_after_k=request.model.wall.ambient_temperature_k,
            )
        if wall_crossed:
            raise PDKDomeBypassStepError(
                "explicit shared-old-wall update reversed a coupled temperature order; "
                "reject and retry at a shorter endpoint"
            )
        wall_residual = math.fsum(
            (wall_storage, -accepted_duty, solid_heat, wall_to_gas, ambient_loss)
        )
        wall_maximum = _relative_residual(
            wall_residual,
            wall_storage,
            accepted_duty,
            solid_heat,
            wall_to_gas,
            ambient_loss,
        )
        wall_ledger = PDSharedWallLedger(
            old_wall_temperature_k=request.prior.wall_temperature_k,
            new_wall_temperature_k=new_wall_temperature,
            accepted_indirect_duty_j=accepted_duty,
            wall_to_bottom_solid_energy_j=solid_heat,
            wall_to_bypass_gas_energy_j=wall_to_gas,
            ambient_loss_energy_j=ambient_loss,
            wall_storage_energy_j=wall_storage,
            common_energy_residual_j=wall_residual,
            no_explicit_equilibrium_crossing=not wall_crossed,
            maximum_relative_residual=wall_maximum,
            relative_limit=request.relative_limit,
        )
        if not wall_ledger.passed:
            raise PDKDomeBypassStepError("shared old-wall common-energy ledger failed")

        outlet, plenum_ledger = _plenum_outlet_and_ledger(
            request,
            saturation_temperature_k=saturation_temperature,
            generated_hexane_kg=generated_hexane,
            generated_water_kg=generated_water,
            generated_hexane_energy_j=generated_energy,
            co_boiling=co_boiling,
            wall_to_bypass_gas_energy_j=wall_to_gas,
            codec=codec,
            condensation=condensation,
        )
        if PD_FLASH_LEVEL0_ARMED and not PD_FOG_ONLY_LEVEL0_ARMED:
            # D12 level 0: the candidate outlet the source bound reads.
            PD_FLASH_LEVEL0_LAST_OUTLET[id(request)] = (outlet, plenum_ledger, co_boiling)
        if PD_FOG_ONLY_LEVEL0_ARMED:
            # D12 fog-only level 0: the candidate outlet's fog, admitted exactly; the
            # source bound is not consulted (the fog takes the excess, nothing is
            # withheld at the source).
            plenum_admitted_supersaturation_pa = _fog_only_level0_admitted(
                outlet,
                label="candidate plenum outlet",
                properties=request.model.properties,
                datum_adapter=codec.component_datum_adapter,
                relative_limit=request.relative_limit,
                ledgers=fog_ledgers,
            )
        _require_single_vapour_phase(
            outlet,
            label="candidate plenum outlet",
            admitted_supersaturation_pa=plenum_admitted_supersaturation_pa,
        )
        pressure_residual = math.fsum(
            (
                request.bypass_inlet.pressure_pa,
                -request.prescribed_bypass_pressure_drop_pa,
                -request.prior.tray.pressure_pa,
            )
        )
        pressure_ledger = PDBypassPressureLedger(
            bypass_inlet_pressure_pa=request.bypass_inlet.pressure_pa,
            prescribed_pressure_drop_pa=request.prescribed_bypass_pressure_drop_pa,
            plenum_pressure_pa=request.prior.tray.pressure_pa,
            pressure_residual_pa=pressure_residual,
        )
        if not pressure_ledger.passed:
            raise PDKDomeBypassStepError("prescribed bypass/plenum pressure algebra failed")

        external_transfers_list: list[tray_host.CellExternalTransfer] = []
        for index, _ in enumerate(request.prior.tray.cells, start=1):
            layer_ledgers = tuple(item for item in packet_ledgers if item.vertical_layer == index)
            layer_generated = math.fsum(item.generated_hexane_kg for item in layer_ledgers)
            layer_generated_water = math.fsum(item.generated_water_kg for item in layer_ledgers)
            layer_condensate = math.fsum(item.condensed_water_kg for item in layer_ledgers)
            layer_condensate_hexane = math.fsum(
                item.condensed_hexane_kg for item in layer_ledgers
            )
            layer_packet_energy = math.fsum(
                item.packet_common_energy_change_j for item in layer_ledgers
            )
            external_transfers_list.append(
                tray_host.CellExternalTransfer(
                    # F4: the co-boiled water leaves the packets exactly as the
                    # n-hexane does, so the host's external water row carries it.
                    # WS-A adds the ARRIVING condensate on the same row, with the
                    # opposite sign; it is exactly 0.0 at the declared default.
                    water_to_cell_kg=math.fsum((-layer_generated_water, layer_condensate)),
                    # CS-1a: the arriving hexane condensate rides the hexane
                    # row exactly as the water condensate rides the water row.
                    # Exactly 0.0 at the declared default, where math.fsum of
                    # one addend reproduces the pre-CS-1 unary minus.
                    hexane_to_cell_kg=(
                        -layer_generated
                        if layer_condensate_hexane == 0.0
                        else math.fsum((-layer_generated, layer_condensate_hexane))
                    ),
                    common_datum_energy_to_cell_j=math.fsum(
                        (layer_packet_energy, wall_storage if index == 1 else 0.0)
                    ),
                    energy_datum_id=request.prior.tray.energy_datum_id,
                )
            )
        external_transfers = tuple(external_transfers_list)
        callback = _AcceptedPDPacketAdvance(
            entries=entries,
            residuals={item.packet_id: item.maximum_relative_residual for item in packet_ledgers},
            # EXPLORATION MODE F4 (Q6): the pinned host's residual pin, raised to an accepted
            # packet residual and logged (host_residual_pin_raised); the limit otherwise.
            maximum_scaled_residual=(
                request.relative_limit
                if PD_EXPLORATION_HOOK is None
                else PD_EXPLORATION_HOOK.host_residual_pin(request.relative_limit, packet_ledgers)
            ),
        )
        host_request = tray_host.KCellTrayStepRequest(
            attempt_id=request.attempt_id,
            prior=request.prior.tray,
            end_time_s=request.end_time_s,
            flow_segments=request.flow_segments,
            external_cell_transfers=external_transfers,
            relative_limit=request.relative_limit,
        )
        host_result = tray_host.evaluate_k_cell_tray_step(
            host_request,
            packet_callback=callback,
            k_planner=k_planner,
            auditor=codec,
        )
        if type(host_result) is tray_host.KCellTrayStepRejection:
            raise PDKDomeBypassStepError(f"K-cell host rejected the PD step: {host_result.reason}")
        if type(host_result) is not tray_host.ValidatedKCellTrayStep:
            raise PDKDomeBypassStepError("K-cell host returned a foreign transaction type")
        candidate = AcceptedPDKTrayState(
            tray=host_result.candidate,
            model_definition_digest=request.prior.model_definition_digest,
            wall_temperature_k=new_wall_temperature,
            plenum_outlet=outlet,
        )
        system_ledger = _system_ledger(
            request,
            host_result,
            wall=wall_ledger,
            plenum=plenum_ledger,
            deposition_bundle=deposition_bundle,
        )
        # EXPLORATION MODE F4 (Q5): admitted only when its failure is the sum of the packet
        # residuals already admitted in this step (its own linked row); a terminal otherwise.
        if not system_ledger.passed and not (
            PD_EXPLORATION_HOOK is not None
            and PD_EXPLORATION_HOOK.admit_system_ledger(system_ledger, packet_ledgers)
        ):
            raise PDKDomeBypassStepError("whole PD mass/common-energy ledger failed")
        return ValidatedPDKTrayStep(
            request=request,
            candidate=candidate,
            host_transaction=host_result,
            saturation_temperature_k=saturation_temperature,
            packet_ledgers=tuple(packet_ledgers),
            inter_k_face=inter_k_face,
            pressure_ledger=pressure_ledger,
            plenum_ledger=plenum_ledger,
            wall_ledger=wall_ledger,
            system_ledger=system_ledger,
            meal_mechanical_work_deposition_ledger=(
                None if deposition_bundle is None else deposition_bundle.ledger
            ),
            meal_mechanical_work_validation_context_digest=(
                None
                if meal_mechanical_work_validation_context is None
                else meal_mechanical_work_validation_context.definition_digest
            ),
            codec_registry_digest=codec.codec_registry_digest,
            auditor_identity_digest=codec.auditor_identity_digest,
            auditor_source_identity=codec.source_identity,
            candidate_digest=candidate.state_digest,
            boiling_authority_id=PD_BOILING_AUTHORITY_ID,
            certified_co_boiling_temperature_k=certified_co_boiling.temperature_k,
            certified_vapour_water_mole_fraction=(
                certified_co_boiling.vapour_water_mole_fraction
            ),
            boiling_temperature_excess_k=math.fsum(
                (saturation_temperature, -certified_co_boiling.temperature_k)
            ),
            condensation_ledger=condensation,
            moisture_pickup=_moisture_pickup_ledger(
                request, host_result, condensation=condensation
            ),
            _seal=_VALIDATED_PD_STEP_SEAL,
        )
    except (
        PDKDomeBypassIntegrationError,
        high_loading.HighLoadingPacketAdapterError,
        component_datum.ComponentEnergyDatumError,
        meal_work.QSCMealMechanicalWorkDepositionError,
    ) as error:
        return PDKTrayStepRejection(
            attempt_id=request.attempt_id,
            reason=str(error),
            rollback_state=request.prior,
        )
    except Exception as error:
        return PDKTrayStepRejection(
            attempt_id=request.attempt_id,
            reason=f"PD engineering callback failed: {type(error).__name__}: {error}",
            rollback_state=request.prior,
        )


def commit_pd_k_dome_bypass_tray_step(
    current: AcceptedPDKTrayState,
    transaction: ValidatedPDKTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner | None = None,
    meal_mechanical_work_validation_context: (
        meal_work.QSCMealMechanicalWorkValidationContext | None
    ) = None,
) -> AcceptedPDKTrayState:
    """Return the sealed in-memory candidate after reversible stale/tamper checks."""

    if type(current) is not AcceptedPDKTrayState:
        raise TypeError("PD commit requires an exact AcceptedPDKTrayState")
    if (
        type(transaction) is not ValidatedPDKTrayStep
        or transaction._seal is not _VALIDATED_PD_STEP_SEAL
    ):
        raise PDKDomeBypassIntegrationError("PD commit requires an evaluator-issued transaction")
    ValidatedPDKTrayStep.__post_init__(transaction)
    deposition_bundle = _validate_meal_mechanical_work_deposition(
        transaction.request,
        meal_mechanical_work_validation_context,
    )
    context_digest = (
        None
        if meal_mechanical_work_validation_context is None
        else meal_mechanical_work_validation_context.definition_digest
    )
    if (
        transaction.meal_mechanical_work_validation_context_digest != context_digest
        or transaction.meal_mechanical_work_deposition_ledger
        != (None if deposition_bundle is None else deposition_bundle.ledger)
    ):
        raise PDKDomeBypassIntegrationError("PD commit meal-work bundle or caller-held pins differ")
    if current is not transaction.request.prior:
        raise StalePDKDomeBypassStateError("PD transaction is stale or foreign")
    if current.state_digest != transaction.request.prior.state_digest:
        raise StalePDKDomeBypassStateError("PD accepted prior digest changed before commit")
    if transaction.candidate.state_digest != transaction.candidate_digest:
        raise PDKDomeBypassIntegrationError("PD candidate digest changed before commit")
    if (
        codec.codec_registry_digest != transaction.codec_registry_digest
        or codec.auditor_identity_digest != transaction.auditor_identity_digest
        or codec.source_identity != transaction.auditor_source_identity
    ):
        raise PDKDomeBypassIntegrationError("PD commit codec violates transaction pins")
    if k_planner is None:
        k_planner = tray_host.NoKRekeyPlanner()
    _validate_pd_k_planner(transaction.request, k_planner)
    reevaluated = evaluate_pd_k_dome_bypass_tray_step(
        transaction.request,
        codec=codec,
        k_planner=k_planner,
        meal_mechanical_work_validation_context=(meal_mechanical_work_validation_context),
    )
    if type(reevaluated) is not ValidatedPDKTrayStep or reevaluated != transaction:
        raise PDKDomeBypassIntegrationError(
            "PD transaction differs from a complete deterministic commit-time reevaluation"
        )
    committed_tray = tray_host.commit_k_cell_tray_step(
        current.tray,
        transaction.host_transaction,
        auditor=codec,
    )
    if committed_tray is not transaction.candidate.tray:
        raise PDKDomeBypassIntegrationError("host and PD candidates diverged at commit")
    if (
        not all(item.passed for item in transaction.packet_ledgers)
        or not transaction.inter_k_face.passed
        or not transaction.pressure_ledger.passed
        or not transaction.plenum_ledger.passed
        or not transaction.wall_ledger.passed
        or not transaction.system_ledger.passed
    ):
        raise PDKDomeBypassIntegrationError("PD acceptance ledgers no longer pass")
    return transaction.candidate


__all__ = (
    "AcceptedPDKTrayState",
    "DEFAULT_RELATIVE_LIMIT",
    "GAS_ENTHALPY_CONVENTION",
    "PD_INTER_K_ALLOCATION_LAW_ID",
    "PD_INTER_K_CONDUCTIVITY_AUTHORITY_ID",
    "PD_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K",
    "PD_K2_LAYER_COUNT",
    "PDBinaryGasStream",
    "PDBypassPressureLedger",
    "PDKDomeBypassConfigurationError",
    "PDKDomeBypassIntegrationError",
    "PDKDomeBypassStepError",
    "PDKTrayEngineeringModel",
    "PDKTrayStepRejection",
    "PDKTrayStepRequest",
    "PDInterKFaceLedger",
    "PDPacketHeatLedger",
    "PDPlenumLedger",
    "PDSharedWallLedger",
    "PDSharedWallParameters",
    "PDSystemLedger",
    "StalePDKDomeBypassStateError",
    "ValidatedPDKTrayStep",
    "commit_pd_k_dome_bypass_tray_step",
    "evaluate_pd_k_dome_bypass_tray_step",
)
