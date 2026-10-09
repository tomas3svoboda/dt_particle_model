"""C8 Tier 1: the sorbed-water arm's declared law for a through-bed layer.

WHAT THIS IS.  ``docs/GT_PS2_C8_TIER1_SORBED_WATER_ARM_DESIGN_2026-09-05.md``
section 1.2 charters one new, unpinned module holding three things and
nothing else:

* the **activity** a layer presents at its interface when its external water
  film is exactly exhausted and its hexane film is present - the accepted
  modified-Luikov law on its qualified domain (delegated to the pinned
  ``props.sorption.water_activity``, never re-derived here), the high-loading
  record's linear ENGINEERING BRIDGE above the fitted ceiling, and exactly
  ``1.0`` at and above the bound-to-free threshold;
* the **threshold** ``x_bf`` at which sorbed water gives way to a free film;
* the **kinetic case** - where a layer's condensate is booked, ``"sorbed"``
  (the working case) or ``"film"`` (today's routing verbatim, the swept
  limit).

Every numeric literal below is provenance class
``DECLARED_ENGINEERING_ASSUMPTION`` and is declared in
``DECLARED_ASSUMPTIONS`` with its criterion, its value and the outcome that
would REJECT it written first, per P-1.  This module mints no isotherm, fits
nothing and reads no measurement.

WHAT IT IS NOT.  It is not a physics decision: D-C8-2 (threshold law and
sweep) and D-C8-3 (kinetic case) were ruled by the owner on 2026-09-05 and
this module only carries them.  It does not touch the frozen Luikov
parameters, it does not extend their qualified domain downward (below
``W_ref`` it REFUSES, exactly as the pinned authority does), and it never
clamps: every out-of-domain input is a typed refusal.

THE CONTINUITY THE READER MUST SEE.  The bridge's left anchor is the
PUBLISHED pair ``(W_cap, 0.799)`` (high-loading record section 3), while the
Luikov correlation evaluated AT ``W_cap`` returns 0.7989999999899889 - the
cap is a ten-decimal published inversion of ``a_w = 0.799``, so the two
differ by 1.00e-11.  The bridge is anchored on the PUBLISHED value because
that is what reproduces the record's own worked number (0.85973 at
``X_w = 0.240`` with ``x_bf = 0.250``); the 1.00e-11 step at the join is
reported by :func:`bridge_join_discontinuity` rather than papered over, and
it is far below any quantity the interface law resolves.

D9-a (owner-ruled 2026-09-13, ``GT_PS2_OWNER_RULINGS_D9_2026-09-13.md``, built
to ``GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md``) adds a SECOND law beside the
C8 one and moves nothing in it.  Condensate is booked to the external film
ALWAYS; the film then drains into the hexane-free shell outside the receding
hexane front, and the arm exchanges with the gas through that shell only,
while no film is present.  The shell is the volume share ``1 - f**3`` of the
sphere whose front sits at dimensionless radius ``f``, which the kernel
already computes.  Its rows are ``d9a_routing``, ``d9a_drain``,
``d9a_exchangeable_water``, ``d9a_uniform_loading`` and the one D9-b constant
``d9b_no_shell_boundary``; its functions are :func:`open_shell_fraction`,
:func:`shell_capacity_kg`, :func:`uptake_capacity_kg`,
:func:`exchangeable_water_kg`, :func:`film_drain_kg`, :func:`drain_split_kg`,
:func:`drain_ledger_residual_kg` and
:func:`drain_bracket_end_for_destination`.  The C8 law id is UNCHANGED -
stamps exist under it - and the routing carries its own
:data:`D9A_WATER_ROUTING_LAW_ID` beside it.  The two existing arm methods
(:meth:`SorbedWaterArm.film_onset_ceiling_kg`,
:meth:`SorbedWaterArm.drawable_floor_kg`) and the activity function are
untouched, so a caller that does not use the new names is byte-identical.
This module still calls neither the kernel nor the cell: the D9-2 edits to
those two byte-pinned files are delivered as proposed diffs, not applied.

physically_qualifying: false.  plant_predictive: false.  No F-gate advances.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from types import MappingProxyType
from typing import ClassVar

from .props import sorption

#: Provenance class of every declared value in this module (D0).
DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"

#: Law id of the arm as a whole.  Bump the revision if a NUMBER or a rule
#: below changes meaning; adding a reported key is additive and keeps the id.
SORBED_WATER_ARM_LAW_ID = "c8-tier1-sorbed-water-arm-v1"

#: The frozen modified-Luikov parameters, read from the pinned authority.
#: This module never constructs its own isotherm constants.
LUIKOV_PARAMS = sorption.LuikovParams()

#: The published high-loading anchor: the last DIRECTLY MEASURED DT-outlet
#: soybean-meal activity, and the ten-decimal loading the frozen correlation
#: inverts it to.
BRIDGE_ANCHOR_ACTIVITY = 0.799
BRIDGE_ANCHOR_LOADING_KG_KG = LUIKOV_PARAMS.W_cap

#: D-C8-2 (owner-ruled 2026-09-05): the bound-to-free threshold.
X_BF_NOMINAL_KG_KG = 0.250
X_BF_SWEEP_BAND_KG_KG = (0.20, 0.32)

#: The closure's own first-tranche declared choice (activity-only hindrance).
WATER_CONDUCTANCE_FACTOR_NOMINAL = 1.0
WATER_CONDUCTANCE_FACTOR_SWEEP_BAND = (0.2, 1.0)

#: D-C8-3 (owner-ruled 2026-09-05): the kinetic case.
CONDENSATE_DESTINATIONS = ("sorbed", "film")
CONDENSATE_DESTINATION_NOMINAL = "sorbed"

#: The record's worked check, reproduced by :meth:`SorbedWaterArm.activity`.
KNOWN_ANSWER_LOADING_KG_KG = 0.240
KNOWN_ANSWER_ACTIVITY = 0.85973

# ---------------------------------------------------------------------------
# D9-a (owner-ruled 2026-09-13) and the one D9-b constant the cell reads.
# A SECOND law beside the C8 one; nothing above this line moves.
# ---------------------------------------------------------------------------

#: D9-a: the id of the water ROUTING law the amended ``condensate_destination``
#: token now selects.  It stands BESIDE :data:`SORBED_WATER_ARM_LAW_ID`, never
#: in place of it: the C8 id stamps evidence that already exists and keeps its
#: exact string.
D9A_WATER_ROUTING_LAW_ID = "d9a-shell-gated-imbibition-v1"

#: D9-a-2: the two declared ends of the imbibition bracket, as the token the
#: caller hands to :func:`film_drain_kg`.  ``"nominal"`` moves the whole opened
#: capacity in one step; ``"none"`` moves nothing.
DRAIN_BRACKET_END_NOMINAL = "nominal"
DRAIN_BRACKET_END_NONE = "none"
DRAIN_BRACKET_ENDS = (DRAIN_BRACKET_END_NOMINAL, DRAIN_BRACKET_END_NONE)

#: D9-a-1: the amended meaning of each ``condensate_destination`` token.
#: ``"sorbed"`` no longer means "book the condensate to the sorbed inventory";
#: it means the D9-a routing - condensate to the film ALWAYS, plus the drain at
#: its nominal.  ``"film"`` is the same routing with no drain, the lower end of
#: the same bracket, which is today's film lane.
CONDENSATE_DESTINATION_DRAIN_BRACKET_END = MappingProxyType(
    {
        CONDENSATE_DESTINATION_NOMINAL: DRAIN_BRACKET_END_NOMINAL,
        "film": DRAIN_BRACKET_END_NONE,
    }
)

#: D9-a-4: the hexane front's dimensionless radius.  ``f = 1`` is a closed
#: shell (the whole particle is behind the front: a film-bearing layer);
#: ``f = 0`` is a fully open particle (a dry-sorbate layer).  Outside this
#: closed interval the routing REFUSES typed; it never clamps.
FRONT_FRACTION_DOMAIN = (0.0, 1.0)

#: D9-b-2 AS AMENDED BY OWNER RULING D9-e (2026-09-13,
#: ``GT_PS2_OWNER_RULINGS_D9E_2026-09-13.md``): the front fraction AT which
#: there is no shell at all, so a film-free layer has nothing to exchange
#: water through and the dry-shell interface has nothing to carry the hexane
#: duty.  The value is 1.0 EXACTLY - the front at the surface - and the 0.99
#: boundary the D9-b build carried, with its (0.95, 0.99) bracket, is
#: WITHDRAWN: it was measured to be an artefact of the water-BLOCKED closure
#: (v), and closure (vi), with the water arm active, has a root at every
#: front position from the layer's own crossing to half radius.  Declared
#: here because the law module is unpinned and the cell and the kernel are
#: not; both read this constant rather than minting their own.
NO_SHELL_FRONT_FRACTION_CEILING = 1.0

#: D9-e, the declared HINDRANCE CONTINUATION that seeds the dry-shell
#: receding-front interface when a direct Newton from the caller's seed
#: refuses.  Both numbers are SEEDING-ONLY: every root the cell reports is the
#: solve at the record's OWN series factor, found and revalidated there by the
#: shipped damped Newton at the requested tolerance, so no root, ledger or
#: reported number depends on either value - what depends on them is only
#: whether the root is REACHED.  Declared here because the law module is
#: unpinned and the cell is not.
#:
#: START is the series factor the ladder begins at; STEPS is the number of
#: geometric intervals from that start down to the record's own factor, so the
#: ladder is ``Phi_k = START (Phi/START)**(k/STEPS)``, k = 0..STEPS, its last
#: rung the caller's record itself.  Fixed count, fixed ratio, no adaptive
#: tolerance and no halving: the same inputs give the same rungs and the same
#: root, bit for bit, on every call.
HINDRANCE_CONTINUATION_START = 0.5
#: The measured working band of the start (see the declared row below).
HINDRANCE_CONTINUATION_START_BRACKET = (0.05, 0.5)
HINDRANCE_CONTINUATION_STEPS = 12
#: The measured working band of the step count.
HINDRANCE_CONTINUATION_STEPS_BRACKET = (4, 24)

#: B1 STAGE 4 (owner ruling 2026-09-19; audit
#: ``docs/evidence/ruling_audit_2026-09-20/B1S4/AUDIT_B1S4.md``): the
#: EXTERNAL-WATER engineering-nil band, whose declared row is
#: ``external_water_nil_band_reclassification`` below.  Declared here because
#: the law module is unpinned and the kernel is not - the same reason the two
#: constants above are declared here.
#:
#: NO NEW NUMBER IS DECLARED.  This is the SAME value as the kernel's
#: ``sp1_k2_law2_tray_integration.NIL_FILM_RECLASSIFICATION_THRESHOLD_KG``,
#: which F-W1 also binds the water landing to as its deposit residual
#: ceiling.  It is written out rather than imported for one reason only: the
#: kernel imports THIS module, so this module cannot import the kernel
#: without a cycle.  The equality is therefore kept the way
#: ``layer_phase_exhaustion.water_depletion_stock`` keeps its mirror of the
#: kernel's own branch - BY TEST, asserted against the kernel's symbol, never
#: by a second decision.
EXTERNAL_WATER_NIL_BAND_KG = 1.0e-4
#: The band's measured bracket, both ends measurements: 0.001 and 0.1
#: monolayer of close-packed liquid water on the tower's smallest layer
#: (PD1 K1), from the shipped 3/R area law at R = 0.885 mm and 600 kg/m3.
EXTERNAL_WATER_NIL_BAND_BRACKET_KG = (2.8886902986632124e-06, 2.8886902986632124e-04)
#: The ACTIVE-SET step the crossing carries, measured on MN2 K1 at
#: W = 0.22009 kg/kg dry: the surface water activity falls from exactly 1.0
#: (the two-liquid law) to the arm's Luikov value 0.78185190540426.  Reported,
#: stamped and declared - never called immaterial.
EXTERNAL_WATER_NIL_BAND_ACTIVITY_STEP = -0.21814809459574003
#: Mirror of the kernel's ``NIL_WATER_FILM_RECLASSIFICATION_LAW_ID``, kept
#: equal by the same test that keeps the band equal.
EXTERNAL_WATER_NIL_BAND_LAW_ID = "b1-stage4-external-water-nil-band-v1"

DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "x_bf": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the loading at which sorbed water in DT meal gives way to a free "
                    "surface film, i.e. the loading above which the surface water "
                    "activity is one and the layer must run the two-liquid law"
                ),
                "value_kg_kg_dry": X_BF_NOMINAL_KG_KG,
                "sweep_kg_kg_dry": list(X_BF_SWEEP_BAND_KG_KG),
                "sources": (
                    "Karnofsky 1985 (20 % wet basis at the toasting trays); "
                    "Aviara 2004 whole soybean 0.20-0.22; Cassini 2006 texturized "
                    "soy 0.28; Zhong and Sun 2000 isolated 11S protein 0.30-0.32 by "
                    "DSC; collected in "
                    "docs/GT_PS2_C8_SORBED_WATER_INTERFACE_SCOPE_AND_LITERATURE_"
                    "2026-09-05.md section 5"
                ),
                "rejecting_outcome": (
                    "the first countercurrent tray's meal overshoots the published "
                    "17-22 % wet-basis outlet band with the threshold held anywhere "
                    "in 0.20-0.32, or it never reaches the band at the top of that "
                    "sweep: then the threshold or the bridge is wrong, not the "
                    "transport group"
                ),
            }
        ),
        "bridge": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "a monotone continuation of the accepted isotherm from its last "
                    "directly measured point to saturation, since no soybean-meal "
                    "measurement exists at 0.20-0.35 kg/kg dry and 80-110 C"
                ),
                "value": (
                    "linear in loading from (W_cap, 0.799) to (x_bf, 1.0); exactly "
                    "1.0 at and above x_bf; the pinned Luikov law below W_cap; a "
                    "typed refusal below W_ref"
                ),
                "anchor_loading_kg_kg_dry": BRIDGE_ANCHOR_LOADING_KG_KG,
                "anchor_activity": BRIDGE_ANCHOR_ACTIVITY,
                "sources": (
                    "docs/GT_PS2_HIGH_LOADING_WATER_SORPTION_2026-09-04.md "
                    "sections 3 and 6 (its own worked value 0.85973 at 0.240 with "
                    "x_bf 0.250 is reproduced here as a known answer)"
                ),
                "rejecting_outcome": (
                    "a measured soybean-meal isotherm above a_w = 0.799 at DT "
                    "temperatures that is not monotone-linear within the band the "
                    "sweep spans; the bridge is then replaced by that measurement "
                    "and this module's middle segment is deleted"
                ),
            }
        ),
        "water_conductance_factor": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the series hindrance the sorbed water phase adds to the "
                    "water-side film conductance, over and above the activity "
                    "depression already carried by the isofugacity row"
                ),
                "value": WATER_CONDUCTANCE_FACTOR_NOMINAL,
                "sweep": list(WATER_CONDUCTANCE_FACTOR_SWEEP_BAND),
                "sources": (
                    "cell_engineering_feasibility.CoreSuppliedWaterClosure docstring "
                    "(first-tranche declared choice: activity-only hindrance); no "
                    "measurement of an internal water resistance for DT meal exists"
                ),
                "rejecting_outcome": (
                    "the condensation duty on the first countercurrent tray misses "
                    "the published 52-76 kg/t band in the SAME direction across the "
                    "whole [0.2, 1.0] sweep: the missing resistance is then not on "
                    "the water side of the film at all"
                ),
            }
        ),
        "condensate_destination": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "where water condensed on a film-free layer is booked within one "
                    "macro step: into the packets' sorbed inventory (instant "
                    "redistribution) or onto the surface as a film (no "
                    "redistribution within the step)"
                ),
                "value": CONDENSATE_DESTINATION_NOMINAL,
                "sweep": list(CONDENSATE_DESTINATIONS),
                "sources": (
                    "Karnofsky 1985 pp. 693-695 observes the flake surface wets "
                    "first and redistributes inward slowly; no source gives the "
                    "redistribution RATE for soybean flakes under steam, so the two "
                    "limits are swept and neither is claimed"
                ),
                "rejecting_outcome": (
                    "the two limits give condensation duties that straddle the "
                    "published band from opposite sides: the kinetics are then "
                    "load-bearing and a rate law must be measured, not declared"
                ),
                # ADDED 2026-09-13, additive: the row above is kept verbatim
                # because stamps exist under it; D9-a-1 amends what its two
                # tokens MEAN, and the amended semantics are the d9a_routing
                # row below.  Nothing in this row's value or sweep moved.
                "amended_by": (
                    "d9a_routing (owner ruling 2026-09-13): the condensate is "
                    "booked to the external film under BOTH tokens; 'sorbed' now "
                    "selects the D9-a drain at its nominal and 'film' selects no "
                    "drain, so the switch is the imbibition bracket rather than "
                    "the destination"
                ),
            }
        ),
        "allocation_basis": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "how one layer's sorbed-water change is divided among the "
                    "packets that make up the layer"
                ),
                "value": "dry matter (every packet's loading moves by the same increment)",
                "sweep": (
                    "capacity-weighted, proportional to (x_bf - X_w), is the named "
                    "alternative (design section 5); not built in Tier 1"
                ),
                "sources": "design record section 1.4",
                "rejecting_outcome": (
                    "the within-layer spread of loadings grows large enough that the "
                    "two bases give materially different tray outlet moisture; the "
                    "basis is then a physics question rather than bookkeeping"
                ),
            }
        ),
        "exhaustion_floor": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the lowest sorbed loading from which the arm may still draw "
                    "water by evaporation"
                ),
                "value_kg_kg_dry": LUIKOV_PARAMS.W_ref,
                "sweep": "none: it is the accepted law's own qualified lower bound",
                "sources": (
                    "props/sorption.py water_activity: 'no low-moisture continuation "
                    "is authorized' (packet A5)"
                ),
                "rejecting_outcome": (
                    "a march needs to cross it: that is the Tier 2 DRY_SURFACE "
                    "branch, and Tier 1 refuses typed rather than continuing"
                ),
            }
        ),
        # -------------------------------------------------------------------
        # D9-a and D9-b (owner-ruled 2026-09-13).  Rejecting outcome first.
        # -------------------------------------------------------------------
        "d9a_routing": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "where water condensing on a layer is booked: liquid water present "
                    "at the surface sets the surface water activity to one, and a shell "
                    "still closed by free hexane cannot take that condensate in, so the "
                    "condensate is liquid water AT the surface and the film is the "
                    "physical state there"
                ),
                "value": (
                    "condensate to the external film ALWAYS; the two "
                    "condensate_destination tokens keep their spellings with amended "
                    "meanings - 'sorbed' is the D9-a routing (film plus the drain of "
                    "row d9a_drain at its nominal) and 'film' is the film with no "
                    "drain, the lower end of the same bracket"
                ),
                "bracket": list(DRAIN_BRACKET_ENDS),
                "sources": (
                    "docs/GT_PS2_FT2A1_CLOSURE_MEASUREMENT_RECORD_2026-09-12.md "
                    "section 4 (at the crossing the surface sits 45 K below the film "
                    "lane and far below the steam dew point under every water closure, "
                    "and the hexane-free shell is 0.77 nm); "
                    "docs/GT_PS2_OWNER_RULINGS_D9_2026-09-13.md section 1; "
                    "docs/GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md row D9-a-1"
                ),
                "rejecting_outcome": (
                    "a per-layer water ledger that does not close on any split: the "
                    "routing is then wrong as bookkeeping before it is wrong as physics"
                ),
            }
        ),
        "d9a_drain": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "the rate at which the surface film imbibes into the hexane-free "
                    "shell: the shell's water capacity opens 5 to 21 times faster than "
                    "condensate arrives, so filling the OPENED capacity within the step "
                    "is the upper limit of the bracket and film-first with no drain is "
                    "its lower limit; no source gives the rate between them"
                ),
                "value": (
                    "drain = min(film_kg, capacity_kg) with capacity_kg = "
                    "max(0, x_bf * dry - retained) * (1 - f**3); booked as a DEPOSIT "
                    "inside the same advance call, on every packet of the layer in "
                    "proportion to its dry matter, after the Law-2 advance and before "
                    "the film-onset check"
                ),
                "bracket": (
                    "nominal: the whole opened capacity in one step; lower limit: "
                    "drain = 0, which is today's film lane"
                ),
                "caloric_treatment": (
                    "the moved water carries its enthalpy at the layer temperature and "
                    "NO sorption heat: the shipped net excess binding enthalpy of "
                    "retained water is exactly 0.0 (props/sorption.py:525-532, "
                    "water_q_net, PHY-021.water); cited, not re-derived here"
                ),
                "sources": (
                    "docs/GT_PS2_PORE_BLOCKING_MEASUREMENT_RECORD_2026-09-12.md "
                    "section 2 (capacity opens at 4.1 to 1.0 kg/s against a mean "
                    "condensation 20.8, 5.3 and 4.9 times smaller; 1 of 561 interval "
                    "pairs had condensate exceed the capacity opened; minimum "
                    "film-to-capacity ratio 1.275 at 28 s); "
                    "docs/GT_PS2_C8_SORBED_WATER_INTERFACE_SCOPE_AND_LITERATURE_"
                    "2026-09-05.md item 6 for the two limiting cases; "
                    "docs/GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md row D9-a-2"
                ),
                "rejecting_outcome": (
                    "the two ends of the bracket giving outlet states on OPPOSITE "
                    "sides of a published band: the rate is then load-bearing and "
                    "unmeasured, and this row cannot be declared at its nominal"
                ),
            }
        ),
        "d9a_exchangeable_water": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "while free hexane fills the pores the matrix water there cannot "
                    "exchange with the gas; only the hexane-free shell outside the "
                    "receding front can, and only while no film is present - behind "
                    "the front the water is inert"
                ),
                "value": (
                    "drawable stock = (retained - W_ref * dry) * (1 - f**3); uptake "
                    "capacity through the shell = (x_bf * dry - retained) * (1 - f**3); "
                    "the arm's exhaustion floor becomes retained - drawable.  The "
                    "activity stays the isotherm's at the layer loading and the water "
                    "conductance factor stays 1.0: the gating is by the shell, not by "
                    "a factor"
                ),
                "bracket": (
                    "none on the form; its one free number is the front fraction f, "
                    "which the declared mass Biot number sets "
                    "(sp1_k2_law2_tray_integration.py:112, D4 held)"
                ),
                "sources": (
                    "docs/GT_PS2_OWNER_RULINGS_D9_2026-09-13.md section 1 items (i) "
                    "and (ii); "
                    "docs/GT_PS2_FT2A1_CLOSURE_MEASUREMENT_RECORD_2026-09-12.md "
                    "section 2 (at the declared Biot number front-gated water IS "
                    "blocked water at the crossing, to five significant figures); "
                    "docs/GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md row D9-a-3"
                ),
                "rejecting_outcome": (
                    "as d9a_routing: a per-layer water ledger that does not close on any split"
                ),
            }
        ),
        "d9a_uniform_loading": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "a layer-level model with no radial water coordinate: the retained "
                    "water loading is ONE number per layer, so the shell's share of it "
                    "can only be the geometric volume share"
                ),
                "value": (
                    "sphere, front at dimensionless radius f, shell share 1 - f**3; "
                    "the layer's drain and draw are divided among its packets on the "
                    "dry-matter basis of the allocation_basis row, so every packet's "
                    "loading moves by the same increment"
                ),
                "bracket": "none at this modelling level",
                "sources": (
                    "docs/GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md row D9-a-4; the "
                    "same sphere geometry the pore-blocking record measured the shell "
                    "with (its section 2)"
                ),
                "rejecting_outcome": (
                    "none at this level: recorded as the abstraction it is.  Pores at "
                    "one radius are not all in one state, which is why a real sealing "
                    "would be partial; a radial water coordinate would replace this "
                    "row rather than reject it"
                ),
            }
        ),
        "d9e_hindrance_continuation_start": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "the series factor at which the dry-shell interface's declared "
                    "seeding continuation begins: a hindrance at which the caller's "
                    "own seed still reaches a root, from which the ladder walks down "
                    "to the record's own factor.  SEEDING ONLY: every reported root "
                    "is the solve at the record's own factor, so this value decides "
                    "whether a root is reached and never what it is"
                ),
                "value": HINDRANCE_CONTINUATION_START,
                "bracket": list(HINDRANCE_CONTINUATION_START_BRACKET),
                "sources": (
                    "measured at D9-e level 0 on nine states - the honest lane's "
                    "step-40 cell, the captured d8dcont cell at four imposed fronts "
                    "and the C8 Tier-2a cell at four - over the starts 0.9, 0.5, 0.3, "
                    "0.1, 0.05 and 0.02: 0.5 down to 0.05 reach the root on all nine, "
                    "0.9 refuses at the first rung on the C8 cell (a sub-nanometre "
                    "shell) and 0.02 refuses at the first rung on the lane cell.  The "
                    "value is the widest working start, which leaves the ladder the "
                    "most room; docs/GT_PS2_D9E_LEVEL0_BUILD_RECORD_2026-09-13.md"
                ),
                "rejecting_outcome": (
                    "a state a lane reaches at which the continuation refuses at its "
                    "FIRST rung: the start is then outside the band this row was "
                    "measured on and must be re-measured, not nudged.  The refusal "
                    "names the rung, so the case is visible rather than silent"
                ),
            }
        ),
        "d9e_hindrance_continuation_steps": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "the number of geometric intervals the continuation takes from "
                    "the declared start to the record's own series factor, the ladder "
                    "re-seeded from each rung's root.  Fixed: no adaptive step, no "
                    "halving, no tolerance that moves - which is what makes the "
                    "reported root bit-for-bit reproducible on repeated calls"
                ),
                "value": HINDRANCE_CONTINUATION_STEPS,
                "bracket": list(HINDRANCE_CONTINUATION_STEPS_BRACKET),
                "sources": (
                    "measured at D9-e level 0 on the same nine states at the declared "
                    "start: 4, 6, 8, 10, 12, 16 and 24 intervals all reach the root on "
                    "all nine, to the same value; 12 is the middle of the measured "
                    "band and costs 13 solves; "
                    "docs/GT_PS2_D9E_LEVEL0_BUILD_RECORD_2026-09-13.md"
                ),
                "rejecting_outcome": (
                    "a state a lane reaches at which the continuation refuses at a "
                    "rung ABOVE the first: the ladder is then too coarse there and the "
                    "count must be re-measured.  The refusal names the rung and its "
                    "series factor, so which end failed is never a guess"
                ),
            }
        ),
        "d9b_no_shell_boundary": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": D9A_WATER_ROUTING_LAW_ID,
                "criterion": (
                    "the front fraction AT which the particle has no hexane-free "
                    "shell at all, so a film-free layer has nothing to exchange "
                    "water through and the dry-shell receding-front interface has "
                    "nothing to carry the hexane duty: the front at the surface, "
                    "f = 1 exactly, which is not a state D9-a's routing produces"
                ),
                "value": NO_SHELL_FRONT_FRACTION_CEILING,
                "bracket": (
                    "WITHDRAWN by owner ruling D9-e (2026-09-13): the 0.99 value "
                    "and its (0.95, 0.99) bracket were artefacts of the "
                    "water-blocked closure (v), and a boundary at exactly one "
                    "carries no bracket because it is the geometry of a front at "
                    "the surface, not a fitted number"
                ),
                "sources": (
                    "docs/GT_PS2_OWNER_RULINGS_D9E_2026-09-13.md (the amendment) on "
                    "docs/GT_PS2_RULING_PACKET_D9E_2026-09-13.md section 3; the "
                    "measurement behind it is closure (vi), "
                    "docs/evidence/ft2a1_closure_2026-09-12/"
                    "FT2A1_CLOSURE_MEASUREMENT_SECTION_8_ADDENDUM.md Table N, which "
                    "converges at every front from the layer's own 0.99999913 to "
                    "half radius with the water condensing into the shell; the "
                    "refuted predecessor is "
                    "docs/GT_PS2_FT2A1_CLOSURE_MEASUREMENT_RECORD_2026-09-12.md "
                    "section 3b with "
                    "docs/GT_PS2_D9_BUILD_DESIGN_BRIEF_2026-09-13.md row D9-b-2"
                ),
                "rejecting_outcome": (
                    "the refusal firing on any D9-a lane (kill criterion K8 as "
                    "D9-e re-defines it): the routing would then have put a "
                    "film-free layer at a closed shell, which its own gate "
                    "forbids, and the drain must be re-measured"
                ),
            }
        ),
        "external_water_nil_band_reclassification": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": EXTERNAL_WATER_NIL_BAND_LAW_ID,
                "criterion": (
                    "the external water mass at or below which a layer no longer "
                    "carries a FILM.  The criterion is geometric, not thermal: at "
                    "most 0.035 monolayers of close-packed liquid water on any "
                    "layer of this tower, read off the shipped 3/R area law "
                    "(src/dtdc_simulator/core2/area.py) at R = 0.885 mm and "
                    "600 kg/m3 - 0.01039 monolayers on MN2 K1, 0.03462 on the "
                    "wettest layer PD1 K1 - which is two orders below a "
                    "continuous film.  At or below it the CELL VIEW books the "
                    "remnant as sorbed water so the C8 Tier 1 arm binds; the "
                    "packets' own payloads are untouched and the same step's "
                    "D9-a drain moves the real mass, so the step is "
                    "mass-conserving.  It is NOT a deadband on the frozen "
                    "PHY-007/023 exact-inventory predicate: the reclassification "
                    "is of the INVENTORY VIEW and every downstream predicate "
                    "still receives, and still compares, an exact zero"
                ),
                "value_kg": EXTERNAL_WATER_NIL_BAND_KG,
                "no_new_number": (
                    "the value IS the kernel's shipped "
                    "NIL_FILM_RECLASSIFICATION_THRESHOLD_KG, read from "
                    "sp1_k2_law2_tray_integration and never restated here, which "
                    "F-W1 already binds the water landing to as its deposit "
                    "residual ceiling, so the two constructions agree by "
                    "construction"
                ),
                "bracket_kg": list(EXTERNAL_WATER_NIL_BAND_BRACKET_KG),
                "bracket_ends": (
                    "both ends are MEASUREMENTS, not choices: 0.001 and 0.1 "
                    "monolayer on the tower's smallest layer (PD1 K1), i.e. "
                    "2.8886902986632124e-06 and 2.8886902986632124e-04 kg"
                ),
                "gate": (
                    "fires only where the layer's model declares a sorbed-water "
                    "arm (SP1KLayerEngineeringModel.declared_sorbed_water_arm is "
                    "not None), which is the containment: the arm is None on "
                    "every certified run and every shipped model, so those runs "
                    "are bit-identical by construction"
                ),
                "active_set_step": EXTERNAL_WATER_NIL_BAND_ACTIVITY_STEP,
                "active_set_step_note": (
                    "the band is picokelvin in MASS (evaporating 1e-4 kg moves "
                    "the measured MN2 layer by 1.36e-05 K) but it is NOT "
                    "immaterial in the ACTIVE SET: at the crossing the surface "
                    "water activity steps from exactly 1.0, the two-liquid law, "
                    "to the arm's Luikov value at the layer's retained loading - "
                    "0.78185190540426 at W = 0.22009, a step of -0.2181, 22 % of "
                    "the water driving force - and the area law flips branch at "
                    "the same instant because external_liquid_present loses its "
                    "last term.  The step is declared and stamped, never called "
                    "nil"
                ),
                "sources": (
                    "owner ruling 2026-09-19 ('I rule towards B1 stage 4 - "
                    "proceed with audit and gated build'); the audit "
                    "docs/evidence/ruling_audit_2026-09-20/B1S4/AUDIT_B1S4.md "
                    "sections 3 and 9 with its scripts c2_monolayer_criterion.py, "
                    "c_nil_band_materiality.py and f_band_reach_census.py; the "
                    "design note "
                    "docs/GT_PS2_B1_STAGE4_WATER_SYMMETRY_DESIGN_NOTE_2026-09-01"
                    ".md, WEAKENED by that audit and kept only for its framing; "
                    "C8 Tier 1 "
                    "docs/GT_PS2_C8_SORBED_WATER_INTERFACE_SCOPE_AND_LITERATURE_"
                    "2026-09-05.md; D9-a "
                    "docs/GT_PS2_OWNER_RULINGS_D9_2026-09-13.md; the build record "
                    "docs/GT_PS2_B1_STAGE4_LEVEL0_BUILD_RECORD_2026-09-20.md"
                ),
                "rejecting_outcome": (
                    "FOUR outcomes reject this row, each written before the "
                    "build. (1) IDENTITY: a certified run - one whose model "
                    "declares no sorbed-water arm - changes any leaf against its "
                    "pre-build twin; the gate would then not be the containment "
                    "it is claimed to be. (2) CONSERVATION: the water ledger of a "
                    "reclassifying step does not close to the kernel's own "
                    "relative limit; the view-only reclassification would then be "
                    "leaking mass, which is the failure mode the F-W1 deposit "
                    "operator was rejected for. (3) CONTINUITY: a reclassified "
                    "packet's water-activity step is LARGER than the isotherm's "
                    "own step between free and sorbed water at that content "
                    "(the stamped -0.2181 at W = 0.22009); the band would then be "
                    "creating a discontinuity rather than crossing the one the "
                    "two-branch surface model already carries. (4) MATERIALITY: "
                    "any of the 60 in-band layer-intervals the audit counted on "
                    "the op200sb fixture corpus (of 588, MN2 K1 26 and MN2 K2 34) "
                    "moves a water ledger by MORE than the mass reclassified on "
                    "that interval; the band would then be doing work beyond "
                    "re-labelling a remnant"
                ),
            }
        ),
    }
)


class SorbedWaterLawError(ValueError):
    """Typed refusal from the C8 Tier 1 sorbed-water arm."""


class SorbedWaterFloorRefusal(SorbedWaterLawError):
    """The loading is below the accepted law's qualified lower bound."""


class SorbedWaterFrontRefusal(SorbedWaterLawError):
    """The hexane front fraction is outside the declared ``[0, 1]`` domain."""


class SorbedWaterStockRefusal(SorbedWaterLawError):
    """A mass handed to the D9-a routing is not a finite, non-negative binary64."""


class SorbedWaterDrainBracketRefusal(SorbedWaterLawError):
    """The drain bracket end, or the token selecting it, is not declared."""


def _require_mass_kg(name: str, value: float) -> float:
    """A stock or a transfer in kilograms: finite, exact binary64, not negative."""

    if type(value) is not float or not math.isfinite(value):
        raise SorbedWaterStockRefusal(f"{name} must be a finite exact binary64")
    if value < 0.0:
        raise SorbedWaterStockRefusal(
            f"{name} {value!r} kg is negative; the D9-a routing has no negative "
            "stock and does not clamp one to zero"
        )
    return value


def _require_x_bf(value: float) -> float:
    """The threshold as the caller supplies it, on the ruled sweep."""

    if type(value) is not float or not math.isfinite(value):
        raise SorbedWaterLawError("x_bf must be a finite exact binary64")
    low, high = X_BF_SWEEP_BAND_KG_KG
    if not low <= value <= high:
        raise SorbedWaterLawError(
            f"x_bf {value!r} kg/kg dry is outside the ruled sweep [{low}, {high}] "
            "(D-C8-2); widen the ruling, do not clamp"
        )
    return value


def open_shell_fraction(*, front_fraction: float) -> float:
    """D9-a-4: the hexane-free volume share of the sphere, ``1 - f**3``.

    ``f = 1`` is a closed shell and returns exactly ``0.0``; ``f = 0`` is the
    whole particle and returns exactly ``1.0``.  Written as three
    multiplications and one subtraction, so both ends are exact and no
    division enters.  Outside ``[0, 1]`` this REFUSES typed.
    """

    if type(front_fraction) is not float or not math.isfinite(front_fraction):
        raise SorbedWaterFrontRefusal("the hexane front fraction must be a finite exact binary64")
    low, high = FRONT_FRACTION_DOMAIN
    if not low <= front_fraction <= high:
        raise SorbedWaterFrontRefusal(
            f"hexane front fraction {front_fraction!r} is outside the declared "
            f"domain [{low}, {high}]; there is no shell share to report - refuse, "
            "do not clamp"
        )
    return 1.0 - front_fraction * front_fraction * front_fraction


def shell_capacity_kg(
    *,
    dry_matter_kg: float,
    retained_water_kg: float,
    front_fraction: float,
    x_bf: float,
) -> float:
    """D9-a-2: the water the OPEN shell can still take, ``max(0, x_bf*dry - retained) * (1 - f**3)``.

    Exactly ``0.0`` on a closed shell (``f = 1``) and on a layer already at or
    above the threshold: those are the definition of the capacity, not a clamp
    of an input - an out-of-domain input is refused above.
    """

    _require_mass_kg("dry_matter_kg", dry_matter_kg)
    _require_mass_kg("retained_water_kg", retained_water_kg)
    _require_x_bf(x_bf)
    shell = open_shell_fraction(front_fraction=front_fraction)
    unfilled_kg = math.fsum((x_bf * dry_matter_kg, -retained_water_kg))
    if unfilled_kg <= 0.0:
        return 0.0
    return unfilled_kg * shell


def uptake_capacity_kg(
    *,
    dry_matter_kg: float,
    retained_water_kg: float,
    front_fraction: float,
    x_bf: float,
) -> float:
    """D9-a-3: the arm's uptake ceiling through the shell - the same number as
    :func:`shell_capacity_kg`, under the name the arm's ceiling reads it by."""

    return shell_capacity_kg(
        dry_matter_kg=dry_matter_kg,
        retained_water_kg=retained_water_kg,
        front_fraction=front_fraction,
        x_bf=x_bf,
    )


def exchangeable_water_kg(
    *,
    dry_matter_kg: float,
    retained_water_kg: float,
    front_fraction: float,
) -> float:
    """D9-a-3: the drawable stock, ``max(0, retained - W_ref*dry) * (1 - f**3)``.

    Exactly ``0.0`` on a closed shell and at or below the accepted law's
    qualified floor.  Behind the front the water is inert, which is what the
    shell factor says.
    """

    _require_mass_kg("dry_matter_kg", dry_matter_kg)
    _require_mass_kg("retained_water_kg", retained_water_kg)
    shell = open_shell_fraction(front_fraction=front_fraction)
    drawable_kg = math.fsum((retained_water_kg, -(LUIKOV_PARAMS.W_ref * dry_matter_kg)))
    if drawable_kg <= 0.0:
        return 0.0
    return drawable_kg * shell


def drain_bracket_end_for_destination(*, condensate_destination: str) -> str:
    """D9-a-1: the amended semantics of the two ``condensate_destination`` tokens.

    ``"sorbed"`` selects the nominal drain, ``"film"`` selects no drain.  Both
    tokens route the condensate itself to the external film: that is the part
    of the routing that is no longer a choice.
    """

    if condensate_destination not in CONDENSATE_DESTINATION_DRAIN_BRACKET_END:
        raise SorbedWaterDrainBracketRefusal(
            f"condensate_destination must be one of {CONDENSATE_DESTINATIONS} "
            "(D-C8-3, amended by D9-a-1)"
        )
    return CONDENSATE_DESTINATION_DRAIN_BRACKET_END[condensate_destination]


def film_drain_kg(*, film_kg: float, capacity_kg: float, bracket_end: str) -> float:
    """D9-a-2: the film-to-shell transfer at one end of the declared bracket.

    ``"nominal"`` moves the whole opened capacity the film can supply,
    ``min(film, capacity)``; ``"none"`` moves exactly ``0.0``.  Nothing between
    them is offered, because no source gives the rate.
    """

    _require_mass_kg("film_kg", film_kg)
    _require_mass_kg("capacity_kg", capacity_kg)
    if bracket_end not in DRAIN_BRACKET_ENDS:
        raise SorbedWaterDrainBracketRefusal(
            f"bracket_end must be one of {DRAIN_BRACKET_ENDS} (D9-a-2)"
        )
    if bracket_end == DRAIN_BRACKET_END_NONE:
        return 0.0
    return min(film_kg, capacity_kg)


def drain_split_kg(
    *,
    film_kg: float,
    retained_water_kg: float,
    drain_kg: float,
) -> tuple[float, float]:
    """The two stocks after the drain, as the kernel would book them.

    ``(film_after, retained_after)``.  The drain may not exceed the film - the
    caller gets it from :func:`film_drain_kg`, which cannot return more - and a
    larger one REFUSES rather than being trimmed.  Both ends of the bracket and
    the full drain (``drain == film``) are EXACTLY conservative: ``film_after``
    is then ``0.0`` exactly and ``retained_after`` is the correctly rounded
    pair sum.  In between, the split closes to binary64 resolution and
    :func:`drain_ledger_residual_kg` reports the residual rather than hiding
    it - measured, not assumed.
    """

    _require_mass_kg("film_kg", film_kg)
    _require_mass_kg("retained_water_kg", retained_water_kg)
    _require_mass_kg("drain_kg", drain_kg)
    if drain_kg > film_kg:
        raise SorbedWaterStockRefusal(
            f"drain {drain_kg!r} kg exceeds the film stock {film_kg!r} kg; the "
            "film cannot supply it - refuse, do not trim"
        )
    return math.fsum((film_kg, -drain_kg)), math.fsum((retained_water_kg, drain_kg))


def drain_ledger_residual_kg(
    *,
    film_kg: float,
    retained_water_kg: float,
    drain_kg: float,
) -> float:
    """The water the split neither kept nor moved, in kilograms.

    Zero at both bracket ends and at the full drain; at most one binary64
    resolution unit of the pair total otherwise.  Reported, never corrected -
    the same treatment :func:`bridge_join_discontinuity` gets.
    """

    film_after, retained_after = drain_split_kg(
        film_kg=film_kg,
        retained_water_kg=retained_water_kg,
        drain_kg=drain_kg,
    )
    return math.fsum((film_after, retained_after)) - math.fsum((film_kg, retained_water_kg))


def bridge_join_discontinuity() -> float:
    """The published anchor minus the Luikov correlation AT the cap.

    Reported, never corrected: see the module docstring.  Measured
    1.00e-11 at the frozen parameters.
    """

    return BRIDGE_ANCHOR_ACTIVITY - sorption.water_activity(
        BRIDGE_ANCHOR_LOADING_KG_KG, LUIKOV_PARAMS
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SorbedWaterArm:
    """One layer's declared sorbed-water arm: threshold, hindrance, routing.

    Carried on ``SP1KLayerEngineeringModel.declared_sorbed_water_arm``.
    ``None`` there - the default and every shipped model - keeps the sealed
    two-film law byte-identical, and a film-free layer keeps refusing with
    the frozen "Law 2 is executable only for TWO_EXTERNAL_LIQUIDS".
    """

    x_bf: float = X_BF_NOMINAL_KG_KG
    water_conductance_factor: float = WATER_CONDUCTANCE_FACTOR_NOMINAL
    condensate_destination: str = CONDENSATE_DESTINATION_NOMINAL
    law_id: str = SORBED_WATER_ARM_LAW_ID

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for name, value in (
            ("x_bf", self.x_bf),
            ("water_conductance_factor", self.water_conductance_factor),
        ):
            if type(value) is not float or not math.isfinite(value):
                raise SorbedWaterLawError(f"{name} must be a finite exact binary64")
        low, high = X_BF_SWEEP_BAND_KG_KG
        if not low <= self.x_bf <= high:
            raise SorbedWaterLawError(
                f"x_bf {self.x_bf!r} kg/kg dry is outside the ruled sweep "
                f"[{low}, {high}] (D-C8-2); widen the ruling, do not clamp"
            )
        if self.x_bf <= LUIKOV_PARAMS.W_ref:
            raise SorbedWaterLawError(
                f"x_bf {self.x_bf!r} kg/kg dry is at or below the accepted law's "
                f"qualified floor {LUIKOV_PARAMS.W_ref}: the arm would have no domain"
            )
        low, high = WATER_CONDUCTANCE_FACTOR_SWEEP_BAND
        if not low <= self.water_conductance_factor <= high:
            raise SorbedWaterLawError(
                f"water_conductance_factor {self.water_conductance_factor!r} is "
                f"outside the declared sweep [{low}, {high}]"
            )
        if self.condensate_destination not in CONDENSATE_DESTINATIONS:
            raise SorbedWaterLawError(
                f"condensate_destination must be one of {CONDENSATE_DESTINATIONS} (D-C8-3)"
            )
        if self.law_id != SORBED_WATER_ARM_LAW_ID:
            raise SorbedWaterLawError("the sorbed-water arm must carry its exact law id")

    @property
    def bridged(self) -> bool:
        """True when a bridge segment of positive width exists.

        ``x_bf <= W_cap`` is the high-loading record's "lower stress case":
        truncate the Luikov branch at the threshold and set ``a_w = 1`` above
        it; do not construct a negative-width bridge (its section 3).
        """

        return self.x_bf > BRIDGE_ANCHOR_LOADING_KG_KG

    def activity(self, loading_kg_kg_dry: float) -> float:
        """The interface water activity at a sorbed loading, ``a_w`` in (0, 1].

        Three segments and one refusal, in the design's own order:
        the pinned Luikov law on its qualified domain, the declared linear
        bridge above the evidence cap, exactly ``1.0`` at and above the
        threshold, and a typed refusal below ``W_ref``.
        """

        loading = loading_kg_kg_dry
        if type(loading) is not float or not math.isfinite(loading):
            raise SorbedWaterLawError("the retained-water loading must be a finite binary64")
        if loading < LUIKOV_PARAMS.W_ref:
            raise SorbedWaterFloorRefusal(
                f"retained-water loading {loading!r} kg/kg dry is below the "
                f"qualified Luikov lower bound {LUIKOV_PARAMS.W_ref}; no "
                "low-moisture continuation is authorized (Tier 2 DRY_SURFACE)"
            )
        if loading >= self.x_bf:
            return 1.0
        if loading <= BRIDGE_ANCHOR_LOADING_KG_KG:
            return sorption.water_activity(loading, LUIKOV_PARAMS)
        # x_bf > W_cap here, because loading < x_bf and loading > W_cap.
        return BRIDGE_ANCHOR_ACTIVITY + (1.0 - BRIDGE_ANCHOR_ACTIVITY) * (
            (loading - BRIDGE_ANCHOR_LOADING_KG_KG) / (self.x_bf - BRIDGE_ANCHOR_LOADING_KG_KG)
        )

    def film_onset_ceiling_kg(self, dry_matter_kg: float) -> float:
        """The retained-water mass at which a packet deposits a free film."""

        return self.x_bf * dry_matter_kg

    def drawable_floor_kg(self, dry_matter_kg: float) -> float:
        """The retained-water mass the arm may never draw below."""

        return LUIKOV_PARAMS.W_ref * dry_matter_kg

    # ----------------------------------------------------------------------
    # D9-a: the shell-gated forms of the same three quantities.  The two
    # methods above are UNTOUCHED - the kernel's existing call sites read
    # them positionally and their behaviour is byte-identical.
    # ----------------------------------------------------------------------

    @property
    def drain_bracket_end(self) -> str:
        """D9-a-1: which end of the imbibition bracket this arm's token selects."""

        return drain_bracket_end_for_destination(condensate_destination=self.condensate_destination)

    def shell_uptake_capacity_kg(
        self,
        *,
        dry_matter_kg: float,
        retained_water_kg: float,
        front_fraction: float,
    ) -> float:
        """D9-a-3: the arm's uptake ceiling through the open shell, at this x_bf."""

        return uptake_capacity_kg(
            dry_matter_kg=dry_matter_kg,
            retained_water_kg=retained_water_kg,
            front_fraction=front_fraction,
            x_bf=self.x_bf,
        )

    def shell_exchangeable_water_kg(
        self,
        *,
        dry_matter_kg: float,
        retained_water_kg: float,
        front_fraction: float,
    ) -> float:
        """D9-a-3: the drawable stock through the open shell."""

        return exchangeable_water_kg(
            dry_matter_kg=dry_matter_kg,
            retained_water_kg=retained_water_kg,
            front_fraction=front_fraction,
        )

    def shell_exhaustion_floor_kg(
        self,
        *,
        dry_matter_kg: float,
        retained_water_kg: float,
        front_fraction: float,
    ) -> float:
        """D9-a-3: the arm's floor WITH the shell - retained minus exchangeable.

        At ``f = 1`` nothing is drawable and the floor is the whole retained
        water; at ``f = 0`` it falls back to :meth:`drawable_floor_kg` to within
        binary64 resolution, which is the continuity the two forms must have.
        A Law-2 water flux below this floor is the typed phase-exhaustion
        refusal the kernel already raises.
        """

        return math.fsum(
            (
                retained_water_kg,
                -exchangeable_water_kg(
                    dry_matter_kg=dry_matter_kg,
                    retained_water_kg=retained_water_kg,
                    front_fraction=front_fraction,
                ),
            )
        )

    def layer_film_drain_kg(
        self,
        *,
        film_kg: float,
        dry_matter_kg: float,
        retained_water_kg: float,
        front_fraction: float,
    ) -> float:
        """D9-a-2: this arm's drain for one layer, at the end its token selects."""

        return film_drain_kg(
            film_kg=film_kg,
            capacity_kg=self.shell_uptake_capacity_kg(
                dry_matter_kg=dry_matter_kg,
                retained_water_kg=retained_water_kg,
                front_fraction=front_fraction,
            ),
            bracket_end=self.drain_bracket_end,
        )

    def as_dict(self) -> dict:
        """The run-record block for this arm, all of it declared."""

        return {
            "law_id": self.law_id,
            "class": DECLARED_ENGINEERING_ASSUMPTION,
            "x_bf_kg_kg_dry": self.x_bf,
            "x_bf_sweep_kg_kg_dry": list(X_BF_SWEEP_BAND_KG_KG),
            "water_conductance_factor": self.water_conductance_factor,
            "water_conductance_factor_sweep": list(WATER_CONDUCTANCE_FACTOR_SWEEP_BAND),
            "condensate_destination": self.condensate_destination,
            "condensate_destinations": list(CONDENSATE_DESTINATIONS),
            "allocation_basis": "dry_matter",
            "bridged": self.bridged,
            "bridge_anchor_loading_kg_kg_dry": BRIDGE_ANCHOR_LOADING_KG_KG,
            "bridge_anchor_activity": BRIDGE_ANCHOR_ACTIVITY,
            "bridge_join_discontinuity": bridge_join_discontinuity(),
            "exhaustion_floor_kg_kg_dry": LUIKOV_PARAMS.W_ref,
            # D9-a: additive reported keys; the C8 law id above is unchanged.
            "d9a_water_routing_law_id": D9A_WATER_ROUTING_LAW_ID,
            "drain_bracket_end": self.drain_bracket_end,
            "drain_bracket_ends": list(DRAIN_BRACKET_ENDS),
            "front_fraction_domain": list(FRONT_FRACTION_DOMAIN),
            "no_shell_front_fraction_ceiling": NO_SHELL_FRONT_FRACTION_CEILING,
            "hindrance_continuation_start": HINDRANCE_CONTINUATION_START,
            "hindrance_continuation_start_bracket": list(HINDRANCE_CONTINUATION_START_BRACKET),
            "hindrance_continuation_steps": HINDRANCE_CONTINUATION_STEPS,
            "hindrance_continuation_steps_bracket": list(HINDRANCE_CONTINUATION_STEPS_BRACKET),
            "declared_assumptions": {
                name: dict(block) for name, block in DECLARED_ASSUMPTIONS.items()
            },
        }


__all__ = [
    "BRIDGE_ANCHOR_ACTIVITY",
    "BRIDGE_ANCHOR_LOADING_KG_KG",
    "CONDENSATE_DESTINATIONS",
    "CONDENSATE_DESTINATION_DRAIN_BRACKET_END",
    "CONDENSATE_DESTINATION_NOMINAL",
    "D9A_WATER_ROUTING_LAW_ID",
    "DECLARED_ASSUMPTIONS",
    "DECLARED_ENGINEERING_ASSUMPTION",
    "DRAIN_BRACKET_ENDS",
    "DRAIN_BRACKET_END_NOMINAL",
    "DRAIN_BRACKET_END_NONE",
    "EXTERNAL_WATER_NIL_BAND_ACTIVITY_STEP",
    "EXTERNAL_WATER_NIL_BAND_BRACKET_KG",
    "EXTERNAL_WATER_NIL_BAND_KG",
    "EXTERNAL_WATER_NIL_BAND_LAW_ID",
    "FRONT_FRACTION_DOMAIN",
    "HINDRANCE_CONTINUATION_START",
    "HINDRANCE_CONTINUATION_START_BRACKET",
    "HINDRANCE_CONTINUATION_STEPS",
    "HINDRANCE_CONTINUATION_STEPS_BRACKET",
    "KNOWN_ANSWER_ACTIVITY",
    "KNOWN_ANSWER_LOADING_KG_KG",
    "LUIKOV_PARAMS",
    "NO_SHELL_FRONT_FRACTION_CEILING",
    "SORBED_WATER_ARM_LAW_ID",
    "SorbedWaterArm",
    "SorbedWaterDrainBracketRefusal",
    "SorbedWaterFloorRefusal",
    "SorbedWaterFrontRefusal",
    "SorbedWaterLawError",
    "SorbedWaterStockRefusal",
    "WATER_CONDUCTANCE_FACTOR_NOMINAL",
    "WATER_CONDUCTANCE_FACTOR_SWEEP_BAND",
    "X_BF_NOMINAL_KG_KG",
    "X_BF_SWEEP_BAND_KG_KG",
    "bridge_join_discontinuity",
    "drain_bracket_end_for_destination",
    "drain_ledger_residual_kg",
    "drain_split_kg",
    "exchangeable_water_kg",
    "film_drain_kg",
    "open_shell_fraction",
    "shell_capacity_kg",
    "uptake_capacity_kg",
]
