"""Tray-type object (S0) -- typed floor families for the C5a tower.

Implements ``docs/GT_PS2_TRAY_TYPE_OBJECT_DESIGN_BASIS_2026-08-11.md`` (D7 build
order item 4, D7 section 2.3: "C5a therefore models a tray as a typed object
carrying its own geometry and loss law, not one generic perforated floor").
Its governing ruling is
``docs/GT_PS2_TRAY_FAMILY_AXIS_RULING_RECORD_2026-08-11.md`` (owner-ratified
2026-08-11, commit ``4c633a2``).

THE FAMILY AXIS IS VAPOUR-PASSAGE GEOMETRY, NEVER HEATING (ruling items 1-4).
A floor family is defined by free area, passage geometry, gas/steam source and
pressure-loss relation. Family 1 is the UNPERFORATED heated disc: its bore area
is IDENTICALLY ZERO -- a defining predicate of the family, not a typical value
-- so vapour bypasses the bed and the contact area is the tray free surface
(the owner's 2026-08-07 DOME_CONTACT ruling). The Schumacher US4619053A
double-bottom, whose bored distance pieces carry 3-30 % bore area, is a
FAMILY 2 member and is that family's HEATED ARCHETYPE; its steam chest is
weld-sealed from its bores, so heating and vapour passage are two separated
fluid systems in one deck. Indirect heating is therefore an ORTHOGONAL
attribute available to any family (``indirect_heating``), with no
discriminating power: the accepted reference machine carries 140-280 kW of
indirect duty on its perforated MN/SP trays. D7:88's family-1 label carries an
erratum; this module follows the ruling.

PRECURSOR, IMPORTED AND NOT REWRITTEN. ``dtdc_stack`` supplies ``VaporContact``
(a pure passage axis), ``TraySpec`` (whose ``indirect_duty_w`` is an
independent field) and the two area laws. This module adds the four things the
precursor does not carry -- free area, passage geometry, gas/steam source and a
pressure-loss SELECTOR -- and nothing else.

HARD ACCEPTANCE: BYTE-EXACT BEHAVIOURAL PARITY. For all six ``REFERENCE_TRAYS``
the typed dispatch reproduces the precursor's gas-side area law bit for bit
(``==`` on binary64, no tolerance). A generalization that changed any
reference-machine number would be a physics change smuggled in as a refactor
and would reopen the T3a tranche, so the parity test carries a tamper probe
showing it bites.

THE FOUR FAMILIES AND THEIR EVIDENTIARY STANDING (Addendum B ``:32-40`` for
families 1-3; family 4's rows are cited to the primary extraction, because
Addendum B does not discuss family 4 at all):

* Family 1 unperforated heated disc -- no dimensioned anchor is NEEDED: zero
  bore is a ruled geometric predicate. Qualitative anchors: imperforate
  double-walled PD decks in the corpus, the owner's 2026-08-07 PD description,
  Kemper section 4.3.1.
* Family 2 porous / staybolt stripping deck -- dimensioned by Schumacher
  US4619053A (bore area 3-30 %, particularly 4.5-20 %; bored-spacer population
  40-200/m2, particularly 60-140/m2; LOWER BORE DIAMETER not more than about
  25 mm; converging-taper bore area ratio f:g about 2.0-1.5; illustrative
  SPACER pitch and diameter e = 80 mm / a = 50 mm against a printed spacer
  outer diameter of 25-75 mm; 12-30 mm steam gap) and French Oil US5992050A
  (3.175 mm bars, 2.032 mm gaps, 20-30 % open, ~30 deg screen arcs).
  THE PRINTED BORE AND THE PRINTED SPACER ARE DIFFERENT DIMENSIONS. The
  register's compressed phrase "pitch/diameter example 80/50 mm" (Addendum B
  ``:34``) is the SPACER pitch and the SPACER diameter (extraction ``:359``,
  "Illustrative spacer pitch and diameter"); a = 50 mm sits inside the printed
  25-75 mm spacer outer-diameter band (``:362``) and the bore is printed
  separately, capped at about 25 mm (``:360``). Reading 50 mm as a bore is a
  transcription error this module previously made and now refuses to make: the
  characteristic opening of the family-2 reference decks is the printed bore
  maximum, not the spacer diameter.
* Family 3 direct-steam sparge floor -- dimensioned by Desmet US11661564B2
  (~4 mm sparge holes, ~2 mm slots or 15-20 mm holes; 5-10 % PDS GRAVITY
  openings, which is a predesolventizing-gravity figure and NOT a sparge-floor
  open fraction).
* Family 4 dryer / cooler air-distribution floor -- **NO ADMITTED DIMENSIONED
  PATENT ANCHOR**, and still the weakest evidentiary standing of the four
  (``evidentiary_standing_rank`` makes that ordering machine-checkable rather
  than prose). The corpus is NOT silent, and an earlier reading of this module
  said it was: US20100088923A1 (corpus member 35, ``B_subsystems``) prints a
  dimensioned DC air floor -- conventional open area about 0.8 % and improved
  about 2.0 % of tray surface, 3-4 mm round holes, 75-300 mm slot length,
  about 0.25 mm slot width (preferred 0.2-1.0 mm), 200-300 mm water tray
  pressure drop, "approximately halved when open area was nearly doubled"
  (extraction ``:562-579``). Every one of those rows is graded LOW and scoped
  to the DC air floor, i.e. strictly weaker than the MEDIUM tray-construction
  rows that dimension families 2-3.

  **AMENDED 2026-08-12 by owner ruling item 3** (``GT_PS2_OWNER_QUEUE_RULING_2026-08-12.md``:
  *"Admit US20100088923A1 as a LOW-match evidence anchor in its own standing
  tier. It does not qualify target-machine dimensions, coefficients or
  reference-machine membership."*). The rows are now admitted through
  ``LOW_MATCH_DC_AIR_FLOOR_US20100088923A1`` in a **fourth standing tier of its
  own**, and they remain PRINTED FOR REFERENCE AND ENFORCED BY NOTHING:
  ``FAMILY4_ANCHORED_VALUES`` stays empty, no family-4 admissibility band exists
  anywhere in this module, and no reference tray becomes a family-4 member.
  Family 4 keeps ``UNANCHORED_DECLARED`` admitted alongside it, so its weakest
  admissible standing is unchanged and family 4 remains strictly the weakest
  family by rank.

  **As the superseded text asked the amendment to note**: the extraction
  (2026-08-10) predates the ruling (2026-08-11), so owner-ratified caveat 5 was
  ratified with these rows already on disk. Caveat 5's "typed with no anchor" is
  superseded for family 4 by ruling item 3, and by nothing else. The withdrawn
  claim ("no dimensioned anchor exists anywhere in the 47-patent corpus") was
  also cited to Addendum B ``:32-40``, which does not discuss family 4 at all.

Every geometry value in this module is provenance class
``DECLARED_ENGINEERING_ASSUMPTION`` (D0), never
``AUTHENTICATED_MACHINE_DATUM``; ``DeclaredGeometryProvenance`` is the typed
carrier and refuses a silent choice.

THE FOUR RULING CAVEATS, BINDING ON THIS OBJECT (ruling ``:119-139``):

1. **Naming trap.** ``DTDCStackResult.exchange_area_m2`` is populated from
   BED-SIDE areas, not the ruled gas-side area. This module does not propagate
   that name: the typed dispatch is ``gas_side_contact_area_m2`` and there is
   deliberately no ``exchange_area_m2`` symbol here.
2. **~610 is a bed-depth ratio, not a family property.** It is never quoted as
   one and no depth, volume, area or ratio is a field of ``TrayTypeSpec``; the
   ratio moves with loaded depth (about 610 at 0.30 m, about 1830 at 0.90 m)
   and belongs to ``TraySpec`` plus the two area laws.
3. **Open-area bracket discipline.** Kemper prints 1-2 % open area for the
   Schumacher staybolt deck; the patent claims 3-30 % (particularly 4.5-20 %).
   Both bands are stored as printed edges, the admitted family-2 band is their
   interval hull, and they are NEVER averaged into one number. The hull is not
   literally their union -- the printed bands are disjoint, so the interior
   window (2 %, 3 %) is admitted with no printed support; that is disclosed
   here rather than hidden, and narrowing the admitted band is an owner ruling.
4. **The DOME_CONTACT film-coefficient basis remains a disclosed, unqualified
   approximation** (``dtdc_stack`` docstring): dome cells reuse the packed-bed
   Coletto film coefficient at the free surface -- area right and dominant,
   coefficient basis not. This module does not touch it and does not repair it.

LIMITATIONS.

* NO PRESSURE-LOSS ARITHMETIC. S0 carries the loss-law SELECTOR only. The laws
  themselves (the nozzle loss of a direct-steam admission, and the family-2
  nozzle contribution named in the next bullet) are the cell-layer packet's
  work. Per-tray loss is in any case corroborated at the LAYER level and
  family-independently -- EP1336426A1 prints 1-5 kPa per tray, bracketing the
  frozen Ergun band 2.0-5.6 kPa and Zuo 2024's 2.5-6.0 kPa -- so what the typed
  object adds is the passage geometry feeding a law, never a new band.
* THE FAMILY-2 SELECTOR IS PROVISIONAL AGAINST THIS MODULE'S OWN GOVERNING
  RULING. ``CANONICAL_LOSS_RELATION[STRIPPING_DECK] = PACKED_BED_ERGUN`` is one
  relation per family, while ruling ``:53-54`` records that the Schumacher
  archetype's "converging-taper bores with area ratio f:g approx. 2.0-1.5" are
  "a designed nozzle requiring its own loss relation" (extraction ``:364``; the
  band is carried as ``FAMILY2_SCHUMACHER_TAPER_BORE_AREA_RATIO_BAND``). A
  bed-plus-nozzle series is therefore NOT expressible here, and
  ``__post_init__`` refuses a family-2 deck carrying a nozzle selector until an
  owner ruling reopens the lattice. Nothing is computed from the selector in
  S0, so this is a disclosed gap, not a wrong number.
* FAMILIES 2-4 ARE NOT GEOMETRICALLY SEPARABLE ON THE DECLARED FIELDS. Only
  family 1 has a numeric predicate (identically zero bore); families 2, 3 and 4
  all pass vapour through the bed, and their printed passage admissions OVERLAP
  (family 2 admits 2.032-25 mm, family 3 admits the printed alternatives 2 mm,
  4 mm or 15-20 mm, family 4 is unbanded), so a field tuple inside the overlap
  -- a 4 mm opening, say -- can be typed into any of them. The
  discriminator is the family-gated ``ADMITTED_PROVENANCE`` lattice -- a
  declared, reviewable, non-defaultable attribution -- not an inferable
  predicate, and no family-3 open-fraction predicate is invented because the
  only Desmet open-area figure is a predesolventizing-gravity number (see
  ``DESMET_PDS_GRAVITY_OPENING_BAND``).
* WHICH BANDS BITE. Enforced fail-closed in ``__post_init__``: the family-2
  open-area bracket, the family-2 characteristic-opening bracket, the family-2
  bored-spacer population band, the family-3 printed passage alternatives, and
  an impossible-geometry bound on the implied open fraction. PRINTED FOR
  REFERENCE AND ENFORCED BY NOTHING, because no typed field carries the
  dimension: the Schumacher spacer pitch / spacer diameter / spacer outer
  diameter / steam gap / taper ratio, the French Oil bar and gap sizes and open
  band, the Desmet predesolventizing-gravity band, and every family-4 DC row.
  Family 4 carries NO admissibility band at all -- not an oversight: no anchor
  is admitted for it, and inventing a band would breach the bracket discipline.
  No declared-vs-implied consistency invariant is imposed either, because the
  implied fraction assumes circular openings while families 2-3 admit printed
  GAPS and SLOTS; a tolerance would have to be invented, so instead the
  divergence is disclosed and only impossible geometry is refused.
* NO SCHEMA ID OR DIGEST DOMAIN. Sibling modules carry a magic string because
  something hashes it. S0 writes no artefact and no receipt, so it deliberately
  carries no magic/version constant; the first consumer that needs one
  introduces it bound to its own digest rather than inheriting an inert
  placeholder.
* Architecture evidence only. No production wiring, no claim-flag move, no
  F-gate advance, no physics reopened (PHY-008 well-mixed default, the Pe
  0.4-6 bracket and PHY-009 mesh independence are untouched).
* Families 3 and 4 have NO member in the accepted reference machine. They are
  constructible and invariant-checked here, but unexercised by plant-anchored
  geometry.

physically_qualifying is False throughout.
"""

from __future__ import annotations

import enum
import math
import struct
from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import ClassVar

from .dtdc_stack import (
    REFERENCE_TRAYS,
    TraySpec,
    VaporContact,
    bed_side_area_m2,
)


#: Provenance class of every geometry value in this module (D0).
GEOMETRY_PROVENANCE_CLASS = "DECLARED_ENGINEERING_ASSUMPTION"


class TrayTypeError(ValueError):
    """Base exception for a refused tray-type construction or dispatch."""


class FamilyInvariantError(TrayTypeError):
    """Raised when a typed floor violates its family's defining predicates."""


class GeometryProvenanceError(TrayTypeError):
    """Raised when a family carries a provenance its evidence cannot support."""


class TrayTypeMismatchError(TrayTypeError):
    """Raised when a typed object disagrees with the reference ``TraySpec``."""


class UnknownReferenceTrayError(TrayTypeError):
    """Raised when a tray is not a member of the accepted reference machine."""


class TrayFloorFamily(enum.Enum):
    """The four mechanically and hydraulically distinct floor families.

    Member values are D7 section 2.3's enumeration numbers, kept because the
    numbering is cited across the ruling record and the design basis. The axis
    is vapour-passage geometry: nothing here keys on heating.
    """

    UNPERFORATED_HEATED_DISC = 1
    STRIPPING_DECK = 2
    DIRECT_STEAM_SPARGE = 3
    AIR_DISTRIBUTION = 4

    @property
    def d7_family_number(self) -> int:
        return self.value

    @property
    def label(self) -> str:
        return {
            1: "unperforated heated disc (bore area identically zero)",
            2: "porous / staybolt stripping deck (heated archetype Schumacher)",
            3: "direct-steam sparge floor",
            4: (
                "dryer / cooler air-distribution floor "
                "(no ADMITTED dimensioned anchor; LOW-match DC rows only)"
            ),
        }[self.value]


class GasSource(enum.Enum):
    """Where a floor's own gas or steam ORIGINATES, if anywhere.

    ORIGINATION, NOT APERTURE GEOMETRY, and the distinction is the frozen
    ledger's own. FROZEN PHY-024 (``release/physics_decisions.yaml:1556``)
    freezes the per-role vapour routing as "the existing BYPASS (PREDESOLV),
    THROUGH_BED (MAIN), and STEAM_SOURCE (SPARGE) vapor routing", so the steam
    SOURCE node sits on the sparge role alone while the countercurrent decks are
    ascent paths. FROZEN PHY-025 (``:1587``) admits steam at "the SP1
    STEAM_SOURCE node in the frozen PHY-024 network, plus the countercurrent-tray
    apertures through which steam ascends" -- it delegates topology to PHY-024,
    and "up through" is an ascent path, not a second source node. Aperture
    geometry is carried by ``free_area_fraction``, ``characteristic_opening_m``
    and ``contact``, which all three family-2 reference decks carry IDENTICALLY;
    only the origination field differs between MN and SP.
    """

    NONE = "none"
    DIRECT_STEAM = "direct_steam"
    AIR = "air"


class PressureLossRelation(enum.Enum):
    """SELECTOR ONLY. S0 implements no pressure-loss arithmetic whatsoever."""

    DOME_BYPASS = "dome_bypass"
    PACKED_BED_ERGUN = "packed_bed_ergun"
    SPARGE_NOZZLE = "sparge_nozzle"
    AIR_PLENUM = "air_plenum"


class DeclaredGeometryProvenance(enum.Enum):
    """Typed geometry provenance; the D0 disposition forbids silent choices.

    Every member is provenance class ``DECLARED_ENGINEERING_ASSUMPTION`` -- a
    declared parameter, never an authenticated machine datum -- but the members
    do NOT have equal evidentiary standing, and the difference is exposed
    rather than smoothed (``evidentiary_standing``).
    """

    #: Zero bore is a ruled geometric predicate (ruling item 2), so family 1
    #: needs no dimensioned anchor. Qualitative anchors: imperforate
    #: double-walled PD decks in the corpus, the owner's 2026-08-07 PD
    #: description, Kemper section 4.3.1 (ruling caveat 5).
    RULED_ZERO_BORE_PREDICATE = "ruled_zero_bore_predicate_v1"

    #: Schumacher US4619053A -- family 2's heated archetype: bored distance
    #: pieces, bore area 3-30 % (particularly 4.5-20 %), bored-spacer population
    #: 40-200/m2 (particularly 60-140/m2), LOWER BORE DIAMETER not more than
    #: about 25 mm, taper-bore area ratio f:g about 2.0-1.5, illustrative SPACER
    #: pitch/diameter 80/50 mm against a 25-75 mm spacer outer diameter, and a
    #: 12-30 mm steam gap. The bore and the spacer are different dimensions.
    SCHUMACHER_US4619053A = "schumacher_us4619053a_v1"

    #: French Oil US5992050A -- porous bar-screen deck: 3.175 mm bars,
    #: 2.032 mm gaps, 20-30 % open, ~30 deg screen arcs.
    FRENCH_OIL_US5992050A = "french_oil_us5992050a_v1"

    #: Desmet US11661564B2 -- sparge / stripping passages: ~4 mm sparge holes,
    #: ~2 mm slots or 15-20 mm holes, 5-10 % PDS gravity openings.
    DESMET_US11661564B2 = "desmet_us11661564b2_v1"

    #: No ADMITTED dimensioned anchor. Retained after the 2026-08-12 amendment
    #: and still admitted for family 4, which is what keeps family 4's weakest
    #: admissible standing at rank 1 and the family strictly weakest by rank.
    UNANCHORED_DECLARED = "unanchored_declared_v1"

    #: US20100088923A1 (corpus member 35, ``B_subsystems``) -- the dimensioned
    #: dryer/cooler AIR-FLOOR rows, every one graded LOW and scoped to the DC air
    #: floor rather than a DT steam tray. Admitted by OWNER RULING item 3
    #: (2026-08-12) as a LOW-match evidence anchor IN ITS OWN STANDING TIER.
    #:
    #: What "admitted" does NOT mean, in the ruling's own words: it "does not
    #: qualify target-machine dimensions, coefficients or reference-machine
    #: membership". So this member grants NO admissibility band --
    #: ``FAMILY4_ANCHORED_VALUES`` is still empty and the ``FAMILY4_DC_*``
    #: constants are still enforced by nothing -- and no reference tray may
    #: become a family-4 member on its strength. It buys exactly one thing: a
    #: family-4 deck may now name where its declaration came from instead of
    #: declaring into the void, and the standing that naming earns is its own
    #: tier, above a bare declaration and below the ruled geometric predicate.
    LOW_MATCH_DC_AIR_FLOOR_US20100088923A1 = "low_match_dc_air_floor_us20100088923a1_v1"

    @property
    def provenance_class(self) -> str:
        return GEOMETRY_PROVENANCE_CLASS

    @property
    def has_dimensioned_patent_anchor(self) -> bool:
        """A MEDIUM-or-better tray-construction anchor that qualifies dimensions.

        Deliberately FALSE for ``LOW_MATCH_DC_AIR_FLOOR_US20100088923A1``. That
        member's rows are dimensioned on the page, but owner ruling item 3
        admits them while withholding exactly this power -- they qualify no
        target-machine dimension -- and this property is what gates dimensional
        enforcement. Returning True here would silently grant family 4 an
        admissibility band the ruling refuses it.
        """

        return self in (
            DeclaredGeometryProvenance.SCHUMACHER_US4619053A,
            DeclaredGeometryProvenance.FRENCH_OIL_US5992050A,
            DeclaredGeometryProvenance.DESMET_US11661564B2,
        )

    @property
    def has_low_match_reference_anchor(self) -> bool:
        """A named LOW-match source that qualifies nothing (ruling item 3)."""

        return self is DeclaredGeometryProvenance.LOW_MATCH_DC_AIR_FLOOR_US20100088923A1

    @property
    def evidentiary_standing(self) -> str:
        """FOUR distinguishable standings, weakest first in the rank below."""

        if self.has_dimensioned_patent_anchor:
            return "dimensioned_patent_anchor"
        if self is DeclaredGeometryProvenance.RULED_ZERO_BORE_PREDICATE:
            return "ruled_geometric_predicate"
        if self.has_low_match_reference_anchor:
            return "low_match_reference_anchor"
        return "declared_without_anchor"

    @property
    def evidentiary_standing_rank(self) -> int:
        """The standing ORDER, machine-checkable: higher is stronger.

        ``evidentiary_standing`` returns bare strings, so "strictly weaker" was
        prose only and no test could assert the ordering the ruling cares about
        (caveat 5). This makes the comparison executable. The ranks are ordinal
        labels of standing, never weights and never anything a number is
        computed from.

        FOUR TIERS since the 2026-08-12 amendment, weakest first::

            1  declared_without_anchor     -- no named source at all
            2  low_match_reference_anchor  -- a named LOW-match source that
                                              qualifies no dimension (item 3)
            3  ruled_geometric_predicate   -- zero bore, a ruled predicate
            4  dimensioned_patent_anchor   -- qualifies dimensions

        THE PLACEMENT OF TIER 2 IS AN ENACTMENT CHOICE, and it is stated rather
        than buried. The ruling says the anchor gets "its own standing tier" and
        the enactment note says family 4 "keeps the weakest position". A named
        LOW-match source cannot coherently rank BELOW no named source at all --
        that would make printed evidence worsen standing -- so tier 2 sits just
        above the bare declaration. Family 4 keeps the weakest position because
        ``UNANCHORED_DECLARED`` remains admitted for it, so its weakest
        admissible rank is still 1, strictly below every other family.

        Note that tier 1 keeps its number: the amendment inserts a tier, it does
        not renumber the floor, so "weakest rank is 1" survives unchanged.
        """

        if self.has_dimensioned_patent_anchor:
            return 4
        if self is DeclaredGeometryProvenance.RULED_ZERO_BORE_PREDICATE:
            return 3
        if self.has_low_match_reference_anchor:
            return 2
        return 1


# --- printed anchor values, as PRINTED EDGES -------------------------------
# Nothing below is averaged, rounded together or pooled. Each entry is either a
# single printed value or a 2-tuple of the printed band edges, with its source.

#: Family 1: the bore area is IDENTICALLY zero (ruling item 2). A predicate.
FAMILY1_BORE_AREA_FRACTION = 0.0

#: Kemper, Solvent Extraction, printed p. 111 section 4.3.2 -- the Schumacher
#: staybolt deck "generally had a 1-2 % open area".
FAMILY2_KEMPER_OPEN_AREA_BAND: tuple[float, float] = (0.01, 0.02)
#: Schumacher US4619053A col. 10 -- bores are 3-30 % of bottom area.
FAMILY2_PATENT_OPEN_AREA_BAND: tuple[float, float] = (0.03, 0.30)
#: ... "particularly" 4.5-20 %.
FAMILY2_PATENT_PREFERRED_OPEN_AREA_BAND: tuple[float, float] = (0.045, 0.20)
#: THE ADMITTED FAMILY-2 BAND IS THE INTERVAL HULL OF THE TWO PRINTED BANDS,
#: carried as a bracket per ruling caveat 3 (``:127-130``). Both readings are
#: nonzero, so the family argument is unaffected by the disagreement -- and the
#: disagreement is NEVER resolved by averaging. A value outside this band is
#: refused, not clipped. Two honest qualifications. (1) Hull, not union: the
#: printed bands are disjoint, so (0.02, 0.03) is admitted with nothing printed
#: behind it. (2) The corpus holds two further readings that fall INSIDE this
#: bracket and are deliberately not stored as anchors -- EP1336426A1's "2 to 5
#: percent for staybolt holes; 5 to 30 percent for larger openings/screens"
#: (extraction :393, whose own limitation says "do not merge into one value")
#: and US5992050A's "prior-art open area 1 to 3 percent" (:194, match LOW, "not
#: a recommended target value"). They neither widen nor narrow the admitted
#: band: caveat 3 fixes it to the Kemper/Schumacher pair, and any change to it
#: is an owner ruling.
FAMILY2_OPEN_AREA_BRACKET: tuple[float, float] = (0.01, 0.30)

#: Bored-spacer population, Schumacher US4619053A. The patent prints a broad
#: band and a "particularly" sub-band in ONE claim row -- extraction :358,
#: "40 to 200 per m2; particularly 60 to 140 per m2" -- exactly the nested
#: structure of the open-area pair above. Addendum B :34 quotes the broad level
#: and the ruling record :52 quotes the preferred level; that is selective
#: quotation of one source, NOT two divergent readings (this module previously
#: named it as a register-vs-ruling disagreement, which does not exist). Both
#: printed levels are carried, neither is averaged, and the ADMITTED band is the
#: broad printed one -- enforced in ``__post_init__``, refused not clipped.
FAMILY2_PATENT_SPACER_DENSITY_PER_M2_BAND: tuple[float, float] = (40.0, 200.0)
FAMILY2_PATENT_PREFERRED_SPACER_DENSITY_PER_M2_BAND: tuple[float, float] = (
    60.0,
    140.0,
)
#: Schumacher's printed SPACER dimensions (extraction :359 "Illustrative spacer
#: pitch and diameter | e = 80 mm; a = 50 mm", corroborated by :362 "Spacer
#: outer diameter | approximately 25 to 75 mm", inside which a = 50 mm sits).
#: NEITHER IS A BORE. Printed for reference: no typed field carries a spacer
#: dimension, so nothing enforces these. Disclosed implication of the pitch,
#: under its printed identity: an 80 mm square lattice is 156.25 spacers/m2 --
#: inside the broad printed population band, outside the preferred sub-band --
#: and at the printed bore maximum below it implies 7.67 % open area, inside
#: both patent bands. Read instead as a 50 mm BORE it would imply 30.68 %, above
#: the patent's own 30 % ceiling (and 39.27 % at 200 spacers/m2), which is
#: independent arithmetic proof that a = 50 mm is not the bore.
FAMILY2_SCHUMACHER_PITCH_M = 0.080
FAMILY2_SCHUMACHER_SPACER_DIAMETER_M = 0.050
FAMILY2_SCHUMACHER_SPACER_OUTER_DIAMETER_M_BAND: tuple[float, float] = (0.025, 0.075)
#: Schumacher's printed BORE, extraction :360 -- "Lower bore diameter | not more
#: than approximately 25 mm". A printed MAXIMUM, not a worked value: declaring
#: it as the family-2 characteristic opening is an upper-edge choice and is said
#: to be one (see ``_COUNTERCURRENT_DECK``).
FAMILY2_SCHUMACHER_LOWER_BORE_DIAMETER_MAX_M = 0.025
#: Extraction :364 "Taper-bore area ratio | f:g about 2.0 to 1.5", which ruling
#: :53-54 reads as "a designed nozzle requiring its own loss relation". Printed
#: for reference and enforced by nothing: S0 has no throat field and one loss
#: selector per family, so the taper is a disclosed gap (see LIMITATIONS). The
#: wide end of the printed bore follows as 25*sqrt(1.5..2.0) = 30.6-35.4 mm --
#: derived, therefore deliberately NOT stored HERE, and NOT admitted as a
#: characteristic opening: the loss-relevant dimension is the throat.
#:
#: OWNER RULING ITEM 4 (2026-08-12) requires that any such wide-end value "shall
#: be stored only as a derived, provenance-bound nozzle input". That obligation
#: is DISCHARGED, and deliberately not discharged here: the cell-layer packet
#: owns ``cell_closure.DerivedWideEndNozzleInput``, which is CONSTRUCTED ON READ
#: from the throat and this printed band (``is_derived_never_stored``), carries
#: the anchor it came from, and sits in the layer that owns the nozzle loss
#: element ruling item 2 created. Adding a second, STORED copy in S0 would
#: contradict that type's central property and put a derived number in the layer
#: that has no throat field to derive it from.
FAMILY2_SCHUMACHER_TAPER_BORE_AREA_RATIO_BAND: tuple[float, float] = (1.5, 2.0)
#: Schumacher steam gap, extraction :363 "approximately 12 to 30 mm". Printed
#: for reference; no typed field carries a steam gap, so nothing enforces it.
FAMILY2_SCHUMACHER_STEAM_GAP_M_BAND: tuple[float, float] = (0.012, 0.030)
#: French Oil US5992050A porous bar-screen deck (extraction :184-196). Printed
#: for reference; the open band is the printed "typically 20 to 30 percent"
#: sub-band as summarized by Addendum B :36. Enforced by nothing.
FAMILY2_FRENCH_OIL_BAR_M = 0.003175
FAMILY2_FRENCH_OIL_GAP_M = 0.002032
FAMILY2_FRENCH_OIL_OPEN_AREA_BAND: tuple[float, float] = (0.20, 0.30)

#: THE ADMITTED FAMILY-2 PASSAGE-DIMENSION BRACKET: the hull of the two printed
#: family-2 passage dimensions, the French Oil bar-screen gap (2.032 mm) and the
#: Schumacher lower bore maximum (25 mm). Printed edges, refused not clipped,
#: never averaged. It exists because a family-2 deck advertising a dimensioned
#: patent anchor must enforce something from that anchor; without it a 5 m
#: "opening" on a 6 m tray was admissible. Hull, not union, and the interior is
#: admitted with no printed support -- disclosed, as with the open-area bracket.
FAMILY2_OPENING_M_BRACKET: tuple[float, float] = (0.002032, 0.025)

#: THE PRINTED FAMILY-2 TRIPLE IS OVERDETERMINED -- open fraction, bore diameter
#: and spacer population are ONE relation in three quantities, so any two fix the
#: third. Overdetermined means REDUNDANT, not contradictory, and this module
#: previously drew the contradiction conclusion in error: it compared an implied
#: 0.1178 against the patent's 4.5 % lower EDGE instead of the band that edge
#: belongs to, and it computed 0.1178 from the spacer diameter mis-transcribed
#: as a bore. At the printed bore maximum the printed record reconciles: 40, 60,
#: 92, 140, 156.25 and 200 spacers/m2 imply 1.96, 2.95, 4.52, 6.87, 7.67 and
#: 9.82 % open area, so populations from about 61/m2 up sit inside the patent's
#: printed 3-30 % and from about 92/m2 up inside its particularly-claimed
#: 4.5-20 %, while the 40/m2 edge lands on 1.96 %, i.e. on Kemper's printed 2 %.
#: Stated without overstatement: the 40 and 60/m2 edges are BELOW the patent's
#: 3 % floor, so not every printed combination is jointly admissible -- but a
#: whole range is, which is the opposite of "not simultaneously satisfiable".
FAMILY2_ANCHOR_TRIPLE_IS_OVERDETERMINED = True
#: ... and therefore the honest reason the reference decks declare no population:
#: the patent prints population BANDS, not a deck value, and the accepted machine
#: has no certified drawing to pick a point from, so the printed record
#: UNDERCONSTRAINS the population. A typed deck declares the open fraction and
#: the characteristic opening; ``opening_count_density_per_m2`` stays None rather
#: than collapsing a printed band to a scalar by back-calculation. The band
#: itself is not suppressed -- it is carried above, un-collapsed, and enforced
#: whenever a deck does declare a population.
FAMILY2_ANCHOR_TRIPLE_UNDERCONSTRAINS_POPULATION = True

#: Desmet US11661564B2 passages: extraction :277 (sparge-hole diameter about
#: 4 mm), :281 (stripping-tray slot width about 2 mm) and :280 (alternative hole
#: diameter 15-20 mm).
FAMILY3_DESMET_SPARGE_HOLE_M = 0.004
FAMILY3_DESMET_SLOT_M = 0.002
FAMILY3_DESMET_HOLE_M_BAND: tuple[float, float] = (0.015, 0.020)
#: THE ADMITTED FAMILY-3 PASSAGE DIMENSIONS, AS DISJOINT PRINTED ALTERNATIVES --
#: not a hull across them, because the source forbids exactly that: extraction
#: :280's own limitation reads "Alternative to slots; do not combine as the same
#: tray geometry". A family-3 deck therefore declares the printed slot width, OR
#: the printed sparge-hole diameter, OR a hole inside the printed 15-20 mm band;
#: anything between those printed alternatives would be an invented dimension and
#: is refused, never clipped and never averaged. Family 3 advertises a
#: dimensioned patent anchor, so it enforces one. Its OPEN FRACTION stays
#: unbanded: the only Desmet open-area figure is the predesolventizing-gravity
#: one below, and inventing a band would breach the bracket discipline.
FAMILY3_DESMET_ADMITTED_OPENING_M_ALTERNATIVES: tuple[tuple[float, float], ...] = (
    (FAMILY3_DESMET_SLOT_M, FAMILY3_DESMET_SLOT_M),
    (FAMILY3_DESMET_SPARGE_HOLE_M, FAMILY3_DESMET_SPARGE_HOLE_M),
    FAMILY3_DESMET_HOLE_M_BAND,
)

#: NOTE THE ATTRIBUTION, AND NOTE THE MISSING FAMILY PREFIX -- deliberate. The
#: 5-10 % figure (extraction :279, "Predesolventizing tray | Gravity opening
#: area") is for PREDESOLVENTIZING GRAVITY openings, not a sparge-floor open
#: fraction, so it is not used as a family-3 admissibility band and must not be
#: swept up by a consumer filtering constants on a ``FAMILY3_`` prefix. Under the
#: ruled taxonomy a perforated tray in PD duty is a family-2 member in PD
#: service, so this datum belongs to no family namespace.
DESMET_PDS_GRAVITY_OPENING_BAND: tuple[float, float] = (0.05, 0.10)

#: Family 4: STILL deliberately empty after the 2026-08-12 amendment, and that
#: is the load-bearing half of owner ruling item 3. An anchor is now ADMITTED
#: for family 4 (``LOW_MATCH_DC_AIR_FLOOR_US20100088923A1``), but the ruling
#: admits it while withholding the power to qualify dimensions -- "It does not
#: qualify target-machine dimensions, coefficients or reference-machine
#: membership" -- so there is still nothing here to enforce and still no
#: family-4 admissibility band anywhere in this module. Admitting an anchor and
#: granting it dimensional force are different acts; only the first happened.
FAMILY4_ANCHORED_VALUES: tuple[float, ...] = ()
#: US20100088923A1 (corpus member 35, ``B_subsystems``) DC air-floor rows,
#: extraction :562-579, every one graded LOW and scoped to the dryer/cooler air
#: floor rather than a DT steam tray. PRINTED FOR REFERENCE AND ENFORCED BY
#: NOTHING -- unchanged by the amendment, which changed the rows' STANDING, not
#: their force. The two printed open areas are separate printed values
#: (conventional and improved) and are never merged into one; the printed dP
#: pair (200-300 mm water, "approximately halved when open area was nearly
#: doubled") is a loss datum and so belongs to the cell-layer packet, not to S0.
FAMILY4_DC_CONVENTIONAL_OPEN_AREA_FRACTION = 0.008
FAMILY4_DC_IMPROVED_OPEN_AREA_FRACTION = 0.020
FAMILY4_DC_ROUND_HOLE_M_BAND: tuple[float, float] = (0.003, 0.004)
FAMILY4_DC_SLOT_LENGTH_M_BAND: tuple[float, float] = (0.075, 0.300)
FAMILY4_DC_SLOT_WIDTH_M = 0.00025
FAMILY4_DC_SLOT_WIDTH_M_PREFERRED_BAND: tuple[float, float] = (0.0002, 0.0010)


# --- the family lattices ---------------------------------------------------
# Each family owns its contact law, its admitted gas/steam sources and its
# pressure-loss selector (D7 section 2.3: "each with its own free area, passage
# geometry, gas/steam source and pressure-loss relation"). The mappings are
# declared and cross-checked rather than inferred, and a member needing a
# different law is an owner ruling, not a code tweak.

#: DOME_CONTACT is admitted for family 1 ALONE: zero bore is the family-1
#: predicate, so a nonzero-bore floor can never carry the dome law (which is
#: exactly the ~2000x area misclassification the ruling quantifies).
CANONICAL_CONTACT: Mapping[TrayFloorFamily, VaporContact] = MappingProxyType(
    {
        TrayFloorFamily.UNPERFORATED_HEATED_DISC: VaporContact.DOME_CONTACT,
        TrayFloorFamily.STRIPPING_DECK: VaporContact.THROUGH_BED,
        TrayFloorFamily.DIRECT_STEAM_SPARGE: VaporContact.THROUGH_BED,
        TrayFloorFamily.AIR_DISTRIBUTION: VaporContact.THROUGH_BED,
    }
)

#: Family 2 admits DIRECT_STEAM because the frozen ledger puts direct-steam
#: apertures ON the countercurrent decks (PHY-008,
#: ``release/physics_decisions.yaml:955``: "countercurrent trays add a 185 C
#: surface plus direct-steam apertures"). That is what makes the reference
#: sparge tray a family-2 member carrying a direct-steam source rather than a
#: family-3 floor (see ``REFERENCE_TRAY_TYPES``).
#:
#: WHICH family-2 deck carries the source is decided by FROZEN PHY-024, not by
#: that aperture sentence: ``release/physics_decisions.yaml:1556`` freezes "the
#: existing BYPASS (PREDESOLV), THROUGH_BED (MAIN), and STEAM_SOURCE (SPARGE)
#: vapor routing", so the origination node is the sparge role alone and the
#: countercurrent decks are ascent paths. FROZEN PHY-025 (``:1587``) injects at
#: "the SP1 STEAM_SOURCE node in the frozen PHY-024 network, plus the
#: countercurrent-tray apertures through which steam ascends" -- it delegates
#: topology to PHY-024, so "up through" is ascent, not a second source. Hence
#: NONE on MN1/MN2 and DIRECT_STEAM on SP1 is ledger-conformant, and the
#: aperture geometry PHY-008 attributes to countercurrent decks is carried by all
#: three of them identically in the geometry fields.
ADMITTED_GAS_SOURCES: Mapping[TrayFloorFamily, tuple[GasSource, ...]] = (
    MappingProxyType(
        {
            TrayFloorFamily.UNPERFORATED_HEATED_DISC: (GasSource.NONE,),
            TrayFloorFamily.STRIPPING_DECK: (
                GasSource.NONE,
                GasSource.DIRECT_STEAM,
            ),
            TrayFloorFamily.DIRECT_STEAM_SPARGE: (GasSource.DIRECT_STEAM,),
            TrayFloorFamily.AIR_DISTRIBUTION: (GasSource.AIR,),
        }
    )
)

#: ONE selector per family, and family 2's is PROVISIONAL against ruling
#: ``:53-54`` (see LIMITATIONS): the Schumacher archetype's converging-taper
#: bores are a designed nozzle in series with the bed, which this one-relation
#: lattice cannot express and ``__post_init__`` refuses to override. Note also
#: that the selector is chosen by the FAMILY, and families 2-4 are not
#: geometrically separable on the declared fields -- the discriminator is the
#: family-gated ``ADMITTED_PROVENANCE`` lattice below.
CANONICAL_LOSS_RELATION: Mapping[TrayFloorFamily, PressureLossRelation] = (
    MappingProxyType(
        {
            TrayFloorFamily.UNPERFORATED_HEATED_DISC:
                PressureLossRelation.DOME_BYPASS,
            TrayFloorFamily.STRIPPING_DECK:
                PressureLossRelation.PACKED_BED_ERGUN,
            TrayFloorFamily.DIRECT_STEAM_SPARGE:
                PressureLossRelation.SPARGE_NOZZLE,
            TrayFloorFamily.AIR_DISTRIBUTION:
                PressureLossRelation.AIR_PLENUM,
        }
    )
)

#: Provenances a family's evidence can support -- and, because families 2-4
#: share one passage predicate, THE ONLY DISCRIMINATOR between them. A typed
#: floor's family is a declared, reviewable attribution carried by this lattice,
#: never something inferred from its numbers.
ADMITTED_PROVENANCE: Mapping[
    TrayFloorFamily, tuple[DeclaredGeometryProvenance, ...]
] = MappingProxyType(
    {
        TrayFloorFamily.UNPERFORATED_HEATED_DISC: (
            DeclaredGeometryProvenance.RULED_ZERO_BORE_PREDICATE,
        ),
        TrayFloorFamily.STRIPPING_DECK: (
            DeclaredGeometryProvenance.SCHUMACHER_US4619053A,
            DeclaredGeometryProvenance.FRENCH_OIL_US5992050A,
        ),
        TrayFloorFamily.DIRECT_STEAM_SPARGE: (
            DeclaredGeometryProvenance.DESMET_US11661564B2,
        ),
        # AMENDED 2026-08-12, owner ruling item 3. Both are admitted, and the
        # order is weakest-first as everywhere else in this module. Keeping
        # UNANCHORED_DECLARED is not vestigial: it is what holds family 4's
        # weakest admissible standing at rank 1, which is the "family 4 keeps
        # the weakest position" half of the ruling's enactment note.
        TrayFloorFamily.AIR_DISTRIBUTION: (
            DeclaredGeometryProvenance.UNANCHORED_DECLARED,
            DeclaredGeometryProvenance.LOW_MATCH_DC_AIR_FLOOR_US20100088923A1,
        ),
    }
)


def _is_negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _require_binary64(name: str, value: float) -> None:
    if type(value) is not float:
        raise TrayTypeError(f"{name} must be an exact binary64 float")
    if value != value or value in (float("inf"), float("-inf")):
        raise TrayTypeError(f"{name} must be finite")
    if _is_negative_zero(value):
        raise TrayTypeError(f"{name} must use canonical positive zero")


@dataclass(frozen=True, slots=True, kw_only=True)
class TrayTypeSpec:
    """One typed floor: passage geometry, gas source and a loss SELECTOR.

    Deliberately absent (ruling caveat 2): bed depth, bed volume, any exchange
    area and any area RATIO. Those belong to ``TraySpec`` and the two area laws,
    because ~610 is a bed-depth ratio and never a family property.

    Every geometry field is a ``DECLARED_ENGINEERING_ASSUMPTION`` value with its
    printed source recorded in the module constants; none is a machine datum.
    """

    family: TrayFloorFamily
    gas_source: GasSource
    free_area_fraction: float
    characteristic_opening_m: float
    #: None where the anchor prints a population BAND rather than a deck value --
    #: which is the reference case: the printed family-2 record underconstrains
    #: the population, so it is left undeclared rather than back-calculated from
    #: the declared open fraction. When it IS declared, the printed band bites.
    opening_count_density_per_m2: float | None
    pressure_loss_relation: PressureLossRelation
    contact: VaporContact
    #: ORTHOGONAL to the family (ruling item 4): no invariant here reads it.
    indirect_heating: bool
    geometry_provenance: DeclaredGeometryProvenance

    architecture_only: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    geometry_values_are_declared_assumptions: ClassVar[bool] = True
    #: S0 carries the selector; the laws are the cell-layer packet's work.
    pressure_loss_arithmetic_implemented: ClassVar[bool] = False
    loss_relation_is_selector_only: ClassVar[bool] = True
    #: The precursor's gas-side area law is reproduced bit for bit, never
    #: generalized (design basis section 4).
    reference_area_law_preserved_exactly: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if type(self.family) is not TrayFloorFamily:
            raise TrayTypeError("family must be a TrayFloorFamily member")
        if type(self.gas_source) is not GasSource:
            raise TrayTypeError("gas_source must be a GasSource member")
        if type(self.pressure_loss_relation) is not PressureLossRelation:
            raise TrayTypeError(
                "pressure_loss_relation must be a PressureLossRelation member"
            )
        if type(self.contact) is not VaporContact:
            raise TrayTypeError("contact must be a VaporContact member")
        if type(self.geometry_provenance) is not DeclaredGeometryProvenance:
            raise TrayTypeError(
                "geometry_provenance must be a DeclaredGeometryProvenance member"
            )
        if type(self.indirect_heating) is not bool:
            raise TrayTypeError("indirect_heating must be an exact bool")

        _require_binary64("free_area_fraction", self.free_area_fraction)
        _require_binary64("characteristic_opening_m", self.characteristic_opening_m)
        if not 0.0 <= self.free_area_fraction < 1.0:
            raise TrayTypeError("free_area_fraction must lie in [0.0, 1.0)")
        if self.characteristic_opening_m < 0.0:
            raise TrayTypeError("characteristic_opening_m must be non-negative")
        if self.opening_count_density_per_m2 is not None:
            _require_binary64(
                "opening_count_density_per_m2", self.opening_count_density_per_m2
            )
            if self.opening_count_density_per_m2 <= 0.0:
                raise TrayTypeError(
                    "opening_count_density_per_m2 must be None or strictly positive"
                )

        # --- the family lattices, fail closed -------------------------------
        expected_contact = CANONICAL_CONTACT[self.family]
        if self.contact is not expected_contact:
            raise FamilyInvariantError(
                f"family {self.family.d7_family_number} requires "
                f"{expected_contact.value}, not {self.contact.value}"
            )
        admitted_gas = ADMITTED_GAS_SOURCES[self.family]
        if self.gas_source not in admitted_gas:
            raise FamilyInvariantError(
                f"family {self.family.d7_family_number} does not admit gas "
                f"source {self.gas_source.value}"
            )
        expected_loss = CANONICAL_LOSS_RELATION[self.family]
        if self.pressure_loss_relation is not expected_loss:
            raise FamilyInvariantError(
                f"family {self.family.d7_family_number} carries "
                f"{expected_loss.value}, not {self.pressure_loss_relation.value}"
            )
        if self.geometry_provenance not in ADMITTED_PROVENANCE[self.family]:
            raise GeometryProvenanceError(
                f"family {self.family.d7_family_number} cannot carry provenance "
                f"{self.geometry_provenance.value}"
            )

        # --- family 1: zero bore is a DEFINING PREDICATE --------------------
        if self.family is TrayFloorFamily.UNPERFORATED_HEATED_DISC:
            if self.free_area_fraction != FAMILY1_BORE_AREA_FRACTION:
                raise FamilyInvariantError(
                    "family 1 bore area is identically zero (ruled predicate); "
                    "a nonzero-bore heated deck is a family-2 member"
                )
            if self.characteristic_opening_m != 0.0:
                raise FamilyInvariantError(
                    "family 1 has no vapour passage, so its characteristic "
                    "opening is exactly zero"
                )
            if self.opening_count_density_per_m2 is not None:
                raise FamilyInvariantError(
                    "family 1 has no opening population to declare"
                )
        else:
            # The complement predicate: families 2-4 pass vapour.
            if self.free_area_fraction <= 0.0:
                raise FamilyInvariantError(
                    f"family {self.family.d7_family_number} must pass vapour, "
                    "so its free area is strictly positive"
                )
            if self.characteristic_opening_m <= 0.0:
                raise FamilyInvariantError(
                    f"family {self.family.d7_family_number} must declare a "
                    "strictly positive characteristic opening"
                )

        # --- family 2: the printed BRACKETS and BAND, never an average -------
        if self.family is TrayFloorFamily.STRIPPING_DECK:
            low, high = FAMILY2_OPEN_AREA_BRACKET
            if not low <= self.free_area_fraction <= high:
                raise FamilyInvariantError(
                    f"family 2 free area {self.free_area_fraction!r} is outside "
                    f"the Kemper/patent bracket [{low}, {high}] -- refused, "
                    "never clipped and never averaged"
                )
            low, high = FAMILY2_OPENING_M_BRACKET
            if not low <= self.characteristic_opening_m <= high:
                raise FamilyInvariantError(
                    f"family 2 characteristic opening "
                    f"{self.characteristic_opening_m!r} m is outside the printed "
                    f"French Oil gap / Schumacher lower-bore bracket "
                    f"[{low}, {high}] -- refused, never clipped and never averaged"
                )
            if self.opening_count_density_per_m2 is not None:
                low, high = FAMILY2_PATENT_SPACER_DENSITY_PER_M2_BAND
                if not low <= self.opening_count_density_per_m2 <= high:
                    raise FamilyInvariantError(
                        f"family 2 opening population "
                        f"{self.opening_count_density_per_m2!r} per m2 is outside "
                        f"the printed bored-spacer band [{low}, {high}] -- "
                        "refused, never clipped and never averaged"
                    )

        # --- family 3: the printed Desmet passage ALTERNATIVES ---------------
        # Family 3 advertises a dimensioned patent anchor, so it enforces one --
        # as disjoint printed alternatives, because the source says "do not
        # combine as the same tray geometry". Its OPEN FRACTION deliberately
        # stays unbanded: the only Desmet open-area figure is the
        # predesolventizing-gravity one, which is not a sparge-floor open
        # fraction, and inventing a band is forbidden.
        if self.family is TrayFloorFamily.DIRECT_STEAM_SPARGE:
            if not any(
                low <= self.characteristic_opening_m <= high
                for low, high in FAMILY3_DESMET_ADMITTED_OPENING_M_ALTERNATIVES
            ):
                raise FamilyInvariantError(
                    f"family 3 characteristic opening "
                    f"{self.characteristic_opening_m!r} m matches none of the "
                    f"printed Desmet alternatives "
                    f"{FAMILY3_DESMET_ADMITTED_OPENING_M_ALTERNATIVES} -- refused, "
                    "never clipped and never averaged"
                )

        # --- impossible geometry, for any family that declares a population --
        # Not a declared-vs-implied consistency check: no tolerance is invented
        # (the implied fraction assumes circular openings while families 2-3
        # admit printed gaps and slots, so a mismatch is often the formula's, not
        # the deck's). This refuses only what cannot exist -- a population whose
        # own openings would cover the whole tray face -- mirroring the
        # free_area_fraction < 1.0 bound. The circular formula UNDERSTATES a slot
        # or gap, so the bound only ever fires on a genuinely impossible deck.
        implied = self.implied_open_area_fraction
        if implied is not None and implied >= 1.0:
            raise FamilyInvariantError(
                f"declared population {self.opening_count_density_per_m2!r} per m2 "
                f"at a {self.characteristic_opening_m!r} m opening implies an open "
                f"fraction of {implied!r} -- impossible geometry, refused"
            )

    @property
    def bore_area_is_identically_zero(self) -> bool:
        """True only for the family-1 predicate."""

        return self.free_area_fraction == 0.0

    @property
    def implied_open_area_fraction(self) -> float | None:
        """Open fraction implied by the declared opening population, or None.

        Disclosed geometry arithmetic (no loss law): CIRCULAR openings of the
        declared characteristic diameter at the declared population. It exposes
        the RELATION among the three printed family-2 quantities -- any two fix
        the third -- and lets a reader check a declared population against the
        printed open-area bands. It is never used to fill in a missing declared
        value, and it is not evidence of a defect when it differs from the
        declared free area: for a printed GAP (French Oil) or SLOT (Desmet) the
        circular formula is the wrong geometry and understates the passage, and
        no passage-shape field exists to distinguish the cases.
        """

        if self.opening_count_density_per_m2 is None:
            return None
        area = math.pi * self.characteristic_opening_m ** 2 / 4.0
        return self.opening_count_density_per_m2 * area


# --- the accepted reference machine, typed ---------------------------------
# PD1-PD3 are family 1: the owner's 2026-08-07 ruling makes them unperforated
# discs, so their bore area is identically zero and vapour bypasses the bed.
#
# MN1, MN2 and SP1 are family 2 on the same declared reference perforation.
#
# THE SP-TRAY DECISION, and why it is family 2 and not family 3: the ruling's
# family axis is vapour-passage GEOMETRY and never duty, and the accepted
# reference machine gives SP1 the very same perforated THROUGH_BED deck as
# MN1/MN2 with no sparge-hole geometry of its own -- so SP1 is typed a family-2
# stripping deck carrying DIRECT_STEAM as its orthogonal gas-source attribute
# (exactly as indirect heating is orthogonal). The frozen ledger speaks the same
# way and decides WHICH deck carries the source: PHY-008 puts "direct-steam
# apertures" ON the countercurrent decks, and FROZEN PHY-024
# (release/physics_decisions.yaml:1556) freezes STEAM_SOURCE on the SPARGE role
# alone with MAIN as THROUGH_BED, which is exactly the MN-NONE / SP-DIRECT_STEAM
# split below (PHY-025:1587 delegates its topology to that network). Family 3
# stays reserved for a floor whose PASSAGE is the dimensioned Desmet
# sparge-nozzle plate that no reference tray declares -- a reservation carried by
# the ADMITTED_PROVENANCE lattice, since families 2-4 share one passage
# predicate and are not geometrically separable on the declared fields.
#
# Consequence, stated rather than hidden: families 3 and 4 have NO member in the
# reference machine, so the six typed trays exercise families 1 and 2 only.
#
# The declared free area is the PRINTED LOWER EDGE of Schumacher's
# particularly-claimed 4.5-20 % bore band -- written as that band's edge rather
# than as a bare 0.045 literal, so the claim "this is a printed edge" cannot
# drift -- and never an average of the Kemper and patent bands. The
# characteristic opening is the PRINTED MAXIMUM of the lower bore diameter
# ("not more than approximately 25 mm"), which is an upper-edge choice within a
# printed one-sided limit, said here rather than implied; the earlier value
# 0.050 was the printed SPACER diameter and was a transcription error. The
# opening population stays None because the printed record underconstrains it
# (see FAMILY2_ANCHOR_TRIPLE_UNDERCONSTRAINS_POPULATION), not because the patent
# contradicts itself -- at this bore the printed triple broadly reconciles.
# None of these three values can move the reference areas: the gas-side law
# reads contact and TraySpec geometry only.
REFERENCE_DECK_FREE_AREA_FRACTION = FAMILY2_PATENT_PREFERRED_OPEN_AREA_BAND[0]

_PREDESOLVENTIZER_DISC = TrayTypeSpec(
    family=TrayFloorFamily.UNPERFORATED_HEATED_DISC,
    gas_source=GasSource.NONE,
    free_area_fraction=FAMILY1_BORE_AREA_FRACTION,
    characteristic_opening_m=0.0,
    opening_count_density_per_m2=None,
    pressure_loss_relation=PressureLossRelation.DOME_BYPASS,
    contact=VaporContact.DOME_CONTACT,
    indirect_heating=True,
    geometry_provenance=DeclaredGeometryProvenance.RULED_ZERO_BORE_PREDICATE,
)

_COUNTERCURRENT_DECK = TrayTypeSpec(
    family=TrayFloorFamily.STRIPPING_DECK,
    gas_source=GasSource.NONE,
    free_area_fraction=REFERENCE_DECK_FREE_AREA_FRACTION,
    characteristic_opening_m=FAMILY2_SCHUMACHER_LOWER_BORE_DIAMETER_MAX_M,
    opening_count_density_per_m2=None,
    pressure_loss_relation=PressureLossRelation.PACKED_BED_ERGUN,
    contact=VaporContact.THROUGH_BED,
    indirect_heating=True,
    geometry_provenance=DeclaredGeometryProvenance.SCHUMACHER_US4619053A,
)

#: The sparge deck differs from the countercurrent deck in ONE typed attribute.
_SPARGE_DECK = replace(_COUNTERCURRENT_DECK, gas_source=GasSource.DIRECT_STEAM)

REFERENCE_TRAY_TYPES: Mapping[str, TrayTypeSpec] = MappingProxyType(
    {
        "PD1": _PREDESOLVENTIZER_DISC,
        "PD2": _PREDESOLVENTIZER_DISC,
        "PD3": _PREDESOLVENTIZER_DISC,
        "MN1": _COUNTERCURRENT_DECK,
        "MN2": _COUNTERCURRENT_DECK,
        "SP1": _SPARGE_DECK,
    }
)

_REFERENCE_SHAPE: Mapping[str, tuple[str, VaporContact]] = MappingProxyType(
    {t.tray_id: (t.role, t.contact) for t in REFERENCE_TRAYS}
)

if set(REFERENCE_TRAY_TYPES) != set(_REFERENCE_SHAPE):  # pragma: no cover
    raise TrayTypeError(
        "the typed reference table must cover exactly the accepted machine"
    )


def reference_tray_type(tray: TraySpec) -> TrayTypeSpec:
    """The typed floor of an accepted reference tray.

    Refuses a tray that is not a member of the accepted machine, and refuses a
    member whose role or vapour contact has been mutated away from the accepted
    machine -- a mutated passage is a different floor and must be typed
    deliberately, not looked up.

    ``indirect_heating`` is read from the tray's own duty, because heating is an
    orthogonal attribute: zeroing a tray's duty changes that flag and changes
    neither its family nor its area.
    """

    if type(tray) is not TraySpec:
        raise TrayTypeError("tray must be an exact TraySpec")
    try:
        typed = REFERENCE_TRAY_TYPES[tray.tray_id]
        role, contact = _REFERENCE_SHAPE[tray.tray_id]
    except KeyError as error:
        raise UnknownReferenceTrayError(
            f"{tray.tray_id!r} is not a tray of the accepted reference machine"
        ) from error
    if tray.role != role or tray.contact is not contact:
        raise TrayTypeMismatchError(
            f"{tray.tray_id!r} no longer matches the accepted machine "
            f"({role}/{contact.value})"
        )
    return replace(typed, indirect_heating=tray.indirect_duty_w > 0.0)


def require_type_matches_tray(tray: TraySpec, tray_type: TrayTypeSpec) -> None:
    """Fail closed when a typed floor disagrees with its tray's passage.

    The two carry the passage axis independently, so a disagreement is a
    smuggled physics change (the ruling quantifies the misclassification error
    class as ~2000x on gas-side area) and is refused rather than reconciled.
    """

    if type(tray) is not TraySpec:
        raise TrayTypeError("tray must be an exact TraySpec")
    if type(tray_type) is not TrayTypeSpec:
        raise TrayTypeError("tray_type must be an exact TrayTypeSpec")
    if tray.contact is not tray_type.contact:
        raise TrayTypeMismatchError(
            f"tray {tray.tray_id!r} is {tray.contact.value} but its typed floor "
            f"(family {tray_type.family.d7_family_number}) is "
            f"{tray_type.contact.value}"
        )


def gas_side_contact_area_m2(tray: TraySpec, tray_type: TrayTypeSpec) -> float:
    """The area vapour actually sees, dispatched by the TYPED floor.

    NAMING (ruling caveat 1): this is the ruled GAS-SIDE area. The precursor's
    ``DTDCStackResult.exchange_area_m2`` field is populated from BED-SIDE areas,
    and that name is deliberately not propagated into this module.

    The dispatch is the typed object's; the two area primitives are the
    precursor's own, imported and unmodified, so the reference numbers are
    reproduced bit for bit rather than re-derived. Unguarded by design: pair it
    with ``require_type_matches_tray`` (``reference_gas_side_contact_area_m2``
    does both), which is what lets a tamper probe show the parity gate bites.
    """

    if type(tray) is not TraySpec:
        raise TrayTypeError("tray must be an exact TraySpec")
    if type(tray_type) is not TrayTypeSpec:
        raise TrayTypeError("tray_type must be an exact TrayTypeSpec")
    if tray_type.contact is VaporContact.DOME_CONTACT:
        # Family 1 only: vapour bypasses the bed, so the area is the tray free
        # surface. The film coefficient used against it stays the disclosed,
        # unqualified packed-bed approximation (ruling caveat 4).
        return tray.cross_section_m2
    return bed_side_area_m2(tray)


def reference_gas_side_contact_area_m2(tray: TraySpec) -> float:
    """Look up the typed floor, cross-check it, then dispatch."""

    tray_type = reference_tray_type(tray)
    require_type_matches_tray(tray, tray_type)
    return gas_side_contact_area_m2(tray, tray_type)


__all__ = (
    "ADMITTED_GAS_SOURCES",
    "ADMITTED_PROVENANCE",
    "CANONICAL_CONTACT",
    "CANONICAL_LOSS_RELATION",
    "DESMET_PDS_GRAVITY_OPENING_BAND",
    "DeclaredGeometryProvenance",
    "FAMILY1_BORE_AREA_FRACTION",
    "FAMILY2_ANCHOR_TRIPLE_IS_OVERDETERMINED",
    "FAMILY2_ANCHOR_TRIPLE_UNDERCONSTRAINS_POPULATION",
    "FAMILY2_FRENCH_OIL_BAR_M",
    "FAMILY2_FRENCH_OIL_GAP_M",
    "FAMILY2_FRENCH_OIL_OPEN_AREA_BAND",
    "FAMILY2_KEMPER_OPEN_AREA_BAND",
    "FAMILY2_OPENING_M_BRACKET",
    "FAMILY2_OPEN_AREA_BRACKET",
    "FAMILY2_PATENT_OPEN_AREA_BAND",
    "FAMILY2_PATENT_PREFERRED_OPEN_AREA_BAND",
    "FAMILY2_PATENT_PREFERRED_SPACER_DENSITY_PER_M2_BAND",
    "FAMILY2_PATENT_SPACER_DENSITY_PER_M2_BAND",
    "FAMILY2_SCHUMACHER_LOWER_BORE_DIAMETER_MAX_M",
    "FAMILY2_SCHUMACHER_PITCH_M",
    "FAMILY2_SCHUMACHER_SPACER_DIAMETER_M",
    "FAMILY2_SCHUMACHER_SPACER_OUTER_DIAMETER_M_BAND",
    "FAMILY2_SCHUMACHER_STEAM_GAP_M_BAND",
    "FAMILY2_SCHUMACHER_TAPER_BORE_AREA_RATIO_BAND",
    "FAMILY3_DESMET_ADMITTED_OPENING_M_ALTERNATIVES",
    "FAMILY3_DESMET_HOLE_M_BAND",
    "FAMILY3_DESMET_SLOT_M",
    "FAMILY3_DESMET_SPARGE_HOLE_M",
    "FAMILY4_ANCHORED_VALUES",
    "FAMILY4_DC_CONVENTIONAL_OPEN_AREA_FRACTION",
    "FAMILY4_DC_IMPROVED_OPEN_AREA_FRACTION",
    "FAMILY4_DC_ROUND_HOLE_M_BAND",
    "FAMILY4_DC_SLOT_LENGTH_M_BAND",
    "FAMILY4_DC_SLOT_WIDTH_M",
    "FAMILY4_DC_SLOT_WIDTH_M_PREFERRED_BAND",
    "FamilyInvariantError",
    "GEOMETRY_PROVENANCE_CLASS",
    "GasSource",
    "GeometryProvenanceError",
    "PressureLossRelation",
    "REFERENCE_DECK_FREE_AREA_FRACTION",
    "REFERENCE_TRAY_TYPES",
    "TrayFloorFamily",
    "TrayTypeError",
    "TrayTypeMismatchError",
    "TrayTypeSpec",
    "UnknownReferenceTrayError",
    "gas_side_contact_area_m2",
    "reference_gas_side_contact_area_m2",
    "reference_tray_type",
    "require_type_matches_tray",
)
