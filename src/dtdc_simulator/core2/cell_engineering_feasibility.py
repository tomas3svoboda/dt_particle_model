"""Bounded carrier-topology cell closures for engineering-feasibility marching.

This module is deliberately isolated from :mod:`cell_closure`.  It implements
the ruled ``BINARY_NO_INERT`` / two-external-liquid branch, the separately
typed corrected ``POSITIVE_CARRIER`` branch, and an explicit shared-old-wall
update.  It is not an F4 implementation, physical qualification, plant
calibration, or production wiring.

The nonlinear block has eight intrinsic-binary unknowns.  Its Jacobian is
assembled by forward analytic differentiation; no finite difference enters
either solver.  Backtracking may shorten a Newton step after an inadmissible
trial, but no trial is clipped, projected, normalised, or switched between
carrier laws.
"""

from __future__ import annotations

import enum
import hashlib
import math
import struct
from dataclasses import InitVar, dataclass, replace
from types import MappingProxyType
from typing import ClassVar

import numpy as np

from . import cell_closure as cc
from . import cell_native_caloric as nc

# D9-b (owner ruling 2026-09-13): the dry-shell receding-front interface reads
# its three declared constants from the modules that own them.  None of the
# three imports this one, so no cycle is created (checked on this tree):
#
#   particle/engineering_foundation.py:47  meal conductivity 0.24 W/(m K),
#                                          declared sensitivity 0.29 at :48
#   sorption_interface.py:134              particle radius 0.885e-3 m
#   through_bed_sorbed_water_law.py        the declared no-shell front ceiling
from . import sorption_interface as _sorption_interface
from . import through_bed_sorbed_water_law as _sorbed_water
from .particle import engineering_foundation as _particle_foundation
from .props import binary_gas
from .props import hexane as hexane_props
from .props import water as water_props

#: D9-b-1: the shell's heat conductance constants, CALLER-FREE because both
#: are already declared elsewhere and read here, never re-declared.
DRY_SHELL_MEAL_CONDUCTIVITY_W_M_K = _particle_foundation.NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K
DRY_SHELL_MEAL_CONDUCTIVITY_SENSITIVITY_W_M_K = (
    _particle_foundation.STRUCTURAL_CONDUCTIVITY_SENSITIVITY_W_M_K
)
DRY_SHELL_PARTICLE_RADIUS_M = _sorption_interface.ParticleSorptionParams().particle_radius_m
#: D9-b-2 as amended by owner ruling D9-e (2026-09-13): the declared front
#: fraction AT which there is no shell, 1.0 exactly.  The withdrawn 0.99
#: boundary lived in the same name and is read from the same module.
NO_SHELL_FRONT_FRACTION_CEILING = _sorbed_water.NO_SHELL_FRONT_FRACTION_CEILING
#: D9-e: the declared hindrance continuation that SEEDS the dry-shell branch
#: when a direct Newton refuses.  Both are read from the law module's rows
#: ``d9e_hindrance_continuation_start`` and ``..._steps``, never re-declared.
DRY_SHELL_HINDRANCE_CONTINUATION_START = _sorbed_water.HINDRANCE_CONTINUATION_START
DRY_SHELL_HINDRANCE_CONTINUATION_STEPS = _sorbed_water.HINDRANCE_CONTINUATION_STEPS

#: D9-e FIX (2026-09-14): the law id of the ONE declared row this fix adds.
D9E_CARRIED_SEED_ROUTING_LAW_ID = "d9e-carried-cell-field-state-seed-width-from-record-v1"
#: D9-e FIX (2026-09-14): the row itself.  It is READ by the kernel that
#: carries the seed (``sp1_k2_law2_tray_integration.D9E_CARRIED_SEED_WIDTH_ROW``)
#: and IMPLEMENTED here, in the entry point and in the evaluator, because the
#: cell owns both widths.  It declares no number: it is a routing rule whose
#: two values are the two unknown-name tuples this module already exports.
D9E_CARRIED_SEED_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "d9e_carried_seed_width_from_record": MappingProxyType(
            {
                "class": "DECLARED_ENGINEERING_ASSUMPTION",
                "law_id": D9E_CARRIED_SEED_ROUTING_LAW_ID,
                "criterion": (
                    "the RECORD decides which interface a cell step solves; the "
                    "seed decides only where the Newton starts.  Owner ruling "
                    "D9-e (2026-09-13) made the ten-unknown dry-shell state an "
                    "ACCEPTED fast-block state, so a tray integration now carries "
                    "it into the next interval as that layer's seed - the shipped "
                    "warm start, unchanged.  A layer can LEAVE the dry-shell "
                    "branch between intervals, its packets regaining an attached "
                    "hexane film, an external water film, or both, and the record "
                    "then presents the shipped eight-unknown state while the "
                    "carried seed is ten wide.  The two extra unknowns - the "
                    "front temperature and the surface hexane activity - are "
                    "properties of the dry-shell FORMULATION and carry no "
                    "inventory, so leaving the branch drops them and keeps the "
                    "eight, which is exactly the state the shipped law would have "
                    "carried; the projection is the exact inverse of the "
                    "extension the branch already makes in the other direction"
                ),
                "value": (
                    "the solved width is len(DRY_SHELL_UNKNOWN_NAMES) when "
                    "presents_dry_shell_receding_front_interface(inputs) is true "
                    "and len(BINARY_NO_INERT_UNKNOWN_NAMES) otherwise, whatever "
                    "the seed's width; a carried ten-unknown seed handed to a "
                    "record that no longer presents the branch is PROJECTED onto "
                    "its first eight components and nothing else about the step "
                    "changes"
                ),
                "bracket": (
                    "none: the row has no free value.  The two widths are the two "
                    "unknown-name tuples already declared in this module, and the "
                    "branch predicate is the one the cell already dispatches on"
                ),
                "rejecting_outcome": (
                    "a step whose record does not present the dry-shell interface "
                    "solving ten rows anyway.  That is how the defect showed on "
                    "the three D9-e honest farm lanes of commit 5fc487d: the "
                    "ten-row branch read the hexane closure's front fraction off "
                    "a record that binds the hexane FILM law, where that closure "
                    "is None, and the march died with an AttributeError.  The "
                    "evaluator's gate below makes that state a typed refusal and "
                    "the entry point's projection makes it unreachable"
                ),
            }
        )
    }
)

UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.314462618

# Ruled molar-mass role split.  The FSG-local value reaches only A_FSG;
# physical state, density, flux/mass conversion and residuals use the canonical
# repository value.
FSG_HEXANE_MOLAR_MASS_KG_MOL = 0.08618
FSG_WATER_MOLAR_MASS_KG_MOL = 0.018015268
CANONICAL_HEXANE_MOLAR_MASS_KG_MOL = hexane_props.M
CANONICAL_WATER_MOLAR_MASS_KG_MOL = water_props.M

LAW2_TEMPERATURE_DOMAIN_K = (330.0, 435.0)

# Frozen single-tray certification box retained as historical evidence.  The
# historical six-layer structural admission is also retained exactly.  The
# separate arbitrary-K admission reaches the 101325 Pa downstream boundary so
# a refined top numerical cell can have less than 1000 Pa of bed drop.  It is a
# noncertifying engineering runtime domain only: it extends neither saturation
# / CERT-01 evidence nor F4 or physical qualification.
FROZEN_SINGLE_TRAY_LAW2_PRESSURE_CERTIFICATION_BOX_PA = (102_325.0, 106_925.0)
ENGINEERING_SIX_LAYER_SERIAL_LAW2_PRESSURE_ADMISSION_PA = (102_325.0, 134_925.0)
ENGINEERING_ARBITRARY_K_SERIAL_LAW2_PRESSURE_ADMISSION_PA = (101_325.0, 134_925.0)
LAW2_PRESSURE_DOMAIN_PA = ENGINEERING_ARBITRARY_K_SERIAL_LAW2_PRESSURE_ADMISSION_PA
ARBITRARY_K_PRESSURE_ADMISSION_EXTENDS_SATURATION_CERT01 = False
ARBITRARY_K_PRESSURE_ADMISSION_EXTENDS_F4 = False
ARBITRARY_K_PRESSURE_ADMISSION_PHYSICALLY_QUALIFYING = False
LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S = (-5.0, 5.0)
LAW2_DIFFUSIVITY_MULTIPLIER_DOMAIN = (0.8, 1.2)

LAW1_TEMPERATURE_DOMAIN_K = LAW2_TEMPERATURE_DOMAIN_K
LAW1_PRESSURE_DOMAIN_PA = ENGINEERING_SIX_LAYER_SERIAL_LAW2_PRESSURE_ADMISSION_PA

POSITIVE_CARRIER_UNKNOWN_NAMES = (
    "gas_temperature_k",
    "gas_hexane_mole_fraction",
    "gas_water_mole_fraction",
    "layer_pressure_pa",
    "interface_temperature_k",
    "external_hexane_saturation",
    "external_water_saturation",
)
POSITIVE_CARRIER_ROW_NAMES = (
    "gas_hexane_balance",
    "gas_water_balance",
    "gas_energy_balance",
    "layer_pressure_balance",
    "interface_energy_balance",
    "film_hexane_identity",
    "film_water_identity",
)
POSITIVE_CARRIER_COLUMN_SCALES = (100.0, 1.0, 1.0, 1.0e4, 100.0, 1.0, 1.0)

BINARY_NO_INERT_UNKNOWN_NAMES = (
    "gas_temperature_k",
    "gas_hexane_mole_fraction",
    "layer_pressure_pa",
    "interface_temperature_k",
    "interface_hexane_mole_fraction",
    "external_hexane_saturation",
    "external_water_saturation",
    "total_molar_flux_mol_m2_s",
)
BINARY_NO_INERT_ROW_NAMES = (
    "gas_hexane_balance",
    "gas_energy_balance",
    "layer_pressure_balance",
    "interface_energy_balance",
    "film_hexane_identity",
    "film_water_identity",
    "hexane_fugacity_equality",
    "water_fugacity_equality",
)
BINARY_NO_INERT_COLUMN_SCALES = (100.0, 1.0, 1.0e4, 100.0, 1.0, 1.0, 1.0, 1.0)

#: D9-b-1 (owner ruling 2026-09-13): the dry-shell receding-front interface.
#: The shipped eight unknowns PLUS the front temperature and the surface
#: hexane activity, in the order ``ft2a1_shell.UNKNOWN_NAMES`` fixes.
DRY_SHELL_UNKNOWN_NAMES = (
    *BINARY_NO_INERT_UNKNOWN_NAMES,
    "front_temperature_k",
    "surface_hexane_activity",
)
#: The ten rows, in the order ``ft2a1_shell.ROW_NAMES`` fixes.  Rows 1-3, 5
#: and 6 are the shipped rows verbatim; row 4 is the surface energy balance in
#: its RECOMBINED form (convection against the front's evaporation duty and
#: the condensing water's latent release); row 7 is the shipped hexane
#: isofugacity row with the surface activity as an UNKNOWN; row 8 is - SINCE
#: OWNER RULING D9-e (2026-09-13) - the shipped water isofugacity row at the
#: shell's SORBED activity, divided through by the pressure; rows 9 and 10 are
#: the shell mass balance and the front energy row carried as Fourier's law
#: solved for the DROP.  Until D9-e row 8 was the binary identity with the
#: water arm blocked (``binary_identity_water_blocked``), which the ruling
#: withdrew.
DRY_SHELL_ROW_NAMES = (
    "gas_hexane_balance",
    "gas_energy_balance",
    "layer_pressure_balance",
    "surface_energy_convection_equals_evaporation_w",
    "film_hexane_identity",
    "film_water_identity",
    "surface_hexane_fugacity_equality",
    "water_fugacity_equality_at_sorbed_activity",
    "shell_mass_balance",
    "front_energy_fourier_drop_k",
)
#: The shipped column scales with the front temperature scaled like the other
#: two temperatures and the surface activity like the other order-one
#: unknowns (``ft2a1_shell.COLUMN_SCALES``).
DRY_SHELL_COLUMN_SCALES = (*BINARY_NO_INERT_COLUMN_SCALES, 100.0, 1.0)
#: D9-e (owner ruling 2026-09-13): the two ends of the water-arm homotopy the
#: dry-shell branch carries.  ``ACTIVE`` is the RULED PHYSICS - the shipped
#: two-film water arm at unit conductance, the water isofugacity row at the
#: shell's sorbed activity, the water latent heat in the surface energy row -
#: and is the value every reported root is solved and revalidated at.
#: ``BLOCKED`` is the WITHDRAWN water-blocked form (closure (v)).  It is no
#: longer any solver's route: the D9-b water-blocking walk that once reached it
#: was RETIRED at D9-e level 0, measured unable to start from the caller's seed
#: on any state it was given - the blocked end is exactly what that seed never
#: reached, which is why the walk ran the other way before the amendment.  The
#: value is kept only so that a test or a diagnostic can evaluate the withdrawn
#: closure and show what was withdrawn.  Neither end is a physical quantity:
#: 1.0 is the closure the cell solves.
DRY_SHELL_WATER_ARM_ACTIVE = 1.0
DRY_SHELL_WATER_ARM_BLOCKED = 0.0

#: The TRIAL domain of the surface hexane activity, READ from
#: ``docs/evidence/ft2a1_closure_2026-09-12/ft2a1_shell.py``
#: (``SURFACE_ACTIVITY_DOMAIN``) and not minted here: the activity is
#: physically at or below one, and the trial domain is opened to 2.0 rather
#: than 1.0 so that a root ABOVE saturation is REPORTED rather than walled off
#: by the trial bound.  A root above one is a measurement to explain, never a
#: state this cell claims.
DRY_SHELL_SURFACE_ACTIVITY_DOMAIN = (0.0, 2.0)


class EngineeringFeasibilityError(cc.CellClosureError):
    """Base class for the bounded engineering-only closure."""


class EngineeringConfigurationError(EngineeringFeasibilityError):
    """Declared inputs contradict the bounded engineering contract."""


class CarrierTopologyConfigurationError(EngineeringConfigurationError):
    """The serialized carrier topology and inlet composition disagree."""


class PositiveCarrierTopologyRefusal(EngineeringConfigurationError):
    """A Law-2-only entry was given the separately typed Law-1 topology."""


class UnsupportedInventoryBranchRefusal(EngineeringConfigurationError):
    """A permanently refused Law-2 inventory branch was selected."""


class UnsupportedStatefulClosureRefusal(EngineeringConfigurationError):
    """Stateful gas or film ownership was requested in this bounded tranche."""


class EngineeringStepRejected(EngineeringFeasibilityError):
    """A trial state left the declared domain; it was not repaired."""


class PositiveCarrierInterfaceCrossoverRefusal(EngineeringStepRejected):
    """The positive-carrier interface reached ``g_I <= 0`` exactly."""


class EngineeringNewtonConvergenceError(EngineeringFeasibilityError):
    """The damped analytic Newton solve did not converge."""


class AcceptedStepIntegrityError(EngineeringFeasibilityError):
    """An accepted fast-step object was forged, substituted, or mismatched."""


class CarrierTopology(enum.Enum):
    POSITIVE_CARRIER = "positive_carrier"
    BINARY_NO_INERT = "binary_no_inert"


class Law2CoefficientModel(enum.Enum):
    COLETTO_NOMINAL = "coletto_nominal"


class Law1ConductanceAuthority(enum.Enum):
    """Authority identity for ``kappa = c*k`` on the positive-carrier law."""

    DECLARED_POSITIVE_CARRIER = "declared_positive_carrier_conductance"


class EngineeringCarrierLaw(enum.Enum):
    """Exact pre-solve carrier-law selection; this enum does not execute a solve."""

    POSITIVE_CARRIER_LAW1 = "positive_carrier_law1"
    BINARY_NO_INERT_LAW2 = "binary_no_inert_law2"


def _require_finite(name: str, value: float) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise EngineeringConfigurationError(f"{name} must be a finite binary64 float")


def _require_positive(name: str, value: float) -> None:
    _require_finite(name, value)
    if value <= 0.0:
        raise EngineeringConfigurationError(f"{name} must be strictly positive")


def _require_nonnegative(name: str, value: float) -> None:
    _require_finite(name, value)
    if value < 0.0:
        raise EngineeringConfigurationError(f"{name} must be nonnegative")


def _require_closed_fraction(name: str, value: float) -> None:
    _require_finite(name, value)
    if not 0.0 <= value <= 1.0:
        raise CarrierTopologyConfigurationError(f"{name} must lie in the closed interval [0, 1]")


def _fsg_prefactor_a() -> float:
    """Return A_FSG using only the ruled FSG-local mass convention."""

    mass_term = (
        1.0 / (FSG_HEXANE_MOLAR_MASS_KG_MOL * 1000.0) + 1.0 / (FSG_WATER_MOLAR_MASS_KG_MOL * 1000.0)
    ) ** 0.5
    diffusion_volume_term = ((6.0 * 16.5 + 14.0 * 1.98) ** (1.0 / 3.0) + 12.7 ** (1.0 / 3.0)) ** 2
    return 1.00e-3 * 1.0e-4 * 101_325.0 * mass_term / diffusion_volume_term


FSG_DIFFUSIVITY_PREFACTOR = _fsg_prefactor_a()


def _binary64_bits(value: float) -> bytes:
    return struct.pack(">d", value)


# --- CAMP-13 transport-correlation authorities (owner batch 2026-08-23) -----
#
# E-8 (CB-3) RULED 2026-08-23: "Yes - bind the existing values, select no new
# number."  The FSG sub-line was RULED verbatim: "the live `_fsg_prefactor_a()`
# function is declared the authority. The `Revision F` annex remains the
# archived derivation and is NOT frozen by this ruling."  The two authority
# declarations below therefore BIND the exact pre-existing executable values.
# Every declared field is proven bit-identical to the live expression at
# construction time, and the K4 declared-engineering model-factory identity
# folds both declarations.  No number is selected, changed, or rounded here.

FSG_DHW_AUTHORITY_DECISION_ID = "E-8-FSG-LIVE-FUNCTION-AUTHORITY-RULED-2026-08-23"
FSG_DHW_AUTHORITY_ID = "live-fsg-prefactor-a-printed-atm-convention-dhw-authority-v1"
#: The bare Fuller-Schettler-Giddings exponent formerly inlined on the K4
#: Law-2 path (``t_interface**1.75``); bound at its existing value by E-8.
FSG_TEMPERATURE_EXPONENT = 1.75
FSG_PRINTED_PREFACTOR = 1.00e-3
FSG_CM2_PER_S_TO_M2_PER_S = 1.0e-4
FSG_REFERENCE_ATMOSPHERE_PA = 101_325.0
FSG_DIFFUSION_VOLUME_CARBON = 16.5
FSG_DIFFUSION_VOLUME_HYDROGEN = 1.98
FSG_DIFFUSION_VOLUME_WATER = 12.7
FSG_HEXANE_CARBON_ATOM_COUNT = 6.0
FSG_HEXANE_HYDROGEN_ATOM_COUNT = 14.0

#: Coletto eq. B.7 packed-bed Nusselt constants and the Chilton-Colburn
#: analogy exponent, formerly bare literals on the K4 execution path.
#:
#: D20, 2026-09-23.  Owner ruling of 2026-09-23, verbatim, "I rule D19 and
#: D20 as recommended"; packet
#: docs/GT_PS2_RULING_PACKET_D20_VOIDAGE_REYNOLDS_2026-09-22.md option 1;
#: build record docs/GT_PS2_D20_LEVEL0_BUILD_RECORD_2026-09-23.md.  The
#: prefactor is now the one the primary source prints: Faner (2008) PhD
#: thesis, printed page 109, Equation 4.27,
#:
#:     Nu_eps = 0.6941 Re_eps^0.579 Pr^(1/3),   R^2 = 0.929,
#:
#: read from an enlargement of the printed line and confirmed independently
#: by the drawn correlation line of the source's own Figure 4.26 (RMS 1.02
#: native pixels over 21 dashes).  The digitized source, its sha256 and the
#: page-by-page reading are pinned by
#: paper/analysis/datasets/faner2008_thesis_4_5_3/DATASET_RECORD.md and its
#: source_reading.md.  No number is selected here; the value is transcribed.
#: WAS, until 2026-09-23: 0.6949, carried from Coletto (2022) eq. B.7 before
#: the thesis was held.  Nothing on the thesis pages prints 0.6949; the
#: difference is 0.115 per cent on the prefactor.
COLETTO_B7_NUSSELT_PREFACTOR = 0.6941
COLETTO_B7_REYNOLDS_EXPONENT = 0.579
COLETTO_B7_PRANDTL_EXPONENT = 1.0 / 3.0
CHILTON_COLBURN_ANALOGY_EXPONENT = 2.0 / 3.0
LAW2_FILM_CORRELATION_DECISION_ID = "E-8-CB3-TRANSPORT-CORRELATION-BINDING-RULED-2026-08-23"
LAW2_FILM_CORRELATION_AUTHORITY_ID = "coletto-b7-nusselt-chilton-colburn-analogy-authority-v1"
#: G-4 carries the open primary-source obligation forward UNCHANGED
#: (docs/evidence/GT_PS2_THESIS_DEPENDENCY_SWEEP_FANER2008_2026-08-22.md).
COLETTO_B7_PRIMARY_SOURCE_CITATION_STATUS = "CITETODO-faner-2008-open"


@dataclass(frozen=True, slots=True, kw_only=True)
class FSGDhwDiffusivityAuthority:
    """E-8-ruled ``D_hw`` authority: the live :func:`_fsg_prefactor_a` law.

    The live function is the declared authority; this record is its typed,
    digest-foldable declaration.  Construction re-derives the prefactor from
    the declared fields in the function's exact operation order and refuses
    unless the result is bit-identical to ``_fsg_prefactor_a()`` and to the
    module constant ``FSG_DIFFUSIVITY_PREFACTOR``.  The archived Revision-F
    annex derivation is reconciled by KAT (printed-convention headline
    ``1.252487038212464e-05 m2/s`` at 334.4803981593741 K / 101325 Pa,
    recorded rel. diff 1.4e-16 in
    ``docs/GT_PS2_DHW_FSG_CONVENTION_VERIFICATION_2026-08-21.md``).
    """

    decision_id: str = FSG_DHW_AUTHORITY_DECISION_ID
    authority_id: str = FSG_DHW_AUTHORITY_ID
    prefactor_a_atm_convention: float = FSG_DIFFUSIVITY_PREFACTOR
    temperature_exponent: float = FSG_TEMPERATURE_EXPONENT
    printed_prefactor: float = FSG_PRINTED_PREFACTOR
    cm2_per_s_to_m2_per_s: float = FSG_CM2_PER_S_TO_M2_PER_S
    reference_atmosphere_pa: float = FSG_REFERENCE_ATMOSPHERE_PA
    diffusion_volume_carbon: float = FSG_DIFFUSION_VOLUME_CARBON
    diffusion_volume_hydrogen: float = FSG_DIFFUSION_VOLUME_HYDROGEN
    diffusion_volume_water: float = FSG_DIFFUSION_VOLUME_WATER
    fsg_hexane_molar_mass_kg_mol: float = FSG_HEXANE_MOLAR_MASS_KG_MOL
    fsg_water_molar_mass_kg_mol: float = FSG_WATER_MOLAR_MASS_KG_MOL

    live_function_is_declared_authority: ClassVar[bool] = True
    revision_f_annex_frozen: ClassVar[bool] = False
    revision_f_annex_is_archived_derivation: ClassVar[bool] = True
    selects_new_number: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        mass_term = (
            1.0 / (self.fsg_hexane_molar_mass_kg_mol * 1000.0)
            + 1.0 / (self.fsg_water_molar_mass_kg_mol * 1000.0)
        ) ** 0.5
        diffusion_volume_term = (
            (
                FSG_HEXANE_CARBON_ATOM_COUNT * self.diffusion_volume_carbon
                + FSG_HEXANE_HYDROGEN_ATOM_COUNT * self.diffusion_volume_hydrogen
            )
            ** (1.0 / 3.0)
            + self.diffusion_volume_water ** (1.0 / 3.0)
        ) ** 2
        rebuilt = (
            self.printed_prefactor
            * self.cm2_per_s_to_m2_per_s
            * self.reference_atmosphere_pa
            * mass_term
            / diffusion_volume_term
        )
        declared_bits = _binary64_bits(self.prefactor_a_atm_convention)
        if (
            self.decision_id != FSG_DHW_AUTHORITY_DECISION_ID
            or self.authority_id != FSG_DHW_AUTHORITY_ID
            or _binary64_bits(rebuilt) != declared_bits
            or _binary64_bits(_fsg_prefactor_a()) != declared_bits
            or _binary64_bits(FSG_DIFFUSIVITY_PREFACTOR) != declared_bits
            or _binary64_bits(self.temperature_exponent) != _binary64_bits(1.75)
            or _binary64_bits(self.fsg_hexane_molar_mass_kg_mol)
            != _binary64_bits(FSG_HEXANE_MOLAR_MASS_KG_MOL)
            or _binary64_bits(self.fsg_water_molar_mass_kg_mol)
            != _binary64_bits(FSG_WATER_MOLAR_MASS_KG_MOL)
        ):
            raise EngineeringConfigurationError(
                "the FSG D_hw authority declaration drifted from the live " "_fsg_prefactor_a() law"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class Law2FilmCorrelationAuthority:
    """E-8-ruled binding of Coletto B.7 and the Chilton-Colburn analogy.

    Binds the packed-bed Nusselt law
    ``Nu_eps = 0.6941 * Re_eps**0.579 * Pr**(1/3)`` and the Chilton-Colburn
    mass-transfer analogy exponent ``(Pr/Sc)**(2/3)`` on the K4 Law-2
    execution path.  No number is selected: under the D20 ruling of
    2026-09-23 the prefactor is transcribed from the primary source, Faner
    (2008) printed page 109 Equation 4.27; it read 0.6949 until then.  The
    CITETODO status string is carried forward unchanged, because closing it
    is a separate decision that this ruling did not take.
    """

    decision_id: str = LAW2_FILM_CORRELATION_DECISION_ID
    authority_id: str = LAW2_FILM_CORRELATION_AUTHORITY_ID
    nusselt_prefactor: float = COLETTO_B7_NUSSELT_PREFACTOR
    reynolds_exponent: float = COLETTO_B7_REYNOLDS_EXPONENT
    prandtl_exponent: float = COLETTO_B7_PRANDTL_EXPONENT
    chilton_colburn_analogy_exponent: float = CHILTON_COLBURN_ANALOGY_EXPONENT
    primary_source_citation_status: str = COLETTO_B7_PRIMARY_SOURCE_CITATION_STATUS

    selects_new_number: ClassVar[bool] = False
    primary_source_citation_open: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if (
            self.decision_id != LAW2_FILM_CORRELATION_DECISION_ID
            or self.authority_id != LAW2_FILM_CORRELATION_AUTHORITY_ID
            or _binary64_bits(self.nusselt_prefactor) != _binary64_bits(0.6941)
            or _binary64_bits(self.reynolds_exponent) != _binary64_bits(0.579)
            or _binary64_bits(self.prandtl_exponent) != _binary64_bits(1.0 / 3.0)
            or _binary64_bits(self.chilton_colburn_analogy_exponent) != _binary64_bits(2.0 / 3.0)
            or self.primary_source_citation_status != COLETTO_B7_PRIMARY_SOURCE_CITATION_STATUS
        ):
            raise EngineeringConfigurationError(
                "the Law-2 film-correlation authority declaration drifted from "
                "the bound Coletto B.7 / Chilton-Colburn constants"
            )


FSG_DHW_DIFFUSIVITY_AUTHORITY = FSGDhwDiffusivityAuthority()
LAW2_FILM_CORRELATION_AUTHORITY = Law2FilmCorrelationAuthority()


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringGasBoundary:
    """Serialized inlet topology and downstream pressure boundary."""

    inlet_molar_flow_mol_s: float
    inlet_hexane_mole_fraction: float
    inlet_water_mole_fraction: float
    inlet_temperature_k: float
    downstream_boundary_pressure_pa: float
    carrier_topology: CarrierTopology = CarrierTopology.POSITIVE_CARRIER

    def __post_init__(self) -> None:
        _require_positive("inlet_molar_flow_mol_s", self.inlet_molar_flow_mol_s)
        _require_closed_fraction("inlet_hexane_mole_fraction", self.inlet_hexane_mole_fraction)
        _require_closed_fraction("inlet_water_mole_fraction", self.inlet_water_mole_fraction)
        _require_positive("inlet_temperature_k", self.inlet_temperature_k)
        _require_positive("downstream_boundary_pressure_pa", self.downstream_boundary_pressure_pa)
        if type(self.carrier_topology) is not CarrierTopology:
            raise CarrierTopologyConfigurationError(
                "carrier_topology must be an exact CarrierTopology member"
            )
        volatile_sum = self.inlet_hexane_mole_fraction + self.inlet_water_mole_fraction
        if self.carrier_topology is CarrierTopology.BINARY_NO_INERT:
            if volatile_sum != 1.0:
                raise CarrierTopologyConfigurationError(
                    "BINARY_NO_INERT requires y_h + y_w == 1 exactly; refused without normalisation"
                )
        elif volatile_sum >= 1.0:
            raise CarrierTopologyConfigurationError(
                "POSITIVE_CARRIER requires y_h + y_w < 1 exactly"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class Law2TransportInputs:
    """Declared Coletto/FSG primitives for the bounded Law-2 path."""

    gas_dynamic_viscosity_pa_s: float
    gas_thermal_conductivity_w_m_k: float
    gas_mass_heat_capacity_j_kg_k: float
    hydraulic_particle_diameter_m: float
    fsg_diffusivity_multiplier: float = 1.0
    coefficient_model: Law2CoefficientModel = Law2CoefficientModel.COLETTO_NOMINAL

    fsg_hexane_molar_mass_kg_mol: ClassVar[float] = FSG_HEXANE_MOLAR_MASS_KG_MOL
    canonical_hexane_molar_mass_kg_mol: ClassVar[float] = CANONICAL_HEXANE_MOLAR_MASS_KG_MOL
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_positive("gas_dynamic_viscosity_pa_s", self.gas_dynamic_viscosity_pa_s)
        _require_positive("gas_thermal_conductivity_w_m_k", self.gas_thermal_conductivity_w_m_k)
        _require_positive("gas_mass_heat_capacity_j_kg_k", self.gas_mass_heat_capacity_j_kg_k)
        _require_positive("hydraulic_particle_diameter_m", self.hydraulic_particle_diameter_m)
        _require_positive("fsg_diffusivity_multiplier", self.fsg_diffusivity_multiplier)
        low, high = LAW2_DIFFUSIVITY_MULTIPLIER_DOMAIN
        if not low <= self.fsg_diffusivity_multiplier <= high:
            raise EngineeringConfigurationError(
                f"fsg_diffusivity_multiplier lies outside [{low}, {high}]"
            )
        if self.coefficient_model is not Law2CoefficientModel.COLETTO_NOMINAL:
            raise EngineeringConfigurationError(
                "only the declared Coletto nominal model is executable in this tranche"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class Law1TransportInputs:
    """Separately declared positive-carrier conductance ``kappa = c*k``.

    This authority is intentionally not a Law-2 coefficient model.  In
    particular, it is never substituted for the FSG/Coletto ``c*k_hw`` used by
    the exact-zero carrier branch.
    """

    positive_carrier_molar_conductance_mol_m2_s: float
    authority: Law1ConductanceAuthority = Law1ConductanceAuthority.DECLARED_POSITIVE_CARRIER

    distinct_from_law2_coefficient_authority: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_positive(
            "positive_carrier_molar_conductance_mol_m2_s",
            self.positive_carrier_molar_conductance_mol_m2_s,
        )
        if self.authority is not Law1ConductanceAuthority.DECLARED_POSITIVE_CARRIER:
            raise EngineeringConfigurationError(
                "the positive-carrier conductance must carry its exact Law-1 authority"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringFeasibilityClaimBoundary:
    engineering_feasibility_only: bool = True
    f4: bool = False
    f_gate_advanced: bool = False
    physically_qualifying: bool = False
    plant_predictive: bool = False
    production_wiring_added: bool = False
    four_branch_obligation_met: bool = False

    def __post_init__(self) -> None:
        expected = (True, False, False, False, False, False, False)
        actual = (
            self.engineering_feasibility_only,
            self.f4,
            self.f_gate_advanced,
            self.physically_qualifying,
            self.plant_predictive,
            self.production_wiring_added,
            self.four_branch_obligation_met,
        )
        if actual != expected:
            raise EngineeringConfigurationError(
                "engineering claim flags are fixed; this module cannot advance qualification"
            )


ENGINEERING_FEASIBILITY_CLAIMS = EngineeringFeasibilityClaimBoundary()


class CoreSuppliedClosureAuthority(enum.Enum):
    """Authority identity for the exhausted-film hexane interface closure."""

    DECLARED_FALLING_RATE_SORPTION = "declared_falling_rate_sorption"


@dataclass(frozen=True, slots=True, kw_only=True)
class CoreSuppliedHexaneClosure:
    """B1 stage-3 CORE_SUPPLIED arm (owner-signed design packet D-B1-2/3).

    Interface closure for a layer whose attached hexane film is exactly
    exhausted while the water film persists (frozen PHY-007 EXTERNAL_WATER
    active set): the hexane liquid reference in the interface isofugacity row
    is scaled by the sorbed-phase activity, and the hexane-side film
    conductance is scaled by the receding-front series-resistance factor from
    the falling-rate law.  Both scalars are CALLER-BOUND from the sorption
    authority and ``qsc_layer_falling_rate_law`` - this record declares them,
    it never re-derives them.
    """

    #: Sorbed-phase hexane activity ``a_h`` in (0, 1]; multiplies the
    #: pure-liquid reference fugacity in the interface equilibrium row.
    hexane_activity: float
    #: Receding-front series factor ``1/(1 + Bi_m (1-f)/f)`` in (0, 1];
    #: multiplies the gas-film hexane conductance in the diffusive flux term.
    hexane_conductance_factor: float
    #: D9-b (owner ruling 2026-09-13): the front's dimensionless radius
    #: ``f = s/R`` in [0, 1], CALLER-BOUND from the kernel's own value at
    #: ``sp1_k2_law2_tray_integration.py:1248`` - this record declares it and
    #: never re-derives it.  ``1.0`` is the front at the surface (no shell),
    #: which is what every pre-D9 caller means, so it is the default and the
    #: shipped constructions stay valid unchanged.
    front_fraction: float = 1.0
    authority: CoreSuppliedClosureAuthority = (
        CoreSuppliedClosureAuthority.DECLARED_FALLING_RATE_SORPTION
    )

    def __post_init__(self) -> None:
        for name, value in (
            ("hexane_activity", self.hexane_activity),
            ("hexane_conductance_factor", self.hexane_conductance_factor),
        ):
            _require_finite(name, value)
            if not 0.0 < value <= 1.0:
                raise EngineeringConfigurationError(
                    f"{name} must lie in (0, 1]; the CORE_SUPPLIED closure can "
                    "only hinder the film law, never amplify it"
                )
        _require_finite("front_fraction", self.front_fraction)
        if not 0.0 <= self.front_fraction <= 1.0:
            raise EngineeringConfigurationError(
                "front_fraction must lie in [0, 1]; outside it there is no "
                "shell geometry to report - refuse, do not clamp"
            )
        if self.authority is not CoreSuppliedClosureAuthority.DECLARED_FALLING_RATE_SORPTION:
            raise EngineeringConfigurationError(
                "the CORE_SUPPLIED closure must carry its exact declared authority"
            )


class CoreSuppliedWaterClosureAuthority(enum.Enum):
    """Authority identity for the exhausted-water-film interface closure."""

    DECLARED_LUIKOV_SORPTION = "declared_luikov_sorption"


@dataclass(frozen=True, slots=True, kw_only=True)
class CoreSuppliedWaterClosure:
    """B1 stage-4 CORE_SUPPLIED water arm (signed plan, water symmetry).

    Interface closure for a layer whose external water film is exactly
    exhausted while the hexane film persists (frozen PHY-007 EXTERNAL_HEXANE
    active set): the water liquid reference in the interface isofugacity row
    is scaled by the retained-water sorption activity (frozen PHY-031
    modified-Luikov, caller-bound via ``tray.luikov_activity``), and the
    water-side film conductance is scaled by a declared series factor.  Both
    scalars are CALLER-BOUND - this record declares them, never re-derives.
    """

    #: Retained-water activity ``a_w`` in (0, 1]; multiplies the pure-liquid
    #: reference fugacity in the interface water equilibrium row (PHY-031).
    water_activity: float
    #: Declared series factor in (0, 1] on the water diffusive film term
    #: (1.0 = activity-only hindrance, the first-tranche declared choice).
    water_conductance_factor: float
    authority: CoreSuppliedWaterClosureAuthority = (
        CoreSuppliedWaterClosureAuthority.DECLARED_LUIKOV_SORPTION
    )

    def __post_init__(self) -> None:
        for name, value in (
            ("water_activity", self.water_activity),
            ("water_conductance_factor", self.water_conductance_factor),
        ):
            _require_finite(name, value)
            if not 0.0 < value <= 1.0:
                raise EngineeringConfigurationError(
                    f"{name} must lie in (0, 1]; the CORE_SUPPLIED closure can "
                    "only hinder the film law, never amplify it"
                )
        if self.authority is not CoreSuppliedWaterClosureAuthority.DECLARED_LUIKOV_SORPTION:
            raise EngineeringConfigurationError(
                "the CORE_SUPPLIED water closure must carry its exact declared authority"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringCellInputs:
    """One bounded zero-inert engineering cell, with no stateful film ownership."""

    grid: cc.GridAuthority
    geometry: cc.CellGeometry
    gas_boundary: EngineeringGasBoundary
    solid: cc.SolidSideAggregates
    hydraulics: cc.LayerHydraulics
    film_area_law: cc.FilmAreaLaw
    transfer: cc.TransferCoefficients
    properties: cc.DeclaredGasProperties | nc.NativeGasProperties
    transport: Law2TransportInputs
    wall: cc.WallNodeParameters
    wall_temperature_k: float
    macro_step_s: float
    gas_closure: cc.GasClosure = cc.GasClosure.QUASI_STEADY
    film_closure: cc.FilmClosure = cc.FilmClosure.QUASI_STEADY_ALGEBRAIC
    accumulation: cc.AccumulationMode = cc.AccumulationMode.TRANSIENT
    require_corroborated_pressure_drop: bool = True
    #: B1 stage 3: ``None`` keeps every sealed two-film behavior byte-identical;
    #: a declared closure admits the EXTERNAL_WATER branch with the attached
    #: hexane film exactly exhausted.
    core_supplied_hexane: CoreSuppliedHexaneClosure | None = None
    #: B1 stage 4 (water symmetry): ``None`` keeps every sealed behavior
    #: byte-identical; a declared closure admits the EXTERNAL_HEXANE branch
    #: with the external water film exactly exhausted.
    core_supplied_water: CoreSuppliedWaterClosure | None = None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for name, expected in (
            ("grid", cc.GridAuthority),
            ("geometry", cc.CellGeometry),
            ("gas_boundary", EngineeringGasBoundary),
            ("solid", cc.SolidSideAggregates),
            ("hydraulics", cc.LayerHydraulics),
            ("film_area_law", cc.FilmAreaLaw),
            ("transfer", cc.TransferCoefficients),
            ("transport", Law2TransportInputs),
            ("wall", cc.WallNodeParameters),
        ):
            if type(getattr(self, name)) is not expected:
                raise EngineeringConfigurationError(f"{name} must be an exact {expected.__name__}")
        # CELL-02b closed caloric-mode union: the declared constant-property
        # closure and the native frozen-law closure, nothing else.
        if type(self.properties) is not cc.DeclaredGasProperties and (
            type(self.properties) is not nc.NativeGasProperties
        ):
            raise EngineeringConfigurationError(
                "properties must be an exact DeclaredGasProperties or NativeGasProperties"
            )
        # CELL-02d W2 forgery gate: the native datum identity may only ride on
        # the native closure type, never on declared constants.
        if type(self.properties) is cc.DeclaredGasProperties and (
            self.properties.energy_datum_id == nc.NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID
        ):
            raise EngineeringConfigurationError(
                "the native zero-gauge datum identity cannot be carried by a "
                "declared constant-property closure"
            )
        _require_positive("wall_temperature_k", self.wall_temperature_k)
        _require_positive("macro_step_s", self.macro_step_s)
        if self.gas_boundary.carrier_topology is not CarrierTopology.BINARY_NO_INERT:
            raise PositiveCarrierTopologyRefusal(
                "the engineering Law-2 entry accepts only declared BINARY_NO_INERT; "
                "positive-carrier Law 1 remains a separate typed-refusal path"
            )
        if self.core_supplied_hexane is not None and self.core_supplied_water is not None:
            # C8 Tier 2a (owner ruling 2026-09-05 item 3): the two arms do NOT
            # assert contradictory active sets when they are declared together
            # - together they ARE the fourth frozen PHY-007 active set, the
            # DRY_SURFACE branch, "dry particle surface with retained-phase
            # constitutive laws".  Both external liquids exactly exhausted,
            # both held by the solid, each presented to the gas through its own
            # sorbed-phase activity.  The mutual exclusion that stood here
            # until 2026-09-05 is lifted; every OTHER refusal below is verbatim.
            if type(self.core_supplied_hexane) is not CoreSuppliedHexaneClosure:
                raise EngineeringConfigurationError(
                    "core_supplied_hexane must be an exact CoreSuppliedHexaneClosure"
                )
            if type(self.core_supplied_water) is not CoreSuppliedWaterClosure:
                raise EngineeringConfigurationError(
                    "core_supplied_water must be an exact CoreSuppliedWaterClosure"
                )
            if self.solid.interface_branch is not cc.InterfaceBranch.DRY_SURFACE:
                raise UnsupportedInventoryBranchRefusal(
                    "the doubly-sorbed CORE_SUPPLIED pair is executable only on the "
                    "DRY_SURFACE branch: attached hexane AND external water both "
                    "exactly exhausted, with both liquids sorbed"
                )
        elif self.core_supplied_hexane is not None:
            if type(self.core_supplied_hexane) is not CoreSuppliedHexaneClosure:
                raise EngineeringConfigurationError(
                    "core_supplied_hexane must be an exact CoreSuppliedHexaneClosure"
                )
            if self.solid.interface_branch is not cc.InterfaceBranch.EXTERNAL_WATER:
                raise UnsupportedInventoryBranchRefusal(
                    "the CORE_SUPPLIED hexane closure is executable only on the "
                    "EXTERNAL_WATER branch: attached hexane exactly exhausted with "
                    "the water film still present"
                )
        elif self.core_supplied_water is not None:
            if type(self.core_supplied_water) is not CoreSuppliedWaterClosure:
                raise EngineeringConfigurationError(
                    "core_supplied_water must be an exact CoreSuppliedWaterClosure"
                )
            if self.solid.interface_branch is not cc.InterfaceBranch.EXTERNAL_HEXANE:
                raise UnsupportedInventoryBranchRefusal(
                    "the CORE_SUPPLIED water closure is executable only on the "
                    "EXTERNAL_HEXANE branch: external water exactly exhausted with "
                    "the hexane film still present"
                )
        else:
            if self.solid.interface_branch is not cc.InterfaceBranch.TWO_EXTERNAL_LIQUIDS:
                raise UnsupportedInventoryBranchRefusal(
                    "Law 2 is executable only for TWO_EXTERNAL_LIQUIDS; the other three "
                    "inventory branches remain permanently refused and F4 remains false"
                )
        if self.gas_closure is not cc.GasClosure.QUASI_STEADY:
            raise UnsupportedStatefulClosureRefusal(
                "dynamic gas storage is outside the bounded Law-2 engineering tranche"
            )
        if self.film_closure is not cc.FilmClosure.QUASI_STEADY_ALGEBRAIC:
            raise UnsupportedStatefulClosureRefusal(
                "transient film ownership remains with the external host"
            )
        if type(self.accumulation) is not cc.AccumulationMode:
            raise EngineeringConfigurationError(
                "accumulation must be an exact AccumulationMode member"
            )
        if type(self.require_corroborated_pressure_drop) is not bool:
            raise EngineeringConfigurationError(
                "require_corroborated_pressure_drop must be an exact bool"
            )
        if self.hydraulics.series.family is not self.geometry.tray_type.family:
            raise EngineeringConfigurationError(
                "the hydraulic series and typed tray must have the same floor family"
            )
        bed_void = self.hydraulics.bed_void_fraction
        if bed_void is None or bed_void != self.film_area_law.bed_void_fraction:
            raise EngineeringConfigurationError(
                "the Coletto path requires one shared, explicit bed void fraction"
            )
        bed_element = self.hydraulics.series.bed_element
        if bed_element is None or (
            bed_element.hydraulic_granule_diameter_m != self.transport.hydraulic_particle_diameter_m
        ):
            raise EngineeringConfigurationError(
                "the Coletto and Ergun paths must use the same hydraulic particle diameter"
            )
        if self.properties.hexane_molar_mass_kg_mol != CANONICAL_HEXANE_MOLAR_MASS_KG_MOL:
            raise EngineeringConfigurationError(
                "physical-state hexane mass must be the canonical repository value; "
                "0.08618 is reserved for the FSG correlation input"
            )
        if self.properties.water_molar_mass_kg_mol != CANONICAL_WATER_MOLAR_MASS_KG_MOL:
            raise EngineeringConfigurationError(
                "physical-state water mass must be the canonical repository value"
            )
        hexane_saturation, water_saturation = declared_external_saturations(self)
        if self.core_supplied_hexane is not None and self.core_supplied_water is not None:
            # C8 Tier 2a: BOTH external saturations exactly zero.  The film-area
            # row then reads area = A_ref * (1 - 0)**n = A_ref exactly, which is
            # the physical statement that a dry surface presents the whole bed
            # area to the gas.  Exact zeros, never a tolerance.
            if hexane_saturation != 0.0 or water_saturation != 0.0:
                raise EngineeringConfigurationError(
                    "the doubly-sorbed CORE_SUPPLIED pair requires BOTH external "
                    "saturations to be exactly zero; a positive film on either "
                    "side must run its own arm or the two-film law"
                )
        elif self.core_supplied_hexane is not None:
            if hexane_saturation != 0.0:
                raise EngineeringConfigurationError(
                    "the CORE_SUPPLIED arm requires the attached-hexane saturation to "
                    "be exactly zero; a positive film must run the two-film law"
                )
            if water_saturation <= 0.0 or water_saturation >= 1.0:
                raise EngineeringConfigurationError(
                    "the CORE_SUPPLIED arm requires one positive water-film saturation "
                    "strictly below one"
                )
        elif self.core_supplied_water is not None:
            if water_saturation != 0.0:
                raise EngineeringConfigurationError(
                    "the CORE_SUPPLIED water arm requires the external-water saturation "
                    "to be exactly zero; a positive film must run the two-film law"
                )
            if hexane_saturation <= 0.0 or hexane_saturation >= 1.0:
                raise EngineeringConfigurationError(
                    "the CORE_SUPPLIED water arm requires one positive hexane-film "
                    "saturation strictly below one"
                )
        else:
            if (
                hexane_saturation <= 0.0
                or water_saturation <= 0.0
                or hexane_saturation + water_saturation >= 1.0
            ):
                raise EngineeringConfigurationError(
                    "the two-free-liquid branch requires two positive component saturations "
                    "whose exact sum is below one"
                )

    @property
    def wall_storage_scale(self) -> float:
        if self.accumulation is cc.AccumulationMode.ZERO_ACCUMULATION:
            return 0.0
        return 1.0 if self.wall.closure is cc.WallClosure.DYNAMIC_NODE else 0.0

    @property
    def inputs_digest(self) -> str:
        payload = (
            self.grid,
            self.geometry,
            self.gas_boundary,
            self.solid,
            self.hydraulics,
            self.film_area_law,
            self.transfer,
            self.properties,
            self.transport,
            self.wall,
            self.wall_temperature_k,
            self.macro_step_s,
            self.gas_closure,
            self.film_closure,
            self.accumulation,
            self.require_corroborated_pressure_drop,
        )
        if self.core_supplied_hexane is not None:
            # Appended only when declared: every sealed two-film digest stays
            # byte-identical, and the 17-element tuple cannot collide with the
            # 16-element form.
            payload = (*payload, self.core_supplied_hexane)
        if self.core_supplied_water is not None:
            # Same conditional pattern; the two closure classes have distinct
            # reprs, so the two arms can never share a digest.
            payload = (*payload, self.core_supplied_water)
        return hashlib.sha256(repr(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierEngineeringCellInputs:
    """One bounded positive-carrier engineering cell.

    The type boundary is deliberate: a Law-1 conductance cannot be smuggled
    into the Law-2 coefficient slot, and topology dispatch happens before a
    residual or Newton iterate exists.
    """

    grid: cc.GridAuthority
    geometry: cc.CellGeometry
    gas_boundary: EngineeringGasBoundary
    solid: cc.SolidSideAggregates
    hydraulics: cc.LayerHydraulics
    film_area_law: cc.FilmAreaLaw
    transfer: cc.TransferCoefficients
    properties: cc.DeclaredGasProperties
    equilibrium: cc.InterfaceEquilibrium
    transport: Law1TransportInputs
    wall: cc.WallNodeParameters
    wall_temperature_k: float
    macro_step_s: float
    gas_closure: cc.GasClosure = cc.GasClosure.QUASI_STEADY
    film_closure: cc.FilmClosure = cc.FilmClosure.QUASI_STEADY_ALGEBRAIC
    accumulation: cc.AccumulationMode = cc.AccumulationMode.TRANSIENT
    require_corroborated_pressure_drop: bool = True

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    def __post_init__(self) -> None:
        # CELL-02b scope boundary: the native caloric closure exists for the
        # water/hexane binary only; the Law-1 inert species has no frozen
        # property law, so the native mode is refused rather than extended.
        if type(self.properties) is nc.NativeGasProperties:
            raise PositiveCarrierTopologyRefusal(
                "the native caloric closure is defined only for the water/hexane "
                "binary; the positive-carrier inert species has no frozen property law"
            )
        for name, expected in (
            ("grid", cc.GridAuthority),
            ("geometry", cc.CellGeometry),
            ("gas_boundary", EngineeringGasBoundary),
            ("solid", cc.SolidSideAggregates),
            ("hydraulics", cc.LayerHydraulics),
            ("film_area_law", cc.FilmAreaLaw),
            ("transfer", cc.TransferCoefficients),
            ("properties", cc.DeclaredGasProperties),
            ("equilibrium", cc.InterfaceEquilibrium),
            ("transport", Law1TransportInputs),
            ("wall", cc.WallNodeParameters),
        ):
            if type(getattr(self, name)) is not expected:
                raise EngineeringConfigurationError(f"{name} must be an exact {expected.__name__}")
        _require_positive("wall_temperature_k", self.wall_temperature_k)
        _require_positive("macro_step_s", self.macro_step_s)
        if self.gas_boundary.carrier_topology is not CarrierTopology.POSITIVE_CARRIER:
            raise CarrierTopologyConfigurationError(
                "the positive-carrier entry accepts only declared POSITIVE_CARRIER"
            )
        if self.solid.interface_branch is not cc.InterfaceBranch.TWO_EXTERNAL_LIQUIDS:
            raise UnsupportedInventoryBranchRefusal(
                "the bounded Law-1 entry is executable only for TWO_EXTERNAL_LIQUIDS; "
                "inactive and sorbed-phase laws are not invented, and the four-branch "
                "obligation remains unmet"
            )
        if self.gas_closure is not cc.GasClosure.QUASI_STEADY:
            raise UnsupportedStatefulClosureRefusal(
                "dynamic gas storage is outside the bounded Law-1 engineering tranche"
            )
        if self.film_closure is not cc.FilmClosure.QUASI_STEADY_ALGEBRAIC:
            raise UnsupportedStatefulClosureRefusal(
                "transient film ownership remains with the external host"
            )
        if type(self.accumulation) is not cc.AccumulationMode:
            raise EngineeringConfigurationError(
                "accumulation must be an exact AccumulationMode member"
            )
        if type(self.require_corroborated_pressure_drop) is not bool:
            raise EngineeringConfigurationError(
                "require_corroborated_pressure_drop must be an exact bool"
            )
        if self.hydraulics.series.family is not self.geometry.tray_type.family:
            raise EngineeringConfigurationError(
                "the hydraulic series and typed tray must have the same floor family"
            )
        bed_void = self.hydraulics.bed_void_fraction
        if bed_void is None or bed_void != self.film_area_law.bed_void_fraction:
            raise EngineeringConfigurationError(
                "the Law-1 hydraulic and film paths require one shared bed void fraction"
            )
        if self.properties.hexane_molar_mass_kg_mol != CANONICAL_HEXANE_MOLAR_MASS_KG_MOL:
            raise EngineeringConfigurationError(
                "physical-state hexane mass must be the canonical repository value"
            )
        if self.properties.water_molar_mass_kg_mol != CANONICAL_WATER_MOLAR_MASS_KG_MOL:
            raise EngineeringConfigurationError(
                "physical-state water mass must be the canonical repository value"
            )
        hexane_saturation, water_saturation = declared_external_saturations(self)
        if (
            hexane_saturation < 0.0
            or water_saturation < 0.0
            or hexane_saturation + water_saturation >= 1.0
        ):
            raise EngineeringConfigurationError(
                "the positive-carrier component saturations must be nonnegative and sum below one"
            )

    @property
    def wall_storage_scale(self) -> float:
        if self.accumulation is cc.AccumulationMode.ZERO_ACCUMULATION:
            return 0.0
        return 1.0 if self.wall.closure is cc.WallClosure.DYNAMIC_NODE else 0.0

    @property
    def inputs_digest(self) -> str:
        payload = (
            self.grid,
            self.geometry,
            self.gas_boundary,
            self.solid,
            self.hydraulics,
            self.film_area_law,
            self.transfer,
            self.properties,
            self.equilibrium,
            self.transport,
            self.wall,
            self.wall_temperature_k,
            self.macro_step_s,
            self.gas_closure,
            self.film_closure,
            self.accumulation,
            self.require_corroborated_pressure_drop,
        )
        return hashlib.sha256(repr(payload).encode("utf-8")).hexdigest()


def declared_external_saturations(
    inputs: EngineeringCellInputs | PositiveCarrierEngineeringCellInputs,
) -> tuple[float, float]:
    volume = inputs.geometry.bed_volume_m3
    dry_holdup = inputs.solid.dry_matter_holdup_kg
    hexane = (
        inputs.solid.attached_hexane_loading_kg_kg
        * dry_holdup
        / inputs.film_area_law.hexane_inventory_capacity_kg(volume)
    )
    water = (
        inputs.solid.external_water_loading_kg_kg
        * dry_holdup
        / inputs.film_area_law.water_inventory_capacity_kg(volume)
    )
    return hexane, water


def doubly_sorbed_total_molar_flux(
    *,
    hexane_conductance_factor: float,
    water_conductance_factor: float,
    interface_hexane_mole_fraction: float,
    hexane_diffusive_flux_mol_m2_s: float,
    water_diffusive_flux_mol_m2_s: float,
) -> float:
    """C8 Tier 2a closed form for the total molar flux on a dry surface.

    With both liquids sorbed, each species carries its own whole-flux
    constitutive scaling,

        n_h = f_h (y_I n_t + D_h),   n_w = f_w ((1 - y_I) n_t + D_w),

    and the binary identity ``n_t = n_h + n_w`` is no longer automatic.
    Closing it gives design section 1 item 3,

        n_t [1 - f_h y_I - f_w (1 - y_I)] = f_h D_h + f_w D_w,

    which DEFINES the total flux whenever the bracket is non-zero.  At
    ``f_h = f_w = 1`` the bracket is exactly zero and, in a binary film where
    ``D_w = -D_h``, so is the right side: the identity degenerates to 0 = 0,
    n_t stays the free unknown the energy rows determine, and the cell does
    not evaluate this function at all.  Asked at that corner it refuses
    rather than divide - the degenerate case is the two-film law, not a
    closed form.

    Q-T2A-1 (owner ruling 2026-09-06).  The cell no longer evaluates this
    function on ANY path: the doubly-sorbed branch below unit conductance is
    a typed refusal there, because releasing the water isofugacity row left
    an interface nothing holds at the two-liquid equilibrium and the Tier 2a
    diagnostic measured it refusing at every deck temperature and every seed
    (build record section 8).  The algebra is kept here, documented and
    tested, so that a later packet which wants hindered conductances on a dry
    surface starts from a written identity rather than a rediscovered one;
    it is NOT a cell path.
    """

    bracket = (
        1.0
        - hexane_conductance_factor * interface_hexane_mole_fraction
        - water_conductance_factor * (1.0 - interface_hexane_mole_fraction)
    )
    if bracket == 0.0:
        raise EngineeringConfigurationError(
            "the doubly-sorbed closed form is degenerate at a zero bracket "
            "(both conductance factors at one): the total flux is the free "
            "unknown the energy rows determine, not a closed form"
        )
    return (
        hexane_conductance_factor * hexane_diffusive_flux_mol_m2_s
        + water_conductance_factor * water_diffusive_flux_mol_m2_s
    ) / bracket


def dispatch_engineering_carrier_topology(
    inputs: EngineeringCellInputs | PositiveCarrierEngineeringCellInputs,
) -> EngineeringCarrierLaw:
    """Select an exact carrier law before solving; never execute or fall back."""

    if type(inputs) is EngineeringCellInputs:
        if inputs.gas_boundary.carrier_topology is not CarrierTopology.BINARY_NO_INERT:
            raise CarrierTopologyConfigurationError(
                "Law-2 inputs and serialized carrier topology disagree"
            )
        return EngineeringCarrierLaw.BINARY_NO_INERT_LAW2
    if type(inputs) is PositiveCarrierEngineeringCellInputs:
        if inputs.gas_boundary.carrier_topology is not CarrierTopology.POSITIVE_CARRIER:
            raise CarrierTopologyConfigurationError(
                "Law-1 inputs and serialized carrier topology disagree"
            )
        return EngineeringCarrierLaw.POSITIVE_CARRIER_LAW1
    raise EngineeringConfigurationError(
        "carrier dispatch requires an exact Law-1 or Law-2 engineering input"
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierCellFieldState:
    gas_temperature_k: float
    gas_hexane_mole_fraction: float
    gas_water_mole_fraction: float
    layer_pressure_pa: float
    interface_temperature_k: float
    external_hexane_saturation: float
    external_water_saturation: float

    def __post_init__(self) -> None:
        for name in POSITIVE_CARRIER_UNKNOWN_NAMES:
            _require_finite(name, getattr(self, name))
        _require_positive_carrier_state_domain(
            self.as_vector(), error_type=EngineeringConfigurationError
        )

    @property
    def gas_carrier_mole_fraction(self) -> float:
        return 1.0 - self.gas_hexane_mole_fraction - self.gas_water_mole_fraction

    def as_vector(self) -> tuple[float, ...]:
        return tuple(getattr(self, name) for name in POSITIVE_CARRIER_UNKNOWN_NAMES)

    @classmethod
    def from_vector(cls, vector: tuple[float, ...]) -> PositiveCarrierCellFieldState:
        if len(vector) != len(POSITIVE_CARRIER_UNKNOWN_NAMES):
            raise EngineeringConfigurationError("the positive-carrier state carries seven unknowns")
        return cls(**dict(zip(POSITIVE_CARRIER_UNKNOWN_NAMES, vector, strict=True)))


def _require_positive_carrier_state_domain(
    vector: tuple[float, ...], *, error_type: type[EngineeringFeasibilityError]
) -> None:
    if len(vector) != 7:
        raise error_type("the positive-carrier state carries seven unknowns")
    if any(type(value) is not float or not math.isfinite(value) for value in vector):
        raise error_type("a positive-carrier state value is not finite binary64")
    t_g, y_h, y_w, pressure, t_interface, s_h, s_w = vector
    t_low, t_high = LAW1_TEMPERATURE_DOMAIN_K
    if not t_low <= t_g <= t_high or not t_low <= t_interface <= t_high:
        raise error_type("a trial temperature left the declared Law-1 domain")
    p_low, p_high = LAW1_PRESSURE_DOMAIN_PA
    if not p_low <= pressure <= p_high:
        raise error_type("trial pressure left the declared Law-1 domain")
    if not 0.0 <= y_h <= 1.0 or not 0.0 <= y_w <= 1.0:
        raise error_type("a trial volatile composition left [0, 1]")
    if 1.0 - y_h - y_w <= 0.0:
        raise error_type("a trial bulk composition left no positive carrier")
    if s_h < 0.0 or s_w < 0.0 or s_h + s_w >= 1.0:
        raise error_type("trial component saturation is inadmissible; refused without projection")


@dataclass(frozen=True, slots=True, kw_only=True)
class BinaryNoInertCellFieldState:
    gas_temperature_k: float
    gas_hexane_mole_fraction: float
    layer_pressure_pa: float
    interface_temperature_k: float
    interface_hexane_mole_fraction: float
    external_hexane_saturation: float
    external_water_saturation: float
    total_molar_flux_mol_m2_s: float

    def __post_init__(self) -> None:
        for name in BINARY_NO_INERT_UNKNOWN_NAMES:
            _require_finite(name, getattr(self, name))
        _require_state_domain(self.as_vector(), error_type=EngineeringConfigurationError)

    @property
    def gas_water_mole_fraction(self) -> float:
        return 1.0 - self.gas_hexane_mole_fraction

    @property
    def interface_water_mole_fraction(self) -> float:
        return 1.0 - self.interface_hexane_mole_fraction

    def as_vector(self) -> tuple[float, ...]:
        return tuple(getattr(self, name) for name in BINARY_NO_INERT_UNKNOWN_NAMES)

    @classmethod
    def from_vector(cls, vector: tuple[float, ...]) -> BinaryNoInertCellFieldState:
        if len(vector) != len(BINARY_NO_INERT_UNKNOWN_NAMES):
            raise EngineeringConfigurationError("the binary Law-2 state carries eight unknowns")
        return cls(**dict(zip(BINARY_NO_INERT_UNKNOWN_NAMES, vector, strict=True)))


def _require_state_domain(
    vector: tuple[float, ...], *, error_type: type[EngineeringFeasibilityError]
) -> None:
    if len(vector) != 8:
        raise error_type("the binary Law-2 state carries eight unknowns")
    if any(type(value) is not float or not math.isfinite(value) for value in vector):
        raise error_type("a binary Law-2 state value is not finite binary64")
    gas_temperature, y_hexane, pressure, interface_temperature, y_interface, s_h, s_w, n_t = vector
    t_low, t_high = LAW2_TEMPERATURE_DOMAIN_K
    if not t_low <= gas_temperature <= t_high:
        raise error_type("trial gas temperature left the declared Law-2 domain")
    if not t_low <= interface_temperature <= t_high:
        raise error_type("trial interface temperature left the declared Law-2 domain")
    p_low, p_high = LAW2_PRESSURE_DOMAIN_PA
    if not p_low <= pressure <= p_high:
        raise error_type("trial pressure left the declared Law-2 domain")
    if not 0.0 <= y_hexane <= 1.0 or not 0.0 <= y_interface <= 1.0:
        raise error_type("trial intrinsic-binary composition left [0, 1]")
    if s_h < 0.0 or s_w < 0.0 or s_h + s_w >= 1.0:
        raise error_type("trial component saturation is inadmissible; refused without projection")
    n_low, n_high = LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S
    if not n_low <= n_t <= n_high:
        raise error_type("trial total molar flux left the declared Law-2 domain")


def _require_dry_shell_state_domain(
    vector: tuple[float, ...], *, error_type: type[EngineeringFeasibilityError]
) -> None:
    """D9-b-1: the ten-unknown dry-shell trial domain.

    The shipped eight-unknown domain applied to the first eight entries -
    including the Law-2 TEMPERATURE FLOOR, which STAYS.  The closure record's
    section 3b measured this interface's root at the captured crossing sitting
    at 291.6 K, BELOW the declared floor; a root below the domain therefore
    refuses typed through the ordinary convergence refusal and is never
    admitted by widening the floor.  The front temperature is held to the same
    declared temperature domain; the surface hexane activity to
    :data:`DRY_SHELL_SURFACE_ACTIVITY_DOMAIN`.
    """

    if len(vector) != len(DRY_SHELL_UNKNOWN_NAMES):
        raise error_type("the dry-shell Law-2 state carries ten unknowns")
    _require_state_domain(tuple(vector[:8]), error_type=error_type)
    front_temperature = vector[8]
    surface_activity = vector[9]
    if any(
        type(value) is not float or not math.isfinite(value)
        for value in (front_temperature, surface_activity)
    ):
        raise error_type("a dry-shell Law-2 state value is not finite binary64")
    t_low, t_high = LAW2_TEMPERATURE_DOMAIN_K
    if not t_low <= front_temperature <= t_high:
        raise error_type("trial front temperature left the declared Law-2 domain")
    activity_low, activity_high = DRY_SHELL_SURFACE_ACTIVITY_DOMAIN
    if not activity_low < surface_activity <= activity_high:
        raise error_type("trial surface hexane activity left the declared domain")


@dataclass(frozen=True, slots=True, kw_only=True)
class DryShellCellFieldState:
    """D9-b-1: the shipped binary state plus the two front unknowns."""

    gas_temperature_k: float
    gas_hexane_mole_fraction: float
    layer_pressure_pa: float
    interface_temperature_k: float
    interface_hexane_mole_fraction: float
    external_hexane_saturation: float
    external_water_saturation: float
    total_molar_flux_mol_m2_s: float
    front_temperature_k: float
    surface_hexane_activity: float

    def __post_init__(self) -> None:
        for name in DRY_SHELL_UNKNOWN_NAMES:
            _require_finite(name, getattr(self, name))
        _require_dry_shell_state_domain(self.as_vector(), error_type=EngineeringConfigurationError)

    @property
    def gas_water_mole_fraction(self) -> float:
        return 1.0 - self.gas_hexane_mole_fraction

    @property
    def interface_water_mole_fraction(self) -> float:
        return 1.0 - self.interface_hexane_mole_fraction

    @property
    def shell_temperature_drop_k(self) -> float:
        """The drop across the dry shell: reported, never solved for."""

        return self.interface_temperature_k - self.front_temperature_k

    def as_vector(self) -> tuple[float, ...]:
        return tuple(getattr(self, name) for name in DRY_SHELL_UNKNOWN_NAMES)

    @classmethod
    def from_vector(cls, vector: tuple[float, ...]) -> DryShellCellFieldState:
        if len(vector) != len(DRY_SHELL_UNKNOWN_NAMES):
            raise EngineeringConfigurationError("the dry-shell Law-2 state carries ten unknowns")
        return cls(**dict(zip(DRY_SHELL_UNKNOWN_NAMES, vector, strict=True)))


class _Dual:
    """Eight-direction forward tangent used as the analytic Jacobian authority.

    D9-b-1 widened it: the tangent width is carried by the OPERANDS rather
    than fixed at eight, so the ten-unknown dry-shell state uses the same
    algebra.  Every default is still eight and the eight-unknown path is
    unchanged bit for bit.
    """

    __slots__ = ("partials", "value")

    def __init__(self, value: float, partials: tuple[float, ...] | None = None) -> None:
        self.value = float(value)
        self.partials = (0.0,) * 8 if partials is None else partials

    @classmethod
    def variable(cls, value: float, index: int, *, width: int = 8) -> _Dual:
        """D9-b-1: ``width`` is the number of unknowns the state carries.

        It defaults to the shipped eight, so every pre-D9 call site is
        unchanged bit for bit; the ten-unknown dry-shell state passes ten and
        the whole tangent algebra follows, because every zero tangent below is
        built from the width of an operand that already carries one.
        """

        tangent = [0.0] * width
        tangent[index] = 1.0
        return cls(value, tuple(tangent))

    def _coerce(self, other: float | _Dual) -> _Dual:
        """``_as_dual`` with the tangent WIDTH of ``self``.

        At width eight this returns exactly what ``_as_dual`` returns - the
        same object for a dual, an identical zero-tangent dual for a float -
        so the eight-unknown path is unchanged; at width ten it keeps the
        strict zips below from ever seeing two different widths.  The same
        pattern ``_Law1Dual._coerce`` already carries.
        """

        if type(other) is _Dual:
            return other
        return _Dual(float(other), (0.0,) * len(self.partials))

    def __add__(self, other: float | _Dual) -> _Dual:
        rhs = self._coerce(other)
        return _Dual(
            self.value + rhs.value,
            tuple(a + b for a, b in zip(self.partials, rhs.partials, strict=True)),
        )

    __radd__ = __add__

    def __neg__(self) -> _Dual:
        return _Dual(-self.value, tuple(-value for value in self.partials))

    def __sub__(self, other: float | _Dual) -> _Dual:
        return self + (-self._coerce(other))

    def __rsub__(self, other: float | _Dual) -> _Dual:
        return self._coerce(other) - self

    def __mul__(self, other: float | _Dual) -> _Dual:
        rhs = self._coerce(other)
        return _Dual(
            self.value * rhs.value,
            tuple(
                left * rhs.value + self.value * right
                for left, right in zip(self.partials, rhs.partials, strict=True)
            ),
        )

    __rmul__ = __mul__

    def __truediv__(self, other: float | _Dual) -> _Dual:
        rhs = self._coerce(other)
        if rhs.value == 0.0:
            raise EngineeringStepRejected("an analytic division encountered zero")
        denominator = rhs.value * rhs.value
        return _Dual(
            self.value / rhs.value,
            tuple(
                (left * rhs.value - self.value * right) / denominator
                for left, right in zip(self.partials, rhs.partials, strict=True)
            ),
        )

    def __rtruediv__(self, other: float | _Dual) -> _Dual:
        return self._coerce(other) / self

    def __pow__(self, exponent: float) -> _Dual:
        if self.value <= 0.0:
            raise EngineeringStepRejected("a fractional-power base is not positive")
        value = self.value**exponent
        factor = exponent * self.value ** (exponent - 1.0)
        return _Dual(value, tuple(factor * partial for partial in self.partials))


def _as_dual(value: float | _Dual) -> _Dual:
    return value if type(value) is _Dual else _Dual(float(value))


def bernoulli_pair(argument: float) -> tuple[float, float]:
    """Return B(x), B'(x), including the exact removable limit at zero."""

    _require_finite("Bernoulli argument", argument)
    if argument == 0.0:
        return 1.0, -0.5
    if abs(argument) < 1.0e-5:
        square = argument * argument
        value = 1.0 - 0.5 * argument + square / 12.0 - square * square / 720.0
        derivative = -0.5 + argument / 6.0 - argument * square / 180.0
        return value, derivative
    denominator = math.expm1(argument)
    value = argument / denominator
    derivative = (denominator - argument * math.exp(argument)) / (denominator * denominator)
    return value, derivative


def _bernoulli(argument: _Dual) -> _Dual:
    value, derivative = bernoulli_pair(argument.value)
    return _Dual(value, tuple(derivative * partial for partial in argument.partials))


class _Law1Dual:
    """Seven-direction forward tangent used by the Law-1 Jacobian authority."""

    __slots__ = ("partials", "value")

    def __init__(self, value: float, partials: tuple[float, ...] | None = None) -> None:
        self.value = float(value)
        self.partials = (0.0,) * 7 if partials is None else partials

    @classmethod
    def variable(cls, value: float, index: int) -> _Law1Dual:
        tangent = [0.0] * 7
        tangent[index] = 1.0
        return cls(value, tuple(tangent))

    def _coerce(self, other: float | _Law1Dual) -> _Law1Dual:
        return other if type(other) is _Law1Dual else _Law1Dual(float(other))

    def __add__(self, other: float | _Law1Dual) -> _Law1Dual:
        rhs = self._coerce(other)
        return _Law1Dual(
            self.value + rhs.value,
            tuple(a + b for a, b in zip(self.partials, rhs.partials, strict=True)),
        )

    __radd__ = __add__

    def __neg__(self) -> _Law1Dual:
        return _Law1Dual(-self.value, tuple(-value for value in self.partials))

    def __sub__(self, other: float | _Law1Dual) -> _Law1Dual:
        return self + (-self._coerce(other))

    def __rsub__(self, other: float | _Law1Dual) -> _Law1Dual:
        return self._coerce(other) - self

    def __mul__(self, other: float | _Law1Dual) -> _Law1Dual:
        rhs = self._coerce(other)
        return _Law1Dual(
            self.value * rhs.value,
            tuple(
                left * rhs.value + self.value * right
                for left, right in zip(self.partials, rhs.partials, strict=True)
            ),
        )

    __rmul__ = __mul__

    def __truediv__(self, other: float | _Law1Dual) -> _Law1Dual:
        rhs = self._coerce(other)
        if rhs.value == 0.0:
            raise EngineeringStepRejected("a Law-1 analytic division encountered zero")
        denominator = rhs.value * rhs.value
        return _Law1Dual(
            self.value / rhs.value,
            tuple(
                (left * rhs.value - self.value * right) / denominator
                for left, right in zip(self.partials, rhs.partials, strict=True)
            ),
        )

    def __rtruediv__(self, other: float | _Law1Dual) -> _Law1Dual:
        return self._coerce(other) / self

    def __pow__(self, exponent: float) -> _Law1Dual:
        if self.value <= 0.0:
            raise EngineeringStepRejected("a Law-1 fractional-power base is not positive")
        value = self.value**exponent
        factor = exponent * self.value ** (exponent - 1.0)
        return _Law1Dual(value, tuple(factor * partial for partial in self.partials))


def _law1_exp(argument: _Law1Dual) -> _Law1Dual:
    try:
        value = math.exp(argument.value)
    except OverflowError as error:
        raise EngineeringStepRejected("a Law-1 exponential is not representable") from error
    if not math.isfinite(value):
        raise EngineeringStepRejected("a Law-1 exponential is not representable")
    return _Law1Dual(value, tuple(value * partial for partial in argument.partials))


def _law1_log(argument: _Law1Dual) -> _Law1Dual:
    if argument.value <= 0.0:
        raise EngineeringStepRejected("a Law-1 logarithm argument is not positive")
    return _Law1Dual(
        math.log(argument.value),
        tuple(partial / argument.value for partial in argument.partials),
    )


def _law1_bernoulli(argument: _Law1Dual) -> _Law1Dual:
    value, derivative = bernoulli_pair(argument.value)
    return _Law1Dual(value, tuple(derivative * partial for partial in argument.partials))


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierFluxVector:
    bulk_carrier_mole_fraction: float
    interface_carrier_mole_fraction: float
    peclet: float
    bernoulli_factor: float
    total_molar_flux_mol_m2_s: float
    hexane_molar_flux_mol_m2_s: float
    water_molar_flux_mol_m2_s: float
    carrier_molar_flux_mol_m2_s: float
    hexane_diffusive_flux_mol_m2_s: float
    water_diffusive_flux_mol_m2_s: float
    carrier_diffusive_flux_mol_m2_s: float
    component_sum_defect_mol_m2_s: float

    laboratory_frame_components: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def positive_carrier_fluxes(
    *,
    bulk_hexane_mole_fraction: float,
    bulk_water_mole_fraction: float,
    interface_hexane_mole_fraction: float,
    interface_water_mole_fraction: float,
    positive_carrier_molar_conductance_mol_m2_s: float,
) -> PositiveCarrierFluxVector:
    """Evaluate the corrected state-driven Law-1 laboratory-frame flux vector."""

    for name, value in (
        ("bulk_hexane_mole_fraction", bulk_hexane_mole_fraction),
        ("bulk_water_mole_fraction", bulk_water_mole_fraction),
        ("interface_hexane_mole_fraction", interface_hexane_mole_fraction),
        ("interface_water_mole_fraction", interface_water_mole_fraction),
    ):
        _require_closed_fraction(name, value)
    _require_positive(
        "positive_carrier_molar_conductance_mol_m2_s",
        positive_carrier_molar_conductance_mol_m2_s,
    )
    g_bulk = 1.0 - bulk_hexane_mole_fraction - bulk_water_mole_fraction
    if g_bulk <= 0.0:
        raise CarrierTopologyConfigurationError(
            "positive-carrier flux evaluation requires g_bulk > 0 exactly"
        )
    g_interface = 1.0 - interface_hexane_mole_fraction - interface_water_mole_fraction
    if g_interface <= 0.0:
        raise PositiveCarrierInterfaceCrossoverRefusal(
            "positive-carrier flux evaluation reached g_interface <= 0; refused without crossover"
        )
    peclet = math.log(g_bulk) - math.log(g_interface)
    bernoulli, _ = bernoulli_pair(peclet)
    kappa = positive_carrier_molar_conductance_mol_m2_s
    n_total = kappa * peclet
    j_hexane = kappa * bernoulli * (interface_hexane_mole_fraction - bulk_hexane_mole_fraction)
    j_water = kappa * bernoulli * (interface_water_mole_fraction - bulk_water_mole_fraction)
    j_carrier = kappa * bernoulli * (g_interface - g_bulk)
    n_hexane = interface_hexane_mole_fraction * n_total + j_hexane
    n_water = interface_water_mole_fraction * n_total + j_water
    n_carrier = g_interface * n_total + j_carrier
    return PositiveCarrierFluxVector(
        bulk_carrier_mole_fraction=g_bulk,
        interface_carrier_mole_fraction=g_interface,
        peclet=peclet,
        bernoulli_factor=bernoulli,
        total_molar_flux_mol_m2_s=n_total,
        hexane_molar_flux_mol_m2_s=n_hexane,
        water_molar_flux_mol_m2_s=n_water,
        carrier_molar_flux_mol_m2_s=n_carrier,
        hexane_diffusive_flux_mol_m2_s=j_hexane,
        water_diffusive_flux_mol_m2_s=j_water,
        carrier_diffusive_flux_mol_m2_s=j_carrier,
        component_sum_defect_mol_m2_s=(n_hexane + n_water + n_carrier - n_total),
    )


def _lift_partials(value: float, dependencies: tuple[tuple[_Dual, float], ...]) -> _Dual:
    """D9-b-1: the tangent width is the width of the dependencies that carry it.

    With at least one dependency - which every shipped call site has - the lift
    returns a tangent of exactly that width, so the same helper serves the
    eight-unknown and the ten-unknown states.  With none there is nothing to
    lift and the shipped eight-wide zero tangent is returned unchanged.
    """

    width = len(dependencies[0][0].partials) if dependencies else 8
    partials = [0.0] * width
    for dependency, derivative in dependencies:
        for index, tangent in enumerate(dependency.partials):
            partials[index] += derivative * tangent
    return _Dual(value, tuple(partials))


def _vapor_fugacity_coefficient(
    temperature: _Dual, pressure: _Dual, y_hexane: _Dual, *, species: str
) -> _Dual:
    try:
        state = binary_gas.state(
            temperature.value,
            pressure.value,
            1.0 - y_hexane.value,
            y_hexane.value,
        )
    except (ValueError, OverflowError) as error:
        raise EngineeringStepRejected(f"binary-gas fugacity state refused: {error}") from error
    y_water = state.y_water
    y_h = state.y_hexane
    d_bmix_d_y = (
        -2.0 * y_water * state.B_ww + 2.0 * (y_water - y_h) * state.B_wh + 2.0 * y_h * state.B_hh
    )
    if species == "hexane":
        coefficient = 2.0 * (y_water * state.B_wh + y_h * state.B_hh) - state.B_mix
        d_coefficient_d_t = (
            2.0 * (y_water * state.dB_wh_dT + y_h * state.dB_hh_dT) - state.dB_mix_dT
        )
        d_coefficient_d_y = 2.0 * (state.B_hh - state.B_wh) - d_bmix_d_y
        phi = state.phi_hexane
    else:
        coefficient = 2.0 * (y_water * state.B_ww + y_h * state.B_wh) - state.B_mix
        d_coefficient_d_t = (
            2.0 * (y_water * state.dB_ww_dT + y_h * state.dB_wh_dT) - state.dB_mix_dT
        )
        d_coefficient_d_y = 2.0 * (state.B_wh - state.B_ww) - d_bmix_d_y
        phi = state.phi_water
    d_phi_d_t = (
        phi
        * (pressure.value / UNIVERSAL_GAS_CONSTANT_J_MOL_K)
        * (d_coefficient_d_t * temperature.value - coefficient)
        / (temperature.value * temperature.value)
    )
    d_phi_d_p = phi * coefficient / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature.value)
    d_phi_d_y = (
        phi
        * pressure.value
        / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature.value)
        * d_coefficient_d_y
    )
    return _lift_partials(
        phi,
        (
            (temperature, d_phi_d_t),
            (pressure, d_phi_d_p),
            (y_hexane, d_phi_d_y),
        ),
    )


def _liquid_reference_fugacity(module: object, temperature: _Dual, pressure: _Dual) -> _Dual:
    """Evaluate the checked Helmholtz saturation map and its IFT tangent."""

    t_value = temperature.value
    p_value = pressure.value
    try:
        saturation = module.saturation(t_value)
        if not saturation.converged:
            raise EngineeringStepRejected(f"{module.__name__}.saturation reports converged=False")
        rho_liquid = saturation.rho_liquid
        rho_vapor = saturation.rho_vapor
        gas_constant = module.R

        def d_pressure_d_temperature(rho: float) -> float:
            delta = rho / module.RHOC
            tau = module.TC / t_value
            residual = module.residual(delta, tau)
            return (
                rho
                * gas_constant
                * (1.0 + delta * residual.phir_d - delta * tau * residual.phir_dt)
            )

        def d_log_fugacity_d_temperature(rho: float) -> float:
            delta = rho / module.RHOC
            tau = module.TC / t_value
            residual = module.residual(delta, tau)
            return 1.0 / t_value - (tau / t_value) * (residual.phir_t + delta * residual.phir_dt)

        d_p_liquid = module._dp_drho(t_value, rho_liquid)
        d_p_vapor = module._dp_drho(t_value, rho_vapor)
        d_f_liquid = d_p_liquid / (rho_liquid * gas_constant * t_value)
        d_f_vapor = d_p_vapor / (rho_vapor * gas_constant * t_value)
        determinant = d_p_liquid * (-d_f_vapor) - (-d_p_vapor) * d_f_liquid
        if determinant == 0.0:
            raise EngineeringStepRejected("the saturation-map IFT Jacobian is singular")
        rhs_pressure = -(d_pressure_d_temperature(rho_liquid) - d_pressure_d_temperature(rho_vapor))
        rhs_fugacity = -(
            d_log_fugacity_d_temperature(rho_liquid) - d_log_fugacity_d_temperature(rho_vapor)
        )
        d_rho_liquid_d_t = (rhs_pressure * (-d_f_vapor) - (-d_p_vapor) * rhs_fugacity) / determinant
        d_saturation_pressure_d_t = (
            d_pressure_d_temperature(rho_liquid) + d_p_liquid * d_rho_liquid_d_t
        )

        log_fugacity = module.ln_fugacity(t_value, rho_liquid)
        poynting = (p_value - saturation.p) / (rho_liquid * gas_constant * t_value)
        fugacity = math.exp(log_fugacity + poynting)
        d_log_fugacity_sat_d_t = (
            d_log_fugacity_d_temperature(rho_liquid) + d_f_liquid * d_rho_liquid_d_t
        )
        numerator = p_value - saturation.p
        denominator = rho_liquid * gas_constant * t_value
        d_numerator_d_t = -d_saturation_pressure_d_t
        d_denominator_d_t = gas_constant * (rho_liquid + t_value * d_rho_liquid_d_t)
        d_poynting_d_t = (d_numerator_d_t * denominator - numerator * d_denominator_d_t) / (
            denominator * denominator
        )
        d_fugacity_d_t = fugacity * (d_log_fugacity_sat_d_t + d_poynting_d_t)
        d_fugacity_d_p = fugacity / denominator
    except EngineeringStepRejected:
        raise
    except (ValueError, OverflowError, ZeroDivisionError) as error:
        raise EngineeringStepRejected(f"pure-liquid fugacity map refused: {error}") from error
    return _lift_partials(
        fugacity,
        ((temperature, d_fugacity_d_t), (pressure, d_fugacity_d_p)),
    )


def _dry_shell_front_saturation_pressure(t_front: _Dual) -> _Dual:
    """Psat_hexane(T_f) with its exact Clapeyron slope, both from the shipped
    hexane property authority (``ft2a1_shell._front_saturation_pressure``)."""

    try:
        saturation = hexane_props.saturation(t_front.value)
    except (ValueError, OverflowError) as error:
        raise EngineeringStepRejected(
            f"hexane saturation refused at the front temperature: {error}"
        ) from error
    if not saturation.converged:
        raise EngineeringStepRejected("hexane saturation did not converge at the front temperature")
    d_volume = 1.0 / saturation.rho_vapor - 1.0 / saturation.rho_liquid
    d_p_d_t = saturation.dh_vap / (t_front.value * d_volume)
    return _lift_partials(saturation.p, ((t_front, d_p_d_t),))


def _dry_shell_front_latent_heat(t_front: _Dual) -> _Dual:
    """The cell's OWN hexane latent-heat authority, read at the FRONT."""

    try:
        law = nc.hexane_latent_heat(t_front.value)
    except (ValueError, OverflowError) as error:
        raise EngineeringStepRejected(
            f"native hexane latent-heat law refused at the front temperature: {error}"
        ) from error
    return _lift_partials(law.value_j_mol, ((t_front, law.d_dT_j_mol_k),))


def _require_dry_shell_is_reachable(core_supplied: CoreSuppliedHexaneClosure | None) -> None:
    """D9-b-2 and R-C5: the two typed guards that stand in front of the branch.

    The FIRST is the no-shell refusal, AS AMENDED BY OWNER RULING D9-e
    (2026-09-13, ``GT_PS2_OWNER_RULINGS_D9E_2026-09-13.md``).  It fires at a
    front fraction of exactly one - the front AT the surface - where the
    particle has no hexane-free shell at all: nothing for the water arm to
    exchange through and nothing for the hexane duty to conduct through, so
    there is no interface to formulate rather than no root to find.  The
    branch is reached only when both arms are core-supplied, which is a layer
    with no film of either species, so "front at the surface" and "no film" is
    the whole state this refusal names.  Under D9-a's routing it does not
    occur, because a layer at a closed shell keeps its film; this refusal is
    the guard that says so, and its firing on any D9-a lane is kill criterion
    K8 as D9-e re-defines it.

    WHAT D9-e WITHDREW, and why it is not silently gone: the D9-b build put
    this guard at a front of 0.99 on the reading that below a thin shell the
    surface falls to the hexane wet-bulb BELOW the Law-2 domain (291.6 K at
    the captured crossing).  That reading was an artefact of the science
    lane's own row D9-b-1, which BLOCKED the water at the dry surface; with
    the water arm active - closure (vi), the addendum's Table N - the same
    interface has a root at every front from 0.99999913 to half radius, and
    the honest lane's death at 0.9914 was the boundary, not the physics.

    The SECOND is the Phi -> 1 guard R-C5 asks for: the shell mass conductance
    is ``k_gas Phi/(1 - Phi)``, unbounded at an unhindered series factor, and a
    closure reporting a receded front with an UNHINDERED film is not a state
    this interface can carry.  Refuse it typed rather than divide.

    THE THIRD, added by the D9-e FIX (2026-09-14, declared row
    ``d9e_carried_seed_width_from_record``): no hexane closure AT ALL.
    Every caller reaches this guard behind a check that the record
    presents the branch, so ``None`` cannot arrive here today; it is
    refused rather than dereferenced so that no future caller can read a
    front fraction off nothing, which is the shape the farm lanes died in.
    """

    if core_supplied is None:
        raise UnsupportedInventoryBranchRefusal(
            "the dry-shell receding-front interface has no hexane closure to "
            "read a front fraction from: this record binds the hexane FILM law "
            "(a positive attached film), which the branch does not carry - "
            "refuse, never read a front off nothing"
        )
    if core_supplied.front_fraction == NO_SHELL_FRONT_FRACTION_CEILING:
        raise UnsupportedInventoryBranchRefusal(
            "no shell to carry the hexane duty (front fraction "
            f"{core_supplied.front_fraction!r} at the declared "
            f"{NO_SHELL_FRONT_FRACTION_CEILING}); the front sits AT the "
            "surface, so there is no hexane-free shell to exchange water "
            "through and none to conduct the hexane duty through; a film-free "
            "layer at a closed shell is outside D9-a's routing"
        )
    if core_supplied.hexane_conductance_factor == 1.0:
        raise UnsupportedInventoryBranchRefusal(
            "the dry-shell receding-front interface has no finite shell mass "
            "conductance at an unhindered series factor (Phi = 1.0 reported "
            f"with front fraction {core_supplied.front_fraction!r}); "
            "refuse, do not divide"
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BinaryNoInertCellEvaluation:
    residual: tuple[float, ...]
    jacobian: tuple[tuple[float, ...], ...]
    row_scales: tuple[float, ...]
    scaled_residual_norm: float
    gas_side_active_area_m2: float
    outlet_molar_flow_mol_s: float
    total_molar_flux_mol_m2_s: float
    hexane_molar_flux_mol_m2_s: float
    water_molar_flux_mol_m2_s: float
    gas_density_kg_m3: float
    superficial_velocity_m_s: float
    layer_pressure_drop_pa: float
    binary_diffusivity_m2_s: float
    schmidt_number: float
    mass_transfer_coefficient_m_s: float
    molar_transport_coefficient_mol_m2_s: float
    heat_transfer_coefficient_w_m2_k: float
    wall_to_gas_w: float
    wall_to_interface_w: float
    solid_to_interface_w: float
    interface_convection_w: float
    independent_water_balance_residual_mol_s: float

    unknowns: ClassVar[tuple[str, ...]] = BINARY_NO_INERT_UNKNOWN_NAMES
    rows: ClassVar[tuple[str, ...]] = BINARY_NO_INERT_ROW_NAMES
    analytic_jacobian: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class DryShellCellEvaluation(BinaryNoInertCellEvaluation):
    """D9-b-1: the same fields, naming the TEN unknowns and rows.

    No field is added: the residual, the Jacobian and the row scales are the
    same tuples, ten long instead of eight, and the two front unknowns are
    read off the state rather than off the evaluation.  The record exists so
    that a caller cannot read ten rows through an eight-row name list.
    """

    unknowns: ClassVar[tuple[str, ...]] = DRY_SHELL_UNKNOWN_NAMES
    rows: ClassVar[tuple[str, ...]] = DRY_SHELL_ROW_NAMES


def _evaluate_binary_no_inert(
    inputs: EngineeringCellInputs,
    vector: tuple[float, ...],
    *,
    water_arm_lambda: float = DRY_SHELL_WATER_ARM_ACTIVE,
) -> BinaryNoInertCellEvaluation:
    """Evaluate the Law-2 cell rows and their analytic Jacobian.

    ``water_arm_lambda`` is consumed ONLY by the ten-unknown dry-shell branch
    and is 1.0 - the water arm ACTIVE, which owner ruling D9-e (2026-09-13)
    made the asked-for closure - everywhere the cell is solved.  At 1.0 the
    water arm is the shipped two-film one at unit conductance, the eighth row
    is the shipped water isofugacity equality held at the SHELL's sorbed
    activity, and the water latent heat is carried in the surface energy row:
    that is closure (vi) of the closure record's section 8 exactly.  At 0.0
    the water arm is blocked and the eighth row is the binary identity: the
    WITHDRAWN closure (v), which D9-e retains only so that a test or a
    diagnostic can evaluate it - no solver route passes through it any more,
    the D9-b water-blocking walk having been RETIRED at D9-e level 0.  Every
    reported dry-shell root is at 1.0, solved and revalidated there by the
    shipped damped Newton.  The eight-unknown path never reads it, at either
    value.
    """

    # D9-b-1 (owner ruling 2026-09-13): ONE evaluator carries the shipped
    # eight-unknown binary state and the ten-unknown dry-shell receding-front
    # state.  The width is read off the trial vector, and every expression the
    # eight-unknown path evaluates is reached unchanged, in the same order and
    # on the same operands - which is what keeps that path bit-identical.
    dry_shell = len(vector) == len(DRY_SHELL_UNKNOWN_NAMES)
    if dry_shell and (inputs.core_supplied_hexane is None or inputs.core_supplied_water is None):
        # D9-e FIX (2026-09-14), declared row
        # ``d9e_carried_seed_width_from_record``: the ten-row branch is the
        # DOUBLY-SORBED state's, and every row it adds is written on the two
        # core-supplied closures - the shell mass row on the hexane closure's
        # front fraction at S4 below, the eighth row on the water closure's
        # sorbed activity.  A ten-wide vector on a record that binds either
        # FILM law is therefore not a state to evaluate but a width that does
        # not belong to this record, and it is refused TYPED here rather than
        # dereferenced.  Reachable only through a caller that bypasses the
        # entry point's projection; the shipped route can no longer produce it.
        raise UnsupportedInventoryBranchRefusal(
            "the ten-unknown dry-shell trial vector was handed to a record "
            "that does not present the dry-shell receding-front interface "
            "(hexane closure "
            f"{'absent' if inputs.core_supplied_hexane is None else 'present'}, "
            "water closure "
            f"{'absent' if inputs.core_supplied_water is None else 'present'}); "
            "the record decides the width, never the seed"
        )
    if dry_shell:
        _require_dry_shell_state_domain(vector, error_type=EngineeringStepRejected)
    else:
        _require_state_domain(vector, error_type=EngineeringStepRejected)
    variables = tuple(
        _Dual.variable(value, index, width=len(vector)) for index, value in enumerate(vector)
    )
    t_g, y_h_bulk, pressure, t_interface, y_h_interface, s_h, s_w, n_total = variables[:8]
    t_front = variables[8] if dry_shell else None
    a_h_surface = variables[9] if dry_shell else None

    properties = inputs.properties
    transport = inputs.transport
    boundary = inputs.gas_boundary
    area_law = inputs.film_area_law
    transfer = inputs.transfer

    s_total = s_h + s_w
    area = inputs.geometry.gas_side_reference_area_m2 * ((1.0 - s_total) ** area_law.exponent.value)
    if area.value <= 0.0:
        raise EngineeringStepRejected("the gas-side active area is not positive")

    molar_mass_bulk = (
        y_h_bulk * properties.hexane_molar_mass_kg_mol
        + (1.0 - y_h_bulk) * properties.water_molar_mass_kg_mol
    )
    molar_mass_interface = (
        y_h_interface * properties.hexane_molar_mass_kg_mol
        + (1.0 - y_h_interface) * properties.water_molar_mass_kg_mol
    )
    density_bulk = pressure * molar_mass_bulk / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * t_g)
    density_interface = (
        pressure * molar_mass_interface / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * t_interface)
    )
    molar_density = pressure / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * t_interface)
    diffusivity = (
        transport.fsg_diffusivity_multiplier
        * FSG_DIFFUSIVITY_PREFACTOR
        * t_interface**FSG_TEMPERATURE_EXPONENT
        / pressure
    )
    schmidt = transport.gas_dynamic_viscosity_pa_s / (density_interface * diffusivity)

    inlet_molar_mass = (
        boundary.inlet_hexane_mole_fraction * properties.hexane_molar_mass_kg_mol
        + boundary.inlet_water_mole_fraction * properties.water_molar_mass_kg_mol
    )
    inlet_mass_flow = boundary.inlet_molar_flow_mol_s * inlet_molar_mass
    bed_void = area_law.bed_void_fraction
    # D20, 2026-09-23 (owner ruling "I rule D19 and D20 as recommended";
    # packet docs/GT_PS2_RULING_PACKET_D20_VOIDAGE_REYNOLDS_2026-09-22.md).
    # The voidage Reynolds number of this correlation divides by the SOLID
    # fraction.  Faner (2008) PhD thesis, printed page 106, Equation 4.23:
    #
    #     Re_eps = D_EE rho_g v_g / [mu_g (1 - eps)]
    #            = D_EE W_g / [mu_g A_t (1 - eps)],
    #
    # which the source's own Table 4.7 confirms: over its 36 rows the
    # voidage implied by the printed Re/Re_eps pair agrees with the voidage
    # implied by the printed Nu/Nu_eps pair (Equation 4.21) to 6.1e-4, while
    # the divide-by-eps reading misses by up to 0.47.  The digitized source
    # is pinned by
    # paper/analysis/datasets/faner2008_thesis_4_5_3/DATASET_RECORD.md.
    # WAS, until 2026-09-23: the same expression divided by ``bed_void``,
    # the convention written in the retired core without the primary source
    # in hand.  At the frozen eps = 0.40 that form returns 1.5 times this
    # Re_eps and, through the 0.579 exponent, a film coefficient 1/0.7908
    # times this one.
    reynolds_void = (
        inlet_mass_flow
        * transport.hydraulic_particle_diameter_m
        / (
            inputs.geometry.cross_section_m2
            * transport.gas_dynamic_viscosity_pa_s
            * (1.0 - bed_void)
        )
    )
    prandtl = (
        transport.gas_mass_heat_capacity_j_kg_k
        * transport.gas_dynamic_viscosity_pa_s
        / transport.gas_thermal_conductivity_w_m_k
    )
    nusselt = (
        COLETTO_B7_NUSSELT_PREFACTOR
        * reynolds_void**COLETTO_B7_REYNOLDS_EXPONENT
        * prandtl**COLETTO_B7_PRANDTL_EXPONENT
    )
    heat_coefficient = (
        nusselt
        * transport.gas_thermal_conductivity_w_m_k
        * (1.0 - bed_void)
        / (transport.hydraulic_particle_diameter_m * bed_void)
    )
    mass_coefficient = (
        heat_coefficient
        / (density_interface * transport.gas_mass_heat_capacity_j_kg_k)
        * (prandtl / schmidt) ** CHILTON_COLBURN_ANALOGY_EXPONENT
    )
    molar_transport = molar_density * mass_coefficient
    peclet = n_total / molar_transport
    flux_kernel = _bernoulli(peclet)
    core_supplied = inputs.core_supplied_hexane
    core_supplied_water = inputs.core_supplied_water
    if core_supplied is not None and core_supplied_water is not None:
        # C8 Tier 2a, the doubly-sorbed flux law.  Each species would carry
        # its own constitutive scaling on its whole flux, in the same
        # whole-flux form B1 stage 3c settled for a single arm:
        #     n_h = f_h * (y_I n_t + D_h)
        #     n_w = f_w * ((1 - y_I) n_t + D_w)
        # with D_w = -D_h in a binary film.  The binary identity
        # n_t = n_h + n_w is then NOT automatic, and closing it gives the
        # closed form of design section 1 item 3, which
        # :func:`doubly_sorbed_total_molar_flux` still carries as algebra:
        #     n_t [1 - f_h y_I - f_w (1 - y_I)] = f_h D_h + f_w D_w.
        # Q-T2A-1 withdrew that below-unit path from the cell; only the
        # factor-one branch is formulated here.
        hexane_factor = core_supplied.hexane_conductance_factor
        water_factor = core_supplied_water.water_conductance_factor
        if hexane_factor == 1.0 and water_factor == 1.0:
            # At f_h = f_w = 1 the bracket is exactly zero AND the right side
            # is exactly zero, so the identity is 0 = 0: it does not define
            # n_t, which stays the free unknown the energy rows determine -
            # exactly as in the two-film law.  The closed form is therefore
            # NOT evaluated at all on this path, and the three expressions
            # below are the two-film branch verbatim, so the working prior
            # (both factors one) is bit-identical to the two-film law.
            hexane_diffusive = molar_transport * flux_kernel * (y_h_interface - y_h_bulk)
            n_hexane = y_h_interface * n_total + hexane_diffusive
            n_water = n_total - n_hexane
        else:
            # D9-b (owner ruling 2026-09-13) replaces the F-T2A-1 refusal that
            # stood here since Q-T2A-1 (2026-09-06).  The state is a layer
            # whose water film has drained, whose hexane front has receded,
            # and whose hexane still evaporates AT the front: the DRY-SHELL
            # RECEDING-FRONT INTERFACE, formulated and measured in
            # docs/evidence/ft2a1_closure_2026-09-12/ft2a1_shell.py and
            # reported in GT_PS2_FT2A1_CLOSURE_MEASUREMENT_RECORD_2026-09-12.md
            # section 3b: the control with the water arm as shipped converges
            # in 3 steps to 338.61 K against the shipped 338.46 K, and a root
            # inside the declared domain exists from f = 0.95 (a 44 micrometre
            # shell) at 357.6 K with the surface hexane activity 0.058.
            #
            # D9-b-2, THE GUARD, first (brief item (f)): a state with no
            # shell, and a closure reporting a receded front with an
            # unhindered film, are both refused typed before anything is
            # solved.  The reasons are in the guard's own docstring; its
            # firing on any D9-a lane is kill criterion K8.
            _require_dry_shell_is_reachable(core_supplied)
            if not dry_shell:
                raise UnsupportedInventoryBranchRefusal(
                    "the dry-shell receding-front interface carries ten "
                    "unknowns and this evaluation was handed the eight-unknown "
                    "binary state; solve it through "
                    "solve_dry_shell_receding_front_fast_block, which extends "
                    "the seed, or through solve_binary_no_inert_fast_block, "
                    "which dispatches to it"
                )
            # D9-b-1, THE INTERFACE (brief item (e)).  Ten unknowns - the
            # shipped eight plus the FRONT temperature and the SURFACE hexane
            # activity - and ten rows.  They are the seven replacements
            # ``ft2a1_shell.REPLACEMENTS`` makes to THIS function, S1 to S7,
            # written out here as a real branch instead of an exec-ed copy:
            #
            #   S1  ten unknowns unpacked, ten-wide domain check (above)
            #   S2  unscaled gas-film flux, the water arm blocked (here)
            #   S3  the surface hexane activity is an UNKNOWN, not the
            #       pinned 1.0
            #   S4  the shell conductances and the front state
            #   S5  surface energy: convection against the front's duty
            #   S6  rows 8, 9 and 10
            #   S7  the three extra row scales
            #
            # S2.  The gas-film whole flux is UNSCALED: the receding-front
            # hindrance is no longer a factor on this flux, it is the shell
            # resistance written out in rows 9 and 10.  The water arm is
            # ACTIVE (owner ruling D9-e, 2026-09-13): the hexane-free shell
            # the front has opened is exactly what row D9-a-3 lets the layer
            # exchange water through, so the water flux is the shipped
            # two-film one at unit conductance, n_w = n_t - n_h, and the
            # eighth row below is the shipped water isofugacity equality at
            # the shell's sorbed activity.  Until D9-e the water was BLOCKED
            # here and n_w was exactly zero; that form is closure (v), which
            # the ruling withdrew after the D9-3 lanes measured it stalling,
            # and it survives only at ``water_arm_lambda = 0`` as the fallback
            # walk's waypoint.
            hexane_diffusive = molar_transport * flux_kernel * (y_h_interface - y_h_bulk)
            n_hexane = y_h_interface * n_total + hexane_diffusive
            n_water = (n_total - n_hexane) * water_arm_lambda
    elif core_supplied_water is not None:
        # B1 stage 4 (water symmetry): the constitutive law is written on
        # the CORE-SUPPLIED species with the series factor on the WHOLE
        # flux (advective + diffusive); hexane is exact binary bookkeeping.
        n_water = core_supplied_water.water_conductance_factor * (
            (1.0 - y_h_interface) * n_total
            + molar_transport * flux_kernel * ((1.0 - y_h_interface) - (1.0 - y_h_bulk))
        )
        n_hexane = n_total - n_water
    elif core_supplied is not None:
        # B1 stage 3c (supersedes the stage-3b note that left the advective
        # term unscaled): in this NON-dilute system the interface row pins
        # y_I near saturation, so an unscaled Stefan/advective term is an
        # unphysical leak AROUND the internal front resistance - measured
        # as a Newton degeneracy at the DT-scale Biot number.  The van Meel
        # normalized-curve reading scales the WHOLE constitutive flux: the
        # hexane reaching the gas side is the factor times what the film
        # law could carry.  Factor 1.0 reproduces the two-film law exactly.
        n_hexane = core_supplied.hexane_conductance_factor * (
            y_h_interface * n_total + molar_transport * flux_kernel * (y_h_interface - y_h_bulk)
        )
        n_water = n_total - n_hexane
    else:
        hexane_diffusive = molar_transport * flux_kernel * (y_h_interface - y_h_bulk)
        n_hexane = y_h_interface * n_total + hexane_diffusive
        n_water = n_total - n_hexane

    inlet_flow = boundary.inlet_molar_flow_mol_s
    outlet_flow = inlet_flow + area * n_total
    if outlet_flow.value <= 0.0:
        raise EngineeringStepRejected("the trial outlet molar flow is not positive")
    superficial_velocity = (
        outlet_flow * molar_mass_bulk / (density_bulk * inputs.geometry.cross_section_m2)
    )
    if density_bulk.value <= 0.0 or superficial_velocity.value < 0.0:
        raise EngineeringStepRejected("a trial gas density or superficial velocity is invalid")
    try:
        pressure_drop = cc.layer_pressure_drop(
            series=inputs.hydraulics.series,
            superficial_velocity_m_s=superficial_velocity.value,
            gas_density_kg_m3=density_bulk.value,
        )
    except cc.CellClosureError as error:
        raise EngineeringStepRejected(f"pressure-loss closure refused: {error}") from error
    pressure_drop_dual = (
        # D9-b-1: the only bare dual in this evaluator, so it is the only one
        # that must be told the width.  ``(0.0,) * 8`` is what the shipped
        # default builds, so the eight-unknown path is unchanged.
        _Dual(pressure_drop.drop_pa, (0.0,) * len(variables))
        + (superficial_velocity - superficial_velocity.value)
        * pressure_drop.d_drop_d_superficial_velocity
        + (density_bulk - density_bulk.value) * pressure_drop.d_drop_d_gas_density
    )

    if type(properties) is cc.DeclaredGasProperties:
        heat_capacity = properties.molar_heat_capacity_j_mol_k
        datum = properties.energy_datum_temperature_k
        hexane_latent = properties.hexane_latent_heat_j_mol
        water_latent = properties.water_latent_heat_j_mol
        # Use the same common-datum gas enthalpy convention as cell_closure:
        # sensible heat plus composition-weighted species offsets at both external
        # boundaries. Omitting these terms while retaining them on the phase source
        # makes the cell solution depend on an arbitrary species datum shift.
        h_gas = (
            heat_capacity * (t_g - datum)
            + y_h_bulk * hexane_latent
            + (1.0 - y_h_bulk) * water_latent
        )
        h_inlet = (
            heat_capacity * (boundary.inlet_temperature_k - datum)
            + boundary.inlet_hexane_mole_fraction * hexane_latent
            + boundary.inlet_water_mole_fraction * water_latent
        )
        h_hexane = heat_capacity * (t_interface - datum) + hexane_latent
        h_water = heat_capacity * (t_interface - datum) + water_latent
        ackermann = n_total * heat_capacity / heat_coefficient
        scale_heat_capacity = heat_capacity
    else:
        # CELL-02b native arm (PHY-057): every caloric quantity is a frozen-law
        # evaluation on the native zero-gauge datum. The energy rows stay
        # datum-invariant by construction - the gas row consumes enthalpy
        # differences and the interface row consumes saturation latents.
        try:
            bulk_enthalpy = nc.mixture_enthalpy(t_g.value, pressure.value, y_h_bulk.value)
            inlet_enthalpy = nc.mixture_enthalpy(
                boundary.inlet_temperature_k,
                pressure.value,
                boundary.inlet_hexane_mole_fraction,
            )
            hexane_partial = nc.hexane_partial_enthalpy(
                t_interface.value, pressure.value, y_h_interface.value
            )
            water_partial = nc.water_partial_enthalpy(
                t_interface.value, pressure.value, y_h_interface.value
            )
            film_heat_capacity = nc.mixture_heat_capacity(
                t_interface.value, pressure.value, y_h_interface.value
            )
            hexane_latent_law = nc.hexane_latent_heat(t_interface.value)
            water_latent_law = nc.water_latent_heat(t_interface.value)
            # Row scaling only: a solve-constant cp at the declared inlet
            # boundary state, so row scales cannot drift between iterates.
            scale_heat_capacity = nc.mixture_heat_capacity(
                boundary.inlet_temperature_k,
                boundary.downstream_boundary_pressure_pa,
                boundary.inlet_hexane_mole_fraction,
            ).value_j_mol_k
        except (ValueError, OverflowError) as error:
            raise EngineeringStepRejected(f"native caloric law refused: {error}") from error
        h_gas = _lift_partials(
            bulk_enthalpy.value_j_mol,
            (
                (t_g, bulk_enthalpy.d_dT_j_mol_k),
                (pressure, bulk_enthalpy.d_dp_j_mol_pa),
                (y_h_bulk, bulk_enthalpy.d_dy_hexane_j_mol),
            ),
        )
        h_inlet = _lift_partials(
            inlet_enthalpy.value_j_mol,
            ((pressure, inlet_enthalpy.d_dp_j_mol_pa),),
        )
        h_hexane = _lift_partials(
            hexane_partial.value_j_mol,
            (
                (t_interface, hexane_partial.d_dT_j_mol_k),
                (pressure, hexane_partial.d_dp_j_mol_pa),
                (y_h_interface, hexane_partial.d_dy_hexane_j_mol),
            ),
        )
        h_water = _lift_partials(
            water_partial.value_j_mol,
            (
                (t_interface, water_partial.d_dT_j_mol_k),
                (pressure, water_partial.d_dp_j_mol_pa),
                (y_h_interface, water_partial.d_dy_hexane_j_mol),
            ),
        )
        hexane_latent = _lift_partials(
            hexane_latent_law.value_j_mol,
            ((t_interface, hexane_latent_law.d_dT_j_mol_k),),
        )
        water_latent = _lift_partials(
            water_latent_law.value_j_mol,
            ((t_interface, water_latent_law.d_dT_j_mol_k),),
        )
        heat_capacity_dual = _lift_partials(
            film_heat_capacity.value_j_mol_k,
            (
                (t_interface, film_heat_capacity.d_dT_j_mol_k2),
                (pressure, film_heat_capacity.d_dp_j_mol_k_pa),
                (y_h_interface, film_heat_capacity.d_dy_hexane_j_mol_k),
            ),
        )
        ackermann = n_total * heat_capacity_dual / heat_coefficient
    heat_kernel = _bernoulli(ackermann)
    gas_to_interface = heat_coefficient * area * (t_g - t_interface) * heat_kernel
    wall_to_gas = transfer.wall_to_gas_ua_w_k * (inputs.wall_temperature_k - t_g)
    wall_to_interface = transfer.wall_to_interface_ua_w_k * (
        inputs.wall_temperature_k - t_interface
    )
    solid_to_interface = transfer.solid_to_interface_ua_w_k * (
        inputs.solid.solid_temperature_k - t_interface
    )

    phi_hexane = _vapor_fugacity_coefficient(t_interface, pressure, y_h_interface, species="hexane")
    phi_water = _vapor_fugacity_coefficient(t_interface, pressure, y_h_interface, species="water")
    liquid_hexane = _liquid_reference_fugacity(hexane_props, t_interface, pressure)
    liquid_water = _liquid_reference_fugacity(water_props, t_interface, pressure)

    declared_hexane, declared_water = declared_external_saturations(inputs)
    if dry_shell:
        # S4 (D9-b-1): the dry shell between f R and R, both conductances
        # referred to the OUTER sphere area, and the front state.
        #
        #   k_heat = k_meal f / (R (1 - f))        [W/(m2 K)]
        #   k_mass = k_gas   f / (Bi_m (1 - f))    [m/s]
        #
        # The mass conductance is NOT taken from the Biot number: this module
        # may not import the kernel that declares it (:112) without a cycle.
        # It is read off the closure's OWN series factor through the identity
        # ``ft2a1_shell.check_phi_identity`` asserts numerically,
        #
        #   Phi = (1/k_gas)/(1/k_gas + 1/k_mass)  =>  k_mass = k_gas Phi/(1-Phi)
        #
        # which is algebraically exact, because Phi/(1 - Phi) = f/(Bi_m (1-f))
        # for Phi = 1/(1 + Bi_m (1-f)/f).  The series over the unhindered gas
        # film therefore reproduces the shipped falling-rate factor by
        # CONSTRUCTION rather than by agreement.
        #
        # CORRECTION to the D9-1 proposal header, measured in this build: the
        # header wrote ``molar_transport * Phi/(1 - Phi)``.  The shell mass row
        # multiplies this conductance by a CONCENTRATION difference in mol/m3
        # to produce a molar flux, so the conductance is a velocity in m/s and
        # the gas-film conductance it is built from is ``mass_coefficient``,
        # which is exactly the ``k_gas`` the harness asserts the identity on.
        # ``molar_transport`` is mol/(m2 s) and would leave the row in the
        # wrong unit by a factor of the molar density.
        front_fraction = core_supplied.front_fraction
        shell_heat_conductance = (
            DRY_SHELL_MEAL_CONDUCTIVITY_W_M_K
            * front_fraction
            / (DRY_SHELL_PARTICLE_RADIUS_M * (1.0 - front_fraction))
        )
        shell_mass_conductance = mass_coefficient * (
            core_supplied.hexane_conductance_factor
            / (1.0 - core_supplied.hexane_conductance_factor)
        )
        front_saturation_pressure_pa = _dry_shell_front_saturation_pressure(t_front)
        front_concentration = front_saturation_pressure_pa / (
            UNIVERSAL_GAS_CONSTANT_J_MOL_K * t_front
        )
        surface_concentration = molar_density * y_h_interface
        front_latent = _dry_shell_front_latent_heat(t_front)
    if core_supplied is None:
        hexane_equilibrium_row = y_h_interface * phi_hexane * pressure - liquid_hexane
    elif dry_shell:
        # S3: the hexane isofugacity row with the SURFACE activity as an
        # unknown.  Behind the front the liquid is at its own saturation; what
        # the gas sees at the surface is the unknown the shell mass row fixes.
        hexane_equilibrium_row = y_h_interface * phi_hexane * pressure - a_h_surface * liquid_hexane
    else:
        # B1 stage 3: the sorbed-phase activity scales the pure-liquid
        # reference - the interface equilibrates with core-held hexane, not
        # with a free film.
        hexane_equilibrium_row = (
            y_h_interface * phi_hexane * pressure - core_supplied.hexane_activity * liquid_hexane
        )
    if core_supplied_water is None:
        water_equilibrium_row = (1.0 - y_h_interface) * phi_water * pressure - liquid_water
    else:
        # B1 stage 4: the PHY-031 Luikov retained-water activity scales the
        # pure-liquid reference - the interface equilibrates with core-held
        # (sorbed) water, not with a free film.
        water_equilibrium_row = (
            1.0 - y_h_interface
        ) * phi_water * pressure - core_supplied_water.water_activity * liquid_water
    if dry_shell:
        # S5.  The surface energy balance, RECOMBINED.  The pair the design
        # brief writes is
        #     (4)   gas_to_interface = A k_heat (T_I - T_f)
        #     (10)  A k_heat (T_I - T_f) = A n_h lambda_h(T_f)
        # and what is carried here is their sum for row 4 and row 10 divided
        # through by ``A k_heat``, so that it reads Fourier's law solved for
        # the DROP.  The root set is identical - one row addition and one row
        # rescale - and the recombination is what makes the block solvable: at
        # a thin shell ``A k_heat`` is of order 1e13 W/K and destroys the
        # Jacobian scaling of the written pair (measurement record, sec. 3b).
        # The WATER LATENT HEAT IS KEPT (owner ruling D9-e): the condensing
        # steam is what heats the dry surface from the co-boiling lock toward
        # the steam temperature as the front recedes, and dropping it was the
        # withdrawn closure (v).  The shipped wall and solid UA terms are
        # DROPPED, which is the primary closure of that record, whose (v-b)
        # variant keeps them and moved the captured root by 0.38 K - a
        # declared modelling level, not a measured one.
        interface_energy_row = gas_to_interface - area * (
            n_hexane * front_latent + water_arm_lambda * n_water * water_latent
        )
        # S6.  Rows 8, 9 and 10: the water isofugacity equality at the shell's
        # sorbed activity (divided through by the pressure, which is its row
        # scale below), the shell mass balance, and the front energy row.  At
        # the withdrawn ``water_arm_lambda = 0`` waypoint the eighth row is
        # the binary identity instead, which is what made n_w exactly zero.
        closing_rows = (
            (1.0 - water_arm_lambda) * (n_total - n_hexane - n_water) / molar_transport
            + water_arm_lambda * water_equilibrium_row / pressure,
            n_hexane - shell_mass_conductance * (front_concentration - surface_concentration),
            (t_interface - t_front) - n_hexane * front_latent / shell_heat_conductance,
        )
    else:
        interface_energy_row = (
            gas_to_interface
            + wall_to_interface
            + solid_to_interface
            - area * (n_hexane * hexane_latent + n_water * water_latent)
        )
        # C8 Tier 2a with Q-T2A-1 (2026-09-06), amended by D9-b: the
        # doubly-sorbed branch below unit conductance now takes the dry-shell
        # branch above rather than refusing, so no EIGHT-unknown path releases
        # this row and the eighth row is the water isofugacity equality on
        # every one of them.
        closing_rows = (water_equilibrium_row,)
    rows = (
        outlet_flow * y_h_bulk - inlet_flow * boundary.inlet_hexane_mole_fraction - area * n_hexane,
        # F1 SIGN FIX (2026-09-02, equations-honesty defect #1) - the
        # ``+ gas_to_interface`` on the fourth line of the gas energy row
        # below.  That row is written (out - in - sources), so a SOURCE into
        # the gas carries a minus (the arriving vapour enthalpy, the wall
        # duty).  ``gas_to_interface`` = h*A*(T_g - T_I)*kernel is heat
        # LEAVING the gas across the interface - a SINK - and must therefore
        # be ADDED, exactly once, opposite to the ``+ gas_to_interface``
        # source it is given in the interface row two rows further down.
        # Until 2026-09-02 it was SUBTRACTED here as well, booking the same
        # q_gI as a source to BOTH control volumes and fabricating 2*q_gI per
        # cell (measured -1,592,790.83 W on this fixture).  That violated
        # registry INVARIANT rows BAL-ENERGY-S (:56), BAL-ENERGY-G (:68) and
        # FLUX-ENERGY-INTERFACE (:75-84), whose declared statement pairs
        # -A*J_E on the solid side with +A*J_E on the gas side.  The certified
        # cell_closure arm has always had it right (cell_closure.py:2944 vs
        # :2950), as do all four native source kernels.  Enforced by
        # tests/test_core2_interface_energy_equal_and_opposite.py.
        outlet_flow * h_gas
        - inlet_flow * h_inlet
        - area * (n_hexane * h_hexane + n_water * h_water)
        + gas_to_interface
        - wall_to_gas,
        pressure - boundary.downstream_boundary_pressure_pa - pressure_drop_dual,
        interface_energy_row,
        s_h - declared_hexane,
        s_w - declared_water,
        hexane_equilibrium_row,
        *closing_rows,
    )
    residual = tuple(row.value for row in rows)
    jacobian = tuple(row.partials for row in rows)
    energy_scale = inlet_flow * scale_heat_capacity * 100.0
    if dry_shell:
        # S7.  Row 8 is already divided through by the pressure inside the
        # row - which is exactly the pressure scale the shipped eighth row
        # carries - so its scale here is one; the shell mass row is scaled by
        # the molar transport, and the front energy row is in kelvin.
        row_scales = (
            inlet_flow,
            energy_scale,
            1.0e4,
            energy_scale,
            1.0,
            1.0,
            pressure.value,
            1.0,
            molar_transport.value,
            1.0,
        )
    else:
        row_scales = (
            inlet_flow,
            energy_scale,
            1.0e4,
            energy_scale,
            1.0,
            1.0,
            pressure.value,
            pressure.value,
        )
    scaled_norm = max(abs(value) / scale for value, scale in zip(residual, row_scales, strict=True))
    independent_water_balance = (
        outlet_flow.value * (1.0 - y_h_bulk.value)
        - inlet_flow * boundary.inlet_water_mole_fraction
        - area.value * n_water.value
    )
    evaluation_type = DryShellCellEvaluation if dry_shell else BinaryNoInertCellEvaluation
    return evaluation_type(
        residual=residual,
        jacobian=jacobian,
        row_scales=row_scales,
        scaled_residual_norm=scaled_norm,
        gas_side_active_area_m2=area.value,
        outlet_molar_flow_mol_s=outlet_flow.value,
        total_molar_flux_mol_m2_s=n_total.value,
        hexane_molar_flux_mol_m2_s=n_hexane.value,
        water_molar_flux_mol_m2_s=n_water.value,
        gas_density_kg_m3=density_bulk.value,
        superficial_velocity_m_s=superficial_velocity.value,
        layer_pressure_drop_pa=pressure_drop.drop_pa,
        binary_diffusivity_m2_s=diffusivity.value,
        schmidt_number=schmidt.value,
        mass_transfer_coefficient_m_s=mass_coefficient.value,
        molar_transport_coefficient_mol_m2_s=molar_transport.value,
        heat_transfer_coefficient_w_m2_k=heat_coefficient,
        wall_to_gas_w=wall_to_gas.value,
        wall_to_interface_w=wall_to_interface.value,
        solid_to_interface_w=solid_to_interface.value,
        interface_convection_w=gas_to_interface.value,
        independent_water_balance_residual_mol_s=independent_water_balance,
    )


#: D9-e: the two field-state types the acceptance surface reconstructs from,
#: as EXACT types.  The dry-shell state is a subclass of the binary one, so an
#: ``isinstance`` here would also admit any future subclass silently; the gate
#: stays an exact-membership test and this tuple is the whole of what it
#: admits.  Nothing else about the check moved.
ACCEPTED_FAST_BLOCK_STATE_TYPES = (BinaryNoInertCellFieldState, DryShellCellFieldState)
#: The two evaluation types the same surface accepts, on the same terms.
ACCEPTED_FAST_BLOCK_EVALUATION_TYPES = (BinaryNoInertCellEvaluation, DryShellCellEvaluation)


def assemble_binary_no_inert_fast_block(
    inputs: EngineeringCellInputs, state: BinaryNoInertCellFieldState | DryShellCellFieldState
) -> BinaryNoInertCellEvaluation | DryShellCellEvaluation:
    """Re-evaluate a fast block from its inputs and its accepted state.

    D9-e widened this by exactly one type.  The ten-unknown dry-shell state
    carries its own ``as_vector()``, and the evaluator dispatches on that
    vector's WIDTH, so the reconstruction of a dry-shell step is the same call
    on the same operands as the reconstruction of a binary one - which is what
    lets the acceptance surface run its substance checks on both without a
    second code path.  Any other type is still refused typed.
    """

    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    if type(state) not in ACCEPTED_FAST_BLOCK_STATE_TYPES:
        raise EngineeringConfigurationError(
            "state must be an exact BinaryNoInertCellFieldState or DryShellCellFieldState"
        )
    return _evaluate_binary_no_inert(inputs, state.as_vector())


#: The unmodified nominal bulk-gas hexane mole fraction guess, kept as a
#: named constant so the Q-F1b interpolation below (and its bit-identity
#: argument) reads directly off the same literal every other branch uses.
_NOMINAL_GAS_HEXANE_MOLE_FRACTION = 0.78


def nominal_binary_no_inert_seed(
    inputs: EngineeringCellInputs,
) -> BinaryNoInertCellFieldState:
    """Return a declared-domain seed, not a fitted or qualified state."""

    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    s_h, s_w = declared_external_saturations(inputs)
    core_supplied = inputs.core_supplied_hexane
    if core_supplied is None:
        gas_hexane_mole_fraction = _NOMINAL_GAS_HEXANE_MOLE_FRACTION
    else:
        # Q-F1b (owner ruling 2026-09-02, batch commit 20f5028).
        #
        # On the CORE_SUPPLIED hexane branch, the whole-flux B1 stage-3c
        # law (:1607-1619 below) scales the ENTIRE constitutive hexane
        # flux by ``hexane_conductance_factor``:
        #   n_hexane = factor * (y_h_interface*n_total + diffusive term)
        # so as factor -> 0, n_hexane -> 0 identically (no hexane crosses
        # the interface at all) - the physical "no-transfer equilibrium".
        # The bulk-gas species row (:1784 below) is
        #   outlet_flow*y_h_bulk - inlet_flow*inlet_hexane_mole_fraction
        #   - area*n_hexane = 0
        # so at that limit it relaxes to
        #   outlet_flow*y_h_bulk = inlet_flow*inlet_hexane_mole_fraction,
        # i.e. the no-transfer limit for the BULK gas composition
        # (``gas_hexane_mole_fraction``, y_h_bulk) is simply the
        # boundary's own inlet hexane mole fraction - there is no local
        # source left to enrich it above that.  The nominal seed's 0.78
        # guess presumes active hexane enrichment; at a strongly
        # hindered factor (e.g. 0.02) that guess sits on the wrong side
        # of the solution fold, which is the measured stiffness this
        # ruling responds to (docs/GT_PS2_F1_SIGN_FIX_CEREMONY_2026-09-02.md
        # Sec 7.3, docs/GT_PS2_F5_PHI_SUPERSESSION_RETRO_SIGN_2026-09-02.md
        # Sec 5.3).
        #
        # Deliberately NOT applied to ``interface_hexane_mole_fraction``:
        # that field is pinned by a SEPARATE isofugacity/activity
        # equilibrium row (``hexane_equilibrium_row``, :1766-1773 below)
        # that depends on ``core_supplied.hexane_activity``, T_I and P -
        # not on ``hexane_conductance_factor`` at all - so its converged
        # value does not move with the factor and no fold crossing is
        # expected there; leaving it at the nominal guess is correct, not
        # an oversight.
        #
        # Interpolated LINEARLY by the factor itself - the same
        # parameter the governing equation itself scales the flux by -
        # between the nominal guess (weight = factor) and the no-transfer
        # limit (weight = 1 - factor).  This is provably bit-identical to
        # the unconditional nominal seed at factor = 1.0 (the previously
        # converging path this ruling must not disturb): in binary64,
        # ``1.0 * nominal + 0.0 * limit`` recombines to exactly
        # ``nominal`` for any finite ``limit`` - asserted, not merely
        # argued, in
        # tests/test_core2_core_supplied_hexane_closure.py::
        # test_nominal_seed_is_bit_identical_at_unhindered_conductance_factor.
        no_transfer_limit = inputs.gas_boundary.inlet_hexane_mole_fraction
        factor = core_supplied.hexane_conductance_factor
        gas_hexane_mole_fraction = (
            factor * _NOMINAL_GAS_HEXANE_MOLE_FRACTION + (1.0 - factor) * no_transfer_limit
        )
    return BinaryNoInertCellFieldState(
        gas_temperature_k=336.0,
        gas_hexane_mole_fraction=gas_hexane_mole_fraction,
        layer_pressure_pa=104_000.0,
        interface_temperature_k=335.5,
        interface_hexane_mole_fraction=0.79,
        external_hexane_saturation=s_h,
        external_water_saturation=s_w,
        total_molar_flux_mol_m2_s=0.002,
    )


def presents_dry_shell_receding_front_interface(inputs: EngineeringCellInputs) -> bool:
    """D9-b-1: does this record present the doubly-sorbed state below unit
    conductance - the state the dry-shell receding-front interface carries?

    Both arms declared and at least one of them hindered.  The factor-one
    corner (both exactly 1.0) stays the shipped two-film-verbatim path and is
    NOT a dry-shell state.
    """

    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    core_supplied = inputs.core_supplied_hexane
    core_supplied_water = inputs.core_supplied_water
    if core_supplied is None or core_supplied_water is None:
        return False
    return not (
        core_supplied.hexane_conductance_factor == 1.0
        and core_supplied_water.water_conductance_factor == 1.0
    )


def dry_shell_surface_activity_seed(state: BinaryNoInertCellFieldState) -> float:
    """The surface activity that makes the isofugacity row EXACTLY the shipped
    one at the seed state: ``a_h,s = y_I phi_h P / f_liquid,h(T_I, P)``.

    ``ft2a1_shell.seed_surface_activity``, which records why: seeding at the
    stamped series factor Phi is wrong by two orders of magnitude - Phi hinders
    the FLUX relative to its driving force, not the surface activity relative
    to saturation - and a Phi seed was MEASURED stalling this closure.
    """

    t_interface = _Dual(state.interface_temperature_k)
    pressure = _Dual(state.layer_pressure_pa)
    y_hexane = _Dual(state.interface_hexane_mole_fraction)
    phi_hexane = _vapor_fugacity_coefficient(t_interface, pressure, y_hexane, species="hexane")
    liquid_hexane = _liquid_reference_fugacity(hexane_props, t_interface, pressure)
    return (y_hexane * phi_hexane * pressure / liquid_hexane).value


def binary_seed_from_dry_shell(seed: DryShellCellFieldState) -> BinaryNoInertCellFieldState:
    """D9-e FIX (2026-09-14): project a carried ten-unknown state back onto
    the eight shipped ones - the exact inverse of
    :func:`dry_shell_seed_from_binary`.

    Declared row ``d9e_carried_seed_width_from_record``.  The ten unknowns
    ARE the eight shipped ones followed by the front temperature and the
    surface hexane activity (``DRY_SHELL_UNKNOWN_NAMES`` is built that
    way), so the projection is a SLICE: no value is recomputed, rescaled
    or averaged, and a round trip through the extension returns the same
    eight components bit for bit.  The two dropped unknowns belong to the
    dry-shell FORMULATION and carry no inventory, so a layer that has left
    the branch loses nothing by dropping them.
    """

    if type(seed) is not DryShellCellFieldState:
        raise EngineeringConfigurationError("seed must be an exact DryShellCellFieldState")
    return BinaryNoInertCellFieldState.from_vector(
        seed.as_vector()[: len(BINARY_NO_INERT_UNKNOWN_NAMES)]
    )


def dry_shell_seed_from_binary(
    inputs: EngineeringCellInputs, seed: BinaryNoInertCellFieldState
) -> DryShellCellFieldState:
    """Extend an eight-unknown seed to the ten-unknown dry-shell state.

    Exactly how ``ft2a1_shell.main`` seeds: the FRONT temperature at the
    interface seed, the surface hexane activity at the isofugacity value of
    the seed state (``vector = tuple(base + [base[3], activity_seed])``).
    """

    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    if type(seed) is not BinaryNoInertCellFieldState:
        raise EngineeringConfigurationError("seed must be an exact BinaryNoInertCellFieldState")
    return DryShellCellFieldState(
        **dict(zip(BINARY_NO_INERT_UNKNOWN_NAMES, seed.as_vector(), strict=True)),
        front_temperature_k=seed.interface_temperature_k,
        surface_hexane_activity=dry_shell_surface_activity_seed(seed),
    )


def nominal_dry_shell_seed(inputs: EngineeringCellInputs) -> DryShellCellFieldState:
    """The declared-domain ten-unknown seed; not a fitted or qualified state."""

    return dry_shell_seed_from_binary(inputs, nominal_binary_no_inert_seed(inputs))


@dataclass(frozen=True, slots=True, kw_only=True)
class BinaryNoInertCellSolve:
    #: D9-b-1: a ``DryShellCellFieldState`` on the ten-unknown dry-shell
    #: branch, the shipped eight-unknown state everywhere else.
    state: BinaryNoInertCellFieldState | DryShellCellFieldState
    evaluation: BinaryNoInertCellEvaluation
    iterations: int
    tolerance: float
    inputs_digest: str
    accepted_iterates: tuple[BinaryNoInertCellFieldState, ...]
    accepted_scaled_residual_norms: tuple[float, ...]
    line_search_rejections: int
    #: D9-e: how many solves the declared HINDRANCE CONTINUATION used to reach
    #: this root, counting the rung on the caller's own record.  ``0`` means
    #: the continuation was not entered at all - either the branch was not the
    #: dry-shell one, or the direct Newton from the caller's seed reached the
    #: root by itself - so a reader can tell the two apart without inference.
    dry_shell_continuation_solves: int = 0

    analytic_jacobian: ClassVar[bool] = True
    rejected_trials_are_never_states: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


DEFAULT_NEWTON_TOLERANCE = 1.0e-10
DEFAULT_NEWTON_MAX_ITERATIONS = 60
DEFAULT_LINE_SEARCH_REDUCTIONS = 32


def _scaled_newton_step(evaluation: BinaryNoInertCellEvaluation) -> tuple[float, ...]:
    row_scale = np.asarray(evaluation.row_scales, dtype=float)
    # D9-b-1: the SAME damped Newton drives both widths; only the declared
    # column scales are selected by the residual length.  At eight this is
    # ``BINARY_NO_INERT_COLUMN_SCALES``, exactly as shipped.
    column_scale = np.asarray(
        DRY_SHELL_COLUMN_SCALES
        if len(evaluation.residual) == len(DRY_SHELL_UNKNOWN_NAMES)
        else BINARY_NO_INERT_COLUMN_SCALES,
        dtype=float,
    )
    residual = np.asarray(evaluation.residual, dtype=float) / row_scale
    jacobian = np.asarray(evaluation.jacobian, dtype=float)
    scaled_jacobian = jacobian * column_scale[np.newaxis, :] / row_scale[:, np.newaxis]
    try:
        scaled_step = np.linalg.solve(scaled_jacobian, -residual)
    except np.linalg.LinAlgError as error:
        raise EngineeringNewtonConvergenceError(
            "the scaled analytic Law-2 Jacobian is singular"
        ) from error
    step = scaled_step * column_scale
    if not np.all(np.isfinite(step)):
        raise EngineeringNewtonConvergenceError("the analytic Newton step is not finite")
    return tuple(float(value) for value in step)


def solve_binary_no_inert_fast_block(
    inputs: EngineeringCellInputs,
    *,
    seed: BinaryNoInertCellFieldState,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
    max_line_search_reductions: int = DEFAULT_LINE_SEARCH_REDUCTIONS,
    water_arm_lambda: float = DRY_SHELL_WATER_ARM_ACTIVE,
) -> BinaryNoInertCellSolve:
    """Solve the bounded 8x8 block by scaled analytic damped Newton.

    D9-b-1: the same driver solves the bounded 10x10 dry-shell block when the
    record presents that state; the width comes from the seed, and nothing
    else about the driver changes.  ``water_arm_lambda`` is the dry-shell
    water-arm homotopy and is 1.0 - the water arm ACTIVE, the closure owner
    ruling D9-e asked for - on every solve whose root is reported.
    """

    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    if type(seed) not in (BinaryNoInertCellFieldState, DryShellCellFieldState):
        raise EngineeringConfigurationError(
            "seed must be an exact BinaryNoInertCellFieldState or DryShellCellFieldState"
        )
    if type(water_arm_lambda) is not float or not 0.0 <= water_arm_lambda <= 1.0:
        raise EngineeringConfigurationError("water_arm_lambda must be a float in [0, 1]")
    if type(seed) is DryShellCellFieldState and not presents_dry_shell_receding_front_interface(
        inputs
    ):
        # D9-e FIX (2026-09-14), declared row
        # ``d9e_carried_seed_width_from_record``: the caller's seed is the
        # PREVIOUS interval's accepted state for this layer, and owner
        # ruling D9-e made the ten-unknown dry-shell state one that may be
        # accepted and carried.  This layer has since LEFT the branch - its
        # packets regained an attached hexane film, an external water film,
        # or both - so the record now presents the shipped eight-unknown
        # state.  Project the carried seed onto those eight and solve the
        # shipped law: the width comes from the RECORD, never from the seed.
        seed = binary_seed_from_dry_shell(seed)
    if type(seed) is BinaryNoInertCellFieldState and presents_dry_shell_receding_front_interface(
        inputs
    ):
        # D9-b-1: the doubly-sorbed state below unit conductance IS the
        # dry-shell receding-front interface, and it carries ten unknowns.
        # The caller's eight-unknown seed is handed to the dry-shell entry
        # point, which extends it and owns the declared seeding route, so
        # every kernel call site reaches the branch without knowing about it.
        return solve_dry_shell_receding_front_fast_block(
            inputs,
            seed=seed,
            tolerance=tolerance,
            max_iterations=max_iterations,
            max_line_search_reductions=max_line_search_reductions,
        )
    state_type = type(seed)
    _require_positive("tolerance", tolerance)
    if type(max_iterations) is not int or max_iterations <= 0:
        raise EngineeringConfigurationError("max_iterations must be a positive int")
    if type(max_line_search_reductions) is not int or max_line_search_reductions <= 0:
        raise EngineeringConfigurationError("max_line_search_reductions must be a positive int")

    vector = seed.as_vector()
    evaluation = _evaluate_binary_no_inert(inputs, vector, water_arm_lambda=water_arm_lambda)
    states = [seed]
    norms = [evaluation.scaled_residual_norm]
    rejected_trials = 0
    iterations = 0
    while evaluation.scaled_residual_norm > tolerance:
        if iterations >= max_iterations:
            raise EngineeringNewtonConvergenceError(
                f"Law 2 did not converge in {max_iterations} accepted Newton steps; "
                f"scaled residual {evaluation.scaled_residual_norm!r}"
            )
        step = _scaled_newton_step(evaluation)
        step_factor = 1.0
        accepted = False
        for _ in range(max_line_search_reductions):
            candidate = tuple(
                float(current + step_factor * increment)
                for current, increment in zip(vector, step, strict=True)
            )
            try:
                candidate_evaluation = _evaluate_binary_no_inert(
                    inputs, candidate, water_arm_lambda=water_arm_lambda
                )
            except EngineeringStepRejected:
                rejected_trials += 1
                step_factor *= 0.5
                continue
            armijo_target = (1.0 - 1.0e-4 * step_factor) * evaluation.scaled_residual_norm
            if (
                candidate_evaluation.scaled_residual_norm <= armijo_target
                or candidate_evaluation.scaled_residual_norm <= tolerance
            ):
                vector = candidate
                evaluation = candidate_evaluation
                state = state_type.from_vector(vector)
                states.append(state)
                norms.append(evaluation.scaled_residual_norm)
                iterations += 1
                accepted = True
                break
            rejected_trials += 1
            step_factor *= 0.5
        if not accepted:
            raise EngineeringNewtonConvergenceError(
                "the analytic Newton direction produced no admissible residual-decreasing "
                "trial; the step is refused without clamping or fallback"
            )

    state = state_type.from_vector(vector)
    final_evaluation = _evaluate_binary_no_inert(
        inputs, state.as_vector(), water_arm_lambda=water_arm_lambda
    )
    if final_evaluation.scaled_residual_norm > tolerance:
        raise EngineeringNewtonConvergenceError(
            "the revalidated Law-2 root does not meet the requested tolerance"
        )
    if inputs.require_corroborated_pressure_drop:
        cc.require_corroborated_layer_drop(final_evaluation.layer_pressure_drop_pa)
    return BinaryNoInertCellSolve(
        state=state,
        evaluation=final_evaluation,
        iterations=iterations,
        tolerance=tolerance,
        inputs_digest=inputs.inputs_digest,
        accepted_iterates=tuple(states),
        accepted_scaled_residual_norms=tuple(norms),
        line_search_rejections=rejected_trials,
    )


def solve_dry_shell_receding_front_fast_block(
    inputs: EngineeringCellInputs,
    *,
    seed: BinaryNoInertCellFieldState | DryShellCellFieldState | None = None,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
    max_line_search_reductions: int = DEFAULT_LINE_SEARCH_REDUCTIONS,
) -> BinaryNoInertCellSolve:
    """D9-b-1: solve the ten-unknown dry-shell receding-front block.

    The named entry point for the branch.  It is the SAME damped Newton
    :func:`solve_binary_no_inert_fast_block` runs - only the width differs -
    and that function dispatches to this state on its own, so a caller that
    already holds an eight-unknown seed need not know the branch exists.  This
    entry point exists for the tests and the diagnostics that want to build
    the ten-unknown seed explicitly.

    Refuses typed if the record does not present the state (both arms
    declared, at least one hindered), so the branch can never be entered by a
    caller that merely guessed.
    """

    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    if not presents_dry_shell_receding_front_interface(inputs):
        raise UnsupportedInventoryBranchRefusal(
            "the dry-shell receding-front interface is the DOUBLY-SORBED state "
            "below unit conductance; this record does not present it"
        )
    _require_dry_shell_is_reachable(inputs.core_supplied_hexane)
    if seed is None:
        seed = nominal_dry_shell_seed(inputs)
    elif type(seed) is BinaryNoInertCellFieldState:
        seed = dry_shell_seed_from_binary(inputs, seed)
    newton = {
        "tolerance": tolerance,
        "max_iterations": max_iterations,
        "max_line_search_reductions": max_line_search_reductions,
    }
    try:
        # D9-e (owner ruling 2026-09-13): the CALLER'S SEED FIRST, at the ruled
        # physics.  With the water arm active the isofugacity seed reaches the
        # root directly at every front the addendum's Table N measured, which
        # is why the continuation below is a fallback and not the route.
        return solve_binary_no_inert_fast_block(inputs, seed=seed, **newton)
    except EngineeringNewtonConvergenceError as direct_refusal:
        # THE DECLARED HINDRANCE CONTINUATION, and the ONLY fallback this entry
        # point has.  Measured at D9-e level 0 on the honest lane's step-40
        # cell: the direct Newton refuses there at EVERY front, because the
        # isofugacity seed puts the surface activity at 1.0 where a series
        # factor of 0.00115 needs about 0.21 and the hexane flux is two and a
        # half orders of magnitude out, so the line search cannot walk the
        # distance.  The root is nonetheless there, and continuing in the
        # HINDRANCE reaches it: solve at a declared start where the caller's
        # own seed still works, then step the series factor geometrically down
        # to the record's own, re-seeding from each rung's root.
        #
        # The front travels WITH the factor, through the mass Biot number the
        # record's own pair already fixes -
        #     Phi = 1/(1 + Bi (1 - f)/f)  =>  Bi = f (1 - Phi)/(Phi (1 - f)),
        # exact, so this module still never imports the kernel that declares
        # Bi - and every rung is therefore a physically consistent state of the
        # SAME particle at a less receded front, not an arbitrary relaxation.
        #
        # Nothing here weakens the answer: the last rung IS the caller's record,
        # so the reported root is the solve at the record's own factor, found
        # and revalidated by the shipped damped Newton at the requested
        # tolerance.  The ladder is a fixed count at a fixed ratio with no
        # adaptive tolerance and no halving, so the same inputs give the same
        # root bit for bit on every call.  If any rung refuses, the refusal is
        # typed and NAMES the rung: there is no second fallback.
        return _dry_shell_hindrance_continuation(inputs, seed, newton, direct_refusal)


def _dry_shell_hindrance_continuation(
    inputs: EngineeringCellInputs,
    seed: DryShellCellFieldState,
    newton: dict[str, object],
    direct_refusal: EngineeringNewtonConvergenceError,
) -> BinaryNoInertCellSolve:
    """D9-e: the declared hindrance continuation - the branch's one fallback.

    Its two numbers are the law module's declared rows
    ``d9e_hindrance_continuation_start`` and ``..._steps``; both are
    SEEDING-ONLY, and the root this returns is the solve on ``inputs`` itself.
    """

    core_supplied = inputs.core_supplied_hexane
    target = core_supplied.hexane_conductance_factor
    front = core_supplied.front_fraction
    start = DRY_SHELL_HINDRANCE_CONTINUATION_START
    steps = DRY_SHELL_HINDRANCE_CONTINUATION_STEPS
    if target >= start:
        raise EngineeringNewtonConvergenceError(
            "the declared hindrance continuation has no room to walk: the "
            f"record's series factor {target!r} is at or above the declared "
            f"start {start!r}, so the ladder would begin at the state that "
            "already refused; the direct refusal stands, unmodified "
            f"({direct_refusal})"
        ) from direct_refusal
    # The mass Biot number the record's own (f, Phi) pair fixes, exactly.
    biot = front * (1.0 - target) / (target * (1.0 - front))
    ratio = (target / start) ** (1.0 / steps)
    state: DryShellCellFieldState = seed
    solve: BinaryNoInertCellSolve | None = None
    for index in range(steps + 1):
        if index == steps:
            # The last rung is the caller's own record, untouched.
            rung = inputs
            rung_factor = target
            rung_front = front
        else:
            rung_factor = start * ratio**index
            rung_front = biot * rung_factor / (biot * rung_factor + 1.0 - rung_factor)
            rung = replace(
                inputs,
                core_supplied_hexane=CoreSuppliedHexaneClosure(
                    hexane_activity=core_supplied.hexane_activity,
                    hexane_conductance_factor=rung_factor,
                    front_fraction=rung_front,
                    authority=core_supplied.authority,
                ),
            )
        try:
            solve = solve_binary_no_inert_fast_block(rung, seed=state, **newton)
        except EngineeringFeasibilityError as rung_refusal:
            raise EngineeringNewtonConvergenceError(
                f"the declared hindrance continuation refused at rung {index} of "
                f"{steps} (series factor {rung_factor!r}, front fraction "
                f"{rung_front!r}, declared start {start!r}): {rung_refusal}"
            ) from rung_refusal
        state = DryShellCellFieldState.from_vector(solve.state.as_vector())
    if solve is None:  # pragma: no cover - the loop always runs at least once
        raise direct_refusal
    return replace(solve, dry_shell_continuation_solves=steps + 1)


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierCellEvaluation:
    residual: tuple[float, ...]
    jacobian: tuple[tuple[float, ...], ...]
    row_scales: tuple[float, ...]
    scaled_residual_norm: float
    gas_side_active_area_m2: float
    outlet_molar_flow_mol_s: float
    interface_hexane_mole_fraction: float
    interface_water_mole_fraction: float
    bulk_carrier_mole_fraction: float
    interface_carrier_mole_fraction: float
    peclet: float
    bernoulli_factor: float
    total_molar_flux_mol_m2_s: float
    laboratory_frame_hexane_molar_flux_mol_m2_s: float
    laboratory_frame_water_molar_flux_mol_m2_s: float
    laboratory_frame_carrier_molar_flux_mol_m2_s: float
    hexane_diffusive_flux_mol_m2_s: float
    water_diffusive_flux_mol_m2_s: float
    carrier_diffusive_flux_mol_m2_s: float
    component_sum_defect_mol_m2_s: float
    independent_carrier_balance_residual_mol_s: float
    gas_density_kg_m3: float
    superficial_velocity_m_s: float
    layer_pressure_drop_pa: float
    blowing_corrected_heat_factor: float
    wall_to_gas_w: float
    wall_to_interface_w: float
    solid_to_interface_w: float
    interface_convection_w: float
    local_mass_entropy_generation_w_m2_k: float
    local_heat_entropy_generation_w_m2_k: float
    local_total_entropy_generation_w_m2_k: float

    unknowns: ClassVar[tuple[str, ...]] = POSITIVE_CARRIER_UNKNOWN_NAMES
    rows: ClassVar[tuple[str, ...]] = POSITIVE_CARRIER_ROW_NAMES
    analytic_jacobian: ClassVar[bool] = True
    laboratory_frame_components_consumed: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def _law1_saturation_pressure(
    law: cc.ClausiusClapeyronSaturation, temperature: _Law1Dual
) -> _Law1Dual:
    exponent = (temperature - law.reference_temperature_k) * law.slope_per_k
    return law.reference_pressure_pa * _law1_exp(exponent)


def _entropy_component(
    *, diffusive_flux_mol_m2_s: float, bulk_fraction: float, interface_fraction: float
) -> float:
    if bulk_fraction == interface_fraction:
        return 0.0
    if bulk_fraction <= 0.0 or interface_fraction <= 0.0:
        raise EngineeringStepRejected(
            "the finite local entropy diagnostic requires positive endpoint fractions"
        )
    return (
        -UNIVERSAL_GAS_CONSTANT_J_MOL_K
        * diffusive_flux_mol_m2_s
        * (math.log(bulk_fraction) - math.log(interface_fraction))
    )


def _evaluate_positive_carrier(
    inputs: PositiveCarrierEngineeringCellInputs, vector: tuple[float, ...]
) -> PositiveCarrierCellEvaluation:
    _require_positive_carrier_state_domain(vector, error_type=EngineeringStepRejected)
    variables = tuple(_Law1Dual.variable(value, index) for index, value in enumerate(vector))
    t_g, y_h_bulk, y_w_bulk, pressure, t_interface, s_h, s_w = variables
    zero = _Law1Dual(0.0)

    properties = inputs.properties
    boundary = inputs.gas_boundary
    transfer = inputs.transfer
    branch = inputs.solid.interface_branch

    s_total = s_h + s_w
    area = inputs.geometry.gas_side_reference_area_m2 * (
        (1.0 - s_total) ** inputs.film_area_law.exponent.value
    )
    if area.value <= 0.0:
        raise EngineeringStepRejected("the Law-1 gas-side active area is not positive")

    y_h_interface = (
        _law1_saturation_pressure(inputs.equilibrium.hexane, t_interface) / pressure
        if branch.hexane_is_active
        else zero
    )
    y_w_interface = (
        _law1_saturation_pressure(inputs.equilibrium.water, t_interface) / pressure
        if branch.water_is_active
        else zero
    )
    g_bulk = 1.0 - y_h_bulk - y_w_bulk
    g_interface = 1.0 - y_h_interface - y_w_interface
    if g_interface.value <= 0.0:
        raise PositiveCarrierInterfaceCrossoverRefusal(
            "the equilibrium-derived interface reached g_interface <= 0; "
            "the exact-zero topology is refused without fallback"
        )

    peclet = _law1_log(g_bulk) - _law1_log(g_interface)
    bernoulli = _law1_bernoulli(peclet)
    kappa = inputs.transport.positive_carrier_molar_conductance_mol_m2_s
    n_total = kappa * peclet
    j_h = kappa * bernoulli * (y_h_interface - y_h_bulk)
    j_w = kappa * bernoulli * (y_w_interface - y_w_bulk)
    j_g = kappa * bernoulli * (g_interface - g_bulk)
    n_h = y_h_interface * n_total + j_h
    n_w = y_w_interface * n_total + j_w
    n_g = g_interface * n_total + j_g

    inlet_flow = boundary.inlet_molar_flow_mol_s
    outlet_flow = inlet_flow + area * n_total
    if outlet_flow.value <= 0.0:
        raise EngineeringStepRejected("the Law-1 trial outlet molar flow is not positive")

    molar_mass_bulk = (
        y_h_bulk * properties.hexane_molar_mass_kg_mol
        + y_w_bulk * properties.water_molar_mass_kg_mol
        + g_bulk * properties.inert_molar_mass_kg_mol
    )
    density_bulk = pressure * molar_mass_bulk / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * t_g)
    superficial_velocity = (
        outlet_flow * molar_mass_bulk / (density_bulk * inputs.geometry.cross_section_m2)
    )
    if density_bulk.value <= 0.0 or superficial_velocity.value < 0.0:
        raise EngineeringStepRejected("a Law-1 gas density or superficial velocity is invalid")
    try:
        pressure_drop = cc.layer_pressure_drop(
            series=inputs.hydraulics.series,
            superficial_velocity_m_s=superficial_velocity.value,
            gas_density_kg_m3=density_bulk.value,
        )
    except cc.CellClosureError as error:
        raise EngineeringStepRejected(f"Law-1 pressure-loss closure refused: {error}") from error
    pressure_drop_dual = (
        _Law1Dual(pressure_drop.drop_pa)
        + (superficial_velocity - superficial_velocity.value)
        * pressure_drop.d_drop_d_superficial_velocity
        + (density_bulk - density_bulk.value) * pressure_drop.d_drop_d_gas_density
    )

    heat_capacity = properties.molar_heat_capacity_j_mol_k
    datum = properties.energy_datum_temperature_k
    hexane_latent = properties.hexane_latent_heat_j_mol
    water_latent = properties.water_latent_heat_j_mol
    h_gas = heat_capacity * (t_g - datum) + y_h_bulk * hexane_latent + y_w_bulk * water_latent
    h_inlet = (
        heat_capacity * (boundary.inlet_temperature_k - datum)
        + boundary.inlet_hexane_mole_fraction * hexane_latent
        + boundary.inlet_water_mole_fraction * water_latent
    )
    h_hexane = heat_capacity * (t_interface - datum) + hexane_latent
    h_water = heat_capacity * (t_interface - datum) + water_latent
    h_carrier = heat_capacity * (t_interface - datum)
    heat_coefficient = transfer.gas_side_heat_coefficient_w_m2_k
    ackermann = n_total * heat_capacity / heat_coefficient
    heat_kernel = _law1_bernoulli(ackermann)
    gas_to_interface = heat_coefficient * area * (t_g - t_interface) * heat_kernel
    wall_to_gas = transfer.wall_to_gas_ua_w_k * (inputs.wall_temperature_k - t_g)
    wall_to_interface = transfer.wall_to_interface_ua_w_k * (
        inputs.wall_temperature_k - t_interface
    )
    solid_to_interface = transfer.solid_to_interface_ua_w_k * (
        inputs.solid.solid_temperature_k - t_interface
    )

    declared_hexane, declared_water = declared_external_saturations(inputs)
    rows = (
        outlet_flow * y_h_bulk - inlet_flow * boundary.inlet_hexane_mole_fraction - area * n_h,
        outlet_flow * y_w_bulk - inlet_flow * boundary.inlet_water_mole_fraction - area * n_w,
        # F1 SIGN FIX (2026-09-02) - the Law-1 twin of the Law-2 correction
        # above, and the same ``+ gas_to_interface`` on the fourth line of
        # the gas energy row below.  Same construction, same defect, same
        # one-character cure: q_gI leaves the gas (a sink in an
        # out-in-sources row) and enters the interface row as a source,
        # exactly once, equal and opposite.  The pre-fix form fabricated
        # 2*q_gI = -1,731,665.24 W on this arm's shipped fixture.
        outlet_flow * h_gas
        - inlet_flow * h_inlet
        - area * (n_h * h_hexane + n_w * h_water + n_g * h_carrier)
        + gas_to_interface
        - wall_to_gas,
        pressure - boundary.downstream_boundary_pressure_pa - pressure_drop_dual,
        gas_to_interface
        + wall_to_interface
        + solid_to_interface
        - area * (n_h * hexane_latent + n_w * water_latent),
        s_h - declared_hexane,
        s_w - declared_water,
    )
    residual = tuple(row.value for row in rows)
    jacobian = tuple(row.partials for row in rows)
    energy_scale = inlet_flow * heat_capacity * 100.0
    row_scales = (inlet_flow, inlet_flow, energy_scale, 1.0e4, energy_scale, 1.0, 1.0)
    scaled_norm = max(abs(value) / scale for value, scale in zip(residual, row_scales, strict=True))

    carrier_balance = (
        outlet_flow.value * g_bulk.value
        - inlet_flow
        * (1.0 - boundary.inlet_hexane_mole_fraction - boundary.inlet_water_mole_fraction)
        - area.value * n_g.value
    )
    component_sum_defect = n_h.value + n_w.value + n_g.value - n_total.value
    mass_entropy = sum(
        (
            _entropy_component(
                diffusive_flux_mol_m2_s=j_h.value,
                bulk_fraction=y_h_bulk.value,
                interface_fraction=y_h_interface.value,
            ),
            _entropy_component(
                diffusive_flux_mol_m2_s=j_w.value,
                bulk_fraction=y_w_bulk.value,
                interface_fraction=y_w_interface.value,
            ),
            _entropy_component(
                diffusive_flux_mol_m2_s=j_g.value,
                bulk_fraction=g_bulk.value,
                interface_fraction=g_interface.value,
            ),
        )
    )
    heat_flux_per_area = heat_coefficient * (t_g.value - t_interface.value) * heat_kernel.value
    heat_entropy = heat_flux_per_area * (1.0 / t_interface.value - 1.0 / t_g.value)
    total_entropy = mass_entropy + heat_entropy
    if mass_entropy < 0.0 or heat_entropy < 0.0 or total_entropy < 0.0:
        raise EngineeringStepRejected("the local Law-1 entropy diagnostic became negative")

    return PositiveCarrierCellEvaluation(
        residual=residual,
        jacobian=jacobian,
        row_scales=row_scales,
        scaled_residual_norm=scaled_norm,
        gas_side_active_area_m2=area.value,
        outlet_molar_flow_mol_s=outlet_flow.value,
        interface_hexane_mole_fraction=y_h_interface.value,
        interface_water_mole_fraction=y_w_interface.value,
        bulk_carrier_mole_fraction=g_bulk.value,
        interface_carrier_mole_fraction=g_interface.value,
        peclet=peclet.value,
        bernoulli_factor=bernoulli.value,
        total_molar_flux_mol_m2_s=n_total.value,
        laboratory_frame_hexane_molar_flux_mol_m2_s=n_h.value,
        laboratory_frame_water_molar_flux_mol_m2_s=n_w.value,
        laboratory_frame_carrier_molar_flux_mol_m2_s=n_g.value,
        hexane_diffusive_flux_mol_m2_s=j_h.value,
        water_diffusive_flux_mol_m2_s=j_w.value,
        carrier_diffusive_flux_mol_m2_s=j_g.value,
        component_sum_defect_mol_m2_s=component_sum_defect,
        independent_carrier_balance_residual_mol_s=carrier_balance,
        gas_density_kg_m3=density_bulk.value,
        superficial_velocity_m_s=superficial_velocity.value,
        layer_pressure_drop_pa=pressure_drop.drop_pa,
        blowing_corrected_heat_factor=heat_kernel.value,
        wall_to_gas_w=wall_to_gas.value,
        wall_to_interface_w=wall_to_interface.value,
        solid_to_interface_w=solid_to_interface.value,
        interface_convection_w=gas_to_interface.value,
        local_mass_entropy_generation_w_m2_k=mass_entropy,
        local_heat_entropy_generation_w_m2_k=heat_entropy,
        local_total_entropy_generation_w_m2_k=total_entropy,
    )


def assemble_positive_carrier_fast_block(
    inputs: PositiveCarrierEngineeringCellInputs, state: PositiveCarrierCellFieldState
) -> PositiveCarrierCellEvaluation:
    if type(inputs) is not PositiveCarrierEngineeringCellInputs:
        raise EngineeringConfigurationError(
            "inputs must be exact PositiveCarrierEngineeringCellInputs"
        )
    if type(state) is not PositiveCarrierCellFieldState:
        raise EngineeringConfigurationError("state must be an exact PositiveCarrierCellFieldState")
    if (
        dispatch_engineering_carrier_topology(inputs)
        is not EngineeringCarrierLaw.POSITIVE_CARRIER_LAW1
    ):
        raise CarrierTopologyConfigurationError("pre-solve topology dispatch did not select Law 1")
    return _evaluate_positive_carrier(inputs, state.as_vector())


def nominal_positive_carrier_seed(
    inputs: PositiveCarrierEngineeringCellInputs,
) -> PositiveCarrierCellFieldState:
    """Return a deterministic declared-domain seed, never a fitted state."""

    if type(inputs) is not PositiveCarrierEngineeringCellInputs:
        raise EngineeringConfigurationError(
            "inputs must be exact PositiveCarrierEngineeringCellInputs"
        )
    s_h, s_w = declared_external_saturations(inputs)
    pressure = inputs.gas_boundary.downstream_boundary_pressure_pa + 3_000.0
    interface_temperature = 335.0
    y_h_interface = (
        inputs.equilibrium.hexane.saturation_pressure_pa(interface_temperature) / pressure
        if inputs.solid.interface_branch.hexane_is_active
        else 0.0
    )
    y_w_interface = (
        inputs.equilibrium.water.saturation_pressure_pa(interface_temperature) / pressure
        if inputs.solid.interface_branch.water_is_active
        else 0.0
    )
    return PositiveCarrierCellFieldState(
        gas_temperature_k=336.0,
        gas_hexane_mole_fraction=(inputs.gas_boundary.inlet_hexane_mole_fraction + y_h_interface)
        / 2.0,
        gas_water_mole_fraction=(inputs.gas_boundary.inlet_water_mole_fraction + y_w_interface)
        / 2.0,
        layer_pressure_pa=pressure,
        interface_temperature_k=interface_temperature,
        external_hexane_saturation=s_h,
        external_water_saturation=s_w,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierCellSolve:
    state: PositiveCarrierCellFieldState
    evaluation: PositiveCarrierCellEvaluation
    iterations: int
    tolerance: float
    inputs_digest: str
    accepted_iterates: tuple[PositiveCarrierCellFieldState, ...]
    accepted_scaled_residual_norms: tuple[float, ...]
    line_search_rejections: int

    analytic_jacobian: ClassVar[bool] = True
    rejected_trials_are_never_states: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def _scaled_positive_carrier_newton_step(
    evaluation: PositiveCarrierCellEvaluation,
) -> tuple[float, ...]:
    row_scale = np.asarray(evaluation.row_scales, dtype=float)
    column_scale = np.asarray(POSITIVE_CARRIER_COLUMN_SCALES, dtype=float)
    residual = np.asarray(evaluation.residual, dtype=float) / row_scale
    jacobian = np.asarray(evaluation.jacobian, dtype=float)
    scaled_jacobian = jacobian * column_scale[np.newaxis, :] / row_scale[:, np.newaxis]
    try:
        scaled_step = np.linalg.solve(scaled_jacobian, -residual)
    except np.linalg.LinAlgError as error:
        raise EngineeringNewtonConvergenceError(
            "the scaled analytic Law-1 Jacobian is singular"
        ) from error
    step = scaled_step * column_scale
    if not np.all(np.isfinite(step)):
        raise EngineeringNewtonConvergenceError("the analytic Law-1 Newton step is not finite")
    return tuple(float(value) for value in step)


def solve_positive_carrier_fast_block(
    inputs: PositiveCarrierEngineeringCellInputs,
    *,
    seed: PositiveCarrierCellFieldState,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
    max_line_search_reductions: int = DEFAULT_LINE_SEARCH_REDUCTIONS,
) -> PositiveCarrierCellSolve:
    """Solve the bounded 7x7 Law-1 block by scaled analytic damped Newton."""

    if type(inputs) is not PositiveCarrierEngineeringCellInputs:
        raise EngineeringConfigurationError(
            "inputs must be exact PositiveCarrierEngineeringCellInputs"
        )
    if type(seed) is not PositiveCarrierCellFieldState:
        raise EngineeringConfigurationError("seed must be an exact PositiveCarrierCellFieldState")
    if (
        dispatch_engineering_carrier_topology(inputs)
        is not EngineeringCarrierLaw.POSITIVE_CARRIER_LAW1
    ):
        raise CarrierTopologyConfigurationError("pre-solve topology dispatch did not select Law 1")
    _require_positive("tolerance", tolerance)
    if type(max_iterations) is not int or max_iterations <= 0:
        raise EngineeringConfigurationError("max_iterations must be a positive int")
    if type(max_line_search_reductions) is not int or max_line_search_reductions <= 0:
        raise EngineeringConfigurationError("max_line_search_reductions must be a positive int")

    vector = seed.as_vector()
    evaluation = _evaluate_positive_carrier(inputs, vector)
    states = [seed]
    norms = [evaluation.scaled_residual_norm]
    rejected_trials = 0
    iterations = 0
    while evaluation.scaled_residual_norm > tolerance:
        if iterations >= max_iterations:
            raise EngineeringNewtonConvergenceError(
                f"Law 1 did not converge in {max_iterations} accepted Newton steps; "
                f"scaled residual {evaluation.scaled_residual_norm!r}"
            )
        step = _scaled_positive_carrier_newton_step(evaluation)
        step_factor = 1.0
        accepted = False
        for _ in range(max_line_search_reductions):
            candidate = tuple(
                float(current + step_factor * increment)
                for current, increment in zip(vector, step, strict=True)
            )
            try:
                candidate_evaluation = _evaluate_positive_carrier(inputs, candidate)
            except EngineeringStepRejected:
                rejected_trials += 1
                step_factor *= 0.5
                continue
            armijo_target = (1.0 - 1.0e-4 * step_factor) * evaluation.scaled_residual_norm
            if (
                candidate_evaluation.scaled_residual_norm <= armijo_target
                or candidate_evaluation.scaled_residual_norm <= tolerance
            ):
                vector = candidate
                evaluation = candidate_evaluation
                state = PositiveCarrierCellFieldState.from_vector(vector)
                states.append(state)
                norms.append(evaluation.scaled_residual_norm)
                iterations += 1
                accepted = True
                break
            rejected_trials += 1
            step_factor *= 0.5
        if not accepted:
            raise EngineeringNewtonConvergenceError(
                "the analytic Law-1 Newton direction produced no admissible "
                "residual-decreasing trial; refused without clamping or fallback"
            )

    state = PositiveCarrierCellFieldState.from_vector(vector)
    final_evaluation = _evaluate_positive_carrier(inputs, state.as_vector())
    if final_evaluation.scaled_residual_norm > tolerance:
        raise EngineeringNewtonConvergenceError(
            "the revalidated Law-1 root does not meet the requested tolerance"
        )
    if inputs.require_corroborated_pressure_drop:
        cc.require_corroborated_layer_drop(final_evaluation.layer_pressure_drop_pa)
    return PositiveCarrierCellSolve(
        state=state,
        evaluation=final_evaluation,
        iterations=iterations,
        tolerance=tolerance,
        inputs_digest=inputs.inputs_digest,
        accepted_iterates=tuple(states),
        accepted_scaled_residual_norms=tuple(norms),
        line_search_rejections=rejected_trials,
    )


_ACCEPTANCE_TOKEN = object()
_POSITIVE_CARRIER_ACCEPTANCE_TOKEN = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedEngineeringFastStep:
    """Immutable owner of the two shared transfers at the accepted old wall."""

    inputs_digest: str
    state: BinaryNoInertCellFieldState
    evaluation: BinaryNoInertCellEvaluation
    old_wall_temperature_k: float
    wall_to_gas_w: float
    wall_to_interface_w: float
    macro_step_s: float
    wall_storage_scale: float
    transfer: cc.TransferCoefficients
    wall_parameters: cc.WallNodeParameters
    claims: EngineeringFeasibilityClaimBoundary
    acceptance_token: InitVar[object] = None

    def __post_init__(self, acceptance_token: object) -> None:
        if acceptance_token is not _ACCEPTANCE_TOKEN:
            raise AcceptedStepIntegrityError(
                "AcceptedEngineeringFastStep must be created by accept_engineering_fast_step"
            )
        expected_gas = self.transfer.wall_to_gas_ua_w_k * (
            self.old_wall_temperature_k - self.state.gas_temperature_k
        )
        expected_interface = self.transfer.wall_to_interface_ua_w_k * (
            self.old_wall_temperature_k - self.state.interface_temperature_k
        )
        if self.wall_to_gas_w != expected_gas or self.wall_to_interface_w != expected_interface:
            raise AcceptedStepIntegrityError(
                "accepted wall transfers do not match the accepted state and old wall"
            )
        if (
            self.evaluation.wall_to_gas_w != self.wall_to_gas_w
            or self.evaluation.wall_to_interface_w != self.wall_to_interface_w
        ):
            raise AcceptedStepIntegrityError(
                "accepted wall transfers were substituted after the fast solve"
            )


def _require_acceptable_tolerance(tolerance: float, *, law: str) -> None:
    if type(tolerance) is not float or not math.isfinite(tolerance) or tolerance <= 0.0:
        raise AcceptedStepIntegrityError(
            f"the {law} solve tolerance must be a finite, strictly positive binary64 float"
        )


def _require_finite_evaluation_core(
    evaluation: BinaryNoInertCellEvaluation | PositiveCarrierCellEvaluation,
    *,
    law: str,
) -> None:
    if type(evaluation.residual) is not tuple or any(
        type(value) is not float or not math.isfinite(value) for value in evaluation.residual
    ):
        raise AcceptedStepIntegrityError(f"the {law} evaluation contains a nonfinite residual")
    if type(evaluation.jacobian) is not tuple or any(
        type(row) is not tuple
        or any(type(value) is not float or not math.isfinite(value) for value in row)
        for row in evaluation.jacobian
    ):
        raise AcceptedStepIntegrityError(f"the {law} evaluation contains a nonfinite Jacobian")
    if (
        type(evaluation.scaled_residual_norm) is not float
        or not math.isfinite(evaluation.scaled_residual_norm)
        or evaluation.scaled_residual_norm < 0.0
    ):
        raise AcceptedStepIntegrityError(
            f"the {law} evaluation contains an invalid scaled residual norm"
        )


def _require_finite_law1_entropy(evaluation: PositiveCarrierCellEvaluation) -> None:
    entropy_values = (
        evaluation.local_mass_entropy_generation_w_m2_k,
        evaluation.local_heat_entropy_generation_w_m2_k,
        evaluation.local_total_entropy_generation_w_m2_k,
    )
    if any(type(value) is not float or not math.isfinite(value) for value in entropy_values):
        raise AcceptedStepIntegrityError("the Law-1 evaluation contains nonfinite local entropy")
    if any(value < 0.0 for value in entropy_values):
        raise AcceptedStepIntegrityError("a negative local Law-1 entropy result cannot be accepted")


def accept_engineering_fast_step(
    inputs: EngineeringCellInputs, solve: BinaryNoInertCellSolve
) -> AcceptedEngineeringFastStep:
    if type(inputs) is not EngineeringCellInputs:
        raise EngineeringConfigurationError("inputs must be an exact EngineeringCellInputs")
    if type(solve) is not BinaryNoInertCellSolve:
        raise AcceptedStepIntegrityError("solve must be an exact BinaryNoInertCellSolve")
    if solve.inputs_digest != inputs.inputs_digest:
        raise AcceptedStepIntegrityError("the solve and engineering inputs have different digests")
    _require_acceptable_tolerance(solve.tolerance, law="Law-2")
    # D9-e: the gate admits the ten-unknown dry-shell evaluation BY NAME, as an
    # exact-type membership exactly as it was written - no isinstance, no
    # removal.  Everything below it then runs on the dry-shell evaluation the
    # same way it runs on the binary one: the reconstruction is the same
    # evaluator call on the same state vector, the equality is the dataclass's
    # own, and the norm bound is the same tolerance.
    if type(solve.evaluation) not in ACCEPTED_FAST_BLOCK_EVALUATION_TYPES:
        raise AcceptedStepIntegrityError("the Law-2 solve contains a foreign evaluation type")
    _require_finite_evaluation_core(solve.evaluation, law="Law-2")
    try:
        recomputed = assemble_binary_no_inert_fast_block(inputs, solve.state)
    except EngineeringFeasibilityError as error:
        raise AcceptedStepIntegrityError(
            "the Law-2 evaluation could not be reconstructed from inputs and solve state"
        ) from error
    _require_finite_evaluation_core(recomputed, law="Law-2")
    if solve.evaluation != recomputed:
        raise AcceptedStepIntegrityError(
            "the supplied Law-2 evaluation does not exactly match inputs and solve state"
        )
    if recomputed.scaled_residual_norm > solve.tolerance:
        raise AcceptedStepIntegrityError("a non-converged fast block cannot be accepted")
    return AcceptedEngineeringFastStep(
        inputs_digest=inputs.inputs_digest,
        state=solve.state,
        evaluation=recomputed,
        old_wall_temperature_k=inputs.wall_temperature_k,
        wall_to_gas_w=recomputed.wall_to_gas_w,
        wall_to_interface_w=recomputed.wall_to_interface_w,
        macro_step_s=inputs.macro_step_s,
        wall_storage_scale=inputs.wall_storage_scale,
        transfer=inputs.transfer,
        wall_parameters=inputs.wall,
        claims=ENGINEERING_FEASIBILITY_CLAIMS,
        acceptance_token=_ACCEPTANCE_TOKEN,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedPositiveCarrierEngineeringFastStep:
    """Immutable owner of accepted Law-1 state and shared old-wall transfers."""

    inputs_digest: str
    state: PositiveCarrierCellFieldState
    evaluation: PositiveCarrierCellEvaluation
    old_wall_temperature_k: float
    wall_to_gas_w: float
    wall_to_interface_w: float
    macro_step_s: float
    wall_storage_scale: float
    transfer: cc.TransferCoefficients
    wall_parameters: cc.WallNodeParameters
    claims: EngineeringFeasibilityClaimBoundary
    acceptance_token: InitVar[object] = None

    def __post_init__(self, acceptance_token: object) -> None:
        if acceptance_token is not _POSITIVE_CARRIER_ACCEPTANCE_TOKEN:
            raise AcceptedStepIntegrityError(
                "AcceptedPositiveCarrierEngineeringFastStep must be created by its Law-1 acceptor"
            )
        expected_gas = self.transfer.wall_to_gas_ua_w_k * (
            self.old_wall_temperature_k - self.state.gas_temperature_k
        )
        expected_interface = self.transfer.wall_to_interface_ua_w_k * (
            self.old_wall_temperature_k - self.state.interface_temperature_k
        )
        if self.wall_to_gas_w != expected_gas or self.wall_to_interface_w != expected_interface:
            raise AcceptedStepIntegrityError(
                "accepted Law-1 wall transfers do not match the accepted state and old wall"
            )
        if (
            self.evaluation.wall_to_gas_w != self.wall_to_gas_w
            or self.evaluation.wall_to_interface_w != self.wall_to_interface_w
        ):
            raise AcceptedStepIntegrityError(
                "accepted Law-1 wall transfers were substituted after the fast solve"
            )


def accept_positive_carrier_engineering_fast_step(
    inputs: PositiveCarrierEngineeringCellInputs,
    solve: PositiveCarrierCellSolve,
) -> AcceptedPositiveCarrierEngineeringFastStep:
    if type(inputs) is not PositiveCarrierEngineeringCellInputs:
        raise EngineeringConfigurationError(
            "inputs must be exact PositiveCarrierEngineeringCellInputs"
        )
    if type(solve) is not PositiveCarrierCellSolve:
        raise AcceptedStepIntegrityError("solve must be an exact PositiveCarrierCellSolve")
    if solve.inputs_digest != inputs.inputs_digest:
        raise AcceptedStepIntegrityError("the Law-1 solve and inputs have different digests")
    _require_acceptable_tolerance(solve.tolerance, law="Law-1")
    if type(solve.evaluation) is not PositiveCarrierCellEvaluation:
        raise AcceptedStepIntegrityError("the Law-1 solve contains a foreign evaluation type")
    _require_finite_evaluation_core(solve.evaluation, law="Law-1")
    _require_finite_law1_entropy(solve.evaluation)
    try:
        recomputed = assemble_positive_carrier_fast_block(inputs, solve.state)
    except EngineeringFeasibilityError as error:
        raise AcceptedStepIntegrityError(
            "the Law-1 evaluation could not be reconstructed from inputs and solve state"
        ) from error
    _require_finite_evaluation_core(recomputed, law="Law-1")
    _require_finite_law1_entropy(recomputed)
    if solve.evaluation != recomputed:
        raise AcceptedStepIntegrityError(
            "the supplied Law-1 evaluation does not exactly match inputs and solve state"
        )
    if recomputed.scaled_residual_norm > solve.tolerance:
        raise AcceptedStepIntegrityError("a non-converged Law-1 fast block cannot be accepted")
    return AcceptedPositiveCarrierEngineeringFastStep(
        inputs_digest=inputs.inputs_digest,
        state=solve.state,
        evaluation=recomputed,
        old_wall_temperature_k=inputs.wall_temperature_k,
        wall_to_gas_w=recomputed.wall_to_gas_w,
        wall_to_interface_w=recomputed.wall_to_interface_w,
        macro_step_s=inputs.macro_step_s,
        wall_storage_scale=inputs.wall_storage_scale,
        transfer=inputs.transfer,
        wall_parameters=inputs.wall,
        claims=ENGINEERING_FEASIBILITY_CLAIMS,
        acceptance_token=_POSITIVE_CARRIER_ACCEPTANCE_TOKEN,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringWallNodeSolve:
    wall_temperature_k: float
    old_wall_temperature_k: float
    steam_side_w: float
    to_ambient_w: float
    shared_to_gas_w: float
    shared_to_interface_w: float
    storage_w: float
    residual_w: float
    row_derivative_w_k: float

    explicit_shared_old_wall: ClassVar[bool] = True
    coupling_uas_excluded_from_implicit_denominator: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def relative_residual(self) -> float:
        scale = max(
            1.0,
            abs(self.storage_w),
            abs(self.steam_side_w),
            abs(self.to_ambient_w),
            abs(self.shared_to_gas_w),
            abs(self.shared_to_interface_w),
        )
        return abs(self.residual_w) / scale


def _advance_shared_old_wall(
    accepted: AcceptedEngineeringFastStep | AcceptedPositiveCarrierEngineeringFastStep,
    *,
    wall: cc.WallNodeParameters,
) -> EngineeringWallNodeSolve:
    """One exact shared-old-wall algebra for the two separately accepted types."""

    storage_conductance = (
        accepted.wall_storage_scale * wall.thermal_capacity_j_k / accepted.macro_step_s
    )
    denominator = storage_conductance + wall.steam_side_ua_w_k + wall.ambient_ua_w_k
    if denominator == 0.0:
        raise EngineeringConfigurationError(
            "the shared-old-wall row has no storage, steam, or ambient conductance"
        )
    numerator = (
        storage_conductance * accepted.old_wall_temperature_k
        + wall.steam_side_ua_w_k * wall.steam_temperature_k
        + wall.ambient_ua_w_k * wall.ambient_temperature_k
        - accepted.wall_to_gas_w
        - accepted.wall_to_interface_w
    )
    temperature = numerator / denominator
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise EngineeringStepRejected("the shared-old-wall update produced an invalid temperature")
    steam_side = wall.steam_side_ua_w_k * (wall.steam_temperature_k - temperature)
    ambient = wall.ambient_ua_w_k * (temperature - wall.ambient_temperature_k)
    storage = storage_conductance * (temperature - accepted.old_wall_temperature_k)
    residual = (
        storage - steam_side + ambient + accepted.wall_to_gas_w + accepted.wall_to_interface_w
    )
    return EngineeringWallNodeSolve(
        wall_temperature_k=temperature,
        old_wall_temperature_k=accepted.old_wall_temperature_k,
        steam_side_w=steam_side,
        to_ambient_w=ambient,
        shared_to_gas_w=accepted.wall_to_gas_w,
        shared_to_interface_w=accepted.wall_to_interface_w,
        storage_w=storage,
        residual_w=residual,
        row_derivative_w_k=denominator,
    )


def advance_engineering_wall_node(
    accepted: AcceptedEngineeringFastStep, *, wall: cc.WallNodeParameters
) -> EngineeringWallNodeSolve:
    """Advance from accepted Law-2 old-wall transfers without recomputing either one."""

    if type(accepted) is not AcceptedEngineeringFastStep:
        raise AcceptedStepIntegrityError("accepted must be an AcceptedEngineeringFastStep")
    if type(wall) is not cc.WallNodeParameters:
        raise EngineeringConfigurationError("wall must be an exact WallNodeParameters")
    if wall != accepted.wall_parameters:
        raise AcceptedStepIntegrityError("the wall parameters differ from the accepted fast step")
    return _advance_shared_old_wall(accepted, wall=wall)


def advance_positive_carrier_engineering_wall_node(
    accepted: AcceptedPositiveCarrierEngineeringFastStep,
    *,
    wall: cc.WallNodeParameters,
) -> EngineeringWallNodeSolve:
    """Advance from accepted Law-1 old-wall transfers without recomputing either one."""

    if type(accepted) is not AcceptedPositiveCarrierEngineeringFastStep:
        raise AcceptedStepIntegrityError(
            "accepted must be an AcceptedPositiveCarrierEngineeringFastStep"
        )
    if type(wall) is not cc.WallNodeParameters:
        raise EngineeringConfigurationError("wall must be an exact WallNodeParameters")
    if wall != accepted.wall_parameters:
        raise AcceptedStepIntegrityError("the wall parameters differ from the accepted fast step")
    return _advance_shared_old_wall(accepted, wall=wall)


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringBoundaryLedger:
    independent_water_balance_residual_mol_s: float
    hexane_fugacity_residual_pa: float
    water_fugacity_residual_pa: float
    wall_to_gas_shared_cancellation_w: float
    wall_to_interface_shared_cancellation_w: float
    wall_row_residual_w: float
    wall_row_relative_residual: float
    claims: EngineeringFeasibilityClaimBoundary


def engineering_boundary_ledger(
    accepted: AcceptedEngineeringFastStep, wall: EngineeringWallNodeSolve
) -> EngineeringBoundaryLedger:
    if type(accepted) is not AcceptedEngineeringFastStep:
        raise AcceptedStepIntegrityError("accepted must be an AcceptedEngineeringFastStep")
    if type(wall) is not EngineeringWallNodeSolve:
        raise EngineeringConfigurationError("wall must be an EngineeringWallNodeSolve")
    return EngineeringBoundaryLedger(
        independent_water_balance_residual_mol_s=(
            accepted.evaluation.independent_water_balance_residual_mol_s
        ),
        hexane_fugacity_residual_pa=accepted.evaluation.residual[6],
        water_fugacity_residual_pa=accepted.evaluation.residual[7],
        wall_to_gas_shared_cancellation_w=(accepted.wall_to_gas_w - wall.shared_to_gas_w),
        wall_to_interface_shared_cancellation_w=(
            accepted.wall_to_interface_w - wall.shared_to_interface_w
        ),
        wall_row_residual_w=wall.residual_w,
        wall_row_relative_residual=wall.relative_residual,
        claims=accepted.claims,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringCellMacroStep:
    fast: BinaryNoInertCellSolve
    accepted: AcceptedEngineeringFastStep
    wall: EngineeringWallNodeSolve
    ledger: EngineeringBoundaryLedger
    claims: EngineeringFeasibilityClaimBoundary


def march_binary_no_inert_engineering_step(
    inputs: EngineeringCellInputs,
    *,
    seed: BinaryNoInertCellFieldState,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
) -> EngineeringCellMacroStep:
    fast = solve_binary_no_inert_fast_block(
        inputs,
        seed=seed,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
    accepted = accept_engineering_fast_step(inputs, fast)
    wall = advance_engineering_wall_node(accepted, wall=inputs.wall)
    ledger = engineering_boundary_ledger(accepted, wall)
    return EngineeringCellMacroStep(
        fast=fast,
        accepted=accepted,
        wall=wall,
        ledger=ledger,
        claims=ENGINEERING_FEASIBILITY_CLAIMS,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierEngineeringBoundaryLedger:
    gas_hexane_balance_residual_mol_s: float
    gas_water_balance_residual_mol_s: float
    independent_carrier_balance_residual_mol_s: float
    gas_energy_balance_residual_w: float
    interface_energy_balance_residual_w: float
    laboratory_frame_total_molar_flux_mol_m2_s: float
    laboratory_frame_hexane_molar_flux_mol_m2_s: float
    laboratory_frame_water_molar_flux_mol_m2_s: float
    laboratory_frame_carrier_molar_flux_mol_m2_s: float
    component_sum_defect_mol_m2_s: float
    local_mass_entropy_generation_w_m2_k: float
    local_heat_entropy_generation_w_m2_k: float
    local_total_entropy_generation_w_m2_k: float
    wall_to_gas_shared_cancellation_w: float
    wall_to_interface_shared_cancellation_w: float
    wall_row_residual_w: float
    wall_row_relative_residual: float
    claims: EngineeringFeasibilityClaimBoundary

    laboratory_frame_components_consumed: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def positive_carrier_engineering_boundary_ledger(
    accepted: AcceptedPositiveCarrierEngineeringFastStep,
    wall: EngineeringWallNodeSolve,
) -> PositiveCarrierEngineeringBoundaryLedger:
    if type(accepted) is not AcceptedPositiveCarrierEngineeringFastStep:
        raise AcceptedStepIntegrityError(
            "accepted must be an AcceptedPositiveCarrierEngineeringFastStep"
        )
    if type(wall) is not EngineeringWallNodeSolve:
        raise EngineeringConfigurationError("wall must be an EngineeringWallNodeSolve")
    evaluation = accepted.evaluation
    return PositiveCarrierEngineeringBoundaryLedger(
        gas_hexane_balance_residual_mol_s=evaluation.residual[0],
        gas_water_balance_residual_mol_s=evaluation.residual[1],
        independent_carrier_balance_residual_mol_s=(
            evaluation.independent_carrier_balance_residual_mol_s
        ),
        gas_energy_balance_residual_w=evaluation.residual[2],
        interface_energy_balance_residual_w=evaluation.residual[4],
        laboratory_frame_total_molar_flux_mol_m2_s=(evaluation.total_molar_flux_mol_m2_s),
        laboratory_frame_hexane_molar_flux_mol_m2_s=(
            evaluation.laboratory_frame_hexane_molar_flux_mol_m2_s
        ),
        laboratory_frame_water_molar_flux_mol_m2_s=(
            evaluation.laboratory_frame_water_molar_flux_mol_m2_s
        ),
        laboratory_frame_carrier_molar_flux_mol_m2_s=(
            evaluation.laboratory_frame_carrier_molar_flux_mol_m2_s
        ),
        component_sum_defect_mol_m2_s=evaluation.component_sum_defect_mol_m2_s,
        local_mass_entropy_generation_w_m2_k=(evaluation.local_mass_entropy_generation_w_m2_k),
        local_heat_entropy_generation_w_m2_k=(evaluation.local_heat_entropy_generation_w_m2_k),
        local_total_entropy_generation_w_m2_k=(evaluation.local_total_entropy_generation_w_m2_k),
        wall_to_gas_shared_cancellation_w=(accepted.wall_to_gas_w - wall.shared_to_gas_w),
        wall_to_interface_shared_cancellation_w=(
            accepted.wall_to_interface_w - wall.shared_to_interface_w
        ),
        wall_row_residual_w=wall.residual_w,
        wall_row_relative_residual=wall.relative_residual,
        claims=accepted.claims,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class PositiveCarrierEngineeringCellMacroStep:
    fast: PositiveCarrierCellSolve
    accepted: AcceptedPositiveCarrierEngineeringFastStep
    wall: EngineeringWallNodeSolve
    ledger: PositiveCarrierEngineeringBoundaryLedger
    claims: EngineeringFeasibilityClaimBoundary


def march_positive_carrier_engineering_step(
    inputs: PositiveCarrierEngineeringCellInputs,
    *,
    seed: PositiveCarrierCellFieldState,
    tolerance: float = DEFAULT_NEWTON_TOLERANCE,
    max_iterations: int = DEFAULT_NEWTON_MAX_ITERATIONS,
) -> PositiveCarrierEngineeringCellMacroStep:
    if (
        dispatch_engineering_carrier_topology(inputs)
        is not EngineeringCarrierLaw.POSITIVE_CARRIER_LAW1
    ):
        raise CarrierTopologyConfigurationError("pre-solve topology dispatch did not select Law 1")
    fast = solve_positive_carrier_fast_block(
        inputs,
        seed=seed,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
    accepted = accept_positive_carrier_engineering_fast_step(inputs, fast)
    wall = advance_positive_carrier_engineering_wall_node(accepted, wall=inputs.wall)
    ledger = positive_carrier_engineering_boundary_ledger(accepted, wall)
    return PositiveCarrierEngineeringCellMacroStep(
        fast=fast,
        accepted=accepted,
        wall=wall,
        ledger=ledger,
        claims=ENGINEERING_FEASIBILITY_CLAIMS,
    )


__all__ = [
    "AcceptedPositiveCarrierEngineeringFastStep",
    "AcceptedEngineeringFastStep",
    "AcceptedStepIntegrityError",
    "ARBITRARY_K_PRESSURE_ADMISSION_EXTENDS_F4",
    "ARBITRARY_K_PRESSURE_ADMISSION_EXTENDS_SATURATION_CERT01",
    "ARBITRARY_K_PRESSURE_ADMISSION_PHYSICALLY_QUALIFYING",
    "BINARY_NO_INERT_COLUMN_SCALES",
    "BINARY_NO_INERT_ROW_NAMES",
    "BINARY_NO_INERT_UNKNOWN_NAMES",
    "ACCEPTED_FAST_BLOCK_EVALUATION_TYPES",
    "ACCEPTED_FAST_BLOCK_STATE_TYPES",
    "DRY_SHELL_COLUMN_SCALES",
    "D9E_CARRIED_SEED_DECLARED_ASSUMPTIONS",
    "D9E_CARRIED_SEED_ROUTING_LAW_ID",
    "DRY_SHELL_HINDRANCE_CONTINUATION_START",
    "DRY_SHELL_HINDRANCE_CONTINUATION_STEPS",
    "DRY_SHELL_MEAL_CONDUCTIVITY_SENSITIVITY_W_M_K",
    "DRY_SHELL_MEAL_CONDUCTIVITY_W_M_K",
    "DRY_SHELL_PARTICLE_RADIUS_M",
    "DRY_SHELL_ROW_NAMES",
    "DRY_SHELL_SURFACE_ACTIVITY_DOMAIN",
    "DRY_SHELL_UNKNOWN_NAMES",
    "DRY_SHELL_WATER_ARM_ACTIVE",
    "DRY_SHELL_WATER_ARM_BLOCKED",
    "DryShellCellEvaluation",
    "DryShellCellFieldState",
    "NO_SHELL_FRONT_FRACTION_CEILING",
    "binary_seed_from_dry_shell",
    "dry_shell_seed_from_binary",
    "dry_shell_surface_activity_seed",
    "nominal_dry_shell_seed",
    "presents_dry_shell_receding_front_interface",
    "solve_dry_shell_receding_front_fast_block",
    "BinaryNoInertCellEvaluation",
    "BinaryNoInertCellFieldState",
    "BinaryNoInertCellSolve",
    "CANONICAL_HEXANE_MOLAR_MASS_KG_MOL",
    "CANONICAL_WATER_MOLAR_MASS_KG_MOL",
    "CHILTON_COLBURN_ANALOGY_EXPONENT",
    "COLETTO_B7_NUSSELT_PREFACTOR",
    "COLETTO_B7_PRANDTL_EXPONENT",
    "COLETTO_B7_PRIMARY_SOURCE_CITATION_STATUS",
    "COLETTO_B7_REYNOLDS_EXPONENT",
    "CarrierTopology",
    "CarrierTopologyConfigurationError",
    "CoreSuppliedClosureAuthority",
    "CoreSuppliedHexaneClosure",
    "CoreSuppliedWaterClosure",
    "CoreSuppliedWaterClosureAuthority",
    "EngineeringCarrierLaw",
    "ENGINEERING_FEASIBILITY_CLAIMS",
    "EngineeringBoundaryLedger",
    "EngineeringCellInputs",
    "EngineeringCellMacroStep",
    "EngineeringConfigurationError",
    "EngineeringFeasibilityClaimBoundary",
    "EngineeringFeasibilityError",
    "EngineeringGasBoundary",
    "EngineeringNewtonConvergenceError",
    "EngineeringStepRejected",
    "EngineeringWallNodeSolve",
    "ENGINEERING_ARBITRARY_K_SERIAL_LAW2_PRESSURE_ADMISSION_PA",
    "ENGINEERING_SIX_LAYER_SERIAL_LAW2_PRESSURE_ADMISSION_PA",
    "FSG_DHW_AUTHORITY_DECISION_ID",
    "FSG_DHW_AUTHORITY_ID",
    "FSG_DHW_DIFFUSIVITY_AUTHORITY",
    "FSG_DIFFUSIVITY_PREFACTOR",
    "FSG_HEXANE_MOLAR_MASS_KG_MOL",
    "FSG_TEMPERATURE_EXPONENT",
    "FSGDhwDiffusivityAuthority",
    "FROZEN_SINGLE_TRAY_LAW2_PRESSURE_CERTIFICATION_BOX_PA",
    "LAW2_FILM_CORRELATION_AUTHORITY",
    "LAW2_FILM_CORRELATION_AUTHORITY_ID",
    "LAW2_FILM_CORRELATION_DECISION_ID",
    "Law2CoefficientModel",
    "Law2FilmCorrelationAuthority",
    "Law2TransportInputs",
    "Law1ConductanceAuthority",
    "Law1TransportInputs",
    "POSITIVE_CARRIER_COLUMN_SCALES",
    "POSITIVE_CARRIER_ROW_NAMES",
    "POSITIVE_CARRIER_UNKNOWN_NAMES",
    "PositiveCarrierCellEvaluation",
    "PositiveCarrierCellFieldState",
    "PositiveCarrierCellSolve",
    "PositiveCarrierEngineeringBoundaryLedger",
    "PositiveCarrierEngineeringCellInputs",
    "PositiveCarrierEngineeringCellMacroStep",
    "PositiveCarrierFluxVector",
    "PositiveCarrierInterfaceCrossoverRefusal",
    "PositiveCarrierTopologyRefusal",
    "UnsupportedInventoryBranchRefusal",
    "UnsupportedStatefulClosureRefusal",
    "accept_engineering_fast_step",
    "accept_positive_carrier_engineering_fast_step",
    "advance_engineering_wall_node",
    "advance_positive_carrier_engineering_wall_node",
    "assemble_binary_no_inert_fast_block",
    "assemble_positive_carrier_fast_block",
    "bernoulli_pair",
    "declared_external_saturations",
    "dispatch_engineering_carrier_topology",
    "engineering_boundary_ledger",
    "march_binary_no_inert_engineering_step",
    "march_positive_carrier_engineering_step",
    "nominal_binary_no_inert_seed",
    "nominal_positive_carrier_seed",
    "positive_carrier_engineering_boundary_ledger",
    "positive_carrier_fluxes",
    "solve_binary_no_inert_fast_block",
    "solve_positive_carrier_fast_block",
]
