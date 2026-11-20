"""Numeric datum bridge between declared Law-2 gas and native packet energy.

The Law-2 engineering cell uses a declared gas enthalpy convention

``h_i(T) = cp_declared * (T - T_datum) + L_i``.

The tag-1 packet instead uses the native Span--Wagner n-hexane and IAPWS-95
water caloric datums.  This module aligns only the arbitrary component zeros:
at one explicit dilute-vapor ``(T_ref, P_ref)``, it derives the additive molar
offset that makes each native vapor enthalpy equal the declared Law-2 value.
The same constant is then applied to every mole of that component in the
packet, irrespective of phase.  Therefore a declared species-offset shift is
carried by both subsystems and cancels from phase-transfer physics.

No heat capacity, latent heat, phase correction, or temperature-dependent fit
is introduced.  The declared Law-2 and native-property caloric curves remain
different away from the reference point, which is one reason this bridge and
its consumers remain engineering-only and nonqualifying.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass
from typing import ClassVar

from . import cell_closure as cc
from . import cell_native_caloric as nc
from .cell_native_caloric import NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID
from .props import hexane, water


DATUM_ADAPTER_DIGEST_DOMAIN = "GT-PS-2/component-energy-datum-adapter/v1"
NATIVE_ZERO_GAUGE_DATUM_ADAPTER_DIGEST_DOMAIN = (
    "GT-PS-2/native-zero-gauge-component-energy-datum-adapter/v1"
)


class ComponentEnergyDatumError(ValueError):
    """Base typed refusal for numeric component-datum work."""


class ComponentEnergyDatumConfigurationError(ComponentEnergyDatumError):
    """The declared reference state or numeric binding is malformed."""


class ComponentEnergyDatumDriftError(ComponentEnergyDatumError):
    """A consumer changed a pinned component, datum, or gas convention."""


def _float_bits(value: float) -> bytes:
    return struct.pack(">d", value)


def _same_float(left: float, right: float) -> bool:
    return _float_bits(left) == _float_bits(right)


def _require_float(name: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise ComponentEnergyDatumConfigurationError(f"{name} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise ComponentEnergyDatumConfigurationError(f"{name} must be strictly positive")


def _require_text(name: str, value: str) -> None:
    if type(value) is not str or not value.strip():
        raise ComponentEnergyDatumConfigurationError(f"{name} must be a nonblank exact str")


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class NumericComponentEnergyDatumAdapter:
    """Exact numeric zero alignment for the four packet component identities."""

    gas_properties: cc.DeclaredGasProperties
    dry_matter_component_id: str
    residual_oil_component_id: str
    hexane_component_id: str
    water_component_id: str
    reference_temperature_k: float
    reference_pressure_pa: float
    dry_matter_offset_j_kg: float
    residual_oil_label_offset_j_kg: float
    hexane_native_to_common_offset_j_mol: float
    water_native_to_common_offset_j_mol: float
    hexane_common_datum_shift_j_mol: float
    water_common_datum_shift_j_mol: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    numeric_zero_alignment_only: ClassVar[bool] = True
    temperature_dependent_caloric_bridge_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.gas_properties) is not cc.DeclaredGasProperties:
            raise ComponentEnergyDatumConfigurationError(
                "gas_properties must be exact DeclaredGasProperties"
            )
        # CELL-02d W2 forgery gate: this adapter is a declared-convention
        # instrument; the native zero-gauge datum identity may never seed it.
        if self.gas_properties.energy_datum_id == NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID:
            raise ComponentEnergyDatumConfigurationError(
                "the native zero-gauge datum identity cannot seed the declared "
                "component-datum adapter"
            )
        for name, value in (
            ("dry-matter component", self.dry_matter_component_id),
            ("residual-oil component", self.residual_oil_component_id),
            ("hexane component", self.hexane_component_id),
            ("water component", self.water_component_id),
            ("common energy datum", self.gas_properties.energy_datum_id),
        ):
            _require_text(name, value)
        for name, value, positive in (
            ("reference temperature", self.reference_temperature_k, True),
            ("reference pressure", self.reference_pressure_pa, True),
            ("dry-matter numeric offset", self.dry_matter_offset_j_kg, False),
            ("residual-oil numeric offset", self.residual_oil_label_offset_j_kg, False),
            ("hexane numeric offset", self.hexane_native_to_common_offset_j_mol, False),
            ("water numeric offset", self.water_native_to_common_offset_j_mol, False),
            ("hexane common-datum shift", self.hexane_common_datum_shift_j_mol, False),
            ("water common-datum shift", self.water_common_datum_shift_j_mol, False),
        ):
            _require_float(name, value, positive=positive)
        if self.dry_matter_offset_j_kg != 0.0 or _float_bits(
            self.dry_matter_offset_j_kg
        ) != _float_bits(0.0):
            raise ComponentEnergyDatumConfigurationError(
                "the bounded bridge pins dry-matter numeric offset to positive zero"
            )
        if self.residual_oil_label_offset_j_kg != 0.0 or _float_bits(
            self.residual_oil_label_offset_j_kg
        ) != _float_bits(0.0):
            raise ComponentEnergyDatumConfigurationError(
                "the bounded bridge pins the non-additive oil-label offset to positive zero"
            )
        if not _same_float(self.gas_properties.hexane_molar_mass_kg_mol, hexane.M):
            raise ComponentEnergyDatumConfigurationError(
                "declared and native n-hexane component molar masses differ"
            )
        if not _same_float(self.gas_properties.water_molar_mass_kg_mol, water.M):
            raise ComponentEnergyDatumConfigurationError(
                "declared and native water component molar masses differ"
            )
        if self.reference_pressure_pa >= min(
            hexane.saturation_pressure(self.reference_temperature_k),
            water.saturation_pressure(self.reference_temperature_k),
        ):
            raise ComponentEnergyDatumConfigurationError(
                "component-datum reference pressure must keep both pure references vapor"
            )
        base_hexane, base_water = self._derived_base_offsets()
        expected_hexane = base_hexane + self.hexane_common_datum_shift_j_mol
        expected_water = base_water + self.water_common_datum_shift_j_mol
        if not _same_float(self.hexane_native_to_common_offset_j_mol, expected_hexane):
            raise ComponentEnergyDatumConfigurationError(
                "declared n-hexane numeric offset does not match the pinned native/Law-2 reference"
            )
        if not _same_float(self.water_native_to_common_offset_j_mol, expected_water):
            raise ComponentEnergyDatumConfigurationError(
                "declared water numeric offset does not match the pinned native/Law-2 reference"
            )

    @classmethod
    def derive(
        cls,
        *,
        gas_properties: cc.DeclaredGasProperties,
        dry_matter_component_id: str,
        residual_oil_component_id: str,
        hexane_component_id: str,
        water_component_id: str,
        reference_temperature_k: float,
        reference_pressure_pa: float,
        hexane_common_datum_shift_j_mol: float = 0.0,
        water_common_datum_shift_j_mol: float = 0.0,
    ) -> NumericComponentEnergyDatumAdapter:
        """Derive both phase-component offsets from existing property authorities."""

        if type(gas_properties) is not cc.DeclaredGasProperties:
            raise ComponentEnergyDatumConfigurationError(
                "gas_properties must be exact DeclaredGasProperties"
            )
        _require_float("reference temperature", reference_temperature_k, positive=True)
        _require_float("reference pressure", reference_pressure_pa, positive=True)
        _require_float("hexane common-datum shift", hexane_common_datum_shift_j_mol)
        _require_float("water common-datum shift", water_common_datum_shift_j_mol)
        native_hexane = (
            hexane.state_Tp(reference_temperature_k, reference_pressure_pa, "vapor").h_mass
            * hexane.M
        )
        native_water = (
            water.state_Tp(reference_temperature_k, reference_pressure_pa, "vapor").h_mass * water.M
        )
        declared_sensible = gas_properties.molar_heat_capacity_j_mol_k * (
            reference_temperature_k - gas_properties.energy_datum_temperature_k
        )
        return cls(
            gas_properties=gas_properties,
            dry_matter_component_id=dry_matter_component_id,
            residual_oil_component_id=residual_oil_component_id,
            hexane_component_id=hexane_component_id,
            water_component_id=water_component_id,
            reference_temperature_k=reference_temperature_k,
            reference_pressure_pa=reference_pressure_pa,
            dry_matter_offset_j_kg=0.0,
            residual_oil_label_offset_j_kg=0.0,
            hexane_native_to_common_offset_j_mol=(
                declared_sensible
                + gas_properties.hexane_latent_heat_j_mol
                - native_hexane
                + hexane_common_datum_shift_j_mol
            ),
            water_native_to_common_offset_j_mol=(
                declared_sensible
                + gas_properties.water_latent_heat_j_mol
                - native_water
                + water_common_datum_shift_j_mol
            ),
            hexane_common_datum_shift_j_mol=hexane_common_datum_shift_j_mol,
            water_common_datum_shift_j_mol=water_common_datum_shift_j_mol,
        )

    @property
    def energy_datum_id(self) -> str:
        return self.gas_properties.energy_datum_id

    @property
    def component_ids(self) -> tuple[str, str, str, str]:
        return (
            self.dry_matter_component_id,
            self.residual_oil_component_id,
            self.hexane_component_id,
            self.water_component_id,
        )

    @property
    def numeric_offsets(self) -> tuple[float, float, float, float]:
        return (
            self.dry_matter_offset_j_kg,
            self.residual_oil_label_offset_j_kg,
            self.hexane_native_to_common_offset_j_mol,
            self.water_native_to_common_offset_j_mol,
        )

    @property
    def common_datum_species_shifts_j_mol(self) -> tuple[float, float]:
        return (
            self.hexane_common_datum_shift_j_mol,
            self.water_common_datum_shift_j_mol,
        )

    @property
    def definition_digest(self) -> str:
        return _framed_digest(
            DATUM_ADAPTER_DIGEST_DOMAIN,
            (
                self.energy_datum_id.encode("utf-8"),
                *(value.encode("utf-8") for value in self.component_ids),
                _float_bits(self.reference_temperature_k),
                _float_bits(self.reference_pressure_pa),
                *(_float_bits(value) for value in self.numeric_offsets),
                *(_float_bits(value) for value in self.common_datum_species_shifts_j_mol),
                _float_bits(self.gas_properties.molar_heat_capacity_j_mol_k),
                _float_bits(self.gas_properties.energy_datum_temperature_k),
                _float_bits(self.gas_properties.hexane_latent_heat_j_mol),
                _float_bits(self.gas_properties.water_latent_heat_j_mol),
                *hexane.caloric_datum_signature(),
                *water.caloric_datum_signature(),
            ),
        )

    def require_gas_properties(self, properties: cc.DeclaredGasProperties) -> None:
        if type(properties) is not cc.DeclaredGasProperties or properties != self.gas_properties:
            raise ComponentEnergyDatumDriftError(
                "Law-2 gas properties drifted from the component-datum adapter"
            )

    def require_component_ids(self, component_ids: tuple[str, str, str, str]) -> None:
        if type(component_ids) is not tuple or component_ids != self.component_ids:
            raise ComponentEnergyDatumDriftError(
                "packet component identities drifted from the component-datum adapter"
            )

    def packet_offset_j(
        self,
        *,
        dry_matter_kg: float,
        residual_oil_label_kg: float,
        total_hexane_kg: float,
        total_water_kg: float,
    ) -> float:
        masses = (
            dry_matter_kg,
            residual_oil_label_kg,
            total_hexane_kg,
            total_water_kg,
        )
        if not all(
            type(value) is float and math.isfinite(value) and value >= 0.0 for value in masses
        ):
            raise ComponentEnergyDatumError(
                "component-datum packet masses must be nonnegative finite binary64"
            )
        return math.fsum(
            (
                dry_matter_kg * self.dry_matter_offset_j_kg,
                residual_oil_label_kg * self.residual_oil_label_offset_j_kg,
                total_hexane_kg
                / self.gas_properties.hexane_molar_mass_kg_mol
                * self.hexane_native_to_common_offset_j_mol,
                total_water_kg
                / self.gas_properties.water_molar_mass_kg_mol
                * self.water_native_to_common_offset_j_mol,
            )
        )

    def native_to_common_energy_j(self, native_energy_j: float, **masses: float) -> float:
        _require_float("native packet energy", native_energy_j)
        return math.fsum((native_energy_j, self.packet_offset_j(**masses)))

    def common_to_native_energy_j(self, common_energy_j: float, **masses: float) -> float:
        _require_float("common-datum packet energy", common_energy_j)
        return math.fsum((common_energy_j, -self.packet_offset_j(**masses)))

    def phase_transfer_offset_j(
        self,
        *,
        hexane_transfer_kg: float,
        water_transfer_kg: float,
    ) -> float:
        for name, value in (
            ("hexane transfer", hexane_transfer_kg),
            ("water transfer", water_transfer_kg),
        ):
            _require_float(name, value)
        return math.fsum(
            (
                hexane_transfer_kg
                / self.gas_properties.hexane_molar_mass_kg_mol
                * self.hexane_native_to_common_offset_j_mol,
                water_transfer_kg
                / self.gas_properties.water_molar_mass_kg_mol
                * self.water_native_to_common_offset_j_mol,
            )
        )

    # -----------------------------------------------------------------
    # F2 (2026-09-02) - THE CODEC CALORIC SEAM, made measurable
    # -----------------------------------------------------------------
    # Every phase change in the engineering lane crosses THIS boundary: a
    # mole of component i leaves a packet as LIQUID, booked in the native
    # frame lifted by ``*_native_to_common_offset_j_mol``, and arrives in a
    # gas control volume as VAPOUR, booked in the declared Law-2 convention
    # ``cp*(T - T_datum) + L_i``.  The crossing conserves energy if and only
    # if the two frames agree on the LIQUID enthalpy at the crossing state,
    # because the vapour enthalpy minus the booked latent IS the implied
    # liquid enthalpy on the gas side.
    #
    # THE DECLARED LATENT AND THE SPECIES SHIFT ARE A PURE GAUGE HERE, and
    # the algebra says so exactly.  Writing ``o_i`` for the offset,
    #
    #     o_i = cp*(T_ref - T_datum) + L_i + shift_i - M_i*h_nat_vap,i(T_ref, P_ref)
    #
    # the seam-consistent latent - the number the interface row would have
    # to consume for the crossing to close - is
    #
    #     L_i^seam(T, P) = [cp*(T - T_datum) + L_i + shift_i]
    #                      - [M_i*h_nat_liq,i(T, P) + o_i]
    #                    = cp*(T - T_ref) + M_i*h_nat_vap,i(T_ref, P_ref)
    #                      - M_i*h_nat_liq,i(T, P)
    #
    # in which BOTH ``L_i`` and ``shift_i`` have cancelled identically.  So "the declared latent is
    # 14 % low" is NOT the defect (wave-1 item 4, sharpened by wave 2): the
    # defect is that the declared arm books the GAUGE constant ``L_i`` where
    # this state function belongs, and the difference is the mis-booking.
    # ``seam_caloric_defect_j_mol`` is exactly that difference, and
    # ``seam_caloric_defect_w`` turns it into the watts it fabricates or
    # destroys at a given crossing rate.
    #
    # The same state function is what the PD dome-bypass boiling branch
    # ALREADY realizes implicitly (it balances the packet's own common-datum
    # energy against the departing vapour's declared enthalpy), which is why
    # F2 and F4 are one cluster and one law with two consumers.
    _SEAM_SPECIES: ClassVar[tuple[str, str]] = ("hexane", "water")

    def _species_module_and_offset(self, species: str) -> tuple[object, float, float, float]:
        """(property module, molar mass, native->common offset, declared species term).

        The declared species term is ``L_i + shift_i`` - exactly the group the
        production stream-enthalpy conventions carry
        (``sp1_k2_law2_tray_integration._bulk_gas_common_datum_enthalpy_j_mol``
        and ``pd_k_dome_bypass_tray_integration._declared_bulk_gas_enthalpy_j_mol``).
        """

        if species == "hexane":
            return (
                hexane,
                self.gas_properties.hexane_molar_mass_kg_mol,
                self.hexane_native_to_common_offset_j_mol,
                math.fsum(
                    (
                        self.gas_properties.hexane_latent_heat_j_mol,
                        self.hexane_common_datum_shift_j_mol,
                    )
                ),
            )
        if species == "water":
            return (
                water,
                self.gas_properties.water_molar_mass_kg_mol,
                self.water_native_to_common_offset_j_mol,
                math.fsum(
                    (
                        self.gas_properties.water_latent_heat_j_mol,
                        self.water_common_datum_shift_j_mol,
                    )
                ),
            )
        raise ComponentEnergyDatumConfigurationError(
            f"seam species must be one of {self._SEAM_SPECIES}, got {species!r}"
        )

    def _require_seam_state(self, temperature_k: float, pressure_pa: float) -> None:
        _require_float("seam temperature", temperature_k, positive=True)
        _require_float("seam pressure", pressure_pa, positive=True)

    def common_datum_liquid_enthalpy_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """The PACKET frame's liquid molar enthalpy, lifted to the common datum.

        This is the number a departing mole actually takes out of a packet's
        ``common_datum_energy_j``.
        """

        self._require_seam_state(temperature_k, pressure_pa)
        module, molar_mass, offset, _ = self._species_module_and_offset(species)
        return math.fsum(
            (module.state_Tp(temperature_k, pressure_pa, "liquid").h_mass * molar_mass, offset)
        )

    def common_datum_vapour_enthalpy_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """The PACKET frame's vapour molar enthalpy, lifted to the common datum.

        DECLARED HAZARD: at a state where the species is supersaturated
        (``P > Psat_i(T)``) this is the metastable vapour root of the frozen
        Helmholtz authority.  It is used only for the VAPOUR-ROW twin of the
        seam measurement, never inside ``seam_latent_j_mol``, which needs the
        liquid root alone.
        """

        self._require_seam_state(temperature_k, pressure_pa)
        module, molar_mass, offset, _ = self._species_module_and_offset(species)
        return math.fsum(
            (module.state_Tp(temperature_k, pressure_pa, "vapor").h_mass * molar_mass, offset)
        )

    def declared_vapour_enthalpy_j_mol(self, species: str, *, temperature_k: float) -> float:
        """The DECLARED convention's pure-species vapour molar enthalpy.

        ``cp*(T - T_datum) + L_i + shift_i`` - the same arithmetic the
        production stream-enthalpy helpers use at a pure-species composition.
        """

        _require_float("seam temperature", temperature_k, positive=True)
        _, _, _, declared_latent = self._species_module_and_offset(species)
        return math.fsum(
            (
                self.gas_properties.molar_heat_capacity_j_mol_k
                * (temperature_k - self.gas_properties.energy_datum_temperature_k),
                declared_latent,
            )
        )

    def seam_latent_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """The latent the interface row must consume for the crossing to close.

        ``declared vapour enthalpy - common-datum liquid enthalpy``.  The
        declared latent and the species shift cancel out of this quantity
        identically (see the block comment above), so it is a gauge-invariant
        state function of ``(T, P)`` and the adapter's reference point.
        """

        return math.fsum(
            (
                self.declared_vapour_enthalpy_j_mol(species, temperature_k=temperature_k),
                -self.common_datum_liquid_enthalpy_j_mol(
                    species, temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
            )
        )

    def seam_caloric_defect_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """``declared latent - seam latent``: the mis-booking per crossing mole.

        Positive means the interface row consumes MORE latent than the
        crossing actually costs, so the seam destroys energy; negative means
        it fabricates energy.  Exactly zero for a seam-consistent arm.
        """

        _, _, _, declared_latent = self._species_module_and_offset(species)
        return math.fsum(
            (
                declared_latent,
                -self.seam_latent_j_mol(
                    species, temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
            )
        )

    def seam_caloric_defect_w(
        self,
        *,
        temperature_k: float,
        pressure_pa: float,
        hexane_molar_rate_mol_s: float,
        water_molar_rate_mol_s: float,
    ) -> float:
        """The seam mis-booking in watts at a declared crossing rate.

        Rates are POSITIVE from the condensed carrier to the gas, matching
        ``release/equation_registry.yaml`` BAL-COMP-S.
        """

        for name, value in (
            ("hexane crossing rate", hexane_molar_rate_mol_s),
            ("water crossing rate", water_molar_rate_mol_s),
        ):
            _require_float(name, value)
        return math.fsum(
            (
                hexane_molar_rate_mol_s
                * self.seam_caloric_defect_j_mol(
                    "hexane", temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
                water_molar_rate_mol_s
                * self.seam_caloric_defect_j_mol(
                    "water", temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
            )
        )

    def _derived_base_offsets(self) -> tuple[float, float]:
        native_hexane = (
            hexane.state_Tp(
                self.reference_temperature_k,
                self.reference_pressure_pa,
                "vapor",
            ).h_mass
            * hexane.M
        )
        native_water = (
            water.state_Tp(
                self.reference_temperature_k,
                self.reference_pressure_pa,
                "vapor",
            ).h_mass
            * water.M
        )
        declared_sensible = self.gas_properties.molar_heat_capacity_j_mol_k * (
            self.reference_temperature_k - self.gas_properties.energy_datum_temperature_k
        )
        return (
            declared_sensible + self.gas_properties.hexane_latent_heat_j_mol - native_hexane,
            declared_sensible + self.gas_properties.water_latent_heat_j_mol - native_water,
        )


# =====================================================================
# Q-F2a - OWNER RULING 2026-09-02 (batch commit ``20f5028``): the
# engineering lane's codec-crossing caloric bookkeeping FLIPS from the
# declared constant-latent arm to the NATIVE caloric arm.
# =====================================================================
# The F246 ceremony record (``docs/GT_PS2_F246_CALORIC_CLUSTER_CEREMONY_
# 2026-09-02.md`` section 3.3) proved that the declared arm's seam defect
# is NOT locally repairable: the crossing enthalpy, the latent and the
# bulk gas enthalpy are three views of ONE caloric curve, so no single
# additive per-species offset can align two phases at once unless the
# declared latent already equals the native one.  The repair the owner
# ruled is (a): run the lane on the native caloric arm, where the packet
# frame and the gas frame are THE SAME frozen-law caloric curve and the
# crossing closes by construction.
#
# The native cell kernels (CELL-02b/02d/02e, ``cell_native_caloric``)
# were already dispatch-capable.  The PACKET CODEC was not: the declared
# bridge above is, by its own W2 forgery gate, forbidden from carrying
# the native zero-gauge datum identity - correctly, because a DECLARED
# constant-property closure claiming the native datum IS a forgery.  The
# gate is NOT relaxed here.  Instead this class is the native sibling of
# the same bridge: it carries NO declared gas properties at all, its four
# component gauges are exact positive zero, and it can carry ONLY the
# native zero-gauge datum identity.  The declared adapter still refuses
# the native identity, and this adapter refuses every other identity, so
# the two types partition the datum space instead of overlapping it.
#
# A ZERO GAUGE HAS NO REFERENCE STATE.  ``reference_temperature_k`` and
# ``reference_pressure_pa`` are therefore exact positive zero sentinels,
# not a physical (T, P): there is no alignment point to record, because
# the packet frame IS the native frame at every state.  They exist only
# because the tag-1 wire schema reserves two slots for them, and 0.0/0.0
# is the honest value to put there.
@dataclass(frozen=True, slots=True, kw_only=True)
class NativeZeroGaugeComponentEnergyDatumAdapter:
    """Identity component bridge for the native zero-gauge packet datum.

    The native-arm counterpart of :class:`NumericComponentEnergyDatumAdapter`.
    Every additive component gauge is exact positive zero, so the packet's
    ``common_datum_energy_j`` IS its native-frame energy and the codec
    crossing conserves energy identically rather than approximately.
    """

    dry_matter_component_id: str
    residual_oil_component_id: str
    hexane_component_id: str
    water_component_id: str

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    numeric_zero_alignment_only: ClassVar[bool] = True
    temperature_dependent_caloric_bridge_implemented: ClassVar[bool] = True
    declared_gas_properties_used: ClassVar[bool] = False
    constant_latent_heat_used: ClassVar[bool] = False
    owner_ruling_id: ClassVar[str] = "Q-F2a"
    owner_ruling_date: ClassVar[str] = "2026-09-02"
    owner_ruling_commit: ClassVar[str] = "20f5028"

    _SEAM_SPECIES: ClassVar[tuple[str, str]] = ("hexane", "water")

    def __post_init__(self) -> None:
        for name, value in (
            ("dry-matter component", self.dry_matter_component_id),
            ("residual-oil component", self.residual_oil_component_id),
            ("hexane component", self.hexane_component_id),
            ("water component", self.water_component_id),
        ):
            _require_text(name, value)

    # -- identity ------------------------------------------------------
    @property
    def energy_datum_id(self) -> str:
        return NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID

    @property
    def reference_temperature_k(self) -> float:
        """Exact positive zero: a zero gauge has no alignment state."""

        return 0.0

    @property
    def reference_pressure_pa(self) -> float:
        """Exact positive zero: a zero gauge has no alignment state."""

        return 0.0

    @property
    def dry_matter_offset_j_kg(self) -> float:
        return 0.0

    @property
    def residual_oil_label_offset_j_kg(self) -> float:
        return 0.0

    @property
    def hexane_native_to_common_offset_j_mol(self) -> float:
        return 0.0

    @property
    def water_native_to_common_offset_j_mol(self) -> float:
        return 0.0

    @property
    def hexane_common_datum_shift_j_mol(self) -> float:
        return 0.0

    @property
    def water_common_datum_shift_j_mol(self) -> float:
        return 0.0

    @property
    def component_ids(self) -> tuple[str, str, str, str]:
        return (
            self.dry_matter_component_id,
            self.residual_oil_component_id,
            self.hexane_component_id,
            self.water_component_id,
        )

    @property
    def numeric_offsets(self) -> tuple[float, float, float, float]:
        return (0.0, 0.0, 0.0, 0.0)

    @property
    def common_datum_species_shifts_j_mol(self) -> tuple[float, float]:
        return (0.0, 0.0)

    @property
    def definition_digest(self) -> str:
        return _framed_digest(
            NATIVE_ZERO_GAUGE_DATUM_ADAPTER_DIGEST_DOMAIN,
            (
                self.energy_datum_id.encode("utf-8"),
                *(value.encode("utf-8") for value in self.component_ids),
                _float_bits(self.reference_temperature_k),
                _float_bits(self.reference_pressure_pa),
                *(_float_bits(value) for value in self.numeric_offsets),
                *(_float_bits(value) for value in self.common_datum_species_shifts_j_mol),
                *hexane.caloric_datum_signature(),
                *water.caloric_datum_signature(),
            ),
        )

    # -- drift gates ---------------------------------------------------
    def require_gas_properties(self, properties: object) -> None:
        """Refuse anything but the exact native caloric closure."""

        if type(properties) is not nc.NativeGasProperties or (
            properties.energy_datum_id != NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID
        ):
            raise ComponentEnergyDatumDriftError(
                "the native zero-gauge component bridge admits only an exact "
                "NativeGasProperties on the native zero-gauge datum"
            )

    def require_component_ids(self, component_ids: tuple[str, str, str, str]) -> None:
        if type(component_ids) is not tuple or component_ids != self.component_ids:
            raise ComponentEnergyDatumDriftError(
                "packet component identities drifted from the native component bridge"
            )

    # -- gauge arithmetic (identity, with the same validation) ----------
    def packet_offset_j(
        self,
        *,
        dry_matter_kg: float,
        residual_oil_label_kg: float,
        total_hexane_kg: float,
        total_water_kg: float,
    ) -> float:
        masses = (dry_matter_kg, residual_oil_label_kg, total_hexane_kg, total_water_kg)
        if not all(
            type(value) is float and math.isfinite(value) and value >= 0.0 for value in masses
        ):
            raise ComponentEnergyDatumError(
                "component-datum packet masses must be nonnegative finite binary64"
            )
        return 0.0

    def native_to_common_energy_j(self, native_energy_j: float, **masses: float) -> float:
        _require_float("native packet energy", native_energy_j)
        self.packet_offset_j(**masses)
        return native_energy_j

    def common_to_native_energy_j(self, common_energy_j: float, **masses: float) -> float:
        _require_float("common-datum packet energy", common_energy_j)
        self.packet_offset_j(**masses)
        return common_energy_j

    def phase_transfer_offset_j(
        self,
        *,
        hexane_transfer_kg: float,
        water_transfer_kg: float,
    ) -> float:
        for name, value in (
            ("hexane transfer", hexane_transfer_kg),
            ("water transfer", water_transfer_kg),
        ):
            _require_float(name, value)
        return 0.0

    # -- the F2 seam, measured on the NATIVE arm ------------------------
    # Same vocabulary as the declared bridge's F2 block, so the seam
    # instrument compares like with like across the Q-F2a flip.
    def _species_module_and_mass(self, species: str) -> tuple[object, float]:
        if species == "hexane":
            return (hexane, hexane.M)
        if species == "water":
            return (water, water.M)
        raise ComponentEnergyDatumConfigurationError(
            f"seam species must be one of {self._SEAM_SPECIES}, got {species!r}"
        )

    def _require_seam_state(self, temperature_k: float, pressure_pa: float) -> None:
        _require_float("seam temperature", temperature_k, positive=True)
        _require_float("seam pressure", pressure_pa, positive=True)

    def common_datum_liquid_enthalpy_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """What a departing mole takes out of a packet's ``common_datum_energy_j``.

        On the zero gauge this IS the native liquid partial enthalpy.
        """

        self._require_seam_state(temperature_k, pressure_pa)
        module, molar_mass = self._species_module_and_mass(species)
        return module.state_Tp(temperature_k, pressure_pa, "liquid").h_mass * molar_mass

    def common_datum_vapour_enthalpy_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """The vapour twin; carries the declared metastable-root hazard."""

        self._require_seam_state(temperature_k, pressure_pa)
        module, molar_mass = self._species_module_and_mass(species)
        return module.state_Tp(temperature_k, pressure_pa, "vapor").h_mass * molar_mass

    def booked_latent_j_mol(self, species: str, *, temperature_k: float) -> float:
        """The latent the NATIVE interface row actually consumes.

        ``cell_engineering_feasibility``'s native arm books the frozen
        SATURATION latent at the interface temperature - the Span-Wagner
        / IAPWS-95 ``dh_vap(T)`` - not a declared constant.
        """

        _require_float("seam temperature", temperature_k, positive=True)
        if species == "hexane":
            return nc.hexane_latent_heat(temperature_k).value_j_mol
        if species == "water":
            return nc.water_latent_heat(temperature_k).value_j_mol
        raise ComponentEnergyDatumConfigurationError(
            f"seam species must be one of {self._SEAM_SPECIES}, got {species!r}"
        )

    def seam_latent_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """The latent the crossing actually costs at ``(T, P)``.

        ``M_i*[h_nat_vap,i(T,P) - h_nat_liq,i(T,P)]``: both roots on the
        SAME frozen law and the SAME zero gauge, so no datum constant
        survives and the quantity is a pure state function.
        """

        return math.fsum(
            (
                self.common_datum_vapour_enthalpy_j_mol(
                    species, temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
                -self.common_datum_liquid_enthalpy_j_mol(
                    species, temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
            )
        )

    def seam_caloric_defect_j_mol(
        self, species: str, *, temperature_k: float, pressure_pa: float
    ) -> float:
        """``booked latent - seam latent``: the residue after the Q-F2a flip.

        On the declared arm this is the constant-versus-state-function
        mis-booking (thousands of J/mol).  On the native arm both terms
        come from the same frozen law and the residue is ONLY the
        saturation-versus-``(T, P)`` root difference: the interface row
        books ``dh_vap(T)`` at ``Psat(T)`` while the crossing happens at
        the tray pressure ``P``.
        """

        return math.fsum(
            (
                self.booked_latent_j_mol(species, temperature_k=temperature_k),
                -self.seam_latent_j_mol(
                    species, temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
            )
        )

    def seam_caloric_defect_w(
        self,
        *,
        temperature_k: float,
        pressure_pa: float,
        hexane_molar_rate_mol_s: float,
        water_molar_rate_mol_s: float,
    ) -> float:
        """The post-flip seam residue in watts at a declared crossing rate."""

        for name, value in (
            ("hexane crossing rate", hexane_molar_rate_mol_s),
            ("water crossing rate", water_molar_rate_mol_s),
        ):
            _require_float(name, value)
        return math.fsum(
            (
                hexane_molar_rate_mol_s
                * self.seam_caloric_defect_j_mol(
                    "hexane", temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
                water_molar_rate_mol_s
                * self.seam_caloric_defect_j_mol(
                    "water", temperature_k=temperature_k, pressure_pa=pressure_pa
                ),
            )
        )


#: The exact pair of component bridges the tag-1 codec admits.  The
#: declared bridge refuses the native datum identity (W2 forgery gate);
#: the native bridge refuses every other identity.  Adding a third type
#: is a Class-B event, not an import.
COMPONENT_ENERGY_DATUM_BRIDGE_UNION = (
    NumericComponentEnergyDatumAdapter,
    NativeZeroGaugeComponentEnergyDatumAdapter,
)


__all__ = (
    "COMPONENT_ENERGY_DATUM_BRIDGE_UNION",
    "ComponentEnergyDatumConfigurationError",
    "ComponentEnergyDatumDriftError",
    "ComponentEnergyDatumError",
    "DATUM_ADAPTER_DIGEST_DOMAIN",
    "NATIVE_ZERO_GAUGE_DATUM_ADAPTER_DIGEST_DOMAIN",
    "NativeZeroGaugeComponentEnergyDatumAdapter",
    "NumericComponentEnergyDatumAdapter",
)
