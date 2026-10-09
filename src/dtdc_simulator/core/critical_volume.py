"""Porosity-dependent critical-volume and phase-ledger kernel.

This module implements the constitutive bookkeeping selected for the
soybean-meal transition loading.  It is intentionally independent of the
production particle and tower solvers so the ledger can be qualified before
it is wired into regime switching, initialization, mass, and energy.

All extensive quantities use one kilogram of composite dry meal as the
reference basis.  The effective critical occupied volume is

``nu_c_eff = kappa_v * epsilon_p / rho_dm,p``.

``nu_c_eff`` is a correlation volume, not a declaration that the measured
dry particle has swollen to that geometry.  The difference from the measured
dry geometric pore volume is exposed as a diagnostic and must be assigned by
a separately approved wet-geometry model before it can alter dimensions.
``kappa_v`` is a material-specific multiplier anchored once at the reference
meal.  It remains fixed when porosity is varied; recalibrating it at every
porosity would erase the approved porosity dependence.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


FANER_SOYBEAN_REFERENCE_POROSITY = 0.141
# 1160 kg/m3 is the rounded literature value.  The executable canonical
# identity uses (1 - 0.141) * 1350 = 1159.65 kg/m3 exactly.
FANER_SOYBEAN_REFERENCE_DRY_MEAL_DENSITY_KG_M3 = 1159.65
FANER_SOYBEAN_REFERENCE_CRITICAL_HEXANE_LOADING = 0.20


def _require_finite(name: str, value: float) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")


def _require_nonnegative(name: str, value: float) -> None:
    _require_finite(name, value)
    if value < 0.0:
        raise ValueError(f"{name} must be non-negative, got {value!r}")


def _require_positive(name: str, value: float) -> None:
    _require_finite(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")


@dataclass(frozen=True)
class CriticalVolumeBasis:
    """Material basis and effective-volume scale at one local state.

    ``dry_meal_envelope_density_kg_m3`` is dry composite-meal mass divided by
    the measured dry particle-envelope volume.  Consequently,
    ``epsilon_p / rho_dm,p`` is the measured dry pore volume per kilogram of
    composite dry meal.
    """

    particle_porosity: float
    dry_meal_envelope_density_kg_m3: float
    critical_volume_multiplier: float

    def __post_init__(self) -> None:
        _require_finite("particle_porosity", self.particle_porosity)
        if not 0.0 < self.particle_porosity < 1.0:
            raise ValueError(
                "particle_porosity must be strictly between zero and one, "
                f"got {self.particle_porosity!r}"
            )
        _require_positive(
            "dry_meal_envelope_density_kg_m3",
            self.dry_meal_envelope_density_kg_m3,
        )
        _require_positive(
            "critical_volume_multiplier",
            self.critical_volume_multiplier,
        )

    @classmethod
    def from_skeletal_density(
        cls,
        *,
        particle_porosity: float,
        dry_solid_skeletal_density_kg_m3: float,
        critical_volume_multiplier: float,
    ) -> CriticalVolumeBasis:
        """Build the canonical basis with ``rho_dm=(1-epsilon)*rho_skel``."""

        _require_positive(
            "dry_solid_skeletal_density_kg_m3",
            dry_solid_skeletal_density_kg_m3,
        )
        _require_finite("particle_porosity", particle_porosity)
        return cls(
            particle_porosity=particle_porosity,
            dry_meal_envelope_density_kg_m3=(
                (1.0 - particle_porosity)
                * dry_solid_skeletal_density_kg_m3
            ),
            critical_volume_multiplier=critical_volume_multiplier,
        )

    @property
    def implied_dry_solid_skeletal_density_kg_m3(self) -> float:
        """Skeletal density implied by the canonical volume identity."""

        return (
            self.dry_meal_envelope_density_kg_m3
            / (1.0 - self.particle_porosity)
        )

    @property
    def dry_geometric_pore_specific_volume_m3_per_kg_dry_meal(self) -> float:
        """Measured dry pore volume; no wet swelling is implied."""

        return (
            self.particle_porosity
            / self.dry_meal_envelope_density_kg_m3
        )

    @property
    def effective_critical_specific_volume_m3_per_kg_dry_meal(self) -> float:
        """Correlation volume available to the critical phase ledger."""

        return (
            self.critical_volume_multiplier
            * self.dry_geometric_pore_specific_volume_m3_per_kg_dry_meal
        )

    @property
    def effective_minus_dry_pore_volume_m3_per_kg_dry_meal(self) -> float:
        """Signed diagnostic volume that has no geometric allocation yet."""

        return (
            self.effective_critical_specific_volume_m3_per_kg_dry_meal
            - self.dry_geometric_pore_specific_volume_m3_per_kg_dry_meal
        )

    @property
    def diagnostic_effective_swelling_volume_m3_per_kg_dry_meal(self) -> float:
        """Positive excess correlation volume, not an actual swelling law."""

        return max(
            self.effective_minus_dry_pore_volume_m3_per_kg_dry_meal,
            0.0,
        )

    @property
    def effective_volume_porosity_derivative_at_fixed_density(
        self,
    ) -> float:
        """Return ``d(nu_c_eff)/d(epsilon_p)`` at fixed dry density."""

        return (
            self.critical_volume_multiplier
            / self.dry_meal_envelope_density_kg_m3
        )

    @property
    def effective_volume_porosity_derivative_at_fixed_skeletal_density(
        self,
    ) -> float:
        """Derivative along ``rho_dm=(1-epsilon_p)*rho_skeletal``."""

        return (
            self.critical_volume_multiplier
            / (
                self.implied_dry_solid_skeletal_density_kg_m3
                * (1.0 - self.particle_porosity) ** 2
            )
        )


@dataclass(frozen=True)
class CriticalFixedPhaseOccupancy:
    """Critical-state phases specified before mobile liquid is recovered.

    Pore gas is represented by its *total mixture volume* and its hexane
    loading separately.  This avoids adding fictitious partial volumes for
    components sharing one gas volume.
    """

    retained_hexane_loading_kg_per_kg_dry_meal: float
    pore_gas_hexane_loading_kg_per_kg_dry_meal: float
    free_water_loading_kg_per_kg_dry_meal: float
    pore_gas_specific_volume_m3_per_kg_dry_meal: float
    liquid_hexane_density_kg_m3: float
    liquid_water_density_kg_m3: float

    def __post_init__(self) -> None:
        for name in (
            "retained_hexane_loading_kg_per_kg_dry_meal",
            "pore_gas_hexane_loading_kg_per_kg_dry_meal",
            "free_water_loading_kg_per_kg_dry_meal",
            "pore_gas_specific_volume_m3_per_kg_dry_meal",
        ):
            _require_nonnegative(name, getattr(self, name))
        _require_positive(
            "liquid_hexane_density_kg_m3",
            self.liquid_hexane_density_kg_m3,
        )
        _require_positive(
            "liquid_water_density_kg_m3",
            self.liquid_water_density_kg_m3,
        )
        if (
            self.pore_gas_hexane_loading_kg_per_kg_dry_meal > 0.0
            and self.pore_gas_specific_volume_m3_per_kg_dry_meal == 0.0
        ):
            raise ValueError(
                "positive pore-gas hexane requires positive total pore-gas volume"
            )

    @property
    def retained_hexane_specific_volume_m3_per_kg_dry_meal(self) -> float:
        """Liquid-like retained partial volume selected under PHY-051."""

        return (
            self.retained_hexane_loading_kg_per_kg_dry_meal
            / self.liquid_hexane_density_kg_m3
        )

    @property
    def free_water_specific_volume_m3_per_kg_dry_meal(self) -> float:
        return (
            self.free_water_loading_kg_per_kg_dry_meal
            / self.liquid_water_density_kg_m3
        )

    @property
    def fixed_occupied_specific_volume_m3_per_kg_dry_meal(self) -> float:
        return (
            self.retained_hexane_specific_volume_m3_per_kg_dry_meal
            + self.free_water_specific_volume_m3_per_kg_dry_meal
            + self.pore_gas_specific_volume_m3_per_kg_dry_meal
        )


@dataclass(frozen=True)
class CriticalPhasePartition:
    """Explicit phase inventory at the critical transition."""

    retained_hexane_loading_kg_per_kg_dry_meal: float
    mobile_liquid_hexane_loading_kg_per_kg_dry_meal: float
    pore_gas_hexane_loading_kg_per_kg_dry_meal: float
    free_water_loading_kg_per_kg_dry_meal: float
    pore_gas_specific_volume_m3_per_kg_dry_meal: float
    liquid_hexane_density_kg_m3: float
    liquid_water_density_kg_m3: float

    def __post_init__(self) -> None:
        fixed = self.fixed_occupancy
        _require_nonnegative(
            "mobile_liquid_hexane_loading_kg_per_kg_dry_meal",
            self.mobile_liquid_hexane_loading_kg_per_kg_dry_meal,
        )
        # Materialize validation of the nested fixed state in __post_init__.
        _ = fixed.fixed_occupied_specific_volume_m3_per_kg_dry_meal

    @property
    def fixed_occupancy(self) -> CriticalFixedPhaseOccupancy:
        return CriticalFixedPhaseOccupancy(
            retained_hexane_loading_kg_per_kg_dry_meal=(
                self.retained_hexane_loading_kg_per_kg_dry_meal
            ),
            pore_gas_hexane_loading_kg_per_kg_dry_meal=(
                self.pore_gas_hexane_loading_kg_per_kg_dry_meal
            ),
            free_water_loading_kg_per_kg_dry_meal=(
                self.free_water_loading_kg_per_kg_dry_meal
            ),
            pore_gas_specific_volume_m3_per_kg_dry_meal=(
                self.pore_gas_specific_volume_m3_per_kg_dry_meal
            ),
            liquid_hexane_density_kg_m3=self.liquid_hexane_density_kg_m3,
            liquid_water_density_kg_m3=self.liquid_water_density_kg_m3,
        )

    @property
    def mobile_liquid_hexane_specific_volume_m3_per_kg_dry_meal(
        self,
    ) -> float:
        return (
            self.mobile_liquid_hexane_loading_kg_per_kg_dry_meal
            / self.liquid_hexane_density_kg_m3
        )

    @property
    def critical_hexane_loading_kg_per_kg_dry_meal(self) -> float:
        """Total internal hexane; attached/interparticle ``X_f`` is absent."""

        return (
            self.retained_hexane_loading_kg_per_kg_dry_meal
            + self.mobile_liquid_hexane_loading_kg_per_kg_dry_meal
            + self.pore_gas_hexane_loading_kg_per_kg_dry_meal
        )

    @property
    def occupied_specific_volume_m3_per_kg_dry_meal(self) -> float:
        return (
            self.fixed_occupancy.fixed_occupied_specific_volume_m3_per_kg_dry_meal
            + self.mobile_liquid_hexane_specific_volume_m3_per_kg_dry_meal
        )


@dataclass(frozen=True)
class DrySideFrontInventory:
    """Hexane inventory immediately outside the dry side of a wet front.

    The pore-gas loading and total gas volume are both explicit.  Their ratio
    is the local hexane mass concentration in the shared gas-mixture volume;
    it is not a pure-component gas volume.
    """

    retained_hexane_loading_kg_per_kg_dry_meal: float
    pore_gas_hexane_loading_kg_per_kg_dry_meal: float
    pore_gas_specific_volume_m3_per_kg_dry_meal: float

    def __post_init__(self) -> None:
        _require_nonnegative(
            "retained_hexane_loading_kg_per_kg_dry_meal",
            self.retained_hexane_loading_kg_per_kg_dry_meal,
        )
        _require_nonnegative(
            "pore_gas_hexane_loading_kg_per_kg_dry_meal",
            self.pore_gas_hexane_loading_kg_per_kg_dry_meal,
        )
        _require_nonnegative(
            "pore_gas_specific_volume_m3_per_kg_dry_meal",
            self.pore_gas_specific_volume_m3_per_kg_dry_meal,
        )
        if (
            self.pore_gas_hexane_loading_kg_per_kg_dry_meal > 0.0
            and self.pore_gas_specific_volume_m3_per_kg_dry_meal == 0.0
        ):
            raise ValueError(
                "positive dry-side pore-gas hexane requires positive gas volume"
            )

    @classmethod
    def from_pore_gas_hexane_density(
        cls,
        *,
        retained_hexane_loading_kg_per_kg_dry_meal: float,
        pore_gas_hexane_density_kg_m3: float,
        pore_gas_specific_volume_m3_per_kg_dry_meal: float,
    ) -> DrySideFrontInventory:
        """Build ``X_h,g=nu_g*c_h,g`` without duplicating the gas mass."""

        _require_nonnegative(
            "pore_gas_hexane_density_kg_m3",
            pore_gas_hexane_density_kg_m3,
        )
        return cls(
            retained_hexane_loading_kg_per_kg_dry_meal=(
                retained_hexane_loading_kg_per_kg_dry_meal
            ),
            pore_gas_hexane_loading_kg_per_kg_dry_meal=(
                pore_gas_hexane_density_kg_m3
                * pore_gas_specific_volume_m3_per_kg_dry_meal
            ),
            pore_gas_specific_volume_m3_per_kg_dry_meal=(
                pore_gas_specific_volume_m3_per_kg_dry_meal
            ),
        )

    @property
    def total_hexane_loading_kg_per_kg_dry_meal(self) -> float:
        return (
            self.retained_hexane_loading_kg_per_kg_dry_meal
            + self.pore_gas_hexane_loading_kg_per_kg_dry_meal
        )

    @property
    def pore_gas_hexane_density_kg_m3(self) -> float | None:
        if self.pore_gas_specific_volume_m3_per_kg_dry_meal == 0.0:
            return None
        return (
            self.pore_gas_hexane_loading_kg_per_kg_dry_meal
            / self.pore_gas_specific_volume_m3_per_kg_dry_meal
        )


@dataclass(frozen=True)
class CriticalPhaseLedger:
    """Closed mass/volume ledger recovered from an effective volume."""

    basis: CriticalVolumeBasis
    partition: CriticalPhasePartition
    phase_volume_residual_m3_per_kg_dry_meal: float

    @property
    def critical_hexane_loading_kg_per_kg_dry_meal(self) -> float:
        return self.partition.critical_hexane_loading_kg_per_kg_dry_meal

    @property
    def critical_loading_porosity_derivative_at_fixed_density(
        self,
    ) -> float:
        """Derivative with fixed non-mobile phase occupancy and properties."""

        return (
            self.basis.effective_volume_porosity_derivative_at_fixed_density
            * self.partition.liquid_hexane_density_kg_m3
        )

    @property
    def critical_loading_porosity_derivative_at_fixed_skeletal_density(
        self,
    ) -> float:
        """Derivative for the canonical porosity-density coupling."""

        return (
            self.basis.effective_volume_porosity_derivative_at_fixed_skeletal_density
            * self.partition.liquid_hexane_density_kg_m3
        )

@dataclass(frozen=True)
class HighLoadingHexanePartition:
    """Internal critical inventory plus separate attached/interparticle liquid."""

    total_hexane_loading_kg_per_kg_dry_meal: float
    internal_critical_hexane_loading_kg_per_kg_dry_meal: float
    attached_interparticle_hexane_loading_kg_per_kg_dry_meal: float
    mass_residual_kg_per_kg_dry_meal: float


@dataclass(frozen=True)
class FrontMobileJumpLedger:
    """Exact wet-minus-dry component jump at the moving front."""

    wet_critical_ledger: CriticalPhaseLedger
    dry_side_inventory: DrySideFrontInventory
    wet_hexane_concentration_kg_m3_particle_envelope: float
    dry_side_hexane_concentration_kg_m3_particle_envelope: float
    mobile_jump_loading_kg_per_kg_dry_meal: float
    mobile_jump_density_kg_m3_particle_envelope: float
    loading_balance_residual_kg_per_kg_dry_meal: float
    concentration_balance_residual_kg_m3_particle_envelope: float

    @property
    def mobile_liquid_reserve_density_kg_m3_particle_envelope(self) -> float:
        """Reserve-density spelling for production radial-state adapters."""

        return self.mobile_jump_density_kg_m3_particle_envelope


def clean_reference_critical_volume_multiplier(
    *,
    liquid_hexane_density_kg_m3: float,
    reference_porosity: float = FANER_SOYBEAN_REFERENCE_POROSITY,
    reference_dry_meal_density_kg_m3: float = (
        FANER_SOYBEAN_REFERENCE_DRY_MEAL_DENSITY_KG_M3
    ),
    reference_critical_hexane_loading: float = (
        FANER_SOYBEAN_REFERENCE_CRITICAL_HEXANE_LOADING
    ),
) -> float:
    """Anchor material-specific ``kappa_v`` in the clean liquid-like limit.

    This clean limiting calibration has no free water or pore gas.  Retained
    and mobile hexane have the same liquid-like partial volume, so their
    internal split does not affect the result.
    """

    basis = CriticalVolumeBasis(
        particle_porosity=reference_porosity,
        dry_meal_envelope_density_kg_m3=reference_dry_meal_density_kg_m3,
        critical_volume_multiplier=1.0,
    )
    _require_positive(
        "liquid_hexane_density_kg_m3",
        liquid_hexane_density_kg_m3,
    )
    _require_nonnegative(
        "reference_critical_hexane_loading",
        reference_critical_hexane_loading,
    )
    if reference_critical_hexane_loading == 0.0:
        raise ValueError("reference critical hexane loading must be positive")
    reference_occupied_volume = (
        reference_critical_hexane_loading
        / liquid_hexane_density_kg_m3
    )
    return (
        reference_occupied_volume
        / basis.dry_geometric_pore_specific_volume_m3_per_kg_dry_meal
    )


def clean_wet_critical_hexane_loading(
    basis: CriticalVolumeBasis,
    *,
    liquid_hexane_density_kg_m3: float,
) -> float:
    """Return the clean wet transition loading at one local liquid density."""

    _require_positive(
        "liquid_hexane_density_kg_m3",
        liquid_hexane_density_kg_m3,
    )
    return (
        basis.effective_critical_specific_volume_m3_per_kg_dry_meal
        * liquid_hexane_density_kg_m3
    )


def clean_wet_critical_hexane_loading_from_skeletal_basis(
    *,
    particle_porosity: float,
    dry_solid_skeletal_density_kg_m3: float,
    critical_volume_multiplier: float,
    liquid_hexane_density_kg_m3: float,
) -> float:
    """Production adapter from canonical particle properties to ``X_c``."""

    return clean_wet_critical_hexane_loading(
        CriticalVolumeBasis.from_skeletal_density(
            particle_porosity=particle_porosity,
            dry_solid_skeletal_density_kg_m3=(
                dry_solid_skeletal_density_kg_m3
            ),
            critical_volume_multiplier=critical_volume_multiplier,
        ),
        liquid_hexane_density_kg_m3=liquid_hexane_density_kg_m3,
    )


def reference_critical_volume_multiplier_from_phase_partition(
    partition: CriticalPhasePartition,
    *,
    reference_porosity: float = FANER_SOYBEAN_REFERENCE_POROSITY,
    reference_dry_meal_density_kg_m3: float = (
        FANER_SOYBEAN_REFERENCE_DRY_MEAL_DENSITY_KG_M3
    ),
    reference_critical_hexane_loading: float = (
        FANER_SOYBEAN_REFERENCE_CRITICAL_HEXANE_LOADING
    ),
    loading_tolerance: float = 1.0e-12,
) -> float:
    """Anchor material-specific ``kappa_v`` from one reference phase state."""

    _require_nonnegative("loading_tolerance", loading_tolerance)
    _require_nonnegative(
        "reference_critical_hexane_loading",
        reference_critical_hexane_loading,
    )
    loading_residual = (
        partition.critical_hexane_loading_kg_per_kg_dry_meal
        - reference_critical_hexane_loading
    )
    if abs(loading_residual) > loading_tolerance:
        raise ValueError(
            "reference phase partition does not recover the specified "
            f"critical loading; residual={loading_residual!r}"
        )
    raw_basis = CriticalVolumeBasis(
        particle_porosity=reference_porosity,
        dry_meal_envelope_density_kg_m3=reference_dry_meal_density_kg_m3,
        critical_volume_multiplier=1.0,
    )
    return (
        partition.occupied_specific_volume_m3_per_kg_dry_meal
        / raw_basis.dry_geometric_pore_specific_volume_m3_per_kg_dry_meal
    )


def solve_critical_phase_ledger(
    basis: CriticalVolumeBasis,
    fixed_occupancy: CriticalFixedPhaseOccupancy,
) -> CriticalPhaseLedger:
    """Fill the remaining effective volume with mobile liquid hexane.

    Retained hexane, free liquid water, and the total pore-gas volume are
    allocated first.  A negative remaining volume is infeasible and is never
    hidden by clipping.
    """

    effective_volume = (
        basis.effective_critical_specific_volume_m3_per_kg_dry_meal
    )
    fixed_volume = (
        fixed_occupancy.fixed_occupied_specific_volume_m3_per_kg_dry_meal
    )
    remaining_volume = effective_volume - fixed_volume
    roundoff_tolerance = 32.0 * math.ulp(
        max(effective_volume, fixed_volume, 1.0e-300)
    )
    if remaining_volume < -roundoff_tolerance:
        raise ValueError(
            "retained hexane, free water, and pore gas exceed the effective "
            f"critical volume by {-remaining_volume!r} m3/kg dry meal"
        )
    mobile_volume = max(remaining_volume, 0.0)
    mobile_loading = (
        mobile_volume * fixed_occupancy.liquid_hexane_density_kg_m3
    )
    partition = CriticalPhasePartition(
        retained_hexane_loading_kg_per_kg_dry_meal=(
            fixed_occupancy.retained_hexane_loading_kg_per_kg_dry_meal
        ),
        mobile_liquid_hexane_loading_kg_per_kg_dry_meal=mobile_loading,
        pore_gas_hexane_loading_kg_per_kg_dry_meal=(
            fixed_occupancy.pore_gas_hexane_loading_kg_per_kg_dry_meal
        ),
        free_water_loading_kg_per_kg_dry_meal=(
            fixed_occupancy.free_water_loading_kg_per_kg_dry_meal
        ),
        pore_gas_specific_volume_m3_per_kg_dry_meal=(
            fixed_occupancy.pore_gas_specific_volume_m3_per_kg_dry_meal
        ),
        liquid_hexane_density_kg_m3=(
            fixed_occupancy.liquid_hexane_density_kg_m3
        ),
        liquid_water_density_kg_m3=(
            fixed_occupancy.liquid_water_density_kg_m3
        ),
    )
    residual = (
        partition.occupied_specific_volume_m3_per_kg_dry_meal
        - effective_volume
    )
    return CriticalPhaseLedger(
        basis=basis,
        partition=partition,
        phase_volume_residual_m3_per_kg_dry_meal=residual,
    )


def solve_front_mobile_jump_ledger(
    *,
    wet_critical_ledger: CriticalPhaseLedger,
    dry_side_inventory: DrySideFrontInventory,
) -> FrontMobileJumpLedger:
    """Return ``Delta C_h=rho_dm*(X_c,wet-X_h,dry)`` exactly.

    The dry-side inventory is independent of any optional pore gas carried in
    the wet critical-state volume ledger.  In particular, the dry-side gas
    contribution may not be inferred from, or hidden inside, the wet mobile
    liquid amount.
    """

    wet_loading = (
        wet_critical_ledger.critical_hexane_loading_kg_per_kg_dry_meal
    )
    dry_loading = (
        dry_side_inventory.total_hexane_loading_kg_per_kg_dry_meal
    )
    mobile_jump_loading = wet_loading - dry_loading
    roundoff_tolerance = 32.0 * math.ulp(
        max(wet_loading, dry_loading, 1.0e-300)
    )
    if mobile_jump_loading < -roundoff_tolerance:
        raise ValueError(
            "dry-side retained plus pore-gas hexane exceeds the wet critical "
            f"loading by {-mobile_jump_loading!r} kg/kg dry meal"
        )
    mobile_jump_loading = max(mobile_jump_loading, 0.0)
    dry_density = (
        wet_critical_ledger.basis.dry_meal_envelope_density_kg_m3
    )
    wet_concentration = dry_density * wet_loading
    dry_concentration = dry_density * dry_loading
    mobile_jump_density = dry_density * mobile_jump_loading
    loading_residual = (
        dry_loading + mobile_jump_loading - wet_loading
    )
    concentration_residual = (
        dry_concentration + mobile_jump_density - wet_concentration
    )
    return FrontMobileJumpLedger(
        wet_critical_ledger=wet_critical_ledger,
        dry_side_inventory=dry_side_inventory,
        wet_hexane_concentration_kg_m3_particle_envelope=wet_concentration,
        dry_side_hexane_concentration_kg_m3_particle_envelope=(
            dry_concentration
        ),
        mobile_jump_loading_kg_per_kg_dry_meal=mobile_jump_loading,
        mobile_jump_density_kg_m3_particle_envelope=mobile_jump_density,
        loading_balance_residual_kg_per_kg_dry_meal=loading_residual,
        concentration_balance_residual_kg_m3_particle_envelope=(
            concentration_residual
        ),
    )


def split_attached_interparticle_hexane(
    *,
    total_hexane_loading_kg_per_kg_dry_meal: float,
    critical_ledger: CriticalPhaseLedger,
) -> HighLoadingHexanePartition:
    """Keep solvent above ``X_c`` as a separate attached inventory ``X_f``."""

    _require_nonnegative(
        "total_hexane_loading_kg_per_kg_dry_meal",
        total_hexane_loading_kg_per_kg_dry_meal,
    )
    critical_loading = (
        critical_ledger.critical_hexane_loading_kg_per_kg_dry_meal
    )
    if total_hexane_loading_kg_per_kg_dry_meal < critical_loading:
        raise ValueError(
            "attached/interparticle partition is only defined for "
            "total hexane loading greater than or equal to X_c"
        )
    attached_loading = (
        total_hexane_loading_kg_per_kg_dry_meal - critical_loading
    )
    mass_residual = (
        critical_loading
        + attached_loading
        - total_hexane_loading_kg_per_kg_dry_meal
    )
    return HighLoadingHexanePartition(
        total_hexane_loading_kg_per_kg_dry_meal=(
            total_hexane_loading_kg_per_kg_dry_meal
        ),
        internal_critical_hexane_loading_kg_per_kg_dry_meal=critical_loading,
        attached_interparticle_hexane_loading_kg_per_kg_dry_meal=(
            attached_loading
        ),
        mass_residual_kg_per_kg_dry_meal=mass_residual,
    )


__all__ = [
    "FANER_SOYBEAN_REFERENCE_CRITICAL_HEXANE_LOADING",
    "FANER_SOYBEAN_REFERENCE_DRY_MEAL_DENSITY_KG_M3",
    "FANER_SOYBEAN_REFERENCE_POROSITY",
    "CriticalFixedPhaseOccupancy",
    "CriticalPhaseLedger",
    "CriticalPhasePartition",
    "CriticalVolumeBasis",
    "DrySideFrontInventory",
    "FrontMobileJumpLedger",
    "HighLoadingHexanePartition",
    "clean_reference_critical_volume_multiplier",
    "clean_wet_critical_hexane_loading",
    "clean_wet_critical_hexane_loading_from_skeletal_basis",
    "reference_critical_volume_multiplier_from_phase_partition",
    "solve_critical_phase_ledger",
    "solve_front_mobile_jump_ledger",
    "split_attached_interparticle_hexane",
]
