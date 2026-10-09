r"""Non-qualifying Gate-1g binary-pore coefficient foundation.

The frozen G1G closure writes the independent pore-gas flux as

``J_w = -J_h = -L_wh_eff grad[(mu_w-mu_h)/(R T)]``.

For an ideal binary gas, the Onsager mobility has the composition degeneracy

``L_wh_eff = c_g D_wh_eff y_w y_h``.

This module evaluates that identity without assigning a soybean-DT transport
coefficient.  The caller must supply a positive, provenance-bearing interval
for the *already pore-effective* binary diffusivity ``D_wh_eff``.  In
particular, porosity is not multiplied into the face mobility a second time
(PHY-020), and no retained-water or surface-film coefficient is inferred.
Every object remains ``physically_qualifying=False``.

The centered face composition factor is the symmetric discrete chain-rule
mean

``m_wh = (y_w,R-y_w,L)/(logit(y_w,R)-logit(y_w,L))``.

At equal states this is exactly ``y_w y_h``.  It is the discrete
thermodynamic face value of that product: for an isothermal ideal binary gas
it makes the log-fugacity flux exactly equal to the centered Fick flux for any
strictly interior pair.  Its pure-component limit is exactly zero, without a
composition floor or clipping.  A face touching exact purity has a zero
mobility and an infinite logit force; callers must evaluate the combined
ideal/Fick limit rather than multiply those two quantities separately.

For a strictly interior face, :meth:`EvaluatedBinaryMobilitySelection.
as_coupled_pore_selection` returns the state-specific positive
``coupled_pore.MobilitySelection`` expected by
``independent_maxwell_stefan_flux``.  Exact-purity states are limiting-oracle
states outside the current ``EquilibriumPoreState`` domain and deliberately
cannot be converted through that logarithmic-force API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dtdc_simulator.core2.particle import coupled_pore


_COMPOSITION_SUM_ABS_TOL = 1.0e-12


class PureComponentMobilityLimitError(ValueError):
    """A split zero-mobility/infinite-force evaluation was requested."""


@dataclass(frozen=True)
class EffectiveBinaryDiffusivityInterval:
    """Caller-supplied positive interval for pore-effective ``D_wh``.

    Endpoints have units of m2/s.  The interval carries numerical/model-form
    provenance only; it is not an identified soybean-DT physical band.
    """

    lower_m2_s: float
    upper_m2_s: float
    provenance: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        endpoints = (self.lower_m2_s, self.upper_m2_s)
        if not all(math.isfinite(value) and value > 0.0 for value in endpoints):
            raise ValueError("effective binary diffusivity endpoints must be positive and finite")
        if self.lower_m2_s > self.upper_m2_s:
            raise ValueError("effective binary diffusivity endpoints must be ordered")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("effective binary diffusivity interval needs explicit provenance")

    def select(self, fraction: float) -> "EffectiveBinaryDiffusivitySelection":
        """Select a declared interval coordinate without clipping it."""

        return EffectiveBinaryDiffusivitySelection(self, fraction)


@dataclass(frozen=True)
class EffectiveBinaryDiffusivitySelection:
    """One explicit endpoint or interior point of a diffusivity interval."""

    interval: EffectiveBinaryDiffusivityInterval
    fraction: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.fraction) or not 0.0 <= self.fraction <= 1.0:
            raise ValueError("effective binary diffusivity fraction must lie in [0, 1]")

    @property
    def value_m2_s(self) -> float:
        """Linearly selected positive effective diffusivity in m2/s."""

        return self.interval.lower_m2_s + self.fraction * (
            self.interval.upper_m2_s - self.interval.lower_m2_s
        )


@dataclass(frozen=True)
class SymmetricBinaryFaceState:
    """Documented centered state used to evaluate ``c_g D y_w y_h``.

    Molar density uses the symmetric arithmetic face mean.  The binary
    composition product uses the symmetric entropy-variable chain-rule mean
    documented in the module header.  Both rules are invariant to swapping
    the left and right cells.  Exact input compositions are retained so a
    pure-component limit can be audited for absence of clipping.
    """

    left_molar_density_mol_m3: float
    right_molar_density_mol_m3: float
    molar_density_mol_m3: float
    left_y_water: float
    left_y_hexane: float
    right_y_water: float
    right_y_hexane: float
    y_water_y_hexane: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        densities = (
            self.left_molar_density_mol_m3,
            self.right_molar_density_mol_m3,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in densities):
            raise ValueError("face molar densities must be positive and finite")
        _validate_binary_composition(self.left_y_water, self.left_y_hexane, "left")
        _validate_binary_composition(self.right_y_water, self.right_y_hexane, "right")

        expected_density = 0.5 * self.left_molar_density_mol_m3 + 0.5 * (
            self.right_molar_density_mol_m3
        )
        expected_product = _binary_chain_rule_product(self.left_y_water, self.right_y_water)
        if (
            not math.isfinite(self.molar_density_mol_m3)
            or self.molar_density_mol_m3 != expected_density
        ):
            raise ValueError("face molar density must equal the symmetric arithmetic mean")
        if not math.isfinite(self.y_water_y_hexane) or self.y_water_y_hexane != expected_product:
            raise ValueError("face composition product must equal the symmetric chain-rule mean")

    @property
    def touches_exact_purity(self) -> bool:
        """Whether either adjacent state is an exact pure component."""

        return any(
            value == 0.0 or value == 1.0 for value in (self.left_y_water, self.right_y_water)
        )

    @property
    def is_uniform_pure_component(self) -> bool:
        """Whether both sides are the same exact pure component."""

        return self.left_y_water == self.right_y_water and (
            self.left_y_water == 0.0 or self.left_y_water == 1.0
        )


@dataclass(frozen=True)
class EvaluatedBinaryMobilityInterval:
    """State-evaluated interval for ``L=c_g D_eff y_w y_h`` in mol/(m s)."""

    diffusivity_interval: EffectiveBinaryDiffusivityInterval
    face: SymmetricBinaryFaceState
    lower_mol_m_s: float
    upper_mol_m_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        scale = self.face.molar_density_mol_m3 * self.face.y_water_y_hexane
        expected = (
            scale * self.diffusivity_interval.lower_m2_s,
            scale * self.diffusivity_interval.upper_m2_s,
        )
        actual = (self.lower_mol_m_s, self.upper_mol_m_s)
        if not all(math.isfinite(value) and value >= 0.0 for value in actual):
            raise ValueError("evaluated binary mobility must be finite and non-negative")
        if actual != expected:
            raise ValueError("evaluated binary mobility endpoints must equal c_g*D_eff*y_w*y_h")

    def select(self, fraction: float) -> "EvaluatedBinaryMobilitySelection":
        """Select the same uncertainty coordinate in ``D_eff`` and ``L``."""

        diffusivity = self.diffusivity_interval.select(fraction)
        return EvaluatedBinaryMobilitySelection(self, diffusivity)


@dataclass(frozen=True)
class EvaluatedBinaryMobilitySelection:
    """One selected diffusivity and its evaluated Onsager face mobility."""

    interval: EvaluatedBinaryMobilityInterval
    diffusivity: EffectiveBinaryDiffusivitySelection
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.diffusivity.interval != self.interval.diffusivity_interval:
            raise ValueError(
                "evaluated mobility and selected diffusivity must use the same interval"
            )

    @property
    def value_mol_m_s(self) -> float:
        """Selected state-dependent Onsager mobility in mol/(m s)."""

        face = self.interval.face
        return face.molar_density_mol_m3 * self.diffusivity.value_m2_s * face.y_water_y_hexane

    def as_coupled_pore_selection(self) -> "coupled_pore.MobilitySelection":
        """Adapt a strictly interior result to ``independent_maxwell_stefan_flux``.

        The existing coupled-pore API intentionally requires a positive
        mobility interval.  At exact purity the product mobility is zero while
        the log-fugacity force is singular, so the split API is not a valid
        numerical representation; use :func:`ideal_binary_fick_limit_flux` for
        that combined limiting oracle.
        """

        if self.interval.face.touches_exact_purity:
            raise PureComponentMobilityLimitError(
                "an exact-purity face has zero mobility and a singular logit force; "
                "evaluate the combined ideal/Fick limit without clipping"
            )

        from dtdc_simulator.core2.particle import coupled_pore

        source = self.interval.diffusivity_interval
        mobility_interval = coupled_pore.MobilityInterval(
            self.interval.lower_mol_m_s,
            self.interval.upper_mol_m_s,
            (
                f"{source.provenance}; state-evaluated "
                "L=c_g*D_wh_eff*y_w*y_h using the symmetric chain-rule face state"
            ),
        )
        return mobility_interval.select(self.diffusivity.fraction)


@dataclass(frozen=True)
class IdealBinaryFickLimitFlux:
    """Combined ideal, isothermal, constant-density Fick limiting oracle."""

    water_diffusive_flux_mol_m2_s: float
    hexane_diffusive_flux_mol_m2_s: float
    diffusivity: EffectiveBinaryDiffusivitySelection
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.water_diffusive_flux_mol_m2_s,
            self.hexane_diffusive_flux_mol_m2_s,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("ideal/Fick limiting fluxes must be finite")
        if self.water_diffusive_flux_mol_m2_s != -self.hexane_diffusive_flux_mol_m2_s:
            raise ValueError("binary ideal/Fick limiting fluxes must be equal and opposite")


def symmetric_binary_face_state(
    left_molar_density_mol_m3: float,
    left_y_water: float,
    left_y_hexane: float,
    right_molar_density_mol_m3: float,
    right_y_water: float,
    right_y_hexane: float,
) -> SymmetricBinaryFaceState:
    """Return the symmetric face state, accepting exact pure compositions.

    Mole fractions are validated but never normalized, floored, or clipped.
    The density interpolation is arithmetic and therefore centered and
    left/right symmetric.  ``y_water_y_hexane`` is the discrete chain-rule
    product rather than a separately averaged storage fraction.
    """

    densities = (left_molar_density_mol_m3, right_molar_density_mol_m3)
    if not all(math.isfinite(value) and value > 0.0 for value in densities):
        raise ValueError("face molar densities must be positive and finite")
    _validate_binary_composition(left_y_water, left_y_hexane, "left")
    _validate_binary_composition(right_y_water, right_y_hexane, "right")

    molar_density = 0.5 * left_molar_density_mol_m3 + 0.5 * right_molar_density_mol_m3
    composition_factor = _binary_chain_rule_product(left_y_water, right_y_water)
    if not math.isfinite(composition_factor) or not 0.0 <= composition_factor <= 0.25:
        raise RuntimeError("binary face composition factor left its analytic bounds")

    return SymmetricBinaryFaceState(
        left_molar_density_mol_m3=left_molar_density_mol_m3,
        right_molar_density_mol_m3=right_molar_density_mol_m3,
        molar_density_mol_m3=molar_density,
        left_y_water=left_y_water,
        left_y_hexane=left_y_hexane,
        right_y_water=right_y_water,
        right_y_hexane=right_y_hexane,
        y_water_y_hexane=composition_factor,
    )


def evaluate_binary_mobility_interval(
    diffusivity: EffectiveBinaryDiffusivityInterval,
    face: SymmetricBinaryFaceState,
) -> EvaluatedBinaryMobilityInterval:
    """Evaluate both positive-``D`` endpoints at one documented face state."""

    scale = face.molar_density_mol_m3 * face.y_water_y_hexane
    lower = scale * diffusivity.lower_m2_s
    upper = scale * diffusivity.upper_m2_s
    if not all(math.isfinite(value) and value >= 0.0 for value in (lower, upper)):
        raise RuntimeError("evaluated binary mobility is negative or non-finite")
    if lower > upper:
        raise RuntimeError("evaluated binary mobility endpoints became unordered")
    return EvaluatedBinaryMobilityInterval(
        diffusivity_interval=diffusivity,
        face=face,
        lower_mol_m_s=lower,
        upper_mol_m_s=upper,
    )


def ideal_binary_fick_limit_flux(
    face: SymmetricBinaryFaceState,
    diffusivity: EffectiveBinaryDiffusivitySelection,
    distance_m: float,
) -> IdealBinaryFickLimitFlux:
    r"""Evaluate the algebraically combined ideal/Fick limit.

    This oracle is deliberately restricted to equal left/right total molar
    density, the condition under which
    ``-c D m_wh Delta(logit(y_w))/Delta r`` reduces exactly to
    ``-c D Delta(y_w)/Delta r``.  It remains finite at exact purity and does
    not introduce a mole-fraction floor.
    """

    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise ValueError("face-state distance must be positive and finite")
    if face.left_molar_density_mol_m3 != face.right_molar_density_mol_m3:
        raise ValueError("ideal/Fick limiting oracle requires exactly constant total molar density")

    gradient = (face.right_y_water - face.left_y_water) / distance_m
    water_flux = -face.molar_density_mol_m3 * diffusivity.value_m2_s * gradient
    if not math.isfinite(water_flux):
        raise RuntimeError("ideal/Fick limiting flux is non-finite")
    return IdealBinaryFickLimitFlux(
        water_diffusive_flux_mol_m2_s=water_flux,
        hexane_diffusive_flux_mol_m2_s=-water_flux,
        diffusivity=diffusivity,
    )


def _validate_binary_composition(
    y_water: float,
    y_hexane: float,
    side: str,
) -> None:
    values = (y_water, y_hexane)
    if not all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in values):
        raise ValueError(f"{side} face mole fractions must lie in [0, 1]")
    if not math.isclose(
        y_water + y_hexane,
        1.0,
        rel_tol=0.0,
        abs_tol=_COMPOSITION_SUM_ABS_TOL,
    ):
        raise ValueError(f"{side} face mole fractions must sum to one")


def _binary_chain_rule_product(left_y_water: float, right_y_water: float) -> float:
    if left_y_water == right_y_water:
        return left_y_water * (1.0 - left_y_water)
    if left_y_water == 0.0 or left_y_water == 1.0 or right_y_water == 0.0 or right_y_water == 1.0:
        return 0.0

    delta_y = right_y_water - left_y_water
    # Direct subtraction of almost-equal logits can round the denominator to
    # zero.  These log1p increments retain the divided-difference limit for
    # adjacent floating-point states without a composition floor.
    delta_logit = math.log1p(delta_y / left_y_water) - math.log1p(-delta_y / (1.0 - left_y_water))
    return delta_y / delta_logit


__all__ = [
    "EffectiveBinaryDiffusivityInterval",
    "EffectiveBinaryDiffusivitySelection",
    "EvaluatedBinaryMobilityInterval",
    "EvaluatedBinaryMobilitySelection",
    "IdealBinaryFickLimitFlux",
    "PureComponentMobilityLimitError",
    "SymmetricBinaryFaceState",
    "evaluate_binary_mobility_interval",
    "ideal_binary_fick_limit_flux",
    "symmetric_binary_face_state",
]
