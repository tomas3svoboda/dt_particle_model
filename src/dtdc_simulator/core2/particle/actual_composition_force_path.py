r"""Actual-``(T, y_h)`` composition forces for a moving dry interface.

This module evaluates only the *composition* part of the two frozen dry-pore
thermodynamic forces.  Along the straight endpoint path

``T(s) = T_L + s (T_R - T_L)``, ``y_h(s) = y_L + s (y_R - y_L)``,

the gas-exchange difference is

``integral_0^1 (partial ln(f_w/f_h)/partial y_h)_(T,P) dy_h/ds ds``.

The ideal-binary contribution is analytic and cancellation-free.  The small
pressure-second-virial contribution uses fixed, deterministic Gauss--Legendre
rules.  The retained-water force follows the unchanged Luikov/cap active set;
cap seams are located and each one-sided branch is integrated separately.
There is no smoothing, clipping, projection, or metastable state storage.

Equal endpoint temperatures are deliberately a separate branch.  It calls
``coupled_pore.isothermal_potential_composition_secant_slopes`` and multiplies
the returned slope by the represented ``delta y_h``.  That preserves the
existing analytic, cancellation-free isothermal operator instead of taking a
difference of nearly equal endpoint potentials.

The result is not a complete non-isothermal chemical-potential difference and
does not invent a Soret coefficient or a standard-state caloric term.  It is a
face reconstruction of the fixed-``(T,P)`` composition force only.  Actual
endpoint equilibrium states are returned so a caller can select mobility and
enthalpy from physical endpoints without constructing common-temperature
virtual states.

The phase guard encloses the *whole* straight path by deterministic adaptive
subdivision.  For the frozen binary activity law, composition monotonicity is
analytic and fixed-``y`` temperature monotonicity follows from

``d ln(a_i)/dT = d ln(phi_i)/dT - d ln(f_i,liq)/dT``.

The certificate uses conservative prospective engineering envelopes for the
negative fixed-composition log-activity temperature derivatives and the
virial composition factor.  The liquid derivative uses each pure authority's
own gas constant; this avoids silently treating ``wa.R*wa.M`` as bit-identical
to the binary-gas ``R``.
Every evaluated actual-path or off-path bounding point must independently
remain inside those envelopes.  Monotone rectangle bounds then prove all
water/hexane inequalities on every accepted leaf; every computed lower/upper
bound is padded outward by 32 representable floats.  A permitted closed hexane
interface endpoint additionally requires a globally negative one-sided path
derivative bound.  Unresolved boxes and cap contacts fail closed.  This is a
continuous segment certificate conditional on the explicitly serialized
engineering envelopes; it is not unconditional interval arithmetic over the
property implementations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Callable

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa


_PHASE_ENCLOSURE_MAX_DEPTH = 28
_PHASE_ENCLOSURE_MIN_FRACTION_WIDTH = 2.0**-28
_CAP_ROOT_MAX_ITERATIONS = 80
_CAP_CONTACT_ULPS = 128.0
_LEFT_HEXANE_CLOSED_MARGIN = 1.0e-10
_DERIVATIVE_CERTIFICATE_MARGIN = 1.0e-12
_ENCLOSURE_ROUNDOFF_ULPS = 32
_TEMPERATURE_ENVELOPE_K = (323.15, 433.0)
_COMPOSITION_ENVELOPE_Y_HEXANE = (0.0, 1.0)
_K_WH_ENVELOPE = (bg.K_WH_MIN, bg.K_WH_MAX)

# Prospective Amendment-20 engineering envelopes.  A deterministic 101 x 5 x
# 101 x 3 audit of the frozen authority over T=[323.15,433] K,
# P=[101,170] kPa, y=[0,1], and k_wh=[0.477,0.511] found
#   d ln(a_w)/dT in [-0.050380,-0.024468] 1/K,
#   d ln(a_h)/dT in [-0.035646,-0.015219] 1/K,
# and V=2 P DeltaB/(R T) in [-0.25727,-0.04492].  The deliberately wider
# values below are the runtime certificate authority.  Inputs outside the
# project pressure range fail closed rather than extrapolating these bounds.
_WATER_LOG_ACTIVITY_DT_MIN_K_INV = -0.060
_WATER_LOG_ACTIVITY_DT_MAX_K_INV = -0.020
_HEXANE_LOG_ACTIVITY_DT_MIN_K_INV = -0.045
_HEXANE_LOG_ACTIVITY_DT_MAX_K_INV = -0.012
_VIRIAL_FACTOR_MIN = -0.30
_VIRIAL_FACTOR_MAX = -0.04

# PART-01 native-cut derivative authority.  These intentionally broad source
# envelopes cover the already admitted A20 project box.  They are not a new
# constitutive law: the local certificate below only encloses the analytic
# identity already evaluated in ``_PathActivities._evaluate_coordinates``.
# The result remains conditional engineering evidence rather than an
# unconditional interval proof of the Helmholtz implementations.
_NATIVE_DERIVATIVE_SOURCE_ENVELOPE_ID = "PART-01/native-linear-path-envelopes/v1"
_NATIVE_DERIVATIVE_SOURCE_AUTHORITY_IDS = (
    "hexane.state_Tp(T,P,'liquid').{h,cp}|Span-Wagner-2003|"
    "docs/GT_PS2_PURE_HEXANE_PROPERTY_AUTHORITY.md",
    "hexane.{h_ideal,cp0,second_virial,second_virial_dT}|"
    "Span-Wagner-2003+Jaeschke-Schley-1995",
    "water.{second_virial,second_virial_dT}|IAPWS-95",
    "binary_gas.cross_second_virial|PHY-053|Plyasunov-Shock-Tsonopoulos",
)
_HEXANE_LIQUID_CP_ENVELOPE_J_MOL_K = (180.0, 300.0)
_HEXANE_IDEAL_CP_ENVELOPE_J_MOL_K = (140.0, 220.0)
_WATER_SECOND_VIRIAL_ENVELOPE_M3_MOL = (-9.0e-4, -2.0e-4)
_CROSS_SECOND_VIRIAL_ENVELOPE_M3_MOL = (-2.1e-4, -5.0e-5)
_HEXANE_SECOND_VIRIAL_ENVELOPE_M3_MOL = (-1.6e-3, -6.5e-4)
_WATER_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K = (1.5e-6, 1.3e-5)
_CROSS_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K = (5.0e-7, 1.7e-6)
_HEXANE_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K = (3.5e-6, 1.4e-5)

# Fixed rules avoid generating nodes inside every nonlinear residual call.
_GAUSS_LEGENDRE_8 = (
    (-0.9602898564975363, 0.1012285362903763),
    (-0.7966664774136267, 0.2223810344533745),
    (-0.5255324099163290, 0.3137066458778873),
    (-0.1834346424956498, 0.3626837833783620),
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
)
_GAUSS_LEGENDRE_16 = (
    (-0.9894009349916499, 0.027152459411754095),
    (-0.9445750230732326, 0.06225352393864789),
    (-0.8656312023878318, 0.09515851168249278),
    (-0.7554044083550030, 0.12462897125553387),
    (-0.6178762444026438, 0.14959598881657673),
    (-0.4580167776572274, 0.16915651939500254),
    (-0.2816035507792589, 0.1826034150449236),
    (-0.09501250983763744, 0.1894506104550685),
    (0.09501250983763744, 0.1894506104550685),
    (0.2816035507792589, 0.1826034150449236),
    (0.4580167776572274, 0.16915651939500254),
    (0.6178762444026438, 0.14959598881657673),
    (0.7554044083550030, 0.12462897125553387),
    (0.8656312023878318, 0.09515851168249278),
    (0.9445750230732326, 0.06225352393864789),
    (0.9894009349916499, 0.027152459411754095),
)


class CompositionForcePathMode(str, Enum):
    """Numerical branch used for one path reconstruction."""

    ISOTHERMAL_EXACT_SECANT = "isothermal_exact_secant"
    NONISOTHERMAL_ACTUAL_PATH = "nonisothermal_actual_path"


class ClosedHexaneDerivativeAuthority(str, Enum):
    """Certificate authority for a closed interface on the linear path.

    The historical broad envelope remains the public default.  Native cut-face
    consumers opt into the local source-envelope enclosure explicitly, so this
    packet cannot silently change archived A20/A22 reference decisions.
    """

    GLOBAL_PROSPECTIVE_ENVELOPE = "global_prospective_envelope"
    NATIVE_LOCAL_ANALYTIC_ENVELOPE = "native_local_analytic_envelope"


class ActualCompositionForceConvergenceError(RuntimeError):
    """A fixed quadrature pair did not agree at the declared tolerance."""


@dataclass(frozen=True)
class GaussLegendreAudit:
    """Deterministic nested-order comparison for one virial contribution."""

    quantity: str
    comparison_order: int
    final_order: int
    comparison_estimate: float
    final_estimate: float
    absolute_change: float
    acceptance_tolerance: float
    converged: bool
    fixed_nodes_and_weights: bool = field(default=True, init=False)

    def __post_init__(self) -> None:
        values = (
            self.comparison_estimate,
            self.final_estimate,
            self.absolute_change,
            self.acceptance_tolerance,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("quadrature audit values must be finite")
        if self.comparison_order <= 0 or self.final_order <= self.comparison_order:
            raise ValueError("quadrature audit orders must be positive and increasing")
        if self.absolute_change < 0.0 or self.acceptance_tolerance <= 0.0:
            raise ValueError("quadrature error and tolerance must be admissible")
        if self.converged != (self.absolute_change <= self.acceptance_tolerance):
            raise ValueError("quadrature convergence flag disagrees with its audit")


@dataclass(frozen=True)
class CapSeamAudit:
    """One deterministically bracketed retained-cap/Luikov path intersection."""

    fraction: float
    water_activity_minus_cap: float
    final_bracket_width: float
    isolation_depth: int
    bisection_iterations: int

    def __post_init__(self) -> None:
        values = (
            self.fraction,
            self.water_activity_minus_cap,
            self.final_bracket_width,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("cap-seam audit values must be finite")
        if not 0.0 <= self.fraction <= 1.0 or self.final_bracket_width < 0.0:
            raise ValueError("cap-seam audit is outside the path")
        if self.isolation_depth < 0 or self.bisection_iterations < 0:
            raise ValueError("cap-seam iteration counts must be nonnegative")


@dataclass(frozen=True)
class ActualPathNodeAudit:
    """One actual-path point consumed by enclosure, seam, or force work."""

    fraction: float
    roles: tuple[str, ...]
    temperature_k: float
    y_hexane: float
    water_activity: float
    hexane_activity: float
    virial_factor: float
    water_log_activity_temperature_derivative_k_inv: float
    hexane_log_activity_temperature_derivative_k_inv: float

    def __post_init__(self) -> None:
        values = (
            self.fraction,
            self.temperature_k,
            self.y_hexane,
            self.water_activity,
            self.hexane_activity,
            self.virial_factor,
            self.water_log_activity_temperature_derivative_k_inv,
            self.hexane_log_activity_temperature_derivative_k_inv,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("actual-path node audit values must be finite")
        if not 0.0 <= self.fraction <= 1.0 or not self.roles:
            raise ValueError("actual-path node needs a valid fraction and role")


@dataclass(frozen=True)
class ActivityEnclosureLeafAudit:
    """One whole-subsegment monotone activity enclosure."""

    start_fraction: float
    end_fraction: float
    depth: int
    water_activity_lower_bound: float
    water_activity_upper_bound: float
    hexane_activity_lower_bound: float
    hexane_activity_upper_bound: float
    hexane_upper_bound_replaced_by_closed_endpoint_derivative: bool

    def __post_init__(self) -> None:
        values = (
            self.start_fraction,
            self.end_fraction,
            self.water_activity_lower_bound,
            self.water_activity_upper_bound,
            self.hexane_activity_lower_bound,
            self.hexane_activity_upper_bound,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("activity-enclosure leaf values must be finite")
        if not 0.0 <= self.start_fraction < self.end_fraction <= 1.0:
            raise ValueError("activity-enclosure leaf must be a nonempty path interval")
        if self.depth < 0:
            raise ValueError("activity-enclosure depth must be nonnegative")
        if self.water_activity_lower_bound > self.water_activity_upper_bound or (
            self.hexane_activity_lower_bound > self.hexane_activity_upper_bound
        ):
            raise ValueError("activity-enclosure bounds must be ordered")


@dataclass(frozen=True)
class NativeClosedHexaneDerivativeCertificate:
    """Analytic full-segment inward derivative enclosure for native cut faces.

    All stored intermediate intervals are recomputed in ``__post_init__`` from
    the physical endpoints and the frozen source envelopes.  A consumer cannot
    replace the final sign bound, its property dependencies, or its provenance
    independently and still construct a certificate.
    """

    left_temperature_k: float
    left_y_hexane: float
    right_temperature_k: float
    right_y_hexane: float
    pressure_pa: float
    k_wh: float
    hexane_liquid_cp_samples_j_mol_k: tuple[float, float, float]
    hexane_ideal_cp_samples_j_mol_k: tuple[float, float, float]
    hexane_liquid_enthalpy_bounds_j_mol: tuple[float, float]
    hexane_ideal_enthalpy_bounds_j_mol: tuple[float, float]
    hexane_partial_residual_enthalpy_bounds_j_mol: tuple[float, float]
    hexane_log_activity_temperature_derivative_bounds_k_inv: tuple[float, float]
    hexane_log_activity_composition_derivative_bounds: tuple[float, float]
    hexane_log_activity_path_derivative_bounds: tuple[float, float]
    source_envelope_id: str
    source_authority_ids: tuple[str, ...]
    hexane_liquid_cp_envelope_j_mol_k: tuple[float, float]
    hexane_ideal_cp_envelope_j_mol_k: tuple[float, float]
    pure_and_cross_second_virial_envelopes_m3_mol: tuple[
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
    ]
    pure_and_cross_second_virial_dt_envelopes_m3_mol_k: tuple[
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
    ]
    analytic_identity: str = field(
        default=(
            "dln(a_h)/ds=dy_h/ds*(1/y_h+V*(1-y_h))"
            "+dT/ds*(h_h_liq-h_h_ideal-h_h_res)/(R_h*T**2)"
        ),
        init=False,
    )
    full_segment_enclosed: bool = field(default=True, init=False)
    finite_difference_used_for_acceptance: bool = field(default=False, init=False)
    conditional_on_frozen_source_envelopes: bool = field(
        default=True,
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)
    industrially_feasible: bool = field(default=False, init=False)
    plant_predictive: bool = field(default=False, init=False)
    production_wired: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.left_temperature_k,
            self.left_y_hexane,
            self.right_temperature_k,
            self.right_y_hexane,
            self.pressure_pa,
            self.k_wh,
            *self.hexane_liquid_cp_samples_j_mol_k,
            *self.hexane_ideal_cp_samples_j_mol_k,
            *self.hexane_liquid_enthalpy_bounds_j_mol,
            *self.hexane_ideal_enthalpy_bounds_j_mol,
            *self.hexane_partial_residual_enthalpy_bounds_j_mol,
            *self.hexane_log_activity_temperature_derivative_bounds_k_inv,
            *self.hexane_log_activity_composition_derivative_bounds,
            *self.hexane_log_activity_path_derivative_bounds,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("native closed-hexane certificate values must be finite")
        if self.pressure_pa <= 0.0 or not (
            0.0 < self.left_y_hexane < 1.0
            and 0.0 < self.right_y_hexane < 1.0
        ):
            raise ValueError("native closed-hexane certificate endpoints are invalid")
        if not _K_WH_ENVELOPE[0] <= self.k_wh <= _K_WH_ENVELOPE[1]:
            raise ValueError("native closed-hexane k_wh is outside its source envelope")
        intervals = (
            self.hexane_liquid_enthalpy_bounds_j_mol,
            self.hexane_ideal_enthalpy_bounds_j_mol,
            self.hexane_partial_residual_enthalpy_bounds_j_mol,
            self.hexane_log_activity_temperature_derivative_bounds_k_inv,
            self.hexane_log_activity_composition_derivative_bounds,
            self.hexane_log_activity_path_derivative_bounds,
        )
        if any(lower > upper for lower, upper in intervals):
            raise ValueError("native closed-hexane certificate interval is unordered")
        if self.source_envelope_id != _NATIVE_DERIVATIVE_SOURCE_ENVELOPE_ID:
            raise ValueError("native closed-hexane source-envelope provenance changed")
        if self.source_authority_ids != _NATIVE_DERIVATIVE_SOURCE_AUTHORITY_IDS:
            raise ValueError("native closed-hexane source authority changed")
        if self.hexane_liquid_cp_envelope_j_mol_k != (
            _HEXANE_LIQUID_CP_ENVELOPE_J_MOL_K
        ) or self.hexane_ideal_cp_envelope_j_mol_k != (
            _HEXANE_IDEAL_CP_ENVELOPE_J_MOL_K
        ):
            raise ValueError("native closed-hexane heat-capacity envelope changed")
        if self.pure_and_cross_second_virial_envelopes_m3_mol != (
            _WATER_SECOND_VIRIAL_ENVELOPE_M3_MOL,
            _CROSS_SECOND_VIRIAL_ENVELOPE_M3_MOL,
            _HEXANE_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        ) or self.pure_and_cross_second_virial_dt_envelopes_m3_mol_k != (
            _WATER_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
            _CROSS_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
            _HEXANE_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
        ):
            raise ValueError("native closed-hexane virial source envelope changed")

        (
            liquid_samples,
            ideal_samples,
            liquid_bounds,
            ideal_bounds,
        ) = _native_hexane_enthalpy_bounds(
            min(self.left_temperature_k, self.right_temperature_k),
            max(self.left_temperature_k, self.right_temperature_k),
            self.pressure_pa,
        )
        expected_residual = _native_hexane_partial_residual_enthalpy_bounds(
            min(self.left_temperature_k, self.right_temperature_k),
            max(self.left_temperature_k, self.right_temperature_k),
            min(self.left_y_hexane, self.right_y_hexane),
            max(self.left_y_hexane, self.right_y_hexane),
            self.pressure_pa,
            self.k_wh,
        )
        expected_temperature = _native_hexane_temperature_derivative_bounds(
            min(self.left_temperature_k, self.right_temperature_k),
            max(self.left_temperature_k, self.right_temperature_k),
            liquid_bounds,
            ideal_bounds,
            expected_residual,
        )
        expected_composition = _composition_derivative_intervals(
            min(self.left_y_hexane, self.right_y_hexane),
            max(self.left_y_hexane, self.right_y_hexane),
        )[1]
        expected_path = _path_derivative_bounds(
            math.fsum((self.right_y_hexane, -self.left_y_hexane)),
            math.fsum((self.right_temperature_k, -self.left_temperature_k)),
            expected_composition,
            expected_temperature,
        )
        actual = (
            self.hexane_liquid_cp_samples_j_mol_k,
            self.hexane_ideal_cp_samples_j_mol_k,
            self.hexane_liquid_enthalpy_bounds_j_mol,
            self.hexane_ideal_enthalpy_bounds_j_mol,
            self.hexane_partial_residual_enthalpy_bounds_j_mol,
            self.hexane_log_activity_temperature_derivative_bounds_k_inv,
            self.hexane_log_activity_composition_derivative_bounds,
            self.hexane_log_activity_path_derivative_bounds,
        )
        expected = (
            liquid_samples,
            ideal_samples,
            liquid_bounds,
            ideal_bounds,
            expected_residual,
            expected_temperature,
            expected_composition,
            expected_path,
        )
        if actual != expected:
            raise ValueError(
                "native closed-hexane certificate changed a source dependency"
            )
        if expected_temperature[0] < _HEXANE_LOG_ACTIVITY_DT_MIN_K_INV or (
            expected_temperature[1] > _HEXANE_LOG_ACTIVITY_DT_MAX_K_INV
        ):
            raise ValueError(
                "native closed-hexane local derivative left the parent envelope"
            )
        if expected_path[1] >= -_DERIVATIVE_CERTIFICATE_MARGIN:
            raise ValueError(
                "native closed-hexane certificate lacks a full-segment inward bound"
            )


@dataclass(frozen=True)
class ActualPathPhaseAudit:
    """Whole-segment activity enclosure and retained-seam audit."""

    evaluated_actual_path_nodes: tuple[ActualPathNodeAudit, ...]
    enclosure_leaves: tuple[ActivityEnclosureLeafAudit, ...]
    cap_seams: tuple[CapSeamAudit, ...]
    minimum_y_hexane: float
    minimum_y_water: float
    minimum_water_activity_minus_luikov_ref: float
    minimum_one_minus_water_activity: float
    minimum_one_minus_hexane_activity: float
    water_activity_at_ref: float
    water_activity_at_cap: float
    actual_path_activity_evaluation_count: int
    off_path_bound_evaluation_count: int
    maximum_enclosure_depth: int
    maximum_cap_isolation_depth: int
    left_hexane_endpoint_treated_as_closed: bool
    closed_endpoint_hexane_log_activity_derivative_upper_bound: float | None
    closed_hexane_derivative_authority: ClosedHexaneDerivativeAuthority
    native_closed_hexane_derivative_certificate: (
        NativeClosedHexaneDerivativeCertificate | None
    )
    virial_factor_envelope: tuple[float, float]
    water_log_activity_temperature_derivative_envelope_k_inv: tuple[float, float]
    hexane_log_activity_temperature_derivative_envelope_k_inv: tuple[float, float]
    project_pressure_envelope_pa: tuple[float, float]
    temperature_envelope_k: tuple[float, float]
    composition_envelope_y_hexane: tuple[float, float]
    k_wh_envelope: tuple[float, float]
    deterministic_guard_passed: bool
    no_clipping_or_projection_used: bool = field(default=True, init=False)
    no_metastable_storage_used: bool = field(default=True, init=False)
    intermediate_storage_states_constructed: bool = field(default=False, init=False)
    endpoint_full_equilibrium_state_count: int = field(default=2, init=False)
    off_path_points_used_only_for_activity_bounds: bool = field(default=True, init=False)
    unresolved_boxes_and_cap_contacts_rejected: bool = field(default=True, init=False)
    continuous_segment_certificate: bool = field(default=True, init=False)
    conditional_on_frozen_engineering_envelopes: bool = field(default=True, init=False)
    unconditional_property_interval_arithmetic_proof: bool = field(
        default=False,
        init=False,
    )
    temperature_monotonicity_identity: str = field(
        default="dln(a_i)/dT=dln(phi_i)/dT-dln(f_i_liq)/dT",
        init=False,
    )
    maximum_allowed_enclosure_depth: int = field(
        default=_PHASE_ENCLOSURE_MAX_DEPTH,
        init=False,
    )
    minimum_allowed_fraction_width: float = field(
        default=_PHASE_ENCLOSURE_MIN_FRACTION_WIDTH,
        init=False,
    )
    outward_roundoff_padding_ulps: int = field(
        default=_ENCLOSURE_ROUNDOFF_ULPS,
        init=False,
    )

    def __post_init__(self) -> None:
        if not self.evaluated_actual_path_nodes or not self.enclosure_leaves:
            raise ValueError("path phase audit needs evaluated nodes and enclosure leaves")
        fractions = tuple(node.fraction for node in self.evaluated_actual_path_nodes)
        if fractions[0] != 0.0 or fractions[-1] != 1.0:
            raise ValueError("path phase audit must include both endpoints")
        if any(right <= left for left, right in zip(fractions, fractions[1:])):
            raise ValueError("path phase audit fractions must be strictly ordered")
        if self.enclosure_leaves[0].start_fraction != 0.0 or (
            self.enclosure_leaves[-1].end_fraction != 1.0
        ):
            raise ValueError("activity-enclosure leaves must cover both path endpoints")
        if any(
            left.end_fraction != right.start_fraction
            for left, right in zip(
                self.enclosure_leaves,
                self.enclosure_leaves[1:],
            )
        ):
            raise ValueError("activity-enclosure leaves must form a contiguous partition")
        margins = (
            self.minimum_y_hexane,
            self.minimum_y_water,
            self.minimum_water_activity_minus_luikov_ref,
            self.minimum_one_minus_water_activity,
            self.minimum_one_minus_hexane_activity,
            self.water_activity_at_ref,
            self.water_activity_at_cap,
        )
        if not all(math.isfinite(value) for value in margins):
            raise ValueError("path phase margins must be finite")
        if min(margins) < 0.0 or self.actual_path_activity_evaluation_count != len(
            self.evaluated_actual_path_nodes
        ):
            raise ValueError("path phase audit contains an invalid margin or count")
        if (
            self.off_path_bound_evaluation_count < 0
            or self.maximum_enclosure_depth < 0
            or self.maximum_cap_isolation_depth < 0
        ):
            raise ValueError("path phase enclosure counts must be nonnegative")
        if max(
            self.maximum_enclosure_depth,
            self.maximum_cap_isolation_depth,
        ) > _PHASE_ENCLOSURE_MAX_DEPTH:
            raise ValueError("path phase enclosure exceeded its frozen depth")
        if self.virial_factor_envelope != (
            _VIRIAL_FACTOR_MIN,
            _VIRIAL_FACTOR_MAX,
        ) or self.water_log_activity_temperature_derivative_envelope_k_inv != (
            _WATER_LOG_ACTIVITY_DT_MIN_K_INV,
            _WATER_LOG_ACTIVITY_DT_MAX_K_INV,
        ) or self.hexane_log_activity_temperature_derivative_envelope_k_inv != (
            _HEXANE_LOG_ACTIVITY_DT_MIN_K_INV,
            _HEXANE_LOG_ACTIVITY_DT_MAX_K_INV,
        ):
            raise ValueError("path phase audit does not serialize the frozen envelopes")
        if (
            self.project_pressure_envelope_pa
            != (bg.PROJECT_PRESSURE_MIN_PA, bg.PROJECT_PRESSURE_MAX_PA)
            or self.temperature_envelope_k != _TEMPERATURE_ENVELOPE_K
            or self.composition_envelope_y_hexane
            != _COMPOSITION_ENVELOPE_Y_HEXANE
            or self.k_wh_envelope != _K_WH_ENVELOPE
        ):
            raise ValueError("path phase audit has an incorrect independent-variable envelope")
        if any(
            not (
                _VIRIAL_FACTOR_MIN <= node.virial_factor <= _VIRIAL_FACTOR_MAX
                and _WATER_LOG_ACTIVITY_DT_MIN_K_INV
                <= node.water_log_activity_temperature_derivative_k_inv
                <= _WATER_LOG_ACTIVITY_DT_MAX_K_INV
                and _HEXANE_LOG_ACTIVITY_DT_MIN_K_INV
                <= node.hexane_log_activity_temperature_derivative_k_inv
                <= _HEXANE_LOG_ACTIVITY_DT_MAX_K_INV
            )
            for node in self.evaluated_actual_path_nodes
        ):
            raise ValueError("actual-path node left a serialized derivative envelope")
        if any(
            not (
                leaf.water_activity_lower_bound > self.water_activity_at_ref
                and leaf.water_activity_upper_bound < 1.0
                and leaf.hexane_activity_lower_bound > 0.0
                and (
                    leaf.hexane_activity_upper_bound < 1.0
                    or leaf.hexane_upper_bound_replaced_by_closed_endpoint_derivative
                )
            )
            for leaf in self.enclosure_leaves
        ):
            raise ValueError("activity-enclosure leaf does not certify its whole segment")
        if self.left_hexane_endpoint_treated_as_closed:
            derivative = self.closed_endpoint_hexane_log_activity_derivative_upper_bound
            if derivative is None or derivative >= -_DERIVATIVE_CERTIFICATE_MARGIN:
                raise ValueError("closed interface needs a strictly inward derivative bound")
            if (
                self.closed_hexane_derivative_authority
                is ClosedHexaneDerivativeAuthority.NATIVE_LOCAL_ANALYTIC_ENVELOPE
            ):
                certificate = self.native_closed_hexane_derivative_certificate
                if certificate is None or derivative != (
                    certificate.hexane_log_activity_path_derivative_bounds[1]
                ):
                    raise ValueError(
                        "native closed interface needs its exact analytic certificate"
                    )
                left = self.evaluated_actual_path_nodes[0]
                right = self.evaluated_actual_path_nodes[-1]
                if (
                    certificate.left_temperature_k != left.temperature_k
                    or certificate.left_y_hexane != left.y_hexane
                    or certificate.right_temperature_k != right.temperature_k
                    or certificate.right_y_hexane != right.y_hexane
                ):
                    raise ValueError(
                        "native closed-interface certificate belongs to another path"
                    )
            elif self.native_closed_hexane_derivative_certificate is not None:
                raise ValueError(
                    "global closed-interface authority cannot carry a native certificate"
                )
        elif self.closed_endpoint_hexane_log_activity_derivative_upper_bound is not None:
            raise ValueError("strict interface cannot carry a closed-endpoint derivative")
        elif self.native_closed_hexane_derivative_certificate is not None:
            raise ValueError("strict interface cannot carry a native derivative certificate")
        if not isinstance(
            self.closed_hexane_derivative_authority,
            ClosedHexaneDerivativeAuthority,
        ):
            raise TypeError("path phase audit needs a named derivative authority")
        if not self.deterministic_guard_passed:
            raise ValueError("a failed path guard cannot construct a result")


@dataclass(frozen=True)
class RetainedPathSegmentAudit:
    """One unsmoothed, constant-active-set retained-force path segment."""

    start_fraction: float
    end_fraction: float
    active_set: cp.RetainedWaterActiveSet
    force_difference_right_minus_left: float
    ideal_mixture_contribution: float
    virial_contribution: float
    quadrature: GaussLegendreAudit | None

    def __post_init__(self) -> None:
        if not 0.0 <= self.start_fraction < self.end_fraction <= 1.0:
            raise ValueError("retained-force segment must be a nonempty path interval")
        if not isinstance(self.active_set, cp.RetainedWaterActiveSet):
            raise TypeError("retained-force segment needs a named active set")
        values = (
            self.force_difference_right_minus_left,
            self.ideal_mixture_contribution,
            self.virial_contribution,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("retained-force segment contributions must be finite")
        if self.active_set is cp.RetainedWaterActiveSet.RETAINED_CAP:
            if values != (0.0, 0.0, 0.0) or self.quadrature is not None:
                raise ValueError("retained-cap segment force must be exactly zero")


@dataclass(frozen=True)
class ActualCompositionForcePath:
    """Composition-force reconstruction and physical endpoint states."""

    mode: CompositionForcePathMode
    left_state: cp.EquilibriumPoreState
    right_state: cp.EquilibriumPoreState
    distance_m: float
    delta_temperature_k: float
    delta_y_hexane: float
    gas_force_difference_right_minus_left: float
    retained_force_difference_right_minus_left: float
    gas_force_gradient_m_inv: float
    retained_force_gradient_m_inv: float
    gas_force_secant_per_y_hexane: float | None
    retained_force_secant_per_y_hexane: float | None
    gas_ideal_mixture_contribution: float | None
    gas_virial_contribution: float | None
    gas_quadrature: GaussLegendreAudit | None
    retained_segments: tuple[RetainedPathSegmentAudit, ...]
    phase_audit: ActualPathPhaseAudit
    evaluation_followed_caller_direction: bool = field(default=True, init=False)
    entropy_ready: bool = field(default=True, init=False)
    actual_endpoint_states_returned: bool = field(default=True, init=False)
    common_temperature_virtual_force_or_storage_states_constructed: bool = field(
        default=False,
        init=False,
    )
    off_path_activity_bound_points_constructed: bool = field(default=True, init=False)
    complete_nonisothermal_chemical_potential: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.mode, CompositionForcePathMode):
            raise TypeError("composition-force path needs a named numerical mode")
        if not math.isfinite(self.distance_m) or self.distance_m <= 0.0:
            raise ValueError("composition-force path distance must be positive and finite")
        values = (
            self.delta_temperature_k,
            self.delta_y_hexane,
            self.gas_force_difference_right_minus_left,
            self.retained_force_difference_right_minus_left,
            self.gas_force_gradient_m_inv,
            self.retained_force_gradient_m_inv,
        )
        optional = (
            self.gas_force_secant_per_y_hexane,
            self.retained_force_secant_per_y_hexane,
            self.gas_ideal_mixture_contribution,
            self.gas_virial_contribution,
        )
        if not all(math.isfinite(value) for value in values) or not all(
            value is None or math.isfinite(value) for value in optional
        ):
            raise ValueError("composition-force result must be finite")
        if self.left_state.pressure_pa != self.right_state.pressure_pa:
            raise ValueError("composition-force endpoint pressures must be identical")
        if not self.retained_segments:
            raise ValueError("composition-force result needs a retained path partition")


@dataclass(frozen=True)
class _ActivityPoint:
    fraction: float
    temperature_k: float
    y_hexane: float
    water_activity: float
    hexane_activity: float
    virial_factor: float
    water_log_activity_temperature_derivative_k_inv: float
    hexane_log_activity_temperature_derivative_k_inv: float


@dataclass(frozen=True)
class _ActivityBox:
    start_fraction: float
    end_fraction: float
    depth: int
    water_activity_lower_bound: float
    water_activity_upper_bound: float
    hexane_activity_lower_bound: float
    hexane_activity_upper_bound: float
    water_path_derivative_bounds: tuple[float, float]
    hexane_path_derivative_bounds: tuple[float, float]


@dataclass(frozen=True)
class _PathCertificate:
    enclosure_leaves: tuple[ActivityEnclosureLeafAudit, ...]
    cap_seams: tuple[CapSeamAudit, ...]
    maximum_cap_isolation_depth: int
    left_hexane_endpoint_treated_as_closed: bool
    closed_endpoint_hexane_log_activity_derivative_upper_bound: float | None
    closed_hexane_derivative_authority: ClosedHexaneDerivativeAuthority
    native_closed_hexane_derivative_certificate: (
        NativeClosedHexaneDerivativeCertificate | None
    )


@dataclass(frozen=True)
class _OrderedEvaluation:
    mode: CompositionForcePathMode
    gas_difference: float
    retained_difference: float
    gas_ideal: float | None
    gas_virial: float | None
    gas_quadrature: GaussLegendreAudit | None
    retained_segments: tuple[RetainedPathSegmentAudit, ...]
    phase_audit: ActualPathPhaseAudit


class _PathActivities:
    """Per-call activity/envelope cache shared by all path consumers."""

    def __init__(
        self,
        left_temperature_k: float,
        left_y_hexane: float,
        right_temperature_k: float,
        right_y_hexane: float,
        pressure_pa: float,
        pore: cp.CoupledPoreParams,
    ) -> None:
        if not bg.PROJECT_PRESSURE_MIN_PA <= pressure_pa <= bg.PROJECT_PRESSURE_MAX_PA:
            raise cp.CoupledPoreTopologyError(
                "actual-path enclosure is authorized only inside the frozen project "
                "pressure envelope"
            )
        if not _K_WH_ENVELOPE[0] <= pore.k_wh <= _K_WH_ENVELOPE[1]:
            raise cp.CoupledPoreTopologyError(
                "actual-path k_wh is outside the frozen enclosure evidence"
            )
        self.left_temperature_k = left_temperature_k
        self.left_y_hexane = left_y_hexane
        self.delta_temperature_k = right_temperature_k - left_temperature_k
        self.delta_y_hexane = right_y_hexane - left_y_hexane
        self.pressure_pa = pressure_pa
        self.pore = pore
        self.activity_at_ref = sp.water_activity(pore.luikov.W_ref, pore.luikov)
        self.activity_at_cap = sp.water_activity(pore.luikov.W_cap, pore.luikov)
        self._points: dict[float, _ActivityPoint] = {}
        self._roles: dict[float, set[str]] = {}
        self._bound_points: dict[tuple[float, float], tuple[float, ...]] = {}

    def _evaluate_coordinates(
        self,
        temperature_k: float,
        y_hexane: float,
    ) -> tuple[float, float, float, float, float]:
        if not _TEMPERATURE_ENVELOPE_K[0] <= temperature_k <= (
            _TEMPERATURE_ENVELOPE_K[1]
        ):
            raise cp.CoupledPoreTopologyError(
                "actual-path temperature left the frozen enclosure evidence"
            )
        water_activity, hexane_activity = cp.binary_gas_component_activities(
            temperature_k,
            self.pressure_pa,
            y_hexane,
            self.pore,
        )
        mixture = bg.state(
            temperature_k,
            self.pressure_pa,
            1.0 - y_hexane,
            y_hexane,
            k_wh=self.pore.k_wh,
        )
        delta_virial = mixture.B_ww - 2.0 * mixture.B_wh + mixture.B_hh
        virial_factor = (
            2.0 * self.pressure_pa * delta_virial / (bg.R * temperature_k)
        )
        liquid_water = wa.state_Tp(temperature_k, self.pressure_pa, "liquid")
        liquid_hexane = hx.state_Tp(temperature_k, self.pressure_pa, "liquid")
        # At fixed (P,y), d ln(phi_i)/dT = -h_i,res/(R_mix T^2), while
        # d ln(f_i,liq)/dT = (h_i,ig-h_i,liq)/(R_i T^2).  Keep the pure
        # authority's R_i explicit: the water and binary-gas molar constants
        # are close but not bit-identical.
        water_liquid_log_fugacity_dt = (
            bg.water_ideal_enthalpy_molar(temperature_k)
            - liquid_water.h_mass * wa.M
        ) / (wa.R * wa.M * temperature_k**2)
        hexane_liquid_log_fugacity_dt = (
            hx.h_ideal(temperature_k) - liquid_hexane.h
        ) / (hx.R * temperature_k**2)
        water_dt = (
            -mixture.partial_residual_enthalpy_water_molar
            / (bg.R * temperature_k**2)
            - water_liquid_log_fugacity_dt
        )
        hexane_dt = (
            -mixture.partial_residual_enthalpy_hexane_molar
            / (bg.R * temperature_k**2)
            - hexane_liquid_log_fugacity_dt
        )
        values = (
            water_activity,
            hexane_activity,
            virial_factor,
            water_dt,
            hexane_dt,
        )
        if not all(math.isfinite(value) for value in values):
            raise cp.CoupledPoreTopologyError(
                "actual-path activity envelope produced a non-finite value"
            )
        if not _VIRIAL_FACTOR_MIN <= virial_factor <= _VIRIAL_FACTOR_MAX:
            raise cp.CoupledPoreTopologyError(
                "actual-path virial factor left its prospective enclosure; fail closed"
            )
        if not (
            _WATER_LOG_ACTIVITY_DT_MIN_K_INV
            <= water_dt
            <= _WATER_LOG_ACTIVITY_DT_MAX_K_INV
            and _HEXANE_LOG_ACTIVITY_DT_MIN_K_INV
            <= hexane_dt
            <= _HEXANE_LOG_ACTIVITY_DT_MAX_K_INV
        ):
            raise cp.CoupledPoreTopologyError(
                "fixed-composition activity temperature derivative left its "
                "prospective enclosure; fail closed"
            )
        return values

    def point(self, fraction: float, *, role: str) -> _ActivityPoint:
        if not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
            raise ValueError("actual-path fraction must lie in [0,1]")
        if not role:
            raise ValueError("actual-path activity evaluation needs an audit role")
        cached = self._points.get(fraction)
        if cached is not None:
            self._roles[fraction].add(role)
            return cached
        temperature_k = math.fsum(
            (self.left_temperature_k, fraction * self.delta_temperature_k)
        )
        y_hexane = math.fsum((self.left_y_hexane, fraction * self.delta_y_hexane))
        if not 0.0 < y_hexane < 1.0:
            raise cp.CoupledPoreTopologyError(
                "actual composition-force path left strict interior composition; "
                "do not clip"
            )
        (
            water_activity,
            hexane_activity,
            virial_factor,
            water_dt,
            hexane_dt,
        ) = self._evaluate_coordinates(temperature_k, y_hexane)
        water_strict = (
            self.activity_at_ref < water_activity < 1.0
        )
        hexane_valid = (
            hexane_activity <= 1.0
            if fraction == 0.0
            else hexane_activity < 1.0
        )
        if not water_strict or not hexane_valid or min(
            water_activity,
            hexane_activity,
        ) <= 0.0:
            raise cp.CoupledPoreTopologyError(
                "actual (T,y) composition-force path leaves the caller-oriented "
                "gas-only domain; only the left interface may close hexane "
                "saturation, and water/Luikov margins remain strict"
            )
        point = _ActivityPoint(
            fraction,
            temperature_k,
            y_hexane,
            water_activity,
            hexane_activity,
            virial_factor,
            water_dt,
            hexane_dt,
        )
        self._points[fraction] = point
        self._roles[fraction] = {role}
        return point

    def bound_point(self, temperature_k: float, y_hexane: float) -> tuple[float, ...]:
        key = (temperature_k, y_hexane)
        cached = self._bound_points.get(key)
        if cached is not None:
            return cached
        if not 0.0 < y_hexane < 1.0:
            raise cp.CoupledPoreTopologyError(
                "activity enclosure corner left strict interior composition"
            )
        evaluated = self._evaluate_coordinates(temperature_k, y_hexane)
        self._bound_points[key] = evaluated
        return evaluated

    def node_audits(self) -> tuple[ActualPathNodeAudit, ...]:
        return tuple(
            ActualPathNodeAudit(
                fraction=point.fraction,
                roles=tuple(sorted(self._roles[fraction])),
                temperature_k=point.temperature_k,
                y_hexane=point.y_hexane,
                water_activity=point.water_activity,
                hexane_activity=point.hexane_activity,
                virial_factor=point.virial_factor,
                water_log_activity_temperature_derivative_k_inv=(
                    point.water_log_activity_temperature_derivative_k_inv
                ),
                hexane_log_activity_temperature_derivative_k_inv=(
                    point.hexane_log_activity_temperature_derivative_k_inv
                ),
            )
            for fraction, point in sorted(self._points.items())
        )

    @property
    def off_path_evaluation_count(self) -> int:
        return len(self._bound_points)


def _quadrature_tolerance(final: float, total: float) -> float:
    scale = max(abs(final), abs(total), math.ulp(0.0))
    return max(1024.0 * math.ulp(scale), 1.0e-24)


def _gauss_integral(
    function: Callable[[float], float],
    start: float,
    end: float,
    rule: tuple[tuple[float, float], ...],
) -> float:
    midpoint = 0.5 * (start + end)
    half_width = 0.5 * (end - start)
    return half_width * math.fsum(
        weight * function(math.fsum((midpoint, half_width * node)))
        for node, weight in rule
    )


def _scalar_times_interval(
    scalar: float,
    interval: tuple[float, float],
) -> tuple[float, float]:
    products = (scalar * interval[0], scalar * interval[1])
    return _outward_lower(min(products)), _outward_upper(max(products))


def _outward_lower(value: float) -> float:
    for _ in range(_ENCLOSURE_ROUNDOFF_ULPS):
        value = math.nextafter(value, -math.inf)
    return value


def _outward_upper(value: float) -> float:
    for _ in range(_ENCLOSURE_ROUNDOFF_ULPS):
        value = math.nextafter(value, math.inf)
    return value


def _interval_add(
    left: tuple[float, float],
    right: tuple[float, float],
) -> tuple[float, float]:
    return (
        _outward_lower(math.fsum((left[0], right[0]))),
        _outward_upper(math.fsum((left[1], right[1]))),
    )


def _interval_subtract(
    left: tuple[float, float],
    right: tuple[float, float],
) -> tuple[float, float]:
    return (
        _outward_lower(math.fsum((left[0], -right[1]))),
        _outward_upper(math.fsum((left[1], -right[0]))),
    )


def _interval_multiply(
    left: tuple[float, float],
    right: tuple[float, float],
) -> tuple[float, float]:
    products = tuple(a * b for a in left for b in right)
    return _outward_lower(min(products)), _outward_upper(max(products))


def _interval_scale(
    scalar: float,
    interval: tuple[float, float],
) -> tuple[float, float]:
    return _interval_multiply((scalar, scalar), interval)


def _native_virial_source_values(
    temperature_k: float,
    k_wh: float,
) -> tuple[float, float, float, float, float, float]:
    cross = bg.cross_second_virial(temperature_k, k_wh)
    values = (
        wa.second_virial(temperature_k),
        cross.B_wh,
        hx.second_virial(temperature_k),
        wa.second_virial_dT(temperature_k),
        cross.dB_wh_dT,
        hx.second_virial_dT(temperature_k),
    )
    envelopes = (
        _WATER_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        _CROSS_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        _HEXANE_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        _WATER_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
        _CROSS_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
        _HEXANE_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
    )
    if any(
        not lower <= value <= upper
        for value, (lower, upper) in zip(values, envelopes, strict=True)
    ):
        raise cp.CoupledPoreTopologyError(
            "analytic virial source left the native derivative envelope; fail closed"
        )
    return values


def _native_hexane_enthalpy_bounds(
    temperature_low_k: float,
    temperature_high_k: float,
    pressure_pa: float,
) -> tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float],
    tuple[float, float],
]:
    midpoint_k = math.fsum(
        (temperature_low_k, 0.5 * (temperature_high_k - temperature_low_k))
    )
    temperatures = (temperature_low_k, midpoint_k, temperature_high_k)
    liquid_states = tuple(hx.state_Tp(T, pressure_pa, "liquid") for T in temperatures)
    liquid_cp = tuple(state.cp for state in liquid_states)
    ideal_cp = tuple(hx.cp0(T) for T in temperatures)
    if any(
        not _HEXANE_LIQUID_CP_ENVELOPE_J_MOL_K[0]
        <= value
        <= _HEXANE_LIQUID_CP_ENVELOPE_J_MOL_K[1]
        for value in liquid_cp
    ) or any(
        not _HEXANE_IDEAL_CP_ENVELOPE_J_MOL_K[0]
        <= value
        <= _HEXANE_IDEAL_CP_ENVELOPE_J_MOL_K[1]
        for value in ideal_cp
    ):
        raise cp.CoupledPoreTopologyError(
            "analytic hexane heat-capacity source left the native derivative "
            "envelope; fail closed"
        )
    radius_k = max(
        midpoint_k - temperature_low_k,
        temperature_high_k - midpoint_k,
    )
    liquid_mid = liquid_states[1].h
    ideal_mid = hx.h_ideal(midpoint_k)
    liquid_radius = _outward_upper(
        _HEXANE_LIQUID_CP_ENVELOPE_J_MOL_K[1] * radius_k
    )
    ideal_radius = _outward_upper(
        _HEXANE_IDEAL_CP_ENVELOPE_J_MOL_K[1] * radius_k
    )
    liquid_bounds = (
        _outward_lower(liquid_mid - liquid_radius),
        _outward_upper(liquid_mid + liquid_radius),
    )
    ideal_bounds = (
        _outward_lower(ideal_mid - ideal_radius),
        _outward_upper(ideal_mid + ideal_radius),
    )
    return liquid_cp, ideal_cp, liquid_bounds, ideal_bounds


def _native_partial_residual_A_bounds(
    y_hexane_bounds: tuple[float, float],
    water_virial: tuple[float, float],
    cross_virial: tuple[float, float],
    hexane_virial: tuple[float, float],
) -> tuple[float, float]:
    y_hexane = y_hexane_bounds
    y_water = (
        _outward_lower(1.0 - y_hexane[1]),
        _outward_upper(1.0 - y_hexane[0]),
    )
    mixture = _interval_add(
        _interval_add(
            _interval_multiply(
                _interval_multiply(y_water, y_water),
                water_virial,
            ),
            _interval_scale(
                2.0,
                _interval_multiply(
                    _interval_multiply(y_water, y_hexane),
                    cross_virial,
                ),
            ),
        ),
        _interval_multiply(
            _interval_multiply(y_hexane, y_hexane),
            hexane_virial,
        ),
    )
    partial = _interval_scale(
        2.0,
        _interval_add(
            _interval_multiply(y_water, cross_virial),
            _interval_multiply(y_hexane, hexane_virial),
        ),
    )
    return _interval_subtract(partial, mixture)


def _native_hexane_partial_residual_enthalpy_bounds(
    temperature_low_k: float,
    temperature_high_k: float,
    y_hexane_low: float,
    y_hexane_high: float,
    pressure_pa: float,
    k_wh: float,
) -> tuple[float, float]:
    midpoint_k = math.fsum(
        (temperature_low_k, 0.5 * (temperature_high_k - temperature_low_k))
    )
    for temperature_k in (temperature_low_k, midpoint_k, temperature_high_k):
        _native_virial_source_values(temperature_k, k_wh)
    y_bounds = (y_hexane_low, y_hexane_high)
    A = _native_partial_residual_A_bounds(
        y_bounds,
        _WATER_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        _CROSS_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        _HEXANE_SECOND_VIRIAL_ENVELOPE_M3_MOL,
    )
    dA_dT = _native_partial_residual_A_bounds(
        y_bounds,
        _WATER_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
        _CROSS_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
        _HEXANE_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
    )
    A_minus_T_dA = _interval_subtract(
        A,
        _interval_multiply(
            (temperature_low_k, temperature_high_k),
            dA_dT,
        ),
    )
    return _interval_scale(pressure_pa, A_minus_T_dA)


def _native_hexane_temperature_derivative_bounds(
    temperature_low_k: float,
    temperature_high_k: float,
    liquid_enthalpy_bounds_j_mol: tuple[float, float],
    ideal_enthalpy_bounds_j_mol: tuple[float, float],
    partial_residual_enthalpy_bounds_j_mol: tuple[float, float],
) -> tuple[float, float]:
    numerator = _interval_subtract(
        _interval_subtract(
            liquid_enthalpy_bounds_j_mol,
            ideal_enthalpy_bounds_j_mol,
        ),
        partial_residual_enthalpy_bounds_j_mol,
    )
    denominator = (
        _outward_lower(hx.R * temperature_low_k * temperature_low_k),
        _outward_upper(hx.R * temperature_high_k * temperature_high_k),
    )
    reciprocal = (
        _outward_lower(1.0 / denominator[1]),
        _outward_upper(1.0 / denominator[0]),
    )
    return _interval_multiply(numerator, reciprocal)


def _native_closed_hexane_derivative_certificate(
    activities: _PathActivities,
) -> NativeClosedHexaneDerivativeCertificate:
    temperature_low_k = min(
        activities.left_temperature_k,
        activities.left_temperature_k + activities.delta_temperature_k,
    )
    temperature_high_k = max(
        activities.left_temperature_k,
        activities.left_temperature_k + activities.delta_temperature_k,
    )
    y_hexane_low = min(
        activities.left_y_hexane,
        activities.left_y_hexane + activities.delta_y_hexane,
    )
    y_hexane_high = max(
        activities.left_y_hexane,
        activities.left_y_hexane + activities.delta_y_hexane,
    )
    (
        liquid_cp,
        ideal_cp,
        liquid_enthalpy,
        ideal_enthalpy,
    ) = _native_hexane_enthalpy_bounds(
        temperature_low_k,
        temperature_high_k,
        activities.pressure_pa,
    )
    partial_residual = _native_hexane_partial_residual_enthalpy_bounds(
        temperature_low_k,
        temperature_high_k,
        y_hexane_low,
        y_hexane_high,
        activities.pressure_pa,
        activities.pore.k_wh,
    )
    temperature_derivative = _native_hexane_temperature_derivative_bounds(
        temperature_low_k,
        temperature_high_k,
        liquid_enthalpy,
        ideal_enthalpy,
        partial_residual,
    )
    composition_derivative = _composition_derivative_intervals(
        y_hexane_low,
        y_hexane_high,
    )[1]
    path_derivative = _path_derivative_bounds(
        activities.delta_y_hexane,
        activities.delta_temperature_k,
        composition_derivative,
        temperature_derivative,
    )
    for fraction in (0.0, 0.5, 1.0):
        point = activities.point(
            fraction,
            role="native_closed_hexane_derivative_source_containment",
        )
        if not temperature_derivative[0] <= (
            point.hexane_log_activity_temperature_derivative_k_inv
        ) <= temperature_derivative[1]:
            raise cp.CoupledPoreTopologyError(
                "analytic hexane temperature derivative escaped its native local "
                "enclosure; fail closed"
            )
    return NativeClosedHexaneDerivativeCertificate(
        left_temperature_k=activities.left_temperature_k,
        left_y_hexane=activities.left_y_hexane,
        right_temperature_k=(
            activities.left_temperature_k + activities.delta_temperature_k
        ),
        right_y_hexane=activities.left_y_hexane + activities.delta_y_hexane,
        pressure_pa=activities.pressure_pa,
        k_wh=activities.pore.k_wh,
        hexane_liquid_cp_samples_j_mol_k=liquid_cp,
        hexane_ideal_cp_samples_j_mol_k=ideal_cp,
        hexane_liquid_enthalpy_bounds_j_mol=liquid_enthalpy,
        hexane_ideal_enthalpy_bounds_j_mol=ideal_enthalpy,
        hexane_partial_residual_enthalpy_bounds_j_mol=partial_residual,
        hexane_log_activity_temperature_derivative_bounds_k_inv=(
            temperature_derivative
        ),
        hexane_log_activity_composition_derivative_bounds=(
            composition_derivative
        ),
        hexane_log_activity_path_derivative_bounds=path_derivative,
        source_envelope_id=_NATIVE_DERIVATIVE_SOURCE_ENVELOPE_ID,
        source_authority_ids=_NATIVE_DERIVATIVE_SOURCE_AUTHORITY_IDS,
        hexane_liquid_cp_envelope_j_mol_k=(
            _HEXANE_LIQUID_CP_ENVELOPE_J_MOL_K
        ),
        hexane_ideal_cp_envelope_j_mol_k=(
            _HEXANE_IDEAL_CP_ENVELOPE_J_MOL_K
        ),
        pure_and_cross_second_virial_envelopes_m3_mol=(
            _WATER_SECOND_VIRIAL_ENVELOPE_M3_MOL,
            _CROSS_SECOND_VIRIAL_ENVELOPE_M3_MOL,
            _HEXANE_SECOND_VIRIAL_ENVELOPE_M3_MOL,
        ),
        pure_and_cross_second_virial_dt_envelopes_m3_mol_k=(
            _WATER_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
            _CROSS_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
            _HEXANE_SECOND_VIRIAL_DT_ENVELOPE_M3_MOL_K,
        ),
    )


def _composition_derivative_intervals(
    lower_y_hexane: float,
    upper_y_hexane: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    water_raw = (
        -1.0 / (1.0 - upper_y_hexane)
        + (-_VIRIAL_FACTOR_MAX) * lower_y_hexane,
        -1.0 / (1.0 - lower_y_hexane)
        + (-_VIRIAL_FACTOR_MIN) * upper_y_hexane,
    )
    hexane_raw = (
        1.0 / upper_y_hexane
        + _VIRIAL_FACTOR_MIN * (1.0 - lower_y_hexane),
        1.0 / lower_y_hexane
        + _VIRIAL_FACTOR_MAX * (1.0 - upper_y_hexane),
    )
    water = (_outward_lower(water_raw[0]), _outward_upper(water_raw[1]))
    hexane = (_outward_lower(hexane_raw[0]), _outward_upper(hexane_raw[1]))
    if water[1] >= 0.0 or hexane[0] <= 0.0:
        raise cp.CoupledPoreTopologyError(
            "prospective activity envelope lost analytic composition monotonicity"
        )
    return water, hexane


def _path_derivative_bounds(
    delta_y_hexane: float,
    delta_temperature_k: float,
    composition_bounds: tuple[float, float],
    temperature_bounds: tuple[float, float],
) -> tuple[float, float]:
    composition = _scalar_times_interval(delta_y_hexane, composition_bounds)
    temperature = _scalar_times_interval(delta_temperature_k, temperature_bounds)
    return (
        _outward_lower(math.fsum((composition[0], temperature[0]))),
        _outward_upper(math.fsum((composition[1], temperature[1]))),
    )


def _activity_box(
    activities: _PathActivities,
    start: float,
    end: float,
    depth: int,
) -> _ActivityBox:
    left = activities.point(start, role="phase_enclosure_endpoint")
    right = activities.point(end, role="phase_enclosure_endpoint")
    temperature_low = min(left.temperature_k, right.temperature_k)
    temperature_high = max(left.temperature_k, right.temperature_k)
    y_low = min(left.y_hexane, right.y_hexane)
    y_high = max(left.y_hexane, right.y_hexane)
    water_lower = _outward_lower(
        activities.bound_point(temperature_high, y_high)[0]
    )
    water_upper = _outward_upper(
        activities.bound_point(temperature_low, y_low)[0]
    )
    hexane_lower = _outward_lower(
        activities.bound_point(temperature_high, y_low)[1]
    )
    hexane_upper = _outward_upper(
        activities.bound_point(temperature_low, y_high)[1]
    )
    water_composition, hexane_composition = _composition_derivative_intervals(
        y_low,
        y_high,
    )
    water_derivative = _path_derivative_bounds(
        activities.delta_y_hexane,
        activities.delta_temperature_k,
        water_composition,
        (
            _WATER_LOG_ACTIVITY_DT_MIN_K_INV,
            _WATER_LOG_ACTIVITY_DT_MAX_K_INV,
        ),
    )
    hexane_derivative = _path_derivative_bounds(
        activities.delta_y_hexane,
        activities.delta_temperature_k,
        hexane_composition,
        (
            _HEXANE_LOG_ACTIVITY_DT_MIN_K_INV,
            _HEXANE_LOG_ACTIVITY_DT_MAX_K_INV,
        ),
    )
    return _ActivityBox(
        start,
        end,
        depth,
        water_lower,
        water_upper,
        hexane_lower,
        hexane_upper,
        water_derivative,
        hexane_derivative,
    )


def _cap_residual(
    activities: _PathActivities,
    fraction: float,
    activity_at_cap: float,
) -> float:
    return math.fsum(
        (
            activities.point(fraction, role="cap_seam_isolation").water_activity,
            -activity_at_cap,
        )
    )


def _bisect_cap_seam(
    activities: _PathActivities,
    lower: float,
    upper: float,
    activity_at_cap: float,
    isolation_depth: int,
) -> CapSeamAudit:
    lower_value = _cap_residual(activities, lower, activity_at_cap)
    upper_value = _cap_residual(activities, upper, activity_at_cap)
    if lower_value == 0.0:
        return CapSeamAudit(lower, 0.0, 0.0, isolation_depth, 0)
    if upper_value == 0.0:
        return CapSeamAudit(upper, 0.0, 0.0, isolation_depth, 0)
    if lower_value * upper_value >= 0.0:
        raise cp.CoupledPoreTopologyError(
            "retained-cap seam locator requires a strict sign bracket"
        )
    best_fraction, best_value = min(
        ((lower, lower_value), (upper, upper_value)),
        key=lambda item: abs(item[1]),
    )
    iterations = 0
    for iteration in range(1, _CAP_ROOT_MAX_ITERATIONS + 1):
        iterations = iteration
        midpoint = lower + 0.5 * (upper - lower)
        if midpoint == lower or midpoint == upper:
            break
        midpoint_value = _cap_residual(activities, midpoint, activity_at_cap)
        if abs(midpoint_value) < abs(best_value):
            best_fraction, best_value = midpoint, midpoint_value
        if midpoint_value == 0.0:
            return CapSeamAudit(
                midpoint,
                0.0,
                upper - lower,
                isolation_depth,
                iterations,
            )
        if lower_value * midpoint_value < 0.0:
            upper = midpoint
            upper_value = midpoint_value
        else:
            lower = midpoint
            lower_value = midpoint_value
    tolerance = _CAP_CONTACT_ULPS * math.ulp(max(1.0, abs(activity_at_cap)))
    if abs(best_value) > tolerance:
        raise cp.CoupledPoreTopologyError(
            "retained-cap seam did not converge to a representable path point; "
            "do not smooth the active-set kink"
        )
    return CapSeamAudit(
        best_fraction,
        best_value,
        upper - lower,
        isolation_depth,
        iterations,
    )


def _deduplicate_cap_seams(seams: list[CapSeamAudit]) -> tuple[CapSeamAudit, ...]:
    unique: list[CapSeamAudit] = []
    for seam in sorted(seams, key=lambda item: item.fraction):
        if unique and abs(seam.fraction - unique[-1].fraction) <= 16.0 * math.ulp(
            max(1.0, abs(seam.fraction))
        ):
            if abs(seam.water_activity_minus_cap) < abs(
                unique[-1].water_activity_minus_cap
            ):
                unique[-1] = seam
            continue
        unique.append(seam)
    return tuple(unique)


def _locate_cap_seams(
    activities: _PathActivities,
    leaves: tuple[ActivityEnclosureLeafAudit, ...],
    activity_at_cap: float,
) -> tuple[tuple[CapSeamAudit, ...], int]:
    seams: list[CapSeamAudit] = []
    maximum_depth = 0
    def isolate(start: float, end: float, depth: int) -> None:
        nonlocal maximum_depth
        maximum_depth = max(maximum_depth, depth)
        box = _activity_box(activities, start, end, depth)
        if not (
            box.water_activity_lower_bound
            <= activity_at_cap
            <= box.water_activity_upper_bound
        ):
            return
        left_value = _cap_residual(activities, start, activity_at_cap)
        right_value = _cap_residual(activities, end, activity_at_cap)
        derivative_excludes_zero = (
            box.water_path_derivative_bounds[1] < -_DERIVATIVE_CERTIFICATE_MARGIN
            or box.water_path_derivative_bounds[0] > _DERIVATIVE_CERTIFICATE_MARGIN
        )
        if derivative_excludes_zero:
            if left_value == 0.0:
                seams.append(CapSeamAudit(start, 0.0, 0.0, depth, 0))
            if left_value * right_value < 0.0:
                seams.append(
                    _bisect_cap_seam(
                        activities,
                        start,
                        end,
                        activity_at_cap,
                        depth,
                    )
                )
            if right_value == 0.0:
                seams.append(CapSeamAudit(end, 0.0, 0.0, depth, 0))
            return
        width = end - start
        if depth >= _PHASE_ENCLOSURE_MAX_DEPTH or (
            width <= _PHASE_ENCLOSURE_MIN_FRACTION_WIDTH
        ):
            raise cp.CoupledPoreTopologyError(
                "retained-cap contact cannot be isolated by the bounded derivative "
                "enclosure; fail closed rather than smooth the seam"
            )
        midpoint = start + 0.5 * width
        activities.point(midpoint, role="cap_seam_subdivision")
        isolate(start, midpoint, depth + 1)
        isolate(midpoint, end, depth + 1)

    for leaf in leaves:
        isolate(leaf.start_fraction, leaf.end_fraction, leaf.depth)
    return _deduplicate_cap_seams(seams), maximum_depth


def _certify_path(
    activities: _PathActivities,
    derivative_authority: ClosedHexaneDerivativeAuthority,
) -> _PathCertificate:
    caller_left = activities.point(0.0, role="phase_enclosure_endpoint")
    activities.point(1.0, role="phase_enclosure_endpoint")
    closed = 1.0 - caller_left.hexane_activity <= _LEFT_HEXANE_CLOSED_MARGIN
    closed_derivative_upper: float | None = None
    native_derivative_certificate: NativeClosedHexaneDerivativeCertificate | None = None
    if closed:
        if (
            derivative_authority
            is ClosedHexaneDerivativeAuthority.NATIVE_LOCAL_ANALYTIC_ENVELOPE
        ):
            native_derivative_certificate = (
                _native_closed_hexane_derivative_certificate(activities)
            )
            closed_derivative_upper = (
                native_derivative_certificate.hexane_log_activity_path_derivative_bounds[
                    1
                ]
            )
        else:
            whole_box = _activity_box(activities, 0.0, 1.0, 0)
            internal_bounds = whole_box.hexane_path_derivative_bounds
            closed_derivative_upper = internal_bounds[1]
        if closed_derivative_upper >= -_DERIVATIVE_CERTIFICATE_MARGIN:
            raise cp.CoupledPoreTopologyError(
                "closed left hexane interface lacks a globally inward analytic "
                "log-activity derivative bound; fail closed"
            )

    leaves: list[ActivityEnclosureLeafAudit] = []

    def enclose(start: float, end: float, depth: int) -> None:
        box = _activity_box(activities, start, end, depth)
        water_enclosed = (
            box.water_activity_lower_bound > activities.activity_at_ref
            and box.water_activity_upper_bound < 1.0
        )
        hexane_enclosed = closed or box.hexane_activity_upper_bound < 1.0
        if water_enclosed and hexane_enclosed:
            leaves.append(
                ActivityEnclosureLeafAudit(
                    start_fraction=start,
                    end_fraction=end,
                    depth=depth,
                    water_activity_lower_bound=box.water_activity_lower_bound,
                    water_activity_upper_bound=box.water_activity_upper_bound,
                    hexane_activity_lower_bound=box.hexane_activity_lower_bound,
                    hexane_activity_upper_bound=box.hexane_activity_upper_bound,
                    hexane_upper_bound_replaced_by_closed_endpoint_derivative=closed,
                )
            )
            return
        width = end - start
        if depth >= _PHASE_ENCLOSURE_MAX_DEPTH or (
            width <= _PHASE_ENCLOSURE_MIN_FRACTION_WIDTH
        ):
            raise cp.CoupledPoreTopologyError(
                "actual path cannot be enclosed inside the caller-oriented phase "
                "domain at the frozen subdivision limit; fail closed"
            )
        midpoint = start + 0.5 * width
        activities.point(midpoint, role="phase_enclosure_subdivision")
        enclose(start, midpoint, depth + 1)
        enclose(midpoint, end, depth + 1)

    enclose(0.0, 1.0, 0)
    leaf_tuple = tuple(leaves)
    seams, maximum_cap_depth = _locate_cap_seams(
        activities,
        leaf_tuple,
        activities.activity_at_cap,
    )
    return _PathCertificate(
        enclosure_leaves=leaf_tuple,
        cap_seams=seams,
        maximum_cap_isolation_depth=maximum_cap_depth,
        left_hexane_endpoint_treated_as_closed=closed,
        closed_endpoint_hexane_log_activity_derivative_upper_bound=(
            closed_derivative_upper
        ),
        closed_hexane_derivative_authority=derivative_authority,
        native_closed_hexane_derivative_certificate=native_derivative_certificate,
    )


def _build_phase_audit(
    activities: _PathActivities,
    certificate: _PathCertificate,
) -> ActualPathPhaseAudit:
    nodes = activities.node_audits()
    return ActualPathPhaseAudit(
        evaluated_actual_path_nodes=nodes,
        enclosure_leaves=certificate.enclosure_leaves,
        cap_seams=certificate.cap_seams,
        minimum_y_hexane=min(node.y_hexane for node in nodes),
        minimum_y_water=min(1.0 - node.y_hexane for node in nodes),
        minimum_water_activity_minus_luikov_ref=min(
            node.water_activity - activities.activity_at_ref for node in nodes
        ),
        minimum_one_minus_water_activity=min(
            1.0 - node.water_activity for node in nodes
        ),
        minimum_one_minus_hexane_activity=min(
            1.0 - node.hexane_activity for node in nodes
        ),
        water_activity_at_ref=activities.activity_at_ref,
        water_activity_at_cap=activities.activity_at_cap,
        actual_path_activity_evaluation_count=len(nodes),
        off_path_bound_evaluation_count=activities.off_path_evaluation_count,
        maximum_enclosure_depth=max(leaf.depth for leaf in certificate.enclosure_leaves),
        maximum_cap_isolation_depth=certificate.maximum_cap_isolation_depth,
        left_hexane_endpoint_treated_as_closed=(
            certificate.left_hexane_endpoint_treated_as_closed
        ),
        closed_endpoint_hexane_log_activity_derivative_upper_bound=(
            certificate.closed_endpoint_hexane_log_activity_derivative_upper_bound
        ),
        closed_hexane_derivative_authority=(
            certificate.closed_hexane_derivative_authority
        ),
        native_closed_hexane_derivative_certificate=(
            certificate.native_closed_hexane_derivative_certificate
        ),
        virial_factor_envelope=(_VIRIAL_FACTOR_MIN, _VIRIAL_FACTOR_MAX),
        water_log_activity_temperature_derivative_envelope_k_inv=(
            _WATER_LOG_ACTIVITY_DT_MIN_K_INV,
            _WATER_LOG_ACTIVITY_DT_MAX_K_INV,
        ),
        hexane_log_activity_temperature_derivative_envelope_k_inv=(
            _HEXANE_LOG_ACTIVITY_DT_MIN_K_INV,
            _HEXANE_LOG_ACTIVITY_DT_MAX_K_INV,
        ),
        project_pressure_envelope_pa=(
            bg.PROJECT_PRESSURE_MIN_PA,
            bg.PROJECT_PRESSURE_MAX_PA,
        ),
        temperature_envelope_k=_TEMPERATURE_ENVELOPE_K,
        composition_envelope_y_hexane=_COMPOSITION_ENVELOPE_Y_HEXANE,
        k_wh_envelope=_K_WH_ENVELOPE,
        deterministic_guard_passed=True,
    )


def _gas_nonisothermal_force(
    activities: _PathActivities,
) -> tuple[float, float, float, GaussLegendreAudit]:
    delta_y = activities.delta_y_hexane
    left_y = activities.left_y_hexane
    ideal = math.fsum(
        (
            math.log1p(-delta_y / (1.0 - left_y)),
            -math.log1p(delta_y / left_y),
        )
    )

    def integrand(fraction: float) -> float:
        return activities.point(fraction, role="gas_quadrature").virial_factor

    comparison = -delta_y * _gauss_integral(
        integrand,
        0.0,
        1.0,
        _GAUSS_LEGENDRE_8,
    )
    final = -delta_y * _gauss_integral(
        integrand,
        0.0,
        1.0,
        _GAUSS_LEGENDRE_16,
    )
    total = math.fsum((ideal, final))
    change = abs(final - comparison)
    tolerance = _quadrature_tolerance(final, total)
    converged = change <= tolerance
    audit = GaussLegendreAudit(
        quantity="gas pressure-second-virial composition-force contribution",
        comparison_order=8,
        final_order=16,
        comparison_estimate=comparison,
        final_estimate=final,
        absolute_change=change,
        acceptance_tolerance=tolerance,
        converged=converged,
    )
    if not converged:
        raise ActualCompositionForceConvergenceError(
            "gas actual-path virial quadrature did not converge at fixed orders 8/16"
        )
    return total, ideal, final, audit


def _active_set_on_segment(
    activities: _PathActivities,
    start: float,
    end: float,
    activity_at_cap: float,
) -> cp.RetainedWaterActiveSet:
    probes = (
        start + 0.5 * (end - start),
        start + 0.25 * (end - start),
        start + 0.75 * (end - start),
    )
    tolerance = _CAP_CONTACT_ULPS * math.ulp(max(1.0, abs(activity_at_cap)))
    signs: list[int] = []
    for fraction in probes:
        residual = _cap_residual(activities, fraction, activity_at_cap)
        if abs(residual) <= tolerance:
            continue
        signs.append(1 if residual > 0.0 else -1)
    if not signs or any(sign != signs[0] for sign in signs[1:]):
        raise cp.CoupledPoreTopologyError(
            "retained-force path segment has an unresolved cap seam; do not smooth"
        )
    if signs[0] > 0:
        return cp.RetainedWaterActiveSet.RETAINED_CAP
    return cp.RetainedWaterActiveSet.LUIKOV_SMOOTH


def _smooth_retained_segment(
    activities: _PathActivities,
    start: float,
    end: float,
) -> RetainedPathSegmentAudit:
    delta_y = activities.delta_y_hexane
    start_point = activities.point(start, role="retained_segment_endpoint")
    end_point = activities.point(end, role="retained_segment_endpoint")
    segment_delta_y = math.fsum((end_point.y_hexane, -start_point.y_hexane))
    ideal = math.log1p(-segment_delta_y / (1.0 - start_point.y_hexane))
    if delta_y == 0.0:
        return RetainedPathSegmentAudit(
            start,
            end,
            cp.RetainedWaterActiveSet.LUIKOV_SMOOTH,
            0.0,
            0.0,
            0.0,
            None,
        )

    def integrand(fraction: float) -> float:
        point = activities.point(
            fraction,
            role="retained_quadrature",
        )
        return point.y_hexane * point.virial_factor

    comparison = -delta_y * _gauss_integral(
        integrand,
        start,
        end,
        _GAUSS_LEGENDRE_8,
    )
    final = -delta_y * _gauss_integral(
        integrand,
        start,
        end,
        _GAUSS_LEGENDRE_16,
    )
    total = math.fsum((ideal, final))
    change = abs(final - comparison)
    tolerance = _quadrature_tolerance(final, total)
    converged = change <= tolerance
    audit = GaussLegendreAudit(
        quantity="smooth retained-water pressure-second-virial force contribution",
        comparison_order=8,
        final_order=16,
        comparison_estimate=comparison,
        final_estimate=final,
        absolute_change=change,
        acceptance_tolerance=tolerance,
        converged=converged,
    )
    if not converged:
        raise ActualCompositionForceConvergenceError(
            "retained actual-path virial quadrature did not converge at fixed orders 8/16"
        )
    return RetainedPathSegmentAudit(
        start,
        end,
        cp.RetainedWaterActiveSet.LUIKOV_SMOOTH,
        total,
        ideal,
        final,
        audit,
    )


def _retained_nonisothermal_force(
    activities: _PathActivities,
    cap_seams: tuple[CapSeamAudit, ...],
) -> tuple[float, tuple[RetainedPathSegmentAudit, ...]]:
    boundaries = [0.0]
    boundaries.extend(
        seam.fraction
        for seam in cap_seams
        if 0.0 < seam.fraction < 1.0
    )
    boundaries.append(1.0)
    activity_at_cap = activities.activity_at_cap
    segments: list[RetainedPathSegmentAudit] = []
    for start, end in zip(boundaries, boundaries[1:]):
        if start == end:
            continue
        active_set = _active_set_on_segment(
            activities,
            start,
            end,
            activity_at_cap,
        )
        if active_set is cp.RetainedWaterActiveSet.RETAINED_CAP:
            segments.append(
                RetainedPathSegmentAudit(
                    start,
                    end,
                    active_set,
                    0.0,
                    0.0,
                    0.0,
                    None,
                )
            )
        else:
            segments.append(_smooth_retained_segment(activities, start, end))
    return (
        math.fsum(segment.force_difference_right_minus_left for segment in segments),
        tuple(segments),
    )


def _retained_isothermal_segments(
    activities: _PathActivities,
    cap_seams: tuple[CapSeamAudit, ...],
    exact_difference: float,
) -> tuple[RetainedPathSegmentAudit, ...]:
    boundaries = [0.0]
    boundaries.extend(
        seam.fraction
        for seam in cap_seams
        if 0.0 < seam.fraction < 1.0
    )
    boundaries.append(1.0)
    segments: list[RetainedPathSegmentAudit] = []
    for start, end in zip(boundaries, boundaries[1:]):
        if start == end:
            continue
        active_set = _active_set_on_segment(
            activities,
            start,
            end,
            activities.activity_at_cap,
        )
        if active_set is cp.RetainedWaterActiveSet.RETAINED_CAP:
            segments.append(
                RetainedPathSegmentAudit(start, end, active_set, 0.0, 0.0, 0.0, None)
            )
            continue
        start_point = activities.point(start, role="retained_segment_endpoint")
        end_point = activities.point(end, role="retained_segment_endpoint")
        segment_delta_y = math.fsum((end_point.y_hexane, -start_point.y_hexane))
        if segment_delta_y == 0.0:
            difference = 0.0
        else:
            local = cp.isothermal_potential_composition_secant_slopes(
                start_point.temperature_k,
                activities.pressure_pa,
                start_point.y_hexane,
                segment_delta_y,
                activities.pore,
            )
            difference = (
                local.retained_water_potential_secant_per_y_hexane
                * segment_delta_y
            )
        segments.append(
            RetainedPathSegmentAudit(
                start,
                end,
                active_set,
                difference,
                difference,
                0.0,
                None,
            )
        )
    smooth_indices = [
        index
        for index, segment in enumerate(segments)
        if segment.active_set is cp.RetainedWaterActiveSet.LUIKOV_SMOOTH
    ]
    if smooth_indices:
        last = smooth_indices[-1]
        other = math.fsum(
            segment.force_difference_right_minus_left
            for index, segment in enumerate(segments)
            if index != last
        )
        corrected = math.fsum((exact_difference, -other))
        segments[last] = replace(
            segments[last],
            force_difference_right_minus_left=corrected,
            ideal_mixture_contribution=corrected,
        )
    elif exact_difference != 0.0:
        raise cp.CoupledPoreTopologyError(
            "isothermal retained-cap path unexpectedly has a nonzero force"
        )
    return tuple(segments)


def _evaluate_ordered(
    left_temperature_k: float,
    left_y_hexane: float,
    right_temperature_k: float,
    right_y_hexane: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
    derivative_authority: ClosedHexaneDerivativeAuthority,
) -> _OrderedEvaluation:
    activities = _PathActivities(
        left_temperature_k,
        left_y_hexane,
        right_temperature_k,
        right_y_hexane,
        pressure_pa,
        pore,
    )
    certificate = _certify_path(activities, derivative_authority)
    delta_y = right_y_hexane - left_y_hexane
    if left_temperature_k == right_temperature_k:
        secants = cp.isothermal_potential_composition_secant_slopes(
            left_temperature_k,
            pressure_pa,
            left_y_hexane,
            delta_y,
            pore,
        )
        gas = secants.gas_exchange_potential_secant_per_y_hexane * delta_y
        retained = secants.retained_water_potential_secant_per_y_hexane * delta_y
        segments = _retained_isothermal_segments(
            activities,
            certificate.cap_seams,
            retained,
        )
        phase_audit = _build_phase_audit(activities, certificate)
        return _OrderedEvaluation(
            mode=CompositionForcePathMode.ISOTHERMAL_EXACT_SECANT,
            gas_difference=gas,
            retained_difference=retained,
            gas_ideal=None,
            gas_virial=None,
            gas_quadrature=None,
            retained_segments=segments,
            phase_audit=phase_audit,
        )
    if delta_y == 0.0:
        _, segments = _retained_nonisothermal_force(
            activities,
            certificate.cap_seams,
        )
        phase_audit = _build_phase_audit(activities, certificate)
        return _OrderedEvaluation(
            mode=CompositionForcePathMode.NONISOTHERMAL_ACTUAL_PATH,
            gas_difference=0.0,
            retained_difference=0.0,
            gas_ideal=0.0,
            gas_virial=0.0,
            gas_quadrature=None,
            retained_segments=segments,
            phase_audit=phase_audit,
        )
    gas, ideal, virial, gas_quadrature = _gas_nonisothermal_force(activities)
    retained, segments = _retained_nonisothermal_force(
        activities,
        certificate.cap_seams,
    )
    phase_audit = _build_phase_audit(activities, certificate)
    return _OrderedEvaluation(
        mode=CompositionForcePathMode.NONISOTHERMAL_ACTUAL_PATH,
        gas_difference=gas,
        retained_difference=retained,
        gas_ideal=ideal,
        gas_virial=virial,
        gas_quadrature=gas_quadrature,
        retained_segments=segments,
        phase_audit=phase_audit,
    )


def evaluate_actual_composition_force_path(
    left_temperature_k: float,
    left_y_hexane: float,
    right_temperature_k: float,
    right_y_hexane: float,
    pressure_pa: float,
    distance_m: float,
    pore: cp.CoupledPoreParams = cp.CoupledPoreParams(),
    *,
    closed_hexane_derivative_authority: ClosedHexaneDerivativeAuthority = (
        ClosedHexaneDerivativeAuthority.GLOBAL_PROSPECTIVE_ENVELOPE
    ),
) -> ActualCompositionForcePath:
    r"""Return the phase-guarded composition force along the actual face path.

    Pressure and material parameters are fixed along this prospective moving-
    interface operator.  The two temperatures and compositions are physical
    endpoints, not values at a fabricated common temperature.  Invalid paths
    raise :class:`~coupled_pore.CoupledPoreTopologyError` without repair.

    Every path fraction follows the caller's orientation: fraction zero is the
    moving interface and fraction one is the dry endpoint.  Symmetric fixed
    Gauss rules and analytic ideal terms preserve reversal antisymmetry without
    obscuring these endpoint roles behind a canonical reordering.
    """

    values = (
        left_temperature_k,
        left_y_hexane,
        right_temperature_k,
        right_y_hexane,
        pressure_pa,
        distance_m,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("actual composition-force path inputs must be finite")
    if pressure_pa <= 0.0 or distance_m <= 0.0:
        raise ValueError("actual composition-force pressure and distance must be positive")
    if not isinstance(
        closed_hexane_derivative_authority,
        ClosedHexaneDerivativeAuthority,
    ):
        raise TypeError("actual composition-force path needs a named derivative authority")
    if not (
        _TEMPERATURE_ENVELOPE_K[0]
        <= left_temperature_k
        <= _TEMPERATURE_ENVELOPE_K[1]
        and _TEMPERATURE_ENVELOPE_K[0]
        <= right_temperature_k
        <= _TEMPERATURE_ENVELOPE_K[1]
    ):
        raise cp.CoupledPoreTopologyError(
            "actual-path endpoints are outside the frozen temperature enclosure evidence"
        )
    if not 0.0 < left_y_hexane < 1.0 or not 0.0 < right_y_hexane < 1.0:
        raise cp.CoupledPoreTopologyError(
            "actual composition-force endpoints require strict interior composition"
        )
    left_state = cp.evaluate_equilibrium(
        left_temperature_k,
        pressure_pa,
        left_y_hexane,
        pore,
    )
    right_state = cp.evaluate_equilibrium(
        right_temperature_k,
        pressure_pa,
        right_y_hexane,
        pore,
    )
    activity_at_ref = sp.water_activity(pore.luikov.W_ref, pore.luikov)
    if not (
        activity_at_ref < left_state.water_activity < 1.0
        and 0.0 < left_state.hexane_activity <= 1.0
    ):
        raise cp.CoupledPoreTopologyError(
            "left moving-interface endpoint requires strict water/Luikov margins; "
            "only hexane saturation may be closed"
        )
    if not (
        activity_at_ref < right_state.water_activity < 1.0
        and 0.0 < right_state.hexane_activity < 1.0
    ):
        raise cp.CoupledPoreTopologyError(
            "right dry endpoint requires strict water, Luikov, and hexane margins"
        )
    ordered = _evaluate_ordered(
        left_temperature_k,
        left_y_hexane,
        right_temperature_k,
        right_y_hexane,
        pressure_pa,
        pore,
        closed_hexane_derivative_authority,
    )
    gas_difference = ordered.gas_difference
    retained_difference = ordered.retained_difference
    delta_y = math.fsum((right_y_hexane, -left_y_hexane))
    if delta_y == 0.0:
        gas_secant = None
        retained_secant = None
    else:
        gas_secant = gas_difference / delta_y
        retained_secant = retained_difference / delta_y
    return ActualCompositionForcePath(
        mode=ordered.mode,
        left_state=left_state,
        right_state=right_state,
        distance_m=distance_m,
        delta_temperature_k=math.fsum(
            (right_temperature_k, -left_temperature_k)
        ),
        delta_y_hexane=delta_y,
        gas_force_difference_right_minus_left=gas_difference,
        retained_force_difference_right_minus_left=retained_difference,
        gas_force_gradient_m_inv=gas_difference / distance_m,
        retained_force_gradient_m_inv=retained_difference / distance_m,
        gas_force_secant_per_y_hexane=gas_secant,
        retained_force_secant_per_y_hexane=retained_secant,
        gas_ideal_mixture_contribution=(
            ordered.gas_ideal
        ),
        gas_virial_contribution=(
            ordered.gas_virial
        ),
        gas_quadrature=ordered.gas_quadrature,
        retained_segments=ordered.retained_segments,
        phase_audit=ordered.phase_audit,
    )


__all__ = [
    "ActivityEnclosureLeafAudit",
    "ActualCompositionForceConvergenceError",
    "ActualCompositionForcePath",
    "ActualPathNodeAudit",
    "ActualPathPhaseAudit",
    "CapSeamAudit",
    "ClosedHexaneDerivativeAuthority",
    "CompositionForcePathMode",
    "GaussLegendreAudit",
    "NativeClosedHexaneDerivativeCertificate",
    "RetainedPathSegmentAudit",
    "evaluate_actual_composition_force_path",
]
