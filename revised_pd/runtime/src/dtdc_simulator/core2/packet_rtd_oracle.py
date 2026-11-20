"""Deterministic equal-weight packet RTD quadrature for PR-02-A1.

This module is an architecture-only manufactured oracle.  It realizes the
frozen tanks-in-series path as cumulative sums of distinct unit-exponential
dwell coordinates.  A power-of-two, unscrambled Sobol digital net supplies a
deterministic product-copula quadrature; each coordinate is Latinized with the
conditional mean of every equal-probability exponential stratum.

The resulting thresholds are reference/proposal coordinates only.  They do
not admit, rekey, discharge, split, or mutate particle state.  In particular,
this module does not authenticate accepted flow, select a production clock
integrator, implement the PR-04/06/07 exit transaction, or confer physical or
plant-predictive qualification.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import InitVar, dataclass
from functools import lru_cache
from typing import ClassVar

import scipy
from scipy.special import gammainc
from scipy.stats import qmc

from .joint_carrier_oracle import DryMatterDischargeBasis, RTDResidenceAuthority


MINIMUM_STAGE_COUNT = 1
MAXIMUM_STAGE_COUNT = 16
MINIMUM_PACKET_COUNT = 2
MAXIMUM_PACKET_COUNT = 65_536
SOBOL_BITS = 30
PACKET_RTD_CONSTRUCTION_ID = "GT-PS-2-PR02-A1/SOBOL-UNSCRAMBLED-BITS30/EXPONENTIAL-STRATUM-CMEAN/V1"
PACKET_RTD_IMPLEMENTATION_ID = f"scipy-{scipy.__version__}/scipy.stats.qmc.Sobol"

_PATH_ISSUE_TOKEN = object()
_ISSUE_TOKEN = object()


def _require_exact_integer(
    name: str,
    value: int,
    *,
    minimum: int,
    maximum: int,
) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an exact integer")
    if not minimum <= value <= maximum:
        raise ValueError(f"{name} must be in [{minimum}, {maximum}]")


def _is_negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _require_binary64(
    name: str,
    value: float,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> None:
    if type(value) is not float:
        raise TypeError(f"{name} must be an exact binary64 float")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if positive and value <= 0.0:
        raise ValueError(f"{name} must be strictly positive")
    if nonnegative and (value < 0.0 or _is_negative_zero(value)):
        raise ValueError(f"{name} must be non-negative with canonical positive zero")


def _same_binary64(left: float, right: float) -> bool:
    return struct.pack(">d", left) == struct.pack(">d", right)


def _require_nonblank(name: str, value: str) -> None:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact string")
    if not value.strip():
        raise ValueError(f"{name} must not be blank")


def _definition_digest(
    *,
    stage_count: int,
    packet_count: int,
    construction_id: str,
    implementation_id: str,
    paths: tuple[PacketRTDPath, ...],
) -> str:
    digest = hashlib.sha256()
    digest.update(b"GT-PS-2/PR02-A1/PACKET-RTD-QUADRATURE/V1\0")
    for payload in (
        stage_count.to_bytes(8, "big"),
        packet_count.to_bytes(8, "big"),
        construction_id.encode("utf-8"),
        implementation_id.encode("utf-8"),
    ):
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    for path in paths:
        digest.update(path.packet_index.to_bytes(8, "big"))
        digest.update(struct.pack(">d", path.packet_fraction))
        for value in path.stratum_indices:
            digest.update(value.to_bytes(8, "big"))
        for value in path.dwell_exposures:
            digest.update(struct.pack(">d", value))
        for value in path.transition_exposures:
            digest.update(struct.pack(">d", value))
    return "sha256:" + digest.hexdigest()


def _power_of_two_exponent(packet_count: int) -> int:
    _require_exact_integer(
        "packet_count",
        packet_count,
        minimum=MINIMUM_PACKET_COUNT,
        maximum=MAXIMUM_PACKET_COUNT,
    )
    if packet_count & (packet_count - 1):
        raise ValueError("packet_count must be a power of two")
    return packet_count.bit_length() - 1


def _exponential_integral_antiderivative(survival_probability: float) -> float:
    """Return an antiderivative of ``-log(v)`` with its exact zero limit."""

    if survival_probability == 0.0:
        return 0.0
    return survival_probability * (1.0 - math.log(survival_probability))


def exponential_stratum_conditional_means(packet_count: int) -> tuple[float, ...]:
    """Return binary64 conditional means for all unit-exponential strata.

    The last value absorbs only the binary64 summation residual.  This makes
    the finite quadrature's dwell sum exactly ``packet_count`` under
    :func:`math.fsum`; it is not a physical correction or a material clamp.
    """

    _power_of_two_exponent(packet_count)
    count_as_float = float(packet_count)
    values = []
    for stratum in range(packet_count):
        upper_survival = (packet_count - stratum) / count_as_float
        lower_survival = (packet_count - stratum - 1) / count_as_float
        conditional_mean = count_as_float * (
            _exponential_integral_antiderivative(upper_survival)
            - _exponential_integral_antiderivative(lower_survival)
        )
        _require_binary64("exponential stratum conditional mean", conditional_mean, positive=True)
        values.append(conditional_mean)
    values[-1] = count_as_float - math.fsum(values[:-1])
    _require_binary64(
        "closure-corrected exponential stratum conditional mean",
        values[-1],
        positive=True,
    )
    if not _same_binary64(math.fsum(values), count_as_float):
        raise ArithmeticError("binary64 exponential-stratum mean closure failed")
    return tuple(values)


def _prefix_sums(values: tuple[float, ...]) -> tuple[float, ...]:
    prefix = tuple(math.fsum(values[:stop]) for stop in range(1, len(values) + 1))
    for value in prefix:
        _require_binary64("cumulative transition exposure", value, positive=True)
    return prefix


@lru_cache(maxsize=64)
def _sobol_stratum_rows(
    *,
    stage_count: int,
    packet_count: int,
) -> tuple[tuple[int, ...], ...]:
    """Return the complete pinned Sobol rank matrix as immutable integers."""

    _require_exact_integer(
        "stage_count",
        stage_count,
        minimum=MINIMUM_STAGE_COUNT,
        maximum=MAXIMUM_STAGE_COUNT,
    )
    exponent = _power_of_two_exponent(packet_count)
    engine = qmc.Sobol(
        d=stage_count,
        scramble=False,
        bits=SOBOL_BITS,
        optimization=None,
    )
    points = engine.random_base2(exponent)
    if points.shape != (packet_count, stage_count):
        raise RuntimeError("Sobol engine returned an unexpected point-array shape")
    rows = tuple(
        tuple(int(float(coordinate) * packet_count) for coordinate in point) for point in points
    )
    if any(not 0 <= stratum < packet_count for row in rows for stratum in row):
        raise RuntimeError("Sobol coordinate did not map to a declared stratum")
    return rows


@dataclass(frozen=True, slots=True)
class PacketRTDPath:
    """One equally weighted manufactured packet path in exposure space."""

    packet_index: int
    packet_fraction: float
    stratum_indices: tuple[int, ...]
    dwell_exposures: tuple[float, ...]
    transition_exposures: tuple[float, ...]
    _issue_token: InitVar[object | None] = None

    def __post_init__(self, _issue_token: object | None) -> None:
        if _issue_token is not _PATH_ISSUE_TOKEN:
            raise TypeError("packet RTD paths are issued only by the oracle builder")
        _require_exact_integer(
            "packet path index",
            self.packet_index,
            minimum=0,
            maximum=MAXIMUM_PACKET_COUNT - 1,
        )
        _require_binary64("packet fraction", self.packet_fraction, positive=True)
        if type(self.stratum_indices) is not tuple:
            raise TypeError("stratum_indices must be an exact tuple")
        if type(self.dwell_exposures) is not tuple:
            raise TypeError("dwell_exposures must be an exact tuple")
        if type(self.transition_exposures) is not tuple:
            raise TypeError("transition_exposures must be an exact tuple")
        if not self.stratum_indices:
            raise ValueError("a packet RTD path needs at least one stage")
        if not (
            len(self.stratum_indices) == len(self.dwell_exposures) == len(self.transition_exposures)
        ):
            raise ValueError("packet RTD path coordinate lengths must agree")
        for stratum in self.stratum_indices:
            _require_exact_integer(
                "packet path stratum",
                stratum,
                minimum=0,
                maximum=MAXIMUM_PACKET_COUNT - 1,
            )
        for value in self.dwell_exposures:
            _require_binary64("packet dwell exposure", value, positive=True)
        for value in self.transition_exposures:
            _require_binary64("packet transition exposure", value, positive=True)
        expected = _prefix_sums(self.dwell_exposures)
        if any(
            not _same_binary64(observed, reference)
            for observed, reference in zip(self.transition_exposures, expected, strict=True)
        ):
            raise ValueError("transition exposures must be binary64 prefix sums of dwell exposures")

    def __copy__(self) -> None:
        raise TypeError("packet RTD path authority views cannot be copied")

    def __deepcopy__(self, memo: object) -> None:
        del memo
        raise TypeError("packet RTD path authority views cannot be deep-copied")

    def __reduce_ex__(self, protocol: int) -> None:
        del protocol
        raise TypeError("packet RTD path authority views cannot be serialized")

    def __replace__(self, /, **changes: object) -> None:
        del changes
        raise TypeError("packet RTD path authority views cannot be replaced")


def _validate_quadrature_paths(
    *,
    stage_count: int,
    packet_count: int,
    paths: tuple[PacketRTDPath, ...],
) -> None:
    if type(paths) is not tuple:
        raise TypeError("quadrature paths must be an exact tuple")
    if len(paths) != packet_count:
        raise ValueError("quadrature path count must equal packet_count")
    if any(type(path) is not PacketRTDPath for path in paths):
        raise TypeError("quadrature paths must be exact PacketRTDPath values")
    if tuple(path.packet_index for path in paths) != tuple(range(packet_count)):
        raise ValueError("quadrature paths must use canonical packet-index order")

    expected_fraction = 1.0 / packet_count
    expected_means = exponential_stratum_conditional_means(packet_count)
    expected_rows = _sobol_stratum_rows(
        stage_count=stage_count,
        packet_count=packet_count,
    )
    expected_strata = set(range(packet_count))
    for path, expected_row in zip(paths, expected_rows, strict=True):
        path.__post_init__(_PATH_ISSUE_TOKEN)
        if len(path.stratum_indices) != stage_count:
            raise ValueError("quadrature path stage count is inconsistent")
        if path.stratum_indices != expected_row:
            raise ValueError("quadrature stratum matrix differs from the pinned Sobol construction")
        if not _same_binary64(path.packet_fraction, expected_fraction):
            raise ValueError("quadrature packet fractions must be exactly equal")
        for stratum, observed_dwell in zip(
            path.stratum_indices,
            path.dwell_exposures,
            strict=True,
        ):
            if stratum >= packet_count:
                raise ValueError("quadrature stratum lies outside packet resolution")
            if not _same_binary64(observed_dwell, expected_means[stratum]):
                raise ValueError("quadrature dwell does not match its exponential stratum")
    for stage_index in range(stage_count):
        observed_strata = {path.stratum_indices[stage_index] for path in paths}
        if observed_strata != expected_strata:
            raise ValueError("every dwell coordinate must use every stratum exactly once")

    if not _same_binary64(math.fsum(path.packet_fraction for path in paths), 1.0):
        raise ValueError("quadrature packet fractions must close exactly to one")
    for stage_index in range(stage_count):
        dwell_mean = math.fsum(path.dwell_exposures[stage_index] for path in paths)
        if not _same_binary64(dwell_mean, float(packet_count)):
            raise ValueError("every finite-Q dwell coordinate must have exact unit mean")
    final_mean = math.fsum(path.transition_exposures[-1] for path in paths) / packet_count
    if abs(final_mean - stage_count) > math.ulp(float(stage_count)):
        raise ValueError("finite-Q final transition mean exceeds one binary64 ulp")


@dataclass(frozen=True, slots=True)
class DeterministicPacketRTDQuadrature:
    """Evaluator-issued immutable PR-02-A1 packet schedule."""

    stage_count: int
    packet_count: int
    construction_id: str
    implementation_id: str
    paths: tuple[PacketRTDPath, ...]
    definition_digest: str
    _issue_token: InitVar[object | None] = None

    architecture_only: ClassVar[bool] = True
    manufactured_quadrature_evidence: ClassVar[bool] = True
    owner_approved_packet_semantics_applied: ClassVar[bool] = True
    final_threshold_is_reference_only: ClassVar[bool] = True
    source_authenticated: ClassVar[bool] = False
    accepted_flow_instrumentation_authenticated: ClassVar[bool] = False
    production_flow_binding_implemented: ClassVar[bool] = False
    production_clock_integrator_selected: ClassVar[bool] = False
    production_rekey_transaction_implemented: ClassVar[bool] = False
    production_exit_selection_implemented: ClassVar[bool] = False
    unequal_packet_weight_quadrature_implemented: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self, _issue_token: object | None) -> None:
        if _issue_token is not _ISSUE_TOKEN:
            raise TypeError("packet RTD quadratures are issued only by the oracle builder")
        _require_exact_integer(
            "stage_count",
            self.stage_count,
            minimum=MINIMUM_STAGE_COUNT,
            maximum=MAXIMUM_STAGE_COUNT,
        )
        _power_of_two_exponent(self.packet_count)
        _require_nonblank("construction_id", self.construction_id)
        _require_nonblank("implementation_id", self.implementation_id)
        if self.construction_id != PACKET_RTD_CONSTRUCTION_ID:
            raise ValueError("quadrature construction identity is not the selected PR-02-A1 engine")
        if self.implementation_id != PACKET_RTD_IMPLEMENTATION_ID:
            raise ValueError("quadrature implementation identity differs from the runtime engine")
        _validate_quadrature_paths(
            stage_count=self.stage_count,
            packet_count=self.packet_count,
            paths=self.paths,
        )
        expected_digest = _definition_digest(
            stage_count=self.stage_count,
            packet_count=self.packet_count,
            construction_id=self.construction_id,
            implementation_id=self.implementation_id,
            paths=self.paths,
        )
        if self.definition_digest != expected_digest:
            raise ValueError("quadrature definition digest does not match its complete schedule")

    def __copy__(self) -> None:
        raise TypeError("packet RTD quadrature authority views cannot be copied")

    def __deepcopy__(self, memo: object) -> None:
        del memo
        raise TypeError("packet RTD quadrature authority views cannot be deep-copied")

    def __reduce_ex__(self, protocol: int) -> None:
        del protocol
        raise TypeError("packet RTD quadrature authority views cannot be serialized")

    def __replace__(self, /, **changes: object) -> None:
        del changes
        raise TypeError("packet RTD quadrature authority views cannot be replaced")


def build_deterministic_packet_rtd_quadrature(
    *,
    stage_count: int,
    packet_count: int,
) -> DeterministicPacketRTDQuadrature:
    """Build the declared-engine PR-02-A1 schedule without moving state."""

    _require_exact_integer(
        "stage_count",
        stage_count,
        minimum=MINIMUM_STAGE_COUNT,
        maximum=MAXIMUM_STAGE_COUNT,
    )
    _power_of_two_exponent(packet_count)
    stratum_rows = _sobol_stratum_rows(
        stage_count=stage_count,
        packet_count=packet_count,
    )
    conditional_means = exponential_stratum_conditional_means(packet_count)
    packet_fraction = 1.0 / packet_count
    paths: list[PacketRTDPath] = []
    for packet_index, strata in enumerate(stratum_rows):
        dwell_exposures = tuple(conditional_means[stratum] for stratum in strata)
        paths.append(
            PacketRTDPath(
                packet_index=packet_index,
                packet_fraction=packet_fraction,
                stratum_indices=strata,
                dwell_exposures=dwell_exposures,
                transition_exposures=_prefix_sums(dwell_exposures),
                _issue_token=_PATH_ISSUE_TOKEN,
            )
        )
    issued_paths = tuple(paths)
    definition_digest = _definition_digest(
        stage_count=stage_count,
        packet_count=packet_count,
        construction_id=PACKET_RTD_CONSTRUCTION_ID,
        implementation_id=PACKET_RTD_IMPLEMENTATION_ID,
        paths=issued_paths,
    )
    return DeterministicPacketRTDQuadrature(
        stage_count=stage_count,
        packet_count=packet_count,
        construction_id=PACKET_RTD_CONSTRUCTION_ID,
        implementation_id=PACKET_RTD_IMPLEMENTATION_ID,
        paths=issued_paths,
        definition_digest=definition_digest,
        _issue_token=_ISSUE_TOKEN,
    )


def validate_deterministic_packet_rtd_quadrature(
    quadrature: DeterministicPacketRTDQuadrature,
) -> None:
    """Revalidate an issued schedule and its complete definition digest."""

    if type(quadrature) is not DeterministicPacketRTDQuadrature:
        raise TypeError("quadrature must be an exact DeterministicPacketRTDQuadrature")
    if quadrature.construction_id != PACKET_RTD_CONSTRUCTION_ID:
        raise ValueError("quadrature construction identity is not the selected PR-02-A1 engine")
    if quadrature.implementation_id != PACKET_RTD_IMPLEMENTATION_ID:
        raise ValueError("quadrature implementation identity differs from the runtime engine")
    _validate_quadrature_paths(
        stage_count=quadrature.stage_count,
        packet_count=quadrature.packet_count,
        paths=quadrature.paths,
    )
    expected_digest = _definition_digest(
        stage_count=quadrature.stage_count,
        packet_count=quadrature.packet_count,
        construction_id=quadrature.construction_id,
        implementation_id=quadrature.implementation_id,
        paths=quadrature.paths,
    )
    if quadrature.definition_digest != expected_digest:
        raise ValueError("quadrature definition digest does not match its complete schedule")


def restore_packet_rtd_path(
    *,
    stage_count: int,
    packet_count: int,
    packet_index: int,
    packet_fraction: float,
    stratum_indices: tuple[int, ...],
    dwell_exposures: tuple[float, ...],
    transition_exposures: tuple[float, ...],
) -> PacketRTDPath:
    """Restore one path only by exact replay of the selected RTD oracle.

    Snapshot decoders must not manufacture the oracle-token-sealed path from wire
    fields.  This entry point rebuilds the complete selected quadrature, then
    returns its issued member only when every supplied integer and binary64 bit
    pattern is identical.  It therefore restores the current deterministic
    variant without broadening the admitted path inventory.
    """

    quadrature = build_deterministic_packet_rtd_quadrature(
        stage_count=stage_count,
        packet_count=packet_count,
    )
    if type(packet_index) is not int or not 0 <= packet_index < packet_count:
        raise ValueError("restored packet RTD index lies outside the selected quadrature")
    issued = quadrature.paths[packet_index]
    if type(packet_fraction) is not float or not _same_binary64(
        packet_fraction,
        issued.packet_fraction,
    ):
        raise ValueError("restored packet RTD fraction differs from the selected oracle")
    if type(stratum_indices) is not tuple or stratum_indices != issued.stratum_indices:
        raise ValueError("restored packet RTD strata differ from the selected oracle")
    if (
        type(dwell_exposures) is not tuple
        or len(dwell_exposures) != len(issued.dwell_exposures)
        or any(
            type(observed) is not float or not _same_binary64(observed, expected)
            for observed, expected in zip(
                dwell_exposures,
                issued.dwell_exposures,
                strict=True,
            )
        )
    ):
        raise ValueError("restored packet RTD dwells differ from the selected oracle")
    if (
        type(transition_exposures) is not tuple
        or len(transition_exposures) != len(issued.transition_exposures)
        or any(
            type(observed) is not float or not _same_binary64(observed, expected)
            for observed, expected in zip(
                transition_exposures,
                issued.transition_exposures,
                strict=True,
            )
        )
    ):
        raise ValueError("restored packet RTD transitions differ from the selected oracle")
    return issued


@dataclass(frozen=True, slots=True)
class PacketRTDMoments:
    """Realized and frozen-reference moments in one declared coordinate."""

    realized_mean: float
    realized_variance: float
    reference_mean: float
    reference_variance: float

    def __post_init__(self) -> None:
        for name, value in (
            ("realized_mean", self.realized_mean),
            ("realized_variance", self.realized_variance),
            ("reference_mean", self.reference_mean),
            ("reference_variance", self.reference_variance),
        ):
            _require_binary64(name, value, positive=True)


def dimensionless_exit_moments(
    quadrature: DeterministicPacketRTDQuadrature,
) -> PacketRTDMoments:
    """Return finite-Q exit moments in dimensionless operational exposure."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    exits = tuple(path.transition_exposures[-1] for path in quadrature.paths)
    mean = math.fsum(exits) / quadrature.packet_count
    variance = math.fsum((value - mean) ** 2 for value in exits) / quadrature.packet_count
    return PacketRTDMoments(
        realized_mean=mean,
        realized_variance=variance,
        reference_mean=float(quadrature.stage_count),
        reference_variance=float(quadrature.stage_count),
    )


def constant_residence_exit_moments(
    quadrature: DeterministicPacketRTDQuadrature,
    residence_authority: RTDResidenceAuthority,
) -> PacketRTDMoments:
    """Scale the dimensionless oracle to one declared constant PHY-009 snapshot."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    if type(residence_authority) is not RTDResidenceAuthority:
        raise TypeError("residence_authority must be an exact RTDResidenceAuthority")
    dimensionless = dimensionless_exit_moments(quadrature)
    mean_residence_time_s = residence_authority.mean_residence_time_s
    scale_s = mean_residence_time_s / quadrature.stage_count
    realized_variance_s2 = dimensionless.realized_variance * scale_s * scale_s
    reference_variance_s2 = mean_residence_time_s * mean_residence_time_s / quadrature.stage_count
    for value in (scale_s, realized_variance_s2, reference_variance_s2):
        _require_binary64("constant-residence derived coordinate", value, positive=True)
    return PacketRTDMoments(
        # Exact finite-Q dimensionless closure makes the mathematical scaled
        # mean equal to the declared PHY-009 mean.  Returning the authority's
        # coordinate avoids inventing a second rounded residence identity.
        realized_mean=mean_residence_time_s,
        realized_variance=realized_variance_s2,
        reference_mean=mean_residence_time_s,
        reference_variance=reference_variance_s2,
    )


@dataclass(frozen=True, slots=True)
class StageExitDistribution:
    """Stage occupancies plus exit probability at one exposure coordinate."""

    exposure: float
    stage_probabilities: tuple[float, ...]
    exit_probability: float
    probability_sum_residual: float

    def __post_init__(self) -> None:
        _require_binary64("distribution exposure", self.exposure, nonnegative=True)
        if type(self.stage_probabilities) is not tuple or not self.stage_probabilities:
            raise TypeError("stage_probabilities must be a nonempty exact tuple")
        for value in self.stage_probabilities:
            _require_binary64("stage probability", value, nonnegative=True)
            if value > 1.0:
                raise ValueError("stage probabilities must not exceed one")
        _require_binary64("exit probability", self.exit_probability, nonnegative=True)
        if self.exit_probability > 1.0:
            raise ValueError("exit probability must not exceed one")
        _require_binary64("probability sum residual", self.probability_sum_residual)


def empirical_stage_exit_distribution(
    quadrature: DeterministicPacketRTDQuadrature,
    exposure: float,
) -> StageExitDistribution:
    """Evaluate the finite packet measure without creating movement proposals."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    _require_binary64("exposure", exposure, nonnegative=True)
    counts = [0] * quadrature.stage_count
    exited = 0
    for path in quadrature.paths:
        crossed = sum(exposure >= threshold for threshold in path.transition_exposures)
        if crossed == quadrature.stage_count:
            exited += 1
        else:
            counts[crossed] += 1
    denominator = float(quadrature.packet_count)
    stages = tuple(count / denominator for count in counts)
    exit_probability = exited / denominator
    residual = math.fsum((*stages, exit_probability)) - 1.0
    return StageExitDistribution(
        exposure=exposure,
        stage_probabilities=stages,
        exit_probability=exit_probability,
        probability_sum_residual=residual,
    )


def erlang_stage_exit_distribution(
    *,
    stage_count: int,
    exposure: float,
) -> StageExitDistribution:
    """Return the analytic tanks-in-series stage/exit marginal reference."""

    _require_exact_integer(
        "stage_count",
        stage_count,
        minimum=MINIMUM_STAGE_COUNT,
        maximum=MAXIMUM_STAGE_COUNT,
    )
    _require_binary64("exposure", exposure, nonnegative=True)
    stage_probabilities: list[float] = []
    probability = math.exp(-exposure)
    for stage_index in range(stage_count):
        if stage_index > 0:
            probability *= exposure / stage_index
        stage_probabilities.append(probability)
    exit_probability = float(gammainc(stage_count, exposure))
    residual = math.fsum((*stage_probabilities, exit_probability)) - 1.0
    return StageExitDistribution(
        exposure=exposure,
        stage_probabilities=tuple(stage_probabilities),
        exit_probability=exit_probability,
        probability_sum_residual=residual,
    )


def _covariance_matrix(rows: tuple[tuple[float, ...], ...]) -> tuple[tuple[float, ...], ...]:
    if not rows:
        raise ValueError("covariance requires at least one row")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("covariance rows must have one common positive width")
    count = len(rows)
    means = tuple(math.fsum(row[column] for row in rows) / count for column in range(width))
    return tuple(
        tuple(
            math.fsum((row[left] - means[left]) * (row[right] - means[right]) for row in rows)
            / count
            for right in range(width)
        )
        for left in range(width)
    )


def dwell_covariance_matrix(
    quadrature: DeterministicPacketRTDQuadrature,
) -> tuple[tuple[float, ...], ...]:
    """Return the finite-Q dwell covariance used by the joint-path gate."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    return _covariance_matrix(tuple(path.dwell_exposures for path in quadrature.paths))


def transition_covariance_matrix(
    quadrature: DeterministicPacketRTDQuadrature,
) -> tuple[tuple[float, ...], ...]:
    """Return covariance of cumulative thresholds for ``min(j,k)`` checks."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    return _covariance_matrix(tuple(path.transition_exposures for path in quadrature.paths))


def accepted_residence_exposure_rate_s_inv(
    *,
    quadrature: DeterministicPacketRTDQuadrature,
    residence_authority: RTDResidenceAuthority,
) -> float:
    """Return declared ``N*q_accepted/M``; this does not authenticate the flow."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    if type(residence_authority) is not RTDResidenceAuthority:
        raise TypeError("residence_authority must be an exact RTDResidenceAuthority")
    rate = (
        quadrature.stage_count
        * residence_authority.accepted_dry_matter_discharge_kg_s
        / residence_authority.dry_matter_holdup_kg
    )
    _require_binary64("accepted residence exposure rate", rate, positive=True)
    return rate


@dataclass(frozen=True, slots=True)
class ManufacturedAcceptedExposureState:
    """Declared exposure inputs, including the exact zero-flow reference.

    This value does not prove that the named flow was accepted by a production
    host.  It lets the manufactured oracle test zero-flow freeze without
    weakening the positive-flow :class:`RTDResidenceAuthority` used by the
    existing moment oracle.
    """

    physical_tray_id: str
    exposure_state_id: str
    accepted_flow_history_id: str
    discharge_basis: DryMatterDischargeBasis
    dry_matter_holdup_kg: float
    accepted_dry_matter_discharge_kg_s: float

    architecture_only: ClassVar[bool] = True
    source_authenticated: ClassVar[bool] = False
    production_flow_binding_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("physical_tray_id", self.physical_tray_id)
        _require_nonblank("exposure_state_id", self.exposure_state_id)
        _require_nonblank("accepted_flow_history_id", self.accepted_flow_history_id)
        if type(self.discharge_basis) is not DryMatterDischargeBasis:
            raise TypeError("discharge_basis must be an exact DryMatterDischargeBasis")
        if self.discharge_basis is not DryMatterDischargeBasis.LOCALLY_ACCEPTED_DISCHARGE:
            raise ValueError("exposure state requires the locally accepted discharge basis")
        _require_binary64(
            "exposure-state dry-matter holdup",
            self.dry_matter_holdup_kg,
            positive=True,
        )
        _require_binary64(
            "exposure-state accepted dry-matter discharge",
            self.accepted_dry_matter_discharge_kg_s,
            nonnegative=True,
        )


def manufactured_accepted_exposure_rate_s_inv(
    *,
    quadrature: DeterministicPacketRTDQuadrature,
    exposure_state: ManufacturedAcceptedExposureState,
) -> float:
    """Return declared ``N*q/M`` with exact zero-flow freeze semantics."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    if type(exposure_state) is not ManufacturedAcceptedExposureState:
        raise TypeError("exposure_state must be an exact ManufacturedAcceptedExposureState")
    rate = (
        quadrature.stage_count
        * exposure_state.accepted_dry_matter_discharge_kg_s
        / exposure_state.dry_matter_holdup_kg
    )
    if exposure_state.accepted_dry_matter_discharge_kg_s == 0.0:
        if rate != 0.0:
            raise ArithmeticError("zero accepted discharge must produce exactly zero exposure")
        return 0.0
    _require_binary64("manufactured accepted residence exposure rate", rate, positive=True)
    return rate


@dataclass(frozen=True, slots=True)
class ManufacturedExposureInterval:
    """One declared piecewise-constant interval for reference integration."""

    start_time_s: float
    end_time_s: float
    exposure_state: ManufacturedAcceptedExposureState

    architecture_only: ClassVar[bool] = True
    production_clock_integrator_selected: ClassVar[bool] = False
    accepted_flow_instrumentation_authenticated: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_binary64("exposure interval start", self.start_time_s, nonnegative=True)
        _require_binary64("exposure interval end", self.end_time_s, nonnegative=True)
        if self.end_time_s <= self.start_time_s:
            raise ValueError("exposure interval end must be strictly after its start")
        if type(self.exposure_state) is not ManufacturedAcceptedExposureState:
            raise TypeError("exposure_state must be an exact ManufacturedAcceptedExposureState")


def integrate_manufactured_piecewise_constant_exposure(
    *,
    quadrature: DeterministicPacketRTDQuadrature,
    intervals: tuple[ManufacturedExposureInterval, ...],
) -> float:
    """Integrate the declared reference exposure over contiguous intervals."""

    validate_deterministic_packet_rtd_quadrature(quadrature)
    if type(intervals) is not tuple or not intervals:
        raise TypeError("intervals must be a nonempty exact tuple")
    if any(type(interval) is not ManufacturedExposureInterval for interval in intervals):
        raise TypeError("intervals must contain exact ManufacturedExposureInterval values")
    physical_tray_id = intervals[0].exposure_state.physical_tray_id
    contributions: list[float] = []
    for index, interval in enumerate(intervals):
        if interval.exposure_state.physical_tray_id != physical_tray_id:
            raise ValueError("piecewise exposure intervals must refer to one physical tray")
        if index > 0 and not _same_binary64(interval.start_time_s, intervals[index - 1].end_time_s):
            raise ValueError("piecewise exposure intervals must be exactly contiguous")
        duration_s = interval.end_time_s - interval.start_time_s
        rate = manufactured_accepted_exposure_rate_s_inv(
            quadrature=quadrature,
            exposure_state=interval.exposure_state,
        )
        contribution = duration_s * rate
        if rate == 0.0:
            if contribution != 0.0:
                raise ArithmeticError("zero exposure rate must contribute exactly zero")
        else:
            _require_binary64("piecewise exposure contribution", contribution, positive=True)
        contributions.append(contribution)
    exposure = math.fsum(contributions)
    _require_binary64("piecewise accepted residence exposure", exposure, nonnegative=True)
    return exposure


__all__ = [
    "MAXIMUM_PACKET_COUNT",
    "MAXIMUM_STAGE_COUNT",
    "MINIMUM_PACKET_COUNT",
    "MINIMUM_STAGE_COUNT",
    "PACKET_RTD_CONSTRUCTION_ID",
    "PACKET_RTD_IMPLEMENTATION_ID",
    "SOBOL_BITS",
    "DeterministicPacketRTDQuadrature",
    "ManufacturedAcceptedExposureState",
    "ManufacturedExposureInterval",
    "PacketRTDMoments",
    "PacketRTDPath",
    "StageExitDistribution",
    "accepted_residence_exposure_rate_s_inv",
    "build_deterministic_packet_rtd_quadrature",
    "constant_residence_exit_moments",
    "dimensionless_exit_moments",
    "dwell_covariance_matrix",
    "empirical_stage_exit_distribution",
    "erlang_stage_exit_distribution",
    "exponential_stratum_conditional_means",
    "integrate_manufactured_piecewise_constant_exposure",
    "manufactured_accepted_exposure_rate_s_inv",
    "restore_packet_rtd_path",
    "transition_covariance_matrix",
    "validate_deterministic_packet_rtd_quadrature",
]
