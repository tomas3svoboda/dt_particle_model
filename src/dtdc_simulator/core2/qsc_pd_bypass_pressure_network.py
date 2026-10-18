"""Small quasi-steady prescribed-pressure authority for the three PD bypasses.

The sealed binary gas path is the four-node serial chain
``PD3 -> PD2 -> PD1 -> DOME``.  Each PD edge applies an externally prescribed
nonnegative pressure drop, then adds local hexane vapour and wall-to-gas heat
in its downstream plenum.  Consequently each upper edge carries every source
admitted below it.

``ZERO_DROP_COLETTO_ORACLE`` is the named negligible-drop limit.
``PRESCRIBED_PROFILE`` is boundary data, not a hydraulic prediction.  This
module has no finite bypass coefficient, Ergun term, sparge-tray pressure
prior, flash model, gas storage, or calibration claim.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import ClassVar


PD_BYPASS_TRAY_ORDER = ("PD3", "PD2", "PD1")
PD_BYPASS_NODE_ORDER = ("PD3", "PD2", "PD1", "DOME")


class QSCPDBypassPressureNetworkError(ValueError):
    """A prescribed profile, gas rate, or local source was inadmissible."""


class QSCPDBypassPressureMode(enum.Enum):
    ZERO_DROP_COLETTO_ORACLE = "ZERO_DROP_COLETTO_ORACLE"
    PRESCRIBED_PROFILE = "PRESCRIBED_PROFILE"


ZERO_DROP_COLETTO_ORACLE = QSCPDBypassPressureMode.ZERO_DROP_COLETTO_ORACLE
PRESCRIBED_PROFILE = QSCPDBypassPressureMode.PRESCRIBED_PROFILE


def _finite(name: str, value: float, *, positive: bool = False, nonnegative: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise QSCPDBypassPressureNetworkError(f"{name} must be a finite binary64 value")
    if positive and value <= 0.0:
        raise QSCPDBypassPressureNetworkError(f"{name} must be strictly positive")
    if nonnegative and value < 0.0:
        raise QSCPDBypassPressureNetworkError(f"{name} must be nonnegative")


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDPrescribedPressureDrops:
    """Pressure losses in physical bottom-to-top order."""

    pd3_pa: float
    pd2_pa: float
    pd1_pa: float

    def __post_init__(self) -> None:
        for tray_id, value in zip(PD_BYPASS_TRAY_ORDER, self.values_pa, strict=True):
            _finite(f"{tray_id} prescribed pressure drop", value, nonnegative=True)

    @property
    def values_pa(self) -> tuple[float, float, float]:
        return self.pd3_pa, self.pd2_pa, self.pd1_pa


@dataclass(frozen=True, slots=True)
class QSCPDUnmatchedPressureDropDiagnostic:
    """Optional literature stress band with unmatched flow architecture."""

    lower_per_tray_pa: ClassVar[float] = 1_000.0
    upper_per_tray_pa: ClassVar[float] = 5_000.0
    applicability: ClassVar[str] = "UNMATCHED_FLOW_ARCHITECTURE_DIAGNOSTIC_ONLY"
    target_hardware_matched: ClassVar[bool] = False
    calibrated: ClassVar[bool] = False
    nominal: ClassVar[bool] = False


UNMATCHED_1_TO_5_KPA_PER_TRAY_DIAGNOSTIC = QSCPDUnmatchedPressureDropDiagnostic()


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDBinaryGasRate:
    """Binary component and common-datum enthalpy rates at one pressure node."""

    pressure_pa: float
    water_molar_flow_mol_s: float
    hexane_molar_flow_mol_s: float
    common_datum_enthalpy_rate_w: float
    energy_datum_id: str

    def __post_init__(self) -> None:
        _finite("gas pressure", self.pressure_pa, positive=True)
        _finite("water molar flow", self.water_molar_flow_mol_s, nonnegative=True)
        _finite("hexane molar flow", self.hexane_molar_flow_mol_s, nonnegative=True)
        _finite("common-datum enthalpy rate", self.common_datum_enthalpy_rate_w)
        if type(self.energy_datum_id) is not str or not self.energy_datum_id.strip():
            raise QSCPDBypassPressureNetworkError("energy datum must be a nonblank string")
        if self.molar_flow_mol_s <= 0.0:
            raise QSCPDBypassPressureNetworkError("upward binary gas flow must be positive")

    @property
    def molar_flow_mol_s(self) -> float:
        return math.fsum((self.water_molar_flow_mol_s, self.hexane_molar_flow_mol_s))

    @property
    def water_mole_fraction(self) -> float:
        return self.water_molar_flow_mol_s / self.molar_flow_mol_s

    @property
    def hexane_mole_fraction(self) -> float:
        return self.hexane_molar_flow_mol_s / self.molar_flow_mol_s


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDLocalGasSource:
    """Hexane and energy entering one downstream PD plenum."""

    tray_id: str
    generated_hexane_molar_flow_mol_s: float
    generated_hexane_common_enthalpy_rate_w: float
    wall_to_gas_heat_w: float
    energy_datum_id: str

    def __post_init__(self) -> None:
        if self.tray_id not in PD_BYPASS_TRAY_ORDER:
            raise QSCPDBypassPressureNetworkError("local source names a non-PD tray")
        _finite(
            "generated hexane molar flow",
            self.generated_hexane_molar_flow_mol_s,
            nonnegative=True,
        )
        _finite(
            "generated hexane common-enthalpy rate",
            self.generated_hexane_common_enthalpy_rate_w,
        )
        _finite("wall-to-gas heat", self.wall_to_gas_heat_w)
        if type(self.energy_datum_id) is not str or not self.energy_datum_id.strip():
            raise QSCPDBypassPressureNetworkError(
                "local-source energy datum must be a nonblank string"
            )
        if (
            self.generated_hexane_molar_flow_mol_s == 0.0
            and self.generated_hexane_common_enthalpy_rate_w != 0.0
        ):
            raise QSCPDBypassPressureNetworkError(
                "zero generated hexane cannot carry generated-vapour enthalpy"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDBypassPressureNetworkRequest:
    """One algebraic network evaluation at an already common time."""

    mode: QSCPDBypassPressureMode
    dome_pressure_pa: float
    pd3_inlet: QSCPDBinaryGasRate
    local_sources: tuple[QSCPDLocalGasSource, ...]
    prescribed_pressure_drops: QSCPDPrescribedPressureDrops | None = None
    diagnostic_metadata: QSCPDUnmatchedPressureDropDiagnostic | None = None

    pressure_loss_predicted: ClassVar[bool] = False
    finite_resistance_coefficient_populated: ClassVar[bool] = False
    ergun_pressure_loss_used: ClassVar[bool] = False
    sparge_pressure_drop_prior_used: ClassVar[bool] = False
    calibrated: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.mode) is not QSCPDBypassPressureMode:
            raise QSCPDBypassPressureNetworkError("pressure mode has a foreign type")
        _finite("dome pressure", self.dome_pressure_pa, positive=True)
        if type(self.pd3_inlet) is not QSCPDBinaryGasRate:
            raise QSCPDBypassPressureNetworkError("PD3 inlet has a foreign type")
        if (
            type(self.local_sources) is not tuple
            or tuple(source.tray_id for source in self.local_sources) != PD_BYPASS_TRAY_ORDER
        ):
            raise QSCPDBypassPressureNetworkError(
                "local sources must be the serial PD3, PD2, PD1 tuple"
            )
        if self.mode is ZERO_DROP_COLETTO_ORACLE:
            if self.prescribed_pressure_drops is not None:
                raise QSCPDBypassPressureNetworkError(
                    "zero-drop oracle does not accept a prescribed profile"
                )
        elif type(self.prescribed_pressure_drops) is not QSCPDPrescribedPressureDrops:
            raise QSCPDBypassPressureNetworkError(
                "prescribed-profile mode requires three explicit drops"
            )
        elif not any(value > 0.0 for value in self.prescribed_pressure_drops.values_pa):
            raise QSCPDBypassPressureNetworkError(
                "an all-zero profile must use the named zero-drop oracle"
            )
        if self.diagnostic_metadata is not None and type(self.diagnostic_metadata) is not (
            QSCPDUnmatchedPressureDropDiagnostic
        ):
            raise QSCPDBypassPressureNetworkError("diagnostic metadata has a foreign type")


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDBypassEdgeResult:
    """One PD request boundary and its downstream-plenum balance."""

    tray_id: str
    inlet_node: str
    outlet_node: str
    bypass_inlet_pressure_pa: float
    prescribed_bypass_pressure_drop_pa: float
    plenum_pressure_pa: float
    inlet: QSCPDBinaryGasRate
    local_source: QSCPDLocalGasSource
    outlet: QSCPDBinaryGasRate

    def __post_init__(self) -> None:
        if self.tray_id not in PD_BYPASS_TRAY_ORDER:
            raise QSCPDBypassPressureNetworkError("edge names a non-PD tray")
        index = PD_BYPASS_TRAY_ORDER.index(self.tray_id)
        if (self.inlet_node, self.outlet_node) != PD_BYPASS_NODE_ORDER[index : index + 2]:
            raise QSCPDBypassPressureNetworkError("edge violates the serial PD node topology")
        if (
            type(self.inlet) is not QSCPDBinaryGasRate
            or type(self.outlet) is not QSCPDBinaryGasRate
            or type(self.local_source) is not QSCPDLocalGasSource
            or self.local_source.tray_id != self.tray_id
        ):
            raise QSCPDBypassPressureNetworkError("edge gas/source types or tray identity differ")
        if (
            self.inlet.pressure_pa != self.bypass_inlet_pressure_pa
            or self.outlet.pressure_pa != self.plenum_pressure_pa
        ):
            raise QSCPDBypassPressureNetworkError("edge stream pressure differs from its node")
        if self.inlet.energy_datum_id != self.outlet.energy_datum_id:
            raise QSCPDBypassPressureNetworkError("edge changes the common-energy datum")
        if not self.balances_exactly:
            raise QSCPDBypassPressureNetworkError("edge pressure/component/enthalpy balance failed")

    @property
    def pressure_residual_pa(self) -> float:
        return math.fsum(
            (
                self.bypass_inlet_pressure_pa,
                -self.prescribed_bypass_pressure_drop_pa,
                -self.plenum_pressure_pa,
            )
        )

    @property
    def balances_exactly(self) -> bool:
        return (
            self.pressure_residual_pa == 0.0
            and self.outlet.water_molar_flow_mol_s == self.inlet.water_molar_flow_mol_s
            and self.outlet.hexane_molar_flow_mol_s
            == math.fsum(
                (
                    self.inlet.hexane_molar_flow_mol_s,
                    self.local_source.generated_hexane_molar_flow_mol_s,
                )
            )
            and self.outlet.common_datum_enthalpy_rate_w
            == math.fsum(
                (
                    self.inlet.common_datum_enthalpy_rate_w,
                    self.local_source.generated_hexane_common_enthalpy_rate_w,
                    self.local_source.wall_to_gas_heat_w,
                )
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDBypassPressureNetworkResult:
    mode: QSCPDBypassPressureMode
    node_pressures_pa: tuple[float, float, float, float]
    edges: tuple[QSCPDBypassEdgeResult, ...]
    dome_outlet: QSCPDBinaryGasRate
    diagnostic_metadata: QSCPDUnmatchedPressureDropDiagnostic | None

    pressure_loss_predicted: ClassVar[bool] = False
    finite_resistance_coefficient_populated: ClassVar[bool] = False
    calibrated: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.mode) is not QSCPDBypassPressureMode:
            raise QSCPDBypassPressureNetworkError("result pressure mode has a foreign type")
        if type(self.node_pressures_pa) is not tuple or len(self.node_pressures_pa) != 4:
            raise QSCPDBypassPressureNetworkError("result must contain exactly four pressure nodes")
        for node_id, pressure_pa in zip(PD_BYPASS_NODE_ORDER, self.node_pressures_pa, strict=True):
            _finite(f"{node_id} result pressure", pressure_pa, positive=True)
        if not (
            self.node_pressures_pa[0]
            >= self.node_pressures_pa[1]
            >= self.node_pressures_pa[2]
            >= self.node_pressures_pa[3]
        ):
            raise QSCPDBypassPressureNetworkError("result pressure nodes reverse direction")
        if (
            type(self.edges) is not tuple
            or len(self.edges) != 3
            or any(type(edge) is not QSCPDBypassEdgeResult for edge in self.edges)
            or tuple(edge.tray_id for edge in self.edges) != PD_BYPASS_TRAY_ORDER
        ):
            raise QSCPDBypassPressureNetworkError("result must contain the three serial PD edges")
        for index, edge in enumerate(self.edges):
            if (
                edge.bypass_inlet_pressure_pa != self.node_pressures_pa[index]
                or edge.plenum_pressure_pa != self.node_pressures_pa[index + 1]
            ):
                raise QSCPDBypassPressureNetworkError("result edge and node pressures differ")
        if (
            self.edges[0].outlet != self.edges[1].inlet
            or self.edges[1].outlet != self.edges[2].inlet
        ):
            raise QSCPDBypassPressureNetworkError("adjacent PD edge gas states do not join")
        if (
            type(self.dome_outlet) is not QSCPDBinaryGasRate
            or self.dome_outlet != self.edges[2].outlet
        ):
            raise QSCPDBypassPressureNetworkError("PD1 outlet and dome outlet differ")
        drops = tuple(edge.prescribed_bypass_pressure_drop_pa for edge in self.edges)
        if self.mode is ZERO_DROP_COLETTO_ORACLE and any(value != 0.0 for value in drops):
            raise QSCPDBypassPressureNetworkError("zero-drop result contains a finite drop")
        if self.mode is PRESCRIBED_PROFILE and not any(value > 0.0 for value in drops):
            raise QSCPDBypassPressureNetworkError("prescribed result aliases the zero-drop oracle")
        if self.diagnostic_metadata is not None and type(self.diagnostic_metadata) is not (
            QSCPDUnmatchedPressureDropDiagnostic
        ):
            raise QSCPDBypassPressureNetworkError("result diagnostic metadata has a foreign type")


def evaluate_qsc_pd_bypass_pressure_network(
    request: QSCPDBypassPressureNetworkRequest,
) -> QSCPDBypassPressureNetworkResult:
    """Propagate pressure, binary components, and common enthalpy without storage."""

    if type(request) is not QSCPDBypassPressureNetworkRequest:
        raise QSCPDBypassPressureNetworkError("network request has a foreign type")
    drops = (
        QSCPDPrescribedPressureDrops(pd3_pa=0.0, pd2_pa=0.0, pd1_pa=0.0)
        if request.mode is ZERO_DROP_COLETTO_ORACLE
        else request.prescribed_pressure_drops
    )
    if drops is None:  # narrowed by the request validation
        raise QSCPDBypassPressureNetworkError("prescribed pressure drops are missing")

    pd1_pressure = math.fsum((request.dome_pressure_pa, drops.pd1_pa))
    pd2_pressure = math.fsum((pd1_pressure, drops.pd2_pa))
    pd3_pressure = math.fsum((pd2_pressure, drops.pd3_pa))
    node_pressures = (pd3_pressure, pd2_pressure, pd1_pressure, request.dome_pressure_pa)
    if not node_pressures[0] >= node_pressures[1] >= node_pressures[2] >= node_pressures[3]:
        raise QSCPDBypassPressureNetworkError("pressure profile reverses the upward gradient")
    if request.pd3_inlet.pressure_pa != node_pressures[0]:
        raise QSCPDBypassPressureNetworkError("PD3 inlet pressure differs from the profile")

    current = request.pd3_inlet
    edges: list[QSCPDBypassEdgeResult] = []
    for index, (tray_id, source, pressure_drop_pa) in enumerate(
        zip(PD_BYPASS_TRAY_ORDER, request.local_sources, drops.values_pa, strict=True)
    ):
        if source.energy_datum_id != current.energy_datum_id:
            raise QSCPDBypassPressureNetworkError("PD source and gas rate use different datums")
        outlet = QSCPDBinaryGasRate(
            pressure_pa=node_pressures[index + 1],
            water_molar_flow_mol_s=current.water_molar_flow_mol_s,
            hexane_molar_flow_mol_s=math.fsum(
                (
                    current.hexane_molar_flow_mol_s,
                    source.generated_hexane_molar_flow_mol_s,
                )
            ),
            common_datum_enthalpy_rate_w=math.fsum(
                (
                    current.common_datum_enthalpy_rate_w,
                    source.generated_hexane_common_enthalpy_rate_w,
                    source.wall_to_gas_heat_w,
                )
            ),
            energy_datum_id=current.energy_datum_id,
        )
        edge = QSCPDBypassEdgeResult(
            tray_id=tray_id,
            inlet_node=PD_BYPASS_NODE_ORDER[index],
            outlet_node=PD_BYPASS_NODE_ORDER[index + 1],
            bypass_inlet_pressure_pa=node_pressures[index],
            prescribed_bypass_pressure_drop_pa=pressure_drop_pa,
            plenum_pressure_pa=node_pressures[index + 1],
            inlet=current,
            local_source=source,
            outlet=outlet,
        )
        edges.append(edge)
        current = outlet

    return QSCPDBypassPressureNetworkResult(
        mode=request.mode,
        node_pressures_pa=node_pressures,
        edges=tuple(edges),
        dome_outlet=current,
        diagnostic_metadata=request.diagnostic_metadata,
    )


__all__ = [
    "PD_BYPASS_NODE_ORDER",
    "PD_BYPASS_TRAY_ORDER",
    "PRESCRIBED_PROFILE",
    "QSCPDBinaryGasRate",
    "QSCPDBypassEdgeResult",
    "QSCPDBypassPressureMode",
    "QSCPDBypassPressureNetworkError",
    "QSCPDBypassPressureNetworkRequest",
    "QSCPDBypassPressureNetworkResult",
    "QSCPDLocalGasSource",
    "QSCPDPrescribedPressureDrops",
    "QSCPDUnmatchedPressureDropDiagnostic",
    "UNMATCHED_1_TO_5_KPA_PER_TRAY_DIAGNOSTIC",
    "ZERO_DROP_COLETTO_ORACLE",
    "evaluate_qsc_pd_bypass_pressure_network",
]
