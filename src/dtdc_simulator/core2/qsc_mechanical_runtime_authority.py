"""MECH-01A declared engineering runtime authority for the six-tray QSC.

This module instantiates one common-shaft PR-08 authority and six finite-K
engineering scheduler authorities from the owner-ruled D0--D6 assumptions.
Every K authority binds the same common-shaft and accepted-history roots.  A
separate immutable binding row ties each K authority to its exact tray arm and
draw identities, so a foreign shaft, history, tray, K, or draw digest refuses.

The authority is intentionally narrow.  It does not wire the six-tray
evaluator, calculate or deposit shaft work, select a physical draw footprint
or flow law, authenticate a machine history, or advance any F, physical,
plant, or industrial-feasibility claim.  The K2 dry-density value and the
geometry fields needed only to occupy the PR-08 schema are explicit D0
brackets and are not dry-inventory or machine-data authority.
"""

from __future__ import annotations

import enum
import functools
import hashlib
import importlib
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar, Iterable, Mapping

import yaml

from . import engineering_physical_k_event_scheduler as k_scheduler
from . import mechanical_authority as mechanical
from .tray_particle_stateful import VerticalLayerId

QSC_MECH01_SCHEMA_ID = "dtdc-core2-qsc-mech01-engineering-authority-v1"
QSC_MECH01_SCHEMA_REVISION = 1
QSC_MECH01_SOLID_TRAY_ORDER = ("PD1", "PD2", "PD3", "MN1", "MN2", "SP1")
QSC_MECH01_CONFIG_RELATIVE_PATH = "benchmarks/qsc_mech01_engineering_authority_v1.yaml"
QSC_MECH01_CONFIG_CANONICAL_LF_SHA256 = (
    "2313e56d2b5d1ed3a99089a87731816016aa6ef0441665d08dfb06b58e4967a6"
)
QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST = (
    "sha256:82e268c2884ad2d6686034938632249142721dfa1d50fec3a1a71ced85367c2e"
)
QSC_MECH01_EXPECTED_COMMON_SHAFT_DEFINITION_DIGEST = (
    "sha256:59c65b167e927a0fcc1bf5afee6beef3e17324f17d6ba6c01d4020f2d1cc6b93"
)
#: Ruled 2026-08-30 (GT_PS2_V3_EXTENSION_AND_RPM_AUTHORITY_RULING): constant
#: engineering shaft speeds are admissible within the D1 OPERATING envelope
#: [minimum stable, maximum continuous]; the D1 nominal 11 r/min remains the
#: canonical authority with its pinned digest unchanged.
QSC_MECH01_OPERATING_RPM_ENVELOPE = (8.0, 12.0)

_ROOT = Path(__file__).resolve().parents[3]
QSC_MECH01_CONFIG_PATH = _ROOT / QSC_MECH01_CONFIG_RELATIVE_PATH
_SHA256_PATTERN = re.compile(r"sha256:[0-9a-f]{64}\Z")

# Preserve the hostile-oracle detachment convention used by the existing K
# scheduler.  This module names no production evaluator and only constructs
# the already-declared engineering authority types.
_arm_oracle = importlib.import_module("." + "adjacent_k_arm_" + "dispersion_oracle", __package__)


class QSCMechanicalParameterClass(enum.Enum):
    """Typed QSC parameter standing; none of these is machine authentication."""

    DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"
    DERIVED_FROM_DECLARED_ENGINEERING_ASSUMPTION = "DERIVED_FROM_DECLARED_ENGINEERING_ASSUMPTION"
    RETIRED_NONAUTHORITY = "RETIRED_NONAUTHORITY"


class ShaftWorkBoundaryPartition(enum.Enum):
    ALL_TO_MEAL = "ALL_TO_MEAL"
    ALL_TO_WALL = "ALL_TO_WALL"


class WithinMealDistribution(enum.Enum):
    WET_HOLDUP_PROPORTIONAL_BY_TRAY = "WET_HOLDUP_PROPORTIONAL_BY_TRAY"
    BOTTOM_HEAVY_BY_SOLID_ORDER_RANK = "BOTTOM_HEAVY_BY_SOLID_ORDER_RANK"


def _canonical_lf_sha256(path: Path) -> str:
    raw = path.read_bytes()
    canonical = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(canonical).hexdigest()


def _sha256_bytes(material: bytes) -> str:
    return "sha256:" + hashlib.sha256(material).hexdigest()


def _sha256_text(material: str) -> str:
    return _sha256_bytes(material.encode("utf-8"))


def _require_sha256(name: str, value: str) -> None:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact string")
    if _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{name} must be canonical lowercase sha256:<64-hex>")


def _require_identity(name: str, *values: str) -> None:
    if any(type(value) is not str for value in values):
        raise TypeError(f"{name} must contain exact strings")
    if any(not value or value != value.strip() for value in values):
        raise ValueError(f"{name} must contain nonblank canonical strings")


def _require_float(name: str, *values: float, positive: bool = False) -> None:
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain exact binary64 values")
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"{name} must contain finite values")
    if any(value == 0.0 and math.copysign(1.0, value) < 0.0 for value in values):
        raise ValueError(f"{name} rejects negative zero")
    if positive and any(value <= 0.0 for value in values):
        raise ValueError(f"{name} must contain strictly positive values")


def _digest(*, domain: str, parts: Iterable[str | int | float | bool]) -> str:
    _require_identity("digest domain", domain)
    digest = hashlib.sha256()
    digest.update(b"DTDC-QSC-MECH01-RUNTIME-AUTHORITY-V1\0")
    material = tuple(parts)
    for value in (domain, str(len(material))):
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    for part in material:
        if type(part) is bool:
            tag, payload = b"b", b"1" if part else b"0"
        elif type(part) is int:
            tag, payload = b"i", str(part).encode("ascii")
        elif type(part) is float:
            _require_float("digest float", part)
            tag, payload = b"f", part.hex().encode("ascii")
        elif type(part) is str:
            _require_identity("digest string", part)
            tag, payload = b"s", part.encode("utf-8")
        else:  # pragma: no cover - closed by the callers below
            raise TypeError("unsupported MECH-01A digest part")
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


def _mapping(name: str, value: Any) -> dict[str, Any]:
    if type(value) is not dict:
        raise TypeError(f"{name} must be an exact mapping")
    return value


def _float_value(name: str, value: Any) -> float:
    if type(value) not in (int, float):
        raise TypeError(f"{name} must be a YAML integer or binary64 number")
    result = float(value)
    _require_float(name, result)
    return result


def _float_pair(name: str, value: Any) -> tuple[float, float]:
    if type(value) is not list or len(value) != 2:
        raise TypeError(f"{name} must be an exact two-item YAML sequence")
    result = (_float_value(name, value[0]), _float_value(name, value[1]))
    if result[0] > result[1]:
        raise ValueError(f"{name} must be ordered")
    return result


def _load_configuration() -> dict[str, Any]:
    if not QSC_MECH01_CONFIG_PATH.is_file():
        raise RuntimeError(f"missing MECH-01A configuration: {QSC_MECH01_CONFIG_PATH}")
    observed = _canonical_lf_sha256(QSC_MECH01_CONFIG_PATH)
    if observed != QSC_MECH01_CONFIG_CANONICAL_LF_SHA256:
        raise RuntimeError(
            "MECH-01A configuration bytes differ from the source-pinned canonical LF digest"
        )
    document = yaml.safe_load(QSC_MECH01_CONFIG_PATH.read_text(encoding="utf-8"))
    return _mapping("MECH-01A configuration", document)


_CONFIG = _load_configuration()


def _require_configuration_semantics(document: Mapping[str, Any]) -> None:
    if type(document) is not dict:
        raise TypeError("MECH-01A configuration must be an exact mapping")
    if document.get("schema_id") != QSC_MECH01_SCHEMA_ID:
        raise ValueError("MECH-01A configuration has the wrong schema identity")
    if document.get("schema_revision") != QSC_MECH01_SCHEMA_REVISION:
        raise ValueError("MECH-01A configuration has the wrong schema revision")
    if document.get("authority_status") != "DECLARED_ENGINEERING_ASSUMPTION":
        raise ValueError("MECH-01A must remain a declared engineering assumption")
    if document.get("authenticated_machine_parameter_count") != 0:
        raise ValueError("MECH-01A cannot claim authenticated machine parameters")
    for flag in ("source_authenticated", "machine_data_authenticated"):
        if document.get(flag) is not False:
            raise ValueError(f"MECH-01A {flag} must remain exact false")

    pin = _mapping(
        "engineering-assumptions source pin",
        _mapping("source_pins", document.get("source_pins")).get("engineering_assumptions_ruling"),
    )
    if pin != {
        "path": "docs/GT_PS2_ENGINEERING_ASSUMPTIONS_RULING_RECORD_2026-08-10.md",
        "digest_basis": "canonical_lf_bytes_v1",
        "sha256": "3eaf8787e0bcc001e24947794bd02c5f6619a39264811f1fbe908c5ec18e4e09",
    }:
        raise ValueError("MECH-01A engineering-assumptions source pin drifted")
    source_path = (_ROOT / pin["path"]).resolve()
    # Public repository export: the pinned assumptions record is not distributed.
    # The pin's path and digest are checked above; the record itself is checked
    # against the digest whenever it is present.
    if not source_path.is_relative_to(_ROOT.resolve()):
        raise ValueError("MECH-01A source pin is outside the repository")
    if source_path.is_file() and _canonical_lf_sha256(source_path) != pin["sha256"]:
        raise ValueError("MECH-01A source pin does not match its canonical LF bytes")

    d1 = _mapping("D1", document.get("d1_common_shaft"))
    if (
        _float_value("D1 nominal", d1.get("nominal_rpm")) != 11.0
        or _float_pair("D1 operating", d1.get("operating_envelope_rpm")) != (8.0, 12.0)
        or _float_pair("D1 design", d1.get("design_envelope_rpm")) != (4.0, 15.0)
        or _float_value("D1 low alarm", d1.get("low_alarm_rpm")) != 8.0
        or _float_value("D1 low-speed trip", d1.get("low_speed_trip_rpm")) != 4.0
        or d1.get("retired_values_rpm") != [3.0, 33.0, 60.0]
        or d1.get("retired_gui_range_rpm") != [0.0, 10.0]
        or d1.get("authenticated_machine_datum") is not False
    ):
        raise ValueError("MECH-01A D1 values differ from the owner ruling")
    speed_states = _mapping(
        "D0 schema-required speed states",
        d1.get("schema_required_speed_states"),
    )
    for name, expected in (
        ("crawl_rpm", (4.0, 1.0, 4.0)),
        ("high_speed_trip_threshold_rpm", (15.0, 15.0, 20.0)),
        ("post_restart_target_rpm", (8.0, 8.0, 12.0)),
    ):
        row = _mapping(f"D0 {name}", speed_states.get(name))
        observed = (
            _float_value(f"D0 {name} nominal", row.get("nominal")),
            *_float_pair(f"D0 {name} bracket", row.get("bracket")),
        )
        if observed != expected:
            raise ValueError(f"MECH-01A D0 {name} bracket drifted")

    d2 = _mapping("D2", document.get("d2_post_gearbox_shaft_power"))
    if (
        _float_value("D2 nominal", d2.get("nominal_w")) != 50_000.0
        or _float_pair("D2 bracket", d2.get("bracket_w")) != (25_000.0, 165_000.0)
        or _float_value("D2 omission ceiling", d2.get("omission_ceiling_w")) != 65_000.0
        or d2.get("omission_authorized") is not False
        or d2.get("ledger_row_required") is not True
        or d2.get("authenticated_machine_datum") is not False
    ):
        raise ValueError("MECH-01A D2 values differ from the owner ruling")

    d3 = _mapping("D3", document.get("d3_dissipation_partition"))
    if (
        d3.get("nominal_boundary_partition") != "ALL_TO_MEAL"
        or d3.get("conservative_boundary_endpoint") != "ALL_TO_WALL"
        or d3.get("nominal_within_meal_distribution") != "WET_HOLDUP_PROPORTIONAL_BY_TRAY"
        or d3.get("sensitivity_within_meal_distribution") != "BOTTOM_HEAVY_BY_SOLID_ORDER_RANK"
        or d3.get("bottom_heavy_rank_multipliers_by_solid_order") != [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        or d3.get("wall_endpoint_is_physical_claim") is not False
        or d3.get("authenticated_machine_datum") is not False
    ):
        raise ValueError("MECH-01A D3 partition or bracket drifted")

    d4 = _mapping("D4", document.get("d4_arm_dispersion"))
    delta = _mapping("D4 delta_z_ex", d4.get("delta_z_ex_m"))
    phi = _mapping("D4 phi_ex", d4.get("phi_ex"))
    state = _mapping("D4 f_state", d4.get("f_state"))
    if (
        _float_value("D4 delta nominal", delta.get("nominal")) != 0.17
        or _float_pair("D4 delta bracket", delta.get("bracket")) != (0.09, 0.33)
        or _float_value("D4 phi nominal", phi.get("nominal")) != 0.5
        or _float_pair("D4 phi bracket", phi.get("bracket")) != (0.3, 0.7)
        or _float_value("D4 state nominal", state.get("nominal")) != 1.0
        or _float_pair("D4 state bracket", state.get("bracket")) != (1.0, 1.0)
        or any(
            d4.get(flag) is not False
            for flag in (
                "directed_drift_selected",
                "saturation_law_selected",
                "state_modifier_identified",
                "authenticated_machine_datum",
            )
        )
    ):
        raise ValueError("MECH-01A D4 closure differs from the ruled bracket")

    d5 = _mapping("D5", document.get("d5_wet_bed_density"))
    if (
        _float_value("D5 nominal", d5.get("nominal_kg_wet_per_m3")) != 600.0
        or _float_pair("D5 bracket", d5.get("bracket_kg_wet_per_m3")) != (480.0, 640.0)
        or _float_value("D5 ceiling", d5.get("tapped_ceiling_kg_wet_per_m3")) != 720.0
        or d5.get("selected_machine_measurement") is not False
        or d5.get("authenticated_machine_datum") is not False
    ):
        raise ValueError("MECH-01A D5 density differs from the owner ruling")

    d6 = _mapping("D6", document.get("d6_discharge_configuration"))
    if (
        d6.get("middle_tray_configuration") != "COMMON_SHAFT_SWEEP_CHUTE"
        or d6.get("sparge_tray_configuration") != "INDEPENDENT_VARIABLE_SPEED_ROTARY_VALVE"
        or d6.get("draw_footprint_selected") is not False
        or d6.get("draw_flow_law_selected") is not False
        or d6.get("authenticated_machine_datum") is not False
    ):
        raise ValueError("MECH-01A D6 configuration differs from the owner ruling")

    trays = document.get("trays")
    if type(trays) is not list or len(trays) != 6:
        raise TypeError("MECH-01A requires six exact tray records")
    if tuple(row.get("tray_id") for row in trays if type(row) is dict) != (
        QSC_MECH01_SOLID_TRAY_ORDER
    ):
        raise ValueError("MECH-01A tray order differs from the six-tray solid order")
    for index, row in enumerate(trays, 1):
        row = _mapping("MECH-01A tray", row)
        if row.get("solid_order_rank") != index:
            raise ValueError("MECH-01A tray ranks must be contiguous in solid order")
        is_sp1 = row.get("tray_id") == "SP1"
        expected = (
            ("ROTARY_AIRLOCK", "INDEPENDENT_DRIVE", "INDEPENDENT_ACTUATOR")
            if is_sp1
            else (
                "INTERTRAY_SWEEP_PORT",
                "COMMON_SHAFT_SWEEP_CONVEYED",
                "NO_SEPARATE_ACTUATOR",
            )
        )
        observed = (
            row.get("draw_device_kind"),
            row.get("conveyance_domain"),
            row.get("actuation_domain"),
        )
        if observed != expected:
            raise ValueError("MECH-01A tray draw domains differ from ruled D6")

    claims = _mapping("MECH-01A claim boundary", document.get("claim_boundary"))
    expected_claims = {
        "f2",
        "f3",
        "f4",
        "f5",
        "f6",
        "f7",
        "f8a",
        "f8b",
        "qsc_10",
        "physically_qualifying",
        "plant_predictive",
        "industry_feasible_operating_envelope",
        "production_wired",
    }
    if set(claims) != expected_claims or any(value is not False for value in claims.values()):
        raise ValueError("MECH-01A claim boundary must be exact and wholly false")

    def visit(value: Any) -> None:
        if type(value) is dict:
            for key, child in value.items():
                if key == "authenticated_machine_datum" and child is not False:
                    raise ValueError("MECH-01A contains a non-false machine-datum flag")
                visit(child)
        elif type(value) is list:
            for child in value:
                visit(child)

    visit(document)


_require_configuration_semantics(_CONFIG)


def qsc_mech01_configuration_errors(
    document: Any,
    *,
    canonical_lf_sha256: str,
) -> tuple[str, ...]:
    """Return fail-closed reasons for an external MECH-01A document."""

    errors: list[str] = []
    if canonical_lf_sha256 != QSC_MECH01_CONFIG_CANONICAL_LF_SHA256:
        errors.append("configuration canonical-LF digest differs from the source pin")
    if type(document) is not dict:
        return (*errors, "configuration is not an exact mapping")
    if document != _CONFIG:
        errors.append("configuration values differ from the instantiated runtime authority")
    try:
        _require_configuration_semantics(document)
    except (KeyError, TypeError, ValueError) as error:
        errors.append(f"{type(error).__name__}: {error}")
    return tuple(errors)


_PIN = _CONFIG["source_pins"]["engineering_assumptions_ruling"]
_PINNED_RULING_DIGEST = "sha256:" + _PIN["sha256"]


def _provenance(section: str, uncertainty_basis_id: str) -> mechanical.MechanicalFieldProvenance:
    return mechanical.MechanicalFieldProvenance(
        evidence_class=mechanical.MechanicalEvidenceClass.PROVISIONAL_ENGINEERING,
        source_id="GT-PS-2-engineering-assumptions-ruling-record-2026-08-10",
        source_content_sha256=_PINNED_RULING_DIGEST,
        source_locator=f"{_PIN['path']}#{section}",
        acquired_at_utc="2026-08-10T00:00:00Z",
        uncertainty_basis_id=uncertainty_basis_id,
    )


def _scalar(
    quantity: mechanical.MechanicalQuantity,
    value: float,
    lower: float,
    upper: float,
    provenance: mechanical.MechanicalFieldProvenance,
) -> mechanical.MechanicalScalarField:
    unit = {
        mechanical.MechanicalQuantity.SHAFT_SPEED_ZERO: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_SPEED_CRAWL: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_SPEED_MINIMUM_STABLE: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_SPEED_NOMINAL: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_SPEED_MAXIMUM_CONTINUOUS: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_SPEED_TRIP: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_SPEED_RESTART: mechanical.MechanicalUnit.RPM,
        mechanical.MechanicalQuantity.SHAFT_RAMP_UP_LIMIT: mechanical.MechanicalUnit.RPM_PER_SECOND,
        mechanical.MechanicalQuantity.SHAFT_RAMP_DOWN_LIMIT: (
            mechanical.MechanicalUnit.RPM_PER_SECOND
        ),
        mechanical.MechanicalQuantity.ARM_INNER_RADIUS: mechanical.MechanicalUnit.METRE,
        mechanical.MechanicalQuantity.ARM_OUTER_RADIUS: mechanical.MechanicalUnit.METRE,
        mechanical.MechanicalQuantity.BLADE_HEIGHT: mechanical.MechanicalUnit.METRE,
        mechanical.MechanicalQuantity.BLADE_WIDTH: mechanical.MechanicalUnit.METRE,
        mechanical.MechanicalQuantity.BLADE_PITCH: mechanical.MechanicalUnit.RADIAN,
        mechanical.MechanicalQuantity.BLADE_RAKE: mechanical.MechanicalUnit.RADIAN,
        mechanical.MechanicalQuantity.TIP_CLEARANCE: mechanical.MechanicalUnit.METRE,
        mechanical.MechanicalQuantity.SWEPT_VOLUME_PER_REVOLUTION: (
            mechanical.MechanicalUnit.CUBIC_METRE_PER_REVOLUTION
        ),
        mechanical.MechanicalQuantity.GEARBOX_RATIO: mechanical.MechanicalUnit.DIMENSIONLESS,
        mechanical.MechanicalQuantity.RATED_POST_GEARBOX_POWER: mechanical.MechanicalUnit.WATT,
        mechanical.MechanicalQuantity.RATED_POST_GEARBOX_TORQUE: (
            mechanical.MechanicalUnit.NEWTON_METRE
        ),
        mechanical.MechanicalQuantity.DRIVETRAIN_EFFICIENCY: (
            mechanical.MechanicalUnit.DIMENSIONLESS
        ),
    }[quantity]
    return mechanical.MechanicalScalarField(
        quantity=quantity,
        unit=unit,
        value=float(value),
        lower_bound=float(lower),
        upper_bound=float(upper),
        provenance=provenance,
    )


def _geometry_row(name: str) -> tuple[float, float, float]:
    row = _mapping(name, _CONFIG["schema_required_geometry_d0"][name])
    lower, upper = _float_pair(f"{name} bracket", row["bracket"])
    return _float_value(f"{name} nominal", row["nominal"]), lower, upper


def _draw_material(tray_id: str, draw_device_id: str) -> tuple[bytes, bytes]:
    geometry = (f"D6|{tray_id}|{draw_device_id}|IDENTITY_ONLY|NO_FOOTPRINT|NO_FLOW_LAW").encode(
        "utf-8"
    )
    uncertainty = (
        f"D0-D6|{tray_id}|DRAW_GEOMETRY_UNKNOWN|REPLACE_WITH_AUTHENTICATED_DRAWING"
    ).encode("utf-8")
    return geometry, uncertainty


def _build_common_shaft_authority() -> mechanical.CommonShaftMechanicalAuthority:
    machine_doc = _CONFIG["machine_identity"]
    d1 = _CONFIG["d1_common_shaft"]
    d2 = _CONFIG["d2_post_gearbox_shaft_power"]
    d0_drive = d2["schema_required_drive_train"]
    history_doc = _CONFIG["accepted_shaft_history"]
    d0_geometry = _CONFIG["schema_required_geometry_d0"]
    d0_speed = d1["schema_required_speed_states"]

    d0_provenance = _provenance("D0", "D0-explicit-schema-occupancy-bracket")
    d1_provenance = _provenance("D1", "D1-owner-ruled-speed-envelope")
    d2_provenance = _provenance("D2", "D2-owner-ruled-power-bracket")
    d3_provenance = _provenance(
        "D3",
        "D3-owner-ruled-dissipation-partition-bed-shear-clearance",
    )
    d6_provenance = _provenance("D6", "D6-declared-discharge-configuration")

    machine = mechanical.MachineConfigurationIdentity(
        unit_id=machine_doc["unit_id"],
        unit_class_id=machine_doc["unit_class_id"],
        serial_id=machine_doc["serial_id"],
        configuration_id=machine_doc["configuration_id"],
        configuration_revision=machine_doc["configuration_revision"],
        provenance=d0_provenance,
    )
    direction = mechanical.ShaftDirectionAuthority(
        unit_id=machine.unit_id,
        configuration_id=machine.configuration_id,
        configuration_revision=machine.configuration_revision,
        direction_reference_frame_id="QSC-REFERENCE-SHAFT-AXIS-Z-UP",
        rotation_sense=mechanical.ShaftRotationSense.POSITIVE_ABOUT_DECLARED_AXIS,
        provenance=d0_provenance,
    )
    history_material = history_doc["history_material"].encode("utf-8")
    history = mechanical.AcceptedShaftHistoryReference(
        unit_id=machine.unit_id,
        configuration_id=machine.configuration_id,
        configuration_revision=machine.configuration_revision,
        history_id=history_doc["history_id"],
        history_status=mechanical.ShaftHistoryStatus.ACCEPTED_ESTIMATED,
        sample_semantics=(
            mechanical.ShaftHistorySampleSemantics.NONNEGATIVE_RPM_MAGNITUDE_IN_DECLARED_SENSE
        ),
        speed_unit=mechanical.MechanicalUnit.RPM,
        direction_authority_digest=direction.definition_digest,
        sensor_channel_id=history_doc["sensor_channel_id"],
        controller_channel_id=history_doc["controller_channel_id"],
        history_schema_id=history_doc["history_schema_id"],
        history_content_length_bytes=len(history_material),
        history_content_sha256=_sha256_bytes(history_material),
        accepted_interval_start_utc=history_doc["interval_start_utc"],
        accepted_interval_end_utc=history_doc["interval_end_utc"],
        provenance=d0_provenance,
    )

    ramp = d1["schema_required_ramp_rate_rpm_per_s"]
    ramp_lower, ramp_upper = _float_pair("D0 ramp bracket", ramp["bracket"])
    ramp_nominal = _float_value("D0 ramp nominal", ramp["nominal"])
    crawl_nominal = _float_value("D0 crawl nominal", d0_speed["crawl_rpm"]["nominal"])
    crawl_lower, crawl_upper = _float_pair("D0 crawl bracket", d0_speed["crawl_rpm"]["bracket"])
    high_trip_nominal = _float_value(
        "D0 high-speed trip nominal",
        d0_speed["high_speed_trip_threshold_rpm"]["nominal"],
    )
    high_trip_lower, high_trip_upper = _float_pair(
        "D0 high-speed trip bracket",
        d0_speed["high_speed_trip_threshold_rpm"]["bracket"],
    )
    restart_nominal = _float_value(
        "D0 restart nominal", d0_speed["post_restart_target_rpm"]["nominal"]
    )
    restart_lower, restart_upper = _float_pair(
        "D0 restart bracket", d0_speed["post_restart_target_rpm"]["bracket"]
    )
    speed = mechanical.ShaftSpeedEnvelope(
        trip_speed_semantics=mechanical.ShaftTripSpeedSemantics.TRIP_THRESHOLD_MAGNITUDE,
        restart_speed_semantics=(
            mechanical.ShaftRestartSpeedSemantics.POST_RESTART_TARGET_MAGNITUDE
        ),
        semantics_provenance=d0_provenance,
        zero_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_ZERO,
            0.0,
            0.0,
            0.0,
            d0_provenance,
        ),
        crawl_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_CRAWL,
            crawl_nominal,
            crawl_lower,
            crawl_upper,
            d0_provenance,
        ),
        minimum_stable_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_MINIMUM_STABLE,
            8.0,
            8.0,
            8.0,
            d1_provenance,
        ),
        nominal_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_NOMINAL,
            11.0,
            11.0,
            11.0,
            d1_provenance,
        ),
        maximum_continuous_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_MAXIMUM_CONTINUOUS,
            12.0,
            12.0,
            12.0,
            d1_provenance,
        ),
        trip_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_TRIP,
            high_trip_nominal,
            high_trip_lower,
            high_trip_upper,
            d0_provenance,
        ),
        restart_rpm=_scalar(
            mechanical.MechanicalQuantity.SHAFT_SPEED_RESTART,
            restart_nominal,
            restart_lower,
            restart_upper,
            d0_provenance,
        ),
        ramp_up_limit_rpm_per_s=_scalar(
            mechanical.MechanicalQuantity.SHAFT_RAMP_UP_LIMIT,
            ramp_nominal,
            ramp_lower,
            ramp_upper,
            d0_provenance,
        ),
        ramp_down_limit_rpm_per_s=_scalar(
            mechanical.MechanicalQuantity.SHAFT_RAMP_DOWN_LIMIT,
            ramp_nominal,
            ramp_lower,
            ramp_upper,
            d0_provenance,
        ),
    )

    power_lower, power_upper = _float_pair("D2 power bracket", d2["bracket_w"])
    power_nominal = _float_value("D2 power nominal", d2["nominal_w"])
    nominal_omega = 2.0 * math.pi * 11.0 / 60.0
    minimum_omega = 2.0 * math.pi * 8.0 / 60.0
    maximum_omega = 2.0 * math.pi * 12.0 / 60.0
    torque_nominal = power_nominal / nominal_omega
    torque_lower = power_lower / maximum_omega
    torque_upper = power_upper / minimum_omega
    ratio_doc = d0_drive["gearbox_ratio_motor_speed_over_post_gearbox_speed"]
    ratio_lower, ratio_upper = _float_pair("D0 gearbox ratio", ratio_doc["bracket"])
    efficiency_doc = d0_drive["efficiency_electrical_input_to_post_gearbox_output"]
    efficiency_lower, efficiency_upper = _float_pair(
        "D0 drive efficiency", efficiency_doc["bracket"]
    )
    drive = mechanical.DriveTrainAuthority(
        unit_id=machine.unit_id,
        configuration_id=machine.configuration_id,
        configuration_revision=machine.configuration_revision,
        motor_id="UNAUTHENTICATED-MOTOR-NOT-ASSIGNED",
        gearbox_id="UNAUTHENTICATED-GEARBOX-NOT-ASSIGNED",
        drive_configuration_id="D0-SCHEMA-OCCUPANCY-NOT-A-NAMEPLATE",
        gearbox_ratio_convention=(
            mechanical.GearboxRatioConvention.MOTOR_SPEED_DIVIDED_BY_POST_GEARBOX_SHAFT_SPEED
        ),
        efficiency_boundary=(
            mechanical.DriveTrainEfficiencyBoundary.ELECTRICAL_INPUT_TO_POST_GEARBOX_MECHANICAL_OUTPUT
        ),
        identity_provenance=d0_provenance,
        convention_provenance=d0_provenance,
        gearbox_ratio=_scalar(
            mechanical.MechanicalQuantity.GEARBOX_RATIO,
            _float_value("D0 gearbox nominal", ratio_doc["nominal"]),
            ratio_lower,
            ratio_upper,
            d0_provenance,
        ),
        rated_post_gearbox_power_w=_scalar(
            mechanical.MechanicalQuantity.RATED_POST_GEARBOX_POWER,
            power_nominal,
            power_lower,
            power_upper,
            d2_provenance,
        ),
        rated_post_gearbox_torque_n_m=_scalar(
            mechanical.MechanicalQuantity.RATED_POST_GEARBOX_TORQUE,
            torque_nominal,
            torque_lower,
            torque_upper,
            d2_provenance,
        ),
        efficiency=_scalar(
            mechanical.MechanicalQuantity.DRIVETRAIN_EFFICIENCY,
            _float_value("D0 efficiency nominal", efficiency_doc["nominal"]),
            efficiency_lower,
            efficiency_upper,
            d0_provenance,
        ),
    )

    arm_count_doc = d0_geometry["arm_count"]
    arm_count_bounds = arm_count_doc["bracket"]
    geometry_quantities = (
        ("arm_inner_radius_m", mechanical.MechanicalQuantity.ARM_INNER_RADIUS),
        ("arm_outer_radius_m", mechanical.MechanicalQuantity.ARM_OUTER_RADIUS),
        ("blade_height_m", mechanical.MechanicalQuantity.BLADE_HEIGHT),
        ("blade_width_m", mechanical.MechanicalQuantity.BLADE_WIDTH),
        ("blade_pitch_rad", mechanical.MechanicalQuantity.BLADE_PITCH),
        ("blade_rake_rad", mechanical.MechanicalQuantity.BLADE_RAKE),
        ("tip_clearance_m", mechanical.MechanicalQuantity.TIP_CLEARANCE),
        (
            "swept_volume_per_revolution_m3",
            mechanical.MechanicalQuantity.SWEPT_VOLUME_PER_REVOLUTION,
        ),
    )
    trays: list[mechanical.TrayMechanicalAuthority] = []
    for tray_doc in _CONFIG["trays"]:
        tray_id = tray_doc["tray_id"]
        scalars = {
            name: _scalar(
                quantity,
                *_geometry_row(name),
                d3_provenance if name == "tip_clearance_m" else d0_provenance,
            )
            for name, quantity in geometry_quantities
        }
        arm = mechanical.ArmGeometryAuthority(
            physical_tray_id=tray_id,
            arm_geometry_id=f"QSC-MECH01-D0-ARM-GEOMETRY-{tray_id}",
            coordinate_frame_id="QSC-REFERENCE-SHAFT-AXIS-Z-UP",
            identity_provenance=d0_provenance,
            arm_count=mechanical.MechanicalIntegerField(
                quantity=mechanical.MechanicalQuantity.ARM_COUNT,
                unit=mechanical.MechanicalUnit.COUNT,
                value=int(arm_count_doc["nominal"]),
                lower_bound=int(arm_count_bounds[0]),
                upper_bound=int(arm_count_bounds[1]),
                provenance=d0_provenance,
            ),
            inner_radius_m=scalars["arm_inner_radius_m"],
            outer_radius_m=scalars["arm_outer_radius_m"],
            blade_height_m=scalars["blade_height_m"],
            blade_width_m=scalars["blade_width_m"],
            blade_pitch_rad=scalars["blade_pitch_rad"],
            blade_rake_rad=scalars["blade_rake_rad"],
            tip_clearance_m=scalars["tip_clearance_m"],
            swept_volume_per_revolution_m3=scalars["swept_volume_per_revolution_m3"],
        )
        kind = mechanical.DrawDeviceKind[tray_doc["draw_device_kind"]]
        conveyance = mechanical.SolidsConveyanceDomain[tray_doc["conveyance_domain"]]
        actuation = mechanical.DrawActuationDomain[tray_doc["actuation_domain"]]
        geometry_material, uncertainty_material = _draw_material(
            tray_id, tray_doc["draw_device_id"]
        )
        independent_conveyance = (
            _sha256_text(f"D6|{tray_id}|INDEPENDENT-ROTARY-CONVEYANCE-NOT-CALIBRATED")
            if conveyance is mechanical.SolidsConveyanceDomain.INDEPENDENT_DRIVE
            else None
        )
        independent_actuation = (
            _sha256_text(f"D6|{tray_id}|INDEPENDENT-ACTUATION-NOT-CALIBRATED")
            if actuation is mechanical.DrawActuationDomain.INDEPENDENT_ACTUATOR
            else None
        )
        draw = mechanical.DrawGeometryReference(
            physical_tray_id=tray_id,
            draw_device_id=tray_doc["draw_device_id"],
            device_kind=kind,
            conveyance_domain=conveyance,
            actuation_domain=actuation,
            draw_geometry_schema_id="qsc-mech01-d6-identity-only-no-footprint-v1",
            draw_geometry_content_length_bytes=len(geometry_material),
            draw_geometry_content_sha256=_sha256_bytes(geometry_material),
            draw_geometry_uncertainty_content_sha256=_sha256_bytes(uncertainty_material),
            independent_conveyance_authority_digest=independent_conveyance,
            independent_actuation_authority_digest=independent_actuation,
            geometry_provenance=d6_provenance,
            conveyance_domain_provenance=d6_provenance,
            actuation_domain_provenance=d6_provenance,
        )
        trays.append(
            mechanical.TrayMechanicalAuthority(
                physical_tray_id=tray_id,
                arm_geometry=arm,
                draw_geometry=draw,
            )
        )

    # PR-08 requires canonical lexical tray order.  The runtime binding below
    # retains the physical solid order independently.
    authority = mechanical.CommonShaftMechanicalAuthority(
        schema_id=mechanical.MECHANICAL_AUTHORITY_SCHEMA_ID,
        authority_revision=1,
        machine=machine,
        shaft_direction=direction,
        accepted_shaft_history=history,
        speed_envelope=speed,
        drive_train=drive,
        trays=tuple(sorted(trays, key=lambda tray: tray.physical_tray_id)),
    )
    mechanical.validate_common_shaft_mechanical_authority(authority)
    return authority


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCTrayMechanicalBinding:
    physical_tray_id: str
    common_shaft_definition_digest: str
    accepted_history_definition_digest: str
    tray_mechanical_definition_digest: str
    arm_geometry_definition_digest: str
    k_authority_definition_digest: str
    draw_geometry_definition_digest: str

    authenticated_machine_data: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("MECH-01A tray binding", self.physical_tray_id)
        for name in (
            "common_shaft_definition_digest",
            "accepted_history_definition_digest",
            "tray_mechanical_definition_digest",
            "arm_geometry_definition_digest",
            "k_authority_definition_digest",
            "draw_geometry_definition_digest",
        ):
            _require_sha256(name, getattr(self, name))

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="qsc-tray-mechanical-binding",
            parts=(
                self.physical_tray_id,
                self.common_shaft_definition_digest,
                self.accepted_history_definition_digest,
                self.tray_mechanical_definition_digest,
                self.arm_geometry_definition_digest,
                self.k_authority_definition_digest,
                self.draw_geometry_definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCTrayShaftWorkFraction:
    physical_tray_id: str
    wet_holdup_proportional_fraction: float
    bottom_heavy_fraction: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("shaft-work tray", self.physical_tray_id)
        _require_float(
            "shaft-work fractions",
            self.wet_holdup_proportional_fraction,
            self.bottom_heavy_fraction,
            positive=True,
        )
        if self.wet_holdup_proportional_fraction > 1.0 or self.bottom_heavy_fraction > 1.0:
            raise ValueError("shaft-work fractions must not exceed one")

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="qsc-tray-shaft-work-fraction",
            parts=(
                self.physical_tray_id,
                self.wet_holdup_proportional_fraction,
                self.bottom_heavy_fraction,
            ),
        )


def _normalized_fractions(values: tuple[float, ...]) -> tuple[float, ...]:
    _require_float("allocation weights", *values, positive=True)
    total = math.fsum(values)
    fractions = tuple(value / total for value in values[:-1])
    closing = 1.0 - math.fsum(fractions)
    return (*fractions, closing)


def _build_k_authorities(
    common: mechanical.CommonShaftMechanicalAuthority,
    *,
    shaft_rpm: float = 11.0,
) -> tuple[k_scheduler.EngineeringKEventAuthority, ...]:
    d4 = _CONFIG["d4_arm_dispersion"]
    density = _geometry_row("k2_reference_dry_matter_density_kg_m3")[0]
    by_tray = {tray.physical_tray_id: tray for tray in common.trays}
    # The canonical 11 r/min segment id (and therefore every canonical digest)
    # must stay byte-identical; ruled variants carry their exact repr.
    rpm_label = "11" if shaft_rpm == 11.0 else repr(shaft_rpm)
    segment = _arm_oracle.ManufacturedAcceptedShaftSegment(
        common_shaft_definition_digest=common.definition_digest,
        accepted_history_definition_digest=common.accepted_shaft_history.definition_digest,
        segment_id=f"QSC-MECH01-D1-CONSTANT-{rpm_label}RPM-ENGINEERING-SEGMENT",
        start_event_id="QSC-MECH01-HISTORY-START",
        end_event_id="QSC-MECH01-HISTORY-END",
        start_event_kind=_arm_oracle.ShaftBoundaryKind.CONTINUOUS,
        end_event_kind=_arm_oracle.ShaftBoundaryKind.CONTINUOUS,
        time_start_s=0.0,
        time_end_s=31_536_000.0,
        rpm_start=shaft_rpm,
        rpm_end=shaft_rpm,
    )
    authorities: list[k_scheduler.EngineeringKEventAuthority] = []
    for tray_doc in _CONFIG["trays"]:
        tray_id = tray_doc["tray_id"]
        depth = _float_value(f"{tray_id} depth", tray_doc["loaded_depth_m"])
        diameter = _float_value(f"{tray_id} diameter", tray_doc["diameter_m"])
        half_depth = 0.5 * depth
        area = math.pi * (0.5 * diameter) ** 2
        half_holdup = density * area * half_depth
        geometry_digest = by_tray[tray_id].arm_geometry.definition_digest
        mesh = _arm_oracle.PhysicalHeightMesh(
            physical_tray_id=tray_id,
            mesh_id=f"QSC-MECH01-{tray_id}-K2-D0-NUMERICAL-MESH",
            tray_geometry_definition_digest=geometry_digest,
            cells=(
                _arm_oracle.PhysicalHeightCell(
                    vertical_layer_id=VerticalLayerId(1),
                    z_lower_m=0.0,
                    z_upper_m=half_depth,
                    dry_matter_holdup_kg=half_holdup,
                ),
                _arm_oracle.PhysicalHeightCell(
                    vertical_layer_id=VerticalLayerId(2),
                    z_lower_m=half_depth,
                    z_upper_m=depth,
                    dry_matter_holdup_kg=half_holdup,
                ),
            ),
            faces=(
                _arm_oracle.PhysicalHeightFace(
                    lower_vertical_layer_id=VerticalLayerId(1),
                    upper_vertical_layer_id=VerticalLayerId(2),
                    face_area_m2=area,
                    dry_matter_bulk_density_kg_m3=density,
                    gradient_distance_m=half_depth,
                ),
            ),
        )
        closure = _arm_oracle.ManufacturedArmClosure(
            closure_id=f"QSC-MECH01-D4-{tray_id}-DECLARED-ENGINEERING-CLOSURE",
            closure_source_digest=_PINNED_RULING_DIGEST,
            tray_geometry_definition_digest=geometry_digest,
            phi_ex=_float_value("D4 phi nominal", d4["phi_ex"]["nominal"]),
            delta_z_ex_m=_float_value("D4 delta nominal", d4["delta_z_ex_m"]["nominal"]),
            f_state=_float_value("D4 state nominal", d4["f_state"]["nominal"]),
        )
        authorities.append(
            k_scheduler.EngineeringKEventAuthority(
                mesh=mesh,
                closure=closure,
                shaft_segments=(segment,),
                mode=_arm_oracle.ArmDispersionMode.ACTIVE_BASELINE,
                scheduler_seed_id=f"QSC-MECH01-{tray_id}-K-EVENT-SEED-V1",
            )
        )
    return tuple(authorities)


def _build_work_fractions() -> tuple[QSCTrayShaftWorkFraction, ...]:
    wet_weights: list[float] = []
    bottom_weights: list[float] = []
    multipliers = _CONFIG["d3_dissipation_partition"][
        "bottom_heavy_rank_multipliers_by_solid_order"
    ]
    for tray_doc, multiplier in zip(_CONFIG["trays"], multipliers, strict=True):
        diameter = _float_value("tray diameter", tray_doc["diameter_m"])
        depth = _float_value("tray loaded depth", tray_doc["loaded_depth_m"])
        wet_volume = math.pi * (0.5 * diameter) ** 2 * depth
        wet_weights.append(wet_volume)
        bottom_weights.append(wet_volume * _float_value("bottom-heavy multiplier", multiplier))
    wet = _normalized_fractions(tuple(wet_weights))
    bottom = _normalized_fractions(tuple(bottom_weights))
    return tuple(
        QSCTrayShaftWorkFraction(
            physical_tray_id=tray_id,
            wet_holdup_proportional_fraction=wet_fraction,
            bottom_heavy_fraction=bottom_fraction,
        )
        for tray_id, wet_fraction, bottom_fraction in zip(
            QSC_MECH01_SOLID_TRAY_ORDER, wet, bottom, strict=True
        )
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMechanicalRuntimeAuthority:
    schema_id: str
    schema_revision: int
    configuration_relative_path: str
    configuration_canonical_lf_sha256: str
    parameter_class: QSCMechanicalParameterClass
    common_shaft: mechanical.CommonShaftMechanicalAuthority
    k_authorities: tuple[k_scheduler.EngineeringKEventAuthority, ...]
    tray_bindings: tuple[QSCTrayMechanicalBinding, ...]
    shaft_work_fractions: tuple[QSCTrayShaftWorkFraction, ...]
    design_envelope_rpm: tuple[float, float]
    low_alarm_rpm: float
    low_speed_trip_rpm: float
    retired_values_rpm: tuple[float, float, float]
    retired_gui_range_rpm: tuple[float, float]
    shaft_power_bracket_w: tuple[float, float, float]
    nominal_boundary_partition: ShaftWorkBoundaryPartition
    conservative_boundary_endpoint: ShaftWorkBoundaryPartition
    nominal_within_meal_distribution: WithinMealDistribution
    sensitivity_within_meal_distribution: WithinMealDistribution
    delta_z_ex_bracket_m: tuple[float, float, float]
    phi_ex_bracket: tuple[float, float, float]
    f_state: float
    wet_bed_density_bracket_kg_m3: tuple[float, float, float, float]
    #: Ruled 2026-08-30: the declared constant operating speed of this
    #: authority's shaft history.  11.0 is the canonical D1 nominal; other
    #: values are admissible only inside the D1 operating envelope.  This is
    #: a projection of the K-authority shaft segments (cross-checked below),
    #: not independent state, so it does not enter the definition digest.
    declared_constant_shaft_rpm: float = 11.0

    mech_01a_runtime_authority_instantiated: ClassVar[bool] = True
    evaluator_wired: ClassVar[bool] = False
    shaft_work_deposition_implemented: ClassVar[bool] = False
    draw_flow_law_selected: ClassVar[bool] = False
    draw_footprint_selected: ClassVar[bool] = False
    authenticated_machine_data: ClassVar[bool] = False
    source_authenticated: ClassVar[bool] = False
    f2: ClassVar[bool] = False
    f3: ClassVar[bool] = False
    f4: ClassVar[bool] = False
    f5: ClassVar[bool] = False
    f6: ClassVar[bool] = False
    f7: ClassVar[bool] = False
    f8a: ClassVar[bool] = False
    f8b: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industry_feasible_operating_envelope: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.schema_id != QSC_MECH01_SCHEMA_ID:
            raise ValueError("runtime authority has the wrong MECH-01A schema")
        if self.schema_revision != QSC_MECH01_SCHEMA_REVISION:
            raise ValueError("runtime authority has the wrong MECH-01A revision")
        if self.configuration_relative_path != QSC_MECH01_CONFIG_RELATIVE_PATH:
            raise ValueError("runtime authority binds another configuration path")
        if self.configuration_canonical_lf_sha256 != (QSC_MECH01_CONFIG_CANONICAL_LF_SHA256):
            raise ValueError("runtime authority binds another configuration digest")
        if self.parameter_class is not (
            QSCMechanicalParameterClass.DECLARED_ENGINEERING_ASSUMPTION
        ):
            raise ValueError("runtime authority must retain declared-assumption standing")
        if type(self.common_shaft) is not mechanical.CommonShaftMechanicalAuthority:
            raise TypeError("runtime common shaft has a foreign exact type")
        if mechanical.validate_common_shaft_mechanical_authority(self.common_shaft) != (
            self.common_shaft.definition_digest
        ):
            raise ValueError("runtime common-shaft identity did not reconstruct")
        if len(self.common_shaft.trays) != 6:
            raise ValueError("runtime common shaft must contain exactly six trays")

        speed = self.common_shaft.speed_envelope
        if (
            speed.nominal_rpm.value != 11.0
            or speed.minimum_stable_rpm.value != 8.0
            or speed.maximum_continuous_rpm.value != 12.0
            or speed.trip_rpm.value != 15.0
            or speed.restart_rpm.value != 8.0
            or self.design_envelope_rpm != (4.0, 15.0)
            or self.low_alarm_rpm != 8.0
            or self.low_speed_trip_rpm != 4.0
            or self.retired_values_rpm != (3.0, 33.0, 60.0)
            or self.retired_gui_range_rpm != (0.0, 10.0)
        ):
            raise ValueError("runtime authority differs from ruled D1 or re-admits retired RPM")
        rpm = self.declared_constant_shaft_rpm
        low_rpm, high_rpm = QSC_MECH01_OPERATING_RPM_ENVELOPE
        if (
            type(rpm) is not float
            or not math.isfinite(rpm)
            or not low_rpm <= rpm <= high_rpm
            or rpm in self.retired_values_rpm
        ):
            raise ValueError(
                "declared constant shaft speed must sit inside the ruled D1 "
                "operating envelope and off the retired values"
            )
        power = self.common_shaft.drive_train.rated_post_gearbox_power_w
        if (
            self.shaft_power_bracket_w != (25_000.0, 50_000.0, 165_000.0)
            or (power.lower_bound, power.value, power.upper_bound) != self.shaft_power_bracket_w
        ):
            raise ValueError("runtime authority differs from ruled D2 power bracket")
        if (
            self.nominal_boundary_partition is not ShaftWorkBoundaryPartition.ALL_TO_MEAL
            or self.conservative_boundary_endpoint is not ShaftWorkBoundaryPartition.ALL_TO_WALL
            or self.nominal_within_meal_distribution
            is not WithinMealDistribution.WET_HOLDUP_PROPORTIONAL_BY_TRAY
            or self.sensitivity_within_meal_distribution
            is not WithinMealDistribution.BOTTOM_HEAVY_BY_SOLID_ORDER_RANK
        ):
            raise ValueError("runtime authority differs from ruled D3 partition brackets")
        if (
            self.delta_z_ex_bracket_m != (0.09, 0.17, 0.33)
            or self.phi_ex_bracket != (0.3, 0.5, 0.7)
            or self.f_state != 1.0
        ):
            raise ValueError("runtime authority differs from ruled D4 closure bracket")
        if self.wet_bed_density_bracket_kg_m3 != (480.0, 600.0, 640.0, 720.0):
            raise ValueError("runtime authority differs from ruled D5 density bracket")

        if type(self.k_authorities) is not tuple or len(self.k_authorities) != 6:
            raise TypeError("runtime authority requires six exact K authorities")
        if any(
            type(authority) is not k_scheduler.EngineeringKEventAuthority
            for authority in self.k_authorities
        ):
            raise TypeError("runtime K authority tuple contains a foreign exact type")
        if tuple(item.mesh.physical_tray_id for item in self.k_authorities) != (
            QSC_MECH01_SOLID_TRAY_ORDER
        ):
            raise ValueError("runtime K authorities differ from solid tray order")
        if type(self.tray_bindings) is not tuple or len(self.tray_bindings) != 6:
            raise TypeError("runtime authority requires six exact tray bindings")
        if any(type(item) is not QSCTrayMechanicalBinding for item in self.tray_bindings):
            raise TypeError("runtime tray binding tuple contains a foreign exact type")
        if tuple(item.physical_tray_id for item in self.tray_bindings) != (
            QSC_MECH01_SOLID_TRAY_ORDER
        ):
            raise ValueError("runtime tray bindings differ from solid tray order")

        common_digest = self.common_shaft.definition_digest
        history_digest = self.common_shaft.accepted_shaft_history.definition_digest
        tray_map = {tray.physical_tray_id: tray for tray in self.common_shaft.trays}
        tray_docs = {row["tray_id"]: row for row in _CONFIG["trays"]}
        reference_density = _geometry_row("k2_reference_dry_matter_density_kg_m3")[0]
        for k_authority, binding in zip(self.k_authorities, self.tray_bindings, strict=True):
            tray_id = binding.physical_tray_id
            tray = tray_map[tray_id]
            expected_binding = QSCTrayMechanicalBinding(
                physical_tray_id=tray_id,
                common_shaft_definition_digest=common_digest,
                accepted_history_definition_digest=history_digest,
                tray_mechanical_definition_digest=tray.definition_digest,
                arm_geometry_definition_digest=tray.arm_geometry.definition_digest,
                k_authority_definition_digest=k_authority.definition_digest,
                draw_geometry_definition_digest=tray.draw_geometry.definition_digest,
            )
            if binding != expected_binding:
                raise ValueError(
                    "runtime tray binding contains a foreign shaft/history/tray/K/draw digest"
                )
            if (
                k_authority.common_shaft_definition_digest != common_digest
                or k_authority.accepted_history_definition_digest != history_digest
            ):
                raise ValueError("runtime K authority binds a foreign shaft or history root")
            if (
                k_authority.mesh.tray_geometry_definition_digest
                != tray.arm_geometry.definition_digest
                or k_authority.closure.tray_geometry_definition_digest
                != tray.arm_geometry.definition_digest
            ):
                raise ValueError("runtime K authority binds a foreign tray geometry")
            if (
                k_authority.closure.phi_ex != 0.5
                or k_authority.closure.delta_z_ex_m != 0.17
                or k_authority.closure.f_state != 1.0
                or k_authority.closure.closure_source_digest != _PINNED_RULING_DIGEST
                or k_authority.mode is not _arm_oracle.ArmDispersionMode.ACTIVE_BASELINE
            ):
                raise ValueError("runtime K authority differs from D4 active baseline")
            segments = k_authority.shaft_segments
            if (
                len(segments) != 1
                or segments[0].rpm_start != self.declared_constant_shaft_rpm
                or segments[0].rpm_end != self.declared_constant_shaft_rpm
                or segments[0].time_start_s != 0.0
                or segments[0].time_end_s != 31_536_000.0
            ):
                raise ValueError(
                    "runtime K authority differs from the declared constant D1 history"
                )
            mesh = k_authority.mesh
            tray_doc = tray_docs[tray_id]
            depth = float(tray_doc["loaded_depth_m"])
            area = math.pi * (0.5 * float(tray_doc["diameter_m"])) ** 2
            if (
                len(mesh.cells) != 2
                or len(mesh.faces) != 1
                or mesh.cells[0].z_lower_m != 0.0
                or mesh.cells[0].z_upper_m != 0.5 * depth
                or mesh.cells[1].z_lower_m != 0.5 * depth
                or mesh.cells[1].z_upper_m != depth
                or mesh.faces[0].face_area_m2 != area
                or mesh.faces[0].dry_matter_bulk_density_kg_m3 != reference_density
                or mesh.faces[0].gradient_distance_m != 0.5 * depth
            ):
                raise ValueError("runtime K authority differs from its D0 K2 numerical mesh")
            draw = tray.draw_geometry
            is_sp1 = tray_id == "SP1"
            expected_draw = (
                (
                    mechanical.DrawDeviceKind.ROTARY_AIRLOCK,
                    mechanical.SolidsConveyanceDomain.INDEPENDENT_DRIVE,
                    mechanical.DrawActuationDomain.INDEPENDENT_ACTUATOR,
                )
                if is_sp1
                else (
                    mechanical.DrawDeviceKind.INTERTRAY_SWEEP_PORT,
                    mechanical.SolidsConveyanceDomain.COMMON_SHAFT_SWEEP_CONVEYED,
                    mechanical.DrawActuationDomain.NO_SEPARATE_ACTUATOR,
                )
            )
            if (draw.device_kind, draw.conveyance_domain, draw.actuation_domain) != expected_draw:
                raise ValueError("runtime draw authority differs from ruled D6")
            geometry_material, uncertainty_material = _draw_material(
                tray_id, tray_doc["draw_device_id"]
            )
            if draw.draw_geometry_content_sha256 != _sha256_bytes(
                geometry_material
            ) or draw.draw_geometry_uncertainty_content_sha256 != _sha256_bytes(
                uncertainty_material
            ):
                raise ValueError("runtime draw authority binds foreign geometry bytes")

        if (
            type(self.shaft_work_fractions) is not tuple
            or len(self.shaft_work_fractions) != 6
            or any(type(item) is not QSCTrayShaftWorkFraction for item in self.shaft_work_fractions)
        ):
            raise TypeError("runtime authority requires six exact shaft-work fractions")
        if tuple(item.physical_tray_id for item in self.shaft_work_fractions) != (
            QSC_MECH01_SOLID_TRAY_ORDER
        ):
            raise ValueError("runtime shaft-work fractions differ from solid tray order")
        for values in (
            tuple(item.wet_holdup_proportional_fraction for item in self.shaft_work_fractions),
            tuple(item.bottom_heavy_fraction for item in self.shaft_work_fractions),
        ):
            if math.fsum(values) != 1.0:
                raise ValueError("runtime shaft-work fractions must close exactly to one")
        if self.shaft_work_fractions != _build_work_fractions():
            raise ValueError("runtime shaft-work fractions differ from D3 bracket definitions")

        # The PR-08 base schema predates the D0 typed class.  Every nested
        # field must therefore use its only compatible non-machine category;
        # the wrapper supplies the stricter typed D0 standing above.
        provenances = [
            self.common_shaft.machine.provenance,
            self.common_shaft.shaft_direction.provenance,
            self.common_shaft.accepted_shaft_history.provenance,
            self.common_shaft.speed_envelope.semantics_provenance,
            self.common_shaft.drive_train.identity_provenance,
            self.common_shaft.drive_train.convention_provenance,
        ]
        speed_fields = (
            self.common_shaft.speed_envelope.zero_rpm,
            self.common_shaft.speed_envelope.crawl_rpm,
            self.common_shaft.speed_envelope.minimum_stable_rpm,
            self.common_shaft.speed_envelope.nominal_rpm,
            self.common_shaft.speed_envelope.maximum_continuous_rpm,
            self.common_shaft.speed_envelope.trip_rpm,
            self.common_shaft.speed_envelope.restart_rpm,
            self.common_shaft.speed_envelope.ramp_up_limit_rpm_per_s,
            self.common_shaft.speed_envelope.ramp_down_limit_rpm_per_s,
        )
        drive_fields = (
            self.common_shaft.drive_train.gearbox_ratio,
            self.common_shaft.drive_train.rated_post_gearbox_power_w,
            self.common_shaft.drive_train.rated_post_gearbox_torque_n_m,
            self.common_shaft.drive_train.efficiency,
        )
        provenances.extend(field.provenance for field in (*speed_fields, *drive_fields))
        for tray in self.common_shaft.trays:
            arm = tray.arm_geometry
            provenances.extend(
                (
                    arm.identity_provenance,
                    arm.arm_count.provenance,
                    arm.inner_radius_m.provenance,
                    arm.outer_radius_m.provenance,
                    arm.blade_height_m.provenance,
                    arm.blade_width_m.provenance,
                    arm.blade_pitch_rad.provenance,
                    arm.blade_rake_rad.provenance,
                    arm.tip_clearance_m.provenance,
                    arm.swept_volume_per_revolution_m3.provenance,
                    tray.draw_geometry.geometry_provenance,
                    tray.draw_geometry.conveyance_domain_provenance,
                    tray.draw_geometry.actuation_domain_provenance,
                )
            )
        if any(
            item.evidence_class is not mechanical.MechanicalEvidenceClass.PROVISIONAL_ENGINEERING
            or item.source_content_sha256 != _PINNED_RULING_DIGEST
            for item in provenances
        ):
            raise ValueError("runtime authority contains non-D0/D1-D6 provenance")

    @property
    def definition_digest(self) -> str:
        # Owner-approved B9 (2026-08-31): the full digest walk is memoized
        # VALUE-KEYED through the frozen dataclass's field-wise hash/eq —
        # a mutated authority hashes differently and re-walks, so the
        # anti-mutation property is preserved while the constant object's
        # repeated walks (measured 56+ per interval) collapse to one.
        return _runtime_authority_definition_digest(self)


@functools.lru_cache(maxsize=64)
def _runtime_authority_definition_digest(authority: QSCMechanicalRuntimeAuthority) -> str:
    return _digest(
        domain="qsc-mechanical-runtime-authority",
        parts=(
            authority.schema_id,
            authority.schema_revision,
            authority.configuration_relative_path,
            authority.configuration_canonical_lf_sha256,
            authority.parameter_class.value,
            authority.common_shaft.definition_digest,
            *(item.definition_digest for item in authority.k_authorities),
            *(binding.definition_digest for binding in authority.tray_bindings),
            *(item.definition_digest for item in authority.shaft_work_fractions),
            *authority.design_envelope_rpm,
            authority.low_alarm_rpm,
            authority.low_speed_trip_rpm,
            *authority.retired_values_rpm,
            *authority.retired_gui_range_rpm,
            *authority.shaft_power_bracket_w,
            authority.nominal_boundary_partition.value,
            authority.conservative_boundary_endpoint.value,
            authority.nominal_within_meal_distribution.value,
            authority.sensitivity_within_meal_distribution.value,
            *authority.delta_z_ex_bracket_m,
            *authority.phi_ex_bracket,
            authority.f_state,
            *authority.wet_bed_density_bracket_kg_m3,
        ),
    )


def _build_runtime_authority(*, shaft_rpm: float = 11.0) -> QSCMechanicalRuntimeAuthority:
    common = _build_common_shaft_authority()
    k_authorities = _build_k_authorities(common, shaft_rpm=shaft_rpm)
    tray_map = {tray.physical_tray_id: tray for tray in common.trays}
    bindings = tuple(
        QSCTrayMechanicalBinding(
            physical_tray_id=tray_id,
            common_shaft_definition_digest=common.definition_digest,
            accepted_history_definition_digest=(common.accepted_shaft_history.definition_digest),
            tray_mechanical_definition_digest=tray_map[tray_id].definition_digest,
            arm_geometry_definition_digest=tray_map[tray_id].arm_geometry.definition_digest,
            k_authority_definition_digest=k_authority.definition_digest,
            draw_geometry_definition_digest=tray_map[tray_id].draw_geometry.definition_digest,
        )
        for tray_id, k_authority in zip(QSC_MECH01_SOLID_TRAY_ORDER, k_authorities, strict=True)
    )
    return QSCMechanicalRuntimeAuthority(
        schema_id=QSC_MECH01_SCHEMA_ID,
        schema_revision=QSC_MECH01_SCHEMA_REVISION,
        configuration_relative_path=QSC_MECH01_CONFIG_RELATIVE_PATH,
        configuration_canonical_lf_sha256=QSC_MECH01_CONFIG_CANONICAL_LF_SHA256,
        parameter_class=QSCMechanicalParameterClass.DECLARED_ENGINEERING_ASSUMPTION,
        common_shaft=common,
        k_authorities=k_authorities,
        tray_bindings=bindings,
        shaft_work_fractions=_build_work_fractions(),
        design_envelope_rpm=(4.0, 15.0),
        low_alarm_rpm=8.0,
        low_speed_trip_rpm=4.0,
        retired_values_rpm=(3.0, 33.0, 60.0),
        retired_gui_range_rpm=(0.0, 10.0),
        shaft_power_bracket_w=(25_000.0, 50_000.0, 165_000.0),
        nominal_boundary_partition=ShaftWorkBoundaryPartition.ALL_TO_MEAL,
        conservative_boundary_endpoint=ShaftWorkBoundaryPartition.ALL_TO_WALL,
        nominal_within_meal_distribution=(WithinMealDistribution.WET_HOLDUP_PROPORTIONAL_BY_TRAY),
        sensitivity_within_meal_distribution=(
            WithinMealDistribution.BOTTOM_HEAVY_BY_SOLID_ORDER_RANK
        ),
        delta_z_ex_bracket_m=(0.09, 0.17, 0.33),
        phi_ex_bracket=(0.3, 0.5, 0.7),
        f_state=1.0,
        wet_bed_density_bracket_kg_m3=(480.0, 600.0, 640.0, 720.0),
        declared_constant_shaft_rpm=shaft_rpm,
    )


# The entire six-tray authority contains exactly one PR-08 common-shaft root.
# Six K authorities refer to it by digest rather than cloning it per tray.
QSC_MECH01_RUNTIME_AUTHORITY = _build_runtime_authority()
QSC_MECH01_COMMON_SHAFT_AUTHORITY = QSC_MECH01_RUNTIME_AUTHORITY.common_shaft


@functools.lru_cache(maxsize=8)
def make_engineering_runtime_authority(shaft_rpm: float = 11.0) -> QSCMechanicalRuntimeAuthority:
    """Config-built MECH-01A authority at a ruled constant operating speed.

    Ruled 2026-08-30 (GT_PS2_V3_EXTENSION_AND_RPM_AUTHORITY_RULING): constant
    shaft speeds other than the D1 nominal are admissible within the D1
    operating envelope; every other pin is unchanged and still enforced.
    """

    if type(shaft_rpm) is not float or not math.isfinite(shaft_rpm):
        raise TypeError("shaft speed must be a finite exact binary64")
    return _build_runtime_authority(shaft_rpm=shaft_rpm)


def validate_qsc_mechanical_runtime_authority(
    authority: QSCMechanicalRuntimeAuthority,
    *,
    preverified_definition_digest: str | None = None,
) -> str:
    """Accept only exact source-pinned, configuration-built MECH-01A graphs.

    The canonical D1 nominal graph must reproduce its pinned digest; a ruled
    constant-speed variant must reproduce, digest for digest, the graph this
    module itself builds from the configuration at its declared speed.

    Owner-approved B9 (2026-08-31): a caller that already ran the FULL
    validation this interval may thread the returned digest back as
    ``preverified_definition_digest``.  The pin branch re-derives the
    authority's digest through the value-keyed memo — a mutated authority
    re-hashes, re-walks, and mismatches — so every field is still proven
    byte-equal to the fully-validated graph; only the redundant structural
    re-walk collapses.  Default ``None`` keeps the sealed full behavior.
    """

    if type(authority) is not QSCMechanicalRuntimeAuthority:
        raise TypeError("authority must be an exact QSCMechanicalRuntimeAuthority")
    if preverified_definition_digest is not None:
        if authority.definition_digest != preverified_definition_digest:
            raise ValueError("runtime authority differs from its interval-preverified digest")
        return preverified_definition_digest
    QSCMechanicalRuntimeAuthority.__post_init__(authority)
    common_digest = authority.common_shaft.definition_digest
    if common_digest != QSC_MECH01_EXPECTED_COMMON_SHAFT_DEFINITION_DIGEST:
        raise ValueError("common shaft is not the exact canonical config-built MECH-01A authority")
    runtime_digest = authority.definition_digest
    if authority.declared_constant_shaft_rpm == 11.0:
        if runtime_digest != QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST:
            raise ValueError("runtime is not the exact canonical config-built MECH-01A authority")
        return runtime_digest
    expected = make_engineering_runtime_authority(authority.declared_constant_shaft_rpm)
    if runtime_digest != expected.definition_digest:
        raise ValueError(
            "runtime is not the exact config-built MECH-01A authority at its declared speed"
        )
    return runtime_digest


__all__ = (
    "QSC_MECH01_COMMON_SHAFT_AUTHORITY",
    "QSC_MECH01_CONFIG_CANONICAL_LF_SHA256",
    "QSC_MECH01_CONFIG_PATH",
    "QSC_MECH01_CONFIG_RELATIVE_PATH",
    "QSC_MECH01_EXPECTED_COMMON_SHAFT_DEFINITION_DIGEST",
    "QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST",
    "QSC_MECH01_OPERATING_RPM_ENVELOPE",
    "QSC_MECH01_RUNTIME_AUTHORITY",
    "QSC_MECH01_SCHEMA_ID",
    "QSC_MECH01_SCHEMA_REVISION",
    "QSC_MECH01_SOLID_TRAY_ORDER",
    "QSCMechanicalParameterClass",
    "QSCMechanicalRuntimeAuthority",
    "QSCTrayMechanicalBinding",
    "QSCTrayShaftWorkFraction",
    "ShaftWorkBoundaryPartition",
    "WithinMealDistribution",
    "make_engineering_runtime_authority",
    "qsc_mech01_configuration_errors",
    "validate_qsc_mechanical_runtime_authority",
)
