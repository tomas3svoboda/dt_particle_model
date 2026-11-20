"""The DECLARED simulation-study lane (GT-PS-2-STUDY-LANE-CHARTER-2026-09-02-01).

The study lane is a variant of the tower that trades bit-identity for speed
under a ratified TRAJECTORY CERTIFICATE: on the declared benchmark runs it must
reproduce the certified lane within per-tray temperature +/-0.1 K, component
inventories +/-0.5% and event times +/-1%.  It is not the certified lane, it
produces no certified evidence, and it never runs unless a human explicitly
opts in and types an acknowledgment string.

The lane's relaxations are independently activatable LEVERS:

``surrogate_props``
    R2 of the charter: Chebyshev surrogates of the pinned property
    authorities' SOLVED outputs.
``gs_pressure_sweep``
    Owner ruling W3-2: Gauss-Seidel (cascade-ordered) sweep of the tower
    face-pressure Picard, on wave-3 audit finding 3.2.

Entry points
------------
``activate_if_opted_in``
    The gate.  Returns ``None`` for an ordinary certified run; raises on a
    half-typed opt-in; installs the selected levers and returns the activation
    (carrying the self-declaration for the run record) on a complete one.
``SurrogateProps``
    The property facade: mirrors the pinned signatures, refuses to evaluate
    outside its declared domain, and falls back to genuinely pinned physics.
``GaussSeidelPressureSweep``
    The W3-2 lever: same fixed point, same acceptance test, cascade sweep
    order; refuses to install over a production Picard it was not derived
    against.
``load_pack``
    Reads a surrogate pack and asserts its props era-binding.

Nothing here is imported by the certified lane.
"""

from __future__ import annotations

from dtdc_simulator.core2.study_lane.activation import (
    ACK_ARGV_PREFIX,
    ACK_ENVIRONMENT_VARIABLE,
    ENABLE_ENVIRONMENT_VARIABLE,
    LEVERS_ARGV_PREFIX,
    LEVERS_ENVIRONMENT_VARIABLE,
    REQUIRED_ACKNOWLEDGMENT,
    STUDY_LANE_ARGV_PREFIXES,
    StudyLaneActivation,
    StudyLaneOptInRefused,
    activate_if_opted_in,
    opt_in_requested,
)
from dtdc_simulator.core2.study_lane.gs_pressure_sweep import (
    GaussSeidelPressureSweep,
    PicardEraMismatch,
)
from dtdc_simulator.core2.study_lane.levers import (
    AVAILABLE_LEVERS,
    DEFAULT_LEVERS,
    GS_PRESSURE_SWEEP_LEVER,
    SURROGATE_PROPS_LEVER,
    LeverStack,
    build_lever_stack,
    parse_lever_list,
    resolve_levers,
)
from dtdc_simulator.core2.study_lane.pack import (
    PackEraMismatch,
    SurrogatePack,
    default_pack_path,
    load_pack,
)
from dtdc_simulator.core2.study_lane.roster import SURROGATE_SPECS, SurrogateSpec
from dtdc_simulator.core2.study_lane.surrogate_props import SurrogateProps

__all__ = [
    "ACK_ARGV_PREFIX",
    "ACK_ENVIRONMENT_VARIABLE",
    "AVAILABLE_LEVERS",
    "DEFAULT_LEVERS",
    "ENABLE_ENVIRONMENT_VARIABLE",
    "GS_PRESSURE_SWEEP_LEVER",
    "LEVERS_ARGV_PREFIX",
    "LEVERS_ENVIRONMENT_VARIABLE",
    "REQUIRED_ACKNOWLEDGMENT",
    "STUDY_LANE_ARGV_PREFIXES",
    "SURROGATE_PROPS_LEVER",
    "SURROGATE_SPECS",
    "GaussSeidelPressureSweep",
    "LeverStack",
    "PackEraMismatch",
    "PicardEraMismatch",
    "StudyLaneActivation",
    "StudyLaneOptInRefused",
    "SurrogatePack",
    "SurrogateProps",
    "SurrogateSpec",
    "activate_if_opted_in",
    "build_lever_stack",
    "default_pack_path",
    "load_pack",
    "opt_in_requested",
    "parse_lever_list",
    "resolve_levers",
]
