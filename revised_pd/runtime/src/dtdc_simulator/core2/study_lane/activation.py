"""The study-lane opt-in gate.

The study lane must be IMPOSSIBLE to enter by accident.  Four independent
conditions must all hold, and every one of them has to be typed by a human who
knows what the lane is:

1. ``DTDC_STUDY_LANE=1`` in the environment;
2. the acknowledgment string ``trajectory-certificate`` — supplied either as
   ``DTDC_STUDY_LANE_ACK`` or as the argv token
   ``--study-lane-ack=trajectory-certificate``;
3. a surrogate pack that exists and whose own digest checks out;
4. that pack's props era-binding matching the installed pinned authorities.

Any one missing and :func:`activate_if_opted_in` returns ``None``: the caller
runs the certified path, unchanged, with no surrogate installed and nothing
patched.  A partial opt-in (the flag without the acknowledgment, or the
acknowledgment misspelled) is a REFUSAL, not a silent downgrade — it raises,
so a typo can never be mistaken for a certified run.

WHICH relaxations run is a separate, optional selection — see
:mod:`dtdc_simulator.core2.study_lane.levers`.  Absent, it is the historical
study lane (``surrogate_props`` alone); present, it must name known levers and
is itself refused when malformed.  Conditions 3 and 4 hold whatever the
selection is: the pack is loaded and era-checked on every activation, because
it is also the record of which pinned-props era the run belongs to.

Activation prints a banner to stderr.  A study-lane run is never quiet.
"""

from __future__ import annotations

import os
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dtdc_simulator.core2.study_lane import levers as levers_module
from dtdc_simulator.core2.study_lane import pack as pack_module

__all__ = [
    "ACK_ARGV_PREFIX",
    "ACK_ENVIRONMENT_VARIABLE",
    "ENABLE_ENVIRONMENT_VARIABLE",
    "LEVERS_ARGV_PREFIX",
    "LEVERS_ENVIRONMENT_VARIABLE",
    "PACK_PATH_ENVIRONMENT_VARIABLE",
    "REQUIRED_ACKNOWLEDGMENT",
    "StudyLaneActivation",
    "StudyLaneOptInRefused",
    "activate_if_opted_in",
    "live_activation",
    "opt_in_requested",
]

ENABLE_ENVIRONMENT_VARIABLE = "DTDC_STUDY_LANE"
ACK_ENVIRONMENT_VARIABLE = "DTDC_STUDY_LANE_ACK"
PACK_PATH_ENVIRONMENT_VARIABLE = "DTDC_STUDY_LANE_PACK"
ACK_ARGV_PREFIX = "--study-lane-ack="
REQUIRED_ACKNOWLEDGMENT = "trajectory-certificate"

#: Re-exported so a caller needs one import to strip both study-lane argv
#: tokens before handing the rest to the march's positional knob parser.
LEVERS_ENVIRONMENT_VARIABLE = levers_module.LEVERS_ENVIRONMENT_VARIABLE
LEVERS_ARGV_PREFIX = levers_module.LEVERS_ARGV_PREFIX

#: Every argv token this package owns.  Anything starting with one of these is
#: ours and must not reach the march.
STUDY_LANE_ARGV_PREFIXES = (ACK_ARGV_PREFIX, LEVERS_ARGV_PREFIX)

#: The live activation of this process, or ``None``: set by
#: :func:`activate_if_opted_in` and cleared by ``deactivate``.  Read through
#: :func:`live_activation` by the exploration mode's launch gate (owner ruling of
#: 2026-10-13, ``docs/GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13.md``), which
#: refuses without it; environment variables alone do not prove the gate ran.
_LIVE_ACTIVATION = None


def live_activation():
    """The study-lane activation live in this process, or ``None``."""

    return _LIVE_ACTIVATION


class StudyLaneOptInRefused(RuntimeError):
    """A study-lane opt-in was attempted but is incomplete or invalid."""


@dataclass(frozen=True)
class StudyLaneActivation:
    """A live study-lane activation and its self-declaration."""

    pack: pack_module.SurrogatePack
    facade: Any
    acknowledgment_source: str
    levers: tuple[str, ...] = levers_module.DEFAULT_LEVERS
    lever_source: str = "default"
    lever_stack: Any = None

    def declaration(self) -> dict[str, Any]:
        record = dict(self.pack.declaration())
        record["acknowledgment"] = REQUIRED_ACKNOWLEDGMENT
        record["acknowledgment_source"] = self.acknowledgment_source
        record["enable_variable"] = ENABLE_ENVIRONMENT_VARIABLE
        record["levers"] = list(self.levers)
        record["lever_source"] = self.lever_source
        record["available_levers"] = list(levers_module.AVAILABLE_LEVERS)
        record["lever_declarations"] = (
            {} if self.lever_stack is None else self.lever_stack.declarations()
        )
        record["installed_targets"] = (
            self.facade.installed_targets()
            if self.lever_stack is None
            else self.lever_stack.installed_targets()
        )
        return record

    def lever_usage(self) -> dict[str, Any]:
        """Per-lever run-time census, keyed by lever name."""

        return {} if self.lever_stack is None else self.lever_stack.usage()

    def deactivate(self) -> None:
        global _LIVE_ACTIVATION
        if _LIVE_ACTIVATION is self:
            _LIVE_ACTIVATION = None
        if self.lever_stack is None:
            self.facade.restore()
            return
        self.lever_stack.restore()


def opt_in_requested(
    environ: dict[str, str] | None = None, argv: Sequence[str] | None = None
) -> bool:
    """True when ANYTHING in the launch asks for the study lane.

    Deliberately generous: it is true for a malformed opt-in too, so the caller
    can distinguish "certified run" from "botched study-lane run" and refuse
    the latter instead of quietly running certified.
    """

    env = os.environ if environ is None else environ
    args = sys.argv if argv is None else argv
    if env.get(ENABLE_ENVIRONMENT_VARIABLE):
        return True
    if env.get(ACK_ENVIRONMENT_VARIABLE):
        return True
    if env.get(LEVERS_ENVIRONMENT_VARIABLE):
        return True
    return any(
        str(token).startswith(STUDY_LANE_ARGV_PREFIXES) for token in args
    )


def _resolve_acknowledgment(
    environ: dict[str, str], argv: Sequence[str]
) -> tuple[str | None, str]:
    for token in argv:
        text = str(token)
        if text.startswith(ACK_ARGV_PREFIX):
            return text[len(ACK_ARGV_PREFIX) :], "argv"
    if ACK_ENVIRONMENT_VARIABLE in environ:
        return environ[ACK_ENVIRONMENT_VARIABLE], ACK_ENVIRONMENT_VARIABLE
    return None, "absent"


def activate_if_opted_in(
    environ: dict[str, str] | None = None,
    argv: Sequence[str] | None = None,
    *,
    banner: bool = True,
) -> StudyLaneActivation | None:
    """Install the surrogate facade iff the full opt-in is present.

    Returns ``None`` when nothing asked for the study lane.  Raises
    :class:`StudyLaneOptInRefused` when something asked but the request is
    incomplete, misspelled, or bound to the wrong props era.
    """

    env = dict(os.environ) if environ is None else dict(environ)
    args = list(sys.argv) if argv is None else list(argv)

    if not opt_in_requested(env, args):
        return None

    enabled = env.get(ENABLE_ENVIRONMENT_VARIABLE, "")
    if enabled != "1":
        raise StudyLaneOptInRefused(
            f"study lane requested but {ENABLE_ENVIRONMENT_VARIABLE}={enabled!r}; "
            f"set {ENABLE_ENVIRONMENT_VARIABLE}=1 to enable it"
        )

    acknowledgment, source = _resolve_acknowledgment(env, args)
    if acknowledgment != REQUIRED_ACKNOWLEDGMENT:
        raise StudyLaneOptInRefused(
            "study lane requested without the required acknowledgment; supply "
            f"{ACK_ENVIRONMENT_VARIABLE}={REQUIRED_ACKNOWLEDGMENT} or the argv "
            f"token {ACK_ARGV_PREFIX}{REQUIRED_ACKNOWLEDGMENT} "
            f"(got {acknowledgment!r} from {source})"
        )

    try:
        selected, lever_source = levers_module.resolve_levers(env, args)
    except ValueError as error:
        raise StudyLaneOptInRefused(str(error)) from error

    override = env.get(PACK_PATH_ENVIRONMENT_VARIABLE)
    loaded = pack_module.load_pack(Path(override) if override else None)

    try:
        stack, facade = levers_module.build_lever_stack(selected, pack=loaded)
    except ValueError as error:
        raise StudyLaneOptInRefused(str(error)) from error
    stack.install()
    activation = StudyLaneActivation(
        pack=loaded,
        facade=facade,
        acknowledgment_source=source,
        levers=tuple(selected),
        lever_source=lever_source,
        lever_stack=stack,
    )
    global _LIVE_ACTIVATION
    _LIVE_ACTIVATION = activation

    if banner:
        adopted = sorted(loaded.adopted())
        surrogates_on = levers_module.SURROGATE_PROPS_LEVER in selected
        print(
            "=" * 78
            + "\nDTDC STUDY LANE ACTIVE - this run is NOT certified evidence.\n"
            + f"  levers         : {', '.join(selected)} (from {lever_source})\n"
            + f"  surrogate pack : {loaded.record}\n"
            + f"  pack digest    : {loaded.pack_digest}\n"
            + f"  certificate    : {loaded.payload['certificate_reference']}\n"
            + (
                f"  adopted        : {len(adopted)} surrogate(s): {', '.join(adopted)}\n"
                if surrogates_on
                else "  adopted        : none installed (surrogate_props lever off)\n"
            )
            + f"  patched        : {len(stack.installed_targets())} target(s)\n"
            + "=" * 78,
            file=sys.stderr,
            flush=True,
        )
    return activation
