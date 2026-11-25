"""The study lane's declared LEVERS, each independently activatable.

The study lane is not one relaxation.  R2 of the charter (Chebyshev surrogates
of the pinned property authorities) was the first and, until now, the only one;
owner ruling **W3-2** adds a second, ``gs_pressure_sweep`` — Gauss-Seidel sweep
ordering for the tower pressure Picard (wave-3 audit finding 3.2).

Each lever is named, selected independently, and stamped into the run
declaration, so a certificate always records WHICH relaxations produced it.  A
declaration that says ``levers: ["surrogate_props"]`` is a different claim from
one that says ``["surrogate_props", "gs_pressure_sweep"]``, and the difference
is visible without reading code.

Selection is a fifth opt-in token, and it is as strict as the four in
:mod:`activation`:

* absent, the selection is :data:`DEFAULT_LEVERS` — surrogates only, the
  historical study lane — so every record written before this module keeps
  meaning exactly what it said;
* present, it must name only known levers, each at most once, and at least
  one.  An unknown name, a duplicate, or an empty list RAISES.  There is no
  "run everything" spelling and no silent "run nothing": a study-lane opt-in
  that installs no lever would produce a run labelled ``study`` that is in fact
  the certified path, which is the one mislabelling this package exists to
  prevent.

Nothing in this module is imported by the certified lane.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Protocol

__all__ = [
    "AVAILABLE_LEVERS",
    "DEFAULT_LEVERS",
    "GS_PRESSURE_SWEEP_LEVER",
    "LEVERS_ARGV_PREFIX",
    "LEVERS_ENVIRONMENT_VARIABLE",
    "SURROGATE_PROPS_LEVER",
    "LeverStack",
    "build_lever_stack",
    "parse_lever_list",
    "resolve_levers",
]

#: R2 of the charter: Chebyshev surrogates over the pinned property solves.
SURROGATE_PROPS_LEVER = "surrogate_props"

#: Owner ruling W3-2: Gauss-Seidel sweep ordering for the tower pressure
#: Picard (wave-3 solver-struggle audit finding 3.2).
GS_PRESSURE_SWEEP_LEVER = "gs_pressure_sweep"

#: Canonical order.  A selection is always reported in this order, whatever
#: order it was typed in, so two declarations naming the same levers compare
#: equal as lists.
AVAILABLE_LEVERS = (SURROGATE_PROPS_LEVER, GS_PRESSURE_SWEEP_LEVER)

#: What an opt-in that does not mention levers means.  Deliberately NOT "all":
#: adding a lever to this package must never change what an existing launch
#: line does.
DEFAULT_LEVERS = (SURROGATE_PROPS_LEVER,)

LEVERS_ENVIRONMENT_VARIABLE = "DTDC_STUDY_LANE_LEVERS"
LEVERS_ARGV_PREFIX = "--study-lane-levers="


class Lever(Protocol):
    """What every lever must provide.  Install/restore must be reversible."""

    name: str

    def install(self) -> None: ...

    def restore(self) -> None: ...

    def installed_targets(self) -> list[str]: ...

    def usage(self) -> dict[str, Any]: ...


def parse_lever_list(text: str) -> tuple[str, ...]:
    """Parse a comma-separated lever selection.  Raises ``ValueError`` on junk."""

    tokens = [token.strip() for token in str(text).split(",")]
    if any(token == "" for token in tokens):
        raise ValueError(
            f"study-lane lever selection {text!r} has an empty entry; spell it as a "
            f"comma-separated list of {list(AVAILABLE_LEVERS)}"
        )
    unknown = [token for token in tokens if token not in AVAILABLE_LEVERS]
    if unknown:
        raise ValueError(
            f"unknown study-lane lever(s) {unknown}; available levers are "
            f"{list(AVAILABLE_LEVERS)}"
        )
    duplicates = sorted({token for token in tokens if tokens.count(token) > 1})
    if duplicates:
        raise ValueError(f"study-lane lever(s) named more than once: {duplicates}")
    return tuple(name for name in AVAILABLE_LEVERS if name in tokens)


def resolve_levers(
    environ: dict[str, str], argv: Sequence[str]
) -> tuple[tuple[str, ...], str]:
    """The selected levers and where the selection came from.

    argv beats the environment, exactly as the acknowledgment does.  Returns
    ``(DEFAULT_LEVERS, "default")`` when nothing named a selection.
    """

    for token in argv:
        text = str(token)
        if text.startswith(LEVERS_ARGV_PREFIX):
            return parse_lever_list(text[len(LEVERS_ARGV_PREFIX) :]), "argv"
    if LEVERS_ENVIRONMENT_VARIABLE in environ:
        return (
            parse_lever_list(environ[LEVERS_ENVIRONMENT_VARIABLE]),
            LEVERS_ENVIRONMENT_VARIABLE,
        )
    return DEFAULT_LEVERS, "default"


class LeverStack:
    """An ordered, reversible set of installed levers."""

    def __init__(self, levers: Sequence[Lever]) -> None:
        self._levers: tuple[Lever, ...] = tuple(levers)
        self._installed = False

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(lever.name for lever in self._levers)

    def lever(self, name: str) -> Any | None:
        for candidate in self._levers:
            if candidate.name == name:
                return candidate
        return None

    def install(self) -> None:
        if self._installed:
            return
        done: list[Lever] = []
        try:
            for lever in self._levers:
                lever.install()
                done.append(lever)
        except BaseException:
            # A half-installed stack is exactly the silent-downgrade failure
            # this package refuses, so unwind what did install before the
            # refusal escapes.
            for lever in reversed(done):
                lever.restore()
            raise
        self._installed = True

    def restore(self) -> None:
        if not self._installed:
            return
        for lever in reversed(self._levers):
            lever.restore()
        self._installed = False

    def installed_targets(self) -> list[str]:
        targets: list[str] = []
        for lever in self._levers:
            targets.extend(lever.installed_targets())
        return targets

    def usage(self) -> dict[str, Any]:
        return {lever.name: lever.usage() for lever in self._levers}

    def declarations(self) -> dict[str, Any]:
        record: dict[str, Any] = {}
        for lever in self._levers:
            declare = getattr(lever, "declaration", None)
            record[lever.name] = {} if declare is None else declare()
        return record


class _SurrogatePropsLever:
    """Adapter giving the R2 property facade the lever protocol."""

    name = SURROGATE_PROPS_LEVER

    def __init__(self, facade: Any) -> None:
        self.facade = facade

    def install(self) -> None:
        self.facade.install()

    def restore(self) -> None:
        self.facade.restore()

    def installed_targets(self) -> list[str]:
        return list(self.facade.installed_targets())

    def usage(self) -> dict[str, Any]:
        return dict(self.facade.usage())

    def declaration(self) -> dict[str, Any]:
        return {
            "record": "GT-PS-2-S2-STUDY-SURROGATE-LANE-2026-09-02-01",
            "relaxation": "Chebyshev surrogates over the pinned property solves",
        }


class _NullFacade:
    """Stand-in for the property facade when the surrogate lever is off.

    The activation record and the runner both report a surrogate census; with
    the lever off the honest census is an empty one, not a missing attribute.
    """

    def install(self) -> None:  # pragma: no cover - never installed
        raise RuntimeError("the null facade is not installable")

    def restore(self) -> None:
        return None

    def installed_targets(self) -> list[str]:
        return []

    def usage(self) -> dict[str, Any]:
        return {}


def build_lever_stack(names: Sequence[str], *, pack: Any) -> tuple[LeverStack, Any]:
    """Construct the selected levers.  Returns ``(stack, property_facade)``.

    ``property_facade`` is the R2 facade when that lever is selected and a
    :class:`_NullFacade` otherwise, so callers that report a surrogate census
    keep working with the lever off.
    """

    selected = tuple(names)
    unknown = [name for name in selected if name not in AVAILABLE_LEVERS]
    if unknown:
        raise ValueError(f"unknown study-lane lever(s) {unknown}")
    if not selected:
        raise ValueError("a study-lane opt-in must select at least one lever")

    built: list[Lever] = []
    facade: Any = _NullFacade()
    for name in AVAILABLE_LEVERS:
        if name not in selected:
            continue
        if name == SURROGATE_PROPS_LEVER:
            from dtdc_simulator.core2.study_lane import surrogate_props

            facade = surrogate_props.SurrogateProps(pack)
            built.append(_SurrogatePropsLever(facade))
        elif name == GS_PRESSURE_SWEEP_LEVER:
            from dtdc_simulator.core2.study_lane import gs_pressure_sweep

            built.append(gs_pressure_sweep.GaussSeidelPressureSweep())
    return LeverStack(built), facade
