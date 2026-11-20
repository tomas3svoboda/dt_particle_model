"""The ``SurrogateProps`` facade: mirrors the pinned props, never extrapolates.

Installation replaces named attributes of the pinned property modules with
wrappers.  Nothing on disk is touched; the pinned sources stay byte-identical
and the certified lane, which never activates, never sees a wrapper.

Three properties are load-bearing and each is enforced structurally rather
than by convention:

**Signature mirroring.**  Every wrapper takes and returns exactly what the
pinned entry point takes and returns — including the frozen ``Saturation``
dataclass, constructed from the pinned module's own class so an ``isinstance``
check downstream cannot tell the lanes apart.

**No extrapolation, ever.**  A wrapper evaluates its series only when every
argument lies inside the DECLARED domain of the adopted surrogate.  A
non-finite argument, an argument past a domain edge, or a phase with no
adopted series all route to the pinned call.  There is no soft edge and no
clamping: the domain test is a closed-interval membership test on the values
the pack declares.

**The fallback is genuinely pinned.**  A naive wrapper would leak: the pinned
``saturation`` calls the module-global ``density``, which is itself patched, so
an out-of-domain "fallback" would still be answered partly by surrogates.  The
facade therefore raises a bypass depth for the whole duration of a fallback
call, and every wrapper short-circuits to its original while that depth is
non-zero.  An out-of-domain input is answered by the pinned authority and
nothing else, all the way down.

The bypass counter is a plain module global.  The study lane runs under the
project's pinned single-worker environment (``OMP_NUM_THREADS=1`` and friends,
asserted by the march script itself), so there is exactly one interpreter
thread walking the property stack; a thread-local would buy nothing and cost a
lookup in the hottest inner loop in the program.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from dtdc_simulator.core2.props import binary_gas, hexane, water
from dtdc_simulator.core2.study_lane import roster
from dtdc_simulator.core2.study_lane.chebyshev import eval_1d, eval_2d
from dtdc_simulator.core2.study_lane.pack import SurrogatePack

__all__ = ["SurrogateProps", "bypass_depth"]

_MODULES = {"water": water, "hexane": hexane}

_PINNED_CACHE_CLEARERS: tuple[tuple[Any, str], ...] = (
    (water, "_residual_cached"),
    (water, "_density_cached"),
    (water, "_saturation_cached"),
    (water, "_saturation_temperature_cached"),
    (water, "_state_Tp_cached"),
    (hexane, "_residual_cached"),
    (hexane, "_state_Tp_cached"),
    (binary_gas, "_state_cached"),
)

_INF = float("inf")

_bypass_depth = 0


def bypass_depth() -> int:
    """Current pinned-fallback nesting depth (0 == surrogates live)."""

    return _bypass_depth


class _Series1D:
    """One adopted one-dimensional Chebyshev output."""

    __slots__ = ("coefficients", "hi", "lo")

    def __init__(self, coefficients: tuple[float, ...], lo: float, hi: float) -> None:
        self.coefficients = coefficients
        self.lo = lo
        self.hi = hi

    def __call__(self, x: float) -> float:
        return eval_1d(self.coefficients, self.lo, self.hi, x)


class _Series2D:
    """One adopted two-dimensional Chebyshev output, plus its inverse transform.

    ``gas_constant`` is non-zero only for a ``compressibility`` series, where
    the fitted quantity is Z and the density is recovered as p / (Z R T).
    """

    __slots__ = ("coefficients", "gas_constant", "x_hi", "x_lo", "y_hi", "y_lo")

    def __init__(
        self,
        coefficients: tuple[tuple[float, ...], ...],
        x_lo: float,
        x_hi: float,
        y_lo: float,
        y_hi: float,
        gas_constant: float = 0.0,
    ) -> None:
        self.coefficients = coefficients
        self.x_lo = x_lo
        self.x_hi = x_hi
        self.y_lo = y_lo
        self.y_hi = y_hi
        self.gas_constant = gas_constant

    def __call__(self, x: float, y: float) -> float:
        value = eval_2d(self.coefficients, self.x_lo, self.x_hi, self.y_lo, self.y_hi, x, y)
        if self.gas_constant == 0.0:
            return value
        return y / (value * self.gas_constant * x)


class SurrogateProps:
    """Install/restore the adopted surrogates of one pack over the pinned props."""

    def __init__(self, pack: SurrogatePack) -> None:
        self.pack = pack
        self._installed = False
        self._originals: dict[tuple[str, str], Callable[..., Any]] = {}
        self._targets: list[str] = []
        self.hits: dict[str, int] = {}
        self.fallbacks: dict[str, int] = {}
        self._adopted = pack.adopted()

    # -- bookkeeping --------------------------------------------------------
    def _count_hit(self, name: str) -> None:
        self.hits[name] = self.hits.get(name, 0) + 1

    def _count_fallback(self, name: str) -> None:
        self.fallbacks[name] = self.fallbacks.get(name, 0) + 1

    def installed_targets(self) -> list[str]:
        return list(self._targets)

    def usage(self) -> dict[str, Any]:
        """Run-time census: surrogate hits and pinned fallbacks per surrogate."""

        names = sorted(set(self.hits) | set(self.fallbacks))
        return {
            name: {
                "surrogate_hits": self.hits.get(name, 0),
                "pinned_fallbacks": self.fallbacks.get(name, 0),
            }
            for name in names
        }

    # -- series access ------------------------------------------------------
    def _series_1d(self, spec: roster.SurrogateSpec) -> dict[str, _Series1D]:
        entry = self._adopted[spec.name]
        axis = entry["domain"][0]
        return {
            output: _Series1D(
                tuple(float(value) for value in payload["coefficients"]),
                float(axis["lo"]),
                float(axis["hi"]),
            )
            for output, payload in entry["outputs"].items()
        }

    def _series_2d(self, spec: roster.SurrogateSpec) -> _Series2D:
        entry = self._adopted[spec.name]
        x_axis, y_axis = entry["domain"]
        payload = entry["outputs"]["rho"]
        coefficients = tuple(
            tuple(float(value) for value in row) for row in payload["coefficients"]
        )
        gas_constant = (
            float(_MODULES[spec.fluid].R)
            if entry.get("transform") == "compressibility"
            else 0.0
        )
        return _Series2D(
            coefficients,
            float(x_axis["lo"]),
            float(x_axis["hi"]),
            float(y_axis["lo"]),
            float(y_axis["hi"]),
            gas_constant,
        )

    # -- installation -------------------------------------------------------
    def install(self) -> None:
        """Patch every adopted surrogate over its pinned entry point."""

        if self._installed:
            return
        builders = {
            "saturation": self._build_saturation,
            "scalar_1d": self._build_scalar_1d,
            "virial_1d": self._build_virial_1d,
        }
        for spec in roster.SURROGATE_SPECS:
            if spec.name not in self._adopted or spec.kind == "density_2d":
                continue
            builders[spec.kind](spec)
        self._build_density()
        _clear_pinned_caches()
        self._installed = True

    def restore(self) -> None:
        """Put every pinned attribute back and drop the caches the lane filled."""

        if not self._installed:
            return
        for (fluid, attribute), original in self._originals.items():
            setattr(_MODULES[fluid], attribute, original)
        self._originals.clear()
        self._targets.clear()
        _clear_pinned_caches()
        self._installed = False

    def _patch(self, fluid: str, attribute: str, wrapper: Callable[..., Any]) -> Callable[..., Any]:
        module = _MODULES[fluid]
        original = getattr(module, attribute)
        self._originals[(fluid, attribute)] = original
        wrapper.__name__ = attribute
        wrapper.__qualname__ = f"study_lane.{fluid}.{attribute}"
        setattr(module, attribute, wrapper)
        self._targets.append(f"{fluid}.{attribute}")
        return original

    # -- wrapper factories --------------------------------------------------
    def _build_saturation(self, spec: roster.SurrogateSpec) -> None:
        series = self._series_1d(spec)
        p_series = series["p"]
        rho_l_series = series["rho_liquid"]
        rho_v_series = series["rho_vapor"]
        dh_series = series["dh_vap"]
        axis = spec.axes[0]
        lo, hi = axis.lo, axis.hi
        name = spec.name
        saturation_class = _MODULES[spec.fluid].Saturation
        original_box: list[Callable[..., Any]] = []

        def wrapper(T: float):
            if _bypass_depth or not (lo <= T <= hi):
                self._count_fallback(name)
                return _pinned(original_box[0], T)
            self._count_hit(name)
            return saturation_class(
                T=T,
                p=p_series(T),
                rho_liquid=rho_l_series(T),
                rho_vapor=rho_v_series(T),
                dh_vap=dh_series(T),
                converged=True,
            )

        original_box.append(self._patch(spec.fluid, spec.patch_targets[0], wrapper))

    def _build_scalar_1d(self, spec: roster.SurrogateSpec) -> None:
        series = self._series_1d(spec)[spec.outputs[0]]
        axis = spec.axes[0]
        lo, hi = axis.lo, axis.hi
        name = spec.name
        original_box: list[Callable[..., Any]] = []

        def wrapper(x: float) -> float:
            if _bypass_depth or not (lo <= x <= hi):
                self._count_fallback(name)
                return _pinned(original_box[0], x)
            self._count_hit(name)
            return series(x)

        original_box.append(self._patch(spec.fluid, spec.patch_targets[0], wrapper))

    def _build_virial_1d(self, spec: roster.SurrogateSpec) -> None:
        series = self._series_1d(spec)
        axis = spec.axes[0]
        lo, hi = axis.lo, axis.hi
        for target in spec.patch_targets:
            self._patch(spec.fluid, target, self._virial_wrapper(spec, target, series[target], lo, hi))

    def _virial_wrapper(
        self,
        spec: roster.SurrogateSpec,
        target: str,
        series: _Series1D,
        lo: float,
        hi: float,
    ) -> Callable[[float], float]:
        name = f"{spec.name}[{target}]"
        module = _MODULES[spec.fluid]
        original = getattr(module, target)

        def wrapper(T: float) -> float:
            if _bypass_depth or not (lo <= T <= hi):
                self._count_fallback(name)
                return _pinned(original, T)
            self._count_hit(name)
            return series(T)

        return wrapper

    def _build_density(self) -> None:
        """One wrapper per fluid, dispatching on the phase's adopted series."""

        for fluid in ("water", "hexane"):
            phase_series: dict[str, tuple[_Series2D, roster.SurrogateSpec]] = {}
            target = ""
            for spec in roster.SURROGATE_SPECS:
                if spec.kind != "density_2d" or spec.fluid != fluid:
                    continue
                target = spec.patch_targets[0]
                if spec.name in self._adopted:
                    phase_series[spec.phase] = (self._series_2d(spec), spec)
            if not phase_series:
                continue
            self._patch(fluid, target, self._density_wrapper(fluid, phase_series))

    def _density_wrapper(
        self, fluid: str, phase_series: dict[str, tuple[_Series2D, roster.SurrogateSpec]]
    ) -> Callable[[float, float, str], float]:
        module = _MODULES[fluid]
        target = "density" if fluid == "hexane" else "_density_uncached"
        original = getattr(module, target)
        bounds = {
            phase: (spec.axes[0].lo, spec.axes[0].hi, spec.axes[1].lo, spec.axes[1].hi)
            for phase, (_, spec) in phase_series.items()
        }
        names = {phase: spec.name for phase, (_, spec) in phase_series.items()}
        series = {phase: entry[0] for phase, entry in phase_series.items()}

        def wrapper(T: float, p: float, phase: str) -> float:
            box = bounds.get(phase)
            if _bypass_depth or box is None:
                self._count_fallback(f"{fluid}.density:{phase}")
                return _pinned(original, T, p, phase)
            t_lo, t_hi, p_lo, p_hi = box
            if not (t_lo <= T <= t_hi and p_lo <= p <= p_hi):
                self._count_fallback(names[phase])
                return _pinned(original, T, p, phase)
            root = series[phase](T, p)
            # Structural sanity, not a tolerance: the pinned Helmholtz algebra
            # takes the logarithm of the reduced density, so a root that is not
            # strictly positive and finite is unusable.  The verified series
            # never produces one inside its declared box; if one ever appears
            # the call is answered by the pinned solve rather than by a value
            # the downstream algebra cannot consume.
            if not (root > 0.0) or root != root or root in (_INF, -_INF):
                self._count_fallback(names[phase])
                return _pinned(original, T, p, phase)
            self._count_hit(names[phase])
            return root

        return wrapper


def _pinned(original: Callable[..., Any], *args: Any) -> Any:
    """Run ``original`` with every surrogate wrapper short-circuited."""

    global _bypass_depth
    _bypass_depth += 1
    try:
        return original(*args)
    finally:
        _bypass_depth -= 1


def _clear_pinned_caches() -> None:
    """Drop every pinned exact-bits memo so no value crosses a lane boundary."""

    for module, attribute in _PINNED_CACHE_CLEARERS:
        cached = getattr(module, attribute, None)
        clear = getattr(cached, "cache_clear", None)
        if clear is not None:
            clear()
