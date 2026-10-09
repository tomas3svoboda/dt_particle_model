"""Scoped actual-path water monotonicity bounds from native EOS expressions.

This is numerical evidence, conditional on four-ULP libm padding. Native
property code is evaluated in isolated interval namespaces; no property or
transport law is modified. The native rectangle enclosure remains fallback.
"""

from contextlib import contextmanager
from dataclasses import replace
import inspect
import math
from functools import lru_cache
from types import SimpleNamespace

from dtdc_simulator.core2.particle import actual_composition_force_path as acfp
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wa


def padded(value, direction):
    for _ in range(4):
        value = math.nextafter(value, direction)
    return value


class Interval:
    def __init__(self, lower, upper=None):
        self.lo = float(lower)
        self.hi = float(lower if upper is None else upper)
        if not math.isfinite(self.lo) or not math.isfinite(self.hi) or self.lo > self.hi:
            raise ValueError("invalid finite interval")

    @staticmethod
    def lift(value):
        return value if isinstance(value, Interval) else Interval(value)

    @staticmethod
    def enclosed(values):
        return Interval(padded(min(values), -math.inf), padded(max(values), math.inf))

    def __add__(self, other):
        other = self.lift(other)
        return self.enclosed((self.lo + other.lo, self.hi + other.hi))

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.lift(other)

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        if other is self:
            return self**2
        return self.enclosed(
            (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.lift(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("interval division crosses zero")
        return self * self.enclosed((1 / other.hi, 1 / other.lo))

    def __rtruediv__(self, other):
        return self.lift(other) / self

    def __pow__(self, exponent):
        if self.lo <= 0 and exponent != int(exponent):
            raise ValueError("fractional power needs a positive interval")
        if exponent < 0 and self.lo <= 0 <= self.hi:
            raise ValueError("negative power crosses zero")
        values = [self.lo**exponent, self.hi**exponent]
        if exponent > 0 and exponent % 2 == 0 and self.lo <= 0 <= self.hi:
            values.append(0.0)
        return self.enclosed(values)

    def __lt__(self, other):
        return self.hi < self.lift(other).lo

    def __le__(self, other):
        return self.hi <= self.lift(other).lo

    def pair(self):
        return self.lo, self.hi


def interval_exp(value):
    value = Interval.lift(value)
    return Interval.enclosed((math.exp(value.lo), math.exp(value.hi)))


def interval_log(value):
    value = Interval.lift(value)
    if value.lo <= 0:
        raise ValueError("interval logarithm is not positive")
    return Interval.enclosed((math.log(value.lo), math.log(value.hi)))


@lru_cache(maxsize=3)
def interval_namespace(module):
    namespace = dict(vars(module))
    namespace["math"] = SimpleNamespace(
        exp=interval_exp,
        log=interval_log,
        sqrt=math.sqrt,
        isfinite=lambda value: isinstance(value, Interval) or math.isfinite(value),
    )
    names = (
        ("_phi0", "_residual_uncached")
        if module is wa
        else ("_residual_uncached",) if module is hx else ("cross_second_virial",)
    )
    for name in names:
        exec(
            compile(
                inspect.getsource(getattr(module, name)),
                f"<native interval {module.__name__}.{name}>",
                "exec",
            ),
            namespace,
        )
    return namespace


def pure_virial(module, temperature):
    tau = module.TC / temperature
    state = interval_namespace(module)["_residual_uncached"](1e-8, tau)
    density = module.RHOC / module.M if module is wa else module.RHOC
    return (
        state.phir_d / density,
        -tau / temperature * state.phir_dt / density,
        tau / temperature**2 * (2 * state.phir_dt + tau * state.phir_dtt) / density,
    )


def water_liquid_enthalpy_bound(temperature, pressure):
    mid = (temperature.lo + temperature.hi) / 2
    points = [wa.state_Tp(t, pressure, "liquid") for t in (temperature.lo, mid, temperature.hi)]
    density = Interval(
        min(p.rho_mass for p in points) - 1e-6, max(p.rho_mass for p in points) + 1e-6
    )
    evaluate = interval_namespace(wa)["_residual_uncached"]
    # Prove this density tube brackets the fixed-pressure liquid root for
    # every T: positive dp/drho and dp/dT make the extrema these two corners.
    res = evaluate(density / wa.RHOC, wa.TC / temperature)
    delta, tau = res.delta, res.tau
    stiffness = 1 + 2 * delta * res.phir_d + delta**2 * res.phir_dd
    thermal = 1 + delta * res.phir_d - delta * tau * res.phir_dt
    if stiffness.lo <= 0 or thermal.lo <= 0:
        raise ValueError("liquid density tube has unresolved monotonicity")

    def pressure_at(t, rho):
        point = evaluate(Interval(rho / wa.RHOC), Interval(wa.TC / t))
        return rho * wa.R * t * (1 + point.delta * point.phir_d)

    if (
        pressure_at(temperature.hi, density.lo).hi >= pressure
        or pressure_at(temperature.lo, density.hi).lo <= pressure
    ):
        raise ValueError("liquid density tube does not enclose the pressure root")
    cp = -wa.R * tau**2 * (res.phi0_tt + res.phir_tt) + wa.R * thermal**2 / stiffness
    if cp.lo <= 0:
        raise ValueError("liquid heat capacity is not certified positive")
    radius = max(abs(cp.lo), abs(cp.hi)) * wa.M * (temperature.hi - temperature.lo) / 2
    center = points[1].h_mass * wa.M
    # Protect the native scalar EOS root and caloric evaluation's final
    # arithmetic too. This is an energy interval, not an energy adjustment.
    radius += 128 * math.ulp(center)
    return Interval(center - radius, center + radius), cp.pair()


def water_path_derivative_bound(activities):
    t0, y0 = activities.left_temperature_k, activities.left_y_hexane
    dt, dy = activities.delta_temperature_k, activities.delta_y_hexane
    temperature = Interval(min(t0, t0 + dt), max(t0, t0 + dt))
    y = Interval(min(y0, y0 + dy), max(y0, y0 + dy))
    if (
        temperature.lo < 323.15
        or temperature.hi > 433
        or temperature.hi - temperature.lo > 0.5
        or y.hi - y.lo > 0.25
    ):
        raise ValueError("outside local certificate domain")
    bw, dbw, d2bw = pure_virial(wa, temperature)
    bh, dbh, d2bh = pure_virial(hx, temperature)
    cross = interval_namespace(bg)["cross_second_virial"](temperature, activities.pore.k_wh)
    bc, dbc = cross.B_wh, cross.dB_wh_dT
    yw = 1 - y
    mix = yw**2 * bw + 2 * yw * y * bc + y**2 * bh
    dmix = yw**2 * dbw + 2 * yw * y * dbc + y**2 * dbh
    aw = 2 * (yw * bw + y * bc) - mix
    daw = 2 * (yw * dbw + y * dbc) - dmix
    residual_h = activities.pressure_pa * (aw - temperature * daw)
    liquid_h, cp = water_liquid_enthalpy_bound(temperature, activities.pressure_pa)
    phi0 = interval_namespace(wa)["_phi0"](1.0, wa.TC / temperature)
    ideal_h = wa.R * wa.M * temperature * (1 + (wa.TC / temperature) * phi0[1])
    qt = -residual_h / (bg.R * temperature**2) - (ideal_h - liquid_h) / (
        wa.R * wa.M * temperature**2
    )
    virial = 2 * activities.pressure_pa * (bw - 2 * bc + bh) / (bg.R * temperature)
    qy = -1 / (1 - y) - virial * y
    derivative = dt * qt + dy * qy
    d2mix = yw**2 * d2bw + 2 * yw * y * cross.d2B_wh_dT2 + y**2 * d2bh
    d2aw = 2 * (yw * d2bw + y * cross.d2B_wh_dT2) - d2mix
    cp_residual = -activities.pressure_pa * temperature * d2aw
    cp_ideal = wa.R * wa.M * (1 - (wa.TC / temperature) ** 2 * phi0[2])
    cp_liquid = Interval(*cp) * wa.M
    qtt = (
        -cp_residual / (bg.R * temperature**2)
        + 2 * residual_h / (bg.R * temperature**3)
        - (cp_ideal - cp_liquid) / (wa.R * wa.M * temperature**2)
        + 2 * (ideal_h - liquid_h) / (wa.R * wa.M * temperature**3)
    )
    virial_dt = (
        2
        * activities.pressure_pa
        / bg.R
        * ((dbw - 2 * dbc + dbh) / temperature - (bw - 2 * bc + bh) / temperature**2)
    )
    second = dt**2 * qtt - 2 * dt * dy * virial_dt * y + dy**2 * (-1 / (1 - y) ** 2 - virial)
    return derivative.pair(), {
        "water_log_temperature_derivative": qt.pair(),
        "water_log_composition_derivative": qy.pair(),
        "liquid_cp_mass_bounds": cp,
        "path_derivative": derivative.pair(),
        "path_second_derivative": second.pair(),
    }


@contextmanager
def monotone_water_path_enclosure(journal):
    original = acfp._activity_box

    def box(activities, start, end, depth):
        raw = original(activities, start, end, depth)
        if (
            raw.water_activity_upper_bound < 1
            and raw.water_activity_lower_bound > activities.activity_at_ref
        ):
            return raw
        cached = getattr(activities, "_study_water_derivative_bound", None)
        if cached is None:
            try:
                cached = water_path_derivative_bound(activities)
            except ValueError, ArithmeticError:
                cached = (None, None)
            activities._study_water_derivative_bound = cached
        derivative, evidence = cached
        if derivative is None:
            return raw
        left = activities.point(start, role="native_interval_water_monotonicity")
        right = activities.point(end, role="native_interval_water_monotonicity")
        if derivative[0] > 0 or derivative[1] < 0:
            method = "whole-segment monotonicity"
            lower = acfp._outward_lower(min(left.water_activity, right.water_activity))
            upper = acfp._outward_upper(max(left.water_activity, right.water_activity))
        else:
            # Linear-interpolation remainder: f(s)-linear(f) lies between
            # -M*w^2/8 and -m*w^2/8 when m <= f'' <= M. Thus even a
            # nonmonotone log-activity can be bounded over the WHOLE segment.
            second_low, second_high = evidence["path_second_derivative"]
            width = end - start
            log_low = (
                min(math.log(left.water_activity), math.log(right.water_activity))
                - max(0.0, second_high) * width**2 / 8
            )
            log_high = (
                max(math.log(left.water_activity), math.log(right.water_activity))
                + max(0.0, -second_low) * width**2 / 8
            )
            lower = acfp._outward_lower(math.exp(log_low))
            upper = acfp._outward_upper(math.exp(log_high))
            method = "whole-segment log-activity interpolation remainder"
        if lower <= activities.activity_at_ref or upper >= 1:
            return raw
        journal.append(
            {
                "left_temperature_k": left.temperature_k,
                "right_temperature_k": right.temperature_k,
                "start_fraction": start,
                "end_fraction": end,
                "water_activity_bounds": (lower, upper),
                "certificate": evidence,
                "libm_padding_ulps": 4,
                "method": method,
            }
        )
        return replace(
            raw,
            water_activity_lower_bound=lower,
            water_activity_upper_bound=upper,
            water_path_derivative_bounds=derivative,
        )

    acfp._activity_box = box
    try:
        yield journal
    finally:
        acfp._activity_box = original


@contextmanager
def monotone_hexane_path_enclosure(journal):
    """Use the existing native derivative envelope to bound a monotone path.

    The independent T/Y corner box can exceed saturation even when the
    actual line segment is strictly admissible. A sign-definite native
    whole-path derivative proves its extremum is an endpoint. No activity,
    derivative authority or physical admission tolerance is modified.
    """
    original = acfp._activity_box

    def box(activities, start, end, depth):
        raw = original(activities, start, end, depth)
        derivative = raw.hexane_path_derivative_bounds
        if raw.hexane_activity_upper_bound < 1.0 or not (
            derivative[0] > 0.0 or derivative[1] < 0.0
        ):
            return raw
        left = activities.point(start, role="native_hexane_monotonicity")
        right = activities.point(end, role="native_hexane_monotonicity")
        lower = acfp._outward_lower(min(left.hexane_activity, right.hexane_activity))
        upper = acfp._outward_upper(max(left.hexane_activity, right.hexane_activity))
        if upper >= 1.0:
            return raw
        journal.append(
            {
                "method": "whole-segment monotonicity using unchanged native derivative envelope",
                "left_temperature_k": left.temperature_k,
                "right_temperature_k": right.temperature_k,
                "start_fraction": start,
                "end_fraction": end,
                "path_derivative_bounds": derivative,
                "hexane_activity_bounds": (lower, upper),
                "activity_clamped": False,
            }
        )
        return replace(raw, hexane_activity_lower_bound=lower, hexane_activity_upper_bound=upper)

    acfp._activity_box = box
    try:
        yield journal
    finally:
        acfp._activity_box = original
