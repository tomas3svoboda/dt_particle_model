"""Existing active-water Faner initializer, relocated from the recovery workspace.

The phase equations, numerical chart and controls are unchanged. Class/codec
identities of the pure-bulk carriers remain in the bulk_endpoint module.
"""

from dataclasses import asdict, replace
import math
from scipy.optimize import least_squares
import finite_film_core2_study as study

core = study.core


def phase_point(water, pressure, pore):
    """Native virial/liquid fugacity equations: a_h=1, a_w=Luikov(W)."""
    target_aw = core.sp.water_activity(water, pore.luikov)

    def equations(x):
        temperature, y = map(float, x)
        gas = core.bg.state(temperature, pressure, 1.0 - y, y, k_wh=pore.k_wh)
        water_liquid = core.wa.state_Tp(temperature, pressure, "liquid")
        hexane_liquid = core.hx.state_Tp(temperature, pressure, "liquid")
        return [
            math.log(gas.fugacity_hexane_pa / hexane_liquid.fugacity),
            math.log(gas.fugacity_water_pa / (target_aw * water_liquid.fugacity)),
        ]

    root = least_squares(
        equations,
        [339.0, 0.96],
        bounds=([330.0, 0.6], [342.0, 1.0 - 1e-12]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=100,
    )
    if not root.success or max(map(abs, equations(root.x))) > 1e-10:
        raise RuntimeError("active-water native phase initializer did not close")
    return tuple(map(float, root.x))


def active_binary_controls(original, config):
    """One numerical temperature chart for every declared moisture case.

    The upper endpoint is the native binary phase point at the existing
    minimum admitted retained-water loading. All moisture cases use this
    same chart. Phase equations and acceptance checks remain unchanged.
    """
    pressure = config.dry.pressure_pa

    minimum_water = original.wet_water_bounds[0]
    endpoint_temperature, endpoint_y = phase_point(minimum_water, pressure, config.dry.pore)
    upper = math.nextafter(endpoint_temperature, -math.inf)
    controls = replace(
        original,
        interface_temperature_bounds_k=(original.interface_temperature_bounds_k[0], upper),
        wet_temperature_bounds_k=(330.0, 345.0),
    )
    audit = {
        "reason": "shared active-water phase roots at low moisture exceed the old 340 K numerical chart",
        "old_interface_temperature_bounds_k": list(original.interface_temperature_bounds_k),
        "new_interface_temperature_bounds_k": list(controls.interface_temperature_bounds_k),
        "old_wet_temperature_bounds_k": list(original.wet_temperature_bounds_k),
        "new_wet_temperature_bounds_k": list(controls.wet_temperature_bounds_k),
        "pressure_pa": pressure,
        "physical_upper_endpoint": "native active-binary phase root at the unchanged minimum retained-water loading",
        "minimum_retained_water_loading": minimum_water,
        "endpoint_temperature_k": endpoint_temperature,
        "endpoint_y_hexane": endpoint_y,
        "old_controls": asdict(original),
        "new_controls": asdict(controls),
        "binary_phase_checks_changed": False,
        "physical_wet_caloric_bounds_changed": False,
        "acceptance_tolerances_changed": False,
        "work_budgets_changed": False,
        "wet_water_bounds_changed": False,
    }
    old_y_bounds = config.interface_composition.y_hexane_bounds
    new_y_bounds = (old_y_bounds[0], max(old_y_bounds[1], (1.0 + endpoint_y) / 2.0))
    config = replace(
        config,
        interface_composition=replace(config.interface_composition, y_hexane_bounds=new_y_bounds),
    )
    audit["old_interface_composition_root_bracket"] = list(old_y_bounds)
    audit["new_interface_composition_root_bracket"] = list(new_y_bounds)
    return controls, config, audit
