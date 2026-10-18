"""Ruled retained-water continuation above the Luikov cap on the Paper-1 particle path.

Owner ruling GT_PS2_OWNER_RULING_PAPER1_FILM_IMBIBITION_BRIDGE_2026-10-17:
(1) the C8 retained-water continuation above W_cap reaches the particle cells'
equilibrium (coupled_pore), the surface node under a film (through the same
equilibrium), the dry-limit branch and the liquid fast path (sp.water_activity
and sp.water_retained no longer refuse above the cap); (2) W_sat, the loading at
a_w = 1, is a declared parameter, nominal 0.55 kg/kg dry, sweep [0.25, 0.97];
(3) the shape is Horrobin's finite-intercept form W_sat - W = c1 (1 - a_w)
joined continuously at (W_cap, 0.799); (4) a film-free surface follows
a_w = P y_w phi / f_liquid(T_s), i.e. the gas activity at the surface
temperature, which the equilibrium below delivers because the retained state
above the cap follows the gas activity instead of a pinned 0.799.

Scope at the n-hexane front (declared reading, evidence 131): the ruling names
the dry cells, the surface node, the dry-limit branch and the fast path; it
does not extend the hexane-wet core, whose retained water stays on the frozen
cap graph (wet_retained_cap, loading <= W_cap with a nonnegative capacity dual).
The zero-volume interface traces of the cut and face-arrival charts
(cut_transport._solve_interface_composition_with_audit) therefore keep the
native retained-water representation: the wet trace sits at the cap with the
KKT dual log(a_w,dry / a_w,cap), and its generalized potential equals the
actual dry-gas water potential exactly as before. The gas side of the trace
(T, y, fugacities) is the same under both laws, and the trace carries no
storage, so no inventory is affected.

There is ONE implementation of the law: ``through_bed_sorbed_water_law.
bridged_water_activity`` / ``bridged_water_retained`` (shared with the tower's
``SorbedWaterArm.activity``). The pinned particle kernels are reached by the
ruled in-memory rebinding pattern (install/uninstall); no particle source byte
carries a second copy. Modes: ``ruled`` (default, W_sat declared), ``former-bridge``
(the C8 tier-1 linear bridge at x_bf 0.25, the same shape at W_sat 0.25, for
replay) and ``plain-cap`` (nothing installed: the former capacity active set).
"""

from contextlib import contextmanager
from dataclasses import replace
import math

from dtdc_simulator.core2 import through_bed_sorbed_water_law as law
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import dry_limit_moisture_branch as dlmb
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa

RULING = "GT_PS2_OWNER_RULING_PAPER1_FILM_IMBIBITION_BRIDGE_2026-10-17"
W_SAT_NOMINAL = 0.55
W_SAT_SWEEP = (0.25, 0.97)
FORMER_BRIDGE_W_SAT = law.X_BF_NOMINAL_KG_KG
MODES = ("ruled", "former-bridge", "plain-cap")


def declared_rows(mode, w_sat):
    """Typed DECLARED_ENGINEERING_ASSUMPTION rows for the record."""
    return [
        {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "authority": RULING,
            "name": "particle_retained_water_continuation_above_cap",
            "mode": mode,
            "saturated_loading_kg_kg_dry": w_sat,
            "sweep_kg_kg_dry": list(W_SAT_SWEEP),
            "shape": "Horrobin finite intercept W_sat - W = c1 (1 - a_w), continuous at (W_cap, 0.799)",
            "implementation": "through_bed_sorbed_water_law.bridged_water_activity / bridged_water_retained",
            "hexane_front_interface_traces": "native cap graph with the KKT dual (wet core outside the ruled routes)",
            "physically_qualifying": False,
        }
    ]


class _Installed:
    active = None
    native_depth = 0


def install(mode="ruled", w_sat=W_SAT_NOMINAL):
    """Rebind the particle kernels to the shared law; returns the uninstall closure."""
    if mode not in MODES:
        raise ValueError(f"unknown imbibition-bridge mode {mode!r}")
    if _Installed.active is not None:
        raise RuntimeError("the imbibition bridge is already installed in this process")
    if mode == "plain-cap":
        _Installed.active = ("plain-cap", None)
        return uninstall
    if mode == "former-bridge":
        w_sat = FORMER_BRIDGE_W_SAT
    if type(w_sat) is not float or not math.isfinite(w_sat):
        raise ValueError("W_sat must be a finite binary64")
    low, high = W_SAT_SWEEP
    if not low <= w_sat <= high:
        raise ValueError(f"W_sat {w_sat!r} is outside the ruled sweep [{low}, {high}]")
    originals = {
        "sp.water_activity": sp.water_activity,
        "sp.water_retained": sp.water_retained,
        "dlmb.dry_limit_water_activity": dlmb.dry_limit_water_activity,
        "dlmb.dry_limit_water_activity_derivative": dlmb.dry_limit_water_activity_derivative,
        "cp._evaluate_equilibrium_cached": cp._evaluate_equilibrium_cached,
        "cut._solve_interface_composition_with_audit": cut._solve_interface_composition_with_audit,
    }
    cap = law.BRIDGE_ANCHOR_LOADING_KG_KG
    native_cap_activity = originals["sp.water_activity"](cap, sp.LuikovParams())
    slope = (1.0 - law.BRIDGE_ANCHOR_ACTIVITY) / (w_sat - cap)

    def water_activity(W_w, p):
        if math.isfinite(W_w) and W_w > p.W_cap:
            if p.W_cap != cap:
                raise ValueError("the ruled continuation is anchored at the frozen W_cap")
            return law.bridged_water_activity(W_w, w_sat)
        return originals["sp.water_activity"](W_w, p)

    def water_retained(a_w, p):
        if math.isfinite(a_w) and a_w >= native_cap_activity and a_w <= 1.0:
            return law.bridged_water_retained(a_w, w_sat)
        return originals["sp.water_retained"](a_w, p)

    def dry_limit_water_activity(loading, luikov=sp.LuikovParams()):
        if math.isfinite(loading) and loading > luikov.W_cap:
            return law.bridged_water_activity(loading, w_sat)
        return originals["dlmb.dry_limit_water_activity"](loading, luikov)

    def dry_limit_water_activity_derivative(loading, luikov=sp.LuikovParams()):
        if math.isfinite(loading) and loading > luikov.W_cap:
            return 0.0 if loading >= w_sat else slope
        return originals["dlmb.dry_limit_water_activity_derivative"](loading, luikov)

    native_cached = originals["cp._evaluate_equilibrium_cached"]

    def evaluate(temperature_k, pressure_pa, y_hexane, params, gas_accessible_fraction, *keys):
        state = native_cached(
            temperature_k, pressure_pa, y_hexane, params, gas_accessible_fraction, *keys
        )
        if _Installed.native_depth:
            return state
        if state.retained_water_active_set is not cp.RetainedWaterActiveSet.RETAINED_CAP:
            return state
        # Above the cap the retained water follows the gas activity along the
        # shared law; every derived field is rebuilt from the native state's own
        # quantities (no second copy of the equilibrium arithmetic).
        retained = law.bridged_water_retained(state.water_activity, w_sat)
        fugacity = state.binary_gas.fugacity_water_pa
        concentration = params.rho_dm_p * retained
        pv_per_concentration = (
            state.retained_water_pv_sensitivity_j_m3 / state.retained_water_concentration_kg_m3
        )
        h_mass = state.retained_water_enthalpy_j_mol / wa.M
        pv = concentration * pv_per_concentration
        retained_energy = concentration * h_mass - params.retained_water_liquid_volume_fraction * pv
        total_water = concentration + state.pore_water_concentration_kg_m3
        return replace(
            state,
            retained_water_active_set=cp.RetainedWaterActiveSet.LUIKOV_SMOOTH,
            retained_water_loading=retained,
            retained_water_fugacity_pa=fugacity,
            retained_water_concentration_kg_m3=concentration,
            total_water_concentration_kg_m3=total_water,
            total_water_concentration_mol_m3=total_water / wa.M,
            retained_water_energy_density_j_m3=retained_energy,
            retained_water_pv_sensitivity_j_m3=pv,
            energy_density_j_m3=state.energy_density_j_m3
            - state.retained_water_energy_density_j_m3
            + retained_energy,
            retained_water_potential_isothermal=math.log(fugacity / cp.REFERENCE_FUGACITY_PA),
        )

    native_interface = originals["cut._solve_interface_composition_with_audit"]

    def interface_composition(*args, **kwargs):
        _Installed.native_depth += 1
        try:
            return native_interface(*args, **kwargs)
        finally:
            _Installed.native_depth -= 1

    evaluate.cache_clear = native_cached.cache_clear
    evaluate.cache_info = native_cached.cache_info
    sp.water_activity = water_activity
    sp.water_retained = water_retained
    dlmb.dry_limit_water_activity = dry_limit_water_activity
    dlmb.dry_limit_water_activity_derivative = dry_limit_water_activity_derivative
    cp._evaluate_equilibrium_cached = evaluate
    cut._solve_interface_composition_with_audit = interface_composition
    native_cached.cache_clear()
    _Installed.active = (mode, w_sat, originals)
    return uninstall


def uninstall():
    active = _Installed.active
    _Installed.active = None
    if active is None or active[0] == "plain-cap":
        return
    originals = active[2]
    sp.water_activity = originals["sp.water_activity"]
    sp.water_retained = originals["sp.water_retained"]
    dlmb.dry_limit_water_activity = originals["dlmb.dry_limit_water_activity"]
    dlmb.dry_limit_water_activity_derivative = originals["dlmb.dry_limit_water_activity_derivative"]
    cp._evaluate_equilibrium_cached = originals["cp._evaluate_equilibrium_cached"]
    cut._solve_interface_composition_with_audit = originals[
        "cut._solve_interface_composition_with_audit"
    ]
    cp._evaluate_equilibrium_cached.cache_clear()


@contextmanager
def installed(mode="ruled", w_sat=W_SAT_NOMINAL):
    install(mode, w_sat)
    try:
        yield
    finally:
        uninstall()


def active():
    return None if _Installed.active is None else _Installed.active[:2]
