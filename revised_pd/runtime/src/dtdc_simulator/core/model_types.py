"""Stable runtime-model contracts shared by the engine and its adapters."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from dtdc_simulator.core import dc, dt_solver


class StageRole(str, Enum):
    PREDESOLV = "PREDESOLV"
    MAIN = "MAIN"
    SPARGE = "SPARGE"
    DRYER = "DRYER"
    COOLER = "COOLER"


class VaporPath(str, Enum):
    BYPASS = "BYPASS"
    THROUGH_BED = "THROUGH_BED"
    STEAM_SOURCE = "STEAM_SOURCE"


class SolidTransferDeviceType(str, Enum):
    PASSIVE_SWEPT_PORT = "PASSIVE_SWEPT_PORT"
    CONTROLLED_GATE = "CONTROLLED_GATE"
    ROTARY_AIRLOCK = "ROTARY_AIRLOCK"


DT_ROLES = {StageRole.PREDESOLV, StageRole.MAIN, StageRole.SPARGE}
DC_ROLES = {StageRole.DRYER, StageRole.COOLER}


@dataclass(frozen=True)
class StageSpec:
    id: str
    role: StageRole
    diameter_m: float
    bed_height_m: float
    # Per-tray sweep-arm turnover/residence multiplier (see schema.StageGeometry): the single
    # central shaft fixes ONE rpm for every tray; this captures the per-tray mixing difference
    # from arm geometry instead. Defaults 1.0 (arm identical to the base tray).
    arm_mixing_factor: float = 1.0
    vapor_path: VaporPath = VaporPath.THROUGH_BED

    @property
    def volume_m3(self) -> float:
        return math.pi / 4.0 * self.diameter_m**2 * self.bed_height_m


@dataclass(frozen=True)
class SolidTransferSpec:
    """Physical meal-transfer boundary below a thermodynamic stage."""

    id: str
    from_stage: str
    to_stage: str | None
    device_type: SolidTransferDeviceType
    controlled: bool
    vapor_seal: bool
    fixed_position_pct: float = 50.0
    capacity_factor: float = 1.0


@dataclass(frozen=True)
class ModelConstants:
    dH_vap_hexane: float
    dH_vap_water: float
    T_boil_hexane: float
    T_boil_water: float  # K, from Antoine(antoine_water) at 1 atm — see config/builder.py
    cp_solid: float
    cp_water_liquid: float
    cp_oil: float
    oil_fraction: float  # kg oil/kg dry solid (X3)
    rho_solid: float  # kg/m3, true particle-solid density (particle physics)
    bed_porosity: float  # eps_b, bulk void fraction (zonal transport physics)
    # --- M3a additions ---
    dt_constants: dt_solver.DTSolverConstants
    # dt_resolve_interval_s is NOT here -- it's a HOT, live-tunable value on
    # `Inputs` now (M3a follow-up "C"), not a cold constant frozen at assembly.
    dt_nz_phz: int
    dt_nz_ftrz: int
    dt_nz_dcz: int
    dt_vapor_feed_water_kg_s: float
    dt_vapor_feed_hex_kg_s: float
    dt_vapor_feed_T: float
    dc_constants: dc.DCConstants
    # M3a follow-up ("A2"): real-time-tuned DT solve convergence settings,
    # deliberately looser than dt_solver.solve_dt()'s own validation-run
    # defaults -- see config/schema.py's own comment for the measured
    # numbers behind this choice.
    dt_outer_tol: float
    dt_outer_max_iter: int
    dt_dcz_inner_max_iter: int
    # Legacy D-basis DT tray inventory. None retains old constructor/config behavior.
    packed_non_oil_dry_solid_density_kg_m3: float | None = None
    dt_outer_relaxation: float = 0.5
    # Steam supply-header conditions for the HMI readout (like the DC air readout),
    # shown for both jacket + sparge steam. Display-only; the physics BCs use
    # dt_constants.T_direct_steam is an independent post-valve boundary.
    steam_supply_barg: float = 0.0
    steam_supply_T: float = 0.0  # K, saturation at the steam supply header pressure


@dataclass(frozen=True)
class OperatingSeed:
    """Initial-condition seed for `init_state` (from operating/disturbance
    defaults). M3a: carries enough of `OperatingDefaults` (feed flow rate,
    per-tray steam duties) to run one real `solve_dt` at assembly time
    (BuildSpec §4: "compute steady-state x0 via initializer at the operating
    defaults") — a one-time setup-phase cost, not tick-budget-constrained."""

    feed_temperature: float
    feed_moisture: float
    feed_hexane: float
    feed_oil: float
    feed_flow_rate: float
    indirect_steam: dict[str, float] = field(default_factory=dict)
    direct_steam: dict[str, float] = field(default_factory=dict)
    sweep_arm_speed: dict[str, float] = field(default_factory=dict)
    # DC air-side conditions, so init_state can pre-solve the DRYER/COOLER to steady
    # state too (chain the air-contact equilibrium from the DT exit). Default flows 0
    # -> the DC seed is a no-op (air_contact_equilibrium returns the meal unchanged),
    # a safe fallback for seeds constructed without these (e.g. unit tests).
    heated_air_temp: float = 360.0
    heated_air_flow: float = 0.0
    ambient_air_temp: float = 298.0
    ambient_air_flow: float = 0.0
    ambient_relative_humidity: float = 0.5


@dataclass
class State:
    """Persistent transient state carried between ticks (BuildSpec §7.12), one
    entry per stage, ordered as `Model.stages`."""

    T: np.ndarray  # K
    X1: np.ndarray  # moisture, kg/kg dry solid
    X2: np.ndarray  # hexane, kg/kg dry solid
    M: np.ndarray  # kg dry solid currently retained (bed holdup)
    # M4 (GUI follow-up): effective per-stage solid discharge (kg/s dry) to the
    # stage below (last entry = product rate) AFTER the gate throttle + any
    # downstream back-pressure rejection -- computed in `step()`, read by
    # `outputs()` for the inter-tray flow arrows + product mass balance. Carried
    # on State (not recomputed in outputs) because the back-pressure cascade
    # makes the NET inter-tray flow differ from a stage's own M*k discharge.
    solid_out: np.ndarray
    # --- M3a: cached DT-solve targets + warm-start, one entry per DT-role stage ---
    dt_target_T: np.ndarray
    dt_target_X1: np.ndarray
    dt_target_X2: np.ndarray
    dt_warm_start_vapor_wV2: float
    dt_warm_start_vapor_flow_kg_s: float
    dt_warm_start_vapor_T: float
    dt_warm_start_T_L_sup: float
    dt_last_solve_sim_time: float
    dt_last_attempt_sim_time: float
    dt_converged: bool
    dt_outer_iterations: int
    # M4 (GUI redesign): the full per-cell axial profile from the last periodic
    # solve_dt() resolve -- visualization-only (the HMI's "profile along the
    # tower"), same cadence/lifecycle as dt_target_T/X1/X2 above (only touched
    # in _apply_dt_result, otherwise carried unchanged tick-to-tick).
    dt_axial_profile: dt_solver.DTAxialProfile

    def copy(self) -> "State":
        return State(
            T=self.T.copy(),
            X1=self.X1.copy(),
            X2=self.X2.copy(),
            M=self.M.copy(),
            solid_out=self.solid_out.copy(),
            dt_target_T=self.dt_target_T.copy(),
            dt_target_X1=self.dt_target_X1.copy(),
            dt_target_X2=self.dt_target_X2.copy(),
            dt_warm_start_vapor_wV2=self.dt_warm_start_vapor_wV2,
            dt_warm_start_vapor_flow_kg_s=self.dt_warm_start_vapor_flow_kg_s,
            dt_warm_start_vapor_T=self.dt_warm_start_vapor_T,
            dt_warm_start_T_L_sup=self.dt_warm_start_T_L_sup,
            dt_last_solve_sim_time=self.dt_last_solve_sim_time,
            dt_last_attempt_sim_time=self.dt_last_attempt_sim_time,
            dt_converged=self.dt_converged,
            dt_outer_iterations=self.dt_outer_iterations,
            dt_axial_profile=self.dt_axial_profile,  # frozen/immutable, safe to share
        )


@dataclass
class Inputs:
    """Hot inputs `u` for one tick — MV effective values + DV values (BuildSpec §5.2)."""

    feed_flow_rate: float  # kg/s dry solid
    feed_temperature: float  # K (static boundary at this fidelity; see DECISIONS.md)
    indirect_steam: dict[str, float] = field(default_factory=dict)  # W, per DT stage
    direct_steam: dict[str, float] = field(default_factory=dict)  # kg/s, per SPARGE stage
    sweep_arm_speed: dict[str, float] = field(default_factory=dict)  # rpm, per stage
    transfer_device_position: dict[str, float] = field(
        default_factory=dict
    )  # 0-100 %, per controlled solids-transfer boundary
    heated_air_temp: float = 380.0  # K
    heated_air_flow: float = 0.0  # kg/s
    ambient_air_temp: float = 298.0  # K
    ambient_air_flow: float = 0.0  # kg/s
    feed_moisture: float = 0.0  # kg/kg dry solid
    feed_hexane: float = 0.0  # kg/kg dry solid
    # M4 (GUI redesign): feed oil (X3) is now a LIVE disturbance, not the frozen
    # `ModelConstants.oil_fraction` -- so an operator can vary incoming oil
    # content. Threaded into the DT solve's `SolidFeed.X3` (`_resolve_dt`); the
    # DC ignores oil (a documented ~1% correction, see `core/dc.py::desorb_hexane`).
    # `init_state` still seeds from the constant (runs before any live DV change).
    feed_oil: float = 0.01  # kg/kg dry solid
    # 0-1 fraction (weather RH, not an absolute humidity ratio) -- converted
    # to the absolute humidity `dc.air_contact_equilibrium` actually needs
    # in `_dc_equilibrium` below, evaluated at `ambient_air_temp` (heating
    # this same air parcel up to `heated_air_temp` for the DRYER doesn't
    # change its ABSOLUTE humidity, only its RH at the new temperature --
    # see `_dc_equilibrium`'s own comment).
    ambient_relative_humidity: float = 0.5  # 0-1
    # M3a follow-up ("C"): HOT, live-tunable -- see config/schema.py's
    # OperatingDefaults.dt_resolve_interval_s for the floor rationale.
    dt_resolve_interval_s: float = 120.0  # s, SIM time


@dataclass(frozen=True)
class MassInventory:
    """Very simple, live plant-wide mass-inventory diagnostic (mass/energy
    balance quality gate follow-up, DECISIONS.md) -- NOT a rigorous
    conservation proof (that's `core/balance.py`'s own job, test-suite-only,
    zero runtime cost); this is a cheap, always-on O(n_stages) signal for
    the real-time engine: is anything ACCUMULATING somewhere in the plant.
    `total_*_holdup_kg` are raw snapshots -- the actual "should read ~0 in
    steady state" signal is the CHANGE in these between consecutive
    `Outputs` (a consumer's own tick-to-tick diff; `Model.step` stays pure,
    no history tracked here). `feed_*`/`product_*` rates are exposed
    alongside for context (their own difference is total evaporation/duty
    -driven mass loss, genuinely nonzero in normal operation -- NOT itself
    a conservation signal, unlike the holdup totals)."""

    total_dry_solid_holdup_kg: float
    total_hexane_holdup_kg: float
    total_water_holdup_kg: float
    feed_dry_solid_kg_s: float
    feed_hexane_kg_s: float
    feed_water_kg_s: float
    product_dry_solid_kg_s: float
    product_hexane_kg_s: float
    product_water_kg_s: float


@dataclass
class Outputs:
    """PVs/KPIs for one tick (BuildSpec §9.1 PV/ node map)."""

    stage_T: dict[str, float]
    stage_X_hex_ppm: dict[str, float]
    stage_X_w_pct: dict[str, float]
    stage_vapor_temp: dict[str, float]
    # Total vapor leaving the TOP face of each DT tray, moving upward. This is
    # sampled from the resolved axial profile and drives the HMI's counterflow
    # arrows; DC trays have no entry.
    stage_vapor_flow_kg_s: dict[str, float]
    stage_level_pct: dict[str, float]
    # DC (DRYER/COOLER) air-outlet readout -- keyed only by DC-role stage
    # ids (see `Model.outputs`'s own diagnostic re-derivation); no entry for
    # DT-role stages, which have no equivalent air-side state at this fidelity.
    stage_air_T_out: dict[str, float]
    stage_air_humidity_out: dict[str, float]
    stage_air_hexane_ppm: dict[
        str, float
    ]  # hexane in the DC exhaust air (mole ppm), vs ~1100 LEL limit
    # M4 (GUI redesign): per-stage solid outflow (kg/s dry) = M/tau, so the HMI
    # can annotate the inter-tray flow arrows with a real number (the same
    # M/tau discharge `step()`'s own holdup balance uses).
    stage_solid_out_kg_s: dict[str, float]
    dt_micro_throughflow_kg_s: float
    kpi_residual_hexane_ppm: float
    kpi_meal_moisture_pct: float
    kpi_steam_consumption_kg_per_t: float
    kpi_throughput_t_per_day: float
    # M4 (GUI redesign): energy + vapor-side KPIs for the process dashboard.
    # All derived from the current tick's inputs + the last DT axial-profile
    # resolve; reuse existing constants (dc.CP_AIR_J_KG_K, dH_vap_water/hexane).
    kpi_exhaust_hexane_ppm: float  # DRYER exhaust-air hexane (mole ppm), vs ~1100 LEL limit
    kpi_direct_steam_kg_s: float  # total sparge (direct) steam
    kpi_indirect_heating_kw: float  # total indirect (jacket) heating duty
    kpi_drying_air_heating_kw: float  # duty to heat ambient air up to the dryer setpoint
    kpi_total_energy_kw: float  # indirect + drying-air + direct-steam latent
    kpi_outlet_vapor_kg_s: float  # vapor leaving the DT top (to condenser)
    kpi_outlet_vapor_hexane_kg_s: float  # hexane in that vapor
    kpi_outlet_vapor_water_kg_s: float  # water in that vapor
    kpi_condenser_duty_kw: float  # latent duty to condense the outlet vapor
    dt_solver_converged: bool
    dt_solver_outer_iterations: int
    mass_inventory: MassInventory
    dt_axial_profile: dt_solver.DTAxialProfile
    dt_last_solve_sim_time: float
    dt_last_attempt_sim_time: float
