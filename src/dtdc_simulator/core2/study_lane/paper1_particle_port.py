"""Study-lane knob: the paper-1 particle laws on the tower's packets (W2-C).

Owner ruling 3 of 2026-10-12 (``docs/GT_PS2_OWNER_RULINGS_2026-10-12_OPEN_ITEMS_BATCH.md``,
record id ``GT-PS-2-OWNER-RULINGS-2026-10-12-OPEN-ITEMS-BATCH-01``): "we NEED to have same
particle model as we proposed it in paper 1"; and the owner's standing intention recorded
under ruling D5 of 2026-10-10 (``docs/GT_PS2_OWNER_RULINGS_2026-10-10_REAUDIT_BATCH.md``,
section 4): the paper's laws are ported as the paper states them, in the paper's convention,
with the tower's bridge documented as a conversion.  The design was the W2-C stage-1 note,
approved by the lead on 2026-10-12 with decisions 1A (the laws on the tower's lumped
carrier, plus the paper's radial face law as a library function), 2a (uniform loading,
carrier distance ``d = R/5``) and 3 (the tower's typed floors at ``W_ref`` kept as
tripwires, the continued potential carried inside the law).

WHAT IS PORTED, BY THE PAPER'S EQUATION LABELS (``paper/01_particle_jfpe/sec_formulation.tex``)

* ``eq:luikov`` with its declared positive-water continuation: the retained-water
  potential ``psi = ln a_w = (1 - A1/W)/A2`` for ``0 < W <= W_cap``, evaluated by the pinned
  ``props.sorption.water_log_activity`` on ``ContinuedPositiveLuikovParams`` (never
  re-derived here);
* ``eq:finitecontact``: ``G_l = (d/L_l + 1/K_l)**-1``, ``N_l = -G_l (psi_s - psi_p)`` with
  ``K_l = 0.32 mol/(m2 s)`` and ``L_l = (1000 - 1) L_wr``, ``L_wr`` the paper run's exact
  reference retained-water mobility; ``psi_s`` from the capacity-bounded retained surface
  state (``a_w = 0.799`` at ``W_cap``), never unit liquid activity; absent with no open dry
  shell and with no external liquid;
* ``eq:liquidstocksupply``: the inward transfer bounded by the film stock after the step's
  own gas exchange;
* the capacity ``W_cap = 0.235670904 kg/kg dry`` as an active set (a film above it).

TWO LAYERS, ONE LAW.

1. THE TOWER'S LUMPED CARRIER (decision 1A, 2a).  A tower packet is a lumped inventory with
   one uniform retained loading ``W`` and a hexane front fraction ``f``.  Per packet and per
   layer macro step ``dt``, at the D9-a drain's place in the step
   (``sp1_k2_law2_tray_integration._packet_entries_after``), the drained mass ``m`` solves the
   backward-Euler balance

       m = dt M_w G_l A_p (psi_s - psi(W + m/dry)),   A_p = 3 dry/(rho_dm R),   d = R/5,

   by bisection over the binary64 bit patterns of ``[0, bound]`` (no tolerance knob; the
   function is monotone and has one root), bounded by ``eq:liquidstocksupply`` (the film
   after the advance) and by the open-shell capacity ``max(0, W_cap dry - retained)(1 - f**3)``;
   ``m = 0`` exactly at ``f = 1`` and at a zero film.  The stocks move through the shipped
   ``through_bed_sorbed_water_law.drain_split_kg``, so the D9-a drain function is never
   called for a packet this law serves (the merge packet's measurement 3.5: never both).
   ``psi`` is taken as zero force at and above ``W_cap`` (the capacity active set).
2. THE PAPER'S RADIAL FACE LAW (decision 1A, reusable under option B).
   :func:`liquid_contact_face_flux` and :func:`supply_limited_surface_transfer` are the paper
   harness's flux and limiter bodies (``paper/analysis/harness/
   experiment_core2_radial_fast_wetting.py``) as library functions, with the same arithmetic,
   and :func:`liquid_connected_path` installs them on an explicitly passed cut-transport
   module for the paper's own particle.  No tower path calls them; ``cut_transport.py`` is
   not edited.

THE CONVENTION CONVERSION (documentary only; the law carries 0.32 in the paper's
convention).  Hemmingsen's contact-only construction gives ``K = rho R I / (3 M_w t)`` with
``I = integral_{W0}^{W_end} dW / (psi_s - psi(W))``, hence ``K_tower = K_paper I_T/I_P`` with
``I_P = 1.485197`` and ``I_T = 0.356052`` (ratio 0.239734): paper 0.253170 at 80 C, 0.317809
at 100 C, adopted 0.32; tower 0.060693 and 0.076190; the conversion of the adopted 0.32 is
0.076715.  Construction-bound (it moves with the endpoint), not an invariant.

WHAT IT IS NOT.  Not the paper's radial particle in the tower (option B, the owner's); not a
fit (no coefficient here was adjusted to any tower or outlet number); not a change to any
certified kernel or its default path; the tower's floors at ``W_ref = 0.023 kg/kg dry`` stay
typed refusals.  ``physically_qualifying: false``.
"""

from __future__ import annotations

import json
import math
import os
import struct
from collections.abc import Mapping
from contextlib import contextmanager
from dataclasses import dataclass, replace
from types import MappingProxyType

from dtdc_simulator.core2 import through_bed_sorbed_water_law as sorbed_water
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa

PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE = "DTDC_PARTICLE_PAPER1_IMBIBITION"
PAPER1_PARTICLE_PORT_DEFAULT = False
PAPER1_PARTICLE_PORT_LAW_ID = "paper1-particle-port-finite-contact-imbibition-lumped-carrier-v1"
RULING_ID = "GT-PS-2-OWNER-RULINGS-2026-10-12-OPEN-ITEMS-BATCH-01"
RULING_RECORD = "docs/GT_PS2_OWNER_RULINGS_2026-10-12_OPEN_ITEMS_BATCH.md"
STANDING_INTENTION_ID = "GT-PS-2-OWNER-RULINGS-2026-10-10-REAUDIT-BATCH-01"
STANDING_INTENTION_RECORD = "docs/GT_PS2_OWNER_RULINGS_2026-10-10_REAUDIT_BATCH.md"
PAPER_EQUATION_LABELS = ("eq:luikov", "eq:finitecontact", "eq:liquidstocksupply")
SIDECAR_FILE_NAME = "GT_PS2_PAPER1_PARTICLE_PORT_ROWS.jsonl"

#: The paper's three transport/capacity inputs (``tab:parameters``) and nothing else.
LIQUID_CONTACT_CONDUCTANCE_MOL_M2_S = 0.32
LIQUID_MOBILITY_TOTAL_FACTOR = 1000.0
#: The paper run's exact reference retained-water mobility (printed rounded as 8.47e-7):
#: ``config.dry.retained_water_mobility.value_mol_m_s`` of every archived checkpoint of
#: ``paper/analysis/results_2026-10-08/dt_contact_0p32`` (r3_b_01, c_dt5).
REFERENCE_RETAINED_WATER_MOBILITY_MOL_M_S = 8.474185510409045e-07
CONTINUED_LUIKOV = sp.ContinuedPositiveLuikovParams()
QUALIFIED_LUIKOV = sp.LuikovParams()
W_CAP_KG_KG = CONTINUED_LUIKOV.W_cap
W_REF_TRIPWIRE_KG_KG = QUALIFIED_LUIKOV.W_ref
#: Decision 2a: the lumped carrier's discretisation distance as a fraction of the radius
#: (the linear-driving-force distance of a uniform sphere; merge packet M5: film end within
#: 0.3 to 7.2 percent of the film window, open loop).  Not a physical coefficient.
CARRIER_DISTANCE_FRACTION_OF_RADIUS = 0.2

#: Documentary conversion (evidence: W2-C ``evidence/s1_conversion_numbers.json``).
CONVERSION = MappingProxyType(
    {
        "formula": "K_tower = K_paper * I_T / I_P, I = integral_{W0}^{W_end} dW/(psi_s - psi(W))",
        "I_P": 1.4851970052952101,
        "I_T": 0.3560515264071491,
        "ratio_I_T_over_I_P": 0.23973353375862572,
        "K_paper_80C_mol_m2_s": 0.25317017840013223,
        "K_paper_100C_mol_m2_s": 0.3178093728852724,
        "K_paper_adopted_mol_m2_s": 0.32,
        "K_tower_80C_mol_m2_s": 0.060693381510165394,
        "K_tower_100C_mol_m2_s": 0.07618956402339912,
        "K_tower_100C_endpoint_bracket_mol_m2_s": (0.06586046451190855, 0.08621846437371039),
        "K_tower_of_adopted_0p32_mol_m2_s": 0.07671473080276023,
        "paper_convention": "psi_s = ln a(W_cap) = ln 0.799, W_end = W_cap - 0.01",
        "tower_convention": "psi_s = 0 at x_bf = 0.250 with the bridge 0.799 to 1, W_end = x_bf - 0.01",
        "carried_value": "0.32 in the paper's convention; the tower number is documentary only",
    }
)


def knob_armed(environ: Mapping[str, str] | None = None) -> bool:
    """True only for the exact string ``"1"`` (the two tower sites read the same rule)."""

    env = os.environ if environ is None else environ
    return env.get(PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE) == "1"


class Paper1PortError(sorbed_water.SorbedWaterLawError):
    """A state or configuration the ported law refuses (typed, never substituted)."""


class Paper1PortArmingRefusal(Paper1PortError):
    """The knob is on but the layer does not carry the paper's convention."""


class Paper1PortFloorRefusal(sorbed_water.SorbedWaterFloorRefusal):
    """The tripwire at ``W_ref`` (decision 3): the tower's floors stay typed refusals."""


PAPER1_PARTICLE_PORT_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "paper1_particle_port": MappingProxyType(
            {
                "class": sorbed_water.DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": PAPER1_PARTICLE_PORT_LAW_ID,
                "status": (
                    "RULED by the owner 2026-10-12 (ruling 3, uniformity required) under the "
                    "standing intention of ruling D5 of 2026-10-10; built 2026-10-12 behind its "
                    "knob on the lead's W2-C brief with the lead's stage-1 decisions 1A, 2a, 3; "
                    "default OFF; no paper-3 result printed until this is the particle the "
                    "tower marches (ruling 3)"
                ),
                "ruling_id": RULING_ID,
                "ruling_record": RULING_RECORD,
                "ruling_words": (
                    "I need uniformity, I don't want to block the work, but we NEED to have "
                    "same particle model as we proposed it in paper 1"
                ),
                "standing_intention_id": STANDING_INTENTION_ID,
                "standing_intention_record": STANDING_INTENTION_RECORD,
                "knob": PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE,
                "default": "0 (off; armed only by the exact string '1')",
                "paper_equation_labels": PAPER_EQUATION_LABELS,
                "paper_inputs": {
                    "K_l_mol_m2_s": LIQUID_CONTACT_CONDUCTANCE_MOL_M2_S,
                    "mobility_total_factor": LIQUID_MOBILITY_TOTAL_FACTOR,
                    "L_wr_mol_m_s": REFERENCE_RETAINED_WATER_MOBILITY_MOL_M_S,
                    "W_cap_kg_kg_dry": W_CAP_KG_KG,
                    "luikov_A1_A2": (CONTINUED_LUIKOV.A1, CONTINUED_LUIKOV.A2),
                },
                "carrier": (
                    "lumped packet, uniform retained loading (row d9a_uniform_loading kept), "
                    "backward Euler per layer macro step, bisection over binary64 bit patterns, "
                    "bounded by the post-advance film (eq:liquidstocksupply) and the open-shell "
                    "capacity max(0, W_cap dry - retained)(1 - f**3); d = R/5"
                ),
                "convention": (
                    "the paper's: capacity W_cap with the Luikov activity up to and including "
                    "W_cap (a_w = 0.799 there), a film above it; the declared arm must carry "
                    "x_bf == W_cap exactly and the destination 'sorbed'"
                ),
                "conversion": dict(CONVERSION),
                "replaces": (
                    "for every packet of a layer with a declared arm: the D9-a nominal drain "
                    "(film_drain_kg is not called) and the film-free closure's arm.activity"
                ),
                "not_ported": (
                    "the paper's radial discretisation and regime A/B/C event topology "
                    "(option B, the owner's), the |beta| <= 3 band, the RULE-2 numerical "
                    "allowances, the prescribed legacy bath, the Faner multiplier 0.146"
                ),
                "tripwires_kept": (
                    "the tower's W_ref = 0.023 kg/kg dry floors (sphere._material_wet_params, "
                    "the D9 exhaustion floors) and this law's own refusal below W_ref"
                ),
                "rejecting_outcomes": (
                    "(i) the closed-loop lumped carrier on the paper's supply history misses "
                    "the paper's film depletion by more than 25 percent of the film window; "
                    "(ii) any tower packet below 0.023 kg/kg dry (refuses typed); "
                    "(iii) a packet with a water film and an open shell on a step that does not "
                    "reach this site (a tray step outside the SP1 K=2 kernel); "
                    "(iv) an implied shell loading above W_cap from the per-packet rows "
                    "(caveat 1 of D5; triggers decision 2 option b)"
                ),
                "physically_qualifying": False,
            }
        )
    }
)


# ---------------------------------------------------------------------------------------
# The paper's convention: arming, potential, surface activity
# ---------------------------------------------------------------------------------------


def require_paper_convention_arm(arm: sorbed_water.SorbedWaterArm | None) -> None:
    """Refuse unless the layer's declared arm carries the paper's convention exactly."""

    if arm is None:
        raise Paper1PortArmingRefusal(
            f"{PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE}=1 needs a declared sorbed-water arm "
            "on every layer it serves; this layer has none (the paper's imbibition law has "
            "no tower carrier without the arm) - refuse, do not fall back"
        )
    if type(arm) is not sorbed_water.SorbedWaterArm:
        raise Paper1PortArmingRefusal("the declared arm must be an exact SorbedWaterArm")
    if arm.x_bf != W_CAP_KG_KG:
        raise Paper1PortArmingRefusal(
            f"the paper's convention places the retained-water capacity at W_cap = "
            f"{W_CAP_KG_KG!r} kg/kg dry; the declared arm carries x_bf = {arm.x_bf!r} - "
            "declare x_bf = W_cap or switch the knob off"
        )
    if arm.condensate_destination != sorbed_water.CONDENSATE_DESTINATION_NOMINAL:
        raise Paper1PortArmingRefusal(
            "the paper's imbibition law replaces the 'sorbed' drain; the declared arm "
            f"carries {arm.condensate_destination!r}"
        )


def retained_potential(loading_kg_kg_dry: float) -> float:
    """``eq:luikov``: ``psi = (1 - A1/W)/A2`` on ``0 < W <= W_cap`` (pinned, continued type)."""

    return sp.water_log_activity(loading_kg_kg_dry, CONTINUED_LUIKOV)


def surface_potential() -> float:
    """``psi_s``: the capacity-bounded retained surface trace, ``ln a_w(W_cap)``."""

    return retained_potential(W_CAP_KG_KG)


def retained_surface_activity(
    arm: sorbed_water.SorbedWaterArm | None, loading_kg_kg_dry: float
) -> float:
    """The film-free layer's interface water activity in the paper's convention.

    The pinned ``water_activity`` on the continued type for ``W_ref <= W <= W_cap``
    (identical in value to the qualified law on that band; 0.799 at ``W_cap``, where the
    shipped arm at ``x_bf = W_cap`` would return 1.0); a typed refusal below ``W_ref``
    (decision 3's tripwire) and above ``W_cap`` (the capacity active set: a film-free
    layer cannot hold water above its capacity in this convention).
    """

    require_paper_convention_arm(arm)
    loading = loading_kg_kg_dry
    if type(loading) is not float or not math.isfinite(loading):
        raise Paper1PortError("the retained-water loading must be a finite binary64")
    if loading < W_REF_TRIPWIRE_KG_KG:
        raise Paper1PortFloorRefusal(
            f"retained-water loading {loading!r} kg/kg dry is below the tower's tripwire "
            f"{W_REF_TRIPWIRE_KG_KG} (decision 3: the tower floors stay typed refusals)"
        )
    if loading > W_CAP_KG_KG:
        raise Paper1PortError(
            f"retained-water loading {loading!r} kg/kg dry is above the paper's capacity "
            f"{W_CAP_KG_KG}: in this convention the excess is a film - refuse, do not clamp"
        )
    return sp.water_activity(loading, CONTINUED_LUIKOV)


# ---------------------------------------------------------------------------------------
# eq:finitecontact on the lumped carrier
# ---------------------------------------------------------------------------------------


def added_liquid_mobility_mol_m_s() -> float:
    """``L_l = (factor - 1) L_wr``, the harness's own arithmetic."""

    return (LIQUID_MOBILITY_TOTAL_FACTOR - 1) * REFERENCE_RETAINED_WATER_MOBILITY_MOL_M_S


def series_conductance_mol_m2_s(distance_m: float) -> float:
    """``eq:finitecontact``: ``G_l = (d/L_l + 1/K_l)**-1``; ``d = 0`` gives ``K_l``."""

    if type(distance_m) is not float or not math.isfinite(distance_m) or distance_m < 0.0:
        raise Paper1PortError("the carrier distance must be a finite non-negative binary64")
    return 1.0 / (
        distance_m / added_liquid_mobility_mol_m_s() + 1.0 / LIQUID_CONTACT_CONDUCTANCE_MOL_M2_S
    )


def carrier_distance_m(radius_m: float) -> float:
    """Decision 2a: ``d = R/5``."""

    return CARRIER_DISTANCE_FRACTION_OF_RADIUS * radius_m


@dataclass(frozen=True, slots=True)
class Paper1Drain:
    """One packet's drain under the ported law, with everything a reader needs to check it."""

    drain_kg: float
    film_kg: float
    capacity_kg: float
    bound_kg: float
    limited_by: str
    conductance_mol_m2_s: float
    outer_area_m2: float
    force_old: float
    force_new: float
    dissipation_w_k: float
    bisection_steps: int


def _bits(value: float) -> int:
    return struct.unpack("<q", struct.pack("<d", value))[0]


def _from_bits(bits: int) -> float:
    return struct.unpack("<d", struct.pack("<q", bits))[0]


def finite_contact_drain(
    *,
    arm: sorbed_water.SorbedWaterArm | None,
    film_kg: float,
    dry_matter_kg: float,
    retained_water_kg: float,
    front_fraction: float,
    macro_step_s: float,
    radius_m: float,
    rho_dm_kg_m3: float,
) -> Paper1Drain:
    """The ported law's film-to-shell transfer for one packet over one macro step."""

    require_paper_convention_arm(arm)
    for name, value in (
        ("film_kg", film_kg),
        ("dry_matter_kg", dry_matter_kg),
        ("retained_water_kg", retained_water_kg),
    ):
        if type(value) is not float or not math.isfinite(value) or value < 0.0:
            raise Paper1PortError(f"{name} must be a finite non-negative binary64")
    for name, value in (
        ("macro_step_s", macro_step_s),
        ("radius_m", radius_m),
        ("rho_dm_kg_m3", rho_dm_kg_m3),
    ):
        if type(value) is not float or not math.isfinite(value) or value <= 0.0:
            raise Paper1PortError(f"{name} must be a finite positive binary64")
    if dry_matter_kg <= 0.0:
        raise Paper1PortError("the packet has no dry matter to carry the law")
    loading = retained_water_kg / dry_matter_kg
    if loading < W_REF_TRIPWIRE_KG_KG:
        raise Paper1PortFloorRefusal(
            f"retained-water loading {loading!r} kg/kg dry is below the tower's tripwire "
            f"{W_REF_TRIPWIRE_KG_KG} (decision 3)"
        )
    capacity_kg = sorbed_water.shell_capacity_kg(
        dry_matter_kg=dry_matter_kg,
        retained_water_kg=retained_water_kg,
        front_fraction=front_fraction,
        x_bf=W_CAP_KG_KG,
    )
    bound_kg = min(film_kg, capacity_kg)
    conductance = series_conductance_mol_m2_s(carrier_distance_m(radius_m))
    area = 3.0 * dry_matter_kg / (rho_dm_kg_m3 * radius_m)
    psi_s = surface_potential()

    def force(mass_kg: float) -> float:
        new_loading = math.fsum((retained_water_kg, mass_kg)) / dry_matter_kg
        if new_loading >= W_CAP_KG_KG:
            return 0.0  # the capacity active set: no uptake at or above W_cap
        return psi_s - retained_potential(new_loading)

    force_old = force(0.0)
    dissipation = bg.R * conductance * force_old * force_old * area
    if bound_kg <= 0.0 or force_old <= 0.0:
        return Paper1Drain(
            0.0, film_kg, capacity_kg, bound_kg,
            "no_open_shell_or_no_film" if bound_kg <= 0.0 else "no_force",
            conductance, area, force_old, force_old, dissipation, 0,
        )  # fmt: skip
    scale = macro_step_s * wa.M * conductance * area

    def residual(mass_kg: float) -> float:
        return mass_kg - scale * force(mass_kg)

    if residual(bound_kg) <= 0.0:
        limited = "film_stock" if bound_kg == film_kg else "shell_capacity"
        return Paper1Drain(
            bound_kg, film_kg, capacity_kg, bound_kg, limited, conductance, area,
            force_old, force(bound_kg), dissipation, 0,
        )  # fmt: skip
    low, high = 0, _bits(bound_kg)
    steps = 0
    while high - low > 1:
        middle = (low + high) // 2
        if residual(_from_bits(middle)) < 0.0:
            low = middle
        else:
            high = middle
        steps += 1
    drained = _from_bits(low)
    return Paper1Drain(
        drained, film_kg, capacity_kg, bound_kg, "rate", conductance, area,
        force_old, force(drained), dissipation, steps,
    )  # fmt: skip


#: In-process tally (report only).
TALLY: dict[str, int] = {
    "site_calls": 0,
    "nonzero_drains": 0,
    "limited_by_rate": 0,
    "limited_by_film_stock": 0,
    "limited_by_shell_capacity": 0,
    "activity_calls": 0,
}


_DECLARATION_STAMPED: set[str] = set()


def _sidecar(row: dict) -> None:
    """Append one report row; the first row of each process and file is the declaration."""

    directory = os.environ.get("DTDC_EVIDENCE_DIR")
    if not directory:
        return
    path = os.path.join(directory, SIDECAR_FILE_NAME)
    with open(path, "a", encoding="utf-8") as handle:
        if path not in _DECLARATION_STAMPED:
            _DECLARATION_STAMPED.add(path)
            declaration = {
                name: dict(block)
                for name, block in PAPER1_PARTICLE_PORT_DECLARED_ASSUMPTIONS.items()
            }
            stamp = {
                "row": "declaration (stamped once per process)",
                "knob": PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE,
                "knob_value": os.environ.get(PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE),
                "declared_assumptions": declaration,
            }
            handle.write(json.dumps(stamp, sort_keys=True, default=list) + "\n")
        handle.write(json.dumps(row, sort_keys=True) + "\n")


def drain_split_at_site(
    *,
    arm: sorbed_water.SorbedWaterArm | None,
    film_kg: float,
    dry_matter_kg: float,
    retained_water_kg: float,
    front_fraction: float,
    macro_step_s: float,
    sphere_params,
    layer: int | None = None,
    packet_id: str | None = None,
) -> tuple[float, float]:
    """The tower site's call: the drain, then the shipped exact split ``(film, retained)``."""

    TALLY["site_calls"] += 1
    result = finite_contact_drain(
        arm=arm,
        film_kg=film_kg,
        dry_matter_kg=dry_matter_kg,
        retained_water_kg=retained_water_kg,
        front_fraction=front_fraction,
        macro_step_s=macro_step_s,
        radius_m=float(sphere_params.particle_radius_m),
        rho_dm_kg_m3=float(sphere_params.wet_core.rho_dm_p),
    )
    if result.drain_kg > 0.0:
        TALLY["nonzero_drains"] += 1
        TALLY[f"limited_by_{result.limited_by}"] += 1
    if film_kg > 0.0 or result.drain_kg > 0.0:
        _sidecar(
            {
                "law_id": PAPER1_PARTICLE_PORT_LAW_ID,
                "row": "attempt (report only; the step may still be rejected and retried)",
                "layer": layer,
                "packet_id": None if packet_id is None else str(packet_id),
                "front_fraction": front_fraction,
                "dry_matter_kg": dry_matter_kg,
                "retained_water_kg": retained_water_kg,
                "macro_step_s": macro_step_s,
                "film_kg": result.film_kg,
                "capacity_kg": result.capacity_kg,
                "drain_kg": result.drain_kg,
                "limited_by": result.limited_by,
                "conductance_mol_m2_s": result.conductance_mol_m2_s,
                "force_old": result.force_old,
                "force_new": result.force_new,
                "dissipation_w_k": result.dissipation_w_k,
                "bisection_steps": result.bisection_steps,
            }
        )
    return sorbed_water.drain_split_kg(
        film_kg=film_kg, retained_water_kg=retained_water_kg, drain_kg=result.drain_kg
    )


def closure_activity_at_site(arm: sorbed_water.SorbedWaterArm | None, loading: float) -> float:
    """The tower site's call for the film-free closure's interface activity."""

    TALLY["activity_calls"] += 1
    return retained_surface_activity(arm, loading)


# ---------------------------------------------------------------------------------------
# The paper's radial face law as library functions (the harness's arithmetic)
# ---------------------------------------------------------------------------------------


def add_retained_transfer(native, addition, enthalpy, dissipation):
    """Book one additional liquid flux in the component, energy and entropy rows."""

    component = replace(
        native.component,
        retained_water_flux_mol_m2_s=native.component.retained_water_flux_mol_m2_s + addition,
        conserved_water_flux_mol_m2_s=native.component.conserved_water_flux_mol_m2_s + addition,
    )
    energy = replace(
        native.energy,
        retained_water_enthalpy_flux_w_m2=native.energy.retained_water_enthalpy_flux_w_m2
        + addition * enthalpy,
        total_energy_flux_w_m2=native.energy.total_energy_flux_w_m2 + addition * enthalpy,
    )
    return replace(
        native,
        component=component,
        energy=energy,
        retained_water_flux_mol_m2_s=native.retained_water_flux_mol_m2_s + addition,
        retained_water_entropy_production_w_m3_k=native.retained_water_entropy_production_w_m3_k
        + dissipation,
    )


def _check_face_parameters(factor, contact_conductance) -> None:
    if not math.isfinite(factor) or factor < 1:
        raise ValueError("fast-wetting multiplier must be finite and >= 1")
    if contact_conductance is not None and (
        not math.isfinite(contact_conductance) or contact_conductance <= 0
    ):
        raise ValueError("liquid contact conductance must be finite and positive")


def liquid_contact_face_flux(
    original,
    args,
    kwargs,
    *,
    factor,
    exterior_liquid_present,
    contact_conductance=None,
    exclude_wet_front=False,
    surface_pairs=None,
):
    """``eq:finitecontact`` on one radial dry face: the native flux plus the liquid path.

    ``original`` is the native ``cut_transport._dry_flux`` and ``args``/``kwargs`` its call.
    Interior faces carry ``(factor - 1) L_wr`` on the actual retained activities; the face
    with the surface trace (``thermodynamic_force_temperature_k`` given) carries the series
    ``G_l``; the wet-front trace is excluded on request.  ``surface_pairs``, when given,
    records the surface face for :func:`supply_limited_surface_transfer`.
    """

    is_wet_front = "temperature_gradient_k_m" in kwargs or "composition_gradient_m_inv" in kwargs
    if factor == 1 or not exterior_liquid_present or (exclude_wet_front and is_wet_front):
        return original(*args, **kwargs)
    config = args[8]
    native = original(*args, **kwargs)
    left = cp.evaluate_equilibrium(args[0], config.dry.pressure_pa, args[1], args[6])
    right = cp.evaluate_equilibrium(args[2], config.dry.pressure_pa, args[3], args[7])
    aw_left = sp.water_activity(left.retained_water_loading, args[6].luikov)
    aw_right = sp.water_activity(right.retained_water_loading, args[7].luikov)
    gradient = math.log(aw_right / aw_left) / args[4]
    mobility = (factor - 1) * config.dry.retained_water_mobility.value_mol_m_s
    if contact_conductance is not None and "thermodynamic_force_temperature_k" in kwargs:
        mobility = args[4] / (args[4] / mobility + 1.0 / contact_conductance)
    addition = -mobility * gradient
    donor = left if addition >= 0 else right
    dissipation = bg.R * mobility * gradient * gradient
    augmented = add_retained_transfer(
        native, addition, donor.retained_water_enthalpy_j_mol, dissipation
    )
    if surface_pairs is not None and "thermodynamic_force_temperature_k" in kwargs:
        surface_pairs[id(augmented)] = (
            augmented,
            native,
            addition,
            donor.retained_water_enthalpy_j_mol,
            dissipation,
        )
    return augmented


def inventory_bounded_transfer_fraction(available_flux, demand_flux, dt_s):
    """``eq:liquidstocksupply`` in mass units, through D9's own ``film_drain_kg``."""

    if available_flux <= 0:
        return 0.0
    requested_mass = demand_flux * dt_s
    drained_mass = sorbed_water.film_drain_kg(
        film_kg=float(available_flux * dt_s),
        capacity_kg=float(requested_mass),
        bracket_end=sorbed_water.DRAIN_BRACKET_END_NOMINAL,
    )
    return drained_mass / requested_mass


def supply_limited_surface_transfer(
    face,
    surface_pairs,
    *,
    factor,
    exterior_liquid_present,
    old_liquid_loading_kg_m2,
    dt_s,
    external_water_flux_kg_m2_s,
    requested_share=None,
):
    """``eq:liquidstocksupply`` at the surface face; returns ``(face, audit or None)``."""

    if factor == 1 or not exterior_liquid_present:
        return face, None
    _, native, addition, enthalpy, dissipation = surface_pairs[id(face)]
    available = (
        old_liquid_loading_kg_m2 / dt_s
        + native.component.conserved_water_flux_mol_m2_s * wa.M
        - external_water_flux_kg_m2_s
    )
    demand = -addition * wa.M
    if demand <= 0:
        return native, None
    bounded_share = inventory_bounded_transfer_fraction(available, demand, dt_s)
    share = bounded_share if requested_share is None else requested_share
    if not 0 <= share <= 1:
        raise ValueError("wet-contact transfer fraction outside [0,1]")
    transferred = demand * share
    audit = {
        "available_water_flux_kg_m2_s": float(available),
        "requested_liquid_uptake_kg_m2_s": float(demand),
        "transferred_liquid_uptake_kg_m2_s": float(transferred),
        "fraction_of_requested_path": float(share),
        "D9_supply_bounded_fraction": float(bounded_share),
        "old_liquid_loading_kg_m2": float(old_liquid_loading_kg_m2),
        "dt_s": float(dt_s),
        "D9_inventory_bounded_transfer_used": True,
    }
    return add_retained_transfer(native, addition * share, enthalpy, dissipation * share), audit


@contextmanager
def liquid_connected_path(
    cut_module,
    factor,
    *,
    exterior_liquid_present,
    contact_conductance=None,
    exclude_wet_front=False,
):
    """The paper particle's liquid path on an EXPLICITLY passed cut-transport module.

    For the paper's own radial particle (its harness and the replay tests); no tower path
    calls it.  Rebinds ``cut_module._dry_flux`` for the duration and restores it; yields the
    surface limiter with a ``last_audit`` attribute, as the harness's context manager does.
    """

    _check_face_parameters(factor, contact_conductance)
    original = cut_module._dry_flux
    surface_pairs: dict = {}

    def flux(*args, **kwargs):
        return liquid_contact_face_flux(
            original,
            args,
            kwargs,
            factor=factor,
            exterior_liquid_present=exterior_liquid_present,
            contact_conductance=contact_conductance,
            exclude_wet_front=exclude_wet_front,
            surface_pairs=surface_pairs,
        )

    def limiter(
        face,
        *,
        old_liquid_loading_kg_m2,
        dt_s,
        external_water_flux_kg_m2_s,
        requested_share=None,
    ):
        result, audit = supply_limited_surface_transfer(
            face,
            surface_pairs,
            factor=factor,
            exterior_liquid_present=exterior_liquid_present,
            old_liquid_loading_kg_m2=old_liquid_loading_kg_m2,
            dt_s=dt_s,
            external_water_flux_kg_m2_s=external_water_flux_kg_m2_s,
            requested_share=requested_share,
        )
        if audit is not None:
            limiter.last_audit = audit
        return result

    limiter.last_audit = None
    cut_module._dry_flux = flux
    try:
        yield limiter
    finally:
        cut_module._dry_flux = original


__all__ = [
    "CARRIER_DISTANCE_FRACTION_OF_RADIUS",
    "CONVERSION",
    "LIQUID_CONTACT_CONDUCTANCE_MOL_M2_S",
    "LIQUID_MOBILITY_TOTAL_FACTOR",
    "PAPER1_PARTICLE_PORT_DECLARED_ASSUMPTIONS",
    "PAPER1_PARTICLE_PORT_DEFAULT",
    "PAPER1_PARTICLE_PORT_ENVIRONMENT_VARIABLE",
    "PAPER1_PARTICLE_PORT_LAW_ID",
    "PAPER_EQUATION_LABELS",
    "Paper1Drain",
    "Paper1PortArmingRefusal",
    "Paper1PortError",
    "Paper1PortFloorRefusal",
    "REFERENCE_RETAINED_WATER_MOBILITY_MOL_M_S",
    "TALLY",
    "W_CAP_KG_KG",
    "add_retained_transfer",
    "carrier_distance_m",
    "closure_activity_at_site",
    "drain_split_at_site",
    "finite_contact_drain",
    "inventory_bounded_transfer_fraction",
    "knob_armed",
    "liquid_connected_path",
    "liquid_contact_face_flux",
    "require_paper_convention_arm",
    "retained_potential",
    "retained_surface_activity",
    "series_conductance_mol_m2_s",
    "supply_limited_surface_transfer",
    "surface_potential",
]
