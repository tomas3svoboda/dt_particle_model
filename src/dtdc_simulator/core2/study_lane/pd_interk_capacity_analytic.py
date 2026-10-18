"""Study-lane knob: the PD inter-K crossing capacity from the analytic layer heat capacity (R3).

Owner ruling of 2026-08-31 (``docs/GT_PS2_FULL_BATCH_RULING_RECORD_2026-08-31.md``, record
id ``GT-PS-2-FULL-BATCH-RULING-RECORD-2026-08-31-01``, decision-index row
``PROC-e7d7424d06``): "B3 - solver/ledger normalization: IMPLEMENT R3 (analytic layer
capacity), removing the class".  R3 is the remedy of
``docs/GT_PS2_PD_INTERK_CROSSING_GUARD_STUDY_2026-08-30.md``; the record of what was built
under B3 (commit 0b417b46) changed the whole-PD ledger instead, and this module is the
missing R3 build, behind its knob, default OFF.  Making it the default replaces the
uniform-layer form the owner made the ruled default on 2026-09-29 and is the owner's.

The guard it serves (``pd_k_dome_bypass_tray_integration``, the inter-K crossing site): the
explicit face heat ``Q = G (T_L - T_U) dt`` must not exceed the room that would carry one
layer onto the other layer's starting temperature.  The former forms measure that room as a
difference of two REBUILT layer energies, ``sum_i [E_i(T_other) - E_i(T_self or stored)]``;
at layers one ULP apart that difference resolves about 3e-7 J out of 1.5e8 J and its codec
rebuild noise (measured up to 2.7e-6 J on the W1-K wall state) flips its sign, so the guard
refuses at every stride a state whose exact two-node evolution is admissible by six orders
of magnitude (``GT_PS2_INTERK_CONDUCTION_WALL_DIAGNOSIS_2026-10-12``).

The ruled form, and nothing else: each layer's room is its ANALYTIC sensible heat capacity
times the separation, ``C_layer |T_L - T_U|``, with

    C_layer = sum_i w_i [ m_dm,i c_core,i(T_i) + m_fw,i du_w,l/dT|P (T_i) ]

``w_i`` the representative particles, ``T_i`` the packet's own decoded temperature, ``m_dm``
and ``m_fw`` the per-particle dry meal and external free water, ``c_core`` the derivative of
the very energy law the packet's codec inverts: tag 2 (film-free, subcritical) the frozen
``wet_core.specific_heat_capacity(T, X_internal, wet)``; tag 1 (film-bearing) the derivative
of ``sphere.combined_high_loading_specific_energy``, which is
``wet_core.specific_heat_capacity(T, X_c(T), wet) + (X_total - X_c(T))
wet_core.hexane_liquid_specific_heat_capacity(T, wet)`` exactly (the two ``dX_c/dT`` terms
cancel because internal and attached n-hexane share one liquid state); the free water
``du/dT|P`` of the same IAPWS liquid state the codec books.  No fit, no constant, no finite
difference.  The datum offsets depend on masses only and carry no temperature derivative.
The crossing capacity is the smaller layer room and the tolerance is the former one, the
ULP sum of the face heat and the three capacities (``_inter_k_crossing_verdict``).  A tag-3
(dry-sorbate) packet has no closed-form derivative here and refuses typed, as the former
level-0 rebuild refuses it.

Report only: when the PD capacity observer is installed, one row per evaluation goes to it
in the existing schema (both former forms, report-only, never raising) with
``form_in_force: "analytic"`` and an ``analytic`` block.  A refusal prints the analytic and
the former capacity.  ``physically_qualifying: false``.
"""

from __future__ import annotations

import math
import os
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from dtdc_simulator.core2 import engineering_subcritical_packet_adapter as tag2
from dtdc_simulator.core2 import high_loading_packet_adapter as tag1
from dtdc_simulator.core2 import pd_k_dome_bypass_tray_integration as pd
from dtdc_simulator.core2.particle import sphere
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.props import critical_volume as cv

PD_INTERK_CAPACITY_ANALYTIC_ENVIRONMENT_VARIABLE = "DTDC_PD_INTERK_CAPACITY_ANALYTIC"
PD_INTERK_CAPACITY_ANALYTIC_DEFAULT = False
RULING_ID = "GT-PS-2-FULL-BATCH-RULING-RECORD-2026-08-31-01"
RULING_DECISION_INDEX_ID = "PROC-e7d7424d06"
RULING_RECORD = "docs/GT_PS2_FULL_BATCH_RULING_RECORD_2026-08-31.md"
PD_INTERK_CAPACITY_ANALYTIC_LAW_ID = (
    "pd-interk-crossing-capacity-analytic-layer-heat-capacity-times-separation-r3-v1"
)
REFUSAL_TEXT = (
    "explicit inter-K conduction would cross the old layer equilibrium; "
    "reject and retry at a shorter endpoint"
)


def knob_armed(environ: Mapping[str, str] | None = None) -> bool:
    """True only for the exact string ``"1"`` (the PD guard site reads the same rule)."""

    env = os.environ if environ is None else environ
    return env.get(PD_INTERK_CAPACITY_ANALYTIC_ENVIRONMENT_VARIABLE) == "1"


class PDInterKAnalyticCapacityRefusal(pd.PDKDomeBypassStepError):
    """A packet whose analytic heat capacity this form cannot give (typed, never substituted)."""


PD_INTERK_CAPACITY_ANALYTIC_DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "pd_interk_capacity_analytic": MappingProxyType(
            {
                "class": pd.DECLARED_ENGINEERING_ASSUMPTION,
                "law_id": PD_INTERK_CAPACITY_ANALYTIC_LAW_ID,
                "status": (
                    "RULED by the owner 2026-08-31 (B3: IMPLEMENT R3, analytic layer "
                    "capacity); built 2026-10-12 behind its knob as a Class-A execution of "
                    "that ruling on the lead's brief (W1-K follow-up); default OFF; making it "
                    "the default replaces the 2026-09-29 ruled uniform-layer default and is "
                    "the owner's (Class B)"
                ),
                "ruling_id": RULING_ID,
                "ruling_decision_index_id": RULING_DECISION_INDEX_ID,
                "ruling_record": RULING_RECORD,
                "ruling_words": (
                    "B3 - solver/ledger normalization: IMPLEMENT R3 (analytic layer "
                    "capacity), removing the class"
                ),
                "knob": PD_INTERK_CAPACITY_ANALYTIC_ENVIRONMENT_VARIABLE,
                "default": "0 (off; armed only by the exact string '1')",
                "former_form": (
                    "crossing capacity = min over the two layers of a rebuilt-and-subtracted "
                    "energy: shipped sum_i [E_i(T_other) - E_i,stored]; uniform-layer (the "
                    "2026-09-29 ruled default) sum_i [E_i(T_other) - E_i(T_self)], each E_i "
                    "re-issued by the codec at the given temperature"
                ),
                "ruled_form": (
                    "crossing capacity = min over the two layers of C_layer |T_L - T_U|, "
                    "C_layer = sum_i w_i [m_dm,i c_core,i(T_i) + m_fw,i du_w,l/dT|P(T_i)], "
                    "c_core the analytic derivative of the energy law the packet's codec "
                    "inverts (tag 2: wet_core.specific_heat_capacity(T, X_internal); tag 1: "
                    "wet_core.specific_heat_capacity(T, X_c(T)) + (X_total - X_c(T)) "
                    "wet_core.hexane_liquid_specific_heat_capacity(T)); no fit, no constant, "
                    "no finite difference; and at the two sites that read that room: site 1, "
                    "the crossing-capacity guard; site 2, the conduction-only order guard's "
                    "uniform layer after the face heat, whose bracket-end room and bisection "
                    "use C_layer (T - T_before) in place of sum_i [E_i(T) - E_i(T_before)]"
                ),
                "measured_case": (
                    "W1-K wall, lane W1I_full step 268 (62.1218 s), PD3: layers 2254.4 kg and "
                    "1286.6 kg dry at the lock 334.4803981611931 K, one ULP apart "
                    "(5.684e-14 K); former (uniform-layer) capacities -5.561e-7 J and "
                    "-3.686e-7 J, rebuilt; true room +4.1e-7 J and +2.3e-7 J; layer heat "
                    "capacities 7.22e6 J/K and 3.96e6 J/K (central difference of the codec "
                    "energy, stable to 1e-6 relative); face heat 1.07e-13 J at 41.7 ms; "
                    "refused 20 times from 41.7 ms to 3.9 ms; explicit crossing stride "
                    "56,543 s.  With site 1 alone on (lane W1M_full5) the same 30 "
                    "evaluations, the wall included, were refused by site 2, whose lower-layer "
                    "bracket-end room is the same rebuilt -5.561e-7 J; the 267 committed rows "
                    "equal W1I_full's"
                ),
                "unchanged": (
                    "the face heat and its conductance, the layer means, the tolerance (the "
                    "ULP sum of the face heat and the three capacities), the refusal text "
                    "and its retry family, the face-heat allocation, the per-packet envelope "
                    "guard, the layer-order guard's tolerance, all-sources reading and lock "
                    "rule, site 2's bracket choice, its 120-halving bisection and its refusal "
                    "text, the lock floor, the shaft-work "
                    "bound, every ledger and digest; knob off, the former guard byte for byte"
                ),
                "linearization": (
                    "the room is C(T_self) |dT|, the first-order term of the exact energy "
                    "difference; the neglected part is about |C'/C| |dT| / 2 of it (5e-4 at "
                    "1 K for the meal), against an explicit stability margin of about 1e5 at "
                    "0.25 s on the measured case"
                ),
                "scope": (
                    "tag-1 and tag-2 packets (and the tag-1-only codec); a tag-3 dry-sorbate "
                    "packet refuses typed (PDInterKAnalyticCapacityRefusal), as the former "
                    "level-0 rebuild refuses it"
                ),
                "rejecting_outcome": (
                    "IDENTITY: knob off, any committed leaf of any lane differing from the "
                    "unedited tree; knob on, any committed leaf differing from the knob-off "
                    "lane before the first evaluation at which the two forms disagree. "
                    "GUARD: a face heat larger than the analytic room plus the unchanged "
                    "tolerance passing; an analytic capacity that is not positive for a "
                    "nonzero separation; any tolerance, budget or retry rule relaxed; a clamp "
                    "in place of the refusal. CALORIC: the analytic layer capacity differing "
                    "from a central difference of the codec's own energy by more than that "
                    "difference's resolution"
                ),
                "physically_qualifying": False,
            }
        )
    }
)


@dataclass(frozen=True)
class AnalyticCrossing:
    """The analytic form's reading of one guard evaluation."""

    lower_heat_capacity_j_k: float | None
    upper_heat_capacity_j_k: float | None
    separation_k: float
    lower_capacity_j: float
    upper_capacity_j: float
    crossing_capacity_j: float
    crossing_bound_j: float
    refuses: bool

    def as_row(self) -> dict[str, object]:
        return {
            "lower_heat_capacity_j_k": self.lower_heat_capacity_j_k,
            "upper_heat_capacity_j_k": self.upper_heat_capacity_j_k,
            "separation_k": self.separation_k,
            "lower_capacity_j": self.lower_capacity_j,
            "upper_capacity_j": self.upper_capacity_j,
            "crossing_capacity_j": self.crossing_capacity_j,
            "crossing_bound_j": self.crossing_bound_j,
            "refuses": self.refuses,
        }


def packet_heat_capacity_j_k(packet: object, *, codec: object) -> float:
    """``dE/dT`` of one weighted packet at its own decoded state, J/K (analytic)."""

    complete = codec.decode(packet.entry.payload_bytes)
    state = complete.state
    temperature = state.temperature_k
    wet = sphere._material_wet_params(
        codec.sphere_params, state.retained_water_loading, state.oil_fraction
    )
    if type(complete) is tag2.CompleteSubcriticalPacket:
        core = wet_core.specific_heat_capacity(temperature, state.internal_hexane_loading, wet)
    elif type(complete) is tag1.CompleteHighLoadingPacket:
        critical = cv.critical_hexane_loading(
            temperature, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
        )
        core = math.fsum(
            (
                wet_core.specific_heat_capacity(temperature, critical, wet),
                (state.total_hexane_loading - critical)
                * wet_core.hexane_liquid_specific_heat_capacity(temperature, wet),
            )
        )
    else:
        raise PDInterKAnalyticCapacityRefusal(
            f"the analytic inter-K capacity has no closed form for a {type(complete).__name__}; "
            "refuse, do not substitute"
        )
    _, water_capacity = wet_core._water_liquid_u_du(temperature, codec.pressure_pa)
    per_particle = math.fsum(
        (
            state.dry_meal_mass_kg * core,
            complete.inventory.external_water_kg * water_capacity,
        )
    )
    capacity = per_particle * packet.entry.weight.representative_particles
    if not (math.isfinite(capacity) and capacity > 0.0):
        raise PDInterKAnalyticCapacityRefusal(
            f"non-positive analytic packet heat capacity {capacity!r} J/K; refuse, do not clamp"
        )
    return capacity


def layer_heat_capacity_j_k(packets: tuple[object, ...], *, codec: object) -> float:
    """The layer's analytic sensible heat capacity, J/K."""

    return math.fsum(packet_heat_capacity_j_k(packet, codec=codec) for packet in packets)


def analytic_crossing(
    bottom_packets: tuple[object, ...],
    upper_packets: tuple[object, ...],
    *,
    lower_temperature_k: float,
    upper_temperature_k: float,
    face_heat_j: float,
    codec: object,
) -> AnalyticCrossing:
    """The ruled form's capacities and verdict; the tolerance is the former one."""

    separation = abs(lower_temperature_k - upper_temperature_k)
    if face_heat_j == 0.0:
        lower_heat = upper_heat = None
        lower_capacity = upper_capacity = 0.0
    else:
        lower_heat = layer_heat_capacity_j_k(bottom_packets, codec=codec)
        upper_heat = layer_heat_capacity_j_k(upper_packets, codec=codec)
        lower_capacity = lower_heat * separation
        upper_capacity = upper_heat * separation
    crossing, bound, refuses = pd._inter_k_crossing_verdict(
        face_heat_j, lower_capacity, upper_capacity
    )
    return AnalyticCrossing(
        lower_heat_capacity_j_k=lower_heat,
        upper_heat_capacity_j_k=upper_heat,
        separation_k=separation,
        lower_capacity_j=lower_capacity,
        upper_capacity_j=upper_capacity,
        crossing_capacity_j=crossing,
        crossing_bound_j=bound,
        refuses=refuses,
    )


def former_crossing(
    bottom_packets: tuple[object, ...],
    upper_packets: tuple[object, ...],
    *,
    lower_temperature_k: float,
    upper_temperature_k: float,
    face_heat_j: float,
    codec: object,
    uniform_layer: bool,
) -> dict[str, object]:
    """One former form's reading, report only: a refusal of its caloric seeds is recorded."""

    try:
        lower, upper = pd._inter_k_crossing_capacities_j(
            bottom_packets,
            upper_packets,
            lower_temperature_k=lower_temperature_k,
            upper_temperature_k=upper_temperature_k,
            face_heat_j=face_heat_j,
            codec=codec,
            uniform_layer=uniform_layer,
        )
        crossing, bound, refuses = pd._inter_k_crossing_verdict(face_heat_j, lower, upper)
    except Exception as error:  # noqa: BLE001 - report only, never raised into the step
        return {"error": f"{type(error).__name__}: {error}"}
    return {
        "lower_capacity_j": lower,
        "upper_capacity_j": upper,
        "crossing_capacity_j": crossing,
        "crossing_bound_j": bound,
        "refuses": refuses,
    }


def require_face_heat_within_analytic_capacity(
    request: object,
    bottom_packets: tuple[object, ...],
    upper_packets: tuple[object, ...],
    *,
    lower_temperature_k: float,
    upper_temperature_k: float,
    face_heat_j: float,
    codec: object,
) -> AnalyticCrossing:
    """The guard under the knob: refuse with the former text, both capacities printed."""

    analytic = analytic_crossing(
        bottom_packets,
        upper_packets,
        lower_temperature_k=lower_temperature_k,
        upper_temperature_k=upper_temperature_k,
        face_heat_j=face_heat_j,
        codec=codec,
    )
    common = dict(
        lower_temperature_k=lower_temperature_k,
        upper_temperature_k=upper_temperature_k,
        face_heat_j=face_heat_j,
        codec=codec,
    )
    former_uniform = pd.PD_INTERK_CAPACITY_UNIFORM_LAYER_ARMED
    former: dict[str, object] | None = None
    observer = pd.PD_INTERK_CAPACITY_OBSERVER
    if observer is not None:
        shipped = former_crossing(bottom_packets, upper_packets, uniform_layer=False, **common)
        uniform = former_crossing(bottom_packets, upper_packets, uniform_layer=True, **common)
        former = uniform if former_uniform else shipped
        prior_tray = getattr(getattr(request, "prior", None), "tray", None)
        observer(
            MappingProxyType(
                {
                    "tray": getattr(prior_tray, "physical_tray_id", None),
                    "t_start_s": getattr(prior_tray, "time_s", None),
                    "t_end_s": getattr(request, "end_time_s", None),
                    "lower_temperature_k": lower_temperature_k,
                    "upper_temperature_k": upper_temperature_k,
                    "face_heat_j": face_heat_j,
                    "form_in_force": "analytic",
                    "shipped": shipped,
                    "uniform_layer": uniform,
                    "analytic": analytic.as_row(),
                }
            )
        )
    if analytic.refuses:
        if former is None:
            former = former_crossing(
                bottom_packets, upper_packets, uniform_layer=former_uniform, **common
            )
        former_text = (
            f"former error {former['error']}"
            if "error" in former
            else f"former capacity {former['crossing_capacity_j']!r} J"
        )
        raise pd.PDKDomeBypassStepError(
            f"{REFUSAL_TEXT} (analytic capacity {analytic.crossing_capacity_j!r} J, "
            f"layer heat capacities {analytic.lower_heat_capacity_j_k!r} and "
            f"{analytic.upper_heat_capacity_j_k!r} J/K; "
            f"{'uniform-layer' if former_uniform else 'shipped'} {former_text}; "
            f"face heat {face_heat_j!r} J)"
        )
    return analytic


def uniform_temperature_after_face_heat_k(
    packets: tuple[object, ...],
    *,
    temperature_before_k: float,
    face_heat_j: float,
    lower_k: float,
    upper_k: float,
    codec: object,
) -> float:
    """Site 2 under the knob: the conduction-only order guard's uniform layer after the face heat.

    The former ``_layer_uniform_temperature_after_face_heat_k`` checks the face heat against
    the rebuilt room ``sum_i [E_i(T) - E_i(T_before)]`` at the bracket ends (the old layer
    envelope) and bisects that room.  At layers one ULP apart the bracket-end room is the very
    quantity the site-1 capacity rebuilds, with the same noise and the same sign (lane
    W1M_full5: the 30 site-2 refusals are the 30 former site-1 refusals, the wall included).
    Here the room is the analytic ``C_layer (T - T_before)``; the bracket (chosen by the caller
    exactly as before), the check, the 120-halving bisection and the refusal text are the
    former ones.  A refusal prints the analytic and the former bracket-end rooms.
    """

    capacity = layer_heat_capacity_j_k(packets, codec=codec)

    def room(temperature_k: float) -> float:
        return capacity * (temperature_k - temperature_before_k)

    if not room(lower_k) <= face_heat_j <= room(upper_k):
        former = []
        for end in (lower_k, upper_k):
            try:
                former.append(
                    repr(
                        pd._layer_uniform_energy_to_temperature_j(
                            packets, end, self_temperature_k=temperature_before_k, codec=codec
                        )
                    )
                )
            except Exception as error:  # noqa: BLE001 - report only
                former.append(f"error {type(error).__name__}: {error}")
        raise pd.PDKDomeBypassStepError(
            f"{REFUSAL_TEXT} (order-guard bracket, analytic room "
            f"[{room(lower_k)!r}, {room(upper_k)!r}] J at layer heat capacity {capacity!r} J/K; "
            f"former room [{former[0]}, {former[1]}] J; face heat {face_heat_j!r} J)"
        )
    lower, upper = lower_k, upper_k
    for _ in range(120):
        middle = 0.5 * (lower + upper)
        if middle == lower or middle == upper:
            break
        if room(middle) < face_heat_j:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


__all__ = (
    "AnalyticCrossing",
    "PDInterKAnalyticCapacityRefusal",
    "PD_INTERK_CAPACITY_ANALYTIC_DECLARED_ASSUMPTIONS",
    "PD_INTERK_CAPACITY_ANALYTIC_DEFAULT",
    "PD_INTERK_CAPACITY_ANALYTIC_ENVIRONMENT_VARIABLE",
    "PD_INTERK_CAPACITY_ANALYTIC_LAW_ID",
    "REFUSAL_TEXT",
    "RULING_ID",
    "analytic_crossing",
    "former_crossing",
    "knob_armed",
    "layer_heat_capacity_j_k",
    "packet_heat_capacity_j_k",
    "require_face_heat_within_analytic_capacity",
    "uniform_temperature_after_face_heat_k",
)
