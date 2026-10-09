"""T3a â€” whole-DTDC tray-stack assembly on the firm reference machine.

Builds the machine of `docs/GT_PS2_T3_REFERENCE_MACHINE_2026-08-07.md`
(owner-accepted 2026-08-07) as ONE marched cell stack and reports
whole-machine mass and energy ledgers.

Architecture, per the frozen ledger and the owner's rulings:

- **Solid descends, vapour rises.**  Feed enters the top predesolventizer
  tray; meal falls PD -> MN -> SP and discharges at the sparge tray.
  Sparge steam enters at SP1 and rises countercurrently.
- **Two vapour-path modes**, because the hardware differs:
  * `THROUGH_BED` (countercurrent + sparge trays) â€” perforated, vapour
    passes through the bed, so the exchange area is the PHY-023 packed-bed
    interfacial area (1-eps)(3/r_p)V.
  * `DOME_CONTACT` (predesolventizer trays) â€” **the owner's ruling of
    2026-08-07**: PD trays are UNPERFORATED DISCS, so vapour cannot pass
    through them, but the header dome is open space and the meal top
    surface is exposed.  The exchange area is therefore the tray free
    surface (the cross-section), NOT the interfacial area.  At the
    reference machine that is 28.27 m2 against ~17250 m2 â€” a ratio near
    1:610, which is the quantitative content of "significantly less
    than countercurrent trays".  The contact fraction is DERIVED from
    geometry, never fitted.
- Each tray carries its own jacket duty and its own holdup (bed depths
  differ), both supported by `FilmDepletionSequence`'s per-cell forms.

DISCLOSED APPROXIMATION (T3a tranche): the DOME_CONTACT cells reuse the
packed-bed Coletto film coefficient at the free surface.  The area is
right and dominates (610x); the coefficient basis is not â€” a free
convection / impingement closure is the proper one and is named as the
next refinement.  Recorded rather than hidden.

physically_qualifying is False throughout.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass

from .film import TrayLayerFilm
from .film_depletion_sequence import (
    FilmDepletionError,
    FilmDepletionSequence,
)
from .tray import GasStream, SolidStream, TrayParams


class VaporContact(enum.Enum):
    """How a tray's vapour meets its meal."""

    THROUGH_BED = "THROUGH_BED"      # perforated: packed-bed interfacial area
    DOME_CONTACT = "DOME_CONTACT"    # unperforated disc: free surface only


@dataclass(frozen=True)
class TraySpec:
    """One physical tray of the reference machine."""

    tray_id: str
    role: str                    # PREDESOLV | MAIN | SPARGE
    contact: VaporContact
    diameter_m: float
    loaded_depth_m: float
    indirect_duty_w: float

    @property
    def cross_section_m2(self) -> float:
        return math.pi * self.diameter_m ** 2 / 4.0

    @property
    def bed_volume_m3(self) -> float:
        return self.cross_section_m2 * self.loaded_depth_m


@dataclass(frozen=True)
class DTDCStackResult:
    """One converged whole-machine solve."""

    tray_ids: tuple[str, ...]
    exchange_area_m2: tuple[float, ...]
    holdup_kg: tuple[float, ...]
    duty_w: tuple[float, ...]
    # boundary streams
    feed: SolidStream
    discharge_x_h: float
    discharge_x_w: float
    discharge_temperature_k: float
    dome_vapor: GasStream
    sparge_steam: GasStream
    # ledgers
    hexane_ledger_residual_rel: float
    water_ledger_residual_rel: float
    energy_ledger_residual_w: float
    energy_ledger_residual_rel: float
    max_layer_energy_residual_w: float
    converged: bool
    macro_steps: int
    energy_screening_tier: bool = False

    @property
    def physically_qualifying(self) -> bool:
        return False


# --- reference machine, 3000 t/d raw soybean -------------------------
# Every number traceable to GT-PS-2-T3-REFERENCE-MACHINE (owner-accepted).
REFERENCE_TRAYS: tuple[TraySpec, ...] = (
    TraySpec("PD1", "PREDESOLV", VaporContact.DOME_CONTACT, 6.0, 0.30, 233.3e3),
    TraySpec("PD2", "PREDESOLV", VaporContact.DOME_CONTACT, 6.0, 0.30, 233.3e3),
    TraySpec("PD3", "PREDESOLV", VaporContact.DOME_CONTACT, 6.0, 0.30, 233.4e3),
    TraySpec("MN1", "MAIN", VaporContact.THROUGH_BED, 6.0, 1.10, 280.0e3),
    TraySpec("MN2", "MAIN", VaporContact.THROUGH_BED, 6.0, 1.00, 280.0e3),
    TraySpec("SP1", "SPARGE", VaporContact.THROUGH_BED, 6.0, 0.90, 140.0e3),
)
REFERENCE_BED_VOID = 0.40
REFERENCE_PARTICLE_RADIUS_M = 0.885e-3
REFERENCE_BULK_DENSITY_KG_M3 = 600.0   # bounded by Kong 720 tapped; disclosed
REFERENCE_CP_GAS_J_KG_K = 1926.0
# QG-6E assumption A-05: shell/ambient loss as a fraction of gross
# incoming thermal duty.  Nominal 1 %, sensitivity interval 0-2 %.
A05_LOSS_FRACTION = 0.01


def exchange_area_m2(tray: TraySpec) -> float:
    """The area vapour actually sees â€” the owner's PD ruling made
    quantitative."""
    if tray.contact is VaporContact.DOME_CONTACT:
        return tray.cross_section_m2
    return bed_side_area_m2(tray)


def bed_side_area_m2(tray: TraySpec) -> float:
    """Bed-side (jacket -> interface) contact area: ALWAYS the packed-bed
    interfacial area.  A predesolventizer bed boils through its whole
    volume against its own jacket regardless of whether the tray is
    perforated; restricting this path as well would starve the bed of
    its own duty and is physically wrong."""
    return (
        (1.0 - REFERENCE_BED_VOID)
        * (3.0 / REFERENCE_PARTICLE_RADIUS_M)
        * tray.bed_volume_m3
    )


def tray_holdup_kg(tray: TraySpec) -> float:
    return tray.bed_volume_m3 * REFERENCE_BULK_DENSITY_KG_M3


class DTDCStack:
    """The whole DT as one marched cell stack (T3a)."""

    physically_qualifying = False

    def __init__(
        self,
        *,
        trays: tuple[TraySpec, ...] = REFERENCE_TRAYS,
        bed_contact_h_w_m2_k: float = 106.5,
        gas_superficial_velocity_m_s: float = 0.354,
        macro_steps_per_residence: int = 15,
        # Plant-scale trays turn over slowly against the 1e-6 discharge
        # stability criterion: a measured single-tray solve needed ~15
        # residences (228 macro-steps), so the 8 that suffice for the
        # small S1/S2 fixtures are not enough here.
        max_residences: float = 30.0,
        pressure_pa: float = 101325.0,
        sparge_contact_temperature_k: float = 373.15,
        energy_screening_tier: bool = False,
    ) -> None:
        if not trays:
            raise FilmDepletionError("the stack needs at least one tray")
        self._trays = trays
        # Cell index 0 is the BOTTOM of the marched stack (gas enters
        # there), so the physical machine is reversed: SP1 is cell 0 and
        # the top predesolventizer tray is the last cell.
        self._cells = tuple(reversed(trays))
        self._contact_h = bed_contact_h_w_m2_k
        self._velocity = gas_superficial_velocity_m_s
        self._steps = macro_steps_per_residence
        self._residences = max_residences
        self._pressure = pressure_pa
        self._sparge_contact_temperature_k = sparge_contact_temperature_k
        self._screening_tier = bool(energy_screening_tier)

    def run(self, *, feed: SolidStream, sparge_steam: GasStream) -> DTDCStackResult:
        cells = self._cells
        areas = tuple(bed_side_area_m2(t) for t in cells)
        gas_areas = tuple(exchange_area_m2(t) for t in cells)
        holdups = tuple(tray_holdup_kg(t) for t in cells)
        duties = list(t.indirect_duty_w for t in cells)

        # --- SPARGE DESUPERHEAT (physically required, energy-exact) ---
        # Design sparge steam is 150-160 C.  Against the ~334.5 K
        # azeotrope-pinned interface that is a relative temperature jump
        # of ~0.246, and the FROZEN reduced film refuses above 0.20 â€”
        # the limitation predicted when the Karnofsky gate was built.
        # Physically the steam does not touch meal at 155 C: it is
        # injected into the vapour space and desuperheats by mixing
        # almost immediately.  So the stream CONTACTS the bed at the
        # local vapour temperature while its superheat is delivered to
        # the sparge tray as heat.  Mass and energy are both preserved
        # exactly; only the delivery path changes.
        contact_t = min(sparge_steam.T, self._sparge_contact_temperature_k)
        steam_mass = sparge_steam.G_w + sparge_steam.G_h
        superheat_w = (
            steam_mass * REFERENCE_CP_GAS_J_KG_K * (sparge_steam.T - contact_t)
        )
        duties[0] += superheat_w
        contact_steam = GasStream(
            G_w=sparge_steam.G_w, G_h=sparge_steam.G_h, T=contact_t
        )
        duties = tuple(duties)
        films = tuple(
            TrayLayerFilm(
                particle_radius_m=REFERENCE_PARTICLE_RADIUS_M,
                bed_void_fraction=REFERENCE_BED_VOID,
                gas_density_kg_m3=0.6,
                gas_superficial_velocity_m_s=self._velocity,
                gas_dynamic_viscosity_pa_s=1.329e-5,
                gas_heat_capacity_j_kg_k=REFERENCE_CP_GAS_J_KG_K,
                gas_thermal_conductivity_w_m_k=0.02371,
                pressure_pa=self._pressure,
                forcing_timescale_s=276.0,
            )
            for _ in cells
        )
        sequence = FilmDepletionSequence(
            films=films,
            interfacial_area_m2=areas,
            gas_side_area_m2=gas_areas,
            bed_contact_h_w_m2_k=(self._contact_h,) * len(cells),
            gas_heat_capacity_j_kg_k=REFERENCE_CP_GAS_J_KG_K,
            layer_holdup_kg=holdups,
            layer_duty_w=duties,
            macro_steps_per_residence=self._steps,
            max_residences=self._residences,
            energy_screening_tier=self._screening_tier,
        )
        result = sequence.run(feed=feed, gas_inlet=contact_steam)
        final = result.steps[-1]

        # --- whole-machine boundary streams -------------------------
        released_h = result.residence_avg_hexane_release_kg_s
        released_w = result.residence_avg_water_release_kg_s
        discharge_x_h = result.discharge_x_h
        discharge_x_w = final.layer_loadings_x_w[0]
        dome = GasStream(
            G_w=sparge_steam.G_w + released_w,
            G_h=sparge_steam.G_h + released_h,
            T=final.films[-1].interface_temperature_k,
        )

        # --- ledgers, per species, across the whole machine ----------
        fed_h = feed.F_dm * feed.X_h
        out_h = feed.F_dm * discharge_x_h + (dome.G_h - sparge_steam.G_h)
        fed_w = feed.F_dm * feed.X_w + sparge_steam.G_w
        out_w = feed.F_dm * discharge_x_w + dome.G_w
        h_res = abs(fed_h - out_h) / max(fed_h, 1.0e-12)
        w_res = abs(fed_w - out_w) / max(fed_w, 1.0e-12)

        # --- machine energy ledger ----------------------------------
        # duty in + steam enthalpy in = solid sensible out + latent of
        # the net release + dome sensible, referenced to the feed state.
        # Machine energy ledger in the form of the QG-6E brief section
        # 4.3 (docs/DTDC_QG6E_QG9E_engineering_assumptions_brief):
        #     Q_ind = H_meal_out + H_vap_out - H_cake_in - m_DS h_DS + Q_loss
        # The direct-steam enthalpy term is the one the first T3a draft
        # omitted, which is why its residual came out at 2.94 relative.
        # Datum: liquid water, liquid hexane and dry solid all at the
        # FEED temperature, so H_cake_in = 0 by construction.
        # Constant-property basis, disclosed: this is an engineering
        # closure check, not the enthalpy-exact node (that arrives with
        # the dynamic wall).  Q_loss is the brief's A-05 term, nominal
        # 1 % of gross incoming thermal duty.
        CP_W_LIQ, CP_W_VAP, LAT_W = 4186.0, 1996.0, 2.34e6
        CP_H_LIQ, LAT_H = 2260.0, 3.35e5
        cp_solid = TrayParams().cp_dry_meal
        # PHYSICAL duty only.  `duties` carries the sparge desuperheat
        # rerouted as tray duty (delivery-path bookkeeping for the
        # sequence), but the LEDGER's steam term below already counts
        # the steam at its full inlet temperature — summing `duties`
        # here double-counted the superheat (~405 kW at reference,
        # measured as part of a 14.4 % machine energy residual).
        duty_total = sum(t.indirect_duty_w for t in cells)
        dT_out = result.discharge_temperature_k - feed.T
        dT_dome = dome.T - feed.T

        steam_mass_in = sparge_steam.G_w + sparge_steam.G_h
        h_steam_in = steam_mass_in * (
            LAT_W + CP_W_VAP * (sparge_steam.T - feed.T)
        )
        q_loss = A05_LOSS_FRACTION * (h_steam_in + duty_total)

        h_meal_out = (
            feed.F_dm * cp_solid * dT_out
            + feed.F_dm * discharge_x_w * CP_W_LIQ * dT_out
            + feed.F_dm * discharge_x_h * CP_H_LIQ * dT_out
        )
        h_vap_out = (
            dome.G_w * (LAT_W + CP_W_VAP * dT_dome)
            + dome.G_h * (LAT_H + CP_H_LIQ * dT_dome)
        )
        energy_residual = (
            duty_total + h_steam_in - h_meal_out - h_vap_out - q_loss
        )
        energy_rel = abs(energy_residual) / max(
            duty_total + h_steam_in, 1.0
        )

        max_layer_energy = max(
            abs(f.energy_residual_w)
            for r in result.steps[-self._steps:]
            for f in r.films
        )
        return DTDCStackResult(
            tray_ids=tuple(t.tray_id for t in cells),
            exchange_area_m2=areas,
            holdup_kg=holdups,
            duty_w=duties,
            feed=feed,
            discharge_x_h=discharge_x_h,
            discharge_x_w=discharge_x_w,
            discharge_temperature_k=result.discharge_temperature_k,
            dome_vapor=dome,
            sparge_steam=sparge_steam,
            hexane_ledger_residual_rel=h_res,
            water_ledger_residual_rel=w_res,
            energy_ledger_residual_w=energy_residual,
            energy_ledger_residual_rel=energy_rel,
            max_layer_energy_residual_w=max_layer_energy,
            converged=result.converged,
            macro_steps=len(result.steps),
            energy_screening_tier=self._screening_tier,
        )
