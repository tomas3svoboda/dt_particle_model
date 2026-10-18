"""Couple three accepted zero-contact PD tray steps to the QSC bypass network.

The tray transactions already own the local hexane-vapour generation and the
wall-to-bypass-gas energy over their common accepted interval.  This adapter
turns those ledger quantities into rates and applies them once, in the causal
``PD3 -> PD2 -> PD1 -> DOME`` order, to the negligible-pressure-drop oracle.

The absolute common-datum enthalpy of the incoming PD3 gas is deliberately a
caller boundary.  A :class:`PDPlenumLedger` stores temperature, composition,
flow and a balance residual, but not the component-datum species offsets needed
to reconstruct that absolute enthalpy.  This module therefore validates the
supplied boundary's pressure, component rates and datum, and derives only the
three relative enthalpy increments from accepted ledgers.  It does not infer a
caloric offset or claim a complete tower coupling.

Q-CS1-R1 (2026-09-05).  "Zero contact" names the bypass: the gas does not flow
THROUGH the bed.  It never meant that the stream cannot touch the bed's top
surface, and CS-1's contact-limited condensation is exactly that touch, carried
on the accepted ledgers like every other local term.  This adapter therefore
admits a packet external-water change ONLY when, and only up to what, that
packet's own condensation ledger accounts for; a packet that received no
condensate keeps the verbatim zero-contact refusal.  The condensate is a sink
on the bypass stream.  The zero-drop oracle is a byte-pinned authority that the
whole certified QSC10 K4 chain fresh-verifies, so it is NOT edited: it keeps
projecting the sources alone, and this adapter subtracts each tray's sink from
that sink-free projection when it audits the accepted plenum rates and when it
states the dome outlet.  ``pressure_network_result`` is therefore the oracle's
sink-free projection and ``dome_outlet`` is the physical one; the two agree
exactly whenever nothing condensed.  Every sink term is exactly 0.0 at the
declared default contact fraction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar

from . import pd_k_dome_bypass_tray_integration as pd
from . import qsc_pd_bypass_pressure_network as pressure_network


PD_CAUSAL_ORDER = pressure_network.PD_BYPASS_TRAY_ORDER


class QSCPDZeroContactNetworkCouplingError(ValueError):
    """The accepted PD transactions cannot form the declared zero-contact chain."""


def _component_rates(stream: pd.PDBinaryGasStream) -> tuple[float, float]:
    """Return water and hexane rates in the pressure-network field order."""

    return (
        stream.molar_flow_mol_s * stream.water_mole_fraction,
        stream.molar_flow_mol_s * stream.hexane_mole_fraction,
    )


def _roundoff_allowance(*values: float) -> float:
    return math.fsum(math.ulp(value) for value in values)


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDCondensateSink:
    """Q-CS1-R1: what one tray's cold top surface took OUT of the bypass stream.

    Rates over the common accepted interval, read from the same plenum ledger
    the local sources come from.  Exactly zero on every tray that did not
    condense, which is every tray at the declared default contact fraction.
    """

    tray_id: str
    condensed_water_molar_flow_mol_s: float
    condensed_hexane_molar_flow_mol_s: float
    condensation_common_enthalpy_rate_w: float

    def __post_init__(self) -> None:
        for label, value in (
            ("condensed water molar flow", self.condensed_water_molar_flow_mol_s),
            ("condensed hexane molar flow", self.condensed_hexane_molar_flow_mol_s),
        ):
            if type(value) is not float or not math.isfinite(value) or value < 0.0:
                raise QSCPDZeroContactNetworkCouplingError(
                    f"{self.tray_id} {label} must be a finite nonnegative binary64 value"
                )
        if type(self.condensation_common_enthalpy_rate_w) is not float or not math.isfinite(
            self.condensation_common_enthalpy_rate_w
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                f"{self.tray_id} condensation common-enthalpy rate must be finite binary64"
            )

    @property
    def is_zero(self) -> bool:
        return (
            self.condensed_water_molar_flow_mol_s == 0.0
            and self.condensed_hexane_molar_flow_mol_s == 0.0
            and self.condensation_common_enthalpy_rate_w == 0.0
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDTransactionRateAudit:
    """Projection of one accepted finite-interval plenum ledger onto gas rates."""

    tray_id: str
    duration_s: float
    actual_inlet_component_rates_mol_s: tuple[float, float]
    network_inlet_component_rates_mol_s: tuple[float, float]
    actual_outlet_component_rates_mol_s: tuple[float, float]
    network_outlet_component_rates_mol_s: tuple[float, float]
    inlet_component_disagreement_bounds_mol_s: tuple[float, float]
    outlet_component_disagreement_bounds_mol_s: tuple[float, float]
    transaction_component_residual_rates_mol_s: tuple[float, float]
    generated_hexane_molar_flow_mol_s: float
    generated_hexane_common_enthalpy_rate_w: float
    wall_to_gas_heat_w: float
    network_common_enthalpy_increment_residual_w: float
    network_common_enthalpy_increment_roundoff_bound_w: float
    transaction_common_energy_residual_rate_w: float
    transaction_ledger_passed: bool
    #: Q-CS1-R1: the condensate sink this transaction projected, 0.0 unless a
    #: condensation ledger accounted for it.
    condensed_water_molar_flow_mol_s: float = 0.0
    condensed_hexane_molar_flow_mol_s: float = 0.0
    condensation_common_enthalpy_rate_w: float = 0.0

    absolute_pd_stream_common_enthalpy_reconstructed: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def component_projection_passed(self) -> bool:
        comparisons = zip(
            (
                *self.actual_inlet_component_rates_mol_s,
                *self.actual_outlet_component_rates_mol_s,
            ),
            (
                *self.network_inlet_component_rates_mol_s,
                *self.network_outlet_component_rates_mol_s,
            ),
            (
                *self.inlet_component_disagreement_bounds_mol_s,
                *self.outlet_component_disagreement_bounds_mol_s,
            ),
            strict=True,
        )
        return all(abs(actual - projected) <= bound for actual, projected, bound in comparisons)

    @property
    def passed(self) -> bool:
        return (
            self.component_projection_passed
            and abs(self.network_common_enthalpy_increment_residual_w)
            <= self.network_common_enthalpy_increment_roundoff_bound_w
            and self.transaction_ledger_passed
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCPDZeroContactNetworkCoupling:
    """Accepted-time ownership and pressure-network projection for three PD steps."""

    accepted_t_n_s: float
    endpoint_time_s: float
    duration_s: float
    transactions: tuple[pd.ValidatedPDKTrayStep, ...]
    pd3_inlet: pressure_network.QSCPDBinaryGasRate
    local_sources: tuple[pressure_network.QSCPDLocalGasSource, ...]
    transaction_rate_audits: tuple[QSCPDTransactionRateAudit, ...]
    #: The pinned oracle's SINK-FREE projection of the sources.
    pressure_network_result: pressure_network.QSCPDBypassPressureNetworkResult
    #: Q-CS1-R1: per-tray condensate sinks, in causal order, and the dome
    #: outlet with every sink subtracted.  ``dome_outlet`` equals
    #: ``pressure_network_result.dome_outlet`` exactly when no tray condensed.
    condensate_sinks: tuple[QSCPDCondensateSink, ...]
    dome_outlet: pressure_network.QSCPDBinaryGasRate

    zero_contact_limit: ClassVar[bool] = True
    #: Q-CS1-R1: a packet external-water change is admitted only when, and
    #: only up to what, its own condensation ledger accounts for.
    ledger_accounted_condensation_admitted: ClassVar[bool] = True
    caller_supplied_absolute_inlet_common_enthalpy: ClassVar[bool] = True
    absolute_pd_stream_common_enthalpy_reconstructed: ClassVar[bool] = False
    complete_contact_wired: ClassVar[bool] = False
    full_tower_coupling_implemented: ClassVar[bool] = False
    mech_qualified: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return all(audit.passed for audit in self.transaction_rate_audits)


def _validate_zero_contact_transactions(
    transactions: tuple[pd.ValidatedPDKTrayStep, ...],
) -> tuple[float, float, float, str]:
    if (
        type(transactions) is not tuple
        or len(transactions) != 3
        or any(type(item) is not pd.ValidatedPDKTrayStep for item in transactions)
    ):
        raise QSCPDZeroContactNetworkCouplingError(
            "coupling requires three exact evaluator-issued PD transactions"
        )
    tray_ids = tuple(item.request.prior.tray.physical_tray_id for item in transactions)
    if tray_ids != PD_CAUSAL_ORDER:
        raise QSCPDZeroContactNetworkCouplingError(
            "PD transactions must be supplied in causal PD3, PD2, PD1 order"
        )

    accepted_t_n_s = transactions[0].request.prior.tray.time_s
    endpoint_time_s = transactions[0].request.end_time_s
    duration_s = endpoint_time_s - accepted_t_n_s
    datum_id = transactions[0].request.prior.tray.energy_datum_id
    for transaction in transactions:
        pd.ValidatedPDKTrayStep.__post_init__(transaction)
        request = transaction.request
        ledger = transaction.plenum_ledger
        if (
            request.prior.tray.time_s != accepted_t_n_s
            or request.end_time_s != endpoint_time_s
            or request.duration_s != duration_s
            or transaction.candidate.tray.time_s != endpoint_time_s
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "all PD transactions must own one common accepted t_n interval"
            )
        if (
            transaction.host_transaction.request.prior is not request.prior.tray
            or transaction.host_transaction.candidate is not transaction.candidate.tray
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "PD transaction is not owned by its accepted t_n host state"
            )
        if (
            request.prior.tray.energy_datum_id != datum_id
            or transaction.candidate.tray.energy_datum_id != datum_id
            or request.model.properties.energy_datum_id != datum_id
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "PD transactions do not share one common-energy datum"
            )
        if (
            request.prescribed_bypass_pressure_drop_pa != 0.0
            or transaction.pressure_ledger.prescribed_pressure_drop_pa != 0.0
            or transaction.pressure_ledger.pressure_residual_pa != 0.0
            or ledger.inlet.pressure_pa != request.prior.tray.pressure_pa
            or ledger.outlet.pressure_pa != request.prior.tray.pressure_pa
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "current coupling admits only exact zero-drop fixed-pressure PD steps"
            )
        if not (
            all(item.passed for item in transaction.packet_ledgers)
            and transaction.inter_k_face.passed
            and transaction.pressure_ledger.passed
            and ledger.passed
            and transaction.wall_ledger.passed
            and transaction.system_ledger.passed
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "a PD transaction contains a failed physical ledger"
            )
        # Q-CS1-R1: a packet's external water may change ONLY by what its own
        # condensation ledger accounts for.  No condensate -> the verbatim
        # zero-contact refusal; condensate -> admitted up to that mass, within
        # the packet ledger's own binary64 roundoff bound, never beyond it.
        for item in transaction.packet_ledgers:
            if item.condensed_water_kg == 0.0:
                if (
                    item.external_water_before_kg != item.external_water_after_kg
                    or item.external_water_residual_kg != 0.0
                ):
                    raise QSCPDZeroContactNetworkCouplingError(
                        "zero-contact PD step changed a packet external-water inventory"
                    )
                continue
            bound = item.mass_binary64_roundoff_bound_kg
            change = math.fsum((item.external_water_after_kg, -item.external_water_before_kg))
            if (
                item.condensed_water_kg < 0.0
                or change > math.fsum((item.condensed_water_kg, bound))
                or abs(item.external_water_residual_kg) > bound
            ):
                raise QSCPDZeroContactNetworkCouplingError(
                    "PD step changed a packet external-water inventory beyond what its "
                    "condensation ledger accounts for"
                )
        system = transaction.system_ledger
        accounted_condensate_kg = math.fsum(
            item.condensed_water_kg for item in transaction.packet_ledgers
        )
        water_bound = system.water_binary64_roundoff_bound_kg
        if accounted_condensate_kg == 0.0:
            if (
                system.packet_water_change_kg != 0.0
                or system.cell_water_change_kg != 0.0
                or abs(system.external_water_to_system_kg) > water_bound
                or abs(system.water_mass_residual_kg) > water_bound
            ):
                raise QSCPDZeroContactNetworkCouplingError(
                    "zero-contact PD step contains a water source beyond binary64 roundoff"
                )
        elif (
            # The condensate lands on packets, never in the cell; the cell's
            # water total may move only by the binary64 noise of carrying the
            # packets' new inventory (measured -1.4e-14 kg on the first
            # condensing march step, against a 1.3e-13 kg bound).
            abs(system.cell_water_change_kg) > water_bound
            or abs(math.fsum((system.packet_water_change_kg, -system.external_water_to_system_kg)))
            > math.fsum((water_bound, water_bound))
            or system.packet_water_change_kg > math.fsum((accounted_condensate_kg, water_bound))
            or abs(system.water_mass_residual_kg) > water_bound
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "PD step moved water beyond what its condensation ledgers account for: "
                f"packet change {system.packet_water_change_kg!r} kg, cell change "
                f"{system.cell_water_change_kg!r} kg, boundary {system.external_water_to_system_kg!r} "
                f"kg, residual {system.water_mass_residual_kg!r} kg, packet-ledger condensate "
                f"{accounted_condensate_kg!r} kg, plenum condensate {ledger.condensed_water_kg!r} "
                f"kg, roundoff bound {water_bound!r} kg"
            )
        packet_generated_mass = math.fsum(
            item.generated_hexane_kg for item in transaction.packet_ledgers
        )
        packet_generated_energy = math.fsum(
            item.generated_hexane_common_energy_j for item in transaction.packet_ledgers
        )
        if (
            packet_generated_mass != ledger.generated_hexane_kg
            or packet_generated_energy != ledger.generated_hexane_common_energy_j
            or transaction.wall_ledger.wall_to_bypass_gas_energy_j
            != ledger.wall_to_bypass_gas_energy_j
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "PD local-source ledgers disagree; a source could be lost or counted twice"
            )
        # Q-CS1-R1: the condensate the packets RECEIVED must be the condensate
        # the plenum stream LOST, component by component, or the sink projected
        # below would not be the mass that landed.
        packet_condensed_water = math.fsum(
            item.condensed_water_kg for item in transaction.packet_ledgers
        )
        packet_condensed_hexane = math.fsum(
            item.condensed_hexane_kg for item in transaction.packet_ledgers
        )
        if abs(
            math.fsum((packet_condensed_water, -ledger.condensed_water_kg))
        ) > _roundoff_allowance(
            *(item.condensed_water_kg for item in transaction.packet_ledgers),
            ledger.condensed_water_kg,
        ) or abs(
            math.fsum((packet_condensed_hexane, -ledger.condensed_hexane_kg))
        ) > _roundoff_allowance(
            *(item.condensed_hexane_kg for item in transaction.packet_ledgers),
            ledger.condensed_hexane_kg,
        ):
            raise QSCPDZeroContactNetworkCouplingError(
                "PD condensation ledgers disagree; condensate could be lost or counted twice"
            )

    if any(
        lower.plenum_ledger.outlet != upper.plenum_ledger.inlet
        for lower, upper in zip(transactions, transactions[1:])
    ):
        raise QSCPDZeroContactNetworkCouplingError(
            "accepted PD gas states do not join exactly in causal order"
        )
    return accepted_t_n_s, endpoint_time_s, duration_s, datum_id


def couple_qsc_pd_zero_contact_network(
    transactions: tuple[pd.ValidatedPDKTrayStep, ...],
    *,
    pd3_inlet: pressure_network.QSCPDBinaryGasRate,
) -> QSCPDZeroContactNetworkCoupling:
    """Apply three accepted PD local-source ledgers once to the zero-drop oracle."""

    accepted_t_n_s, endpoint_time_s, duration_s, datum_id = _validate_zero_contact_transactions(
        transactions
    )
    if type(pd3_inlet) is not pressure_network.QSCPDBinaryGasRate:
        raise QSCPDZeroContactNetworkCouplingError(
            "PD3 inlet boundary must be an exact QSC binary-gas rate"
        )
    actual_pd3_inlet = transactions[0].plenum_ledger.inlet
    actual_water_rate, actual_hexane_rate = _component_rates(actual_pd3_inlet)
    if (
        pd3_inlet.pressure_pa != actual_pd3_inlet.pressure_pa
        or pd3_inlet.water_molar_flow_mol_s != actual_water_rate
        or pd3_inlet.hexane_molar_flow_mol_s != actual_hexane_rate
        or pd3_inlet.energy_datum_id != datum_id
    ):
        raise QSCPDZeroContactNetworkCouplingError(
            "caller PD3 boundary pressure, component rates, or datum differs from the "
            "validated PD3 plenum inlet"
        )

    local_sources: list[pressure_network.QSCPDLocalGasSource] = []
    condensate_sinks: list[QSCPDCondensateSink] = []
    for transaction in transactions:
        ledger = transaction.plenum_ledger
        tray_id = transaction.request.prior.tray.physical_tray_id
        local_sources.append(
            pressure_network.QSCPDLocalGasSource(
                tray_id=tray_id,
                generated_hexane_molar_flow_mol_s=(ledger.generated_hexane_mol / duration_s),
                generated_hexane_common_enthalpy_rate_w=(
                    ledger.generated_hexane_common_energy_j / duration_s
                ),
                wall_to_gas_heat_w=ledger.wall_to_bypass_gas_energy_j / duration_s,
                energy_datum_id=datum_id,
            )
        )
        # Q-CS1-R1: the condensate sink, from the same plenum ledger the
        # sources come from.  0.0 / duration is exactly 0.0.
        condensate_sinks.append(
            QSCPDCondensateSink(
                tray_id=tray_id,
                condensed_water_molar_flow_mol_s=ledger.condensed_water_mol / duration_s,
                condensed_hexane_molar_flow_mol_s=ledger.condensed_hexane_mol / duration_s,
                condensation_common_enthalpy_rate_w=(
                    ledger.condensation_common_energy_j / duration_s
                ),
            )
        )
    sources = tuple(local_sources)
    sinks = tuple(condensate_sinks)
    network_result = pressure_network.evaluate_qsc_pd_bypass_pressure_network(
        pressure_network.QSCPDBypassPressureNetworkRequest(
            mode=pressure_network.ZERO_DROP_COLETTO_ORACLE,
            dome_pressure_pa=transactions[-1].candidate.tray.pressure_pa,
            pd3_inlet=pd3_inlet,
            local_sources=sources,
        )
    )

    audits: list[QSCPDTransactionRateAudit] = []
    cumulative_hexane_residual = 0.0
    cumulative_water_residual = 0.0
    # Q-CS1-R1: the oracle projects the sources alone; every sink admitted
    # below a node is subtracted from the projection at that node here.
    upstream_water_sink = 0.0
    upstream_hexane_sink = 0.0
    for transaction, source, sink, edge in zip(
        transactions, sources, sinks, network_result.edges, strict=True
    ):
        ledger = transaction.plenum_ledger
        inlet_actual = _component_rates(ledger.inlet)
        outlet_actual = _component_rates(ledger.outlet)
        inlet_network = (
            math.fsum((edge.inlet.water_molar_flow_mol_s, -upstream_water_sink)),
            math.fsum((edge.inlet.hexane_molar_flow_mol_s, -upstream_hexane_sink)),
        )
        inclusive_water_sink = math.fsum(
            (upstream_water_sink, sink.condensed_water_molar_flow_mol_s)
        )
        inclusive_hexane_sink = math.fsum(
            (upstream_hexane_sink, sink.condensed_hexane_molar_flow_mol_s)
        )
        outlet_network = (
            math.fsum((edge.outlet.water_molar_flow_mol_s, -inclusive_water_sink)),
            math.fsum((edge.outlet.hexane_molar_flow_mol_s, -inclusive_hexane_sink)),
        )
        inlet_water_bound = math.fsum(
            (
                cumulative_water_residual,
                _roundoff_allowance(*inlet_actual, *inlet_network, upstream_water_sink),
            )
        )
        inlet_hexane_bound = math.fsum(
            (
                cumulative_hexane_residual,
                _roundoff_allowance(*inlet_actual, *inlet_network, upstream_hexane_sink),
            )
        )
        water_residual_rate = abs(ledger.water_molar_residual) / duration_s
        hexane_residual_rate = abs(ledger.hexane_molar_residual) / duration_s
        cumulative_water_residual = math.fsum((cumulative_water_residual, water_residual_rate))
        cumulative_hexane_residual = math.fsum((cumulative_hexane_residual, hexane_residual_rate))
        outlet_water_bound = math.fsum(
            (
                cumulative_water_residual,
                _roundoff_allowance(*outlet_actual, *outlet_network, inclusive_water_sink),
            )
        )
        outlet_hexane_bound = math.fsum(
            (
                cumulative_hexane_residual,
                _roundoff_allowance(*outlet_actual, *outlet_network, inclusive_hexane_sink),
            )
        )
        upstream_water_sink = inclusive_water_sink
        upstream_hexane_sink = inclusive_hexane_sink
        common_increment_residual = math.fsum(
            (
                edge.outlet.common_datum_enthalpy_rate_w,
                -edge.inlet.common_datum_enthalpy_rate_w,
                -source.generated_hexane_common_enthalpy_rate_w,
                -source.wall_to_gas_heat_w,
            )
        )
        common_increment_roundoff_bound = _roundoff_allowance(
            edge.outlet.common_datum_enthalpy_rate_w,
            edge.inlet.common_datum_enthalpy_rate_w,
            source.generated_hexane_common_enthalpy_rate_w,
            source.wall_to_gas_heat_w,
            common_increment_residual,
        )
        audit = QSCPDTransactionRateAudit(
            tray_id=source.tray_id,
            duration_s=duration_s,
            actual_inlet_component_rates_mol_s=inlet_actual,
            network_inlet_component_rates_mol_s=inlet_network,
            actual_outlet_component_rates_mol_s=outlet_actual,
            network_outlet_component_rates_mol_s=outlet_network,
            inlet_component_disagreement_bounds_mol_s=(
                inlet_water_bound,
                inlet_hexane_bound,
            ),
            outlet_component_disagreement_bounds_mol_s=(
                outlet_water_bound,
                outlet_hexane_bound,
            ),
            transaction_component_residual_rates_mol_s=(
                ledger.water_molar_residual / duration_s,
                ledger.hexane_molar_residual / duration_s,
            ),
            generated_hexane_molar_flow_mol_s=(source.generated_hexane_molar_flow_mol_s),
            generated_hexane_common_enthalpy_rate_w=(
                source.generated_hexane_common_enthalpy_rate_w
            ),
            wall_to_gas_heat_w=source.wall_to_gas_heat_w,
            condensed_water_molar_flow_mol_s=sink.condensed_water_molar_flow_mol_s,
            condensed_hexane_molar_flow_mol_s=sink.condensed_hexane_molar_flow_mol_s,
            condensation_common_enthalpy_rate_w=sink.condensation_common_enthalpy_rate_w,
            network_common_enthalpy_increment_residual_w=common_increment_residual,
            network_common_enthalpy_increment_roundoff_bound_w=(common_increment_roundoff_bound),
            transaction_common_energy_residual_rate_w=(
                ledger.common_energy_residual_j / duration_s
            ),
            transaction_ledger_passed=ledger.passed,
        )
        if not audit.passed:
            raise QSCPDZeroContactNetworkCouplingError(
                f"{source.tray_id} transaction does not project onto the serial rate network"
            )
        audits.append(audit)

    # Q-CS1-R1: the physical dome outlet is the oracle's sink-free dome outlet
    # less every sink admitted below it; a sink larger than the stream it
    # leaves is refused typed by the gas-rate's own nonnegativity, and by the
    # audit above before that.  Exactly the oracle's own object when nothing
    # condensed: fsum with zero terms is the identity.
    oracle_dome = network_result.dome_outlet
    dome_outlet = (
        oracle_dome
        if all(sink.is_zero for sink in sinks)
        else pressure_network.QSCPDBinaryGasRate(
            pressure_pa=oracle_dome.pressure_pa,
            water_molar_flow_mol_s=math.fsum(
                (oracle_dome.water_molar_flow_mol_s, -upstream_water_sink)
            ),
            hexane_molar_flow_mol_s=math.fsum(
                (oracle_dome.hexane_molar_flow_mol_s, -upstream_hexane_sink)
            ),
            common_datum_enthalpy_rate_w=math.fsum(
                (
                    oracle_dome.common_datum_enthalpy_rate_w,
                    *(-sink.condensation_common_enthalpy_rate_w for sink in sinks),
                )
            ),
            energy_datum_id=oracle_dome.energy_datum_id,
        )
    )
    result = QSCPDZeroContactNetworkCoupling(
        accepted_t_n_s=accepted_t_n_s,
        endpoint_time_s=endpoint_time_s,
        duration_s=duration_s,
        transactions=transactions,
        pd3_inlet=pd3_inlet,
        local_sources=sources,
        transaction_rate_audits=tuple(audits),
        pressure_network_result=network_result,
        condensate_sinks=sinks,
        dome_outlet=dome_outlet,
    )
    if not result.passed:
        raise QSCPDZeroContactNetworkCouplingError(
            "zero-contact PD coupling contains a failed transaction-rate audit"
        )
    return result


__all__ = (
    "PD_CAUSAL_ORDER",
    "QSCPDCondensateSink",
    "QSCPDTransactionRateAudit",
    "QSCPDZeroContactNetworkCoupling",
    "QSCPDZeroContactNetworkCouplingError",
    "couple_qsc_pd_zero_contact_network",
)
