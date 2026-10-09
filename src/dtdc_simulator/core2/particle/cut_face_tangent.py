r"""Zero-volume right-limit transition at an inner master face.

This module supplies the singular seam between an accepted exact-face arrival
and the ordinary positive-volume face-departure chart.  It does not create a
small dry cell.  Instead, it evaluates the dry-side interface flux from
one-sided gradient primitives at exactly zero newborn volume.

The locally equilibrated dry phase has one composition degree of freedom.
Consequently the gas-exchange and retained-water thermodynamic-force
gradients are not independent unknowns.  The analytic fixed-``(T,P)``
composition derivative of the authoritative equilibrium state maps the
gas-exchange-force gradient to the retained-water-force gradient.  On the
smooth Luikov branch this map has a positive finite ratio; on the retained-
water cap branch the retained-force derivative and retained flux are exactly
zero.  The latter is an active constitutive branch, not a null mode.  Summed
component continuity then recovers the interface total Stefan flux
algebraically,

``N_t,Gamma = (N_w + N_h)_outer - N_w,ret,Gamma``.

All three component/energy flux-continuity rows are reported.  After the
summed recovery, the water row is algebraically dependent on the hexane row;
the square residual therefore consists of three Rankine--Hugoniot rows plus
the independent hexane and energy continuity rows.  There are five unknowns:
``(T_Gamma, s_dot, grad(T), grad(phi_wh), N_t,outer)``.

The outer dry flux is the zero-newborn-thickness limit of an ordinary
internal dry face.  It therefore uses the same symmetric endpoint-temperature
thermodynamic-force reconstruction as that face.  The moving-interface flux
remains a distinct one-sided ``T_Gamma`` trace.  Keeping those roles separate
makes the tangent and strict-cut operators continuous at the topology handoff
without adding a Soret term or changing a governing balance.

``continue_to_positive_sweep`` maps a certified tangent to an ordinary
strict-cut departure candidate at a finite positive swept volume.  It is a
predictor/overlap audit, not an accepted time step: the ordinary departure
nonlinear solve still owns acceptance and cumulative-ledger commitment.
There is no epsilon cell, clipping, pinning, tolerance relaxation, or contact
fallback anywhere in this module.  Every result remains a numerical oracle
with ``physically_qualifying=False``.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np
from scipy import optimize

from dtdc_simulator.core2.particle import actual_composition_force_path as acfp
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_face_departure as departure
from dtdc_simulator.core2.particle import cut_face_event_integrator as face_solver
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx


class FaceTangentTopologyError(ValueError):
    """The proposed right-limit state is outside the exact-face chart."""


class FaceTangentStepError(RuntimeError):
    """Rejected tangent/continuation carrying the exact immutable rollback."""

    def __init__(
        self,
        message: str,
        rollback_state: face_solver.FaceArrivalCommittedState,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


class FaceTangentNoAdmissibleRootError(FaceTangentStepError):
    """The exact-face equations have no resolved ``q_h >= 0`` root.

    This is deliberately distinct from iteration exhaustion or an
    inadmissible nonlinear trial.  It reports a converged bounded
    reconstruction whose residual cannot close on the frozen
    monotonic-primary-drainage branch, and it carries the exact immutable
    input as its rollback token.
    """


@dataclass(frozen=True)
class FaceTangentLayout:
    """Five-variable rank after exact summed-component elimination."""

    unknown_blocks: tuple[tuple[str, int], ...] = (
        ("interface_temperature", 1),
        ("front_speed", 1),
        ("dry_temperature_gradient", 1),
        ("dry_gas_force_gradient", 1),
        ("outer_total_stefan_flux", 1),
    )
    independent_residual_blocks: tuple[tuple[str, int], ...] = (
        ("rh_water", 1),
        ("rh_hexane", 1),
        ("rh_energy", 1),
        ("dry_hexane_flux_continuity", 1),
        ("dry_energy_flux_continuity", 1),
    )
    diagnostic_residual_blocks: tuple[tuple[str, int], ...] = (
        ("rh_water", 1),
        ("rh_hexane", 1),
        ("rh_energy", 1),
        ("dry_water_flux_continuity", 1),
        ("dry_hexane_flux_continuity", 1),
        ("dry_energy_flux_continuity", 1),
    )
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def unknown_count(self) -> int:
        return sum(size for _, size in self.unknown_blocks)

    @property
    def independent_residual_count(self) -> int:
        return sum(size for _, size in self.independent_residual_blocks)

    @property
    def diagnostic_residual_count(self) -> int:
        return sum(size for _, size in self.diagnostic_residual_blocks)

    @property
    def is_square(self) -> bool:
        return self.unknown_count == self.independent_residual_count == 5


@dataclass(frozen=True)
class FaceTangentUnknowns:
    """Right-limit algebraic variables in the declared rank ordering."""

    interface_temperature_k: float
    front_speed_m_s: float
    dry_temperature_gradient_k_m: float
    dry_gas_force_gradient_m_inv: float
    outer_total_stefan_flux_mol_m2_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def vector(self) -> tuple[float, ...]:
        return (
            self.interface_temperature_k,
            self.front_speed_m_s,
            self.dry_temperature_gradient_k_m,
            self.dry_gas_force_gradient_m_inv,
            self.outer_total_stefan_flux_mol_m2_s,
        )

    @classmethod
    def from_vector(cls, values: Sequence[float]) -> "FaceTangentUnknowns":
        vector = tuple(values)
        if len(vector) != 5:
            raise ValueError(f"tangent vector has length {len(vector)}, expected 5")
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("tangent candidate must be finite")
        return cls(*vector)


@dataclass(frozen=True)
class FaceTangentForceMap:
    """Analytic active-set map between the two dry force gradients."""

    gas_force_derivative_per_y: float
    retained_force_derivative_per_y: float
    retained_per_gas_gradient_ratio: float
    retained_water_active_set: cp.RetainedWaterActiveSet
    derivative_source: str = field(
        default="analytic frozen-virial activity derivative",
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.gas_force_derivative_per_y,
            self.retained_force_derivative_per_y,
            self.retained_per_gas_gradient_ratio,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("face-tangent force map must be finite")
        if self.gas_force_derivative_per_y >= 0.0:
            raise ValueError("gas force must decrease with y_hexane")
        if not isinstance(self.retained_water_active_set, cp.RetainedWaterActiveSet):
            raise TypeError("face-tangent force map requires a retained-water active set")
        if self.retained_water_active_set is cp.RetainedWaterActiveSet.LUIKOV_SMOOTH:
            if self.retained_force_derivative_per_y >= 0.0:
                raise ValueError("smooth retained-water force must decrease with y_hexane")
            if self.retained_per_gas_gradient_ratio <= 0.0:
                raise ValueError("smooth force-gradient ratio must be positive")
        elif self.retained_water_active_set is cp.RetainedWaterActiveSet.RETAINED_CAP:
            if self.retained_force_derivative_per_y != 0.0:
                raise ValueError("retained-cap force derivative must be exactly zero")
            if self.retained_per_gas_gradient_ratio != 0.0:
                raise ValueError("retained-cap force-gradient ratio must be exactly zero")
        else:  # pragma: no cover - exhaustive guard against a future active set
            raise ValueError("unrecognized retained-water active set")


@dataclass(frozen=True)
class FaceTangentDirectFlux:
    """Dry interface flux evaluated directly from non-singular gradients."""

    temperature_gradient_k_m: float
    gas_force_gradient_m_inv: float
    retained_force_gradient_m_inv: float
    binary_mobility_mol_m_s: float
    retained_water_mobility_mol_m_s: float
    independent_water_flux_mol_m2_s: float
    retained_water_flux_mol_m2_s: float
    conductive_heat_flux_w_m2: float
    recovered_total_stefan_flux_mol_m2_s: float
    component: cp.ComponentMolarFluxes
    energy: cp.ComponentEnergyFlux
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceTangentResidualBlocks:
    """Three RH rows and all three zero-volume dry continuity diagnostics."""

    rh_water_mol_s: float
    rh_hexane_mol_s: float
    rh_energy_w: float
    dry_water_flux_continuity_mol_s: float
    dry_hexane_flux_continuity_mol_s: float
    dry_energy_flux_continuity_w: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def independent_vector(self) -> tuple[float, ...]:
        return (
            self.rh_water_mol_s,
            self.rh_hexane_mol_s,
            self.rh_energy_w,
            self.dry_hexane_flux_continuity_mol_s,
            self.dry_energy_flux_continuity_w,
        )

    @property
    def diagnostic_vector(self) -> tuple[float, ...]:
        return (
            self.rh_water_mol_s,
            self.rh_hexane_mol_s,
            self.rh_energy_w,
            self.dry_water_flux_continuity_mol_s,
            self.dry_hexane_flux_continuity_mol_s,
            self.dry_energy_flux_continuity_w,
        )


@dataclass(frozen=True)
class FaceTangentLedger:
    """Rank, algebraic recovery, zero-volume, and rollback identities."""

    interface_area_m2: float
    newborn_old_volume_m3: float
    newborn_current_volume_m3: float
    newborn_inventory_change_mol_or_j: tuple[float, float, float]
    summed_component_recovery_error_mol_s: float
    dependent_water_row_identity_error_mol_s: float
    interface_total_minus_components_error_mol_m2_s: float
    front_hexane_jump_mol_m2_s: float
    used_coalescing_state_subtraction: bool
    total_stefan_flux_was_fitted: bool
    retained_force_was_independent_unknown: bool
    independent_rank_formula: str
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceTangentAssembly:
    """Complete exact-face right-limit residual without state commitment."""

    before: face_solver.FaceArrivalCommittedState
    candidate: FaceTangentUnknowns
    layout: FaceTangentLayout
    interface: cut.CutInterfaceState
    wet_cell: ww.WetWaterCellState
    wet_interface_flux: ww.WetWaterFaceFlux
    outer_dry_cell: cp.EquilibriumPoreState
    outer_dry_flux: cut.DryFaceFlux
    force_map: FaceTangentForceMap
    direct_interface_flux: FaceTangentDirectFlux
    residuals: FaceTangentResidualBlocks
    ledger: FaceTangentLedger
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square:
            raise RuntimeError("face tangent lost its square independent rank")
        if len(self.candidate.vector()) != self.layout.unknown_count:
            raise RuntimeError("face tangent candidate does not match its rank")
        if (
            len(self.residuals.independent_vector)
            != self.layout.independent_residual_count
        ):
            raise RuntimeError("face tangent independent residual lost rank")
        if (
            len(self.residuals.diagnostic_vector)
            != self.layout.diagnostic_residual_count
        ):
            raise RuntimeError("face tangent diagnostic residual lost a row")


@dataclass(frozen=True)
class FaceTangentOverlap:
    """Finite positive-sweep predictor evaluated by the ordinary chart."""

    before: face_solver.FaceArrivalCommittedState
    tangent: FaceTangentAssembly
    dt_s: float
    candidate: departure.FaceDepartureUnknowns
    ordinary_assembly: departure.FaceDepartureAssembly
    newborn_center_distance_m: float
    newborn_temperature_k: float
    newborn_y_hexane: float
    newborn_trace_gradients: departure.NewbornDryTraceGradients
    interface_water_flux_difference_mol_m2_s: float
    interface_hexane_flux_difference_mol_m2_s: float
    interface_energy_flux_difference_w_m2: float
    outer_water_flux_difference_mol_m2_s: float
    outer_hexane_flux_difference_mol_m2_s: float
    outer_energy_flux_difference_w_m2: float
    requires_ordinary_nonlinear_solve: bool = field(default=True, init=False)
    committed: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def maximum_flux_overlap_error(self) -> float:
        """Raw mixed-unit diagnostic; consumers should inspect typed fields."""

        return max(
            abs(value)
            for value in (
                self.interface_water_flux_difference_mol_m2_s,
                self.interface_hexane_flux_difference_mol_m2_s,
                self.interface_energy_flux_difference_w_m2,
                self.outer_water_flux_difference_mol_m2_s,
                self.outer_hexane_flux_difference_mol_m2_s,
                self.outer_energy_flux_difference_w_m2,
            )
        )


@dataclass(frozen=True)
class FaceTangentRankAudit:
    """Numerical local-rank certificate in independently scaled coordinates."""

    perturbations: tuple[float, ...]
    physical_jacobian: tuple[tuple[float, ...], ...]
    scaled_jacobian: tuple[tuple[float, ...], ...]
    scaled_singular_values: tuple[float, ...]
    numerical_rank: int
    expected_rank: int
    condition_proxy: float
    rank_tolerance: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def full_rank(self) -> bool:
        return self.numerical_rank == self.expected_rank


_LAYOUT = FaceTangentLayout()
_ROOT_RESIDUAL_TOLERANCE = 2.0e-11


def _stock_candidate_from_committed_face(
    before: face_solver.FaceArrivalCommittedState,
    bound: float,
    report: dict,
) -> FaceTangentUnknowns:
    """Construct the exact, flux-matched right-limit candidate (the stock solve).

    For fixed ``(T_Gamma, q_h)``, where ``q_h >= 0`` is the conserved
    outward n-hexane flux, the outer total Stefan flux and the remaining
    component/energy variables are recovered algebraically.  This coordinate
    makes the frozen monotonic-primary-drainage sign condition part of the
    nonlinear domain instead of rejecting an unconstrained root afterwards.
    The symmetric future-internal-face reconstruction differs from the
    moving-interface reconstruction on the committed arrival side, so the
    committed ``(T_Gamma, N_t)`` pair is only a seed.  A safeguarded two-row
    reconstruction closes the remaining water and energy RH equations.

    This reconstructs an exact zero-volume tangent seed; it is not an
    accepted finite continuation and commits no state or ledger.

    B-F1: ``bound`` is the residual bound in force (the literal
    ``_ROOT_RESIDUAL_TOLERANCE`` unless ``DTDC_FACE_TANGENT_RESIDUAL_BOUND``
    declares another); ``report`` receives report-only fields and never feeds
    back into the solve.
    """

    try:
        _validate_before(before)
        initial_data = _face_data(
            before,
            before.interface_temperature_k,
            before.dry_total_stefan_fluxes_mol_m2_s[0],
        )
        initial_hexane_flux = (
            initial_data[4].component.conserved_hexane_flux_mol_m2_s
        )
        hexane_flux_scale = max(
            abs(initial_hexane_flux),
            1.0e-6,
        )
        report["seed_hexane_flux_mol_m2_s"] = float(initial_hexane_flux)
        initial = _flux_matched_candidate_from_hexane_flux(
            before,
            before.interface_temperature_k,
            max(initial_hexane_flux, 0.0),
        )
        report["seed_front_speed_m_s"] = float(initial.front_speed_m_s)
        initial_assembly = _assemble_face_tangent(before, initial)
        water_scale, energy_scale = _rh_reconstruction_scales(initial_assembly)

        def remaining_rh(coordinates: np.ndarray) -> np.ndarray:
            candidate = _flux_matched_candidate_from_hexane_flux(
                before,
                float(coordinates[0]),
                float(coordinates[1]) * hexane_flux_scale,
            )
            residuals = _assemble_face_tangent(before, candidate).residuals
            return np.asarray(
                (
                    residuals.rh_water_mol_s / water_scale,
                    residuals.rh_energy_w / energy_scale,
                ),
                dtype=float,
            )

        temperature_bounds = (
            before.config.dry.conditioned_temperature_domain.solver_bounds_k
        )
        lower_temperature = math.nextafter(
            temperature_bounds[0],
            temperature_bounds[1],
        )
        upper_temperature = math.nextafter(
            temperature_bounds[1],
            temperature_bounds[0],
        )
        solution = optimize.least_squares(
            remaining_rh,
            np.asarray(
                (
                    before.interface_temperature_k,
                    max(initial_hexane_flux, 0.0) / hexane_flux_scale,
                ),
                dtype=float,
            ),
            method="trf",
            bounds=(
                np.asarray((lower_temperature, 0.0), dtype=float),
                np.asarray((upper_temperature, math.inf), dtype=float),
            ),
            ftol=1.0e-13,
            xtol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=100,
            x_scale="jac",
            diff_step=1.0e-5,
        )
        candidate = _flux_matched_candidate_from_hexane_flux(
            before,
            float(solution.x[0]),
            float(solution.x[1]) * hexane_flux_scale,
        )
        scaled_residuals = remaining_rh(solution.x)
        maximum_scaled_residual = float(np.max(np.abs(scaled_residuals)))
        report["stock_reconstruction_q_h_mol_m2_s"] = float(solution.x[1]) * hexane_flux_scale
        report["stock_two_row_scaled_residual"] = maximum_scaled_residual
        if not solution.success:
            raise RuntimeError(
                "face-tangent two-row reconstruction did not converge: "
                f"{solution.message}"
            )
        if solution.active_mask[0] != 0:
            raise FaceTangentTopologyError(
                "face-tangent two-row reconstruction touched a temperature guard"
            )
        if maximum_scaled_residual > bound:
            raise FaceTangentNoAdmissibleRootError(
                "face-tangent bounded reconstruction found no admissible "
                "primary-drainage root: "
                f"q_h={float(solution.x[1]) * hexane_flux_scale:.16e}, "
                f"maximum scaled residual={maximum_scaled_residual:.3e} > "
                f"{bound:.3e}",
                before,
            )
        _finish_root_report(
            report,
            before,
            candidate,
            hexane_flux=float(solution.x[1]) * hexane_flux_scale,
            two_row=maximum_scaled_residual,
        )
        return candidate
    except FaceTangentNoAdmissibleRootError:
        raise
    except FaceTangentStepError:
        raise
    except Exception as exc:
        if not isinstance(before, face_solver.FaceArrivalCommittedState):
            raise
        raise FaceTangentStepError(
            f"face-tangent seed rejected with exact rollback: {exc}",
            before,
        ) from exc


# ---------------------------------------------------------------------------
# B-F1 (2026-09-29): the inward-root seed, the declared residual bound and the
# per-call tangent report.  Build record (lead-placed):
# GT_PS2_FACE_TANGENT_INWARD_SEED_BUILD_RECORD_2026-09-29, measurement basis
# M-F1 of the same day.
#
# WHAT IT IS.  The stock bounded solve above seeds its two-row reconstruction
# at the committed arrival flux ``q_seed``.  M-F1 measured that every refused
# arrival state of the refinement ladder has an inward (``q_h > 0``,
# ``s_dot < 0``) root of the SAME exact-face rows beside an outward one, and
# that the stock trust-region walk reaches the inward root only when it lies
# within about 1.9 times ``q_seed``; beyond that it descends to the ``q_h = 0``
# bound and refuses.  When ``DTDC_FACE_TANGENT_INWARD_SEED=1`` and the stock
# solve refuses, a bracketing continuation in the flux ratio
# ``x = q_h / q_scale`` locates the inward root on the energy-row curve,
# polishes it, and hands it to the SAME stock reconstruction
# (``_flux_matched_candidate_from_hexane_flux``: ``q_h >= 0`` and
# ``s_dot <= 0`` raise), the SAME assembly (``_assemble_face_tangent``:
# ``s_dot > 0`` raises) and the SAME five-row check
# (``_maximum_scaled_root_residual`` against the residual bound in force).
# The downstream ``continue_to_positive_sweep`` repeats its own five-row check
# unchanged.  The outward root is never admitted: the scan runs over
# ``x > 0`` only and every admitted candidate passes the stock sign guards.
# A state with no admissible inward root re-raises the stock refusal object
# unchanged, text included.
#
# THE SEARCH MAP IS EVALUATION ONLY.  While scanning, the rows are evaluated
# without the sign guards (``_unguarded_*`` below, the M-F1 prototype's
# functions, whose rows and candidates M-F1 measured bit-equal to the stock
# assembly on every certified state and at 45 of 45 inward roots).  Nothing
# evaluated there is committed, returned or admitted; the prototype's exact
# grid, bracketing, polish and ordering are kept so that the knob reproduces
# the measured lane S roots bit for bit.
#
# THE RESIDUAL BOUND.  ``DTDC_FACE_TANGENT_RESIDUAL_BOUND`` replaces the
# literal ``_ROOT_RESIDUAL_TOLERANCE = 2.0e-11`` at every place this module
# compares against it (the stock two-row check, the seed's five-row check,
# ``continue_to_positive_sweep``).  Unset, the literal is in force and every
# comparison is the pre-B-F1 comparison.  Any other value is an owner-ruled
# declaration; the value in force is written on every tangent-call report,
# and a caller that records a job states ``residual_bound_declared``.
#
# THE FLOOR, AND WHY THERE IS NO EXTRA POLISH STAGE.  At the inward root the
# five-row residual is floored by the interface-composition bisection lattice
# (the 2026-08-22 tangent composition-floor record): the scaled water row
# jumps by about 1.7e-10 across one lattice step of y_Gamma and is smooth to
# about 2e-13 within a step, so the best point lies anywhere in [0, 8.5e-11]
# depending on where the root falls.  B-F1 step 2b measured four Newton
# polishes (relative steps 1e-7, 1e-9, 1e-11 and a regression Jacobian) on
# the 15 floor states: none moves the residual off the lattice value, so no
# polish stage is built; the owner-ruled bound is the declared answer there.
#
# Neither knob clips, projects, relaxes a guard or changes an equation.  Both
# are off/literal on disk.
# ---------------------------------------------------------------------------

FACE_TANGENT_INWARD_SEED_ENVIRONMENT_VARIABLE = "DTDC_FACE_TANGENT_INWARD_SEED"
FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE = "DTDC_FACE_TANGENT_RESIDUAL_BOUND"
#: The continuation grid in ``x = q_h / q_scale`` (``numpy.logspace`` arguments):
#: 1e-4 to 10**2.5 on 40 points, the M-F1 lane S grid, kept for bit identity
#: with the measured roots.  M-F1 found every inward root at 1.09 to 12.4 times
#: the committed flux and no sign change of the water row below it.
INWARD_SEED_FLUX_RATIO_GRID = (-4.0, 2.5, 40)
#: Bisection steps on the bracketed sign change (M-F1 lane S).
INWARD_SEED_BISECTION_STEPS = 80
#: Newton-polish iterations on the two-row reduction: UNDERIVED, the M-F1 prototype's
#: setting (2026-09-29), kept for bit identity with the measured lane S roots.
INWARD_SEED_NEWTON_ITERATIONS = 40

_FACE_TANGENT_CALL_REPORTS: list[dict] = []


def face_tangent_inward_seed_armed() -> bool:
    """Return whether ``DTDC_FACE_TANGENT_INWARD_SEED`` arms the inward seed.

    Unset or ``"0"`` is off; ``"1"`` is on; anything else is refused, so a
    mistyped knob never silently selects a path.
    """

    raw = os.environ.get(FACE_TANGENT_INWARD_SEED_ENVIRONMENT_VARIABLE)
    if raw is None or raw == "0":
        return False
    if raw == "1":
        return True
    raise ValueError(
        f"{FACE_TANGENT_INWARD_SEED_ENVIRONMENT_VARIABLE} must be '0' or '1', got {raw!r}"
    )


def face_tangent_residual_bound() -> float:
    """Return the root-residual bound in force (the literal when unset)."""

    raw = os.environ.get(FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE)
    if raw is None:
        return _ROOT_RESIDUAL_TOLERANCE
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(
            f"{FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE} is not a number: {raw!r}"
        ) from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(
            f"{FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE} must be finite and "
            f"positive, got {raw!r}"
        )
    return value


def face_tangent_residual_bound_declared() -> bool:
    """True when the bound in force is not the literal 2.0e-11."""

    return face_tangent_residual_bound() != _ROOT_RESIDUAL_TOLERANCE


def face_tangent_call_reports() -> tuple[dict, ...]:
    """Return copies of every tangent-call report recorded in this process."""

    return tuple(dict(report) for report in _FACE_TANGENT_CALL_REPORTS)


def drain_face_tangent_call_reports() -> tuple[dict, ...]:
    """Return and clear the tangent-call reports (report-only bookkeeping)."""

    reports = face_tangent_call_reports()
    _FACE_TANGENT_CALL_REPORTS.clear()
    return reports


def _plain_float(value: object) -> float | None:
    if value is None:
        return None
    return float(value)  # type: ignore[arg-type]


def _new_call_report(
    before: face_solver.FaceArrivalCommittedState,
    *,
    armed: bool,
    bound: float,
) -> dict:
    report: dict = {
        "call_index": len(_FACE_TANGENT_CALL_REPORTS),
        "inward_seed_armed": armed,
        "residual_bound": bound,
        "residual_bound_declared": bound != _ROOT_RESIDUAL_TOLERANCE,
        "path": None,
        "seed_hexane_flux_mol_m2_s": None,
        "seed_front_speed_m_s": None,
        "stock_refusal": None,
        "stock_reconstruction_q_h_mol_m2_s": None,
        "stock_two_row_scaled_residual": None,
        "root_hexane_flux_mol_m2_s": None,
        "root_flux_ratio_over_seed": None,
        "root_interface_temperature_k": None,
        "root_front_speed_m_s": None,
        "front_speed_ratio_root_over_seed": None,
        "root_five_row_scaled_residual": None,
        "root_two_row_scaled_residual_seed_scales": None,
        "seed_search_evaluations": None,
        "seed_search_candidates": None,
        "physically_qualifying": False,
    }
    try:
        report["time_s"] = float(before.time_s)
        report["arrival_face_index"] = int(before.arrival_face_index)
        report["cells"] = int(before.geometry.master_grid.n)
    except Exception:  # noqa: BLE001 - report only; the solve owns validation
        report["time_s"] = None
        report["arrival_face_index"] = None
        report["cells"] = None
    return report


def _finish_root_report(
    report: dict,
    before: face_solver.FaceArrivalCommittedState,
    candidate: FaceTangentUnknowns,
    *,
    hexane_flux: float,
    two_row: float | None,
    five_row: float | None = None,
) -> None:
    """Fill the root fields of a report (report-only; never gates)."""

    if five_row is None:
        try:
            five_row = _maximum_scaled_root_residual(_assemble_face_tangent(before, candidate))
        except Exception:  # noqa: BLE001 - report only
            five_row = None
    seed_flux = report.get("seed_hexane_flux_mol_m2_s")
    seed_speed = report.get("seed_front_speed_m_s")
    report["root_hexane_flux_mol_m2_s"] = float(hexane_flux)
    report["root_flux_ratio_over_seed"] = (
        float(hexane_flux) / seed_flux if seed_flux is not None and seed_flux > 0.0 else None
    )
    report["root_interface_temperature_k"] = float(candidate.interface_temperature_k)
    report["root_front_speed_m_s"] = float(candidate.front_speed_m_s)
    report["front_speed_ratio_root_over_seed"] = (
        float(candidate.front_speed_m_s) / seed_speed
        if seed_speed is not None and seed_speed != 0.0
        else None
    )
    report["root_five_row_scaled_residual"] = _plain_float(five_row)
    report["root_two_row_scaled_residual_seed_scales"] = _plain_float(two_row)


def candidate_from_committed_face(
    before: face_solver.FaceArrivalCommittedState,
) -> FaceTangentUnknowns:
    """Construct the exact, flux-matched right-limit candidate.

    The stock bounded reconstruction (``_stock_candidate_from_committed_face``)
    runs first and unchanged.  Only when it refuses and
    ``DTDC_FACE_TANGENT_INWARD_SEED=1`` does the inward-root seed run; it can
    only return a candidate that the stock reconstruction, the stock sign
    guards and the stock five-row check accept, and otherwise re-raises the
    stock refusal unchanged.  Every call appends one report-only row to the
    process's tangent-call reports.
    """

    armed = face_tangent_inward_seed_armed()
    bound = face_tangent_residual_bound()
    report = _new_call_report(before, armed=armed, bound=bound)
    _FACE_TANGENT_CALL_REPORTS.append(report)
    try:
        candidate = _stock_candidate_from_committed_face(before, bound, report)
    except FaceTangentStepError as stock_error:
        report["stock_refusal"] = str(stock_error)
        if not armed:
            report["path"] = "stock_refused_seed_off"
            raise
        seeded = _inward_root_seed(before, bound, report)
        if seeded is None:
            report["path"] = "inward_seed_no_admissible_root"
            raise
        report["path"] = "inward_seed_admitted"
        return seeded
    report["path"] = "stock_certified"
    return candidate


class _SeedMapUndefined(Exception):
    """A search-map evaluation left a closure's domain (evaluation only)."""


def _unguarded_outer_total_stefan(
    before: face_solver.FaceArrivalCommittedState,
    temperature_k: float,
    hexane_flux_mol_m2_s: float,
) -> float:
    """``_outer_total_stefan_from_hexane_flux`` without the ``q_h >= 0`` guard.

    Search-map evaluation only (B-F1); never admits a candidate.
    """

    interface, _wet, _wet_flux, outer_cell, zero_flux, _force_map, _mobility = _face_data(
        before, temperature_k, 0.0
    )
    diffusive = zero_flux.component.conserved_hexane_flux_mol_m2_s
    advective = math.fsum((hexane_flux_mol_m2_s, -diffusive))
    if advective > 0.0:
        upwind = interface.y_hexane
    elif advective < 0.0:
        upwind = outer_cell.y_hexane
    else:
        return 0.0
    if not math.isfinite(upwind) or upwind <= 0.0:
        raise _SeedMapUndefined("upwind y_hexane not positive")
    return advective / upwind


def _unguarded_flux_matched_candidate(
    before: face_solver.FaceArrivalCommittedState,
    temperature_k: float,
    outer_total_stefan_flux_mol_m2_s: float,
) -> tuple:
    """``_flux_matched_candidate`` algebra without the ``s_dot <= 0`` guard.

    Search-map evaluation only (B-F1); never admits a candidate.
    """

    data = _face_data(before, temperature_k, outer_total_stefan_flux_mol_m2_s)
    interface, _wet_cell, _wet_flux, _outer_cell, outer, force_map, mobility = data
    y_hexane = interface.y_hexane
    retained_mobility = before.config.dry.retained_water_mobility.value_mol_m_s
    total_outer = math.fsum(
        (
            outer.component.conserved_water_flux_mol_m2_s,
            outer.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    coefficient = math.fsum(
        (mobility, y_hexane * retained_mobility * force_map.retained_per_gas_gradient_ratio)
    )
    if not math.isfinite(coefficient) or coefficient <= 0.0:
        raise _SeedMapUndefined("force map cannot resolve continuity")
    gas_gradient = (
        outer.component.conserved_hexane_flux_mol_m2_s - y_hexane * total_outer
    ) / coefficient
    provisional = _direct_flux(
        before,
        interface,
        force_map,
        outer,
        temperature_gradient_k_m=0.0,
        gas_force_gradient_m_inv=gas_gradient,
    )
    conductivity = before.config.dry.thermal_conductivity.value_w_m_k
    temperature_gradient = (
        -(outer.energy.total_energy_flux_w_m2 - provisional.energy.total_energy_flux_w_m2)
        / conductivity
    )
    direct = _direct_flux(
        before,
        interface,
        force_map,
        outer,
        temperature_gradient_k_m=temperature_gradient,
        gas_force_gradient_m_inv=gas_gradient,
    )
    wet_hexane_concentration = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[before.arrival_face_index - 1]
        / hx.M
    )
    jump = math.fsum((wet_hexane_concentration, -interface.dry.total_hexane_concentration_mol_m3))
    if not math.isfinite(jump) or jump <= 0.0:
        raise _SeedMapUndefined("hexane capacity jump not positive")
    speed = -direct.component.conserved_hexane_flux_mol_m2_s / jump
    return (
        temperature_k,
        speed,
        temperature_gradient,
        gas_gradient,
        outer_total_stefan_flux_mol_m2_s,
    )


def _unguarded_rh_water_energy(
    before: face_solver.FaceArrivalCommittedState,
    vector: tuple,
) -> tuple[float, float, float, float]:
    """RH water and energy rows and their fixed physical scales at ``vector``.

    The same arithmetic as ``_assemble_face_tangent`` and
    ``_rh_reconstruction_scales`` without the ``s_dot <= 0`` guard.
    Search-map evaluation only (B-F1).
    """

    temperature, speed, temperature_gradient, gas_gradient, outer_total = vector
    interface, wet_cell, wet_flux, _outer_cell, outer_flux, force_map, _mobility = _face_data(
        before, temperature, outer_total
    )
    direct = _direct_flux(
        before,
        interface,
        force_map,
        outer_flux,
        temperature_gradient_k_m=temperature_gradient,
        gas_force_gradient_m_inv=gas_gradient,
    )
    area = before.geometry.master_grid.areas[before.arrival_face_index]
    rh_water = area * math.fsum(
        (
            direct.component.conserved_water_flux_mol_m2_s,
            -wet_flux.retained_water_flux_mol_m2_s,
            -speed * interface.dry.total_water_concentration_mol_m3,
            speed * wet_cell.retained_water_concentration_mol_m3,
        )
    )
    rh_energy = area * math.fsum(
        (
            direct.energy.total_energy_flux_w_m2,
            -wet_flux.total_energy_flux_w_m2,
            -speed * interface.dry.energy_density_j_m3,
            speed * wet_cell.energy_density_j_m3,
        )
    )
    water_scale = max(
        area * abs(direct.component.conserved_water_flux_mol_m2_s),
        area * abs(wet_flux.retained_water_flux_mol_m2_s),
        area
        * abs(speed)
        * max(
            interface.dry.total_water_concentration_mol_m3,
            interface.wet.retained_water_concentration_mol_m3,
        ),
        abs(rh_water),
        1.0e-300,
    )
    energy_scale = max(
        area * abs(direct.energy.total_energy_flux_w_m2),
        area * abs(wet_flux.total_energy_flux_w_m2),
        area
        * abs(speed)
        * max(
            abs(interface.dry.energy_density_j_m3),
            abs(interface.wet.energy_density_j_m3),
        ),
        abs(rh_energy),
        1.0e-300,
    )
    return rh_water, rh_energy, water_scale, energy_scale


class _InwardSeedReduction:
    """The stock two-row reduction in ``(T_Gamma, x = q_h / q_scale)``.

    Row scales are frozen at the committed seed exactly as the stock solve
    freezes them.  The water row is followed along the curve on which the
    energy row vanishes.  This is the M-F1 prototype's ``Reduced`` class,
    kept operation for operation so that the knob reproduces the measured
    roots bit for bit.  Evaluation only (B-F1).
    """

    def __init__(self, before: face_solver.FaceArrivalCommittedState) -> None:
        self.before = before
        self.t0 = before.interface_temperature_k
        data = _face_data(before, self.t0, before.dry_total_stefan_fluxes_mol_m2_s[0])
        self.q0 = data[4].component.conserved_hexane_flux_mol_m2_s
        self.q_scale = max(abs(self.q0), 1.0e-6)
        seed_vector = _unguarded_flux_matched_candidate(
            before,
            self.t0,
            _unguarded_outer_total_stefan(before, self.t0, max(self.q0, 0.0)),
        )
        _water, _energy, self.water_scale, self.energy_scale = _unguarded_rh_water_energy(
            before, seed_vector
        )
        lower, upper = before.config.dry.conditioned_temperature_domain.solver_bounds_k
        self.t_lower = math.nextafter(lower, upper)
        self.t_upper = math.nextafter(upper, lower)
        self.evaluations = 0

    def vector(self, temperature, ratio) -> tuple:
        return _unguarded_flux_matched_candidate(
            self.before,
            temperature,
            _unguarded_outer_total_stefan(self.before, temperature, ratio * self.q_scale),
        )

    def rows(self, temperature, ratio):
        self.evaluations += 1
        vector = self.vector(temperature, ratio)
        rh_water, rh_energy, _ws, _es = _unguarded_rh_water_energy(self.before, vector)
        return (
            np.array([rh_water / self.water_scale, rh_energy / self.energy_scale]),
            vector,
        )

    def safe_rows(self, temperature, ratio):
        try:
            return self.rows(temperature, ratio)
        except Exception:  # noqa: BLE001 - a map point outside the closures' domain
            return None, None

    def temperature_on_energy_row(self, ratio, count=33):
        """Scan, bracket nearest the committed ``T_Gamma``, then Brent."""

        temperatures = np.concatenate(
            [
                np.linspace(self.t_lower, self.t_upper, count),
                [self.t0],
                self.t0 + np.linspace(-8.0, 8.0, 161),
            ]
        )
        temperatures = np.unique(
            temperatures[(temperatures >= self.t_lower) & (temperatures <= self.t_upper)]
        )
        values = []
        for temperature in temperatures:
            rows, _ = self.safe_rows(temperature, ratio)
            values.append(np.nan if rows is None else rows[1])
        values = np.array(values)
        brackets = []
        for index in range(len(temperatures) - 1):
            left, right = values[index], values[index + 1]
            if np.isfinite(left) and np.isfinite(right) and left * right <= 0:
                brackets.append((temperatures[index], temperatures[index + 1]))
        if not brackets:
            return None, 0
        brackets.sort(key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - self.t0))
        left, right = brackets[0]
        temperature = optimize.brentq(
            lambda value: self.rows(value, ratio)[0][1],
            left,
            right,
            xtol=1e-13,
            rtol=4 * np.finfo(float).eps,
            maxiter=200,
        )
        return temperature, len(brackets)

    def temperature_on_energy_row_warm(self, ratio, guess):
        """Warm-started bracket around a nearby curve point; full scan otherwise."""

        if guess is not None:
            at_guess, _ = self.safe_rows(guess, ratio)
            if at_guess is not None:
                for half_width in (0.02, 0.1, 0.5, 2.0):
                    for left, right in (
                        (guess - half_width, guess),
                        (guess, guess + half_width),
                    ):
                        left, right = max(left, self.t_lower), min(right, self.t_upper)
                        at_left, _ = self.safe_rows(left, ratio)
                        at_right, _ = self.safe_rows(right, ratio)
                        if (
                            at_left is not None
                            and at_right is not None
                            and at_left[1] * at_right[1] <= 0
                        ):
                            temperature = optimize.brentq(
                                lambda value: self.rows(value, ratio)[0][1],
                                left,
                                right,
                                xtol=1e-13,
                                rtol=4 * np.finfo(float).eps,
                                maxiter=200,
                            )
                            return temperature, -1
        return self.temperature_on_energy_row(ratio)

    def water_on_curve(self, ratio, guess=None):
        try:
            temperature, brackets = self.temperature_on_energy_row_warm(ratio, guess)
        except Exception:  # noqa: BLE001 - a curve point outside the domain
            return np.nan, None, -1
        if temperature is None:
            return np.nan, None, 0
        rows, _vector = self.safe_rows(temperature, ratio)
        if rows is None:
            return np.nan, temperature, brackets
        return rows[0], temperature, brackets

    def newton_polish(self, temperature, ratio, iterations=INWARD_SEED_NEWTON_ITERATIONS):
        point = np.array([temperature, ratio], dtype=float)
        best = (np.inf, point.copy())
        for _ in range(iterations):
            rows, _vector = self.safe_rows(point[0], point[1])
            if rows is None:
                break
            norm = float(np.max(np.abs(rows)))
            if norm < best[0]:
                best = (norm, point.copy())
            if norm == 0.0:
                break
            jacobian = np.empty((2, 2))
            finite = True
            for column in range(2):
                step = max(abs(point[column]), 1e-3 if column == 0 else 1e-9) * 1e-7
                upper, lower = point.copy(), point.copy()
                upper[column] += step
                lower[column] -= step
                at_upper, _ = self.safe_rows(*upper)
                at_lower, _ = self.safe_rows(*lower)
                if at_upper is None or at_lower is None:
                    finite = False
                    break
                jacobian[:, column] = (at_upper - at_lower) / (2 * step)
            if not finite:
                break
            try:
                newton_step = np.linalg.solve(jacobian, -rows)
            except np.linalg.LinAlgError:
                break
            damping = 1.0
            improved = False
            for _search in range(40):
                trial = point + damping * newton_step
                at_trial, _ = self.safe_rows(trial[0], trial[1])
                if at_trial is not None and np.max(np.abs(at_trial)) < norm:
                    point = trial
                    improved = True
                    break
                damping *= 0.5
            if not improved:
                break
        return best


def _inward_seed_root_candidates(
    reduction: _InwardSeedReduction,
) -> list[tuple[str, float, float]]:
    """Bracket every water-row sign change on the energy curve at ``x > 0``.

    Returns ``(how, T_Gamma, x)`` polish results in the prototype's order:
    the Newton polish, then the stock least-squares settings re-seeded at the
    bracketed root.  Evaluation only.  Budgets (80 bisection steps, 40 Newton
    iterations, ``max_nfev=100`` for the re-seeded least squares, the stock
    solve's own value) are UNDERIVED: the M-F1 prototype's settings of
    2026-09-29, kept for bit identity with the measured lane S roots.
    """

    ratios = np.logspace(*INWARD_SEED_FLUX_RATIO_GRID)
    points = []
    guess = None
    for ratio in ratios:
        water, temperature, _brackets = reduction.water_on_curve(float(ratio), guess)
        if temperature is not None:
            guess = temperature
        points.append((float(ratio), water, temperature))
    found: list[tuple[str, float, float]] = []
    for (ratio_a, water_a, temperature_a), (ratio_b, water_b, _tb) in zip(points[:-1], points[1:]):
        if not (np.isfinite(water_a) and np.isfinite(water_b)) or water_a * water_b > 0:
            continue
        low, high, water_low = ratio_a, ratio_b, water_a
        for _ in range(INWARD_SEED_BISECTION_STEPS):
            middle = 0.5 * (low + high)
            if middle in (low, high):
                break
            water_middle, _t, _ = reduction.water_on_curve(middle, temperature_a)
            if not np.isfinite(water_middle):
                break
            if (water_middle < 0) == (water_low < 0):
                low, water_low = middle, water_middle
            else:
                high = middle
        ratio_root = 0.5 * (low + high)
        _water, temperature_root, _ = reduction.water_on_curve(ratio_root, temperature_a)
        _norm, polished = reduction.newton_polish(temperature_root, ratio_root)
        found.append(("newton", float(polished[0]), float(polished[1])))
        try:
            solution = optimize.least_squares(
                lambda c: reduction.rows(float(c[0]), float(c[1]))[0],
                np.array([temperature_root, ratio_root]),
                method="trf",
                bounds=(
                    np.array([reduction.t_lower, 0.0]),
                    np.array([reduction.t_upper, np.inf]),
                ),
                ftol=1e-13,
                xtol=1e-13,
                gtol=1e-13,
                max_nfev=100,
                x_scale="jac",
                diff_step=1e-5,
            )
            found.append(("stock_lsq_reseeded", float(solution.x[0]), float(solution.x[1])))
        except Exception:  # noqa: BLE001 - the re-seeded polish is optional
            pass
    return found


def _inward_root_seed(
    before: face_solver.FaceArrivalCommittedState,
    bound: float,
    report: dict,
) -> FaceTangentUnknowns | None:
    """Return the admitted inward root, or ``None`` (the caller re-raises).

    Admission, in order, for every polished candidate: the stock
    reconstruction ``_flux_matched_candidate_from_hexane_flux`` at
    ``(T_Gamma, q_h = x * q_scale)`` must succeed (``q_h >= 0``,
    ``s_dot <= 0``) and reproduce the searched vector bit for bit; ``q_h > 0``
    and ``s_dot < 0`` (inward, primary drainage); ``T_Gamma`` strictly inside
    the solver bounds; the stock assembly must succeed; and the stock five-row
    metric must be ``<= bound``.  Among admissible candidates the one with
    the smallest ``max(two-row, five-row)`` is returned (the prototype's
    order).  Nothing here clips, projects or relaxes.
    """

    candidates_report: list[dict] = []
    report["seed_search_candidates"] = candidates_report
    try:
        reduction = _InwardSeedReduction(before)
        found = _inward_seed_root_candidates(reduction)
    except Exception as exc:  # noqa: BLE001 - the stock refusal stands
        report["seed_search_error"] = f"{type(exc).__name__}: {exc}"
        return None
    admissible: list[tuple[float, int, FaceTangentUnknowns, float, float]] = []
    for order, (how, temperature, ratio) in enumerate(found):
        row: dict = {"how": how, "x": ratio, "interface_temperature_k": temperature}
        candidates_report.append(row)
        try:
            scaled, vector = reduction.rows(temperature, ratio)
            two_row = float(np.max(np.abs(scaled)))
            hexane_flux = float(ratio * reduction.q_scale)
            row.update(
                q_h_mol_m2_s=hexane_flux,
                front_speed_m_s=float(vector[1]),
                two_row_scaled_residual_seed_scales=two_row,
            )
            stock = _flux_matched_candidate_from_hexane_flux(
                before, temperature, ratio * reduction.q_scale
            )
            if tuple(stock.vector()) != tuple(float(value) for value in vector):
                row["rejected"] = "stock reconstruction differs from the searched vector"
                continue
            if not (hexane_flux > 0.0 and stock.front_speed_m_s < 0.0):
                row["rejected"] = "not an inward primary-drainage root"
                continue
            if not reduction.t_lower < temperature < reduction.t_upper:
                row["rejected"] = "interface temperature touches a solver bound"
                continue
            five_row = _maximum_scaled_root_residual(_assemble_face_tangent(before, stock))
            row["five_row_scaled_residual"] = float(five_row)
            if not five_row <= bound:
                row["rejected"] = f"five-row residual {five_row:.3e} > bound {bound:.3e}"
                continue
        except Exception as exc:  # noqa: BLE001 - an inadmissible candidate
            row["rejected"] = f"{type(exc).__name__}: {exc}"
            continue
        row["admissible"] = True
        admissible.append((max(two_row, float(five_row)), order, stock, two_row, five_row))
    report["seed_search_evaluations"] = reduction.evaluations
    if not admissible:
        return None
    admissible.sort(key=lambda item: (item[0], item[1]))
    _key, order, chosen, two_row, five_row = admissible[0]
    candidates_report[order]["selected"] = True
    report["seed_search_selected"] = candidates_report[order]["how"]
    _finish_root_report(
        report,
        before,
        chosen,
        hexane_flux=candidates_report[order]["q_h_mol_m2_s"],
        two_row=two_row,
        five_row=five_row,
    )
    return chosen


def _outer_total_stefan_from_hexane_flux(
    before: face_solver.FaceArrivalCommittedState,
    temperature_k: float,
    outward_hexane_flux_mol_m2_s: float,
) -> float:
    """Recover the outer-face ``N_t`` from one admissible hexane flux.

    The ordinary dry-face Maxwell--Stefan contribution is independent of
    ``N_t``.  Its advective contribution is piecewise linear because the
    existing face operator upwinds on the sign of ``N_t``.  Selecting the
    upwind mole fraction from the sign of the required advective remainder
    therefore gives the exact algebraic inverse without fitting another flux.
    """

    if not math.isfinite(outward_hexane_flux_mol_m2_s):
        raise ValueError("face-tangent hexane flux must be finite")
    if outward_hexane_flux_mol_m2_s < 0.0:
        raise FaceTangentTopologyError(
            "primary-drainage face tangent requires outward q_h >= 0"
        )
    zero_data = _face_data(before, temperature_k, 0.0)
    interface, _wet, _wet_flux, outer_cell, zero_flux, _force_map, _mobility = (
        zero_data
    )
    diffusive_hexane_flux = (
        zero_flux.component.conserved_hexane_flux_mol_m2_s
    )
    advective_hexane_flux = math.fsum(
        (outward_hexane_flux_mol_m2_s, -diffusive_hexane_flux)
    )
    if advective_hexane_flux > 0.0:
        upwind_y_hexane = interface.y_hexane
    elif advective_hexane_flux < 0.0:
        upwind_y_hexane = outer_cell.y_hexane
    else:
        return 0.0
    if not math.isfinite(upwind_y_hexane) or upwind_y_hexane <= 0.0:
        raise FaceTangentTopologyError(
            "outer-face algebraic Stefan recovery requires positive y_hexane"
        )
    return advective_hexane_flux / upwind_y_hexane


def _flux_matched_candidate_from_hexane_flux(
    before: face_solver.FaceArrivalCommittedState,
    temperature_k: float,
    outward_hexane_flux_mol_m2_s: float,
) -> FaceTangentUnknowns:
    """Close a tangent from admissible ``(T_Gamma, q_h >= 0)`` coordinates."""

    outer_total_stefan_flux = _outer_total_stefan_from_hexane_flux(
        before,
        temperature_k,
        outward_hexane_flux_mol_m2_s,
    )
    candidate = _flux_matched_candidate(
        before,
        temperature_k,
        outer_total_stefan_flux,
    )
    if candidate.front_speed_m_s > 0.0:
        raise RuntimeError(
            "admissible hexane-flux coordinate recovered reverse-core motion"
        )
    return candidate


def _flux_matched_candidate(
    before: face_solver.FaceArrivalCommittedState,
    temperature_k: float,
    outer_total_stefan_flux_mol_m2_s: float,
) -> FaceTangentUnknowns:
    """Algebraically close dry continuity and hexane RH for one seed pair."""

    data = _face_data(
        before,
        temperature_k,
        outer_total_stefan_flux_mol_m2_s,
    )
    interface, _wet_cell, _wet_flux, _outer_cell, outer, force_map, mobility = data
    y_hexane = interface.y_hexane
    retained_mobility = before.config.dry.retained_water_mobility.value_mol_m_s
    total_outer = math.fsum(
        (
            outer.component.conserved_water_flux_mol_m2_s,
            outer.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    coefficient = math.fsum(
        (
            mobility,
            y_hexane
            * retained_mobility
            * force_map.retained_per_gas_gradient_ratio,
        )
    )
    if not math.isfinite(coefficient) or coefficient <= 0.0:
        raise FaceTangentTopologyError(
            "local-equilibrium force map cannot resolve component continuity"
        )
    gas_gradient = (
        outer.component.conserved_hexane_flux_mol_m2_s
        - y_hexane * total_outer
    ) / coefficient
    provisional = _direct_flux(
        before,
        interface,
        force_map,
        outer,
        temperature_gradient_k_m=0.0,
        gas_force_gradient_m_inv=gas_gradient,
    )
    enthalpy_only = provisional.energy.total_energy_flux_w_m2
    conductivity = before.config.dry.thermal_conductivity.value_w_m_k
    temperature_gradient = -(
        outer.energy.total_energy_flux_w_m2 - enthalpy_only
    ) / conductivity
    direct = _direct_flux(
        before,
        interface,
        force_map,
        outer,
        temperature_gradient_k_m=temperature_gradient,
        gas_force_gradient_m_inv=gas_gradient,
    )
    wet_hexane_concentration = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[before.arrival_face_index - 1]
        / hx.M
    )
    jump = math.fsum(
        (
            wet_hexane_concentration,
            -interface.dry.total_hexane_concentration_mol_m3,
        )
    )
    if not math.isfinite(jump) or jump <= 0.0:
        raise FaceTangentTopologyError(
            "right-limit hexane capacity jump must remain positive"
        )
    speed = -direct.component.conserved_hexane_flux_mol_m2_s / jump
    if speed > 0.0:
        raise FaceTangentTopologyError(
            "right-limit candidate requests reverse-core motion"
        )
    return FaceTangentUnknowns(
        interface_temperature_k=temperature_k,
        front_speed_m_s=speed,
        dry_temperature_gradient_k_m=temperature_gradient,
        dry_gas_force_gradient_m_inv=gas_gradient,
        outer_total_stefan_flux_mol_m2_s=outer_total_stefan_flux_mol_m2_s,
    )


def _rh_reconstruction_scales(
    assembly: FaceTangentAssembly,
) -> tuple[float, float]:
    """Fixed physical scales for the two exact-tangent consistency rows."""

    area = assembly.ledger.interface_area_m2
    candidate = assembly.candidate
    interface = assembly.interface
    wet_flux = assembly.wet_interface_flux
    direct = assembly.direct_interface_flux
    water_scale = max(
        area * abs(direct.component.conserved_water_flux_mol_m2_s),
        area * abs(wet_flux.retained_water_flux_mol_m2_s),
        area
        * abs(candidate.front_speed_m_s)
        * max(
            interface.dry.total_water_concentration_mol_m3,
            interface.wet.retained_water_concentration_mol_m3,
        ),
        abs(assembly.residuals.rh_water_mol_s),
        1.0e-300,
    )
    energy_scale = max(
        area * abs(direct.energy.total_energy_flux_w_m2),
        area * abs(wet_flux.total_energy_flux_w_m2),
        area
        * abs(candidate.front_speed_m_s)
        * max(
            abs(interface.dry.energy_density_j_m3),
            abs(interface.wet.energy_density_j_m3),
        ),
        abs(assembly.residuals.rh_energy_w),
        1.0e-300,
    )
    return water_scale, energy_scale


def _maximum_scaled_root_residual(assembly: FaceTangentAssembly) -> float:
    """Return the five-row root residual in fixed physical scales."""

    water_scale, energy_scale = _rh_reconstruction_scales(assembly)
    area = assembly.ledger.interface_area_m2
    hexane_flux = (
        assembly.direct_interface_flux.component.conserved_hexane_flux_mol_m2_s
    )
    hexane_scale = max(
        area * abs(hexane_flux),
        abs(assembly.residuals.rh_hexane_mol_s),
        abs(assembly.residuals.dry_hexane_flux_continuity_mol_s),
        1.0e-300,
    )
    residuals = assembly.residuals
    return max(
        abs(residuals.rh_water_mol_s) / water_scale,
        abs(residuals.rh_hexane_mol_s) / hexane_scale,
        abs(residuals.rh_energy_w) / energy_scale,
        abs(residuals.dry_hexane_flux_continuity_mol_s) / hexane_scale,
        abs(residuals.dry_energy_flux_continuity_w) / energy_scale,
    )


def assemble_face_tangent(
    before: face_solver.FaceArrivalCommittedState,
    candidate: FaceTangentUnknowns,
) -> FaceTangentAssembly:
    """Assemble the exact zero-volume tangent residual or roll back exactly."""

    try:
        return _assemble_face_tangent(before, candidate)
    except FaceTangentStepError:
        raise
    except Exception as exc:
        if not isinstance(before, face_solver.FaceArrivalCommittedState):
            raise
        raise FaceTangentStepError(
            f"face-tangent residual rejected with exact rollback: {exc}",
            before,
        ) from exc


def audit_face_tangent_rank(
    assembly: FaceTangentAssembly,
    *,
    rank_tolerance: float = 1.0e-8,
) -> FaceTangentRankAudit:
    """Finite-difference and independently scale the exact 5-by-5 Jacobian.

    The perturbations are resolution probes, not nonlinear tolerances.  Every
    perturbed point is reassembled through the exact topology checks; an
    inadmissible probe fails with the original committed state as rollback.
    """

    before = assembly.before if isinstance(assembly, FaceTangentAssembly) else None
    try:
        if not isinstance(assembly, FaceTangentAssembly):
            raise TypeError("rank audit requires FaceTangentAssembly")
        if not math.isfinite(rank_tolerance) or rank_tolerance <= 0.0:
            raise ValueError("rank tolerance must be positive and finite")
        candidate = assembly.candidate
        values = np.asarray(candidate.vector(), dtype=float)
        perturbations = np.asarray(
            (
                1.0e-4,
                max(abs(candidate.front_speed_m_s) * 1.0e-4, 1.0e-12),
                max(abs(candidate.dry_temperature_gradient_k_m) * 1.0e-5, 1.0e-2),
                max(abs(candidate.dry_gas_force_gradient_m_inv) * 1.0e-5, 1.0e-2),
                max(
                    abs(candidate.outer_total_stefan_flux_mol_m2_s) * 1.0e-5,
                    1.0e-9,
                ),
            ),
            dtype=float,
        )
        jacobian = np.empty((5, 5), dtype=float)
        for column, delta in enumerate(perturbations):
            upper = values.copy()
            lower = values.copy()
            upper[column] += delta
            lower[column] -= delta
            upper_residual = assemble_face_tangent(
                assembly.before,
                FaceTangentUnknowns.from_vector(upper),
            ).residuals.independent_vector
            lower_residual = assemble_face_tangent(
                assembly.before,
                FaceTangentUnknowns.from_vector(lower),
            ).residuals.independent_vector
            jacobian[:, column] = (
                np.asarray(upper_residual, dtype=float)
                - np.asarray(lower_residual, dtype=float)
            ) / (2.0 * delta)
        row_scale = np.maximum(np.max(np.abs(jacobian), axis=1), 1.0e-300)
        scaled = jacobian / row_scale[:, None]
        column_scale = np.maximum(np.max(np.abs(scaled), axis=0), 1.0e-300)
        scaled = scaled / column_scale[None, :]
        singular_values = np.linalg.svd(scaled, compute_uv=False)
        numerical_rank = int(np.linalg.matrix_rank(scaled, tol=rank_tolerance))
        condition_proxy = float(singular_values[0] / singular_values[-1])
        return FaceTangentRankAudit(
            perturbations=tuple(float(value) for value in perturbations),
            physical_jacobian=tuple(
                tuple(float(value) for value in row) for row in jacobian
            ),
            scaled_jacobian=tuple(
                tuple(float(value) for value in row) for row in scaled
            ),
            scaled_singular_values=tuple(
                float(value) for value in singular_values
            ),
            numerical_rank=numerical_rank,
            expected_rank=_LAYOUT.unknown_count,
            condition_proxy=condition_proxy,
            rank_tolerance=rank_tolerance,
        )
    except FaceTangentStepError:
        raise
    except Exception as exc:
        if before is None:
            raise
        raise FaceTangentStepError(
            f"face-tangent rank audit rejected with exact rollback: {exc}",
            before,
        ) from exc


def _assemble_face_tangent(
    before: face_solver.FaceArrivalCommittedState,
    candidate: FaceTangentUnknowns,
) -> FaceTangentAssembly:
    _validate_before(before)
    if not isinstance(candidate, FaceTangentUnknowns):
        raise TypeError("candidate must be FaceTangentUnknowns")
    if not all(math.isfinite(value) for value in candidate.vector()):
        raise ValueError("tangent candidate must remain finite")
    if candidate.front_speed_m_s > 0.0:
        raise FaceTangentTopologyError(
            "primary-drainage face tangent requires s_dot <= 0"
        )
    data = _face_data(
        before,
        candidate.interface_temperature_k,
        candidate.outer_total_stefan_flux_mol_m2_s,
    )
    interface, wet_cell, wet_flux, outer_cell, outer_flux, force_map, _mobility = (
        data
    )
    direct = _direct_flux(
        before,
        interface,
        force_map,
        outer_flux,
        temperature_gradient_k_m=candidate.dry_temperature_gradient_k_m,
        gas_force_gradient_m_inv=candidate.dry_gas_force_gradient_m_inv,
    )
    grid = before.geometry.master_grid
    k = before.arrival_face_index
    area = grid.areas[k]
    speed = candidate.front_speed_m_s
    wet_hexane_concentration = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[k - 1]
        / hx.M
    )

    # O9a/O10a ALE FRONT-DONOR CORRECTION AND ITS ENERGY COMPLETION, ruled
    # 2026-08-21 in docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md
    # and docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md, applied
    # here as ONE COHERENT FAMILY TREATMENT because Stage 1b measured that
    # moving the mass donor alone breaks the water-datum gauge covariance of
    # the energy row - the two rows would convect different material states.
    #
    # PHYSICS GROUND, as recorded in the O10a ruling: two-phase moving-boundary
    # (Rankine-Hugoniot) jump conditions must convect mass and enthalpy OF THE
    # SAME MATERIAL STATE - an internal-consistency requirement of the balance
    # laws, not an empirical claim - and WHICH state it is, the bulk undrained
    # one, is the O9-measured fact (Peclet 3.4e3 forbids pre-drainage).  These
    # three per-area rows ARE those jump conditions.
    #
    # THE CONSUMED MATERIAL IS THIS CHART'S OWN.  At the exact tangent the
    # front sits on face ``k`` moving inward at ``s_dot <= 0``, so the material
    # it consumes is wet cell ``k - 1`` - which is precisely the ``wet_cell``
    # this module already evaluates from the committed state and already uses
    # as the left state of the front-face diffusive flux, and whose own
    # historical hexane loading the ``rh_hexane`` row below already convects.
    #
    # ROLE SEPARATION IS UNCHANGED: the zero-volume equilibrium trace keeps its
    # potential/flux role in full as the right state of ``wet_flux``.
    wet_water_donor_concentration = wet_cell.retained_water_concentration_mol_m3
    wet_energy_donor_density = wet_cell.energy_density_j_m3

    rh_water_per_area = math.fsum(
        (
            direct.component.conserved_water_flux_mol_m2_s,
            -wet_flux.retained_water_flux_mol_m2_s,
            -speed * interface.dry.total_water_concentration_mol_m3,
            # O9a, carried here by the O10a family treatment.
            speed * wet_water_donor_concentration,
        )
    )
    rh_hexane_per_area = math.fsum(
        (
            direct.component.conserved_hexane_flux_mol_m2_s,
            -speed * interface.dry.total_hexane_concentration_mol_m3,
            speed * wet_hexane_concentration,
        )
    )
    rh_energy_per_area = math.fsum(
        (
            direct.energy.total_energy_flux_w_m2,
            -wet_flux.total_energy_flux_w_m2,
            -speed * interface.dry.energy_density_j_m3,
            # O10a: the SAME material state as ``rh_water_per_area``'s donor
            # above, so the water caloric datum stays an exact gauge freedom of
            # this row.
            speed * wet_energy_donor_density,
        )
    )
    water_continuity_per_area = math.fsum(
        (
            direct.component.conserved_water_flux_mol_m2_s,
            -outer_flux.component.conserved_water_flux_mol_m2_s,
        )
    )
    hexane_continuity_per_area = math.fsum(
        (
            direct.component.conserved_hexane_flux_mol_m2_s,
            -outer_flux.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    energy_continuity_per_area = math.fsum(
        (
            direct.energy.total_energy_flux_w_m2,
            -outer_flux.energy.total_energy_flux_w_m2,
        )
    )
    residuals = FaceTangentResidualBlocks(
        rh_water_mol_s=area * rh_water_per_area,
        rh_hexane_mol_s=area * rh_hexane_per_area,
        rh_energy_w=area * rh_energy_per_area,
        dry_water_flux_continuity_mol_s=area * water_continuity_per_area,
        dry_hexane_flux_continuity_mol_s=area * hexane_continuity_per_area,
        dry_energy_flux_continuity_w=area * energy_continuity_per_area,
    )
    target_sum = math.fsum(
        (
            outer_flux.component.conserved_water_flux_mol_m2_s,
            outer_flux.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    direct_sum = math.fsum(
        (
            direct.component.conserved_water_flux_mol_m2_s,
            direct.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    recovery_error = area * math.fsum((direct_sum, -target_sum))
    dependent_identity = math.fsum(
        (
            residuals.dry_water_flux_continuity_mol_s,
            residuals.dry_hexane_flux_continuity_mol_s,
            -recovery_error,
        )
    )
    component_total_error = math.fsum(
        (
            direct.component.gas_water_flux_mol_m2_s,
            direct.component.gas_hexane_flux_mol_m2_s,
            -direct.recovered_total_stefan_flux_mol_m2_s,
        )
    )
    ledger = FaceTangentLedger(
        interface_area_m2=area,
        newborn_old_volume_m3=0.0,
        newborn_current_volume_m3=0.0,
        newborn_inventory_change_mol_or_j=(0.0, 0.0, 0.0),
        summed_component_recovery_error_mol_s=recovery_error,
        dependent_water_row_identity_error_mol_s=dependent_identity,
        interface_total_minus_components_error_mol_m2_s=component_total_error,
        front_hexane_jump_mol_m2_s=rh_hexane_per_area,
        used_coalescing_state_subtraction=False,
        total_stefan_flux_was_fitted=False,
        retained_force_was_independent_unknown=False,
        independent_rank_formula="3 RH + (hexane, energy) dry continuity = 5",
    )
    return FaceTangentAssembly(
        before=before,
        candidate=candidate,
        layout=_LAYOUT,
        interface=interface,
        wet_cell=wet_cell,
        wet_interface_flux=wet_flux,
        outer_dry_cell=outer_cell,
        outer_dry_flux=outer_flux,
        force_map=force_map,
        direct_interface_flux=direct,
        residuals=residuals,
        ledger=ledger,
    )


def continue_to_positive_sweep(
    tangent: FaceTangentAssembly,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
) -> FaceTangentOverlap:
    """Map the tangent to a finite strict-cut ordinary departure candidate.

    The newborn composition is the unique open-domain root whose discrete
    gas thermodynamic-force difference equals the tangent gradient.  No
    primitive is clipped or projected.  The returned ordinary assembly owns
    the exact ALE/GCL and material ledgers, but it is deliberately uncommitted
    and still requires the safeguarded ordinary nonlinear solve.
    """

    before = tangent.before if isinstance(tangent, FaceTangentAssembly) else None
    try:
        if not isinstance(tangent, FaceTangentAssembly):
            raise TypeError("continuation requires a FaceTangentAssembly")
        before = tangent.before
        if not isinstance(
            surface_boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("continuation requires a supported pore boundary")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("positive-sweep duration must be positive and finite")
        candidate_tangent = tangent.candidate
        outward_hexane_flux = (
            tangent.direct_interface_flux.component.conserved_hexane_flux_mol_m2_s
        )
        if outward_hexane_flux < 0.0:
            raise FaceTangentTopologyError(
                "positive-sweep continuation requires tangent q_h >= 0"
            )
        maximum_scaled_root_residual = _maximum_scaled_root_residual(tangent)
        residual_bound = face_tangent_residual_bound()
        if maximum_scaled_root_residual > residual_bound:
            raise FaceTangentTopologyError(
                "positive-sweep continuation requires an exact face-tangent root: "
                f"maximum scaled residual={maximum_scaled_root_residual:.3e} > "
                f"{residual_bound:.3e}"
            )
        face_radius = before.geometry.front.radius_m
        new_radius = math.fsum(
            (face_radius, candidate_tangent.front_speed_m_s * dt_s)
        )
        k = before.arrival_face_index
        inner_face = before.geometry.master_grid.faces[k - 1]
        if not inner_face < new_radius < face_radius:
            raise FaceTangentTopologyError(
                "finite tangent continuation must remain strictly inside cell k-1"
            )
        geometry = cg.partition_master_grid(
            before.geometry.master_grid,
            radius_m=new_radius,
        )
        # ``FaceDepartureUnknowns`` persists the front in volume coordinate
        # ``z`` and the ordinary residual reconstructs geometry from that
        # value.  Canonicalize through the same stored coordinate before
        # deriving either newborn trace.  Otherwise a radius -> z -> radius
        # round trip can change the center distance by one ulp and break the
        # exact ``trace = interface + gradient * distance`` identity even
        # though the represented physical front is unchanged.
        geometry = cg.partition_master_grid(
            before.geometry.master_grid,
            z=geometry.front.z,
        )
        if geometry.cut_cell_index != k - 1:
            raise FaceTangentTopologyError(
                "positive-sweep geometry did not enter the ordinary departure chart"
            )
        newborn_distance = math.fsum(
            (
                cut._dry_piece_center(  # noqa: SLF001
                    geometry.cells[k - 1],
                    geometry.front.radius_m,
                ),
                -geometry.front.radius_m,
            )
        )
        if newborn_distance <= 0.0:
            raise RuntimeError("strict departure lost its positive newborn distance")
        newborn_temperature = math.fsum(
            (
                candidate_tangent.interface_temperature_k,
                candidate_tangent.dry_temperature_gradient_k_m
                * newborn_distance,
            )
        )
        newborn_y, newborn_composition_gradient = (
            _newborn_composition_from_force_gradient(
            before,
            tangent.interface,
            newborn_temperature,
            newborn_distance,
            candidate_tangent.dry_gas_force_gradient_m_inv,
        )
        )
        trace_gradients = departure.NewbornDryTraceGradients(
            candidate_tangent.dry_temperature_gradient_k_m,
            newborn_composition_gradient,
        )
        direct = tangent.direct_interface_flux
        candidate = departure.FaceDepartureUnknowns(
            wet_temperatures_k=before.wet_temperatures_k,
            wet_retained_water_loadings=before.wet_retained_water_loadings,
            dry_temperatures_k=(newborn_temperature, *before.dry_temperatures_k),
            dry_y_hexane=(newborn_y, *before.dry_y_hexane),
            dry_total_stefan_fluxes_mol_m2_s=(
                direct.recovered_total_stefan_flux_mol_m2_s,
                candidate_tangent.outer_total_stefan_flux_mol_m2_s,
                *before.dry_total_stefan_fluxes_mol_m2_s[1:],
            ),
            front_z=geometry.front.z,
            interface_temperature_k=candidate_tangent.interface_temperature_k,
        )
        state = departure.adapt_face_arrival_committed_state(before)
        ordinary = departure.assemble_face_departure(
            state,
            candidate,
            dt_s,
            surface_boundary,
            trace_gradients,
            enforce_reduced_film_thresholds=False,
        )
        interface_flux = ordinary.dry_face_fluxes[0]
        outer_flux = ordinary.dry_face_fluxes[1]
        target_outer = tangent.outer_dry_flux
        return FaceTangentOverlap(
            before=before,
            tangent=tangent,
            dt_s=dt_s,
            candidate=candidate,
            ordinary_assembly=ordinary,
            newborn_center_distance_m=newborn_distance,
            newborn_temperature_k=newborn_temperature,
            newborn_y_hexane=newborn_y,
            newborn_trace_gradients=trace_gradients,
            interface_water_flux_difference_mol_m2_s=math.fsum(
                (
                    interface_flux.component.conserved_water_flux_mol_m2_s,
                    -direct.component.conserved_water_flux_mol_m2_s,
                )
            ),
            interface_hexane_flux_difference_mol_m2_s=math.fsum(
                (
                    interface_flux.component.conserved_hexane_flux_mol_m2_s,
                    -direct.component.conserved_hexane_flux_mol_m2_s,
                )
            ),
            interface_energy_flux_difference_w_m2=math.fsum(
                (
                    interface_flux.energy.total_energy_flux_w_m2,
                    -direct.energy.total_energy_flux_w_m2,
                )
            ),
            outer_water_flux_difference_mol_m2_s=math.fsum(
                (
                    outer_flux.component.conserved_water_flux_mol_m2_s,
                    -target_outer.component.conserved_water_flux_mol_m2_s,
                )
            ),
            outer_hexane_flux_difference_mol_m2_s=math.fsum(
                (
                    outer_flux.component.conserved_hexane_flux_mol_m2_s,
                    -target_outer.component.conserved_hexane_flux_mol_m2_s,
                )
            ),
            outer_energy_flux_difference_w_m2=math.fsum(
                (
                    outer_flux.energy.total_energy_flux_w_m2,
                    -target_outer.energy.total_energy_flux_w_m2,
                )
            ),
        )
    except FaceTangentStepError:
        raise
    except Exception as exc:
        if before is None:
            raise
        raise FaceTangentStepError(
            f"face-tangent positive-sweep continuation rejected with exact rollback: {exc}",
            before,
        ) from exc


def _validate_before(before: face_solver.FaceArrivalCommittedState) -> None:
    if not isinstance(before, face_solver.FaceArrivalCommittedState):
        raise TypeError("face tangent requires FaceArrivalCommittedState")
    k = before.arrival_face_index
    if not 0 < k < before.geometry.master_grid.n:
        raise FaceTangentTopologyError("face tangent requires an interior master face")
    if before.geometry.cut_cell_index is not None or not before.geometry.front.at_master_face:
        raise FaceTangentTopologyError("face tangent requires exact no-cut geometry")
    if len(before.dry_total_stefan_fluxes_mol_m2_s) != len(
        before.dry_temperatures_k
    ) + 1:
        raise FaceTangentTopologyError(
            "committed dry Stefan-flux array does not span the dry domain"
        )


def _face_data(
    before: face_solver.FaceArrivalCommittedState,
    interface_temperature_k: float,
    outer_total_stefan_flux_mol_m2_s: float,
) -> tuple[
    cut.CutInterfaceState,
    ww.WetWaterCellState,
    ww.WetWaterFaceFlux,
    cp.EquilibriumPoreState,
    cut.DryFaceFlux,
    FaceTangentForceMap,
    float,
]:
    k = before.arrival_face_index
    config = before.config
    interface = cut.evaluate_interface_state(
        interface_temperature_k,
        config,
        before.historical_hexane_loadings[k - 1],
        before.oil_fraction_labels[k - 1],
    )
    wet_cell = ww.evaluate_cell(
        before.wet_temperatures_k[-1],
        before.wet_retained_water_loadings[-1],
        before.historical_hexane_loadings[k - 1],
        before.oil_fraction_labels[k - 1],
        config.wet,
    )
    wet_center = 0.5 * math.fsum(
        (
            before.geometry.master_grid.faces[k - 1],
            before.geometry.master_grid.faces[k],
        )
    )
    wet_flux = cut._wet_flux(  # noqa: SLF001
        wet_cell.temperature_k,
        wet_cell.retained_water_loading,
        interface.temperature_k,
        interface.wet.retained_water_loading,
        before.geometry.front.radius_m - wet_center,
        config,
        "wet cell k-1",
        "exact-face interface",
        left_capacity_dual_over_rt=(
            wet_cell.retained_water_capacity_dual_over_rt
        ),
        right_capacity_dual_over_rt=(
            interface.retained_water_trace_capacity_dual_over_rt
        ),
    )
    outer_pore = replace(config.dry.pore, w_o=before.oil_fraction_labels[k])
    cut._validate_open_dry_band(  # noqa: SLF001
        before.dry_temperatures_k[0],
        before.dry_y_hexane[0],
        config,
        pore=outer_pore,
    )
    outer_cell = cp.evaluate_equilibrium(
        before.dry_temperatures_k[0],
        config.dry.pressure_pa,
        before.dry_y_hexane[0],
        outer_pore,
    )
    outer_center = 0.5 * math.fsum(
        (
            before.geometry.master_grid.faces[k],
            before.geometry.master_grid.faces[k + 1],
        )
    )
    newborn_pore = replace(
        config.dry.pore,
        w_o=before.oil_fraction_labels[k - 1],
    )
    outer_flux = cut._dry_flux(  # noqa: SLF001
        interface.temperature_k,
        interface.y_hexane,
        outer_cell.temperature_k,
        outer_cell.y_hexane,
        outer_center - before.geometry.front.radius_m,
        outer_total_stefan_flux_mol_m2_s,
        newborn_pore,
        outer_pore,
        config,
    )
    force_map = _one_sided_force_map(interface, newborn_pore, config)
    mobility = _interface_binary_mobility(interface, config)
    return (
        interface,
        wet_cell,
        wet_flux,
        outer_cell,
        outer_flux,
        force_map,
        mobility,
    )


def _one_sided_force_map(
    interface: cut.CutInterfaceState,
    pore: cp.CoupledPoreParams,
    config: cut.CutTransportConfig,
) -> FaceTangentForceMap:
    temperature = interface.temperature_k
    y_hexane = interface.y_hexane
    derivatives = cp.isothermal_potential_composition_derivatives(
        temperature,
        config.dry.pressure_pa,
        y_hexane,
        pore,
    )
    if derivatives.retained_water_active_set is not interface.dry.retained_water_active_set:
        raise FaceTangentTopologyError(
            "analytic tangent active set disagrees with the interface equilibrium state"
        )
    seam = cp.smooth_luikov_composition_lower_bound(
        temperature,
        config.dry.pressure_pa,
        pore,
    )
    if y_hexane == seam:
        raise FaceTangentTopologyError(
            "exact retained-cap/Luikov seam requires an unambiguous one-sided "
            "active-set direction; no branch smoothing is authorized"
        )
    gas_derivative = (
        derivatives.gas_exchange_potential_derivative_per_y_hexane
    )
    retained_derivative = (
        derivatives.retained_water_potential_derivative_per_y_hexane
    )
    if not math.isfinite(gas_derivative) or gas_derivative == 0.0:
        raise FaceTangentTopologyError(
            "analytic gas thermodynamic-force derivative must be finite and nonzero"
        )
    if not math.isfinite(retained_derivative):
        raise FaceTangentTopologyError(
            "analytic retained-water thermodynamic-force derivative must be finite"
        )
    if derivatives.retained_water_active_set is cp.RetainedWaterActiveSet.RETAINED_CAP:
        if retained_derivative != 0.0:
            raise FaceTangentTopologyError(
                "retained-cap tangent must have an exactly zero force derivative"
            )
        ratio = 0.0
    else:
        ratio = retained_derivative / gas_derivative
    if not math.isfinite(ratio) or ratio < 0.0:
        raise FaceTangentTopologyError(
            "local-equilibrium force gradients lost their common direction"
        )
    return FaceTangentForceMap(
        gas_force_derivative_per_y=gas_derivative,
        retained_force_derivative_per_y=retained_derivative,
        retained_per_gas_gradient_ratio=ratio,
        retained_water_active_set=derivatives.retained_water_active_set,
    )


def _interface_binary_mobility(
    interface: cut.CutInterfaceState,
    config: cut.CutTransportConfig,
) -> float:
    dry = interface.dry
    face = tc.symmetric_binary_face_state(
        dry.binary_gas.molar_density_mol_m3,
        dry.y_water,
        dry.y_hexane,
        dry.binary_gas.molar_density_mol_m3,
        dry.y_water,
        dry.y_hexane,
    )
    result = tc.evaluate_binary_mobility_interval(
        config.dry.binary_diffusivity.interval,
        face,
    ).select(config.dry.binary_diffusivity.fraction)
    mobility = result.value_mol_m_s
    if not math.isfinite(mobility) or mobility <= 0.0:
        raise FaceTangentTopologyError(
            "exact-face binary mobility must be positive and finite"
        )
    return mobility


def _direct_flux(
    before: face_solver.FaceArrivalCommittedState,
    interface: cut.CutInterfaceState,
    force_map: FaceTangentForceMap,
    outer_flux: cut.DryFaceFlux,
    *,
    temperature_gradient_k_m: float,
    gas_force_gradient_m_inv: float,
) -> FaceTangentDirectFlux:
    values = (temperature_gradient_k_m, gas_force_gradient_m_inv)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("direct interface gradients must be finite")
    config = before.config
    binary_mobility = _interface_binary_mobility(interface, config)
    retained_mobility = config.dry.retained_water_mobility.value_mol_m_s
    retained_gradient = (
        force_map.retained_per_gas_gradient_ratio * gas_force_gradient_m_inv
    )
    independent_water = -binary_mobility * gas_force_gradient_m_inv
    retained_water = -retained_mobility * retained_gradient
    outer_total_components = math.fsum(
        (
            outer_flux.component.conserved_water_flux_mol_m2_s,
            outer_flux.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    recovered_stefan = math.fsum((outer_total_components, -retained_water))
    component = cp.compose_component_fluxes(
        interface.dry.y_water,
        interface.dry.y_hexane,
        recovered_stefan,
        independent_water,
        retained_water,
    )
    conductive = (
        -config.dry.thermal_conductivity.value_w_m_k
        * temperature_gradient_k_m
    )
    energy = cp.component_energy_flux(
        component,
        conductive,
        interface.dry.water_gas_partial_enthalpy_j_mol,
        interface.dry.hexane_gas_partial_enthalpy_j_mol,
        interface.dry.retained_water_enthalpy_j_mol,
    )
    return FaceTangentDirectFlux(
        temperature_gradient_k_m=temperature_gradient_k_m,
        gas_force_gradient_m_inv=gas_force_gradient_m_inv,
        retained_force_gradient_m_inv=retained_gradient,
        binary_mobility_mol_m_s=binary_mobility,
        retained_water_mobility_mol_m_s=retained_mobility,
        independent_water_flux_mol_m2_s=independent_water,
        retained_water_flux_mol_m2_s=retained_water,
        conductive_heat_flux_w_m2=conductive,
        recovered_total_stefan_flux_mol_m2_s=recovered_stefan,
        component=component,
        energy=energy,
    )


def _newborn_composition_from_actual_path_force_gradient(
    before: face_solver.FaceArrivalCommittedState,
    interface: cut.CutInterfaceState,
    newborn_temperature_k: float,
    distance_m: float,
    gas_force_gradient_m_inv: float,
) -> tuple[float, float]:
    """Invert the finite actual-path gas force without a virtual ``T_Gamma`` state.

    The exact one-sided tangent supplies the local predictor.  A deterministic
    expansion follows only the predictor's Newton direction and every sampled
    point must remain in the same admitted actual-path component.  An invalid
    sample is a topology failure, not a signal to search another component or
    reconstruction.  The unchanged first internal dry face remains symmetric,
    so its common-temperature endpoint domain still participates in the open
    newborn interval.
    """

    config = before.config
    k = before.arrival_face_index
    pore = replace(config.dry.pore, w_o=before.oil_fraction_labels[k - 1])
    if not all(
        math.isfinite(value)
        for value in (
            newborn_temperature_k,
            distance_m,
            gas_force_gradient_m_inv,
        )
    ) or distance_m <= 0.0:
        raise FaceTangentTopologyError(
            "actual-path newborn inversion requires finite inputs and positive distance"
        )
    if not before.dry_temperatures_k:
        raise FaceTangentTopologyError(
            "actual-path newborn inversion requires the first outer dry cell"
        )

    internal_face_temperature = 0.5 * math.fsum(
        (newborn_temperature_k, before.dry_temperatures_k[0])
    )
    storage_interval = cp.gas_only_composition_interval(
        newborn_temperature_k,
        config.dry.pressure_pa,
        pore,
    )
    internal_interval = cp.gas_only_composition_interval(
        internal_face_temperature,
        config.dry.pressure_pa,
        pore,
    )
    lower = max(
        storage_interval.lower_y_hexane,
        internal_interval.lower_y_hexane,
    )
    upper = min(
        storage_interval.upper_y_hexane,
        internal_interval.upper_y_hexane,
    )
    if not lower < upper:
        raise FaceTangentTopologyError(
            "actual-path newborn storage and first internal face have no common "
            "open gas-only interval"
        )

    force_map = _one_sided_force_map(interface, pore, config)
    force_derivative = force_map.gas_force_derivative_per_y
    predictor_gradient = gas_force_gradient_m_inv / force_derivative
    predictor_y = math.fsum(
        (interface.y_hexane, predictor_gradient * distance_m)
    )
    if not lower < predictor_y < upper:
        raise FaceTangentTopologyError(
            "actual-path newborn tangent predictor left its open storage/internal-face "
            "domain; do not clip"
        )

    evaluations: dict[float, tuple[float, float]] = {}

    def evaluate_gradient(composition_gradient_m_inv: float) -> tuple[float, float]:
        if not math.isfinite(composition_gradient_m_inv):
            raise FaceTangentTopologyError(
                "actual-path newborn bracket produced a non-finite gradient"
            )
        cached = evaluations.get(composition_gradient_m_inv)
        if cached is not None:
            return cached
        composition = math.fsum(
            (
                interface.y_hexane,
                composition_gradient_m_inv * distance_m,
            )
        )
        if not lower < composition < upper:
            raise FaceTangentTopologyError(
                "actual-path newborn bracket left its open storage/internal-face domain"
            )
        try:
            path = acfp.evaluate_actual_composition_force_path(
                interface.temperature_k,
                interface.y_hexane,
                newborn_temperature_k,
                composition,
                config.dry.pressure_pa,
                distance_m,
                pore,
                closed_hexane_derivative_authority=(
                    acfp.ClosedHexaneDerivativeAuthority.NATIVE_LOCAL_ANALYTIC_ENVELOPE
                ),
            )
        except (
            acfp.ActualCompositionForceConvergenceError,
            cp.CoupledPoreTopologyError,
        ) as exc:
            raise FaceTangentTopologyError(
                "actual-path newborn force bracket left its admitted path component"
            ) from exc
        residual = math.fsum(
            (
                path.gas_force_gradient_m_inv,
                -gas_force_gradient_m_inv,
            )
        )
        if not math.isfinite(residual):
            raise FaceTangentTopologyError(
                "actual-path newborn force residual became non-finite"
            )

        # The exact tangent derivative fixes the local force-map orientation.
        # Reject a sampled reversal rather than accepting an ambiguous second
        # root.  This is a finite numerical audit, not a continuous proof.
        for other_gradient, (_, other_residual) in evaluations.items():
            if other_gradient == composition_gradient_m_inv:
                continue
            orientation = (
                (residual - other_residual)
                * (composition_gradient_m_inv - other_gradient)
                * force_derivative
            )
            if orientation < 0.0:
                raise FaceTangentTopologyError(
                    "actual-path newborn force map has an ambiguous sampled orientation"
                )
        result = (composition, residual)
        evaluations[composition_gradient_m_inv] = result
        return result

    _, predictor_residual = evaluate_gradient(predictor_gradient)
    if predictor_residual == 0.0:
        gradient = predictor_gradient
    else:
        newton_step = -predictor_residual / force_derivative
        if not math.isfinite(newton_step) or newton_step == 0.0:
            raise FaceTangentTopologyError(
                "actual-path newborn predictor has no representable bracket direction"
            )
        direction = math.copysign(1.0, newton_step)
        radius = abs(newton_step)
        previous_absolute_residual = abs(predictor_residual)
        bracket: tuple[float, float] | None = None
        exact_expansion_root: float | None = None
        for _ in range(100):
            trial_gradient = math.fsum(
                (predictor_gradient, direction * radius)
            )
            if trial_gradient == predictor_gradient:
                radius *= 2.0
                if not math.isfinite(radius):
                    break
                continue
            trial_y = math.fsum(
                (interface.y_hexane, trial_gradient * distance_m)
            )
            if not lower < trial_y < upper:
                break
            _, trial_residual = evaluate_gradient(trial_gradient)
            if trial_residual == 0.0:
                exact_expansion_root = trial_gradient
                break
            if (trial_residual < 0.0) != (predictor_residual < 0.0):
                bracket = tuple(sorted((predictor_gradient, trial_gradient)))
                break
            if abs(trial_residual) >= previous_absolute_residual:
                raise FaceTangentTopologyError(
                    "actual-path newborn expansion did not approach one local force root"
                )
            previous_absolute_residual = abs(trial_residual)
            radius *= 2.0
            if not math.isfinite(radius):
                break
        if exact_expansion_root is not None:
            gradient = exact_expansion_root
        elif bracket is None:
            raise FaceTangentTopologyError(
                "actual-path newborn force continuation has no bracketed local root"
            )
        else:
            gradient, root_result = optimize.brentq(
                lambda value: evaluate_gradient(value)[1],
                bracket[0],
                bracket[1],
                xtol=1.0e-10,
                rtol=1.0e-15,
                maxiter=100,
                full_output=True,
                disp=False,
            )
            if not root_result.converged:
                raise FaceTangentTopologyError(
                    "actual-path newborn bracketed force root did not converge"
                )

    root = math.fsum((interface.y_hexane, gradient * distance_m))
    if not lower < root < upper:
        raise FaceTangentTopologyError(
            "actual-path newborn force root left its open storage/internal-face domain"
        )
    cut._validate_open_dry_band(  # noqa: SLF001
        newborn_temperature_k,
        root,
        config,
        pore=pore,
    )
    try:
        final_path = acfp.evaluate_actual_composition_force_path(
            interface.temperature_k,
            interface.y_hexane,
            newborn_temperature_k,
            root,
            config.dry.pressure_pa,
            distance_m,
            pore,
            closed_hexane_derivative_authority=(
                acfp.ClosedHexaneDerivativeAuthority.NATIVE_LOCAL_ANALYTIC_ENVELOPE
            ),
        )
    except (
        acfp.ActualCompositionForceConvergenceError,
        cp.CoupledPoreTopologyError,
    ) as exc:
        raise FaceTangentTopologyError(
            "actual-path newborn force root failed final path certification"
        ) from exc
    final_residual = math.fsum(
        (
            final_path.gas_force_gradient_m_inv,
            -gas_force_gradient_m_inv,
        )
    )
    if not math.isfinite(final_residual) or abs(final_residual) > abs(
        predictor_residual
    ):
        raise FaceTangentTopologyError(
            "actual-path newborn force root failed its final residual improvement"
        )
    return root, gradient


def _newborn_composition_from_force_gradient(
    before: face_solver.FaceArrivalCommittedState,
    interface: cut.CutInterfaceState,
    newborn_temperature_k: float,
    distance_m: float,
    gas_force_gradient_m_inv: float,
) -> tuple[float, float]:
    """Recover the finite newborn composition from the exact gas-force trace.

    The root spans the complete common gas-only interval at the interface
    force temperature and the newborn storage temperature.  It does not push
    a retained-cap tangent onto the smooth Luikov side.  If a finite root
    crosses the retained-cap kink, the public piecewise secant supplies its
    unique endpoint force difference; the already-rejected exact seam tangent
    is never regularized here.
    """

    if before.config.moving_interface_composition_force_authority is (
        cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
    ):
        return _newborn_composition_from_actual_path_force_gradient(
            before,
            interface,
            newborn_temperature_k,
            distance_m,
            gas_force_gradient_m_inv,
        )

    config = before.config
    k = before.arrival_face_index
    pore = replace(config.dry.pore, w_o=before.oil_fraction_labels[k - 1])
    face_temperature = interface.temperature_k
    face_interval = cp.gas_only_composition_interval(
        face_temperature,
        config.dry.pressure_pa,
        pore,
    )
    storage_interval = cp.gas_only_composition_interval(
        newborn_temperature_k,
        config.dry.pressure_pa,
        pore,
    )
    lower = max(
        face_interval.lower_y_hexane,
        storage_interval.lower_y_hexane,
    )
    upper = min(
        face_interval.upper_y_hexane,
        storage_interval.upper_y_hexane,
    )
    if not lower < upper:
        raise FaceTangentTopologyError(
            "newborn force root has no common gas-only interval"
        )

    def is_common_admissible(composition: float) -> bool:
        try:
            cp.evaluate_equilibrium(
                face_temperature,
                config.dry.pressure_pa,
                composition,
                pore,
            )
            cp.evaluate_equilibrium(
                newborn_temperature_k,
                config.dry.pressure_pa,
                composition,
                pore,
            )
        except cp.CoupledPoreTopologyError:
            return False
        return True

    def inward_representable_endpoint(boundary: float, toward: float) -> float:
        composition = boundary
        for _ in range(128):
            composition = math.nextafter(composition, toward)
            if not lower < composition < upper:
                break
            if is_common_admissible(composition):
                return composition

        # An activity root can be more than 128 composition ulps outside its
        # admissible side when its local activity slope is very small.  Contract
        # toward a known strict interior point; this only constructs a root
        # bracket and never returns a projected physical candidate.
        strict_interior = boundary + 0.5 * (toward - boundary)
        if not lower < strict_interior < upper or not is_common_admissible(
            strict_interior
        ):
            raise FaceTangentTopologyError(
                "newborn force root has no representable common gas-only interior"
            )
        rejected = boundary
        accepted = strict_interior
        for _ in range(96):
            trial = rejected + 0.5 * (accepted - rejected)
            if trial == rejected or trial == accepted:
                break
            if is_common_admissible(trial):
                accepted = trial
            else:
                rejected = trial
        if lower < accepted < upper and is_common_admissible(accepted):
            return accepted
        raise FaceTangentTopologyError(
            "newborn force root has no representable common gas-only endpoint"
        )

    lower_open = inward_representable_endpoint(lower, upper)
    upper_open = inward_representable_endpoint(upper, lower)
    if not lower_open < upper_open:
        raise FaceTangentTopologyError(
            "newborn force root has no common open gas-only interval"
        )
    lower_gradient = (lower_open - interface.y_hexane) / distance_m
    upper_gradient = (upper_open - interface.y_hexane) / distance_m

    def residual(composition_gradient_m_inv: float) -> float:
        delta_y = composition_gradient_m_inv * distance_m
        secant = cp.isothermal_potential_composition_secant_slopes(
            face_temperature,
            config.dry.pressure_pa,
            interface.y_hexane,
            delta_y,
            pore,
        )
        return math.fsum(
            (
                secant.gas_exchange_potential_secant_per_y_hexane
                * composition_gradient_m_inv,
                -gas_force_gradient_m_inv,
            )
        )

    f_lower = residual(lower_gradient)
    f_upper = residual(upper_gradient)
    if f_lower == 0.0:
        gradient = lower_gradient
    elif f_upper == 0.0:
        gradient = upper_gradient
    elif f_lower * f_upper >= 0.0:
        raise FaceTangentTopologyError(
            "newborn thermodynamic-force continuation has no bracketed root"
        )
    else:
        gradient = optimize.brentq(
            residual,
            lower_gradient,
            upper_gradient,
            xtol=1.0e-10,
            rtol=1.0e-15,
            maxiter=100,
        )
    root = math.fsum((interface.y_hexane, gradient * distance_m))
    if not lower < root < upper or not is_common_admissible(root):
        raise FaceTangentTopologyError(
            "newborn thermodynamic-force root left the common open gas-only interval"
        )
    force_root = cp.evaluate_equilibrium(
        face_temperature,
        config.dry.pressure_pa,
        root,
        pore,
    )
    if force_root.retained_water_active_set is not interface.dry.retained_water_active_set:
        raise FaceTangentTopologyError(
            "newborn thermodynamic-force root left the exact tangent's one-sided "
            "retained-water active set"
        )
    cut._validate_open_dry_band(  # noqa: SLF001
        newborn_temperature_k,
        root,
        config,
        pore=pore,
    )
    return root, gradient


__all__ = [
    "FaceTangentAssembly",
    "FaceTangentDirectFlux",
    "FaceTangentForceMap",
    "FaceTangentLayout",
    "FaceTangentLedger",
    "FaceTangentOverlap",
    "FaceTangentRankAudit",
    "FaceTangentResidualBlocks",
    "FaceTangentNoAdmissibleRootError",
    "FaceTangentStepError",
    "FaceTangentTopologyError",
    "FaceTangentUnknowns",
    "FACE_TANGENT_INWARD_SEED_ENVIRONMENT_VARIABLE",
    "FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE",
    "INWARD_SEED_BISECTION_STEPS",
    "INWARD_SEED_FLUX_RATIO_GRID",
    "INWARD_SEED_NEWTON_ITERATIONS",
    "assemble_face_tangent",
    "audit_face_tangent_rank",
    "candidate_from_committed_face",
    "continue_to_positive_sweep",
    "drain_face_tangent_call_reports",
    "face_tangent_call_reports",
    "face_tangent_inward_seed_armed",
    "face_tangent_residual_bound",
    "face_tangent_residual_bound_declared",
]
