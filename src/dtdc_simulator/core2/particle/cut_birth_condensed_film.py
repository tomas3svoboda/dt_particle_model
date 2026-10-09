r"""BC-4a: the attached-hexane film channel on the birth committed state.

AUTHORITY.  A1, RATIFIED 2026-08-22
(``docs/GT_PS2_BC3_CONDENSED_SURFACE_BIRTH_2026-08-22.md``, Status block +
section 3), under the ratified M2 design
(``docs/GT_PS2_RR2_BIRTH_CONTINUATION_DESIGN_2026-08-21.md``): the ratified
D2 deposit path extends INTO the birth committed state - the cut state gains
the attached-hexane film channel, the intensive twin of
``surface_active_set.ExternalSurfaceInventories.attached_hexane_kg_per_kg_dry``,
with the CC-2a accumulator discipline (exact ``math.fsum`` identity,
non-negative cone, typed refusals), so the mesh-certified condensed-surface
birth root can COMMIT with the film's mass and energy booked.

THE CONSTRUCTION (beside; the sealed-companion pattern).  Nothing on
``cut_integrator.CutIntegratorState`` changes.  The film channel is a sealed
companion record - :class:`CommittedBirthAttachedFilm` - BOUND TO the
committed state BY IDENTITY (``record.committed_state is step.after``), built
from an ACCEPTED :class:`cut_birth_integrator.BirthIntegratorStep` whose
boundary is the condensed-solvent variant carrying an armed
:class:`condensed_solvent_reduced_film.CondensedRootCommitChannel` with
exactly one recorded root admission.  The booking is construction-exact
against the birth solve's own condensed-root evidence, never an estimate:

* MASS - the condensate is the birth step's own surface hexane throughput,
  ``dt * A * N_h`` with every factor the accepted assembly's own; the closure
  identity against the integrator ledger's ``boundary_hexane_out_mol`` is
  enforced EXACTLY (same expression, bit for bit).  Physically this is A1's
  measured statement: the newborn shell's evaporative load re-condenses on
  its own sub-dew surface, so the hexane crossing the newborn surface lands
  in the attached external film rather than the bulk.  The committed
  particle-side ledgers are untouched - they correctly book the same moles
  as boundary outflow; this record says where that outflow LIVES.
* ENERGY - the convected hexane enthalpy across the same face over the same
  step, ``dt * A * (N_h * h_h)`` with ``h_h`` the accepted root's own
  surface-equilibrium gas partial enthalpy (the assembly's
  ``surface_film_audit.hexane_gas_partial_enthalpy_j_mol``); the identity
  ``energy.gas_hexane_enthalpy_flux_w_m2 == N_h * h_h`` is enforced EXACTLY
  (that is literally how the frozen surface energy constructs it), and the
  float-association residual between the two lawful bookings of the same
  quantity is REPORTED, never absorbed.
* EVIDENCE - the channel's typed root admission is cross-checked EXACTLY
  against the accepted assembly (root surface temperature vs the surface
  film audit, cell temperature vs the committed newborn dry temperature,
  pressure and boundary composition vs the boundary), and the sweep
  journal's own condensed-root count corroborates.

The drainage latch is honored by A1's amendment scope: the film on the
newborn state is the pre-recession epoch's continuation (the D3 birth gate
required the epoch film at EXACTLY zero, enforced here on ``film_before``),
and recession onset still requires the film at zero - film re-evaporation is
the newborn shell's first task (BC-4b), measured, not assumed.

Every claim flag reachable from this module is False.  Nothing here makes
the particle model physically qualifying, plant predictive, or production
wired.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.particle import condensed_solvent_reduced_film as csrf
from dtdc_simulator.core2.particle import cut_birth_integrator as cbi
from dtdc_simulator.core2.particle import surface_active_set as sas
from dtdc_simulator.core2.props import hexane as hx


class BirthCondensedFilmRefusalCode(Enum):
    NOT_A_BIRTH_STEP = "not_a_birth_step"
    UNACCEPTED_STEP = "unaccepted_step"
    FILM_BEFORE_NOT_TYPED = "film_before_not_typed"
    FILM_BEFORE_NOT_ZERO = "film_before_not_zero"
    NOT_CONDENSED_VARIANT = "not_condensed_variant"
    NO_COMMIT_CHANNEL = "no_commit_channel"
    NO_CONDENSED_ROOT_ADMISSION = "no_condensed_root_admission"
    JOURNAL_ROOT_EVIDENCE_MISSING = "journal_root_evidence_missing"
    SURFACE_AUDIT_MISSING = "surface_audit_missing"
    ADMISSION_STATE_MISMATCH = "admission_state_mismatch"
    NONPOSITIVE_DEPOSIT = "nonpositive_deposit"
    MASS_CLOSURE_BROKEN = "mass_closure_broken"
    ENERGY_IDENTITY_BROKEN = "energy_identity_broken"
    NONFINITE_INPUT = "nonfinite_input"


class BirthCondensedFilmError(ValueError):
    """Typed refusal for the birth attached-hexane film channel; fail closed."""

    def __init__(self, code: BirthCondensedFilmRefusalCode, message: str) -> None:
        super().__init__(f"{code.value}: {message}")
        self.code = code


@dataclass(frozen=True)
class BirthCondensedFilmAccumulator:
    """The attached-hexane film booked by one committed condensed-surface birth.

    The exact mirror of CC-2a's ``FreeWaterFilmAccumulator`` and BC-1's
    ``QSCHexaneFilmAccumulator``: the intensive twin is carried on the
    particle's own dry-mass basis, the accumulator identity holds
    construction-exactly through ``math.fsum``, and leaving the non-negative
    cone is a typed refusal, never a clamp.  A birth film step may not move
    the free water.  ``rebasing_residual_kg`` is the one float rounding the
    intensive re-basing introduces (``(mol*M)/basis*basis`` vs ``mol*M``),
    reported as evidence rather than hidden; the PRIMARY mass closure is in
    moles and exact.
    """

    before: sas.ExternalSurfaceInventories
    after: sas.ExternalSurfaceInventories
    basis_kg_dry: float
    duration_s: float
    deposit_kg_per_kg_dry: float
    deposit_kg: float
    deposit_mol: float
    deposit_energy_j: float
    molar_enthalpy_j_mol: float
    rebasing_residual_kg: float
    exact_accumulator_identity: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        for payload in (self.before, self.after):
            if type(payload) is not sas.ExternalSurfaceInventories:
                raise BirthCondensedFilmError(
                    BirthCondensedFilmRefusalCode.FILM_BEFORE_NOT_TYPED,
                    "the film accumulator needs typed ExternalSurfaceInventories",
                )
        values = (
            self.basis_kg_dry,
            self.duration_s,
            self.deposit_kg_per_kg_dry,
            self.deposit_kg,
            self.deposit_mol,
            self.deposit_energy_j,
            self.molar_enthalpy_j_mol,
            self.rebasing_residual_kg,
        )
        if not all(math.isfinite(value) for value in values):
            raise BirthCondensedFilmError(
                BirthCondensedFilmRefusalCode.NONFINITE_INPUT,
                "film accumulator values must be finite",
            )
        if self.basis_kg_dry <= 0.0 or self.duration_s <= 0.0:
            raise BirthCondensedFilmError(
                BirthCondensedFilmRefusalCode.NONFINITE_INPUT,
                "the film accumulator basis and duration must be strictly positive",
            )
        if self.after.free_water_kg_per_kg_dry != self.before.free_water_kg_per_kg_dry:
            raise BirthCondensedFilmError(
                BirthCondensedFilmRefusalCode.NONFINITE_INPUT,
                "a birth hexane film step may not move the free water",
            )
        if self.after.attached_hexane_kg_per_kg_dry != math.fsum(
            (self.before.attached_hexane_kg_per_kg_dry, self.deposit_kg_per_kg_dry)
        ):
            raise BirthCondensedFilmError(
                BirthCondensedFilmRefusalCode.NONFINITE_INPUT,
                "the film accumulator identity must hold construction-exactly",
            )
        if self.deposit_kg != self.deposit_kg_per_kg_dry * self.basis_kg_dry:
            raise BirthCondensedFilmError(
                BirthCondensedFilmRefusalCode.NONFINITE_INPUT,
                "the extensive film deposit must be the intensive gain times basis",
            )
        if not self.deposit_kg_per_kg_dry > 0.0:
            raise BirthCondensedFilmError(
                BirthCondensedFilmRefusalCode.NONPOSITIVE_DEPOSIT,
                "a committed condensed-surface birth books a strictly positive "
                "condensate (the receding branch certifies a positive outward "
                "surface hexane flux); a nonpositive deposit is a contradiction",
            )
        # The non-negative cone is ExternalSurfaceInventories' own; with a
        # positive deposit onto a lawful before-inventory it holds by
        # construction, and the typed inventory constructor enforces it.

    @property
    def topology_after(self) -> sas.SurfaceTopology:
        return self.after.topology

    def journal_record(self) -> dict:
        return {
            "accumulator": "birth-condensed-attached-hexane-film",
            "attached_hexane_before_kg_per_kg_dry": (
                self.before.attached_hexane_kg_per_kg_dry
            ),
            "attached_hexane_after_kg_per_kg_dry": (
                self.after.attached_hexane_kg_per_kg_dry
            ),
            "free_water_kg_per_kg_dry": self.before.free_water_kg_per_kg_dry,
            "basis_kg_dry": self.basis_kg_dry,
            "duration_s": self.duration_s,
            "deposit_kg_per_kg_dry": self.deposit_kg_per_kg_dry,
            "deposit_kg": self.deposit_kg,
            "deposit_mol": self.deposit_mol,
            "deposit_energy_j": self.deposit_energy_j,
            "molar_enthalpy_j_mol": self.molar_enthalpy_j_mol,
            "rebasing_residual_kg": self.rebasing_residual_kg,
            "topology_after": self.after.topology.value,
            "physically_qualifying": False,
        }


@dataclass(frozen=True)
class CommittedBirthAttachedFilm:
    """Sealed companion record: the newborn state's attached-hexane film.

    Bound to the committed :class:`cut_integrator.CutIntegratorState` BY
    IDENTITY - consumers gate on :meth:`is_bound_to`, never on equality - so
    the record can neither be re-attached to another state nor survive a
    state it does not describe.  All closure identities were enforced at
    construction; the residual fields carry the evidence.
    """

    committed_state: object
    film: BirthCondensedFilmAccumulator
    admission: csrf.CondensedRootCommitAdmission
    condensed_root_evaluations: int
    condensed_sweep_samples: int
    coldest_condensed_root_temperature_k: float
    surface_hexane_flux_mol_m2_s: float
    surface_gas_hexane_enthalpy_flux_w_m2: float
    hexane_gas_partial_enthalpy_j_mol: float
    #: EXACT (enforced == 0.0): dt*A*N_h vs the integrator ledger's
    #: boundary_hexane_out_mol - the same expression, bit for bit.
    mass_closure_residual_mol: float
    #: EXACT (enforced == 0.0): the assembly's gas-hexane enthalpy flux vs
    #: the recomputed N_h * h_h - the frozen construction, bit for bit.
    enthalpy_flux_identity_residual_w_m2: float
    #: Reported evidence only: (dt*A)*(N_h*h_h) vs (dt*A*N_h)*h_h - the two
    #: lawful float associations of one exact quantity.
    energy_association_residual_j: float
    physically_qualifying: bool = field(default=False, init=False)

    def is_bound_to(self, state: object) -> bool:
        return state is self.committed_state

    def journal_record(self) -> dict:
        return {
            "record": "committed-birth-attached-hexane-film",
            "film": self.film.journal_record(),
            "admission": self.admission.journal_record(),
            "condensed_root_evaluations": self.condensed_root_evaluations,
            "condensed_sweep_samples": self.condensed_sweep_samples,
            "coldest_condensed_root_temperature_k": (
                self.coldest_condensed_root_temperature_k
            ),
            "surface_hexane_flux_mol_m2_s": self.surface_hexane_flux_mol_m2_s,
            "surface_gas_hexane_enthalpy_flux_w_m2": (
                self.surface_gas_hexane_enthalpy_flux_w_m2
            ),
            "hexane_gas_partial_enthalpy_j_mol": (
                self.hexane_gas_partial_enthalpy_j_mol
            ),
            "mass_closure_residual_mol": self.mass_closure_residual_mol,
            "enthalpy_flux_identity_residual_w_m2": (
                self.enthalpy_flux_identity_residual_w_m2
            ),
            "energy_association_residual_j": self.energy_association_residual_j,
            "physically_qualifying": False,
        }


def require_pre_birth_film_admissible(
    film_before: sas.ExternalSurfaceInventories,
) -> None:
    """Refuse a pre-birth external inventory the D3 gate could not have passed.

    The typed birth gate admits only with the epoch's attached-hexane film at
    EXACTLY zero (float equality - the drainage latch's and the radial
    activation endpoint's own requirement), so the newborn film is entirely
    the birth step's own condensate.  Enforced independently here so the
    companion record can never smuggle a pre-existing film into the booking.
    """

    if type(film_before) is not sas.ExternalSurfaceInventories:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.FILM_BEFORE_NOT_TYPED,
            "film_before must be typed ExternalSurfaceInventories",
        )
    if film_before.attached_hexane_kg_per_kg_dry != 0.0:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.FILM_BEFORE_NOT_ZERO,
            "the D3 birth gate admits only with the epoch film at EXACTLY zero "
            f"(got {film_before.attached_hexane_kg_per_kg_dry!r} kg/kg dry); the "
            "newborn film must be entirely the birth step's own condensate",
        )


def commit_birth_attached_film(
    step: cbi.BirthIntegratorStep,
    *,
    film_before: sas.ExternalSurfaceInventories,
) -> CommittedBirthAttachedFilm:
    """Book the accepted condensed-surface birth's condensate, exactly.

    ``film_before`` is the pre-birth external inventory - the epoch's ending
    state, whose attached hexane the D3 birth gate certified at EXACTLY zero
    (float equality; enforced again here so the newborn film is entirely the
    birth step's own condensate).  Every quantity is the accepted assembly's
    own; every closure identity is enforced exactly; nothing is estimated.
    """

    if type(step) is not cbi.BirthIntegratorStep:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.NOT_A_BIRTH_STEP,
            "the film channel books exactly one accepted BirthIntegratorStep; "
            f"got {type(step).__name__}",
        )
    require_pre_birth_film_admissible(film_before)
    if not step.ledger.accepted:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.UNACCEPTED_STEP,
            "only an accepted birth step can book its condensate",
        )
    boundary = step.boundary
    if type(boundary) is not csrf.CondensedSolventReducedFilmPoreBoundary:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.NOT_CONDENSED_VARIANT,
            "the film channel requires the exact condensed-solvent reduced-film "
            f"boundary variant; got {type(boundary).__name__}: a Dirichlet or "
            "plain reduced-film birth has no condensed-root evidence to book",
        )
    channel = boundary.condensed_root_commit_channel
    if channel is None:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.NO_COMMIT_CHANNEL,
            "the boundary carries no condensed-root commit channel: an accepted "
            "condensed root cannot have committed without one (the gate refuses "
            "channel-less condensed roots), and a gas-only birth has no "
            "condensate to book",
        )
    admission = channel.admission
    if admission is None:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.NO_CONDENSED_ROOT_ADMISSION,
            "the commit channel is empty: the accepted root evaluated gas-only "
            "(above the dew locus), so there is no condensate to book and no "
            "film record to build; commit the gas-only birth without one",
        )
    journal = boundary.sweep_journal
    if journal.condensed_root_evaluations < 1:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.JOURNAL_ROOT_EVIDENCE_MISSING,
            "the sweep journal carries no condensed-root evaluation: the "
            "channel admission and the journal must corroborate one another",
        )
    assembly = step.assembly
    surface = assembly.dry_face_fluxes[-1]
    audit = surface.surface_film_audit
    if audit is None:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.SURFACE_AUDIT_MISSING,
            "the accepted assembly's surface flux carries no reduced-film "
            "audit; the condensate's enthalpy reference is missing",
        )
    mismatches = []
    if admission.surface_temperature_k != audit.surface_temperature_k:
        mismatches.append(
            "root surface temperature "
            f"{admission.surface_temperature_k!r} != audit "
            f"{audit.surface_temperature_k!r}"
        )
    if admission.cell_temperature_k != assembly.candidate.dry_temperature_k:
        mismatches.append(
            f"root cell temperature {admission.cell_temperature_k!r} != "
            f"newborn dry temperature {assembly.candidate.dry_temperature_k!r}"
        )
    if admission.pressure_pa != assembly.config.dry.pressure_pa:
        mismatches.append(
            f"admission pressure {admission.pressure_pa!r} != chart pressure "
            f"{assembly.config.dry.pressure_pa!r}"
        )
    if admission.boundary_y_hexane != boundary.y_hexane:
        mismatches.append(
            f"admission boundary composition {admission.boundary_y_hexane!r} != "
            f"boundary {boundary.y_hexane!r}"
        )
    if mismatches:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.ADMISSION_STATE_MISMATCH,
            "the channel admission does not describe this accepted assembly: "
            + "; ".join(mismatches),
        )
    flux = surface.component.conserved_hexane_flux_mol_m2_s
    if not flux > 0.0:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.NONPOSITIVE_DEPOSIT,
            "the accepted surface hexane flux is not strictly positive "
            f"({flux!r} mol/m2/s); the receding branch cannot have certified",
        )
    dt_s = assembly.dt_s
    area = assembly.before.grid.areas[-1]
    # EXACTLY the integrator ledger's expression (cut_birth_integrator
    # ``_accept_birth``: ``dt_s * area * surface.component
    # .conserved_hexane_flux_mol_m2_s``), bit for bit.
    deposit_mol = dt_s * area * flux
    mass_closure = math.fsum((deposit_mol, -step.ledger.boundary_hexane_out_mol))
    if mass_closure != 0.0:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.MASS_CLOSURE_BROKEN,
            "the condensate does not close against the birth ledger's own "
            f"boundary hexane outflow: dt*A*N_h = {deposit_mol!r} mol vs "
            f"ledger {step.ledger.boundary_hexane_out_mol!r} mol "
            f"(residual {mass_closure!r} mol)",
        )
    h_mol = audit.hexane_gas_partial_enthalpy_j_mol
    flux_e = surface.energy.gas_hexane_enthalpy_flux_w_m2
    # EXACTLY the frozen surface-energy construction
    # (coupled_transport.reduced_film_surface_energy_flux:
    # ``conserved_hexane_flux_mol_m2_s * hexane_gas_partial_enthalpy_j_mol``).
    enthalpy_identity = math.fsum((flux_e, -(flux * h_mol)))
    if enthalpy_identity != 0.0:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.ENERGY_IDENTITY_BROKEN,
            "the surface gas-hexane enthalpy flux does not recompute from the "
            f"accepted root's own partial enthalpy: {flux_e!r} W/m2 vs "
            f"N_h*h_h = {flux * h_mol!r} W/m2",
        )
    deposit_energy_j = dt_s * area * flux_e
    energy_association = math.fsum((deposit_energy_j, -(deposit_mol * h_mol)))
    # The particle's own dry-mass basis: the SAME fsum expression the epoch
    # state's ``basis_kg_dry`` evaluates on the same grid and authority.
    basis_kg_dry = math.fsum(
        volume * assembly.config.wet.wet.rho_dm_p
        for volume in assembly.before.grid.volumes
    )
    deposit_kg_raw = deposit_mol * hx.M
    deposit_intensive = deposit_kg_raw / basis_kg_dry
    deposit_kg = deposit_intensive * basis_kg_dry
    rebasing_residual = math.fsum((deposit_kg_raw, -deposit_kg))
    try:
        film_after = sas.ExternalSurfaceInventories(
            math.fsum((film_before.attached_hexane_kg_per_kg_dry, deposit_intensive)),
            film_before.free_water_kg_per_kg_dry,
        )
    except ValueError as exc:
        raise BirthCondensedFilmError(
            BirthCondensedFilmRefusalCode.NONFINITE_INPUT,
            f"the external surface inventories refused their own cone: {exc}",
        ) from exc
    film = BirthCondensedFilmAccumulator(
        before=film_before,
        after=film_after,
        basis_kg_dry=basis_kg_dry,
        duration_s=dt_s,
        deposit_kg_per_kg_dry=deposit_intensive,
        deposit_kg=deposit_kg,
        deposit_mol=deposit_mol,
        deposit_energy_j=deposit_energy_j,
        molar_enthalpy_j_mol=h_mol,
        rebasing_residual_kg=rebasing_residual,
    )
    return CommittedBirthAttachedFilm(
        committed_state=step.after,
        film=film,
        admission=admission,
        condensed_root_evaluations=journal.condensed_root_evaluations,
        condensed_sweep_samples=journal.condensed_sample_count,
        coldest_condensed_root_temperature_k=(
            journal.coldest_condensed_root_temperature_k
        ),
        surface_hexane_flux_mol_m2_s=flux,
        surface_gas_hexane_enthalpy_flux_w_m2=flux_e,
        hexane_gas_partial_enthalpy_j_mol=h_mol,
        mass_closure_residual_mol=mass_closure,
        enthalpy_flux_identity_residual_w_m2=enthalpy_identity,
        energy_association_residual_j=energy_association,
    )


__all__ = [
    "BirthCondensedFilmAccumulator",
    "BirthCondensedFilmError",
    "BirthCondensedFilmRefusalCode",
    "CommittedBirthAttachedFilm",
    "commit_birth_attached_film",
    "require_pre_birth_film_admissible",
]
