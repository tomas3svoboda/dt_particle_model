"""Study-lane lever W3-2: Gauss-Seidel sweep ordering for the tower pressure Picard.

Owner ruling **W3-2** (``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md``,
2026-09-03), on wave-3 solver-struggle audit finding **3.2**
(``docs/GT_PS2_WAVE3_SOLVER_STRUGGLE_AUDIT_2026-09-03.md``).

The problem
-----------
``qsc_six_tray_k2_mechanism_interval._close_through_pressures`` solves two
unknown gas faces of the through-bed cascade:

* ``sp1_mn2`` — the SP1 top / MN2 inlet face;
* ``mn2_mn1`` — the MN2 top / MN1 inlet face.

The chain must be EVALUATED bottom-up (SP1 → MN2 → MN1), because composition,
molar flow and temperature propagate upward.  But the pressure level propagates
DOWNWARD from the one fixed boundary, ``MN1_TOP_PRESSURE_PA``: MN1's inlet
pressure is the top pressure plus MN1's own drop, and MN2's inlet pressure is
``mn2_mn1`` plus MN2's own drop.  As a system in the two faces the coupling is
therefore triangular, and the measured Jacobian at the operating point is

    J = [[5.379e-03, 9.652e-01],
         [6.105e-03, -1.691e-04]]

— an upper-triangular structure to three digits, with the strong entry
``J01 = 0.9652`` being the pass-through of a change in the upper face into the
lower face's inlet estimate.

The shipped update is ``face_pressures = proposed``: BOTH faces are advanced
from the same stale iterate.  That is a Jacobi sweep of a triangular system,
and it throws away the one number the same chain evaluation already produced —
the new upper face.  Re-derived spectra (audit §3.2, independently recomputed
there from J): Jacobi ρ = **0.079418**; sweeping in cascade order ρ =
**0.014545**; the ratio is **5.460×**.

The lever
---------
Sweep the cascade in its own order — upper face first, lower face consuming the
already-updated upper value — realised as an exact forward substitution inside
the SAME chain evaluation, at zero extra evaluation cost::

    mn2_mn1 <- proposed[1]                       # the just-computed upper face
    sp1_mn2 <- proposed[0] + (proposed[1] - mn2_mn1_old)
             = proposed[0] + residuals[1]

The substitution coefficient is **1.0**, not the measured 0.9652: raising the
top boundary of a tray raises its inlet pressure one-for-one at leading order
(the deviation is only the drop's own weak pressure dependence, through gas
density and superficial velocity).  A unit coefficient is derivative-free,
needs no calibration, is era-independent, and — decisively — is EXACT at the
fixed point, where ``residuals[1]`` is 0.0.  Using the fitted 0.9652 would buy
ρ = 0.011189 instead of 0.014545 (audit §4.4), i.e. 1.3× on top of 5.46×, in
exchange for a fitted constant that would have to be re-measured every era.
Not taken; named in the record as an unbuilt improvement.

The fixed point is IDENTICAL by construction
--------------------------------------------
Write ``P(x)`` for the map "evaluate the chain at the iterate ``x`` and read
the two proposed faces", and ``r(x) = P(x) - x``.  The two updates are

    U_jacobi(x) = P(x)
    U_sweep(x)  = P(x) + (r_1(x), 0)

so ``U_sweep(x) = U_jacobi(x) + (r_1(x), 0)``: they differ by a multiple of a
component of the residual, and by nothing else.

*Every fixed point of one is a fixed point of the other.*  If ``P(x*) = x*``
then ``r(x*) = 0`` and both updates return ``x*``.  Conversely, if
``U_sweep(x*) = x*``, its second component gives ``P_1(x*) = x*_1``, hence
``r_1(x*) = 0``, and its first component then gives ``P_0(x*) + 0 = x*_0``.
So ``Fix(U_sweep) = Fix(U_jacobi) = {x : P(x) = x}`` — the same set, with no
tolerance argument anywhere.  In binary64 the statement is bit-exact as well:
at ``r_1 = 0.0`` the update evaluates ``proposed[0] + 0.0``, and IEEE-754
addition of a zero to a finite non-zero float returns that float unchanged.

*The acceptance test is untouched.*  ``residuals``, ``relative``, the
normalizers, the limit ``request.pressure_face_relative_limit`` and the
returned ``QSCMechanismPressureSolveEvidence`` are computed exactly as the
production function computes them, from the same ``P(x)`` at the same ``x``.
Only the NEXT iterate differs.

*What does differ, and by how much.*  Neither lane returns the fixed point; both
return the first iterate that satisfies ``relative(x) <= tol``, together with
the transactions evaluated AT that iterate.  For a contraction of factor ρ,
``|x_accepted - x*| <= |r(x_accepted)| / (1 - ρ)``, so with tol = 1e-10, a face
near 1.05e5 Pa and ρ <= 0.0794 the accepted faces of the two lanes each lie
within **1.2e-05 Pa** of the common root and therefore within 2.4e-05 Pa
(2.3e-10 relative) of each other.  That is why this is a study-lane lever under
the trajectory certificate and not an exact-lane edit: the answers agree to the
gate, not to the bit.  The audit's independent statement of the same fact is
"roots agree to 4.4e-11 relative — inside the 1e-10 gate but not bit-identical".

Era-binding
-----------
The lever replaces a private function of a certified module, so it carries the
same discipline the surrogate pack carries: the source of the production
function is digested at install time and must match
:data:`PRODUCTION_PICARD_SOURCE_DIGEST`.  If the production Picard is edited —
a different residual, a different acceptance, a third face — the lever REFUSES
to install rather than silently substituting a body derived against a function
that no longer exists.
"""

from __future__ import annotations

import hashlib
import inspect
from typing import Any

__all__ = [
    "PATCH_TARGET",
    "PRODUCTION_PICARD_FUNCTION",
    "PRODUCTION_PICARD_SOURCE_DIGEST",
    "GaussSeidelPressureSweep",
    "PicardEraMismatch",
    "production_picard_source_digest",
]

#: The certified module whose face-pressure Picard this lever re-sweeps.
MECHANISM_MODULE = "dtdc_simulator.core2.qsc_six_tray_k2_mechanism_interval"

#: The private function replaced, and the target name stamped into the run
#: declaration's ``installed_targets``.
PRODUCTION_PICARD_FUNCTION = "_close_through_pressures"
PATCH_TARGET = f"qsc_six_tray_k2_mechanism_interval.{PRODUCTION_PICARD_FUNCTION}"

#: ``sha256`` over ``inspect.getsource`` of the production function, at
#: HEAD ``f5a70bb`` (2026-09-03).  The lever's body is a copy of that source
#: with exactly one statement changed — the iterate update — so a change to the
#: original must re-derive the copy, not be absorbed silently.
#:
#: Re-keyed 2026-09-25: HEAD ``f6cd1c65``.  The production function gained a
#: report-only refusal instrument on the ordered-structural-domain raise (a
#: ``+ _ordered_structural_domain_detail(...)`` appended to the same message,
#: naming the evaluation index, both iterates, the admitted band and the
#: failing clause; ``tests/test_core2_mechanism_ordered_domain_refusal_detail.py``).
#: No residual, acceptance, face count or control flow changed — only the
#: source text moved, so only this pin moved.  Re-derived two independent
#: ways (``production_picard_source_digest()`` against the live import, and an
#: independent ``ast.get_source_segment`` extraction of the same function from
#: the file text) and confirmed to agree; the superseded literal above was
#: confirmed to match ``inspect.getsource`` of the pre-instrument function at
#: HEAD ``f6cd1c65`` (the module's committed state immediately before the
#: instrument was added to the working tree), so the move is explained by
#: that edit alone.  Holder re-keyed per the repository rule: a stale pin is
#: repaired by re-keying the HOLDER, never by editing the pin TARGET.
PRODUCTION_PICARD_SOURCE_DIGEST = (
    "sha256:76528e03fb9e58f357efd8da45c786fc62991a6de511238f4e007a58a3643357"
)

#: The forward-substitution coefficient.  1.0 by derivation, not by fit; see
#: the module docstring.
FORWARD_SUBSTITUTION_COEFFICIENT = 1.0


class PicardEraMismatch(RuntimeError):
    """The production Picard is not the function this lever was derived against."""


def _mechanism_module() -> Any:
    from dtdc_simulator.core2 import qsc_six_tray_k2_mechanism_interval

    return qsc_six_tray_k2_mechanism_interval


def production_picard_source_digest(function: Any | None = None) -> str:
    """``sha256`` over the production Picard's own source text."""

    target = (
        getattr(_mechanism_module(), PRODUCTION_PICARD_FUNCTION)
        if function is None
        else function
    )
    source = inspect.getsource(target)
    return "sha256:" + hashlib.sha256(source.encode("utf-8")).hexdigest()


class GaussSeidelPressureSweep:
    """Install/restore the cascade-ordered sweep over the production Picard."""

    name = "gs_pressure_sweep"

    def __init__(self) -> None:
        self._installed = False
        self._original: Any = None
        self.solves = 0
        self.chain_evaluations = 0
        self.sweep_updates = 0
        self.domain_refusals = 0
        self.evaluation_limit_refusals = 0
        self.max_forward_substitution_pa = 0.0

    # -- lever protocol -----------------------------------------------------
    def install(self) -> None:
        if self._installed:
            return
        module = _mechanism_module()
        original = getattr(module, PRODUCTION_PICARD_FUNCTION)
        measured = production_picard_source_digest(original)
        if measured != PRODUCTION_PICARD_SOURCE_DIGEST:
            raise PicardEraMismatch(
                "the production face-pressure Picard has changed since this lever "
                f"was derived: {PATCH_TARGET} digests {measured}, expected "
                f"{PRODUCTION_PICARD_SOURCE_DIGEST}.  Re-derive the sweep against "
                "the current function before re-enabling gs_pressure_sweep."
            )
        self._original = original
        replacement = self._build()
        replacement.__name__ = PRODUCTION_PICARD_FUNCTION
        replacement.__qualname__ = f"study_lane.gs_pressure_sweep.{PATCH_TARGET}"
        setattr(module, PRODUCTION_PICARD_FUNCTION, replacement)
        self._installed = True

    def restore(self) -> None:
        if not self._installed:
            return
        setattr(_mechanism_module(), PRODUCTION_PICARD_FUNCTION, self._original)
        self._original = None
        self._installed = False

    def installed_targets(self) -> list[str]:
        return [PATCH_TARGET] if self._installed else []

    def usage(self) -> dict[str, Any]:
        return {
            "picard_solves": self.solves,
            "chain_evaluations": self.chain_evaluations,
            "sweep_updates": self.sweep_updates,
            "domain_refusals": self.domain_refusals,
            "evaluation_limit_refusals": self.evaluation_limit_refusals,
            "max_forward_substitution_pa": self.max_forward_substitution_pa,
            "mean_evaluations_per_solve": (
                self.chain_evaluations / self.solves if self.solves else 0.0
            ),
        }

    def declaration(self) -> dict[str, Any]:
        return {
            "record": "GT-PS-2-W32-GS-SWEEP-STUDY-LEVER-2026-09-03-01",
            "ruling": "W3-2",
            "finding": "GT_PS2_WAVE3_SOLVER_STRUGGLE_AUDIT_2026-09-03.md #3.2",
            "relaxation": (
                "cascade-ordered (Gauss-Seidel) sweep of the tower face-pressure "
                "Picard, by exact forward substitution of the just-computed upper "
                "face into the lower face update"
            ),
            "fixed_point": "identical by construction (see module docstring)",
            "acceptance_test": "unchanged (same residual, normalizers and limit)",
            "forward_substitution_coefficient": FORWARD_SUBSTITUTION_COEFFICIENT,
            "patch_target": PATCH_TARGET,
            "production_picard_source_digest": PRODUCTION_PICARD_SOURCE_DIGEST,
            "derived_spectral_radius_jacobi": 0.079418,
            "derived_spectral_radius_sweep": 0.014545,
        }

    # -- the replacement ----------------------------------------------------
    def _build(self) -> Any:
        """The production loop, with exactly one statement changed.

        Every helper, constant and exception is reached through the mechanism
        module at call time, so the copy tracks the original everywhere except
        at the iterate update; the source digest asserted in :meth:`install`
        is what makes "everywhere except" a checked claim rather than a hope.
        """

        lever = self
        module = _mechanism_module()

        def _close_through_pressures(
            request: Any,
            *,
            codec: Any,
            planners: dict[str, Any],
            bundles: dict[str, Any],
            bundle_digests: dict[str, str] | None = None,
            work_context: Any,
        ) -> tuple[Any, Any]:
            lever.solves += 1
            face_pressures = (
                request.through_bed_states[1].cell_field_states[0].layer_pressure_pa,
                request.through_bed_states[0].cell_field_states[0].layer_pressure_pa,
            )
            prepared_cache: dict[str, Any] | None = (
                {} if module.through_bed.HOIST_ENABLED else None
            )
            low, high = module.ef.ENGINEERING_SIX_LAYER_SERIAL_LAW2_PRESSURE_ADMISSION_PA
            for evaluation in range(1, request.pressure_root_max_evaluations + 1):
                sp1_mn2, mn2_mn1 = face_pressures
                if (
                    not low <= mn2_mn1 < sp1_mn2 <= high
                    or not module.MN1_TOP_PRESSURE_PA < mn2_mn1
                ):
                    lever.domain_refusals += 1
                    raise module.QSCK2MechanismIntervalError(
                        "pressure Picard iterate left the ordered structural domain"
                    )
                lever.chain_evaluations += 1
                transactions = module._evaluate_through_chain(
                    request,
                    sp1_mn2_pressure_pa=sp1_mn2,
                    mn2_mn1_pressure_pa=mn2_mn1,
                    codec=codec,
                    planners=planners,
                    bundles=bundles,
                    bundle_digests=bundle_digests,
                    work_context=work_context,
                    prepared_cache=prepared_cache,
                )
                mn1, mn2, _sp1 = transactions
                proposed = (
                    mn2.accepted_fast_steps[0].state.layer_pressure_pa,
                    mn1.accepted_fast_steps[0].state.layer_pressure_pa,
                )
                residuals = (proposed[0] - sp1_mn2, proposed[1] - mn2_mn1)
                relative = max(
                    abs(residuals[0]) / max(abs(proposed[0]), abs(sp1_mn2), 1.0),
                    abs(residuals[1]) / max(abs(proposed[1]), abs(mn2_mn1), 1.0),
                )
                if relative <= request.pressure_face_relative_limit:
                    evidence = module.QSCMechanismPressureSolveEvidence(
                        evaluations=evaluation,
                        sp1_mn2_face_pressure_pa=sp1_mn2,
                        mn2_mn1_face_pressure_pa=mn2_mn1,
                        sp1_mn2_residual_pa=residuals[0],
                        mn2_mn1_residual_pa=residuals[1],
                        maximum_relative_residual=relative,
                        relative_limit=request.pressure_face_relative_limit,
                    )
                    return transactions, evidence
                # ---- W3-2: the ONE changed statement -----------------------
                # Production: ``face_pressures = proposed`` (Jacobi).  Here the
                # upper face is taken as just computed and forward-substituted
                # into the lower face, which is the cascade's own sweep order.
                correction = FORWARD_SUBSTITUTION_COEFFICIENT * residuals[1]
                if abs(correction) > lever.max_forward_substitution_pa:
                    lever.max_forward_substitution_pa = abs(correction)
                lever.sweep_updates += 1
                face_pressures = (proposed[0] + correction, proposed[1])
                # ------------------------------------------------------------
            lever.evaluation_limit_refusals += 1
            raise module.QSCK2MechanismIntervalError(
                "pressure Picard iteration exceeded its evaluation limit"
            )

        return _close_through_pressures
