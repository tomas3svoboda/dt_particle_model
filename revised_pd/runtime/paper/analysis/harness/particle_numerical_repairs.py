"""Explicit numerical coordinates and pure-bulk film closure error budgets.

Physical residuals and stored states are unchanged. Endpoint closure uses a
prospective 128-epsilon arithmetic-scale budget, not a composition clamp.
"""

from contextlib import contextmanager, nullcontext
import math
import os
import sys
import time
from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded
import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize._numdiff import group_columns
from dtdc_simulator.core2.particle import external_film_reduced as efr


# Owner ruling GT_PS2_OWNER_RULING_PAPER1_NUMERICAL_GATE_HYGIENE_2026-10-17 item 3:
# the dense certificate Jacobian taken at an iterate that already meets gtol
# (it owns the native condition gate) is counted outside max_nfev, the search
# budget, which is unchanged. False restores the former counting for replay.
CERTIFICATE_JACOBIAN_OUTSIDE_SEARCH_BUDGET = True

# Lead decision (Class A) under the hygiene principle of
# GT_PS2_OWNER_RULING_PAPER1_NUMERICAL_GATE_HYGIENE_2026-10-17, recorded in
# docs/GT_PS2_PAPER1_LINE_SEARCH_MERIT_BUILD_2026-10-18.md (evidence 149, 150):
# with a per-row acceptance clause, the line search's sufficient-decrease test
# is taken on the gate-normalized residual r_i / g_i, g the clause's own row
# gates at the current iterate from the Jacobian the Newton already holds
# there (no residual call is added), so a step is judged in the units that
# accept the iterate and rounding noise below a row's gate cannot refuse it.
# The norm of the normalized vector is declared here ("2" or "max"); without
# an acceptance clause the former unweighted 2-norm is used, bit for bit.
GATE_NORMALIZED_MERIT_NORM = "max"
GATE_NORMALIZED_MERIT_NORMS = ("2", "max")

# Kernel-envelope JACOBIAN_BROYDEN_UPDATES and its companion JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES
# (the lead's Class-A numerics of 2026-10-31, record JACOBIAN_LEVERS_RECORD, measured in
# JACOBIAN_LEVERS_MEASUREMENT sections 3.2, 3.3 and 5). Armed only, and only inside a
# JacobianLeverStage scope: the front driver enters one for its front stage, the runner one for its
# dry stage; everywhere else, and with the knobs absent, domain_aware_newton runs its former
# statements, bit for bit. Armed, the search Jacobian is held across the iterations of a solve and
# updated after every accepted step by Broyden's rank-one formula J + (dr - J s) s^T / (s^T s) on
# the unscaled step s and residual change dr; it is rebuilt as a fresh dense difference Jacobian
# at the next iterate after an accepted step whose merit ratio merit_after / merit_before exceeds
# the declared contraction ratio below, and at the same iterate when its direction finds no
# admissible decreasing step (the original step and failure test are then taken on the fresh
# Jacobian). With the companion, in the dry stage only, a solve starts from the previous converged
# solve's dense certificate Jacobian of the same layout (the calling site, the unknown count, the
# acceptance class and the difference mode), then Broyden. Unchanged: the convergence test (gtol
# and the acceptance clause), max_nfev and its counting, the dense certificate Jacobian at every
# returned or provisional iterate, the line search and its merit, the damping ladder and the bounds
# handling. Each solve is stamped (JacobianLeverStage.drain): the Jacobian source of every
# iteration (fresh, fresh_certificate, carried, broyden_updated), the rebuild and update counts.
# The ratio is the one measured (no sweep; the measurement's section 6): an accepted step that
# removes less than half of the merit marks the held Jacobian as stale.
BROYDEN_REBUILD_CONTRACTION_RATIO = 0.5  # WAS absent until 2026-10-31 (no Broyden update)
JACOBIAN_LEVERS_RECORD = "docs/GT_PS2_DTDC_KERNEL_JACOBIAN_LEVERS_BUILD_2026-10-31.md"
JACOBIAN_LEVERS_MEASUREMENT = "docs/GT_PS2_DTDC_KERNEL_JACOBIAN_LEVER_MEASURED_2026-10-31.md"
JACOBIAN_LEVER_STAGES = ("front", "dry", "birth")  # WAS ("front", "dry") until 2026-11-02
# Kernel-envelope JACOBIAN_LEVERS_AT_BIRTH (the companion of JACOBIAN_BROYDEN_UPDATES; the lead's
# Class-A numerics of 2026-11-02, record BIRTH_COST_LEVERS_RECORD): the "birth" stage, the birth's
# one domain-aware Newton solve (the Paper-1 runner's fresh_dt_birth and the cohort tower's born
# job enter it), takes the Broyden lever exactly as the front stage does; the carried companion is
# refused there (a birth is one solve). With the knob absent jacobian_lever_stage("birth") is None
# and the birth runs its former statements.
BIRTH_COST_LEVERS_RECORD = "docs/GT_PS2_DTDC_KERNEL_BIRTH_COST_LEVERS_2026-11-02.md"
# The stage scope in force; None outside every scope (the former Newton).
_ACTIVE_JACOBIAN_LEVERS = None
# Kernel-envelope JACOBIAN_COLUMN_GROUPS (the lead's Class-A numerics of 2026-11-02, record
# particle_jacobian_column_groups.RECORD; a VALUED knob, armed only): a declared site passes its
# declared structure as column_groups (a particle_jacobian_column_groups.ColumnGroupPlan) and the
# solve takes its finite-difference Jacobians by column groups, exact bit for bit where the structure
# holds (the same step, sign, subtraction and division per column, the last probe the former last
# probe); a probe moving a row outside the structure falls back to the one-at-a-time Jacobian and
# disables the layout; the first solves of a process are bit-checked. The search Jacobian of the
# dry-rowsum value assembles the dry stage's dense rows from the inventory terms (not bit-identical;
# the certificate stays exact). Without column_groups (every caller with the knob absent) the solve
# runs its former statements, bit for bit.
# Kernel-envelope DRY_FACE_FLUX_MEMO (particle_face_flux_memo, armed only): the installed memo, whose
# frame scopes each domain-aware Newton solve; None (the knob absent) runs the solve as before.
_FACE_FLUX_MEMO = None
# Kernel-envelope DEFERRED_CERTIFICATE (a VALUED knob armed only by exactly "1"; record
# DEFERRED_CERTIFICATE_RECORD; the lead's Class-A numerics, a measurement, nothing adopted). Inside a
# DeferredCertificateScope (the front driver's run enters one when the knob is armed), a solve called
# with deferrable=True whose iterate meets gtol returns that iterate at once with its dense
# certificate pending (result.jac None, result.deferred_certificate a PendingCertificate). The
# certificate is the same dense finite-difference Jacobian at the same iterate by the same
# statements, computed later: at the next deferrable solve's certificate point (the in-process
# stand-in for a certifier working while the driver's next step searches) or at a drain the driver
# asks for, inside the residual's process-global hooks as they were at the return (restored for its
# columns, then put back), outside max_nfev as the hygiene ruling's item 3 counts it, its proof
# points charged to the step's own allowance. Its consumers (the film node's condition gate and the
# condition number the accepted step's ledger stores) run when it is in. A certificate whose
# consumers refuse after the driver moved on raises DeferredCertificateRollback and the run is
# replayed with that certificate taken at once (today's path). Not deferred in this build: an
# iterate returned on an acceptance clause's candidate test, a solve whose certificate is charged
# to max_nfev (CERTIFICATE_JACOBIAN_OUTSIDE_SEARCH_BUDGET False), a carried-Jacobian scope. Without
# the scope (the knob absent) every solve runs its former statements, bit for bit.
DEFERRED_CERTIFICATE_RECORD = "docs/GT_PS2_DTDC_KERNEL_DEFERRED_CERTIFICATE_PROOF_2026-11-03.md"
_DEFERRED_CERTIFICATES = None
# The process-global hooks a deferred certificate's residual may read: the module attributes and
# class attributes of the repository's own modules (the harness and dtdc_simulator.core2) that are
# callables or classes (the scopes of a step install their closures there). Others are left alone
# and, when they differ at the certificate, stamped.
_HOOK_ROOTS = (
    os.path.normcase(os.path.dirname(os.path.abspath(__file__))),
    os.path.normcase(
        os.path.join(
            os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            ),
            "src",
            "dtdc_simulator",
            "core2",
        )
    ),
)


class CertificateChargeMismatch(RuntimeError):
    """A deferred certificate took another number of residual calls than was charged at its return
    (a retreating probe): the run is replayed with it taken at once."""


class DeferredCertificateRollback(BaseException):
    """A deferred certificate refused after the driver moved on; the run replays it at once.

    A BaseException, so that no step handler (``except Exception``) takes it for its own refusal.
    """

    def __init__(self, index, cause):
        super().__init__(f"deferred certificate {index} refused after the driver moved on: {cause}")
        self.index = index
        self.cause = cause


def _hook_snapshot():
    entries = []
    for module in list(sys.modules.values()):
        path = getattr(module, "__file__", None)
        if not path:
            continue
        path = os.path.normcase(os.path.abspath(path))
        if not path.startswith(_HOOK_ROOTS):
            continue
        name = getattr(module, "__name__", None)
        for key, value in list(vars(module).items()):
            entries.append((module, key, value))
            if isinstance(value, type) and getattr(value, "__module__", None) == name:
                for ckey, cvalue in list(vars(value).items()):
                    if ckey.startswith("__") and ckey.endswith("__"):
                        continue
                    entries.append((value, ckey, cvalue))
    return entries


_MISSING = object()


@contextmanager
def _restored_hooks(entries, stamp):
    swapped = []
    for owner, key, then in entries:
        now = vars(owner).get(key, _MISSING)
        if now is then:
            continue
        if callable(then) or isinstance(then, (staticmethod, classmethod, property)):
            swapped.append((owner, key, now))
            setattr(owner, key, then)
        else:
            stamp.setdefault("differing_noncallables_left", set()).add(
                f"{getattr(owner, '__name__', owner)}.{key}"
            )
    stamp["hooks_restored"] = sorted(
        {f"{getattr(owner, '__name__', owner)}.{key}" for owner, key, _ in swapped}
    )
    try:
        yield
    finally:
        for owner, key, now in reversed(swapped):
            if now is _MISSING:
                delattr(owner, key)
            else:
                setattr(owner, key, now)


class PendingCertificate:
    """A returned iterate's dense certificate, not yet computed (DEFERRED_CERTIFICATE)."""

    def __init__(self, index, compute, hooks):
        self.index = index
        self._compute = compute
        self._hooks = hooks
        self.after = []  # callables(jac) after the certificate (physical_chart_solve's pullback)
        self.verdict = None  # the consumer's gate and fills (the film node), set by the consumer
        self.done = False
        self.jac = None
        self.stamp = {"index": index}

    def certify(self):
        if self.done:
            return self.jac
        started = time.perf_counter()
        memo = _FACE_FLUX_MEMO
        with (
            _restored_hooks(self._hooks, self.stamp),
            memo.frame() if memo is not None else nullcontext(),
        ):
            jac = self._compute()
            for hook in self.after:
                hook(jac)
        self.stamp["certificate_wall_s"] = time.perf_counter() - started
        self.done = True
        self.jac = jac
        return jac


class DeferredCertificateScope:
    """The deferred certificates of one run (DEFERRED_CERTIFICATE): their order, drains, stamps.

    eager: the indices (in deferral order) taken at once, as today; the replay after a rollback.
    """

    def __init__(self, eager=()):
        self.eager = frozenset(eager)
        self.count = 0
        self.pending = []
        self.rolling_back = False
        self.on_drained = []
        self.by_result = {}  # id(result) -> (result, its PendingCertificate)
        self.watchers = []  # callables returning the containers a sentinel may have been copied to
        # per-step counter dicts (the driver's): a late certificate's increments are taken out of
        # the step then open and kept on the certificate for its own step's row
        self.step_counters = []
        self.resolved = []
        self.totals = {
            "certificate_points": 0,
            "deferred": 0,
            "taken_at_once_by_replay": 0,
            "resolved_inside_their_step": 0,
            "certified_late": 0,
            "certificate_wall_s": 0.0,
            "drains": 0,
            "hooks_restored": [],
            "differing_noncallables_left": [],
        }

    @contextmanager
    def active(self):
        global _DEFERRED_CERTIFICATES
        previous = _DEFERRED_CERTIFICATES
        _DEFERRED_CERTIFICATES = self
        try:
            yield self
        finally:
            _DEFERRED_CERTIFICATES = previous

    def register(self, result, pending):
        self.pending.append(pending)
        self.by_result[id(result)] = (result, pending)

    def pending_for(self, result):
        entry = self.by_result.get(id(result))
        return None if entry is None or entry[0] is not result else entry[1]

    def next_index(self):
        index = self.count
        self.count += 1
        self.totals["certificate_points"] += 1
        return index

    def defers(self, index):
        if index in self.eager:
            self.totals["taken_at_once_by_replay"] += 1
            return False
        self.totals["deferred"] += 1
        return True

    def _account(self, pending):
        self.totals["certificate_wall_s"] += pending.stamp.get("certificate_wall_s", 0.0)
        for key in ("hooks_restored", "differing_noncallables_left"):
            for name in sorted(pending.stamp.get(key, ())):
                if name not in self.totals[key]:
                    self.totals[key].append(name)

    def resolve_now(self, pending, *, verdict):
        """Inside its own step (an exception after the return, or a consumer needing the matrix
        now): the certificate at once and, when the step had reached its gate, the verdict, whose
        refusal is the step's own error as today."""
        if pending in self.pending:
            self.pending.remove(pending)
        try:
            pending.certify()
        except CertificateChargeMismatch as error:
            self.rolling_back = True
            raise DeferredCertificateRollback(pending.index, f"{type(error).__name__}: {error}")
        self._account(pending)
        self.totals["resolved_inside_their_step"] += 1
        if verdict and pending.verdict is not None:
            pending.verdict()

    def drain(self):
        """Every pending certificate in order; a refusal now is a rollback of the run."""
        if self.rolling_back or not self.pending:
            return
        self.totals["drains"] += 1
        while self.pending:
            pending = self.pending.pop(0)
            before = [dict(counter) for counter in self.step_counters]
            try:
                pending.certify()
                if pending.verdict is not None:
                    pending.verdict()
            except Exception as error:
                self.rolling_back = True
                raise DeferredCertificateRollback(
                    pending.index, f"{type(error).__name__}: {error}"
                ) from error
            pending.counter_deltas = []
            for counter, then in zip(self.step_counters, before):
                delta = {key: counter[key] - then.get(key, 0) for key in counter}
                for key, change in delta.items():
                    counter[key] -= change
                pending.counter_deltas.append(delta)
            self.resolved.append(pending)
            self._account(pending)
            self.totals["certified_late"] += 1
        for callback in self.on_drained:
            callback()

    def watched(self):
        containers = []
        for watcher in self.watchers:
            containers.extend(watcher())
        return containers

    def fill(self, sentinel, value, containers):
        """Replace the sentinel (bit for bit) by the certificate's value in the containers."""
        for container in containers:
            _replace_float(container, sentinel, value)


def _replace_float(node, sentinel, value):
    if isinstance(node, dict):
        for key, item in node.items():
            if isinstance(item, float) and item == sentinel:
                node[key] = value
            elif isinstance(item, (dict, list)):
                _replace_float(item, sentinel, value)
    elif isinstance(node, list):
        for position, item in enumerate(node):
            if isinstance(item, float) and item == sentinel:
                node[position] = value
            elif isinstance(item, (dict, list)):
                _replace_float(item, sentinel, value)


def deferrable_kwargs():
    """{"deferrable": True} inside a DEFERRED_CERTIFICATE scope; {} otherwise (the former call)."""
    return {} if _DEFERRED_CERTIFICATES is None else {"deferrable": True}


def pending_certificate(result):
    """The result's pending certificate inside a DEFERRED_CERTIFICATE scope; None otherwise."""
    scope = _DEFERRED_CERTIFICATES
    return None if scope is None else scope.pending_for(result)


def deferred_certificate_sentinel(index):
    """A float no certificate value takes (negative, about -1e301), unique per deferral."""
    return -math.ldexp(1.0, 1000) * (1.0 + index * 2.0**-30)


class JacobianLeverStage:
    """The armed Jacobian levers of one stage of one run, its carried Jacobians and its stamps.

    Built by jacobian_lever_stage (None when JACOBIAN_BROYDEN_UPDATES is not armed). Every
    domain_aware_newton call made inside active() runs with the levers; drain() returns the
    per-solve stamps since the previous drain (one history row's stamp); totals counts the stage.
    """

    def __init__(self, stage, *, carried_across_advances):
        if stage not in JACOBIAN_LEVER_STAGES:
            raise ValueError(f"undeclared Jacobian lever stage {stage!r}")
        if carried_across_advances and stage != "dry":
            raise ValueError("the carried Jacobian is declared for the dry stage only")
        self.stage = stage
        self.carried_across_advances = bool(carried_across_advances)
        self.carried = {}
        self.pending = []
        self.totals = {
            "solves": 0,
            "solves_not_converged": 0,
            "iterations": 0,
            "fresh_search_jacobians": 0,
            "fresh_certificate_jacobians": 0,
            "carried_jacobians_used": 0,
            "broyden_updated_jacobians_used": 0,
            "broyden_updates": 0,
            "rebuilds_after_stalled_contraction": 0,
            "rebuilds_after_failed_direction": 0,
        }

    @contextmanager
    def active(self):
        global _ACTIVE_JACOBIAN_LEVERS
        previous = _ACTIVE_JACOBIAN_LEVERS
        _ACTIVE_JACOBIAN_LEVERS = self
        try:
            yield self
        finally:
            _ACTIVE_JACOBIAN_LEVERS = previous

    def declared(self):
        """The stage's declared settings (the history record's stamp)."""
        return {
            "stage": self.stage,
            "broyden_updates": True,
            "carried_across_advances": self.carried_across_advances,
            "rebuild_contraction_ratio": BROYDEN_REBUILD_CONTRACTION_RATIO,
            "record": JACOBIAN_LEVERS_RECORD,
            "physically_qualifying": False,
        }

    def drain(self):
        """Every solve's stamp since the previous drain (refused attempts included), then reset."""
        solves, self.pending = self.pending, []
        return {"stage": self.stage, "solves": solves}

    def _open(self, site, unknowns):
        stamp = {
            "site": site,
            "unknowns": unknowns,
            "jacobian_sources": [],
            "broyden_updates": 0,
            "rebuilds": 0,
            "rebuilds_after_stalled_contraction": 0,
            "rebuilds_after_failed_direction": 0,
            "outcome": None,
        }
        self.pending.append(stamp)
        self.totals["solves"] += 1
        return stamp

    def _source(self, stamp, source):
        stamp["jacobian_sources"].append(source)
        self.totals["iterations"] += 1
        self.totals[
            {
                "fresh": "fresh_search_jacobians",
                "fresh_certificate": "fresh_certificate_jacobians",
                "carried": "carried_jacobians_used",
                "broyden_updated": "broyden_updated_jacobians_used",
            }[source]
        ] += 1

    def _count(self, stamp, name):
        stamp[name] += 1
        if name.startswith("rebuilds_"):
            stamp["rebuilds"] += 1
        self.totals[name] += 1

    def _close(self, stamp, outcome):
        stamp["outcome"] = outcome
        if outcome != "converged":
            self.totals["solves_not_converged"] += 1


def jacobian_lever_stage(stage, environ=None, *, carried_across_advances=None):
    """The stage's armed Jacobian levers, or None when JACOBIAN_BROYDEN_UPDATES is not armed.

    The carried Jacobian follows its knob in the dry stage. It is refused for the front stage:
    asking for it there (carried_across_advances=True) raises; there the carried matrix is a poorer
    first Jacobian than a fresh one and costs wall (the measurement's section 3.3).
    """
    import particle_kernel_envelope as envelope

    if stage not in JACOBIAN_LEVER_STAGES:
        raise ValueError(f"undeclared Jacobian lever stage {stage!r}")
    knobs = envelope.validate_environment(environ)
    if carried_across_advances and stage != "dry":
        raise envelope.EnvelopeKnobRefused(
            f"{envelope.PREFIX + envelope.JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES} is refused for the "
            f"{stage} stage: the carried Jacobian is declared for the dry stage only (in the front "
            "the front coordinate and the cut cell move between advances; "
            f"{JACOBIAN_LEVERS_MEASUREMENT} section 3.3)"
        )
    if not knobs[envelope.JACOBIAN_BROYDEN_UPDATES]:
        return None
    if stage == "birth" and not knobs[envelope.JACOBIAN_LEVERS_AT_BIRTH]:
        # JACOBIAN_LEVERS_AT_BIRTH absent: the birth keeps its former Newton
        return None
    carried = knobs[envelope.JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES] and stage == "dry"
    if carried_across_advances is not None:
        carried = bool(carried_across_advances) and carried
    return JacobianLeverStage(stage, carried_across_advances=carried)


def _jacobian_lever_layout(size, acceptance, difference_mode):
    """The solve's calling site outside this module, and its layout key."""
    frame = sys._getframe(1)
    while frame is not None and frame.f_code.co_filename.endswith("particle_numerical_repairs.py"):
        frame = frame.f_back
    where = (
        ("?", "?", 0)
        if frame is None
        else (os.path.basename(frame.f_code.co_filename), frame.f_code.co_name, frame.f_lineno)
    )
    return "%s:%s:%d" % where, (*where, int(size), type(acceptance).__name__, difference_mode)


def gate_normalized_merit(residual, gates, norm=None):
    """The line-search merit in gate units: the norm of residual_i / gate_i."""
    norm = GATE_NORMALIZED_MERIT_NORM if norm is None else norm
    if norm not in GATE_NORMALIZED_MERIT_NORMS:
        raise ValueError(f"undeclared gate-normalized merit norm {norm!r}")
    scaled = np.asarray(residual, dtype=float) / gates
    return float(np.linalg.norm(scaled) if norm == "2" else np.max(np.abs(scaled)))


def merit_row_gates(acceptance, x, jac, size):
    """The acceptance clause's row gates at x (its gates(x, jac)), validated."""
    gates = np.asarray(acceptance.gates(x, jac)[0], dtype=float)
    if gates.shape != (size,) or not np.all(np.isfinite(gates)) or not np.all(gates > 0):
        raise ValueError("acceptance row gates must be positive, finite and one per residual row")
    return gates


def search_column_scale(jac):
    """The column 2-norm scaling of the domain-aware Newton's least-squares step.

    domain_aware_newton solves (jac * search_column_scale(jac)) for its step;
    owner ruling GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17 point 3
    takes the birth condition certificate on exactly this matrix.
    """
    return 1 / np.maximum(np.linalg.norm(jac, axis=0), np.finfo(float).tiny)


def searched_matrix_condition(jac):
    """np.linalg.cond of the column-2-norm-scaled matrix the Newton step solves."""
    matrix = np.asarray(jac, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return math.inf
    return float(np.linalg.cond(matrix * search_column_scale(matrix)))


def domain_aware_newton(fun, x0, **kwargs):
    """The domain-aware Newton (_domain_aware_newton, whose docstring is the contract).

    Kernel-envelope DRY_FACE_FLUX_MEMO (armed only): the solve runs inside its own frame of the
    installed face-flux memo (particle_face_flux_memo); absent, the solve itself, unchanged.
    """
    face_memo = _FACE_FLUX_MEMO
    if face_memo is None:
        return _domain_aware_newton(fun, x0, **kwargs)
    with face_memo.frame():
        return _domain_aware_newton(fun, x0, **kwargs)


def _domain_aware_newton(
    fun, x0, *, bounds, max_nfev=400, gtol=1e-11, difference_mode="forward",
    search_jacobian_sparsity=None, acceptance=None, column_groups=None, deferrable=False,
    **unused
):
    """Backtrack undefined physical trials; never replace their residuals.

    Jacobians use finite differences with one-sided retreat when a probe is
    outside the physical domain. Every real residual call counts in nfev.

    acceptance (default None: the former code path, bit for bit) is a declared
    per-row acceptance clause whose row gates are never tighter than gtol:
    an object with candidate(residual) (every row with a fixed gate passes it)
    and accepts(x, residual, jac) (every row passes its gate, some gates
    depending on the dense Jacobian at x). An iterate it accepts is returned
    as converged although its maximum exceeds gtol. The dense Jacobian at a
    candidate iterate is its certificate: counted outside max_nfev when the
    iterate is returned (the hygiene ruling's item 3), charged to max_nfev
    when the clause then refuses it (it is the next step's search Jacobian).
    With a clause, the line search's sufficient-decrease merit is the
    gate-normalized residual (gate_normalized_merit with the clause's
    gates(x, jac) at the current iterate, from the Jacobian already held there:
    no residual call is added); the merit and each accepted step's damping,
    length and retreat count are returned as line_search. Without a clause the
    merit is the former unweighted 2-norm and the result carries no new key.

    Inside a JacobianLeverStage scope (kernel-envelope JACOBIAN_BROYDEN_UPDATES, armed only) the
    search Jacobian is held and Broyden-updated between rebuilds; see that section above. The
    result's keys are unchanged.
    """
    # Kernel-envelope JACOBIAN_BROYDEN_UPDATES (inside a JacobianLeverStage scope only): the held
    # search Jacobian, its source and the solve's stamp; None outside every scope (the former Newton).
    levers = _ACTIVE_JACOBIAN_LEVERS
    stamp = held = layout = held_source = None
    force_fresh = False
    if levers is not None:
        site, layout = _jacobian_lever_layout(np.size(x0), acceptance, difference_mode)
        stamp = levers._open(site, int(np.size(x0)))
        carried = levers.carried.get(layout) if levers.carried_across_advances else None
        if carried is not None and carried.shape == (np.size(x0), np.size(x0)):
            held = carried.copy()
            held_source = "carried"
    x = np.array(x0, dtype=float, copy=True)
    lower, upper = map(np.asarray, bounds)
    work = 0
    certificate = {"active": False, "calls": 0}
    rejected = 0
    rejection_causes = {}
    if difference_mode not in ("forward", "central"):
        raise ValueError("unknown numerical derivative mode")
    pattern = None
    if search_jacobian_sparsity is not None:
        if difference_mode != "forward":
            raise ValueError("grouped search currently requires forward differences")
        pattern = np.asarray(search_jacobian_sparsity, dtype=bool)
        if pattern.ndim != 2 or pattern.shape[1] != x.size or not np.all(np.any(pattern, axis=0)):
            raise ValueError("invalid declared search Jacobian structure")
        groups = group_columns(pattern)
        for group in np.unique(groups):
            if np.any(np.sum(pattern[:, groups == group], axis=1) > 1):
                raise ValueError("declared Jacobian colors collide in a residual row")
    # Kernel-envelope JACOBIAN_COLUMN_GROUPS (armed only; see the declaration above)
    cg_stamp = None
    if column_groups is not None:
        import particle_jacobian_column_groups as cgj

        if pattern is not None or difference_mode != "forward":
            raise ValueError("column groups take forward differences and no grouped search")
        cg_stamp = cgj.open_solve(column_groups)

    def evaluate(values):
        nonlocal work, rejected
        if certificate["active"]:
            certificate["calls"] += 1
        else:
            if work >= max_nfev:
                raise RuntimeError("domain-aware Newton residual-call budget exhausted")
            work += 1
        try:
            result = np.asarray(fun(values), dtype=float)
            if not np.all(np.isfinite(result)):
                raise ValueError("physical residual is nonfinite")
            return result
        except StepBudgetExceeded:
            raise
        except (ValueError, RuntimeError) as exc:
            if "budget" in str(exc).lower() or "source" in str(exc).lower():
                raise
            rejected += 1
            label = f"{type(exc).__name__}: {exc}"
            rejection_causes[label] = rejection_causes.get(label, 0) + 1
            raise

    def jacobian(values, residual):
        columns = []
        for column in range(values.size):
            h = (
                np.cbrt(np.finfo(float).eps)
                if difference_mode == "central"
                else math.sqrt(np.finfo(float).eps)
            ) * max(1.0, abs(values[column]))
            for retreat in range(16):
                found = False
                samples = []
                for sign in (1.0, -1.0):
                    trial = values.copy()
                    trial[column] += sign * h
                    if not lower[column] < trial[column] < upper[column]:
                        continue
                    try:
                        sample = evaluate(trial)
                    except StepBudgetExceeded:
                        raise
                    except (ValueError, RuntimeError) as exc:
                        if "budget" in str(exc).lower() or "source" in str(exc).lower():
                            raise
                        continue
                    represented_h = trial[column] - values[column]
                    if represented_h == 0:
                        continue
                    samples.append((represented_h, sample))
                    if difference_mode == "forward":
                        break
                if len(samples) == 2:
                    columns.append(
                        (samples[0][1] - samples[1][1]) / (samples[0][0] - samples[1][0])
                    )
                    found = True
                elif len(samples) == 1 and (difference_mode == "forward" or retreat >= 4):
                    columns.append((samples[0][1] - residual) / samples[0][0])
                    found = True
                if found:
                    break
                h *= 0.5
            else:
                raise RuntimeError("no admissible one-sided physical Jacobian probe")
        return np.column_stack(columns)

    def one_column(values, residual, column, call=None):
        # Kernel-envelope JACOBIAN_COLUMN_GROUPS: one column by the one-at-a-time loop's statements
        call = evaluate if call is None else call
        h = (
            np.cbrt(np.finfo(float).eps)
            if difference_mode == "central"
            else math.sqrt(np.finfo(float).eps)
        ) * max(1.0, abs(values[column]))
        for retreat in range(16):
            samples = []
            for sign in (1.0, -1.0):
                trial = values.copy()
                trial[column] += sign * h
                if not lower[column] < trial[column] < upper[column]:
                    continue
                try:
                    sample = call(trial)
                except StepBudgetExceeded:
                    raise
                except (ValueError, RuntimeError) as exc:
                    if "budget" in str(exc).lower() or "source" in str(exc).lower():
                        raise
                    continue
                represented_h = trial[column] - values[column]
                if represented_h == 0:
                    continue
                samples.append((represented_h, sample))
                if difference_mode == "forward":
                    break
            if len(samples) == 2:
                return (samples[0][1] - samples[1][1]) / (samples[0][0] - samples[1][0])
            if len(samples) == 1 and (difference_mode == "forward" or retreat >= 4):
                return (samples[0][1] - residual) / samples[0][0]
            h *= 0.5
        raise RuntimeError("no admissible one-sided physical Jacobian probe")

    def unchecked_residual(values):
        # Kernel-envelope JACOBIAN_COLUMN_GROUPS' bit check: a residual call outside the solve's count
        result = np.asarray(fun(values), dtype=float)
        if not np.all(np.isfinite(result)):
            raise ValueError("physical residual is nonfinite")
        return result

    def column_group_charge(extra):
        # Kernel-envelope JACOBIAN_COLUMN_GROUPS: a group probe's other columns booked as the former
        # loop books them, each call checked against the residual-call budget in turn, then the
        # residual's own counters (the site's residual_charge), then counted
        nonlocal work
        for _ in range(extra):
            if not certificate["active"] and work >= max_nfev:
                raise RuntimeError("domain-aware Newton residual-call budget exhausted")
            if column_groups.residual_charge is not None:
                column_groups.residual_charge(1)
            if certificate["active"]:
                certificate["calls"] += 1
            else:
                work += 1

    def column_group_jacobian(values, residual, search):
        # Kernel-envelope JACOBIAN_COLUMN_GROUPS (armed only; see the declaration above)
        try:
            matrix = cgj.grouped_jacobian(
                column_groups,
                cg_stamp,
                values,
                residual,
                lower,
                upper,
                evaluate,
                one_column,
                search=search,
                charge=column_group_charge,
            )
        except cgj._StructureViolation as violation:
            cgj.violation(column_groups, cg_stamp, violation.args[0])
            cgj._count(cg_stamp, "dense_fallback_jacobians")
            return jacobian(values, residual)
        cgj._count(cg_stamp, "grouped_jacobians")
        if cg_stamp["bit_checked"]:
            uncounted = column_groups.residual_uncounted
            with uncounted() if uncounted is not None else nullcontext():
                dense = np.column_stack(
                    [
                        one_column(values, residual, column, unchecked_residual)
                        for column in range(values.size)
                    ]
                )
            cgj._count(cg_stamp, "bit_checks")
            exact = np.ones(matrix.shape, dtype=bool)
            if search and column_groups.rowsum is not None:
                rows = list(column_groups.rowsum.dense_rows)
                exact[np.ix_(rows, list(column_groups.rowsum.far_columns))] = False
                assembled = matrix[~exact]
                former = dense[~exact]
                cg_stamp["rowsum"][-1]["checked_maximum_relative_difference"] = (
                    float(np.max(np.abs(assembled - former) / np.maximum(np.abs(former), 1e-300)))
                    if assembled.size and cg_stamp["rowsum"][-1]["assembled"]
                    else 0.0
                )
            if matrix[exact].tobytes() != dense[exact].tobytes():
                raise cgj.ColumnGroupJacobianMismatch(
                    f"{column_groups.site}: grouped and one-at-a-time Jacobians differ at "
                    f"{np.argwhere((matrix != dense) & exact)[:8].tolist()}"
                )
            cg_stamp["bit_checks_equal"] += 1
        return matrix

    def grouped_search_jacobian(values, residual):
        matrix = np.zeros(pattern.shape)
        for group in np.unique(groups):
            columns = np.flatnonzero(groups == group)
            trial = values.copy()
            steps = math.sqrt(np.finfo(float).eps) * np.maximum(1.0, np.abs(values[columns]))
            signs = np.where(values[columns] + steps < upper[columns], 1.0, -1.0)
            trial[columns] += signs * steps
            if np.any(trial[columns] <= lower[columns]) or np.any(trial[columns] >= upper[columns]):
                return jacobian(values, residual)
            try:
                sample = evaluate(trial)
            except StepBudgetExceeded:
                raise
            except (ValueError, RuntimeError) as exc:
                if "budget" in str(exc).lower() or "source" in str(exc).lower():
                    raise
                return jacobian(values, residual)
            for column in columns:
                rows = pattern[:, column]
                matrix[rows, column] = (sample[rows] - residual[rows]) / (trial[column] - values[column])
        return matrix

    residual = evaluate(x)
    best = (float(np.max(np.abs(residual))), x.copy())
    line_search = None
    if acceptance is not None:
        line_search = {
            "merit": "gate-normalized residual r_i / g_i, g the acceptance clause's row gates "
            "at the current iterate from the Newton's Jacobian there",
            "norm": GATE_NORMALIZED_MERIT_NORM,
            "merit_residual_calls": 0,
            "steps": [],
        }
    try:
        for iteration in range(max_nfev):
            maximum = float(np.max(np.abs(residual)))
            # Dense original-coordinate derivatives own final conditioning.
            # Grouping only accelerates intermediate Newton search directions.
            provisional = (
                acceptance is not None and maximum > gtol and bool(acceptance.candidate(residual))
            )
            deferred = _DEFERRED_CERTIFICATES
            if (
                deferred is not None
                and deferrable
                and maximum <= gtol
                and CERTIFICATE_JACOBIAN_OUTSIDE_SEARCH_BUDGET
                and difference_mode == "forward"
                and not (levers is not None and levers.carried_across_advances)
            ):
                # Kernel-envelope DEFERRED_CERTIFICATE (inside its scope only): the previous
                # deferred certificate is taken now, after this solve's search; this iterate is
                # returned with its own certificate pending (the same statements, later)
                deferred.drain()
                index = deferred.next_index()
                if deferred.defers(index):
                    if levers is not None:
                        force_fresh = False
                        levers._source(stamp, "fresh_certificate")
                        levers._close(stamp, "converged")
                    result = OptimizeResult(
                        x=x,
                        fun=residual,
                        jac=None,
                        success=True,
                        active_mask=np.zeros(x.size),
                        nfev=work,
                        optimality=None,
                        message="domain-aware Newton converged",
                        certificate_jacobian_calls_outside_search_budget=None,
                        rejected_trial_evaluations=rejected,
                        **({} if line_search is None else {"line_search": line_search}),
                    )

                    def compute(result=result, x=x, residual=residual, iteration=iteration):
                        calls_before = certificate["calls"]
                        certificate["active"] = True
                        try:
                            if (
                                cg_stamp is not None
                                and column_groups.layout not in cgj._STATE["disabled"]
                            ):
                                jac = column_group_jacobian(x, residual, False)
                            else:
                                jac = jacobian(x, residual)
                        finally:
                            certificate["active"] = False
                        if pattern is not None and iteration == 0 and np.any(jac[~pattern] != 0.0):
                            raise ValueError(
                                "declared grouped Jacobian misses measured dependencies"
                            )
                        result["jac"] = jac
                        result["optimality"] = float(np.max(np.abs(jac.T @ residual)))
                        result["certificate_jacobian_calls_outside_search_budget"] = certificate[
                            "calls"
                        ]
                        result["rejected_trial_evaluations"] = rejected
                        if certificate["calls"] - calls_before != x.size:
                            # the caller charged one forward probe per unknown at the return
                            raise CertificateChargeMismatch(
                                f"certificate took {certificate['calls'] - calls_before} residual "
                                f"calls; {x.size} were charged at the return"
                            )
                        return jac

                    pending = PendingCertificate(index, compute, _hook_snapshot())
                    # the residual calls the former certificate makes (one forward probe per
                    # unknown, its retreats aside), charged by the caller at the return
                    pending.predicted_calls = int(x.size)
                    deferred.register(result, pending)
                    return result
            certificate["active"] = (
                maximum <= gtol or provisional
            ) and CERTIFICATE_JACOBIAN_OUTSIDE_SEARCH_BUDGET
            calls_before = certificate["calls"]
            # JACOBIAN_BROYDEN_UPDATES: the held Jacobian is the search Jacobian unless this iterate
            # is a certificate candidate or a rebuild is due (never outside a lever scope).
            reused = held is not None and not force_fresh and not (maximum <= gtol or provisional)
            try:
                if reused:
                    jac = held
                elif cg_stamp is not None and column_groups.layout not in cgj._STATE["disabled"]:
                    # Kernel-envelope JACOBIAN_COLUMN_GROUPS (armed only)
                    jac = column_group_jacobian(x, residual, not (maximum <= gtol or provisional))
                else:
                    jac = (
                        jacobian(x, residual)
                        if pattern is None or iteration == 0 or maximum <= gtol or provisional
                        else grouped_search_jacobian(x, residual)
                    )
            finally:
                certificate["active"] = False
            if levers is not None:
                force_fresh = False
                if reused:
                    levers._source(stamp, held_source)
                else:
                    held, held_source = jac, "fresh"
                    levers._source(
                        stamp,
                        "fresh_certificate" if maximum <= gtol or provisional else "fresh",
                    )
            if pattern is not None and iteration == 0 and np.any(jac[~pattern] != 0.0):
                raise ValueError("declared grouped Jacobian misses measured dependencies")
            accepted_by_rows = provisional and bool(acceptance.accepts(x, residual, jac))
            if provisional and not accepted_by_rows:
                # Not returned: this Jacobian is the next step's search Jacobian.
                charged = certificate["calls"] - calls_before
                certificate["calls"] = calls_before
                work += charged
                if work > max_nfev:
                    raise RuntimeError("domain-aware Newton residual-call budget exhausted")
            if maximum <= gtol or accepted_by_rows:
                if levers is not None:
                    levers._close(stamp, "converged")
                    if levers.carried_across_advances:
                        # JACOBIAN_CARRIED_ACROSS_DRY_ADVANCES: this certificate Jacobian starts
                        # the next solve of the same layout
                        levers.carried[layout] = np.array(jac, copy=True)
                return OptimizeResult(
                    x=x,
                    fun=residual,
                    jac=jac,
                    success=True,
                    active_mask=np.zeros(x.size),
                    nfev=work,
                    optimality=float(np.max(np.abs(jac.T @ residual))),
                    message=(
                        "domain-aware Newton converged"
                        if maximum <= gtol
                        else "domain-aware Newton converged under the declared per-row acceptance"
                    ),
                    certificate_jacobian_calls_outside_search_budget=certificate["calls"],
                    rejected_trial_evaluations=rejected,
                    **({} if line_search is None else {"line_search": line_search}),
                )
            column_scale = search_column_scale(jac)
            scaled_jac = jac * column_scale
            if acceptance is None:
                old_norm = np.linalg.norm(residual)
            else:
                merit_gates = merit_row_gates(acceptance, x, jac, residual.size)
                old_norm = gate_normalized_merit(residual, merit_gates, line_search["norm"])
            advanced = False
            for damping in (0.0, 1e-6, 1e-4, 0.01, 1.0, 100.0):
                matrix = np.vstack((scaled_jac, math.sqrt(damping) * np.eye(x.size)))
                rhs = np.r_[-residual, np.zeros(x.size)]
                step = np.linalg.lstsq(matrix, rhs, rcond=1e-14)[0] * column_scale
                # Tangent search directions near a bound. This changes a
                # search direction only; no stored coordinate is projected.
                numerical_distance = 32 * np.finfo(float).eps * np.maximum(1.0, np.abs(x))
                blocked = ((step > 0) & ((upper - x) < numerical_distance)) | (
                    (step < 0) & ((x - lower) < numerical_distance)
                )
                if np.any(blocked):
                    allowed = ~blocked
                    step[:] = 0.0
                    step[allowed] = (
                        np.linalg.lstsq(matrix[:, allowed], rhs, rcond=1e-14)[0]
                        * column_scale[allowed]
                    )
                alpha = 1.0
                for index, change in enumerate(step):
                    if change > 0:
                        alpha = min(alpha, 0.99 * (upper[index] - x[index]) / change)
                    elif change < 0:
                        alpha = min(alpha, 0.99 * (lower[index] - x[index]) / change)
                for retreat in range(20):
                    trial = x + alpha * step
                    try:
                        candidate = evaluate(trial)
                    except StepBudgetExceeded:
                        raise
                    except (ValueError, RuntimeError) as exc:
                        if "budget" in str(exc).lower() or "source" in str(exc).lower():
                            raise
                        alpha *= 0.5
                        continue
                    merit = (
                        np.linalg.norm(candidate)
                        if acceptance is None
                        else gate_normalized_merit(candidate, merit_gates, line_search["norm"])
                    )
                    if merit < (1 - 1e-4 * alpha) * old_norm:
                        if line_search is not None:
                            line_search["steps"].append(
                                {
                                    "damping": damping,
                                    "step_length": float(alpha),
                                    "retreats": retreat,
                                    "merit_before": float(old_norm),
                                    "merit_after": float(merit),
                                }
                            )
                        if levers is not None:
                            # JACOBIAN_BROYDEN_UPDATES: the rank-one update on the accepted step,
                            # then a rebuild at the next iterate after a stalled contraction
                            s = trial - x
                            denominator = float(s @ s)
                            if denominator > 0.0 and np.isfinite(denominator):
                                y = candidate - residual
                                held = held + np.outer(y - held @ s, s) / denominator
                                held_source = "broyden_updated"
                                levers._count(stamp, "broyden_updates")
                            if merit > BROYDEN_REBUILD_CONTRACTION_RATIO * old_norm:
                                force_fresh = True
                                levers._count(stamp, "rebuilds_after_stalled_contraction")
                        x, residual = trial, candidate
                        if np.max(np.abs(residual)) < best[0]:
                            best = (float(np.max(np.abs(residual))), x.copy())
                        advanced = True
                        break
                    alpha *= 0.5
                if advanced:
                    break
            if not advanced and reused:
                # JACOBIAN_BROYDEN_UPDATES: the held (updated or carried) direction found no
                # admissible decreasing step; the dense Jacobian is rebuilt at this same iterate
                # and the original step and failure test are taken on it
                force_fresh = True
                levers._count(stamp, "rebuilds_after_failed_direction")
                continue
            if not advanced:
                if levers is not None:
                    levers._close(stamp, "no admissible decreasing step")
                return OptimizeResult(
                    x=x,
                    fun=residual,
                    jac=jac,
                    success=False,
                    active_mask=np.zeros(x.size),
                    nfev=work,
                    optimality=float(np.max(np.abs(jac.T @ residual))),
                    message="domain-aware Newton has no admissible decreasing step",
                    rejected_trial_evaluations=rejected,
                    **({} if line_search is None else {"line_search": line_search}),
                )
    except Exception as exc:
        if levers is not None:
            levers._close(stamp, f"raised {type(exc).__name__}")
        exc.optimizer_debug = {
            "best_maximum_scaled_residual": best[0],
            "best_coordinates": best[1].tolist(),
            "residual_calls": work,
            "rejected_trials": rejected,
            "rejection_causes": rejection_causes,
        }
        if line_search is not None:
            exc.optimizer_debug["line_search"] = line_search
        raise
    if levers is not None:
        levers._close(stamp, "iteration limit")
    raise RuntimeError("domain-aware Newton iteration limit exhausted")


def physical_chart_solve(
    fun, x0, *, decode, encode, physical_bounds, temperature_indices, physical_fun=None, **kwargs
):
    """Pull back the unchanged native residual into physical primitives.

    The encode callback still enforces every original native solver bound.
    Temperature origins are translations, with unit slopes. The returned
    Jacobian is explicitly in physical coordinates, for a conditioning audit.
    """
    initial = np.asarray(decode(x0), dtype=float)
    origin = np.zeros(initial.size)
    origin[list(temperature_indices)] = initial[list(temperature_indices)]
    lower, upper = (np.asarray(v, dtype=float) for v in physical_bounds)
    native_lower, native_upper = (np.asarray(v) for v in kwargs.pop("bounds"))

    def native(values):
        encoded = np.asarray(encode(values), dtype=float)
        if np.any(encoded <= native_lower) or np.any(encoded >= native_upper):
            raise ValueError("physical trial left the original open native solver bounds")
        return encoded

    def residual(values):
        physical = values + origin
        encoded = native(physical)
        return fun(encoded) if physical_fun is None else physical_fun(physical, encoded)

    try:
        result = domain_aware_newton(
            residual, initial - origin, bounds=(lower - origin, upper - origin), **kwargs
        )
    except Exception as exc:
        debug = getattr(exc, "optimizer_debug", None)
        if debug is not None:
            physical = np.asarray(debug["best_coordinates"]) + origin
            debug["temperature_origins_k"] = origin.tolist()
            debug["best_physical_coordinates"] = physical.tolist()
            debug["best_native_coordinates"] = native(physical).tolist()
        raise
    result.physical_coordinates = result.x + origin
    result.x = native(result.physical_coordinates)
    if result.success:
        # Restore the native coordinate Jacobian for its original condition
        # gate. Physical increments remain the nonlinear search coordinates.
        # A raw matrix in K, kg/kg and seconds is not the native dimensionless
        # barrier matrix, particularly near a vanishing event duration.
        pullback = []
        for column in range(result.x.size):
            h = 1e-5 * max(1.0, abs(result.x[column]))
            for retreat in range(16):
                left = result.x.copy()
                right = result.x.copy()
                left[column] -= h
                right[column] += h
                try:
                    if np.any(left <= native_lower) or np.any(right >= native_upper):
                        raise ValueError("Jacobian pullback probe outside native bounds")
                    change = (np.asarray(decode(right)) - np.asarray(decode(left))) / (
                        right[column] - left[column]
                    )
                    pullback.append(change)
                    break
                except StepBudgetExceeded:
                    raise
                except ValueError, RuntimeError:
                    h *= 0.5
            else:
                raise RuntimeError("native Jacobian pullback has no admissible probes")
        pending = pending_certificate(result)
        if pending is None:
            result.native_jac = result.jac @ np.column_stack(pullback)
        else:
            # DEFERRED_CERTIFICATE: the pullback (the decode map only, no residual call) now, its
            # product with the certificate when the certificate is in
            pullback_matrix = np.column_stack(pullback)
            result.native_jac = None

            def native_product(jac, result=result, pullback_matrix=pullback_matrix):
                result.native_jac = jac @ pullback_matrix

            pending.after.append(native_product)
    result.physical_coordinate_audit = {
        "jacobian_coordinates": "physical primitives; translated temperatures",
        "temperature_origins_k": origin.tolist(),
        "native_solver_bounds_preserved": True,
        "physical_residual_changed": False,
        "native_jacobian_available": bool(result.success),
    }
    return result


def centered_temperature_solve(solver, fun, x0, *, temperature_index, **kwargs):
    """Solve an affine, unit-slope temperature increment; return physical x."""
    initial = np.array(x0, dtype=float, copy=True)
    origin = float(initial[temperature_index])
    initial[temperature_index] = 0.0
    lower, upper = (np.array(v, dtype=float, copy=True) for v in kwargs["bounds"])
    lower[temperature_index] -= origin
    upper[temperature_index] -= origin
    kwargs["bounds"] = (lower, upper)

    def residual(coordinates):
        physical = np.array(coordinates, copy=True)
        physical[temperature_index] += origin
        return fun(physical)

    result = solver(residual, initial, **kwargs)
    result.x = np.array(result.x, copy=True)
    result.x[temperature_index] += origin
    # The unit-slope affine map leaves every Jacobian column and its units
    # unchanged; the finite-difference probes are now around increments.
    result.temperature_coordinate_audit = {
        "origin_k": origin,
        "temperature_index": temperature_index,
        "affine_slope": 1.0,
        "physical_residual_changed": False,
    }
    return result


def audit_film_profile(
    film,
    *,
    interface_hexane_mass_fraction,
    bulk_hexane_mass_fraction,
    fluxes,
    closure_residual_budget_kg_m2_s=0.0,
):
    """Keep native auditing, with a bounded pure-bulk endpoint fallback.

    Only the native reconstructed-profile range failure can enter this
    branch. The logarithm signs, monotonic direction, flux-coordinate closure
    and reconstructed endpoint must all pass independent error budgets.
    """
    wi, wb = interface_hexane_mass_fraction, bulk_hexane_mass_fraction
    if not math.isfinite(closure_residual_budget_kg_m2_s) or closure_residual_budget_kg_m2_s < 0:
        raise ValueError("finite nonnegative prospective film closure error budget required")
    try:
        efr.conditioned_ideal_binary_mass_coordinate_residual(
            film, interface_hexane_mass_fraction=wi, bulk_hexane_mass_fraction=wb, fluxes=fluxes
        )
        return None
    except efr.FilmProfileAdmissibilityError as exc:
        if (
            str(exc) != "candidate implies a mass-fraction profile outside [0,1]"
            or wb not in (0.0, 1.0)
            or not 0.0 < wi < 1.0
        ):
            raise
        original = exc
    total = fluxes.total_kg_m2_s
    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    pe = total / conductance
    numerator = fluxes.hexane_kg_m2_s - wb * total
    denominator = fluxes.hexane_kg_m2_s - wi * total
    if total == 0.0 or numerator * denominator <= 0.0:
        raise original
    phi = efr._inverse_exprel(pe)
    diffusive = conductance * (wi - wb) * phi
    residual = math.fsum((fluxes.hexane_kg_m2_s, -wi * total, -diffusive))
    flux_scale = math.fsum((abs(fluxes.hexane_kg_m2_s), abs(wi * total), abs(diffusive)))
    flux_budget = 128.0 * np.finfo(float).eps * flux_scale + closure_residual_budget_kg_m2_s
    profile_term = denominator / conductance * efr._exprel(pe)
    implied = wi - profile_term
    endpoint_budget = 128.0 * np.finfo(float).eps * math.fsum((1.0, abs(wi), abs(profile_term)))
    endpoint_budget += closure_residual_budget_kg_m2_s / conductance * abs(efr._exprel(pe))
    # d w_h/dz has sign -denominator throughout the exact exponential
    # profile. Require the direction of the prescribed physical endpoints.
    monotone = denominator * (wi - wb) > 0.0
    if (
        not monotone
        or not math.isfinite(implied)
        or abs(residual) > flux_budget
        or abs(implied - wb) > endpoint_budget
    ):
        raise original
    return {
        "branch": "pure_bulk_endpoint_numerical_closure",
        "bulk_hexane_mass_fraction": wb,
        "unmodified_implied_bulk_hexane_mass_fraction": implied,
        "conditioned_residual_kg_m2_s": residual,
        "flux_error_budget_kg_m2_s": flux_budget,
        "prospective_closure_residual_budget_kg_m2_s": closure_residual_budget_kg_m2_s,
        "closure_error_amplification_to_bulk_fraction": abs(efr._exprel(pe)) / conductance,
        "endpoint_error_budget": endpoint_budget,
        "epsilon_multiplier": 128,
        "profile_monotonicity_checked": True,
        "composition_clamped": False,
    }
