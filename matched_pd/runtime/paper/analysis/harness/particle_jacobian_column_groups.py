"""Kernel-envelope JACOBIAN_COLUMN_GROUPS: finite-difference Jacobians by column groups.

Record: docs/GT_PS2_DTDC_KERNEL_JACOBIAN_COLUMN_GROUPS_2026-11-02.md (the lead's Class-A numerics of
2026-11-02, on the owner's direction "We need to find opportunities to speed up our algorithms").

The domain-aware Newton (particle_numerical_repairs) takes every finite-difference Jacobian one column
at a time: one complete residual evaluation per unknown. The measured structure of the particle
residual is local (record section 2): a radial chain of cells, a cut/front border and the surface
node. Where a residual row depends on at most one unknown of a group, perturbing the whole group at
once and reading each column from the rows it owns gives the same difference quotient, bit for bit:
the same step h = sqrt(eps) max(1, |x_j|), the same sign (the first of +h and -h inside the open
bounds, as the one-at-a-time loop takes it), the same subtraction and division, a row outside the
column's declared structure the same signed zero 0.0 / h.

The structure is DECLARED, never learned from one Jacobian (a measured pattern misses entries that
are numerically zero at the difference step; the record found two at the birth): the src's own
equation-structural supersets (cut_sparsity.build_cut_jacobian_sparsity for the front's cut layout,
cut_birth_integrator.birth_sparsity for the birth, both imported, the pinned files unedited) with the
surface node's three unknowns and three rows taken as dense border blocks, and for the dry stage the
cell stencil with the conservative surface elimination's dense rows (the outer cell's three rows and
the surface rows depend on every cell's T and y through the whole-particle inventory sums of
experiment_core2_dt_dry.advance's residual).

Declared values (a VALUED knob, off when unset, any other value refused):

* ``1``: every Jacobian of a declared site exact by column groups (search and certificate);
* ``dry-rowsum``: as ``1``, and in the dry stage's film-bearing and film-free solves the SEARCH
  Jacobian's dense rows assembled from the per-cell inventory-change terms the residual's sums add
  (chain 257's Broyden class: the iterate path is not bit-identical; the dense certificate at every
  returned or provisional iterate stays the exact Jacobian).

Counting (the former rate): a group probe of k columns is charged as the k residual calls of the
one-at-a-time loop, to the Newton's residual-call count (max_nfev and the certificate count, each
call checked in turn as the former loop checks it) and to the residual's own assembly counters
through the site's residual_charge (the birth's assembly budget, the dry advance's assembly
counter and half-step threshold, the front's work count), so every budget decision is the former
one; only the composition-path proof's points, consumed inside each actual residual call, are
consumed by fewer calls. The bit check's one-at-a-time residual calls are counted by neither (the
site's residual_uncounted context).

Guards (never a silent fill): a group probe that moves a row outside its columns' declared structure
discards the grouped Jacobian, takes the one-at-a-time Jacobian at that iterate and disables the
layout for the rest of the process (stamped); a group probe the physics refuses (a domain error)
gives its columns the one-at-a-time procedure (exact); in the first BIT_CHECK_SOLVES solves of a
process every exact grouped Jacobian is compared bit for bit with the one-at-a-time Jacobian at the
same iterate (its residual calls outside the solve's count) and a difference raises
ColumnGroupJacobianMismatch. The depletion-event solve, the hexane-exhausted vertex, the nested
three-unknown surface solves and every other site stay one-at-a-time.

physically_qualifying: false.
"""

from __future__ import annotations

import math

import numpy as np

RECORD = "docs/GT_PS2_DTDC_KERNEL_JACOBIAN_COLUMN_GROUPS_2026-11-02.md"
# The number of solves per process in which every exact grouped Jacobian is checked bit for bit
# against the one-at-a-time Jacobian at the same iterate (the brief's "first N solves", N declared).
# WAS absent until 2026-11-02 (no grouped Jacobian).
BIT_CHECK_SOLVES = 3
EXACT = "1"
DRY_ROWSUM = "dry-rowsum"
VALUES = (EXACT, DRY_ROWSUM)
# dry-rowsum: the least-squares fit of a dense row's two sum coefficients needs at least this many
# far groups (two unknowns); with fewer the search Jacobian is the exact grouped one (stamped).
ROWSUM_MINIMUM_FAR_GROUPS = 2
# The front's cut structure: cut_sparsity's two widenings, the face-limit basis (an O(1) seam
# widening) and the face-conditioned composition (a radius-two stencil). Both are needed: without
# the second a bit-checked positive-target solve (29 unknowns) moved three rows outside the
# structure (record section 4; the guard fell back to the one-at-a-time Jacobian), with both every
# checked Jacobian of the record's runs was bit-equal.
FRONT_FACE_LIMIT_REGULARIZED = True
FRONT_FACE_CONDITIONED_COMPOSITION = True

_STATE = {"checks_left": BIT_CHECK_SOLVES, "disabled": set(), "pending": [], "totals": None}


class ColumnGroupJacobianMismatch(RuntimeError):
    """An exact grouped Jacobian differed from the one-at-a-time Jacobian (a structure defect)."""

    def __init__(self, message):
        # "source" keeps the domain-aware Newton from reading this as a rejected trial
        super().__init__(f"column-group source structure defect: {message}")


class _StructureViolation(Exception):
    pass


class Counting:
    """A site's switch for its own residual-call counters: on, except inside uncounted() (the bit
    check's one-at-a-time calls)."""

    def __init__(self):
        self.on = True

    def uncounted(self):
        from contextlib import contextmanager

        @contextmanager
        def scope():
            previous, self.on = self.on, False
            try:
                yield
            finally:
                self.on = previous

        return scope()


def _fresh_totals():
    return {
        "solves": 0,
        "checked_solves": 0,
        "grouped_jacobians": 0,
        "rowsum_jacobians": 0,
        "dense_fallback_jacobians": 0,
        "group_evaluations": 0,
        "solo_column_evaluations": 0,
        "refused_group_probes": 0,
        "structure_violations": 0,
        "bit_checks": 0,
    }


_STATE["totals"] = _fresh_totals()


def reset():
    """Forget the process state (the bit-check count, disabled layouts, stamps)."""
    _STATE.update(checks_left=BIT_CHECK_SOLVES, disabled=set(), pending=[], totals=_fresh_totals())


def mode(environ=None):
    """The armed declared value ("1" or "dry-rowsum"), or None when the knob is unset."""
    import particle_kernel_envelope as envelope

    return envelope.valued(envelope.JACOBIAN_COLUMN_GROUPS, environ)


def drain():
    """Every solve stamp since the previous drain and the process totals; None when nothing ran."""
    solves, _STATE["pending"] = _STATE["pending"], []
    if not solves and _STATE["totals"]["solves"] == 0:
        return None
    return {"record": RECORD, "solves": solves, "totals": dict(_STATE["totals"])}


def greedy_colors(structure, extra_conflicts=()):
    """Deterministic greedy coloring of the column-intersection graph (largest degree first), as
    cut_sparsity._greedy_column_coloring orders it; extra_conflicts adds declared column pairs."""
    structure = np.asarray(structure, dtype=bool)
    n = structure.shape[1]
    neighbours = [set() for _ in range(n)]
    for row in structure:
        columns = np.flatnonzero(row).tolist()
        for c in columns:
            neighbours[c].update(columns)
    for a, b in extra_conflicts:
        neighbours[a].add(b)
        neighbours[b].add(a)
    for c in range(n):
        neighbours[c].discard(c)
    order = sorted(range(n), key=lambda c: (-len(neighbours[c]), c))
    colors = [-1] * n
    for c in order:
        used = {colors[k] for k in neighbours[c] if colors[k] >= 0}
        color = 0
        while color in used:
            color += 1
        colors[c] = color
    return colors


def _groups_from_colors(colors, last_singleton):
    n = len(colors)
    groups = {}
    for c in range(n):
        if last_singleton and c == n - 1:
            continue
        groups.setdefault(colors[c], []).append(c)
    ordered = [tuple(groups[k]) for k in sorted(groups)]
    if last_singleton:
        # the last probe is the one-at-a-time loop's last probe, x + h e_(n-1), bit for bit, so any
        # state a caller reads from the residual's last call is the former one
        ordered.append((n - 1,))
    return tuple(ordered)


class RowsumSpec:
    """dry-rowsum: which dense rows are assembled, which columns are far, and the per-cell terms.

    terms(values) returns (water_per_cell, hexane_per_cell) as float arrays from the search
    coordinates, the terms the residual's inventory sums add (each local to its own cell)."""

    def __init__(self, dense_rows, far_columns, cell_of_column, terms):
        self.dense_rows = tuple(int(r) for r in dense_rows)
        self.far_columns = tuple(int(c) for c in far_columns)
        self.cell_of_column = dict(cell_of_column)
        self.terms = terms


class ColumnGroupPlan:
    """One solve's declared structure (rows x unknowns, a superset of every derivative), its groups
    and, for dry-rowsum, the search structure and its groups."""

    def __init__(
        self, structure, site, *, rowsum=None, residual_charge=None, residual_uncounted=None
    ):
        structure = np.asarray(structure, dtype=bool)
        # the site's own assembly counters: residual_charge(k) books k residual calls as the
        # residual itself would (its checks in order, each raising as the residual raises);
        # residual_uncounted() a context in which the residual books nothing (the bit check)
        self.residual_charge = residual_charge
        self.residual_uncounted = residual_uncounted
        if structure.ndim != 2 or not structure.any(axis=0).all():
            raise ValueError("a declared column-group structure needs one supported row per column")
        self.structure = structure
        self.site = site
        self.layout = (site, structure.shape, structure.tobytes())
        self.groups = _groups_from_colors(greedy_colors(structure), True)
        self.rowsum = rowsum
        self.search_structure = None
        self.search_groups = None
        if rowsum is not None:
            search = structure.copy()
            dense = list(rowsum.dense_rows)
            far = list(rowsum.far_columns)
            search[np.ix_(dense, far)] = False
            # near columns (touching a dense row exactly) are never grouped with far columns, so a
            # dense row's change in a group holding far columns comes from the sums alone
            near = [c for c in range(structure.shape[1]) if search[dense, c].any()]
            self.search_structure = search
            self.search_groups = _groups_from_colors(
                greedy_colors(search, [(a, b) for a in near for b in far]), True
            )

    def stamp_static(self):
        out = {
            "site": self.site,
            "unknowns": int(self.structure.shape[1]),
            "rows": int(self.structure.shape[0]),
            "declared_nonzeros": int(self.structure.sum()),
            "groups": len(self.groups),
        }
        if self.rowsum is not None:
            out["search_groups"] = len(self.search_groups)
            out["assembled_rows"] = list(self.rowsum.dense_rows)
            out["far_columns"] = list(self.rowsum.far_columns)
        return out


def open_solve(plan):
    """The solve's stamp and whether its exact grouped Jacobians are bit-checked; None when the
    layout was disabled by a violation earlier in this process."""
    if plan.layout in _STATE["disabled"]:
        return None
    checked = _STATE["checks_left"] > 0
    if checked:
        _STATE["checks_left"] -= 1
        _STATE["totals"]["checked_solves"] += 1
    _STATE["totals"]["solves"] += 1
    stamp = {
        **plan.stamp_static(),
        "bit_checked": checked,
        "grouped_jacobians": 0,
        "rowsum_jacobians": 0,
        "dense_fallback_jacobians": 0,
        "group_evaluations": 0,
        "solo_column_evaluations": 0,
        "refused_group_probes": 0,
        "structure_violations": [],
        "bit_checks": 0,
        "bit_checks_equal": 0,
        "rowsum": [],
    }
    _STATE["pending"].append(stamp)
    return stamp


def _count(stamp, key, k=1):
    stamp[key] += k
    _STATE["totals"][key] += k


def _first_step(values, lower, upper, column):
    """The one-at-a-time loop's first admissible trial value of a forward-difference column (retreat
    0), or None when the loop would need a retreat or a probe the bounds refuse at h itself."""
    h = math.sqrt(np.finfo(float).eps) * max(1.0, abs(values[column]))
    for sign in (1.0, -1.0):
        trial = values.copy()
        trial[column] += sign * h
        if not lower[column] < trial[column] < upper[column]:
            continue
        represented = trial[column] - values[column]
        if represented == 0:
            return None
        return trial[column], represented
    return None


def grouped_jacobian(
    plan, stamp, values, residual, lower, upper, evaluate, dense_column, *, search, charge
):
    """The Jacobian at values by column groups (see the module docstring); raises
    _StructureViolation when a group probe moves a row outside its columns' declared structure."""
    rowsum = plan.rowsum if search else None
    structure = plan.search_structure if rowsum is not None else plan.structure
    groups = plan.search_groups if rowsum is not None else plan.groups
    m, n = plan.structure.shape
    if residual.shape != (m,):
        raise ValueError("declared column-group structure does not match the residual rows")
    first = [_first_step(values, lower, upper, j) for j in range(n)]
    matrix = np.empty((m, n))
    solo = []
    samples = {}
    for group in groups:
        columns = [j for j in group if first[j] is not None]
        solo.extend(j for j in group if first[j] is None)
        if not columns:
            continue
        trial = values.copy()
        for j in columns:
            trial[j] = first[j][0]
        try:
            sample = evaluate(trial)
        except (ValueError, RuntimeError) as exc:
            text = str(exc).lower()
            if "budget" in text or "source" in text:
                raise
            _count(stamp, "refused_group_probes")
            solo.extend(columns)
            continue
        _count(stamp, "group_evaluations")
        # the former rate: the other columns' residual calls booked as the former loop books them
        charge(len(columns) - 1)
        allowed = structure[:, columns].any(axis=1)
        far_in_group = rowsum is not None and any(j in rowsum.far_columns for j in columns)
        if far_in_group:
            allowed[list(rowsum.dense_rows)] = True
        moved = sample != residual
        if np.any(moved & ~allowed):
            rows = np.flatnonzero(moved & ~allowed).tolist()
            raise _StructureViolation({"rows": rows[:12], "columns": list(columns)})
        for j in columns:
            difference = sample - residual
            difference[~structure[:, j]] = 0.0
            matrix[:, j] = difference / first[j][1]
        samples[tuple(columns)] = (trial, sample)
    if rowsum is not None:
        # far columns the fit cannot assemble take the one-at-a-time procedure (exact)
        solo.extend(_assemble_rowsum(plan, stamp, values, residual, matrix, samples, first, solo))
    for j in sorted(set(solo)):
        _count(stamp, "solo_column_evaluations")
        matrix[:, j] = dense_column(values, residual, j)
    return matrix


def _assemble_rowsum(plan, stamp, values, residual, matrix, samples, first, solo):
    rowsum = plan.rowsum
    base_w, base_h = (np.asarray(v, dtype=float) for v in rowsum.terms(values))
    far = set(rowsum.far_columns)
    sums, changes, per_column = [], [], {}
    for columns, (trial, sample) in samples.items():
        held = [j for j in columns if j in far]
        if not held:
            continue
        w, h = (np.asarray(v, dtype=float) for v in rowsum.terms(trial))
        dw, dh = w - base_w, h - base_h
        sums.append((math.fsum(dw), math.fsum(dh)))
        changes.append(sample[list(rowsum.dense_rows)] - residual[list(rowsum.dense_rows)])
        for j in held:
            cell = rowsum.cell_of_column[j]
            per_column[j] = (dw[cell] / first[j][1], dh[cell] / first[j][1])
    record = {"far_groups": len(sums), "assembled": False}
    if len(sums) < ROWSUM_MINIMUM_FAR_GROUPS:
        record["reason"] = "too few far groups; the far columns one at a time"
        stamp["rowsum"].append(record)
        return [j for j in per_column if j not in solo]
    a = np.asarray(sums)
    b = np.asarray(changes)
    coefficients, _, rank, _ = np.linalg.lstsq(a, b, rcond=None)
    fitted = a @ coefficients
    scale = np.maximum(np.max(np.abs(b), axis=0), np.finfo(float).tiny)
    record.update(
        assembled=True,
        rank=int(rank),
        fit_relative_residual=float(np.max(np.abs(fitted - b) / scale)),
    )
    for j, (gw, gh) in per_column.items():
        matrix[list(rowsum.dense_rows), j] = coefficients[0] * gw + coefficients[1] * gh
    stamp["rowsum"].append(record)
    _count(stamp, "rowsum_jacobians")
    return []


def violation(plan, stamp, detail):
    """A group probe moved a row outside the declared structure: the layout is disabled for the rest
    of this process and the violation stamped (the caller then takes the one-at-a-time Jacobian)."""
    _STATE["disabled"].add(plan.layout)
    stamp["structure_violations"].append(detail)
    stamp["structure_violation_count"] = len(stamp["structure_violations"])
    _STATE["totals"]["structure_violations"] += 1


# ---------------------------------------------------------------- the declared structures


def _spans(blocks):
    out, start = {}, 0
    for name, size in blocks:
        out[name] = list(range(start, start + size))
        start += size
    return out


def surface_extension(native, unknown_blocks, residual_blocks):
    """The native structure with the surface node's three unknowns and three rows appended.

    The surface node (temperature, film coordinate, liquid or event duration) closes the outer
    face of the outermost dry piece. Its three rows depend on the outer two dry pieces' T and y,
    the outer two Stefan fluxes, the two border unknowns (front coordinate, interface temperature)
    and its own three unknowns; its three unknowns reach the outer two dry pieces' three rows, the
    three Rankine-Hugoniot rows (the face-limit basis can carry a single dry piece's rows there)
    and its own three rows. Every extra entry is a declared over-approximation of the measured
    surface coupling (record section 2: the surface rows measured on the outer dry piece, its
    outer Stefan flux, the front coordinate and the surface unknowns)."""
    unknown = _spans(unknown_blocks)
    residual = _spans(residual_blocks)
    m, n = native.shape
    structure = np.zeros((m + 3, n + 3), dtype=bool)
    structure[:m, :n] = native
    surface_columns = [n, n + 1, n + 2]
    surface_rows = [m, m + 1, m + 2]
    reached = []
    for name in ("dry_temperature", "dry_y_hexane", "dry_total_stefan_flux"):
        reached += unknown[name][-2:]
    reached += unknown["front_z"] + unknown["interface_temperature"] + surface_columns
    structure[np.ix_(surface_rows, reached)] = True
    rows = []
    for name in ("dry_water", "dry_hexane", "dry_energy"):
        rows += residual[name][-2:]
    rows += residual["rh_water"] + residual["rh_hexane"] + residual["rh_energy"] + surface_rows
    structure[np.ix_(rows, surface_columns)] = True
    return structure


def front_plan(
    layout,
    *,
    stall_slot=None,
    stall_row=None,
    duration_column=False,
    site="front",
    residual_charge=None,
    residual_uncounted=None,
):
    """The front's cut layout (core2_water_film_node's physical chart: the native unknowns then the
    surface node's three, the native rows then its three): cut_sparsity's structural superset with
    both widenings (FRONT_FACE_LIMIT_REGULARIZED, FRONT_FACE_CONDITIONED_COMPOSITION) and the surface
    node by
    surface_extension (its third unknown a dense column when it is the film-depletion event's
    duration); the stalled front's slot and row deleted as the solve does."""
    from dtdc_simulator.core2.particle import cut_sparsity

    native = cut_sparsity.build_cut_jacobian_sparsity(
        layout,
        face_conditioned_composition=FRONT_FACE_CONDITIONED_COMPOSITION,
        face_limit_regularized=FRONT_FACE_LIMIT_REGULARIZED,
    ).jac_sparsity.toarray()
    structure = surface_extension(native, layout.unknown_blocks, layout.residual_blocks)
    if duration_column:
        # the film-depletion event: the surface node's third unknown is the event duration, which
        # every row's backward-Euler step reads (a dense column; record section 4: a bit-checked
        # stalled event moved every wet row)
        structure[:, -1] = True
    if stall_slot is not None:
        structure = np.delete(structure, stall_slot, axis=1)
    if stall_row is not None:
        structure = np.delete(structure, stall_row, axis=0)
    return ColumnGroupPlan(
        structure, site, residual_charge=residual_charge, residual_uncounted=residual_uncounted
    )


def birth_plan(layout, *, site="birth", residual_charge=None, residual_uncounted=None):
    """The birth layout: cut_birth_integrator.birth_sparsity with the surface node by
    surface_extension."""
    from dtdc_simulator.core2.particle import cut_birth_integrator as bi

    native = bi.birth_sparsity(layout).matrix.toarray().astype(bool)
    structure = surface_extension(native, layout.unknown_blocks, layout.residual_blocks)
    return ColumnGroupPlan(
        structure, site, residual_charge=residual_charge, residual_uncounted=residual_uncounted
    )


def dry_structure(cells, unknowns, rows):
    """experiment_core2_dt_dry.advance's residual: unknowns T (cells), q (cells), the Stefan fluxes
    (cells; flux k on the face between cells k and k + 1, the last on the surface face) and the
    surface unknowns; rows water, hexane, energy (cells) and the surface rows. A cell row depends on
    the T and q of its neighbours and on its two adjacent Stefan fluxes; the outer cell's three rows
    and every surface row depend on every cell's T and q (the conservative elimination's inventory
    sums), on the two outer Stefan fluxes and on every surface unknown; every surface unknown is a
    dense column."""
    n = cells
    structure = np.zeros((rows, unknowns), dtype=bool)
    for block in range(3):
        for c in range(n):
            row = block * n + c
            for k in range(max(0, c - 1), min(n, c + 2)):
                structure[row, k] = True
                structure[row, n + k] = True
            for k in (c - 1, c):
                if 0 <= k < n:
                    structure[row, 2 * n + k] = True
    dense_rows = [n - 1, 2 * n - 1, 3 * n - 1, *range(3 * n, rows)]
    for row in dense_rows:
        structure[row, : 2 * n] = True
        structure[row, 2 * n + n - 2 : 3 * n] = True
    structure[:, 3 * n :] = True
    return structure, dense_rows


def dry_plan(
    cells,
    unknowns,
    rows,
    *,
    rowsum_terms=None,
    site="dry",
    residual_charge=None,
    residual_uncounted=None,
):
    structure, dense_rows = dry_structure(cells, unknowns, rows)
    rowsum = None
    if rowsum_terms is not None:
        n = cells
        far_cells = range(0, n - 2)
        far = [c for c in far_cells] + [n + c for c in far_cells]
        cell_of = {c: c for c in far_cells}
        cell_of.update({n + c: c for c in far_cells})
        rowsum = RowsumSpec(dense_rows, far, cell_of, rowsum_terms)
    return ColumnGroupPlan(
        structure,
        site,
        rowsum=rowsum,
        residual_charge=residual_charge,
        residual_uncounted=residual_uncounted,
    )
