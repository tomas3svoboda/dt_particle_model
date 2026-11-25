r"""Equation-structural Jacobian sparsity for the Gate-1g cut residual.

This module describes which residual equations *may* depend on which
cut-transport unknowns.  It does not evaluate derivatives and does not claim
an analytic Jacobian or solver-performance qualification.  The returned
SciPy CSR matrix is suitable for ``scipy.optimize.least_squares``
``jac_sparsity`` because fixed residual scales do not alter structural zeros.
The optional face-conditioned composition chart is triangular rather than
diagonal: physical ``y_i`` then also depends on the temperatures that define
its two adjacent symmetric face states.  Its exact radius-one coordinate
dependence is expanded explicitly below.

Regional finite-volume equations have nearest-neighbour state stencils and
two adjacent face-flux unknowns.  The three Rankine--Hugoniot rows touch only
the two interface traces and the innermost dry Stefan flux.  Front coordinate
``z`` and interface temperature ``T_Gamma`` are conservatively treated as
dense border columns.  That deliberate over-approximation prevents false
sparsity while retaining O(N) nonzeros and a bounded column coloring.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse

from dtdc_simulator.core2.particle import cut_transport as cut


@dataclass(frozen=True)
class StructuralBlock:
    """Named half-open span in the documented residual or unknown ordering."""

    name: str
    start: int
    stop: int

    def __post_init__(self) -> None:
        if not self.name or self.start < 0 or self.stop <= self.start:
            raise ValueError("a structural block needs a name and positive span")

    @property
    def size(self) -> int:
        return self.stop - self.start

    @property
    def slice(self) -> slice:
        return slice(self.start, self.stop)

    def at(self, local_index: int) -> int:
        if not 0 <= local_index < self.size:
            raise IndexError(f"{self.name} local index is outside its block")
        return self.start + local_index


@dataclass(frozen=True)
class CutJacobianSparsity:
    """Immutable row supports and a deterministic valid column coloring."""

    layout: cut.CutTransportLayout
    unknown_blocks: tuple[StructuralBlock, ...]
    residual_blocks: tuple[StructuralBlock, ...]
    row_columns: tuple[tuple[int, ...], ...]
    column_colors: tuple[int, ...]
    border_column_names: tuple[str, ...] = (
        "front_z",
        "interface_temperature",
    )
    border_row_names: tuple[str, ...] = (
        "rh_water",
        "rh_hexane",
        "rh_energy",
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if len(self.row_columns) != self.layout.residual_count:
            raise ValueError("one structural support is required per residual row")
        if len(self.column_colors) != self.layout.unknown_count:
            raise ValueError("one color is required per unknown column")
        if any(not columns for columns in self.row_columns):
            raise ValueError("a residual row cannot have empty structural support")
        if any(
            column < 0 or column >= self.layout.unknown_count
            for columns in self.row_columns
            for column in columns
        ):
            raise ValueError("structural support contains an invalid column")
        if any(color < 0 for color in self.column_colors):
            raise ValueError("column colors must be non-negative")
        for columns in self.row_columns:
            colors = tuple(self.column_colors[column] for column in columns)
            if len(colors) != len(set(colors)):
                raise ValueError("columns sharing a residual row must have distinct colors")

    @property
    def shape(self) -> tuple[int, int]:
        return self.layout.residual_count, self.layout.unknown_count

    @property
    def nnz(self) -> int:
        return sum(len(columns) for columns in self.row_columns)

    @property
    def color_count(self) -> int:
        return max(self.column_colors, default=-1) + 1

    @property
    def column_groups(self) -> tuple[tuple[int, ...], ...]:
        return tuple(
            tuple(
                column
                for column, column_color in enumerate(self.column_colors)
                if column_color == color
            )
            for color in range(self.color_count)
        )

    @property
    def jac_sparsity(self) -> sparse.csr_matrix:
        """Return a fresh boolean CSR matrix accepted by SciPy least-squares."""

        indptr = np.empty(self.layout.residual_count + 1, dtype=np.int64)
        indptr[0] = 0
        for row, columns in enumerate(self.row_columns, start=1):
            indptr[row] = indptr[row - 1] + len(columns)
        indices = np.fromiter(
            (column for columns in self.row_columns for column in columns),
            dtype=np.int64,
            count=self.nnz,
        )
        data = np.ones(self.nnz, dtype=bool)
        return sparse.csr_matrix((data, indices, indptr), shape=self.shape)

    def supports(self, residual_row: int, unknown_column: int) -> bool:
        """Return whether a derivative is admitted by the structural pattern."""

        if not 0 <= residual_row < self.layout.residual_count:
            raise IndexError("residual row is outside the pattern")
        if not 0 <= unknown_column < self.layout.unknown_count:
            raise IndexError("unknown column is outside the pattern")
        return unknown_column in self.row_columns[residual_row]


def build_cut_jacobian_sparsity(
    layout: cut.CutTransportLayout,
    *,
    face_conditioned_composition: bool = False,
    face_limit_regularized: bool = False,
) -> CutJacobianSparsity:
    """Build the conservative O(N) sparsity pattern for one fixed cut chart."""

    if not isinstance(layout, cut.CutTransportLayout):
        raise TypeError("cut Jacobian sparsity requires a CutTransportLayout")
    if not layout.is_square:
        raise ValueError("cut Jacobian sparsity requires the square cut residual")
    if not isinstance(face_conditioned_composition, bool):
        raise TypeError("face-conditioned-composition flag must be boolean")
    if not isinstance(face_limit_regularized, bool):
        raise TypeError("face-limit-regularized flag must be boolean")

    unknown_blocks = _block_spans(layout.unknown_blocks)
    residual_blocks = _block_spans(layout.residual_blocks)
    unknown = {block.name: block for block in unknown_blocks}
    residual = {block.name: block for block in residual_blocks}
    supports: list[set[int]] = [set() for _ in range(layout.residual_count)]

    wet_equations = (residual["wet_water"], residual["wet_energy"])
    for equation in wet_equations:
        for piece in range(layout.wet_piece_count):
            row = equation.at(piece)
            for neighbour in range(
                max(0, piece - 1),
                min(layout.wet_piece_count, piece + 2),
            ):
                supports[row].add(unknown["wet_temperature"].at(neighbour))
                supports[row].add(unknown["wet_retained_water"].at(neighbour))

    dry_equations = (
        residual["dry_water"],
        residual["dry_hexane"],
        residual["dry_energy"],
    )
    for equation in dry_equations:
        for piece in range(layout.dry_piece_count):
            row = equation.at(piece)
            for neighbour in range(
                max(0, piece - 1),
                min(layout.dry_piece_count, piece + 2),
            ):
                supports[row].add(unknown["dry_temperature"].at(neighbour))
                supports[row].add(unknown["dry_y_hexane"].at(neighbour))
            supports[row].add(unknown["dry_total_stefan_flux"].at(piece))
            supports[row].add(unknown["dry_total_stefan_flux"].at(piece + 1))

    interface_columns = (
        unknown["wet_temperature"].at(layout.wet_piece_count - 1),
        unknown["wet_retained_water"].at(layout.wet_piece_count - 1),
        unknown["dry_temperature"].at(0),
        unknown["dry_y_hexane"].at(0),
        unknown["dry_total_stefan_flux"].at(0),
    )
    for name in ("rh_water", "rh_hexane", "rh_energy"):
        supports[residual[name].start].update(interface_columns)

    if face_limit_regularized:
        # The sparse production integrator replaces the first dry rows by the
        # exact conservative face-limit basis
        #
        #   C_w = B_w,cut + B_d,0 + RH_w,
        #   C_h = B_d,0 + RH_h,
        #   C_E = B_E,w,cut + (B_E,d,0 - H_REF B_h,d,0)
        #         + (RH_E - H_REF RH_h).
        #
        # Admit the union of every source-row support before coloring.  This
        # is an O(1) seam widening and therefore preserves the O(N) radial
        # structure.  Original-basis consumers keep their independently
        # verified narrower pattern.
        wet_water_cut = residual["wet_water"].at(layout.wet_piece_count - 1)
        wet_energy_cut = residual["wet_energy"].at(layout.wet_piece_count - 1)
        dry_water_first = residual["dry_water"].at(0)
        dry_hexane_first = residual["dry_hexane"].at(0)
        dry_energy_first = residual["dry_energy"].at(0)
        rh_water = residual["rh_water"].start
        rh_hexane = residual["rh_hexane"].start
        rh_energy = residual["rh_energy"].start
        supports[dry_water_first].update(supports[wet_water_cut])
        supports[dry_water_first].update(supports[rh_water])
        original_dry_hexane_support = set(supports[dry_hexane_first])
        supports[dry_hexane_first].update(supports[rh_hexane])
        supports[dry_energy_first].update(supports[wet_energy_cut])
        supports[dry_energy_first].update(original_dry_hexane_support)
        supports[dry_energy_first].update(supports[rh_energy])
        supports[dry_energy_first].update(supports[rh_hexane])

    if face_conditioned_composition:
        # The exact nonlinear chart has the triangular coordinate relation
        #
        #   y_i = y_i(eta_i, T_{i-1}, T_i, T_{i+1}, T_Gamma),
        #
        # with the absent outer neighbour replaced by the fixed boundary
        # temperature.  Expand every physical y-column dependency through
        # that relation before coloring.  This grows a regional radius-one
        # physical stencil to at most radius two and therefore remains O(N).
        dry_y = unknown["dry_y_hexane"]
        dry_temperature = unknown["dry_temperature"]
        for row_support in supports:
            composition_pieces = tuple(
                column - dry_y.start for column in row_support if dry_y.start <= column < dry_y.stop
            )
            for piece in composition_pieces:
                for neighbour in range(
                    max(0, piece - 1),
                    min(layout.dry_piece_count, piece + 2),
                ):
                    row_support.add(dry_temperature.at(neighbour))

    # Geometry and interface-state effects are deliberately dense border
    # columns.  This includes harmless structural overestimates, but it cannot
    # omit a derivative as the cut approaches either side of its master cell.
    border_columns = (
        unknown["front_z"].start,
        unknown["interface_temperature"].start,
    )
    for row_support in supports:
        row_support.update(border_columns)

    row_columns = tuple(tuple(sorted(columns)) for columns in supports)
    colors = _greedy_column_coloring(layout.unknown_count, row_columns)
    return CutJacobianSparsity(
        layout=layout,
        unknown_blocks=unknown_blocks,
        residual_blocks=residual_blocks,
        row_columns=row_columns,
        column_colors=colors,
    )


def _block_spans(blocks: tuple[tuple[str, int], ...]) -> tuple[StructuralBlock, ...]:
    result: list[StructuralBlock] = []
    start = 0
    for name, size in blocks:
        result.append(StructuralBlock(name, start, start + size))
        start += size
    return tuple(result)


def _greedy_column_coloring(
    column_count: int,
    row_columns: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    """Color the column-intersection graph in deterministic degree order."""

    neighbours = [set() for _ in range(column_count)]
    for columns in row_columns:
        for position, column in enumerate(columns):
            neighbours[column].update(columns[:position])
            neighbours[column].update(columns[position + 1 :])
    order = sorted(range(column_count), key=lambda column: (-len(neighbours[column]), column))
    colors = [-1] * column_count
    for column in order:
        forbidden = {
            colors[neighbour] for neighbour in neighbours[column] if colors[neighbour] >= 0
        }
        color = 0
        while color in forbidden:
            color += 1
        colors[column] = color
    return tuple(colors)


__all__ = [
    "CutJacobianSparsity",
    "StructuralBlock",
    "build_cut_jacobian_sparsity",
]
