"""DCZ domain construction and tray-duty apportionment."""

from __future__ import annotations

import math

from dtdc_simulator.core.dt_types import DTSolverConstants, DTTray
from dtdc_simulator.core.zones import dcz


def _tray_q_Iv(tray: DTTray, A_bed_m2: float) -> float:
    return tray.Q_indirect_w / (A_bed_m2 * tray.bed_height_m) if tray.bed_height_m > 0.0 else 0.0


def _average_q_Iv_over_depth(trays: list[DTTray], depth_m: float, A_bed_m2: float) -> float:
    """Length/volume-average jacket heat density over a zone from its top.

    This lets the thin FTRZ cross a real-tray boundary without inheriting the
    first tray's duty density over its entire length.  The FTRZ free-boundary
    iteration calls it with its current length guess.
    """
    if depth_m <= 0.0:
        return _tray_q_Iv(trays[0], A_bed_m2)
    remaining_m = depth_m
    integral_w_m2 = 0.0
    covered_m = 0.0
    for tray in trays:
        share_m = min(max(remaining_m, 0.0), tray.bed_height_m)
        integral_w_m2 += _tray_q_Iv(tray, A_bed_m2) * share_m
        covered_m += share_m
        remaining_m -= share_m
        if remaining_m <= 1.0e-12:
            break
    if remaining_m > 1.0e-9:
        raise ValueError(
            f"FTRZ length ({depth_m:.4f} m) exceeds the remaining DT bed height ({covered_m:.4f} m)"
        )
    return integral_w_m2 / max(covered_m, 1.0e-12)


def _build_dcz_domain(
    remaining: list[DTTray],
    L_FTRZ_m: float,
    A_bed_m2: float,
    hQ: float,
    hM: float,
    aV: float,
    c: DTSolverConstants,
    nz_dcz: int,
    axial_cell_fractions: tuple[float, ...] | None = None,
) -> tuple[dcz.DCZConstants, tuple[float, ...]]:
    total_remaining_m = sum(tray.bed_height_m for tray in remaining)
    if L_FTRZ_m >= total_remaining_m:
        raise ValueError(
            f"FTRZ length ({L_FTRZ_m:.4f} m) leaves no DCZ domain within "
            f"the remaining DT bed ({total_remaining_m:.4f} m)"
        )

    segments: list[tuple[float, float, float]] = []
    ftrz_left_m = L_FTRZ_m
    z = 0.0
    for tray in remaining:
        ftrz_share_m = min(max(ftrz_left_m, 0.0), tray.bed_height_m)
        dcz_share_m = tray.bed_height_m - ftrz_share_m
        ftrz_left_m -= ftrz_share_m
        if dcz_share_m <= 1.0e-12:
            continue
        q_Iv_tray = _tray_q_Iv(tray, A_bed_m2)
        segments.append((z, z + dcz_share_m, q_Iv_tray))
        z += dcz_share_m
    L_DCZ_total_m = z

    fractions = (
        tuple(1.0 / nz_dcz for _ in range(nz_dcz))
        if axial_cell_fractions is None
        else axial_cell_fractions
    )
    if len(fractions) != nz_dcz or any(value <= 0.0 for value in fractions):
        raise ValueError("DCZ axial cell fractions must be positive and length nz_dcz")
    if not math.isclose(sum(fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("DCZ axial cell fractions must sum to one")
    cell_widths = tuple(L_DCZ_total_m * value for value in fractions)
    cell_faces = [0.0]
    for width in cell_widths:
        cell_faces.append(cell_faces[-1] + width)
    profile: list[float] = []
    for j in range(nz_dcz):
        cell_start = cell_faces[j]
        cell_end = cell_faces[j + 1]
        # Conservative cell average of the piecewise-constant tray duty.
        # Sampling at the cell centre made the total applied jacket heat
        # depend on whether a tray boundary happened to fall to the left or
        # right of a centre.  Exact overlap integration preserves
        # sum(q_j * A_bed * dz) for every axial mesh, including cells that
        # straddle two real trays.
        integrated_q_w_m2 = 0.0
        for start, end, seg_q in segments:
            overlap_m = max(0.0, min(cell_end, end) - max(cell_start, start))
            integrated_q_w_m2 += seg_q * overlap_m
        profile.append(integrated_q_w_m2 / cell_widths[j])

    dcz_c = dcz.DCZConstants(
        diameter_m=remaining[0].diameter_m,
        bed_height_m=L_DCZ_total_m,
        hM=hM,
        hQ=hQ,
        aV=aV,
        D_ax=c.D_ax,
        k_mixL=c.k_mixL,
        rho_V=c.rho_V,
        cp_V=c.cp_V,
        alpha_V=c.ftrz.bed_porosity,
        alpha_L=1.0 - c.ftrz.bed_porosity,
        particle=c.particle,
        dH_vap_water=c.ftrz.vapor_enthalpy_ref.dH_vap_water,
        antoine_water=c.ftrz.antoine_water,
        antoine_hexane=c.antoine_hexane,
        luikov=c.luikov,
        water_diffusivity=c.water_diffusivity,
        vapor_enthalpy_ref=c.ftrz.vapor_enthalpy_ref,
        pressure_pa=c.pressure_pa,
    )
    return dcz_c, tuple(profile)
