"""``core2`` — clean re-implementation of the DTDC core against the frozen
GT-PS-2 physics ledger (``release/physics_decisions.yaml``, 53/53 FROZEN).

This package is intentionally separate from the legacy ``core/`` so both stay
runnable for A/B comparison. The legacy transient DAE and the ``model.py``
lag-toward-steady surrogate are being retired; nothing here imports them.

Build order (each phase gated by owner review — see
``docs/CORE2_ARCHITECTURE.md``):

    props/      PHY-049 Span-Wagner n-hexane, PHY-030 IAPWS-95 water   (Phase 0)
    particle    resolved wet core + Fickian dry shell + RH front       (Phase 1)
    zones       PHZ / FTRZ / DCZ floating regimes, wall node, pressure (Phase 2)
    assembly    whole-DTDC backward-Euler DAE, RTD, discharge, DC       (Phase 3)
    solver      semismooth-Newton active-set + analytic sparse Jacobian (Phase 4)
    shim        RuntimeFacade-compatible Model.step adapter             (Phase 5)

Only Phase 0 (the two property authorities) plus a single quasi-steady tray
ledger (``core2.tray``) are implemented in the first review slice.
"""

from __future__ import annotations

__all__ = ["props", "tray"]
