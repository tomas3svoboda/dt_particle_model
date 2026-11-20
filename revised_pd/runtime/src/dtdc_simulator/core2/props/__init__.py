"""Frozen pure-fluid property authorities for GT-PS-2.

- ``hexane`` — PHY-049 published Span-Wagner (2003) 12-term Helmholtz EOS.
- ``water``  — PHY-030 IAPWS-95 (Wagner & Pruss 2002) reference Helmholtz EOS.

Both are self-contained (no CoolProp/REFPROP at runtime; those are oracles
only) and expose one common caloric datum per fluid whose numerical zero is
arbitrary and must cancel under datum shifts (PHY-021 datum-invariance rule).

Design contract shared by both modules
--------------------------------------
Each module builds its dimensionless Helmholtz energy ``a/(RT) = phi0 + phir``
in reduced coordinates ``delta = rho/rho_c``, ``tau = Tc/T`` and returns a
:class:`HelmholtzResult` carrying ``phir`` and its first/second derivatives.
The public phase-state basis is explicit: n-hexane uses molar
``PhaseState`` and mass-native IAPWS-95 uses ``MassPhaseState``; both expose
named conversions to the other basis.
All measurable properties are derived from that single potential, so pressure,
fugacity, saturation, density, caloric states and latent heat are mutually
consistent by construction (the inconsistency the frozen decisions call out in
the legacy Antoine + constant-cp treatment).
"""

from __future__ import annotations

from dtdc_simulator.core2.props.state import HelmholtzResult, MassPhaseState, PhaseState

__all__ = [
    "HelmholtzResult",
    "MassPhaseState",
    "PhaseState",
    "binary_gas",
    "hexane",
    "water",
]
