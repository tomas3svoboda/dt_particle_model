"""Shared result containers for the pure-fluid Helmholtz authorities.

These carry molar quantities in SI (J/mol, mol/m3, Pa, K). Per-mass views are
provided through the molar mass so callers can work in whichever basis their
balance uses; the frozen ledger closes in extensive SI either way.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HelmholtzResult:
    """Dimensionless residual Helmholtz energy and its reduced derivatives.

    ``phir`` = alpha^r(delta, tau); subscripts denote partial derivatives with
    respect to ``delta`` and ``tau`` (e.g. ``phir_dt`` = d2 alpha^r/ddelta dtau).
    ``phi0*`` mirror the ideal-gas part where a module supplies it.
    """

    delta: float
    tau: float
    phir: float
    phir_d: float
    phir_dd: float
    phir_t: float
    phir_tt: float
    phir_dt: float
    #: Third mixed derivative d3 alpha^r/(ddelta dtau^2); optional because only
    #: the CELL-02 native caloric derivatives consume it.
    phir_dtt: float = 0.0
    #: Fourth mixed derivative d4 alpha^r/(ddelta dtau^3); CELL-02a2.
    phir_dttt: float = 0.0
    # Ideal-gas part (optional; caloric absolutes need it, saturation does not).
    phi0: float = 0.0
    phi0_t: float = 0.0
    phi0_tt: float = 0.0

    @property
    def compressibility(self) -> float:
        """Z = p / (rho R T) = 1 + delta * alpha^r_delta."""
        return 1.0 + self.delta * self.phir_d


@dataclass(frozen=True)
class PhaseState:
    """A fully resolved single-phase point in molar SI units."""

    T: float  # K
    rho: float  # mol/m3
    p: float  # Pa
    Z: float  # -
    h: float  # J/mol   (ideal datum + residual; datum arbitrary)
    u: float  # J/mol
    s: float  # J/mol/K
    cv: float  # J/mol/K
    cp: float  # J/mol/K
    ln_fugacity: float  # ln(f/Pa)  — pure-fluid fugacity
    molar_mass: float  # kg/mol

    @property
    def rho_mass(self) -> float:
        """kg/m3."""
        return self.rho * self.molar_mass

    @property
    def h_mass(self) -> float:
        """J/kg."""
        return self.h / self.molar_mass

    @property
    def u_mass(self) -> float:
        """J/kg."""
        return self.u / self.molar_mass

    @property
    def fugacity(self) -> float:
        """Pa."""
        import math

        return math.exp(self.ln_fugacity)


@dataclass(frozen=True)
class MassPhaseState:
    """A fully resolved single-phase point in mass-specific SI units.

    IAPWS-95 is formulated with mass density and a specific gas constant.  A
    separate container keeps that native basis explicit instead of storing
    kg/m3 and J/kg in fields documented as molar.  The named molar views are
    available for mixture work, but there are intentionally no ambiguous
    ``rho``/``h``/``u`` aliases.
    """

    T: float  # K
    rho_mass: float  # kg/m3
    p: float  # Pa
    Z: float  # -
    h_mass: float  # J/kg (datum arbitrary)
    u_mass: float  # J/kg
    s_mass: float  # J/(kg K)
    cv_mass: float  # J/(kg K)
    cp_mass: float  # J/(kg K)
    ln_fugacity: float  # ln(f/Pa)
    molar_mass: float  # kg/mol

    @property
    def rho_molar(self) -> float:
        """mol/m3."""
        return self.rho_mass / self.molar_mass

    @property
    def h_molar(self) -> float:
        """J/mol."""
        return self.h_mass * self.molar_mass

    @property
    def u_molar(self) -> float:
        """J/mol."""
        return self.u_mass * self.molar_mass

    @property
    def s_molar(self) -> float:
        """J/(mol K)."""
        return self.s_mass * self.molar_mass

    @property
    def cv_molar(self) -> float:
        """J/(mol K)."""
        return self.cv_mass * self.molar_mass

    @property
    def cp_molar(self) -> float:
        """J/(mol K)."""
        return self.cp_mass * self.molar_mass

    @property
    def fugacity(self) -> float:
        """Pa."""
        import math

        return math.exp(self.ln_fugacity)
