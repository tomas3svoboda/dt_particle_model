"""The declared roster of study-lane surrogates.

One module owns the roster so the build script and the run-time facade can
never disagree about what is surrogated, on which domain, or against which
pinned entry point.  Everything numeric in a pack — coefficients, degrees,
measured errors — is produced by the builder; everything structural — names,
domains, patch targets, output fields — is declared HERE and copied into the
pack for the record.

Two construction shapes appear, and the split is a deliberate engineering
judgement, recorded in
``docs/GT_PS2_S2_STUDY_SURROGATE_LANE_RECORD_2026-09-02.md``:

* **One-dimensional solved outputs are fitted directly.**  ``saturation(T)``
  and ``saturation_temperature(p)`` return few numbers, are analytic in their
  single argument away from the critical point, and are the calls whose
  iteration counts dominate (a 35-step bisection, a 2x2 Newton seeded by two
  density solves).  Fitting the outputs deletes the whole iteration.

* **The two-dimensional solve is replaced at the ROOT, not at the outputs.**
  ``density(T, p, phase)`` is one scalar Newton root; ``state_Tp`` then derives
  ten mutually consistent properties from that root through the PINNED
  Helmholtz derivative table.  Surrogating the root and leaving the algebra
  pinned keeps every Maxwell relation between the ten fields exact by
  construction, needs one series instead of ten, and still deletes the
  ten-to-twenty residual evaluations the Newton walk costs.  The residual
  evaluation that remains is the reason this is a ~5x call, not a ~60x one —
  an honest trade of speed for structural consistency.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "SURROGATE_SPECS",
    "VAPOUR_BRANCH_EXCLUSION",
    "Axis",
    "SurrogateSpec",
    "spec_by_name",
]


@dataclass(frozen=True)
class Axis:
    """One declared domain axis of a surrogate."""

    variable: str
    unit: str
    lo: float
    hi: float

    def contains(self, value: float) -> bool:
        return self.lo <= value <= self.hi

    def as_dict(self) -> dict[str, object]:
        return {"variable": self.variable, "unit": self.unit, "lo": self.lo, "hi": self.hi}


@dataclass(frozen=True)
class SurrogateSpec:
    """A declared surrogate: what is replaced, where, and with what outputs."""

    name: str
    fluid: str
    kind: str
    #: Attribute(s) of the pinned module replaced when this surrogate installs.
    #: Several targets share one spec where they read the same solved object
    #: (the four virial derivatives all come from one residual evaluation).
    patch_targets: tuple[str, ...]
    axes: tuple[Axis, ...]
    outputs: tuple[str, ...]
    #: ``liquid``/``vapor`` for the density specs; empty otherwise.
    phase: str = ""
    #: What is actually fitted.  ``identity`` fits the solved output itself.
    #: ``compressibility`` fits Z = p / (rho R T) and reconstructs
    #: rho = p / (Z R T) -- the right choice on a vapour branch, where rho
    #: spans nearly three decades across the declared box while Z stays within
    #: a couple of percent of unity.  A polynomial that is near-minimax in
    #: ABSOLUTE error on rho is useless in RELATIVE error at the thin end of
    #: that range; on Z the two measures coincide.
    transform: str = "identity"
    degree_ladder: tuple[int, ...] = (10, 12, 14, 16, 20, 24)
    #: Second-axis ladder for two-dimensional specs.
    second_degree_ladder: tuple[int, ...] = (6, 8, 10, 12, 14)
    #: Relative-error target the builder aims for on the ladder.
    relative_target: float = 1.0e-10
    #: Acceptance threshold; a measured maximum relative error above this makes
    #: the surrogate NOT ADOPTED and the facade falls back to the pinned call.
    relative_acceptance: float = 1.0e-8
    #: Optional absolute acceptance for quantities whose reference resolution is
    #: itself coarse (a bisection's terminal bracket).
    absolute_acceptance: float | None = None
    notes: str = ""
    metadata: dict[str, object] = field(default_factory=dict)

    @property
    def dimension(self) -> int:
        return len(self.axes)

    def contains(self, *values: float) -> bool:
        return all(axis.contains(value) for axis, value in zip(self.axes, values, strict=True))


_WATER_SATURATION_NOTE = (
    "IAPWS-95 coexistence: a 2x2 Newton on (rho_l, rho_v) seeded by two full "
    "density solves.  The upper domain bound covers the bisection excursion "
    "that saturation_temperature's first midpoint (460.078 K) drives into it."
)

_TSAT_NOTE = (
    "The pinned inversion is a 35-halving bisection on a [273.16, 646.996] K "
    "bracket returning the bracket midpoint, so its own resolution is a "
    "staircase of step <= 5e-9 K.  Its ACCURACY is coarser than that by four "
    "orders: every bracket decision reads saturation_pressure, which carries "
    "the D1 non-convergence noise, so a flipped comparison displaces the "
    "bracket.  The reference, not the Chebyshev decay, therefore sets the "
    "measured sup-norm, and the acceptance is stated ABSOLUTELY and justified "
    "by MATERIALITY against the certificate's 0.1 K temperature bound rather "
    "than by the reference's precision.  The pack also records the inversion "
    "residual - how well Psat(T_surrogate) reproduces the requested p against "
    "how well Psat(T_pinned) does - which is the direction-free comparison."
)

_DENSITY_NOTE = (
    "Newton root of p(T, rho) - p on the declared branch.  Only the root is "
    "surrogated; every property derived from it stays on the pinned Helmholtz "
    "table, so the ten state_Tp fields remain mutually consistent."
)

#: WHY THERE IS NO VAPOUR-BRANCH DENSITY SURROGATE (measured, 2026-09-02).
#:
#: A vapour surrogate needs a rectangular (T, p) box, and any box wide enough
#: to be useful reaches deep into p >> Psat(T), where the PINNED solve does not
#: return a root at all.  Measured on ``water.density(T, p, 'vapor')``:
#: at 290 K / 250 kPa (p/Psat = 130) it runs its full 200-iteration cap in
#: 42 ms and returns a wandering iterate of 98 kg/m3 -- neither the vapour root
#: (0.0019 kg/m3) nor the liquid one (998 kg/m3); at 330 K / 250 kPa
#: (p/Psat = 14.5) the same, 44 ms; only at 370 K / 250 kPa (p/Psat = 2.8) does
#: it converge, in 1.2 ms.  There is no function there to interpolate, so a
#: surrogate could only launder a non-convergence into a smooth-looking number.
#:
#: The cost of leaving it pinned is nil.  In the profiled two-interval march
#: 534 of 3,415 water density calls were vapour-branch and every one of them
#: came from inside ``saturation(T)`` at p = Psat_seed(T) -- which the
#: saturation surrogate deletes outright -- while ``state_Tp(..., 'vapor')``
#: accounted for 4 of 19,377 water and 4 of 33,232 hexane calls.
#:
#: WORDING REFRESHED 2026-09-03 (owed by S2 record Addendum A.6, discharged in
#: Addendum B).  The measurement above is left verbatim because it is what was
#: seen on 2026-09-02, but it now describes a failure mode the pinned authority
#: no longer exhibits.  The D2 vapour-density ceremony
#: (``docs/GT_PS2_D2_VAPOUR_DENSITY_REFUSAL_CEREMONY_2026-09-03.md``, commit
#: ``0f5faa7``) gave ``water.density(T, p, 'vapor')`` and its ``hexane.density``
#: twin a converge-or-raise discipline: a ``'vapor'`` request that converges onto
#: ``rho >= RHOC`` below ``TC`` is the liquid root, so the call now RAISES with
#: "the request is above the vapour spinodal" instead of returning a wandering
#: iterate.  The D3 fail-open ceremony
#: (``docs/GT_PS2_D3_PROPS_FAIL_OPEN_CEREMONY_2026-09-03.md``, commit
#: ``0992d5d``) added the symmetric liquid mirror.  The exclusion is therefore
#: unchanged but strictly stronger: it is no longer only this roster that
#: declines to fit the vapour branch above the spinodal -- the pinned surface
#: itself refuses to answer there.
VAPOUR_BRANCH_EXCLUSION = (
    "vapour-branch density is NOT surrogated: for p >> Psat(T) there is no "
    "vapour root to fit, and the pinned solve now says so -- since the D2 "
    "ceremony (2026-09-03) it RAISES above the vapour spinodal rather than "
    "returning the non-physical iterate measured on 2026-09-02 (200-iteration "
    "cap, 42-52 ms), and D3 added the symmetric liquid mirror; the measured "
    "call census shows the study lane loses nothing by leaving it pinned"
)

_VIRIAL_NOTE = (
    "The zero-density limit of alpha^r_delta and its analytic T-derivatives. "
    "Not a solve, but one full 56-term (water) residual evaluation per new "
    "temperature; the four derivatives share that single evaluation."
)


SURROGATE_SPECS: tuple[SurrogateSpec, ...] = (
    SurrogateSpec(
        name="water.saturation",
        fluid="water",
        kind="saturation",
        patch_targets=("_saturation_uncached",),
        axes=(Axis("T", "K", 290.0, 480.0),),
        outputs=("p", "rho_liquid", "rho_vapor", "dh_vap"),
        degree_ladder=(10, 12, 14, 16, 20, 24, 28),
        notes=_WATER_SATURATION_NOTE,
    ),
    SurrogateSpec(
        name="hexane.saturation",
        fluid="hexane",
        kind="saturation",
        patch_targets=("saturation",),
        axes=(Axis("T", "K", 290.0, 420.0),),
        outputs=("p", "rho_liquid", "rho_vapor", "dh_vap"),
        degree_ladder=(10, 12, 14, 16, 20, 24, 28),
        notes=(
            "Span-Wagner (2003) coexistence.  Unlike the water twin this call "
            "carries NO exact-bits memo in the pinned module, so every call is "
            "a full solve and the surrogate bites on all of them."
        ),
    ),
    SurrogateSpec(
        name="water.saturation_temperature",
        fluid="water",
        kind="scalar_1d",
        patch_targets=("_saturation_temperature_uncached",),
        axes=(Axis("p", "Pa", 2.0e4, 4.0e5),),
        outputs=("T",),
        degree_ladder=(10, 12, 14, 16, 20, 24),
        relative_target=1.0e-11,
        relative_acceptance=1.0e-9,
        absolute_acceptance=1.0e-3,
        notes=_TSAT_NOTE,
    ),
    SurrogateSpec(
        name="water.second_virial",
        fluid="water",
        kind="virial_1d",
        patch_targets=(
            "second_virial",
            "second_virial_dT",
            "second_virial_dT2",
            "second_virial_dT3",
        ),
        axes=(Axis("T", "K", 290.0, 420.0),),
        outputs=("second_virial", "second_virial_dT", "second_virial_dT2", "second_virial_dT3"),
        notes=_VIRIAL_NOTE,
    ),
    SurrogateSpec(
        name="hexane.second_virial",
        fluid="hexane",
        kind="virial_1d",
        patch_targets=(
            "second_virial",
            "second_virial_dT",
            "second_virial_dT2",
            "second_virial_dT3",
        ),
        axes=(Axis("T", "K", 290.0, 420.0),),
        outputs=("second_virial", "second_virial_dT", "second_virial_dT2", "second_virial_dT3"),
        notes=_VIRIAL_NOTE,
    ),
    SurrogateSpec(
        name="water.density:liquid",
        fluid="water",
        kind="density_2d",
        phase="liquid",
        patch_targets=("_density_uncached",),
        axes=(Axis("T", "K", 290.0, 370.0), Axis("p", "Pa", 500.0, 2.5e5)),
        outputs=("rho",),
        notes=_DENSITY_NOTE,
    ),
    SurrogateSpec(
        name="hexane.density:liquid",
        fluid="hexane",
        kind="density_2d",
        phase="liquid",
        patch_targets=("density",),
        axes=(Axis("T", "K", 290.0, 370.0), Axis("p", "Pa", 500.0, 2.5e5)),
        outputs=("rho",),
        notes=_DENSITY_NOTE,
    ),
)


def spec_by_name(name: str) -> SurrogateSpec:
    for spec in SURROGATE_SPECS:
        if spec.name == name:
            return spec
    raise KeyError(f"no surrogate spec named {name!r}")
