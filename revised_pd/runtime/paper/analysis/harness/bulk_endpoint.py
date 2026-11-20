"""Exact pure-hexane *bulk* carrier for the shared finite binary film.

This does not alter a pore-gas domain, water law, surface composition or flux.
Native film partition/profile equations already permit bulk mass fractions 0
and 1. The native fast-mass boundary requires an interior bulk because it also
uses bulk composition as pore-surface composition. The study adapter resolves
the latter independently, so the former can be exactly pure.
"""

from dataclasses import dataclass, fields
from contextlib import contextmanager
import math

import finite_film_core2_study as study

core = study.core
INHERITED_FANER_COMMON_EXPOSURE_MULTIPLIER = 0.14578543636088787


@contextmanager
def declared_faner_film_multiplier(multiplier):
    """Admit a stated common exposure input; keep the source pair unchanged.

    Reuses the predecessor's named Faner boundary extension in a scoped study.
    Native B.7--B.10 regeneration, joint scaling, and validity-ratio checks are
    untouched. Only this stated multiplier joins the existing allowed set.
    """
    if isinstance(multiplier, bool) or not math.isfinite(multiplier) or not 0 < multiplier <= 4:
        raise ValueError("declared Faner common film multiplier must lie in (0,4]")
    original_private = core.efr._ALLOWED_STRESS_MULTIPLIERS
    original_public = core.efr.ALLOWED_STRESS_MULTIPLIERS
    if (
        multiplier not in original_private
        and multiplier != INHERITED_FANER_COMMON_EXPOSURE_MULTIPLIER
    ):
        raise ValueError(
            "only the inherited owner-approved Faner common multiplier can extend the native whitelist"
        )
    admitted = original_private | frozenset((float(multiplier),))
    audit = {
        "declared_multiplier": multiplier,
        "source_pair_edited": False,
        "native_joint_heat_mass_scaling": True,
        "native_pair_regeneration_check_preserved": True,
        "native_validity_ratio_check_preserved": True,
        "approved_identified_boundary_input": multiplier
        == INHERITED_FANER_COMMON_EXPOSURE_MULTIPLIER,
        "falling_rate_fit_performed": False,
    }
    core.efr._ALLOWED_STRESS_MULTIPLIERS = admitted
    core.efr.ALLOWED_STRESS_MULTIPLIERS = admitted
    try:
        yield audit
    finally:
        core.efr._ALLOWED_STRESS_MULTIPLIERS = original_private
        core.efr.ALLOWED_STRESS_MULTIPLIERS = original_public
        audit["native_multiplier_sets_restored"] = (
            core.efr._ALLOWED_STRESS_MULTIPLIERS is original_private
            and core.efr.ALLOWED_STRESS_MULTIPLIERS is original_public
        )


def native_validation_probe(boundary):
    """Validate every native field without evaluating a surrogate bath.

    The interior composition belongs only to this unused structural-validation
    object. Actual objects and every physical calculation retain y_bulk = 1.
    """
    kwargs = {
        field.name: getattr(boundary, field.name)
        for field in fields(core.ct.ReducedFilmPoreBoundary)
        if field.init
    }
    kwargs["y_hexane"] = 0.5
    return core.ct.ReducedFilmPoreBoundary(**kwargs)


@dataclass(frozen=True)
class PureHexaneBulkFilmBoundary(core.ct.ReducedFilmPoreBoundary):
    def __post_init__(self):
        if self.y_hexane != 1.0:
            raise ValueError("this named bulk endpoint requires exactly y_hexane=1")
        native_validation_probe(self)


@dataclass(frozen=True)
class PureHexaneFiniteFilmBoundary(study.FiniteFilmStudyBoundary):
    def __post_init__(self):
        if self.y_hexane != 1.0:
            raise ValueError("this named study carrier requires exactly y_hexane=1")
        if not isinstance(self.film_boundary, PureHexaneBulkFilmBoundary):
            raise TypeError("the exact pure bulk film carrier is required")
        # Reuse the existing study's complete structural validation. The probe
        # does not flow into any physics; self retains exact endpoint fields.
        study.FiniteFilmStudyBoundary(
            self.temperature_k,
            0.5,
            self.label,
            physically_qualifying=self.physically_qualifying,
            film_boundary=native_validation_probe(self.film_boundary),
            accepted_beta_limit=self.accepted_beta_limit,
            actual_surface_path=self.actual_surface_path,
            numerical_cap_seam=self.numerical_cap_seam,
        )
        if self.temperature_k != self.film_boundary.temperature_k:
            raise ValueError("bulk film and integration carrier temperatures differ")
        if self.physically_qualifying or self.film_boundary.physically_qualifying:
            raise ValueError("the finite-film integration carrier is not qualifying")
        if not math.isfinite(self.temperature_k):
            raise ValueError("bulk temperature must be finite")


@contextmanager
def pure_bulk_domain_adapter():
    """Admit only the named external bulk, preserving native pore checks.

    Native cut assembly normally treats its supplied boundary as a pore state.
    In this study the named object is the external bath, and the existing film
    adapter supplies a separately admitted binary surface. Check bulk gas and
    carrier fields here; route every other boundary to the unchanged validator.
    """
    original = core.cut._validate_dry_boundary
    audit = {
        "named_bulk_validations": 0,
        "native_pore_boundary_validations": 0,
        "bulk_y_hexane": 1.0,
        "surface_domain_relaxed": False,
    }

    def validate_bulk_or_native(boundary, config):
        if type(boundary) is not PureHexaneFiniteFilmBoundary:
            audit["native_pore_boundary_validations"] += 1
            return original(boundary, config)
        boundary.__post_init__()
        if not boundary.actual_surface_path:
            raise ValueError("pure bulk requires the resolved actual binary surface path")
        if boundary.film_boundary.binary_gas_interaction_k_wh != config.dry.pore.k_wh:
            raise ValueError("pure bulk and native pore interaction parameters differ")
        # Evaluate the actual free-gas bath, not a sorbing pore or a fictitious
        # interior-composition replacement for the experimental pure inlet.
        core.bg.state(
            boundary.temperature_k,
            config.dry.pressure_pa,
            0.0,
            1.0,
            k_wh=config.dry.pore.k_wh,
        )
        audit["named_bulk_validations"] += 1

    core.cut._validate_dry_boundary = validate_bulk_or_native
    try:
        yield audit
    finally:
        core.cut._validate_dry_boundary = original
        audit["native_validator_restored"] = core.cut._validate_dry_boundary is original


def make_pure_bulk_boundary(
    pair, *, temperature_k, forcing_timescale_s, beta_limit=1.0, stress_multiplier=1.0
):
    """Use native film parameters, replacing only its integration carrier."""
    template = core.reduced_film_boundary(
        temperature_k=temperature_k,
        y_hexane=0.5,
        pair=pair,
        forcing_timescale_s=forcing_timescale_s,
        label="native parameter validation only; no physical evaluation",
    )
    kwargs = {
        field.name: getattr(template, field.name)
        for field in fields(core.ct.ReducedFilmPoreBoundary)
        if field.init
    }
    kwargs.update(
        y_hexane=1.0,
        stress_multiplier=stress_multiplier,
        label="EXACT pure-hexane bulk; active binary interior/surface water",
    )
    with declared_faner_film_multiplier(stress_multiplier):
        scaled_mass = stress_multiplier * pair.binary_mass_transfer_m_s
        ratio = pair.state.bulk_binary_diffusivity_m2_s / scaled_mass**2 / forcing_timescale_s
        if ratio <= 0.01:
            kwargs["fidelity_band"] = core.efr.FidelityBand.NOMINAL
        elif ratio <= 0.10:
            kwargs["fidelity_band"] = core.efr.FidelityBand.BOUNDED_ANOMALY
        else:
            raise core.efr.ReducedFilmValidityError(
                "declared Faner scaled film exceeds native validity-ratio band"
            )
        bulk = PureHexaneBulkFilmBoundary(**kwargs)
        return PureHexaneFiniteFilmBoundary(
            temperature_k,
            1.0,
            "exact pure bulk with resolved finite binary surface",
            film_boundary=bulk,
            accepted_beta_limit=beta_limit,
            actual_surface_path=True,
            numerical_cap_seam=True,
        )
