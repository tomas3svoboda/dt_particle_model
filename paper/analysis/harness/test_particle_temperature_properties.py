"""Caloric consistency, actual face-law wiring, refusal and restoration checks."""

from dataclasses import replace
import math

import pytest
import particle_temperature_properties as law
import positive_water_continuation as continuation
import finite_film_core2_study as study

core = study.core


@pytest.mark.parametrize("oil", [0.012, 0.0195, 0.03])
@pytest.mark.parametrize("temperature", [333.15, 378.15, 413.15])
def test_integrated_composite_caloric_derivative(oil, temperature):
    h = 0.001
    derivative = (
        law.composite_energy(temperature + h, 273.15, oil)
        - law.composite_energy(temperature - h, 273.15, oil)
    ) / (2 * h)
    assert derivative == pytest.approx(law.composite_capacity(temperature, oil), rel=1e-9)
    if oil == core.sp.W_O_REF:
        assert law.composite_energy(temperature, 273.15, oil) == core.sp.CP_NATIVE_DRY_MEAL * (
            temperature - 273.15
        )


def test_native_wet_capacity_inverse_cache_and_restoration():
    p = replace(core.wet_core.WetCoreParams(), w_o=0.013740458015267173)
    temperature, loading = 378.15, 0.2017
    before = core.wet_core.specific_energy(temperature, loading, p)
    function = core.wet_core.specific_heat_capacity
    with law.installed(4e-10):
        after = core.wet_core.specific_energy(temperature, loading, p)
        assert after != before
        h = 0.001
        derivative = (
            core.wet_core.specific_energy(temperature + h, loading, p)
            - core.wet_core.specific_energy(temperature - h, loading, p)
        ) / (2 * h)
        assert derivative == pytest.approx(
            core.wet_core.specific_heat_capacity(temperature, loading, p), rel=1e-9
        )
        recovered = core.wet_core.temperature_from_energy_density(p.rho_dm_p * after, loading, p)
        assert recovered == pytest.approx(temperature, abs=1e-7)
    assert core.wet_core.specific_heat_capacity is function
    assert core.wet_core.specific_energy(temperature, loading, p) == before


def test_real_cut_and_fully_dry_faces_use_local_temperature_once():
    original = core.q.install_engineering_case(
        "db_max_dapp_max", stress_multiplier=1, structural_conductivity_w_m_k=0.24
    ).configuration.cut
    reference = continuation.continued_configuration(original)
    T = 393.15
    p = reference.dry.pore
    args = (T, 0.2, T, 0.3, 0.001, 0.0, p, p, reference)
    old = core.cut._dry_flux(*args)
    grid = core.backend.uniform_grid(2, 0.0015)
    states = tuple(
        core.cp.evaluate_equilibrium(T, reference.dry.pressure_pa, y, p) for y in (0.2, 0.3)
    )
    boundary = core.ct.DirichletPoreBoundary(T, 0.4, "property wiring check")
    old_dry = core.ct._face_constitutive_data(grid, states, boundary, reference.dry)
    with law.installed(4e-10):
        config = continuation.continued_configuration(original)
        new = core.cut._dry_flux(*args[:-1], config)
        expected = law.effective_diffusivity(T, reference.dry.pressure_pa, 4e-10) / 6e-10
        assert (
            new.independent_water_flux_mol_m2_s / old.independent_water_flux_mol_m2_s
            == pytest.approx(expected, rel=1e-12)
        )
        assert new.retained_water_flux_mol_m2_s == old.retained_water_flux_mol_m2_s
        assert new.gas_entropy_production_w_m3_k >= 0
        reverse = core.cut._dry_flux(T, 0.3, T, 0.2, 0.001, 0.0, p, p, config)
        assert reverse.independent_water_flux_mol_m2_s == pytest.approx(
            -new.independent_water_flux_mol_m2_s, rel=1e-12
        )
        new_dry = core.ct._face_constitutive_data(grid, states, boundary, config.dry)
        for a, b in zip(
            old_dry.independent_water_fluxes_mol_m2_s[1:],
            new_dry.independent_water_fluxes_mol_m2_s[1:],
        ):
            assert b / a == pytest.approx(expected, rel=1e-12)
        assert new_dry.retained_water_fluxes_mol_m2_s == old_dry.retained_water_fluxes_mol_m2_s


def test_bad_property_inputs_refuse():
    for bad in (0.0, -1.0, math.inf, math.nan):
        with pytest.raises(ValueError):
            law.effective_diffusivity(bad, 101325, 4e-10)
        with pytest.raises(ValueError):
            law.oil_heat_capacity(bad)
    assert law.oil_heat_capacity(378.15) == pytest.approx(2158.45)
