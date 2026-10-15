"""Collect the paper harness's physical-property tests in the standard suite."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "paper/analysis/harness"))
from test_particle_temperature_properties import (  # noqa: E402
    test_bad_property_inputs_refuse as test_bad_property_inputs_refuse,
    test_integrated_composite_caloric_derivative as test_integrated_composite_caloric_derivative,
    test_native_wet_capacity_inverse_cache_and_restoration as test_native_wet_capacity_inverse_cache_and_restoration,
    test_real_cut_and_fully_dry_faces_use_local_temperature_once as test_real_cut_and_fully_dry_faces_use_local_temperature_once,
)
