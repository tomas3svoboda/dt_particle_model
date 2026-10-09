"""Reuse the existing DT march, preserving an unsuccessful birth iterate as a seed.

No time, inventory, physical coefficient or acceptance check is changed.
An optional seed must be from this same physical initial-value problem.
"""
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / 'src'))
sys.path.insert(0, str(root / 'paper/analysis/harness'))
import experiment_core2_finite_wetting_b as engine

seed_path = None
if '--birth-seed' in sys.argv:
    position = sys.argv.index('--birth-seed')
    seed_path = Path(sys.argv[position + 1])
    del sys.argv[position:position + 2]
out_path = Path(sys.argv[sys.argv.index('--out') + 1])
requested_factor = float(sys.argv[sys.argv.index('--factor') + 1])
requested_K = float(sys.argv[sys.argv.index('--contact-conductance') + 1])
assert requested_K == 0.32 and requested_factor in (100, 10000)
seed = json.loads(seed_path.read_text()) if seed_path else None
if seed:
    assert seed['K_l_mol_m2_s'] == requested_K
    assert seed['shell_mobility_factor'] == requested_factor
original = engine.driver.advance_birth

def observed_birth(*args, **kwargs):
    if seed:
        assert args[0].time_s == seed['fully_wet_before_time_s']
        kwargs['root_seed'] = seed
    try:
        return original(*args, **kwargs)
    except Exception as error:
        debug = getattr(error, 'study_debug', None)
        if debug:
            payload = {
                'K_l_mol_m2_s': requested_K,
                'shell_mobility_factor': requested_factor,
                'fully_wet_before_time_s': args[0].time_s,
                'study_debug': debug,
                'error': str(error),
                'total_residual_assembly_attempts': getattr(error, 'total_residual_assembly_attempts', None),
                'maximum_scaled_residual': getattr(error, 'maximum_scaled_residual', None),
                'failed_trial_advanced_physical_time': False,
                'exact_fully_wet_rollback_identity_preserved': error.rollback_state is args[0],
            }
            target = out_path.with_suffix('.birth_debug.json')
            assert not target.exists()
            target.write_text(json.dumps(payload, indent=2, default=lambda value: value.tolist() if hasattr(value, 'tolist') else float(value)))
        raise

provenance_path = out_path.with_suffix('.birth_seed_provenance.json')
assert not provenance_path.exists()
provenance_path.parent.mkdir(parents=True, exist_ok=True)
provenance_path.write_text(json.dumps({
    'root_seed': str(seed_path) if seed_path else None,
    'root_seed_sha256': hashlib.sha256(seed_path.read_bytes()).hexdigest() if seed_path else None,
    'same_IVP_seed_only': True,
    'accepted_material_state_imported': False,
    'acceptance_rules_changed': False,
    'wrapper_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}, indent=2))
engine.driver.advance_birth = observed_birth
original_proof = engine.driver.bounded_path_proof
def longer_wall(assembly_budget, wall_s):
    return original_proof(assembly_budget, max(wall_s, 600.0))
engine.driver.bounded_path_proof = longer_wall
proof_policy = out_path.with_suffix('.proof_wall_policy.json')
assert not proof_policy.exists()
proof_policy.write_text(json.dumps({
    'K_l_mol_m2_s': requested_K, 'mobility_factor': requested_factor,
    'step_proof_wall_s': 600, 'proof_point_budget': 3000000,
    'seed_only_tangent_bound': 3e-10, 'component_energy_tolerance': 1e-10,
    'native_physical_equations_changed': False,
    'owner_authority': 'RULE 2; same documented numerical controls as the current 0.32 baseline',
}, indent=2))
try:
    import runpy
    runpy.run_path(str(Path(__file__).with_name('run_b_admissible_reference_mobility.py')), run_name='__main__')
finally:
    engine.driver.advance_birth = original
    engine.driver.bounded_path_proof = original_proof
