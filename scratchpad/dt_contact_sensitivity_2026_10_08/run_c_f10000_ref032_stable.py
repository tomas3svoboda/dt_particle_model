"""Continue C using the same cancellation-resistant liquid force as B."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "paper/analysis/harness"))
sys.path.insert(0, str(Path(__file__).parent))
import experiment_core2_radial_fast_wetting as engine
from liquid_gradient_f10000_ref032 import installed_function

out = Path(sys.argv[sys.argv.index("--out") + 1])
original = engine.liquid_connected_fast_path
replacement, policy = installed_function(engine)
policy_path = out.with_suffix(".stable_gradient_policy.json")
assert not policy_path.exists()
policy_path.write_text(json.dumps(policy, indent=2))
engine.liquid_connected_fast_path = replacement
try:
    engine.main()
finally:
    engine.liquid_connected_fast_path = original
