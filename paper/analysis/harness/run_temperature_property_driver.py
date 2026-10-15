"""Run an existing native paper CLI under an explicit temperature-property law."""

import argparse
import hashlib
import json
from pathlib import Path
import runpy
import sys

import particle_temperature_properties as properties


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Db-ref", type=float, required=True)
    parser.add_argument("--constant-oil", action="store_true")
    parser.add_argument("--constant-D", action="store_true")
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    args.driver = args.driver.resolve()
    argv = args.arguments
    if argv and argv[0] == "--":
        argv = argv[1:]
    law = properties.declaration(args.Db_ref, not args.constant_oil, not args.constant_D)
    properties.effective_diffusivity(
        properties.REFERENCE_T_K, properties.REFERENCE_P_PA, args.Db_ref
    )
    if "--resume" in argv:
        checkpoint = Path(argv[argv.index("--resume") + 1])
        parent_path = checkpoint.with_name(checkpoint.name.replace(".accepted.pkl", ".json"))
        parent = json.loads(parent_path.read_text(encoding="utf-8"))
        if json.dumps(parent.get("temperature_property_variant"), sort_keys=True) != json.dumps(
            law, sort_keys=True
        ):
            raise ValueError(
                "checkpoint property law differs; a full-feed case requires the same declared law on every segment"
            )
    out = Path(argv[argv.index("--out") + 1])
    provenance = out.with_name("properties_" + out.name)
    if out.exists() or provenance.exists():
        raise ValueError("refusing to overwrite study evidence")
    out.parent.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__).resolve(), Path(properties.__file__).resolve(), args.driver.resolve()]
    entry = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    with properties.installed(
        args.Db_ref, oil_law=not args.constant_oil, temperature_law=not args.constant_D
    ) as law:
        provenance.write_text(
            json.dumps({**law, "source_sha256": entry, "argv": argv}, indent=2), encoding="utf-8"
        )
        sys.argv = [str(args.driver), *argv]
        try:
            runpy.run_path(str(args.driver), run_name="__main__")
        finally:
            if out.is_file():
                data = json.loads(out.read_text(encoding="utf-8"))
                data["temperature_property_variant"] = law
                data["property_variant_source_sha256"] = entry
                data["property_variant_sources_unchanged"] = all(
                    hashlib.sha256(p.read_bytes()).hexdigest() == entry[str(p)] for p in sources
                )
                out.write_text(json.dumps(data, indent=2, allow_nan=True), encoding="utf-8")


if __name__ == "__main__":
    main()
