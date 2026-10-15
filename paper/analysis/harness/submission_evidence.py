"""Restore the paper's compressed evidence without changing archived bytes."""

import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[3]
_ready = False


def ensure_evidence():
    global _ready
    if _ready:
        return
    description = json.loads((ROOT / "data/evidence.json").read_text(encoding="utf-8"))
    path = ROOT / description["archive"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == description["sha256"], "Evidence archive hash mismatch"
    with tarfile.open(path, "r:xz") as archive:
        manifest = archive.extractfile("EVIDENCE_MANIFEST.sha256").read().decode("utf-8")
        expected = dict((name, digest) for digest, name in (line.split("  ", 1) for line in manifest.splitlines()))
        members = {item.name: item for item in archive.getmembers() if item.name != "EVIDENCE_MANIFEST.sha256"}
        assert set(members) == set(expected)
        for name, member in members.items():
            target = (ROOT / name).resolve()
            if not member.isfile() or not target.is_relative_to(ROOT.resolve()):
                raise ValueError("Evidence member is not a regular file inside the repository")
            if target.exists():
                if hashlib.sha256(target.read_bytes()).hexdigest() != expected[name]:
                    raise ValueError(f"Preserve the modified extracted file before retrying: {name}")
                continue
            content = archive.extractfile(member).read()
            assert hashlib.sha256(content).hexdigest() == expected[name], name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
    _ready = True


if __name__ == "__main__":
    ensure_evidence()
    print("Paper evidence extracted and verified; extracted files are Git-ignored")
