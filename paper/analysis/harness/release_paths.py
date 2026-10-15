"""Resolve archived producer paths inside the public particle repository.

Raw run records and their SHA-256 links are preserved byte for byte. An old
absolute producer path is provenance, not an instruction to read that host.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def resolve_record_path(value):
    text = str(value).replace("\\", "/")
    if "paper/analysis/" in text:
        text = "paper/analysis/" + text.rsplit("paper/analysis/", 1)[1]
    path = Path(text)
    candidate = path if path.is_absolute() else ROOT / path
    candidate = candidate.resolve()
    if not candidate.is_relative_to(ROOT.resolve()):
        raise ValueError(f"Record path is outside the public repository: {value}")
    return candidate
