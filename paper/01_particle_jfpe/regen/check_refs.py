r"""Check that every \ref/\eqref target in the main text and the supplement
resolves to a \label defined in the same document, and (2026-10-01) that every
cross-document \ref* resolves in the other one: \ref*{S-...} in the main text
against the supplement's labels, \ref*{M-...} in the supplement against the
main text's (the xr prefixes of main.tex and supplementary.tex).

Run:
  .venv\Scripts\python.exe paper\01_particle_jfpe\regen\check_refs.py
"""

from __future__ import annotations

import pathlib
import re
import argparse

HERE = pathlib.Path(__file__).resolve().parent
PAPER = HERE.parent

MAIN_FILES = [
    "main.tex",
    "sec_intro_reshape.tex",
    "sec_formulation.tex",
    "sec_discretization_reshape.tex",
    "sec_verification_reshape.tex",
    "sec_validation_reshape.tex",
    "sec_discussion_reshape.tex",
]
MOVED = [
    "moved_full/" + n
    for n in (
        "sec_intro.tex",
        "sec_discretization.tex",
        "sec_verification.tex",
        "sec_validation.tex",
        "sec_results.tex",
        "sec_discussion.tex",
    )
]
SUPP_FILES = [
    "supplementary_submission.tex",
    "supplementary_submission_detail.tex",
    "supplementary_submission_protocol.tex",
    "supplementary_submission_evidence.tex",
    "supplementary_submission_dtmarch.tex",
    "supplementary_submission_finite_contact.tex",
]
ARCHIVE_FILES = ["supplementary.tex", *MOVED]

LABEL = re.compile(r"\\label\{([^}]*)\}")
REF = re.compile(r"\\(?:eq)?ref\{([^}]*)\}")
XREF = re.compile(r"\\ref\*\{([SM])-([^}]*)\}")


def scan(names: list[str]) -> tuple[set[str], list[tuple[str, str]]]:
    labels: set[str] = set()
    refs: list[tuple[str, str]] = []
    for name in names:
        text = (PAPER / name).read_text(encoding="utf-8")
        text = re.sub(r"(?<!\\)%.*", "", text)
        labels.update(LABEL.findall(text))
        refs.extend((name, r) for r in REF.findall(text))
    return labels, refs


def report(title: str, names: list[str], fallback: list[str] | None = None) -> int:
    labels, refs = scan(names)
    fallback_labels = scan(fallback)[0] if fallback else set()
    missing = [(f, r) for f, r in refs if r not in labels | fallback_labels]
    print(f"{title}: {len(labels)} labels, {len(refs)} references, " f"{len(missing)} unresolved")
    for f, r in missing:
        print(f"  UNRESOLVED  {r}   (in {f})")
    return len(missing)


def cross(title: str, names: list[str], prefix: str, other: list[str]) -> int:
    labels, _ = scan(other)
    refs = []
    for name in names:
        text = re.sub(r"(?<!\\)%.*", "", (PAPER / name).read_text(encoding="utf-8"))
        refs.extend((name, k) for p, k in XREF.findall(text) if p == prefix)
    missing = [(f, r) for f, r in refs if r not in labels]
    print(f"{title}: {len(refs)} references, {len(missing)} unresolved")
    for f, r in missing:
        print(f"  UNRESOLVED  {prefix}-{r}   (in {f})")
    return len(missing)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-archives",
        action="store_true",
        help="also inspect historical documents that are not built",
    )
    args = parser.parse_args()
    bad = report("main document", MAIN_FILES)
    bad += report("submitted supplement", SUPP_FILES, MAIN_FILES)
    bad += cross("main document into the supplement (S-)", MAIN_FILES, "S", SUPP_FILES)
    bad += cross("supplement into the main document (M-)", SUPP_FILES, "M", MAIN_FILES)
    if args.include_archives:
        bad += report("internal full-length supplement", ARCHIVE_FILES, MAIN_FILES)
        bad += cross(
            "internal full-length supplement into the main document (M-)",
            ARCHIVE_FILES,
            "M",
            MAIN_FILES,
        )
        bad += cross(
            "archived passages into the full-length supplement (S-)", MOVED, "S", ARCHIVE_FILES
        )
    if bad == 0:
        print("all checked references resolve")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
