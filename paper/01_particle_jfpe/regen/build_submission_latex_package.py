"""Build the minimal LaTeX package of the submitted manuscript and supplement.

The package holds exactly what the two documents compile from: the two master
files, the section files they input, the shared preamble and macros, the
bibliography pruned to the cited entries (internal ``annote`` and ``file`` fields
removed), the figures the documents include, and a build script. Full-line
LaTeX comments are dropped from the copied sources. Nothing else is shipped.

Usage, from paper/01_particle_jfpe:
    python regen/build_submission_latex_package.py --out submission_latex_<date>
The package is then compiled in place (pdflatex + bibtex) to check that it is
self-contained, and zipped beside the directory.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
SHARED = PAPER.parent / "shared"

MAIN = "main.tex"
SUPP = "supplementary_submission.tex"
CRLF = "\r\n"

CITE_RE = re.compile(r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\]){0,2}\{([^}]*)\}")
INPUT_RE = re.compile(r"\\input\{([^}]*)\}")
GRAPHICS_RE = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}")


def read(path: Path) -> str:
    return io.open(path, encoding="utf-8", newline="").read()


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    io.open(path, "w", encoding="utf-8", newline="").write(text)


def strip_comment_lines(text: str) -> str:
    """Drop lines whose first non-blank character is %, keeping line endings."""
    nl = CRLF if CRLF in text else "\n"
    kept = [line for line in text.split(nl) if not line.lstrip().startswith("%")]
    return nl.join(kept)


def flatten_paths(text: str) -> str:
    text = text.replace(r"\input{../shared/preamble.tex}", r"\input{preamble.tex}")
    text = text.replace(r"\input{../shared/macros.tex}", r"\input{macros.tex}")
    text = text.replace(r"\bibliography{../shared/refs}", r"\bibliography{refs}")
    return text


def closure(master: Path) -> list[Path]:
    """Every .tex file reached by \\input from the master, in order, master first."""
    seen: list[Path] = []

    def visit(path: Path) -> None:
        if path in seen:
            return
        seen.append(path)
        for name in INPUT_RE.findall(read(path)):
            target = (path.parent / name).resolve()
            if target.exists():
                visit(target)

    visit(master.resolve())
    return seen


def bib_entries(text: str) -> list[tuple[str, str, str]]:
    """(type, key, full entry text) for every @entry, by brace matching."""
    out = []
    for m in re.finditer(r"@(\w+)\s*\{", text):
        kind = m.group(1).lower()
        i = m.end()
        depth = 1
        while depth and i < len(text):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
            i += 1
        body = text[m.start():i]
        if kind in ("comment",):
            continue
        key = ""
        if kind not in ("string", "preamble"):
            key = body[m.end() - m.start():].split(",", 1)[0].strip()
        out.append((kind, key, body))
    return out


def drop_fields(entry: str, names=("annote", "file")) -> str:
    """Remove the named fields (brace-delimited values) from one entry."""
    for name in names:
        while True:
            m = re.search(r",\s*%s\s*=\s*\{" % name, entry, re.I)
            if not m:
                break
            i = m.end()
            depth = 1
            while depth and i < len(entry):
                if entry[i] == "{":
                    depth += 1
                elif entry[i] == "}":
                    depth -= 1
                i += 1
            entry = entry[: m.start()] + entry[i:]
    return entry


def pruned_bib(cited: set[str]) -> tuple[str, list[str]]:
    text = read(SHARED / "refs.bib")
    nl = CRLF if CRLF in text else "\n"
    entries = bib_entries(text)
    by_key = {k: (t, b) for t, k, b in entries if k}
    wanted = set(cited)
    # follow crossref fields once
    for key in list(wanted):
        if key in by_key:
            m = re.search(r"crossref\s*=\s*\{([^}]*)\}", by_key[key][1], re.I)
            if m:
                wanted.add(m.group(1).strip())
    missing = sorted(k for k in wanted if k not in by_key)
    parts = [b for t, k, b in entries if t in ("string", "preamble")]
    for t, k, b in entries:
        if k in wanted:
            parts.append(drop_fields(b))
    return (nl + nl).join(parts) + nl, missing


BUILD_BAT = """@echo off
rem Build the manuscript (main.pdf) and the Supplementary Material
rem (supplementary_submission.pdf). The two documents cross-reference each
rem other through the xr package, so each is compiled after the other's .aux
rem exists. Requires pdflatex, bibtex and the apacite, natbib, siunitx and xr
rem packages (any current TeX Live or MiKTeX).
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
bibtex supplementary_submission
pdflatex -interaction=nonstopmode supplementary_submission.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
"""

BUILD_SH = """#!/bin/sh
# Same sequence as build.bat, for Unix-like systems.
set -e
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
bibtex supplementary_submission
pdflatex -interaction=nonstopmode supplementary_submission.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode supplementary_submission.tex
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args()
    out = args.out if args.out.is_absolute() else (PAPER / args.out)
    if out.exists():
        raise FileExistsError('Use a fresh output folder; preserve existing submission evidence')
    out.mkdir(parents=True)

    tex_files: list[Path] = []
    for master in (PAPER / MAIN, PAPER / SUPP):
        for path in closure(master):
            if path not in tex_files:
                tex_files.append(path)

    cited: set[str] = set()
    figures: set[str] = set()
    shipped: list[Path] = []
    for src in tex_files:
        text = read(src)
        for group in CITE_RE.findall(text):
            cited.update(k.strip() for k in group.split(",") if k.strip())
        figures.update(GRAPHICS_RE.findall(text))
        if src.parent == SHARED:
            dest = out / src.name
            write(dest, strip_comment_lines(text))
        else:
            dest = out / src.name
            write(dest, flatten_paths(strip_comment_lines(text)))
        shipped.append(dest)

    bib, missing = pruned_bib(cited)
    if missing:
        raise SystemExit(f"cited keys missing from refs.bib: {missing}")
    write(out / "refs.bib", bib)
    shipped.append(out / "refs.bib")

    for fig in sorted(figures):
        src = PAPER / fig
        if not src.exists():
            raise SystemExit(f"figure not found: {src}")
        dest = out / fig
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        shipped.append(dest)

    write(out / "build.bat", BUILD_BAT.replace("\n", CRLF))
    write(out / "build.sh", BUILD_SH)
    shipped += [out / "build.bat", out / "build.sh"]

    readme = [
        "Submission LaTeX package",
        "",
        "main.tex                      manuscript (inputs the sec_*.tex files)",
        "supplementary_submission.tex  Supplementary Material (inputs the",
        "                              supplementary_submission_*.tex files)",
        "preamble.tex, macros.tex      shared preamble and macros",
        "core2_particle_results.tex    result macros used by both documents",
        "refs.bib                      bibliography (cited entries only)",
        "figures/                      the figures the documents include",
        "build.bat, build.sh           compile sequence (pdflatex + bibtex; the",
        "                              documents cross-reference each other with xr)",
        "",
        "Cited keys: %d. Figures: %d." % (len(cited), len(figures)),
    ]
    write(out / "README.txt", CRLF.join(readme) + CRLF)
    shipped.append(out / "README.txt")

    manifest = []
    for path in sorted(shipped):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest.append(f"{digest}  {path.relative_to(out).as_posix()}")
    write(out / "MANIFEST.sha256", CRLF.join(manifest) + CRLF)

    print(f"package: {out}")
    print(f"  tex files {len(tex_files)}, cited keys {len(cited)}, figures {len(figures)}")

    if not args.no_build:
        result = subprocess.run(["cmd", "/c", "build.bat"], cwd=out, capture_output=True, text=True)
        log = out / "build.log.txt"
        write(log, result.stdout + result.stderr)
        undefined = re.findall(r"(?:Reference|Citation) `[^']*' .*undefined", result.stdout)
        errors = re.findall(r"^! .*", result.stdout, re.M)
        print(f"  build exit {result.returncode}; errors {len(errors)}; undefined in all passes {len(undefined)}")
        for pdf in ("main.pdf", "supplementary_submission.pdf"):
            p = out / pdf
            print(f"  {pdf}: {'ok' if p.exists() else 'MISSING'}")
        # remove build by-products so the package stays minimal; keep the two PDFs
        for p in out.iterdir():
            if p.suffix in (".aux", ".bbl", ".blg", ".log", ".out", ".toc") or p.name == "build.log.txt":
                p.unlink()

    zip_path = out.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(out.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(out.parent).as_posix())
    print(f"  zip: {zip_path} ({zip_path.stat().st_size // 1024} kB)")


if __name__ == "__main__":
    main()
