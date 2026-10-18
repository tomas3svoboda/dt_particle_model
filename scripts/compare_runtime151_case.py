"""Compare a fresh qualification case directory with its runtime-151 reference.

Every JSON and JSONL record is compared value by value (floats exactly, no tolerance).
Differences are classified so that only the expected provenance fields are excused:

  wall_clock   wall-clock durations and timestamps (keys naming wall, elapsed, perf, utc, ...)
  path         strings that name the producer's or the rerun's own directories
  git_stamp    the era stamp of the checkout (repo_head, repo_branch)
  source_hash  source-digest maps (the shipped tree is a subset of the runtime tree)
  other        everything else -- a reproduction defect unless explained

Pickled checkpoints that are not byte-identical are unpickled (classes from src/ and
the harness) and compared value by value under the same classes.

    python scripts/compare_runtime151_case.py <fresh case dir> <reference case dir> [--json out.json]
"""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

WALL = re.compile(r"(wall|elapsed|perf_counter|_utc|utc_|timestamp|started|finished|seconds_since|clock_s|cpu_s|runtime_s)", re.I)
HASHKEY = re.compile(r"(sha256|source_hash|hashes|source_pins)", re.I)
GITSTAMP = re.compile(r"^repo_(head|branch|dirty|status)$")
PATHLIKE = re.compile(r"([A-Za-z]:[\\/]|[\\/](paper|reruns|scratchpad)[\\/]|\.pkl$|\.json$|\.jsonl$)")


def classify(path, a, b):
    joined = "/".join(map(str, path))
    if path and isinstance(path[-1], str) and GITSTAMP.match(path[-1]):
        return "git_stamp"
    if any(HASHKEY.search(str(k)) for k in path if isinstance(k, str)):
        return "source_hash"
    if any(WALL.search(str(k)) for k in path if isinstance(k, str)):
        return "wall_clock"
    if isinstance(a, str) and isinstance(b, str) and (PATHLIKE.search(a) or PATHLIKE.search(b)):
        return "path"
    if any(isinstance(k, str) and PATHLIKE.search(k) for k in path):
        return "source_hash" if "sha" in joined.lower() else "path"
    return "other"


def same_float(a, b):
    return a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b))


def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                cls = classify(path + [k], a.get(k), b.get(k))
                out.append((cls, path + [k], "missing" if k not in a else a[k], "missing" if k not in b else b[k]))
            else:
                walk(a[k], b[k], path + [k], out)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append((classify(path, a, b), path + ["<len>"], len(a), len(b)))
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, path + [i], out)
        return
    if type(a) is not type(b) or not same_float(a, b):
        out.append((classify(path, a, b), path, a, b))


def plain(obj, depth=0):
    """Unpickled object -> JSON-like tree (dataclasses, slots, numpy arrays, enums)."""
    if depth > 60:
        return repr(obj)
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, bytes):
        return {"<bytes>": hashlib.sha256(obj).hexdigest()}
    if isinstance(obj, dict):
        return {str(k): plain(v, depth + 1) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set, frozenset)):
        items = sorted(obj, key=repr) if isinstance(obj, (set, frozenset)) else obj
        return [plain(v, depth + 1) for v in items]
    if hasattr(obj, "dtype") and hasattr(obj, "tolist"):
        return {"<array>": str(obj.dtype), "shape": list(getattr(obj, "shape", ())), "values": plain(obj.tolist(), depth + 1)}
    if hasattr(obj, "value") and type(obj).__class__.__name__ == "EnumType":
        return {"<enum>": f"{type(obj).__name__}.{obj.name}"}
    state = {}
    if hasattr(obj, "__dict__"):
        state.update(vars(obj))
    for cls in type(obj).__mro__:
        for slot in getattr(cls, "__slots__", ()):
            if isinstance(slot, str) and hasattr(obj, slot):
                state[slot] = getattr(obj, slot)
    return {"<type>": f"{type(obj).__module__}.{type(obj).__qualname__}", **{k: plain(v, depth + 1) for k, v in state.items()}}


def load_pickle(p):
    import pickle
    return plain(pickle.loads(p.read_bytes()))


def load_jsonl(p):
    return [json.loads(line) for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]


def compare(fresh, ref):
    report = {"fresh": str(fresh), "reference": str(ref), "files": {}}
    names = sorted({p.name for p in fresh.iterdir() if p.is_file()} | {p.name for p in ref.iterdir() if p.is_file()})
    for name in names:
        fa, fb = fresh / name, ref / name
        if not fa.is_file() or not fb.is_file():
            report["files"][name] = {"status": "missing_in_fresh" if not fa.is_file() else "missing_in_reference"}
            continue
        entry = {"bytes_identical": fa.read_bytes() == fb.read_bytes(),
                 "sha256_fresh": hashlib.sha256(fa.read_bytes()).hexdigest(),
                 "sha256_reference": hashlib.sha256(fb.read_bytes()).hexdigest()}
        if name.endswith((".json", ".jsonl", ".pkl")) and not entry["bytes_identical"] or name.endswith((".json", ".jsonl")):
            loader = load_pickle if name.endswith(".pkl") else (load_jsonl if name.endswith(".jsonl") else lambda p: json.loads(p.read_text(encoding="utf-8")))
            a, b = loader(fa), loader(fb)
            diffs = []
            walk(a, b, [], diffs)
            counts = {}
            for cls, *_ in diffs:
                counts[cls] = counts.get(cls, 0) + 1
            entry["differences_by_class"] = counts
            entry["other_examples"] = [{"path": "/".join(map(str, p)), "fresh": repr(x)[:200], "reference": repr(y)[:200]}
                                       for cls, p, x, y in diffs if cls == "other"][:40]
            entry["excused_examples"] = [{"class": cls, "path": "/".join(map(str, p)), "fresh": repr(x)[:120], "reference": repr(y)[:120]}
                                         for cls, p, x, y in diffs if cls != "other"][:12]
            if name.endswith(".jsonl"):
                entry["records"] = [len(a), len(b)]
        report["files"][name] = entry
    report["other_differences"] = sum(f.get("differences_by_class", {}).get("other", 0) for f in report["files"].values())
    report["pickles_bytes_identical"] = all(v.get("bytes_identical", False) for k, v in report["files"].items() if k.endswith(".pkl"))
    report["missing_files"] = [k for k, v in report["files"].items() if "status" in v]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fresh", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    report = compare(args.fresh, args.reference)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    for name, f in report["files"].items():
        print(f"{name:34s} identical={f.get('bytes_identical', f.get('status'))} diffs={f.get('differences_by_class', '')}")
    print(json.dumps({k: report[k] for k in ("other_differences", "pickles_bytes_identical", "missing_files")}))
    return 0 if report["other_differences"] == 0 and not report["missing_files"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
