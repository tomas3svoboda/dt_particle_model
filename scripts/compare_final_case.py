"""Compare a fresh qualification case directory with its archived original.

Every JSON and JSONL record that the archive holds is compared value by value (floats
exactly, no tolerance).  Differences are classified so that only the expected provenance
fields are excused:

  wall_clock   wall-clock durations and timestamps (keys naming wall, elapsed, perf, utc, ...)
  path         strings that name the producer's or the rerun's own directories
  git_stamp    the era stamp of the checkout (repo_head, repo_branch)
  source_hash  source-digest maps (the shipped tree is a subset of the producing tree)
  knob_stamp   the stamp blocks of the declared driver options (KNOB_STAMPS below); an archived case
               produced on a byte mirror with the option as a mirror variant carries no stamp or a
               different one, a rerun from this tree carries the option's own stamp
  knob_registry  the kernel-envelope registry of declared options in force (its 'armed' list and
               'declared_rows' texts): descriptive rows, not values the solver computes; they differ where
               the archived case ran with an option as a mirror variant and the rerun with the landed option
  source_hash  also the count of native source pins checked (the public tree ships the import closure of src,
               a subset of the producing tree)
  other        everything else -- a reproduction defect unless explained

The independent whole-history audit is produced by a separate auditor, not by the runner; a rerun
does not write it, so its absence in the fresh directory is listed, not a failure.

The archive carries no pickled checkpoints, run logs or temporary files; such files in the
fresh directory are listed as "not_archived" and are not a failure.  A record the archive
holds that the fresh run did not write IS a failure.

    python scripts/compare_final_case.py <fresh case dir> <archived case dir> [--json out.json]
"""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

WALL = re.compile(r"(wall|elapsed|perf_counter|_utc|utc_|timestamp|started|finished|seconds_since|clock_s|cpu_s|runtime_s)", re.I)
HASHKEY = re.compile(r"(sha256|source_hash|hashes|source_pins)", re.I)
GITSTAMP = re.compile(r"^(repo_(head|branch|dirty|status)|worktree_dirty)$")
PATHLIKE = re.compile(r"([A-Za-z]:[\/]|[\/](paper|reruns|scratchpad)[\/]|\.pkl$|\.json$|\.jsonl$)")
NOT_ARCHIVED = (".pkl", ".tmp", ".log")
SEPARATELY_PRODUCED = ("independent_whole_history_audit.json",)
KNOB_STAMPS = ("faner_measured_vapour_trace", "front_stage_noise_aware_energy_gate", "kernel_envelope_birth_seed_bound",
               "declared_dry_advance_ceiling")


def classify(path, a, b):
    joined = "/".join(map(str, path))
    if any(isinstance(k, str) and k in KNOB_STAMPS for k in path):
        return "knob_stamp"
    if (len(path) >= 2 and path[0] in ("kernel_envelope", "kernel_envelope_hygiene")
            and path[1] in ("armed", "declared_rows")):
        return "knob_registry"
    if path and path[-1] == "native_source_pins_checked":
        return "source_hash"
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


def load_jsonl(p):
    return [json.loads(line) for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]


def load_json(p):
    return json.loads(p.read_text(encoding="utf-8"))


def compare(fresh, ref):
    report = {"fresh": str(fresh), "reference": str(ref), "files": {}, "not_archived": []}
    fresh_names = {p.name for p in fresh.iterdir() if p.is_file()}
    ref_names = {p.name for p in ref.iterdir() if p.is_file()}
    report["not_archived"] = sorted(n for n in fresh_names - ref_names if n.endswith(NOT_ARCHIVED))
    for name in sorted((fresh_names | ref_names) - set(report["not_archived"])):
        fa, fb = fresh / name, ref / name
        if not fa.is_file() or not fb.is_file():
            report["files"][name] = {"status": "missing_in_fresh" if not fa.is_file() else "missing_in_reference"}
            continue
        entry = {"bytes_identical": fa.read_bytes() == fb.read_bytes(),
                 "sha256_fresh": hashlib.sha256(fa.read_bytes()).hexdigest(),
                 "sha256_reference": hashlib.sha256(fb.read_bytes()).hexdigest()}
        if name.endswith((".json", ".jsonl")):
            loader = load_jsonl if name.endswith(".jsonl") else load_json
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
    report["missing_files"] = [k for k, v in report["files"].items() if "status" in v and k not in SEPARATELY_PRODUCED]
    report["separately_produced_not_compared"] = [k for k in SEPARATELY_PRODUCED if k in report["files"]]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fresh", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    report = compare(args.fresh, args.reference)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    for name, f in report["files"].items():
        print(f"{name:34s} identical={f.get('bytes_identical', f.get('status'))} diffs={f.get('differences_by_class', '')}")
    print(json.dumps({k: report[k] for k in ("other_differences", "missing_files", "not_archived")}))
    return 0 if report["other_differences"] == 0 and not report["missing_files"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
