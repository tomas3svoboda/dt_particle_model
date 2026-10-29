"""Existing audited Faner history reader; functions reused without changing their definitions.
Original reader SHA-256: 796bca7c3e910e96731b08624c122c52b27e2eea4e14af1ea02335ac611096e0.
Only path configuration and the main-program wrapper were replaced for distribution.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
MAIN = Path(__file__).resolve().parents[1]
RADIUS = 0.0015
HOLDER_CLOCK = 235.058
MEAS_FALL = 25.350441949123166
OUTLET = 240.0

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_case(folder):
    case = json.loads((folder / "case.json").read_text(encoding="utf-8"))
    birth = json.loads((folder / "birth.json").read_text(encoding="utf-8"))
    rows, betas, birth_betas = [], [], []
    for r in birth["feed_rows"]:
        rows.append({"t": r["t_s"], "H": r["hexane_kg_kg_dry"], "W": r["retained_water_kg_kg_dry"],
                     "T": r["temperature_k"], "front_z": 1.0, "stage": "A"})
        if r.get("beta") is not None:
            betas.append(abs(r["beta"]))
    feed_end = rows[-1]["t"]
    core_end = None
    flux = []
    area = 4 * math.pi * RADIUS ** 2
    counts = {}
    for name, stage in (("front.history.jsonl", "B"), ("dry.history.jsonl", "C")):
        with (folder / name).open(encoding="utf-8") as stream:
            for line in stream:
                r = json.loads(line)
                counts[f"{stage}:{r['kind']}"] = counts.get(f"{stage}:{r['kind']}", 0) + 1
                if r["kind"] == "initial" and name == "front.history.jsonl":
                    s0 = r.get("surface") or {}
                    if s0.get("beta") is not None:
                        birth_betas.append(abs(s0["beta"]))
                if r["kind"] in ("refusal", "initial") or r["time_s"] <= rows[-1]["t"]:
                    continue
                surface = r.get("surface") or {}
                if surface.get("beta") is not None:
                    betas.append(abs(surface["beta"]))
                rows.append({"t": r["time_s"], "H": r["hexane_kg_kg_dry"], "W": r["particle_water_kg_kg_dry"],
                             "T": r["temperature_k"],
                             "front_z": r["front_z"] if r.get("front_z") is not None else (0.0 if stage == "C" else None),
                             "stage": stage})
                if r.get("front_z") is not None:
                    core_end = r["time_s"]
                if stage == "C":
                    ledger = r["native_report"]["ledger"]
                    total = ledger["stefan_recovery"]["total_stefan_fluxes_mol_m2_s"][-1]
                    upstream = r["y_hexane_cells"][-1] if total >= 0 else surface["surface_y_hexane"]
                    net = ledger["boundary_hexane_out_mol"] / (area * ledger["dt_s"])
                    motion = upstream * total
                    flux.append({"t": r["time_s"], "diffusion_percent": 100 * (net - motion) / net})
    t = [r["t"] for r in rows]
    assert all(b > a for a, b in zip(t, t[1:])), folder
    assert abs(t[-1] - OUTLET) < 1e-9, (folder, t[-1])
    return {"case": case, "birth_inputs": birth["inputs"], "rows": rows, "feed_end": feed_end,
            "core_end": core_end, "flux": flux, "betas": betas, "birth_betas": birth_betas,
            "row_kinds": counts}


def arr(c, key):
    return np.array([np.nan if r[key] is None else r[key] for r in c["rows"]], dtype=float)


def at(c, key, clock):
    return float(np.interp(clock, arr(c, "t"), arr(c, key)))


def crossing(c, level):
    t, x = arr(c, "t"), arr(c, "H")
    for i in range(len(t) - 1):
        if x[i] >= level >= x[i + 1]:
            return float(t[i] + (t[i + 1] - t[i]) * (x[i] - level) / (x[i] - x[i + 1]))
    raise ValueError("crossing missing")


def measured(quantity):
    base = MAIN / "paper/analysis/datasets/faner2019_v2"
    with (base / "adjudication_per_point.csv").open(newline="", encoding="utf-8") as stream:
        bands = {r["point_id"]: float(r["uncertainty_abs"]) for r in csv.DictReader(stream)}
    with (base / "curves.csv").open(newline="", encoding="utf-8") as stream:
        return sorted(
            (float(r["time_s"]), float(r["value"]), bands[r["point_id"]])
            for r in csv.DictReader(stream)
            if r["figure"] == "2" and r["material"] == "soybean" and r["temperature_c"] == "120"
            and r["quantity"] == quantity and r["series"] == "experimental" and not r["superseded_by"]
        )


def source_scores(c):
    out = {}
    for quantity, key, offset in (("solvent_loading", "H", 0.0), ("particle_temperature", "T", 273.15)):
        pts = [p for p in measured(quantity) if p[0] <= OUTLET]
        dev = np.array([at(c, key, tm) - offset - v for tm, v, _ in pts])
        unc = np.array([u for _, _, u in pts])
        out[quantity] = {
            "points": len(pts),
            "rms_deviation": float(np.sqrt(np.mean(dev ** 2))),
            "mean_deviation": float(np.mean(dev)),
            "max_abs_deviation": float(np.max(np.abs(dev))),
            "max_abs_at_clock_s": float(pts[int(np.argmax(np.abs(dev)))][0]),
            "within_uncertainty": int(np.sum(np.abs(dev) <= unc)),
            "deviations": [float(d) for d in dev],
            "clocks_s": [p[0] for p in pts],
        }
    return out


def observables(c):
    c20, c10 = crossing(c, 0.20), crossing(c, 0.10)
    W, H = at(c, "W", OUTLET), at(c, "H", OUTLET)
    diff = [f["diffusion_percent"] for f in c["flux"]]
    return {
        "FanerFreeEnd_s": c["feed_end"],
        "FanerDryTime_s": c["core_end"],
        "FanerCrossTwenty_s": c20,
        "FanerCrossTen_s": c10,
        "FanerFallingInterval_s": c10 - c20,
        "FanerFallingShortfall_pct": 100 * (1 - (c10 - c20) / MEAS_FALL),
        "FanerExitHexane_kg_kg_dry": H,
        "FanerExitWater_kg_kg_dry": W,
        "FanerExitTemperature_C": at(c, "T", OUTLET) - 273.15,
        "FanerLastHolderClockTemperature_C": at(c, "T", HOLDER_CLOCK) - 273.15,
        "outlet_moisture_pct_wb": 100 * W / (1 + W + H),
        "outlet_hexane_ppm_wb": 1e6 * H / (1 + W + H),
        "FanerMaxBeta": max(c["betas"] + c["birth_betas"]),
        "max_abs_beta_after_birth": max(c["betas"]),
        "regime_C_advances": len(c["flux"]),
        "advances_with_diffusion_opposing": sum(d < 0 for d in diff),
        "diffusion_percent_range": [min(diff), max(diff)] if diff else None,
    }


def three_rung(v_c, v_m, v_f, ratio=2.0):
    """Unchanged from the final ladder's numbers_from_cases.py."""
    e1, e2 = v_m - v_c, v_f - v_m
    res = {"values": [v_c, v_m, v_f], "e_coarse": e1, "e_fine": e2}
    if e1 == 0 or e2 == 0 or (e1 > 0) != (e2 > 0):
        res.update({"order": None, "reason": "differences change sign (or vanish)", "envelope": max(abs(e1), abs(e2))})
        return res
    p = math.log(e1 / e2) / math.log(ratio)
    res["order"] = p
    if p <= 0:
        res.update({"reason": "same sign but not decreasing (p <= 0)", "envelope": max(abs(e1), abs(e2))})
        return res
    corr = e2 / (ratio ** p - 1)
    res.update({"richardson": v_f + corr, "band": abs(corr), "gci_fine": 1.25 * abs(corr),
                "band_relative": abs(corr) / abs(v_f), "gci_fine_relative": 1.25 * abs(corr) / abs(v_f)})
    if p < 0.5:
        res.update({"reason": "order below 0.5: differences barely decrease, Richardson value not credible",
                    "envelope": max(abs(e1), abs(e2))})
    return res


def history_array(c):
    W, H, T = arr(c, "W"), arr(c, "H"), arr(c, "T")
    return np.column_stack([arr(c, "t"), 1e6 * H / (1 + W + H), 100 * W / (1 + W + H), T, H])


def whole_history(a, b):
    """a, b: history arrays; a minus b on the union of clocks; hexane relative with b's outlet floor."""
    clocks = np.union1d(a[:, 0], b[:, 0])
    out = {}
    for i, label in ((1, "hexane_ppm_wb"), (2, "moisture_points"), (3, "temperature_k"), (4, "hexane_kg_kg_dry")):
        d = np.interp(clocks, a[:, 0], a[:, i]) - np.interp(clocks, b[:, 0], b[:, i])
        k = int(np.argmax(np.abs(d)))
        out[label] = {"maximum_absolute_difference": float(abs(d[k])), "signed": float(d[k]),
                      "at_clock_s": float(clocks[k])}
    floor = abs(b[-1, 1])
    ah, bh = np.interp(clocks, a[:, 0], a[:, 1]), np.interp(clocks, b[:, 0], b[:, 1])
    rel = 100 * np.abs(ah - bh) / np.maximum(np.abs(bh), floor)
    k = int(np.argmax(rel))
    out["hexane_relative_outlet_floor_pct"] = {"maximum": float(rel[k]), "at_clock_s": float(clocks[k]),
                                                "floor_ppm": float(floor)}
    return out
