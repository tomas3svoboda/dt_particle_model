"""Faner et al. (2019) figure 2, soybean, 120 C: the measured vapour (gas) temperature as the Faner
march's bath, behind kernel-envelope FANER_MEASURED_VAPOUR_TRACE (off by default).

The owner's ruling, reported with
docs/GT_PS2_PAPER1_FANER_MEASURED_VAPOUR_TRACE_AND_DT_BATH_CONSTANTS_2026-10-28.md: Faner's measured
vapour temperature trace is the base-case bath of the Faner march, not a sensitivity. Built for
landing by the three-knob build of 2026-10-28/29 (record
docs/GT_PS2_PAPER1_THREE_HARNESS_KNOBS_BUILD_2026-10-29.md) from that record's measurement variant
(its records/variant_trace_A.diff, module faner_measured_trace.py), the same points, rule and
arithmetic.

The trace: the 17 digitized experimental gas_temperature points of SOURCE (figure 2, soybean,
series experimental, superseded rows excluded), 117.87 C at 0 s to 117.82 C at 235.058 s, maximum
120.706 C at 54.981 s; numpy.interp in time between the points, constant before the first (0 s)
and after the last (235.058 s). boundary_at(base, t) replaces only the bulk temperature of the
boundary and of its film boundary (the rule of faner_postcore_measured_bath.py); no other
coefficient changes. The file is read on first use, never at import, so a run with the knob absent
does not read it. physically_qualifying false.
"""

import csv
import hashlib
from dataclasses import replace
from pathlib import Path

import numpy as np

import particle_kernel_envelope as envelope

SOURCE_RELATIVE = "paper/analysis/datasets/faner2019_v2/curves.csv"
SOURCE = Path(__file__).resolve().parents[1] / "datasets" / "faner2019_v2" / "curves.csv"
POINT_COUNT = 17
LAST_CLOCK_S = 235.058
RULE = (
    "linear interpolation in time between the 17 measured points; constant before the first point "
    "(0 s) and after the last (235.058 s); only the bulk temperature of the boundary and its film "
    "boundary replaced"
)
SELECTION = {
    "figure": "2",
    "material": "soybean",
    "series": "experimental",
    "quantity": "gas_temperature",
    "superseded_by": "",
}
_TRACE = {}


class MeasuredTraceRefused(ValueError):
    """The declared trace file does not give the declared 17 points from 0 s."""


def _load():
    if _TRACE:
        return _TRACE
    data = SOURCE.read_bytes()
    with SOURCE.open(encoding="utf-8") as stream:
        rows = [
            r
            for r in csv.DictReader(stream)
            if r["figure"] == "2"
            and r["material"] == "soybean"
            and r["series"] == "experimental"
            and r["quantity"] == "gas_temperature"
            and not r["superseded_by"]
        ]
    rows.sort(key=lambda r: float(r["time_s"]))
    clocks = tuple(float(r["time_s"]) for r in rows)
    temperatures = tuple(float(r["value"]) + 273.15 for r in rows)
    if len(clocks) != POINT_COUNT or clocks[0] != 0.0 or clocks[-1] != LAST_CLOCK_S:
        raise MeasuredTraceRefused(
            f"the measured trace must be the {POINT_COUNT} points from 0 s to {LAST_CLOCK_S} s"
        )
    _TRACE.update(
        clocks_s=clocks,
        temperatures_k=temperatures,
        point_ids=tuple(r["point_id"] for r in rows),
        sha256=hashlib.sha256(data).hexdigest(),
        sha256_lf=hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest(),
    )
    return _TRACE


def armed(environ=None):
    """True when DTDC_KERNEL_ENVELOPE_FANER_MEASURED_VAPOUR_TRACE is exactly "1"."""
    return envelope.armed(envelope.FANER_MEASURED_VAPOUR_TRACE, environ)


def temperature_at(clock_s):
    """The measured gas temperature (K) at clock_s, by the declared rule."""
    trace = _load()
    return float(np.interp(clock_s, trace["clocks_s"], trace["temperatures_k"]))


def boundary_at(base, clock_s):
    """base with the bulk temperature of the boundary and its film boundary at the trace."""
    t = temperature_at(clock_s)
    return replace(
        base, temperature_k=t, film_boundary=replace(base.film_boundary, temperature_k=t)
    )


def bath_boundary(base, clock_s):
    """The Faner bath at clock_s: the trace when the knob is armed, else base itself."""
    return boundary_at(base, clock_s) if armed() else base


def declared_row():
    """The trace's stamp: the file by path, its sha256, the points and the rule."""
    trace = _load()
    return {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "knob": envelope.PREFIX + envelope.FANER_MEASURED_VAPOUR_TRACE,
        "record": envelope.THREE_KNOBS_RECORD,
        "measurement": envelope.FANER_TRACE_MEASUREMENT,
        "physically_qualifying": False,
        "source": SOURCE_RELATIVE,
        "source_sha256": trace["sha256"],
        "source_sha256_lf": trace["sha256_lf"],
        "selection": dict(SELECTION),
        "point_ids": list(trace["point_ids"]),
        "clocks_s": list(trace["clocks_s"]),
        "temperatures_k": list(trace["temperatures_k"]),
        "rule": RULE,
        "replaced": "the bulk temperature of the boundary and of its film boundary only",
    }
