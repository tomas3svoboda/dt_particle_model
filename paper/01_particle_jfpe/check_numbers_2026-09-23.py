#!/usr/bin/env python
"""Mechanical numbers audit of the particle manuscript.

Part B of the referee audit (2026-09-23), extended on 2026-09-24 so that it
scores the repaired manuscript: every assertion now carries an ``expect``
field, "record" for a number that must be found in its named record and
"absent" for a superseded literal that must no longer appear anywhere in the
live manuscript sources, and the run ends with an explicit count of
mismatch / stale / unrecorded.

Extended again later on 2026-09-24, with ``N15``--``N30``, for the numbers the
fold-in of the re-characterized fatal limits put on the page: where the
measured diffusivity law is evaluated (limit 2), the placement and the reach
of the two domain extensions (limit 3), and the resolution bound measured
against the ladder's departure-tangent refusals (limit 1, whose record lives
under ``docs/evidence/`` and not under ``paper/analysis/``).

Extended once more on 2026-09-24, with ``Q1``--``Q17``, for the D22 mesh-flux
fold-in: the observed orders of the boundary-consistent reconstruction of the
arbitrary Lagrangian--Eulerian mesh flux against the shipped upwinded one, its
ledgers and its identity at zero mesh speed (``item_20_mms_mesh_flux``), and
the measurement that no trajectory reported in the manuscript marches the
operator it repairs (``item_21_mesh_flux_trajectory_delta``).

Extended on 2026-09-25 twice: with ``O1``--``O18`` and ``X15``--``X22`` for the
oil arm's source-pinned activity coefficient (``item_23_oil_arm_activity``),
and with ``J1``--``J28`` for the journal-condition sensitivity sweep folded
into Sec. S10.9 --- the thesis-run band it is scored against, the swept box,
the two inert axes, the per-axis spans in band widths, the per-trace table and
the three cross-checks (``item_22_faner_journal_sensitivity``).

Extended on 2026-09-26 for the independent-review pass, and only where a
number moved location: ``B1``--``B7`` for the dry-meal reading of the printed
sample mass (0.51-0.61 duty, 1.45-1.75 crossing) now printed beside the
charge-reading headline, the 95 per cent prediction factor 3.17 now printed in
Sec. 5.2, and the space-resolved window residuals 0.770 and 0.697 now printed
in Sec. 5.2; ``N11`` re-keyed to the 30.62 that Sec. S5.1 prints after the
corrected voidage Reynolds number left the main text; and an optional
``page`` key in the coverage pass for a literal a float cannot carry (the
trailing zero of 0.770).

Extended again on 2026-09-26 for the measurement fold-in of the independent
review's priorities 1 and 3, with ``K1``--``K30`` against the two refinement
records (``item_24_faner_refinement_and_band``, the outer-cell read refined
and the 95 per cent band marched; ``item_27_boundary_value_surface_read``, the
law read at the surface boundary value, whose numbers the paper now prints)
and ``Z1``--``Z7`` for the sentences that measurement withdrew. The printed
9.9 and 20.6 per cent stay under ``F9``/``F10``: the paper still prints them,
named as the superseded reported-resolution values.

Extended a third time on 2026-09-26 for the uptake fold-in (placeholder P5),
with ``W1``--``W33`` against ``item_26_uptake_refinement`` (the refined
factor and retardation, the refined converted count and band, the orders,
the largest change of simulated uptake, the identity count, the ledgers, the
withheld sub-step and the stated departure, the Table S27/S28 cells and the
stopped isotherm branch) and ``Z8``--``Z9`` for the placeholder and the
Sec. 5 sentence that said the uptake comparison was not refined.

Extended a fourth time on 2026-09-26 for the fourth referee read
(``REFEREE_AUDIT_4_2026-09-26.md``), with ``A1``--``A28`` against
``item_28_prediction_factor_boundary_read`` (the 95 per cent factor re-derived
at the boundary read, 3.140 and 3.132, and the band re-marched, at most 5.7 to
47.4 and 10.4 to 61.2 per cent), ``item_26`` (the figure-3.17-free 7.6 to 18.6,
18.4 refined, and the refined values now parenthesised), the last-tray total
floors and the surface-edge residuals beside the space-resolved ones; ``B3`` and
``B4`` re-keyed to the record row (1.74, not the 1.75 of the prose);
``Z10``--``Z16`` for the superseded sentences.

Extended a fifth time on 2026-09-27 for the owner's ruling "Adopt" on the
particle geometry of the superheated-hexane comparison
(``docs/GT_PS2_RULING_FANER_COMPARISON_GEOMETRY_2026-09-27.md``), with
``G1``--``G50`` against ``item_30_aris_sphere_ladder_and_population`` (the
declared volume-to-surface sphere: the ladder, the continuum values, the
orders, the 95 per cent band, the population, the conventions, the
adversarial Fourier number), ``H1``--``H20`` against
``item_29_faner_geometry_and_temperature_sensitivity`` (the radius sweep, the
D over R squared identity, the temperature alternatives), ``N18`` re-keyed to
the 11.53 that Sec. S15 prints after the main-text 7.0-to-11.5 sentence was
replaced by a pointer, and ``Z17``--``Z19`` for the withdrawn sentences.

Extended a sixth time on 2026-09-27 for the fifth referee read
(``REFEREE_AUDIT_5_2026-09-27.md``), with ``I1``--``I68`` against
``item_31_species_own_geometry`` (each thesis sample's own volume-to-surface
sphere, now the declared geometry: the radii, the ladder, the continuum values,
the orders, the 95 per cent band, the temperature alternatives and the D over R
squared identity at those radii, the population scaled to each sample, the
adversarial sphericity), ``G27`` re-keyed to the thesis diameter 1.95 mm that
Sec. 5.2 now prints in place of the 0.975 mm radius, and ``Z20``--``Z30`` for
the withdrawn wording ("a derivation, not a fit", "the source particle's", the
sunflower-inside-the-band claim at the declared sphere, the population's
"up when the mass-mean is", "halves the soybean miss", "calibrated", the bold
declared row, the one-sphericity sentence, "either way", the 2 per cent
diameter sentence and "at a fixed radius").

Extended a seventh time on 2026-09-27 for the sixth referee read
(``REFEREE_AUDIT_6_2026-09-27.md``), with ``L1``--``L13`` (the borrowed
sphericity's effect from item_29's radius sweep and item_31's adversarial
paragraph, the thesis volume-equivalent soybean residual from item_30's
conventions table, the faster temperature edge's signed residual from
item_31's temperature table, both 95 per cent factors) and ``Z31``--``Z38``
for the withdrawn wording ("anywhere inside", "because its mass-mean is",
"a 1.0 mm class", the unlabelled abstract comparator, the unqualified
"both directions improve", the single factor 3.14, the misattached caveat and
the singular sphere).

Extended an eighth time on 2026-09-26 (the machine date; the day's earlier
extensions are stamped 2026-09-27, one day ahead) for the space-resolved solve
marched at the declared spheres and refined
(``item_32_space_resolved_declared_geometry``), with ``S1``--``S40`` (the
declared-sphere residuals and refusals, the 120-cell refusals of Table S40,
the stride-ceiling extension, the refusal state, the ledgers, the identity)
and ``Z39``--``Z41`` for the withdrawn wording ("the one reading no evaluation
choice affects", "lies farther from both", "has not been marched on
Whitaker's").

Extended a ninth time on 2026-09-27 for the owner's strategic clarity pass
(``STRATEGIC_PASS_2026-09-27.md``), with ``T1``--``T6`` against
``item_33_refinement_ladder_d20_tree`` (the ladder re-run on the current tree,
printed in S14.3 and pointed to from Sec. 6.1) and ``P1``--``P12`` for the
wording the pass withdrew (the defect-log framing of the abstract, the
unnamed "source" phrases, the disagreement labels of Table 2, the question
headings). Numbers moved from the main text to the supplement keep their
assertions, which search both.

Extended a tenth time on 2026-09-27 for the corrections fold-in
(``CORRECTIONS_FOLDIN_2026-09-27.md``): version 2 of the three digitizations.
``REKEY_2026_09_27`` re-keys in place every assertion whose printed number
moved, to the version-2 value and the re-scoring record that carries it
(``item_35_faner2008_v2_rescoring``, ``item_36_cardarelli_v2_rescoring``,
``item_37_faner2019_v2_rescoring``); ``WITHDRAWN_2026_09_27`` turns the
vessel-gas floor (U8, U10, O10; ``item_34_published_vapour_floors``) into
expected absences; ``C1``--``C45`` bind the numbers the fold-in put on the
page (the floors at published compositions, the oil-arm intervals, the
uptake re-fit, the mechanism paragraph from ``item_38_steam_implications_note``)
and ``Z42``--``Z58`` the withdrawn wording.

Extended an eleventh time on 2026-09-27 for the final polish
(``READER_CHECK_RESPONSE_2026-09-27.md``): ``REKEY_POLISH_2026_09_27`` re-keys
F20, N7 and B4 to the journal conditions on version 2
(``item_39_journal_ranges_v2_and_outcomes``); ``P1``--``P20`` bind the
numbers item_39 put on the page (the crossing and dry-crossing ranges, the
ratio to the demand, the Table S7 journal columns, S13.2, S18, the item_29
S2 re-evaluation, the bracketing criterion's 3.7 and 4.2 per cent) and
``Z59``--``Z68`` the withdrawn ranges and wording.

Extended a twelfth time on 2026-09-27 for the condensation of the main text
(``CONDENSATION_PLAN_2026-09-27.md``): no number changed and no assertion was
added or removed.  ``RELOCATED_CONDENSATION_2026_09_27`` re-labels the
``where`` of the twelve record assertions whose number left the main text
and is now printed only in the Supplementary Material, in the verbatim moved
paragraphs; the coverage pass, which searches both documents, checks them
there.  ``F13`` gains the ``page`` literal 0.410, the regression's printed
value, which stays in Sec. 5.2 (its float token 0.41 had been matching an
unrelated activity that the condensation moved to S9.11).

Extended a thirteenth time on 2026-09-28 for the condensation of the
Supplementary Material (``CONDENSATION_SUPPLEMENT_PLAN_2026-09-28.md``): the
119-page supplement became the condensed ``supplementary.tex`` and, unchanged
in content, the extended technical supplement ``supplementary_extended.tex``.
No number changed and no assertion was added, removed or weakened.  The
absence checks now run on both supplements (``EXT_SOURCES``).  The coverage
pass first searches the submitted documents (the main text and
``supplementary.tex``); an assertion found only in the extended report counts
as reaching the page only when it is named in
``RELOCATED_SUPPLEMENT_CONDENSATION_2026_09_28`` with the extended-report
section that prints it, and it is reported there and counted separately.

Extended a fourteenth time on 2026-09-30 for the paper-1 rewrite R-P1
(``REFEREE_AUDIT_7_2026-09-30.md``): ``Y1``--``Y61`` bind the refinement
ladder re-earned on the F36 farm campaign (two labelled tiers, the chains,
the Richardson orders and grid-convergence indices, the declared bound's
rejecting outcomes, the hosts); ``E1``--``E11`` the face-departure mechanism
and the ledger reach (the build record, M-F1, the PART-01 and PART-02
records); ``M1``--``M35`` the comparator C1 of the surface read against the
weighted mean of Crank and the verified audit items; ``Z69``--``Z84`` the
withdrawn wording.  C41 is withdrawn from the page (audit A10).  Four
evidence folders had not landed when the block was written: an assertion on
them reports AWAITING-LANDED-EVIDENCE, counted apart, unless the environment
variable named in ``PENDING_ROOTS`` points at the staging folder.

Extended a fifteenth time on 2026-09-30 for the eighth referee read, part 1
(``REFEREE_AUDIT_8_2026-09-30.md``): ``TR1``--``TR86`` the two-regime
derivation of Sec. 5.3, Table 4 and S9.15; ``RM1``--``RM61`` the regime map
as the Fickian shell's, its pressure-feasible form, the resistance table and
the constant-gas march at the band state; ``FE1``--``FE7`` the coupled march's
admissible limit and cost; ``HA1``--``HA24`` the held window at three meshes;
``FC1``--``FC7`` the film axis's lower corner and the F37b axis; ``TC1``--``TC9``
and ``SG1``--``SG3`` the 96-cell time chain and the pending items 04 and 10;
``Z90``--``Z108`` the withdrawn wording.

Extended a sixteenth time on 2026-09-30 for the eighth referee read, part 2
(``REFEREE_AUDIT_8_2026-09-30.md``, section 8.M): ``MF1``--``MF159`` the
moving-front march of the journal runs of Faner et al. (2019) printed in the
abstract, Secs. 5, 5.4, 7.2, 7.3, the conclusions, Tables 2 and 5 and S9.16
(the emergency march records, step 1 and step 2, and their sidecars; step 2
staged through ``PENDING_ROOTS`` until it lands); ``MZ1``--``MZ10`` the
withdrawn wording that said the moving front was not marched.

Extended a seventeenth time on 2026-09-30 for the eighth referee read, part 3
(``REFEREE_AUDIT_8_2026-09-30.md``, section 8.N): ``EN1``--``EN*`` the soybean
size-ensemble march of Sec. 5.4 (Table 5, Figure 3, the abstract, Secs. 5,
7.2, 7.3, the conclusions), the candidates of S9.16 (Table S13) and the
ensemble in full with the sunflower run of S9.17 (Tables S15 and S16), bound
to the step-3 and step-4 records and sidecars, both placed in the tree when
the block was run (``PENDING_ROOTS`` entries kept, as for the blocks above);
``EZ1``--``EZ8`` the
withdrawn wording (the three unmeasured candidates, the unexplained factor of
two, the old titles).  The two records join the corpus as files.

Extended an eighteenth time on 2026-09-30 for the eighth referee read, part 4
(``REFEREE_AUDIT_8_2026-09-30.md``, section 8.O, the owner's instruction that
the paper carry no certification or other internal process vocabulary):
``PZ1``--``PZ*`` assert that vocabulary ABSENT from the printed text of the
main text and both supplements.  They search with ``printed_only``: LaTeX
comments and the arguments of the url, label and ref commands are removed
first, because a repository path in the data statement or a label key is not
printed prose.  No numeric assertion changed.

Extended a nineteenth time on 2026-10-01 for the eighth referee read, part 6
(``REFEREE_AUDIT_8_2026-09-30.md``, section 8.Q, the whole-manuscript referee
read): ``P6_1``--``P6_*`` bind the numbers the pass prints (the sample clock of
the soybean ensemble beside the crossing-matched 0.85, the plateau end by one
rule on both curves, the rate exponent across the classes, the dry-shell limit
on the declared spheres) to the landed step-4 superposition, the landed
two-regime derivation, and the part-6 derived record (staged through
``PENDING_ROOTS`` as ``DTDC_PAPER_PART6_DERIVED`` until it lands);
``P6Z1``--``P6Z*`` assert the withdrawn wording absent.  An absence check may
now carry ``sources``, the manuscript files it searches (default: all three
documents), so that wording withdrawn from the submitted text is not asserted
absent from the deposited extended report.

Extended a twentieth time on 2026-10-01 for the eighth referee read, part 7
(``REFEREE_AUDIT_8_2026-09-30.md``, section 8.R): ``C14`` is re-keyed in place
to the version-2 re-score of Table 3's sunflower cell of the "each 2019 meal's
diameter" row (11.6 per cent, the record
``GT_PS2_PAPER1_TABLE3_SUNFLOWER_RESCORE_2026-10-01``, its evidence staged
through ``PENDING_ROOTS`` as ``DTDC_PAPER_TABLE3_RESCORE`` until it lands);
``P7_1``--``P7_*`` bind the companions (13 of 14, Table S21's 0.116, the
version-1 9.1 and 0.091 kept as a pointer, the radius 0.670 mm);
``P7Z1``--``P7Z*`` assert the withdrawn wording absent, the abstract's
unqualified clock, and every literal supplement pointer ("Sec.~S17.5") in the
main text, which now refers to the supplement by label through xr.

Three passes:

  Pass 1 (coverage)   every numeric literal printed in the manuscript sources
                      that is a *result* (filtered by an exclusion list of
                      structural, LaTeX and declared-parameter numbers) is
                      looked up in the curated record corpus at its own printed
                      precision.  A number found nowhere is reported as
                      UNRECORDED-BY-SEARCH, which is a pointer for a human
                      check, not a verdict on its own.

  Pass 2 (assertions) a curated table of load-bearing numbers, each with the
                      record file and the field that should carry it, checked
                      to the printed digit.  This is the pass whose verdicts
                      enter the report.

  Pass 3 (internal)   every number printed in more than one place in the
                      manuscript is checked for agreement between those places.

Read-only.  Writes nothing.  Run from the repository root:

  $env:PYTHONHASHSEED='1'; $env:OMP_NUM_THREADS='1'
  $env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'
  $env:NUMEXPR_NUM_THREADS='1'; $env:PYTHONDONTWRITEBYTECODE='1'
  $env:PYTHONIOENCODING='utf-8'
  .\.venv\Scripts\python.exe paper\01_particle_jfpe\check_numbers_2026-09-23.py
"""

from __future__ import annotations

import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
ANALYSIS = PAPER / "analysis"

MAIN_SOURCES = [
    "main.tex",
     "sec_intro_reshape.tex",
    "sec_formulation.tex",
     "sec_discretization_reshape.tex",
     "sec_verification_reshape.tex",
     "sec_validation_reshape.tex",
     "sec_results_reshape.tex",
     "sec_discussion_reshape.tex",
]
ARCHIVED_SUPP_SOURCES = ["supplementary.tex"] + [
     "moved_full/" + name for name in (
          "sec_intro.tex", "sec_discretization.tex", "sec_verification.tex",
          "sec_validation.tex", "sec_results.tex", "sec_discussion.tex",
     )
]
#: 2026-09-28: the full, pre-condensation supplement, deposited and not submitted.
SUPP_SOURCES = [
     "supplementary_submission.tex", "supplementary_submission_detail.tex",
     "supplementary_submission_protocol.tex", "supplementary_submission_evidence.tex",
     "supplementary_submission_dtmarch.tex",
]
EXT_SOURCES = ["supplementary_extended.tex", *ARCHIVED_SUPP_SOURCES]

# ---------------------------------------------------------------------------
# 0.  Source reading
# ---------------------------------------------------------------------------

#: lines that carry no result number at all
SKIP_LINE = re.compile(
    r"\\(documentclass|usepackage|input|newcommand|renewcommand|setcounter"
    r"|setlength|includegraphics|label|bibliography|captionsetup|begin\{tabular\}"
    r"|cmidrule|toprule|midrule|bottomrule|addlinespace|vspace|hspace)"
)

#: numeric literals that belong to the typesetting, not to the physics
STRUCTURAL = re.compile(
    r"(0\.\d+\\textwidth"            # column widths
    r"|table-format=[^\]]*"          # siunitx table formats
    r"|\[\d+pt\]"                    # skips
    r"|\d+pt\b"                      # lengths
    r"|S\[[^\]]*\]"                  # siunitx column specs
    r"|width=[\d.]+"
    r"|\\the[a-z]+)"
)

NUMBER = re.compile(
    r"(?<![A-Za-z0-9_.])"
    r"(\d+\.\d+(?:[eE][-+]?\d+)?|\d+(?:[eE][-+]?\d+)?)"
    r"(?![0-9])"
)


@dataclass
class Printed:
    value: float
    text: str
    file: str
    line: int
    context: str


def strip_structural(line: str) -> str:
    return STRUCTURAL.sub(" ", line)


def read_printed(files: list[str]) -> list[Printed]:
    out: list[Printed] = []
    for name in files:
        path = HERE / name
        if not path.exists():
            continue
        for i, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            if SKIP_LINE.search(line):
                continue
            clean = strip_structural(line)
            for m in NUMBER.finditer(clean):
                tok = m.group(1)
                try:
                    val = float(tok)
                except ValueError:
                    continue
                out.append(
                    Printed(val, tok, name, i, clean[max(0, m.start() - 60) : m.end() + 60])
                )
    return out


# ---------------------------------------------------------------------------
# 1.  The curated record corpus
# ---------------------------------------------------------------------------

#: raw per-step solver dumps are excluded: they carry millions of floats and
#: would make every printed number "match" something by accident.
CORPUS_EXCLUDE_DIRS = {"cells", "logs", "runs", "out", "__pycache__"}
CORPUS_SUFFIXES = {".md", ".csv", ".json", ".txt"}
CORPUS_MAX_BYTES = 3_000_000


def corpus_files() -> list[Path]:
    files: list[Path] = []
    # 2026-09-30 (R-P1): the rewrite's numbers live in evidence outside
    # paper/analysis; each extra root is read when it has landed or is staged.
    for rec in globals().get("EXTRA_CORPUS_ROOTS", []):
        root, _ = resolve_record(rec)
        if root is None or not root.exists():
            continue
        # 2026-09-30 (REFEREE_AUDIT_8 part 3): a file root (a record) is one file
        if root.is_file():
            if root.suffix.lower() in CORPUS_SUFFIXES and root.stat().st_size <= CORPUS_MAX_BYTES:
                files.append(root)
            continue
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in CORPUS_SUFFIXES:
                continue
            if EXTRA_CORPUS_EXCLUDE & set(p.name for p in path.parents):
                continue
            if path.stat().st_size > CORPUS_MAX_BYTES:
                continue
            files.append(path)
    for path in ANALYSIS.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in CORPUS_SUFFIXES:
            continue
        if CORPUS_EXCLUDE_DIRS & set(p.name for p in path.parents):
            continue
        if path.stat().st_size > CORPUS_MAX_BYTES:
            continue
        files.append(path)
    return sorted(files)


@dataclass
class Corpus:
    files: list[Path] = field(default_factory=list)
    #: rounded-string -> set of file indexes
    index: dict[str, set[int]] = field(default_factory=dict)
    values: list[tuple[float, int]] = field(default_factory=list)


def build_corpus() -> Corpus:
    c = Corpus(files=corpus_files())
    for fi, path in enumerate(c.files):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for m in NUMBER.finditer(text):
            try:
                v = float(m.group(1))
            except ValueError:
                continue
            c.values.append((v, fi))
    return c


def significant_digits(tok: str) -> int:
    t = tok.lower().split("e")[0].lstrip("-+0.")
    return len(t.replace(".", "").rstrip()) or 1


def matches_at_precision(printed: float, tok: str, candidate: float) -> bool:
    """True when `candidate` rounds to `printed` at the precision printed."""
    if printed == 0:
        return candidate == 0
    sig = significant_digits(tok)
    if candidate == 0:
        return False
    if printed * candidate < 0:
        return False
    try:
        exp = math.floor(math.log10(abs(printed)))
    except ValueError:
        return False
    half = 0.5 * 10 ** (exp - sig + 1)
    return abs(abs(candidate) - abs(printed)) <= half * 1.000001


def lookup(c: Corpus, p: Printed, limit: int = 4) -> list[str]:
    hits: list[str] = []
    seen: set[int] = set()
    for v, fi in c.values:
        if fi in seen:
            continue
        if matches_at_precision(p.value, p.text, v):
            seen.add(fi)
            # 2026-09-30 (R-P1): extra corpus roots may be staged outside paper/
            try:
                hits.append(str(c.files[fi].relative_to(PAPER)))
            except ValueError:
                hits.append(str(c.files[fi]))
            if len(hits) >= limit:
                break
    return hits


# ---------------------------------------------------------------------------
# 2.  The curated assertion table
# ---------------------------------------------------------------------------
# Each entry: (id, printed value, unit, where printed, record path, how to find
# it in that record, expected record value or None when the check is a search).

#: 2026-09-27: the two records behind the geometry ruling.
_I29 = "analysis/results_2026-09-19/item_29_faner_geometry_and_temperature_sensitivity/RESULT.md"
_I30D = "analysis/results_2026-09-19/item_30_aris_sphere_ladder_and_population"
_I30 = _I30D + "/RESULT.md"
#: 2026-09-27, fifth referee read: each sample's own volume-to-surface sphere.
_I31 = "analysis/results_2026-09-19/item_31_species_own_geometry/RESULT.md"
#: 2026-09-26 (machine date): the space-resolved solve at the declared spheres.
_I32D = "analysis/results_2026-09-19/item_32_space_resolved_declared_geometry"
_I32 = _I32D + "/RESULT.md"
#: 2026-09-27, strategic pass: the refinement ladder re-run on the current tree.
_I33 = "analysis/results_2026-09-19/item_33_refinement_ladder_d20_tree/RESULT.md"

ASSERTIONS: list[dict] = [
    # --- superheated-hexane cluster -----------------------------------------
    dict(id="F1", printed=68.7145, unit="degC", where="sec_validation.tex constant-rate temperature",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"68\.7145"),
    dict(id="F2", printed=0.137, unit="degC", where="sec_validation.tex / tab:basis",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.137"),
    dict(id="F3", printed=0.093, unit="degC", where="sec_validation.tex / tab:basis",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.093"),
    dict(id="F4", printed=0.19935, unit="kg/kg", where="sec_validation.tex critical loading",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.19935"),
    dict(id="F5", printed=1.257, unit="-", where="sec_validation.tex GAB activity at Xc",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"1\.257"),
    dict(id="F6", printed=45.9, unit="-", where="sec_validation.tex fixed Re_eps",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig", grep=r"45\.9"),
    dict(id="F7", printed=131.3, unit="W/m2/K", where="sec_validation.tex borrowed film",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"131\.3"),
    dict(id="F8", printed=103.7, unit="W/m2/K", where="sec_formulation.tex / sec_validation.tex",
         record="analysis", grep=r"103\.7"),
    dict(id="F9", printed=9.9, unit="percent", where="tab:fanermarch sunflower window RMS",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.09[89]|9\.9"),
    dict(id="F10", printed=20.6, unit="percent", where="tab:fanermarch soybean window RMS",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.20[56]|20\.6"),
    dict(id="F11", printed=1.165e-9, unit="m2/s", where="sec_validation.tex identified scalar",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"1\.165"),
    dict(id="F12", printed=2.91, unit="-", where="sec_validation.tex identified/fixed",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"2\.91"),
    dict(id="F13", printed=0.410, unit="ln D", where="sec_validation.tex regression sigma",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.410"),
    dict(id="F14", printed=1.507, unit="-", where="sec_validation.tex sigma as a factor",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"1\.50[67]"),
    dict(id="F15", printed=3.513e-9, unit="m2/s", where="eq:defflaw prefactor",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"3\.513"),
    dict(id="F16", printed=0.959, unit="-", where="eq:defflaw beta",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.959"),
    # --- uptake / tail ------------------------------------------------------
    # 2026-09-25 (third referee read, M5): the benchmark retardation and the
    # surface-state capacity are the 50 C set of item_13's 34 simulated rows,
    # a different set from the 34 tabulated states of item_05 that carry
    # 1.19-18.98.  Both are now printed in Sec. S10.5 as well as in 5.2.
    dict(id="U1", printed=8.53, unit="-",
         where="sec_validation.tex / S10.5 benchmark retardation, low end",
         record="analysis/results_2026-09-19/item_13_cardarelli_uptake_benchmark", grep=r"8\.5[23]"),
    dict(id="U2", printed=19.54, unit="-",
         where="sec_validation.tex / S10.5 benchmark retardation, high end",
         record="analysis/results_2026-09-19/item_13_cardarelli_uptake_benchmark", grep=r"19\.5"),
    dict(id="U2a", printed=8.22, unit="-",
         where="sec_validation.tex / S10.5 surface-state capacity, low end",
         record="analysis/results_2026-09-19/item_13_cardarelli_uptake_benchmark", grep=r"8\.22"),
    dict(id="U2b", printed=37.36, unit="-",
         where="sec_validation.tex / S10.5 surface-state capacity, high end",
         record="analysis/results_2026-09-19/item_13_cardarelli_uptake_benchmark", grep=r"37\.36"),
    dict(id="U2c", printed=18.98, unit="-",
         where="sec_validation.tex / S10.4 tabulated-state retardation, high end",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21",
         grep=r"18\.98"),
    # U3 removed 2026-09-24 (second referee read, Pass 2b): the bare "1216"
    # this once tracked is gone by design, superseded by the two floors the
    # repair separated at every occurrence, 1216.9 ppm (U7) and 1216.358 /
    # 1216.4 ppm (U8) --- printing an unqualified "1216" again would undo
    # that separation.
    dict(id="U4", printed=0.9031, unit="-", where="sec_validation.tex structural factor",
         record="analysis/results_2026-09-19/item_06_cardarelli_tail", grep=r"0\.903"),
    # the measured final loading is a digitized datum, not an item_06 output
    dict(id="U5", printed=7.415e-3, unit="kg/kg", where="sec_validation.tex measured final loading",
         record="analysis/datasets/faner2019/curves.csv", grep=r"0\.007415"),
    dict(id="U6", printed=0.040, unit="-", where="sec_validation.tex form-mapping RMS low end",
         record="../docs/evidence/form_mapping_2026-09-15", grep=r"0\.040"),
    # 2026-09-25 (third referee read, M2): the pure-vapour floor left the
    # abstract, which now carries the vessel-gas floor alone; 5.3 and the
    # conclusions still print it with its construction.
    dict(id="U7", printed=1216.88, unit="ppm",
         where="sec_validation.tex / sec_discussion.tex pure-vapour floor",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21",
         grep=r"1216\.88"),
    dict(id="U8", printed=1216.358, unit="ppm", where="sec_validation.tex vessel-composition floor",
         record="analysis/results_2026-09-19/item_06_cardarelli_tail", grep=r"1216\.3"),
    # --- verification / numerics --------------------------------------------
    dict(id="V1", printed=1.8e-12, unit="-", where="tab:verification energy ledger",
         record="analysis/results_2026-09-19/item_16_mms_coupled_thermal", grep=r"1\.78[0-9]*e-12|1\.787"),
    dict(id="V2", printed=0.998, unit="order", where="tab:verification temporal order",
         record="analysis/results_2026-09-19/item_16_mms_coupled_thermal", grep=r"0\.998"),
    dict(id="V3", printed=3.4228875e-10, unit="-", where="sec_verification.tex datum invariance",
         record="analysis/results_2026-09-19/item_15_datum_invariance_cold_cache", grep=r"3\.42288"),
    dict(id="V4", printed=0.72, unit="order", where="tab:verification moving-front spatial order",
         record="analysis/results_2026-09-19/item_16_mms_coupled_thermal", grep=r"0\.71[89]|0\.72"),
    dict(id="V5", printed=0.75, unit="order", where="tab:verification moving-front spatial order",
         record="analysis/results_2026-09-19/item_16_mms_coupled_thermal", grep=r"0\.74[67]|0\.75"),
    # --- results ------------------------------------------------------------
    # the ladder observables are recorded in item_01, not item_10
    dict(id="R1", printed=8.161, unit="percent", where="sec_results.tex worst ladder observable",
         record="analysis/results_2026-09-19/item_01_refinement_ladder", grep=r"8\.161"),
    dict(id="R2", printed=0.4529, unit="K", where="sec_results.tex worst temperature spread",
         record="analysis/results_2026-09-19/item_01_refinement_ladder", grep=r"0\.4529"),
    # the sunflower particle-area demand: the paper prints the 2026-09-19 value
    dict(id="R9", printed=15.98, unit="W/m2/K", where="supplementary.tex S13.2, sunflower",
         record="analysis/results_2026-09-19/item_05_faner_comparison/item05_constant_rate.csv",
         grep=r"15\.97"),
    dict(id="R10", printed=17.228, unit="W/m2/K", where="the live value S13.2 should carry",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig",
         grep=r"17\.22"),
    # printed value updated 2026-09-24 (Pass 2b): the length pass tightened
    # the page from 2.81e-11 to the record's own 2.809e-11; same record, same
    # digits, now checked at the precision actually printed.
    dict(id="R3", printed=2.809e-11, unit="-", where="sec_results.tex departure-tangent residual",
         record="analysis/results_2026-09-19/item_10_dynamic_stress", grep=r"2\.80[89]"),
    dict(id="R4", printed=414, unit="points", where="sec_results.tex regime grid",
         record="analysis/results_2026-09-19/item_07_regime_map", grep=r"414"),
    dict(id="R5", printed=375, unit="points", where="sec_results.tex receding points",
         record="analysis/results_2026-09-19/item_07_regime_map", grep=r"375"),
    dict(id="R6", printed=29.8, unit="percent", where="sec_results.tex boundary energy",
         record="analysis/results_2026-09-19/item_08_coefficient_sensitivity", grep=r"29\.8|0\.298"),
    dict(id="R7", printed=1.41, unit="K", where="sec_results.tex peak temperature",
         record="analysis/results_2026-09-19/item_08_coefficient_sensitivity", grep=r"1\.41"),
    # NB: a literal "344.6" is NOT in the record; the record carries 344.5609751681902.
    # A literal regex missed it on the first run. This is the tool's own failure mode:
    # a printed number rounded from a record value needs a rounding-aware pattern.
    dict(id="R8", printed=344.6, unit="-", where="sec_formulation.tex Xc/Xe at the top of the band",
         record="analysis/results_2026-09-19/item_07_regime_map", grep=r"344\.5[6-9]|344\.6"),
    # the ladder job counts: the live generation is 36 jobs / 5 complete
    dict(id="R12", printed=36, word=True, unit="jobs", where="sec_results.tex / supplementary.tex S14.3",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/rerun_2026-09-20",
         grep=r"\"jobs_seen\":\s*36"),
    dict(id="R14", printed=5, word=True, unit="jobs", where="sec_results.tex completed count",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/rerun_2026-09-20",
         grep=r"\"jobs_completed\":\s*5"),
    dict(id="R15", printed=30.26, unit="percent", where="sec_results.tex second pair, worst relative",
         record="analysis/results_2026-09-19/item_01_refinement_ladder",
         grep=r"30\.26|0\.302647"),
    dict(id="R16", printed=0.9362, unit="K", where="sec_results.tex second pair, worst temperature",
         record="analysis/results_2026-09-19/item_01_refinement_ladder",
         grep=r"0\.9362|0\.936151"),
    dict(id="R17", printed=8.09e-4, unit="-", where="sec_results.tex / S14.1 mesh triple",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/RESULT.md",
         grep=r"8\.09e-4"),
    dict(id="R18", printed=2.23e-3, unit="-", where="sec_results.tex / S14.1 wide-pair triple",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/RESULT.md",
         grep=r"2\.23e-3"),
    dict(id="R19", printed=8.13e-11, unit="-", where="sec_results.tex ladder tangent residual, low",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/cells",
         grep=r"8\.133e-11", raw=True),
    dict(id="R20", printed=5.03e-2, unit="-", where="sec_results.tex ladder tangent residual, high",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/cells",
         grep=r"5\.026e-02", raw=True),
    dict(id="R21", printed=25, word=True, unit="gates", where="tab:verification acceptance criteria",
         record="analysis/results_2026-09-19/item_10_dynamic_stress",
         grep=r"twenty-five boolean acceptance"),
    # --- the superseded literals that must no longer be printed -------------
    dict(id="X1", printed=27, unit="jobs", where="the 2026-09-19 ladder job count",
         expect="absent", grep=r"twenty-seven jobs"),
    dict(id="X2", printed=23, unit="jobs", where="the 2026-09-19 refusal count",
         expect="absent", grep=r"other twenty-three"),
    dict(id="X3", printed=47, unit="tests", where="'47 acceptance tests'",
         expect="absent", grep=r"47 acceptance tests"),
    dict(id="X4", printed=3.819, unit="percent", where="the superseded finite-pressure seam",
         expect="absent", grep=r"\$3\.819\\,\\%\$ discrepancy"),
    dict(id="X5", printed=7.9e-4, unit="-", where="the superseded mesh triple",
         expect="absent", grep=r"\$7\.9\$, \$9\.6\$ and"),
    dict(id="X6", printed=15.98, unit="W/m2/K", where="the superseded S13.2 demand",
         expect="absent", grep=r"\\SI\{15\.98\}"),
    dict(id="X7", printed=1706.8, unit="W/kg", where="the superseded S13.2 duty",
         expect="absent", grep=r"1706\.8"),
    dict(id="X8", printed=3.2184, unit="m2/kg", where="the superseded S13.2 second area",
         expect="absent", grep=r"3\.2184"),
    dict(id="X9", printed=25, unit="-", where="the hard-coded 'eq. (25) of the main text'",
         expect="absent", grep=r"eq\.~\(25\)"),
    dict(id="X10", printed=8, unit="-", where="the main text's rounded 8-to-19 ratio",
         expect="absent", grep=r"coefficient \$8\$\s*\n?to \$19\$"),
    dict(id="X11", printed=1216, unit="ppm", where="the unqualified 1216 ppm label",
         expect="absent", grep=r"\\SI\{1216\}\{\\ppm\}|\\SI\{2655\}\{\\ppm\}"),
    dict(id="X12", printed=1.49, unit="bar", where="the superseded ceiling arithmetic",
         expect="absent", grep=r"\\SI\{1\.49\}\{\\bar\}"),
    dict(id="X13", printed=16, unit="W/m2/K", where="the superseded 'about 16' particle demand",
         expect="absent", grep=r"demands about \\SI\{16\}"),
    dict(id="X14", printed=0, unit="-", where="the -ise spelling",
         expect="absent", grep=r"memoised"),
    # --- 2026-09-25: the ideal-solution oil arm the source pin superseded ---
    dict(id="X15", printed=1.863e-3, unit="kg/kg", where="the ideal-solution floor at 120 C",
         expect="absent", grep=r"1\.863\\times10\^\{-3\}"),
    dict(id="X16", printed=3.883e-3, unit="kg/kg", where="the ideal-solution floor at 100 C",
         expect="absent", grep=r"3\.883\\times10\^\{-3\}"),
    dict(id="X17", printed=3.98, unit="-", where="the ideal-solution measured-over-floor pair",
         expect="absent", grep=r"\$3\.98\$ and \$1\.91\$"),
    dict(id="X18", printed=23.6, unit="percent", where="the ideal-solution gap closure, supplement",
         expect="absent", grep=r"\$23\.6\$ and \\SI\{37\.0\}"),
    dict(id="X19", printed=24, unit="percent", where="the ideal-solution gap closure, main text",
         expect="absent", grep=r"closing \$24\$ to \\SI\{37\}"),
    dict(id="X20", printed=17, unit="percent", where="the ideal-solution double-count exposure",
         expect="absent", grep=r"\$17\$ to \\SI\{53\}"),
    dict(id="X21", printed=1.5, unit="-", where="the declared ideal-to-1.5 coefficient interval",
         expect="absent", grep=r"bracketed over \$1\.0\$ to \$1\.5\$"),
    dict(id="X22", printed=0, unit="-", where="the withdrawn 'none is held' oil-arm sentence",
         expect="absent", grep=r"coefficient that would close it nor a fatty-acid"),
    # --- derived trend / envelope ------------------------------------------
    dict(id="D1", printed=105.840, unit="degC", where="tab:basis / sec_validation.tex discharge",
         record="analysis/results_2026-09-19/item_12_published_envelope", grep=r"105\.84"),
    dict(id="D2", printed=108.351, unit="degC", where="sec_validation.tex M11 pairing",
         record="analysis/results_2026-09-19/item_18_water_branch_placement", grep=r"108\.35"),
    dict(id="D3", printed=4.82, unit="K", where="sec_validation.tex isotherm spread",
         record="analysis/results_2026-09-19/item_18_water_branch_placement", grep=r"4\.82"),
    dict(id="D4", printed=105.915, unit="degC", where="sec_validation.tex marched discharge",
         record="analysis/results_2026-09-19/item_17_derived_trend_temperature", grep=r"105\.91"),
    # D5 removed 2026-09-24 (second referee read, Pass 2b): the length pass
    # cut the sentence that stated this offset as its own number. It is not
    # lost -- Table~S10 (T13) still prints the two values it comes from,
    # 136.00000 and 135.97613 degC, on the page.
    dict(id="D6", printed=18.240, unit="percent", where="sec_validation.tex band entry moisture",
         record="analysis/results_2026-09-19/item_12_published_envelope", grep=r"18\.24"),
    dict(id="D7", printed=41.4, unit="percent", where="sec_validation.tex band admissibility",
         record="analysis/results_2026-09-19/item_12_published_envelope", grep=r"41\.4|0\.414"),
    dict(id="D8", printed=61.4, unit="percent", where="sec_validation.tex band admissibility",
         record="analysis/results_2026-09-19/item_12_published_envelope", grep=r"61\.4|0\.614"),
    dict(id="D9", printed=1.24, unit="kJ/mol", where="sec_formulation.tex isosteric heat low end",
         record="analysis/results_2026-09-19/item_11_adsorption_isotherm", grep=r"1\.24"),
    dict(id="D10", printed=4.39, unit="kJ/mol", where="sec_formulation.tex isosteric heat high end",
         record="analysis/results_2026-09-19/item_11_adsorption_isotherm", grep=r"4\.39"),
    dict(id="D11", printed=10.3, unit="K", where="sec_discussion.tex isotherm gap",
         record="analysis/results_2026-09-19/item_11_adsorption_isotherm", grep=r"10\.3"),
    dict(id="D12", printed=367.863, unit="K", where="sec_discussion.tex warmest water measurement",
         record="analysis/results_2026-09-19/item_11_adsorption_isotherm", grep=r"367\.863"),
    dict(id="D13", printed=1.2376, unit="kJ/mol", where="sec_formulation.tex source column, low end",
         record="analysis/results_2026-09-19/item_11_adsorption_isotherm", grep=r"1\.2376"),
    dict(id="D14", printed=4.3885, unit="kJ/mol", where="sec_formulation.tex source column, high end",
         record="analysis/results_2026-09-19/item_11_adsorption_isotherm", grep=r"4\.3885"),
    dict(id="D15", printed=90.27, unit="-", where="sec_formulation.tex Xc/Xe at the band top",
         record="analysis/results_2026-09-19/item_07_regime_map", grep=r"90\.274"),
    dict(id="D16", printed=105.832, unit="degC", where="sec_validation.tex warmest placed entry",
         record="analysis/results_2026-09-19/item_18_water_branch_placement", grep=r"105\.83"),
    # --- the two floors, separated -----------------------------------------
    dict(id="U9", printed=1216.9, unit="ppm", where="abstract / 5.3: pure-vapour total floor",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21",
         grep=r"1216\.88"),
    dict(id="U10", printed=1216.4, unit="ppm", where="tab:basis / 5.3 / 5.5 / 7.4: vessel-gas floor",
         record="analysis/results_2026-09-19/item_06_cardarelli_tail", grep=r"1216\.3"),
    dict(id="U11", printed=2654.6, unit="ppm", where="supplementary.tex S9.9 pure-vapour floor",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21",
         grep=r"2654\.60|0\.0026546"),
    dict(id="U12", printed=148793, unit="Pa", where="supplementary.tex S9.6 stage-entry pressure",
         record="analysis/results_2026-09-19/item_06_cardarelli_tail", grep=r"148 ?793|148793"),
    # --- the form-mapping set ----------------------------------------------
    # printed value updated 2026-09-24 (Pass 2b): the record's 7.62 is
    # declared on the page as 7.6, so the assertion now tracks what is
    # actually printed; the record grep is unchanged and still matches.
    dict(id="U13", printed=7.6, unit="-", where="sec_validation.tex ratio low end",
         record="../docs/evidence/form_mapping_2026-09-15", grep=r"7\.62"),
    dict(id="U14", printed=18.9, unit="-", where="sec_validation.tex ratio high end",
         record="../docs/evidence/form_mapping_2026-09-15", grep=r"18\.9"),
    dict(id="U15", printed=18.6, unit="-", where="sec_validation.tex ratio without figure 3.17",
         record="../docs/evidence/form_mapping_2026-09-15", grep=r"18\.61"),
    # --- the constant-rate closure -----------------------------------------
    dict(id="F17", printed=16.3, unit="percent", where="sec_validation.tex earlier soybean window",
         record="analysis/results_2026-09-19/item_05_faner_comparison/RESULT_2026-09-22.md",
         grep=r"16\.3|0\.163"),
    dict(id="F18", printed=9.7, unit="percent", where="sec_validation.tex best sunflower window",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig",
         grep=r"9\.7 % with 15 of 15"),
    dict(id="F19", printed=18.65, unit="W/m2/K", where="supplementary.tex S13.2 demand, high",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig",
         grep=r"18\.652"),
    dict(id="F20", printed=1671.4, unit="W/kg", where="supplementary.tex S13.2 duty, low",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig",
         grep=r"1671\.43"),
    dict(id="F21", printed=2.9231, unit="m2/kg", where="supplementary.tex S13.2 single area",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig",
         grep=r"2\.92314"),
    dict(id="F22", printed=6.18, unit="-", where="supplementary.tex S13.1 soybean third set",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"6\.18"),
    # --- section 5.4 scope corrections --------------------------------------
    dict(id="D17", printed=0.29, unit="percent", where="sec_validation.tex plateau miss, low",
         record="analysis/results_2026-09-19/item_14_derived_trend_benchmark", grep=r"0\.29"),
    dict(id="D18", printed=0.78, unit="percent", where="sec_validation.tex plateau miss, high",
         record="analysis/results_2026-09-19/item_14_derived_trend_benchmark", grep=r"0\.78"),
    dict(id="D19", printed=6.43, unit="K", where="sec_validation.tex elevation at that activity",
         record="analysis/results_2026-09-19/item_18_water_branch_placement", grep=r"6\.430"),
    # D20/D21 re-pinned 2026-09-25: the oil arm now runs at a source-pinned
    # activity coefficient (item_23), so the closure moved from 23.6/37.0 to
    # 33.2/44.7 and its record moved with it.
    dict(id="D20", printed=33.2, unit="percent", where="supplementary.tex S9.11 oil-arm gap closure, low",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity", grep=r"33\.2 %"),
    dict(id="D21", printed=44.7, unit="percent", where="supplementary.tex S9.11 oil-arm gap closure, high",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity", grep=r"44\.7 %"),
    # --- 2026-09-25: the oil arm's source-pinned activity coefficient -------
    dict(id="O1", printed=5.2, unit="-", where="supplementary.tex S9.11 pinned coefficient",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"\*\*5\.2\*\*, bracket"),
    dict(id="O2", printed=4.6, unit="-", where="supplementary.tex S9.11 coefficient bracket, low",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"bracket \*\*4\.6 to 5\.9\*\*"),
    dict(id="O3", printed=5.9, unit="-", where="supplementary.tex S9.11 coefficient bracket, high",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"4\.6 to 5\.9"),
    dict(id="O4", printed=22.8, unit="percent", where="S9.11 / sec_discussion four-method spread",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"22\.8"),
    dict(id="O5", printed=1.46, unit="MPa", where="supplementary.tex S9.11 Henry constant at 105 C",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"1\.46 MPa"),
    dict(id="O6", printed=2.216e-3, unit="kg/kg", where="S9.11 / sec_validation floor with the arm, 120 C",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"2\.216"),
    dict(id="O7", printed=4.201e-3, unit="kg/kg", where="S9.11 / sec_validation floor with the arm, 100 C",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"4\.201"),
    dict(id="O8", printed=3.35, unit="-", where="supplementary.tex S9.11 measured over floor, 120 C",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"\*\*3\.35\*\*"),
    dict(id="O9", printed=1.77, unit="-", where="supplementary.tex S9.11 measured over floor, 100 C",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"\*\*1\.77\*\*"),
    dict(id="O10", printed=2085.5, unit="ppm", where="supplementary.tex S9.11 vessel floor with the arm",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"2085\.5"),
    dict(id="O11", printed=8.7, unit="percent", where="S9.11 Henry cross-check against 14.2 atm",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"8\.7 % against 14\.2"),
    dict(id="O12", printed=13.5, unit="percent", where="S9.11 Henry cross-check against 15.5 atm",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"13\.5 % against 15\.5"),
    dict(id="O13", printed=100.9, unit="degC", where="S9.11 the first printed Henry constant",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"100\.9"),
    dict(id="O14", printed=102.2, unit="degC", where="S9.11 the second printed Henry constant",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"102\.2"),
    dict(id="O15", printed=0.67, unit="-", where="S9.11 pinned arm over the measured equilibrium",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"\*\*0\.67\*\*"),
    dict(id="O16", printed=2.89, unit="-", where="S9.11 mole-basis arm over the measured equilibrium",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"\*\*2\.89\*\*"),
    dict(id="O17", printed=1.5, unit="-", where="S9.11 lowest ratio of condition to Henry limit",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"1\.48"),
    dict(id="O18", printed=2.5, unit="-", where="S9.11 highest ratio of condition to Henry limit",
         record="analysis/results_2026-09-19/item_23_oil_arm_activity",
         grep=r"2\.53"),
    # --- 2026-09-25: the journal-condition sensitivity sweep (item_22), the
    #     numbers Sec. S10.9 and the two main-text sentences put on the page --
    dict(id="J1", printed=0.099450, unit="-", where="S10.9 width of the thesis-run band",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.099450"),
    dict(id="J2", printed=0.117872, unit="-", where="S10.9 thesis-run band, low end",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.117872"),
    dict(id="J3", printed=0.018422, unit="-", where="S10.9 thesis-run band, high end",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.018422"),
    dict(id="J4", printed=960, unit="points", where="S10.9 scored grid points",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"\"grid_points_n\": 960"),
    dict(id="J5", printed=240, unit="conditions", where="S10.9 distinct scored conditions",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"240 distinct"),
    dict(id="J6", printed=0.2388, unit="kg/kg", where="S10.9 window lower bound",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.2388"),
    dict(id="J7", printed=0.490628, unit="kg/kg", where="S10.9 soybean curve's first point",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.490628"),
    dict(id="J8", printed=0.499279, unit="kg/kg", where="S10.9 sunflower curve's first point",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.499279"),
    dict(id="J9", printed=11.2, unit="band widths", where="S10.9 bed-depth span, soybean",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"11\.2 and 9\.7"),
    dict(id="J10", printed=9.7, unit="band widths", where="S10.9 bed-depth span, sunflower",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"11\.2 and 9\.7"),
    dict(id="J11", printed=6.9, unit="band widths", where="S10.9 gas-velocity span, soybean",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"6\.9 and 6\.2"),
    dict(id="J12", printed=6.2, unit="band widths", where="S10.9 gas-velocity span, sunflower",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"6\.9 and 6\.2"),
    dict(id="J13", printed=1.7, unit="band widths", where="S10.9 initial-loading span, soybean",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"1\.7 and 0\.3"),
    dict(id="J14", printed=0.3, unit="band widths", where="S10.9 initial-loading span, sunflower",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"1\.7 and 0\.3"),
    dict(id="J15", printed=0.7551, unit="-", where="Table S19 soybean residual span, low",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.7551"),
    dict(id="J16", printed=1.5716, unit="-", where="Table S19 soybean residual span, high",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"1\.5716"),
    dict(id="J17", printed=0.7840, unit="-", where="Table S19 sunflower residual span, low",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"-0\.783955"),
    dict(id="J18", printed=0.117, unit="m/s",
         where="S10.9 thesis printed approach velocity, low",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.117 to 0\.129"),
    dict(id="J19", printed=0.0437, unit="-", where="Table S19 soybean best printed coordinate",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.043704"),
    dict(id="J20", printed=0.0040, unit="-", where="Table S19 sunflower best printed coordinate",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"-0\.004006"),
    dict(id="J21", printed=0.0939, unit="-", where="Table S19 marched soybean residual",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"-0\.09389039"),
    dict(id="J22", printed=0.1767, unit="-", where="Table S19 marched sunflower residual",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"-0\.1766786"),
    dict(id="J23", printed=0.1243, unit="m/s", where="S10.9 marched soybean on the velocity axis",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.1243"),
    dict(id="J24", printed=0.1176, unit="m/s", where="S10.9 marched sunflower on the velocity axis",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.1176"),
    dict(id="J25", printed=1.038, unit="-", where="S10.9 implied footprint over the printed holder, soybean",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"1\.038"),
    dict(id="J26", printed=0.965, unit="-", where="S10.9 implied footprint over the printed holder, sunflower",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.965"),
    dict(id="J27", printed=0.129, unit="m/s",
         where="S10.9 thesis printed approach velocity, high",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"0\.117 to 0\.129"),
    dict(id="J28", printed=593.2, unit="-", where="S10.9 voidage Reynolds number, high",
         record="analysis/results_2026-09-19/item_22_faner_journal_sensitivity",
         grep=r"593\.2"),
    # --- 2026-09-24 second referee read: the numbers the repair itself put on
    #     the page, which the 2026-09-23 table could not have covered ---------
    dict(id="N1", printed=9.43, unit="percent", where="tab:verification / 6.3 ablation difference, low",
         record="analysis/results_2026-09-19/item_09_ablation_comparison", grep=r"9\.43"),
    dict(id="N2", printed=10.93, unit="percent", where="tab:verification / 6.3 ablation difference, high",
         record="analysis/results_2026-09-19/item_09_ablation_comparison", grep=r"10\.93"),
    dict(id="N3", printed=25.71, unit="percent", where="sec_results.tex wet volume fraction, low",
         record="analysis/results_2026-09-19/item_09_ablation_comparison", grep=r"25\.71"),
    dict(id="N4", printed=29.33, unit="percent", where="sec_results.tex wet volume fraction, high",
         record="analysis/results_2026-09-19/item_09_ablation_comparison", grep=r"29\.33"),
    dict(id="N5", printed=3.799, unit="percent", where="sec_verification.tex / Table S4 finite-pressure seam",
         record="analysis/results_2026-09-19/item_09_ablation_comparison", grep=r"3\.799|0\.037992"),
    dict(id="N6", printed=6.34e-16, unit="-", where="tab:verification ideal-binary reduction",
         record="analysis/results_2026-09-19/item_09_ablation_comparison", grep=r"6\.3423"),
    dict(id="N7", printed=16.18, unit="W/m2/K", where="sec_validation.tex / S13.2 particle demand, low",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig",
         grep=r"16\.177"),
    dict(id="N8", printed=9.55e-4, unit="-", where="sec_results.tex / S14.1 mesh triple, middle",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/RESULT.md",
         grep=r"9\.55e-4"),
    dict(id="N9", printed=9.03e-4, unit="-", where="sec_results.tex / S14.1 mesh triple, last",
         record="analysis/results_2026-09-19/item_01_refinement_ladder/RESULT.md",
         grep=r"9\.03e-4"),
    dict(id="N10", printed=0.6941, unit="-", where="sec_formulation.tex source film prefactor",
         record="../docs/evidence/reeps_definition_2026-09-22", grep=r"0\.6941"),
    # 2026-09-26 (independent review 1, length pass): the corrected voidage
    # Reynolds number left the main text's parameter paragraph; it is printed
    # in Sec. S5.1 at the record's own four digits, 30.62.  Same record.
    dict(id="N11", printed=30.62, unit="-",
         where="supplementary.tex S5.1 corrected voidage Reynolds (main text until 2026-09-25)",
         record="../docs/evidence/reeps_definition_2026-09-22", grep=r"30\.6"),
    # --- 2026-09-26, independent review 1: numbers that moved location ----
    # The dry-meal reading of the printed sample mass, marched at Whitaker's
    # correlation: duty ratio 0.51-0.61 and crossing-time ratio 1.45-1.75
    # (record table "lit_whitaker on the dry-meal area"), now printed beside
    # the charge-reading headline in the abstract, Table 2, Sec. 5.1 and
    # Sec. S10.8.
    dict(id="B1", printed=0.51, unit="-", where="abstract / tab:basis / Sec. 5.1 / S10.8 dry-reading duty, low",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig/DETAIL.md",
         grep=r"lit_whitaker on the dry-meal area \| 0\.57 / 1\.55 \| 0\.51 / 1\.74"),
    dict(id="B2", printed=0.61, unit="-", where="abstract / tab:basis / Sec. 5.1 / S10.8 dry-reading duty, high",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig/DETAIL.md",
         grep=r"\| 0\.61 / 1\.45 \|"),
    # 2026-09-26, fourth referee read minor 1: B3/B4 re-keyed from the
    # RESULT_2026-09-22 prose (which carries the 1.75 slip of the
    # source_central row) to the record row itself, lit_whitaker on the
    # dry-meal area: crossing ratios 1.55, 1.74, 1.45, 1.68.
    dict(id="B3", printed=1.45, unit="-", where="tab:basis / Sec. 5.1 / S10.8 dry-reading crossing, low",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig/DETAIL.md",
         grep=r"\| 0\.61 / 1\.45 \|"),
    dict(id="B4", printed=1.74, unit="-", where="tab:basis / Sec. 5.1 / S10.8 dry-reading crossing, high",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-22_rig/DETAIL.md",
         grep=r"lit_whitaker on the dry-meal area \| 0\.57 / 1\.55 \| 0\.51 / 1\.74"),
    dict(id="B5", printed=3.17, unit="-", where="Sec. 5.2 / S10.4 / S15 95 per cent prediction factor at the surface state",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface/RESULT.md",
         grep=r"3\.17 to 3\.18 at the surface state"),
    dict(id="B6", printed=0.770, page="0.770", unit="-", where="Sec. 5.2 / S10.7 space-resolved window residual, soybean",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.770|0\.77[0-9]"),
    dict(id="B7", printed=0.697, unit="-", where="Sec. 5.2 / S10.7 space-resolved window residual, sunflower",
         record="analysis/results_2026-09-19/item_05_faner_comparison", grep=r"0\.697|0\.69[67]"),
    # --- 2026-09-26, measurement fold-in: item_24 and item_27 -----------
    dict(id="K1", printed=24.4, unit="percent",
         where="abstract / tab:basis / Sec. 5.2 / tab:fanermarch / S18 sunflower window RMS, boundary read, as printed",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"\*\*0\.244\*\*"),
    dict(id="K2", printed=36.6, unit="percent",
         where="abstract / tab:basis / Sec. 5.2 / tab:fanermarch / S18 soybean window RMS, boundary read, as printed",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"\*\*0\.366\*\*"),
    dict(id="K3", printed=24.0, unit="percent",
         where="Sec. 5.2 / S18 converged sunflower window RMS",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"converged values 24\.0 and 35\.8"),
    dict(id="K4", printed=35.8, unit="percent",
         where="Sec. 5.2 / S18 converged soybean window RMS",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"converged values 24\.0 and 35\.8"),
    dict(id="K5", printed=1.96, unit="order",
         where="Sec. 5.2 / S18 observed mesh order, low",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"mesh 1\.96 to 2\.03"),
    dict(id="K6", printed=2.03, unit="order",
         where="Sec. 5.2 / S18 observed mesh order, high",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"mesh 1\.96 to 2\.03"),
    dict(id="K7", printed=1.05, unit="order",
         where="Sec. 5.2 / S18 observed time order, low",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"time 1\.05 to 1\.16"),
    dict(id="K8", printed=1.16, unit="order",
         where="Sec. 5.2 / S18 observed time order, high",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"time 1\.05 to 1\.16"),
    dict(id="K9", printed=0.8, unit="points",
         where="Sec. 5.2 / S18 reported level above continuum (0.4 and 0.8)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"within 0\.4 and 0\.8 points"),
    dict(id="K10", printed=2, unit="count",
         where="Sec. 5.2 / tab:fanermarch / S18 band counts as printed (2 and 0 of 15)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"2 of 15 and 0 of 15 points inside the reading band as printed"),
    dict(id="K11", printed=13.2, unit="percent",
         where="Sec. 5.2 / tab:fanermarch / S10.4 / S18 one-sigma, sunflower, fast edge",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"33\.9 %, 13\.2 %"),
    dict(id="K12", printed=33.9, unit="percent",
         where="Sec. 5.2 / tab:fanermarch / S10.4 / S18 one-sigma, sunflower, slow edge",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"33\.9 %, 13\.2 %"),
    dict(id="K13", printed=24.2, unit="percent",
         where="Sec. 5.2 / tab:fanermarch / S10.4 / S18 one-sigma, soybean, fast edge",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"46\.6 %, 24\.2 %"),
    dict(id="K14", printed=46.6, unit="percent",
         where="Sec. 5.2 / tab:fanermarch / S10.4 / S18 one-sigma, soybean, slow edge",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"46\.6 %, 24\.2 %"),
    dict(id="K15", printed=47.6, unit="percent",
         where="tab:basis / Sec. 5.2 / S10.4 / S18 95 per cent range, sunflower (5.7-47.6)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"at most 5\.7 to 47\.6 %"),
    dict(id="K16", printed=61.4, unit="percent",
         where="tab:basis / Sec. 5.2 / S10.4 / S18 95 per cent range, soybean (10.7-61.4)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"at most 10\.7 to 61\.4 %"),
    dict(id="K17", printed=14.7, unit="percent",
         where="S10.4 / S18 95 per cent fast edge, sunflower (8 of 15)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"14\.7 % \(8/15\)"),
    dict(id="K18", printed=5.7, unit="percent",
         where="S10.4 / S18 lowest marched in band, sunflower (15 of 15 at F^0.75)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"5\.7 % at D x F\^0\.75 \(15/15\)"),
    dict(id="K19", printed=10.7, unit="percent",
         where="S10.4 / S18 95 per cent fast edge, soybean (12 of 15)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"10\.7 % \(12/15\)"),
    dict(id="K20", printed=0.42, unit="-",
         where="S10.4 95 per cent widths over the span (0.42 and 0.51)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"Widths at least 0\.42 and 0\.51"),
    dict(id="K21", printed=2.78, unit="percent",
         where="sec_formulation.tex threshold / S18 switch width at 60 cells",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"2\.78,\s*3\.26, 3\.36 and 3\.38 per cent"),
    dict(id="K22", printed=3.38, unit="percent",
         where="sec_formulation.tex threshold / S18 switch width at 480 cells",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"2\.78,\s*3\.26, 3\.36 and 3\.38 per cent"),
    dict(id="K23", printed=9.8, unit="percent",
         where="sec_formulation.tex threshold, superseded read (9.7 to 9.8)",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"9\.7\s+to 9\.8 per cent"),
    dict(id="K24", printed=54658, page="54\\,658", unit="checks",
         where="S18 identity with the committed read switched in",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/RESULT.md", grep=r"54,658 of 54,658"),
    dict(id="K25", printed=66.6, unit="percent",
         where="S9.10 table note, journal soybean, boundary read",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/DETAIL.md", grep=r"\| 0\.6663 \|"),
    dict(id="K26", printed=60.2, unit="percent",
         where="S9.10 table note, journal sunflower, boundary read",
         record="analysis/results_2026-09-19/item_27_boundary_value_surface_read/DETAIL.md", grep=r"\| 0\.6018 \|"),
    dict(id="K27", printed=1.51, unit="-",
         where="Sec. 5.2 / S15 / S17 outer-cell over face loading at the reported mesh",
         record="analysis/results_2026-09-19/item_24_faner_refinement_and_band/RESULT.md", grep=r"1\.51, 1\.28 and 1\.15"),
    dict(id="K28", printed=19.7, unit="percent",
         where="Sec. 5.2 / S17 outer-cell read at the finest level, sunflower",
         record="analysis/results_2026-09-19/item_24_faner_refinement_and_band/RESULT.md", grep=r"\*\*0\.197\*\*"),
    dict(id="K29", printed=30.8, unit="percent",
         where="Sec. 5.2 / S17 outer-cell read at the finest level, soybean",
         record="analysis/results_2026-09-19/item_24_faner_refinement_and_band/RESULT.md", grep=r"\*\*0\.308\*\*"),
    dict(id="K30", printed=490, unit="checks",
         where="S17 identity of the committed marches",
         record="analysis/results_2026-09-19/item_24_faner_refinement_and_band/RESULT.md", grep=r"490 of 490 checks exact"),
    # sentences the measurement withdrew
    dict(id="Z1", printed=0, unit="-", where="the placeholder sentence of the threshold paragraph (P1)",
         expect="absent", grep=r"pending the refinement measurement"),
    dict(id="Z2", printed=0, unit="-", where="the withdrawn claim that the fitted scalar does not improve on the law (S10.6)",
         expect="absent", grep=r"neither direction improves on the"),
    dict(id="Z3", printed=0, unit="-", where="the withdrawn conclusion that the fitted scalar does no better",
         expect="absent", grep=r"does no better on the other"),
    dict(id="Z4", printed=0, unit="-", where="the placeholder refinement sentence of Sec. 5 (P4)",
         expect="absent", grep=r"No marched trajectory\s+compared here"),
    dict(id="Z5", printed=0, unit="-", where="the four measurement placeholders P1-P4",
         expect="absent", grep=r"PLACEHOLDER P[1-4]"),
    dict(id="Z6", printed=0, unit="-", where="the withdrawn Sec. 5.2 claim",
         expect="absent", grep=r"describes the two thesis windows better than"),
    dict(id="Z7", printed=0, unit="-", where="the withdrawn S9.10 reading",
         expect="absent", grep=r"the formulation reproduces the thesis falling-rate leg"),
    # --- 2026-09-26, the uptake fold-in: item_26 (placeholder P5) -------
    #     The printed uptake numbers stay at the reported level; the
    #     refined values, orders, ledgers, the withheld run and the
    #     stopped isotherm branch are asserted against item_26.
    dict(id="W1", printed=7.6, unit="-",
         where="Sec. 5.3 / S9.5 / S19 / tab:basis / Sec. 8 refined form-mapping factor, low (7.6 to 18.7)",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 7\.6–18\.7"),
    dict(id="W2", printed=18.7, unit="-",
         where="Sec. 5.3 / S9.5 / S19 / tab:basis / Sec. 8 refined form-mapping factor, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 7\.6–18\.7"),
    dict(id="W3", printed=8.51, unit="-",
         where="Sec. 5.3 / S10.5 / S19 refined retardation, low",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 8\.51–19\.27"),
    dict(id="W4", printed=19.27, unit="-",
         where="Sec. 5.3 / S10.5 / S19 refined retardation, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 8\.51–19\.27"),
    dict(id="W5", printed=34, unit="-",
         where="Sec. 5.3 / S19 converted in-box count, refined (34 to 35)",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\*\*34 to 35\*\* \(R6\)"),
    dict(id="W6", printed=4.0, page="4.0", unit="-",
         where="Sec. 5.3 / S19 converted band, refined (3.9 to 4.0)",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\*\*3\.9 to 4\.0\*\*"),
    dict(id="W7", printed=0.0053, unit="-",
         where="Sec. 5.3 / Sec. 7.3 / S19 largest change of simulated uptake at a measured time",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"0\.0053 \(converted, fixed D\*\)"),
    dict(id="W8", printed=2.0, page="2.00", unit="-",
         where="Sec. 5.3 / S19 mesh order of the marched uptake times, low",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\| 2\.00–2\.01 \| 0\.99–1\.06 \|"),
    dict(id="W9", printed=2.01, unit="-",
         where="Sec. 5.3 / S19 mesh order of the marched uptake times, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\| 2\.00–2\.01 \| 0\.99–1\.06 \|"),
    dict(id="W10", printed=0.97, unit="-",
         where="Sec. 5.3 / S19 time order of the marched uptake times, low",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"0\.97–1\.05 \(34/34\)"),
    dict(id="W11", printed=1.06, unit="-",
         where="Sec. 5.3 / S19 time order of the marched uptake times, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\| 2\.00–2\.01 \| 0\.99–1\.06 \|"),
    dict(id="W12", printed=128514, page="128\,514", unit="-",
         where="S19 identity of the reported level",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\*\*128,514 checks, 0 failures\*\*"),
    dict(id="W13", printed=515, unit="-",
         where="S19 marches over seven levels",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"515 marches \(59 at the reported level"),
    dict(id="W14", printed=1.7e-10, unit="-",
         where="S19 largest per-step ledger",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"largest 1\.70e-10 of the equilibrium inventory"),
    dict(id="W15", printed=3.71e-09, unit="-",
         where="S19 cumulative ledger",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"cumulative at most 3\.71e-9"),
    dict(id="W16", printed=0.034, page="0.0340", unit="-",
         where="S19 the refused sub-stride",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"sub-stride 0\.0340 s"),
    dict(id="W17", printed=74, unit="-",
         where="S19 the m4k4 runs not withheld (stated departure)",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"not to the other 74 runs"),
    dict(id="W18", printed=0.0045, unit="-",
         where="S19 largest uptake change, declared coefficients",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"is 0\.0045 \(declared coefficients\)"),
    dict(id="W19", printed=0.0047, unit="-",
         where="S19 largest uptake change, converted pipeline",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"0\.0047 \(converted, pipeline\)"),
    dict(id="W20", printed=0.0031, unit="-",
         where="S19 largest master-curve change",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"move by at most 0\.0031"),
    dict(id="W21", printed=2.6, unit="-",
         where="S19 largest printed-ratio change, per cent",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"largest 2\.6 per\s+cent"),
    dict(id="W22", printed=18.654, unit="-",
         where="tab:item26 finest form-mapping factor, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"7\.594–18\.654"),
    dict(id="W23", printed=18.53, unit="-",
         where="tab:item26 formal-order estimate of the factor, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"7\.592–18\.53"),
    dict(id="W24", printed=18.425, unit="-",
         where="tab:item26 finest factor without figure 3.17",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"7\.594–18\.425"),
    dict(id="W25", printed=17.683, unit="-",
         where="tab:item26 finest soybean factor",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"7\.594–17\.683"),
    dict(id="W26", printed=19.271, unit="-",
         where="tab:item26 finest retardation, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"8\.513–19\.271"),
    dict(id="W27", printed=31.79, unit="-",
         where="tab:item26 fixed-coefficient band factor, finest",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"\| 31\.77 \| 31\.79 \| 31\.80 \|"),
    dict(id="W28", printed=3.954, unit="-",
         where="tab:item26 converted band, finest pipeline",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"3\.954 \(pipeline\); 3\.948 \(fixed\)"),
    dict(id="W29", printed=2.1832, unit="-",
         where="tab:item26 converted time ratio, finest fixed, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"0\.5529–2\.1832 \(fixed\)"),
    dict(id="W30", printed=7.173e-09, unit="-",
         where="tab:item26 finest converted coefficient, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"2\.608e-10–7\.173e-9"),
    dict(id="W31", printed=1.87, unit="-",
         where="tab:item26orders fitted-factor time chain, no asymptotic order",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"0\.11–1\.87 on 12, one not contracting"),
    dict(id="W32", printed=2.27, unit="-",
         where="tab:item26orders converted residual, diagonal, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"0\.99–2\.27"),
    dict(id="W33", printed=0.0, word=True, unit="-",
         where="S19 / Sec. 5.3 isotherm fit uncertainty not printed by the source (branch stopped)",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"no standard\s+error, confidence interval or fit statistic"),
    dict(id="Z8", printed=0, unit="-", where="the uptake placeholder P5",
         expect="absent", grep=r"PLACEHOLDER P5"),
    dict(id="Z9", printed=0, unit="-", where="the Sec. 5 sentence that the uptake comparison is not refined",
         expect="absent", grep=r"uptake comparison is carried at\s+the discretization"),
    dict(id="N12", printed=0.253, unit="-", where="S9.6 / S9.9 pure-vapour hexane activity (left Sec. 5.4 on 2026-09-26, uptake length pass)",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21",
         grep=r"0\.2534"),
    # The vessel-gas hexane activity 0.222 of Sec. 5.3 and Sec. S9.6 is not
    # itself in any record: it is now printed as its own construction,
    # 0.612 (the record's own mole fraction) times 0.3625 (the record's own
    # pressure-feasible activity ceiling at 378.15 K).  This asserts the
    # ceiling factor of that construction against the record that computes it.
    dict(id="N13", printed=0.3625, unit="-",
         where="sec_validation.tex / S9.6 pressure-feasible activity ceiling factor",
         record="analysis/results_2026-09-19/item_06_cardarelli_tail", grep=r"activity 0\.3625"),
    dict(id="N14", printed=373.0, unit="K", where="sec_formulation.tex n-hexane isotherm domain top",
         record="analysis/results_2026-09-19/item_07_regime_map",
         grep=r"isotherm_measured_range_K\"?:\s*\[\s*323\.0,\s*373\.0"),
    # --- 2026-09-24, the fold-in of the re-characterized fatal limits -------
    #     limit 2: where the diffusivity law is actually evaluated
    dict(id="N15", printed=96.8, unit="percent",
         where="sec_validation.tex / S15 window duration inside the measured range",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface",
         grep=r"96\.8|0\.9678"),
    dict(id="N16", printed=0.74, unit="-",
         where="sec_validation.tex / 7.3 / S15 surface-state leverage, high end",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface",
         grep=r"0\.29 to 0\.74"),
    dict(id="N17", printed=0.30, unit="-",
         where="sec_validation.tex / 7.3 / S15 largest fitted leverage",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface",
         grep=r"0\.304594"),
    dict(id="N18", printed=11.53, unit="-",
         where="S15 volume mean over the measured maximum (re-keyed 2026-09-27)",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface",
         grep=r"11\.53"),
    dict(id="N19", printed=0.72, unit="-",
         where="sec_validation.tex / 7.3 / S15 temperature extrapolation, fraction of the fitted span",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface",
         grep=r"0\.7196"),
    dict(id="N20", printed=3.18, unit="-",
         where="S15 prediction factor at the surface state",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface",
         grep=r"3\.17 to 3\.18|3\.180"),
    #     limit 3: the two domain extensions
    dict(id="N21", printed=0.070, unit="-",
         where="sec_validation.tex / 7.3 / S16 lowest measured water activity held",
         record="analysis/results_2026-09-19/item_19_domain_extensions_evidence",
         grep=r"0\.070"),
    dict(id="N22", printed=0.046973, unit="-",
         where="sec_validation.tex / S16 continuation upper edge",
         record="analysis/results_2026-09-19/item_19_domain_extensions_evidence",
         grep=r"0\.046973"),
    dict(id="N23", printed=9.6919e-5, unit="kg/kg",
         where="S16 underflow of the model's retained-water activity",
         record="analysis/results_2026-09-19/item_19_domain_extensions_evidence",
         grep=r"9\.6919"),
    dict(id="N24", printed=17.9, unit="percent",
         where="S16 hexane activity convention spread in retained loading",
         record="analysis/results_2026-09-19/item_19_domain_extensions_evidence",
         grep=r"17\.9"),
    dict(id="N25", printed=0.693147, unit="-",
         where="S16 fixed-ratio limit of the scaled Maxwell-Stefan flux",
         record="analysis/results_2026-09-19/item_19_domain_extensions_evidence",
         grep=r"0\.693147"),
    dict(id="N26", printed=1.21, unit="-",
         where="S16 branch over measurement, high end of the low-activity straddle",
         record="analysis/results_2026-09-19/item_19_domain_extensions_evidence",
         grep=r"1\.210|1\.21"),
    #     limit 1: the departure-tangent refusals against the arrival's own
    #     certified resolution (Class-A measurement, outside paper/analysis)
    dict(id="N27", printed=1.938e-9, unit="-",
         where="S14.3 tightest defensible resolution bound",
         record="../docs/evidence/d21_measurement_2026-09-24",
         grep=r"1\.938e-09"),
    dict(id="N28", printed=1.759e11, unit="1/K",
         where="S14.3 temperature derivative a 1e-2 floor would need",
         record="../docs/evidence/d21_measurement_2026-09-24",
         grep=r"1\.759e\+11"),
    dict(id="N29", printed=5.73, unit="1/K",
         where="S14.3 largest measured row derivative",
         record="../docs/evidence/d21_measurement_2026-09-24",
         grep=r"5\.73"),
    dict(id="N30", printed=4.25, unit="-",
         where="S14.3 ceiling of achieved residual over the resolution floor",
         record="../docs/evidence/d21_measurement_2026-09-24",
         grep=r"4\.252"),
    # --- 2026-09-24, the D22 mesh-flux fold-in ----------------------------
    #     item_20 builds and verifies a boundary-consistent ALE mesh flux
    #     beside the shipped wet-core energy step; item_21 measures what it
    #     changes in the campaign's own result items, which is nothing.
    dict(id="Q1", printed=2.397, unit="-",
         where="Table 1 / Sec. 4 / S7.2 reconstructed order, 12 percent, Linf",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"2\.3974"),
    dict(id="Q2", printed=2.030, unit="-",
         where="S7.2 reconstructed order, 12 percent, L2",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"2\.0297"),
    dict(id="Q3", printed=2.354, unit="-",
         where="Table 1 / Sec. 4 / S7.2 reconstructed order, 1.2 percent, Linf",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"2\.3541"),
    dict(id="Q4", printed=0.7339, unit="-",
         where="Table S5 shipped control at 1024 steps, Linf",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"0\.7339"),
    dict(id="Q5", printed=0.8365, unit="-",
         where="Table S5 shipped control at 1024 steps, L2",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"0\.8365"),
    dict(id="Q6", printed=0.0399, unit="-",
         where="Table S5 shipped mesh flux at 1.2 percent, Linf",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"0\.0399"),
    dict(id="Q7", printed=1.9963, unit="-",
         where="S7.2 fixed-radius order, re-measured on both boundary kinds",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"1\.9963"),
    dict(id="Q8", printed=1.157e-12, unit="-",
         where="S7.2 worst accumulated relative energy ledger, 128-step ladders",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"1\.157e-12"),
    dict(id="Q9", printed=0.219, unit="-",
         where="S7.2 worst cell Peclet number of the mesh flux",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"2\.1866e-1"),
    dict(id="Q10", printed=23.9, unit="-",
         where="S7.2 error ratio at 128 cells, 12 percent, 1024 steps",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"23\.9"),
    dict(id="Q11", printed=5.97, unit="-",
         where="S7.2 error ratio at 256 cells, 1.2 percent",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"5\.97"),
    dict(id="Q12", printed=5.455e-3, unit="K",
         where="S7.2 moving-grid backward-Euler floor at 48 cells",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"5\.455"),
    dict(id="Q13", printed=1.836e-3, unit="K",
         where="S7.2 anchored boundary-cell defect at 256 cells",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"1\.836e-3"),
    dict(id="Q14", printed=6.098e-4, unit="K",
         where="S7.2 linear boundary-cell defect at 256 cells",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"6\.098e-4"),
    dict(id="Q15", printed=3.01, unit="-",
         where="S7.2 ratio of the two boundary-cell residues",
         record="analysis/results_2026-09-19/item_20_mms_mesh_flux", grep=r"3\.01"),
    dict(id="Q16", printed=37.0, unit="-", word=True,
         where="S7.2 scored observables per mesh, spelled out on the page",
         record="analysis/results_2026-09-19/item_21_mesh_flux_trajectory_delta", grep=r"37 observables per mesh"),
    dict(id="Q17", printed=414.0, unit="-",
         where="Sec. 6.2 / S7.2 regime-map grid points, byte identical in both arms",
         record="analysis/results_2026-09-19/item_21_mesh_flux_trajectory_delta", grep=r"414 grid points"),
    # --- 2026-09-26, fourth referee read (REFEREE_AUDIT_4): M1, M2 and the
    # minor items.  A1-A12 against item_28 (the 95 per cent prediction factor
    # re-derived at the boundary read and the band re-marched), A13-A21 for
    # the uptake numbers now parenthesised and the figure-3.17-free headline
    # (item_26), A22-A26 for the last-tray total floors and the isotherm refit,
    # A27-A28 for the surface-edge residuals beside the space-resolved ones,
    # and Z10-Z16 for the superseded sentences.
    dict(id="A1", printed=3.140, page="3.140", unit="-",
         where="Sec. 5.2 / tab:fanermarch / S10.4 / S15 / S18 prediction factor at the boundary read, sunflower",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"\*\*F = 3\.140\*\*|\*\*3\.1401\*\*"),
    dict(id="A2", printed=3.132, unit="-",
         where="Sec. 5.2 / tab:fanermarch / S10.4 / S15 / S18 prediction factor at the boundary read, soybean",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"\*\*3\.1323\*\*"),
    dict(id="A3", printed=3.14, unit="-", where="Sec. 5.2 / Sec. 7.3 factor at the boundary read, rounded",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/factor_item28.json", grep=r"3\.14011"),
    dict(id="A4", printed=47.4, unit="percent",
         where="tab:basis / Sec. 5.2 / tab:fanermarch / S10.4 / tab:item28 95 per cent band, sunflower, slow edge",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/band_item28.json", grep=r"0\.4744"),
    dict(id="A5", printed=14.4, unit="percent", where="S10.4 / tab:item28 95 per cent band, sunflower, fast edge",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/band_item28.json", grep=r"0\.1439"),
    dict(id="A6", printed=61.2, unit="percent",
         where="tab:basis / Sec. 5.2 / tab:fanermarch / S10.4 / tab:item28 95 per cent band, soybean, slow edge",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/band_item28.json", grep=r"0\.6121"),
    dict(id="A7", printed=10.4, unit="percent",
         where="tab:basis / Sec. 5.2 / tab:fanermarch / S10.4 / tab:item28 95 per cent band, soybean, fast edge",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/band_item28.json", grep=r"0\.1043"),
    dict(id="A8", printed=0.20, page="0.20", unit="-",
         where="Sec. 5.2 / Sec. 7.3 / S15 / S18 leverage at the boundary read, low",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"\*\*0\.20-0\.69\*\*"),
    dict(id="A9", printed=0.69, unit="-",
         where="Sec. 5.2 / Sec. 7.3 / S15 / S18 leverage at the boundary read, high",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"\*\*0\.20-0\.69\*\*"),
    dict(id="A10", printed=1.03, unit="-",
         where="Sec. 5.2 / S18 lowest boundary loading over the measured minimum",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"\| 1\.032 \|"),
    dict(id="A11", printed=33248, page="33\\,248", unit="-", where="S18 identity checks of the re-marched edges",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"33,248 of 33,248"),
    dict(id="A12", printed=7.92e-4, unit="kg/kg", where="S18 lowest boundary loading in the sunflower window",
         record="analysis/results_2026-09-19/item_28_prediction_factor_boundary_read/RESULT.md", grep=r"7\.920e-4"),
    dict(id="A13", printed=18.6, unit="-",
         where="abstract / tab:basis / Sec. 5.3 / S9.5 / Sec. 8 form-mapping factor without figure 3.17",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"7\.6 to 18\.6 \(without 3\.17\)"),
    dict(id="A14", printed=18.4, unit="-",
         where="abstract / tab:basis / Sec. 5.3 / S9.5 / Sec. 8 the same, refined",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 7\.6–18\.4"),
    dict(id="A15", printed=0.114, unit="-", where="Sec. 5.3 fit residual, refined, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 0\.041–0\.114"),
    dict(id="A16", printed=0.56, unit="-", where="Sec. 5.3 converted time ratio, refined, low",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 0\.56–2\.20"),
    dict(id="A17", printed=2.20, page="2.20", unit="-", where="Sec. 5.3 converted time ratio, refined, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 0\.56–2\.20"),
    dict(id="A18", printed=0.57, unit="-", where="Sec. 5.3 declared (fixed) time ratio, refined, low",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 0\.57–18\.2"),
    dict(id="A19", printed=18.2, unit="-", where="Sec. 5.3 declared (fixed) time ratio, refined, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 0\.57–18\.2"),
    dict(id="A20", printed=24.3, unit="-", where="Sec. 5.3 declared (tabulated) time ratio, refined, high",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"prints 6\.1–24\.3"),
    dict(id="A21", printed=12, unit="percent", where="Sec. 5.3 / S19 isotherm refit, C differs by 12 to 32 per cent",
         record="analysis/results_2026-09-19/item_26_uptake_refinement/RESULT.md", grep=r"C by\s+12 to 32 per cent"),
    dict(id="A22", printed=4.5, unit="ppm", where="Sec. 5.4 / Sec. 8 / S9.9 last-tray total floor, low (110 C, y 0.003)",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21/floor_at_desolventizer_conditions.csv", grep=r"4\.5079"),
    dict(id="A23", printed=66.0, page="66.0", unit="ppm", where="Sec. 5.4 / Sec. 8 / S9.9 last-tray total floor, high (100 C, y 0.03)",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21/floor_at_desolventizer_conditions.csv", grep=r"66\.021"),
    dict(id="A24", printed=152, unit="ppm", where="Sec. 5.4 / Sec. 8 / S9.9 total floor at y 0.1, low (110 C)",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21/floor_at_desolventizer_conditions.csv", grep=r"152\.22"),
    dict(id="A25", printed=222, unit="ppm", where="Sec. 5.4 / Sec. 8 / S9.9 total floor at y 0.1, high (100 C)",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-21/floor_at_desolventizer_conditions.csv", grep=r"222\.45"),
    dict(id="A26", printed=60, word=True, unit="cells", where="Sec. 5.5 the wet marches' discretization",
         record="analysis/results_2026-09-19/item_14_derived_trend_benchmark/run_item14_with_water.py",
         grep=r"TAIL_CELLS = 60", raw=True),
    dict(id="A27", printed=0.141, unit="-", where="S10.7 (left Sec. 5.2 in the item_32 pass) surface edge on the borrowed leg, soybean",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface/RESULT.md", grep=r"0\.1413 and 0\.1102"),
    dict(id="A28", printed=0.110, page="0.110", unit="-", where="S10.7 (left Sec. 5.2 in the item_32 pass) surface edge on the borrowed leg, sunflower",
         record="analysis/results_2026-09-19/item_05_faner_comparison/rerun_2026-09-24_surface/RESULT.md", grep=r"0\.1413 and 0\.1102"),
    dict(id="Z10", printed=0, unit="-", where="the dry-reading crossing slip 1.75 (minor 1)",
         expect="absent", grep=r"1\.45\$--\$1\.75|1\.45\$ to \$1\.75"),
    dict(id="Z11", printed=0, unit="-", where="'the law has to be evaluated at one place' (M1)",
         expect="absent", grep=r"has to be evaluated at one place"),
    dict(id="Z12", printed=0, unit="-", where="'marched only against the outer-cell read' of the space-resolved solve (M1)",
         expect="absent", grep=r"space-resolved solve, marched only against"),
    dict(id="Z13", printed=0, unit="-", where="the main-text band at the superseded factor (M2)",
         expect="absent", grep=r"SIrange\{5\.7\}\{47\.6\}|SIrange\{10\.7\}\{61\.4\}|5\.7\$--\$47\.6"),
    dict(id="Z14", printed=0, unit="-", where="'Extended validation detail' and 'the validation experiments' (minor 8)",
         expect="absent", grep=r"Extended validation detail|validation experiments"),
    dict(id="Z15", printed=0, unit="-", where="the bold superseded column of Table S8 (minor 7)",
         expect="absent", grep=r"mathbf\{9\.9\}|mathbf\{16\.3\}"),
    dict(id="Z16", printed=0, unit="-", where="'tens to low hundreds' of ppm (minor 5)",
         expect="absent", grep=r"tens to low hundreds"),
    # --- 2026-09-27, the owner's ruling on the comparison's geometry -------
    # G: item_30 (the declared volume-to-surface sphere R_d = 0.74 x 0.885 mm).
    dict(id="G1", printed=7.9, unit="percent",
         where="abstract / tab:basis / Sec. 5.2 / tab:fanermarch / Sec. 8 / S9.10 window RMS at R_d, sunflower",
         record=_I30, grep=r"7\.9 and 18\.7 per cent"),
    dict(id="G2", printed=18.7, unit="percent",
         where="abstract / tab:basis / Sec. 5.2 / tab:fanermarch / Sec. 8 / S9.10 window RMS at R_d, soybean",
         record=_I30, grep=r"7\.9 and 18\.7 per cent"),
    dict(id="G3", printed=0.0791, unit="-", where="Table S32 reported cell, sunflower",
         record=_I30, grep=r"\*\*0\.0791\*\*"),
    dict(id="G4", printed=0.1868, unit="-", where="Table S32 / S21 R2 reported cell, soybean",
         record=_I30, grep=r"\*\*0\.1868\*\*"),
    dict(id="G5", printed=7.8, unit="percent", where="Sec. 5.2 / tab:fanermarch converged, sunflower",
         record=_I30, grep=r"converged values 7\.8 and 17\.8 per cent"),
    dict(id="G6", printed=17.8, unit="percent", where="Sec. 5.2 / tab:fanermarch converged, soybean",
         record=_I30, grep=r"converged values 7\.8 and 17\.8 per cent"),
    dict(id="G7", printed=0.0777, unit="-", where="S21 continuum, sunflower, low",
         record=_I30, grep=r"\*\*0\.0777 to 0\.0779\*\*"),
    dict(id="G8", printed=0.1786, unit="-", where="S21 continuum, soybean, high",
         record=_I30, grep=r"\*\*0\.1784 to 0\.1786\*\*"),
    dict(id="G9", printed=2.13, unit="order", where="Sec. 5.2 / S21 mesh order, high",
         record=_I30, grep=r"2\.13 \(sunflower\) and 2\.01, 2\.03, 2\.05"),
    dict(id="G10", printed=1.10, page="1.10", unit="order", where="Sec. 5.2 / S21 time order, low",
         record=_I30, grep=r"1\.10 to 1\.14 \(soybean\)"),
    dict(id="G11", printed=1.33, unit="order", where="Sec. 5.2 / S21 time order, high",
         record=_I30, grep=r"time 1\.29 to 1\.33"),
    dict(id="G12", printed=5.2, unit="percent", where="tab:basis / Sec. 5.2 / tab:fanermarch / S10.4 / S21 band at R_d, sunflower, low",
         record=_I30, grep=r"at most 5\.2 to 40\.0 %"),
    dict(id="G13", printed=40.0, page="40.0", unit="percent", where="band at R_d, sunflower, high",
         record=_I30, grep=r"at most 5\.2 to 40\.0 %"),
    dict(id="G14", printed=9.1, unit="percent", where="band at R_d, soybean, low",
         record=_I30, grep=r"at most 9\.1 to 49\.6 %"),
    dict(id="G15", printed=49.6, unit="percent", where="band at R_d, soybean, high",
         record=_I30, grep=r"at most 9\.1 to 49\.6 %"),
    dict(id="G16", printed=1.20, page="1.20", unit="-", where="Sec. 5.2 / S21 band minimum, D multiplier, sunflower",
         record=_I30, grep=r"0\.050 at D x 1\.20"),
    dict(id="G17", printed=1.51, unit="-", where="Sec. 5.2 / S21 band minimum, D multiplier, soybean",
         record=_I30, grep=r"0\.090 at D x 1\.51"),
    dict(id="G18", printed=0.655, unit="mm", where="Sec. 5.2 / S21 declared radius",
         record=_I30, grep=r"= 0\.655 mm"),
    dict(id="G19", printed=0.6549, unit="mm", where="S21 declared radius, four digits",
         record=_I30, grep=r"0\.6549 mm"),
    dict(id="G20", printed=1.310, page="1.310", unit="mm", where="tab:fanermarch / Table S35 declared diameter",
         record=_I30D + "/conventions_table.csv", grep=r"1\.3098"),
    dict(id="G21", printed=29.1, unit="percent", where="tab:fanermarch thesis sphere, sunflower",
         record=_I30, grep=r"thesis spheres 29\.1 and 37\.5"),
    dict(id="G22", printed=37.5, unit="percent", where="tab:fanermarch thesis sphere, soybean",
         record=_I30, grep=r"thesis spheres 29\.1 and 37\.5"),
    dict(id="G23", printed=13.3, unit="percent", where="tab:fanermarch psi times thesis sphere, sunflower",
         record=_I30, grep=r"\| 1\.443 \| 0\.133 \(10\)"),
    dict(id="G24", printed=19.7, unit="percent", where="tab:fanermarch psi times thesis sphere, soybean",
         record=_I30, grep=r"\| 1\.332 \| 0\.197 \(6\)"),
    dict(id="G25", printed=1.443, unit="mm", where="tab:fanermarch / Table S35 psi times thesis sphere, sunflower",
         record=_I30, grep=r"\| 1\.443 \|"),
    dict(id="G26", printed=1.332, unit="mm", where="tab:fanermarch / Table S35 psi times thesis sphere, soybean",
         record=_I30, grep=r"\| 1\.332 \|"),
    # G27 re-keyed 2026-09-27 (fifth referee read): Sec. 5.2 now prints the
    # thesis diameter 1.95 mm, not the 0.975 mm radius.
    dict(id="G27", printed=1.95, unit="mm", where="Sec. 5.2 / tab:fanermarch thesis diameter, sunflower",
         record=_I30D + "/conventions_table.csv", grep=r"1\.95,0\.975"),
    dict(id="G28", printed=0.90, page="0.90", unit="mm", where="Sec. 5.2 thesis sphere radius, soybean",
         record=_I30D + "/conventions_table.csv", grep=r"1\.8,0\.9,"),
    dict(id="G29", printed=0.32, unit="-", where="Sec. 5.2 / Sec. 7.3 / S21 Fourier number at the last scored sample",
         record=_I30, grep=r"about 0\.32\s+and 0\.33"),
    dict(id="G30", printed=0.033, unit="-", where="Sec. 5.2 / S21 largest population shift",
         record=_I30, grep=r"-0\.033 to \+0\.027"),
    dict(id="G31", printed=0.027, unit="-", where="S21 population shift up",
         record=_I30, grep=r"-0\.033 to \+0\.027"),
    dict(id="G32", printed=33330, page="33\\,330", unit="checks", where="S21 identity checks",
         record=_I30, grep=r"33,330 of 33,330"),
    dict(id="G33", printed=68, unit="marches", where="S21 marches run",
         record=_I30, grep=r"68 marches"),
    dict(id="G34", printed=0.0013, unit="-", where="S21 numerical uncertainty, sunflower",
         record=_I30, grep=r"\*\*\+0\.0013 to \+0\.0014\*\*"),
    dict(id="G35", printed=0.0083, unit="-", where="S21 claim, numerical uncertainty, soybean",
         record=_I30, grep=r"0\.0013 and 0\.0083\s+above"),
    dict(id="G36", printed=0.35, unit="-", where="S21 band width, sunflower",
         record=_I30, grep=r"0\.35 and 0\.40 of the span"),
    dict(id="G37", printed=3.137, unit="-", where="S21 factor at R_d, for information",
         record=_I30, grep=r"F = 3\.137"),
    dict(id="G38", printed=1.036, unit="-", where="S21 lowest boundary loading over the measured minimum",
         record=_I30, grep=r"lowest 1\.036 and 1\.048"),
    dict(id="G39", printed=0.029, unit="-", where="S21 R2 distance",
         record=_I30, grep=r"within 0\.029 of"),
    dict(id="G40", printed=268.6, unit="s", where="S21 time to the last measured loading, sunflower",
         record=_I30, grep=r"268\.6 and 292\.9"),
    dict(id="G41", printed=0.053, unit="-", where="S21 continuum mean signed residual, sunflower",
         record=_I30, grep=r"\+0\.053 and \+0\.167"),
    dict(id="G42", printed=0.0546, unit="-", where="Table S34 population (a1), sunflower",
         record=_I30, grep=r"0\.0546 \(15\)"),
    dict(id="G43", printed=0.2136, unit="-", where="Table S34 population (b), soybean",
         record=_I30, grep=r"0\.2136 \(4\)"),
    dict(id="G44", printed=0.00037, unit="-", where="S21 largest 60-to-240 change, sunflower",
         record=_I30, grep=r"0\.00037 \(sunflower\) and 0\.00130"),
    dict(id="G45", printed=0.00199, unit="-", where="S21 time-step change at 240 cells, soybean",
         record=_I30, grep=r"0\.00070 and 0\.00199"),
    dict(id="G46", printed=5.5e-10, unit="m2/s", where="S21 window geometric-mean coefficient at R_d",
         record=_I30, grep=r"5\.5e-10 and 5\.7e-10"),
    dict(id="G47", printed=2.25, unit="-", where="S20 / S21 slab against sphere slowest mode",
         record=_I30, grep=r"a factor 2\.25 apart"),
    dict(id="G48", printed=0.676, unit="-", where="S21 band minimum as a radius factor, sunflower",
         record=_I30, grep=r"0\.74/sqrt\(1\.20\) = 0\.676"),
    dict(id="G49", printed=0.0524, unit="-", where="Table S33 lowest marched inside the band, sunflower",
         record=_I30, grep=r"0\.0524 \(15/15\)"),
    dict(id="G50", printed=0.367, unit="-", where="Table S33 slow edge, sunflower",
         record=_I30, grep=r"0\.367 \(2/15\)"),
    # H: item_29 (the radius and the temperature dependence of the law).
    dict(id="H1", printed=16520, page="16\\,520", unit="checks", where="S20 identity checks",
         record=_I29, grep=r"16,520 of 16,520"),
    dict(id="H2", printed=0.677, unit="-", where="S20 residual minimum over the radius, sunflower",
         record=_I29, grep=r"0\.677 x 0\.885 = 0\.599 mm"),
    dict(id="H3", printed=0.601, unit="-", where="S20 residual minimum over the radius, soybean",
         record=_I29, grep=r"0\.601 x 0\.885 = 0\.532 mm"),
    dict(id="H4", printed=44.81, unit="kJ/mol", where="S20 / Table S31 implied activation energy, exact",
         record=_I29, grep=r"Ea = 44\.81 kJ/mol"),
    dict(id="H5", printed=44.8, unit="kJ/mol", where="Sec. 5.2 implied activation energy",
         record=_I29, grep=r"implies Ea = 44\.8\s+kJ/mol"),
    dict(id="H6", printed=4.294, unit="-", where="Table S31 factor at 136 C under the implied Ea",
         record=_I29, grep=r"\*\*4\.294\*\*"),
    dict(id="H7", printed=1.27, unit="-", where="S20 / Table S31 position in the band",
         record=_I29, grep=r"\*\*\+1\.27\*\*"),
    dict(id="H8", printed=32.65, unit="kJ/mol", where="S20 superseded hand value",
         record=_I29, grep=r"32\.65 kJ/mol"),
    dict(id="H9", printed=2.15, unit="-", where="S20 constant factor minimising the residual, sunflower",
         record=_I29, grep=r"\*\*2\.15 \(sunflower\) and 2\.78"),
    dict(id="H10", printed=2.78, unit="-", where="S20 constant factor minimising the residual, soybean",
         record=_I29, grep=r"\*\*2\.15 \(sunflower\) and 2\.78"),
    dict(id="H11", printed=0.198, unit="-", where="S20 / Table S31 Ea + 16.1, sunflower",
         record=_I29, grep=r"between 0\.198 and 0\.282"),
    dict(id="H12", printed=0.395, unit="-", where="S20 / Table S31 Ea - 16.1, soybean",
         record=_I29, grep=r"0\.326 and 0\.395"),
    dict(id="H13", printed=0.0093, unit="-", where="S20 cross-grid interpolation agreement",
         record=_I29, grep=r"at most 0\.0093"),
    dict(id="H14", printed=3.5e-13, unit="-", where="S20 loading path agreement under the similarity",
         record=_I29, grep=r"3\.5e-13"),
    dict(id="H15", printed=7.607e-4, unit="kg/kg", where="S20 boundary datum loading at 136 C",
         record=_I29, grep=r"7\.607e-4"),
    dict(id="H16", printed=0.065, unit="-", where="S20 slope per unit log of the Ea lever",
         record=_I29, grep=r"-0\.080\s+\(sunflower\) and -0\.065"),
    dict(id="H17", printed=0.49, unit="-", where="S20 slope per unit log of the radius",
         record=_I29, grep=r"0\.49 to 0\.55 of the span"),
    dict(id="H18", printed=1.102, unit="-", where="S20 thesis sphere radius factors",
         record=_I29, grep=r"factors 1\.102 and 1\.017"),
    dict(id="H19", printed=0.458, unit="-", where="Table S30 radius 1.2, soybean",
         record=_I29, grep=r"\| 0\.458 \| 0 \|"),
    dict(id="H20", printed=0.093, unit="-", where="Table S30 / S20 soybean minimum",
         record=_I29, grep=r"\*\*0\.093\*\* \| \*\*14\*\*"),
    dict(id="Z17", printed=0, unit="-", where="'a scalar fitted to one trace does better on both' (geometry decision)",
         expect="absent", grep=r"does better on\s+both"),
    dict(id="Z18", printed=0, unit="-", where="'both closer than the law' (geometry decision)",
         expect="absent", grep=r"both closer than the law"),
    dict(id="Z19", printed=0, unit="-", where="'Three further uncertainties attach' (geometry decision)",
         expect="absent", grep=r"Three further uncertainties attach"),
    # --- 2026-09-27, fifth referee read: each sample's own volume-to-surface
    # sphere, the declared geometry (item_31), and the sentences it withdrew.
    # I: item_31 (R_d = 0.74 x d_thesis / 2 = 0.7215 and 0.666 mm).
    dict(id="I1", printed=13.3, unit="percent",
         where="abstract / tab:basis / Sec. 5.2 / tab:fanermarch / Sec. 8 / S9.10 window RMS at the declared sphere, sunflower",
         record=_I31, grep=r"13\.3 and 19\.7 per cent of\s+the window span"),
    dict(id="I2", printed=19.7, unit="percent",
         where="abstract / tab:basis / Sec. 5.2 / tab:fanermarch / Sec. 8 / S9.10 window RMS at the declared sphere, soybean",
         record=_I31, grep=r"13\.3 and 19\.7 per cent of\s+the window span"),
    dict(id="I3", printed=0.7215, unit="mm", where="Sec. 5.2 / S22 declared radius, sunflower",
         record=_I31, grep=r"\*\*R_d,sun = 0\.74 x 1\.95 mm / 2 = 0\.7215 mm\*\*"),
    dict(id="I4", printed=0.666, unit="mm", where="Sec. 5.2 / S22 declared radius, soybean",
         record=_I31, grep=r"\*\*R_d,soy = 0\.74 x 1\.80 mm / 2 = 0\.6660 mm\*\*"),
    dict(id="I5", printed=13.0, page="13.0", unit="percent", where="Sec. 5.2 / tab:fanermarch converged, sunflower",
         record=_I31, grep=r"converged 13\.0\s+and 18\.9 per cent"),
    dict(id="I6", printed=18.9, unit="percent", where="Sec. 5.2 / tab:fanermarch converged, soybean",
         record=_I31, grep=r"converged 13\.0\s+and 18\.9 per cent"),
    dict(id="I7", printed=1.91, unit="order", where="Sec. 5.2 / S22 mesh order, low",
         record=_I31, grep=r"orders 1\.91 to 1\.98"),
    dict(id="I8", printed=2.04, unit="order", where="Sec. 5.2 / S22 mesh order, high",
         record=_I31, grep=r"2\.01, 2\.02, 2\.04 \(soybean\)"),
    dict(id="I9", printed=1.19, unit="order", where="Sec. 5.2 / S22 time order, high",
         record=_I31, grep=r"1\.16 to 1\.19 \(soybean\)"),
    dict(id="I10", printed=5.1, unit="percent", where="tab:basis / Sec. 5.2 / tab:fanermarch / S10.4 / S22 band, sunflower, low",
         record=_I31, grep=r"at most 5\.1 to 40\.4 %"),
    dict(id="I11", printed=40.4, unit="percent", where="band, sunflower, high",
         record=_I31, grep=r"at most 5\.1 to 40\.4 %"),
    dict(id="I12", printed=8.9, unit="percent", where="band, soybean, low",
         record=_I31, grep=r"at most 8\.9 to 50\.3 %"),
    dict(id="I13", printed=50.3, unit="percent", where="band, soybean, high",
         record=_I31, grep=r"at most 8\.9 to 50\.3 %"),
    dict(id="I14", printed=1.45, unit="-", where="Sec. 5.2 / S22 band minimum, D multiplier, sunflower",
         record=_I31, grep=r"0\.050 at D x 1\.45"),
    dict(id="I15", printed=1.56, unit="-", where="Sec. 5.2 / S22 band minimum, D multiplier, soybean",
         record=_I31, grep=r"0\.089 at D x 1\.56"),
    dict(id="I16", printed=1.034, unit="-", where="S22 lowest boundary loading over the measured minimum (main text: at or above 1.03)",
         record=_I31, grep=r"lowest 1\.034 and 1\.047"),
    dict(id="I17", printed=10.0, page="10.0", unit="percent", where="Sec. 5.2 Ea inside the half-width, sunflower, low",
         record=_I31, grep=r"10\.0 to 17\.7 per cent\s+\(sunflower\) and 18\.0 to 22\.4"),
    dict(id="I18", printed=17.7, unit="percent", where="Sec. 5.2 Ea inside the half-width, sunflower, high",
         record=_I31, grep=r"10\.0 to 17\.7 per cent"),
    dict(id="I19", printed=18.0, page="18.0", unit="percent", where="Sec. 5.2 Ea inside the half-width, soybean, low",
         record=_I31, grep=r"18\.0 to 22\.4 per cent"),
    dict(id="I20", printed=22.4, unit="percent", where="Sec. 5.2 Ea inside the half-width, soybean, high",
         record=_I31, grep=r"18\.0 to 22\.4 per cent"),
    dict(id="I21", printed=2.85, unit="-", where="Sec. 5.2 transported identification over the law at 136 C, sunflower",
         record=_I31, grep=r"factor 2\.854 and 2\.431"),
    dict(id="I22", printed=2.43, unit="-", where="Sec. 5.2 transported identification over the law at 136 C, soybean",
         record=_I31, grep=r"factor 2\.854 and 2\.431"),
    dict(id="I23", printed=32.3, unit="kJ/mol", where="Sec. 5.2 / S22 Ea implied by the transported identification, sunflower",
         record=_I31, grep=r"\*\*32\.3 and 27\.4 kJ/mol\*\*"),
    dict(id="I24", printed=27.4, unit="kJ/mol", where="Sec. 5.2 / S22 Ea implied by the transported identification, soybean",
         record=_I31, grep=r"\*\*32\.3 and 27\.4 kJ/mol\*\*"),
    dict(id="I25", printed=0.0329, unit="-", where="Table S39 / Sec. 5.2 largest population shift down (0.033)",
         record=_I31, grep=r"\*\*-0\.0329\*\*"),
    dict(id="I26", printed=0.0269, unit="-", where="Table S39 / Sec. 5.2 population shift up (0.027)",
         record=_I31, grep=r"\*\*\+0\.0269\*\*"),
    dict(id="I27", printed=0.3, unit="-", where="Sec. 5.2 / S22 Fourier number at the last scored sample",
         record=_I31, grep=r"Fourier number of about 0\.3"),
    dict(id="I28", printed=17043, page="17\\,043", unit="checks", where="S22 identity checks",
         record=_I31, grep=r"17,043 of 17,043 checks exact over 2 marches"),
    dict(id="I29", printed=82, unit="marches", where="S22 marches run",
         record=_I31, grep=r"82 marches"),
    dict(id="I30", printed=0.81525, unit="-", where="S22 radius factor, sunflower",
         record=_I31, grep=r"radius factor 0\.81525"),
    dict(id="I31", printed=0.75254, unit="-", where="S22 radius factor, soybean",
         record=_I31, grep=r"radius factor 0\.75254"),
    dict(id="I32", printed=1.9e-12, unit="-", where="S22 whole-march ledger",
         record=_I31, grep=r"every whole march to at most 1\.9e-12"),
    dict(id="I33", printed=0.1330, page="0.1330", unit="-", where="Table S36 reported cell, sunflower",
         record=_I31, grep=r"\*\*0\.1330\*\* \| \*\*10\*\*"),
    dict(id="I34", printed=0.1971, unit="-", where="Table S36 reported cell, soybean",
         record=_I31, grep=r"\*\*0\.1971\*\* \| \*\*6\*\*"),
    dict(id="I35", printed=0.1308, unit="-", where="Table S36 finest cell, full-step, sunflower",
         record=_I31, grep=r"\*\*0\.1308\*\*"),
    dict(id="I36", printed=0.1898, unit="-", where="Table S36 finest cell, full-step, soybean",
         record=_I31, grep=r"\*\*0\.1898\*\*"),
    dict(id="I37", printed=0.1304, unit="-", where="S22 continuum, sunflower",
         record=_I31, grep=r"\*\*0\.1304 to 0\.1305\*\*"),
    dict(id="I38", printed=0.1892, unit="-", where="S22 continuum, soybean",
         record=_I31, grep=r"\*\*0\.1889 to 0\.1892\*\*"),
    dict(id="I39", printed=0.0026, unit="-", where="S22 numerical uncertainty, sunflower (main text: 0.2 points)",
         record=_I31, grep=r"\*\*\+0\.0024 to \+0\.0026\*\*"),
    dict(id="I40", printed=0.0082, unit="-", where="S22 numerical uncertainty, soybean (main text: 0.8 points)",
         record=_I31, grep=r"\*\*\+0\.0079 to \+0\.0082\*\*"),
    dict(id="I41", printed=0.00132, unit="-", where="S22 largest 60-to-240 change, soybean",
         record=_I31, grep=r"0\.00070 \(sunflower\) and 0\.00132"),
    dict(id="I42", printed=0.00097, unit="-", where="S22 time-step change at 240 cells, sunflower",
         record=_I31, grep=r"0\.00097 and 0\.00200"),
    dict(id="I43", printed=0.111, unit="-", where="S22 continuum mean signed residual, sunflower",
         record=_I31, grep=r"\+0\.111 and \+0\.178"),
    dict(id="I44", printed=0.0508, unit="-", where="Table S37 lowest marched inside the band, sunflower",
         record=_I31, grep=r"0\.0508 \(15/15\)"),
    dict(id="I45", printed=0.0896, unit="-", where="Table S37 lowest marched inside the band, soybean",
         record=_I31, grep=r"0\.0896 \(14/15\)"),
    dict(id="I46", printed=0.404, unit="-", where="Table S37 slow edge, sunflower",
         record=_I31, grep=r"0\.404 \(1/15\)"),
    dict(id="I47", printed=0.334, unit="-", where="Table S37 fast edge, soybean",
         record=_I31, grep=r"0\.334 \(5/15\)"),
    dict(id="I48", printed=40.6, unit="percent", where="Table S37 note, decimated band, sunflower",
         record=_I31, grep=r"decimated: 5\.0 to\s+40\.6 and 9\.3 to 50\.5 %"),
    dict(id="I49", printed=3.138, unit="-", where="S22 factor at the declared radii, for information",
         record=_I31, grep=r"F = 3\.138 and 3\.130"),
    dict(id="I50", printed=0.072, unit="-", where="S22 Ea-pair slope, sunflower",
         record=_I31, grep=r"Ea pair -0\.072 \(sunflower\) and -0\.042"),
    dict(id="I51", printed=0.243, unit="-", where="S22 matched constant-pair slope, soybean",
         record=_I31, grep=r"matched constant pair -0\.190 and -0\.243"),
    dict(id="I52", printed=0.38, unit="-", where="S22 Ea over constant, sunflower",
         record=_I31, grep=r"Ea over constant:\s+0\.38 and 0\.17"),
    dict(id="I53", printed=5.5e-14, unit="-", where="S22 D/R^2 identity, soybean",
         record=_I31, grep=r"1\.2e-14 \(sunflower\) and 5\.5e-14 \(soybean\)"),
    dict(id="I54", printed=5.8e-16, unit="-", where="S22 D/R^2 identity under Ea + 16.1, soybean",
         record=_I31, grep=r"2\.1e-16 and 5\.8e-16"),
    dict(id="I55", printed=4.293, unit="-", where="Table S38 factor at 136 C under 44.81 kJ/mol (item_31's rounding)",
         record=_I31, grep=r"\| 44\.81 \| 4\.293 \| \+1\.27 \|"),
    dict(id="I56", printed=7.74e-10, unit="m2/s", where="S22 transported identification, sunflower",
         record=_I31, grep=r"7\.74e-10 \(sunflower\) and\s+6\.60e-10 m2/s"),
    dict(id="I57", printed=32.34, unit="kJ/mol", where="Table S38 transported implied Ea, sunflower",
         record=_I31, grep=r"32\.34 \(sun\), 27\.45 \(soy\)"),
    dict(id="I58", printed=2.854, unit="-", where="Table S38 transported factor, sunflower",
         record=_I31, grep=r"2\.854 \(sun\), 2\.431 \(soy\)"),
    dict(id="I59", printed=0.100, page="0.100", unit="-", where="Table S38 / S22 Ea + 16.1 at the declared radius, sunflower",
         record=_I31, grep=r"\*\*0\.100 \(12\)\*\*"),
    dict(id="I60", printed=0.224, unit="-", where="Table S38 / S22 Ea - 16.1 at the declared radius, soybean",
         record=_I31, grep=r"0\.224 \(3\)"),
    dict(id="I61", printed=0.073, unit="-", where="Table S38 constant matched to Ea + 16.1, sunflower",
         record=_I31, grep=r"\| 0\.073 \(15\) \|"),
    dict(id="I62", printed=0.355, unit="-", where="Table S38 constant matched to Ea - 16.1, soybean",
         record=_I31, grep=r"\| 0\.355 \(0\) \|"),
    dict(id="I63", printed=0.1025, unit="-", where="Table S39 population (a1), sunflower",
         record=_I31, grep=r"0\.1025 \(14\)"),
    dict(id="I64", printed=0.2240, page="0.2240", unit="-", where="Table S39 population (b), soybean",
         record=_I31, grep=r"0\.2240 \(3\)"),
    dict(id="I65", printed=1.1017, unit="-", where="Table S39 note / S22 class scaling, sunflower",
         record=_I31, grep=r"1\.1017\s+sunflower, 1\.0169 soybean"),
    dict(id="I66", printed=0.017, unit="-", where="S22 unscaled population (b), soybean",
         record=_I31, grep=r"\*\*\+0\.017\*\*"),
    dict(id="I67", printed=0.67, unit="-", where="S22 adversarial sunflower sphericity reproducing 7.9 per cent",
         record=_I31, grep=r"a sunflower sphericity of 0\.67"),
    dict(id="I68", printed=1.95, unit="mm", where="Sec. 5.2 / S22 thesis sunflower diameter",
         record=_I31, grep=r"1\.95 mm \(sunflower\) and 1\.80 mm \(soybean\)"),
    dict(id="Z20", printed=0, unit="-", where="'a derivation, not a fit' (fifth referee read, M3)",
         expect="absent", grep=r"derivation, not a fit"),
    dict(id="Z21", printed=0, unit="-", where="'the source particle's volume-to-surface sphere' (M1, minor 1)",
         expect="absent", grep=r"source particle's"),
    dict(id="Z22", printed=0, unit="-", where="'reproduces the sunflower trace inside the digitization band' at the declared sphere (M1)",
         expect="absent", grep=r"reproduces the\s+sunflower trace inside"),
    dict(id="Z23", printed=0, unit="-", where="'up when the mass-mean is 2.0 mm' (minor 2)",
         expect="absent", grep=r"up when\s+the mass-mean is"),
    dict(id="Z24", printed=0, unit="-", where="'halves the soybean miss' / 'miss is halved' (M1)",
         expect="absent", grep=r"halves the soybean miss|miss is\s+halved"),
    dict(id="Z25", printed=0, unit="-", where="'nothing in it is calibrated to the observation' (M3)",
         expect="absent", grep=r"nothing in it is\s+calibrated"),
    dict(id="Z26", printed=0, unit="-", where="the bold declared row of Table 3 (M3)",
         expect="absent", grep=r"textbf\{Volume-to-surface|textbf\{\\SI\{7\.9\}"),
    dict(id="Z27", printed=0, unit="-", where="'one mean sphericity for both thesis samples' without the diameter (M1)",
         expect="absent", grep=r"from one mean sphericity for both thesis samples"),
    dict(id="Z28", printed=0, unit="-", where="'at most $0.033$ of the span either way' (minor 2)",
         expect="absent", grep=r"at\s+most \$0\.033\$ of the span either way"),
    dict(id="Z29", printed=0, unit="-", where="'within 2 per cent of each other' for both species (minor 8)",
         expect="absent", grep=r"report within \\SI\{2\}\{\\percent\} of each other\. What"),
    dict(id="Z30", printed=0, unit="-", where="'monodisperse and spherical at a fixed radius' (minor 7)",
         expect="absent", grep=r"spherical at a fixed radius"),
    # --- 2026-09-27, sixth referee read (REFEREE_AUDIT_6_2026-09-27.md): the
    # borrowed sphericity quantified (M1), the comparator labelled (M2), the
    # temperature edges' signed residual and both factors (minors 1, 2).
    dict(id="L1", printed=0.49, unit="-", where="Sec. 5.2 residual per unit ln psi (per unit ln R at 0.885 mm), low",
         record=_I29, grep=r"radius \(0\.49 to 0\.55 of the span\)"),
    dict(id="L2", printed=0.55, unit="-", where="Sec. 5.2 residual per unit ln psi, high",
         record=_I29, grep=r"radius \(0\.49 to 0\.55 of the span\)"),
    dict(id="L3", printed=11.4, unit="percent", where="Sec. 5.2 / S22 soybean at the 0.575 mm sphere bracketing psi 0.67",
         record=_I29, grep=r"\| 0\.65\* \| 0\.575 \| 0\.057 \| 15 \| 0\.114 \| 12 \|"),
    dict(id="L4", printed=15.3, unit="percent", where="Sec. 5.2 / S22 soybean at the 0.620 mm sphere bracketing psi 0.67",
         record=_I29, grep=r"\| 0\.70 \| 0\.620 \| \*\*0\.055\*\* \| \*\*15\*\* \| 0\.153 \| 9 \|"),
    dict(id="L5", printed=0.575, unit="mm", where="Sec. 5.2 / S22 marched sphere below psi 0.67 (soybean)",
         record=_I29, grep=r"\| 0\.65\* \| 0\.575 \|"),
    dict(id="L6", printed=0.620, page="0.620", unit="mm", where="Sec. 5.2 / S22 marched sphere above psi 0.67 (soybean)",
         record=_I29, grep=r"\| 0\.70 \| 0\.620 \|"),
    dict(id="L7", printed=29.1, unit="percent", where="Sec. 5.2 / Sec. 8 / tab:basis sunflower at psi 1 (thesis volume-equivalent sphere)",
         record=_I31, grep=r"psi = 1 gives the thesis volume-equivalent\s+sphere's 29\.1 per cent"),
    dict(id="L8", printed=37.5, unit="percent", where="Sec. 5.2 / Sec. 8 / tab:basis soybean at psi 1 (thesis volume-equivalent sphere)",
         record=_I30D + "/conventions_table.csv",
         grep=r"faner2008_soybean_136C,\"volume-equivalent sphere, thesis\",1\.8,0\.9,1\.0169491525423728,0\.3750982"),
    dict(id="L9", printed=7.9, unit="percent", where="Sec. 5.2 / Sec. 8 sunflower at psi 0.67",
         record=_I31, grep=r"sunflower sphericity of 0\.67\s+would reproduce item 30's 7\.9 per cent"),
    dict(id="L10", printed=0.069, unit="-", where="Sec. 5.2 mean signed residual at Ea + 16.1, sunflower (minor 1)",
         record="analysis/results_2026-09-19/item_31_species_own_geometry/temperature_table.csv",
         grep=r"faner2008_sunflower_136C,arrh_ea_plus_half95,.*,0\.0689243"),
    dict(id="L11", printed=0.137, unit="-", where="Sec. 5.2 mean signed residual at Ea + 16.1, soybean (minor 1)",
         record="analysis/results_2026-09-19/item_31_species_own_geometry/temperature_table.csv",
         grep=r"faner2008_soybean_136C,arrh_ea_plus_half95,.*,0\.1373089"),
    dict(id="L12", printed=3.13, unit="-", where="Sec. 5.2 soybean 95 per cent factor (minor 2)",
         record=_I31, grep=r"soybean, F = 3\.1323"),
    dict(id="L13", printed=3.14, unit="-", where="Sec. 5.2 sunflower 95 per cent factor",
         record=_I31, grep=r"sunflower, F = 3\.1401"),
    dict(id="Z31", printed=0, unit="-", where="'anywhere inside the source's 95 % half-width' (minor 1)",
         expect="absent", grep=r"anywhere\s+inside the source's"),
    dict(id="Z32", printed=0, unit="-", where="S21 'because its mass-mean is 2.0 mm' (minor 6)",
         expect="absent", grep=r"because its mass-mean is"),
    dict(id="Z33", printed=0, unit="-", where="'whenever a 1.0 mm class carries' (minor 5)",
         expect="absent", grep=r"whenever a \\SI\{1\.0\}"),
    dict(id="Z34", printed=0, unit="-", where="abstract '(volume-equivalent: 24.4 and 36.6 %)' unlabelled (M2)",
         expect="absent", grep=r"\(volume-equivalent: \$24\.4\$ and \\SI\{36\.6\}"),
    dict(id="Z35", printed=0, unit="-", where="S10.6 'And both directions improve on the' unqualified (minor 7)",
         expect="absent", grep=r"And both directions improve on the"),
    dict(id="Z36", printed=0, unit="-", where="'$3.14$ at the boundary read' for both samples (minor 2)",
         expect="absent", grep=r"\$3\.14\$ at the boundary read"),
    dict(id="Z37", printed=0, unit="-", where="Sec. 8 caveat attached to the surface read (minor 3)",
         expect="absent", grep=r"surface read, a convention proposed"),
    dict(id="Z38", printed=0, unit="-", where="abstract 'on a declared volume-to-surface sphere' singular (minor 4)",
         expect="absent", grep=r"on a\s+declared volume-to-surface sphere chosen"),
    # --- 2026-09-26 (machine date; the day's earlier rows are stamped
    # 2026-09-27): the space-resolved solve marched at the declared spheres
    # and refined (item_32), Sec. 5.2, Sec. 7.3, Sec. 8, S10.7 and new S23.
    dict(id="S1", printed=0.736, unit="-", where="Sec. 5.2 / S23 space-resolved residual, sunflower, declared sphere, Whitaker leg",
         record=_I32, grep=r"\| Whitaker \| sunflower \| 0\.7215 mm \(own\) \| \*\*0\.7355\*\*, 2 of 15 \|"),
    dict(id="S2", printed=0.875, unit="-", where="Sec. 5.2 / S23 space-resolved residual, sunflower, declared sphere, borrowed leg",
         record=_I32, grep=r"\| borrowed \| sunflower \| 0\.7215 mm \(own\) \| \*\*0\.8748\*\*, 0 of 15 \|"),
    dict(id="S3", printed=88.7, unit="s", where="S23 soybean declared-sphere refusal, 60 cells, borrowed",
         record=_I32, grep=r"\| borrowed \| soybean \| 0\.6660 mm \(own\) \| \*\*refused at 88\.7 s\*\*"),
    dict(id="S4", printed=147.1, unit="s", where="S23 soybean declared-sphere refusal, 60 cells, Whitaker",
         record=_I32, grep=r"\| Whitaker \| soybean \| 0\.6660 mm \(own\) \| \*\*refused at 147\.1 s\*\*"),
    dict(id="S5", printed=3732, unit="-", where="S23 identity checks",
         record=_I32D + "/identity.json", grep=r"\"checks\": 3732"),
    dict(id="S6", printed=131.34, unit="W/m2/K", where="S23 borrowed film coefficient",
         record=_I32D + "/DETAIL.md", grep=r"h = 131\.34 W/m2K"),
    dict(id="S7", printed=96.89, unit="W/m2/K", where="S23 Whitaker film coefficient, sunflower",
         record=_I32D + "/DETAIL.md", grep=r"h = 96\.89 and 75\.71"),
    dict(id="S8", printed=75.71, unit="W/m2/K", where="S23 Whitaker film coefficient, soybean",
         record=_I32D + "/DETAIL.md", grep=r"h = 96\.89 and 75\.71"),
    dict(id="S9", printed=300.1, unit="s", where="S23 last scored sample after the start of drying",
         record=_I32, grep=r"last scored sample \(300\.1 s after"),
    dict(id="S10", printed=0.564, unit="-", where="S23 space-resolved residual, sunflower, 0.885 mm, Whitaker leg",
         record=_I32, grep=r"\| Whitaker \| sunflower \| 0\.885 mm \| 0\.5642, 5 of 15 \|"),
    dict(id="S11", printed=0.538, unit="-", where="S23 space-resolved residual, soybean, 0.885 mm, Whitaker leg",
         record=_I32, grep=r"\| Whitaker \| soybean \| 0\.885 mm \| 0\.5384, 6 of 15 \|"),
    dict(id="S12", printed=58.4, unit="s", where="Table S40 120-cell refusal, borrowed sunflower 0.885 mm",
         record=_I32, grep=r"\| borrowed \| sunflower \| 0\.885 mm \|.*\| refused at 58\.4 s \|"),
    dict(id="S13", printed=59.6, unit="s", where="Table S40 120-cell refusal, borrowed sunflower 0.7215 mm",
         record=_I32, grep=r"\| borrowed \| sunflower \| 0\.7215 mm \(own\) \|.*\| refused at 59\.6 s \|"),
    dict(id="S14", printed=52.2, unit="s", where="Table S40 120-cell refusal, borrowed soybean 0.885 mm",
         record=_I32, grep=r"\| borrowed \| soybean \| 0\.885 mm \|.*\| refused at 52\.2 s \|"),
    dict(id="S15", printed=40.5, unit="s", where="Table S40 120-cell refusal, borrowed soybean 0.666 mm",
         record=_I32, grep=r"\| borrowed \| soybean \| 0\.6660 mm \(own\) \|.*\| refused at 40\.5 s \|"),
    dict(id="S16", printed=81.9, unit="s", where="Table S40 120-cell refusal, Whitaker sunflower 0.885 mm",
         record=_I32, grep=r"\| Whitaker \| sunflower \| 0\.885 mm \|.*\| refused at 81\.9 s \|"),
    dict(id="S17", printed=82.8, unit="s", where="Table S40 120-cell refusal, Whitaker sunflower 0.7215 mm",
         record=_I32, grep=r"\| Whitaker \| sunflower \| 0\.7215 mm \(own\) \|.*\| refused at 82\.8 s \|"),
    dict(id="S18", printed=85.6, unit="s", where="Table S40 120-cell refusal, Whitaker soybean 0.885 mm",
         record=_I32, grep=r"\| Whitaker \| soybean \| 0\.885 mm \|.*\| refused at 85\.6 s \|"),
    dict(id="S19", printed=70.0, page="70.0", unit="s", where="Table S40 120-cell refusal, Whitaker soybean 0.666 mm",
         record=_I32, grep=r"\| Whitaker \| soybean \| 0\.6660 mm \(own\) \|.*\| refused at 70\.0 s \|"),
    dict(id="S20", printed=-0.37, unit="-", where="S23 mean signed residual where the march covers the window, one end",
         record=_I32, grep=r"negative \(-0\.37 to -0\.82 of the span\)"),
    dict(id="S21", printed=-0.82, unit="-", where="S23 mean signed residual where the march covers the window, other end",
         record=_I32, grep=r"negative \(-0\.37 to -0\.82 of the span\)"),
    dict(id="S22", printed=0.95, unit="-", where="S23 volume fraction above the measured loading range at the window start",
         record=_I32, grep=r"falls from 0\.95 to 1\.0 at the window start to 0 to 0\.42 at"),
    dict(id="S23", printed=0.42, unit="-", where="S23 volume fraction above the measured loading range at the window end",
         record=_I32, grep=r"falls from 0\.95 to 1\.0 at the window start to 0 to 0\.42 at"),
    dict(id="S24", printed=97, unit="-", where="S23 within-particle coefficient ratio in the window, low",
         record=_I32, grep=r"coefficient ratio reaches 97 to 101"),
    dict(id="S25", printed=101, unit="-", where="S23 within-particle coefficient ratio in the window, high",
         record=_I32, grep=r"coefficient ratio reaches 97 to 101"),
    dict(id="S26", printed=40, unit="s", where="S23 120-cell refusals, earliest (Sec. 5.2 says 'inside the window')",
         record=_I32, grep=r"Every 120-cell march\s+refuses 40 to 86 s"),
    dict(id="S27", printed=86, unit="s", where="S23 120-cell refusals, latest",
         record=_I32, grep=r"Every 120-cell march\s+refuses 40 to 86 s"),
    dict(id="S28", printed=45.6, unit="s", where="S23 extension refusals at the 5/4 s ceiling, earliest",
         record=_I32, grep=r"between 45\.6 and\s+144\.0 s"),
    dict(id="S29", printed=144.0, page="144.0", unit="s", where="S23 extension refusals at the 5/4 s ceiling, latest",
         record=_I32, grep=r"between 45\.6 and\s+144\.0 s"),
    dict(id="S30", printed=32.5, unit="s", where="S23 extension refusals at the 5/16 s ceiling, earliest",
         record=_I32, grep=r"between 32\.5 and 57\.4 s"),
    dict(id="S31", printed=57.4, unit="s", where="S23 extension refusals at the 5/16 s ceiling, latest",
         record=_I32, grep=r"between 32\.5 and 57\.4 s"),
    dict(id="S32", printed=0.003, unit="kg/kg", where="S23 outer-cell loading at the last accepted stride, low",
         record=_I32D + "/DETAIL.md", grep=r"outer cell at 0\.003 to 0\.012\s+kg/kg"),
    dict(id="S33", printed=0.012, unit="kg/kg", where="S23 outer-cell loading at the last accepted stride, high",
         record=_I32D + "/DETAIL.md", grep=r"outer cell at 0\.003 to 0\.012\s+kg/kg"),
    dict(id="S34", printed=15, unit="-", where="S23 coefficient ratio at the last accepted stride, low",
         record=_I32D + "/DETAIL.md", grep=r"coefficient varying\s+by 15 to 61"),
    dict(id="S35", printed=61, unit="-", where="S23 coefficient ratio at the last accepted stride, high",
         record=_I32D + "/DETAIL.md", grep=r"coefficient varying\s+by 15 to 61"),
    dict(id="S36", printed=9.6e-15, unit="-", where="S23 per-stride inventory ledger",
         record=_I32, grep=r"at most 9\.6e-15 of the initial inventory"),
    dict(id="S37", printed=7.7e-15, unit="-", where="S23 between-stride remap ledger",
         record=_I32, grep=r"to at most 7\.7e-15"),
    dict(id="S38", printed=6.7e-13, unit="-", where="S23 whole-march ledger",
         record=_I32, grep=r"to at most 6\.7e-13"),
    dict(id="S39", printed=0.770, page="0.770", unit="-", where="S23 paper configuration through the declared code path, soybean",
         record=_I32D + "/analysis.json", grep=r"\"paper_configuration_via_declared_code_path\": \{\s+\"faner2008_soybean_136C\": 0\.7700807507969353"),
    dict(id="S40", printed=0.697, unit="-", where="S23 paper configuration through the declared code path, sunflower",
         record=_I32D + "/analysis.json", grep=r"\"faner2008_sunflower_136C\": 0\.6967799577128457\s+\}"),
    dict(id="Z39", printed=0, unit="-", where="Sec. 5.2 'the one reading no evaluation choice affects' (item_32)",
         expect="absent", grep=r"the one reading no\s+evaluation choice affects"),
    dict(id="Z40", printed=0, unit="-", where="Sec. 7.3 / Sec. 8 'lies farther from both' (item_32: unrefined)",
         expect="absent", grep=r"lies farther\s+from both"),
    dict(id="Z41", printed=0, unit="-", where="S10.7 'has not been marched on Whitaker's' (item_32 marched it)",
         expect="absent", grep=r"has not been marched on Whitaker's"),
    # ----- 2026-09-27, owner's strategic clarity pass (STRATEGIC_PASS_2026-09-27.md)
    # T: the refinement ladder re-run on the current tree (item_33), now printed
    # in S14.3 and pointed to from Sec. 6.1.
    dict(id="T1", printed=36, unit="-", where="S14.3 / Sec. 6.1 re-run ladder, jobs",
         record=_I33, grep=r"The 36 jobs of"),
    dict(id="T2", printed=4, unit="-", word=True, where="S14.3 / Sec. 6.1 re-run ladder, four of 36 complete",
         record=_I33, grep=r"So \*\*four trajectories complete\*\*"),
    dict(id="T3", printed=15, unit="-", where="S14.3 re-run ladder, observables meeting their criterion on the one pair",
         record=_I33, grep=r"15 of the 22 observables meet their criterion"),
    dict(id="T4", printed=22, unit="-", where="S14.3 re-run ladder, observables scored on the one pair",
         record=_I33, grep=r"15 of the 22 observables meet their criterion"),
    dict(id="T5", printed=3, unit="-", word=True, where="S14.3 re-run ladder, cells that now complete",
         record=_I33, grep=r"refused pre-D20, \*\*completes now\*\* \| 3"),
    dict(id="T6", printed=4, unit="-", word=True, where="S14.3 re-run ladder, cells that now refuse",
         record=_I33, grep=r"completed pre-D20, \*\*refuses now\*\* \| 4"),
    # P: wording the strategic pass withdrew from every live source.
    dict(id="P1", printed=0, unit="-", where="abstract / Sec. 8 'Two outcomes disagree'",
         expect="absent", grep=r"Two outcomes disagree"),
    dict(id="P2", printed=0, unit="-", where="abstract 'the source is silent'",
         expect="absent", grep=r"the source is\s+silent"),
    dict(id="P3", printed=0, unit="-", where="abstract / Sec. 8 'another source's coefficients'",
         expect="absent", grep=r"another source's coefficients"),
    dict(id="P4", printed=0, unit="-", where="abstract / Sec. 8 'an independent bed correlation'",
         expect="absent", grep=r"an independent\s+bed correlation"),
    dict(id="P5", printed=0, unit="-", where="Table 2 'a disagreement, identified as a form mapping'",
         expect="absent", grep=r"a disagreement, identified as a form"),
    dict(id="P6", printed=0, unit="-", where="Table 2 / Sec. 8 'the source does not say which'",
         expect="absent", grep=r"the source does not say"),
    dict(id="P7", printed=0, unit="-", where="abstract 'No moving-front trajectory meets the refinement criteria'",
         expect="absent", grep=r"No moving-front trajectory meets the\s+refinement criteria"),
    dict(id="P8", printed=0, unit="-", where="Sec. 1 'A definitional conflict compounds it'",
         expect="absent", grep=r"definitional conflict"),
    dict(id="P9", printed=0, unit="-", where="Table 2 / Sec. 8 'the temperature source's own moisture'",
         expect="absent", grep=r"the temperature source's own"),
    dict(id="P10", printed=0, unit="-", where="Sec. 5.1 / Sec. 7.3 / S10.8 'until the experimental authors clarify'",
         expect="absent", grep=r"until the experimental authors clarify"),
    dict(id="P11", printed=0, unit="-", where="Secs. 5 and 6 headings phrased as questions",
         expect="absent", grep=r"\\subsection\{(Does|Can|Is|Where is|Which|How fast)[^}]*\?\}"),
    dict(id="P12", printed=0, unit="-", where="Table 2 'against 100--500 ppm: a disagreement'",
         expect="absent", grep=r"\\SIrange\{100\}\{500\}\{\\ppm\}: a disagreement"),
]

# ---------------------------------------------------------------------------
# 2026-09-27, corrections fold-in (CORRECTIONS_FOLDIN_2026-09-27.md): version 2
# of the three digitizations.  Every assertion whose printed number moved is
# re-keyed in place to the version-2 value and to the re-scoring record that
# carries it (item_35: Faner (2008); item_36: Cardarelli (1998); item_37:
# Faner et al. (2019)); the vessel-gas floor is withdrawn (item_34) and its
# assertions now expect absence; new C-series assertions bind the numbers the
# fold-in put on the page (floors at published compositions, oil-arm
# intervals, the mechanism paragraph); Z42-Z58 are the withdrawn wordings.
# ---------------------------------------------------------------------------
_R = "analysis/results_2026-09-19/"
_I34 = _R + "item_34_published_vapour_floors/RESULT.md"
_I35 = _R + "item_35_faner2008_v2_rescoring/RESULT.md"
_I35D = _R + "item_35_faner2008_v2_rescoring/DETAIL.md"
_I35P = _R + "item_35_faner2008_v2_rescoring/outputs/printed_numbers_v1_v2.csv"
_I35O = _R + "item_35_faner2008_v2_rescoring/outputs"
_I36 = _R + "item_36_cardarelli_v2_rescoring/RESULT.md"
_I36P = _R + "item_36_cardarelli_v2_rescoring/outputs/printed_numbers.csv"
_I36O = _R + "item_36_cardarelli_v2_rescoring/outputs"
_I37 = _R + "item_37_faner2019_v2_rescoring/RESULT.md"
_I38 = _R + "item_38_steam_implications_note/NOTE.md"
_F19V2 = "analysis/datasets/faner2019_v2/DATASET_RECORD.md"
_F08V2 = "analysis/datasets/faner2008_thesis_95_108_v2/DATASET_RECORD.md"
_C98V2 = "analysis/datasets/cardarelli1998_thesis_ch3_v2/DATASET_RECORD.md"

#: id -> (printed on the page at version 2, record, grep)
REKEY_2026_09_27 = {
    # constant-rate leg (item_35 section 5 and 6)
    "B1": (0.53, _I35P, r'"0\.53, 0\.61"'),
    "B4": (1.68, _I35P, r'"1\.45, 1\.68"'),
    "F9": (12.7, _I35P, r'"12\.7, 16\.4"'),
    "F10": (16.4, _I35P, r'"12\.7, 16\.4"'),
    "F17": (14.0, _I35P, r'F17,16\.3,[^\n]*,16\.3,14\.0,14\.0,'),
    "F18": (11.8, _I35O + "/item05_tail_residuals.csv", r"frozen_rig_state,law_log_loading_surface,[^\n]*0\.1175288"),
    "F19": (17.82, _I35P, r'"16\.18, 17\.82"'),
    "A27": (0.135, _I35P, r'"0\.1354, 0\.0994"'),
    "A28": (0.099, _I35P, r'"0\.1354, 0\.0994"'),
    "B6": (0.682, _I35P, r'"0\.682, 0\.899"'),
    "B7": (0.899, _I35P, r'"0\.682, 0\.899"'),
    "S39": (0.682, _I35P, r'"0\.682, 0\.899"'),
    "S40": (0.899, _I35P, r'"0\.682, 0\.899"'),
    # the 0.885 mm sphere (items 27, 28)
    "K1": (31.7, _I35P, r'"31\.7, 30\.3"'),
    "K2": (30.3, _I35P, r'"31\.7, 30\.3"'),
    "K3": (31.2, _I35P, r'"31\.2, 29\.7"'),
    "K4": (29.7, _I35P, r'"31\.2, 29\.7"'),
    "K5": (1.97, _I35P, r'"1\.97, 2\.02"'),
    "K6": (2.02, _I35P, r'"1\.97, 2\.02"'),
    "K7": (1.02, _I35P, r'"1\.02, 1\.16"'),
    "K9": (0.5, _I35P, r'"0\.6, 0\.5"'),
    "K11": (17.0, _I35P, r'"44\.1, 17\.0"'),
    "K12": (44.1, _I35P, r'"44\.1, 17\.0"'),
    "K13": (19.7, _I35P, r'"39\.0, 19\.7"'),
    "K14": (39.0, _I35P, r'"39\.0, 19\.7"'),
    "K15": (61.9, _I35P, r'"61\.9, 18\.3, 6\.1"'),
    "K16": (51.8, _I35P, r'"51\.8, 9\.9"'),
    "K17": (18.3, _I35P, r'"61\.9, 18\.3, 6\.1"'),
    "K18": (6.1, _I35P, r'"61\.9, 18\.3, 6\.1"'),
    "K19": (9.9, _I35P, r'"51\.8, 9\.9"'),
    "K20": (0.56, _I35P, r'"0\.56, 0\.43"'),
    "K28": (25.4, _I35P, r'"25\.4, 25\.6"'),
    "K29": (25.6, _I35P, r'"25\.4, 25\.6"'),
    "A4": (61.7, _I35P, r'"61\.7, 17\.9, 51\.7, 9\.6"'),
    "A5": (17.9, _I35P, r'"61\.7, 17\.9, 51\.7, 9\.6"'),
    "A6": (51.7, _I35P, r'"61\.7, 17\.9, 51\.7, 9\.6"'),
    "A7": (9.6, _I35P, r'"61\.7, 17\.9, 51\.7, 9\.6"'),
    # radius sweep and temperature (item_29)
    "H2": (0.672, _I35P, r'"0\.672, 0\.612"'),
    "H3": (0.612, _I35P, r'"0\.672, 0\.612"'),
    "H9": (2.18, _I35P, r'"2\.18, 2\.64"'),
    "H10": (2.64, _I35P, r'"2\.18, 2\.64"'),
    "H11": (0.257, _I35P, r'"0\.257, 0\.366, 0\.267, 0\.329"'),
    "H12": (0.329, _I35P, r'"0\.257, 0\.366, 0\.267, 0\.329"'),
    "H13": (0.0146, _I35P, r'H13,0\.0093,[^\n]*,0\.0093,0\.0146,'),
    "H16": (0.058, _I35P, r'"-0\.103, -0\.058"'),
    "H17": (0.47, _I35P, r'"0\.47, 0\.66"'),
    "L1": (0.47, _I35P, r'"0\.47, 0\.66"'),
    "L2": (0.66, _I35P, r'"0\.47, 0\.66"'),
    "H19": (0.383, _I35P, r'"0\.383, 0\.077, 14\.000"'),
    "H20": (0.077, _I35P, r'"0\.383, 0\.077, 14\.000"'),
    "L3": (8.7, _I35P, r'"8\.7, 11\.9"'),
    "L4": (11.9, _I35P, r'"8\.7, 11\.9"'),
    # the 2019 meal's diameter (item_30)
    "G1": (9.9, _I35P, r'"0\.0993, 0\.1473"'),
    "G2": (14.7, _I35P, r'"0\.0993, 0\.1473"'),
    "G3": (0.0993, _I35P, r'"0\.0993, 0\.1473"'),
    "G4": (0.1473, _I35P, r'"0\.0993, 0\.1473"'),
    "G5": (9.4, _I35P, r'"0\.0940, 0\.0943, 0\.1428, 0\.1429"'),
    "G6": (14.3, _I35P, r'"0\.0940, 0\.0943, 0\.1428, 0\.1429"'),
    "G7": (0.0940, _I35P, r'"0\.0940, 0\.0943, 0\.1428, 0\.1429"'),
    "G8": (0.1429, _I35P, r'"0\.0940, 0\.0943, 0\.1428, 0\.1429"'),
    "G9": (2.11, _I35P, r'"2\.01, 2\.11"'),
    "G10": (1.07, _I35P, r'"1\.07, 1\.16"'),
    "G11": (1.16, _I35P, r'"1\.07, 1\.16"'),
    "G12": (5.4, _I35P, r'"5\.4, 51\.4, 7\.6, 41\.6"'),
    "G13": (51.4, _I35P, r'"5\.4, 51\.4, 7\.6, 41\.6"'),
    "G14": (7.6, _I35P, r'"5\.4, 51\.4, 7\.6, 41\.6"'),
    "G15": (41.6, _I35P, r'"5\.4, 51\.4, 7\.6, 41\.6"'),
    "G16": (1.21, _I35P, r'"1\.209, 1\.451, 0\.673"'),
    "G17": (1.45, _I35P, r'"1\.209, 1\.451, 0\.673"'),
    "G48": (0.673, _I35P, r'"1\.209, 1\.451, 0\.673"'),
    "G21": (37.9, _I35P, r'"37\.9, 31\.0"'),
    "G22": (31.0, _I35P, r'"37\.9, 31\.0"'),
    "L7": (37.9, _I35P, r'"37\.9, 31\.0"'),
    "L8": (31.0, _I35P, r'"37\.9, 31\.0"'),
    "G23": (17.2, _I35P, r'"17\.2, 15\.6"'),
    "G24": (15.6, _I35P, r'"17\.2, 15\.6"'),
    "G30": (0.035, _I35P, r'"-0\.035, 0\.033"'),
    "G31": (0.033, _I35P, r'"-0\.035, 0\.033"'),
    "G34": (0.0050, _I35P, r'"0\.0050, 0\.0044"'),
    "G35": (0.0044, _I35P, r'"0\.0050, 0\.0044"'),
    "G36": (0.46, _I35P, r'"0\.46, 0\.34"'),
    "G39": (0.033, _I35P, r'G39,0\.029,[^\n]*,0\.029,0\.033,'),
    "G41": (0.082, _I35P, r'"0\.082, 0\.124"'),
    "G42": (0.0639, _I35P, r'"0\.0639, 0\.1708"'),
    "G43": (0.1708, _I35P, r'"0\.0639, 0\.1708"'),
    "G44": (0.00100, _I35P, r'"0\.00100, 0\.00081, 0\.00155, 0\.00138"'),
    "G45": (0.00138, _I35P, r'"0\.00100, 0\.00081, 0\.00155, 0\.00138"'),
    "G49": (0.0542, _I35P, r'"0\.0542, 0\.4771"'),
    "G50": (0.477, _I35P, r'"0\.0542, 0\.4771"'),
    "L9": (9.9, _I35D, r"the R_d value itself is 9\.9"),
    # the declared spheres (item_31)
    "I1": (17.2, _I35, r"\*\*17\.2 % \(9 of 14\)\*\*"),
    "I2": (15.6, _I35, r"\*\*15\.6 % \(8 of 15\)\*\*"),
    "I33": (0.172, _I35P, r'"17\.2, 15\.6"'),
    "I34": (0.1562, _I35D, r"0\.1562\(8\)"),
    "I35": (0.1675, _I35P, r'"0\.1675, 0\.1523"'),
    "I36": (0.1523, _I35P, r'"0\.1675, 0\.1523"'),
    "I5": (16.7, _I35P, r'"0\.1667, 0\.1669, 0\.1517, 0\.1519"'),
    "I6": (15.2, _I35P, r'"0\.1667, 0\.1669, 0\.1517, 0\.1519"'),
    "I37": (0.1667, _I35P, r'"0\.1667, 0\.1669, 0\.1517, 0\.1519"'),
    "I38": (0.1519, _I35P, r'"0\.1667, 0\.1669, 0\.1517, 0\.1519"'),
    "I39": (0.0053, _I35P, r'"0\.0051, 0\.0053, 0\.0043, 0\.0045"'),
    "I40": (0.0045, _I35P, r'"0\.0051, 0\.0053, 0\.0043, 0\.0045"'),
    "I7": (2.00, _I35P, r'"2\.00, 2\.00, 2\.05, 2\.12"'),
    "I8": (2.12, _I35P, r'"2\.00, 2\.00, 2\.05, 2\.12"'),
    "I9": (1.16, _I35P, r'"1\.10, 1\.11, 1\.15, 1\.16"'),
    "I41": (0.00084, _I35P, r'"0\.00118, 0\.00084, 0\.00161, 0\.00140"'),
    "I42": (0.00161, _I35P, r'"0\.00118, 0\.00084, 0\.00161, 0\.00140"'),
    "I43": (0.158, _I35P, r'"0\.158, 0\.134"'),
    "I10": (5.4, _I35P, r'"5\.4, 52\.6, 7\.7, 42\.2"'),
    "I11": (52.6, _I35P, r'"5\.4, 52\.6, 7\.7, 42\.2"'),
    "I12": (7.7, _I35P, r'"5\.4, 52\.6, 7\.7, 42\.2"'),
    "I13": (42.2, _I35P, r'"5\.4, 52\.6, 7\.7, 42\.2"'),
    "I48": (52.8, _I35P, r'"5\.5, 52\.8, 7\.7, 42\.4"'),
    "I14": (1.46, _I35P, r'"1\.46, 1\.50"'),
    "I15": (1.50, _I35P, r'"1\.46, 1\.50"'),
    "I44": (0.0537, _I35P, r'"0\.0537, 0\.0773, 0\.5262, 0\.3013"'),
    "I45": (0.0773, _I35P, r'"0\.0537, 0\.0773, 0\.5262, 0\.3013"'),
    "I46": (0.526, _I35P, r'"0\.0537, 0\.0773, 0\.5262, 0\.3013"'),
    "I47": (0.301, _I35P, r'"0\.0537, 0\.0773, 0\.5262, 0\.3013"'),
    "I17": (13.0, _I35P, r'"0\.130, 0\.227, 0\.142, 0\.182"'),
    "I18": (22.7, _I35P, r'"0\.130, 0\.227, 0\.142, 0\.182"'),
    "I19": (14.2, _I35P, r'"0\.130, 0\.227, 0\.142, 0\.182"'),
    "I20": (18.2, _I35P, r'"0\.130, 0\.227, 0\.142, 0\.182"'),
    "I59": (0.130, _I35P, r'"0\.130, 0\.227, 0\.142, 0\.182"'),
    "I60": (0.182, _I35P, r'"0\.130, 0\.227, 0\.142, 0\.182"'),
    "L10": (0.097, _I35P, r'"0\.097, 0\.098"'),
    "L11": (0.098, _I35P, r'"0\.097, 0\.098"'),
    "I61": (0.084, _I35P, r'"0\.084, 0\.293"'),
    "I62": (0.293, _I35P, r'"0\.084, 0\.293"'),
    "I50": (0.091, _I35P, r'"-0\.091, -0\.039, -0\.258, -0\.196, 0\.355, 0\.197"'),
    "I51": (0.196, _I35P, r'"-0\.091, -0\.039, -0\.258, -0\.196, 0\.355, 0\.197"'),
    "I52": (0.35, _I35P, r'"-0\.091, -0\.039, -0\.258, -0\.196, 0\.355, 0\.197"'),
    "I53": (4.8e-14, _I35O + "/shadow_v2_rule/item_31_species_own_geometry/similarity_table.csv", r"4\.8\d*e-14"),
    "I25": (0.041, _I35P, r'"-0\.0406, 0\.0322"'),
    "I26": (0.032, _I35P, r'"-0\.0406, 0\.0322"'),
    "I63": (0.1315, _I35P, r'"0\.1315, 0\.1797"'),
    "I64": (0.1797, _I35P, r'"0\.1315, 0\.1797"'),
    "I66": (0.015, _I35P, r'I66,\+0\.017,[^\n]*,0\.017,0\.015,'),
    # the space-resolved solve (item_32 re-scored in item_35)
    "S1": (0.956, _I35P, r'"0\.956, 1\.129"'),
    "S2": (1.129, _I35P, r'"0\.956, 1\.129"'),
    "S10": (0.734, _I35P, r'"0\.734, 0\.473"'),
    "S11": (0.473, _I35P, r'"0\.734, 0\.473"'),
    "S20": (-0.34, _I35P, r'"-0\.34, -1\.07"'),
    "S21": (-1.07, _I35P, r'"-0\.34, -1\.07"'),
    "S22": (0.90, _I35O + "/spaceresolved_rescore.csv", r"0\.903"),
    # the uptake benchmark (item_36)
    "U1": (8.52, _I36, r"8\.52 to 19\.54"),
    "U5": (0.0105, _F19V2, r"0\.0105"),
    "A18": (0.59, _I36, r"0\.589–17\.51"),
    "W22": (20.053, _I36P, r"20\.05306"),
    "W25": (19.844, _I36P, r"19\.84437"),
    "W28": (2.459, _I36P, r"2\.45859"),
    # the journal traces and the oil arm (item_37)
    "D20": (31.8, _I37, r"31\.8 % \(121 to 24\.5\)"),
    "D21": (33.4, _I37, r"33\.4 % \(44\.2 to 29\.6\)"),
    "O8": (3.61, _I37, r"\*\*3\.61 \(0\.90 to 6\.32\)\*\*"),
    "O9": (2.50, _I37, r"\*\*2\.50 \(1\.79 to 2\.98\)\*\*"),
    "J7": (0.492364, _I37, r"0\.490628 → 0\.492364"),
    "J8": (0.50164, _I37, r"0\.499279 → 0\.50164"),
    "J9": (1.117, _I37, r"layer\s+depth 1\.116 → 1\.117"),
    "J10": (0.957, _I37, r"depth 0\.968 → 0\.957"),
    "J11": (0.702, _I37, r"gas velocity 0\.688 → 0\.702"),
    "J12": (0.610, _I37, r"gas velocity 0\.615 → 0\.610"),
    "J13": (0.157, _I37, r"initial loading 0\.170 → 0\.157"),
    "J14": (0.034, _I37, r"initial loading 0\.034 → 0\.034"),
    "J15": (0.7506, _I37, r"−0\.7551 to 1\.5716 → −0\.7506 to 1\.5716"),
    "J17": (0.7864, _I37, r"−0\.7864 to 0\.9194"),
    "J19": (0.0264, _I37, r"\*\*−0\.0264\*\*, now in the signed band"),
    "J20": (0.0109, _I37, r"\*\*−0\.0109\*\*"),
    "J21": (0.908, _I37, r"0\.906 → 0\.908"),
    "J22": (0.815, _I37, r"0\.823 → 0\.815"),
}

#: the vessel-gas floor and its oil-arm companion are withdrawn (item_34)
WITHDRAWN_2026_09_27 = {
    "U8": r"1216\.4|1216\.358",
    "U10": r"1216\.4",
    "O10": r"2085\.5",
}

for _a in ASSERTIONS:
    _k = _a["id"]
    if _k in REKEY_2026_09_27:
        _p, _rec, _g = REKEY_2026_09_27[_k]
        _a.update(printed=_p, record=_rec, grep=_g)
        _a.pop("page", None)
        _a["where"] = _a["where"] + " [re-keyed to version 2, 2026-09-27]"
    if _k in WITHDRAWN_2026_09_27:
        _a.update(expect="absent", grep=WITHDRAWN_2026_09_27[_k], printed=0)
        _a["where"] = _a["where"] + " [withdrawn, item_34, 2026-09-27]"

ASSERTIONS += [
    # C: numbers the fold-in put on the page
    dict(id="C1", printed=1.06, unit="-", where="abstract / Table 2 / Sec. 5.1 duty ratio, high",
         record=_I35P, grep=r'"0\.82, 1\.06"'),
    dict(id="C2", printed=0.86, unit="-", where="abstract / Table 2 / Sec. 5.1 crossing ratio, low",
         record=_I35P, grep=r'"0\.86, 1\.12"'),
    dict(id="C3", printed=0.92, unit="-", where="Sec. 5.1 duty ratio, 2008 soybean",
         record=_I35P, grep=r'"1\.06, 0\.92"'),
    dict(id="C4", printed=1.07, unit="-", where="Sec. 5.1 Whitaker over the demand, high",
         record=_I35P, grep=r'"0\.87, 1\.07"'),
    dict(id="C5", printed=2.7, unit="-", where="Sec. 5.1 Bird and Bradshaw-Myers over the demand, high",
         record=_I35P, grep=r'"1\.78, 2\.74"'),
    dict(id="C6", printed=90.2, unit="W/m2/K", where="S10.3 demand, 2008 sunflower",
         record=_I35P, grep=r'"90\.2, 81\.7"'),
    dict(id="C7", printed=81.7, unit="W/m2/K", where="S10.3 demand, 2008 soybean",
         record=_I35P, grep=r'"90\.2, 81\.7"'),
    dict(id="C8", printed=1.074, unit="-", where="S10.3 Whitaker over the demand, 2008 sunflower",
         record=_I35P, grep=r'"1\.074, 0\.927"'),
    dict(id="C9", printed=156.5, unit="W/m2/K", where="S10.3 dry-reading demand, 2008 sunflower",
         record=_I35P, grep=r'"141\.8, 156\.5"'),
    dict(id="C10", printed=3292.2, unit="W/kg", where="S13.2 specific duty, 2008 sunflower",
         record=_I35P, grep=r'"3292\.2, 3532\.1"'),
    dict(id="C11", printed=12.6, unit="percent", where="S9.14 / S22 15-point sunflower window",
         record=_I35, grep=r"12\.6 % \(9 of 15\)"),
    dict(id="C12", printed=0.20027, unit="kg/kg", where="S9.14 the 60.04 s sunflower point",
         record=_F08V2, grep=r"X = 0\.20027"),
    dict(id="C13", printed=70.11, unit="s", where="S9.14 / S10.2 sunflower window start",
         record=_F08V2, grep=r"14 points from 70\.11 s"),
    dict(id="C14", printed=9.1, unit="percent", where="Table 3 / S21 sunflower at 0.670 mm (version 1 points)",
         record=_I37, grep=r"0\.670 mm for sunflower \(9\.1 per cent, 14 of 15\)"),
    dict(id="C15", printed=8.4, unit="-", where="abstract / Table 2 / Sec. 5.3 form-mapping factor, low",
         record=_I36, grep=r"\*\*8\.4 to 20\.1\*\*"),
    dict(id="C16", printed=20.1, unit="-", where="abstract / Table 2 / Sec. 5.3 form-mapping factor, high",
         record=_I36, grep=r"\*\*8\.4 to 20\.1\*\*"),
    dict(id="C17", printed=20.0, unit="-", where="Table 2 / Sec. 5.3 form-mapping factor, refined, high",
         record=_I36, grep=r"8\.4 to 20\.0 refined"),
    dict(id="C18", printed=124, unit="points", where="Sec. 5.3 / S10.5 independent points",
         record=_I36, grep=r"12 series, 124 points"),
    dict(id="C19", printed=0.098, unit="-", where="Sec. 5.3 fit rms, high",
         record=_I36, grep=r"rms 0\.013 to 0\.098"),
    dict(id="C20", printed=0.085, unit="-", where="Sec. 5.3 analytic sphere rms, high",
         record=_I36, grep=r"0\.015 to 0\.085"),
    dict(id="C21", printed=0.146, unit="-", where="Sec. 5.3 converted rms, high",
         record=_I36, grep=r"0\.033 to 0\.146"),
    dict(id="C22", printed=2.12, unit="-", where="Sec. 5.3 converted time ratio, high",
         record=_I36, grep=r"0\.86 to 2\.12"),
    dict(id="C23", printed=2.5, unit="-", where="Sec. 5.3 converted band",
         record=_I36, grep=r"contracts to 2\.5"),
    dict(id="C24", printed=52, unit="points", where="Sec. 5.3 converted in box",
         record=_I36, grep=r"52 of the 124 points"),
    dict(id="C25", printed=22, unit="points", where="Sec. 5.3 fixed coefficient in box",
         record=_I36, grep=r"22 and 3 of the 124"),
    dict(id="C26", printed=11.6, unit="ppm", where="Sec. 5.4 / Sec. 7 / S9.11 sparged-section entry floor",
         record=_I34, grep=r"\*\*11\.6\*\* \| \*\*20\.3\*\*"),
    dict(id="C27", printed=20.3, unit="ppm", where="Sec. 5.4 / S9.11 the same with the arm",
         record=_I34, grep=r"\*\*11\.6\*\* \| \*\*20\.3\*\*"),
    dict(id="C28", printed=47.1, unit="ppm", where="Table 2 / Sec. 5.4 Fig. 5.3 vapours, with the arm, high",
         record=_I34, grep=r"47\.1 / 21\.7 / 12\.0"),
    dict(id="C29", printed=6.8, unit="ppm", where="Table 2 / Sec. 5.4 Fig. 5.3 vapours, arm-free, low",
         record=_I34, grep=r"26\.8 / 12\.4 / 6\.8"),
    dict(id="C30", printed=26.8, unit="ppm", where="Sec. 8 / S9.8 published gas, arm-free, high",
         record=_I37, grep=r"3\.9 to\s+26\.8 ppm without the arm"),
    dict(id="C31", printed=6.9, unit="ppm", where="Sec. 5.4 model vapour of 2002, with the arm",
         record=_I34, grep=r"it floors at 3\.9 to 6\.9 ppm"),
    dict(id="C32", printed=6.57, unit="-", where="Sec. 5.4 / S9.11 soybean final loading over the arm-free floor",
         record=_I37, grep=r"\*\*6\.57 \(1\.64 to 11\.50\)\*\*"),
    dict(id="C33", printed=3.96, unit="-", where="Sec. 5.4 / S9.11 sunflower final loading over the arm-free floor",
         record=_I37, grep=r"\*\*3\.96 \(2\.83 to 4\.71\)\*\*"),
    dict(id="C34", printed=0.0080, page="0.0080", unit="kg/kg", where="Sec. 5.4 / S9.9 / S9.11 soybean final loading",
         record=_F19V2, grep=r"soybean 0\.0080 ± 0\.006"),
    dict(id="C35", printed=2.35e-3, unit="kg/kg", where="S9.11 soybean ratio bounded away from one above",
         record=_I37, grep=r"m > 2\.35e-3"),
    dict(id="C36", printed=9.86e-3, unit="kg/kg", where="S9.11 sunflower mole-basis floor",
         record=_I37, grep=r"that floor is 9\.86e-3"),
    dict(id="C37", printed=61.3, unit="degC", where="Sec. 7.2 mechanism paragraph, co-boiling point",
         record=_I38, grep=r"334\.48 K \(61\.33 °C\)"),
    dict(id="C38", printed=0.054, unit="kg/kg", where="Sec. 7.2 mechanism paragraph, water load at the co-boiling point",
         record=_I38, grep=r"0\.054 kg of water per kg of hexane"),
    dict(id="C39", printed=0.129, unit="kg/kg", where="Sec. 7.2 mechanism paragraph, water load at 75 C",
         record=_I38, grep=r"0\.129 at 75 °C"),
    dict(id="C40", printed=105.84, unit="degC", where="Sec. 7.2 mechanism paragraph, condensation stops",
         record=_I38, grep=r"105\.84 °C at 19\.0 per cent"),
    dict(id="C41", printed=98, unit="percent", where="Sec. 7.2 mechanism paragraph, share of the duty (tower lane, measured)",
         record=_I38, grep=r"than 98 per cent of the duty came from the sparge steam"),
    dict(id="C42", printed=0.908, unit="-", where="Table S19 march on version 2, soybean",
         record=_I37, grep=r"Whitaker predicted / measured \| 0\.823 → 0\.815 \| 0\.906 → 0\.908"),
    dict(id="C43", printed=0.9194, unit="-", where="Table S19 sunflower residual span, high",
         record=_I37, grep=r"−0\.7864 to 0\.9194"),
    dict(id="C44", printed=0.94, unit="deg", where="S9.14 page rotation, low end",
         record=_C98V2, grep=r"\*\*-0\.94\*\*"),
    dict(id="C45", printed=0.028, unit="kg/kg", where="S9.14 largest version-2 move of a Figure 4.21 point",
         record=_F08V2, grep=r"\+0\.0080 to \+0\.0281"),
    # Z: wording the fold-in withdrew from every live source
    dict(id="Z42", printed=0, unit="-", where="abstract 'plant validation remains open' (the author's closing sentence)",
         expect="absent", grep=r"plant\s+validation remains open"),
    dict(id="Z43", printed=0, unit="-", where="Sec. 5.4 / S9.6 'the one vapour composition ... report for the industrial vessel' (item_34)",
         expect="absent", grep=r"one vapour\s+composition"),
    dict(id="Z44", printed=0, unit="-", where="S9.1 'logged every 2 s' (faner2019_v2)",
         expect="absent", grep=r"logged every"),
    dict(id="Z45", printed=0, unit="-", where="Sec. 5.2 'borrows both factors from the soybean meal' (item_37)",
         expect="absent", grep=r"borrows both\s+factors"),
    dict(id="Z46", printed=0, unit="-", where="sphericity 'measured on industrial soybean meal' (item_37)",
         expect="absent", grep=r"measured (by \\citet\{faner2019kinetics\} )?on\s+industrial\s+soybean\s+meal(:| and)"),
    dict(id="Z47", printed=0, unit="-", where="the version-1 headline '13.3 and 19.7'",
         expect="absent", grep=r"13\.3\$ and \\SI\{19\.7\}|13\.3\}\{\\percent\} and \\SI\{19\.7"),
    dict(id="Z48", printed=0, unit="-", where="the version-1 duty range '0.82 to 0.98'",
         expect="absent", grep=r"0\.82\$ ?(to|--) ?\$0\.98"),
    dict(id="Z49", printed=0, unit="-", where="the version-1 uptake factor '7.6 to 18.6'",
         expect="absent", grep=r"7\.6\$\s*(to|--)\s*\$18\.6"),
    dict(id="Z50", printed=0, unit="-", where="'the three bracket the demand at all four conditions' (item_35 M2)",
         expect="absent", grep=r"bracket\s+the coefficient the measured duty\s+demands at all four"),
    dict(id="Z51", printed=0, unit="-", where="'set-of-spheres figures' (item_36 sample wording)",
         expect="absent", grep=r"set-of-spheres"),
    dict(id="Z52", printed=0, unit="-", where="'the soybean set is a packed cylinder read as a set of spheres'",
         expect="absent", grep=r"packed cylinder read as a set"),
    dict(id="Z53", printed=0, unit="-", where="'thirty-five of the hundred' (version-1 uptake count)",
         expect="absent", grep=r"thirty-five of the hundred"),
    dict(id="Z54", printed=0, unit="-", where="'The vessel-gas floor and the last-tray floors are different quantities'",
         expect="absent", grep=r"vessel-gas floor and the last-tray floors"),
    dict(id="Z55", printed=0, unit="-", where="'at the vessel's reported gas' (item_34)",
         expect="absent", grep=r"vessel's (own )?reported"),
    dict(id="Z56", printed=0, unit="-", where="the version-1 final loading 7.415e-3 as a current value (item_37)",
         expect="absent", grep=r"7\.415\\times10\^\{-3\}"),
    dict(id="Z57", printed=0, unit="-", where="Sec. 5.2 'closer on sunflower than on soybean' (item_35: the traces trade places)",
         expect="absent", grep=r"closer on\s+sunflower than on soybean"),
    dict(id="Z58", printed=0, unit="-", where="the version-1 declared-sphere counts '10 and 6 of 15'",
         expect="absent", grep=r"\$10\$ and\s+\$6\$ of \$15\$"),
]

# ---------------------------------------------------------------------------
# 2026-09-27, final polish (READER_CHECK_RESPONSE_2026-09-27.md): the journal
# conditions' stage 1 on version 2 (item_39 section 1), item_29's S2
# re-evaluated (item_39 section 2), the pre-registered bracketing criterion
# quoted and declared departed from (item_39 section 3; owner ruling A,
# LIMITS_AUDIT addendum 13).
# ---------------------------------------------------------------------------
_I39 = _R + "item_39_journal_ranges_v2_and_outcomes/RESULT.md"
_I39C = _R + "item_39_journal_ranges_v2_and_outcomes/outputs/ranges_before_after.csv"
_I39L = _R + "item_39_journal_ranges_v2_and_outcomes/outputs/log_journal_ranges_v2.txt"

#: id -> (printed on the page after the polish, record, grep)
REKEY_POLISH_2026_09_27 = {
    "F20": (1690.8, _I39C, r"1690\.7739563620464"),
    "N7": (16.47, _I39C, r"16\.468894992239445"),
    "B4": (1.70, _I39C, r"1\.6954969156220498"),
}
for _a in ASSERTIONS:
    _k = _a["id"]
    if _k in REKEY_POLISH_2026_09_27:
        _p, _rec, _g = REKEY_POLISH_2026_09_27[_k]
        _a.update(printed=_p, record=_rec, grep=_g)
        _a.pop("page", None)
        _a["where"] = _a["where"] + " [re-keyed to item_39, 2026-09-27]"
        if _k == "B4":
            _a["page"] = "1.70"

ASSERTIONS += [
    dict(id="P1", printed=1.13, unit="-", where="abstract / Table 2 / Sec. 5.1 / Sec. 8 / S10.8 / S16 crossing ratio, high",
         record=_I39, grep=r"0\.86 to 1\.12 -> 0\.86 to 1\.13"),
    dict(id="P2", printed=0.86, unit="-", where="Sec. 5.1 / S10.3 Whitaker over the demand, low",
         record=_I39, grep=r"0\.87 to 1\.07 -> 0\.86 to 1\.07"),
    dict(id="P3", printed=92.9, unit="W/m2/K", where="S10.3 / S13.2 / Table S7 demand, high (2019 sunflower)",
         record=_I39, grep=r"91\.6 -> 92\.9"),
    dict(id="P4", printed=88.3, unit="W/m2/K", where="Table S7 demand, 2019 soybean",
         record=_I39L, grep=r"86\.52025633777734 -> 88\.30631706632722"),
    dict(id="P5", printed=131.8, unit="W/m2/K", where="Table S7 dry-reading demand, 2019 soybean",
         record=_I39L, grep=r"128\.96951666426838 -> 131\.78516856237235"),
    dict(id="P6", printed=139.5, unit="W/m2/K", where="Table S7 dry-reading demand, 2019 sunflower",
         record=_I39L, grep=r"137\.348884292203 -> 139\.54563251774033"),
    dict(id="P7", printed=16.53, unit="W/m2/K", where="S13.2 particle-area demand, 2019 soybean",
         record=_I39L, grep=r"16\.177352583302923 -> 16\.530535216576247"),
    dict(id="P8", printed=17.50, page="17.50", unit="W/m2/K", where="S13.2 particle-area demand, 2019 sunflower",
         record=_I39L, grep=r"17\.228422541912526 -> 17\.503972699038204"),
    dict(id="P9", printed=2431.8, unit="W/kg", where="S13.2 specific duty, 2019 soybean",
         record=_I39L, grep=r"2433\.024768739569 -> 2431\.7899835114126"),
    dict(id="P10", printed=0.969, unit="-", where="S18 crossing ratio, 2019 soybean",
         record=_I39L, grep=r"0\.9734866726619311 -> 0\.9694699569988898"),
    dict(id="P11", printed=1.129, unit="-", where="S18 crossing ratio, 2019 sunflower",
         record=_I39L, grep=r"1\.1215754069848691 -> 1\.1290967979156452"),
    dict(id="P12", printed=3.7, unit="percent", where="Sec. 5.1 / S10.3 tabulated coefficient below Whitaker, 2019 soybean",
         record=_I39, grep=r"\*\*75\.73 below 78\.65\*\* \(3\.7 %\)"),
    dict(id="P13", printed=4.2, unit="percent", where="Sec. 5.1 / S10.3 tabulated coefficient below Whitaker, 2019 sunflower",
         record=_I39, grep=r"\*\*76\.67 below 80\.01\*\* \(4\.2 %\)"),
    dict(id="P14", printed=0.0146, unit="-", where="S20 S2 cross-grid difference, sunflower (fires)",
         record=_I39, grep=r"\*\*0\.01463\*\*"),
    dict(id="P15", printed=0.0065, page="0.0065", unit="-", where="S20 S2 cross-grid difference, soybean",
         record=_I39, grep=r"0\.00652 \| 0\.00652"),
    dict(id="P16", printed=3e-14, unit="-", where="S20 / S9.14 exact-partner agreement",
         record=_I39, grep=r"holds to 3e-14"),
    dict(id="P17", printed=5.4e-13, unit="-", where="S20 loading-path agreement, printed as 5e-13",
         record=_I39, grep=r"agrees to 5\.4e-13"),
    dict(id="P18", printed=0.65, unit="-", where="S20 radius factor of the S2 firing",
         record=_I39, grep=r"at f_R = 0\.65"),
    dict(id="P19", printed=1.06, unit="-", where="Sec. 5.1 / S10.3 Whitaker's marched duty over the measured, 2008 sunflower",
         record=_I39C, grep=r"0\.9815780148887441,1\.0568546803333532"),
    dict(id="P20", printed=1.074, unit="-", where="S10.3 Whitaker over the demand, 2008 sunflower (reading ii)",
         record=_I39, grep=r"Whitaker at 1\.074 times"),
    # Z: ranges and wording the polish withdrew
    dict(id="Z59", printed=0, unit="-", where="the crossing range '0.86 to 1.12' (item_39: 1.13)",
         expect="absent", grep=r"0\.86\$\s*(to|--)\s*\$1\.12"),
    dict(id="Z60", printed=0, unit="-", where="the dry crossing range '1.45 to 1.68' (item_39: 1.70)",
         expect="absent", grep=r"1\.45\$\s*(to|--)\s*\$1\.68"),
    dict(id="Z61", printed=0, unit="-", where="Whitaker over the demand '0.87 to 1.07' (item_39: 0.86)",
         expect="absent", grep=r"0\.87\$\s*(to|--)\s*\$1\.07"),
    dict(id="Z62", printed=0, unit="-", where="'charge-mass basis' (READER_CHECK A2, S9)",
         expect="absent", grep=r"charge-mass\s+basis"),
    dict(id="Z63", printed=0, unit="-", where="'the stage-1 posing ... fires' (decision A: declared departure)",
         expect="absent", grep=r"stage-1 posing(,| \(the constant-rate leg\),) (therefore )?fires"),
    dict(id="Z64", printed=0, unit="-", where="the demand range '81.7 to 91.6' (item_39: 92.9)",
         expect="absent", grep=r"\\SIrange\{81\.7\}\{91\.6\}"),
    dict(id="Z65", printed=0, unit="-", where="Table S7 'inputs of version~1 of theirs' (item_39)",
         expect="absent", grep=r"inputs of version~1 of theirs"),
    dict(id="Z66", printed=0, unit="-", where="S20 'the interpolation's own error' (item_39 section 2)",
         expect="absent", grep=r"the interpolation's own error"),
    dict(id="Z67", printed=0, unit="-", where="S18 journal crossings '0.973 and 1.122' (item_39)",
         expect="absent", grep=r"\$0\.973\$ and \$1\.122\$"),
    dict(id="Z68", printed=0, unit="-", where="'that source's Figure' (READER_CHECK source item 1)",
         expect="absent", grep=r"that source's Figure"),
]

# ---------------------------------------------------------------------------
# 2026-09-27, the condensation of the main text (CONDENSATION_PLAN_2026-09-27.md).
# No number changed and no record moved.  Twelve record assertions whose
# number the main text no longer prints are now checked where the
# Supplementary Material prints it: the paragraphs the condensation shortened
# stand there in full, verbatim, under "Moved here from Sec. X of the main text
# (condensation of 2026-09-27)".  The coverage pass already searches the main
# text and the supplement together; this map only re-labels ``where`` so that
# the report says which document carries the number now.
# ---------------------------------------------------------------------------
#: id -> where the number is printed since the condensation
RELOCATED_CONDENSATION_2026_09_27 = {
    "F1": "S9.2 (moved paragraph of Sec. 5.1); Table 2 prints the deviations",
    "U2a": "S10.5 (moved paragraph of Sec. 5.3)",
    "U2b": "S10.5 (moved paragraph of Sec. 5.3)",
    "U7": "S9.11 (moved paragraph of Sec. 5.4) and S9.9",
    "D3": "S9.13 (moved paragraph of the former Sec. 5.6)",
    "U9": "S9.11 (moved paragraph of Sec. 5.4)",
    "U11": "S9.9 and S9.11 (moved paragraph of Sec. 5.4)",
    "O6": "S9.11 (moved paragraph of Sec. 5.4)",
    "O7": "S9.11 (moved paragraph of Sec. 5.4)",
    "L5": "S22 (moved paragraph of Sec. 5.2)",
    "L6": "S22 (moved paragraph of Sec. 5.2)",
    "S1": "S23 and S22 (moved paragraph of Sec. 5.2)",
}
for _a in ASSERTIONS:
    _k = _a["id"]
    if _k in RELOCATED_CONDENSATION_2026_09_27:
        _a["where"] = (_a["where"] + " [main text condensed 2026-09-27: now printed in "
                       + RELOCATED_CONDENSATION_2026_09_27[_k] + "]")
    # F13 (the regression's 0.410 in ln D) stays in Sec. 5.2.  Its coverage
    # token was the float 0.41, which the page's 0.410 cannot match and which
    # the condensed main text no longer prints elsewhere (the activity 0.41 of
    # the experiments' own gas moved to S9.11); the literal is named here so
    # that coverage finds the regression's own printed value.
    if _k == "F13":
        _a["page"] = "0.410"

# ---------------------------------------------------------------------------
# 2026-09-28, the condensation of the Supplementary Material
# (CONDENSATION_SUPPLEMENT_PLAN_2026-09-28.md).  No number changed and no
# record moved.  The 77 record assertions below were printed, before the
# condensation, only in passages of the supplement that now stand in the
# extended technical supplement alone (supplementary_extended.tex, deposited
# with the data and code).  Each is re-pointed to the section of that report
# where its number is first printed (E<n> is S<n> of the 119-page supplement,
# found by searching the assertion's own page tokens).  Two of them are still
# printed in the submitted supplement at another precision that the coverage
# tokens do not match: F5 as 1.2570 (S10.2) and V3 as
# 3.422887528969909e-10 (S7).  The coverage pass accepts an assertion printed
# only in the extended report when, and only when, it is named here.
# ---------------------------------------------------------------------------
#: id -> extended-report section printing the number since 2026-09-28
RELOCATED_SUPPLEMENT_CONDENSATION_2026_09_28 = {
    "F5": "E17", "V3": "E7", "R9": "E9.6", "R8": "E9.8", "R19": "E14.4", "R20": "E14.4",
    "U12": "E9.6", "U14": "E10.4", "O14": "E9.11", "J6": "E10.9", "J7": "E10.9", "J8": "E10.9",
    "J24": "E10.9", "J26": "E10.9", "J28": "E10.9", "K21": "E7.2", "W2": "E19", "W9": "E19",
    "W14": "E22", "W15": "E19", "W17": "E18", "W19": "E19", "W20": "E19", "W29": "E20",
    "W31": "E19", "N13": "E9.6", "N15": "E15", "N25": "E16.2", "N28": "E14.4", "N29": "E14.4",
    "Q2": "E7.1", "Q10": "E7.2", "Q11": "E7.2", "Q12": "E7.2", "Q13": "E7.2", "Q14": "E7.2",
    "Q15": "E7.2", "A10": "E17", "A11": "E18", "A12": "E18", "A20": "E22", "G3": "E21",
    "G9": "E21", "G32": "E21", "G33": "E21", "G35": "E9.12", "G41": "E21", "G42": "E21",
    "G47": "E20", "G48": "E21", "G49": "E21", "H6": "E20", "H8": "E20", "H9": "E20",
    "H10": "E20", "H15": "E20", "H16": "E20", "I16": "E17", "I21": "E22", "I22": "E22",
    "I31": "E22", "I42": "E22", "I43": "E22", "I51": "E22", "I53": "E22", "I56": "E22",
    "S22": "E3.4", "S24": "E23", "S27": "E23", "S28": "E10.8", "S29": "E23", "S35": "E23",
    "S38": "E23", "C36": "E9.11", "C45": "E9.14", "P10": "E17", "P15": "E20",
}
for _a in ASSERTIONS:
    _k = _a["id"]
    if _k in RELOCATED_SUPPLEMENT_CONDENSATION_2026_09_28:
        _a["where"] = (_a["where"] + " [supplement condensed 2026-09-28: now printed in the extended report, section "
                       + RELOCATED_SUPPLEMENT_CONDENSATION_2026_09_28[_k] + "]")


# ---------------------------------------------------------------------------
# 2026-09-30, the paper-1 rewrite R-P1 (REFEREE_AUDIT_7_2026-09-30.md): the
# refinement ladder re-founded on the F36 farm campaign in two labelled tiers
# (owner ruling D1 of 2026-09-30), the audit items A1-A10 and B1-B8, the
# comparator C1 of the surface read against the weighted mean of Crank.
#
# Records.  Four evidence folders had not landed in the tree when this block
# was written; their assumed landed paths are the constants below, each also
# listed in REFEREE_AUDIT_7_2026-09-30.md.  An assertion whose record lies
# under a root that does not exist yet is reported AWAITING-LANDED-EVIDENCE
# and counted apart; setting the named environment variable to the staging
# folder checks it now (the rewrite was checked that way before the landing).
# ---------------------------------------------------------------------------
_F36 = "../docs/evidence/farm_f36_2026-09-29"          # the F36 collection's evidence/
_F36R = "../docs/GT_PS2_F36_LADDER_COLLECTION_RECORD_2026-09-30.md"  # the F36 collection record
_C1 = "../docs/evidence/paper1_crank_comparator_2026-09-29"  # scratchpad crank_comparator/
_C1R = "../docs/GT_PS2_PAPER1_CRANK_COMPARATOR_RECORD_2026-09-29.md"  # the comparator record (its RESULT.md)
_VN = "../docs/evidence/paper1_referee_audit_2026-09-29/verify_N"  # scratchpad paper_audit/verify_N/
_VP = "../docs/evidence/paper1_referee_audit_2026-09-29/verify_P"  # scratchpad paper_audit/verify_P/
_SEED = "../docs/evidence/face_tangent_inward_seed_2026-09-29"  # landed
_SEEDR = "../docs/GT_PS2_FACE_TANGENT_INWARD_SEED_BUILD_RECORD_2026-09-29.md"  # landed
_MF1 = "../docs/evidence/face_departure_m_f1_2026-09-29"  # landed
_P02 = "../docs/GT_PS2_PART02_EXTINCTION_LOCALIZER_AND_DRY_CONTINUATION_2026-08-22.md"  # landed
_P01 = "../docs/GT_PS2_PART01_RECONVENED_CEREMONY_PACKET_2026-08-22.md"  # landed
_I34D = "analysis/results_2026-09-19/item_34_published_vapour_floors/DETAIL.md"
_I03P = "analysis/results_2026-09-19/item_03_coupled_water_pore/rerun_2026-09-20/pytest_coupled_suites.txt"

#: root -> environment variable naming its staging folder until it lands
PENDING_ROOTS = {
    _F36: "DTDC_PAPER_F36_EVIDENCE",
    _F36R: "DTDC_PAPER_F36_RECORD",
    _C1: "DTDC_PAPER_C1_EVIDENCE",
    _VN: "DTDC_PAPER_VERIFY_N",
    _VP: "DTDC_PAPER_VERIFY_P",
}

ASSERTIONS += [
    # --- the ladder on the present formulation (F36), Sec. 6.1, S14.1, S14.3, Fig. S1
    dict(id="Y1", printed=8.17e-4, unit="-", where="Sec. 6.1 / S14.1 front, correlated 12 against 24, first time (declared bound)",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.0008166409969210554"),
    dict(id="Y2", printed=9.73e-4, unit="-", where="Sec. 6.1 / S14.1 same pair, second time",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.000972755243434695"),
    dict(id="Y3", printed=9.38e-4, unit="-", where="Sec. 6.1 / S14.1 same pair, third time",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.0009377434504371801"),
    dict(id="Y4", printed=2.26e-3, unit="-", where="Sec. 6.1 / S14.1 front, composition 12 against 96, first time",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.0022639572711690985"),
    dict(id="Y5", printed=1.98e-3, unit="-", where="Sec. 6.1 / S14.1 same pair, second time",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.0019754792798941196"),
    dict(id="Y6", printed=1.57e-3, unit="-", where="Sec. 6.1 / S14.1 same pair, third time",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.0015722256011982206"),
    dict(id="Y7", printed=13, word=True, unit="cells", where="Sec. 6.1 / S14.3 complete under every criterion at the literal",
         record=_F36 + "/ladder/summary.json", grep=r'"gate_complete_by_arm": \{\s*"A": 13'),
    dict(id="Y8", printed=25, word=True, unit="cells", where="Sec. 6.1 / S14.3 complete under every criterion with the declared bound (13 + 9 + 3)",
         record=_F36 + "/ladder/summary.json", grep=r'"A": 13,\s*"B": 9,\s*"C": 3'),
    dict(id="Y9", printed=0.235, unit="-", where="Sec. 6.1 / S14.3 front ratio, literal mesh chain",
         record=_F36 + "/ladder/chains.csv", grep=r"literal,combined,mesh,0\.0375,.*,0\.23460694192631468,"),
    dict(id="Y10", printed=0.441, unit="-", where="Sec. 6.1 / S14.3 front ratio, literal temperature step chain",
         record=_F36 + "/ladder/chains.csv", grep=r"literal,temperature_only,time,48,.*,0\.4412757772605264,"),
    dict(id="Y11", printed=0.442, unit="-", where="Sec. 6.1 / S14.3 front ratio, literal correlated step chain",
         record=_F36 + "/ladder/chains.csv", grep=r"literal,combined,time,48,.*,0\.44187196026224324,"),
    dict(id="Y12", printed=2.09, unit="order", where="Sec. 6.1 / S14.3 observed order in space, literal",
         record=_F36 + "/ladder/chains.csv", grep=r"MONOTONE,2\.0916823921155885"),
    dict(id="Y13", printed=1.18, unit="order", where="Sec. 6.1 / S14.3 observed order in time, literal",
         record=_F36 + "/ladder/chains.csv", grep=r"MONOTONE,1\.180247538566519"),
    dict(id="Y14", printed=0.012, page="0.012", unit="percent", where="Sec. 6.1 / S14.3 GCI on the finest level, literal mesh chain",
         record=_F36 + "/ladder/chain_observables.csv", grep=r"final_front_position_over_radius,.*,0\.00012400882813959775,True"),
    dict(id="Y15", printed=0.05, page="0.05", unit="percent", where="Sec. 6.1 / S14.3 GCI on the finest level, literal step chains",
         record=_F36 + "/ladder/chain_observables.csv", grep=r"0\.0005122682278630624,True"),
    dict(id="Y16", printed=3.46, unit="order", where="Sec. 6.1 / S14.3 observed order in space, declared tier, high",
         record=_F36 + "/ladder/chains.csv", grep=r"MONOTONE,3\.4645216925972853"),
    dict(id="Y17", printed=0.97, unit="order", where="Sec. 6.1 / S14.3 observed order in time, declared tier, low",
         record=_F36 + "/ladder/chains.csv", grep=r"MONOTONE,0\.966301434564368"),
    dict(id="Y18", printed=2.75, unit="order", where="Sec. 6.1 / S14.3 observed order in time, declared tier, high",
         record=_F36 + "/ladder/chains.csv", grep=r"MONOTONE,2\.752509093032342"),
    dict(id="Y19", printed=0.004, page="0.004", unit="percent", where="Sec. 6.1 / S14.3 GCI, declared tier, low (3.54e-5)",
         record=_F36 + "/ladder/chain_observables.csv", grep=r"3\.543171914333638e-05,True"),
    dict(id="Y20", printed=0.052, page="0.052", unit="percent", where="Sec. 6.1 / S14.3 GCI, declared tier, high (5.16e-4)",
         record=_F36 + "/ladder/chain_observables.csv", grep=r"0\.0005159685535594883,True"),
    dict(id="Y21", printed=0.20, page="0.20", unit="percent", where="Sec. 6.1 / S14.3 worst finest front pair, declared tier",
         record=_F36 + "/ladder/chains.csv", grep=r"7\.280727949542739,0\.0019955985271276153"),
    dict(id="Y22", printed=14.8, unit="percent", where="Sec. 6.1 / S14.3 worst finest pulse water-flux pair, N <= 48",
         record=_F36 + "/ladder/observable_summary.csv", grep=r"pulse_water_flux_time_integral_mol_s_m2,10,0,0,0,0\.148"),
    dict(id="Y23", printed=16.7, unit="percent", where="S14.3 worst finest pulse water-flux pair, N = 96 chains",
         record=_F36 + "/ladder/observable_summary.csv", grep=r"n96ext,pulse_water_flux_time_integral_mol_s_m2,3,0,0,0,0\.167"),
    dict(id="Y24", printed=7.28, unit="-", where="Sec. 6.1 / S14.3 front ratio, composition 12-cell step chain (fires)",
         record=_F36 + "/ladder/summary.json", grep=r"ratio 7\.28"),
    dict(id="Y25", printed=2.39, unit="-", where="Sec. 6.1 / S14.3 front ratio, 24-48-96 composition (fires)",
         record=_F36 + "/ladder/summary.json", grep=r"ratio 2\.39"),
    dict(id="Y26", printed=6.84, unit="-", where="Sec. 6.1 / S14.3 front ratio, 24-48-96 temperature (fires)",
         record=_F36 + "/ladder/summary.json", grep=r"ratio 6\.84"),
    dict(id="Y27", printed=8.70, page="8.70", unit="-", where="Sec. 6.1 / S14.3 front ratio, 24-48-96 correlated (fires)",
         record=_F36 + "/ladder/chains.csv", grep=r"8\.696768281759088"),
    dict(id="Y28", printed=0.074, page="0.074", unit="percent", where="Sec. 6.1 / S14.3 finest 48-96 front pair, high",
         record=_F36 + "/ladder/chains.csv", grep=r"8\.696768281759088,0\.0007400868063275764"),
    dict(id="Y29", printed=0.048, page="0.048", unit="percent", where="S14.3 finest 48-96 front pair, low",
         record=_F36 + "/ladder/chains.csv", grep=r"2\.3907209594519387,0\.0004796415571099159"),
    dict(id="Y30", printed=7.8e-13, unit="-", where="Sec. 6.1 / Table 1 largest ladder ledger (7.773e-13 in S14.3)",
         record=_F36 + "/ladder/summary.json", grep=r'"largest_ledger": 7\.772904082778828e-13'),
    dict(id="Y31", printed=7.773e-13, unit="-", where="S14.3 largest ladder ledger",
         record=_F36 + "/ladder/summary.json", grep=r'"largest_ledger": 7\.772904082778828e-13'),
    dict(id="Y32", printed=47, unit="-", where="S14.3 observables meeting on the finest pair, literal (47 of 66)",
         record=_F36 + "/ladder/summary.json", grep=r'"observables_final_pair_meeting": "47 of 66"'),
    dict(id="Y33", printed=154, unit="-", where="S14.3 observables meeting on the finest pair, declared tier (154 of 220)",
         record=_F36 + "/ladder/summary.json", grep=r'"observables_final_pair_meeting": "154 of 220"'),
    dict(id="Y34", printed=44, unit="-", where="S14.3 observables meeting on the finest pair, N = 96 chains (44 of 66)",
         record=_F36 + "/ladder/summary.json", grep=r'"observables_final_pair_meeting": "44 of 66"'),
    dict(id="Y35", printed=8.027, unit="percent", where="S14.3 worst relative, correlated 12 against 24 (F36)",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.08027281145424683"),
    dict(id="Y36", printed=0.3673, unit="K", where="S14.3 worst temperature, same pair (F36)",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.36729138683023166"),
    dict(id="Y37", printed=29.97, unit="percent", where="S14.3 worst relative, composition 12 against 96 (F36)",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.2996594982604303"),
    dict(id="Y38", printed=0.7588, unit="K", where="S14.3 worst temperature, same pair (F36)",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.7587597160173232"),
    dict(id="Y39", printed=4.779, unit="-", where="S14.3 face-speed jump, composition N = 48, dt 0.075",
         record=_F36 + "/ladder/jump_by_mesh.csv", grep=r"4\.778588555565851"),
    dict(id="Y40", printed=1.036, unit="-", where="S14.3 face-speed jump, N = 96, same physical face",
         record=_F36 + "/ladder/jump_by_mesh.csv", grep=r"1\.0363889018374868"),
    dict(id="Y41", printed=1.454, unit="-", where="S14.3 face-speed jump, composition N = 12, dt 0.075",
         record=_F36 + "/ladder/jump_by_mesh.csv", grep=r"1\.4540128432834523"),
    dict(id="Y42", printed=2.656, unit="-", where="S14.3 face-speed jump, temperature N = 48, dt 0.075",
         record=_F36 + "/ladder/jump_by_mesh.csv", grep=r"temperature_only,0\.075,.*,2\.65[56]"),
    dict(id="Y43", printed=8.51e-11, unit="-", where="S14.3 original-solve residual on the present ladder, low",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"8\.514893087468349e-11"),
    dict(id="Y44", printed=5.24e-2, unit="-", where="S14.3 original-solve residual on the present ladder, high",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"0\.052367305311826"),
    dict(id="Y45", printed=2.14e-11, unit="-", where="S14.3 best inward candidate, low",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"2\.1356592544897533e-11"),
    dict(id="Y46", printed=7.68e-11, unit="-", where="S14.3 best inward candidate, high",
         record=_F36 + "/ladder/printed_numbers_ladder.json", grep=r"7\.683620997145214e-11"),
    dict(id="Y47", printed=2.9e-11, unit="-", where="S14.3 best inward candidate of the ten refused N <= 48 cells, low (2.932e-11)",
         record=_F36 + "/checks.json", grep=r'"f36-033-f36an12d0375temp": 2\.932170603533585e-11'),
    dict(id="Y48", printed=7.3e-11, unit="-", where="S14.3 same, high (7.270e-11)",
         record=_F36 + "/checks.json", grep=r'"f36-023-f36an24d0375comp": 7\.270207881872931e-11'),
    dict(id="Y49", printed=5.295e-11, unit="-", where="S14.3 96-cell bootstrap residual at the literal",
         record=_F36 + "/jobs/f36-001-f36an96d01875comb/farm_stderr.txt", grep=r"maximum scaled residual=5\.295e-11 > 2\.000e-11"),
    dict(id="Y50", printed=65, unit="-", where="S14.3 in-process determinism audit bit-exact, 65 of 65",
         record=_F36 + "/headline.json", grep=r'"True\|mesh48": 11,\s*"True\|mesh24": 20,\s*"True\|mesh12": 34'),
    dict(id="Y51", printed=1.71e-11, unit="-", where="Sec. 6.1 (1.7e-11) / S14.3 cross-host checked-observable difference, low",
         record=_F36 + "/determinism/aggregate.json", grep=r"1\.705358016217157e-11"),
    dict(id="Y52", printed=6.27e-10, unit="-", where="Sec. 6.1 (6.3e-10) / S14.3 cross-host checked-observable difference, high",
         record=_F36 + "/determinism/aggregate.json", grep=r"6\.272659519058966e-10"),
    dict(id="Y53", printed=9.29e-10, unit="K", where="S14.3 cross-host temperature difference, max",
         record=_F36 + "/determinism/aggregate.json", grep=r"9\.289919944421854e-10"),
    dict(id="Y54", printed=28, unit="-", where="S14.3 completed cross-host pairs",
         record=_F36 + "/determinism/aggregate.json", grep=r'"both_completed": 28'),
    dict(id="Y55", printed=11, unit="-", where="S14.3 farm-internal bit-identical 12-cell pairs",
         record=_F36 + "/determinism/summary_stdout.txt", grep=r"farm-internal pairs 14 bit-identical 11"),
    dict(id="Y56", printed=9474, unit="-", where="S14.3 host statement, AMD EPYC 9474F, CPython 3.14.5 (server); CPython 3.14.7 (workstation)",
         record=_F36R, grep=r"AMD EPYC 9474F, CPython 3\.14\.5"),
    dict(id="Y57", printed=0.279, unit="K", where="S14.5 conductivity axis, peak temperature, composition only",
         record=_F36 + "/items_08_21/ITEMS_08_21.md", grep=r"moves 0\.279 K \(composition-only\)"),
    dict(id="Y58", printed=0.266, unit="K", where="S14.5 conductivity axis, peak temperature, temperature only; Sec. 6.3 prints 0.27 to 0.28",
         record=_F36 + "/items_08_21/ITEMS_08_21.md", grep=r"and 0\.266 K \(temperature-only\)"),
    dict(id="Y59", printed=0.94, unit="percent", where="S14.5 rectangle pulse water half-width on the F36 point set",
         record=_F36R, grep=r"reads 0\.94 percent on F36"),
    dict(id="Y60", printed=0.98, unit="percent", where="S14.5 committed record on the same two corners",
         record=_F36R, grep=r"the committed record gives 0\.98 percent"),
    dict(id="Y61", printed=7, word=True, unit="levels", where="S14.5 seven of twenty-one item-08 scenario levels complete",
         record=_F36R, grep=r"7 of 21 scenario levels complete"),
    # --- the mechanism (M-F1 and the build record, landed)
    dict(id="E1", printed=2.886e-11, unit="-", where="Sec. 6.1 (2.9e-11) / S3.4 polished residual, low",
         record=_SEED + "/polish/step2b_table.md", grep=r"2\.886e-11"),
    dict(id="E2", printed=8.172e-11, unit="-", where="Sec. 6.1 (8.2e-11) / S3.4 polished residual, high",
         record=_SEEDR, grep=r"2\.886e-11 to 8\.172e-11"),
    dict(id="E3", printed=1.67e-10, unit="-", where="Sec. 6.1 (1.7e-10) / S3.4 one lattice step on the water row, low",
         record=_SEEDR, grep=r"jumps 1\.67e-10 to 1\.73e-10"),
    dict(id="E4", printed=1.73e-10, unit="-", where="S3.4 one lattice step on the water row, high",
         record=_SEEDR, grep=r"1\.67e-10 to 1\.73e-10"),
    dict(id="E5", printed=1.88, unit="-", where="S14.3 inward root found only within 1.88 times the seed",
         record=_MF1 + "/M_F1_RESULT_head.md", grep=r"within 1\.88 times that flux"),
    dict(id="E6", printed=1.96, unit="-", where="S14.3 walks to the bound from 1.96 times",
         record=_MF1 + "/M_F1_RESULT_head.md", grep=r"from 1\.96"),
    dict(id="E7", printed=24, unit="-", where="Sec. 6.1 / S14.3 decisive arrival states with an inward root, 24 of 24",
         record="../docs/GT_PS2_FACE_DEPARTURE_M_F1_RECORD_2026-09-29.md", grep=r"24 of 24 decisive states"),
    dict(id="E8", printed=15, word=True, unit="-", where="Sec. 6.1 / S3.4 floor states above the literal, 15 of 24",
         record="../docs/GT_PS2_FACE_DEPARTURE_M_F1_RECORD_2026-09-29.md", grep=r"meets 2\.0e-11 on 9 states and sits at 2\.1e-11 to"),
    dict(id="E9", printed=4.732, unit="s", where="S2.6 first face arrival of the two-cell project run",
         record=_P01, grep=r"crossed z = 0\.125 at t = 4\.732 s"),
    dict(id="E10", printed=5.1e-13, unit="-", where="S2.6 macro ledgers of the two-cell project run",
         record=_P01, grep=r"macro ledgers <= 5\.1e-13"),
    dict(id="E11", printed=99.9999977, unit="percent", where="S2.6 core consumed by the dry continuation before it stopped",
         record=_P02, grep=r"99\.9999977 % consumed by volume"),
    # --- the comparator C1 (the weighted mean of Crank against the surface read), Sec. 5.2 and S15
    dict(id="M1", printed=1.411, unit="-", where="Sec. 5.2 (1.41) / S15 weighted-mean RMS, sunflower",
         record=_C1 + "/score.out", grep=r"CX  declared   faner2008_sunflower_136C n 60 t_max    t_end  2400\.0 rms_v2 1\.4112 band 0"),
    dict(id="M2", printed=0.979, unit="-", where="Sec. 5.2 (0.98) / S15 weighted-mean RMS, soybean",
         record=_C1 + "/score.out", grep=r"CX  declared   faner2008_soybean_136C   n 60 t_max    t_end  2400\.0 rms_v2 0\.9786 band 1"),
    dict(id="M3", printed=0.289, unit="-", where="Sec. 5.2 (0.29) / S15 weighted-mean drying time, sunflower",
         record=_C1 + "/score.out", grep=r"rms_v2 1\.4112 band 0 signed -1\.3465 \| dry 0\.289"),
    dict(id="M4", printed=0.302, unit="-", where="Sec. 5.2 (0.30) / S15 weighted-mean drying time, soybean",
         record=_C1 + "/score.out", grep=r"rms_v2 0\.9786 band 1 signed -0\.8934 \| dry 0\.302"),
    dict(id="M5", printed=75, unit="-", where="Sec. 5.2 / S15 weighted mean over the surface value, 5 to 75",
         record=_C1R, grep=r"It runs from about 5 to 75 across the window"),
    dict(id="M6", printed=45, unit="-", where="S15 time-weighted geometric mean 40 to 45",
         record=_C1R, grep=r"time-weighted geometric mean of 40 to 45"),
    dict(id="M7", printed=14.6, unit="-", where="S15 the surface read closer by a factor of 2.9 to 14.6 in RMS",
         record=_C1R, grep=r"by a factor of 2\.9 to 14\.6 in RMS"),
    dict(id="M8", printed=676, unit="-", where="S15 identity checks of the comparator",
         record=_C1 + "/score.out", grep=r"checks 676, mismatches 0"),
    dict(id="M9", printed=9.6e-15, unit="-", where="S15 largest stride ledger of the comparator marches",
         record=_C1R, grep=r"at most 9\.6e-15 of the initial inventory"),
    dict(id="M10", printed=0.634, unit="-", where="S15 pore-gas-weighted variant, sunflower",
         record=_C1R, grep=r"CC Crank in pore-gas c \| sunflower \| 0\.634"),
    dict(id="M11", printed=0.465, unit="-", where="S15 pore-gas-weighted variant, soybean",
         record=_C1R, grep=r"CC Crank in pore-gas c \| soybean \| 0\.465"),
    dict(id="M12", printed=1.410, page="1.410", unit="-", where="S15 current-profile variant, sunflower",
         record=_C1R, grep=r"CXp Crank, current profile \| sunflower \| 1\.410"),
    dict(id="M13", printed=0.977, unit="-", where="S15 current-profile variant, soybean",
         record=_C1R, grep=r"CXp Crank, current profile \| soybean \| 0\.977"),
    # --- the verified audit items (derived numbers, verify_N and verify_P)
    dict(id="M14", printed=10.20, page="10.20", unit="percent", where="S14.6 (Sec. 6.1 and Table 1: 10.2) ablation gap at t = 0",
         record=_VN + "/verify_numbers.txt", grep=r"oracle s/R=0\.5445 solver s/R=0\.6000 rel=10\.20%"),
    dict(id="M15", printed=0.0110, page="0.0110", unit="kg/kg", where="S14.6 initial shell loading if the core is at Xc",
         record=_VN + "/verify_numbers.txt", grep=r"shell loading if core at Xc 0\.01101 vs Xe 0\.02336"),
    dict(id="M16", printed=5.684e-14, unit="K", where="S14.4 the pin miss, one ULP of 343.48 K",
         record=_VN + "/verify_numbers.txt", grep=r"difference 5\.684e-14 K, ulp 5\.684e-14 K"),
    dict(id="M17", printed=168, unit="tests", where="S14.4 tests passing on the coupled path",
         record=_I03P, grep=r"1 failed, 168 passed"),
    dict(id="M18", printed=9.86, unit="-", where="S14.5 the law's 95 per cent band at a fixed state, 9.86 to one",
         record=_VN + "/verify_numbers.txt", grep=r"end-to-end band at fixed state 9\.86:1"),
    dict(id="M19", printed=1.2e-3, page="1.2e-3", unit="kg/kg", where="S12.4 loading at which the law returns 4.0e-10 at 343 K",
         record=_VN + "/verify_numbers.txt", grep=r"4\.0e-10 at 343 K: 1\.162e-03 kg/kg"),
    dict(id="M20", printed=0.880, page="0.880", unit="-", where="Sec. 5.1 duty x crossing product, journal soybean",
         record=_VP + "/verify_P.out", grep=r"journal soybean charge       0\.9076 x 0\.9695 = 0\.8799"),
    dict(id="M21", printed=0.920, page="0.920", unit="-", where="Sec. 5.1 duty x crossing product, journal sunflower",
         record=_VP + "/verify_P.out", grep=r"journal sunflower charge     0\.8152 x 1\.1291 = 0\.9204"),
    dict(id="M22", printed=2.653, unit="nm", where="S15 (Sec. 5.2: 0.5 to 2.7) Knudsen pore at 50 C, tau 10",
         record=_VP + "/verify_P.out", grep=r"T  50\.0 C tau 10\.0: d = 2\.653 nm"),
    dict(id="M23", printed=0.531, unit="nm", where="S15 Knudsen pore at 50 C, tau 2",
         record=_VP + "/verify_P.out", grep=r"T  50\.0 C tau  2\.0: d = 0\.531 nm"),
    dict(id="M24", printed=2.357, unit="nm", where="S15 Knudsen pore at 136 C, tau 10 (printed 2.36)",
         record=_VP + "/verify_P.out", grep=r"T 136\.0 C tau 10\.0: d = 2\.357 nm"),
    dict(id="M25", printed=184, unit="-", where="Sec. 5.2 / S15 gas diffusion over the law, 37 to 184",
         record=_VP + "/VERIFY_P.md", grep=r"1\.29e-7 to 6\.47e-7 m2/s, is 37 to 184 times the law"),
    dict(id="M26", printed=1.584, page="1.6", unit="kPa", where="Sec. 5.5 band top reached 1.6 kPa above the dome",
         record=_VP + "/verify_P.out", grep=r"reached at P = 100\.884 kPa, 1\.584 kPa above the dome"),
    dict(id="M27", printed=0.292, page="0.29", unit="K/kPa", where="Sec. 5.5 (0.29) T* per kPa",
         record=_VP + "/VERIFY_P.md", grep=r"1\.584 kPa above the dome, 0\.292 K per kPa"),
    dict(id="M28", printed=0.687, page="0.69", unit="K", where="Sec. 7.3 (0.69) T* shift at Pixton's adsorption activity",
         record=_VP + "/verify_P.out", grep=r"Pixton adsorption aw 0\.8170: 105\.153 C; shift 0\.687 K"),
    dict(id="M29", printed=0.071, unit="kg/kg", where="Sec. 7.2 (0.07) plateau end, soybean, Faner et al. (2019) Fig. 2",
         record=_VP + "/verify_P.out", grep=r"last sample within 0\.6 K 84\.9 s \(X 0\.071\)"),
    dict(id="M30", printed=0.097, unit="kg/kg", where="Sec. 7.2 (0.10) plateau end, sunflower",
         record=_VP + "/verify_P.out", grep=r"last sample within 0\.6 K 109\.0 s \(X 0\.097\)"),
    dict(id="M31", printed=75, unit="percent", where="Sec. 7.2 direct steam share of the heat, Kemper (2013) p. 111",
         record=_VP + "/kemper.txt", grep=r"provides approximately 75%"),
    dict(id="M32", printed=19.4, page="19", unit="-", where="S10.7 (19) weighted mean over the surface value, span 52",
         record=_VP + "/verify_P.out", grep=r"Dbar/D\(surface\) = 19\.4"),
    dict(id="M33", printed=37.2, page="37", unit="-", where="S10.7 (37) weighted mean over the surface value, span 101",
         record=_VP + "/verify_P.out", grep=r"Dbar/D\(surface\) = 37\.2"),
    dict(id="M34", printed=112.4, unit="ppm", where="abstract / Table 2 / Sec. 5.4 / Sec. 8 last-tray floor with the oil arm, high",
         record=_I34D, grep=r"With the arm the same\s+rows are 8\.1 to 112\.4"),
    dict(id="M35", printed=8.1, unit="ppm", where="Table 2 / Sec. 5.4 / Sec. 8 last-tray floor with the oil arm, low",
         record=_I34D, grep=r"8\.1 to 112\.4"),
    # --- wording withdrawn by the rewrite (absent from the main text and both supplements)
    dict(id="Z69", printed=0, unit="-", where="'a limit of the two-region model and not of any tolerance' (M-F1, F36)",
         expect="absent", grep=r"a limit of the two-region model and not of any\s+tolerance"),
    dict(id="Z70", printed=0, unit="-", where="'the obstacle lies in the two-region description rather than in a tolerance'",
         expect="absent", grep=r"obstacle lies in\s+the two-region\s+description rather than in a\s+tolerance"),
    dict(id="Z71", printed=0, unit="-", where="'governed by (its|the) region of lowest diffusivity' (audit A7, comparator C1)",
         expect="absent", grep=r"governed by\s+(its|the)\s+region\s+of\s+lowest\s+diffusivity"),
    dict(id="Z72", printed=0, unit="-", where="'the physically indicated' surface edge (audit A7)",
         expect="absent", grep=r"physically\s+indicated"),
    dict(id="Z73", printed=0, unit="-", where="'No adsorption-branch water isotherm exists for soybean meal at any temperature' (audit A8)",
         expect="absent", grep=r"No adsorption-branch water isotherm exists"),
    dict(id="Z74", printed=0, unit="-", where="'more than 98 percent of the duty' (audit A10; replaces C41)",
         expect="absent", grep=r"\\SI\{98\}\{\\percent\} of the duty"),
    dict(id="Z75", printed=0, unit="-", where="'168 assertions' (audit A5: 168 tests)",
         expect="absent", grep=r"168 assertions"),
    dict(id="Z76", printed=0, unit="-", where="the pin miss 'about 3.6e-13 K' (audit A5: 5.684e-14 K)",
         expect="absent", grep=r"3\.6e-13\}\{\\kelvin\}"),
    dict(id="Z77", printed=0, unit="-", where="'has not been executed' for the 360-trajectory campaign (audit B8)",
         expect="absent", grep=r"has not been\s+executed"),
    dict(id="Z78", printed=0, unit="-", where="'No other tolerance in this work carries an exception' (audit A4)",
         expect="absent", grep=r"No other\s+tolerance in this\s+work carries an\s+exception"),
    dict(id="Z79", printed=0, unit="-", where="'declared safety factors over measured floors' (audit A4)",
         expect="absent", grep=r"declared safety factors over measured\s+floors"),
    dict(id="Z80", printed=0, unit="-", where="'coupled trajectory through activation and extinction' (audit A6)",
         expect="absent", grep=r"through activation and extinction"),
    dict(id="Z81", printed=0, unit="-", where="'all twenty-five boolean acceptance and contract gates are true / hold' (18 of 25 true)",
         expect="absent", grep=r"twenty-five\s+boolean acceptance and contract gates (hold|are true)"),
    dict(id="Z82", printed=0, unit="-", where="'place the prior inside a measured band' (audit A3)",
         expect="absent", grep=r"place the prior inside a measured band"),
    dict(id="Z83", printed=0, unit="-", where="the abstract's bare 'Results compared with experiment are mesh-converged;' (audit B3)",
         expect="absent", grep=r"Results compared with experiment are mesh-converged;"),
    dict(id="Z84", printed=0, unit="-", where="'retained only as that limit and as a verification oracle' (audit B1)",
         expect="absent", grep=r"retained only as that limit and as a verification oracle"),
]

# C41 (the tower lane's 98 per cent) is withdrawn by audit item A10; its
# wording is checked absent by Z74, its record stays where it was.
for _a in ASSERTIONS:
    if _a["id"] == "C41":
        _a["where"] = _a["where"] + " [withdrawn from the page 2026-09-30, audit A10: see Z74]"
        _a["withdrawn_from_page"] = True
    # R6/R7 (item 08) and Q17 (item 21) are not re-earned on the F36 tree; the
    # page says so where it prints them (Sec. 6.3, S14.5, S7.2).
    if _a["id"] in ("R6", "R7"):
        _a["where"] = _a["where"] + " [2026-09-20 tree; not re-earned on F36, stated in Sec. 6.3 and S14.5]"
    if _a["id"] == "Q17":
        _a["where"] = _a["where"] + " [2026-09-20 tree; not re-earned on F36 (packaging defect), stated in S7.2]"
    # R1, R2, R15-R18: the 2026-09-20 values, still printed beside the F36 ones.
    if _a["id"] in ("R1", "R2", "R15", "R16", "R17", "R18", "N8", "N9"):
        _a["where"] = _a["where"] + " [2026-09-20 tree, printed beside its F36 value in S14.1 / S14.3]"

#: extra corpus roots for Pass 1 (the numbers of the rewrite live in evidence
#: outside paper/analysis); each is read only if it exists or is staged.
EXTRA_CORPUS_ROOTS = [_F36, _C1, _VN, _VP, _SEED, _MF1]
EXTRA_CORPUS_EXCLUDE = {"jobs", "extract", "runs", "queue"}

# ---------------------------------------------------------------------------
# 2026-09-30, the verification pass (REFEREE_AUDIT_7_VERIFY_2026-09-30.md):
# the held ablation window of item 09 (Sec. 6.1, S14.6), from the adverse-
# challenge repair and return-hold build record, section 4.3, and its
# evidence item09/; the t = 0 ablation gap on the item-09 measure (derived,
# the derivation's output printed in the verify file); five withdrawn
# wordings.  The build record and its evidence had not landed when this block
# was written: their assumed landed paths are the two constants below, each
# with a staging variable in PENDING_ROOTS, as in the R-P1 block above.
# ---------------------------------------------------------------------------
_ADV = "../docs/evidence/adverse_repair_return_hold_2026-09-30"   # scratchpad adverse_build/evidence/
_ADVR = "../docs/GT_PS2_ADVERSE_CHALLENGE_REPAIR_AND_RETURN_HOLD_BUILD_RECORD_2026-09-30.md"  # scratchpad adverse_build/docs/
_VER = "01_particle_jfpe/REFEREE_AUDIT_7_VERIFY_2026-09-30.md"
PENDING_ROOTS[_ADV] = "DTDC_PAPER_ADVERSE_EVIDENCE"
PENDING_ROOTS[_ADVR] = "DTDC_PAPER_ADVERSE_RECORD"
EXTRA_CORPUS_ROOTS.append(_ADV)

_H12 = _ADV + "/item09/holdB_n12_1800.csv"
ASSERTIONS += [
    dict(id="AH1", printed=2.025, unit="s", where="Sec. 6.1 / S14.6 12-cell hold at the declared bound accepted to 2.025 s",
         record=_H12, grep=r"return_hold,23,2\.0249999999999995,"),
    dict(id="AH2", printed=23.9, unit="percent", where="Sec. 6.1 / S14.6 excess over the oracle at 0.9 s (relative to the oracle)",
         record=_H12, grep=r"return_hold,8,0\.8999999999999998,[^\n]*,23\.90896505262101,"),
    dict(id="AH3", printed=163, unit="percent", where="Sec. 6.1 / S14.6 excess over the oracle at 1.65 s",
         record=_H12, grep=r"return_hold,18,1\.6499999999999995,[^\n]*,163\.0767123030493,"),
    dict(id="AH4", printed=1.725, unit="s", where="Sec. 6.1 / S14.6 first held step with the oracle front at the center",
         record=_H12, grep=r"return_hold,19,1\.7249999999999994,[^\n]*ZeroDivisionError"),
    dict(id="AH5", printed=0.42, page="0.42", unit="-", where="Sec. 6.1 / S14.6 resolved s/R at 1.725 s",
         record=_H12, grep=r"return_hold,19,1\.7249999999999994,337\.71683316903943,0\.4232437070107339,"),
    dict(id="AH6", printed=0.39, page="0.39", unit="-", where="Sec. 6.1 / S14.6 resolved s/R at 2.025 s",
         record=_H12, grep=r"return_hold,23,2\.0249999999999995,337\.709349389664,0\.3859358468704856,"),
    dict(id="AH7", printed=0.13, page="0.13", unit="points", where="Sec. 6.1 / S14.6 12 against 24 cells at 0.9 s (0.132)",
         record=_ADVR, grep=r"0\.036 to 0\.132\s+points over the hold"),
    dict(id="AH8", printed=0.05, page="0.05", unit="points", where="S14.6 12 against 24 cells over the three segments, 0.02 to 0.05",
         record=_ADVR, grep=r"0\.021 to 0\.051\s+points over the three segments"),
    dict(id="AH9", printed=4.51e-11, unit="-", where="S14.6 held face arrival resolved at the declared bound",
         record=_ADVR, grep=r"4\.513023587547526e-11"),
    dict(id="AH10", printed=5.44e-12, unit="-", where="S14.6 held face arrival resolved below the literal",
         record=_ADVR, grep=r"5\.439436022034718e-12"),
    dict(id="AH11", printed=2.93e-13, unit="-", where="S14.6 cumulative ledger of the 12-cell hold",
         record=_ADVR, grep=r"cumulative 2\.9271954109817904e-13, every hold flag true"),
    dict(id="AH12", printed=1.05, unit="s", where="S14.6 the literal hold stops at the next face arrival",
         record=_ADVR, grep=r"\*\*24 of 24 to 2\.025 s\*\* \(the literal stopped at 1\.05 s\)"),
    dict(id="AH13", printed=0.675, unit="s", where="S14.6 the 24-cell hold, 9 of 9 to 0.9 s",
         record=_ADVR, grep=r"0\.675 s \(9 steps, horizon 0\.9 s\)"),
    dict(id="AH14", printed=1.8, page="1.8", unit="s", where="S14.6 the 12-cell hold, 24 steps to 2.025 s",
         record=_ADVR, grep=r"1\.8 s \(24 steps, horizon 2\.025 s\)"),
    dict(id="AH15", printed=9.25, unit="percent", where="S14.6 the t = 0 gap on the item-09 measure (derived, 1 - 0.544482/0.6)",
         record=_VER, grep=r"gap/solver=9\.253%"),
    # --- wording withdrawn by the verification pass
    dict(id="Z85", printed=0, unit="-", where="'falls to 9.43' (the t = 0 gap was on another measure; it grows)",
         expect="absent", grep=r"falls to\s*\$9\.43\$"),
    dict(id="Z86", printed=0, unit="-", where="'schedule extension not made' (the held window is measured)",
         expect="absent", grep=r"schedule extension not made"),
    dict(id="Z87", printed=0, unit="-", where="'the same jobs differ by' (six comparisons change path)",
         expect="absent", grep=r"the same\s+jobs differ by"),
    dict(id="Z88", printed=0, unit="-", where="'do not identify it without identifying it' (E12.4 garble)",
         expect="absent", grep=r"identify it\s+without identifying it"),
    dict(id="Z89", printed=0, unit="-", where="item 04 'passes under the declared bound at' (F37a local: no candidate)",
         expect="absent", grep=r"passes under the declared bound at"),
]


# ---------------------------------------------------------------------------
# 2026-09-30, REFEREE_AUDIT_8 part 1 (REFEREE_AUDIT_8_2026-09-30.md): the
# regime map relabelled as the Fickian shell's and the formulation's own
# recession at one band state (the regime-map check); the new Sec. 5.3 and
# S9.15 (the two-regime derivation, its 37 entries TR1-TR37 as the derivation
# wrote them, re-labelled from its draft, and TR38-TR84); the coupled march's
# admissible limit and cost (the feasibility record); the held window at three
# meshes (F38a); the film axis's lower corner and the F37b axis; the 96-cell
# time chain and the pending items 04 and 10; the withdrawn wording Z90-Z108.
# Every root below has landed in the tree (staged or committed); the pending
# variables are kept in the pattern of the blocks above.
# ---------------------------------------------------------------------------
_TR = "../docs/evidence/paper1_two_regime_derivation_2026-09-30"
_TRR = "../docs/GT_PS2_PAPER1_TWO_REGIME_DERIVATION_2026-09-30.md"
_RM = "../docs/evidence/paper1_regime_map_check_2026-09-30"
_RMR = "../docs/GT_PS2_REGIME_MAP_CHECK_2026-09-30.md"
_FEAS = "../docs/GT_PS2_FANER_MARCH_FEASIBILITY_2026-09-30.md"
_F38A = "../docs/evidence/farm_f38a_2026-09-30"
_F38AR = "../docs/GT_PS2_F38A_COLLECTION_RECORD_2026-09-30.md"
_F37B = "../docs/evidence/farm_f37b_2026-09-30"
_B410 = "../docs/GT_PS2_ITEM04_ITEM10_COMPLETION_BUILD_RECORD_2026-09-30.md"
_S1 = "../docs/GT_PS2_ITEM04_STAGE1_OUTCOME_2026-09-30.md"
_I07 = "analysis/results_2026-09-19/item_07_regime_map/item07_regime_map.json"
_I27RUN = "analysis/results_2026-09-19/item_27_boundary_value_surface_read/runs/band_faner2008_soybean_136C_plus_one_sigma.json"
PENDING_ROOTS[_TR] = "DTDC_PAPER_TWO_REGIME_EVIDENCE"
PENDING_ROOTS[_TRR] = "DTDC_PAPER_TWO_REGIME_RECORD"
PENDING_ROOTS[_RM] = "DTDC_PAPER_REGIME_MAP_EVIDENCE"
EXTRA_CORPUS_ROOTS += [_TR, _RM, _F38A, _F37B]
ASSERTIONS += [
    dict(id="TR1", printed=0.179, unit="kg/kg", where="Sec. 5.3 / Table 4 mobile liquid X_c - W_eq(1, T_b)",
         record=_TR + "/two_regime.json", grep=r'"mobile_liquid": 0\.17886238534399065'),
    dict(id="TR2", printed=0.0205, unit="kg/kg", where="Sec. 5.3 / Table 4 W_eq(a_h = 1, T_b), eq. (gab)",
         record=_TR + "/two_regime.json", grep=r'"W_eq_a1_Tb_GAB": 0\.020487614656009356'),
    dict(id="TR3", printed=242.1, unit="s", where="Sec. 5.3 / Table 4 thesis soybean, crossing to the last measured loading 0.1006",
         record=_TR + "/two_regime.json", grep=r'"measured_s": 242\.08470213991097'),
    dict(id="TR4", printed=240.0, unit="s", where="Sec. 5.3 / Table 4 thesis sunflower, crossing to the last measured loading 0.1042",
         record=_TR + "/two_regime.json", grep=r'"measured_s": 239\.97917627119546'),
    dict(id="TR5", printed=25.4, unit="s", where="Sec. 5.3 / Table 4 journal soybean, crossing to 0.10",
         record=_TR + "/two_regime.json", grep=r'"measured_s": 25\.350441949123166'),
    dict(id="TR6", printed=42.5, unit="s", where="Sec. 5.3 / Table 4 journal sunflower, crossing to 0.10",
         record=_TR + "/two_regime.json", grep=r'"measured_s": 42\.50008976874932'),
    dict(id="TR7", printed=9.6, unit="s", where="Sec. 5.3 / Table 4 thesis soybean heat-limited to 0.1006 at r_c",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rc": 9\.639322528152599'),
    dict(id="TR8", printed=17.2, unit="s", where="Sec. 5.3 / Table 4 thesis soybean heat-limited to 0.1006 at r_pre",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rpre": 17\.227356634679467'),
    dict(id="TR9", printed=10.0, unit="s", where="Sec. 5.3 / Table 4 thesis sunflower heat-limited to 0.1042 at r_c",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rc": 9\.973015840038675'),
    dict(id="TR10", printed=18.8, unit="s", where="Sec. 5.3 / Table 4 thesis sunflower heat-limited to 0.1042 at r_pre",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rpre": 18\.792353640412102'),
    dict(id="TR11", printed=14.0, unit="s", where="Sec. 5.3 / Table 4 journal soybean heat-limited to 0.10 at r_c",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rc": 13\.996091624664562'),
    dict(id="TR12", printed=19.4, unit="s", where="Sec. 5.3 / Table 4 journal soybean heat-limited to 0.10 at r_pre",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rpre": 19\.43218357986262'),
    dict(id="TR13", printed=20.1, unit="s", where="Sec. 5.3 / Table 4 journal sunflower heat-limited to 0.10 at r_c",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rc": 20\.08126815886464'),
    dict(id="TR14", printed=32.6, unit="s", where="Sec. 5.3 / Table 4 journal sunflower heat-limited to 0.10 at r_pre",
         record=_TR + "/two_regime.json", grep=r'"heat_limited_s_rpre": 32\.57086398880304'),
    dict(id="TR15", printed=25, unit="-", where="Sec. 5.3 / Table 4 thesis soybean measured over heat-limited, at r_c (at least)",
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_heat_limited_rc": 25\.114285929626128'),
    dict(id="TR16", printed=13, unit="-", where="Sec. 5.3 / Table 4 thesis sunflower measured over heat-limited, at r_pre (at least)",
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_heat_limited_rpre": 12\.770043649834854'),
    dict(id="TR17", printed=1.3, unit="-", where="Sec. 5.3 / Table 4 journal soybean measured over heat-limited, at r_pre",
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_heat_limited_rpre": 1\.3045596160070028'),
    dict(id="TR18", printed=2.1, unit="-", where="Sec. 5.3 / Table 4 journal sunflower measured over heat-limited, at r_c",
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_heat_limited_rc": 2\.116404672878598'),
    dict(id="TR19", printed=0.80, unit="-", where="Sec. 5.3 / Table 4 thesis soybean measured over the Fick shell at 0.666 mm",
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_fick": 0\.7972259354606936'),
    dict(id="TR20", printed=0.73, unit="-", where="Sec. 5.3 / Table 4 thesis sunflower measured over the Fick shell at 0.7215 mm",
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_fick": 0\.7341519868141062'),
    dict(id="TR21", printed=541, unit="s", where="Sec. 5.3 / Table 4 journal soybean Fick shell (0.885 mm) crossing to 0.10",
         record=_TR + "/two_regime.json", grep=r'"fick_s": 541\.1079073315648'),
    dict(id="TR22", printed=528, unit="s", where="Sec. 5.3 / Table 4 journal sunflower Fick shell (0.885 mm) crossing to 0.10",
         record=_TR + "/two_regime.json", grep=r'"fick_s": 527\.5980457378759'),
    dict(id="TR23", printed=41.2, unit="s", where="Sec. 5.3 / Table 4 journal soybean plateau end (+1 K, interpolated) after the crossing",
         record=_TR + "/two_regime.json", grep=r'"plateau_end_measured_interp_plus1K_minus_tc_s": 41\.22426880936834'),
    dict(id="TR24", printed=45.6, unit="s", where="Sec. 5.3 / Table 4 journal sunflower plateau end (+1 K, interpolated) after the crossing",
         record=_TR + "/two_regime.json", grep=r'"plateau_end_measured_interp_plus1K_minus_tc_s": 45\.594488081137484'),
    dict(id="TR25", printed=22.7, unit="s", where="Sec. 5.3 / Table 4 journal soybean predicted plateau end, particle mean +1 K, r_c",
         record=_TR + "/two_regime.json", grep=r'22\.670559128192604'),
    dict(id="TR26", printed=33.0, unit="s", where="Sec. 5.3 / Table 4 journal soybean predicted plateau end, r_pre",
         record=_TR + "/two_regime.json", grep=r'32\.95521018918573'),
    dict(id="TR27", printed=34.2, unit="s", where="Sec. 5.3 / Table 4 journal sunflower predicted plateau end, r_c",
         record=_TR + "/two_regime.json", grep=r'34\.18243207652007'),
    dict(id="TR28", printed=58.0, page='58.0', unit="s", where="Sec. 5.3 / Table 4 journal sunflower predicted plateau end, r_pre",
         record=_TR + "/two_regime.json", grep=r'58\.008021436324135'),
    dict(id="TR29", printed=34.4, unit="s", where="Sec. 5.3 / Table 4 thesis soybean plateau left before the crossing (+1 K, interpolated; negative in the record)",
         record=_TR + "/two_regime.json", grep=r'"plateau_end_measured_interp_plus1K_minus_tc_s": \-34\.37507059280935'),
    dict(id="TR30", printed=18.4, unit="s", where="Sec. 5.3 / Table 4 thesis sunflower plateau left before the crossing (negative in the record)",
         record=_TR + "/two_regime.json", grep=r'"plateau_end_measured_interp_plus1K_minus_tc_s": \-18\.40639242020449'),
    dict(id="TR31", printed=11.0, unit="K", where="Sec. 5.3 / Table 4 thesis soybean Tp at the crossing above T_b",
         record=_TR + "/two_regime.json", grep=r'"Tp_at_crossing_minus_Tb_K": 10\.960188367259633'),
    dict(id="TR32", printed=9.8, unit="K", where="Sec. 5.3 / Table 4 thesis sunflower Tp at the crossing above T_b",
         record=_TR + "/two_regime.json", grep=r'"Tp_at_crossing_minus_Tb_K": 9\.767556735813159'),
    dict(id="TR33", printed=0.4, unit="K", where="Sec. 5.3 / Table 4 journal soybean Tp at the crossing, below T_b (negative in the record)",
         record=_TR + "/two_regime.json", grep=r'"Tp_at_crossing_minus_Tb_K": \-0\.37659958231071755'),
    dict(id="TR34", printed=0.1, unit="K", where="Sec. 5.3 / Table 4 journal sunflower Tp at the crossing, below T_b (negative in the record)",
         record=_TR + "/two_regime.json", grep=r'"Tp_at_crossing_minus_Tb_K": \-0\.12974693501004708'),
    dict(id="TR35", printed=3.2, unit="s", where="Sec. 5.3 / Table 4 thesis soybean Fick-shell lumped temperature +1 K after its crossing",
         record=_TR + "/two_regime.json", grep=r'"plateau_end_pred_fick_plus1K_s": 3\.2446130218057703'),
    dict(id="TR36", printed=7.2, unit="s", where="Sec. 5.3 / Table 4 journal sunflower Fick-shell lumped temperature +1 K after its crossing",
         record=_TR + "/two_regime.json", grep=r'"plateau_end_pred_fick_plus1K_s": 7\.16123768119828'),
    dict(id="TR37", printed=0.43, unit="kg/kg", where="Sec. 5.3 / Table 4 thesis soybean X when the thermocouple is 1 K above its plateau",
         record=_TR + "/two_regime.json", grep=r'"X_at_above_plateau_plus_1\.0K": 0\.4261972766132316'),
    dict(id="TR38", printed=0.53, page="0.53", unit="-", where='Sec. 5.3 rate over the 20 s before the crossing over r_c, low (thesis sunflower)',
         record=_TR + "/derived_checks.json", grep='"rate_pre_over_const":\\ 0\\.5250206059745008'),
    dict(id="TR39", printed=0.72, page="0.72", unit="-", where='Sec. 5.3 rate over the 20 s before the crossing over r_c, high (journal soybean)',
         record=_TR + "/derived_checks.json", grep='"rate_pre_over_const":\\ 0\\.7170100660124655'),
    dict(id="TR40", printed=80, unit="percent", where='Sec. 5.3 mobile liquid gone when the mean temperature is 1 K above T_b, low (thesis sunflower, r_c), percent',
         record=_TR + "/derived_checks.json", grep='"mobile_fraction_gone_at_meanT_plus1K_rc":\\ 0\\.7963075856135521'),
    dict(id="TR41", printed=96, unit="percent", where='Sec. 5.3 the same, high (journal sunflower, r_pre), percent',
         record=_TR + "/derived_checks.json", grep='"mobile_fraction_gone_at_meanT_plus1K_rpre":\\ 0\\.9559320074070002'),
    dict(id="TR42", printed=3.9, unit="-", where='Sec. 5.3 thesis soybean Tp(t_c) - T_b in reading half-widths',
         record=_TR + "/derived_checks.json", grep='"Tp_tc_minus_Tb_over_max_halfwidth_near_tc":\\ 3\\.947924445826702'),
    dict(id="TR43", printed=4.4, unit="-", where='Sec. 5.3 / S9.15 thesis sunflower Tp(t_c) - T_b in reading half-widths',
         record=_TR + "/derived_checks.json", grep='"Tp_tc_minus_Tb_over_max_halfwidth_near_tc":\\ 4\\.401426083424129'),
    dict(id="TR44", printed=1.6, page="1.6", unit="percent", where='Sec. 5.3 shell term over film term, low (journal soybean), percent',
         record=_TR + "/derived_checks.json", grep='"shell_over_film_closed_form":\\ 0\\.01637743930353334'),
    dict(id="TR45", printed=2.1, page="2.1", unit="percent", where='Sec. 5.3 shell term over film term, high (thesis sunflower), percent',
         record=_TR + "/derived_checks.json", grep='"shell_over_film_closed_form":\\ 0\\.020511872203710466'),
    dict(id="TR46", printed=66.43, page="2.3", unit="C", where='Sec. 5.3 thesis soybean early plateau (2.3 K below T_b = 68.71 C, derived)',
         record=_TR + "/two_regime.json", grep='"plateau_mean_t_le_30s_C":\\ 66\\.43357333333333'),
    dict(id="TR47", printed=1.257, unit="-", where='Sec. 5.3 start activity of the dry-shell storage in the marches of Sec. 5.2',
         record=_I27RUN, grep='"tail_start_activity":\\s+1\\.257002957166127'),
    dict(id="TR48", printed=303.7, unit="s", where='Table 4 dry-shell limit, thesis soybean, to the last measured loading',
         record=_TR + "/two_regime.json", grep='"fick_s":\\ 303\\.65883919721864'),
    dict(id="TR49", printed=326.9, unit="s", where='Table 4 dry-shell limit, thesis sunflower',
         record=_TR + "/two_regime.json", grep='"fick_s":\\ 326\\.87942085752377'),
    dict(id="TR50", printed=0.047, page="0.047", unit="-", where='Table 4 journal soybean measured over the dry-shell limit',
         record=_TR + "/two_regime.json", grep='"ratio_meas_over_fick":\\ 0\\.04684914340678751'),
    dict(id="TR51", printed=0.081, page="0.081", unit="-", where='Table 4 journal sunflower measured over the dry-shell limit',
         record=_TR + "/two_regime.json", grep='"ratio_meas_over_fick":\\ 0\\.0805539181050425'),
    dict(id="TR52", printed=14.5, unit="s", where='Table 4 thesis soybean predicted plateau end, r_c',
         record=_TR + "/two_regime.json", grep='"t_s":\\ 14\\.547259287806787'),
    dict(id="TR53", printed=28.8, unit="s", where='Table 4 thesis soybean predicted plateau end, r_pre',
         record=_TR + "/two_regime.json", grep='"t_s":\\ 28\\.776436826542042'),
    dict(id="TR54", printed=15.3, unit="s", where='Table 4 thesis sunflower predicted plateau end, r_c',
         record=_TR + "/two_regime.json", grep='"t_s":\\ 15\\.343609641091765'),
    dict(id="TR55", printed=32.4, unit="s", where='Table 4 thesis sunflower predicted plateau end, r_pre',
         record=_TR + "/two_regime.json", grep='"t_s":\\ 32\\.430937481959994'),
    dict(id="TR56", printed=3.4, page="3.4", unit="s", where='Table 4 thesis sunflower dry-shell lumped temperature +1 K',
         record=_TR + "/two_regime.json", grep='"plateau_end_pred_fick_plus1K_s":\\ 3\\.3575018591186883'),
    dict(id="TR57", printed=5.2, page="5.2", unit="s", where='Table 4 journal soybean dry-shell lumped temperature +1 K',
         record=_TR + "/two_regime.json", grep='"plateau_end_pred_fick_plus1K_s":\\ 5\\.183139246190109'),
    dict(id="TR58", printed=0.1006, unit="-", where='Table 4 thesis soybean last measured loading',
         record=_TR + "/two_regime.json", grep='"X":\\ 0\\.10059'),
    dict(id="TR59", printed=0.1042, unit="-", where='Table 4 thesis sunflower last measured loading',
         record=_TR + "/two_regime.json", grep='"X":\\ 0\\.10418'),
    dict(id="TR60", printed=0.3, page="0.30", unit="-", where='Sec. 5.3 thesis sunflower loading when the thermocouple is 1 K above its plateau',
         record=_TR + "/two_regime.json", grep='"X_at_above_plateau_plus_1\\.0K":\\ 0\\.30204285570011524'),
    dict(id="TR61", printed=17.3, unit="s", where='S9.15 core extinction without sensible heat, thesis soybean',
         record=_TR + "/two_regime.json", grep='"core_extinction_pred_no_sensible_s":\\ 17\\.27869686130633'),
    dict(id="TR62", printed=18.6, unit="s", where='S9.15 the same, thesis sunflower',
         record=_TR + "/two_regime.json", grep='"core_extinction_pred_no_sensible_s":\\ 18\\.56952708742658'),
    dict(id="TR63", printed=25.0, page="25.0", unit="s", where='S9.15 the same, journal soybean',
         record=_TR + "/two_regime.json", grep='"core_extinction_pred_no_sensible_s":\\ 25\\.03809755985096'),
    dict(id="TR64", printed=36.1, unit="s", where='S9.15 the same, journal sunflower',
         record=_TR + "/two_regime.json", grep='"core_extinction_pred_no_sensible_s":\\ 36\\.12287794871382'),
    dict(id="TR65", printed=22.6, unit="s", where='Table S (criterion) thesis_soybean_136C t_above_plateau_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_plateau_plus_5\\.0K_minus_t_cross":\\ \\-22\\.605046611762944'),
    dict(id="TR66", printed=14.4, unit="s", where='Table S (criterion) thesis_soybean_136C t_above_Tb_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_Tb_plus_5\\.0K_minus_t_cross":\\ \\-14\\.397634152386878'),
    dict(id="TR67", printed=8.9, unit="s", where='Table S (criterion) thesis_sunflower_136C t_above_plateau_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_plateau_plus_5\\.0K_minus_t_cross":\\ \\-8\\.946132685759956'),
    dict(id="TR68", printed=8.4, unit="s", where='Table S (criterion) thesis_sunflower_136C t_above_Tb_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_Tb_plus_5\\.0K_minus_t_cross":\\ \\-8\\.404795872920992'),
    dict(id="TR69", printed=39.2, unit="s", where='Table S (criterion) journal_soybean_120C t_above_Tb_plus_1.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_Tb_plus_1\\.0K_minus_t_cross":\\ 39\\.166936780709705'),
    dict(id="TR70", printed=54.5, unit="s", where='Table S (criterion) journal_soybean_120C t_above_plateau_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_plateau_plus_5\\.0K_minus_t_cross":\\ 54\\.53793007300497'),
    dict(id="TR71", printed=63.7, unit="s", where='Table S (criterion) journal_sunflower_100C t_above_plateau_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_plateau_plus_5\\.0K_minus_t_cross":\\ 63\\.66819156452938'),
    dict(id="TR72", printed=64.7, unit="s", where='Table S (criterion) journal_sunflower_100C t_above_Tb_plus_5.0K_minus_t_cross',
         record=_TR + "/two_regime.json", grep='"t_above_Tb_plus_5\\.0K_minus_t_cross":\\ 64\\.65835568772611'),
    dict(id="TR73", printed=20.0, page="20.0", unit="-", where='S9.15 2019 Fig. 3/4 soybean 136 C t_cross_to_0.10',
         record=_TR + "/two_regime.json", grep='"t_cross_to_0\\.10":\\ 20\\.047185856855943'),
    dict(id="TR74", printed=24.0, page="24.0", unit="-", where='S9.15 2019 Fig. 3/4 soybean 136 C t_cross_to_0.10',
         record=_TR + "/two_regime.json", grep='"t_cross_to_0\\.10":\\ 24\\.00069063157309'),
    dict(id="TR75", printed=8.5, page="8.5", unit="-", where='S9.15 2019 Fig. 3/4 soybean 136 C film_only_heat_limited_s',
         record=_TR + "/two_regime.json", grep='"film_only_heat_limited_s":\\ 8\\.459663973323071'),
    dict(id="TR76", printed=9.9, page="9.9", unit="-", where='S9.15 2019 Fig. 3/4 soybean 136 C film_only_heat_limited_s',
         record=_TR + "/two_regime.json", grep='"film_only_heat_limited_s":\\ 9\\.925365908568283'),
    dict(id="TR77", printed=2.37, page="2.37", unit="-", where='S9.15 2019 Fig. 3/4 soybean 136 C ratio',
         record=_TR + "/two_regime.json", grep='"ratio":\\ 2\\.3697378430246485'),
    dict(id="TR78", printed=2.68, page="2.68", unit="-", where='S9.15 2019 Fig. 3/4 soybean 136 C ratio',
         record=_TR + "/two_regime.json", grep='"ratio":\\ 2\\.6838340633701296'),
    dict(id="TR79", printed=18.5, page="18.5", unit="-", where='S9.15 2019 Fig. 3/4 sunflower 136 C t_cross_to_0.10',
         record=_TR + "/two_regime.json", grep='"t_cross_to_0\\.10":\\ 18\\.473238129954048'),
    dict(id="TR80", printed=22.7, page="22.7", unit="-", where='S9.15 2019 Fig. 3/4 sunflower 136 C t_cross_to_0.10',
         record=_TR + "/two_regime.json", grep='"t_cross_to_0\\.10":\\ 22\\.689897826837004'),
    dict(id="TR81", printed=9.4, page="9.4", unit="-", where='S9.15 2019 Fig. 3/4 sunflower 136 C film_only_heat_limited_s',
         record=_TR + "/two_regime.json", grep='"film_only_heat_limited_s":\\ 9\\.435541385736888'),
    dict(id="TR82", printed=11.0, page="11.0", unit="-", where='S9.15 2019 Fig. 3/4 sunflower 136 C film_only_heat_limited_s',
         record=_TR + "/two_regime.json", grep='"film_only_heat_limited_s":\\ 10\\.971559750856846'),
    dict(id="TR83", printed=1.8, page="1.80", unit="-", where='S9.15 2019 Fig. 3/4 sunflower 136 C ratio',
         record=_TR + "/two_regime.json", grep='"ratio":\\ 1\\.8044789808567328'),
    dict(id="TR84", printed=2.31, page="2.31", unit="-", where='S9.15 2019 Fig. 3/4 sunflower 136 C ratio',
         record=_TR + "/two_regime.json", grep='"ratio":\\ 2\\.3085375839366726'),
    dict(id="RM1", printed=7800.0, page='7.80', unit="-", where="Sec. 6.2 / S11.1 Lambda minimum, Fickian shell (the check's recomputation)",
         record=_RM + "/lambda_points.json", grep='"lambda_min":\\ 7804\\.738606368701'),
    dict(id="RM2", printed=62900.0, unit="-", where='Sec. 6.2 / S11.1 Lambda maximum, Fickian shell',
         record=_RM + "/lambda_points.json", grep='"lambda_max":\\ 62914\\.76131270393'),
    dict(id="RM3", printed=7800.0, page='7.80', unit="-", where='Sec. 6.2 / S11.1 Lambda minimum in the committed item-07 record',
         record=_I07, grep='"lambda_min":\\ 7804\\.738606368701'),
    dict(id="RM4", printed=62900.0, unit="-", where='Sec. 6.2 / S11.1 Lambda maximum in the committed item-07 record',
         record=_I07, grep='"lambda_max":\\ 62914\\.76131270393'),
    dict(id="RM5", printed=3.12e-06, unit="-", where='S11.1 D_eff for Lambda = 1, low (Fickian shell)',
         record=_I07, grep='3\\.1218954425474805e\\-06'),
    dict(id="RM6", printed=2.52e-05, unit="-", where='S11.1 D_eff for Lambda = 1, high',
         record=_I07, grep='2\\.5165904525081572e\\-05'),
    dict(id="RM7", printed=3.81e-06, unit="-", where='S11.1 k_dry for Lambda = 1, low',
         record=_I07, grep='3\\.814685059474882e\\-06'),
    dict(id="RM8", printed=3.08e-05, unit="-", where='S11.1 k_dry for Lambda = 1, high',
         record=_I07, grep='3\\.075054939113002e\\-05'),
    dict(id="RM9", printed=14200.0, unit="-", where='S11.1 pressure-feasible Lambda, low',
         record=_RM + "/lambda_points.json", grep='14171\\.518666644564'),
    dict(id="RM10", printed=3390000.0, unit="-", where='S11.1 pressure-feasible Lambda, high',
         record=_RM + "/lambda_points.json", grep='3391276\\.592736881'),
    dict(id="RM11", printed=400, page="400", unit="-", where='S11.1 Fick supply overstated by up to 400',
         record=_RM + "/lambda_points.json", grep='400\\.39517427628897'),
    dict(id="RM12", printed=306, unit="-", where='S11.1 states with c_sat(T_g) above the pressure ceiling',
         record=_RM + "/lambda_points.json", grep='"points_with_T_g_above_hexane_Tb_at_P\\ \\(c_sat\\(T_g\\)\\ above\\ the\\ pressure\\ ceiling\\)":\\ 306'),
    dict(id="RM13", printed=260.0, page="261", unit="-", where="Sec. 6.2 / S11.1 march flux over the Fick map's capacity, low (261)",
         record=_RM + "/march_compare.json", grep='260\\.7400327678003'),
    dict(id="RM14", printed=1800.0, page="1831", unit="-", where='Sec. 6.2 / S11.1 the same, high (1831)',
         record=_RM + "/march_compare.json", grep='1831\\.0572985728527'),
    dict(id="RM15", printed=1090.0, unit="-", where='Sec. 6.2 / S11.1 shell Peclet number, low',
         record=_RM + "/march_compare.json", grep='1088\\.0743351799933'),
    dict(id="RM16", printed=2930.0, unit="-", where='Sec. 6.2 / S11.1 shell Peclet number, high',
         record=_RM + "/march_compare.json", grep='2927\\.0868933431248'),
    dict(id="RM17", printed=0.653, unit="-", where='Sec. 6.2 / S11.1 loading removed at film 0.6264 over film 1, low (0.65)',
         record=_RM + "/march_compare.json", grep='0\\.6525984223690261'),
    dict(id="RM18", printed=0.698, unit="-", where='Sec. 6.2 / S11.1 the same, high (0.70)',
         record=_RM + "/march_compare.json", grep='0\\.6977145907280035'),
    dict(id="RM19", printed=21000.0, page='2.1\\times10^{4}', unit="-", where="S11.1 flux over a Fickian shell on the march's own states, low",
         record=_RM + "/march_compare.json", grep='21062\\.531316308978'),
    dict(id="RM20", printed=38000.0, page='3.8\\times10^{4}', unit="-", where='S11.1 the same, high',
         record=_RM + "/march_compare.json", grep='38045\\.35753491091'),
    dict(id="RM21", printed=0.22, page="0.22", unit="-", where='S11.1 film share of the heat resistance, low',
         record=_RM + "/march_compare.json", grep='0\\.2193178502540356'),
    dict(id="RM22", printed=0.96, page="0.96", unit="-", where='S11.1 film share of the heat resistance, high',
         record=_RM + "/march_compare.json", grep='0\\.9573277630804378'),
    dict(id="RM23", printed=337.36, unit="K", where='S11.1 interface temperature, low',
         record=_RM + "/march_compare.json", grep='337\\.35533197217813'),
    dict(id="RM24", printed=337.78, unit="K", where='S11.1 interface temperature, high',
         record=_RM + "/march_compare.json", grep='337\\.7777895732413'),
    dict(id="RM25", printed=336.78, unit="K", where="S11.1 the map's interface temperature at the marched state",
         record=_RM + "/march_compare.json", grep='"item07_T_eq_K":\\ 336\\.7789926584653'),
    dict(id="RM26", printed=0.025, page="0.025", unit="-", where="S11.1 flux over the map's conductive capacity, low",
         record=_RM + "/march_compare.json", grep='0\\.024789675485870397'),
    dict(id="RM27", printed=0.17, page="0.17", unit="-", where="S11.1 flux over the map's conductive capacity, high",
         record=_RM + "/march_compare.json", grep='0\\.1740864866273847'),
    dict(id="RM28", printed=0.01064, unit="-", where='Table S (march at the band state) s/R~0.85 hexane_flux_out_kg_m2_s',
         record=_RM + "/march_compare.json", grep='"hexane_flux_out_kg_m2_s":\\ 0\\.01063704557776373'),
    dict(id="RM29", printed=0.0069, page='6.90\\times10^{-3}', unit="-", where='Table S (march at the band state) s/R~0.55 hexane_flux_out_kg_m2_s',
         record=_RM + "/march_compare.json", grep='"hexane_flux_out_kg_m2_s":\\ 0\\.006899574637263331'),
    dict(id="RM30", printed=0.00332, unit="-", where='Table S (march at the band state) s/R~0.3 hexane_flux_out_kg_m2_s',
         record=_RM + "/march_compare.json", grep='"hexane_flux_out_kg_m2_s":\\ 0\\.003316378925695945'),
    dict(id="RM31", printed=346, unit="-", where='Table S (march at the band state) s/R~0.85 ratio_flux_over_fick_map',
         record=_RM + "/march_compare.json", grep='"ratio_flux_over_fick_map":\\ 346\\.4014752396197'),
    dict(id="RM32", printed=1042, unit="-", where='Table S (march at the band state) s/R~0.55 ratio_flux_over_fick_map',
         record=_RM + "/march_compare.json", grep='"ratio_flux_over_fick_map":\\ 1041\\.7028606638435'),
    dict(id="RM33", printed=1382, unit="-", where='Table S (march at the band state) s/R~0.3 ratio_flux_over_fick_map',
         record=_RM + "/march_compare.json", grep='"ratio_flux_over_fick_map":\\ 1382\\.1197476282246'),
    dict(id="RM34", printed=1380, unit="-", where='Table S (march at the band state) s/R~0.85 shell_peclet_N_L_over_cDb',
         record=_RM + "/march_compare.json", grep='"shell_peclet_N_L_over_cDb":\\ 1379\\.7434154303405'),
    dict(id="RM35", printed=2921, unit="-", where='Table S (march at the band state) s/R~0.55 shell_peclet_N_L_over_cDb',
         record=_RM + "/march_compare.json", grep='"shell_peclet_N_L_over_cDb":\\ 2921\\.0252406466866'),
    dict(id="RM36", printed=2404, unit="-", where='Table S (march at the band state) s/R~0.3 shell_peclet_N_L_over_cDb',
         record=_RM + "/march_compare.json", grep='"shell_peclet_N_L_over_cDb":\\ 2404\\.0665908755554'),
    dict(id="RM37", printed=0.936, unit="-", where='Table S (march at the band state) s/R~0.85 film_share_of_heat_resistance',
         record=_RM + "/march_compare.json", grep='"film_share_of_heat_resistance":\\ 0\\.9357822950763027'),
    dict(id="RM38", printed=0.759, unit="-", where='Table S (march at the band state) s/R~0.55 film_share_of_heat_resistance',
         record=_RM + "/march_compare.json", grep='"film_share_of_heat_resistance":\\ 0\\.7586344247313951'),
    dict(id="RM39", printed=0.532, unit="-", where='Table S (march at the band state) s/R~0.3 film_share_of_heat_resistance',
         record=_RM + "/march_compare.json", grep='"film_share_of_heat_resistance":\\ 0\\.5324193102385351'),
    dict(id="RM40", printed=0.848, unit="-", where='Table S (march at the band state) s/R~0.85 s_over_R',
         record=_RM + "/march_compare.json", grep='"s_over_R":\\ 0\\.8478904024251663'),
    dict(id="RM41", printed=0.546, unit="-", where='Table S (march at the band state) s/R~0.55 s_over_R',
         record=_RM + "/march_compare.json", grep='"s_over_R":\\ 0\\.5459334955628685'),
    dict(id="RM42", printed=0.303, unit="-", where='Table S (march at the band state) s/R~0.3 s_over_R',
         record=_RM + "/march_compare.json", grep='"s_over_R":\\ 0\\.3034137153616382'),
    dict(id="RM43", printed=337.72, unit="-", where='Table S (march at the band state) s/R~0.85 T_interface_K',
         record=_RM + "/march_compare.json", grep='"T_interface_K":\\ 337\\.7170786205579'),
    dict(id="RM44", printed=0.282, unit="-", where='S4.6 spherical thermal crossing at h 106.5',
         record=_RM + "/resistances.json", grep='"thermal_crossing_s_over_R_sphere_h=106\\.5":\\ 0\\.28197994031458257'),
    dict(id="RM45", printed=0.277, unit="-", where='S4.6 spherical thermal crossing at h 103.74',
         record=_RM + "/resistances.json", grep='"thermal_crossing_s_over_R_sphere_h=103\\.73636577329322":\\ 0\\.2766872646537152'),
    dict(id="RM46", printed=2.25, unit="mm", where='S4.6 slab thermal crossing shell, mm',
         record=_RM + "/resistances.json", grep='"thermal_slab_crossing_shell_mm_h=106\\.5":\\ 2\\.2535211267605635'),
    dict(id="RM47", printed=0.039, page="0.039", unit="-", where='Table S resistances, thermal slab at 0.90 R',
         record=_RM + "/resistances.json", grep='"thermal_slab_h_L_over_k":\\ 0\\.03927187499999999'),
    dict(id="RM48", printed=0.196, unit="-", where='Table S resistances, thermal slab at 0.50 R',
         record=_RM + "/resistances.json", grep='"thermal_slab_h_L_over_k":\\ 0\\.196359375'),
    dict(id="RM49", printed=0.044, page="0.044", unit="-", where='Table S resistances, thermal sphere at 0.90 R',
         record=_RM + "/resistances.json", grep='"thermal_sphere_h_R_L_over_k_s":\\ 0\\.04363541666666666'),
    dict(id="RM50", printed=0.393, unit="-", where='Table S resistances, thermal sphere at 0.50 R',
         record=_RM + "/resistances.json", grep='"thermal_sphere_h_R_L_over_k_s":\\ 0\\.39271875'),
    dict(id="RM51", printed=0.275, unit="-", where='Table S resistances, thermal slab at 0.30 R',
         record=_RM + "/resistances.json", grep='"thermal_slab_h_L_over_k":\\ 0\\.274903125'),
    dict(id="RM52", printed=0.916, unit="-", where='Table S resistances, thermal sphere at 0.30 R',
         record=_RM + "/resistances.json", grep='"thermal_sphere_h_R_L_over_k_s":\\ 0\\.9163437500000002'),
    dict(id="RM53", printed=170000.0, page='1.7\\times10^{5}', unit="-", where='Table S resistances, Fickian mass slab at 0.30 R',
         record=_RM + "/resistances.json", grep='"mass_slab_km_L_over_D":\\ 168813\\.75000000003'),
    dict(id="RM54", printed=1.0003, unit="-", where='Sec. 6.3 / S14.5 hexane removed, factor 3 in D_b (composition only)',
         record=_RM + "/item08_front_travel.json", grep='"hexane_out_ratio":\\ 1\\.0002868727712702'),
    dict(id="RM55", printed=1.171, unit="-", where='Sec. 6.3 / S14.5 hexane removed, film 0.5 to 1.0 (correlated)',
         record=_RM + "/item08_front_travel.json", grep='"hexane_out_ratio_hi_over_lo":\\ 1\\.1709085855485382'),
    dict(id="RM56", printed=7.05, unit="-", where='Sec. 6.2 / S11.1 R5 reaches 7.05 s',
         record=_FEAS, grep='\\|\\s+R5_front090_y085_dt0075\\s+\\|\\s+N=12,\\s+dt\\s+0\\.075,\\s+y\\s+0\\.85,\\s+film\\s+1\\.0,\\s+front\\s+0\\.9\\s+\\|\\s+refused\\s+at\\s+7\\.05\\s+s'),
    dict(id="RM57", printed=0.066, unit="-", where='S11.1 R5 front at the end',
         record=_FEAS, grep='\\|\\s+93\\s+\\|\\s+7\\.05\\s+\\|\\s+0\\.900\\s+\\->\\s+0\\.066\\s+\\|'),
    dict(id="RM58", printed=7.2e-11, unit="-", where='S11.1 R5 largest scaled residual',
         record=_FEAS, grep='\\|\\s+16\\.3\\s+\\|\\s+7\\.22e\\-11\\s+\\|\\s+1\\.2e\\-12\\s+/\\s+1\\.8e\\-12\\s+\\|'),
    dict(id="RM59", printed=1.8e-12, unit="-", where='S11.1 R5 cumulative ledger',
         record=_FEAS, grep='\\|\\s+16\\.3\\s+\\|\\s+7\\.22e\\-11\\s+\\|\\s+1\\.2e\\-12\\s+/\\s+1\\.8e\\-12\\s+\\|'),
    dict(id="RM60", printed=2.55, unit="-", where='S11.1 R9 reaches 2.55 s',
         record=_FEAS, grep='\\|\\s+R9_front090_y085_dt0075_film06264\\s+\\|\\s+N=12,\\s+dt\\s+0\\.075,\\s+y\\s+0\\.85,\\s+film\\s+0\\.6264,\\s+front\\s+0\\.9\\s+\\|\\s+refused\\s+at\\s+2\\.55\\s+s'),
    dict(id="RM61", printed=0.15, page="0.150", unit="-", where='Sec. 6.2 / S11.1 the start loading',
         record=_FEAS, grep='\\|\\s+start\\s+loading\\s+\\|\\s+0\\.1499\\s+kg/kg,\\s+front\\s+0\\.90\\s+R\\s+\\(R5\\)'),
    dict(id="FE1", printed=0.857, unit="-", where='S9.15 gas-only chart y_h at most 0.857 at 409.15 K',
         record=_FEAS, grep='admits\\s+at\\s+most\\s+y_h\\s+=\\s+0\\.857\\s+at\\s+136\\s+C\\s+\\(at\\s+least\\s+14\\.3\\s+mol%\\s+water'),
    dict(id="FE2", printed=0.143, unit="-", where='S9.15 water-vapour mole fraction at least 0.143',
         record=_FEAS, grep='\\(at\\s+least\\s+14\\.3\\s+mol%\\s+water'),
    dict(id="FE3", printed=0.047, unit="-", where='S9.15 isotherm qualified down to a_w 0.047',
         record=_FEAS, grep='\\(a_w\\s+>=\\s+0\\.047\\)'),
    dict(id="FE4", printed=0.023, unit="-", where='S9.15 retained water 0.023 kg/kg',
         record=_FEAS, grep='qualified\\s+only\\s+down\\s+to\\s+0\\.023\\s+kg/kg\\s+retained\\s+water'),
    dict(id="FE5", printed=54, page="54", unit="-", where='S9.15 wall per physical second between face events, low',
         record=_FEAS, grep='Between\\s+face\\s+events\\s+54\\s+to\\s+84\\s+s\\s+of\\s+wall\\s+per\\s+physical'),
    dict(id="FE6", printed=119, page="119", unit="-", where='S9.15 wall per physical second between face events, high (N = 48)',
         record=_FEAS, grep='119\\s+at\\s+N\\s+=\\s+48\\s+\\(R10\\)'),
    dict(id="FE7", printed=110, page="110", unit="-", where='S9.15 each face crossing about 50 to 110 s',
         record=_FEAS, grep='each\\s+master\\-face\\s+crossing\\s+costs\\s+about\\s+50\\s+to\\s+110\\s+s'),
    dict(id="HA1", printed=11.86, unit="-", where='S14.6 Table S (held window, server) N = 12 excess at 0.225 s',
         record=_F38A + "/science_table.csv", grep='f38an12h1800comp,12,1\\.8,composition_only,24/24,2\\.0249999999999995,11\\.858604953566028,'),
    dict(id="HA2", printed=11.82, unit="-", where='S14.6 Table S (held window, server) N = 24 excess at 0.225 s',
         record=_F38A + "/science_table.csv", grep='f38an24h1800comp,24,1\\.8,composition_only,24/24,2\\.025,11\\.822115551810377,'),
    dict(id="HA3", printed=11.79, unit="-", where='S14.6 Table S (held window, server) N = 48 excess at 0.225 s',
         record=_F38A + "/science_table.csv", grep='f38an48h1800comp,48,1\\.8,composition_only,24/24,2\\.024999999999999,11\\.793817566383598,'),
    dict(id="HA4", printed=23.91, unit="-", where='S14.6 Table S (held window, server) N = 12 excess at 0.9 s',
         record=_F38A + "/science_table.csv", grep='f38an12h1800comp,12,1\\.8,composition_only,24/24,2\\.0249999999999995,11\\.858604953566028,23\\.908965051592407,'),
    dict(id="HA5", printed=23.78, unit="-", where='S14.6 Table S (held window, server) N = 24 excess at 0.9 s',
         record=_F38A + "/science_table.csv", grep='f38an24h1800comp,24,1\\.8,composition_only,24/24,2\\.025,11\\.822115551810377,23\\.777335741856053,'),
    dict(id="HA6", printed=23.66, unit="-", where='S14.6 Table S (held window, server) N = 48 excess at 0.9 s',
         record=_F38A + "/science_table.csv", grep='f38an48h1800comp,48,1\\.8,composition_only,24/24,2\\.024999999999999,11\\.793817566383598,23\\.65526960425258,'),
    dict(id="HA7", printed=163.1, unit="-", where='S14.6 Table S (held window, server) N = 12 excess at 1.65 s',
         record=_F38A + "/science_table.csv", grep='f38an12h1800comp,12,1\\.8,composition_only,24/24,2\\.0249999999999995,11\\.858604953566028,23\\.908965051592407,163\\.07671232671115,'),
    dict(id="HA8", printed=153.6, unit="-", where='S14.6 Table S (held window, server) N = 24 excess at 1.65 s',
         record=_F38A + "/science_table.csv", grep='f38an24h1800comp,24,1\\.8,composition_only,24/24,2\\.025,11\\.822115551810377,23\\.777335741856053,153\\.60965080627662,'),
    dict(id="HA9", printed=146.6, unit="-", where='S14.6 Table S (held window, server) N = 48 excess at 1.65 s',
         record=_F38A + "/science_table.csv", grep='f38an48h1800comp,48,1\\.8,composition_only,24/24,2\\.024999999999999,11\\.793817566383598,23\\.65526960425258,146\\.57587575874336,'),
    dict(id="HA10", printed=0.423, unit="-", where='S14.6 Table S (held window, server) N = 12 s/R at 1.725 s',
         record=_F38A + "/science_table.csv", grep='f38an12h1800comp,12,1\\.8,composition_only,24/24,2\\.0249999999999995,11\\.858604953566028,23\\.908965051592407,163\\.07671232671115,1\\.7249999999999994,0\\.42324370701081804,'),
    dict(id="HA11", printed=0.424, unit="-", where='S14.6 Table S (held window, server) N = 24 s/R at 1.725 s',
         record=_F38A + "/science_table.csv", grep='f38an24h1800comp,24,1\\.8,composition_only,24/24,2\\.025,11\\.822115551810377,23\\.777335741856053,153\\.60965080627662,1\\.7249999999999999,0\\.4237437184206638,'),
    dict(id="HA12", printed=0.424, unit="-", where='S14.6 Table S (held window, server) N = 48 s/R at 1.725 s',
         record=_F38A + "/science_table.csv", grep='f38an48h1800comp,48,1\\.8,composition_only,24/24,2\\.024999999999999,11\\.793817566383598,23\\.65526960425258,146\\.57587575874336,1\\.7249999999999992,0\\.4242380553334914,'),
    dict(id="HA13", printed=0.386, unit="-", where='S14.6 Table S (held window, server) N = 12 s/R at 2.025 s',
         record=_F38A + "/science_table.csv", grep='f38an12h1800comp,12,1\\.8,composition_only,24/24,2\\.0249999999999995,11\\.858604953566028,23\\.908965051592407,163\\.07671232671115,1\\.7249999999999994,0\\.42324370701081804,0\\.38593584686987614,'),
    dict(id="HA14", printed=0.387, unit="-", where='S14.6 Table S (held window, server) N = 48 s/R at 2.025 s',
         record=_F38A + "/science_table.csv", grep='f38an48h1800comp,48,1\\.8,composition_only,24/24,2\\.024999999999999,11\\.793817566383598,23\\.65526960425258,146\\.57587575874336,1\\.7249999999999992,0\\.4242380553334914,0\\.38730148080271354,'),
    dict(id="HA15", printed=1.09e-12, unit="-", where='S14.6 largest cumulative ledger over the 18 jobs',
         record=_F38A + "/jobs_table.csv", grep='1\\.0911484609925196e\\-12'),
    dict(id="HA16", printed=0.13, page="0.13", unit="-", where='S14.6 / Sec. 6.1 12 against 24 cells at 0.9 s, composition',
         record=_F38A + "/analysis/dependence.json", grep='"signed_at_0p9":\\s+0\\.13162930973635412'),
    dict(id="HA17", printed=0.12, page="0.12", unit="-", where='S14.6 / Sec. 6.1 24 against 48 cells at 0.9 s, composition',
         record=_F38A + "/analysis/dependence.json", grep='"signed_at_0p9":\\s+0\\.12206613760347196'),
    dict(id="HA18", printed=1.0, page="1.0", unit="-", where='S14.6 composition minus correlated at 0.9 s, 12 cells',
         record=_F38A + "/analysis/dependence.json", grep='"diff_points":\\s+1\\.0040380639086095'),
    dict(id="HA19", printed=7.1e-12, unit="-", where='S14.6 twin front agreement, 12 cells',
         record=_F38A + "/analysis/twin_check.json", grep='"front_max_rel":\\s+7\\.121319245511469e\\-12'),
    dict(id="HA20", printed=2.4e-08, unit="-", where='S14.6 twin excess agreement, 12 cells, points',
         record=_F38A + "/analysis/twin_check.json", grep='"excess_max_abs_points":\\s+2\\.3661840486965957e\\-08'),
    dict(id="HA21", printed=1.8, page="1.8", unit="-", where='S14.6 correlated oracle at the center from 1.8 s',
         record=_F38A + "/item09/item09_f38an12h1800comb.csv", grep='return_hold,20,1\\.7999999999999994,[^\\n]*ZeroDivisionError'),
    dict(id="HA22", printed=1.17, unit="-", where='S14.6 successive-difference ratio at 0.9 s, high (temperature)',
         record=_F38AR, grep='1\\.08\\s+\\(composition\\),\\s+1\\.06\\s+\\(combined\\),\\s+1\\.17\\s+\\(temperature\\)'),
    dict(id="HA23", printed=0.19, page="0.19", unit="-", where='S14.6 resolved front between successive meshes at most 0.19 percent',
         record=_F38AR, grep='differs\\s+by\\s+at\\s+most\\s+0\\.19\\s+percent\\s+of\\s+its\\s+position'),
    dict(id="HA24", printed=3.7, page="3.7", unit="-", where='S14.6 oracle front differs by 3.7 and 2.9 percent at 1.65 s',
         record=_F38AR, grep='differs\\s+by\\s+3\\.7\\s+and\\s+2\\.9\\s+percent\\s+at\\s+1\\.65\\s+s'),
    dict(id="FC1", printed=0.6264, unit="-", where='Sec. 6.2 / 6.3 / S14.2 / S14.5 declared film corner',
         record=_B410, grep='\\*\\*declared\\s+corner\\*\\*\\s+\\|\\s+\\*\\*0\\.6264\\*\\*'),
    dict(id="FC2", printed=0.62631, unit="-", where='S14.5 smallest admitted multiplier',
         record=_B410, grep='\\*\\*0\\.6263060583204717\\s+on\\s+19\\s+of\\s+19\\s+states\\*\\*'),
    dict(id="FC3", printed=0.09997, unit="-", where='S14.5 ratio at the declared corner',
         record=_B410, grep='ratio\\s+0\\.09997000809857089'),
    dict(id="FC4", printed=0.157, unit="-", where='S14.5 ratio at 0.5',
         record=_B410, grep='\\|\\s+ratio\\s+at\\s+0\\.5\\s+\\(former\\s+corner\\)\\s+\\|\\s+0\\.15690371147557045,\\s+refused\\s+\\|'),
    dict(id="FC5", printed=59.46, unit="-", where='S14.5 film 1.0 against 2.0, cumulative boundary energy half-width (workstation)',
         record=_F37B + "/item08/twin_ranges_by_axis.csv", grep='film_multiplier,combined,cumulative_boundary_energy_out_j,.*,0\\.5945859436187543\\b'),
    dict(id="FC6", printed=2.125, unit="-", where='Sec. 6.3 / S14.5 film 1.0 against 2.0, peak particle temperature (workstation)',
         record=_F37B + "/item08/twin_ranges_by_axis.csv", grep='film_multiplier,combined,peak_particle_temperature_k,.*,2\\.1252519076634258,'),
    dict(id="FC7", printed=0.2185, unit="-", where='S14.5 conductivity axis, temperature only, workstation (Sec. 6.3 prints 0.22)',
         record=_F37B + "/item08/twin_ranges_by_axis.csv", grep='structural_conductivity,temperature_only,peak_particle_temperature_k,.*,0\\.2185235276444928,'),
    dict(id="TC1", printed=1.44, page="1.440", unit="-", where='Sec. 6.1 / S14.3 N = 96 like-for-like time chain at 0.225 s',
         record=_B410, grep='\\|\\s+0\\.225\\s+\\(t\\*\\)\\s+\\|\\s+yes\\s+\\|'),
    dict(id="TC2", printed=0.000433, unit="-", where='S14.3 finest pair of that chain',
         record=_B410, grep='The\\s+finest\\s+pair\\s+differs\\s+by\\s+4\\.33e\\-4\\s+relative\\s+at\\s+0\\.225\\s+s'),
    dict(id="TC3", printed=0.517, unit="-", where='S14.3 water inventory ratio at t*',
         record=_B410, grep='contract\\s+at\\s+t\\*\\s+\\(ratios\\s+0\\.517\\s+and\\s+0\\.378\\)'),
    dict(id="TC4", printed=3.165, unit="-", where='S14.3 interface temperature ratio at t*',
         record=_B410, grep='the\\s+interface\\s+temperature\\s+\\(3\\.165\\)\\s+do\\s+not'),
    dict(id="TC5", printed=2.154, unit="-", where='S14.3 the 12-cell composition chain at 0.15 s',
         record=_B410, grep='reads\\s+0\\.566\\s+at\\s+0\\.075\\s+s\\s+and\\s+2\\.154\\s+at\\s+0\\.15\\s+s'),
    dict(id="TC6", printed=15.7, page="15.7", unit="-", where='S14.3 the earlier reading at a face-split end',
         record=_B410, grep="15\\.675\\s+\\(the\\s+diagnosis's\\s+R\\-4,\\s+15\\.7\\)"),
    dict(id="TC7", printed=1.78e-13, unit="-", where='S14.3 ledgers of the refined N = 96 cells',
         record=_B410, grep='ledgers\\s+at\\s+most\\s+1\\.78e\\-13'),
    dict(id="TC8", printed=8.49e-11, unit="-", where="S14.3 the scalar probe's N = 48 root",
         record=_B410, grep='8\\.485064453787296e\\-11'),
    dict(id="TC9", printed=9.5e-14, unit="-", where='S14.3 the repaired adverse challenge, cumulative ledger',
         record=_ADVR, grep='9\\.505487584485193e\\-14'),
    dict(id="SG1", printed=69, page="69", unit="-", where='Sec. 7.3 / S14.2 stage 1: all 69 jobs complete',
         record=_S1, grep='all\\s+69\\s+screen\\s+jobs\\s+completed'),
    dict(id="SG2", printed=1.59, unit="-", where='S14.2 baseline hexane flux integral in space',
         record=_S1, grep='1\\.59\\s+percent\\s+in\\s+space\\s+and\\s+2\\.88\\s+percent\\s+in\\s+time'),
    dict(id="SG3", printed=2.88, unit="-", where='S14.2 baseline hexane flux integral in time',
         record=_S1, grep='1\\.59\\s+percent\\s+in\\s+space\\s+and\\s+2\\.88\\s+percent\\s+in\\s+time'),
    dict(id="TR85", printed=12.8, unit="-", where='Sec. 5 lead / Sec. 5.3 thesis runs at least 12.8 times the heat-limited time (sunflower, r_pre)',
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_heat_limited_rpre": 12\.770043649834854'),
    dict(id="TR86", printed=25.1, unit="-", page="25$ times", where='Sec. 5.3 thesis runs up to 25.1 times the heat-limited time (soybean, r_c), printed 12.8 to 25 since part 6 (S12)',
         record=_TR + "/two_regime.json", grep=r'"ratio_meas_over_heat_limited_rc": 25\.114285929626128'),
    dict(id="Z90", printed=0, unit="-", where="'the throttled end holds throughout' (regime-map check: no throttle on the total flux)",
         expect="absent", grep='throttled end\\s+holds\\s+throughout'),
    dict(id="Z91", printed=0, unit="-", where="'limited throughout by transport through the dry shell' (introduction)",
         expect="absent", grep='limited throughout by transport\\s+through'),
    dict(id="Z92", printed=0, unit="-", where="'mass path is internally' (Sec. 2.1 and its E copies)",
         expect="absent", grep='mass path is\\s+internally'),
    dict(id="Z93", printed=0, unit="-", where="'reduce to the scalar Fick description' (withdrawn; E16.2)",
         expect="absent", grep='reduce to the scalar Fick\\s+description'),
    dict(id="Z94", printed=0, unit="-", where="'recovering the isothermal description' (S5)",
         expect="absent", grep='recovering the isothermal\\s+description'),
    dict(id="Z95", printed=0, unit="-", where="'recession is transport-limited throughout' (Sec. 6 lead and its E copy)",
         expect="absent", grep='recession is transport-limited\\s+throughout'),
    dict(id="Z96", printed=0, unit="-", where="'isothermal description ... is recovered' (E5)",
         expect="absent", grep='isothermal description of\\s+\\\\citet\\{cardarelli2002modeling\\} is recovered'),
    dict(id="Z97", printed=0, unit="-", where="'two ends of one expression' (the formulation has no transport end)",
         expect="absent", grep='two ends of one\\s+expression'),
    dict(id="Z98", printed=0, unit="-", where="'sit on opposite sides of one inversion'",
         expect="absent", grep='sit on opposite\\s+sides of one\\s+inversion'),
    dict(id="Z99", printed=0, unit="-", where="the slab table's 0.200 at 0.50 R (slab value 0.196)",
         expect="absent", grep='0\\.50 & 443 & 0\\.200'),
    dict(id="Z100", printed=0, unit="-", where="'throttled by a growing dry shell' (regime-B definition)",
         expect="absent", grep='throttled by a growing\\s+dry shell'),
    dict(id="Z101", printed=0, unit="-", where="'eq.~(22) of the main' (the diffusivity law is eq. 20; by label now)",
         expect="absent", grep='eq\\.~\\(22\\) of the main'),
    dict(id="Z102", printed=0, unit="-", where='the literal main-text equation numbers 21 and 22 in the extended report',
         expect="absent", grep='\\[equation~\\(2[12]\\) of the main text\\]'),
    dict(id="Z103", printed=0, unit="-", where="'a separate obstacle under diagnosis' (the refined route is built)",
         expect="absent", grep='a separate obstacle under\\s+diagnosis'),
    dict(id="Z104", printed=0, unit="-", where="'re-scored probe not yet returned' (the repaired challenge passes)",
         expect="absent", grep='re-scored probe not yet\\s+returned'),
    dict(id="Z105", printed=0, unit="-", where="'scale inversely with' (the resistance ratios scale with the film)",
         expect="absent", grep='scale inversely\\s+with'),
    dict(id="Z106", printed=0, unit="-", where="'heat-to-transport transition is the crossing' (Sec. 6.2)",
         expect="absent", grep='heat-to-transport\\s+transition is the crossing'),
    dict(id="Z107", printed=0, unit="-", where="item 04 'stopping at the 96-cell bootstraps' (stage 1 now completes all 69 jobs)",
         expect="absent", grep='stopping at the 96-cell bootstraps'),
    dict(id="Z108", printed=0, unit="-", where="'for lack of a cut candidate, under diagnosis' (Sec. 6.1)",
         expect="absent", grep='cut candidate, under\\s+diagnosis'),
]
for _a in ASSERTIONS:
    if _a["id"] in ("R6", "R7"):
        _a["where"] = _a["where"] + " [film 0.5 refuses at any bound on this tree; the ruled corner 0.6264 not yet run; stated in Sec. 6.3 and S14.5 (audit 8)]"
    if _a["id"] == "Y58":
        _a["where"] = _a["where"] + " [the server's value; the workstation reads 0.2185 on another path, FC7 (audit 8)]"


# ---------------------------------------------------------------------------
# 2026-09-30, REFEREE_AUDIT_8 part 2: the moving-front march of the journal
# runs of Faner et al. (2019), Sec. 5.4, Table 5, Figure 3, the abstract,
# Secs. 5, 7.2, 7.3, the conclusions and S9.16.  Records: the emergency march,
# step 1 (landed) and step 2 (its record and sidecar land in the chain after
# this block was written: a staging variable in PENDING_ROOTS, as in the
# blocks above).  MF1-MF159 bind the printed numbers; MZ1-MZ10 the withdrawn
# wording.  Every artifact behind them is TEMPORARILY_UNCERTIFIED_EMERGENCY.
# ---------------------------------------------------------------------------
_MS1 = "../docs/evidence/paper1_faner2019_march_emergency_2026-09-30"
_MS1R = "../docs/GT_PS2_FANER2019_MOVING_FRONT_MARCH_EMERGENCY_2026-09-30.md"
_MS2 = "../docs/evidence/paper1_faner2019_march_step2_2026-09-30"     # scratchpad faner2019_march/
_MS2R = "../docs/GT_PS2_FANER2019_MARCH_STEP2_2026-09-30.md"          # scratchpad faner2019_march/docs/
PENDING_ROOTS[_MS1] = "DTDC_PAPER_MARCH_STEP1_EVIDENCE"
PENDING_ROOTS[_MS1R] = "DTDC_PAPER_MARCH_STEP1_RECORD"
PENDING_ROOTS[_MS2] = "DTDC_PAPER_MARCH_STEP2_EVIDENCE"
PENDING_ROOTS[_MS2R] = "DTDC_PAPER_MARCH_STEP2_RECORD"
EXTRA_CORPUS_ROOTS += [_MS1, _MS2]
ASSERTIONS += [
    dict(id='MF1', printed=14.15, unit="-", where='Abstract/Sec. 5.4/Table 5 soybean marched crossing to 0.10, 24 cells, dt 0.075 (L3)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 14\\.14666075013411'),
    dict(id='MF2', printed=20.35, unit="-", where='Sec. 5.4/Table 5 sunflower marched crossing to 0.10, 24 cells, dt 0.075 (L4)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 20\\.350317295584063'),
    dict(id='MF3', printed=25.35, unit="-", where='Sec. 5.4/Table 5 soybean measured crossing to 0.10',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"measured_crossing_to_0p10_s": 25\\.350441949123166'),
    dict(id='MF4', printed=42.5, page='42.50', unit="-", where='Sec. 5.4/Table 5 sunflower measured crossing to 0.10',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"measured_crossing_to_0p10_s": 42\\.50008976874932'),
    dict(id='MF5', printed=0.56, unit="-", where='Abstract/Sec. 5/5.4/Table 2/5/conclusions marched over measured, soybean (0.558 at 24 cells)',
         record=_MS2R + '', grep='0\\.558 \\(soybean\\) and 0\\.479\\s+\\(sunflower\\)'),
    dict(id='MF6', printed=0.48, unit="-", where='Abstract/Sec. 5/5.4/Table 2/5/conclusions marched over measured, sunflower (0.479 at 24 cells)',
         record=_MS2R + '', grep='0\\.558 \\(soybean\\) and 0\\.479\\s+\\(sunflower\\)'),
    dict(id='MF7', printed=27, unit="-", where='Sec. 5.4/S9.16 face events in the ladder (26 completed + 1 refused)',
         record=_MS2R + '', grep='exercised on\\s+27 faces'),
    dict(id='MF8', printed=1.4e-11, page='1.4', unit="-", where='Sec. 5.4 largest face-event scaled residual',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 1\\.3895058786824721e-11'),
    dict(id='MF9', printed=2.6e-13, page='2.6', unit="-", where='Sec. 5.4 largest face-event step ledger',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"step_ledger": 2\\.6370518049044265e-13'),
    dict(id='MF10', printed=15.15, unit="-", where='Sec. 5.4/S9.16 refused march (L6) end time',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"t_end_s": 15\\.149999999999903'),
    dict(id='MF11', printed=0.1245, unit="-", where='Sec. 5.4/S9.16 refused march (L6) loading',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"X_end": 0\\.12453241245237281'),
    dict(id='MF12', printed=0.0014, page='1.4', unit="-", where='Sec. 5.4/Table 5 soybean spread 12 against 24 cells, s',
         record=_MS2R + '', grep='they move by 1\\.4e-3 s and 1\\.9e-3 s'),
    dict(id='MF13', printed=0.0019, page='1.9e-3', unit="-", where='Sec. 5 lead/Sec. 5.4/Table 5 sunflower spread 12 against 24 cells, s',
         record=_MS2R + '', grep='they move by 1\\.4e-3 s and 1\\.9e-3 s'),
    dict(id='MF14', printed=0.00029, page='2.9e-4', unit="-", where='Sec. 5.4/Table 5 largest spread over the two steps, s',
         record=_MS2R + '', grep='Halving dt moves them by at most\\s+2\\.9e-4 s'),
    dict(id='MF15', printed=68.75, unit="-", where='Sec. 5.4/Table 5 marched interface temperature, C',
         record=_MS1 + '/evidence/plateau.json', grep='"marched_interface_temperature_c_range": \\[\\s*68\\.75335616347832'),
    dict(id='MF16', printed=69.1, unit="-", where='Sec. 5.4/Table 5 measured plateau, soybean, C',
         record=_MS1 + '/evidence/plateau.json', grep='"measured_plateau_mean_c_to_crossing": 69\\.11449999999999'),
    dict(id='MF17', printed=68.5, unit="-", where='Sec. 5.4/Table 5 measured plateau, sunflower, C',
         record=_MS1 + '/evidence/plateau.json', grep='"measured_plateau_mean_c_to_crossing": 68\\.51'),
    dict(id='MF18', printed=0.52, unit="-", where='Sec. 5.4 surface above the interface at 0.10, soybean, K',
         record=_MS2R + '', grep='0\\.52 and 0\\.43 K\\s+above the interface'),
    dict(id='MF19', printed=0.43, unit="-", where='Sec. 5.4 surface above the interface at 0.10, sunflower, K',
         record=_MS2R + '', grep='0\\.52 and 0\\.43 K\\s+above the interface'),
    dict(id='MF20', printed=0.063, unit="-", where='Sec. 5.4/Table 5 measured plateau end loading, soybean',
         record=_MS1 + '/evidence/plateau.json', grep='"measured_plateau_end_loading": 0\\.06325968872941157'),
    dict(id='MF21', printed=0.085, unit="-", where='Sec. 5.4/Table 5 measured plateau end loading, sunflower',
         record=_MS1 + '/evidence/plateau.json', grep='"measured_plateau_end_loading": 0\\.08540557964660069'),
    dict(id='MF22', printed=0.57, unit="-", where='Sec. 5.4/Table 5 marched mean rise at the plateau-end loading, soybean (2 cells), K',
         record=_MS1 + '/evidence/plateau.json', grep='"marched_mean_T_rise_at_that_loading_K": 0\\.5656988495848623'),
    dict(id='MF23', printed=0.27, unit="-", where='Sec. 5.4/Table 5 marched mean rise at the plateau-end loading, sunflower (2 cells), K',
         record=_MS1 + '/evidence/plateau.json', grep='"marched_mean_T_rise_at_that_loading_K": 0\\.26638279702984846'),
    dict(id='MF24', printed=1.1, unit="-", where='Sec. 5.4 marched mean rise by 0.041 kg/kg, at most (soybean +1.10, sunflower +0.92), K',
         record=_MS1R + '', grep='the rise is \\+1\\.10 K and \\+0\\.92 K'),
    dict(id='MF25', printed=0.041, unit="-", where='Sec. 5.4 end loading of the 2-cell marches',
         record=_MS1 + '/evidence/plateau.json', grep='"march_end_loading": 0\\.04074283405161792'),
    dict(id='MF26', printed=13.99, unit="-", where='Sec. 5.4/S9.16 film coefficient with the Ackermann factor, soybean, W/(m2 K)',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"h_with_model_ackermann_w_m2_k": 13\\.991080612041394'),
    dict(id='MF27', printed=16.41, unit="-", where='Sec. 5.4/S9.16 film coefficient with the Ackermann factor, sunflower, W/(m2 K)',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"h_with_model_ackermann_w_m2_k": 16\\.41150781648396'),
    dict(id='MF28', printed=0.1995, unit="-", where='Sec. 5.4/S9.16 start loading at the front 0.999 R',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"X_start": 0\\.19953429803323314'),
    dict(id='MF29', printed=0.05, unit="-", where='Sec. 5.4/S9.16 retained water held (declared, inside the qualified Luikov span)',
         record=_MS1R + '', grep='default 0\\.05 kg/kg, inside the'),
    dict(id='MF30', printed=16.18, unit="-", where='Sec. 5.4 (16.2)/Table 5/Table S soybean, film without the Ackermann factor (P6, 4 cells)',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 16\\.18064951129601'),
    dict(id='MF31', printed=22.14, unit="-", where='Sec. 5.4 (22.1)/Table 5/Table S sunflower, film without the Ackermann factor (4 cells)',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 22\\.136711776895908'),
    dict(id='MF32', printed=14.12, unit="-", where='Sec. 5.4/Table 5/Table S soybean at conductivity 0.29 (4 cells)',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 14\\.120077532089681'),
    dict(id='MF33', printed=20.31, unit="-", where='Sec. 5.4/Table 5/Table S sunflower at conductivity 0.29 (4 cells)',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 20\\.31134859045187'),
    dict(id='MF34', printed=0.034, unit="-", where='Sec. 5.4/Table S heat Biot number, soybean',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"heat_biot_h_eff_R_over_k": 0\\.03387462309263205'),
    dict(id='MF35', printed=0.045, unit="-", where='Sec. 5.4/Table S heat Biot number, sunflower',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"heat_biot_h_eff_R_over_k": 0\\.04531163467457226'),
    dict(id='MF36', printed=0.275, unit="-", where='Sec. 5.4/Table 5 onset loading, soybean (0.27465 in S9.16)',
         record=_MS2 + '/evidence/onset_reading.json', grep='"onset_loading_kg_kg": 0\\.27465229762394155'),
    dict(id='MF37', printed=0.243, unit="-", where='Sec. 5.4/Table 5 onset loading, sunflower (0.24282 in S9.16)',
         record=_MS2 + '/evidence/onset_reading.json', grep='"onset_loading_kg_kg": 0\\.2428'),
    dict(id='MF38', printed=9.43, unit="-", where='Sec. 5.4 (9.4)/Table 5/Table S soybean from the onset, 12 cells',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 9\\.428063023652967'),
    dict(id='MF39', printed=15.26, unit="-", where='Sec. 5.4 (15.3)/Table 5/Table S sunflower from the onset, 12 cells',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 15\\.25679176434619'),
    dict(id='MF40', printed=0.986, unit="-", where='Sec. 5.4/Table 5/S9.16 model rate over measured, 0.3 to 3 s, soybean',
         record=_MS2R + '', grep='It is 0\\.986 and 0\\.982 over 0\\.3 to 3\\.0 s'),
    dict(id='MF41', printed=0.982, unit="-", where='Sec. 5.4/Table 5/S9.16 model rate over measured, 0.3 to 3 s, sunflower',
         record=_MS2R + '', grep='It is 0\\.986 and 0\\.982 over 0\\.3 to 3\\.0 s'),
    dict(id='MF42', printed=0.039, unit="-", where="Sec. 5.4/S9.16 wet start below the model's interface root, K",
         record=_MS2R + '', grep='which is 0\\.039 K'),
    dict(id='MF43', printed=1.28, unit="-", where='Sec. 5.4 sorption-heat upper bound on the time factor (1.278 in S9.16)',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"upper_bound_time_factor_if_all_desorption_heat_charged_in_window": 1\\.2777050254748477'),
    dict(id='MF44', printed=1.027, unit="-", where='Sec. 5.4/7.2/Table S marched over the pure film limit',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"marched_over_film_limit": 1\\.02723647098005'),
    dict(id='MF45', printed=0.72, unit="-", where='Sec. 5.4/7.2 measured rate before the crossing over r_c, soybean',
         record=_TR + '/derived_checks.json', grep='"rate_pre_over_const":\\ 0\\.7170100660124655'),
    dict(id='MF46', printed=0.61, unit="-", where='Sec. 5.4/7.2 measured rate before the crossing over r_c, sunflower',
         record=_TR + '/derived_checks.json', grep='"rate_pre_over_const":\\ 0\\.613316830126657'),
    dict(id='MF47', printed=0.57, unit="-", where='Sec. 5.4/Table 5 RMS over the covered window, soybean (L3)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"rms_over_span": 0\\.5726658800735105'),
    dict(id='MF48', printed=0.98, unit="-", where='Sec. 5.4/Table 5 RMS over the covered window, sunflower (L4)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"rms_over_span": 0\\.9818992710437606'),
    dict(id='MF49', printed=0.4, unit="-", where='Table 2 interface within 0.4 K of the measured plateaux (0.36 below, 0.24 above)',
         record=_MS1R + '', grep='0\\.36 K below and 0\\.24 K above'),
    dict(id='MF50', printed=13.77, unit="-", where='Table 5/Table S pure film limit, soybean, s',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"film_limit_crossing_to_0p10_s": 13\\.772960204404244'),
    dict(id='MF51', printed=19.81, unit="-", where='Table 5/Table S pure film limit, sunflower, s',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"film_limit_crossing_to_0p10_s": 19\\.809239752211866'),
    dict(id='MF52', printed=18.08, unit="-", where='Table 5/Table S sorption-heat upper bound, soybean, s',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"upper_bound_crossing_to_0p10_s": 18\\.07708190588503'),
    dict(id='MF53', printed=26.0, page='26.00', unit="-", where='Table 5/Table S sorption-heat upper bound, sunflower, s',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"upper_bound_crossing_to_0p10_s": 26\\.004150360623594'),
    dict(id='MF54', printed=382, unit="-", where='S9.16 option-off identity, recorded cells (0 differ)',
         record=_MS2 + '/evidence/knob_off_identity_K2_F1to5.json', grep='"cells_compared": 382,\\s+"cells_differing": 0'),
    dict(id='MF55', printed=149, unit="-", where='S9.16 face-event test files passed with every switch off',
         record=_MS2 + '/runs/U1_spot_tests_F1_F5.stdout.txt', grep='149 passed'),
    dict(id='MF56', printed=1.6e-15, page='1.6', unit="-", where='S9.16 water ledger at or below',
         record=_MS1R + '', grep='stayed at or below 1\\.6e-15'),
    dict(id='MF57', printed=341.9034, unit="-", where="S9.16 the pure vapor's interface root, K",
         record=_MS1R + '', grep='interface root, 341\\.9034 K'),
    dict(id='MF58', printed=350.0, page='350.0', unit="-", where='S9.16 interface-temperature chart upper bound, K',
         record=_MS1R + '', grep='\\(334\\.5, 350\\.0\\) K'),
    dict(id='MF59', printed=341.0, page='341.0', unit="-", where="S9.16 binary chart's upper bound, K",
         record=_MS1R + '', grep='\\(334\\.5, 341\\.0\\) K'),
    dict(id='MF60', printed=0.00726, page='7.26', unit="-", where='S9.16 measured constant rate, soybean, 1/s',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"measured_constant_rate_1_s": 0\\.007260603277429244'),
    dict(id='MF61', printed=0.00505, page='5.05e-3', unit="-", where='S9.16 measured constant rate, sunflower, 1/s',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"measured_constant_rate_1_s": 0\\.005048149310668734'),
    dict(id='MF62', printed=341.8645, unit="-", where='S9.16 model boiling point and wet start, K',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"model_boiling_point_k": 341\\.86451253613916'),
    dict(id='MF63', printed=334.93, unit="-", where='S9.16 enthalpy of vaporization, kJ/kg',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"model_dh_vap_j_kg": 334929\\.46668371576'),
    dict(id='MF64', printed=47.42, unit="-", where='S9.16 film supply h a, soybean, W/(kg K)',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"h_times_area_per_kg_w_kg_k": 47\\.41672749479113'),
    dict(id='MF65', printed=54.04, unit="-", where='S9.16 film supply h a, sunflower, W/(kg K)',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"h_times_area_per_kg_w_kg_k": 54\\.04339498673724'),
    dict(id='MF66', printed=0.11966, unit="-", where='S9.16 film multiplier with the Ackermann factor, soybean',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"multiplier_with_ackermann": 0\\.11965745545344508'),
    dict(id='MF67', printed=0.14517, unit="-", where='S9.16 film multiplier with the Ackermann factor, sunflower',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"multiplier_with_ackermann": 0\\.145168314706609'),
    dict(id='MF68', printed=0.1044, page='0.10440', unit="-", where='S9.16 film multiplier without the Ackermann factor, soybean',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"multiplier_no_ackermann": 0\\.10439987738892059'),
    dict(id='MF69', printed=0.13332, unit="-", where='S9.16 film multiplier without the Ackermann factor, sunflower',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"multiplier_no_ackermann": 0\\.13332389336134887'),
    dict(id='MF70', printed=27.55, unit="-", where='S9.16 forcing time X_c/r_c, soybean, s',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"forcing_timescale_s_Xc_over_rate": 27\\.545920408808488'),
    dict(id='MF71', printed=39.62, unit="-", where='S9.16 forcing time X_c/r_c, sunflower, s',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"forcing_timescale_s_Xc_over_rate": 39\\.61847950442373'),
    dict(id='MF72', printed=0.00587, unit="-", where='S9.16 relaxation-to-forcing ratio, soybean',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"relaxation_to_forcing_ratio_with_ackermann": 0\\.0058713333845248'),
    dict(id='MF73', printed=0.00297, unit="-", where='S9.16 relaxation-to-forcing ratio, sunflower',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"relaxation_to_forcing_ratio_with_ackermann": 0\\.0029668922576557614'),
    dict(id='MF74', printed=23.0, page='23.0', unit="-", where='S9.16 voidage Reynolds number of the film correlation, soybean',
         record=_MS1 + '/evidence/derive_film_2019_soybean.json', grep='"coletto_pair_reynolds_voidage": 23\\.04094134754'),
    dict(id='MF75', printed=25.0, page='25.0', unit="-", where='S9.16 voidage Reynolds number of the film correlation, sunflower',
         record=_MS1 + '/evidence/derive_film_2019_sunflower.json', grep='"coletto_pair_reynolds_voidage": 24\\.961019793168333'),
    dict(id='MF76', printed=0.586, unit="-", where='S9.16 soybean RMS matched at the crossing (0.571 matched at the start loading)',
         record=_MS1R + '', grep='soybean RMS 0\\.571 and\\s+0\\.586'),
    dict(id='MF77', printed=3.064e-10, page='3.064', unit="-", where='S9.16 bootstrap refusal, 12 cells',
         record=_MS2R + '', grep='3\\.064e-10 > 1e-10'),
    dict(id='MF78', printed=5.576e-10, page='5.576', unit="-", where='S9.16 bootstrap refusal, 24 cells',
         record=_MS2R + '', grep='\\(N = 12\\) and 5\\.576e-10 \\(N = 24\\)'),
    dict(id='MF79', printed=14.148, unit="-", where='Table S soybean 12 cells, dt 0.075 and 0.0375 (L1, L5)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 14\\.148087035320884'),
    dict(id='MF80', printed=14.147, unit="-", where='Table S soybean 24 cells, dt 0.075 (L3)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 14\\.14666075013411'),
    dict(id='MF81', printed=14.146, unit="-", where='Table S soybean 24 cells, dt 0.0375 (L7)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 14\\.146433748759534'),
    dict(id='MF82', printed=20.352, unit="-", where='Table S sunflower 12 cells, dt 0.075 (L2)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 20\\.352232981911754'),
    dict(id='MF83', printed=20.35, page='20.350', unit="-", where='Table S sunflower 24 cells (L4, L8)',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"marched_crossing_to_0p10_s": 20\\.35003100539005'),
    dict(id='MF84', printed=2.9e-11, page='2.90', unit="-", where='Table S largest scaled residual, L1',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 2\\.9023237938306313e-11'),
    dict(id='MF85', printed=7.68e-11, page='7.68', unit="-", where='Table S largest scaled residual, L5',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 7\\.683888529901792e-11'),
    dict(id='MF86', printed=3.53e-11, page='3.53', unit="-", where='Table S largest scaled residual, L3',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 3\\.533897346395623e-11'),
    dict(id='MF87', printed=7.9e-11, page='7.90', unit="-", where='Table S largest scaled residual, L7',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 7\\.896085626781921e-11'),
    dict(id='MF88', printed=8.17e-11, page='8.17', unit="-", where='Table S largest scaled residual, L2',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 8\\.167888949550835e-11'),
    dict(id='MF89', printed=8.6e-11, page='8.60', unit="-", where='Table S largest scaled residual, L6',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 8\\.599077739146427e-11'),
    dict(id='MF90', printed=3.97e-11, page='3.97', unit="-", where='Table S largest scaled residual, L4',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 3\\.970398305152482e-11'),
    dict(id='MF91', printed=7.51e-11, page='7.51', unit="-", where='Table S largest scaled residual, L8',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_scaled_residual": 7\\.506080279455198e-11'),
    dict(id='MF92', printed=9.87e-13, page='9.87', unit="-", where='Table S step ledger, L1',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 9\\.869431328253435e-13'),
    dict(id='MF93', printed=1.33e-12, page='1.33', unit="-", where='Table S cumulative ledger, L1',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 1\\.3331450656753897e-12'),
    dict(id='MF94', printed=9.97e-12, page='9.97', unit="-", where='Table S step ledger, L5',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 9\\.971231335234605e-12'),
    dict(id='MF95', printed=9.7e-11, page='9.70', unit="-", where='Table S cumulative ledger, L5',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 9\\.703322886446678e-11'),
    dict(id='MF96', printed=1.03e-12, page='1.03', unit="-", where='Table S step ledger, L3',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 1\\.0250135474964306e-12'),
    dict(id='MF97', printed=1.57e-12, page='1.57', unit="-", where='Table S cumulative ledger, L3',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 1\\.5686761239648605e-12'),
    dict(id='MF98', printed=1.12e-12, page='1.12', unit="-", where='Table S step ledger, L7',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 1\\.1165623651284726e-12'),
    dict(id='MF99', printed=2.86e-12, page='2.86', unit="-", where='Table S cumulative ledger, L7',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 2\\.860900312206035e-12'),
    dict(id='MF100', printed=4.39e-12, page='4.39', unit="-", where='Table S step ledger, L2',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 4\\.390614265697943e-12'),
    dict(id='MF101', printed=3.02e-11, page='3.02', unit="-", where='Table S cumulative ledger, L2',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 3\\.019353796594835e-11'),
    dict(id='MF102', printed=1.59e-11, page='1.59', unit="-", where='Table S step ledger, L6',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 1\\.587993579808176e-11'),
    dict(id='MF103', printed=9.9996e-11, page='9.9996', unit="-", where='Table S cumulative ledger, L6',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 9\\.999602603747037e-11'),
    dict(id='MF104', printed=6.93e-13, page='6.93', unit="-", where='Table S step ledger, L4',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 6\\.929408282826315e-13'),
    dict(id='MF105', printed=1.16e-12, page='1.16', unit="-", where='Table S cumulative ledger, L4',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 1\\.1551364555108297e-12'),
    dict(id='MF106', printed=4.26e-12, page='4.26', unit="-", where='Table S step ledger, L8',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_step_ledger": 4\\.2644748542357325e-12'),
    dict(id='MF107', printed=8.08e-11, page='8.08', unit="-", where='Table S cumulative ledger, L8',
         record=_MS2 + '/evidence/step2/ladder.json', grep='"max_cumulative_ledger": 8\\.084025256418155e-11'),
    dict(id='MF108', printed=12.75, page='12.750', unit="-", where='Table S last face stride end, L3',
         record=_MS2R + '', grep='169, 12\\.750 \\(19\\)'),
    dict(id='MF109', printed=18.262, unit="-", where='Table S last face stride end, L8',
         record=_MS2R + '', grep='485, 18\\.262 \\(19\\)'),
    dict(id='MF110', printed=0.0058, unit="-", where='S9.16 shortest face event, s',
         record=_MS2R + '', grep='Events of 0\\.0058 to 0\\.0717 s'),
    dict(id='MF111', printed=0.0666, unit="-", where='S9.16 longest departure remainder, s',
         record=_MS2R + '', grep='departure remainders 0\\.0033 to 0\\.0666 s'),
    dict(id='MF112', printed=1.39e-11, page='1.39', unit="-", where='S9.16 face-event scaled residual, largest',
         record=_MS2R + '', grep='Maximum scaled residual 2\\.6e-14 to 1\\.39e-11'),
    dict(id='MF113', printed=2.64e-13, page='2.64', unit="-", where='S9.16 face-event step ledger, largest',
         record=_MS2R + '', grep='step ledger at most 2\\.64e-13'),
    dict(id='MF114', printed=2400000.0, page='2.4', unit="-", where='S9.16 condition proxy, largest',
         record=_MS2R + '', grep='condition proxy 8\\.4e5 to 2\\.4e6'),
    dict(id='MF115', printed=341.9033561634786, unit="-", where='S9.16 interface temperature on every stride, upper, K',
         record=_MS2R + '', grep='341\\.9033561634781 to 341\\.9033561634786 K'),
    dict(id='MF116', printed=4.6e-12, page='4.6', unit="-", where='S9.16 energy cumulative ledger, at most',
         record=_MS2R + '', grep='energy cumulative ledger at most 4\\.6e-12'),
    dict(id='MF117', printed=9.84e-11, page='9.84', unit="-", where='S9.16 L6 cumulative ledger at 2.25 s',
         record=_MS2R + '', grep='already 9\\.84e-11 at 2\\.25 s'),
    dict(id='MF118', printed=98, unit="-", where='S9.16 share of the ledger carried by 28 strides above 1e-11, percent',
         record=_MS2R + '', grep='28 strides accepted above 1e-11 carry 98 percent'),
    dict(id='MF119', printed=4e-14, page='4.0', unit="-", where='S9.16 median accepted scaled residual, upper',
         record=_MS2R + '', grep='median accepted scaled residual is 2\\.2e-14 to 4\\.0e-14'),
    dict(id='MF120', printed=134, unit="-", where='S9.16 strides still needed to 0.10 (derived)',
         record=_MS2R + '', grep='about 134 more strides'),
    dict(id='MF121', printed=2.5e-13, page='2.5', unit="-", where='S9.16 observed mean addition per stride',
         record=_MS2R + '', grep='1\\.9e-13 from the start instead of the observed 2\\.5e-13'),
    dict(id='MF122', printed=4.3e-11, page='4.3', unit="-", where='S9.16 residual cap that would suffice (derived)',
         record=_MS2R + '', grep='at about 4\\.3e-11'),
    dict(id='MF123', printed=12.675, unit="-", where='S9.16 onset march, sunflower, 4 cells, refusal time',
         record=_MS2R + '', grep='12\\.675 s, X 0\\.1798'),
    dict(id='MF124', printed=0.1798, unit="-", where='S9.16 onset march, sunflower, 4 cells, refusal loading',
         record=_MS2R + '', grep='12\\.675 s, X 0\\.1798'),
    dict(id='MF125', printed=3e-11, page='3.0', unit="-", where='S9.16 ledgers at the end of the 0.075 s runs, upper',
         record=_MS2R + '', grep='end at 1\\.2e-12 to 3\\.0e-11'),
    dict(id='MF126', printed=14.16, unit="-", where='Table S baseline soybean, 4 cells',
         record=_MS2R + '', grep='14\\.16 s \\(0\\.559\\)'),
    dict(id='MF127', printed=0.559, unit="-", where='Table S baseline soybean ratio',
         record=_MS2R + '', grep='14\\.16 s \\(0\\.559\\)'),
    dict(id='MF128', printed=20.37, unit="-", where='Table S baseline sunflower, 4 cells',
         record=_MS2R + '', grep='20\\.37 s \\(0\\.479\\)'),
    dict(id='MF129', printed=9.45, unit="-", where='Table S onset, soybean, 4 cells',
         record=_MS2 + '/evidence/step2/hypotheses.json', grep='"marched_crossing_to_0p10_s": 9\\.448496842782063'),
    dict(id='MF130', printed=0.373, unit="-", where='Table S onset, soybean ratio',
         record=_MS2R + '', grep='9\\.45 s \\(0\\.373\\)'),
    dict(id='MF131', printed=0.359, unit="-", where='Table S onset, sunflower ratio, 12 cells',
         record=_MS2R + '', grep='15\\.26 s \\(0\\.359\\)'),
    dict(id='MF132', printed=0.638, unit="-", where='Table S film without Ackermann, soybean ratio',
         record=_MS2R + '', grep='16\\.18 s \\(0\\.638'),
    dict(id='MF133', printed=0.521, unit="-", where='Table S film without Ackermann, sunflower ratio',
         record=_MS2R + '', grep='22\\.14 s \\(0\\.521\\)'),
    dict(id='MF134', printed=18, unit="-", where='Table S share of the gap closed by (b), soybean, percent',
         record=_MS2R + '', grep='closes 18 and 8 percent of the gap'),
    dict(id='MF135', printed=0.557, unit="-", where='Table S conductivity 0.29, soybean ratio',
         record=_MS2R + '', grep='14\\.12 s \\(0\\.557\\)'),
    dict(id='MF136', printed=0.478, unit="-", where='Table S conductivity 0.29, sunflower ratio',
         record=_MS2R + '', grep='20\\.31 s \\(0\\.478\\)'),
    dict(id='MF137', printed=0.713, unit="-", where='Table S sorption heat upper bound, soybean ratio',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"upper_bound_ratio_to_measured": 0\\.7130874460557596'),
    dict(id='MF138', printed=0.612, unit="-", where='Table S sorption heat upper bound, sunflower ratio',
         record=_MS2 + '/evidence/step2/derived_bounds.json', grep='"upper_bound_ratio_to_measured": 0\\.6118610690498981'),
    dict(id='MF139', printed=0.27465, unit="-", where='S9.16 onset loading, soybean',
         record=_MS2 + '/evidence/onset_reading.json', grep='"onset_loading_kg_kg": 0\\.27465229762394155'),
    dict(id='MF140', printed=30.62, unit="-", where='S9.16 onset time, soybean, s',
         record=_MS2 + '/evidence/onset_reading.json', grep='"onset_time_s": 30\\.61634981002494'),
    dict(id='MF141', printed=0.24282, unit="-", where='S9.16 onset loading, sunflower',
         record=_MS2R + '', grep='Onset X 0\\.24282 at 51\\.29 s'),
    dict(id='MF142', printed=51.29, unit="-", where='S9.16 onset time, sunflower, s',
         record=_MS2R + '', grep='Onset X 0\\.24282 at 51\\.29 s'),
    dict(id='MF143', printed=0.020466, unit="-", where='S9.16 sorbed remainder at unit activity, kg/kg',
         record=_MS2R + '', grep='0\\.020466 kg/kg at a = 1'),
    dict(id='MF144', printed=1.0009, unit="-", where='S9.16 film heat over the derived supply, soybean, strides 1 to 40',
         record=_MS2R + '', grep='strides 1 to 40 \\(to 3\\.0 s\\) \\| 1\\.0009'),
    dict(id='MF145', printed=1.0001, unit="-", where='S9.16 film heat over the derived supply, sunflower, strides 1 to 40',
         record=_MS2R + '', grep='sunflower, strides 1 to 40 \\| 1\\.0001'),
    dict(id='MF146', printed=0.9829, unit="-", where='S9.16 latent of the mobile liquid, soybean',
         record=_MS2R + '', grep='\\| 0\\.9829 \\| 0\\.0127 \\|'),
    dict(id='MF147', printed=0.9767, unit="-", where='S9.16 latent of the mobile liquid, sunflower',
         record=_MS2R + '', grep='\\| 0\\.9767 \\| 0\\.0170 \\|'),
    dict(id='MF148', printed=0.0127, unit="-", where='S9.16 sensible share, soybean',
         record=_MS2R + '', grep='\\| 0\\.9829 \\| 0\\.0127 \\|'),
    dict(id='MF149', printed=0.017, page='0.0170', unit="-", where='S9.16 sensible share, sunflower',
         record=_MS2R + '', grep='\\| 0\\.9767 \\| 0\\.0170 \\|'),
    dict(id='MF150', printed=0.0018, unit="-", where='S9.16 desorption share, soybean (0.0015 + 0.0003, derived sum)',
         record=_MS2R + '', grep='\\| 0\\.0015 \\+ 0\\.0003 \\|'),
    dict(id='MF151', printed=0.0019, unit="-", where='S9.16 desorption share, sunflower (0.0015 + 0.0004, derived sum)',
         record=_MS2R + '', grep='\\| 0\\.0015 \\+ 0\\.0004 \\|'),
    dict(id='MF152', printed=0.885, unit="-", where='S9.16 model rate over measured on the first stride',
         record=_MS2R + '', grep="The model's rate is 0\\.885 of the measured on the first stride"),
    dict(id='MF153', printed=3.5e-15, page='3.5', unit="-", where='S9.16 binding deficits equal, relative',
         record=_MS2R + '', grep='3\\.5e-15 relative at T_Gamma'),
    dict(id='MF154', printed=2446.37, unit="-", where='S9.16 binding deficit at unit activity, J/kg dry',
         record=_MS2R + '', grep='2446\\.37 J/kg dry'),
    dict(id='MF155', printed=9301, unit="-", where='S9.16 desorption heat of the remainder, J/kg dry',
         record=_MS2R + '', grep='9301 J/kg dry, of which 26 percent is binding'),
    dict(id='MF156', printed=33493, unit="-", where='S9.16 mobile latent heat from the crossing to 0.10, J/kg dry',
         record=_MS2R + '', grep='against 33493 J/kg'),
    dict(id='MF157', printed=1.278, unit="-", where='S9.16 upper-bound time factor',
         record=_MS2R + '', grep='by the factor 1\\.278'),
    dict(id='MF158', printed=59.92, unit="-", where='S9.16 covered samples, soybean, s',
         record=_MS2R + '', grep='50\\.04 and 59\\.92 s for soybean, 70\\.12 and 80\\.00 s for sunflower'),
    dict(id='MF159', printed=80.0, page='80.00', unit="-", where='S9.16 covered samples, sunflower, s',
         record=_MS2R + '', grep='50\\.04 and 59\\.92 s for soybean, 70\\.12 and 80\\.00 s for sunflower'),
    dict(id='MZ1', expect="absent", where='Sec. 5 lead: the march is now reported (Sec. 5.4)', grep='A march of the moving-front\\s+system against the measured curves is not reported here'),
    dict(id='MZ2', expect="absent", where='Sec. 5.3 end: replaced by the pointer to Sec. 5.4', grep='A march of that system through these records is not reported here'),
    dict(id='MZ3', expect="absent", where="Sec. 7.3: replaced by the march's status", grep='none comes from a march of the moving front'),
    dict(id='MZ4', expect="absent", where='conclusions: replaced by the march', grep='compared in\\s+closed form and not marched'),
    dict(id='MZ5', expect="absent", where='Sec. 5 lead: withdrawn', grep='No marched number compared with experiment comes from the moving-front'),
    dict(id='MZ6', expect="absent", where='Sec. 5 lead: a seventh comparison now exists', grep='classification rather than a seventh comparison'),
    dict(id='MZ7', expect="absent", where='Sec. 5 lead: seven comparisons', grep='The six comparisons of Table'),
    dict(id='MZ8', expect="absent", where="S9.15: the default configuration's refusal", grep='Why the coupled system is not marched here'),
    dict(id='MZ9', expect="absent", where='abstract: moved (part-1 note, Sec. 2.2 calls eq. (9) a separate description)', grep="These comparisons exercise the formulation's two limits"),
    dict(id='MZ10', expect="absent", where='intro: replaced by the march', grep='closed form, not marched\\s*\\(Section'),
]


# ---------------------------------------------------------------------------
# 2026-09-30, REFEREE_AUDIT_8 part 3 (the owner's decision
# docs/GT_PS2_OWNER_DECISION_PAPER1_GO_SOYBEAN_ENSEMBLE_2026-09-30.md): the
# soybean size-ensemble march as the main-text benchmark (Sec. 5.4, Table 5,
# Figure 3, the abstract, Secs. 5, 7.2, 7.3, the conclusions), the single
# sphere and the candidates in S9.16, the ensemble in full and the sunflower
# run in S9.17.  Records: the march's step 3 and step 4, both placed in the
# tree (untracked until the next chain) when this block was run; the staging
# variables in PENDING_ROOTS are kept as in the blocks above.  EN1-EN* bind the printed
# numbers; EZ1-EZ* the withdrawn wording.  Every artifact behind them is
# TEMPORARILY_UNCERTIFIED_EMERGENCY_2026-09-30.
# ---------------------------------------------------------------------------
_MS3 = "../docs/evidence/paper1_faner2019_march_step3_2026-09-30"
_MS3R = "../docs/GT_PS2_FANER2019_MARCH_STEP3_2026-09-30.md"
_MS4 = "../docs/evidence/paper1_faner2019_march_step4_2026-09-30"     # scratchpad faner2019_march/ (evidence/step4)
_MS4R = "../docs/GT_PS2_FANER2019_MARCH_STEP4_2026-09-30.md"          # scratchpad faner2019_march/docs/
PENDING_ROOTS[_MS3] = "DTDC_PAPER_MARCH_STEP3_EVIDENCE"
PENDING_ROOTS[_MS3R] = "DTDC_PAPER_MARCH_STEP3_RECORD"
PENDING_ROOTS[_MS4] = "DTDC_PAPER_MARCH_STEP4_EVIDENCE"
PENDING_ROOTS[_MS4R] = "DTDC_PAPER_MARCH_STEP4_RECORD"
# the two records themselves join the corpus (a file root is read as one file)
EXTRA_CORPUS_ROOTS += [_MS3, _MS4, _MS3R, _MS4R]

_E4 = _MS4 + "/evidence/step4/"
_SUP = _E4 + "superpose_lognormal.json"
_DT = _E4 + "derived_table.md"
_SC = _E4 + "scores/lognormal_soy__sample_clock.json"


def _en(i, printed, where, record, grep, page=None):
    if record.endswith(".md"):
        grep = grep.replace(" ", r"\s+")  # a record wraps its lines anywhere
    d = dict(id=f"EN{i}", printed=printed, unit="-", where=where, record=record, grep=grep)
    if page is not None:
        d["page"] = page
    return d


_EN = [
    # --- Sec. 5.4 lead, Table 5, abstract, Secs. 5, 7, conclusions: the soybean ensemble
    (21.59, 'Sec. 5.4/Table 5/Table S13 soybean ensemble marched, crossing to 0.10, s', _SUP, r'"crossing_to_0p10_s": 21\.593478954118133'),
    (0.85, 'Abstract/Sec. 5/5.4/Table 2/5/intro/7.2/conclusions soybean ensemble over measured (0.8518)', _SUP, r'"ratio_to_measured": 0\.8517989152794639'),
    (0.852, 'Table S13/S16 soybean ensemble marched over measured', _SUP, r'"ratio_to_measured": 0\.8517989152794639'),
    (0.75, 'Abstract/5.4/Table 2/5/intro/7.2/conclusions bracket low end (uniform, number basis, 0.752)', _DT, r'\| soybean \| uniform_1_3 \| None \| 1-3 \| number \| declared \| closed_form \| hold \| 1\.446 \| 19\.06 \(0\.752\)'),
    (0.752, 'S16 uniform number basis, soybean', _DT, r'\| soybean \| uniform_1_3 \| None \| 1-3 \| number \| declared \| closed_form \| hold \| 1\.446 \| 19\.06 \(0\.752\)'),
    (1.04, 'Abstract/5.4/Table 2/5/intro/7.2/conclusions bracket high end (uniform, mass basis, 1.036)', _DT, r'\| soybean \| uniform_1_3 \| None \| 1-3 \| mass \| declared \| closed_form \| hold \| 1\.011 \| 26\.27 \(1\.036\)'),
    (1.036, 'S16 uniform mass basis, soybean', _DT, r'\| soybean \| uniform_1_3 \| None \| 1-3 \| mass \| declared \| closed_form \| hold \| 1\.011 \| 26\.27 \(1\.036\)'),
    (0.70, 'Sec. 5/5.4/Table 2/7.2/conclusions sunflower ensemble over measured (0.698)', _SUP, r'"ratio_to_measured": 0\.6979931127546416', '0.70'),
    (0.698, 'S9.17/Table S13/S16 sunflower ensemble over measured', _SUP, r'"ratio_to_measured": 0\.6979931127546416'),
    (29.66, 'S9.17/Table S13 sunflower ensemble marched, crossing to 0.10, s', _SUP, r'"crossing_to_0p10_s": 29\.664769950041034'),
    (0.5, 'Sec. 5.4/S9.17 the size statement, lower end of the range, mm', _MS4R, r'ranged from 0\.5 to 4 10-3 m'),
    (95.1, 'Sec. 5.4/S9.17 soybean share within 1 to 3 mm, percent', _MS4R, r'between 1 and 3 10-3 m was 95\.1%'),
    (97.1, 'Sec. 5.4/S9.17 sunflower share within 1 to 3 mm, percent', _MS4R, r'approximately 97\.1% of the particles'),
    (1.77, 'Sec. 5.4/S9.17 soybean mean equivalent diameter, mm', _MS4R, r'average equivalent diameter of 1\.77 10-3 m'),
    (1.81, 'S9.17 sunflower mean equivalent diameter, mm', _MS4R, r'average equivalent diameter of 1\.81 10-3 m'),
    (0.74, 'Sec. 5.4/S9.17 sphericity', _MS4R, r'sphericity shape factor \(psi\) was 0\.74'),
    (0.280, 'Sec. 5.4 sigma_ln d, soybean (0.27996)', _E4 + "derive_ensemble.json", r'"sigma_ln": 0\.27996199199685146', '0.280'),
    (0.2800, 'S9.17 sigma_ln d, soybean', _DT, r'\| 1\.50 / 2\.15 / 3\.03 \| 0\.2800 \|', '0.2800'),
    (0.2521, 'S9.17 sigma_ln d, sunflower', _DT, r'\| 1\.54 / 2\.12 / 2\.91 \| 0\.2521 \|'),
    (0.5335, 'S9.17 mu_ln d, soybean (ln mm)', _MS4R, r'sigma_ln 0\.2800 \(mu_ln 0\.5335'),
    (0.5623, 'S9.17 mu_ln d, sunflower (ln mm)', _MS4R, r'sunflower 0\.2521 \(0\.5623\)'),
    (0.79, 'Sec. 5.4/Table 5 log-normal readings, low end (0.788, the other statement)', _DT, r'\| soybean \| lognormal \| 0\.971 \| 0\.5-4 \| number \| declared \| closed_form \| hold \| 1\.176 \| 19\.96 \(0\.788\)'),
    (0.788, 'S16 the other statement, soybean', _DT, r'\| soybean \| lognormal \| 0\.971 \| 0\.5-4 \| number \| declared \| closed_form \| hold \| 1\.176 \| 19\.96 \(0\.788\)'),
    (0.726, 'S9.17/S16 the other statement, sunflower', _DT, r'\| sunflower \| lognormal \| 0\.951 \| 0\.5-4 \| number \| declared \| closed_form \| hold \| 1\.224 \| 30\.85 \(0\.726\)'),
    (0.850, 'S16 untruncated, soybean', _DT, r'\| soybean \| lognormal \| 0\.951 \| none \| number \| declared \| closed_form \| hold \| 1\.224 \| 21\.55 \(0\.850\)', '0.850'),
    (0.685, 'S16 untruncated, sunflower', _DT, r'\| sunflower \| lognormal \| 0\.971 \| none \| number \| declared \| closed_form \| hold \| 1\.181 \| 29\.11 \(0\.685\)'),
    (0.686, 'S16 mass basis, sunflower', _DT, r'\| sunflower \| lognormal \| 0\.971 \| 0\.5-4 \| mass \| declared \| closed_form \| hold \| 0\.902 \| 29\.15 \(0\.686\)'),
    (0.646, 'S16 uniform number basis, sunflower (0.65 in S9.17)', _DT, r'\| sunflower \| uniform_1_3 \| None \| 1-3 \| number \| declared \| closed_form \| hold \| 1\.411 \| 27\.44 \(0\.646\)'),
    (0.65, 'S9.17 sunflower bracket low end (0.646)', _DT, r'27\.44 \(0\.646\)'),
    (0.894, 'S16 uniform mass basis, sunflower', _DT, r'\| sunflower \| uniform_1_3 \| None \| 1-3 \| mass \| declared \| closed_form \| hold \| 0\.988 \| 38\.01 \(0\.894\)'),
    (0.89, 'S9.17 sunflower bracket high end (0.894)', _DT, r'38\.01 \(0\.894\)'),
    (0.68, 'S9.17 sunflower log-normal readings low end (0.682 primary)', _DT, r'\| 0\.17108 \(1\.178\) \| 28\.97 \| 0\.682 \|'),
    (0.524, 'Sec. 5.4 smallest class radius, mm (0.5241)', _MS4R, r'\| soybean \| 1 \| 0\.5241 \| 0\.1428 \| 18\.1 \|'),
    (1.219, 'Sec. 5.4 largest class radius, mm (1.2190)', _MS4R, r'\| soybean \| 7 \| 1\.2190 \| 0\.1429 \| 42\.2 \|'),
    (0.14579, 'Sec. 5.4/Table 5/S9.17 soybean ensemble multiplier', _SUP, r'"multiplier": 0\.14578543636088787'),
    (0.17108, 'S9.17 sunflower ensemble multiplier', _MS4R, r'\| sunflower \| 0\.14517 \| \*\*0\.17108\*\* \(1\.178\)'),
    (0.11966, 'Sec. 5.4/Table 5/S9.17 single-sphere multiplier, soybean', _MS4R, r'\| soybean \| 0\.11966 \| \*\*0\.14579\*\* \(1\.218\)'),
    (0.14517, 'S9.17 single-sphere multiplier, sunflower', _MS4R, r'\| sunflower \| 0\.14517 \| \*\*0\.17108\*\*'),
    (1.218, 'S9.17 M over M1, soybean', _MS4R, r'\*\*0\.14579\*\* \(1\.218\)'),
    (1.178, 'S9.17 M over M1, sunflower', _MS4R, r'\*\*0\.17108\*\* \(1\.178\)'),
    (0.14508, 'S9.17 initial-rate rule multiplier, soybean', _MS4R, r'\(1\.446\) \| 0\.14508 \|'),
    (0.17014, 'S9.17 initial-rate rule multiplier, sunflower', _MS4R, r'\(1\.411\) \| 0\.17014 \|'),
    (0.006, 'Sec. 5/5.4/Table 5 class 4 on 12 cells moves the soybean ensemble time, s (derived 21.5935 - 21.5875)', _MS4R, r'Class 4 at N = 12 moves the ensemble time by 0\.006 and 0\.009 s'),
    (21.588, 'S9.17 soybean ensemble time with class 4 on 12 cells', _E4 + "mesh_check.json", r'"crossing_to_0p10_s": 21\.58751328179912'),
    (21.593, 'S9.17 soybean ensemble time with class 4 on 4 cells', _E4 + "mesh_check.json", r'"crossing_to_0p10_s": 21\.593478954118133'),
    (29.665, 'S9.17 sunflower ensemble time, class 4 on 4 cells', _E4 + "mesh_check.json", r'"crossing_to_0p10_s": 29\.664769950041034'),
    (29.656, 'S9.17 sunflower ensemble time, class 4 on 12 cells', _E4 + "mesh_check.json", r'"crossing_to_0p10_s": 29\.655938857700264'),
    (21.28, 'Sec. 5.4/Table 5/S9.17 soybean ensemble, closed form, s', _DT, r'\| 0\.14579 \(1\.218\) \| 21\.28 \| 0\.840 \|'),
    (0.84, 'Sec. 5.4/Table 5 soybean ensemble closed form over measured (0.840)', _DT, r'\| 0\.14579 \(1\.218\) \| 21\.28 \| 0\.840 \|'),
    (0.840, 'S16 primary, soybean', _DT, r'\| 0\.14579 \(1\.218\) \| 21\.28 \| 0\.840 \|', '0.840'),
    (28.97, 'S9.17 sunflower ensemble, closed form, s', _DT, r'\| 0\.17108 \(1\.178\) \| 28\.97 \| 0\.682 \|'),
    (0.682, 'S9.17/S16 sunflower ensemble closed form over measured', _DT, r'\| 0\.17108 \(1\.178\) \| 28\.97 \| 0\.682 \|'),
    (21.50, 'S9.17 seven-class closed form, soybean, s', _E4 + "classes.json", r'"crossing_to_0p10_s": 21\.50161251206662', '21.50'),
    (29.30, 'S9.17 seven-class closed form, sunflower, s', _E4 + "classes.json", r'"crossing_to_0p10_s": 29\.30444352547402', '29.30'),
    (0.848, 'S16 seven-class closed form over measured, soybean', _MS4R, r'21\.50 \(0\.848\)'),
    (0.690, 'S16 seven-class closed form over measured, sunflower', _MS4R, r'29\.30 \(0\.690\)', '0.690'),
    (0.4, 'Sec. 5.4 each class\'s own falling rate adds 0.1 to 0.4 s', _MS4R, r'front recedes adds 0\.1 to 0\.4 s'),
    (87.0, 'Sec. 5.4/Table 5 ensemble mean temperature leaves the plateau by 2 K, s', _SUP, r'"mean_T_plateau_end_2K_s": 87\.00710924436292', '87.0'),
    (0.045, 'Sec. 5.4/Table 5 loading there', _SUP, r'"mean_T_plateau_end_loading": 0\.04480312655949622'),
    (89.7, 'Sec. 5.4/Table 5 measured plateau end, s', _SUP, r'"measured_plateau_end_s": 89\.71539763431348'),
    (0.063, 'Sec. 5.4/Table 5 measured plateau-end loading', _SUP, r'"measured_plateau_end_loading": 0\.06325968872941157'),
    (0.71, 'Sec. 5.4 boiling mass fraction at 0.10, soybean', _SUP, r'"boiling_mass_fraction_at_0p10": 0\.7143592690248014'),
    (0.74, 'Sec. 5.4 ensemble rate at the crossing over the measured constant rate (0.742)', _SUP, r'"rate_at_crossing_over_measured_constant_rate": 0\.7415906135450228'),
    (0.97, 'Sec. 5.4/7.2 marched ensemble mean rate over the 20 s before the crossing (0.973)', _E4 + "rpre_marched_lognormal.json", r'"r_pre_over_r_c": 0\.9729910846403282'),
    (0.973, 'Table 5/S9.17 marched r_pre over r_c, soybean', _E4 + "rpre_marched_lognormal.json", r'"r_pre_over_r_c": 0\.9729910846403282'),
    (0.974, 'S9.17 marched r_pre over r_c, sunflower', _E4 + "rpre_marched_lognormal.json", r'"r_pre_over_r_c": 0\.9738097704463864'),
    (0.199, 'Table 5 marched ensemble onset by the step-2 rule, kg/kg', _E4 + "rpre_marched_lognormal.json", r'"loading": 0\.1992360137709006'),
    (0.718, 'Table 5/S9.17 measured r_pre over r_c, soybean', _DT, r'\| soybean \| measured \| - \| 25\.35 \| 1 \| 0\.718 \(r_pre / r_c\) \| 0\.275 \|'),
    (0.72, 'Sec. 5.4/7.2 measured r_pre over r_c, soybean (0.718)', _DT, r'0\.718 \(r_pre / r_c\)'),
    (0.614, 'S9.17 measured r_pre over r_c, sunflower', _DT, r'\| sunflower \| measured \| - \| 42\.50 \| 1 \| 0\.614 \(r_pre / r_c\) \| 0\.243 \|'),
    (0.275, 'Table 5/S9.17 measured onset, soybean', _DT, r'0\.718 \(r_pre / r_c\) \| 0\.275 \|'),
    (0.243, 'S9.17 measured onset, sunflower', _DT, r'0\.614 \(r_pre / r_c\) \| 0\.243 \|'),
    (0.946, 'Table 5/S9.17 closed-form ensemble r_pre, soybean', _DT, r'\| 0\.840 \| 0\.833; 0\.946 \| 0\.215; 0\.315 \| 0\.41 \| 0\.31 \|'),
    (0.943, 'S9.17 closed-form ensemble r_pre, sunflower', _DT, r'\| 0\.682 \| 0\.862; 0\.943 \| 0\.216; 0\.305 \| 0\.56 \| 0\.38 \|'),
    (0.215, 'Table 5/S9.17 closed-form onset, soybean', _DT, r'0\.833; 0\.946 \| 0\.215; 0\.315'),
    (0.216, 'S9.17 closed-form onset, sunflower', _DT, r'0\.862; 0\.943 \| 0\.216; 0\.305'),
    (0.41, 'Table 5 closed-form boiling mass at the measured plateau-end loading, soybean', _DT, r'0\.215; 0\.315 \| 0\.41 \| 0\.31 \|'),
    (0.31, 'S9.17 closed-form boiling mass at the measured plateau-end time, soybean', _DT, r'0\.215; 0\.315 \| 0\.41 \| 0\.31 \|'),
    (0.38, 'S9.17 closed-form boiling mass at the measured plateau-end time, sunflower', _DT, r'0\.216; 0\.305 \| 0\.56 \| 0\.38 \|'),
    (0.43, 'Table 5 marched boiling mass at the measured plateau-end loading, soybean', _SUP, r'"boiling_mass_fraction_at_measured_plateau_end_loading": 0\.42868051751655206'),
    (4.7, 'Sec. 5.4 the ensemble crosses 0.20 before the measurement, s (4.70)', _SUP, r'"shift_s": 4\.695825474706858'),
    (4.70, 'S9.17 the same, soybean', _MS4R, r'4\.70 s and\s+4\.59 s', '4.70'),
    (4.59, 'S9.17 the same, sunflower', _MS4R, r'4\.70 s and\s+4\.59 s'),
    (0.092, 'Sec. 5.4/Table 2/5 RMS over the rule window, soybean ensemble, of the span', _SUP, r'"rms_over_span": 0\.09249449164217872'),
    (0.062, 'Sec. 5.4 crossing-matched RMS, soybean ensemble', _SUP, r'"rms_over_span": 0\.06222668475804202'),
    (150, 'Sec. 5.4 last sample outside the band, s (150.116)', _SC, r'"t_s": 150\.116'),
    (160, 'Sec. 5.4/Table 5 first in-band sample, s', _SC, r'"t_s": 160\.0'),
    (230, 'Sec. 5.4/Table 5 last in-band sample, s (230.116)', _SC, r'"t_s": 230\.116'),
    (0.0205, 'Sec. 5.4/S9.16 sorbed share, the model isotherm (S1)', _MS3R, r'\| S1, the model\'s \(Cardarelli and Crapiste 1996\) \| 0\.0205 \|'),
    (0.0417, 'Sec. 5.4/S9.16 sorbed share, oil-arm upper bound (S5a)', _MS3R, r'\| S5a, S1 plus the \(2\.24\) oil arm, upper bound \| 0\.0417 \|'),
    (0.0259, 'S9.16 sorbed share, Tabla 2.6 soybean (S3a)', _MS3R, r'\| S3a, Cardarelli 1998 Tabla 2\.6 \| 0\.0259 \|'),
    (0.0096, 'S9.16 sorbed share, sunflower GAB (S2)', _MS3R, r'\| S2, the sunflower meal\'s own 1996 GAB \| 0\.0096 \|'),
    (0.0163, 'S9.16 sorbed share, Tabla 2.6 sunflower (S3b)', _MS3R, r'\| S3b, Tabla 2\.6 \| 0\.0163 \|'),
    (0.0292, 'S9.16 sorbed share, sunflower oil-arm upper bound (S5b)', _MS3R, r'\| S5b, S2 plus the oil arm, upper bound \| 0\.0292 \|'),
    (0.558, 'Sec. 5.4/Table S13 partition, soybean ratio low end', _MS3R, r'stays at 0\.558 to\s+0\.563 \(soybean\) and 0\.478 to 0\.481 \(sunflower\)'),
    (0.563, 'Sec. 5.4/Table S13 partition, soybean ratio high end', _MS3R, r'stays at 0\.558 to\s+0\.563 \(soybean\)'),
    (0.478, 'Table S13 partition, sunflower ratio low end', _MS3R, r'0\.478 to 0\.481 \(sunflower\)'),
    (0.481, 'Table S13 partition, sunflower ratio high end', _MS3R, r'0\.478 to 0\.481 \(sunflower\)'),
    (14.27, 'Table S13 partition, soybean time high end (M07)', _MS3R, r'\| M07 \| soybean \| 4 \| S5a, 0\.0417 \| 338 \| 3 \| 14\.273 \(0\.563\)'),
    (20.30, 'Table S13 partition, sunflower time low end (M04)', _MS3R, r'\| M04 \| sunflower \| 12 \| S2, 0\.0096 \| 574 \| 11 \| 20\.300 \(0\.478\)', '20.30'),
    (20.43, 'Table S13 partition, sunflower time high end (M10)', _MS3R, r'\| M10 \| sunflower \| 4 \| S5b, 0\.0292 \| 509 \| 3 \| 20\.426 \(0\.481\)'),
    (0.12, 'S9.16 the partition moves the time by at most, s', _MS3R, r'moves by at most 0\.12 s \(0\.9 percent\)'),
    (26.8, 'S9.16 core end after the crossing, low end, s', _MS3R, r'26\.8 to\s+43\.3 s after the crossing'),
    (43.3, 'S9.16 core end after the crossing, high end, s', _MS3R, r'26\.8 to\s+43\.3 s after the crossing'),
    (0.136, 'S9.16 the value the data would need, soybean, low end (reader only)', _MS3R, r'X_s = 0\.136 to 0\.169\s+\(soybean\)'),
    (0.169, 'S9.16 the value the data would need, soybean, high end (reader only)', _MS3R, r'X_s = 0\.136 to 0\.169\s+\(soybean\)'),
    (3.3, 'S9.16 the needed value over the sources, low end', _MS3R, r'That is 3\.3 to\s+16 times the source values'),
    (16, 'S9.16 the needed value over the sources, high end', _MS3R, r'That is 3\.3 to\s+16 times the source values'),
    (17, 'S9.16 gas cooling across the layer, high end, K', _MS3R, r'cools the gas by 7 to 17 K at the constant rate'),
    (1.125, 'Table S13 (g) soybean coefficient ratio low end', _MS3R, r'\| 0\.214 to 0\.336 \| 1\.125 to 1\.219 \|'),
    (1.219, 'Table S13 (g) soybean coefficient ratio high end', _MS3R, r'\| 0\.214 to 0\.336 \| 1\.125 to 1\.219 \|'),
    (1.129, 'Table S13 (g) sunflower coefficient ratio low end', _MS3R, r'\| 0\.219 to 0\.348 \| 1\.129 to 1\.229 \|'),
    (1.229, 'Table S13 (g) sunflower coefficient ratio high end', _MS3R, r'\| 0\.219 to 0\.348 \| 1\.129 to 1\.229 \|'),
    (2.0, 'Sec. 5.4/S9.16/Table S13 vapor-flow resistance, soybean, percent of the driving force', _MS3R, r'1/51\.3 = 2\.0 percent \(soybean\)', '2.0'),
    (3.2, 'S9.16/Table S13 vapor-flow resistance, sunflower, percent', _MS3R, r'1/31\.3 = 3\.2 percent'),
    (14.5, 'S9.16 front elevation the falling rate would need, soybean, K', _MS3R, r'sit 14\.5 K and 12\.1 K above T_b'),
    (12.1, 'S9.16 front elevation the falling rate would need, sunflower, K', _MS3R, r'sit 14\.5 K and 12\.1 K above T_b'),
    (15, 'Sec. 5.4/7.2/conclusions the time left unexplained, soybean, percent (1 - 0.852)', _MS4R, r'it leaves 15 and 30\s+percent of the time'),
    (0.0011, 'S9.17 untruncated share above 4 mm, soybean', _MS4R, r'with 0\.0011 and 0\.0005 of the particles above 4 mm'),
    (0.0005, 'S9.17 untruncated share above 4 mm, sunflower', _MS4R, r'with 0\.0011 and 0\.0005 of the particles above 4 mm'),
    (3.03, 'S9.17 mass-weighted d90, soybean, mm', _MS4R, r'1\.50 / 2\.15 / 3\.03 mm \(mass mean 2\.21 mm\)'),
    (2.21, 'S9.17 mass mean, soybean, mm', _MS4R, r'\(mass mean 2\.21 mm\)'),
    (2.91, 'S9.17 mass-weighted d90, sunflower, mm', _MS4R, r'sunflower 1\.54 / 2\.12 / 2\.91 mm \(2\.18 mm\)'),
    (2.18, 'S9.17 mass mean, sunflower, mm', _MS4R, r'sunflower 1\.54 / 2\.12 / 2\.91 mm \(2\.18 mm\)'),
    (2.00, 'S9.17 uniform mean, mm', _MS4R, r'its mean is 2\.00 mm', '2.00'),
    (18.1, 'S9.17 smallest marched Reynolds number', _MS4R, r'Every marched class lies at Re 18\.1 to 43\.0'),
    (43.0, 'S9.17 largest marched Reynolds number', _MS4R, r'Every marched class lies at Re 18\.1 to 43\.0', '43.0'),
    (11.8, 'S9.17 the 0.5 mm particle, soybean Reynolds number', _MS4R, r'sits at Re 11\.8 \(soybean\) and 13\.4'),
    (13.4, 'S9.17 the 0.5 mm particle, sunflower Reynolds number', _MS4R, r'sits at Re 11\.8 \(soybean\) and 13\.4'),
    (0.0016, 'S9.17 relaxation-to-forcing ratio, low end', _MS4R, r'ratio at every marched class is 0\.0016 to 0\.0055'),
    (0.0055, 'S9.17 relaxation-to-forcing ratio, high end', _MS4R, r'ratio at every marched class is 0\.0016 to 0\.0055'),
    (0.8725, 'S9.17 Ackermann factor, soybean', _MS4R, r'Ackermann factor \(0\.8725 soybean, 0\.9184 sunflower\)'),
    (0.9184, 'S9.17 Ackermann factor, sunflower', _MS4R, r'Ackermann factor \(0\.8725 soybean, 0\.9184 sunflower\)'),
    (0.2388, 'S9.17 the window rule threshold, kg/kg', _MS4R, r'the samples with X > 0\.2388'),
    (20.6, 'S9.17 the first class finishes inside the window, soybean, s', _MS4R, r'the first at 20\.6 s and 33\.5 s'),
    (33.5, 'S9.17 the first class finishes inside the window, sunflower, s', _MS4R, r'the first at 20\.6 s and 33\.5 s'),
    (2.278e-10, 'S9.17 start probe, 2-cell root', _MS4R, r'the N = 2 root at 2\.278e-10, 1\.905e-10 and 1\.470e-09', '2.278'),
    (1.905e-10, 'S9.17 start probe, 2-cell root', _MS4R, r'2\.278e-10, 1\.905e-10 and 1\.470e-09', '1.905'),
    (1.470e-09, 'S9.17 start probe, 2-cell root', _MS4R, r'1\.905e-10 and 1\.470e-09', '1.470'),
    (1.405e-10, 'S9.17 start probe, soybean class 7 bootstrap', _MS4R, r'N = 4 bootstrap at 1\.405e-10', '1.405'),
    (1.855e-09, 'S9.17 bootstrap refusal, soybean class 2', _MS4R, r'\(1\.855e-09,', '1.855'),
    (1.173e-09, 'S9.17 bootstrap refusal, soybean class 5', _MS4R, r'1\.173e-09, 7\.976e-10, 6\.002e-10', '1.173'),
    (7.976e-10, 'S9.17 bootstrap refusal, soybean class 6', _MS4R, r'1\.173e-09, 7\.976e-10, 6\.002e-10', '7.976'),
    (6.002e-10, 'S9.17 bootstrap refusal, soybean class 7', _MS4R, r'1\.173e-09, 7\.976e-10, 6\.002e-10', '6.002'),
    (5.382e-10, 'S9.17 sunflower class 5, 2-cell root', _MS4R, r'Sunflower class 5 refused at the N = 2 root \(5\.382e-10\)', '5.382'),
    (1.140, 'S9.17/Table S15 r4 over r5, sunflower', _MS4R, r'r_4 / r_5 = 1\.140', '1.140'),
    (0.005, 'S9.17 class 5 carried by class 6 instead, s', _MS4R, r'the ensemble time to 0\.10 moves by 0\.005 s'),
    (39.075, 'S9.17/Table S15 sunflower class 7 refusal time on 12 cells, s', _MS4R, r'class 7 at N = 12 refused at 39\.075 s \(X 0\.0960, front 0\.751 R'),
    (0.0960, 'S9.17/Table S15 sunflower class 7 loading there', _MS4R, r'refused at 39\.075 s \(X 0\.0960', '0.0960'),
    (0.751, 'S9.17 sunflower class 7 front there, R', _MS4R, r'\(X 0\.0960, front 0\.751 R'),
    (9.99965e-11, 'S9.17/Table S15 sunflower class 7 cumulative ledger', _MS4R, r'cumulative\s+ledger 9\.99965e-11', '9.99965'),
    (149.0, 'S9.17 sunflower coverage end on the sample clock, s', _MS4R, r'clock only to 149\.0 s', '149.0'),
    (9.72e-11, 'S9.17 largest scaled residual of any class', _MS4R, r'at most 9\.72e-11 scaled', '9.72'),
    (341.90336, 'S9.17 interface on every stride, K', _MS4R, r'the interface at\s+341\.90336 K'),
    (1.7e-15, 'S9.17 cumulative water ledger at most', _MS4R, r'the cumulative water ledger at most 1\.7e-15', '1.7'),
    (1.1e-17, 'S9.17 ledger margin, low end', _MS4R, r'the 1e-10 gate by 1\.1e-17 to 8\.0e-16', '1.1'),
    (8.0e-16, 'S9.17 ledger margin, high end', _MS4R, r'the 1e-10 gate by 1\.1e-17 to 8\.0e-16', '8.0'),
    (32.40, 'S9.17 soybean class 4 on 12 cells reaches the center cell, s', _E4 + "mesh_check.json", r'"march_end_s": 32\.399999999999814', '32.40'),
    (1.06e-11, 'S9.17 soybean class 4 on 12 cells, cumulative ledger', _E4 + "mesh_check.json", r'"max_cum_hexane_ledger": 1\.059413704197396e-11', '1.06'),
    (0.018, 'S9.17 soybean class 4 time to 0.10 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.1": -0\.01841261863848942'),
    (0.078, 'S9.17 soybean class 4 time to 0.03 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.03": -0\.07777025705449958'),
    (0.039, 'S9.17 soybean class 4 time to 0.05 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.05": -0\.03948628703286161'),
    (0.0007, 'S9.17 soybean class 4 time to 0.15 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.15": -0\.000670243713230434'),
    (40.80, 'S9.17 sunflower class 4 on 12 cells refuses, s', _MS4R, r'\*\*refused\*\* at 40\.80 s \(X 0\.0238, front 0\.336 R, 7 faces\)', '40.80'),
    (0.0238, 'S9.17 sunflower class 4 on 12 cells refusal loading', _MS4R, r'at 40\.80 s \(X 0\.0238, front 0\.336 R'),
    (0.336, 'S9.17 sunflower class 4 on 12 cells refusal front, R', _MS4R, r'\(X 0\.0238, front 0\.336 R'),
    (0.0008, 'S9.17 sunflower class 4 time to 0.15 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.15": -0\.0008241200047081065'),
    (0.024, 'S9.17 sunflower class 4 time to 0.10 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.1": -0\.02358859735037555'),
    (0.054, 'S9.17 sunflower class 4 time to 0.05 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.05": -0\.05431233754905662'),
    (0.106, 'S9.17 sunflower class 4 time to 0.03 moves, s', _E4 + "mesh_check.json", r'"diff_N12_minus_N4_t_to_0\.03": -0\.10564151139844569'),
    (0.012, 'S9.17 the tails or the basis move the log-normal result by at most', _MS4R, r'by at most 0\.012 with the tails'),
    (1.67, 'S9.17 a third of the mass below, uniform mass basis, mm', _MS4R, r'puts 1/3 of the mass below\s+1\.67 mm'),
    (1.60, 'S9.17 a seventh of the mass below, soybean log-normal, mm', _MS4R, r'1/7 of the mass below 1\.60 mm \(soybean\) and 1\.62 mm', '1.60'),
    (1.62, 'S9.17 a seventh of the mass below, sunflower log-normal, mm', _MS4R, r'1\.60 mm \(soybean\) and 1\.62 mm'),
    (0.831, 'S9.17/S16 march-shaped classes, soybean', _MS4R, r'gives 0\.831 and 0\.678'),
    (0.678, 'S9.17/S16 march-shaped classes, sunflower', _MS4R, r'gives 0\.831 and 0\.678'),
    (0.813, 'S16 Crank tail, soybean, low end', _MS4R, r'ratios 0\.813 to 0\.820 soybean, 0\.662 to 0\.669 sunflower'),
    (0.820, 'S16 Crank tail, soybean, high end', _MS4R, r'ratios 0\.813 to 0\.820 soybean', '0.820'),
    (0.662, 'S16 Crank tail, sunflower, low end', _MS4R, r'0\.662 to 0\.669 sunflower'),
    (0.669, 'S16 Crank tail, sunflower, high end', _MS4R, r'0\.662 to 0\.669 sunflower'),
    (0.8, 'S9.17 the Crank tail shortens the time by at most, s', _MS4R, r'shortens the time to 0\.10 by 0\.5 to 0\.8 s'),
    (0.543, 'S16 single sphere closed form, soybean', _DT, r'\| soybean \| single declared sphere, closed form \| 0\.11966 \(1\) \| 13\.77 \| 0\.543 \|'),
    (0.466, 'S16 single sphere closed form, sunflower', _DT, r'\| sunflower \| single declared sphere, closed form \| 0\.14517 \(1\) \| 19\.81 \| 0\.466 \|'),
    (0.240, 'S9.17 RMS over the covered window, sunflower ensemble', _SUP, r'"rms_over_span": 0\.23976791601946934', '0.240'),
    (0.167, 'S9.17 crossing-matched RMS, sunflower ensemble', _SUP, r'"rms_over_span": 0\.16715037880338834'),
    (0.870, 'S9.17 sunflower ensemble rate at the crossing', _SUP, r'"rate_at_crossing_over_measured_constant_rate": 0\.8695850540756426', '0.870'),
    (142.2, 'S9.17 sunflower ensemble plateau end, s', _SUP, r'"mean_T_plateau_end_2K_s": 142\.18670406515898'),
    (0.030, 'S9.17 sunflower ensemble plateau-end loading', _SUP, r'"mean_T_plateau_end_loading": 0\.030106871958711943', '0.030'),
    (117.4, 'S9.17 measured plateau end, sunflower, s', _SUP, r'"measured_plateau_end_s": 117\.36475679816168'),
    (0.085, 'S9.17 measured plateau-end loading, sunflower', _SUP, r'"measured_plateau_end_loading": 0\.08540557964660069'),
    (0.1428, 'S9.17/Table S15 class mass weight, low end', _MS4R, r'\| soybean \| 1 \| 0\.5241 \| 0\.1428 \|'),
    (0.1430, 'S9.17/Table S15 class mass weight, high end (sunflower class 5)', _MS4R, r'\| sunflower \| 5 \| 0\.9266 \| 0\.1430 \|', '0.1430'),
    (0.74, 'S9.17 sunflower log-normal readings, high end (0.739, the other statement on the mass basis)', _DT, r'\| sunflower \| lognormal \| 0\.951 \| 0\.5-4 \| mass \| declared \| closed_form \| hold \| 0\.884 \| 31\.42 \(0\.739\)'),
    (0.666, 'Sec. 5.4 the declared sphere, soybean, mm (self-check W1)', _MS4R, r'0\.666 mm \(soybean\) and\s+0\.7215 mm \(sunflower\)'),
]

# Table S15 (the classes): every row of the record's per-class table, one
# assertion per printed cell that is not already bound above.
_ROWS = {
    ("soybean", 1): ("0.5241", "18.1", "0.01243", "23.55", "217", "16.275", "0.0185", "16.65", "0.0157", "9.16e-11", "1.6543561e-11"),
    ("soybean", 2): ("0.6474", "22.4", "0.00921", "31.79", "297", "22.275", "0.0179", "22.80", "0.0153", "5.67e-11", "1.8908571e-11"),
    ("soybean", 3): ("0.7296", "25.2", "0.00777", "37.68", "354", "26.550", "0.0178", "27.23", "0.0151", "8.47e-11", "6.6507671e-11"),
    ("soybean", 4): ("0.8079", "27.9", "0.00672", "43.56", "412", "30.900", "0.0176", "31.70", "0.0150", "9.72e-11", "9.9999204e-11"),
    ("soybean", 5): ("0.8941", "30.9", "0.00582", "50.31", "478", "35.850", "0.0176", "36.84", "0.0148", "7.06e-11", "9.9999953e-11"),
    ("soybean", 6): ("1.0053", "34.8", "0.00493", "59.42", "569", "42.675", "0.0174", "43.89", "0.0146", "9.45e-11", "9.9999989e-11"),
    ("soybean", 7): ("1.2190", "42.2", "0.00375", "78.15", "758", "56.850", "0.0170", "58.56", "0.0143", "9.39e-11", "9.9999970e-11"),
    ("sunflower", 1): ("0.5712", "19.8", "0.00829", "36.44", "322", "24.150", "0.0192", "24.73", "0.0165", "5.61e-11", "1.3460115e-11"),
    ("sunflower", 2): ("0.6915", "23.9", "0.00632", "47.80", "427", "32.025", "0.0188", "32.83", "0.0162", "9.64e-11", "9.9868381e-11"),
    ("sunflower", 3): ("0.7704", "26.7", "0.00542", "55.75", "524", "39.300", "0.0150", "39.34", "0.0149", "9.32e-11", "2.5818322e-11"),
    ("sunflower", 4): ("0.8450", "29.2", "0.00475", "63.56", "574", "43.050", "0.0185", "44.19", "0.0159", "9.44e-11", "9.9999942e-11"),
    ("sunflower", 5): ("0.9266", "32.1", "0.00417", "72.46", None, "49.075", "0.0185", "50.38", "0.0159", None, None),
    ("sunflower", 6): ("1.0317", "35.7", "0.00358", "84.41", "808", "60.600", "0.0145", "60.69", "0.0144", "6.09e-11", "9.9929407e-11"),
    ("sunflower", 7): ("1.2427", "43.0", "0.00275", "109.97", "521", "39.075", "0.0960", None, None, "9.41e-11", "9.9996462e-11"),
}
# what the page prints for the record's longer ledger values
_PAGE_LEDGER = {"1.6543561e-11": "1.65", "1.8908571e-11": "1.89", "6.6507671e-11": "6.65",
                "1.3460115e-11": "1.35", "9.9868381e-11": "9.9868", "2.5818322e-11": "2.58",
                "9.9929407e-11": "9.993", "9.9996462e-11": "9.99965"}
_n = len(_EN)
for (trace, k), cells in _ROWS.items():
    rad = cells[0]
    for pos, tok in enumerate(cells):
        if tok is None:
            continue
        page = _PAGE_LEDGER.get(tok, tok)
        if "e" in page:
            page = page.split("e")[0]
        grep = r"\| " + trace + r" \| " + str(k) + r" \| " + re.escape(rad) + r" \|"
        if pos:
            grep += r"[^\n]*" + re.escape(tok)
        _EN.append((float(tok), f"Table S15 {trace} class {k}, column {pos + 1}", _MS4R, grep, page))

ASSERTIONS += [_en(i + 1, *row) for i, row in enumerate(_EN)]

ASSERTIONS += [
    dict(id='EZ1', expect="absent", where='Sec. 7.2: the candidates are now measured (steps 3 and 4)', grep='Three candidates, none measured'),
    dict(id='EZ2', expect="absent", where='S9.16: the candidates are now measured', grep='Three candidates are named and not measured'),
    dict(id='EZ3', expect="absent", where='abstract: the benchmark is the soybean ensemble', grep='mesh-invariant factor of two still unexplained'),
    dict(id='EZ4', expect="absent", where='conclusions: the factor is carried mostly by the size spread', grep='none of the five candidates\\s+measured explains'),
    dict(id='EZ5', expect="absent", where='conclusions: the candidates measured', grep='evaporation zone are the candidates for its factor of two'),
    dict(id='EZ6', expect="absent", where='S9.16 retitled for the single sphere', grep='The moving-front march of the journal runs in full'),
    dict(id='EZ7', expect="absent", where='Sec. 5 lead and intro: the benchmark is 0.85', grep='dries in\\s+about half the measured time'),
    dict(id='EZ8', expect="absent", where='Sec. 5.4 retitled for the soybean run', grep='The moving-front system marched through the journal runs\\}'),
]

# ---------------------------------------------------------------------------
# 2026-09-30, REFEREE_AUDIT_8 part 4 (section 8.O): the process vocabulary is
# absent from the printed text of the main text and both supplements (the
# owner's instruction: a process-modeling paper, not a software-development
# one).  printed_only: comments and \url/\label/\ref arguments removed first.
# ---------------------------------------------------------------------------
_PZ = [
    ("uncertified / non-certified / temporarily uncertified", r"(?i)\b(un|non-)certified\b"),
    ("certified, certification, certify, re-certification (not 'certificate')", r"(?i)certif(y|ied|ies|ying|ication)"),
    ("emergency (lane, contract, records)", r"(?i)\bemergency\b"),
    ("contract change(s)", r"(?i)\bcontract\s+changes?\b"),
    ("debt list", r"(?i)\bdebts?\b"),
    ("gate in the numerical sense (gate, gates, gated, gate-complete, subgate)", r"(?i)(\bgat(e|es|ed|ing)\b|subgate)"),
    ("knob", r"(?i)\bknobs?\b"),
    ("flag, flags, flagged", r"(?i)\bflag(s|ged|ging)?\b"),
    ("frozen (gate, contract, value, set, isotherm, branch)", r"(?i)\bfrozen\b"),
    ("typed (refusal, outcome, construction, inapplicability)", r"(?i)\btyped\b"),
    ("committed (stride, read, record, march)", r"(?i)\bcommit(s|ted|ting)?\b"),
    ("owner", r"(?i)\bowner('s)?\b"),
    ("ruling, ruled", r"(?i)\brul(ing|ings|ed)\b"),
    ("sidecar", r"(?i)\bsidecars?\b"),
    ("farm", r"(?i)\bfarm\b"),
    ("fail-closed", r"(?i)fail[- ]closed"),
    ("licence, license (to subdivide)", r"(?i)\blicen[cs](e|es|ed|ing)\b"),
    ("re-earned", r"(?i)re-earn"),
    ("unscored", r"(?i)\bunscored\b"),
    ("shipped (kernel, scheme)", r"(?i)\bshipped\b"),
    ("production (path, route, gate, operator)", r"(?i)\bproduction\s+(path|route|gate|residual|operator)"),
    ("contract in the process sense (residual, failed-leaf, thermal, campaign, preselection, evidence, event-chart, marching, A/B/C, own)",
     r"(?i)(residual|failed-leaf|thermal|campaign|preselection|evidence|event-chart|marching|own)\s+contracts?\b|\bcontracts?\s+(gates|[ABC]\b|itself|stays|family)|\bthe\s+contract('s)?\b"),
    ("the pin, regression pin, source-pinned", r"(?i)\b(the|regression)\s+pin\b|source-pinned"),
    ("authority (property authority, the certified authority)", r"(?i)\bauthority\b"),
    ("switch in the software sense (behind a switch, switches of the feasibility marches)", r"(?i)behind\s+a\s+switch|face-tangent\s+switches|every\s+switch"),
    ("product path, stock tangent solve, scorer, execution adapter", r"(?i)product\s+path|stock\s+tangent|\bscorer\b|execution\s+adapter"),
    ("the code state called a tree", r"(?i)\b(current|that|source)\s+tree\b"),
]
ASSERTIONS += [
    dict(id=f"PZ{i + 1}", expect="absent", printed_only=True,
         where=f"process vocabulary absent from the printed text: {label}", grep=pat)
    for i, (label, pat) in enumerate(_PZ)
]

# ---------------------------------------------------------------------------
# 2026-10-01, REFEREE_AUDIT_8 part 6 (section 8.Q, the whole-manuscript referee
# read): the numbers the pass prints.  Landed records: the step-4
# superposition (_SUP) and the two-regime derivation (_TR).  The part-6
# derived record (derived/plateau_rule.json and derived/part6_derived.json of
# the pass's folder) is staged through DTDC_PAPER_PART6_DERIVED until it lands
# at _P6.
# ---------------------------------------------------------------------------
_P6 = "../docs/evidence/paper1_referee_audit8_part6_2026-10-01"
PENDING_ROOTS[_P6] = "DTDC_PAPER_PART6_DERIVED"
EXTRA_CORPUS_ROOTS += [_P6]
_P6D = _P6 + "/part6_derived.json"
_P6R = _P6 + "/plateau_rule.json"
_TRJ = _TR + "/two_regime.json"
_P6L = [
    # M1: the sample clock beside the crossing-matched ratio
    (41.17, 'Table 5/S9.17 ensemble crossing of 0.20 on the sample clock, s', _SUP, r'"ensemble_crossing_0p20_s": 41\.167954798392195'),
    (41.2, 'Sec. 5.4 ensemble crossing of 0.20 on the sample clock, s', _SUP, r'"ensemble_crossing_0p20_s": 41\.167954798392195'),
    (45.86, 'Table 5/S9.17 measured crossing of 0.20 on the sample clock, s', _SUP, r'"measured_crossing_s": 45\.86378027309905'),
    (45.9, 'Sec. 5.4 measured crossing on the sample clock, s', _SUP, r'"measured_crossing_s": 45\.86378027309905'),
    (62.8, 'Sec. 5.4 ensemble reaches 0.10 on the sample clock, s', _SUP, r'"ensemble_0p10_s": 62\.76143375251033'),
    (71.2, 'Sec. 5.4 measured curve reaches 0.10 on the sample clock, s', _SUP, r'"measured_time_s": 71\.21422222222222'),
    (16.90, 'Table 5 ensemble 0.10 after the measured crossing, sample clock, s', _SUP, r'"marched_crossing_to_level_s": 16\.897653479411275', '16.90'),
    (16.9, 'Sec. 5.4/7.2 ensemble 0.10 after the measured crossing, s', _SUP, r'"marched_crossing_to_level_s": 16\.897653479411275'),
    (21.59, 'Sec. 5.4/Table 5 the headline on the crossing-matched clock (scores.crossing_matched), s', _SUP, r'"marched_crossing_to_level_s": 21\.593478954118126'),
    (25.35, 'Sec. 5.4/Table 5 the measured time from the measured crossing (scores.sample_clock and crossing_matched), s', _SUP, r'"measured_crossing_to_level_s": 25\.350441949123166'),
    (0.67, 'Table 2/5, Sec. 5.4/7.2 sample-clock ratio 16.90/25.35 (derived)', _P6D, r'"sample_clock_ratio_to_measured_25p35": 0\.6665624809746459'),
    (60.47, 'S9.17 sunflower ensemble crossing on the sample clock, s', _P6R, r'"crossing_0p20_s": 60\.4695208037018'),
    (65.06, 'S9.17 sunflower measured crossing on the sample clock, s', _P6R, r'"measured_crossing_s": 65\.05777158948685'),
    # M3: the rate exponent and the first class's core end
    (1.42, 'Sec. 5.4 rate per unit mass across the classes as R^-1.42 (derived from Table S15)', _P6D, r'"rate_per_mass_scales_as_R_to_minus": 1\.4210000000000003'),
    (40, 'Sec. 5.4 the first class reaches its remainder about 40 s after the first sample (23.55 + 16.65, derived)', _P6D, r'"sum_s": 40\.203263288775936'),
    # M5: the plateau by one rule
    (69.3, 'Sec. 5.4 measured early plateau, soybean, C (69.33)', _TRJ, r'"plateau_mean_t_le_30s_C": 69\.32933333333334'),
    (69.33, 'S9.15/S9.16/Table S14 measured early plateau, soybean, C', _TRJ, r'"plateau_mean_t_le_30s_C": 69\.32933333333334'),
    (68.35, 'S9.15/S9.16/Table S14 measured early plateau, sunflower, C', _TRJ, r'"plateau_mean_t_le_30s_C": 68\.352'),
    (60.9, 'Sec. 5.4/S9.17 ensemble mean temperature 1 K above its early plateau, sample clock, s', _P6R, r'"end_plus_1K_s": 60\.88336132736254'),
    (15.0, 'Table 2/5, Sec. 5.4 the same after the measured crossing, s', _P6R, r'"end_plus_1K_minus_measured_crossing_s": 15\.01958105426349'),
    (0.107, 'Table 5, Sec. 5.4/S9.17 loading there', _P6R, r'"X_at_end_plus_1K": 0\.10683335059862901'),
    (87.1, 'Sec. 5.4/S9.17 thermocouple 1 K above its early plateau, sample clock, s', _TRJ, r'"t_above_plateau_plus_1\.0K": 87\.08804908246739'),
    (0.067, 'Table 5, Secs. 5.4/7.2, S9.15-S9.17 measured plateau-end loading by the rule, soybean', _TRJ, r'"X_at_above_plateau_plus_1\.0K": 0\.06746227681163401'),
    (0.095, 'Sec. 7.2, S9.15-S9.17 measured plateau-end loading by the rule, sunflower', _TRJ, r'"X_at_above_plateau_plus_1\.0K": 0\.09464818704236118'),
    (26, 'Secs. 5.4/7.2 the ensemble leaves the plateau at least 26 s before the thermocouple (derived)', _P6D, r'"measured_minus_ensemble_end_1K_s": 26\.20468775510485'),
    (85.7, 'S9.17 ensemble 2 K end on the early-plateau base, s', _P6R, r'"end_plus_2K_s": 85\.68909578186263'),
    (90.4, 'S9.17 measured 2 K end on the early-plateau base, s', _TRJ, r'"t_above_plateau_plus_2\.0K": 90\.43427946053505'),
    (25, 'S9.17 the ensemble 1 K and 2 K ends 25 s apart (derived)', _P6D, r'"ensemble_2K_minus_1K_s": 24\.80573445450009'),
    (3, 'S9.17 the measured 1 K and 2 K ends 3 s apart (derived)', _P6D, r'"measured_2K_minus_1K_s": 3\.346230378067659'),
    (93.4, 'S9.17 sunflower ensemble 1 K end, s', _P6R, r'"end_plus_1K_s": 93\.35095712602062'),
    (0.092, 'S9.17 sunflower ensemble loading at its 1 K end', _P6R, r'"X_at_end_plus_1K": 0\.09189576225178712'),
    (110.7, 'S9.17 sunflower thermocouple 1 K end, s', _TRJ, r'"t_above_plateau_plus_1\.0K": 110\.65225967062433'),
    # M8: the dry-shell limit on the declared spheres
    (306, 'S9.16 dry-shell limit on the declared soybean sphere, s (derived, Deff/R^2 similarity)', _P6D, r'"soybean_declared_s": 306\.37913243320884'),
    (351, 'S9.16 dry-shell limit on the declared sunflower sphere, s (derived)', _P6D, r'"sunflower_declared_s": 350\.92964090778514'),
    (12, 'Secs. 5.4/8, S9.16 dry-shell limit over the measured time, soybean declared sphere (derived)', _P6D, r'"soybean_over_measured_25p35": 12\.085961831684767'),
    (8, 'S9.16 dry-shell limit over the measured time, sunflower declared sphere (derived)', _P6D, r'"sunflower_over_measured_42p50": 8\.25716802135965'),
]


def _p6(i, printed, where, record, grep, page=None):
    d = dict(id=f"P6_{i}", printed=printed, unit="-", where=where, record=record, grep=grep)
    if page is not None:
        d["page"] = page
    if isinstance(printed, int) and printed < 100:
        d["word"] = True  # a small integer: the page check cannot tell it from other text
    return d


ASSERTIONS += [_p6(i + 1, *row) for i, row in enumerate(_P6L)]

_MS = MAIN_SOURCES
_MSS = _MS + SUPP_SOURCES + ARCHIVED_SUPP_SOURCES
_P6ZL = [
    ("M2: 'without fitted parameters' for the march", r"(?i)without\s+fitted\s+parameters", _MSS),
    ("M5/claim 5: 'holds the boiling plateau'", r"(?i)holds\s+the\s+boiling\s+plateau", _MS),
    ("M5/claim 5: 'holds the measured plateau temperature'", r"(?i)holds\s+the\s+measured\s+plateau", _MS),
    ("M5: the 2 K rule as the main-text reading", r"by\s+the\s+measured\s+rule\s+of\s+a", _MS),
    ("M1: 'from the measured crossing of the critical loading in 21.59'", r"from\s+the\s+measured\s+crossing\s+of\s+the\s+critical\s+loading\s+in", _MS),
    ("M1: Table 5 caption 'Times from the measured crossing'", r"in\s+closed\s+form\s+and\s+marched\.\s+Times\s+from\s+the\s+measured", _MS),
    ("M4a: 'monodisperse sphere'", r"(?i)monodisperse\s+sphere", _MS),
    ("M4c: 'no size-distribution robustness ... anywhere'", r"no\s+size-distribution\s+robustness", _MSS),
    ("M4d: 'no rapeseed or sunflower meal has been simulated'", r"no\s+rapeseed\s+or\s+sunflower\s+meal", _MSS),
    ("M4d: the rapeseed keyword", r"oilseed\s+meal;\s+rapeseed;", _MS),
    ("M8: 'more than twenty times too slow'", r"more\s+than\s+twenty\s+times", _MSS),
    ("M8: 'more than ten times too slow'", r"more\s+than\s+ten\s+times\s+too\s+slow", _MSS),
    ("M10: 'failed-leaf'", r"(?i)failed-leaf", _MSS),
    ("M10: 'formulation of 2026-09-20'", r"formulation\s+of\s+2026-09-20", _MSS),
    ("M10: 'can be made available on reasonable request'", r"can\s+be\s+made\s+available", _MS),
    ("M10: 'in the author's repository'", r"in\s+the\s+author's\s+repository", _MS),
    ("M10: 'our radius sweep'", r"our\s+radius\s+sweep", _MS),
    ("M10: 'AMD EPYC' in the main text", r"EPYC", _MS),
    ("M10: 'server' and 'workstation' in the main text", r"(?i)\bserver\b|workstation", _MS),
    ("M10: 'prolongation bootstrap'", r"(?i)prolongation\s+bootstrap", _MSS),
    ("M10: 'pinned to four published determinations'", r"pinned\s+to\s+four", _MS),
    ("M10: 'still to be added'", r"still\s+to\s+be\s+added", _MSS),
    ("S15: 'used for nothing'", r"used\s+for\s+nothing", _MSS),
    ("M10: the queued uniform marches", r"queued\s+behind", _MSS),
    ("M10: dated edit note in S9", r"Moved\s+here\s+from\s+Section~5", _MSS),
    ("M10: dated edit note in S9.10", r"moved\s+there\s+on\s+2026", _MSS),
    ("M10: 'as the manuscript printed it before'", r"as\s+the\s+manuscript\s+printed\s+it\s+before", _MSS),
    ("M10: acquisition dates in S5.1", r"acquired\s+on\s+(11|22)\s+September", _MSS),
    ("M10: internal item and file names in the printed supplement", r"item\\_39|item\\_28|LIMITS\\_AUDIT|w95p|\bw96\b", _MSS),
    ("M10: 'the banked 2026-09-03 value'", r"banked", _MSS),
    ("M10: software narrative in S14", r"rollback\s+indicator|job's\s+copy\s+of\s+the\s+code|interface\s+that\s+would\s+execute", _MSS),
    ("claim 2: 'kinetics set discharge residual'", r"kinetics\s+set\s+discharge", _MS),
    ("claim 4: 'verified reference for equipment-scale'", r"verified\s+reference\s+for\s+equipment-scale", _MS),
]
ASSERTIONS += [
    dict(id=f"P6Z{i + 1}", expect="absent", printed_only=True, sources=src,
         where=f"withdrawn wording absent (part 6): {label}", grep=pat)
    for i, (label, pat, src) in enumerate(_P6ZL)
]

# ---------------------------------------------------------------------------
# 2026-10-01, REFEREE_AUDIT_8 part 7 (section 8.R): Table 3's sunflower cell of
# the "each 2019 meal's diameter" row (the sphere of the 2019 sunflower meal's
# own diameter, radius 0.670 mm) re-scored on version 2 of the thesis
# digitization by the paper's own scorer over the rule window
# (GT_PS2_PAPER1_TABLE3_SUNFLOWER_RESCORE_2026-10-01), staged through
# DTDC_PAPER_TABLE3_RESCORE until it lands at _P7.  C14 is re-keyed in place;
# the version-1 value stays printed as a pointer with its item_37 binding.
# ---------------------------------------------------------------------------
_P7 = "../docs/evidence/paper1_table3_sunflower_rescore_2026-10-01"
PENDING_ROOTS[_P7] = "DTDC_PAPER_TABLE3_RESCORE"
EXTRA_CORPUS_ROOTS += [_P7]
_P7J = _P7 + "/rescore_output.json"
_P7RMS = r'"rms_over_span_as_printed": 0\.11584103586403083'
for _a in ASSERTIONS:
    if _a["id"] == "C14":
        _a.update(printed=11.6, record=_P7J, grep=_P7RMS)
        _a.pop("page", None)
        _a["where"] = ("Table 3 sunflower at 0.670 mm, version 2, rule window"
                       " [re-keyed to the version-2 re-score, 2026-10-01]")
_P7L = [
    (0.116, 'Table S21 the 0.670 mm sunflower row on version 2, fraction of the span', _P7J, _P7RMS),
    (13, 'Table 3 / S21 the 0.670 mm sunflower row on version 2: points in the reading band', _P7J,
     r'"in_band_as_printed": 13,'),
    (14, 'Table 3 the 0.670 mm sunflower row on version 2: points in the rule window', _P7J,
     r'"window_points": 14\b'),
    (9.1, 'Table 3 footnote a / extended report: the version-1 score of the 0.670 mm row, kept as a pointer', _I37,
     r"0\.670 mm for sunflower \(9\.1 per cent, 14 of 15\)"),
    (0.091, 'S21 note: the version-1 score of the 0.670 mm row, fraction of the span', _P7J,
     r'"rms_normalized_by_window_span": 0\.09106158238364007'),
    (0.670, 'Table 3 footnote a / S21 note: the radius of the 2019 sunflower meal\'s sphere, mm', _I37,
     r"0\.670 mm for sunflower", '0.670'),
]


def _p7(i, printed, where, record, grep, page=None):
    d = dict(id=f"P7_{i}", printed=printed, unit="-", where=where, record=record, grep=grep)
    if page is not None:
        d["page"] = page
    if isinstance(printed, int) and printed < 100:
        d["word"] = True  # a small integer: the page check cannot tell it from other text
    return d


ASSERTIONS += [_p7(i + 1, *row) for i, row in enumerate(_P7L)]

_ALL3 = _MSS + ["supplementary_extended.tex"]
_P7ZL = [
    ("Table 3: the 'version-1 points' cell", r"version-1\s+points", _ALL3),
    ("Table 3 footnote a: 'Scored against version~1 ... and not re-scored'",
     r"Scored\s+against\s+version~1\s+of\s+the\s+thesis\s+digitization\s+and\s+not\s+re-scored", _ALL3),
    ("Table S21: 'not re-scored; of fifteen'", r"not\s+re-scored;\s+of\s+fifteen", _ALL3),
    ("extended report: the 0.6697 mm march 'not re-scored on version 2'",
     r"digitization\s+and\s+not\s+re-scored\s+on\s+version\s+2", _ALL3),
    ("abstract: the clock without 'its own crossing'", r"passes\s+from\s+the\s+critical\s+loading", _MS),
    ("main text: a literal supplement pointer (S<n>) in place of a label",
     r"(?<![A-Za-z0-9\\{])S[0-9]+(?:\.[0-9]+)*", _MS),
]
ASSERTIONS += [
    dict(id=f"P7Z{i + 1}", expect="absent", printed_only=True, sources=src,
         where=f"withdrawn wording absent (part 7): {label}", grep=pat)
    for i, (label, pat, src) in enumerate(_P7ZL)
]


# Current reduced continuation; older step-4 assertions still bind the
# preserved historical reports, not the revised reader-facing trajectory.
_TC = "analysis/results_2026-10-01/faner2019_thermal_continuation/continuation.json"
EXTRA_CORPUS_ROOTS += ["analysis/results_2026-10-01/faner2019_thermal_continuation"]
_TCL = [
     (20.87, "continued soybean crossing interval, s", r'"crossing_matched_interval_s": 20\.8720977373767'),
     (0.82, "continued soybean crossing-matched ratio", r'"crossing_matched_ratio": 0\.8233425586530507'),
     (0.64, "continued soybean sample-clock ratio", r'"sample_clock_ratio_after_measured_crossing": 0\.6358297844240106'),
     (41.1, "continued soybean crossing, s", r'"model_crossing_0p20_s": 41\.11024857528673'),
     (41.11, "continued soybean crossing, s", r'"model_crossing_0p20_s": 41\.11024857528673'),
     (62.0, "continued soybean target on sample clock, s", r'"model_0p10_s": 61\.98234631266343'),
     (61.98, "continued soybean target on sample clock, s", r'"model_0p10_s": 61\.98234631266343'),
     (12.28, "continued soybean temperature RMS, K", r'"rms": 12\.277843960666063'),
     (12.3, "continued soybean temperature RMS, K", r'"rms": 12\.277843960666063'),
     (16.01, "historical held soybean temperature RMS, K", r'"rms": 16\.01170311276723'),
     (12.6, "continued soybean normalized loading RMS, percent", r'"window_loading_rms_normalized": 0\.12606535803776966'),
     (0, "continued soybean loading points in reading band", r'"window_in_reading_band": 0'),
     (0.0128, "step-halving temperature change, K", r'"nominal_minus_fine_max_temperature_k": 0\.012781230626359275'),
     (43.14, "continued model plateau end, s", r'"model_plateau_end_s": 43\.14482688395394'),
     (43.1, "continued model plateau end, s", r'"model_plateau_end_s": 43\.14482688395394'),
     (87.09, "measured plateau end, s", r'"measured_plateau_end_s": 87\.08804908246739'),
]
ASSERTIONS += [dict(id=f"TC{i+1}", printed=v, unit="-", where=w,
                       record=_TC, grep=g, word=(v == 0))
                  for i, (v, w, g) in enumerate(_TCL)]


# Direct binary histories remain separate from all scalar uptake records.
_BK = "analysis/results_2026-10-01/binary_particle_kinetics/summary.json"
EXTRA_CORPUS_ROOTS += ["analysis/results_2026-10-01/binary_particle_kinetics"]
ASSERTIONS += [
    dict(id="BK1", printed=6.16e-6, unit="kg/kg", where="binary half-step water difference",
         record=_BK, grep=r'"water_kg_kg_dry": 6\.1558925439908485e-06'),
    dict(id="BK2", printed=1.28e-3, unit="kg/kg", where="binary water-richer mesh difference",
         record=_BK, grep=r'"water_kg_kg_dry": 0\.001277034879409264'),
    dict(id="BK3", printed=700, unit="s", where="binary fine-grid last accepted time",
         record=_BK, grep=r'"last_accepted_time_s": 700\.0'),
    dict(id="BK4", printed=0.12118387, unit="kg/kg", where="binary initial water inventory",
         record=_BK, grep=r'"water_kg_kg_dry": 0\.12118386956873105'),
    dict(id="BK5", printed=0.01064392, unit="kg/kg", where="binary initial hexane inventory",
         record=_BK, grep=r'"hexane_kg_kg_dry": 0\.01064391940026349'),
]


# Separate nitrogen/one-sorbate identifications (2026-10-01).  Preserve all
# historical scalar-uptake assertions.  Each Table 8 row binds its activity,
# pore mobility, fitted count and fractional RMS in the same named JSON row;
# the manuscript's RMS column is 100 * fit_rms, not an apparent diffusivity.
_CK = "analysis/results_2026-10-01/cardarelli_carrier_kinetics/fits.json"
EXTRA_CORPUS_ROOTS += [_CK]
_CK_ROWS = [
    ("water", "0.218", "7.808065728139263e-08", 7, "0.05036091562476088",
     (("CK1", 0.218, "activity"), ("CK2", 7.81e-8, "D"),
      ("CK3", 7, "points"), ("CK4", 5.04, "RMS percent"))),
    ("water", "0.422", "5.7012535857401553e-08", 8, "0.002591486613072366",
     (("CK5", 0.422, "activity"), ("CK6", 5.70e-8, "D"),
      ("CK7", 8, "points"), ("CK8", 0.259, "RMS percent"))),
    ("water", "0.624", "5.174079723488737e-08", 6, "0.044294648833983215",
     (("CK9", 0.624, "activity"), ("CK10", 5.17e-8, "D"),
      ("CK11", 6, "points"), ("CK12", 4.43, "RMS percent"))),
    ("water", "0.806", "2.7288571775353137e-08", 9, "0.06923057017863649",
     (("CK13", 0.806, "activity"), ("CK14", 2.73e-8, "D"),
      ("CK15", 9, "points"), ("CK16", 6.92, "RMS percent"))),
    ("hexane", "0.212", "1.5537357336562158e-10", 9, "0.07124694775305075",
     (("CK17", 0.212, "activity"), ("CK18", 1.55e-10, "D"),
      ("CK19", 9, "points"), ("CK20", 7.12, "RMS percent"))),
    ("hexane", "0.423", "2.463825795080323e-10", 9, "0.05902710309466754",
     (("CK21", 0.423, "activity"), ("CK22", 2.46e-10, "D"),
      ("CK23", 9, "points"), ("CK24", 5.90, "RMS percent"))),
    ("hexane", "0.621", "3.6726765925702894e-10", 9, "0.0862215356814302",
     (("CK25", 0.621, "activity"), ("CK26", 3.67e-10, "D"),
      ("CK27", 9, "points"), ("CK28", 8.62, "RMS percent"))),
    ("hexane", "0.813", "5.893139688739871e-10", 9, "0.09104244822708095",
     (("CK29", 0.813, "activity"), ("CK30", 5.89e-10, "D"),
      ("CK31", 9, "points"), ("CK32", 9.10, "RMS percent"))),
]
for _sorbate, _activity, _diffusivity, _count, _rms, _bindings in _CK_ROWS:
    # Stop at the next mobility key: a changed cell cannot match a later
    # series's count or RMS.  The config identifies the sorbate and activity.
    _ck_row = (
        r'"config":\s*\{\s*"sorbate": "' + _sorbate
        + r'",\s*"activity_boundary": ' + re.escape(_activity) + r','
        + r'(?:(?!"diffusivity_m2_s").)*"diffusivity_m2_s": '
        + re.escape(_diffusivity) + r',\s*"fit_method": "[^"\n]+",'
        + r'\s*"fitted_point_count": ' + str(_count)
        + r',\s*"fit_rms": ' + re.escape(_rms) + r','
    )
    for _id, _value, _column in _bindings:
        ASSERTIONS.append(dict(
            id=_id, printed=_value, unit="-",
            where=f"Table 8 {_sorbate} a={_activity}, {_column}",
            record=_CK, grep="(?s)" + _ck_row,
            page={"CK6": "5.70", "CK24": "5.90", "CK32": "9.10"}.get(_id),
        ))

ASSERTIONS += [
    dict(id="CK33", printed=0.00158, unit="-", where="Sec. 5.1 fixed-D mesh sensitivity, minimum (water a=0.422)",
         record=_CK, grep=r'"mesh_max_loading_difference": 0\.0015834161825498405,'),
    dict(id="CK34", printed=0.00741, unit="-", where="Sec. 5.1 fixed-D mesh sensitivity, maximum (hexane a=0.212)",
         record=_CK, grep=r'"mesh_max_loading_difference": 0\.007413630386539791,'),
    dict(id="CK35", printed=24, unit="cells", where="Sec. 5.1 nominal carrier-fit mesh",
         record=_CK, grep=r'"cells": 24,\s*"refined_cells": 48,'),
    dict(id="CK36", printed=48, unit="cells", where="Sec. 5.1 refined carrier-fit mesh",
         record=_CK, grep=r'"cells": 24,\s*"refined_cells": 48,'),
    # The carrier harness uses this existing soybean Luikov prior, not a
    # newly fitted sunflower isotherm.  Bind its parameters to the already
    # landed binary configuration; do not treat them as fit outputs.
    dict(id="CK37", printed=0.880, page="0.880", unit="kg/kg", where="Sec. 5.1 soybean water prior A1",
         record=_BK, grep=r'"luikov":\s*\{\s*"A1": 0\.88,'),
    dict(id="CK38", printed=12.184, unit="-", where="Sec. 5.1 soybean water prior A2",
         record=_BK, grep=r'"luikov":\s*\{\s*"A1": 0\.88,\s*"A2": 12\.184,'),
    dict(id="CK39", printed=0.023, unit="kg/kg", where="Sec. 5.1 soybean water prior lower anchor",
         record=_BK, grep=r'"W_cap": 0\.235670904,\s*"W_ref": 0\.023'),
    dict(id="CK40", printed=0.2356709040, page="0.2356709040", unit="kg/kg", where="Sec. 5.1 soybean water prior evidence cap",
         record=_BK, grep=r'"W_cap": 0\.235670904,\s*"W_ref": 0\.023'),
]


def resolve_record(rec: str):
    """(path or None, pending environment variable or None) for a record.

    A record under a root of ``PENDING_ROOTS`` that has not landed resolves to
    its staging folder when the named variable is set, else to None."""
    import os
    path = PAPER / rec
    if path.exists():
        return path, None
    for root, env in sorted(PENDING_ROOTS.items(), key=lambda kv: -len(kv[0])):
        if rec == root or rec.startswith(root + "/"):
            staged = os.environ.get(env)
            if staged:
                rest = rec[len(root):].lstrip("/")
                return (Path(staged) / rest if rest else Path(staged)), None
            return None, env
    return path, None


#: Pass-1 literals that are not results, each with the reason it is excluded.
#: They are declared here rather than left as bare "unrecorded" pointers.
PASS1_DECLARED = {
     "101.325": "declared parameter: the bed pressure, Table S2",
     "6.16e-6": "binary water half-step maximum, rounded; BK1 binds the exact JSON value",
    "223.25": "derived quotient of the two declared diffusivities, printed with"
              " its basis in Sec. S4.3",
    "3.4228875": "mantissa of 3.4228875e-10; the tool cannot fold a LaTeX"
                 " \\times10^{n} into its mantissa (assertion V3 finds it)",
    "254000": "Table 1 feed hexane, 254252 ppm rounded to the thousand; DT1 binds the record",
    "216000": "Table 1 after-PD hexane at 0.50 mm, 215915 ppm rounded; DT4 binds the record",
    "225000": "Table 1 after-PD hexane at 0.885 mm, 225410 ppm rounded; DT5 binds the record",
    "232000": "Table 1 after-PD hexane at 1.5 mm, 232405 ppm rounded; DT6 binds the record",
    "3200": "Table 1 after-strip hexane at 0.885 and 1.5 mm, 3215 and 3244 ppm rounded; DT11, DT12",
    "3600": "Table 1 after-strip hexane at 0.50 mm, 3649 ppm rounded; DT10 binds the record",
}


# ---------------------------------------------------------------------------
# 2026-10-02/03: the desolventizer march (Sec. 5.1, Table 1, Fig. 2) and the
# face-heated layer march (Sec. 5.4, Table 2, Fig. 3).  Records:
# v5_stations.json (station values read from the six dt_R*_v5 runs by
# paper/analysis/harness/dtdc_v5_stations.py) and layer_march.json.
_V5 = "analysis/results_2026-10-01/dtdc_particle_destiny/v5_stations.json"
_LM = "analysis/results_2026-10-01/faner2019_layer_march/layer_march.json"
_ENS = "analysis/results_2026-10-01/dtdc_particle_destiny/fick_ensemble_exit_residual.json"
EXTRA_CORPUS_ROOTS += [_V5, _LM, _ENS,
                       "analysis/results_2026-10-01/faner2019_layer_march"]
_DTL = [
    (254252, "feed hexane, ppm wet", r'"hexane_ppm_wb": 254252\.5194988399'),
    (49, "feed temperature, C", r'"T_c": 49\.02011022717181'),
    (8.1, "feed moisture, percent wet", r'"moisture_wb_pct": 8\.108735078765733'),
    (215915, "after PD hexane 0.50 mm, ppm", r'"hexane_ppm_wb": 215914\.5'),
    (225410, "after PD hexane 0.885 mm, ppm", r'"hexane_ppm_wb": 225409\.8'),
    (232405, "after PD hexane 1.5 mm, ppm", r'"hexane_ppm_wb": 232404\.8'),
    (0.78, "strip duration 0.50 mm, s", r'"strip_duration_s": 0\.7777580965951643'),
    (1.8, "strip duration 0.885 mm, s", r'"strip_duration_s": 1\.7977580965942366'),
    (4.2, "strip duration 1.5 mm, s", r'"strip_duration_s": 4\.217758096592036'),
    (3649, "after strip hexane 0.50 mm, ppm", r'"hexane_ppm_wb": 3649\.355234888234'),
    (3215, "after strip hexane 0.885 mm, ppm", r'"hexane_ppm_wb": 321[45]\.'),
    (3244, "after strip hexane 1.5 mm, ppm", r'"hexane_ppm_wb": 324[34]\.'),
    (18.9, "exit moisture at 105 C, percent", r'"moisture_wb_pct": 18\.86'),
    (14.1, "exit moisture at 110 C, percent", r'"moisture_wb_pct": 14\.1[23]'),
    (17, "exit hexane 0.50 mm at 105 C, ppm", r'"hexane_ppm_wb": 1[67]\.'),
    (31, "exit hexane 0.885 mm at 105 C, ppm", r'"hexane_ppm_wb": 3[01]\.'),
    (259, "exit hexane 1.5 mm at 105 C, ppm", r'"hexane_ppm_wb": 259\.'),
    (14, "exit hexane 0.50 mm at 110 C, ppm", r'"hexane_ppm_wb": 14\.'),
    (23, "exit hexane 0.885 mm at 110 C, ppm", r'"hexane_ppm_wb": 23\.'),
    (197, "exit hexane 1.5 mm at 110 C, ppm", r'"hexane_ppm_wb": 197\.'),
    (18.6, "PD removal 0.50 mm, percent", r'"pd_hexane_removed_fraction": 0\.185'),
    (14.1, "PD removal 0.885 mm, percent", r'"pd_hexane_removed_fraction": 0\.1407'),
    (10.7, "PD removal 1.5 mm, percent", r'"pd_hexane_removed_fraction": 0\.107'),
    (74, "implied jacket 0.50 mm at 110 C, W/kg", r'"implied_jacket_w_per_kg_dry": 73\.9'),
    (81, "implied jacket 1.5 mm at 110 C, W/kg", r'"implied_jacket_w_per_kg_dry": 80\.8'),
    (-49, "implied jacket 1.5 mm at 105 C, W/kg", r'"implied_jacket_w_per_kg_dry": -48\.8'),
    (-55, "implied jacket 0.885 mm at 105 C, W/kg", r'"implied_jacket_w_per_kg_dry": -54\.5'),
    (0.52, "Ackermann maximum in the strip (B, 1.5 mm)", r'"B": 0\.51[6-7]'),
    (1.09, "Ackermann maximum at the toaster hand-off", r'"toaster": 1\.08[5-6]'),
]
ASSERTIONS += [dict(id=f"DT{i+1}", printed=v, unit="-", where=w, record=_V5, grep=g)
               for i, (v, w, g) in enumerate(_DTL)]
_ENSL = [
    (36, "ensemble residual at 105 C, D 2.6e-10, ppm", r'"mass_weighted_ppm": 3[56]\.'),
    (306, "ensemble residual at 105 C, lowest D, ppm", r'"mass_weighted_ppm": 30[5-6]\.'),
    (18, "ensemble residual at 105 C, highest D, ppm", r'"mass_weighted_ppm": 1[78]\.'),
    (23, "ensemble residual at 110 C, D 2.6e-10, ppm", r'"mass_weighted_ppm": 2[23]\.'),
    (198, "ensemble residual at 110 C, lowest D, ppm", r'"mass_weighted_ppm": 19[78]\.'),
    (15, "ensemble residual at 110 C, highest D, ppm", r'"mass_weighted_ppm": 1[45]\.'),
]
ASSERTIONS += [dict(id=f"EN{i+1}", printed=v, unit="ppm", where=w, record=_ENS, grep=g)
               for i, (v, w, g) in enumerate(_ENSL)]
_LML = [
    (5.4, "layer depth, mm", r'"layer_depth_m": 0\.00540'),
    (92, "face film bare h, W/m2K", r'"bare_h_w_m2_k": 92\.37'),
    (81, "face film effective h, W/m2K", r'"effective_h_w_m2_k": 80\.58'),
    (0.057, "bed conductivity, W/mK", r'"effective_w_m_k": 0\.0567'),
    (0.040, "Zehner-Schluender part, W/mK", r'"zehner_schluender_w_m_k": 0\.0404'),
    (0.016, "radiative part, W/mK", r'"radiative_w_m_k_at_film_mean_T": 0\.016'),
    (43.3, "wicking ends, s", r'"stage1_end_s": 43\.3'),
    (78.4, "layer plateau end, s", r'"model_plateau_end_s": 78\.38'),
    (87.1, "measured plateau end, s", r'"measured_plateau_end_s": 87\.088'),
    (4.1, "layer temperature RMS, K", r'"rms": 4\.085'),
    (37.1, "layer crossing-to-0.10 interval, s", r'"crossing_matched_interval_s": 37\.1'),
    (1.46, "layer crossing-matched ratio", r'"crossing_matched_ratio": 1\.46'),
    (9.5, "layer window loading RMS, percent", r'"window_loading_rms_normalized": 0\.095'),
    (25.4, "measured crossing-to-0.10 interval, s", r'"measured_0p10_s": 71\.2142'),
    (71.0, "plateau end at 2k, s", r'"model_plateau_end_s": 71\.0'),
    (2.6, "temperature RMS at 2k, K", r'"rms": 2\.597'),
    (1.13, "crossing-matched ratio at 2k", r'"crossing_matched_ratio": 1\.128'),
    (6.4, "window RMS at 2k, percent", r'"window_loading_rms_normalized": 0\.0639'),
    (91.2, "plateau end at k/2, s", r'"model_plateau_end_s": 91\.22'),
    (108, "plateau end, planar face, s", r'"model_plateau_end_s": 107\.68'),
    (109, "0.10 reached, planar face, s", r'"model_0p10_s": 108\.96'),
    (36, "0.10 reached, two-sided, s", r'"model_0p10_s": 35\.68'),
    (61, "plateau end, two-sided, s", r'"model_plateau_end_s": 61\.1'),
    (53, "plateau end, probe at quarter depth, s", r'"model_plateau_end_s": 52\.6'),
    (139, "plateau end, probe at three-quarter depth, s", r'"model_plateau_end_s": 138\.5'),
    (96, "plateau end, exposed half diameter, s", r'"model_plateau_end_s": 95\.5'),
    (59, "plateau end, exposed 1.5 diameters, s", r'"model_plateau_end_s": 59\.39'),
    (77, "plateau end, 12 sublayers, s", r'"model_plateau_end_s": 77\.2'),
    (81, "plateau end, 48 sublayers, s", r'"model_plateau_end_s": 80\.8'),
]
ASSERTIONS += [dict(id=f"LM{i+1}", printed=v, unit="-", where=w, record=_LM, grep=g)
               for i, (v, w, g) in enumerate(_LML)]
# the main-text wording withdrawn with the ensemble as the Sec. 5.4 result
_LMZ = [
    ("the 0.82 crossing-matched ratio as the Sec. 5.4 result", r"\$0\.82\$ of the measured time"),
    ("the 0.64 sample-clock ratio", r"\$0\.64\$ of the observed interval"),
    ("the ensemble's 0.71 still-boiling mass", r"\$0\.71\$ of the sample mass"),
    ("ladder as the word for the bed temperature history", r"(?i)\bbed ladders?\b|ladder to \\SI\{1[01][05]\}"),
]
ASSERTIONS += [dict(id=f"LMZ{i+1}", expect="absent", printed_only=True,
                    where=f"withdrawn wording absent: {label}", grep=pat)
               for i, (label, pat) in enumerate(_LMZ)]


def live_sources_text(names: list[str] | None = None) -> str:
    """The sources searched.  Default (2026-09-28): the main text and both
    supplements, so that every absence check runs on the submitted and on the
    extended supplement alike."""
    parts = []
    for name in (MAIN_SOURCES + SUPP_SOURCES + EXT_SOURCES if names is None else names):
        path = HERE / name
        if path.exists():
            parts.append(path.read_text(encoding="utf-8-sig"))
    return "\n".join(parts)


def printed_only_text(text: str) -> str:
    '''The printed prose of LaTeX sources (2026-09-30, part 4): whole-line and
    inline comments removed, and the arguments of the url, label, ref, ref*,
    eqref and autoref commands dropped, since a repository path or a label key
    is not prose.  Line count is preserved so that a hit keeps its line
    number.'''
    out = []
    for raw in text.split("\n"):
        if raw.lstrip().startswith("%"):
            out.append("")
            continue
        line = re.split(r"(?<![\\])%", raw, maxsplit=1)[0]
        line = re.sub(r"\\(?:url|label|ref\*?|eqref|autoref)\{[^}]*\}", " ", line)
        out.append(line)
    return "\n".join(out)


def run_assertions() -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    live = live_sources_text()
    printed = printed_only_text(live)
    for a in ASSERTIONS:
        if a.get("expect") == "absent":
            if a.get("sources"):
                src = live_sources_text(a["sources"])
                m = re.search(a["grep"], printed_only_text(src) if a.get("printed_only") else src)
            else:
                src = live
                m = re.search(a["grep"], printed if a.get("printed_only") else live)
            if m:
                ln = src[: m.start()].count("\n") + 1
                rows.append((a["id"], "STALE-STILL-PRINTED", f"live source line {ln}"))
            else:
                rows.append((a["id"], "ABSENT-OK", a["where"]))
            continue
        root, pending = resolve_record(a["record"])
        if root is None:
            if os.environ.get("DTDC_PUBLIC_PACKAGE"):
                rows.append((a["id"], "RECORD-NOT-SHIPPED", a["record"]))
            else:
                rows.append((a["id"], "AWAITING-LANDED-EVIDENCE", f"{a['record']} (set {pending})"))
            continue
        # the public reproducibility package ships the records the submitted
        # documents rest on, not the internal process records; there, a missing
        # record is reported as such and is not a mismatch
        if not root.exists() and os.environ.get("DTDC_PUBLIC_PACKAGE"):
            rows.append((a["id"], "RECORD-NOT-SHIPPED", a["record"]))
            continue
        pat = re.compile(a["grep"])
        hits: list[str] = []
        if root.exists():
            raw = bool(a.get("raw"))
            targets = [root] if root.is_file() else [
                p for p in root.rglob("*")
                if p.is_file()
                and p.suffix.lower() in CORPUS_SUFFIXES
                and (raw or not (CORPUS_EXCLUDE_DIRS & set(q.name for q in p.parents)))
                and (raw or p.stat().st_size <= CORPUS_MAX_BYTES)
            ]
            for p in targets:
                try:
                    text = p.read_text(encoding="utf-8", errors="replace")
                except OSError:
                    continue
                m = pat.search(text)
                if m:
                    ln = text[: m.start()].count("\n") + 1
                    try:
                        shown = p.relative_to(PAPER)
                    except ValueError:
                        shown = p
                    hits.append(f"{shown}:{ln}")
                    if len(hits) >= 3:
                        break
        verdict = "FOUND" if hits else "NOT-FOUND-IN-NAMED-RECORD"
        rows.append((a["id"], verdict, "; ".join(hits) if hits else a["record"]))
    return rows


# ---------------------------------------------------------------------------
# 2b.  Coverage drift (added 2026-09-24)
# ---------------------------------------------------------------------------
# An assertion binds a printed number to a record.  After a repair pass an
# assertion can keep passing while the manuscript no longer prints the number
# it was written for, or prints it at a coarser precision: the assertion then
# guards the record and not the page.  This pass searches each record-expecting
# assertion's own pattern in the live manuscript sources and says which of them
# still reach the page.  Patterns that are record-structural (JSON keys, prose)
# are reported as not applicable rather than as failures.

def page_tokens(x: float) -> set[str]:
    """The decimal strings a manuscript might print for ``x``.

    A record keeps 1216.358141953133 and the page prints 1216.4; a record keeps
    9.03e-4 and the page prints ``$9.03\\times10^{-4}$``.  So the value at its
    own precision, the same value rounded down to three significant figures,
    and its mantissa all count as the page's form of it.  Nothing shorter than
    three digits is generated unless the value itself is shorter, because a
    one- or two-digit token matches the page by accident.
    """
    if x == 0:
        return {"0"}
    ax = abs(float(x))
    nat = len([c for c in ("%r" % ax).split("e")[0] if c.isdigit()])
    nat = max(1, min(nat, 10))
    out: set[str] = set()
    for sig in range(nat, 2, -1):
        s = f"%.{sig}g" % ax
        if "e" not in s and "E" not in s:
            out.add(s)
            if "." in s:
                out.add(s.rstrip("0").rstrip("."))
    mant = f"%.{nat}g" % (ax / 10.0 ** math.floor(math.log10(ax)))
    if "e" not in mant and len([c for c in mant if c.isdigit()]) >= 3:
        out.add(mant)
        out.add(mant.rstrip("0").rstrip(".") if "." in mant else mant)
    out = {t for t in out if len([c for c in t if c.isdigit()]) >= 3}
    if not out:
        s = f"%.{nat}g" % ax
        out = {s} if "e" not in s else {f"%.{nat}g" % (ax / 10.0 ** math.floor(math.log10(ax)))}
    return out


def run_coverage() -> list[tuple[str, str, str]]:
    live = live_sources_text(MAIN_SOURCES + SUPP_SOURCES)
    ext = live_sources_text(EXT_SOURCES)
    rows: list[tuple[str, str, str]] = []
    for a in ASSERTIONS:
        if a.get("expect") == "absent":
            continue
        if a.get("word"):
            rows.append((a["id"], "N/A-SPELLED-AS-A-WORD", a["where"]))
            continue
        if a.get("withdrawn_from_page"):
            rows.append((a["id"], "WITHDRAWN-FROM-PAGE", a["where"]))
            continue
        # ``page`` (added 2026-09-26) names the literal the page prints when
        # a float cannot carry it, e.g. a trailing zero such as 0.770.
        toks = set(page_tokens(a["printed"])) | ({a["page"]} if a.get("page") else set())
        def _hit(text: str) -> bool:
            return any(
                re.search(r"(?<![0-9.])" + re.escape(tok) + r"(?![0-9.])", text)
                for tok in toks
            )
        if _hit(live):
            verdict = "PRINTED"
        elif _hit(ext):
            verdict = "PRINTED-IN-INTERNAL-ARCHIVE"
        else:
            verdict = "NOT-PRINTED"
        rows.append((a["id"], verdict, a["where"]))
    return rows


# ---------------------------------------------------------------------------
# 3.  Internal consistency: the same number printed twice
# ---------------------------------------------------------------------------

def internal_pairs(main: list[Printed], supp: list[Printed]) -> list[str]:
    """Numbers printed in the main text that appear nowhere in the supplement."""
    supp_tokens = {p.text for p in supp}
    supp_vals = [p.value for p in supp]
    out: list[str] = []
    for p in main:
        if p.text in supp_tokens:
            continue
        if any(matches_at_precision(p.value, p.text, v) for v in supp_vals):
            continue
        out.append(f"{p.file}:{p.line}  {p.text}   | {p.context[:110]}")
    return out


# ---------------------------------------------------------------------------

def main() -> int:
    main_printed = read_printed(MAIN_SOURCES)
    supp_printed = read_printed(SUPP_SOURCES)
    print(f"printed numeric literals: main {len(main_printed)}, supplement {len(supp_printed)}")

    corpus = build_corpus()
    print(f"record corpus: {len(corpus.files)} curated files, {len(corpus.values)} numeric tokens")

    unmatched_main = []
    for p in main_printed:
        if not lookup(corpus, p, limit=1):
            unmatched_main.append(p)
    undeclared = [p for p in unmatched_main if p.text not in PASS1_DECLARED]
    print(f"\nPASS 1  main-text literals with no corpus hit at printed precision: "
          f"{len(unmatched_main)} of {len(main_printed)}")
    for p in unmatched_main:
        note = PASS1_DECLARED.get(p.text)
        tag = f"DECLARED ({note})" if note else "UNRECORDED"
        print(f"  {p.file}:{p.line}  {p.text}  {tag}")

    print("\nPASS 2  curated assertions")
    rows = run_assertions()
    for aid, verdict, where in rows:
        print(f"  {aid:4s} {verdict:26s} {where}")

    print("\nPASS 2b  assertion coverage: submitted versus preserved internal sources")
    cov = run_coverage()
    for aid, verdict, where in cov:
        if verdict not in ("PRINTED", "PRINTED-IN-INTERNAL-ARCHIVE"):
            print(f"  {aid:4s} {verdict:30s} {where}")
    print(f"  ({sum(1 for _, v, _ in cov if v == 'PRINTED')} printed, "
          f"{sum(1 for _, v, _ in cov if v == 'PRINTED-IN-INTERNAL-ARCHIVE')} internal archive only, "
          f"{sum(1 for _, v, _ in cov if v == 'NOT-PRINTED')} in neither source, "
          f"{sum(1 for _, v, _ in cov if v == 'N/A-SPELLED-AS-A-WORD')} spelled as a word)")

    print("\nPASS 3  main-text numbers absent from the supplement")
    for row in internal_pairs(main_printed, supp_printed):
        print(f"  {row}")

    mismatch = sum(1 for _, v, _ in rows if v == "NOT-FOUND-IN-NAMED-RECORD")
    stale = sum(1 for _, v, _ in rows if v == "STALE-STILL-PRINTED")
    unrecorded = len(undeclared)
    print("\nCOUNTS")
    print(f"  assertions checked against a record : "
          f"{sum(1 for a in ASSERTIONS if a.get('expect') != 'absent')}")
    print(f"  superseded literals checked absent  : "
          f"{sum(1 for a in ASSERTIONS if a.get('expect') == 'absent')}")
    print(f"  mismatch                            : {mismatch}")
    print(f"  awaiting landed evidence            : "
          f"{sum(1 for _, v, _ in rows if v == 'AWAITING-LANDED-EVIDENCE')}")
    if os.environ.get("DTDC_PUBLIC_PACKAGE"):
        print(f"  record not shipped in the package   : "
              f"{sum(1 for _, v, _ in rows if v == 'RECORD-NOT-SHIPPED')}")
    print(f"  stale                               : {stale}")
    print(f"  unrecorded                          : {unrecorded}")
    print(f"  assertions printed in submitted sources: "
          f"{sum(1 for _, v, _ in cov if v == 'PRINTED')}")
    print(f"  assertions only in preserved internal archive: "
          f"{sum(1 for _, v, _ in cov if v == 'PRINTED-IN-INTERNAL-ARCHIVE')}")
    print(f"  assertions no longer in any preserved source: "
          f"{sum(1 for _, v, _ in cov if v == 'NOT-PRINTED')}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
