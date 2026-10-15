"""JFPE graphical abstract for the particle manuscript.

Canvas is 400 x 400 px at 72 dpi, the square size allowed by the journal.
Type is Arial at 10–12 pt. The prose is three sentences. Preferred upload
files are the PDF and the TIFF written beside this script's output.

Numbers are the printed v6 exits: strip 0.8–4.2 s; moisture 18.9% at
105 °C and 14.1% at 110 °C; hexane 16/28/369 ppm and 14/25/387 ppm.
"""

from __future__ import annotations

import os
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
OUT = os.path.join(PAPER, "collaterals")

for name in ("arial.ttf", "arialbd.ttf"):
    path = os.path.join(os.environ.get("WINDIR", "C:/Windows"), "Fonts", name)
    if os.path.isfile(path):
        font_manager.fontManager.addfont(path)
plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEXANE = "#1B4F72"
MEAL = "#F3EEE6"
EDGE = "#1C1C1C"
WATER = "#1A7A72"
STEAM = "#A15C38"
INK = "#1A1A1A"
RULE = "#C8C2B8"

SENTENCES = [
    "With a prescribed front temperature, condensation heat at the surface removes free n-hexane in 0.8–4.2 s; moisture reaches 18.9% at 105 °C (14.1% at 110 °C).",
    "After the front vanishes, the 1.5 mm sphere retains 369 and 387 ppm, inside the published 100–500 ppm class for that fraction, not a plant average.",
    "The 0.50 mm sphere reaches the exit isotherm at 16 and 14 ppm.",
]


def _words(text: str) -> int:
    return len(text.split())


def draw_particle(ax, cx, cy, radius, kind):
    if kind == "A":
        ax.add_patch(Circle((cx, cy), radius, facecolor=HEXANE, edgecolor=EDGE, lw=1.1, zorder=2))
        for dx, dy in ((-0.28, 0.22), (0.30, 0.05), (-0.05, -0.32)):
            ax.add_patch(
                Circle((cx + dx * radius, cy + dy * radius), 0.16 * radius,
                       facecolor=WATER, edgecolor="none", zorder=3)
            )
    elif kind == "B":
        ax.add_patch(Circle((cx, cy), radius, facecolor=MEAL, edgecolor=EDGE, lw=1.1, zorder=2))
        ax.add_patch(Circle((cx, cy), 0.52 * radius, facecolor=HEXANE, edgecolor=EDGE, lw=0.9, zorder=3))
        for dx, dy in ((0.72, 0.28), (-0.70, 0.35), (0.62, -0.48)):
            ax.add_patch(
                Circle((cx + dx * radius, cy + dy * radius), 0.11 * radius,
                       facecolor=WATER, edgecolor="none", zorder=3)
            )
    elif kind == "C":
        ax.add_patch(Circle((cx, cy), radius, facecolor=MEAL, edgecolor=EDGE, lw=1.1, zorder=2))
        for dx, dy, color in (
            (-0.35, 0.40, HEXANE),
            (0.42, 0.22, HEXANE),
            (0.05, -0.38, HEXANE),
            (-0.48, -0.18, WATER),
            (0.38, -0.42, WATER),
            (0.10, 0.48, WATER),
        ):
            ax.add_patch(
                Circle((cx + dx * radius, cy + dy * radius), 0.12 * radius,
                       facecolor=color, edgecolor="none", zorder=3)
            )
    else:
        raise ValueError(kind)


def main() -> None:
    prose = " ".join(SENTENCES)
    n_words = _words(prose)
    if len(SENTENCES) != 3 or n_words > 80:
        raise SystemExit(f"prose is {len(SENTENCES)} sentences and {n_words} words")

    os.makedirs(OUT, exist_ok=True)
    side = 400 / 72
    fig = plt.figure(figsize=(side, side), dpi=72)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    centers = (17.5, 50.0, 82.5)
    labels = ("wet", "front", "diffusion")
    kinds = ("A", "B", "C")
    cy, radius = 79.0, 9.2
    for cx, kind, letter, label in zip(centers, kinds, "ABC", labels):
        draw_particle(ax, cx, cy, radius, kind)
        ax.text(cx, cy - radius - 3.2, letter, ha="center", va="top",
                fontsize=12, fontweight="bold", color=INK)
        ax.text(cx, cy - radius - 8.0, label, ha="center", va="top",
                fontsize=10, color=INK)

    for x0, x1, caption in ((29.2, 38.4, "0.8–4.2 s"), (61.6, 70.8, None)):
        ax.add_patch(
            FancyArrowPatch(
                (x0, cy), (x1, cy),
                arrowstyle="-|>", mutation_scale=9, lw=1.0,
                color=STEAM if caption else INK, zorder=1,
            )
        )
        if caption:
            ax.text((x0 + x1) / 2, cy + 3.3, caption, ha="center", va="bottom",
                    fontsize=10, color=STEAM)

    ax.plot([8, 92], [60.5, 60.5], color=RULE, lw=0.6, solid_capstyle="butt")

    ax.text(50, 57.6, "Exit hexane, ppm", ha="center", va="top",
            fontsize=11, fontweight="bold", color=INK)
    ax.text(46, 53.2, "105 °C", ha="center", va="top", fontsize=10,
            fontweight="bold", color=INK)
    ax.text(68, 53.2, "110 °C", ha="center", va="top", fontsize=10,
            fontweight="bold", color=INK)

    rows = (
        ("0.50 mm", "16", "14", "isotherm"),
        ("0.885 mm", "28", "25", ""),
        ("1.50 mm", "369", "387", ""),
    )
    y = 49.4
    for radius_label, a, b, note in rows:
        weight = "bold" if radius_label.startswith("0.50") or radius_label.startswith("1.50") else "regular"
        ax.text(22, y, radius_label, ha="left", va="top", fontsize=10,
                fontweight=weight, color=INK)
        ax.text(46, y, a, ha="center", va="top", fontsize=10,
                fontweight=weight, color=INK)
        ax.text(68, y, b, ha="center", va="top", fontsize=10,
                fontweight=weight, color=INK)
        if note:
            ax.text(78, y, note, ha="left", va="top", fontsize=10, color=WATER)
        y -= 4.6

    wrapped = textwrap.fill(prose, width=64)
    ax.text(8, 33.2, wrapped, ha="left", va="top", fontsize=10, color=INK,
            linespacing=1.25)

    pdf = os.path.join(OUT, "graphical_abstract.pdf")
    tif = os.path.join(OUT, "graphical_abstract.tif")
    png = os.path.join(OUT, "graphical_abstract.png")
    fig.savefig(pdf, dpi=72)
    fig.savefig(png, dpi=72)
    fig.savefig(tif, dpi=72)
    plt.close(fig)

    sidecar = os.path.join(OUT, "graphical_abstract_text.txt")
    with open(sidecar, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(prose + "\n")
        handle.write(f"\n{n_words} words, {len(SENTENCES)} sentences.\n")
    print(f"words={n_words} sentences={len(SENTENCES)}")
    print(pdf)


if __name__ == "__main__":
    main()
