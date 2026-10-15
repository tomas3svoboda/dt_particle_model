"""Current graphical abstract from completed native particle experiments."""

import json
import textwrap
from pathlib import Path

import fig_graphical_abstract as style

PAPER = Path(__file__).resolve().parents[1]


def main():
    audit = json.loads((PAPER / "core2_particle_numbers_audit.json").read_text(encoding="utf-8"))
    assert audit["passed"]
    numbers = audit["macros"]
    sentences = [
        "One conservative particle model couples hexane removal, water uptake and heat transfer through three phase regimes.",
        f"The prescribed DT particle exits with {numbers['DTFiniteMoisture']}% moisture and {numbers['DTFiniteHexaneWet']} ppm hexane on a wet basis.",
        f"The Faner falling interval is 25.4 s measured and {numbers['FanerFallingInterval']} s predicted.",
    ]
    prose = " ".join(sentences)
    assert len(prose.split()) <= 80
    fig = style.plt.figure(figsize=(400 / 72, 400 / 72), dpi=72)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 100), ylim=(0, 100))
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.text(50, 96, "One particle model, two prescribed baths", ha="center", va="top", fontsize=12)
    for x, kind, label in zip((18, 50, 82), "ABC", ("wet", "moving front", "radial dry tail")):
        style.draw_particle(ax, x, 79, 9.2, kind)
        ax.text(x, 66, label, ha="center", va="top", fontsize=10)
    for x in (30, 62):
        ax.annotate(
            "", xy=(x + 9, 79), xytext=(x, 79), arrowprops={"arrowstyle": "->", "color": style.INK}
        )
    ax.text(
        50,
        55,
        "Faner: prescribed hexane bath\nDT: prescribed steam-rich bath",
        ha="center",
        va="top",
        fontsize=11,
        linespacing=1.4,
    )
    for y, sentence in zip((43, 28, 13), sentences):
        ax.text(
            8, y, textwrap.fill(sentence, width=58),
            ha="left", va="top", fontsize=10, linespacing=1.25,
        )
    out = PAPER / "collaterals"
    out.mkdir(exist_ok=True)
    for suffix in ("pdf", "png", "tif"):
        fig.savefig(out / f"graphical_abstract_core2.{suffix}", dpi=72)
    style.plt.close(fig)
    (out / "graphical_abstract_core2_text.txt").write_text(prose + "\n", encoding="utf-8")
    print(f"{len(prose.split())} words, three sentences")


if __name__ == "__main__":
    main()
