"""Shared drawing style for the current graphical abstract."""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle

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
