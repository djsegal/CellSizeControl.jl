#!/usr/bin/env python3
"""Shared publication-figure style for the cell-size-control analysis figures.

Ports the science-space PubPlots conventions to matplotlib so every plot_*.py
script draws from one place: the Okabe-Ito colorblind-safe palette and a clean,
journal-ready rcParams (sans-serif, top/right spines off, ticks out, frameless
legends, 300 dpi, tight bounding box).

Usage:
    from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE
    apply_style()
"""
from __future__ import annotations

import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe categorical palette (Okabe & Ito 2008, "Color
# Universal Design"). Named constants + ordered list, matching PubPlots.
BLACK = "#000000"
ORANGE = "#e69f00"
SKY = "#56b4e9"
GREEN = "#009e73"
YELLOW = "#f0e442"
BLUE = "#0072b2"
VERM = "#d55e00"
REDPURPLE = "#cc79a7"

OKABE = [BLACK, ORANGE, SKY, GREEN, YELLOW, BLUE, VERM, REDPURPLE]


def apply_style() -> None:
    """Apply the shared publication rcParams (mirrors PubPlots' nature theme)."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "legend.fontsize": 9,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )
