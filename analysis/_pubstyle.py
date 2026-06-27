#!/usr/bin/env python3
"""Shared publication-figure style for the cell-size-control analysis figures.

A matplotlib port of the science-space **PubPlots** philosophy (its
``_BADNESS_WEIGHTS`` vocabulary, in priority order):

  * accessibility (weight 5): the Okabe-Ito colorblind-safe categorical palette;
    sequential = viridis, diverging = vik (Crameri) -- both perceptually uniform
    and grayscale-survivable. NEVER rainbow/jet/magma/RdBu/tab.
  * redundant encoding: pair colour with a marker shape (MARKER_CYCLE) so groups
    read in grayscale and for colourblind readers.
  * no occlusion (weight 4): legends framed opaque in a clear corner
    (``opaque_legend``); in-plot text haloed so gridlines/data can't cut it
    (``halo``); no text sitting on data.
  * legible craft: font floors, no faint/hairline/tiny marks, big legend markers.
  * clean arrows: ``pub_arrow`` draws a filled-triangle annotation arrow
    (PubPlots ``flow_arrow!``), not matplotlib's thin caret default.

``pub_audit(fig)`` is the matplotlib analog of PubPlots ``publication_audit`` --
it returns the violations (text/legend over data, sub-floor fonts, non-CB-safe
colormaps) so a plot script can assert a figure is clean before saving.

Usage:
    from _pubstyle import (apply_style, pub_arrow, halo, opaque_legend,
                           SEQ_CMAP, DIV_CMAP, MARKER_CYCLE, pub_audit,
                           BLUE, VERM, GREEN, OKABE)
    apply_style()
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.patheffects as _pe

# --- Okabe-Ito colorblind-safe categorical palette (Okabe & Ito 2008) ---------
BLACK = "#000000"
ORANGE = "#e69f00"
SKY = "#56b4e9"
GREEN = "#009e73"
YELLOW = "#f0e442"
BLUE = "#0072b2"
VERM = "#d55e00"
REDPURPLE = "#cc79a7"

OKABE = [BLACK, ORANGE, SKY, GREEN, YELLOW, BLUE, VERM, REDPURPLE]

# Redundant shape encoding, same order as OKABE (PubPlots MARKER_CYCLE): so a
# series is told apart by BOTH colour and glyph -- robust in grayscale / for CVD.
MARKER_CYCLE = ["o", "s", "^", "D", "v", "*", "P", "X"]

# Colormaps that survive grayscale + colour-vision deficiency (PubPlots
# SEQ_COLORMAP=:viridis, DIV_COLORMAP=:vik). vik needs cmcrameri; fall back to a
# CB-safe built-in diverging if it is not installed.
SEQ_CMAP = "viridis"
try:  # Crameri's vik, the portfolio diverging standard
    import cmcrameri.cm as _cmc  # noqa: F401

    DIV_CMAP = _cmc.vik
except Exception:  # CB-safe blue-white-red fallback (NOT RdBu/seismic)
    DIV_CMAP = "coolwarm"

# colormaps allowed by the audit (perceptually-uniform / single-hue, CB-safe)
_GOOD_CMAPS = {
    "viridis", "cividis", "vik", "Blues", "Greys", "gray", "Greens", "Oranges",
    "Purples", "coolwarm", "bwr_cb",
}


def apply_style() -> None:
    """Apply the shared publication rcParams (mirrors PubPlots' nature theme)."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.prop_cycle": plt.cycler(color=OKABE[1:]),  # skip black for lines
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            # opaque, legible legends by default (occlusion guard)
            "legend.frameon": True,
            "legend.framealpha": 0.93,
            "legend.facecolor": "white",
            "legend.edgecolor": "0.8",
            "legend.fontsize": 11,
            "legend.markerscale": 1.6,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )


def pub_arrow(ax, xy, xytext, text="", *, color="0.25", lw=1.6, scale=14,
              shrinkA=3, shrinkB=4, connectionstyle=None, **kw):
    """A clean filled-triangle annotation arrow (PubPlots ``flow_arrow!``).

    Points FROM ``xytext`` (where the label sits) TO ``xy`` (the target). Uses
    ``arrowstyle='-|>'`` (a solid triangle head) at a visible ``scale``, instead
    of matplotlib's default thin two-stroke caret. ``shrinkB`` keeps the head off
    the target marker so it isn't buried under it.
    """
    aprops = dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=scale,
                  shrinkA=shrinkA, shrinkB=shrinkB)
    if connectionstyle:
        aprops["connectionstyle"] = connectionstyle
    return ax.annotate(text, xy=xy, xytext=xytext, arrowprops=aprops,
                       color=color, **kw)


def halo(artist, lw=2.2, fg="white"):
    """White stroke behind text so gridlines/data can't cut the glyphs
    (PubPlots ``halo_text!`` == matplotlib ``path_effects.withStroke``)."""
    arts = artist if isinstance(artist, (list, tuple)) else [artist]
    for a in arts:
        if a is not None:
            a.set_path_effects([_pe.withStroke(linewidth=lw, foreground=fg)])
    return artist


def opaque_legend(ax, *, loc="best", markerscale=1.7, **kw):
    """Framed, opaque legend with legible markers, masking whatever sits behind
    it (PubPlots ``legend_opaque!``). Default ``loc='best'`` lets matplotlib pick
    the emptiest spot; pass an explicit corner to override."""
    leg = ax.legend(loc=loc, frameon=True, framealpha=0.93, facecolor="white",
                    edgecolor="0.8", markerscale=markerscale, **kw)
    if leg is not None:
        leg.set_zorder(6)
    return leg


def clean_hexbin(ax, x, y, **kwargs):
    """A ``hexbin`` that tiles cleanly: each hexagon's edge in its own face
    colour fills matplotlib's antialiasing seams. Any kwarg passes through."""
    kwargs.setdefault("edgecolors", "face")
    kwargs.setdefault("linewidths", 0.4)
    kwargs.setdefault("cmap", SEQ_CMAP)
    return ax.hexbin(x, y, **kwargs)


def pub_audit(fig, *, min_font=8.0, occlude_frac=0.12):
    """matplotlib analog of PubPlots ``publication_audit`` -- return the list of
    violations (empty == clean): in-plot text or legend occluding data, sub-floor
    fonts, and non-CB-safe colormaps. Call after the figure is fully built."""
    import numpy as np

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    issues = []
    for ax in fig.axes:
        # data points in display coords (lines + collection offsets + bar corners)
        chunks = []
        for ln in ax.get_lines():
            xy = ln.get_xydata()
            if len(xy) >= 2:  # densify so a curve CROSSING a box (no vertex in it) is caught
                seg = np.linspace(0, 1, 8)[:, None]
                dense = np.concatenate(
                    [xy[i] + seg * (xy[i + 1] - xy[i]) for i in range(len(xy) - 1)])
                chunks.append(ax.transData.transform(dense))
            elif len(xy):
                chunks.append(ax.transData.transform(xy))
        for col in ax.collections:
            try:
                off = col.get_offsets()
                if len(off):
                    chunks.append(ax.transData.transform(off))
            except Exception:
                pass
        for p in getattr(ax, "patches", []):
            try:
                chunks.append(p.get_window_extent(r).get_points())
            except Exception:
                pass
        pts = np.vstack(chunks) if chunks else np.empty((0, 2))

        # candidate occluders: in-plot text + the legend box
        boxes = [("text '%s'" % t.get_text()[:24], t) for t in ax.texts
                 if t.get_text().strip()]
        leg = ax.get_legend()
        if leg is not None:
            boxes.append(("legend", leg))
        for name, art in boxes:
            try:
                bb = art.get_window_extent(r)
            except Exception:
                continue
            if len(pts):
                frac = float(np.mean(
                    (pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1)
                    & (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1)))
                if frac > occlude_frac:
                    issues.append("%s: %s occludes %.0f%% of data"
                                  % (ax.get_title()[:24] or "axes", name, 100 * frac))
        # font floors
        for t in ax.texts:
            if t.get_text().strip() and t.get_fontsize() < min_font:
                issues.append("sub-floor font %.0fpt: '%s'"
                              % (t.get_fontsize(), t.get_text()[:24]))
        # colormaps. NB: a contour drawn with an explicit ``colors=`` argument still carries a
        # default placeholder colormap named "from_list" -- it is not a real colour choice (the
        # lines are the given colours), so it is exempt; only a genuinely-set cmap is audited.
        for im in list(ax.images) + list(ax.collections):
            cm = getattr(im, "get_cmap", lambda: None)()
            nm = getattr(cm, "name", "") or ""
            if nm == "from_list":  # placeholder from an explicit-colour contour / LineCollection
                continue
            if nm and nm not in _GOOD_CMAPS and not nm.startswith("cmc"):
                issues.append("non-CB-safe colormap '%s' (use viridis/vik)" % nm)
    return issues
