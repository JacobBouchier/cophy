#!/usr/bin/env python3
"""
iTOL annotation file builders for CoPhy outputs.

Produces:
  * DATASET_COLORSTRIP files - one per taxonomic rank, one per cooccur
    protein (presence yes/no). Each unique label gets a deterministic
    color from a well-spaced HSL hue rotation.
  * DATASET_GRADIENT files - one per (primary x cooccur) pair, with the
    kb distance for each species.

All files are native iTOL text format: a header declaring the dataset
type, then a DATA block. iTOL accepts these via direct drag-and-drop
onto a tree visualization.

Pure-function design: palette generation and text-building have no I/O
so they're unit-testable. File writes happen in build_itol.py.
"""

from __future__ import annotations

import colorsys


# --------------------------------------------------------------------------
# Color palette
# --------------------------------------------------------------------------

def hsl_palette(n: int, saturation: float = 0.65,
                lightness: float = 0.55) -> list[str]:
    """
    Generate n visually distinct hex colors by evenly rotating through
    HSL hue space.

    Works for any n — at 6 colors they're red/yellow/green/cyan/blue/magenta;
    at 60 they're close-but-distinguishable; at 600 they're crowded but
    still deterministically ordered.

    Returns list of hex strings like ['#A2FF5D', ...] in uppercase, with
    a leading '#', so they paste straight into iTOL files.

    Deterministic: the same n always produces the same colors in the same
    order, so reruns don't reshuffle the visualization.
    """
    if n <= 0:
        return []
    out: list[str] = []
    for i in range(n):
        hue = (i / n) % 1.0
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        out.append(f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}")
    return out


def assign_colors(labels: list[str]) -> dict[str, str]:
    """
    Map each distinct label in `labels` to a unique hex color from the
    HSL palette. Preserves first-seen order for determinism.

    'unclassified' (case-insensitive) gets a neutral grey, so that gap
    in the tree reads visually as 'no data' rather than as 'a particular
    group that happens to look grey'.
    """
    seen: list[str] = []
    for lab in labels:
        if lab not in seen:
            seen.append(lab)
    # Pull 'unclassified' out and assign it grey; colour everything else.
    colored = [lab for lab in seen if lab.lower() != "unclassified"]
    palette = hsl_palette(len(colored))
    mapping = dict(zip(colored, palette))
    for lab in seen:
        if lab.lower() == "unclassified":
            mapping[lab] = "#BDBDBD"
    return mapping


# --------------------------------------------------------------------------
# iTOL DATASET_COLORSTRIP builder
# --------------------------------------------------------------------------

def build_colorstrip(dataset_label: str,
                     rows: list[tuple[str, str, str]],
                     legend_color_to_label: dict[str, str] | None = None,
                     ) -> str:
    """
    Build a DATASET_COLORSTRIP text file for iTOL.

    rows: list of (tree_node_id, color_hex, label) tuples.

    legend_color_to_label: optional dict mapping color -> label for the
    legend block. If omitted, derived from rows. Pass explicitly when you
    want labels not appearing in rows (e.g. an empty result set) still
    listed.

    The output is a single string (one file's worth of text) which the
    caller writes to disk with whatever filename.
    """
    if legend_color_to_label is None:
        legend_color_to_label = {}
        for _node, color, label in rows:
            legend_color_to_label.setdefault(color, label)

    # Build legend arrays in stable order.
    legend_pairs = sorted(legend_color_to_label.items(), key=lambda x: x[1])
    legend_colors = [c for c, _ in legend_pairs]
    legend_labels = [l for _, l in legend_pairs]
    # iTOL wants LEGEND_SHAPES as an integer per legend entry; 1=square.
    legend_shapes = ["1"] * len(legend_pairs)

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{dataset_label}",
        "COLOR\t#000000",
        "STRIP_WIDTH\t25",
        "MARGIN\t0",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t0",
    ]
    if legend_pairs:
        lines.append(f"LEGEND_TITLE\t{dataset_label}")
        lines.append("LEGEND_SHAPES\t" + "\t".join(legend_shapes))
        lines.append("LEGEND_COLORS\t" + "\t".join(legend_colors))
        lines.append("LEGEND_LABELS\t" + "\t".join(legend_labels))
    lines.append("DATA")
    for node, color, label in rows:
        # Lines after DATA: <id> <tab> <color> <tab> <label>
        lines.append(f"{node}\t{color}\t{label}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# iTOL DATASET_GRADIENT builder
# --------------------------------------------------------------------------

def build_gradient(dataset_label: str,
                   rows: list[tuple[str, str]],
                   color_min: str = "#FFFFFF",
                   color_max: str = "#1976D2",
                   ) -> str:
    """
    Build a DATASET_GRADIENT text file for iTOL.

    rows: list of (tree_node_id, value_string) tuples. value_string may
    be empty - iTOL renders no shading for those nodes, which is what we
    want for species lacking distance data.

    color_min and color_max define the gradient endpoints.
    """
    lines = [
        "DATASET_GRADIENT",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{dataset_label}",
        "COLOR\t#000000",
        f"COLOR_MIN\t{color_min}",
        f"COLOR_MAX\t{color_max}",
        "STRIP_WIDTH\t25",
        "MARGIN\t0",
        "SHOW_INTERNAL\t0",
        "DATA",
    ]
    for node, value in rows:
        if value == "" or value is None:
            # iTOL allows skipping rows with no data — we just omit them.
            continue
        lines.append(f"{node}\t{value}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Helpers for the driver
# --------------------------------------------------------------------------

def encode_species_id(species_name: str) -> str:
    """
    Match the encoding used in tree leaf labels: spaces -> underscores.
    Mirrors the pipeline's `encode_species` so iTOL files use the same
    leaf names the tree was built with.
    """
    return species_name.replace(" ", "_")


def parse_kb_value(s: str) -> str | None:
    """
    Take a value from a *_dist_kb column and return either a numeric
    string suitable for the iTOL gradient, or None.

    The kb column can be:
        '' (blank)         -> None (no data; skipped from gradient)
        '0.5', '3.0', etc. -> '0.5'
        '>100'             -> '100'  (capped: iTOL still gets a value)
    """
    s = (s or "").strip()
    if not s:
        return None
    if s.startswith(">"):
        s = s[1:].strip()
    try:
        float(s)
    except ValueError:
        return None
    return s
