#!/usr/bin/env python3
"""
Build iTOL annotation files from the cooccurrence CSV + the species tree.

Produces in <outdir>/itol/:
  * <rank>_color_strip.txt        one per phylum/class/order/family
  * <cooccur_id>_color_strip.txt  one per cooccur protein (presence yes/no)
  * <cooccur>_via_<primary>_distance_color_gradient.txt
                                  one per (primary, cooccur) pair, kb distance
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itol import (    # noqa: E402
    assign_colors, build_colorstrip, build_gradient,
    encode_species_id, parse_kb_value,
)


TAXONOMY_RANKS = ("phylum", "class", "order", "family")


def read_tree_leaves(newick_path: str) -> list[str]:
    """
    Extract leaf node IDs from a Newick tree file. Tokens are anything
    that isn't a control char or whitespace, before the first ':' (which
    starts a branch length) or ')' or ','.
    """
    with open(newick_path) as fh:
        text = fh.read()
    # Strip branch lengths / internal-node support labels and walk the
    # token list. A leaf name is any name not immediately preceded by ')'.
    tokens = re.findall(r"[A-Za-z0-9_.\-]+", text.split(";", 1)[0])
    # The Newick parser approach: leaf names appear immediately after '(' or ','
    # In our trees they are alphanumeric/underscored.
    # Simpler robust approach: match identifiers, then exclude those that
    # appear right after a ')' (those are internal-node support labels).
    leaves: list[str] = []
    i = 0
    cleaned = text.split(";", 1)[0]
    n = len(cleaned)
    while i < n:
        c = cleaned[i]
        if c in "(,":
            # The next identifier (skipping whitespace) is a leaf name
            # IF the character after it is not '(' (which would mean it
            # is an internal subtree, not a leaf).
            j = i + 1
            while j < n and cleaned[j].isspace():
                j += 1
            if j < n and cleaned[j] != "(":
                m = re.match(r"[A-Za-z0-9_.\-]+", cleaned[j:])
                if m:
                    leaves.append(m.group(0))
        i += 1
    return leaves


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True,
                    help="Newick tree file (species_tree.nwk)")
    ap.add_argument("--cooccurrence-csv", required=True)
    ap.add_argument("--collapsed-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cooccur-ids", nargs="*", default=[],
                    help="cooccur ids in the order they should be processed")
    ap.add_argument("--primary-groups", nargs="*", default=[])
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    tree_leaves = read_tree_leaves(args.tree)
    sys.stderr.write(f"itol: {len(tree_leaves)} tree leaves\n")

    # Read the collapsed CSV (one row per species, taxonomy + presence + dist).
    with open(args.collapsed_csv) as fh:
        collapsed_rows = list(csv.DictReader(fh))
    sys.stderr.write(f"itol: {len(collapsed_rows)} species in collapsed CSV\n")

    # Species name -> tree node id mapping. Tree leaves use underscores.
    sp_to_node: dict[str, str] = {}
    for sp in (r["species"] for r in collapsed_rows):
        sp_to_node[sp] = encode_species_id(sp)

    # --- 1. Taxonomy color strips, one per rank. ---
    for rank in TAXONOMY_RANKS:
        if rank not in (collapsed_rows[0].keys() if collapsed_rows else []):
            sys.stderr.write(f"itol: skipping {rank} (no column in CSV)\n")
            continue
        labels_per_species = [r.get(rank, "unclassified") for r in collapsed_rows]
        color_map = assign_colors(labels_per_species)
        strip_rows: list[tuple[str, str, str]] = []
        for r in collapsed_rows:
            label = r.get(rank, "unclassified")
            node = sp_to_node[r["species"]]
            strip_rows.append((node, color_map[label], label))
        text = build_colorstrip(
            dataset_label=rank.capitalize(),
            rows=strip_rows,
            legend_color_to_label={v: k for k, v in color_map.items()},
        )
        out_path = os.path.join(args.out_dir, f"{rank.capitalize()}_color_strip.txt")
        with open(out_path, "w") as fh:
            fh.write(text)
        sys.stderr.write(f"itol: wrote {out_path} "
                         f"({len(color_map)} unique labels)\n")

    # --- 2. Presence color strips, one per primary + one per cooccur id. ---
    # The collapsed CSV has yes/no presence columns for both primaries and
    # cooccur ids, all in the same simple format, so the same loop handles
    # both. Each one gets its own <id>_color_strip.txt file.
    presence_ids = list(args.primary_groups) + list(args.cooccur_ids)
    for pid in presence_ids:
        if pid not in (collapsed_rows[0].keys() if collapsed_rows else []):
            sys.stderr.write(f"itol: skipping {pid} presence (no column)\n")
            continue
        # Two-label palette: yes (saturated), no (grey-white).
        color_map = {"yes": "#1976D2", "no": "#E0E0E0"}
        strip_rows = []
        for r in collapsed_rows:
            presence = r.get(pid, "no").strip().lower() or "no"
            strip_rows.append((sp_to_node[r["species"]],
                               color_map.get(presence, "#E0E0E0"),
                               presence))
        text = build_colorstrip(
            dataset_label=f"{pid} presence",
            rows=strip_rows,
            legend_color_to_label={v: k for k, v in color_map.items()},
        )
        out_path = os.path.join(args.out_dir, f"{pid}_color_strip.txt")
        with open(out_path, "w") as fh:
            fh.write(text)
        sys.stderr.write(f"itol: wrote {out_path}\n")

    # --- 3. Distance gradients, one per (primary, cooccur) pair. ---
    primary_groups = args.primary_groups
    for cid in args.cooccur_ids:
        for pg in primary_groups:
            col = f"{cid}_dist_kb_via_{pg}"
            if not collapsed_rows or col not in collapsed_rows[0].keys():
                sys.stderr.write(f"itol: skipping {col} (no column)\n")
                continue
            # Build rows for every species in the collapsed CSV. iTOL ignores
            # nodes without data, so blank values mean "no shading".
            grad_rows: list[tuple[str, str]] = []
            n_valued = 0
            for r in collapsed_rows:
                v = parse_kb_value(r.get(col, ""))
                grad_rows.append((sp_to_node[r["species"]], v or ""))
                if v:
                    n_valued += 1
            text = build_gradient(
                dataset_label=f"{cid} kb from {pg}",
                rows=grad_rows,
            )
            out_path = os.path.join(
                args.out_dir,
                f"{cid}_via_{pg}_distance_color_gradient.txt")
            with open(out_path, "w") as fh:
                fh.write(text)
            sys.stderr.write(
                f"itol: wrote {out_path} ({n_valued} species with distance)\n")

    sys.stderr.write("itol: done\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
