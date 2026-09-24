#!/usr/bin/env python3
"""
Assemble the co-occurrence CSV (step 9) from:
  * the per-primary-group members.tsv files (species + ortholog accession)
  * the per-cooccurrence-id deduped member files (which species contain it)

For each co-occurrence id we run the SAME search+consensus+dedup machinery
(via build_group.py upstream), producing a members file per cooccur id whose
'species' column is exactly the set of species containing that protein.

When [neighborhood] neighborhood=true, also computes the signed gene-order
distance between each primary ortholog and each present co-occurrence
ortholog, adding a <cooccur_id>_dist column.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cooccurrence import (PrimaryMember, build_rows, write_csv,    # noqa: E402
                          collapse_by_species, write_collapsed_csv)
from fetch import decode_species                                   # noqa: E402


def read_members(path: str):
    """Return list of (species_display, accession) from a members.tsv.
    Species column was written in encoded (underscored) form; decode for
    human-readable CSV output."""
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out.append((decode_species(f[0]), f[1]))
    return out


def read_species_set(path: str):
    return {sp for sp, _ in read_members(path)}


def read_species_to_acc(path: str) -> dict[str, str]:
    """For neighborhood: species_display -> accession of that species' hit."""
    return {sp: acc for sp, acc in read_members(path)}


def read_species_to_taxid(path: str) -> dict[str, str]:
    """For taxonomy: species_display -> taxon_id (column 3) for each row.
    Returns {} if the file lacks a taxon_id column."""
    out: dict[str, str] = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2].strip():
                out[decode_species(f[0])] = f[2].strip()
    return out


def _compute_distances(primary_members, cooccur_id_to_species_acc,
                       cooccur_ids, cfg, max_dist):
    """
    Returns {(species, primary_group, cooccur_id): formatted_distance_str}
    using neighborhood.closest_signed_distance + format_distance.
    Network-heavy; called only when neighborhood=true.
    """
    from fetch_neighborhood import (resolve_positions, configure,    # noqa: E402
                                     emit_warn_summary, prefetch_idmapping)
    from neighborhood import (closest_signed_distance, format_distance,
                              format_distance_kb)

    configure(api_key=cfg.get("ncbi_api_key", ""),
              email=cfg.get("ncbi_email", ""),
              max_lookups=cfg.get("neighborhood_max_lookups_per_run", 0))

    cache_dir = cfg["cache_dir"]
    # Collect all accessions whose positions we'll need.
    needed = {m.ortholog_acc for m in primary_members}
    for cid in cooccur_ids:
        needed.update(cooccur_id_to_species_acc.get(cid, {}).values())
    sys.stderr.write(
        f"neighborhood: resolving positions for {len(needed)} accessions "
        f"via NCBI EFetch\n")

    # Batch idmapping pre-pass: recover UniProt -> RefSeq mappings the
    # inline xrefs miss. Non-fatal if it fails.
    prefetch_idmapping(sorted(needed), cache_dir)

    positions: dict[str, list] = {}
    for i, acc in enumerate(sorted(needed)):
        if i and i % 25 == 0:
            sys.stderr.write(f"neighborhood: {i}/{len(needed)} resolved\n")
        try:
            positions[acc] = resolve_positions(acc, cache_dir)
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"WARN: resolve_positions {acc}: {exc}\n")
            positions[acc] = []

    n_with_pos = sum(1 for p in positions.values() if p)
    sys.stderr.write(
        f"neighborhood: {n_with_pos}/{len(positions)} accessions resolved "
        f"to coordinates\n")
    emit_warn_summary()

    out: dict[tuple[str, str, str], tuple[str, str]] = {}
    for m in primary_members:
        for cid in cooccur_ids:
            cooccur_acc = cooccur_id_to_species_acc.get(cid, {}).get(m.species)
            if not cooccur_acc:
                out[(m.species, m.primary_group, cid)] = ("", "")
                continue
            d = closest_signed_distance(positions.get(m.ortholog_acc, []),
                                        positions.get(cooccur_acc, []))
            out[(m.species, m.primary_group, cid)] = (
                format_distance(d, max_dist),
                format_distance_kb(d, max_dist),
            )
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--primary-members", nargs="+", required=True,
                    help="members.tsv per primary group; basename = query id")
    ap.add_argument("--cooccur-members", nargs="*", default=[],
                    help="members.tsv per cooccur id (legacy search path)")
    ap.add_argument("--cooccur-rbh", nargs="*", default=[],
                    help="RBH TSV per cooccur id (cooccur_rbh/<cid>.tsv)")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-collapsed-csv", default="",
                    help="Optional path for the one-row-per-species pivot CSV")
    args = ap.parse_args()

    from config import load_config
    cfg = load_config(args.config)

    # Primary members -> PrimaryMember rows, and per-primary species sets
    # (mirrors the cooccur_species structure so build_rows can emit primary
    # presence columns).
    primary: list[PrimaryMember] = []
    primary_ids: list[str] = []
    primary_species: dict[str, set] = {}
    for path in args.primary_members:
        group = os.path.basename(path).split(".")[0]
        if group not in primary_ids:
            primary_ids.append(group)
            primary_species[group] = set()
        for species, acc in read_members(path):
            primary.append(PrimaryMember(species=species,
                                         primary_group=group,
                                         ortholog_acc=acc))
            primary_species[group].add(species)

    # Cooccurrence presence sets + species->acc map per cooccur id.
    cooccur_ids: list[str] = []
    cooccur_species: dict[str, set] = {}
    cooccur_species_acc: dict[str, dict[str, str]] = {}

    # Build a taxon_id -> species_display map by scanning primary members.
    # We need this to convert RBH TSVs (taxon_id-keyed) back to species
    # strings (used as the join key in the rest of the pipeline).
    taxid_to_species: dict[str, str] = {}
    for path in args.primary_members:
        for sp, tid in read_species_to_taxid(path).items():
            taxid_to_species.setdefault(tid, sp)

    # Path A: classic per-tool search (cooccur_members points at members.tsv
    # files in results/groups/, same format as primaries).
    for path in args.cooccur_members:
        cid = os.path.basename(path).split(".")[0]
        cooccur_ids.append(cid)
        sp_to_acc = read_species_to_acc(path)
        cooccur_species[cid] = set(sp_to_acc.keys())
        cooccur_species_acc[cid] = sp_to_acc

    # Path B: RBH TSV (cooccur_rbh points at results/cooccur_rbh/<cid>.tsv).
    # Schema: taxon_id\treciprocal_tools\tforward_hit\tforward_evalue
    # Species are called present when reciprocal_tools is non-empty.
    for path in getattr(args, "cooccur_rbh", None) or []:
        cid = os.path.basename(path).split(".")[0]
        if cid in cooccur_ids:
            # Both files given for the same id; RBH wins (richer info).
            cooccur_species[cid].clear()
            cooccur_species_acc[cid].clear()
        else:
            cooccur_ids.append(cid)
            cooccur_species[cid] = set()
            cooccur_species_acc[cid] = {}
        with open(path) as fh:
            fh.readline()  # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                tid, recip_tools, fwd_hit, _fwd_e = f[:4]
                if not recip_tools.strip():
                    continue  # not present in this species
                sp = taxid_to_species.get(tid.strip())
                if sp is None:
                    # Species not in tree (shouldn't normally happen,
                    # but be defensive). Skip.
                    continue
                cooccur_species[cid].add(sp)
                cooccur_species_acc[cid][sp] = fwd_hit.strip()

    rows = build_rows(primary, cooccur_ids, cooccur_species,
                      primary_ids=primary_ids,
                      primary_species=primary_species)

    # Build species -> taxon_id map by scanning every members file.
    # We use the first taxon_id we see for a species (they should all match).
    species_to_taxid: dict[str, str] = {}
    for path in list(args.primary_members) + list(args.cooccur_members):
        for sp, tid in read_species_to_taxid(path).items():
            species_to_taxid.setdefault(sp, tid)

    # Optionally look up phylum/class/order/family for each species.
    taxonomy_added = False
    if cfg.get("taxonomy", True) and species_to_taxid:
        try:
            from taxonomy import fetch_lineages, RANKS, EMPTY_LINEAGE
            unique_taxids = sorted({t for t in species_to_taxid.values() if t})
            sys.stderr.write(
                f"taxonomy: resolving lineages for {len(unique_taxids)} "
                f"unique taxon_ids\n")
            taxid_to_lineage = fetch_lineages(unique_taxids, cfg["cache_dir"])
            for row in rows:
                tid = species_to_taxid.get(row["species"], "")
                ranks = taxid_to_lineage.get(tid, dict(EMPTY_LINEAGE))
                for r in RANKS:
                    row[r] = ranks.get(r, "unclassified")
            taxonomy_added = True
            sys.stderr.write(
                f"taxonomy: added {list(RANKS)} columns to {len(rows)} rows\n")
        except Exception as exc:  # noqa: BLE001 — must not break the rule
            sys.stderr.write(f"WARN: taxonomy lookup failed: {exc}\n")

    # Optionally compute distances and merge into rows.
    header_extra: list[str] = []
    if cfg.get("neighborhood") and cooccur_ids:
        dist_map = _compute_distances(
            primary, cooccur_species_acc, cooccur_ids,
            cfg, cfg["neighborhood_max_dist"])
        for row in rows:
            for cid in cooccur_ids:
                key = (row["species"], row["primary_query_id"], cid)
                bp, kb = dist_map.get(key, ("", ""))
                row[f"{cid}_dist"] = bp
                row[f"{cid}_dist_kb"] = kb
        header_extra = []
        for cid in cooccur_ids:
            header_extra += [f"{cid}_dist", f"{cid}_dist_kb"]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    # Interleave each cooccur id with its _dist and _dist_kb columns.
    if header_extra:
        import csv as _csv
        header = ["species"]
        if taxonomy_added:
            from taxonomy import RANKS as _TAX_RANKS
            header += list(_TAX_RANKS)
        header += ["primary_query_id", "ortholog_accession"]
        # Primary presence yes/no columns (one per primary id).
        header += list(primary_ids)
        for cid in cooccur_ids:
            header += [cid, f"{cid}_dist", f"{cid}_dist_kb"]
        with open(args.out_csv, "w", newline="", encoding="utf-8") as fh:
            w = _csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    else:
        # No neighborhood — but we may still have taxonomy. Build header manually.
        import csv as _csv
        header = ["species"]
        if taxonomy_added:
            from taxonomy import RANKS as _TAX_RANKS
            header += list(_TAX_RANKS)
        header += ["primary_query_id", "ortholog_accession"]
        header += list(primary_ids) + list(cooccur_ids)
        with open(args.out_csv, "w", newline="", encoding="utf-8") as fh:
            w = _csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    sys.stderr.write(f"cooccurrence: {len(rows)} rows, "
                     f"{len(cooccur_ids)} cooccur columns"
                     f"{' + distances' if header_extra else ''} "
                     f"-> {args.out_csv}\n")

    # Optionally write the collapsed pivot CSV (one row per species).
    if args.out_collapsed_csv:
        # Determine the primary group order from the input file basenames so
        # the columns appear in the same order the user listed them.
        primary_groups: list[str] = []
        seen: set[str] = set()
        for path in args.primary_members:
            pg = os.path.basename(path).split(".")[0]
            if pg not in seen:
                primary_groups.append(pg)
                seen.add(pg)
        collapsed = collapse_by_species(rows, primary_groups, cooccur_ids)
        os.makedirs(os.path.dirname(os.path.abspath(args.out_collapsed_csv)),
                    exist_ok=True)
        write_collapsed_csv(collapsed, primary_groups, cooccur_ids,
                            args.out_collapsed_csv)
        sys.stderr.write(f"cooccurrence (collapsed): {len(collapsed)} species "
                         f"-> {args.out_collapsed_csv}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
