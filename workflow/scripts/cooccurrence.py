#!/usr/bin/env python3
"""
Build the co-occurrence CSV (step 9).

Inputs (in memory):
  * primary_members: for each primary query group, the deduped list of
    ConsensusCall representatives (one per species) that made it into the
    tree. We use these to enumerate (species, group, ortholog-accession).
  * cooccur_species: for each co-occurrence query id, the SET of species
    keys in which that protein has an ortholog (by the same search+consensus
    method). Presence is membership in that set.

Output rows (one per (species, primary_group)):
    species_name, primary_query_id, ortholog_accession,
    <cooccur_id_1>, <cooccur_id_2>, ...        # each "yes"/"no"

Row granularity and the two identifier columns match the agreed spec.
Pure functions; fully unit-tested. CSV writing uses the stdlib csv module
so quoting/escaping is correct for species names containing commas.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass


@dataclass
class PrimaryMember:
    species: str          # species key (row identity)
    primary_group: str    # the query UniProt id defining the group
    ortholog_acc: str     # the representative hit's UniProt accession


def build_rows(primary_members: list[PrimaryMember],
               cooccur_ids: list[str],
               cooccur_species: dict[str, set[str]],
               primary_ids: list[str] | None = None,
               primary_species: dict[str, set[str]] | None = None,
               ) -> list[dict[str, str]]:
    """
    Produce one row per (species, primary_group). The yes/no columns report
    whether each co-occurrence protein is present in that row's species.

    If primary_ids and primary_species are given, additional yes/no columns
    are emitted for each primary query id (using the same presence convention
    as the cooccur columns). This lets users see, on a row for species X under
    primary group A, whether species X also has primary B, C, ... — useful
    when there are multiple primaries and you want a complete presence matrix.

    Output is sorted by (species, primary_group) for determinism.
    """
    primary_ids = primary_ids or []
    primary_species = primary_species or {}
    rows: list[dict[str, str]] = []
    for m in sorted(primary_members,
                    key=lambda x: (x.species, x.primary_group, x.ortholog_acc)):
        row = {
            "species": m.species,
            "primary_query_id": m.primary_group,
            "ortholog_accession": m.ortholog_acc,
        }
        # Primary presence columns (one per primary id).
        for pid in primary_ids:
            present = m.species in primary_species.get(pid, set())
            row[pid] = "yes" if present else "no"
        # Cooccur presence columns (one per cooccur id).
        for cid in cooccur_ids:
            present = m.species in cooccur_species.get(cid, set())
            row[cid] = "yes" if present else "no"
        rows.append(row)
    return rows


def write_csv(rows: list[dict[str, str]], cooccur_ids: list[str],
              path: str) -> None:
    header = ["species", "primary_query_id", "ortholog_accession"] + list(cooccur_ids)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def collapse_by_species(
    rows: list[dict[str, str]],
    primary_groups: list[str],
    cooccur_ids: list[str],
) -> list[dict[str, str]]:
    r"""
    Pivot per-row data into one row per species (wide format).

    Input rows are the detailed format: one row per (species, primary_group).
    Each row carries primary_query_id, ortholog_accession, and per cooccur id
    both a yes/no column and (when neighborhood is enabled) an _dist column.

    Output row layout (one per species, columns in this order):
        species
        <primary_group>_ortholog                 \   one pair per primary
        ...                                      /
        <cooccur_id>                             \   one block per cooccur:
        <cooccur_id>_dist_via_<primary_group>     |  presence OR'd across rows,
        <cooccur_id>_dist_via_<primary_group>     |  per-primary distance from
        ...                                      /   that primary's row
        ...

    Collapse rules:
      * <primary>_ortholog: the accession from the row whose primary_query_id
        matches that primary, or empty if the species has no such row.
      * <cooccur>: 'yes' if any row for this species had 'yes', else 'no'.
      * <cooccur>_dist_via_<primary>: distance from this species' row for
        that primary (empty if the species has no row for that primary,
        OR if that row's distance was already blank).
    """
    # Group rows by species, preserving first-seen order for stable output.
    by_species: dict[str, list[dict[str, str]]] = {}
    order: list[str] = []
    for r in rows:
        sp = r["species"]
        if sp not in by_species:
            by_species[sp] = []
            order.append(sp)
        by_species[sp].append(r)

    out: list[dict[str, str]] = []
    # Recognise taxonomy columns if present (assume all rows for a species
    # carry the same lineage; pick from the first).
    tax_keys = ("phylum", "class", "order", "family")
    for sp in order:
        sp_rows = by_species[sp]
        # Index this species' rows by their primary group, for lookup.
        by_primary = {r["primary_query_id"]: r for r in sp_rows}

        new: dict[str, str] = {"species": sp}

        # Taxonomy columns (carried through from any row for this species).
        first = sp_rows[0]
        for k in tax_keys:
            if k in first:
                new[k] = first[k]

        # One ortholog column per primary group, plus a yes/no presence
        # column (parallel to the cooccur columns). Presence is OR'd across
        # any rows for this species (in practice a species has at most one
        # row per primary, so this is just "did any row for this primary
        # mark its primary as yes").
        for pg in primary_groups:
            r = by_primary.get(pg)
            new[f"{pg}_ortholog"] = r["ortholog_accession"] if r else ""
            # Primary presence: 'yes' if any row for this species said yes
            # for this primary id (added by build_rows when primary_ids
            # were passed). Falls back to inferring from ortholog_accession
            # being non-empty.
            presence_pg = "no"
            for rr in sp_rows:
                if rr.get(pg, "").strip().lower() == "yes":
                    presence_pg = "yes"
                    break
            else:
                # Inference fallback: if this species has an actual ortholog
                # row for this primary, treat it as 'yes' even without an
                # explicit column. This keeps behavior sensible for older
                # detailed CSVs that lack the presence columns.
                if r is not None and r.get("ortholog_accession", "").strip():
                    presence_pg = "yes"
            new[pg] = presence_pg

        # For each cooccur id: presence OR'd across rows, then one distance
        # column per primary group (read from that primary's row), plus a
        # corresponding kb column.
        for cid in cooccur_ids:
            presence = "no"
            for r in sp_rows:
                if r.get(cid, "no").strip().lower() == "yes":
                    presence = "yes"
                    break
            new[cid] = presence
            dist_key = f"{cid}_dist"
            dist_kb_key = f"{cid}_dist_kb"
            for pg in primary_groups:
                r = by_primary.get(pg)
                new[f"{cid}_dist_via_{pg}"] = (
                    r.get(dist_key, "") if r else ""
                )
                new[f"{cid}_dist_kb_via_{pg}"] = (
                    r.get(dist_kb_key, "") if r else ""
                )
        out.append(new)
    return out


def write_collapsed_csv(
    collapsed: list[dict[str, str]],
    primary_groups: list[str],
    cooccur_ids: list[str],
    path: str,
) -> None:
    """Write the one-row-per-species (wide) CSV with the agreed column order.
    If the rows carry taxonomy columns (phylum/class/order/family), they
    appear right after `species`."""
    header = ["species"]
    # Detect taxonomy columns by sampling the first row.
    if collapsed:
        for k in ("phylum", "class", "order", "family"):
            if k in collapsed[0]:
                header.append(k)
    for pg in primary_groups:
        header.append(f"{pg}_ortholog")
        header.append(pg)  # primary presence yes/no
    for cid in cooccur_ids:
        header.append(cid)
        for pg in primary_groups:
            header.append(f"{cid}_dist_via_{pg}")
            header.append(f"{cid}_dist_kb_via_{pg}")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        for row in collapsed:
            writer.writerow(row)
