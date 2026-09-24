#!/usr/bin/env python3
"""
Fetch query data from UniProt / EBI and parse it.

Design notes
------------
* The network layer (`_http_get_json`, `_http_get_text`) is deliberately
  thin and isolated so that the *parsing* functions below can be unit
  tested against saved JSON fixtures with no network access.
* Everything is cached to `cache_dir` keyed by accession, so re-runs and
  Snakemake retries don't re-hit the API.
* We never trust the network to be up: failures raise a clear FetchError
  that the calling rule surfaces, rather than producing silent empties.

What we extract per query UniProt accession:
    - protein sequence (for aa mode and as the mmseqs2/HMMER query)
    - organism scientific name + NCBI taxon id
    - reviewed (Swiss-Prot) flag  -> used by the dedup tie-break
    - EMBL/ENA CDS cross-reference ids (for nt mode)
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict
from typing import Any


UNIPROT_REST = "https://rest.uniprot.org"
EBI_PROTEINS = "https://www.ebi.ac.uk/proteins/api"
USER_AGENT = "CoPhy/1.0 (https://github.com/<your-username>/cophy; research pipeline)"


class FetchError(RuntimeError):
    """Raised when a remote resource cannot be retrieved or parsed."""


@dataclass
class QueryEntry:
    accession: str
    scientific_name: str
    taxon_id: int
    reviewed: bool
    sequence: str
    embl_cds_ids: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# --------------------------------------------------------------------------
# Network layer (isolated, retried)
# --------------------------------------------------------------------------

def _http_get(url: str, *, accept: str, retries: int = 4,
              backoff: float = 2.0, timeout: int = 60) -> bytes:
    last_err: Exception | None = None
    for attempt in range(retries):
        req = urllib.request.Request(
            url, headers={"Accept": accept, "User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            # 4xx (except 429) won't be fixed by retrying.
            if exc.code != 429 and 400 <= exc.code < 500:
                raise FetchError(f"HTTP {exc.code} for {url}") from exc
            last_err = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last_err = exc
        sleep_for = backoff * (2 ** attempt)
        time.sleep(sleep_for)
    raise FetchError(f"Failed to GET {url} after {retries} attempts: {last_err}")


def _http_get_json(url: str, **kw) -> Any:
    raw = _http_get(url, accept="application/json", **kw)
    try:
        return json.loads(raw.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise FetchError(f"Invalid JSON from {url}: {exc}") from exc


def _http_get_text(url: str, **kw) -> str:
    return _http_get(url, accept="text/plain", **kw).decode("utf-8")


def _parse_link_header_next(link_header: str | None) -> str | None:
    """
    Extract the rel="next" URL from an HTTP Link header.

    UniProt's paginated REST endpoints (including idmapping/results) put
    the next-page cursor here, not in the response body. Format:

        Link: <https://rest.uniprot.org/idmapping/results/{job}?cursor=X&size=500>; rel="next"

    Returns None if no Link header or no next URL.
    """
    if not link_header:
        return None
    # The header can contain multiple comma-separated entries; we want
    # the one with rel="next".
    for part in link_header.split(","):
        part = part.strip()
        if 'rel="next"' not in part and "rel=next" not in part:
            continue
        # The URL is wrapped in <...>
        lt = part.find("<")
        gt = part.find(">", lt + 1)
        if lt >= 0 and gt > lt:
            return part[lt + 1:gt]
    return None


def _http_get_json_with_link(url: str, **kw) -> tuple[Any, str | None]:
    """Like _http_get_json but also returns the next-page URL from Link header."""
    last_err: Exception | None = None
    for attempt in range(4):
        req = urllib.request.Request(
            url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                body = resp.read()
                link_header = resp.headers.get("Link")
                payload = json.loads(body.decode("utf-8"))
                return payload, _parse_link_header_next(link_header)
        except urllib.error.HTTPError as exc:
            if exc.code != 429 and 400 <= exc.code < 500:
                raise FetchError(f"HTTP {exc.code} for {url}") from exc
            last_err = exc
        except (urllib.error.URLError, TimeoutError, ValueError) as exc:
            last_err = exc
        time.sleep(2.0 * (2 ** attempt))
    raise FetchError(f"Failed to GET {url}: {last_err}")


# --------------------------------------------------------------------------
# Parsing layer (pure functions; unit-testable against fixtures)
# --------------------------------------------------------------------------

def parse_uniprot_entry(obj: dict[str, Any]) -> QueryEntry:
    """
    Turn a UniProtKB entry JSON object into a QueryEntry.
    Tolerant of missing optional fields; raises FetchError on missing
    essentials (accession or sequence).
    """
    accession = obj.get("primaryAccession")
    if not accession:
        raise FetchError("UniProt entry missing primaryAccession")

    entry_type = obj.get("entryType", "")
    et = entry_type.lower()
    # NB: "unreviewed" contains the substring "reviewed", so test the
    # unreviewed/TrEMBL case explicitly and treat everything else as
    # reviewed only when it actually says swiss-prot/reviewed.
    if "unreviewed" in et or "trembl" in et:
        reviewed = False
    else:
        reviewed = ("swiss-prot" in et) or ("reviewed" in et)

    organism = obj.get("organism", {}) or {}
    scientific_name = organism.get("scientificName", "") or ""
    taxon_id = organism.get("taxonId")
    if taxon_id is None:
        taxon_id = -1

    seq_obj = obj.get("sequence", {}) or {}
    sequence = seq_obj.get("value", "") or ""
    if not sequence:
        raise FetchError(f"{accession}: no sequence in UniProt entry")

    embl_cds: list[str] = []
    for xref in obj.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "EMBL":
            continue
        # Each EMBL xref may carry properties; the CDS protein id is in a
        # property with key "ProteinId". The molecule/CDS linkage is what
        # we map later via the ID-mapping endpoint.
        for prop in xref.get("properties", []) or []:
            if prop.get("key") == "ProteinId":
                val = prop.get("value", "")
                if val and val != "-":
                    embl_cds.append(val)

    # De-dup while preserving order.
    seen: set[str] = set()
    embl_cds = [x for x in embl_cds if not (x in seen or seen.add(x))]

    return QueryEntry(
        accession=accession,
        scientific_name=scientific_name,
        taxon_id=int(taxon_id),
        reviewed=bool(reviewed),
        sequence=sequence,
        embl_cds_ids=embl_cds,
    )


def species_key(scientific_name: str, rank: str) -> str:
    """
    Collapse an organism scientific name to a de-duplication key.
      rank=binomial  -> 'Genus species' (first two whitespace tokens)
      rank=subspecies-> keep up to three tokens (genus species subsp.)
      rank=none      -> the full name unchanged
    Falls back gracefully on single-token or empty names.
    """
    name = " ".join(scientific_name.split())  # normalise whitespace
    if rank == "none" or not name:
        return name
    tokens = name.split(" ")
    if rank == "binomial":
        return " ".join(tokens[:2])
    if rank == "subspecies":
        return " ".join(tokens[:3])
    return name


# FASTA headers must be whitespace-free because every common alignment and
# tree tool (MAFFT, ClipKIT, FastTree, etc.) splits the header at the first
# whitespace and uses only the first token as the sequence id. If we wrote
# headers as "Genus species" the species epithet would silently disappear
# and 113 distinct Corynebacterium species would all collide into a single
# row keyed "Corynebacterium". Encode to underscores when writing FASTAs and
# decode back when presenting to the user (CSV, logs).
_FASTA_UNSAFE_TO_UNDERSCORE = str.maketrans({
    " ": "_", "\t": "_", "(": "_", ")": "_", ",": "_", ";": "_",
    ":": "_", "'": "_", '"': "_",
})


def encode_species(name: str) -> str:
    """Turn a display species name into a FASTA-safe token."""
    return name.translate(_FASTA_UNSAFE_TO_UNDERSCORE)


def decode_species(encoded: str) -> str:
    """Best-effort inverse of encode_species for display.
    We can't perfectly recover which underscore was originally a space vs
    something else, but a single underscore between two alpha tokens is
    overwhelmingly a space — that handles "Genus_species" cleanly. Multiple
    consecutive underscores are left alone so we don't mangle anything."""
    out = []
    i = 0
    while i < len(encoded):
        ch = encoded[i]
        if (ch == "_" and 0 < i < len(encoded) - 1
                and encoded[i - 1].isalpha() and encoded[i + 1].isalpha()):
            out.append(" ")
        else:
            out.append(ch)
        i += 1
    return "".join(out)


# --------------------------------------------------------------------------
# Cached high-level fetchers
# --------------------------------------------------------------------------

def _cache_path(cache_dir: str, name: str) -> str:
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, name)


def fetch_entry(accession: str, cache_dir: str) -> QueryEntry:
    """Fetch + parse a UniProtKB entry, using a JSON cache."""
    cache = _cache_path(cache_dir, f"uniprot_{accession}.json")
    if os.path.exists(cache):
        with open(cache, "r", encoding="utf-8") as fh:
            obj = json.load(fh)
    else:
        url = f"{UNIPROT_REST}/uniprotkb/{urllib.parse.quote(accession)}.json"
        obj = _http_get_json(url)
        with open(cache, "w", encoding="utf-8") as fh:
            json.dump(obj, fh)
    return parse_uniprot_entry(obj)


def fetch_structure(accession: str, cache_dir: str) -> str | None:
    """
    Fetch the AlphaFold predicted structure (PDB) for an accession, for
    foldseek queries. Returns the path to the cached .pdb, or None if no
    model exists (many entries have none).

    Uses the AlphaFold API to discover the current file URL rather than
    hardcoding a model version (which has changed: v4 → v5 over time, and
    will keep changing). The API endpoint returns a JSON with a 'pdbUrl'
    field pointing at the canonical current model file.
    """
    cache = _cache_path(cache_dir, f"AF-{accession}.pdb")
    if os.path.exists(cache):
        return cache

    # Step 1: query the AlphaFold prediction API to get the current pdbUrl.
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    try:
        payload = _http_get_json(api_url)
    except FetchError:
        return None
    if not payload or not isinstance(payload, list) or not payload[0]:
        return None
    pdb_url = payload[0].get("pdbUrl")
    if not pdb_url:
        # Fall back to the historical v4 URL pattern in case the API
        # changes its response shape unexpectedly.
        pdb_url = (f"https://alphafold.ebi.ac.uk/files/"
                   f"AF-{accession}-F1-model_v4.pdb")

    # Step 2: fetch the actual structure file.
    try:
        data = _http_get(pdb_url, accept="text/plain")
    except FetchError:
        return None
    with open(cache, "wb") as fh:
        fh.write(data)
    return cache


def map_uniprot_to_embl_cds(accessions: list[str], cache_dir: str) -> dict[str, list[str]]:
    """
    Use the UniProt ID-mapping endpoint to map UniProtKB AC -> EMBL CDS.
    Returns {accession: [cds_id, ...]}. Used for nt mode. This is a thin
    wrapper that submits a job, polls, then reads results.
    """
    cache = _cache_path(cache_dir, "uniprot_to_embl_cds.json")
    if os.path.exists(cache):
        with open(cache, "r", encoding="utf-8") as fh:
            cached = json.load(fh)
        if set(accessions).issubset(cached.keys()):
            return {a: cached[a] for a in accessions}

    run_url = f"{UNIPROT_REST}/idmapping/run"
    data = urllib.parse.urlencode({
        "from": "UniProtKB_AC-ID",
        "to": "EMBL-GenBank-DDBJ_CDS",
        "ids": ",".join(accessions),
    }).encode()
    req = urllib.request.Request(
        run_url, data=data,
        headers={"User-Agent": USER_AGENT,
                 "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            job = json.loads(resp.read().decode())
    except (urllib.error.URLError, ValueError) as exc:
        raise FetchError(f"ID-mapping submit failed: {exc}") from exc

    job_id = job.get("jobId")
    if not job_id:
        raise FetchError("ID-mapping returned no jobId")

    status_url = f"{UNIPROT_REST}/idmapping/status/{job_id}"
    for _ in range(60):
        st = _http_get_json(status_url)
        if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
            break
        time.sleep(3)

    results_url = f"{UNIPROT_REST}/idmapping/results/{job_id}?size=500"
    mapping: dict[str, list[str]] = {a: [] for a in accessions}
    next_url: str | None = results_url
    while next_url:
        payload, next_url = _http_get_json_with_link(next_url)
        for row in payload.get("results", []) or []:
            frm = row.get("from")
            to = row.get("to")
            if frm in mapping and to:
                mapping[frm].append(to)

    with open(cache, "w", encoding="utf-8") as fh:
        json.dump(mapping, fh)
    return mapping


def map_uniprot_to_refseq(accessions: list[str], cache_dir: str
                          ) -> dict[str, list[str]]:
    """
    Batch-map UniProt accessions -> RefSeq protein ids via the UniProt
    idmapping service (UniProtKB_AC-ID -> RefSeq_Protein).

    Used by the neighborhood feature as a fallback for accessions whose
    cached UniProt JSON has no inline RefSeq cross-reference. Many TrEMBL
    entries have no inline xref but DO have a RefSeq mapping known to
    UniProt — the idmapping service exposes those.

    Returns {accession: [refseq_id, ...]}. Accessions with no mapping get
    an empty list. The full result is cached on disk so re-runs are free.

    Defensive: if the idmapping service fails (network error, no jobId,
    timeout, malformed result), this raises FetchError. Callers should
    catch and fall through to a per-accession path; idmapping is a
    speedup-and-recovery tool, not a critical path.
    """
    cache = _cache_path(cache_dir, "uniprot_to_refseq_idmapping.json")
    if os.path.exists(cache):
        with open(cache, "r", encoding="utf-8") as fh:
            cached = json.load(fh)
        if set(accessions).issubset(cached.keys()):
            return {a: cached[a] for a in accessions}

    run_url = f"{UNIPROT_REST}/idmapping/run"
    data = urllib.parse.urlencode({
        "from": "UniProtKB_AC-ID",
        "to": "RefSeq_Protein",
        "ids": ",".join(accessions),
    }).encode()
    req = urllib.request.Request(
        run_url, data=data,
        headers={"User-Agent": USER_AGENT,
                 "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            job = json.loads(resp.read().decode())
    except (urllib.error.URLError, ValueError) as exc:
        raise FetchError(f"ID-mapping submit (RefSeq) failed: {exc}") from exc

    job_id = job.get("jobId")
    if not job_id:
        raise FetchError("ID-mapping (RefSeq) returned no jobId")

    status_url = f"{UNIPROT_REST}/idmapping/status/{job_id}"
    for _ in range(120):  # up to ~6 minutes; bigger batches take longer
        try:
            st = _http_get_json(status_url)
        except FetchError:
            time.sleep(3)
            continue
        if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
            break
        time.sleep(3)

    results_url = f"{UNIPROT_REST}/idmapping/results/{job_id}?size=500"
    mapping: dict[str, list[str]] = {a: [] for a in accessions}
    next_url: str | None = results_url
    pages = 0
    while next_url and pages < 2000:  # safety cap on pagination
        payload, next_url = _http_get_json_with_link(next_url)
        for row in payload.get("results", []) or []:
            frm = row.get("from")
            to = row.get("to")
            if frm in mapping and to:
                mapping[frm].append(to)
        pages += 1

    # Merge with any prior cache so re-runs don't lose previously-mapped IDs.
    if os.path.exists(cache):
        with open(cache, "r", encoding="utf-8") as fh:
            prior = json.load(fh)
        prior.update(mapping)
        mapping = prior

    with open(cache, "w", encoding="utf-8") as fh:
        json.dump(mapping, fh)
    return {a: mapping.get(a, []) for a in accessions}
