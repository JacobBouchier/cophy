#!/usr/bin/env python3
"""
Resolve UniProt accessions -> GenePos (assembly, contig, ordinal) via NCBI.

EBI's coordinates service is calculated from GRC/Ensembl assemblies and
does not cover bacterial genomes, which is what most CoPhy runs care
about. NCBI EFetch on the protein and nuccore databases does — every
bacterial RefSeq protein record carries a /coded_by= qualifier giving the
contig and start coordinate, and the contig's feature table gives every
CDS's coordinates so we can convert to ordinals.

Strategy per UniProt accession:
  1. UniProt accession -> RefSeq protein id
     a. Look in cached UniProt JSON for a 'RefSeq' xref (free).
     b. If absent, call NCBI esearch on the protein DB for the accession.
  2. EFetch the RefSeq protein in GenBank format -> parse /coded_by=...
     into (contig, start, strand).
  3. EFetch each unique contig's feature table (ft format) once -> sorted
     [(start, protein_id)] for the whole contig. Cached.
  4. Convert (contig, start) -> ordinal using that list.

Network calls go through one rate-limited helper. Parsing is pure.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch import _cache_path, fetch_entry, FetchError  # noqa: E402
from neighborhood import GenePos  # noqa: E402


ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "CoPhy/1.0"


# --------------------------------------------------------------------------
# Rate-limited network layer + run-wide call counter
# --------------------------------------------------------------------------

class _RunBudget:
    """Tracks total NCBI calls in this run; enforces the optional cap."""
    def __init__(self, limit: int):
        self.limit = limit  # 0 = unlimited
        self.used = 0
        self.over = False  # latches once we exceed; we keep warning quietly

    def consume(self) -> bool:
        """Return True if a call is allowed; False if the cap was hit."""
        if self.limit and self.used >= self.limit:
            if not self.over:
                sys.stderr.write(
                    f"WARN: NCBI lookup cap reached "
                    f"({self.limit}); further accessions will be blank\n")
                self.over = True
            return False
        self.used += 1
        return True


_BUDGET = _RunBudget(limit=0)
_API_KEY = ""
_EMAIL = ""
_LAST_CALL = [0.0]   # mutable to avoid global keyword

# Warn-event counters: per-call sites tally into here instead of writing to
# stderr. resolve_positions_summary() emits one consolidated block at the
# end of a run. This trims tens of thousands of stderr writes that were
# eating real I/O time on the cluster.
_WARN_COUNTS: dict[str, int] = {}


def _warn(kind: str, detail: str = "") -> None:
    """Tally a warning. detail is kept only for the first few of each kind."""
    n = _WARN_COUNTS.get(kind, 0)
    _WARN_COUNTS[kind] = n + 1
    # Keep up to 3 example details per kind for the summary.
    if detail and n < 3:
        examples = _WARN_COUNTS.setdefault(f"_examples_{kind}", [])
        examples.append(detail)  # type: ignore[union-attr]


def emit_warn_summary() -> None:
    """Flush the WARN tally as a small summary block to stderr."""
    keys = [k for k in _WARN_COUNTS if not k.startswith("_examples_")]
    if not keys:
        return
    sys.stderr.write("\n=== neighborhood WARN summary ===\n")
    for kind in sorted(keys):
        n = _WARN_COUNTS[kind]
        sys.stderr.write(f"  {kind}: {n}\n")
        ex = _WARN_COUNTS.get(f"_examples_{kind}")
        if ex:
            sys.stderr.write(f"    examples: {', '.join(ex)}\n")
    sys.stderr.write("=================================\n")


def reset_warn_counts() -> None:
    """For tests."""
    _WARN_COUNTS.clear()


def configure(api_key: str, email: str, max_lookups: int) -> None:
    """Call once before resolving anything."""
    global _BUDGET, _API_KEY, _EMAIL
    _BUDGET = _RunBudget(limit=max_lookups)
    _API_KEY = api_key
    _EMAIL = email
    reset_warn_counts()


def _entrez_get(path: str, params: dict[str, str]) -> str:
    """Rate-limited GET against NCBI Entrez. Returns text body."""
    if not _BUDGET.consume():
        raise FetchError("NCBI run budget exhausted")
    # Rate limit: 10/s with API key, 3/s without. Sleep to stay under.
    min_interval = 0.11 if _API_KEY else 0.34
    elapsed = time.time() - _LAST_CALL[0]
    if elapsed < min_interval:
        time.sleep(min_interval - elapsed)

    full_params = dict(params)
    full_params["email"] = _EMAIL
    if _API_KEY:
        full_params["api_key"] = _API_KEY
    url = f"{ENTREZ_BASE}/{path}?" + urllib.parse.urlencode(full_params)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    # Retry on transient 429/5xx.
    last: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                _LAST_CALL[0] = time.time()
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504):
                last = exc
                time.sleep(1.5 * (2 ** attempt))
                continue
            raise FetchError(f"HTTP {exc.code} for {url}") from exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            time.sleep(1.5 * (2 ** attempt))
    raise FetchError(f"NCBI {path} failed after retries: {last}")


# --------------------------------------------------------------------------
def _normalize_accession(raw: str) -> str:
    """
    Normalize an NCBI accession for comparison.

    NCBI's feature table puts database-prefixed ids like 'emb|CCO08829.1|',
    'gb|NP_414542.1|', 'ref|WP_035850696.1|' in the protein_id qualifier,
    while EFetch returns the bare accession 'CCO08829.1'. We strip the
    db|...| wrapper. We also strip the trailing version suffix (e.g. '.1')
    defensively, since some sources omit or change it.

    Returns the bare accession with no version, lowercased — both inputs
    pass through the same function, so any consistent rule works.
    """
    s = raw.strip()
    # Strip db|...| wrapper: take the inner middle piece if there are pipes.
    if "|" in s:
        parts = [p for p in s.split("|") if p]
        if parts:
            # 'emb|CCO08829.1|' -> parts ['emb','CCO08829.1']; take the
            # rightmost non-empty piece, which is the accession.
            s = parts[-1]
    # Strip trailing version suffix '.<digits>'.
    if "." in s:
        head, _, tail = s.rpartition(".")
        if tail.isdigit():
            s = head
    return s.lower()


# Pure parsers (testable without network)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class RawCoords:
    contig: str
    start: int       # 1-based, on the FORWARD strand of the contig
    strand: str      # '+' or '-'


_CODED_BY_RE = re.compile(
    r'/coded_by="(?P<body>[^"]+)"',
    re.IGNORECASE,
)


def parse_coded_by(text: str) -> list[RawCoords]:
    """
    Pull every /coded_by= qualifier out of a GenBank-format text dump and
    parse it into RawCoords entries.

    /coded_by examples we must handle:
        /coded_by="NC_000913.3:4174554..4174916"
        /coded_by="complement(NC_000913.3:4174554..4174916)"
        /coded_by="join(NC_000913.3:1..50,NC_000913.3:100..200)"
        /coded_by="complement(join(...))"

    We return one RawCoords per /coded_by qualifier, taking the *lowest*
    start coordinate across joined segments and reporting the *outer*
    complement as the strand.

    Multi-line /coded_by qualifiers (which span continuation lines in
    GenBank text) are unwrapped before this function sees them.
    """
    out: list[RawCoords] = []
    for m in _CODED_BY_RE.finditer(text):
        body = m.group("body")
        out_one = _parse_coded_by_body(body)
        if out_one is not None:
            out.append(out_one)
    return out


def _parse_coded_by_body(body: str) -> RawCoords | None:
    """Inner parser for the string inside coded_by="...". Handles
    complement() and join() wrappers and picks the minimum start."""
    body = body.strip()
    strand = "+"
    if body.lower().startswith("complement(") and body.endswith(")"):
        strand = "-"
        body = body[len("complement("):-1].strip()
    if body.lower().startswith("join(") and body.endswith(")"):
        body = body[len("join("):-1].strip()

    # Now body is a comma-separated list of "<contig>:<start>..<end>" segments.
    contig: str | None = None
    min_start: int | None = None
    for seg in body.split(","):
        seg = seg.strip()
        if not seg:
            continue
        # Expected: "<contig>:<start>..<end>" (with optional > or < for fuzzy
        # endpoints, which we strip).
        if ":" not in seg:
            continue
        c, _, rng = seg.partition(":")
        c = c.strip()
        # Range may itself be wrapped in complement(...) per-segment.
        rng = rng.strip()
        if rng.lower().startswith("complement(") and rng.endswith(")"):
            strand = "-"
            rng = rng[len("complement("):-1].strip()
        rng = rng.replace("<", "").replace(">", "")
        if ".." not in rng:
            continue
        s_str, _, _e_str = rng.partition("..")
        try:
            s = int(s_str)
        except ValueError:
            continue
        if contig is None:
            contig = c
        if min_start is None or s < min_start:
            min_start = s
    if contig is None or min_start is None:
        return None
    return RawCoords(contig=contig, start=min_start, strand=strand)


def parse_genbank_for_coords(text: str) -> list[RawCoords]:
    """Given a GenBank record (protein DB), pull /coded_by qualifiers out.
    GenBank wraps long qualifiers across multiple indented lines; we
    unwrap them first."""
    unwrapped: list[str] = []
    for line in text.splitlines():
        if line.startswith("                     ") and unwrapped:
            unwrapped[-1] = unwrapped[-1] + line.strip()
        else:
            unwrapped.append(line)
    joined = "\n".join(unwrapped)
    return parse_coded_by(joined)


def parse_refseq_xref_from_entry(obj: dict) -> list[str]:
    """Return RefSeq protein ids from a cached UniProt JSON entry."""
    out: list[str] = []
    for xref in obj.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "RefSeq":
            continue
        rsid = xref.get("id", "")
        if rsid and rsid != "-":
            out.append(rsid)
    return out


def parse_locus_tags_from_entry(obj: dict) -> list[str]:
    """
    Return locus_tag / ORF name strings from a UniProt JSON entry.

    UniProt's JSON puts these under `genes[*].orfNames[*].value` and
    `genes[*].orderedLocusNames[*].value`. Locus_tags are the names
    NCBI uses to index CDS features in nucleotide records (e.g.
    'FDA94_30825'), so when a TrEMBL accession has no inline RefSeq
    cross-reference, searching NCBI for the locus_tag is often the
    fastest way to find the corresponding protein record.

    Returns a deduped list, preserving first-seen order.
    """
    seen: list[str] = []
    for gene in obj.get("genes", []) or []:
        for tag in (gene.get("orfNames") or []) + (gene.get("orderedLocusNames") or []):
            val = tag.get("value", "").strip()
            if val and val not in seen:
                seen.append(val)
    return seen


def parse_ipg_tsv(text: str) -> list[RawCoords]:
    """
    Parse NCBI's IPG (Identical Protein Groups) tab-separated response.

    Format (header row + one row per (protein, genome) instance):
        Id  Source  Nucleotide Accession  Start  Stop  Strand  Protein  ...
        33217001  RefSeq  NZ_AHLT01000014.1  13020  13304  -  WP_000134546.1  ...

    WP_ accessions don't have /coded_by qualifiers in their GenBank
    records (they're deliberately non-redundant summaries). IPG is NCBI's
    workaround: it lists every genome where that protein is annotated.
    We return one RawCoords per row; resolve_positions will use the
    closest-pair logic to find a usable (contig, start) combination.

    Returns empty list on empty input or unrecognised headers.
    """
    out: list[RawCoords] = []
    lines = text.splitlines()
    if not lines:
        return out
    header = [c.strip() for c in lines[0].split("\t")]
    # Be defensive about column naming variations across IPG versions.
    def find_col(*names: str) -> int:
        lowered = [h.lower() for h in header]
        for n in names:
            for i, h in enumerate(lowered):
                if h == n.lower():
                    return i
        return -1

    i_nucl = find_col("Nucleotide Accession", "Nucleotide_Accession",
                       "nucleotide accession")
    i_start = find_col("Start")
    i_strand = find_col("Strand")
    if i_nucl < 0 or i_start < 0:
        return out

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= max(i_nucl, i_start):
            continue
        contig = parts[i_nucl].strip()
        if not contig or contig == "-":
            continue
        try:
            start = int(parts[i_start].strip())
        except (ValueError, IndexError):
            continue
        strand = "+"
        if 0 <= i_strand < len(parts):
            s = parts[i_strand].strip()
            if s == "-":
                strand = "-"
        out.append(RawCoords(contig=contig, start=start, strand=strand))
    return out


# --------------------------------------------------------------------------
# Feature table parser (for contig gene order)
# --------------------------------------------------------------------------

# NCBI feature table format ('rettype=ft') looks like:
#
#   >Feature gnl|NCBI|NC_000913.3
#   190	255	gene
#               gene  thrL
#               locus_tag  b0001
#   190	255	CDS
#               product  thr operon leader peptide
#               protein_id  NP_414542.1
#   complement(337..2799)	CDS
#               protein_id  NP_414543.1
#
# Feature-line columns: start TAB end TAB feature_name; complement() may
# wrap the coordinate column instead of using a sign.

def parse_feature_table_for_cds(text: str) -> list[tuple[int, str]]:
    """
    Return sorted [(start, protein_id), ...] for every CDS on the contig
    whose feature table we just fetched. Start is the lower numeric
    coordinate regardless of strand.
    """
    genes: list[tuple[int, str]] = []
    current_kind: str | None = None
    current_start: int | None = None
    current_protein: str | None = None

    def _flush() -> None:
        nonlocal current_kind, current_start, current_protein
        if (current_kind == "CDS" and current_start is not None
                and current_protein):
            genes.append((current_start, current_protein))
        current_kind = None
        current_start = None
        current_protein = None

    for raw in text.splitlines():
        if not raw or raw.startswith(">Feature"):
            _flush()
            continue
        # Feature line: "<start>\t<end>\t<kind>" — start/end may carry
        # fuzzy markers (<, >) but NCBI's rettype=ft doesn't use
        # complement(); reverse-strand features have start > end.
        if raw[0].isdigit() or raw[0] in "<>":
            _flush()
            parts = raw.split("\t")
            if len(parts) < 3:
                continue
            kind = parts[2].strip()
            if kind != "CDS":
                continue
            start_field = parts[0].strip().replace("<", "").replace(">", "")
            end_field = parts[1].strip().replace("<", "").replace(">", "")
            try:
                if "^" in start_field:
                    start_field = start_field.split("^", 1)[0]
                s = int(start_field)
                e = int(end_field) if end_field else s
            except ValueError:
                continue
            # On reverse strand, NCBI writes end < start. Normalise to the
            # lower coordinate so the gene-order sort is consistent.
            current_kind = "CDS"
            current_start = min(s, e)
        elif raw.startswith("\t\t\t") or raw.startswith("   "):
            # Qualifier line: "\t\t\tkey\tvalue" (or 12-space prefix).
            stripped = raw.strip()
            # Form is "<key>\t<value>".
            if "\t" in stripped:
                k, _, v = stripped.partition("\t")
            else:
                k, _, v = stripped.partition(" ")
            k = k.strip().lower()
            v = v.strip()
            if k == "protein_id":
                # Strip db|...| wrappers and version suffix for matching.
                current_protein = _normalize_accession(v)
    _flush()
    genes.sort(key=lambda x: x[0])
    return genes


# --------------------------------------------------------------------------
# Cached high-level resolvers
# --------------------------------------------------------------------------

def _read_uniprot_cache(accession: str, cache_dir: str) -> dict | None:
    """Read the cached UniProt JSON if present (written by fetch_entry)."""
    path = _cache_path(cache_dir, f"uniprot_{accession}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def prefetch_idmapping(accessions: list[str], cache_dir: str) -> int:
    """
    Batch pre-pass: use UniProt idmapping to resolve UniProt -> RefSeq for
    accessions that don't already have a RefSeq id from inline cross-refs.

    This populates per-accession ncbi_refseq_for_<acc>.json files so the
    regular uniprot_to_refseq() function picks the mappings up
    transparently. Accessions that already have a non-empty cached
    mapping (from xref or a prior esearch) are skipped.

    Defensive: any failure here is non-fatal — we just log it via _warn
    and return 0. The per-accession path will run as before for whatever
    didn't get pre-resolved.

    Returns: count of accessions newly resolved.
    """
    from fetch import map_uniprot_to_refseq  # local import: same module dir

    # Filter: which accessions don't yet have a mapping?
    needs_lookup: list[str] = []
    for acc in accessions:
        cache_file = _cache_path(cache_dir, f"ncbi_refseq_for_{acc}.json")
        if os.path.exists(cache_file):
            try:
                with open(cache_file) as fh:
                    if json.load(fh):
                        continue  # already has a non-empty mapping
            except (json.JSONDecodeError, OSError):
                pass  # treat unreadable cache as missing
        # Also try the inline xref before adding to the batch list — free.
        obj = _read_uniprot_cache(acc, cache_dir)
        if obj is not None and parse_refseq_xref_from_entry(obj):
            continue
        needs_lookup.append(acc)

    if not needs_lookup:
        sys.stderr.write(
            "neighborhood: idmapping pre-pass not needed (all accessions "
            "already mapped or have inline RefSeq xrefs)\n")
        return 0

    sys.stderr.write(
        f"neighborhood: idmapping pre-pass for {len(needs_lookup)} accessions "
        f"(batch UniProt -> RefSeq)\n")

    try:
        mapping = map_uniprot_to_refseq(needs_lookup, cache_dir)
    except Exception as exc:  # noqa: BLE001 — must not break the rule
        _warn("idmapping_failed", str(exc)[:200])
        sys.stderr.write(
            f"neighborhood: idmapping batch failed ({exc}); "
            f"continuing with per-accession path\n")
        return 0

    # Write per-accession cache files so uniprot_to_refseq() finds them.
    n_resolved = 0
    for acc in needs_lookup:
        refseqs = mapping.get(acc, [])
        if refseqs:
            n_resolved += 1
        # Always write (even empty) so we don't repeatedly retry esearch
        # for accessions idmapping confirmed have no RefSeq.
        cache_file = _cache_path(cache_dir, f"ncbi_refseq_for_{acc}.json")
        with open(cache_file, "w") as fh:
            json.dump(refseqs, fh)

    sys.stderr.write(
        f"neighborhood: idmapping recovered {n_resolved}/{len(needs_lookup)} "
        f"accessions ({100*n_resolved//max(len(needs_lookup),1)}%)\n")
    return n_resolved


def _esearch_protein_to_accessions(term: str) -> list[str]:
    """
    Run a single NCBI esearch on the protein DB and return up to 5
    accessions (with version) for the matching hits.

    Returns [] on no hits or on any failure (failures are tallied by the
    caller via _warn, not here). Raises FetchError for the caller to
    catch and tally.
    """
    text = _entrez_get("esearch.fcgi", {
        "db": "protein", "term": term, "retmode": "json"})
    payload = json.loads(text)
    ids = payload.get("esearchresult", {}).get("idlist", []) or []
    if not ids:
        return []
    summary = _entrez_get("esummary.fcgi", {
        "db": "protein", "id": ",".join(ids[:5]), "retmode": "json"})
    spayload = json.loads(summary)
    result = spayload.get("result", {})
    out: list[str] = []
    for uid in ids[:5]:
        item = result.get(uid, {})
        acc = item.get("accessionversion") or item.get("caption")
        if acc:
            out.append(acc)
    return out


def uniprot_to_refseq(accession: str, cache_dir: str) -> list[str]:
    """
    Map UniProt accession -> RefSeq / GenBank protein id(s).

    Strategy, in order (each step skipped if the prior gave a hit):
      1. UniProt cross-references: cached, no network.
      2. Locus_tag esearch: take the UniProt entry's orfNames /
         orderedLocusNames (e.g. 'FDA94_30825') and search the NCBI
         protein DB by Gene Name. This recovers most TrEMBL A0A...
         accessions whose entries do have a locus_tag from the parent
         assembly submission, even when UniProt has no inline RefSeq
         xref. Searching by locus_tag is more specific than searching
         by the bare UniProt accession.
      3. Bare-accession esearch: last-ditch fallback for entries that
         have neither RefSeq xref nor a locus_tag.

    Empty result is cached so we don't re-search forever. A non-empty
    result is also cached.
    """
    cache = _cache_path(cache_dir, f"ncbi_refseq_for_{accession}.json")
    if os.path.exists(cache):
        with open(cache) as fh:
            return json.load(fh)

    refs: list[str] = []
    obj = _read_uniprot_cache(accession, cache_dir)

    # 1. UniProt cross-references (free).
    if obj is not None:
        refs = parse_refseq_xref_from_entry(obj)

    # 2. Locus_tag-first esearch.
    if not refs and obj is not None:
        locus_tags = parse_locus_tags_from_entry(obj)
        for tag in locus_tags:
            try:
                # Use [Gene Name] field — NCBI indexes locus_tags there
                # for the protein DB. Quote in case the tag contains
                # characters NCBI treats specially.
                hits = _esearch_protein_to_accessions(f'"{tag}"[Gene Name]')
            except (FetchError, ValueError) as exc:
                _warn("esearch_locus_tag_failed", f"{accession}/{tag}: {exc}")
                hits = []
            if hits:
                refs.extend(hits)
                break  # one good locus_tag is enough; don't keep searching

    # 3. Bare-accession esearch (last resort).
    if not refs:
        try:
            refs = _esearch_protein_to_accessions(accession)
        except (FetchError, ValueError) as exc:
            _warn("esearch_failed", f"{accession}: {exc}")

    with open(cache, "w") as fh:
        json.dump(refs, fh)
    return refs


def fetch_ipg_coords(refseq_id: str, cache_dir: str) -> list[RawCoords]:
    """EFetch the IPG entry for a protein -> list of (contig, start) per
    genome where it's annotated. Used when the protein's GenBank record
    lacks /coded_by (typical for WP_ non-redundant summaries)."""
    cache = _cache_path(cache_dir, f"ncbi_ipg_{refseq_id}.tsv")
    if os.path.exists(cache):
        with open(cache) as fh:
            text = fh.read()
    else:
        try:
            text = _entrez_get("efetch.fcgi", {
                "db": "ipg", "id": refseq_id,
                "rettype": "ipg", "retmode": "text"})
        except FetchError as exc:
            _warn("efetch_ipg_failed", f"{refseq_id}: {exc}")
            text = ""
        with open(cache, "w") as fh:
            fh.write(text)
    return parse_ipg_tsv(text) if text else []


def fetch_protein_coords(refseq_id: str, cache_dir: str) -> list[RawCoords]:
    """EFetch a RefSeq protein in GenBank -> parsed RawCoords list.

    Falls back to the IPG database when the GenBank record lacks any
    /coded_by qualifier (typical for WP_ non-redundant accessions, which
    are summaries with no single genome location)."""
    cache = _cache_path(cache_dir, f"ncbi_protein_{refseq_id}.gb")
    if os.path.exists(cache):
        with open(cache) as fh:
            text = fh.read()
    else:
        try:
            text = _entrez_get("efetch.fcgi", {
                "db": "protein", "id": refseq_id,
                "rettype": "gp", "retmode": "text"})
        except FetchError as exc:
            _warn("efetch_protein_failed", f"{refseq_id}: {exc}")
            text = ""
        with open(cache, "w") as fh:
            fh.write(text)

    coords = parse_genbank_for_coords(text) if text else []
    if not coords:
        # WP_ records and other non-redundant entries have no /coded_by;
        # IPG lists all the genomes that contain this protein.
        coords = fetch_ipg_coords(refseq_id, cache_dir)
    return coords


def fetch_contig_gene_order(contig: str, cache_dir: str
                            ) -> list[tuple[int, str]]:
    """EFetch a nuccore feature table -> sorted [(start, protein_id)]."""
    safe = contig.replace("/", "_")
    cache = _cache_path(cache_dir, f"ncbi_contig_{safe}.ft")
    if os.path.exists(cache):
        with open(cache) as fh:
            text = fh.read()
    else:
        try:
            text = _entrez_get("efetch.fcgi", {
                "db": "nuccore", "id": contig,
                "rettype": "ft", "retmode": "text"})
        except FetchError as exc:
            _warn("efetch_contig_failed", f"{contig}: {exc}")
            text = ""
        with open(cache, "w") as fh:
            fh.write(text)
    return parse_feature_table_for_cds(text) if text else []


def resolve_positions(accession: str, cache_dir: str) -> list[GenePos]:
    """
    Public entry: UniProt accession -> list of GenePos (one per
    assembly/contig where it's annotated).

    For each RefSeq protein id we get from UniProt, we fetch either the
    GenBank record (whose /coded_by tells us contig + start) OR the IPG
    table (which tells us the same thing across multiple genomes). We
    convert each result directly into a GenePos with nucleotide-level
    start coordinate. No feature-table lookup is needed because distance
    is now measured in nucleotides between start coordinates rather than
    gene ordinals.

    All NCBI calls are cached on disk.
    """
    refseq_ids = uniprot_to_refseq(accession, cache_dir)
    if not refseq_ids:
        _warn("no_refseq_protein", accession)
        return []

    out: list[GenePos] = []
    for rsid in refseq_ids:
        coords_list = fetch_protein_coords(rsid, cache_dir)
        if not coords_list:
            _warn("no_coords_any_source", f"{rsid} (from {accession})")
            continue
        for c in coords_list:
            # We use contig accession as the assembly key. Bacterial
            # genomes usually have one chromosome and the contig accession
            # uniquely identifies it. For multi-replicon assemblies,
            # different contigs correctly fail to share (assembly, contig).
            out.append(GenePos(assembly=c.contig, contig=c.contig,
                               start=c.start, strand=c.strand))

    # Dedup identical positions.
    seen: dict[tuple, GenePos] = {}
    for p in out:
        seen[(p.assembly, p.contig, p.start, p.strand)] = p
    return list(seen.values())
