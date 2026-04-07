"""
NCBI / PubMed query helpers built on top of ``Bio.Entrez``.

All network-bound work lives here so the Streamlit layer stays thin.
"""

from __future__ import annotations

import re
import time
import textwrap
from dataclasses import dataclass, field
from typing import Optional

from Bio import Entrez, SeqIO


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class CDSResult:
    """Holds everything extracted from a GenBank CDS feature."""

    accession: str = ""
    description: str = ""
    cds_dna: str = ""           # raw origin sub-sequence for CDS
    cds_start: int = 0          # 1-based inclusive
    cds_end: int = 0            # 1-based inclusive
    translation: str = ""
    full_sequence: str = ""     # entire mRNA record origin


@dataclass
class PrimerHit:
    """One primer found in a PubMed article."""

    direction: str = ""         # "Forward" or "Reverse"
    sequence: str = ""
    position: Optional[int] = None
    strand: str = ""


@dataclass
class PubMedPrimerResult:
    """Aggregates primer data from a single PubMed article."""

    pmid: str = ""
    pmcid: str = ""
    title: str = ""
    citation: str = ""
    primers: list[PrimerHit] = field(default_factory=list)
    raw_text: str = ""          # the text that was searched


# ── Entrez helpers ────────────────────────────────────────────────────


def search_nucleotide(gene: str, organism: str = "Homo sapiens", retmax: int = 20) -> list[dict]:
    """
    Search NCBI Nucleotide for *gene* in *organism*.

    Returns a list of dicts with keys ``id``, ``title``.
    """
    if organism.lower() in ("human", "homo sapiens"):
        org_query = '("Homo sapiens"[Organism] OR "human gene")'
    elif organism.lower() in ("mouse", "mus musculus"):
        org_query = '("Mus musculus"[Organism] OR "mouse gene")'
    else:
        org_query = f'("{organism}"[Organism])'

    query = f"{gene} AND {org_query}"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        return []

    # Fetch summaries so we can display titles
    handle = Entrez.esummary(db="nucleotide", id=",".join(ids), retmax=retmax)
    summaries = Entrez.read(handle)
    handle.close()

    results = []
    for s in summaries:
        results.append({"id": str(s["Id"]), "title": s.get("Title", ""), "accession": s.get("AccessionVersion", s.get("Caption", ""))})
    return results


def _pick_best_record(results: list[dict], gene: str) -> Optional[dict]:
    """
    Heuristically pick the record most likely to be
    "transcript variant 1, mRNA" for the given gene.
    """
    gene_upper = gene.upper()
    scored: list[tuple[int, dict]] = []
    for r in results:
        title = r["title"].upper()
        score = 0
        # Must mention the gene
        if gene_upper in title:
            score += 10
        else:
            continue
        # Prefer mRNA records
        if "MRNA" in title:
            score += 5
        # Prefer transcript variant 1
        if "TRANSCRIPT VARIANT 1" in title:
            score += 8
        elif "TRANSCRIPT VARIANT" in title:
            score += 2
        # Prefer Homo sapiens
        if "HOMO SAPIENS" in title:
            score += 3
        # Penalise predicted / model records
        if "PREDICTED" in title:
            score -= 6
        scored.append((score, r))
    if not scored:
        return results[0] if results else None
    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1]


def fetch_genbank_record(accession_or_id: str) -> SeqIO.SeqRecord:
    """Download and parse a single GenBank flat-file record."""
    handle = Entrez.efetch(db="nucleotide", id=accession_or_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def extract_cds(record: SeqIO.SeqRecord) -> CDSResult:
    """
    Pull the first CDS feature from a SeqRecord and build a ``CDSResult``.

    Raises ``ValueError`` if the record has no CDS annotation.
    """
    cds_feature = None
    for feat in record.features:
        if feat.type == "CDS":
            cds_feature = feat
            break
    if cds_feature is None:
        raise ValueError("No CDS feature found in this GenBank record.")

    start = int(cds_feature.location.start)   # 0-based
    end = int(cds_feature.location.end)        # 0-based exclusive
    cds_dna = str(record.seq[start:end])

    translation = cds_feature.qualifiers.get("translation", [""])[0]

    accession = record.id
    description = record.description

    return CDSResult(
        accession=accession,
        description=description,
        cds_dna=cds_dna,
        cds_start=start + 1,   # 1-based
        cds_end=end,           # 1-based inclusive
        translation=translation,
        full_sequence=str(record.seq),
    )


def lookup_gene(gene: str, organism: str = "Homo sapiens") -> tuple[list[dict], Optional[CDSResult]]:
    """
    High-level convenience: search → pick best → fetch → extract CDS.

    Returns ``(summaries, cds_result | None)``.
    """
    summaries = search_nucleotide(gene, organism)
    if not summaries:
        return summaries, None

    best = _pick_best_record(summaries, gene)
    if best is None:
        return summaries, None

    record = fetch_genbank_record(best["id"])
    cds = extract_cds(record)
    return summaries, cds


# ── PubMed primer search ─────────────────────────────────────────────


def search_pubmed_for_primers(gene: str, organism: str = "human", retmax: int = 15) -> list[dict]:
    """
    Search PubMed for articles likely to contain PCR primers for *gene*.

    Returns basic metadata (PMID, title, source, authors, date).
    """
    query = f'("{gene}" AND (primer OR "PCR" OR "RT-PCR") AND ("{organism}"))'
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="xml")
    records = Entrez.read(handle)
    handle.close()

    articles: list[dict] = []
    for article in records.get("PubmedArticle", []):
        medline = article.get("MedlineCitation", {})
        pmid = str(medline.get("PMID", ""))
        art = medline.get("Article", {})
        title = art.get("ArticleTitle", "")

        # Build a rough citation
        journal = art.get("Journal", {})
        journal_title = journal.get("ISOAbbreviation", journal.get("Title", ""))
        pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
        year = pub_date.get("Year", "")
        month = pub_date.get("Month", "")
        volume = journal.get("JournalIssue", {}).get("Volume", "")
        issue = journal.get("JournalIssue", {}).get("Issue", "")
        pages = art.get("Pagination", {}).get("MedlinePgn", "")

        # Authors
        author_list = art.get("AuthorList", [])
        authors_str = ""
        if author_list:
            names = []
            for a in author_list:
                ln = a.get("LastName", "")
                ini = a.get("Initials", "")
                if ln:
                    names.append(f"{ln} {ini}".strip())
            authors_str = ", ".join(names)

        # Abstract
        abstract_parts = art.get("Abstract", {}).get("AbstractText", [])
        abstract = " ".join(str(p) for p in abstract_parts)

        # DOI / PMCID
        doi = ""
        pmcid = ""
        id_list = article.get("PubmedData", {}).get("ArticleIdList", [])
        for aid in id_list:
            attrs = aid.attributes if hasattr(aid, "attributes") else {}
            if attrs.get("IdType") == "doi":
                doi = str(aid)
            elif attrs.get("IdType") == "pmc":
                pmcid = str(aid)

        citation_parts = []
        if authors_str:
            citation_parts.append(authors_str + ".")
        if title:
            citation_parts.append(title)
        journal_part = journal_title
        if year:
            journal_part += f" {year}"
        if month:
            journal_part += f" {month}"
        if volume:
            journal_part += f";{volume}"
        if issue:
            journal_part += f"({issue})"
        if pages:
            journal_part += f":{pages}"
        journal_part += "."
        citation_parts.append(journal_part)
        if doi:
            citation_parts.append(f"doi: {doi};")
        if pmcid:
            citation_parts.append(f"PMCID: {pmcid}")

        citation = " ".join(citation_parts)

        articles.append({
            "pmid": pmid,
            "pmcid": pmcid,
            "title": title,
            "abstract": abstract,
            "citation": citation,
            "doi": doi,
            "authors": authors_str,
        })

    return articles


def fetch_pmc_fulltext(pmcid: str) -> str:
    """
    Attempt to fetch full-text XML from PMC and return it as raw text
    (tags stripped).  Returns empty string on failure.
    """
    try:
        pmc_id = pmcid.replace("PMC", "")
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        raw = handle.read()
        handle.close()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        # Crude tag stripping – good enough for regex primer search
        text = re.sub(r"<[^>]+>", " ", raw)
        text = re.sub(r"\s+", " ", text)
        return text
    except Exception:
        return ""


def extract_primers_from_article(article: dict, gene: str) -> list[PrimerHit]:
    """
    Try to extract primer sequences from an article's abstract and,
    if a PMCID exists, from the full text.

    Returns a list of ``PrimerHit`` (direction not yet assigned).
    """
    from utils import extract_primers_from_text

    texts_to_search: list[str] = []
    if article.get("abstract"):
        texts_to_search.append(article["abstract"])

    pmcid = article.get("pmcid", "")
    if pmcid:
        ft = fetch_pmc_fulltext(pmcid)
        if ft:
            texts_to_search.append(ft)

    primers: list[PrimerHit] = []
    seen: set[str] = set()
    for txt in texts_to_search:
        for seq in extract_primers_from_text(txt):
            if seq not in seen:
                seen.add(seq)
                primers.append(PrimerHit(sequence=seq))
    return primers
