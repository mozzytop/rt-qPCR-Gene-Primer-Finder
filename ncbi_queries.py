"""
NCBI / PubMed / PMC query helpers built on ``Bio.Entrez``.

**Key upgrade**: primer search now targets PubMed Central (PMC) full-text
instead of PubMed abstracts. Primers are extracted with direction-aware
regex and programmatically verified against the CDS before being returned.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Optional

from Bio import Entrez, SeqIO

from utils import (
    extract_primers_with_direction,
    find_primer_on_either_strand,
    reverse_complement,
    verify_primer_pair,
    ExtractedPrimer,
)


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class CDSResult:
    """Holds everything extracted from a GenBank CDS feature."""
    accession: str = ""
    description: str = ""
    cds_dna: str = ""
    cds_start: int = 0
    cds_end: int = 0
    translation: str = ""
    full_sequence: str = ""


@dataclass
class VerifiedPrimerPair:
    """A primer pair that has been verified against the CDS."""
    forward: str = ""
    reverse: str = ""
    reverse_comp: str = ""
    fwd_position: Optional[int] = None
    fwd_strand: str = ""
    rev_position: Optional[int] = None
    rev_strand: str = ""
    source_pmcid: str = ""
    source_pmid: str = ""
    citation: str = ""
    title: str = ""
    doi: str = ""


# ── NCBI Nucleotide helpers ──────────────────────────────────────────


def parse_ncbi_url(url: str) -> Optional[str]:
    """Extract an accession or GI from an NCBI nuccore URL.

    Handles forms like:
      https://www.ncbi.nlm.nih.gov/nuccore/NM_005269.3
      https://www.ncbi.nlm.nih.gov/nuccore/AF003837.1
      https://www.ncbi.nlm.nih.gov/nuccore/12345678
    """
    m = re.search(r"nuccore/([A-Za-z0-9_.]+)", url)
    return m.group(1) if m else None


def search_nucleotide(gene: str, organism: str = "Homo sapiens", retmax: int = 100) -> list[dict]:
    """Search NCBI Nucleotide for *gene* in *organism*.

    Returns a list of dicts with: id, title, accession, length, link, score.
    Results are scored and sorted so the best "transcript variant 1, mRNA"
    record appears first, but ALL results are returned for user selection.
    """
    gene = gene.strip()
    if organism.lower() in ("human", "homo sapiens"):
        organism_name = "Homo sapiens"
    elif organism.lower() in ("mouse", "mus musculus"):
        organism_name = "Mus musculus"
    else:
        organism_name = organism.strip()

    org_query = f'"{organism_name}"[Organism]'

    # Query titles first so gene symbols like HES1 do not float unrelated
    # transcripts to the top. Fall back to the broader legacy query only if
    # the title-constrained search is sparse.
    query_variants = [
        f'("{gene}"[Title]) AND {org_query} AND mRNA[Title]',
        f'("{gene}"[Title]) AND {org_query}',
        f'("{gene}"[Gene Name]) AND {org_query} AND mRNA[Title]',
        f"{gene} AND {org_query}",
    ]

    ids: list[str] = []
    seen_ids: set[str] = set()
    for query in query_variants:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
        record = Entrez.read(handle)
        handle.close()
        for found_id in record.get("IdList", []):
            if found_id not in seen_ids:
                seen_ids.add(found_id)
                ids.append(found_id)
            if len(ids) >= retmax:
                break
        if len(ids) >= retmax:
            break

    if not ids:
        return []

    handle = Entrez.esummary(db="nucleotide", id=",".join(ids), retmax=retmax)
    summaries = Entrez.read(handle)
    handle.close()

    gene_upper = gene.upper()
    results = []
    for s in summaries:
        acc = s.get("AccessionVersion", s.get("Caption", ""))
        title = s.get("Title", "")
        length = int(s.get("Length", 0))
        link = f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"

        # Score for sorting
        title_upper = title.upper()
        score = 0
        if gene_upper in title_upper:
            score += 10
        if "MRNA" in title_upper:
            score += 5
        if "TRANSCRIPT VARIANT 1" in title_upper:
            score += 8
        elif "TRANSCRIPT VARIANT" in title_upper:
            score += 2
        if "HOMO SAPIENS" in title_upper or "MUS MUSCULUS" in title_upper:
            score += 3
        if "PREDICTED" in title_upper:
            score -= 6

        results.append({
            "id": str(s["Id"]),
            "title": title,
            "accession": acc,
            "length": length,
            "link": link,
            "score": score,
        })

    results.sort(key=lambda r: r["score"], reverse=True)
    return results


def _pick_best_record(results: list[dict], gene: str) -> Optional[dict]:
    """Heuristically pick the best "transcript variant 1, mRNA" record."""
    gene_upper = gene.upper()
    scored: list[tuple[int, dict]] = []
    for r in results:
        title = r["title"].upper()
        score = 0
        if gene_upper in title:
            score += 10
        else:
            continue
        if "MRNA" in title:
            score += 5
        if "TRANSCRIPT VARIANT 1" in title:
            score += 8
        elif "TRANSCRIPT VARIANT" in title:
            score += 2
        if "HOMO SAPIENS" in title:
            score += 3
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
    """Pull the first CDS feature from a SeqRecord."""
    cds_feature = None
    for feat in record.features:
        if feat.type == "CDS":
            cds_feature = feat
            break
    if cds_feature is None:
        raise ValueError("No CDS feature found in this GenBank record.")

    start = int(cds_feature.location.start)
    end = int(cds_feature.location.end)
    cds_dna = str(record.seq[start:end])
    translation = cds_feature.qualifiers.get("translation", [""])[0]

    return CDSResult(
        accession=record.id,
        description=record.description,
        cds_dna=cds_dna,
        cds_start=start + 1,
        cds_end=end,
        translation=translation,
        full_sequence=str(record.seq),
    )


def lookup_gene(gene: str, organism: str = "Homo sapiens") -> tuple[list[dict], Optional[CDSResult]]:
    """High-level: search → pick best → fetch → extract CDS."""
    summaries = search_nucleotide(gene, organism)
    if not summaries:
        return summaries, None
    best = _pick_best_record(summaries, gene)
    if best is None:
        return summaries, None
    record = fetch_genbank_record(best["id"])
    cds = extract_cds(record)
    return summaries, cds


def fetch_cds_by_accession(accession: str) -> CDSResult:
    """Fetch a GenBank record by accession and extract CDS."""
    record = fetch_genbank_record(accession)
    return extract_cds(record)


# ═════════════════════════════════════════════════════════════════════
# PMC FULL-TEXT PRIMER SEARCH (upgraded)
# ═════════════════════════════════════════════════════════════════════


def _fetch_pmc_articles_from_query(query: str, retmax: int) -> list[dict]:
    """Run a PMC query and fetch article metadata for the matching IDs."""
    handle = Entrez.esearch(db="pmc", term=query, retmax=retmax, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        return []

    articles: list[dict] = []
    for pmc_id in ids:
        try:
            article_info = _fetch_pmc_summary(pmc_id)
            if article_info:
                articles.append(article_info)
        except Exception:
            continue
        time.sleep(0.35)  # respect rate limit

    return articles


def search_pmc_for_primers(gene: str, organism: str = "human", retmax: int = 20) -> list[dict]:
    """
    Search **PubMed Central** (not regular PubMed) for open-access articles
    likely containing PCR primer sequences for *gene*.

    PMC gives us full-text access which is where primer sequences actually live
    (Materials & Methods, supplementary tables, etc).

    Returns basic metadata: PMCID, PMID, title, authors, citation.
    """
    query = (
        f'"{gene}" AND (primer OR "PCR" OR "RT-PCR" OR "qPCR" OR "real-time PCR")'
        f' AND ("{organism}") AND open access[filter]'
    )
    articles = _fetch_pmc_articles_from_query(query, retmax)
    if articles:
        return articles

    # Fallback: try without the open-access filter and slightly shorter PCR terms.
    query_fallback = (
        f'"{gene}" AND (primer OR "PCR" OR "RT-PCR" OR "qPCR")'
        f' AND ("{organism}")'
    )
    return _fetch_pmc_articles_from_query(query_fallback, retmax)


def search_pmc_for_gene_mentions(gene: str, retmax: int = 100) -> list[dict]:
    """
    Broader PMC fallback for cases where a paper contains primers in a table or
    supplement but does not rank well for PCR-specific queries.
    """
    query = f'"{gene}" AND open access[filter]'
    articles = _fetch_pmc_articles_from_query(query, retmax)
    if articles:
        return articles
    return _fetch_pmc_articles_from_query(f'"{gene}"', retmax)


def _fetch_pmc_summary(pmc_id: str) -> Optional[dict]:
    """Fetch metadata for a single PMC article using efetch XML."""
    try:
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        raw = handle.read()
        handle.close()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")

        # Extract metadata from XML with regex (lightweight, no lxml needed)
        title = _xml_extract(raw, "article-title") or "Unknown Title"
        pmid = _xml_extract(raw, "article-id", attr_filter='pub-id-type="pmid"') or ""
        pmcid = f"PMC{pmc_id}"
        doi = _xml_extract(raw, "article-id", attr_filter='pub-id-type="doi"') or ""

        # Authors
        author_matches = re.findall(
            r"<surname>([^<]+)</surname>\s*<given-names>([^<]+)</given-names>",
            raw,
        )
        authors = ", ".join(f"{sn} {gn[0]}" for sn, gn in author_matches[:6]) if author_matches else ""
        if len(author_matches) > 6:
            authors += " et al."

        # Journal info
        journal = _xml_extract(raw, "journal-title") or _xml_extract(raw, "abbrev-journal-title") or ""
        year = _xml_extract(raw, "year") or ""
        volume = _xml_extract(raw, "volume") or ""
        issue = _xml_extract(raw, "issue") or ""
        fpage = _xml_extract(raw, "fpage") or ""
        lpage = _xml_extract(raw, "lpage") or ""
        pages = f"{fpage}-{lpage}" if fpage and lpage else fpage

        # Build citation
        cite_parts = []
        if authors:
            cite_parts.append(f"{authors}.")
        if title:
            cite_parts.append(title)
        journal_part = journal
        if year:
            journal_part += f" {year}"
        if volume:
            journal_part += f";{volume}"
        if issue:
            journal_part += f"({issue})"
        if pages:
            journal_part += f":{pages}."
        else:
            journal_part += "."
        cite_parts.append(journal_part)
        if doi:
            cite_parts.append(f"doi: {doi};")
        cite_parts.append(f"PMCID: {pmcid}")

        citation = " ".join(cite_parts)

        return {
            "pmc_id": pmc_id,
            "pmcid": pmcid,
            "pmid": pmid,
            "title": title,
            "citation": citation,
            "doi": doi,
            "authors": authors,
            "raw_xml": raw,  # keep for full-text primer extraction
        }
    except Exception:
        return None


def _xml_extract(xml: str, tag: str, attr_filter: str = "") -> Optional[str]:
    """Quick-and-dirty single-value XML tag extractor."""
    if attr_filter:
        pattern = rf"<{tag}\s+{attr_filter}[^>]*>([^<]+)</{tag}>"
    else:
        pattern = rf"<{tag}[^>]*>([^<]+)</{tag}>"
    m = re.search(pattern, xml)
    return m.group(1).strip() if m else None


def fetch_pmc_fulltext(pmc_id: str, raw_xml: str = "") -> str:
    """
    Get stripped plain text from a PMC article.

    If *raw_xml* is already available (from the search step), reuses it.
    Otherwise fetches from NCBI.
    """
    if not raw_xml:
        try:
            handle = Entrez.efetch(db="pmc", id=pmc_id.replace("PMC", ""), rettype="xml")
            raw_xml = handle.read()
            handle.close()
            if isinstance(raw_xml, bytes):
                raw_xml = raw_xml.decode("utf-8", errors="replace")
        except Exception:
            return ""

    # Strip XML tags → plain text
    text = re.sub(r"<[^>]+>", " ", raw_xml)
    text = re.sub(r"\s+", " ", text)
    return text


def extract_and_verify_primers(
    article: dict,
    gene: str,
    cds_sequence: str,
) -> list[VerifiedPrimerPair]:
    """
    The core upgrade: extract primers from a PMC full-text article,
    assign directions, and programmatically verify them against the CDS.

    Steps:
      1. Get full text from the article's raw_xml.
      2. Run direction-aware regex extraction.
      3. For each candidate pair (forward + reverse), verify against CDS.
      4. Return only verified pairs.
    """
    raw_xml = article.get("raw_xml", "")
    pmc_id = article.get("pmc_id", "")
    full_text = fetch_pmc_fulltext(pmc_id, raw_xml)

    if not full_text:
        return []

    # ── Extract primers with direction labels ─────────────────────────
    candidates = extract_primers_with_direction(full_text)

    if not candidates:
        return []

    # ── Separate into forward / reverse / unknown ─────────────────────
    forwards = [c for c in candidates if c.direction == "forward"]
    reverses = [c for c in candidates if c.direction == "reverse"]
    unknowns = [c for c in candidates if c.direction == ""]

    # ── Try to match unknown primers by CDS strand mapping ────────────
    # Forward primers typically match the sense strand;
    # Reverse primers typically match the antisense strand.
    for unk in unknowns:
        pos_sense = find_primer_on_either_strand(unk.sequence, cds_sequence)
        if pos_sense[1] == "sense":
            forwards.append(unk)
        elif pos_sense[1] == "antisense":
            reverses.append(unk)

    # ── Check all possible fwd × rev combinations ────────────────────
    verified_pairs: list[VerifiedPrimerPair] = []

    for fwd_candidate in forwards:
        for rev_candidate in reverses:
            result = verify_primer_pair(fwd_candidate.sequence, rev_candidate.sequence, cds_sequence)
            if result["both_map"]:
                pair = VerifiedPrimerPair(
                    forward=result["forward_seq"],
                    reverse=result["reverse_seq"],
                    reverse_comp=result["reverse_complement"],
                    fwd_position=result["fwd_sense_pos"],
                    fwd_strand="sense",
                    rev_position=result["reverse_antisense_binding_pos"],
                    rev_strand="antisense",
                    source_pmcid=article.get("pmcid", ""),
                    source_pmid=article.get("pmid", ""),
                    citation=article.get("citation", ""),
                    title=article.get("title", ""),
                    doi=article.get("doi", ""),
                )
                verified_pairs.append(pair)

    # If we had no labelled fwd/rev but have ≥2 unknowns that map,
    # try pairing them by strand (sense = fwd, antisense = rev)
    if not verified_pairs and len(unknowns) >= 2:
        sense_hits = []
        antisense_hits = []
        for unk in unknowns:
            pos, strand = find_primer_on_either_strand(unk.sequence, cds_sequence)
            if strand == "sense":
                sense_hits.append((unk, pos))
            elif strand == "antisense":
                antisense_hits.append((unk, pos))
        for fwd_unk, fwd_pos in sense_hits:
            for rev_unk, rev_pos in antisense_hits:
                result = verify_primer_pair(fwd_unk.sequence, rev_unk.sequence, cds_sequence)
                if result["both_map"]:
                    pair = VerifiedPrimerPair(
                        forward=result["forward_seq"],
                        reverse=result["reverse_seq"],
                        reverse_comp=result["reverse_complement"],
                        fwd_position=result["fwd_sense_pos"],
                        fwd_strand="sense",
                        rev_position=result["reverse_antisense_binding_pos"],
                        rev_strand="antisense",
                        source_pmcid=article.get("pmcid", ""),
                        source_pmid=article.get("pmid", ""),
                        citation=article.get("citation", ""),
                        title=article.get("title", ""),
                        doi=article.get("doi", ""),
                    )
                    verified_pairs.append(pair)

    return verified_pairs


# ── Legacy PubMed search (fallback) ──────────────────────────────────


def search_pubmed_for_primers(gene: str, organism: str = "human", retmax: int = 15) -> list[dict]:
    """Fallback PubMed abstract search."""
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

        journal = art.get("Journal", {})
        journal_title = journal.get("ISOAbbreviation", journal.get("Title", ""))
        pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
        year = pub_date.get("Year", "")
        volume = journal.get("JournalIssue", {}).get("Volume", "")
        issue = journal.get("JournalIssue", {}).get("Issue", "")
        pages = art.get("Pagination", {}).get("MedlinePgn", "")

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

        abstract_parts = art.get("Abstract", {}).get("AbstractText", [])
        abstract = " ".join(str(p) for p in abstract_parts)

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
