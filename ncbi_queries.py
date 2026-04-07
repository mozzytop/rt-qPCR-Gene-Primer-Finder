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


# ── Species detection ────────────────────────────────────────────────

# NCBI Taxonomy IDs
_TXID_HUMAN = "9606"
_TXID_MOUSE = "10090"

_MOUSE_KEYWORDS = re.compile(
    r"\b(mouse|mus|murine|mus\s*musculus)\b", re.IGNORECASE
)


def detect_species(text: str) -> tuple[str, str, str]:
    """Auto-detect species from the user's input string.

    Scans for mouse-related keywords (mouse, mus, murine, mus musculus).
    Defaults to Homo sapiens if none are found.

    Returns (organism_name, taxonomy_id, short_label).
    """
    if _MOUSE_KEYWORDS.search(text):
        return "Mus musculus", _TXID_MOUSE, "mouse"
    return "Homo sapiens", _TXID_HUMAN, "human"


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


def search_nucleotide(gene: str, organism: str = "Homo sapiens", retmax: int = 20) -> list[dict]:
    """Search NCBI Nucleotide for *gene* in *organism*.

    Returns a list of dicts with: id, title, accession, length, link, score.
    Results are scored and sorted so the best "transcript variant 1, mRNA"
    record appears first, but ALL results are returned for user selection.
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
    """Pick the highest-scored result (already sorted by search_nucleotide)."""
    if not results:
        return None
    # Results are already sorted by score in search_nucleotide.
    # Return the first one that contains the gene name.
    gene_upper = gene.upper()
    for r in results:
        if gene_upper in r["title"].upper():
            return r
    return results[0]


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


def search_pmc_for_primers(
    gene: str,
    organism: str = "human",
    taxonomy_id: str = "",
    retmax: int = 20,
) -> list[dict]:
    """
    Search **PubMed Central** for open-access articles containing PCR
    primer sequences for *gene*, filtered by species.

    Uses taxonomy-aware organism name to prevent species bleeding.
    """
    # Resolve taxonomy ID and organism name
    if not taxonomy_id:
        _, taxonomy_id, _ = detect_species(organism)

    org_name = "Mus musculus" if taxonomy_id == _TXID_MOUSE else "Homo sapiens"

    query = (
        f'"{gene}" AND (primer OR "PCR" OR "RT-PCR" OR "qPCR" OR "real-time PCR")'
        f' AND ("{org_name}") AND open access[filter]'
    )
    handle = Entrez.esearch(db="pmc", term=query, retmax=retmax, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    ids = record.get("IdList", [])
    if not ids:
        # Fallback: try without open access filter
        query_fallback = (
            f'"{gene}" AND (primer OR "PCR" OR "RT-PCR" OR "qPCR")'
            f' AND ("{org_name}")'
        )
        handle = Entrez.esearch(db="pmc", term=query_fallback, retmax=retmax, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        ids = record.get("IdList", [])

    if not ids:
        return []

    # Fetch article summaries from PMC
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
    accession: str = "",
) -> list[VerifiedPrimerPair]:
    """
    Triple-Check Verification of primers from a PMC full-text article.

    The selected CDS sequence is the **biological anchor** (ground truth).
    A primer pair is only accepted if ALL THREE checks pass:

      Check 1 (Forward):  The Forward primer exists as an exact substring
                          of the CDS (either strand).
      Check 2 (Reverse):  The Reverse *Complement* of the Reverse primer
                          exists as an exact substring of the CDS.
      Check 3 (Species):  If the accession indicates a species (NM_ = human
                          RefSeq, NM_ + title context), reject primers whose
                          surrounding text mentions only the wrong species.
                          In practice, Checks 1+2 already enforce this:
                          mouse primers can't physically match a human CDS.

    Steps:
      1. Get full text from the article's raw_xml.
      2. Run direction-aware regex extraction.
      3. For each candidate pair, apply the Triple-Check.
      4. Return only fully verified pairs.
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

    # ── Normalise CDS for matching (biological anchor) ────────────────
    cds_upper = cds_sequence.upper().replace("U", "T")

    def _species_from_text(text: str) -> set[str]:
        labels: set[str] = set()
        txt = text.lower()
        if re.search(r"\b(homo\s+sapiens|human)\b", txt):
            labels.add("human")
        if re.search(r"\b(mus\s+musculus|mouse|murine)\b", txt):
            labels.add("mouse")
        return labels

    def _infer_expected_species() -> str:
        # Best effort from accession first, then gene symbol style, then article title.
        # RefSeq prefixes alone are not species-specific, so this remains conservative.
        acc = accession.strip().upper()
        if "HOMO SAPIENS" in acc or acc.startswith("NC_000"):
            return "human"
        if "MUS MUSCULUS" in acc or acc.startswith("NC_000067"):
            return "mouse"

        gene_clean = re.sub(r"[^A-Za-z0-9]", "", gene)
        if gene_clean and any(ch.isalpha() for ch in gene_clean):
            if gene_clean.upper() == gene_clean:
                return "human"
            if gene_clean[0:1].isupper() and gene_clean[1:].lower() == gene_clean[1:]:
                return "mouse"

        title_species = _species_from_text(article.get("title", ""))
        if len(title_species) == 1:
            return next(iter(title_species))
        return ""

    expected_species = _infer_expected_species()

    def _species_consistent(fwd_ctx: str, rev_ctx: str) -> bool:
        if not expected_species:
            return True
        pair_species = _species_from_text(fwd_ctx) | _species_from_text(rev_ctx)
        if not pair_species:
            pair_species = _species_from_text(article.get("title", ""))
        if not pair_species:
            return True
        if expected_species == "human":
            return "human" in pair_species and "mouse" not in pair_species
        return "mouse" in pair_species and "human" not in pair_species

    def _triple_check(fwd: ExtractedPrimer, rev: ExtractedPrimer) -> dict | None:
        """Return match metadata if all three checks pass, else None."""
        fwd_up = fwd.sequence.upper().replace("U", "T")
        rev_up = rev.sequence.upper().replace("U", "T")
        rev_comp_up = reverse_complement(rev_up).upper().replace("U", "T")

        # Check 1: candidate Forward exists exactly in selected CDS.
        fwd_pos = cds_upper.find(fwd_up)
        if fwd_pos == -1:
            return None

        # Check 2: reverse-complement(candidate Reverse) exists exactly in CDS.
        rev_pos = cds_upper.find(rev_comp_up)
        if rev_pos == -1:
            return None

        # Check 3: species consistency from local primer context.
        if not _species_consistent(fwd.context, rev.context):
            return None

        return {
            "fwd_pos": fwd_pos + 1,  # 1-indexed
            "fwd_strand": "sense",
            "rev_pos": rev_pos + 1,
            "rev_strand": "antisense",
            "rev_comp": rev_comp_up,
        }

    # ── Separate into forward / reverse / unknown ─────────────────────
    forwards = [c for c in candidates if c.direction == "forward"]
    reverses = [c for c in candidates if c.direction == "reverse"]
    unknowns = [c for c in candidates if c.direction == ""]

    # ── Try to assign unknown primers by CDS strand mapping ───────────
    for unk in unknowns:
        unk_up = unk.sequence.upper().replace("U", "T")
        if cds_upper.find(unk_up) != -1:
            forwards.append(unk)  # maps to sense strand → likely forward
        elif cds_upper.find(reverse_complement(unk_up).upper()) != -1:
            reverses.append(unk)  # reverse-complement maps to sense

    # ── Check all possible fwd × rev combinations ────────────────────
    verified_pairs: list[VerifiedPrimerPair] = []

    for fwd_candidate in forwards:
        for rev_candidate in reverses:
            match = _triple_check(fwd_candidate, rev_candidate)
            if match is not None:
                pair = VerifiedPrimerPair(
                    forward=fwd_candidate.sequence.upper(),
                    reverse=rev_candidate.sequence.upper(),
                    reverse_comp=match["rev_comp"],
                    fwd_position=match["fwd_pos"],
                    fwd_strand=match["fwd_strand"],
                    rev_position=match["rev_pos"],
                    rev_strand=match["rev_strand"],
                    source_pmcid=article.get("pmcid", ""),
                    source_pmid=article.get("pmid", ""),
                    citation=article.get("citation", ""),
                    title=article.get("title", ""),
                    doi=article.get("doi", ""),
                )
                verified_pairs.append(pair)

    # If no labelled fwd/rev but ≥2 unknowns map, try pairing by strand
    if not verified_pairs and len(unknowns) >= 2:
        sense_hits = []
        antisense_hits = []
        for unk in unknowns:
            unk_up = unk.sequence.upper().replace("U", "T")
            if cds_upper.find(unk_up) != -1:
                sense_hits.append((unk, cds_upper.find(unk_up) + 1))
            else:
                unk_rc = reverse_complement(unk_up).upper().replace("U", "T")
                rc_pos = cds_upper.find(unk_rc)
                if rc_pos != -1:
                    antisense_hits.append((unk, rc_pos + 1))
        for fwd_unk, _ in sense_hits:
            for rev_unk, _ in antisense_hits:
                match = _triple_check(fwd_unk, rev_unk)
                if match is not None:
                    pair = VerifiedPrimerPair(
                        forward=fwd_unk.sequence.upper(),
                        reverse=rev_unk.sequence.upper(),
                        reverse_comp=match["rev_comp"],
                        fwd_position=match["fwd_pos"],
                        fwd_strand="sense",
                        rev_position=match["rev_pos"],
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
