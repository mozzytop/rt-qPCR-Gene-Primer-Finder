"""
Gene Primer Lookup Tool – Streamlit Application (v3)
=====================================================

Features:
  - Tabbed UI: Main Workflow | Filter DNA | Sequence Calculator
  - Flexible input: Gene Name, NCBI URL, or raw pasted sequence
  - Variant selection from multiple NCBI results
  - Clickable DOI/PMCID reference links
  - Built-in Filter DNA utility (strips non-DNA characters)
  - Built-in Sequence Calculator (Reverse / Complement / RC)

Run with:
    python3 -m streamlit run app.py
"""

from __future__ import annotations

import re
import textwrap
from html import escape as html_escape
import streamlit as st
from Bio import Entrez

from utils import (
    build_report_docx,
    build_report_pdf,
    filter_dna,
    reverse_complement,
    complement_only,
    reverse_only,
    find_primer_in_sequence,
    format_origin_block,
    format_filtered_dna,
    find_primer_on_either_strand,
    extract_primers_from_text,
    extract_primers_with_direction,
    verify_primer_pair,
)
from ncbi_queries import (
    lookup_gene,
    search_nucleotide,
    parse_ncbi_url,
    fetch_cds_by_accession,
    fetch_genbank_record,
    extract_cds,
    search_pmc_for_primers,
    search_pmc_for_gene_mentions,
    search_pubmed_for_primers,
    extract_and_verify_primers,
    CDSResult,
    VerifiedPrimerPair,
)

# ── Page config ──────────────────────────────────────────────────────

st.set_page_config(
    page_title="Gene Primer Lookup Tool",
    layout="wide",
)

# ── Custom CSS ───────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

:root {
    --bg-primary: #0f1117;
    --bg-card: #1a1d29;
    --accent-1: #6c63ff;
    --accent-2: #00d2ff;
    --accent-3: #7c3aed;
    --text-primary: #e4e4e7;
    --text-secondary: #a1a1aa;
    --border-color: #27273a;
    --success: #22c55e;
    --warning: #f59e0b;
    --error: #ef4444;
}
html, body, [class*="css"] { font-family: 'Inter', sans-serif; }

.main-header {
    background: linear-gradient(135deg, #6c63ff 0%, #00d2ff 50%, #7c3aed 100%);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent;
    background-clip: text; font-weight: 700; font-size: 2.4rem;
    letter-spacing: -0.02em; margin-bottom: 0.1rem;
}
.sub-header { color: var(--text-secondary); font-size: 1.05rem; font-weight: 300; margin-bottom: 1.2rem; }

.card {
    background: var(--bg-card); border: 1px solid var(--border-color);
    border-radius: 16px; padding: 1.5rem; margin-bottom: 1rem;
    transition: border-color 0.3s ease, box-shadow 0.3s ease;
}
.card:hover { border-color: var(--accent-1); box-shadow: 0 0 20px rgba(108,99,255,0.08); }
.card-title {
    font-size: 1.1rem; font-weight: 600; color: var(--text-primary);
    margin-bottom: 0.8rem; display: flex; align-items: center; gap: 0.5rem;
}

.sequence-box {
    background: #111318; border: 1px solid #2a2d3a; border-radius: 10px;
    padding: 1rem 1.2rem; font-family: 'JetBrains Mono', monospace;
    font-size: 0.78rem; line-height: 1.6; color: #c4c4cc;
    overflow-x: auto; white-space: pre; max-height: 400px; overflow-y: auto;
}

.primer-badge {
    display: inline-block; padding: 0.35rem 0.85rem; border-radius: 8px;
    font-family: 'JetBrains Mono', monospace; font-size: 0.82rem;
    font-weight: 500; letter-spacing: 0.03em;
}
.primer-fwd { background: rgba(34,197,94,0.12); border: 1px solid rgba(34,197,94,0.3); color: #22c55e; }
.primer-fwd-verification { background: rgba(168,85,247,0.16); border: 1px solid rgba(168,85,247,0.42); color: #d8b4fe; }
.primer-rev { background: rgba(14,165,233,0.12); border: 1px solid rgba(14,165,233,0.3); color: #0ea5e9; }
.primer-rc { background: rgba(245,158,11,0.12); border: 1px solid rgba(245,158,11,0.3); color: #f59e0b; }
.primer-warning { background: rgba(255, 8, 32, 0.33); border: 1px solid rgba(255, 8, 32, 0.95); color: #f50b0b; }

.status-pill {
    display: inline-flex; align-items: center; gap: 0.4rem;
    padding: 0.25rem 0.7rem; border-radius: 999px; font-size: 0.78rem; font-weight: 500;
}
.status-success { background: rgba(34,197,94,0.12); color: #22c55e; }
.status-warning { background: rgba(245,158,11,0.12); color: #f59e0b; }
.status-error   { background: rgba(239,68,68,0.12); color: #ef4444; }

.step-indicator {
    display: inline-flex; align-items: center; justify-content: center;
    width: 28px; height: 28px; border-radius: 50%;
    background: linear-gradient(135deg, var(--accent-1), var(--accent-2));
    color: white; font-size: 0.8rem; font-weight: 700; margin-right: 0.6rem;
}

.verified-pair-box {
    background: linear-gradient(135deg, rgba(34,197,94,0.06), rgba(108,99,255,0.06));
    border: 1px solid rgba(34,197,94,0.25); border-radius: 14px;
    padding: 1.3rem 1.5rem; margin: 0.8rem 0;
    font-family: 'JetBrains Mono', monospace; font-size: 0.82rem; line-height: 2;
}
.verified-pair-box .label {
    color: #a1a1aa; font-family: 'Inter', sans-serif; font-size: 0.72rem;
    text-transform: uppercase; letter-spacing: 0.08em; font-weight: 600;
}

.report-box {
    background: #0c0e14; border: 1px solid #6c63ff44; border-radius: 12px;
    padding: 1.5rem; font-family: 'JetBrains Mono', monospace; font-size: 0.76rem;
    line-height: 1.65; color: #d4d4dc; overflow-x: auto; white-space: normal; overflow-wrap: anywhere;
}
.report-box strong { color: #ffffff; font-weight: 700; }
.report-pre {
    white-space: pre-wrap; word-break: normal; overflow-wrap: anywhere;
    margin: 0.2rem 0 1rem 0;
}
.report-line { margin-bottom: 0.6rem; }
.report-reference { white-space: pre-wrap; overflow-wrap: anywhere; }
.report-highlight-fwd {
    background: rgba(34,197,94,0.22); color: #9df0b7; border-radius: 4px; padding: 0 1px;
}
.report-highlight-rc {
    background: rgba(108,99,255,0.24); color: #c5b8ff; border-radius: 4px; padding: 0 1px;
}

div.stButton > button {
    background: linear-gradient(135deg, #6c63ff, #7c3aed) !important;
    color: white !important; border: none !important; border-radius: 10px !important;
    padding: 0.55rem 1.8rem !important; font-weight: 600 !important;
    transition: all 0.3s ease !important;
}
div.stButton > button:hover {
    box-shadow: 0 4px 20px rgba(108,99,255,0.35) !important;
    transform: translateY(-1px) !important;
}
.stTextArea textarea {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.8rem !important; border-radius: 10px !important;
}

.gradient-divider {
    height: 2px; background: linear-gradient(90deg, transparent, #6c63ff, #00d2ff, transparent);
    border: none; margin: 2rem 0; border-radius: 1px;
}

/* Utility result cards */
.util-result {
    background: #111318; border: 1px solid #2a2d3a; border-radius: 12px;
    padding: 1rem 1.3rem; margin: 0.6rem 0;
}
.util-result-label {
    font-size: 0.72rem; font-weight: 600; text-transform: uppercase;
    letter-spacing: 0.08em; margin-bottom: 0.4rem;
}
.util-result-seq {
    font-family: 'JetBrains Mono', monospace; font-size: 0.8rem;
    line-height: 1.6; word-break: break-all; white-space: pre-wrap;
}
.label-reverse { color: #f59e0b; }
.label-complement { color: #00d2ff; }
.label-rc { color: #6c63ff; }
.label-filtered { color: #22c55e; }

/* Variant table styling */
.variant-table {
    width: 100%; border-collapse: separate; border-spacing: 0;
    font-size: 0.82rem;
}
.variant-table th {
    background: #1a1d29; color: #a1a1aa; font-weight: 600; font-size: 0.72rem;
    text-transform: uppercase; letter-spacing: 0.06em; padding: 0.7rem 1rem;
    border-bottom: 1px solid #27273a; text-align: left;
}
.variant-table td {
    padding: 0.6rem 1rem; border-bottom: 1px solid #1a1d29; color: #d4d4dc;
}
.variant-table tr:hover td { background: #1a1d2944; }
.variant-table a { color: #6c63ff; text-decoration: none; }
.variant-table a:hover { color: #00d2ff; text-decoration: underline; }
</style>
""", unsafe_allow_html=True)

# ── Header ───────────────────────────────────────────────────────────

st.markdown('<p class="main-header">Gene Primer Lookup Tool</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">'
    "Automated NCBI gene lookup · CDS extraction · PMC primer search · Sequence utilities"
    "</p>",
    unsafe_allow_html=True,
)
st.markdown(
    '<p class="sub-header">'
    "MANUALLY VERIFY SEQUENCES AND SOURCES AS THIS PROGRAM CAN MAKE MISTAKES!"
    "</p>",
    unsafe_allow_html=True,
)


def _has_secret_ncbi_credentials() -> bool:
    """Return True when Streamlit secrets contain an NCBI email."""
    if "ncbi" in st.secrets:
        ncbi_secrets = st.secrets["ncbi"]
        return bool(
            str(ncbi_secrets.get("ncbi_email") or ncbi_secrets.get("email") or "").strip()
        )
    return bool(str(st.secrets.get("ncbi_email", "") or "").strip())


# ── Sidebar config ───────────────────────────────────────────────────

with st.sidebar:
    st.markdown("### Configuration")
    if _has_secret_ncbi_credentials():
        email = ""
        api_key = ""
        st.success("Logged in as ESFCOM Lab Member")
    else:
        email = st.text_input(
            "NCBI Email (required)",
            placeholder="you@example.com",
            help="NCBI requires an email for Entrez API usage.",
        )
        api_key = st.text_input(
            "NCBI API Key (optional)",
            type="password",
            help="An API key raises the rate limit from 3 to 10 requests/sec.",
        )
    organism = st.selectbox("Organism", ["Homo sapiens", "Mus musculus"], index=0)

    st.markdown("---")
    st.markdown(
        "**How it works:**\n\n"
        "1. Enter a gene name, NCBI URL, or sequence\n"
        "2. Select a transcript variant\n"
        "3. Auto-search PMC for primer pairs\n"
        "4. Verify primers & generate report\n\n"
        "**Utilities:**\n"
        "- Filter DNA — clean raw sequences\n"
        "- Sequence Calculator — RC, complement"
    )

# ── Session state init ──────────────────────────────────────────────

for key, default in {
    "cds": None,
    "gene": "",
    "input_type": "gene_name",
    "summaries": [],
    "variant_result_limit": 10,
    "selected_accession": "",
    "pmc_articles": [],
    "verified_pairs": [],
    "final_report": "",
    "final_report_html": "",
    "fwd_primer": "",
    "rev_primer": "",
    "reference": "",
    "auto_search_done": False,
}.items():
    if key not in st.session_state:
        st.session_state[key] = default

# ── Helpers ──────────────────────────────────────────────────────────


def _get_ncbi_credentials() -> tuple[str, str]:
    """Prefer Streamlit secrets when present, otherwise use manual sidebar input."""
    secret_email = ""
    secret_api_key = ""

    if "ncbi" in st.secrets:
        ncbi_secrets = st.secrets["ncbi"]
        secret_email = str(
            ncbi_secrets.get("ncbi_email") or ncbi_secrets.get("email") or ""
        ).strip()
        secret_api_key = str(
            ncbi_secrets.get("ncbi_api_key") or ncbi_secrets.get("api_key") or ""
        ).strip()
    else:
        secret_email = str(st.secrets.get("ncbi_email", "") or "").strip()
        secret_api_key = str(
            st.secrets.get("ncbi_api_key", "") or st.secrets.get("api_key", "") or ""
        ).strip()

    if secret_email:
        return secret_email, secret_api_key

    manual_email = email.strip()
    manual_api_key = api_key.strip()
    if manual_email:
        return manual_email, manual_api_key

    return "", ""


def _set_entrez():
    resolved_email, resolved_api_key = _get_ncbi_credentials()

    if not resolved_email:
        st.error("Please enter your NCBI email in the sidebar or add it to Streamlit secrets.")
        st.stop()
    Entrez.email = resolved_email
    Entrez.api_key = resolved_api_key or None



def _detect_input_type(text: str) -> str:
    """Detect whether input is a URL, raw sequence, or gene name."""
    text = text.strip()
    if re.match(r"https?://", text):
        return "url"
    # If >50% of characters are ATCG, treat as raw sequence
    dna_chars = sum(1 for c in text.upper() if c in "ATCGU")
    if len(text) > 20 and dna_chars / max(len(text.replace(" ", "").replace("\n", "")), 1) > 0.7:
        return "raw_sequence"
    return "gene_name"


def _make_ref_links(pair: VerifiedPrimerPair) -> str:
    """Build clickable reference link string for a verified primer pair."""
    parts = []
    if pair.citation:
        parts.append(pair.citation)
    links = []
    if pair.doi:
        doi_url = f"https://doi.org/{pair.doi}" if not pair.doi.startswith("http") else pair.doi
        links.append(f"[DOI]({doi_url})")
    if pair.source_pmcid:
        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pair.source_pmcid}/"
        links.append(f"[{pair.source_pmcid}]({pmc_url})")
    if pair.source_pmid:
        pm_url = f"https://pubmed.ncbi.nlm.nih.gov/{pair.source_pmid}/"
        links.append(f"[PMID:{pair.source_pmid}]({pm_url})")
    if links:
        parts.append(" · ".join(links))
    return "\n\n".join(parts) if parts else ""


def _scan_pmc_articles_for_pairs(
    articles: list[dict],
    gene_name: str,
    filtered_cds: str,
    progress_label: str,
) -> list[VerifiedPrimerPair]:
    """Extract and verify primers across a batch of PMC articles."""
    verified_pairs: list[VerifiedPrimerPair] = []
    progress = st.progress(0, text=progress_label)

    for i, art in enumerate(articles):
        title = art.get("title", "Untitled article")
        progress.progress(
            (i + 1) / len(articles),
            text=f"Scanning article {i + 1}/{len(articles)}: {title[:60]}…",
        )
        try:
            pairs = extract_and_verify_primers(art, gene_name, filtered_cds)
            verified_pairs.extend(pairs)
        except Exception:
            continue

    progress.empty()
    return verified_pairs


def _find_all_occurrences(sequence: str, motif: str) -> list[tuple[int, int]]:
    """Return all non-overlapping motif spans within sequence."""
    spans: list[tuple[int, int]] = []
    if not motif:
        return spans
    start = 0
    while True:
        idx = sequence.find(motif, start)
        if idx == -1:
            break
        spans.append((idx, idx + len(motif)))
        start = idx + len(motif)
    return spans


def _build_highlight_spans(sequence: str, fwd: str, rev_comp: str) -> list[tuple[int, int, str]]:
    """Return sequence spans for forward and reverse-complement highlights."""
    spans: list[tuple[int, int, str]] = []
    occupied: set[int] = set()
    for start, end in _find_all_occurrences(sequence, fwd.lower()):
        spans.append((start, end, "report-highlight-fwd"))
        occupied.update(range(start, end))
    for start, end in _find_all_occurrences(sequence, rev_comp.lower()):
        if any(idx in occupied for idx in range(start, end)):
            continue
        spans.append((start, end, "report-highlight-rc"))
    spans.sort(key=lambda item: item[0])
    return spans


def _apply_text_highlights(sequence: str, spans: list[tuple[int, int, str]]) -> str:
    """Uppercase highlighted sequence spans for plain-text exports."""
    chars = list(sequence.lower())
    for start, end, _ in spans:
        for idx in range(start, min(end, len(chars))):
            chars[idx] = chars[idx].upper()
    return "".join(chars)


def _format_origin_preserving_case(sequence: str, line_width: int = 60, block_size: int = 10) -> str:
    """Format an ORIGIN block without forcing the sequence to lowercase."""
    lines: list[str] = []
    for i in range(0, len(sequence), line_width):
        chunk = sequence[i : i + line_width]
        blocks = [chunk[j : j + block_size] for j in range(0, len(chunk), block_size)]
        lines.append(f"{i + 1:>9} {' '.join(blocks)}")
    return "\n".join(lines)


def _render_sequence_html(
    sequence: str,
    spans: list[tuple[int, int, str]],
    *,
    line_width: int = 60,
    block_size: int | None = None,
    show_origin_numbers: bool = False,
) -> str:
    """Render a highlighted DNA sequence block as HTML."""
    ordered = sorted(spans, key=lambda item: item[0])

    def render_range(start: int, end: int) -> str:
        parts: list[str] = []
        cursor = start
        for span_start, span_end, css_class in ordered:
            if span_end <= start or span_start >= end:
                continue
            seg_start = max(start, span_start)
            seg_end = min(end, span_end)
            if cursor < seg_start:
                parts.append(html_escape(sequence[cursor:seg_start]))
            parts.append(
                f'<span class="{css_class}">{html_escape(sequence[seg_start:seg_end])}</span>'
            )
            cursor = seg_end
        if cursor < end:
            parts.append(html_escape(sequence[cursor:end]))
        return "".join(parts)

    lines: list[str] = []
    for line_start in range(0, len(sequence), line_width):
        line_end = min(line_start + line_width, len(sequence))
        line_html_parts: list[str] = []
        if show_origin_numbers:
            line_html_parts.append(f"{line_start + 1:>9} ")
        if block_size:
            for block_start in range(line_start, line_end, block_size):
                block_end = min(block_start + block_size, line_end)
                line_html_parts.append(render_range(block_start, block_end))
                if block_end < line_end:
                    line_html_parts.append(" ")
        else:
            line_html_parts.append(render_range(line_start, line_end))
        lines.append("".join(line_html_parts))
    return "<br>".join(lines)


def _build_report(
    gene: str, cds: CDSResult, fwd: str, rev: str, reference: str, organism_label: str,
) -> tuple[str, str]:
    """Build plain-text and HTML report outputs."""
    filtered = filter_dna(cds.cds_dna)
    rev_comp = reverse_complement(rev)
    link = f"https://www.ncbi.nlm.nih.gov/nuccore/{cds.accession}?from={cds.cds_start}&to={cds.cds_end}"
    org_label = "Human" if "homo" in organism_label.lower() else "Mouse"
    spans = _build_highlight_spans(filtered.lower(), fwd.lower(), rev_comp.lower())
    highlighted_sequence = _apply_text_highlights(filtered, spans)
    origin_block = _format_origin_preserving_case(highlighted_sequence)
    filtered_block = format_filtered_dna(highlighted_sequence)
    origin_block_html = _render_sequence_html(
        filtered.lower(),
        spans,
        line_width=60,
        block_size=10,
        show_origin_numbers=True,
    )
    filtered_block_html = _render_sequence_html(filtered.lower(), spans, line_width=60)
    wrapped_reference = (
        textwrap.fill(reference, width=90, break_long_words=True, break_on_hyphens=True)
        if reference else ""
    )

    lines: list[str] = []
    lines.append(f"Gene: {gene.upper()}")
    lines.append(f"Link: {link}")
    lines.append("")
    lines.append("ORIGIN")
    lines.append(origin_block)
    lines.append("")

    trans = cds.translation
    if trans:
        lines.append('translation="')
        lines.extend(textwrap.wrap(trans, width=60))
        lines[-1] += '"'
        lines.append("")

    lines.append("Filter DNA results:")
    lines.append(f">filtered DNA sequence consisting of {len(filtered)} bases.")
    lines.append(filtered_block)
    lines.append("")

    lines.append(f"{org_label} {gene.upper()} Primer Sequence: Forward : 5'-{fwd}-3'")
    lines.append(f"{org_label} {gene.upper()} Primer Sequence: Reverse : 5'-{rev}-3'")
    lines.append(f"{org_label} {gene.upper()} Primer Sequence: Reverse Comp. : 5'-{rev_comp}-3'")

    if wrapped_reference:
        lines.append("")
        lines.append("Reference:")
        lines.append(wrapped_reference)

    html_parts: list[str] = [
        f'<div class="report-line"><strong>Gene:</strong> {html_escape(gene.upper())}</div>',
        f'<div class="report-line"><strong>Link:</strong> {html_escape(link)}</div>',
        '<div class="report-line"><strong>ORIGIN</strong></div>',
        f'<div class="report-pre">{origin_block_html}</div>',
    ]
    if trans:
        trans_html = "<br>".join(html_escape(chunk) for chunk in textwrap.wrap(trans, width=60))
        html_parts.extend([
            '<div class="report-line"><strong>translation</strong></div>',
            f'<div class="report-pre">&quot;{trans_html}&quot;</div>',
        ])
    html_parts.extend([
        '<div class="report-line"><strong>Filter DNA results:</strong></div>',
        f'<div class="report-line">&gt; filtered DNA sequence consisting of {len(filtered)} bases.</div>',
        f'<div class="report-pre">{filtered_block_html}</div>',
        f'<div class="report-line"><strong>{html_escape(org_label)} {html_escape(gene.upper())} Primer Sequence: Forward :</strong> 5&#39;-{html_escape(fwd)}-3&#39;</div>',
        f'<div class="report-line">{html_escape(org_label)} {html_escape(gene.upper())} Primer Sequence: Reverse : 5&#39;-{html_escape(rev)}-3&#39;</div>',
        f'<div class="report-line"><strong>{html_escape(org_label)} {html_escape(gene.upper())} Primer Sequence: Reverse Comp. :</strong> 5&#39;-{html_escape(rev_comp)}-3&#39;</div>',
    ])
    if wrapped_reference:
        html_parts.extend([
            '<div class="report-line"><strong>Reference:</strong></div>',
            f'<div class="report-reference">{html_escape(wrapped_reference)}</div>',
        ])

    return "\n".join(lines), "".join(html_parts)


# ═══════════════════════════════════════════════════════════════════
# TABS
# ═══════════════════════════════════════════════════════════════════

tab_main, tab_filter, tab_calc = st.tabs([
    "Primer Finder",
    "Filter DNA",
    "Sequence Calculator",
])

# ═══════════════════════════════════════════════════════════════════
# TAB 1 — MAIN WORKFLOW
# ═══════════════════════════════════════════════════════════════════

with tab_main:

    # ── STEP 1 — Flexible Input ──────────────────────────────────────

    st.markdown(
        '<div class="card"><div class="card-title">'
        '<span class="step-indicator">1</span> Gene Input &amp; Variant Selection'
        "</div></div>",
        unsafe_allow_html=True,
    )

    st.caption(
        "Enter a **gene name** (e.g. GLI1), a **direct NCBI URL** "
        "(e.g. https://www.ncbi.nlm.nih.gov/nuccore/NM_005269.3), "
        "or paste a **raw DNA sequence**."
    )

    with st.form(key="search_form"):
        col_input, col_btn = st.columns([4, 1])

        with col_input:
            user_input = st.text_input(
                "Input",
                placeholder="GLI1  or  https://www.ncbi.nlm.nih.gov/nuccore/NM_005269.3  or  ATGTTCAACT...",
                label_visibility="collapsed",
            )
        with col_btn:
            search_btn = st.form_submit_button("Search", use_container_width=True)

    if search_btn and user_input.strip():
        _set_entrez()
        inp = user_input.strip()
        input_type = _detect_input_type(inp)

        # Reset downstream state
        for k in ["pmc_articles", "verified_pairs", "final_report", "final_report_html", "fwd_primer",
                   "rev_primer", "reference", "auto_search_done", "selected_accession"]:
            st.session_state[k] = "" if isinstance(st.session_state.get(k), str) else (
                [] if isinstance(st.session_state.get(k), list) else False
            )
        st.session_state.cds = None

        if input_type == "url":
            # ── Direct NCBI URL ──────────────────────────────────────
            accession = parse_ncbi_url(inp)
            if accession:
                st.session_state.gene = accession
                st.session_state.input_type = "url"
                with st.spinner(f"Fetching GenBank record for {accession}…"):
                    try:
                        cds = fetch_cds_by_accession(accession)
                        st.session_state.cds = cds
                        st.session_state.selected_accession = accession
                    except Exception as e:
                        st.error(f"Failed to fetch record: {e}")
            else:
                st.error("Could not parse an accession from that URL.")

        elif input_type == "raw_sequence":
            # ── Raw pasted sequence ──────────────────────────────────
            st.session_state.input_type = "raw_sequence"
            clean = filter_dna(inp)
            if len(clean) < 30:
                st.error("The pasted sequence is too short for CDS analysis.")
            else:
                st.session_state.gene = "USER_SEQ"
                cds = CDSResult(
                    accession="USER_INPUT",
                    description="User-pasted sequence",
                    cds_dna=clean,
                    cds_start=1,
                    cds_end=len(clean),
                    translation="",
                    full_sequence=clean,
                )
                st.session_state.cds = cds

        else:
            # ── Gene name search ─────────────────────────────────────
            st.session_state.gene = inp
            st.session_state.input_type = "gene_name"
            with st.spinner("Querying NCBI Nucleotide…"):
                summaries = search_nucleotide(inp, organism, retmax=100)
                st.session_state.summaries = summaries

    # ── Variant selection table ──────────────────────────────────────

    if (st.session_state.input_type == "gene_name"
            and st.session_state.summaries
            and st.session_state.cds is None):

        display_limit = int(st.session_state.variant_result_limit)
        results = st.session_state.summaries[:display_limit]
        st.markdown(
            f"Showing **{len(results)}** of **{len(st.session_state.summaries)}** results. Select a variant:"
        )

        # Build HTML table
        table_html = '<table class="variant-table"><thead><tr>'
        table_html += "<th></th><th>Accession</th><th>Description</th><th>Length (bp)</th><th>NCBI Link</th>"
        table_html += "</tr></thead><tbody>"
        for i, r in enumerate(results):
            table_html += f"<tr>"
            table_html += f'<td style="text-align:center;">{i + 1}</td>'
            table_html += f'<td><code>{r["accession"]}</code></td>'
            table_html += f'<td>{r["title"]}</td>'
            table_html += f'<td style="text-align:right;">{r["length"]:,}</td>'
            table_html += f'<td><a href="{r["link"]}" target="_blank">View</a></td>'
            table_html += "</tr>"
        table_html += "</tbody></table>"
        st.markdown(table_html, unsafe_allow_html=True)

        st.selectbox(
            "How many gene input results should be shown?",
            options=[10, 50, 100],
            key="variant_result_limit",
            help="This adjusts the Gene Input & Variant Selection result count.",
        )

        # Selection UI
        options = [f'{r["accession"]}  —  {r["title"][:80]}  ({r["length"]:,} bp)' for r in results]
        with st.form(key="variant_selection_form"):
            selected = st.selectbox("Choose a variant", options, index=0)
            select_btn = st.form_submit_button("Use Selected Variant", use_container_width=True)

        if select_btn:
            _set_entrez()
            idx = options.index(selected)
            chosen = results[idx]
            with st.spinner(f'Fetching CDS for {chosen["accession"]}…'):
                try:
                    cds = fetch_cds_by_accession(chosen["accession"])
                    st.session_state.cds = cds
                    st.session_state.selected_accession = chosen["accession"]
                except Exception as e:
                    st.error(f"Failed to extract CDS: {e}")

    # ── CDS display ──────────────────────────────────────────────────

    cds: CDSResult | None = st.session_state.cds

    if cds is not None:
        filtered = filter_dna(cds.cds_dna)
        link = f"https://www.ncbi.nlm.nih.gov/nuccore/{cds.accession}?from={cds.cds_start}&to={cds.cds_end}"

        st.markdown(
            f'<span class="status-pill status-success">● Found</span> &nbsp; '
            f"**{cds.description}**",
            unsafe_allow_html=True,
        )
        if cds.accession != "USER_INPUT":
            st.markdown(f"[NCBI Link]({link})")
        st.markdown(f"CDS span: **{cds.cds_start}..{cds.cds_end}** &nbsp;|&nbsp; Filtered bases: **{len(filtered)}**")

        with st.expander("ORIGIN (GenBank-formatted CDS)", expanded=False):
            st.markdown(f'<div class="sequence-box">{format_origin_block(filtered)}</div>', unsafe_allow_html=True)

        if cds.translation:
            with st.expander("Amino-Acid Translation", expanded=False):
                st.markdown(f'<div class="sequence-box">{cds.translation}</div>', unsafe_allow_html=True)

        with st.expander("Filtered DNA (continuous)", expanded=False):
            st.markdown(f'<div class="sequence-box">{format_filtered_dna(filtered)}</div>', unsafe_allow_html=True)

        st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

        # ── STEP 2 — PMC Primer Search ───────────────────────────────

        st.markdown(
            '<div class="card"><div class="card-title">'
            '<span class="step-indicator">2</span> PMC Full-Text Primer Search &amp; Auto-Verification'
            "</div></div>",
            unsafe_allow_html=True,
        )

        st.caption(
            "Searches PubMed Central for open-access papers with PCR primers. "
            "Extracts primer sequences from full-text, detects Forward/Reverse, "
            "and verifies them against the CDS."
        )

        pmc_btn = st.button("Search PMC for Primers", use_container_width=True)

        if pmc_btn:
            _set_entrez()
            org_short = "human" if "homo" in organism.lower() else "mouse"
            gene_name = st.session_state.gene

            with st.spinner("Searching PubMed Central…"):
                articles = search_pmc_for_primers(gene_name, org_short, retmax=20)

            all_articles = list(articles)
            all_verified: list[VerifiedPrimerPair] = []

            if not articles:
                with st.spinner("No focused PMC hits found. Expanding to a broader PMC gene search…"):
                    all_articles = search_pmc_for_gene_mentions(gene_name, retmax=100)

                if not all_articles:
                    st.session_state.pmc_articles = []
                    st.session_state.verified_pairs = []
                    st.session_state.auto_search_done = True
                    st.warning("No PMC articles found. Try manual primer input in Step 3.")
                else:
                    st.info(
                        f"Focused PMC search returned no hits. Scanning **{len(all_articles)}** "
                        f"broader full-text articles mentioning **{gene_name.upper()}**…"
                    )
                    all_verified = _scan_pmc_articles_for_pairs(
                        all_articles,
                        gene_name,
                        filtered,
                        "Scanning broader PMC full-text results…",
                    )

                    st.session_state.pmc_articles = all_articles
                    st.session_state.verified_pairs = all_verified
                    st.session_state.auto_search_done = True

                    if all_verified:
                        best = all_verified[0]
                        st.session_state.fwd_primer = best.forward
                        st.session_state.rev_primer = best.reverse
                        st.session_state.reference = best.citation
            else:
                st.info(f"Found **{len(articles)}** PMC articles. Scanning for primers…")
                all_verified = _scan_pmc_articles_for_pairs(
                    articles,
                    gene_name,
                    filtered,
                    "Extracting & verifying primers…",
                )

                if not all_verified:
                    with st.spinner("Expanding to a broader PMC gene search…"):
                        broader_articles = search_pmc_for_gene_mentions(gene_name, retmax=100)

                    seen_article_ids = {
                        art.get("pmcid") or art.get("pmc_id") or art.get("title", "")
                        for art in all_articles
                    }
                    extra_articles = [
                        art for art in broader_articles
                        if (art.get("pmcid") or art.get("pmc_id") or art.get("title", ""))
                        not in seen_article_ids
                    ]

                    if extra_articles:
                        st.info(
                            f"No verified pairs found in the focused PMC search. "
                            f"Scanning **{len(extra_articles)}** additional full-text articles mentioning "
                            f"**{gene_name.upper()}**…"
                        )
                        all_articles.extend(extra_articles)
                        all_verified.extend(
                            _scan_pmc_articles_for_pairs(
                                extra_articles,
                                gene_name,
                                filtered,
                                "Scanning broader PMC full-text results…",
                            )
                        )

                st.session_state.pmc_articles = all_articles
                st.session_state.verified_pairs = all_verified
                st.session_state.auto_search_done = True

                if all_verified:
                    best = all_verified[0]
                    st.session_state.fwd_primer = best.forward
                    st.session_state.rev_primer = best.reverse
                    st.session_state.reference = best.citation

        # ── Display verified results ─────────────────────────────────

        if st.session_state.verified_pairs:
            pairs = st.session_state.verified_pairs
            st.success(f"Found **{len(pairs)}** verified primer pair(s)!")

            org_label = "Human" if "homo" in organism.lower() else "Mouse"
            gene_upper = st.session_state.gene.upper()

            for idx, pair in enumerate(pairs):
                exp_title = f"Pair {idx + 1}  —  {pair.source_pmcid}  |  {pair.title[:70]}"
                with st.expander(exp_title, expanded=(idx == 0)):
                    st.markdown(
                        f"""<div class="verified-pair-box">
<span class="label">Forward Primer</span><br>
<span class="primer-badge primer-fwd">{org_label} {gene_upper} Primer Sequence: Forward : 5'-{pair.forward}-3'</span>
<span class="status-pill status-success">● CDS {pair.fwd_strand}, pos {pair.fwd_position}</span><br><br>
<span class="label">Reverse Primer</span><br>
<span class="primer-badge primer-rev">{org_label} {gene_upper} Primer Sequence: Reverse : 5'-{pair.reverse}-3'</span>
<span class="status-pill status-success">● CDS {pair.rev_strand}, pos {pair.rev_position}</span><br><br>
<span class="label">Reverse Complement</span><br>
<span class="primer-badge primer-rc">{org_label} {gene_upper} Primer Sequence: Reverse Comp. : 5'-{pair.reverse_comp}-3'</span>
</div>""",
                        unsafe_allow_html=True,
                    )

                    # Clickable reference links
                    ref_md = _make_ref_links(pair)
                    if ref_md:
                        st.markdown(f"**Reference:** {ref_md}")

                    if st.button(f"Select this pair for report (Double Check Species in PMC Report)", key=f"use_pair_{idx}"):
                        st.session_state.fwd_primer = pair.forward
                        st.session_state.rev_primer = pair.reverse
                        st.session_state.reference = pair.citation
                        st.rerun()

        elif st.session_state.auto_search_done and not st.session_state.verified_pairs:
            st.warning(
                "No verified primer pairs found in PMC full text. "
                "Use **Step 3** below to paste primers manually."
            )

        if st.session_state.pmc_articles:
            with st.expander("All PMC articles scanned", expanded=False):
                for art in st.session_state.pmc_articles:
                    pmcid = art.get("pmcid", "")
                    doi = art.get("doi", "")
                    links_parts = []
                    if pmcid:
                        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                        links_parts.append(f"[{pmcid}]({pmc_url})")
                    if doi:
                        doi_url = f"https://doi.org/{doi}" if not doi.startswith("http") else doi
                        links_parts.append(f"[DOI]({doi_url})")
                    link_str = " · ".join(links_parts) if links_parts else ""
                    st.markdown(f"- **{pmcid}** | {art['title']} {link_str}")

        st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

        # ── STEP 3 — Manual Primer Input ─────────────────────────────

        st.markdown(
            '<div class="card"><div class="card-title">'
            '<span class="step-indicator">3</span> Primer Input &amp; Verification'
            "</div>"
            '<span style="color:#a1a1aa;font-size:0.82rem;">'
            "Auto-populated from Step 2 if a verified pair was found. Edit or enter manually."
            "</span></div>",
            unsafe_allow_html=True,
        )

        pcol1, pcol2 = st.columns(2)
        with pcol1:
            fwd_input = st.text_input(
                "Forward Primer",
                value=st.session_state.fwd_primer,
                placeholder="e.g. AGCGTGAGCCTGAATCTGTG",
            )
        with pcol2:
            rev_input = st.text_input(
                "Reverse Primer",
                value=st.session_state.rev_primer,
                placeholder="e.g. CAGCATGTACTGGGCTTTGAA",
            )

        ref_input = st.text_area(
            "Reference Citation",
            value=st.session_state.reference,
            height=80,
            placeholder="Paste the full citation here…",
        )

        with st.expander("Or paste text containing 5'-...-3' primer notation", expanded=False):
            raw_snippet = st.text_area(
                "Paste abstract / methods text",
                height=120,
                placeholder="The tool will regex-extract primer sequences from this text…",
                key="raw_snippet",
            )
            if raw_snippet:
                found_directed = extract_primers_with_direction(raw_snippet)
                if found_directed:
                    st.markdown("**Extracted sequences:**")
                    for ep in found_directed:
                        dir_label = f" [{ep.direction}]" if ep.direction else ""
                        if ep.direction == "forward":
                            pos = find_primer_in_sequence(ep.sequence, filtered)
                            tag = (
                                f'<span class="status-pill status-success">● Matches sense strand, pos {pos}</span>'
                                if pos is not None
                                else '<span class="status-pill status-warning">● No valid forward-primer CDS match</span>'
                            )
                        elif ep.direction == "reverse":
                            pos = find_primer_in_sequence(ep.sequence, reverse_complement(filtered))
                            tag = (
                                f'<span class="status-pill status-success">● Matches antisense strand, pos {pos}</span>'
                                if pos is not None
                                else '<span class="status-pill status-warning">● No valid reverse-primer CDS match</span>'
                            )
                        else:
                            pos, strand = find_primer_on_either_strand(ep.sequence, filtered)
                            tag = (
                                f'<span class="status-pill status-success">● Maps ({strand}, pos {pos})</span>'
                                if pos is not None
                                else '<span class="status-pill status-warning">● No CDS match</span>'
                            )
                        st.markdown(f'`{ep.sequence}`{dir_label} {tag}', unsafe_allow_html=True)

                    fwd_c = [ep for ep in found_directed if ep.direction == "forward"]
                    rev_c = [ep for ep in found_directed if ep.direction == "reverse"]
                    if fwd_c and rev_c:
                        result = verify_primer_pair(fwd_c[0].sequence, rev_c[0].sequence, filtered)
                        if result["both_map"]:
                            st.success("Verified pair found in pasted text!")
                            if st.button("Use this pair"):
                                st.session_state.fwd_primer = result["forward_seq"]
                                st.session_state.rev_primer = result["reverse_seq"]
                                st.rerun()
                else:
                    st.warning("No primer patterns found in the pasted text.")

        if fwd_input and rev_input:
            fwd_clean = filter_dna(fwd_input.strip()).upper()
            rev_clean = filter_dna(rev_input.strip()).upper()

            st.markdown("#### Verification")
            result = verify_primer_pair(fwd_clean, rev_clean, filtered)

            vcol1, vcol2 = st.columns(2)
            with vcol1:
                if result["fwd_maps"]:
                    st.markdown(
                        f'<span class="primer-badge primer-fwd-verification">FWD: {fwd_clean}</span><br>'
                        f'<span class="status-pill status-success">● sense strand, pos {result["fwd_sense_pos"]}</span>',
                        unsafe_allow_html=True,
                    )
                elif result["fwd_antisense_pos"] is not None:
                    st.markdown(
                        f'<span class="primer-badge primer-fwd-verification">FWD: {fwd_clean}</span><br>'
                        '<span class="status-pill status-error">● Found only on antisense strand; forward primers must match the CDS sense strand</span>',
                        unsafe_allow_html=True,
                    )
                else:
                    st.markdown(
                        f'<span class="primer-badge primer-fwd-verification">FWD: {fwd_clean}</span><br>'
                        f'<span class="status-pill status-error">● Not found in CDS</span>',
                        unsafe_allow_html=True,
                    )
            with vcol2:
                if result["rev_maps"]:
                    st.markdown(
                        f'<span class="primer-badge primer-rev">REV: {rev_clean}</span><br>'
                        f'<span class="status-pill status-success">● antisense strand, pos {result["reverse_antisense_binding_pos"]}</span>',
                        unsafe_allow_html=True,
                    )
                elif result["rev_sense_pos"] is not None:
                    st.markdown(
                        f'<span class="primer-badge primer-rev">REV: {rev_clean}</span><br>'
                        '<span class="status-pill status-error">● Found only on the CDS sense strand; reverse primers must bind the antisense strand</span>',
                        unsafe_allow_html=True,
                    )
                else:
                    st.markdown(
                        f'<span class="primer-badge primer-rev">REV: {rev_clean}</span><br>'
                        f'<span class="status-pill status-error">● Not found in CDS</span>',
                        unsafe_allow_html=True,
                    )

            rev_comp = reverse_complement(rev_clean)
            st.markdown(
                f'<span class="primer-badge primer-rc">Reverse Complement: {rev_comp}</span>',
                unsafe_allow_html=True,
            )

            if result["both_map"]:
                st.markdown(
                    '<span class="primer-badge primer-warning">REPLACE THIS WARNING TEXT WITH YOUR MESSAGE</span>',
                    unsafe_allow_html=True,
                )
                st.success("Both primers verified against the CDS!")

        st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

        # ── STEP 4 — Generate Report ─────────────────────────────────

        st.markdown(
            '<div class="card"><div class="card-title">'
            '<span class="step-indicator">4</span> Generate Formatted Report'
            "</div></div>",
            unsafe_allow_html=True,
        )

        generate_btn = st.button("Generate Report", use_container_width=True)
        if generate_btn:
            if not fwd_input or not rev_input:
                st.error("Please enter both primers before generating a report.")
            else:
                fwd_clean = filter_dna(fwd_input.strip()).upper()
                rev_clean = filter_dna(rev_input.strip()).upper()
                verification = verify_primer_pair(fwd_clean, rev_clean, filtered)
                if not verification["both_map"]:
                    st.error("Both primers must be verified against the CDS before generating a report.")
                else:
                    report_text, report_html = _build_report(
                        gene=st.session_state.gene,
                        cds=cds,
                        fwd=fwd_clean,
                        rev=rev_clean,
                        reference=ref_input.strip(),
                        organism_label=organism,
                    )
                    st.session_state.final_report = report_text
                    st.session_state.final_report_html = report_html

        if st.session_state.final_report:
            st.markdown(
                f'<div class="report-box">{st.session_state.final_report_html}</div>',
                unsafe_allow_html=True,
            )
            report_basename = f"{st.session_state.gene.upper()}_primer_report"
            pdf_bytes = build_report_pdf(st.session_state.final_report)
            docx_bytes = build_report_docx(st.session_state.final_report)
            dcol1, dcol2, dcol3 = st.columns(3)
            with dcol1:
                st.download_button(
                    label="Download TXT Report",
                    data=st.session_state.final_report,
                    file_name=f"{report_basename}.txt",
                    mime="text/plain",
                    use_container_width=True,
                )
            with dcol2:
                st.download_button(
                    label="Download PDF Report",
                    data=pdf_bytes,
                    file_name=f"{report_basename}.pdf",
                    mime="application/pdf",
                    use_container_width=True,
                )
            with dcol3:
                st.download_button(
                    label="Download Word Report",
                    data=docx_bytes,
                    file_name=f"{report_basename}.docx",
                    mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                    use_container_width=True,
                )

    elif st.session_state.gene and st.session_state.input_type != "gene_name":
        st.warning("No CDS found. Check the input and try again.")


# ═══════════════════════════════════════════════════════════════════
# TAB 2 — FILTER DNA UTILITY
# ═══════════════════════════════════════════════════════════════════

with tab_filter:
    st.markdown(
        '<div class="card"><div class="card-title">'
        "Filter DNA Utility"
        "</div></div>",
        unsafe_allow_html=True,
    )

    st.caption(
        "Paste raw, messy sequence data and get a clean, continuous DNA string. "
        "Removes all spaces, line breaks, numbers, headers, and non-DNA characters. "
        "Equivalent to [bioinformatics.org/sms2/filter_dna.html](http://www.bioinformatics.org/sms2/filter_dna.html)."
    )

    filter_input = st.text_area(
        "Paste raw sequence",
        height=200,
        placeholder=(
            ">gi|12345|ref|NM_005269.3| Homo sapiens GLI1\n"
            "   1 atgttcaact cgatgacccc accaccaatc agtagctatg\n"
            "  61 gcgagccctg ctgtctccgg ..."
        ),
        key="filter_dna_input",
    )

    fcol1, fcol2 = st.columns([1, 1])
    with fcol1:
        case_option = st.radio("Output case", ["lowercase", "UPPERCASE", "As-is"], horizontal=True)
    with fcol2:
        line_width_filter = st.number_input("Characters per line (0 = no wrapping)", min_value=0, value=60, step=10)

    if filter_input:
        cleaned = filter_dna(filter_input)
        if case_option == "lowercase":
            cleaned = cleaned.lower()
        elif case_option == "UPPERCASE":
            cleaned = cleaned.upper()

        base_count = len(cleaned)

        # Wrap if requested
        if line_width_filter > 0:
            display_seq = "\n".join(
                cleaned[i : i + line_width_filter]
                for i in range(0, len(cleaned), line_width_filter)
            )
        else:
            display_seq = cleaned

        st.markdown(
            f'<div class="util-result">'
            f'<div class="util-result-label label-filtered">Filtered DNA — {base_count:,} bases</div>'
            f'<div class="util-result-seq">{display_seq}</div>'
            f'</div>',
            unsafe_allow_html=True,
        )

        # Base composition
        comp_upper = cleaned.upper()
        a_count = comp_upper.count("A")
        t_count = comp_upper.count("T")
        g_count = comp_upper.count("G")
        c_count = comp_upper.count("C")
        gc_pct = (g_count + c_count) / max(base_count, 1) * 100

        st.markdown(
            f"**Base composition:** "
            f"A={a_count:,} · T={t_count:,} · G={g_count:,} · C={c_count:,} · "
            f"GC%={gc_pct:.1f}%"
        )

        st.download_button(
            label="Download Filtered Sequence",
            data=display_seq,
            file_name="filtered_dna.txt",
            mime="text/plain",
            use_container_width=True,
        )

    else:
        st.info("Paste a sequence above to filter it.")


# ═══════════════════════════════════════════════════════════════════
# TAB 3 — SEQUENCE CALCULATOR
# ═══════════════════════════════════════════════════════════════════

with tab_calc:
    st.markdown(
        '<div class="card"><div class="card-title">'
        "Sequence Calculator"
        "</div></div>",
        unsafe_allow_html=True,
    )

    st.caption(
        "Paste a DNA sequence and get three transformations: "
        "**Reverse**, **Complement**, and **Reverse Complement**."
    )

    calc_input = st.text_area(
        "Paste DNA sequence",
        height=150,
        placeholder="ATGTTCAACTCGATGACCCCACCACCAATC...",
        key="calc_input",
    )

    calc_case = st.radio(
        "Output case", ["Preserve original", "lowercase", "UPPERCASE"],
        horizontal=True, key="calc_case",
    )

    if calc_input:
        clean_seq = filter_dna(calc_input)
        if not clean_seq:
            st.error("No valid DNA characters found in the input.")
        else:
            seq_reverse = reverse_only(clean_seq)
            seq_complement = complement_only(clean_seq)
            seq_rc = reverse_complement(clean_seq)

            if calc_case == "lowercase":
                seq_reverse = seq_reverse.lower()
                seq_complement = seq_complement.lower()
                seq_rc = seq_rc.lower()
                clean_seq = clean_seq.lower()
            elif calc_case == "UPPERCASE":
                seq_reverse = seq_reverse.upper()
                seq_complement = seq_complement.upper()
                seq_rc = seq_rc.upper()
                clean_seq = clean_seq.upper()

            # Input echo
            st.markdown(
                f'<div class="util-result">'
                f'<div class="util-result-label" style="color:#d4d4dc;">Input Sequence — {len(clean_seq):,} bases</div>'
                f'<div class="util-result-seq">{clean_seq}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Reverse
            st.markdown(
                f'<div class="util-result">'
                f'<div class="util-result-label label-reverse">Reverse</div>'
                f'<div class="util-result-seq">{seq_reverse}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Complement
            st.markdown(
                f'<div class="util-result">'
                f'<div class="util-result-label label-complement">Complement</div>'
                f'<div class="util-result-seq">{seq_complement}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Reverse Complement
            st.markdown(
                f'<div class="util-result">'
                f'<div class="util-result-label label-rc">Reverse Complement</div>'
                f'<div class="util-result-seq">{seq_rc}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Download all
            all_output = (
                f"Input ({len(clean_seq)} bases):\n{clean_seq}\n\n"
                f"Reverse:\n{seq_reverse}\n\n"
                f"Complement:\n{seq_complement}\n\n"
                f"Reverse Complement:\n{seq_rc}\n"
            )
            st.download_button(
                label="Download All Results",
                data=all_output,
                file_name="sequence_calculations.txt",
                mime="text/plain",
                use_container_width=True,
            )

    else:
        st.info("Paste a sequence above to calculate transformations.")
