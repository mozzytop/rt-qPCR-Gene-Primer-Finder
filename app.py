"""
Gene Primer Lookup Tool – Streamlit Application
=================================================

A bioinformatics automation tool that:
1. Queries NCBI Nucleotide for a gene's mRNA / CDS annotation.
2. Extracts the CDS origin sequence and amino-acid translation.
3. Searches PubMed for PCR primer sequences in published literature.
4. Verifies primers against the CDS, computes reverse complement.
5. Produces a strictly formatted report.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import textwrap
import streamlit as st

from Bio import Entrez

from utils import (
    filter_dna,
    reverse_complement,
    format_origin_block,
    format_filtered_dna,
    format_translation,
    find_primer_on_either_strand,
    extract_primers_from_text,
)
from ncbi_queries import (
    lookup_gene,
    search_pubmed_for_primers,
    extract_primers_from_article,
    fetch_genbank_record,
    extract_cds,
    search_nucleotide,
    CDSResult,
    PrimerHit,
)

# ── Page config ──────────────────────────────────────────────────────

st.set_page_config(
    page_title="Gene Primer Lookup Tool",
    page_icon="🧬",
    layout="wide",
)

# ── Custom CSS ───────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

/* Root variables */
:root {
    --bg-primary: #0f1117;
    --bg-card: #1a1d29;
    --bg-card-hover: #1f2335;
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

html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

/* Header styling */
.main-header {
    background: linear-gradient(135deg, #6c63ff 0%, #00d2ff 50%, #7c3aed 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    font-weight: 700;
    font-size: 2.4rem;
    letter-spacing: -0.02em;
    margin-bottom: 0.1rem;
}

.sub-header {
    color: var(--text-secondary);
    font-size: 1.05rem;
    font-weight: 300;
    margin-bottom: 2rem;
}

/* Card containers */
.card {
    background: var(--bg-card);
    border: 1px solid var(--border-color);
    border-radius: 16px;
    padding: 1.5rem;
    margin-bottom: 1rem;
    transition: border-color 0.3s ease, box-shadow 0.3s ease;
}
.card:hover {
    border-color: var(--accent-1);
    box-shadow: 0 0 20px rgba(108, 99, 255, 0.08);
}

.card-title {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 0.8rem;
    display: flex;
    align-items: center;
    gap: 0.5rem;
}

/* Sequence display */
.sequence-box {
    background: #111318;
    border: 1px solid #2a2d3a;
    border-radius: 10px;
    padding: 1rem 1.2rem;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.78rem;
    line-height: 1.6;
    color: #c4c4cc;
    overflow-x: auto;
    white-space: pre;
    max-height: 400px;
    overflow-y: auto;
}

/* Primer badge */
.primer-badge {
    display: inline-block;
    padding: 0.35rem 0.85rem;
    border-radius: 8px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.82rem;
    font-weight: 500;
    letter-spacing: 0.03em;
}
.primer-fwd {
    background: rgba(34, 197, 94, 0.12);
    border: 1px solid rgba(34, 197, 94, 0.3);
    color: #22c55e;
}
.primer-rev {
    background: rgba(239, 68, 68, 0.12);
    border: 1px solid rgba(239, 68, 68, 0.3);
    color: #ef4444;
}
.primer-rc {
    background: rgba(108, 99, 255, 0.12);
    border: 1px solid rgba(108, 99, 255, 0.3);
    color: #6c63ff;
}

/* Status indicators */
.status-pill {
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    padding: 0.25rem 0.7rem;
    border-radius: 999px;
    font-size: 0.78rem;
    font-weight: 500;
}
.status-success {
    background: rgba(34, 197, 94, 0.12);
    color: #22c55e;
}
.status-warning {
    background: rgba(245, 158, 11, 0.12);
    color: #f59e0b;
}
.status-error {
    background: rgba(239, 68, 68, 0.12);
    color: #ef4444;
}

/* Step indicators */
.step-indicator {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 28px;
    height: 28px;
    border-radius: 50%;
    background: linear-gradient(135deg, var(--accent-1), var(--accent-2));
    color: white;
    font-size: 0.8rem;
    font-weight: 700;
    margin-right: 0.6rem;
}

/* Report output */
.report-box {
    background: #0c0e14;
    border: 1px solid #6c63ff44;
    border-radius: 12px;
    padding: 1.5rem;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.76rem;
    line-height: 1.65;
    color: #d4d4dc;
    overflow-x: auto;
    white-space: pre-wrap;
    word-break: break-all;
}

/* Streamlit overrides */
.stTextInput > div > div > input,
.stSelectbox > div > div {
    border-radius: 10px !important;
}
div.stButton > button {
    background: linear-gradient(135deg, #6c63ff, #7c3aed) !important;
    color: white !important;
    border: none !important;
    border-radius: 10px !important;
    padding: 0.55rem 1.8rem !important;
    font-weight: 600 !important;
    transition: all 0.3s ease !important;
}
div.stButton > button:hover {
    box-shadow: 0 4px 20px rgba(108, 99, 255, 0.35) !important;
    transform: translateY(-1px) !important;
}

.stTextArea textarea {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.8rem !important;
    border-radius: 10px !important;
}

/* Divider */
.gradient-divider {
    height: 2px;
    background: linear-gradient(90deg, transparent, #6c63ff, #00d2ff, transparent);
    border: none;
    margin: 2rem 0;
    border-radius: 1px;
}
</style>
""", unsafe_allow_html=True)

# ── Header ───────────────────────────────────────────────────────────

st.markdown('<p class="main-header">🧬 Gene Primer Lookup Tool</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">'
    "Automated NCBI gene lookup · CDS extraction · PubMed primer search · Reverse complement · Formatted report"
    "</p>",
    unsafe_allow_html=True,
)
st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

# ── Sidebar config ───────────────────────────────────────────────────

with st.sidebar:
    st.markdown("### ⚙️ Configuration")
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
        "**Tip:** If automatic primer extraction fails, you can paste "
        "primer sequences manually in Step 3."
    )

# ── Session state init ──────────────────────────────────────────────

if "cds" not in st.session_state:
    st.session_state.cds = None
if "gene" not in st.session_state:
    st.session_state.gene = ""
if "summaries" not in st.session_state:
    st.session_state.summaries = []
if "pubmed_articles" not in st.session_state:
    st.session_state.pubmed_articles = []
if "final_report" not in st.session_state:
    st.session_state.final_report = ""
if "fwd_primer" not in st.session_state:
    st.session_state.fwd_primer = ""
if "rev_primer" not in st.session_state:
    st.session_state.rev_primer = ""
if "reference" not in st.session_state:
    st.session_state.reference = ""

# ── Helpers ──────────────────────────────────────────────────────────

def _set_entrez():
    if not email:
        st.error("⚠️  Please enter your NCBI email in the sidebar.")
        st.stop()
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def _build_report(
    gene: str,
    cds: CDSResult,
    fwd: str,
    rev: str,
    reference: str,
    organism_label: str,
) -> str:
    """Build the final formatted output string."""
    filtered = filter_dna(cds.cds_dna)
    origin_block = format_origin_block(filtered)
    filtered_block = format_filtered_dna(filtered)
    rev_comp = reverse_complement(rev)
    link = f"https://www.ncbi.nlm.nih.gov/nuccore/{cds.accession}?from={cds.cds_start}&to={cds.cds_end}"

    org_label = "Human" if "homo" in organism_label.lower() else "Mouse"

    lines: list[str] = []
    lines.append(f"Gene: {gene.upper()}")
    lines.append(f"Link: {link}")
    lines.append("ORIGIN")
    lines.append(origin_block)
    lines.append("//")

    # Translation block
    trans = cds.translation
    lines.append(f'/translation="{trans[0:60]}')
    for i in range(60, len(trans), 60):
        lines.append(trans[i : i + 60])
    # close with quote on last line
    if lines[-1] and not lines[-1].endswith('"'):
        lines[-1] += '"'

    lines.append("Filter DNA results:")
    lines.append(f">filtered DNA sequence consisting of {len(filtered)} bases.")
    lines.append(filtered_block)

    lines.append(f"{org_label} {gene.upper()} Primer Sequence: Forward : 5'-{fwd}-3'")
    lines.append(f"{org_label} {gene.upper()} Primer Sequence: Reverse : 5'-{rev}-3'")
    lines.append(f"{org_label} {gene.upper()} Primer Sequence: Reverse Comp. : 5'-{rev_comp}-3'")

    if reference:
        lines.append(f"Reference: {reference}")

    return "\n".join(lines)

# ═══════════════════════════════════════════════════════════════════
# STEP 1 — Gene Lookup
# ═══════════════════════════════════════════════════════════════════

st.markdown(
    '<div class="card"><div class="card-title">'
    '<span class="step-indicator">1</span> Gene Lookup &amp; CDS Extraction'
    "</div></div>",
    unsafe_allow_html=True,
)

col1, col2 = st.columns([3, 1])
with col1:
    gene_input = st.text_input(
        "Gene Symbol",
        value=st.session_state.gene,
        placeholder="e.g. GLI1, TP53, BRCA1",
        label_visibility="collapsed",
    )
with col2:
    search_btn = st.button("🔍  Search NCBI", use_container_width=True)

if search_btn and gene_input:
    _set_entrez()
    st.session_state.gene = gene_input.strip()
    with st.spinner("Querying NCBI Nucleotide…"):
        summaries, cds = lookup_gene(gene_input.strip(), organism)
        st.session_state.summaries = summaries
        st.session_state.cds = cds

cds: CDSResult | None = st.session_state.cds

if cds is not None:
    filtered = filter_dna(cds.cds_dna)
    link = f"https://www.ncbi.nlm.nih.gov/nuccore/{cds.accession}?from={cds.cds_start}&to={cds.cds_end}"

    st.markdown(
        f'<span class="status-pill status-success">● Found</span> &nbsp; '
        f"**{cds.description}**",
        unsafe_allow_html=True,
    )
    st.markdown(f"[🔗 NCBI Link]({link})")
    st.markdown(f"CDS span: **{cds.cds_start}..{cds.cds_end}** &nbsp;|&nbsp; Filtered bases: **{len(filtered)}**")

    with st.expander("📄 ORIGIN (GenBank-formatted CDS)", expanded=False):
        st.markdown(
            f'<div class="sequence-box">{format_origin_block(filtered)}</div>',
            unsafe_allow_html=True,
        )

    with st.expander("🧪 Amino-Acid Translation", expanded=False):
        st.markdown(
            f'<div class="sequence-box">{cds.translation}</div>',
            unsafe_allow_html=True,
        )

    with st.expander("🔬 Filtered DNA (continuous)", expanded=False):
        st.markdown(
            f'<div class="sequence-box">{format_filtered_dna(filtered)}</div>',
            unsafe_allow_html=True,
        )

    st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

    # ═════════════════════════════════════════════════════════════════
    # STEP 2 — PubMed Primer Search
    # ═════════════════════════════════════════════════════════════════

    st.markdown(
        '<div class="card"><div class="card-title">'
        '<span class="step-indicator">2</span> PubMed Primer Search'
        "</div></div>",
        unsafe_allow_html=True,
    )

    pubmed_btn = st.button("🔎  Search PubMed for Primers", use_container_width=True)

    if pubmed_btn:
        _set_entrez()
        org_short = "human" if "homo" in organism.lower() else "mouse"
        with st.spinner("Searching PubMed…"):
            articles = search_pubmed_for_primers(st.session_state.gene, org_short)
            st.session_state.pubmed_articles = articles

    if st.session_state.pubmed_articles:
        st.markdown(f"Found **{len(st.session_state.pubmed_articles)}** candidate articles.")
        for i, art in enumerate(st.session_state.pubmed_articles):
            with st.expander(f"PMID {art['pmid']} — {art['title'][:90]}…" if len(art['title']) > 90 else f"PMID {art['pmid']} — {art['title']}", expanded=False):
                st.markdown(f"**Citation:** {art['citation']}")
                if art.get("pmcid"):
                    st.markdown(f"**PMCID:** {art['pmcid']}")

                # Try automatic extraction
                _set_entrez()
                from ncbi_queries import extract_primers_from_article
                primers = extract_primers_from_article(art, st.session_state.gene)
                if primers:
                    st.markdown("**Auto-detected primer sequences:**")
                    for p in primers:
                        pos, strand = find_primer_on_either_strand(p.sequence, filtered)
                        match_str = (
                            f'<span class="status-pill status-success">● Maps to CDS ({strand} strand, pos {pos})</span>'
                            if pos is not None
                            else '<span class="status-pill status-warning">● Not found in CDS</span>'
                        )
                        st.markdown(
                            f'<span class="primer-badge primer-fwd">{p.sequence}</span> {match_str}',
                            unsafe_allow_html=True,
                        )
                else:
                    st.info("No primer sequences auto-detected in abstract / full text.")

                # Use-this-article button
                if st.button(f"Use this reference", key=f"use_{i}"):
                    st.session_state.reference = art["citation"]

    st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

    # ═════════════════════════════════════════════════════════════════
    # STEP 3 — Manual Primer Input & Verification
    # ═════════════════════════════════════════════════════════════════

    st.markdown(
        '<div class="card"><div class="card-title">'
        '<span class="step-indicator">3</span> Primer Input &amp; Verification'
        "</div></div>",
        unsafe_allow_html=True,
    )

    st.caption(
        "Enter the forward and reverse primer sequences (DNA only, no 5'/3' markers). "
        "The tool will verify they map to the CDS and compute the reverse complement."
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

    # Optional: paste raw text to auto-extract
    with st.expander("📋 Or paste text containing 5'-…-3' primer notation", expanded=False):
        raw_snippet = st.text_area(
            "Paste abstract / methods text",
            height=120,
            placeholder="The tool will regex-extract primer sequences from this text…",
            key="raw_snippet",
        )
        if raw_snippet:
            found = extract_primers_from_text(raw_snippet)
            if found:
                st.markdown("**Extracted sequences:**")
                for seq in found:
                    pos, strand = find_primer_on_either_strand(seq, filtered)
                    tag = (
                        f'<span class="status-pill status-success">● Maps ({strand}, pos {pos})</span>'
                        if pos is not None
                        else '<span class="status-pill status-warning">● No CDS match</span>'
                    )
                    st.markdown(f'`{seq}` {tag}', unsafe_allow_html=True)
            else:
                st.warning("No 5'-…-3' patterns found in the pasted text.")

    # Verify primers
    if fwd_input and rev_input:
        fwd_clean = filter_dna(fwd_input.strip())
        rev_clean = filter_dna(rev_input.strip())

        st.markdown("#### Verification")
        vcol1, vcol2 = st.columns(2)
        with vcol1:
            fwd_pos, fwd_strand = find_primer_on_either_strand(fwd_clean, filtered)
            if fwd_pos is not None:
                st.markdown(
                    f'<span class="primer-badge primer-fwd">FWD: {fwd_clean}</span><br>'
                    f'<span class="status-pill status-success">● Found on {fwd_strand} strand at position {fwd_pos}</span>',
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    f'<span class="primer-badge primer-fwd">FWD: {fwd_clean}</span><br>'
                    f'<span class="status-pill status-error">● Not found in CDS</span>',
                    unsafe_allow_html=True,
                )
        with vcol2:
            rev_pos, rev_strand = find_primer_on_either_strand(rev_clean, filtered)
            if rev_pos is not None:
                st.markdown(
                    f'<span class="primer-badge primer-rev">REV: {rev_clean}</span><br>'
                    f'<span class="status-pill status-success">● Found on {rev_strand} strand at position {rev_pos}</span>',
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

    st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

    # ═════════════════════════════════════════════════════════════════
    # STEP 4 — Generate Report
    # ═════════════════════════════════════════════════════════════════

    st.markdown(
        '<div class="card"><div class="card-title">'
        '<span class="step-indicator">4</span> Generate Formatted Report'
        "</div></div>",
        unsafe_allow_html=True,
    )

    generate_btn = st.button("📝  Generate Report", use_container_width=True)
    if generate_btn:
        if not fwd_input or not rev_input:
            st.error("Please enter both forward and reverse primers before generating a report.")
        else:
            fwd_clean = filter_dna(fwd_input.strip())
            rev_clean = filter_dna(rev_input.strip())
            report = _build_report(
                gene=st.session_state.gene,
                cds=cds,
                fwd=fwd_clean,
                rev=rev_clean,
                reference=ref_input.strip(),
                organism_label=organism,
            )
            st.session_state.final_report = report

    if st.session_state.final_report:
        st.markdown(
            f'<div class="report-box">{st.session_state.final_report}</div>',
            unsafe_allow_html=True,
        )
        st.download_button(
            label="⬇️  Download Report (.txt)",
            data=st.session_state.final_report,
            file_name=f"{st.session_state.gene.upper()}_primer_report.txt",
            mime="text/plain",
            use_container_width=True,
        )

elif st.session_state.gene:
    st.warning("No results found. Try a different gene name or organism.")
