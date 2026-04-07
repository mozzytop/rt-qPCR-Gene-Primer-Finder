"""
Gene Primer Lookup Tool – Streamlit Application (v2)
=====================================================

Upgraded workflow:
  Step 1 — NCBI gene lookup + CDS extraction
  Step 2 — PMC full-text primer search → auto-extract → verify → display
  Step 3 — Manual primer input fallback + verification
  Step 4 — Formatted report generation + download

Run with:
    python3 -m streamlit run app.py
"""

from __future__ import annotations

import streamlit as st
from Bio import Entrez

from utils import (
    filter_dna,
    reverse_complement,
    format_origin_block,
    format_filtered_dna,
    find_primer_on_either_strand,
    extract_primers_from_text,
    extract_primers_with_direction,
    verify_primer_pair,
)
from ncbi_queries import (
    lookup_gene,
    search_pmc_for_primers,
    search_pubmed_for_primers,
    extract_and_verify_primers,
    CDSResult,
    VerifiedPrimerPair,
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

html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

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

.verified-pair-box {
    background: linear-gradient(135deg, rgba(34,197,94,0.06), rgba(108,99,255,0.06));
    border: 1px solid rgba(34, 197, 94, 0.25);
    border-radius: 14px;
    padding: 1.3rem 1.5rem;
    margin: 0.8rem 0;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.82rem;
    line-height: 2;
}
.verified-pair-box .label {
    color: #a1a1aa;
    font-family: 'Inter', sans-serif;
    font-size: 0.72rem;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    font-weight: 600;
}

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
    "Automated NCBI gene lookup · CDS extraction · PMC full-text primer search · Reverse complement · Formatted report"
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
        "**How it works:**\n\n"
        "1. Enter a gene symbol and search NCBI\n"
        "2. The tool searches PMC full-text for primers\n"
        "3. Primers are auto-verified against the CDS\n"
        "4. If auto-search fails, paste primers manually\n"
        "5. Generate a formatted report"
    )

# ── Session state init ──────────────────────────────────────────────

for key, default in {
    "cds": None,
    "gene": "",
    "summaries": [],
    "pmc_articles": [],
    "verified_pairs": [],
    "final_report": "",
    "fwd_primer": "",
    "rev_primer": "",
    "reference": "",
    "auto_search_done": False,
}.items():
    if key not in st.session_state:
        st.session_state[key] = default

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

    trans = cds.translation
    lines.append(f'/translation="{trans[0:60]}')
    for i in range(60, len(trans), 60):
        lines.append(trans[i : i + 60])
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
    # Reset downstream state on new search
    st.session_state.pmc_articles = []
    st.session_state.verified_pairs = []
    st.session_state.final_report = ""
    st.session_state.fwd_primer = ""
    st.session_state.rev_primer = ""
    st.session_state.reference = ""
    st.session_state.auto_search_done = False
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
        st.markdown(f'<div class="sequence-box">{format_origin_block(filtered)}</div>', unsafe_allow_html=True)

    with st.expander("🧪 Amino-Acid Translation", expanded=False):
        st.markdown(f'<div class="sequence-box">{cds.translation}</div>', unsafe_allow_html=True)

    with st.expander("🔬 Filtered DNA (continuous)", expanded=False):
        st.markdown(f'<div class="sequence-box">{format_filtered_dna(filtered)}</div>', unsafe_allow_html=True)

    st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

    # ═════════════════════════════════════════════════════════════════
    # STEP 2 — PMC Full-Text Primer Search + Auto-Verification
    # ═════════════════════════════════════════════════════════════════

    st.markdown(
        '<div class="card"><div class="card-title">'
        '<span class="step-indicator">2</span> PMC Full-Text Primer Search &amp; Auto-Verification'
        "</div></div>",
        unsafe_allow_html=True,
    )

    st.caption(
        "Searches PubMed Central for open-access papers with PCR primers. "
        "Extracts primer sequences from full-text (Methods/Supplementary), "
        "detects Forward/Reverse labels, and verifies them against the CDS."
    )

    pmc_btn = st.button("🔎  Search PMC for Primers", use_container_width=True)

    if pmc_btn:
        _set_entrez()
        org_short = "human" if "homo" in organism.lower() else "mouse"
        gene_name = st.session_state.gene

        with st.spinner("Searching PubMed Central for full-text articles…"):
            articles = search_pmc_for_primers(gene_name, org_short)
            st.session_state.pmc_articles = articles

        if not articles:
            st.warning("No PMC articles found. Try the manual primer input in Step 3.")
        else:
            st.info(f"Found **{len(articles)}** PMC articles. Scanning full text for primers…")

            all_verified: list[VerifiedPrimerPair] = []
            progress = st.progress(0, text="Extracting & verifying primers…")

            for i, art in enumerate(articles):
                progress.progress(
                    (i + 1) / len(articles),
                    text=f"Scanning article {i + 1}/{len(articles)}: {art['title'][:60]}…",
                )
                try:
                    pairs = extract_and_verify_primers(art, gene_name, filtered)
                    all_verified.extend(pairs)
                except Exception:
                    continue

            progress.empty()
            st.session_state.verified_pairs = all_verified
            st.session_state.auto_search_done = True

            if all_verified:
                # Auto-populate Step 3 with the first verified pair
                best = all_verified[0]
                st.session_state.fwd_primer = best.forward
                st.session_state.rev_primer = best.reverse
                st.session_state.reference = best.citation

    # ── Display verified results ──────────────────────────────────────

    if st.session_state.verified_pairs:
        pairs = st.session_state.verified_pairs
        st.success(f"✅ Found **{len(pairs)}** verified primer pair(s) that map to the CDS!")

        org_label = "Human" if "homo" in organism.lower() else "Mouse"
        gene_upper = st.session_state.gene.upper()

        for idx, pair in enumerate(pairs):
            with st.expander(
                f"✅ Verified Pair {idx + 1}  —  {pair.source_pmcid}  |  {pair.title[:70]}…"
                if len(pair.title) > 70
                else f"✅ Verified Pair {idx + 1}  —  {pair.source_pmcid}  |  {pair.title}",
                expanded=(idx == 0),
            ):
                st.markdown(
                    f"""<div class="verified-pair-box">
<span class="label">Forward Primer</span><br>
<span class="primer-badge primer-fwd">{org_label} {gene_upper} Primer Sequence: Forward : 5'-{pair.forward}-3'</span>
<span class="status-pill status-success">● CDS {pair.fwd_strand} strand, pos {pair.fwd_position}</span><br><br>
<span class="label">Reverse Primer</span><br>
<span class="primer-badge primer-rev">{org_label} {gene_upper} Primer Sequence: Reverse : 5'-{pair.reverse}-3'</span>
<span class="status-pill status-success">● CDS {pair.rev_strand} strand, pos {pair.rev_position}</span><br><br>
<span class="label">Reverse Complement</span><br>
<span class="primer-badge primer-rc">{org_label} {gene_upper} Primer Sequence: Reverse Comp. : 5'-{pair.reverse_comp}-3'</span>
</div>""",
                    unsafe_allow_html=True,
                )
                st.markdown(f"**Reference:** {pair.citation}")

                if st.button(f"⬆️ Use this pair for report", key=f"use_pair_{idx}"):
                    st.session_state.fwd_primer = pair.forward
                    st.session_state.rev_primer = pair.reverse
                    st.session_state.reference = pair.citation
                    st.rerun()

    elif st.session_state.auto_search_done and not st.session_state.verified_pairs:
        st.warning(
            "⚠️ No verified primer pairs found in PMC full text. "
            "This can happen when primers are in supplementary PDFs or non-standard formatting. "
            "Use **Step 3** below to paste primers manually."
        )

    # Show individual articles for inspection
    if st.session_state.pmc_articles:
        with st.expander("📚 All PMC articles scanned", expanded=False):
            for art in st.session_state.pmc_articles:
                st.markdown(f"- **{art.get('pmcid', '')}** | {art['title']}")

    st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

    # ═════════════════════════════════════════════════════════════════
    # STEP 3 — Manual Primer Input & Verification (fallback)
    # ═════════════════════════════════════════════════════════════════

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

    # Optional: paste raw text to auto-extract
    with st.expander("📋 Or paste text containing 5'-…-3' primer notation", expanded=False):
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
                    pos, strand = find_primer_on_either_strand(ep.sequence, filtered)
                    dir_label = f" [{ep.direction}]" if ep.direction else ""
                    tag = (
                        f'<span class="status-pill status-success">● Maps ({strand}, pos {pos})</span>'
                        if pos is not None
                        else '<span class="status-pill status-warning">● No CDS match</span>'
                    )
                    st.markdown(f'`{ep.sequence}`{dir_label} {tag}', unsafe_allow_html=True)

                # Auto-assign if we got a clear fwd+rev
                fwd_candidates = [ep for ep in found_directed if ep.direction == "forward"]
                rev_candidates = [ep for ep in found_directed if ep.direction == "reverse"]
                if fwd_candidates and rev_candidates:
                    result = verify_primer_pair(fwd_candidates[0].sequence, rev_candidates[0].sequence, filtered)
                    if result["both_map"]:
                        st.success("✅ Verified pair found in pasted text!")
                        if st.button("⬆️ Use this pair"):
                            st.session_state.fwd_primer = result["forward_seq"]
                            st.session_state.rev_primer = result["reverse_seq"]
                            st.rerun()
            else:
                st.warning("No primer patterns found in the pasted text.")

    # Verify manually entered primers
    if fwd_input and rev_input:
        fwd_clean = filter_dna(fwd_input.strip()).upper()
        rev_clean = filter_dna(rev_input.strip()).upper()

        st.markdown("#### Verification")
        result = verify_primer_pair(fwd_clean, rev_clean, filtered)

        vcol1, vcol2 = st.columns(2)
        with vcol1:
            if result["fwd_maps"]:
                fwd_pos, fwd_strand = find_primer_on_either_strand(fwd_clean, filtered)
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
            if result["rev_maps"]:
                rev_pos, rev_strand = find_primer_on_either_strand(rev_clean, filtered)
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

        if result["both_map"]:
            st.success("✅ Both primers verified against the CDS!")
        elif result["fwd_maps"] or result["rev_maps"]:
            st.warning("⚠️ Only one primer maps to the CDS. The other may bind outside the CDS region.")

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
            fwd_clean = filter_dna(fwd_input.strip()).upper()
            rev_clean = filter_dna(rev_input.strip()).upper()
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
