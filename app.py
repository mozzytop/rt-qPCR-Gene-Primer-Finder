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
import streamlit as st
from Bio import Entrez

from utils import (
    filter_dna,
    reverse_complement,
    complement_only,
    reverse_only,
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
    search_pubmed_for_primers,
    extract_and_verify_primers,
    detect_species,
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
.primer-rev { background: rgba(239,68,68,0.12); border: 1px solid rgba(239,68,68,0.3); color: #ef4444; }
.primer-rc  { background: rgba(108,99,255,0.12); border: 1px solid rgba(108,99,255,0.3); color: #6c63ff; }

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
    line-height: 1.65; color: #d4d4dc; overflow-x: auto; white-space: pre-wrap; word-break: break-all;
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

st.markdown('<p class="main-header">🧬 Gene Primer Lookup Tool</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">'
    "Automated NCBI gene lookup · CDS extraction · PMC primer search · Sequence utilities"
    "</p>",
    unsafe_allow_html=True,
)

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
        "1. Enter a gene name, NCBI URL, or sequence\n"
        "2. Select a transcript variant\n"
        "3. Auto-search PMC for primer pairs\n"
        "4. Verify primers & generate report\n\n"
        "**Utilities:**\n"
        "- 🧹 Filter DNA — clean raw sequences\n"
        "- 🔄 Sequence Calculator — RC, complement"
    )

# ── Session state init ──────────────────────────────────────────────

for key, default in {
    "cds": None,
    "gene": "",
    "input_type": "gene_name",
    "summaries": [],
    "selected_accession": "",
    "detected_organism": "Homo sapiens",
    "detected_taxid": "9606",
    "detected_species_label": "human",
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


def _build_report(
    gene: str, cds: CDSResult, fwd: str, rev: str, reference: str, organism_label: str,
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
# TABS
# ═══════════════════════════════════════════════════════════════════

tab_main, tab_filter, tab_calc = st.tabs([
    "🧬 Primer Finder",
    "🧹 Filter DNA",
    "🔄 Sequence Calculator",
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

    col_input, col_btn = st.columns([4, 1])
    with col_input:
        user_input = st.text_area(
            "Input",
            height=80,
            placeholder="GLI1  or  https://www.ncbi.nlm.nih.gov/nuccore/NM_005269.3  or  ATGTTCAACT...",
            label_visibility="collapsed",
        )
    with col_btn:
        search_btn = st.button("🔍 Search", use_container_width=True)

    if search_btn and user_input.strip():
        _set_entrez()
        inp = user_input.strip()
        input_type = _detect_input_type(inp)

        # Reset downstream state
        for k in ["pmc_articles", "verified_pairs", "final_report", "fwd_primer",
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

            # Auto-detect species from input text
            det_org, det_taxid, det_label = detect_species(inp)
            st.session_state.detected_organism = det_org
            st.session_state.detected_taxid = det_taxid
            st.session_state.detected_species_label = det_label

            # Use sidebar organism override if it disagrees with auto-detect
            # (sidebar is the explicit user choice, takes precedence)
            if "mus" in organism.lower():
                st.session_state.detected_organism = "Mus musculus"
                st.session_state.detected_taxid = "10090"
                st.session_state.detected_species_label = "mouse"
            elif "homo" in organism.lower():
                st.session_state.detected_organism = "Homo sapiens"
                st.session_state.detected_taxid = "9606"
                st.session_state.detected_species_label = "human"

            with st.spinner(f"Querying NCBI Nucleotide (txid{st.session_state.detected_taxid})…"):
                summaries = search_nucleotide(
                    inp,
                    organism=st.session_state.detected_organism,
                    taxonomy_id=st.session_state.detected_taxid,
                )
                st.session_state.summaries = summaries

    # ── Variant selection table ──────────────────────────────────────

    if (st.session_state.input_type == "gene_name"
            and st.session_state.summaries
            and st.session_state.cds is None):

        results = st.session_state.summaries[:10]  # show top 10
        det_label = st.session_state.detected_species_label
        det_taxid = st.session_state.detected_taxid
        species_icon = "🐭" if det_label == "mouse" else "🧬"
        st.markdown(
            f'{species_icon} Species filter: '
            f'<span class="status-pill status-success">● {st.session_state.detected_organism} (txid{det_taxid})</span>'
            f' &nbsp; **{len(results)}** results found — select a variant:',
            unsafe_allow_html=True,
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
            table_html += f'<td><a href="{r["link"]}" target="_blank">🔗 View</a></td>'
            table_html += "</tr>"
        table_html += "</tbody></table>"
        st.markdown(table_html, unsafe_allow_html=True)

        # Selection UI
        options = [f'{r["accession"]}  —  {r["title"][:80]}  ({r["length"]:,} bp)' for r in results]
        selected = st.selectbox("Choose a variant", options, index=0)
        select_btn = st.button("✅ Use Selected Variant", use_container_width=True)

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
            st.markdown(f"[🔗 NCBI Link]({link})")
        st.markdown(f"CDS span: **{cds.cds_start}..{cds.cds_end}** &nbsp;|&nbsp; Filtered bases: **{len(filtered)}**")

        with st.expander("📄 ORIGIN (GenBank-formatted CDS)", expanded=False):
            st.markdown(f'<div class="sequence-box">{format_origin_block(filtered)}</div>', unsafe_allow_html=True)

        if cds.translation:
            with st.expander("🧪 Amino-Acid Translation", expanded=False):
                st.markdown(f'<div class="sequence-box">{cds.translation}</div>', unsafe_allow_html=True)

        with st.expander("🔬 Filtered DNA (continuous)", expanded=False):
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

        pmc_btn = st.button("🔎  Search PMC for Primers", use_container_width=True)

        if pmc_btn:
            _set_entrez()
            gene_name = st.session_state.gene
            taxid = st.session_state.detected_taxid

            with st.spinner("Searching PubMed Central…"):
                articles = search_pmc_for_primers(gene_name, taxonomy_id=taxid)
                st.session_state.pmc_articles = articles

            if not articles:
                st.warning("No PMC articles found. Try manual primer input in Step 3.")
            else:
                st.info(f"Found **{len(articles)}** PMC articles. Scanning for primers…")

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
                    best = all_verified[0]
                    st.session_state.fwd_primer = best.forward
                    st.session_state.rev_primer = best.reverse
                    st.session_state.reference = best.citation

        # ── Display verified results ─────────────────────────────────

        if st.session_state.verified_pairs:
            pairs = st.session_state.verified_pairs
            st.success(f"✅ Found **{len(pairs)}** verified primer pair(s)!")

            org_label = "Human" if "homo" in organism.lower() else "Mouse"
            gene_upper = st.session_state.gene.upper()

            for idx, pair in enumerate(pairs):
                exp_title = f"✅ Pair {idx + 1}  —  {pair.source_pmcid}  |  {pair.title[:70]}"
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

                    if st.button(f"⬆️ Use this pair for report", key=f"use_pair_{idx}"):
                        st.session_state.fwd_primer = pair.forward
                        st.session_state.rev_primer = pair.reverse
                        st.session_state.reference = pair.citation
                        st.rerun()

        elif st.session_state.auto_search_done and not st.session_state.verified_pairs:
            st.warning(
                "⚠️ No verified primer pairs found in PMC full text. "
                "Use **Step 3** below to paste primers manually."
            )

        if st.session_state.pmc_articles:
            with st.expander("📚 All PMC articles scanned", expanded=False):
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

                    fwd_c = [ep for ep in found_directed if ep.direction == "forward"]
                    rev_c = [ep for ep in found_directed if ep.direction == "reverse"]
                    if fwd_c and rev_c:
                        result = verify_primer_pair(fwd_c[0].sequence, rev_c[0].sequence, filtered)
                        if result["both_map"]:
                            st.success("✅ Verified pair found in pasted text!")
                            if st.button("⬆️ Use this pair"):
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
                    fwd_pos, fwd_strand = find_primer_on_either_strand(fwd_clean, filtered)
                    st.markdown(
                        f'<span class="primer-badge primer-fwd">FWD: {fwd_clean}</span><br>'
                        f'<span class="status-pill status-success">● {fwd_strand} strand, pos {fwd_pos}</span>',
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
                        f'<span class="status-pill status-success">● {rev_strand} strand, pos {rev_pos}</span>',
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

        st.markdown('<div class="gradient-divider"></div>', unsafe_allow_html=True)

        # ── STEP 4 — Generate Report ─────────────────────────────────

        st.markdown(
            '<div class="card"><div class="card-title">'
            '<span class="step-indicator">4</span> Generate Formatted Report'
            "</div></div>",
            unsafe_allow_html=True,
        )

        generate_btn = st.button("📝  Generate Report", use_container_width=True)
        if generate_btn:
            if not fwd_input or not rev_input:
                st.error("Please enter both primers before generating a report.")
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

    elif st.session_state.gene and st.session_state.input_type != "gene_name":
        st.warning("No CDS found. Check the input and try again.")


# ═══════════════════════════════════════════════════════════════════
# TAB 2 — FILTER DNA UTILITY
# ═══════════════════════════════════════════════════════════════════

with tab_filter:
    st.markdown(
        '<div class="card"><div class="card-title">'
        "🧹 Filter DNA Utility"
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
            label="⬇️  Download Filtered Sequence",
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
        "🔄 Sequence Calculator"
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
                f'<div class="util-result-label label-reverse">↩️ Reverse</div>'
                f'<div class="util-result-seq">{seq_reverse}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Complement
            st.markdown(
                f'<div class="util-result">'
                f'<div class="util-result-label label-complement">🔀 Complement</div>'
                f'<div class="util-result-seq">{seq_complement}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Reverse Complement
            st.markdown(
                f'<div class="util-result">'
                f'<div class="util-result-label label-rc">🔄 Reverse Complement</div>'
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
                label="⬇️  Download All Results",
                data=all_output,
                file_name="sequence_calculations.txt",
                mime="text/plain",
                use_container_width=True,
            )

    else:
        st.info("Paste a sequence above to calculate transformations.")
