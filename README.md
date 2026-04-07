# 🧬 Gene Primer Lookup Tool

A Streamlit-powered bioinformatics automation tool that replaces a manual
data-gathering workflow for gene sequence retrieval, PCR primer matching, and
formatted report generation.

## Features

| Step | What it does |
|------|-------------|
| **1 — Gene Lookup** | Queries NCBI Nucleotide for a gene symbol, auto-selects the best "transcript variant 1, mRNA" record, and extracts the CDS origin sequence + amino-acid translation. |
| **2 — PubMed Search** | Searches PubMed for papers containing PCR primers for the gene. Attempts automatic regex extraction of `5'-…-3'` primer sequences from abstracts and PMC full text. |
| **3 — Primer Verification** | Lets you enter (or paste) forward/reverse primers and verifies they map to the CDS on either strand. Computes the reverse complement of the reverse primer. |
| **4 — Report Generation** | Produces a strictly formatted text report matching the provided template, ready for download. |

## Prerequisites

- Python 3.10+
- An internet connection (for NCBI Entrez API calls)
- An email address (NCBI API requirement)
- *(Optional)* An NCBI API key to raise the rate limit from 3 → 10 requests/sec.
  Register at <https://www.ncbi.nlm.nih.gov/account/>

## Installation

```bash
cd gene-primer-tool
pip install -r requirements.txt
```

## Running the App

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

## Usage

1. Enter your NCBI email in the sidebar.
2. Type a gene symbol (e.g. `GLI1`, `TP53`, `BRCA1`) and click **Search NCBI**.
3. Click **Search PubMed for Primers** to find literature with primer sequences.
4. Enter the forward and reverse primer sequences (from literature or manually).
5. Click **Generate Report** to produce the formatted output.
6. Download the report as a `.txt` file.

## Project Structure

```
gene-primer-tool/
├── app.py              # Streamlit UI
├── ncbi_queries.py     # Bio.Entrez query logic
├── utils.py            # DNA filtering, reverse complement, formatting
├── requirements.txt    # Python dependencies
└── README.md           # This file
```

## Notes on Primer Extraction

Fully automated extraction of primer sequences from published PDFs is inherently
brittle — primers may be in supplementary tables, figures, or non-standard
notation. The tool uses a two-pronged approach:

1. **Automatic**: Regex searches abstracts and PMC Open Access full text for
   `5'-SEQUENCE-3'` patterns.
2. **Manual fallback**: A text-paste box lets you drop in any snippet of text
   (e.g. copied from a paper) for regex extraction, or you can type primer
   sequences directly.

Both approaches feed into the same verification pipeline that maps primers
against the CDS.
