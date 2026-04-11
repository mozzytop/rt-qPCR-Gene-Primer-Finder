[README.md](https://github.com/user-attachments/files/26578760/README.md)
# 🧬 RT-qPCR Gene Primer Finder

A Streamlit-powered bioinformatics automation tool that replaces a manual
data-gathering workflow for gene sequence retrieval, PCR primer matching, and
formatted report generation.

## Features

| Tab | What it does |
|-----|--------------|
| **Primer Finder** | Full workflow: NCBI gene lookup → variant selection → PMC primer search → CDS verification → formatted report |
| **Filter DNA** | Standalone utility: strips non-DNA characters from messy sequences, shows base composition + GC% |
| **Sequence Calculator** | Standalone utility: computes Reverse, Complement, and Reverse Complement |

### Primer Finder Workflow

| Step | Details |
|------|---------|
| **1 — Flexible Input** | Accepts a gene name (e.g. `GLI1`), a direct NCBI URL (e.g. `https://www.ncbi.nlm.nih.gov/nuccore/NM_005269.3`), or a raw pasted DNA sequence |
| **2 — Variant Selection** | Returns 5–10 NCBI results with accession, description, bp length, and clickable links. You choose which variant to use |
| **3 — PMC Primer Search** | Searches PubMed Central full-text for primers, extracts with direction-aware regex, verifies against the CDS |
| **4 — Manual Input** | Fallback: paste primers manually or extract from text snippets |
| **5 — Report Generation** | Produces a formatted text report with clickable DOI/PMCID reference links |

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
2. Type a gene symbol (e.g. `GLI1`, `TP53`, `BRCA1`) and click **Search**.
3. Select a transcript variant from the results table.
4. Click **Search PMC for Primers** to auto-extract and verify primers.
5. Review verified pairs or enter primers manually.
6. Click **Generate Report** to produce the formatted output.
7. Download the report as a `.txt` file.

## Project Structure

```
gene-primer-tool/
├── .github/
│   └── workflows/
│       └── keep_awake.yml      # GitHub Action: pings app every 8h
├── app.py                      # Streamlit UI (tabbed: Primer Finder | Filter DNA | Seq Calc)
├── ncbi_queries.py             # Bio.Entrez query logic (Nucleotide, PMC, PubMed)
├── utils.py                    # DNA filtering, reverse complement, primer extraction
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

## Example Images
<img width="1452" height="731" alt="image" src="https://github.com/user-attachments/assets/271acd25-a6e8-4ba0-9097-b3d72d7c54e6" />
<img width="1403" height="715" alt="image" src="https://github.com/user-attachments/assets/afcadbd5-980e-4d92-b0ea-e751156c425d" />
<img width="1424" height="752" alt="image" src="https://github.com/user-attachments/assets/83114e62-ca9a-4678-ab25-07032ac67c39" />
<img width="1469" height="793" alt="image" src="https://github.com/user-attachments/assets/917389a5-9e95-4b7c-bdea-28c0655a0e04" />
<img width="1465" height="748" alt="image" src="https://github.com/user-attachments/assets/1ab07136-17e4-4936-9b80-23169839f738" />
<img width="1468" height="795" alt="image" src="https://github.com/user-attachments/assets/6bc47ae6-a916-493e-a1ab-2e836a10f867" />


## Notes on Primer Extraction

Fully automated extraction of primer sequences from published PDFs is inherently
brittle — primers may be in supplementary tables, figures, or non-standard
notation. The tool uses a multi-pronged approach:

1. **PMC Full-Text**: Searches PubMed Central for open-access papers, extracts
   primers from Methods/Supplementary sections using direction-aware regex.
2. **Direction Inference**: Scans surrounding text for "forward", "reverse",
   "sense", "antisense" keywords to label primers automatically.
3. **CDS Verification**: All candidate pairs are programmatically tested against
   the CDS sequence before being accepted.
4. **Manual Fallback**: A text-paste box for regex extraction from any snippet,
   or direct primer entry.
