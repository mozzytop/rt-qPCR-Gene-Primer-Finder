[README.md](https://github.com/user-attachments/files/26643855/README.md)
---
title: RT-qPCR Gene Primer Finder
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: streamlit
sdk_version: 1.42.0
app_file: app.py
pinned: true
---

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

## Usage

1. Enter your NCBI email in the sidebar (or pre-configured for lab members).
2. Type a gene symbol (e.g. `GLI1`, `TP53`, `BRCA1`) and click **Search**.
3. Select a transcript variant from the results table.
4. Click **Search PMC for Primers** to auto-extract and verify primers.
5. Review verified pairs or enter primers manually.
6. Click **Generate Report** to produce the formatted output.
7. Download the report as `.txt`, `.pdf`, or `.docx`.
