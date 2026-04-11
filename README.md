---
title: rt-qPCR Gene Primer Finder
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 8501
pinned: false
license: gpl-3.0
---

# rt-qPCR Gene Primer Finder

A Streamlit-powered bioinformatics tool that automates gene sequence retrieval, PCR primer matching, and formatted report generation — deployed via Docker on Hugging Face Spaces.

## Deployment: Hugging Face Spaces

This repo is configured to run as a **Docker**-based Space on Hugging Face. The Streamlit SDK was deprecated on HF Spaces, so it now uses a Dockerfile that installs and runs Streamlit directly.

### Setup on Hugging Face
1. Create a [Hugging Face account](https://huggingface.co/join).
2. Navigate to **Spaces** and click **Create New Space**.
3. Name your Space (e.g., `PCR-Primer-Finder`).
4. Select **Docker** as the SDK.
5. Choose your preferred visibility (Public or Private).
6. **Create Space**.

### Syncing the Repository
You can deploy the code by either manually uploading the files or by setting up a GitHub Action to sync your repo automatically.

#### Option A: Automatic Sync (Recommended)
Updates Hugging Face Space every time you push to the `HuggingFace` branch on GitHub.

1. **Generate a Hugging Face Token:**
   - Go to your [Hugging Face Settings > Access Tokens](https://huggingface.co/settings/tokens).
   - Create a new token with **Write** access.
2. **Add the Token to GitHub Secrets:**
   - In this GitHub repository, go to **Settings > Secrets and variables > Actions**.
   - Create a **New repository secret** named `HF_TOKEN`.
   - Paste your Hugging Face token here.
3. **Configure the Workflow:**
   - The workflow file at `.github/workflows/sync_to_hf.yml` handles this automatically.

### Dependencies
This app requires several Python libraries. The `requirements.txt` file contains all necessary dependencies:

```text
streamlit>=1.35.0
biopython>=1.83
pandas
requests
```

Check out the configuration reference at https://huggingface.co/docs/hub/spaces-config-reference
