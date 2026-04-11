# Deployment: Hugging Face Spaces

I built this repo to run as a **Streamlit** SDK on Hugging Face Spaces. Follow these steps to host your own version of the Gene Primer Finder.

## 1. Space Configuration (YAML Metadata)
Hugging Face requires a specific metadata block at the very top of the `README.md` file in the root directory. **Ensure this block is the first thing in your file:**

```yaml
---
title: rt-qPCR Gene Primer Finder
colorFrom: blue
colorTo: green
sdk: streamlit
app_file: app.py
pinned: false
license: mit
---
```

## 2. Setup on Hugging Face
1. Create a [Hugging Face account](https://huggingface.co/join).
2. Navigate to **Spaces** and click **Create New Space**.
3. Name your Space (e.g., `gene-primer-finder`).
4. Select **Streamlit** as the SDK.
5. Choose your preferred visibility (Public or Private).
6. **Create Space**.

## 3. Syncing the Repository
You can deploy the code by either manually uploading the files or by setting up a GitHub Action to sync your the repo automatically.

### Option A: Automatic Sync (Recommended)
Updates Hugging Face Space every time you push to the `HuggingFace` branch on GitHub.

1. **Generate a Hugging Face Token:**
   - Go to your [Hugging Face Settings > Access Tokens](https://huggingface.co/settings/tokens).
   - Create a new token with **Write** access.
2. **Add the Token to GitHub Secrets:**
   - In this GitHub repository, go to **Settings > Secrets and variables > Actions**.
   - Create a **New repository secret** named `HF_TOKEN`.
   - Paste your Hugging Face token here.
3. **Configure the Workflow:**
   - Create a folder named `.github/workflows` (if it doesn't exist).
   - Create a file inside called `sync_to_hf.yml` and paste the following:

```yaml
name: Sync to Hugging Face Hub
on:
  push:
    branches: [HuggingFace]
  workflow_dispatch:

jobs:
  sync-to-hub:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
        with:
          fetch-depth: 0
          lfs: true
      - name: Push to hub
        env:
          HF_TOKEN: ${{ secrets.HF_TOKEN }}
        run: git push --force https://YOUR_HF_USERNAME:$HF_TOKEN@huggingface.co/spaces/YOUR_HF_USERNAME/YOUR_SPACE_NAME HuggingFace:main
```
*Note: Replace `YOUR_HF_USERNAME` and `YOUR_SPACE_NAME` with your actual details.*

## 4. Dependencies
This app requires several Python libraries to function. Double check that your `requirements.txt` file is present in the root directory and contains _at least_ the following:

```text
streamlit
pandas
biopython
requests
```
