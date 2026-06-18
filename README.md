# 🧬 Protein Sequence Analyzer

An interactive web app that analyzes a protein sequence and returns its
biochemical properties, amino-acid composition and predicted membrane topology —
with a **hybrid architecture**: the heavy numerical work runs in **C++**
(exposed to Python via **pybind11**), while the interface and the interactive
charts are built in **Python** (**Streamlit + Plotly**).

> Paste a sequence, upload a FASTA file, or fetch one directly from **UniProt**
> by accession — and get the full analysis in seconds.

## ✨ Features

- **Provenance** parsed from the FASTA header (UniProt / NCBI): protein name,
  organism, database, accession (with a link to UniProt), gene.
- **Biochemical metrics**: length, molecular weight, isoelectric point (solved
  numerically by bisection), net charge at pH 7, peptide bonds, cysteines /
  disulfide bonds, **GRAVY**, **aliphatic index**, **molar extinction
  coefficient at 280 nm** and **A₂₈₀ of a 1 g/L solution** (to estimate protein
  concentration from a spectrophotometer reading).
- **Titration curve** (net charge vs pH) — the zero crossing matches the pI.
- **Amino-acid composition**: frequency bar chart and physico-chemical class
  pie chart, with an interactive hover-highlight effect.
- **Membrane domains**: Kyte-Doolittle hydrophobicity profile with automatic
  detection of transmembrane segments, highlighted on the sequence.
- **Downloadable report** (self-contained HTML) and **CSV** export.
- **6 languages**: English, Italiano, Español, Deutsch, Français, 中文.

## 🧱 Architecture

| Layer | Tech | Role |
|-------|------|------|
| Compute core | **C++** (`seq_core.cpp`) via **pybind11** | heavy math: pI by bisection, titration curve, hydrophobicity profile, transmembrane detection |
| Frontend | **Python** — **Streamlit** + **Plotly** | UI, interactive charts, FASTA/UniProt input, report generation |

The C++ module is compiled to a native Python extension (`seq_core`) with
`setup.py build_ext`.

## 🚀 Run locally

```bash
# 1. install dependencies
pip install -r requirements.txt

# 2. compile the C++ core
python setup.py build_ext --inplace

# 3. launch the app
streamlit run app.py
```

Then open http://localhost:8501.

## ☁️ Deploy (Streamlit Community Cloud)

The repo is ready to deploy:

- `requirements.txt` — Python dependencies (incl. `pybind11`, `setuptools`).
- `packages.txt` — system packages (`build-essential`, i.e. the C++ compiler).
- `app.py` compiles `seq_core` automatically on first start if the binary
  isn't present, so the native module is built on the server.

Point Streamlit Community Cloud at this repo with `app.py` as the main file.

## ⚠️ Note

Isoelectric point and transmembrane-domain predictions are approximate
theoretical estimates (they depend on the pKa set and the hydropathy threshold).

## 📂 Files

```
app.py          # Streamlit frontend (UI, charts, i18n, report)
seq_core.cpp    # C++ compute core (pybind11)
setup.py        # builds the C++ extension
analizza.py     # command-line version of the analysis
KCNQ1.fasta     # example sequence (the KCNQ1 potassium channel)
requirements.txt / packages.txt   # deploy config
```
