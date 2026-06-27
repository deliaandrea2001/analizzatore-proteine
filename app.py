#!/usr/bin/env python3
# =====================================================================
#  app.py  -  Web app interattiva (Streamlit) dell'analizzatore proteico.
# ---------------------------------------------------------------------
#  Incolla una sequenza (o un FASTA), premi "Analizza/Analyze" e vedi:
#  provenienza (dall'intestazione FASTA), metriche biochimiche e grafici.
#  Il calcolo gira sul core C++. Interfaccia bilingue (English/Italiano),
#  lingua di default: inglese.
#
#  Avvio:   streamlit run app.py
# =====================================================================

import os
import re
import json
import html
import sys
import subprocess
import urllib.request
import urllib.error
import urllib.parse
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components


def _carica_core():
    """Importa il core C++. Se il modulo compilato non c'è ancora (es. al primo
    avvio su un server di deploy), lo compila al volo con setup.py e riprova.
    In locale, dove il .so esiste già, importa e basta."""
    try:
        import seq_core
        return seq_core
    except ImportError:
        cartella = os.path.dirname(os.path.abspath(__file__))
        subprocess.run([sys.executable, "setup.py", "build_ext", "--inplace"],
                       cwd=cartella, check=True)
        import seq_core
        return seq_core


seq_core = _carica_core()        # <-- il nostro core C++

AA_VALIDI = set("ACDEFGHIKLMNPQRSTVWY")
NOMI_AA = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
}

# Nome scientifico e nome comune (per lingua) degli organismi UniProt comuni.
ORG_SCI = {
    "HUMAN": "Homo sapiens", "MOUSE": "Mus musculus", "RAT": "Rattus norvegicus",
    "BOVIN": "Bos taurus", "PIG": "Sus scrofa", "CHICK": "Gallus gallus",
    "YEAST": "Saccharomyces cerevisiae", "ECOLI": "Escherichia coli",
    "DROME": "Drosophila melanogaster", "CAEEL": "Caenorhabditis elegans",
    "XENLA": "Xenopus laevis", "DANRE": "Danio rerio",
}
ORG_COMUNE = {
    "HUMAN": {"it": "uomo", "en": "human"}, "MOUSE": {"it": "topo", "en": "mouse"},
    "RAT": {"it": "ratto", "en": "rat"}, "BOVIN": {"it": "bovino", "en": "cattle"},
    "PIG": {"it": "maiale", "en": "pig"}, "CHICK": {"it": "pollo", "en": "chicken"},
    "YEAST": {"it": "lievito", "en": "yeast"}, "ECOLI": {"it": "", "en": ""},
    "DROME": {"it": "moscerino", "en": "fruit fly"}, "CAEEL": {"it": "", "en": ""},
    "XENLA": {"it": "rana", "en": "frog"}, "DANRE": {"it": "zebrafish", "en": "zebrafish"},
}
DB_DESC = {
    "sp": {"it": "Swiss-Prot (UniProt, revisionata manualmente)",
           "en": "Swiss-Prot (UniProt, manually reviewed)"},
    "tr": {"it": "TrEMBL (UniProt, annotazione automatica)",
           "en": "TrEMBL (UniProt, automatically annotated)"},
}

# Etichette delle classi chimico-fisiche per lingua (chiave = nome canonico dal core C++).
CLASSI_LABEL = {
    "Idrofobici":       {"it": "Idrofobici", "en": "Hydrophobic", "es": "Hidrofóbicos",
                         "de": "Hydrophob", "fr": "Hydrophobes", "zh": "疏水性"},
    "Polari":           {"it": "Polari", "en": "Polar", "es": "Polares",
                         "de": "Polar", "fr": "Polaires", "zh": "极性"},
    "Carichi positivi": {"it": "Carichi positivi", "en": "Positively charged", "es": "Carga positiva",
                         "de": "Positiv geladen", "fr": "Chargés positivement", "zh": "带正电"},
    "Carichi negativi": {"it": "Carichi negativi", "en": "Negatively charged", "es": "Carga negativa",
                         "de": "Negativ geladen", "fr": "Chargés négativement", "zh": "带负电"},
}
# Verdetti topologici restituiti dal core C++ (in italiano) -> traduzione.
VERDETTO_LABEL = {
    "Sequenza troppo corta per l'analisi":
        {"it": "Sequenza troppo corta per l'analisi", "en": "Sequence too short for analysis",
         "es": "Secuencia demasiado corta para el análisis", "de": "Sequenz zu kurz für die Analyse",
         "fr": "Séquence trop courte pour l'analyse", "zh": "序列太短，无法分析"},
    "Proteina globulare solubile (nessun dominio di membrana)":
        {"it": "Proteina globulare solubile (nessun dominio di membrana)",
         "en": "Soluble globular protein (no membrane domain)",
         "es": "Proteína globular soluble (sin dominio de membrana)",
         "de": "Lösliches globuläres Protein (keine Membrandomäne)",
         "fr": "Protéine globulaire soluble (aucun domaine membranaire)",
         "zh": "可溶性球状蛋白（无膜结构域）"},
    "Proteina di membrana (single-pass)":
        {"it": "Proteina di membrana (single-pass)", "en": "Membrane protein (single-pass)",
         "es": "Proteína de membrana (single-pass)", "de": "Membranprotein (single-pass)",
         "fr": "Protéine membranaire (single-pass)", "zh": "膜蛋白（单次跨膜）"},
    "Proteina di membrana (multi-pass)":
        {"it": "Proteina di membrana (multi-pass)", "en": "Membrane protein (multi-pass)",
         "es": "Proteína de membrana (multi-pass)", "de": "Membranprotein (multi-pass)",
         "fr": "Protéine membranaire (multi-pass)", "zh": "膜蛋白（多次跨膜）"},
}

# Idrofobicità di Kyte-Doolittle per residuo (KD): negativa = idrofilico,
# positiva = idrofobico. Usata per colorare le barre della composizione in modo
# coerente con l'analisi dei domini di membrana (anch'essa basata su Kyte-Doolittle).
HYDRO_KD = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5,
}
# Scala di colori divergente dell'idrofobicità: teal (idrofilico) → sabbia
# (neutro) → arancio (idrofobico). Gli estremi coincidono coi colori delle classi
# "Polari" e "Idrofobici", così barre e torta condividono la stessa palette.
HYDRO_SCALE = [[0.0, "#3aa6a0"], [0.5, "#ece6d6"], [1.0, "#e07a3f"]]

# Scale per-residuo per il profilo del grafico (la finestra le media). La scala
# Kyte-Doolittle è anche quella usata dal core C++ per il rilevamento dei domini
# TM. 'verso': +1 = valori alti idrofobici, -1 = valori alti idrofilici.
HOPP_WOODS = {                            # Hopp & Woods (1981): alti = idrofilico
    'R': 3.0, 'D': 3.0, 'E': 3.0, 'K': 3.0, 'S': 0.3, 'N': 0.2, 'Q': 0.2,
    'G': 0.0, 'P': 0.0, 'T': -0.4, 'A': -0.5, 'H': -0.5, 'C': -1.0, 'M': -1.3,
    'V': -1.5, 'I': -1.8, 'L': -1.8, 'Y': -2.3, 'F': -2.5, 'W': -3.4,
}
EISENBERG = {                             # Eisenberg (1984), consenso normalizzato
    'I': 1.38, 'F': 1.19, 'V': 1.08, 'L': 1.06, 'W': 0.81, 'M': 0.64, 'A': 0.62,
    'G': 0.48, 'C': 0.29, 'Y': 0.26, 'P': 0.12, 'T': -0.05, 'S': -0.18, 'H': -0.40,
    'E': -0.74, 'N': -0.78, 'Q': -0.85, 'D': -0.90, 'K': -1.50, 'R': -2.53,
}
SCALE_IDRO = {                            # chiave -> (nome, tabella, è_idrofobicità)
    "kd": ("Kyte-Doolittle", HYDRO_KD, True),
    "hopp": ("Hopp-Woods", HOPP_WOODS, False),
    "eisenberg": ("Eisenberg", EISENBERG, True),
}

# Colori delle classi (chiave = nome canonico), armonizzati con HYDRO_SCALE.
# Servono anche al JS per attenuare le fette della torta.
COLORI_CLASSI = {
    "Idrofobici": "#e07a3f",        # arancio  (estremo idrofobico della scala)
    "Polari": "#3aa6a0",            # teal     (estremo idrofilico della scala)
    "Carichi positivi": "#3b7dd8",  # blu
    "Carichi negativi": "#d1495b",  # rosso
}

# ===================== TRADUZIONI INTERFACCIA =====================
TESTI = {
    "en": {
        "page_title": "Protein sequence analyzer",
        "title": "🧬 Protein sequence analyzer",
        "caption": "Paste a protein sequence or a FASTA and press *Analyze*.",
        "lang_label": "Language",
        "input_label": "Paste your protein sequence here (plain text or FASTA format)",
        "input_placeholder": "Paste a protein sequence, or a full FASTA (with the > line)",
        "btn_example": "🧪 Use KCNQ1",
        "btn_analyze": "🔬 Analyze",
        "err_no_aa": "No valid amino acid found. Paste a protein sequence.",
        "warn_discarded": "{n} non-standard characters ignored.",
        "prov_header": "📍 Sequence provenance",
        "lbl_protein": "Protein", "lbl_organism": "Organism", "lbl_db": "Database",
        "lbl_accession": "Accession", "lbl_gene": "Gene", "lbl_taxid": "Tax ID",
        "info_no_header": ("ℹ️ No FASTA header recognized: this is a raw sequence, so its "
                           "provenance can't be determined. Paste a FASTA (with the `>` line) "
                           "to see protein and organism."),
        "res_header": "📊 Results",
        "sidebar_title": "Key data",
        "m_length": "Length", "m_mw": "Molecular weight", "m_pi": "Isoelectric point",
        "m_charge": "Charge at pH 7", "m_bonds": "Peptide bonds", "m_disulfide": "Disulfide bonds",
        "unit_aa": "aa", "disulfide_val": "max {n}", "cys_delta": "{n} Cys",
        "m_gravy": "GRAVY", "gravy_hydrophobic": "hydrophobic", "gravy_hydrophilic": "hydrophilic",
        "m_aliphatic": "Aliphatic index", "aliphatic_delta": "↑ = more thermostable",
        "m_abs": "Abs 280 nm (1 g/L)", "abs_delta": "ε ≈ {v} M⁻¹cm⁻¹",
        "exp_props": "ℹ️ What these properties mean",
        "props_md": (
            "- **GRAVY** (Grand Average of Hydropathy): mean hydrophobicity "
            "(Kyte-Doolittle scale) over all residues. **Negative** → likely "
            "hydrophilic/soluble protein; **positive** → hydrophobic. Here: `{gravy}` ({idro}).\n"
            "- **Aliphatic index**: relative volume occupied by aliphatic side chains "
            "(Ala, Val, Ile, Leu). Higher values are associated with greater thermostability. "
            "Here: `{ali}`.\n"
            "- **Molar extinction coefficient at 280 nm** and **Abs of a 1 g/L solution** "
            "(ProtParam's *Abs 0.1%*): used to obtain protein concentration from a "
            "spectrophotometer reading (`conc [g/L] = A₂₈₀ / Abs(1 g/L)`). Computed from "
            "Trp, Tyr and disulfide bonds.\n"
            "  - Reduced cysteines: ε ≈ `{est_red}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_red}`)\n"
            "  - Cysteines in disulfide bonds: ε ≈ `{est_ox}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_ox}`)"
        ),
        "chart_freq_title": "Amino acid frequency",
        "chart_classes_title": "Physico-chemical classes",
        "col_code": "Code", "col_name": "Name", "col_count": "Count", "col_pct": "Percentage",
        "hv_count": "Count", "hv_pct": "Percentage", "hv_aa": "amino acids",
        "exp_grouping": "📖 How amino acids are grouped",
        "grouping_md": (
            "- **Hydrophobic:** A, V, L, I, P, F, M, W\n"
            "- **Polar:** G, S, T, C, Y, N, Q\n"
            "- **Positively charged (+):** K, R, H\n"
            "- **Negatively charged (−):** D, E"
        ),
        "exp_freq_table": "📋 Full frequency table",
        "btn_csv": "📥 Download frequency table (.csv)",
        "csv_prefix": "frequencies_", "default_name": "sequence",
        "mem_header": "🧬 Membrane domains — Kyte-Doolittle hydrophobicity",
        "m_tm": "Transmembrane domains",
        "profile_title": "Hydrophobicity profile (window 19 aa)",
        "struct3d_header": "🧬 3D structure (AlphaFold)",
        "scale_label": "Scale", "window_label": "Window (aa)",
        "scale_tm_note": "Transmembrane regions are detected with Kyte-Doolittle.",
        "crosshl_hint": "💡 Hover the profile to light up the matching residues in 3D — and hover the 3D model to mark the position on the profile.",
        "struct3d_spinner": "Loading 3D structure…",
        "struct3d_plddt": "Colour = pLDDT confidence: blue = high, red = low.",
        "struct3d_source": "Source: AlphaFold DB · {id}",
        "struct3d_notfound": "No AlphaFold model available for {acc}.",
        "struct3d_need_acc": "Load the protein via a UniProt accession (UniProt tab, or a FASTA with a UniProt header) to see its 3D structure.",
        "ax_pos": "position (residue)", "ax_hydro": "average hydrophobicity",
        "hydro_legend": "Hydrophobicity (KD)",
        "tm_threshold": "TM threshold (1.6)",
        "exp_tm_detail": "🔍 Hydrophobic regions detail and sequence map",
        "tm_seg": "Segment", "tm_res": "Residues", "tm_len": "Length",
        "tm_avg": "Avg. hydrophobicity", "tm_sub": "Subsequence",
        "no_tm": "No transmembrane segment detected.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>blue</span> "
                       "= transmembrane residues &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>green</span> = hydrophilic residues (water-exposed)"),
        "btn_report": "⬇️ Download report (HTML)", "report_prefix": "report_",
        "exp_export_header": "🧩 Export charts (vector SVG / PDF)",
        "exp_export_hint": "High-resolution vector files — ideal for a thesis or exam report.",
        "exp_export_format": "Format", "exp_export_prepare": "Prepare files",
        "exp_export_spinner": "Generating files…", "exp_export_ready": "Download:",
        "exp_export_error": "Export failed. Is the 'kaleido' package installed?",
        "caption_pi": ("ℹ️ The isoelectric point is an approximate theoretical value: "
                       "it depends on the pKa set used."),
        # --- report HTML ---
        "r_title": "🧬 Protein analysis report",
        "r_provenance": "Provenance", "r_metrics": "Metrics",
        "r_length": "Length", "r_mw": "Molecular weight", "r_pi": "Isoelectric point",
        "r_charge": "Charge at pH 7",
        "r_bonds": "Peptide bonds", "r_cys": "Cysteines", "r_disulfide": "max {n} disulfide bonds",
        "r_gravy": "GRAVY (average hydrophobicity)", "r_aliphatic": "Aliphatic index",
        "r_ext": "Molar extinction at 280 nm", "r_ext_ox": "oxidized Cys", "r_abs": "Abs 1 g/L",
        "r_membrane": "Membrane domains (Kyte-Doolittle)",
        "r_verdict": "Verdict:", "r_tm_count": "transmembrane domain(s)",
        "r_seq_map": "Sequence map",
        "r_legend": ('<span style="background:#2b6cb0;color:#fff;padding:1px 5px">blue</span> = '
                     'transmembrane residues (hydrophobic) &nbsp; '
                     '<span style="color:#2f855a">green</span> = hydrophilic residues (water-exposed)'),
        "r_freq": "Amino acid frequency", "r_classes": "Physico-chemical classes",
        "r_footer": "Isoelectric point and membrane domains = approximate theoretical predictions.",
        "r_report_title": "Protein report",
        "upload_label": "Upload a FASTA file (.fasta / .txt)",
        "accession_label": "UniProt accession",
        "accession_ph": "Paste a UniProt accession",
        "btn_fetch": "Fetch from UniProt",
        "fetch_error": "Could not fetch '{acc}' from UniProt. Check the accession.",
        "fetch_spinner": "Fetching from UniProt…",
        "tab_paste": "📋 Paste sequence", "tab_file": "📁 Upload file",
        "tab_uniprot": "🔎 Search UniProt",
        "titration_header": "📈 Titration curve (net charge vs pH)",
        "ax_ph": "pH", "ax_charge": "net charge", "titration_pi": "pI = {v}",
        "r_titration": "Titration curve (net charge vs pH)",
        "ss_header": "🌀 Secondary structure composition (Chou-Fasman)",
        "ss_helix": "Helix", "ss_sheet": "Sheet", "ss_coil": "Coil / turn",
        "ss_note": "Indicative estimate from Chou-Fasman propensities — not a real structure prediction.",
        "ss_chart_title": "Predicted secondary structure",
        "compare_header": "🔬 Compare with a second protein",
        "compare_label": "Paste the second sequence (plain text or FASTA)",
        "compare_ph": "Paste a second protein sequence, or a full FASTA (with the > line)",
        "btn_compare": "Compare",
        "compare_metric": "Metric", "compare_a": "Sequence A", "compare_b": "Sequence B",
        "compare_err": "No valid amino acid in the second sequence.",
        "compare_chart_title": "Class composition: A vs B (%)",
        "compare_y": "% of residues",
    },
    "it": {
        "page_title": "Analizzatore di sequenze proteiche",
        "title": "🧬 Analizzatore di sequenze proteiche",
        "caption": "Incolla una sequenza proteica o un FASTA e premi *Analizza*.",
        "lang_label": "Lingua",
        "input_label": "Incolla qui la tua sequenza proteica (testo o formato FASTA)",
        "input_placeholder": "Incolla una sequenza proteica, oppure un FASTA completo (con la riga >)",
        "btn_example": "🧪 Usa KCNQ1",
        "btn_analyze": "🔬 Analizza",
        "err_no_aa": "Nessun amminoacido valido trovato. Incolla una sequenza proteica.",
        "warn_discarded": "{n} caratteri non standard ignorati.",
        "prov_header": "📍 Provenienza della sequenza",
        "lbl_protein": "Proteina", "lbl_organism": "Organismo", "lbl_db": "Banca dati",
        "lbl_accession": "Accession", "lbl_gene": "Gene", "lbl_taxid": "Tax ID",
        "info_no_header": ("ℹ️ Nessuna intestazione FASTA riconosciuta: è una sequenza grezza, "
                           "quindi non posso dirne la provenienza. Incolla un FASTA (con la riga `>`) "
                           "per vedere proteina e organismo."),
        "res_header": "📊 Risultati",
        "sidebar_title": "Dati chiave",
        "m_length": "Lunghezza", "m_mw": "Peso molecolare", "m_pi": "Punto isoelettrico",
        "m_charge": "Carica a pH 7", "m_bonds": "Legami peptidici", "m_disulfide": "Ponti disolfuro",
        "unit_aa": "aa", "disulfide_val": "max {n}", "cys_delta": "{n} Cys",
        "m_gravy": "GRAVY", "gravy_hydrophobic": "idrofobica", "gravy_hydrophilic": "idrofilica",
        "m_aliphatic": "Indice alifatico", "aliphatic_delta": "↑ = più termostabile",
        "m_abs": "Abs 280 nm (1 g/L)", "abs_delta": "ε ≈ {v} M⁻¹cm⁻¹",
        "exp_props": "ℹ️ Cosa significano queste proprietà",
        "props_md": (
            "- **GRAVY** (Grand Average of Hydropathy): media dell'idrofobicità "
            "(scala Kyte-Doolittle) su tutti i residui. **Negativo** → proteina "
            "tendenzialmente idrofilica/solubile; **positivo** → idrofobica. Qui: `{gravy}` ({idro}).\n"
            "- **Indice alifatico**: volume relativo occupato dalle catene laterali "
            "alifatiche (Ala, Val, Ile, Leu). Valori più alti si associano a una maggiore "
            "termostabilità. Qui: `{ali}`.\n"
            "- **Coefficiente di estinzione molare a 280 nm** e **Abs di una soluzione "
            "da 1 g/L** (la *Abs 0.1%* di ProtParam): servono a ricavare la concentrazione "
            "della proteina da una lettura allo spettrofotometro "
            "(`conc [g/L] = A₂₈₀ / Abs(1 g/L)`). Calcolati da Trp, Tyr e ponti disolfuro.\n"
            "  - Cisteine ridotte: ε ≈ `{est_red}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_red}`)\n"
            "  - Cisteine in ponti disolfuro: ε ≈ `{est_ox}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_ox}`)"
        ),
        "chart_freq_title": "Frequenza degli amminoacidi",
        "chart_classes_title": "Classi chimico-fisiche",
        "col_code": "Codice", "col_name": "Nome", "col_count": "Conteggio", "col_pct": "Percentuale",
        "hv_count": "Conteggio", "hv_pct": "Percentuale", "hv_aa": "amminoacidi",
        "exp_grouping": "📖 Come sono raggruppati gli amminoacidi",
        "grouping_md": (
            "- **Idrofobici:** A, V, L, I, P, F, M, W\n"
            "- **Polari:** G, S, T, C, Y, N, Q\n"
            "- **Carichi positivi (+):** K, R, H\n"
            "- **Carichi negativi (−):** D, E"
        ),
        "exp_freq_table": "📋 Tabella completa delle frequenze",
        "btn_csv": "📥 Scarica tabella frequenze (.csv)",
        "csv_prefix": "frequenze_", "default_name": "sequenza",
        "mem_header": "🧬 Domini di membrana — idrofobicità di Kyte-Doolittle",
        "m_tm": "Domini transmembrana",
        "profile_title": "Profilo di idrofobicità (finestra 19 aa)",
        "struct3d_header": "🧬 Struttura 3D (AlphaFold)",
        "scale_label": "Scala", "window_label": "Finestra (aa)",
        "scale_tm_note": "Le regioni transmembrana sono rilevate con Kyte-Doolittle.",
        "crosshl_hint": "💡 Passa il mouse sul profilo per illuminare i residui corrispondenti nel 3D — e viceversa, passa sul modello 3D per segnare la posizione sul profilo.",
        "struct3d_spinner": "Caricamento struttura 3D…",
        "struct3d_plddt": "Colore = confidenza pLDDT: blu = alta, rosso = bassa.",
        "struct3d_source": "Fonte: AlphaFold DB · {id}",
        "struct3d_notfound": "Nessun modello AlphaFold disponibile per {acc}.",
        "struct3d_need_acc": "Carica la proteina tramite accession UniProt (tab UniProt, o un FASTA con intestazione UniProt) per vederne la struttura 3D.",
        "ax_pos": "posizione (residuo)", "ax_hydro": "idrofobicità media",
        "hydro_legend": "Idrofobicità (KD)",
        "tm_threshold": "soglia TM (1.6)",
        "exp_tm_detail": "🔍 Dettaglio regioni idrofobiche e mappa della sequenza",
        "tm_seg": "Segmento", "tm_res": "Residui", "tm_len": "Lunghezza",
        "tm_avg": "Idrofob. media", "tm_sub": "Sottosequenza",
        "no_tm": "Nessun segmento transmembrana rilevato.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>azzurro</span> "
                       "= residui transmembrana &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>verde</span> = residui idrofilici (esposti all'acqua)"),
        "btn_report": "⬇️ Scarica il report (HTML)", "report_prefix": "report_",
        "exp_export_header": "🧩 Esporta grafici (vettoriale SVG / PDF)",
        "exp_export_hint": "File vettoriali ad alta risoluzione — ideali per una tesi o una relazione d'esame.",
        "exp_export_format": "Formato", "exp_export_prepare": "Prepara i file",
        "exp_export_spinner": "Generazione dei file…", "exp_export_ready": "Scarica:",
        "exp_export_error": "Esportazione non riuscita. Il pacchetto 'kaleido' è installato?",
        "caption_pi": ("ℹ️ Il punto isoelettrico è un valore teorico approssimato: "
                       "dipende dal set di pKa usato."),
        # --- report HTML ---
        "r_title": "🧬 Report di analisi proteica",
        "r_provenance": "Provenienza", "r_metrics": "Metriche",
        "r_length": "Lunghezza", "r_mw": "Peso molecolare", "r_pi": "Punto isoelettrico",
        "r_charge": "Carica a pH 7",
        "r_bonds": "Legami peptidici", "r_cys": "Cisteine", "r_disulfide": "max {n} ponti disolfuro",
        "r_gravy": "GRAVY (idrofobicità media)", "r_aliphatic": "Indice alifatico",
        "r_ext": "Estinzione molare a 280 nm", "r_ext_ox": "Cys ossidate", "r_abs": "Abs 1 g/L",
        "r_membrane": "Domini di membrana (Kyte-Doolittle)",
        "r_verdict": "Verdetto:", "r_tm_count": "dominio/i transmembrana",
        "r_seq_map": "Mappa della sequenza",
        "r_legend": ('<span style="background:#2b6cb0;color:#fff;padding:1px 5px">azzurro</span> = '
                     'residui transmembrana (idrofobici) &nbsp; '
                     '<span style="color:#2f855a">verde</span> = residui idrofilici (esposti all\'acqua)'),
        "r_freq": "Frequenza degli amminoacidi", "r_classes": "Classi chimico-fisiche",
        "r_footer": "Punto isoelettrico e domini di membrana = predizioni teoriche approssimate.",
        "r_report_title": "Report proteico",
        "upload_label": "Carica un file FASTA (.fasta / .txt)",
        "accession_label": "Accession UniProt",
        "accession_ph": "Incolla un accession UniProt",
        "btn_fetch": "Scarica da UniProt",
        "fetch_error": "Impossibile scaricare '{acc}' da UniProt. Controlla l'accession.",
        "fetch_spinner": "Scaricamento da UniProt…",
        "tab_paste": "📋 Incolla sequenza", "tab_file": "📁 Carica file",
        "tab_uniprot": "🔎 Cerca su UniProt",
        "titration_header": "📈 Curva di titolazione (carica netta vs pH)",
        "ax_ph": "pH", "ax_charge": "carica netta", "titration_pi": "pI = {v}",
        "r_titration": "Curva di titolazione (carica netta vs pH)",
        "ss_header": "🌀 Composizione in struttura secondaria (Chou-Fasman)",
        "ss_helix": "Elica", "ss_sheet": "Foglietto β", "ss_coil": "Coil / turn",
        "ss_note": "Stima indicativa dalle propensioni di Chou-Fasman — non è una vera predizione strutturale.",
        "ss_chart_title": "Struttura secondaria stimata",
        "compare_header": "🔬 Confronta con una seconda proteina",
        "compare_label": "Incolla la seconda sequenza (testo o FASTA)",
        "compare_ph": "Incolla una seconda sequenza proteica, oppure un FASTA completo (con la riga >)",
        "btn_compare": "Confronta",
        "compare_metric": "Metrica", "compare_a": "Sequenza A", "compare_b": "Sequenza B",
        "compare_err": "Nessun amminoacido valido nella seconda sequenza.",
        "compare_chart_title": "Composizione in classi: A vs B (%)",
        "compare_y": "% dei residui",
    },
    "es": {
        "page_title": "Analizador de secuencias proteicas",
        "title": "🧬 Analizador de secuencias proteicas",
        "caption": "Pega una secuencia proteica o un FASTA y pulsa *Analizar*.",
        "lang_label": "Idioma",
        "input_label": "Pega aquí tu secuencia proteica (texto plano o formato FASTA)",
        "input_placeholder": "Pega una secuencia proteica, o un FASTA completo (con la línea >)",
        "btn_example": "🧪 Usar KCNQ1",
        "btn_analyze": "🔬 Analizar",
        "err_no_aa": "No se encontró ningún aminoácido válido. Pega una secuencia proteica.",
        "warn_discarded": "{n} caracteres no estándar ignorados.",
        "prov_header": "📍 Procedencia de la secuencia",
        "lbl_protein": "Proteína", "lbl_organism": "Organismo", "lbl_db": "Base de datos",
        "lbl_accession": "Accession", "lbl_gene": "Gen", "lbl_taxid": "Tax ID",
        "info_no_header": ("ℹ️ No se reconoció ninguna cabecera FASTA: es una secuencia sin "
                           "procesar, por lo que no se puede determinar su procedencia. Pega un "
                           "FASTA (con la línea `>`) para ver la proteína y el organismo."),
        "res_header": "📊 Resultados",
        "sidebar_title": "Datos clave",
        "m_length": "Longitud", "m_mw": "Peso molecular", "m_pi": "Punto isoeléctrico",
        "m_charge": "Carga a pH 7", "m_bonds": "Enlaces peptídicos", "m_disulfide": "Puentes disulfuro",
        "unit_aa": "aa", "disulfide_val": "máx {n}", "cys_delta": "{n} Cys",
        "m_gravy": "GRAVY", "gravy_hydrophobic": "hidrofóbica", "gravy_hydrophilic": "hidrofílica",
        "m_aliphatic": "Índice alifático", "aliphatic_delta": "↑ = más termoestable",
        "m_abs": "Abs 280 nm (1 g/L)", "abs_delta": "ε ≈ {v} M⁻¹cm⁻¹",
        "exp_props": "ℹ️ Qué significan estas propiedades",
        "props_md": (
            "- **GRAVY** (Grand Average of Hydropathy): hidrofobicidad media "
            "(escala Kyte-Doolittle) sobre todos los residuos. **Negativo** → proteína "
            "probablemente hidrofílica/soluble; **positivo** → hidrofóbica. Aquí: `{gravy}` ({idro}).\n"
            "- **Índice alifático**: volumen relativo ocupado por las cadenas laterales "
            "alifáticas (Ala, Val, Ile, Leu). Valores más altos se asocian con mayor "
            "termoestabilidad. Aquí: `{ali}`.\n"
            "- **Coeficiente de extinción molar a 280 nm** y **Abs de una solución de 1 g/L** "
            "(la *Abs 0.1%* de ProtParam): sirven para obtener la concentración de la proteína "
            "a partir de una lectura en el espectrofotómetro (`conc [g/L] = A₂₈₀ / Abs(1 g/L)`). "
            "Calculados a partir de Trp, Tyr y puentes disulfuro.\n"
            "  - Cisteínas reducidas: ε ≈ `{est_red}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_red}`)\n"
            "  - Cisteínas en puentes disulfuro: ε ≈ `{est_ox}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_ox}`)"
        ),
        "chart_freq_title": "Frecuencia de aminoácidos",
        "chart_classes_title": "Clases fisicoquímicas",
        "col_code": "Código", "col_name": "Nombre", "col_count": "Recuento", "col_pct": "Porcentaje",
        "hv_count": "Recuento", "hv_pct": "Porcentaje", "hv_aa": "aminoácidos",
        "exp_grouping": "📖 Cómo se agrupan los aminoácidos",
        "grouping_md": (
            "- **Hidrofóbicos:** A, V, L, I, P, F, M, W\n"
            "- **Polares:** G, S, T, C, Y, N, Q\n"
            "- **Carga positiva (+):** K, R, H\n"
            "- **Carga negativa (−):** D, E"
        ),
        "exp_freq_table": "📋 Tabla completa de frecuencias",
        "btn_csv": "📥 Descargar tabla de frecuencias (.csv)",
        "csv_prefix": "frecuencias_", "default_name": "secuencia",
        "mem_header": "🧬 Dominios de membrana — hidrofobicidad de Kyte-Doolittle",
        "m_tm": "Dominios transmembrana",
        "profile_title": "Perfil de hidrofobicidad (ventana 19 aa)",
        "struct3d_header": "🧬 Estructura 3D (AlphaFold)",
        "scale_label": "Escala", "window_label": "Ventana (aa)",
        "scale_tm_note": "Las regiones transmembrana se detectan con Kyte-Doolittle.",
        "crosshl_hint": "💡 Pasa el ratón por el perfil para iluminar los residuos correspondientes en 3D — y al revés, pasa por el modelo 3D para marcar la posición en el perfil.",
        "struct3d_spinner": "Cargando estructura 3D…",
        "struct3d_plddt": "Color = confianza pLDDT: azul = alta, rojo = baja.",
        "struct3d_source": "Fuente: AlphaFold DB · {id}",
        "struct3d_notfound": "No hay modelo AlphaFold disponible para {acc}.",
        "struct3d_need_acc": "Carga la proteína mediante un accession de UniProt (pestaña UniProt, o un FASTA con encabezado UniProt) para ver su estructura 3D.",
        "ax_pos": "posición (residuo)", "ax_hydro": "hidrofobicidad media",
        "hydro_legend": "Hidrofobicidad (KD)",
        "tm_threshold": "umbral TM (1.6)",
        "exp_tm_detail": "🔍 Detalle de regiones hidrofóbicas y mapa de la secuencia",
        "tm_seg": "Segmento", "tm_res": "Residuos", "tm_len": "Longitud",
        "tm_avg": "Hidrofob. media", "tm_sub": "Subsecuencia",
        "no_tm": "No se detectó ningún segmento transmembrana.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>azul</span> "
                       "= residuos transmembrana &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>verde</span> = residuos hidrofílicos (expuestos al agua)"),
        "btn_report": "⬇️ Descargar el informe (HTML)", "report_prefix": "informe_",
        "exp_export_header": "🧩 Exportar gráficos (vectorial SVG / PDF)",
        "exp_export_hint": "Archivos vectoriales de alta resolución — ideales para una tesis o un informe.",
        "exp_export_format": "Formato", "exp_export_prepare": "Preparar archivos",
        "exp_export_spinner": "Generando archivos…", "exp_export_ready": "Descargar:",
        "exp_export_error": "Error de exportación. ¿Está instalado el paquete 'kaleido'?",
        "caption_pi": ("ℹ️ El punto isoeléctrico es un valor teórico aproximado: "
                       "depende del conjunto de pKa utilizado."),
        "r_title": "🧬 Informe de análisis proteico",
        "r_provenance": "Procedencia", "r_metrics": "Métricas",
        "r_length": "Longitud", "r_mw": "Peso molecular", "r_pi": "Punto isoeléctrico",
        "r_charge": "Carga a pH 7",
        "r_bonds": "Enlaces peptídicos", "r_cys": "Cisteínas", "r_disulfide": "máx {n} puentes disulfuro",
        "r_gravy": "GRAVY (hidrofobicidad media)", "r_aliphatic": "Índice alifático",
        "r_ext": "Extinción molar a 280 nm", "r_ext_ox": "Cys oxidadas", "r_abs": "Abs 1 g/L",
        "r_membrane": "Dominios de membrana (Kyte-Doolittle)",
        "r_verdict": "Veredicto:", "r_tm_count": "dominio(s) transmembrana",
        "r_seq_map": "Mapa de la secuencia",
        "r_legend": ('<span style="background:#2b6cb0;color:#fff;padding:1px 5px">azul</span> = '
                     'residuos transmembrana (hidrofóbicos) &nbsp; '
                     '<span style="color:#2f855a">verde</span> = residuos hidrofílicos (expuestos al agua)'),
        "r_freq": "Frecuencia de aminoácidos", "r_classes": "Clases fisicoquímicas",
        "r_footer": "Punto isoeléctrico y dominios de membrana = predicciones teóricas aproximadas.",
        "r_report_title": "Informe proteico",
        "upload_label": "Sube un archivo FASTA (.fasta / .txt)",
        "accession_label": "Accession de UniProt",
        "accession_ph": "Pega un accession de UniProt",
        "btn_fetch": "Descargar de UniProt",
        "fetch_error": "No se pudo descargar '{acc}' de UniProt. Comprueba el accession.",
        "fetch_spinner": "Descargando de UniProt…",
        "tab_paste": "📋 Pegar secuencia", "tab_file": "📁 Subir archivo",
        "tab_uniprot": "🔎 Buscar en UniProt",
        "titration_header": "📈 Curva de titulación (carga neta vs pH)",
        "ax_ph": "pH", "ax_charge": "carga neta", "titration_pi": "pI = {v}",
        "r_titration": "Curva de titulación (carga neta vs pH)",
        "ss_header": "🌀 Composición de estructura secundaria (Chou-Fasman)",
        "ss_helix": "Hélice", "ss_sheet": "Lámina β", "ss_coil": "Coil / giro",
        "ss_note": "Estimación indicativa a partir de las propensiones de Chou-Fasman — no es una predicción estructural real.",
        "ss_chart_title": "Estructura secundaria estimada",
        "compare_header": "🔬 Comparar con una segunda proteína",
        "compare_label": "Pega la segunda secuencia (texto o FASTA)",
        "compare_ph": "Pega una segunda secuencia proteica, o un FASTA completo (con la línea >)",
        "btn_compare": "Comparar",
        "compare_metric": "Métrica", "compare_a": "Secuencia A", "compare_b": "Secuencia B",
        "compare_err": "Ningún aminoácido válido en la segunda secuencia.",
        "compare_chart_title": "Composición de clases: A vs B (%)",
        "compare_y": "% de residuos",
    },
    "de": {
        "page_title": "Proteinsequenz-Analysator",
        "title": "🧬 Proteinsequenz-Analysator",
        "caption": "Füge eine Proteinsequenz oder ein FASTA ein und drücke *Analysieren*.",
        "lang_label": "Sprache",
        "input_label": "Füge hier deine Proteinsequenz ein (reiner Text oder FASTA-Format)",
        "input_placeholder": "Füge eine Proteinsequenz ein oder ein vollständiges FASTA (mit der >-Zeile)",
        "btn_example": "🧪 KCNQ1 verwenden",
        "btn_analyze": "🔬 Analysieren",
        "err_no_aa": "Keine gültige Aminosäure gefunden. Füge eine Proteinsequenz ein.",
        "warn_discarded": "{n} nicht standardmäßige Zeichen ignoriert.",
        "prov_header": "📍 Herkunft der Sequenz",
        "lbl_protein": "Protein", "lbl_organism": "Organismus", "lbl_db": "Datenbank",
        "lbl_accession": "Accession", "lbl_gene": "Gen", "lbl_taxid": "Tax-ID",
        "info_no_header": ("ℹ️ Kein FASTA-Header erkannt: Dies ist eine reine Sequenz, daher kann "
                           "die Herkunft nicht bestimmt werden. Füge ein FASTA (mit der `>`-Zeile) "
                           "ein, um Protein und Organismus zu sehen."),
        "res_header": "📊 Ergebnisse",
        "sidebar_title": "Kerndaten",
        "m_length": "Länge", "m_mw": "Molekulargewicht", "m_pi": "Isoelektrischer Punkt",
        "m_charge": "Ladung bei pH 7", "m_bonds": "Peptidbindungen", "m_disulfide": "Disulfidbrücken",
        "unit_aa": "AS", "disulfide_val": "max {n}", "cys_delta": "{n} Cys",
        "m_gravy": "GRAVY", "gravy_hydrophobic": "hydrophob", "gravy_hydrophilic": "hydrophil",
        "m_aliphatic": "Aliphatischer Index", "aliphatic_delta": "↑ = thermostabiler",
        "m_abs": "Abs 280 nm (1 g/L)", "abs_delta": "ε ≈ {v} M⁻¹cm⁻¹",
        "exp_props": "ℹ️ Was diese Eigenschaften bedeuten",
        "props_md": (
            "- **GRAVY** (Grand Average of Hydropathy): mittlere Hydrophobizität "
            "(Kyte-Doolittle-Skala) über alle Reste. **Negativ** → wahrscheinlich "
            "hydrophiles/lösliches Protein; **positiv** → hydrophob. Hier: `{gravy}` ({idro}).\n"
            "- **Aliphatischer Index**: relatives Volumen der aliphatischen Seitenketten "
            "(Ala, Val, Ile, Leu). Höhere Werte sind mit größerer Thermostabilität verbunden. "
            "Hier: `{ali}`.\n"
            "- **Molarer Extinktionskoeffizient bei 280 nm** und **Abs einer 1-g/L-Lösung** "
            "(ProtParams *Abs 0.1%*): dienen dazu, die Proteinkonzentration aus einer "
            "Spektrophotometer-Messung zu ermitteln (`conc [g/L] = A₂₈₀ / Abs(1 g/L)`). "
            "Berechnet aus Trp, Tyr und Disulfidbrücken.\n"
            "  - Reduzierte Cysteine: ε ≈ `{est_red}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_red}`)\n"
            "  - Cysteine in Disulfidbrücken: ε ≈ `{est_ox}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_ox}`)"
        ),
        "chart_freq_title": "Aminosäure-Häufigkeit",
        "chart_classes_title": "Physikalisch-chemische Klassen",
        "col_code": "Code", "col_name": "Name", "col_count": "Anzahl", "col_pct": "Prozent",
        "hv_count": "Anzahl", "hv_pct": "Prozent", "hv_aa": "Aminosäuren",
        "exp_grouping": "📖 Wie die Aminosäuren gruppiert werden",
        "grouping_md": (
            "- **Hydrophob:** A, V, L, I, P, F, M, W\n"
            "- **Polar:** G, S, T, C, Y, N, Q\n"
            "- **Positiv geladen (+):** K, R, H\n"
            "- **Negativ geladen (−):** D, E"
        ),
        "exp_freq_table": "📋 Vollständige Häufigkeitstabelle",
        "btn_csv": "📥 Häufigkeitstabelle herunterladen (.csv)",
        "csv_prefix": "haeufigkeit_", "default_name": "sequenz",
        "mem_header": "🧬 Membrandomänen — Kyte-Doolittle-Hydrophobizität",
        "m_tm": "Transmembrandomänen",
        "profile_title": "Hydrophobizitätsprofil (Fenster 19 AS)",
        "struct3d_header": "🧬 3D-Struktur (AlphaFold)",
        "scale_label": "Skala", "window_label": "Fenster (AS)",
        "scale_tm_note": "Transmembranregionen werden mit Kyte-Doolittle erkannt.",
        "crosshl_hint": "💡 Fahre über das Profil, um die passenden Reste in 3D hervorzuheben — und über das 3D-Modell, um die Position im Profil zu markieren.",
        "struct3d_spinner": "3D-Struktur wird geladen…",
        "struct3d_plddt": "Farbe = pLDDT-Konfidenz: blau = hoch, rot = niedrig.",
        "struct3d_source": "Quelle: AlphaFold DB · {id}",
        "struct3d_notfound": "Kein AlphaFold-Modell für {acc} verfügbar.",
        "struct3d_need_acc": "Lade das Protein über eine UniProt-Accession (UniProt-Tab oder FASTA mit UniProt-Header), um die 3D-Struktur zu sehen.",
        "ax_pos": "Position (Rest)", "ax_hydro": "mittlere Hydrophobizität",
        "hydro_legend": "Hydrophobizität (KD)",
        "tm_threshold": "TM-Schwelle (1.6)",
        "exp_tm_detail": "🔍 Details zu hydrophoben Regionen und Sequenzkarte",
        "tm_seg": "Segment", "tm_res": "Reste", "tm_len": "Länge",
        "tm_avg": "Mittl. Hydrophob.", "tm_sub": "Teilsequenz",
        "no_tm": "Kein Transmembransegment erkannt.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>blau</span> "
                       "= Transmembranreste &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>grün</span> = hydrophile Reste (wasserexponiert)"),
        "btn_report": "⬇️ Bericht herunterladen (HTML)", "report_prefix": "bericht_",
        "exp_export_header": "🧩 Diagramme exportieren (vektor SVG / PDF)",
        "exp_export_hint": "Hochauflösende Vektordateien — ideal für Abschlussarbeit oder Bericht.",
        "exp_export_format": "Format", "exp_export_prepare": "Dateien vorbereiten",
        "exp_export_spinner": "Dateien werden erzeugt…", "exp_export_ready": "Herunterladen:",
        "exp_export_error": "Export fehlgeschlagen. Ist das Paket 'kaleido' installiert?",
        "caption_pi": ("ℹ️ Der isoelektrische Punkt ist ein ungefährer theoretischer Wert: "
                       "er hängt vom verwendeten pKa-Satz ab."),
        "r_title": "🧬 Proteinanalyse-Bericht",
        "r_provenance": "Herkunft", "r_metrics": "Kennzahlen",
        "r_length": "Länge", "r_mw": "Molekulargewicht", "r_pi": "Isoelektrischer Punkt",
        "r_charge": "Ladung bei pH 7",
        "r_bonds": "Peptidbindungen", "r_cys": "Cysteine", "r_disulfide": "max {n} Disulfidbrücken",
        "r_gravy": "GRAVY (mittlere Hydrophobizität)", "r_aliphatic": "Aliphatischer Index",
        "r_ext": "Molare Extinktion bei 280 nm", "r_ext_ox": "oxidierte Cys", "r_abs": "Abs 1 g/L",
        "r_membrane": "Membrandomänen (Kyte-Doolittle)",
        "r_verdict": "Befund:", "r_tm_count": "Transmembrandomäne(n)",
        "r_seq_map": "Sequenzkarte",
        "r_legend": ('<span style="background:#2b6cb0;color:#fff;padding:1px 5px">blau</span> = '
                     'Transmembranreste (hydrophob) &nbsp; '
                     '<span style="color:#2f855a">grün</span> = hydrophile Reste (wasserexponiert)'),
        "r_freq": "Aminosäure-Häufigkeit", "r_classes": "Physikalisch-chemische Klassen",
        "r_footer": "Isoelektrischer Punkt und Membrandomänen = ungefähre theoretische Vorhersagen.",
        "r_report_title": "Proteinbericht",
        "upload_label": "FASTA-Datei hochladen (.fasta / .txt)",
        "accession_label": "UniProt-Accession",
        "accession_ph": "Füge eine UniProt-Accession ein",
        "btn_fetch": "Von UniProt laden",
        "fetch_error": "'{acc}' konnte nicht von UniProt geladen werden. Prüfe die Accession.",
        "fetch_spinner": "Lade von UniProt…",
        "tab_paste": "📋 Sequenz einfügen", "tab_file": "📁 Datei hochladen",
        "tab_uniprot": "🔎 UniProt durchsuchen",
        "titration_header": "📈 Titrationskurve (Nettoladung vs. pH)",
        "ax_ph": "pH", "ax_charge": "Nettoladung", "titration_pi": "pI = {v}",
        "r_titration": "Titrationskurve (Nettoladung vs. pH)",
        "ss_header": "🌀 Sekundärstruktur-Zusammensetzung (Chou-Fasman)",
        "ss_helix": "Helix", "ss_sheet": "Faltblatt", "ss_coil": "Coil / Turn",
        "ss_note": "Indikative Schätzung aus Chou-Fasman-Propensitäten — keine echte Strukturvorhersage.",
        "ss_chart_title": "Geschätzte Sekundärstruktur",
        "compare_header": "🔬 Mit einem zweiten Protein vergleichen",
        "compare_label": "Zweite Sequenz einfügen (Text oder FASTA)",
        "compare_ph": "Füge eine zweite Proteinsequenz ein oder ein vollständiges FASTA (mit der >-Zeile)",
        "btn_compare": "Vergleichen",
        "compare_metric": "Kennzahl", "compare_a": "Sequenz A", "compare_b": "Sequenz B",
        "compare_err": "Keine gültige Aminosäure in der zweiten Sequenz.",
        "compare_chart_title": "Klassenzusammensetzung: A vs B (%)",
        "compare_y": "% der Reste",
    },
    "fr": {
        "page_title": "Analyseur de séquences protéiques",
        "title": "🧬 Analyseur de séquences protéiques",
        "caption": "Collez une séquence protéique ou un FASTA et appuyez sur *Analyser*.",
        "lang_label": "Langue",
        "input_label": "Collez ici votre séquence protéique (texte brut ou format FASTA)",
        "input_placeholder": "Collez une séquence protéique, ou un FASTA complet (avec la ligne >)",
        "btn_example": "🧪 Utiliser KCNQ1",
        "btn_analyze": "🔬 Analyser",
        "err_no_aa": "Aucun acide aminé valide trouvé. Collez une séquence protéique.",
        "warn_discarded": "{n} caractères non standard ignorés.",
        "prov_header": "📍 Provenance de la séquence",
        "lbl_protein": "Protéine", "lbl_organism": "Organisme", "lbl_db": "Base de données",
        "lbl_accession": "Accession", "lbl_gene": "Gène", "lbl_taxid": "Tax ID",
        "info_no_header": ("ℹ️ Aucun en-tête FASTA reconnu : c'est une séquence brute, sa "
                           "provenance ne peut donc pas être déterminée. Collez un FASTA (avec la "
                           "ligne `>`) pour voir la protéine et l'organisme."),
        "res_header": "📊 Résultats",
        "sidebar_title": "Données clés",
        "m_length": "Longueur", "m_mw": "Poids moléculaire", "m_pi": "Point isoélectrique",
        "m_charge": "Charge à pH 7", "m_bonds": "Liaisons peptidiques", "m_disulfide": "Ponts disulfure",
        "unit_aa": "aa", "disulfide_val": "max {n}", "cys_delta": "{n} Cys",
        "m_gravy": "GRAVY", "gravy_hydrophobic": "hydrophobe", "gravy_hydrophilic": "hydrophile",
        "m_aliphatic": "Indice aliphatique", "aliphatic_delta": "↑ = plus thermostable",
        "m_abs": "Abs 280 nm (1 g/L)", "abs_delta": "ε ≈ {v} M⁻¹cm⁻¹",
        "exp_props": "ℹ️ Ce que signifient ces propriétés",
        "props_md": (
            "- **GRAVY** (Grand Average of Hydropathy) : hydrophobicité moyenne "
            "(échelle Kyte-Doolittle) sur tous les résidus. **Négatif** → protéine "
            "probablement hydrophile/soluble ; **positif** → hydrophobe. Ici : `{gravy}` ({idro}).\n"
            "- **Indice aliphatique** : volume relatif occupé par les chaînes latérales "
            "aliphatiques (Ala, Val, Ile, Leu). Des valeurs plus élevées sont associées à une "
            "plus grande thermostabilité. Ici : `{ali}`.\n"
            "- **Coefficient d'extinction molaire à 280 nm** et **Abs d'une solution à 1 g/L** "
            "(l'*Abs 0.1%* de ProtParam) : servent à obtenir la concentration de la protéine à "
            "partir d'une lecture au spectrophotomètre (`conc [g/L] = A₂₈₀ / Abs(1 g/L)`). "
            "Calculés à partir de Trp, Tyr et ponts disulfure.\n"
            "  - Cystéines réduites : ε ≈ `{est_red}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_red}`)\n"
            "  - Cystéines en ponts disulfure : ε ≈ `{est_ox}` M⁻¹cm⁻¹ (Abs 1 g/L = `{abs_ox}`)"
        ),
        "chart_freq_title": "Fréquence des acides aminés",
        "chart_classes_title": "Classes physico-chimiques",
        "col_code": "Code", "col_name": "Nom", "col_count": "Nombre", "col_pct": "Pourcentage",
        "hv_count": "Nombre", "hv_pct": "Pourcentage", "hv_aa": "acides aminés",
        "exp_grouping": "📖 Comment les acides aminés sont regroupés",
        "grouping_md": (
            "- **Hydrophobes :** A, V, L, I, P, F, M, W\n"
            "- **Polaires :** G, S, T, C, Y, N, Q\n"
            "- **Chargés positivement (+) :** K, R, H\n"
            "- **Chargés négativement (−) :** D, E"
        ),
        "exp_freq_table": "📋 Tableau complet des fréquences",
        "btn_csv": "📥 Télécharger le tableau des fréquences (.csv)",
        "csv_prefix": "frequences_", "default_name": "sequence",
        "mem_header": "🧬 Domaines membranaires — hydrophobicité de Kyte-Doolittle",
        "m_tm": "Domaines transmembranaires",
        "profile_title": "Profil d'hydrophobicité (fenêtre 19 aa)",
        "struct3d_header": "🧬 Structure 3D (AlphaFold)",
        "scale_label": "Échelle", "window_label": "Fenêtre (aa)",
        "scale_tm_note": "Les régions transmembranaires sont détectées avec Kyte-Doolittle.",
        "crosshl_hint": "💡 Survolez le profil pour illuminer les résidus correspondants en 3D — et survolez le modèle 3D pour marquer la position sur le profil.",
        "struct3d_spinner": "Chargement de la structure 3D…",
        "struct3d_plddt": "Couleur = confiance pLDDT : bleu = élevée, rouge = faible.",
        "struct3d_source": "Source : AlphaFold DB · {id}",
        "struct3d_notfound": "Aucun modèle AlphaFold disponible pour {acc}.",
        "struct3d_need_acc": "Chargez la protéine via une accession UniProt (onglet UniProt, ou un FASTA avec en-tête UniProt) pour voir sa structure 3D.",
        "ax_pos": "position (résidu)", "ax_hydro": "hydrophobicité moyenne",
        "hydro_legend": "Hydrophobicité (KD)",
        "tm_threshold": "seuil TM (1.6)",
        "exp_tm_detail": "🔍 Détail des régions hydrophobes et carte de la séquence",
        "tm_seg": "Segment", "tm_res": "Résidus", "tm_len": "Longueur",
        "tm_avg": "Hydrophob. moy.", "tm_sub": "Sous-séquence",
        "no_tm": "Aucun segment transmembranaire détecté.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>bleu</span> "
                       "= résidus transmembranaires &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>vert</span> = résidus hydrophiles (exposés à l'eau)"),
        "btn_report": "⬇️ Télécharger le rapport (HTML)", "report_prefix": "rapport_",
        "exp_export_header": "🧩 Exporter les graphiques (vectoriel SVG / PDF)",
        "exp_export_hint": "Fichiers vectoriels haute résolution — idéaux pour un mémoire ou un rapport.",
        "exp_export_format": "Format", "exp_export_prepare": "Préparer les fichiers",
        "exp_export_spinner": "Génération des fichiers…", "exp_export_ready": "Télécharger :",
        "exp_export_error": "Échec de l'export. Le paquet « kaleido » est-il installé ?",
        "caption_pi": ("ℹ️ Le point isoélectrique est une valeur théorique approximative : "
                       "elle dépend du jeu de pKa utilisé."),
        "r_title": "🧬 Rapport d'analyse protéique",
        "r_provenance": "Provenance", "r_metrics": "Métriques",
        "r_length": "Longueur", "r_mw": "Poids moléculaire", "r_pi": "Point isoélectrique",
        "r_charge": "Charge à pH 7",
        "r_bonds": "Liaisons peptidiques", "r_cys": "Cystéines", "r_disulfide": "max {n} ponts disulfure",
        "r_gravy": "GRAVY (hydrophobicité moyenne)", "r_aliphatic": "Indice aliphatique",
        "r_ext": "Extinction molaire à 280 nm", "r_ext_ox": "Cys oxydées", "r_abs": "Abs 1 g/L",
        "r_membrane": "Domaines membranaires (Kyte-Doolittle)",
        "r_verdict": "Verdict :", "r_tm_count": "domaine(s) transmembranaire(s)",
        "r_seq_map": "Carte de la séquence",
        "r_legend": ('<span style="background:#2b6cb0;color:#fff;padding:1px 5px">bleu</span> = '
                     'résidus transmembranaires (hydrophobes) &nbsp; '
                     '<span style="color:#2f855a">vert</span> = résidus hydrophiles (exposés à l\'eau)'),
        "r_freq": "Fréquence des acides aminés", "r_classes": "Classes physico-chimiques",
        "r_footer": "Point isoélectrique et domaines membranaires = prédictions théoriques approximatives.",
        "r_report_title": "Rapport protéique",
        "upload_label": "Téléverser un fichier FASTA (.fasta / .txt)",
        "accession_label": "Accession UniProt",
        "accession_ph": "Collez une accession UniProt",
        "btn_fetch": "Récupérer depuis UniProt",
        "fetch_error": "Impossible de récupérer '{acc}' depuis UniProt. Vérifiez l'accession.",
        "fetch_spinner": "Récupération depuis UniProt…",
        "tab_paste": "📋 Coller la séquence", "tab_file": "📁 Téléverser un fichier",
        "tab_uniprot": "🔎 Rechercher sur UniProt",
        "titration_header": "📈 Courbe de titration (charge nette vs pH)",
        "ax_ph": "pH", "ax_charge": "charge nette", "titration_pi": "pI = {v}",
        "r_titration": "Courbe de titration (charge nette vs pH)",
        "ss_header": "🌀 Composition en structure secondaire (Chou-Fasman)",
        "ss_helix": "Hélice", "ss_sheet": "Feuillet β", "ss_coil": "Coil / coude",
        "ss_note": "Estimation indicative à partir des propensions de Chou-Fasman — pas une vraie prédiction structurale.",
        "ss_chart_title": "Structure secondaire estimée",
        "compare_header": "🔬 Comparer avec une deuxième protéine",
        "compare_label": "Collez la deuxième séquence (texte ou FASTA)",
        "compare_ph": "Collez une deuxième séquence protéique, ou un FASTA complet (avec la ligne >)",
        "btn_compare": "Comparer",
        "compare_metric": "Métrique", "compare_a": "Séquence A", "compare_b": "Séquence B",
        "compare_err": "Aucun acide aminé valide dans la deuxième séquence.",
        "compare_chart_title": "Composition en classes : A vs B (%)",
        "compare_y": "% des résidus",
    },
    "zh": {
        "page_title": "蛋白质序列分析器",
        "title": "🧬 蛋白质序列分析器",
        "caption": "粘贴一条蛋白质序列或 FASTA，然后点击 *分析*。",
        "lang_label": "语言",
        "input_label": "在此粘贴你的蛋白质序列（纯文本或 FASTA 格式）",
        "input_placeholder": "粘贴一条蛋白质序列，或完整的 FASTA（含 > 行）",
        "btn_example": "🧪 使用 KCNQ1",
        "btn_analyze": "🔬 分析",
        "err_no_aa": "未找到有效的氨基酸。请粘贴一条蛋白质序列。",
        "warn_discarded": "已忽略 {n} 个非标准字符。",
        "prov_header": "📍 序列来源",
        "lbl_protein": "蛋白质", "lbl_organism": "物种", "lbl_db": "数据库",
        "lbl_accession": "登录号", "lbl_gene": "基因", "lbl_taxid": "Tax ID",
        "info_no_header": ("ℹ️ 未识别到 FASTA 头部：这是一条原始序列，因此无法确定其来源。"
                           "请粘贴 FASTA（含 `>` 行）以查看蛋白质和物种。"),
        "res_header": "📊 结果",
        "sidebar_title": "关键数据",
        "m_length": "长度", "m_mw": "分子量", "m_pi": "等电点",
        "m_charge": "pH 7 时的电荷", "m_bonds": "肽键", "m_disulfide": "二硫键",
        "unit_aa": "aa", "disulfide_val": "最多 {n}", "cys_delta": "{n} 个 Cys",
        "m_gravy": "GRAVY", "gravy_hydrophobic": "疏水", "gravy_hydrophilic": "亲水",
        "m_aliphatic": "脂肪族指数", "aliphatic_delta": "↑ = 更耐热",
        "m_abs": "Abs 280 nm (1 g/L)", "abs_delta": "ε ≈ {v} M⁻¹cm⁻¹",
        "exp_props": "ℹ️ 这些属性的含义",
        "props_md": (
            "- **GRAVY**（总平均亲水性）：所有残基的平均疏水性（Kyte-Doolittle 标度）。"
            "**负值** → 可能为亲水/可溶蛋白；**正值** → 疏水。此处：`{gravy}`（{idro}）。\n"
            "- **脂肪族指数**：脂肪族侧链（Ala、Val、Ile、Leu）所占的相对体积。"
            "数值越高通常热稳定性越好。此处：`{ali}`。\n"
            "- **280 nm 摩尔消光系数** 和 **1 g/L 溶液的吸光度**（ProtParam 的 *Abs 0.1%*）："
            "用于根据分光光度计读数计算蛋白质浓度（`conc [g/L] = A₂₈₀ / Abs(1 g/L)`）。"
            "由 Trp、Tyr 和二硫键计算得出。\n"
            "  - 还原型半胱氨酸：ε ≈ `{est_red}` M⁻¹cm⁻¹（Abs 1 g/L = `{abs_red}`）\n"
            "  - 形成二硫键的半胱氨酸：ε ≈ `{est_ox}` M⁻¹cm⁻¹（Abs 1 g/L = `{abs_ox}`）"
        ),
        "chart_freq_title": "氨基酸频率",
        "chart_classes_title": "理化分类",
        "col_code": "代码", "col_name": "名称", "col_count": "数量", "col_pct": "百分比",
        "hv_count": "数量", "hv_pct": "百分比", "hv_aa": "个氨基酸",
        "exp_grouping": "📖 氨基酸如何分组",
        "grouping_md": (
            "- **疏水性：** A, V, L, I, P, F, M, W\n"
            "- **极性：** G, S, T, C, Y, N, Q\n"
            "- **带正电 (+)：** K, R, H\n"
            "- **带负电 (−)：** D, E"
        ),
        "exp_freq_table": "📋 完整频率表",
        "btn_csv": "📥 下载频率表 (.csv)",
        "csv_prefix": "frequencies_", "default_name": "sequence",
        "mem_header": "🧬 膜结构域 — Kyte-Doolittle 疏水性",
        "m_tm": "跨膜结构域",
        "profile_title": "疏水性分布（窗口 19 aa）",
        "struct3d_header": "🧬 三维结构 (AlphaFold)",
        "scale_label": "标度", "window_label": "窗口 (aa)",
        "scale_tm_note": "跨膜区域使用 Kyte-Doolittle 检测。",
        "crosshl_hint": "💡 将鼠标悬停在曲线上可在 3D 中高亮对应残基；悬停在 3D 模型上可在曲线上标记对应位置。",
        "struct3d_spinner": "正在加载三维结构…",
        "struct3d_plddt": "颜色 = pLDDT 置信度：蓝色=高，红色=低。",
        "struct3d_source": "来源：AlphaFold DB · {id}",
        "struct3d_notfound": "未找到 {acc} 的 AlphaFold 模型。",
        "struct3d_need_acc": "通过 UniProt 登录号加载蛋白质（UniProt 标签页，或带 UniProt 头的 FASTA）以查看其三维结构。",
        "ax_pos": "位置（残基）", "ax_hydro": "平均疏水性",
        "hydro_legend": "疏水性 (KD)",
        "tm_threshold": "TM 阈值 (1.6)",
        "exp_tm_detail": "🔍 疏水区域详情与序列图",
        "tm_seg": "片段", "tm_res": "残基", "tm_len": "长度",
        "tm_avg": "平均疏水性", "tm_sub": "子序列",
        "no_tm": "未检测到跨膜片段。",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>蓝色</span> "
                       "= 跨膜残基 &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>绿色</span> = 亲水残基（暴露于水）"),
        "btn_report": "⬇️ 下载报告 (HTML)", "report_prefix": "report_",
        "exp_export_header": "🧩 导出图表（矢量 SVG / PDF）",
        "exp_export_hint": "高分辨率矢量文件 —— 适合论文或考试报告。",
        "exp_export_format": "格式", "exp_export_prepare": "准备文件",
        "exp_export_spinner": "正在生成文件…", "exp_export_ready": "下载：",
        "exp_export_error": "导出失败。是否已安装 'kaleido' 包？",
        "caption_pi": "ℹ️ 等电点是一个近似的理论值：取决于所用的 pKa 集合。",
        "r_title": "🧬 蛋白质分析报告",
        "r_provenance": "来源", "r_metrics": "指标",
        "r_length": "长度", "r_mw": "分子量", "r_pi": "等电点", "r_charge": "pH 7 时的电荷",
        "r_bonds": "肽键", "r_cys": "半胱氨酸", "r_disulfide": "最多 {n} 个二硫键",
        "r_gravy": "GRAVY（平均疏水性）", "r_aliphatic": "脂肪族指数",
        "r_ext": "280 nm 摩尔消光系数", "r_ext_ox": "氧化型 Cys", "r_abs": "Abs 1 g/L",
        "r_membrane": "膜结构域 (Kyte-Doolittle)",
        "r_verdict": "判定：", "r_tm_count": "个跨膜结构域",
        "r_seq_map": "序列图",
        "r_legend": ('<span style="background:#2b6cb0;color:#fff;padding:1px 5px">蓝色</span> = '
                     '跨膜残基（疏水）&nbsp; '
                     '<span style="color:#2f855a">绿色</span> = 亲水残基（暴露于水）'),
        "r_freq": "氨基酸频率", "r_classes": "理化分类",
        "r_footer": "等电点和膜结构域 = 近似的理论预测。",
        "r_report_title": "蛋白质报告",
        "upload_label": "上传 FASTA 文件 (.fasta / .txt)",
        "accession_label": "UniProt 登录号",
        "accession_ph": "粘贴一个 UniProt 登录号",
        "btn_fetch": "从 UniProt 获取",
        "fetch_error": "无法从 UniProt 获取 '{acc}'。请检查登录号。",
        "fetch_spinner": "正在从 UniProt 获取…",
        "tab_paste": "📋 粘贴序列", "tab_file": "📁 上传文件",
        "tab_uniprot": "🔎 搜索 UniProt",
        "titration_header": "📈 滴定曲线（净电荷 vs pH）",
        "ax_ph": "pH", "ax_charge": "净电荷", "titration_pi": "pI = {v}",
        "r_titration": "滴定曲线（净电荷 vs pH）",
        "ss_header": "🌀 二级结构组成（Chou-Fasman）",
        "ss_helix": "α-螺旋", "ss_sheet": "β-折叠", "ss_coil": "无规卷曲 / 转角",
        "ss_note": "基于 Chou-Fasman 倾向性的指示性估计 —— 并非真正的结构预测。",
        "ss_chart_title": "预测的二级结构",
        "compare_header": "🔬 与第二条蛋白质比较",
        "compare_label": "粘贴第二条序列（纯文本或 FASTA）",
        "compare_ph": "粘贴第二条蛋白质序列，或完整的 FASTA（含 > 行）",
        "btn_compare": "比较",
        "compare_metric": "指标", "compare_a": "序列 A", "compare_b": "序列 B",
        "compare_err": "第二条序列中没有有效的氨基酸。",
        "compare_chart_title": "类别组成：A vs B（%）",
        "compare_y": "残基百分比",
    },
}


def chiave(k):
    return k if isinstance(k, str) else chr(k)


def parse_fasta_header(testo: str, lang: str):
    """Legge la prima riga '>' e ne ricava la provenienza. None se assente."""
    riga = next((l.strip() for l in testo.splitlines() if l.startswith(">")), None)
    if not riga:
        return None
    h = riga[1:].strip()
    info = {"intestazione": h}

    # Formato UniProt:  sp|P51787|KCNQ1_HUMAN  Descrizione OS=... OX=... GN=...
    m = re.match(r"(sp|tr)\|([^|]+)\|(\S+)\s*(.*)", h)
    if m:
        db, acc, entry, resto = m.groups()
        info["banca_dati"] = (DB_DESC.get(db, {}).get(lang)
                              or DB_DESC.get(db, {}).get("en", db))
        info["accession"] = acc
        info["entry"] = entry
        info["url"] = f"https://www.uniprot.org/uniprotkb/{acc}"
        # campi OS= (organismo), OX= (tax id), GN= (gene)
        for tag, campo in [("OS", "organismo"), ("OX", "tax_id"), ("GN", "gene")]:
            mm = re.search(rf"{tag}=(.+?)(?:\s+\w\w=|$)", resto)
            if mm:
                info[campo] = mm.group(1).strip()
        # nome proteina = testo prima del primo TAG=
        nome = re.split(r"\s+\w\w=", resto)[0].strip()
        if nome:
            info["proteina"] = nome
        # organismo dal suffisso dell'entry (es. _HUMAN) se OS= manca
        if "organismo" not in info and "_" in entry:
            suff = entry.split("_")[-1]
            if suff in ORG_SCI:
                comune = ORG_COMUNE.get(suff, {}).get(lang, "")
                info["organismo"] = (f"{ORG_SCI[suff]} ({comune})"
                                     if comune else ORG_SCI[suff])
        return info

    # Formato NCBI:  NP_000209.2 descrizione [Homo sapiens]
    m = re.match(r"(\S+)\s+(.*)", h)
    if m:
        acc, resto = m.groups()
        info["accession"] = acc
        org = re.search(r"\[(.+?)\]", resto)
        if org:
            info["organismo"] = org.group(1)
            resto = resto[:org.start()].strip()
        info["proteina"] = resto or h
        info["banca_dati"] = ("NCBI / GenBank (probabile)" if lang == "it"
                              else "NCBI / GenBank (likely)")
        return info

    info["proteina"] = h
    return info


def valida(seq: str):
    righe = [r for r in seq.splitlines() if not r.startswith(">")]
    seq = "".join(righe).upper()
    seq = "".join(c for c in seq if not c.isspace())
    pulita = "".join(c for c in seq if c in AA_VALIDI)
    return pulita, len(seq) - len(pulita)


def carica_esempio():
    """Restituisce il FASTA KCNQ1 COMPLETO (con intestazione) per la demo."""
    path = os.path.join(os.path.dirname(__file__), "KCNQ1.fasta")
    if os.path.isfile(path):
        with open(path) as f:
            return f.read()
    return ">esempio\nMAAASSPPRAERKRWGWGRLPGARRGSAGLAKKCPFSLELAEGGPAGGALYAPIAPGAPGP"


def fetch_uniprot(acc: str):
    """Scarica il FASTA di un accession da UniProt. Gestisce anche gli accession
    secondari/obsoleti (es. P62158 → P0DP23 Calmodulina): il .fasta diretto in quei
    casi risponde 200 ma vuoto, quindi si ricade sull'API di ricerca (sec_acc).
    Risolvendo all'entry corrente si sistema anche il viewer 3D (l'header avrà
    l'accession primario). None se fallisce del tutto."""
    acc = acc.strip()

    def _get(url):
        req = urllib.request.Request(url, headers={"User-Agent": "analizzatore-proteine"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            return resp.read().decode("utf-8", "ignore")

    try:
        # 1) accession diretto
        testo = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
        if testo.lstrip().startswith(">"):
            return testo
        # 2) accession secondario/obsoleto → risolvi all'entry corrente
        q = urllib.parse.quote(f"sec_acc:{acc}")
        testo = _get("https://rest.uniprot.org/uniprotkb/search"
                     f"?query={q}&format=fasta&size=1")
        return testo if testo.lstrip().startswith(">") else None
    except Exception:
        return None


# UA simil-browser: l'API di AlphaFold rifiuta (403) gli User-Agent generici.
_UA_BROWSER = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
               "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124 Safari/537.36")


@st.cache_data(show_spinner=False)
def fetch_alphafold(acc: str):
    """Scarica il modello predetto da AlphaFold DB per un accession UniProt.
    Ritorna {'pdb', 'id'} oppure None (nessun modello / rete non disponibile).
    In cache: lo stesso accession non viene riscaricato a ogni rerun."""
    acc = (acc or "").strip()
    if not acc:
        return None

    def _get(url, timeout):
        req = urllib.request.Request(url, headers={"User-Agent": _UA_BROWSER})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", "ignore")

    try:
        meta = json.loads(_get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", 20))
        if not meta:
            return None
        entry = meta[0]
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            return None
        pdb = _get(pdb_url, 30)
        if "ATOM" not in pdb:
            return None
        return {"pdb": pdb, "id": entry.get("entryId", "")}
    except Exception:
        return None


def visualizzatore_3d(pdb_text, altezza=360):
    """Mostra una struttura PDB con 3Dmol.js in un widget interattivo (rotazione,
    zoom). Catena 'cartoon' colorata per pLDDT (B-factor): blu = confidenza alta,
    rosso = bassa. La struttura è iniettata come testo (niente fetch dal browser,
    quindi nessun problema di CORS). Sfondo trasparente: si adatta al tema."""
    html_v = """
<div id="vwrap" style="width:100%%;height:%dpx;position:relative;"></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
(function(){
  function start(){
    if(typeof $3Dmol === 'undefined'){ return setTimeout(start, 80); }
    const el = document.getElementById('vwrap');
    const v = $3Dmol.createViewer(el, {backgroundAlpha: 0});
    v.addModel(%s, 'pdb');
    v.setStyle({}, {cartoon: {colorscheme: {prop:'b', gradient:'roygb', min:50, max:90}}});
    v.zoomTo();
    v.render();
    window.addEventListener('resize', function(){ v.resize(); });
  }
  start();
})();
</script>
""" % (altezza, json.dumps(pdb_text))
    components.html(html_v, height=altezza + 12)


_CROSS_TEMPLATE = """
<div id="cwrap" style="display:flex;gap:10px;width:100%;height:__ALT__px;
 font-family:-apple-system,'Segoe UI',sans-serif;">
  <div id="cplot" style="flex:3;min-width:0;height:100%;"></div>
  <div id="cmol" style="flex:2;min-width:0;height:100%;position:relative;">
    <div style="position:absolute;top:6px;left:0;right:0;text-align:center;
     font-size:13px;font-weight:600;opacity:.75;pointer-events:none;z-index:5;
     color:__FONTC__;">__TITOLO3D__</div>
  </div>
</div>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
(function(){
  const D = __PAYLOAD__;
  const HL = '#ffd400';                       // colore di evidenziazione condiviso
  function baseShapes(){
    const sh = [];
    if(D.soglia !== null){
      sh.push({type:'line',xref:'paper',x0:0,x1:1,y0:D.soglia,y1:D.soglia,
               line:{color:'#c0504d',dash:'dash',width:1.5}});
    } else {
      sh.push({type:'line',xref:'paper',x0:0,x1:1,y0:0,y1:0,
               line:{color:'#888',dash:'dot',width:1}});
    }
    D.seg.forEach(function(s){
      sh.push({type:'rect',xref:'x',yref:'paper',x0:s[0],x1:s[1],y0:0,y1:1,
               fillcolor:'#2b6cb0',opacity:0.18,line:{width:0},layer:'below'});
    });
    return sh;
  }
  const fontc = D.scuro ? '#fafafa' : '#262730';
  const layout = {title:{text:D.titolo,font:{size:14}},margin:{t:42,b:44,l:56,r:12},
    xaxis:{title:D.asseX,zeroline:false},yaxis:{title:D.asseY,zeroline:false},
    shapes:baseShapes(),paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',
    font:{color:fontc},hovermode:'x',showlegend:false};
  const trace = {x:D.x,y:D.y,type:'scatter',mode:'lines',line:{color:D.linea,width:2},
    hovertemplate:'%{x}<br>%{y:.2f}<extra></extra>'};

  let viewer=null, baseStyle=null;
  function set3D(lo,hi){
    if(!viewer) return;
    viewer.setStyle({}, baseStyle);
    if(lo!=null){
      const arr=[]; for(let r=lo;r<=hi;r++) arr.push(r);
      viewer.setStyle({resi:arr},{cartoon:{color:HL},stick:{color:HL,radius:0.25}});
    }
    viewer.render();
  }
  function markPlot(resi){
    const sh = baseShapes();
    if(resi!=null){
      sh.push({type:'line',xref:'x',yref:'paper',x0:resi,x1:resi,y0:0,y1:1,
               line:{color:HL,width:2}});
    }
    Plotly.relayout('cplot',{shapes:sh});
  }

  Plotly.newPlot('cplot',[trace],layout,{displayModeBar:false,responsive:true})
   .then(function(){
    const gd = document.getElementById('cplot');
    gd.on('plotly_hover',function(ev){
      const x = Math.round(ev.points[0].x);
      let lo=x, hi=x;                          // se dentro un TM, evidenzio tutto il dominio
      for(const s of D.seg){ if(x>=s[0] && x<=s[1]){ lo=s[0]; hi=s[1]; break; } }
      set3D(lo,hi);
    });
    gd.on('plotly_unhover',function(){ set3D(null,null); });
  });

  function startMol(){
    if(typeof $3Dmol==='undefined' || typeof Plotly==='undefined'){ return setTimeout(startMol,80); }
    viewer = $3Dmol.createViewer(document.getElementById('cmol'),{backgroundAlpha:0});
    viewer.addModel(D.pdb,'pdb');
    baseStyle = {cartoon:{colorscheme:{prop:'b',gradient:'roygb',min:50,max:90}}};
    viewer.setStyle({}, baseStyle);
    viewer.zoomTo();
    viewer.render();
    viewer.setHoverable({}, true,
      function(atom){ if(atom){ markPlot(atom.resi);
        viewer.setStyle({}, baseStyle);
        viewer.setStyle({resi:[atom.resi]},{cartoon:{color:HL},stick:{color:HL,radius:0.25}});
        viewer.render(); } },
      function(){ markPlot(null); set3D(null,null); });
    viewer.render();
    window.addEventListener('resize', function(){ viewer.resize(); });
  }
  startMol();
})();
</script>
"""


def visualizzatore_combinato(pos, sco, segmenti, pdb_text, *, asse_x, asse_y,
                             titolo, titolo_3d, soglia, linea, scuro, altezza=400):
    """Profilo 2D (Plotly.js) + struttura 3D (3Dmol.js) nello STESSO iframe, con
    cross-highlighting: passando sul grafico si illuminano i residui corrispondenti
    sul modello 3D (e i domini TM per intero), e viceversa il passaggio sul 3D
    segna la posizione sul grafico. Stare nello stesso frame evita la fragile
    comunicazione cross-iframe."""
    data = {
        "x": list(pos), "y": [round(v, 4) for v in sco],
        "seg": [[int(s.inizio), int(s.fine)] for s in segmenti],
        "asseX": asse_x, "asseY": asse_y, "titolo": titolo,
        "soglia": soglia, "linea": linea, "scuro": bool(scuro), "pdb": pdb_text,
    }
    html_w = (_CROSS_TEMPLATE.replace("__ALT__", str(altezza))
              .replace("__FONTC__", "#fafafa" if scuro else "#262730")
              .replace("__TITOLO3D__", html.escape(str(titolo_3d)))
              .replace("__PAYLOAD__", json.dumps(data)))
    components.html(html_w, height=altezza + 16)


def profilo_finestra(seq, tabella, finestra):
    """Profilo a media mobile di una scala per-residuo, con lo stesso allineamento
    del core C++ (finestra a sinistra [i, i+w-1], valore al residuo centrale,
    1-based). Così la curva combacia col rilevamento TM. Ritorna (pos, punteggi)."""
    N = len(seq)
    pos, sco = [], []
    if N < finestra:
        return pos, sco
    for i in range(0, N - finestra + 1):
        s = sum(tabella.get(seq[k], 0.0) for k in range(i, i + finestra))
        sco.append(s / finestra)
        pos.append(i + finestra // 2 + 1)
    return pos, sco


def heatmap_sequenza(seq, segmenti):
    """Sequenza intera in HTML: residui transmembrana evidenziati in azzurro."""
    in_tm = [False] * len(seq)
    for s in segmenti:                       # inizio/fine sono 1-based
        for i in range(s.inizio - 1, min(s.fine, len(seq))):
            in_tm[i] = True
    pezzi = []
    for i, a in enumerate(seq):
        if in_tm[i]:
            pezzi.append(f"<span style='background:#2b6cb0;color:#fff'>{a}</span>")
        else:
            pezzi.append(f"<span style='color:#2f855a'>{a}</span>")
    corpo = "".join(pezzi)
    return (f"<div style='font-family:ui-monospace,Menlo,monospace;font-size:13px;"
            f"line-height:1.9;word-break:break-all;letter-spacing:1px'>{corpo}</div>")


def costruisci_report_html(meta, r, prof, fig_bar, fig_pie, fig_idro, fig_tit, heatmap, T, lang):
    """Report HTML autonomo con i grafici Plotly interattivi (per il download)."""
    # Il report ha sfondo bianco: i grafici vanno resi sempre in tema chiaro, a
    # prescindere dal tema dell'app (in dark mode testo/linee chiari sparirebbero
    # sul bianco). Riallineo qui tutte le figure prima di serializzarle.
    for _f in (fig_bar, fig_pie, fig_idro, fig_tit):
        _f.update_layout(template="plotly_white", paper_bgcolor="white",
                         plot_bgcolor="white", font_color="#222")
    fig_idro.update_traces(line_color="#1d2733")     # linea scura, leggibile sul bianco

    prov = ""
    if meta:
        voci = [(T["lbl_protein"], meta.get("proteina")),
                (T["lbl_organism"], meta.get("organismo")),
                (T["lbl_db"], meta.get("banca_dati")),
                (T["lbl_accession"], meta.get("accession")),
                (T["lbl_gene"], meta.get("gene"))]
        prov = "".join(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>"
                       for k, v in voci if v)
        prov = f"<h2>{T['r_provenance']}</h2><table>{prov}</table>"
    kpi = f"""
      <ul>
        <li>{T['r_length']}: <b>{r.lunghezza}</b> {T['unit_aa']}</li>
        <li>{T['r_mw']}: <b>{r.peso_molecolare:,.1f} Da</b> ({r.peso_molecolare/1000:.1f} kDa)</li>
        <li>{T['r_pi']}: <b>{r.punto_isoelettrico:.2f}</b></li>
        <li>{T['r_charge']}: <b>{r.carica_a_pH7:+.2f}</b></li>
        <li>{T['r_bonds']}: <b>{r.legami_peptidici}</b> &middot; {T['r_cys']}: <b>{r.cisteine}</b> ({T['r_disulfide'].format(n=r.ponti_disolfuro_max)})</li>
        <li>{T['r_gravy']}: <b>{r.gravy:+.3f}</b> &middot; {T['r_aliphatic']}: <b>{r.indice_alifatico:.1f}</b></li>
        <li>{T['r_ext']}: <b>{r.estinzione_ox:,.0f}</b> M&middot;cm<sup>-1</sup> ({T['r_ext_ox']}) &middot; {T['r_abs']}: <b>{r.abs280_ox:.3f}</b></li>
      </ul>"""
    verdetto = VERDETTO_LABEL.get(prof.verdetto, {}).get(lang, prof.verdetto)
    righe_tm = "".join(
        f"<tr><td>TM{i}</td><td>{s.inizio}-{s.fine}</td><td>{s.lunghezza} {T['unit_aa']}</td>"
        f"<td>{s.idro_media:.2f}</td><td style='font-family:monospace'>{s.sottosequenza}</td></tr>"
        for i, s in enumerate(prof.segmenti, 1)
    )
    membrana = f"""
      <h2>{T['r_membrane']}</h2>
      <p><b>{T['r_verdict']}</b> {verdetto} &middot; <b>{prof.n_domini}</b> {T['r_tm_count']}</p>
      {fig_idro.to_html(full_html=False, include_plotlyjs=False)}
      <table><tr><th>{T['tm_seg']}</th><th>{T['tm_res']}</th><th>{T['tm_len']}</th><th>{T['tm_avg']}</th><th>{T['tm_sub']}</th></tr>{righe_tm}</table>
      <h3>{T['r_seq_map']}</h3>
      <p style="font-size:12px;color:#666">{T['r_legend']}</p>
      {heatmap}
    """

    bar = fig_bar.to_html(full_html=False, include_plotlyjs="cdn")
    pie = fig_pie.to_html(full_html=False, include_plotlyjs=False)
    tit = fig_tit.to_html(full_html=False, include_plotlyjs=False)
    return f"""<!DOCTYPE html><html lang="{lang}"><head><meta charset="utf-8">
<title>{T['r_report_title']}</title>
<style>body{{font-family:-apple-system,Segoe UI,sans-serif;max-width:900px;
margin:20px auto;color:#1d2733}} table{{border-collapse:collapse;width:100%}}
td,th{{padding:4px 12px;border-bottom:1px solid #eee;text-align:left;font-size:13px}}
h1{{color:#3b7dd8}}</style>
</head><body>
<h1>{T['r_title']}</h1>
{prov}
<h2>{T['r_metrics']}</h2>{kpi}
<h2>{T['r_titration']}</h2>{tit}
{membrana}
<h2>{T['r_freq']}</h2>{bar}
<h2>{T['r_classes']}</h2>{pie}
<p style="color:#888;font-size:12px">{T['r_footer']}</p>
</body></html>"""


def _tema_grafico():
    """Colori coerenti col tema (chiaro/scuro) di Streamlit. Lo sfondo è
    sempre trasparente, così il grafico eredita lo sfondo dell'app."""
    # Un eventuale tema forzato in .streamlit/config.toml è quello realmente
    # renderizzato: ha la precedenza. Senza config si segue il tema attivo nel
    # frontend (menu "Settings" → Chiaro/Scuro/Sistema), riportato da st.context.
    base = st.get_option("theme.base")
    if base in ("light", "dark"):
        tipo = base
    else:
        try:
            tipo = st.context.theme.type
        except Exception:
            tipo = None
        tipo = tipo or "light"
    scuro = (tipo == "dark")
    return {"sfondo": "rgba(0,0,0,0)",
            "testo": "#fafafa" if scuro else "#262730",
            "template": "plotly_dark" if scuro else "plotly_white"}


def _hex_rgba(hex_c, alpha):
    h = hex_c.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def grafico_con_hover(fig, tipo, div_id, colori=None, height=430):
    """Rende un grafico Plotly in cui l'elemento sotto il cursore resta in
    primo piano e gli altri si attenuano (listener hover/unhover in JS)."""
    tema = _tema_grafico()
    fig.update_layout(paper_bgcolor=tema["sfondo"], plot_bgcolor=tema["sfondo"],
                      font_color=tema["testo"], template=tema["template"],
                      margin=dict(t=60, b=20, l=10, r=10))

    if tipo == "bar":
        # le barre sono punti di un'unica traccia: attenuo via marker.opacity.
        # 'reset' è ritardato (debounce) per non lampeggiare tra barre adiacenti.
        post = (
            "var gd=document.getElementById('%s');"
            "var n=gd.data[0].x.length;var t=null;"
            "function evidenzia(pn){var op=[];"
            "  for(var i=0;i<n;i++){op.push(pn<0?1.0:(i===pn?1.0:0.18));}"
            "  Plotly.restyle(gd,{'marker.opacity':[op]});}"
            "gd.on('plotly_hover',function(d){"
            "  if(t){clearTimeout(t);t=null;}evidenzia(d.points[0].pointNumber);});"
            "gd.on('plotly_unhover',function(){"
            "  if(t)clearTimeout(t);t=setTimeout(function(){evidenzia(-1);},110);});"
        ) % div_id
    else:  # pie / torta: attenuo cambiando i colori delle fette non puntate
        pieni = colori
        sbiaditi = [_hex_rgba(c, 0.16) for c in colori]
        post = (
            "var gd=document.getElementById('%s');"
            "var pieni=%s;var sbiaditi=%s;var t=null;"
            "function evidenzia(pn){"
            "  var c=pieni.map(function(x,i){return (pn<0||i===pn)?pieni[i]:sbiaditi[i];});"
            "  Plotly.restyle(gd,{'marker.colors':[c]});}"
            "gd.on('plotly_hover',function(d){"
            "  if(t){clearTimeout(t);t=null;}evidenzia(d.points[0].pointNumber);});"
            "gd.on('plotly_unhover',function(){"
            "  if(t)clearTimeout(t);t=setTimeout(function(){evidenzia(-1);},110);});"
        ) % (div_id, json.dumps(pieni), json.dumps(sbiaditi))

    inner = fig.to_html(include_plotlyjs="cdn", full_html=False, div_id=div_id,
                        config={"displayModeBar": True, "displaylogo": False,
                                "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                                "responsive": True},
                        post_script=post)
    # transizioni CSS: l'opacità delle barre e il colore/opacità delle fette
    # cambiano in modo fluido invece che a scatti.
    stile = (
        "<style>"
        ".js-plotly-plot .barlayer .point path,"
        ".js-plotly-plot .point path{transition:opacity .18s ease;}"
        ".js-plotly-plot .pielayer .slice path,"
        ".js-plotly-plot .slice path{"
        "transition:fill .22s ease,fill-opacity .22s ease;}"
        "</style>"
    )
    page = (f"<!DOCTYPE html><html><head><meta charset='utf-8'>{stile}</head>"
            f"<body style='margin:0;background:{tema['sfondo']}'>{inner}</body></html>")
    st.iframe(page, height=height)


def _fig_per_stampa(fig, line_scuro=False):
    """Copia 'pronta per la stampa' di una figura Plotly: tema chiaro, sfondo
    bianco e dimensioni fisse — adatta all'export vettoriale (SVG/PDF) da inserire
    in una tesi o relazione. Non altera la figura mostrata a schermo."""
    f = go.Figure(fig)
    f.update_layout(template="plotly_white", paper_bgcolor="white",
                    plot_bgcolor="white", font_color="#222",
                    width=1000, height=520, margin=dict(t=60, b=60, l=70, r=40))
    if line_scuro:                       # profilo idrofobicità: linea scura sul bianco
        f.update_traces(line_color="#1d2733")
    return f


def griglia_metriche(cards, ncol=3):
    """Rende un elenco di metriche come griglia di card 'dashboard':
    bordo leggero, leggera ombreggiatura, icona colorata dedicata.
    Ogni elemento di `cards` è un dict con chiavi:
        icon  -> emoji/icona
        color -> colore esadecimale dell'icona (#rrggbb)
        label -> etichetta della metrica
        value -> valore principale
        sub   -> riga secondaria opzionale (unità / nota)
    Il layout è responsivo: `ncol` colonne su desktop, 2 sotto i 760px, 1 su
    mobile. I colori NON sono cotti in Python: il CSS ha varianti chiara e scura,
    e il theme-watcher JS sceglie quale applicare mettendo la classe `app-dark`
    su <html>. Così al cambio tema le card si aggiornano subito, senza rerun."""
    css = (
        "<style>"
        ".mcards{display:grid;grid-template-columns:repeat(%d,minmax(0,1fr));"
        "gap:14px;margin:.25rem 0 1.1rem;}"
        ".mcards .mc{background:#ffffff;border:1px solid #e7eaf0;border-radius:14px;"
        "padding:15px 17px;box-shadow:0 1px 2px rgba(16,24,40,.04),0 2px 6px rgba(16,24,40,.06);"
        "display:flex;align-items:flex-start;gap:13px;"
        "transition:box-shadow .18s ease,transform .18s ease,background .18s ease,"
        "border-color .18s ease;}"
        ".mcards .mc:hover{transform:translateY(-2px);box-shadow:0 6px 18px rgba(16,24,40,.12);}"
        ".mcards .mc .ico{flex:0 0 auto;width:40px;height:40px;border-radius:11px;"
        "display:flex;align-items:center;justify-content:center;font-size:20px;"
        "line-height:1;}"
        ".mcards .mc .bd{display:flex;flex-direction:column;gap:1px;min-width:0;}"
        ".mcards .mc .lab{font-size:.72rem;font-weight:600;letter-spacing:.04em;"
        "text-transform:uppercase;color:#667085;white-space:nowrap;overflow:hidden;"
        "text-overflow:ellipsis;}"
        ".mcards .mc .val{font-size:1.5rem;font-weight:700;line-height:1.2;color:#1f2430;}"
        ".mcards .mc .sub{font-size:.78rem;color:#98a2b3;}"
        # --- variante scura (classe 'app-dark' su <html>, messa dal watcher) ---
        "html.app-dark .mcards .mc{background:#1b2230;border-color:rgba(255,255,255,0.10);"
        "box-shadow:0 1px 3px rgba(0,0,0,.45);}"
        "html.app-dark .mcards .mc:hover{box-shadow:0 6px 18px rgba(0,0,0,.55);}"
        "html.app-dark .mcards .mc .lab{color:#9aa0a6;}"
        "html.app-dark .mcards .mc .val{color:#fafafa;}"
        "html.app-dark .mcards .mc .sub{color:#8a9099;}"
        "@media(max-width:760px){.mcards{grid-template-columns:repeat(2,minmax(0,1fr));}}"
        "@media(max-width:480px){.mcards{grid-template-columns:1fr;}}"
        "</style>"
    ) % ncol

    items = []
    for c in cards:
        ico_bg = _hex_rgba(c["color"], 0.14)
        sub = (f'<div class="sub">{html.escape(str(c["sub"]))}</div>'
               if c.get("sub") else "")
        items.append(
            '<div class="mc">'
            f'<div class="ico" style="background:{ico_bg};color:{c["color"]}">'
            f'{c["icon"]}</div>'
            '<div class="bd">'
            f'<div class="lab">{html.escape(str(c["label"]))}</div>'
            f'<div class="val">{html.escape(str(c["value"]))}</div>'
            f'{sub}</div></div>'
        )
    st.markdown(css + '<div class="mcards">' + "".join(items) + "</div>",
                unsafe_allow_html=True)


# CSS statico del popover 'Dati chiave'. Colori in due varianti: chiara (default)
# e scura (selettore 'html.app-dark', attivato dal theme-watcher JS). Così il
# popover segue il tema in tempo reale, senza dover ricalcolare nulla in Python.
_KP_CSS = """
#keypanel-root .kp-trigger{position:fixed;left:14px;top:50%;
 transform:translateY(-50%);width:48px;height:48px;border-radius:14px;display:flex;
 align-items:center;justify-content:center;font-size:24px;cursor:pointer;
 z-index:100000;background:#ffffff;border:1px solid #e7eaf0;
 box-shadow:0 4px 16px rgba(0,0,0,.18);
 transition:box-shadow .2s,transform .2s,background .2s,border-color .2s;}
html.app-dark #keypanel-root .kp-trigger{background:#1b2230;
 border-color:rgba(255,255,255,0.12);box-shadow:0 4px 16px rgba(0,0,0,.45);}
#keypanel-root .kp-trigger:hover{transform:translateY(-50%) scale(1.05);
 box-shadow:0 4px 16px rgba(0,0,0,.25),0 0 18px rgba(59,125,216,.55);}
#keypanel-root .kp-panel{position:fixed;left:74px;top:50%;
 transform-origin:left center;
 transform:translateY(-50%) translateX(-14px) scale(.94);width:252px;z-index:100000;
 background:rgba(255,255,255,0.97);backdrop-filter:blur(12px);
 -webkit-backdrop-filter:blur(12px);border:1px solid #e7eaf0;border-radius:16px;
 padding:14px 16px;box-shadow:0 12px 38px rgba(16,24,40,.18);opacity:0;visibility:hidden;
 transition:opacity .3s ease,transform .3s cubic-bezier(.16,1,.3,1),
 visibility 0s linear .3s,background .25s,border-color .25s;
 font-family:-apple-system,'Segoe UI',sans-serif;}
html.app-dark #keypanel-root .kp-panel{background:rgba(27,34,48,0.92);
 border-color:rgba(255,255,255,0.12);box-shadow:0 12px 38px rgba(0,0,0,.45);}
#keypanel-root .kp-panel.open{opacity:1;visibility:visible;
 transform:translateY(-50%) translateX(0) scale(1);
 transition:opacity .3s ease,transform .3s cubic-bezier(.16,1,.3,1),
 visibility 0s,background .25s,border-color .25s;}
#keypanel-root .kp-head{display:flex;align-items:center;justify-content:space-between;}
#keypanel-root .kp-title{font-weight:700;font-size:.98rem;color:#1f2430;}
html.app-dark #keypanel-root .kp-title{color:#fafafa;}
#keypanel-root .kp-close{cursor:pointer;color:#667085;font-size:1.25rem;
 line-height:1;opacity:.65;}#keypanel-root .kp-close:hover{opacity:1;}
html.app-dark #keypanel-root .kp-close{color:#aab2bd;}
#keypanel-root .kp-name{font-size:.82rem;font-weight:600;color:#1f2430;
 opacity:.75;margin:2px 0 10px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
html.app-dark #keypanel-root .kp-name{color:#fafafa;}
#keypanel-root .kp-row{display:flex;align-items:center;gap:8px;padding:6px 0;
 border-top:1px solid #e7eaf0;}
html.app-dark #keypanel-root .kp-row{border-top-color:rgba(255,255,255,0.12);}
#keypanel-root .kp-row:first-of-type{border-top:none;}
#keypanel-root .kp-dot{width:9px;height:9px;border-radius:50%;flex:0 0 auto;}
#keypanel-root .kp-l{font-size:.78rem;color:#667085;flex:1;}
html.app-dark #keypanel-root .kp-l{color:#aab2bd;}
#keypanel-root .kp-v{font-size:.9rem;font-weight:700;color:#1f2430;}
html.app-dark #keypanel-root .kp-v{color:#fafafa;}
"""


def prepara_layout_pannello():
    """Nasconde la sidebar nativa di Streamlit (non più usata) e riserva un po'
    di spazio a sinistra per l'icona-trigger del popover 'Dati chiave'."""
    st.markdown(
        "<style>"
        "[data-testid='stSidebar'],[data-testid='stSidebarCollapseButton'],"
        "[data-testid='stSidebarCollapsedControl'],[data-testid='collapsedControl']"
        "{display:none!important;}"
        "section[data-testid='stMain'] .block-container{padding-left:64px;}"
        "</style>",
        unsafe_allow_html=True,
    )


def inietta_theme_watcher():
    """Tiene i componenti custom (card, popover, icona) allineati al tema attivo
    di Streamlit SENZA rerun né refresh. Il cambio tema dal menu è lato browser:
    qui rilevo chiaro/scuro dal colore di sfondo dell'app e metto/tolgo la classe
    'app-dark' su <html> — il CSS (varianti chiara/scura) fa il resto all'istante.
    Un breve polling reagisce ai cambi di tema in tempo reale. Idempotente: il
    flag su window evita di creare più intervalli a ogni rerun."""
    components.html("""
<script>
const doc = window.parent.document;
const win = doc.defaultView;
if(!win.__themeWatcher){
  win.__themeWatcher = true;
  const root = doc.documentElement;
  function bg(el){ return (getComputedStyle(el).backgroundColor || '').match(/[\\d.]+/g); }
  function isDark(){
    let m = null;
    const cand = [doc.querySelector('[data-testid="stApp"]'), doc.body, root];
    for(const el of cand){
      if(!el) continue;
      const c = bg(el);
      if(c && (c.length < 4 || parseFloat(c[3]) > 0)){ m = c; break; }
    }
    if(!m) return false;
    const lum = 0.2126*(+m[0]) + 0.7152*(+m[1]) + 0.0722*(+m[2]);
    return lum < 128;
  }
  function apply(){ root.classList.toggle('app-dark', isDark()); }
  apply();
  win.setInterval(apply, 350);
}
</script>
""", height=0)


def inietta_invio_submit():
    """Nella text_area della sequenza (dentro un st.form) fa sì che il semplice
    Invio avvii l'analisi, mentre Shift+Invio inserisce un a-capo (utile per i
    FASTA multiriga). Di default Streamlit richiederebbe Ctrl/⌘+Invio.

    st.markdown non esegue JS: inietto un listener nel documento padre via
    components.html (iframe invisibile)."""
    components.html("""
<script>
const doc = window.parent.document;
doc.querySelectorAll('[data-testid="stForm"]').forEach(function(form){
  const ta = form.querySelector('textarea');
  if(ta && !ta.dataset.enterBound){
    ta.dataset.enterBound = "1";
    ta.addEventListener('keydown', function(e){
      if(e.key === 'Enter' && !e.shiftKey && !e.isComposing){
        e.preventDefault();
        const btn = form.querySelector('[data-testid="stFormSubmitButton"] button');
        if(btn) btn.click();
      }
    });
  }
});
</script>
""", height=0)


def pannello_dati_chiave(T, meta, r):
    """Popover 'Dati chiave' sovrapposto, nascosto di default e mostrato al
    passaggio del mouse su un'icona fissa a sinistra (event listener JS
    onmouseover/onmouseout), con pulsante di chiusura (×).

    st.markdown non esegue <script>, quindi inietto HTML+JS via components.html:
    dall'iframe (height=0, invisibile) lo script costruisce il popover dentro
    window.parent.document, così fluttua sopra la pagina principale. I colori NON
    sono cotti qui: il tema lo gestisce il CSS (classe 'app-dark' su <html>)."""
    voci = [
        ("#3b7dd8", T["m_length"], f"{r.lunghezza} {T['unit_aa']}"),
        ("#7e57c2", T["m_mw"], f"{r.peso_molecolare/1000:.1f} kDa"),
        ("#c0504d", T["m_pi"], f"{r.punto_isoelettrico:.2f}"),
        ("#e8a33d", T["m_charge"], f"{r.carica_a_pH7:+.1f}"),
        ("#17a2b8", T["m_gravy"], f"{r.gravy:+.3f}"),
    ]
    righe = "".join(
        f'<div class="kp-row"><span class="kp-dot" style="background:{c}"></span>'
        f'<span class="kp-l">{html.escape(str(l))}</span>'
        f'<span class="kp-v">{html.escape(str(v))}</span></div>'
        for c, l, v in voci
    )
    nome = ""
    if meta and (meta.get("proteina") or meta.get("accession")):
        n = meta.get("proteina") or meta.get("accession")
        nome = f'<div class="kp-name">{html.escape(str(n))}</div>'
    titolo = html.escape(str(T["sidebar_title"]))

    inner = (
        f"<style>{_KP_CSS}</style>"
        '<div class="kp-trigger" id="kp-trigger" title="Dati chiave">🧬</div>'
        '<div class="kp-panel" id="kp-panel">'
        f'<div class="kp-head"><span class="kp-title">{titolo}</span>'
        '<span class="kp-close" id="kp-close">&times;</span></div>'
        f"{nome}{righe}</div>"
    )

    js = """
const doc = window.parent.document;
let root = doc.getElementById('keypanel-root');
if(!root){ root = doc.createElement('div'); root.id='keypanel-root'; doc.body.appendChild(root); }
root.innerHTML = %s;
const trigger = doc.getElementById('kp-trigger');
const panel   = doc.getElementById('kp-panel');
const closeb  = doc.getElementById('kp-close');
let t = null;
const show = () => { if(t){clearTimeout(t); t=null;} panel.classList.add('open'); };
const hide = () => { t = setTimeout(() => panel.classList.remove('open'), 200); };
trigger.addEventListener('mouseover', show);
trigger.addEventListener('mouseout',  hide);
panel.addEventListener('mouseover', show);
panel.addEventListener('mouseout',  hide);
closeb.addEventListener('click', () => { if(t) clearTimeout(t); panel.classList.remove('open'); });
""" % json.dumps(inner)

    components.html(f"<script>{js}</script>", height=0)


# ===================== INTERFACCIA =====================
st.set_page_config(page_title="Protein sequence analyzer",
                   page_icon="🧬", layout="wide",
                   initial_sidebar_state="collapsed")

# --- Selettore di lingua (default: inglese), persistito nella URL ---
LINGUE = {"English": "en", "Italiano": "it", "Español": "es",
          "Deutsch": "de", "Français": "fr", "中文": "zh"}
codici = list(LINGUE.values())
lang_url = st.query_params.get("lang", "en")
top_l, top_r = st.columns([5, 1])
with top_r:
    scelta = st.selectbox("Language / 语言", list(LINGUE.keys()),
                          index=codici.index(lang_url) if lang_url in codici else 0,
                          key="lang_sel")
lang = LINGUE[scelta]
if st.query_params.get("lang") != lang:
    st.query_params["lang"] = lang
T = TESTI[lang]

with top_l:
    st.title(T["title"])
    st.caption(T["caption"])

# Nasconde la sidebar nativa e prepara lo spazio per il popover 'Dati chiave'.
prepara_layout_pannello()
# Allinea i componenti custom al tema attivo in tempo reale (senza rerun).
inietta_theme_watcher()

# --- Ripristino dopo il refresh della pagina ---
# Un refresh del browser crea una sessione nuova e azzera lo stato in memoria.
# Per non perdere l'analisi, la sequenza viene salvata nella URL (query param):
# al primo caricamento la recuperiamo e rilanciamo l'analisi automaticamente.
if "sessione_avviata" not in st.session_state:
    st.session_state.sessione_avviata = True
    seq_salvata = st.query_params.get("seq", "")
    if seq_salvata:
        st.session_state.seq_input = seq_salvata
        st.session_state.autorun = True

# --- Sorgenti alternative: carica un file FASTA o recupera da UniProt ---
def _imposta_sequenza(testo):
    # gira PRIMA della creazione del text_area, quindi può scrivere su seq_input
    st.session_state.seq_input = testo
    st.session_state.autorun = True        # rilancia l'analisi in automatico

# Tre metodi di input in tab intercambiabili: mostra solo i campi necessari.
tab_paste, tab_file, tab_uniprot = st.tabs(
    [T["tab_paste"], T["tab_file"], T["tab_uniprot"]])

# I metodi 'file' e 'UniProt' vengono eseguiti PRIMA di creare la text_area,
# perché scrivono su st.session_state.seq_input (vedi _imposta_sequenza). L'ordine
# del codice è indipendente dall'ordine visivo dei tab definito in st.tabs().
with tab_file:
    up = st.file_uploader(T["upload_label"], type=["fasta", "fa", "faa", "txt", "seq"])
    if up is not None:
        fid = (up.name, up.size)
        if st.session_state.get("_ultimo_file") != fid:   # carica solo i file nuovi
            st.session_state._ultimo_file = fid
            _imposta_sequenza(up.getvalue().decode("utf-8", "ignore"))
with tab_uniprot:
    # Form: premendo Invio nel campo accession parte il fetch (oltre al pulsante).
    with st.form("form_uniprot", border=False):
        acc = st.text_input(T["accession_label"], placeholder=T["accession_ph"],
                            key="acc_input")
        fetch_clic = st.form_submit_button(T["btn_fetch"])
    if fetch_clic and acc.strip():
        with st.spinner(T["fetch_spinner"]):
            scaricato = fetch_uniprot(acc)
        if scaricato:
            _imposta_sequenza(scaricato)
        else:
            st.error(T["fetch_error"].format(acc=acc.strip()))

# 'Incolla sequenza' è il tab predefinito (il primo): contiene la text_area.
with tab_paste:
    # Form: nella text_area Ctrl/⌘+Invio avvia l'analisi (oltre al pulsante).
    with st.form("form_paste", border=False):
        seq_input = st.text_area(
            T["input_label"], height=130, key="seq_input",
            placeholder=T["input_placeholder"],
        )
        analizza = st.form_submit_button(T["btn_analyze"], type="primary")

# Nella casella sequenza: Invio = avvia analisi, Shift+Invio = a-capo.
inietta_invio_submit()

# autorun = analisi rilanciata in automatico dopo un refresh (vedi sopra)
autorun = st.session_state.pop("autorun", False)

# Memorizzo qual è la sequenza analizzata. Così, se cambio lingua (rerun senza
# clic su Analizza), l'analisi viene ri-renderizzata da sola nella nuova lingua,
# senza dover ripremere il pulsante.
if analizza or autorun:
    st.session_state.seq_analizzata = st.session_state.get("seq_input", "")

testo = st.session_state.get("seq_analizzata", "")
if testo:
    meta = parse_fasta_header(testo, lang)
    pulita, scartati = valida(testo)
    if not pulita:
        if "seq" in st.query_params:        # input non valido: ripulisco la URL
            del st.query_params["seq"]
        st.session_state.pop("seq_analizzata", None)
        st.error(T["err_no_aa"])
        st.stop()
    if scartati:
        st.warning(T["warn_discarded"].format(n=scartati))

    # Salvo la sequenza nella URL: così sopravvive al refresh della pagina.
    # (Le sequenze molto lunghe non entrano in una URL: in quel caso non la
    # salvo, l'analisi resta valida ma non verrà ripristinata dopo il refresh.)
    if len(testo) <= 6000:
        st.query_params["seq"] = testo
    elif "seq" in st.query_params:
        del st.query_params["seq"]

    # --- PROVENIENZA (dall'intestazione FASTA) ---
    if meta and (meta.get("proteina") or meta.get("accession")):
        st.subheader(T["prov_header"])
        p1, p2, p3 = st.columns(3)
        p1.markdown(f"**{T['lbl_protein']}**\n\n{meta.get('proteina','—')}")
        org = meta.get("organismo", "—")
        p2.markdown(f"**{T['lbl_organism']}**\n\n{org}")
        p3.markdown(f"**{T['lbl_db']}**\n\n{meta.get('banca_dati','—')}")
        dettagli = []
        if meta.get("accession"):
            acc = meta["accession"]
            if meta.get("url"):
                dettagli.append(f"{T['lbl_accession']}: [`{acc}`]({meta['url']}) ↗")
            else:
                dettagli.append(f"{T['lbl_accession']}: `{acc}`")
        if meta.get("gene"):
            dettagli.append(f"{T['lbl_gene']}: **{meta['gene']}**")
        if meta.get("tax_id"):
            dettagli.append(f"{T['lbl_taxid']}: `{meta['tax_id']}`")
        if dettagli:
            st.markdown(" &nbsp;•&nbsp; ".join(dettagli))
    else:
        st.info(T["info_no_header"])

    # >>> calcolo pesante nel core C++ <<<
    r = seq_core.analizza(pulita)

    # Popover 'Dati chiave' (icona a sinistra, si apre al passaggio del mouse).
    pannello_dati_chiave(T, meta, r)

    st.subheader(T["res_header"])
    # Metriche principali come griglia di card 'dashboard' (3×2).
    griglia_metriche([
        {"icon": "📏", "color": "#3b7dd8", "label": T["m_length"],
         "value": f"{r.lunghezza}", "sub": T["unit_aa"]},
        {"icon": "⚖️", "color": "#7e57c2", "label": T["m_mw"],
         "value": f"{r.peso_molecolare/1000:.1f} kDa",
         "sub": f"{r.peso_molecolare:,.0f} Da"},
        {"icon": "🎯", "color": "#c0504d", "label": T["m_pi"],
         "value": f"{r.punto_isoelettrico:.2f}", "sub": "pH"},
        {"icon": "⚡", "color": "#e8a33d", "label": T["m_charge"],
         "value": f"{r.carica_a_pH7:+.1f}", "sub": "pH 7"},
        {"icon": "🔗", "color": "#2a9d8f", "label": T["m_bonds"],
         "value": f"{r.legami_peptidici}", "sub": ""},
        {"icon": "🧬", "color": "#4e9a5f", "label": T["m_disulfide"],
         "value": T["disulfide_val"].format(n=r.ponti_disolfuro_max),
         "sub": T["cys_delta"].format(n=r.cisteine)},
    ], ncol=3)

    # --- Proprietà biochimiche aggiuntive (stile ProtParam, calcolo in C++) ---
    idro = T["gravy_hydrophobic"] if r.gravy > 0 else T["gravy_hydrophilic"]
    griglia_metriche([
        {"icon": "💧", "color": "#17a2b8", "label": T["m_gravy"],
         "value": f"{r.gravy:+.3f}", "sub": idro},
        {"icon": "🌡️", "color": "#d9822b", "label": T["m_aliphatic"],
         "value": f"{r.indice_alifatico:.1f}", "sub": T["aliphatic_delta"]},
        {"icon": "🔦", "color": "#5c6bc0", "label": T["m_abs"],
         "value": f"{r.abs280_ox:.2f}",
         "sub": T["abs_delta"].format(v=f"{r.estinzione_ox:,.0f}")},
    ], ncol=3)
    with st.expander(T["exp_props"]):
        st.markdown(T["props_md"].format(
            gravy=f"{r.gravy:+.3f}", idro=idro, ali=f"{r.indice_alifatico:.1f}",
            est_red=f"{r.estinzione_red:,.0f}", abs_red=f"{r.abs280_red:.3f}",
            est_ox=f"{r.estinzione_ox:,.0f}", abs_ox=f"{r.abs280_ox:.3f}"))

    # --- Curva di titolazione (carica netta vs pH) — dati dal core C++ ---
    fig_tit = px.line(x=list(r.titolazione_pH), y=list(r.titolazione_carica),
                      labels={"x": T["ax_ph"], "y": T["ax_charge"]},
                      title=T["r_titration"])
    fig_tit.update_traces(line_color="#3b7dd8")
    fig_tit.add_hline(y=0, line_dash="dot", line_color="#888")
    fig_tit.add_vline(x=r.punto_isoelettrico, line_dash="dash", line_color="#c0504d",
                      annotation_text=T["titration_pi"].format(v=f"{r.punto_isoelettrico:.2f}"),
                      annotation_position="top")
    with st.expander(T["titration_header"]):
        st.plotly_chart(fig_tit)

    # --- Struttura secondaria (Chou-Fasman, conteggi dal core C++) ---
    ss_tot = r.ss_helix + r.ss_sheet + r.ss_coil
    with st.expander(T["ss_header"]):
        if ss_tot > 0:
            s1, s2, s3 = st.columns(3)
            s1.metric(T["ss_helix"], f"{100*r.ss_helix/ss_tot:.0f}%",
                      f"{r.ss_helix} {T['unit_aa']}", delta_color="off")
            s2.metric(T["ss_sheet"], f"{100*r.ss_sheet/ss_tot:.0f}%",
                      f"{r.ss_sheet} {T['unit_aa']}", delta_color="off")
            s3.metric(T["ss_coil"], f"{100*r.ss_coil/ss_tot:.0f}%",
                      f"{r.ss_coil} {T['unit_aa']}", delta_color="off")
            fig_ss = px.pie(names=[T["ss_helix"], T["ss_sheet"], T["ss_coil"]],
                            values=[r.ss_helix, r.ss_sheet, r.ss_coil],
                            title=T["ss_chart_title"], hole=0.4,
                            color_discrete_sequence=["#e15759", "#4e79a7", "#bab0ac"])
            fig_ss.update_traces(textinfo="label+percent", sort=False)
            tema_ss = _tema_grafico()
            fig_ss.update_layout(paper_bgcolor=tema_ss["sfondo"],
                                 font_color=tema_ss["testo"], template=tema_ss["template"])
            st.plotly_chart(fig_ss)
        st.caption(T["ss_note"])

    # DataFrame delle frequenze con colonne nella lingua scelta.
    cc, cn, cq, cp = T["col_code"], T["col_name"], T["col_count"], T["col_pct"]
    conteggi = {chiave(k): v for k, v in r.conteggi.items()}
    df = pd.DataFrame(
        [(aa, NOMI_AA[aa], conteggi.get(aa, 0)) for aa in sorted(AA_VALIDI)],
        columns=[cc, cn, cq],
    )
    df[cp] = (100 * df[cq] / df[cq].sum()).round(2)
    ch = T["hydro_legend"]
    df[ch] = df[cc].map(HYDRO_KD)         # idrofobicità Kyte-Doolittle per residuo
    classi = {chiave(k): v for k, v in r.classi.items() if v > 0}

    d = df.sort_values(cq, ascending=False)
    fig = px.bar(d, x=cc, y=cq, title=T["chart_freq_title"],
                 hover_data=[cn, cp, ch])
    # Barre colorate per idrofobicità (Kyte-Doolittle): stessa scala divergente
    # della torta delle classi. Range simmetrico, così il neutro (KD≈0) sta al centro.
    fig.update_traces(
        marker=dict(color=d[ch], colorscale=HYDRO_SCALE, cmin=-4.5, cmax=4.5,
                    showscale=True,
                    colorbar=dict(title=dict(text=ch, side="right"),
                                  thickness=12, len=0.85, outlinewidth=0)),
        hovertemplate="<b>%{customdata[0]} (%{x})</b><br>"
                      + T["hv_count"] + ": %{y}<br>"
                      + T["hv_pct"] + ": %{customdata[1]}%<br>"
                      + ch + ": %{customdata[2]:+.1f}<extra></extra>",
    )
    fig.update_layout(hoverlabel=dict(font_size=14, font_family="ui-monospace"))

    # Classi: etichette nella lingua scelta, colori per indice (chiave canonica).
    classi_nomi = [CLASSI_LABEL.get(k, {}).get(lang, k) for k in classi.keys()]
    classi_colori = [COLORI_CLASSI.get(k, "#888888") for k in classi.keys()]
    color_map = {nome: col for nome, col in zip(classi_nomi, classi_colori)}
    tema = _tema_grafico()
    fig2 = px.pie(names=classi_nomi, values=list(classi.values()),
                  title=T["chart_classes_title"], hole=0.4,
                  color=classi_nomi, color_discrete_map=color_map)
    # fette separate da un sottile spazio del colore di sfondo + tooltip prominente
    fig2.update_traces(
        textinfo="label+percent",
        marker=dict(line=dict(color=tema["sfondo"], width=2)),
        pull=[0.03] * len(classi),
        sort=False,
        hoverlabel=dict(font_size=15),
        hovertemplate="<b>%{label}</b><br>%{value} " + T["hv_aa"]
                      + " (%{percent})<extra></extra>",
    )

    g1, g2 = st.columns([3, 2])
    with g1:
        grafico_con_hover(fig, "bar", "graf_freq")
    with g2:
        grafico_con_hover(fig2, "pie", "graf_classi", colori=classi_colori)

    # --- Legenda delle classi chimico-fisiche (consultazione rapida) ---
    with st.expander(T["exp_grouping"]):
        st.markdown(T["grouping_md"])

    with st.expander(T["exp_freq_table"]):
        st.dataframe(d, hide_index=True)
        base = (meta.get("accession") if meta else None) or T["default_name"]
        st.download_button(
            T["btn_csv"],
            data=d.to_csv(index=False).encode("utf-8"),
            file_name=f"{T['csv_prefix']}{base}.csv",
            mime="text/csv",
        )

    # --- RILEVATORE DI DOMINI DI MEMBRANA (calcolo nel core C++) ---
    st.subheader(T["mem_header"])

    # Controlli: scala del profilo + finestra della media mobile.
    ctl1, ctl2 = st.columns([1, 1])
    scala_key = ctl1.selectbox(
        T["scale_label"], list(SCALE_IDRO.keys()),
        format_func=lambda k: SCALE_IDRO[k][0], key="idro_scala")
    finestra = ctl2.select_slider(
        T["window_label"], options=[7, 9, 11, 13, 15, 17, 19, 21, 23, 25],
        value=19, key="idro_finestra")

    # Il rilevamento TM (verdetto/segmenti/heatmap) resta su Kyte-Doolittle nel
    # core C++, con la finestra scelta dall'utente.
    prof = seq_core.idrofobicita(pulita, finestra=finestra)
    verdetto = VERDETTO_LABEL.get(prof.verdetto, {}).get(lang, prof.verdetto)

    v1, v2 = st.columns([2, 1])
    if prof.n_domini == 0:
        v1.success(f"**{verdetto}**")
    else:
        v1.info(f"**{verdetto}**")
    v2.metric(T["m_tm"], prof.n_domini)

    # Curva del profilo: calcolata in Python con la scala scelta (così posso
    # offrire scale alternative oltre a quella KD del core C++).
    nome_scala, tabella_scala, scala_idrofoba = SCALE_IDRO[scala_key]
    pos_idro, sco_idro = profilo_finestra(pulita, tabella_scala, finestra)
    titolo_prof = f"{nome_scala} — {T['window_label']} {finestra}"

    tema_idro = _tema_grafico()
    # linea leggibile su entrambi i temi: azzurro chiaro sullo scuro, navy sul chiaro
    linea_idro = "#9ec1f0" if tema_idro["template"] == "plotly_dark" else "#1d2733"
    fig_idro = px.line(x=pos_idro, y=sco_idro,
                       labels={"x": T["ax_pos"], "y": nome_scala},
                       title=titolo_prof)
    fig_idro.update_traces(line_color=linea_idro, line_width=2)
    if scala_key == "kd":                 # soglia TM solo per Kyte-Doolittle
        fig_idro.add_hline(y=1.6, line_dash="dash", line_color="#c0504d",
                           annotation_text=T["tm_threshold"])
    else:                                 # linea di riferimento allo zero
        fig_idro.add_hline(y=0, line_dash="dot", line_color="#888")
    for i, s in enumerate(prof.segmenti, 1):
        fig_idro.add_vrect(x0=s.inizio, x1=s.fine, fillcolor="#2b6cb0",
                           opacity=0.18, line_width=0,
                           annotation_text=f"TM{i}", annotation_position="top")

    # Profilo 2D + struttura 3D AlphaFold. Se c'è l'accession e il modello, uso
    # il widget COMBINATO con cross-highlighting; altrimenti solo il profilo.
    acc_3d = (meta.get("accession") if meta else None)
    af = None
    if acc_3d:
        with st.spinner(T["struct3d_spinner"]):
            af = fetch_alphafold(acc_3d)

    if af:
        soglia_w = 1.6 if scala_key == "kd" else None
        visualizzatore_combinato(
            pos_idro, sco_idro, prof.segmenti, af["pdb"],
            asse_x=T["ax_pos"], asse_y=nome_scala, titolo=titolo_prof,
            titolo_3d=T["struct3d_header"], soglia=soglia_w, linea=linea_idro,
            scuro=(tema_idro["template"] == "plotly_dark"), altezza=400)
        st.caption(T["crosshl_hint"])
        note = [T["struct3d_plddt"], T["struct3d_source"].format(id=af.get("id", ""))]
        if prof.segmenti and scala_key != "kd":
            note.append(T["scale_tm_note"])
        st.caption("  ·  ".join(note))
    else:
        col_prof, col_3d = st.columns([3, 2])
        with col_prof:
            st.plotly_chart(fig_idro, width="stretch")
            if prof.segmenti and scala_key != "kd":   # TM = rilevate con KD
                st.caption(T["scale_tm_note"])
        with col_3d:
            st.markdown(f"**{T['struct3d_header']}**")
            if acc_3d:
                st.info(T["struct3d_notfound"].format(acc=acc_3d))
            else:
                st.info(T["struct3d_need_acc"])

    heat = heatmap_sequenza(pulita, prof.segmenti)
    with st.expander(T["exp_tm_detail"], expanded=True):
        if prof.segmenti:
            tm_df = pd.DataFrame(
                [(f"TM{i}", f"{s.inizio}–{s.fine}", f"{s.lunghezza} {T['unit_aa']}",
                  round(s.idro_media, 2), s.sottosequenza)
                 for i, s in enumerate(prof.segmenti, 1)],
                columns=[T["tm_seg"], T["tm_res"], T["tm_len"], T["tm_avg"], T["tm_sub"]],
            )
            st.dataframe(tm_df, hide_index=True)
        else:
            st.write(T["no_tm"])
        st.markdown(T["legend_seq"], unsafe_allow_html=True)
        st.markdown(heat, unsafe_allow_html=True)

    # --- DOWNLOAD del report ---
    nome_file = (meta.get("accession") if meta else None) or "report"
    report = costruisci_report_html(meta, r, prof, fig, fig2, fig_idro, fig_tit, heat, T, lang)
    st.download_button(T["btn_report"], data=report,
                       file_name=f"{T['report_prefix']}{nome_file}.html", mime="text/html")

    # --- ESPORTAZIONE VETTORIALE DEI GRAFICI (SVG / PDF, per tesi/relazioni) ---
    with st.expander(T["exp_export_header"]):
        st.caption(T["exp_export_hint"])
        fmt = st.radio(T["exp_export_format"], ["SVG", "PDF"],
                       horizontal=True, key="export_fmt")
        figure_export = [
            (T["chart_freq_title"], "frequenza_aa", _fig_per_stampa(fig)),
            (T["chart_classes_title"], "classi_chimiche", _fig_per_stampa(fig2)),
            (titolo_prof, "profilo_idrofobicita",
             _fig_per_stampa(fig_idro, line_scuro=True)),
            (T["r_titration"], "titolazione", _fig_per_stampa(fig_tit)),
        ]
        if ss_tot > 0:
            figure_export.append(
                (T["ss_chart_title"], "struttura_secondaria", _fig_per_stampa(fig_ss)))

        seq_key = f"{len(pulita)}:{pulita[:24]}"   # per invalidare export di altre sequenze
        if st.button(T["exp_export_prepare"]):
            ext = fmt.lower()
            preparati = {}
            try:
                with st.spinner(T["exp_export_spinner"]):
                    for label, base, f in figure_export:
                        preparati[base] = (label, f.to_image(format=ext, scale=2))
                st.session_state["export_files"] = {"ext": ext, "seq": seq_key,
                                                     "items": preparati}
            except Exception:
                st.session_state.pop("export_files", None)
                st.error(T["exp_export_error"])

        exp = st.session_state.get("export_files")
        if exp and exp.get("seq") == seq_key and exp.get("ext") == fmt.lower():
            mime = "image/svg+xml" if exp["ext"] == "svg" else "application/pdf"
            st.markdown(f"**{T['exp_export_ready']}**")
            for base, (label, data) in exp["items"].items():
                st.download_button(f"⬇️ {label} (.{exp['ext']})", data=data,
                                   file_name=f"{base}.{exp['ext']}", mime=mime,
                                   key=f"dl_{base}_{exp['ext']}")

    st.caption(T["caption_pi"])

    # --- Confronto con una seconda proteina ---
    with st.expander(T["compare_header"]):
        st.text_area(T["compare_label"], key="seq_b_input",
                     placeholder=T["compare_ph"], height=110)
        if st.button(T["btn_compare"]):
            st.session_state.compare_seq = st.session_state.get("seq_b_input", "")
        seq_b = st.session_state.get("compare_seq", "")
        if seq_b:
            meta_b = parse_fasta_header(seq_b, lang)
            pulita_b, _ = valida(seq_b)
            if not pulita_b:
                st.error(T["compare_err"])
            else:
                rb = seq_core.analizza(pulita_b)
                nome_a = ((meta.get("proteina") or meta.get("accession")) if meta else None) or T["compare_a"]
                nome_b = ((meta_b.get("proteina") or meta_b.get("accession")) if meta_b else None) or T["compare_b"]
                nome_a, nome_b = str(nome_a)[:28], str(nome_b)[:28]
                ssa = (r.ss_helix + r.ss_sheet + r.ss_coil) or 1
                ssb = (rb.ss_helix + rb.ss_sheet + rb.ss_coil) or 1
                righe = [
                    (T["m_length"], f"{r.lunghezza}", f"{rb.lunghezza}"),
                    (T["m_mw"], f"{r.peso_molecolare/1000:.1f} kDa", f"{rb.peso_molecolare/1000:.1f} kDa"),
                    (T["m_pi"], f"{r.punto_isoelettrico:.2f}", f"{rb.punto_isoelettrico:.2f}"),
                    (T["m_charge"], f"{r.carica_a_pH7:+.1f}", f"{rb.carica_a_pH7:+.1f}"),
                    (T["m_gravy"], f"{r.gravy:+.3f}", f"{rb.gravy:+.3f}"),
                    (T["m_aliphatic"], f"{r.indice_alifatico:.1f}", f"{rb.indice_alifatico:.1f}"),
                    (T["m_disulfide"], f"{r.ponti_disolfuro_max}", f"{rb.ponti_disolfuro_max}"),
                    (T["ss_helix"], f"{100*r.ss_helix/ssa:.0f}%", f"{100*rb.ss_helix/ssb:.0f}%"),
                    (T["ss_sheet"], f"{100*r.ss_sheet/ssa:.0f}%", f"{100*rb.ss_sheet/ssb:.0f}%"),
                ]
                df_cmp = pd.DataFrame(righe, columns=[T["compare_metric"], nome_a, nome_b])
                st.dataframe(df_cmp, hide_index=True)

                # composizione in classi a confronto (% dei residui)
                ca = {chiave(k): v for k, v in r.classi.items()}
                cb = {chiave(k): v for k, v in rb.classi.items()}
                ta = sum(ca.values()) or 1
                tb = sum(cb.values()) or 1
                rows = []
                for k in CLASSI_LABEL:
                    lab = CLASSI_LABEL[k].get(lang, k)
                    rows.append({"c": lab, "s": nome_a, "p": 100 * ca.get(k, 0) / ta})
                    rows.append({"c": lab, "s": nome_b, "p": 100 * cb.get(k, 0) / tb})
                fig_cmp = px.bar(pd.DataFrame(rows), x="c", y="p", color="s", barmode="group",
                                 title=T["compare_chart_title"],
                                 labels={"p": T["compare_y"], "c": "", "s": ""},
                                 color_discrete_sequence=["#4e79a7", "#f28e2b"])
                tema_c = _tema_grafico()
                fig_cmp.update_layout(paper_bgcolor=tema_c["sfondo"], plot_bgcolor=tema_c["sfondo"],
                                      font_color=tema_c["testo"], template=tema_c["template"],
                                      legend_title_text="")
                st.plotly_chart(fig_cmp)
