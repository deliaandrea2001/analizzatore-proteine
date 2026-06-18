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
import sys
import subprocess
import urllib.request
import urllib.error
import pandas as pd
import plotly.express as px
import streamlit as st


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

# Colori delle classi (chiave = nome canonico). Servono al JS per attenuare le fette.
COLORI_CLASSI = {
    "Idrofobici": "#9ecae1",        # azzurro chiaro
    "Polari": "#2b6cb0",            # blu
    "Carichi positivi": "#fcae91",  # rosa/salmone
    "Carichi negativi": "#ef3b2c",  # rosso
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
        "ax_pos": "position (residue)", "ax_hydro": "average hydrophobicity",
        "tm_threshold": "TM threshold (1.6)",
        "exp_tm_detail": "🔍 Hydrophobic regions detail and sequence map",
        "tm_seg": "Segment", "tm_res": "Residues", "tm_len": "Length",
        "tm_avg": "Avg. hydrophobicity", "tm_sub": "Subsequence",
        "no_tm": "No transmembrane segment detected.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>blue</span> "
                       "= transmembrane residues &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>green</span> = hydrophilic residues (water-exposed)"),
        "btn_report": "⬇️ Download report (HTML)", "report_prefix": "report_",
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
        "titration_header": "📈 Titration curve (net charge vs pH)",
        "ax_ph": "pH", "ax_charge": "net charge", "titration_pi": "pI = {v}",
        "r_titration": "Titration curve (net charge vs pH)",
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
        "ax_pos": "posizione (residuo)", "ax_hydro": "idrofobicità media",
        "tm_threshold": "soglia TM (1.6)",
        "exp_tm_detail": "🔍 Dettaglio regioni idrofobiche e mappa della sequenza",
        "tm_seg": "Segmento", "tm_res": "Residui", "tm_len": "Lunghezza",
        "tm_avg": "Idrofob. media", "tm_sub": "Sottosequenza",
        "no_tm": "Nessun segmento transmembrana rilevato.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>azzurro</span> "
                       "= residui transmembrana &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>verde</span> = residui idrofilici (esposti all'acqua)"),
        "btn_report": "⬇️ Scarica il report (HTML)", "report_prefix": "report_",
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
        "titration_header": "📈 Curva di titolazione (carica netta vs pH)",
        "ax_ph": "pH", "ax_charge": "carica netta", "titration_pi": "pI = {v}",
        "r_titration": "Curva di titolazione (carica netta vs pH)",
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
        "ax_pos": "posición (residuo)", "ax_hydro": "hidrofobicidad media",
        "tm_threshold": "umbral TM (1.6)",
        "exp_tm_detail": "🔍 Detalle de regiones hidrofóbicas y mapa de la secuencia",
        "tm_seg": "Segmento", "tm_res": "Residuos", "tm_len": "Longitud",
        "tm_avg": "Hidrofob. media", "tm_sub": "Subsecuencia",
        "no_tm": "No se detectó ningún segmento transmembrana.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>azul</span> "
                       "= residuos transmembrana &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>verde</span> = residuos hidrofílicos (expuestos al agua)"),
        "btn_report": "⬇️ Descargar el informe (HTML)", "report_prefix": "informe_",
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
        "titration_header": "📈 Curva de titulación (carga neta vs pH)",
        "ax_ph": "pH", "ax_charge": "carga neta", "titration_pi": "pI = {v}",
        "r_titration": "Curva de titulación (carga neta vs pH)",
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
        "ax_pos": "Position (Rest)", "ax_hydro": "mittlere Hydrophobizität",
        "tm_threshold": "TM-Schwelle (1.6)",
        "exp_tm_detail": "🔍 Details zu hydrophoben Regionen und Sequenzkarte",
        "tm_seg": "Segment", "tm_res": "Reste", "tm_len": "Länge",
        "tm_avg": "Mittl. Hydrophob.", "tm_sub": "Teilsequenz",
        "no_tm": "Kein Transmembransegment erkannt.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>blau</span> "
                       "= Transmembranreste &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>grün</span> = hydrophile Reste (wasserexponiert)"),
        "btn_report": "⬇️ Bericht herunterladen (HTML)", "report_prefix": "bericht_",
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
        "titration_header": "📈 Titrationskurve (Nettoladung vs. pH)",
        "ax_ph": "pH", "ax_charge": "Nettoladung", "titration_pi": "pI = {v}",
        "r_titration": "Titrationskurve (Nettoladung vs. pH)",
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
        "ax_pos": "position (résidu)", "ax_hydro": "hydrophobicité moyenne",
        "tm_threshold": "seuil TM (1.6)",
        "exp_tm_detail": "🔍 Détail des régions hydrophobes et carte de la séquence",
        "tm_seg": "Segment", "tm_res": "Résidus", "tm_len": "Longueur",
        "tm_avg": "Hydrophob. moy.", "tm_sub": "Sous-séquence",
        "no_tm": "Aucun segment transmembranaire détecté.",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>bleu</span> "
                       "= résidus transmembranaires &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>vert</span> = résidus hydrophiles (exposés à l'eau)"),
        "btn_report": "⬇️ Télécharger le rapport (HTML)", "report_prefix": "rapport_",
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
        "titration_header": "📈 Courbe de titration (charge nette vs pH)",
        "ax_ph": "pH", "ax_charge": "charge nette", "titration_pi": "pI = {v}",
        "r_titration": "Courbe de titration (charge nette vs pH)",
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
        "ax_pos": "位置（残基）", "ax_hydro": "平均疏水性",
        "tm_threshold": "TM 阈值 (1.6)",
        "exp_tm_detail": "🔍 疏水区域详情与序列图",
        "tm_seg": "片段", "tm_res": "残基", "tm_len": "长度",
        "tm_avg": "平均疏水性", "tm_sub": "子序列",
        "no_tm": "未检测到跨膜片段。",
        "legend_seq": ("<span style='background:#2b6cb0;color:#fff;padding:1px 6px;border-radius:3px'>蓝色</span> "
                       "= 跨膜残基 &nbsp;&nbsp; "
                       "<span style='color:#2f855a'>绿色</span> = 亲水残基（暴露于水）"),
        "btn_report": "⬇️ 下载报告 (HTML)", "report_prefix": "report_",
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
        "titration_header": "📈 滴定曲线（净电荷 vs pH）",
        "ax_ph": "pH", "ax_charge": "净电荷", "titration_pi": "pI = {v}",
        "r_titration": "滴定曲线（净电荷 vs pH）",
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
    """Scarica il FASTA di un accession da UniProt. None se fallisce."""
    acc = acc.strip()
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "analizzatore-proteine"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            testo = resp.read().decode("utf-8", "ignore")
        return testo if testo.lstrip().startswith(">") else None
    except Exception:
        return None


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
    tipo = None
    try:
        tipo = st.context.theme.type        # Streamlit recenti: tema attivo reale
    except Exception:
        tipo = st.get_option("theme.base")  # fallback
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
                        config={"displayModeBar": False, "responsive": True},
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


# ===================== INTERFACCIA =====================
st.set_page_config(page_title="Protein sequence analyzer",
                   page_icon="🧬", layout="wide")

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

src1, src2 = st.columns(2)
with src1:
    up = st.file_uploader(T["upload_label"], type=["fasta", "fa", "faa", "txt", "seq"])
    if up is not None:
        fid = (up.name, up.size)
        if st.session_state.get("_ultimo_file") != fid:   # carica solo i file nuovi
            st.session_state._ultimo_file = fid
            _imposta_sequenza(up.getvalue().decode("utf-8", "ignore"))
with src2:
    acc = st.text_input(T["accession_label"], placeholder=T["accession_ph"], key="acc_input")
    if st.button(T["btn_fetch"]) and acc.strip():
        with st.spinner(T["fetch_spinner"]):
            scaricato = fetch_uniprot(acc)
        if scaricato:
            _imposta_sequenza(scaricato)
        else:
            st.error(T["fetch_error"].format(acc=acc.strip()))

seq_input = st.text_area(
    T["input_label"], height=130, key="seq_input",
    placeholder=T["input_placeholder"],
)

analizza = st.button(T["btn_analyze"], type="primary")
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

    st.subheader(T["res_header"])
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric(T["m_length"], f"{r.lunghezza}", T["unit_aa"])
    c2.metric(T["m_mw"], f"{r.peso_molecolare/1000:.1f} kDa")
    c3.metric(T["m_pi"], f"{r.punto_isoelettrico:.2f}")
    c4.metric(T["m_charge"], f"{r.carica_a_pH7:+.1f}")
    c5.metric(T["m_bonds"], f"{r.legami_peptidici}")
    c6.metric(T["m_disulfide"], T["disulfide_val"].format(n=r.ponti_disolfuro_max),
              T["cys_delta"].format(n=r.cisteine))

    # --- Proprietà biochimiche aggiuntive (stile ProtParam, calcolo in C++) ---
    b1, b2, b3 = st.columns(3)
    idro = T["gravy_hydrophobic"] if r.gravy > 0 else T["gravy_hydrophilic"]
    b1.metric(T["m_gravy"], f"{r.gravy:+.3f}", idro, delta_color="off")
    b2.metric(T["m_aliphatic"], f"{r.indice_alifatico:.1f}",
              T["aliphatic_delta"], delta_color="off")
    b3.metric(T["m_abs"], f"{r.abs280_ox:.2f}",
              T["abs_delta"].format(v=f"{r.estinzione_ox:,.0f}"), delta_color="off")
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

    # DataFrame delle frequenze con colonne nella lingua scelta.
    cc, cn, cq, cp = T["col_code"], T["col_name"], T["col_count"], T["col_pct"]
    conteggi = {chiave(k): v for k, v in r.conteggi.items()}
    df = pd.DataFrame(
        [(aa, NOMI_AA[aa], conteggi.get(aa, 0)) for aa in sorted(AA_VALIDI)],
        columns=[cc, cn, cq],
    )
    df[cp] = (100 * df[cq] / df[cq].sum()).round(2)
    classi = {chiave(k): v for k, v in r.classi.items() if v > 0}

    d = df.sort_values(cq, ascending=False)
    fig = px.bar(d, x=cc, y=cq, color=cq,
                 color_continuous_scale="Blues", title=T["chart_freq_title"],
                 hover_data=[cn, cp])
    fig.update_layout(coloraxis_showscale=False,
                      hoverlabel=dict(font_size=14, font_family="ui-monospace"))
    # tooltip prominente sulla barra sotto il cursore
    fig.update_traces(
        hovertemplate="<b>%{customdata[0]} (%{x})</b><br>"
                      + T["hv_count"] + ": %{y}<br>"
                      + T["hv_pct"] + ": %{customdata[1]}%<extra></extra>",
    )

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
    prof = seq_core.idrofobicita(pulita)
    verdetto = VERDETTO_LABEL.get(prof.verdetto, {}).get(lang, prof.verdetto)
    st.subheader(T["mem_header"])

    v1, v2 = st.columns([2, 1])
    if prof.n_domini == 0:
        v1.success(f"**{verdetto}**")
    else:
        v1.info(f"**{verdetto}**")
    v2.metric(T["m_tm"], prof.n_domini)

    # grafico del profilo con soglia e regioni TM evidenziate
    fig_idro = px.line(x=list(prof.posizioni), y=list(prof.punteggi),
                       labels={"x": T["ax_pos"], "y": T["ax_hydro"]},
                       title=T["profile_title"])
    fig_idro.update_traces(line_color="#1d2733")
    fig_idro.add_hline(y=1.6, line_dash="dash", line_color="#c0504d",
                       annotation_text=T["tm_threshold"])
    for i, s in enumerate(prof.segmenti, 1):
        fig_idro.add_vrect(x0=s.inizio, x1=s.fine, fillcolor="#2b6cb0",
                           opacity=0.18, line_width=0,
                           annotation_text=f"TM{i}", annotation_position="top")
    st.plotly_chart(fig_idro)

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

    st.caption(T["caption_pi"])
