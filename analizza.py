#!/usr/bin/env python3
# =====================================================================
#  analizza.py  -  Parte Python dell'analizzatore di sequenze proteiche.
# ---------------------------------------------------------------------
#  Progetto IBRIDO C++ + Python:
#    * C++  (seq_core) -> calcolo: conteggi, peso molecolare, punto
#                         isoelettrico (bisezione), carica netta.
#    * Python (qui)    -> validazione input, Pandas per le percentuali,
#                         grafici matplotlib e report HTML automatico.
#
#  Compila prima il core:   python setup.py build_ext --inplace
#  Poi lancia:
#      python analizza.py KCNQ1.fasta          # da file FASTA
#      python analizza.py MKTAYIAKQR...        # sequenza diretta
#      python analizza.py                      # usa KCNQ1.fasta se presente
# =====================================================================

import sys
import os
import io
import base64

import pandas as pd
import matplotlib
matplotlib.use("Agg")              # backend senza finestre (solo file)
import matplotlib.pyplot as plt

import seq_core                    # <-- il nostro core C++

AA_VALIDI = set("ACDEFGHIKLMNPQRSTVWY")
NOMI_AA = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
}


def leggi_sequenza(arg: str) -> tuple[str, str]:
    """Restituisce (nome, sequenza) da un file FASTA o da una stringa."""
    if os.path.isfile(arg):
        nome, righe = os.path.basename(arg), []
        with open(arg) as f:
            for r in f:
                if r.startswith(">"):
                    nome = r[1:].strip()
                else:
                    righe.append(r.strip())
        return nome, "".join(righe)
    return "sequenza inserita", arg


def valida(seq: str) -> str:
    """Pulisce la sequenza tenendo solo i 20 amminoacidi standard."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    pulita = "".join(c for c in seq if c in AA_VALIDI)
    scartati = len(seq) - len(pulita)
    if scartati:
        print(f"  Attenzione: {scartati} caratteri non standard ignorati.")
    if not pulita:
        sys.exit("Errore: nessun amminoacido valido nella sequenza.")
    return pulita


def chiave(k) -> str:
    """Le chiavi dei dict da C++ possono arrivare come str o int: normalizzo."""
    return k if isinstance(k, str) else chr(k)


def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def grafico_frequenze(df: pd.DataFrame) -> str:
    d = df.sort_values("conteggio", ascending=False)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(d["amminoacido"], d["conteggio"], color="#3b7dd8")
    ax.set_ylabel("conteggio")
    ax.set_title("Frequenza degli amminoacidi")
    ax.tick_params(axis="x", rotation=0)
    return fig_to_base64(fig)


def grafico_classi(classi: dict) -> str:
    etichette = [k for k, v in classi.items() if v > 0]
    valori = [classi[k] for k in etichette]
    colori = ["#e07b39", "#3b7dd8", "#4caf50", "#c0504d"]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(valori, labels=etichette, autopct="%1.1f%%",
           colors=colori[:len(etichette)], startangle=90)
    ax.set_title("Classi chimico-fisiche")
    return fig_to_base64(fig)


def kpi_html(titolo, valore, sotto=""):
    return f"""
    <div class="kpi">
      <div class="kpi-val">{valore}</div>
      <div class="kpi-tit">{titolo}</div>
      <div class="kpi-sub">{sotto}</div>
    </div>"""


def genera_report(nome, seq, r, df, b64_freq, b64_classi, out_path):
    righe_tab = "".join(
        f"<tr><td>{row.amminoacido}</td><td>{NOMI_AA[row.amminoacido]}</td>"
        f"<td>{row.conteggio}</td><td>{row.percentuale:.2f}%</td></tr>"
        for row in df.sort_values("conteggio", ascending=False).itertuples()
    )
    html = f"""<!DOCTYPE html>
<html lang="it"><head><meta charset="utf-8">
<title>Analisi proteica - {nome}</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, sans-serif;
         margin: 0; background: #f4f6fa; color: #1d2733; }}
  header {{ background: #1d2733; color: #fff; padding: 26px 40px; }}
  header h1 {{ margin: 0 0 6px; font-size: 22px; }}
  header p {{ margin: 0; color: #9fb0c3; font-size: 13px; }}
  main {{ max-width: 980px; margin: 0 auto; padding: 28px 40px 60px; }}
  .kpis {{ display: flex; flex-wrap: wrap; gap: 14px; margin: 20px 0 30px; }}
  .kpi {{ background: #fff; border-radius: 12px; padding: 18px 22px;
         flex: 1 1 150px; box-shadow: 0 1px 4px rgba(0,0,0,.07); }}
  .kpi-val {{ font-size: 26px; font-weight: 700; color: #3b7dd8; }}
  .kpi-tit {{ font-size: 13px; font-weight: 600; margin-top: 4px; }}
  .kpi-sub {{ font-size: 11px; color: #7a8aa0; }}
  .card {{ background: #fff; border-radius: 12px; padding: 22px;
          margin-bottom: 22px; box-shadow: 0 1px 4px rgba(0,0,0,.07); }}
  .card h2 {{ margin-top: 0; font-size: 16px; }}
  img {{ max-width: 100%; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th, td {{ text-align: left; padding: 6px 10px; border-bottom: 1px solid #eef1f5; }}
  th {{ color: #7a8aa0; font-weight: 600; }}
  .seq {{ font-family: ui-monospace, Menlo, monospace; font-size: 11px;
         word-break: break-all; color: #54657a; line-height: 1.7; }}
  footer {{ text-align: center; color: #9fb0c3; font-size: 11px; padding: 20px; }}
</style></head><body>
<header>
  <h1>Analisi di sequenza proteica</h1>
  <p>{nome}</p>
</header>
<main>
  <div class="kpis">
    {kpi_html("Lunghezza", r.lunghezza, "amminoacidi")}
    {kpi_html("Peso molecolare", f"{r.peso_molecolare/1000:.1f} kDa", f"{r.peso_molecolare:,.0f} Da")}
    {kpi_html("Punto isoelettrico", f"{r.punto_isoelettrico:.2f}", "pI (bisezione, C++)")}
    {kpi_html("Carica a pH 7", f"{r.carica_a_pH7:+.1f}", "netta")}
    {kpi_html("Legami peptidici", r.legami_peptidici, "L - 1")}
    {kpi_html("Ponti disolfuro", f"max {r.ponti_disolfuro_max}", f"{r.cisteine} cisteine")}
  </div>

  <div class="card">
    <h2>Frequenza degli amminoacidi</h2>
    <img src="data:image/png;base64,{b64_freq}">
  </div>

  <div class="card">
    <h2>Composizione chimico-fisica</h2>
    <img src="data:image/png;base64,{b64_classi}">
  </div>

  <div class="card">
    <h2>Tabella delle frequenze</h2>
    <table>
      <tr><th>Codice</th><th>Nome</th><th>Conteggio</th><th>Percentuale</th></tr>
      {righe_tab}
    </table>
  </div>

  <div class="card">
    <h2>Sequenza analizzata</h2>
    <div class="seq">{seq}</div>
  </div>
</main>
<footer>
  Report generato da un'architettura ibrida: core di calcolo in C++ (pybind11),
  analisi e visualizzazione in Python (Pandas + matplotlib).
</footer>
</body></html>"""
    with open(out_path, "w") as f:
        f.write(html)


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else "KCNQ1.fasta"
    nome, grezza = leggi_sequenza(arg)
    print(f"Analizzo: {nome}")
    seq = valida(grezza)

    # --- calcolo pesante: tutto in C++ ---
    r = seq_core.analizza(seq)

    # --- analisi comoda in Python con Pandas ---
    conteggi = {chiave(k): v for k, v in r.conteggi.items()}
    df = pd.DataFrame(
        [(aa, conteggi.get(aa, 0)) for aa in sorted(AA_VALIDI)],
        columns=["amminoacido", "conteggio"],
    )
    df["percentuale"] = 100 * df["conteggio"] / df["conteggio"].sum()
    classi = {chiave(k): v for k, v in r.classi.items()}

    # --- riepilogo a schermo ---
    print(f"  Lunghezza:          {r.lunghezza} aa")
    print(f"  Peso molecolare:    {r.peso_molecolare:,.1f} Da ({r.peso_molecolare/1000:.1f} kDa)")
    print(f"  Punto isoelettrico: {r.punto_isoelettrico:.2f}")
    print(f"  Carica a pH 7:      {r.carica_a_pH7:+.2f}")
    print(f"  Cisteine:           {r.cisteine} (max {r.ponti_disolfuro_max} ponti disolfuro)")

    # --- grafici + report HTML ---
    b64_freq = grafico_frequenze(df)
    b64_classi = grafico_classi(classi)
    out = "report.html"
    genera_report(nome, seq, r, df, b64_freq, b64_classi, out)
    print(f"\nReport salvato in '{out}' (aprilo nel browser; da lì 'Stampa -> PDF').")


if __name__ == "__main__":
    main()
