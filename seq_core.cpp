// =====================================================================
//  seq_core.cpp  -  Core C++ per l'analisi di una sequenza proteica.
// ---------------------------------------------------------------------
//  Riceve una sequenza di amminoacidi (stringa) e calcola, alla velocità
//  dell'hardware:
//    * conteggio dei 20 amminoacidi standard
//    * lunghezza e numero di legami peptidici (L-1)
//    * numero di cisteine (-> potenziali ponti disolfuro)
//    * peso molecolare (somma delle masse dei residui + acqua)
//    * carica netta a un dato pH
//    * punto isoelettrico (pI) risolto NUMERICAMENTE per bisezione
//
//  Il pI è il motivo per cui ha senso usare il C++: non è una semplice
//  conta, ma la ricerca dello zero di una funzione (carica netta = 0)
//  fatta iterando decine di volte. È il tipo di calcolo "pesante" che
//  conviene tenere fuori da Python.
//
//  Esposto a Python con pybind11 come modulo "seq_core".
// =====================================================================

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>      // per convertire std::map / std::string in dict / str Python
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>

namespace py = pybind11;

// --- Dati biochimici dei 20 amminoacidi standard --------------------
// Massa media del RESIDUO (amminoacido - acqua), in Dalton.
static double massa_residuo(char a) {
    switch (a) {
        case 'A': return 71.0788;  case 'R': return 156.1875;
        case 'N': return 114.1038; case 'D': return 115.0886;
        case 'C': return 103.1388; case 'E': return 129.1155;
        case 'Q': return 128.1307; case 'G': return 57.0519;
        case 'H': return 137.1411; case 'I': return 113.1594;
        case 'L': return 113.1594; case 'K': return 128.1741;
        case 'M': return 131.1926; case 'F': return 147.1766;
        case 'P': return 97.1167;  case 'S': return 87.0782;
        case 'T': return 101.1051; case 'W': return 186.2132;
        case 'Y': return 163.1760; case 'V': return 99.1326;
        default:  return 0.0;
    }
}

// Dichiarazione anticipata: la scala Kyte-Doolittle è definita più sotto,
// ma serve già qui per calcolare il GRAVY dentro analizza().
static double kyte_doolittle(char a);

// Parametri conformazionali di Chou-Fasman (×100): propensione di ciascun
// amminoacido a elica (Pa), foglietto beta (Pb) e turn/coil (Pt).
// Classifichiamo ogni residuo nella conformazione con propensione massima:
// è una stima INDICATIVA della composizione, non una vera predizione strutturale.
static void chou_fasman(char a, double& Pa, double& Pb, double& Pt) {
    switch (a) {
        case 'A': Pa=142; Pb=83;  Pt=66;  break;  case 'R': Pa=98;  Pb=93;  Pt=95;  break;
        case 'N': Pa=67;  Pb=89;  Pt=156; break;  case 'D': Pa=101; Pb=54;  Pt=146; break;
        case 'C': Pa=70;  Pb=119; Pt=119; break;  case 'Q': Pa=111; Pb=110; Pt=98;  break;
        case 'E': Pa=151; Pb=37;  Pt=74;  break;  case 'G': Pa=57;  Pb=75;  Pt=156; break;
        case 'H': Pa=100; Pb=87;  Pt=95;  break;  case 'I': Pa=108; Pb=160; Pt=47;  break;
        case 'L': Pa=121; Pb=130; Pt=59;  break;  case 'K': Pa=114; Pb=74;  Pt=101; break;
        case 'M': Pa=145; Pb=105; Pt=60;  break;  case 'F': Pa=113; Pb=138; Pt=60;  break;
        case 'P': Pa=57;  Pb=55;  Pt=152; break;  case 'S': Pa=77;  Pb=75;  Pt=143; break;
        case 'T': Pa=83;  Pb=119; Pt=96;  break;  case 'W': Pa=108; Pb=137; Pt=96;  break;
        case 'Y': Pa=69;  Pb=147; Pt=114; break;  case 'V': Pa=106; Pb=170; Pt=50;  break;
        default:  Pa=0;   Pb=0;   Pt=0;   break;
    }
}

// Contributo di carica di una catena laterale a un dato pH.
// Gruppi basici (positivi):  carica = +1 / (1 + 10^(pH - pKa))
// Gruppi acidi (negativi):   carica = -1 / (1 + 10^(pKa - pH))
static double carica_netta(const std::map<char,int>& conteggi,
                           bool ha_residui, double pH) {
    if (!ha_residui) return 0.0;
    auto get = [&](char c){ auto it = conteggi.find(c); return it==conteggi.end()?0:it->second; };

    // pKa (set classico, valori approssimati)
    const double pKa_Nterm = 9.69, pKa_Cterm = 2.34;
    const double pKa_K = 10.5, pKa_R = 12.5, pKa_H = 6.0;       // basici
    const double pKa_D = 3.65, pKa_E = 4.25, pKa_C = 8.3, pKa_Y = 10.07; // acidi

    auto pos = [&](double pKa, int n){ return n / (1.0 + std::pow(10.0, pH - pKa)); };
    auto neg = [&](double pKa, int n){ return n / (1.0 + std::pow(10.0, pKa - pH)); };

    double q = 0.0;
    q += pos(pKa_Nterm, 1);                                  // N-terminale (uno solo)
    q += pos(pKa_K, get('K')) + pos(pKa_R, get('R')) + pos(pKa_H, get('H'));
    q -= neg(pKa_Cterm, 1);                                  // C-terminale (uno solo)
    q -= neg(pKa_D, get('D')) + neg(pKa_E, get('E'))
       + neg(pKa_C, get('C')) + neg(pKa_Y, get('Y'));
    return q;
}

// Risultato dell'analisi, esposto a Python come oggetto con attributi.
struct Risultato {
    int lunghezza = 0;
    int legami_peptidici = 0;
    int cisteine = 0;
    int ponti_disolfuro_max = 0;
    double peso_molecolare = 0.0;
    double carica_a_pH7 = 0.0;
    double punto_isoelettrico = 0.0;
    double gravy = 0.0;               // idrofobicità media (Kyte-Doolittle)
    double indice_alifatico = 0.0;    // indice alifatico (Ikai): stabilità termica
    double estinzione_red = 0.0;      // ε molare a 280 nm, Cys ridotte (M^-1 cm^-1)
    double estinzione_ox = 0.0;       // ε molare a 280 nm, Cys in ponti disolfuro
    double abs280_red = 0.0;          // assorbanza di una soluzione 1 g/L (Cys ridotte)
    double abs280_ox = 0.0;           // assorbanza di una soluzione 1 g/L (Cys ossidate)
    std::vector<double> titolazione_pH;     // curva di titolazione: asse pH
    std::vector<double> titolazione_carica; // curva di titolazione: carica netta
    int ss_helix = 0;                 // residui classificati come elica (Chou-Fasman)
    int ss_sheet = 0;                 // residui classificati come foglietto beta
    int ss_coil = 0;                  // residui classificati come coil/turn
    std::map<char,int> conteggi;      // amminoacido -> numero
    std::map<std::string,int> classi; // classe chimica -> numero
};

Risultato analizza(const std::string& seq_in) {
    Risultato r;
    std::map<char,int> c;

    // Conteggio + peso molecolare (un solo passaggio sulla stringa).
    double massa = 0.0;
    for (char ch : seq_in) {
        char a = std::toupper(static_cast<unsigned char>(ch));
        double m = massa_residuo(a);
        if (m == 0.0) continue;        // ignora caratteri non standard
        c[a]++;
        massa += m;
        r.lunghezza++;
        // struttura secondaria: assegno il residuo alla conformazione dominante
        double Pa, Pb, Pt; chou_fasman(a, Pa, Pb, Pt);
        if (Pa >= Pb && Pa >= Pt)      r.ss_helix++;
        else if (Pb >= Pt)             r.ss_sheet++;
        else                           r.ss_coil++;
    }
    if (r.lunghezza > 0) massa += 18.01528;   // aggiungo una molecola d'acqua

    r.conteggi = c;
    r.legami_peptidici = r.lunghezza > 0 ? r.lunghezza - 1 : 0;
    r.cisteine = c.count('C') ? c['C'] : 0;
    r.ponti_disolfuro_max = r.cisteine / 2;
    r.peso_molecolare = massa;

    // Classi chimico-fisiche (per il grafico a torta).
    auto somma = [&](const std::string& set){
        int t = 0; for (char a : set) if (c.count(a)) t += c[a]; return t;
    };
    r.classi["Idrofobici"]       = somma("AVLIMFWPG");
    r.classi["Polari"]           = somma("STCYNQ");
    r.classi["Carichi positivi"] = somma("KRH");
    r.classi["Carichi negativi"] = somma("DE");

    bool ha = r.lunghezza > 0;
    r.carica_a_pH7 = carica_netta(c, ha, 7.0);

    // Punto isoelettrico: bisezione sul pH per cui la carica netta = 0.
    double lo = 0.0, hi = 14.0;
    for (int it = 0; it < 100 && ha; ++it) {
        double mid = 0.5 * (lo + hi);
        double q = carica_netta(c, ha, mid);
        if (q > 0.0) lo = mid; else hi = mid;   // troppo positivo -> alzo il pH
    }
    r.punto_isoelettrico = ha ? 0.5 * (lo + hi) : 0.0;

    // --- Proprietà biochimiche aggiuntive (stile ProtParam) -----------
    auto cnt = [&](char a){ auto it = c.find(a); return it==c.end()?0:it->second; };

    // GRAVY: media dei valori di idrofobicità Kyte-Doolittle su tutti i residui.
    // Negativo = idrofilico/solubile, positivo = idrofobico.
    if (ha) {
        double s = 0.0;
        for (auto& kv : c) s += kyte_doolittle(kv.first) * kv.second;
        r.gravy = s / r.lunghezza;
    }

    // Indice alifatico (Ikai 1980): volume relativo delle catene alifatiche
    // (Ala, Val, Ile, Leu). Più alto = proteina tendenzialmente più termostabile.
    if (ha) {
        auto molperc = [&](char a){ return 100.0 * cnt(a) / r.lunghezza; };
        r.indice_alifatico = molperc('A') + 2.9 * molperc('V')
                           + 3.9 * (molperc('I') + molperc('L'));
    }

    // Coefficiente di estinzione molare a 280 nm (Pace 1995):
    //   ε = nTyr*1490 + nTrp*5500  (+ 125 per ogni ponte disolfuro)
    // e assorbanza di una soluzione da 1 g/L (la "Abs 0.1%" di ProtParam),
    // utile per ricavare la concentrazione da una lettura allo spettrofotometro.
    r.estinzione_red = cnt('Y') * 1490.0 + cnt('W') * 5500.0;
    r.estinzione_ox  = r.estinzione_red + r.ponti_disolfuro_max * 125.0;
    if (r.peso_molecolare > 0.0) {
        r.abs280_red = r.estinzione_red / r.peso_molecolare;
        r.abs280_ox  = r.estinzione_ox  / r.peso_molecolare;
    }

    // Curva di titolazione: carica netta a 141 valori di pH (0.0, 0.1, ... 14.0),
    // con la stessa funzione usata per il pI (lo zero della curva = punto isoelettrico).
    const int NP = 141;
    for (int i = 0; i < NP; ++i) {
        double pH = 14.0 * i / (NP - 1);
        r.titolazione_pH.push_back(pH);
        r.titolazione_carica.push_back(ha ? carica_netta(c, ha, pH) : 0.0);
    }

    return r;
}

// =====================================================================
//  RILEVATORE DI DOMINI DI MEMBRANA (idrofobicità di Kyte-Doolittle)
// ---------------------------------------------------------------------
//  Calcola il profilo di idrofobicità con la sliding window e individua
//  matematicamente i segmenti transmembrana: tratti in cui la media
//  resta sopra una soglia (default 1.6) per almeno N residui di fila
//  (19-25 aa = lunghezza fisica per attraversare il doppio strato
//  lipidico). Restituisce a Python il profilo, i segmenti trovati e il
//  verdetto topologico.
// =====================================================================

// Scala di idrofobicità di Kyte & Doolittle (1982): valori alti = idrofobico.
static double kyte_doolittle(char a) {
    switch (a) {
        case 'A': return  1.8; case 'R': return -4.5; case 'N': return -3.5;
        case 'D': return -3.5; case 'C': return  2.5; case 'Q': return -3.5;
        case 'E': return -3.5; case 'G': return -0.4; case 'H': return -3.2;
        case 'I': return  4.5; case 'L': return  3.8; case 'K': return -3.9;
        case 'M': return  1.9; case 'F': return  2.8; case 'P': return -1.6;
        case 'S': return -0.8; case 'T': return -0.7; case 'W': return -0.9;
        case 'Y': return -1.3; case 'V': return  4.2;
        default:  return  0.0;
    }
}

struct Segmento {
    int inizio = 0;            // primo residuo (1-based)
    int fine = 0;              // ultimo residuo (1-based)
    int lunghezza = 0;
    double idro_media = 0.0;
    std::string sottosequenza;
};

struct ProfiloIdro {
    std::vector<double> punteggi;   // media della finestra
    std::vector<int> posizioni;     // residuo centrale di ciascuna finestra (1-based)
    std::vector<Segmento> segmenti; // domini transmembrana predetti
    std::string verdetto;           // responso biologico
    int n_domini = 0;
};

ProfiloIdro idrofobicita(const std::string& seq_in,
                         int finestra = 19, double soglia = 1.6,
                         int lung_min = 19) {
    ProfiloIdro out;

    // Tengo solo i residui standard (in maiuscolo). Nessuno dei 20 ha KD = 0,
    // quindi il valore 0.0 identifica i caratteri non standard da scartare.
    std::string seq;
    for (char ch : seq_in) {
        char a = std::toupper(static_cast<unsigned char>(ch));
        if (kyte_doolittle(a) != 0.0) seq.push_back(a);
    }
    const int N = static_cast<int>(seq.size());
    if (N < finestra) { out.verdetto = "Sequenza troppo corta per l'analisi"; return out; }

    // Profilo: media della finestra allineata a sinistra [i, i+finestra-1].
    std::vector<double> smoothed(N - finestra + 1, 0.0);
    for (int i = 0; i <= N - finestra; ++i) {
        double s = 0.0;
        for (int k = 0; k < finestra; ++k) s += kyte_doolittle(seq[i + k]);
        smoothed[i] = s / finestra;
        out.punteggi.push_back(smoothed[i]);
        out.posizioni.push_back(i + finestra / 2 + 1);   // residuo centrale, 1-based
    }

    // Cerco i tratti in cui la media resta >= soglia (potenziali eliche TM).
    // Raccolgo gli intervalli grezzi di residui [r_start, r_end] (0-based).
    std::vector<std::pair<int,int>> grezzi;
    int p = 0;
    const int M = static_cast<int>(smoothed.size());
    while (p < M) {
        if (smoothed[p] >= soglia) {
            int q = p;
            while (q + 1 < M && smoothed[q + 1] >= soglia) ++q;
            grezzi.push_back({p, q + finestra - 1});   // residui coperti dalle finestre
            p = q + 1;
        } else {
            ++p;
        }
    }

    // Fondo gli intervalli che si sovrappongono: due finestre adiacenti che
    // restano alte descrivono la STESSA elica, non due domini distinti.
    std::vector<std::pair<int,int>> uniti;
    for (auto& g : grezzi) {
        if (!uniti.empty() && g.first <= uniti.back().second + 1)
            uniti.back().second = std::max(uniti.back().second, g.second);
        else
            uniti.push_back(g);
    }

    // Trasformo in segmenti, scartando quelli troppo corti per attraversare il bilayer.
    for (auto& u : uniti) {
        int r_start = u.first, r_end = u.second;
        int len = r_end - r_start + 1;
        if (len < lung_min) continue;
        Segmento seg;
        seg.inizio = r_start + 1;              // 1-based per la visualizzazione
        seg.fine = r_end + 1;
        seg.lunghezza = len;
        double s = 0.0;
        for (int k = r_start; k <= r_end; ++k) s += kyte_doolittle(seq[k]);
        seg.idro_media = s / len;
        seg.sottosequenza = seq.substr(r_start, len);
        out.segmenti.push_back(seg);
    }

    // Verdetto topologico in base al numero di domini trovati.
    out.n_domini = static_cast<int>(out.segmenti.size());
    if (out.n_domini == 0)
        out.verdetto = "Proteina globulare solubile (nessun dominio di membrana)";
    else if (out.n_domini == 1)
        out.verdetto = "Proteina di membrana (single-pass)";
    else
        out.verdetto = "Proteina di membrana (multi-pass)";

    return out;
}

PYBIND11_MODULE(seq_core, m) {
    m.doc() = "Core C++ per l'analisi di sequenze proteiche";

    py::class_<Risultato>(m, "Risultato")
        .def_readonly("lunghezza", &Risultato::lunghezza)
        .def_readonly("legami_peptidici", &Risultato::legami_peptidici)
        .def_readonly("cisteine", &Risultato::cisteine)
        .def_readonly("ponti_disolfuro_max", &Risultato::ponti_disolfuro_max)
        .def_readonly("peso_molecolare", &Risultato::peso_molecolare)
        .def_readonly("carica_a_pH7", &Risultato::carica_a_pH7)
        .def_readonly("punto_isoelettrico", &Risultato::punto_isoelettrico)
        .def_readonly("gravy", &Risultato::gravy)
        .def_readonly("indice_alifatico", &Risultato::indice_alifatico)
        .def_readonly("estinzione_red", &Risultato::estinzione_red)
        .def_readonly("estinzione_ox", &Risultato::estinzione_ox)
        .def_readonly("abs280_red", &Risultato::abs280_red)
        .def_readonly("abs280_ox", &Risultato::abs280_ox)
        .def_readonly("titolazione_pH", &Risultato::titolazione_pH)
        .def_readonly("titolazione_carica", &Risultato::titolazione_carica)
        .def_readonly("ss_helix", &Risultato::ss_helix)
        .def_readonly("ss_sheet", &Risultato::ss_sheet)
        .def_readonly("ss_coil", &Risultato::ss_coil)
        .def_readonly("conteggi", &Risultato::conteggi)
        .def_readonly("classi", &Risultato::classi);

    m.def("analizza", &analizza, py::arg("sequenza"),
          "Analizza una sequenza proteica e restituisce un oggetto Risultato");

    py::class_<Segmento>(m, "Segmento")
        .def_readonly("inizio", &Segmento::inizio)
        .def_readonly("fine", &Segmento::fine)
        .def_readonly("lunghezza", &Segmento::lunghezza)
        .def_readonly("idro_media", &Segmento::idro_media)
        .def_readonly("sottosequenza", &Segmento::sottosequenza);

    py::class_<ProfiloIdro>(m, "ProfiloIdro")
        .def_readonly("punteggi", &ProfiloIdro::punteggi)
        .def_readonly("posizioni", &ProfiloIdro::posizioni)
        .def_readonly("segmenti", &ProfiloIdro::segmenti)
        .def_readonly("verdetto", &ProfiloIdro::verdetto)
        .def_readonly("n_domini", &ProfiloIdro::n_domini);

    m.def("idrofobicita", &idrofobicita,
          py::arg("sequenza"), py::arg("finestra") = 19,
          py::arg("soglia") = 1.6, py::arg("lung_min") = 19,
          "Profilo di Kyte-Doolittle e rilevamento dei domini transmembrana");
}
