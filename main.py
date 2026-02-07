"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.02
Système : Phénotypage, Scanner IA, Nutrition DZ, Santé & Laboratoire Génomique
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import os
from datetime import datetime, date, timedelta
from Bio import Align  
from Bio.Seq import Seq

# ============================================================================
# BLOC 1 : GESTION DE LA BASE DE DONNÉES
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_pro.db"):
        self.db_path = db_path
        if not os.path.exists('data'): os.makedirs('data')
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row

    def execute_query(self, query: str, params: tuple = ()):
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except sqlite3.Error as e:
            st.error(f"Erreur SQL: {e}")
            return None

    def fetch_all_as_df(self, query: str, params: tuple = ()):
        return pd.read_sql_query(query, self.conn, params=params)

def init_database(db: DatabaseManager):
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE,
            nom TEXT, race TEXT, poids REAL, hauteur REAL, longueur REAL, 
            largeur_bassin REAL, circ_canon REAL, note_mamelle INTEGER, 
            profondeur_mamelle REAL, attache_arriere REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS genotypes (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            gene_nom TEXT, score_homologie REAL, classement TEXT, date_test DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, qte_matin REAL, qte_soir REAL, taux_gras REAL
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, maladie_suspectee TEXT, traitement TEXT, rappel_date DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# BLOC 2 : MOTEUR BIOINFORMATIQUE EXPERT (VOTRE LABO)
# ============================================================================

class BioInfoEngine:
    GENES_QUALITE = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Muscle/Viande)": "AAGCTTGATTAGCAGGTTCCCGG",
        "CAST (Tendreté)": "TGGGGCCCAAGTCGATTGCAGAA",
        "LALBA (Lait)": "GCTAGCTAGCTAGCTGATCGATG"
    }
    GENES_SANTE = {
        "Scrapie_ARR (Résistance)": "TGGTACCCATAATCAGTGGAACA",
        "Scrapie_VRQ (Sensibilité)": "TGGTAGCCATAATCAGTGGAACA",
        "Arachnomélie": "CCGTAGCTAGCTGATCGATCGTA"
    }

    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'

    def filtrer_sequence(self, seq):
        if ">" in seq: seq = "".join(seq.split('\n')[1:])
        return seq.upper().strip().replace(" ", "").replace("\n", "").replace("\r", "")

    def extraire_multi_fasta(self, raw_text):
        sequences = {}
        current_id = None
        for line in raw_text.split('\n'):
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line.upper().replace(" ", "")
        return sequences if sequences else {"Individu_Unique": raw_text.upper().replace(" ", "")}

    def calculer_heterozygotie(self, sequences_dict):
        if len(sequences_dict) < 2: return 0.0
        seqs = list(sequences_dict.values())
        distances = []
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                score = self.aligner.score(seqs[i], seqs[j])
                distances.append(1 - (score / max(len(seqs[i]), len(seqs[j]))))
        return round(np.mean(distances) * 100, 2)

    def traduire(self, dna_seq):
        try:
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            return str(Seq(clean_dna).translate(to_stop=True))
        except: return "Séquence invalide"

# ============================================================================
# BLOC 3 : INTERFACE GÉNOMIQUE (LABO)
# ============================================================================

def bloc_genomique(db, genomique):
    st.title("🧬 Laboratoire Génomique Intégré")
    dna_input = st.text_area("Collez vos séquences ADN (NCBI / FASTA)", height=150)
    
    if dna_input:
        data_dict = genomique.extraire_multi_fasta(dna_input)
        tab_ciblage, tab_diversite, tab_traduction = st.tabs([
            "🎯 Génotypage Ciblé", "📊 Analyse de Population", "🔬 Séquençage Moléculaire"
        ])

        with tab_ciblage:
            st.subheader("Diagnostic de marqueurs par animal")
            c1, c2 = st.columns([1, 2])
            with c1:
                target_id = st.selectbox("Sélectionner l'individu", list(data_dict.keys()))
                tous_m = {**genomique.GENES_QUALITE, **genomique.GENES_SANTE}
                choix_m = st.multiselect("Marqueurs à cibler", list(tous_m.keys()), default=["FecB (Prolificité)"])
                btn_run = st.button("Lancer le diagnostic")
            with c2:
                if btn_run:
                    seq_test = data_dict[target_id]
                    res_list = []
                    for m_nom in choix_m:
                        ref = tous_m[m_nom]
                        score = round((genomique.aligner.score(seq_test, ref) / len(ref)) * 100, 2)
                        if "ARR" in m_nom: verdict = "💎 R1 (TRÈS RÉSISTANT)" if score > 88 else "CLASSIQUE"
                        elif "VRQ" in m_nom: verdict = "❌ SENSIBLE" if score > 88 else "SAIN"
                        elif score > 85: verdict = "✅ ÉLITE / POSITIF"
                        else: verdict = "➖ NÉGATIF"
                        
                        db.execute_query("INSERT INTO genotypes (brebis_id, gene_nom, score_homologie, classement, date_test) VALUES (?,?,?,?,?)",
                                         (target_id, m_nom, score, verdict, date.today()))
                        res_list.append({"Marqueur": m_nom, "Homologie": f"{score}%", "Verdict": verdict})
                    st.table(pd.DataFrame(res_list))

        with tab_diversite:
            if len(data_dict) > 1:
                score_h = genomique.calculer_heterozygotie(data_dict)
                st.metric("Indice d'Hétérozygotie", f"{score_h}%")
                if score_h < 12: st.error("⚠️ Risque de consanguinité !")
                else: st.success("✅ Bonne variabilité génétique.")

        with tab_traduction:
            st.subheader("Séquence Protéique")
            id_tr = st.selectbox("Traduire :", list(data_dict.keys()), key="tr_sel")
            st.code(genomique.traduire(data_dict[id_tr]))

# ============================================================================
# BLOC 4 : NUTRITION DZ FLEXIBLE
# ============================================================================

def bloc_nutrition():
    st.title("🌾 Nutrition & Rations DZ")
    with st.expander("💰 Prix du Marché (DA/Quintal)", expanded=True):
        c1, c2, c3 = st.columns(3)
        p_orge = c1.number_input("Orge (Chaïr)", value=5500)
        p_son = c2.number_input("Son (Nkhala)", value=2500)
        p_foin = c3.number_input("Foin/Luzerne", value=4500)

    stade = st.selectbox("Stade physiologique", ["Entretien", "Gestation", "Allaitement", "Engraissement"])
    rations = {"Entretien": [0.4, 0.2, 1.0], "Gestation": [0.7, 0.3, 1.2], "Allaitement": [1.0, 0.4, 1.5], "Engraissement": [0.8, 0.2, 0.5]}
    r = rations[stade]
    cout = (p_orge/100*r[0]) + (p_son/100*r[1]) + (p_foin/100*r[2])
    st.metric("Coût journalier par tête", f"{round(cout, 2)} DA")

# ============================================================================
# BLOC 5 : PHÉNOTYPAGE & SCANNER IA
# ============================================================================

def bloc_inscription(db):
    st.title("📝 Phénotypage & Scanner IA")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Scanner")
        etalon = st.selectbox("Étalon", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
        img = st.file_uploader("Importer l'image de la brebis", type=['jpg','png'])
        if img: st.image(img, width=300)
    with col2:
        st.subheader("Fiche")
        with st.form("f_ins"):
            uid = st.text_input("ID Boucle")
            h = st.number_input("Hauteur (cm)")
            l = st.number_input("Longueur (cm)")
            b = st.number_input("Bassin (cm)")
            c = st.number_input("Canon (cm)")
            prof = st.number_input("Profondeur Mamelle (cm)")
            if st.form_submit_button("Enregistrer"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, hauteur, longueur, largeur_bassin, circ_canon, profondeur_mamelle, created_at) VALUES (?,?,?,?,?,?,?)", 
                                 (uid, h, l, b, c, prof, date.today()))
                st.success("Données enregistrées !")

# ============================================================================
# BLOC 6 : SANTÉ & CALENDRIER IA
# ============================================================================

def bloc_sante():
    st.title("🩺 Santé & Calendrier IA")
    symp = st.multiselect("Symptômes détectés :", ["Toux", "Boiterie", "Lésions buccales"])
    if st.button("Analyser"):
        if "Boiterie" in symp: st.error("Suspicion de Piétin. Traitement : Sulfate de Zinc.")
    
    st.subheader("📅 Rappels Vaccinaux")
    cal = [{"Mois": "Mars", "Vaccin": "Pestevax (Entéro)"}, {"Mois": "Avril", "Vaccin": "Clavelée"}]
    st.table(pd.DataFrame(cal))

# ============================================================================
# MAIN : NAVIGATION
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🧬")
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    db, genomique = st.session_state.db, st.session_state.genomique
    st.sidebar.title("🐑 EXPERT OVIN DZ")
    menu = ["📊 Dashboard", "📝 Phénotype/Scanner", "🧬 Laboratoire Génomique", "🌾 Nutrition DZ", "🩺 Santé IA"]
    choice = st.sidebar.radio("Navigation", menu)

    if choice == "📊 Dashboard": st.dataframe(db.fetch_all_as_df("SELECT * FROM brebis"))
    elif choice == "📝 Phénotype/Scanner": bloc_inscription(db)
    elif choice == "🧬 Laboratoire Génomique": bloc_genomique(db, genomique)
    elif choice == "🌾 Nutrition DZ": bloc_nutrition()
    elif choice == "🩺 Santé IA": bloc_sante()

if __name__ == "__main__":
    main()
