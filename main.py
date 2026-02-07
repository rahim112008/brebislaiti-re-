"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.02
Système Modulaire : Phénotypage, Lait, Santé & Laboratoire Génomique
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
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
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, poids REAL, note_mamelle INTEGER, 
            tour_poitrine REAL, longueur REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS genotypes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, gene_nom TEXT, resultat TEXT, 
            score_homologie REAL, classement TEXT, date_test DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, quantite_lait REAL
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# BLOC 2 : MOTEUR BIOINFORMATIQUE EXPERT
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
# BLOC 3 : FONCTIONS D'INTERFACE (PAR MODULE)
# ============================================================================

def bloc_dashboard(db):
    st.title("📊 Performances du Cheptel")
    t1, t2, t3 = st.tabs(["📋 Effectifs", "🧬 Résultats Génomiques", "📈 Courbes Lait"])
    with t1:
        st.dataframe(db.fetch_all_as_df("SELECT * FROM brebis"), width='stretch')
    with t2:
        st.dataframe(db.fetch_all_as_df("SELECT * FROM genotypes ORDER BY date_test DESC"), width='stretch')

def bloc_inscription(db):
    st.title("📝 Enregistrement Phénotypique")
    with st.form("form_ins"):
        c1, c2 = st.columns(2)
        uid = c1.text_input("ID Boucle")
        nom = c1.text_input("Nom")
        race = c2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune", "Autre"])
        poids = c2.number_input("Poids (kg)", 10.0, 150.0, 50.0)
        note_m = st.slider("Note Mamelle", 1, 10, 5)
        if st.form_submit_button("Enregistrer"):
            db.execute_query("INSERT INTO brebis (identifiant_unique, nom, race, poids, note_mamelle, created_at) VALUES (?,?,?,?,?,?)",
                             (uid, nom, race, poids, note_m, date.today()))
            st.success("Animal enregistré !")

def bloc_genomique(db, genomique):
    st.title("🧬 Laboratoire Génomique Intégré")
    dna_input = st.text_area("Collez vos séquences ADN (NCBI / FASTA)", height=150)
    if dna_input:
        data_dict = genomique.extraire_multi_fasta(dna_input)
        tab1, tab2, tab3 = st.tabs(["🎯 Génotypage Ciblé", "📊 Diversité", "🔬 Traduction"])
        
        with tab1:
            c1, c2 = st.columns([1, 2])
            with c1:
                target_id = st.selectbox("Individu", list(data_dict.keys()))
                tous_m = {**genomique.GENES_QUALITE, **genomique.GENES_SANTE}
                choix_m = st.multiselect("Marqueurs", list(tous_m.keys()), default=["FecB (Prolificité)"])
                btn = st.button("Lancer le diagnostic")
            with c2:
                if btn:
                    seq_test = data_dict[target_id]
                    res_list = []
                    for m_nom in choix_m:
                        ref = tous_m[m_nom]
                        score = round((genomique.aligner.score(seq_test, ref) / len(ref)) * 100, 2)
                        verdict = "💎 ÉLITE" if score > 85 else "➖ NÉGATIF"
                        db.execute_query("INSERT INTO genotypes (brebis_id, gene_nom, score_homologie, classement, date_test) VALUES (?,?,?,?,?)",
                                         (target_id, m_nom, score, verdict, date.today()))
                        res_list.append({"Marqueur": m_nom, "Homologie": f"{score}%", "Verdict": verdict})
                    st.table(pd.DataFrame(res_list))
        
        with tab2:
            if len(data_dict) > 1:
                score_h = genomique.calculer_heterozygotie(data_dict)
                st.metric("Hétérozygotie", f"{score_h}%")
            else: st.info("Besoin de 2 séquences minimum.")
            
        with tab3:
            id_tr = st.selectbox("Traduire :", list(data_dict.keys()))
            st.code(genomique.traduire(data_dict[id_tr]))

# ============================================================================
# BLOC MAIN : POINT D'ENTRÉE
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    db, genomique = st.session_state.db, st.session_state.genomique

    st.sidebar.title("🐑 EXPERT OVIN DZ")
    menu = ["📊 Dashboard Élite", "📝 Inscription", "📷 Scanner IA 1m", "🧬 Laboratoire Génomique", "🥛 Lait", "🩺 Santé"]
    choice = st.sidebar.radio("Navigation", menu)

    if choice == "📊 Dashboard Élite": bloc_dashboard(db)
    elif choice == "📝 Inscription": bloc_inscription(db)
    elif choice == "📷 Scanner IA 1m": st.camera_input("Scanner l'animal")
    elif choice == "🧬 Laboratoire Génomique": bloc_genomique(db, genomique)
    elif choice == "🥛 Lait": st.title("🥛 Suivi Laitier")
    elif choice == "🩺 Santé": st.title("🩺 Carnet de Santé")

if __name__ == "__main__":
    main()
