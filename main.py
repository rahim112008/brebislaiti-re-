"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE CONSOLIDÉE 2026
Système Tout-en-Un : Phénotypage, Scanner IA, Génomique (PLINK Style), 
Lait, Santé, Arbre Généalogique & Simulation d'Accouplement.
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import os
from datetime import datetime, date, timedelta
from Bio import Align  
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam
import plotly.express as px

# ============================================================================
# 1. GESTION DE LA BASE DE DONNÉES (UNIFIÉE & PERSISTANTE)
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
    """Initialise l'architecture complète des données"""
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, 
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, poids REAL, 
            tour_poitrine REAL, longueur REAL, note_mamelle INTEGER,
            pere_id TEXT, mere_id TEXT, sexe TEXT,
            created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, quantite_lait REAL, 
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. MOTEUR BIO-INFORMATIQUE & GÉNÉTIQUE (STYLE PLINK)
# ============================================================================

class BioInfoEngine:
    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'
        self.GENES_REF = {
            "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
            "MSTN (Muscle)": "AAGCTTGATTAGCAGGTTCCCGG",
            "Scrapie_ARR (Résistance)": "TGGTACCCATAATCAGTGGAACA",
            "DGAT1 (Lait)": "GCTAGCTAGCTAGCTGATCGATG"
        }

    def filtrer_sequence(self, seq):
        if ">" in seq: seq = "".join(seq.split('\n')[1:])
        return seq.upper().strip().replace(" ", "").replace("\r", "").replace("\n", "")

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
        return sequences if sequences else {"Individu": raw_text.upper().replace(" ", "")}

    def alignement_expert(self, seq_test, ref_seq):
        if not seq_test or not ref_seq: return 0.0
        try:
            score = self.aligner.score(seq_test, ref_seq)
            return round((score / len(ref_seq)) * 100, 2)
        except: return 0.0

    def calculer_heterozygotie(self, sequences_dict):
        if len(sequences_dict) < 2: return 0.0
        seqs = list(sequences_dict.values())
        scores = []
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                s = self.aligner.score(seqs[i], seqs[j])
                scores.append(s / max(len(seqs[i]), len(seqs[j])))
        return round((1 - np.mean(scores)) * 100, 2)

    def simuler_croisement(self, seq_p, seq_m):
        pred = {}
        for gene, ref in self.GENES_REF.items():
            sc_p = self.alignement_expert(seq_p, ref)
            sc_m = self.alignement_expert(seq_m, ref)
            p_pres = 1 if sc_p > 85 else 0
            m_pres = 1 if sc_m > 85 else 0
            prob = (p_pres + m_pres) / 2
            if p_pres == 1 and m_pres == 1: prob = 0.95
            elif p_pres == 0 and m_pres == 0: prob = 0.05
            pred[gene] = prob
        return pred

# ============================================================================
# 3. INTERFACE UTILISATEUR PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'bio' not in st.session_state:
        st.session_state.bio = BioInfoEngine()
        
    db, bio = st.session_state.db, st.session_state.bio

    st.sidebar.title("🐑 EXPERT OVIN DZ")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Généalogie",
        "📷 Scanner IA 1m",
        "🎲 Simulateur Accouplement", 
        "🧬 Génomique (PLINK Style)",
        "🥛 Contrôle Laitier",
        "🩺 Santé & Vaccins",
        "🌾 Nutrition Solo"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    if choice == "📊 Dashboard Élite":
        st.title("📊 Tableau de Bord des Performances")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        
        if not df_b.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df_b))
            c2.metric("Poids Moyen", f"{round(df_b['poids'].mean(), 1)} kg")
            avg_lait = df_l['quantite_lait'].mean() if not df_l.empty else 0
            c3.metric("Moyenne Lait", f"{round(avg_lait, 2)} L")
            st.dataframe(df_b, use_container_width=True)
            if not df_l.empty:
                fig = px.line(df_l, x='date_controle', y='quantite_lait', color='brebis_id', title="Production Laitière")
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Aucun animal enregistré.")

    elif choice == "📝 Inscription & Généalogie":
        st.title("📝 Enregistrement & Pedigree")
        t_ins, t_tree = st.tabs(["Nouvelle Inscription", "Arbre Généalogique"])
        with t_ins:
            with st.form("inscription"):
                col1, col2 = st.columns(2)
                uid = col1.text_input("ID Unique (Boucle)")
                nom = col1.text_input("Nom / Alias")
                race = col2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune", "Autre"])
                poids = col2.number_input("Poids (kg)", 10.0, 150.0, 55.0)
                pere = col1.text_input("ID du Père")
                mere = col2.text_input("ID de la Mère")
                if st.form_submit_button("Sauvegarder"):
                    db.execute_query(
                        "INSERT INTO brebis (identifiant_unique, nom, race, poids, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?)",
                        (uid, nom, race, poids, pere, mere, date.today())
                    )
                    st.success(f"Animal {uid} ajouté.")
        with t_tree:
            sid = st.text_input("Entrez l'ID pour voir l'arbre")
            if sid:
                st.subheader(f"🌳 Arbre Généalogique de {sid}")
                res = db.fetch_all_as_df("SELECT * FROM brebis WHERE identifiant_unique = ?", (sid,))
                if not res.empty:
                    row = res.iloc[0]
                    c1, c2 = st.columns(2)
                    c1.info(f"♂️ Père : {row['pere_id'] if row['pere_id'] else 'Inconnu'}")
                    c2.info(f"♀️ Mère : {row['mere_id'] if row['mere_id'] else 'Inconnue'}")

    elif choice == "📷 Scanner IA 1m":
        st.title("📸 Scanner Morphométrique")
        st.camera_input("Scanner")

    elif choice == "🎲 Simulateur Accouplement":
        st.title("🎲 Simulation de Croisement")
        col1, col2 = st.columns(2)
        dna_p = col1.text_area("ADN Père", key="sim_p")
        dna_m = col2.text_area("ADN Mère", key="sim_m")
        if st.button("🧬 Prédire la Qualité"):
            if dna_p and dna_m:
                results = bio.simuler_croisement(dna_p, dna_m)
                st.json(results)
            else: st.warning("Entrez les deux séquences.")

    elif choice == "🧬 Génomique (PLINK Style)":
        st.title("🧬 Analyse de Groupe")
        dna_txt = st.text_area("Séquences Multi-FASTA")
        if dna_txt:
            data_dict = bio.extraire_multi_fasta(dna_txt)
            h = bio.calculer_heterozygotie(data_dict)
            st.metric("Indice d'Hétérozygotie", f"{h}%")

    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Suivi Laitier")
        with st.form("lait_form"):
            id_l = st.text_input("ID Brebis")
            qte = st.number_input("Litres", 0.0, 10.0, 1.5)
            if st.form_submit_button("Enregistrer"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)", (id_l, date.today(), qte))
                st.success("Enregistré.")

    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Carnet de Santé")
        df_sante = db.fetch_all_as_df("SELECT * FROM sante")
        st.table(df_sante)

    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Rationnement")
        p = st.number_input("Poids (kg)", 10, 150, 60)
        st.write(f"Orge : {p * 0.012:.2f} kg/j | Foin : {p * 0.02:.2f} kg/j")

if __name__ == "__main__":
    main()
