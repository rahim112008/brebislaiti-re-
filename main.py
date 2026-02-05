"""
EXPERT OVIN DZ PRO - VERSION ULTRA EXPERT (SQLITE INTEGRATED)
Fusion : Database, Génomique SNP, Biochimie & Morphométrie
Auteur : rahim LABORATOIRE GenApAgiE 
"""

"""
EXPERT OVIN DZ PRO - VERSION INTÉGRALE CUMULATIVE 2026
Système Expert : SQL / Génomique / Biochimie / Morphométrie
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
from datetime import datetime, date
from typing import Dict, List, Any

# ============================================================================
# 1. GESTIONNAIRE DE BASE DE DONNÉES (SQLITE) - MODULE DATA
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_manager.db"):
        self.db_path = db_path
        if not os.path.exists('data'):
            os.makedirs('data')
        self.conn = None
        self.connect()
    
    def connect(self):
        try:
            self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
            self.conn.row_factory = sqlite3.Row
            return True
        except sqlite3.Error as e:
            st.error(f"Erreur connexion DB: {e}")
            return False
    
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

def init_database(db_manager: DatabaseManager):
    """Initialisation complète de toutes vos tables SQL"""
    tables = [
        # Table brebis
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, date_naissance DATE, race TEXT, sexe TEXT,
            statut TEXT DEFAULT 'active', notes TEXT, created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )""",
        # Table gestations
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT,
            date_insemination DATE, date_mise_bas_prevu DATE, statut TEXT DEFAULT 'en_cours',
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table caracteres_morpho
        """CREATE TABLE IF NOT EXISTS caracteres_morpho (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_mesure DATE,
            hauteur REAL, longueur REAL, ial REAL, yield_est REAL, 
            attache_arriere INTEGER, sillon_median INTEGER,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table biochimie
        """CREATE TABLE IF NOT EXISTS biochimie_lait (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_analyse DATE,
            fat REAL, prot REAL, bhb REAL, ratio_tbtp REAL, diagnostic TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table sequences_genetiques
        """CREATE TABLE IF NOT EXISTS sequences_genetiques (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT,
            accession_number TEXT UNIQUE, sequence_type TEXT, longueur INTEGER,
            date_sequencage DATE, labo TEXT, sequence_dna TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )"""
    ]
    for table_sql in tables:
        db_manager.execute_query(table_sql)

# ============================================================================
# 2. MODULE GÉNOMIQUE & SNP - MODULE BIOINFO
# ============================================================================

class IntegrationGenomique:
    def __init__(self, email: str = "labo@expert-ovin.dz"):
        self.email = email
    
    def analyser_snp(self, sequence_ref: str, sequence_stu: str) -> Dict:
        if len(sequence_ref) != len(sequence_stu):
            return {"erreur": "Les séquences doivent avoir la même longueur"}
        snps = []
        for i, (ref, etu) in enumerate(zip(sequence_ref, sequence_stu)):
            if ref != etu:
                snps.append({
                    'position': i + 1, 'reference': ref, 'etudie': etu,
                    'type': self._determiner_type_mutation(ref, etu)
                })
        return {'total_snps': len(snps), 'frequence': len(snps)/len(sequence_ref), 'snps': snps}

    def _determiner_type_mutation(self, ref: str, etu: str) -> str:
        transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
        if (ref, etu) in transitions: return 'transition'
        return 'transversion'

    def rechercher_genes_candidats(self, race: str) -> List[Dict]:
        genes = {
            'Lacaune': [{'gene': 'CSN1S1', 'fonction': 'Caséine alpha-S1', 'chromosome': '6'}],
            'Ouled Djellal': [{'gene': 'GDF9', 'fonction': 'Fécondité', 'chromosome': 'X'}]
        }
        return genes.get(race, [{'gene': 'GENERIC', 'fonction': 'Standard', 'chromosome': 'NA'}])

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ PRO", layout="wide")

    # Initialisation DB
    if 'db_manager' not in st.session_state:
        db_m = DatabaseManager()
        init_database(db_m)
        st.session_state.db_manager = db_m

    db = st.session_state.db_manager
    genomique = IntegrationGenomique()

    # --- AUTHENTIFICATION ---
    if 'auth' not in st.session_state: st.session_state.auth = False
    if not st.session_state.auth:
        st.title("🔐 Accès Expert Ovin DZ Pro")
        pwd = st.text_input("Mot de passe", type="password")
        if st.button("Connexion"):
            if pwd == "admin123":
                st.session_state.auth = True
                st.rerun()
        return

    # --- SIDEBAR NAVIGATION ---
    st.sidebar.title("🐑 Menu Intégral")
    menu = ["📊 Dashboard", "📝 Inscription Animal", "📷 Morphométrie IA", "🧪 Biochimie Laitière", "🧬 Génomique & SNP"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE 1: DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title("📊 Tableau de Bord Central")
        df_brebis = db.fetch_all_as_df("SELECT * FROM brebis")
        if df_brebis.empty:
            st.info("Aucune donnée enregistrée.")
        else:
            col1, col2 = st.columns(2)
            col1.metric("Total Brebis", len(df_brebis))
            st.dataframe(df_brebis, use_container_width=True)

    # --- MODULE 2: INSCRIPTION ---
    elif choice == "📝 Inscription Animal":
        st.title("📝 Enregistrement Permanent (SQL)")
        with st.form("new_brebis"):
            identifiant = st.text_input("ID Unique (Boucle)*")
            race = st.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra"])
            date_n = st.date_input("Date de Naissance")
            notes = st.text_area("Observations")
            if st.form_submit_button("Sauvegarder en Base"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, race, date_naissance, notes) VALUES (?,?,?,?)",
                                (identifiant, race, date_n, notes))
                st.success(f"Animal {identifiant} enregistré définitivement.")

    # --- MODULE 3: MORPHOMETRIE IA ---
    elif choice == "📷 Morphométrie IA":
        st.title("📐 Morphométrie & Scanner IA")
        
        df_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        if not df_list.empty:
            with st.form("scanner_form"):
                target = st.selectbox("Brebis", df_list['identifiant_unique'])
                h = st.number_input("Hauteur au garrot (cm)", 50.0, 100.0, 70.0)
                l = st.number_input("Longueur de corps (cm)", 50.0, 120.0, 80.0)
                attache = st.slider("Score Attache Arrière (1-9)", 1, 9, 5)
                sillon = st.slider("Score Sillon Médian (1-9)", 1, 9, 5)
                if st.form_submit_button("💾 Enregistrer Mesures"):
                    ial = (attache * 0.6) + (sillon * 0.4)
                    db.execute_query("INSERT INTO caracteres_morpho (brebis_id, date_mesure, hauteur, longueur, ial) VALUES (?,?,?,?,?)",
                                    (target, date.today(), h, l, ial))
                    st.success("Mesures et pointage mamelle archivés.")
                    

    # --- MODULE 4: BIOCHIMIE ---
    elif choice == "🧪 Biochimie Laitière":
        st.title("🧪 Analyse Biochimique & Métabolique")
        df_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        if not df_list.empty:
            with st.form("biochem_form"):
                target = st.selectbox("Brebis", df_list['identifiant_unique'])
                fat = st.number_input("TB (g/L)", 20.0, 95.0, 45.0)
                prot = st.number_input("TP (g/L)", 20.0, 75.0, 38.0)
                bhb = st.number_input("BHB (mmol/L)", 0.0, 5.0, 0.5)
                if st.form_submit_button("Analyser & Sauvegarder"):
                    ratio = fat / prot
                    diag = "Normal" if bhb < 1.2 else "Cétose Subclinique ⚠️"
                    db.execute_query("INSERT INTO biochimie_lait (brebis_id, date_analyse, fat, prot, bhb, ratio_tbtp, diagnostic) VALUES (?,?,?,?,?,?,?)",
                                    (target, date.today(), fat, prot, bhb, ratio, diag))
                    st.success(f"Analyse terminée. Diagnostic : {diag}")
                    

    # --- MODULE 5: GÉNOMIQUE & SNP ---
    elif choice == "🧬 Génomique & SNP":
        st.title("🧬 Bioinformatique & GBLUP")
        df_list = db.fetch_all_as_df("SELECT identifiant_unique, race FROM brebis")
        if not df_list.empty:
            target = st.selectbox("Sélectionner la brebis", df_list['identifiant_unique'])
            race_sel = df_list[df_list['identifiant_unique'] == target]['race'].values[0]
            
            tab1, tab2 = st.tabs(["Analyse SNP", "Gènes Candidats"])
            with tab1:
                seq_ref = st.text_area("Séquence Référence NCBI", "ATGCGTACGTAGCTAGCTAGCGATCGATCGATCGA")
                seq_stu = st.text_area("Séquence Animal", "ATGCGTACGTGGCTAGCTAGCCATCGATCGATCGA")
                if st.button("Lancer l'analyse SNP"):
                    res = genomique.analyser_snp(seq_ref, seq_stu)
                    st.json(res)
                    
            with tab2:
                st.write(f"Gènes prioritaires pour la race {race_sel} :")
                st.table(genomique.rechercher_genes_candidats(race_sel))

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.auth = False
        st.rerun()

if __name__ == "__main__":
    main()
    main()
