"""
PROJET : EXPERT OVIN DZ PRO (VERSION INTÉGRALE 2026)
Domaine : Sélection génétique, Génomique, Morphométrie & Gestion Laitière
Auteur : rahim LABORATOIRE GenApAgiE 
"""

"""
EXPERT OVIN DZ PRO - VERSION ULTRA EXPERT (SQLITE INTEGRATED)
Fusion : Database, Génomique SNP, Biochimie & Morphométrie
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date
from typing import Dict, List, Any

# ============================================================================
# 1. GESTIONNAIRE DE BASE DE DONNÉES (SQLITE)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_manager.db"):
        self.db_path = db_path
        # Création du dossier data si inexistant
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
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, date_naissance DATE, race TEXT, sexe TEXT,
            statut TEXT DEFAULT 'active', created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )""",
        """CREATE TABLE IF NOT EXISTS caracteres_morpho (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, date_mesure DATE, 
            hauteur REAL, longueur REAL, ial REAL, yield_est REAL,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS sequences_genetiques (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, accession_number TEXT UNIQUE,
            sequence_type TEXT, longueur INTEGER, date_sequencage DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )"""
    ]
    for table_sql in tables:
        db_manager.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE MÉTIER (SNP & BIOCHIMIE)
# ============================================================================

class IntegrationGenomique:
    def analyser_snp(self, sequence_reference: str, sequence_etudiee: str) -> Dict:
        if len(sequence_reference) != len(sequence_etudiee):
            return {"erreur": "Séquences de longueurs différentes"}
        snps = []
        for i, (ref, etu) in enumerate(zip(sequence_reference, sequence_etudiee)):
            if ref != etu:
                snps.append({'pos': i + 1, 'ref': ref, 'alt': etu})
        return {'total': len(snps), 'details': snps}

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin SQL Pro", layout="wide")

    # Initialisation DB
    if 'db_manager' not in st.session_state:
        db_m = DatabaseManager()
        init_database(db_m)
        st.session_state.db_manager = db_m

    db = st.session_state.db_manager

    # Sidebar Navigation
    st.sidebar.title("🧬 Ovin Manager SQL")
    menu = ["📊 Dashboard", "📝 Enregistrement", "📷 Morphométrie", "🧪 Biochimie", "🧬 Génomique SNP"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- ENREGISTREMENT ---
    if choice == "📝 Enregistrement":
        st.title("🐑 Inscription d'un nouvel animal")
        with st.form("new_animal"):
            identifiant = st.text_input("Identifiant Unique (Boucle)*")
            nom = st.text_input("Nom/Surnom")
            race = st.selectbox("Race", ["Ouled Djellal", "Rembi", "Lacaune", "Hamra"])
            date_n = st.date_input("Date de naissance", date(2022, 1, 1))
            if st.form_submit_button("Ajouter à la base SQL"):
                query = "INSERT INTO brebis (identifiant_unique, nom, race, date_naissance) VALUES (?, ?, ?, ?)"
                db.execute_query(query, (identifiant, nom, race, date_n))
                st.success(f"Animal {identifiant} ajouté à la base de données !")

    # --- DASHBOARD ---
    elif choice == "📊 Dashboard":
        st.title("📊 Vue d'ensemble du Troupeau")
        df_brebis = db.fetch_all_as_df("SELECT * FROM brebis")
        st.dataframe(df_brebis, use_container_width=True)

    # --- MORPHOMÉTRIE (SCANNER) ---
    elif choice == "📷 Morphométrie":
        st.title("📏 Mesures Morphométriques")
        df_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        if df_list.empty:
            st.warning("Aucun animal en base.")
        else:
            with st.form("morpho_entry"):
                target = st.selectbox("Sélectionner la brebis", df_list['identifiant_unique'])
                h = st.number_input("Hauteur (cm)", 50.0, 100.0, 70.0)
                l = st.number_input("Longueur (cm)", 50.0, 120.0, 80.0)
                # Calcul automatique
                poids = (h**2 * l) / 10800
                if st.form_submit_button("Sauvegarder les mesures"):
                    query = "INSERT INTO caracteres_morpho (brebis_id, date_mesure, hauteur, longueur, yield_est) VALUES (?, ?, ?, ?, ?)"
                    db.execute_query(query, (target, date.today(), h, l, poids))
                    st.success("Mesures enregistrées dans SQL.")

    # --- GÉNOMIQUE SNP ---
    elif choice == "🧬 Génomique SNP":
        st.title("🧬 Laboratoire Génomique")
        engine = IntegrationGenomique()
        
        df_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        target_dna = st.selectbox("Animal pour analyse", df_list['identifiant_unique'])
        
        seq_ref = st.text_area("Séquence Référence NCBI", "ATGCGTACGTAGCTAGCTAGCGATCGATCGATCGA")
        seq_stu = st.text_area("Séquence Séquencée", "ATGCGTACGTGGCTAGCTAGCCATCGATCGATCGA")
        
        if st.button("Lancer Analyse SNP"):
            res = engine.analyser_snp(seq_ref, seq_stu)
            st.json(res)
            # Optionnel : Sauvegarde auto dans SQL
            query = "INSERT INTO sequences_genetiques (brebis_id, sequence_type, longueur, date_sequencage) VALUES (?, ?, ?, ?)"
            db.execute_query(query, (target_dna, "ADN", len(seq_stu), date.today()))
            st.info("Données de séquençage archivées.")

if __name__ == "__main__":
    main()
