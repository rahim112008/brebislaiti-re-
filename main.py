"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.02.09
Système : Bio-Informatique, API Ensembl/EBI & Moteur de Prédiction de Valeur (Éco/Santé).
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
import re
import requests
from datetime import datetime, date, timedelta
import random

# ============================================================================
# 1. DATABASE & CONFIG
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
            return None

    def fetch_all_as_df(self, query: str, params: tuple = ()):
        return pd.read_sql_query(query, self.conn, params=params)

def init_database(db: DatabaseManager):
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, sexe TEXT, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS genomique (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, marqueur TEXT, homologie REAL, prediction TEXT, categorie TEXT, date_test DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_controle DATE, quantite_lait REAL
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. MOTEUR IA DE PRÉDICTION ÉCONOMIQUE & SANITAIRE
# ============================================================================

class BioAIEngine:
    # Dictionnaire des gènes d'intérêt (Économique et Sanitaire)
    GENE_CATALOG = {
        "CAST": {"nom": "Calpastatine", "cat": "Économique", "effet": "Tendreté de la viande et croissance."},
        "DGAT1": {"nom": "DGAT1", "cat": "Économique", "effet": "Teneur en gras et rendement laitier."},
        "MSTN": {"nom": "Myostatine", "cat": "Économique", "effet": "Développement musculaire (double muscle)."},
        "PrP": {"nom": "Prion Protein", "cat": "Sanitaire", "effet": "Résistance génétique à la Tremblante (Scrapie)."},
        "FecB": {"nom": "Booroola", "cat": "Sanitaire/Prod", "effet": "Taux de gémellité et santé reproductive."},
        "LEP": {"nom": "Leptine", "cat": "Économique", "effet": "Efficacité alimentaire et réserves de graisses."}
    }

    @staticmethod
    def predire_impact(gene_symbol, score_homologie):
        gene_info = BioAIEngine.GENE_CATALOG.get(gene_symbol, {"nom": "Inconnu", "cat": "Général", "effet": "Inconnu"})
        
        if score_homologie >= 95:
            statut = "Optimal"
            impact = f"Haut potentiel : {gene_info['effet']}"
        elif score_homologie >= 80:
            statut = "Standard"
            impact = f"Performance moyenne pour {gene_info['nom']}."
        else:
            statut = "Risque/Faible"
            impact = f"Performance réduite ou vulnérabilité ({gene_info['cat']})."
            
        return gene_info['cat'], statut, impact

    @staticmethod
    def comparer_sequences(seq_ref, seq_user):
        seq_ref = seq_ref.upper().strip()[:100]
        seq_user = seq_user.upper().strip()[:100]
        if not seq_ref or not seq_user: return 0.0
        matches = sum(1 for a, b in zip(seq_ref, seq_user) if a == b)
        return round((matches / max(len(seq_ref), 1)) * 100, 2)

class WebBioAPI:
    SERVER = "https://rest.ensembl.org"
    @classmethod
    def get_gene_data(cls, symbol):
        ext = f"/lookup/symbol/ovis_aries/{symbol}?"
        r = requests.get(cls.SERVER+ext, headers={"Content-Type": "application/json"})
        if r.ok:
            data = r.json()
            ext_seq = f"/sequence/id/{data['id']}?type=genomic"
            r_seq = requests.get(cls.SERVER+ext_seq, headers={"Content-Type": "text/plain"})
            data['sequence'] = r_seq.text if r_seq.ok else ""
            return data
        return None

# ============================================================================
# 3. INTERFACE STREAMLIT
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ - GENOMIQUE", layout="wide")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    engine = BioAIEngine()
    api_web = WebBioAPI()

    st.sidebar.title("🧬 Ovin-Bio-Predictor")
    menu = ["📊 Tableau de Bord", "🧬 Analyse & Prédiction API", "📉 Études de Valeur"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE : ANALYSE & PRÉDICTION (COEUR DU SYSTÈME) ---
    if choice == "🧬 Analyse & Prédiction API":
        st.title("🧬 Analyse Génomique Comparative")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("Paramètres d'analyse")
            gene_target = st.selectbox("Sélectionner le gène cible", list(engine.GENE_CATALOG.keys()))
            info = engine.GENE_CATALOG[gene_target]
            st.info(f"**Catégorie :** {info['cat']}  \n**Description :** {info['effet']}")
            
            fasta_in = st.text_area("Coller la séquence FASTA de l'animal", height=200)
            animal_id = st.selectbox("Sélectionner l'animal", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")['identifiant_unique'])

        with col2:
            st.subheader("Résultats de l'IA")
            if st.button("Lancer la Prédiction"):
                with st.spinner("Analyse des bases Ensembl..."):
                    remote_data = api_web.get_gene_data(gene_target)
                    if remote_data and fasta_in:
                        # 1. Comparaison
                        score = engine.comparer_sequences(remote_data['sequence'], fasta_in)
                        # 2. Prédiction
                        cat, statut, impact = engine.predire_impact(gene_target, score)
                        
                        # Affichage
                        st.metric("Score d'Homologie", f"{score}%", delta=f"{statut}")
                        
                        if statut == "Optimal":
                            st.balloons()
                            st.success(f"🌟 **VALEUR ÉLEVÉE** : {impact}")
                        else:
                            st.warning(f"⚠️ **NOTE** : {impact}")
                        
                        # Sauvegarde
                        db.execute_query("""INSERT INTO genomique 
                            (brebis_id, marqueur, homologie, prediction, categorie, date_test) 
                            VALUES (?,?,?,?,?,?)""",
                            (animal_id, gene_target, score, impact, cat, date.today()))
                    else:
                        st.error("Données manquantes pour l'analyse.")

    # --- MODULE : ÉTUDES DE VALEUR ---
    elif choice == "📉 Études de Valeur":
        st.title("📉 Synthèse Économique et Sanitaire")
        df_gen = db.fetch_all_as_df("SELECT * FROM genomique")
        
        if not df_gen.empty:
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("Répartition par Catégorie")
                st.plotly_chart(px.bar(df_gen, x="categorie", color="prediction", barmode="group"))
            with c2:
                st.subheader("Top Animaux par Gène")
                st.dataframe(df_gen.sort_values(by="homologie", ascending=False))
        else:
            st.info("Aucun test génomique enregistré.")

    # --- MODULE : DASHBOARD ---
    elif choice == "📊 Tableau de Bord":
        st.title("📊 Vue d'ensemble")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        if df_b.empty:
            if st.button("Initialiser données démo"):
                for i in range(5):
                    db.execute_query("INSERT INTO brebis (identifiant_unique, race, poids, created_at) VALUES (?,?,?,?)",
                                     (f"DZ-OVIS-{i}", "Ouled Djellal", 65+i, date.today()))
                st.rerun()
        else:
            st.write(f"Nombre d'animaux suivis : {len(df_b)}")
            st.dataframe(df_b)

if __name__ == "__main__":
    main()
