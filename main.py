"""
EXPERT OVIN DZ PRO - VERSION ELITE MASTER 2026
Système Intégral : Phénotypage, Bio-Informatique (Ensembl/EBI),
Moteur de Prédiction de Valeur (Éco/Santé) & Scanner Morphométrique.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
import re
import requests
from datetime import datetime, date, timedelta
import random
from typing import Dict, Optional, Any

# ============================================================================
# 1. DATABASE MASTER (ARCHITECTURE ROBUSTE)
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
            if "duplicate column name" not in str(e).lower():
                st.error(f"Erreur SQL: {e}")
            return None

    def fetch_all_as_df(self, query: str, params: tuple = ()):
        return pd.read_sql_query(query, self.conn, params=params)

def init_database(db: DatabaseManager):
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, sexe TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS genomique (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, marqueur TEXT, homologie REAL, 
            prediction TEXT, categorie TEXT, date_test DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, quantite_lait REAL
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)
    # Migration de sécurité
    db.execute_query("ALTER TABLE brebis ADD COLUMN sexe TEXT DEFAULT 'Femelle'")

# ============================================================================
# 2. MOTEURS IA & CONNECTEURS BIO-INFORMATIQUES PRO
# ============================================================================

class WebBioAPI:
    """Interface professionnelle avec Ensembl REST API."""
    BASE_URL = "https://rest.ensembl.org"
    SPECIES = "ovis_aries"

    @classmethod
    @st.cache_data(ttl=3600)
    def fetch_gene_metadata(cls, symbol: str) -> Optional[Dict[str, Any]]:
        endpoint = f"{cls.BASE_URL}/lookup/symbol/{cls.SPECIES}/{symbol}?"
        try:
            response = requests.get(endpoint, headers={"Content-Type": "application/json"}, timeout=10)
            return response.json() if response.status_code == 200 else None
        except Exception as e:
            st.error(f"Erreur API Metadata: {e}")
            return None

    @classmethod
    @st.cache_data(ttl=86400)
    def fetch_genomic_sequence(cls, gene_id: str) -> Optional[str]:
        endpoint = f"{cls.BASE_URL}/sequence/id/{gene_id}?type=genomic"
        try:
            response = requests.get(endpoint, headers={"Content-Type": "text/plain"}, timeout=15)
            return response.text if response.status_code == 200 else None
        except Exception as e:
            st.error(f"Erreur API Sequence: {e}")
            return None

class AIEngine:
    """Moteur d'IA pour le phénotypage et la valeur génétique."""
    
    GENE_CATALOG = {
        "CAST": {"nom": "Calpastatine", "cat": "Économique", "effet": "Tendreté de la viande."},
        "DGAT1": {"nom": "DGAT1", "cat": "Économique", "effet": "Richesse et rendement laitier."},
        "MSTN": {"nom": "Myostatine", "cat": "Économique", "effet": "Développement musculaire."},
        "PrP": {"nom": "Prion", "cat": "Sanitaire", "effet": "Résistance à la Tremblante."},
        "FecB": {"nom": "Booroola", "cat": "Sanitaire", "effet": "Prolificité/Fécondité."}
    }

    @staticmethod
    def calculate_genetic_homology(seq_ref: str, seq_user: str) -> Dict[str, Any]:
        s1 = re.sub(r'[^ATGC]', '', seq_ref.upper())
        s2 = re.sub(r'[^ATGC]', '', seq_user.upper())
        min_len = min(len(s1), len(s2))
        if min_len == 0: return {"homology": 0.0, "status": "Invalide"}
        
        matches = sum(1 for a, b in zip(s1[:min_len], s2[:min_len]) if a == b)
        score = round((matches / min_len) * 100, 2)
        
        status = "Optimal" if score >= 95 else "Standard" if score >= 85 else "Mutation"
        return {"homology": score, "status": status, "matches": matches}

    @staticmethod
    def nutrition_recommandee(poids: float) -> Dict[str, Any]:
        return {
            "Orge (kg)": round(poids * 0.012, 2),
            "Luzerne (kg)": round(poids * 0.025, 2),
            "CMV (g)": 35,
            "Eau (L)": round(poids * 0.1, 1)
        }

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    ia = AIEngine()
    api_web = WebBioAPI()

    st.sidebar.title("🐑 Bio-Master v2026")
    menu = [
        "📊 Dashboard Élite", "🧬 Analyse & Prédiction", "📈 GWAS & PLINK Pro", 
        "🌐 Recherche Bio-Web", "📝 Inscription", "📷 Scanner IA", 
        "🥛 Contrôle Laitier", "🌾 Nutrition Solo", "📈 Statistiques"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE : ANALYSE & PRÉDICTION (COEUR) ---
    if choice == "🧬 Analyse & Prédiction":
        st.title("🧬 Analyse Génomique & Valeur Prédite")
        c1, c2 = st.columns([1, 1])
        
        with c1:
            gene_target = st.selectbox("Gène cible", list(ia.GENE_CATALOG.keys()))
            st.info(f"**Impact :** {ia.GENE_CATALOG[gene_target]['effet']}")
            fasta_in = st.text_area("Séquence FASTA de l'animal", height=200)
            ani_df = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
            target_id = st.selectbox("Assigner à", ani_df['identifiant_unique'] if not ani_df.empty else ["Aucun"])

        with c2:
            if st.button("Lancer l'Alignement & Prédiction"):
                with st.spinner("Analyse Ensembl en cours..."):
                    meta = api_web.fetch_gene_metadata(gene_target)
                    if meta and fasta_in:
                        ref_seq = api_web.fetch_genomic_sequence(meta['id'])
                        res = ia.calculate_genetic_homology(ref_seq, fasta_in)
                        
                        st.metric("Homologie Génomique", f"{res['homology']}%", delta=res['status'])
                        
                        # Déduction de la valeur
                        cat = ia.GENE_CATALOG[gene_target]['cat']
                        impact = f"Haut potentiel {cat}" if res['homology'] > 90 else f"Valeur {cat} moyenne"
                        
                        if res['homology'] > 95: st.success(f"🌟 Individu Élite pour {gene_target}")
                        else: st.warning(f"Individu standard pour {gene_target}")
                        
                        db.execute_query("""INSERT INTO genomique 
                            (brebis_id, marqueur, homologie, prediction, categorie, date_test) 
                            VALUES (?,?,?,?,?,?)""",
                            (target_id, gene_target, res['homology'], impact, cat, date.today()))
                    else:
                        st.error("Erreur : Gène introuvable ou séquence vide.")

    # --- MODULE : RECHERCHE BIO-WEB ---
    elif choice == "🌐 Recherche Bio-Web":
        st.title("🌐 Explorateur Bio-Informatique")
        search_gene = st.text_input("Symbole du gène", "CAST")
        if st.button("Rechercher sur Ensembl"):
            data = api_web.fetch_gene_metadata(search_gene)
            if data:
                st.write(f"**Nom:** {data['display_name']} | **Chromosome:** {data['seq_region_name']}")
                st.json(data)
                
            else: st.error("Non trouvé.")

    # --- MODULE : DASHBOARD ---
    elif choice == "📊 Dashboard Élite":
        st.title("📊 Performance Globale du Troupeau")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df_b.empty:
            st.plotly_chart(px.pie(df_b, names='race', title="Composition du Troupeau"))
            df_g = db.fetch_all_as_df("SELECT * FROM genomique")
            if not df_g.empty:
                st.subheader("Analyse de Valeur Économique vs Sanitaire")
                st.plotly_chart(px.scatter(df_g, x="homologie", y="marqueur", color="categorie", size="homology"))
        else:
            st.info("Veuillez inscrire des animaux.")

    # --- MODULE : INSCRIPTION ---
    elif choice == "📝 Inscription":
        st.title("📝 Enregistrement Phénotypique")
        with st.form("reg"):
            uid = st.text_input("ID Unique (DZ-XXX)")
            race = st.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune"])
            poids = st.number_input("Poids (kg)", 10, 150, 60)
            if st.form_submit_button("Sauvegarder"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, race, poids, created_at) VALUES (?,?,?,?)", 
                                 (uid, race, poids, date.today()))
                st.success("Animal enregistré.")

    # --- MODULE : NUTRITION ---
    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Rationnement IA de Précision")
        p = st.number_input("Poids de l'animal (kg)", 20, 150, 65)
        st.table(pd.DataFrame([ia.nutrition_recommandee(p)]))

    # --- MODULE : SCANNER ---
    elif choice == "📷 Scanner IA":
        st.title("📷 Morphométrie par Image")
        st.info("Utilisez l'étalon de 1 mètre pour la calibration.")
        st.camera_input("Prendre une photo de profil")

    # --- MODULE : STATISTIQUES ---
    elif choice == "📈 Statistiques":
        st.title("📈 Analyse Data")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            st.plotly_chart(px.box(df, x="race", y="poids", points="all"))

if __name__ == "__main__":
    main()
