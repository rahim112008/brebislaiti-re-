"""
EXPERT OVIN DZ PRO - VERSION ELITE 2026.02
Système : Bio-Informatique, API Ensembl/EBI & Diagnostic de Valeur.
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
from typing import Dict, Optional, Any

# ============================================================================
# 1. DATABASE & CONFIGURATION
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
            brebis_id TEXT, marqueur TEXT, homologie REAL, 
            prediction TEXT, categorie TEXT, date_test DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. MOTEURS IA & CONNECTEURS BIO-INFORMATIQUES PRO
# ============================================================================

class WebBioAPI:
    BASE_URL = "https://rest.ensembl.org"
    SPECIES = "ovis_aries"

    @classmethod
    @st.cache_data(ttl=3600)
    def fetch_gene_metadata(cls, symbol: str) -> Optional[Dict[str, Any]]:
        endpoint = f"{cls.BASE_URL}/lookup/symbol/{cls.SPECIES}/{symbol}?"
        try:
            r = requests.get(endpoint, headers={"Content-Type": "application/json"}, timeout=10)
            return r.json() if r.ok else None
        except: return None

    @classmethod
    @st.cache_data(ttl=86400)
    def fetch_genomic_sequence(cls, gene_id: str) -> Optional[str]:
        endpoint = f"{cls.BASE_URL}/sequence/id/{gene_id}?type=genomic"
        try:
            r = requests.get(endpoint, headers={"Content-Type": "text/plain"}, timeout=15)
            return r.text if r.ok else None
        except: return None

class AIEngine:
    GENE_CATALOG = {
        "CAST": {"nom": "Calpastatine", "cat": "Économique", "effet": "Tendreté de la viande."},
        "DGAT1": {"nom": "DGAT1", "cat": "Économique", "effet": "Richesse laitière (gras)."},
        "MSTN": {"nom": "Myostatine", "cat": "Économique", "effet": "Développement musculaire."},
        "PrP": {"nom": "Prion", "cat": "Sanitaire", "effet": "Résistance à la Tremblante."},
        "FecB": {"nom": "Booroola", "cat": "Sanitaire", "effet": "Taux de gémellité."}
    }

    @staticmethod
    def calculate_homology(seq_ref: str, seq_user: str) -> Dict[str, Any]:
        s1 = re.sub(r'[^ATGC]', '', seq_ref.upper())
        s2 = re.sub(r'[^ATGC]', '', seq_user.upper())
        min_len = min(len(s1), len(s2))
        if min_len == 0: return {"score": 0.0, "status": "Invalide"}
        
        matches = sum(1 for a, b in zip(s1[:min_len], s2[:min_len]) if a == b)
        score = round((matches / min_len) * 100, 2)
        return {"score": score, "analyzed": min_len}

# ============================================================================
# 3. INTERFACE UTILISATEUR
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
        "📊 Dashboard Élite", 
        "🧬 Génomique & Alignement", # Restauré ici
        "🌐 Recherche Bio-Web", 
        "📝 Inscription", 
        "📷 Scanner IA", 
        "🌾 Nutrition Solo"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE : GÉNOMIQUE & ALIGNEMENT (RESTAURÉ) ---
    if choice == "🧬 Génomique & Alignement":
        st.title("🧬 Diagnostic Génomique & Alignement")
        
        col_in, col_res = st.columns([1, 1])
        
        with col_in:
            gene_choice = st.selectbox("Sélectionner le gène d'intérêt", list(ia.GENE_CATALOG.keys()))
            st.info(f"**Cible :** {ia.GENE_CATALOG[gene_choice]['nom']} | **Effet :** {ia.GENE_CATALOG[gene_choice]['effet']}")
            
            fasta_input = st.text_area("Séquence ADN (Format FASTA)", height=200, help="Collez ici la séquence brute ATGC")
            
            ani_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
            target_id = st.selectbox("Assigner l'analyse à :", ani_list['identifiant_unique'] if not ani_list.empty else ["Aucun animal"])

        with col_res:
            if st.button("Lancer l'Alignement Professionnel"):
                with st.spinner("Interrogation d'Ensembl REST API..."):
                    # 1. Chercher métadonnées
                    meta = api_web.fetch_gene_metadata(gene_choice)
                    if meta:
                        # 2. Chercher séquence de référence
                        ref_seq = api_web.fetch_genomic_sequence(meta['id'])
                        if ref_seq and fasta_input:
                            # 3. Calculer Homologie
                            res = ia.calculate_homology(ref_seq, fasta_input)
                            
                            st.subheader("Résultats de l'Analyse")
                            st.metric("Homologie Génomique", f"{res['score']}%")
                            
                            # Prédiction IA
                            cat = ia.GENE_CATALOG[gene_choice]['cat']
                            prediction = "Élite" if res['score'] > 98 else "Standard"
                            
                            if res['score'] > 95:
                                st.success(f"✅ Potentiel {cat} confirmé.")
                            else:
                                st.warning(f"⚠️ Variations détectées par rapport à la référence.")
                            
                            # Sauvegarde
                            db.execute_query("""INSERT INTO genomique 
                                (brebis_id, marqueur, homologie, prediction, categorie, date_test) 
                                VALUES (?,?,?,?,?,?)""",
                                (target_id, gene_choice, res['score'], prediction, cat, date.today()))
                            
                            
                        else:
                            st.error("Séquence de référence ou utilisateur introuvable.")
                    else:
                        st.error("Gène introuvable sur les serveurs internationaux.")

    # --- AUTRES MODULES ---
    elif choice == "🌐 Recherche Bio-Web":
        st.title("🌐 Consultation des Bases Mondiales")
        gene_name = st.text_input("Symbole (ex: CAST)", "CAST")
        if st.button("Chercher"):
            data = api_web.fetch_gene_metadata(gene_name)
            if data: st.json(data)

    elif choice == "📊 Dashboard Élite":
        st.title("📊 Statut du Troupeau")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df_b.empty:
            st.dataframe(df_b)
            df_g = db.fetch_all_as_df("SELECT * FROM genomique")
            if not df_g.empty:
                st.plotly_chart(px.bar(df_g, x="brebis_id", y="homologie", color="marqueur", title="Profil Génomique"))
        else:
            st.info("Base de données vide. Utilisez l'onglet Inscription.")

    elif choice == "📝 Inscription":
        st.title("📝 Nouvel Animal")
        with st.form("f"):
            id_u = st.text_input("ID Unique")
            race = st.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra"])
            poids = st.number_input("Poids (kg)", 10, 150, 60)
            if st.form_submit_button("Ajouter"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, race, poids, created_at) VALUES (?,?,?,?)", (id_u, race, poids, date.today()))
                st.success("Ajouté !")

    elif choice == "📷 Scanner IA":
        st.title("📷 Scanner Morphométrique")
        st.camera_input("Scanner (Calibration 1m)")

    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Rationnement")
        p = st.number_input("Poids", 20, 150, 60)
        st.write(f"Orge: {p*0.012}kg | Luzerne: {p*0.025}kg")

if __name__ == "__main__":
    main()
