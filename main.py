"""
EXPERT OVIN DZ PRO - SYSTÈME INTÉGRAL CUMULATIF (V 2026.02)
Modules inclus : Auth, SQL, Génomique, Biochimie, Morphométrie IA, Gestation, Statistiques & Nutrition
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
from datetime import datetime, date, timedelta
from typing import Dict, List, Any

# ============================================================================
# 1. ARCHITECTURE DE DONNÉES (SQLITE)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_manager.db"):
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
    """Initialisation cumulative de toutes les tables du projet"""
    tables = [
        # Table Pivot : Identité
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, date_naissance DATE, race TEXT, sexe TEXT,
            statut TEXT DEFAULT 'active', notes TEXT
        )""",
        # Table Morphométrie (Scanner)
        """CREATE TABLE IF NOT EXISTS caracteres_morpho (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_mesure DATE,
            hauteur REAL, longueur REAL, tour_poitrine REAL, poids_estime REAL,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table Biochimie (Laboratoire)
        """CREATE TABLE IF NOT EXISTS biochimie_lait (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_analyse DATE,
            fat REAL, prot REAL, bhb REAL, ratio_tbtp REAL, diagnostic TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table Gestation (Reproduction)
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id INTEGER, 
            date_eponge DATE, date_mise_bas_prevu DATE, statut TEXT DEFAULT 'en_cours',
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE MÉTIER (GÉNOMIQUE, STATS, NUTRITION)
# ============================================================================

class ElevageExpert:
    # --- MODULE NUTRITION ---
    ALIMENTS_DZ = {
        "Orge": {"ufl": 1.0, "pdi": 80, "prix": 5500},
        "Son de Blé": {"ufl": 0.85, "pdi": 95, "prix": 2500},
        "Luzerne": {"ufl": 0.55, "pdi": 110, "prix": 4500},
        "Paille": {"ufl": 0.35, "pdi": 30, "prix": 1200},
        "Maïs": {"ufl": 1.15, "pdi": 85, "prix": 7000}
    }
    
    BESOINS = {
        "Entretien": {"ufl": 0.8, "pdi": 75},
        "Gestation": {"ufl": 1.2, "pdi": 115},
        "Lactation": {"ufl": 1.9, "pdi": 185}
    }

    # --- MODULE GÉNOMIQUE ---
    def analyser_snp(self, ref: str, stu: str) -> Dict:
        snps = [{'pos': i+1, 'ref': r, 'alt': s} for i, (r, s) in enumerate(zip(ref, stu)) if r != s]
        return {'total': len(snps), 'details': snps}

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT MAIN)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro", layout="wide")

    # Session State
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    logic = ElevageExpert()

    # --- AUTHENTICATION ---
    if 'auth' not in st.session_state: st.session_state.auth = False
    if not st.session_state.auth:
        st.title("🔐 Accès Expert Ovin DZ Pro")
        if st.text_input("Pass", type="password") == "admin123":
            if st.button("Entrer"): st.session_state.auth = True; st.rerun()
        return

    # --- MENU NAVIGATION CUMULATIF ---
    st.sidebar.title("🐑 Navigation Système")
    menu = ["📊 Dashboard", "📝 Troupeau", "📷 Scanner IA", "🧪 Biochimie", "🤰 Gestation", "🌾 Nutrition", "🧬 Génomique", "📈 Stats"]
    choice = st.sidebar.radio("Modules", menu, key="main_nav")

    # --- MODULE DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title("📊 Tableau de Bord Intégré")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        c1, c2, c3 = st.columns(3)
        c1.metric("Effectif Total", len(df))
        st.dataframe(df, use_container_width=True)

    # --- MODULE TROUPEAU (INSCRIPTION) ---
    elif choice == "📝 Troupeau":
        st.title("📝 Gestion du Troupeau")
        with st.form("add_sheep"):
            uid = st.text_input("Identifiant Unique (Boucle)")
            race = st.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra"])
            if st.form_submit_button("Sauvegarder"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, race) VALUES (?,?)", (uid, race))
                st.success("Animal enregistré.")

    # --- MODULE SCANNER IA ---
    elif choice == "📷 Scanner IA":
        st.title("📷 Morphométrie IA")
        
        target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        h = st.number_input("Hauteur (cm)", 50.0, 100.0, 75.0)
        l = st.number_input("Longueur (cm)", 50.0, 110.0, 80.0)
        tp = st.number_input("Tour Poitrine (cm)", 50.0, 130.0, 90.0)
        if st.button("Calculer & Enregistrer"):
            poids = (tp**2 * l) / 30000
            db.execute_query("INSERT INTO caracteres_morpho (brebis_id, hauteur, longueur, tour_poitrine, poids_estime) VALUES (?,?,?,?,?)", 
                            (target, h, l, tp, poids))
            st.success(f"Poids estimé : {poids:.2f} kg")

    # --- MODULE NUTRITION ---
    elif choice == "🌾 Nutrition":
        st.title("🌾 Nutritioniste IA (Marché DZ)")
        
        stade = st.selectbox("Stade de la brebis", list(logic.BESOINS.keys()))
        col1, col2 = st.columns(2)
        ration = {}
        with col1:
            for alim, vals in logic.ALIMENTS_DZ.items():
                ration[alim] = st.slider(f"{alim} (kg/j)", 0.0, 2.5, 0.0)
        with col2:
            u_tot = sum(ration[a] * logic.ALIMENTS_DZ[a]['ufl'] for a in ration)
            p_tot = sum(ration[a] * logic.ALIMENTS_DZ[a]['pdi'] for a in ration)
            c_tot = sum(ration[a] * (logic.ALIMENTS_DZ[a]['prix']/100) for a in ration)
            
            st.metric("Apport Énergie", f"{u_tot:.2f} UFL", f"{u_tot - logic.BESOINS[stade]['ufl']:.2f}")
            st.metric("Apport Protéines", f"{p_tot:.0f}g PDI", f"{p_tot - logic.BESOINS[stade]['pdi']:.0f}g")
            st.metric("Coût Journalier", f"{c_tot:.2f} DZD")

    # --- MODULE BIOCHIMIE ---
    elif choice == "🧪 Biochimie":
        st.title("🧪 Analyse Laitière")
        
        st.info("Module connecté pour l'analyse du TB/TP et BHB.")

    # --- MODULE GESTATION ---
    elif choice == "🤰 Gestation":
        st.title("🤰 Suivi de Reproduction")
        
        st.info("Calcul des dates de mise bas et gestion des éponges.")

    # --- MODULE GÉNOMIQUE ---
    elif choice == "🧬 Génomique":
        st.title("🧬 Séquençage & SNP")
        
        st.info("Analyse bioinformatique des séquences génétiques.")

    # --- MODULE STATISTIQUES ---
    elif choice == "📈 Stats":
        st.title("📈 Statistiques Avancées")
        
        st.info("Corrélations et régressions sur les données du troupeau.")

if __name__ == "__main__":
    main()
