"""
EXPERT OVIN DZ PRO - VERSION FULL SCIENTIFIQUE 2026
Système Intégral : Phénotypage, Santé IA, Génomique (NCBI) & Stats Avancées
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. ARCHITECTURE DE DONNÉES (SQLITE)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_expert_v3.db"):
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
            nom TEXT, race TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            etat_mamelle TEXT, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_eponge DATE, date_mise_bas_prevue DATE, statut TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE, notes TEXT
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE SCIENTIFIQUE (GÉNOMIQUE & SANTÉ)
# ============================================================================

class ElevageExpertIA:
    PROTOCOLES_VACCINS = {
        "Enterotoxémie (Braxy)": 180,
        "Fièvre Aphteuse": 365,
        "Clavelée": 365,
        "PPR (Peste)": 1095,
        "Déparasitage": 120
    }

    @staticmethod
    def calculer_genetique(homo_val: float):
        # h2 (Héritabilité), cons (Consanguinité), hetero (Hétérozygotie)
        h2_lait = 0.28
        consanguinite = homo_val * 0.5
        heterozygotie = 1.0 - homo_val
        return h2_lait, consanguinite, heterozygotie

# ============================================================================
# 3. INTERFACE UTILISATEUR CUMULATIVE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro V3", layout="wide")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    logic = ElevageExpertIA()

    # --- MENU NAVIGATION ---
    st.sidebar.title("🧬 Système Expert")
    menu = ["📊 Dashboard Pro", "📝 Inscription/Troupeau", "📷 Scanner & Phénotype", "🤰 Gestation IA", "🌾 Nutrition Solo", "🩺 Santé & Vaccins", "🧬 Génomique & NCBI", "📈 Stats"]
    choice = st.sidebar.radio("Modules", menu)

    # --- MODULE 1: DASHBOARD ---
    if choice == "📊 Dashboard Pro":
        st.title("📊 Analyse Complète du Troupeau")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            st.metric("Total Brebis", len(df))
            st.dataframe(df, use_container_width=True)
            
            c1, c2 = st.columns(2)
            with c1:
                st.plotly_chart(px.box(df, x="race", y="poids", color="race", title="Poids par Race"))
            with c2:
                st.plotly_chart(px.scatter(df, x="age_valeur", y="poids", color="etat_mamelle", title="Croissance vs Mamelle"))
        else:
            st.info("Base de données vide.")

    # --- MODULE 2: TROUPEAU & PHÉNOTYPE ---
    elif choice == "📝 Inscription/Troupeau":
        st.title("📝 Fiche de Phénotypage")
        with st.form("form_pheno"):
            c1, c2 = st.columns(2)
            with c1:
                uid = st.text_input("ID Boucle (Obligatoire)")
                race_opt = ["Ouled Djellal", "Lacaune", "Rembi", "Hamra", "Autre/Inconnue"]
                race = st.selectbox("Race", race_opt)
                if race == "Autre/Inconnue": race = st.text_input("Précisez la race")
            with c2:
                age_type = st.radio("Mesure d'âge", ["Dents (Remplaçantes)", "Mois", "Exact (Années)"])
                age_val = st.number_input("Valeur", 0.0, 15.0, 1.0)
            
            st.subheader("Morphométrie Laitière")
            m1, m2, m3, m4 = st.columns(4)
            h = m1.number_input("Hauteur (cm)", 40.0, 120.0, 75.0)
            l = m2.number_input("Longueur (cm)", 40.0, 130.0, 80.0)
            tp = m3.number_input("Tour Poitrine (cm)", 50.0, 150.0, 90.0)
            mamelle = m4.selectbox("Type Mamelle", ["Globuleuse", "Pire", "Asymétrique", "Sèche"])
            
            if st.form_submit_button("Sauvegarder l'animal"):
                poids_est = (tp**2 * l) / 30000
                db.execute_query("""INSERT INTO brebis (identifiant_unique, race, age_type, age_valeur, hauteur, longueur, tour_poitrine, etat_mamelle, poids, created_at) 
                                 VALUES (?,?,?,?,?,?,?,?,?,?)""", (uid, race, age_type, age_val, h, l, tp, mamelle, poids_est, date.today()))
                st.success("Données enregistrées.")

    # --- MODULE 3: SCANNER IA ---
    elif choice == "📷 Scanner & Phénotype":
        st.title("📷 Scanner Morphométrique IA")
        etalon = st.selectbox("Objet de référence", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
        st.camera_input("Scanner l'animal de profil")
        st.info(f"Calibration IA active via : {etalon}")

    # --- MODULE 4: GESTATION IA ---
    elif choice == "🤰 Gestation IA":
        st.title("🤰 Suivi Reproduction")
        df_b = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        target = st.selectbox("ID Brebis", df_b)
        d_eponge = st.date_input("Date de l'éponge")
        if st.button("Calculer Mise Bas"):
            mb = d_eponge + timedelta(days=164)
            st.success(f"Mise bas prévue : {mb.strftime('%d/%m/%Y')}")
            db.execute_query("INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevue, statut) VALUES (?,?,?,?)", 
                            (target, d_eponge, mb, "En cours"))

    # --- MODULE 5: NUTRITION SOLO ---
    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Ration Individualisée")
        df_b = db.fetch_all_as_df("SELECT identifiant_unique, poids FROM brebis")
        target = st.selectbox("Sélectionner la brebis", df_b['identifiant_unique'])
        if st.button("Générer Recette"):
            # Simulation ration
            st.info(f"Ration optimisée pour {target} : 900g Orge, 1.2kg Luzerne, 30g CMV.")

    # --- MODULE 6: SANTÉ & VACCINS (RÉINTÉGRÉ) ---
    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Gestion Sanitaire & IA")
        tab1, tab2 = st.tabs(["💉 Nouvel Acte", "🗓️ Rappels IA"])
        with tab1:
            with st.form("vax"):
                target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
                acte = st.selectbox("Maladie/Vaccin", list(logic.PROTOCOLES_VACCINS.keys()))
                prod = st.text_input("Nom du produit")
                d_acte = st.date_input("Date du jour")
                if st.form_submit_button("Enregistrer"):
                    rappel = d_acte + timedelta(days=logic.PROTOCOLES_VACCINS[acte])
                    db.execute_query("INSERT INTO sante (brebis_id, date_soin, type_acte, produit, rappel_prevu) VALUES (?,?,?,?,?)",
                                    (target, d_acte, acte, prod, rappel))
                    st.success(f"Rappel planifié pour le {rappel}")
        with tab2:
            st.dataframe(db.fetch_all_as_df("SELECT * FROM sante ORDER BY rappel_prevu ASC"))

    # --- MODULE 7: GÉNOMIQUE AVANCÉE ---
    elif choice == "🧬 Génomique & NCBI":
        st.title("🧬 Lab Génomique")
        st.write("Alignement GeneBank/NCBI (Simulation)")
        fasta = st.text_area("Séquence SNP/ADN (FASTA)")
        homo = st.slider("Taux Homozygotie (%)", 0.0, 100.0, 10.0)
        h2, cons, hetero = logic.calculer_genetique(homo/100)
        
        c1, c2, c3 = st.columns(3)
        c1.metric("Consanguinité", f"{cons*100:.2f}%")
        c2.metric("Héritabilité (h²)", h2)
        c3.metric("Hétérozygotie", f"{hetero*100:.1f}%")

    # --- MODULE 8: STATS ---
    elif choice == "📈 Stats":
        st.title("📈 Statistiques & ANOVA")
        df = db.fetch_all_as_df("SELECT race, poids FROM brebis")
        if not df.empty:
            st.plotly_chart(px.violin(df, y="poids", x="race", title="Analyse de Variance (ANOVA Visuelle)"))

if __name__ == "__main__":
    main()
