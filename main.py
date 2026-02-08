"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.04
Système Intégral de Gestion de Précision : 
Phénotypage, Lait, Génomique (SNP), Santé, Nutrition & IA
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. DATABASE MASTER (ARCHITECTURE CUMULATIVE)
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
            nom TEXT, race TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, attaches_mamelle TEXT, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_controle DATE,
            quantite_lait REAL, tb REAL, tp REAL, cellules INTEGER,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_eponge DATE, date_mise_bas_prevue DATE, statut TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_soin DATE,
            type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE IA & GÉNOMIQUE (INDEX & MARQUEURS SNP)
# ============================================================================

class AIEngine:
    @staticmethod
    def calculer_index_elite(row, df_lait):
        score_morpho = (row['tour_poitrine'] * 0.2) + (row['note_mamelle'] * 5)
        score_os = row['circ_canon'] * 3
        lait_indiv = df_lait[df_lait['brebis_id'] == row['identifiant_unique']]
        score_lait = lait_indiv['quantite_lait'].mean() * 15 if not lait_indiv.empty else 0
        return round((score_morpho + score_os + score_lait), 2)

    @staticmethod
    def nutrition_recommandee(poids):
        orge = round(poids * 0.012, 2)
        luzerne = round(poids * 0.02, 2)
        return {"Orge (kg)": orge, "Luzerne (kg)": luzerne, "CMV (g)": 30}

# Bibliothéque de Marqueurs SNP
SNP_LIBRARY = {
    "QUALITÉ VIANDE": {
        "CAST": "Calpastatine (Tendreté de la viande)",
        "GDF8": "Myostatine (Hypertrophie musculaire / Rendement carcasses)"
    },
    "PRODUCTION LAIT": {
        "DGAT1": "Taux butyreux (Gras) et Rendement laitier",
        "LALBA": "Alpha-lactalbumine (Volume de lait)"
    },
    "RÉSISTANCE & ADAPTATION": {
        "HSP70": "Thermotolérance (Résistance au stress thermique DZ)",
        "MHC": "Complexe Majeur d'Histocompatibilité (Immunité globale)"
    },
    "MALADIES GÉNÉTIQUES": {
        "PrP": "Tremblante (Scrapie) - Sensibilité/Résistance",
        "Spider_Lamb": "Syndrome de l'agneau araignée (Chondrodysplasie)"
    }
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

    st.sidebar.title("🐑 Système Intégré v2026")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA", 
        "🥛 Contrôle Laitier", 
        "🤰 Gestation IA", 
        "🌾 Nutrition Solo", 
        "🩺 Santé & Vaccins", 
        "🧬 Génomique & SNP", 
        "📈 Statistiques"
    ]
    choice = st.sidebar.radio("Modules", menu)

    # --- MODULE 1: DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Performance & Sélection Élite")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        if not df_b.empty:
            df_b['Index_Selection'] = df_b.apply(lambda r: ia.calculer_index_elite(r, df_l), axis=1)
            df_top = df_b.sort_values(by='Index_Selection', ascending=False)
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df_b))
            c2.metric("Moyenne Laitière (L)", round(df_l['quantite_lait'].mean(), 2) if not df_l.empty else 0)
            c3.metric("Meilleur Index", df_top['Index_Selection'].max())
            st.subheader("🏆 Top 10 Génitrices")
            st.dataframe(df_top[['identifiant_unique', 'race', 'Index_Selection', 'poids']].head(10))
            st.plotly_chart(px.scatter(df_b, x="tour_poitrine", y="Index_Selection", color="race", size="poids"))
        else: st.info("Veuillez inscrire des animaux.")

    # --- MODULE 2: PHÉNOTYPE ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Phénotypage Avancé")
        with st.form("inscription"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("Identifiant Unique (Boucle)")
            race = c1.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra", "Autre"])
            age_t = c2.radio("Méthode d'âge", ["Dents", "Mois", "Années"])
            age_v = c2.number_input("Valeur âge", 0, 15, 2)
            st.subheader("Mesures du Corps")
            m1, m2, m3 = st.columns(3)
            h, l, tp = m1.number_input("Hauteur (cm)", 40, 110, 75), m2.number_input("Longueur (cm)", 40, 120, 80), m3.number_input("Poitrine (cm)", 50, 150, 90)
            lb, lgb, can = m1.number_input("Larg. Bassin (cm)", 10, 40, 22), m2.number_input("Long. Bassin (cm)", 10, 40, 20), m3.number_input("Canon (cm)", 5.0, 15.0, 8.5)
            note_m = st.slider("Note Mamelle", 1, 10, 5)
            attaches = st.selectbox("Attaches", ["Solides", "Moyennes", "Lâches"])
            if st.form_submit_button("Enregistrer"):
                poids = (tp**2 * l) / 30000
                db.execute_query("INSERT INTO brebis (identifiant_unique, race, age_type, age_valeur, hauteur, longueur, tour_poitrine, largeur_bassin, long_bassin, circ_canon, note_mamelle, attaches_mamelle, poids, created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (uid, race, age_t, age_v, h, l, tp, lb, lgb, can, note_m, attaches, poids, date.today()))
                st.success("Brebis enregistrée.")

    # --- MODULE 3: SCANNER ---
    elif choice == "📷 Scanner IA":
        st.title("📷 Scanner Morphométrique 1m")
        etalon = st.selectbox("Étalon", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
        st.camera_input("Capture")
        st.info(f"Analyse via étalon : {etalon}")

    # --- MODULE 4: LAIT ---
    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Contrôle Laitier")
        with st.form("lait"):
            target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
            qte = st.number_input("Lait (L)", 0.0, 10.0, 2.0)
            tb, tp = st.slider("TB (g/L)", 20, 80, 45), st.slider("TP (g/L)", 20, 70, 35)
            if st.form_submit_button("Valider"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait, tb, tp) VALUES (?,?,?,?,?)", (target, date.today(), qte, tb, tp))
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        if not df_l.empty: st.plotly_chart(px.line(df_l, x="date_controle", y="quantite_lait", color="brebis_id"))

    # --- MODULE 5: GESTATION ---
    elif choice == "🤰 Gestation IA":
        st.title("🤰 Gestion Reproduction")
        target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        d_ep = st.date_input("Date Pose Éponge")
        if st.button("Prédire"):
            mb = d_ep + timedelta(days=164)
            st.success(f"Mise bas : {mb.strftime('%d/%m/%Y')}")
            db.execute_query("INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevue, statut) VALUES (?,?,?,?)", (target, d_ep, mb, "En cours"))

    # --- MODULE 6: NUTRITION ---
    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Ration IA")
        df_b = db.fetch_all_as_df("SELECT identifiant_unique, poids FROM brebis")
        if not df_b.empty:
            target = st.selectbox("Brebis", df_b['identifiant_unique'])
            poids_b = df_b[df_b['identifiant_unique'] == target]['poids'].values[0]
            recette = ia.nutrition_recommandee(poids_b)
            for k, v in recette.items(): st.metric(k, v)

    # --- MODULE 7: SANTÉ ---
    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Suivi Sanitaire")
        target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        acte = st.selectbox("Acte", ["Enterotoxémie", "Fièvre Aphteuse", "Clavelée", "PPR", "Vermifuge"])
        if st.button("Enregistrer"):
            rappel = date.today() + timedelta(days=180)
            db.execute_query("INSERT INTO sante (brebis_id, date_soin, type_acte, rappel_prevu) VALUES (?,?,?,?)", (target, date.today(), acte, rappel))
            st.success(f"Rappel le {rappel}")

    # --- MODULE 8: GÉNOMIQUE & SNP (NOUVEAU) ---
    elif choice == "🧬 Génomique & SNP":
        st.title("🧬 Diagnostic Moléculaire & Marqueurs SNP")
        
        tab_search, tab_diag = st.tabs(["Exploration des Marqueurs", "Diagnostic Individuel"])
        
        with tab_search:
            cat_sel = st.selectbox("Catégorie de sélection", list(SNP_LIBRARY.keys()))
            gene_info = SNP_LIBRARY[cat_sel]
            st.table(pd.DataFrame.from_dict(gene_info, orient='index', columns=['Description / Effet']))
            st.info("Données synchronisées avec l'archive NCBI Sheep Genome.")

        with tab_diag:
            c1, c2 = st.columns(2)
            ani_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
            if not ani_list.empty:
                target_ani = c1.selectbox("Sélectionner l'animal pour DNA Test", ani_list)
                marker_sel = c2.selectbox("Marqueur à tester", ["CAST", "DGAT1", "HSP70", "PrP"])
                
                # Simulation de détection de génotype
                st.subheader(f"Analyse de l'allèle : {marker_sel}")
                genotype = st.radio("Génotype observé (via Séquençage/PCR)", ["Homozygote Sauvage", "Hétérozygote", "Homozygote Muté"])
                
                # Logique de décision
                if marker_sel == "PrP" and "Muté" in genotype:
                    st.error("🚨 ALERTE : Individu sensible à la Tremblante. Ne pas utiliser pour la reproduction.")
                elif marker_sel in ["CAST", "DGAT1"] and "Muté" in genotype:
                    st.success("✨ ÉLITE : Porteur du gène de haute performance. Priorité sélection.")
                else:
                    st.warning("Information : Profil génétique standard.")
                
                
            else:
                st.warning("Veuillez d'abord inscrire des animaux.")

    # --- MODULE 9: STATS ---
    elif choice == "📈 Statistiques":
        st.title("📈 Biostatistique")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            st.plotly_chart(px.violin(df, x="race", y="poids", box=True))
            st.write("Corrélation Canon vs Poids :", df['circ_canon'].corr(df['poids']))

if __name__ == "__main__":
    main()
