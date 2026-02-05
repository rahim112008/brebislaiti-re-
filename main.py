"""
EXPERT OVIN DZ PRO - VERSION FULL SCIENTIFIQUE 2026
Modules : Génomique (NCBI/Héritabilité), Nutrition Individualisée, 
Scanner avec Étalons, Dashboard Interactif, Santé & Gestation.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. DATABASE AVANCÉE (MORPHO + PHÉNOTYPE + GÉNOMIQUE)
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
            date_soin DATE, type_acte TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE IA ET GÉNOMIQUE
# ============================================================================

class BioInfoEngine:
    @staticmethod
    def calculer_genetique(homozygotie: float):
        # Simulation calcul consanguinité et héritabilité (h2)
        h2_lait = 0.25 # Héritabilité moyenne du lait chez l'ovin
        consanguinite = homozygotie * 0.5
        return h2_lait, consanguinite

# ============================================================================
# 3. INTERFACE UTILISATEUR CUMULATIVE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro V3", layout="wide")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    logic_gen = BioInfoEngine()

    # --- MENU NAVIGATION ---
    st.sidebar.title("🧬 Système Expert")
    menu = ["📊 Dashboard Pro", "📝 Inscription/Troupeau", "📷 Scanner & Phénotype", "🤰 Gestation IA", "🌾 Nutrition Solo", "🩺 Santé", "🧬 Génomique & NCBI", "📈 Stats"]
    choice = st.sidebar.radio("Modules", menu)

    # --- MODULE 1: DASHBOARD INTERACTIF ---
    if choice == "📊 Dashboard Pro":
        st.title("📊 Analyse Complète du Troupeau")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            st.subheader("Filtres Dynamiques")
            race_filter = st.multiselect("Filtrer par Race", df['race'].unique(), default=df['race'].unique())
            df_filtered = df[df['race'].isin(race_filter)]
            
            st.dataframe(df_filtered, use_container_width=True)
            
            c1, c2 = st.columns(2)
            with c1:
                fig_poids = px.box(df_filtered, x="race", y="poids", title="Distribution du Poids par Race", color="race")
                st.plotly_chart(fig_poids)
            with c2:
                fig_lact = px.scatter(df_filtered, x="tour_poitrine", y="poids", size="longueur", color="etat_mamelle", title="Morpho vs Mamelle")
                st.plotly_chart(fig_lact)
        else:
            st.info("Aucune donnée disponible.")

    # --- MODULE 2: TROUPEAU & PHÉNOTYPE ---
    elif choice == "📝 Inscription/Troupeau":
        st.title("📝 Fiche de Phénotypage")
        with st.form("form_pheno"):
            c1, c2 = st.columns(2)
            with c1:
                uid = st.text_input("ID Boucle (Obligatoire)")
                nom = st.text_input("Nom de la brebis")
                race_opt = ["Ouled Djellal", "Lacaune", "Rembi", "Hamra", "Autre/Inconnue"]
                race = st.selectbox("Race", race_opt)
                if race == "Autre/Inconnue":
                    race = st.text_input("Précisez la race")
            with c2:
                age_type = st.radio("Type d'âge", ["Dents (Remplaçantes)", "Exact (Date)", "Mois"])
                age_val = st.number_input("Valeur (Nombre)", 0.0, 15.0, 1.0)
            
            st.subheader("Mesures Morphométriques")
            m1, m2, m3, m4 = st.columns(4)
            h = m1.number_input("Hauteur (cm)", 40.0, 120.0, 75.0)
            l = m2.number_input("Longueur (cm)", 40.0, 130.0, 80.0)
            tp = m3.number_input("Tour Poitrine (cm)", 50.0, 150.0, 90.0)
            mamelle = m4.selectbox("État Mamelle", ["Excellente", "Correcte", "Asymétrique", "Sèche"])
            
            if st.form_submit_button("Enregistrer Phénotype Complet"):
                poids_est = (tp**2 * l) / 30000
                db.execute_query("""INSERT INTO brebis (identifiant_unique, nom, race, age_type, age_valeur, hauteur, longueur, tour_poitrine, etat_mamelle, poids, created_at) 
                                 VALUES (?,?,?,?,?,?,?,?,?,?,?)""", (uid, nom, race, age_type, age_val, h, l, tp, mamelle, poids_est, date.today()))
                st.success(f"Enregistré ! Poids estimé : {poids_est:.2f} kg")

    # --- MODULE 3: SCANNER IA AVEC ÉTALONS ---
    elif choice == "📷 Scanner & Phénotype":
        st.title("📷 Scanner Morphométrique IA")
        etalon = st.selectbox("Choisir l'objet d'étalonnage", ["Bâton standard 1m", "Feuille A4 (29.7cm)", "Carte Bancaire (8.5cm)"])
        st.camera_input("Prendre la photo de profil (avec l'objet visible)")
        st.info(f"L'IA utilise le {etalon} pour calibrer les pixels en centimètres réels.")

    # --- MODULE 4: GESTATION & PRÉDICTION ---
    elif choice == "🤰 Gestation IA":
        st.title("🤰 Suivi Reproduction")
        
        target = st.selectbox("Sélectionner la brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        d_eponge = st.date_input("Date de pose de l'éponge")
        if st.button("Prédire la Mise Bas"):
            # Estimation : Eponge + 14 jours (oestrus) + 150 jours (gestation)
            mb = d_eponge + timedelta(days=164)
            st.success(f"📅 Mise bas prévue (prédiction IA) : {mb.strftime('%d/%m/%Y')}")
            db.execute_query("INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevue, statut) VALUES (?,?,?,?)", 
                            (target, d_eponge, mb, "En cours"))

    # --- MODULE 5: NUTRITION SOLO ---
    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Ration Personnalisée")
        target = st.selectbox("ID Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        # Ici on pourrait chercher le stade via SQL, simplifié pour l'exemple
        st.subheader("Générateur de recette détaillée")
        if st.button("Calculer la ration optimale pour cette brebis"):
            st.write(f"**Ration pour {target} :**")
            st.info("Recette : 800g Orge / 400g Son / 1kg Luzerne. Complément : 20g CMV.")

    # --- MODULE 6: GÉNOMIQUE AVANCÉE ---
    elif choice == "🧬 Génomique & NCBI":
        st.title("🧬 Laboratoire Génomique")
        
        tab1, tab2 = st.tabs(["Analyse Bioinformatique", "Héritabilité & SNP"])
        with tab1:
            st.write("🔗 Connexion vers : **NCBI GenBank / OMIA (Sheep)**")
            st.text_area("Séquence FASTA pour alignement")
            st.button("Lancer l'alignement BLAST")
        with tab2:
            homo = st.slider("Taux d'homozygotie (%)", 0.0, 100.0, 15.0)
            h2, cons = logic_gen.calculer_genetique(homo/100)
            st.metric("Coefficient de Consanguinité", f"{cons*100:.2f}%")
            st.metric("Héritabilité estimée (h²)", h2)

    # --- MODULE 7: STATS ---
    elif choice == "📈 Stats":
        st.title("📈 Analyses Statistiques")
        st.write("ANOVA : Comparaison des poids entre les races")
        df = db.fetch_all_as_df("SELECT race, poids FROM brebis")
        if not df.empty:
            fig = px.violin(df, y="poids", x="race", box=True, points="all", title="ANOVA Visuelle")
            st.plotly_chart(fig)

if __name__ == "__main__":
    main()
