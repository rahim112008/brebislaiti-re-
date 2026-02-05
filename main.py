"""
EXPERT OVIN DZ PRO - VERSION ULTIMATE MASTER 2026
------------------------------------------------
Auteur: Expert Ovin DZ & Gemini
Modules: Génétique, Biochimie, Nutrition, Finance, Stocks, Simulation
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date

# ============================================================================
# 1. ARCHITECTURE DE LA BASE DE DONNÉES
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_ultimate.db"):
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
        "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT)",
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE NOT NULL,
            owner_id TEXT, race TEXT, type_animal TEXT, poids REAL, note_mamelle REAL, 
            pere_id TEXT, mere_id TEXT, snp_gene TEXT, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS lait_biochimie (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            qte_lait REAL, tb REAL, tp REAL, esd REAL, date_controle DATE, owner_id TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS stocks (
            owner_id TEXT, aliment TEXT, quantite_q REAL, prix_q REAL, PRIMARY KEY(owner_id, aliment)
        )"""
    ]
    for t in tables: db.execute_query(t)
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert')")
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('eleveur1', 'ovin2026', 'Eleveur')")

# ============================================================================
# 2. MOTEUR SCIENTIFIQUE (NUTRITION & GÉNÉTIQUE)
# ============================================================================

TABLE_VALEURS = {
    "Orge": {"UFL": 1.0, "PDI": 80},
    "Son de blé": {"UFL": 0.85, "PDI": 95},
    "Foin de luzerne": {"UFL": 0.65, "PDI": 90},
    "Foin d'avoine": {"UFL": 0.55, "PDI": 50},
    "Maïs grain": {"UFL": 1.15, "PDI": 75},
    "Paille de blé": {"UFL": 0.35, "PDI": 30}
}

class ScienceEngine:
    @staticmethod
    def calculer_isg(row):
        """Index de Sélection Global (Base 100)"""
        m = row['note_mamelle'] if not pd.isna(row.get('note_mamelle')) else 5.0
        p = row['poids'] if not pd.isna(row.get('poids')) else 50.0
        return round((m * 6) + (p * 0.4), 1)

    @staticmethod
    def predire_agnelle(agnelle_row, df_troupeau):
        pere = df_troupeau[df_troupeau['identifiant_unique'] == agnelle_row['pere_id']]
        mere = df_troupeau[df_troupeau['identifiant_unique'] == agnelle_row['mere_id']]
        score_p = 50 if pere.empty else ScienceEngine.calculer_isg(pere.iloc[0])
        score_m = 50 if mere.empty else ScienceEngine.calculer_isg(mere.iloc[0])
        return round((score_p * 0.5) + (score_m * 0.5), 1)

    @staticmethod
    def calculer_esd(tb, tp):
        return round((tp * 0.85) + (tb * 0.1) + 5.5, 2)

# ============================================================================
# 3. INTERFACE UTILISATEUR
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin Master Pro", layout="wide", page_icon="🐑")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    
    db = st.session_state.db

    # --- Authentification ---
    if 'auth' not in st.session_state:
        st.title("🛡️ Station Master Ovin DZ")
        with st.form("login_form"):
            u = st.text_input("Identifiant")
            p = st.text_input("Mot de passe", type="password")
            if st.form_submit_button("Entrer dans le système"):
                res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
                if res:
                    st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                    st.rerun()
                else: st.error("Accès refusé.")
        return

    user, role = st.session_state.username, st.session_state.role

    # --- Barre Latérale & Navigation ---
    st.sidebar.title(f"👤 {user} ({role})")
    
    if role == "Expert":
        eleveurs = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
        view_user = st.sidebar.selectbox("📂 Dossier Éleveur :", eleveurs)
        if st.sidebar.button("🧪 Injecter Données Démo"):
            inject_demo_data(db, view_user)
            st.sidebar.success("Données injectées !")
    else: view_user = user

    menu = ["📊 Dashboard Élite", "🧬 Sélection & Agnelles", "🥛 Labo Biochimie", "🍲 Nutrition & Ration", "📦 Stocks", "📝 Registre"]
    choice = st.sidebar.radio("Navigation Master", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title(f"📊 Dashboard - {view_user}")
        df = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df.empty:
            df['ISG'] = df.apply(ScienceEngine.calculer_isg, axis=1)
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df))
            c2.metric("ISG Moyen", round(df['ISG'].mean(), 1))
            c3.metric("Poids Moyen", f"{round(df['poids'].mean(), 1)} kg")
            
            st.plotly_chart(px.scatter(df, x="poids", y="note_mamelle", color="race", trendline="ols", title="Corrélation Poids/Morphologie"))
            st.subheader("🏆 Top 10 Élite")
            st.dataframe(df.sort_values(by="ISG", ascending=False).head(10)[['identifiant_unique', 'race', 'type_animal', 'ISG']])
        else: st.info("Aucune donnée. Utilisez le module Registre.")

    # --- 2. GÉNETIQUE & AGNELLES ---
    elif choice == "🧬 Sélection & Agnelles":
        st.title("🧬 Sélection Prédictive")
        df_all = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        agnelles = df_all[df_all['type_animal'] == 'Agnelle']
        if not agnelles.empty:
            res = []
            for _, row in agnelles.iterrows():
                pot = ScienceEngine.predire_agnelle(row, df_all)
                res.append({"ID": row['identifiant_unique'], "Père": row['pere_id'], "Mère": row['mere_id'], "Index Potentiel": pot})
            st.table(pd.DataFrame(res).sort_values(by="Index Potentiel", ascending=False))
        else: st.warning("Pas d'agnelles avec parents identifiés.")

    # --- 3. LABO BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Analyses Biochimiques")
        
        with st.form("lait"):
            ids = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=? AND type_animal='Brebis Adulte'", (view_user,))
            target = st.selectbox("Brebis", ids['identifiant_unique']) if not ids.empty else None
            c1, c2, c3 = st.columns(3)
            qte = c1.number_input("Quantité (L)", 0.0, 10.0, 1.5)
            tb = c2.number_input("TB (g/L)", 20.0, 90.0, 45.0)
            tp = c3.number_input("TP (g/L)", 20.0, 90.0, 38.0)
            if st.form_submit_button("Enregistrer Labo"):
                esd = ScienceEngine.calculer_esd(tb, tp)
                db.execute_query("INSERT INTO lait_biochimie (brebis_id, qte_lait, tb, tp, esd, date_controle, owner_id) VALUES (?,?,?,?,?,?,?)",
                                (target, qte, tb, tp, esd, date.today(), view_user))
                st.success(f"Analyse enregistrée. ESD: {esd}")

    # --- 4. NUTRITION & RATION ---
    elif choice == "🍲 Nutrition & Ration":
        st.title("🍲 Rationnement et Rentabilité")
        
        px_lait = st.sidebar.number_input("Prix Lait (DA/L)", 50, 200, 100)
        
        st.subheader("⚙️ Configurer Prix Aliments (DA/100kg)")
        prix = {}
        cols = st.columns(3)
        for i, al in enumerate(TABLE_VALEURS.keys()):
            prix[al] = cols[i%3].number_input(f"{al}", 500, 15000, 4500)
        
        st.divider()
        choix = st.multiselect("Aliments du mélange", list(TABLE_VALEURS.keys()), default=["Orge", "Foin de luzerne"])
        qtes = {a: st.number_input(f"Kg de {a}/jour", 0.0, 5.0, 0.5) for a in choix}
        
        if choix:
            ufl = sum(qtes[a] * TABLE_VALEURS[a]['UFL'] for a in choix)
            pdi = sum(qtes[a] * TABLE_VALEURS[a]['PDI'] for a in choix)
            cout = sum(qtes[a] * (prix[a]/100) for a in choix)
            st.metric("Coût Journalier / Tête", f"{round(cout, 2)} DA")
            st.info(f"Apport : {round(ufl, 2)} UFL | {round(pdi, 0)} g PDI")

    # --- 5. STOCKS ---
    elif choice == "📦 Stocks":
        st.title("📦 Inventaire & Autonomie")
        
        al = st.selectbox("Aliment", list(TABLE_VALEURS.keys()))
        q = st.number_input("Quantité en stock (Quintaux)", 0.0, 1000.0)
        if st.button("Mettre à jour Stock"):
            db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (view_user, al, q))
        
        st.dataframe(db.fetch_all_as_df("SELECT aliment, quantite_q FROM stocks WHERE owner_id=?", (view_user,)))

    # --- 6. REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Registre du Troupeau")
        with st.form("add"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle")
            cat = c1.selectbox("Type", ["Brebis Adulte", "Agnelle", "Bélier"])
            rac = c1.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi"])
            pds = c2.number_input("Poids (kg)", 10.0, 150.0, 60.0)
            mam = c2.slider("Note Mamelle", 1.0, 10.0, 5.0)
            p_id = c2.text_input("ID Père")
            m_id = c2.text_input("ID Mère")
            if st.form_submit_button("Enregistrer"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, note_mamelle, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?,?)",
                                (uid, view_user, rac, cat, pds, mam, p_id, m_id, date.today()))
                st.success("Animal ajouté.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

# ============================================================================
# 4. FONCTION DE SIMULATION
# ============================================================================

def inject_demo_data(db, user):
    data = [
        ("D_B01", "Brebis Adulte", "Ouled Djellal", 72, 8.5, "P_A", "M_A"),
        ("D_B02", "Brebis Adulte", "Lacaune", 55, 9.5, "P_B", "M_B"),
        ("D_A01", "Agnelle", "Lacaune", 24, 5.0, "D_B03", "D_B02"),
        ("D_B03", "Bélier", "Lacaune", 90, 7.0, None, None)
    ]
    for uid, cat, rac, pds, mam, p, m in data:
        db.execute_query("INSERT OR IGNORE INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, note_mamelle, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?,?)",
                        (uid, user, rac, cat, pds, mam, p, m, date.today()))
    db.execute_query("INSERT OR IGNORE INTO lait_biochimie (brebis_id, qte_lait, tb, tp, esd, date_controle, owner_id) VALUES (?,?,?,?,?,?,?)",
                    ("D_B02", 3.5, 62, 48, 15.5, date.today(), user))

if __name__ == "__main__":
    main()
