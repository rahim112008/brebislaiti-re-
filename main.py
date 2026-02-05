"""
EXPERT OVIN DZ PRO - VERSION ULTIMATE MASTER V16.FIX
----------------------------------------------------
FIX: Nettoyage complet de l'indentation et migration DB
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date

# ============================================================================
# 1. DATABASE ENGINE
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_v16_final.db"):
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
        try:
            return pd.read_sql_query(query, self.conn, params=params)
        except:
            return pd.DataFrame()

def init_database(db: DatabaseManager):
    tables = [
        "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT)",
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE NOT NULL,
            owner_id TEXT, race TEXT, type_animal TEXT, poids REAL, note_mamelle REAL, 
            pere_id TEXT, mere_id TEXT, created_at DATE
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
# 2. MOTEUR SCIENTIFIQUE
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
        m = row['note_mamelle'] if row.get('note_mamelle') else 5.0
        p = row['poids'] if row.get('poids') else 50.0
        return round((m * 6) + (p * 0.4), 1)

    @staticmethod
    def predire_agnelle(agnelle_row, df_troupeau):
        if df_troupeau.empty: return 50.0
        pere = df_troupeau[df_troupeau['identifiant_unique'] == agnelle_row['pere_id']]
        mere = df_troupeau[df_troupeau['identifiant_unique'] == agnelle_row['mere_id']]
        score_p = 50 if pere.empty else ScienceEngine.calculer_isg(pere.iloc[0])
        score_m = 50 if mere.empty else ScienceEngine.calculer_isg(mere.iloc[0])
        return round((score_p * 0.5) + (score_m * 0.5), 1)

# ============================================================================
# 3. INTERFACE PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin Ultimate", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    
    db = st.session_state.db

    if 'auth' not in st.session_state:
        st.title("🛡️ Station Master Ovin DZ")
        u = st.text_input("Username")
        p = st.text_input("Password", type="password")
        if st.button("Connexion"):
            res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role
    u_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
    view_user = st.sidebar.selectbox("📁 Éleveur", u_list) if (role == "Expert" and u_list) else user

    if role == "Expert":
        if st.sidebar.button("🧪 ACTION : Injecter Données Démo"):
            inject_demo_data(db, view_user)
            st.rerun()

    menu = ["📊 Dashboard", "🍲 Nutrition & Ration", "📦 Stocks & Autonomie", "🧬 Génétique Agnelles", "🥛 Labo Biochimie", "📝 Registre"]
    choice = st.sidebar.radio("Modules Master", menu)

    if choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - {view_user}")
        df = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df.empty:
            df['ISG'] = df.apply(ScienceEngine.calculer_isg, axis=1)
            c1, c2, c3 = st.columns(3)
            c1.metric("Têtes", len(df))
            c2.metric("Score Moyen", f"{round(df['ISG'].mean(), 1)}/100")
            c3.metric("Poids Moyen", f"{round(df['poids'].mean(), 1)} kg")
            st.plotly_chart(px.scatter(df, x="poids", y="note_mamelle", color="race", title="Analyse Troupeau"))
        else:
            st.info("Troupeau vide. Cliquez sur 'Injecter Données' ou allez dans 'Registre'.")

    elif choice == "🍲 Nutrition & Ration":
        st.title("🍲 Nutrition & Rentabilité")
        
        cols = st.columns(3)
        prix = {al: cols[i%3].number_input(f"{al} (DA/100kg)", 50, 15000, 4500) for i, al in enumerate(TABLE_VALEURS.keys())}
        st.divider()
        choix = st.multiselect("Composer le mélange", list(TABLE_VALEURS.keys()), default=["Orge"])
        qtes = {a: st.number_input(f"Kg de {a}/animal", 0.0, 5.0, 0.5) for a in choix}
        if choix:
            cout = sum(qtes[a] * (prix[a]/100) for a in choix)
            st.metric("Coût Journalier", f"{round(cout, 2)} DA")

    elif choice == "📦 Stocks & Autonomie":
        st.title("📦 Gestion des Stocks")
        
        al = st.selectbox("Aliment", list(TABLE_VALEURS.keys()))
        q = st.number_input("Quantité (Quintaux)", 0.0, 1000.0)
        if st.button("Mettre à jour"):
            db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (view_user, al, q))
        st.dataframe(db.fetch_all_as_df("SELECT aliment, quantite_q FROM stocks WHERE owner_id=?", (view_user,)))

    elif choice == "🧬 Génétique Agnelles":
        st.title("🧬 Sélection Agnelles")
        df_all = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df_all.empty and 'type_animal' in df_all.columns:
            agnelles = df_all[df_all['type_animal'] == 'Agnelle']
            if not agnelles.empty:
                res = [{"ID": r['identifiant_unique'], "Potentiel": ScienceEngine.predire_agnelle(r, df_all)} for _, r in agnelles.iterrows()]
                st.dataframe(pd.DataFrame(res))
            else:
                st.info("Aucune agnelle enregistrée.")
        else:
            st.warning("Données insuffisantes.")

    elif choice == "📝 Registre":
        st.title("📝 Registre")
        with st.form("reg"):
            uid = st.text_input("ID")
            cat = st.selectbox("Type", ["Brebis Adulte", "Agnelle", "Bélier"])
            rac = st.selectbox("Race", ["Ouled Djellal", "Lacaune"])
            pds = st.number_input("Poids", 10, 150, 60)
            p_id = st.text_input("ID Père")
            m_id = st.text_input("ID Mère")
            if st.form_submit_button("Inscrire"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?)",
                                (uid, view_user, rac, cat, pds, p_id, m_id, date.today()))
                st.success("OK")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

def inject_demo_data(db, user):
    db.execute_query("INSERT OR IGNORE INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, created_at) VALUES (?,?,?,?,?,?)",
                    ("DEMO_01", user, "Lacaune", "Brebis Adulte", 65, date.today()))
    db.execute_query("INSERT OR IGNORE INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?)",
                    ("AG_DEMO", user, "Lacaune", "Agnelle", 25, "B_ELITE", "DEMO_01", date.today()))
    db.execute_query("INSERT OR REPLACE INTO stocks VALUES (?,?,100,4500)", (user, "Orge"))

if __name__ == "__main__":
    main()
