"""
EXPERT OVIN DZ PRO - VERSION ULTIMATE MASTER V15.FINAL
----------------------------------------------------
Modules inclus : Génétique, Biochimie, Nutrition IA, Finance, Stocks, Outils Expert
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
    def __init__(self, db_path: str = "data/ovin_master_ultimate_v15.db"):
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
# 2. CALCULS SCIENTIFIQUES (NUTRITION & GÉNÉTIQUE)
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

# ============================================================================
# 3. INTERFACE STREAMLIT
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin Ultimate", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    
    db = st.session_state.db

    if 'auth' not in st.session_state:
        st.title("🛡️ Connexion Expert Ovin")
        with st.form("login"):
            u = st.text_input("Username")
            p = st.text_input("Password", type="password")
            if st.form_submit_button("Se connecter"):
                res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
                if res:
                    st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                    st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role
    
    # --- Barre Latérale ---
    st.sidebar.title(f"📍 {user}")
    u_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
    view_user = st.sidebar.selectbox("📁 Éleveur", u_list) if role == "Expert" else user

    if role == "Expert":
        if st.sidebar.button("🧪 Injecter Données Démo"):
            inject_demo_data(db, view_user)
            st.sidebar.success("Données de test injectées !")

    # NOUVEAU MENU COMPLET
    menu = ["📊 Dashboard", "🍲 Nutrition & Ration", "📦 Stocks & Autonomie", "🧬 Génétique Agnelles", "🥛 Labo Biochimie", "📝 Registre"]
    choice = st.sidebar.radio("Modules", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - {view_user}")
        df = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df.empty:
            df['ISG'] = df.apply(ScienceEngine.calculer_isg, axis=1)
            c1, c2, c3 = st.columns(3)
            c1.metric("Têtes", len(df))
            c2.metric("Score Moyen", f"{round(df['ISG'].mean(), 1)}/100")
            c3.metric("Poids Moyen", f"{round(df['poids'].mean(), 1)} kg")
            st.plotly_chart(px.scatter(df, x="poids", y="note_mamelle", color="race", title="Analyse du Troupeau"))
        else:
            st.info("Le troupeau est vide. Enregistrez des animaux dans le module 'Registre' ou utilisez le bouton démo.")

    # --- 2. NUTRITION & RATION ---
    elif choice == "🍲 Nutrition & Ration":
        st.title("🍲 Nutrition de Précision")
        st.subheader("Configuration des prix (DA/Quintal)")
        cols = st.columns(3)
        prix = {al: cols[i%3].number_input(f"{al}", 500, 15000, 4500) for i, al in enumerate(TABLE_VALEURS.keys())}
        
        st.divider()
        choix = st.multiselect("Composer le mélange", list(TABLE_VALEURS.keys()), default=["Orge", "Foin de luzerne"])
        qtes = {a: st.number_input(f"Kg de {a}/jour/animal", 0.0, 5.0, 0.5) for a in choix}
        
        if choix:
            ufl = sum(qtes[a] * TABLE_VALEURS[a]['UFL'] for a in choix)
            cout = sum(qtes[a] * (prix[a]/100) for a in choix)
            st.metric("Coût Journalier par Animal", f"{round(cout, 2)} DA")
            st.info(f"Énergie totale du mélange : {round(ufl, 2)} UFL")

    # --- 3. STOCKS ---
    elif choice == "📦 Stocks & Autonomie":
        st.title("📦 Gestion des Stocks")
        with st.form("stk"):
            al = st.selectbox("Aliment", list(TABLE_VALEURS.keys()))
            q = st.number_input("Quantité en stock (Quintaux)", 0.0, 1000.0)
            if st.form_submit_button("Mettre à jour"):
                db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (view_user, al, q))
        
        df_s = db.fetch_all_as_df("SELECT aliment, quantite_q FROM stocks WHERE owner_id=?", (view_user,))
        st.table(df_s)

    # --- 4. GÉNÉTIQUE ---
    elif choice == "🧬 Génétique Agnelles":
        st.title("🧬 Prédiction Génétique des Agnelles")
        df_all = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        agnelles = df_all[df_all['type_animal'] == 'Agnelle']
        if not agnelles.empty:
            res = [{"ID": r['identifiant_unique'], "Potentiel": ScienceEngine.predire_agnelle(r, df_all)} for _, r in agnelles.iterrows()]
            st.dataframe(pd.DataFrame(res).sort_values(by="Potentiel", ascending=False))
        else: st.warning("Pas d'agnelles avec pedigree complet (père/mère).")

    # --- 5. LABO BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Labo : Qualité du Lait")
        # Logique de saisie TB/TP/ESD simplifiée
        st.info("Module d'analyse des composants du lait.")

    # --- 6. REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Inscription des Sujets")
        with st.form("reg"):
            uid = st.text_input("ID Boucle")
            cat = st.selectbox("Catégorie", ["Brebis Adulte", "Agnelle", "Bélier"])
            rac = st.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi"])
            pds = st.number_input("Poids (kg)", 10, 150, 60)
            mam = st.slider("Note Mamelle", 1.0, 10.0, 5.0)
            p_id = st.text_input("ID Père (pour les agnelles)")
            m_id = st.text_input("ID Mère (pour les agnelles)")
            if st.form_submit_button("Inscrire"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, note_mamelle, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?,?)",
                                (uid, view_user, rac, cat, pds, mam, p_id, m_id, date.today()))
                st.success("Animal enregistré.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

def inject_demo_data(db, user):
    db.execute_query("INSERT OR IGNORE INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, note_mamelle, created_at) VALUES (?,?,?,?,?,?,?)",
                    ("D_B01", user, "Ouled Djellal", "Brebis Adulte", 75, 8.0, date.today()))
    db.execute_query("INSERT OR IGNORE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (user, "Orge", 50))

if __name__ == "__main__":
    main()
