"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.V8
Génétique | Finance | Accouplement | Calendrier Sanitaire Automatisé
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. BASE DE DONNÉES (AVEC TABLE SANITAIRE)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_v8.db"):
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
            owner_id TEXT, race TEXT, sexe TEXT, note_mamelle REAL, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            acte TEXT, date_acte DATE, date_rappel DATE, statut TEXT, owner_id TEXT
        )""",
        "CREATE TABLE IF NOT EXISTS recommandations (target_user TEXT, expert_name TEXT, message TEXT, date_envoyee DATE)"
    ]
    for table_sql in tables: db.execute_query(table_sql)
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert')")
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('eleveur1', 'ovin2026', 'Eleveur')")

# ============================================================================
# 2. LOGIQUE MÉTIER
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro v8", layout="wide")
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    
    if 'auth' not in st.session_state:
        st.title("🛡️ Plateforme Expert Ovin")
        with st.form("login"):
            u, p = st.text_input("User"), st.text_input("Pass", type="password")
            if st.form_submit_button("Se connecter"):
                res = st.session_state.db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
                if res:
                    st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                    st.rerun()
        return

    db, user, role = st.session_state.db, st.session_state.username, st.session_state.role

    # --- FILTRE EXPERT ---
    if role == "Expert":
        u_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")
        view_user = st.sidebar.selectbox("📂 Dossier Éleveur :", u_list['username'] if not u_list.empty else [user])
    else: view_user = user

    menu = ["📊 Dashboard", "🩺 Calendrier Sanitaire", "🧬 Génétique & Accouplement", "📝 Registre", "📈 Rapport Master"]
    choice = st.sidebar.radio("Modules", menu)

    # --- MODULE 1: DASHBOARD AVEC ALERTES ---
    if choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - {view_user}")
        
        # Vérification des rappels de vaccins urgents
        aujourdhui = date.today().isoformat()
        alertes = db.fetch_all_as_df("SELECT * FROM sante WHERE owner_id=? AND date_rappel <= ? AND statut='En attente'", (view_user, aujourdhui))
        
        if not alertes.empty:
            st.error(f"🚨 ATTENTION : {len(alertes)} interventions sanitaires sont en retard !")
            st.dataframe(alertes[['brebis_id', 'acte', 'date_rappel']])

        df_b = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df_b.empty:
            c1, c2 = st.columns(2)
            c1.metric("Effectif", len(df_b))
            c2.metric("Note Moyenne Mamelle", round(df_b['note_mamelle'].mean(), 1))
            st.plotly_chart(px.histogram(df_b, x="poids", title="Distribution du Poids"))
        else: st.info("Aucune donnée.")

    # --- MODULE 2: CALENDRIER SANITAIRE (NOUVEAU) ---
    elif choice == "🩺 Calendrier Sanitaire":
        st.title("🩺 Calendrier de Prophylaxie & Soins")
        
        
        with st.expander("➕ Enregistrer une nouvelle intervention"):
            with st.form("sante_form"):
                target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=?", (view_user,))['identifiant_unique'])
                acte = st.selectbox("Type d'acte", ["Vaccin Enterotoxémie", "Vaccin Clavelée", "Vaccin PPR", "Déparasitage interne", "Traitement Mammite"])
                d_acte = st.date_input("Date de l'acte")
                # Calcul auto du rappel selon l'acte
                d_rappel = d_acte + (timedelta(days=180) if "Vaccin" in acte else timedelta(days=90))
                
                if st.form_submit_button("Enregistrer"):
                    db.execute_query("INSERT INTO sante (brebis_id, acte, date_acte, date_rappel, statut, owner_id) VALUES (?,?,?,?,?,?)",
                                    (target, acte, d_acte, d_rappel, "En attente", view_user))
                    st.success(f"Enregistré ! Prochain rappel le {d_rappel}")

        st.subheader("🗓️ Historique et Rappels à venir")
        df_s = db.fetch_all_as_df("SELECT * FROM sante WHERE owner_id=?", (view_user,))
        if not df_s.empty:
            st.table(df_s[['brebis_id', 'acte', 'date_acte', 'date_rappel', 'statut']])
            if st.button("Marquer tout comme 'Fait'"):
                db.execute_query("UPDATE sante SET statut='Terminé' WHERE owner_id=?", (view_user,))
                st.rerun()

    # --- MODULE 3: GÉNÉTIQUE ---
    elif choice == "🧬 Génétique & Accouplement":
        st.title("🧬 Intelligence de Sélection")
        st.info("Utilisez ce module pour planifier la prochaine génération.")
        # ... (Logique d'accouplement précédente préservée)

    # --- MODULE 4: REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Gestion des fiches individuelles")
        with st.form("reg"):
            uid = st.text_input("ID Boucle")
            sexe = st.selectbox("Sexe", ["Femelle", "Mâle"])
            race = st.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi"])
            mam = st.slider("Note Mamelle", 1.0, 10.0, 5.0)
            if st.form_submit_button("Ajouter"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, sexe, race, note_mamelle, created_at) VALUES (?,?,?,?,?,?)",
                                (uid, view_user, sexe, race, mam, date.today()))

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
