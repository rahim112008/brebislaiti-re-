"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE REPARÉE 2026
Logiciel Master : Génétique, Santé, Finance et Rapports PDF
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. BASE DE DONNÉES ET STRUCTURE
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_final.db"):
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
            owner_id TEXT, race TEXT, sexe TEXT, note_mamelle REAL DEFAULT 5.0, 
            poids REAL DEFAULT 60.0, pere_id TEXT, mere_id TEXT, created_at DATE
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
# 2. LOGIQUE MÉTIER ET INTERFACE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro", layout="wide")
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

    # --- BARRE LATERALE ---
    st.sidebar.title(f"👤 {user}")
    if role == "Expert":
        u_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")
        view_user = st.sidebar.selectbox("📂 Dossier Éleveur :", u_list['username'] if not u_list.empty else [user])
    else: view_user = user

    menu = ["📊 Dashboard", "🩺 Calendrier Sanitaire", "🧬 Génétique & Accouplement", "📝 Registre", "📈 Rapport Master"]
    choice = st.sidebar.radio("Modules", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - {view_user}")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df_b.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df_b))
            c2.metric("Note Mamelle Moy.", round(df_b['note_mamelle'].mean(), 1))
            c3.metric("Poids Moyen", f"{round(df_b['poids'].mean(), 1)} kg")
            st.plotly_chart(px.histogram(df_b, x="poids", color="race", title="Analyse de la structure du troupeau"))
        else: st.info("Le troupeau est vide. Enregistrez des animaux dans le module 'Registre'.")

    # --- 2. CALENDRIER SANITAIRE ---
    elif choice == "🩺 Calendrier Sanitaire":
        st.title("🩺 Suivi Sanitaire & Prophylaxie")
        with st.expander("➕ Ajouter un soin/vaccin"):
            with st.form("sante"):
                ids = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=?", (view_user,))
                target = st.selectbox("Brebis", ids['identifiant_unique']) if not ids.empty else None
                acte = st.selectbox("Intervention", ["Vaccin Entero", "Vaccin PPR", "Déparasitage", "Traitement Mammite"])
                date_j = st.date_input("Date de l'acte")
                if st.form_submit_button("Enregistrer"):
                    rappel = date_j + timedelta(days=180)
                    db.execute_query("INSERT INTO sante (brebis_id, acte, date_acte, date_rappel, statut, owner_id) VALUES (?,?,?,?,?,?)",
                                    (target, acte, date_j, rappel, "En attente", view_user))
                    st.success("Soin enregistré !")
        
        df_s = db.fetch_all_as_df("SELECT * FROM sante WHERE owner_id=?", (view_user,))
        st.dataframe(df_s, use_container_width=True)

    # --- 3. GÉNÉTIQUE & ACCOUPLEMENT (REPARÉ) ---
    elif choice == "🧬 Génétique & Accouplement":
        st.title("🧬 Intelligence de Sélection")
        df_all = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        
        if not df_all.empty:
            femelles = df_all[df_all['sexe'] == 'Femelle']
            males = df_all[df_all['sexe'] == 'Mâle']
            
            st.subheader("Suggestions d'accouplement pour améliorer le troupeau")
            if not males.empty and not femelles.empty:
                results = []
                meilleur_belier = males.sort_values(by='note_mamelle', ascending=False).iloc[0]
                for _, f in femelles.iterrows():
                    action = "Amélioration Mamelle" if f['note_mamelle'] < 6 else "Maintien"
                    results.append({"Femelle": f['identifiant_unique'], "Note": f['note_mamelle'], "Bélier suggéré": meilleur_belier['identifiant_unique'], "Objectif": action})
                st.table(pd.DataFrame(results))
            else:
                st.warning("Pour simuler un accouplement, vous devez avoir enregistré au moins une Femelle et un Mâle.")
        else:
            st.info("Aucune donnée génétique disponible.")

    # --- 4. REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Registre des Animaux")
        with st.form("add_ov"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("Identifiant (Boucle)")
            sex = c1.selectbox("Sexe", ["Femelle", "Mâle"])
            rac = c2.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra"])
            mam = c2.slider("Note Mamelle (1-10)", 1.0, 10.0, 5.0)
            if st.form_submit_button("Inscrire l'animal"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, sexe, note_mamelle, created_at) VALUES (?,?,?,?,?,?)",
                                (uid, view_user, rac, sex, mam, date.today()))
                st.success("Animal ajouté au registre.")

    # --- 5. RAPPORT MASTER (REPARÉ) ---
    elif choice == "📈 Rapport Master":
        st.title("📈 Génération du Rapport d'Expertise")
        st.write("Ce module compile toutes les données pour créer une synthèse professionnelle.")
        
        df_rep = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df_rep.empty:
            st.subheader(f"Synthèse pour l'éleveur : {view_user}")
            st.write(f"- Nombre d'animaux analysés : {len(df_rep)}")
            st.write(f"- Qualité génétique moyenne : {round(df_rep['note_mamelle'].mean(), 2)} / 10")
            
            if st.button("📄 Télécharger le Rapport Complet"):
                csv = df_rep.to_csv(index=False).encode('utf-16')
                st.download_button("Confirmer le téléchargement (CSV/Excel)", csv, f"Rapport_{view_user}.csv", "text/csv")
        else:
            st.error("Impossible de générer un rapport : le troupeau est vide.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
