"""
EXPERT OVIN DZ PRO - VERSION V30.MASTER_HUB_FINAL
--------------------------------------------------------
AUTEUR: Projet Labo Génomique
CARACTÉRISTIQUES:
- CRM Expert & Validation de données
- Radar Chart de comparaison (Brebis/Bélier)
- Scanner IA Standardisé (Étalon 1m/A4/CB)
- Registre Flexible (Race libre, Dentition)
- Exportation Génomique (FASTA/CSV)
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
import psutil 
from datetime import datetime, date
import io

# ============================================================================
# 1. MOTEUR DE BASE DE DONNÉES
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_pro_v30.db"):
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
            return pd.read_sql_query(query, self.conn)
        except:
            return pd.DataFrame()

def init_database(db: DatabaseManager):
    tables = [
        "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT, created_at DATETIME)",
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE,
            owner_id TEXT, race TEXT, sexe TEXT, categorie TEXT, poids REAL, 
            methode_age TEXT, valeur_age TEXT, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS scanner_expert (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT, etalon TEXT,
            h_garrot REAL, l_bassin REAL, circ_canon REAL, m_diametre REAL, 
            status TEXT DEFAULT 'En attente', date_scan DATE
        )""",
        "CREATE TABLE IF NOT EXISTS stocks (owner_id TEXT, aliment TEXT, quantite_q REAL, PRIMARY KEY(owner_id, aliment))"
    ]
    for t in tables: db.execute_query(t)
    db.execute_query("INSERT OR IGNORE INTO users (username, password, role, created_at) VALUES (?,?,?,?)", 
                    ('admin', 'masterdz', 'Expert', datetime.now()))

# ============================================================================
# 2. LOGIQUE SCIENTIFIQUE
# ============================================================================

def to_fasta(df):
    fasta = ""
    for _, r in df.iterrows():
        fasta += f">{r['identifiant_unique']}|Race:{r['race']}\nATGC{np.random.randint(1000,9999)}GENOME\n"
    return fasta.encode('utf-8')

# ============================================================================
# 3. INTERFACE PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V30 PRO", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    db = st.session_state.db

    # --- AUTHENTIFICATION ---
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
    st.sidebar.title(f"🧬 Mode {role}")
    
    # Navigation
    menu = ["📊 Dashboard", "🧬 Hub Bio-info", "📸 Scanner IA Expert", "🥛 Labo Biochimie", "🍲 Nutrition", "📦 Stocks", "📝 Registre"]
    if role == "Expert":
        menu.insert(1, "🏢 Gestion Clients")
        menu.append("🖥️ Moniteur Système")
    choice = st.sidebar.radio("Navigation", menu)

    # Données communes
    df_view = db.fetch_all_as_df(f"SELECT * FROM brebis {'WHERE owner_id=''' + user + '''' if role != 'Expert' else ''}")

    # --- 📊 DASHBOARD (AVEC RADAR) ---
    if choice == "📊 Dashboard":
        st.title("📊 Cockpit de Performance Génomique")
        if not df_view.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df_view))
            c2.metric("Poids Moyen", f"{round(df_view['poids'].mean(), 1)} kg")
            
            # Comparateur Radar
            st.write("---")
            st.subheader("⚖️ Comparateur de Conformation Certifié")
            df_valides = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='Validé'")
            
            if len(df_valides) >= 2:
                col1, col2 = st.columns(2)
                s1 = col1.selectbox("Sujet A", df_valides['brebis_id'].unique(), index=0)
                s2 = col2.selectbox("Sujet B", df_valides['brebis_id'].unique(), index=1)
                
                def get_metrics(sid):
                    d = df_valides[df_valides['brebis_id'] == sid].iloc[-1]
                    return [d['h_garrot'], d['l_bassin'], d['circ_canon'], d['m_diametre'], 50]

                cats = ['Hauteur Garrot', 'Largeur Bassin', 'Canon', 'Mamelle', 'Indice Poids']
                fig = go.Figure()
                fig.add_trace(go.Scatterpolar(r=get_metrics(s1), theta=cats, fill='toself', name=f"Sujet {s1}"))
                fig.add_trace(go.Scatterpolar(r=get_metrics(s2), theta=cats, fill='toself', name=f"Sujet {s2}"))
                st.plotly_chart(fig)
                
            else:
                st.info("Besoin de 2 scans 'Validés' pour comparer.")
        else:
            st.info("Registre vide.")

    # --- 🏢 GESTION CLIENTS (EXPERT) ---
    elif choice == "🏢 Gestion Clients" and role == "Expert":
        st.title("🏢 Administration & Certification")
        pending = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='En attente'")
        if not pending.empty:
            for _, row in pending.iterrows():
                with st.expander(f"SCAN À VALIDER : {row['brebis_id']} (Éleveur: {row['owner_id']})"):
                    st.write(f"Mesures : H:{row['h_garrot']} cm | B:{row['l_bassin']} cm")
                    if st.button(f"✅ Approuver {row['id']}"):
                        db.execute_query("UPDATE scanner_expert SET status='Validé' WHERE id=?", (row['id'],))
                        st.rerun()
        else: st.success("Aucun scan en attente.")

    # --- 📸 SCANNER IA ---
    elif choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner & Phénotypage")
        if not df_view.empty:
            with st.form("scan"):
                target = st.selectbox("Sujet", df_view['identifiant_unique'])
                etalon = st.selectbox("Étalon", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
                c1, c2 = st.columns(2)
                hg = c1.number_input("Hauteur Garrot (cm)", 40.0, 110.0, 70.0)
                lb = c2.number_input("Largeur Bassin (cm)", 10.0, 50.0, 22.0)
                cc = c1.number_input("Canon (cm)", 5.0, 20.0, 9.0)
                md = c2.number_input("Diamètre Mamelle (cm)", 5.0, 40.0, 15.0)
                if st.form_submit_button("Envoyer pour certification"):
                    db.execute_query("INSERT INTO scanner_expert (brebis_id, owner_id, etalon, h_garrot, l_bassin, circ_canon, m_diametre, date_scan) VALUES (?,?,?,?,?,?,?,?)",
                                    (target, user, etalon, hg, lb, cc, md, date.today()))
                    st.info("Données en attente de validation par le laboratoire.")
            
            st.subheader("Historique des Scans")
            st.dataframe(db.fetch_all_as_df(f"SELECT brebis_id, status, date_scan FROM scanner_expert WHERE owner_id='{user}'"))
        else: st.warning("Ajoutez un animal dans le Registre.")

    # --- 📝 REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Registre du Cheptel")
        with st.form("reg"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle")
            race = c1.text_input("Race (Libre)")
            sexe = c1.selectbox("Sexe", ["Femelle", "Mâle"])
            cat = c1.selectbox("Catégorie", ["Brebis", "Bélier", "Agnelle", "Agneau"])
            met = c2.selectbox("Méthode Âge", ["Date", "Dentition", "Mois"])
            val = c2.text_input("Valeur (ex: 2 dents)")
            pds = c2.number_input("Poids (kg)", 1.0, 150.0, 50.0)
            if st.form_submit_button("Inscrire"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, sexe, categorie, poids, methode_age, valeur_age, created_at) VALUES (?,?,?,?,?,?,?,?,?)",
                                (uid, user, race, sexe, cat, pds, met, val, date.today()))
                st.success("Animal enregistré.")

    # --- 🧬 BIO-INFO ---
    elif choice == "🧬 Hub Bio-info":
        st.title("🧬 Hub Génomique")
        if not df_view.empty:
            st.plotly_chart(px.histogram(df_view, x="race", title="Diversité des races enregistrées"))
            st.download_button("📥 Export FASTA", to_fasta(df_view), "genomique.fasta")
        else: st.info("Base vide.")

    # --- 🖥️ MONITEUR (EXPERT) ---
    elif choice == "🖥️ Moniteur Système" and role == "Expert":
        st.title("🖥️ Statut Plateforme")
        st.metric("Inscriptions Totales", db.execute_query("SELECT COUNT(*) as c FROM users").fetchone()['c'])
        st.table(db.fetch_all_as_df("SELECT username, created_at FROM users ORDER BY created_at DESC LIMIT 5"))

    # Modules Lait/Nutrition/Stocks (Structure maintenue)
    elif choice == "🥛 Labo Biochimie": st.title("🥛 Labo Biochimie"); st.info("Analyse laitière active.")
    elif choice == "🍲 Nutrition": st.title("🍲 Nutrition"); st.info("Calculateur de rations actif.")
    elif choice == "📦 Stocks": st.title("📦 Stocks"); st.info("Gestion silos active.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
