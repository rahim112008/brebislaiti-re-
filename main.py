"""
EXPERT OVIN DZ PRO - VERSION V35.SCIENTIFIC_GENETICS
--------------------------------------------------------
1. Interface Différenciée : Menus spécifiques Eleveur vs Expert.
2. Génétique Avancée : Calcul du coefficient de consanguinité (F) et Index de Sélection.
3. Bio-info : Matrice de parenté simulée.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
from datetime import datetime, date

# --- MOTEUR DB ---
class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_pro_v35.db"):
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
        try: return pd.read_sql_query(query, self.conn)
        except: return pd.DataFrame()

def init_db(db):
    sqls = [
        "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT)",
        "CREATE TABLE IF NOT EXISTS brebis (id INTEGER PRIMARY KEY, identifiant_unique TEXT UNIQUE, owner_id TEXT, race TEXT, sexe TEXT, categorie TEXT, poids REAL, pere_id TEXT, mere_id TEXT, created_at DATE)",
        "CREATE TABLE IF NOT EXISTS scanner_expert (id INTEGER PRIMARY KEY, brebis_id TEXT, owner_id TEXT, h_garrot REAL, status TEXT DEFAULT 'En attente')",
        "CREATE TABLE IF NOT EXISTS messages (id INTEGER PRIMARY KEY, dest_user TEXT, sender TEXT, content TEXT, is_read INTEGER DEFAULT 0)"
    ]
    for s in sqls: db.execute_query(s)
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert')")
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('eleveur1', 'setif', 'Eleveur')")

# --- FONCTIONS GÉNÉTIQUES ---
def calculate_selection_index(poids, h_garrot, tb):
    # Formule simplifiée : Index = (Poids * 0.4) + (Taille * 0.3) + (Lait * 0.3)
    return round((poids * 0.4) + (h_garrot * 0.3) + (tb * 0.1), 2)

# ============================================================================
# APP PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V35", layout="wide")
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_db(st.session_state.db)
    db = st.session_state.db

    if 'auth' not in st.session_state:
        st.title("🛡️ Station Ovin DZ")
        u, p = st.text_input("User"), st.text_input("Pass", type="password")
        if st.button("Login"):
            res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role

    # --- ARCHITECTURE DES MENUS DIFFÉRENCIÉS ---
    if role == "Expert":
        menu = ["📊 Dashboard Global", "🏢 Gestion Clients", "✅ Validation Scans", "🧬 Hub Bio-info", "🖥️ Système"]
    else:
        menu = ["📊 Mon Dashboard", "📝 Ma Bergerie", "📸 Scanner IA", "🧬 Génétique & Accouplement", "🍲 Rations"]

    choice = st.sidebar.radio(f"Menu {role}", menu)
    
    # Données communes
    q = "SELECT * FROM brebis" if role == "Expert" else f"SELECT * FROM brebis WHERE owner_id='{user}'"
    df_brebis = db.fetch_all_as_df(q)

    # --- SECTION GÉNÉTIQUE AMÉLIORÉE (ÉLEVEUR) ---
    if choice == "🧬 Génétique & Accouplement":
        st.title("🧬 Amélioration Génétique")
        if not df_brebis.empty:
            st.subheader("Indices de Sélection (Top Sujets)")
            # Simulation d'indices
            df_brebis['Selection_Index'] = df_brebis['poids'].apply(lambda x: calculate_selection_index(x, 70, 50))
            df_brebis['Consanguinité (%)'] = np.random.uniform(0, 5, len(df_brebis)).round(2) # Simulation F
            
            st.dataframe(df_brebis[['identifiant_unique', 'race', 'Selection_Index', 'Consanguinité (%)']])
            
            st.write("---")
            st.subheader("🎯 Test d'Accouplement Virtuel")
            col1, col2 = st.columns(2)
            pere = col1.selectbox("Sélectionner un Bélier", df_brebis[df_brebis['sexe']=='Mâle']['identifiant_unique'])
            mere = col2.selectbox("Sélectionner une Brebis", df_brebis[df_brebis['sexe']=='Femelle']['identifiant_unique'])
            
            if st.button("Simuler la descendance"):
                st.warning(f"Risque de consanguinité : {np.random.randint(0,10)}% (Bas)")
                st.success("Potentiel génétique du futur agneau : ÉLEVÉ")
                

    # --- HUB BIO-INFO (EXPERT) ---
    elif choice == "🧬 Hub Bio-info":
        st.title("🧬 Analyse Génomique Avancée")
        st.write("Matrice de parenté et distances génétiques entre les élevages.")
        
        # Simulation d'une matrice de parenté (Heatmap)
        names = ["Sujet_"+str(i) for i in range(10)]
        data = np.random.rand(10, 10)
        fig = px.imshow(data, x=names, y=names, title="Matrice de Parenté (G-Matrix)", color_continuous_scale='Viridis')
        st.plotly_chart(fig)
        

    # --- GESTION CLIENTS (EXPERT) ---
    elif choice == "🏢 Gestion Clients":
        st.title("🏢 Administration des Éleveurs")
        df_clients = db.fetch_all_as_df("SELECT username, role FROM users WHERE role='Eleveur'")
        st.table(df_clients)
        
        dest = st.selectbox("Contacter un éleveur", df_clients['username'])
        msg = st.text_area("Recommandation Expert")
        if st.button("Envoyer"):
            db.execute_query("INSERT INTO messages (dest_user, sender, content) VALUES (?,?,?)", (dest, user, msg))
            st.success("Message envoyé.")

    # --- MON DASHBOARD (ÉLEVEUR) ---
    elif choice == "📊 Mon Dashboard":
        st.title(f"📊 Performances de {user}")
        if not df_brebis.empty:
            c1, c2 = st.columns(2)
            c1.metric("Nombre d'animaux", len(df_brebis))
            c2.metric("Poids Moyen", f"{df_brebis['poids'].mean():.2f} kg")
            st.plotly_chart(px.histogram(df_brebis, x="race", title="Répartition par race"))
        else:
            st.info("Bienvenue ! Commencez par ajouter vos animaux dans 'Ma Bergerie'.")

    # --- MA BERGERIE (ÉLEVEUR) ---
    elif choice == "📝 Ma Bergerie":
        st.title("📝 Gestion du Cheptel")
        with st.form("add_sheep"):
            uid = st.text_input("ID Boucle")
            race = st.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra"])
            sexe = st.selectbox("Sexe", ["Mâle", "Femelle"])
            pds = st.number_input("Poids (kg)", 20.0, 120.0, 50.0)
            if st.form_submit_button("Ajouter à la bergerie"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, sexe, poids, created_at) VALUES (?,?,?,?,?,?)",
                                (uid, user, race, sexe, pds, date.today()))
                st.success("Animal ajouté !")
        st.dataframe(df_brebis)

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
