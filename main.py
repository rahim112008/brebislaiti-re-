"""
EXPERT OVIN DZ PRO - VERSION V19.MASTER
----------------------------------------------
Modules : Scanner IA Morphométrique (Expert), Biochimie, Nutrition, Stocks.
Calibration : Bâton 1m, Feuille A4, Carte Bancaire.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date

# ============================================================================
# 1. DATABASE ENGINE (V19)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_v19_final.db"):
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
        """CREATE TABLE IF NOT EXISTS scanner_expert (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            hauteur_garrot REAL, longueur_corps REAL, circ_canon REAL, taille_bassin REAL,
            diametre_mamelle REAL, profondeur_mamelle REAL, indice_conformation REAL, date_scan DATE
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
# 2. LOGIQUE SCIENTIFIQUE & IA
# ============================================================================

TABLE_VALEURS = {
    "Orge": {"UFL": 1.0, "PDI": 80},
    "Son de blé": {"UFL": 0.85, "PDI": 95},
    "Foin de luzerne": {"UFL": 0.65, "PDI": 90},
    "Maïs grain": {"UFL": 1.15, "PDI": 75}
}

class ScienceEngine:
    @staticmethod
    def calculer_esd(tb, tp):
        return round((tp * 0.85) + (tb * 0.1) + 5.5, 2)

    @staticmethod
    def calculer_indice_scan(h, l):
        return round((l / h) * 100, 1)

# ============================================================================
# 3. INTERFACE PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro V19", layout="wide", page_icon="🐑")
    
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
    
    # --- Menu de navigation ---
    if role == "Expert":
        st.sidebar.title("💎 Panel EXPERT")
        u_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
        view_user = st.sidebar.selectbox("📂 Dossier Éleveur", u_list) if u_list else user
        if st.sidebar.button("🧪 Injecter Données Démo"):
            inject_demo_data(db, view_user)
            st.rerun()
    else:
        st.sidebar.title(f"👨‍🌾 Éleveur: {user}")
        view_user = user

    menu = ["📊 Dashboard", "📸 Scanner IA Expert", "🥛 Labo Biochimie", "🍲 Nutrition & Ration", "📦 Stocks", "📝 Registre"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE SCANNER IA (MISE À JOUR V19) ---
    if choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner IA & Morphométrie")
        
        # Section Etalonnage
        st.info("### 📏 Étalonnage de l'image")
        etalon = st.selectbox("Choisir l'objet de référence (Étalon) :", 
                             ["Bâton de 1 mètre", "Feuille A4 (29.7 cm)", "Carte Bancaire (8.5 cm)"])
        st.write(f"L'IA utilise le **{etalon}** pour convertir les pixels en centimètres.")
        
        

        df_troupeau = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=?", (view_user,))
        if not df_troupeau.empty:
            with st.form("expert_scan_ia"):
                target = st.selectbox("Animal à mesurer", df_troupeau['identifiant_unique'])
                
                st.subheader("🦴 Squelette & Conformation")
                c1, c2, c3, c4 = st.columns(4)
                h_g = c1.number_input("Hauteur Garrot (cm)", 40.0, 110.0, 70.0)
                l_c = c2.number_input("Longueur Corps (cm)", 40.0, 130.0, 85.0)
                c_c = c3.number_input("Circ. Canon (cm)", 5.0, 20.0, 10.0)
                t_b = c4.number_input("Taille Bassin (cm)", 10.0, 45.0, 22.0)
                
                st.subheader("🥛 Morphologie Mamelle")
                m1, m2 = st.columns(2)
                m_d = m1.number_input("Diamètre Mamelle (cm)", 5.0, 40.0, 18.0)
                m_p = m2.number_input("Profondeur Mamelle (cm)", 5.0, 40.0, 15.0)
                
                if st.form_submit_button("Lancer Analyse & Enregistrer"):
                    ind = ScienceEngine.calculer_indice_scan(h_g, l_c)
                    db.execute_query("""INSERT INTO scanner_expert 
                        (brebis_id, hauteur_garrot, longueur_corps, circ_canon, taille_bassin, diametre_mamelle, profondeur_mamelle, indice_conformation, date_scan) 
                        VALUES (?,?,?,?,?,?,?,?,?)""", 
                        (target, h_g, l_c, c_c, t_b, m_d, m_p, ind, date.today()))
                    st.success(f"Analyse terminée. Indice de compacité : {ind}%")
            
            st.subheader("Historique Morphométrique")
            st.dataframe(db.fetch_all_as_df("SELECT * FROM scanner_expert ORDER BY date_scan DESC"))
        else:
            st.warning("Veuillez d'abord enregistrer des animaux dans le Registre.")

    # --- DASHBOARD ---
    elif choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - {view_user}")
        df = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df))
            c2.metric("Poids Moyen", f"{round(df['poids'].mean(), 1)} kg")
            st.plotly_chart(px.bar(df, x="identifiant_unique", y="poids", color="race", title="Poids du Troupeau"))
        else: st.info("Aucune donnée disponible.")

    # --- BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Qualité du Lait & Biochimie")
        
        df_adultes = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=? AND type_animal='Brebis Adulte'", (view_user,))
        if not df_adultes.empty:
            with st.form("labo"):
                target = st.selectbox("Brebis", df_adultes['identifiant_unique'])
                q, tb, tp = st.columns(3)
                q_v = q.number_input("Lait (L)", 0.1, 10.0, 1.5)
                tb_v = tb.number_input("TB (g/L)", 20.0, 90.0, 45.0)
                tp_v = tp.number_input("TP (g/L)", 20.0, 90.0, 38.0)
                if st.form_submit_button("Calculer ESD"):
                    esd = ScienceEngine.calculer_esd(tb_v, tp_v)
                    db.execute_query("INSERT INTO lait_biochimie (brebis_id, qte_lait, tb, tp, esd, date_controle, owner_id) VALUES (?,?,?,?,?,?,?)",
                                    (target, q_v, tb_v, tp_v, esd, date.today(), view_user))
                    st.success(f"ESD Calculé : {esd}")
            st.dataframe(db.fetch_all_as_df("SELECT * FROM lait_biochimie WHERE owner_id=?", (view_user,)))

    # --- NUTRITION ---
    elif choice == "🍲 Nutrition & Ration":
        st.title("🍲 Calculateur de Ration")
        # Logique de ration simple
        choix = st.multiselect("Aliments disponibles", list(TABLE_VALEURS.keys()), default=["Orge"])
        if choix:
            qtes = {a: st.number_input(f"Kg de {a}", 0.0, 5.0, 0.5) for a in choix}
            total_ufl = sum(qtes[a] * TABLE_VALEURS[a]["UFL"] for a in choix)
            st.metric("Total Énergie (UFL)", f"{round(total_ufl, 2)}")

    # --- STOCKS ---
    elif choice == "📦 Stocks":
        st.title("📦 Inventaire Silos")
        al = st.selectbox("Aliment", list(TABLE_VALEURS.keys()))
        q = st.number_input("Quantité (Qx)", 0.0, 500.0)
        if st.button("Mettre à jour Stock"):
            db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (view_user, al, q))
        st.dataframe(db.fetch_all_as_df("SELECT aliment, quantite_q FROM stocks WHERE owner_id=?", (view_user,)))

    # --- REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Gestion du Registre")
        with st.form("reg"):
            uid = st.text_input("ID Boucle")
            cat = st.selectbox("Type", ["Brebis Adulte", "Agnelle", "Bélier"])
            rac = st.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi"])
            pds = st.number_input("Poids (kg)", 10.0, 150.0, 60.0)
            if st.form_submit_button("Enregistrer"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, created_at) VALUES (?,?,?,?,?,?)",
                                (uid, view_user, rac, cat, pds, date.today()))
                st.success("Animal ajouté.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

def inject_demo_data(db, user):
    db.execute_query("INSERT OR IGNORE INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, created_at) VALUES (?,?,?,?,?,?)",
                    ("TEST_01", user, "Ouled Djellal", "Brebis Adulte", 75, date.today()))

if __name__ == "__main__":
    main()
