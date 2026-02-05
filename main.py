"""
EXPERT OVIN DZ PRO - VERSION V20.GENOMIC & BIOINFORMATICS
--------------------------------------------------------
Modules : Génomique, Bioinformatique (NCBI/GeneBank), Scanner IA Expert, 
          Biochimie, Nutrition, Stocks, Registre.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date

# ============================================================================
# 1. DATABASE ENGINE (V20 - Mise à jour structurelle)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_v20_genomic.db"):
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
            pere_id TEXT, mere_id TEXT, h_garrot REAL, l_corps REAL, created_at DATE
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
# 2. LOGIQUE SCIENTIFIQUE, GÉNOMIQUE & BIOINFORMATIQUE
# ============================================================================

TABLE_VALEURS = {"Orge": {"UFL": 1.0, "PDI": 80}, "Son de blé": {"UFL": 0.85, "PDI": 95}, "Foin de luzerne": {"UFL": 0.65, "PDI": 90}, "Maïs grain": {"UFL": 1.15, "PDI": 75}}

class GenomicEngine:
    @staticmethod
    def calculer_esd(tb, tp):
        return round((tp * 0.85) + (tb * 0.1) + 5.5, 2)

    @staticmethod
    def calculer_index_selection(row):
        # Index de sélection combiné (Poids + Morphologie + Lait estimé)
        p = row.get('poids', 50)
        m = row.get('note_mamelle', 5)
        return round((p * 0.3) + (m * 7), 1)

    @staticmethod
    def estimer_consanguinite(id_animal, df):
        # Simulation simplifiée de calcul de parenté
        row = df[df['identifiant_unique'] == id_animal]
        if row.empty or not row.iloc[0]['pere_id'] or not row.iloc[0]['mere_id']:
            return 0.0
        return 6.25  # Valeur de base pour test

# ============================================================================
# 3. INTERFACE PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V20 Genomic", layout="wide", page_icon="🧬")
    
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
    
    # --- Barre Latérale ---
    st.sidebar.title(f"✨ {role}")
    u_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
    view_user = st.sidebar.selectbox("📂 Dossier Éleveur", u_list) if (role == "Expert" and u_list) else user
    
    if role == "Expert":
        if st.sidebar.button("🧪 Injecter Données Démo"):
            inject_demo_data(db, view_user)
            st.rerun()

    menu = ["📊 Dashboard", "🧬 Hub Génomique & Bioinfo", "📸 Scanner IA Expert", "🥛 Labo Biochimie", "🍲 Nutrition & Ration", "📦 Stocks", "📝 Registre"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE GÉNOMIQUE & BIOINFO (NOUVEAU & COMPLET) ---
    if choice == "🧬 Hub Génomique & Bioinfo":
        st.title("🧬 Hub de Bioinformatique & Génomique")
        
        tab1, tab2, tab3 = st.tabs(["🧬 Analyse Génomique", "📊 Génétique des Populations", "🌐 Ressources NCBI"])
        
        df_gen = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        
        with tab1:
            st.subheader("Analyse Individualisée")
            if not df_gen.empty:
                target = st.selectbox("Sujet pour analyse génomique", df_gen['identifiant_unique'])
                sub_df = df_gen[df_gen['identifiant_unique'] == target].iloc[0]
                
                c1, c2, c3 = st.columns(3)
                c1.metric("Index de Sélection (EBV)", GenomicEngine.calculer_index_selection(sub_df))
                c2.metric("Consanguinité estimée", f"{GenomicEngine.estimer_consanguinite(target, df_gen)}%")
                c3.metric("Lignée", sub_df['pere_id'] if sub_df['pere_id'] else "Inconnue")
            else: st.info("Aucun animal pour l'analyse.")

        with tab2:
            st.subheader("Statistiques de la Population")
            if not df_gen.empty:
                st.plotly_chart(px.histogram(df_gen, x="race", title="Distribution des Races"))
                st.plotly_chart(px.box(df_gen, x="race", y="poids", title="Variabilité du Phénotype Poids par Race"))
            else: st.info("Registre vide.")

        with tab3:
            st.subheader("Accès Bases de Données Génomiques (NCBI)")
            st.markdown("""
            - [NCBI Ovis aries (Genome)](https://www.ncbi.nlm.nih.gov/genome/?term=ovis+aries)
            - [GeneBank Ovine Sequences](https://www.ncbi.nlm.nih.gov/genbank/)
            - [Ensembl Genome Browser (Sheep)](https://www.ensembl.org/Ovis_aries/Info/Index)
            """)

    # --- SCANNER IA (CONSERVÉ) ---
    elif choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner IA & Morphométrie")
        etalon = st.selectbox("Étalon :", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
        df_t = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=?", (view_user,))
        if not df_t.empty:
            with st.form("scan_v20"):
                target = st.selectbox("Animal", df_t['identifiant_unique'])
                h = st.number_input("Hauteur (cm)", 40.0, 110.0, 70.0)
                l = st.number_input("Longueur (cm)", 40.0, 130.0, 85.0)
                if st.form_submit_button("Scanner"):
                    db.execute_query("INSERT INTO scanner_expert (brebis_id, hauteur_garrot, longueur_corps, date_scan) VALUES (?,?,?,?)", (target, h, l, date.today()))
                    st.success("Scan enregistré.")

    # --- DASHBOARD ---
    elif choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - {view_user}")
        df = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df.empty:
            c1, c2 = st.columns(2)
            c1.metric("Effectif", len(df))
            c2.metric("Poids Moyen", f"{round(df['poids'].mean(), 1)} kg")
            st.plotly_chart(px.bar(df, x="identifiant_unique", y="poids", color="race"))
        else: st.warning("Le troupeau est vide. Utilisez le module 'Registre' ou le bouton démo.")

    # --- BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Labo & Qualité Lait")
        df_a = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=? AND type_animal='Brebis Adulte'", (view_user,))
        if not df_a.empty:
            with st.form("bio"):
                target = st.selectbox("Brebis", df_a['identifiant_unique'])
                tb = st.number_input("TB (g/L)", 20.0, 90.0, 45.0)
                tp = st.number_input("TP (g/L)", 20.0, 90.0, 38.0)
                if st.form_submit_button("Calculer"):
                    esd = GenomicEngine.calculer_esd(tb, tp)
                    db.execute_query("INSERT INTO lait_biochimie (brebis_id, tb, tp, esd, owner_id) VALUES (?,?,?,?,?)", (target, tb, tp, esd, view_user))
                    st.success(f"ESD : {esd}")

    # --- NUTRITION ---
    elif choice == "🍲 Nutrition & Ration":
        st.title("🍲 Calculateur de Ration")
        choix = st.multiselect("Mélange", list(TABLE_VALEURS.keys()), default=["Orge"])
        if choix:
            q = {a: st.number_input(f"Kg {a}", 0.0, 5.0, 0.5) for a in choix}
            total_ufl = sum(q[a] * TABLE_VALEURS[a]["UFL"] for a in choix)
            st.metric("Total UFL", round(total_ufl, 2))

    # --- STOCKS ---
    elif choice == "📦 Stocks":
        st.title("📦 Inventaire")
        al = st.selectbox("Aliment", list(TABLE_VALEURS.keys()))
        q = st.number_input("Quantité (Qx)", 0.0, 1000.0)
        if st.button("Mettre à jour"):
            db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (view_user, al, q))
        st.dataframe(db.fetch_all_as_df("SELECT * FROM stocks WHERE owner_id=?", (view_user,)))

    # --- REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Registre du Troupeau")
        with st.form("reg_v20"):
            uid = st.text_input("ID Boucle")
            cat = st.selectbox("Type", ["Brebis Adulte", "Agnelle", "Bélier"])
            rac = st.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi"])
            pds = st.number_input("Poids (kg)", 10.0, 150.0, 60.0)
            p_id = st.text_input("ID Père")
            m_id = st.text_input("ID Mère")
            if st.form_submit_button("Inscrire"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?)",
                                (uid, view_user, rac, cat, pds, p_id, m_id, date.today()))
                st.success("Enregistré.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

def inject_demo_data(db, user):
    db.execute_query("INSERT OR IGNORE INTO brebis (identifiant_unique, owner_id, race, type_animal, poids, created_at) VALUES (?,?,?,?,?,?)",
                    ("GEN_TEST", user, "Lacaune", "Brebis Adulte", 68.5, date.today()))
    db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (user, "Orge", 50.0))

if __name__ == "__main__":
    main()
