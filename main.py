"""
EXPERT OVIN DZ PRO - VERSION V24.GENOMIC_MASTER_FULL
--------------------------------------------------------
- Scanner : Morphométrie complète (Bassin, Canon, Mamelle) + Étallons
- Registre : Race libre, Age (Dentition/Date), Sexe, Catégories
- Modules : Bio-informatique, Labo, Nutrition et Stocks restaurés
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date
import io

# ============================================================================
# 1. DATABASE ENGINE (V24)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_v24_full.db"):
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
        "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT)",
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE,
            owner_id TEXT, race TEXT, sexe TEXT, categorie TEXT, poids REAL, 
            methode_age TEXT, valeur_age TEXT, pere_id TEXT, mere_id TEXT, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS scanner_expert (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, etalon TEXT,
            h_garrot REAL, l_corps REAL, l_bassin REAL, circ_canon REAL,
            m_diametre REAL, m_profondeur REAL, m_forme TEXT, date_scan DATE
        )""",
        """CREATE TABLE IF NOT EXISTS lait_biochimie (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            tb REAL, tp REAL, esd REAL, qte REAL, date_controle DATE, owner_id TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS stocks (
            owner_id TEXT, aliment TEXT, quantite_q REAL, prix_q REAL, PRIMARY KEY(owner_id, aliment)
        )"""
    ]
    for t in tables: db.execute_query(t)
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert')")

# ============================================================================
# 2. LOGIQUE SCIENTIFIQUE
# ============================================================================

TABLE_VALEURS = {"Orge": 1.0, "Son": 0.85, "Foin Luzerne": 0.65, "Maïs": 1.15}

class GenomicEngine:
    @staticmethod
    def calculer_esd(tb, tp):
        return round((tp * 0.85) + (tb * 0.1) + 5.5, 2)

    @staticmethod
    def safe_ebv(row):
        try:
            p = float(row['poids']) if row['poids'] else 50.0
            return round(p * 0.4, 1) # Index simplifié
        except: return 0.0

# ============================================================================
# 3. INTERFACE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V24 PRO", layout="wide", page_icon="🧬")
    
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
    st.sidebar.title(f"🧬 {role}")
    
    menu = ["📊 Dashboard", "🧬 Hub Bio-info", "📸 Scanner IA Expert", "🥛 Labo Biochimie", "🍲 Nutrition", "📦 Stocks", "📝 Registre"]
    choice = st.sidebar.radio("Navigation", menu)
    
    df_all = db.fetch_all_as_df(f"SELECT * FROM brebis WHERE owner_id='{user}'")

    # --- REGISTRE (Modifié selon vos demandes) ---
    if choice == "📝 Registre":
        st.title("📝 Enregistrement du Sujet")
        with st.form("reg_form"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle (Unique)")
            race = c1.text_input("Race (Optionnel - Laisser vide si inconnu)")
            sexe = c1.selectbox("Sexe", ["Femelle", "Mâle"])
            cat = c1.selectbox("Catégorie", ["Brebis Adulte", "Bélier", "Agnelle", "Agneau"])
            
            methode_age = c2.selectbox("Détermination de l'âge", ["Date de naissance", "Dentition (Remplacements)", "Mois approximatifs"])
            valeur_age = c2.text_input("Valeur (ex: 15/05/2024 ou 2 dents)")
            pds = c2.number_input("Poids initial (kg)", 5.0, 150.0, 45.0)
            
            p_id = st.text_input("ID Père (Bio-info)")
            m_id = st.text_input("ID Mère (Bio-info)")
            
            if st.form_submit_button("Inscrire au Registre"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, sexe, categorie, poids, methode_age, valeur_age, pere_id, mere_id, created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                                (uid, user, race, sexe, cat, pds, methode_age, valeur_age, p_id, m_id, date.today()))
                st.success(f"Sujet {uid} enregistré avec succès.")

    # --- SCANNER IA EXPERT (Détaillé) ---
    elif choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner Morphométrique & Phénotypique")
        if not df_all.empty:
            with st.form("scan_expert"):
                target = st.selectbox("Sélectionner le sujet", df_all['identifiant_unique'])
                etalon = st.radio("Étalon de mesure", ["Feuille A4", "Bâton 1 mètre", "Carte Bancaire"], horizontal=True)
                
                st.subheader("📏 Mensurations de Conformation")
                c1, c2, c3 = st.columns(3)
                h_g = c1.number_input("Hauteur au Garrot (cm)", 30.0, 120.0, 70.0)
                l_c = c2.number_input("Longueur du corps (cm)", 30.0, 150.0, 80.0)
                l_b = c3.number_input("Largeur du Bassin (cm)", 10.0, 50.0, 22.0)
                c_c = c1.number_input("Circonférence du Canon (cm)", 5.0, 20.0, 9.0)
                
                st.subheader("🥛 Caractères Phénotypiques Mamelle")
                m_d = c1.number_input("Diamètre Mamelle (cm)", 0.0, 40.0, 15.0)
                m_p = c2.number_input("Profondeur Mamelle (cm)", 0.0, 40.0, 12.0)
                m_f = c3.selectbox("Forme des Trayons", ["Cylindrique", "Conique", "Asymétrique", "Petit/Absent"])
                
                if st.form_submit_button("Enregistrer les données du Scanner"):
                    db.execute_query("INSERT INTO scanner_expert (brebis_id, etalon, h_garrot, l_corps, l_bassin, circ_canon, m_diametre, m_profondeur, m_forme, date_scan) VALUES (?,?,?,?,?,?,?,?,?,?)",
                                    (target, etalon, h_g, l_c, l_b, c_c, m_d, m_p, m_f, date.today()))
                    st.success("Analyse morphométrique sauvegardée.")
        else: st.warning("Veuillez d'abord remplir le Registre.")

    # --- HUB BIO-INFORMATIQUE (Sécurisé) ---
    elif choice == "🧬 Hub Bio-info":
        st.title("🧬 Laboratoire Bio-informatique")
        if df_all.empty:
            st.info("La base de données est vide. Les outils bio-informatiques s'activeront après l'ajout de sujets.")
        else:
            t1, t2 = st.tabs(["📊 Génétique", "🌐 NCBI"])
            with t1:
                st.dataframe(df_all[['identifiant_unique', 'race', 'sexe', 'pere_id', 'mere_id']])
                st.plotly_chart(px.pie(df_all, names='sexe', title="Répartition par Sexe"))
            with t2:
                st.link_button("Accéder à GeneBank (NCBI)", "https://www.ncbi.nlm.nih.gov/genbank/")

    # --- LABO BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Analyse Biochimique du Lait")
        if not df_all.empty:
            with st.form("labo"):
                target = st.selectbox("Sujet", df_all[df_all['sexe']=="Femelle"]['identifiant_unique'])
                tb = st.number_input("Taux Butyreux (g/L)", 10.0, 100.0, 45.0)
                tp = st.number_input("Taux Protéique (g/L)", 10.0, 100.0, 38.0)
                qte = st.number_input("Production Journalière (Litres)", 0.0, 10.0, 1.5)
                if st.form_submit_button("Valider l'analyse"):
                    esd = GenomicEngine.calculer_esd(tb, tp)
                    db.execute_query("INSERT INTO lait_biochimie (brebis_id, tb, tp, esd, qte, date_controle, owner_id) VALUES (?,?,?,?,?,?,?)",
                                    (target, tb, tp, esd, qte, date.today(), user))
                    st.success(f"Analyse terminée. ESD: {esd}")
        else: st.info("Aucun sujet femelle trouvé.")

    # --- NUTRITION ---
    elif choice == "🍲 Nutrition":
        st.title("🍲 Formulation de la Ration")
        choix = st.multiselect("Sélectionner les aliments", list(TABLE_VALEURS.keys()))
        if choix:
            st.write("Détails de la ration :")
            total_ufl = 0
            for a in choix:
                q = st.number_input(f"Quantité de {a} (kg)", 0.1, 10.0, 0.5)
                total_ufl += q * TABLE_VALEURS[a]
            st.metric("Énergie Totale (UFL)", round(total_ufl, 2))

    # --- STOCKS ---
    elif choice == "📦 Stocks":
        st.title("📦 Gestion des Stocks Silos")
        with st.form("stock_add"):
            al = st.selectbox("Aliment", list(TABLE_VALEURS.keys()))
            quant = st.number_input("Quantité en Quintaux (Qx)", 0.0, 500.0, 10.0)
            if st.form_submit_button("Mettre à jour le stock"):
                db.execute_query("INSERT OR REPLACE INTO stocks (owner_id, aliment, quantite_q) VALUES (?,?,?)", (user, al, quant))
        
        st.subheader("Inventaire Actuel")
        st.dataframe(db.fetch_all_as_df(f"SELECT aliment, quantite_q FROM stocks WHERE owner_id='{user}'"))

    # --- DASHBOARD ---
    elif choice == "📊 Dashboard":
        st.title(f"📊 Dashboard Expert - {user}")
        if not df_all.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df_all))
            c2.metric("Poids Moyen", f"{round(df_all['poids'].mean(), 1)} kg")
            st.plotly_chart(px.bar(df_all, x="identifiant_unique", y="poids", color="categorie"))
        else: st.info("Bienvenue. Commencez par enregistrer vos animaux dans le Registre.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
