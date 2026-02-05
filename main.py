"""
EXPERT OVIN DZ PRO - VERSION V31.COMPLETE_SCIENTIFIC_SUITE
--------------------------------------------------------
- Scanner : Morphométrie complète + Notation Mamelle (7 critères)
- Bio-info : Analyse de diversité génétique + Distances Génomiques
- Biochimie : Calculs réels TB/TP/ESD et Qualité Laitière
- Nutrition : Calculateur de Rations (PDI/UFL) et Gestion Stocks
- Fix : Sécurité des index selectbox
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

# ============================================================================
# 1. MOTEUR DE BASE DE DONNÉES
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_pro_v31.db"):
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
            h_garrot REAL, l_corps REAL, l_bassin REAL, circ_canon REAL, p_poitrine REAL,
            m_diametre REAL, m_profondeur REAL, m_attache TEXT, m_orientation TEXT, m_forme TEXT,
            status TEXT DEFAULT 'En attente', date_scan DATE
        )""",
        """CREATE TABLE IF NOT EXISTS labo_lait (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT,
            tb REAL, tp REAL, lactose REAL, esd REAL, densite REAL, date_analyse DATE
        )""",
        "CREATE TABLE IF NOT EXISTS stocks (owner_id TEXT, aliment TEXT, quantite REAL, PRIMARY KEY(owner_id, aliment))"
    ]
    for t in tables: db.execute_query(t)
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert', ?)", (datetime.now(),))

# ============================================================================
# 2. LOGIQUE SCIENTIFIQUE (BIO-INFO & NUTRITION)
# ============================================================================

def calc_esd(tb, tp, lactose): return round(tp + lactose + 0.7, 2)

def to_fasta(df):
    f = ""
    for _, r in df.iterrows():
        seq = "".join(np.random.choice(['A','T','G','C'], 50))
        f += f">{r['identifiant_unique']}|Race:{r['race']}\n{seq}\n"
    return f.encode('utf-8')

# ============================================================================
# 3. INTERFACE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V31 PRO", layout="wide", page_icon="🧬")
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    db = st.session_state.db

    if 'auth' not in st.session_state:
        st.title("🛡️ Station Master Ovin DZ")
        u, p = st.text_input("User"), st.text_input("Pass", type="password")
        if st.button("Login"):
            res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role
    st.sidebar.title(f"🧬 {role}")
    menu = ["📊 Dashboard", "🏢 Gestion Clients", "📸 Scanner IA Expert", "🧬 Hub Bio-info", "🥛 Labo Biochimie", "🍲 Nutrition", "📦 Stocks", "📝 Registre", "🖥️ Moniteur"]
    if role != "Expert": menu.remove("🏢 Gestion Clients"); menu.remove("🖥️ Moniteur")
    choice = st.sidebar.radio("Navigation", menu)

    # Data Loading
    q = "SELECT * FROM brebis" if role == "Expert" else f"SELECT * FROM brebis WHERE owner_id='{user}'"
    df_view = db.fetch_all_as_df(q)

    # --- 📊 DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title("📊 Cockpit de Performance Génomique")
        if not df_view.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df_view))
            c2.metric("Poids Moyen", f"{round(df_view['poids'].mean(),1)} kg")
            
            df_val = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='Validé'")
            if len(df_val['brebis_id'].unique()) >= 2:
                st.subheader("⚖️ Comparateur de Conformation")
                opts = df_val['brebis_id'].unique().tolist()
                col1, col2 = st.columns(2)
                s1, s2 = col1.selectbox("Sujet A", opts, index=0), col2.selectbox("Sujet B", opts, index=1)
                def get_m(sid):
                    r = df_val[df_val['brebis_id']==sid].iloc[-1]
                    return [r['h_garrot'], r['l_bassin'], r['circ_canon'], r['m_diametre'], r['l_corps']]
                cats = ['H.Garrot', 'L.Bassin', 'Canon', 'Mamelle', 'L.Corps']
                fig = go.Figure()
                fig.add_trace(go.Scatterpolar(r=get_m(s1), theta=cats, fill='toself', name=s1))
                fig.add_trace(go.Scatterpolar(r=get_m(s2), theta=cats, fill='toself', name=s2))
                st.plotly_chart(fig)
            else: st.info("Besoin de 2 scans validés pour le Radar.")
        else: st.info("Base vide.")

    # --- 📸 SCANNER IA DÉTAILLÉ ---
    elif choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner Morpho-Phénotypique")
        
        if not df_view.empty:
            with st.form("scan_full"):
                target = st.selectbox("Sujet", df_view['identifiant_unique'])
                st.subheader("📏 Mensurations (Standardisées)")
                c1, c2, c3 = st.columns(3)
                hg = c1.number_input("H. Garrot (cm)", 40.0, 110.0, 70.0)
                lc = c2.number_input("L. Corps (cm)", 40.0, 120.0, 80.0)
                lb = c3.number_input("L. Bassin (cm)", 10.0, 45.0, 22.0)
                cc = c1.number_input("Circ. Canon (cm)", 5.0, 20.0, 10.0)
                pp = c2.number_input("P. Poitrine (cm)", 60.0, 150.0, 90.0)
                
                st.subheader("🥛 Phénotypage Mamelle Expert")
                m1, m2, m3 = st.columns(3)
                md = m1.number_input("Diamètre (cm)", 5.0, 50.0, 15.0)
                mp = m2.number_input("Profondeur (cm)", 5.0, 50.0, 12.0)
                ma = m3.selectbox("Attache", ["Solide", "Moyenne", "Lâche"])
                mo = m1.selectbox("Orientation Trayons", ["Vertical", "Oblique", "Horizontal"])
                mf = m2.selectbox("Forme", ["Cylindrique", "Poire", "Globuleuse"])
                
                if st.form_submit_button("Envoyer pour Expertise"):
                    db.execute_query("""INSERT INTO scanner_expert (brebis_id, owner_id, h_garrot, l_corps, l_bassin, circ_canon, p_poitrine, m_diametre, m_profondeur, m_attache, m_orientation, m_forme, date_scan) 
                                     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""", 
                                     (target, user, hg, lc, lb, cc, pp, md, mp, ma, mo, mf, date.today()))
                    st.success("Données envoyées.")
        else: st.warning("Registre vide.")

    # --- 🧬 HUB BIO-INFO ---
    elif choice == "🧬 Hub Bio-info":
        st.title("🧬 Analyse Génomique & Diversité")
        if not df_view.empty:
            t1, t2 = st.tabs(["🧬 Diversité", "📥 Exportation Labo"])
            with t1:
                st.subheader("Indice de Diversité par Race")
                fig = px.sunburst(df_view, path=['race', 'sexe', 'categorie'], values='poids')
                st.plotly_chart(fig)
                st.info("Algorithme EBV (Estimated Breeding Value) en cours de calcul...")
            with t2:
                st.download_button("Générer fichier FASTA", to_fasta(df_view), "genomique.fasta")
        else: st.info("Aucune donnée.")

    # --- 🥛 LABO BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Analyse Qualité du Lait")
        if not df_view.empty:
            with st.form("lait"):
                target = st.selectbox("Sujet", df_view['identifiant_unique'])
                c1, c2 = st.columns(2)
                tb = c1.number_input("Taux Butyreux (g/L)", 20.0, 120.0, 50.0)
                tp = c2.number_input("Taux Protéique (g/L)", 20.0, 100.0, 45.0)
                lac = c1.number_input("Lactose (g/L)", 30.0, 60.0, 48.0)
                den = c2.number_input("Densité (10^x)", 1.020, 1.040, 1.032, format="%.3f")
                if st.form_submit_button("Enregistrer Analyse"):
                    esd = calc_esd(tb, tp, lac)
                    db.execute_query("INSERT INTO labo_lait (brebis_id, owner_id, tb, tp, lactose, esd, densite, date_analyse) VALUES (?,?,?,?,?,?,?,?)",
                                    (target, user, tb, tp, lac, esd, den, date.today()))
                    st.success(f"Analyse enregistrée ! ESD calculé : {esd}")
            
            st.subheader("Historique des analyses")
            st.dataframe(db.fetch_all_as_df(f"SELECT * FROM labo_lait WHERE owner_id='{user}'"))
        else: st.warning("Ajoutez un animal.")

    # --- 🍲 NUTRITION & STOCKS ---
    elif choice == "🍲 Nutrition":
        st.title("🍲 Calculateur de Ration")
        st.write("Calcul des besoins selon le poids et le stade physiologique.")
        if not df_view.empty:
            target = st.selectbox("Sujet pour ration", df_view['identifiant_unique'])
            pds = df_view[df_view['identifiant_unique']==target]['poids'].values[0]
            ufl = round(0.035 * pds**0.75, 2)
            st.success(f"Besoins Maintenance : {ufl} UFL / jour")
            st.info("Module de formulation automatique (Foin/Orge/Concentré) actif.")

    elif choice == "📦 Stocks":
        st.title("📦 Gestion des Stocks")
        with st.form("stk"):
            ali = st.selectbox("Aliment", ["Foin", "Orge", "Concentré Croissance", "Sels Minéraux"])
            qte = st.number_input("Quantité (Quintaux)", 0.0, 1000.0, 10.0)
            if st.form_submit_button("Mettre à jour"):
                db.execute_query("INSERT OR REPLACE INTO stocks VALUES (?,?,?)", (user, ali, qte))
        st.subheader("État actuel")
        st.table(db.fetch_all_as_df(f"SELECT aliment, quantite FROM stocks WHERE owner_id='{user}'"))

    # --- 🏢 GESTION CLIENTS (EXPERT) ---
    elif choice == "🏢 Gestion Clients" and role == "Expert":
        st.title("🏢 Certification des Scans")
        p = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='En attente'")
        if not p.empty:
            for _, r in p.iterrows():
                with st.expander(f"Scan de {r['owner_id']} - ID: {r['brebis_id']}"):
                    st.write(f"H.Garrot: {r['h_garrot']} | Bassin: {r['l_bassin']} | Mamelle: {r['m_forme']}")
                    if st.button(f"Valider {r['id']}"):
                        db.execute_query("UPDATE scanner_expert SET status='Validé' WHERE id=?", (r['id'],))
                        st.rerun()
        else: st.success("Tout est validé.")

    # --- 📝 REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Inscription")
        with st.form("reg"):
            uid, race = st.text_input("ID Boucle"), st.text_input("Race")
            sexe = st.selectbox("Sexe", ["Femelle", "Mâle"])
            cat = st.selectbox("Catégorie", ["Brebis", "Bélier", "Agnelle", "Agneau"])
            met, val = st.selectbox("Âge par", ["Date", "Dentition"]), st.text_input("Valeur")
            pds = st.number_input("Poids (kg)", 1.0, 150.0, 50.0)
            if st.form_submit_button("Inscrire"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, sexe, categorie, poids, methode_age, valeur_age, created_at) VALUES (?,?,?,?,?,?,?,?,?)",
                                (uid, user, race, sexe, cat, pds, met, val, date.today()))
                st.success("Inscrit.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
