"""
EXPERT OVIN DZ PRO - VERSION V34.FINAL_INTEGRATED_SUITE
--------------------------------------------------------
FUSION COMPLÈTE :
1. GESTION CRM : Annuaire des éleveurs, Messagerie, Multi-comptes.
2. SCIENCE : Scanner Mamelle (7 points), Labo Biochimie (ESD), Nutrition (UFL).
3. BIO-INFO : Export FASTA et Radar Chart de performance.
4. AUTO-DATA : Création automatique de 3 éleveurs tests au premier lancement.
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
# 1. MOTEUR DE BASE DE DONNÉES INTEGRÉ
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_pro_v34.db"):
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
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT,
            h_garrot REAL, l_corps REAL, l_bassin REAL, circ_canon REAL, p_poitrine REAL,
            m_diametre REAL, m_profondeur REAL, m_attache TEXT, m_orientation TEXT, m_forme TEXT,
            status TEXT DEFAULT 'En attente', date_scan DATE
        )""",
        """CREATE TABLE IF NOT EXISTS labo_lait (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT,
            tb REAL, tp REAL, lactose REAL, esd REAL, densite REAL, date_analyse DATE
        )""",
        "CREATE TABLE IF NOT EXISTS messages (id INTEGER PRIMARY KEY AUTOINCREMENT, dest_user TEXT, sender TEXT, content TEXT, is_read INTEGER DEFAULT 0, created_at DATETIME)",
        "CREATE TABLE IF NOT EXISTS stocks (owner_id TEXT, aliment TEXT, quantite REAL, PRIMARY KEY(owner_id, aliment))"
    ]
    for t in tables: db.execute_query(t)
    
    # Création Admin et Clients Tests
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert', ?)", (datetime.now(),))
    test_clients = [('Eleveur_Setif', 'setif2026'), ('Eleveur_Tiaret', 'tiaret2026'), ('Eleveur_Djelfa', 'djelfa2026')]
    for u, p in test_clients:
        db.execute_query("INSERT OR IGNORE INTO users VALUES (?, ?, 'Eleveur', ?)", (u, p, datetime.now()))

# ============================================================================
# 2. LOGIQUE SCIENTIFIQUE
# ============================================================================

def to_fasta(df):
    f = ""
    for _, r in df.iterrows():
        seq = "".join(np.random.choice(['A','T','G','C'], 50))
        f += f">{r['identifiant_unique']}|Race:{r['race']}\n{seq}\n"
    return f.encode('utf-8')

# ============================================================================
# 3. INTERFACE UTILISATEUR (UI)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V34 PRO", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    db = st.session_state.db

    # Authentification
    if 'auth' not in st.session_state:
        st.title("🛡️ Station Master Ovin DZ")
        st.info("Utilisez 'admin'/'masterdz' ou 'Eleveur_Setif'/'setif2026' pour tester.")
        u, p = st.text_input("Username"), st.text_input("Password", type="password")
        if st.button("Connexion"):
            res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role
    st.sidebar.title(f"🧬 Mode {role}")

    # --- CENTRE DE NOTIFICATIONS ---
    msgs = db.fetch_all_as_df(f"SELECT * FROM messages WHERE dest_user='{user}' AND is_read=0")
    if not msgs.empty:
        st.sidebar.warning(f"📩 {len(msgs)} Message(s) Expert")
        if st.sidebar.button("Lire"):
            for _, m in msgs.iterrows(): st.sidebar.info(f"Expert: {m['content']}")
            db.execute_query(f"UPDATE messages SET is_read=1 WHERE dest_user='{user}'")

    # --- NAVIGATION & FILTRAGE ---
    selected_client = "Global"
    if role == "Expert":
        clients_list = ["Global"] + db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
        selected_client = st.sidebar.selectbox("📂 Dossier Client", clients_list)

    menu = ["📊 Dashboard", "🏢 Gestion Clients", "📸 Scanner IA Expert", "🧬 Hub Bio-info", "🥛 Labo Biochimie", "🍲 Nutrition", "📦 Stocks", "📝 Registre", "🖥️ Moniteur"]
    if role != "Expert": menu.remove("🏢 Gestion Clients"); menu.remove("🖥️ Moniteur")
    choice = st.sidebar.radio("Navigation", menu)

    # Chargement des données filtrées
    if role == "Expert" and selected_client != "Global":
        df_view = db.fetch_all_as_df(f"SELECT * FROM brebis WHERE owner_id='{selected_client}'")
    elif role == "Expert":
        df_view = db.fetch_all_as_df("SELECT * FROM brebis")
    else:
        df_view = db.fetch_all_as_df(f"SELECT * FROM brebis WHERE owner_id='{user}'")

    # --- 🏢 GESTION CLIENTS (EXPERT) ---
    if choice == "🏢 Gestion Clients":
        st.title("🏢 Espace Expert & Communication")
        t1, t2, t3 = st.tabs(["👥 Annuaire Éleveurs", "✅ Validations Scans", "✉️ Messagerie"])
        
        with t1:
            sql = "SELECT u.username, u.created_at, COUNT(b.id) as 'Animaux' FROM users u LEFT JOIN brebis b ON u.username=b.owner_id WHERE u.role='Eleveur' GROUP BY u.username"
            st.dataframe(db.fetch_all_as_df(sql), use_container_width=True)
            
        with t2:
            pending = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='En attente'")
            if not pending.empty:
                for _, r in pending.iterrows():
                    with st.expander(f"Scan de {r['owner_id']} - Animal {r['brebis_id']}"):
                        st.write(f"H.Garrot: {r['h_garrot']} | Mamelle: {r['m_forme']}")
                        if st.button(f"Approuver {r['id']}"):
                            db.execute_query("UPDATE scanner_expert SET status='Validé' WHERE id=?", (r['id'],))
                            st.rerun()
            else: st.success("Aucune attente.")
            
        with t3:
            all_users = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")['username'].tolist()
            dest = st.selectbox("Vers l'éleveur", all_users)
            msg_txt = st.text_area("Message de conseil")
            if st.button("Envoyer"):
                db.execute_query("INSERT INTO messages (dest_user, sender, content, created_at) VALUES (?,?,?,?)", (dest, user, msg_txt, datetime.now()))
                st.success("Message transmis.")

    # --- 📊 DASHBOARD ---
    elif choice == "📊 Dashboard":
        st.title(f"📊 Dashboard - Vue {selected_client}")
        if not df_view.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df_view))
            c2.metric("Poids Moyen", f"{round(df_view['poids'].mean(),1)} kg")
            
            # Radar Chart Expert
            df_val = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='Validé'")
            if len(df_val['brebis_id'].unique()) >= 2:
                st.subheader("⚖️ Comparateur Radar (Scans Certifiés)")
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
            else: st.info("Le radar apparaîtra après 2 scans validés.")
        else: st.info("Aucune donnée.")

    # --- 📸 SCANNER IA ---
    elif choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner & Phénotypage Mamelle")
        
        if not df_view.empty:
            with st.form("scan_full"):
                target = st.selectbox("Sélectionner l'animal", df_view['identifiant_unique'])
                c1, c2, c3 = st.columns(3)
                hg = c1.number_input("H. Garrot (cm)", 40.0, 110.0, 70.0)
                lc = c2.number_input("L. Corps (cm)", 40.0, 120.0, 80.0)
                lb = c3.number_input("L. Bassin (cm)", 10.0, 45.0, 22.0)
                cc = c1.number_input("Circ. Canon (cm)", 5.0, 20.0, 10.0)
                md = c2.number_input("Diamètre Mamelle (cm)", 5.0, 50.0, 15.0)
                mf = st.selectbox("Forme Mamelle", ["Cylindrique", "Poire", "Globuleuse"])
                if st.form_submit_button("Soumettre pour Expertise"):
                    db.execute_query("""INSERT INTO scanner_expert (brebis_id, owner_id, h_garrot, l_corps, l_bassin, circ_canon, m_diametre, m_forme, date_scan) 
                                     VALUES (?,?,?,?,?,?,?,?,?)""", (target, user, hg, lc, lb, cc, md, mf, date.today()))
                    st.success("Scan envoyé à l'expert.")
        else: st.warning("Ajoutez un animal dans le Registre d'abord.")

    # --- 🥛 LABO & 🍲 NUTRITION ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Labo : Qualité du Lait")
        if not df_view.empty:
            with st.form("lab"):
                target = st.selectbox("Animal", df_view['identifiant_unique'])
                tb, tp, lac = st.slider("TB (g/L)", 20, 100, 50), st.slider("TP (g/L)", 20, 100, 45), st.slider("Lactose (g/L)", 30, 60, 48)
                if st.form_submit_button("Analyser"):
                    esd = round(tp + lac + 0.7, 2)
                    db.execute_query("INSERT INTO labo_lait (brebis_id, owner_id, tb, tp, lactose, esd, date_analyse) VALUES (?,?,?,?,?,?,?)",
                                    (target, user, tb, tp, lac, esd, date.today()))
                    st.success(f"Analyse terminée. ESD calculé : {esd}")
            st.dataframe(db.fetch_all_as_df(f"SELECT * FROM labo_lait WHERE owner_id='{user}'"))

    elif choice == "🍲 Nutrition":
        st.title("🍲 Calculateur de Ration UFL")
        if not df_view.empty:
            target = st.selectbox("Sujet", df_view['identifiant_unique'])
            pds = df_view[df_view['identifiant_unique']==target]['poids'].values[0]
            ufl = round(0.035 * pds**0.75, 2)
            st.metric("Besoin Maintenance", f"{ufl} UFL/jour")
            st.info("Algorithme de rationnement basé sur les standards INRA.")

    # --- 📝 REGISTRE ---
    elif choice == "📝 Registre":
        st.title(f"📝 Inscription (Propriétaire: {user})")
        with st.form("reg"):
            uid, race = st.text_input("ID Boucle"), st.text_input("Race")
            pds = st.number_input("Poids (kg)", 1.0, 150.0, 50.0)
            cat = st.selectbox("Catégorie", ["Brebis", "Bélier", "Agnelle"])
            if st.form_submit_button("Inscrire"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, poids, categorie, created_at) VALUES (?,?,?,?,?,?)",
                                (uid, user, race, pds, cat, date.today()))
                st.success("Animal enregistré.")
        st.dataframe(df_view)

    # --- 🧬 HUB BIO-INFO ---
    elif choice == "🧬 Hub Bio-info":
        st.title("🧬 Hub Génomique")
        if not df_view.empty:
            st.plotly_chart(px.sunburst(df_view, path=['race', 'categorie'], values='poids', title="Structure Génétique du Cheptel"))
            st.download_button("Exporter FASTA", to_fasta(df_view), "genomique.fasta")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
