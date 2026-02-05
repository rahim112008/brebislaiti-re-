"""
EXPERT OVIN DZ PRO - VERSION V33.EXPERT_COMMUNICATION_HUB
--------------------------------------------------------
Nouveautés :
1. Système de Notifications : L'Expert envoie des messages aux éleveurs.
2. Modules Scientifiques Complets : Biochimie, Nutrition et Stocks activés.
3. Annuaire Clients : Vue 360° du cheptel par éleveur.
4. Correctif Sécurité : Gestion des index et des dossiers clients.
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
# 1. MOTEUR DE BASE DE DONNÉES (V33)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_pro_v33.db"):
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
        "CREATE TABLE IF NOT EXISTS brebis (id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE, owner_id TEXT, race TEXT, sexe TEXT, categorie TEXT, poids REAL, created_at DATE)",
        "CREATE TABLE IF NOT EXISTS scanner_expert (id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT, h_garrot REAL, l_corps REAL, l_bassin REAL, circ_canon REAL, m_diametre REAL, m_forme TEXT, status TEXT DEFAULT 'En attente', date_scan DATE)",
        "CREATE TABLE IF NOT EXISTS labo_lait (id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT, tb REAL, tp REAL, lactose REAL, esd REAL, date_analyse DATE)",
        "CREATE TABLE IF NOT EXISTS messages (id INTEGER PRIMARY KEY AUTOINCREMENT, dest_user TEXT, sender TEXT, content TEXT, is_read INTEGER DEFAULT 0, created_at DATETIME)",
        "CREATE TABLE IF NOT EXISTS stocks (owner_id TEXT, aliment TEXT, quantite REAL, PRIMARY KEY(owner_id, aliment))"
    ]
    for t in tables: db.execute_query(t)
    
    # Comptes par défaut
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert', ?)", (datetime.now(),))
    for u, p in [('Eleveur_Setif', 'setif2026'), ('Eleveur_Tiaret', 'tiaret2026')]:
        db.execute_query("INSERT OR IGNORE INTO users VALUES (?, ?, 'Eleveur', ?)", (u, p, datetime.now()))

# ============================================================================
# 2. INTERFACE PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin V33 PRO", layout="wide", page_icon="📡")
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_database(st.session_state.db)
    db = st.session_state.db

    if 'auth' not in st.session_state:
        st.title("🛡️ Station Master Ovin DZ")
        u, p = st.text_input("Username"), st.text_input("Password", type="password")
        if st.button("Connexion"):
            res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role
    st.sidebar.title(f"🧬 Mode {role}")

    # --- CENTRE DE MESSAGERIE (NOTIFICATION) ---
    msgs = db.fetch_all_as_df(f"SELECT * FROM messages WHERE dest_user='{user}' AND is_read=0")
    if not msgs.empty:
        st.sidebar.warning(f"📩 {len(msgs)} Nouveau(x) message(s) Expert")
        if st.sidebar.button("Lire les messages"):
            for _, m in msgs.iterrows():
                st.sidebar.info(f"De {m['sender']}: {m['content']}")
            db.execute_query(f"UPDATE messages SET is_read=1 WHERE dest_user='{user}'")

    # --- NAVIGATION ---
    menu = ["📊 Dashboard", "🏢 Gestion Clients", "📸 Scanner IA Expert", "🧬 Hub Bio-info", "🥛 Labo Biochimie", "🍲 Nutrition", "📦 Stocks", "📝 Registre"]
    if role != "Expert": menu.remove("🏢 Gestion Clients")
    choice = st.sidebar.radio("Navigation", menu)

    # Filtrage Data
    q = "SELECT * FROM brebis" if role == "Expert" else f"SELECT * FROM brebis WHERE owner_id='{user}'"
    df_view = db.fetch_all_as_df(q)

    # --- 🏢 GESTION CLIENTS & MESSAGERIE (EXPERT) ---
    if choice == "🏢 Gestion Clients" and role == "Expert":
        st.title("🏢 Espace Expert & Communication")
        t1, t2, t3 = st.tabs(["👥 Clients", "✅ Validations", "✉️ Envoyer Message"])
        
        with t1:
            df_c = db.fetch_all_as_df("SELECT u.username, COUNT(b.id) as 'Animaux' FROM users u LEFT JOIN brebis b ON u.username=b.owner_id WHERE u.role='Eleveur' GROUP BY u.username")
            st.table(df_c)
        
        with t2:
            p = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='En attente'")
            if not p.empty:
                for _, r in p.iterrows():
                    if st.button(f"Valider {r['brebis_id']} ({r['owner_id']})"):
                        db.execute_query("UPDATE scanner_expert SET status='Validé' WHERE id=?", (r['id'],))
                        st.rerun()
            else: st.success("Aucun scan en attente.")
            
        with t3:
            st.subheader("Communiquer avec un éleveur")
            dest = st.selectbox("Destinataire", df_c['username'].tolist())
            txt = st.text_area("Message (ex: Votre scan est rejeté, l'étalon n'est pas visible)")
            if st.button("Envoyer la notification"):
                db.execute_query("INSERT INTO messages (dest_user, sender, content, created_at) VALUES (?,?,?,?)",
                                (dest, user, txt, datetime.now()))
                st.success("Message envoyé !")

    # --- 📊 DASHBOARD ---
    elif choice == "📊 Dashboard":
        st.title("📊 Cockpit de Pilotage")
        if not df_view.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Cheptel", len(df_view))
            c2.metric("Poids Moyen", f"{round(df_view['poids'].mean(), 1)} kg")
            
            # Radar Chart
            df_val = db.fetch_all_as_df("SELECT * FROM scanner_expert WHERE status='Validé'")
            if len(df_val['brebis_id'].unique()) >= 2:
                st.subheader("⚖️ Comparateur de Conformation")
                opts = df_val['brebis_id'].unique().tolist()
                s1 = st.selectbox("Sujet A", opts, index=0)
                s2 = st.selectbox("Sujet B", opts, index=1)
                
                def get_m(sid):
                    r = df_val[df_val['brebis_id']==sid].iloc[-1]
                    return [r['h_garrot'], r['l_bassin'], r['circ_canon'], r['m_diametre'], 50]
                
                cats = ['H.Garrot', 'L.Bassin', 'Canon', 'Mamelle', 'Poids']
                fig = go.Figure()
                fig.add_trace(go.Scatterpolar(r=get_m(s1), theta=cats, fill='toself', name=s1))
                fig.add_trace(go.Scatterpolar(r=get_m(s2), theta=cats, fill='toself', name=s2))
                st.plotly_chart(fig)
            else: st.info("Attente de scans validés pour le Radar.")
        else: st.info("Base vide.")

    # --- 📸 SCANNER IA DÉTAILLÉ ---
    elif choice == "📸 Scanner IA Expert":
        st.title("📸 Scanner IA & Phénotypage")
        
        if not df_view.empty:
            with st.form("scan"):
                target = st.selectbox("Sujet", df_view['identifiant_unique'])
                c1, c2 = st.columns(2)
                hg = c1.number_input("H. Garrot (cm)", 40.0, 110.0, 70.0)
                lb = c2.number_input("L. Bassin (cm)", 10.0, 45.0, 22.0)
                cc = c1.number_input("Canon (cm)", 5.0, 20.0, 10.0)
                md = c2.number_input("Mamelle Diamètre (cm)", 5.0, 50.0, 15.0)
                if st.form_submit_button("Envoyer pour Expertise"):
                    db.execute_query("INSERT INTO scanner_expert (brebis_id, owner_id, h_garrot, l_bassin, circ_canon, m_diametre, date_scan) VALUES (?,?,?,?,?,?,?)",
                                    (target, user, hg, lb, cc, md, date.today()))
                    st.success("Données envoyées.")
        else: st.warning("Registre vide.")

    # --- 🥛 LABO BIOCHIMIE ---
    elif choice == "🥛 Labo Biochimie":
        st.title("🥛 Analyse Qualité Laitière")
        if not df_view.empty:
            with st.form("lab"):
                target = st.selectbox("Animal", df_view['identifiant_unique'])
                tb = st.number_input("Taux Butyreux (g/L)", 30.0, 100.0, 50.0)
                tp = st.number_input("Taux Protéique (g/L)", 30.0, 100.0, 45.0)
                lac = st.number_input("Lactose (g/L)", 30.0, 60.0, 48.0)
                if st.form_submit_button("Calculer & Enregistrer"):
                    esd = round(tp + lac + 0.7, 2)
                    db.execute_query("INSERT INTO labo_lait (brebis_id, owner_id, tb, tp, lactose, esd, date_analyse) VALUES (?,?,?,?,?,?,?)",
                                    (target, user, tb, tp, lac, esd, date.today()))
                    st.success(f"ESD calculé : {esd}")
            st.dataframe(db.fetch_all_as_df(f"SELECT * FROM labo_lait WHERE owner_id='{user}'"))

    # --- 🍲 NUTRITION & 📦 STOCKS ---
    elif choice == "🍲 Nutrition":
        st.title("🍲 Rationnement")
        if not df_view.empty:
            target = st.selectbox("Sujet", df_view['identifiant_unique'])
            p = df_view[df_view['identifiant_unique']==target]['poids'].values[0]
            st.metric("Besoin Maintenance (UFL)", round(0.035 * p**0.75, 2))
            st.info("Recommandation : 60% Fourrage / 40% Concentré")

    elif choice == "📦 Stocks":
        st.title("📦 Magasin")
        with st.form("st"):
            ali = st.selectbox("Aliment", ["Foin", "Orge", "Son"])
            qte = st.number_input("Quantité (Q)", 0.0, 500.0, 10.0)
            if st.form_submit_button("MAJ Stock"):
                db.execute_query("INSERT OR REPLACE INTO stocks VALUES (?,?,?)", (user, ali, qte))
        st.table(db.fetch_all_as_df(f"SELECT * FROM stocks WHERE owner_id='{user}'"))

    # --- 📝 REGISTRE ---
    elif choice == "📝 Registre":
        st.title("📝 Registre du Cheptel")
        with st.form("reg"):
            uid, rac = st.text_input("ID Boucle"), st.text_input("Race")
            pds = st.number_input("Poids (kg)", 10.0, 150.0, 55.0)
            if st.form_submit_button("Ajouter"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, poids, created_at) VALUES (?,?,?,?,?)",
                                (uid, user, rac, pds, date.today()))
                st.success("Animal ajouté.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
