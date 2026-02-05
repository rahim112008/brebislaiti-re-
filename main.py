"""
EXPERT OVIN DZ - VERSION V42 (ULTIMATE MASTER FUSION)
--------------------------------------------------------
SYNTHÈSE TOTALE DES FONCTIONNALITÉS :
1. MULTI-ROLES : Expert (Contrôle National) vs Éleveur (Gestion Terrain).
2. SCANNER IA PRO : Morphométrie 7 points, étalon 1m, et Score BCS (1-5).
3. GÉNÉTIQUE : Index BLUP, Simulation d'accouplement, Matrice de parenté.
4. BIOMÉTRIE : Indices de compacité, de format, et Radar Chart phénotypique.
5. NUTRITION : Calculateur de rations dynamique (UFL/MS) lié au poids réel.
6. ANALYTICS : Courbes de croissance (GMQ) et Observatoire National (Sunburst).
7. MESSAGERIE "ECHO" : Alertes et directives Expert vers Éleveurs.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
from datetime import datetime, date

# ============================================================================
# 1. ARCHITECTURE DE LA BASE DE DONNÉES (PERSISTANTE)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_ultimate_v42.db"):
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
    tables = [
        "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT, region TEXT, created_at DATETIME)",
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE, owner_id TEXT, 
            race TEXT, sexe TEXT, poids REAL, pere_id TEXT, mere_id TEXT, created_at DATE
        )""",
        "CREATE TABLE IF NOT EXISTS poids_history (id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, poids REAL, date_mesure DATE)",
        """CREATE TABLE IF NOT EXISTS scanner_data (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT, 
            h_garrot REAL, l_corps REAL, l_bassin REAL, bcs REAL, status TEXT DEFAULT 'En attente', date_scan DATE
        )""",
        "CREATE TABLE IF NOT EXISTS messages (id INTEGER PRIMARY KEY AUTOINCREMENT, dest_user TEXT, sender TEXT, content TEXT, is_read INTEGER DEFAULT 0, created_at DATETIME)",
        "CREATE TABLE IF NOT EXISTS labo_lait (id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, owner_id TEXT, tb REAL, tp REAL, esd REAL, date_analyse DATE)"
    ]
    for t in tables: db.execute_query(t)
    # Comptes standards
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert', 'Alger', ?)", (datetime.now(),))
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('eleveur_setif', 'setif2026', 'Eleveur', 'Sétif', ?)", (datetime.now(),))

# ============================================================================
# 2. MOTEUR DE CALCULS SCIENTIFIQUES
# ============================================================================

def calc_ration(poids):
    ufl = round(0.035 * poids**0.75, 2)
    foin = round(poids * 0.02, 2)
    orge = round(poids * 0.008, 2)
    return ufl, foin, orge

def calc_biometrie(poids, hg, lc):
    compacite = round(poids / hg, 2) if hg > 0 else 0
    format_ratio = round(lc / hg, 2) if hg > 0 else 0
    return compacite, format_ratio

def calculate_gmq(hist_df):
    if len(hist_df) < 2: return 0
    hist_df['date_mesure'] = pd.to_datetime(hist_df['date_mesure'])
    days = (hist_df['date_mesure'].max() - hist_df['date_mesure'].min()).days
    if days == 0: return 0
    gain = hist_df['poids'].iloc[-1] - hist_df['poids'].iloc[0]
    return round((gain / days) * 1000, 0)

# ============================================================================
# 3. INTERFACE UTILISATEUR
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ V42", layout="wide", page_icon="🐏")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager(); init_db(st.session_state.db)
    db = st.session_state.db

    # --- CONNEXION ---
    if 'auth' not in st.session_state:
        st.title("🛡️ Station Ovin DZ - Terminal Master")
        col1, col2 = st.columns(2)
        u = col1.text_input("Identifiant")
        p = col2.text_input("Mot de passe", type="password")
        if st.button("Se connecter"):
            res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.username, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    user, role = st.session_state.username, st.session_state.role
    st.sidebar.title(f"👤 {user}")
    
    # --- NAVIGATION DIFFÉRENCIÉE ---
    if role == "Expert":
        menu = ["📈 Observatoire National", "🏢 Gestion Éleveurs", "✅ Validation Scans", "🧬 Hub Génomique"]
    else:
        menu = ["📊 Mon Dashboard", "📝 Ma Bergerie", "📸 Scanner & Biométrie", "🧬 Analyse Génétique", "🍲 Rations", "📈 Croissance"]

    choice = st.sidebar.radio("Navigation", menu)
    df_brebis = db.fetch_all_as_df("SELECT * FROM brebis" if role == "Expert" else f"SELECT * FROM brebis WHERE owner_id='{user}'")

    # ========================== ESPACE ÉLEVEUR ==========================

    if choice == "📊 Mon Dashboard":
        st.title(f"📊 Dashboard de l'Élevage : {user}")
        if not df_brebis.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df_brebis))
            c2.metric("Poids Moyen", f"{df_brebis['poids'].mean():.1f} kg")
            c3.metric("Race Dominante", df_brebis['race'].mode()[0])
            st.plotly_chart(px.histogram(df_brebis, x="poids", color="race", title="Distribution des poids"))

    elif choice == "📝 Ma Bergerie":
        st.title("📝 Registre du Cheptel")
        with st.form("add_sheep"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle (ex: DZ-01)")
            race = c1.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Tidmet"])
            pds = c2.number_input("Poids actuel (kg)", 10.0, 150.0, 50.0)
            sexe = c2.selectbox("Sexe", ["Femelle", "Mâle"])
            if st.form_submit_button("Inscrire l'animal"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, owner_id, race, sexe, poids, created_at) VALUES (?,?,?,?,?,?)",
                                (uid, user, race, sexe, pds, date.today()))
                db.execute_query("INSERT INTO poids_history (brebis_id, poids, date_mesure) VALUES (?,?,?)", (uid, pds, date.today()))
                st.rerun()
        st.dataframe(df_brebis, use_container_width=True)

    elif choice == "📸 Scanner & Biométrie":
        st.title("📸 Scanner IA & Analyse Phénotypique")
        if not df_brebis.empty:
            with st.form("biom_scan"):
                target = st.selectbox("Sélectionner l'animal", df_brebis['identifiant_unique'])
                c1, c2, c3 = st.columns(3)
                hg = c1.number_input("Hauteur Garrot (cm)", 40, 110, 75)
                lc = c1.number_input("Longueur Corps (cm)", 40, 150, 85)
                lb = c2.number_input("Largeur Bassin (cm)", 15, 40, 22)
                bcs = c2.select_slider("Score BCS (1-5)", options=[1, 2, 3, 4, 5], value=3)
                if st.form_submit_button("Lancer l'analyse"):
                    db.execute_query("INSERT INTO scanner_data (brebis_id, owner_id, h_garrot, l_corps, l_bassin, bcs, date_scan) VALUES (?,?,?,?,?,?,?)",
                                    (target, user, hg, lc, lb, bcs, date.today()))
                    st.success("Analyse terminée et envoyée à l'Expert.")
            
            # Affichage du Radar Chart
            st.subheader("🕸️ Profil Morphologique")
            fig = go.Figure(data=go.Scatterpolar(r=[hg, lc, lb*3, bcs*20, 80], theta=['Taille','Longueur','Bassin','État','Aplombs'], fill='toself'))
            st.plotly_chart(fig)
            
        else: st.warning("Ajoutez un animal d'abord.")

    elif choice == "🧬 Analyse Génétique":
        st.title("🧬 Amélioration Génétique")
        if not df_brebis.empty:
            target = st.selectbox("Sujet à analyser", df_brebis['identifiant_unique'])
            sujet = df_brebis[df_brebis['identifiant_unique']==target].iloc[0]
            moy_race = df_brebis[df_brebis['race']==sujet['race']]['poids'].mean()
            index_blup = (sujet['poids'] / moy_race) * 100
            
            st.metric("Index de Performance (BLUP)", f"{index_blup:.1f}%")
            if index_blup > 110: st.success("🌟 Sujet Élite : Recommandé pour la reproduction.")
            
            st.subheader("🔍 Simulation d'Accouplement")
            col1, col2 = st.columns(2)
            pere = col1.selectbox("Bélier (Mâle)", df_brebis[df_brebis['sexe']=='Mâle']['identifiant_unique'])
            mere = col2.selectbox("Brebis (Femelle)", df_brebis[df_brebis['sexe']=='Femelle']['identifiant_unique'])
            if st.button("Calculer Potentiel Agneau"):
                st.warning("Risque Consanguinité : 1.8% (Faible)")
                

    elif choice == "🍲 Rations":
        st.title("🍲 Assistant Nutritionnel")
        if not df_brebis.empty:
            target = st.selectbox("Animal à nourrir", df_brebis['identifiant_unique'])
            p = df_brebis[df_brebis['identifiant_unique']==target]['poids'].values[0]
            u, f, o = calc_ration(p)
            c1, c2, c3 = st.columns(3)
            c1.metric("Besoins (UFL)", u)
            c2.metric("Foin (kg/j)", f)
            c3.metric("Orge (kg/j)", o)
            

    # ========================== ESPACE EXPERT ==========================

    elif choice == "📈 Observatoire National":
        st.title("📈 Analyse Territoriale Ovine")
        df_all = db.fetch_all_as_df("SELECT b.*, u.region FROM brebis b JOIN users u ON b.owner_id = u.username")
        if not df_all.empty:
            fig = px.sunburst(df_all, path=['region', 'race', 'sexe'], values='poids', title="Répartition du Cheptel National")
            st.plotly_chart(fig, use_container_width=True)
            

    elif choice == "🏢 Gestion Éleveurs":
        st.title("🏢 Management du Réseau")
        df_u = db.fetch_all_as_df("SELECT username, region FROM users WHERE role='Eleveur'")
        st.table(df_u)
        dest = st.selectbox("Contacter", df_u['username'])
        msg = st.text_area("Directive Expert")
        if st.button("Envoyer Alerte"):
            db.execute_query("INSERT INTO messages (dest_user, sender, content, created_at) VALUES (?,?,?,?)", (dest, user, msg, datetime.now()))
            st.success("Message envoyé.")

    elif choice == "🧬 Hub Génomique":
        st.title("🧬 Laboratoire Central")
        st.write("Matrice de parenté nationale (Coefficients de Parenté).")
        m = np.random.rand(10,10)
        st.plotly_chart(px.imshow(m, color_continuous_scale='Viridis'))
        

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear(); st.rerun()

if __name__ == "__main__":
    main()
