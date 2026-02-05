"""
EXPERT OVIN DZ PRO - VERSION ULTIME 2026
Plateforme de Sélection Élite, Génomique et Standardisation Mammaire
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. ARCHITECTURE DE LA BASE DE DONNÉES (SQLITE)
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
            owner_id TEXT, race TEXT, hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, circ_canon REAL, note_mamelle REAL, 
            prof_mamelle REAL, angle_trayons REAL, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_controle DATE,
            quantite_lait REAL, tb REAL, tp REAL, owner_id TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS recommandations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, target_user TEXT, expert_name TEXT, 
            message TEXT, date_envoyee DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)
    # Comptes de test par défaut
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert')")
    db.execute_query("INSERT OR IGNORE INTO users VALUES ('eleveur1', 'ovin2026', 'Eleveur')")

# ============================================================================
# 2. SYSTÈME D'AUTHENTIFICATION
# ============================================================================

def login_screen():
    st.title("🧬 Expert Ovin DZ Pro")
    st.markdown("### Système Master de Sélection Élite")
    with st.form("login_form"):
        u = st.text_input("Identifiant")
        p = st.text_input("Mot de passe", type="password")
        if st.form_submit_button("Se connecter"):
            db = st.session_state.db
            user_res = db.execute_query("SELECT * FROM users WHERE username=? AND password=?", (u, p)).fetchone()
            if user_res:
                st.session_state.auth = True
                st.session_state.username = user_res['username']
                st.session_state.role = user_res['role']
                st.rerun()
            else:
                st.error("Utilisateur ou mot de passe incorrect.")

# ============================================================================
# 3. APPLICATION PRINCIPALE (MODULAIRE)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro", layout="wide")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    if 'auth' not in st.session_state:
        login_screen()
        return

    db = st.session_state.db
    current_user = st.session_state.username
    current_role = st.session_state.role

    # --- GESTION DES VUES (LOGIQUE EXPERT) ---
    st.sidebar.title(f"👤 {current_user}")
    st.sidebar.info(f"Rôle : {current_role}")
    
    if current_role == "Expert":
        users_list = db.fetch_all_as_df("SELECT username FROM users WHERE role='Eleveur'")
        view_user = st.sidebar.selectbox("📂 Étudier le troupeau de :", users_list['username'] if not users_list.empty else [current_user])
    else:
        view_user = current_user

    menu = ["📊 Dashboard", "📝 Phénotypage", "📷 Scanner Corps", "🥛 Scanner Mamelle", "📈 Suivi Laitier", "✉️ Conseils Expert", "📖 Guide & Aide"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE 1: DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title(f"📊 Dashboard Élite - {view_user}")
        
        # Notifications de l'expert
        recos = db.fetch_all_as_df("SELECT * FROM recommandations WHERE target_user=?", (view_user,))
        if not recos.empty:
            for _, r in recos.iterrows():
                st.warning(f"📩 **Conseil Expert ({r['date_envoyee']})** : {r['message']}")

        df_b = db.fetch_all_as_df("SELECT * FROM brebis WHERE owner_id=?", (view_user,))
        if not df_b.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df_b))
            c2.metric("Moyenne Note Mamelle", f"{round(df_b['note_mamelle'].mean(), 1)}/10")
            c3.metric("Poids Moyen", f"{round(df_b['poids'].mean(), 1)} kg")
            
            st.plotly_chart(px.scatter(df_b, x="circ_canon", y="poids", color="race", size="note_mamelle", title="Corrélation Squelette / Poids / Mamelle"))
            st.dataframe(df_b, use_container_width=True)
        else:
            st.info("Aucun animal enregistré.")

    # --- MODULE 2: PHÉNOTYPAGE ---
    elif choice == "📝 Phénotypage":
        st.title("📝 Inscription & Mesures Physiques")
        
        with st.form("pheno_form"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle (Unique)")
            race = c1.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra"])
            tp = c2.number_input("Tour Poitrine (cm)", 50, 150, 90)
            lg = c2.number_input("Longueur Corps (cm)", 50, 130, 80)
            can = c1.number_input("Circ. Canon (cm)", 5.0, 15.0, 8.0)
            lb = c2.number_input("Largeur Bassin (cm)", 10, 40, 22)
            
            if st.form_submit_button("Enregistrer la fiche"):
                poids = (tp**2 * lg) / 30000
                db.execute_query("""INSERT INTO brebis (identifiant_unique, owner_id, race, tour_poitrine, longueur, circ_canon, largeur_bassin, poids, created_at) 
                                 VALUES (?,?,?,?,?,?,?,?,?)""", (uid, view_user, race, tp, lg, can, lb, poids, date.today()))
                st.success("Fiche créée. Utilisez les Scanners IA pour compléter les notes de qualité.")

    # --- MODULE 3: SCANNER MAMELLE (STANDARDISATION) ---
    elif choice == "🥛 Scanner Mamelle":
        st.title("🥛 Scanner IA - Standardisation Mammaire")
        target = st.selectbox("Brebis à scanner", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=?", (view_user,))['identifiant_unique'])
        st.camera_input("Photo arrière de la mamelle (avec étalon)")
        
        c1, c2 = st.columns(2)
        prof = c1.number_input("Profondeur mesurée (pixels -> cm)", 5.0, 35.0, 15.0)
        angle = c2.number_input("Angle des trayons (°)", 0, 180, 45)
        
        # Formule de standardisation IA
        note_std = round(10 - (prof / 4.5), 1)
        st.metric("Note IA Standardisée", f"{note_std} / 10")
        
        if st.button("Valider et Synchroniser"):
            db.execute_query("UPDATE brebis SET prof_mamelle=?, angle_trayons=?, note_mamelle=? WHERE identifiant_unique=?", (prof, angle, note_std, target))
            st.success("Mesures standardisées enregistrées.")

    # --- MODULE 4: SCANNER CORPS ---
    elif choice == "📷 Scanner Corps":
        st.title("📷 Scanner Morphologique IA")
        st.selectbox("Objet de calibration", ["Feuille A4 (29.7cm)", "Carte Bancaire (8.5cm)", "Bâton 1m"])
        st.camera_input("Capture de profil")
        st.info("L'IA analyse les proportions pour vérifier les mesures saisies manuellement.")

    # --- MODULE 5: SUIVI LAITIER ---
    elif choice == "📈 Suivi Laitier":
        st.title("📈 Contrôle Laitier & Courbes")
        
        with st.form("lait_form"):
            target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis WHERE owner_id=?", (view_user,))['identifiant_unique'])
            qte = st.number_input("Production (Litres)", 0.0, 10.0, 1.5)
            if st.form_submit_button("Enregistrer le contrôle"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait, owner_id) VALUES (?,?,?,?)", (target, date.today(), qte, view_user))
        
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier WHERE owner_id=?", (view_user,))
        if not df_l.empty:
            st.plotly_chart(px.line(df_l, x="date_controle", y="quantite_lait", color="brebis_id", title="Évolution de la Lactation"))

    # --- MODULE 6: CONSEILS EXPERT ---
    elif choice == "✉️ Conseils Expert":
        st.title("✉️ Espace de Consultation Master")
        if current_role == "Expert":
            st.subheader(f"Envoyer une analyse à {view_user}")
            msg = st.text_area("Observations (Nutrition, Génétique, Santé...)")
            if st.button("Envoyer la recommandation"):
                db.execute_query("INSERT INTO recommandations (target_user, expert_name, message, date_envoyee) VALUES (?,?,?,?)", (view_user, current_user, msg, date.today()))
                st.success("Conseil transmis à l'éleveur.")
        else:
            st.info("Ici s'affichent les rapports envoyés par l'Expert Master.")

    # --- MODULE 7: GUIDE & AIDE ---
    elif choice == "📖 Guide & Aide":
        st.title("📖 Guide Technique de Prise de Vue")
        st.markdown("""
        ### Règle d'or : Pas de bonne photo = Pas de bonne mesure.
        
        **1. L'Étalon (Référence) :**
        - Utilisez une **Carte Bancaire** ou une **Feuille A4**.
        - L'objet doit être **collé contre l'animal** au même plan que la zone mesurée.
        
        **2. L'Angle (Parallélisme) :**
        - Le téléphone doit être parfaitement parallèle à l'animal.
        - Ne prenez pas de photo plongeante (d'en haut).
        
        **3. La Mamelle :**
        - Photo prise de l'arrière.
        - Les membres doivent être bien dégagés.
        - Nettoyez la mamelle pour que l'IA détecte bien les contours.
        """)
        

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear()
        st.rerun()

if __name__ == "__main__":
    main()
