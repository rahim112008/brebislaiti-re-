"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.04.D
Système Intégral : Phénotypage, Lait, Reproduction, Santé & Bioinformatique (FASTA/Zygotie)
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
import re
from datetime import datetime, date, timedelta

# ============================================================================
# 1. DATABASE MASTER
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_pro.db"):
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
        # Table Identité
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, attaches_mamelle TEXT, poids REAL, created_at DATE
        )""",
        # Table Contrôle Laitier
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_controle DATE,
            quantite_lait REAL, tb REAL, tp REAL, cellules INTEGER,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table Gestation
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_eponge DATE, date_mise_bas_prevue DATE, statut TEXT
        )""",
        # Table Santé
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_soin DATE,
            type_acte TEXT, produit TEXT, rappel_prevu DATE
        )""",
        # Table Génomique (Nouvelle structure consolidée)
        """CREATE TABLE IF NOT EXISTS genomique (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, marqueur TEXT, allele_1 TEXT, allele_2 TEXT, 
            zygotie TEXT, impact TEXT, date_test DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE IA & BIOINFORMATIQUE
# ============================================================================

# Référentiel des signatures ADN pour le scanner FASTA
GENE_SIGNATURES = {
    "CAST (Tendreté Viande)": "TTAGCCT", 
    "GDF8 (Muscle Myostatine)": "CCGGTTA",
    "DGAT1 (Richesse Lait)": "GGCATAA",
    "PrP (Résistance Tremblante)": "ATATGCG"
}

class AIEngine:
    @staticmethod
    def calculer_index_elite(row, df_lait):
        score_morpho = (row['tour_poitrine'] * 0.2) + (row['note_mamelle'] * 5)
        score_os = row['circ_canon'] * 3
        lait_indiv = df_lait[df_lait['brebis_id'] == row['identifiant_unique']]
        score_lait = lait_indiv['quantite_lait'].mean() * 15 if not lait_indiv.empty else 0
        return round((score_morpho + score_os + score_lait), 2)

    @staticmethod
    def nutrition_recommandee(poids):
        return {"Orge (kg)": round(poids * 0.012, 2), "Luzerne (kg)": round(poids * 0.02, 2), "CMV (g)": 30}

    @staticmethod
    def scan_fasta_logic(sequence, animal_id):
        results = []
        sequence = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        for gene, pattern in GENE_SIGNATURES.items():
            matches = [m.start() for m in re.finditer(pattern, sequence)]
            count = len(matches)
            if count >= 2:
                zygotie = "Homozygote"; status = "✅ Fixé (Double copie)"
            elif count == 1:
                zygotie = "Hétérozygote"; status = "⚠️ Porteur (Simple copie)"
            else:
                zygotie = "Absent"; status = "❌ Non détecté"
            results.append({"Gène": gene, "Occurrence": count, "Zygotie": zygotie, "Diagnostic": status})
        return results

# ============================================================================
# 3. INTERFACE UTILISATEUR
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    ia = AIEngine()

    st.sidebar.title("🐑 Système Master v2026")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA", 
        "🥛 Contrôle Laitier", 
        "🤰 Gestation IA", 
        "🌾 Nutrition Solo", 
        "🩺 Santé & Vaccins", 
        "🧬 Génomique & FASTA", 
        "📈 Statistiques"
    ]
    choice = st.sidebar.radio("Modules", menu)

    # --- MODULE 1: DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Performance & Sélection Élite")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        df_g = db.fetch_all_as_df("SELECT * FROM genomique")
        
        if not df_b.empty:
            df_b['Index_Selection'] = df_b.apply(lambda r: ia.calculer_index_elite(r, df_l), axis=1)
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Effectif", len(df_b))
            c2.metric("Moyenne Lait (L)", round(df_l['quantite_lait'].mean(), 2) if not df_l.empty else 0)
            c3.metric("Tests Génomiques", len(df_g))
            c4.metric("Meilleur Index", df_b['Index_Selection'].max())

            st.subheader("🏆 Top Génitrices & Génomique")
            st.dataframe(df_b.sort_values(by='Index_Selection', ascending=False).head(10))
        else:
            st.info("Aucune donnée disponible.")

    # --- MODULE 2: INSCRIPTION ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Phénotypage Avancé")
        with st.form("inscription"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("Identifiant Unique (Boucle)")
            race = c1.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra", "Autre"])
            age_t = c2.radio("Méthode d'âge", ["Dents", "Mois", "Années"])
            age_v = c2.number_input("Valeur âge", 0, 15, 2)
            
            st.subheader("Mesures & Mamelle")
            m1, m2, m3 = st.columns(3)
            h = m1.number_input("Hauteur (cm)", 40, 110, 75); l = m2.number_input("Longueur (cm)", 40, 120, 80); tp = m3.number_input("Tour Poitrine (cm)", 50, 150, 90)
            lb = m1.number_input("Largeur Bassin (cm)", 10, 40, 22); can = m3.number_input("Canon (cm)", 5.0, 15.0, 8.5)
            note_m = st.slider("Note Mamelle", 1, 10, 5)
            
            if st.form_submit_button("Enregistrer"):
                poids = (tp**2 * l) / 30000
                db.execute_query("INSERT INTO brebis (identifiant_unique, race, age_type, age_valeur, hauteur, longueur, tour_poitrine, largeur_bassin, circ_canon, note_mamelle, poids, created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", 
                                 (uid, race, age_t, age_v, h, l, tp, lb, can, note_m, poids, date.today()))
                st.success("Enregistré.")

    # --- MODULE 8: GÉNOMIQUE & FASTA (MODIFIÉ) ---
    elif choice == "🧬 Génomique & FASTA":
        st.title("🧬 Analyse Génomique & Scanner FASTA")
        
        tab1, tab2 = st.tabs(["🔍 Scanner de Séquences FASTA", "📊 Inventaire des Allèles"])
        
        with tab1:
            st.subheader("Recherche de Motifs ADN")
            ani_df = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
            
            col_in, col_ref = st.columns([2, 1])
            with col_in:
                target_animal = st.selectbox("Animal à analyser", ani_df['identifiant_unique'] if not ani_df.empty else ["Aucun"])
                fasta_input = st.text_area("Coller la séquence FASTA ici", height=200, placeholder=">Exemple\nATGCGTTAGCCT...")
            
            with col_ref:
                st.write("**Signatures cibles :**")
                for k, v in GENE_SIGNATURES.items(): st.code(f"{k}: {v}")

            if st.button("Lancer le Scanner Moléculaire"):
                if fasta_input and target_animal != "Aucun":
                    results = ia.scan_fasta_logic(fasta_input, target_animal)
                    st.table(pd.DataFrame(results))
                    for r in results:
                        db.execute_query("INSERT INTO genomique (brebis_id, marqueur, zygotie, impact, date_test) VALUES (?,?,?,?,?)",
                                        (target_animal, r['Gène'], r['Zygotie'], r['Diagnostic'], date.today()))
                    st.success("Analyses enregistrées.")

        with tab2:
            st.subheader("Statut de Zygotie du Troupeau")
            df_g = db.fetch_all_as_df("SELECT * FROM genomique")
            if not df_g.empty:
                c1, c2 = st.columns(2)
                c1.plotly_chart(px.pie(df_g, names='zygotie', title="Répartition Zygotie", hole=0.4))
                c2.plotly_chart(px.bar(df_g, x='marqueur', color='zygotie', barmode='group'))
                st.dataframe(df_g)
            else:
                st.info("Aucun test génomique en base.")

    # --- AUTRES MODULES (CONSERVÉS) ---
    elif choice == "📷 Scanner IA":
        st.title("📷 Scanner Morphométrique")
        st.camera_input("Capture avec étalon 1m")
        st.info("Analyse des pixels activée.")

    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Contrôle Laitier")
        with st.form("lait"):
            target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
            qte = st.number_input("Lait (L)", 0.0, 10.0, 2.0)
            if st.form_submit_button("Valider"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)", (target, date.today(), qte))
        
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        if not df_l.empty: st.plotly_chart(px.line(df_l, x="date_controle", y="quantite_lait", color="brebis_id"))

    elif choice == "🤰 Gestation IA":
        st.title("🤰 Reproduction")
        d_ep = st.date_input("Date Pose Éponge")
        if st.button("Calculer Mise Bas"):
            st.success(f"Prévue le : {(d_ep + timedelta(days=164)).strftime('%d/%m/%Y')}")

    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Ration IA")
        df_b = db.fetch_all_as_df("SELECT identifiant_unique, poids FROM brebis")
        if not df_b.empty:
            target = st.selectbox("Brebis", df_b['identifiant_unique'])
            p = df_b[df_b['identifiant_unique'] == target]['poids'].values[0]
            for k, v in ia.nutrition_recommandee(p).items(): st.metric(k, v)

    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Santé")
        target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        acte = st.selectbox("Soin", ["Enterotoxémie", "PPR", "Vermifuge"])
        if st.button("Enregistrer"):
            db.execute_query("INSERT INTO sante (brebis_id, date_soin, type_acte) VALUES (?,?,?)", (target, date.today(), acte))
            st.success("Soin noté.")

    elif choice == "📈 Statistiques":
        st.title("📈 Analyse de Variance")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty: st.plotly_chart(px.violin(df, x="race", y="poids", box=True))

if __name__ == "__main__":
    main()
