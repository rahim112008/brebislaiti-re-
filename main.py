"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.04.I
Système Intégral : Phénotypage, Lait, GWAS, PLINK & Accouplement IA
Inclus : Générateur de Base de Données de Démonstration
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
import random

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
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, sexe TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, attaches_mamelle TEXT, poids REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_controle DATE,
            quantite_lait REAL, tb REAL, tp REAL, cellules INTEGER,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS genomique (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, marqueur TEXT, zygotie TEXT, impact TEXT, date_test DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_eponge DATE, date_mise_bas_prevue DATE, statut TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

def seed_data_demo(db: DatabaseManager):
    """Génère des données fictives pour tester toutes les fonctionnalités"""
    races = ["Ouled Djellal", "Rembi", "Hamra", "Lacaune"]
    # Création de 12 animaux (10 femelles, 2 mâles)
    for i in range(1, 13):
        uid = f"DZ-2026-{100+i}"
        sexe = "Mâle" if i > 10 else "Femelle"
        race = random.choice(races)
        db.execute_query("""INSERT OR IGNORE INTO brebis 
            (identifiant_unique, nom, race, sexe, age_type, age_valeur, hauteur, longueur, tour_poitrine, circ_canon, note_mamelle, poids, created_at) 
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (uid, f"Animal_{i}", race, sexe, "Années", random.randint(2, 5), 75, 80, 95, 8.5, random.randint(4, 9), random.randint(55, 85), date.today()))

        # Ajout de données génomiques pour chaque animal
        markers = [("CAST (Viande)", "Homozygote"), ("DGAT1 (Lait)", "Hétérozygote")]
        for m, z in markers:
            db.execute_query("INSERT OR IGNORE INTO genomique (brebis_id, marqueur, zygotie, impact, date_test) VALUES (?,?,?,?,?)",
                             (uid, m, z, "Analyse Démo", date.today()))

    # Ajout de contrôles laitiers pour les femelles
    for i in range(1, 11):
        uid = f"DZ-2026-{100+i}"
        for d in range(5):
            db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)",
                             (uid, date.today() - timedelta(days=d*7), round(random.uniform(1.0, 3.5), 2)))

# ============================================================================
# 2. MOTEURS IA & LOGIQUE
# ============================================================================

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
            zygotie = "Homozygote" if count >= 2 else "Hétérozygote" if count == 1 else "Absent"
            status = "✅ Fixé" if count >= 2 else "⚠️ Porteur" if count == 1 else "❌ Non détecté"
            results.append({"Gène": gene, "Occurrence": count, "Zygotie": zygotie, "Diagnostic": status})
        return results

class RelationshipEngine:
    @staticmethod
    def calculer_coefficient_parente(id1, id2, db):
        g1 = db.fetch_all_as_df("SELECT marqueur, zygotie FROM genomique WHERE brebis_id=?", (id1,))
        g2 = db.fetch_all_as_df("SELECT marqueur, zygotie FROM genomique WHERE brebis_id=?", (id2,))
        if g1.empty or g2.empty: return 0.0
        communs = pd.merge(g1, g2, on="marqueur")
        if communs.empty: return 0.0
        matches = sum(communs['zygotie_x'] == communs['zygotie_y'])
        return round((matches / len(communs)) * 0.5, 3)

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
    rel_engine = RelationshipEngine()

    # --- SIDEBAR & DEMO ---
    st.sidebar.title("🐑 Système Master v2026")
    if st.sidebar.button("🚀 Charger Données Démo"):
        seed_data_demo(db)
        st.sidebar.success("Base de données initialisée !")

    menu = [
        "📊 Dashboard Élite", "📝 Inscription & Phénotype", "📷 Scanner IA", 
        "🧬 Génomique & FASTA", "📈 GWAS & PLINK", "⚤ Accouplement IA",
        "🥛 Contrôle Laitier", "🤰 Gestation IA", "🌾 Nutrition Solo", 
        "🩺 Santé & Vaccins", "📈 Statistiques"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE : DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Statut Global & Sélection")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        if not df_b.empty:
            df_b['Index_Selection'] = df_b.apply(lambda r: ia.calculer_index_elite(r, df_l), axis=1)
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df_b))
            c2.metric("Moyenne Lait (L)", round(df_l['quantite_lait'].mean(), 2) if not df_l.empty else 0)
            c3.metric("Meilleur Index", df_b['Index_Selection'].max())
            st.plotly_chart(px.bar(df_b.sort_values('Index_Selection', ascending=False).head(10), 
                                   x='identifiant_unique', y='Index_Selection', color='race', title="Top 10 Individus"))
        else:
            st.info("Base vide. Utilisez le bouton 'Charger Données Démo' pour tester.")

    # --- MODULE : INSCRIPTION ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Enregistrement Phénotypique")
        with st.form("inscription"):
            c1, c2, c3 = st.columns(3)
            uid = c1.text_input("ID Unique")
            sexe = c1.selectbox("Sexe", ["Femelle", "Mâle"])
            race = c2.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra"])
            poids_manuel = c3.number_input("Poids (kg)", 10.0, 150.0, 60.0)
            
            st.subheader("Morphométrie (cm)")
            m1, m2, m3 = st.columns(3)
            h = m1.number_input("Hauteur", 40, 110, 75); l = m2.number_input("Longueur", 40, 120, 80); tp = m3.number_input("Tour Poitrine", 50, 150, 90)
            can = m1.number_input("Canon", 5.0, 15.0, 8.5); note_m = m2.slider("Note Mamelle", 1, 10, 5)

            if st.form_submit_button("Sauvegarder"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, sexe, race, hauteur, longueur, tour_poitrine, circ_canon, note_mamelle, poids, created_at) VALUES (?,?,?,?,?,?,?,?,?,?)", 
                                 (uid, sexe, race, h, l, tp, can, note_m, poids_manuel, date.today()))
                st.success("Enregistré.")

    # --- MODULE : GÉNOMIQUE & FASTA ---
    elif choice == "🧬 Génomique & FASTA":
        st.title("🧬 Scanner Moléculaire FASTA")
        ani_df = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        target = st.selectbox("Assigner à :", ani_df['identifiant_unique'] if not ani_df.empty else ["Aucun"])
        fasta_data = st.text_area("Coller Séquence FASTA", height=150)
        
        if st.button("Analyser ADN") and fasta_data and target != "Aucun":
            results = ia.scan_fasta_logic(fasta_data, target)
            st.table(pd.DataFrame(results))
            for r in results:
                db.execute_query("INSERT INTO genomique (brebis_id, marqueur, zygotie, impact, date_test) VALUES (?,?,?,?,?)",
                                 (target, r['Gène'], r['Zygotie'], r['Diagnostic'], date.today()))

    # --- MODULE : GWAS & PLINK ---
    elif choice == "📈 GWAS & PLINK":
        st.title("📈 GWAS (Genome-Wide Association Study)")
        query = "SELECT g.marqueur, g.zygotie, l.quantite_lait FROM genomique g JOIN controle_laitier l ON g.brebis_id = l.brebis_id"
        df_gwas = db.fetch_all_as_df(query)
        if not df_gwas.empty:
            df_plot = df_gwas.groupby(['marqueur', 'zygotie'])['quantite_lait'].mean().reset_index()
            df_plot['significativite'] = np.random.uniform(1, 5, len(df_plot)) 
            st.plotly_chart(px.scatter(df_plot, x="marqueur", y="significativite", size="quantite_lait", color="zygotie", title="Manhattan Plot (Simulation)"))
            
            st.download_button("Exporter pour PLINK (.PED)", df_gwas.to_csv(), "ovin_plink_export.ped")
        else:
            st.warning("Données croisées (ADN + Lait) insuffisantes.")

    # --- MODULE : ACCOUPLEMENT IA ---
    elif choice == "⚤ Accouplement IA":
        st.title("⚤ Planificateur d'Accouplement")
        ani_df = db.fetch_all_as_df("SELECT identifiant_unique, sexe FROM brebis")
        males = ani_df[ani_df['sexe'] == 'Mâle']
        femelles = ani_df[ani_df['sexe'] == 'Femelle']
        
        c1, c2 = st.columns(2)
        belier = c1.selectbox("Bélier", males['identifiant_unique'] if not males.empty else ["Aucun"])
        brebis = c2.selectbox("Brebis", femelles['identifiant_unique'] if not femelles.empty else ["Aucune"])
        
        if st.button("Calculer Consanguinité"):
            if belier != "Aucun" and brebis != "Aucune":
                coef = rel_engine.calculer_coefficient_parente(belier, brebis, db)
                fig = go.Figure(go.Indicator(mode="gauge+number", value=coef*100, title={'text': "Risque %"},
                                            gauge={'axis': {'range': [0, 50]}, 'steps': [{'range': [0, 6.25], 'color': "green"}, {'range': [12.5, 50], 'color': "red"}]}))
                st.plotly_chart(fig)
                
                if coef > 0.125: st.error("☢️ Risque élevé (Inbreeding) !")
                else: st.success("✅ Accouplement sécurisé.")

    # --- MODULES DE GESTION ---
    elif choice == "📷 Scanner IA":
        st.title("📷 Scanner Morphométrique")
        st.camera_input("Étallon 1m obligatoire")

    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Suivi Production")
        with st.form("lait"):
            target = st.selectbox("Animal", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")['identifiant_unique'])
            qte = st.number_input("Litres", 0.0, 10.0, 2.0)
            if st.form_submit_button("Valider"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)", (target, date.today(), qte))

    elif choice == "🤰 Gestation IA":
        st.title("🤰 Suivi Reproduction")
        d_ep = st.date_input("Date Pose Éponge")
        st.info(f"Mise bas estimée : {d_ep + timedelta(days=164)}")

    elif choice == "📈 Statistiques":
        st.title("📈 Analyse du Troupeau")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty: st.plotly_chart(px.violin(df, x="race", y="poids", box=True, points="all"))

if __name__ == "__main__":
    main()
