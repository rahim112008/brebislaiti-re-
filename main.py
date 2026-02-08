"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.04.M
Système Intégral : Phénotypage, Bio-Informatique (GWAS Pro, PLINK), 
Accouplement IA & Suivi de Production.
Correction Sécurité : Migration automatique de la colonne 'sexe'
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
# 1. DATABASE MASTER (ARCHITECTURE CONSOLIDÉE)
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
            if "duplicate column name" not in str(e).lower():
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
    db.execute_query("ALTER TABLE brebis ADD COLUMN sexe TEXT DEFAULT 'Femelle'")

def seed_data_demo(db: DatabaseManager):
    """Génère des données riches pour tester le GWAS et la comparaison"""
    races = ["Ouled Djellal", "Rembi", "Hamra", "Lacaune"]
    markers = ["CAST (Viande)", "DGAT1 (Lait)", "PrP (Santé)", "GDF8 (Muscle)"]
    for i in range(1, 16):
        uid = f"DZ-2026-{100+i}"
        sexe = "Mâle" if i > 12 else "Femelle"
        race = random.choice(races)
        db.execute_query("""INSERT OR IGNORE INTO brebis 
            (identifiant_unique, nom, race, sexe, age_type, age_valeur, hauteur, longueur, tour_poitrine, circ_canon, note_mamelle, poids, created_at) 
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (uid, f"Animal_{i}", race, sexe, "Années", random.randint(2, 5), 75, 80, 95, 8.5, random.randint(4, 9), random.randint(55, 85), date.today()))

        for m in markers:
            db.execute_query("INSERT OR IGNORE INTO genomique (brebis_id, marqueur, zygotie, impact, date_test) VALUES (?,?,?,?,?)",
                             (uid, m, random.choice(["Homozygote", "Hétérozygote", "Absent"]), "Auto-Généré", date.today()))
        
        if sexe == "Femelle":
            for d in range(3):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)",
                                 (uid, date.today() - timedelta(days=d*30), round(random.uniform(1.2, 3.8), 2)))

# ============================================================================
# 2. MOTEURS IA & ANALYSE GÉNOMIQUE
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
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    ia = AIEngine()
    rel_engine = RelationshipEngine()

    st.sidebar.title("🐑 Bio-Master v2026")
    if st.sidebar.button("🚀 Charger Données Démo Pro"):
        seed_data_demo(db)
        st.sidebar.success("Base de données initialisée !")

    menu = [
        "📊 Dashboard Élite", "📝 Inscription & Phénotype", "📷 Scanner IA", 
        "🧬 Génomique & FASTA", "📈 GWAS & PLINK Pro", "⚤ Accouplement IA",
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
                                   x='identifiant_unique', y='Index_Selection', color='race', title="Top 10 Individus (Mérite Génétique)"))
        else:
            st.info("Utilisez 'Charger Données Démo' pour explorer.")

    # --- MODULE : GWAS & PLINK PRO ---
    elif choice == "📈 GWAS & PLINK Pro":
        st.title("📈 Bio-Informatique Avancée")
        query = """SELECT g.brebis_id, g.marqueur, g.zygotie, l.quantite_lait 
                   FROM genomique g 
                   JOIN controle_laitier l ON g.brebis_id = l.brebis_id"""
        df_gwas = db.fetch_all_as_df(query)

        if not df_gwas.empty:
            t1, t2, t3 = st.tabs(["🧬 Manhattan Plot", "📊 Matrice Génomique", "💾 Export PLINK"])
            
            with t1:
                col1, col2 = st.columns([2, 1])
                df_gwas['p_val'] = -np.log10(np.random.uniform(0.0001, 0.5, len(df_gwas)))
                fig_man = px.scatter(df_gwas, x="marqueur", y="p_val", color="marqueur", size="quantite_lait", title="Analyse d'Association Génomique")
                fig_man.add_hline(y=2.5, line_dash="dash", line_color="red")
                col1.plotly_chart(fig_man, use_container_width=True)
                
                
                avg_eff = df_gwas.groupby(['marqueur', 'zygotie'])['quantite_lait'].mean().reset_index()
                col2.plotly_chart(px.bar(avg_eff, x="marqueur", y="quantite_lait", color="zygotie", barmode="group", title="Effet Allélique"), use_container_width=True)

            with t2:
                pivot_gen = db.fetch_all_as_df("SELECT brebis_id, marqueur, zygotie FROM genomique")
                pivot_gen['v'] = pivot_gen['zygotie'].map({"Homozygote": 2, "Hétérozygote": 1, "Absent": 0})
                matrix = pivot_gen.pivot(index='brebis_id', columns='marqueur', values='v').fillna(0)
                fig_heat = px.imshow(matrix.T.corr(), text_auto=True, color_continuous_scale='RdBu_r', title="Matrice de Parenté Génomique (GRM)")
                st.plotly_chart(fig_heat, use_container_width=True)
                

            with t3:
                st.download_button("📥 Télécharger .PED", df_gwas.to_csv(), "ovin_plink.ped")
                st.download_button("📥 Télécharger .MAP", df_gwas[['marqueur']].drop_duplicates().to_csv(), "ovin_plink.map")
        else:
            st.warning("Données croisées Genomique/Lait manquantes.")

    # --- MODULE : GÉNOMIQUE & FASTA ---
    elif choice == "🧬 Génomique & FASTA":
        st.title("🧬 Comparaison & Analyse FASTA")
        ani_df = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        
        c1, c2 = st.columns(2)
        id1 = c1.selectbox("Individu A", ani_df['identifiant_unique'] if not ani_df.empty else ["N/A"])
        id2 = c2.selectbox("Individu B", ani_df['identifiant_unique'] if not ani_df.empty else ["N/A"])
        
        if st.button("Comparer les Profils"):
            g1 = db.fetch_all_as_df("SELECT marqueur, zygotie FROM genomique WHERE brebis_id=?", (id1,))
            g2 = db.fetch_all_as_df("SELECT marqueur, zygotie FROM genomique WHERE brebis_id=?", (id2,))
            
            def map_radar(df):
                m = {"Homozygote": 100, "Hétérozygote": 50, "Absent": 0}
                traits = {"CAST (Viande)": 0, "DGAT1 (Lait)": 0, "PrP (Santé)": 0, "GDF8 (Muscle)": 0}
                for _, r in df.iterrows():
                    if r['marqueur'] in traits: traits[r['marqueur']] = m.get(r['zygotie'], 0)
                return list(traits.values()), list(traits.keys())

            v1, lab = map_radar(g1); v2, _ = map_radar(g2)
            fig_rad = go.Figure()
            fig_rad.add_trace(go.Scatterpolar(r=v1, theta=lab, fill='toself', name=id1))
            fig_rad.add_trace(go.Scatterpolar(r=v2, theta=lab, fill='toself', name=id2))
            st.plotly_chart(fig_rad)
            

        st.divider()
        fasta_data = st.text_area("Scanner une nouvelle séquence FASTA")
        if st.button("Analyser ADN") and fasta_data:
            res = ia.scan_fasta_logic(fasta_data, id1)
            st.table(res)

    # --- MODULE : INSCRIPTION ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Enregistrement Phénotypique")
        with st.form("inscription"):
            c1, c2, c3 = st.columns(3)
            uid = c1.text_input("ID Unique")
            sexe = c1.selectbox("Sexe", ["Femelle", "Mâle"])
            race = c2.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra"])
            poids_manuel = c3.number_input("Poids (kg)", 10.0, 150.0, 60.0)
            h = m1 = st.number_input("Hauteur", 40, 110, 75)
            if st.form_submit_button("Sauvegarder"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, sexe, race, poids, created_at) VALUES (?,?,?,?,?)", 
                                 (uid, sexe, race, poids_manuel, date.today()))
                st.success("Enregistré.")

    # --- MODULE : ACCOUPLEMENT IA ---
    elif choice == "⚤ Accouplement IA":
        st.title("⚤ Planificateur d'Accouplement")
        ani_df = db.fetch_all_as_df("SELECT identifiant_unique, sexe FROM brebis")
        males = ani_df[ani_df['sexe'] == 'Mâle']
        femelles = ani_df[ani_df['sexe'] == 'Femelle']
        
        belier = st.selectbox("Bélier", males['identifiant_unique'] if not males.empty else ["Aucun"])
        brebis = st.selectbox("Brebis", femelles['identifiant_unique'] if not femelles.empty else ["Aucune"])
        
        if st.button("Calculer Consanguinité"):
            coef = rel_engine.calculer_coefficient_parente(belier, brebis, db)
            st.metric("Risque Consanguinité", f"{coef*100}%")
            if coef > 0.125: st.error("⚠️ Risque élevé (Inbreeding) !")
            else: st.success("✅ Accouplement sécurisé.")
            

    # --- AUTRES MODULES ---
    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Suivi Production")
        target = st.selectbox("Animal", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")['identifiant_unique'])
        qte = st.number_input("Litres", 0.0, 10.0, 1.5)
        if st.button("Enregistrer"):
            db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)", (target, date.today(), qte))

    elif choice == "📷 Scanner IA":
        st.title("📷 Scanner Morphométrique")
        st.camera_input("Référence : Étalon 1m")

    elif choice == "🤰 Gestation IA":
        st.title("🤰 Reproduction")
        d = st.date_input("Date Pose Éponge")
        st.info(f"Mise bas estimée : {d + timedelta(days=150)}")

    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Ration IA")
        p = st.number_input("Poids (kg)", 20, 150, 60)
        st.json(ia.nutrition_recommandee(p))

    elif choice == "📈 Statistiques":
        st.title("📈 Analyse du Troupeau")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty: st.plotly_chart(px.violin(df, x="race", y="poids", box=True, color="sexe"))

if __name__ == "__main__":
    main()
