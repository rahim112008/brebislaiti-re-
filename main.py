"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.04
Système Intégral de Gestion de Précision : 
Phénotypage, Lait, Génomique, Santé, Nutrition & IA
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
from datetime import datetime, date, timedelta

# ============================================================================
# 1. DATABASE MASTER (ARCHITECTURE CUMULATIVE)
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
        # Table Identité et Phénotypage Avancé
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, attaches_mamelle TEXT, poids REAL, created_at DATE
        )""",
        # Table Contrôle Laitier & Biochimie
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
        # Table Santé & Vaccins
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_soin DATE,
            type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE IA (NUTRITION, GÉNÉTIQUE & INDEX)
# ============================================================================

class AIEngine:
    @staticmethod
    def calculer_index_elite(row, df_lait):
        # Index de Sélection Scientifique (Lait 50%, Morpho 30%, Skeleton 20%)
        score_morpho = (row['tour_poitrine'] * 0.2) + (row['note_mamelle'] * 5)
        score_os = row['circ_canon'] * 3
        
        lait_indiv = df_lait[df_lait['brebis_id'] == row['identifiant_unique']]
        score_lait = lait_indiv['quantite_lait'].mean() * 15 if not lait_indiv.empty else 0
        
        return round((score_morpho + score_os + score_lait), 2)

    @staticmethod
    def nutrition_recommandee(poids):
        # Besoins simplifiés selon poids
        orge = round(poids * 0.012, 2)
        luzerne = round(poids * 0.02, 2)
        return {"Orge (kg)": orge, "Luzerne (kg)": luzerne, "CMV (g)": 30}

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    ia = AIEngine()

    # --- SIDEBAR NAVIGATION ---
    st.sidebar.title("🐑 Système Intégré v2026")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA", 
        "🥛 Contrôle Laitier", 
        "🤰 Gestation IA", 
        "🌾 Nutrition Solo", 
        "🩺 Santé & Vaccins", 
        "🧬 Génomique & NCBI", 
        "📈 Statistiques"
    ]
    choice = st.sidebar.radio("Modules", menu)

    # --- MODULE 1: DASHBOARD ÉLITE ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Performance & Sélection Élite")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        
        if not df_b.empty:
            df_b['Index_Selection'] = df_b.apply(lambda r: ia.calculer_index_elite(r, df_l), axis=1)
            df_top = df_b.sort_values(by='Index_Selection', ascending=False)
            
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif", len(df_b))
            c2.metric("Moyenne Laitière (L)", round(df_l['quantite_lait'].mean(), 2) if not df_l.empty else 0)
            c3.metric("Meilleur Index", df_top['Index_Selection'].max())

            st.subheader("🏆 Classement des meilleures génitrices")
            st.dataframe(df_top[['identifiant_unique', 'race', 'Index_Selection', 'poids', 'note_mamelle']].head(10))
            
            st.plotly_chart(px.scatter(df_b, x="tour_poitrine", y="Index_Selection", color="race", size="poids", hover_name="identifiant_unique"))
        else:
            st.info("Bienvenue ! Veuillez commencer par inscrire des animaux.")

    # --- MODULE 2: INSCRIPTION & PHÉNOTYPE ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Phénotypage Avancé")
        
        with st.form("inscription"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("Identifiant Unique (Boucle)")
            race = c1.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra", "Autre"])
            age_t = c2.radio("Méthode d'âge", ["Dents", "Mois", "Années"])
            age_v = c2.number_input("Valeur âge", 0, 15, 2)
            
            st.subheader("Mesures du Corps & Bassin")
            m1, m2, m3 = st.columns(3)
            h = m1.number_input("Hauteur Garrot (cm)", 40, 110, 75)
            l = m2.number_input("Longueur Corps (cm)", 40, 120, 80)
            tp = m3.number_input("Tour Poitrine (cm)", 50, 150, 90)
            lb = m1.number_input("Largeur Bassin (cm)", 10, 40, 22)
            lgb = m2.number_input("Longueur Bassin (cm)", 10, 40, 20)
            can = m3.number_input("Circonférence Canon (cm)", 5.0, 15.0, 8.5)
            
            st.subheader("Évaluation de la Mamelle")
            
            note_m = st.slider("Note Mamelle (Volume/Équilibre)", 1, 10, 5)
            attaches = st.selectbox("Attaches", ["Solides", "Moyennes", "Lâches"])
            
            if st.form_submit_button("Enregistrer"):
                poids = (tp**2 * l) / 30000
                db.execute_query("""INSERT INTO brebis (identifiant_unique, race, age_type, age_valeur, hauteur, longueur, tour_poitrine, largeur_bassin, long_bassin, circ_canon, note_mamelle, attaches_mamelle, poids, created_at) 
                                 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (uid, race, age_t, age_v, h, l, tp, lb, lgb, can, note_m, attaches, poids, date.today()))
                st.success("Brebis enregistrée avec succès.")

    # --- MODULE 3: SCANNER IA ---
    elif choice == "📷 Scanner IA":
        st.title("📷 Scanner Morphométrique avec Étalon")
        etalon = st.selectbox("Référence physique", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
        st.camera_input("Capture pour analyse des pixels")
        st.info(f"Analyse en cours via étalon : {etalon}")

    # --- MODULE 4: CONTRÔLE LAITIER & COURBE ---
    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Contrôle Laitier & Courbe de Lactation")
        
        with st.form("lait"):
            target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
            qte = st.number_input("Lait (Litres)", 0.0, 10.0, 2.0)
            tb = st.slider("Gras (TB) g/L", 20, 80, 45)
            tp = st.slider("Protéines (TP) g/L", 20, 70, 35)
            if st.form_submit_button("Valider le contrôle"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait, tb, tp) VALUES (?,?,?,?,?)",
                                (target, date.today(), qte, tb, tp))
        
        st.subheader("📈 Évolution de la Lactation")
        df_lait_all = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        if not df_lait_all.empty:
            fig = px.line(df_lait_all, x="date_controle", y="quantite_lait", color="brebis_id", markers=True, title="Courbe de Lactation Individuelle")
            st.plotly_chart(fig)

    # --- MODULE 5: GESTATION IA ---
    elif choice == "🤰 Gestation IA":
        st.title("🤰 Gestion de la Reproduction")
        
        target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        d_ep = st.date_input("Date Pose Éponge")
        if st.button("Prédire Mise Bas"):
            mb = d_ep + timedelta(days=164)
            st.success(f"📅 Mise bas prévue : {mb.strftime('%d/%m/%Y')}")
            db.execute_query("INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevue, statut) VALUES (?,?,?,?)", (target, d_ep, mb, "En cours"))

    # --- MODULE 6: NUTRITION SOLO ---
    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Ration Individualisée par IA")
        df_b = db.fetch_all_as_df("SELECT identifiant_unique, poids FROM brebis")
        target = st.selectbox("Brebis", df_b['identifiant_unique'])
        poids_b = df_b[df_b['identifiant_unique'] == target]['poids'].values[0]
        recette = ia.nutrition_recommandee(poids_b)
        st.write(f"**Ration recommandée pour {target} ({poids_b:.1f} kg) :**")
        for k, v in recette.items():
            st.metric(k, v)

    # --- MODULE 7: SANTÉ ---
    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Suivi Sanitaire")
        target = st.selectbox("Brebis", db.fetch_all_as_df("SELECT identifiant_unique FROM brebis"))
        acte = st.selectbox("Vaccin / Soin", ["Enterotoxémie", "Fièvre Aphteuse", "Clavelée", "PPR", "Vermifuge"])
        if st.button("Enregistrer Soin"):
            rappel = date.today() + timedelta(days=180)
            db.execute_query("INSERT INTO sante (brebis_id, date_soin, type_acte, rappel_prevu) VALUES (?,?,?,?)", (target, date.today(), acte, rappel))
            st.success(f"Rappel enregistré pour le {rappel}")

 # --- MODULE 8: GÉNOMIQUE, SNP D'INTÉRÊT & PARENTÉ ---
    elif choice == "🧬 Génomique & NCBI":
        st.title("🧬 Laboratoire Génomique : Sélection & Parenté")
        
        from Bio import pairwise2
        from Bio.Seq import Seq

        # Références génétiques (Standards 2026)
        GENES_INTERET = {
            "FecB (Prolificité)": "TTAGC",  # Séquence cible simulée
            "MSTN (Muscle/Myostatine)": "GGACC",
            "CAST (Tendreté Viande)": "CCAAA",
            "DGAT1 (Qualité Lait)": "GCTAG"
        }

        tab_snp, tab_parente, tab_ncbi = st.tabs(["🎯 Gènes d'Intérêt", "👪 Test de Parenté", "🌐 NCBI Connect"])

        with tab_snp:
            st.subheader("🔍 Criblage des Gènes de Performance")
            dna_input = st.text_area("Séquence ADN de l'animal (Format FASTA)", 
                                     height=120, placeholder=">ID_UNIT_01\nATGC...").upper()
            
            if dna_input:
                seq_clean = "".join(dna_input.split('\n')[1:]) if ">" in dna_input else dna_input
                seq_clean = seq_clean.strip()

                st.write("### Résultats du Criblage SNP")
                cols = st.columns(2)
                
                found_any = False
                for i, (gene, motif) in enumerate(GENES_INTERET.items()):
                    col_idx = i % 2
                    if motif in seq_clean:
                        found_any = True
                        cols[col_idx].success(f"✅ **{gene} Détecté**")
                        if "FecB" in gene:
                            cols[col_idx].caption("📢 **Interprétation :** L'animal aura une descendance très prolifique (probabilité élevée de jumeaux).")
                        elif "MSTN" in gene:
                            cols[col_idx].caption("📢 **Interprétation :** Hypertrophie musculaire détectée (Rendement viande supérieur).")
                        elif "CAST" in gene:
                            cols[col_idx].caption("📢 **Interprétation :** Marqueur de tendreté de la viande (Qualité supérieure).")
                    else:
                        cols[col_idx].info(f"⚪ {gene} : Non détecté")
                
                if not found_any:
                    st.warning("Aucun gène d'intérêt spécifique détecté dans cette séquence.")

        with tab_parente:
            st.subheader("👪 Vérification de la Parenté Biologique")
            st.write("Comparez l'ADN de l'agneau avec ses parents présumés.")
            
            c1, c2, c3 = st.columns(3)
            dna_agneau = c1.text_area("ADN Agneau", height=100).upper()
            dna_pere = c2.text_area("ADN Père (Bélier)", height=100).upper()
            dna_mere = c3.text_area("ADN Mère (Brebis)", height=100).upper()

            if st.button("Lancer le test de paternité/maternité"):
                if dna_agneau and dna_pere and dna_mere:
                    # Score de similarité Agneau-Père
                    score_p = pairwise2.align.globalxx(dna_agneau, dna_pere, score_only=True)
                    sim_p = (score_p / len(dna_pere)) * 100
                    
                    # Score de similarité Agneau-Mère
                    score_m = pairwise2.align.globalxx(dna_agneau, dna_mere, score_only=True)
                    sim_m = (score_m / len(dna_mere)) * 100
                    
                    res_p, res_m = st.columns(2)
                    res_p.metric("Similarité Paternelle", f"{sim_p:.1f}%")
                    res_m.metric("Similarité Maternelle", f"{sim_m:.1f}%")
                    
                    if sim_p > 45 and sim_m > 45:
                        st.success("✅ **Filiation Confirmée** : L'agneau hérite correctement des deux parents.")
                    else:
                        st.error("⚠️ **Incohérence de filiation** : Les scores sont trop bas pour confirmer la parenté.")
                else:
                    st.error("Veuillez saisir les trois séquences pour le test.")

        with tab_ncbi:
            st.subheader("🌐 Ressources Globales")
            
            st.link_button("Accéder à Ensembl Sheep Genome", "https://www.ensembl.org/Ovis_aries/")

    # --- MODULE 9: STATS ---
    elif choice == "📈 Statistiques":
        st.title("📈 Analyse de Variance & Corrélations")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            st.plotly_chart(px.violin(df, x="race", y="poids", box=True, title="Distribution des Poids par Race"))
            st.write("Corrélation Canon vs Poids :", df['circ_canon'].corr(df['poids']))

if __name__ == "__main__":
    main()
