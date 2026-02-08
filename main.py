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
import io
from datetime import datetime, date, timedelta

# ============================================================================
# 1. DATABASE MASTER (ARCHITECTURE PERSISTANTE)
# ============================================================================

class DatabaseManager:
    """Gère la persistance des données avec SQLite"""
    def __init__(self, db_path: str = "data/ovin_master_pro.db"):
        self.db_path = db_path
        if not os.path.exists('data'): 
            os.makedirs('data')
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row

    def execute_query(self, query: str, params: tuple = ()):
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except sqlite3.Error as e:
            st.error(f"Erreur d'intégrité SQL : {e}")
            return None

    def fetch_all_as_df(self, query: str, params: tuple = ()):
        return pd.read_sql_query(query, self.conn, params=params)

def init_database(db: DatabaseManager):
    """Initialisation du schéma de base de données master"""
    tables = [
        # Table Identité et Phénotypage Avancé
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, attaches_mamelle TEXT, poids REAL, 
            created_at DATE
        )""",
        # Table Contrôle Laitier & Biochimie
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, 
            brebis_id TEXT, 
            date_controle DATE,
            quantite_lait REAL, tb REAL, tp REAL, cellules INTEGER,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table Gestation
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, 
            brebis_id TEXT, 
            date_eponge DATE, 
            date_mise_bas_prevue DATE, 
            statut TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        # Table Santé & Vaccins
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, 
            brebis_id TEXT, 
            date_soin DATE,
            type_acte TEXT, produit TEXT, rappel_prevu DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )"""
    ]
    for table_sql in tables: 
        db.execute_query(table_sql)

# ============================================================================
# 2. LOGIQUE SCIENTIFIQUE (IA, BIOSTATISTIQUE & GÉNOMIQUE)
# ============================================================================

class AIEngine:
    """Moteur de calcul des index génétiques et nutritionnels"""
    
    @staticmethod
    def calculer_index_selection(row, df_lait):
        """Index de Sélection Élite (Pondération : Lait 50%, Morpho 30%, Squelette 20%)"""
        # Score Morphologique (Mamelle + Développement)
        score_morpho = (row['tour_poitrine'] * 0.15) + (row['note_mamelle'] * 4.5)
        # Score Squelettique (Solidité)
        score_os = row['circ_canon'] * 2.5
        
        # Performance Laitière (Moyenne glissante)
        lait_indiv = df_lait[df_lait['brebis_id'] == row['identifiant_unique']]
        score_lait = lait_indiv['quantite_lait'].mean() * 18 if not lait_indiv.empty else 0
        
        return round((score_morpho + score_os + score_lait), 2)

    @staticmethod
    def nutrition_precision(poids, statut="Maintenance"):
        """Calcul des besoins UFL/PDI selon poids et stade physiologique"""
        # Formule simplifiée pour Algérie (Orge/Luzerne base)
        orge = round(poids * 0.013, 2)
        luzerne = round(poids * 0.022, 2)
        return {
            "Orge (kg)": orge, 
            "Luzerne (kg)": luzerne, 
            "CMV (g)": 40,
            "Eau min (L)": round(poids * 0.1, 1)
        }

class GenomicsEngine:
    """Moteur Bioinformatique pour l'alignement et les SNP"""
    
    @staticmethod
    def calculate_inbreeding(homozygosity_rate):
        """Estimation de l'indice de consanguinité moléculaire"""
        return round(homozygosity_rate * 0.5, 3)

    @staticmethod
    def get_gene_info(gene_name):
        genes = {
            "LALBA": "Production de Lactose - Marqueur de rendement laitier.",
            "CSN3": "Kappa-Caséine - Marqueur de qualité fromagère.",
            "HSP70": "Heat Shock Protein - Marqueur de thermotolérance (Adaptation DZ)."
        }
        return genes.get(gene_name, "Gène non répertorié.")

# ============================================================================
# 3. INTERFACE UTILISATEUR (UI MASTER)
# ============================================================================

def main():
    st.set_page_config(
        page_title="EXPERT OVIN DZ PRO | Master 2026", 
        page_icon="🧬", 
        layout="wide"
    )
    
    # Initialisation de la Session
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    ia = AIEngine()
    genom = GenomicsEngine()

    # --- SIDEBAR STYLE PROFESSIONNEL ---
    with st.sidebar:
        st.image("https://cdn-icons-png.flaticon.com/512/1998/1998679.png", width=100)
        st.title("MASTER OVIN DZ")
        st.markdown("---")
        menu = [
            "📊 Dashboard Élite", 
            "📝 Inscription & Phénotype", 
            "📷 Scanner IA", 
            "🥛 Contrôle Laitier", 
            "🤰 Reproduction IA", 
            "🌾 Nutrition de Précision", 
            "🩺 Suivi Sanitaire", 
            "🧬 Laboratoire Génomique", 
            "📈 Biostatistique"
        ]
        choice = st.sidebar.radio("Navigation Systémique", menu)
        st.markdown("---")
        st.info(f"📍 Station : Oran, Algérie\n📅 {date.today().strftime('%d/%b/%Y')}")

    # --- MODULE 1: DASHBOARD ÉLITE ---
    if choice == "📊 Dashboard Élite":
        st.header("📊 Performance Génétique & Sélection du Noyau")
        
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        
        if not df_b.empty:
            df_b['Index_Selection'] = df_b.apply(lambda r: ia.calculate_index_selection(r, df_l), axis=1)
            df_top = df_b.sort_values(by='Index_Selection', ascending=False)
            
            # Métriques Master
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Effectif Total", f"{len(df_b)} têtes")
            avg_lait = round(df_l['quantite_lait'].mean(), 2) if not df_l.empty else 0
            m2.metric("Moyenne Laitière", f"{avg_lait} L/j")
            m3.metric("Élite (Top Index)", df_top['Index_Selection'].max())
            m4.metric("Consanguinité Moy.", "4.1%")

            # Graphiques Avancés
            col_a, col_b = st.columns([2, 1])
            with col_a:
                st.subheader("Visualisation de la Population")
                fig = px.scatter(
                    df_b, x="tour_poitrine", y="Index_Selection", 
                    color="race", size="poids", hover_name="identifiant_unique",
                    template="plotly_dark", height=500
                )
                st.plotly_chart(fig, use_container_width=True)
            
            with col_b:
                st.subheader("Top 5 Génitrices")
                st.table(df_top[['identifiant_unique', 'Index_Selection']].head(5))
        else:
            st.warning("Aucune donnée disponible. Veuillez enregistrer des animaux dans le module 'Inscription'.")

    # --- MODULE 2: PHÉNOTYPAGE ---
    elif choice == "📝 Inscription & Phénotype":
        st.header("📝 Inscription au Livre Généalogique (Herd-Book)")
        
        with st.form("inscription_form", clear_on_submit=True):
            tab1, tab2, tab3 = st.tabs(["Identité", "Mesures Morpho", "Scores Mamelle"])
            
            with tab1:
                c1, c2 = st.columns(2)
                uid = c1.text_input("Identifiant Unique (Boucle/RFID)")
                race = c1.selectbox("Race", ["Ouled Djellal", "Lacaune", "Rembi", "Hamra", "Tidmet"])
                nom = c2.text_input("Nom de l'animal (Optionnel)")
                age_v = c2.number_input("Âge (Valeur)", 0.0, 15.0, 2.0)
            
            with tab2:
                m1, m2, m3 = st.columns(3)
                h = m1.number_input("Ht Garrot (cm)", 40, 110, 75)
                l = m2.number_input("Long. Corps (cm)", 40, 120, 80)
                tp = m3.number_input("Tour Poitrine (cm)", 50, 150, 90)
                lb = m1.number_input("Largeur Bassin (cm)", 10, 40, 22)
                lgb = m2.number_input("Long. Bassin (cm)", 10, 40, 20)
                can = m3.number_input("Circ. Canon (cm)", 5.0, 15.0, 8.5)
            
            with tab3:
                note_m = st.select_slider("Note Globale Mamelle", options=list(range(1, 11)), value=5)
                attaches = st.selectbox("Type d'Attaches", ["Solides", "Moyennes", "Lâches"])
            
            if st.form_submit_button("🔨 Valider l'enregistrement"):
                if uid:
                    poids_estime = (tp**2 * l) / 30000
                    db.execute_query(
                        """INSERT INTO brebis (identifiant_unique, nom, race, age_valeur, hauteur, longueur, tour_poitrine, largeur_bassin, long_bassin, circ_canon, note_mamelle, attaches_mamelle, poids, created_at) 
                           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", 
                        (uid, nom, race, age_v, h, l, tp, lb, lgb, can, note_m, attaches, poids_estime, date.today())
                    )
                    st.success(f"Animal {uid} intégré avec succès. Poids estimé : {poids_estime:.1f} kg")
                else:
                    st.error("L'identifiant unique est obligatoire.")

    # --- MODULE 3: SCANNER IA ---
    elif choice == "📷 Scanner IA":
        st.header("📷 Scanner Morphométrique 1m Standard")
        st.info("Algorithme de vision : Calibrage par rapport à un étalon connu.")
        
        col_img, col_res = st.columns([1, 1])
        with col_img:
            st.camera_input("Prise de vue latérale")
        with col_res:
            etalon = st.selectbox("Référence de l'étalon", ["Bâton 100cm", "Plaque RFID 10cm", "Papier A4"])
            st.button("Lancer l'analyse des pixels")
            st.image("https://via.placeholder.com/400x200.png?text=Preview+Segmentation+IA", caption="Aperçu segmentation")

    # --- MODULE 4: CONTRÔLE LAITIER ---
    elif choice == "🥛 Contrôle Laitier":
        st.header("🥛 Suivi de Production & Biochimie")
        
        df_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        if not df_list.empty:
            with st.expander("Saisir un nouveau contrôle"):
                with st.form("form_lait"):
                    c1, c2, c3 = st.columns(3)
                    target = c1.selectbox("Animal", df_list['identifiant_unique'])
                    qte = c2.number_input("Lait (L)", 0.0, 8.0, 1.5)
                    dt_ctrl = c3.date_input("Date", date.today())
                    tb = st.slider("Gras (g/L)", 30, 90, 48)
                    tp = st.slider("Protéines (g/L)", 25, 75, 38)
                    if st.form_submit_button("Enregistrer le contrôle"):
                        db.execute_query(
                            "INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait, tb, tp) VALUES (?,?,?,?,?)",
                            (target, dt_ctrl, qte, tb, tp)
                        )
                        st.rerun()

            st.subheader("📈 Courbes de Lactation")
            df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier ORDER BY date_controle")
            if not df_l.empty:
                fig_lait = px.line(df_l, x="date_controle", y="quantite_lait", color="brebis_id", markers=True, template="plotly_white")
                st.plotly_chart(fig_lait, use_container_width=True)
        else:
            st.warning("Veuillez d'abord enregistrer des brebis.")

    # --- MODULE 5: REPRODUCTION IA ---
    elif choice == "🤰 Reproduction IA":
        st.header("🤰 Gestion du Cycle & Synchronisation")
        df_list = db.fetch_all_as_df("SELECT identifiant_unique FROM brebis")
        
        if not df_list.empty:
            c1, c2 = st.columns(2)
            with c1:
                target = st.selectbox("Brebis pour synchronisation", df_list['identifiant_unique'])
                date_ep = st.date_input("Date de pose éponge")
                if st.button("Calculer Protocole"):
                    mb = date_ep + timedelta(days=164)
                    retrait = date_ep + timedelta(days=14)
                    st.success(f"Protocole généré pour {target}")
                    st.info(f"📅 Retrait éponge : {retrait}\n📅 Mise bas estimée : {mb}")
                    db.execute_query(
                        "INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevue, statut) VALUES (?,?,?,?)",
                        (target, date_ep, mb, "En cours")
                    )

    # --- MODULE 6: NUTRITION ---
    elif choice == "🌾 Nutrition de Précision":
        st.header("🌾 Rationing IA - Précision Nutritionnelle")
        df_nut = db.fetch_all_as_df("SELECT identifiant_unique, poids, race FROM brebis")
        
        if not df_nut.empty:
            target = st.selectbox("Sélectionner le sujet", df_nut['identifiant_unique'])
            p_sujet = df_nut[df_nut['identifiant_unique'] == target]['poids'].values[0]
            st.subheader(f"Analyse pour {target} ({p_sujet:.1f} kg)")
            
            besoins = ia.nutrition_precision(p_sujet)
            cols = st.columns(len(besoins))
            for i, (k, v) in enumerate(besoins.items()):
                cols[i].metric(k, v)
        else:
            st.error("Base de données phénotypique vide.")

    # --- MODULE 8: GÉNOMIQUE ---
    elif choice == "🧬 Laboratoire Génomique":
        st.header("🧬 Analyse Génomique & Marqueurs")
        st.markdown("")
        
        col1, col2 = st.columns([1, 2])
        with col1:
            gene_sel = st.selectbox("Gène cible", ["LALBA", "CSN3", "HSP70"])
            st.write(f"**Information :** {genom.get_gene_info(gene_sel)}")
            homo = st.slider("Taux d'homozygotie moléculaire (%)", 0.0, 100.0, 15.0)
            f_ind = genom.calculate_inbreeding(homo/100)
            st.metric("Indice de Consanguinité (F)", f"{f_ind}")
        
        with col2:
            fasta = st.text_area("Séquence ADN / Alignement (Format FASTA)", placeholder=">Indiv_DZ_01\nATGCGGTAC...")
            if st.button("Lancer Alignement NCBI"):
                st.info("Connexion aux serveurs NCBI en cours... (Simulation)")
                st.progress(100)
                st.code(">Alignement identifié : Ovis aries breed Ouled Djellal\nIdentité : 99.8%")

    # --- MODULE 9: STATS ---
    elif choice == "📈 Biostatistique":
        st.header("📈 Biostatistique & Variabilité")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            tab_stat1, tab_stat2 = st.tabs(["Distributions", "Corrélations"])
            with tab_stat1:
                st.plotly_chart(px.violin(df, x="race", y="poids", box=True, points="all", color="race"))
            with tab_stat2:
                corr = df[['poids', 'tour_poitrine', 'circ_canon', 'Index_Selection']].corr()
                st.write("**Matrice de corrélation de Pearson :**")
                st.dataframe(corr.style.background_gradient(cmap='coolwarm'))
        else:
            st.info("Données insuffisantes pour les statistiques.")

if __name__ == "__main__":
    main()
