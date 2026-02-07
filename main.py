"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.02
Système : Phénotypage, Lait, Santé & Génotypage Sélectif
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date, timedelta
from Bio import Align  
from Bio.Seq import Seq

# ============================================================================
# 1. GESTION DE LA BASE DE DONNÉES
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
            nom TEXT, race TEXT, poids REAL, note_mamelle INTEGER, 
            tour_poitrine REAL, longueur REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS genotypes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, gene_nom TEXT, resultat TEXT, 
            score_homologie REAL, classement TEXT, date_test DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, quantite_lait REAL
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. MOTEUR BIOINFORMATIQUE (ALIGNEUR HAUTE PERFORMANCE)
# ============================================================================

class BioInfoEngine:
    GENES_QUALITE = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Muscle/Viande)": "AAGCTTGATTAGCAGGTTCCCGG",
        "CAST (Tendreté)": "TGGGGCCCAAGTCGATTGCAGAA",
        "LALBA (Lait)": "GCTAGCTAGCTAGCTGATCGATG"
    }
    
    GENES_SANTE = {
        "Scrapie_ARR (Résistance)": "TGGTACCCATAATCAGTGGAACA",
        "Scrapie_VRQ (Sensibilité)": "TGGTAGCCATAATCAGTGGAACA",
        "Arachnomélie": "CCGTAGCTAGCTGATCGATCGTA"
    }

    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'

    def filtrer_sequence(self, seq):
        if ">" in seq:
            seq = "".join(seq.split('\n')[1:])
        return seq.upper().strip().replace(" ", "").replace("\n", "").replace("\r", "")

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🧬")
    
    # Initialisation DB et Moteur
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    db = st.session_state.db
    genomique = st.session_state.genomique

    # Navigation latérale
    st.sidebar.title("🐑 EXPERT OVIN DZ")
    st.sidebar.markdown(f"**Date:** {date.today()}")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA 1m", 
        "🧬 Génotypage Sélectif", 
        "🥛 Contrôle Laitier", 
        "🩺 Santé & Vaccins"
    ]
    choice = st.sidebar.radio("Menu Principal", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Performances et Analyses du Cheptel")
        t1, t2 = st.tabs(["📈 Production & Effectif", "🧬 Archive Génotypique"])
        
        with t1:
            df_b = db.fetch_all_as_df("SELECT * FROM brebis")
            if not df_b.empty:
                st.dataframe(df_b, width='stretch')
            else: st.info("Aucun animal enregistré.")
            
        with t2:
            df_g = db.fetch_all_as_df("SELECT * FROM genotypes ORDER BY date_test DESC")
            if not df_g.empty:
                st.dataframe(df_g, width='stretch')
            else: st.info("Aucune donnée génomique disponible.")

    # --- 2. INSCRIPTION ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Enregistrement Phénotypique")
        with st.form("inscription"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle")
            nom = c1.text_input("Nom de l'animal")
            race = c2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Autre"])
            poids = c2.number_input("Poids (kg)", 10.0, 150.0, 50.0)
            if st.form_submit_button("Sauvegarder"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, nom, race, poids, created_at) VALUES (?,?,?,?,?)",
                                 (uid, nom, race, poids, date.today()))
                st.success("Animal ajouté !")

    # --- 3. SCANNER IA ---
    elif choice == "📷 Scanner IA 1m":
        st.title("📸 Scanner IA Morphométrique")
        st.info("Placez l'étalon de 1m parallèlement à l'animal.")
        st.camera_input("Capturer")

    # --- 4. GÉNOTYPAGE SÉLECTIF (LE NOUVEAU BLOC) ---
    elif choice == "🧬 Génotypage Sélectif":
        st.title("🧬 Laboratoire de Génotypage à la Carte")
        st.markdown("---")
        
        col_in, col_out = st.columns([1, 1.2])

        with col_in:
            st.subheader("⚙️ Paramètres")
            id_animal = st.text_input("ID de l'animal à analyser")
            dna_raw = st.text_area("Séquence ADN (NCBI / FASTA)", height=200, help="Collez ici la séquence texte.")
            
            tous_marqueurs = {**genomique.GENES_QUALITE, **genomique.GENES_SANTE}
            choix = st.multiselect(
                "Marqueurs à cibler :",
                options=list(tous_marqueurs.keys()),
                default=["FecB (Prolificité)", "Scrapie_ARR (Résistance)"]
            )
            
            bouton = st.button("🚀 Lancer l'Analyse Moléculaire")

        with col_out:
            st.subheader("📋 Diagnostic Final")
            if bouton and dna_raw and id_animal:
                with st.spinner("Alignement des séquences..."):
                    seq_clean = genomique.filtrer_sequence(dna_raw)
                    resultats_final = []

                    for m_nom in choix:
                        ref_seq = tous_marqueurs[m_nom]
                        score = round((genomique.aligner.score(seq_clean, ref_seq) / len(ref_seq)) * 100, 2)
                        
                        # Logique de Classification
                        if "ARR" in m_nom:
                            statut = "🏆 R1 (TRÈS RÉSISTANT)" if score > 88 else "CLASSIQUE"
                        elif "VRQ" in m_nom:
                            statut = "❌ SENSIBLE (À ÉLIMINER)" if score > 88 else "SAIN"
                        elif score > 85:
                            statut = "💎 ÉLITE / PRÉSENT"
                        elif score > 60:
                            statut = "🔸 PORTEUR SIMPLE"
                        else:
                            statut = "➖ NÉGATIF"

                        # Enregistrement DB
                        db.execute_query(
                            "INSERT INTO genotypes (brebis_id, gene_nom, resultat, score_homologie, classement, date_test) VALUES (?,?,?,?,?,?)",
                            (id_animal, m_nom, "CIBLÉ", score, statut, date.today())
                        )
                        
                        resultats_final.append({
                            "Marqueur": m_nom,
                            "Homologie": f"{score}%",
                            "Diagnostic": statut
                        })

                    st.table(pd.DataFrame(resultats_final))
                    st.balloons()
            else:
                st.info("Veuillez entrer l'ID de l'animal et sa séquence ADN pour voir le diagnostic.")

    # --- 5. LAIT & 6. SANTÉ (STUBS) ---
    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Suivi Laitier")
        st.write("Enregistrez vos traites ici.")
    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Carnet de Santé")
        st.write("Gestion des vaccins et traitements.")

if __name__ == "__main__":
    main()
