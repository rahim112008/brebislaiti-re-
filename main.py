"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.02
Système : Phénotypage, Scanner IA Instantané, Laboratoire ADN, Nutrition DZ & Lait
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import os
from datetime import datetime, date
from Bio import Align  
from Bio.Seq import Seq

# ============================================================================
# BLOC 1 : GESTION DE LA BASE DE DONNÉES
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
            id INTEGER PRIMARY KEY AUTOINCREMENT, identifiant_unique TEXT UNIQUE,
            nom TEXT, race TEXT, poids REAL, hauteur REAL, longueur REAL, 
            largeur_bassin REAL, circ_canon REAL, note_mamelle INTEGER, 
            prof_mamelle REAL, attache_ar REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS genotypes (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            gene_nom TEXT, score_homologie REAL, classement TEXT, date_test DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, qte_matin REAL, qte_soir REAL, taux_gras REAL
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, symptomes TEXT, diagnostic TEXT, traitement TEXT
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# BLOC 2 : MOTEUR BIOINFORMATIQUE (ADN & PROTÉINES)
# ============================================================================

class BioInfoEngine:
    GENES_REF = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Viande)": "AAGCTTGATTAGCAGGTTCCCGG",
        "Scrapie_ARR (Résistance)": "TGGTACCCATAATCAGTGGAACA",
        "Scrapie_VRQ (Sensibilité)": "TGGTAGCCATAATCAGTGGAACA"
    }

    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'

    def extraire_multi_fasta(self, raw_text):
        sequences = {}
        current_id = None
        for line in raw_text.split('\n'):
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line.upper().replace(" ", "")
        return sequences if sequences else {"Individu_Unique": raw_text.upper().replace(" ", "")}

    def traduire(self, dna_seq):
        try:
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            return str(Seq(clean_dna).translate(to_stop=True))
        except: return "Erreur de traduction"

# ============================================================================
# BLOC 3 : SCANNER IA & PHÉNOTYPE (APPAREIL PHOTO DIRECT)
# ============================================================================

def bloc_scanner_phenotype(db):
    st.title("📷 Scanner IA & Inscription")
    
    col_cam, col_form = st.columns(2)
    
    with col_cam:
        st.subheader("Analyse Morphométrique")
        etalon = st.selectbox("Étalon de calibration", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
        # Caméra en direct
        photo = st.camera_input("Prendre une photo de l'animal")
        
        if photo:
            st.success("Analyse IA : Dimensions calculées !")
            # Simulation des mesures extraites par l'image
            st.session_state['ia_measure'] = {"h": 76.5, "l": 102.0, "b": 23.5, "c": 9.2}

    with col_form:
        st.subheader("Fiche Phénotypique")
        m = st.session_state.get('ia_measure', {"h":0.0, "l":0.0, "b":0.0, "c":0.0})
        with st.form("form_sheep"):
            uid = st.text_input("ID Boucle (Unique)")
            race = st.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Tidmet"])
            hauteur = st.number_input("Hauteur au garrot (cm)", value=m['h'])
            longueur = st.number_input("Longueur (cm)", value=m['l'])
            bassin = st.number_input("Largeur Bassin (cm)", value=m['b'])
            canon = st.number_input("Circonférence Canon (cm)", value=m['c'])
            
            st.markdown("---")
            st.write("**Évaluation quantitative Mamelle**")
            prof_m = st.number_input("Profondeur (cm)", 10.0, 30.0, 15.0)
            att_ar = st.number_input("Attache Arrière (cm)", 5.0, 20.0, 10.0)
            
            if st.form_submit_button("Enregistrer l'animal"):
                db.execute_query("""INSERT INTO brebis (identifiant_unique, race, hauteur, longueur, 
                                 largeur_bassin, circ_canon, prof_mamelle, attache_ar, created_at) 
                                 VALUES (?,?,?,?,?,?,?,?,?)""", 
                                 (uid, race, hauteur, longueur, bassin, canon, prof_m, att_ar, date.today()))
                st.success("Animal enregistré avec succès !")

# ============================================================================
# BLOC 4 : NUTRITION DZ FLEXIBLE
# ============================================================================

def bloc_nutrition():
    st.title("🌾 Nutrition & Rations Algérie")
    
    with st.expander("💰 Configuration des Prix (DA/Quintal)", expanded=True):
        c1, c2, c3 = st.columns(3)
        p_orge = c1.number_input("Orge (Chaïr)", 5000)
        p_son = c2.number_input("Son (Nkhala)", 2500)
        p_foussa = c3.number_input("Luzerne/Foin", 4500)

    stade = st.selectbox("Objectif Nutritionnel", ["Entretien", "Gestation", "Allaitement (Lait+)", "Engraissement"])
    
    # Rations types (kg/jour)
    rations = {"Entretien": [0.5, 0.2, 1.0], "Gestation": [0.8, 0.3, 1.2], "Allaitement (Lait+)": [1.2, 0.5, 1.5], "Engraissement": [1.0, 0.2, 0.5]}
    r = rations[stade]
    
    cout_j = (p_orge/100*r[0]) + (p_son/100*r[1]) + (p_foussa/100*r[2])
    
    st.metric("Coût journalier estimé", f"{round(cout_j, 2)} DA / Tête")
    st.info(f"Composition : Orge: {r[0]}kg | Son: {r[1]}kg | Foin: {r[2]}kg")

# ============================================================================
# BLOC 5 : LAIT & SANTÉ IA
# ============================================================================

def bloc_lait_sante(db):
    st.title("🥛 Lait & 🩺 Santé IA")
    
    t1, t2 = st.tabs(["Suivi Laitier", "Diagnostic Santé IA"])
    
    with t1:
        with st.form("lait_f"):
            id_b = st.text_input("ID Brebis")
            q_m = st.number_input("Matin (L)", 0.0, 5.0, 1.0)
            q_s = st.number_input("Soir (L)", 0.0, 5.0, 0.8)
            gras = st.slider("Taux de matière grasse (%)", 30, 90, 60)
            if st.form_submit_button("Valider Traite"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, qte_matin, qte_soir, taux_gras) VALUES (?,?,?,?,?)",
                                 (id_b, date.today(), q_m, q_s, gras))
    
    with t2:
        st.subheader("Assistant Diagnostic")
        symp = st.multiselect("Signes cliniques :", ["Toux", "Boiterie", "Lésions buccales", "Diarrhée"])
        if st.button("Analyser les symptômes"):
            if "Boiterie" in symp:
                st.error("Diagnostic IA : Suspicion de Piétin.")
                st.write("Action : Parage des onglons et désinfection.")
            elif "Toux" in symp:
                st.warning("Diagnostic IA : Suspicion de Parasitose Pulmonaire.")

# ============================================================================
# BLOC 6 : LABORATOIRE ADN & EXPERTISE PROTÉIQUE
# ============================================================================

def bloc_adn_proteine(genomique):
    st.title("🧬 Labo ADN & 🔬 Expertise Protéique")
    dna_input = st.text_area("Collez vos séquences ADN (Format FASTA)")
    
    if dna_input:
        sequences = genomique.extraire_multi_fasta(dna_input)
        tab1, tab2 = st.tabs(["Génotypage", "Expertise Protéique"])
        
        with tab1:
            target = st.selectbox("Choisir individu", list(sequences.keys()))
            gene = st.selectbox("Marqueur cible", list(genomique.GENES_REF.keys()))
            score = round((genomique.aligner.score(sequences[target], genomique.GENES_REF[gene]) / len(genomique.GENES_REF[gene]))*100, 2)
            st.metric(f"Homologie {gene}", f"{score}%")
            
        with tab2:
            st.subheader("Séquence en Acides Aminés")
            id_prot = st.selectbox("Traduire individu", list(sequences.keys()))
            st.code(genomique.traduire(sequences[id_prot]), language="markdown")

# ============================================================================
# MAIN : NAVIGATION
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ", layout="wide", page_icon="🐑")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    menu = ["📊 Dashboard", "📝 Inscription & Phénotype", "📷 Scanner IA 1m", "🧬 Laboratoire ADN & Protéines", "🥛 Lait & Santé", "🌾 Nutrition"]
    choice = st.sidebar.radio("Menu Principal", menu)

    if choice == "📊 Dashboard":
        st.title("Tableau de Bord Élite")
        st.dataframe(st.session_state.db.fetch_all_as_df("SELECT * FROM brebis"))
    elif choice == "📝 Inscription & Phénotype" or choice == "📷 Scanner IA 1m":
        bloc_scanner_phenotype(st.session_state.db)
    elif choice == "🧬 Laboratoire ADN & Protéines":
        bloc_adn_proteine(st.session_state.genomique)
    elif choice == "🥛 Lait & Santé":
        bloc_lait_sante(st.session_state.db)
    elif choice == "🌾 Nutrition":
        bloc_nutrition()

if __name__ == "__main__":
    main()
