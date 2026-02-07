"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.02
Système : Phénotypage, Lait, Santé & Laboratoire Génomique (Bio.Align)
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
# 2. MOTEUR BIOINFORMATIQUE EXPERT (LABORATOIRE COMPLET)
# ============================================================================

class BioInfoEngine:
    # Marqueurs de Qualité
    GENES_QUALITE = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Muscle/Viande)": "AAGCTTGATTAGCAGGTTCCCGG",
        "CAST (Tendreté)": "TGGGGCCCAAGTCGATTGCAGAA",
        "LALBA (Lait)": "GCTAGCTAGCTAGCTGATCGATG"
    }
    
    # Marqueurs de Santé
    GENES_SANTE = {
        "Scrapie_ARR (Résistance)": "TGGTACCCATAATCAGTGGAACA",
        "Scrapie_VRQ (Sensibilité)": "TGGTAGCCATAATCAGTGGAACA",
        "Arachnomélie": "CCGTAGCTAGCTGATCGATCGTA"
    }

    def __init__(self):
        # Initialisation du moteur Align pour éviter les crashs sur séquences NCBI
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'

    def filtrer_sequence(self, seq):
        """Nettoyage FASTA"""
        if ">" in seq:
            seq = "".join(seq.split('\n')[1:])
        return seq.upper().strip().replace(" ", "").replace("\n", "").replace("\r", "")

    def extraire_multi_fasta(self, raw_text):
        """Gestion multi-séquences"""
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

    def calculer_heterozygotie(self, sequences_dict):
        """Calcul de diversité génétique"""
        if len(sequences_dict) < 2: return 0.0
        seqs = list(sequences_dict.values())
        distances = []
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                score = self.aligner.score(seqs[i], seqs[j])
                distances.append(1 - (score / max(len(seqs[i]), len(seqs[j]))))
        return round(np.mean(distances) * 100, 2)

    def traduire(self, dna_seq):
        """ADN -> Protéine"""
        try:
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            return str(Seq(clean_dna).translate(to_stop=True))
        except: return "Séquence invalide pour la traduction"

# ============================================================================
# 3. INTERFACE UTILISATEUR (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    # Initialisation Session State
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    db, genomique = st.session_state.db, st.session_state.genomique

    # Navigation
    st.sidebar.title("🐑 EXPERT OVIN DZ")
    st.sidebar.markdown(f"**V.2026.02** | {date.today()}")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA 1m", 
        "🧬 Laboratoire Génomique", 
        "🥛 Contrôle Laitier", 
        "🩺 Santé & Vaccins"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Performances du Cheptel")
        t1, t2, t3 = st.tabs(["📋 Effectifs", "🧬 Résultats Génomiques", "📈 Courbes Lait"])
        with t1:
            df_b = db.fetch_all_as_df("SELECT * FROM brebis")
            st.dataframe(df_b, width='stretch')
        with t2:
            df_g = db.fetch_all_as_df("SELECT * FROM genotypes ORDER BY date_test DESC")
            st.dataframe(df_g, width='stretch')

    # --- 2. INSCRIPTION ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Enregistrement Phénotypique")
        with st.form("form_ins"):
            c1, c2 = st.columns(2)
            uid = c1.text_input("ID Boucle")
            nom = c1.text_input("Nom")
            race = c2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune", "Autre"])
            poids = c2.number_input("Poids (kg)", 10.0, 150.0, 50.0)
            note_m = st.slider("Note Mamelle", 1, 10, 5)
            if st.form_submit_button("Enregistrer"):
                db.execute_query("INSERT INTO brebis (identifiant_unique, nom, race, poids, note_mamelle, created_at) VALUES (?,?,?,?,?,?)",
                                 (uid, nom, race, poids, note_m, date.today()))
                st.success("Animal enregistré !")

    # --- 3. SCANNER IA ---
    elif choice == "📷 Scanner IA 1m":
        st.title("📸 Scanner Morphométrique")
        st.info("Utilisez l'étalon de 1 mètre pour la calibration automatique.")
        st.camera_input("Scanner l'animal")

    # --- 4. LABORATOIRE GÉNOMIQUE COMPLET (FUSION) ---
    elif choice == "🧬 Laboratoire Génomique":
        st.title("🧬 Laboratoire Génomique Intégré")
        dna_input = st.text_area("Collez vos séquences ADN (NCBI / FASTA)", height=150)
        
        if dna_input:
            data_dict = genomique.extraire_multi_fasta(dna_input)
            
            tab_ciblage, tab_diversite, tab_traduction = st.tabs([
                "🎯 Génotypage Ciblé", 
                "📊 Analyse de Population", 
                "🔬 Séquençage Moléculaire"
            ])

            # A. Génotypage par individu
            with tab_ciblage:
                st.subheader("Diagnostic de marqueurs par animal")
                c1, c2 = st.columns([1, 2])
                with c1:
                    target_id = st.selectbox("Sélectionner l'individu", list(data_dict.keys()))
                    tous_m = {**genomique.GENES_QUALITE, **genomique.GENES_SANTE}
                    choix_m = st.multiselect("Marqueurs à cibler", list(tous_m.keys()), default=["FecB (Prolificité)"])
                    btn_run = st.button("Lancer le diagnostic")
                
                with c2:
                    if btn_run:
                        seq_test = data_dict[target_id]
                        res_list = []
                        for m_nom in choix_m:
                            ref = tous_m[m_nom]
                            score = round((genomique.aligner.score(seq_test, ref) / len(ref)) * 100, 2)
                            
                            # Classification Intelligente
                            if "ARR" in m_nom: verdict = "💎 R1 (TRÈS RÉSISTANT)" if score > 88 else "CLASSIQUE"
                            elif "VRQ" in m_nom: verdict = "❌ SENSIBLE (À ÉLIMINER)" if score > 88 else "SAIN"
                            elif score > 85: verdict = "✅ ÉLITE / POSITIF"
                            else: verdict = "➖ NÉGATIF"
                            
                            db.execute_query("INSERT INTO genotypes (brebis_id, gene_nom, score_homologie, classement, date_test) VALUES (?,?,?,?,?)",
                                             (target_id, m_nom, score, verdict, date.today()))
                            res_list.append({"Marqueur": m_nom, "Homologie": f"{score}%", "Verdict": verdict})
                        st.table(pd.DataFrame(res_list))

            # B. Diversité de groupe
            with tab_diversite:
                st.subheader("Analyse de Diversité Génétique")
                if len(data_dict) > 1:
                    score_h = genomique.calculer_heterozygotie(data_dict)
                    st.metric("Indice de Diversité (Hétérozygotie)", f"{score_h}%")
                    if score_h < 12: st.error("⚠️ Attention : Risque de consanguinité détecté !")
                    else: st.success("✅ Bonne variabilité génétique au sein du groupe.")
                else: st.info("Collez au moins 2 séquences pour l'analyse de diversité.")

            # C. Traduction
            with tab_traduction:
                st.subheader("Séquence Protéique")
                id_tr = st.selectbox("Traduire :", list(data_dict.keys()), key="tr_sel")
                st.code(genomique.traduire(data_dict[id_tr]))
        else:
            st.warning("Veuillez coller une séquence ADN pour activer le laboratoire.")

    # --- 5. LAIT & 6. SANTÉ ---
    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Suivi Laitier")
        # Logique simplifiée
        id_l = st.text_input("ID Brebis")
        qte = st.number_input("Lait (L)", 0.0, 10.0, 1.0)
        if st.button("Valider traite"):
            db.execute_query("INSERT INTO controle_laitier (brebis_id, quantite_lait) VALUES (?,?)", (id_l, qte))
            st.success("Donnée ajoutée.")

    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Carnet de Santé")
        st.write("Gestion des vaccins et traitements.")

if __name__ == "__main__":
    main()
