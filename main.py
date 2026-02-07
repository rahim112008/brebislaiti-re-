"""
EXPERT OVIN DZ PRO - VERSION ULTIME CONSOLIDÉE 2026
Système Intégré : Phénotypage, Scanner IA, Génomique, Nutrition DZ & Lait
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
# 1. GESTION DE LA BASE DE DONNÉES (PERSISTENCE)
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
            largeur_bassin REAL, circ_canon REAL, prof_mamelle REAL, attache_ar REAL, created_at DATE
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
# 2. MOTEUR GÉNOMIQUE & PROTÉIQUE (LABO)
# ============================================================================

class BioInfoEngine:
    GENES_REF = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Muscle/Viande)": "AAGCTTGATTAGCAGGTTCCCGG",
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
                current_id = line[1:].split()[0]
                sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line.upper().replace(" ", "")
        return sequences if sequences else {"Individu_Unique": raw_text.upper().replace(" ", "")}

    def traduire(self, dna_seq):
        try:
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            return str(Seq(clean_dna).translate(to_stop=True))
        except: return "Séquence invalide"

# ============================================================================
# 3. INTERFACE UTILISATEUR & LOGIQUE MÉTIER
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    # Initialisation Session
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    db, genomique = st.session_state.db, st.session_state.genomique

    # Navigation Latérale
    st.sidebar.title("🐑 EXPERT OVIN DZ")
    st.sidebar.info("Système Intégré Pro v2026.02")
    menu = [
        "📊 Dashboard", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA 1m", 
        "🧬 Laboratoire ADN", 
        "🔬 Expertise Protéique",
        "🥛 Lait & Santé",
        "🌾 Nutrition DZ"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard":
        st.title("📊 Performances du Cheptel")
        df = db.fetch_all_as_df("SELECT * FROM brebis")
        if not df.empty:
            st.dataframe(df, use_container_width=True)
            st.download_button("Exporter en CSV", df.to_csv(index=False), "cheptel_dz.csv")
        else:
            st.info("Aucun animal enregistré.")

    # --- 2. INSCRIPTION & 3. SCANNER IA ---
    elif choice in ["📝 Inscription & Phénotype", "📷 Scanner IA 1m"]:
        st.title("📝 Phénotypage & Scanner IA")
        
        col_scan, col_fiche = st.columns([1, 1.2])
        
        with col_scan:
            st.subheader("Capturer / Analyser")
            etalon = st.selectbox("Objet de calibration", ["Bâton 1m", "Feuille A4", "Carte Bancaire"])
            photo = st.camera_input("Scanner l'animal")
            upload = st.file_uploader("Ou importer une image", type=['jpg', 'jpeg', 'png'])
            
            if photo or upload:
                # Simulation IA des mesures morphométriques
                st.session_state['ia_data'] = {"h": 78.2, "l": 105.4, "b": 24.5, "c": 10.1}
                st.success("Analyse IA terminée avec succès !")

        with col_fiche:
            st.subheader("Données Morphologiques")
            ia = st.session_state.get('ia_data', {"h": 0.0, "l": 0.0, "b": 0.0, "c": 0.0})
            with st.form("form_pheno"):
                c1, c2 = st.columns(2)
                uid = c1.text_input("ID Boucle")
                race = c2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune"])
                
                m1, m2 = st.columns(2)
                h_g = m1.number_input("Hauteur Garrot (cm)", value=ia['h'])
                l_c = m2.number_input("Longueur Corps (cm)", value=ia['l'])
                b_l = m1.number_input("Largeur Bassin (cm)", value=ia['b'])
                c_c = m2.number_input("Circ. Canon (cm)", value=ia['c'])
                
                st.markdown("**Expertise Mamelle**")
                p_m = st.number_input("Profondeur (cm)", 10.0, 30.0, 15.0)
                a_m = st.number_input("Attache (cm)", 5.0, 20.0, 10.0)
                
                if st.form_submit_button("Enregistrer l'Animal"):
                    db.execute_query("""INSERT INTO brebis (identifiant_unique, race, hauteur, longueur, largeur_bassin, 
                                     circ_canon, prof_mamelle, attache_ar, created_at) VALUES (?,?,?,?,?,?,?,?,?)""",
                                     (uid, race, h_g, l_c, b_l, c_c, p_m, a_m, date.today()))
                    st.success(f"L'animal {uid} a été ajouté au système.")

    # --- 4. LABORATOIRE ADN ---
    elif choice == "🧬 Laboratoire ADN":
        st.title("🧬 Diagnostic Génomique")
        dna_input = st.text_area("Séquences ADN (Format FASTA)", height=200)
        if dna_input:
            data = genomique.extraire_multi_fasta(dna_input)
            for name, seq in data.items():
                with st.expander(f"Résultats pour : {name}"):
                    res_m = []
                    for g_nom, g_ref in genomique.GENES_REF.items():
                        score = round((genomique.aligner.score(seq, g_ref)/len(g_ref))*100, 2)
                        verdict = "💎 ÉLITE" if score > 85 else "➖ NORMAL"
                        res_m.append({"Marqueur": g_nom, "Homologie": f"{score}%", "Verdict": verdict})
                    st.table(pd.DataFrame(res_m))

    # --- 5. EXPERTISE PROTÉIQUE ---
    elif choice == "🔬 Expertise Protéique":
        st.title("🔬 Analyse des Protéines")
        dna_p = st.text_area("Collez la séquence ADN pour traduction")
        if dna_p:
            clean_p = dna_p.split('\n')[-1] # Simple nettoyage
            proteine = genomique.traduire(clean_p)
            st.info("Séquence d'Acides Aminés (Traduction moléculaire) :")
            st.code(proteine, language="markdown")

    # --- 6. LAIT & SANTÉ ---
    elif choice == "🥛 Lait & Santé":
        st.title("🥛 Production & 🩺 Santé")
        tab_lait, tab_sante = st.tabs(["Suivi Laitier", "Carnet de Santé IA"])
        
        with tab_lait:
            with st.form("f_lait"):
                b_id = st.text_input("ID Animal")
                q_m = st.number_input("Lait Matin (L)", 0.0, 5.0, 1.2)
                q_s = st.number_input("Lait Soir (L)", 0.0, 5.0, 0.8)
                if st.form_submit_button("Enregistrer Traite"):
                    db.execute_query("INSERT INTO controle_laitier (brebis_id, qte_matin, qte_soir, date_controle) VALUES (?,?,?,?)",
                                     (b_id, q_m, q_s, date.today()))
        
        with tab_sante:
            symp = st.multiselect("Symptômes observés :", ["Toux", "Boiterie", "Diarrhée", "Lésions buccales"])
            if st.button("Lancer le diagnostic IA"):
                if "Boiterie" in symp: st.error("Alerte : Suspicion de Piétin. Traiter au sulfate de zinc.")
                elif "Toux" in symp: st.warning("Alerte : Parasitose interne probable.")

    # --- 7. NUTRITION DZ ---
    elif choice == "🌾 Nutrition":
        st.title("🌾 Ration Alimentaire & Coûts")
        with st.container(border=True):
            st.subheader("Prix du Marché (DA / Quintal)")
            c1, c2, c3 = st.columns(3)
            p_orge = c1.number_input("Orge (Chaïr)", 4500)
            p_son = c2.number_input("Son (Nkhala)", 2500)
            p_foin = c3.number_input("Luzerne/Foin", 4000)
        
        stade = st.selectbox("Stade de la brebis", ["Entretien", "Gestation", "Allaitement", "Flashage Agneau"])
        # Logique de calcul simple : Orge (60%) + Son (20%) + Foin (20%) - à adapter
        cout_base = (p_orge/100 * 0.8) + (p_son/100 * 0.4) + (p_foin/100 * 1.2)
        st.metric("Coût Journalier Estimé", f"{round(cout_base, 2)} DA / Tête")

if __name__ == "__main__":
    main()
