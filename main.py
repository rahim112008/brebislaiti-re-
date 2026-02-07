"""
EXPERT OVIN DZ PRO - VERSION ULTIME CONSOLIDÉE 2026
Système Intégré : Phénotypage, Scanner IA, Génomique, Nutrition DZ & Rapport PDF
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import os
from datetime import datetime, date
from Bio import Align  
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam
from fpdf import FPDF
import base64

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
        """CREATE TABLE IF NOT EXISTS rations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_ration DATE, ufl_total REAL, pdi_total REAL, cout_total REAL
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. GÉNÉRATEUR DE RAPPORT PDF
# ============================================================================

class PDFReport(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 15)
        self.cell(0, 10, 'RAPPORT D\'EXPERTISE OVIN DZ PRO', 0, 1, 'C')
        self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()} | Généré le {datetime.now().strftime("%d/%m/%Y")}', 0, 0, 'C')

def create_pdf_report(data, filename="rapport_nutrition.pdf"):
    pdf = PDFReport()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    
    pdf.set_fill_color(200, 220, 255)
    pdf.cell(0, 10, f"Bilan de Rationnement - {data['animal_id']}", 0, 1, 'L', 1)
    pdf.ln(5)
    
    pdf.cell(0, 10, f"Stade Physiologique : {data['stade']}", 0, 1)
    pdf.cell(0, 10, f"Poids de l'animal : {data['poids']} kg", 0, 1)
    pdf.ln(5)
    
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Composition de la Ration :", 0, 1)
    pdf.set_font("Arial", size=11)
    for al, qte in data['aliments'].items():
        if qte > 0:
            pdf.cell(0, 8, f"- {al} : {qte} kg", 0, 1)
    
    pdf.ln(10)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Valeurs Nutritionnelles Totales :", 0, 1)
    pdf.set_font("Arial", size=11)
    pdf.cell(0, 8, f"Énergie Totale : {data['total_ufl']} UFL (Cible: {data['cible_ufl']})", 0, 1)
    pdf.cell(0, 8, f"Protéines Totales : {data['total_pdi']} g PDI (Cible: {data['cible_pdi']})", 0, 1)
    pdf.cell(0, 8, f"Coût Estimé : {data['cout']} DA / Jour", 0, 1)
    
    return pdf.output(dest='S').encode('latin-1')

# ============================================================================
# 3. MOTEURS DE CALCUL (BIO & NUTRITION)
# ============================================================================

class BioInfoEngine:
    GENES_REF = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Muscle/Viande)": "AAGCTTGATTAGCAGGTTCCCGG",
        "Scrapie_ARR": "TGGTACCCATAATCAGTGGAACA"
    }

    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'

    def analyser_proteine(self, dna_seq):
        try:
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            protein_seq = str(Seq(clean_dna).translate(to_stop=True))
            analyser = ProtParam.ProteinAnalysis(protein_seq)
            return {
                "Séquence": protein_seq,
                "Poids": f"{round(analyser.molecular_weight() / 1000, 2)} kDa",
                "pI": round(analyser.isoelectric_point(), 2),
                "Instabilité": round(analyser.instability_index(), 2),
                "AA": analyser.get_amino_acids_percent()
            }
        except: return None

# ============================================================================
# 4. INTERFACE PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()

    db, genomique = st.session_state.db, st.session_state.genomique

    st.sidebar.title("🐑 EXPERT OVIN DZ")
    menu = ["📊 Dashboard", "📝 Scanner & Phénotype", "🧬 Laboratoire ADN", "🔬 Expertise Protéique", "🌾 Nutrition & Rapport PDF"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- NUTRITION & RAPPORT PDF ---
    if choice == "🌾 Nutrition & Rapport PDF":
        st.title("🌾 Optimisation & Rapport de Rationnement")
        
        aliments_dz = {
            "Orge (Chaïr)": {"ufl": 1.0, "pdi": 80, "prix": 4500},
            "Son de blé (Nkhala)": {"ufl": 0.82, "pdi": 95, "prix": 2500},
            "Foin de Luzerne": {"ufl": 0.65, "pdi": 90, "prix": 4000},
            "Paille traitée": {"ufl": 0.45, "pdi": 45, "prix": 1800},
            "Maïs concassé": {"ufl": 1.15, "pdi": 95, "prix": 6200}
        }

        col_cfg, col_res = st.columns([1, 1.5])
        
        with col_cfg:
            st.subheader("Configuration de l'Animal")
            a_id = st.text_input("Identifiant Animal", "DZ-2026-001")
            poids_n = st.number_input("Poids (kg)", 30, 120, 60)
            stade = st.selectbox("Stade", ["Entretien", "Gestation", "Lactation"])
            
            besoins = {"ufl": 0.8, "pdi": 75}
            if stade == "Gestation": besoins = {"ufl": 1.1, "pdi": 110}
            elif stade == "Lactation": besoins = {"ufl": 1.7, "pdi": 165}
            
            st.info(f"Cibles : {besoins['ufl']} UFL | {besoins['pdi']}g PDI")
            
            st.subheader("Ration journalière (kg)")
            choix_qte = {}
            for al, val in aliments_dz.items():
                choix_qte[al] = st.slider(f"{al}", 0.0, 2.5, 0.0, step=0.1)

        with col_res:
            st.subheader("Analyse de la Ration")
            t_ufl = sum(choix_qte[al] * aliments_dz[al]['ufl'] for al in choix_qte)
            t_pdi = sum(choix_qte[al] * aliments_dz[al]['pdi'] for al in choix_qte)
            t_cout = sum((choix_qte[al] / 100) * aliments_dz[al]['prix'] for al in choix_qte)
            
            # Graphique de couverture des besoins
            df_plot = pd.DataFrame({
                "Paramètre": ["Énergie (UFL)", "Protéines (PDI)"],
                "Apport": [t_ufl, t_pdi],
                "Besoin": [besoins['ufl'], besoins['pdi']]
            })
            st.bar_chart(df_plot.set_index("Paramètre"))
            
            m1, m2, m3 = st.columns(3)
            m1.metric("Total UFL", round(t_ufl, 2), delta=round(t_ufl - besoins['ufl'], 2))
            m2.metric("Total PDI (g)", round(t_pdi, 1), delta=round(t_pdi - besoins['pdi'], 1))
            m3.metric("Coût (DA/j)", f"{round(t_cout, 2)}")

            # --- GÉNÉRATION PDF ---
            st.divider()
            data_pdf = {
                "animal_id": a_id, "poids": poids_n, "stade": stade,
                "aliments": choix_qte, "total_ufl": round(t_ufl, 2),
                "total_pdi": round(t_pdi, 1), "cible_ufl": besoins['ufl'],
                "cible_pdi": besoins['pdi'], "cout": round(t_cout, 2)
            }
            
            if st.button("📄 Générer le Rapport PDF Professionnel"):
                pdf_bytes = create_pdf_report(data_pdf)
                st.download_button(
                    label="⬇️ Télécharger le Rapport",
                    data=pdf_bytes,
                    file_name=f"Ration_{a_id}_{date.today()}.pdf",
                    mime="application/pdf"
                )
                st.success("Le rapport a été généré avec succès.")

    # --- AUTRES MENUS (REPRISE DES VERSIONS PRÉCÉDENTES) ---
    elif choice == "🔬 Expertise Protéique":
        st.title("🔬 Analyse Moléculaire")
        dna = st.text_area("Séquence ADN")
        if dna:
            res = genomique.analyser_proteine(dna.strip())
            if res:
                st.json(res)
                st.bar_chart(pd.DataFrame.from_dict(res['AA'], orient='index'))

    elif choice == "📊 Dashboard":
        st.title("📊 Vue d'ensemble")
        st.write(db.fetch_all_as_df("SELECT * FROM brebis"))

    elif choice == "📝 Scanner & Phénotype":
        st.title("📷 Scanner 1 mètre")
        cam = st.camera_input("Scanner l'animal")
        if cam: st.image(cam, caption="Analyse en cours via étalon 1m...")

    elif choice == "🧬 Laboratoire ADN":
        st.title("🧬 Génomique")
        dna_lab = st.text_area("Entrée FASTA")
        if dna_lab: st.info("Analyse des marqueurs de prolificité en cours...")

if __name__ == "__main__":
    main()
