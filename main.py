import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta
from PIL import Image
from contextlib import contextmanager
from Bio import Entrez, SeqIO, Align, SeqUtils
from scipy import stats
import io

# ==========================================
# CONFIGURATION SCIENTIFIQUE
# ==========================================
st.set_page_config(page_title="OVIN-GENOMICS PRO", layout="wide", page_icon="🧬")

Entrez.email = "votre.recherche@institut.com"
DB_NAME = "ovin_research_v9.db"

# Pondérations pour l'Index de Sélection (Modifiables par le généticien)
WEIGHTS = {"TB": 0.45, "TP": 0.35, "MORPHO": 0.20}

# ==========================================
# GESTION DES DONNÉES (SQL RELATIONNEL)
# ==========================================
@contextmanager
def get_db_connection():
    conn = sqlite3.connect(DB_NAME, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()

def init_database():
    with get_db_connection() as conn:
        conn.executescript('''
            CREATE TABLE IF NOT EXISTS individus (
                id TEXT PRIMARY KEY, boucle TEXT UNIQUE, race TEXT, date_naiss DATE, sexe TEXT);
            
            CREATE TABLE IF NOT EXISTS lab_genomics (
                id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_analyse DATE,
                tb REAL, tp REAL, dgat1_genotype TEXT, csn3_genotype TEXT, gc_content REAL,
                FOREIGN KEY(animal_id) REFERENCES individus(id));

            CREATE TABLE IF NOT EXISTS phenotypage (
                id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_mesure DATE,
                hg REAL, pm REAL, ic_score REAL,
                FOREIGN KEY(animal_id) REFERENCES individus(id));
        ''')

# ==========================================
# MODULE BIOINFORMATIQUE (Expert)
# ==========================================
class BioInfoEngine:
    @staticmethod
    def analyze_sequence(record):
        """Calcule les paramètres biophysiques de la séquence récupérée"""
        details = {
            "ID": record.id,
            "Description": record.description,
            "Longueur": len(record.seq),
            "GC%": round(SeqUtils.gc_fraction(record.seq) * 100, 2),
            "Masse Moléculaire": round(SeqUtils.molecular_weight(record.seq, seq_type="DNA"), 2)
        }
        return details

    @staticmethod
    def find_motifs(sequence, motif):
        """Recherche de sites de restriction ou d'amorces"""
        return [m.start() for m in SeqUtils.nt_search(str(sequence), motif)[1:]]

# ==========================================
# VUES INTERFACE (UI)
# ==========================================

def view_genomics_lab():
    st.title("🧬 Laboratoire de Génomique Moléculaire")
    
    tab1, tab2 = st.tabs(["Extraction & Analyse Séquence", "Alignement de Variants"])
    
    with tab1:
        acc_id = st.text_input("ID Accession NCBI (ex: NM_001009378.1)", "NM_001009378.1")
        motif_search = st.text_input("Rechercher un motif (ex: TATA, CCGG)", "ATGC")
        
        if st.button("Lancer l'analyse biophysique"):
            with st.spinner("Accès aux serveurs NCBI..."):
                try:
                    with Entrez.efetch(db="nucleotide", id=acc_id, rettype="fasta", retmode="text") as h:
                        record = SeqIO.read(h, "fasta")
                        stats = BioInfoEngine.analyze_sequence(record)
                        
                        # Affichage des métriques expertes
                        c1, c2, c3, c4 = st.columns(4)
                        c1.metric("GC Content", f"{stats['GC%']}%")
                        c2.metric("Longueur (pb)", stats['Longueur'])
                        c3.metric("Masse (Da)", stats['Masse Moléculaire'])
                        
                        motifs_found = BioInfoEngine.find_motifs(record.seq, motif_search)
                        c4.metric(f"Motifs '{motif_search}'", len(motifs_found))
                        
                        st.text_area("Séquence FASTA complète", str(record.seq), height=250)
                except Exception as e:
                    st.error(f"Erreur d'accès : {e}")

def view_genetic_selection():
    st.title("📊 Génétique Quantitative & Sélection")
    
    with get_db_connection() as conn:
        df = pd.read_sql('''
            SELECT i.boucle, l.tb, l.tp, p.hg, p.pm 
            FROM individus i 
            JOIN lab_genomics l ON i.id = l.animal_id 
            JOIN phenotypage p ON i.id = p.animal_id''', conn)
    
    if not df.empty:
        # Calcul de l'Index de Sélection Scientifique
        # Normalisation Z-score pour comparer des unités différentes
        for col in ['tb', 'tp', 'hg', 'pm']:
            df[f'z_{col}'] = (df[col] - df[col].mean()) / df[col].std()
        
        df['Selection_Index'] = (df['z_tb'] * WEIGHTS['TB']) + \
                                (df['z_tp'] * WEIGHTS['TP']) + \
                                (df['z_hg'] * WEIGHTS['MORPHO'])
        
        st.subheader("🏆 Classement par Valeur Génétique Estimée (EBV)")
        st.dataframe(df.sort_values(by='Selection_Index', ascending=False), use_container_width=True)
        
        # Distribution de la population
        fig = px.histogram(df, x="Selection_Index", nbins=20, title="Distribution de l'Index de Sélection du Troupeau",
                           marginal="rug", color_discrete_sequence=['#2e7d32'])
        st.plotly_chart(fig)
        
        

def main():
    # Initialisation DB (Automatique)
    init_database()
    
    st.sidebar.title("🔬 OVIN-GENOMICS V9")
    st.sidebar.markdown("---")
    menu = st.sidebar.selectbox("Expertise", 
        ["Génomique Moléculaire", "Génétique Quantitative", "Morphométrie IA", "Base de Données"])
    
    if menu == "Génomique Moléculaire": view_genomics_lab()
    elif menu == "Génétique Quantitative": view_genetic_selection()
    elif menu == "Base de Données":
        st.title("📂 Gestion du Registre")
        # Logique d'ajout d'animaux ici...

if __name__ == "__main__":
    main()
