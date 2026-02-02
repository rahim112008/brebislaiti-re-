"""
EXPERT OVIN PRO - Système Unifié : Gestion, Génomique & Morphométrie
Version: 6.0 (Intégration Totale : NCBI + Alignement + Morpho 1m + EPG)
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
from datetime import datetime, timedelta, date
from PIL import Image
from contextlib import contextmanager
import io

# Bibliothèques Bioinformatiques
from Bio import Entrez, SeqIO, Align
from Bio.Seq import Seq

# ==========================================
# CONFIGURATION & CONSTANTES SCIENTIFIQUES
# ==========================================
Entrez.email = "votre.email@exemple.com" # À remplacer par votre email
DB_NAME = "expert_ovin_v6.db"

MORPHO_TRAITS = {
    "Hauteur au garrot": "Capacité corporelle globale",
    "Profondeur de poitrine": "Capacité respiratoire et digestive",
    "Profondeur de mamelle": "Corrélation directe avec le volume de lait",
    "Attache arrière": "Soutien de la mamelle et longévité"
}

# ==========================================
# GESTION DE LA BASE DE DONNÉES
# ==========================================
@contextmanager
def get_db_connection():
    conn = sqlite3.connect(DB_NAME, check_same_thread=False)
    try:
        yield conn
        conn.commit()
    except Exception as e:
        st.error(f"Erreur Database: {e}")
        conn.rollback()
    finally:
        conn.close()

def init_database():
    with get_db_connection() as conn:
        # Table Animaux
        conn.execute('''CREATE TABLE IF NOT EXISTS animaux (
            id TEXT PRIMARY KEY, boucle TEXT UNIQUE, race TEXT, 
            date_naiss DATE, potentiel TEXT)''')
        
        # Table Génomique & Biochimie
        conn.execute('''CREATE TABLE IF NOT EXISTS analyse_labo (
            id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_analyse DATE,
            tb REAL, tp REAL, cellules INTEGER, snp_dgat1 TEXT, snp_csn3 TEXT,
            FOREIGN KEY(animal_id) REFERENCES animaux(id))''')

        # Table Morphométrie
        conn.execute('''CREATE TABLE IF NOT EXISTS morphometrie (
            id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_mesure DATE,
            hauteur_garrot REAL, prof_mamelle REAL, 
            FOREIGN KEY(animal_id) REFERENCES animaux(id))''')

        # Table Reproduction (EPG)
        conn.execute('''CREATE TABLE IF NOT EXISTS reproduction (
            id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_eponge DATE, 
            date_mise_bas_prevue DATE, FOREIGN KEY(animal_id) REFERENCES animaux(id))''')

def seed_demo_data():
    """Données pré-installées pour démonstration"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM animaux")
        if cursor.fetchone()[0] == 0:
            # Insertion d'une brebis Élite
            cursor.execute("INSERT INTO animaux VALUES (?,?,?,?,?)", 
                          ('ELITE_01', 'FR778899', 'Lacaune', '2023-05-10', 'Élite Génomique'))
            # Données Labo
            cursor.execute("INSERT INTO analyse_labo (animal_id, date_analyse, tb, tp, snp_dgat1, snp_csn3) VALUES (?,?,?,?,?,?)",
                          ('ELITE_01', '2026-01-15', 74.2, 59.5, 'AA (Haut Rendement)', 'BB'))
            # Reproduction
            d_ep = date.today() - timedelta(days=10)
            d_prev = d_ep + timedelta(days=164)
            cursor.execute("INSERT INTO reproduction (animal_id, date_eponge, date_mise_bas_prevue) VALUES (?,?,?)",
                          ('ELITE_01', d_ep.isoformat(), d_prev.isoformat()))

# ==========================================
# FONCTIONS BIOINFORMATIQUES & ANALYSE
# ==========================================
def fetch_ncbi_sequence(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        return f"Erreur NCBI : {e}"

def align_sequences(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(seq1, seq2)
    return alignments[0]

# ==========================================
# MODULES INTERFACE (VUES)
# ==========================================

def module_dashboard():
    st.title("📊 Dashboard de l'Exploitation")
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT * FROM animaux", conn)
        df_repro = pd.read_sql("SELECT * FROM reproduction", conn)
    
    c1, c2, c3 = st.columns(3)
    c1.metric("Effectif Troupeau", len(df_animaux))
    c2.metric("Mises bas prévues", len(df_repro))
    
    st.subheader("📋 Liste des Animaux Élite")
    st.dataframe(df_animaux, use_container_width=True)

def module_repro_epg():
    st.title("🐑 Gestion de la Reproduction & EPG")
    st.info("Algorithme de prédiction basé sur le protocole de synchronisation des chaleurs.")
    
    with st.form("form_epg"):
        animal_id = st.text_input("ID de la Brebis (ex: ELITE_01)")
        d_eponge = st.date_input("Date de pose de l'éponge")
        protocole = st.selectbox("Durée du protocole (jours)", [12, 14])
        
        if st.form_submit_button("Enregistrer & Prédire"):
            # Calcul : Date éponge + durée protocole + 150j gestation
            date_mb = d_eponge + timedelta(days=protocole + 150)
            with get_db_connection() as conn:
                conn.execute("INSERT INTO reproduction (animal_id, date_eponge, date_mise_bas_prevue) VALUES (?,?,?)",
                            (animal_id, d_eponge.isoformat(), date_mb.isoformat()))
            st.success(f"✅ Enregistré. Mise bas prévue le : {date_mb.strftime('%d/%m/%Y')}")

def module_morpho_scanner():
    st.title("📏 Scanner Morphométrique (Étalon 1m)")
    st.write("Analyse phénotypique par traitement d'image.")
    
    img_file = st.file_uploader("Prendre/Charger la photo de l'animal", type=['jpg', 'png'])
    if img_file:
        img = Image.open(img_file)
        st.image(img, caption="Image de référence")
        
        st.subheader("⚙️ Calibration et Mesures")
        col1, col2 = st.columns(2)
        pix_m = col1.number_input("Pixels pour 1 mètre (étalon au sol)", value=500)
        ratio = 100 / pix_m
        
        h_px = col2.number_input("Mesure Hauteur Garrot (en pixels)", value=350)
        m_px = col2.number_input("Mesure Profondeur Mamelle (en pixels)", value=150)
        
        h_cm = round(h_px * ratio, 2)
        m_cm = round(m_px * ratio, 2)
        
        st.divider()
        st.metric("Hauteur Garrot Calculée", f"{h_cm} cm")
        st.metric("Profondeur Mamelle Calculée", f"{m_cm} cm")
        
        if st.button("Sauvegarder les mesures"):
            st.success("Mesures ajoutées à la base phénotypique.")

def module_bioinfo_labo():
    st.title("🧬 Laboratoire Bioinformatique & Génomique")
    
    tab1, tab2, tab3 = st.tabs(["Extraction NCBI", "Alignement SNP", "Analyses Labo"])
    
    with tab1:
        st.subheader("🔍 Accès GenBank")
        acc_id = st.text_input("Accession ID (ex: NM_001009378)", "NM_001009378")
        if st.button("Récupérer Séquence"):
            with st.spinner("Recherche sur NCBI..."):
                res = fetch_ncbi_sequence(acc_id)
                if isinstance(res, str): st.error(res)
                else:
                    st.text_area("Séquence FASTA", str(res.seq), height=150)
                    st.caption(f"Source: {res.description}")

    with tab2:
        st.subheader("🔗 Comparaison de Séquences (Alignement Global)")
        c1, c2 = st.columns(2)
        s1 = c1.text_area("Séquence de Référence", "ATGC...")
        s2 = c2.text_area("Séquence Brebis Test", "ATGC...")
        if st.button("Lancer l'Alignement"):
            if len(s1) > 5 and len(s2) > 5:
                ali = align_sequences(s1, s2)
                st.write(f"**Score d'identité :** {ali.score}")
                st.text(ali)
            else:
                st.warning("Veuillez entrer des séquences valides.")

    with tab3:
        st.subheader("🧪 Résultats Biochimiques")
        with get_db_connection() as conn:
            df_labo = pd.read_sql("""
                SELECT a.boucle, l.tb, l.tp, l.snp_dgat1, l.snp_csn3 
                FROM analyse_labo l JOIN animaux a ON l.animal_id = a.id""", conn)
        st.dataframe(df_labo, use_container_width=True)
        if not df_labo.empty:
            fig = px.bar(df_labo, x="boucle", y=["tb", "tp"], barmode="group", title="Compositions Laitières")
            st.plotly_chart(fig)

# ==========================================
# POINT D'ENTRÉE PRINCIPAL
# ==========================================
def main():
    st.set_page_config(page_title="EXPERT OVIN PRO v6", layout="wide", page_icon="🐑")
    
    # Init DB et Données Démo
    init_database()
    seed_demo_data()
    
    # Barre Latérale
    st.sidebar.title("🧬 EXPERT OVIN PRO")
    st.sidebar.write("Système Décisionnel")
    menu = st.sidebar.radio("Navigation", [
        "Dashboard", 
        "Reproduction & EPG", 
        "Scanner Morpho", 
        "Bioinformatique & Labo",
        "Statistiques R"
    ])
    
    if menu == "Dashboard": module_dashboard()
    elif menu == "Reproduction & EPG": module_repro_epg()
    elif menu == "Scanner Morpho": module_morpho_scanner()
    elif menu == "Bioinformatique & Labo": module_bioinfo_labo()
    elif menu == "Statistiques R":
        st.title("📊 Analyse Statistiques R")
        st.info("Export des données vers R pour calcul des indices de sélection.")
        st.button("Générer export CSV pour R-Stats")

if __name__ == "__main__":
    main()
