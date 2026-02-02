import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta, date
from PIL import Image
from contextlib import contextmanager
from Bio import Entrez, SeqIO, Align
import scipy.stats as stats
import io

# ==========================================
# CONFIGURATION & DESIGN UI
# ==========================================
st.set_page_config(
    page_title="EXPERT OVIN PRO | Intelligence Zootechnique",
    page_icon="🧬",
    layout="wide"
)

# Style CSS personnalisé pour un rendu professionnel
st.markdown("""
    <style>
    .main { background-color: #f8f9fa; }
    .stMetric { background-color: #ffffff; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
    .sidebar .sidebar-content { background-image: linear-gradient(#2e7d32, #1b5e20); color: white; }
    h1, h2, h3 { color: #1b5e20; font-family: 'Segoe UI', sans-serif; }
    </style>
    """, unsafe_allow_html=True)

# Configuration API NCBI
Entrez.email = "votre.email@expert-ovin.com" 
DB_NAME = "expert_ovin_v7.db"

# ==========================================
# GESTION DES DONNÉES (DAL - Data Access Layer)
# ==========================================
@contextmanager
def get_db_connection():
    conn = sqlite3.connect(DB_NAME, check_same_thread=False)
    conn.row_factory = sqlite3.Row # Accès aux colonnes par nom
    try:
        yield conn
        conn.commit()
    except Exception as e:
        st.error(f"⚠️ Erreur Critique Database : {e}")
        conn.rollback()
    finally:
        conn.close()

def init_database():
    """Schéma de base de données relationnel complet"""
    with get_db_connection() as conn:
        conn.executescript('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY, boucle TEXT UNIQUE, race TEXT, date_naiss DATE, potentiel TEXT);
            
            CREATE TABLE IF NOT EXISTS analyse_labo (
                id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_analyse DATE,
                tb REAL, tp REAL, cellules INTEGER, snp_dgat1 TEXT, snp_csn3 TEXT,
                FOREIGN KEY(animal_id) REFERENCES animaux(id) ON DELETE CASCADE);

            CREATE TABLE IF NOT EXISTS morphometrie (
                id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_mesure DATE,
                hauteur_garrot REAL, prof_mamelle REAL, prof_poitrine REAL,
                FOREIGN KEY(animal_id) REFERENCES animaux(id) ON DELETE CASCADE);

            CREATE TABLE IF NOT EXISTS reproduction (
                id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, date_eponge DATE, 
                date_mise_bas_prevue DATE, statut TEXT DEFAULT 'En attente',
                FOREIGN KEY(animal_id) REFERENCES animaux(id) ON DELETE CASCADE);
        ''')

def seed_demo_data():
    """Injection de données de recherche pour test immédiat"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM animaux")
        if cursor.fetchone()[0] == 0:
            cursor.execute("INSERT INTO animaux VALUES (?,?,?,?,?)", 
                          ('ELITE_01', 'FR778899', 'Lacaune', '2023-05-10', 'Élite Génomique'))
            cursor.execute("INSERT INTO analyse_labo (animal_id, date_analyse, tb, tp, snp_dgat1, snp_csn3) VALUES (?,?,?,?,?,?)",
                          ('ELITE_01', '2026-01-15', 74.2, 59.5, 'AA (Haut Rendement)', 'BB'))
            cursor.execute("INSERT INTO morphometrie (animal_id, date_mesure, hauteur_garrot, prof_mamelle) VALUES (?,?,?,?)",
                          ('ELITE_01', '2026-01-20', 72.5, 28.4))

# ==========================================
# MODULES SCIENTIFIQUES & BIOINFO
# ==========================================
def fetch_ncbi_sequence(accession_id):
    try:
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
            return SeqIO.read(handle, "fasta")
    except Exception as e:
        return f"Erreur NCBI : {str(e)}"

# ==========================================
# INTERFACE UTILISATEUR (UI)
# ==========================================

def view_dashboard():
    st.title("📊 Dashboard Analytique")
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT * FROM animaux", conn)
        df_repro = pd.read_sql("SELECT * FROM reproduction WHERE statut='En attente'", conn)
    
    m1, m2, m3 = st.columns(3)
    m1.metric("Effectif Total", len(df_animaux))
    m2.metric("Gestation en cours", len(df_repro))
    m3.metric("Indice de Sélection Moyen", "112.5") # Placeholder

    st.divider()
    st.subheader("📋 Registre des Brebis Élites")
    st.dataframe(df_animaux, use_container_width=True)

def view_reproduction():
    st.title("🐑 Planning de Reproduction & EPG")
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("Nouveau Protocole")
        with st.form("repro_form"):
            a_id = st.text_input("ID Animal")
            d_ep = st.date_input("Date de pose éponge")
            prot = st.selectbox("Protocole", [12, 14])
            if st.form_submit_button("Calculer Prédiction"):
                d_mb = d_ep + timedelta(days=prot + 150)
                with get_db_connection() as conn:
                    conn.execute("INSERT INTO reproduction (animal_id, date_eponge, date_mise_bas_prevue) VALUES (?,?,?)",
                                (a_id, d_ep, d_mb))
                st.success(f"📅 Mise bas : {d_mb.strftime('%d/%m/%Y')}")
    
    with col2:
        st.subheader("Calendrier des mises bas")
        with get_db_connection() as conn:
            df = pd.read_sql("SELECT * FROM reproduction", conn)
            st.table(df)

def view_morphometrie():
    st.title("📏 Phénotypage par Image (Standard 1m)")
    
    
    uploaded_file = st.file_uploader("Prendre/Charger photo", type=['jpg', 'png'])
    if uploaded_file:
        img = Image.open(uploaded_file)
        st.image(img, caption="Analyse morphométrique v2.0")
        
        c1, c2 = st.columns(2)
        pix_etalon = c1.number_input("Pixels pour 1 mètre (Étalon au sol)", value=500)
        ratio = 100 / pix_etalon
        
        h_px = c2.number_input("Mesure Hauteur Garrot (Pixels)", value=350)
        st.info(f"📏 Taille réelle : **{round(h_px * ratio, 2)} cm**")

def view_bioinfo():
    st.title("🧬 Laboratoire Bioinformatique")
    tab1, tab2 = st.tabs(["Séquençage NCBI", "Alignement SNP"])
    
    with tab1:
        acc_id = st.text_input("Accession ID GenBank", "NM_001009378")
        if st.button("Analyser Gène"):
            record = fetch_ncbi_sequence(acc_id)
            if hasattr(record, 'seq'):
                st.success(f"Gène : {record.description}")
                st.text_area("Séquence", str(record.seq[:500]) + "...", height=150)
            else: st.error(record)

def view_statistics():
    st.title("📊 Intelligence Statistique & R-Stats")
    
    with get_db_connection() as conn:
        query = """
            SELECT a.boucle, m.hauteur_garrot, m.prof_mamelle, l.tb, l.tp
            FROM animaux a
            LEFT JOIN morphometrie m ON a.id = m.animal_id
            LEFT JOIN analyse_labo l ON a.id = l.animal_id
        """
        df = pd.read_sql(query, conn)

    if not df.dropna().empty:
        # Graphique Radar : Comparaison Élite vs Moyenne
        st.subheader("🕸️ Profil Morpho-Laitier (Comparateur Radar)")
        
        # Exemple de radar chart pour la première brebis
        categories = ['Hauteur Garrot', 'Prof. Mamelle', 'Taux Butyreux', 'Taux Protéique']
        values = [df.iloc[0]['hauteur_garrot']/80, df.iloc[0]['prof_mamelle']/40, 
                  df.iloc[0]['tb']/80, df.iloc[0]['tp']/60]
        
        fig = go.Figure()
        fig.add_trace(go.Scatterpolar(r=values, theta=categories, fill='toself', name=df.iloc[0]['boucle']))
        fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 1])), showlegend=True)
        st.plotly_chart(fig)

        # Export R
        st.divider()
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("📥 Télécharger CSV pour R-Stats", csv, "export_r.csv", "text/csv")
    else:
        st.warning("Données insuffisantes pour les statistiques avancées.")

# ==========================================
# MAIN APP FLOW
# ==========================================
def main():
    init_database()
    seed_demo_data()
    
    st.sidebar.title("🐑 EXPERT OVIN PRO")
    st.sidebar.caption("Système Expert de Sélection v7.0")
    
    menu = st.sidebar.radio("Navigation", 
        ["Dashboard", "Reproduction & EPG", "Scanner Morpho", "Bioinformatique", "Statistiques R"])
    
    if menu == "Dashboard": view_dashboard()
    elif menu == "Reproduction & EPG": view_reproduction()
    elif menu == "Scanner Morpho": view_morphometrie()
    elif menu == "Bioinformatique": view_bioinfo()
    elif menu == "Statistiques R": view_statistics()

if __name__ == "__main__":
    main()
