"""
OVINMASTER MASTER-V50 : ULTIMATE GENOMIC & BUSINESS EDITION
-----------------------------------------------------------
Système Expert Intégré pour la Filière Ovine Algérienne.
Dépendances : streamlit, pandas, biopython, plotly
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import io
import sqlite3
from datetime import datetime, date

# ============================================================================
# 1. DATABASE & SECURITY
# ============================================================================

class Database:
    def __init__(self, db_path: str = "ovin_v50_master.db"):
        self.conn = sqlite3.connect(db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self.create_tables()

    def create_tables(self):
        sqls = [
            "CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, role TEXT, region TEXT)",
            "CREATE TABLE IF NOT EXISTS brebis (id_u TEXT PRIMARY KEY, owner TEXT, race TEXT, poids REAL, created_at DATE)",
            "CREATE TABLE IF NOT EXISTS sante (id_u TEXT, acte TEXT, produit TEXT, date_acte DATE)",
            "CREATE TABLE IF NOT EXISTS stocks (owner TEXT, aliment TEXT, qte REAL, PRIMARY KEY(owner, aliment))"
        ]
        for s in sqls: self.conn.execute(s)
        self.conn.execute("INSERT OR IGNORE INTO users VALUES ('admin', 'masterdz', 'Expert', 'Alger')")
        self.conn.execute("INSERT OR IGNORE INTO users VALUES ('eleveur', 'dz2026', 'Eleveur', 'Sétif')")
        self.conn.commit()

db = Database()

# ============================================================================
# 2. SCIENTIFIC & GENOMIC ENGINES
# ============================================================================

@st.cache_data
def run_genomic_pipeline(fasta_str):
    """Analyse Biopython des séquences ADN"""
    fasta_io = io.StringIO(fasta_str)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    results = []
    for rec in records:
        gc = gc_fraction(rec.seq) * 100
        results.append({
            "Sujet": rec.id,
            "Longeur_pb": len(rec.seq),
            "GC_Percent": round(gc, 2),
            "Sequence": str(rec.seq[:50]) + "..."
        })
    return pd.DataFrame(results)

def get_market_metrics(poids, bcs):
    """Calcul Business & Rendement"""
    rendement = round(42 + (max(0, bcs - 2.5) * 4), 1)
    valeur_est = round(poids * 1250 * (1 + (bcs - 3) * 0.1), 0)
    return rendement, valeur_est

# ============================================================================
# 3. USER INTERFACE (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="OvinMaster Pro V50", layout="wide", page_icon="🧬")
    
    # Custom CSS for Mobile
    st.markdown("<style>.stButton>button {width:100%; border-radius:10px;}</style>", unsafe_allow_html=True)

    if 'auth' not in st.session_state:
        st.title("🐑 OvinMaster Pro V50")
        st.subheader("Station de Génomique & Biométrie Algérienne")
        col1, col2 = st.columns(2)
        u = col1.text_input("Login")
        p = col2.text_input("Pass", type="password")
        if st.button("🚀 Accéder au Système"):
            res = db.conn.execute("SELECT * FROM users WHERE username=? AND password=?", (u,p)).fetchone()
            if res:
                st.session_state.auth, st.session_state.user, st.session_state.role = True, res['username'], res['role']
                st.rerun()
        return

    # --- SIDEBAR & NAV ---
    st.sidebar.title(f"👤 {st.session_state.user}")
    mode = st.sidebar.radio("Navigation", ["📊 Dashboard", "📸 Scanner & BCS", "🧬 Labo ADN", "💉 Santé & Soins", "🍲 Stocks"])

    # --- MODULES ---
    
    if mode == "📊 Dashboard":
        st.title("📊 Tableau de Bord Stratégique")
        df = pd.read_sql(f"SELECT * FROM brebis WHERE owner='{st.session_state.user}'", db.conn)
        if not df.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df))
            c2.metric("Poids Moyen (kg)", round(df['poids'].mean(), 1))
            c3.metric("Race Dominante", df['race'].mode()[0])
            st.plotly_chart(px.sunburst(df, path=['race', 'id_u'], values='poids', title="Structure du Troupeau"))
            
        else:
            st.info("Aucun animal enregistré. Commencez par le Scanner.")

    elif mode == "📸 Scanner & BCS":
        st.title("📸 Scanner Biométrique IA")
        with st.form("scan_form"):
            col1, col2 = st.columns(2)
            id_u = col1.text_input("ID Animal (Boucle)")
            race = col1.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra"])
            poids = col2.number_input("Poids (kg)", 5.0, 150.0, 55.0)
            bcs = col2.select_slider("Score BCS (1-5)", options=[1, 2, 3, 4, 5], value=3)
            
            # Paramètres 1 mètre
            st.divider()
            h_garrot = st.slider("Hauteur au Garrot (cm) - Réf 1m", 40, 110, 75)
            l_corps = st.slider("Longueur de Corps (cm)", 40, 130, 85)
            
            if st.form_submit_button("🎯 Analyser & Sauvegarder"):
                db.conn.execute("INSERT OR REPLACE INTO brebis VALUES (?,?,?,?,?)", (id_u, st.session_state.user, race, poids, date.today()))
                db.conn.commit()
                rend, val = get_market_metrics(poids, bcs)
                st.success(f"Analyse terminée ! Rendement : {rend}% | Valeur : {val:,} DZD")
                
                # Radar Chart
                fig = go.Figure(data=go.Scatterpolar(r=[h_garrot, l_corps, bcs*20, 80], theta=['Taille','Format','Muscle','Aplombs'], fill='toself'))
                st.plotly_chart(fig)
                
                

    elif mode == "🧬 Labo ADN":
        st.title("🧬 Laboratoire de Génomique In Silico")
        st.write("Analyseur de séquences basé sur Biopython.")
        up = st.file_uploader("Charger fichier FASTA", type=["fasta", "fa"])
        if up:
            fasta_str = up.getvalue().decode("utf-8")
            df_gen = run_genomic_pipeline(fasta_str)
            st.dataframe(df_gen, use_container_width=True)
            st.plotly_chart(px.bar(df_gen, x="Sujet", y="GC_Percent", title="Stabilité Génomique (GC%)"))
            st.info("💡 Note : Un taux de GC > 45% indique une forte adaptation thermique chez l'ovin.")
            

    elif mode == "💉 Santé & Soins":
        st.title("💉 Carnet de Santé Numérique")
        
        # Interface de saisie simplifiée ici...

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.clear()
        st.rerun()

if __name__ == "__main__":
    main()
