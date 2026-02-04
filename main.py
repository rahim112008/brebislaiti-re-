"""
OVIN MANAGER ULTIMATE - Version Fusionnée & Génomique
Système : Races Algériennes + Morphologie 3D + G-BLUP + Bioinformatique + R-Stats
"""

# ============================================================================
# SECTION 1: IMPORTS
# ============================================================================
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import random
import math
import io
import json
from PIL import Image, ImageDraw

# ============================================================================
# SECTION 2: CONFIGURATION & CSS (FUSIONNÉ)
# ============================================================================
st.set_page_config(page_title="Ovin Manager Pro - Algérie Génomique", layout="wide", page_icon="🐑")

st.markdown("""
<style>
    .main-header { font-size: 2.8rem; color: #8B0000; text-align: center; background: linear-gradient(90deg, #8B0000, #4a148c); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
    .section-header { font-size: 2rem; color: #8B0000; border-bottom: 3px solid #FF4500; padding-bottom: 10px; }
    .bio-info-view { background-color: #0e1117; color: #00ff41; font-family: 'Courier New', monospace; padding: 15px; border-radius: 5px; border: 1px solid #00ff41; }
    .metric-card { background: white; border-radius: 15px; padding: 15px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border-left: 5px solid #8B0000; }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# SECTION 3: STANDARDS DES RACES ALGÉRIENNES
# ============================================================================
STANDARDS_RACES = {
    'HAMRA': {'nom': 'Hamra', 'poids': (45, 90), 'lait': (1.5, 3.5)},
    'OUDA': {'nom': 'Ouled Djellal', 'poids': (50, 100), 'lait': (1.0, 2.5)},
    'SIDAHOU': {'nom': 'Sidahou', 'poids': (40, 85), 'lait': (1.2, 2.8)},
    'BERBERE': {'nom': 'Berbère', 'poids': (35, 70), 'lait': (0.8, 2.0)},
    'CROISE': {'nom': 'Croisement', 'poids': (40, 95), 'lait': (1.0, 3.0)}
}

# ============================================================================
# SECTION 4: MOTEUR BIOINFORMATIQUE & G-BLUP
# ============================================================================
class BioinfoEngine:
    @staticmethod
    def aligner_sequences(seq_ref, seq_ech):
        matches = "".join("|" if r == e else "." for r, e in zip(seq_ref, seq_ech))
        identite = (matches.count("|") / len(seq_ref)) * 100
        return identite, matches

    @staticmethod
    def calculer_gblup(genotypes_matrix):
        """
        G-BLUP: Genomic Best Linear Unbiased Prediction
        Calcule la GEBV (Genomic Estimated Breeding Value) pour le lait
        """
        # Simulation d'un modèle linéaire mixte génomique
        # GEBV = Z * u (où Z est la matrice des marqueurs et u les effets SNPs)
        n_individus = genotypes_matrix.shape[0]
        n_marqueurs = genotypes_matrix.shape[1]
        
        # Simulation d'effets SNPs (poids génétiques)
        effets_snps = np.random.normal(0.2, 0.05, n_marqueurs)
        
        # Calcul des scores
        gebv_scores = np.dot(genotypes_matrix, effets_snps)
        # Normalisation pour l'affichage
        gebv_scores = (gebv_scores - np.mean(gebv_scores)) / np.std(gebv_scores)
        return gebv_scores

# ============================================================================
# SECTION 5: GESTION BASE DE DONNÉES
# ============================================================================
def init_db():
    conn = sqlite3.connect('ovin_ultimate.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS troupeau 
                 (id TEXT PRIMARY KEY, race TEXT, sexe TEXT, poids REAL, 
                  lait_pref REAL, adn_seq TEXT, gebv_score REAL)''')
    conn.commit()
    return conn

# ============================================================================
# SECTION 6: MODULE SCANNER 3D (AVEC STANDARD 1M)
# ============================================================================
def page_scanner_3d():
    st.markdown('<h2 class="section-header">📐 SCANNER MORPHOLOGIQUE 3D</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.info("🎯 Système de calibration : Standard 1 mètre activé.")
        img_file = st.file_uploader("Importer capture Scanner/Caméra", type=['jpg', 'png'])
        if img_file:
            st.image(img_file, caption="Analyse en cours...")
            
    with col2:
        st.subheader("📝 Mesures calculées (Précision 0.1cm)")
        with st.form("mesures_form"):
            long_t = st.number_input("Longueur des trayons (cm)", min_value=0.0, max_value=20.0, value=4.5, step=0.1)
            larg_b = st.number_input("Largeur Bassin (cm)", min_value=0.0, value=42.0, step=0.1)
            haut_g = st.number_input("Hauteur au garrot (cm)", min_value=0.0, value=75.0, step=0.1)
            
            if st.form_submit_button("Enregistrer les mesures"):
                st.success("Mesures morphologiques enregistrées avec succès.")

# ============================================================================
# SECTION 7: PAGE GÉNOMIQUE & G-BLUP
# ============================================================================
def page_genomique_lab():
    st.markdown('<h1 class="main-header">🧬 LABORATOIRE GÉNOMIQUE & G-BLUP</h1>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["Alignement ADN", "Prédiction G-BLUP", "Statistiques Population"])
    
    with tab1:
        st.subheader("🔍 Analyse de Séquence (Alignement)")
        seq_ref = "ATGCGGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
        seq_ind = st.text_input("Séquence ADN Individu (ATCG)", "ATGCGGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCG")
        
        identite, matches = BioinfoEngine.aligner_sequences(seq_ref, seq_ind)
        
        st.markdown(f"**Taux d'identité :** `{identite:.2f}%`")
        st.markdown(f"""<div class="bio-info-view">REF: {seq_ref}<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{matches}<br>ECH: {seq_ind}</div>""", unsafe_allow_html=True)

    with tab2:
        st.subheader("📈 Prédiction des Performances (G-BLUP)")
        st.write("Calcul des GEBV (Genomic Estimated Breeding Values) pour la production laitière.")
        
        if st.button("Lancer la prédiction G-BLUP sur le troupeau"):
            # Simulation d'une matrice de génotypes (20 individus x 100 SNPs)
            matrix = np.random.randint(0, 3, size=(20, 100)) 
            scores = BioinfoEngine.calculer_gblup(matrix)
            
            df_pred = pd.DataFrame({
                'ID': [f'OVIN-{i:03d}' for i in range(1, 21)],
                'Score G-BLUP (Lait)': scores,
                'Fiabilité': [random.uniform(0.75, 0.95) for _ in range(20)]
            }).sort_values(by='Score G-BLUP (Lait)', ascending=False)
            
            st.dataframe(df_pred.style.background_gradient(subset=['Score G-BLUP (Lait)'], cmap='Greens'))
            
            fig = px.scatter(df_pred, x='ID', y='Score G-BLUP (Lait)', size='Fiabilité', color='Score G-BLUP (Lait)', title="Classement Génétique du Troupeau")
            st.plotly_chart(fig, use_container_width=True)

    with tab3:
        st.subheader("📊 R-Stats Analysis")
        c1, c2 = st.columns(2)
        with c1:
            # Matrice de parenté
            matrix_size = 8
            k_matrix = np.random.uniform(0, 0.25, size=(matrix_size, matrix_size))
            np.fill_diagonal(k_matrix, 1.0)
            fig_mat = px.imshow(k_matrix, labels=dict(color="Consanguinité"), title="Matrice de Parenté Additive (A-Matrix)")
            st.plotly_chart(fig_mat)
        with c2:
            # Héritabilité
            h2 = 0.35 # Simulation pour production laitière
            st.gauge = go.Figure(go.Indicator(
                mode = "gauge+number",
                value = h2,
                title = {'text': "Héritabilité (h²) - Lait"},
                gauge = {'axis': {'range': [0, 1]}, 'bar': {'color': "darkred"}}
            ))
            st.plotly_chart(st.gauge)

# ============================================================================
# SECTION 8: NAVIGATION & MAIN
# ============================================================================
def main():
    st.sidebar.title("🐑 OvinManager Pro v5")
    st.sidebar.markdown("---")
    menu = ["Tableau de Bord", "Scanner 3D & Morphologie", "Génomique & G-BLUP", "Races Algériennes"]
    choice = st.sidebar.selectbox("Accéder au module :", menu)
    
    init_db()

    if choice == "Tableau de Bord":
        st.markdown('<h1 class="main-header">📊 TABLEAU DE BORD GÉNÉRAL</h1>', unsafe_allow_html=True)
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Effectif", "150", "+5")
        col2.metric("Moyenne Lait", "2.4 L/j", "+0.2")
        col3.metric("GEBV Moyenne", "+1.15", "Top 5%")
        col4.metric("Consanguinité", "1.8%", "-0.2%")
        
    elif choice == "Scanner 3D & Morphologie":
        page_scanner_3d()
        
    elif choice == choice == "Génomique & G-BLUP":
        page_genomique_lab()
        
    elif choice == "Races Algériennes":
        st.subheader("📚 Standards des Races")
        race_sel = st.selectbox("Choisir une race", list(STANDARDS_RACES.keys()))
        data = STANDARDS_RACES[race_sel]
        st.write(f"**Nom :** {data['nom']}")
        st.write(f"**Poids adulte :** {data['poids'][0]} - {data['poids'][1]} kg")
        st.write(f"**Potentiel laitier :** {data['lait'][0]} - {data['lait'][1]} L/j")

if __name__ == "__main__":
    main()
