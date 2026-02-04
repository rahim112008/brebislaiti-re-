"""
OVIN MANAGER PRO - Version Intégrale 2026
Système de Gestion, Scanner 3D (Standard 1m) et Analyse Génétique
Races Algériennes : Hamra, Ouda, Sidahou, Berbère
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import json
import random
import math
import io
from PIL import Image, ImageDraw
from scipy import stats

# ========== CONFIGURATION PAGE ==========

st.set_page_config(
    page_title="Ovin Manager Pro - Expert Algérie",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ========== CSS PERSONNALISÉ (STYLISATION AVANCÉE) ==========

st.markdown("""
<style>
    .main-header {
        font-size: 2.8rem;
        color: #8B0000;
        text-align: center;
        margin-bottom: 1rem;
        background: linear-gradient(90deg, #8B0000, #FF4500);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-weight: bold;
    }
    .metric-card {
        background: linear-gradient(135deg, #fdfbfb 0%, #ebedee 100%);
        padding: 20px;
        border-radius: 15px;
        box-shadow: 2px 4px 10px rgba(0,0,0,0.1);
        border-left: 5px solid #8B0000;
        text-align: center;
    }
    .scanner-container {
        background-color: #000000;
        color: #00FF00;
        padding: 20px;
        border-radius: 15px;
        font-family: 'Courier New', Courier, monospace;
        border: 2px solid #333;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #f0f2f6;
        border-radius: 10px 10px 0px 0px;
        gap: 1px;
        padding-top: 10px;
        padding-bottom: 10px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #8B0000 !important;
        color: white !important;
    }
</style>
""", unsafe_allow_html=True)

# ========== STANDARDS DES RACES ALGÉRIENNES ==========

STANDARDS_RACES = {
    'HAMRA': {
        'nom_complet': 'Hamra (Rousse)',
        'poids_adulte': {'femelle': (45, 65), 'male': (65, 90)},
        'mensurations': {'longueur_cm': (95, 125), 'hauteur_cm': (65, 85), 'tour_poitrine_cm': (95, 120)},
        'couleur_pref': '#8B4513'
    },
    'OUDA': {
        'nom_complet': 'Ouled Djellal (Ouda)',
        'poids_adulte': {'femelle': (50, 70), 'male': (70, 100)},
        'mensurations': {'longueur_cm': (100, 130), 'hauteur_cm': (70, 90), 'tour_poitrine_cm': (100, 130)},
        'couleur_pref': '#F5F5DC'
    },
    'SIDAHOU': {
        'nom_complet': 'Sidahou',
        'poids_adulte': {'femelle': (40, 60), 'male': (60, 85)},
        'mensurations': {'longueur_cm': (90, 120), 'hauteur_cm': (60, 80), 'tour_poitrine_cm': (90, 115)},
        'couleur_pref': '#2F4F4F'
    },
    'BERBERE': {
        'nom_complet': 'Brebis Berbère',
        'poids_adulte': {'femelle': (35, 50), 'male': (50, 70)},
        'mensurations': {'longueur_cm': (80, 110), 'hauteur_cm': (55, 75), 'tour_poitrine_cm': (85, 110)},
        'couleur_pref': '#A0522D'
    }
}

# ========== INITIALISATION BASE DE DONNÉES ==========

def init_db():
    conn = sqlite3.connect('ovin_manager_v2.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table Principale
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant TEXT UNIQUE,
            nom TEXT,
            race TEXT,
            sexe TEXT,
            poids FLOAT,
            date_naissance DATE,
            score_condition INTEGER,
            longueur_cm FLOAT,
            hauteur_cm FLOAT,
            statut TEXT DEFAULT 'Actif'
        )
    ''')
    
    # Table Génétique
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS genotypage (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            marqueur TEXT,
            genotype TEXT,
            trait_associe TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Peuplement initial si vide
    cursor.execute("SELECT COUNT(*) FROM brebis")
    if cursor.fetchone()[0] == 0:
        races = list(STANDARDS_RACES.keys())
        for i in range(15):
            race = random.choice(races)
            sexe = 'F' if random.random() > 0.3 else 'M'
            p_min, p_max = STANDARDS_RACES[race]['poids_adulte']['femelle' if sexe == 'F' else 'male']
            cursor.execute('''
                INSERT INTO brebis (identifiant, nom, race, sexe, poids, date_naissance, score_condition)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            ''', (f"DZ-{race[:2]}-{i:03d}", f"Animal {i}", race, sexe, random.uniform(p_min, p_max), 
                  (date.today() - timedelta(days=random.randint(300, 1000))).isoformat(), random.randint(2, 4)))
    
    conn.commit()
    return conn

db_conn = init_db()

# ========== LOGIQUE SCANNER 3D (AVEC STANDARD 1 MÈTRE) ==========

class ScannerOvin:
    @staticmethod
    def simuler_scan(distance_camera=1.0):
        """
        Simule un scan avec un étalonnage basé sur la distance (Standard 1m).
        """
        # Le facteur de correction simule la précision optique à 1 mètre
        precision_factor = 1.0 if distance_camera == 1.0 else (1.0 / distance_camera)
        
        n_points = 1200
        # Génération d'un nuage de points en forme d'ellipsoïde pour l'ovin
        t = np.linspace(0, 2*np.pi, n_points)
        u = np.linspace(0, np.pi, n_points)
        
        # Dimensions moyennes d'une brebis (en cm) corrigées par la distance
        x = 40 * np.outer(np.cos(t), np.sin(u)) * precision_factor
        y = 100 * np.outer(np.sin(t), np.sin(u)) * precision_factor
        z = 60 * np.outer(np.ones(np.size(t)), np.cos(u)) * precision_factor
        
        return x.flatten(), y.flatten(), z.flatten()

# ========== INTERFACE UTILISATEUR STREAMLIT ==========

def main():
    st.markdown('<h1 class="main-header">OVIN MANAGER PRO v2.0</h1>', unsafe_allow_html=True)
    
    tabs = st.tabs(["📊 Dashboard", "🔍 Scanner 3D AI", "🧬 Génétique & QTL", "📋 Inventaire", "🔊 Echo-Assist"])

    # --- TAB 1: DASHBOARD ---
    with tabs[0]:
        df = pd.read_sql_query("SELECT * FROM brebis", db_conn)
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown(f'<div class="metric-card"><h3>Effectif Total</h3><h2>{len(df)}</h2><p>Têtes</p></div>', unsafe_allow_html=True)
        with col2:
            st.markdown(f'<div class="metric-card"><h3>Poids Moyen</h3><h2>{df["poids"].mean():.1f}</h2><p>kg</p></div>', unsafe_allow_html=True)
        with col3:
            st.markdown(f'<div class="metric-card"><h3>Race Dominante</h3><h2>{df["race"].mode()[0]}</h2><p>Majorité</p></div>', unsafe_allow_html=True)
        with col4:
            st.markdown(f'<div class="metric-card"><h3>Santé Troupeau</h3><h2>94%</h2><p>Indice Global</p></div>', unsafe_allow_html=True)
            
        st.markdown("---")
        c1, c2 = st.columns(2)
        with c1:
            fig_race = px.pie(df, names='race', title="Répartition par Race", hole=0.4, color_discrete_sequence=px.colors.sequential.RdBu)
            st.plotly_chart(fig_race, use_container_width=True)
        with c2:
            fig_weight = px.histogram(df, x='poids', color='race', title="Distribution des Poids par Race", marginal="rug")
            st.plotly_chart(fig_weight, use_container_width=True)

    # --- TAB 2: SCANNER 3D AI (STANDARD 1M) ---
    with tabs[1]:
        st.subheader("🛰️ Analyse Morphologique 3D par LiDAR")
        
        col_s1, col_s2 = st.columns([1, 2])
        
        with col_s1:
            st.info("Configuration du capteur")
            dist = st.slider("Distance de scan (mètres)", 0.5, 3.0, 1.0, step=0.1)
            if dist == 1.0:
                st.success("✅ Standard 1 mètre respecté - Précision Optimale")
            else:
                st.warning("⚠️ Étalonnage requis hors standard 1m")
                
            mode = st.selectbox("Mode de Scan", ["Morphologie complète", "Volume carcasse", "Estimation état gras"])
            
            if st.button("LANCER LE SCAN AI", use_container_width=True):
                with st.spinner("Reconstruction du nuage de points..."):
                    st.session_state['scan_data'] = ScannerOvin.simuler_scan(dist)
                    st.toast("Scan terminé avec succès !")

        with col_s2:
            if 'scan_data' in st.session_state:
                x, y, z = st.session_state['scan_data']
                fig_3d = go.Figure(data=[go.Scatter3d(
                    x=x[::10], y=y[::10], z=z[::10],
                    mode='markers',
                    marker=dict(size=2, color=z[::10], colorscale='Viridis', opacity=0.8)
                )])
                fig_3d.update_layout(scene=dict(aspectmode='data'), margin=dict(l=0, r=0, b=0, t=0))
                st.plotly_chart(fig_3d, use_container_width=True)
                
                # Résultats de l'analyse AI
                st.markdown("""
                <div class="scanner-container">
                    > ANALYSE AI TERMINEE...<br>
                    > LONGUEUR ESTIME : 104.2 cm<br>
                    > HAUTEUR GARROT : 72.1 cm<br>
                    > VOLUME CORPOREL : 0.082 m3<br>
                    > POIDS ESTIME VIA SCAN : 64.5 kg (+/- 0.8)
                </div>
                """, unsafe_allow_html=True)
            else:
                st.image("https://img.freepik.com/vecteurs-premium/concept-technologie-numerique-radar-ligne-balayage-bleu_356415-18.jpg", caption="Prêt pour le scan")

    # --- TAB 3: GÉNÉTIQUE & QTL ---
    with tabs[2]:
        st.subheader("🧬 Laboratoire de Sélection Génomique")
        
        qtl_traits = ["Production Laitière", "Vitesse Croissance", "Résistance Parasitaire", "Prolificité"]
        selected_qtl = st.selectbox("Sélectionner un QTL à analyser", qtl_traits)
        
        col_g1, col_g2 = st.columns(2)
        
        with col_g1:
            # Simulation d'un graphique de Manhattan
            positions = np.arange(100)
            p_values = -np.log10(np.random.uniform(0, 1, 100))
            p_values[random.randint(0,99)] = 7.5 # Signal fort
            
            fig_qtl = px.scatter(x=positions, y=p_values, title=f"Analyse GWAS : {selected_qtl}",
                                labels={'x': 'Chromosome / Position', 'y': '-log10(p-value)'})
            fig_qtl.add_hline(y=5, line_dash="dash", line_color="red", annotation_text="Seuil Significativité")
            st.plotly_chart(fig_qtl, use_container_width=True)
            
        with col_g2:
            st.markdown("### Top Géniteurs (Index Génomique)")
            geno_df = df[df['race'] == 'OUDA'].head(5)
            st.table(geno_df[['identifiant', 'race', 'poids', 'score_condition']])
            st.button("Générer Certificat de Lignée")

    # --- TAB 4: INVENTAIRE ---
    with tabs[3]:
        st.subheader("📋 Gestion de la Base de Données")
        edited_df = st.data_editor(df, num_rows="dynamic")
        if st.button("Sauvegarder les modifications"):
            edited_df.to_sql('brebis', db_conn, if_exists='replace', index=False)
            st.success("Base de données mise à jour !")

    # --- TAB 5: ECHO-ASSIST (ASSISTANT VOCAL SIMULÉ) ---
    with tabs[4]:
        st.subheader("🔊 Echo-Assist : Interface Vocale")
        st.info("Posez une question sur votre troupeau (Simulation)")
        user_voice = st.text_input("Commande vocale (ex: 'Quel est l'état de la Hamra 05 ?')")
        
        if user_voice:
            st.write("🤖 **Echo répond :**")
            st.write(f"D'après les dernières analyses du {date.today()}, l'animal présente un score de condition de 3.5. Sa croissance est supérieure de 12% par rapport à la moyenne de la race {df['race'].iloc[0]}.")

# Lancement de l'application
if __name__ == "__main__":
    main()
