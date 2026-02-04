"""
EXPERT OVIN DZ PRO - Version Cloud Ultime
Fonctionne sans OpenCV ni dépendances système
"""

# Imports de base garantis
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime, timedelta
import io
import base64
import json
from dataclasses import dataclass
from typing import Dict, List, Optional

# Désactivé sur cloud
OPENCV_AVAILABLE = False
RPY2_AVAILABLE = False
FIREBASE_AVAILABLE = False

# Tentative imports optionnels (sans blocage)
try:
    from PIL import Image, ImageDraw, ImageFont
    PIL_AVAILABLE = True
except:
    PIL_AVAILABLE = False

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate, Table, Paragraph
    from reportlab.lib.styles import getSampleStyleSheet
    REPORTLAB_AVAILABLE = True
except:
    REPORTLAB_AVAILABLE = False

try:
    import openpyxl
    OPENPYXL_AVAILABLE = True
except:
    OPENPYXL_AVAILABLE = False

try:
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    SKLEARN_AVAILABLE = True
except:
    SKLEARN_AVAILABLE = False

# Configuration Streamlit
st.set_page_config(
    page_title="Expert Ovin DZ Pro",
    page_icon="🐑",
    layout="wide"
)

# =============================================================================
# APPLICATION SIMPLIFIÉE MAIS COMPLÈTE
# =============================================================================

def main():
    st.title("🐑 Expert Ovin DZ Pro")
    st.caption("Version Cloud - Mode sans OpenCV activé")
    
    # Sidebar navigation simplifiée
    menu = st.sidebar.radio(
        "Navigation",
        ["🏠 Accueil", "📏 Saisie Manuelle", "⚖️ Comparatif", "📊 Excel", "📈 Stats"]
    )
    
    if menu == "🏠 Accueil":
        render_home()
    elif menu == "📏 Saisie Manuelle":
        render_saisie()
    elif menu == "⚖️ Comparatif":
        render_comparatif()
    elif menu == "📊 Excel":
        render_excel()
    elif menu == "📈 Stats":
        render_stats()

def render_home():
    """Page d'accueil avec démo des fonctionnalités"""
    
    st.markdown("""
    ## Bienvenue dans Expert Ovin DZ Pro
    
    Cette version cloud fonctionne sans :
    - ❌ OpenCV (traitement image désactivé)
    - ❌ R (export R désactivé, CSV disponible)
    - ❌ Firebase (stockage local SQLite uniquement)
    
    **Fonctionnalités disponibles :**
    - ✅ Saisie manuelle complète avec calculs
    - ✅ Mode comparatif 2 animaux
    - ✅ Export Excel professionnel
    - ✅ Statistiques de base + ML (clustering)
    """)
    
    # Démonstration visuelle avec données simulées
    st.subheader("📊 Démonstration - Troupeau exemple")
    
    np.random.seed(42)
    demo_data = pd.DataFrame({
        'Animal_ID': [f'MOUTON_{i:03d}' for i in range(1, 21)],
        'Race': np.random.choice(['Lacaune', 'Manech', 'Basco'], 20),
        'Hauteur_Garrot': np.random.normal(68, 3, 20).round(1),
        'Longueur_Corps': np.random.normal(78, 4, 20).round(1),
        'Score_IAL': np.random.normal(75, 8, 20).round(1)
    })
    
    col1, col2, col3 = st.columns(3)
    col1.metric("Total animaux", len(demo_data))
    col2.metric("IAL moyen", f"{demo_data['Score_IAL'].mean():.1f}")
    col3.metric("Race majoritaire", demo_data['Race'].mode()[0])
    
    # Graphique démo
    fig = px.scatter(demo_data, x='Hauteur_Garrot', y='Longueur_Corps', 
                    color='Score_IAL', size='Score_IAL',
                    hover_data=['Animal_ID', 'Race'],
                    title="Nuage de points morphométrique (données démo)")
    st.plotly_chart(fig, use_container_width=True)
    
    st.dataframe(demo_data, use_container_width=True)

def render_saisie():
    """Formulaire saisie manuelle complet"""
    
    st.header("📏 Saisie Manuelle au Ruban")
    
    with st.form("saisie_complete"):
        
        # Section 1: ID
        st.subheader("🐑 Identification")
        col1, col2, col3 = st.columns(3)
        with col1:
            animal_id = st.text_input("ID Animal", f"MOUTON_{datetime.now().strftime('%H%M%S')}")
            race = st.selectbox("Race", ["Lacaune", "Manech", "Basco-Béarnaise", "Autre"])
        with col2:
            date_mesure = st.date_input("Date", datetime.now())
            age_mois = st.number_input("Âge (mois)", 8, 180, 24)
        with col3:
            operateur = st.text_input("Opérateur", "Tech")
            lactation = st.number_input("N° Lactation", 0, 10, 1)
        
        # Section 2: Mensurations
        st.subheader("📐 Mensurations Corporelles (cm)")
        col1, col2, col3 = st.columns(3)
        with col1:
            h_garrot = st.number_input("Hauteur garrot", 50.0, 90.0, 68.0, 0.5)
            l_corps = st.number_input("Longueur corps", 60.0, 100.0, 78.0, 0.5)
        with col2:
            l_bassin = st.number_input("Largeur bassin", 15.0, 30.0, 21.0, 0.5)
            t_poitrine = st.number_input("Tour poitrine", 70.0, 120.0, 92.0, 0.5)
        with col3:
            p_poitrine = st.number_input("Profondeur poitrine", 25.0, 45.0, 34.0, 0.5)
            l_epaules = st.number_input("Largeur épaules", 15.0, 30.0, 22.0, 0.5)
        
        # Section 3: Mamelle
        st.subheader("🥛 Mensurations Mamelle (cm)")
        col1, col2, col3 = st.columns(3)
        with col1:
            l_mamelle = st.number_input("Longueur mamelle", 10.0, 30.0, 18.0, 0.5)
            lr_mamelle = st.number_input("Largeur mamelle", 8.0, 25.0, 15.0, 0.5)
        with col2:
            pr_mamelle = st.number_input("Profondeur mamelle", 8.0, 25.0, 16.0, 0.5)
            ecart_t = st.number_input("Écart tétines", 5.0, 20.0, 11.0, 0.5)
        with col3:
            cirf_g = st.number_input("Circonf. tétine G", 2.0, 5.0, 3.2, 0.1)
            cirf_d = st.number_input("Circonf. tétine D", 2.0, 5.0, 3.2, 0.1)
        
        # Section 4: Scores
        st.subheader("⭐ Scores Subjectifs (1-9)")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            sc_attache = st.slider("Attache", 1, 9, 7)
        with col2:
            sc_prof = st.slider("Profondeur", 1, 9, 7)
        with col3:
            sc_sym = st.slider("Symétrie", 1, 9, 8)
        with col4:
            sc_ec = st.slider("État corporel", 1, 9, 5)
        
        # Bouton validation
        submitted = st.form_submit_button("💾 Calculer et Enregistrer", type="primary")
        
        if submitted:
            # Calculs
            indice_conf = (l_corps / h_garrot) * 100
            poids_est = (t_poitrine ** 2) * l_corps / 10800
            volume_mam = (4/3) * 3.14159 * (l_mamelle/2) * (lr_mamelle/2) * (pr_mamelle/2)
            score_global = (sc_attache * 0.3 + sc_prof * 0.25 + sc_sym * 0.25 + (10-abs(sc_ec-5))*0.2)
            
            # Stockage session
            if 'mesures' not in st.session_state:
                st.session_state.mesures = []
            
            st.session_state.mesures.append({
                'id': animal_id,
                'race': race,
                'date': date_mesure,
                'h_garrot': h_garrot,
                'l_corps': l_corps,
                'indice_conf': round(indice_conf, 2),
                'poids_kg': round(poids_est, 2),
                'volume_mam': round(volume_mam, 2),
                'score_global': round(score_global, 2)
            })
            
            # Résultats
            st.success(f"✅ {animal_id} enregistré!")
            
            col_r1, col_r2, col_r3, col_r4 = st.columns(4)
            col_r1.metric("Indice Conf.", f"{indice_conf:.1f}")
            col_r2.metric("Poids estimé", f"{poids_est:.1f} kg")
            col_r3.metric("Volume mamelle", f"{volume_mam:.0f} cm³")
            col_r4.metric("Score global", f"{score_global:.1f}/9")
    
    # Affichage historique
    if st.session_state.get('mesures'):
        st.subheader("📋 Historique session")
        st.dataframe(pd.DataFrame(st.session_state.mesures), use_container_width=True)

def render_comparatif():
    """Mode comparatif 2 animaux"""
    
    st.header("⚖️ Mode Comparatif")
    
    mesures = st.session_state.get('mesures', [])
    
    if len(mesures) < 2:
        st.warning("Il faut au moins 2 animaux saisis. Allez dans 'Saisie Manuelle' d'abord.")
        
        # Données démo pour test
        if st.checkbox("Utiliser données démo"):
            mesures = [
                {'id': 'DEMO_A', 'race': 'Lacaune', 'h_garrot': 70.0, 'l_corps': 80.0, 
                 'volume_mam': 2500, 'score_global': 85.0, 'indice_conf': 114.3},
                {'id': 'DEMO_B', 'race': 'Lacaune', 'h_garrot': 68.0, 'l_corps': 76.0,
                 'volume_mam': 2100, 'score_global': 78.0, 'indice_conf': 111.8},
            ]
        else:
            return
    
    # Sélection
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Animal A (Référence)")
        choix_a = st.selectbox("Choisir A", [m['id'] for m in mesures], key='sel_a')
        a = next(m for m in mesures if m['id'] == choix_a)
    
    with col2:
        st.subheader("Animal B (À évaluer)")
        choix_b = st.selectbox("Choisir B", [m['id'] for m in mesures if m['id'] != choix_a], key='sel_b')
        b = next(m for m in mesures if m['id'] == choix_b)
    
    # Comparaison visuelle
    st.markdown("---")
    
    # Scores côte à côte
    col_a, col_mid, col_b = st.columns([2, 1, 2])
    
    with col_a:
        st.markdown(f"""
        <div style='background-color: #E3F2FD; padding: 20px; border-radius: 10px; text-align: center;'>
            <h2>{a['id']}</h2>
            <h1 style='color: #1976D2;'>{a['score_global']:.1f}</h1>
            <p>Score Global</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col_mid:
        diff = b['score_global'] - a['score_global']
        emoji = "🟢" if diff > 0 else "🔴" if diff < 0 else "⚪"
        st.markdown(f"""
        <div style='text-align: center; padding-top: 40px;'>
            <h1>{emoji}</h1>
            <h2>{diff:+.1f}</h2>
            <p>points</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col_b:
        color_b = '#4CAF50' if diff > 0 else '#F44336'
        st.markdown(f"""
        <div style='background-color: {'#E8F5E9' if diff > 0 else '#FFEBEE'}; padding: 20px; border-radius: 10px; text-align: center;'>
            <h2>{b['id']}</h2>
            <h1 style='color: {color_b};'>{b['score_global']:.1f}</h1>
            <p>{'Meilleur' if diff > 0 else 'Inférieur'}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Graphique radar
    st.subheader("🕸️ Profil Comparatif")
    
    categories = ['Hauteur', 'Longueur', 'Volume\nMamelle', 'Score\nGlobal', 'Indice\nConf.']
    
    # Normalisation 0-10
    def norm(val, min_v, max_v):
        return min(10, max(0, (val - min_v) / (max_v - min_v) * 10))
    
    values_a = [
        norm(a['h_garrot'], 60, 75),
        norm(a['l_corps'], 70, 85),
        norm(a['volume_mam'], 1500, 3000),
        a['score_global'],
        norm(a['indice_conf'], 100, 120)
    ]
    
    values_b = [
        norm(b['h_garrot'], 60, 75),
        norm(b['l_corps'], 70, 85),
        norm(b['volume_mam'], 1500, 3000),
        b['score_global'],
        norm(b['indice_conf'], 100, 120)
    ]
    
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=values_a + [values_a[0]],
        theta=categories + [categories[0]],
        fill='toself',
        name=f"{a['id']} (Réf)",
        line_color='#2196F3'
    ))
    fig.add_trace(go.Scatterpolar(
        r=values_b + [values_b[0]],
        theta=categories + [categories[0]],
        fill='toself',
        name=f"{b['id']} (Test)",
        line_color='#4CAF50' if diff > 0 else '#F44336'
    ))
    fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 10])))
    st.plotly_chart(fig, use_container_width=True)
    
    # Recommandation
    st.subheader("🎯 Recommandation")
    if diff > 5:
        st.success(f"✅ **Sélectionner {b['id']}** - Supériorité significative")
    elif diff > 0:
        st.info(f"ℹ️ **Légère préférence {b['id']}** - Avantage marginal")
    elif diff < -5:
        st.error(f"❌ **Écarter {b['id']}** - Infériorité significative")
    else:
        st.warning(f"⚖️ **Équivalence** - Critères secondaires à considérer")

def render_excel():
    """Import/Export Excel"""
    
    st.header("📊 Import/Export Excel")
    
    if not OPENPYXL_AVAILABLE:
        st.error("Module Excel non disponible")
        return
    
    tabs = st.tabs(["📥 Import", "📤 Export", "📄 Template"])
    
    with tabs[0]:
        st.subheader("Importer fichier Excel")
        uploaded = st.file_uploader("Fichier .xlsx", type=['xlsx', 'xls'])
        
        if uploaded:
            try:
                df = pd.read_excel(uploaded)
                st.success(f"{len(df)} lignes importées")
                st.dataframe(df.head())
                
                if st.button("Ajouter aux données"):
                    st.session_state['imported_data'] = df
                    st.balloons()
            except Exception as e:
                st.error(f"Erreur: {e}")
    
    with tabs[1]:
        st.subheader("Exporter vers Excel")
        
        mesures = st.session_state.get('mesures', [])
        if not mesures:
            st.warning("Aucune donnée à exporter")
            return
        
        df_export = pd.DataFrame(mesures)
        
        # Création Excel en mémoire
        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            df_export.to_excel(writer, sheet_name='Données', index=False)
            
            # Feuille stats
            stats = pd.DataFrame({
                'Statistique': ['Moyenne', 'Min', 'Max', 'Écart-type'],
                'Score Global': [
                    df_export['score_global'].mean(),
                    df_export['score_global'].min(),
                    df_export['score_global'].max(),
                    df_export['score_global'].std()
                ]
            })
            stats.to_excel(writer, sheet_name='Stats', index=False)
        
        st.download_button(
            "⬇️ Télécharger Excel",
            buffer.getvalue(),
            f"Export_Ovin_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
        
        st.dataframe(df_export, use_container_width=True)
    
    with tabs[2]:
        st.subheader("Générer Template vierge")
        
        if st.button("📄 Créer template"):
            buffer = io.BytesIO()
            
            # Template avec openpyxl
            from openpyxl import Workbook
            from openpyxl.styles import Font, PatternFill
            
            wb = Workbook()
            ws = wb.active
            ws.title = "Saisie"
            
            headers = ['ID_Animal', 'Date', 'Race', 'Age_mois', 'Hauteur_Garrot', 
                      'Longueur_Corps', 'Score_Attache', 'Score_Symetrie']
            
            for col, h in enumerate(headers, 1):
                cell = ws.cell(1, col, h)
                cell.font = Font(bold=True, color="FFFFFF")
                cell.fill = PatternFill(start_color="366092", end_color="366092", fill_type="solid")
            
            wb.save(buffer)
            
            st.download_button(
                "⬇️ Télécharger Template",
                buffer.getvalue(),
                "Template_Saisie_Ovin.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

def render_stats():
    """Statistiques et ML basique"""
    
    st.header("📈 Statistiques & Analyse")
    
    mesures = st.session_state.get('mesures', [])
    
    if len(mesures) < 3:
        st.warning(f"Il faut au moins 3 animaux pour les stats (actuel: {len(mesures)})")
        
        # Générer données démo
        if st.button("Générer données démo (20 animaux)"):
            np.random.seed(42)
            for i in range(20):
                st.session_state.mesures.append({
                    'id': f'DEMO_{i+1:03d}',
                    'race': np.random.choice(['Lacaune', 'Manech']),
                    'h_garrot': np.random.normal(68, 3),
                    'l_corps': np.random.normal(78, 4),
                    'volume_mam': np.random.normal(2200, 400),
                    'score_global': np.random.normal(76, 6),
                    'indice_conf': np.random.normal(113, 4)
                })
            st.rerun()
        return
    
    df = pd.DataFrame(mesures)
    
    # Stats descriptives
    st.subheader("Statistiques Descriptives")
    st.dataframe(df.describe(), use_container_width=True)
    
    # Graphiques
    col1, col2 = st.columns(2)
    
    with col1:
        fig = px.histogram(df, x='score_global', nbins=10, 
                          title="Distribution Scores",
                          color_discrete_sequence=['#2E5090'])
        st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        fig = px.scatter(df, x='volume_mam', y='score_global', 
                        color='race', size='h_garrot',
                        title="Volume Mamelle vs Score")
        st.plotly_chart(fig, use_container_width=True)
    
    # Clustering si sklearn disponible
    if SKLEARN_AVAILABLE and len(mesures) >= 5:
        st.subheader("🧬 Clustering Automatique (K-means)")
        
        n_clusters = st.slider("Nombre groupes", 2, 5, 3)
        
        # Préparation données
        X = df[['h_garrot', 'l_corps', 'volume_mam', 'score_global']].values
        X_scaled = StandardScaler().fit_transform(X)
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        df['cluster'] = kmeans.fit_predict(X_scaled)
        
        fig = px.scatter(df, x='volume_mam', y='score_global', 
                        color='cluster.astype(str)',
                        hover_data=['id'],
                        title=f"Classification en {n_clusters} groupes")
        st.plotly_chart(fig, use_container_width=True)
        
        # Profils clusters
        st.write("Profils des groupes:")
        st.dataframe(df.groupby('cluster')[['h_garrot', 'score_global']].mean())

if __name__ == "__main__":
    main()
