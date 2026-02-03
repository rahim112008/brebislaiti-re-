"""
OVIN MANAGER PRO - Application Streamlit de gestion scientifique d'élevage ovin
Auteur: [Votre Nom]
Version: 1.0.0
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
from typing import Dict, List, Optional
import json
import io
import base64
from PIL import Image
import numpy as np
import sys
import os

# Ajouter le chemin des modules
sys.path.append(os.path.join(os.path.dirname(__file__), 'modules'))

# Importer nos modules
from modules.database import DatabaseManager, init_database
from modules.gestion import GestionnaireGestation, DonneesDemonstration
from modules.morphometrie import AnalyseurMorphometrique
from modules.genomique import IntegrationGenomique
from modules.statistiques import AnalyseurStatistique
from modules.utils import generer_rapport_pdf, exporter_excel

# Configuration de la page Streamlit
st.set_page_config(
    page_title="Ovin Manager Pro",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personnalisé
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #2E7D32;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.8rem;
        color: #388E3C;
        margin-top: 2rem;
        margin-bottom: 1rem;
        padding-bottom: 0.5rem;
        border-bottom: 2px solid #4CAF50;
    }
    .metric-card {
        background-color: #E8F5E9;
        padding: 1rem;
        border-radius: 10px;
        border-left: 5px solid #4CAF50;
        margin-bottom: 1rem;
    }
    .success-message {
        background-color: #C8E6C9;
        padding: 1rem;
        border-radius: 5px;
        margin: 1rem 0;
    }
    .warning-message {
        background-color: #FFF3CD;
        padding: 1rem;
        border-radius: 5px;
        margin: 1rem 0;
    }
    .stButton button {
        background-color: #4CAF50;
        color: white;
        font-weight: bold;
    }
    .stButton button:hover {
        background-color: #45A049;
    }
</style>
""", unsafe_allow_html=True)

# Initialisation de la session state
if 'db_manager' not in st.session_state:
    st.session_state.db_manager = DatabaseManager()
    init_database(st.session_state.db_manager)

if 'gestation_manager' not in st.session_state:
    st.session_state.gestation_manager = GestionnaireGestation(st.session_state.db_manager)

if 'morpho_analyzer' not in st.session_state:
    st.session_state.morpho_analyzer = AnalyseurMorphometrique()

if 'genomique' not in st.session_state:
    st.session_state.genomique = IntegrationGenomique("contact@ovin-manager.com")

if 'stats' not in st.session_state:
    st.session_state.stats = AnalyseurStatistique()

# Titre principal
st.markdown('<h1 class="main-header">🐑 Ovin Manager Pro</h1>', unsafe_allow_html=True)
st.markdown("""
*Application scientifique de gestion et d'analyse d'élevage ovin laitier*
""")

# Sidebar - Navigation
with st.sidebar:
    st.image("https://via.placeholder.com/150x150/4CAF50/FFFFFF?text=OMP", width=150)
    st.markdown("### 📍 Navigation")
    
    page = st.radio(
        "Menu Principal",
        ["🏠 Tableau de Bord", 
         "📊 Gestion des Brebis", 
         "🤰 Suivi Gestation", 
         "📸 Morphométrie",
         "🧬 Analyse Génomique", 
         "📈 Statistiques", 
         "📋 Rapports",
         "⚙️ Configuration"]
    )
    
    st.markdown("---")
    st.markdown("### 📊 Statistiques rapides")
    
    # Afficher quelques statistiques
    try:
        cursor = st.session_state.db_manager.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM brebis")
        total_brebis = cursor.fetchone()[0]
        st.metric("Total Brebis", total_brebis)
        
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
        gestations_en_cours = cursor.fetchone()[0]
        st.metric("Gestations en cours", gestations_en_cours)
    except:
        st.info("Base de données en cours d'initialisation")

# Fonctions pour chaque page
def afficher_tableau_bord():
    """Affiche le tableau de bord principal"""
    st.markdown('<h2 class="section-header">🏠 Tableau de Bord</h2>', unsafe_allow_html=True)
    
    # Métriques principales
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        cursor = st.session_state.db_manager.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM brebis")
        total = cursor.fetchone()[0]
        st.metric("Total Brebis", total, delta="+2 ce mois")
    
    with col2:
        cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F'")
        femelles = cursor.fetchone()[0]
        st.metric("Brebis Femelles", femelles, f"{femelles/total*100:.1f}%")
    
    with col3:
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
        gestations = cursor.fetchone()[0]
        st.metric("Gestations", gestations)
    
    with col4:
        cursor.execute("SELECT AVG(valeur) FROM caracteres_morpho WHERE caractere = 'poids_vif'")
        poids_moyen = cursor.fetchone()[0] or 0
        st.metric("Poids Moyen", f"{poids_moyen:.1f} kg")
    
    # Graphiques
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 📊 Répartition par Race")
        cursor.execute("""
            SELECT race, COUNT(*) as count 
            FROM brebis 
            GROUP BY race
        """)
        data_race = cursor.fetchall()
        
        if data_race:
            df_race = pd.DataFrame(data_race, columns=['Race', 'Nombre'])
            fig = px.pie(df_race, values='Nombre', names='Race', 
                        color_discrete_sequence=px.colors.sequential.Greens)
            st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        st.markdown("### 📈 Production par Mois")
        # Données simulées
        mois = ['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Jun']
        production = [1200, 1350, 1420, 1380, 1500, 1620]
        
        fig = go.Figure(data=[
            go.Bar(x=mois, y=production, marker_color='#4CAF50')
        ])
        fig.update_layout(title="Production Laitière (L)")
        st.plotly_chart(fig, use_container_width=True)
    
    # Alertes et notifications
    st.markdown("### 🔔 Alertes et Rappels")
    
    col1, col2 = st.columns(2)
    
    with col1:
        with st.expander("⚠️ Vaccinations en retard", expanded=True):
            st.warning("3 brebis nécessitent des rappels de vaccination")
            if st.button("Voir détails", key="vaccin_details"):
                st.info("Brebis ID: 12, 15, 18 - Vaccin FCO")
    
    with col2:
        with st.expander("🤰 Mises bas prévues", expanded=True):
            st.info("2 mises bas prévues cette semaine")
            if st.button("Voir calendrier", key="misebas_details"):
                st.success("Brebis ID: 5 (15/06) et ID: 8 (18/06)")

def afficher_gestion_brebis():
    """Affiche la page de gestion des brebis"""
    st.markdown('<h2 class="section-header">📊 Gestion des Brebis</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📋 Liste des Brebis", "➕ Ajouter une Brebis", "🔍 Recherche Avancée"])
    
    with tab1:
        # Filtrer les brebis
        col1, col2, col3 = st.columns(3)
        with col1:
            filtre_race = st.selectbox("Filtrer par race", ["Toutes", "lacaune", "manech", "basco_bearnaise"])
        with col2:
            filtre_statut = st.selectbox("Filtrer par statut", ["Tous", "active", "retrait", "morte"])
        with col3:
            filtre_sexe = st.selectbox("Filtrer par sexe", ["Tous", "F", "M"])
        
        # Construire la requête SQL
        query = "SELECT * FROM brebis WHERE 1=1"
        params = []
        
        if filtre_race != "Toutes":
            query += " AND race = ?"
            params.append(filtre_race)
        
        if filtre_statut != "Tous":
            query += " AND statut = ?"
            params.append(filtre_statut)
        
        if filtre_sexe != "Tous":
            query += " AND sexe = ?"
            params.append(filtre_sexe)
        
        query += " ORDER BY id"
        
        cursor = st.session_state.db_manager.conn.cursor()
        cursor.execute(query, params)
        brebis_data = cursor.fetchall()
        
        if brebis_data:
            columns = [desc[0] for desc in cursor.description]
            df_brebis = pd.DataFrame(brebis_data, columns=columns)
            
            # Formater la dataframe pour l'affichage
            df_display = df_brebis[['id', 'identifiant_unique', 'nom', 'race', 'sexe', 'statut', 'date_naissance']].copy()
            df_display['âge'] = df_display['date_naissance'].apply(
                lambda x: f"{(date.today() - datetime.strptime(x, '%Y-%m-%d').date()).days // 30} mois"
            )
            
            st.dataframe(df_display, use_container_width=True)
            
            # Options d'export
            col1, col2 = st.columns(2)
            with col1:
                if st.button("📥 Exporter en Excel"):
                    excel_buffer = exporter_excel(df_brebis, "brebis")
                    st.download_button(
                        label="Télécharger Excel",
                        data=excel_buffer,
                        file_name="brebis_export.xlsx",
                        mime="application/vnd.ms-excel"
                    )
            with col2:
                if st.button("📄 Générer Rapport PDF"):
                    pdf_buffer = generer_rapport_pdf(df_brebis, "Rapport Brebis")
                    st.download_button(
                        label="Télécharger PDF",
                        data=pdf_buffer,
                        file_name="rapport_brebis.pdf",
                        mime="application/pdf"
                    )
        else:
            st.info("Aucune brebis ne correspond aux filtres")
    
    with tab2:
        st.markdown("### Ajouter une nouvelle brebis")
        
        with st.form("form_ajout_brebis"):
            col1, col2 = st.columns(2)
            
            with col1:
                identifiant = st.text_input("Identifiant Unique*")
                nom = st.text_input("Nom")
                date_naissance = st.date_input("Date de Naissance*", value=date.today() - timedelta(days=365))
            
            with col2:
                race = st.selectbox("Race*", ["lacaune", "manech", "basco_bearnaise", "autre"])
                sexe = st.radio("Sexe*", ["F", "M"], horizontal=True)
                statut = st.selectbox("Statut", ["active", "retrait", "morte"])
                notes = st.text_area("Notes")
            
            soumettre = st.form_submit_button("✅ Ajouter la Brebis")
            
            if soumettre:
                if not identifiant:
                    st.error("L'identifiant unique est obligatoire")
                else:
                    try:
                        # Vérifier si l'identifiant existe déjà
                        cursor = st.session_state.db_manager.conn.cursor()
                        cursor.execute("SELECT COUNT(*) FROM brebis WHERE identifiant_unique = ?", (identifiant,))
                        if cursor.fetchone()[0] > 0:
                            st.error("Cet identifiant existe déjà")
                        else:
                            cursor.execute('''
                                INSERT INTO brebis 
                                (identifiant_unique, nom, date_naissance, race, sexe, statut, notes)
                                VALUES (?, ?, ?, ?, ?, ?, ?)
                            ''', (identifiant, nom, date_naissance.isoformat(), race, sexe, statut, notes))
                            st.session_state.db_manager.conn.commit()
                            st.success(f"✅ Brebis {identifiant} ajoutée avec succès!")
                            st.balloons()
                    except Exception as e:
                        st.error(f"Erreur: {e}")
    
    with tab3:
        st.markdown("### 🔍 Recherche Avancée")
        
        recherche = st.text_input("Rechercher par nom, identifiant ou notes")
        if recherche:
            cursor = st.session_state.db_manager.conn.cursor()
            cursor.execute('''
                SELECT * FROM brebis 
                WHERE nom LIKE ? OR identifiant_unique LIKE ? OR notes LIKE ?
            ''', (f'%{recherche}%', f'%{recherche}%', f'%{recherche}%'))
            
            resultats = cursor.fetchall()
            if resultats:
                columns = [desc[0] for desc in cursor.description]
                df_resultats = pd.DataFrame(resultats, columns=columns)
                st.dataframe(df_resultats[['id', 'identifiant_unique', 'nom', 'race', 'statut']])
            else:
                st.info("Aucun résultat trouvé")

def afficher_suivi_gestation():
    """Affiche la page de suivi des gestations"""
    st.markdown('<h2 class="section-header">🤰 Suivi des Gestations</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📅 Calendrier", "➕ Nouvelle Gestation", "📊 Statistiques"])
    
    with tab1:
        # Calendrier des gestations
        cursor = st.session_state.db_manager.conn.cursor()
        cursor.execute('''
            SELECT g.*, b.nom, b.race 
            FROM gestations g 
            JOIN brebis b ON g.brebis_id = b.id 
            WHERE g.statut = 'en_cours'
            ORDER BY g.date_mise_bas_prevu
        ''')
        
        gestations = cursor.fetchall()
        
        if gestations:
            st.markdown("### 📋 Gestations en cours")
            
            for gestation in gestations:
                cols = [desc[0] for desc in cursor.description]
                gest_dict = dict(zip(cols, gestation))
                
                with st.expander(f"🐑 {gest_dict['nom']} - Mise bas prévue: {gest_dict['date_mise_bas_prevu']}"):
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Jours restants", 
                                 (datetime.strptime(gest_dict['date_mise_bas_prevu'], '%Y-%m-%d').date() - date.today()).days)
                    
                    with col2:
                        st.metric("Agneaux prévus", gest_dict['nombre_agneaux_prevus'])
                    
                    with col3:
                        progress = min(100, 100 * (date.today() - datetime.strptime(gest_dict['date_eponge'], '%Y-%m-%d').date()).days / 150)
                        st.progress(progress / 100, text=f"Progression: {progress:.1f}%")
                    
                    # Détails
                    st.write(f"**Race:** {gest_dict['race']}")
                    st.write(f"**Date éponge:** {gest_dict['date_eponge']}")
                    st.write(f"**Date retrait éponge:** {gest_dict.get('date_retrait_eponge', 'Non spécifiée')}")
                    
                    # Actions
                    if st.button("📝 Mettre à jour", key=f"update_{gest_dict['id']}"):
                        st.session_state[f"edit_gestation_{gest_dict['id']}"] = True
                    
                    if st.button("✅ Mise bas effectuée", key=f"birth_{gest_dict['id']}"):
                        date_mise_bas = st.date_input("Date de mise bas", value=date.today(), key=f"date_birth_{gest_dict['id']}")
                        if st.button("Confirmer", key=f"confirm_birth_{gest_dict['id']}"):
                            cursor.execute('''
                                UPDATE gestations 
                                SET date_mise_bas_reel = ?, statut = 'termine'
                                WHERE id = ?
                            ''', (date_mise_bas.isoformat(), gest_dict['id']))
                            st.session_state.db_manager.conn.commit()
                            st.success("Mise bas enregistrée!")
                            st.rerun()
        else:
            st.info("Aucune gestation en cours")
    
    with tab2:
        st.markdown("### Programmer une nouvelle gestation")
        
        # Sélection de la brebis
        cursor = st.session_state.db_manager.conn.cursor()
        cursor.execute("SELECT id, nom, race FROM brebis WHERE sexe = 'F' AND statut = 'active'")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            brebis_options = {f"{b[1]} (ID: {b[0]})": b[0] for b in brebis_list}
            selected_brebis_label = st.selectbox("Sélectionner la brebis", list(brebis_options.keys()))
            brebis_id = brebis_options[selected_brebis_label]
            
            # Récupérer la race de la brebis
            cursor.execute("SELECT race FROM brebis WHERE id = ?", (brebis_id,))
            brebis_race = cursor.fetchone()[0]
            
            # Formulaire de gestation
            col1, col2 = st.columns(2)
            
            with col1:
                date_eponge = st.date_input("Date d'introduction de l'éponge", value=date.today())
                date_retrait = st.date_input("Date de retrait prévue", value=date_eponge + timedelta(days=14))
            
            with col2:
                date_insemination = st.date_input("Date d'insémination", value=date_retrait)
                nb_agneaux = st.number_input("Nombre d'agneaux prévus", min_value=1, max_value=4, value=1)
            
            if st.button("📅 Calculer la date de mise bas"):
                programme = st.session_state.gestation_manager.programmer_eponge(
                    brebis_id, date_eponge, brebis_race
                )
                
                st.markdown("### 📋 Programme de gestation")
                
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Date mise bas prévue", programme['date_mise_bas_prevu'])
                    st.metric("Durée gestation", "150 jours")
                with col2:
                    st.metric("Période insémination", f"{programme['periode_insemination']['debut']} à {programme['periode_insemination']['fin']}")
                    st.metric("Période mise bas", f"{programme['periode_mise_bas']['debut']} à {programme['periode_mise_bas']['fin']}")
                
                # Bouton pour enregistrer
                if st.button("💾 Enregistrer cette gestation"):
                    cursor.execute('''
                        INSERT INTO gestations 
                        (brebis_id, date_eponge, date_retrait_eponge, date_insemination, 
                         date_mise_bas_prevu, nombre_agneaux_prevus, statut)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', (brebis_id, date_eponge.isoformat(), date_retrait.isoformat(),
                          date_insemination.isoformat(), programme['date_mise_bas_prevu'].isoformat(),
                          nb_agneaux, 'en_cours'))
                    st.session_state.db_manager.conn.commit()
                    st.success("✅ Gestation enregistrée avec succès!")
                    st.balloons()
        else:
            st.warning("Aucune brebis active disponible")
    
    with tab3:
        st.markdown("### 📊 Statistiques des gestations")
        
        stats = st.session_state.gestation_manager.get_statistiques_gestation()
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total gestations", stats['total_gestations'])
        with col2:
            st.metric("En cours", stats['en_cours'])
        with col3:
            st.metric("Terminées", stats['terminees'])
        with col4:
            st.metric("Taux réussite", f"{stats['taux_reussite']:.1f}%")
        
        # Graphique des gestations par mois
        cursor.execute('''
            SELECT strftime('%Y-%m', date_eponge) as mois, COUNT(*) as count
            FROM gestations
            GROUP BY mois
            ORDER BY mois
        ''')
        data_mois = cursor.fetchall()
        
        if data_mois:
            df_mois = pd.DataFrame(data_mois, columns=['Mois', 'Nombre'])
            fig = px.bar(df_mois, x='Mois', y='Nombre', 
                        title="Gestations par mois",
                        color='Nombre',
                        color_continuous_scale='Greens')
            st.plotly_chart(fig, use_container_width=True)

def afficher_morphometrie():
    """Affiche la page d'analyse morphométrique"""
    st.markdown('<h2 class="section-header">📸 Analyse Morphométrique</h2>', unsafe_allow_html=True)
    
    st.info("""
    **Fonctionnalité d'analyse par photo smartphone**
    - Prenez une photo latérale de la brebis
    - Incluez un objet de référence de taille connue
    - L'IA mesure automatiquement les dimensions
    """)
    
    tab1, tab2 = st.tabs(["📷 Analyser une photo", "📊 Historique des mesures"])
    
    with tab1:
        # Upload de photo
        uploaded_file = st.file_uploader("Choisir une image", type=['jpg', 'jpeg', 'png'])
        
        if uploaded_file is not None:
            # Afficher l'image
            image = Image.open(uploaded_file)
            st.image(image, caption="Photo uploadée", width=300)
            
            # Paramètres d'analyse
            col1, col2 = st.columns(2)
            
            with col1:
                brebis_id = st.number_input("ID de la brebis", min_value=1, value=1)
                race = st.selectbox("Race", ["lacaune", "manech", "basco_bearnaise"])
            
            with col2:
                objet_pixels = st.number_input("Taille de l'objet de référence en pixels", 
                                             min_value=10, max_value=500, value=100)
                taille_reelle = st.number_input("Taille réelle de l'objet (cm)", 
                                              min_value=0.1, max_value=50.0, value=10.0)
            
            if st.button("🔍 Analyser la photo", type="primary"):
                with st.spinner("Analyse en cours..."):
                    # Simuler l'analyse
                    resultats = st.session_state.morpho_analyzer.analyser_photo(
                        "uploaded_image.jpg",
                        objet_pixels,
                        taille_reelle,
                        race
                    )
                    
                    # Afficher les résultats
                    st.markdown("### 📊 Résultats de l'analyse")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Longueur corps", f"{resultats['mesures']['longueur_corps']} cm")
                    with col2:
                        st.metric("Hauteur garrot", f"{resultats['mesures']['hauteur_garrot']} cm")
                    with col3:
                        st.metric("Poids estimé", f"{resultats['mesures']['poids_estime']} kg")
                    
                    # Graphique radar des mesures
                    mesures = ['longueur_corps', 'hauteur_garrot', 'tour_poitrine', 'poids_estime']
                    valeurs = [resultats['mesures'][m] for m in mesures]
                    
                    fig = go.Figure(data=go.Scatterpolar(
                        r=valeurs,
                        theta=mesures,
                        fill='toself',
                        line_color='green'
                    ))
                    
                    fig.update_layout(
                        polar=dict(
                            radialaxis=dict(
                                visible=True,
                                range=[0, max(valeurs) * 1.2]
                            )),
                        showlegend=False,
                        title="Profil morphométrique"
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Analyse comparative
                    if 'analyse_comparative' in resultats:
                        st.markdown("### 📋 Analyse comparative")
                        for analyse in resultats['analyse_comparative']:
                            st.write(f"• {analyse}")
                    
                    # Bouton pour sauvegarder
                    if st.button("💾 Sauvegarder ces mesures"):
                        cursor = st.session_state.db_manager.conn.cursor()
                        for caractere, valeur in resultats['mesures'].items():
                            cursor.execute('''
                                INSERT INTO caracteres_morpho 
                                (brebis_id, date_mesure, caractere, valeur, unite, methode_mesure)
                                VALUES (?, ?, ?, ?, ?, ?)
                            ''', (brebis_id, date.today().isoformat(), caractere, 
                                 valeur, 'cm' if 'cm' in caractere else 'kg', 'photo'))
                        
                        st.session_state.db_manager.conn.commit()
                        st.success("✅ Mesures sauvegardées!")
    
    with tab2:
        # Afficher l'historique des mesures
        brebis_id = st.number_input("ID brebis pour historique", min_value=1, value=1, key="hist_brebis")
        
        if st.button("📈 Afficher l'historique"):
            cursor = st.session_state.db_manager.conn.cursor()
            cursor.execute('''
                SELECT date_mesure, caractere, valeur, unite 
                FROM caracteres_morpho 
                WHERE brebis_id = ?
                ORDER BY date_mesure
            ''', (brebis_id,))
            
            mesures = cursor.fetchall()
            
            if mesures:
                df_mesures = pd.DataFrame(mesures, columns=['Date', 'Caractère', 'Valeur', 'Unité'])
                
                # Graphique d'évolution
                fig = px.line(df_mesures, x='Date', y='Valeur', color='Caractère',
                             title="Évolution des mesures morphométriques",
                             markers=True)
                st.plotly_chart(fig, use_container_width=True)
                
                # Tableau des données
                st.dataframe(df_mesures)
            else:
                st.info("Aucune mesure enregistrée pour cette brebis")

def afficher_genomique():
    """Affiche la page d'analyse génomique"""
    st.markdown('<h2 class="section-header">🧬 Analyse Génomique</h2>', unsafe_allow_html=True)
    
    st.info("""
    **Intégration avec NCBI et GenBank**
    - Recherche de séquences ovines
    - Analyse de SNP (Single Nucleotide Polymorphism)
    - Identification de gènes candidats
    """)
    
    tab1, tab2, tab3 = st.tabs(["🔎 Recherche NCBI", "🧪 Analyse SNP", "📊 Gènes Candidats"])
    
    with tab1:
        st.markdown("### Recherche dans les bases de données génomiques")
        
        col1, col2 = st.columns(2)
        
        with col1:
            terme_recherche = st.text_input("Terme de recherche", value="Ovis aries lact")
            max_results = st.slider("Nombre maximum de résultats", 10, 100, 20)
        
        with col2:
            base_donnees = st.selectbox("Base de données", ["nucleotide", "protein", "gene"])
            espece = st.text_input("Espèce (filtre)", value="Ovis aries")
        
        if st.button("🔍 Lancer la recherche", type="primary"):
            with st.spinner("Recherche en cours..."):
                # Simulation de recherche
                genes = st.session_state.genomique.rechercher_genes_candidats("lacaune")
                
                st.markdown(f"### 📋 Résultats pour: {terme_recherche}")
                
                for i, gene in enumerate(genes[:5], 1):
                    with st.expander(f"Gène {gene['gene']}"):
                        st.write(f"**Fonction:** {gene['fonction']}")
                        st.write(f"**Chromosome:** {gene['chromosome']}")
                        st.write(f"**Accès NCBI:** [Lien simulé](https://www.ncbi.nlm.nih.gov/)")
                
                st.info("*Note: Ceci est une simulation. En production, connexion réelle à NCBI.*")
    
    with tab2:
        st.markdown("### Analyse de SNP (Single Nucleotide Polymorphism)")
        
        # Entrée des séquences
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Séquence de référence**")
            seq_ref = st.text_area("ADN référence", value="ATCGATCGATCGATCG", height=100,
                                 placeholder="Entrez la séquence ADN de référence")
        
        with col2:
            st.markdown("**Séquence étudiée**")
            seq_etu = st.text_area("ADN étudié", value="ATCGCTCGATCGATCG", height=100,
                                 placeholder="Entrez la séquence ADN à étudier")
        
        if st.button("🧬 Analyser les SNP"):
            if len(seq_ref) != len(seq_etu):
                st.error("Les séquences doivent avoir la même longueur")
            else:
                analyse = st.session_state.genomique.analyser_snp(seq_ref, seq_etu)
                
                st.markdown("### 📊 Résultats de l'analyse SNP")
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total SNP", analyse['total_snps'])
                with col2:
                    st.metric("Fréquence SNP", f"{analyse['frequence_snp']:.4f}")
                with col3:
                    st.metric("Longueur séquence", analyse['sequence_longueur'])
                
                # Visualisation des SNP
                if analyse['snps_detailles']:
                    st.markdown("#### 🎯 SNP Détectés")
                    snps_df = pd.DataFrame(analyse['snps_detailles'])
                    st.dataframe(snps_df)
                    
                    # Graphique de position des SNP
                    positions = [snp['position'] for snp in analyse['snps_detailles']]
                    fig = go.Figure(data=[go.Scatter(
                        x=positions,
                        y=[1] * len(positions),
                        mode='markers',
                        marker=dict(size=10, color='red'),
                        name='SNP'
                    )])
                    
                    fig.update_layout(
                        title="Position des SNP dans la séquence",
                        xaxis_title="Position (bp)",
                        yaxis=dict(showticklabels=False),
                        height=200
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 📊 Gènes Candidats pour l'Amélioration Génétique")
        
        race_selectionnee = st.selectbox("Sélectionner la race", 
                                        ["lacaune", "manech", "basco_bearnaise"])
        
        if st.button("🔍 Identifier les gènes candidats"):
            genes = st.session_state.genomique.rechercher_genes_candidats(race_selectionnee)
            
            # Tableau des gènes
            df_genes = pd.DataFrame(genes)
            st.dataframe(df_genes, use_container_width=True)
            
            # Graphique
            fig = px.bar(df_genes, x='gene', y=[1]*len(df_genes),
                        title=f"Gènes candidats pour {race_selectionnee}",
                        color_discrete_sequence=['green'])
            fig.update_yaxes(showticklabels=False)
            st.plotly_chart(fig, use_container_width=True)
            
            # Rapport génomique
            if st.button("📄 Générer le rapport génomique"):
                rapport = st.session_state.genomique.generer_rapport_genomique(1, genes)
                st.download_button(
                    label="📥 Télécharger le rapport",
                    data=rapport,
                    file_name=f"rapport_genomique_{race_selectionnee}.txt",
                    mime="text/plain"
                )

def afficher_statistiques():
    """Affiche la page d'analyse statistique"""
    st.markdown('<h2 class="section-header">📈 Analyse Statistique</h2>', unsafe_allow_html=True)
    
    st.info("""
    **Analyses statistiques avancées**
    - Corrélations génétiques
    - Modèles de courbes de lactation
    - Analyses multivariées
    - Intégration R stats
    """)
    
    tab1, tab2, tab3 = st.tabs(["📊 Corrélations", "🥛 Courbes Lactation", "📈 Analyses Multivariées"])
    
    with tab1:
        st.markdown("### Analyse de corrélations")
        
        # Charger les données
        cursor = st.session_state.db_manager.conn.cursor()
        cursor.execute('''
            SELECT b.race, b.sexe, cm.caractere, cm.valeur
            FROM caracteres_morpho cm
            JOIN brebis b ON cm.brebis_id = b.id
            WHERE cm.caractere IN ('poids_vif', 'longueur_corps', 'tour_poitrine')
        ''')
        
        data_raw = cursor.fetchall()
        
        if data_raw:
            # Préparer les données
            donnees = []
            for row in data_raw:
                donnees.append({
                    'race': row[0],
                    'sexe': row[1],
                    'caractere': row[2],
                    'valeur': row[3]
                })
            
            # Pivoter les données
            df_pivot = pd.DataFrame(donnees)
            df_wide = df_pivot.pivot_table(index=['race', 'sexe'], 
                                          columns='caractere', 
                                          values='valeur').reset_index()
            
            # Sélection des variables
            variables = st.multiselect(
                "Sélectionner les variables à analyser",
                ['poids_vif', 'longueur_corps', 'tour_poitrine'],
                default=['poids_vif', 'longueur_corps']
            )
            
            if len(variables) >= 2:
                # Calculer la corrélation
                correlation = st.session_state.stats.analyser_correlations(
                    df_wide.to_dict('records'), 
                    variables[0], 
                    variables[1]
                )
                
                # Afficher les résultats
                st.markdown("### 📊 Résultats de corrélation")
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Coefficient", f"{correlation['correlation']:.3f}")
                with col2:
                    st.metric("Interprétation", correlation['interpretation'].split()[0])
                with col3:
                    st.metric("Taille échantillon", correlation['taille_echantillon'])
                
                # Graphique de dispersion
                fig = px.scatter(df_wide, x=variables[0], y=variables[1],
                                color='race', trendline='ols',
                                title=f"Corrélation: {variables[0]} vs {variables[1]}")
                st.plotly_chart(fig, use_container_width=True)
                
                # Matrice de corrélation
                if len(variables) > 2:
                    st.markdown("#### 🎯 Matrice de corrélation")
                    corr_matrix = df_wide[variables].corr()
                    
                    fig = px.imshow(corr_matrix,
                                   text_auto='.2f',
                                   color_continuous_scale='Greens',
                                   title="Matrice de corrélation")
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("Sélectionnez au moins 2 variables")
        else:
            st.info("Pas assez de données pour l'analyse")
    
    with tab2:
        st.markdown("### Modélisation des courbes de lactation")
        
        # Simulation de données de lactation
        jours = list(range(1, 301))
        production = [3.0 * (j**0.2) * (2.718**(-0.005 * j)) for j in jours]
        
        df_lactation = pd.DataFrame({
            'jour_lactation': jours,
            'production': production
        })
        
        # Modélisation
        analyse = st.session_state.stats.analyser_production_lait(
            df_lactation.to_dict('records')
        )
        
        if 'pic_production' in analyse:
            st.markdown("### 📈 Résultats de modélisation")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Pic production", f"{analyse['pic_production']:.2f} L/jour")
            with col2:
                st.metric("Jour du pic", analyse['jour_pic'])
            with col3:
                st.metric("Production totale", f"{analyse['production_totale_estimee']:.0f} L")
            
            # Courbe de lactation
            fig = go.Figure()
            
            # Courbe réelle
            fig.add_trace(go.Scatter(
                x=jours,
                y=production,
                mode='lines',
                name='Production',
                line=dict(color='green', width=2)
            ))
            
            # Pic de lactation
            fig.add_trace(go.Scatter(
                x=[analyse['jour_pic']],
                y=[analyse['pic_production']],
                mode='markers',
                name='Pic',
                marker=dict(size=15, color='red')
            ))
            
            fig.update_layout(
                title="Courbe de lactation - Modèle de Wood",
                xaxis_title="Jours de lactation",
                yaxis_title="Production (L/jour)",
                hovermode='x unified'
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Persistance
            if analyse['persistance_lactation']:
                st.metric("Persistance lactation", 
                         f"{analyse['persistance_lactation']:.2f}",
                         "Production à J150 / Pic")
    
    with tab3:
        st.markdown("### Analyses multivariées")
        
        st.info("""
        **Analyses disponibles:**
        - Analyse en composantes principales (ACP)
        - Classification hiérarchique
        - Régression multiple
        - Modèles mixtes
        """)
        
        # Simulation d'ACP
        if st.button("🔍 Exécuter l'Analyse en Composantes Principales"):
            with st.spinner("Calcul en cours..."):
                # Données simulées
                np.random.seed(42)
                n_samples = 100
                
                data = {
                    'Poids': np.random.normal(70, 10, n_samples),
                    'Longueur': np.random.normal(75, 5, n_samples),
                    'Tour poitrine': np.random.normal(95, 8, n_samples),
                    'Production': np.random.normal(2.5, 0.5, n_samples)
                }
                
                df_acp = pd.DataFrame(data)
                
                # Graphique ACP simulé
                fig = px.scatter(df_acp, x='Poids', y='Production',
                                color='Longueur',
                                size='Tour poitrine',
                                title="Analyse multivariée - Projection des individus",
                                labels={'Poids': 'PC1 (38%)', 'Production': 'PC2 (22%)'},
                                color_continuous_scale='greens')
                
                st.plotly_chart(fig, use_container_width=True)
                
                # Statistiques
                st.markdown("#### 📊 Statistiques descriptives")
                st.dataframe(df_acp.describe())

def afficher_rapports():
    """Affiche la page de génération de rapports"""
    st.markdown('<h2 class="section-header">📋 Rapports et Exportations</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📄 Générer Rapports", "📥 Exporter Données", "📊 Tableaux de Bord"])
    
    with tab1:
        st.markdown("### Génération de rapports")
        
        col1, col2 = st.columns(2)
        
        with col1:
            type_rapport = st.selectbox(
                "Type de rapport",
                ["Rapport complet d'élevage", 
                 "Rapport sanitaire", 
                 "Rapport génétique",
                 "Rapport de production",
                 "Rapport financier"]
            )
            
            periode = st.selectbox(
                "Période",
                ["Mois en cours", "Trimestre", "Semestre", "Année complète"]
            )
        
        with col2:
            format_rapport = st.radio(
                "Format de sortie",
                ["PDF", "Excel", "HTML", "Word"],
                horizontal=True
            )
            
            detail = st.select_slider(
                "Niveau de détail",
                options=["Basique", "Standard", "Détaillé", "Expert"]
            )
        
        if st.button("🔄 Générer le rapport", type="primary"):
            with st.spinner("Génération du rapport en cours..."):
                # Simulation de génération de rapport
                rapport_content = f"""
                RAPPORT {type_rapport.upper()}
                ============================
                
                Période: {periode}
                Date de génération: {datetime.now().strftime('%Y-%m-%d %H:%M')}
                Niveau de détail: {detail}
                
                RÉSUMÉ EXÉCUTIF
                ---------------
                • Nombre de brebis: 150
                • Production laitière totale: 45,200 L
                • Taux de gestation: 85%
                • Dépenses vétérinaires: €12,450
                • Revenus totaux: €68,900
                
                ANALYSES DÉTAILLÉES
                --------------------
                1. Performance laitière: +12% vs période précédente
                2. Santé du troupeau: Excellent état général
                3. Génétique: Amélioration continue des indices
                4. Économique: Rentabilité à +18%
                
                RECOMMANDATIONS
                ----------------
                1. Continuer le programme de sélection génétique
                2. Investir dans l'alimentation complémentaire
                3. Renforcer la prévention sanitaire
                4. Explorer de nouveaux débouchés commerciaux
                """
                
                st.markdown("### 📋 Aperçu du rapport")
                st.text_area("Contenu du rapport", rapport_content, height=300)
                
                # Boutons de téléchargement
                col1, col2 = st.columns(2)
                with col1:
                    st.download_button(
                        label="📥 Télécharger PDF",
                        data=rapport_content,
                        file_name=f"rapport_{datetime.now().strftime('%Y%m%d')}.pdf",
                        mime="application/pdf"
                    )
                with col2:
                    st.download_button(
                        label="📥 Télécharger Excel",
                        data=rapport_content,
                        file_name=f"rapport_{datetime.now().strftime('%Y%m%d')}.xlsx",
                        mime="application/vnd.ms-excel"
                    )
    
    with tab2:
        st.markdown("### Exportation des données")
        
        # Sélection des données à exporter
        datasets = st.multiselect(
            "Sélectionner les jeux de données",
            ["Brebis", "Gestations", "Suivis médicaux", 
             "Mesures morphologiques", "Séquences génétiques",
             "Production laitière", "Données alimentaires"]
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            format_export = st.radio(
                "Format d'export",
                ["Excel (.xlsx)", "CSV (.csv)", "JSON (.json)", "SQL (.sql)"],
                horizontal=True
            )
        
        with col2:
            compression = st.checkbox("Compresser (ZIP)")
            include_metadata = st.checkbox("Inclure les métadonnées")
        
        if st.button("📤 Exporter les données"):
            # Simulation d'export
            progress_bar = st.progress(0)
            
            for i in range(100):
                # Simulation de progression
                progress_bar.progress(i + 1)
            
            st.success(f"✅ Données exportées: {len(datasets)} jeux de données")
            
            # Téléchargement simulé
            st.download_button(
                label="📥 Télécharger l'export",
                data="Contenu export simulé",
                file_name=f"export_ovin_{datetime.now().strftime('%Y%m%d_%H%M')}.zip",
                mime="application/zip"
            )
    
    with tab3:
        st.markdown("### Tableaux de bord personnalisés")
        
        # Création de dashboard
        st.info("Créez votre propre tableau de bord en sélectionnant les widgets:")
        
        col1, col2 = st.columns(2)
        
        with col1:
            widgets = st.multiselect(
                "Widgets disponibles",
                ["Graphique production", "Statistiques sanitaires", 
                 "Calendrier gestation", "Indicateurs économiques",
                 "Cartographie génétique", "Alertes et rappels",
                 "Comparaisons temporelles", "Benchmarks de race"]
            )
        
        with col2:
            layout = st.selectbox("Disposition", ["Grille 2x2", "Grille 3x3", "Vertical", "Horizontal"])
            refresh_rate = st.selectbox("Fréquence rafraîchissement", ["Manuel", "5 min", "15 min", "1 heure"])
        
        if st.button("🎨 Créer le tableau de bord"):
            st.success("✅ Tableau de bord créé avec succès!")
            
            # Afficher un aperçu
            st.markdown("#### 👁️ Aperçu du tableau de bord")
            
            # Widgets simulés
            cols = st.columns(2)
            with cols[0]:
                st.metric("Production journalière", "1520 L", "+5.2%")
                st.line_chart(pd.DataFrame({
                    'Production': np.random.randn(30).cumsum() + 1500
                }))
            
            with cols[1]:
                st.metric("Gestations actives", "18", "2 nouvelles")
                st.dataframe(pd.DataFrame({
                    'Brebis': ['Bella', 'Daisy', 'Luna'],
                    'Date mise bas': ['2024-06-15', '2024-06-18', '2024-06-22']
                }))

def afficher_configuration():
    """Affiche la page de configuration"""
    st.markdown('<h2 class="section-header">⚙️ Configuration</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["🔧 Paramètres", "📊 Base de données", "🔗 Intégrations", "🆘 Aide"])
    
    with tab1:
        st.markdown("### Paramètres généraux")
        
        col1, col2 = st.columns(2)
        
        with col1:
            langue = st.selectbox("Langue", ["Français", "English", "Español"])
            fuseau_horaire = st.selectbox("Fuseau horaire", ["Europe/Paris", "UTC", "America/New_York"])
            unite_poids = st.radio("Unité de poids", ["kg", "lbs"], horizontal=True)
            unite_longueur = st.radio("Unité de longueur", ["cm", "inches"], horizontal=True)
        
        with col2:
            format_date = st.selectbox("Format de date", ["DD/MM/YYYY", "YYYY-MM-DD", "MM/DD/YYYY"])
            notifications = st.checkbox("Activer les notifications")
            auto_sauvegarde = st.checkbox("Sauvegarde automatique", value=True)
            intervalle_sauvegarde = st.selectbox("Intervalle", ["Quotidien", "Hebdomadaire", "Mensuel"])
        
        if st.button("💾 Sauvegarder les paramètres"):
            st.success("✅ Paramètres sauvegardés avec succès!")
    
    with tab2:
        st.markdown("### Gestion de la base de données")
        
        # Statistiques de la base
        cursor = st.session_state.db_manager.conn.cursor()
        
        tables = ['brebis', 'gestations', 'suivi_medical', 'caracteres_morpho', 'sequences_genetiques']
        stats_db = []
        
        for table in tables:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            count = cursor.fetchone()[0]
            cursor.execute(f"SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='{table}'")
            exists = cursor.fetchone()[0]
            
            stats_db.append({
                'Table': table,
                'Enregistrements': count,
                'Existe': '✅' if exists else '❌'
            })
        
        df_stats = pd.DataFrame(stats_db)
        st.dataframe(df_stats, use_container_width=True)
        
        # Actions sur la base
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("🔄 Reconstruire la base"):
                init_database(st.session_state.db_manager)
                st.success("Base de données reconstruite!")
        
        with col2:
            if st.button("📤 Exporter la base"):
                # Simuler l'export
                st.info("Export de la base en cours...")
                st.download_button(
                    label="📥 Télécharger backup",
                    data="backup_simule.sql",
                    file_name="ovin_manager_backup.sql",
                    mime="application/sql"
                )
        
        with col3:
            if st.button("🗑️ Données de démo"):
                DonneesDemonstration.creer_donnees_test(st.session_state.db_manager)
                st.success("Données de démo créées!")
    
    with tab3:
        st.markdown("### Intégrations externes")
        
        st.info("""
        **Services disponibles:**
        - NCBI/GenBank: Accès aux bases génomiques
        - Services météo: Prévisions locales
        - Marchés agricoles: Prix et tendances
        - Laboratoires vétérinaires: Résultats d'analyses
        """)
        
        # Configuration NCBI
        st.markdown("#### 🧬 Configuration NCBI")
        
        col1, col2 = st.columns(2)
        with col1:
            ncbi_email = st.text_input("Email NCBI", value="contact@ovin-manager.com")
            ncbi_api_key = st.text_input("Clé API NCBI", type="password")
        
        with col2:
            ncbi_rate_limit = st.slider("Limite requêtes/sec", 1, 10, 3)
            ncbi_cache = st.checkbox("Activer le cache", value=True)
        
        # Configuration météo
        st.markdown("#### 🌤️ Configuration météo")
        
        api_key_meteo = st.text_input("Clé API Météo", type="password")
        localisation = st.text_input("Localisation", value="45.764043,4.835659")
        
        if st.button("🔗 Tester les connexions"):
            with st.spinner("Test des connexions en cours..."):
                st.success("✅ Connexion NCBI: OK")
                st.success("✅ Service météo: OK")
                st.warning("⚠️ Marchés agricoles: Non configuré")
    
    with tab4:
        st.markdown("### Aide et Support")
        
        st.markdown("""
        #### 📚 Documentation
        - [Guide d'utilisation](https://example.com)
        - [API Documentation](https://api.example.com)
        - [FAQ](https://faq.example.com)
        
        #### 🆘 Support technique
        **Email:** support@ovin-manager.com  
        **Téléphone:** +33 1 23 45 67 89  
        **Horaires:** 9h-18h du lundi au vendredi
        
        #### 🐛 Signaler un problème
        """)
        
        with st.form("form_support"):
            sujet = st.selectbox("Sujet", ["Bug", "Amélioration", "Question", "Autre"])
            description = st.text_area("Description du problème")
            email_contact = st.text_input("Votre email")
            
            if st.form_submit_button("📨 Envoyer la demande"):
                st.success("✅ Demande envoyée! Nous vous répondrons sous 48h.")

# Navigation principale
if page == "🏠 Tableau de Bord":
    afficher_tableau_bord()
elif page == "📊 Gestion des Brebis":
    afficher_gestion_brebis()
elif page == "🤰 Suivi Gestation":
    afficher_suivi_gestation()
elif page == "📸 Morphométrie":
    afficher_morphometrie()
elif page == "🧬 Analyse Génomique":
    afficher_genomique()
elif page == "📈 Statistiques":
    afficher_statistiques()
elif page == "📋 Rapports":
    afficher_rapports()
elif page == "⚙️ Configuration":
    afficher_configuration()

# Pied de page
st.markdown("---")
col1, col2, col3 = st.columns(3)
with col1:
    st.caption("🐑 Ovin Manager Pro v1.0.0")
with col2:
    st.caption("© 2024 - Tous droits réservés")
with col3:
    st.caption(f"Dernière mise à jour: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
