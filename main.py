"""
OVIN MANAGER PRO - Version simplifiée pour Streamlit Cloud
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import io

# Configuration de la page
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
    }
</style>
""", unsafe_allow_html=True)

# Initialisation de la base de données
def init_database():
    """Initialise la base de données SQLite"""
    conn = sqlite3.connect('ovin_manager.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table des brebis
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT,
            date_naissance DATE,
            race TEXT,
            sexe TEXT,
            statut TEXT DEFAULT 'active',
            notes TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Table gestations
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_eponge DATE,
            date_mise_bas_prevu DATE,
            nombre_agneaux_prevus INTEGER DEFAULT 1,
            statut TEXT DEFAULT 'en_cours',
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
    ''')
    
    conn.commit()
    return conn

# Connexion à la base de données
conn = init_database()

# Titre principal
st.markdown('<h1 class="main-header">🐑 Ovin Manager Pro</h1>', unsafe_allow_html=True)
st.markdown("""
*Application scientifique de gestion et d'analyse d'élevage ovin laitier*
""")

# Sidebar - Navigation
with st.sidebar:
    st.markdown("### 📍 Navigation")
    
    page = st.radio(
        "Menu Principal",
        ["🏠 Tableau de Bord", 
         "📊 Gestion des Brebis", 
         "🤰 Suivi Gestation", 
         "📊 Statistiques"]
    )
    
    st.markdown("---")
    st.markdown("### 📊 Statistiques rapides")
    
    # Statistiques rapides
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM brebis")
    total_brebis = cursor.fetchone()[0]
    st.metric("Total Brebis", total_brebis)
    
    cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
    gestations_en_cours = cursor.fetchone()[0]
    st.metric("Gestations en cours", gestations_en_cours)

# Fonction pour afficher le tableau de bord
def afficher_tableau_bord():
    st.markdown('<h2 class="section-header">🏠 Tableau de Bord</h2>', unsafe_allow_html=True)
    
    # Métriques
    col1, col2, col3 = st.columns(3)
    
    with col1:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM brebis")
        total = cursor.fetchone()[0]
        st.metric("Total Brebis", total)
    
    with col2:
        cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F'")
        femelles = cursor.fetchone()[0]
        st.metric("Brebis Femelles", femelles)
    
    with col3:
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
        gestations = cursor.fetchone()[0]
        st.metric("Gestations en cours", gestations)
    
    # Graphique
    st.markdown("### 📊 Répartition par Race")
    cursor.execute("SELECT race, COUNT(*) as count FROM brebis GROUP BY race")
    data_race = cursor.fetchall()
    
    if data_race:
        df_race = pd.DataFrame(data_race, columns=['Race', 'Nombre'])
        fig = px.pie(df_race, values='Nombre', names='Race', 
                    color_discrete_sequence=px.colors.sequential.Greens)
        st.plotly_chart(fig, use_container_width=True)

# Fonction pour gérer les brebis
def afficher_gestion_brebis():
    st.markdown('<h2 class="section-header">📊 Gestion des Brebis</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["📋 Liste des Brebis", "➕ Ajouter une Brebis"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM brebis ORDER BY id")
        brebis_data = cursor.fetchall()
        
        if brebis_data:
            columns = [desc[0] for desc in cursor.description]
            df_brebis = pd.DataFrame(brebis_data, columns=columns)
            st.dataframe(df_brebis[['id', 'nom', 'race', 'sexe', 'statut']])
    
    with tab2:
        with st.form("form_ajout_brebis"):
            identifiant = st.text_input("Identifiant Unique*")
            nom = st.text_input("Nom")
            date_naissance = st.date_input("Date de Naissance", value=date.today() - timedelta(days=365))
            race = st.selectbox("Race", ["lacaune", "manech", "basco_bearnaise"])
            sexe = st.radio("Sexe", ["F", "M"], horizontal=True)
            statut = st.selectbox("Statut", ["active", "retrait"])
            
            if st.form_submit_button("✅ Ajouter la Brebis"):
                if identifiant:
                    try:
                        cursor = conn.cursor()
                        cursor.execute('''
                            INSERT INTO brebis (identifiant_unique, nom, date_naissance, race, sexe, statut)
                            VALUES (?, ?, ?, ?, ?, ?)
                        ''', (identifiant, nom, date_naissance.isoformat(), race, sexe, statut))
                        conn.commit()
                        st.success(f"✅ Brebis {identifiant} ajoutée!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Erreur: {e}")

# Fonction pour le suivi des gestations
def afficher_suivi_gestation():
    st.markdown('<h2 class="section-header">🤰 Suivi des Gestations</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["📅 Calendrier", "➕ Nouvelle Gestation"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT g.*, b.nom, b.race 
            FROM gestations g 
            JOIN brebis b ON g.brebis_id = b.id 
            WHERE g.statut = 'en_cours'
        ''')
        
        gestations = cursor.fetchall()
        
        if gestations:
            for gestation in gestations:
                cols = ['id', 'brebis_id', 'date_eponge', 'date_mise_bas_prevu', 
                       'nombre_agneaux_prevus', 'statut', 'nom', 'race']
                gest_dict = dict(zip(cols, gestation))
                
                with st.expander(f"🐑 {gest_dict['nom']} - Mise bas: {gest_dict['date_mise_bas_prevu']}"):
                    st.write(f"**Race:** {gest_dict['race']}")
                    st.write(f"**Date éponge:** {gest_dict['date_eponge']}")
                    st.write(f"**Agneaux prévus:** {gest_dict['nombre_agneaux_prevus']}")
    
    with tab2:
        cursor = conn.cursor()
        cursor.execute("SELECT id, nom FROM brebis WHERE sexe = 'F' AND statut = 'active'")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            with st.form("form_gestation"):
                brebis_id = st.selectbox("Brebis", [f"{b[1]} (ID: {b[0]})" for b in brebis_list])
                date_eponge = st.date_input("Date d'éponge", value=date.today())
                nb_agneaux = st.number_input("Nombre d'agneaux prévus", min_value=1, max_value=4, value=1)
                
                if st.form_submit_button("📅 Ajouter la gestation"):
                    # Extraire l'ID de la brebis
                    brebis_id_num = int(brebis_id.split("ID: ")[1].rstrip(")"))
                    date_mise_bas = date_eponge + timedelta(days=150)
                    
                    cursor.execute('''
                        INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevu, nombre_agneaux_prevus)
                        VALUES (?, ?, ?, ?)
                    ''', (brebis_id_num, date_eponge.isoformat(), date_mise_bas.isoformat(), nb_agneaux))
                    conn.commit()
                    st.success("✅ Gestation ajoutée!")
                    st.rerun()

# Fonction pour les statistiques
def afficher_statistiques():
    st.markdown('<h2 class="section-header">📈 Statistiques</h2>', unsafe_allow_html=True)
    
    cursor = conn.cursor()
    
    # Statistiques générales
    cursor.execute("SELECT COUNT(*) FROM brebis")
    total = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F'")
    femelles = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
    gestations = cursor.fetchone()[0]
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total brebis", total)
    with col2:
        st.metric("Femelles", femelles)
    with col3:
        st.metric("Gestations", gestations)
    
    # Graphique de répartition
    st.markdown("### 📊 Répartition par statut")
    cursor.execute("SELECT statut, COUNT(*) as count FROM brebis GROUP BY statut")
    data_statut = cursor.fetchall()
    
    if data_statut:
        df_statut = pd.DataFrame(data_statut, columns=['Statut', 'Nombre'])
        fig = px.bar(df_statut, x='Statut', y='Nombre', 
                    color='Statut', color_discrete_sequence=px.colors.sequential.Greens)
        st.plotly_chart(fig, use_container_width=True)

# Navigation principale
if page == "🏠 Tableau de Bord":
    afficher_tableau_bord()
elif page == "📊 Gestion des Brebis":
    afficher_gestion_brebis()
elif page == "🤰 Suivi Gestation":
    afficher_suivi_gestation()
elif page == "📊 Statistiques":
    afficher_statistiques()

# Pied de page
st.markdown("---")
st.caption("🐑 Ovin Manager Pro v1.0.0 | © 2024")
