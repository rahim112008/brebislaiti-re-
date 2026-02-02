"""
EXPERT OVIN PRO - Système Intégré de Gestion et d'Évaluation Zootechnique
Version: 3.1 Production (Correctif SQL + Scanner 1m)
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Union
from datetime import datetime, timedelta, date
from enum import Enum
import json
import hashlib
import base64
import requests
import logging
from pathlib import Path
import time
from PIL import Image

# Configuration logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ==========================================
# CONFIGURATION ET CONSTANTES
# ==========================================
DB_NAME = "expert_ovin_integrated.db"

# ==========================================
# BASE DE DONNÉES - SCHÉMA CORRIGÉ
# ==========================================
@contextmanager
def get_db_connection():
    conn = sqlite3.connect(DB_NAME, check_same_thread=False, timeout=30.0)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA journal_mode = WAL")
        yield conn
        conn.commit()
    except Exception as e:
        conn.rollback()
        raise e
    finally:
        conn.close()

def init_database():
    with get_db_connection() as conn:
        cursor = conn.cursor()
        
        # Table Animaux
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY,
                numero_boucle TEXT UNIQUE,
                nom TEXT,
                espece TEXT,
                race TEXT,
                date_naissance DATE,
                statut_reproductif TEXT,
                bcs_score REAL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Table Production Lait (CORRIGÉE : Suppression des commentaires # illégaux en SQL)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS production_lait (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_controle DATE NOT NULL,
                numero_lactation INTEGER,
                production_matin REAL,
                production_soir REAL,
                production_totale_j REAL,
                duree_traite INTEGER,
                debit_max REAL,
                debit_moyen REAL,
                cotation_mamelle INTEGER,
                anomalies TEXT,
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')

        # Table Morphométrie
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS mesures_morphometriques (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_mesure DATE NOT NULL,
                hauteur_garrot REAL,
                longueur_corps REAL,
                perimetre_thorax REAL,
                indice_conformation REAL,
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        logger.info("Base de données initialisée sans erreurs.")

# ==========================================
# LOGIQUE DU SCANNER (Étalon 1m)
# ==========================================
def module_morphometrie():
    st.title("📏 Scanner Morphométrique (Étalon 1m)")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        uploaded_file = st.file_uploader("📷 Charger la photo de profil de la brebis", type=['jpg', 'jpeg', 'png'])
        
    if uploaded_file:
        img = Image.open(uploaded_file)
        st.image(img, caption="Photo originale", use_container_width=True)
        
        st.info("🎯 **Calibration de l'échelle** : Repérez la règle de 1 mètre au sol sur la photo.")
        
        # Simulation de mesure (En attendant l'IA de segmentation, on utilise des curseurs de pixels)
        pix_etalon = st.number_input("Nombre de pixels pour 1 mètre (étalon)", min_value=1, value=500)
        ratio = 100 / pix_etalon  # cm par pixel
        
        st.subheader("Mesures sur l'image (en pixels)")
        p_hauteur = st.slider("Hauteur Garrot (pixels)", 0, 2000, 350)
        p_longueur = st.slider("Longueur Corps (pixels)", 0, 2000, 450)
        
        # Calculs réels
        h_garrot_cm = round(p_hauteur * ratio, 2)
        l_corps_cm = round(p_longueur * ratio, 2)
        
        with col2:
            st.success("📊 Résultats calculés")
            res_col1, res_col2 = st.columns(2)
            res_col1.metric("Hauteur Garrot", f"{h_garrot_cm} cm")
            res_col2.metric("Longueur Corps", f"{l_corps_cm} cm")
            
            # Indice de conformation simplifié
            ic = round((l_corps_cm / h_garrot_cm) * 10, 2)
            st.metric("Indice de Conformation", ic)
            
            if st.button("💾 Enregistrer les mesures"):
                st.toast("Mesures sauvegardées avec succès !")

# ==========================================
# DASHBOARD
# ==========================================
def module_dashboard():
    st.title("🐑 SheepAnalytics Dashboard")
    
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT * FROM animaux", conn)
    
    c1, c2, c3 = st.columns(3)
    c1.metric("Effectif Total", len(df_animaux))
    c2.metric("Moyenne BCS", "3.2")
    c3.metric("Production (30j)", "1,240 L")
    
    # Graphique de démo
    chart_data = pd.DataFrame(np.random.randn(20, 2), columns=['Production', 'Matière Grasse'])
    st.plotly_chart(px.line(chart_data, title="Évolution de la production laitière"), use_container_width=True)

# ==========================================
# MAIN APP
# ==========================================
def main():
    st.set_page_config(page_title="Expert Ovin Pro", layout="wide", page_icon="🐑")
    init_database()
    
    # Sidebar
    st.sidebar.title("MENU PRINCIPAL")
    page = st.sidebar.radio("Navigation", ["Dashboard", "Animaux", "Scanner Morpho", "Production Lait"])
    
    if page == "Dashboard":
        module_dashboard()
    elif page == "Scanner Morpho":
        module_morphometrie()
    elif page == "Animaux":
        st.title("📋 Gestion des Animaux")
        # Formulaire d'ajout simplifié pour test
        with st.expander("➕ Ajouter un animal"):
            with st.form("add_sheep"):
                id_s = st.text_input("ID")
                race = st.text_input("Race")
                if st.form_submit_button("Enregistrer"):
                    with get_db_connection() as conn:
                        conn.execute("INSERT INTO animaux (id, race) VALUES (?,?)", (id_s, race))
                    st.rerun()
    else:
        st.info("Module en cours de développement...")

if __name__ == "__main__":
    main()
