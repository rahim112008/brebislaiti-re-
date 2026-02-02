"""
EXPERT OVIN PRO - Système Intégré de Gestion et d'Évaluation Zootechnique
Version: 3.1 (Correction SQL + Scanner 1m)
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

class ConstantesReproduction:
    DUREE_CYCLE_ESTRAL = 17
    DUREE_GESTATION = 150
    DUREE_PROTOCOLE_EPG = 14
    DUREE_EFFET_EPG = 48
    SCORES_CORP_BCS = {'maigre': (1.0, 2.0), 'optimal': (2.5, 3.5), 'surgras': (4.0, 5.0)}

class GenesLaitiers(Enum):
    DGAT1 = {"chrom": "OAR9", "desc": "Diacylglycerol acyltransferase", "impact": "Matières grasses +15-20%"}
    LALBA = {"chrom": "OAR3", "desc": "Alpha-lactalbumine", "impact": "Quantité et qualité protéines"}
    CSN1S1 = {"chrom": "OAR3", "desc": "Caséine alpha-s1", "impact": "Rendement fromager"}
    CSN3 = {"chrom": "OAR6", "desc": "Caséine kappa", "impact": "Coagulation et texture"}

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
    """Initialisation avec correction des commentaires SQL (# -> --)"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        
        # 1. TABLE ANIMAUX
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY,
                numero_boucle TEXT UNIQUE,
                nom TEXT,
                espece TEXT CHECK(espece IN ('Bélier', 'Brebis', 'Agneau/elle')),
                race TEXT NOT NULL,
                date_naissance DATE,
                statut_reproductif TEXT,
                bcs_score REAL,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # 2. TABLE PRODUCTION LAITIÈRE (CORRIGÉE : Suppression des #)
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

        # 3. AUTRES TABLES (Suivi, Morpho, Génétique...)
        cursor.execute('CREATE TABLE IF NOT EXISTS suivi_medical (id INTEGER PRIMARY KEY, animal_id TEXT, date_intervention DATE, type_intervention TEXT, produit TEXT, FOREIGN KEY(animal_id) REFERENCES animaux(id))')
        cursor.execute('CREATE TABLE IF NOT EXISTS suivi_reproductif (id INTEGER PRIMARY KEY, animal_id TEXT, date_eponge_pose DATE, date_mise_bas_prevue DATE, FOREIGN KEY(animal_id) REFERENCES animaux(id))')
        
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

# ==========================================
# CALCULS SCIENTIFIQUES & GENOMIQUE
# ==========================================
class CalculateurZootechnique:
    @staticmethod
    def indice_conformation(perimetre_thorax: float, canon: float, hauteur_garrot: float) -> float:
        if canon <= 0 or hauteur_garrot <= 0: return 0
        return round((perimetre_thorax / (canon * hauteur_garrot)) * 1000, 2)

class NCBIConnector:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    @staticmethod
    def search_snp(gene_symbol: str):
        return [{"rs_id": "rs12345", "gene": gene_symbol, "status": "Simulated NCBI Data"}]

# ==========================================
# MODULES INTERFACE
# ==========================================

def module_dashboard():
    st.title("🏠 Dashboard")
    with get_db_connection() as conn:
        count = pd.read_sql("SELECT COUNT(*) as total FROM animaux", conn).iloc[0]['total']
    st.metric("Total Animaux", count)

def module_morphometrie():
    """Scanner Morphométrique avec étalon de 1m"""
    st.title("📏 Scanner Morphométrique (Étalon 1m)")
    
    with get_db_connection() as conn:
        df_brebis = pd.read_sql("SELECT id FROM animaux", conn)
    
    if df_brebis.empty:
        st.warning("Veuillez d'abord ajouter un animal.")
        return

    animal_id = st.selectbox("Sélectionner l'animal", df_brebis['id'])
    uploaded_file = st.file_uploader("Charger la photo (avec règle de 1m visible)", type=['jpg', 'png'])

    if uploaded_file:
        img = Image.open(uploaded_file)
        st.image(img, caption="Référence pour mesures")
        
        col1, col2 = st.columns(2)
        with col1:
            st.info("🎯 **Calibration**")
            pix_etalon = st.number_input("Pixels correspondant à l'étalon (1m)", min_value=1, value=500)
            ratio = 100 / pix_etalon # cm/pixel
            
        with col2:
            st.info("📐 **Mesures (en pixels sur l'image)**")
            p_hauteur = st.number_input("Hauteur Garrot (px)", value=350)
            p_longueur = st.number_input("Longueur Corps (px)", value=450)
            
        h_cm = round(p_hauteur * ratio, 2)
        l_cm = round(p_longueur * ratio, 2)
        
        st.success(f"Résultats : Hauteur = **{h_cm} cm** | Longueur = **{l_cm} cm**")
        
        if st.button("Sauvegarder les mesures"):
            with get_db_connection() as conn:
                conn.execute("INSERT INTO mesures_morphometriques (animal_id, date_mesure, hauteur_garrot, longueur_corps) VALUES (?,?,?,?)",
                            (animal_id, date.today().isoformat(), h_cm, l_cm))
            st.success("Enregistré !")

def module_gestion_animaux():
    st.title("📋 Gestion des Animaux")
    with st.form("add"):
        c1, c2 = st.columns(2)
        id_a = c1.text_input("ID")
        num = c1.text_input("Boucle")
        race = c2.selectbox("Race", ["Lacaune", "Manech", "Assaf"])
        espece = c2.selectbox("Espèce", ["Brebis", "Bélier"])
        if st.form_submit_button("Ajouter"):
            with get_db_connection() as conn:
                conn.execute("INSERT INTO animaux (id, numero_boucle, race, espece) VALUES (?,?,?,?)", (id_a, num, race, espece))
            st.rerun()

# ==========================================
# MAIN
# ==========================================
def main():
    st.set_page_config(page_title="Expert Ovin Pro", layout="wide")
    init_database()
    
    menu = st.sidebar.radio("Navigation", 
        ["Dashboard", "Gestion Animaux", "Scanner Morpho", "Reproduction", "Génomique"])
    
    if menu == "Dashboard": module_dashboard()
    elif menu == "Gestion Animaux": module_gestion_animaux()
    elif menu == "Scanner Morpho": module_morphometrie()
    else: st.info(f"Module {menu} en attente de données...")

if __name__ == "__main__":
    main()
