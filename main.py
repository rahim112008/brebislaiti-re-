"""
EXPERT OVIN PRO - Système Intégré de Gestion et d'Évaluation Zootechnique
Version: 3.2 (Correction SQL + Données de démo + Scanner 1m)
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta, date
from enum import Enum
import logging
from PIL import Image

# Configuration logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ==========================================
# CONFIGURATION ET CONSTANTES SCIENTIFIQUES
# ==========================================
DB_NAME = "expert_ovin_integrated.db"

class ConstantesReproduction:
    DUREE_CYCLE_ESTRAL = 17
    DUREE_GESTATION = 150
    DUREE_PROTOCOLE_EPG = 14
    SCORES_CORP_BCS = {'maigre': (1.0, 2.0), 'optimal': (2.5, 3.5), 'surgras': (4.0, 5.0)}

class GenesLaitiers(Enum):
    DGAT1 = {"chrom": "OAR9", "impact": "Matières grasses +15-20%"}
    LALBA = {"chrom": "OAR3", "impact": "Qualité protéines"}

# ==========================================
# BASE DE DONNÉES - SCHÉMA ET SEEDING
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
    """Initialisation avec correction des erreurs de syntaxe SQL"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        
        # Table Animaux
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY,
                numero_boucle TEXT UNIQUE,
                nom TEXT,
                espece TEXT,
                race TEXT NOT NULL,
                date_naissance DATE,
                statut_reproductif TEXT,
                bcs_score REAL,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Table Production Lait (CORRIGÉE : Suppression des commentaires # illégaux)
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
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')

def seed_data():
    """Ajoute des données de démonstration si la base est vide"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        # Vérifier si la base est vide
        cursor.execute("SELECT COUNT(*) FROM animaux")
        if cursor.fetchone()[0] == 0:
            # Insertion Brebis
            sample_sheep = [
                ('BR_001', 'FR12345', 'Belle', 'Brebis', 'Lacaune', '2023-01-15', 'Lactante', 3.5),
                ('BR_002', 'FR67890', 'Noiraude', 'Brebis', 'Manech', '2022-05-20', 'Gestante', 3.0)
            ]
            cursor.executemany("INSERT INTO animaux (id, numero_boucle, nom, espece, race, date_naissance, statut_reproductif, bcs_score) VALUES (?,?,?,?,?,?,?,?)", sample_sheep)
            
            # Insertion Production
            sample_prod = [
                ('BR_001', '2026-01-20', 1, 1.2, 1.0, 2.2, 300, 450, 400, 4, 'RAS'),
                ('BR_001', '2026-02-01', 1, 1.3, 1.1, 2.4, 280, 460, 410, 5, 'Excellente traite')
            ]
            cursor.executemany("INSERT INTO production_lait (animal_id, date_controle, numero_lactation, production_matin, production_soir, production_totale_j, duree_traite, debit_max, debit_moyen, cotation_mamelle, anomalies) VALUES (?,?,?,?,?,?,?,?,?,?,?)", sample_prod)
            logger.info("Données de démonstration injectées.")

# ==========================================
# MODULES INTERFACE
# ==========================================

def module_dashboard():
    st.title("🏠 Dashboard de l'Exploitation")
    
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT * FROM animaux", conn)
        df_prod = pd.read_sql("SELECT * FROM production_lait", conn)

    c1, c2, c3 = st.columns(3)
    c1.metric("Total Animaux", len(df_animaux))
    if not df_prod.empty:
        c2.metric("Moyenne Production (L)", round(df_prod['production_totale_j'].mean(), 2))
        c3.metric("Dernière Mesure", df_prod['date_controle'].max())

    if not df_prod.empty:
        st.subheader("📈 Évolution de la production")
        fig = px.line(df_prod, x='date_controle', y='production_totale_j', color='animal_id', markers=True)
        st.plotly_chart(fig, use_container_width=True)

def module_morphometrie():
    st.title("📏 Scanner Morphométrique (Standard 1m)")
    st.info("Utilisez une photo avec une règle de 1 mètre au sol pour calibrer les mesures.")
    
    with get_db_connection() as conn:
        df_brebis = pd.read_sql("SELECT id FROM animaux", conn)
    
    if df_brebis.empty:
        st.warning("Aucun animal disponible.")
        return

    animal_id = st.selectbox("Animal à mesurer", df_brebis['id'])
    file = st.file_uploader("Charger la photo de profil", type=['jpg', 'png'])

    if file:
        img = Image.open(file)
        st.image(img, caption="Aperçu pour analyse")
        
        col1, col2 = st.columns(2)
        pix_etalon = col1.number_input("Pixels pour 1 mètre (étalon)", min_value=1, value=500)
        ratio = 100 / pix_etalon
        
        p_hauteur = col2.number_input("Mesure hauteur (pixels)", value=400)
        h_cm = round(p_hauteur * ratio, 2)
        
        st.success(f"📏 Hauteur calculée : **{h_cm} cm**")
        
        if st.button("Enregistrer la mesure"):
            with get_db_connection() as conn:
                conn.execute("INSERT INTO mesures_morphometriques (animal_id, date_mesure, hauteur_garrot, longueur_corps) VALUES (?,?,?,?)",
                            (animal_id, date.today().isoformat(), h_cm, 0))
            st.success("Mesure sauvegardée en base !")

# ==========================================
# MAIN APP
# ==========================================
def main():
    st.set_page_config(page_title="Expert Ovin Pro", layout="wide", page_icon="🐑")
    
    # Initialisation
    init_database()
    seed_data()
    
    # Navigation
    with st.sidebar:
        st.title("🐑 EXPERT OVIN")
        menu = st.radio("Navigation", ["Dashboard", "Gestion Animaux", "Scanner Morpho", "Production Laitière"])
        st.divider()
        st.caption("Mode Démonstration Actif")

    if menu == "Dashboard":
        module_dashboard()
    elif menu == "Scanner Morpho":
        module_morphometrie()
    elif menu == "Gestion Animaux":
        st.title("📋 Registre Animalier")
        with get_db_connection() as conn:
            df = pd.read_sql("SELECT * FROM animaux", conn)
            st.dataframe(df, use_container_width=True)
    elif menu == "Production Laitière":
        st.title("🥛 Contrôle Laitier")
        with get_db_connection() as conn:
            df = pd.read_sql("SELECT * FROM production_lait", conn)
            st.dataframe(df, use_container_width=True)

if __name__ == "__main__":
    main()
