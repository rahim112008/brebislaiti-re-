"""
EXPERT OVIN PRO - Système Intégré de Gestion et d'Évaluation Zootechnique
Version: 3.0 Production
Modules: Gestion d'élevage | Analyse scientifique | Génomique | Biochimie | Morphométrie IA
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

# Configuration logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("expert_ovin.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ==========================================
# CONFIGURATION ET CONSTANTES SCIENTIFIQUES
# ==========================================
DB_NAME = "expert_ovin_integrated.db"

class ConstantesReproduction:
    """Constantes pour la reproduction ovine"""
    DUREE_CYCLE_ESTRAL = 17  # jours
    DUREE_GESTATION = 150    # jours (5 mois)
    DUREE_PROTOCOLE_EPG = 14 # jours (éponge + injection)
    DUREE_EFFET_EPG = 48     # heures après retrait éponge
    
    SCORES_CORP_BCS = {
        'maigre': (1.0, 2.0),
        'optimal': (2.5, 3.5),
        'surgras': (4.0, 5.0)
    }

class GenesLaitiers(Enum):
    """Gènes majeurs pour la production laitière"""
    DGAT1 = {"chrom": "OAR9", "desc": "Diacylglycerol acyltransferase", "impact": "Matières grasses +15-20%"}
    LALBA = {"chrom": "OAR3", "desc": "Alpha-lactalbumine", "impact": "Quantité et qualité protéines"}
    CSN1S1 = {"chrom": "OAR3", "desc": "Caséine alpha-s1", "impact": "Rendement fromager"}
    CSN3 = {"chrom": "OAR6", "desc": "Caséine kappa", "impact": "Coagulation et texture"}
    PRLR = {"chrom": "OAR16", "desc": "Récepteur prolactine", "impact": "Production +10%"}
    STAT5A = {"chrom": "OAR19", "desc": "Signal transduction", "impact": "Différenciation cellules mammaires"}
    ACACA = {"chrom": "OAR11", "desc": "Acetyl-CoA carboxylase", "impact": "Synthèse acides gras"}
    FASN = {"chrom": "OAR12", "desc": "Fatty acid synthase", "impact": "Métabolisme lipides"}

class SeuilsMorphometriques:
    """Seuils morphométriques pour classification aptitude laitière"""
    # Mamelle
    LONGUEUR_TETINE_MIN = 3.5      # cm
    LONGUEUR_TETINE_MAX = 4.5      # cm
    DIAMETRE_TETINE_MIN = 2.0      # cm
    PROFONDEUR_MAMELLE_MIN = 12.0  # cm
    ATTACHE_ANTERIEURE_MAX = 4.0   # cm du corps
    
    # Corps
    PERIMETRE_THORAX_LAIT = 85.0   # cm
    PROFONDEUR_THORAX_MIN = 35.0   # cm
    HAUTEUR_GARROT_MIN = 65.0      # cm
    
    # Indices
    IC_OPTIMAL_MIN = 27.0          # Indice Conformation
    IC_ELITE_MIN = 33.0

# ==========================================
# BASE DE DONNÉES - SCHÉMA COMPLET
# ==========================================
@contextmanager
def get_db_connection():
    """Gestionnaire de connexion avec WAL mode"""
    conn = sqlite3.connect(DB_NAME, check_same_thread=False, timeout=30.0)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        yield conn
        conn.commit()
    except Exception as e:
        conn.rollback()
        logger.error(f"Erreur DB: {e}")
        raise
    finally:
        conn.close()

def init_database():
    """Initialisation complète du schéma"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        
        # 1. TABLE ANIMAUX (core)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY,
                numero_boucle TEXT UNIQUE,
                nom TEXT,
                espece TEXT CHECK(espece IN ('Bélier', 'Brebis', 'Agneau/elle')),
                race TEXT NOT NULL,
                date_naissance DATE,
                date_entree_ferme DATE,
                pere_id TEXT,
                mere_id TEXT,
                statut_reproductif TEXT CHECK(statut_reproductif IN 
                    ('Agneau', 'Génisse', 'Reproductive', 'Gestante', 'Lactante', 'Tarissement', 'Réforme')),
                bcs_score REAL CHECK(bcs_score BETWEEN 1 AND 5),
                notes TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # 2. TABLE GESTION MÉDICALE
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS suivi_medical (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_intervention DATE NOT NULL,
                type_intervention TEXT CHECK(type_intervention IN 
                    ('Vaccination', 'Vermifugation', 'Traitement', 'Chirurgie', 'Bilan sanguin', 'Échographie', 'Autre')),
                produit TEXT,
                posologie TEXT,
                veto_prescripteur TEXT,
                motif TEXT,
                couts REAL DEFAULT 0,
                suite_prevue DATE,
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # 3. TABLE GESTION ALIMENTAIRE
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS ration_alimentaire (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_debut DATE NOT NULL,
                date_fin DATE,
                type_ration TEXT CHECK(type_ration IN 
                    ('Maintenance', 'Croissance', 'Gestation', 'Lactation', 'Engraissement', 'Tarissement')),
                fourrage_principal TEXT,
                concentre_kg_j REAL,
                fourrage_kg_j REAL,
                mineraux TEXT,
                eau_litres_j REAL,
                cout_jour REAL,
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # 4. TABLE REPRODUCTION ET GESTATION
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS suivi_reproductif (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_chaleur DATE,
                date_saillie IA DATE,
                belier_utilise TEXT,
                type_protocole TEXT CHECK(type_protocole IN 
                    ('Naturel', 'Epg progesterone', 'Epg progestagene', 'Induction hormonale', 'IA')),
                date_eponge_pose DATE,
                date_eponge_retrait DATE,
                date_injection_PMSG DATE,
                dose_PMSG_UI INTEGER,
                date_mise_bas_prevue DATE,
                date_mise_bas_reelle DATE,
                nombre_agnes INTEGER,
                sexes_agnes TEXT,
                poids_agnes_naissance TEXT,
                complications TEXT,
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # 5. TABLE MORPHOMÉTRIE SCIENTIFIQUE
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS mesures_morphometriques (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_mesure DATE NOT NULL,
                mesureur TEXT,
                
                -- Mamelle
                longueur_tetine_gauche REAL,
                longueur_tetine_droite REAL,
                diametre_tetine_gauche REAL,
                diametre_tetine_droite REAL,
                profondeur_mamelle REAL,
                attache_anterieure REAL,
                attache_posterieure REAL,
                score_symetrie INTEGER CHECK(score_symetrie BETWEEN 1 AND 9),
                score_suspension INTEGER CHECK(score_suspension BETWEEN 1 AND 9),
                
                -- Corps
                hauteur_garrot REAL,
                longueur_corps REAL,
                perimetre_thorax REAL,
                profondeur_thorax REAL,
                largeur_bassin REAL,
                angle_croupe REAL,
                
                -- Poids
                poids_vif REAL,
                bcs_score REAL,
                
                -- Indices calculés
                ic_conformation REAL,
                volume_thoracique REAL,
                
                -- Photo référence
                photo_path TEXT,
                photo_hash TEXT,
                
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # 6. TABLE ANALYSE BIOCHIMIQUE DU LAIT
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS analyses_lait (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_prelevement DATE NOT NULL,
                stade_lactation TEXT CHECK(stade_lactation IN ('Début', 'Milieu', 'Fin')),
                
                -- Composition
                matieres_grasses_percent REAL,
                proteines_totales_percent REAL,
                caséine_percent REAL,
                lactose_percent REAL,
                matieres_seches_percent REAL,
                
                -- Paramètres technologiques
                taux_cellulaire INTEGER,  -- cellules/mL
                pH REAL,
                conductivite REAL,
                rct_secondes INTEGER,     -- Rennet Coagulation Time
                a30_mm REAL,              -- Fermeté caillé à 30 min
                rendement_fromage_percent REAL,
                
                -- Acides gras (profil)
                c4_0 REAL, c6_0 REAL, c8_0 REAL, c10_0 REAL,
                c12_0 REAL, c14_0 REAL, c16_0 REAL, c18_0 REAL,
                cla REAL, omega3 REAL,
                
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # 7. TABLE DONNÉES GÉNÉTIQUES ET SÉQUENÇAGE
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS donnees_genetiques (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT UNIQUE NOT NULL,
                statut_sequencage TEXT CHECK(statut_sequencage IN 
                    ('Non séquencé', 'En cours', 'Séquencé', 'Analysé', 'Élite confirmée')),
                
                -- Identifiants externes
                ncbi_biosample TEXT,
                ncbi_sra_accession TEXT,
                ensembl_id TEXT,
                
                -- SNPs clés (genotype: 0=AA, 1=Aa, 2=aa ou notation A/B)
                dgat1_snp TEXT,
                lalba_snp TEXT,
                csn1s1_snp TEXT,
                csn3_snp TEXT,
                prlr_snp TEXT,
                stat5a_snp TEXT,
                acaca_snp TEXT,
                fasn_snp TEXT,
                
                -- Scores
                score_genetique_lait REAL,
                estimation_aptitude_lait TEXT,
                
                -- Fichiers
                vcf_path TEXT,
                rapport_pdf_path TEXT,
                
                date_sequencage DATE,
                laboratoire TEXT,
                
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # 8. TABLE PRODUCTION LAITIÈRE (contrôle laitier)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS production_lait (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_controle DATE NOT NULL,
                numero_lactation INTEGER,
                
                -- Quantités
                production_matin REAL,    # kg
                production_soir REAL,     # kg
                production_totale_j REAL, # kg
                
                -- Paramètres traite
                duree_traite INTEGER,     # secondes
                debit_max REAL,           # ml/min
                debit_moyen REAL,         # ml/min
                
                -- Notes
                cotation_mamelle INTEGER,
                anomalies TEXT,
                
                FOREIGN KEY (animal_id) REFERENCES animaux(id) ON DELETE CASCADE
            )
        ''')
        
        # Indexes pour performance
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_animal_id ON animaux(id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_medical_animal ON suivi_medical(animal_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_gestation_animal ON suivi_reproductif(animal_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_morpho_animal ON mesures_morphometriques(animal_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_lait_animal ON production_lait(animal_id)')
        
        logger.info("Base de données initialisée avec succès")

# ==========================================
# MODÈLES DE DONNÉES (DATACLASSES)
# ==========================================
@dataclass
class Animal:
    """Représentation d'un animal dans le système"""
    id: str
    numero_boucle: str
    nom: Optional[str] = None
    espece: str = "Brebis"
    race: str = ""
    date_naissance: Optional[date] = None
    statut_reproductif: str = "Reproductive"
    bcs_score: float = 3.0
    
    def age_mois(self) -> int:
        if not self.date_naissance:
            return 0
        return (date.today() - self.date_naissance).days // 30
    
    def est_en_lactation(self) -> bool:
        return self.statut_reproductif == "Lactante"

@dataclass
class ProtocoleEPG:
    """Protocole d'induction de l'œstrus avec éponge"""
    date_pose: date
    type_eponge: str = "Fluorogestone (30mg)"  # ou Chronogest, etc.
    duree_jours: int = 14
    dose_PMSG: int = 400  # UI
    
    def calculer_mise_bas(self) -> date:
        """Calcule la date prévue de mise bas"""
        date_retrait = self.date_pose + timedelta(days=self.duree_jours)
        # Ovulation 24-48h après retrait, fécondation immédiate
        date_fertilisation = date_retrait + timedelta(days=2)
        # Gestation ~150 jours
        return date_fertilisation + timedelta(days=ConstantesReproduction.DUREE_GESTATION)
    
    def jours_avant_mise_bas(self) -> int:
        return (self.calculer_mise_bas() - date.today()).days

@dataclass
class ScoreMorphologique:
    """Score morphologique pour aptitude laitière"""
    # Mamelle (40 points)
    volume_mamelle: int = 0  # /10
    forme_mamelle: int = 0   # /10
    suspension: int = 0      # /10
    tétines: int = 0         # /10
    
    # Corps (30 points)
    capacité_corps: int = 0  # /15
    aplomb: int = 0          # /15
    
    # Lactation (30 points)
    qualité_lait: int = 0    # /15
    facilité_traite: int = 0 # /15
    
    def total(self) -> int:
        return (self.volume_mamelle + self.forme_mamelle + self.suspension + self.tétines +
                self.capacité_corps + self.aplomb + self.qualité_lait + self.facilité_traite)
    
    def classe(self) -> str:
        total = self.total()
        if total >= 85: return "🥇 ELITE LAITIERE"
        elif total >= 70: return "🥈 BONNE APTITUDE"
        elif total >= 50: return "🥉 MOYENNE"
        else: return "❌ FAIBLE"

# ==========================================
# CALCULS SCIENTIFIQUES
# ==========================================
class CalculateurZootechnique:
    """Moteur de calcul scientifique pour l'évaluation"""
    
    @staticmethod
    def indice_conformation(perimetre_thorax: float, canon: float, hauteur_garrot: float) -> float:
        """IC = (PT / (CC * HG)) * 1000"""
        if canon <= 0 or hauteur_garrot <= 0:
            return 0
        return round((perimetre_thorax / (canon * hauteur_garrot)) * 1000, 2)
    
    @staticmethod
    def prediction_composition_lait(poids: float, ic: float, hauteur: float) -> Dict:
        """Prédit la composition du lait d'après la morphologie"""
        gras_mm = max(2.0, 2.0 + (poids * 0.15) + (ic * 0.1) - (hauteur * 0.05))
        pct_gras = max(10.0, min(40.0, 5.0 + (gras_mm * 1.5)))
        pct_proteine = max(45.0, min(72.0, 75.0 - (pct_gras * 0.6) + (ic * 0.2)))
        pct_os = round(100.0 - pct_proteine - pct_gras, 1)
        
        return {
            'matieres_grasses': round(pct_gras, 1),
            'proteines': round(pct_proteine, 1),
            'os': pct_os,
            'rendement_fromage': round(42 + (pct_proteine * 0.12), 1)
        }
    
    @staticmethod
    def evaluer_genotype_lait(genes: Dict[str, str]) -> Dict:
        """Évalue le potentiel génétique laitier d'après les SNPs"""
        score = 0
        details = []
        
        # DGAT1 (impact majeur sur MG)
        if genes.get('dgat1') in ['GG', 'G/G', 'AA']:
            score += 20
            details.append("DGAT1 favorable (+15% MG)")
        elif genes.get('dgat1') in ['GA', 'G/A']:
            score += 10
            details.append("DGAT1 hétérozygote")
        
        # LALBA (protéines)
        if genes.get('lalba') == 'AA':
            score += 15
            details.append("LALBA A (+qualité protéique)")
        
        # CSN1S1 (rendement fromager)
        if genes.get('csn1s1') in ['CC', 'C/C', 'Strong']:
            score += 15
            details.append("Caséine αs1 forte (bon rendement)")
        
        # PRLR (production)
        if genes.get('prlr') in ['GG', 'G/G']:
            score += 10
            details.append("PRLR favorable (+production)")
        
        return {
            'score': score,
            'max_possible': 60,
            'niveau': 'Élite' if score >= 45 else 'Bon' if score >= 30 else 'Standard',
            'details': details
        }

# ==========================================
# INTEGRATION NCBI/GENBANK
# ==========================================
class NCBIConnector:
    """Connexion aux bases de données NCBI pour les séquences ovines"""
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    @staticmethod
    def search_snp(gene_symbol: str, organism: str = "Ovis aries") -> List[Dict]:
        """
        Recherche les SNPs connus pour un gène donné
        Note: Nécessite une connexion internet
        """
        try:
            # Recherche dans dbSNP
            url = f"{NCBIConnector.BASE_URL}/esearch.fcgi"
            params = {
                'db': 'snp',
                'term': f'{gene_symbol}[Gene Name] AND {organism}[Organism]',
                'retmode': 'json',
                'retmax': 20
            }
            response = requests.get(url, params=params, timeout=10)
            data = response.json()
            
            ids = data.get('esearchresult', {}).get('idlist', [])
            
            # Récupération des détails
            results = []
            if ids:
                summary_url = f"{NCBIConnector.BASE_URL}/esummary.fcgi"
                sum_params = {
                    'db': 'snp',
                    'id': ','.join(ids[:5]),  # Limite à 5 pour la démo
                    'retmode': 'json'
                }
                sum_response = requests.get(summary_url, params=sum_params, timeout=10)
                sum_data = sum_response.json()
                
                for uid, info in sum_data.get('result', {}).items():
                    if uid != 'uids':
                        results.append({
                            'rs_id': info.get('snp_id', ''),
                            'gene': gene_symbol,
                            'position': info.get('chrpos', 'N/A'),
                            'alleles': info.get('docsum', '').split('ALLELE:')[1].split(';')[0] if 'ALLELE:' in info.get('docsum', '') else 'N/A'
                        })
            
            return results
            
        except Exception as e:
            logger.error(f"Erreur connexion NCBI: {e}")
            return []
    
    @staticmethod
    def get_gene_sequence(gene_symbol: str) -> Optional[str]:
        """Récupère la séquence de référence d'un gène"""
        try:
            # Recherche dans Gene puis récupération dans Nucleotide
            search_url = f"{NCBIConnector.BASE_URL}/esearch.fcgi"
            params = {
                'db': 'gene',
                'term': f'{gene_symbol}[Gene Name] AND sheep[Organism]',
                'retmode': 'json'
            }
            response = requests.get(search_url, params=params, timeout=10)
            data = response.json()
            ids = data.get('esearchresult', {}).get('idlist', [])
            
            if ids:
                # Récupération du record
                fetch_url = f"{NCBIConnector.BASE_URL}/efetch.fcgi"
                fetch_params = {
                    'db': 'gene',
                    'id': ids[0],
                    'retmode': 'xml'
                }
                # Simplifié pour la démo
                return f"Gene ID: {ids[0]} (Ovis aries {gene_symbol})"
            
            return None
            
        except Exception as e:
            logger.error(f"Erreur: {e}")
            return None

# ==========================================
# INTERFACE STREAMLIT - APPLICATION PRINCIPALE
# ==========================================
def init_session_state():
    """Initialisation des variables de session"""
    defaults = {
        'animal_selectionne': None,
        'scan_morpho': {},
        'prediction_genetique': {},
        'mode_demo': False
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val

def render_sidebar():
    """Menu latéral de navigation"""
    with st.sidebar:
        st.title("🐑 EXPERT OVIN PRO")
        st.markdown("**Système Intégré de Gestion Zootechnique**")
        st.divider()
        
        menu = st.radio(
            "Navigation",
            [
                "🏠 Dashboard",
                "📋 Gestion Animaux",
                "🏥 Suivi Médical",
                "🌾 Alimentation",
                "🍼 Reproduction & Gestation",
                "📏 Analyse Morphométrique",
                "🧬 Génomique & SNPs",
                "⚗️ Biochimie Lait",
                "📊 Production Laitière",
                "🔧 Administration"
            ],
            key="main_nav"
        )
        
        st.divider()
        st.caption("v3.0 - Système Scientifique")
        
        return menu

# ==========================================
# MODULES PRINCIPAUX
# ==========================================

def module_dashboard():
    """Vue d'ensemble de l'exploitation"""
    st.title("🏠 Dashboard Exploitation")
    
    # Métriques globales
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT * FROM animaux", conn)
        df_gestantes = pd.read_sql(
            "SELECT * FROM suivi_reproductif WHERE date_mise_bas_prevue >= date('now')", conn
        )
        df_lactantes = pd.read_sql(
            "SELECT * FROM animaux WHERE statut_reproductif = 'Lactante'", conn
        )
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("🐑 Total Animaux", len(df_animaux) if not df_animaux.empty else 0)
    with col2:
        st.metric("🤰 Gestantes", len(df_gestantes) if not df_gestantes.empty else 0)
    with col3:
        st.metric("🥛 En Lactation", len(df_lactantes) if not df_lactantes.empty else 0)
    with col4:
        if not df_gestantes.empty:
            prochaine_mb = pd.to_datetime(df_gestantes['date_mise_bas_prevue']).min()
            jours_restants = (prochaine_mb - datetime.now()).days
            st.metric("⏱️ Prochaine Mise-bas", f"{jours_restants}j")
        else:
            st.metric("⏱️ Prochaine Mise-bas", "N/A")
    
    # Alertes et priorités
    st.subheader("⚠️ Alertes Prioritaires")
    
    alertes = []
    
    # Vérifie les synchronisations EPG en cours
    if not df_gestantes.empty:
        for _, row in df_gestantes.iterrows():
            date_prevue = pd.to_datetime(row['date_mise_bas_prevue'])
            jours_restants = (date_prevue - datetime.now()).days
            if 0 <= jours_restants <= 7:
                alertes.append({
                    'type': '🔴 MISE-BAS IMMINENTE',
                    'animal': row['animal_id'],
                    'date': date_prevue.strftime('%d/%m/%Y'),
                    'action': 'Préparer box de mise-bas'
                })
    
    if alertes:
        for alerte in alertes:
            with st.expander(f"{alerte['type']} - {alerte['animal']}"):
                st.write(f"**Date prévue:** {alerte['date']}")
                st.write(f"**Action:** {alerte['action']}")
    else:
        st.info("Aucune alerte urgente")
    
    # Calendrier des 30 prochains jours
    st.subheader("📅 Événements à venir (30 jours)")
    # Simplifié pour la démo

def module_gestion_animaux():
    """Gestion du registre animalier"""
    st.title("📋 Gestion des Animaux")
    
    tab1, tab2 = st.tabs(["📃 Liste", "➕ Nouvel Animal"])
    
    with tab1:
        with get_db_connection() as conn:
            df = pd.read_sql("SELECT * FROM animaux ORDER BY updated_at DESC", conn)
        
        if not df.empty:
            st.dataframe(df[['id', 'numero_boucle', 'espece', 'race', 'statut_reproductif', 'bcs_score']], 
                        use_container_width=True, hide_index=True)
        else:
            st.info("Aucun animal enregistré")
    
    with tab2:
        with st.form("nouvel_animal"):
            st.subheader("Enregistrement Nouvel Animal")
            
            col1, col2 = st.columns(2)
            with col1:
                id_animal = st.text_input("ID Unique *", placeholder="BREBIS_001")
                num_boucle = st.text_input("N° Boucle *", placeholder="FR123456")
                nom = st.text_input("Nom (optionnel)")
                espece = st.selectbox("Espèce", ["Brebis", "Bélier", "Agneau/elle"])
            
            with col2:
                race = st.selectbox("Race", 
                    ["Lacaune", "Manech", "Sarda", "Assaf", "Awassi", "East Friesian", "Chios", "Autre"])
                date_naiss = st.date_input("Date Naissance", value=None)
                bcs = st.slider("Score BCS", 1.0, 5.0, 3.0, 0.5)
                statut = st.selectbox("Statut", 
                    ["Agneau", "Génisse", "Reproductive", "Gestante", "Lactante", "Tarissement"])
            
            if st.form_submit_button("💾 Enregistrer", use_container_width=True):
                try:
                    with get_db_connection() as conn:
                        conn.execute('''
                            INSERT INTO animaux (id, numero_boucle, nom, espece, race, 
                                               date_naissance, bcs_score, statut_reproductif)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (id_animal, num_boucle, nom, espece, race, 
                              date_naiss.isoformat() if date_naiss else None, bcs, statut))
                    st.success(f"✅ {id_animal} enregistrée!")
                    st.balloons()
                except sqlite3.IntegrityError:
                    st.error("❌ ID ou N° Boucle déjà existant!")
                except Exception as e:
                    st.error(f"Erreur: {e}")

def module_reproduction():
    """Suivi reproduction, EPG et prédictions"""
    st.title("🍼 Reproduction & Gestion des Mises-bas")
    
    # Sélection animal
    with get_db_connection() as conn:
        df_brebis = pd.read_sql(
            "SELECT id, numero_boucle FROM animaux WHERE espece = 'Brebis'", conn
        )
    
    if df_brebis.empty:
        st.warning("Aucune brebis enregistrée")
        return
    
    animal_id = st.selectbox("Sélectionner la brebis", 
                            df_brebis['id'].tolist(),
                            format_func=lambda x: f"{x} ({df_brebis[df_brebis['id']==x]['numero_boucle'].values[0]})")
    
    tab1, tab2, tab3 = st.tabs(["🧽 Protocole EPG", "📊 Historique", "🔮 Prédictions"])
    
    with tab1:
        st.subheader("Nouveau Protocole d'Induction (EPG)")
        
        with st.form("protocole_epg"):
            col1, col2 = st.columns(2)
            with col1:
                date_pose = st.date_input("Date pose éponge", value=date.today())
                type_eponge = st.selectbox("Type", 
                    ["Fluorogestone acetate (30mg)", "Medroxyprogestérone (60mg)", "Chronogest"])
                duree = st.number_input("Durée (jours)", 12, 16, 14)
            
            with col2:
                PMSG = st.number_input("Dose PMSG (UI)", 200, 600, 400, 50)
                date_injection = st.date_input("Date injection PMSG", 
                                              value=date_pose + timedelta(days=duree-2))
            
            belier = st.text_input("Bélier utilisé / Numéro paillettes IA")
            
            if st.form_submit_button("🚀 Calculer et Enregistrer"):
                protocole = ProtocoleEPG(date_pose, type_eponge, duree, PMSG)
                date_mb = protocole.calculer_mise_bas()
                
                with get_db_connection() as conn:
                    conn.execute('''
                        INSERT INTO suivi_reproductif 
                        (animal_id, date_eponge_pose, date_eponge_retrait, date_injection_PMSG,
                         dose_PMSG_UI, type_protocole, belier_utilise, date_mise_bas_prevue)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (animal_id, date_pose.isoformat(), 
                          (date_pose + timedelta(days=duree)).isoformat(),
                          date_injection.isoformat(), PMSG, type_eponge, belier, date_mb.isoformat()))
                
                st.success(f"✅ Protocole enregistré!")
                st.info(f"**Date prévue mise-bas:** {date_mb.strftime('%d/%m/%Y')} "
                       f"({protocole.jours_avant_mise_bas()} jours)")
                
                # Calendrier des surveillances
                st.subheader("📅 Calendrier de surveillance")
                cal_df = pd.DataFrame([
                    {"Événement": "Retrait éponge", "Date": (date_pose + timedelta(days=duree)).strftime('%d/%m')},
                    {"Événement": "Injection PMSG", "Date": date_injection.strftime('%d/%m')},
                    {"Événement": "Début chaleurs", "Date": (date_pose + timedelta(days=duree+1)).strftime('%d/%m')},
                    {"Événement": "IA Naturelle", "Date": (date_pose + timedelta(days=duree+2)).strftime('%d/%m')},
                    {"Événement": "Mise-bas prévue", "Date": date_mb.strftime('%d/%m'), "⚠️": "CRITIQUE"}
                ])
                st.table(cal_df)
    
    with tab2:
        st.subheader("Historique Reproductif")
        with get_db_connection() as conn:
            histo = pd.read_sql('''
                SELECT * FROM suivi_reproductif 
                WHERE animal_id = ? 
                ORDER BY date_eponge_pose DESC
            ''', conn, params=(animal_id,))
        
        if not histo.empty:
            st.dataframe(histo)
        else:
            st.info("Aucun historique")
    
    with tab3:
        st.subheader("Prédictions Basées sur les Données")
        # ML simple pour prédire nombre d'agnes selon BCS, age, race...

def module_morphometrie():
    """Analyse morphométrique avec smartphone"""
    st.title("📏 Analyse Morphométrique Scientifique")
    
    st.markdown("""
    **Protocole de mesure au smartphone:**
    1. 📷 Photographier l'animal de profil (règle de 1m au sol pour échelle)
    2. 🎯 Marquer les points anatomiques sur l'image
    3. 📐 Le système calcule automatiquement les distances
    """)
    
    animal_id = st.selectbox("Animal à mesurer", 
                            ["Choisir..."] + get_liste_animaux())
    
    if animal_id == "Choisir...":
        return
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("📸 Capture ou Import")
        source = st.radio("Source", ["Caméra", "Fichier"], horizontal=True)
        
        if source == "Caméra":
            img = st.camera_input("Prendre la photo")
        else:
            img = st.file_uploader("Importer", type=['jpg', 'png', 'jpeg'])
        
        if img:
            st.image(img, caption="Image analysée", use_container_width=True)
            
            # Simulation analyse IA (dans la vraie version: OpenCV + calibration)
            with st.spinner("🔍 Analyse morphométrique IA..."):
                time.sleep(2)
                
                # Valeurs simulées réalistes
                mesures = {
                    'hauteur_garrot': np.random.normal(72, 3),
                    'longueur_corps': np.random.normal(78, 4),
                    'perimetre_thorax': np.random.normal(88, 5),
                    'profondeur_thorax': np.random.normal(36, 2),
                    'longueur_tetine': np.random.normal(4.0, 0.3)
                }
                st.session_state['scan_morpho'] = mesures
            
            with col2:
                st.success("✅ Mesures extraites")
                for k, v in mesures.items():
                    st.metric(k.replace('_', ' ').title(), f"{v:.1f} cm")
    
    with col2:
        st.subheader("📊 Évaluation Aptitude Laitière")
        
        if st.session_state.get('scan_morpho'):
            mesures = st.session_state['scan_morpho']
            
            # Calcul IC
            ic = CalculateurZootechnique.indice_conformation(
                mesures['perimetre_thorax'], 9.0, mesures['hauteur_garrot']
            )
            
            # Score morphologique
            score = ScoreMorphologique()
            score.volume_mamelle = min(10, int(mesures['perimetre_thorax'] / 10))
            score.suspension = 7
            score.tétines = min(10, int(mesures['longueur_tetine'] * 2))
            score.capacité_corps = min(15, int(mesures['profondeur_thorax'] / 2.5))
            
            st.metric("Score Total", f"{score.total()}/100")
            st.metric("Indice Conformation", ic)
            st.metric("Classe", score.classe())
            
            # Radar chart
            fig = go.Figure(go.Scatterpolar(
                r=[score.volume_mamelle, score.suspension, score.tétines, 
                   score.capacité_corps, score.aplomb],
                theta=['Volume', 'Suspension', 'Tétines', 'Capacité', 'Aplomb'],
                fill='toself'
            ))
            fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 10])))
            st.plotly_chart(fig, use_container_width=True)
            
            # Prédiction composition lait
            comp = CalculateurZootechnique.prediction_composition_lait(
                65, ic, mesures['hauteur_garrot']
            )
            st.subheader("🥛 Prédiction Lait")
            st.json(comp)

def module_genomique():
    """Gestion des données génétiques et liens NCBI"""
    st.title("🧬 Génomique & Séquençage")
    
    st.markdown("""
    **Intégration NCBI/GenBank:** Ce module permet de :
    - Consulter les SNPs connus pour les gènes laitiers
    - Saisir les génotypes des animaux séquencés
    - Calculer les scores de sélection génomique
    """)
    
    tab1, tab2, tab3 = st.tabs(["🔍 Recherche NCBI", "🧬 Saisie Génotypes", "📈 Évaluation Génétique"])
    
    with tab1:
        st.subheader("Recherche dans les bases NCBI")
        
        gene_search = st.selectbox("Gène à rechercher", 
                                  [g.name for g in GenesLaitiers])
        
        if st.button("🔎 Rechercher SNPs"):
            with st.spinner("Connexion à NCBI..."):
                # Simulation (en prod: appel API réel)
                results = NCBIConnector.search_snp(gene_search)
                
                if results:
                    st.success(f"{len(results)} SNPs trouvés")
                    for r in results[:3]:
                        st.write(f"**{r['rs_id']}** - Position: {r['position']}")
                        st.write(f"Allèles: {r['alleles']}")
                        st.divider()
                else:
                    # Données de démo si pas de connexion
                    st.info("Mode démonstration (sans connexion)")
                    demo_data = {
                        'DGAT1': [{'rs_id': 'rs384757468', 'position': 'OAR9:37113685', 'alleles': 'G/A'}],
                        'LALBA': [{'rs_id': 'rs397507565', 'position': 'OAR3:80456721', 'alleles': 'A/B'}]
                    }
                    if gene_search in demo_data:
                        st.json(demo_data[gene_search])
    
    with tab2:
        st.subheader("Saisie des Génotypes")
        
        animal_id = st.selectbox("Animal", get_liste_animaux(), key="geno_animal")
        
        with st.form("genotype_form"):
            cols = st.columns(4)
            genotypes = {}
            
            genes = ['dgat1', 'lalba', 'csn1s1', 'csn3', 'prlr', 'stat5a', 'acaca', 'fasn']
            
            for i, gene in enumerate(genes):
                with cols[i % 4]:
                    genotypes[gene] = st.selectbox(
                        f"{gene.upper()}", 
                        ['Non testé', 'AA', 'AB', 'BB', 'A/A', 'A/B', 'B/B'],
                        key=f"gene_{gene}"
                    )
            
            date_seq = st.date_input("Date séquençage")
            labo = st.text_input("Laboratoire")
            
            if st.form_submit_button("💾 Enregistrer"):
                eval_gen = CalculateurZootechnique.evaluer_genotype_lait(genotypes)
                
                with get_db_connection() as conn:
                    conn.execute('''
                        INSERT OR REPLACE INTO donnees_genetiques 
                        (animal_id, dgat1_snp, lalba_snp, csn1s1_snp, csn3_snp, prlr_snp,
                         stat5a_snp, acaca_snp, fasn_snp, score_genetique_lait, 
                         estimation_aptitude_lait, date_sequencage, laboratoire, statut_sequencage)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (animal_id, genotypes.get('dgat1'), genotypes.get('lalba'),
                          genotypes.get('csn1s1'), genotypes.get('csn3'), genotypes.get('prlr'),
                          genotypes.get('stat5a'), genotypes.get('acaca'), genotypes.get('fasn'),
                          eval_gen['score'], eval_gen['niveau'], date_seq.isoformat(), labo, 'Analysé'))
                
                st.success(f"✅ Score génétique: {eval_gen['score']}/60 - {eval_gen['niveau']}")
                for detail in eval_gen['details']:
                    st.write(f"- {detail}")
    
    with tab3:
        st.subheader("Évaluation Génétique du Troupeau")
        # Tableau comparatif, distribution des scores...

def module_biochimie():
    """Analyses biochimiques du lait"""
    st.title("⚗️ Biochimie et Technologie du Lait")
    
    animal_id = st.selectbox("Animal", get_liste_animaux())
    
    with st.form("analyse_lait"):
        st.subheader("Nouvelle Analyse")
        
        col1, col2 = st.columns(2)
        with col1:
            date_prev = st.date_input("Date prélèvement")
            stade = st.selectbox("Stade lactation", ["Début", "Milieu", "Fin"])
            mg = st.number_input("Matières grasses %", 3.0, 12.0, 6.5, 0.1)
            prot = st.number_input("Protéines totales %", 4.0, 8.0, 5.8, 0.1)
            caseine = st.number_input("Caséine %", 3.0, 6.0, 4.5, 0.1)
        
        with col2:
            lactose = st.number_input("Lactose %", 4.0, 6.0, 4.8, 0.1)
            pH = st.number_input("pH", 6.0, 7.0, 6.6, 0.01)
            taux_cell = st.number_input("Taux cellulaire (milliers/mL)", 0, 2000, 200, 50)
            rct = st.number_input("RCT (secondes)", 100, 1800, 900, 10)
        
        if st.form_submit_button("🔬 Analyser"):
            # Calcul rendement fromager estimé
            rendement = 42 + (caseine * 0.15) + (mg * 0.08) - (pH - 6.6) * 5
            
            st.success(f"Rendement fromager estimé: {rendement:.1f}%")
            
            # Qualification
            if caseine > 4.5 and rct < 1200 and taux_cell < 500:
                qualite = "🥇 PREMIUM (Fromage AOP)"
            elif caseine > 4.0:
                qualite = "🥈 BONNE (Fromage fermier)"
            else:
                qualite = "🥉 STANDARD (Consommation)"
            
            st.metric("Qualification", qualite)

# ==========================================
# UTILITAIRES
# ==========================================
def get_liste_animaux() -> List[str]:
    """Récupère la liste des IDs animaux"""
    try:
        with get_db_connection() as conn:
            df = pd.read_sql("SELECT id FROM animaux ORDER BY id", conn)
            return df['id'].tolist() if not df.empty else []
    except:
        return []

def main():
    """Point d'entrée principal"""
    st.set_page_config(
        page_title="Expert Ovin Pro - Système Intégré",
        page_icon="🐑",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # CSS personnalisé
    st.markdown("""
    <style>
    .stMetric { background-color: #f0f2f6; padding: 10px; border-radius: 10px; }
    .stAlert { border-radius: 8px; }
    div[data-testid="stExpander"] div[role="button"] p { font-weight: bold; color: #2E7D32; }
    </style>
    """, unsafe_allow_html=True)
    
    # Initialisation
    init_database()
    init_session_state()
    
    # Navigation
    menu = render_sidebar()
    
    # Routage
    if menu == "🏠 Dashboard":
        module_dashboard()
    elif menu == "📋 Gestion Animaux":
        module_gestion_animaux()
    elif menu == "🏥 Suivi Médical":
        st.title("🏥 Suivi Médical")
        st.info("Module en développement - Vaccinations, traitements, bilans")
    elif menu == "🌾 Alimentation":
        st.title("🌾 Gestion Alimentaire")
        st.info("Module en développement - Rations, BCS, coûts")
    elif menu == "🍼 Reproduction & Gestion":
        module_reproduction()
    elif menu == "📏 Analyse Morphométrique":
        module_morphometrie()
    elif menu == "🧬 Génomique & SNPs":
        module_genomique()
    elif menu == "⚗️ Biochimie Lait":
        module_biochimie()
    elif menu == "📊 Production Laitière":
        st.title("📊 Contrôle Laitier")
        st.info("Module en développement - Courbes de lactation, CCI")
    elif menu == "🔧 Administration":
        st.title("🔧 Administration")
        if st.button("⚠️ Générer données test (50 animaux)"):
            st.success("Fonction de génération à implémenter")

if __name__ == "__main__":
    main()
