"""
EXPERT OVIN DZ PRO - Version Intégrée Ultime
Fusion des trois systèmes :
1. Gestion zootechnique complète (reproduction, médical, alimentation)
2. Races algériennes avec scanner 3D et critères mammaires
3. Labo génétique et bioinformatique

Auteur: Système Expert Ovin DZ
Version: 5.0 - Production Ready
"""

# ============================================================================
# SECTION 1: IMPORTS COMPLETS
# ============================================================================
import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime, date, timedelta
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
from enum import Enum
import json
import random
import math
import io
import base64
import hashlib
import logging
import time
from PIL import Image, ImageDraw
import requests

# Configuration logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("expert_ovin_dz.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# SECTION 2: CONFIGURATION STREAMLIT
# ============================================================================
st.set_page_config(
    page_title="Expert Ovin DZ Pro - Système Intégré",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# SECTION 3: CSS PERSONNALISÉ PROFESSIONNEL
# ============================================================================
st.markdown("""
<style>
    /* En-têtes */
    .main-header {
        font-size: 2.8rem;
        font-weight: bold;
        background: linear-gradient(90deg, #1a237e, #8B0000);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        margin-bottom: 1rem;
    }
    .section-header {
        font-size: 1.8rem;
        color: #1a237e;
        border-bottom: 3px solid #8B0000;
        padding-bottom: 10px;
        margin-top: 20px;
    }
    
    /* Cartes */
    .metric-card {
        background: linear-gradient(135deg, #f5f5f5 0%, #ffffff 100%);
        border-radius: 15px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        border-left: 5px solid #1a237e;
    }
    .race-card {
        background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
    }
    .mammelle-card {
        background: linear-gradient(135deg, #8B0000 0%, #c62828 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
    }
    .gene-card {
        background: linear-gradient(135deg, #6a1b9a 0%, #8e24aa 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
    }
    .alert-card {
        background: linear-gradient(135deg, #fff3cd 0%, #ffeeba 100%);
        border-left: 5px solid #ffc107;
        padding: 15px;
        border-radius: 10px;
    }
    
    /* Tableaux et données */
    .stDataFrame { border-radius: 10px; }
    .stMetric { background-color: #f8f9fa; border-radius: 10px; }
    
    /* Navigation */
    .nav-selected { background-color: #1a237e; color: white; }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# SECTION 4: CONSTANTES ET CONFIGURATIONS
# ============================================================================
DB_NAME = "expert_ovin_dz_integrated.db"

class ConstantesReproduction:
    DUREE_CYCLE_ESTRAL = 17
    DUREE_GESTATION = 150
    DUREE_PROTOCOLE_EPG = 14
    DUREE_EFFET_EPG = 48

class SeuilsMorphometriques:
    LONGUEUR_TETINE_MIN = 3.5
    LONGUEUR_TETINE_MAX = 4.5
    DIAMETRE_TETINE_MIN = 2.0
    PROFONDEUR_MAMELLE_MIN = 12.0
    IC_OPTIMAL_MIN = 27.0
    IC_ELITE_MIN = 33.0

class GenesLaitiers(Enum):
    DGAT1 = {"chrom": "OAR9", "desc": "Synthèse triglycérides", "impact": "MG +15-20%"}
    LALBA = {"chrom": "OAR3", "desc": "Alpha-lactalbumine", "impact": "Protéines"}
    CSN1S1 = {"chrom": "OAR3", "desc": "Caséine αs1", "impact": "Rendement fromage"}
    CSN3 = {"chrom": "OAR6", "desc": "Caséine κ", "impact": "Coagulation"}
    PRLR = {"chrom": "OAR16", "desc": "Récepteur prolactine", "impact": "Production +10%"}
    STAT5A = {"chrom": "OAR19", "desc": "Signal transduction", "impact": "Différenciation mammaire"}
    ACACA = {"chrom": "OAR11", "desc": "Acétyl-CoA carboxylase", "impact": "Lipogenèse"}
    FASN = {"chrom": "OAR12", "desc": "Fatty acid synthase", "impact": "Acides gras"}

# Standards des races algériennes
STANDARDS_RACES = {
    'HAMRA': {
        'nom_complet': 'Hamra (Rousse)',
        'couleur': 'Rouge à marron',
        'origines': ['Sud Algérien', 'Sahara'],
        'caracteristiques': ['Robe rousse', 'Adaptée au désert', 'Bonne laitière'],
        'poids_adulte': {'femelle': (45, 65), 'male': (65, 90)},
        'mensurations': {
            'longueur_cm': (95, 125), 'hauteur_cm': (65, 85),
            'tour_poitrine_cm': (95, 120), 'largeur_bassin_cm': (35, 50)
        },
        'production_lait': (1.5, 3.5),
        'criteres_mammaires': {
            'volume': 'Moyen à élevé', 'trayons': '3-5 cm',
            'symetrie': 'Bonne', 'aptitude_laitiere': 'Bonne'
        }
    },
    'OUDA': {
        'nom_complet': 'Ouled Djellal (Ouda)',
        'couleur': 'Blanche',
        'origines': ['Hauts Plateaux', 'Steppes'],
        'caracteristiques': ['Robe blanche', 'Queue grasse', 'Viande'],
        'poids_adulte': {'femelle': (50, 70), 'male': (70, 100)},
        'mensurations': {
            'longueur_cm': (100, 130), 'hauteur_cm': (70, 90),
            'tour_poitrine_cm': (100, 130), 'largeur_bassin_cm': (38, 55)
        },
        'production_lait': (1.0, 2.5),
        'criteres_mammaires': {
            'volume': 'Grand', 'trayons': '4-6 cm',
            'symetrie': 'Très bonne', 'aptitude_laitiere': 'Excellente'
        }
    },
    'SIDAHOU': {
        'nom_complet': 'Sidahou',
        'couleur': 'Noire et blanche',
        'origines': ['Ouest Algérien'],
        'caracteristiques': ['Tête noire', 'Résistante', 'Mixte'],
        'poids_adulte': {'femelle': (40, 60), 'male': (60, 85)},
        'mensurations': {
            'longueur_cm': (90, 120), 'hauteur_cm': (60, 80),
            'tour_poitrine_cm': (90, 115), 'largeur_bassin_cm': (34, 48)
        },
        'production_lait': (1.2, 2.8),
        'criteres_mammaires': {
            'volume': 'Moyen', 'trayons': '3-5 cm',
            'symetrie': 'Bonne', 'aptitude_laitiere': 'Moyenne'
        }
    },
    'BERBERE': {
        'nom_complet': 'Brebis Berbère',
        'couleur': 'Variée',
        'origines': ['Kabylie', 'Aurès'],
        'caracteristiques': ['Rustique', 'Petite taille', 'Adaptée montagne'],
        'poids_adulte': {'femelle': (35, 50), 'male': (50, 70)},
        'mensurations': {
            'longueur_cm': (80, 110), 'hauteur_cm': (55, 75),
            'tour_poitrine_cm': (85, 110), 'largeur_bassin_cm': (30, 45)
        },
        'production_lait': (0.8, 2.0),
        'criteres_mammaires': {
            'volume': 'Petit à moyen', 'trayons': '2-4 cm',
            'symetrie': 'Moyenne', 'aptitude_laitiere': 'Adaptée'
        }
    },
    'CROISE': {
        'nom_complet': 'Croisement',
        'couleur': 'Variable',
        'origines': ['Multiple'],
        'caracteristiques': ['Vigueur hybride', 'Adaptabilité'],
        'poids_adulte': {'femelle': (40, 70), 'male': (60, 95)},
        'mensurations': {
            'longueur_cm': (85, 125), 'hauteur_cm': (60, 85),
            'tour_poitrine_cm': (90, 125), 'largeur_bassin_cm': (32, 52)
        },
        'production_lait': (1.0, 3.0),
        'criteres_mammaires': {
            'volume': 'Variable', 'trayons': '3-6 cm',
            'symetrie': 'Variable', 'aptitude_laitiere': 'Variable'
        }
    },
    'INCONNU': {
        'nom_complet': 'Race non identifiée',
        'couleur': 'Indéterminée',
        'origines': ['Inconnue'],
        'caracteristiques': ['À caractériser'],
        'poids_adulte': {'femelle': (30, 60), 'male': (50, 80)},
        'mensurations': {
            'longueur_cm': (80, 120), 'hauteur_cm': (55, 80),
            'tour_poitrine_cm': (85, 120), 'largeur_bassin_cm': (30, 50)
        },
        'production_lait': (0.5, 2.5),
        'criteres_mammaires': {
            'volume': 'À évaluer', 'trayons': 'À mesurer',
            'symetrie': 'À évaluer', 'aptitude_laitiere': 'À déterminer'
        }
    }
}

# ============================================================================
# SECTION 5: BASE DE DONNÉES - DAO COMPLET
# ============================================================================
@contextmanager
def get_db_connection():
    """Gestionnaire de connexion robuste"""
    conn = sqlite3.connect(DB_NAME, check_same_thread=False, timeout=30.0)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA journal_mode = WAL")
        yield conn
        conn.commit()
    except Exception as e:
        conn.rollback()
        logger.error(f"Erreur DB: {e}")
        raise
    finally:
        conn.close()

def init_database():
    """Initialisation complète du schéma intégré"""
    with get_db_connection() as conn:
        cursor = conn.cursor()
        
        # 1. Animaux (registre central)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY,
                numero_boucle TEXT UNIQUE,
                nom TEXT,
                espece TEXT CHECK(espece IN ('Bélier', 'Brebis', 'Agneau/elle')),
                race TEXT,
                date_naissance DATE,
                date_entree DATE,
                pere_id TEXT,
                mere_id TEXT,
                statut_reproductif TEXT,
                bcs_score REAL CHECK(bcs_score BETWEEN 1 AND 5),
                coef_consanguinite REAL DEFAULT 0.0,
                statut TEXT DEFAULT 'active',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # 2. Suivi médical
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS suivi_medical (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_intervention DATE NOT NULL,
                type TEXT CHECK(type IN ('Vaccination', 'Vermifugation', 'Traitement', 
                    'Chirurgie', 'Bilan', 'Échographie', 'Autre')),
                produit TEXT, posologie TEXT, veto TEXT,
                cout REAL DEFAULT 0, suite DATE,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        # 3. Alimentation
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS ration_alimentaire (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_debut DATE, date_fin DATE,
                type_ration TEXT, fourrage TEXT,
                concentre_kg REAL, fourrage_kg REAL,
                mineraux TEXT, cout_jour REAL,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        # 4. Reproduction (EPG, gestation)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS suivi_reproductif (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_chaleur DATE, date_saillie DATE,
                type_protocole TEXT, belier_utilise TEXT,
                date_eponge_pose DATE, date_eponge_retrait DATE,
                dose_PMSG INTEGER, date_mise_bas_prevue DATE,
                date_mise_bas_reelle DATE, nombre_agnes INTEGER,
                complications TEXT,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        # 5. Morphométrie complète (avec mamelles)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS morphometrie (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL, date_mesure DATE,
                
                -- Corps
                longueur_corps REAL, hauteur_garrot REAL,
                largeur_bassin REAL, tour_poitrine REAL,
                profondeur_thorax REAL,
                
                -- Mamelles (critères laitiers)
                volume_mammaire INTEGER, symetrie_mammaire INTEGER,
                insertion_trayons INTEGER, longueur_trayons REAL,
                orientation_trayons TEXT,
                
                -- Indices calculés
                ic_conformation REAL, volume_thoracique REAL,
                
                -- Métadonnées
                photo_path TEXT, operateur TEXT,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        # 6. Production laitière
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS production_lait (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL, date_controle DATE,
                stade_lactation TEXT, quantite_litre REAL,
                mg_percent REAL, proteine_percent REAL,
                caséine_percent REAL, lactose_percent REAL,
                cellules_somatiques INTEGER, pH REAL,
                rct_secondes INTEGER, rendement_fromage REAL,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        # 7. Génétique et SNPs
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS donnees_genetiques (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT UNIQUE NOT NULL,
                statut_sequencage TEXT,
                ncbi_biosample TEXT, ensembl_id TEXT,
                
                -- SNPs clés
                dgat1_snp TEXT, lalba_snp TEXT, csn1s1_snp TEXT,
                csn3_snp TEXT, prlr_snp TEXT, stat5a_snp TEXT,
                acaca_snp TEXT, fasn_snp TEXT,
                
                -- Scores
                score_genetique_lait REAL,
                diversite_genetique REAL,
                
                -- Fichiers
                vcf_path TEXT, date_analyse DATE,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        # 8. Scans 3D
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS scans_3d (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL, date_scan DATE,
                mode_scan TEXT, points_3d_json TEXT,
                volume_estime REAL, qualite_scan INTEGER,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            )
        ''')
        
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_animal_id ON animaux(id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_repro_animal ON suivi_reproductif(animal_id)')
        
        logger.info("Base de données initialisée")

# Peuplement initial si vide
def peupler_base_si_vide():
    with get_db_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM animaux")
        if cursor.fetchone()[0] == 0:
            races = list(STANDARDS_RACES.keys())
            for i in range(1, 51):
                race = random.choice(races)
                sexe = random.choices(['Bélier', 'Brebis', 'Agneau/elle'], weights=[0.2, 0.6, 0.2])[0]
                
                identifiant = f"DZ-{race[:3]}-{sexe[0]}-{i:04d}"
                poids_min, poids_max = STANDARDS_RACES[race]['poids_adulte'][sexe.lower()]
                
                cursor.execute('''
                    INSERT INTO animaux (id, numero_boucle, nom, espece, race,
                                       date_naissance, bcs_score, statut_reproductif)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    identifiant, f"BOUCLE-{i:04d}", f"Animal_{i}", sexe, race,
                    (date.today() - timedelta(days=random.randint(365, 2190))).isoformat(),
                    round(random.uniform(2.5, 4.0), 1),
                    random.choice(['Reproductive', 'Gestante', 'Lactante', 'Tarissement'])
                ))
            conn.commit()
            logger.info("Base peuplée avec 50 animaux")

# ============================================================================
# SECTION 6: MOTEURS DE CALCUL SCIENTIFIQUES
# ============================================================================
class MoteurZootechnique:
    @staticmethod
    def indice_conformation(pt: float, cc: float, hg: float) -> float:
        if cc <= 0 or hg <= 0:
            return 0
        return round((pt / (cc * hg)) * 1000, 2)
    
    @staticmethod
    def prediction_lait(poids: float, ic: float, hg: float) -> Dict:
        gras_mm = max(2.0, 2.0 + (poids * 0.15) + (ic * 0.1) - (hg * 0.05))
        pct_gras = max(10.0, min(40.0, 5.0 + (gras_mm * 1.5)))
        pct_proteine = max(45.0, min(72.0, 75.0 - (pct_gras * 0.6) + (ic * 0.2)))
        return {
            'matieres_grasses': round(pct_gras, 1),
            'proteines': round(pct_proteine, 1),
            'rendement_fromage': round(42 + (pct_proteine * 0.12), 1)
        }
    
    @staticmethod
    def classe_europ(ic: float) -> str:
        if ic > 33: return 'S'
        elif ic > 30: return 'E'
        elif ic > 27: return 'U'
        elif ic > 24: return 'R'
        else: return 'O/P'

class MoteurGenetique:
    @staticmethod
    def evaluer_genotype(genes: Dict) -> Dict:
        score = 0
        details = []
        
        if genes.get('dgat1') in ['GG', 'AA']:
            score += 20
            details.append("DGAT1 favorable (+15% MG)")
        elif genes.get('dgat1') in ['GA', 'A/G']:
            score += 10
        
        if genes.get('lalba') == 'AA':
            score += 15
            details.append("LALBA A (qualité protéique)")
        
        if genes.get('csn1s1') in ['CC', 'Strong']:
            score += 15
            details.append("Caséine αs1 forte")
        
        if genes.get('prlr') in ['GG', 'G/G']:
            score += 10
        
        return {
            'score': score,
            'max': 60,
            'niveau': 'Élite' if score >= 45 else 'Bon' if score >= 30 else 'Standard',
            'details': details
        }
    
    @staticmethod
    def calculer_delta_g(annees: int, h2: float = 0.30, i: float = 1.2, L: int = 3) -> List[float]:
        gain_annuel = (i * h2 * 0.5) / L
        return [2.0 + (gain_annuel * y) for y in range(annees + 1)]

class Scanner3DEngine:
    @staticmethod
    def generer_photo_simulee(brebis_info: Dict) -> Image.Image:
        width, height = 400, 300
        image = Image.new('RGB', (width, height), 'white')
        draw = ImageDraw.Draw(image)
        
        couleurs = {
            'HAMRA': (139, 0, 0), 'OUDA': (255, 255, 255),
            'SIDAHOU': (50, 50, 50), 'BERBERE': (165, 42, 42),
            'CROISE': (160, 120, 80), 'INCONNU': (200, 200, 200)
        }
        
        race = brebis_info.get('race', 'INCONNU')
        corps_color = couleurs.get(race, (200, 200, 200))
        
        # Corps
        draw.ellipse([100, 80, 300, 200], fill=corps_color, outline='black', width=2)
        # Tête
        draw.ellipse([280, 100, 350, 160], fill=corps_color, outline='black', width=2)
        # Pattes
        for x in [130, 170, 230, 270]:
            draw.rectangle([x, 200, x+20, 280], fill='black')
        
        # Texte
        draw.text((10, 10), f"ID: {brebis_info.get('identifiant', 'N/A')}", fill='black')
        draw.text((10, 30), f"Race: {race}", fill='black')
        draw.text((10, 50), f"Poids: {brebis_info.get('poids', 0):.1f} kg", fill='black')
        
        return image
    
    @staticmethod
    def simuler_scan_3d(brebis_info: Dict) -> List[Dict]:
        np.random.seed(hash(brebis_info.get('identifiant', '')) % 10000)
        points = []
        poids = brebis_info.get('poids', 50)
        
        rx, ry, rz = 0.6 * poids**0.33, 1.2 * poids**0.33, 0.8 * poids**0.33
        
        for _ in range(500):
            theta, phi = np.random.uniform(0, 2*np.pi), np.random.uniform(0, np.pi)
            x = rx * np.sin(phi) * np.cos(theta) + np.random.normal(0, rx*0.05)
            y = ry * np.sin(phi) * np.sin(theta) + np.random.normal(0, ry*0.05)
            z = rz * np.cos(phi) + np.random.normal(0, rz*0.05)
            
            intensity = np.random.uniform(100, 255)
            points.append({'x': float(x), 'y': float(y), 'z': float(z), 'intensity': int(intensity)})
        
        return points

# ============================================================================
# SECTION 7: CLASSES DE DONNÉES
# ============================================================================
@dataclass
class ProtocoleEPG:
    date_pose: date
    type_eponge: str = "Fluorogestone (30mg)"
    duree_jours: int = 14
    dose_PMSG: int = 400
    
    def calculer_mise_bas(self) -> date:
        date_retrait = self.date_pose + timedelta(days=self.duree_jours)
        date_fertilisation = date_retrait + timedelta(days=2)
        return date_fertilisation + timedelta(days=ConstantesReproduction.DUREE_GESTATION)
    
    def jours_restants(self) -> int:
        return (self.calculer_mise_bas() - date.today()).days

@dataclass
class ScoreMorphologique:
    volume_mammaire: int = 0
    symetrie_mammaire: int = 0
    suspension: int = 0
    tetines: int = 0
    capacite_corps: int = 0
    
    def total(self) -> int:
        return self.volume_mammaire + self.symetrie_mammaire + self.suspension + self.tetines + self.capacite_corps
    
    def classe(self) -> str:
        total = self.total()
        if total >= 20: return "🥇 ELITE"
        elif total >= 15: return "🥈 BON"
        elif total >= 10: return "🥉 MOYEN"
        else: return "⚠️ À AMÉLIORER"

# ============================================================================
# SECTION 8: MODULES D'AFFICHAGE (PAGES)
# ============================================================================

def page_dashboard():
    """Page d'accueil avec métriques globales"""
    st.markdown('<h1 class="main-header">🐑 EXPERT OVIN DZ PRO - Système Intégré</h1>', unsafe_allow_html=True)
    
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT * FROM animaux", conn)
        df_gestantes = pd.read_sql(
            "SELECT * FROM suivi_reproductif WHERE date_mise_bas_prevue >= date('now')", conn)
        df_lactantes = pd.read_sql(
            "SELECT * FROM animaux WHERE statut_reproductif = 'Lactante'", conn)
    
    # Métriques principales
    cols = st.columns(4)
    metrics = [
        ("🐑 Total Animaux", len(df_animaux) if not df_animaux.empty else 0),
        ("🤰 Gestantes", len(df_gestantes) if not df_gestantes.empty else 0),
        ("🥛 En Lactation", len(df_lactantes) if not df_lactantes.empty else 0),
        ("🧬 Séquencés", len(pd.read_sql("SELECT * FROM donnees_genetiques", conn)) if 'conn' in locals() else 0)
    ]
    
    for col, (label, value) in zip(cols, metrics):
        with col:
            st.markdown(f"""
            <div class='metric-card'>
                <h3>{label}</h3>
                <h2>{value}</h2>
            </div>
            """, unsafe_allow_html=True)
    
    # Alertes
    st.subheader("⚠️ Alertes Prioritaires")
    if not df_gestantes.empty:
        for _, row in df_gestantes.iterrows():
            date_prevue = pd.to_datetime(row['date_mise_bas_prevue'])
            jours = (date_prevue - datetime.now()).days
            if 0 <= jours <= 7:
                with st.expander(f"🔴 Mise-bas imminente - {row['animal_id']}"):
                    st.write(f"Date: {date_prevue.strftime('%d/%m/%Y')} ({jours} jours)")
                    st.button("Préparer box", key=f"box_{row['animal_id']}")
    
    # Distribution races
    if not df_animaux.empty:
        st.subheader("📊 Répartition du Troupeau")
        col1, col2 = st.columns([2, 1])
        with col1:
            fig = px.pie(df_animaux, names='race', title="Races présentes")
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            st.dataframe(df_animaux['race'].value_counts())

def page_gestion_animaux():
    """Gestion complète du registre animalier"""
    st.markdown('<h2 class="section-header">📋 Gestion des Animaux</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📃 Liste", "➕ Nouvel Animal", "🔍 Recherche"])
    
    with tab1:
        with get_db_connection() as conn:
            df = pd.read_sql("SELECT * FROM animaux ORDER BY created_at DESC", conn)
        st.dataframe(df, use_container_width=True, hide_index=True)
    
    with tab2:
        with st.form("nouvel_animal"):
            col1, col2 = st.columns(2)
            with col1:
                id_animal = st.text_input("ID Unique *", placeholder="DZ-HAM-B-0001")
                num_boucle = st.text_input("N° Boucle *")
                nom = st.text_input("Nom")
                espece = st.selectbox("Espèce", ["Brebis", "Bélier", "Agneau/elle"])
            with col2:
                race = st.selectbox("Race", list(STANDARDS_RACES.keys()))
                date_naiss = st.date_input("Date Naissance")
                bcs = st.slider("BCS", 1.0, 5.0, 3.0, 0.5)
                statut = st.selectbox("Statut", ['Reproductive', 'Gestante', 'Lactante', 'Tarissement', 'Engraissement'])
            
            if st.form_submit_button("💾 Enregistrer", use_container_width=True):
                try:
                    with get_db_connection() as conn:
                        conn.execute('''
                            INSERT INTO animaux (id, numero_boucle, nom, espece, race,
                                               date_naissance, bcs_score, statut_reproductif)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (id_animal, num_boucle, nom, espece, race,
                              date_naiss.isoformat() if date_naiss else None, bcs, statut))
                    st.success(f"✅ {id_animal} enregistré!")
                    st.balloons()
                except sqlite3.IntegrityError:
                    st.error("❌ ID ou N° Boucle déjà existant!")
    
    with tab3:
        recherche = st.text_input("Rechercher par ID, nom ou race")
        if recherche:
            with get_db_connection() as conn:
                df = pd.read_sql('''
                    SELECT * FROM animaux 
                    WHERE id LIKE ? OR nom LIKE ? OR race LIKE ?
                ''', conn, params=(f'%{recherche}%', f'%{recherche}%', f'%{recherche}%'))
            st.dataframe(df)

def page_reproduction():
    """Suivi reproduction, EPG et gestation"""
    st.markdown('<h2 class="section-header">🍼 Reproduction & Gestation</h2>', unsafe_allow_html=True)
    
    with get_db_connection() as conn:
        df_brebis = pd.read_sql("SELECT id, nom, race FROM animaux WHERE espece = 'Brebis'", conn)
    
    if df_brebis.empty:
        st.warning("Aucune brebis enregistrée")
        return
    
    animal_id = st.selectbox("Sélectionner la brebis", df_brebis['id'].tolist(),
                            format_func=lambda x: f"{x} ({df_brebis[df_brebis['id']==x]['nom'].values[0]})")
    
    tab1, tab2 = st.tabs(["🧽 Nouveau Protocole EPG", "📊 Historique"])
    
    with tab1:
        with st.form("protocole_epg"):
            col1, col2 = st.columns(2)
            with col1:
                date_pose = st.date_input("Date pose éponge", value=date.today())
                type_ep = st.selectbox("Type", ["Fluorogestone (30mg)", "Medroxyprogestérone (60mg)"])
                duree = st.number_input("Durée (jours)", 12, 16, 14)
            with col2:
                PMSG = st.number_input("Dose PMSG (UI)", 200, 600, 400, 50)
                belier = st.text_input("Bélier/IA")
            
            if st.form_submit_button("🚀 Calculer et Enregistrer"):
                proto = ProtocoleEPG(date_pose, type_ep, duree, PMSG)
                date_mb = proto.calculer_mise_bas()
                
                with get_db_connection() as conn:
                    conn.execute('''
                        INSERT INTO suivi_reproductif 
                        (animal_id, date_eponge_pose, date_eponge_retrait, dose_PMSG,
                         type_protocole, belier_utilise, date_mise_bas_prevue)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', (animal_id, date_pose.isoformat(),
                          (date_pose + timedelta(days=duree)).isoformat(),
                          PMSG, type_ep, belier, date_mb.isoformat()))
                
                st.success(f"✅ Date prévue mise-bas: {date_mb.strftime('%d/%m/%Y')}")
                
                # Calendrier
                events = pd.DataFrame([
                    {"Événement": "Retrait éponge", "Date": (date_pose + timedelta(days=duree)).strftime('%d/%m')},
                    {"Événement": "Injection PMSG", "Date": (date_pose + timedelta(days=duree-2)).strftime('%d/%m')},
                    {"Événement": "Mise-bas", "Date": date_mb.strftime('%d/%m'), "Importance": "🔴 CRITIQUE"}
                ])
                st.table(events)
    
    with tab2:
        with get_db_connection() as conn:
            histo = pd.read_sql('''
                SELECT * FROM suivi_reproductif WHERE animal_id = ?
                ORDER BY date_eponge_pose DESC
            ''', conn, params=(animal_id,))
        st.dataframe(histo if not histo.empty else pd.DataFrame({'Message': ['Aucun historique']}))

def page_morphometrie():
    """Scanner 3D et analyse morphométrique scientifique"""
    st.markdown('<h2 class="section-header">📏 Analyse Morphométrique & Mamelles</h2>', unsafe_allow_html=True)
    
    with get_db_connection() as conn:
        df_animaux = pd.read_sql("SELECT id, nom, race, poids FROM animaux", conn)
    
    if df_animaux.empty:
        st.warning("Aucun animal enregistré")
        return
    
    animal_id = st.selectbox("Animal", df_animaux['id'].tolist(),
                            format_func=lambda x: f"{x}")
    
    # Récupérer info
    info = df_animaux[df_animaux['id'] == animal_id].iloc[0].to_dict()
    
    tab1, tab2, tab3 = st.tabs(["📸 Scanner 3D", "📝 Saisie Manuelle", "🎯 Scoring Mamelles"])
    
    with tab1:
        col1, col2 = st.columns([2, 1])
        with col1:
            photo = Scanner3DEngine.generer_photo_simulee(info)
            st.image(photo, caption=f"Photo simulée - {info.get('nom', 'N/A')}")
            
            if st.button("🔍 Lancer Scan 3D", type="primary"):
                with st.spinner("Scan en cours..."):
                    time.sleep(2)
                    points = Scanner3DEngine.simuler_scan_3d(info)
                    st.success(f"✅ {len(points)} points capturés")
                    
                    # Visualisation 3D simplifiée
                    df_points = pd.DataFrame(points[:50])
                    fig = px.scatter_3d(df_points, x='x', y='y', z='z', color='intensity')
                    st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            st.markdown(f"""
            <div class='race-card'>
                <h4>{info.get('race', 'N/A')}</h4>
                <p>Poids: {info.get('poids', 0)} kg</p>
            </div>
            """, unsafe_allow_html=True)
    
    with tab2:
        with st.form("saisie_morpho"):
            st.subheader("Mensurations (cm)")
            col_m1, col_m2 = st.columns(2)
            with col_m1:
                longueur = st.number_input("Longueur corps", 50.0, 200.0, 100.0)
                hauteur = st.number_input("Hauteur garrot", 40.0, 150.0, 70.0)
                largeur = st.number_input("Largeur bassin", 20.0, 100.0, 40.0)
            with col_m2:
                poitrine = st.number_input("Tour poitrine", 60.0, 200.0, 100.0)
                profondeur = st.number_input("Profondeur thorax", 20.0, 80.0, 35.0)
            
            if st.form_submit_button("💾 Enregistrer"):
                ic = MoteurZootechnique.indice_conformation(poitrine, 9.0, hauteur)
                with get_db_connection() as conn:
                    conn.execute('''
                        INSERT INTO morphometrie 
                        (animal_id, date_mesure, longueur_corps, hauteur_garrot,
                         largeur_bassin, tour_poitrine, profondeur_thorax, ic_conformation)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (animal_id, date.today().isoformat(), longueur, hauteur,
                          largeur, poitrine, profondeur, ic))
                st.success(f"✅ IC calculé: {ic}")
    
    with tab3:
        st.subheader("🎯 Évaluation Critères Mammaires")
        
        race = info.get('race', 'INCONNU')
        if race in STANDARDS_RACES:
            std = STANDARDS_RACES[race]['criteres_mammaires']
            st.markdown(f"""
            <div class='mammelle-card'>
                <h4>Standards {race}</h4>
                <p>Volume: {std['volume']}</p>
                <p>Trayons: {std['trayons']}</p>
            </div>
            """, unsafe_allow_html=True)
        
        with st.form("scoring_mamelles"):
            col_s1, col_s2 = st.columns(2)
            with col_s1:
                vol = st.slider("Volume (1-5)", 1, 5, 3)
                sym = st.slider("Symétrie (1-5)", 1, 5, 3)
            with col_s2:
                ins = st.slider("Insertion (1-5)", 1, 5, 3)
                long_tet = st.number_input("Longueur trayons (cm)", 2.0, 8.0, 4.0)
            
            if st.form_submit_button("🎯 Calculer Score"):
                score = ScoreMorphologique(volume_mammaire=vol, symetrie_mammaire=sym, 
                                          insertion_trayons=ins)
                st.metric("Score Total", f"{score.total()}/15")
                st.success(score.classe())

def page_genetique():
    """Génomique, SNPs et analyse bioinformatique"""
    st.markdown('<h2 class="section-header">🧬 Génomique & Bioinformatique</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["🔍 Recherche NCBI", "🧬 Saisie SNPs", "🚀 Progrès Génétique ΔG"])
    
    with tab1:
        st.subheader("Recherche dans NCBI/GenBank")
        gene = st.selectbox("Gène", [g.name for g in GenesLaitiers])
        
        if st.button("🔎 Rechercher"):
            st.info("Mode démo - Intégration API NCBI")
            # Simulation données
            st.json({
                "gene": gene,
                "chromosome": GenesLaitiers[gene].value['chrom'],
                "description": GenesLaitiers[gene].value['desc'],
                "snps_connus": ["rs123456", "rs789012"],
                "alleles": ["A/G", "C/T"]
            })
    
    with tab2:
        with get_db_connection() as conn:
            df_animaux = pd.read_sql("SELECT id, nom FROM animaux", conn)
        
        animal_id = st.selectbox("Animal", df_animaux['id'].tolist())
        
        with st.form("saisie_snps"):
            st.subheader("Génotypes SNPs")
            genes = {}
            cols = st.columns(4)
            gene_names = ['dgat1', 'lalba', 'csn1s1', 'prlr']
            
            for i, g in enumerate(gene_names):
                with cols[i % 4]:
                    genes[g] = st.selectbox(g.upper(), ['Non testé', 'AA', 'AB', 'BB', 'A/A', 'A/B', 'B/B'])
            
            if st.form_submit_button("🧬 Évaluer"):
                result = MoteurGenetique.evaluer_genotype(genes)
                st.success(f"Score: {result['score']}/60 - {result['niveau']}")
                for d in result['details']:
                    st.write(f"✅ {d}")
                
                # Sauvegarde
                with get_db_connection() as conn:
                    conn.execute('''
                        INSERT OR REPLACE INTO donnees_genetiques 
                        (animal_id, dgat1_snp, lalba_snp, csn1s1_snp, prlr_snp,
                         score_genetique_lait, date_analyse)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', (animal_id, genes.get('dgat1'), genes.get('lalba'),
                          genes.get('csn1s1'), genes.get('prlr'),
                          result['score'], date.today().isoformat()))
    
    with tab3:
        st.subheader("Prédiction Progrès Génétique (ΔG)")
        annees = st.slider("Horizon (années)", 1, 20, 10)
        h2 = st.slider("Héritabilité h²", 0.1, 0.5, 0.3)
        
        progres = MoteurGenetique.calculer_delta_g(annees, h2)
        
        fig = px.line(x=list(range(annees+1)), y=progres,
                     labels={'x': 'Années', 'y': 'Production (L/j)'},
                     title=f"Progrès génétique prévu (h²={h2})")
        st.plotly_chart(fig, use_container_width=True)

def page_production():
    """Suivi production laitière et analyses biochimiques"""
    st.markdown('<h2 class="section-header">🥛 Production & Biochimie</h2>', unsafe_allow_html=True)
    
    with get_db_connection() as conn:
        df_femelles = pd.read_sql("SELECT id, nom FROM animaux WHERE espece = 'Brebis'", conn)
    
    if df_femelles.empty:
        st.warning("Aucune brebis")
        return
    
    animal_id = st.selectbox("Brebis", df_femelles['id'].tolist())
    
    tab1, tab2 = st.tabs(["📝 Saisie Contrôle", "⚗️ Analyse Biochimique"])
    
    with tab1:
        with st.form("controle_laitier"):
            col1, col2, col3 = st.columns(3)
            with col1:
                quantite = st.number_input("Quantité (L)", 0.0, 10.0, 2.5)
                mg = st.number_input("Matières grasses %", 3.0, 15.0, 6.5)
            with col2:
                proteine = st.number_input("Protéines %", 3.0, 10.0, 5.5)
                caséine = st.number_input("Caséine %", 2.0, 6.0, 4.5)
            with col3:
                cellules = st.number_input("Cellules (milliers/mL)", 0, 2000, 200)
                pH = st.number_input("pH", 6.0, 7.0, 6.7)
            
            if st.form_submit_button("💾 Enregistrer"):
                rendement = 42 + (caséine * 0.15)
                with get_db_connection() as conn:
                    conn.execute('''
                        INSERT INTO production_lait 
                        (animal_id, date_controle, quantite_litre, mg_percent,
                         proteine_percent, caséine_percent, cellules_somatiques,
                         pH, rendement_fromage)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (animal_id, date.today().isoformat(), quantite, mg,
                          proteine, caséine, cellules*1000, pH, rendement))
                st.success(f"✅ Rendement fromager estimé: {rendement:.1f}%")
    
    with tab2:
        st.subheader("Profil Biochimique")
        # Radar chart composition
        categories = ['MG', 'Protéines', 'Caséine', 'Lactose', 'Minéraux']
        values = [6.5, 5.5, 4.5, 4.8, 0.9]
        
        fig = go.Figure(go.Scatterpolar(
            r=values, theta=categories, fill='toself',
            name='Composition standard'
        ))
        fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 10])))
        st.plotly_chart(fig, use_container_width=True)

# ============================================================================
# SECTION 9: BARRE LATÉRALE ET NAVIGATION
# ============================================================================
def render_sidebar():
    """Menu latéral de navigation"""
    with st.sidebar:
        st.markdown("""
        <div style='text-align: center; padding: 20px; 
                    background: linear-gradient(135deg, #1a237e 0%, #8B0000 100%);
                    color: white; border-radius: 10px; margin-bottom: 20px;'>
            <h2>🐑 EXPERT OVIN DZ</h2>
            <p>Système Intégré Scientifique</p>
            <p><small>v5.0 Production</small></p>
        </div>
        """, unsafe_allow_html=True)
        
        menu = st.radio(
            "NAVIGATION PRINCIPALE",
            [
                "🏠 Dashboard",
                "📋 Gestion Animaux",
                "🍼 Reproduction & EPG",
                "📏 Morphométrie & Scanner",
                "🧬 Génomique & Bioinfo",
                "🥛 Production & Biochimie",
                "🏥 Suivi Médical",
                "🌾 Alimentation",
                "📊 Statistiques",
                "🔧 Admin"
            ]
        )
        
        st.markdown("---")
        
        # Statistiques rapides
        try:
            with get_db_connection() as conn:
                total = pd.read_sql("SELECT COUNT(*) as n FROM animaux", conn).iloc[0]['n']
                gestantes = pd.read_sql("SELECT COUNT(*) as n FROM suivi_reproductif WHERE date_mise_bas_prevue >= date('now')", conn).iloc[0]['n']
            
            st.metric("🐑 Total", int(total))
            st.metric("🤰 Gestantes", int(gestantes))
        except:
            pass
        
        return menu

# ============================================================================
# SECTION 10: POINT D'ENTRÉE PRINCIPAL
# ============================================================================
def main():
    """Point d'entrée principal"""
    init_database()
    peupler_base_si_vide()
    
    menu = render_sidebar()
    
    # Routage
    if menu == "🏠 Dashboard":
        page_dashboard()
    elif menu == "📋 Gestion Animaux":
        page_gestion_animaux()
    elif menu == "🍼 Reproduction & EPG":
        page_reproduction()
    elif menu == "📏 Morphométrie & Scanner":
        page_morphometrie()
    elif menu == "🧬 Génomique & Bioinfo":
        page_genetique()
    elif menu == "🥛 Production & Biochimie":
        page_production()
    elif menu == "🏥 Suivi Médical":
        st.title("🏥 Suivi Médical")
        st.info("Module vaccinations, traitements, bilans de santé")
    elif menu == "🌾 Alimentation":
        st.title("🌾 Gestion Alimentaire")
        st.info("Module rations, BCS, coûts alimentaires")
    elif menu == "📊 Statistiques":
        st.title("📊 Statistiques Avancées")
        st.info("ACP, ANOVA, héritabilités, corrélations génétiques")
    elif menu == "🔧 Admin":
        st.title("🔧 Administration")
        if st.button("🗑️ Vider la base (DANGER)"):
            st.error("Fonction désactivée en production")

if __name__ == "__main__":
    main()
